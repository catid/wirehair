#!/usr/bin/env python3
"""Deterministic pure-Python selftests for the WH2 all-K survivor-mask holdout.

These tests never execute the codec binary: they exercise the sealed task
ledger, mask validation, the reducer's accept/reject gates on synthetic
receipts, thermal-receipt tamper rejection with the control-stage-only CPU
busy floor, canonical-JSON round trips, and the environment-ownership marker
lifecycle.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import time
import unittest
from unittest import mock
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import wh2_row_mask_allk_holdout as subject


MOMENT = datetime(2026, 7, 21, 0, 0, 0, tzinfo=timezone.utc)


class CanonicalJsonTest(unittest.TestCase):
    def test_round_trip_and_stability(self):
        sample = {"b": [1, 2], "a": {"z": "x", "k": None}, "n": 7}
        encoded = subject.canonical_json(sample)
        self.assertTrue(encoded.endswith(b"\n"))
        self.assertNotIn(b" ", encoded)
        decoded = json.loads(encoded.decode("ascii"))
        self.assertEqual(decoded, sample)
        self.assertEqual(subject.canonical_json(decoded), encoded)

    def test_nan_rejected(self):
        with self.assertRaises(ValueError):
            subject.canonical_json({"bad": float("nan")})

    def test_load_json_rejects_nonfinite_and_ambiguous_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "receipt.json"
            for payload in (
                    b'{"bad":NaN}\n',
                    b'{"bad":Infinity}\n',
                    b'{"bad":1e400}\n',
                    b'{"same":1,"same":2}\n'):
                with self.subTest(payload=payload):
                    path.write_bytes(payload)
                    with self.assertRaises(subject.CampaignError):
                        subject.load_json(path)

    def test_json_lines_and_read_jsonl(self):
        rows = [{"K": 2, "cause": "q>H"}, {"K": 3, "cause": "field_shortfall"}]
        data = subject.json_lines(rows)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rows.jsonl"
            path.write_bytes(data)
            self.assertEqual(subject.read_jsonl(path), rows)
            path.write_bytes(b'{"K": 2}\n')
            with self.assertRaises(subject.CampaignError):
                subject.read_jsonl(path)
            path.write_bytes(b'{"K":2}')
            with self.assertRaises(subject.CampaignError):
                subject.read_jsonl(path)
            path.write_bytes(b'{"K":NaN}\n')
            with self.assertRaises(subject.CampaignError):
                subject.read_jsonl(path)

    def test_write_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "immutable.json"
            subject.write_once(path, subject.canonical_json({"v": 1}))
            subject.write_once(path, subject.canonical_json({"v": 1}))
            with self.assertRaises(subject.CampaignError):
                subject.write_once(path, subject.canonical_json({"v": 2}))


class ArmCatalogTest(unittest.TestCase):
    def test_sealed_catalog(self):
        arms = subject.arms_catalog()
        self.assertEqual(
            [arm["name"] for arm in arms],
            ["p48_r3_pfx007", "p48_r3_sfx380", "p32_r7_pfx07f", "p32_r7_sfx3f8"])
        by_name = {arm["name"]: arm for arm in arms}
        self.assertEqual(by_name["p48_r3_pfx007"]["mask"], 0x007)
        self.assertEqual(by_name["p48_r3_sfx380"]["mask"], 0x380)
        self.assertEqual(by_name["p32_r7_pfx07f"]["mask"], 0x07F)
        self.assertEqual(by_name["p32_r7_sfx3f8"]["mask"], 0x3F8)
        self.assertEqual(by_name["p32_r7_pfx07f"]["mask_hex"], "0x07f")
        for arm in arms:
            self.assertEqual(bin(arm["mask"]).count("1"), arm["rows"])
            self.assertLess(arm["mask"], 1 << 10)
        self.assertEqual(sum(arm["canonical_suffix"] for arm in arms), 2)
        for arm in arms:
            if arm["role"] == "candidate":
                partner = by_name[arm["comparator"]]
                self.assertEqual(partner["role"], "comparator")
                self.assertEqual(partner["panel"], arm["panel"])
                self.assertTrue(partner["canonical_suffix"])
                self.assertFalse(arm["canonical_suffix"])

    def test_mask_validation_rejects(self):
        def mutated(index, **overrides):
            arms = [dict(arm) for arm in subject.ARMS]
            arms[index].update(overrides)
            return arms
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(0, mask=0x00F))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(0, mask=1 << 10))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(0, mask=0))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(1, mask=0x1C0))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(0, comparator="p32_r7_sfx3f8"))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(mutated(0, name="p48_r3_sfx380"))
        with self.assertRaises(subject.CampaignError):
            subject.validate_arms(list(subject.ARMS)[:3])

    def test_sealed_comparator_totals(self):
        self.assertEqual(subject.SOURCE_ARM_FAILURES, {
            "p48_r3": 685, "p32_r7": 687, "p48_r0": 735, "prod244": 691,
        })
        self.assertEqual(subject.CHUNK_MAX, 250)
        self.assertEqual(subject.DEFAULT_WORKERS, 126)
        self.assertEqual((subject.K_LO, subject.K_HI), (2, 64000))
        self.assertEqual(subject.CPU_BUSY_FLOOR_STAGES, ("control",))
        self.assertEqual(subject.CPU_BUSY_FLOOR, Decimal("90"))


class TaskLedgerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, cls.hard, _ = subject.synthetic_reduction_fixture()
        cls.tasks = subject.build_tasks(
            "/frozen/wirehair_v2_bench", cls.hard, k_lo=2, k_hi=501)

    def test_identity_and_cardinality(self):
        tasks = self.tasks
        self.assertEqual([task["job"] for task in tasks], list(range(len(tasks))))
        self.assertEqual(len({task["output_name"] for task in tasks}), len(tasks))
        self.assertEqual(
            sum(len(task["Ks"]) for task in tasks), 4 * (4 + 500 * 9))
        cartesian = {
            (task["arm"], task["stage"], K, task["seed_index"], task["schedule"])
            for task in tasks for K in task["Ks"]
        }
        self.assertEqual(len(cartesian), 4 * (4 + 500 * 9))

    def test_stage_order_and_chunk_bounds(self):
        stages = [task["stage"] for task in self.tasks]
        self.assertEqual(stages, sorted(stages, key=("hard", "control").index))
        for task in self.tasks:
            self.assertTrue(task["Ks"])
            self.assertLessEqual(len(task["Ks"]), subject.CHUNK_MAX)
        chunks = [
            task for task in self.tasks
            if task["stage"] == "control" and task["arm"] == "p48_r3_pfx007" and
            task["seed_index"] == 0 and task["schedule"] == "burst"
        ]
        self.assertEqual([len(task["Ks"]) for task in chunks], [250, 250])
        self.assertEqual(chunks[0]["Ks"][0], 2)
        self.assertEqual(chunks[1]["Ks"][-1], 501)

    def test_hard_stage_covers_exactly_the_dev_cells(self):
        expected = {subject.cell_tuple(cell) for cell in self.hard}
        for arm in subject.arms_catalog():
            observed = {
                (K, task["seed_index"], task["schedule"])
                for task in self.tasks
                if task["stage"] == "hard" and task["arm"] == arm["name"]
                for K in task["Ks"]
            }
            self.assertEqual(observed, expected)

    def test_argv_shape(self):
        by_arm = {}
        for task in self.tasks:
            by_arm.setdefault(task["arm"], task)
        for arm in subject.arms_catalog():
            task = by_arm[arm["name"]]
            argv = task["argv"]
            def value_of(flag):
                return argv[argv.index(flag) + 1]
            self.assertEqual(argv[1], "precodefail")
            self.assertEqual(value_of("--mixed-grouped-gf256-rows"), str(arm["rows"]))
            self.assertEqual(
                value_of("--mixed-grouped-gf256-row-mask"),
                "0x{:x}".format(arm["mask"]))
            self.assertEqual(value_of("--mixed-period"), str(arm["period"]))
            self.assertEqual(value_of("--mixed-geometry"), "shared-x")
            self.assertEqual(value_of("--loss"), "0.50")
            self.assertEqual(value_of("--bb-list"), "64")
            self.assertEqual(value_of("--trials"), "1")
            self.assertEqual(value_of("--seed"), task["seed"])
            self.assertEqual(value_of("--schedule"), task["schedule"])
            self.assertEqual(
                value_of("--N"), ",".join(str(K) for K in task["Ks"]))
            self.assertIn("--binary-dense-two-anchor", argv)

    def test_duplicate_hard_cell_rejected(self):
        duplicated = list(self.hard) + [dict(self.hard[0])]
        with self.assertRaises(subject.CampaignError):
            subject.build_tasks(
                "/frozen/wirehair_v2_bench", duplicated, k_lo=2, k_hi=50)


class BenchmarkCsvTest(unittest.TestCase):
    def setUp(self):
        self.task = {
            "period": 48, "rows": 3, "mask": 0x007,
            "seed_index": 0, "schedule": "adversarial", "Ks": [257],
        }
        preamble = " ".join(
            "{}={}".format(key, value)
            for key, value in subject.expected_preamble(self.task).items())
        self.prefix = (
            "# precodefail: " + preamble + "\n" +
            ",".join(subject.CSV_FIELDS) + "\n"
        )
        self.row = {field: "0" for field in subject.CSV_FIELDS}
        for field in subject.CSV_MILLI_FIELDS:
            self.row[field] = "0.000"
        self.row.update({
            "N": "257", "bb": "64", "heavy_family": "periodic",
            "mix_count": "2", "overhead": "0", "trials": "1",
            "success": "1", "rank_fail": "0", "error": "0",
            "fail_rate": "0.00000000",
            "inact_mu": "0.000", "inact_max": "0",
            "binary_def_mu": "0.000", "binary_def_max": "0",
            "heavy_gain_mu": "0.000", "heavy_gain_min": "0",
            "heavy_shortfall": "0", "first_rank_fail": "-1",
            "binary_def_hist": "0:1", "heavy_gain_hist": "0:1",
            "failure_trials": "",
            "active_packet_peel_seed_xor": "0x0",
        })

    def payload(self, fields):
        return (self.prefix + ",".join(fields) + "\n").encode("ascii")

    def test_exact_row_width_accepted(self):
        values = [self.row[field] for field in subject.CSV_FIELDS]
        rows = subject.parse_benchmark_csv(self.payload(values), self.task)
        self.assertEqual(len(rows), 1)

    def test_missing_and_extra_row_columns_rejected(self):
        values = [self.row[field] for field in subject.CSV_FIELDS]
        with self.assertRaises(subject.CampaignError):
            subject.parse_benchmark_csv(self.payload(values[:-1]), self.task)
        with self.assertRaises(subject.CampaignError):
            subject.parse_benchmark_csv(
                self.payload(values + ["unexpected"]), self.task)

    def test_reduction_metrics_are_validated_before_publication(self):
        for field, value in (
                ("N", "0257"),
                ("seed_attempt", "+0"),
                ("seed_attempt", "256"),
                ("block_xors_mu", "NaN"),
                ("block_muladds_mu", "-1.000"),
                ("block_muladds_mu", "-0.000"),
                ("block_xors_mu", "1.0000"),
                ("block_xors_mu", "1.001"),
                ("solve_ms_mu", "1e0")):
            row = dict(self.row)
            row[field] = value
            with self.subTest(field=field, value=value), \
                    self.assertRaises(subject.CampaignError):
                subject.parse_benchmark_csv(
                    self.payload([
                        row[name] for name in subject.CSV_FIELDS
                    ]), self.task)
        fractional_timing = dict(self.row)
        fractional_timing["solve_ms_mu"] = "1.234"
        subject.parse_benchmark_csv(
            self.payload([
                fractional_timing[name] for name in subject.CSV_FIELDS
            ]), self.task)

    def test_single_trial_aggregate_contradictions_rejected(self):
        for field, value in (
                ("fail_rate", "1.00000000"),
                ("binary_def_max", "13"),
                ("heavy_shortfall", "1"),
                ("first_rank_fail", "0"),
                ("failure_trials", "0"),
                ("binary_def_hist", "13:1")):
            row = dict(self.row)
            row[field] = value
            values = [row[name] for name in subject.CSV_FIELDS]
            with self.subTest(field=field), \
                    self.assertRaises(subject.CampaignError):
                subject.parse_benchmark_csv(
                    self.payload(values), self.task)

        forged_success = dict(self.row)
        forged_success.update({
            "inact_mu": "13.000", "inact_max": "13",
            "binary_def_mu": "13.000", "binary_def_max": "13",
            "binary_def_hist": "13:1",
        })
        values = [
            forged_success[name] for name in subject.CSV_FIELDS
        ]
        with self.assertRaises(subject.CampaignError):
            subject.parse_benchmark_csv(
                self.payload(values), self.task)
        for inactivated, binary_deficit, heavy_gain in (
                (0, 13, 0),
                (12, 1, 2)):
            row = dict(self.row)
            row.update({
                "inact_mu": "{}.000".format(inactivated),
                "inact_max": str(inactivated),
                "binary_def_mu": "{}.000".format(binary_deficit),
                "binary_def_max": str(binary_deficit),
                "binary_def_hist": "{}:1".format(binary_deficit),
                "heavy_gain_mu": "{}.000".format(heavy_gain),
                "heavy_gain_min": str(heavy_gain),
                "heavy_gain_hist": "{}:1".format(heavy_gain),
            })
            values = [row[name] for name in subject.CSV_FIELDS]
            with self.subTest(
                    inactivated=inactivated,
                    binary_deficit=binary_deficit,
                    heavy_gain=heavy_gain), \
                    self.assertRaises(subject.CampaignError):
                subject.parse_benchmark_csv(
                    self.payload(values), self.task)

    def test_consistent_single_trial_rank_failures_accepted(self):
        for binary_deficit, heavy_gain, heavy_shortfall, outcome in (
                (13, 0, 0, "q>H"),
                (12, 11, 1, "field_shortfall")):
            row = dict(self.row)
            row.update({
                "success": "0", "rank_fail": "1",
                "fail_rate": "1.00000000",
                "inact_mu": "{}.000".format(binary_deficit),
                "inact_max": str(binary_deficit),
                "binary_def_mu": "{}.000".format(binary_deficit),
                "binary_def_max": str(binary_deficit),
                "heavy_gain_mu": "{}.000".format(heavy_gain),
                "heavy_gain_min": str(heavy_gain),
                "heavy_shortfall": str(heavy_shortfall),
                "first_rank_fail": "0",
                "binary_def_hist": "{}:1".format(binary_deficit),
                "heavy_gain_hist": "{}:1".format(heavy_gain),
                "failure_trials": "0",
            })
            values = [row[name] for name in subject.CSV_FIELDS]
            with self.subTest(outcome=outcome):
                parsed = subject.parse_benchmark_csv(
                    self.payload(values), self.task)
                self.assertEqual(subject.classify(parsed[0]), outcome)


class SourceSummaryTest(unittest.TestCase):
    @staticmethod
    def make_summary():
        def fs_record(K, seed_index, schedule):
            return {"K": K, "seed_index": seed_index, "schedule": schedule,
                    "cause": "field_shortfall", "binary_def": 12, "heavy_gain": 11}

        def qh_record(K, seed_index, schedule):
            return {"K": K, "seed_index": seed_index, "schedule": schedule,
                    "cause": "q>H", "binary_def": 13, "heavy_gain": 0}

        arms = {}
        fs_by_arm = {
            "p48_r3": [fs_record(100 + index, index % 3,
                                 subject.SCHEDULES[index % 3])
                       for index in range(7)],
            "p32_r7": [fs_record(200 + index, index % 3,
                                 subject.SCHEDULES[index % 3])
                       for index in range(9)],
            "p48_r0": [],
            "prod244": [],
        }
        for offset, arm in enumerate(subject.SOURCE_ARMS):
            total = subject.SOURCE_ARM_FAILURES[arm]
            records = list(fs_by_arm[arm])
            index = 0
            base = 1000 + offset * 10000
            while len(records) < total:
                records.append(qh_record(
                    base + index // 9,
                    (index // 3) % 3,
                    subject.SCHEDULES[index % 3]))
                index += 1
            arms[arm] = {"failures": total, "failure_records": records}
        return {
            "schema": 3, "head": subject.SOURCE_HEAD,
            "K_range": [subject.K_LO, subject.K_HI],
            "K_count": subject.K_HI - subject.K_LO + 1,
            "cells_per_arm": (subject.K_HI - subject.K_LO + 1) * 9,
            "validation_issue_count": 0, "timing_promotional": False,
            "arms": arms,
        }

    def test_accept_and_hard_union(self):
        summary = self.make_summary()
        records = subject.validate_source_summary(summary)
        self.assertEqual(len(records["p48_r3"]), 685)
        self.assertEqual(len(records["p32_r7"]), 687)
        self.assertEqual(len(records["p48_r0"]), 735)
        self.assertEqual(len(records["prod244"]), 691)
        hard = subject.derive_hard_keys(records)
        self.assertEqual(len(hard), 16)
        self.assertEqual(
            sum(cell["source_arms"] == ["p48_r3"] for cell in hard), 7)
        self.assertEqual(
            sum(cell["source_arms"] == ["p32_r7"] for cell in hard), 9)

    def test_reject_total_changed(self):
        summary = self.make_summary()
        summary["arms"]["p48_r3"]["failure_records"].pop()
        summary["arms"]["p48_r3"]["failures"] = 684
        with self.assertRaises(subject.CampaignError):
            subject.validate_source_summary(summary)

    def test_reject_head_changed(self):
        summary = self.make_summary()
        summary["head"] = "0" * 40
        with self.assertRaises(subject.CampaignError):
            subject.validate_source_summary(summary)

    def test_reject_duplicate_record(self):
        summary = self.make_summary()
        records = summary["arms"]["prod244"]["failure_records"]
        records[-1] = dict(records[0])
        with self.assertRaises(subject.CampaignError):
            subject.validate_source_summary(summary)

    def test_reject_field_shortfall_receipt_changed(self):
        summary = self.make_summary()
        summary["arms"]["p48_r3"]["failure_records"][0]["binary_def"] = 11
        with self.assertRaises(subject.CampaignError):
            subject.validate_source_summary(summary)

    def test_reject_overlapping_union(self):
        summary = self.make_summary()
        p48 = summary["arms"]["p48_r3"]["failure_records"][0]
        target = summary["arms"]["p32_r7"]["failure_records"][0]
        target.update({key: p48[key] for key in ("K", "seed_index", "schedule")})
        records = subject.validate_source_summary(summary)
        with self.assertRaises(subject.CampaignError):
            subject.derive_hard_keys(records)


class ReducerTest(unittest.TestCase):
    def setUp(self):
        self.outcomes, self.hard, self.source_records = \
            subject.synthetic_reduction_fixture()

    def reduce(self, outcomes):
        return subject.reduce_outcomes(
            outcomes, self.hard, self.source_records, k_lo=2, k_hi=13)

    def test_accept_path(self):
        reduction = self.reduce(self.outcomes)
        self.assertTrue(reduction["accepted"])
        self.assertEqual(reduction["codec_errors"], 0)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        p32 = reduction["candidates"]["p32_r7_pfx07f"]
        self.assertEqual(
            p48["acceptance"]["raw_field_shortfalls"],
            {"candidate": 1, "comparator": 2})
        self.assertEqual(
            p32["acceptance"]["raw_field_shortfalls"],
            {"candidate": 0, "comparator": 2})
        self.assertEqual(
            p48["acceptance"]["headline_field_shortfalls"],
            {"candidate": 0, "comparator": 0})
        self.assertEqual(
            p48["acceptance"]["headline_total_failures"],
            {"candidate": 0, "comparator": 1})
        self.assertTrue(
            p48["acceptance"]["gate_actual_recovery_improved"])
        self.assertTrue(p48["acceptance"]["accepted"])
        self.assertTrue(p32["acceptance"]["accepted"])
        # h1 is selection material, while q1 independently proves a headline
        # repair outside the four development cells.
        self.assertEqual(p48["repairs_vs_comparator"], 2)
        self.assertEqual(p48["headline_repairs"], 1)
        self.assertEqual(p48["introductions_vs_comparator"], 0)
        self.assertEqual(p32["repairs_vs_comparator"], 3)
        self.assertEqual(p48["headline_common_success"], 103)
        self.assertEqual(p48["acceptance"]["common_success_xor_ratio"], "1")
        self.assertEqual(p48["acceptance"]["common_success_muladd_ratio"], "1")
        fix = p48["fix_confirmation"]
        self.assertEqual(fix["panel_source_failures_fixed"], 1)
        self.assertEqual(fix["panel_source_failures_residual"], 1)
        self.assertEqual(len(fix["cells"]), len(self.hard))
        arm = reduction["arms"]["p48_r3_sfx380"]
        self.assertEqual(arm["failures"], 3)
        self.assertEqual(arm["field_shortfalls"], 2)
        self.assertEqual(arm["causes"], {"field_shortfall": 2, "q>H": 1})
        self.assertEqual(arm["weak_K"], 3)
        self.assertEqual(
            arm["weak_K_multiplicity"],
            {"3": 1, "5": 1, "11": 1})
        self.assertEqual(
            [record["K"] for record in arm["failure_records"]], [3, 5, 11])

    def test_reject_headline_introduction(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "control", 4, 0, "burst")] = \
            subject.synth_cell("q>H")
        reduction = self.reduce(outcomes)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        self.assertFalse(reduction["accepted"])
        self.assertFalse(p48["acceptance"]["gate_zero_headline_introductions"])
        self.assertEqual(p48["headline_introductions"], 1)
        self.assertEqual(
            p48["headline_introduction_keys"],
            [{"K": 4, "seed_index": 0, "schedule": "burst"}])
        self.assertTrue(p48["acceptance"]["gate_field_shortfall_improved"])

    def test_reject_cost_above_tolerance(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "control", 6, 0, "burst")] = \
            subject.synth_cell("success", xors=Decimal("1100"))
        reduction = self.reduce(outcomes)
        self.assertFalse(reduction["accepted"])
        self.assertFalse(
            reduction["candidates"]["p48_r3_pfx007"]
            ["acceptance"]["gate_cost_within_tolerance"])

    def test_cost_boundary_is_inclusive(self):
        # 103 common successes at 1000 XORs each; +51.5 lands exactly on
        # the sealed 1.0005 ceiling, which is accepted (<=).
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "control", 6, 0, "burst")] = \
            subject.synth_cell("success", xors=Decimal("1051.5"))
        reduction = self.reduce(outcomes)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        self.assertTrue(p48["acceptance"]["gate_cost_within_tolerance"])
        self.assertEqual(
            Decimal(p48["acceptance"]["common_success_xor_ratio"]),
            Decimal("1.0005"))

    def test_reject_field_shortfall_tie(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "control", 4, 2, "adversarial")] = \
            subject.synth_cell("field_shortfall")
        reduction = self.reduce(outcomes)
        self.assertFalse(reduction["accepted"])
        self.assertFalse(
            reduction["candidates"]["p48_r3_pfx007"]
            ["acceptance"]["gate_field_shortfall_improved"])

    def test_reject_failure_cause_relabel_without_recovery(self):
        outcomes = dict(self.outcomes)
        hard_keys = {subject.cell_tuple(cell) for cell in self.hard}
        for key in subject.domain_iter(2, 13):
            comparator = outcomes[("p48_r3_sfx380", "control", *key)]
            outcomes[("p48_r3_pfx007", "control", *key)] = dict(comparator)
            if key in hard_keys:
                outcomes[("p48_r3_pfx007", "hard", *key)] = dict(comparator)
            if comparator["outcome"] != "field_shortfall":
                continue
            outcomes[("p48_r3_pfx007", "control", *key)] = \
                subject.synth_cell("q>H")
            if key in hard_keys:
                outcomes[("p48_r3_pfx007", "hard", *key)] = \
                    subject.synth_cell("q>H")
        reduction = self.reduce(outcomes)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        self.assertFalse(reduction["accepted"])
        self.assertTrue(
            p48["acceptance"]["gate_field_shortfall_improved"])
        self.assertFalse(
            p48["acceptance"]["gate_actual_recovery_improved"])
        self.assertEqual(p48["repairs_vs_comparator"], 0)

    def test_reject_hard_relabel_even_with_independent_recovery(self):
        outcomes = dict(self.outcomes)
        panel_keys = {
            subject.cell_tuple(record)
            for record in self.source_records["p48_r3"]
            if record["cause"] == "field_shortfall"
        }
        for key in panel_keys:
            outcomes[("p48_r3_pfx007", "control", *key)] = \
                subject.synth_cell("q>H")
            outcomes[("p48_r3_pfx007", "hard", *key)] = \
                subject.synth_cell("q>H")
        reduction = self.reduce(outcomes)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        self.assertFalse(reduction["accepted"])
        self.assertTrue(
            p48["acceptance"]["gate_field_shortfall_improved"])
        self.assertTrue(
            p48["acceptance"]["gate_actual_recovery_improved"])
        self.assertFalse(
            p48["acceptance"]["gate_dev_fix_confirmed"])
        self.assertEqual(
            p48["acceptance"]["panel_source_failures_fixed"], 0)
        self.assertEqual(
            p48["fix_confirmation"]["panel_source_failures_fixed"], 0)

    def test_reject_dev_union_introductions_that_erase_raw_gain(self):
        outcomes = dict(self.outcomes)
        for key in (
                (7, 2, "repair-only"),
                (9, 0, "burst")):
            outcomes[("p48_r3_pfx007", "control", *key)] = \
                subject.synth_cell("q>H")
            outcomes[("p48_r3_pfx007", "hard", *key)] = \
                subject.synth_cell("q>H")
        reduction = self.reduce(outcomes)
        p48 = reduction["candidates"]["p48_r3_pfx007"]
        self.assertFalse(reduction["accepted"])
        self.assertTrue(
            p48["acceptance"]["gate_actual_recovery_improved"])
        self.assertTrue(
            p48["acceptance"]["gate_dev_fix_confirmed"])
        self.assertFalse(
            p48["acceptance"]["gate_raw_total_recovery_improved"])
        self.assertEqual(
            p48["acceptance"]["raw_total_failures"],
            {"candidate": 3, "comparator": 3})

    def test_reject_negative_reduced_receipts(self):
        for field in ("xors", "muladds"):
            for value in (Decimal("-1"), Decimal("-0")):
                outcomes = dict(self.outcomes)
                outcomes[("p48_r3_pfx007", "control", 6, 0, "burst")] = \
                    subject.synth_cell("success", **{field: value})
                with self.subTest(field=field, value=value), \
                        self.assertRaises(subject.CampaignError):
                    self.reduce(outcomes)
        for field in ("seed_attempt", "inactivated", "binary_deficit"):
            outcomes = dict(self.outcomes)
            outcomes[("p48_r3_pfx007", "control", 6, 0, "burst")] = \
                subject.synth_cell("success", **{field: -1})
            with self.subTest(field=field), self.assertRaises(
                    subject.CampaignError):
                self.reduce(outcomes)

    def test_reject_comparator_source_deviation(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_sfx380", "control", 3, 0, "burst")] = \
            subject.synth_cell("success")
        outcomes[("p48_r3_sfx380", "hard", 3, 0, "burst")] = \
            subject.synth_cell("success")
        with self.assertRaises(subject.CampaignError):
            self.reduce(outcomes)

    def test_reject_hard_stage_divergence(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "hard", 3, 0, "burst")] = \
            subject.synth_cell("field_shortfall")
        with self.assertRaises(subject.CampaignError):
            self.reduce(outcomes)

    def test_reject_missing_and_extra_cells(self):
        missing = dict(self.outcomes)
        del missing[("p32_r7_pfx07f", "control", 2, 0, "burst")]
        with self.assertRaises(subject.CampaignError):
            self.reduce(missing)
        extra = dict(self.outcomes)
        extra[("p32_r7_pfx07f", "control", 999, 0, "burst")] = \
            subject.synth_cell("success")
        with self.assertRaises(subject.CampaignError):
            self.reduce(extra)

    def test_paired_receipt_deltas_are_reported_not_rejected(self):
        outcomes = dict(self.outcomes)
        outcomes[("p48_r3_pfx007", "control", 8, 2, "adversarial")] = \
            subject.synth_cell("success", seed_attempt=3, inactivated=9)
        reduction = self.reduce(outcomes)
        self.assertTrue(reduction["accepted"])
        deltas = reduction["candidates"]["p48_r3_pfx007"] \
            ["paired_receipt_deltas_vs_comparator"]
        self.assertEqual(deltas["seed_attempt"], {
            "cells": 1, "candidate_minus_comparator_sum": 2, "max_abs": 2})
        self.assertEqual(deltas["inactivated"], {
            "cells": 1, "candidate_minus_comparator_sum": 4, "max_abs": 4})
        self.assertEqual(deltas["binary_deficit"]["cells"], 0)


class ThermalGateTest(unittest.TestCase):
    def validate(self, root, design, tasks):
        return subject.validate_thermal(root, design, tasks)

    def test_live_thermal_snapshot_retries_only_incomplete_tail(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            path = base / "thermal.csv"
            complete = subject.make_thermal_bytes(("100.0",))
            header = complete.split(b"\n", 1)[0] + b"\n"

            for raw in (b"", header[:-1], complete + b"torn"):
                path.write_bytes(raw)
                self.assertIsNone(subject.read_thermal_snapshot(
                    path, allow_incomplete=True))
                with self.assertRaises(subject.CampaignError):
                    subject.read_thermal_rows(path)

            path.write_bytes(complete)
            rows = subject.read_thermal_snapshot(
                path, allow_incomplete=True)
            self.assertIsNotNone(rows)
            self.assertEqual(len(rows), 1)

            path.write_bytes(complete + b"malformed\n")
            with self.assertRaises(subject.CampaignError):
                subject.read_thermal_snapshot(
                    path, allow_incomplete=True)

            path.write_bytes(complete)
            replacement = base / "replacement.csv"
            replacement.write_bytes(complete)
            real_read = os.read
            replaced = False

            def replace_during_read(descriptor, size):
                nonlocal replaced
                chunk = real_read(descriptor, size)
                if chunk and not replaced:
                    replaced = True
                    os.replace(replacement, path)
                return chunk

            with mock.patch.object(
                    subject.os, "read", side_effect=replace_during_read), \
                    self.assertRaises(subject.CampaignError):
                subject.read_thermal_snapshot(
                    path, allow_incomplete=True)
            self.assertTrue(replaced)

            path.write_bytes(complete)
            changed_mode = False

            def chmod_during_read(descriptor, size):
                nonlocal changed_mode
                chunk = real_read(descriptor, size)
                if chunk and not changed_mode:
                    changed_mode = True
                    path.chmod(0o400)
                return chunk

            with mock.patch.object(
                    subject.os, "read", side_effect=chmod_during_read), \
                    self.assertRaises(subject.CampaignError):
                subject.read_thermal_snapshot(
                    path, allow_incomplete=True)
            self.assertTrue(changed_mode)

    def test_live_thermal_reader_retries_a_torn_observation(self):
        rows = [{"monotonic_s": "100.0"}]
        with mock.patch.object(
                subject, "read_thermal_snapshot",
                side_effect=(None, rows)) as snapshots, \
                mock.patch.object(subject.time, "sleep") as sleep:
            self.assertIs(
                subject.read_live_thermal_rows(Path("/unused")), rows)
        self.assertEqual(snapshots.call_count, 2)
        sleep.assert_called_once_with(0.005)

    def test_thermal_snapshot_rejects_oversized_stream_boundedly(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(b"x" * 33)
            real_read = os.read
            requests = []

            def record_read(descriptor, size):
                requests.append(size)
                return real_read(descriptor, size)

            with mock.patch.object(
                    subject, "MAX_THERMAL_STREAM_BYTES", 32), \
                    mock.patch.object(
                        subject.os, "read", side_effect=record_read), \
                    self.assertRaises(subject.CampaignError):
                subject.read_thermal_snapshot(
                    path, allow_incomplete=True)
            self.assertEqual(requests, [33])

    def test_terminal_thermal_hashes_are_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            csv_path = base / "thermal.csv"
            stderr_path = base / "thermal.stderr"
            csv_path.write_bytes(b"x" * 33)
            stderr_path.write_bytes(b"y" * 17)
            with mock.patch.object(
                    subject, "MAX_THERMAL_STREAM_BYTES", 32), \
                    self.assertRaises(subject.CampaignError):
                subject.bounded_thermal_sha256(csv_path)
            with mock.patch.object(
                    subject, "MAX_THERMAL_STDERR_BYTES", 16), \
                    self.assertRaises(subject.CampaignError):
                subject.bounded_thermal_stderr_sha256(stderr_path)

    def test_bounded_thermal_snapshot_close_error_is_attempted_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(subject.make_thermal_bytes(("100.0",)))
            real_close = os.close
            closes = []

            def close_then_report_error(descriptor):
                closes.append(descriptor)
                real_close(descriptor)
                raise OSError("injected close report")

            with mock.patch.object(
                    subject.os, "close",
                    side_effect=close_then_report_error), \
                    self.assertRaises(subject.CampaignError):
                subject.read_thermal_rows(path)
            self.assertEqual(len(closes), 1)

    def test_success_and_busy_floor_stage_scoping(self):
        with tempfile.TemporaryDirectory() as temporary:
            success = Path(temporary) / "success"
            success.mkdir()
            tasks, design = subject.make_success_fixture(success)
            receipt = self.validate(success, design, tasks)
            self.assertEqual(receipt["successful_segments"], 1)
            self.assertEqual(receipt["samples"], 3)

            hard_burst = Path(temporary) / "hard-burst"
            hard_burst.mkdir()
            tasks, design = subject.make_success_fixture(hard_burst, busy="20.0")
            receipt = self.validate(hard_burst, design, tasks)
            self.assertEqual(receipt["successful_segments"], 1)

            control_floor = Path(temporary) / "control-floor"
            control_floor.mkdir()
            tasks, design = subject.make_success_fixture(
                control_floor, stage="control", busy=str(subject.CPU_BUSY_FLOOR))
            receipt = self.validate(control_floor, design, tasks)
            self.assertEqual(receipt["successful_segments"], 1)

            control_low = Path(temporary) / "control-low"
            control_low.mkdir()
            tasks, design = subject.make_success_fixture(
                control_low, stage="control", busy="89.9")
            with self.assertRaises(subject.CampaignError):
                self.validate(control_low, design, tasks)

    def test_busy_floor_is_duration_weighted_over_benchmark_window(self):
        monotonic = [
            Decimal(value)
            for value in ("0", "1", "3.5", "4.5", "7", "8")
        ]
        busy = ("100", "100", "80", "100", "80", "100")
        rows = [{"cpu_busy_pct": value} for value in busy]
        # The unweighted six-row mean exceeds 90%, but the two long 80%
        # intervals dominate the actual benchmark window.
        self.assertGreaterEqual(
            sum(Decimal(value) for value in busy) / len(busy),
            subject.CPU_BUSY_FLOOR)
        with self.assertRaises(subject.CampaignError):
            subject.validate_successful_thermal_coverage(
                rows, {"monotonic": monotonic},
                start=Decimal("0.1"), end=Decimal("7.9"),
                stage="control", segment=0,
            )
        intervals = subject.validate_successful_thermal_coverage(
            rows, {"monotonic": monotonic},
            start=Decimal("0.1"), end=Decimal("7.9"),
            stage="hard", segment=0,
        )
        self.assertEqual(
            sum((duration for _, duration in intervals), Decimal(0)),
            Decimal("7.8"))

    def test_tampered_thermal_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Path(temporary) / "fixture"
            fixture.mkdir()
            tasks, design = subject.make_success_fixture(fixture)

            def tampered(name, monotonic, **kwargs):
                target = Path(temporary) / name
                shutil.copytree(fixture, target)
                csv_path = target / "thermal" / "segment000.csv"
                csv_path.write_bytes(subject.make_thermal_bytes(monotonic, **kwargs))
                final_path = subject.segment_path(target, 0, "final")
                final = subject.load_json(final_path)
                final["thermal_csv_sha256"] = subject.sha256_file(csv_path)
                final_path.write_bytes(subject.canonical_json(final))
                return target

            cadence = tampered("cadence", ("100.0", "101.0", "104.0"))
            with self.assertRaises(subject.CampaignError):
                self.validate(cadence, design, tasks)
            edac = tampered("edac", ("100.0", "101.0", "102.0"), edac_ce="1")
            with self.assertRaises(subject.CampaignError):
                self.validate(edac, design, tasks)
            blank = tampered("blank", ("100.0", "101.0", "102.0"), blank_dimm=True)
            with self.assertRaises(subject.CampaignError):
                self.validate(blank, design, tasks)
            blank_busy = tampered(
                "blank-busy", ("100.0", "101.0", "102.0"),
                blank_busy_indices=(0, 1))
            with self.assertRaises(subject.CampaignError):
                self.validate(blank_busy, design, tasks)
            hot_cpu = tampered(
                "hot-cpu", ("100.0", "101.0", "102.0"),
                cpu_temperature="9999")
            with self.assertRaises(subject.CampaignError):
                self.validate(hot_cpu, design, tasks)
            cold_dimm = tampered(
                "cold-dimm", ("100.0", "101.0", "102.0"),
                dimm_temperature="-9999")
            with self.assertRaises(subject.CampaignError):
                self.validate(cold_dimm, design, tasks)
            for name, overrides in (
                    ("offset-utc", {"utc": "2026-07-21T01:00:00.000+01:00"}),
                    ("bad-monotonic", {"monotonic_s": "1e2"}),
                    ("signed-zero-busy", {"cpu_busy_pct": "-0"}),
                    ("blank-mhz", {"cpu_avg_mhz": ""}),
                    ("signed-zero-load", {"load1": "-0"}),
                    ("noncanonical-counter", {"edac_ce": "+0"})):
                target = tampered(
                    name, ("100.0", "101.0", "102.0"),
                    row_overrides=overrides)
                with self.subTest(name=name), \
                        self.assertRaises(subject.CampaignError):
                    self.validate(target, design, tasks)
            with self.assertRaises(subject.CampaignError):
                subject.validate_thermal_rows([], segment=0)

    def test_failed_ready_segment_preserves_bounded_diagnostic_evidence(self):
        for diagnostic in ("torn_csv", "stderr"):
            with self.subTest(diagnostic=diagnostic), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / diagnostic
                root.mkdir()
                tasks, design = subject.make_success_fixture(root)
                intent0 = subject.load_json(
                    subject.segment_path(root, 0, "intent"))
                ready0 = subject.load_json(
                    subject.segment_path(root, 0, "ready"))
                final0_path = subject.segment_path(root, 0, "final")
                final0 = subject.load_json(final0_path)

                intent1 = dict(intent0)
                intent1["segment"] = 1
                subject.write_once(
                    subject.segment_path(root, 1, "intent"),
                    subject.canonical_json(intent1))
                ready1 = dict(ready0)
                ready1["segment"] = 1
                ready1["intent_sha256"] = subject.sha256_file(
                    subject.segment_path(root, 1, "intent"))
                subject.write_once(
                    subject.segment_path(root, 1, "ready"),
                    subject.canonical_json(ready1))
                shutil.copytree(
                    root / "attempts" / "segment000",
                    root / "attempts" / "segment001")
                manifest_sha256, file_count = \
                    subject.seal_attempt_manifest(
                        root, 1, {0}, require_complete=True)
                csv_path = root / "thermal" / "segment000.csv"
                stderr_path = root / "thermal" / "segment000.stderr"
                retry_csv = root / "thermal" / "segment001.csv"
                retry_stderr = root / "thermal" / "segment001.stderr"
                retry_csv.write_bytes(csv_path.read_bytes())
                retry_stderr.write_bytes(stderr_path.read_bytes())
                receipt_path = subject.job_receipt_path(root, 0)
                receipt = subject.load_json(receipt_path)
                receipt["segment"] = 1
                receipt_path.write_bytes(subject.canonical_json(receipt))
                final1 = dict(final0)
                final1.update({
                    "segment": 1,
                    "intent_sha256": subject.sha256_file(
                        subject.segment_path(root, 1, "intent")),
                    "ready_sha256": subject.sha256_file(
                        subject.segment_path(root, 1, "ready")),
                    "thermal_csv": retry_csv.name,
                    "thermal_csv_sha256":
                        subject.bounded_thermal_sha256(retry_csv),
                    "thermal_stderr": retry_stderr.name,
                    "thermal_stderr_sha256":
                        subject.bounded_thermal_stderr_sha256(
                            retry_stderr),
                    "attempt_manifest":
                        "segment001.attempts.sha256",
                    "attempt_manifest_sha256": manifest_sha256,
                    "attempt_file_count": file_count,
                })
                subject.write_once(
                    subject.segment_path(root, 1, "final"),
                    subject.canonical_json(final1))

                if diagnostic == "torn_csv":
                    with csv_path.open("ab") as stream:
                        stream.write(b"torn")
                else:
                    stderr_path.write_bytes(b"sampler traceback\n")
                final0.update({
                    "state": "failed",
                    "failures": [{"status": "sampler_failed"}],
                    "rolled_back_jobs": [0],
                    "published_jobs": [],
                    "sampler_returncode": 1,
                    "thermal_csv_sha256":
                        subject.bounded_thermal_sha256(csv_path),
                    "thermal_stderr_sha256":
                        subject.bounded_thermal_stderr_sha256(stderr_path),
                })
                final0_path.write_bytes(subject.canonical_json(final0))
                receipt = self.validate(root, design, tasks)
                self.assertEqual(receipt["failed_segments"], 1)
                self.assertEqual(receipt["successful_segments"], 1)
                self.assertEqual(receipt["samples"], 3)

    def test_interrupted_segment_reconciliation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "interrupted"
            subject.make_private_directory(root)
            for directory in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / directory)
            task = {"job": 0, "stage": "hard", "output_name": "retry.csv",
                    "argv": ["/bin/true"]}
            intent = {
                "schema": subject.SCHEMA + ".segment_intent", "segment": 0,
                "boot_id": subject.required_boot_id(),
                "controller_pid": 99999998,
                "controller_start_ticks": 1,
                "stage": "hard", "jobs": [0],
                "jobs_sha256": subject.sha256_bytes(subject.json_lines([0])),
                "workers": 1,
                "retry_policy": "stage-atomic non-selective retry",
            }
            subject.write_once(
                subject.segment_path(root, 0, "intent"),
                subject.canonical_json(intent))
            subject.make_private_directory(root / "attempts" / "segment000")
            (root / "raw" / "retry.csv").write_bytes(b"failed output\n")
            (root / "stderr" / "retry.csv.stderr").write_bytes(b"failure\n")
            (root / "exit" / "retry.csv.exit").write_bytes(b"1\n")
            reconciled = subject.reconcile_incomplete_segments(root, [task])
            self.assertEqual(len(reconciled), 1)
            self.assertEqual(reconciled[0]["state"], "interrupted")
            self.assertEqual(subject.completed_jobs(root, [task]), set())
            for directory, name in (
                    ("raw", "retry.csv"), ("stderr", "retry.csv.stderr"),
                    ("exit", "retry.csv.exit")):
                self.assertFalse((root / directory / name).exists())
            self.assertEqual(
                subject.reconcile_incomplete_segments(root, [task]), [])

    def test_interrupted_attempt_accepts_atomic_output_prefix_after_pid(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "partial-attempt"
            subject.make_private_directory(root)
            subject.make_private_directory(root / "segments")
            attempt = root / "attempts" / "segment000"
            subject.make_private_directory(attempt, parents=True)
            (attempt / "job00000.pid").write_bytes(subject.canonical_json({
                "pid": 99999999, "start_ticks": 1,
            }))
            (attempt / "job00000.stdout").write_bytes(b"partial output\n")
            manifest_sha256, file_count = subject.seal_attempt_manifest(
                root, 0, {0}, require_complete=False)
            final = {
                "attempt_manifest": "segment000.attempts.sha256",
                "attempt_manifest_sha256": manifest_sha256,
                "attempt_file_count": file_count,
            }
            subject.validate_attempt_manifest(
                root, 0, final, {0}, require_complete=False)
            with self.assertRaises(subject.CampaignError):
                subject.validate_attempt_manifest(
                    root, 0, final, {0}, require_complete=True)

    def test_attempt_output_without_pid_is_rejected(self):
        for suffixes in (("stderr",), ("pid", "stderr"), ("pid", "exit")):
            with self.subTest(suffixes=suffixes), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "orphan-output"
                subject.make_private_directory(root)
                subject.make_private_directory(root / "segments")
                attempt = root / "attempts" / "segment000"
                subject.make_private_directory(attempt, parents=True)
                for suffix in suffixes:
                    payload = (
                        subject.canonical_json({
                            "pid": 99999999, "start_ticks": 1,
                        }) if suffix == "pid" else
                        b"1\n" if suffix == "exit" else b"orphan\n")
                    (attempt / ("job00000." + suffix)).write_bytes(payload)
                with self.assertRaises(subject.CampaignError):
                    subject.seal_attempt_manifest(
                        root, 0, {0}, require_complete=False)

    def test_attempt_exit_receipt_must_be_canonical_and_bounded(self):
        invalid_signal = "{}\n".format(-subject.signal.NSIG).encode("ascii")
        for payload in (
                b"+1\n", b"01\n", b" 1\n", b"1", b"256\n",
                invalid_signal):
            with self.subTest(payload=payload), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "bad-exit"
                subject.make_private_directory(root)
                subject.make_private_directory(root / "segments")
                attempt = root / "attempts" / "segment000"
                subject.make_private_directory(attempt, parents=True)
                (attempt / "job00000.pid").write_bytes(subject.canonical_json({
                    "pid": 99999999, "start_ticks": 1,
                }))
                (attempt / "job00000.stdout").write_bytes(b"failed output\n")
                (attempt / "job00000.stderr").write_bytes(b"failure\n")
                (attempt / "job00000.exit").write_bytes(payload)
                manifest_sha256, file_count = subject.seal_attempt_manifest(
                    root, 0, {0}, require_complete=False)
                final = {
                    "attempt_manifest": "segment000.attempts.sha256",
                    "attempt_manifest_sha256": manifest_sha256,
                    "attempt_file_count": file_count,
                }
                with self.assertRaises(subject.CampaignError):
                    subject.validate_attempt_manifest(
                        root, 0, final, {0}, require_complete=False)

    def test_attempt_exit_receipt_accepts_maximum_signal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "signal-exit"
            subject.make_private_directory(root)
            subject.make_private_directory(root / "segments")
            attempt = root / "attempts" / "segment000"
            subject.make_private_directory(attempt, parents=True)
            (attempt / "job00000.pid").write_bytes(subject.canonical_json({
                "pid": 99999999, "start_ticks": 1,
            }))
            (attempt / "job00000.stdout").write_bytes(b"")
            (attempt / "job00000.stderr").write_bytes(b"")
            (attempt / "job00000.exit").write_bytes(
                "{}\n".format(-(subject.signal.NSIG - 1)).encode("ascii"))
            manifest_sha256, file_count = subject.seal_attempt_manifest(
                root, 0, {0}, require_complete=False)
            subject.validate_attempt_manifest(
                root, 0, {
                    "attempt_manifest": "segment000.attempts.sha256",
                    "attempt_manifest_sha256": manifest_sha256,
                    "attempt_file_count": file_count,
                }, {0}, require_complete=False)

    def test_attempt_pid_receipt_must_be_exact_canonical_schema(self):
        payloads = (
            b'{"pid":99999999,"start_ticks":1,"extra":0}\n',
            b'{\n  "pid": 99999999,\n  "start_ticks": 1\n}\n',
            b'{"pid":99999999,"start_ticks":1}',
        )
        for payload in payloads:
            with self.subTest(payload=payload), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "bad-pid"
                subject.make_private_directory(root)
                subject.make_private_directory(root / "segments")
                attempt = root / "attempts" / "segment000"
                subject.make_private_directory(attempt, parents=True)
                (attempt / "job00000.pid").write_bytes(payload)
                manifest_sha256, file_count = subject.seal_attempt_manifest(
                    root, 0, {0}, require_complete=False)
                with self.assertRaises(subject.CampaignError):
                    subject.validate_attempt_manifest(
                        root, 0, {
                            "attempt_manifest":
                                "segment000.attempts.sha256",
                            "attempt_manifest_sha256": manifest_sha256,
                            "attempt_file_count": file_count,
                        }, {0}, require_complete=False)

    def test_prior_boot_reconciliation_never_uses_pid_signal_authority(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "prior-boot"
            subject.make_private_directory(root)
            for directory in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / directory)
            current_boot = subject.required_boot_id()
            old_boot = "00000000-0000-0000-0000-000000000000"
            if old_boot == current_boot:
                old_boot = "11111111-1111-1111-1111-111111111111"
            task = {
                "job": 0, "stage": "hard",
                "output_name": "retry.csv", "argv": ["/bin/true"],
            }
            intent = {
                "schema": subject.SCHEMA + ".segment_intent",
                "segment": 0, "boot_id": old_boot,
                "controller_pid": os.getpid(),
                "controller_start_ticks":
                    subject.process_start_ticks(os.getpid()),
                "stage": "hard", "jobs": [0],
                "jobs_sha256":
                    subject.sha256_bytes(subject.json_lines([0])),
                "workers": 1,
                "retry_policy": "stage-atomic non-selective retry",
            }
            subject.write_once(
                subject.segment_path(root, 0, "intent"),
                subject.canonical_json(intent))
            subject.make_private_directory(
                root / "attempts" / "segment000")
            pid_path = root / "thermal" / "segment000.pid"
            pid_path.write_text("{}\n".format(os.getpid()), encoding="ascii")
            with mock.patch.object(
                    subject, "terminate_owned_segment_processes") as cleanup, \
                    mock.patch.object(
                        subject, "other_samplers",
                        return_value=[{
                            "pid": os.getpid(),
                            "command": "exact current sampler",
                            "i2c_devices": ["/dev/i2c-1"],
                        }]), \
                    self.assertRaises(subject.CampaignError):
                subject.reconcile_incomplete_segments(root, [task])
            cleanup.assert_not_called()
            self.assertTrue(pid_path.exists())
            self.assertFalse(
                subject.segment_path(root, 0, "final").exists())
            with mock.patch.object(
                    subject, "terminate_owned_segment_processes") as cleanup, \
                    mock.patch.object(
                        subject, "other_samplers",
                        return_value=[]) as samplers, \
                    mock.patch.object(subject, "process_alive") as alive, \
                    mock.patch.object(
                        subject, "discover_owned_processes") as discover, \
                    mock.patch.object(
                        subject, "terminate_verified_process") as terminate:
                reconciled = subject.reconcile_incomplete_segments(
                    root, [task])
            cleanup.assert_not_called()
            samplers.assert_called_once_with()
            alive.assert_not_called()
            discover.assert_not_called()
            terminate.assert_not_called()
            self.assertFalse(pid_path.exists())
            self.assertEqual(
                reconciled[0]["process_actions"][0]["action"],
                "boot_changed")

    def test_nonfinite_completed_job_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "nonfinite"
            root.mkdir()
            tasks, _ = subject.make_success_fixture(root)
            path = subject.job_receipt_path(root, int(tasks[0]["job"]))
            receipt = json.loads(path.read_text(encoding="ascii"))
            receipt["started_monotonic_s"] = float("nan")
            receipt["ended_monotonic_s"] = float("nan")
            path.write_text(
                json.dumps(receipt, separators=(",", ":")) + "\n",
                encoding="ascii")
            with self.assertRaises(subject.CampaignError):
                subject.completed_jobs(root, tasks)

    def test_boolean_job_identity_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "boolean-job"
            root.mkdir()
            tasks, _ = subject.make_success_fixture(root)
            path = subject.job_receipt_path(root, int(tasks[0]["job"]))
            receipt = subject.load_json(path)
            receipt["job"] = False
            receipt["returncode"] = False
            path.write_bytes(subject.canonical_json(receipt))
            with self.assertRaises(subject.CampaignError):
                subject.completed_jobs(root, tasks)

    def test_ready_pid_bound_and_boolean_segment_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            fixture = base / "fixture"
            fixture.mkdir()
            tasks, design = subject.make_success_fixture(fixture)

            oversized = base / "oversized"
            shutil.copytree(fixture, oversized)
            ready_path = subject.segment_path(oversized, 0, "ready")
            final_path = subject.segment_path(oversized, 0, "final")
            ready = subject.load_json(ready_path)
            ready["sampler_pid"] = subject.MAX_PROCESS_PID + 1
            ready_path.write_bytes(subject.canonical_json(ready))
            final = subject.load_json(final_path)
            final["ready_sha256"] = subject.sha256_file(ready_path)
            final_path.write_bytes(subject.canonical_json(final))
            with self.assertRaises(subject.CampaignError):
                self.validate(oversized, design, tasks)

            boolean_segment = base / "boolean-segment"
            shutil.copytree(fixture, boolean_segment)
            intent_path = subject.segment_path(
                boolean_segment, 0, "intent")
            ready_path = subject.segment_path(
                boolean_segment, 0, "ready")
            final_path = subject.segment_path(
                boolean_segment, 0, "final")
            intent = subject.load_json(intent_path)
            intent["segment"] = False
            intent["previously_complete"] = False
            intent_path.write_bytes(subject.canonical_json(intent))
            ready = subject.load_json(ready_path)
            ready["segment"] = False
            ready["intent_sha256"] = subject.sha256_file(intent_path)
            ready_path.write_bytes(subject.canonical_json(ready))
            final = subject.load_json(final_path)
            final["segment"] = False
            final["intent_sha256"] = subject.sha256_file(intent_path)
            final["ready_sha256"] = subject.sha256_file(ready_path)
            final["published_jobs"] = [False]
            final_path.write_bytes(subject.canonical_json(final))
            with self.assertRaises(subject.CampaignError):
                self.validate(boolean_segment, design, tasks)

    def test_final_stage_and_retry_policy_are_bound_to_intent(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            fixture = base / "fixture"
            fixture.mkdir()
            tasks, design = subject.make_success_fixture(fixture)
            for name, field, value in (
                    ("stage", "stage", "control"),
                    ("retry", "retry_policy", "selective retry")):
                with self.subTest(field=field):
                    target = base / name
                    shutil.copytree(fixture, target)
                    final_path = subject.segment_path(target, 0, "final")
                    final = subject.load_json(final_path)
                    final[field] = value
                    final_path.write_bytes(subject.canonical_json(final))
                    with self.assertRaises(subject.CampaignError):
                        self.validate(target, design, tasks)

    def test_rollback_rejects_symlinked_runtime_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "campaign"
            subject.make_private_directory(root)
            for directory in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / directory)
            outside = base / "outside"
            subject.make_private_directory(outside)
            victim = outside / "job.csv"
            victim.write_bytes(b"keep me\n")
            (root / "raw").rmdir()
            (root / "raw").symlink_to(outside, target_is_directory=True)
            task = {"job": 0, "output_name": "job.csv"}
            with self.assertRaises(subject.CampaignError):
                subject.rollback_segment_outputs(root, [task])
            self.assertEqual(victim.read_bytes(), b"keep me\n")


class OwnerMarkerTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        base = Path(self.temporary.name)
        self.marker_path = base / "owner.json"
        self.boot_path = base / "boot_id"
        self.boot_path.write_text("boot-1\n", encoding="ascii")
        self.entry = {"role": "controller", "pid": 1234, "start_ticks": 99}

    def publish(self, phase="executing", protected=None, ttl=48, now=MOMENT):
        if protected is None:
            protected = [self.entry]
        return subject.write_owner_marker(
            Path("/tmp/example-root"), phase, protected, ttl,
            marker_path=self.marker_path, boot_id_path=self.boot_path, now=now)

    def load(self, now=MOMENT):
        return subject.load_active_owner_marker(
            marker_path=self.marker_path, boot_id_path=self.boot_path, now=now)

    def test_publish_schema_and_canonical_bytes(self):
        marker = self.publish()
        self.assertEqual(
            self.marker_path.read_bytes(), subject.canonical_json(marker))
        self.assertEqual(marker["schema"], "wirehair.environment_owner.v1")
        self.assertEqual(marker["campaign_root"], "/tmp/example-root")
        self.assertEqual(marker["phase"], "executing")
        self.assertEqual(marker["boot_id"], "boot-1")
        self.assertEqual(marker["created_utc"], "2026-07-21T00:00:00+00:00")
        self.assertEqual(marker["expires_utc"], "2026-07-23T00:00:00+00:00")
        self.assertEqual(marker["protected"], [self.entry])
        self.assertRegex(
            marker["expires_utc"],
            r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\+00:00$")
        self.assertFalse(
            self.marker_path.with_name(self.marker_path.name + ".partial").exists())

    def test_lifecycle_active_then_complete(self):
        self.publish()
        active = self.load()
        self.assertIsNotNone(active)
        self.assertEqual(active["campaign_root"], "/tmp/example-root")
        self.assertEqual(subject.active_owner_live_identities(
            active, alive_probe=lambda pid: pid == 1234,
            start_ticks_probe=lambda pid: 99), [self.entry])
        self.assertEqual(subject.active_owner_live_identities(
            active, alive_probe=lambda pid: pid == 1234,
            start_ticks_probe=lambda pid: 100), [])
        alive = iter((True, False))
        self.assertEqual(subject.active_owner_live_identities(
            active, alive_probe=lambda _pid: next(alive),
            start_ticks_probe=lambda _pid: None), [])
        self.publish(phase="complete", protected=[])
        self.assertIsNone(self.load())

    def test_expiry_and_boot_inertness(self):
        self.publish(ttl=48)
        self.assertIsNotNone(self.load(now=MOMENT + timedelta(hours=47)))
        self.assertIsNone(self.load(now=MOMENT + timedelta(hours=49)))
        self.boot_path.write_text("boot-2\n", encoding="ascii")
        self.assertIsNone(self.load())

    def test_missing_marker_is_inert(self):
        self.assertIsNone(self.load())

    def test_malformed_marker_fails_closed(self):
        self.marker_path.write_text("{malformed", encoding="ascii")
        self.marker_path.chmod(0o644)
        with self.assertRaises(subject.CampaignError):
            self.load()
        self.marker_path.write_text(
            '{"schema":"wirehair.environment_owner.v1",'
            '"campaign_root":"/tmp/example-root",'
            '"phase":"executing","phase":"complete",'
            '"boot_id":"boot-1",'
            '"created_utc":"2026-07-21T00:00:00+00:00",'
            '"expires_utc":"2026-07-23T00:00:00+00:00",'
            '"protected":[]}\n',
            encoding="ascii")
        with self.assertRaises(subject.CampaignError):
            self.load()
        self.marker_path.write_bytes(subject.canonical_json({
            "schema": subject.OWNER_MARKER_SCHEMA,
            "campaign_root": "/tmp/example-root",
            "phase": "complete",
            "boot_id": "boot-1",
            "created_utc": "2026-07-21T00:00:00+00:00",
            "expires_utc": "2026-07-23T00:00:00+00:00",
            "protected": [self.entry],
        }))
        self.assertIsNone(self.load())

    def test_stale_scratch_does_not_block_marker_publication(self):
        stale = self.marker_path.with_name(
            self.marker_path.name + ".partial")
        stale.write_text("stale\n", encoding="ascii")
        marker = self.publish()
        self.assertEqual(self.load(), marker)
        self.assertEqual(stale.read_text(encoding="ascii"), "stale\n")
        self.marker_path.write_bytes(
            subject.canonical_json({"schema": "other.schema.v1"}))
        with self.assertRaises(subject.CampaignError):
            self.load()
        self.marker_path.write_bytes(subject.canonical_json({
            "schema": subject.OWNER_MARKER_SCHEMA,
            "campaign_root": "/tmp/example-root",
            "phase": "complete",
            "boot_id": "boot-1",
            "created_utc": "2026-07-21T00:00:00+00:00",
            "expires_utc": "2026-07-23T00:00:00+00:00",
            "protected": [
                {"role": "reader", "pid": True, "start_ticks": 1}],
        }))
        with self.assertRaises(subject.CampaignError):
            self.load()

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(subject.CampaignError):
            self.publish(phase="running")
        with self.assertRaises(subject.CampaignError):
            self.publish(ttl=0)
        with self.assertRaises(subject.CampaignError):
            self.publish(protected=[{"role": "controller", "pid": 1234}])
        with self.assertRaises(subject.CampaignError):
            self.publish(protected=[
                {"role": "controller", "pid": 1, "start_ticks": 0}])
        self.assertFalse(self.marker_path.exists())

    def test_nonblocking_owner_lock(self):
        lock_path = Path(self.temporary.name) / "owner.lock"
        first = subject.acquire_environment_lock(lock_path)
        self.assertIsNotNone(first)
        self.assertIsNone(subject.acquire_environment_lock(lock_path))
        try:
            subject.os.close(first)
            first = None
            second = subject.acquire_environment_lock(lock_path)
            self.assertIsNotNone(second)
            subject.os.close(second)
        finally:
            if first is not None:
                subject.os.close(first)

    def test_successful_reduce_never_mutates_environment_owner(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in (
                    "raw", "stderr", "exit", "segments", "thermal",
                    "attempts", "job_receipts", "frozen"):
                subject.make_private_directory(root / name)
            task = {
                "job": 0, "arm": "synthetic", "stage": "control",
                "seed_index": 0, "schedule": "burst",
                "output_name": "job.csv",
            }
            (root / "raw" / "job.csv").write_bytes(b"synthetic csv\n")
            (root / "stderr" / "job.csv.stderr").write_bytes(b"")
            (root / "exit" / "job.csv.exit").write_bytes(b"0\n")
            args = subject.argparse.Namespace(
                root=str(root), expected_receipts_sha256="0" * 64)
            owner_lookup = mock.Mock(return_value={
                "campaign_root": str(root), "phase": "executing",
                "protected": [{
                    "role": "controller", "pid": os.getpid(),
                    "start_ticks": subject.process_start_ticks(os.getpid()),
                }],
            })
            owner_mutation = mock.Mock()
            with mock.patch.object(
                    subject, "verify_root",
                    return_value=({"total_cells": 1}, [task])), \
                    mock.patch.object(
                        subject, "completed_jobs", return_value={0}), \
                    mock.patch.object(
                        subject, "validate_thermal", return_value={}), \
                    mock.patch.object(
                        subject, "parse_benchmark_csv",
                        return_value=[{
                            "N": "2", "seed_attempt": "0",
                            "inact_max": "0", "binary_def_max": "0",
                        }]), \
                    mock.patch.object(
                        subject, "classify", return_value="success"), \
                    mock.patch.object(
                        subject, "decimal_field", return_value=Decimal(0)), \
                    mock.patch.object(subject, "read_jsonl", return_value=[]), \
                    mock.patch.object(subject, "load_json", return_value={}), \
                    mock.patch.object(
                        subject, "validate_source_summary", return_value={}), \
                    mock.patch.object(
                        subject, "reduce_outcomes",
                        return_value={"accepted": True, "candidates": {}}), \
                    mock.patch.object(
                        subject, "load_active_owner_marker", owner_lookup), \
                    mock.patch.object(
                        subject, "write_owner_marker", owner_mutation):
                subject.command_reduce(args)
            owner_lookup.assert_not_called()
            owner_mutation.assert_not_called()
            self.assertTrue((root / "validated_summary.json").is_file())


class ProcessGuardTest(unittest.TestCase):
    def test_resume_finalizes_owner_after_last_segment_crash_window(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in (
                    "segments", "job_receipts", "raw", "stderr", "exit",
                    "attempts"):
                subject.make_private_directory(root / name)
            subject.segment_path(root, 0, "final").write_bytes(b"{}\n")
            task = {"job": 0, "stage": "control"}
            design = {"worker_count": 1, "owner_ttl_hours": 168}
            owner = {
                "campaign_root": str(root), "phase": "executing",
                "protected": [{
                    "role": "controller", "pid": 99999999,
                    "start_ticks": 1,
                }],
            }
            args = subject.argparse.Namespace(
                expected_receipts_sha256="0" * 64, preflight_only=False,
                resume=True, stage="all",
            )
            finalize = mock.Mock()
            validate = mock.Mock(return_value={})
            with mock.patch.object(
                    subject, "verify_root",
                    side_effect=((design, [task]), (design, [task]))), \
                    mock.patch.object(
                        subject, "load_active_owner_marker",
                        return_value=owner), \
                    mock.patch.object(
                        subject, "active_owner_live_identities",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "segment_indices", return_value=[0]), \
                    mock.patch.object(
                        subject, "reconcile_incomplete_segments",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "completed_jobs", return_value={0}), \
                    mock.patch.object(
                        subject, "other_samplers", return_value=[]), \
                    mock.patch.object(
                        subject, "validate_thermal", validate), \
                    mock.patch.object(
                        subject, "write_owner_marker", finalize):
                subject.command_launch_locked(args, root)
            validate.assert_called_once_with(root, design, [task])
            finalize.assert_called_once_with(root, "complete", [], 168)

    def test_complete_owner_finalization_rechecks_for_foreign_race(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in (
                    "segments", "job_receipts", "raw", "stderr", "exit",
                    "attempts"):
                subject.make_private_directory(root / name)
            subject.segment_path(root, 0, "final").write_bytes(b"{}\n")
            task = {"job": 0, "stage": "control"}
            design = {"worker_count": 1, "owner_ttl_hours": 168}
            initial_owner = {
                "campaign_root": str(root), "phase": "executing",
                "protected": [],
            }
            foreign_owner = {
                "campaign_root": "/tmp/foreign", "phase": "executing",
                "protected": [],
            }
            args = subject.argparse.Namespace(
                expected_receipts_sha256="0" * 64, preflight_only=False,
                resume=True, stage="all",
            )
            finalize = mock.Mock()
            with mock.patch.object(
                    subject, "verify_root",
                    side_effect=((design, [task]), (design, [task]))), \
                    mock.patch.object(
                        subject, "load_active_owner_marker",
                        side_effect=(initial_owner, foreign_owner)), \
                    mock.patch.object(
                        subject, "active_owner_live_identities",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "segment_indices", return_value=[0]), \
                    mock.patch.object(
                        subject, "reconcile_incomplete_segments",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "completed_jobs", return_value={0}), \
                    mock.patch.object(
                        subject, "other_samplers", return_value=[]), \
                    mock.patch.object(
                        subject, "write_owner_marker", finalize), \
                    self.assertRaises(subject.CampaignError):
                subject.command_launch_locked(args, root)
            finalize.assert_not_called()

    def test_complete_owner_finalization_refuses_live_protected_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in (
                    "segments", "job_receipts", "raw", "stderr", "exit",
                    "attempts"):
                subject.make_private_directory(root / name)
            subject.segment_path(root, 0, "final").write_bytes(b"{}\n")
            task = {"job": 0, "stage": "control"}
            design = {"worker_count": 1, "owner_ttl_hours": 168}
            owner = {
                "campaign_root": str(root), "phase": "executing",
                "protected": [{
                    "role": "reader", "pid": 99999999,
                    "start_ticks": 1,
                }],
            }
            args = subject.argparse.Namespace(
                expected_receipts_sha256="0" * 64, preflight_only=False,
                resume=True, stage="all",
            )
            finalize = mock.Mock()
            with mock.patch.object(
                    subject, "verify_root",
                    side_effect=((design, [task]), (design, [task]))), \
                    mock.patch.object(
                        subject, "load_active_owner_marker",
                        side_effect=(owner, owner)), \
                    mock.patch.object(
                        subject, "active_owner_live_identities",
                        side_effect=([], [{
                            "role": "reader", "pid": 99999999,
                            "start_ticks": 1,
                        }])), \
                    mock.patch.object(
                        subject, "segment_indices", return_value=[0]), \
                    mock.patch.object(
                        subject, "reconcile_incomplete_segments",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "completed_jobs", return_value={0}), \
                    mock.patch.object(
                        subject, "other_samplers", return_value=[]), \
                    mock.patch.object(
                        subject, "write_owner_marker", finalize), \
                    self.assertRaises(subject.CampaignError):
                subject.command_launch_locked(args, root)
            finalize.assert_not_called()

    def test_multi_entry_cleanup_attempts_every_child(self):
        processes = [object(), object(), object()]
        seen = []

        def cleanup(job, process):
            seen.append((job, process))
            if job in (0, 1):
                raise subject.CampaignError(
                    "injected cleanup failure {}".format(job))
            return "terminated"

        results, errors = subject.attempt_all_cleanups(
            list(enumerate(processes)), cleanup)
        self.assertEqual(seen, list(enumerate(processes)))
        self.assertEqual(results, [(2, processes[2], "terminated")])
        self.assertEqual(len(errors), 2)
        self.assertIn("job0:", errors[0])
        self.assertIn("job1:", errors[1])

    def test_identical_worker_exceptions_retain_each_job_identity(self):
        futures = []
        for _ in range(4):
            future = subject.concurrent.futures.Future()
            future.set_exception(subject.CampaignError("same failure"))
            futures.append(future)
        failures = [
            subject.benchmark_future_result(future, job)
            for job, future in zip((5, 14, 23, 32), futures)
        ]
        self.assertEqual(
            [failure["job"] for failure in failures], [5, 14, 23, 32])
        self.assertTrue(all(
            failure["status"] == "worker_exception" and
            failure["detail"] == "CampaignError('same failure')"
            for failure in failures))
        self.assertEqual(len(failures), 4)

    def test_worker_result_mismatched_job_is_bound_to_expected_future(self):
        future = subject.concurrent.futures.Future()
        future.set_result({"job": 99, "status": "success"})
        self.assertEqual(subject.benchmark_future_result(future, 7), {
            "job": 7, "status": "worker_exception",
            "detail": "worker returned a mismatched job identity",
        })

    def test_cancelled_worker_future_retains_job_identity(self):
        future = subject.concurrent.futures.Future()
        self.assertTrue(future.cancel())
        self.assertEqual(subject.benchmark_future_result(future, 8), {
            "job": 8, "status": "future_cancelled",
        })

    def test_discovered_cleanup_attempts_every_identity_before_failing(self):
        identities = [
            ("first", 701, 10, ["first"]),
            ("second", 702, 20, ["second"]),
            ("third", 703, 30, ["third"]),
        ]
        seen = []

        def terminate(pid, start_ticks, tokens):
            seen.append((pid, start_ticks, tokens))
            if pid in (701, 702):
                raise subject.CampaignError(
                    "injected cleanup failure {}".format(pid))
            return "terminated"

        with mock.patch.object(
                subject, "discover_owned_processes",
                return_value=identities), \
                mock.patch.object(
                    subject, "terminate_verified_process",
                    side_effect=terminate), \
                self.assertRaises(subject.CampaignError):
            subject.terminate_discovered_owned_processes(
                [("unused", lambda tokens: False)])
        self.assertEqual(seen, [
            (701, 10, ["first"]),
            (702, 20, ["second"]),
            (703, 30, ["third"]),
        ])

    def test_discovered_cleanup_rejects_same_process_command_change(self):
        identity = ("owned", 701, 10, ["owned"])
        with mock.patch.object(
                subject, "discover_owned_processes",
                return_value=[identity]), \
                mock.patch.object(
                    subject, "terminate_verified_process",
                    return_value="identity_changed"), \
                self.assertRaises(subject.CampaignError):
            subject.terminate_discovered_owned_processes(
                [("unused", lambda tokens: False)])

    def test_discovery_identity_probe_race_to_exit_is_inert(self):
        proc = mock.Mock()
        entry = mock.Mock()
        entry.name = "701"
        proc.iterdir.return_value = [entry]
        with mock.patch.object(
                subject, "Path", return_value=proc), \
                mock.patch.object(
                    subject, "process_tokens",
                    side_effect=(["owned"], None)), \
                mock.patch.object(
                    subject, "process_alive",
                    side_effect=(True, False)), \
                mock.patch.object(
                    subject, "process_start_ticks", return_value=None):
            self.assertEqual(subject.discover_owned_processes([
                ("owned", lambda tokens: tokens == ["owned"]),
            ]), [])

    def test_benchmark_binding_accepts_terminal_before_exec_observation(self):
        process = mock.Mock(pid=701)
        wrapper = ["python", "-c", "wrapper"]
        benchmark = ["/bin/true"]
        with mock.patch.object(
                subject, "process_start_ticks", return_value=20), \
                mock.patch.object(
                    subject, "process_tokens", return_value=wrapper), \
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=True):
            identity = subject.bind_benchmark_exec_or_terminal(
                process, 9, wrapper, benchmark, subject.threading.Event())
        self.assertEqual(identity, (20, False))

    def test_benchmark_binding_accepts_exit_during_live_exec_reproof(self):
        process = mock.Mock(pid=701)
        wrapper = ["python", "-c", "wrapper"]
        benchmark = ["/bin/true"]
        with mock.patch.object(
                subject, "process_start_ticks",
                side_effect=(20, 20, 20)), \
                mock.patch.object(
                    subject, "process_tokens",
                    side_effect=(benchmark, benchmark)), \
                mock.patch.object(
                    subject, "pidfd_has_exited",
                    side_effect=(False, True)):
            identity = subject.bind_benchmark_exec_or_terminal(
                process, 9, wrapper, benchmark, subject.threading.Event())
        self.assertEqual(identity, (20, True))

    def test_benchmark_binding_rejects_live_unexpected_command(self):
        process = mock.Mock(pid=701)
        wrapper = ["python", "-c", "wrapper"]
        benchmark = ["/bin/true"]
        with mock.patch.object(
                subject, "process_start_ticks", return_value=20), \
                mock.patch.object(
                    subject, "process_tokens", return_value=["unexpected"]), \
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=False), \
                self.assertRaises(subject.CampaignError):
            subject.bind_benchmark_exec_or_terminal(
                process, 9, wrapper, benchmark, subject.threading.Event())

    def test_benchmark_binding_rejects_lost_start_even_when_terminal(self):
        process = mock.Mock(pid=701)
        wrapper = ["python", "-c", "wrapper"]
        benchmark = ["/bin/true"]
        with mock.patch.object(
                subject, "process_start_ticks",
                side_effect=(20, None)), \
                mock.patch.object(
                    subject, "process_tokens", return_value=[]), \
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=True), \
                self.assertRaises(subject.CampaignError):
            subject.bind_benchmark_exec_or_terminal(
                process, 9, wrapper, benchmark, subject.threading.Event())

    def test_real_fast_benchmark_binding_stress(self):
        parent_start = subject.process_start_ticks(os.getpid())
        self.assertIsNotNone(parent_start)
        benchmark = ["/bin/true"]
        for _ in range(32):
            wrapper = [
                str(subject.PYTHON_PATH), "-I", "-c",
                subject.PDEATH_EXEC_CODE, str(os.getpid()),
                str(parent_start), *benchmark,
            ]
            process = subject.subprocess.Popen(
                wrapper, stdout=subject.subprocess.PIPE,
                stderr=subject.subprocess.PIPE, start_new_session=True)
            pidfd = os.pidfd_open(process.pid, 0)
            try:
                start_ticks, _ = subject.bind_benchmark_exec_or_terminal(
                    process, pidfd, wrapper, benchmark,
                    subject.threading.Event())
                self.assertGreaterEqual(start_ticks, 0)
                stdout, stderr = process.communicate(timeout=5.0)
                self.assertEqual((process.returncode, stdout, stderr),
                                 (0, b"", b""))
            finally:
                os.close(pidfd)
                if process.poll() is None:
                    process.kill()
                process.communicate(timeout=5.0)

    def test_terminal_before_exec_observation_publishes_full_success(self):
        task = {
            "job": 0, "stage": "hard", "arm": "p48_r3_pfx007",
            "output_name": "job00000.csv",
            "argv": ["/fake/wirehair_v2_bench", "precodefail"],
            "period": 48, "rows": 3, "mask": 0x007,
            "seed_index": 0, "schedule": "adversarial", "Ks": [257],
        }
        preamble = " ".join(
            "{}={}".format(key, value)
            for key, value in subject.expected_preamble(task).items())
        row = {field: "0" for field in subject.CSV_FIELDS}
        for field in subject.CSV_MILLI_FIELDS:
            row[field] = "0.000"
        row.update({
            "N": "257", "bb": "64", "heavy_family": "periodic",
            "mix_count": "2", "overhead": "0", "trials": "1",
            "success": "1", "rank_fail": "0", "error": "0",
            "fail_rate": "0.00000000",
            "inact_mu": "0.000", "inact_max": "0",
            "binary_def_mu": "0.000", "binary_def_max": "0",
            "heavy_gain_mu": "0.000", "heavy_gain_min": "0",
            "heavy_shortfall": "0", "first_rank_fail": "-1",
            "binary_def_hist": "0:1", "heavy_gain_hist": "0:1",
            "failure_trials": "",
            "active_packet_peel_seed_xor": "0x0",
        })
        payload = (
            "# precodefail: " + preamble + "\n" +
            ",".join(subject.CSV_FIELDS) + "\n" +
            ",".join(row[field] for field in subject.CSV_FIELDS) + "\n"
        ).encode("ascii")

        class FakeSampler:
            pid = 700
            returncode = None

            def poll(self):
                return self.returncode

        class FastBenchmark:
            pid = 701
            returncode = None

            def poll(self):
                return self.returncode

            def communicate(self, timeout=None):
                self.returncode = 0
                return payload, b""

        sampler = FakeSampler()
        benchmark = FastBenchmark()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            thermal_csv = root / "thermal" / "segment000.csv"
            thermal_bytes = subject.make_thermal_bytes(
                ("100.0", "101.0", "102.0"))
            pidfd = os.open("/dev/null", os.O_RDONLY)

            def launch_sampler(_command, _error_stream):
                thermal_csv.write_bytes(thermal_bytes)
                return sampler

            def wait_ready(
                    process, pid_file, sampler_path, csv_path,
                    controller_identity, wrapper_identity, identity_sink):
                identity_sink.append((702, 30, ["sampler-child"]))
                return 702, subject.read_thermal_rows(csv_path)

            def start_ticks(pid):
                return {700: 10, 701: 20}.get(pid, 40)

            def tokens(pid):
                if pid == 700:
                    return ["sampler-wrapper"]
                if pid == 701:
                    return [
                        str(subject.PYTHON_PATH), "-I", "-c",
                        subject.PDEATH_EXEC_CODE, str(os.getpid()), "40",
                        *task["argv"],
                    ]
                return ["other"]

            patches = (
                mock.patch.object(
                    subject.subprocess, "Popen",
                    return_value=benchmark),
                mock.patch.object(
                    subject, "launch_thermal_sampler",
                    side_effect=launch_sampler),
                mock.patch.object(subject, "validate_trusted_root_tool"),
                mock.patch.object(
                    subject, "process_start_ticks", side_effect=start_ticks),
                mock.patch.object(
                    subject, "process_tokens", side_effect=tokens),
                mock.patch.object(
                    subject.os, "pidfd_open", return_value=pidfd),
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=True),
                mock.patch.object(
                    subject, "wait_sampler_ready", side_effect=wait_ready),
                mock.patch.object(
                    subject, "validate_thermal_rows", return_value={}),
                mock.patch.object(
                    subject, "read_live_thermal_rows",
                    side_effect=lambda path:
                        subject.read_thermal_rows(path)),
                mock.patch.object(
                    subject, "wait_thermal_coverage", return_value=True),
                mock.patch.object(
                    subject, "validate_successful_thermal_coverage",
                    return_value=[]),
                mock.patch.object(subject, "assert_process_identity"),
                mock.patch.object(
                    subject, "sole_sampler_competitors", return_value={}),
                mock.patch.object(
                    subject, "assert_current_controller_owner"),
                mock.patch.object(subject, "write_owner_marker"),
                mock.patch.object(
                    subject, "stop_launched_sampler",
                    return_value=(0, [])),
            )
            for patcher in patches:
                patcher.start()
            try:
                final = subject.run_segment(
                    root, {"worker_count": 1, "owner_ttl_hours": 1},
                    [task], "hard", 0, False)
            finally:
                for patcher in reversed(patches):
                    patcher.stop()
                try:
                    os.close(pidfd)
                except OSError:
                    pass
            self.assertEqual(final["state"], "success", final)
            self.assertEqual(final["failures"], [])
            self.assertEqual(
                sorted(path.suffix for path in
                       (root / "attempts" / "segment000").iterdir()),
                [".exit", ".pid", ".stderr", ".stdout"])
            self.assertEqual(
                subject.load_json(
                    root / "attempts" / "segment000" / "job00000.pid"),
                {"pid": 701, "start_ticks": 20})
            self.assertEqual(
                (root / "attempts" / "segment000" /
                 "job00000.stdout").read_bytes(), payload)
            self.assertTrue(subject.job_receipt_path(root, 0).is_file())

    def exercise_benchmark_cleanup(self, transient):
        class FakeSampler:
            pid = 700
            returncode = None

            def poll(self):
                return self.returncode

        class FakeBenchmark:
            pid = 701
            returncode = None

            def poll(self):
                return self.returncode

            @staticmethod
            def communicate(timeout=None):
                return b"", b""

        sampler = FakeSampler()
        benchmark = FakeBenchmark()
        calls = []
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            task = {
                "job": 0, "stage": "hard",
                "output_name": "job00000.csv",
                "argv": ["/fake/wirehair_v2_bench", "precodefail"],
            }
            design = {"worker_count": 1, "owner_ttl_hours": 1}
            rows = [{"monotonic_s": str(time.monotonic())}]
            real_atomic_result = subject.atomic_result
            pidfd = os.open("/dev/null", os.O_RDONLY)
            receipt_injection = {"hit": False}

            def wait_ready(
                    process, pid_file, sampler_path, csv_path,
                    controller_identity,
                    wrapper_identity,
                    identity_sink):
                identity_sink.append((702, 30, ["sampler-child"]))
                return 702, rows

            def start_ticks(pid):
                return {700: 10, 701: 20}.get(pid, 40)

            def tokens(pid):
                if pid == 700:
                    return ["sampler-wrapper"]
                if pid == 701:
                    return list(task["argv"])
                return ["other"]

            def atomic_result(path, data):
                if path.name == "job00000.pid":
                    receipt_injection["hit"] = True
                    raise subject.CampaignError(
                        "injected PID receipt write failure")
                return real_atomic_result(path, data)

            def terminate(process, owned_pidfd):
                calls.append((process.pid, owned_pidfd))
                if not transient or len(calls) == 1:
                    raise subject.CampaignError(
                        "injected cleanup failure")
                process.returncode = -15
                return "terminated"

            patches = (
                mock.patch.object(
                    subject.subprocess, "Popen",
                    side_effect=[sampler, benchmark]),
                mock.patch.object(subject, "validate_trusted_root_tool"),
                mock.patch.object(
                    subject, "process_start_ticks", side_effect=start_ticks),
                mock.patch.object(
                    subject, "process_tokens", side_effect=tokens),
                mock.patch.object(
                    subject.os, "pidfd_open", return_value=pidfd),
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=False),
                mock.patch.object(
                    subject, "wait_sampler_ready", side_effect=wait_ready),
                mock.patch.object(
                    subject, "validate_thermal_rows", return_value={}),
                mock.patch.object(
                    subject, "read_live_thermal_rows", return_value=rows),
                mock.patch.object(subject, "assert_process_identity"),
                mock.patch.object(
                    subject, "sole_sampler_competitors", return_value={}),
                mock.patch.object(
                    subject, "assert_current_controller_owner"),
                mock.patch.object(subject, "write_owner_marker"),
                mock.patch.object(
                    subject, "atomic_result", side_effect=atomic_result),
                mock.patch.object(
                    subject, "terminate_direct_child_by_pidfd",
                    side_effect=terminate),
                mock.patch.object(
                    subject, "stop_launched_sampler",
                    return_value=(0, [])),
            )
            for patcher in patches:
                patcher.start()
            try:
                try:
                    final = subject.run_segment(
                        root, design, [task], "hard", 0, False)
                    outcome = "returned:" + str(final["state"])
                except BaseException as exc:
                    outcome = "raised:" + repr(exc)
                result = {
                    "outcome": outcome,
                    "benchmark_poll": benchmark.poll(),
                    "termination_calls": len(calls),
                    "receipt_injection_hit": receipt_injection["hit"],
                    "final_exists": subject.segment_path(
                        root, 0, "final").exists(),
                }
            finally:
                for patcher in reversed(patches):
                    patcher.stop()
                try:
                    os.close(pidfd)
                except OSError:
                    pass
            return result

    def test_persistent_benchmark_cleanup_fault_stays_incomplete(self):
        result = self.exercise_benchmark_cleanup(False)
        self.assertTrue(result["outcome"].startswith("raised:"))
        self.assertIsNone(result["benchmark_poll"])
        self.assertGreaterEqual(result["termination_calls"], 2)
        self.assertTrue(result["receipt_injection_hit"])
        self.assertFalse(result["final_exists"])

    def test_transient_benchmark_cleanup_fault_retries_before_final(self):
        result = self.exercise_benchmark_cleanup(True)
        self.assertEqual(result["outcome"], "returned:failed")
        self.assertEqual(result["benchmark_poll"], -15)
        self.assertEqual(result["termination_calls"], 2)
        self.assertTrue(result["receipt_injection_hit"])
        self.assertTrue(result["final_exists"])

    def test_success_revalidates_sealed_thermal_stream_after_sampler_stop(self):
        class FakeSampler:
            pid = 700
            returncode = None

            def poll(self):
                return self.returncode

        real_validate = subject.validate_thermal_rows
        real_coverage = subject.validate_successful_thermal_coverage
        for failure_kind in ("semantic", "coverage"):
            with self.subTest(failure_kind=failure_kind), \
                    tempfile.TemporaryDirectory() as temporary:
                sampler = FakeSampler()
                root = Path(temporary) / "campaign"
                subject.make_private_directory(root)
                for name in subject.RUNTIME_DIRECTORY_NAMES:
                    subject.make_private_directory(root / name)
                subject.make_private_directory(root / "frozen")
                thermal_csv = root / "thermal" / "segment000.csv"
                thermal_bytes = subject.make_thermal_bytes(
                    ("100.0", "101.0", "102.0"))
                validation_calls = []
                coverage_calls = []

                def launch_sampler(_command, _error_stream):
                    thermal_csv.write_bytes(thermal_bytes)
                    return sampler

                def wait_ready(
                        process, pid_file, sampler_path, csv_path,
                        controller_identity, wrapper_identity, identity_sink):
                    identity_sink.append((702, 30, ["sampler-child"]))
                    return 702, subject.read_thermal_rows(thermal_csv)

                def validate_rows(rows, *, segment):
                    validation_calls.append(len(rows))
                    if len(validation_calls) == 3:
                        return real_validate(rows, segment=segment)
                    return {}

                def validate_coverage(
                        rows, thermal, *, start, end, stage, segment):
                    coverage_calls.append((start, end))
                    if failure_kind == "coverage" and \
                            len(coverage_calls) == 2:
                        return real_coverage(
                            rows, thermal, start=start, end=end,
                            stage=stage, segment=segment)
                    return []

                def stop_and_corrupt(*_args):
                    if failure_kind == "semantic":
                        invalid = subject.make_thermal_bytes(
                            ("103.0",), edac_ce="1")
                        with thermal_csv.open("ab") as stream:
                            stream.write(invalid.split(b"\n", 1)[1])
                            stream.flush()
                            os.fsync(stream.fileno())
                    sampler.returncode = 0
                    return 0, []

                patches = (
                    mock.patch.object(
                        subject, "launch_thermal_sampler",
                        side_effect=launch_sampler),
                    mock.patch.object(subject, "validate_trusted_root_tool"),
                    mock.patch.object(
                        subject, "process_start_ticks", return_value=10),
                    mock.patch.object(
                        subject, "process_tokens",
                        return_value=["sampler-wrapper"]),
                    mock.patch.object(
                        subject, "wait_sampler_ready",
                        side_effect=wait_ready),
                    mock.patch.object(
                        subject, "validate_thermal_rows",
                        side_effect=validate_rows),
                    mock.patch.object(
                        subject, "read_live_thermal_rows",
                        side_effect=lambda _path:
                            subject.read_thermal_rows(thermal_csv)),
                    mock.patch.object(
                        subject, "wait_thermal_coverage",
                        return_value=True),
                    mock.patch.object(
                        subject, "validate_successful_thermal_coverage",
                        side_effect=validate_coverage),
                    mock.patch.object(subject, "assert_process_identity"),
                    mock.patch.object(
                        subject, "sole_sampler_competitors",
                        return_value={}),
                    mock.patch.object(
                        subject, "assert_current_controller_owner"),
                    mock.patch.object(subject, "write_owner_marker"),
                    mock.patch.object(
                        subject, "stop_launched_sampler",
                        side_effect=stop_and_corrupt),
                )
                for patcher in patches:
                    patcher.start()
                try:
                    final = subject.run_segment(
                        root, {"worker_count": 1, "owner_ttl_hours": 1},
                        [], "hard", 0, False)
                finally:
                    for patcher in reversed(patches):
                        patcher.stop()
                self.assertEqual(final["state"], "failed")
                self.assertEqual(
                    [failure["status"] for failure in final["failures"]],
                    ["sealed_thermal_validation_error"])
                self.assertEqual(final["rolled_back_jobs"], [])
                self.assertEqual(
                    final["thermal_csv_sha256"],
                    subject.bounded_thermal_sha256(thermal_csv))
                self.assertEqual(len(validation_calls), 3)
                self.assertEqual(len(coverage_calls), 1 if
                                 failure_kind == "semantic" else 2)

    def test_success_final_rechecks_validated_thermal_digest_and_empty_stderr(
            self):
        for mutation in ("csv", "stderr"):
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "campaign"
                subject.make_private_directory(root)
                for name in subject.RUNTIME_DIRECTORY_NAMES:
                    subject.make_private_directory(root / name)
                intent = {
                    "stage": "hard", "jobs": [],
                    "jobs_sha256": subject.sha256_bytes(
                        subject.json_lines([])),
                }
                subject.write_once(
                    subject.segment_path(root, 0, "intent"),
                    subject.canonical_json(intent))
                subject.make_private_directory(
                    root / "attempts" / "segment000")
                thermal_csv = root / "thermal" / "segment000.csv"
                thermal_stderr = \
                    root / "thermal" / "segment000.stderr"
                thermal_csv.write_bytes(subject.make_thermal_bytes(
                    ("100.0", "101.0", "102.0")))
                thermal_stderr.write_bytes(b"")
                validated = subject.bounded_thermal_sha256(thermal_csv)
                if mutation == "csv":
                    thermal_csv.write_bytes(subject.make_thermal_bytes(
                        ("200.0", "201.0", "202.0")))
                else:
                    thermal_stderr.write_bytes(b"late diagnostic\n")
                with self.assertRaises(subject.CampaignError):
                    subject.write_segment_final(
                        root, 0, intent, "success", failures=[],
                        rolled_back=[], process_actions=[],
                        jobs_ended_monotonic_s=1.0,
                        sampler_returncode=0,
                        validated_thermal_csv_sha256=validated)
                self.assertFalse(
                    subject.segment_path(root, 0, "final").exists())

    def test_registration_fallback_cleanup_failure_stays_incomplete(self):
        class FakeSampler:
            pid = 700
            returncode = None

            def poll(self):
                return self.returncode

        class FakeBenchmark:
            def __init__(self, pid):
                self.pid = pid
                self.returncode = None

            def poll(self):
                return self.returncode

            def communicate(self, timeout=None):
                if self.returncode is None:
                    self.returncode = 0
                return b"", b""

        sampler = FakeSampler()
        benchmarks = [FakeBenchmark(701), FakeBenchmark(702)]
        popen_values = [sampler, *benchmarks]
        popen_lock = subject.threading.Lock()
        registered = {"pid": None}
        fallback_started = subject.threading.Event()
        fallback_attempts = []
        pidfd = os.open("/dev/null", os.O_RDONLY)

        def popen(*args, **kwargs):
            with popen_lock:
                return popen_values.pop(0)

        def pidfd_open(pid, flags):
            if pid in (701, 702):
                registered["pid"] = pid
                if not fallback_started.wait(5.0):
                    raise subject.CampaignError(
                        "duplicate registration did not run")
            return pidfd

        def direct_cleanup(process):
            if process.pid != registered["pid"]:
                fallback_attempts.append(process.pid)
                fallback_started.set()
                raise subject.CampaignError(
                    "injected persistent fallback cleanup failure")
            process.returncode = -15
            return "terminated"

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            task = {
                "job": 0, "stage": "hard",
                "output_name": "job00000.csv",
                "argv": ["/fake/wirehair_v2_bench", "precodefail"],
            }
            rows = [{"monotonic_s": str(time.monotonic())}]

            def wait_ready(
                    process, pid_file, sampler_path, csv_path,
                    controller_identity, wrapper_identity, identity_sink):
                identity_sink.append((703, 30, ["sampler-child"]))
                return 703, rows

            patches = (
                mock.patch.object(
                    subject.subprocess, "Popen", side_effect=popen),
                mock.patch.object(subject, "validate_trusted_root_tool"),
                mock.patch.object(
                    subject, "process_start_ticks",
                    side_effect=lambda pid:
                        {700: 10, 701: 20, 702: 21}.get(pid, 30)),
                mock.patch.object(
                    subject, "process_tokens",
                    side_effect=lambda pid:
                        ["sampler-wrapper"] if pid == 700 else
                        list(task["argv"])),
                mock.patch.object(
                    subject.os, "pidfd_open", side_effect=pidfd_open),
                mock.patch.object(
                    subject, "pidfd_has_exited", return_value=False),
                mock.patch.object(
                    subject, "wait_sampler_ready", side_effect=wait_ready),
                mock.patch.object(
                    subject, "validate_thermal_rows", return_value={}),
                mock.patch.object(
                    subject, "read_live_thermal_rows", return_value=rows),
                mock.patch.object(subject, "assert_process_identity"),
                mock.patch.object(
                    subject, "sole_sampler_competitors", return_value={}),
                mock.patch.object(
                    subject, "assert_current_controller_owner"),
                mock.patch.object(subject, "write_owner_marker"),
                mock.patch.object(
                    subject, "terminate_direct_unreaped_child",
                    side_effect=direct_cleanup),
                mock.patch.object(
                    subject, "stop_launched_sampler",
                    return_value=(0, [])),
            )
            for patcher in patches:
                patcher.start()
            try:
                with self.assertRaises(subject.CampaignError):
                    subject.run_segment(
                        root, {"worker_count": 2, "owner_ttl_hours": 1},
                        [dict(task), dict(task)], "hard", 0, False)
            finally:
                for patcher in reversed(patches):
                    patcher.stop()
                try:
                    os.close(pidfd)
                except OSError:
                    pass
            self.assertGreaterEqual(len(fallback_attempts), 2)
            self.assertFalse(subject.segment_path(root, 0, "final").exists())

    def test_benchmark_wrapper_dies_with_launching_parent(self):
        launcher_code = (
            "import os,subprocess,sys\n"
            "raw=open('/proc/{}/stat'.format(os.getpid()),'rb').read()\n"
            "start=int(raw.rsplit(b') ',1)[1].split()[19])\n"
            "child=subprocess.Popen([sys.argv[2],'-I','-c',sys.argv[1],"
            "str(os.getpid()),str(start),'/bin/sleep','60'])\n"
            "print(child.pid,flush=True)\n"
        )
        launcher = subject.subprocess.Popen(
            [
                str(subject.PYTHON_PATH), "-I", "-c", launcher_code,
                subject.PDEATH_EXEC_CODE, str(subject.PYTHON_PATH),
            ],
            stdout=subject.subprocess.PIPE,
            text=True,
        )
        child_pid = None
        try:
            assert launcher.stdout is not None
            child_pid = int(launcher.stdout.readline().strip())
            launcher.stdout.close()
            launcher.wait(timeout=5.0)
            deadline = time.monotonic() + 2.0
            while time.monotonic() < deadline and \
                    subject.process_alive(child_pid):
                time.sleep(0.01)
            self.assertFalse(subject.process_alive(child_pid))
        finally:
            if launcher.poll() is None:
                launcher.kill()
                launcher.wait(timeout=5.0)
            if child_pid is not None and subject.process_alive(child_pid):
                os.kill(child_pid, 9)

    def test_privileged_sampler_supervisor_binds_and_reaps(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            forbidden = base / "forbidden"
            current_start = subject.process_start_ticks(os.getpid())
            self.assertIsNotNone(current_start)
            mismatch = subject.subprocess.run(
                [
                    str(subject.SUDO_PATH), "-n", str(subject.PYTHON_PATH),
                    "-I", "-c", subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
                    str(os.getpid()), str(int(current_start) + 1),
                    subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
                    str(subject.PYTHON_PATH), "-I", "-c",
                    "import pathlib,sys;"
                    "pathlib.Path(sys.argv[1]).write_text('bad')",
                    str(forbidden),
                ],
                stdout=subject.subprocess.PIPE,
                stderr=subject.subprocess.PIPE,
                timeout=5.0,
            )
            self.assertEqual(mismatch.returncode, 125)
            self.assertEqual((mismatch.stdout, mismatch.stderr), (b"", b""))
            self.assertFalse(forbidden.exists())

            marker = base / "child.pid"
            launcher_code = (
                "import os,subprocess,sys,time\n"
                "supervisor=bytes.fromhex(sys.argv[1]).decode()\n"
                "pdeath=bytes.fromhex(sys.argv[2]).decode()\n"
                "marker=sys.argv[3]\n"
                "raw=open('/proc/self/stat','rb').read()\n"
                "start=int(raw.rsplit(b') ',1)[1].split()[19])\n"
                "child_code=\"import os,pathlib,sys,time;"
                "pathlib.Path(sys.argv[1]).write_text(str(os.getpid()));"
                "time.sleep(60)\"\n"
                "command=['/usr/bin/sudo','-n','/usr/bin/python3.12',"
                "'-I','-c',supervisor,str(os.getpid()),str(start),"
                "pdeath.encode().hex(),'/usr/bin/python3.12','-I','-c',"
                "child_code,marker]\n"
                "wrapper=subprocess.Popen(command)\n"
                "deadline=time.monotonic()+5\n"
                "while time.monotonic()<deadline and "
                "not os.path.exists(marker): time.sleep(.01)\n"
                "if not os.path.exists(marker): raise SystemExit(2)\n"
                "print(wrapper.pid,open(marker).read(),flush=True)\n"
            )
            launcher = subject.subprocess.Popen(
                [
                    str(subject.PYTHON_PATH), "-I", "-c", launcher_code,
                    subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE.encode(
                        "utf-8").hex(),
                    subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
                    str(marker),
                ],
                stdout=subject.subprocess.PIPE,
                text=True,
            )
            wrapper_pid = child_pid = None
            try:
                assert launcher.stdout is not None
                wrapper_pid, child_pid = (
                    int(value)
                    for value in launcher.stdout.readline().split())
                launcher.stdout.close()
                launcher.wait(timeout=5.0)
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline and (
                        subject.process_alive(wrapper_pid) or
                        subject.process_alive(child_pid)):
                    time.sleep(0.05)
                self.assertFalse(subject.process_alive(wrapper_pid))
                self.assertFalse(subject.process_alive(child_pid))
            finally:
                if launcher.poll() is None:
                    launcher.kill()
                    launcher.wait(timeout=5.0)
                for pid in (child_pid, wrapper_pid):
                    if pid is not None and subject.process_alive(pid):
                        subject.subprocess.run([
                            str(subject.SUDO_PATH), "-n",
                            str(subject.KILL_PATH), "-KILL", str(pid),
                        ], check=False)

    def test_privileged_sampler_supervisor_honors_pre_spawn_stop(self):
        current_start = subject.process_start_ticks(os.getpid())
        self.assertIsNotNone(current_start)
        popen = mock.Mock()

        def install(signum, handler):
            if signum == subject.signal.SIGTERM:
                handler(signum, None)

        argv = [
            "supervisor", str(os.getpid()), str(current_start),
            subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
            "/bin/true",
        ]
        with mock.patch.object(subject.sys, "argv", argv), \
                mock.patch.object(
                    subject.signal, "signal", side_effect=install), \
                mock.patch.object(subject.subprocess, "Popen", popen), \
                self.assertRaises(SystemExit) as caught:
            exec(subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE, {})
        self.assertEqual(caught.exception.code, 125)
        popen.assert_not_called()

    def test_privileged_sampler_supervisor_kills_group_descendant(self):
        with tempfile.TemporaryDirectory() as temporary:
            marker = Path(temporary) / "group.pids"
            current_start = subject.process_start_ticks(os.getpid())
            self.assertIsNotNone(current_start)
            child_code = (
                "import os,pathlib,subprocess,sys;"
                "child=subprocess.Popen(['/bin/sleep','60']);"
                "pathlib.Path(sys.argv[1]).write_text("
                "str(os.getpid())+' '+str(child.pid))"
            )
            wrapper = subject.subprocess.Popen(
                [
                    str(subject.SUDO_PATH), "-n", str(subject.PYTHON_PATH),
                    "-I", "-c",
                    subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
                    str(os.getpid()), str(current_start),
                    subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
                    str(subject.PYTHON_PATH), "-I", "-c",
                    child_code, str(marker),
                ],
                stdout=subject.subprocess.PIPE,
                stderr=subject.subprocess.PIPE,
            )
            owned = {}
            try:
                deadline = time.monotonic() + 5.0
                values = []
                while time.monotonic() < deadline:
                    try:
                        values = marker.read_text(encoding="ascii").split()
                    except FileNotFoundError:
                        values = []
                    if len(values) == 2:
                        break
                    time.sleep(0.01)
                self.assertEqual(len(values), 2)
                for pid in (int(value) for value in values):
                    owned[pid] = subject.process_start_ticks(pid)
                stdout, stderr = wrapper.communicate(timeout=20.0)
                self.assertEqual(wrapper.returncode, 0)
                self.assertEqual((stdout, stderr), (b"", b""))
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline and any(
                        subject.process_alive(pid) for pid in owned):
                    time.sleep(0.05)
                self.assertTrue(owned)
                self.assertFalse(any(
                    subject.process_alive(pid) for pid in owned))
            finally:
                if wrapper.poll() is None:
                    subject.terminate_privileged_direct_unreaped_child(
                        wrapper)
                for pid, start_ticks in owned.items():
                    if start_ticks is not None and \
                            subject.process_start_ticks(pid) == start_ticks:
                        subject.subprocess.run([
                            str(subject.SUDO_PATH), "-n",
                            str(subject.KILL_PATH), "-KILL", str(pid),
                        ], check=False)

    def test_unrepresentable_pid_rejected_before_os_call(self):
        with mock.patch.object(subject.os, "kill") as kill:
            with self.assertRaises(subject.CampaignError):
                subject.process_alive(10 ** 100)
        kill.assert_not_called()

    def test_i2c_parser_and_sole_reader_invariant(self):
        self.assertEqual(subject.parse_fuser_i2c_result(
            Path("/dev/i2c-1"), 0, b" 41 17 41",
            b"/dev/i2c-1:         \n"), (17, 41))
        self.assertEqual(subject.parse_fuser_i2c_result(
            Path("/dev/i2c-2"), 1, b"", b""), ())
        with self.assertRaises(subject.CampaignError):
            subject.parse_fuser_i2c_result(
                Path("/dev/i2c-1"), 0, b"41\n",
                b"/dev/i2c-1:         \n")
        self.assertEqual(subject.validate_sole_sampler_inventory({
            71: ["/dev/i2c-1", "/dev/i2c-2"],
            91: ["/dev/i2c-1"],
        }, 71, 70), {91: ["/dev/i2c-1"]})
        self.assertEqual(subject.validate_sole_sampler_inventory({
            70: ["/dev/i2c-1"],
            71: ["/dev/i2c-1", "/dev/i2c-2"],
        }, 71, 70), {70: ["/dev/i2c-1"]})
        with self.assertRaises(subject.CampaignError):
            subject.validate_sole_sampler_inventory({
                71: ["/dev/i2c-1"],
            }, 71, 70)

    def test_sampler_command_positions(self):
        self.assertTrue(subject.is_wirehair_sampler_command([
            "sudo", "-n", "python3", "-I",
            "/tmp/wirehair_expo_thermal_sampler_hardened.py",
            "--csv=/tmp/t.csv", "--pid-file=/tmp/t.pid",
        ]))
        self.assertTrue(subject.is_wirehair_sampler_command([
            "bash", "-c",
            "python3 /tmp/wirehair_expo_thermal_sampler.py "
            "--csv /tmp/t.csv --pid-file /tmp/t.pid",
        ]))
        self.assertFalse(subject.is_wirehair_sampler_command([
            "python3", "/tmp/review.py",
            "/tmp/wirehair_expo_thermal_sampler.py",
            "--csv", "notes", "--pid-file", "fixture",
        ]))
        sampler = Path("/campaign/frozen/wirehair_expo_thermal_sampler.py")
        csv_path = Path("/campaign/thermal/segment000.csv")
        pid_path = Path("/campaign/thermal/segment000.pid")
        controller_identity = (1234, 99)
        sampled = subject.campaign_sampled_command(
            sampler, csv_path, pid_path)
        supervisor = [
            str(subject.PYTHON_PATH), "-I", "-c",
            subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
            str(controller_identity[0]), str(controller_identity[1]),
            subject.PDEATH_EXEC_CODE.encode("utf-8").hex(), *sampled,
        ]
        self.assertTrue(subject.campaign_sampler_command_matches([
            str(subject.SUDO_PATH), "-n", *supervisor,
        ], sampler, csv_path, pid_path, controller_identity))
        self.assertTrue(subject.campaign_sampler_command_matches(
            supervisor, sampler, csv_path, pid_path, controller_identity))
        self.assertFalse(subject.campaign_sampler_command_matches(
            sampled, sampler, csv_path, pid_path, controller_identity))
        self.assertFalse(subject.campaign_sampler_command_matches([
            "python3", "/tmp/review.py", str(sampler),
            "--csv", str(csv_path), "--pid-file", str(pid_path),
        ], sampler, csv_path, pid_path, controller_identity))

    def test_sampler_launch_is_detached_from_callers_terminal(self):
        launched = object()
        error_stream = object()
        command = [str(subject.SUDO_PATH), "-n", "/bin/true"]
        with mock.patch.object(
                subject.subprocess, "Popen",
                return_value=launched) as popen:
            self.assertIs(
                subject.launch_thermal_sampler(command, error_stream),
                launched)
        popen.assert_called_once_with(
            command,
            stdin=subject.subprocess.DEVNULL,
            stdout=subject.subprocess.DEVNULL,
            stderr=error_stream,
            start_new_session=True,
        )

    def test_sampler_pid_publication_waits_for_canonical_newline(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_path = base / "thermal.pid"
            descriptor = os.open(
                str(pid_path),
                os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o444)
            os.fchmod(descriptor, 0o444)
            publication = None
            try:
                pid, publication = subject.read_sampler_pid_publication(
                    pid_path)
                self.assertIsNone(pid)
                self.assertIsNotNone(publication)

                os.write(descriptor, b"701")
                os.fsync(descriptor)
                pid, publication = subject.read_sampler_pid_publication(
                    pid_path, publication)
                self.assertIsNone(pid)
                self.assertEqual(publication[3], b"701")
                with mock.patch.object(
                        subject.os, "pidfd_open") as pidfd_open:
                    self.assertIsNone(subject.bind_sampler_identity(
                        pid_path, base / subject.SAMPLER_NAME,
                        base / "thermal.csv", (1234, 99),
                        (700, 10, ["/unused"])))
                pidfd_open.assert_not_called()

                os.write(descriptor, b"\n")
                os.fsync(descriptor)
                pid, publication = subject.read_sampler_pid_publication(
                    pid_path, publication)
                self.assertEqual(pid, 701)
                self.assertEqual(publication[3], b"701\n")
            finally:
                subject.close_sampler_pid_publication(publication)
                os.close(descriptor)

    def test_sampler_pid_publication_rejects_mutation_and_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            for case in (
                    "mutation", "metadata", "replacement",
                    "disappearance"):
                with self.subTest(case=case):
                    pid_path = base / (case + ".pid")
                    descriptor = os.open(
                        str(pid_path),
                        os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                        0o444)
                    os.fchmod(descriptor, 0o444)
                    os.write(descriptor, b"70")
                    os.fsync(descriptor)
                    _pid, publication = \
                        subject.read_sampler_pid_publication(pid_path)
                    self.assertIsNotNone(publication)
                    if case == "mutation":
                        os.ftruncate(descriptor, 0)
                        os.lseek(descriptor, 0, os.SEEK_SET)
                        os.write(descriptor, b"71")
                        os.fsync(descriptor)
                    elif case != "metadata":
                        os.close(descriptor)
                        descriptor = -1
                        pid_path.unlink()
                        if case == "replacement":
                            replacement = os.open(
                                str(pid_path),
                                os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                                os.O_CLOEXEC, 0o444)
                            os.fchmod(replacement, 0o444)
                            os.write(replacement, b"70")
                            os.fsync(replacement)
                            os.close(replacement)
                    try:
                        if case == "metadata":
                            real_read = os.read
                            changed_mode = False

                            def chmod_during_read(value, size):
                                nonlocal changed_mode
                                chunk = real_read(value, size)
                                if chunk and not changed_mode:
                                    changed_mode = True
                                    os.fchmod(descriptor, 0o666)
                                return chunk

                            with mock.patch.object(
                                    subject.os, "read",
                                    side_effect=chmod_during_read), \
                                    self.assertRaises(subject.CampaignError):
                                subject.read_sampler_pid_publication(
                                    pid_path, publication)
                            self.assertTrue(changed_mode)
                        else:
                            with self.assertRaises(subject.CampaignError):
                                subject.read_sampler_pid_publication(
                                    pid_path, publication)
                    finally:
                        subject.close_sampler_pid_publication(publication)
                        if descriptor >= 0:
                            os.close(descriptor)

    def test_incomplete_pid_publication_close_error_is_not_retried(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_path = base / "thermal.pid"
            descriptor = os.open(
                str(pid_path),
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o444)
            os.fchmod(descriptor, 0o444)
            os.write(descriptor, b"70")
            os.fsync(descriptor)
            os.close(descriptor)
            real_close = os.close
            closes = []

            def close_then_report_error(value):
                closes.append(value)
                real_close(value)
                raise OSError("injected close report")

            with mock.patch.object(
                    subject.os, "close",
                    side_effect=close_then_report_error), \
                    self.assertRaises(subject.CampaignError):
                subject.bind_sampler_identity(
                    pid_path, base / subject.SAMPLER_NAME,
                    base / "thermal.csv", (1234, 99),
                    (700, 10, ["/unused"]))
            self.assertEqual(len(closes), 1)

    def test_sampler_pid_publication_rejects_noncanonical_content(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            for index, raw in enumerate((
                    b"0", b"01", b"+7\n", b"7 \n",
                    str(subject.MAX_PROCESS_PID + 1).encode("ascii"))):
                with self.subTest(raw=raw):
                    path = base / "receipt{}.pid".format(index)
                    descriptor = os.open(
                        str(path),
                        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                        0o444)
                    os.fchmod(descriptor, 0o444)
                    os.write(descriptor, raw)
                    os.fsync(descriptor)
                    os.close(descriptor)
                    with self.assertRaises(subject.CampaignError):
                        subject.read_sampler_pid_publication(path)

    def test_sampler_readiness_bounds_incomplete_pid_publication(self):
        class FakeSampler:
            @staticmethod
            def poll():
                return None

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_path = base / "thermal.pid"
            descriptor = os.open(
                str(pid_path),
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o444)
            os.fchmod(descriptor, 0o444)
            identities = []
            try:
                with mock.patch.object(
                        subject.time, "monotonic",
                        side_effect=(0.0, 0.1, 13.0)), \
                        mock.patch.object(subject.time, "sleep"), \
                        mock.patch.object(
                            subject.os, "pidfd_open") as pidfd_open, \
                        self.assertRaises(subject.CampaignError):
                    subject.wait_sampler_ready(
                        FakeSampler(), pid_path,
                        base / subject.SAMPLER_NAME,
                        base / "thermal.csv", (1234, 99),
                        (700, 10, ["/unused"]), identities)
                pidfd_open.assert_not_called()
                self.assertEqual(identities, [])
            finally:
                os.close(descriptor)

    def test_sampler_pid_receipt_requires_exact_child_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = base / subject.SAMPLER_NAME
            csv_path = base / "thermal.csv"
            pid_path = base / "thermal.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            expected = subject.campaign_sampled_command(
                sampler, csv_path, pid_path)
            with mock.patch.object(
                    subject, "process_alive", return_value=True), \
                    mock.patch.object(
                        subject, "process_tokens", return_value=expected):
                self.assertEqual(
                    subject.sampler_identity(
                        pid_path, sampler, csv_path), 701)
            for changed in (
                    [*expected, "--extra"],
                    [
                        *expected[:expected.index("--pid-file") + 1],
                        str(base / "other.pid"),
                        *expected[expected.index("--pid-file") + 2:],
                    ],
            ):
                with self.subTest(changed=changed), \
                        mock.patch.object(
                            subject, "process_alive", return_value=True), \
                        mock.patch.object(
                            subject, "process_tokens", return_value=changed), \
                        self.assertRaises(subject.CampaignError):
                    subject.sampler_identity(
                        pid_path, sampler, csv_path)

    def test_sampler_identity_cmdline_race_to_exit_is_inert(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_path = base / "thermal.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            alive = iter((True, False))
            with mock.patch.object(
                    subject, "process_alive",
                    side_effect=lambda _pid: next(alive)), \
                    mock.patch.object(
                        subject, "process_tokens", return_value=[]):
                self.assertIsNone(subject.sampler_identity(
                    pid_path, base / subject.SAMPLER_NAME,
                    base / "thermal.csv"))

    def test_sampler_binding_rejects_reused_pid_before_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = base / subject.SAMPLER_NAME
            csv_path = base / "thermal.csv"
            pid_path = base / "thermal.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            controller_identity = (1234, 99)
            child_tokens = subject.campaign_sampled_command(
                sampler, csv_path, pid_path)
            supervisor_tokens = [
                str(subject.PYTHON_PATH), "-I", "-c",
                subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
                str(controller_identity[0]), str(controller_identity[1]),
                subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
                *child_tokens,
            ]
            descriptors = [
                os.open("/dev/null", os.O_RDONLY),
                os.open("/dev/null", os.O_RDONLY),
            ]
            child_start_reads = iter((20, 21))

            def start_ticks(pid):
                if pid == 701:
                    return next(child_start_reads)
                return 10

            def tokens(pid):
                return child_tokens if pid == 701 else supervisor_tokens

            try:
                with mock.patch.object(
                        subject.os, "pidfd_open",
                        side_effect=descriptors), \
                        mock.patch.object(
                            subject, "pidfd_has_exited",
                            return_value=False), \
                        mock.patch.object(
                            subject, "process_start_ticks",
                            side_effect=start_ticks), \
                        mock.patch.object(
                            subject, "process_parent_pid",
                            return_value=700), \
                        mock.patch.object(
                            subject, "process_tokens",
                            side_effect=tokens), \
                        self.assertRaises(subject.CampaignError):
                    subject.bind_sampler_identity(
                        pid_path, sampler, csv_path,
                        controller_identity,
                        (700, 10, supervisor_tokens))
            finally:
                for descriptor in descriptors:
                    try:
                        os.close(descriptor)
                    except OSError:
                        pass

    def test_sampler_binding_requires_exact_sudo_supervisor_ancestry(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = base / subject.SAMPLER_NAME
            csv_path = base / "thermal.csv"
            pid_path = base / "thermal.pid"
            pid_path.write_text("702\n", encoding="ascii")
            pid_path.chmod(0o444)
            controller_identity = (1234, 99)
            child_tokens = subject.campaign_sampled_command(
                sampler, csv_path, pid_path)
            supervisor_tokens = [
                str(subject.PYTHON_PATH), "-I", "-c",
                subject.PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
                str(controller_identity[0]), str(controller_identity[1]),
                subject.PDEATH_EXEC_CODE.encode("utf-8").hex(),
                *child_tokens,
            ]
            wrapper_tokens = [
                str(subject.SUDO_PATH), "-n", *supervisor_tokens,
            ]

            for supervisor_parent, accepted in ((700, True), (699, False)):
                descriptors = [
                    os.open("/dev/null", os.O_RDONLY)
                    for _ in range(3)
                ]
                starts = {700: 10, 701: 20, 702: 30}
                parents = {701: supervisor_parent, 702: 701}
                tokens = {
                    700: wrapper_tokens,
                    701: supervisor_tokens,
                    702: child_tokens,
                }
                try:
                    patches = (
                        mock.patch.object(
                            subject.os, "pidfd_open",
                            side_effect=descriptors),
                        mock.patch.object(
                            subject, "pidfd_has_exited",
                            return_value=False),
                        mock.patch.object(
                            subject, "process_start_ticks",
                            side_effect=lambda pid: starts.get(pid)),
                        mock.patch.object(
                            subject, "process_parent_pid",
                            side_effect=lambda pid: parents.get(pid)),
                        mock.patch.object(
                            subject, "process_tokens",
                            side_effect=lambda pid: tokens.get(pid)),
                    )
                    for patcher in patches:
                        patcher.start()
                    try:
                        if accepted:
                            self.assertEqual(
                                subject.bind_sampler_identity(
                                    pid_path, sampler, csv_path,
                                    controller_identity,
                                    (700, 10, wrapper_tokens)),
                                (702, 30, child_tokens))
                        else:
                            with self.assertRaises(subject.CampaignError):
                                subject.bind_sampler_identity(
                                    pid_path, sampler, csv_path,
                                    controller_identity,
                                    (700, 10, wrapper_tokens))
                    finally:
                        for patcher in reversed(patches):
                            patcher.stop()
                finally:
                    for descriptor in descriptors:
                        try:
                            os.close(descriptor)
                        except OSError:
                            pass

    def test_sampler_binding_accepts_direct_exec_and_closes_all_pidfds(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = base / subject.SAMPLER_NAME
            csv_path = base / "thermal.csv"
            pid_path = base / "thermal.pid"
            controller_identity = (1234, 99)
            child_tokens = subject.campaign_sampled_command(
                sampler, csv_path, pid_path)
            supervisor_tokens, wrapper_tokens = \
                subject.campaign_sampler_supervisor_commands(
                    sampler, csv_path, pid_path, controller_identity)

            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            direct_tokens = {
                700: supervisor_tokens,
                701: child_tokens,
            }
            direct_closes = []
            real_close = os.close

            def close_direct(descriptor):
                if descriptor in (900, 901):
                    direct_closes.append(descriptor)
                else:
                    real_close(descriptor)

            with mock.patch.object(
                    subject.os, "pidfd_open",
                    side_effect=(900, 901)), \
                    mock.patch.object(
                        subject, "pidfd_has_exited",
                        return_value=False), \
                    mock.patch.object(
                        subject, "process_start_ticks",
                        side_effect=lambda pid: {700: 10, 701: 30}.get(pid)), \
                    mock.patch.object(
                        subject, "process_parent_pid",
                        side_effect=lambda pid: {700: 1, 701: 700}.get(pid)), \
                    mock.patch.object(
                        subject, "process_tokens",
                        side_effect=lambda pid: direct_tokens.get(pid)), \
                    mock.patch.object(
                        subject.os, "close",
                        side_effect=close_direct):
                self.assertEqual(
                    subject.bind_sampler_identity(
                        pid_path, sampler, csv_path, controller_identity,
                        (700, 10, wrapper_tokens)),
                    (701, 30, child_tokens))
            self.assertEqual(direct_closes, [900, 901])

            pid_path.chmod(0o644)
            pid_path.write_text("702\n", encoding="ascii")
            pid_path.chmod(0o444)
            forked_tokens = {
                700: wrapper_tokens,
                701: supervisor_tokens,
                702: child_tokens,
            }
            forked_closes = []

            def fail_first_close(descriptor):
                if descriptor not in (910, 911, 912):
                    real_close(descriptor)
                    return
                forked_closes.append(descriptor)
                if len(forked_closes) == 1:
                    raise OSError("injected close failure")

            with mock.patch.object(
                    subject.os, "pidfd_open",
                    side_effect=(910, 911, 912)), \
                    mock.patch.object(
                        subject, "pidfd_has_exited",
                        return_value=False), \
                    mock.patch.object(
                        subject, "process_start_ticks",
                        side_effect=lambda pid:
                            {700: 10, 701: 20, 702: 30}.get(pid)), \
                    mock.patch.object(
                        subject, "process_parent_pid",
                        side_effect=lambda pid:
                            {700: 1, 701: 700, 702: 701}.get(pid)), \
                    mock.patch.object(
                        subject, "process_tokens",
                        side_effect=lambda pid: forked_tokens.get(pid)), \
                    mock.patch.object(
                        subject.os, "close",
                        side_effect=fail_first_close), \
                    self.assertRaises(subject.CampaignError):
                subject.bind_sampler_identity(
                    pid_path, sampler, csv_path, controller_identity,
                    (700, 10, wrapper_tokens))
            self.assertEqual(forked_closes, [912, 911, 910])

    def test_benchmark_direct_and_parent_death_commands_match(self):
        expected = ["/campaign/wirehair_v2_bench", "precodefail"]
        self.assertTrue(subject.benchmark_command_matches(
            expected, expected))
        self.assertFalse(subject.benchmark_command_matches(
            expected, expected, (1234, 99)))
        self.assertTrue(subject.benchmark_command_matches([
            str(subject.PYTHON_PATH), "-I", "-c",
            subject.PDEATH_EXEC_CODE, "1234", "99", *expected,
        ], expected, (1234, 99)))
        self.assertFalse(subject.benchmark_command_matches([
            str(subject.PYTHON_PATH), "-I", "-c",
            subject.PDEATH_EXEC_CODE, "1234", "100", *expected,
        ], expected, (1234, 99)))
        self.assertFalse(subject.benchmark_command_matches([
            str(subject.PYTHON_PATH), "-I", "-c",
            subject.PDEATH_EXEC_CODE, "1234", "99",
            "/other/wirehair_v2_bench", "precodefail",
        ], expected, (1234, 99)))

    def test_pre_ready_sampler_is_reconciled_by_exact_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            pid_path = root / "thermal" / "segment000.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            state = {"alive": True}

            def terminate_discovered(matchers):
                labels = [label for label, _ in matchers]
                if labels == ["thermal"]:
                    state["alive"] = False
                    return [{
                        "kind": "thermal", "pid": 701,
                        "action": "terminated",
                    }]
                return []

            with mock.patch.object(
                    subject, "process_alive",
                    side_effect=lambda pid: state["alive"]), \
                    mock.patch.object(
                        subject, "process_start_ticks", return_value=20), \
                    mock.patch.object(
                        subject, "terminate_discovered_owned_processes",
                        side_effect=terminate_discovered), \
                    mock.patch.object(
                        subject, "discover_owned_processes", return_value=[]):
                actions = subject.terminate_owned_segment_processes(
                    root, 0, [], (1234, 99))
            self.assertEqual(actions, [{
                "kind": "thermal", "pid": 701, "action": "terminated",
            }])
            self.assertFalse(pid_path.exists())

    def test_incomplete_sampler_publication_is_removed_only_after_no_reader(self):
        for blocker in ("none", "direct", "inventory", "wrapper"):
            with self.subTest(blocker=blocker), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "campaign"
                subject.make_private_directory(root)
                for name in subject.RUNTIME_DIRECTORY_NAMES:
                    subject.make_private_directory(root / name)
                subject.make_private_directory(root / "frozen")
                pid_path = root / "thermal" / "segment000.pid"
                descriptor = os.open(
                    str(pid_path),
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                    0o444)
                os.fchmod(descriptor, 0o444)
                os.write(descriptor, b"70")
                os.fsync(descriptor)
                os.close(descriptor)
                direct = [(
                    "unreceipted_thermal", 701, 20,
                    subject.campaign_sampled_command(
                        root / "frozen" / subject.SAMPLER_NAME,
                        root / "thermal" / "segment000.csv", pid_path),
                )] if blocker == "direct" else []
                discoveries = iter((direct, []))

                def terminate_wrappers(_matchers):
                    if blocker == "wrapper" and \
                            terminate_wrappers.calls == 0:
                        terminate_wrappers.calls += 1
                        raise subject.CampaignError("injected wrapper failure")
                    terminate_wrappers.calls += 1
                    return []

                terminate_wrappers.calls = 0
                with mock.patch.object(
                        subject, "terminate_discovered_owned_processes",
                        side_effect=terminate_wrappers), \
                        mock.patch.object(
                            subject, "discover_owned_processes",
                            side_effect=lambda _matchers:
                                next(discoveries)), \
                        mock.patch.object(
                            subject, "other_samplers",
                            return_value=(
                                [{"pid": 702, "command": "foreign",
                                  "i2c_devices": ["/dev/i2c-1"]}]
                                if blocker == "inventory" else [])):
                    if blocker != "none":
                        with self.assertRaises(subject.CampaignError):
                            subject.terminate_owned_segment_processes(
                                root, 0, [], (1234, 99))
                    else:
                        self.assertEqual(
                            subject.terminate_owned_segment_processes(
                                root, 0, [], (1234, 99)),
                            [])
                self.assertEqual(pid_path.exists(), blocker != "none")

    def test_reused_unbound_complete_sampler_receipt_is_never_signaled(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            pid_path = root / "thermal" / "segment000.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            start = mock.Mock()
            terminate = mock.Mock()
            with mock.patch.object(
                    subject, "terminate_discovered_owned_processes",
                    return_value=[]), \
                    mock.patch.object(
                        subject, "discover_owned_processes",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "other_samplers", return_value=[]), \
                    mock.patch.object(
                        subject, "process_alive", return_value=True), \
                    mock.patch.object(
                        subject, "process_start_ticks", start), \
                    mock.patch.object(
                        subject, "terminate_verified_process", terminate):
                self.assertEqual(
                    subject.terminate_owned_segment_processes(
                        root, 0, [], (1234, 99)),
                    [])
            start.assert_not_called()
            terminate.assert_not_called()
            self.assertFalse(pid_path.exists())

    def test_receipted_sampler_post_stop_probe_race_to_exit_is_inert(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            pid_path = root / "thermal" / "segment000.pid"
            pid_path.write_text("701\n", encoding="ascii")
            pid_path.chmod(0o444)
            subject.segment_path(root, 0, "ready").write_bytes(
                subject.canonical_json({
                    "sampler_pid": 701, "sampler_start_ticks": 20,
                }))
            alive = iter((True, True, False))
            with mock.patch.object(
                    subject, "terminate_discovered_owned_processes",
                    return_value=[]), \
                    mock.patch.object(
                        subject, "discover_owned_processes",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "other_samplers", return_value=[]), \
                    mock.patch.object(
                        subject, "process_alive",
                        side_effect=lambda _pid: next(alive)), \
                    mock.patch.object(
                        subject, "process_start_ticks", return_value=None), \
                    mock.patch.object(
                        subject, "terminate_verified_process",
                        return_value="terminated"):
                actions = subject.terminate_owned_segment_processes(
                    root, 0, [], (1234, 99))
            self.assertEqual(actions, [{
                "kind": "thermal", "pid": 701, "action": "terminated",
            }])
            self.assertFalse(pid_path.exists())

    def test_segment_cleanup_attempts_other_categories_after_corruption(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "campaign"
            subject.make_private_directory(root)
            for name in subject.RUNTIME_DIRECTORY_NAMES:
                subject.make_private_directory(root / name)
            subject.make_private_directory(root / "frozen")
            attempt = root / "attempts" / "segment000"
            subject.make_private_directory(attempt)
            thermal_pid = root / "thermal" / "segment000.pid"
            thermal_pid.write_text("malformed\n", encoding="ascii")
            thermal_pid.chmod(0o444)
            (attempt / "job00000.pid").write_text(
                "not json\n", encoding="ascii")
            (attempt / "job00001.pid").write_bytes(subject.canonical_json({
                "pid": 702, "start_ticks": 20,
            }))
            tasks = [
                {"job": 0, "argv": ["/bin/false"]},
                {"job": 1, "argv": ["/bin/true"]},
            ]
            discoveries = []
            terminated = []

            def discover_wrappers(matchers):
                discoveries.append([label for label, _ in matchers])
                if len(discoveries) == 1:
                    raise subject.CampaignError(
                        "injected thermal discovery failure")
                return []

            def terminate(pid, start_ticks, tokens):
                terminated.append((pid, start_ticks, tokens))
                return "terminated"

            with mock.patch.object(
                    subject, "terminate_discovered_owned_processes",
                    side_effect=discover_wrappers), \
                    mock.patch.object(
                        subject, "discover_owned_processes",
                        return_value=[]), \
                    mock.patch.object(
                        subject, "process_alive",
                        side_effect=lambda pid: pid == 702), \
                    mock.patch.object(
                        subject, "process_start_ticks",
                        side_effect=lambda pid: 20 if pid == 702 else None), \
                    mock.patch.object(
                        subject, "process_tokens",
                        side_effect=lambda pid: ["/bin/true"]), \
                    mock.patch.object(
                        subject, "terminate_verified_process",
                        side_effect=terminate), \
                    self.assertRaises(subject.CampaignError):
                subject.terminate_owned_segment_processes(
                    root, 0, tasks, (1234, 99))
            self.assertEqual(discoveries, [
                ["thermal"],
                ["benchmark_job_0", "benchmark_job_1"],
            ])
            self.assertEqual(terminated, [(702, 20, ["/bin/true"])])

    def test_live_receipt_with_unreadable_start_fails_closed(self):
        for kind in ("thermal", "benchmark"):
            with self.subTest(kind=kind), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "campaign"
                subject.make_private_directory(root)
                for name in subject.RUNTIME_DIRECTORY_NAMES:
                    subject.make_private_directory(root / name)
                subject.make_private_directory(root / "frozen")
                tasks = []
                protected_path = (
                    root / "thermal" / "segment000.pid")
                if kind == "thermal":
                    protected_path.write_text("701\n", encoding="ascii")
                    protected_path.chmod(0o444)
                    subject.segment_path(root, 0, "ready").write_bytes(
                        subject.canonical_json({
                            "sampler_pid": 701,
                            "sampler_start_ticks": 20,
                        }))
                else:
                    attempt = root / "attempts" / "segment000"
                    subject.make_private_directory(attempt)
                    protected_path = attempt / "job00000.pid"
                    protected_path.write_bytes(subject.canonical_json({
                        "pid": 701, "start_ticks": 20,
                    }))
                    tasks = [{"job": 0, "argv": ["/bin/true"]}]
                with mock.patch.object(
                        subject, "terminate_discovered_owned_processes",
                        return_value=[]), \
                        mock.patch.object(
                            subject, "discover_owned_processes",
                            return_value=[]), \
                        mock.patch.object(
                            subject, "process_alive", return_value=True), \
                        mock.patch.object(
                            subject, "process_start_ticks",
                            return_value=None), \
                        self.assertRaises(subject.CampaignError):
                    subject.terminate_owned_segment_processes(
                        root, 0, tasks, (1234, 99))
                self.assertTrue(protected_path.exists())

    def test_pid_reuse_after_term_never_gets_kill(self):
        state = {"start_ticks": 99}
        signals = []

        def send(signal_name):
            signals.append(signal_name)
            if signal_name == "TERM":
                state["start_ticks"] = 100

        action = subject.terminate_verified_process(
            1234, 99, ["owned", "worker"],
            alive_probe=lambda pid: True,
            start_ticks_probe=lambda pid: state["start_ticks"],
            tokens_probe=lambda pid: ["owned", "worker"],
            signal_sender=send,
        )
        self.assertEqual(action, "identity_reused")
        self.assertEqual(signals, ["TERM"])

    def test_signal_helper_race_to_natural_exit_is_stopped(self):
        alive = iter((True, True, False))
        action = subject.terminate_verified_process(
            1234, 99, ["owned", "worker"],
            alive_probe=lambda pid: next(alive),
            start_ticks_probe=lambda pid: 99,
            tokens_probe=lambda pid: ["owned", "worker"],
            signal_sender=lambda signal_name: False,
        )
        self.assertEqual(action, "terminated")

    def test_identity_probe_race_to_exit_is_not_misclassified_as_corruption(self):
        for field in ("start", "tokens_none", "tokens_empty"):
            with self.subTest(field=field):
                alive_values = iter((True, False))
                sender = mock.Mock()
                action = subject.terminate_verified_process(
                    1234, 99, ["owned", "worker"],
                    alive_probe=lambda _pid: next(alive_values),
                    start_ticks_probe=(
                        (lambda _pid: None) if field == "start" else
                        (lambda _pid: 99)),
                    tokens_probe=(
                        (lambda _pid: [])
                        if field == "tokens_empty" else
                        (lambda _pid: None)),
                    signal_sender=sender,
                )
                self.assertEqual(action, "already_exited")
                sender.assert_not_called()

    def test_unreadable_identity_that_remains_live_still_fails_closed(self):
        for field in ("start", "tokens"):
            with self.subTest(field=field), \
                    self.assertRaises(subject.CampaignError):
                subject.terminate_verified_process(
                    1234, 99, ["owned", "worker"],
                    alive_probe=lambda _pid: True,
                    start_ticks_probe=(
                        (lambda _pid: None) if field == "start" else
                        (lambda _pid: 99)),
                    tokens_probe=lambda _pid: None,
                    signal_sender=lambda _signal: True,
                )

    def test_sampler_stop_rejects_identity_loss(self):
        class FakeSampler:
            pid = 700

            @staticmethod
            def wait(timeout):
                return 0

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_file = base / "sampler.pid"
            sampler_path = base / "sampler.py"
            csv_path = base / "sampler.csv"
            identity = (701, 99, ["python3", str(sampler_path)])
            for action in ("identity_reused", "identity_changed"):
                with self.subTest(action=action), \
                        mock.patch.object(
                            subject, "terminate_verified_process",
                            return_value=action), \
                        mock.patch.object(subject, "other_samplers", return_value=[]), \
                        self.assertRaises(subject.CampaignError):
                    subject.stop_launched_sampler(
                        FakeSampler(), pid_file, sampler_path, csv_path,
                        identity)

    def test_sampler_stop_exception_still_waits_and_checks_readers(self):
        class FakeSampler:
            pid = 700

            def __init__(self):
                self.wait_calls = []

            def wait(self, timeout):
                self.wait_calls.append(timeout)
                return 0

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = FakeSampler()
            check = mock.Mock(return_value=[])
            pid_file = base / "sampler.pid"
            pid_file.write_text("701\n", encoding="ascii")
            with mock.patch.object(
                    subject, "terminate_verified_process",
                    side_effect=subject.CampaignError("injected")), \
                    mock.patch.object(subject, "other_samplers", check), \
                    self.assertRaises(subject.CampaignError):
                subject.stop_launched_sampler(
                    sampler, pid_file, base / "sampler.py",
                    base / "sampler.csv",
                    (701, 99, ["python3", str(base / "sampler.py")]))
            self.assertEqual(sampler.wait_calls, [15.0])
            check.assert_called_once_with()
            self.assertTrue(pid_file.exists())

    def test_sampler_stop_keeps_pid_when_reader_check_fails(self):
        class FakeSampler:
            pid = 700

            @staticmethod
            def wait(timeout):
                return 0

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pid_file = base / "sampler.pid"
            pid_file.write_text("701\n", encoding="ascii")
            with mock.patch.object(
                    subject, "terminate_verified_process",
                    return_value="terminated"), \
                    mock.patch.object(
                        subject, "other_samplers",
                        side_effect=subject.CampaignError("injected")), \
                    self.assertRaises(subject.CampaignError):
                subject.stop_launched_sampler(
                    FakeSampler(), pid_file, base / "sampler.py",
                    base / "sampler.csv",
                    (701, 99, ["python3", str(base / "sampler.py")]))
            self.assertTrue(pid_file.exists())

    def test_sampler_stop_exception_still_stops_timed_out_wrapper(self):
        class FakeSampler:
            pid = 700

            def __init__(self):
                self.wait_calls = []

            def wait(self, timeout):
                self.wait_calls.append(timeout)
                if len(self.wait_calls) == 1:
                    raise subject.subprocess.TimeoutExpired("wait", timeout)
                return -15

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            sampler = FakeSampler()
            wrapper_stop = mock.Mock(return_value="terminated")
            check = mock.Mock(return_value=[])
            with mock.patch.object(
                    subject, "terminate_verified_process",
                    side_effect=subject.CampaignError("injected")), \
                    mock.patch.object(
                        subject, "terminate_privileged_direct_unreaped_child",
                        wrapper_stop), \
                    mock.patch.object(subject, "other_samplers", check), \
                    self.assertRaises(subject.CampaignError):
                subject.stop_launched_sampler(
                    sampler, base / "sampler.pid", base / "sampler.py",
                    base / "sampler.csv",
                    (701, 99, ["python3", str(base / "sampler.py")]))
            self.assertEqual(sampler.wait_calls, [15.0, 5.0])
            wrapper_stop.assert_called_once_with(sampler)
            check.assert_called_once_with()

    def test_failed_binding_cleanup_checks_readers_after_stop_exception(self):
        class FakeSampler:
            pid = 700

        sampler = FakeSampler()
        wrapper_stop = mock.Mock(
            side_effect=subject.CampaignError("injected"))
        check = mock.Mock(return_value=[])
        with mock.patch.object(
                subject, "terminate_privileged_direct_unreaped_child",
                wrapper_stop), \
                mock.patch.object(subject, "other_samplers", check), \
                self.assertRaises(subject.CampaignError):
            subject.cleanup_failed_sampler_binding(sampler)
        wrapper_stop.assert_called_once_with(sampler)
        check.assert_called_once_with()

    def test_cleanup_action_must_prove_stopped(self):
        for action in subject.STOPPED_PROCESS_ACTIONS:
            subject.require_stopped_action(action, "test")
        for action in ("identity_reused", "identity_changed",
                       "identity_replaced_before_term"):
            with self.subTest(action=action), self.assertRaises(
                    subject.CampaignError):
                subject.require_stopped_action(action, "test")


class SelftestCommandTest(unittest.TestCase):
    def test_selftest_command_passes(self):
        subject.command_selftest(None)


if __name__ == "__main__":
    unittest.main()
