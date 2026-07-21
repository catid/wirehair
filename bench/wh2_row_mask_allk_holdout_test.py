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
import shutil
import sys
import tempfile
import unittest
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
        self.assertTrue(p48["acceptance"]["accepted"])
        self.assertTrue(p32["acceptance"]["accepted"])
        # h1 is repaired but is selection material, so it is a repair in the
        # raw ledger and absent from the headline ledger.
        self.assertEqual(p48["repairs_vs_comparator"], 1)
        self.assertEqual(p48["headline_repairs"], 0)
        self.assertEqual(p48["introductions_vs_comparator"], 0)
        self.assertEqual(p32["repairs_vs_comparator"], 2)
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
        self.assertEqual(arm["weak_K_multiplicity"], {"3": 1, "5": 1, "11": 1})
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
        outcomes[("p48_r3_pfx007", "control", 3, 0, "burst")] = \
            subject.synth_cell("field_shortfall")
        outcomes[("p48_r3_pfx007", "hard", 3, 0, "burst")] = \
            subject.synth_cell("field_shortfall")
        reduction = self.reduce(outcomes)
        self.assertFalse(reduction["accepted"])
        self.assertFalse(
            reduction["candidates"]["p48_r3_pfx007"]
            ["acceptance"]["gate_field_shortfall_improved"])

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
            with self.assertRaises(subject.CampaignError):
                subject.validate_thermal_rows([], segment=0)

    def test_interrupted_segment_reconciliation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "interrupted"
            for directory in (
                    "raw", "stderr", "exit", "thermal", "segments",
                    "attempts", "job_receipts"):
                (root / directory).mkdir(parents=True, exist_ok=True)
            task = {"job": 0, "stage": "hard", "output_name": "retry.csv",
                    "argv": ["/bin/true"]}
            intent = {
                "schema": subject.SCHEMA + ".segment_intent", "segment": 0,
                "stage": "hard", "jobs": [0],
                "jobs_sha256": subject.sha256_bytes(subject.json_lines([0])),
                "workers": 1,
                "retry_policy": "stage-atomic non-selective retry",
            }
            subject.write_once(
                subject.segment_path(root, 0, "intent"),
                subject.canonical_json(intent))
            (root / "attempts" / "segment000").mkdir()
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
        with self.assertRaises(subject.CampaignError):
            self.load()
        self.marker_path.write_bytes(
            subject.canonical_json({"schema": "other.schema.v1"}))
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


class SelftestCommandTest(unittest.TestCase):
    def test_selftest_command_passes(self):
        subject.command_selftest(None)


if __name__ == "__main__":
    unittest.main()
