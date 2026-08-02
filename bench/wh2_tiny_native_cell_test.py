#!/usr/bin/env python3
"""Bounded CLI/identity smoke tests for the native tiny-WH2 cell worker."""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import unittest


SCHEMA = "wirehair.wh2.tiny-native-cell.v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class TinyNativeCellTest(unittest.TestCase):
    worker = None

    @classmethod
    def setUpClass(cls):
        if not hasattr(os, "sched_getaffinity"):
            raise unittest.SkipTest("Linux CPU affinity is required")
        affinity = os.sched_getaffinity(0)
        if not affinity:
            raise unittest.SkipTest("no allowed CPU is available")
        cls.cpu = min(affinity)

    def command(self, **changes):
        values = {
            "k": 5,
            "block_bytes": 64,
            "layout": "disabled",
            "scope": "solve",
            "left": "off",
            "right": "on",
            "construction_attempt": 3,
            "trace_seed": 1,
            "replicate": 0,
            "cpu": self.cpu,
            "invocations_per_slot": 2,
        }
        values.update(changes)
        return [
            str(self.worker),
            "--k", str(values["k"]),
            "--block-bytes", str(values["block_bytes"]),
            "--layout", str(values["layout"]),
            "--scope", str(values["scope"]),
            "--left", str(values["left"]),
            "--right", str(values["right"]),
            "--construction-attempt", str(values["construction_attempt"]),
            "--trace-seed", str(values["trace_seed"]),
            "--replicate", str(values["replicate"]),
            "--cpu", str(values["cpu"]),
            "--invocations-per-slot", str(values["invocations_per_slot"]),
        ]

    def run_raw(self, command):
        return subprocess.run(
            command, text=True, capture_output=True, check=False)

    def run_cell(self, **changes):
        result = self.run_raw(self.command(**changes))
        self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
        self.assertEqual(result.stderr, "")
        self.assertTrue(result.stdout.endswith("\n"))
        record = json.loads(result.stdout)
        canonical = json.dumps(record, sort_keys=True, separators=(",", ":"))
        self.assertEqual(result.stdout.rstrip("\n"), canonical)
        self.assert_record_shape(record)
        return record

    def assert_record_shape(self, record):
        self.assertEqual(set(record), {
            "arms", "panel", "request", "request_sha256", "schema",
            "source_sha256", "trace",
        })
        self.assertEqual(record["schema"], SCHEMA)
        self.assertRegex(record["request_sha256"], SHA256)
        self.assertRegex(record["source_sha256"], SHA256)
        self.assertRegex(record["trace"]["sha256"], SHA256)
        request_json = json.dumps(
            record["request"], sort_keys=True, separators=(",", ":"))
        self.assertEqual(
            record["request_sha256"],
            hashlib.sha256(request_json.encode("utf-8")).hexdigest())
        self.assertEqual(record["trace"]["packet_count"],
                         record["request"]["K"] + 256)
        self.assertEqual(record["trace"]["loss_seed"],
                         record["request"]["trace_seed"])
        self.assertEqual(record["trace"]["requested_seed"],
                         record["request"]["trace_seed"])
        self.assertEqual(len(record["panel"]["slots"]), 8)
        self.assertEqual(
            [slot["index"] for slot in record["panel"]["slots"]],
            list(range(8)))
        for arm in record["arms"].values():
            self.assertRegex(arm["equations_sha256"], SHA256)
            self.assertRegex(arm["construction_sha256"], SHA256)
            for summary in arm["direct"].values():
                self.assertEqual(summary["calls"], 5)
                self.assertTrue(summary["identity_ok"])

    def assert_zero_direct(self, arm):
        for summary in arm["direct"].values():
            self.assertEqual(summary["attempts_min"], 0)
            self.assertEqual(summary["attempts_max"], 0)
            self.assertEqual(summary["attempts_total"], 0)
            self.assertEqual(summary["completions_total"], 0)
            self.assertEqual(summary["fallbacks_total"], 0)

    def assert_complete_outcome(self, record, comparable):
        panel = record["panel"]
        self.assertEqual(panel["status"], "complete")
        self.assertEqual(panel["diagnostic"], "")
        self.assertEqual(panel["comparable"], comparable)
        self.assertIsNotNone(panel["left_preflight"])
        self.assertIsNotNone(panel["right_preflight"])
        for slot in panel["slots"]:
            self.assertIn(slot["side"], ("left", "right"))
            if comparable:
                self.assertEqual(slot["disposition"], "success")
                self.assertGreater(slot["elapsed_ns"], 0)
            else:
                self.assertEqual(slot["disposition"], "preflight_failure")
                self.assertIsNone(slot["elapsed_ns"])

    def test_parser_rejects_malformed_and_forbidden_requests(self):
        help_result = self.run_raw([str(self.worker), "--help"])
        self.assertEqual(help_result.returncode, 0)
        self.assertIn("--scope <solve|receive|encoder>", help_result.stdout)
        self.assertEqual(help_result.stderr, "")

        cases = []
        cases.append(self.command(k=9))
        cases.append(self.command(invocations_per_slot=32769))
        cases.append(self.command(left="wh1"))
        cases.append(self.command()[:-2])
        duplicate = self.command()
        duplicate[5] = "--k"
        cases.append(duplicate)
        leading_zero = self.command()
        leading_zero[2] = "05"
        cases.append(leading_zero)
        for command in cases:
            with self.subTest(command=command):
                result = self.run_raw(command)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(result.stdout, "")
                self.assertIn("wirehair_wh2_tiny_native_cell:", result.stderr)
                self.assertIn("usage:", result.stderr)

    def test_requested_weak_attempt_is_not_retried(self):
        record = self.run_cell(construction_attempt=0)
        self.assertEqual(record["request"]["construction_attempt"], 0)
        self.assert_complete_outcome(record, False)
        self.assertEqual(record["arms"]["left"]["equations_sha256"],
                         record["arms"]["right"]["equations_sha256"])
        self.assertEqual(record["arms"]["left"]["construction_sha256"],
                         record["arms"]["right"]["construction_sha256"])
        for side in ("left_preflight", "right_preflight"):
            self.assertEqual(record["panel"][side]["outcome"],
                             "bad_peel_seed")
        for slot in record["panel"]["slots"]:
            self.assertEqual(slot["outcome"], "bad_peel_seed")
        self.assert_zero_direct(record["arms"]["left"])
        direct = record["arms"]["right"]["direct"]
        self.assertEqual(direct["construction"]["attempts_min"], 1)
        self.assertEqual(direct["construction"]["attempts_max"], 1)
        self.assertEqual(direct["construction"]["attempts_total"], 5)
        self.assertEqual(direct["construction"]["completions_total"], 5)
        self.assertEqual(direct["operation"]["attempts_total"], 0)

    def test_known_attempt_three_succeeds_and_seed_is_not_requalified(self):
        even = self.run_cell(construction_attempt=3, replicate=0)
        odd = self.run_cell(construction_attempt=3, replicate=1)
        for record, order in ((even, "ABBA"), (odd, "BAAB")):
            self.assertEqual(record["request"]["construction_attempt"], 3)
            self.assertEqual(record["panel"]["order"], order)
            self.assert_complete_outcome(record, True)
            self.assertEqual(record["arms"]["left"]["equations_sha256"],
                             record["arms"]["right"]["equations_sha256"])
            self.assertEqual(record["arms"]["left"]["construction_sha256"],
                             record["arms"]["right"]["construction_sha256"])
            self.assert_zero_direct(record["arms"]["left"])
            direct = record["arms"]["right"]["direct"]
            for phase in ("construction", "operation"):
                self.assertEqual(direct[phase]["attempts_min"], 1)
                self.assertEqual(direct[phase]["attempts_max"], 1)
                self.assertEqual(direct[phase]["attempts_total"], 5)
                self.assertEqual(direct[phase]["completions_total"], 5)
                self.assertEqual(direct[phase]["fallbacks_total"], 0)
        self.assertEqual(even["source_sha256"], odd["source_sha256"])
        self.assertEqual(even["trace"]["sha256"], odd["trace"]["sha256"])
        self.assertEqual(even["trace"]["derived_seed"],
                         odd["trace"]["derived_seed"])

    def test_attempt_zero_all_k_and_layouts_has_exact_direct_identity(self):
        for k in range(2, 9):
            for layout in ("disabled", "two07"):
                with self.subTest(k=k, layout=layout):
                    record = self.run_cell(
                        k=k, layout=layout, construction_attempt=0)
                    self.assertEqual(record["request"]["construction_attempt"],
                                     0)
                    panel = record["panel"]
                    self.assertEqual(panel["status"], "complete")
                    self.assertEqual(panel["left_preflight"]["outcome"],
                                     panel["right_preflight"]["outcome"])
                    self.assertEqual(
                        record["arms"]["left"]["equations_sha256"],
                        record["arms"]["right"]["equations_sha256"])
                    self.assertEqual(
                        record["arms"]["left"]["construction_sha256"],
                        record["arms"]["right"]["construction_sha256"])
                    self.assert_zero_direct(record["arms"]["left"])

    def test_encoder_smoke_k5_k8_both_layouts(self):
        attempts = {
            (5, "disabled"): 3,
            (5, "two07"): 2,
            (8, "disabled"): 0,
            (8, "two07"): 0,
        }
        for (k, layout), attempt in attempts.items():
            with self.subTest(k=k, layout=layout):
                record = self.run_cell(
                    k=k, layout=layout, scope="encoder", left="on",
                    right="wh1", construction_attempt=attempt, replicate=1)
                self.assert_complete_outcome(record, True)
                self.assertEqual(record["arms"]["left"]["mode"], "on")
                self.assertEqual(record["arms"]["right"]["mode"], "wh1")
                self.assert_zero_direct(record["arms"]["right"])
                construction = record["arms"]["left"]["direct"][
                    "construction"]
                operation = record["arms"]["left"]["direct"]["operation"]
                self.assertEqual(construction["attempts_total"], 0)
                self.assertGreater(operation["attempts_total"], 0)
                self.assertEqual(operation["attempts_total"],
                                 operation["completions_total"])
                self.assertEqual(operation["fallbacks_total"], 0)
                for slot in record["panel"]["slots"]:
                    self.assertIsNone(slot["decoded_overhead"])

    def test_receive_smoke_k5_k8_both_layouts(self):
        attempts = {
            (5, "disabled"): 3,
            (5, "two07"): 2,
            (8, "disabled"): 0,
            (8, "two07"): 0,
        }
        for (k, layout), attempt in attempts.items():
            with self.subTest(k=k, layout=layout):
                record = self.run_cell(
                    k=k, layout=layout, scope="receive", left="off",
                    right="on", construction_attempt=attempt, replicate=1)
                self.assert_complete_outcome(record, True)
                self.assert_zero_direct(record["arms"]["left"])
                for side in ("left_preflight", "right_preflight"):
                    self.assertEqual(record["panel"][side][
                        "decoded_overhead"], 0)
                direct = record["arms"]["right"]["direct"]
                self.assertGreater(direct["construction"]["attempts_total"], 0)
                self.assertGreater(direct["operation"]["attempts_total"], 0)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: wh2_tiny_native_cell_test.py WORKER")
    worker = Path(sys.argv[1]).resolve()
    if not worker.is_file():
        raise SystemExit("worker executable does not exist: {}".format(worker))
    TinyNativeCellTest.worker = worker
    unittest.main(argv=[sys.argv[0]], verbosity=2)


if __name__ == "__main__":
    main()
