#!/usr/bin/env python3
"""Neutral synthetic reducer/controller tests; no actual mapping or timing run."""
import copy
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location("mulrows_r0_test", Path(__file__).with_name("Wh2ThueMorseMulRowsR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def identity():
    canonical = b"neutral synthetic identity receipt"
    derived = dict.fromkeys(M.R.H.TARGET_DERIVED_KEYS, 0)
    derived.update(ccd_id=6, complex_id=6, core_id=50, family=26, full_apic_id=100,
                   initial_apic_id_8=100, logical_processors_per_package=128, model=8,
                   package_id=0, stepping=1, thread_id=0, threads_per_core=2)
    return dict(requested_cpu=50, before_cpu=50, after_cpu=50,
                raw_capture_complete=True, semantic_validation_passed=True,
                singleton_affinity_verified=True, canonical_bytes=len(canonical),
                canonical_hex=canonical.hex(), canonical_sha256=M.C.sha(canonical),
                capabilities=dict.fromkeys(M.R.H.TARGET_CAPABILITY_KEYS, 0), derived=derived)


class ScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lookup = bytes((i * 17 + i // 37) & 255 for i in range(39936))
        rows = M.expected_rows(cls.lookup)
        cls.base = dict(schema=M.PROTOCOL + ".raw.v1", protocol=M.PROTOCOL,
                        outcome="COMPLETE", lookup_sha256=M.LOOKUP_SHA,
                        warm_callbacks=576, measured_callbacks=2304, row_calls=11796480,
                        target_cpu=50, calls_per_callback=4096, preflight_row_calls=256,
                        checked_rows=184320, output_stride=38, output_address=65536,
                        baseline_function_address=4096, candidate_function_address=8192,
                        expected_rows_hex=rows.hex(), expected_rows_sha256=M.C.sha(rows),
                        checksum=M.output_checksum(rows), elapsed_ns=10000000,
                        identity_before=identity(), identity_after=identity(),
                        runtime_before=dict(polynomial=333, address=1048576, ssse3=1, avx2=1, gfni=1, avx512=1),
                        runtime_after=dict(polynomial=333, address=1048576, ssse3=1, avx2=1, gfni=1, avx512=1),
                        panels=[])
        for coord in M.roster():
            _, _, comparison, order = coord
            sides = [0, 1, 1, 0, 1, 0, 0, 1]
            times = [900 if comparison == 2 and (side ^ order) == 1 else 1000 for side in sides]
            cls.base["panels"].append(list(coord) + [times])

    def reduce(self, raw):
        return M.reduce_raw(raw, self.lookup)

    def test_complete_pass_and_direction(self):
        result = self.reduce(self.base)
        self.assertEqual(result["outcome"], "PASS")
        self.assertEqual(len(result["statistics"]), 24)
        for item in result["statistics"]:
            self.assertAlmostEqual(item["estimate"]["geometric_mean"], .9 if item["comparison"] == 2 else 1)
        self.assertFalse(result["WH1_compared"])
        self.assertFalse(result["production_promotion_claimed"])

    def test_control_failure_cannot_be_pooled(self):
        raw = copy.deepcopy(self.base)
        for row in raw["panels"]:
            if row[2] == 0:
                row[4] = [1050, 1000, 1000, 1050, 1000, 1050, 1050, 1000]
        result = self.reduce(raw)
        self.assertEqual(result["outcome"], "CONTROL_FAIL")
        self.assertEqual(len(result["failed_controls"]), 8)

    def test_no_gain_and_local_regression_fail(self):
        for family, value in ((3, 1000), (0, 1030), (1, 1030), (2, 1030)):
            raw = copy.deepcopy(self.base)
            for row in raw["panels"]:
                if row[1] == family and row[2] == 2:
                    sides = [0, 1, 1, 0, 1, 0, 0, 1]
                    row[4] = [value if (side ^ row[3]) == 1 else 1000 for side in sides]
            with self.subTest(family=family):
                self.assertEqual(self.reduce(raw)["outcome"], "SHORT_FAIL")

    def test_one_percent_gain_is_acceptable(self):
        raw = copy.deepcopy(self.base)
        for row in raw["panels"]:
            if row[2] == 2:
                row[4] = [990 if value == 900 else value for value in row[4]]
        self.assertEqual(self.reduce(raw)["outcome"], "PASS")

    def test_roster_and_status_forgeries(self):
        mutations = [(["outcome"], "PASS"), (["panels", 0, 0], True),
                     (["panels", 0, 4, 0], 0), (["panels", 0, 4, 0], True),
                     (["panels", 0, 4, 0], 10000001), (["elapsed_ns"], 1),
                     (["row_calls"], 1), (["checked_rows"], 1), (["warm_callbacks"], False),
                     (["checksum"], 0), (["expected_rows_sha256"], "0" * 64),
                     (["expected_rows_hex"], "00" * 1536), (["output_stride"], 6),
                     (["output_address"], (1 << 64) - 1), (["candidate_function_address"], 4096),
                     (["runtime_after", "address"], 42), (["identity_before", "derived", "full_apic_id"], 50),
                     (["identity_before", "canonical_sha256"], "0" * 64)]
        for path, replacement in mutations:
            raw = copy.deepcopy(self.base)
            target = raw
            for field in path[:-1]:
                target = target[field]
            target[path[-1]] = replacement
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.reduce(raw)
        for raw in (dict(self.base, extra=1), dict(self.base, panels=self.base["panels"][:-1])):
            with self.assertRaises(ValueError):
                self.reduce(raw)

    def test_family_depth_and_slots(self):
        for family in range(4):
            values = M.ids(family)
            self.assertEqual(len(set(values)), 64)
            for value in values:
                self.assertEqual(sum(bool(x) for x in (value >> 24, (value >> 17) & 127,
                                                       (value >> 10) & 127)), family)
        self.assertEqual(len(set(M.roster())), 288)

    def test_field_multiply_neutral_examples(self):
        for x in range(256):
            self.assertEqual(M.multiply(x, 0), 0)
            self.assertEqual(M.multiply(x, 1), x)
            self.assertEqual(M.multiply(x, 2), ((x << 1) ^ (333 if x & 128 else 0)))

    def test_build_never_contains_worker_invocation(self):
        commands = M.build_commands(Path("/tmp/neutral-mulrows-build"))
        self.assertEqual(len(commands), 6)
        self.assertTrue(all("--worker" not in command for command in commands))
        for command in commands[:4]:
            self.assertIn("-std=gnu++11", command)
            self.assertNotIn("-flto", command)

    def test_build_manifest_rejects_stale_inputs_and_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "source", Path(temp) / "object"
            M.C.write_new(source, b"neutral source")
            M.C.write_new(output, b"neutral object")
            manifest = dict(inputs=[M.N.pin(source)], outputs=[M.N.pin(output)])
            M.validate_build_manifest(manifest, {source}, {output})
            for key in ("inputs", "outputs"):
                changed = copy.deepcopy(manifest)
                changed[key][0]["sha256"] = "0" * 64
                with self.assertRaisesRegex(ValueError, "stale build"):
                    M.validate_build_manifest(changed, {source}, {output})
            with self.assertRaisesRegex(ValueError, "dependency target"):
                M.dependency_names(b"wrong.o: source.cpp\n", "right.o")

    def test_expired_deadline_is_terminal_before_oracle(self):
        with mock.patch.object(M, "authenticated_lookup") as lookup, self.assertRaises(TimeoutError):
            M.reduce_raw(self.base, deadline=0)
        lookup.assert_not_called()
        with self.assertRaises(TimeoutError):
            M.expected_rows(self.lookup, deadline=0)
        with self.assertRaises(TimeoutError):
            M.output_checksum(bytes(1536), deadline=0)

    def test_observed_worker_limits(self):
        for elapsed in (30.01, .0001):
            with tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / "claimed"
                receipt_path = Path(temp) / "receipt.json"
                receipt = dict(build_directory="/tmp/neutral-unused", worker_argv=["never-executed", "--worker"])
                M.C.write_new(receipt_path, M.C.canonical(receipt))
                observation = dict(stdout=M.C.canonical(self.base) + b"\n", stderr=b"",
                                   failure=None, returncode=0, elapsed_seconds=elapsed)
                with mock.patch.object(M, "OUTPUT", output), \
                        mock.patch.object(M, "current_receipt", return_value=receipt), \
                        mock.patch.object(M.C, "capture_worker", return_value=observation), \
                        mock.patch.object(M, "reduce_raw") as reducer, \
                        mock.patch("sys.stdout", new_callable=io.StringIO):
                    self.assertEqual(M.run(receipt_path), 1)
                    reducer.assert_not_called()
                summary = M.C.strict_json(M.C.read_regular(output / "summary.json", 65536))
                self.assertEqual(summary["outcome"], "INVALID")

    def test_create_only_terminal_and_invalid_post_pin(self):
        for post_error in (False, True):
            with tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / "claimed"
                receipt_path = Path(temp) / "receipt.json"
                receipt = dict(build_directory="/tmp/neutral-unused", worker_argv=["never-executed", "--worker"])
                M.C.write_new(receipt_path, M.C.canonical(receipt))
                observation = dict(stdout=M.C.canonical(self.base) + b"\n", stderr=b"",
                                   failure=None, returncode=0, elapsed_seconds=.1)
                result = self.reduce(self.base)
                with mock.patch.object(M, "OUTPUT", output), \
                        mock.patch.object(M, "current_receipt", side_effect=[receipt, ValueError("post") if post_error else receipt]), \
                        mock.patch.object(M.C, "capture_worker", return_value=observation) as capture, \
                        mock.patch.object(M, "reduce_raw", return_value=result), \
                        mock.patch("sys.stdout", new_callable=io.StringIO):
                    self.assertEqual(M.run(receipt_path), int(post_error))
                    self.assertEqual(capture.call_count, 1)
                manifest = M.C.strict_json(M.C.read_regular(output / "COMPLETE.json", 65536))
                self.assertEqual(manifest["outcome"], "INVALID" if post_error else "PASS")
                self.assertEqual(set(manifest["files"]), {"CLAIM.json", "raw.json", "stderr.txt", "analysis.json", "summary.json"})
                for name, info in manifest["files"].items():
                    data = M.C.read_regular(output / name, 2 * 1024 * 1024)
                    self.assertEqual(info, dict(bytes=len(data), sha256=M.C.sha(data)))
                with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current_receipt", return_value=receipt), \
                        mock.patch.object(M.C, "capture_worker") as capture, self.assertRaises(FileExistsError):
                    M.run(receipt_path)
                capture.assert_not_called()


if __name__ == "__main__":
    unittest.main()
