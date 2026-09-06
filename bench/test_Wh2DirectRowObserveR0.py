#!/usr/bin/env python3
"""Neutral .66 schema/math/launcher tests; no fixed data or observation run."""
import copy
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock


SPEC = importlib.util.spec_from_file_location("direct_observe_neutral_test", Path(__file__).with_name("Wh2DirectRowObserveR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
HISTORICAL_READER = M.historical_receipt


def identity():
    data = b"neutral synthetic identity, not a hardware capture"
    derived = dict(ccd_id=6, complex_id=6, core_id=50, family=26,
                   full_apic_id=100, initial_apic_id_8=100,
                   logical_processors_per_package=128, model=8, package_id=0,
                   stepping=1, thread_id=0, threads_per_core=2)
    return dict(requested_cpu=50, before_cpu=50, after_cpu=50,
                raw_capture_complete=True, semantic_validation_passed=True,
                singleton_affinity_verified=True, canonical_bytes=len(data),
                canonical_hex=data.hex(), canonical_sha256=M.C.sha(data),
                capabilities=dict.fromkeys(M.R.H.TARGET_CAPABILITY_KEYS, 0), derived=derived)


def fixture(lookup):
    rows = M.expected_rows(lookup)
    runtime = dict(polynomial=333, address=1048576, ssse3=1, avx2=1, gfni=1, avx512=1)
    raw = dict(protocol=M.PROTOCOL, schema=M.PROTOCOL + ".raw.v1", outcome="DIAGNOSTIC_COMPLETE",
               target_cpu=50, lookup_sha256=M.M.LOOKUP_SHA, expected_rows_hex=rows.hex(),
               expected_rows_sha256=M.C.sha(rows), output_address=65536, output_stride=38,
               baseline_function_address=4096, identity_before=identity(), identity_after=identity(),
               runtime_before=runtime, runtime_after=dict(runtime), monotonic_resolution_ns=1,
               thread_resolution_ns=1, preflight_row_calls=64, callbacks=240, row_calls=983040,
               checked_rows=15360, checksum=M.output_checksum(rows), sum_row_wall_ns=0,
               elapsed_ns=10000000, records=[], failure=None)
    for index, rep, order, position in M.roster():
        m0, c0, wall = 10**12 + index * 20000, 10**6 + index * 4000, 1000 + index
        ru0 = [index // 30, 0, index // 50, index // 70]
        ru1 = [(index + 1) // 30, 0, (index + 1) // 50, (index + 1) // 70]
        raw["records"].append([index, rep, order, position, m0, c0, m0 + 100,
                               m0 + 100 + wall, c0 + wall + 200, m0 + wall + 500, ru0, ru1])
        raw["sum_row_wall_ns"] += wall
    return raw


class ObservationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lookup = bytes((i * 23 + i // 17) & 255 for i in range(39936))
        cls.raw = fixture(cls.lookup)

    def setUp(self):
        clock = mock.patch.object(M.time, "monotonic", return_value=100.0)
        clock.start()
        self.addCleanup(clock.stop)
        # Any accidental fixed lookup/campaign/build entry is a test failure.
        for owner, name in ((M.M, "authenticated_lookup"), (M.M, "run"), (M.M, "reduce_raw"),
                            (M.M, "build_only"), (M, "historical_receipt"),
                            (M.N, "checked"), (M.C, "capture_worker")):
            patch = mock.patch.object(owner, name, side_effect=AssertionError("forbidden actual work"))
            patch.start()
            self.addCleanup(patch.stop)

    def reduce(self, raw):
        return M.reduce_raw(raw, lookup=self.lookup)

    def test_complete_roster_and_descriptive_bins(self):
        result = self.reduce(self.raw)
        self.assertEqual(result["outcome"], "DIAGNOSTIC_COMPLETE")
        self.assertEqual(len(result["records"]), 240)
        self.assertEqual(len(result["bins"]), 24)
        self.assertEqual(sum(row["warmup"] for row in result["records"]), 48)
        self.assertEqual(result["bins"][0]["metrics"]["row_wall_ns"], dict(min=1000, median=1004.5, max=1009))
        self.assertEqual(result["bins"][-1]["metrics"]["row_wall_ns"], dict(min=1230, median=1234.5, max=1239))
        self.assertAlmostEqual(result["first_over_last_bin_ratios"]["row_wall_ns"]["median"], 1004.5 / 1234.5)
        self.assertEqual(sum(row["counter_deltas"]["minflt"] for row in result["records"]), 8)
        self.assertEqual(sum(row["counter_deltas"]["nvcsw"] for row in result["records"]), 4)
        self.assertEqual(sum(row["counter_deltas"]["nivcsw"] for row in result["records"]), 3)
        self.assertTrue(result["all_samples_retained"])
        for key in ("speed_claimed", "AA_equivalence_claimed", "production_promotion_claimed", "WH1_compared"):
            self.assertFalse(result[key])

    def test_different_clock_epochs_and_cpu_above_inner_wall(self):
        raw = copy.deepcopy(self.raw)
        for record in raw["records"]:
            # The CPU epoch may even exceed the monotonic epoch. Never compare
            # their absolute timestamps, only within each clock's domain.
            record[5] += 10**15
            record[8] += 10**15
        result = self.reduce(raw)
        for row in result["records"]:
            self.assertEqual(row["thread_span_ns"], row["row_wall_ns"] + 200)
            self.assertEqual(row["instrumentation_envelope_ns"], 500)
            self.assertNotIn("off_cpu_ns", row)

    def test_zero_denominator_is_null_not_clamped(self):
        raw = copy.deepcopy(self.raw)
        for row in raw["records"]:
            row[8] = row[5]
        self.assertIsNone(self.reduce(raw)["first_over_last_bin_ratios"]["thread_span_ns"]["median"])

    def test_all_samples_retained_including_large_inner_interval(self):
        raw = copy.deepcopy(self.raw)
        # A last-record excursion has no subsequent clock record to move.
        raw["records"][-1][7] += 1000000
        raw["records"][-1][8] += 1000000
        raw["records"][-1][9] += 1000000
        raw["sum_row_wall_ns"] += 1000000
        result = self.reduce(raw)
        self.assertEqual(len(result["records"]), 240)
        self.assertEqual(result["bins"][-1]["metrics"]["row_wall_ns"]["max"], 1001239)

    def test_strict_schema_and_numeric_forgeries(self):
        mutations = [(["outcome"], "PASS"), (["failure"], "error"), (["callbacks"], True),
                     (["records", 0, 0], False), (["records", 0, 2], 1),
                     (["records", 0, 4], True), (["records", 0, 5], 1.0),
                     (["records", 0, 10, 0], False), (["records", 0, 11, 0], 1.5),
                     (["records", 0, 9], None), (["records", 0, 10], [0]),
                     (["row_calls"], 1), (["checked_rows"], 0), (["preflight_row_calls"], 256),
                     (["sum_row_wall_ns"], self.raw["sum_row_wall_ns"] + 1),
                     (["sum_row_wall_ns"], 100000001), (["elapsed_ns"], 1),
                     (["elapsed_ns"], 10000000001), (["thread_resolution_ns"], 0),
                     (["monotonic_resolution_ns"], 1.0), (["output_stride"], 6),
                     (["output_address"], (1 << 64) - 1), (["baseline_function_address"], 0),
                     (["expected_rows_hex"], "00" * 384), (["expected_rows_sha256"], "0" * 64),
                     (["checksum"], 1), (["runtime_after", "address"], 9),
                     (["identity_after", "derived", "full_apic_id"], 50)]
        for path, value in mutations:
            raw = copy.deepcopy(self.raw)
            target = raw
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.reduce(raw)
        for raw in (dict(self.raw, extra=1), dict(self.raw, records=self.raw["records"][:-1])):
            with self.assertRaises(ValueError):
                self.reduce(raw)
        raw = dict(self.raw)
        del raw["records"]
        with self.assertRaises(ValueError):
            self.reduce(raw)

    def test_clock_and_counter_monotonicity(self):
        replacements = [(0, 6, self.raw["records"][0][4] - 1),
                        (0, 7, self.raw["records"][0][6]),
                        (0, 8, self.raw["records"][0][5] - 1),
                        (0, 9, self.raw["records"][0][7] - 1),
                        (1, 4, self.raw["records"][0][9] - 1),
                        (1, 5, self.raw["records"][0][8] - 1),
                        (30, 10, [0, 0, 0, 0]),
                        (30, 11, [0, 0, 0, 0])]
        for index, column, value in replacements:
            raw = copy.deepcopy(self.raw)
            raw["records"][index][column] = value
            with self.subTest(index=index, column=column), self.assertRaises(ValueError):
                self.reduce(raw)
        raw = copy.deepcopy(self.raw)
        raw["records"][-1][5] = raw["records"][0][5] + raw["elapsed_ns"]
        raw["records"][-1][8] = raw["records"][-1][5] + 100
        with self.assertRaisesRegex(ValueError, "thread record span"):
            self.reduce(raw)

    def test_deadline_precedes_any_fixed_oracle(self):
        with self.assertRaises(TimeoutError):
            M.reduce_raw(self.raw, deadline=0)
        with self.assertRaises(TimeoutError):
            M.expected_rows(self.lookup, deadline=0)
        with self.assertRaises(TimeoutError):
            M.output_checksum(bytes(384), deadline=0)

    def test_low_only_oracle_and_independent_checksum(self):
        rows = M.expected_rows(self.lookup)
        expected = bytes(self.lookup[78 * j + offset] for j in range(64) for offset in range(6))
        self.assertEqual(rows, expected)
        value = 1469598103934665603
        for byte in expected * 240:
            value = ((value * 1099511628211) & 0xffffffffffffffff) ^ byte
        self.assertEqual(M.output_checksum(rows), value)
        self.assertEqual(len(set(M.roster())), 240)

    def test_build_contains_only_baseline_and_new_driver(self):
        commands = M.build_commands(Path("/tmp/neutral-observe-build"))
        self.assertEqual(len(commands), 3)
        self.assertEqual([Path(command[command.index("-c") + 1]).name for command in commands[:2]],
                         ["Wh2ThueMorseNativeCodec.cpp", "Wh2DirectRowObserveR0.cpp"])
        for command in commands:
            self.assertNotIn("--worker", command)
            self.assertFalse(any("Wh2ThueMorseMulRowsR0.o" in part for part in command))
            self.assertNotIn("-flto", command)
        self.assertEqual(sum(part.endswith("libwirehair.a") for part in commands[-1]), 1)
        for command in commands[:2]:
            self.assertIn("-std=gnu++11", command)

    def test_manifest_retains_exact_input_output_pins(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "neutral-source", Path(temp) / "neutral-object"
            M.C.write_new(source, b"source")
            M.C.write_new(output, b"object")
            manifest = dict(inputs=[M.N.pin(source)], outputs=[M.N.pin(output)])
            M.M.validate_build_manifest(manifest, {source}, {output})
            for key in ("inputs", "outputs"):
                changed = copy.deepcopy(manifest)
                changed[key][0]["sha256"] = "0" * 64
                with self.assertRaises(ValueError):
                    M.M.validate_build_manifest(changed, {source}, {output})

    def test_historical_interpreter_is_authenticated_separately(self):
        old = dict(source_commit="historical-commit", build_directory="/tmp/neutral-old",
                   interpreter=dict(path="/tmp/neutral-python312", bytes=10, sha256="1" * 64),
                   python_version=[3, 12, 3], preserved="exact-build-data")
        current = dict(old, source_commit="new-commit",
                       interpreter=dict(path="/tmp/neutral-python38", bytes=20, sha256="2" * 64),
                       python_version=[3, 8, 20])
        data = M.C.canonical(old)
        with mock.patch.object(M.C, "read_regular", return_value=data), \
                mock.patch.object(M, "HISTORICAL_SHA", M.C.sha(data)), \
                mock.patch.object(M.M, "current_receipt", side_effect=lambda *a, **k: copy.deepcopy(current)), \
                mock.patch.object(M.N, "pin", return_value=old["interpreter"]) as pin:
            self.assertEqual(HISTORICAL_READER(), old)
            pin.assert_called_once_with(old["interpreter"]["path"], None, installed=True)
            current["preserved"] = "changed-build-data"
            with self.assertRaisesRegex(ValueError, "historical build/source"):
                HISTORICAL_READER()
            current["preserved"] = old["preserved"]
            pin.return_value = dict(old["interpreter"], sha256="3" * 64)
            with self.assertRaisesRegex(ValueError, "historical interpreter"):
                HISTORICAL_READER()

    def run_synthetic(self, output, receipt_path, raw, duration=.1, post_error=False):
        receipt = dict(build_directory="/tmp/neutral-unused", worker_argv=["never-executed", "--worker"])
        M.C.write_new(receipt_path, M.C.canonical(receipt))
        raw_bytes = M.C.canonical(raw) + b"\n"
        observation = dict(stdout=raw_bytes, stderr=b"", failure=None, returncode=0, elapsed_seconds=duration)
        reducer = M.reduce_raw
        with mock.patch.object(M, "OUTPUT", output), \
                mock.patch.object(M, "current_receipt", side_effect=[receipt, ValueError("post-pin") if post_error else receipt]), \
                mock.patch.object(M.C, "capture_worker", return_value=observation) as capture, \
                mock.patch.object(M, "reduce_raw", side_effect=lambda value, **kw: reducer(value, self.lookup, **kw)), \
                mock.patch("sys.stdout", new_callable=io.StringIO):
            code = M.run(receipt_path)
            self.assertEqual(capture.call_count, 1)
        summary = M.C.strict_json(M.C.read_regular(output / "summary.json", 65536))
        self.assertEqual(M.C.read_regular(output / "raw.json", 1048576), raw_bytes)
        return code, summary, receipt

    def test_success_publication_and_create_only_namespace(self):
        with tempfile.TemporaryDirectory() as temp:
            output, receipt_path = Path(temp) / "out", Path(temp) / "receipt"
            code, summary, receipt = self.run_synthetic(output, receipt_path, self.raw)
            self.assertEqual(code, 0)
            self.assertEqual(summary["outcome"], "DIAGNOSTIC_COMPLETE")
            manifest = M.C.strict_json(M.C.read_regular(output / "COMPLETE.json", 65536))
            self.assertEqual(set(manifest["files"]), {"CLAIM.json", "raw.json", "stderr.txt", "analysis.json", "summary.json"})
            for name, pin in manifest["files"].items():
                data = M.C.read_regular(output / name, 1048576)
                self.assertEqual(pin, dict(bytes=len(data), sha256=M.C.sha(data)))
            with mock.patch.object(M, "OUTPUT", output), mock.patch.object(M, "current_receipt", return_value=receipt), \
                    mock.patch.object(M.C, "capture_worker") as capture, self.assertRaises(FileExistsError):
                M.run(receipt_path)
            capture.assert_not_called()

    def test_invalid_prefix_is_preserved_without_success_analysis(self):
        raw = dict(self.raw, outcome="INVALID", failure="neutral failed clock",
                   records=[self.raw["records"][0], [1, 0, 0, 1] + [None] * 8], callbacks=1)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "out"
            code, summary, _ = self.run_synthetic(output, Path(temp) / "receipt", raw)
            self.assertEqual(code, 1)
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertFalse((output / "analysis.json").exists())

    def test_observed_overrun_internal_interval_and_postpins(self):
        for duration, post_error in ((10.01, False), (.0001, False), (.1, True)):
            with self.subTest(duration=duration, post_error=post_error), tempfile.TemporaryDirectory() as temp:
                code, summary, _ = self.run_synthetic(Path(temp) / "out", Path(temp) / "receipt", self.raw,
                                                       duration=duration, post_error=post_error)
                self.assertEqual(code, 1)
                self.assertEqual(summary["outcome"], "INVALID")


if __name__ == "__main__":
    unittest.main()
