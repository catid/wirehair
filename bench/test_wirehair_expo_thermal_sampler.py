#!/usr/bin/env python3
"""Deterministic tests for fail-closed SPD5118 plausibility accounting."""

import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock


SOURCE = Path(__file__).with_name("wirehair_expo_thermal_sampler.py")
SPEC = importlib.util.spec_from_file_location("wh2_thermal_sampler", SOURCE)
SAMPLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SAMPLER)

WAVE94_RAW_FIXTURE = (
    b"2026-08-11T12:38:19.520Z,2487309.324575,62.70802406438003,"
    b"2912.8568203125,69.0,40.75,42.5,44.0,44.0,43.0,44.5,45.0,"
    b"42.75,0,75.044921875,73.958984375,68.51904296875,0,0\n"
    b"2026-08-11T12:38:20.520Z,2487310.324568,62.52050941479803,"
    b"2953.5854140625,69.0,121.0,42.5,44.0,43.75,43.0,44.5,45.0,"
    b"42.75,0,75.044921875,73.958984375,68.51904296875,0,0\n"
    b"2026-08-11T12:38:21.520Z,2487311.324574,62.65333229158528,"
    b"3035.928890625,69.0,40.75,42.5,44.25,43.75,43.0,44.5,45.0,"
    b"42.75,0,75.044921875,73.958984375,68.51904296875,0,0\n"
)
WAVE94_RAW_FIXTURE_SHA256 = (
    "9696f19f45f49c634e3f8720825edee4c88c2129b1ed5317082bdfec74fca2fc")


def temperatures(value):
    return {sensor: float(value) for sensor in SAMPLER.DIMMS}


def no_attempt_errors():
    return {sensor: 0 for sensor in SAMPLER.DIMMS}


class PlausibilityMonitorTests(unittest.TestCase):
    def observe(self, monitor, when, values, attempts=None, ce=0, ue=0):
        return monitor.observe(
            float(when), values,
            no_attempt_errors() if attempts is None else attempts,
            ce, ue)

    def test_steady_samples_remain_valid(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        decisions = []
        for index, value in enumerate((40.0, 40.25, 40.0, 40.5, 40.25)):
            result = self.observe(monitor, index, temperatures(value))
            decisions.append(result["decision"])
            self.assertEqual(result["fault_count"], 0)
        self.assertEqual(decisions, ["continue"] * 5)
        summary = monitor.summary()
        self.assertEqual(summary["dimm_invalid_samples_total"], 0)
        self.assertEqual(summary["dimm_read_error_samples_total"], 0)
        self.assertEqual(summary["sample_count"], 5)

    def test_legitimate_ramp_below_safety_threshold_remains_valid(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        for index, value in enumerate((40.0, 44.0, 48.0, 52.0, 56.0)):
            result = self.observe(monitor, index, temperatures(value))
            self.assertEqual(result["decision"], "continue")
            self.assertEqual(result["fault_count"], 0)
        self.assertEqual(
            monitor.summary()["dimm_invalid_samples_total"], 0)

    def test_sustained_real_overheat_aborts_after_three_coherent_samples(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        results = [
            self.observe(monitor, index, temperatures(value))
            for index, value in enumerate((80.0, 86.0, 90.0, 91.0, 92.0))
        ]
        self.assertEqual(
            [result["decision"] for result in results],
            ["continue", "continue", "continue", "continue", "thermal_abort"])
        self.assertEqual(results[-1]["fault_count"], 0)
        self.assertEqual(len(results[-1]["hot_sensors"]), len(SAMPLER.DIMMS))
        self.assertEqual(
            monitor.summary()["dimm_invalid_samples_total"], 0)

    def test_sustained_hot_jump_is_confirmed_instead_of_dismissed_forever(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        self.observe(monitor, 0, temperatures(40.0))
        results = [
            self.observe(monitor, index, temperatures(121.0))
            for index in (1, 2, 3)
        ]
        self.assertEqual(
            [result["decision"] for result in results],
            ["continue", "continue", "thermal_abort"])
        self.assertEqual(
            monitor.summary()["dimm_invalid_samples_total"],
            3 * len(SAMPLER.DIMMS))

    def test_repeated_above_plausibility_cap_still_thermal_aborts(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        self.observe(monitor, 0, temperatures(40.0))
        results = [
            self.observe(monitor, index, temperatures(150.0))
            for index in (1, 2, 3)
        ]
        self.assertEqual(results[-1]["decision"], "thermal_abort")
        first_sensor = SAMPLER.dimm_name(SAMPLER.DIMMS[0])
        self.assertEqual(
            results[-1]["sensors"][first_sensor]["reason"],
            "absolute_range")

    def test_wave94_single_spike_is_preserved_and_flagged_exactly_once(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        target = SAMPLER.DIMMS[0]
        self.assertEqual(
            hashlib.sha256(WAVE94_RAW_FIXTURE).hexdigest(),
            WAVE94_RAW_FIXTURE_SHA256)
        results = []
        for row in csv.reader(
                WAVE94_RAW_FIXTURE.decode("ascii").splitlines()):
            values = {
                sensor: float(value)
                for sensor, value in zip(SAMPLER.DIMMS, row[5:13])
            }
            results.append(self.observe(monitor, float(row[1]), values))
        initial_result, spike_result, recovered_result = results

        target_name = SAMPLER.dimm_name(target)
        self.assertEqual(initial_result["fault_count"], 0)
        self.assertEqual(spike_result["fault_count"], 1)
        self.assertEqual(spike_result["sensors"][target_name]["raw_c"], 121.0)
        self.assertEqual(
            spike_result["sensors"][target_name]["reason"], "jump+rate")
        self.assertFalse(spike_result["sensors"][target_name]["valid"])
        self.assertEqual(spike_result["decision"], "continue")
        self.assertEqual(recovered_result["fault_count"], 0)
        self.assertEqual(recovered_result["decision"], "continue")
        self.assertEqual(
            monitor.summary()["dimm_invalid_samples_total"], 1)

    def test_read_failure_and_attempt_errors_are_explicit_then_recover(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        target = SAMPLER.DIMMS[3]
        self.observe(monitor, 0, temperatures(42.0))
        missing = temperatures(42.0)
        del missing[target]
        attempts = no_attempt_errors()
        attempts[target] = 5
        result = self.observe(monitor, 1, missing, attempts)
        target_result = result["sensors"][SAMPLER.dimm_name(target)]
        self.assertEqual(result["read_error_count"], 1)
        self.assertEqual(target_result["reason"], "read_error")
        self.assertEqual(target_result["attempt_errors"], 5)
        self.assertEqual(result["decision"], "continue")
        recovered = self.observe(monitor, 2, temperatures(42.25))
        self.assertEqual(recovered["consecutive_fault_rows"], 0)
        summary = monitor.summary()
        self.assertEqual(summary["dimm_attempt_errors_total"], 5)
        self.assertEqual(summary["dimm_read_error_samples_total"], 1)

    def test_sustained_read_failure_fails_closed(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        target = SAMPLER.DIMMS[-1]
        self.observe(monitor, 0, temperatures(42.0))
        decisions = []
        for index in range(1, SAMPLER.TELEMETRY_FAULT_ABORT_SAMPLES + 1):
            missing = temperatures(42.0)
            del missing[target]
            decisions.append(self.observe(
                monitor, index, missing)["decision"])
        self.assertEqual(
            decisions[:-1],
            ["continue"] * (SAMPLER.TELEMETRY_FAULT_ABORT_SAMPLES - 1))
        self.assertEqual(decisions[-1], "telemetry_abort")

    def test_edac_increase_aborts_and_counter_decrease_is_rejected(self):
        monitor = SAMPLER.DimmPlausibilityMonitor(10, 20)
        result = self.observe(monitor, 0, temperatures(40), ce=11, ue=20)
        self.assertEqual(result["decision"], "edac_abort")
        self.assertEqual(result["edac_ce_delta"], 1)
        decreasing = SAMPLER.DimmPlausibilityMonitor(10, 20)
        self.observe(decreasing, 0, temperatures(40), ce=10, ue=20)
        with self.assertRaisesRegex(ValueError, "must not decrease"):
            self.observe(decreasing, 1, temperatures(40), ce=9, ue=20)

    def test_replay_validation_bytes_are_deterministic(self):
        fixture = (
            (0.0, 40.75),
            (1.0, 121.0),
            (2.0, 40.75),
        )
        outputs = []
        for _ in range(2):
            monitor = SAMPLER.DimmPlausibilityMonitor()
            records = []
            for when, target_value in fixture:
                values = temperatures(40.75)
                values[SAMPLER.DIMMS[0]] = target_value
                records.append(SAMPLER.canonical_json_bytes(
                    self.observe(monitor, when, values)))
            outputs.append(b"".join(records))
        self.assertEqual(outputs[0], outputs[1])

    def test_nonfinite_and_nonmonotonic_inputs_are_rejected(self):
        monitor = SAMPLER.DimmPlausibilityMonitor()
        self.observe(monitor, 1, temperatures(40))
        with self.assertRaisesRegex(ValueError, "must increase"):
            self.observe(monitor, 1, temperatures(40))
        bad = temperatures(40)
        bad[SAMPLER.DIMMS[0]] = float("nan")
        with self.assertRaisesRegex(ValueError, "must be finite"):
            self.observe(
                SAMPLER.DimmPlausibilityMonitor(), 0, bad)
        with self.assertRaisesRegex(ValueError, "must be finite"):
            SAMPLER.DimmPlausibilityMonitor().observe(
                True, temperatures(40), no_attempt_errors(), 0, 0)
        with self.assertRaises(ValueError):
            SAMPLER.DimmPlausibilityMonitor(True, 0)

    def test_malformed_final_sensor_leaves_monitor_state_unchanged(self):
        monitor = SAMPLER.DimmPlausibilityMonitor(7, 9)
        before = copy.deepcopy(monitor.summary())
        attempts = no_attempt_errors()
        attempts[SAMPLER.DIMMS[-1]] = -1
        with self.assertRaisesRegex(ValueError, "nonnegative integer"):
            self.observe(
                monitor, 0, temperatures(40), attempts, ce=7, ue=9)
        self.assertEqual(monitor.summary(), before)

    def test_simultaneous_thermal_and_edac_abort_prioritizes_thermal(self):
        monitor = SAMPLER.DimmPlausibilityMonitor(4, 5)
        self.observe(monitor, 0, temperatures(40), ce=4, ue=5)
        self.observe(monitor, 1, temperatures(95), ce=4, ue=5)
        self.observe(monitor, 2, temperatures(95), ce=4, ue=5)
        result = self.observe(monitor, 3, temperatures(95), ce=5, ue=5)
        self.assertEqual(result["decision"], "thermal_abort")
        self.assertEqual(result["edac_ce_delta"], 1)


class AcquisitionTests(unittest.TestCase):
    def test_i2c_failures_are_counted_and_retried(self):
        target_fd = 101
        target_address = 0x50
        failures = 0

        def fake_read(fd, address):
            nonlocal failures
            if fd == target_fd and address == target_address and failures < 2:
                failures += 1
                raise OSError("synthetic I2C failure")
            return 40.0 + (address - 0x50) / 4.0

        with mock.patch.object(
                SAMPLER, "read_spd5118_temperature", side_effect=fake_read), \
                mock.patch.object(SAMPLER.time, "sleep") as sleep:
            values, pending, attempts = SAMPLER.read_dimm_temperatures(
                {1: target_fd, 2: 102}, 5, 0.01)
        self.assertEqual(pending, [])
        self.assertEqual(attempts[(1, 0x50)], 2)
        self.assertEqual(values[(1, 0x50)], 40.0)
        self.assertEqual(sleep.call_count, 2)

    def test_implausible_decoded_value_is_not_retried_or_dropped(self):
        def fake_read(fd, address):
            if fd == 101 and address == 0x50:
                return 121.0
            return 40.0

        with mock.patch.object(
                SAMPLER, "read_spd5118_temperature", side_effect=fake_read) \
                as read:
            values, pending, attempts = SAMPLER.read_dimm_temperatures(
                {1: 101, 2: 102}, 5, 0.01)
        self.assertEqual(read.call_count, len(SAMPLER.DIMMS))
        self.assertEqual(values[(1, 0x50)], 121.0)
        self.assertEqual(pending, [])
        self.assertEqual(sum(attempts.values()), 0)


class EvidenceTests(unittest.TestCase):
    RAW_FIXTURE = (
        b"monotonic_s,dimm_i2c1_50_c\n"
        b"0.000000,40.75\n"
        b"1.000000,121.0\n"
        b"2.000000,40.75\n"
    )

    def write_sealed(self, path, payload):
        with SAMPLER.open_exclusive_evidence(path, newline="") as output:
            output.write(payload.decode("ascii"))
            binding = SAMPLER.seal_and_hash_output(output)
            SAMPLER.validate_path_binding(path, binding)
            return binding

    def test_terminal_error_records_are_nonempty_bounded_and_capped(self):
        class BrokenTextError(Exception):
            def __str__(self):
                raise RuntimeError("synthetic string conversion failure")

        errors = []
        SAMPLER.append_terminal_error(errors, "empty", MemoryError())
        SAMPLER.append_terminal_error(
            errors, "long", OSError("x" * 5000))
        SAMPLER.append_terminal_error(errors, "", BrokenTextError())
        self.assertEqual(errors[0]["message"], "MemoryError")
        self.assertEqual(len(errors[1]["message"]), 4000)
        self.assertEqual(errors[2]["message"], "BrokenTextError")
        self.assertEqual(errors[2]["phase"], "unknown")
        capped = []
        for index in range(SAMPLER.MAX_TERMINAL_ERRORS + 3):
            SAMPLER.append_terminal_error(capped, f"phase-{index}", OSError(index))
        self.assertEqual(len(capped), SAMPLER.MAX_TERMINAL_ERRORS)
        self.assertEqual(capped[-1]["phase"], "phase-63")

    def test_raw_fixture_is_immutable_and_hash_reproducible(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-evidence-") as temp:
            first = os.path.join(temp, "first.csv")
            second = os.path.join(temp, "second.csv")
            first_binding = self.write_sealed(first, self.RAW_FIXTURE)
            second_binding = self.write_sealed(second, self.RAW_FIXTURE)
            expected = hashlib.sha256(self.RAW_FIXTURE).hexdigest()
            self.assertEqual(first_binding["sha256"], expected)
            self.assertEqual(second_binding["sha256"], expected)
            self.assertEqual(first_binding["mode"], 0o444)
            self.assertEqual(Path(first).read_bytes(), self.RAW_FIXTURE)

    def test_terminal_receipt_is_canonical_sealed_and_self_bound(self):
        fixture = {
            "outcome": "stopped",
            "raw_samples_preserved": True,
            "schema": SAMPLER.SAMPLER_SCHEMA,
            "summary": {"dimm_invalid_samples_total": 1},
            "thresholds": SAMPLER.threshold_metadata(),
        }
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-receipt-") as temp:
            paths = [os.path.join(temp, name) for name in ("a.json", "b.json")]
            bindings = []
            for path in paths:
                bindings.append(SAMPLER.write_terminal_receipt(
                    path, copy.deepcopy(fixture)))
            self.assertEqual(bindings[0]["sha256"], bindings[1]["sha256"])
            receipt = json.loads(Path(paths[0]).read_text(encoding="ascii"))
            self_hash = receipt.pop("self_sha256_excluding_field")
            self.assertEqual(
                self_hash,
                hashlib.sha256(SAMPLER.canonical_json_bytes(receipt)).hexdigest())
            self.assertEqual(os.stat(paths[0]).st_mode & 0o7777, 0o444)

    def test_intent_marker_is_bound_then_replaced_in_the_same_receipt_inode(self):
        fixture = {
            "outcome": "thermal_abort",
            "schema": SAMPLER.SAMPLER_SCHEMA,
        }

        class NonzeroSizeSpy:
            def __init__(self, stream):
                self.stream = stream
                self.observed_sizes = []

            def __getattr__(self, name):
                return getattr(self.stream, name)

            def _observe(self):
                self.observed_sizes.append(os.fstat(self.stream.fileno()).st_size)

            def flush(self):
                result = self.stream.flush()
                self._observe()
                return result

            def seek(self, *args):
                result = self.stream.seek(*args)
                self._observe()
                return result

            def truncate(self, *args):
                result = self.stream.truncate(*args)
                self._observe()
                return result

            def write(self, *args):
                result = self.stream.write(*args)
                self._observe()
                return result

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-intent-") as temp:
            path = os.path.join(temp, "receipt.json")
            destination = SAMPLER.prepare_evidence_destination(path)
            output = SAMPLER.open_exclusive_evidence(
                path, newline="", destination=destination)
            try:
                reservation = SAMPLER.binding_from_fd(output.fileno())
                SAMPLER.write_terminal_intent_marker(output)
                self.assertEqual(Path(path).read_bytes(), b"!")
                with self.assertRaisesRegex(
                        RuntimeError, "receipt inode is not empty"):
                    SAMPLER.write_terminal_intent_marker(output)
                with self.assertRaisesRegex(
                        RuntimeError, "receipt content differs"):
                    SAMPLER.validate_terminal_intent_state(
                        path, output, reservation, False, destination)
                marked = SAMPLER.validate_terminal_intent_state(
                    path, output, reservation, True, destination)
                self.assertEqual(marked["size"], 1)
                self.assertEqual(marked["inode"], reservation["inode"])
                size_spy = NonzeroSizeSpy(output)
                final = SAMPLER.write_terminal_receipt(
                    path, fixture, destination=destination, output=size_spy)
                self.assertEqual(final["inode"], reservation["inode"])
                self.assertGreater(final["size"], 1)
                self.assertTrue(size_spy.observed_sizes)
                self.assertGreater(min(size_spy.observed_sizes), 0)
            finally:
                output.close()
                SAMPLER.close_evidence_destination(destination)
            receipt_bytes = Path(path).read_bytes()
            self.assertTrue(receipt_bytes.startswith(b"{"))
            self.assertTrue(receipt_bytes.endswith(b"\n"))
            self.assertNotIn(b"!", receipt_bytes)
            receipt = json.loads(receipt_bytes)
            self_hash = receipt.pop("self_sha256_excluding_field")
            self.assertEqual(
                self_hash,
                hashlib.sha256(SAMPLER.canonical_json_bytes(receipt)).hexdigest())

    def test_destinations_must_be_absent_and_distinct(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-paths-") as temp:
            same = os.path.join(temp, "same")
            with self.assertRaisesRegex(ValueError, "must be distinct"):
                SAMPLER.prepare_output_paths([same, same])
            Path(same).write_text("occupied", encoding="ascii")
            with self.assertRaises(FileExistsError):
                SAMPLER.prepare_output_paths([same])

    def test_output_parent_must_match_externally_trusted_owner_uid(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-owner-") as temp:
            path = os.path.join(temp, "evidence")
            wrong_uid = 1 if os.geteuid() == 0 else 0
            with self.assertRaisesRegex(
                    PermissionError, "externally trusted UID"):
                SAMPLER.prepare_evidence_destination(
                    path, expected_owner_uid=wrong_uid)

    def test_path_replacement_cannot_validate_as_the_sealed_inode(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-binding-") as temp:
            path = os.path.join(temp, "evidence")
            binding = self.write_sealed(path, b"original\n")
            os.unlink(path)
            Path(path).write_bytes(b"replacement\n")
            with self.assertRaisesRegex(RuntimeError, "binding or content"):
                SAMPLER.validate_path_binding(path, binding)

    def test_same_size_write_through_open_fd_breaks_content_binding(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-content-") as temp:
            path = os.path.join(temp, "evidence")
            output = SAMPLER.open_exclusive_evidence(path, newline="")
            writable_duplicate = os.dup(output.fileno())
            try:
                output.write("original\n")
                binding = SAMPLER.seal_and_hash_output(output)
                output.close()
                os.pwrite(writable_duplicate, b"tampered\n", 0)
                os.fsync(writable_duplicate)
                with self.assertRaisesRegex(RuntimeError, "binding or content"):
                    SAMPLER.validate_path_binding(path, binding)
            finally:
                if not output.closed:
                    output.close()
                os.close(writable_duplicate)

    def test_open_then_path_swap_is_detected_after_hashing(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-open-swap-") as temp:
            path = os.path.join(temp, "evidence")
            binding = self.write_sealed(path, b"original\n")
            real_binding_from_fd = SAMPLER.binding_from_fd

            def hash_then_swap(fd):
                observed = real_binding_from_fd(fd)
                os.unlink(path)
                Path(path).write_bytes(b"replacement\n")
                return observed

            with mock.patch.object(
                    SAMPLER, "binding_from_fd", side_effect=hash_then_swap), \
                    self.assertRaisesRegex(
                        RuntimeError,
                        "changed after validation opened"):
                SAMPLER.validate_path_binding(path, binding)
            self.assertEqual(Path(path).read_bytes(), b"replacement\n")

    def test_fifo_replacement_validation_is_bounded_and_fails_closed(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-fifo-") as temp:
            path = os.path.join(temp, "evidence")
            binding = self.write_sealed(path, b"original\n")
            os.unlink(path)
            os.mkfifo(path, 0o444)
            program = (
                "import importlib.util,json,sys\n"
                "spec=importlib.util.spec_from_file_location('sampler',sys.argv[1])\n"
                "module=importlib.util.module_from_spec(spec)\n"
                "spec.loader.exec_module(module)\n"
                "try:\n"
                " module.validate_path_binding(sys.argv[2],json.loads(sys.argv[3]))\n"
                "except RuntimeError as exc:\n"
                " print(type(exc).__name__+': '+str(exc))\n"
                " raise SystemExit(0)\n"
                "raise SystemExit(2)\n"
            )
            completed = subprocess.run(
                [sys.executable, "-c", program, str(SOURCE), path,
                 json.dumps(binding, sort_keys=True)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=2.0)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn(b"not a regular file", completed.stdout)

    def test_hard_link_breaks_single_link_evidence_invariant(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-nlink-") as temp:
            path = os.path.join(temp, "evidence")
            link = os.path.join(temp, "second-name")
            binding = self.write_sealed(path, b"original\n")
            self.assertEqual(binding["uid"], os.geteuid())
            self.assertEqual(binding["nlink"], 1)
            os.link(path, link)
            with self.assertRaisesRegex(RuntimeError, "binding or content"):
                SAMPLER.validate_path_binding(path, binding)


class RuntimeTests(unittest.TestCase):
    def freeze_source(self, directory):
        path = Path(directory, "frozen_sampler.py")
        path.write_bytes(SOURCE.read_bytes())
        path.chmod(0o444)
        return str(path)

    def arguments(self, directory, interval=0.001):
        return types.SimpleNamespace(
            csv=os.path.join(directory, "raw.csv"),
            dimm_attempts=5,
            dimm_retry_delay=0.0,
            expected_output_owner_uid=os.geteuid(),
            expected_source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            interval=interval,
            pid_file=os.path.join(directory, "sampler.pid"),
            receipt=os.path.join(directory, "receipt.json"),
            validation_jsonl=os.path.join(directory, "validation.jsonl"),
        )

    def test_synthetic_runtime_aborts_hot_and_seals_all_artifacts(self):
        real_os_open = os.open
        real_open_exclusive = SAMPLER.open_exclusive_evidence
        real_write_intent = SAMPLER.write_terminal_intent_marker
        sample_number = 0
        cpu_total = 1000
        observed_raw_modes = []
        observed_intent_bytes = []
        observed_validation_at_intent = []

        class ReceiptCloseFailure:
            def __init__(self, stream):
                self.stream = stream

            def __getattr__(self, name):
                return getattr(self.stream, name)

            def close(self):
                self.stream.close()
                raise OSError("synthetic post-publication close failure")

        def synthetic_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                return real_os_open("/dev/null", os.O_RDWR)
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        def synthetic_dimms(_fds, _attempts, _retry_delay):
            nonlocal sample_number
            observed_raw_modes.append(os.stat(args.csv).st_mode & 0o7777)
            value = 40.0 if sample_number == 0 else 95.0
            sample_number += 1
            return temperatures(value), [], no_attempt_errors()

        def synthetic_cpu_stat():
            nonlocal cpu_total
            cpu_total += 100
            return cpu_total, cpu_total // 10

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-runtime-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)

            def open_with_receipt_close_failure(path, *positional, **keywords):
                output = real_open_exclusive(path, *positional, **keywords)
                if os.path.abspath(path) == os.path.abspath(args.receipt):
                    return ReceiptCloseFailure(output)
                return output

            def observe_intent(output):
                real_write_intent(output)
                observed_intent_bytes.append(Path(args.receipt).read_bytes())
                observed_validation_at_intent.append(
                    Path(args.validation_jsonl).read_bytes())

            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=synthetic_open), \
                    mock.patch.object(
                        SAMPLER, "open_exclusive_evidence",
                        side_effect=open_with_receipt_close_failure), \
                    mock.patch.object(
                        SAMPLER, "write_terminal_intent_marker",
                        side_effect=observe_intent) as intent_write, \
                    mock.patch.object(SAMPLER.signal, "signal"), \
                    mock.patch.object(
                        SAMPLER, "find_tctl_path", return_value="/fake/tctl"), \
                    mock.patch.object(
                        SAMPLER, "read_text", return_value="50000"), \
                    mock.patch.object(
                        SAMPLER, "discover_edac_paths",
                        side_effect=lambda counter: (f"/fake/mc0/{counter}",)), \
                    mock.patch.object(SAMPLER, "sum_edac", return_value=0), \
                    mock.patch.object(
                        SAMPLER, "read_cpu_stat", side_effect=synthetic_cpu_stat), \
                    mock.patch.object(
                        SAMPLER, "average_cpu_mhz", return_value=4000.0), \
                    mock.patch.object(
                        SAMPLER, "read_dimm_temperatures",
                        side_effect=synthetic_dimms):
                return_code = SAMPLER.run_sampler(args)

            self.assertEqual(return_code, SAMPLER.EXIT_THERMAL_ABORT)
            self.assertEqual(intent_write.call_count, 1)
            self.assertEqual(observed_intent_bytes, [b"!"])
            validation_at_intent = observed_validation_at_intent[0].splitlines()
            self.assertEqual(
                json.loads(validation_at_intent[-1])["decision"],
                "thermal_abort")
            self.assertEqual(set(observed_raw_modes), {0o600})
            self.assertFalse(os.path.exists(args.pid_file))
            raw = Path(args.csv).read_bytes()
            rows = raw.decode("ascii").splitlines()
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows[0], (
                "utc,monotonic_s,cpu_busy_pct,cpu_avg_mhz,cpu_tctl_c,"
                "dimm_i2c1_50_c,dimm_i2c1_51_c,dimm_i2c1_52_c,"
                "dimm_i2c1_53_c,dimm_i2c2_50_c,dimm_i2c2_51_c,"
                "dimm_i2c2_52_c,dimm_i2c2_53_c,dimm_read_errors,"
                "load1,load5,load15,edac_ce,edac_ue"
            ))
            self.assertEqual(len(rows[0].split(",")), 19)
            self.assertEqual([len(row.split(",")) for row in rows[1:]], [19] * 4)
            self.assertIn(",95.0,95.0,95.0,95.0,95.0,95.0,95.0,95.0,", raw.decode("ascii"))
            validation_lines = Path(args.validation_jsonl).read_text(
                encoding="ascii").splitlines()
            self.assertEqual(len(validation_lines), 5)
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "thermal_abort")
            self.assertEqual(receipt["exit_code"], SAMPLER.EXIT_THERMAL_ABORT)
            self.assertEqual(receipt["summary"]["sample_count"], 4)
            self.assertEqual(
                receipt["raw_csv"]["binding"]["sha256"],
                hashlib.sha256(raw).hexdigest())
            self.assertEqual(receipt["thresholds"], SAMPLER.threshold_metadata())
            self.assertEqual(
                receipt["sampler_source"]["binding"],
                receipt["sampler_source"]["binding_finished"])
            for path in (args.csv, args.validation_jsonl, args.receipt):
                self.assertEqual(os.stat(path).st_mode & 0o7777, 0o444)

    def test_post_admission_final_timestamp_failure_preserves_paired_streams(self):
        real_os_open = os.open
        sample_number = 0
        cpu_total = 1000

        def synthetic_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                return real_os_open("/dev/null", os.O_RDWR)
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        def synthetic_dimms(_fds, _attempts, _retry_delay):
            nonlocal sample_number
            value = 40.0 if sample_number == 0 else 95.0
            sample_number += 1
            return temperatures(value), [], no_attempt_errors()

        def synthetic_cpu_stat():
            nonlocal cpu_total
            cpu_total += 100
            return cpu_total, cpu_total // 10

        with tempfile.TemporaryDirectory(
                prefix="wh2-thermal-post-admission-final-time-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)
            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=synthetic_open), \
                    mock.patch.object(SAMPLER.signal, "signal"), \
                    mock.patch.object(
                        SAMPLER, "find_tctl_path", return_value="/fake/tctl"), \
                    mock.patch.object(
                        SAMPLER, "read_text", return_value="50000"), \
                    mock.patch.object(
                        SAMPLER, "discover_edac_paths",
                        side_effect=lambda counter: (f"/fake/mc0/{counter}",)), \
                    mock.patch.object(SAMPLER, "sum_edac", return_value=0), \
                    mock.patch.object(
                        SAMPLER, "read_cpu_stat",
                        side_effect=synthetic_cpu_stat), \
                    mock.patch.object(
                        SAMPLER, "average_cpu_mhz", return_value=4000.0), \
                    mock.patch.object(
                        SAMPLER, "read_dimm_temperatures",
                        side_effect=synthetic_dimms), \
                    mock.patch.object(
                        SAMPLER.time, "monotonic_ns",
                        side_effect=(1_000_000_000,
                                     OSError("synthetic final clock failure"))):
                return_code = SAMPLER.run_sampler(args)

            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            raw = Path(args.csv).read_bytes()
            validation = Path(args.validation_jsonl).read_bytes()
            receipt = json.loads(Path(args.receipt).read_bytes())
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertEqual(receipt["exit_code"], SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertEqual(
                [error["phase"] for error in receipt["errors"]],
                ["build_terminal_timestamps"],
            )
            self.assertIsNone(receipt["finished_monotonic_ns"])
            self.assertIsNone(receipt["finished_utc"])
            self.assertTrue(receipt["raw_samples_preserved"])
            self.assertIsNotNone(receipt["raw_csv"]["binding"])
            self.assertIsNotNone(receipt["validation_jsonl"]["binding"])
            self.assertEqual(
                receipt["raw_csv"]["binding"]["sha256"],
                hashlib.sha256(raw).hexdigest(),
            )
            self.assertEqual(
                receipt["validation_jsonl"]["binding"]["sha256"],
                hashlib.sha256(validation).hexdigest(),
            )
            raw_sample_count = len(raw.splitlines()) - 1
            validation_records = [
                json.loads(line)
                for line in validation.splitlines()[1:]
            ]
            self.assertEqual(len(validation_records), raw_sample_count)
            self.assertEqual(validation_records[-1]["decision"], "thermal_abort")
            self.assertEqual(receipt["summary"]["sample_count"], raw_sample_count)
            preimage = dict(receipt)
            observed_self_hash = preimage.pop("self_sha256_excluding_field")
            self.assertEqual(
                observed_self_hash,
                hashlib.sha256(
                    SAMPLER.canonical_json_bytes(preimage)
                ).hexdigest(),
            )

    def test_setup_failure_still_emits_fail_closed_terminal_receipt(self):
        real_os_open = os.open
        real_write_intent = SAMPLER.write_terminal_intent_marker
        observed_intent_bytes = []

        def failing_i2c_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                raise OSError("synthetic I2C unavailable")
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-runtime-error-") \
                as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)

            def observe_intent(output):
                real_write_intent(output)
                observed_intent_bytes.append(Path(args.receipt).read_bytes())

            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=failing_i2c_open), \
                    mock.patch.object(
                        SAMPLER, "write_terminal_intent_marker",
                        side_effect=observe_intent) as intent_write, \
                    mock.patch.object(SAMPLER.signal, "signal"):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertEqual(intent_write.call_count, 1)
            self.assertEqual(observed_intent_bytes, [b"!"])
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertEqual(receipt["exit_code"], SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertFalse(receipt["raw_samples_preserved"])
            self.assertEqual(receipt["errors"][0]["phase"], "open_i2c")
            self.assertIsNone(receipt["raw_csv"]["binding"])
            self.assertEqual(os.stat(args.receipt).st_mode & 0o7777, 0o444)

    def test_final_timestamp_failure_is_retained_as_canonical_sampler_error(self):
        real_os_open = os.open

        def failing_i2c_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                raise OSError("synthetic I2C unavailable")
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        with tempfile.TemporaryDirectory(
                prefix="wh2-thermal-final-time-error-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)
            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=failing_i2c_open), \
                    mock.patch.object(SAMPLER.signal, "signal"), \
                    mock.patch.object(
                        SAMPLER.time, "monotonic_ns",
                        side_effect=(1_000_000_000,
                                     OSError("synthetic final clock failure"))):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            receipt_bytes = Path(args.receipt).read_bytes()
            receipt = json.loads(receipt_bytes)
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertEqual(receipt["started_monotonic_ns"], 1_000_000_000)
            self.assertIsNone(receipt["finished_monotonic_ns"])
            self.assertIsNone(receipt["finished_utc"])
            self.assertIn(
                "build_terminal_timestamps",
                [error["phase"] for error in receipt["errors"]],
            )
            preimage = dict(receipt)
            observed_self_hash = preimage.pop("self_sha256_excluding_field")
            expected_self_hash = hashlib.sha256(
                SAMPLER.canonical_json_bytes(preimage)
            ).hexdigest()
            self.assertEqual(observed_self_hash, expected_self_hash)

    def test_stop_signal_publishes_intent_before_cleanup(self):
        real_os_open = os.open
        real_write_intent = SAMPLER.write_terminal_intent_marker
        observed_intent_bytes = []

        def synthetic_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                return real_os_open("/dev/null", os.O_RDWR)
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        def install_and_stop(signum, handler):
            if signum == SAMPLER.signal.SIGTERM:
                handler(signum, None)

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-stop-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)

            def observe_intent(output):
                real_write_intent(output)
                observed_intent_bytes.append(Path(args.receipt).read_bytes())

            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=synthetic_open), \
                    mock.patch.object(
                        SAMPLER.signal, "signal", side_effect=install_and_stop), \
                    mock.patch.object(
                        SAMPLER, "write_terminal_intent_marker",
                        side_effect=observe_intent) as intent_write, \
                    mock.patch.object(
                        SAMPLER, "find_tctl_path", return_value="/fake/tctl"), \
                    mock.patch.object(
                        SAMPLER, "discover_edac_paths",
                        side_effect=lambda counter: (f"/fake/mc0/{counter}",)), \
                    mock.patch.object(SAMPLER, "sum_edac", return_value=0), \
                    mock.patch.object(
                        SAMPLER, "read_cpu_stat", return_value=(1000, 100)):
                return_code = SAMPLER.run_sampler(args)

            self.assertEqual(return_code, SAMPLER.EXIT_OK)
            self.assertEqual(intent_write.call_count, 1)
            self.assertEqual(observed_intent_bytes, [b"!"])
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "stopped")
            self.assertEqual(receipt["signal"], "SIGTERM")
            self.assertEqual(receipt["exit_code"], SAMPLER.EXIT_OK)

    def test_signal_install_failure_is_receipted(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-signal-error-") \
                as temp:
            args = self.arguments(temp)
            with mock.patch.object(
                    SAMPLER.signal, "signal",
                    side_effect=RuntimeError("synthetic non-main thread")):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertEqual(receipt["errors"][0]["phase"],
                             "install_signal_handlers")

    def test_output_collision_after_receipt_reservation_is_receipted(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-collision-") as temp:
            args = self.arguments(temp)
            Path(args.csv).write_bytes(b"preexisting\n")
            with mock.patch.object(SAMPLER.signal, "signal"):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertEqual(Path(args.csv).read_bytes(), b"preexisting\n")
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["errors"][0]["phase"],
                             "prepare_raw_csv_destination")

    def test_wrong_output_parent_owner_is_receipted(self):
        real_directory_binding = SAMPLER.directory_binding_from_fd
        with tempfile.TemporaryDirectory(
                prefix="wh2-thermal-owner-receipt-") as temp:
            root = Path(temp)
            receipt_parent = root / "trusted"
            wrong_parent = root / "wrong-owner"
            receipt_parent.mkdir(mode=0o700)
            wrong_parent.mkdir(mode=0o700)
            args = self.arguments(str(receipt_parent))
            args.csv = str(wrong_parent / "raw.csv")
            wrong_info = wrong_parent.stat()

            def bind_with_wrong_owner(fd):
                binding = real_directory_binding(fd)
                if (binding["device"], binding["inode"]) == \
                        (wrong_info.st_dev, wrong_info.st_ino):
                    binding = dict(binding)
                    binding["uid"] = (
                        args.expected_output_owner_uid + 1
                    ) % (SAMPLER.MAX_UID + 1)
                return binding

            with mock.patch.object(
                    SAMPLER, "directory_binding_from_fd",
                    side_effect=bind_with_wrong_owner), \
                    mock.patch.object(SAMPLER.signal, "signal"):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["expected_output_owner_uid"], os.geteuid())
            self.assertEqual(receipt["errors"][0]["phase"],
                             "prepare_raw_csv_destination")
            self.assertIn("externally trusted UID",
                          receipt["errors"][0]["message"])

    def test_expected_source_sha_mismatch_is_receipted(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-source-sha-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)
            args.expected_source_sha256 = "0" * 64
            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(SAMPLER.signal, "signal"):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["errors"][0]["phase"], "bind_source")
            self.assertEqual(
                receipt["sampler_source"]["expected_sha256"], "0" * 64)
            self.assertNotEqual(
                receipt["sampler_source"]["sha256"], "0" * 64)

    def test_fifo_source_open_is_bounded_and_receipted(self):
        with tempfile.TemporaryDirectory(prefix="wh2-thermal-source-fifo-") \
                as temp:
            fifo = os.path.join(temp, "sampler-source-fifo")
            os.mkfifo(fifo, 0o444)
            program = (
                "import importlib.util,os,sys,types\n"
                "spec=importlib.util.spec_from_file_location('sampler',sys.argv[1])\n"
                "module=importlib.util.module_from_spec(spec)\n"
                "spec.loader.exec_module(module)\n"
                "module.__file__=sys.argv[2]\n"
                "root=sys.argv[3]\n"
                "args=types.SimpleNamespace(csv=os.path.join(root,'raw.csv'),"
                "pid_file=os.path.join(root,'sampler.pid'),"
                "validation_jsonl=os.path.join(root,'validation.jsonl'),"
                "receipt=os.path.join(root,'receipt.json'),"
                "expected_output_owner_uid=os.geteuid(),"
                "expected_source_sha256='0'*64,interval=1.0,"
                "dimm_attempts=1,dimm_retry_delay=0.0)\n"
                "raise SystemExit(module.run_sampler(args))\n"
            )
            completed = subprocess.run(
                [sys.executable, "-c", program, str(SOURCE), fifo, temp],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=2.0)
            self.assertEqual(completed.returncode, SAMPLER.EXIT_SAMPLER_ERROR)
            receipt = json.loads(
                Path(temp, "receipt.json").read_text(encoding="ascii"))
            self.assertEqual(receipt["errors"][0]["phase"], "bind_source")
            self.assertIn("not a regular file",
                          receipt["errors"][0]["message"])

    def test_pid_seal_failure_is_receipted_and_unverified_pid_is_retained(self):
        real_os_open = os.open
        real_seal = SAMPLER.seal_and_hash_output
        seal_calls = 0

        def synthetic_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                return real_os_open("/dev/null", os.O_RDWR)
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        def fail_first_seal(output):
            nonlocal seal_calls
            seal_calls += 1
            if seal_calls == 1:
                raise RuntimeError("synthetic PID seal failure")
            return real_seal(output)

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-pid-seal-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)
            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=synthetic_open), \
                    mock.patch.object(SAMPLER.signal, "signal"), \
                    mock.patch.object(
                        SAMPLER, "seal_and_hash_output",
                        side_effect=fail_first_seal):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertTrue(os.path.exists(args.pid_file))
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertEqual(
                [error["phase"] for error in receipt["errors"]],
                ["create_pid_file", "unlink_pid_file"])
            self.assertIsNone(receipt["pid_file"]["binding"])

    def test_replaced_pid_path_is_receipted_and_not_unlinked(self):
        real_os_open = os.open
        replacement = b"do-not-unlink-replacement\n"

        def synthetic_open(path, flags, mode=0o777, *, dir_fd=None):
            if str(path).startswith("/dev/i2c-"):
                return real_os_open("/dev/null", os.O_RDWR)
            return real_os_open(path, flags, mode, dir_fd=dir_fd)

        with tempfile.TemporaryDirectory(prefix="wh2-thermal-pid-race-") as temp:
            args = self.arguments(temp)
            frozen_source = self.freeze_source(temp)

            def replace_pid_path():
                os.unlink(args.pid_file)
                Path(args.pid_file).write_bytes(replacement)
                return None

            with mock.patch.object(SAMPLER, "__file__", frozen_source), \
                    mock.patch.object(
                        SAMPLER.os, "open", side_effect=synthetic_open), \
                    mock.patch.object(SAMPLER.signal, "signal"), \
                    mock.patch.object(
                        SAMPLER, "find_tctl_path", side_effect=replace_pid_path):
                return_code = SAMPLER.run_sampler(args)
            self.assertEqual(return_code, SAMPLER.EXIT_SAMPLER_ERROR)
            self.assertEqual(Path(args.pid_file).read_bytes(), replacement)
            receipt = json.loads(Path(args.receipt).read_text(encoding="ascii"))
            self.assertEqual(receipt["outcome"], "sampler_error")
            self.assertIn(
                "unlink_pid_file",
                [error["phase"] for error in receipt["errors"]])


class ContractTests(unittest.TestCase):
    def test_thresholds_and_exit_codes_are_sealed_constants(self):
        self.assertEqual(SAMPLER.threshold_metadata(), {
            "dimm_safety_c_inclusive": 90.0,
            "hot_confirm_samples": 3,
            "max_dimm_jump_c": 12.0,
            "max_dimm_rate_c_per_s": 6.0,
            "max_plausible_dimm_c_exclusive": 130.0,
            "min_plausible_dimm_c_exclusive": 0.0,
            "telemetry_fault_abort_samples": 8,
        })
        self.assertEqual(
            SAMPLER.decision_exit_code("thermal_abort"),
            SAMPLER.EXIT_THERMAL_ABORT)
        self.assertEqual(
            SAMPLER.decision_exit_code("telemetry_abort"),
            SAMPLER.EXIT_TELEMETRY_ABORT)
        self.assertEqual(
            SAMPLER.decision_exit_code("edac_abort"),
            SAMPLER.EXIT_EDAC_ABORT)

    def test_invalid_sampling_arguments_fail(self):
        for interval in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(interval=interval), self.assertRaises(ValueError):
                SAMPLER.validate_sampling_arguments(interval, 1, 0.0)
        with self.assertRaises(ValueError):
            SAMPLER.validate_sampling_arguments(1.0, 0, 0.0)
        with self.assertRaises(ValueError):
            SAMPLER.validate_sampling_arguments(1.0, 1, float("nan"))

    def test_output_owner_uid_argument_is_canonical_and_bounded(self):
        base = [
            "--csv", "/tmp/raw.csv", "--pid-file", "/tmp/pid",
            "--validation-jsonl", "/tmp/validation.jsonl",
            "--receipt", "/tmp/receipt.json",
            "--expected-source-sha256", "0" * 64,
        ]
        for value in ("", "+1", "01", "-1", "1.0",
                      str(SAMPLER.MAX_UID + 1)):
            with self.subTest(value=value), self.assertRaises(SystemExit):
                SAMPLER.parse_arguments(
                    [*base, "--expected-output-owner-uid", value])
        for value in ("0", str(SAMPLER.MAX_UID)):
            parsed = SAMPLER.parse_arguments(
                [*base, "--expected-output-owner-uid", value])
            self.assertEqual(parsed.expected_output_owner_uid, int(value))


if __name__ == "__main__":
    unittest.main(verbosity=2)
