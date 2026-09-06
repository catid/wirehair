#!/usr/bin/env python3
"""Neutral receipt/publication tests, with no fresh roots or candidate scores."""

import contextlib
import importlib.util
import io
import os
from pathlib import Path
import stat
import sys
import tempfile
import time
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


W = sibling("_recovery_wrapper_test", "Wh2ThueMorseRecoveryRunR0.py")
HISTORY = sibling("_recovery_controller_tests", "test_Wh2NoncommutingRadixRunR0.py")
HISTORY.M = W.C


class WrapperTests(unittest.TestCase):
    def test_namespace_bounds_and_source_closure(self):
        self.assertEqual(W.C.PROTOCOL, "wirehair.wh2.thue-morse-recovery-r0")
        self.assertEqual(str(W.C.OUTPUT), "/var/tmp/wh2-thue-morse-recovery-r0")
        self.assertEqual(len(W.C.SOURCES), len(set(W.C.SOURCES)))
        self.assertEqual(W.C.SOURCES[0], "bench/Wh2ThueMorseRecoveryR0.py")
        required = (
            "Wh2ThueMorseRecoveryR0.py", "Wh2ThueMorseRecoveryRunR0.py",
            "test_Wh2ThueMorseRecoveryR0.py", "Wh2ThueMorseRecoveryHistoryR0.py",
            "test_Wh2ThueMorseRecoveryHistoryR0.py", "test_Wh2ThueMorseRecoveryRunR0.py",
            "Wh2ThueMorseR0.py", "Wh2NoncommutingRadixR0.py",
            "Wh2NoncommutingRadixRunR0.py", "test_Wh2ThueMorseR0.py",
            "test_Wh2NoncommutingRadixR0.py", "test_Wh2NoncommutingRadixRunR0.py",
            "Wh2FrozenTrace.cpp", "Wh2FrozenTrace.h", "Wh2FrozenTraceTest.cpp",
            "wh2_benchmark_contract_v4.json",
            "wh2_mix2_seed_repair_contract.py", "wh2_precodefail_work_screen.py",
        )
        self.assertEqual(set(W.C.SOURCES), {"bench/" + name for name in required})
        self.assertEqual(W.C.MAX_STDOUT, 4 * 1024 * 1024)
        self.assertEqual(W.C.MAX_STDERR, 1024 * 1024)
        self.assertEqual(W.C.MAX_RECEIPTS, 8 * 1024 * 1024)

    def test_import_is_inert_and_historical_controller_unchanged(self):
        with mock.patch.object(W.C.subprocess, "Popen") as spawn, \
                mock.patch.object(os, "open", side_effect=AssertionError("artifact read")):
            fresh = sibling("_fresh_recovery_wrapper", "Wh2ThueMorseRecoveryRunR0.py")
            old = sibling("_unchanged_radix_controller", "Wh2NoncommutingRadixRunR0.py")
        spawn.assert_not_called()
        self.assertIsNot(fresh.C, W.C)
        self.assertEqual(fresh.C.PROTOCOL, W.C.PROTOCOL)
        self.assertEqual(old.PROTOCOL, "wirehair.wh2.noncommuting-radix-r0")
        self.assertEqual(str(old.OUTPUT), "/var/tmp/wh2-noncommuting-radix-r0")
        self.assertEqual(old.strict_json(b'{"outcome":"EXHAUSTED"}')["outcome"], "EXHAUSTED")

    def test_receipt_binds_historical_inputs_and_deadline(self):
        c = W.C
        data = b"neutral committed source"
        history = {"neutral_sha256": "1" * 64}
        deadline = time.monotonic() + 5
        with mock.patch.object(c, "git", side_effect=["", "\n".join(c.SOURCES), "commit", "commit"]), \
                mock.patch.object(c, "read_regular", return_value=data), \
                mock.patch.object(c, "git_bytes", return_value=data), \
                mock.patch.object(W.H, "input_receipt", return_value=history) as inputs, \
                mock.patch.object(c, "capture_worker") as capture:
            receipt = c.current_receipt(deadline=deadline)
        inputs.assert_called_once_with(deadline=deadline)
        capture.assert_not_called()
        self.assertEqual(receipt["historical_inputs"], history)
        self.assertEqual(set(receipt["sources"]), set(c.SOURCES))
        self.assertEqual(receipt["worker_argv"][1:],
                         ["-I", "-B", "-S", str(HERE / "Wh2ThueMorseRecoveryR0.py"), "--worker"])
        self.assertEqual(receipt["worker_seconds"], 60)
        self.assertEqual(receipt["outer_seconds"], 70)

    def test_bad_historical_inputs_fail_before_launch(self):
        with mock.patch.object(W, "_current_receipt", return_value={}), \
                mock.patch.object(W.H, "input_receipt", side_effect=ValueError("input drift")), \
                mock.patch.object(W.C, "capture_worker") as capture:
            with self.assertRaisesRegex(ValueError, "input drift"):
                W.current_receipt()
        capture.assert_not_called()

    def test_receipt_mode_never_launches(self):
        receipt = {"neutral": True}
        with mock.patch.object(sys, "argv", ["wrapper", "--make-receipt", "/neutral/receipt"]), \
                mock.patch.object(W.C, "current_receipt", return_value=receipt), \
                mock.patch.object(W.C, "write_new") as publish, \
                mock.patch.object(W.C, "capture_worker") as capture:
            self.assertEqual(W.C.main(), 0)
        publish.assert_called_once_with(Path("/neutral/receipt"), W.C.canonical(receipt))
        capture.assert_not_called()


class PublicationTests(HISTORY.PublicationTests):
    def observation(self, result=None, **changes):
        if result is None:
            result = {"protocol": W.C.PROTOCOL, "outcome": "PASS"}
        return super().observation(result=result, **changes)

    def test_success_manifest_and_no_overwrite(self):
        c = W.C
        with tempfile.TemporaryDirectory(prefix="wh2-recovery-publish-test-") as root:
            code, output, capture = self.run_fake(root, self.observation())
            self.assertEqual(code, 0)
            capture.assert_called_once()
            manifest = c.strict_json(c.read_regular(output / "COMPLETE.json", 65536))
            self.assertEqual(manifest["outcome"], "PASS")
            for name, entry in manifest["files"].items():
                path = output / name
                data = c.read_regular(path, c.MAX_RECEIPTS)
                self.assertEqual(len(data), entry["bytes"])
                self.assertEqual(c.sha(data), entry["sha256"])
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o400)
            with mock.patch.object(c, "OUTPUT", output), \
                    mock.patch.object(c, "current_receipt", return_value=self.fake_receipt()), \
                    mock.patch.object(c, "capture_worker") as capture, \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(FileExistsError):
                    c.run(Path(root) / "receipt.json", time.monotonic())
                capture.assert_not_called()

    def test_historical_input_drift_after_worker_is_invalid(self):
        before = dict(self.fake_receipt(), historical_inputs={"sha256": "a" * 64})
        after = dict(self.fake_receipt(), historical_inputs={"sha256": "b" * 64})
        with tempfile.TemporaryDirectory(prefix="wh2-recovery-drift-test-") as root, \
                mock.patch.object(self, "fake_receipt", return_value=before):
            code, output, _ = self.run_fake(root, self.observation(), receipts=[before, after])
            self.assertEqual(code, 1)
            summary = W.C.strict_json(W.C.read_regular(output / "summary.json", 65536))
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertEqual(summary["worker"]["post_pins"], "failed_or_unobserved")

    def test_exhausted_preserves_raw_but_is_invalid(self):
        observation = self.observation(result={"protocol": W.C.PROTOCOL, "outcome": "EXHAUSTED"})
        raw = observation["stdout"]
        with tempfile.TemporaryDirectory(prefix="wh2-recovery-exhausted-test-") as root:
            code, output, _ = self.run_fake(root, observation)
            self.assertEqual(code, 1)
            self.assertEqual(W.C.read_regular(output / "raw.json", 65536), raw)


def load_tests(loader, tests, pattern):
    return unittest.TestSuite((tests, loader.loadTestsFromTestCase(HISTORY.FileTests),
                              loader.loadTestsFromTestCase(HISTORY.CaptureTests)))


if __name__ == "__main__":
    unittest.main()
