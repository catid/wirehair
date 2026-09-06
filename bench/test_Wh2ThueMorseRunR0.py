#!/usr/bin/env python3
"""Neutral .62 wrapper tests plus the unchanged controller regression suite.

No frozen polynomial, parameter selection, matrix candidate, or codec is run.
"""

import contextlib
import importlib.util
import io
from pathlib import Path
import stat
import sys
import tempfile
import time
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def load_sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


W = load_sibling("_thue_morse_controller_wrapper", "Wh2ThueMorseRunR0.py")
HISTORY = load_sibling("_thue_morse_controller_regressions",
                       "test_Wh2NoncommutingRadixRunR0.py")
HISTORY.M = W.C


class WrapperTests(unittest.TestCase):
    def test_exact_namespace_and_transitive_source_pins(self):
        self.assertEqual(W.C.PROTOCOL, "wirehair.wh2.thue-morse-r0")
        self.assertEqual(str(W.C.OUTPUT), "/var/tmp/wh2-thue-morse-r0")
        self.assertEqual(W.C.SOURCES, (
            "bench/Wh2ThueMorseR0.py",
            "bench/Wh2ThueMorseRunR0.py",
            "bench/test_Wh2ThueMorseR0.py",
            "bench/test_Wh2ThueMorseRunR0.py",
            "bench/Wh2NoncommutingRadixR0.py",
            "bench/Wh2NoncommutingRadixRunR0.py",
            "bench/test_Wh2NoncommutingRadixR0.py",
            "bench/test_Wh2NoncommutingRadixRunR0.py",
        ))
        self.assertEqual(W.C.MAX_STDOUT, 4 * 1024 * 1024)
        self.assertEqual(W.C.MAX_STDERR, 1024 * 1024)
        self.assertEqual(W.C.MAX_RECEIPTS, 8 * 1024 * 1024)

    def test_import_is_inert_and_historical_namespace_is_unchanged(self):
        with mock.patch.object(W.C.subprocess, "Popen") as spawn:
            fresh = load_sibling("_thue_morse_second_wrapper", "Wh2ThueMorseRunR0.py")
            old = load_sibling("_radix_original_controller", "Wh2NoncommutingRadixRunR0.py")
        spawn.assert_not_called()
        self.assertIsNot(fresh.C, W.C)
        self.assertEqual(fresh.C.PROTOCOL, W.C.PROTOCOL)
        self.assertEqual(old.PROTOCOL, "wirehair.wh2.noncommuting-radix-r0")
        self.assertEqual(str(old.OUTPUT), "/var/tmp/wh2-noncommuting-radix-r0")
        self.assertEqual(old.SOURCES[0], "bench/Wh2NoncommutingRadixR0.py")
        self.assertEqual(old.strict_json(b'{"protocol":"wirehair.wh2.noncommuting-radix-r0","outcome":"EXHAUSTED"}')["outcome"],
                         "EXHAUSTED")

    def test_receipt_uses_new_worker_and_all_pins(self):
        c = W.C
        data = b"neutral committed source"
        with mock.patch.object(c, "git", side_effect=["", "\n".join(c.SOURCES), "commit", "commit"]), \
                mock.patch.object(c, "read_regular", return_value=data), \
                mock.patch.object(c, "git_bytes", return_value=data), \
                mock.patch.object(c, "capture_worker") as capture:
            receipt = c.current_receipt()
        capture.assert_not_called()
        self.assertEqual(receipt["protocol"], c.PROTOCOL)
        self.assertEqual(set(receipt["sources"]), set(c.SOURCES))
        self.assertEqual(receipt["worker_argv"][1:],
                         ["-I", "-B", "-S", str(HERE / "Wh2ThueMorseR0.py"), "--worker"])
        self.assertEqual(receipt["worker_seconds"], 60)
        self.assertEqual(receipt["outer_seconds"], 70)
        self.assertEqual(receipt["environment"],
                         {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C", "TZ": "UTC"})

    def test_receipt_mode_does_not_launch(self):
        receipt = {"neutral": True}
        with mock.patch.object(sys, "argv", ["wrapper", "--make-receipt", "/neutral/receipt"]), \
                mock.patch.object(W.C, "current_receipt", return_value=receipt), \
                mock.patch.object(W.C, "write_new") as publish, \
                mock.patch.object(W.C, "capture_worker") as capture:
            self.assertEqual(W.C.main(), 0)
        publish.assert_called_once_with(Path("/neutral/receipt"), W.C.canonical(receipt))
        capture.assert_not_called()

    def test_run_mode_delegates_without_extra_work(self):
        with mock.patch.object(sys, "argv", ["wrapper", "--run", "--receipt", "/neutral/receipt"]), \
                mock.patch.object(W.C, "run", return_value=7) as run:
            self.assertEqual(W.C.main(), 7)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args[0][0], Path("/neutral/receipt"))


class PublicationTests(HISTORY.PublicationTests):
    def observation(self, result=None, **changes):
        if result is None:
            result = {"protocol": W.C.PROTOCOL, "outcome": "PASS"}
        return super().observation(result=result, **changes)

    def test_success_manifest_and_no_overwrite(self):
        c = W.C
        with tempfile.TemporaryDirectory(prefix="wh2-thue-publish-test-") as root:
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

    def test_exhaustion_is_invalid_with_raw_evidence_preserved(self):
        observation = self.observation(result={"protocol": W.C.PROTOCOL,
                                               "outcome": "EXHAUSTED"})
        raw = observation["stdout"]
        with tempfile.TemporaryDirectory(prefix="wh2-thue-exhaustion-test-") as root:
            code, output, _ = self.run_fake(root, observation)
            self.assertEqual(code, 1)
            self.assertEqual(W.C.read_regular(output / "raw.json", 65536), raw)
            summary = W.C.strict_json(W.C.read_regular(output / "summary.json", 65536))
            self.assertEqual(summary["outcome"], "INVALID")
            self.assertIn("EXHAUSTED", summary["failure"])


def load_tests(loader, tests, pattern):
    # Each inherited test resolves M in its original module globals. Point it
    # at this wrapper's separately loaded controller, not the historical one.
    return unittest.TestSuite((tests,
                              loader.loadTestsFromTestCase(HISTORY.FileTests),
                              loader.loadTestsFromTestCase(HISTORY.CaptureTests)))


if __name__ == "__main__":
    unittest.main()
