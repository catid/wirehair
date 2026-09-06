#!/usr/bin/env python3
"""Neutral launcher tests: no candidate bytes, GF search, or codec workloads."""

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
SPEC = importlib.util.spec_from_file_location(
    "noncommuting_launcher", str(HERE / "Wh2NoncommutingRadixRunR0.py"))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
PYTHON = str(Path(sys.executable).resolve())


class FileTests(unittest.TestCase):
    def test_expired_shared_deadline(self):
        with self.assertRaises(TimeoutError):
            M.time_left(time.monotonic() - 1)

    def test_clean_status_cannot_hide_wrong_commit_bytes(self):
        with mock.patch.object(M, "git", side_effect=["", "\n".join(M.SOURCES), "commit", "commit"]), \
                mock.patch.object(M, "read_regular", return_value=b"worktree bytes"), \
                mock.patch.object(M, "git_bytes", return_value=b"different committed bytes"):
            with self.assertRaisesRegex(ValueError, "committed blob"):
                M.current_receipt()

    def test_canonical_and_strict_json(self):
        data = M.canonical({"b": [1, True], "a": "text"})
        self.assertEqual(data, b'{"a":"text","b":[1,true]}')
        self.assertEqual(M.strict_json(data), {"a": "text", "b": [1, True]})
        for bad in (b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":Infinity}',
                    b'\xff', b'{'):
            with self.assertRaises((ValueError, UnicodeError)):
                M.strict_json(bad)
        with self.assertRaises(ValueError):
            M.canonical({"a": float("nan")})

    def test_regular_exclusive_readonly_and_cap(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-file-test-") as root:
            path = Path(root) / "owned"
            M.write_new(path, b"abc")
            self.assertEqual(M.read_regular(path, 3), b"abc")
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o400)
            with self.assertRaises(FileExistsError):
                M.write_new(path, b"replacement")
            self.assertEqual(M.read_regular(path, 3), b"abc")
            with self.assertRaises(ValueError):
                M.read_regular(path, 2)
            with self.assertRaises(ValueError):
                M.read_regular(Path(root), 100000)

    def test_hardlink_boundary_and_symlink(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-link-test-") as root:
            path = Path(root) / "code"
            M.write_new(path, b"same installed code")
            linked = Path(root) / "linked"
            os.link(str(path), str(linked))
            with self.assertRaises(ValueError):
                M.read_regular(path, 100)
            self.assertEqual(M.read_regular(path, 100, single_link=False),
                             b"same installed code")
            symbolic = Path(root) / "symbolic"
            symbolic.symlink_to(path)
            with self.assertRaises(ValueError):
                M.read_regular(symbolic, 100, single_link=False)

    def test_fifo_does_not_block(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-fifo-test-") as root:
            path = Path(root) / "fifo"
            os.mkfifo(str(path))
            with self.assertRaises(ValueError):
                M.read_regular(path, 100)

    def test_changed_descriptor_rejected(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-drift-test-") as root:
            path = Path(root) / "source"
            M.write_new(path, b"abc")
            real = os.fstat
            count = [0]

            def changed(fd):
                info = real(fd)
                count[0] += 1
                if count[0] == 2:
                    return mock.Mock(st_dev=info.st_dev, st_ino=info.st_ino + 1)
                return info

            with mock.patch.object(M.os, "fstat", side_effect=changed):
                with self.assertRaises(ValueError):
                    M.read_regular(path, 100)


class CaptureTests(unittest.TestCase):
    def capture(self, code, **kwargs):
        return M.capture_worker([PYTHON, "-I", "-B", "-S", "-c", code],
                                time.monotonic() + 5, **kwargs)

    def test_success_and_returncode(self):
        result = self.capture("import sys;sys.stdout.write('ok')")
        self.assertEqual(result["stdout"], b"ok")
        self.assertEqual(result["stderr"], b"")
        self.assertEqual(result["returncode"], 0)
        self.assertIsNone(result["failure"])
        result = self.capture("raise SystemExit(7)")
        self.assertEqual(result["returncode"], 7)

    def test_live_stdout_cap(self):
        result = self.capture("import os;os.write(1,b'x'*10000)", stdout_cap=32)
        self.assertEqual(result["stdout"], b"x" * 32)
        self.assertEqual(result["failure"], "stdout cap")

    def test_live_stderr_cap(self):
        result = self.capture("import os;os.write(2,b'x'*10000)", stderr_cap=32)
        self.assertEqual(result["stderr"], b"x" * 32)
        self.assertEqual(result["failure"], "stderr cap")

    def test_exact_cap_is_accepted(self):
        result = self.capture("import os;os.write(1,b'x'*32)", stdout_cap=32)
        self.assertEqual(result["stdout"], b"x" * 32)
        self.assertIsNone(result["failure"])

    def test_deadline_and_reap(self):
        result = M.capture_worker(
            [PYTHON, "-I", "-B", "-S", "-c", "import time;time.sleep(5)"],
            time.monotonic() + 0.15)
        self.assertEqual(result["failure"], "worker deadline")
        self.assertIsInstance(result["returncode"], int)
        self.assertLess(result["elapsed_seconds"], 3)

    def test_selector_failure_never_spawns(self):
        with mock.patch.object(M.selectors, "DefaultSelector", side_effect=OSError("EMFILE")), \
                mock.patch.object(M.subprocess, "Popen") as spawn:
            with self.assertRaises(OSError):
                self.capture("not executed")
            spawn.assert_not_called()

    def test_spawn_failure_closes_selector(self):
        selector = mock.Mock()
        with mock.patch.object(M.selectors, "DefaultSelector", return_value=selector), \
                mock.patch.object(M.subprocess, "Popen", side_effect=OSError("spawn failure")):
            with self.assertRaises(OSError):
                self.capture("not executed")
        selector.close.assert_called_once()

    def test_post_spawn_setup_failure_reaps_child(self):
        selector = mock.Mock()
        selector.register.side_effect = OSError("registration failed")
        child = mock.Mock()
        child.poll.return_value = None
        with mock.patch.object(M.selectors, "DefaultSelector", return_value=selector), \
                mock.patch.object(M.subprocess, "Popen", return_value=child), \
                mock.patch.object(M.os, "set_blocking"):
            with self.assertRaises(OSError):
                self.capture("not executed")
        child.kill.assert_called_once()
        child.wait.assert_called_once()
        child.stdout.close.assert_called_once()
        child.stderr.close.assert_called_once()

    def test_selector_close_failure_still_reaps_child(self):
        selector = mock.Mock()
        selector.register.side_effect = OSError("registration failed")
        selector.close.side_effect = OSError("close failed")
        child = mock.Mock()
        child.poll.return_value = None
        with mock.patch.object(M.selectors, "DefaultSelector", return_value=selector), \
                mock.patch.object(M.subprocess, "Popen", return_value=child), \
                mock.patch.object(M.os, "set_blocking"):
            with self.assertRaises(OSError):
                self.capture("not executed")
        child.kill.assert_called_once()
        child.wait.assert_called_once()
        child.stdout.close.assert_called_once()
        child.stderr.close.assert_called_once()


class PublicationTests(unittest.TestCase):
    def fake_receipt(self):
        return {"protocol": M.PROTOCOL, "worker_argv": ["never-executed"]}

    def observation(self, result=None, **changes):
        if result is None:
            result = {"protocol": M.PROTOCOL, "outcome": "EXHAUSTED"}
        data = {"stdout": M.canonical(result) + b"\n", "stderr": b"",
                "returncode": 0, "elapsed_seconds": 0.001, "failure": None}
        data.update(changes)
        return data

    def run_fake(self, root, observation, receipts=None, age=0):
        receipt = self.fake_receipt()
        receipt_path = Path(root) / "receipt.json"
        M.write_new(receipt_path, M.canonical(receipt))
        output = Path(root) / "bundle"
        with mock.patch.object(M, "OUTPUT", output), \
                mock.patch.object(M, "current_receipt",
                                  side_effect=receipts or [receipt, receipt]), \
                mock.patch.object(M, "capture_worker", return_value=observation) as capture, \
                contextlib.redirect_stdout(io.StringIO()):
            code = M.run(receipt_path, time.monotonic() - age)
        return code, output, capture

    def test_success_manifest_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-publish-test-") as root:
            code, output, capture = self.run_fake(root, self.observation())
            self.assertEqual(code, 0)
            capture.assert_called_once()
            manifest = M.strict_json(M.read_regular(output / "COMPLETE.json", 65536))
            self.assertEqual(manifest["outcome"], "EXHAUSTED")
            for name, entry in manifest["files"].items():
                path = output / name
                data = M.read_regular(path, M.MAX_RECEIPTS)
                self.assertEqual(len(data), entry["bytes"])
                self.assertEqual(M.sha(data), entry["sha256"])
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o400)
            with mock.patch.object(M, "OUTPUT", output), \
                    mock.patch.object(M, "current_receipt", return_value=self.fake_receipt()), \
                    mock.patch.object(M, "capture_worker") as capture:
                with self.assertRaises(FileExistsError):
                    M.run(Path(root) / "receipt.json", time.monotonic())
                capture.assert_not_called()

    def test_invalid_worker_observations_preserved(self):
        cases = [self.observation(stderr=b"diagnostic"),
                 self.observation(returncode=1),
                 self.observation(failure="stdout cap"),
                 self.observation(stdout=b'{"outcome":"PASS","outcome":"FAIL"}'),
                 self.observation(stdout=b'{}'),
                 self.observation(elapsed_seconds=60.001),
                 self.observation(result={"protocol": "wrong", "outcome": "PASS"}),
                 self.observation(result={"protocol": M.PROTOCOL, "outcome": "MAYBE"})]
        for observation in cases:
            with self.subTest(observation=observation), \
                    tempfile.TemporaryDirectory(prefix="wh2-radix-invalid-test-") as root:
                expected_raw = observation["stdout"]
                code, output, _ = self.run_fake(root, observation)
                self.assertEqual(code, 1)
                self.assertEqual(M.read_regular(output / "raw.json", M.MAX_STDOUT),
                                 expected_raw)
                summary = M.strict_json(M.read_regular(output / "summary.json", 65536))
                self.assertEqual(summary["outcome"], "INVALID")
                self.assertTrue(summary["failure"])

    def test_source_drift_is_invalid(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-source-test-") as root:
            code, output, _ = self.run_fake(root, self.observation(),
                                          receipts=[self.fake_receipt(), {}])
            self.assertEqual(code, 1)
            summary = M.strict_json(M.read_regular(output / "summary.json", 65536))
            self.assertIn("identity changed", summary["failure"])

    def test_late_start_never_launches(self):
        with tempfile.TemporaryDirectory(prefix="wh2-radix-late-test-") as root:
            receipt_path = Path(root) / "receipt.json"
            M.write_new(receipt_path, M.canonical(self.fake_receipt()))
            output = Path(root) / "bundle"
            with mock.patch.object(M, "OUTPUT", output), \
                    mock.patch.object(M, "capture_worker") as capture:
                with self.assertRaises(TimeoutError):
                    M.run(receipt_path, time.monotonic() - 11)
                capture.assert_not_called()
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
