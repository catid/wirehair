#!/usr/bin/env python3
"""Regression tests for the WH2 thermal reader janitor.

The scenarios mirror the 2026-07-18/2026-07-21 incidents: a root-owned
sampler that a non-root ``kill -0`` misreported as dead, duplicate I2C
readers on the same host, orphaned load fillers, and PID reuse racing a
cleanup signal.
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest

import wh2_thermal_reader_janitor as janitor

SAMPLER_CMD = (
    "/usr/bin/python3",
    "/tmp/frozen/wirehair_expo_thermal_sampler.py",
    "--csv", "/tmp/wirehair-thermal.csv",
    "--pid-file", "/tmp/wirehair-thermal.pid",
)
FILLER_CMD = ("/usr/bin/bash", "-c", "while :; do :; done")


def completed(returncode: int) -> subprocess.CompletedProcess:
    return subprocess.CompletedProcess(args=(), returncode=returncode)


class FakeProc:
    """A /proc facsimile with root-owned process semantics."""

    def __init__(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.root = Path(self._directory.name)
        self.root_owned: set[int] = set()
        self.sudo_kill_log: list[tuple[str, ...]] = []
        self.sudo_available = True

    def cleanup(self) -> None:
        self._directory.cleanup()

    def add(self, pid: int, cmdline: tuple[str, ...], *,
            start_ticks: int, uid: int = 0) -> None:
        entry = self.root / str(pid)
        entry.mkdir()
        (entry / "cmdline").write_bytes(
            b"\0".join(token.encode() for token in cmdline) + b"\0")
        (entry / "stat").write_text(
            "%d (comm with (parens) and spaces) S 1 1 1 0 -1 0 0 0 0 0 "
            "0 0 0 0 20 0 1 0 %d 0 0" % (pid, start_ticks),
            encoding="ascii")
        (entry / "status").write_text(
            "Name:\tfake\nUid:\t%d\t%d\t%d\t%d\n" % (uid, uid, uid, uid),
            encoding="ascii")
        if uid == 0:
            self.root_owned.add(pid)

    def remove(self, pid: int) -> None:
        entry = self.root / str(pid)
        for child in entry.iterdir():
            child.unlink()
        entry.rmdir()
        self.root_owned.discard(pid)

    def kill_fn(self, pid: int, signal_number: int) -> None:
        assert signal_number == 0
        if not (self.root / str(pid)).exists():
            raise ProcessLookupError(pid)
        if pid in self.root_owned:
            raise PermissionError(pid)

    def run_fn(self, command, **_ignored) -> subprocess.CompletedProcess:
        command = tuple(command)
        if command == ("sudo", "-n", "true"):
            return completed(0 if self.sudo_available else 1)
        assert command[:3] == ("sudo", "-n", "kill")
        if not self.sudo_available:
            return completed(1)
        self.sudo_kill_log.append(command)
        pid = int(command[-1])
        if command[3] == "-0":
            return completed(0 if (self.root / str(pid)).exists() else 1)
        # Delivered TERM/KILL: the fake process exits immediately.
        if (self.root / str(pid)).exists():
            self.remove(pid)
            return completed(0)
        return completed(1)

    def affinity_fn(self, pid: int) -> set:
        if not (self.root / str(pid)).exists():
            raise OSError("no such process")
        return {127}


class FakeProcTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.proc = FakeProc()
        self.addCleanup(self.proc.cleanup)

    def run_janitor(self, mode: str, **overrides):
        parameters = dict(
            proc_root=self.proc.root,
            kill_fn=self.proc.kill_fn,
            run_fn=self.proc.run_fn,
            affinity_fn=self.proc.affinity_fn,
            sleep_fn=lambda _seconds: None,
            monotonic_fn=lambda: 0.0,
        )
        parameters.update(overrides)
        return janitor.run_janitor(mode, **parameters)


class LivenessTest(FakeProcTestCase):
    def test_root_owned_stale_pid_is_reported_alive(self) -> None:
        """The original bug: EPERM from kill -0 was misread as dead."""
        self.proc.add(3028028, SAMPLER_CMD, start_ticks=100, uid=0)
        state = janitor.pid_liveness(
            3028028, kill_fn=self.proc.kill_fn, run_fn=self.proc.run_fn,
            proc_root=self.proc.root)
        self.assertEqual(state, "alive")

    def test_exited_pid_is_dead(self) -> None:
        state = janitor.pid_liveness(
            424242, kill_fn=self.proc.kill_fn, run_fn=self.proc.run_fn,
            proc_root=self.proc.root)
        self.assertEqual(state, "dead")

    def test_eperm_without_sudo_answer_fails_closed(self) -> None:
        """EPERM plus an unanswerable sudo must never look dead."""
        self.proc.add(555, SAMPLER_CMD, start_ticks=7, uid=0)
        self.proc.sudo_available = False
        state = janitor.pid_liveness(
            555, kill_fn=self.proc.kill_fn, run_fn=self.proc.run_fn,
            proc_root=self.proc.root)
        self.assertEqual(state, "unknown")


class ScanTest(FakeProcTestCase):
    def test_scan_classifies_and_receipts_identity(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.add(200, FILLER_CMD, start_ticks=22, uid=0)
        self.proc.add(300, ("/usr/bin/vim", "notes.txt"),
                      start_ticks=33, uid=1000)
        records = janitor.scan_processes(
            proc_root=self.proc.root, affinity_fn=self.proc.affinity_fn)
        self.assertEqual(
            [(r["pid"], r["role"]) for r in records],
            [(100, "reader"), (200, "filler")])
        reader = records[0]
        self.assertEqual(reader["uid"], 0)
        self.assertEqual(reader["start_ticks"], 11)
        self.assertEqual(reader["affinity"], (127,))
        self.assertEqual(reader["cmdline"], list(SAMPLER_CMD))


class EnforceTest(FakeProcTestCase):
    def write_pid_file(self, pid: int) -> Path:
        path = self.proc.root / "expected.pid"
        path.write_text("%d\n" % pid, encoding="ascii")
        return path

    def test_sole_expected_reader_passes(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        receipt = self.run_janitor(
            "enforce", expect_pid_file=self.write_pid_file(100))
        self.assertEqual(receipt["violations"], [])
        self.assertEqual(receipt["expected_pid"], 100)

    def test_duplicate_reader_is_rejected(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.add(101, SAMPLER_CMD, start_ticks=12, uid=0)
        receipt = self.run_janitor(
            "enforce", expect_pid_file=self.write_pid_file(100))
        self.assertEqual(
            receipt["violations"],
            ["unexpected extra I2C reader PID 101"])

    def test_unexpected_reader_without_expectation_is_rejected(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        receipt = self.run_janitor("enforce")
        self.assertEqual(
            receipt["violations"],
            ["no I2C reader expected but PID 100 is live"])

    def test_stale_expected_pid_is_rejected(self) -> None:
        receipt = self.run_janitor(
            "enforce", expect_pid_file=self.write_pid_file(100))
        self.assertEqual(
            receipt["violations"],
            ["expected sole I2C reader PID 100 is not a live sampler"])

    def test_lingering_fillers_are_rejected(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.add(200, FILLER_CMD, start_ticks=22, uid=0)
        receipt = self.run_janitor(
            "enforce", expect_pid_file=self.write_pid_file(100))
        self.assertEqual(
            receipt["violations"],
            ["load filler PID 200 is still present"])

    def test_noncanonical_pid_file_is_rejected(self) -> None:
        path = self.proc.root / "expected.pid"
        path.write_text("0100\n", encoding="ascii")
        with self.assertRaises(janitor.JanitorError):
            self.run_janitor("enforce", expect_pid_file=path)

    def test_unclassifiable_protected_pid_fails_closed(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.sudo_available = False
        receipt = self.run_janitor(
            "enforce", expect_pid_file=self.write_pid_file(100))
        self.assertIn(
            "cannot determine liveness of protected PID 100",
            receipt["violations"])


class StopOrphansTest(FakeProcTestCase):
    def test_orphans_are_terminated_through_sudo(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.add(200, FILLER_CMD, start_ticks=22, uid=0)
        receipt = self.run_janitor("stop-orphans")
        self.assertEqual(receipt["violations"], [])
        self.assertEqual(
            [(a["pid"], a["outcome"]) for a in receipt["actions"]],
            [(100, "stopped"), (200, "stopped")])
        self.assertIn(
            ("sudo", "-n", "kill", "-TERM", "100"),
            self.proc.sudo_kill_log)

    def test_reused_pid_is_never_signalled(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        record = janitor.scan_processes(
            proc_root=self.proc.root,
            affinity_fn=self.proc.affinity_fn)[0]
        self.proc.remove(100)
        self.proc.add(100, SAMPLER_CMD, start_ticks=99, uid=0)
        action = janitor.stop_process(
            record, proc_root=self.proc.root, run_fn=self.proc.run_fn,
            sleep_fn=lambda _seconds: None, monotonic_fn=lambda: 0.0)
        self.assertEqual(action["outcome"], "pid-reused")
        self.assertEqual(action["signals"], [])
        self.assertEqual(self.proc.sudo_kill_log, [])

    def test_without_sudo_orphans_are_reported_not_leaked(self) -> None:
        self.proc.add(100, SAMPLER_CMD, start_ticks=11, uid=0)
        self.proc.sudo_available = False
        receipt = self.run_janitor("stop-orphans")
        self.assertEqual(receipt["actions"], [])
        self.assertIn(
            "sudo -n is unavailable; cannot stop root-owned orphans",
            receipt["violations"])


class RealProcSmokeTest(unittest.TestCase):
    def test_audit_scan_of_real_proc_is_read_only_and_clean(self) -> None:
        records = janitor.scan_processes()
        for record in records:
            self.assertIn(record["role"], ("reader", "filler"))


if __name__ == "__main__":
    unittest.main()
