#!/usr/bin/env python3
"""Regression and tamper tests for the sealed P32 timing harness."""

from __future__ import annotations

import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

import wh2_p32_dispatch_timing as timing


def thermal_bytes(times=(99.5, 100.5, 101.5)):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=timing.THERMAL_FIELDS,
                            lineterminator="\n")
    writer.writeheader()
    for index, monotonic in enumerate(times):
        row = {field: "1" for field in timing.THERMAL_FIELDS}
        row.update({
            "utc": "2026-07-19T00:00:%02d.000Z" % index,
            "monotonic_s": "%.6f" % monotonic,
            "cpu_busy_pct": "100.0", "cpu_avg_mhz": "4500.0",
            "cpu_tctl_c": "60.0", "dimm_read_errors": "0",
            "load1": "128.0", "load5": "128.0", "load15": "128.0",
            "edac_ce": "0", "edac_ue": "0",
        })
        for field in timing.DIMM_FIELDS:
            row[field] = "45.0"
        writer.writerow(row)
    return stream.getvalue().encode("ascii")


class CanonicalReceiptTests(unittest.TestCase):
    def test_sealed_record_detects_tamper(self):
        value = timing.sealed_record("test.schema", {"a": 1})
        self.assertEqual(timing.verify_sealed(value, "test.schema", "test"), value)
        value["a"] = 2
        with self.assertRaisesRegex(timing.TimingError, "self hash mismatch"):
            timing.verify_sealed(value, "test.schema", "test")

    def test_task_cartesian_product_is_exact_and_unique(self):
        tasks = timing.generate_tasks()
        self.assertEqual(len(tasks), 13 * 6 * 3 * 3 * 2)
        self.assertEqual([task["job"] for task in tasks], list(range(1404)))
        self.assertEqual(len({task["task_id"] for task in tasks}), 1404)
        coordinates = {(task["K"], task["bb"], task["seed_index"],
                        task["schedule"], task["cache_state"])
                       for task in tasks}
        self.assertEqual(len(coordinates), 1404)

    def test_automatic_route_boundaries(self):
        self.assertEqual(timing.expected_rhs_route("p32_r7", 3199, 4096),
                         "streamed")
        self.assertEqual(timing.expected_rhs_route("p32_r7", 3200, 4096),
                         "joint-delta")
        self.assertEqual(timing.expected_rhs_route("p32_r7", 9999, 1280),
                         "streamed")
        self.assertEqual(timing.expected_rhs_route("p32_r7", 10000, 1280),
                         "joint-delta")
        self.assertEqual(timing.expected_rhs_route("prod244", 10000, 4096),
                         "streamed")

    def test_child_environment_is_explicit_and_sanitized(self):
        previous = os.environ.get("ANTHROPIC_API_KEY")
        os.environ["ANTHROPIC_API_KEY"] = "must-not-leak"
        try:
            environment = timing.sanitized_environment(Path("/empty"), allocator=True)
        finally:
            if previous is None:
                del os.environ["ANTHROPIC_API_KEY"]
            else:
                os.environ["ANTHROPIC_API_KEY"] = previous
        self.assertEqual(set(environment), {
            "HOME", "PATH", "LC_ALL", "LANG", "TZ", "PYTHONDONTWRITEBYTECODE",
            "MALLOC_MMAP_THRESHOLD_", "MALLOC_TRIM_THRESHOLD_"})
        self.assertNotIn("ANTHROPIC_API_KEY", environment)


class BoundedChildTests(unittest.TestCase):
    def test_normal_child_output(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        rc, stdout, stderr = timing.run_bounded(
            (sys.executable, "-c", "print('ok')"), environment, 5,
            stdout_limit=64, stderr_limit=64)
        self.assertEqual((rc, stdout, stderr), (0, b"ok\n", b""))

    def test_oversized_stdout_is_killed_at_bound(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with self.assertRaisesRegex(timing.BoundedProcessError, "stdout-limit"):
            timing.run_bounded(
                (sys.executable, "-c", "import os; os.write(1, b'x' * 65536)"),
                environment, 5, stdout_limit=1024, stderr_limit=64)

    def test_timeout_is_bounded(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with self.assertRaisesRegex(timing.BoundedProcessError, "timeout"):
            timing.run_bounded(
                (sys.executable, "-c", "import time; time.sleep(60)"),
                environment, 0.05, stdout_limit=64, stderr_limit=64)

    def test_timeout_still_applies_after_child_closes_pipes(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        script = "import os,time; os.close(1); os.close(2); time.sleep(60)"
        with self.assertRaisesRegex(timing.BoundedProcessError, "timeout"):
            timing.run_bounded(
                (sys.executable, "-c", script), environment, 0.05,
                stdout_limit=64, stderr_limit=64)

    def test_timeout_still_applies_when_descendant_holds_pipes(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        script = ("import os,time; child=os.fork(); "
                  "os._exit(0) if child else time.sleep(60)")
        with self.assertRaisesRegex(timing.BoundedProcessError, "timeout"):
            timing.run_bounded(
                (sys.executable, "-c", script), environment, 0.05,
                stdout_limit=64, stderr_limit=64)


class ProcessIdentityTests(unittest.TestCase):
    def make_fixture(self, root: Path):
        proc = root / "proc"
        process = proc / "123"
        process.mkdir(parents=True)
        # fields[19] after the closing comm parenthesis is Linux stat field 22.
        rest = [b"S"] + [b"0"] * 18 + [b"4242"]
        (process / "stat").write_bytes(b"123 (thermal sampler) " + b" ".join(rest) + b"\n")
        (process / "cmdline").write_bytes(
            b"/usr/bin/python3\0/frozen/sampler.py\0--csv\0/tmp/t.csv\0")
        (process / "status").write_text(
            "Name:\tpython3\nCpus_allowed_list:\t127\n", encoding="ascii")
        boot = root / "boot_id"
        boot.write_text("01234567-89ab-cdef-0123-456789abcdef\n", encoding="ascii")
        csv_path = root / "thermal.csv"
        csv_path.write_bytes(thermal_bytes())
        return proc, boot, csv_path

    def test_pid_reuse_start_tick_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            proc, boot, csv_path = self.make_fixture(root)
            identity = timing.capture_process_identity(
                123, 127, csv_path, proc_root=proc, boot_id_path=boot)
            self.assertTrue(timing.process_identity_matches(
                identity, 127, csv_path, proc_root=proc, boot_id_path=boot))
            stat_path = proc / "123/stat"
            stat_path.write_bytes(stat_path.read_bytes().replace(b"4242", b"4243"))
            self.assertFalse(timing.process_identity_matches(
                identity, 127, csv_path, proc_root=proc, boot_id_path=boot))

    def test_csv_inode_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            proc, boot, csv_path = self.make_fixture(root)
            identity = timing.capture_process_identity(
                123, 127, csv_path, proc_root=proc, boot_id_path=boot)
            replacement = root / "replacement.csv"
            replacement.write_bytes(thermal_bytes())
            os.replace(replacement, csv_path)
            self.assertFalse(timing.process_identity_matches(
                identity, 127, csv_path, proc_root=proc, boot_id_path=boot))

    def test_public_cmdline_redacts_ephemeral_pidfile(self):
        identity = {
            "pid": 123, "cmdline": ["python3", "sampler.py", "--pid-file",
                                      "/tmp/bootstrap.pid"],
        }
        public = timing._public_sampler_identity(identity)
        self.assertEqual(public["cmdline"][-1], "<ephemeral-bootstrap-only>")
        self.assertNotIn("bootstrap.pid", json.dumps(public))


class ThermalSealingTests(unittest.TestCase):
    def test_post_end_coverage_and_health(self):
        summary = timing.validate_thermal_interval(thermal_bytes(), 100.0, 101.0)
        self.assertEqual(summary["sample_count"], 3)
        self.assertGreaterEqual(summary["post_end_margin_s"], 0.0)

    def test_missing_post_end_sample_is_rejected(self):
        with self.assertRaisesRegex(timing.TimingError, "bracket benchmark end"):
            timing.validate_thermal_interval(
                thermal_bytes((99.5, 100.5)), 100.0, 101.0)

    def test_csv_append_after_terminal_receipt_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            original = thermal_bytes()
            path.write_bytes(original)
            stat = path.stat()
            receipt = {
                "thermal_csv_sha256": timing.sha256_bytes(original),
                "thermal_csv_size": len(original),
                "thermal_csv_device": stat.st_dev,
                "thermal_csv_inode": stat.st_ino,
                "timing_start_monotonic_s": 100.0,
                "benchmark_end_monotonic_s": 101.0,
                "thermal_summary": timing.validate_thermal_interval(
                    original, 100.0, 101.0),
            }
            timing.verify_terminal_thermal(path, receipt)
            with path.open("ab") as stream:
                stream.write(thermal_bytes((102.5,)).splitlines(keepends=True)[1])
            with self.assertRaisesRegex(timing.TimingError,
                                        "terminal thermal receipt mismatch"):
                timing.verify_terminal_thermal(path, receipt)

    def test_intermediate_edac_change_is_rejected(self):
        raw = thermal_bytes()
        lines = raw.decode("ascii").splitlines()
        reader = list(csv.DictReader(io.StringIO("\n".join(lines) + "\n")))
        reader[1]["edac_ce"] = "1"
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=timing.THERMAL_FIELDS,
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(reader)
        with self.assertRaisesRegex(timing.TimingError, "EDAC counters changed"):
            timing.validate_thermal_interval(
                stream.getvalue().encode("ascii"), 100.0, 101.0)


class PidfdStopTests(unittest.TestCase):
    def test_pidfd_stop_validates_identity_before_signal(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        process = subprocess.Popen(
            (sys.executable, "-c", "import time; time.sleep(60)"),
            env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            deadline = time.monotonic() + 5.0
            while True:
                cmdline = Path("/proc/%d/cmdline" % process.pid).read_bytes()
                if b"time.sleep(60)" in cmdline:
                    break
                if time.monotonic() >= deadline:
                    self.fail("test child did not finish exec")
                time.sleep(0.005)
            stat_raw = Path("/proc/%d/stat" % process.pid).read_bytes()
            tick = timing._process_start_tick(stat_raw)
            rc, stdout, stderr = timing.run_bounded(
                (sys.executable, "-c", timing.PIDFD_STOP_PROGRAM,
                 str(process.pid), str(tick), timing.sha256_bytes(cmdline)),
                environment, 5.0, 1024, 1024)
            self.assertEqual((rc, stdout, stderr), (0, b"", b""))
            self.assertEqual(process.wait(timeout=5), -15)
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=5)

    def test_pidfd_stop_rejects_wrong_cmdline_hash(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        process = subprocess.Popen(
            (sys.executable, "-c", "import time; time.sleep(60)"),
            env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            tick = timing._process_start_tick(
                Path("/proc/%d/stat" % process.pid).read_bytes())
            rc, _stdout, _stderr = timing.run_bounded(
                (sys.executable, "-c", timing.PIDFD_STOP_PROGRAM,
                 str(process.pid), str(tick), "0" * 64),
                environment, 5.0, 1024, 1024)
            self.assertEqual(rc, 73)
            self.assertIsNone(process.poll())
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=5)


class SegmentReceiptTests(unittest.TestCase):
    def make_fixture(self, root: Path):
        (root / "segments").mkdir()
        (root / "interrupted").mkdir()
        for name in ("design.json", "prepare_receipt.json",
                     "prelaunch_receipt.json", "tasks_manifest.jsonl"):
            (root / name).write_bytes((name + "\n").encode("ascii"))
        design = {
            "controller_core": 126, "thermal_core": 127, "task_count": 1,
            "tools": {"python3": {"path": "/usr/bin/python3"}},
        }
        raw = thermal_bytes()
        thermal = root / "segments/segment0000.thermal.csv"
        thermal.write_bytes(raw)
        stat = thermal.stat()
        identity = {
            "pid": 1234, "start_tick": 5678,
            "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
            "cmdline": timing._sampler_cmdline(root, design, 0, public=True),
            "cmdline_sha256": "a" * 64, "affinity": [127],
            "csv_device": stat.st_dev, "csv_inode": stat.st_ino,
        }
        terminal = timing.sealed_record(
            "wirehair.wh2.p32_dispatch.segment_terminal.v2", {
                "segment": 0, "status": "complete",
                "started_utc": "2026-07-19T00:00:00.000Z",
                "ended_utc": "2026-07-19T00:00:01.000Z",
                "timing_start_monotonic_s": 100.0,
                "benchmark_end_monotonic_s": 101.0, "duration_s": 1.0,
                "controller_pid": 4321, "controller_affinity": [126],
                "design_sha256": timing.sha256_file(root / "design.json"),
                "prepare_receipt_sha256": timing.sha256_file(
                    root / "prepare_receipt.json"),
                "prelaunch_receipt_sha256": timing.sha256_file(
                    root / "prelaunch_receipt.json"),
                "tasks_manifest_sha256": timing.sha256_file(
                    root / "tasks_manifest.jsonl"),
                "resumed_jobs_before": [],
                "completed_tasks": [{"job": 0, "receipt_sha256": "b" * 64,
                                     "retry_count": 0}],
                "recovered_interrupted_transactions": [], "failure": None,
                "sampler_identity_start": identity,
                "sampler_identity_end": identity,
                "graceful_stop": True,
                "graceful_stop_mechanism":
                    "sudo-python-pidfd-identity-verified-sigterm",
                "post_end_sample": True,
                "thermal_csv_name": thermal.name,
                "thermal_csv_sha256": timing.sha256_bytes(raw),
                "thermal_csv_size": len(raw),
                "thermal_csv_device": stat.st_dev,
                "thermal_csv_inode": stat.st_ino,
                "thermal_summary": timing.validate_thermal_interval(
                    raw, 100.0, 101.0),
            })
        terminal_path = root / "segments/segment0000.terminal.json"
        terminal_path.write_bytes(timing.canonical_json(terminal))
        entries = {0: {"receipt_sha256": "b" * 64,
                       "receipt": {"attempt": 0}}}
        return design, entries, terminal, terminal_path

    @staticmethod
    def reseal(terminal):
        payload = dict(terminal)
        schema = payload.pop("schema")
        payload.pop("self_sha256_excluding_field")
        return timing.sealed_record(schema, payload)

    def test_exact_terminal_receipt_replays(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, _terminal, _path = self.make_fixture(root)
            records = timing.verify_segments(root, design, entries)
            self.assertEqual(len(records), 1)

    def test_terminal_resumed_set_must_match_prior_segments(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, terminal, path = self.make_fixture(root)
            terminal["resumed_jobs_before"] = [0]
            path.write_bytes(timing.canonical_json(self.reseal(terminal)))
            with self.assertRaisesRegex(timing.TimingError,
                                        "segment task ledger malformed"):
                timing.verify_segments(root, design, entries)

    def test_ephemeral_pidfile_never_enters_terminal_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, _terminal, _path = self.make_fixture(root)
            (root / "segments/.segment0000.bootstrap.pid").write_text(
                "1234\n", encoding="ascii")
            with self.assertRaisesRegex(timing.TimingError,
                                        "ephemeral or unknown"):
                timing.verify_segments(root, design, entries)


class TransactionTests(unittest.TestCase):
    def make_root(self, root: Path):
        for name in ("ledger", ".transactions", "interrupted"):
            (root / name).mkdir()

    def test_interrupted_task_is_archived_not_completed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_root(root)
            task = timing.generate_tasks()[0]
            staging = timing.begin_task_transaction(root, task)
            (staging / "partial.stdout").write_bytes(b"partial")
            recovered = timing.recover_interrupted_transactions(root)
            self.assertEqual(len(recovered), 1)
            self.assertFalse(staging.exists())
            self.assertFalse((root / "ledger/0000").exists())
            timing._verify_interrupted_archive(root, recovered)

    def test_atomic_commit_has_exact_immutable_file_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_root(root)
            task = timing.generate_tasks()[0]
            receipt = timing.sealed_record(
                "wirehair.wh2.p32_dispatch.task.v2", {"job": 0})
            final = timing.commit_task_transaction(
                root, task, b"stdout\n", b"", receipt)
            self.assertEqual(sorted(path.name for path in final.iterdir()),
                             ["intent.json", "receipt.json", "stderr.bin",
                              "stdout.csv"])
            self.assertFalse(any((root / ".transactions").iterdir()))

    def test_unbound_completed_task_is_archived_for_rerun(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_root(root)
            task = timing.generate_tasks()[0]
            receipt = timing.sealed_record(
                "wirehair.wh2.p32_dispatch.task.v2", {"job": 0})
            timing.commit_task_transaction(root, task, b"stdout\n", b"", receipt)
            recovered = timing.recover_unbound_task_entries(root, {})
            self.assertEqual(len(recovered), 1)
            self.assertFalse((root / "ledger/0000").exists())
            timing._verify_interrupted_archive(root, recovered)

    def test_symlink_transaction_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_root(root)
            (root / ".transactions/0000.part").symlink_to(root / "ledger")
            with self.assertRaisesRegex(timing.TimingError, "unsafe"):
                timing.recover_interrupted_transactions(root)


if __name__ == "__main__":
    unittest.main()
