#!/usr/bin/env python3
"""Regression and tamper tests for the sealed P32 timing harness."""

from __future__ import annotations

import csv
from contextlib import nullcontext
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

import wh2_p32_dispatch_timing as timing
import wirehair_expo_thermal_sampler as thermal_sampler


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


def root_owned_sealed_stat(path):
    actual = path.lstat()
    return mock.Mock(
        st_mode=(actual.st_mode & ~0o777) | 0o444,
        st_uid=0, st_nlink=1, st_dev=actual.st_dev, st_ino=actual.st_ino,
        st_size=actual.st_size)


class CanonicalReceiptTests(unittest.TestCase):
    @staticmethod
    def root_receipts():
        campaign = timing.sealed_record(
            "wirehair.wh2.p32_dispatch.campaign.v2", {
                "completed_utc": "2026-07-20T00:00:00.000Z",
                "design_sha256": "a" * 64,
                "prepare_receipt_sha256": "b" * 64,
                "prelaunch_receipt_sha256": "c" * 64,
                "tasks_manifest_sha256": "d" * 64,
                "task_count": 0, "task_receipts": [],
                "terminal_receipts": [],
            })
        summary = timing.sealed_record(
            "wirehair.wh2.p32_dispatch.summary.v2", {
                "created_utc": "2026-07-20T00:00:01.000Z",
                "design_sha256": "a" * 64,
                "campaign_receipt_sha256": "e" * 64,
                "task_count": 0, "segment_count": 0,
                "comparison": "p32_r7_vs_prod244",
                "timing_evidence_promotional": True,
                "architecture_promotion_ready": False,
                "architecture_promotion_blocker":
                    "independent-cross-payload-recovery-required",
                "dispatch_speed_policy": timing.DISPATCH_SPEED_POLICY,
                "seed_fixes": "none", "aggregates": {},
            })
        return {
            "campaign_receipt.json": campaign,
            "validated_summary.json": summary,
        }

    def test_root_atomic_receipts_recover_before_and_after_link_crashes(self):
        for name, receipt in self.root_receipts().items():
            raw = timing.canonical_json(receipt)
            for state in ("before-link", "after-link"):
                with self.subTest(name=name, state=state), \
                        tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    part = root / (name + ".part")
                    part.write_bytes(raw)
                    if state == "after-link":
                        (root / name).write_bytes(raw)
                    timing.recover_root_atomic_receipts(root)
                    self.assertEqual((root / name).read_bytes(), raw)
                    self.assertFalse(part.exists())

    def test_root_atomic_receipt_mismatched_link_remnant_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = self.root_receipts()["campaign_receipt.json"]
            (root / "campaign_receipt.json").write_bytes(
                timing.canonical_json(receipt))
            (root / "campaign_receipt.json.part").write_bytes(b"changed\n")
            with self.assertRaisesRegex(
                    timing.TimingError, "link remnant changed"):
                timing.recover_root_atomic_receipts(root)

    def test_campaign_lock_rejects_concurrent_recovery_or_verification(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock = root / "controller.lock"
            lock.write_bytes(b"")
            design = {"controller_lock_sha256": timing.sha256_bytes(b"")}
            with lock.open("rb") as owner:
                timing.fcntl.flock(
                    owner.fileno(), timing.fcntl.LOCK_EX |
                    timing.fcntl.LOCK_NB)
                with mock.patch.object(
                        timing, "load_design", return_value=design), \
                        self.assertRaisesRegex(
                            timing.TimingError, "another campaign controller"):
                    with timing._exclusive_campaign_lock(root):
                        self.fail("concurrent lock unexpectedly acquired")

    def test_eviction_receipt_is_bound_to_exact_command_value(self):
        timing._validate_evict_bytes_receipt("4096", 4096)
        for value, expected in (("8192", 4096), ("4096", 8192),
                                ("04096", 4096), ("4096", True),
                                ("4096", 4095)):
            with self.subTest(value=value, expected=expected), \
                    self.assertRaises(timing.TimingError):
                timing._validate_evict_bytes_receipt(value, expected)

    def test_seed_and_result_receipts_are_bounded_to_cpp_domains(self):
        self.assertEqual(timing.parse_hex_uint(
            "0xffffffffffffffff", "matrix seed", (1 << 64) - 1),
            (1 << 64) - 1)
        self.assertEqual(timing.parse_hex_uint(
            "0xffffffff", "peel seed", (1 << 32) - 1), (1 << 32) - 1)
        self.assertEqual(timing.parse_sint(
            str(-(1 << 31)), "result", -(1 << 31), (1 << 31) - 1),
            -(1 << 31))
        for value, maximum in (("0x10000000000000000", (1 << 64) - 1),
                               ("0x100000000", (1 << 32) - 1),
                               ("0x00", (1 << 32) - 1)):
            with self.subTest(value=value), self.assertRaises(timing.TimingError):
                timing.parse_hex_uint(value, "seed", maximum)
        for value in (str(-(1 << 31) - 1), str(1 << 31)):
            with self.subTest(value=value), self.assertRaises(timing.TimingError):
                timing.parse_sint(
                    value, "result", -(1 << 31), (1 << 31) - 1)

    def test_runtime_tool_receipt_is_exact_and_executable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tool = root / "tool"
            tool.write_bytes(b"tool\n")
            tool.chmod(0o755)
            receipt = {
                "path": str(tool.resolve()),
                "sha256": timing.sha256_file(tool),
            }
            timing._verify_runtime_tool_receipt("tool", receipt)
            invalid = (
                {**receipt, "extra": True},
                {**receipt, "path": tool},
                {**receipt, "path": "tool"},
                {**receipt, "sha256": receipt["sha256"].upper()},
            )
            for mutation in invalid:
                with self.subTest(mutation=mutation), self.assertRaises(
                        timing.TimingError):
                    timing._verify_runtime_tool_receipt("tool", mutation)
            link = root / "link"
            link.symlink_to(tool)
            with self.assertRaises(timing.TimingError):
                timing._verify_runtime_tool_receipt(
                    "tool", {**receipt, "path": str(link)})
            tool.chmod(0o644)
            with self.assertRaises(timing.TimingError):
                timing._verify_runtime_tool_receipt("tool", receipt)

    @staticmethod
    def make_speed_entries():
        work_fields = (
            "inactivated", "binary_def", "heavy_gain", "block_xors",
            "block_muladds", "joint_source_xors", "joint_marginal_xors",
            "joint_marginal_copies", "joint_active_deltas",
            "joint_scratch_bytes", "dual_source_columns",
            "intermediate_bytes",
        )
        tasks = []
        entries = {}
        for seed_index, seed in enumerate(timing.SEEDS):
            for schedule in timing.SCHEDULES:
                trace = "%d:%s" % (seed_index, schedule)
                for cache_state in timing.CACHE_STATES:
                    job = len(tasks)
                    task = {
                        "job": job, "K": 3200, "bb": 64,
                        "seed_index": seed_index, "seed": seed,
                        "schedule": schedule, "cache_state": cache_state,
                    }
                    rows = tuple({"arm": arm, **{field: "1" for field in work_fields}}
                                 for arm in ("control", "candidate"))
                    parsed = timing.ParsedOutput(
                        preamble={"candidate_preflight_rhs_route": "streamed"},
                        rows=rows, stdout_sha256="a" * 64,
                        trace_sha256=trace, cell_class="common-success",
                        common_success=True, timed_control_ns=100,
                        timed_candidate_ns=90,
                        work_signatures={"control": ("control",),
                                         "candidate": ("candidate",)},
                        contaminations=(),
                    )
                    tasks.append(task)
                    entries[job] = {"parsed": parsed}
        return tasks, entries

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

    def test_active_runner_must_match_frozen_runner_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "runner.py"
            source.write_bytes(b"exact frozen runner\n")
            design = {"immutable_files": {
                "frozen/wh2_p32_dispatch_timing.py":
                    timing.sha256_file(source),
            }}
            with mock.patch.object(timing, "__file__", str(source)):
                timing.verify_active_runner_identity(design)
                source.write_bytes(b"stale active runner\n")
                with self.assertRaisesRegex(
                        timing.TimingError, "active runner differs"):
                    timing.verify_active_runner_identity(design)

    def test_committed_source_blob_rejects_active_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            repo = Path(temporary)
            source = repo / "bench/runner.py"
            source.parent.mkdir()
            source.write_bytes(b"committed bytes\n")
            with mock.patch.object(
                    timing, "run_bounded",
                    return_value=(0, b"committed bytes\n", b"")):
                relative, blob = timing._committed_source_blob(
                    Path("/frozen/git"), repo, "a" * 40, source)
                self.assertEqual((relative, blob),
                                 ("bench/runner.py", b"committed bytes\n"))
                source.write_bytes(b"assume-unchanged mismatch\n")
                with self.assertRaisesRegex(
                        timing.TimingError, "differs from committed HEAD"):
                    timing._committed_source_blob(
                        Path("/frozen/git"), repo, "a" * 40, source)

    def test_controller_termination_signal_is_deferred_and_restored(self):
        previous = signal.getsignal(signal.SIGTERM)
        with self.assertRaisesRegex(timing.ControllerTermination, "SIGTERM"):
            with timing.DeferredTermination() as guard:
                os.kill(os.getpid(), signal.SIGTERM)
                deadline = time.monotonic() + 1.0
                while guard.error() is None and time.monotonic() < deadline:
                    time.sleep(0.001)
                with self.assertRaisesRegex(
                        timing.ControllerTermination, "SIGTERM"):
                    guard.raise_if_requested()
        self.assertIs(signal.getsignal(signal.SIGTERM), previous)

    def test_pending_termination_is_not_lost_at_guard_exit(self):
        guard = timing.DeferredTermination()
        with self.assertRaisesRegex(timing.ControllerTermination, "SIGTERM"):
            with guard:
                guard._record(signal.SIGTERM, None)

    def test_clustered_bootstrap_is_deterministic_and_pairs_cache_states(self):
        values = []
        for seed_index, timings in enumerate(((100, 90), (200, 205))):
            for cache_state in ("cold", "warm"):
                values.append({
                    "K": 3200, "bb": 64, "seed_index": seed_index,
                    "schedule": "burst", "cache_state": cache_state,
                    "control_ns": timings[0], "candidate_ns": timings[1],
                })
        first = timing._clustered_bootstrap(
            values, "unit-test", repetitions=200)
        second = timing._clustered_bootstrap(
            values, "unit-test", repetitions=200)
        self.assertEqual(first, second)
        self.assertEqual(first["cluster_count"], 2)
        self.assertLess(float(first["lower_95"]), float(first["upper_95"]))

    def test_dispatch_speed_gate_is_fail_closed(self):
        summary = {
            "cell_count": 18,
            "outcome_counts": {
                "common-success": 18, "control-only": 0,
                "candidate-only": 0, "common-failure": 0,
            },
            "ratio_of_sums": 0.999,
            "bootstrap": {"upper_one_sided_95": "1.005"},
        }
        self.assertTrue(timing._dispatch_speed_gate(summary)["speed_eligible"])
        for change in (
                {"ratio_of_sums": 1.0001},
                {"bootstrap": {"upper_one_sided_95": "1.010"}},
                {"bootstrap": None},
                {"outcome_counts": {
                    "common-success": 17, "control-only": 1,
                    "candidate-only": 0, "common-failure": 0,
                }}):
            candidate = dict(summary)
            candidate.update(change)
            self.assertFalse(
                timing._dispatch_speed_gate(candidate)["speed_eligible"])

    def test_reducer_separates_timing_evidence_from_architecture_promotion(self):
        tasks, entries = self.make_speed_entries()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "design.json").write_bytes(b"design\n")
            (root / "campaign_receipt.json").write_bytes(b"campaign\n")
            with mock.patch.object(timing, "KS", (3200,)), \
                    mock.patch.object(timing, "WIDTHS", (64,)), \
                    mock.patch.object(
                        timing, "_verify_campaign_data_locked",
                        return_value=({}, tasks, entries, [])), \
                    mock.patch.object(
                        timing, "_exclusive_campaign_lock",
                        return_value=nullcontext({})), \
                    mock.patch.object(timing, "recover_root_atomic_receipts"), \
                    mock.patch("sys.stdout", new=io.StringIO()):
                timing.reduce_campaign(mock.Mock(result_dir=str(root)))
            summary = timing.load_canonical(
                root / "validated_summary.json", "unit summary")
        self.assertIs(summary["timing_evidence_promotional"], True)
        self.assertIs(summary["architecture_promotion_ready"], False)
        self.assertEqual(
            summary["architecture_promotion_blocker"],
            "independent-cross-payload-recovery-required")
        self.assertTrue(
            summary["aggregates"]["overall_speed_gate"]["speed_eligible"])
        self.assertEqual(summary["aggregates"]["speed_eligible_K_bb"],
                         [{"K": 3200, "bb": 64}])


class BoundedChildTests(unittest.TestCase):
    @staticmethod
    def process_is_live(pid):
        try:
            state, _group, _session, _tick = \
                timing._process_stat_session_identity(
                    Path("/proc/%d/stat" % pid).read_bytes())
        except (OSError, timing.TimingError):
            return False
        return state != b"Z"

    def test_normal_child_output(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        rc, stdout, stderr = timing.run_bounded(
            (sys.executable, "-c", "print('ok')"), environment, 5,
            stdout_limit=64, stderr_limit=64)
        self.assertEqual((rc, stdout, stderr), (0, b"ok\n", b""))

    def test_nonfinite_or_boolean_child_timeouts_are_rejected_before_launch(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        for timeout in (float("nan"), float("inf"), True, 0.0):
            with self.subTest(timeout=timeout), self.assertRaises(
                    timing.TimingError):
                timing.run_bounded(
                    (sys.executable, "-c", "raise SystemExit(99)"),
                    environment, timeout, stdout_limit=64, stderr_limit=64)
            with self.subTest(privileged_timeout=timeout), self.assertRaises(
                    timing.TimingError):
                timing.run_privileged_bounded(
                    Path("/sudo"), Path("/timeout"), ("/true",), environment,
                    helper_timeout_s=timeout)

    def test_exec_failure_preserves_original_error_without_signaling(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with mock.patch.object(
                timing, "_terminate_bounded_process_group_fallback") as fallback, \
                mock.patch.object(
                    timing, "_terminate_bounded_process_session") as session:
            with self.assertRaises(FileNotFoundError):
                timing.run_bounded(
                    ("/definitely/no-such-wh2-executable",), environment, 5,
                    stdout_limit=64, stderr_limit=64)
        fallback.assert_not_called()
        session.assert_not_called()

    def test_oversized_stdout_is_killed_at_bound(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with self.assertRaisesRegex(timing.BoundedProcessError, "stdout-limit"):
            timing.run_bounded(
                (sys.executable, "-c", "import os; os.write(1, b'x' * 65536)"),
                environment, 5, stdout_limit=1024, stderr_limit=64)

    def test_baseexception_reaps_detached_child(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pid"

            def interrupt(_selector, _timeout=None):
                deadline = time.monotonic() + 2.0
                while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                        time.monotonic() < deadline:
                    time.sleep(0.005)
                raise KeyboardInterrupt()

            script = (
                "import os,time\n"
                "child = os.fork()\n"
                "if child:\n"
                "    os._exit(0)\n"
                "open(%r, 'w').write(str(os.getpid()))\n"
                "time.sleep(60)\n" % str(pid_path))
            with mock.patch.object(
                    timing.selectors.DefaultSelector, "select", new=interrupt):
                with self.assertRaises(KeyboardInterrupt):
                    timing.run_bounded(
                        (sys.executable, "-c", script), environment, 30,
                        stdout_limit=64, stderr_limit=64)
            self.assertTrue(pid_path.is_file())
            pid = int(pid_path.read_text(encoding="ascii"))
            self.assertFalse(self.process_is_live(pid))

    def test_setup_baseexception_reaps_started_child(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pid"
            calls = 0

            def fail_second_set_blocking(_fd, _blocking):
                nonlocal calls
                calls += 1
                if calls == 2:
                    deadline = time.monotonic() + 2.0
                    while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                            time.monotonic() < deadline:
                        time.sleep(0.005)
                    raise KeyboardInterrupt()

            script = (
                "import os,time\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "os.fchmod(fd, 0o644)\n"
                "os.write(fd, str(os.getpid()).encode('ascii')); os.close(fd)\n"
                "time.sleep(60)\n" % str(pid_path))
            with mock.patch.object(
                    timing.os, "set_blocking",
                    side_effect=fail_second_set_blocking):
                with self.assertRaises(KeyboardInterrupt):
                    timing.run_bounded(
                        (sys.executable, "-c", script), environment, 30,
                        stdout_limit=64, stderr_limit=64)
            self.assertTrue(pid_path.is_file())
            pid = int(pid_path.read_text(encoding="ascii"))
            self.assertFalse(self.process_is_live(pid))

    def test_identity_capture_baseexception_reaps_whole_session(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pids"
            script = (
                "import os,time\n"
                "child=os.fork()\n"
                "if child: time.sleep(60)\n"
                "os.setpgid(0,0)\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "text='%%d %%d' %% (os.getppid(),os.getpid())\n"
                "os.write(fd, text.encode('ascii')); os.close(fd)\n"
                "time.sleep(60)\n" % str(pid_path))
            original = timing._process_stat_session_identity
            calls = 0

            def interrupt_first_capture(raw):
                nonlocal calls
                calls += 1
                if calls == 1:
                    deadline = time.monotonic() + 2.0
                    while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                            time.monotonic() < deadline:
                        time.sleep(0.005)
                    raise KeyboardInterrupt()
                return original(raw)

            pids = ()
            try:
                with mock.patch.object(
                        timing, "_process_stat_session_identity",
                        side_effect=interrupt_first_capture):
                    with self.assertRaises(KeyboardInterrupt):
                        timing.run_bounded(
                            (sys.executable, "-c", script), environment, 30,
                            stdout_limit=64, stderr_limit=64)
                self.assertTrue(pid_path.is_file())
                pids = tuple(map(
                    int, pid_path.read_text(encoding="ascii").split()))
                self.assertEqual(len(pids), 2)
                self.assertTrue(all(not self.process_is_live(pid) for pid in pids))
            finally:
                if not pids and pid_path.is_file() and pid_path.stat().st_size:
                    pids = tuple(map(
                        int, pid_path.read_text(encoding="ascii").split()))
                for pid in pids:
                    if self.process_is_live(pid):
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    def test_popen_constructor_baseexception_reaps_whole_session(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pids"
            script = (
                "import os,time\n"
                "child=os.fork()\n"
                "if child: time.sleep(60)\n"
                "os.setpgid(0,0)\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "text='%%d %%d' %% (os.getppid(),os.getpid())\n"
                "os.write(fd, text.encode('ascii')); os.close(fd)\n"
                "time.sleep(60)\n" % str(pid_path))
            original = timing.subprocess.Popen._execute_child

            def interrupt_after_fork(process, *args, **kwargs):
                original(process, *args, **kwargs)
                deadline = time.monotonic() + 2.0
                while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                        time.monotonic() < deadline:
                    time.sleep(0.005)
                raise KeyboardInterrupt()

            pids = ()
            try:
                with mock.patch.object(
                        timing.subprocess.Popen, "_execute_child",
                        new=interrupt_after_fork):
                    with self.assertRaises(KeyboardInterrupt):
                        timing.run_bounded(
                            (sys.executable, "-c", script), environment, 30,
                            stdout_limit=64, stderr_limit=64)
                self.assertTrue(pid_path.is_file())
                pids = tuple(map(
                    int, pid_path.read_text(encoding="ascii").split()))
                self.assertEqual(len(pids), 2)
                self.assertTrue(all(not self.process_is_live(pid) for pid in pids))
            finally:
                if not pids and pid_path.is_file() and pid_path.stat().st_size:
                    pids = tuple(map(
                        int, pid_path.read_text(encoding="ascii").split()))
                for pid in pids:
                    if self.process_is_live(pid):
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    def test_signal_after_kernel_fork_is_deferred_until_exact_cleanup(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pids"
            script = (
                "import os,time\n"
                "child=os.fork()\n"
                "if child: time.sleep(60)\n"
                "os.setpgid(0,0)\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "text='%%d %%d' %% (os.getppid(),os.getpid())\n"
                "os.write(fd, text.encode('ascii')); os.close(fd)\n"
                "time.sleep(60)\n" % str(pid_path))
            fork_owner = timing.subprocess if hasattr(
                timing.subprocess, "_fork_exec") else \
                timing.subprocess._posixsubprocess
            fork_name = "_fork_exec" if fork_owner is timing.subprocess else \
                "fork_exec"
            original = getattr(fork_owner, fork_name)

            def signal_after_fork(*args):
                pid = original(*args)
                deadline = time.monotonic() + 2.0
                while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                        time.monotonic() < deadline:
                    time.sleep(0.005)
                os.kill(os.getpid(), signal.SIGINT)
                return pid

            pids = ()
            try:
                with mock.patch.object(
                        fork_owner, fork_name,
                        new=signal_after_fork):
                    with self.assertRaisesRegex(
                            timing.ControllerTermination, "SIGINT"):
                        timing.run_bounded(
                            (sys.executable, "-c", script), environment, 30,
                            stdout_limit=64, stderr_limit=64)
                self.assertTrue(pid_path.is_file())
                pids = tuple(map(
                    int, pid_path.read_text(encoding="ascii").split()))
                self.assertEqual(len(pids), 2)
                self.assertTrue(all(not self.process_is_live(pid) for pid in pids))
            finally:
                if not pids and pid_path.is_file() and pid_path.stat().st_size:
                    pids = tuple(map(
                        int, pid_path.read_text(encoding="ascii").split()))
                for pid in pids:
                    if self.process_is_live(pid):
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    def test_successful_leader_cannot_leave_closed_pipe_descendant(self):
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "pid"
            script = (
                "import os,time\n"
                "child=os.fork()\n"
                "if child: os._exit(0)\n"
                "os.setpgid(0,0)\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "os.write(fd, str(os.getpid()).encode('ascii')); os.close(fd)\n"
                "os.close(1); os.close(2)\n"
                "time.sleep(60)\n" % str(pid_path))
            with self.assertRaisesRegex(
                    timing.BoundedProcessError, "descendant-survival"):
                timing.run_bounded(
                    (sys.executable, "-c", script), environment, 5,
                    stdout_limit=64, stderr_limit=64)
            self.assertTrue(pid_path.is_file())
            pid = int(pid_path.read_text(encoding="ascii"))
            self.assertFalse(Path("/proc/%d" % pid).exists())

    def test_privileged_timeout_cleans_root_child_on_controller_exception(self):
        sudo = Path("/usr/bin/sudo")
        timeout = Path("/usr/bin/timeout")
        if not sudo.is_file() or not timeout.is_file():
            self.skipTest("sudo/timeout unavailable")
        try:
            probe = subprocess.run(
                (str(sudo), "-n", "/usr/bin/true"),
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=2, check=False)
        except (OSError, subprocess.SubprocessError):
            self.skipTest("passwordless sudo unavailable")
        if probe.returncode != 0:
            self.skipTest("passwordless sudo unavailable")
        environment = timing.sanitized_environment(Path("/tmp"), allocator=False)
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "root-pid"
            script = (
                "import os,time\n"
                "child=os.fork()\n"
                "if child: time.sleep(60)\n"
                "os.setpgid(0,0)\n"
                "fd=os.open(%r, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0o600)\n"
                "os.fchmod(fd, 0o644)\n"
                "text='%%d %%d' %% (os.getppid(),os.getpid())\n"
                "os.write(fd, text.encode('ascii')); os.close(fd)\n"
                "os.close(1); os.close(2)\n"
                "time.sleep(60)\n" % str(pid_path))

            def interrupt():
                deadline = time.monotonic() + 2.0
                while (not pid_path.exists() or pid_path.stat().st_size == 0) and \
                        time.monotonic() < deadline:
                    time.sleep(0.005)
                raise KeyboardInterrupt()

            pids = ()
            try:
                with self.assertRaises(KeyboardInterrupt):
                    timing.run_privileged_bounded(
                        sudo, timeout, (sys.executable, "-c", script), environment,
                        helper_timeout_s=0.25, stdout_limit=64, stderr_limit=64,
                        poll_callback=interrupt)
                self.assertTrue(pid_path.is_file())
                pids = tuple(map(
                    int, pid_path.read_text(encoding="ascii").split()))
                self.assertEqual(len(pids), 2)
                self.assertTrue(all(not self.process_is_live(pid) for pid in pids))
            finally:
                if not pids and pid_path.is_file() and pid_path.stat().st_size:
                    pids = tuple(map(
                        int, pid_path.read_text(encoding="ascii").split()))
                live = tuple(pid for pid in pids if self.process_is_live(pid))
                if live:
                    subprocess.run(
                        (str(sudo), "-n", "/usr/bin/kill", "-KILL",
                         *(str(pid) for pid in live)),
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                        timeout=2, check=False)

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

    def test_partial_worktree_add_is_always_removed(self):
        for failure in ("nonzero", "exception"):
            with self.subTest(failure=failure), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                repo = root / "repo"
                repo.mkdir()
                workspace = root / "workspace"
                staging = root / "staging"
                staging.mkdir()
                source = workspace / "source"
                calls = []

                def bounded(argv, *_args, **_kwargs):
                    calls.append(tuple(argv))
                    if "add" in argv:
                        source.mkdir(parents=True)
                        if failure == "exception":
                            raise KeyboardInterrupt()
                        return 1, b"", b"partial add"
                    if "remove" in argv:
                        source.rmdir()
                        return 0, b"", b""
                    raise AssertionError("unexpected bounded command")

                expected = KeyboardInterrupt if failure == "exception" else \
                    timing.TimingError
                with mock.patch.object(
                        timing, "run_bounded", side_effect=bounded), \
                        mock.patch.object(
                            timing, "_git_worktree_registered",
                            return_value=False) as registered:
                    with self.assertRaises(expected):
                        timing._build_frozen_binary(
                            repo, "a" * 40, workspace, staging,
                            {"git": Path("/git")}, 1,
                            Path("/cc"), Path("/cxx"))
                self.assertTrue(any("add" in argv for argv in calls))
                self.assertTrue(any("remove" in argv for argv in calls))
                registered.assert_called_once_with(Path("/git"), repo, source)
                self.assertFalse(source.exists())

    def test_worktree_cleanup_failure_preserves_primary_exception(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo = root / "repo"
            repo.mkdir()
            workspace = root / "workspace"
            staging = root / "staging"
            staging.mkdir()
            source = workspace / "source"

            def bounded(argv, *_args, **_kwargs):
                if "add" in argv:
                    source.mkdir(parents=True)
                    raise KeyboardInterrupt("primary")
                if "remove" in argv:
                    raise timing.BoundedProcessError("cleanup timeout")
                raise AssertionError("unexpected bounded command")

            with mock.patch.object(
                    timing, "run_bounded", side_effect=bounded), \
                    mock.patch.object(
                        timing, "_git_worktree_registered", return_value=True):
                with self.assertRaisesRegex(
                        timing.TimingError,
                        "frozen build failed primary.*cleanup timeout") as caught:
                    timing._build_frozen_binary(
                        repo, "a" * 40, workspace, staging,
                        {"git": Path("/git")}, 1,
                        Path("/cc"), Path("/cxx"))
            self.assertIsInstance(caught.exception.__cause__, KeyboardInterrupt)

    def test_prepare_sudo_preflight_is_privileged_and_bounded(self):
        failures = ((1, b"", b""), (0, b"unexpected", b""),
                    (0, b"", b"unexpected"))
        for outcome in failures:
            with self.subTest(outcome=outcome), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                args = mock.Mock(
                    result_dir=str(root / "result"), repo=str(root),
                    build_jobs=1, core=0, controller_core=1, thermal_core=2,
                    numa_node=0, evict_bytes=4096,
                    c_compiler=None, cxx_compiler=None)
                with mock.patch.object(
                        timing, "resolve_tool",
                        side_effect=lambda name: Path("/usr/bin") / name), \
                        mock.patch.object(timing.shutil, "which", return_value=None), \
                        mock.patch.object(
                            timing, "run_privileged_bounded",
                            return_value=outcome) as bounded:
                    with self.assertRaisesRegex(
                            timing.TimingError,
                            "passwordless sudo preflight failed"):
                        timing.prepare_campaign(args)
                call = bounded.call_args
                self.assertEqual(call.args[:3], (
                    Path("/usr/bin/sudo"), Path("/usr/bin/timeout"),
                    ("/usr/bin/true",)))
                self.assertEqual(call.args[3]["LC_ALL"], "C")


class ThermalSamplerTests(unittest.TestCase):
    def test_evidence_open_preserves_primary_failure_when_close_also_fails(self):
        with mock.patch.object(thermal_sampler.os, "open", return_value=77), \
                mock.patch.object(thermal_sampler.os, "fchmod"), \
                mock.patch.object(
                    thermal_sampler.os, "fdopen",
                    side_effect=RuntimeError("fdopen-primary")), \
                mock.patch.object(
                    thermal_sampler.os, "close",
                    side_effect=OSError("close-secondary")) as close:
            with self.assertRaisesRegex(RuntimeError, "fdopen-primary"):
                thermal_sampler.open_exclusive_evidence("evidence")
        close.assert_called_once_with(77)

    def test_evidence_files_are_readonly_even_under_umask_zero(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence"
            old_umask = os.umask(0)
            try:
                with thermal_sampler.open_exclusive_evidence(path) as stream:
                    self.assertEqual(os.fstat(stream.fileno()).st_mode & 0o777,
                                     0o444)
                    stream.write("root-open-fd-remains-writable\n")
                    stream.flush()
                    os.fsync(stream.fileno())
            finally:
                os.umask(old_umask)
            self.assertEqual(path.stat().st_mode & 0o777, 0o444)
            self.assertEqual(path.read_text(encoding="ascii"),
                             "root-open-fd-remains-writable\n")

    def test_sampler_argument_bounds_reject_nan_and_infinity(self):
        thermal_sampler.validate_sampling_arguments(1.0, 5, 0.01)
        invalid = (
            (float("nan"), 5, 0.01), (float("inf"), 5, 0.01),
            (0.0, 5, 0.01), (1.0, 0, 0.01),
            (1.0, 5, float("nan")), (1.0, 5, float("inf")),
            (1.0, 5, -0.01),
        )
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(ValueError):
                thermal_sampler.validate_sampling_arguments(*values)

    def test_tctl_discovery_is_python38_compatible(self):
        globs = {
            "/sys/class/hwmon/hwmon*/name": ["/sys/class/hwmon/hwmon7/name"],
            "/sys/class/hwmon/hwmon7/temp*_label": [
                "/sys/class/hwmon/hwmon7/temp2_label"],
        }
        values = {
            "/sys/class/hwmon/hwmon7/name": "k10temp",
            "/sys/class/hwmon/hwmon7/temp2_label": "Tctl",
        }
        with mock.patch.object(
                thermal_sampler.glob, "glob",
                side_effect=lambda pattern: globs.get(pattern, [])), \
                mock.patch.object(
                    thermal_sampler, "read_text",
                    side_effect=lambda path: values.get(path)):
            self.assertEqual(
                thermal_sampler.find_tctl_path(),
                "/sys/class/hwmon/hwmon7/temp2_input")

    def test_cpu_stat_does_not_double_count_guest_time(self):
        raw = "cpu 10 20 30 40 5 6 7 8 100 50\ncpu0 0 0 0 0\n"
        with mock.patch.object(thermal_sampler, "read_text", return_value=raw):
            total, idle = thermal_sampler.read_cpu_stat()
        self.assertEqual(total, sum((10, 20, 30, 40, 5, 6, 7, 8)))
        self.assertEqual(idle, 45)
        for malformed in ("intr 1 2 3\n", "cpu 1 2 3 4 5 6 7\n",
                          "cpu 1 2 3 4 5 6 7 bad\n",
                          "cpu 1 2 3 4 5 6 7 -1\n"):
            with self.subTest(malformed=malformed), mock.patch.object(
                    thermal_sampler, "read_text", return_value=malformed), \
                    self.assertRaises(RuntimeError):
                thermal_sampler.read_cpu_stat()

    def test_bootstrap_distinguishes_header_append_from_invalid_complete_row(self):
        raw = thermal_bytes()
        header, first_row = raw.splitlines(keepends=True)[:2]
        self.assertEqual(timing._bootstrap_thermal_csv_state(header), "header")
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(header + first_row[:-3]),
            "header")
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(header + first_row), "row")
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(b"utc,partial"), "incomplete")
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(b"utc,wrong\n"),
            "invalid-row")
        invalid = bytearray(header + first_row)
        invalid[0:3] = b"bad"
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(bytes(invalid)), "invalid-row")
        rows = raw.decode("ascii").splitlines()
        fields = rows[1].split(",")
        fields[timing.THERMAL_FIELDS.index(timing.DIMM_FIELDS[0])] = ""
        complete_bad_row = (rows[0] + "\n" + ",".join(fields) + "\n").encode(
            "ascii")
        self.assertEqual(
            timing._bootstrap_thermal_csv_state(complete_bad_row), "row")
        with self.assertRaises(timing.TimingError):
            timing._parse_thermal_csv(complete_bad_row)

    def test_graceful_sampler_seals_csv_before_controller_receives_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            with path.open("xb") as output:
                output.write(b"sample\n")
                thermal_sampler.seal_csv_output(output)
                self.assertEqual(os.fstat(output.fileno()).st_mode & 0o777, 0o444)
            self.assertEqual(path.stat().st_mode & 0o777, 0o444)

    def test_edac_sum_requires_stable_complete_nonnegative_inventory(self):
        paths = (
            "/sys/devices/system/edac/mc/mc0/ce_count",
            "/sys/devices/system/edac/mc/mc1/ce_count",
        )
        with mock.patch.object(
                thermal_sampler.glob, "glob", return_value=list(reversed(paths))), \
                mock.patch.object(
                    thermal_sampler, "read_text", side_effect=["2", "3"]):
            self.assertEqual(thermal_sampler.sum_edac("ce_count", paths), 5)
        for values in (("2", None), ("2", "bad"), ("2", "-1"),
                       (str((1 << 64) - 1), "1")):
            with self.subTest(values=values), mock.patch.object(
                    thermal_sampler.glob, "glob", return_value=list(paths)), \
                    mock.patch.object(
                        thermal_sampler, "read_text", side_effect=values), \
                    self.assertRaises(RuntimeError):
                thermal_sampler.sum_edac("ce_count", paths)

    def test_edac_inventory_empty_or_changed_is_rejected(self):
        path = "/sys/devices/system/edac/mc/mc0/ce_count"
        with mock.patch.object(thermal_sampler.glob, "glob", return_value=[]), \
                self.assertRaisesRegex(RuntimeError, "inventory is empty"):
            thermal_sampler.discover_edac_paths("ce_count")
        with mock.patch.object(
                thermal_sampler.glob, "glob", return_value=[path]), \
                self.assertRaisesRegex(RuntimeError, "inventory changed"):
            thermal_sampler.sum_edac("ce_count", (path, path + ".other"))

    def test_implausible_spd_read_is_retried(self):
        dimm = (1, 0x50)
        with mock.patch.object(thermal_sampler, "DIMMS", [dimm]), \
                mock.patch.object(
                    thermal_sampler, "read_spd5118_temperature",
                    side_effect=[240.75, 48.25]) as read:
            temperatures, pending = thermal_sampler.read_dimm_temperatures(
                {1: object()}, attempts=2, retry_delay=0.0)
        self.assertEqual(temperatures, {dimm: 48.25})
        self.assertEqual(pending, [])
        self.assertEqual(read.call_count, 2)

    def test_repeated_implausible_spd_reads_are_reported_failed(self):
        dimm = (1, 0x50)
        with mock.patch.object(thermal_sampler, "DIMMS", [dimm]), \
                mock.patch.object(
                    thermal_sampler, "read_spd5118_temperature",
                    side_effect=[-40.0, 130.0]):
            temperatures, pending = thermal_sampler.read_dimm_temperatures(
                {1: object()}, attempts=2, retry_delay=0.0)
        self.assertEqual(temperatures, {})
        self.assertEqual(pending, [dimm])

    def test_plausible_spd_boundary_values_and_oserror_retry(self):
        dimms = [(1, 0x50), (1, 0x51)]
        with mock.patch.object(thermal_sampler, "DIMMS", dimms), \
                mock.patch.object(
                    thermal_sampler, "read_spd5118_temperature",
                    side_effect=[-39.75, OSError("transient"), 129.75]):
            temperatures, pending = thermal_sampler.read_dimm_temperatures(
                {1: object()}, attempts=2, retry_delay=0.0)
        self.assertEqual(temperatures, {
            (1, 0x50): -39.75,
            (1, 0x51): 129.75,
        })
        self.assertEqual(pending, [])


class PostTaskBindingTests(unittest.TestCase):
    @staticmethod
    def owner(root: Path):
        process = mock.Mock(pid=123)
        return timing.SamplerOwner(
            process, 456, 789, {"pid": 456}, root / "thermal.csv.part",
            root / "bootstrap.pid", 0)

    @staticmethod
    def design():
        return {
            "thermal_core": 127,
            "tools": {
                "fuser": {"path": "/frozen/fuser"},
                "sudo": {"path": "/frozen/sudo"},
                "timeout": {"path": "/frozen/timeout"},
            },
        }

    def test_task_is_not_bound_when_post_task_environment_fails(self):
        for failure in ("filler", "identity"):
            with self.subTest(failure=failure), \
                    tempfile.TemporaryDirectory() as temporary:
                completed = []
                filler = (999,) if failure == "filler" else ()
                identity_ok = failure != "identity"
                with mock.patch.object(
                        timing, "filler_pids", return_value=filler), \
                        mock.patch.object(
                            timing, "process_identity_matches",
                            return_value=identity_ok), \
                        mock.patch.object(
                            timing, "sole_i2c_readers", return_value=(456,)):
                    with self.assertRaises(timing.TimingError):
                        timing._bind_task_after_environment_checks(
                            completed, {"job": 0}, mock.Mock(),
                            self.owner(Path(temporary)), self.design())
                self.assertEqual(completed, [])


class ThermalArchiveTests(unittest.TestCase):
    @staticmethod
    def make_root(root: Path):
        (root / "segments").mkdir()
        (root / "interrupted").mkdir()
        return {"tools": {
            "fuser": {"path": "/frozen/fuser"},
            "sudo": {"path": "/frozen/sudo"},
            "timeout": {"path": "/frozen/timeout"},
        }}

    def test_invalid_thermal_part_is_archived_and_verifiable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            csv_part = root / "segments/.segment0000.thermal.csv.part"
            pid_file = root / "segments/.segment0000.bootstrap.pid"
            csv_part.write_bytes(thermal_bytes())
            pid_file.write_text("123\n", encoding="ascii")
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                record = timing.archive_invalid_thermal_segment(
                    root, design, 0, "unit-invalid-thermal")
            self.assertIsNotNone(record)
            self.assertFalse(csv_part.exists())
            self.assertFalse(pid_file.exists())
            archive = root / "interrupted" / record["archive"]
            timing._verify_thermal_failure_archive(archive)
            timing._verify_interrupted_archive(root, [record])
            adopted = timing.adopt_unbound_interrupted_archives(root)
            timing._verify_interrupted_archive(root, adopted)
            (archive / "thermal.csv.part").write_bytes(b"tampered\n")
            with self.assertRaisesRegex(
                    timing.TimingError, "archive changed|does not replay"):
                timing._verify_interrupted_archive(root, [record])

    def test_incomplete_archive_is_completed_on_next_invocation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            (root / "segments/.segment0000.thermal.csv.part").write_bytes(
                thermal_bytes())
            staging = root / "interrupted/.thermal-segment0000.invalid.123.part"
            staging.mkdir()
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                recovered = timing.recover_orphaned_thermal_segments(
                    root, design)
                self.assertEqual(
                    timing.recover_orphaned_thermal_segments(root, design), [])
            self.assertEqual(len(recovered), 1)
            archive = root / "interrupted" / recovered[0]["archive"]
            timing._verify_thermal_failure_archive(archive)
            self.assertFalse(any((root / "segments").iterdir()))

    def test_atomic_intent_and_failure_receipts_recover_before_and_after_link(self):
        for receipt_name in ("intent", "thermal_failure"):
            for after_link in (False, True):
                with self.subTest(receipt=receipt_name, after_link=after_link), \
                        tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    design = self.make_root(root)
                    staging = root / \
                        "interrupted/.thermal-segment0000.invalid.123.part"
                    staging.mkdir()
                    intent = timing.sealed_record(
                        "wirehair.wh2.p32_dispatch.thermal_archive_intent.v2", {
                            "created_utc": "2026-07-20T00:00:00.000Z",
                            "segment": 0, "reason": "atomic-unit-test",
                        })
                    intent_raw = timing.canonical_json(intent)
                    if receipt_name == "intent":
                        final_path = staging / "intent.json"
                        part_path = staging / "intent.json.part"
                        part_path.write_bytes(intent_raw)
                        if after_link:
                            os.link(part_path, final_path)
                    else:
                        (staging / "intent.json").write_bytes(intent_raw)
                        artifacts = timing._thermal_archive_file_receipts(staging)
                        failure = timing.sealed_record(
                            "wirehair.wh2.p32_dispatch.thermal_failure.v2", {
                                "archived_utc": "2026-07-20T00:00:01.000Z",
                                "segment": 0, "reason": "atomic-unit-test",
                                "intent_sha256": timing.sha256_file(
                                    staging / "intent.json"),
                                "artifacts": artifacts,
                            })
                        final_path = staging / "thermal_failure.json"
                        part_path = staging / "thermal_failure.json.part"
                        part_path.write_bytes(timing.canonical_json(failure))
                        if after_link:
                            os.link(part_path, final_path)
                    with mock.patch.object(
                            timing, "sole_i2c_readers", return_value=()):
                        record = timing._complete_thermal_archive(
                            root, design, staging)
                    archive = root / "interrupted" / record["archive"]
                    self.assertFalse((archive / "intent.json.part").exists())
                    self.assertFalse(
                        (archive / "thermal_failure.json.part").exists())
                    timing._verify_thermal_failure_archive(archive)

    def test_invalid_atomic_receipt_partials_are_preserved_and_bound(self):
        for receipt_name in ("intent", "thermal_failure"):
            with self.subTest(receipt=receipt_name), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                design = self.make_root(root)
                staging = root / \
                    "interrupted/.thermal-segment0000.invalid.123.part"
                staging.mkdir()
                if receipt_name == "thermal_failure":
                    intent = timing.sealed_record(
                        "wirehair.wh2.p32_dispatch.thermal_archive_intent.v2", {
                            "created_utc": "2026-07-20T00:00:00.000Z",
                            "segment": 0, "reason": "partial-unit-test",
                        })
                    (staging / "intent.json").write_bytes(
                        timing.canonical_json(intent))
                (staging / (receipt_name + ".json.part")).write_bytes(
                    b"truncated evidence")
                with mock.patch.object(
                        timing, "sole_i2c_readers", return_value=()):
                    record = timing._complete_thermal_archive(
                        root, design, staging)
                archive = root / "interrupted" / record["archive"]
                remnants = [path.name for path in archive.iterdir()
                             if ".interrupted." in path.name]
                self.assertEqual(len(remnants), 1)
                self.assertTrue(remnants[0].startswith(receipt_name + "."))
                timing._verify_thermal_failure_archive(archive)

    def test_atomic_remnant_filename_must_bind_its_content_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            staging = root / \
                "interrupted/.thermal-segment0000.invalid.123.part"
            staging.mkdir()
            intent = timing.sealed_record(
                "wirehair.wh2.p32_dispatch.thermal_archive_intent.v2", {
                    "created_utc": "2026-07-20T00:00:00.000Z",
                    "segment": 0, "reason": "remnant-name-test",
                })
            (staging / "intent.json").write_bytes(
                timing.canonical_json(intent))
            (staging / ("intent.interrupted.%s.part" % ("0" * 64))).write_bytes(
                b"content with a different hash")
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                with self.assertRaisesRegex(
                        timing.TimingError, "remnant name hash changed"):
                    timing._complete_thermal_archive(root, design, staging)

    def test_unbound_final_thermal_and_terminal_partial_are_archived(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            (root / "segments/segment0000.thermal.csv").write_bytes(
                thermal_bytes())
            (root / "segments/segment0000.terminal.json.part").write_bytes(
                b"partial terminal evidence")
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                recovered = timing.recover_orphaned_thermal_segments(
                    root, design)
            self.assertEqual(len(recovered), 1)
            archive = root / "interrupted" / recovered[0]["archive"]
            self.assertEqual(
                {"thermal.csv.unbound-final", "terminal.json.part"},
                {"thermal.csv.unbound-final", "terminal.json.part"} &
                {path.name for path in archive.iterdir()})
            self.assertFalse(any((root / "segments").iterdir()))
            timing._verify_thermal_failure_archive(archive)

    def test_terminal_atomic_link_remnant_is_removed_only_when_identical(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            terminal = root / "segments/segment0000.terminal.json"
            raw = timing.canonical_json(timing.sealed_record(
                "wirehair.wh2.p32_dispatch.segment_terminal.v2",
                {"segment": 0}))
            terminal.write_bytes(raw)
            os.link(terminal, terminal.with_name(terminal.name + ".part"))
            (root / "segments/segment0000.thermal.csv").write_bytes(
                thermal_bytes())
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                self.assertEqual(
                    timing.recover_orphaned_thermal_segments(root, design), [])
            self.assertFalse(terminal.with_name(terminal.name + ".part").exists())

    def test_terminal_without_final_thermal_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            terminal = timing.sealed_record(
                "wirehair.wh2.p32_dispatch.segment_terminal.v2",
                {"segment": 0})
            (root / "segments/segment0000.terminal.json").write_bytes(
                timing.canonical_json(terminal))
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                with self.assertRaisesRegex(
                        timing.TimingError, "lacks its final thermal CSV"):
                    timing.recover_orphaned_thermal_segments(root, design)

    def test_outer_thermal_archive_receipt_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            (root / "segments/.segment0000.thermal.csv.part").write_bytes(
                thermal_bytes())
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                record = timing.archive_invalid_thermal_segment(
                    root, design, 0, "outer-exactness")
            mutations = []
            extra = json.loads(json.dumps(record))
            extra["extra"] = True
            mutations.append(extra)
            source = json.loads(json.dumps(record))
            source["source"] = "prior-controller-interruption"
            mutations.append(source)
            reason = json.loads(json.dumps(record))
            reason["reason"] = "changed"
            mutations.append(reason)
            file_extra = json.loads(json.dumps(record))
            first = next(iter(file_extra["files"].values()))
            first["extra"] = True
            mutations.append(file_extra)
            for mutation in mutations:
                with self.subTest(mutation=mutation), self.assertRaises(
                        timing.TimingError):
                    timing._verify_interrupted_archive(root, [mutation])
            for malformed in (None, 7, "record", []):
                with self.subTest(malformed=malformed), self.assertRaises(
                        timing.TimingError):
                    timing._verify_interrupted_archive(root, [malformed])

    def test_adoption_semantically_rejects_tampered_thermal_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            (root / "segments/.segment0000.thermal.csv.part").write_bytes(
                thermal_bytes())
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                record = timing.archive_invalid_thermal_segment(
                    root, design, 0, "adoption-tamper")
            archive = root / "interrupted" / record["archive"]
            (archive / "thermal.csv.part").write_bytes(b"tampered")
            with self.assertRaises(timing.TimingError):
                timing.adopt_unbound_interrupted_archives(root)

    def test_live_i2c_reader_prevents_thermal_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design = self.make_root(root)
            csv_part = root / "segments/.segment0000.thermal.csv.part"
            csv_part.write_bytes(thermal_bytes())
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=(999,)):
                with self.assertRaisesRegex(
                        timing.TimingError, "I2C reader lives"):
                    timing.archive_invalid_thermal_segment(
                        root, design, 0, "must-not-move-live-evidence")
            self.assertTrue(csv_part.is_file())


class ProcessIdentityTests(unittest.TestCase):
    def setUp(self):
        patcher = mock.patch.object(
            timing, "_root_readonly_single_link_stat",
            side_effect=lambda path, _description: root_owned_sealed_stat(path))
        patcher.start()
        self.addCleanup(patcher.stop)

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
            "Name:\tpython3\nUid:\t0\t0\t0\t0\n"
            "Cpus_allowed_list:\t127\n", encoding="ascii")
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

    def test_nonroot_sampler_identity_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            proc, boot, csv_path = self.make_fixture(root)
            status = proc / "123/status"
            status.write_text(status.read_text(encoding="ascii").replace(
                "Uid:\t0\t0\t0\t0", "Uid:\t1000\t1000\t1000\t1000"),
                encoding="ascii")
            with self.assertRaisesRegex(timing.TimingError, "exact root process"):
                timing.capture_process_identity(
                    123, 127, csv_path, proc_root=proc, boot_id_path=boot)

    def test_fuser_return_code_and_pid_output_are_exact(self):
        device = Path("/dev/i2c-1")
        valid = ((1, b"", b"", ()),
                 (0, b" 42", b"/dev/i2c-1:         \n", (42,)),
                 (0, b" 42 7 42", b"/dev/i2c-1: \n", (7, 42)),
                 (0, b"  3212   7", b"/dev/i2c-1:\n", (7, 3212)))
        for returncode, stdout, stderr, expected in valid:
            with self.subTest(returncode=returncode, stdout=stdout):
                self.assertEqual(timing._parse_fuser_device_result(
                    device, returncode, stdout, stderr), expected)
        invalid = (
            (0, b"", b"/dev/i2c-1: \n"),
            (0, b"garbage", b"/dev/i2c-1: \n"),
            (0, b" 0", b"/dev/i2c-1: \n"),
            (0, b"42", b"/dev/i2c-1: \n"),
            (0, b" 42\n", b"/dev/i2c-1: \n"),
            (0, b" 42", b"fatal diagnostic\n"),
            (1, b" 42", b""), (1, b"", b"missing device\n"),
            (2, b"", b""),
        )
        for returncode, stdout, stderr in invalid:
            with self.subTest(returncode=returncode, stdout=stdout,
                              stderr=stderr), self.assertRaises(
                                  timing.TimingError):
                timing._parse_fuser_device_result(
                    device, returncode, stdout, stderr)

    def test_fuser_queries_each_required_i2c_device_independently(self):
        outcomes = ((0, b" 42", b"/dev/i2c-1: \n"),
                    (0, b" 7 42", b"/dev/i2c-2: \n"))
        with mock.patch.object(timing.Path, "is_symlink", return_value=False), \
                mock.patch.object(
                    timing.Path, "is_char_device", return_value=True), \
                mock.patch.object(
                    timing, "run_privileged_bounded",
                    side_effect=outcomes) as bounded:
            self.assertEqual(timing.sole_i2c_readers(
                Path("/fuser"), Path("/sudo"), Path("/timeout")), (7, 42))
        commands = [call.args[2] for call in bounded.call_args_list]
        self.assertEqual(commands, [
            ("/fuser", "/dev/i2c-1"), ("/fuser", "/dev/i2c-2")])

    def test_public_cmdline_redacts_ephemeral_pidfile(self):
        identity = {
            "pid": 123, "cmdline": ["python3", "sampler.py", "--pid-file",
                                      "/tmp/bootstrap.pid"],
        }
        public = timing._public_sampler_identity(identity)
        self.assertEqual(public["cmdline"][-1], "<ephemeral-bootstrap-only>")
        self.assertNotIn("bootstrap.pid", json.dumps(public))


class ThermalSealingTests(unittest.TestCase):
    @staticmethod
    def with_temperature(raw: bytes, field: str, value: str) -> bytes:
        rows = raw.decode("ascii").splitlines()
        header = rows[0].split(",")
        for index in range(1, len(rows)):
            fields = rows[index].split(",")
            fields[header.index(field)] = value
            rows[index] = ",".join(fields)
        return ("\n".join(rows) + "\n").encode("ascii")

    def test_post_end_coverage_and_health(self):
        summary = timing.validate_thermal_interval(thermal_bytes(), 100.0, 101.0)
        self.assertEqual(summary["sample_count"], 3)
        self.assertGreaterEqual(summary["post_end_margin_s"], 0.0)

    def test_controller_renames_root_sealed_csv_without_chmod(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "segments").mkdir()
            part = root / "segments/.segment0000.thermal.csv.part"
            part.write_bytes(thermal_bytes())
            part.chmod(0o444)
            real_stat = part.stat()
            sealed_stat = mock.Mock(
                st_mode=(real_stat.st_mode & ~0o777) | 0o444,
                st_uid=0, st_nlink=1, st_dev=real_stat.st_dev,
                st_ino=real_stat.st_ino)
            identity = {"pid": 456, "csv_device": real_stat.st_dev,
                        "csv_inode": real_stat.st_ino}
            owner = timing.SamplerOwner(
                mock.Mock(), 456, 789, identity, part,
                root / "bootstrap.pid", 0)
            stop = timing.SamplerStop(
                identity_end=identity, graceful=True,
                mechanism="sudo-timeout-python-pidfd-identity-verified-sigterm",
                forced_reason=None)
            with mock.patch.object(
                    timing, "process_identity_matches", return_value=True), \
                    mock.patch.object(
                        timing, "latest_thermal_time", return_value=101.5), \
                    mock.patch.object(
                        timing, "_stop_owned_sampler", return_value=stop), \
                    mock.patch.object(
                        timing, "_root_sealed_thermal_stat",
                        return_value=sealed_stat), \
                    mock.patch.object(
                        timing.os, "chmod",
                        side_effect=PermissionError("not owner")) as chmod:
                result = timing.stop_sampler(
                    owner, root, {"thermal_core": 127}, 100.0, 101.0)
            chmod.assert_not_called()
            self.assertFalse(part.exists())
            self.assertTrue((root / "segments/segment0000.thermal.csv").is_file())
            self.assertEqual(result["thermal_csv_inode"], real_stat.st_ino)

    def test_root_sealed_csv_requires_owner_mode_and_single_link(self):
        path = mock.Mock()
        for uid, mode, links in ((1000, 0o100444, 1), (0, 0o100644, 1),
                                 (0, 0o100444, 2), (0, 0o040444, 1)):
            path.lstat.return_value = mock.Mock(
                st_uid=uid, st_mode=mode, st_nlink=links)
            with self.subTest(uid=uid, mode=oct(mode), links=links), \
                    self.assertRaisesRegex(
                        timing.TimingError, "root sampler did not seal"):
                timing._root_sealed_thermal_stat(path)

    def test_missing_post_end_sample_is_rejected(self):
        with self.assertRaisesRegex(timing.TimingError, "bracket benchmark end"):
            timing.validate_thermal_interval(
                thermal_bytes((99.5, 100.5)), 100.0, 101.0)

    def test_nonpositive_cpu_temperature_is_rejected(self):
        for value in ("0.0", "-0.25"):
            with self.subTest(value=value):
                rows = thermal_bytes().decode("ascii").splitlines()
                header = rows[0].split(",")
                fields = rows[1].split(",")
                fields[header.index("cpu_tctl_c")] = value
                fields[header.index(timing.DIMM_FIELDS[0])] = "-39.75"
                rows[1] = ",".join(fields)
                raw = ("\n".join(rows) + "\n").encode("ascii")
                with self.assertRaisesRegex(
                        timing.TimingError, "thermal temperature implausible"):
                    timing.validate_thermal_interval(raw, 100.0, 101.0)

    def test_negative_dimm_temperature_remains_in_its_own_domain(self):
        rows = thermal_bytes().decode("ascii").splitlines()
        header = rows[0].split(",")
        for index in range(1, len(rows)):
            fields = rows[index].split(",")
            fields[header.index(timing.DIMM_FIELDS[0])] = "-39.75"
            rows[index] = ",".join(fields)
        raw = ("\n".join(rows) + "\n").encode("ascii")
        summary = timing.validate_thermal_interval(raw, 100.0, 101.0)
        self.assertEqual(summary["dimm_max_c"][timing.DIMM_FIELDS[0]], -39.75)

    def test_live_dimm_limit_passes_below_and_aborts_at_equality(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            now = time.monotonic()
            recent = (now - 1.5, now - 0.75, now - 0.1)
            path.write_bytes(self.with_temperature(
                thermal_bytes(recent), timing.DIMM_FIELDS[0], "83.75"))
            live = timing.enforce_live_thermal_safety(path, now_s=now)
            self.assertEqual(live["dimm_max_c"], 83.75)
            for value in ("84.0", "84.25"):
                with self.subTest(value=value):
                    path.write_bytes(self.with_temperature(
                        thermal_bytes(recent), timing.DIMM_FIELDS[0], value))
                    with self.assertRaisesRegex(
                            timing.TimingError, "abort threshold"):
                        timing.enforce_live_thermal_safety(path, now_s=now)
                    with self.assertRaisesRegex(
                            timing.TimingError, "temperature gate exceeded"):
                        timing.validate_thermal_interval(
                            path.read_bytes(), recent[0] + 0.1, recent[-1] - 0.05)

    def test_live_thermal_safety_rejects_stale_and_gapped_streams(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            now = time.monotonic()
            path.write_bytes(thermal_bytes((now - 6.0, now - 5.0, now - 4.0)))
            with self.assertRaisesRegex(timing.TimingError, "sample is stale"):
                timing.enforce_live_thermal_safety(path, now_s=now)
            path.write_bytes(thermal_bytes((now - 4.0, now - 0.2)))
            with self.assertRaisesRegex(timing.TimingError, "cadence gap"):
                timing.enforce_live_thermal_safety(path, now_s=now)

    def test_live_recency_clock_is_sampled_after_stable_csv_read(self):
        rows = timing._parse_thermal_csv(thermal_bytes((100.0, 101.0, 102.0)))
        events = []

        def read_rows(_path):
            events.append("read")
            return rows

        def read_clock():
            events.append("clock")
            return 102.25

        with mock.patch.object(
                timing, "stable_thermal_rows", side_effect=read_rows), \
                mock.patch.object(
                    timing.time, "monotonic", side_effect=read_clock):
            live = timing.enforce_live_thermal_safety(Path("unused"))
        self.assertEqual(events, ["read", "clock"])
        self.assertEqual(live["latest_age_s"], 0.25)

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
                "thermal_csv_uid": 0, "thermal_csv_mode": 0o444,
                "thermal_csv_nlink": 1,
                "timing_start_monotonic_s": 100.0,
                "benchmark_end_monotonic_s": 101.0,
                "thermal_summary": timing.validate_thermal_interval(
                    original, 100.0, 101.0),
            }
            with mock.patch.object(
                    timing, "_root_sealed_thermal_stat",
                    side_effect=root_owned_sealed_stat):
                timing.verify_terminal_thermal(path, receipt)
                with self.assertRaisesRegex(
                        timing.TimingError, "terminal thermal receipt mismatch"):
                    timing.verify_terminal_thermal(
                        path, {**receipt, "thermal_csv_uid": False})
                with path.open("ab") as stream:
                    stream.write(
                        thermal_bytes((102.5,)).splitlines(keepends=True)[1])
                with self.assertRaisesRegex(
                        timing.TimingError, "terminal thermal receipt mismatch"):
                    timing.verify_terminal_thermal(path, receipt)

    def test_offline_thermal_verifier_rejects_mode_and_hardlink_mutation(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            raw = thermal_bytes()
            path.write_bytes(raw)
            path.chmod(0o444)
            actual = path.stat()
            receipt = {
                "thermal_csv_sha256": timing.sha256_bytes(raw),
                "thermal_csv_size": len(raw),
                "thermal_csv_device": actual.st_dev,
                "thermal_csv_inode": actual.st_ino,
                "thermal_csv_uid": 0, "thermal_csv_mode": 0o444,
                "thermal_csv_nlink": 1,
                "timing_start_monotonic_s": 100.0,
                "benchmark_end_monotonic_s": 101.0,
                "thermal_summary": timing.validate_thermal_interval(
                    raw, 100.0, 101.0),
            }

            def strict_root_stat(candidate):
                current = candidate.lstat()
                if (current.st_mode & 0o777) != 0o444 or current.st_nlink != 1:
                    raise timing.TimingError("root sampler did not seal")
                return root_owned_sealed_stat(candidate)

            with mock.patch.object(
                    timing, "_root_sealed_thermal_stat",
                    side_effect=strict_root_stat):
                timing.verify_terminal_thermal(path, receipt)
                path.chmod(0o644)
                with self.assertRaisesRegex(timing.TimingError, "did not seal"):
                    timing.verify_terminal_thermal(path, receipt)
                path.chmod(0o444)
                os.link(path, Path(temporary) / "extra-link.csv")
                with self.assertRaisesRegex(timing.TimingError, "did not seal"):
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

    def test_bad_post_end_evidence_still_stops_owned_sampler(self):
        for failure in ("malformed", "timeout"):
            with self.subTest(failure=failure), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                csv_path = root / "thermal.csv.part"
                csv_path.write_bytes(thermal_bytes())
                pid_file = root / "bootstrap.pid"
                process = mock.Mock(returncode=0)
                process.pid = 123
                process.communicate.return_value = (b"", b"")
                identity = {"pid": 123, "start_tick": 4242,
                            "process_group": 123, "session_id": 123,
                            "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
                            "cmdline_sha256": "a" * 64}
                owner = timing.SamplerOwner(
                    process, 123, 4242, identity, csv_path, pid_file, 0)
                design = {
                    "thermal_core": 127,
                    "tools": {
                        "sudo": {"path": "/frozen/sudo"},
                        "fuser": {"path": "/frozen/fuser"},
                        "timeout": {"path": "/frozen/timeout"},
                        "env": {"path": "/frozen/env"},
                        "python3": {"path": "/frozen/python3"},
                    },
                }
                latest = (mock.Mock(side_effect=timing.TimingError("bad row"))
                          if failure == "malformed" else
                          mock.Mock(return_value=100.5))
                monotonic = ([100.0, 116.0, 200.0]
                             if failure == "timeout" else [100.0, 200.0])
                with mock.patch.object(
                        timing, "process_identity_matches",
                        side_effect=[True, False]), \
                        mock.patch.object(
                            timing, "capture_process_identity",
                            return_value=identity), \
                        mock.patch.object(
                            timing, "sole_i2c_readers",
                            side_effect=[(123,), ()]), \
                        mock.patch.object(
                            timing, "latest_thermal_time", latest), \
                        mock.patch.object(
                            timing, "run_privileged_bounded",
                            return_value=(0, b"", b"")) as stop, \
                        mock.patch.object(
                            timing.time, "monotonic", side_effect=monotonic):
                    with self.assertRaisesRegex(
                            timing.TimingError,
                            "coverage check failed|did not emit"):
                        timing.stop_sampler(
                            owner, root, design, 100.0, 101.0)
                stop.assert_called_once()
                process.communicate.assert_called_once_with(timeout=15)

    def test_bootstrap_timeout_prefers_identity_then_falls_back_to_session(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            csv_path = root / "thermal.csv.part"
            csv_path.write_bytes(thermal_bytes())
            pid_file = root / "bootstrap.pid"
            pid_file.write_text("123\n", encoding="ascii")
            process = mock.Mock(pid=999)
            argv = ["/frozen/python3", "/frozen/sampler.py", "--csv",
                    str(csv_path)]
            boot_id = "01234567-89ab-cdef-0123-456789abcdef"
            identity = {"pid": 123, "cmdline": argv,
                        "process_group": 999, "session_id": 999,
                        "boot_id": boot_id}
            design = {
                "thermal_core": 127,
                "tools": {
                    "sudo": {"path": "/frozen/sudo"},
                    "fuser": {"path": "/frozen/fuser"},
                    "timeout": {"path": "/frozen/timeout"},
                },
            }
            with mock.patch.object(
                    timing, "capture_process_identity",
                    return_value=identity), \
                    mock.patch.object(
                        timing, "sole_i2c_readers", return_value=(123,)), \
                    mock.patch.object(timing, "_stop_owned_sampler") as stop, \
                    mock.patch.object(
                        timing, "_kill_owned_process_session") as kill:
                timing._cleanup_timed_out_sampler(
                    root, design, process, csv_path, pid_file, 7, argv,
                    launcher_start_tick=4242, boot_id=boot_id)
                stop.assert_called_once()
                owner = stop.call_args.args[0]
                self.assertEqual((owner.pid, owner.segment), (123, 7))
                self.assertEqual(owner.launcher_start_tick, 4242)
                kill.assert_not_called()

                bad_identity = dict(identity, cmdline=["unexpected"])
                with mock.patch.object(
                        timing, "capture_process_identity",
                        return_value=bad_identity), \
                        mock.patch.object(
                            timing, "sole_i2c_readers",
                            side_effect=[(123,), ()]):
                    timing._cleanup_timed_out_sampler(
                        root, design, process, csv_path, pid_file, 8, argv,
                        launcher_start_tick=4242, boot_id=boot_id)
                self.assertEqual(stop.call_count, 1)
                kill.assert_called_once()

                pid_file.unlink()
                kill.reset_mock()
                with mock.patch.object(
                        timing, "sole_i2c_readers", return_value=()):
                    timing._cleanup_timed_out_sampler(
                        root, design, process, csv_path, pid_file, 9, argv,
                        launcher_start_tick=4242, boot_id=boot_id)
                kill.assert_called_once()

    def test_hung_sampler_uses_forced_session_cleanup_and_is_not_promotional(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            process = mock.Mock(pid=123, returncode=-9)
            identity = {
                "pid": 456, "start_tick": 789, "process_group": 123,
                "session_id": 123,
                "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
                "cmdline_sha256": "a" * 64,
            }
            owner = timing.SamplerOwner(
                process, 456, 321, identity, root / "thermal.csv.part",
                root / "bootstrap.pid", 0)
            design = {
                "thermal_core": 127,
                "tools": {
                    "sudo": {"path": "/frozen/sudo"},
                    "fuser": {"path": "/frozen/fuser"},
                    "timeout": {"path": "/frozen/timeout"},
                    "env": {"path": "/frozen/env"},
                    "python3": {"path": "/frozen/python3"},
                },
            }
            with mock.patch.object(
                    timing, "capture_process_identity", return_value=identity), \
                    mock.patch.object(
                        timing, "sole_i2c_readers", side_effect=[(456,), ()]), \
                    mock.patch.object(
                        timing, "run_privileged_bounded",
                        return_value=(0, b"", b"")), \
                    mock.patch.object(
                        timing, "process_identity_matches", return_value=True), \
                    mock.patch.object(
                        timing, "_kill_owned_sampler_session",
                        return_value=(b"", b"")) as kill, \
                    mock.patch.object(
                        timing.time, "monotonic", side_effect=[0.0, 16.0]):
                stop = timing._stop_owned_sampler(owner, root, design)
            self.assertFalse(stop.graceful)
            self.assertIn("ignored graceful stop", str(stop.forced_reason))
            kill.assert_called_once_with(owner, root, design)

    def test_fuser_failure_after_identity_proof_still_cleans_session(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            process = mock.Mock(pid=123)
            identity = {
                "pid": 456, "process_group": 123, "session_id": 123,
                "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
            }
            owner = timing.SamplerOwner(
                process, 456, 321, identity, root / "thermal.csv.part",
                root / "bootstrap.pid", 0)
            design = {
                "thermal_core": 127,
                "tools": {
                    "sudo": {"path": "/frozen/sudo"},
                    "fuser": {"path": "/frozen/fuser"},
                    "timeout": {"path": "/frozen/timeout"},
                },
            }
            with mock.patch.object(
                    timing, "capture_process_identity", return_value=identity), \
                    mock.patch.object(
                        timing, "sole_i2c_readers",
                        side_effect=timing.TimingError("fuser timeout")), \
                    mock.patch.object(
                        timing, "_kill_owned_sampler_session",
                        return_value=(b"", b"")) as kill:
                with self.assertRaisesRegex(
                        timing.TimingError, "inspection failed.*fuser timeout"):
                    timing._stop_owned_sampler(owner, root, design)
            kill.assert_called_once_with(owner, root, design)

    def test_coverage_and_shutdown_failures_are_both_reported(self):
        owner = timing.SamplerOwner(
            mock.Mock(), 123, 4242, {}, Path("thermal.csv.part"),
            Path("bootstrap.pid"), 0)
        with mock.patch.object(
                timing, "process_identity_matches", return_value=True), \
                mock.patch.object(
                    timing, "latest_thermal_time",
                    side_effect=timing.TimingError("bad row")), \
                mock.patch.object(
                    timing, "_stop_owned_sampler",
                    side_effect=timing.TimingError("stop failed")):
            with self.assertRaisesRegex(
                    timing.TimingError,
                    "thermal evidence failure.*bad row.*shutdown failure.*stop failed"):
                timing.stop_sampler(
                    owner, Path("/tmp"), {"thermal_core": 127}, 1.0, 2.0)


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
    def setUp(self):
        patcher = mock.patch.object(
            timing, "_root_sealed_thermal_stat",
            side_effect=root_owned_sealed_stat)
        patcher.start()
        self.addCleanup(patcher.stop)

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
            "process_group": 1234, "session_id": 1234,
            "boot_id": "01234567-89ab-cdef-0123-456789abcdef",
            "cmdline": timing._sampler_cmdline(root, design, 0, public=True),
            "cmdline_sha256": timing.sha256_bytes(
                b"\0".join(os.fsencode(value) for value in
                            timing._sampler_cmdline(
                                root, design, 0, public=False)) + b"\0"),
            "uid": 0, "affinity": [127],
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
                    "sudo-timeout-python-pidfd-identity-verified-sigterm",
                "post_end_sample": True,
                "thermal_csv_name": thermal.name,
                "thermal_csv_sha256": timing.sha256_bytes(raw),
                "thermal_csv_size": len(raw),
                "thermal_csv_device": stat.st_dev,
                "thermal_csv_inode": stat.st_ino,
                "thermal_csv_uid": 0, "thermal_csv_mode": 0o444,
                "thermal_csv_nlink": 1,
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

    def test_failed_terminal_can_leave_postcheck_task_unbound_for_rerun(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, terminal, path = self.make_fixture(root)
            terminal["status"] = "failed"
            terminal["completed_tasks"] = []
            terminal["failure"] = {
                "class": "TimingError",
                "message": "post-task environment changed",
            }
            path.write_bytes(timing.canonical_json(self.reseal(terminal)))
            with self.assertRaisesRegex(
                    timing.TimingError,
                    "task ledger is not exactly terminal-bound"):
                timing.verify_segments(root, design, entries)
            records = timing.verify_segments(
                root, design, entries, allowed_unbound_job=0)
            self.assertEqual(len(records), 1)

            design["task_count"] = 2
            entries[1] = {"receipt_sha256": "c" * 64,
                          "receipt": {"attempt": 0}}
            with self.assertRaisesRegex(
                    timing.TimingError,
                    "task ledger is not exactly terminal-bound"):
                timing.verify_segments(
                    root, design, entries, allowed_unbound_job=0)

    def test_unbound_task_allowance_does_not_relax_archive_ledger(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, terminal, path = self.make_fixture(root)
            terminal["status"] = "failed"
            terminal["completed_tasks"] = []
            terminal["failure"] = {
                "class": "TimingError",
                "message": "post-task environment changed",
            }
            path.write_bytes(timing.canonical_json(self.reseal(terminal)))
            (root / "interrupted/unbound-extra").mkdir()
            with self.assertRaisesRegex(
                    timing.TimingError,
                    "interrupted archive ledger is not exact"):
                timing.verify_segments(
                    root, design, entries, allowed_unbound_job=0)

    def test_campaign_path_may_contain_bootstrap_pid_substring(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "valid-bootstrap.pid-campaign"
            root.mkdir()
            design, entries, _terminal, _path = self.make_fixture(root)
            self.assertEqual(len(timing.verify_segments(root, design, entries)), 1)

    def test_terminal_replays_private_cmdline_hash_and_boot_identity(self):
        for field, value in (("cmdline_sha256", "a" * 64),
                             ("boot_id", "A1234567-89ab-cdef-0123-456789abcdef")):
            with self.subTest(field=field), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                design, entries, terminal, path = self.make_fixture(root)
                start = dict(terminal["sampler_identity_start"])
                start[field] = value
                terminal["sampler_identity_start"] = start
                terminal["sampler_identity_end"] = start
                path.write_bytes(timing.canonical_json(self.reseal(terminal)))
                with self.assertRaisesRegex(
                        timing.TimingError,
                        "terminal cmdline did not redact ephemeral PID file"):
                    timing.verify_segments(root, design, entries)

    def test_terminal_identity_integer_fields_reject_booleans(self):
        for field, value in (("start_tick", True), ("uid", False)):
            with self.subTest(field=field), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                design, entries, terminal, path = self.make_fixture(root)
                start = dict(terminal["sampler_identity_start"])
                start[field] = value
                terminal["sampler_identity_start"] = start
                terminal["sampler_identity_end"] = start
                path.write_bytes(timing.canonical_json(self.reseal(terminal)))
                with self.assertRaises(timing.TimingError):
                    timing.verify_segments(root, design, entries)

    def test_terminal_segment_zero_rejects_boolean_alias(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, terminal, path = self.make_fixture(root)
            terminal["segment"] = False
            path.write_bytes(timing.canonical_json(self.reseal(terminal)))
            with self.assertRaisesRegex(
                    timing.TimingError, "identity/status malformed"):
                timing.verify_segments(root, design, entries)

    def test_terminal_symlink_is_rejected_offline(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, _terminal, path = self.make_fixture(root)
            target = root / "terminal-target.json"
            path.rename(target)
            path.symlink_to(target)
            with self.assertRaisesRegex(timing.TimingError, "missing or unsafe"):
                timing.verify_segments(root, design, entries)

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

    def test_terminal_can_bind_only_same_index_thermal_failure_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            design, entries, terminal, path = self.make_fixture(root)
            design["tools"].update({
                "fuser": {"path": "/frozen/fuser"},
                "sudo": {"path": "/frozen/sudo"},
                "timeout": {"path": "/frozen/timeout"},
            })
            (root / "segments/.segment0007.thermal.csv.part").write_bytes(
                thermal_bytes())
            with mock.patch.object(
                    timing, "sole_i2c_readers", return_value=()):
                record = timing.archive_invalid_thermal_segment(
                    root, design, 7, "wrong-terminal-index")
            terminal["recovered_interrupted_transactions"] = [record]
            path.write_bytes(timing.canonical_json(self.reseal(terminal)))
            with self.assertRaisesRegex(
                    timing.TimingError,
                    "thermal interrupted receipt identity changed"):
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
