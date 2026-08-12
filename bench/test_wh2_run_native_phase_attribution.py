#!/usr/bin/env python3
"""Adversarial tests for the sealed native n16 phase runner."""

from __future__ import annotations

import copy
from contextlib import ExitStack
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_run_native_phase_attribution as subject
import test_wh2_phase_attribution as phase_fixture


class PhaseNativeRunnerTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        phase_fixture.PhaseAttributionTest.setUpClass()
        cls.fixture = phase_fixture.PhaseAttributionTest
        cls.contract = cls.fixture.contract
        cls.description = dict(cls.fixture.description)
        cls.description["resolved_path"] = "/fixture/phase-worker"

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        recovery_dir = self.root / "recovery"
        recovery_dir.mkdir()
        witness = {"fixture": "timing-proxy"}
        witness_sha256 = contract_api.sha256_json(witness)
        binary = self.description["binary_sha256"]
        source = self.description["source_git_commit"]
        recovery_freeze_sha256 = "b" * 64
        recovery_result_sha256 = "c" * 64
        recovery_receipt_sha256 = "d" * 64
        self.campaign = {
            "directory": str(recovery_dir.resolve()),
            "directory_identity": (
                recovery_dir.stat().st_dev, recovery_dir.stat().st_ino),
            "candidate_id": "two07",
            "candidate_arm": subject.phase_api.LEFT_ARM,
            "summary": {
                "source_git_commit": source,
                "worker_binary_sha256": binary,
                "timing_proxy_witness_sha256": witness_sha256,
                "summary_sha256": "a" * 64,
                "recovery_freeze_sha256": recovery_freeze_sha256,
                "recovery_result_sha256": recovery_result_sha256,
                "recovery_execution_receipt_sha256":
                    recovery_receipt_sha256,
            },
            "freeze": {
                "source_git_commit": source,
                "timing_proxy_witness_sha256": witness_sha256,
                "arms": [
                    {"binary_sha256": binary},
                    {},
                    {},
                    {
                        "arm": subject.phase_api.LEFT_ARM,
                        "arm_descriptor_sha256":
                            subject.phase_api.EXPECTED_ARMS[2][2],
                    },
                ],
            },
            "receipt": {
                "freeze_manifest_sha256": recovery_freeze_sha256,
                "result_stream_sha256": recovery_result_sha256,
                "receipt_sha256": recovery_receipt_sha256,
            },
            "timing_proxy_witness": witness,
        }
        self.recovery_dir = recovery_dir
        self.binding = subject._campaign_binding(self.campaign)
        self.bundle_counter = 0
        self.thermal = {
            "pid": 9000,
            "cpu": 127,
            "process_start_ticks": 1,
            "script_path": "/fixture/sampler.py",
            "script_sha256": "e" * 64,
            "csv_path": "/fixture/thermal.csv",
            "csv_device": 1,
            "csv_inode": 2,
            "window_csv_sha256": "f" * 64,
            "window_start_monotonic_ns": self.fixture.window_start_ns,
            "window_end_monotonic_ns": self.fixture.window_end_ns,
            "sample_count": 2,
            "cpu_tctl_max_millic": 65000,
            "dimm_max_millic": 43000,
            "dimm_read_errors": 0,
            "edac_ce_max": 0,
            "edac_ue_max": 0,
            "terminal_status": "ok",
        }

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _sampler(self):
        return {
            "schema": native_api.SAMPLER_SCHEMA,
            "pid": self.thermal["pid"],
            "cpu": self.thermal["cpu"],
            "process_start_ticks": self.thermal["process_start_ticks"],
            "script_path": self.thermal["script_path"],
            "script_sha256": self.thermal["script_sha256"],
            "csv_path": self.thermal["csv_path"],
            "csv_device": self.thermal["csv_device"],
            "csv_inode": self.thermal["csv_inode"],
            "window_start_monotonic_ns":
                self.thermal["window_start_monotonic_ns"],
            "window_end_monotonic_ns":
                self.thermal["window_end_monotonic_ns"],
            "terminal_status": "ok",
        }

    def _trace_worker(self, body: str, name: str) -> Path:
        path = self.root / name
        path.write_text(
            "#!{}\n{}\n".format(sys.executable, body), encoding="utf-8")
        path.chmod(0o755)
        return path

    @staticmethod
    def _host():
        return {
            "controller_cpu": 8,
            "hostname": "fixture-host",
            "kernel_release": "fixture-kernel",
            "machine": "x86_64",
            "system": "Linux",
        }

    @staticmethod
    def _write_object(path: Path, value) -> bytes:
        data = (contract_api.canonical_json(value) + "\n").encode("utf-8")
        path.write_bytes(data)
        return data

    def _bundle(self, include_summary: bool = True):
        self.bundle_counter += 1
        output = self.root / "phase-{}-{}".format(
            "complete" if include_summary else "prospective",
            self.bundle_counter)
        output.mkdir()
        trace_data = self.fixture.trace_data
        records = self.fixture.build_records()
        native_data = phase_fixture.canonical_jsonl(records)
        description = subject._public_phase_description(self.description)
        with mock.patch.object(
                subject.runner_api, "_host_identity",
                return_value=self._host()), mock.patch.object(
                subject.phase_api, "_deterministic_message_sha256",
                side_effect=phase_fixture.fake_message_sha256):
            _traces, trace_manifest_sha256 = \
                subject.phase_api.validate_trace_manifest(
                    self.contract, description, trace_data)
            freeze = subject._phase_freeze(
                self.contract, self.description, self.binding, list(range(8)), 8,
                trace_manifest_sha256, subject._hash_bytes(trace_data))
            payloads, metadata = subject.phase_api.validate_native_records(
                self.contract, description, trace_data, native_data,
                list(range(8)) * 3, self.fixture.runtime_workers,
                self.fixture.window_start_ns, self.fixture.window_end_ns)
            analysis = subject.phase_api.build_phase_analysis(payloads)
        freeze_bytes = self._write_object(
            output / subject.PHASE_SOURCE_NAMES["freeze"], freeze)
        (output / subject.PHASE_SOURCE_NAMES["trace"]).write_bytes(trace_data)
        (output / subject.PHASE_SOURCE_NAMES["native"]).write_bytes(native_data)
        analysis_bytes = self._write_object(
            output / subject.PHASE_SOURCE_NAMES["analysis"], analysis)
        sampler_bytes = self._write_object(
            output / subject.PHASE_SOURCE_NAMES["sampler"], self._sampler())
        execution = subject._execution_receipt(
            self.contract, freeze, freeze_bytes, trace_data, native_data,
            analysis_bytes, sampler_bytes, metadata, self.thermal, 8)
        execution_bytes = self._write_object(
            output / subject.PHASE_SOURCE_NAMES["execution"], execution)
        summary = subject._phase_summary(
            self.contract, output.resolve(), freeze, execution, analysis,
            self.binding, freeze_bytes, trace_data, native_data,
            analysis_bytes, execution_bytes, sampler_bytes)
        if include_summary:
            self._write_object(
                output / subject.PHASE_SOURCE_NAMES["summary"], summary)
        return output, freeze, execution, analysis, summary

    def _loader_patches(self, campaign=None):
        selected = self.campaign if campaign is None else campaign
        return (
            mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=selected),
            mock.patch.object(
                subject.phase_api, "_deterministic_message_sha256",
                side_effect=phase_fixture.fake_message_sha256),
            mock.patch.object(
                subject, "_validate_thermal", return_value=self.thermal),
            mock.patch.object(subject, "_validate_physical_layout"),
        )

    def test_exact_command_and_schedule_geometry(self) -> None:
        waves = subject._phase_job_waves()
        jobs = [job for wave in waves for job in wave]
        self.assertEqual(len(waves), 3)
        self.assertTrue(all(len(wave) == 8 for wave in waves))
        self.assertEqual([job.ordinal for job in jobs], list(range(0, 48, 2)))
        self.assertEqual([job.command() for job in jobs], [
            "P {} 0\n".format(value).encode("ascii") for value in range(24)
        ])
        with mock.patch.object(
                subject.runner_api, "_host_identity",
                return_value=self._host()):
            freeze = subject._phase_freeze(
                self.contract, self.description, self.binding,
                list(range(8)), 8, "1" * 64, "2" * 64)
        self.assertEqual(freeze["coordinate_cpus"], list(range(8)) * 3)
        self.assertEqual(freeze["commands"]["job_ordinals"],
                         list(range(0, 48, 2)))

    def test_raw_response_cap_precedes_json_parser(self) -> None:
        process = mock.Mock()
        process.poll.return_value = None
        process.stdout.fileno.return_value = 17
        worker = subject.runner_api.PersistentWorker(3, process, 1)
        job = subject.PhaseJob(0)
        selector = mock.Mock()
        selector.select.return_value = [
            (mock.Mock(data=worker), subject.runner_api.selectors.EVENT_READ),
        ]
        parser = mock.Mock(side_effect=AssertionError(
            "oversized phase response reached JSON parser"))
        line = b'{"value":"oversized"}\n'
        with mock.patch.object(
                subject.runner_api.selectors, "DefaultSelector",
                return_value=selector), mock.patch.object(
                subject.runner_api, "_write_worker"), mock.patch.object(
                subject.runner_api.os, "read", return_value=line), \
                mock.patch.object(
                    subject.runner_api, "_parse_canonical_line", parser), \
                self.assertRaisesRegex(
                    subject.runner_api.RunnerError, "exact byte cap"):
            subject.runner_api.run_job_batch(
                [worker], [job], 0, mock.Mock(), time.monotonic() + 1.0,
                mock.Mock(), maximum_response_bytes=16)
        parser.assert_not_called()

    def test_trace_stdout_cap_precedes_manifest_parser_and_reaps(self) \
            -> None:
        worker = self._trace_worker(
            "import os\nos.write(1, b'x' * {})".format(
                subject.phase_api.MAX_TRACE_BYTES + 4096),
            "trace-stdout-overflow.py")
        parser = mock.Mock(side_effect=AssertionError(
            "overflow reached phase manifest parser"))
        processes = []
        real_popen = subject.subprocess.Popen

        def capture(*args, **kwargs):
            process = real_popen(*args, **kwargs)
            processes.append(process)
            return process

        with mock.patch.object(
                subject.subprocess, "Popen", side_effect=capture), \
                mock.patch.object(
                    subject.phase_api, "validate_trace_manifest", parser), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError,
                    "stdout exceeds its exact byte cap"):
            subject._emit_phase_traces(
                self.contract, {**self.description,
                                "resolved_path": str(worker)},
                self.root, time.monotonic() + 3.0)
        parser.assert_not_called()
        self.assertEqual(len(processes), 1)
        self.assertIsNotNone(processes[0].poll())
        self.assertFalse((self.root / "phase-traces.jsonl").exists())

    def test_trace_stderr_flood_timeout_and_nonzero_are_reaped(self) -> None:
        cases = (
            (
                "stderr-cap",
                "import os\nos.write(2, b'e' * {})".format(
                    subject.MAX_PHASE_TRACE_STDERR_BYTES + 4096),
                3.0, "stderr exceeds its exact byte cap"),
            (
                "timeout", "import time\ntime.sleep(10)",
                0.1, "hard wall expired"),
            (
                "nonzero", "import os\nos.write(2, b'no')\nraise SystemExit(7)",
                3.0, "exited 7: no"),
        )
        real_popen = subject.subprocess.Popen
        for label, body, seconds, diagnostic in cases:
            with self.subTest(label=label):
                worker = self._trace_worker(
                    body, "trace-{}.py".format(label))
                processes = []

                def capture(*args, **kwargs):
                    process = real_popen(*args, **kwargs)
                    processes.append(process)
                    return process

                with mock.patch.object(
                        subject.subprocess, "Popen", side_effect=capture), \
                        self.assertRaisesRegex(
                            subject.PhaseNativeRunnerError, diagnostic):
                    subject._run_capped_phase_command(
                        [str(worker)], time.monotonic() + seconds,
                        "{} fixture".format(label))
                self.assertEqual(len(processes), 1)
                self.assertIsNotNone(processes[0].poll())

    def test_trace_selector_handoff_interrupt_owns_group_and_pipes(self) \
            -> None:
        class HandoffInterrupt(BaseException):
            pass

        primary = HandoffInterrupt("injected trace selector handoff")
        process = mock.Mock()
        process.pid = 4242
        process.returncode = None
        process.stdin = None
        process.stdout = mock.Mock()
        process.stderr = mock.Mock()
        process.wait.return_value = -15
        process.poll.return_value = -15
        signals = []

        def signal_group(pid, sig):
            signals.append((pid, sig))
            if sig == subject.signal.SIGTERM:
                raise OSError("injected TERM cleanup failure")

        with mock.patch.object(
                subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject.selectors, "DefaultSelector",
                    side_effect=primary), mock.patch.object(
                    subject.os, "killpg", side_effect=signal_group), \
                mock.patch.object(subject.time, "sleep"), self.assertRaises(
                        HandoffInterrupt) as raised:
            subject._run_capped_phase_command(
                ["/fixture/worker"], time.monotonic() + 1.0,
                "handoff fixture")
        self.assertIs(raised.exception, primary)
        self.assertIsInstance(raised.exception.__cause__,
                              subject.runner_api.RunnerError)
        self.assertEqual(signals, [
            (4242, subject.signal.SIGTERM),
            (4242, subject.signal.SIGKILL),
        ])
        self.assertEqual(process.wait.call_count, 1)
        process.stdout.close.assert_called_once_with()
        process.stderr.close.assert_called_once_with()

    def test_trace_cleanup_never_signals_reaped_pid_and_waits_last(self) \
            -> None:
        reaped = mock.Mock()
        reaped.pid = 5000
        reaped.returncode = 7
        reaped.stdin = None
        reaped.stdout = mock.Mock()
        reaped.stderr = mock.Mock()
        reaped.poll.return_value = 7
        with mock.patch.object(subject.os, "killpg") as killpg, \
                mock.patch.object(subject.time, "sleep") as sleep:
            subject._finish_capped_phase_process(
                reaped, None, successful=False, reaped=True, primary=None)
        killpg.assert_not_called()
        sleep.assert_not_called()
        reaped.wait.assert_not_called()

        live = mock.Mock()
        live.pid = 5001
        live.returncode = None
        live.stdin = None
        live.stdout = mock.Mock()
        live.stderr = mock.Mock()
        live.wait.return_value = -9
        live.poll.return_value = -9
        events = []

        def signal_group(_pid, sig):
            events.append(sig)

        def grace(_seconds):
            events.append("grace")

        def wait(*_args, **_kwargs):
            events.append("wait")
            return -9

        live.wait.side_effect = wait
        with mock.patch.object(
                subject.os, "killpg", side_effect=signal_group), \
                mock.patch.object(subject.time, "sleep", side_effect=grace):
            subject._finish_capped_phase_process(
                live, None, successful=False, reaped=False, primary=None)
        self.assertEqual(events, [
            subject.signal.SIGTERM, "grace", subject.signal.SIGKILL, "wait",
        ])
        self.assertEqual(live.wait.call_count, 1)

        handoff = mock.Mock()
        handoff.pid = 5002
        handoff.returncode = 7
        handoff.stdin = None
        handoff.stdout = mock.Mock()
        handoff.stderr = mock.Mock()
        handoff.poll.return_value = 7
        with mock.patch.object(subject.os, "killpg") as killpg, \
                mock.patch.object(subject.time, "sleep") as sleep:
            subject._finish_capped_phase_process(
                handoff, None, successful=False, reaped=False, primary=None)
        killpg.assert_not_called()
        sleep.assert_not_called()
        handoff.wait.assert_not_called()

    def test_cross_lane_wave_overlap_is_rejected(self) -> None:
        records = self.fixture.build_records()
        by_coordinate = {
            row["payload"]["coordinate_ordinal"]: row for row in records
        }
        first_lane_end = by_coordinate[0]["finished_monotonic_ns"]
        last_lane_end = max(
            by_coordinate[index]["finished_monotonic_ns"]
            for index in range(8))
        by_coordinate[8]["started_monotonic_ns"] = first_lane_end
        by_coordinate[8]["finished_monotonic_ns"] = first_lane_end + 1
        self.assertLess(
            by_coordinate[8]["started_monotonic_ns"], last_lane_end)
        with self.assertRaisesRegex(
                subject.PhaseNativeRunnerError, "wave barrier"):
            subject._validate_phase_wave_barriers(
                phase_fixture.canonical_jsonl(records))

    def test_worker_constructor_handoff_interrupt_is_cleaned(self) -> None:
        class HandoffInterrupt(BaseException):
            pass

        primary = HandoffInterrupt("injected phase worker handoff interrupt")
        process = mock.Mock()
        cleanup = []
        with mock.patch.object(
                subject.runner_api.subprocess, "Popen",
                return_value=process), mock.patch.object(
                subject.runner_api, "PersistentWorker",
                side_effect=primary), mock.patch.object(
                subject.runner_api, "terminate_workers",
                side_effect=lambda workers: cleanup.append(list(workers))), \
                self.assertRaises(HandoffInterrupt) as raised:
            subject.spawn_phase_workers(
                {"resolved_path": "/fixture/worker"}, [17],
                time.monotonic() + 1.0)
        self.assertIs(raised.exception, primary)
        self.assertEqual(len(cleanup), 1)
        self.assertEqual(len(cleanup[0]), 1)
        self.assertIs(cleanup[0][0].process, process)

    def test_cleanup_exhausts_workers_and_affinity_preserving_control(self) \
            -> None:
        class CleanupInterrupt(BaseException):
            pass

        primary = CleanupInterrupt("primary control")
        events = []
        with mock.patch.object(
                subject.runner_api, "terminate_workers",
                side_effect=lambda _workers: events.append("terminate")), \
                mock.patch.object(
                    subject.runner_api, "_restore_controller_affinity",
                    side_effect=lambda _affinity: events.append("restore")), \
                self.assertRaises(CleanupInterrupt) as raised:
            try:
                raise primary
            finally:
                subject._finish_phase_cleanup(
                    [mock.Mock()], False, True, {7}, sys.exc_info()[1])
        self.assertIs(raised.exception, primary)
        self.assertEqual(events, ["terminate", "restore"])

        events.clear()
        with mock.patch.object(
                subject.runner_api, "terminate_workers",
                side_effect=subject.runner_api.RunnerError("terminate")), \
                mock.patch.object(
                    subject.runner_api, "_restore_controller_affinity",
                    side_effect=lambda _affinity: events.append("restore")), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError, "terminate"):
            subject._finish_phase_cleanup(
                [mock.Mock()], False, True, {7}, None)
        self.assertEqual(events, ["restore"])

    def test_strict_loader_double_derives_and_accepts_complete_bundle(self) \
            -> None:
        output, _freeze, _execution, analysis, summary = self._bundle()
        original = subject._validate_phase_snapshots
        semantic = mock.Mock(wraps=original)
        with ExitStack() as stack:
            for patcher in self._loader_patches():
                stack.enter_context(patcher)
            stack.enter_context(mock.patch.object(
                subject, "_validate_phase_snapshots", semantic))
            loaded = subject.load_completed_phase_screen(
                self.contract, output, self.recovery_dir)
        self.assertEqual(semantic.call_count, 2)
        self.assertEqual(loaded["summary"], summary)
        self.assertEqual(loaded["analysis"], analysis)
        self.assertEqual(loaded["freeze"]["coordinate_cpus"],
                         list(range(8)) * 3)

    def test_loader_rejects_terminal_recovery_input_substitution(self) -> None:
        output, _freeze, _execution, _analysis, _summary = self._bundle()
        substituted = copy.deepcopy(self.campaign)
        substituted["summary"]["recovery_result_sha256"] = "9" * 64
        substituted["receipt"]["result_stream_sha256"] = "9" * 64
        with mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                side_effect=(self.campaign, substituted)), mock.patch.object(
                subject.phase_api, "_deterministic_message_sha256",
                side_effect=phase_fixture.fake_message_sha256), \
                mock.patch.object(
                    subject, "_validate_thermal", return_value=self.thermal), \
                mock.patch.object(subject, "_validate_physical_layout"), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError,
                    "recovery input changed"):
            subject.load_completed_phase_screen(
                self.contract, output, self.recovery_dir)

    def test_loader_rejects_extra_artifact_name(self) -> None:
        output, _freeze, _execution, _analysis, _summary = self._bundle()
        (output / "unbound-extra.json").write_bytes(b"{}\n")
        with mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError,
                    "too many artifacts"):
            subject.load_completed_phase_screen(
                self.contract, output, self.recovery_dir)

    def test_exact_roster_scan_stops_at_expected_plus_one_and_closes(self) \
            -> None:
        expected = sorted(subject.PHASE_DEPENDENCY_NAMES.values())

        class Entry:
            def __init__(self, name):
                self.name = name

            def is_file(self, follow_symlinks=True):
                self_test.assertFalse(follow_symlinks)
                return True

        class BoundedScan:
            def __init__(self):
                self.values = [Entry(name) for name in expected] + [
                    Entry("first-extra"), Entry("must-not-be-read")]
                self.index = 0
                self.closed = False

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                self.closed = True

            def __iter__(self):
                return self

            def __next__(self):
                if self.index >= len(self.values):
                    raise StopIteration
                if self.index > len(expected):
                    raise AssertionError("roster scan consumed beyond overflow")
                value = self.values[self.index]
                self.index += 1
                return value

        self_test = self
        scan = BoundedScan()
        with mock.patch.object(
                subject.os, "scandir", return_value=scan), mock.patch.object(
                subject.os, "listdir",
                side_effect=AssertionError("unbounded listdir used")), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError, "too many artifacts"):
            subject._require_exact_phase_names(17, complete=False)
        self.assertEqual(scan.index, len(expected) + 1)
        self.assertTrue(scan.closed)

    def test_loader_rejects_worker_siblings_and_controller_overlap(self) \
            -> None:
        output, _freeze, _execution, _analysis, _summary = self._bundle()

        def sibling_topology(cpu, _root):
            if cpu in (0, 1):
                return (0, 0)
            return (0, cpu)

        with mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), mock.patch.object(
                subject.runner_api, "_cpu_topology",
                side_effect=sibling_topology), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError, "isolated physical cores"):
            subject.load_completed_phase_screen(
                self.contract, output, self.recovery_dir)

        def controller_overlap(cpu, _root):
            return (0, 0) if cpu in (0, 8) else (0, cpu)

        with mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), mock.patch.object(
                subject.runner_api, "_cpu_topology",
                side_effect=controller_overlap), \
                self.assertRaisesRegex(
                    subject.PhaseNativeRunnerError, "isolated physical cores"):
            subject.load_completed_phase_screen(
                self.contract, output, self.recovery_dir)

    def test_commit_happy_path_and_prelink_dependency_race(self) -> None:
        output, _freeze, _execution, _analysis, summary = self._bundle(
            include_summary=False)
        validated = {"summary": summary}
        events = []
        real_roster = subject._require_exact_phase_names
        real_link = subject.runner_api._link_unnamed_completion_marker

        def roster(*args, **kwargs):
            events.append("roster")
            return real_roster(*args, **kwargs)

        def link(*args, **kwargs):
            events.append("link")
            return real_link(*args, **kwargs)

        finish = mock.Mock(side_effect=lambda: events.append("finish"))
        with mock.patch.object(
                subject, "_validate_phase_snapshots",
                return_value=validated) as semantic, mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), mock.patch.object(
                subject.runner_api, "_git_head",
                return_value=self.description["source_git_commit"]), \
                mock.patch.object(
                    subject, "_require_exact_phase_names",
                    side_effect=roster), mock.patch.object(
                    subject.runner_api, "_link_unnamed_completion_marker",
                    side_effect=link):
            subject._commit_completed_phase_screen(
                self.contract, output, self.recovery_dir, self.binding,
                summary, self.description["source_git_commit"],
                time.monotonic() + 5.0, finish)
        self.assertEqual(semantic.call_count, 2)
        finish.assert_called_once_with()
        self.assertEqual(events[-3:], ["finish", "roster", "link"])
        self.assertEqual(
            subject._parse_object_bytes(
                (output / "run-summary.json").read_bytes(), "summary"),
            summary)

        raced, _freeze, _execution, _analysis, raced_summary = self._bundle(
            include_summary=False)
        real_check = subject._require_pinned_phase_unchanged

        def mutate_then_check(names, directory_fd, descriptors, fingerprints):
            (raced / subject.PHASE_SOURCE_NAMES["analysis"]).write_bytes(
                b"substituted\n")
            return real_check(names, directory_fd, descriptors, fingerprints)

        with mock.patch.object(
                subject, "_validate_phase_snapshots",
                return_value={"summary": raced_summary}), mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), mock.patch.object(
                subject.runner_api, "_git_head",
                return_value=self.description["source_git_commit"]), \
                mock.patch.object(
                    subject, "_require_pinned_phase_unchanged",
                    side_effect=mutate_then_check), self.assertRaisesRegex(
                        subject.PhaseNativeRunnerError, "changed"):
            subject._commit_completed_phase_screen(
                self.contract, raced, self.recovery_dir, self.binding,
                raced_summary, self.description["source_git_commit"],
                time.monotonic() + 5.0, mock.Mock())
        self.assertFalse((raced / "run-summary.json").exists())

    def test_postlink_lost_ack_preserves_published_summary(self) -> None:
        output, _freeze, _execution, _analysis, summary = self._bundle(
            include_summary=False)
        real_link = subject.runner_api._link_unnamed_completion_marker

        class LostAck(BaseException):
            pass

        lost = LostAck("injected post-link lost acknowledgement")

        def link_then_interrupt(*args, **kwargs):
            real_link(*args, **kwargs)
            raise lost

        with mock.patch.object(
                subject, "_validate_phase_snapshots",
                return_value={"summary": summary}), mock.patch.object(
                subject.recovery_api, "load_completed_campaign",
                return_value=self.campaign), mock.patch.object(
                subject.runner_api, "_git_head",
                return_value=self.description["source_git_commit"]), \
                mock.patch.object(
                    subject.runner_api, "_link_unnamed_completion_marker",
                    side_effect=link_then_interrupt), self.assertRaises(
                        LostAck) as raised:
            subject._commit_completed_phase_screen(
                self.contract, output, self.recovery_dir, self.binding,
                summary, self.description["source_git_commit"],
                time.monotonic() + 5.0, mock.Mock())
        self.assertIs(raised.exception, lost)
        self.assertTrue((output / "run-summary.json").is_file())

    def test_run_happy_path_restores_affinity_before_terminal_checks(self) \
            -> None:
        output = self.root / "run-output"
        records = self.fixture.build_records()
        native_data = phase_fixture.canonical_jsonl(records)
        trace_data = self.fixture.trace_data
        runtime_rows = [{
            "cpu": cpu,
            "pid": identity[0],
            "process_start_ticks": identity[1],
        } for cpu, identity in self.fixture.runtime_workers.items()]
        workers = [SimpleNamespace(cpu=cpu) for cpu in range(8)]
        args = SimpleNamespace(
            deadline_seconds=60.0,
            contract=contract_api.DEFAULT_CONTRACT,
            recovery_dir=self.recovery_dir,
            worker=Path("/fixture/worker"),
            sampler_pid=9000,
            sampler_cpu=127,
            sampler_script=Path("/fixture/sampler.py"),
            sampler_csv=Path("/fixture/thermal.csv"),
            cpus="0,1,2,3,4,5,6,7",
            controller_cpu=8,
            output_dir=output,
        )
        events = []
        git_calls = 0
        recovery_calls = 0

        def git_head(_deadline):
            nonlocal git_calls
            git_calls += 1
            events.append("git{}".format(git_calls))
            return self.description["source_git_commit"]

        def load_binding(_contract, _path):
            nonlocal recovery_calls
            recovery_calls += 1
            events.append("recovery{}".format(recovery_calls))
            return self.campaign, self.binding

        def emit_traces(_contract, _description, target, _deadline):
            path = target / subject.PHASE_SOURCE_NAMES["trace"]
            path.write_bytes(trace_data)
            _rows, digest = subject.phase_api.validate_trace_manifest(
                self.contract,
                subject._public_phase_description(self.description),
                trace_data)
            return path, trace_data, digest

        def run_jobs(_contract, _freeze, _description, _trace, _workers,
                     target, _window_start, _deadline):
            path = target / subject.PHASE_SOURCE_NAMES["native"]
            path.write_bytes(native_data)
            return path, self.fixture.window_end_ns - 100

        def sampler_attestation(*call_args):
            self._write_object(Path(call_args[-1]), self._sampler())

        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(
                subject.contract_api, "load_contract",
                return_value=self.contract))
            stack.enter_context(mock.patch.object(
                subject, "_load_recovery_binding", side_effect=load_binding))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "describe_worker",
                return_value=self.description))
            stack.enter_context(mock.patch.object(
                subject.os, "sched_getaffinity", return_value=set(range(10))))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "select_cpu_layout",
                return_value=(list(range(8)), 8)))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_preflight_sampler"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_git_head", side_effect=git_head))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_require_worker_source_commit"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_create_output_dir",
                side_effect=lambda path: (
                    path.mkdir(), path.resolve(strict=True))[1]))
            stack.enter_context(mock.patch.object(
                subject, "_emit_phase_traces", side_effect=emit_traces))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_host_identity",
                return_value=self._host()))
            stack.enter_context(mock.patch.object(
                subject, "spawn_phase_workers", return_value=workers))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_pin_controller",
                side_effect=lambda _cpu: events.append("pin")))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "choose_new_sampler_start",
                return_value=self.fixture.window_start_ns))
            stack.enter_context(mock.patch.object(
                subject, "_run_phase_jobs", side_effect=run_jobs))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_wait_for_sampler_sample",
                return_value=self.fixture.window_end_ns))
            stack.enter_context(mock.patch.object(
                subject.native_api, "write_sampler_attestation",
                side_effect=sampler_attestation))
            stack.enter_context(mock.patch.object(
                subject, "_runtime_worker_map",
                return_value=(self.fixture.runtime_workers, runtime_rows)))
            stack.enter_context(mock.patch.object(
                subject.phase_api, "_deterministic_message_sha256",
                side_effect=phase_fixture.fake_message_sha256))
            stack.enter_context(mock.patch.object(
                subject, "_validate_thermal", return_value=self.thermal))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "quit_workers",
                side_effect=lambda _workers, _deadline: events.append("quit")))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_restore_controller_affinity",
                side_effect=lambda _affinity: events.append("restore")))
            commit = stack.enter_context(mock.patch.object(
                subject, "_commit_completed_phase_screen",
                side_effect=lambda *_args: events.append("commit")))
            summary = subject.run_phase_screen(args)
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(recovery_calls, 2)
        self.assertEqual(git_calls, 3)
        self.assertLess(events.index("quit"), events.index("restore"))
        self.assertLess(events.index("restore"), events.index("git3"))
        self.assertLess(events.index("git3"), events.index("recovery2"))
        self.assertLess(events.index("recovery2"), events.index("commit"))
        commit.assert_called_once()


if __name__ == "__main__":
    unittest.main()
