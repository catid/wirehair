#!/usr/bin/env python3
"""Bounded tests for the fixed native WH2 short-screen controller."""

from __future__ import annotations

import copy
from contextlib import ExitStack, redirect_stderr
import csv
import dis
import errno
import hashlib
import inspect
from io import FileIO, StringIO
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_run_native_short_screen as subject


FAKE_WORKER = r'''#!/usr/bin/env python3
import hashlib
import json
import os
from pathlib import Path
import sys
import time

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))

def digest(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()

def event(state, command, cpu):
    path = os.environ.get("WH2_FAKE_WORKER_LOG")
    if not path:
        return
    value = canonical({"command": command, "cpu": cpu, "state": state}) + "\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    try:
        os.write(descriptor, value.encode("utf-8"))
    finally:
        os.close(descriptor)

binary = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
arms = [
    {"arm": "wirehair2_head", "codec": "wirehair2_certified",
     "arm_descriptor_sha256": digest("head")},
    {"arm": "wirehair1", "codec": "wirehair1",
     "arm_descriptor_sha256": digest("wh1")},
    {"arm": "wirehair2_dense_two07_basis_v1",
     "codec": "wirehair2_experiment",
     "arm_descriptor_sha256":
         "9527f200ad38c7eec6502b2f768fdd67"
         "b92787fb227eed3d7616274ffc2df388"},
]
if sys.argv[1:] == ["--describe"]:
    print(canonical({"arms": arms, "binary_sha256": binary,
                     "source_git_commit": "1111111111111111111111111111111111111111",
                     "schema": "wirehair.wh2.native-worker-description.v1"}),
          flush=True)
    raise SystemExit(0)
if len(sys.argv) == 3 and sys.argv[1] == "--emit-traces":
    print(canonical({"kind": sys.argv[2]}), flush=True)
    raise SystemExit(0)
if len(sys.argv) != 3 or sys.argv[1] != "--worker":
    raise SystemExit(2)
cpu = int(sys.argv[2])
os.sched_setaffinity(0, {cpu})
for raw in sys.stdin:
    command = raw.rstrip("\n")
    if command == "Q":
        raise SystemExit(0)
    event("start", command, cpu)
    if command == os.environ.get("WH2_FAKE_WORKER_FAIL"):
        sys.stderr.write("injected worker failure\n")
        sys.stderr.flush()
        raise SystemExit(7)
    time.sleep(float(os.environ.get("WH2_FAKE_WORKER_DELAY", "0.002")))
    event("done", command, cpu)
    print(canonical({
        "command": command,
        "cpu": cpu,
        "finished_monotonic_ns": time.monotonic_ns(),
        "pid": os.getpid(),
    }), flush=True)
raise SystemExit(9)
'''


def sampler_row(monotonic_s: float) -> list:
    return [
        "2026-08-02T00:00:00Z", "{:.6f}".format(monotonic_s),
        "75.0", "4200.0", "60.0",
        "40.0", "40.0", "40.0", "40.0",
        "40.0", "40.0", "40.0", "40.0",
        "0", "1.0", "1.0", "1.0", "0", "0",
    ]


class _ManagedPauseFixture:
    """Small exact-owner stand-in for full-run cleanup control-flow tests."""

    def __init__(self) -> None:
        self.records = [{
            "pid": 101,
            "process_start_ticks": 1001,
            "cpu_affinity": [9],
        }]
        self.pause_calls = 0
        self.resume_calls = 0
        self.stopped = False
        self.resumed = False

    def pause(self, _pids, _sibling_cpus, _deadline) -> None:
        self.pause_calls += 1
        self.stopped = True

    def resume(self, _deadline) -> None:
        self.resume_calls += 1
        self.stopped = False
        self.resumed = True

    @property
    def cleanup_complete(self) -> bool:
        return self.resumed


class NativeRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.qualification_counter = 0

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run_timing_cleanup_failure(
            self, managed_paused: _ManagedPauseFixture,
            body_failure: BaseException,
            restore_affinity: mock.Mock,
            finish_hard_wall=None):
        """Reach T with mocks, own a stop, then fail the first guard check."""
        contract = contract_api.load_contract()
        qualification = mock.Mock(map_sha256="a" * 64)
        description = {
            "binary_sha256": "b" * 64,
            "source_git_commit": "1" * 40,
        }
        guard = mock.Mock()
        guard.check.side_effect = body_failure
        args = mock.Mock(
            deadline_seconds=60.0, contract=None, worker=Path("worker"),
            sampler_cpu=99, sampler_pid=123,
            sampler_script=Path("s.py"),
            sampler_csv=Path("thermal.csv"), cpus=None,
            controller_cpu=None, pause_sibling_pids="101",
            output_dir=self.root / "cleanup-control-flow")
        patches = (
            mock.patch.object(
                contract_api, "load_contract", return_value=contract),
            mock.patch.object(
                subject, "describe_worker", return_value=description),
            mock.patch.object(
                subject.os, "sched_getaffinity",
                return_value=set(range(12))),
            mock.patch.object(
                subject, "_cpu_topology",
                side_effect=lambda cpu, _root: (0, cpu)),
            mock.patch.object(subject, "_development_timing_shape"),
            mock.patch.object(subject, "_preflight_sampler"),
            mock.patch.object(subject, "_git_head", return_value="1" * 40),
            mock.patch.object(subject, "_require_worker_source_commit"),
            mock.patch.object(
                subject, "_create_output_dir", return_value=args.output_dir),
            mock.patch.object(
                subject, "_timing_qualification_controls", return_value=[]),
            mock.patch.object(
                subject, "choose_new_sampler_start", return_value=1000000000),
            mock.patch.object(subject, "spawn_workers", return_value=[]),
            mock.patch.object(
                subject, "run_timing_qualification",
                return_value=(self.root / "qualification-native.jsonl",
                              1500000000)),
            mock.patch.object(
                native_api, "assemble_timing_qualification",
                return_value=(qualification, {}, "c" * 64)),
            mock.patch.object(subject, "quit_workers"),
            mock.patch.object(
                subject, "_wait_for_sampler_sample",
                return_value=2000000000),
            mock.patch.object(native_api, "write_sampler_attestation"),
            mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification),
            mock.patch.object(
                subject, "select_cpu_layout",
                return_value=(list(range(8)), 8)),
            mock.patch.object(
                native_api, "publish_timing_trace_manifest",
                return_value="d" * 64),
            mock.patch.object(
                subject, "write_development_timing_freeze",
                return_value={"trace_manifest_sha256": "d" * 64}),
            mock.patch.object(subject, "_pin_controller"),
            mock.patch.object(
                subject, "_restore_controller_affinity",
                side_effect=restore_affinity),
            mock.patch.object(
                native_api, "timing_sibling_cpus", return_value=[9]),
            mock.patch.object(
                subject, "ManagedPausedSiblingProcesses",
                return_value=managed_paused),
            mock.patch.object(
                subject, "SiblingIsolationGuard", return_value=guard),
        )
        with ExitStack() as stack:
            for patcher in patches:
                stack.enter_context(patcher)
            return subject.run_short_screen(args, finish_hard_wall)

    def _inner_cleanup_entry_offset(self) -> tuple:
        """Return and structurally verify the direct timing-cleanup entry."""
        lines, first_line = inspect.getsourcelines(subject.run_short_screen)
        capture_index = next(
            index for index, line in enumerate(lines)
            if "cleanup_primary = sys.exc_info()[1]" in line)
        capture_line = first_line + capture_index
        instructions = list(dis.get_instructions(subject.run_short_screen))
        candidates = [
            (index, instruction)
            for index, instruction in enumerate(instructions)
            if instruction.starts_line == capture_line and
            instruction.opname == "LOAD_GLOBAL" and
            instruction.argval == "sys"
        ]
        self.assertTrue(candidates)
        index, target = max(
            candidates, key=lambda value: value[1].offset)
        if sys.version_info >= (3, 11):
            # The exception-path copy must begin immediately after the finally
            # unwinder state.  A nested try would insert the unsafe NOP here.
            self.assertEqual(instructions[index - 1].opname, "PUSH_EXC_INFO")
            self.assertTrue(any(
                entry.start <= target.offset < entry.end
                for entry in dis.Bytecode(
                    subject.run_short_screen).exception_entries))
        else:
            # Python 3.8 uses CALL_FINALLY and jumps directly to this opcode.
            self.assertTrue(target.is_jump_target)
        return target.offset, target.opname

    def _topology(self, values: Mapping[int, tuple]) -> Path:
        root = self.root / "sysfs"
        for cpu, value in values.items():
            package, core = value[:2]
            llc = value[2] if len(value) == 3 else core
            path = root / "cpu{}".format(cpu) / "topology"
            path.mkdir(parents=True)
            (path / "physical_package_id").write_text(
                str(package) + "\n", encoding="ascii")
            (path / "core_id").write_text(
                str(core) + "\n", encoding="ascii")
            # Deliberately use index7 so the reader cannot assume index3.
            for index, level, cache_type, cache_id in (
                    (0, 1, "Data", core), (7, 3, "Unified", llc)):
                cache = root / "cpu{}".format(cpu) / "cache" / \
                    "index{}".format(index)
                cache.mkdir(parents=True)
                (cache / "level").write_text(
                    str(level) + "\n", encoding="ascii")
                (cache / "type").write_text(
                    cache_type + "\n", encoding="ascii")
                (cache / "id").write_text(
                    str(cache_id) + "\n", encoding="ascii")
        return root

    def _timing_contract(self, repetitions: int) -> Mapping[str, Any]:
        contract = copy.deepcopy(contract_api.load_contract())
        domain = contract["timing"]["domains"]["development"]
        original_repetitions = domain["paired_repetitions"]
        cells_per_replicate = domain["expected_cells"] // \
            original_repetitions
        domain["paired_repetitions"] = repetitions
        domain["independent_rounds"] = \
            repetitions // subject.TIMING_WORKER_COUNT
        domain["expected_cells"] = cells_per_replicate * repetitions
        return contract

    def _fake_worker(self) -> Path:
        path = self.root / "fake_worker.py"
        path.write_text(FAKE_WORKER, encoding="utf-8")
        path.chmod(0o755)
        return path

    def _description(self, path: Path) -> Mapping[str, Any]:
        return subject.describe_worker(path, time.monotonic() + 5.0)

    def _qualification(
            self, contract: Optional[Mapping[str, Any]] = None,
            retry_offsets: Optional[Sequence[int]] = None,
            description: Optional[Mapping[str, Any]] = None,
            ) -> contract_api.TimingQualification:
        candidate = contract_api.load_contract() if contract is None else \
            contract
        base_cells = list(contract_api.iter_timing_base_cells(
            candidate, "development"))
        offsets = [0] * len(base_cells) if retry_offsets is None else \
            list(retry_offsets)
        if len(offsets) != len(base_cells):
            raise AssertionError("qualification fixture offset cardinality")
        audit_rows: List[Dict[str, Any]] = []
        selected_traces: List[str] = []
        for ordinal, (base_cell, selected) in enumerate(zip(
                base_cells, offsets)):
            for retry in range(selected + 1):
                terminal = retry == selected
                trace = hashlib.sha256(
                    "runner qualification:{}:{}".format(
                        ordinal, retry).encode("ascii")).hexdigest()
                audit_rows.append({
                    "ordinal": ordinal,
                    "base_cell_sha256":
                        contract_api.sha256_json(base_cell),
                    "loss_retry_offset": retry,
                    "loss_seed": contract_api._qualified_timing_loss_seed(
                        base_cell["base_loss_seed"], retry),
                    "trace_sha256": trace,
                    "wirehair2_head_outcome":
                        "success" if terminal else "need_more_at_bound",
                    "wirehair2_head_decoded_extra":
                        4 if terminal else None,
                    "wirehair1_outcome": "success",
                    "wirehair1_decoded_extra": 0,
                })
                if terminal:
                    selected_traces.append(trace)
        self.qualification_counter += 1
        stem = "qualification-{}".format(self.qualification_counter)
        audit_path = self.root / (stem + ".jsonl")
        audit_path.write_bytes("".join(
            contract_api.canonical_json(row) + "\n" for row in audit_rows
        ).encode("utf-8"))
        described = {
            value["arm"]: value for value in
            (description or {}).get("arms", [])
        }
        binary = (description or {}).get(
            "binary_sha256", hashlib.sha256(b"runner binary").hexdigest())
        controls = []
        for arm, scope in (
                ("wirehair2_head", "decoder_solve"),
                ("wirehair1", "receive_to_success")):
            controls.append({
                "arm": arm,
                "scope": scope,
                "binary_sha256": binary,
                "arm_descriptor_sha256": described.get(arm, {}).get(
                    "arm_descriptor_sha256",
                    hashlib.sha256(("descriptor:" + arm).encode(
                        "ascii")).hexdigest()),
                "construction_policy":
                    "not_applicable" if arm == "wirehair1" else "raw_base",
                "repair_map_sha256": "0" * 64,
            })
        value = {
            "schema": contract_api.TIMING_QUALIFICATION_MAP_SCHEMA,
            "contract_sha256": contract_api.contract_sha256(candidate),
            "phase": "development",
            "source_git_commit": "1" * 40,
            "base_domain_sha256": candidate["timing"]["domains"]
                ["development"]["base_domain_sha256"],
            "qualified_domain_sha256":
                contract_api._timing_domain_sha256_from_offsets(
                    candidate, "development", offsets),
            "entry_kind": contract_api.TIMING_QUALIFICATION_ENTRY_KIND,
            "controls": controls,
            "qualification_audit_sha256":
                contract_api.timing_qualification_audit_sha256(
                    candidate, "development", audit_path),
            "selected_trace_roster_sha256":
                contract_api.timing_selected_trace_roster_sha256(
                    selected_traces),
            "retry_offsets": offsets,
        }
        map_path = self.root / (stem + ".json")
        map_path.write_bytes(
            (contract_api.canonical_json(value) + "\n").encode("utf-8"))
        return contract_api.load_timing_qualification_map(
            candidate, "development", map_path, audit_path,
            contract_api.timing_qualification_map_sha256(value))

    def _worker_cpus(self) -> list:
        values = sorted(os.sched_getaffinity(0))
        if len(values) < 2:
            self.skipTest("two logical CPUs are required for worker fixture")
        return values[:2]

    def _validator(self, value, _line, worker, job) -> int:
        self.assertEqual(set(value), {
            "command", "cpu", "finished_monotonic_ns", "pid",
        })
        self.assertEqual(value["command"], job.command().decode("ascii").strip())
        self.assertEqual(value["cpu"], worker.cpu)
        self.assertEqual(value["pid"], worker.process.pid)
        self.assertIs(type(value["finished_monotonic_ns"]), int)
        return value["finished_monotonic_ns"]

    def _completed_timing_fixture(self, stem: str = "completed-timing"):
        contract = contract_api.load_contract()
        qualification = self._qualification(contract)
        directory = self.root / stem
        directory.mkdir()
        binary = "a" * 64
        trace_sha256 = "b" * 64
        freeze = {
            "source_git_commit": "1" * 40,
            "arm_roster": [value[0] for value in subject.EXPECTED_ARMS],
            "trace_manifest_sha256": trace_sha256,
            "host_identity": {"controller_cpu": 8},
            "arms": [
                {
                    "arm": arm,
                    "codec": codec,
                    "binary_sha256": binary,
                    "arm_descriptor_sha256":
                        subject.TIMING_CANDIDATE_DESCRIPTOR_SHA256
                        if index == 2 else hashlib.sha256(
                            arm.encode("ascii")).hexdigest(),
                }
                for index, (arm, codec) in enumerate(subject.EXPECTED_ARMS)
            ],
        }
        timing_summary = {
            "schema": contract_api.SCHEMA + ".timing-summary.v1",
            "phase": "development",
        }
        thermal = {
            "sample_count": 2,
            "cpu_tctl_max_millic": 60000,
            "dimm_max_millic": 40000,
        }
        qualification_thermal = {
            "sample_count": 3,
            "cpu_tctl_max_millic": 61000,
            "dimm_max_millic": 41000,
        }
        execution_receipt = {
            "record_count": 25344,
            "worker_cpus": list(range(8)),
            "qualification_worker_cpus": list(range(12)),
            "freeze_manifest_sha256": "c" * 64,
            "result_stream_sha256": "d" * 64,
            "receipt_sha256": "e" * 64,
            "validator_summary_sha256":
                contract_api.sha256_json(timing_summary),
            "timing_qualification_execution_receipt_sha256": "f" * 64,
            "thermal": thermal,
            "qualification_thermal": qualification_thermal,
        }
        offsets = list(qualification.retry_offsets)
        run_summary = {
            "schema": subject.SUMMARY_SCHEMA,
            "status": "complete",
            "output_dir": str(directory.resolve()),
            "source_git_commit": freeze["source_git_commit"],
            "contract_sha256": contract_api.contract_sha256(contract),
            "worker_binary_sha256": binary,
            "controller_cpu": 8,
            "worker_cpus": execution_receipt["worker_cpus"],
            "qualification_worker_cpus":
                execution_receipt["qualification_worker_cpus"],
            "timing_qualification_map_sha256": qualification.map_sha256,
            "timing_qualification_execution_receipt_sha256": "f" * 64,
            "timing_qualified_domain_sha256":
                qualification.qualified_domain_sha256,
            "qualification_attempt_count": len(offsets) + sum(offsets),
            "qualification_retried_cell_count": 0,
            "qualification_max_retry_offset": 0,
            "qualification_sum_retry_offsets": 0,
            "timing_records": execution_receipt["record_count"],
            "timing_trace_manifest_sha256": trace_sha256,
            "timing_freeze_sha256":
                execution_receipt["freeze_manifest_sha256"],
            "timing_architecture_artifact_sha256":
                contract_api.architecture_artifact_sha256(freeze),
            "timing_result_sha256":
                execution_receipt["result_stream_sha256"],
            "timing_execution_receipt_sha256":
                execution_receipt["receipt_sha256"],
            "timing_validator_summary_sha256":
                execution_receipt["validator_summary_sha256"],
            "thermal_samples": thermal["sample_count"],
            "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "dimm_max_millic": thermal["dimm_max_millic"],
            "timing_thermal_samples": thermal["sample_count"],
            "timing_cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "timing_dimm_max_millic": thermal["dimm_max_millic"],
            "qualification_thermal_samples":
                qualification_thermal["sample_count"],
            "qualification_cpu_tctl_max_millic":
                qualification_thermal["cpu_tctl_max_millic"],
            "qualification_dimm_max_millic":
                qualification_thermal["dimm_max_millic"],
            "overall_cpu_tctl_max_millic": 61000,
            "overall_dimm_max_millic": 41000,
        }
        run_summary["summary_sha256"] = \
            contract_api.sha256_json(run_summary)
        names = (
            "run-summary.json", "timing-freeze.json",
            "timing-traces.jsonl", "timing-native-results.jsonl",
            "timing-results.jsonl", "timing-execution.json",
            "sampler-attestation.json", "timing-qualification-map.json",
            "timing-qualification-audit.jsonl",
            "timing-qualification-native.jsonl",
            "qualification-sampler-attestation.json",
            "timing-qualification-execution.json",
        )
        for name in names:
            (directory / name).write_bytes(b"{}\n")
        (directory / "run-summary.json").write_bytes(
            (contract_api.canonical_json(run_summary) + "\n").encode("utf-8"))
        validated = {
            "summary": timing_summary,
            "execution_receipt": execution_receipt,
        }
        return (contract, qualification, directory, freeze, run_summary,
                validated)

    def test_default_worker_name_matches_cmake_target(self) -> None:
        stderr = StringIO()
        with redirect_stderr(stderr):
            with mock.patch.object(sys, "argv", [
                    "runner", "--output-dir", "out", "--sampler-pid", "1",
                    "--sampler-cpu", "2", "--sampler-script", "s.py",
                    "--sampler-csv", "s.csv"]), \
                    mock.patch.object(subject, "run_short_screen",
                                      side_effect=subject.RunnerError("stop")) \
                    as run:
                self.assertEqual(subject.main(), 1)
        self.assertIn("stop", stderr.getvalue())
        self.assertEqual(
            run.call_args.args[0].worker,
            Path("build/wirehair_wh2_contract_worker"))

    def test_hard_wall_interrupts_synchronous_work(self) -> None:
        started = time.monotonic()
        with self.assertRaisesRegex(subject.RunnerError, "hard wall expired"):
            with subject._hard_wall(0.05):
                time.sleep(2.0)
        self.assertLess(time.monotonic() - started, 1.0)

    def test_hard_wall_finish_leaves_no_fallible_context_teardown(self) -> None:
        marker = self.root / "finished-marker"
        timer_calls = []

        def set_timer(which, seconds, interval=0.0):
            timer_calls.append((which, seconds, interval))
            if len(timer_calls) > 2:
                raise OSError(errno.EIO, "injected post-finish teardown")
            return (0.0, 0.0)

        def finish_inside_commit(_args, finish_hard_wall):
            marker.write_bytes(b"complete\n")
            finish_hard_wall()
            return {"status": "complete"}

        argv = [
            "--output-dir", str(self.root / "run"),
            "--sampler-pid", "1", "--sampler-cpu", "2",
            "--sampler-script", "s.py", "--sampler-csv", "s.csv",
        ]
        with mock.patch.object(
                subject.signal, "getsignal", return_value=signal.SIG_DFL), \
                mock.patch.object(subject.signal, "signal"), \
                mock.patch.object(
                    subject.signal, "setitimer", side_effect=set_timer), \
                mock.patch.object(
                    subject, "run_short_screen", side_effect=finish_inside_commit), \
                mock.patch("builtins.print"):
            self.assertEqual(subject.main(argv), 0)
        self.assertTrue(marker.exists())
        self.assertEqual(len(timer_calls), 2)

    def test_inner_cleanup_entry_interrupt_runs_full_outer_cleanup(self) \
            -> None:
        managed = _ManagedPauseFixture()
        restore_affinity = mock.Mock()
        body_failure = subject.RunnerError("injected timing body failure")
        cleanup_interrupt = subject.RunnerError(
            "injected first inner-finally opcode interrupt")
        target_offset, target_opname = self._inner_cleanup_entry_offset()
        injected = [False]
        if hasattr(sys, "monitoring"):
            monitoring = sys.monitoring
            tool_id = next(
                candidate for candidate in range(5, -1, -1)
                if monitoring.get_tool(candidate) is None)
            monitoring.use_tool_id(tool_id, "wh2 cleanup-entry test")

            def instruction(_code, offset):
                if offset == target_offset and not injected[0]:
                    injected[0] = True
                    raise cleanup_interrupt

            monitoring.register_callback(
                tool_id, monitoring.events.INSTRUCTION, instruction)
            monitoring.set_local_events(
                tool_id, subject.run_short_screen.__code__,
                monitoring.events.INSTRUCTION)
            try:
                with self.assertRaises(subject.RunnerError) as raised:
                    self._run_timing_cleanup_failure(
                        managed, body_failure, restore_affinity)
            finally:
                monitoring.set_local_events(
                    tool_id, subject.run_short_screen.__code__, 0)
                monitoring.register_callback(
                    tool_id, monitoring.events.INSTRUCTION, None)
                monitoring.free_tool_id(tool_id)
        else:
            previous_trace = sys.gettrace()

            def trace(frame, event, _arg):
                if frame.f_code is subject.run_short_screen.__code__:
                    frame.f_trace_opcodes = True
                    if (event == "opcode" and
                            frame.f_lasti == target_offset and
                            not injected[0]):
                        injected[0] = True
                        raise cleanup_interrupt
                return trace

            sys.settrace(trace)
            try:
                with self.assertRaises(subject.RunnerError) as raised:
                    self._run_timing_cleanup_failure(
                        managed, body_failure, restore_affinity)
            finally:
                sys.settrace(previous_trace)
        self.assertTrue(injected[0], target_opname)
        self.assertFalse(managed.stopped)
        self.assertTrue(managed.cleanup_complete)
        self.assertGreaterEqual(managed.resume_calls, 1)
        restore_affinity.assert_called()
        self.assertIn("injected timing body failure", str(raised.exception))
        self.assertIn(
            "injected first inner-finally opcode interrupt",
            str(raised.exception))

    def test_cleanup_error_preserves_authoritative_keyboard_interrupt(
            self) -> None:
        managed = _ManagedPauseFixture()
        restore_affinity = mock.Mock()
        primary = KeyboardInterrupt("authoritative timing interrupt")
        cleanup_error = subject.RunnerError(
            "injected cleanup-entry failure after resume")
        real_finish = subject._finish_timing_cleanup
        finish_calls = [0]

        def finish_then_fail(*args, **kwargs):
            finish_calls[0] += 1
            if finish_calls[0] == 1:
                managed.resume(args[-1])
                raise cleanup_error
            return real_finish(*args, **kwargs)

        with mock.patch.object(
                subject, "_finish_timing_cleanup",
                side_effect=finish_then_fail), \
                self.assertRaises(KeyboardInterrupt) as raised:
            self._run_timing_cleanup_failure(
                managed, primary, restore_affinity)
        self.assertIs(raised.exception, primary)
        self.assertFalse(managed.stopped)
        self.assertTrue(managed.cleanup_complete)
        self.assertGreaterEqual(finish_calls[0], 2)
        restore_affinity.assert_called()
        self.assertIsInstance(primary.__cause__, subject.RunnerError)
        self.assertIn(
            "injected cleanup-entry failure after resume",
            str(primary.__cause__))

    def test_interrupted_primary_capture_recovers_body_failure(self) -> None:
        managed = _ManagedPauseFixture()
        restore_affinity = mock.Mock()
        body_failure = subject.RunnerError(
            "injected body failure before capture")
        capture_interrupt = subject.RunnerError(
            "injected primary-capture interrupt")
        real_exc_info = sys.exc_info
        injected = [False]

        def interrupt_capture_once():
            if managed.stopped and not injected[0]:
                injected[0] = True
                raise capture_interrupt
            return real_exc_info()

        with mock.patch.object(
                subject.sys, "exc_info", side_effect=interrupt_capture_once), \
                self.assertRaises(subject.RunnerError) as raised:
            self._run_timing_cleanup_failure(
                managed, body_failure, restore_affinity)
        self.assertTrue(injected[0])
        self.assertFalse(managed.stopped)
        self.assertTrue(managed.cleanup_complete)
        restore_affinity.assert_called()
        self.assertIn(
            "injected body failure before capture", str(raised.exception))
        self.assertIn(
            "injected primary-capture interrupt", str(raised.exception))

    def test_cpu_parser_and_physical_core_selection(self) -> None:
        topology = self._topology({
            0: (0, 0, 0),
            1: (0, 1, 0), 2: (0, 2, 0),
            3: (0, 3, 1), 4: (0, 4, 1),
            5: (0, 5, 2), 6: (0, 6, 2),
            7: (0, 7, 3), 8: (0, 8, 3),
            9: (0, 9, 4),
            11: (0, 1, 0),
        })
        self.assertEqual(subject.parse_cpu_list("0,3,7"), [0, 3, 7])
        self.assertEqual(subject.parse_pid_list("1,3,7"), [1, 3, 7])
        for invalid in ("", "00", "0,0", "1,0", "0, 1", "+1"):
            with self.assertRaises(subject.RunnerError):
                subject.parse_cpu_list(invalid)
        for invalid in (
                "", "0", "01", "1,1", "2,1", "1, 2", "+1",
                str(subject.MAX_LINUX_PID + 1), "9" * 10000):
            with self.assertRaises(subject.RunnerError):
                subject.parse_pid_list(invalid)
        affinity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_controller=9, affinity=affinity,
            sysfs_root=topology), (list(range(1, 9)), 9))
        diverse_workers = [1, 2, 3, 4, 5, 6, 7, 9]
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_workers=diverse_workers, affinity=affinity,
            sysfs_root=topology), (diverse_workers, 8))
        secondary_smt = [2, 3, 4, 5, 6, 7, 8, 11]
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_workers=secondary_smt, explicit_controller=9,
            affinity=affinity, sysfs_root=topology), (secondary_smt, 9))
        self.assertEqual(subject.select_cpu_layout(
            0, affinity=affinity, sysfs_root=topology),
            ([1, 2, 3, 4, 5, 6, 7, 9], 8))
        invalid_rosters = (
            list(range(1, 8)),
            list(range(1, 10)),
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 2, 3, 4, 5, 6, 7, 11],
            [2, 3, 4, 5, 6, 7, 8, 10],
        )
        for invalid in invalid_rosters:
            with self.assertRaises(subject.RunnerError):
                subject.select_cpu_layout(
                    0, explicit_workers=invalid, affinity=affinity,
                    sysfs_root=topology)
        with self.assertRaises(subject.RunnerError):
            subject.select_cpu_layout(
                0, explicit_workers=list(range(1, 9)),
                explicit_controller=8, affinity=affinity,
                sysfs_root=topology)
        for invalid_affinity in ([True] + affinity, [-1] + affinity):
            with self.assertRaisesRegex(
                    subject.RunnerError, "affinity contains an invalid"):
                subject.select_cpu_layout(
                    0, affinity=invalid_affinity, sysfs_root=topology)

    def test_default_qualification_affinity_rejects_sampler_logical_or_core(
            self) -> None:
        contract = contract_api.load_contract()
        for sampler_cpu, topology in (
                (1, {0: (0, 0), 1: (0, 1)}),
                (17, {0: (0, 0), 1: (0, 17), 17: (0, 17)})):
            args = mock.Mock(
                deadline_seconds=1.0, contract=None, worker=Path("worker"),
                sampler_cpu=sampler_cpu)

            def cpu_topology(cpu, _root):
                return topology[cpu]

            with self.subTest(sampler_cpu=sampler_cpu), \
                    mock.patch.object(
                        contract_api, "load_contract", return_value=contract), \
                    mock.patch.object(
                        subject, "describe_worker", return_value={}), \
                    mock.patch.object(
                        os, "sched_getaffinity", return_value={0, 1}), \
                    mock.patch.object(
                        subject, "_cpu_topology", side_effect=cpu_topology), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "sampler physical core"):
                subject.run_short_screen(args)

    def test_llc_selection_maximizes_coverage_then_round_robin_fills(self) \
            -> None:
        topology = self._topology({
            0: (0, 0, 9),
            1: (0, 1, 0), 2: (0, 2, 0), 3: (0, 3, 0),
            4: (0, 4, 1), 5: (0, 5, 1),
            6: (0, 6, 2),
            7: (0, 7, 3), 8: (0, 8, 3), 9: (0, 9, 3),
            10: (0, 10, 4),
        })
        affinity = list(range(11))
        expected = ([1, 2, 3, 4, 5, 6, 7, 8], 10)
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_controller=10, affinity=affinity,
            sysfs_root=topology), expected)
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_controller=10, affinity=reversed(affinity),
            sysfs_root=topology), expected)

    def test_explicit_worker_roster_rejects_avoidable_llc_reuse(self) -> None:
        values = {0: (0, 0, 9), 10: (0, 10, 9)}
        values.update({1: (0, 1, 0), 2: (0, 2, 0)})
        for cpu in range(3, 10):
            values[cpu] = (0, cpu, cpu - 2)
        topology = self._topology(values)
        affinity = list(range(11))
        with self.assertRaisesRegex(
                subject.RunnerError, "maximize available LLC coverage"):
            subject.select_cpu_layout(
                0, explicit_workers=list(range(1, 9)),
                explicit_controller=10, affinity=affinity,
                sysfs_root=topology)
        diverse = [1, 3, 4, 5, 6, 7, 8, 9]
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_workers=diverse, explicit_controller=10,
            affinity=affinity, sysfs_root=topology), (diverse, 10))

    def test_llc_identity_includes_physical_package(self) -> None:
        values = {
            0: (0, 0, 8), 10: (0, 10, 9), 9: (0, 9, 0),
        }
        for cpu in range(1, 5):
            values[cpu] = (0, cpu, cpu - 1)
        for cpu in range(5, 9):
            values[cpu] = (1, cpu - 4, cpu - 5)
        topology = self._topology(values)
        affinity = list(range(11))
        under_diverse = [1, 2, 3, 4, 5, 6, 7, 9]
        with self.assertRaisesRegex(
                subject.RunnerError, "maximize available LLC coverage"):
            subject.select_cpu_layout(
                0, explicit_workers=under_diverse, explicit_controller=10,
                affinity=affinity, sysfs_root=topology)
        diverse = list(range(1, 9))
        self.assertEqual(subject.select_cpu_layout(
            0, explicit_workers=diverse, explicit_controller=10,
            affinity=affinity, sysfs_root=topology), (diverse, 10))

    def test_llc_identity_includes_cache_level(self) -> None:
        topology = self._topology({0: (0, 0, 0), 1: (0, 1, 0)})
        (topology / "cpu1" / "cache" / "index7" / "level").write_text(
            "4\n", encoding="ascii")
        self.assertEqual(subject._cpu_llc_identity(0, topology), (3, 0))
        self.assertEqual(subject._cpu_llc_identity(1, topology), (4, 0))

    def test_llc_topology_rejects_missing_and_malformed_metadata(self) -> None:
        values = {cpu: (0, cpu, cpu) for cpu in range(10)}
        topology = self._topology(values)
        affinity = list(range(10))
        cache_id = topology / "cpu1" / "cache" / "index7" / "id"
        cache_id.unlink()
        with self.assertRaisesRegex(subject.RunnerError, "cache ID"):
            subject.select_cpu_layout(
                0, explicit_controller=9, affinity=affinity,
                sysfs_root=topology)
        cache_id.write_text("1\n", encoding="ascii")
        level = topology / "cpu1" / "cache" / "index7" / "level"
        level.write_text("03\n", encoding="ascii")
        with self.assertRaisesRegex(subject.RunnerError, "malformed cache"):
            subject.select_cpu_layout(
                0, explicit_controller=9, affinity=affinity,
                sysfs_root=topology)

    def test_worker_description_is_exact_and_hashes_executable(self) -> None:
        worker = self._fake_worker()
        description = self._description(worker)
        self.assertEqual(
            description["binary_sha256"],
            hashlib.sha256(worker.read_bytes()).hexdigest())
        self.assertEqual(
            [arm["arm"] for arm in description["arms"]],
            [value[0] for value in subject.EXPECTED_ARMS])
        self.assertEqual(description["source_git_commit"], "1" * 40)
        subject._require_worker_source_commit(description, "1" * 40)
        with self.assertRaises(subject.RunnerError):
            subject._require_worker_source_commit(description, "2" * 40)

        candidate_digest_prefix = (
            "9527f200ad38c7eec6502b2f768fdd67")
        self.assertIn(candidate_digest_prefix, FAKE_WORKER)
        wrong_worker = self.root / "wrong_candidate_worker.py"
        wrong_worker.write_text(
            FAKE_WORKER.replace(candidate_digest_prefix, "0" * 32),
            encoding="utf-8")
        wrong_worker.chmod(0o755)
        with self.assertRaisesRegex(
                subject.RunnerError, "wrong timing candidate"):
            self._description(wrong_worker)

    def test_controller_affinity_is_singleton_and_restored(self) -> None:
        original = set(os.sched_getaffinity(0))
        if not original:
            self.skipTest("controller affinity is empty")
        selected = min(original)
        try:
            subject._pin_controller(selected)
            self.assertEqual(os.sched_getaffinity(0), {selected})
        finally:
            subject._restore_controller_affinity(original)
        self.assertEqual(os.sched_getaffinity(0), original)

    def test_sibling_guard_rejects_exact_affinity_occupant_and_attests_clean(
            self) -> None:
        values = {cpu: (0, cpu, cpu) for cpu in range(10)}
        values[11] = (0, 1, 1)
        topology = self._topology(values)
        workers = list(range(1, 9))
        foreign = [{
            "pid": 1234, "tid": 1234, "process_start_ticks": 99,
            "cpu_affinity": [11], "state": "R",
        }]
        guard = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: foreign,
            sysfs_root=topology)
        self.assertEqual(guard.sibling_cpus, [11])
        with self.assertRaisesRegex(
                subject.RunnerError, "foreign exact-affinity task"):
            guard.check("before-worker-spawn")

        controller_thread = [{
            "pid": os.getpid(), "tid": os.getpid() + 1,
            "process_start_ticks": 77,
            "cpu_affinity": [11], "state": "R",
        }]
        guard = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: controller_thread,
            sysfs_root=topology)
        with self.assertRaisesRegex(
                subject.RunnerError, "foreign exact-affinity task"):
            guard.check("before-worker-spawn")

        misplaced_controller = [{
            "pid": os.getpid(), "tid": os.getpid(),
            "process_start_ticks": 77,
            "cpu_affinity": [workers[0]], "state": "R",
        }]
        guard = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: misplaced_controller,
            sysfs_root=topology)
        with self.assertRaisesRegex(
                subject.RunnerError, "foreign exact-affinity task"):
            guard.check("before-timing-wave-0")

        pinned_controller = copy.deepcopy(misplaced_controller)
        pinned_controller[0]["cpu_affinity"] = [9]
        guard = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: pinned_controller,
            sysfs_root=topology)
        guard.check("before-timing-wave-0")

        worker_processes = []
        for cpu in workers:
            worker = mock.Mock()
            worker.cpu = cpu
            worker.process.pid = 3000 + cpu
            worker.start_ticks = 7000 + cpu
            worker_processes.append(worker)
        unexpected_worker_thread = [{
            "pid": worker_processes[0].process.pid,
            "tid": worker_processes[0].process.pid + 100,
            "process_start_ticks": worker_processes[0].start_ticks,
            "cpu_affinity": [worker_processes[0].cpu], "state": "R",
        }]
        guard = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: unexpected_worker_thread,
            sysfs_root=topology)
        guard.bind_workers(worker_processes)
        with self.assertRaisesRegex(
                subject.RunnerError, "foreign exact-affinity task"):
            guard.check("before-timing-wave-0")

        clean = subject.SiblingIsolationGuard(
            workers, 9, scan=lambda _cpus: [], sysfs_root=topology)
        clean.check("before-worker-spawn")
        first = clean._checks[0]["monotonic_ns"]
        clean.check("after-timing-interval")
        last = clean._checks[-1]["monotonic_ns"]
        attestation = clean.attestation()
        validated = native_api.validate_sibling_isolation(
            attestation, first, last, workers, 9, topology)
        self.assertEqual(validated, attestation)
        bool_cpu = copy.deepcopy(attestation)
        bool_cpu["worker_cpus"][0] = True
        with self.assertRaises(native_api.NativeEvidenceError):
            native_api.validate_sibling_isolation(
                bool_cpu, first, last, workers, 9, topology)
        bool_ordinal = copy.deepcopy(attestation)
        bool_ordinal["checks"][0]["ordinal"] = False
        bool_ordinal["checks_sha256"] = contract_api.sha256_json(
            bool_ordinal["checks"])
        with self.assertRaises(native_api.NativeEvidenceError):
            native_api.validate_sibling_isolation(
                bool_ordinal, first, last, workers, 9, topology)

    def test_exact_affinity_scan_fails_closed_on_pid_reuse(self) -> None:
        proc_root = self.root / "proc"
        process_root = proc_root / "4242"
        process_root.mkdir(parents=True)
        (process_root / "exe").write_bytes(b"executable")
        with mock.patch.object(
                subject, "_process_start_ticks",
                side_effect=(111, 222)), mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(4242, "R")]), mock.patch.object(
                        os, "sched_getaffinity", return_value={7}), \
                self.assertRaisesRegex(
                    subject.RunnerError,
                    "changed identity during SMT isolation scan"):
            subject._scan_exact_affinity_tasks({7}, proc_root)

        with mock.patch.object(
                subject, "_process_start_ticks",
                side_effect=(111, 111)), mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(4242, "R")]), mock.patch.object(
                        os, "sched_getaffinity", return_value={7}):
            self.assertEqual(
                subject._scan_exact_affinity_tasks({7}, proc_root),
                [{"pid": 4242, "tid": 4242,
                  "process_start_ticks": 111,
                  "cpu_affinity": [7], "state": "R"}])

    def test_exact_affinity_scan_retains_observed_process_on_exit(self) \
            -> None:
        proc_root = self.root / "proc-exit"
        process_root = proc_root / "4242"
        process_root.mkdir(parents=True)
        (process_root / "exe").write_bytes(b"executable")
        with mock.patch.object(
                subject, "_process_start_ticks",
                side_effect=(111, ProcessLookupError())), mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(4242, "R")]), mock.patch.object(
                        os, "sched_getaffinity", return_value={7}):
            self.assertEqual(
                subject._scan_exact_affinity_tasks({7}, proc_root),
                [{"pid": 4242, "tid": 4242,
                  "process_start_ticks": 111,
                  "cpu_affinity": [7], "state": "R"}])

    def test_managed_sibling_process_is_resumed_after_guard_window(self) -> None:
        if (not hasattr(os, "pidfd_open") or
                not hasattr(signal, "pidfd_send_signal")):
            self.skipTest("Linux pidfd signal support is required")
        allowed = sorted(os.sched_getaffinity(0))
        if not allowed:
            self.skipTest("one schedulable CPU is required")
        cpu = allowed[0]
        process = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
        manager = None
        try:
            os.sched_setaffinity(process.pid, {cpu})
            manager = subject.ManagedPausedSiblingProcesses()
            manager.pause([process.pid], [cpu], time.monotonic() + 5.0)
            states = subject._process_task_states(process.pid)
            self.assertTrue(states)
            self.assertTrue(all(state in ("T", "t")
                                for _, state in states))
            manager.resume(time.monotonic() + 5.0)
            states = subject._process_task_states(process.pid)
            self.assertTrue(states)
            self.assertTrue(all(state not in ("T", "t")
                                for _, state in states))
            self.assertIsNone(process.poll())
        finally:
            if manager is not None:
                manager.resume(time.monotonic() + 5.0)
            if process.poll() is None:
                process.terminate()
            process.wait(timeout=5.0)

    def test_managed_sibling_pause_rollback_resumes_earlier_process(self) \
            -> None:
        if (not hasattr(os, "pidfd_open") or
                not hasattr(signal, "pidfd_send_signal")):
            self.skipTest("Linux pidfd signal support is required")
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < 2:
            self.skipTest("two schedulable CPUs are required")
        processes = [
            subprocess.Popen(
                [sys.executable, "-c", "import time; time.sleep(60)"],
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            for _ in range(2)
        ]
        try:
            os.sched_setaffinity(processes[0].pid, {allowed[0]})
            os.sched_setaffinity(processes[1].pid, {allowed[1]})
            manager = subject.ManagedPausedSiblingProcesses()
            with self.assertRaisesRegex(
                    subject.RunnerError, "not pinned only to T siblings"):
                manager.pause(
                    [process.pid for process in processes], [allowed[0]],
                    time.monotonic() + 5.0)
            states = subject._process_task_states(processes[0].pid)
            self.assertTrue(states)
            self.assertTrue(all(state not in ("T", "t")
                                for _, state in states))
            self.assertIsNone(processes[0].poll())
        finally:
            for process in processes:
                try:
                    os.kill(process.pid, signal.SIGCONT)
                except ProcessLookupError:
                    pass
                if process.poll() is None:
                    process.terminate()
                process.wait(timeout=5.0)

    def test_stop_deadline_rollback_uses_fresh_resume_grace(self) -> None:
        if (not hasattr(os, "pidfd_open") or
                not hasattr(signal, "pidfd_send_signal")):
            self.skipTest("Linux pidfd signal support is required")
        allowed = sorted(os.sched_getaffinity(0))
        if not allowed:
            self.skipTest("one schedulable CPU is required")
        cpu = allowed[0]
        process = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
        original_states = subject._process_task_states
        state_calls = 0

        def staged_states(pid: int, *args, **kwargs):
            nonlocal state_calls
            state_calls += 1
            if state_calls == 2:
                return [(pid, "R")]
            if state_calls == 3:
                return [(pid, "T")]
            return original_states(pid, *args, **kwargs)

        try:
            os.sched_setaffinity(process.pid, {cpu})
            manager = subject.ManagedPausedSiblingProcesses()
            with mock.patch.object(
                    subject, "_process_task_states",
                    side_effect=staged_states), self.assertRaisesRegex(
                        subject.RunnerError,
                        "hard wall expired during stopping managed") as raised:
                manager.pause(
                    [process.pid], [cpu], time.monotonic() - 1.0)
            self.assertNotIn("rollback failed", str(raised.exception))
            states = original_states(process.pid)
            self.assertTrue(states)
            self.assertTrue(all(state not in ("T", "t")
                                for _, state in states))
        finally:
            try:
                os.kill(process.pid, signal.SIGCONT)
            except ProcessLookupError:
                pass
            if process.poll() is None:
                process.terminate()
            process.wait(timeout=5.0)

    def test_pid_reuse_rollback_resumes_exact_signaled_process(self) -> None:
        if (not hasattr(os, "pidfd_open") or
                not hasattr(signal, "pidfd_send_signal")):
            self.skipTest("Linux pidfd signal support is required")
        allowed = sorted(os.sched_getaffinity(0))
        if not allowed:
            self.skipTest("one schedulable CPU is required")
        cpu = allowed[0]
        process = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
        real_start_ticks = subject._process_start_ticks
        real_pidfd_signal = signal.pidfd_send_signal
        start_calls = 0
        observed_signals = []

        def changed_identity(pid: int, *args, **kwargs):
            nonlocal start_calls
            start_calls += 1
            value = real_start_ticks(pid, *args, **kwargs)
            return value + 1 if start_calls == 2 else value

        def record_signal(pidfd, sig, siginfo=None, flags=0):
            observed_signals.append(sig)
            return real_pidfd_signal(pidfd, sig, siginfo, flags)

        try:
            os.sched_setaffinity(process.pid, {cpu})
            manager = subject.ManagedPausedSiblingProcesses()
            with mock.patch.object(
                    subject, "_process_start_ticks",
                    side_effect=changed_identity), mock.patch.object(
                        signal, "pidfd_send_signal",
                        side_effect=record_signal), self.assertRaisesRegex(
                            subject.RunnerError,
                            "changed identity while stopping"):
                manager.pause(
                    [process.pid], [cpu], time.monotonic() + 5.0)
            self.assertIn(signal.SIGSTOP, observed_signals)
            self.assertIn(signal.SIGCONT, observed_signals)
            states = subject._process_task_states(process.pid)
            self.assertTrue(states)
            self.assertTrue(all(state not in ("T", "t")
                                for _, state in states))
        finally:
            try:
                os.kill(process.pid, signal.SIGCONT)
            except ProcessLookupError:
                pass
            if process.poll() is None:
                process.terminate()
            process.wait(timeout=5.0)

    def test_partial_resume_failure_retains_every_exact_handle_for_retry(
            self) -> None:
        records = [
            {"pid": 101, "process_start_ticks": 1001,
             "cpu_affinity": [11]},
            {"pid": 202, "process_start_ticks": 2002,
             "cpu_affinity": [12]},
        ]
        manager = subject.ManagedPausedSiblingProcesses(records)
        first = FileIO(os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        second = FileIO(os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(first.close)
        self.addCleanup(second.close)
        first_fd = first.fileno()
        second_fd = second.fileno()
        manager._pidfds = {101: first, 202: second}
        failed_once = False
        signals = []

        def signal_with_one_failure(pidfd, sig, siginfo=None, flags=0):
            nonlocal failed_once
            signals.append((pidfd, sig))
            if pidfd == second_fd and not failed_once:
                failed_once = True
                raise OSError(errno.EIO, "injected resume failure")

        def start_ticks(pid, *args, **kwargs):
            return 1001 if pid == 101 else 2002

        with mock.patch.object(
                signal, "pidfd_send_signal",
                side_effect=signal_with_one_failure,
                create=True), mock.patch.object(
                    subject, "_process_start_ticks",
                    side_effect=start_ticks), mock.patch.object(
                    subject, "_process_task_states",
                        return_value=[(1, "R")]):
            manager.resume(time.monotonic() + 5.0)
            self.assertTrue(manager._resumed)
            self.assertEqual(manager._pidfds, {})
            self.assertEqual(
                signals,
                [(first_fd, signal.SIGCONT),
                 (second_fd, signal.SIGCONT),
                 (second_fd, signal.SIGCONT)])
            self.assertTrue(first.closed)
            self.assertTrue(second.closed)

    def test_resume_defers_control_flow_until_every_pid_is_continued(
            self) -> None:
        class ResumeInterrupt(BaseException):
            pass

        records = [
            {"pid": 101, "process_start_ticks": 1001,
             "cpu_affinity": [11]},
            {"pid": 202, "process_start_ticks": 2002,
             "cpu_affinity": [12]},
        ]
        manager = subject.ManagedPausedSiblingProcesses(records)
        first = FileIO(os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        second = FileIO(os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(first.close)
        self.addCleanup(second.close)
        first_fd = first.fileno()
        second_fd = second.fileno()
        manager._pidfds = {101: first, 202: second}
        control = ResumeInterrupt("injected one-shot resume interrupt")
        signals = []

        def signal_with_interrupt(pidfd, sig, siginfo=None, flags=0):
            signals.append((pidfd, sig))
            if (pidfd == second_fd and
                    signals.count((second_fd, sig)) == 1):
                raise control

        def start_ticks(pid, *args, **kwargs):
            return 1001 if pid == 101 else 2002

        with mock.patch.object(
                signal, "pidfd_send_signal",
                side_effect=signal_with_interrupt,
                create=True), mock.patch.object(
                    subject, "_process_start_ticks",
                    side_effect=start_ticks), mock.patch.object(
                    subject, "_process_task_states",
                        return_value=[(1, "R")]), self.assertRaises(
                                ResumeInterrupt) as raised:
            manager.resume(time.monotonic() + 5.0)
        self.assertIs(raised.exception, control)
        self.assertTrue(manager._resumed)
        self.assertEqual(manager._pidfds, {})
        self.assertEqual(
            signals,
            [(first_fd, signal.SIGCONT),
             (second_fd, signal.SIGCONT),
             (second_fd, signal.SIGCONT)])
        self.assertTrue(first.closed)
        self.assertTrue(second.closed)
        manager.resume(time.monotonic() + 5.0)

    def test_pidfd_close_interrupt_leaves_resume_state_converged(self) -> None:
        class CloseInterrupt(BaseException):
            pass

        records = [{
            "pid": 101, "process_start_ticks": 1001,
            "cpu_affinity": [11],
        }, {
            "pid": 202, "process_start_ticks": 2002,
            "cpu_affinity": [12],
        }]
        manager = subject.ManagedPausedSiblingProcesses(records)
        control = CloseInterrupt("injected pidfd close interrupt")

        class InterruptAfterClose(FileIO):
            armed = True

            def close(self):
                super().close()
                if self.armed:
                    self.armed = False
                    raise control

        first = InterruptAfterClose(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        second = FileIO(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(first.close)
        self.addCleanup(second.close)
        manager._pidfds = {101: first, 202: second}

        with mock.patch.object(
                signal, "pidfd_send_signal", create=True), mock.patch.object(
                    subject, "_process_start_ticks",
                    side_effect=lambda pid: 1001 if pid == 101 else 2002), \
                mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(1, "R")]), \
                self.assertRaises(CloseInterrupt) as raised:
            manager.resume(time.monotonic() + 5.0)
        self.assertIs(raised.exception, control)
        self.assertTrue(first.closed)
        self.assertTrue(second.closed)
        self.assertTrue(manager._resumed)
        self.assertEqual(manager._pidfds, {})
        manager.resume(time.monotonic() + 5.0)

    def test_resume_state_store_interrupt_keeps_exact_close_ownership(
            self) -> None:
        class StateInterrupt(BaseException):
            pass

        control = StateInterrupt("injected resume-state store interrupt")

        class InterruptingManager(
                subject.ManagedPausedSiblingProcesses):
            armed = False

            def __setattr__(self, name, value):
                if name == "_resumed" and value is True and self.armed:
                    object.__setattr__(self, "armed", False)
                    raise control
                super().__setattr__(name, value)

        manager = InterruptingManager([{
            "pid": 101, "process_start_ticks": 1001,
            "cpu_affinity": [11],
        }])
        handle = FileIO(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(handle.close)
        manager._pidfds = {101: handle}
        manager.armed = True
        with mock.patch.object(
                signal, "pidfd_send_signal", create=True), mock.patch.object(
                    subject, "_process_start_ticks", return_value=1001), \
                mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(1, "R")]), self.assertRaises(
                            StateInterrupt) as raised:
            manager.resume(time.monotonic() + 5.0)
        self.assertIs(raised.exception, control)
        self.assertTrue(manager._resumed)
        self.assertEqual(manager._pidfds, {})
        self.assertTrue(handle.closed)
        manager.resume(time.monotonic() + 5.0)

    def test_pre_drain_interrupt_retries_persistent_close_roster(self) -> None:
        class DrainInterrupt(BaseException):
            pass

        control = DrainInterrupt("injected pre-drain interrupt")

        class InterruptingManager(
                subject.ManagedPausedSiblingProcesses):
            armed = True

            def _drain_resumed_pidfds(self):
                if self.armed:
                    self.armed = False
                    raise control
                return super()._drain_resumed_pidfds()

        manager = InterruptingManager([{
            "pid": 101, "process_start_ticks": 1001,
            "cpu_affinity": [11],
        }])
        handle = FileIO(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(handle.close)
        manager._pidfds = {101: handle}
        with mock.patch.object(
                signal, "pidfd_send_signal", create=True), mock.patch.object(
                    subject, "_process_start_ticks", return_value=1001), \
                mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(1, "R")]), self.assertRaises(
                            DrainInterrupt) as raised:
            manager.resume(time.monotonic() + 5.0)
        self.assertIs(raised.exception, control)
        self.assertTrue(manager._resumed)
        self.assertEqual(manager._pidfds, {})
        self.assertTrue(handle.closed)
        manager.resume(time.monotonic() + 5.0)

    def test_post_close_interrupt_never_retries_reused_descriptor(self) \
            -> None:
        class PostCloseInterrupt(BaseException):
            pass

        control = PostCloseInterrupt("injected interrupt after real close")
        descriptor = os.open(os.devnull, os.O_RDONLY)

        class InterruptAfterClose(FileIO):
            armed = True

            def close(self):
                super().close()
                if self.armed:
                    self.armed = False
                    raise control

        handle = InterruptAfterClose(descriptor, "rb", closefd=True)
        self.addCleanup(lambda: FileIO.close(handle))
        manager = subject.ManagedPausedSiblingProcesses([{
            "pid": 101, "process_start_ticks": 1001,
            "cpu_affinity": [11],
        }])
        manager._pidfds = {101: handle}

        replacement = -1
        try:
            with mock.patch.object(
                    signal, "pidfd_send_signal", create=True), \
                    mock.patch.object(
                        subject, "_process_start_ticks",
                        return_value=1001), mock.patch.object(
                            subject, "_process_task_states",
                            return_value=[(1, "R")]), \
                    self.assertRaises(PostCloseInterrupt) as raised:
                manager.resume(time.monotonic() + 5.0)
            self.assertIs(raised.exception, control)
            self.assertEqual(manager._pidfds, {})
            self.assertTrue(handle.closed)
            replacement = os.open(os.devnull, os.O_RDONLY)
            self.assertEqual(replacement, descriptor)
            manager.resume(time.monotonic() + 5.0)
            os.fstat(replacement)
        finally:
            if replacement >= 0:
                os.close(replacement)

    def test_pre_close_interrupt_retains_exact_pidfd_owner(self) -> None:
        class PreCloseInterrupt(BaseException):
            pass

        control = PreCloseInterrupt("injected interrupt before FileIO close")

        class InterruptBeforeClose(FileIO):
            armed = True

            def close(self):
                if self.armed:
                    self.armed = False
                    raise control
                return super().close()

        handle = InterruptBeforeClose(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(lambda: FileIO.close(handle))
        manager = subject.ManagedPausedSiblingProcesses()
        manager._resumed = True
        manager._pidfds = {101: handle}
        with self.assertRaises(PreCloseInterrupt) as raised:
            manager._drain_resumed_pidfds()
        self.assertIs(raised.exception, control)
        self.assertIs(manager._pidfds.get(101), handle)
        self.assertFalse(handle.closed)
        self.assertFalse(manager.cleanup_complete)
        manager._drain_resumed_pidfds()
        self.assertTrue(handle.closed)
        self.assertTrue(manager.cleanup_complete)

    def test_pidfd_drain_identity_pop_preserves_swapped_owner(self) -> None:
        manager = subject.ManagedPausedSiblingProcesses()
        manager._resumed = True
        replacement = FileIO(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(replacement.close)

        class SwapOwnerAfterClose(FileIO):
            armed = True

            def close(self):
                super().close()
                if self.armed:
                    self.armed = False
                    manager._pidfds[101] = replacement

        original = SwapOwnerAfterClose(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(lambda: FileIO.close(original))
        manager._pidfds = {101: original}
        manager._drain_resumed_pidfds()
        self.assertTrue(original.closed)
        self.assertIs(manager._pidfds.get(101), replacement)
        self.assertFalse(replacement.closed)
        manager._drain_resumed_pidfds()
        self.assertTrue(replacement.closed)
        self.assertTrue(manager.cleanup_complete)

    def test_pidfd_drain_persistent_failure_is_bounded(self) -> None:
        class PersistentDrainFailure(subject.ManagedPausedSiblingProcesses):
            drain_calls = 0

            def _drain_resumed_pidfds(self):
                self.drain_calls += 1
                raise subject.RunnerError("injected persistent drain failure")

        for future_deadline, expected_calls in ((False, 2), (True, 8)):
            with self.subTest(future_deadline=future_deadline):
                manager = PersistentDrainFailure()
                manager._resumed = True
                handle = FileIO(
                    os.open(os.devnull, os.O_RDONLY),
                    "rb", closefd=True)
                self.addCleanup(handle.close)
                manager._pidfds = {101: handle}
                deadline = time.monotonic() + (60.0 if future_deadline else -1.0)
                with self.assertRaisesRegex(
                        subject.RunnerError, "did not converge"):
                    manager.resume(deadline)
                self.assertEqual(manager.drain_calls, expected_calls)
                self.assertIs(manager._pidfds.get(101), handle)
                self.assertFalse(handle.closed)
                self.assertFalse(manager.cleanup_complete)
                handle.close()
                manager._pidfds.clear()

    def test_timing_cleanup_gives_paused_load_a_fresh_grace_period(self) \
            -> None:
        manager = mock.Mock()
        manager.cleanup_complete = False

        def resume(_deadline):
            manager.cleanup_complete = True

        manager.resume.side_effect = resume
        subject._finish_timing_cleanup(
            [], True, False, set(), None, manager, False,
            time.monotonic() - 60.0)
        self.assertGreater(
            manager.resume.call_args.args[0], time.monotonic() + 4.0)

    def test_timing_cleanup_retries_pre_body_resume_failure(self) -> None:
        primary = subject.RunnerError("injected campaign failure")
        interruption = subject.RunnerError(
            "injected pre-body resume interruption")

        class InterruptedManager(
                subject.ManagedPausedSiblingProcesses):
            resume_calls = 0

            def resume(self, deadline):
                self.resume_calls += 1
                if self.resume_calls == 1:
                    raise interruption
                return super().resume(deadline)

        manager = InterruptedManager([{
            "pid": 101, "process_start_ticks": 1001,
            "cpu_affinity": [11],
        }])
        handle = FileIO(
            os.open(os.devnull, os.O_RDONLY), "rb", closefd=True)
        self.addCleanup(handle.close)
        descriptor = handle.fileno()
        manager._pidfds = {101: handle}
        signals = []

        def send_signal(pidfd, sig, siginfo=None, flags=0):
            signals.append((pidfd, sig))

        with mock.patch.object(
                signal, "pidfd_send_signal", side_effect=send_signal,
                create=True), mock.patch.object(
                    subject, "_process_start_ticks", return_value=1001), \
                mock.patch.object(
                    subject, "_process_task_states",
                    return_value=[(1, "R")]), self.assertRaises(
                            subject.RunnerError) as raised:
            subject._finish_timing_cleanup(
                [], True, False, set(), primary, manager, False,
                time.monotonic() - 60.0)
        self.assertIs(raised.exception.__cause__, primary)
        self.assertIn("pre-body resume interruption", str(raised.exception))
        self.assertEqual(manager.resume_calls, 2)
        self.assertEqual(signals, [(descriptor, signal.SIGCONT)])
        self.assertTrue(manager.cleanup_complete)
        self.assertEqual(manager._pidfds, {})
        self.assertTrue(handle.closed)

    def test_live_exact_affinity_sibling_occupant_fails_before_timing(self) \
            -> None:
        if os.geteuid() != 0:
            self.skipTest("root is required for a fail-closed all-/proc scan")
        allowed = sorted(os.sched_getaffinity(0))
        by_core: Dict[tuple, List[int]] = {}
        for cpu in allowed:
            by_core.setdefault(
                subject._cpu_topology(cpu, subject.CPU_SYSFS_ROOT), []).append(cpu)
        pair = next((sorted(cpus) for cpus in by_core.values()
                     if len(cpus) >= 2), None)
        if pair is None:
            self.skipTest("an allowed SMT sibling pair is required")
        worker_cpu, occupied_sibling = pair[:2]
        controller_cpu = next((
            cpu for cpu in allowed
            if subject._cpu_topology(cpu, subject.CPU_SYSFS_ROOT) !=
                subject._cpu_topology(worker_cpu, subject.CPU_SYSFS_ROOT)), None)
        if controller_cpu is None:
            self.skipTest("a separate controller core is required")
        process = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)
        original_affinity = set(os.sched_getaffinity(0))
        try:
            os.sched_setaffinity(process.pid, {occupied_sibling})
            subject._pin_controller(controller_cpu)
            guard = subject.SiblingIsolationGuard(
                [worker_cpu], controller_cpu)
            self.assertIn(occupied_sibling, guard.sibling_cpus)
            with self.assertRaisesRegex(
                    subject.RunnerError,
                    "pid={}.*cpus={}".format(
                        process.pid, occupied_sibling)):
                guard.check("before-worker-spawn")
        finally:
            subject._restore_controller_affinity(original_affinity)
            if process.poll() is None:
                process.terminate()
            process.wait(timeout=5.0)

    def test_timing_freeze_is_published_before_result_sink(self) -> None:
        worker = self._fake_worker()
        description = self._description(worker)
        contract = contract_api.load_contract()
        output = self.root / "output"
        output.mkdir()
        trace_hash = hashlib.sha256(b"timing trace").hexdigest()
        qualification = self._qualification(contract, description=description)
        freeze = subject.write_development_timing_freeze(
            contract, description, [0], 1, "1" * 40, trace_hash, output,
            qualification)
        path = output / "timing-freeze.json"
        self.assertTrue(path.is_file())
        self.assertFalse((output / "recovery-freeze.json").exists())
        self.assertEqual(freeze["trace_manifest_sha256"], trace_hash)
        self.assertEqual(freeze["arm_roster"], [
            "wirehair2_head", "wirehair1",
            "wirehair2_dense_two07_basis_v1"])
        self.assertEqual(freeze["host_identity"]["controller_cpu"], 1)
        self.assertEqual(
            [value["construction_policy"] for value in freeze["arms"]],
            ["raw_base", "not_applicable", "raw_base"])
        sink = subject.AtomicLineSink(output / "later-native-results.jsonl")
        try:
            self.assertTrue(path.exists())
        finally:
            sink.abort()

    def test_completed_timing_loader_reopens_exact_terminal_bundle(self) \
            -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture()
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated) as validate, mock.patch.object(
                contract_api, "load_freeze_manifest",
                return_value=freeze):
            loaded = subject.load_completed_timing_screen(
                contract, directory)
        self.assertEqual(set(loaded), {
            "directory", "directory_identity", "run_summary", "freeze",
            "summary", "execution_receipt", "timing_qualification",
        })
        self.assertEqual(loaded["run_summary"], run_summary)
        self.assertEqual(loaded["freeze"], freeze)
        self.assertEqual(loaded["summary"], validated["summary"])
        self.assertEqual(
            loaded["execution_receipt"], validated["execution_receipt"])
        self.assertIs(loaded["timing_qualification"], qualification)
        self.assertFalse(validate.call_args.kwargs["verify_live_sampler"])
        self.assertEqual(validate.call_args.args[1:3], (
            "timing", "development"))

    def test_completed_timing_loader_rejects_provenance_substitution(self) \
            -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture()

        def publish(value):
            unsigned = {
                key: item for key, item in value.items()
                if key != "summary_sha256"
            }
            value["summary_sha256"] = contract_api.sha256_json(unsigned)
            (directory / "run-summary.json").write_bytes(
                (contract_api.canonical_json(value) + "\n").encode("utf-8"))

        mutations = (
            ("result hash", "timing_result_sha256", "1" * 64),
            ("execution hash", "timing_execution_receipt_sha256", "2" * 64),
            ("qualification execution hash",
             "timing_qualification_execution_receipt_sha256", "3" * 64),
            ("validator hash", "timing_validator_summary_sha256", "4" * 64),
            ("source commit", "source_git_commit", "2" * 40),
        )
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest",
                return_value=freeze):
            for name, field, replacement in mutations:
                changed = copy.deepcopy(run_summary)
                changed[field] = replacement
                publish(changed)
                with self.subTest(name=name), self.assertRaises(
                        subject.RunnerError):
                    subject.load_completed_timing_screen(contract, directory)
        publish(copy.deepcopy(run_summary))

    def test_completed_timing_loader_fails_closed_on_host_and_symlink(self) \
            -> None:
        (contract, qualification, directory, freeze, _run_summary,
         validated) = self._completed_timing_fixture()
        invalid_freeze = copy.deepcopy(freeze)
        invalid_freeze["host_identity"] = {"name": "missing controller"}
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest",
                return_value=invalid_freeze), self.assertRaisesRegex(
                    subject.RunnerError, "valid controller CPU"):
            subject.load_completed_timing_screen(contract, directory)

        d12_freeze = copy.deepcopy(freeze)
        d12_freeze["arm_roster"].insert(
            2, "wirehair2_raw_d12_h12_periodic")
        d12_freeze["arms"].insert(2, {
            "arm": "wirehair2_raw_d12_h12_periodic",
            "codec": "wirehair2_experiment",
            "binary_sha256": "a" * 64,
            "arm_descriptor_sha256": "9" * 64,
        })
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest",
                return_value=d12_freeze), self.assertRaisesRegex(
                    subject.RunnerError, "wrong exact arm roster"):
            subject.load_completed_timing_screen(contract, directory)

        target = directory / "timing-results-target.jsonl"
        target.write_bytes(b"{}\n")
        result = directory / "timing-results.jsonl"
        result.unlink()
        result.symlink_to(target.name)
        with self.assertRaisesRegex(
                subject.RunnerError, "cannot read completed timing result"):
            subject.load_completed_timing_screen(contract, directory)

    def test_completed_read_rejects_named_entry_replacement(self) -> None:
        directory = self.root / "entry-replacement"
        directory.mkdir()
        artifact = directory / "artifact.json"
        replacement = directory / "replacement.json"
        artifact.write_bytes(b"old-bytes\n")
        replacement.write_bytes(b"new-bytes\n")
        directory_fd = os.open(
            str(directory), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        real_open = os.open
        replaced = False

        def replace_after_open(path, *args, **kwargs):
            nonlocal replaced
            descriptor = real_open(path, *args, **kwargs)
            if (path == artifact.name and kwargs.get("dir_fd") == directory_fd
                    and not replaced):
                os.replace(str(replacement), str(artifact))
                replaced = True
            return descriptor

        try:
            with mock.patch.object(
                    subject.os, "open", side_effect=replace_after_open), \
                    self.assertRaisesRegex(
                        subject.RunnerError,
                        "single-link|changed while it was being read"):
                subject._read_completed_regular_bytes(
                    Path(artifact.name), "replacement artifact", directory_fd,
                    1024)
        finally:
            os.close(directory_fd)
        self.assertTrue(replaced)

    def test_completed_loader_rejects_hardlinked_dependency(self) -> None:
        (contract, _qualification, directory, _freeze, _run_summary,
         _validated) = self._completed_timing_fixture("hardlinked-dependency")
        result = directory / "timing-results.jsonl"
        os.link(str(result), str(directory / "result-hardlink.jsonl"))
        with self.assertRaisesRegex(
                subject.RunnerError, "regular non-symlink file"):
            subject.load_completed_timing_screen(contract, directory)

    def test_completed_loader_rejects_mutation_during_semantic_validation(
            self) -> None:
        (contract, qualification, directory, freeze, _run_summary,
         validated) = self._completed_timing_fixture("semantic-mutation")
        replacement = self.root / "replacement-results.jsonl"
        replacement.write_bytes(b'{"mutated":true}\n')

        def mutate(*_args, **_kwargs):
            os.replace(
                str(replacement), str(directory / "timing-results.jsonl"))
            return validated

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=mutate), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaisesRegex(
                    subject.RunnerError, "changed during semantic validation"):
            subject.load_completed_timing_screen(contract, directory)

    def test_completed_loader_rejects_parent_replacement_during_validation(
            self) -> None:
        (contract, qualification, directory, freeze, _run_summary,
         validated) = self._completed_timing_fixture("parent-replacement")
        moved = self.root / "moved-parent-replacement"

        def replace_parent(*_args, **_kwargs):
            directory.rename(moved)
            directory.mkdir()
            return validated

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=replace_parent), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaisesRegex(
                    subject.RunnerError, "directory changed before terminal"):
            subject.load_completed_timing_screen(contract, directory)

    def test_completed_bounds_cover_full_qualification_retry_envelope(
            self) -> None:
        self.assertEqual(
            subject.MAX_QUALIFICATION_ATTEMPTS, 2304 * 256)
        self.assertEqual(
            subject.COMPLETED_ARTIFACT_BYTE_LIMITS[
                "timing-qualification-native.jsonl"],
            2304 * 256 * subject.MAX_QUALIFICATION_NATIVE_RECORD_BYTES)
        self.assertGreater(
            subject.COMPLETED_ARTIFACT_BYTE_LIMITS[
                "timing-qualification-native.jsonl"],
            int(540.6 * 1024 * 1024))
        self.assertEqual(
            subject.MAX_COMPLETED_BUNDLE_BYTES,
            sum(subject.COMPLETED_ARTIFACT_BYTE_LIMITS.values()))

    def test_producer_applies_loader_cap_to_every_dependency(self) -> None:
        (_contract, _qualification, directory, _freeze, _run_summary,
         _validated) = self._completed_timing_fixture("producer-bounds")
        (directory / "run-summary.json").unlink()
        parent_fd, _parent_identity = subject._open_completion_parent(
            directory / "run-summary.json")
        try:
            for name in subject.COMPLETED_DEPENDENCY_NAMES.values():
                descriptors = {}
                with self.subTest(name=name), mock.patch.dict(
                        subject.COMPLETED_ARTIFACT_BYTE_LIMITS,
                        {name: 1}), self.assertRaisesRegex(
                            subject.RunnerError, "bounded artifact size"):
                    subject._open_pinned_completed_bundle(
                        subject.COMPLETED_DEPENDENCY_NAMES, parent_fd,
                        descriptors)
                self.assertEqual(descriptors, {})
            holder = [-1]
            with mock.patch.dict(
                    subject.COMPLETED_ARTIFACT_BYTE_LIMITS,
                    {"run-summary.json": 1}), self.assertRaisesRegex(
                        subject.RunnerError, "publication byte cap"):
                subject._open_unnamed_completion_marker(
                    parent_fd, b"{}\n", holder)
            self.assertEqual(holder, [-1])
        finally:
            os.close(parent_fd)

    def test_validator_failure_leaves_no_marker_or_staging_name(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         _validated) = self._completed_timing_fixture("validator-failure")
        marker = directory / "run-summary.json"
        marker.unlink()
        names_before = {path.name for path in directory.iterdir()}

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=subject.RunnerError(
                    "injected semantic validator failure")), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaisesRegex(
                    subject.RunnerError, "semantic validator failure"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze, {}, {},
                qualification)
        self.assertFalse(marker.exists())
        self.assertEqual(
            {path.name for path in directory.iterdir()}, names_before)

    def test_strict_loader_cannot_observe_marker_during_validation(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("private-validation")
        marker = directory / "run-summary.json"
        marker.unlink()
        observations = []

        def validate_while_private(*_args, **_kwargs):
            visible = marker.exists()
            observations.append(visible)
            if not visible:
                with self.assertRaisesRegex(
                        subject.RunnerError,
                        "cannot read completed timing summary"):
                    subject.load_completed_timing_screen(contract, directory)
            return validated

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=validate_while_private), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(observations, [False, True])
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_hard_wall_failure_leaves_no_marker(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("hard-wall-failure")
        marker = directory / "run-summary.json"
        marker.unlink()

        def fail_finish() -> None:
            raise OSError(errno.EIO, "injected hard-wall finish failure")

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaisesRegex(OSError, "hard-wall finish failure"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification, fail_finish)
        self.assertFalse(marker.exists())

    def test_dependency_mutation_during_validation_fails_prelink(self) -> None:
        for index, name in enumerate((
                "timing-freeze.json", "timing-results.jsonl")):
            (contract, qualification, directory, freeze, run_summary,
             validated) = self._completed_timing_fixture(
                 "dependency-mutation-{}".format(index))
            marker = directory / "run-summary.json"
            marker.unlink()
            replacement = self.root / "replacement-{}".format(index)
            replacement.write_bytes(b'{"mutated":true}\n')
            mutated = False

            def mutate_dependency(*_args, _name=name, **_kwargs):
                nonlocal mutated
                if not mutated:
                    os.replace(str(replacement), str(directory / _name))
                    mutated = True
                return validated

            with self.subTest(name=name), mock.patch.object(
                    contract_api, "load_timing_qualification_map",
                    return_value=qualification), mock.patch.object(
                    native_api, "validate_execution_receipt",
                    side_effect=mutate_dependency), mock.patch.object(
                    contract_api, "load_freeze_manifest",
                    return_value=freeze), self.assertRaisesRegex(
                        subject.RunnerError,
                        "single-link|changed during semantic validation"):
                subject._commit_completed_timing_screen(
                    contract, directory, run_summary, freeze,
                    validated["summary"], validated["execution_receipt"],
                    qualification)
            self.assertTrue(mutated)
            self.assertFalse(marker.exists())

    def test_finish_mutation_is_caught_by_post_finish_guard(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("finish-mutation")
        marker = directory / "run-summary.json"
        marker.unlink()
        result = directory / "timing-results.jsonl"
        replacement = self.root / "finish-replacement.jsonl"
        replacement.write_bytes(b'{"mutated_after_reread":true}\n')

        def mutate_at_finish() -> None:
            os.replace(str(replacement), str(result))

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaisesRegex(
                    subject.RunnerError, "changed before publication"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification, mutate_at_finish)
        self.assertFalse(marker.exists())

    def test_parent_replacement_during_validation_fails_prelink(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("parent-replacement-commit")
        marker = directory / "run-summary.json"
        marker.unlink()
        moved = self.root / "moved-parent-replacement-commit"

        def replace_parent(*_args, **_kwargs):
            directory.rename(moved)
            directory.mkdir()
            return validated

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=replace_parent), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                self.assertRaises(subject.RunnerError):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertFalse(marker.exists())
        self.assertFalse((moved / "run-summary.json").exists())

    def test_unsupported_otmpfile_leaves_no_name(self) -> None:
        directory = self.root / "unsupported-otmpfile"
        directory.mkdir()
        marker = directory / "run-summary.json"
        parent_fd, _identity = subject._open_completion_parent(marker)
        holder = [-1]
        try:
            with mock.patch.object(subject.os, "O_TMPFILE", 0), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "requires O_TMPFILE"):
                subject._open_unnamed_completion_marker(
                    parent_fd, b"{}\n", holder)
        finally:
            os.close(parent_fd)
        self.assertEqual(holder, [-1])
        self.assertEqual(list(directory.iterdir()), [])

    def test_unsupported_proc_link_leaves_no_name(self) -> None:
        directory = self.root / "unsupported-proc-link"
        directory.mkdir()
        marker = directory / "run-summary.json"
        parent_fd, parent_identity = subject._open_completion_parent(marker)
        holder = [-1]
        try:
            subject._open_unnamed_completion_marker(
                parent_fd, b"{}\n", holder)
            with mock.patch.object(
                    subject.os, "link",
                    side_effect=OSError(errno.EOPNOTSUPP, "no proc link")), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "cannot publish"):
                subject._link_unnamed_completion_marker(
                    marker, holder[0], parent_fd, parent_identity)
        finally:
            if holder[0] >= 0:
                os.close(holder[0])
            os.close(parent_fd)
        self.assertEqual(list(directory.iterdir()), [])

    def test_unnamed_stage_return_interrupt_is_closed_without_a_name(
            self) -> None:
        directory = self.root / "unnamed-return-interrupt"
        directory.mkdir()
        marker = directory / "run-summary.json"
        helper_code = subject._open_unnamed_completion_marker.__code__
        captured = []

        class ReturnInterrupt(BaseException):
            pass

        def interrupt_return(frame, event, _argument):
            if frame.f_code is helper_code and event == "return":
                descriptor = frame.f_locals["descriptor_holder"][0]
                captured.append(descriptor)
                self.assertEqual(os.fstat(descriptor).st_nlink, 0)
                raise ReturnInterrupt("injected unnamed-stage return interrupt")
            return interrupt_return

        previous_trace = sys.gettrace()
        sys.settrace(interrupt_return)
        try:
            with self.assertRaisesRegex(
                    ReturnInterrupt, "unnamed-stage return interrupt"):
                subject._commit_completed_timing_screen(
                    {}, directory, {"status": "complete"}, {}, {}, {},
                    mock.Mock())
        finally:
            sys.settrace(previous_trace)
        self.assertEqual(len(captured), 1)
        with self.assertRaises(OSError):
            os.fstat(captured[0])
        self.assertFalse(marker.exists())
        self.assertEqual(list(directory.iterdir()), [])

    def test_parent_fd_return_interrupt_is_closed_without_a_marker(self) -> None:
        directory = self.root / "parent-return-interrupt"
        directory.mkdir()
        marker = directory / "run-summary.json"
        helper_code = subject._open_completion_parent.__code__
        captured = []

        class ReturnInterrupt(BaseException):
            pass

        def interrupt_return(frame, event, _argument):
            if frame.f_code is helper_code and event == "return":
                descriptor = frame.f_locals["descriptor_holder"][0]
                captured.append(descriptor)
                self.assertTrue(stat.S_ISDIR(os.fstat(descriptor).st_mode))
                raise ReturnInterrupt("injected parent-fd return interrupt")
            return interrupt_return

        previous_trace = sys.gettrace()
        sys.settrace(interrupt_return)
        try:
            with self.assertRaisesRegex(
                    ReturnInterrupt, "parent-fd return interrupt"):
                subject._commit_completed_timing_screen(
                    {}, directory, {"status": "complete"}, {}, {}, {},
                    mock.Mock())
        finally:
            sys.settrace(previous_trace)
        self.assertEqual(len(captured), 1)
        with self.assertRaises(OSError):
            os.fstat(captured[0])
        self.assertFalse(marker.exists())

    def test_final_link_interrupt_leaves_strictly_loadable_marker(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("link-interrupt")
        marker = directory / "run-summary.json"
        marker.unlink()
        real_link = os.link

        class LinkInterrupt(BaseException):
            pass

        def link_then_interrupt(source, destination, *args, **kwargs):
            real_link(source, destination, *args, **kwargs)
            raise LinkInterrupt("injected after successful final link")

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject.os, "link", side_effect=link_then_interrupt), \
                self.assertRaisesRegex(LinkInterrupt, "successful final link"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertTrue(marker.is_file())
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_final_link_oserror_leaves_strictly_loadable_marker(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("link-oserror")
        marker = directory / "run-summary.json"
        marker.unlink()
        real_link = os.link

        def link_then_error(source, destination, *args, **kwargs):
            real_link(source, destination, *args, **kwargs)
            raise OSError(errno.EIO, "injected OSError after final link")

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject.os, "link", side_effect=link_then_error), \
                self.assertRaisesRegex(
                    subject.RunnerError, "OSError after final link"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertTrue(marker.is_file())
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_post_link_directory_fsync_failure_keeps_loadable_marker(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("link-fsync-failure")
        marker = directory / "run-summary.json"
        marker.unlink()
        real_fsync = os.fsync
        failed = False

        def fail_linked_directory_fsync(descriptor: int) -> None:
            nonlocal failed
            if (not failed and marker.exists() and
                    stat.S_ISDIR(os.fstat(descriptor).st_mode)):
                failed = True
                raise OSError(errno.EIO, "injected linked directory fsync")
            real_fsync(descriptor)

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject.os, "fsync",
                    side_effect=fail_linked_directory_fsync), \
                self.assertRaisesRegex(OSError, "linked directory fsync"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertTrue(failed)
        self.assertTrue(marker.is_file())
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_post_link_close_interrupt_drains_all_and_keeps_loadable_marker(
            self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("close-interrupt")
        marker = directory / "run-summary.json"
        marker.unlink()
        real_close = os.close
        close_interrupt = type(
            "CloseInterrupt", (BaseException,), {})(
                "injected after a successful descriptor close")
        linked_closes = []

        def close_then_interrupt(descriptor: int) -> None:
            linked = marker.exists()
            real_close(descriptor)
            if linked:
                linked_closes.append(descriptor)
                if len(linked_closes) == 1:
                    raise close_interrupt

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject.os, "close", side_effect=close_then_interrupt), \
                self.assertRaises(type(close_interrupt)) as raised:
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertIs(raised.exception, close_interrupt)
        self.assertIsInstance(
            raised.exception.__cause__, subject.RunnerError)
        self.assertEqual(
            len(linked_closes), len(subject.COMPLETED_DEPENDENCY_NAMES) + 2)
        self.assertTrue(marker.is_file())
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_preexisting_marker_is_preserved(self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("preexisting-marker")
        marker = directory / "run-summary.json"
        original = marker.read_bytes()
        with self.assertRaisesRegex(
                subject.RunnerError, "refusing to replace"):
            subject._commit_completed_timing_screen(
                contract, directory, run_summary, freeze,
                validated["summary"], validated["execution_receipt"],
                qualification)
        self.assertEqual(marker.read_bytes(), original)

    def test_public_loader_keeps_two_stable_bundle_reads(self) -> None:
        (contract, qualification, directory, freeze, _run_summary,
         validated) = self._completed_timing_fixture("two-loader-reads")
        real_read = subject._read_completed_bundle
        calls = []

        def record_read(source_names, directory_fd):
            calls.append(tuple(source_names.items()))
            return real_read(source_names, directory_fd)

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject, "_read_completed_bundle",
                    side_effect=record_read):
            subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(calls, [
            tuple(subject.COMPLETED_SOURCE_NAMES.items()),
            tuple(subject.COMPLETED_SOURCE_NAMES.items()),
        ])

    def test_public_loader_close_interrupt_drains_both_directory_fds(
            self) -> None:
        (contract, qualification, directory, freeze, run_summary,
         validated) = self._completed_timing_fixture("loader-close-interrupt")
        real_close = os.close
        real_open = os.open
        close_interrupt = type(
            "LoaderCloseInterrupt", (BaseException,), {})(
                "injected loader directory close interrupt")
        loader_directory_fds = []
        directory_closes = []

        def track_loader_directory_open(path, flags, *args, **kwargs):
            descriptor = real_open(path, flags, *args, **kwargs)
            if (str(path) == str(directory) and
                    flags & getattr(os, "O_DIRECTORY", 0)):
                loader_directory_fds.append(descriptor)
            return descriptor

        def close_then_interrupt(descriptor: int) -> None:
            is_loader_directory = (
                len(loader_directory_fds) == 2 and
                descriptor in loader_directory_fds)
            real_close(descriptor)
            if is_loader_directory:
                directory_closes.append(descriptor)
                if len(directory_closes) == 1:
                    raise close_interrupt

        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze), \
                mock.patch.object(
                    subject.os, "open",
                    side_effect=track_loader_directory_open), \
                mock.patch.object(
                    subject.os, "close", side_effect=close_then_interrupt), \
                self.assertRaises(type(close_interrupt)) as raised:
            subject.load_completed_timing_screen(contract, directory)
        self.assertIs(raised.exception, close_interrupt)
        self.assertEqual(len(directory_closes), 2)
        self.assertIsInstance(
            raised.exception.__cause__, subject.RunnerError)
        with mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification), mock.patch.object(
                native_api, "validate_execution_receipt",
                return_value=validated), mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=freeze):
            loaded = subject.load_completed_timing_screen(contract, directory)
        self.assertEqual(loaded["run_summary"], run_summary)

    def test_completion_publication_has_no_named_rollback_machinery(self) -> None:
        source = inspect.getsource(subject)
        for forbidden in (
                "_CompletionMarkerTransaction", "_PublishedCompletionMarker",
                "_rollback_completion_marker", "_rename_noreplace_at",
                "renameat2", "quarantine"):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)

    def test_runner_has_no_explicit_failure_point_after_summary_commit(
            self) -> None:
        source = inspect.getsource(subject.run_short_screen)
        tail = source.split("_commit_completed_timing_screen(", 1)[1]
        before_return = tail.split("return summary", 1)[0]
        self.assertNotIn("_remaining(", before_return)
        self.assertNotIn("_git_head(", before_return)
        self.assertNotIn("fail(", before_return)

    def test_timing_cleanup_restores_affinity_after_terminate_oserror(
            self) -> None:
        events = []
        primary = subject.RunnerError("injected primary failure")

        def terminate(_workers):
            events.append("terminate")
            raise OSError(errno.EIO, "injected terminate failure")

        def restore(_affinity):
            events.append("restore")
            raise subject.RunnerError("injected restore failure")

        with mock.patch.object(
                subject, "terminate_workers", side_effect=terminate), \
                mock.patch.object(
                    subject, "_restore_controller_affinity",
                    side_effect=restore), self.assertRaises(
                        subject.RunnerError) as raised:
            subject._finish_timing_cleanup(
                [mock.Mock()], False, True, {7}, primary)
        self.assertEqual(events, ["terminate", "restore"])
        self.assertIs(raised.exception.__cause__, primary)
        message = str(raised.exception)
        self.assertIn("injected primary failure", message)
        self.assertIn("injected terminate failure", message)
        self.assertIn("injected restore failure", message)

    def test_timing_cleanup_does_not_swallow_control_flow(self) -> None:
        class CustomInterrupt(BaseException):
            pass

        for control in (KeyboardInterrupt("interrupt"), SystemExit(7),
                        CustomInterrupt("custom interrupt")):
            events = []

            def terminate(_workers, _control=control):
                events.append("terminate")
                raise _control

            def restore(_affinity):
                events.append("restore")

            with self.subTest(control=type(control).__name__), \
                    mock.patch.object(
                        subject, "terminate_workers", side_effect=terminate), \
                    mock.patch.object(
                        subject, "_restore_controller_affinity",
                        side_effect=restore), self.assertRaises(
                            type(control)) as raised:
                subject._finish_timing_cleanup(
                    [mock.Mock()], False, True, {7}, None)
            self.assertEqual(events, ["terminate", "restore"])
            self.assertIs(raised.exception, control)
            self.assertIsInstance(
                raised.exception.__cause__, subject.RunnerError)

    def test_spawn_cleanup_preserves_primary_control_after_terminate_oserror(
            self) -> None:
        class SpawnInterrupt(BaseException):
            pass

        primary = SpawnInterrupt("injected second-spawn interrupt")
        first_process = mock.Mock()
        cleanup_calls = []

        def terminate(workers):
            cleanup_calls.append(list(workers))
            raise OSError(errno.EIO, "injected partial-pool cleanup failure")

        with mock.patch.object(
                subject.subprocess, "Popen",
                side_effect=(first_process, primary)), mock.patch.object(
                subject, "terminate_workers", side_effect=terminate), \
                self.assertRaises(SpawnInterrupt) as raised:
            subject.spawn_workers(
                {"resolved_path": "/injected/worker"}, [3, 5],
                time.monotonic() + 1.0)
        self.assertIs(raised.exception, primary)
        self.assertEqual(len(cleanup_calls), 1)
        self.assertEqual(len(cleanup_calls[0]), 1)
        self.assertIs(cleanup_calls[0][0].process, first_process)
        self.assertIsInstance(
            raised.exception.__cause__, subject.RunnerError)
        self.assertIn(
            "partial-pool cleanup failure", str(raised.exception.__cause__))

    def test_spawn_constructor_interrupt_cleans_provisional_child(self) -> None:
        class ConstructorInterrupt(BaseException):
            pass

        primary = ConstructorInterrupt("injected worker constructor interrupt")
        process = mock.Mock()
        cleanup_calls = []

        def terminate(workers):
            cleanup_calls.append(list(workers))

        with mock.patch.object(
                subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject, "PersistentWorker", side_effect=primary), \
                mock.patch.object(
                    subject, "terminate_workers", side_effect=terminate), \
                self.assertRaises(ConstructorInterrupt) as raised:
            subject.spawn_workers(
                {"resolved_path": "/injected/worker"}, [3],
                time.monotonic() + 1.0)
        self.assertIs(raised.exception, primary)
        self.assertEqual(len(cleanup_calls), 1)
        self.assertEqual(len(cleanup_calls[0]), 1)
        self.assertIs(cleanup_calls[0][0].process, process)

    def test_spawn_first_post_assignment_interrupt_cleans_child_once(
            self) -> None:
        class AssignmentInterrupt(BaseException):
            pass

        primary = AssignmentInterrupt(
            "injected first post-Popen-assignment interrupt")
        process = mock.Mock()
        cleanup_calls = []
        code = subject.spawn_workers.__code__

        def terminate(workers):
            cleanup_calls.append(list(workers))

        def interrupt_after_assignment(frame, event, _argument):
            provisional = frame.f_locals.get("provisional_worker")
            if (frame.f_code is code and event == "line" and
                    provisional is not None and
                    provisional.process is process and
                    not frame.f_locals.get("workers")):
                raise primary
            return interrupt_after_assignment

        previous_trace = sys.gettrace()
        with mock.patch.object(
                subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject, "terminate_workers", side_effect=terminate):
            sys.settrace(interrupt_after_assignment)
            try:
                with self.assertRaises(AssignmentInterrupt) as raised:
                    subject.spawn_workers(
                        {"resolved_path": "/injected/worker"}, [3],
                        time.monotonic() + 1.0)
            finally:
                sys.settrace(previous_trace)
        self.assertIs(raised.exception, primary)
        self.assertEqual(len(cleanup_calls), 1)
        self.assertEqual(len(cleanup_calls[0]), 1)
        self.assertIs(cleanup_calls[0][0].process, process)

    def test_spawn_append_handoff_interrupt_cleans_child_once(self) -> None:
        class AppendInterrupt(BaseException):
            pass

        primary = AppendInterrupt("injected worker append handoff interrupt")
        process = mock.Mock()
        cleanup_calls = []
        code = subject.spawn_workers.__code__

        def terminate(workers):
            cleanup_calls.append(list(workers))

        def interrupt_handoff(frame, event, _argument):
            if (frame.f_code is code and event == "line" and
                    len(frame.f_locals.get("workers", ())) == 1 and
                    frame.f_locals.get("provisional_worker") is not None):
                raise primary
            return interrupt_handoff

        previous_trace = sys.gettrace()
        with mock.patch.object(
                subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject, "terminate_workers", side_effect=terminate):
            sys.settrace(interrupt_handoff)
            try:
                with self.assertRaises(AppendInterrupt) as raised:
                    subject.spawn_workers(
                        {"resolved_path": "/injected/worker"}, [3],
                        time.monotonic() + 1.0)
            finally:
                sys.settrace(previous_trace)
        self.assertIs(raised.exception, primary)
        self.assertEqual(len(cleanup_calls), 1)
        self.assertEqual(len(cleanup_calls[0]), 1)
        self.assertIs(cleanup_calls[0][0].process, process)

    def test_qualification_cleanup_control_outranks_ordinary_primary(
            self) -> None:
        class CleanupInterrupt(BaseException):
            pass

        primary = subject.RunnerError("injected qualification failure")
        cleanup = CleanupInterrupt("injected qualification cleanup interrupt")
        with mock.patch.object(
                subject, "terminate_workers", side_effect=cleanup), \
                self.assertRaises(CleanupInterrupt) as raised:
            subject._finish_worker_pool_cleanup(
                [mock.Mock()], False,
                "qualification worker cleanup failed", primary)
        self.assertIs(raised.exception, cleanup)
        self.assertIsInstance(
            raised.exception.__cause__, subject.RunnerError)
        diagnostic = str(raised.exception.__cause__)
        self.assertIn("injected qualification failure", diagnostic)
        self.assertIn("qualification cleanup interrupt", diagnostic)

    def test_qualification_pool_is_quit_and_reaped_before_exact_eight_spawn(
            self) -> None:
        contract = contract_api.load_contract()
        qualification_offsets = [0] * len(list(
            contract_api.iter_timing_base_cells(contract, "development")))
        qualification_offsets[7] = 1
        qualification_offsets[-1] = 2
        qualification = self._qualification(
            contract, retry_offsets=qualification_offsets)

        class Process:
            def __init__(self) -> None:
                self.status = None

            def poll(self):
                return self.status

        qualification_workers = [
            mock.Mock(cpu=cpu, process=Process()) for cpu in range(12)
        ]
        timing_workers = [
            mock.Mock(cpu=cpu, process=Process()) for cpu in range(8)
        ]
        spawn_count = 0
        lifecycle = []

        def spawn(_description, cpus, _deadline):
            nonlocal spawn_count
            spawn_count += 1
            if spawn_count == 1:
                self.assertEqual(cpus, list(range(12)))
                lifecycle.append("qualification_spawn")
                return qualification_workers
            self.assertEqual(spawn_count, 2)
            self.assertTrue(all(
                worker.process.poll() is not None
                for worker in qualification_workers))
            self.assertEqual(cpus, list(range(8)))
            lifecycle.append("timing_spawn")
            return timing_workers

        def quit_pool(workers, _deadline):
            for worker in workers:
                worker.process.status = 0
            lifecycle.append(
                "qualification_reaped" if workers is qualification_workers
                else "timing_reaped")

        def assemble_qualification(*_args, **_kwargs):
            self.assertTrue(all(
                worker.process.poll() is None
                for worker in qualification_workers))
            lifecycle.append("qualification_receipt_published_live")
            return qualification, qualification_metadata, "f" * 64

        def restore_affinity(_affinity):
            lifecycle.append("affinity_restored")

        published_summary = {}

        def commit_summary(_contract, _directory, value, *_args, **_kwargs):
            self.assertTrue(all(
                worker.process.poll() is not None
                for worker in qualification_workers + timing_workers))
            self.assertIn("affinity_restored", lifecycle)
            lifecycle.append("summary_published")
            published_summary.update(value)

        description = {
            "binary_sha256": "a" * 64,
            "source_git_commit": "1" * 40,
        }
        qualification_metadata = {
            "qualification_attempt_count":
                len(qualification_offsets) + sum(qualification_offsets),
            "qualification_worker_cpus": list(range(12)),
        }
        receipt = {
            "record_count": 1,
            "freeze_manifest_sha256": "b" * 64,
            "result_stream_sha256": "c" * 64,
            "receipt_sha256": "d" * 64,
            "validator_summary_sha256":
                contract_api.sha256_json({"validated": True}),
            "thermal": {
                "sample_count": 2,
                "cpu_tctl_max_millic": 60000,
                "dimm_max_millic": 40000,
            },
        }
        frozen = {"trace_manifest_sha256": "e" * 64}
        def assembled(_contract, kind, *_args, **_kwargs):
            self.assertEqual(kind, "timing")
            value = dict(receipt)
            value["qualification_thermal"] = {
                "sample_count": 3,
                "cpu_tctl_max_millic": 70000,
                "dimm_max_millic": 50000,
            }
            return {
                "summary": {"validated": True},
                "execution_receipt": value,
            }

        validation_modes = []

        def validate(*args, **kwargs):
            live = kwargs.get("verify_live_sampler")
            validation_modes.append(live)
            if live is False:
                self.assertIn("timing_reaped", lifecycle)
                self.assertIn("affinity_restored", lifecycle)
                lifecycle.append("post_shutdown_validation")
            return assembled(*args, **kwargs)

        args = mock.Mock(
            deadline_seconds=60.0, contract=None, worker=Path("worker"),
            sampler_cpu=99, sampler_pid=123, sampler_script=Path("s.py"),
            sampler_csv=Path("thermal.csv"), cpus=None, controller_cpu=None,
            pause_sibling_pids=None,
            output_dir=self.root / "run")
        managed_paused = mock.Mock(records=[])
        managed_paused.resume.side_effect = \
            lambda _deadline: lifecycle.append("sibling_load_resumed")
        sibling_guard = mock.Mock()
        sibling_guard.attestation.return_value = {"terminal_status": "clean"}
        patches = (
            mock.patch.object(
                contract_api, "load_contract", return_value=contract),
            mock.patch.object(
                contract_api, "load_timing_qualification_map",
                return_value=qualification),
            mock.patch.object(
                os, "sched_getaffinity", return_value=set(range(12))),
            mock.patch.object(
                subject, "_cpu_topology",
                side_effect=lambda cpu, _root: (0, cpu)),
            mock.patch.object(
                subject, "describe_worker", return_value=description),
            mock.patch.object(subject, "_development_timing_shape"),
            mock.patch.object(subject, "_preflight_sampler"),
            mock.patch.object(subject, "_git_head", return_value="1" * 40),
            mock.patch.object(subject, "_require_worker_source_commit"),
            mock.patch.object(
                subject, "_timing_qualification_controls", return_value=[]),
            mock.patch.object(
                subject, "choose_new_sampler_start",
                side_effect=(1000000000, 3000000000)),
            mock.patch.object(subject, "spawn_workers", side_effect=spawn),
            mock.patch.object(
                subject, "run_timing_qualification",
                return_value=(self.root / "qualification-native.jsonl",
                              1500000000)),
            mock.patch.object(
                native_api, "assemble_timing_qualification",
                side_effect=assemble_qualification),
            mock.patch.object(subject, "quit_workers", side_effect=quit_pool),
            mock.patch.object(
                subject, "_wait_for_sampler_sample",
                side_effect=(2000000000, 9000000000)),
            mock.patch.object(native_api, "write_sampler_attestation"),
            mock.patch.object(
                subject, "select_cpu_layout",
                return_value=(list(range(8)), 8)),
            mock.patch.object(
                native_api, "timing_sibling_cpus", return_value=[]),
            mock.patch.object(
                subject, "ManagedPausedSiblingProcesses",
                return_value=managed_paused),
            mock.patch.object(
                subject, "SiblingIsolationGuard",
                return_value=sibling_guard),
            mock.patch.object(
                native_api, "publish_timing_trace_manifest",
                return_value="e" * 64),
            mock.patch.object(
                subject, "write_development_timing_freeze",
                return_value=frozen),
            mock.patch.object(subject, "_pin_controller"),
            mock.patch.object(
                subject, "_restore_controller_affinity",
                side_effect=restore_affinity),
            mock.patch.object(
                subject, "_run_timing_jobs",
                return_value=(self.root / "timing-native.jsonl",
                              8000000000)),
            mock.patch.object(
                native_api, "assemble_results", side_effect=assembled),
            mock.patch.object(
                native_api, "validate_execution_receipt",
                side_effect=validate),
            mock.patch.object(
                contract_api, "load_freeze_manifest", return_value=frozen),
            mock.patch.object(
                contract_api, "architecture_artifact_sha256",
                return_value="9" * 64),
            mock.patch.object(
                subject, "_commit_completed_timing_screen",
                side_effect=commit_summary),
        )
        with ExitStack() as stack:
            for patcher in patches:
                stack.enter_context(patcher)
            result = subject.run_short_screen(args)
        self.assertEqual(published_summary, result)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["cpu_tctl_max_millic"], 60000)
        self.assertEqual(result["dimm_max_millic"], 40000)
        self.assertEqual(result["qualification_cpu_tctl_max_millic"], 70000)
        self.assertEqual(result["qualification_dimm_max_millic"], 50000)
        self.assertEqual(result["overall_cpu_tctl_max_millic"], 70000)
        self.assertEqual(result["overall_dimm_max_millic"], 50000)
        self.assertEqual(
            result["qualification_attempt_count"],
            len(qualification_offsets) + 3)
        self.assertEqual(result["qualification_retried_cell_count"], 2)
        self.assertEqual(result["qualification_max_retry_offset"], 2)
        self.assertEqual(result["qualification_sum_retry_offsets"], 3)
        self.assertEqual(spawn_count, 2)
        self.assertEqual(validation_modes, [True, False])
        self.assertTrue(all(
            worker.process.poll() is not None for worker in timing_workers))
        self.assertLess(
            lifecycle.index("qualification_receipt_published_live"),
            lifecycle.index("qualification_reaped"))
        self.assertLess(
            lifecycle.index("timing_reaped"),
            lifecycle.index("affinity_restored"))
        self.assertLess(
            lifecycle.index("affinity_restored"),
            lifecycle.index("sibling_load_resumed"))
        self.assertLess(
            lifecycle.index("sibling_load_resumed"),
            lifecycle.index("post_shutdown_validation"))
        self.assertLess(
            lifecycle.index("post_shutdown_validation"),
            lifecycle.index("summary_published"))

    def test_affinity_restore_failure_never_publishes_complete_summary(
            self) -> None:
        contract = contract_api.load_contract()
        qualification = self._qualification(contract)

        class Process:
            def __init__(self) -> None:
                self.status = None

            def poll(self):
                return self.status

        for message in ("restore failure", "hard wall during restore"):
            qualification_workers = [
                mock.Mock(cpu=cpu, process=Process()) for cpu in range(12)
            ]
            timing_workers = [
                mock.Mock(cpu=cpu, process=Process()) for cpu in range(8)
            ]
            pools = iter((qualification_workers, timing_workers))

            def quit_pool(workers, _deadline):
                for worker in workers:
                    worker.process.status = 0

            description = {
                "binary_sha256": "a" * 64,
                "source_git_commit": "1" * 40,
            }
            receipt = {
                "record_count": 1,
                "freeze_manifest_sha256": "b" * 64,
                "result_stream_sha256": "c" * 64,
                "receipt_sha256": "d" * 64,
                "validator_summary_sha256":
                    contract_api.sha256_json({"validated": True}),
                "thermal": {
                    "sample_count": 2,
                    "cpu_tctl_max_millic": 60000,
                    "dimm_max_millic": 40000,
                },
                "qualification_thermal": {
                    "sample_count": 2,
                    "cpu_tctl_max_millic": 61000,
                    "dimm_max_millic": 41000,
                },
            }
            assembled_evidence = {
                "summary": {"validated": True},
                "execution_receipt": receipt,
            }
            args = mock.Mock(
                deadline_seconds=60.0, contract=None, worker=Path("worker"),
                sampler_cpu=99, sampler_pid=123,
                sampler_script=Path("s.py"),
                sampler_csv=Path("thermal.csv"), cpus=None,
                controller_cpu=None,
                pause_sibling_pids=None,
                output_dir=self.root / ("run-" + message.replace(" ", "-")))
            summary_writer = mock.Mock()
            managed_paused = mock.Mock(records=[])
            sibling_guard = mock.Mock()
            sibling_guard.attestation.return_value = {
                "terminal_status": "clean"}
            patches = (
                mock.patch.object(
                    contract_api, "load_contract", return_value=contract),
                mock.patch.object(
                    contract_api, "load_timing_qualification_map",
                    return_value=qualification),
                mock.patch.object(
                    os, "sched_getaffinity", return_value=set(range(12))),
                mock.patch.object(
                    subject, "_cpu_topology",
                    side_effect=lambda cpu, _root: (0, cpu)),
                mock.patch.object(
                    subject, "describe_worker", return_value=description),
                mock.patch.object(subject, "_development_timing_shape"),
                mock.patch.object(subject, "_preflight_sampler"),
                mock.patch.object(
                    subject, "_git_head", return_value="1" * 40),
                mock.patch.object(subject, "_require_worker_source_commit"),
                mock.patch.object(
                    subject, "_create_output_dir", return_value=args.output_dir),
                mock.patch.object(
                    subject, "_timing_qualification_controls",
                    return_value=[]),
                mock.patch.object(
                    subject, "choose_new_sampler_start",
                    side_effect=(1000000000, 3000000000)),
                mock.patch.object(
                    subject, "spawn_workers",
                    side_effect=lambda *_args: next(pools)),
                mock.patch.object(
                    subject, "run_timing_qualification",
                    return_value=(Path("qualification-native"), 1500000000)),
                mock.patch.object(
                    native_api, "assemble_timing_qualification",
                    return_value=(qualification, {}, "f" * 64)),
                mock.patch.object(
                    subject, "quit_workers", side_effect=quit_pool),
                mock.patch.object(
                    subject, "_wait_for_sampler_sample",
                    side_effect=(2000000000, 9000000000)),
                mock.patch.object(native_api, "write_sampler_attestation"),
                mock.patch.object(
                    subject, "select_cpu_layout",
                    return_value=(list(range(8)), 8)),
                mock.patch.object(
                    native_api, "timing_sibling_cpus", return_value=[]),
                mock.patch.object(
                    subject, "ManagedPausedSiblingProcesses",
                    return_value=managed_paused),
                mock.patch.object(
                    subject, "SiblingIsolationGuard",
                    return_value=sibling_guard),
                mock.patch.object(
                    native_api, "publish_timing_trace_manifest",
                    return_value="e" * 64),
                mock.patch.object(
                    subject, "write_development_timing_freeze",
                    return_value={"trace_manifest_sha256": "e" * 64}),
                mock.patch.object(subject, "_pin_controller"),
                mock.patch.object(
                    subject, "_restore_controller_affinity",
                    side_effect=subject.RunnerError(message)),
                mock.patch.object(
                    subject, "_run_timing_jobs",
                    return_value=(Path("t"), 8000000000)),
                mock.patch.object(
                    native_api, "assemble_results",
                    return_value=assembled_evidence),
                mock.patch.object(
                    native_api, "validate_execution_receipt",
                    return_value=assembled_evidence),
                mock.patch.object(
                    subject, "_atomic_write_object", summary_writer),
            )
            with self.subTest(message=message), ExitStack() as stack:
                for patcher in patches:
                    stack.enter_context(patcher)
                with self.assertRaisesRegex(subject.RunnerError, message):
                    subject.run_short_screen(args)
            summary_writer.assert_not_called()

    def test_strict_timing_validator_precomputes_cell_indexes_once(self) -> None:
        worker = self._fake_worker()
        description = self._description(worker)
        contract = contract_api.load_contract()
        output = self.root / "output"
        output.mkdir()
        qualification = self._qualification(contract, description=description)
        freeze = subject.write_development_timing_freeze(
            contract, description, [0], 1, "1" * 40,
            "3" * 64, output,
            qualification)
        with mock.patch.object(
                contract_api, "_timing_cell_indexes",
                wraps=contract_api._timing_cell_indexes) as cell_indexes:
            validator = subject._strict_response_validator(
                contract, freeze, "timing", description, 0,
                qualification)
            self.assertTrue(callable(validator))
            self.assertEqual(cell_indexes.call_count, 1)

    def test_post_link_publish_failure_preserves_visible_commit(self) -> None:
        destination = self.root / "atomic-output.json"
        real_fsync = os.fsync

        def fail_directory_fsync(descriptor: int) -> None:
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise OSError("injected directory fsync failure")
            real_fsync(descriptor)

        with mock.patch.object(os, "fsync", side_effect=fail_directory_fsync):
            with self.assertRaises(subject.RunnerError):
                subject._atomic_write_bytes(destination, b"complete\n")
        self.assertEqual(destination.read_bytes(), b"complete\n")

    def test_atomic_write_fdopen_failure_closes_and_unlinks_staging(
            self) -> None:
        destination = self.root / "atomic-fdopen-failure.json"
        descriptors = []

        def reject_fdopen(descriptor: int, _mode: str):
            descriptors.append(descriptor)
            raise OSError("injected fdopen failure")

        with mock.patch.object(
                subject.os, "fdopen", side_effect=reject_fdopen), \
                self.assertRaisesRegex(OSError, "injected fdopen failure"):
            subject._atomic_write_bytes(destination, b"complete\n")
        self.assertEqual(len(descriptors), 1)
        with self.assertRaises(OSError) as raised:
            os.fstat(descriptors[0])
        self.assertEqual(raised.exception.errno, errno.EBADF)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_atomic_line_sink_fdopen_failure_closes_and_unlinks_staging(
            self) -> None:
        destination = self.root / "sink-fdopen-failure.jsonl"
        descriptors = []

        def reject_fdopen(descriptor: int, _mode: str):
            descriptors.append(descriptor)
            raise OSError("injected fdopen failure")

        with mock.patch.object(
                subject.os, "fdopen", side_effect=reject_fdopen), \
                self.assertRaisesRegex(OSError, "injected fdopen failure"):
            subject.AtomicLineSink(destination)
        self.assertEqual(len(descriptors), 1)
        with self.assertRaises(OSError) as raised:
            os.fstat(descriptors[0])
        self.assertEqual(raised.exception.errno, errno.EBADF)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_atomic_line_sink_field_failure_closes_and_unlinks_staging(
            self) -> None:
        class InjectedInterrupt(BaseException):
            pass

        class FailingSink(subject.AtomicLineSink):
            def __setattr__(self, name, value):
                if name == "sha256":
                    raise InjectedInterrupt()
                object.__setattr__(self, name, value)

        destination = self.root / "sink-field-failure.jsonl"
        descriptors = []
        real_fdopen = os.fdopen

        def capture_fdopen(descriptor: int, mode: str):
            descriptors.append(descriptor)
            return real_fdopen(descriptor, mode)

        with mock.patch.object(
                subject.os, "fdopen", side_effect=capture_fdopen), \
                self.assertRaises(InjectedInterrupt):
            FailingSink(destination)
        self.assertEqual(len(descriptors), 1)
        with self.assertRaises(OSError) as raised:
            os.fstat(descriptors[0])
        self.assertEqual(raised.exception.errno, errno.EBADF)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_post_link_deadline_preserves_visible_commit(self) -> None:
        destination = self.root / "deadline-output.json"
        error = subject.RunnerError("injected hard-wall deadline")
        real_fsync = os.fsync

        def fail_directory_fsync(descriptor: int) -> None:
            if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise error
            real_fsync(descriptor)

        with mock.patch.object(os, "fsync", side_effect=fail_directory_fsync):
            with self.assertRaisesRegex(subject.RunnerError, "hard-wall"):
                subject._atomic_write_bytes(destination, b"complete\n")
        self.assertEqual(destination.read_bytes(), b"complete\n")

    def test_atomic_callers_preserve_post_publish_commit_by_inode(self) \
            -> None:
        class InjectedInterrupt(BaseException):
            pass

        for writer in ("atomic", "sink"):
            for replacement in (False, True):
                with self.subTest(writer=writer, replacement=replacement):
                    destination = self.root / "{}-{}.jsonl".format(
                        writer, "replacement" if replacement else "own")
                    real_publish = subject._publish_staged

                    def publish_then_interrupt(
                            staged, target, expected_identity=None):
                        real_publish(staged, target, expected_identity)
                        if replacement:
                            target.unlink()
                            target.write_bytes(b"replacement winner\n")
                        raise InjectedInterrupt()

                    with mock.patch.object(
                            subject, "_publish_staged",
                            side_effect=publish_then_interrupt):
                        if writer == "atomic":
                            with self.assertRaises(InjectedInterrupt):
                                subject._atomic_write_bytes(
                                    destination, b"new artifact\n")
                        else:
                            sink = subject.AtomicLineSink(destination)
                            try:
                                sink.write(b"new artifact\n")
                                with self.assertRaises(InjectedInterrupt):
                                    sink.publish()
                            finally:
                                sink.abort()
                    if replacement:
                        self.assertEqual(
                            destination.read_bytes(), b"replacement winner\n")
                    else:
                        self.assertEqual(
                            destination.read_bytes(), b"new artifact\n")

        for writer in ("atomic", "sink"):
            with self.subTest(writer=writer, preexisting=True):
                destination = self.root / (writer + "-preexisting.jsonl")
                destination.write_bytes(b"preexisting winner\n")
                if writer == "atomic":
                    with self.assertRaises(subject.RunnerError):
                        subject._atomic_write_bytes(
                            destination, b"new artifact\n")
                else:
                    sink = subject.AtomicLineSink(destination)
                    try:
                        sink.write(b"new artifact\n")
                        with self.assertRaises(subject.RunnerError):
                            sink.publish()
                    finally:
                        sink.abort()
                self.assertEqual(
                    destination.read_bytes(), b"preexisting winner\n")

    def test_exact_recovery_and_homogeneous_timing_wave_domains(self) -> None:
        contract = contract_api.load_contract()
        recovery = subject._recovery_jobs(contract, 3)
        self.assertEqual(len(recovery), 1080)
        self.assertEqual(
            {job.ordinal for job in recovery}, set(range(1080)))
        panels = contract_api.timing_panels(
            contract, [value[0] for value in subject.EXPECTED_ARMS])
        for repetitions in (16, 24):
            with self.subTest(repetitions=repetitions):
                candidate = self._timing_contract(repetitions)
                qualification = self._qualification(candidate)
                bounded_cells = candidate["timing"]["domains"][
                    "development"]["expected_cells"]
                with mock.patch.object(
                        subject, "FROZEN_TIMING_BASE_CELLS", bounded_cells):
                    waves = subject._timing_job_waves(
                        candidate, len(panels), qualification)
                cells = list(contract_api.iter_timing_cells(
                    candidate, "development", qualification))
                round_count = repetitions // subject.TIMING_WORKER_COUNT
                self.assertEqual(len(waves), 264 * round_count)
                self.assertTrue(all(
                    len(jobs) == subject.TIMING_WORKER_COUNT
                    for _, jobs in waves))
                self.assertEqual(
                    {job.ordinal for _, jobs in waves for job in jobs},
                    set(range(len(cells) * len(panels))))
                commands = [job.command() for _, jobs in waves for job in jobs]
                self.assertEqual(len(commands), len(set(commands)))

                worker_counts = {
                    replicate: [0] * subject.TIMING_WORKER_COUNT
                    for replicate in range(repetitions)
                }
                for cohort_index in range(264):
                    stable_index, panel = divmod(cohort_index, len(panels))
                    width_index, k_index = divmod(
                        stable_index, len(contract_api.EXPECTED_TIMING_SHORT_K))
                    for round_index in range(round_count):
                        rotation, jobs = waves[
                            round_index * 264 + cohort_index]
                        replicates = [
                            cells[job.cell]["replicate"] for job in jobs
                        ]
                        self.assertEqual(replicates, list(range(
                            round_index * subject.TIMING_WORKER_COUNT,
                            (round_index + 1) * subject.TIMING_WORKER_COUNT)))
                        self.assertEqual({job.item for job in jobs}, {panel})
                        self.assertEqual(
                            {cells[job.cell]["block_bytes"] for job in jobs},
                            {candidate["timing"]["domains"]["development"]
                             ["block_bytes"][width_index]})
                        self.assertEqual(
                            {cells[job.cell]["K"] for job in jobs},
                            {contract_api.EXPECTED_TIMING_SHORT_K[k_index]})
                        orders = [contract_api.timing_order(
                            panels[panel], replicate)
                            for replicate in replicates]
                        self.assertEqual(orders.count("ABBA"), 4)
                        self.assertEqual(orders.count("BAAB"), 4)
                        expected_rotation = contract_api.timing_worker_slot(
                            candidate, "development",
                            [value[0] for value in subject.EXPECTED_ARMS],
                            cohort_index,
                            replicates[0])
                        self.assertEqual(rotation, expected_rotation)
                        for position, replicate in enumerate(replicates):
                            slot = (position + rotation) % \
                                subject.TIMING_WORKER_COUNT
                            self.assertEqual(slot,
                                contract_api.timing_worker_slot(
                                    candidate, "development",
                                    [value[0] for value in
                                     subject.EXPECTED_ARMS],
                                    cohort_index, replicate))
                            worker_counts[replicate][slot] += 1
                self.assertTrue(all(
                    counts == [33] * subject.TIMING_WORKER_COUNT
                    for counts in worker_counts.values()))

                with mock.patch.object(
                        subject, "FROZEN_TIMING_BASE_CELLS", bounded_cells):
                    rebuilt = subject._timing_job_waves(
                        candidate, len(panels), qualification)
                self.assertEqual(
                    [(rotation, [job.command() for job in jobs])
                     for rotation, jobs in waves],
                    [(rotation, [job.command() for job in jobs])
                     for rotation, jobs in rebuilt])

    def test_timing_waves_pack_each_cells_own_selected_retry(self) -> None:
        contract = self._timing_contract(16)
        base_cells = list(contract_api.iter_timing_base_cells(
            contract, "development"))
        stable_fields = [
            field for field in contract["timing"]["cell_key"]
            if field not in (
                "replicate", "base_loss_seed", "base_cell_sha256",
                "loss_retry_offset", "loss_seed")
        ]
        first_identity = {
            field: base_cells[0][field] for field in stable_fields
        }
        first_cohort = [
            ordinal for ordinal, cell in enumerate(base_cells)
            if all(cell[field] == value
                   for field, value in first_identity.items())
        ]
        self.assertEqual(len(first_cohort), 16)
        offsets = [0] * len(base_cells)
        expected_first = [0, 1, 2, 3, 0, 1, 2, 3]
        for ordinal, retry in zip(first_cohort[:8], expected_first):
            offsets[ordinal] = retry
        qualification = self._qualification(contract, offsets)
        panels = contract_api.timing_panels(
            contract, [value[0] for value in subject.EXPECTED_ARMS])
        cells = list(contract_api.iter_timing_cells(
            contract, "development", qualification))
        with mock.patch.object(
                subject, "FROZEN_TIMING_BASE_CELLS",
                contract["timing"]["domains"]["development"][
                    "expected_cells"]):
            waves = subject._timing_job_waves(
                contract, len(panels), qualification)
        first_wave = waves[0][1]
        self.assertEqual(
            [job.retry_offset for job in first_wave], expected_first)
        self.assertGreater(len({job.retry_offset for job in first_wave}), 1)
        for _, jobs in waves:
            for job in jobs:
                self.assertEqual(
                    job.retry_offset,
                    cells[job.cell]["loss_retry_offset"])
                prefix, cell, packed = job.command().decode("ascii").split()
                self.assertEqual(prefix, "T")
                self.assertEqual(int(cell), job.cell)
                panel, retry = divmod(int(packed), 256)
                self.assertEqual(panel, job.item)
                self.assertEqual(retry, job.retry_offset)

    def test_timing_wave_builder_fails_closed_on_domain_mutations(self) -> None:
        contract = self._timing_contract(16)
        qualification = self._qualification(contract)
        panel_count = len(contract_api.timing_panels(
            contract, [value[0] for value in subject.EXPECTED_ARMS]))
        original = list(contract_api.iter_timing_cells(
            contract, "development", qualification))
        mutations = {}
        mutations["missing"] = original[:-1]
        mutations["extra"] = original + [dict(original[-1])]
        duplicate = [dict(cell) for cell in original]
        duplicate[-1] = dict(duplicate[0])
        mutations["duplicate"] = duplicate
        identity_drift = [dict(cell) for cell in original]
        identity_drift[13]["K"] += 1
        mutations["identity_drift"] = identity_drift
        bad_replicate = [dict(cell) for cell in original]
        bad_replicate[0]["replicate"] = 16
        mutations["bad_replicate"] = bad_replicate
        for name, cells in mutations.items():
            with self.subTest(name=name), mock.patch.object(
                    subject, "FROZEN_TIMING_BASE_CELLS", len(original)), \
                    mock.patch.object(
                        contract_api, "iter_timing_cells",
                        side_effect=lambda _contract, _phase, _qualification,
                        values=cells:
                            iter(values)):
                with self.assertRaises(subject.RunnerError):
                    subject._timing_job_waves(
                        contract, panel_count, qualification)
        for bad_panel_count in (True, panel_count - 1, panel_count + 1):
            with self.subTest(panel_count=bad_panel_count), mock.patch.object(
                    subject, "FROZEN_TIMING_BASE_CELLS", len(original)):
                with self.assertRaises(subject.RunnerError):
                    subject._timing_job_waves(
                        contract, bad_panel_count, qualification)

    def test_timing_wave_builder_rejects_multi_candidate_panel_drift(self) \
            -> None:
        contract = contract_api.load_contract()
        qualification = self._qualification(contract)
        expanded_arms = subject.EXPECTED_ARMS + (
            ("candidate_two", "wirehair2_experiment"),
        )
        with mock.patch.object(subject, "EXPECTED_ARMS", expanded_arms):
            roster = [value[0] for value in expanded_arms]
            panel_count = len(contract_api.timing_panels(contract, roster))
            self.assertEqual(panel_count, 18)
            with self.assertRaisesRegex(
                    subject.RunnerError, "panel cardinality differs"):
                subject._timing_job_waves(
                    contract, panel_count, qualification)

    def test_development_timing_shape_rejects_artifact_bound_drift(
            self) -> None:
        contract = contract_api.load_contract()
        panel_count = len(contract_api.timing_panels(
            contract, [value[0] for value in subject.EXPECTED_ARMS]))
        self.assertEqual(
            subject._development_timing_shape(contract, panel_count)[2],
            subject.FROZEN_TIMING_BASE_CELLS)
        for expected_cells in (
                subject.FROZEN_TIMING_BASE_CELLS - 1,
                subject.FROZEN_TIMING_BASE_CELLS + 1):
            candidate = copy.deepcopy(contract)
            candidate["timing"]["domains"]["development"][
                "expected_cells"] = expected_cells
            with self.subTest(expected_cells=expected_cells), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "cell cardinality differs"):
                subject._development_timing_shape(candidate, panel_count)
        for candidate_panels in (
                subject.FROZEN_TIMING_PANEL_COUNT - 1,
                subject.FROZEN_TIMING_PANEL_COUNT + 1):
            with self.subTest(panel_count=candidate_panels), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "panel cardinality differs"):
                subject._development_timing_shape(contract, candidate_panels)

    def test_timing_worker_count_rejects_before_native_dispatch(self) -> None:
        contract = self._timing_contract(16)
        qualification = self._qualification(contract)
        for count in (0, 7, 9):
            with self.subTest(count=count), self.assertRaisesRegex(
                    subject.RunnerError, "exactly eight timing workers"):
                subject._run_timing_jobs(
                    contract, {}, qualification, {}, [object()] * count,
                    self.root, 0, time.monotonic() + 1.0)
        for repetitions in (0, 7, 12):
            candidate = self._timing_contract(16)
            domain = candidate["timing"]["domains"]["development"]
            domain["paired_repetitions"] = repetitions
            with self.subTest(repetitions=repetitions), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "positive multiple of eight"):
                subject._timing_job_waves(
                    candidate, 11, qualification)

    def test_sampler_fixture_waits_for_an_exact_new_sample(self) -> None:
        path = self.root / "thermal.csv"
        baseline = 1000.0
        with path.open("w", encoding="ascii", newline="") as output:
            writer = csv.writer(output, lineterminator="\n")
            writer.writerow(native_api.THERMAL_HEADER)
            writer.writerow(sampler_row(baseline))

        def append() -> None:
            time.sleep(0.05)
            with path.open("a", encoding="ascii", newline="") as output:
                csv.writer(output, lineterminator="\n").writerow(
                    sampler_row(baseline + 1.0))

        thread = threading.Thread(target=append)
        thread.start()
        try:
            selected = subject.choose_new_sampler_start(
                path, time.monotonic() + 2.0)
        finally:
            thread.join()
        self.assertEqual(selected, 1001000000000)
        self.assertEqual(subject._wait_for_sampler_sample(
            path, time.monotonic() + 1.0,
            at_or_after_ns=1000500000000), selected)

    def test_sampler_never_selects_an_unterminated_endpoint(self) -> None:
        path = self.root / "partial-thermal.csv"
        header = ",".join(native_api.THERMAL_HEADER)
        first = ",".join(sampler_row(1000.0))
        partial = ",".join(sampler_row(1001.0))
        path.write_bytes((header + "\n" + first + "\n" + partial).encode(
            "ascii"))
        self.assertEqual(subject._valid_sampler_samples(path), [1000000000000])
        with path.open("ab") as output:
            output.write(b"\n")
        self.assertEqual(subject._valid_sampler_samples(path), [
            1000000000000, 1001000000000])

    def test_root_owned_sampler_denial_requests_sudo_without_weakening(self) \
            -> None:
        fake_stat = type("FakeStat", (), {"st_uid": 0})()
        script = self.root / "sampler.py"
        thermal = self.root / "thermal.csv"
        script.write_text("pass\n")
        thermal.write_text("fixture\n")
        with mock.patch.object(
                native_api, "_parse_proc_start_ticks", return_value=1), \
                mock.patch.object(
                    native_api, "_verify_live_sampler_process",
                    side_effect=native_api.NativeEvidenceError("access denied")), \
                mock.patch.object(Path, "stat", return_value=fake_stat), \
                mock.patch.object(os, "geteuid", return_value=1000):
            with self.assertRaisesRegex(
                    subject.RunnerError,
                    r"sudo -n python3\.12; validation is not weakened"):
                subject._preflight_sampler(
                    3320493, 127, script, thermal)

    def test_persistent_workers_obey_barriers_and_exact_q(self) -> None:
        worker_path = self._fake_worker()
        description = self._description(worker_path)
        cpus = self._worker_cpus()
        log = self.root / "worker.log"
        old = os.environ.get("WH2_FAKE_WORKER_LOG")
        os.environ["WH2_FAKE_WORKER_LOG"] = str(log)
        workers = []
        sink = subject.AtomicLineSink(self.root / "native.jsonl")
        try:
            workers = subject.spawn_workers(
                description, cpus, time.monotonic() + 5.0)
            first = [subject.Job("recovery", value, 0, value)
                     for value in range(6)]
            second = [subject.Job("timing", value, 0, value)
                      for value in range(6)]
            _, used0 = subject.run_job_batch(
                workers, first, 0, sink, time.monotonic() + 5.0,
                self._validator)
            _, used1 = subject.run_job_batch(
                workers, second, 1, sink, time.monotonic() + 5.0,
                self._validator)
            self.assertEqual(used0, set(cpus))
            self.assertEqual(used1, set(cpus))
            sink.publish()
            subject.quit_workers(workers, time.monotonic() + 5.0)
            workers = []
        finally:
            sink.abort()
            subject.terminate_workers(workers)
            if old is None:
                os.environ.pop("WH2_FAKE_WORKER_LOG", None)
            else:
                os.environ["WH2_FAKE_WORKER_LOG"] = old
        events = [json.loads(line) for line in log.read_text().splitlines()]
        last_recovery_done = max(
            index for index, event in enumerate(events)
            if event["state"] == "done" and event["command"].startswith("R "))
        first_timing_start = min(
            index for index, event in enumerate(events)
            if event["state"] == "start" and event["command"].startswith("T "))
        self.assertLess(last_recovery_done, first_timing_start)
        self.assertEqual(len((self.root / "native.jsonl").read_text().splitlines()),
                         12)

    def test_job_batch_rejects_exact_raw_cap_before_json_parse(self) -> None:
        process = mock.Mock()
        process.poll.return_value = None
        process.stdout.fileno.return_value = 17
        worker = subject.PersistentWorker(3, process, 1)
        job = subject.Job("recovery", 1, 0, 0)
        selector = mock.Mock()
        selector.select.return_value = [
            (mock.Mock(data=worker), subject.selectors.EVENT_READ),
        ]
        parser = mock.Mock(side_effect=AssertionError(
            "oversized raw response reached JSON parser"))
        validator = mock.Mock()
        sink = mock.Mock()
        line = b'{"value":"oversized"}\n'
        self.assertGreater(len(line), 16)

        with mock.patch.object(
                subject.selectors, "DefaultSelector",
                return_value=selector), mock.patch.object(
                subject, "_write_worker"), mock.patch.object(
                subject.os, "read", return_value=line), mock.patch.object(
                subject, "_parse_canonical_line", parser), \
                self.assertRaisesRegex(
                    subject.RunnerError, "exact byte cap"):
            subject.run_job_batch(
                [worker], [job], 0, sink, time.monotonic() + 1.0,
                validator, maximum_response_bytes=16)
        parser.assert_not_called()
        validator.assert_not_called()
        sink.write.assert_not_called()
        selector.close.assert_called_once_with()

    def test_job_batch_rejects_invalid_exact_raw_caps(self) -> None:
        for invalid in (
                None, False, 0, -1, 1.0, subject.MAX_RESPONSE_BYTES + 1):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(
                    subject.RunnerError, "response byte cap"):
                subject.run_job_batch(
                    [], [], 0, mock.Mock(), time.monotonic() + 1.0,
                    mock.Mock(), maximum_response_bytes=invalid)

    def test_timing_runner_dispatches_only_t_commands_and_no_recovery_file(
            self) -> None:
        cpus = sorted(os.sched_getaffinity(0))[:subject.TIMING_WORKER_COUNT]
        if len(cpus) != subject.TIMING_WORKER_COUNT:
            self.skipTest("eight logical CPUs are required")
        worker_path = self._fake_worker()
        description = self._description(worker_path)
        log = self.root / "timing-worker.log"
        old = os.environ.get("WH2_FAKE_WORKER_LOG")
        os.environ["WH2_FAKE_WORKER_LOG"] = str(log)
        workers = []
        jobs = [
            subject.Job("timing", ordinal, 0, ordinal)
            for ordinal in range(subject.TIMING_WORKER_COUNT)
        ]
        try:
            workers = subject.spawn_workers(
                description, cpus, time.monotonic() + 5.0)
            sibling_guard = mock.Mock()
            with mock.patch.object(
                    contract_api, "timing_panels", return_value=[]), \
                    mock.patch.object(
                        subject, "_timing_job_waves",
                        return_value=[(0, jobs)]), \
                    mock.patch.object(
                        subject, "_strict_response_validator",
                        return_value=self._validator):
                path, _ = subject._run_timing_jobs(
                    {}, {}, mock.Mock(), description, workers,
                    self.root, 0, time.monotonic() + 5.0, sibling_guard)
            subject.quit_workers(workers, time.monotonic() + 5.0)
            workers = []
        finally:
            subject.terminate_workers(workers)
            if old is None:
                os.environ.pop("WH2_FAKE_WORKER_LOG", None)
            else:
                os.environ["WH2_FAKE_WORKER_LOG"] = old
        commands = [json.loads(line)["command"]
                    for line in log.read_text().splitlines()]
        self.assertTrue(commands)
        self.assertTrue(all(command.startswith("T ") for command in commands))
        self.assertEqual(path, self.root / "timing-native-results.jsonl")
        self.assertFalse(
            (self.root / "recovery-native-results.jsonl").exists())
        self.assertEqual(
            [call.args[0] for call in sibling_guard.check.call_args_list],
            ["before-timing-wave-0", "after-timing-wave-0"])

    def test_qualification_adaptive_tail_retries_only_one_unresolved_cell(
            self) -> None:
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < 9:
            self.skipTest("more than eight logical CPUs are required")
        cpus = allowed[:min(12, len(allowed))]
        contract = self._timing_contract(8)
        expected_cells = contract["timing"]["domains"]["development"][
            "expected_cells"]
        target_cell = expected_cells - 1
        worker_path = self._fake_worker()
        description = self._description(worker_path)
        log = self.root / "qualification-worker.log"
        output = self.root / "qualification-native.jsonl"
        old_log = os.environ.get("WH2_FAKE_WORKER_LOG")
        os.environ["WH2_FAKE_WORKER_LOG"] = str(log)
        workers = []
        observed = []

        def validate(_contract, _description, _cells, _window_start,
                     _value, _worker, job):
            observed.append((job.cell, job.retry_offset))
            success = not (
                job.cell == target_cell and job.retry_offset == 0)
            return time.monotonic_ns(), success

        try:
            workers = subject.spawn_workers(
                description, cpus, time.monotonic() + 10.0)
            with mock.patch.object(
                    subject, "_strict_qualification_response",
                    side_effect=validate):
                path, maximum_end = subject.run_timing_qualification(
                    contract, description, workers, output, 0,
                    time.monotonic() + 10.0)
            self.assertEqual(path, output)
            self.assertGreater(maximum_end, 0)
            subject.quit_workers(workers, time.monotonic() + 10.0)
            workers = []
        finally:
            subject.terminate_workers(workers)
            if old_log is None:
                os.environ.pop("WH2_FAKE_WORKER_LOG", None)
            else:
                os.environ["WH2_FAKE_WORKER_LOG"] = old_log
        self.assertEqual(len(observed), expected_cells + 1)
        self.assertEqual(
            sorted(value for value in observed if value[0] == target_cell),
            [(target_cell, 0), (target_cell, 1)])
        self.assertTrue(all(
            retry == 0 for cell, retry in observed if cell != target_cell))
        self.assertEqual(
            len(output.read_text(encoding="utf-8").splitlines()),
            expected_cells + 1)
        commands = [
            json.loads(line)["command"]
            for line in log.read_text(encoding="utf-8").splitlines()
            if json.loads(line)["state"] == "start"
        ]
        self.assertEqual(commands.count(
            "L {} 0".format(target_cell)), 1)
        self.assertEqual(commands.count(
            "L {} 1".format(target_cell)), 1)
        self.assertFalse(any(command.endswith(" 2")
                             for command in commands))

    def test_qualification_selector_register_failure_closes_and_aborts(
            self) -> None:
        contract = {
            "timing": {"domains": {"development": {"expected_cells": 1}}},
        }
        process = mock.Mock()
        process.stdout = mock.Mock()
        worker = mock.Mock(
            cpu=3, process=process, pending=None, buffer=b"")
        sink = mock.Mock()
        selector = mock.Mock()
        selector.register.side_effect = OSError(
            "injected selector registration failure")
        with mock.patch.object(
                contract_api, "iter_timing_base_cells",
                return_value=iter(({},))), \
                mock.patch.object(
                    subject, "AtomicLineSink", return_value=sink), \
                mock.patch.object(
                    subject.selectors, "DefaultSelector",
                    return_value=selector), \
                self.assertRaisesRegex(
                    OSError, "selector registration failure"):
            subject.run_timing_qualification(
                contract, {}, [worker], self.root / "qualification.jsonl",
                0, time.monotonic() + 1.0)
        selector.close.assert_called_once_with()
        sink.abort.assert_called_once_with()
        sink.publish.assert_not_called()

    def test_timing_sink_constructor_failure_has_no_recovery_side_effect(self) \
            -> None:
        construction_paths = []

        def construct(path, *_args):
            construction_paths.append(path)
            raise OSError("injected timing sink construction failure")

        workers = [mock.Mock(cpu=cpu) for cpu in range(8)]
        qualification = mock.Mock()
        with mock.patch.object(
                contract_api, "timing_panels", return_value=[]), \
                mock.patch.object(
                    subject, "_timing_job_waves", return_value=[]), \
                mock.patch.object(
                    subject, "AtomicLineSink", side_effect=construct), \
                self.assertRaisesRegex(
                    OSError, "timing sink construction failure"):
            subject._run_timing_jobs(
                {}, {}, qualification, {}, workers, self.root, 0,
                time.monotonic() + 1.0)
        self.assertEqual(
            construction_paths,
            [self.root / "timing-native-results.jsonl"])

    def test_worker_error_terminates_complete_pool(self) -> None:
        worker_path = self._fake_worker()
        description = self._description(worker_path)
        cpus = self._worker_cpus()
        old = os.environ.get("WH2_FAKE_WORKER_FAIL")
        os.environ["WH2_FAKE_WORKER_FAIL"] = "R 1 0"
        workers = []
        sink = subject.AtomicLineSink(self.root / "failed.jsonl")
        try:
            workers = subject.spawn_workers(
                description, cpus, time.monotonic() + 5.0)
            jobs = [subject.Job("recovery", value, 0, value)
                    for value in range(6)]
            with self.assertRaises(subject.RunnerError):
                subject.run_job_batch(
                    workers, jobs, 0, sink, time.monotonic() + 5.0,
                    self._validator)
        finally:
            sink.abort()
            subject.terminate_workers(workers)
            if old is None:
                os.environ.pop("WH2_FAKE_WORKER_FAIL", None)
            else:
                os.environ["WH2_FAKE_WORKER_FAIL"] = old
        self.assertTrue(workers)
        self.assertTrue(all(worker.process.poll() is not None
                            for worker in workers))
        self.assertFalse((self.root / "failed.jsonl").exists())

    def test_terminate_workers_reports_stubborn_survivor_after_full_cleanup(
            self) -> None:
        class StubbornProcess:
            def __init__(self) -> None:
                self.stdin = mock.Mock()
                self.stdout = mock.Mock()
                self.stderr = mock.Mock()
                self.terminate_calls = 0
                self.kill_calls = 0
                self.wait_timeouts: List[float] = []

            def poll(self):
                return None

            def terminate(self) -> None:
                self.terminate_calls += 1

            def kill(self) -> None:
                self.kill_calls += 1

            def wait(self, timeout=None):
                self.wait_timeouts.append(timeout)
                raise subprocess.TimeoutExpired("stubborn-worker", timeout)

        process = StubbornProcess()
        worker = subject.PersistentWorker(17, process, 1)
        with self.assertRaisesRegex(
                subject.RunnerError, "unreaped CPUs 17"):
            subject.terminate_workers([worker])
        self.assertEqual(process.terminate_calls, 1)
        self.assertEqual(process.kill_calls, 1)
        self.assertEqual(len(process.wait_timeouts), 2)
        for stream in (process.stdin, process.stdout, process.stderr):
            stream.close.assert_called_once_with()

    def test_terminate_workers_defers_first_failure_and_cleans_later_worker(
            self) -> None:
        class PollInterrupt(BaseException):
            pass

        poll_interrupt = PollInterrupt("injected first-worker poll interrupt")

        class Process:
            def __init__(self, first: bool) -> None:
                self.stdin = mock.Mock()
                self.stdout = mock.Mock()
                self.stderr = mock.Mock()
                self.first = first
                self.status = None
                self.poll_calls = 0
                self.terminate_calls = 0
                self.kill_calls = 0
                self.wait_calls = 0

            def poll(self):
                self.poll_calls += 1
                if self.first and self.poll_calls == 1:
                    raise poll_interrupt
                return self.status

            def terminate(self) -> None:
                self.terminate_calls += 1

            def kill(self) -> None:
                self.kill_calls += 1

            def wait(self, timeout=None):
                self.wait_calls += 1
                if self.wait_calls == 1:
                    if self.first:
                        raise OSError(
                            errno.EIO, "injected first-worker wait failure")
                    raise subprocess.TimeoutExpired(
                        "injected later stubborn worker", timeout)
                self.status = -9
                return self.status

        first = Process(True)
        later = Process(False)
        workers = (
            subject.PersistentWorker(11, first, 1),
            subject.PersistentWorker(13, later, 1),
        )
        with self.assertRaises(PollInterrupt) as raised:
            subject.terminate_workers(workers)
        self.assertIs(raised.exception, poll_interrupt)
        self.assertIsInstance(
            raised.exception.__cause__, subject.RunnerError)
        diagnostic = str(raised.exception.__cause__)
        self.assertIn("first-worker poll interrupt", diagnostic)
        self.assertIn("first-worker wait failure", diagnostic)
        for process in (first, later):
            self.assertEqual(process.terminate_calls, 1)
            self.assertEqual(process.kill_calls, 1)
            self.assertEqual(process.wait_calls, 2)
            self.assertEqual(process.status, -9)
            for stream in (process.stdin, process.stdout, process.stderr):
                stream.close.assert_called_once_with()

    def test_git_head_rejects_relevant_dirty_or_untracked_source(self) -> None:
        repo = self.root / "repo"
        repo.mkdir()
        for directory in ("codec", "bench", "include"):
            (repo / directory).mkdir()
        (repo / "CMakeLists.txt").write_text("cmake_minimum_required(VERSION 3.10)\n")
        (repo / "codec" / "codec.cpp").write_text("int codec = 1;\n")
        (repo / "README.md").write_text("fixture\n")
        subprocess.run(["git", "init", "-q"], cwd=str(repo), check=True)
        subprocess.run(["git", "config", "user.email", "test@example.com"],
                       cwd=str(repo), check=True)
        subprocess.run(["git", "config", "user.name", "Test"],
                       cwd=str(repo), check=True)
        subprocess.run(["git", "add", "CMakeLists.txt", "codec/codec.cpp",
                        "README.md"],
                       cwd=str(repo), check=True)
        subprocess.run(["git", "commit", "-qm", "fixture"], cwd=str(repo),
                       check=True)
        head = subject._git_head(time.monotonic() + 5.0, repo)
        self.assertRegex(head, r"^[0-9a-f]{40}$")
        # Unrelated experiment/result artifacts do not falsify source identity.
        (repo / "preferred-route.csv").write_text("user artifact\n")
        self.assertEqual(subject._git_head(time.monotonic() + 5.0, repo), head)
        (repo / "bench" / "untracked.py").write_text("pass\n")
        with self.assertRaises(subject.RunnerError):
            subject._git_head(time.monotonic() + 5.0, repo)
        (repo / "bench" / "untracked.py").unlink()
        (repo / "new_codec.h").write_text("int untracked_header;\n")
        with self.assertRaises(subject.RunnerError):
            subject._git_head(time.monotonic() + 5.0, repo)
        (repo / "new_codec.h").unlink()
        (repo / "README.md").write_text("tracked unrelated change\n")
        with self.assertRaises(subject.RunnerError):
            subject._git_head(time.monotonic() + 5.0, repo)
        subprocess.run(["git", "checkout", "--", "README.md"], cwd=str(repo),
                       check=True)
        (repo / "codec" / "codec.cpp").write_text("int codec = 2;\n")
        with self.assertRaises(subject.RunnerError):
            subject._git_head(time.monotonic() + 5.0, repo)


if __name__ == "__main__":
    unittest.main()
