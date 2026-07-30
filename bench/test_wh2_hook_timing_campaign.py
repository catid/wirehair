#!/usr/bin/env python3
"""Tests for the paired WH2 hook-overhead campaign."""

import argparse
import ast
import contextlib
import dataclasses
import fcntl
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import signal
import stat
import sys
import tempfile
import textwrap
import threading
import time
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
campaign = importlib.import_module("wh2_hook_timing_campaign")


MOCK_TIMING = r'''#!/usr/bin/python3
import argparse
import ctypes
import hashlib
import json
import os
import signal
import sys
import time

p = argparse.ArgumentParser()
p.add_argument("--K", type=int, required=True)
p.add_argument("--bb", type=int, required=True)
p.add_argument("--construction-seed", type=int, required=True)
p.add_argument("--loss-seed", type=int, required=True)
p.add_argument("--loss-ppm", type=int, required=True)
p.add_argument("--schedule", required=True)
p.add_argument("--scope", required=True)
p.add_argument("--warmup-reps", type=int, required=True)
p.add_argument("--inner-reps", type=int, required=True)
p.add_argument("--max-working-mib", type=int, required=True)
p.add_argument("--context-sha256", required=True)
p.add_argument("--ready-fd", type=int)
p.add_argument("--go-fd", type=int)
a = p.parse_args()
mode = os.environ.get("MOCK_MODE", "normal")
actual_hooks = __ACTUAL_HOOKS__
hooks = (
    1 - actual_hooks
    if mode == "wrong-hooks" or (
        mode == "wrong-hooks-on" and actual_hooks == 1)
    else actual_hooks
)
cpu = ctypes.CDLL(None).sched_getcpu()

if mode == "exit":
    sys.stderr.write("mock failure\n")
    raise SystemExit(1)
if mode == "exit-on" and actual_hooks == 1:
    raise SystemExit(7)
if mode == "leader-exit-descendant" and actual_hooks == 1:
    child = os.fork()
    if child == 0:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        time.sleep(30)
        os._exit(0)
    pid_dir = os.environ.get("MOCK_PID_DIR")
    if pid_dir:
        with open(os.path.join(pid_dir, "descendant"), "w") as f:
            f.write(str(child))
    os._exit(7)
if mode == "stderr":
    sys.stderr.write("surprise\n")
if mode in ("sleep", "ignore-term-sleep"):
    if mode == "ignore-term-sleep":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    child = os.fork()
    if child == 0:
        time.sleep(30)
        os._exit(0)
    pid_dir = os.environ.get("MOCK_PID_DIR")
    if pid_dir:
        with open(os.path.join(pid_dir, str(os.getpid())), "w") as f:
            f.write(str(child))
    time.sleep(30)

def h(tag):
    material = "|".join(map(str, (
        tag, a.K, a.bb, a.construction_seed, a.loss_seed,
        a.loss_ppm, a.schedule,
        hooks if mode == "semantic-drift" else 0,
    )))
    return hashlib.sha256(material.encode("ascii")).hexdigest()

weak = mode == "weak" or (
    mode == "weak-disagreement" and actual_hooks == 1)
weak_result = int(os.environ.get("MOCK_WEAK_RESULT", "4"))
status = "weak-root" if weak else "success"
header = {
    "record": "header",
    "schema": "wirehair.wh2.hook-timing.v1",
    "profile": "dispatch-v1",
    "hooks_compiled": hooks,
    "K": a.K,
    "block_bytes": a.bb,
    "construction_seed": a.construction_seed,
    "loss_seed": a.loss_seed,
    "loss_ppm": a.loss_ppm,
    "schedule": a.schedule,
    "scope_request": a.scope,
    "warmup_reps": a.warmup_reps,
    "inner_reps": a.inner_reps,
    "max_working_mib": a.max_working_mib,
    "context_sha256": a.context_sha256,
    "cpu": cpu,
    "clock": "CLOCK_MONOTONIC",
    "start_barrier":
        "none" if a.ready_fd is None else "ready-go-pipe-v1",
}
semantic = {
    "record": "semantic",
    "schema": "wirehair.wh2.hook-timing.v1",
    "profile": "dispatch-v1",
    "status": status,
    "encoder_result": weak_result if weak else 0,
    "decoder_result": None if weak else 0,
    "recover_result": None if weak else 0,
    "direct_result": None if weak else 0,
    "trace_packets": 2 * a.K + 512,
    "delivered_packets": 0 if weak else a.K,
    "overhead_packets": None if weak else 0,
    "decoder_solve_attempts": 0 if weak else 1,
    "profile_sha256": h("profile"),
    "system_sha256": h("system"),
    "coefficients_sha256": h("coefficients"),
    "trace_sha256": h("trace"),
    "row_sha256": h("row"),
    "message_sha256": h("message"),
    "payload_sha256": "not_applicable" if weak else h("payload"),
    "intermediate_sha256": "not_applicable" if weak else h("intermediate"),
    "recovered_sha256": "not_applicable" if weak else h("recovered"),
    "encoder_stats_sha256": "not_applicable" if weak else h("encoder_stats"),
    "decoder_stats_sha256": "not_applicable" if weak else h("decoder_stats"),
    "direct_stats_sha256": "not_applicable" if weak else h("direct_stats"),
    "semantic_sha256": h("semantic"),
}
records = [header, semantic]
if a.scope != "semantic":
    if a.ready_fd is None or a.go_fd is None:
        raise SystemExit(3)
    ready = b"X" if mode == "bad-ready" else (
        b"RR" if mode == "extra-ready" else b"R")
    os.write(a.ready_fd, ready)
    os.close(a.ready_fd)
    go = os.read(a.go_fd, 2)
    extra = os.read(a.go_fd, 1)
    os.close(a.go_fd)
    if go != b"G" or extra != b"":
        raise SystemExit(3)
    if mode == "term-proof":
        signal.signal(signal.SIGTERM, lambda unused_s, unused_f: os._exit(19))
        pid_dir = os.environ.get("MOCK_PID_DIR")
        if pid_dir:
            with open(
                    os.path.join(pid_dir, str(actual_hooks)), "w") as f:
                f.write(str(os.getpid()))
        time.sleep(30)
if not weak and a.scope != "semantic":
    lifecycle = {
        "row": "caller-row-buffer-reuse-v1",
        "encoder":
            "fresh-first-then-transactional-reinitialize-including-prior-release-v1",
        "decoder-feed": "distinct-preinitialized-endpoints-v1",
        "decoder-full":
            "fresh-first-then-transactional-reinitialize-including-prior-release-v1",
        "direct": "fresh-first-then-transactional-output-reuse-v1",
    }[a.scope]
    scale = 1000000 if mode == "slow" else 100000000
    elapsed = a.inner_reps * scale
    records.append({
        "record": "timing",
        "schema": "wirehair.wh2.hook-timing.v1",
        "profile": "dispatch-v1",
        "hooks_compiled": hooks,
        "scope": a.scope,
        "lifecycle": lifecycle,
        "semantic_sha256": semantic["semantic_sha256"],
        "unit": "ns",
        "clock": "CLOCK_MONOTONIC",
        "warmup_reps": a.warmup_reps,
        "measured_reps": 1,
        "inner_reps": a.inner_reps,
        "work_items_per_rep": a.K,
        "elapsed_ns": elapsed,
        "min_ns": elapsed,
        "max_ns": elapsed,
        "minor_faults": actual_hooks,
        "major_faults": 0,
        "voluntary_context_switches": 1,
        "involuntary_context_switches": 2,
        "sink": 123,
        "result": 0,
    })
lines = [
    json.dumps(row, separators=(",", ":"), ensure_ascii=True) + "\n"
    for row in records
]
prefix = "".join(lines)
done_hash = hashlib.sha256(prefix.encode("ascii")).hexdigest()
if mode == "hash-mismatch":
    done_hash = "0" * 64
done = {
    "record": "done",
    "schema": "wirehair.wh2.hook-timing.v1",
    "status": status,
    "records_before_done": len(records),
    "stream_sha256": done_hash,
}
lines.append(json.dumps(done, separators=(",", ":")) + "\n")
output = "".join(lines)
if mode == "truncated":
    output = "".join(lines[:-1])
elif mode == "malformed":
    output = output.replace(
        '"record":"header"', '"record":"header","record":"header"', 1)
sys.stdout.write(output)
'''


def _write_executable(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)


def _mutate_authenticated_stream(payload, record_index, key, value):
    records = [
        json.loads(line)
        for line in payload.decode("ascii").splitlines()]
    records[record_index][key] = value
    lines = [
        json.dumps(row, separators=(",", ":"), ensure_ascii=True) + "\n"
        for row in records[:-1]
    ]
    records[-1]["stream_sha256"] = hashlib.sha256(
        "".join(lines).encode("ascii")).hexdigest()
    lines.append(json.dumps(
        records[-1], separators=(",", ":"), ensure_ascii=True) + "\n")
    return "".join(lines).encode("ascii")


def _sibling_pair():
    if not hasattr(os, "sched_getaffinity"):
        return None
    allowed = set(os.sched_getaffinity(0))
    root = Path("/sys/devices/system/cpu")
    for cpu in sorted(allowed):
        try:
            text = (
                root / "cpu{}".format(cpu) /
                "topology/thread_siblings_list").read_text().strip()
            siblings = campaign._parse_cpu_list(text)
        except (OSError, campaign.CampaignError):
            continue
        for sibling in sorted(siblings):
            if sibling != cpu and sibling in allowed:
                return cpu, sibling
    return None


class MockFixture:
    def __init__(self, root):
        self.root = Path(root)
        self.paths = {}
        self.bindings = {}
        for treatment in campaign.TREATMENTS:
            path = self.root / treatment / treatment
            hooks = campaign.TREATMENT_STATE[treatment]
            _write_executable(
                path, MOCK_TIMING.replace("__ACTUAL_HOOKS__", str(hooks)))
            current = path.stat()
            payload = path.read_bytes()
            mode = stat.S_IMODE(current.st_mode)
            self.paths[treatment] = path
            self.bindings[treatment] = campaign.BinaryBinding(
                treatment=treatment,
                hook_state=hooks,
                replicate=treatment[-1],
                path=str(path),
                bytes=len(payload),
                sha256=hashlib.sha256(payload).hexdigest(),
                build_id="11112222" if "hook-on" in treatment
                else "33334444",
                device=current.st_dev,
                inode=current.st_ino,
                mode=mode,
                execution_fd=campaign._sealed_executable_memfd(
                    "mock-{}".format(treatment), payload),
            )

    def close(self):
        campaign.close_execution_bindings(self.bindings)

    def expected(self, treatment, cpu, scope="row", inner=1):
        return campaign.InvocationExpected(
            treatment=treatment,
            K=64,
            block_bytes=2,
            construction_seed=1,
            loss_seed=2,
            loss_ppm=350000,
            schedule="repair-only",
            scope=scope,
            inner_reps=inner,
            max_working_mib=1024,
            context_sha256="0" * 64,
            cpu=cpu,
        )

    def specs(self, cpus, scope="row", inner=1):
        output = []
        for index, treatment in enumerate(
                ("hook-on-a", "hook-off-a")):
            expected = self.expected(
                treatment, cpus[index], scope=scope, inner=inner)
            output.append(campaign.InvocationSpec(
                ordinal=index,
                pair_ordinal=0,
                phase="test",
                case_id="test",
                block=None,
                round_index=None,
                treatment=treatment,
                cpu=cpus[index],
                executable_fd=self.bindings[treatment].execution_fd,
                expected=expected,
                command=campaign.invocation_command(
                    self.bindings[treatment], expected),
            ))
        return tuple(output)


class StrictJsonTests(unittest.TestCase):
    def test_duplicate_key_noncanonical_and_done_hash_rejected(self):
        with self.assertRaisesRegex(campaign.CampaignError, "duplicate"):
            campaign.strict_json_line(b'{"a":1,"a":2}\n', "duplicate")
        with self.assertRaisesRegex(campaign.CampaignError, "compact"):
            campaign.strict_json_line(b'{"a": 1}\n', "spaced")

    def test_source_contains_no_duplicate_dict_literal_keys(self):
        tree = ast.parse(
            Path(campaign.__file__).read_text(encoding="utf-8"))
        duplicates = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Dict):
                continue
            keys = []
            for key in node.keys:
                if isinstance(key, ast.Constant) and isinstance(key.value, str):
                    keys.append(key.value)
            repeated = sorted({key for key in keys if keys.count(key) > 1})
            if repeated:
                duplicates.append((node.lineno, repeated))
        self.assertEqual(duplicates, [])

    def test_integer_and_nesting_work_are_bounded(self):
        with self.assertRaisesRegex(campaign.CampaignError, "integer"):
            campaign.strict_json_document(
                ("{" + '"n":' + "9" * 5000 + "}").encode("ascii"), "huge")
        with self.assertRaisesRegex(campaign.CampaignError, "nesting"):
            campaign.strict_json_document(
                ("[" * 10000 + "]" * 10000).encode("ascii"), "deep")
        quoted = {
            "brackets": "[[[{{{",
            "escaped_quote": '\\"still-in-string[[[',
        }
        payload = campaign.canonical_json_bytes(quoted)
        self.assertEqual(
            campaign.strict_json_document(payload, "quoted"), quoted)


@unittest.skipUnless(_sibling_pair(), "requires allowed SMT siblings")
class StreamAndProcessTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = MockFixture(self.temporary.name)
        self.cpus = _sibling_pair()

    def tearDown(self):
        self.fixture.close()
        self.temporary.cleanup()

    def _captures(self, mode="normal", scope="row", inner=1,
                  timeout=5):
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": mode,
        }
        specs = self.fixture.specs(self.cpus, scope=scope, inner=inner)
        captures = campaign.run_simultaneous_pair(
            specs, timeout, environment=environment)
        return specs, captures

    def test_valid_success_and_semantic_streams(self):
        for scope in ("row", "semantic"):
            specs, captures = self._captures(scope=scope)
            parsed = [
                campaign.parse_timing_stream(
                    capture.stdout, spec.expected)
                for spec, capture in zip(specs, captures)
            ]
            self.assertTrue(all(capture.returncode == 0 for capture in captures))
            self.assertTrue(all(not capture.stderr for capture in captures))
            self.assertEqual({row.status for row in parsed}, {"success"})
            self.assertEqual(
                [row.timing is None for row in parsed],
                [scope == "semantic", scope == "semantic"])
            self.assertEqual(parsed[0].semantic, parsed[1].semantic)

    def test_truncation_malformed_and_done_hash_mismatch(self):
        for mode, error in (
                ("truncated", "truncated"),
                ("malformed", "duplicate"),
                ("hash-mismatch", "done")):
            with self.subTest(mode=mode):
                specs, captures = self._captures(mode=mode)
                with self.assertRaisesRegex(
                        campaign.CampaignError,
                        "done" if mode == "truncated" else error):
                    campaign.parse_timing_stream(
                        captures[0].stdout, specs[0].expected)

    def test_wrong_hooks_rejected(self):
        specs, captures = self._captures(mode="wrong-hooks")
        with self.assertRaisesRegex(campaign.CampaignError, "header"):
            campaign.parse_timing_stream(
                captures[0].stdout, specs[0].expected)

    def test_hash_valid_boolean_numeric_aliases_are_rejected(self):
        specs, captures = self._captures()
        for record_index, key, value in (
                (0, "hooks_compiled", True),
                (1, "encoder_result", False),
                (2, "measured_reps", True),
                (-1, "records_before_done", True)):
            with self.subTest(record=record_index, key=key):
                mutated = _mutate_authenticated_stream(
                    captures[0].stdout, record_index, key, value)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "integer"):
                    campaign.parse_timing_stream(
                        mutated, specs[0].expected)

    def test_ready_barrier_rejects_bad_extra_and_early_exit(self):
        for mode, error in (
                ("bad-ready", "ready"),
                ("extra-ready", "ready"),
                ("exit-on", "start barrier|ready pipe")):
            with self.subTest(mode=mode):
                started = time.monotonic()
                with self.assertRaisesRegex(campaign.CampaignError, error):
                    self._captures(mode=mode, timeout=2)
                self.assertLess(time.monotonic() - started, 1.5)

    def test_invocation_spec_is_immutable_and_barrier_args_do_not_accumulate(self):
        specs = self.fixture.specs(self.cpus)
        original_commands = tuple(spec.command for spec in specs)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            specs[0].command = ("changed",)
        first = campaign.run_simultaneous_pair(
            specs, 5, environment={
                "LC_ALL": "C", "LANG": "C", "TZ": "UTC"})
        second = campaign.run_simultaneous_pair(
            specs, 5, environment={
                "LC_ALL": "C", "LANG": "C", "TZ": "UTC"})
        self.assertEqual(
            tuple(spec.command for spec in specs), original_commands)
        for captures in (first, second):
            for capture in captures:
                self.assertEqual(capture.launched_command.count("--ready-fd"), 1)
                self.assertEqual(capture.launched_command.count("--go-fd"), 1)

    def test_pipe_and_selector_construction_failures_restore_resources(self):
        specs = self.fixture.specs(self.cpus)
        baseline_fds = set(os.listdir("/proc/self/fd"))
        original_mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        with mock.patch.object(
                campaign.selectors, "DefaultSelector",
                side_effect=RuntimeError("selector injection")):
            with self.assertRaisesRegex(RuntimeError, "selector injection"):
                campaign.run_simultaneous_pair(specs, 5)
        self.assertEqual(
            signal.pthread_sigmask(signal.SIG_BLOCK, set()), original_mask)
        self.assertEqual(set(os.listdir("/proc/self/fd")), baseline_fds)

        real_pipe2 = campaign.os.pipe2
        calls = []

        def fail_second_pipe(flags):
            calls.append(flags)
            if len(calls) == 2:
                raise OSError("pipe injection")
            return real_pipe2(flags)

        with mock.patch.object(campaign.os, "pipe2", side_effect=fail_second_pipe):
            with self.assertRaisesRegex(OSError, "pipe injection"):
                campaign.run_simultaneous_pair(specs, 5)
        self.assertEqual(set(os.listdir("/proc/self/fd")), baseline_fds)

    def test_second_spawn_failure_reaps_first_and_leaks_no_fds(self):
        specs = self.fixture.specs(self.cpus)
        baseline_fds = set(os.listdir("/proc/self/fd"))
        real_popen = campaign.subprocess.Popen
        spawned = []

        def fail_second_spawn(*args, **kwargs):
            if spawned:
                raise OSError("spawn injection")
            child = real_popen(*args, **kwargs)
            spawned.append(child.pid)
            return child

        with mock.patch.object(
                campaign.subprocess, "Popen", side_effect=fail_second_spawn):
            with self.assertRaisesRegex(OSError, "spawn injection"):
                campaign.run_simultaneous_pair(specs, 5)
        self.assertEqual(len(spawned), 1)
        self.assertFalse(campaign._group_exists(spawned[0]))
        self.assertEqual(set(os.listdir("/proc/self/fd")), baseline_fds)

    def test_spawn_signal_race_reaps_tracked_child(self):
        specs = self.fixture.specs(self.cpus)
        real_popen = campaign.subprocess.Popen
        spawned = []

        def interrupt_after_spawn(*args, **kwargs):
            child = real_popen(*args, **kwargs)
            spawned.append(child.pid)
            signal.raise_signal(signal.SIGTERM)
            return child

        previous = campaign._install_signal_handlers()
        try:
            with mock.patch.object(
                    campaign.subprocess, "Popen",
                    side_effect=interrupt_after_spawn):
                with self.assertRaises(campaign.CampaignInterrupted):
                    campaign.run_simultaneous_pair(specs, 5)
        finally:
            campaign._restore_signal_handlers(previous)
        self.assertEqual(len(spawned), 1)
        self.assertFalse(campaign._group_exists(spawned[0]))

    def test_child_restores_original_mask_and_accepts_sigterm(self):
        pid_dir = Path(self.temporary.name) / "term-proof"
        pid_dir.mkdir()
        specs = self.fixture.specs(self.cpus)
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": "term-proof",
            "MOCK_PID_DIR": str(pid_dir),
        }
        result = []

        def invoke():
            try:
                result.append(campaign.run_simultaneous_pair(
                    specs, 5, environment=environment))
            except BaseException as error:
                result.append(error)

        worker = threading.Thread(target=invoke)
        worker.start()
        deadline = time.monotonic() + 2
        while len(list(pid_dir.iterdir())) != 2 and \
                time.monotonic() < deadline:
            time.sleep(0.01)
        pids = [int(path.read_text()) for path in pid_dir.iterdir()]
        self.assertEqual(len(pids), 2)
        for pid in pids:
            os.kill(pid, signal.SIGTERM)
        worker.join(2)
        if worker.is_alive():
            for pid in pids:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            worker.join(2)
        self.assertFalse(worker.is_alive())
        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[0], tuple)
        self.assertEqual(
            [capture.returncode for capture in result[0]], [19, 19])

    def test_double_sigterm_during_cleanup_cannot_escape_reaping(self):
        pid_dir = Path(self.temporary.name) / "double-term"
        pid_dir.mkdir()
        specs = self.fixture.specs(self.cpus)
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": "ignore-term-sleep",
            "MOCK_PID_DIR": str(pid_dir),
        }
        real_terminate = campaign._terminate_children

        def signal_twice_then_reap(children, pgids):
            signal.raise_signal(signal.SIGTERM)
            signal.raise_signal(signal.SIGTERM)
            return real_terminate(children, pgids)

        previous = campaign._install_signal_handlers()
        try:
            with mock.patch.object(
                    campaign, "_terminate_children",
                    side_effect=signal_twice_then_reap):
                with self.assertRaises(campaign.CampaignInterrupted):
                    campaign.run_simultaneous_pair(
                        specs, 0.2, environment=environment)
        finally:
            campaign._restore_signal_handlers(previous)
        child_pids = [
            int(path.read_text()) for path in pid_dir.iterdir()]
        self.assertEqual(len(child_pids), 2)
        for pid in child_pids:
            deadline = time.monotonic() + 2
            while (Path("/proc") / str(pid)).exists() and \
                    time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertFalse((Path("/proc") / str(pid)).exists())

    def test_sealed_execution_is_pinned_when_original_changes(self):
        binding = self.fixture.bindings["hook-on-a"]
        self.assertEqual(
            fcntl.fcntl(binding.execution_fd, campaign.F_GET_SEALS) &
            campaign.EXECUTION_SEALS,
            campaign.EXECUTION_SEALS)
        self.assertEqual(
            stat.S_IMODE(os.fstat(binding.execution_fd).st_mode),
            campaign.EXECUTION_MODE)
        with self.assertRaises(OSError):
            os.write(binding.execution_fd, b"x")
        with self.assertRaises(OSError):
            os.fchmod(binding.execution_fd, 0)
        with self.assertRaises(OSError):
            os.fchmod(binding.execution_fd, campaign.EXECUTION_MODE | 0o010)
        original = Path(binding.path).read_bytes()
        Path(binding.path).write_bytes(b"X" * len(original))
        specs, captures = self._captures()
        parsed = campaign.parse_timing_stream(
            captures[0].stdout, specs[0].expected)
        self.assertEqual(parsed.status, "success")
        self.assertEqual(
            captures[0].launched_command[0],
            "/proc/self/fd/{}".format(binding.execution_fd))
        with self.assertRaisesRegex(campaign.CampaignError, "changed"):
            campaign._reverify_binaries(self.fixture.bindings)

    def test_early_leader_exit_descendant_is_killed(self):
        pid_dir = Path(self.temporary.name) / "leader-pid"
        pid_dir.mkdir()
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": "leader-exit-descendant",
            "MOCK_PID_DIR": str(pid_dir),
        }
        specs = self.fixture.specs(self.cpus)
        with self.assertRaisesRegex(
                campaign.CampaignError, "before start barrier"):
            campaign.run_simultaneous_pair(
                specs, 3, environment=environment)
        descendant = int((pid_dir / "descendant").read_text())
        deadline = time.monotonic() + 2
        while (Path("/proc") / str(descendant)).exists() and \
                time.monotonic() < deadline:
            time.sleep(0.01)
        self.assertFalse((Path("/proc") / str(descendant)).exists())

    def test_weak_encoder_results_exact_set(self):
        for result in (1, 4, 8):
            with self.subTest(result=result):
                environment = {
                    "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
                    "MOCK_MODE": "weak",
                    "MOCK_WEAK_RESULT": str(result),
                }
                specs = self.fixture.specs(self.cpus, scope="row")
                captures = campaign.run_simultaneous_pair(
                    specs, 5, environment=environment)
                parsed = campaign.parse_timing_stream(
                    captures[0].stdout, specs[0].expected)
                self.assertEqual(parsed.status, "weak-root")
                self.assertIsNone(parsed.timing)
        for result in (3, 9):
            with self.subTest(invalid=result):
                environment = {
                    "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
                    "MOCK_MODE": "weak",
                    "MOCK_WEAK_RESULT": str(result),
                }
                specs = self.fixture.specs(self.cpus, scope="row")
                captures = campaign.run_simultaneous_pair(
                    specs, 5, environment=environment)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "weak-root"):
                    campaign.parse_timing_stream(
                        captures[0].stdout, specs[0].expected)

    def test_timeout_reaps_mock_process_groups(self):
        pid_dir = Path(self.temporary.name) / "pids"
        pid_dir.mkdir()
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": "sleep",
            "MOCK_PID_DIR": str(pid_dir),
        }
        specs = self.fixture.specs(self.cpus, scope="row")
        started = time.monotonic()
        with self.assertRaisesRegex(campaign.CampaignError, "timed out"):
            campaign.run_simultaneous_pair(
                specs, 0.25, environment=environment)
        self.assertLess(time.monotonic() - started, 1.0)
        deadline = time.monotonic() + 2
        child_pids = []
        while time.monotonic() < deadline:
            child_pids = [
                int(path.read_text()) for path in pid_dir.iterdir()]
            if len(child_pids) == 2:
                break
            time.sleep(0.01)
        self.assertEqual(len(child_pids), 2)
        for pid in child_pids:
            deadline = time.monotonic() + 2
            while Path("/proc/{}".format(pid)).exists() and \
                    time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertFalse(Path("/proc/{}".format(pid)).exists())


class DesignAndMathTests(unittest.TestCase):
    def test_frozen_coverage_and_period_cpu_balance(self):
        campaign.validate_frozen_design()
        self.assertEqual(len(campaign.frozen_preflight_cases()), 14 * 6)
        timed = campaign.frozen_timed_cases()
        self.assertEqual(
            sum(row["kind"] == "primary" for row in timed), 14 * 5)
        self.assertEqual(
            sum(row["kind"] == "anchor" for row in timed), 3 * 5)
        by_period = [
            {"on-aa": 0, "off-aa": 0, "cross": 0}
            for unused in range(4)]
        cpu_counts = {
            treatment: [0, 0] for treatment in campaign.TREATMENTS}
        exact_orientations = {}
        cpu_rows = {
            cpu: {row: 0 for row in range(4)}
            for cpu in ("cpu_a", "cpu_b")}
        cpu_period = {
            cpu: [
                {treatment: 0 for treatment in campaign.TREATMENTS}
                for unused in range(4)]
            for cpu in ("cpu_a", "cpu_b")}
        carryovers = {cpu: {} for cpu in ("cpu_a", "cpu_b")}
        self.assertEqual(campaign.WILLIAMS_ROWS, (
            ("hook-on-a", "hook-on-b", "hook-off-b", "hook-off-a"),
            ("hook-on-b", "hook-off-a", "hook-on-a", "hook-off-b"),
            ("hook-off-a", "hook-off-b", "hook-on-b", "hook-on-a"),
            ("hook-off-b", "hook-on-a", "hook-off-a", "hook-on-b"),
        ))
        for block in range(campaign.MEASURED_BLOCKS):
            rows = campaign.williams_block(block)
            sequences = {"cpu_a": [], "cpu_b": []}
            for cpu in ("cpu_a", "cpu_b"):
                row_ids = {
                    row["{}_row".format(cpu)] for row in rows}
                self.assertEqual(len(row_ids), 1)
                cpu_rows[cpu][row_ids.pop()] += 1
            for row in rows:
                pair = (
                    row["cpu_a_treatment"], row["cpu_b_treatment"])
                exact = tuple(sorted(pair))
                directions = exact_orientations.setdefault(
                    exact, {exact: 0, exact[::-1]: 0})
                directions[pair] += 1
                states = {campaign.TREATMENT_STATE[value] for value in pair}
                kind = "cross"
                if states == {1}:
                    kind = "on-aa"
                elif states == {0}:
                    kind = "off-aa"
                by_period[row["round"]][kind] += 1
                cpu_counts[pair[0]][0] += 1
                cpu_counts[pair[1]][1] += 1
                for cpu, treatment in (
                        ("cpu_a", pair[0]), ("cpu_b", pair[1])):
                    cpu_period[cpu][row["round"]][treatment] += 1
                    sequences[cpu].append(treatment)
            for cpu, sequence in sequences.items():
                for transition in zip(sequence, sequence[1:]):
                    carryovers[cpu][transition] = (
                        carryovers[cpu].get(transition, 0) + 1)
        self.assertEqual(
            by_period,
            [{"on-aa": 4, "off-aa": 4, "cross": 8}] * 4)
        self.assertTrue(all(value == [16, 16]
                            for value in cpu_counts.values()))
        self.assertEqual(len(exact_orientations), 4)
        self.assertEqual(set(exact_orientations), {
            tuple(sorted(pair))
            for pair in (
                ("hook-on-a", "hook-on-b"),
                ("hook-off-a", "hook-off-b"),
                ("hook-on-b", "hook-off-a"),
                ("hook-on-a", "hook-off-b"),
            )
        })
        self.assertTrue(all(
            sorted(directions.values()) == [8, 8]
            for directions in exact_orientations.values()))
        self.assertTrue(all(
            set(rows.values()) == {4} for rows in cpu_rows.values()))
        self.assertTrue(all(
            set(period.values()) == {4}
            for periods in cpu_period.values() for period in periods))
        directed_nonself = {
            (left, right)
            for left in campaign.TREATMENTS
            for right in campaign.TREATMENTS
            if left != right
        }
        self.assertTrue(all(
            set(rows) == directed_nonself and set(rows.values()) == {4}
            for rows in carryovers.values()))
        contract = campaign.pre_results_contract(
            "0" * 64, 1, 2, 1024, 0, 1, "1" * 64,
            runner_python_sha256="2" * 64)
        self.assertEqual(contract["measurement"]["blocks"], 16)
        self.assertEqual(
            contract["measurement"]
            ["scheduled_invocations_per_treatment"], 32)
        self.assertEqual(
            contract["measurement"]
            ["successful_timing_observations_per_treatment"], 32)
        self.assertEqual(
            contract["measurement"]
            ["weak_root_timing_observations_per_treatment"], 0)
        self.assertIsNone(
            contract["measurement"]["weak_root_summary"])
        self.assertEqual(contract["timed_anchors"], [
            {
                **anchor,
                "scopes": list(campaign.SCOPES),
                "loss_ppm": campaign.PRIMARY_LOSS_PPM,
            }
            for anchor in campaign.TIMED_ANCHORS
        ])
        self.assertIn(
            "16 Williams-block inferential sampling units",
            contract["statistics"]["effect_sampling_unit"])
        self.assertIn(
            "assumes independence",
            contract["statistics"]["ci_assumption"])

    def test_cpu_list_endpoints_are_canonical_and_bounded(self):
        self.assertEqual(campaign._parse_cpu_list("0-2,4"), {0, 1, 2, 4})
        for text in ("00", "01-2", "0-01", "9" * 4096):
            with self.subTest(text=text[:20]):
                with self.assertRaises(campaign.CampaignError):
                    campaign._parse_cpu_list(text)

    def test_ci_formula_constant_and_known_variance(self):
        same = campaign._mean_ci([0.25] * 24)
        self.assertEqual(same["mean_log_ratio"], 0.25)
        self.assertEqual(same["sample_sd_log_ratio"], 0.0)
        self.assertEqual(same["ci95_low_log_ratio"], 0.25)
        self.assertEqual(same["ci95_high_log_ratio"], 0.25)
        values = [float(value) for value in range(1, 25)]
        result = campaign._mean_ci(values)
        mean = 12.5
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / 23)
        half = campaign.T_CRITICAL_975[23] * sd / math.sqrt(24)
        self.assertAlmostEqual(result["mean_log_ratio"], mean, places=15)
        self.assertAlmostEqual(
            result["ci95_low_log_ratio"], mean - half, places=14)
        self.assertAlmostEqual(
            result["ci95_high_log_ratio"], mean + half, places=14)

    def test_cross_effect_ci_clusters_two_pairs_within_each_block(self):
        pairs = []
        ordinal = 0
        for block in range(campaign.MEASURED_BLOCKS):
            cross_index = 0
            for row in campaign.williams_block(block):
                treatments = (
                    row["cpu_a_treatment"], row["cpu_b_treatment"])
                states = {
                    campaign.TREATMENT_STATE[value]
                    for value in treatments}
                elapsed = {treatment: 1 for treatment in treatments}
                if states == {0, 1}:
                    exponent = block if cross_index == 0 else -block
                    cross_index += 1
                    on = next(
                        value for value in treatments
                        if campaign.TREATMENT_STATE[value])
                    off = next(
                        value for value in treatments
                        if not campaign.TREATMENT_STATE[value])
                    if exponent >= 0:
                        elapsed[on] = 1 << exponent
                    else:
                        elapsed[off] = 1 << (-exponent)
                outcomes = []
                for cpu, treatment in enumerate(treatments):
                    expected = campaign.InvocationExpected(
                        treatment=treatment,
                        K=64,
                        block_bytes=2,
                        construction_seed=1,
                        loss_seed=2,
                        loss_ppm=350000,
                        schedule="repair-only",
                        scope="row",
                        inner_reps=1,
                        max_working_mib=1,
                        context_sha256="0" * 64,
                        cpu=cpu,
                    )
                    spec = campaign.InvocationSpec(
                        ordinal=ordinal,
                        pair_ordinal=ordinal // 2,
                        phase="measured",
                        case_id="cluster",
                        block=block,
                        round_index=row["round"],
                        treatment=treatment,
                        cpu=cpu,
                        executable_fd=0,
                        expected=expected,
                        command=("mock",),
                    )
                    ordinal += 1
                    timing = {
                        "elapsed_ns": elapsed[treatment],
                        "minor_faults": 0,
                        "major_faults": 0,
                        "voluntary_context_switches": 0,
                        "involuntary_context_switches": 0,
                    }
                    parsed = campaign.ParsedTimingStream(
                        header={},
                        semantic={},
                        timing=timing,
                        done={},
                        status="success",
                    )
                    capture = campaign.ProcessCapture(
                        returncode=0,
                        stdout=b"",
                        stderr=b"",
                        launched_command=("mock",),
                    )
                    outcomes.append(campaign.InvocationOutcome(
                        spec=spec, capture=capture, parsed=parsed))
                pairs.append(tuple(outcomes))
            self.assertEqual(cross_index, 2)
        summary = campaign.summarize_measurements(pairs)
        effect = summary["hook_on_over_hook_off"]
        self.assertEqual(effect["raw_cross_pairs"], 32)
        self.assertEqual(effect["effect_blocks"], 16)
        self.assertEqual(effect["n"], 16)
        self.assertEqual(effect["degrees_of_freedom"], 15)
        self.assertEqual(effect["contrast_components"], [
            {
                "numerator": "hook-on-a",
                "denominator": "hook-off-b",
                "pairs": 16,
                "numerator_cpu_a": 8,
                "numerator_cpu_b": 8,
            },
            {
                "numerator": "hook-on-b",
                "denominator": "hook-off-a",
                "pairs": 16,
                "numerator_cpu_a": 8,
                "numerator_cpu_b": 8,
            },
        ])
        self.assertEqual(
            summary["hook_on_a_over_b"]["contrast"], {
                "numerator": "hook-on-a",
                "denominator": "hook-on-b",
                "pairs": 16,
                "numerator_cpu_a": 8,
                "numerator_cpu_b": 8,
            })
        self.assertEqual(
            summary["hook_off_a_over_b"]["contrast"], {
                "numerator": "hook-off-a",
                "denominator": "hook-off-b",
                "pairs": 16,
                "numerator_cpu_a": 8,
                "numerator_cpu_b": 8,
            })
        self.assertAlmostEqual(effect["mean_log_ratio"], 0.0, places=15)
        self.assertAlmostEqual(
            effect["sample_sd_log_ratio"], 0.0, places=15)
        corrupt = list(pairs)
        block_zero_cross = [
            index for index, pair in enumerate(corrupt[:4])
            if {
                campaign.TREATMENT_STATE[pair[0].spec.treatment],
                campaign.TREATMENT_STATE[pair[1].spec.treatment],
            } == {0, 1}
        ]
        self.assertEqual(len(block_zero_cross), 2)
        corrupt[block_zero_cross[1]] = corrupt[block_zero_cross[0]]
        with self.assertRaisesRegex(campaign.CampaignError, "duplicate"):
            campaign.summarize_measurements(corrupt)

    def test_reserved_integrations_fail_closed(self):
        parser = campaign._parser()
        common = [
            "--build-receipt", "/tmp/a",
            "--build-receipt-sha256", "0" * 64,
            "--raw-log", "/tmp/b", "--output", "/tmp/c",
            "--cpu-a", "0", "--cpu-b", "64",
            "--raw-root", "1", "--loss-root", "2",
            "--max-working-mib", "1",
        ]
        for option in ("--filler-controller", "--thermal-window-controller"):
            args = parser.parse_args(common + [option, "/tmp/controller"])
            with self.assertRaisesRegex(
                    campaign.CampaignError, "not implemented"):
                campaign.run_campaign(args)

    def test_cli_uint_conversion_is_length_bounded(self):
        for text in ("1" * 21, "9" * 100000):
            with self.subTest(digits=len(text)):
                with self.assertRaises(argparse.ArgumentTypeError):
                    campaign._canonical_uint(text)

    def test_signal_handler_raises_and_restores_without_external_signal(self):
        previous = campaign._install_signal_handlers()
        try:
            with self.assertRaises(campaign.CampaignInterrupted):
                signal.raise_signal(signal.SIGTERM)
        finally:
            campaign._restore_signal_handlers(previous)


class BuildReceiptTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.readelf = self.root / "readelf"
        _write_executable(
            self.readelf,
            "#!/usr/bin/python3\n"
            "import sys\n"
            "data=open(sys.argv[-1],encoding='utf-8').read()\n"
            "v='11112222' if '# on' in data else '33334444'\n"
            "print('Build ID: '+v)\n")
        self.binary_rows = {}
        for treatment in campaign.TREATMENTS:
            path = self.root / treatment
            state = "on" if "hook-on" in treatment else "off"
            _write_executable(
                path, "#!/bin/sh\nexit 0\n# {}\n".format(state))
            payload = path.read_bytes()
            current = path.stat()
            self.binary_rows[treatment] = {
                "path": str(path),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "mode": stat.S_IMODE(current.st_mode),
                "build_id": {
                    "value":
                        "11112222" if state == "on" else "33334444"},
            }
        self.receipt = self._receipt()

    def tearDown(self):
        self.temporary.cleanup()

    def _tool(self, path):
        payload = Path(path).read_bytes()
        return {
            "configured_path": str(path),
            "resolved_path": str(Path(path).resolve()),
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "version": {"stdout": {}, "stderr": {}},
        }

    def _receipt(self):
        python = Path(sys.executable).resolve()
        return {
            "schema": campaign.BUILD_SCHEMA,
            "source": {},
            "provenance_helper": {},
            "hook_macro": {
                "name": campaign.HOOK_MACRO,
                "off_contract": "macro-absent",
                "on_token": "-D{}=1".format(campaign.HOOK_MACRO),
            },
            "configuration": {
                "target": campaign.BUILD_TARGET,
                "binary_relative": "codec/wirehair_v2_hook_timing",
                "build_type": "Release",
                "build_order": [
                    "hook-on-a", "hook-off-a",
                    "hook-on-b", "hook-off-b"],
                "generator": "cmake-default",
                "definitions": [],
                "flags": {
                    "c_flags": "", "cxx_flags": "",
                    "exe_linker_flags": "",
                    "module_linker_flags": "",
                    "shared_linker_flags": "",
                },
                "build_shared": "OFF",
                "static_pic": "ON",
                "march_native": "OFF",
                "lto": "OFF",
                "pgo": "OFF",
                "jobs": 1,
                "build_environment": {},
            },
            "tools": {
                "python": self._tool(python),
                "git": {},
                "cmake": {},
                "readelf": self._tool(self.readelf),
            },
            "builds": {
                treatment: {
                    "name": treatment,
                    "hook_state":
                        "on" if "hook-on" in treatment else "off",
                    "replicate": treatment[-1],
                    "binary": self.binary_rows[treatment],
                }
                for treatment in campaign.TREATMENTS
            },
            "comparability": {
                "cache_normalized_sha256_identical": True,
                "compile_commands_normalized_sha256_identical": True,
                "link_normalized_sha256_identical": True,
                "toolchains_identical": True,
                "hook_on_ab_binary_identical": True,
                "hook_off_ab_binary_identical": True,
                "hook_on_off_binary_identical": False,
                "only_allowed_compile_difference":
                    "-D{}=1".format(campaign.HOOK_MACRO),
            },
            "timing_executed": False,
        }

    def _publish(self, receipt=None):
        receipt = self.receipt if receipt is None else receipt
        payload = campaign.canonical_json_bytes(receipt)
        path = self.root / "receipt.json"
        path.write_bytes(payload)
        return path, hashlib.sha256(payload).hexdigest()

    def test_valid_receipt_authenticates_binary_hashes_and_build_ids(self):
        path, digest = self._publish()
        unused, bindings, binding = campaign.authenticate_build_receipt(
            path, digest)
        try:
            self.assertEqual(set(bindings), set(campaign.TREATMENTS))
            self.assertEqual(binding["sha256"], digest)
        finally:
            campaign.close_execution_bindings(bindings)

    def test_authentication_rejects_final_component_symlinks(self):
        receipt_path, digest = self._publish()
        receipt_link = self.root / "receipt-link.json"
        receipt_link.symlink_to(receipt_path)
        with self.assertRaisesRegex(campaign.CampaignError, "cannot open"):
            campaign.authenticate_build_receipt(receipt_link, digest)

        for kind in ("binary", "readelf"):
            with self.subTest(kind=kind):
                mutated = json.loads(json.dumps(self.receipt))
                if kind == "binary":
                    target = Path(
                        mutated["builds"]["hook-on-a"]["binary"]["path"])
                    link = self.root / "binary-link"
                    link.symlink_to(target)
                    mutated["builds"]["hook-on-a"]["binary"]["path"] = \
                        str(link)
                else:
                    target = Path(
                        mutated["tools"]["readelf"]["resolved_path"])
                    link = self.root / "readelf-link"
                    link.symlink_to(target)
                    mutated["tools"]["readelf"]["resolved_path"] = str(link)
                path, mutated_digest = self._publish(mutated)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "cannot open"):
                    campaign.authenticate_build_receipt(
                        path, mutated_digest)
                link.unlink()

    def test_receipt_hash_and_configuration_mutations_fail(self):
        path, digest = self._publish()
        with self.assertRaisesRegex(campaign.CampaignError, "mismatch"):
            campaign.authenticate_build_receipt(path, "f" * 64)
        for key, value in (
                ("build_type", "Debug"),
                ("binary_relative", "wrong"),
                ("build_order", list(campaign.TREATMENTS)),
                ("definitions", [{"name": "X", "value": "1"}]),
                ("build_shared", "ON"),
                ("static_pic", "OFF"),
                ("march_native", "ON"),
                ("lto", "ON"),
                ("pgo", "GENERATE")):
            with self.subTest(key=key):
                mutated = json.loads(json.dumps(self.receipt))
                mutated["configuration"][key] = value
                payload = campaign.canonical_json_bytes(mutated)
                mutated_path = self.root / "mutated-{}.json".format(key)
                mutated_path.write_bytes(payload)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "configuration"):
                    campaign.authenticate_build_receipt(
                        mutated_path, hashlib.sha256(payload).hexdigest())

    def test_self_inconsistent_identical_on_off_binding_fails(self):
        # Make the off files byte-identical to on, update every authenticated
        # field except retain the distinct claimed off build ID, and lie
        # consistently in comparability.  Independent byte authentication must
        # reject the experiment without trusting the claimed build IDs.
        on_payload = Path(
            self.binary_rows["hook-on-a"]["path"]).read_bytes()
        mutated = json.loads(json.dumps(self.receipt))
        for treatment in ("hook-off-a", "hook-off-b"):
            path = Path(mutated["builds"][treatment]["binary"]["path"])
            path.write_bytes(on_payload)
            path.chmod(0o755)
            mutated["builds"][treatment]["binary"]["bytes"] = len(on_payload)
            mutated["builds"][treatment]["binary"]["sha256"] = \
                hashlib.sha256(on_payload).hexdigest()
            mutated["builds"][treatment]["binary"]["build_id"]["value"] = \
                "33334444"
        mutated["comparability"]["hook_on_off_binary_identical"] = False
        payload = campaign.canonical_json_bytes(mutated)
        path = self.root / "identical.json"
        path.write_bytes(payload)
        with self.assertRaisesRegex(campaign.CampaignError, "identical"):
            campaign.authenticate_build_receipt(
                path, hashlib.sha256(payload).hexdigest())

    def test_same_state_mode_mismatch_is_rejected(self):
        treatment = "hook-on-b"
        binary = self.receipt["builds"][treatment]["binary"]
        path = Path(binary["path"])
        path.chmod(0o700)
        binary["mode"] = 0o700
        receipt_path, digest = self._publish()
        with self.assertRaisesRegex(campaign.CampaignError, "A/B"):
            campaign.authenticate_build_receipt(receipt_path, digest)

    def test_cross_state_mode_mismatch_is_rejected(self):
        for treatment in ("hook-on-a", "hook-on-b"):
            binary = self.receipt["builds"][treatment]["binary"]
            path = Path(binary["path"])
            path.chmod(0o700)
            binary["mode"] = 0o700
        receipt_path, digest = self._publish()
        with self.assertRaisesRegex(campaign.CampaignError, "modes"):
            campaign.authenticate_build_receipt(receipt_path, digest)

    def test_owner_execute_bit_is_required(self):
        for treatment in campaign.TREATMENTS:
            binary = self.receipt["builds"][treatment]["binary"]
            path = Path(binary["path"])
            path.chmod(0o610)
            binary["mode"] = 0o610
        receipt_path, digest = self._publish()
        with self.assertRaisesRegex(
                campaign.CampaignError, "owner-executable"):
            campaign.authenticate_build_receipt(receipt_path, digest)

    def test_partial_auth_failure_clears_sink_before_fd_reuse(self):
        mutated = json.loads(json.dumps(self.receipt))
        mutated["builds"]["hook-on-b"]["binary"]["sha256"] = "0" * 64
        receipt_path, digest = self._publish(mutated)
        sink = {}
        created = []
        real_create = campaign._sealed_executable_memfd

        def capture_fd(*args, **kwargs):
            descriptor = real_create(*args, **kwargs)
            created.append(descriptor)
            return descriptor

        with mock.patch.object(
                campaign, "_sealed_executable_memfd",
                side_effect=capture_fd):
            with self.assertRaisesRegex(campaign.CampaignError, "binding"):
                campaign.authenticate_build_receipt(
                    receipt_path, digest, sink)
        self.assertEqual(len(created), 1)
        self.assertEqual(sink, {})
        replacement = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        try:
            self.assertEqual(replacement, created[0])
            campaign.close_execution_bindings(sink)
            os.fstat(replacement)
        finally:
            os.close(replacement)


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.output = self.root / "receipt.json"
        self.payload = b'{"ok":true}\n'

    def tearDown(self):
        self.temporary.cleanup()

    def test_success_and_prelink_failure_boundaries(self):
        publication = campaign._publish_no_replace(
            self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(publication["post_commit_warnings"], [])
        self.assertEqual(publication["post_commit_interruptions"], [])
        self.assertEqual(self.output.read_bytes(), self.payload)

        second = self.root / "second.json"
        with mock.patch.object(
                campaign.os, "link",
                side_effect=OSError("prelink injection")):
            with self.assertRaisesRegex(
                    campaign.CampaignError, "prelink injection"):
                campaign._publish_no_replace(second, self.payload)
        self.assertFalse(second.exists())

    def test_precommit_interruption_leaves_no_output(self):
        with mock.patch.object(
                campaign.os, "fsync",
                side_effect=campaign.CampaignInterrupted(
                    "write interruption")):
            with self.assertRaisesRegex(
                    campaign.CampaignInterrupted, "write interruption"):
                campaign._publish_no_replace(self.output, self.payload)
        self.assertFalse(self.output.exists())

    def test_concurrent_output_creation_is_preserved(self):
        real_link = campaign.os.link
        real_open = campaign.os.open
        real_close = campaign.os.close
        competing = b"concurrent owner\n"

        def create_output_then_link(source, destination, **kwargs):
            descriptor = real_open(
                destination,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o644,
                dir_fd=kwargs["dst_dir_fd"],
            )
            try:
                os.write(descriptor, competing)
                os.fsync(descriptor)
            finally:
                real_close(descriptor)
            return real_link(source, destination, **kwargs)

        with mock.patch.object(
                campaign.os, "link", side_effect=create_output_then_link):
            with self.assertRaisesRegex(
                    campaign.CampaignError, "File exists"):
                campaign._publish_no_replace(self.output, self.payload)
        self.assertEqual(self.output.read_bytes(), competing)

    def test_publication_links_the_held_anonymous_inode(self):
        real_link = campaign.os.link
        linked_inode = []

        def record_source_inode(source, destination, **kwargs):
            source_fd = int(source.rsplit("/", 1)[1])
            linked_inode.append(os.fstat(source_fd).st_ino)
            return real_link(source, destination, **kwargs)

        stale_named_temporary = self.output.with_name(
            ".{}.{}.tmp".format(self.output.name, os.getpid()))
        stale_named_temporary.write_bytes(b"forged\n")
        with mock.patch.object(
                campaign.os, "link", side_effect=record_source_inode):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertEqual(self.output.stat().st_ino, linked_inode[0])
        self.assertEqual(stale_named_temporary.read_bytes(), b"forged\n")

    def test_exception_after_real_link_is_committed(self):
        real_link = campaign.os.link

        def link_then_raise(*args, **kwargs):
            real_link(*args, **kwargs)
            raise MemoryError("post-link injection")

        with mock.patch.object(
                campaign.os, "link", side_effect=link_then_raise):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertTrue(any(
            "MemoryError" in warning
            for warning in publication["post_commit_warnings"]))

    def test_stat_failure_after_real_link_is_committed(self):
        real_stat = campaign.os.stat
        output_stats = []
        before_fds = set(os.listdir("/proc/self/fd"))

        def fail_postlink_stat(path, *args, **kwargs):
            if (path == self.output.name and
                    kwargs.get("dir_fd") is not None):
                output_stats.append(path)
                if len(output_stats) == 2:
                    raise OSError("post-link stat injection")
            return real_stat(path, *args, **kwargs)

        with mock.patch.object(
                campaign.os, "stat", side_effect=fail_postlink_stat):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertEqual(before_fds, set(os.listdir("/proc/self/fd")))
        self.assertTrue(any(
            "verification failed after commit" in warning
            for warning in publication["post_commit_warnings"]))

    def test_link_and_verification_failure_is_explicitly_ambiguous(self):
        real_stat = campaign.os.stat
        output_stats = []

        def fail_verification(path, *args, **kwargs):
            if (path == self.output.name and
                    kwargs.get("dir_fd") is not None):
                output_stats.append(path)
                if len(output_stats) == 2:
                    raise OSError("verification injection")
            return real_stat(path, *args, **kwargs)

        with mock.patch.object(
                campaign.os, "link",
                side_effect=OSError("link outcome injection")), \
                mock.patch.object(
                    campaign.os, "stat", side_effect=fail_verification):
            with self.assertRaisesRegex(
                    campaign.CampaignError, "outcome is ambiguous"):
                campaign._publish_no_replace(
                    self.output, self.payload)
        self.assertFalse(self.output.exists())

    def test_postcommit_directory_sync_is_a_warning(self):
        fsync_calls = []
        real_fsync = campaign.os.fsync

        def fail_directory_sync(descriptor):
            fsync_calls.append(descriptor)
            if len(fsync_calls) == 2:
                raise OSError("directory fsync injection")
            return real_fsync(descriptor)

        with mock.patch.object(
                campaign.os, "fsync", side_effect=fail_directory_sync):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertTrue(any(
            "directory sync" in warning
            for warning in publication["post_commit_warnings"]))

    def test_parent_drift_does_not_skip_held_directory_fsync(self):
        parent = self.root / "publish-parent"
        renamed = self.root / "renamed-parent"
        parent.mkdir()
        output = parent / "receipt.json"
        real_stat = campaign.os.stat
        real_fsync = campaign.os.fsync
        fsync_calls = []
        moved = []

        def move_parent_before_identity_check(path, *args, **kwargs):
            if Path(path) == parent and not moved:
                try:
                    real_stat(output)
                except FileNotFoundError:
                    pass
                else:
                    parent.rename(renamed)
                    parent.mkdir()
                    moved.append(True)
            return real_stat(path, *args, **kwargs)

        def track_fsync(descriptor):
            fsync_calls.append(descriptor)
            return real_fsync(descriptor)

        with mock.patch.object(
                campaign.os, "stat",
                side_effect=move_parent_before_identity_check), \
                mock.patch.object(
                    campaign.os, "fsync", side_effect=track_fsync):
            publication = campaign._publish_no_replace(
                output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(len(fsync_calls), 2)
        self.assertFalse(output.exists())
        self.assertEqual(
            (renamed / output.name).read_bytes(), self.payload)
        self.assertTrue(any(
            "parent validation failed" in warning
            for warning in publication["post_commit_warnings"]))

    def test_postcommit_signal_mask_restore_failure_is_a_warning(self):
        real_pthread_sigmask = campaign.signal.pthread_sigmask

        def restore_then_raise(how, mask):
            result = real_pthread_sigmask(how, mask)
            if how == signal.SIG_SETMASK:
                raise OSError("mask restore injection")
            return result

        with mock.patch.object(
                campaign.signal, "pthread_sigmask",
                side_effect=restore_then_raise):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertTrue(any(
            "signal-mask restoration" in warning
            for warning in publication["post_commit_warnings"]))
        self.assertEqual(publication["post_commit_interruptions"], [])

    def test_real_pending_signal_is_reported_after_commit(self):
        original_mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
        if signal.SIGTERM in original_mask:
            self.skipTest("SIGTERM is blocked by the test runner")
        real_link = campaign.os.link

        def link_then_signal(*args, **kwargs):
            result = real_link(*args, **kwargs)
            os.kill(os.getpid(), signal.SIGTERM)
            self.assertIn(signal.SIGTERM, signal.sigpending())
            return result

        previous_handlers = campaign._install_signal_handlers()
        try:
            with mock.patch.object(
                    campaign.os, "link", side_effect=link_then_signal):
                publication = campaign._publish_no_replace(
                    self.output, self.payload)
        finally:
            campaign._restore_signal_handlers(previous_handlers)
            signal.pthread_sigmask(signal.SIG_SETMASK, original_mask)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertTrue(publication["post_commit_interruptions"])
        self.assertTrue(any(
            "interrupted after receipt link" in warning
            for warning in publication["post_commit_warnings"]))

    def test_postcommit_close_interruption_is_classified(self):
        real_close = campaign.os.close
        interrupted = []
        before_fds = set(os.listdir("/proc/self/fd"))

        def close_then_interrupt(descriptor):
            real_close(descriptor)
            if self.output.exists() and not interrupted:
                interrupted.append(descriptor)
                raise campaign.CampaignInterrupted(
                    "close interruption injection")

        with mock.patch.object(
                campaign.os, "close", side_effect=close_then_interrupt):
            publication = campaign._publish_no_replace(
                self.output, self.payload)
        self.assertTrue(publication["committed"])
        self.assertEqual(self.output.read_bytes(), self.payload)
        self.assertEqual(before_fds, set(os.listdir("/proc/self/fd")))
        self.assertTrue(any(
            "close interruption injection" in message
            for message in publication["post_commit_interruptions"]))

    def test_initial_signal_block_exception_restores_exact_mask(self):
        real_pthread_sigmask = campaign.signal.pthread_sigmask
        initial_mask = real_pthread_sigmask(signal.SIG_BLOCK, set())
        changed = []

        def block_then_raise(how, mask):
            result = real_pthread_sigmask(how, mask)
            if (how == signal.SIG_BLOCK and
                    set(mask) == set(campaign.MANAGED_SIGNALS) and
                    not changed):
                changed.append(True)
                raise OSError("signal block injection")
            return result

        try:
            with mock.patch.object(
                    campaign.signal, "pthread_sigmask",
                    side_effect=block_then_raise):
                with self.assertRaisesRegex(
                        campaign.CampaignError, "signal block injection"):
                    campaign._publish_no_replace(
                        self.output, self.payload)
            self.assertEqual(
                real_pthread_sigmask(signal.SIG_BLOCK, set()),
                initial_mask)
            self.assertFalse(self.output.exists())
        finally:
            real_pthread_sigmask(signal.SIG_SETMASK, initial_mask)

    def test_final_component_symlink_is_not_followed(self):
        target = self.root / "redirected.json"
        self.output.symlink_to(target)
        with self.assertRaisesRegex(campaign.CampaignError, "exists"):
            campaign._publish_no_replace(self.output, self.payload)
        self.assertTrue(self.output.is_symlink())
        self.assertFalse(target.exists())

    def test_cli_publication_line_surfaces_postcommit_warnings(self):
        publication = {
            "path": str(self.output),
            "sha256": "1" * 64,
            "committed": True,
            "post_commit_interruptions": [],
            "post_commit_warnings": ["directory sync warning"],
        }
        fake_parser = mock.Mock()
        fake_parser.parse_args.return_value = object()
        output = io.StringIO()
        with mock.patch.object(campaign, "_parser", return_value=fake_parser), \
                mock.patch.object(
                    campaign, "run_campaign",
                    return_value=({}, publication)), \
                contextlib.redirect_stdout(output):
            campaign.main([])
        row = json.loads(output.getvalue())
        self.assertEqual(row["post_commit_warnings"], [
            "directory sync warning"])
        self.assertTrue(row["committed"])

    def test_cli_surfaces_commit_before_late_interrupt(self):
        publication = {
            "path": str(self.output),
            "sha256": "1" * 64,
            "committed": True,
            "post_commit_interruptions": ["late SIGTERM"],
            "post_commit_warnings": ["late SIGTERM"],
        }
        fake_parser = mock.Mock()
        fake_parser.parse_args.return_value = object()
        output = io.StringIO()
        with mock.patch.object(campaign, "_parser", return_value=fake_parser), \
                mock.patch.object(
                    campaign, "run_campaign",
                    return_value=({}, publication)), \
                contextlib.redirect_stdout(output):
            with self.assertRaisesRegex(
                    campaign.CampaignInterrupted, "after the receipt commit"):
                campaign.main([])
        row = json.loads(output.getvalue())
        self.assertTrue(row["committed"])
        self.assertEqual(row["post_commit_interruptions"], ["late SIGTERM"])


class CampaignCleanupTests(unittest.TestCase):
    def _populate_sink(self, sink, descriptors):
        payload = b"#!/bin/sh\nexit 0\n"
        digest = hashlib.sha256(payload).hexdigest()
        for treatment in campaign.TREATMENTS:
            descriptor = campaign._sealed_executable_memfd(
                "cleanup-{}".format(treatment), payload)
            descriptors.append(descriptor)
            sink[treatment] = campaign.BinaryBinding(
                treatment=treatment,
                hook_state=campaign.TREATMENT_STATE[treatment],
                replicate=treatment[-1],
                path="/unused/{}".format(treatment),
                bytes=len(payload),
                sha256=digest,
                build_id="11112222",
                device=1,
                inode=len(descriptors),
                mode=0o755,
                execution_fd=descriptor,
            )

    def test_run_campaign_closes_execution_fds_once_on_success_and_error(self):
        for failure in (False, True):
            with self.subTest(failure=failure):
                descriptors = []
                close_counts = {}
                real_close = campaign.os.close

                def fake_unclosed(unused_args, sink):
                    self._populate_sink(sink, descriptors)
                    if failure:
                        raise campaign.CampaignError("injected campaign error")
                    return {"ok": True}, {"committed": True}

                def tracking_close(descriptor):
                    if descriptor in descriptors:
                        close_counts[descriptor] = (
                            close_counts.get(descriptor, 0) + 1)
                    return real_close(descriptor)

                with mock.patch.object(
                        campaign, "_run_campaign_unclosed",
                        side_effect=fake_unclosed), \
                        mock.patch.object(
                            campaign.os, "close",
                            side_effect=tracking_close):
                    if failure:
                        with self.assertRaisesRegex(
                                campaign.CampaignError, "injected"):
                            campaign.run_campaign(object())
                    else:
                        campaign.run_campaign(object())
                self.assertEqual(len(descriptors), 4)
                self.assertEqual(
                    close_counts,
                    {descriptor: 1 for descriptor in descriptors})
                for descriptor in descriptors:
                    with self.assertRaises(OSError):
                        os.fstat(descriptor)

    def test_raw_close_failure_is_surfaced_or_preserves_primary(self):
        class StubRaw:
            def __init__(self, error=None):
                self.closed = False
                self.error = error
                self.calls = 0

            def close(self):
                self.calls += 1
                if self.error is not None:
                    raise self.error
                self.closed = True

        success = StubRaw()
        campaign._close_raw_log(success, False)
        self.assertTrue(success.closed)
        self.assertEqual(success.calls, 1)

        only_failure = StubRaw(OSError("fsync injection"))
        with self.assertRaisesRegex(campaign.CampaignError, "fsync injection"):
            campaign._close_raw_log(only_failure, False)
        self.assertEqual(only_failure.calls, 1)

        primary = StubRaw(OSError("cleanup injection"))
        campaign._close_raw_log(primary, True)
        self.assertEqual(primary.calls, 1)

    def test_real_raw_close_error_marks_descriptor_closed(self):
        with tempfile.TemporaryDirectory() as root:
            raw = campaign.AppendOnlyRawLog(Path(root) / "raw.jsonl")
            descriptor = raw.descriptor
            with mock.patch.object(
                    campaign.os, "fsync",
                    side_effect=OSError("fsync injection")):
                with self.assertRaisesRegex(
                        campaign.CampaignError, "fsync injection"):
                    raw.close()
            self.assertTrue(raw.closed)
            self.assertEqual(raw.descriptor, -1)
            with self.assertRaises(OSError):
                os.fstat(descriptor)
            campaign._close_raw_log(raw, True)

    def test_raw_parent_fsync_failure_prevents_output_publication(self):
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            args = argparse.Namespace(
                build_receipt=root / "build.json",
                build_receipt_sha256="0" * 64,
                raw_log=root / "raw.jsonl",
                output=root / "output.json",
                cpu_a=0,
                cpu_b=1,
                raw_root=1,
                loss_root=2,
                max_working_mib=1024,
                timeout_seconds=5,
                calibration_max_inner_reps=10,
                calibration_max_steps=2,
                filler_controller=None,
                thermal_window_controller=None,
            )
            descriptors = []

            def fake_auth(unused_path, unused_sha, sink):
                self._populate_sink(sink, descriptors)
                return {}, sink, {}

            class FakeRunner:
                def __init__(self, **unused_kwargs):
                    self.next_invocation = 0
                    self.next_pair = 0
                    self.execution_order = []

                def run_semantic_preflight(self):
                    return [{} for unused in range(84)]

                def run_timed_case(self, case):
                    return dict(case)

            real_sync = campaign._fsync_directory
            raw_parent_syncs = []

            def fail_final_raw_parent_sync(path, description):
                if description == "raw-log parent":
                    raw_parent_syncs.append(str(path))
                    if len(raw_parent_syncs) == 2:
                        raise campaign.CampaignError(
                            "raw parent fsync injection")
                return real_sync(path, description)

            with mock.patch.object(
                    campaign, "authenticate_build_receipt",
                    side_effect=fake_auth), \
                    mock.patch.object(
                        campaign, "validate_smt_siblings",
                        return_value={}), \
                    mock.patch.object(campaign, "CampaignRunner", FakeRunner), \
                    mock.patch.object(
                        campaign, "_reverify_binaries",
                        return_value=None), \
                    mock.patch.object(
                        campaign, "_fsync_directory",
                        side_effect=fail_final_raw_parent_sync), \
                    mock.patch.object(
                        campaign, "_publish_no_replace") as publish:
                with self.assertRaisesRegex(
                        campaign.CampaignError,
                        "raw parent fsync injection"):
                    campaign.run_campaign(args)
            self.assertEqual(len(raw_parent_syncs), 2)
            publish.assert_not_called()
            self.assertFalse(args.output.exists())

    def test_raw_and_output_parents_must_exist(self):
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            missing = root / "missing"
            with self.assertRaisesRegex(
                    campaign.CampaignError, "parent"):
                campaign.AppendOnlyRawLog(missing / "raw.jsonl")
            with self.assertRaisesRegex(
                    campaign.CampaignError, "parent"):
                campaign._publish_no_replace(
                    missing / "output.json", b"{}\n")

            redirected = root / "redirected-raw.jsonl"
            raw_link = root / "raw-link.jsonl"
            raw_link.symlink_to(redirected)
            with self.assertRaisesRegex(
                    campaign.CampaignError, "cannot create raw log"):
                campaign.AppendOnlyRawLog(raw_link)
            self.assertTrue(raw_link.is_symlink())
            self.assertFalse(redirected.exists())

            raw_path = root / "sync-failure.jsonl"
            real_close = campaign.os.close

            def close_then_raise(descriptor):
                real_close(descriptor)
                raise OSError("cleanup close injection")

            with mock.patch.object(
                    campaign, "_fsync_directory",
                    side_effect=campaign.CampaignError(
                        "primary directory sync failure")), \
                    mock.patch.object(
                        campaign.os, "close",
                        side_effect=close_then_raise):
                with self.assertRaisesRegex(
                        campaign.CampaignError,
                        "primary directory sync failure"):
                    campaign.AppendOnlyRawLog(raw_path)
            self.assertFalse(raw_path.exists())


@unittest.skipUnless(_sibling_pair(), "requires allowed SMT siblings")
class RunnerFailureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = MockFixture(self.temporary.name)
        self.cpus = _sibling_pair()

    def tearDown(self):
        self.fixture.close()
        self.temporary.cleanup()

    def _runner(self, mode, inner_cap=1000000):
        raw = campaign.AppendOnlyRawLog(
            Path(self.temporary.name) / "raw-{}.jsonl".format(mode))
        environment = {
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
            "MOCK_MODE": mode,
        }
        runner = campaign.CampaignRunner(
            bindings=self.fixture.bindings,
            raw_log=raw,
            context_sha256="0" * 64,
            raw_root=1,
            loss_root=2,
            max_working_mib=1024,
            cpu_a=self.cpus[0],
            cpu_b=self.cpus[1],
            timeout_seconds=5,
            calibration_max_inner_reps=inner_cap,
            calibration_max_steps=4,
            process_environment=environment,
        )
        return runner, raw

    def _case(self, scope="row"):
        return {
            "case_id": "test",
            "kind": "primary",
            "K": 64,
            "block_bytes": 2,
            "schedule": "repair-only",
            "loss_ppm": 350000,
            "scope": scope,
        }

    def test_semantic_drift_and_weak_disagreement_fail(self):
        for mode, error in (
                ("semantic-drift", "semantic drift"),
                ("weak-disagreement", "semantic drift")):
            with self.subTest(mode=mode):
                runner, raw = self._runner(mode)
                try:
                    with self.assertRaisesRegex(
                            campaign.CampaignError, error):
                        runner.execute_pair(
                            "hook-on-a", "hook-off-a",
                            "test", self._case(scope="semantic"), 1)
                finally:
                    binding = raw.close()
                records = [
                    json.loads(line)
                    for line in Path(binding["path"]).read_text().splitlines()]
                self.assertEqual(len(records), 2)
                self.assertTrue(all(
                    row["validation"]["ok"] is False for row in records))
                self.assertTrue(all(
                    "semantic drift" in row["validation"]["error"]
                    for row in records))

    def test_one_bad_stream_marks_valid_peer_as_pair_invalid(self):
        runner, raw = self._runner("wrong-hooks-on")
        try:
            with self.assertRaises(campaign.CampaignError):
                runner.execute_pair(
                    "hook-on-a", "hook-off-a",
                    "test", self._case(scope="semantic"), 1)
        finally:
            binding = raw.close()
        records = [
            json.loads(line)
            for line in Path(binding["path"]).read_text().splitlines()]
        self.assertEqual(len(records), 2)
        self.assertFalse(records[0]["validation"]["ok"])
        self.assertIn("header", records[0]["validation"]["error"])
        self.assertFalse(records[1]["validation"]["ok"])
        self.assertIn("peer", records[1]["validation"]["error"])

    def test_sub_100ms_calibration_hits_cap_fail_closed(self):
        runner, raw = self._runner("slow", inner_cap=10)
        try:
            with self.assertRaisesRegex(
                    campaign.CampaignError, "above cap"):
                runner.calibrate(self._case())
        finally:
            binding = raw.close()
        self.assertEqual(binding["records"], 4)

    def test_complete_timed_case_has_frozen_invocation_counts(self):
        runner, raw = self._runner("normal")
        try:
            result = runner.run_timed_case(self._case())
        finally:
            binding = raw.close()
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["calibration"]["inner_reps"], 1)
        self.assertEqual(runner.next_invocation, 148)
        self.assertEqual(runner.next_pair, 74)
        self.assertEqual(binding["records"], 148)
        self.assertEqual(len(result["warmup_pairs"]), 8)
        self.assertEqual(len(result["measured_pairs"]), 64)
        self.assertEqual(
            result["summary"]["hook_on_over_hook_off"]["n"], 16)
        self.assertEqual(
            result["summary"]["hook_on_over_hook_off"]
            ["degrees_of_freedom"], 15)
        self.assertEqual(
            result["summary"]["hook_on_over_hook_off"]
            ["raw_cross_pairs"], 32)
        self.assertEqual(
            result["summary"]["hook_on_over_hook_off"]
            ["effect_blocks"], 16)
        counts = {
            treatment: sum(
                row["treatment"] == treatment
                for row in runner.execution_order)
            for treatment in campaign.TREATMENTS
        }
        self.assertEqual(set(counts.values()), {37})

    def test_append_many_preparation_failure_does_not_advance_state(self):
        path = Path(self.temporary.name) / "raw-prepare.jsonl"
        raw = campaign.AppendOnlyRawLog(path)
        try:
            with self.assertRaises(campaign.CampaignError):
                raw.append_many([{
                    "raw_ordinal": 0,
                    "bad": object(),
                }])
            self.assertEqual(raw.records, 0)
            self.assertEqual(raw.bytes, 0)
        finally:
            binding = raw.close()
        self.assertEqual(binding["records"], 0)
        self.assertEqual(path.read_bytes(), b"")


if __name__ == "__main__":
    unittest.main()
