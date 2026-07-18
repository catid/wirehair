#!/usr/bin/env python3
"""Unit tests for the isolated preferred-attempt timing layer.

All process execution is mocked.  These tests exercise the exact CLI/output
contract, retry/evidence behavior, and rational gates without running a timing
campaign or allocating the LLC-eviction buffer.
"""

from dataclasses import replace
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import wh2_preferred_attempt_timing as timing


ROUTE_CONTEXT_SHA256 = "b" * 64
TRACE_SHA256 = "c" * 64
THERMAL_SHA256 = "d" * 64
PERFORMANCE_SHA256 = "e" * 64


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def canonical_json(value):
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ensure_ascii=True,
    ) + "\n").encode("ascii")


def panel_spec(metric="solve", K=4096, width=64, schedule="burst"):
    return timing.TimingPanelSpec(
        K=K, block_bytes=width, preferred_attempt=2, metric=metric,
        seed=timing.timing_panel_seed(metric, K, width, schedule),
        schedule=schedule,
    )


def timing_stdout(
    spec,
    route_cache_sha256="a" * 64,
    trace_sha256=TRACE_SHA256,
    candidate_ns=90,
    control_ns=100,
    cpu_by_cycle=None,
    major_fault_by_cycle=None,
    cycle_index=None,
    intermediate_blocks=3,
    canonical_attempt=0,
):
    cpu_by_cycle = cpu_by_cycle or {}
    major_fault_by_cycle = major_fault_by_cycle or {}
    preamble = (
        "# preferredtiming: schema=v1 policy=h12-q0-adaptive "
        "metric=%s cycles=%s order=CTTCTCCT discard_cycle=0 "
        "cycle_mode=%s cycle_index=%s "
        "overhead=%s payload=%s payload_alignment=%s "
        "payload_prefaulted=%s route_cache_sha256=%s "
        "route_context_sha256=%s trace_sha256=%s" % (
            spec.metric,
            "4" if cycle_index is None else "1",
            "full" if cycle_index is None else "replacement",
            "all" if cycle_index is None else str(cycle_index),
            "4" if spec.metric == "solve" else "none",
            "distinct-zero-v1" if spec.metric == "solve" else "none",
            "64" if spec.metric == "solve" else "none",
            "1" if spec.metric == "solve" else "none",
            route_cache_sha256, ROUTE_CONTEXT_SHA256, trace_sha256,
        )
    )
    lines = [preamble, ",".join(timing.CSV_FIELDS)]
    cycles = range(4) if cycle_index is None else (cycle_index,)
    for cycle in cycles:
        cpu_before, cpu_after = cpu_by_cycle.get(cycle, (4, 4))
        for slot, marker in enumerate(timing.ORDER):
            candidate = marker == "T"
            distinct_candidate = (
                candidate and canonical_attempt != spec.preferred_attempt)
            if spec.metric == "solve":
                inactivated = 6 if distinct_candidate else 7
                binary_def = 2 if distinct_candidate else 3
                heavy_gain = 1
                block_xors = 1001 if distinct_candidate else 1003
                block_muladds = 71 if distinct_candidate else 73
                source_bytes = spec.K * spec.block_bytes
                intermediate_bytes = (
                    spec.K + intermediate_blocks) * spec.block_bytes
            else:
                inactivated = binary_def = heavy_gain = 0
                block_xors = block_muladds = 0
                source_bytes = intermediate_bytes = 0
            elapsed = candidate_ns if candidate else control_ns
            if isinstance(elapsed, dict):
                elapsed = elapsed[cycle]
            attempt = spec.preferred_attempt if candidate else canonical_attempt
            values = (
                spec.K, spec.block_bytes, spec.metric, cycle, slot,
                "candidate" if candidate else "control", attempt, 0,
                elapsed, 0, cpu_before, cpu_after, 0,
                major_fault_by_cycle.get(cycle, 0), inactivated, binary_def,
                heavy_gain, block_xors, block_muladds, source_bytes,
                intermediate_bytes,
            )
            lines.append(",".join(str(value) for value in values))
    return ("\n".join(lines) + "\n").encode("ascii")


def mutate_csv_field(raw, data_row, field, value):
    lines = raw.decode("ascii").splitlines()
    fields = lines[data_row + 2].split(",")
    fields[timing.CSV_FIELDS.index(field)] = str(value)
    lines[data_row + 2] = ",".join(fields)
    return ("\n".join(lines) + "\n").encode("ascii")


def mutate_arm_field(raw, arm, field, value):
    lines = raw.decode("ascii").splitlines()
    arm_index = timing.CSV_FIELDS.index("arm")
    field_index = timing.CSV_FIELDS.index(field)
    for line_index in range(2, len(lines)):
        fields = lines[line_index].split(",")
        if fields[arm_index] == arm:
            fields[field_index] = str(value)
            lines[line_index] = ",".join(fields)
    return ("\n".join(lines) + "\n").encode("ascii")


def clean_evidence(**changes):
    values = {
        "core": 4,
        "numa_node": 0,
        "exclusive_core": True,
        "load_workers_stopped": True,
        "sibling_busy_ticks": 0,
        "cpu_temperature_millic": 61_000,
        "dimm_temperature_millic": 52_000,
        "dimm_read_errors": 0,
        "edac_ce_delta": 0,
        "edac_ue_delta": 0,
        "throttled": False,
        "aperf_mperf_excursion": False,
        "telemetry_gap": False,
        "thermal_interval_sha256": THERMAL_SHA256,
        "performance_interval_sha256": PERFORMANCE_SHA256,
    }
    values.update(changes)
    return timing.EnvironmentEvidence(**values)


class EvidenceQueue:
    def __init__(self, values):
        self.values = list(values)
        self.starts = []
        self.finishes = []

    def start(self, spec, attempt_index, cycle_index):
        token = (spec.key(), attempt_index, cycle_index)
        self.starts.append(token)
        return token

    def finish(self, token):
        self.finishes.append(token)
        return self.values[len(self.finishes) - 1]


class ExecutorQueue:
    def __init__(self, values):
        self.values = list(values)
        self.calls = []

    def __call__(self, command, **kwargs):
        self.calls.append((tuple(command), kwargs))
        value = self.values[len(self.calls) - 1]
        if isinstance(value, BaseException):
            raise value
        if isinstance(value, subprocess.CompletedProcess):
            return value
        return subprocess.CompletedProcess(command, 0, stdout=value, stderr=b"")


def make_config(directory):
    binary = directory / "wh2_wirehair_test"
    binary.write_bytes(b"frozen timing binary\n")
    binary.chmod(0o755)
    route_cache = directory / "routes.json"
    route_cache.write_bytes(b'{"schema":"test"}\n')
    isolation = timing.LinuxIsolation(
        core=4, sibling_cpus=(5,), numa_node=0, llc_bytes=32 * 1024 * 1024,
        eviction_bytes=timing.MIN_EVICTION_BYTES,
    )
    numactl = directory / "numactl"
    taskset = directory / "taskset"
    for executable in (numactl, taskset):
        executable.write_bytes(b"#!/bin/sh\nexit 1\n")
        executable.chmod(0o755)
    return timing.TimingRunnerConfig(
        binary=binary, binary_sha256=sha256(binary.read_bytes()),
        route_cache=route_cache,
        route_cache_sha256=sha256(route_cache.read_bytes()),
        route_context_sha256=ROUTE_CONTEXT_SHA256,
        isolation=isolation,
        launcher=(
            str(numactl.resolve()), "--physcpubind=4", "--membind=0",
            str(taskset.resolve()), "-c", "4"),
        timeout_seconds=30,
    )


def make_systemd_config(directory):
    config = make_config(directory)
    sudo = directory / "sudo"
    systemd_run = directory / "systemd-run"
    for executable in (sudo, systemd_run):
        executable.write_bytes(b"#!/bin/sh\nexit 1\n")
        executable.chmod(0o755)
    return replace(
        config,
        launcher=(
            str(sudo.resolve()), "-n", str(systemd_run.resolve()),
            "--scope", "--quiet", "--collect",
            "--slice=wirehair-wh2-timing.slice",
            "--property=AllowedCPUs=4", "--",
            config.launcher[0], "--physcpubind=4", "--membind=0",
            config.launcher[3], "-c", "4",
        ),
        launcher_sha256=(),
    )


def mock_command(spec, route_hash, cycle_index=None):
    command = (
        "/mock/numactl", "--physcpubind=4", "--membind=0",
        "/mock/taskset", "-c", "4", "/mock/wh2_wirehair_test",
        "preferredtiming",
        "--N", str(spec.K), "--bb", str(spec.block_bytes),
        "--preferred-attempt", str(spec.preferred_attempt),
        "--evict-bytes", str(timing.MIN_EVICTION_BYTES),
        "--metric", spec.metric, "--route-cache", "/mock/routes.json",
        "--route-cache-sha256", route_hash,
        "--route-context-sha256", ROUTE_CONTEXT_SHA256,
        "--loss", "0.50", "--seed", str(spec.seed),
        "--schedule", spec.schedule,
    )
    if cycle_index is not None:
        command += ("--cycle-index", str(cycle_index))
    return command


def parsed_panel(spec, candidate_ns, control_ns=1000):
    route_hash = "a" * 64
    raw = timing_stdout(
        spec, route_cache_sha256=route_hash,
        candidate_ns=candidate_ns, control_ns=control_ns)
    parsed = timing.parse_timing_output(
        raw, spec, route_hash, ROUTE_CONTEXT_SHA256)
    command = mock_command(spec, route_hash)
    completed = subprocess.CompletedProcess(
        command, 0, stdout=raw, stderr=b"")
    contamination = {cycle: () for cycle in range(4)}
    record = timing._build_attempt_record(
        0, None, command, completed, parsed,
        clean_evidence(), contamination, timing.TIMED_CYCLES)
    record.validate()
    return timing.TimingPanelResult(
        spec=spec, trace_sha256=parsed.trace_sha256,
        canonical_attempt=parsed.canonical_attempt,
        cycles=tuple(
            timing.TimingCycle(cycle, 0, parsed.cycle_rows(cycle))
            for cycle in timing.TIMED_CYCLES),
        attempts=(record,),
    )


def parsed_panel_with_one_retry(spec, candidate_ns, control_ns=1000):
    route_hash = "a" * 64
    full_raw = timing_stdout(
        spec, route_cache_sha256=route_hash, candidate_ns=candidate_ns,
        control_ns=control_ns, cpu_by_cycle={1: (4, 5)})
    retry_raw = timing_stdout(
        spec, route_cache_sha256=route_hash, candidate_ns=candidate_ns,
        control_ns=control_ns, cycle_index=1)
    full = timing.parse_timing_output(
        full_raw, spec, route_hash, ROUTE_CONTEXT_SHA256)
    retry = timing.parse_timing_output(
        retry_raw, spec, route_hash, ROUTE_CONTEXT_SHA256, 1)
    full_command = mock_command(spec, route_hash)
    retry_command = mock_command(spec, route_hash, 1)
    full_record = timing._build_attempt_record(
        0, None, full_command,
        subprocess.CompletedProcess(
            full_command, 0, stdout=full_raw, stderr=b""),
        full, clean_evidence(),
        {0: (), 1: ("migration",), 2: (), 3: ()}, (2, 3))
    retry_record = timing._build_attempt_record(
        1, 1, retry_command,
        subprocess.CompletedProcess(
            retry_command, 0, stdout=retry_raw, stderr=b""),
        retry, clean_evidence(), {1: ()}, (1,))
    return timing.TimingPanelResult(
        spec=spec, trace_sha256=full.trace_sha256,
        canonical_attempt=full.canonical_attempt,
        cycles=(
            timing.TimingCycle(1, 1, retry.cycle_rows(1)),
            timing.TimingCycle(2, 0, full.cycle_rows(2)),
            timing.TimingCycle(3, 0, full.cycle_rows(3)),
        ),
        attempts=(full_record, retry_record),
    )


def complete_campaign(sample, solve_candidate=900, setup_candidate=1005):
    panels = []
    for K in sample:
        for width in timing.WIDTHS:
            for schedule in timing.SCHEDULES:
                value = (
                    solve_candidate(K, width, schedule)
                    if callable(solve_candidate) else solve_candidate)
                panels.append(parsed_panel(
                    panel_spec("solve", K, width, schedule), value))
            value = (
                setup_candidate(K, width)
                if callable(setup_candidate) else setup_candidate)
            panels.append(parsed_panel(
                panel_spec("setup", K, width, "burst"), value))
    return panels


class SelectionAndParsingTests(unittest.TestCase):
    def test_selection_is_deterministic_stratified_and_ascending(self):
        routed = list(range(4096, 64001, 997))
        selected = timing.select_timing_sample(routed)
        self.assertEqual(32, len(selected))
        self.assertEqual(tuple(sorted(selected)), selected)
        self.assertEqual(selected, timing.select_timing_sample(reversed(routed)))
        for low, high in timing.BANDS:
            self.assertEqual(8, sum(low <= K <= high for K in selected))
        domain = b"wirehair.wh2.h12-preferred-attempt.timing-sample.v1|"
        expected = []
        for low, high in timing.BANDS:
            expected.extend(sorted(
                (K for K in routed if low <= K <= high),
                key=lambda K: (
                    hashlib.sha256(
                        domain + str(K).encode("ascii")).digest(), K),
            )[:8])
        self.assertEqual(tuple(sorted(expected)), selected)

    def test_selection_rejects_underpowered_duplicate_and_out_of_domain(self):
        with self.assertRaises(timing.TimingError):
            timing.select_timing_sample([4096, 4097, 4098, 4099])
        with self.assertRaises(timing.TimingError):
            timing.select_timing_sample([4096, 4097, 4098, 4099, 4099])
        with self.assertRaises(timing.TimingError):
            timing.select_timing_sample([4096, 4097, 4098, 4099, 64001])

    def test_panel_spec_rejects_numeric_width_aliases(self):
        with self.assertRaisesRegex(timing.TimingError, "width"):
            replace(panel_spec(), block_bytes=64.0)

    def test_panel_spec_enforces_frozen_trace_seed(self):
        spec = panel_spec()
        payload = b"\0".join((
            timing.TRACE_SEED_DOMAIN, b"solve", b"4096", b"64", b"burst"))
        expected = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
        self.assertEqual(expected, spec.seed)
        with self.assertRaisesRegex(timing.TimingError, "seed"):
            replace(spec, seed=spec.seed ^ 1)

    def test_strict_parser_accepts_exact_solve_and_setup_fixtures(self):
        for metric in ("solve", "setup"):
            spec = panel_spec(metric)
            raw = timing_stdout(spec)
            parsed = timing.parse_timing_output(
                raw, spec, "a" * 64, ROUTE_CONTEXT_SHA256)
            self.assertEqual(32, len(parsed.rows))
            self.assertEqual(0, parsed.canonical_attempt)
            self.assertEqual(TRACE_SHA256, parsed.trace_sha256)
            self.assertEqual(sha256(raw), parsed.stdout_sha256)
            replacement = timing_stdout(spec, cycle_index=2)
            parsed_replacement = timing.parse_timing_output(
                replacement, spec, "a" * 64, ROUTE_CONTEXT_SHA256, 2)
            self.assertEqual(8, len(parsed_replacement.rows))
            self.assertEqual((2,), tuple(sorted(set(
                row.cycle for row in parsed_replacement.rows))))

    def test_strict_parser_rejects_contract_mutations(self):
        spec = panel_spec()
        valid = timing_stdout(spec)
        mutations = (
            valid.replace(b"cycles=4", b"cycles=04", 1),
            valid.replace(b"cycle,slot", b"slot,cycle", 1),
            mutate_csv_field(valid, 0, "result", 1),
            mutate_csv_field(valid, 1, "attempt", 1),
            mutate_csv_field(valid, 0, "source_bytes", spec.K * 64 - 64),
            mutate_csv_field(
                valid, 0, "intermediate_bytes", spec.K * spec.block_bytes),
            mutate_csv_field(valid, 0, "intermediate_bytes", spec.K * 64 + 1),
            valid.replace(TRACE_SHA256.encode("ascii"), b"C" * 64, 1),
            valid[:-1],
            valid + b"extra\n",
        )
        for raw in mutations:
            with self.subTest(raw=raw[:100]):
                with self.assertRaises(timing.TimingError):
                    timing.parse_timing_output(
                        raw, spec, "a" * 64, ROUTE_CONTEXT_SHA256)
        replacement = timing_stdout(spec, cycle_index=2)
        with self.assertRaises(timing.TimingError):
            timing.parse_timing_output(
                replacement, spec, "a" * 64, ROUTE_CONTEXT_SHA256, 1)
        with self.assertRaises(timing.TimingError):
            timing.parse_timing_output(
                replacement, spec, "a" * 64, ROUTE_CONTEXT_SHA256)

    def test_wide_neutral_alias_is_allowed_but_bb64_must_remain_direct(self):
        for width in (1280, 4096):
            spec = panel_spec(width=width)
            raw = timing_stdout(
                spec, canonical_attempt=spec.preferred_attempt)
            parsed = timing.parse_timing_output(
                raw, spec, "a" * 64, ROUTE_CONTEXT_SHA256)
            self.assertEqual(spec.preferred_attempt, parsed.canonical_attempt)
            self.assertEqual(1, len({
                row.work_signature() for row in parsed.rows}))
        spec64 = panel_spec(width=64)
        raw64 = timing_stdout(
            spec64, canonical_attempt=spec64.preferred_attempt)
        with self.assertRaisesRegex(timing.TimingError, "eligible"):
            timing.parse_timing_output(
                raw64, spec64, "a" * 64, ROUTE_CONTEXT_SHA256)

    def test_wide_zero_attempt_alias_matches_native_cli_domain(self):
        spec = replace(panel_spec(width=1280), preferred_attempt=0)
        parsed = timing.parse_timing_output(
            timing_stdout(spec, canonical_attempt=0), spec,
            "a" * 64, ROUTE_CONTEXT_SHA256)
        self.assertEqual(0, parsed.canonical_attempt)
        with self.assertRaisesRegex(timing.TimingError, "eligible domain"):
            replace(panel_spec(width=64), preferred_attempt=0)

    def test_wide_neutral_alias_rejects_different_work_between_labels(self):
        spec = panel_spec(width=1280)
        raw = timing_stdout(
            spec, canonical_attempt=spec.preferred_attempt)
        raw = mutate_arm_field(raw, "candidate", "block_xors", 999)
        with self.assertRaisesRegex(timing.TimingError, "neutral alias"):
            timing.parse_timing_output(
                raw, spec, "a" * 64, ROUTE_CONTEXT_SHA256)

    def test_solve_intermediate_domain_must_exceed_source_domain(self):
        spec = panel_spec()
        with self.assertRaisesRegex(timing.TimingError, "byte count"):
            timing.parse_timing_output(
                timing_stdout(spec, intermediate_blocks=0), spec,
                "a" * 64, ROUTE_CONTEXT_SHA256)


class IsolationTests(unittest.TestCase):
    def test_linux_isolation_discovers_siblings_numa_llc_and_eviction(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            cpu = root / "cpu4"
            (cpu / "topology").mkdir(parents=True)
            (cpu / "topology/thread_siblings_list").write_text(
                "4-5\n", encoding="ascii")
            (cpu / "node0").mkdir()
            for index, level, cache_type, size, shared in (
                (0, "1", "Data", "32K", "4"),
                (1, "3", "Unified", "32M", "0-7"),
            ):
                cache = cpu / "cache" / ("index%d" % index)
                cache.mkdir(parents=True)
                (cache / "level").write_text(level + "\n", encoding="ascii")
                (cache / "type").write_text(cache_type + "\n", encoding="ascii")
                (cache / "size").write_text(size + "\n", encoding="ascii")
                (cache / "shared_cpu_list").write_text(
                    shared + "\n", encoding="ascii")
            result = timing.inspect_linux_isolation(4, root)
            self.assertEqual((5,), result.sibling_cpus)
            self.assertEqual(0, result.numa_node)
            self.assertEqual(32 * 1024 * 1024, result.llc_bytes)
            self.assertEqual(256 * 1024 * 1024, result.eviction_bytes)

    def test_runner_config_rejects_nonfrozen_eviction_size(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            with self.assertRaises(timing.TimingError):
                replace(config.isolation, eviction_bytes=123)

    def test_runner_config_rejects_missing_isolation_launcher(self):
        with tempfile.TemporaryDirectory() as name:
            config = make_config(Path(name))
            with self.assertRaisesRegex(timing.TimingError, "schema"):
                replace(config, launcher=())

    def test_systemd_launcher_schema_and_hashes_are_frozen(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_systemd_config(directory)
            self.assertEqual(4, len(config.launcher_sha256))
            Path(config.launcher[2]).write_bytes(b"changed systemd-run\n")
            with self.assertRaisesRegex(timing.TimingError, "not frozen"):
                timing.run_timing_panel(
                    panel_spec(), config, EvidenceQueue([clean_evidence()]),
                    directory / "attempts", executor=ExecutorQueue([]))

    def test_systemd_launcher_round_trips_through_panel_seal(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_systemd_config(directory)
            spec = panel_spec()
            executor = ExecutorQueue([
                timing_stdout(spec, config.route_cache_sha256)])
            attempts = directory / "attempts"
            result = timing.run_timing_panel(
                spec, config, EvidenceQueue([clean_evidence()]), attempts,
                executor=executor)
            self.assertEqual(result, timing.load_timing_panel_result(
                attempts, spec, config))
            command = executor.calls[0][0]
            self.assertEqual("sudo", Path(command[0]).name)
            self.assertEqual("systemd-run", Path(command[2]).name)
            self.assertEqual("numactl", Path(command[9]).name)
            self.assertEqual(
                str(config.binary.resolve()), command[len(config.launcher)])

    def test_performance_interval_hash_is_mandatory(self):
        with self.assertRaisesRegex(timing.TimingError, "performance interval"):
            clean_evidence(performance_interval_sha256="missing")


class RunnerTests(unittest.TestCase):
    def test_default_executor_reaps_timeout_process_group(self):
        started = time.monotonic()
        with self.assertRaises(subprocess.TimeoutExpired):
            timing._execute_timing_process(
                ("/bin/sh", "-c", "sleep 5 & wait"),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=0.05)
        self.assertLess(time.monotonic() - started, 1.0)

    def test_probe_start_failure_is_retained_as_unsealed_evidence(self):
        class FailingProbe:
            def start(self, _spec, _attempt_index, _cycle_index):
                raise RuntimeError("probe unavailable")

        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            attempts = directory / "attempts"
            with self.assertRaisesRegex(timing.TimingError, "failed to start"):
                timing.run_timing_panel(
                    panel_spec(), config, FailingProbe(), attempts,
                    executor=ExecutorQueue([]))
            payload = json.loads(
                (attempts / "attempt-00.json").read_text())
            self.assertIn("probe start failed", payload["failure"])
            self.assertFalse((attempts / timing.PANEL_RESULT_NAME).exists())

    def test_clean_panel_uses_one_attempt_and_persists_hashed_evidence(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            executor = ExecutorQueue([
                timing_stdout(spec, config.route_cache_sha256)])
            probe = EvidenceQueue([clean_evidence()])
            attempts = directory / "attempts"
            result = timing.run_timing_panel(
                spec, config, probe, attempts, executor=executor)
            self.assertEqual((0, 0, 0), tuple(
                cycle.attempt_index for cycle in result.cycles))
            self.assertEqual(1, len(result.attempts))
            self.assertEqual(Fraction(9, 10), result.panel_ratio())
            self.assertEqual(5, len(list(attempts.iterdir())))
            payload = json.loads((attempts / "attempt-00.json").read_text())
            record_hash = payload.pop("record_sha256")
            self.assertEqual(record_hash, sha256(canonical_json(payload)))
            self.assertEqual(
                payload["stdout_sha256"],
                sha256((attempts / "attempt-00.stdout").read_bytes()))
            self.assertEqual(
                PERFORMANCE_SHA256,
                payload["environment"]["performance_interval_sha256"])
            command, kwargs = executor.calls[0]
            self.assertEqual(
                str(config.binary.resolve()), command[len(config.launcher)])
            self.assertIn("--evict-bytes", command)
            self.assertIn(str(timing.MIN_EVICTION_BYTES), command)
            self.assertEqual(30, kwargs["timeout"])
            panel_raw = (attempts / timing.PANEL_RESULT_NAME).read_bytes()
            panel = json.loads(panel_raw)
            self.assertEqual(
                config.binary_sha256,
                panel["runtime_binding"]["binary_sha256"])
            self.assertEqual(
                config.route_cache_sha256,
                panel["runtime_binding"]["route_cache_sha256"])
            self_hash = panel.pop("self_sha256_excluding_field")
            self.assertEqual(self_hash, sha256(canonical_json(panel)))
            self.assertEqual(
                (sha256(panel_raw) + "  panel.json\n").encode("ascii"),
                (attempts / timing.PANEL_RESULT_SIDECAR_NAME).read_bytes())

    def test_run_or_resume_reconstructs_clean_panel_without_invocation(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            attempts = directory / "attempts"
            original = timing.run_timing_panel(
                spec, config, EvidenceQueue([clean_evidence()]), attempts,
                executor=ExecutorQueue([
                    timing_stdout(spec, config.route_cache_sha256)]))
            executor = ExecutorQueue([])
            probe = EvidenceQueue([])
            resumed = timing.run_or_resume_timing_panel(
                spec, config, probe, attempts, executor=executor)
            self.assertEqual(original, resumed)
            self.assertEqual(Fraction(9, 10), resumed.panel_ratio())
            self.assertEqual([], executor.calls)
            self.assertEqual([], probe.starts)

    def test_failed_or_partial_panel_is_never_rerun_by_controller(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            attempts = directory / "attempts"
            with self.assertRaises(timing.TimingError):
                timing.run_timing_panel(
                    spec, config, EvidenceQueue([clean_evidence()]), attempts,
                    executor=ExecutorQueue([b"malformed\n"]))
            executor = ExecutorQueue([
                timing_stdout(spec, config.route_cache_sha256)])
            with self.assertRaisesRegex(
                    timing.TimingError, "controller-level rerun is forbidden"):
                timing.run_or_resume_timing_panel(
                    spec, config, EvidenceQueue([clean_evidence()]), attempts,
                    executor=executor)
            self.assertEqual([], executor.calls)

    def test_resume_rejects_tampered_retained_stdout(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            attempts = directory / "attempts"
            timing.run_timing_panel(
                spec, config, EvidenceQueue([clean_evidence()]), attempts,
                executor=ExecutorQueue([
                    timing_stdout(spec, config.route_cache_sha256)]))
            stdout = attempts / "attempt-00.stdout"
            stdout.write_bytes(stdout.read_bytes() + b"tamper\n")
            with self.assertRaisesRegex(timing.TimingError, "stream hash"):
                timing.load_timing_panel_result(attempts, spec, config)

    def test_resume_recomputes_contamination_from_resealed_stdout(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            attempts = directory / "attempts"
            timing.run_timing_panel(
                spec, config, EvidenceQueue([clean_evidence()]), attempts,
                executor=ExecutorQueue([
                    timing_stdout(spec, config.route_cache_sha256)]))

            stdout_path = attempts / "attempt-00.stdout"
            stdout = mutate_csv_field(
                stdout_path.read_bytes(), 0, "cpu_after", 5)
            stdout_path.write_bytes(stdout)
            attempt_path = attempts / "attempt-00.json"
            attempt = json.loads(attempt_path.read_text())
            attempt["stdout_sha256"] = sha256(stdout)
            attempt.pop("record_sha256")
            attempt["record_sha256"] = sha256(canonical_json(attempt))
            attempt_path.write_bytes(canonical_json(attempt))

            panel_path = attempts / timing.PANEL_RESULT_NAME
            panel = json.loads(panel_path.read_text())
            panel["attempt_record_sha256"][0] = attempt["record_sha256"]
            panel.pop("self_sha256_excluding_field")
            panel["self_sha256_excluding_field"] = sha256(
                canonical_json(panel))
            panel_raw = canonical_json(panel)
            panel_path.write_bytes(panel_raw)
            (attempts / timing.PANEL_RESULT_SIDECAR_NAME).write_bytes(
                (sha256(panel_raw) + "  panel.json\n").encode("ascii"))
            with self.assertRaisesRegex(
                    timing.TimingError, "contamination changed"):
                timing.load_timing_panel_result(attempts, spec, config)

    def test_retry_replaces_only_contaminated_cycles(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            first = timing_stdout(
                spec, config.route_cache_sha256,
                cpu_by_cycle={1: (4, 5)})
            second = timing_stdout(
                spec, config.route_cache_sha256, cycle_index=1)
            executor = ExecutorQueue([first, second])
            probe = EvidenceQueue([clean_evidence(), clean_evidence()])
            result = timing.run_timing_panel(
                spec, config, probe, directory / "attempts", executor=executor)
            self.assertEqual((1, 0, 0), tuple(
                cycle.attempt_index for cycle in result.cycles))
            self.assertEqual((2, 3), result.attempts[0].accepted_cycles)
            self.assertEqual((1,), result.attempts[1].accepted_cycles)
            self.assertEqual(("migration",),
                             result.attempts[0].cycle_contamination[1])
            self.assertEqual((None, 1), tuple(
                attempt.cycle_index for attempt in result.attempts))
            self.assertNotIn("--cycle-index", executor.calls[0][0])
            self.assertEqual(("--cycle-index", "1"),
                             executor.calls[1][0][-2:])

    def test_environment_and_row_contamination_get_two_bounded_retries(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            full = timing_stdout(spec, config.route_cache_sha256)
            cycle_one = timing_stdout(
                spec, config.route_cache_sha256, cycle_index=1)
            cycle_two_fault = timing_stdout(
                spec, config.route_cache_sha256,
                major_fault_by_cycle={2: 1}, cycle_index=2)
            cycle_three = timing_stdout(
                spec, config.route_cache_sha256, cycle_index=3)
            cycle_two = timing_stdout(
                spec, config.route_cache_sha256, cycle_index=2)
            executor = ExecutorQueue([
                full, cycle_one, cycle_two_fault, cycle_three, cycle_two])
            probe = EvidenceQueue([
                clean_evidence(cpu_temperature_millic=85_000),
                clean_evidence(), clean_evidence(), clean_evidence(),
                clean_evidence(),
            ])
            result = timing.run_timing_panel(
                spec, config, probe, directory / "attempts", executor=executor)
            self.assertEqual((1, 4, 3), tuple(
                cycle.attempt_index for cycle in result.cycles))
            self.assertEqual((None, 1, 2, 3, 2), tuple(
                attempt.cycle_index for attempt in result.attempts))
            self.assertEqual(5, len(executor.calls))
            self.assertEqual(17, len(list((directory / "attempts").iterdir())))

    def test_retry_exhaustion_retains_full_plus_six_replacements(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            outputs = [timing_stdout(spec, config.route_cache_sha256)]
            outputs.extend(
                timing_stdout(
                    spec, config.route_cache_sha256, cycle_index=cycle)
                for _round in range(2) for cycle in timing.TIMED_CYCLES)
            executor = ExecutorQueue(outputs)
            probe = EvidenceQueue([
                clean_evidence(telemetry_gap=True) for _value in outputs])
            with self.assertRaisesRegex(timing.TimingError, "two retries"):
                timing.run_timing_panel(
                    spec, config, probe, directory / "attempts",
                    executor=executor)
            self.assertEqual(7, len(executor.calls))
            self.assertEqual(21, len(list((directory / "attempts").iterdir())))

    def test_replacement_work_signature_must_match_full_attempt(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            full = timing_stdout(
                spec, config.route_cache_sha256,
                cpu_by_cycle={1: (4, 5)})
            replacement = timing_stdout(
                spec, config.route_cache_sha256, cycle_index=1,
                intermediate_blocks=4)
            with self.assertRaisesRegex(timing.TimingError, "work changed"):
                timing.run_timing_panel(
                    spec, config,
                    EvidenceQueue([clean_evidence(), clean_evidence()]),
                    directory / "attempts",
                    executor=ExecutorQueue([full, replacement]))
            self.assertEqual(6, len(list((directory / "attempts").iterdir())))

    def test_edac_change_is_fatal_and_retained(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            executor = ExecutorQueue([
                timing_stdout(spec, config.route_cache_sha256)])
            probe = EvidenceQueue([clean_evidence(edac_ce_delta=1)])
            with self.assertRaisesRegex(timing.TimingError, "EDAC"):
                timing.run_timing_panel(
                    spec, config, probe, directory / "attempts",
                    executor=executor)
            self.assertEqual(1, len(executor.calls))
            self.assertEqual(3, len(list((directory / "attempts").iterdir())))
            payload = json.loads(
                (directory / "attempts/attempt-00.json").read_text())
            self.assertEqual(
                "wirehair.wh2.h12_preferred_attempt.failed_attempt.v1",
                payload["schema"])

    def test_isolation_loss_is_fatal_not_retryable_contamination(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            executor = ExecutorQueue([
                timing_stdout(spec, config.route_cache_sha256)])
            probe = EvidenceQueue([
                clean_evidence(load_workers_stopped=False)])
            with self.assertRaisesRegex(timing.TimingError, "isolation failed"):
                timing.run_timing_panel(
                    spec, config, probe, directory / "attempts",
                    executor=executor)
            self.assertEqual(1, len(executor.calls))
            payload = json.loads(
                (directory / "attempts/attempt-00.json").read_text())
            self.assertEqual(
                "wirehair.wh2.h12_preferred_attempt.failed_attempt.v1",
                payload["schema"])

    def test_parse_failure_and_timeout_are_retained_as_failed_attempts(self):
        for value in (
                b"malformed\n",
                subprocess.TimeoutExpired(("mock",), 30, output=b"partial")):
            with self.subTest(value=type(value).__name__):
                with tempfile.TemporaryDirectory() as name:
                    directory = Path(name)
                    config = make_config(directory)
                    spec = panel_spec()
                    executor = ExecutorQueue([value])
                    probe = EvidenceQueue([clean_evidence()])
                    with self.assertRaises(timing.TimingError):
                        timing.run_timing_panel(
                            spec, config, probe, directory / "attempts",
                            executor=executor)
                    payload = json.loads(
                        (directory / "attempts/attempt-00.json").read_text())
                    self.assertEqual(
                        "wirehair.wh2.h12_preferred_attempt.failed_attempt.v1",
                        payload["schema"])
                    self.assertEqual(3, len(list(
                        (directory / "attempts").iterdir())))

    def test_symlink_input_is_rejected_before_attempt_directory_creation(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            link = directory / "binary-link"
            link.symlink_to(config.binary)
            config = replace(config, binary=link)
            attempts = directory / "attempts"
            with self.assertRaisesRegex(timing.TimingError, "symlink"):
                timing.run_timing_panel(
                    panel_spec(), config, EvidenceQueue([clean_evidence()]),
                    attempts, executor=ExecutorQueue([]))
            self.assertFalse(attempts.exists())

    def test_route_cache_change_during_process_is_fatal_and_retained(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            config = make_config(directory)
            spec = panel_spec()
            raw = timing_stdout(spec, config.route_cache_sha256)

            def mutate_cache(command, **_kwargs):
                config.route_cache.write_bytes(b"changed route cache\n")
                return subprocess.CompletedProcess(
                    command, 0, stdout=raw, stderr=b"")

            attempts = directory / "attempts"
            with self.assertRaisesRegex(timing.TimingError, "changed during"):
                timing.run_timing_panel(
                    spec, config, EvidenceQueue([clean_evidence()]),
                    attempts, executor=mutate_cache)
            payload = json.loads((attempts / "attempt-00.json").read_text())
            self.assertEqual(
                "wirehair.wh2.h12_preferred_attempt.failed_attempt.v1",
                payload["schema"])


class AnalysisTests(unittest.TestCase):
    SAMPLE = (4096, 8192, 16384, 32768, 64000)

    def test_exact_medians_sign_gate_thresholds_and_realized_counts_pass(self):
        report = timing.analyze_timing(
            self.SAMPLE, complete_campaign(self.SAMPLE))
        self.assertTrue(report["accepted"])
        self.assertEqual(
            {"numerator": 9, "denominator": 10},
            report["solve"]["global_ratio"])
        self.assertEqual(
            {"numerator": 201, "denominator": 200},
            report["setup"]["global_ratio"])
        self.assertEqual(5, report["solve"]["sign"]["wins"])
        self.assertEqual(1, report["solve"]["sign"]["tail_numerator"])
        self.assertEqual(32, report["solve"]["sign"]["denominator"])
        self.assertEqual(45, report["realized"]["solve_panels"])
        self.assertEqual(15, report["realized"]["setup_panels"])
        self.assertEqual(1920, report["realized"]["accepted_total_invocations"])
        self.assertEqual(1440, report["realized"]["accepted_timed_invocations"])
        self.assertEqual(60, report["realized"][
            "physical_process_attempts_including_retries"])

    def test_physical_accounting_counts_replacement_as_eight_invocations(self):
        panels = complete_campaign(self.SAMPLE)
        panels[0] = parsed_panel_with_one_retry(panels[0].spec, 900)
        report = timing.analyze_timing(self.SAMPLE, panels)
        self.assertTrue(report["accepted"])
        self.assertEqual(60, report["realized"][
            "physical_full_process_attempts"])
        self.assertEqual(1, report["realized"][
            "physical_replacement_process_attempts"])
        self.assertEqual(1928, report["realized"][
            "physical_invocations_including_retries"])

    def test_panel_validation_rejects_repeated_clean_full_run(self):
        panel = parsed_panel(panel_spec(), 900)
        duplicate = replace(
            panel.attempts[0], attempt_index=1, accepted_cycles=(),
            record_sha256="0" * 64)
        duplicate = replace(
            duplicate,
            record_sha256=sha256(canonical_json(duplicate.core_record())))
        repeated = replace(panel, attempts=(panel.attempts[0], duplicate))
        with self.assertRaisesRegex(timing.TimingError, "repeated a full run"):
            repeated.validate()

    def test_sign_gate_rejects_four_wins_and_one_loss_despite_fast_median(self):
        slow_K = self.SAMPLE[-1]
        panels = complete_campaign(
            self.SAMPLE,
            solve_candidate=lambda K, _width, _schedule: (
                1100 if K == slow_K else 900))
        report = timing.analyze_timing(self.SAMPLE, panels)
        self.assertFalse(report["accepted"])
        self.assertTrue(report["solve"]["global_strictly_below_one"])
        self.assertEqual(4, report["solve"]["sign"]["wins"])
        self.assertEqual(1, report["solve"]["sign"]["losses"])
        self.assertEqual(6, report["solve"]["sign"]["tail_numerator"])
        self.assertFalse(report["solve"]["sign"][
            "twenty_tail_le_denominator"])

    def test_setup_per_width_gate_rejects_1_006_even_if_global_is_1_005(self):
        panels = complete_campaign(
            self.SAMPLE,
            setup_candidate=lambda _K, width: (
                1006 if width == 4096 else 1005))
        report = timing.analyze_timing(self.SAMPLE, panels)
        self.assertFalse(report["accepted"])
        self.assertTrue(report["setup"]["global_at_most_1_005"])
        self.assertFalse(report["setup"]["per_width_at_most_1_005"])

    def test_analysis_rejects_missing_duplicate_and_noncanonical_coverage(self):
        panels = complete_campaign(self.SAMPLE)
        with self.assertRaises(timing.TimingError):
            timing.analyze_timing(self.SAMPLE, panels[:-1])
        with self.assertRaises(timing.TimingError):
            timing.analyze_timing(self.SAMPLE, panels + [panels[0]])
        with self.assertRaises(timing.TimingError):
            timing.analyze_timing(tuple(reversed(self.SAMPLE)), panels)

    def test_analysis_rejects_one_K_using_different_attempt_across_widths(self):
        panels = complete_campaign(self.SAMPLE)
        changed = []
        for panel in panels:
            if panel.spec.K == self.SAMPLE[0] and panel.spec.block_bytes == 1280:
                spec = replace(panel.spec, preferred_attempt=3)
                candidate_ns = 900 if spec.metric == "solve" else 1005
                changed.append(parsed_panel(spec, candidate_ns))
            else:
                changed.append(panel)
        with self.assertRaisesRegex(timing.TimingError, "across widths"):
            timing.analyze_timing(self.SAMPLE, changed)


if __name__ == "__main__":
    unittest.main()
