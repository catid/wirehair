#!/usr/bin/env python3
"""Bounded tests for the fixed native WH2 short-screen controller."""

from __future__ import annotations

import copy
from contextlib import redirect_stderr
import csv
import hashlib
from io import StringIO
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any, Mapping
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
    {"arm": "wirehair2_identity", "codec": "wirehair2_experiment",
     "arm_descriptor_sha256": digest("identity")},
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


class NativeRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

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
        for invalid in ("", "00", "0,0", "1,0", "0, 1", "+1"):
            with self.assertRaises(subject.RunnerError):
                subject.parse_cpu_list(invalid)
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

    def test_freezes_are_published_before_result_sink(self) -> None:
        worker = self._fake_worker()
        description = self._description(worker)
        contract = contract_api.load_contract()
        output = self.root / "output"
        output.mkdir()
        trace_hashes = {
            "recovery": hashlib.sha256(b"recovery trace").hexdigest(),
            "timing": hashlib.sha256(b"timing trace").hexdigest(),
        }
        freezes = subject.write_development_freezes(
            contract, description, [0], 1, "1" * 40, trace_hashes, output)
        for kind in ("recovery", "timing"):
            path = output / "{}-freeze.json".format(kind)
            self.assertTrue(path.is_file())
            self.assertEqual(freezes[kind]["trace_manifest_sha256"],
                             trace_hashes[kind])
            self.assertEqual(freezes[kind]["arm_roster"], [
                "wirehair2_head", "wirehair1", "wirehair2_identity"])
            self.assertEqual(freezes[kind]["host_identity"]["controller_cpu"],
                             1)
            self.assertEqual(
                [value["construction_policy"]
                 for value in freezes[kind]["arms"]],
                ["raw_base", "not_applicable", "raw_base"])
        sink = subject.AtomicLineSink(output / "later-native-results.jsonl")
        try:
            self.assertTrue(all(
                (output / "{}-freeze.json".format(kind)).exists()
                for kind in ("recovery", "timing")))
        finally:
            sink.abort()

    def test_strict_timing_validator_precomputes_cell_indexes_once(self) -> None:
        worker = self._fake_worker()
        description = self._description(worker)
        contract = contract_api.load_contract()
        output = self.root / "output"
        output.mkdir()
        freezes = subject.write_development_freezes(
            contract, description, [0], 1, "1" * 40,
            {"recovery": "2" * 64, "timing": "3" * 64}, output)
        with mock.patch.object(
                contract_api, "_timing_cell_indexes",
                wraps=contract_api._timing_cell_indexes) as cell_indexes:
            validator = subject._strict_response_validator(
                contract, freezes["timing"], "timing", description, 0)
            self.assertTrue(callable(validator))
            self.assertEqual(cell_indexes.call_count, 1)

    def test_post_link_publish_failure_removes_own_artifact(self) -> None:
        destination = self.root / "atomic-output.json"
        with mock.patch.object(
                subject, "_fsync_directory",
                side_effect=OSError("injected directory fsync failure")):
            with self.assertRaises(subject.RunnerError):
                subject._atomic_write_bytes(destination, b"complete\n")
        self.assertFalse(destination.exists())

    def test_post_link_deadline_removes_own_artifact(self) -> None:
        destination = self.root / "deadline-output.json"
        error = subject.RunnerError("injected hard-wall deadline")
        with mock.patch.object(
                subject, "_fsync_directory", side_effect=error):
            with self.assertRaisesRegex(subject.RunnerError, "hard-wall"):
                subject._atomic_write_bytes(destination, b"complete\n")
        self.assertFalse(destination.exists())

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
                waves = subject._timing_job_waves(candidate, len(panels))
                cells = list(contract_api.iter_timing_cells(
                    candidate, "development"))
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

                rebuilt = subject._timing_job_waves(candidate, len(panels))
                self.assertEqual(
                    [(rotation, [job.command() for job in jobs])
                     for rotation, jobs in waves],
                    [(rotation, [job.command() for job in jobs])
                     for rotation, jobs in rebuilt])

    def test_timing_wave_builder_fails_closed_on_domain_mutations(self) -> None:
        contract = self._timing_contract(16)
        panel_count = len(contract_api.timing_panels(
            contract, [value[0] for value in subject.EXPECTED_ARMS]))
        original = list(contract_api.iter_timing_cells(
            contract, "development"))
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
                    contract_api, "iter_timing_cells",
                    side_effect=lambda _contract, _phase, values=cells:
                        iter(values)):
                with self.assertRaises(subject.RunnerError):
                    subject._timing_job_waves(contract, panel_count)
        for bad_panel_count in (True, panel_count - 1, panel_count + 1):
            with self.subTest(panel_count=bad_panel_count):
                with self.assertRaises(subject.RunnerError):
                    subject._timing_job_waves(contract, bad_panel_count)

    def test_timing_wave_builder_derives_multi_candidate_cohort_count(self) \
            -> None:
        contract = self._timing_contract(16)
        expanded_arms = subject.EXPECTED_ARMS + (
            ("candidate_two", "wirehair2_experiment"),
        )
        with mock.patch.object(subject, "EXPECTED_ARMS", expanded_arms):
            roster = [value[0] for value in expanded_arms]
            panel_count = len(contract_api.timing_panels(contract, roster))
            self.assertEqual(panel_count, 18)
            waves = subject._timing_job_waves(contract, panel_count)
            self.assertEqual(
                len(waves),
                contract_api.timing_cohort_count(
                    contract, "development", roster) * 2)
            self.assertEqual(
                sum(len(jobs) for _, jobs in waves), 384 * panel_count)

    def test_timing_worker_count_rejects_before_native_dispatch(self) -> None:
        contract = self._timing_contract(16)
        for count in (0, 7, 9):
            with self.subTest(count=count), self.assertRaisesRegex(
                    subject.RunnerError, "exactly eight timing workers"):
                subject._run_native_jobs(
                    contract, {}, {}, [object()] * count,
                    self.root, 0, time.monotonic() + 1.0)
        for repetitions in (0, 7, 12):
            candidate = self._timing_contract(16)
            domain = candidate["timing"]["domains"]["development"]
            domain["paired_repetitions"] = repetitions
            with self.subTest(repetitions=repetitions), \
                    self.assertRaisesRegex(
                        subject.RunnerError, "positive multiple of eight"):
                subject._timing_job_waves(candidate, 11)

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
