#!/usr/bin/env python3
"""Lightweight tests for the fixed precodefail work/rank sidecar."""

from __future__ import annotations

import concurrent.futures
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_precodefail_work_screen as subject


FAKE_WORKER = r'''#!/usr/bin/env python3
import os
import sys

HEADER = (
    "N,bb,heavy_family,mix_count,staircase,binary_dense_rows,"
    "gf256_heavy_rows,source_hits,dense_identity_corner,"
    "overhead,trials,success,rank_fail,error,"
    "fail_rate,inact_mu,inact_max,binary_def_mu,binary_def_max,heavy_gain_mu,"
    "heavy_gain_min,heavy_shortfall,solve_ms_mu,build_ms_mu,peel_ms_mu,"
    "project_ms_mu,residual_ms_mu,backsub_ms_mu,seed_attempt,block_xors_mu,"
    "block_muladds_mu,first_rank_fail,binary_def_hist,heavy_gain_hist,"
    "failure_trials,active_packet_peel_seed_xor,precode_attempt,"
    "packet_attempt,attempt_mode,construction_seed_basis,"
    "seed_schedule_sha256,effective_precode_seed,effective_packet_seed"
)

RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_SHA256 = (
    "90a98a3db207852dabdf5fb27573ef48b"
    "ce52e0228cee4e291d96fa44ed509a7"
)
RAW_PRECODE_BASE_SEED = 0x487468302aad7105
RAW_PACKET_BASE_SEED = 0x4ec72102
RAW_PRECODE_ATTEMPT_STRIDE = 0x9e3779b97f4a7c15
RAW_PACKET_ATTEMPT_STRIDE = 0x9e3779b9
STAIRCASE_BY_K = {
    2: 2, 3: 3, 4: 3, 5: 5, 6: 6, 8: 6, 16: 13, 32: 13, 64: 19,
    100: 26, 101: 26, 128: 30, 256: 30, 512: 38, 513: 38,
    1000: 50, 1001: 50, 2048: 54, 4096: 70, 5000: 66, 5001: 66,
    8192: 78, 10000: 86, 10001: 86, 16384: 114, 20000: 134,
    20001: 134, 32768: 194, 49152: 378, 64000: 346,
}

args = sys.argv[1:]
if not args or args[0] != "precodefail":
    raise SystemExit(2)
flags = set(("--paired-overhead-stream", "--full-payload-solve"))
values = {}
i = 1
while i < len(args):
    token = args[i]
    if token in flags:
        values[token] = True
        i += 1
    else:
        values[token] = args[i + 1]
        i += 2

trial = int(values["--exact-precode-attempt"])
schedule = values["--schedule"]
loss = float(values["--loss"])
dense = int(values["--binary-dense-rows"])
heavy = int(values["--gf256-heavy-rows"])
dense_anchor_layout = values["--dense-anchors"]
seed_basis = values["--construction-seed-basis"]
if seed_basis != RAW_SEED_BASIS:
    raise SystemExit(2)
effective_precode_seed = (
    RAW_PRECODE_BASE_SEED + trial * RAW_PRECODE_ATTEMPT_STRIDE) & ((1 << 64) - 1)
effective_packet_seed = (
    RAW_PACKET_BASE_SEED + trial * RAW_PACKET_ATTEMPT_STRIDE) & ((1 << 32) - 1)
metadata = [
    "trials=1", "threads=1", "loss={:.17g}".format(loss),
    "seed={}".format(values["--seed"]), "source_hits_override=0",
    "packet_peel_seed_xor=0x0",
    "binary_dense_rows_override={}".format(dense),
    "gf256_heavy_rows_override={}".format(heavy),
    "dense_anchor_layout={}".format(dense_anchor_layout),
    "odd_packet_peel_seed_xor=0x0", "packet_row_seed_multiplier=0x1",
    "packet_row_seed_avalanche=0", "seed_block_bytes_override=0",
    "overhead_stream=paired", "full_payload_solve=1",
    "schedule={}".format(schedule), "cold_solve_wide_xor=policy",
    "exact_attempt_mode=1",
    "exact_precode_attempt={}".format(trial),
    "exact_packet_attempt={}".format(trial),
    "construction_seed_basis={}".format(seed_basis),
    "seed_schedule_sha256={}".format(RAW_SEED_SCHEDULE_SHA256),
    "source_git_commit=" + "1" * 40,
]
corruption = os.environ.get("WH2_FAKE_PRECODEFAIL_CORRUPTION", "")
if corruption == "stdout_overflow" and schedule == "iid" and trial == 0 and \
        dense_anchor_layout == "disabled":
    import time
    sys.stdout.write("x" * (2 * 1024 * 1024 + 1))
    sys.stdout.flush()
    time.sleep(10)
    raise SystemExit(9)
if corruption == "metadata" and schedule == "iid" and trial == 0 and \
        dense_anchor_layout == "disabled":
    metadata[0] = "trials=2"
print("# precodefail: " + " ".join(metadata))
print(HEADER)
schedule_index = ("iid", "burst", "adversarial", "repair-only").index(schedule)
for K in [int(value) for value in values["--N"].split(",")]:
    for overhead in [int(value) for value in values["--overhead"].split(",")]:
        inact = 4 + (K + trial + schedule_index) % 17
        binary_def = (K + dense + trial) % 4
        rank_fail = int(binary_def > 0 and
                        (K + trial + schedule_index + dense + heavy) %
                        (7 + overhead) == 0)
        gain = binary_def - rank_fail
        success = 1 - rank_fail
        shortfall = int(rank_fail and binary_def <= heavy and gain < binary_def)
        block_xors = K * 3 + dense + overhead
        muladds = K + heavy * 2 + overhead
        staircase = STAIRCASE_BY_K[K]
        source_hits = 3 if K >= 10000 else 2
        row = [
            str(K), "2", "periodic", "3", str(staircase), str(dense),
            str(heavy), str(source_hits), "0", str(overhead), "1",
            str(success), str(rank_fail), "0",
            "{:.8f}".format(float(rank_fail)), "{}.000".format(inact),
            str(inact), "{}.000".format(binary_def), str(binary_def),
            "{}.000".format(gain), str(gain), str(shortfall),
            "0.010", "0.001", "0.002", "0.003", "0.002", "0.002", "",
            "{}.000".format(block_xors), "{}.000".format(muladds),
            "0" if rank_fail else "-1", "{}:1".format(binary_def),
            "{}:1".format(gain), "0" if rank_fail else "", "0x0",
            str(trial), str(trial), "exact", seed_basis,
            RAW_SEED_SCHEDULE_SHA256,
            "0x{:016x}".format(effective_precode_seed),
            "0x{:08x}".format(effective_packet_seed),
        ]
        if corruption == "error" and schedule == "iid" and trial == 0 and \
                dense_anchor_layout == "disabled" and K == 2 and overhead == 0:
            row[11:15] = ["0", "0", "1", "1.00000000"]
        print(",".join(row))
'''


class PrecodefailWorkScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = contract_api.load_contract()
        cls.cell_map, cls.trace_map, cls.proof = subject.build_frozen_domain(
            cls.contract)

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.worker = self.root / "fake_precodefail.py"
        self.worker.write_text(FAKE_WORKER, encoding="utf-8")
        self.worker.chmod(0o755)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _one_output(self) -> tuple:
        invocation = subject.make_invocations()[0]
        command = invocation.argv(self.worker)
        import subprocess
        completed = subprocess.run(command, check=True, capture_output=True)
        return invocation, command, completed.stdout, completed.stderr

    def test_frozen_trace_identity_matches_cpp_oracle(self) -> None:
        self.assertEqual(len(self.cell_map), 360)
        self.assertEqual(len(self.trace_map), 360)
        self.assertEqual(
            self.proof["trace_aggregate_sha256"],
            subject.EXPECTED_FROZEN_TRACE_AGGREGATE_SHA256)
        self.assertEqual(
            self.proof["total_attempted_candidates"], 5345254)
        self.assertEqual(
            self.proof["maximum_attempted_candidates"], 128334)
        self.assertEqual(self.proof["root_stratum_case_count"], 12)
        self.assertEqual(self.proof["nested_prefix_count"], 1440)
        for (_trial, _schedule, K), trace in self.trace_map.items():
            self.assertEqual(len(trace.prefix_attempted_candidates), 4)
            for overhead, attempts in zip(
                    subject.OVERHEADS, trace.prefix_attempted_candidates):
                self.assertLessEqual(
                    attempts, (K + overhead) * 256 + 65536)
        self.assertEqual(
            self.proof["status"], "audited_source_equivalence")
        self.assertFalse(self.proof["runtime_packet_ids_observed"])
        self.assertEqual(
            self.proof["authoritative_trace_evidence"],
            "native_frozen_trace_ledger")

    def test_raw_seed_schedule_goldens_and_provenance(self) -> None:
        self.assertEqual(
            [subject._effective_precode_seed(value)
             for value in (0, 1, 2, 255)],
            ["0x487468302aad7105", "0xe6abe1e9a9f7ed1a",
             "0x84e35ba32942692f", "0xe1b6a7f5f5df09f0"])
        self.assertEqual(
            [subject._effective_packet_seed(value)
             for value in (0, 1, 2, 255)],
            ["0x4ec72102", "0xecfe9abb", "0x8b361474", "0xe8096049"])
        provenance = subject._source_provenance(self.worker)
        self.assertEqual(
            provenance["wirehair_v2_raw_seed_source_sha256"],
            hashlib.sha256(subject.RAW_SEED_SOURCE.read_bytes()).hexdigest())
        self.assertEqual(subject.ARM_DESCRIPTOR_SHA256, {
            "wirehair2_raw_d12_h12_periodic":
                "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11",
            "wirehair2_dense_two07_basis_v1":
                "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388",
        })
        for arm in subject.ARM_TRANSFORMS:
            self.assertEqual(subject._arm_descriptor_sha256(arm),
                             subject.ARM_DESCRIPTOR_SHA256[arm])

    def test_one_invocation_is_exactly_120_bound_rows(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        result = subject.parse_invocation_output(
            invocation, command, stdout, stderr, self.cell_map,
            self.trace_map, hashlib.sha256(self.worker.read_bytes()).hexdigest(),
            "1" * 40)
        self.assertEqual(len(result.rows), 120)
        self.assertEqual(len({(row["K"], row["overhead"])
                              for row in result.rows}), 120)
        first = result.rows[0]
        self.assertEqual(
            command[-2:], ["--construction-seed-basis",
                           subject.RAW_SEED_BASIS])
        wide_mode_index = command.index("--cold-solve-wide-xor")
        self.assertEqual(command[wide_mode_index + 1], "policy")
        self.assertEqual(first["heavy_family"], "periodic-cauchy")
        self.assertEqual(first["arm"],
                         "wirehair2_raw_d12_h12_periodic")
        self.assertEqual(first["precode_attempt"], 0)
        self.assertEqual(first["packet_attempt"], 0)
        self.assertEqual(first["construction_seed_basis"],
                         subject.RAW_SEED_BASIS)
        self.assertEqual(first["seed_schedule_sha256"],
                         subject.RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(first["effective_precode_seed"],
                         "0x487468302aad7105")
        self.assertEqual(first["effective_packet_seed"], "0x4ec72102")
        self.assertEqual(first["staircase"], 2)
        self.assertEqual(first["binary_dense_rows"], 12)
        self.assertEqual(first["gf256_heavy_rows"], 12)
        self.assertEqual(first["dense_anchor_layout"], "disabled")
        self.assertEqual(first["source_hits"], 2)
        self.assertIs(first["dense_identity_corner"], False)
        self.assertEqual(first["arm_descriptor_sha256"],
                         subject.ARM_DESCRIPTOR_SHA256[first["arm"]])
        self.assertEqual(first["packet_count"], first["K"])
        self.assertEqual(set(first).isdisjoint({
            "solve_ms_mu", "build_ms_mu", "peel_ms_mu",
        }), True)

    def test_metadata_tamper_is_rejected(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        stdout = stdout.replace(b"trials=1", b"trials=2", 1)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "metadata trials"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

        invocation, command, stdout, stderr = self._one_output()
        stdout = stdout.replace(b"1" * 40, b"2" * 40, 1)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "metadata source_git_commit"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

    def test_raw_seed_schedule_tamper_is_rejected(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        stdout = stdout.replace(
            subject.RAW_SEED_SCHEDULE_SHA256.encode("ascii"),
            b"0" * 64, 1)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "metadata seed_schedule_sha256"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

        invocation, command, stdout, stderr = self._one_output()
        lines = stdout.decode("ascii").splitlines()
        fields = lines[2].split(",")
        source_hits_index = subject.CSV_HEADER.index("source_hits")
        fields[source_hits_index] = "3"
        lines[2] = ",".join(fields)
        stdout = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "realized raw structure"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

        invocation, command, stdout, stderr = self._one_output()
        lines = stdout.decode("ascii").splitlines()
        fields = lines[2].split(",")
        seed_index = subject.CSV_HEADER.index("effective_precode_seed")
        fields[seed_index] = "0x487468302aad7106"
        lines[2] = ",".join(fields)
        stdout = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "field effective_precode_seed"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

    def test_worker_error_row_is_rejected(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        lines = stdout.decode("ascii").splitlines()
        fields = lines[2].split(",")
        success_index = subject.CSV_HEADER.index("success")
        fields[success_index:success_index + 4] = \
            ["0", "0", "1", "1.00000000"]
        lines[2] = ",".join(fields)
        stdout = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "error or invalid terminal"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64, "1" * 40)

    def test_running_stdout_cap_kills_worker_immediately(self) -> None:
        invocation = subject.make_invocations()[0]
        registry = subject._ProcessRegistry()
        started = time.monotonic()
        with mock.patch.dict(
                os.environ,
                {"WH2_FAKE_PRECODEFAIL_CORRUPTION": "stdout_overflow"}):
            with self.assertRaisesRegex(subject.WorkScreenError,
                                        "stdout exceeded"):
                subject._run_invocation(
                    invocation, self.worker, time.monotonic() + 20.0,
                    registry, self.cell_map, self.trace_map, "0" * 64,
                    "1" * 40)
        self.assertLess(time.monotonic() - started, 5.0)
        self.assertTrue(registry.cancelled.is_set())
        self.assertFalse(registry.processes)

    def test_invocations_execute_the_pinned_worker_inode(self) -> None:
        invocation = subject.make_invocations()[0]
        resolved, worker_fd = subject._open_pinned_worker(self.worker)
        pinned_sha256 = subject._sha256_fd(worker_fd, "test worker")
        replacement = self.root / "replacement-worker.py"
        replacement.write_text(
            "#!/usr/bin/env python3\nraise SystemExit(91)\n",
            encoding="utf-8")
        replacement.chmod(0o755)
        os.replace(str(replacement), str(self.worker))
        self.assertNotEqual(
            pinned_sha256, hashlib.sha256(self.worker.read_bytes()).hexdigest())
        registry = subject._ProcessRegistry()
        try:
            result = subject._run_invocation(
                invocation, resolved, time.monotonic() + 20.0,
                registry, self.cell_map, self.trace_map, pinned_sha256,
                "1" * 40, worker_fd)
        finally:
            os.close(worker_fd)
        self.assertEqual(len(result.rows), 120)
        self.assertTrue(all(
            row["worker_binary_sha256"] == pinned_sha256
            for row in result.rows))

    def test_partial_submission_failure_shuts_down_executor(self) -> None:
        class FailingExecutor:
            instance = None

            def __init__(self, max_workers):
                self.max_workers = max_workers
                self.futures = []
                self.shutdown_called = False
                FailingExecutor.instance = self

            def submit(self, *args, **kwargs):
                del args, kwargs
                if self.futures:
                    raise RuntimeError("injected submit failure")
                future = concurrent.futures.Future()
                self.futures.append(future)
                return future

            def shutdown(self, wait=True):
                self.shutdown_called = wait

        provenance = {
            "source_git_commit": "1" * 40,
            "tracked_worktree_clean": True,
            "tracked_status_sha256": hashlib.sha256(b"").hexdigest(),
            "untracked_source_clean": True,
            "untracked_source_sha256": hashlib.sha256(b"").hexdigest(),
            "worker_source_git_commit": "1" * 40,
            "worker_binary_sha256": hashlib.sha256(
                self.worker.read_bytes()).hexdigest(),
            "controller_source_sha256": "4" * 64,
            "wirehair_v2_bench_source_sha256": "2" * 64,
            "wirehair_v2_raw_seed_source_sha256": "5" * 64,
            "wh2_frozen_trace_source_sha256": "3" * 64,
        }
        with mock.patch.object(subject, "_source_provenance",
                               return_value=provenance), \
                mock.patch.object(
                    subject, "build_frozen_domain",
                    return_value=(self.cell_map, self.trace_map, self.proof)), \
                mock.patch.object(
                    subject.concurrent.futures, "ThreadPoolExecutor",
                    FailingExecutor):
            with self.assertRaisesRegex(RuntimeError,
                                        "injected submit failure"):
                subject.run_campaign(
                    self.worker, self.root / "submit-failure", jobs=2,
                    deadline_seconds=20.0)
        executor = FailingExecutor.instance
        self.assertIsNotNone(executor)
        self.assertTrue(executor.shutdown_called)
        self.assertTrue(executor.futures[0].cancelled())

    def test_dirty_tracked_worktree_is_rejected_before_launch(self) -> None:
        provenance = {
            "source_git_commit": "1" * 40,
            "tracked_worktree_clean": False,
            "tracked_status_sha256": hashlib.sha256(b" M dirty").hexdigest(),
            "untracked_source_clean": True,
            "untracked_source_sha256": hashlib.sha256(b"").hexdigest(),
            "worker_source_git_commit": "1" * 40,
            "worker_binary_sha256": hashlib.sha256(
                self.worker.read_bytes()).hexdigest(),
            "controller_source_sha256": "4" * 64,
            "wirehair_v2_bench_source_sha256": "2" * 64,
            "wirehair_v2_raw_seed_source_sha256": "5" * 64,
            "wh2_frozen_trace_source_sha256": "3" * 64,
        }
        with mock.patch.object(subject, "_source_provenance",
                               return_value=provenance), \
                mock.patch.object(
                    subject, "build_frozen_domain",
                    return_value=(self.cell_map, self.trace_map, self.proof)):
            with self.assertRaisesRegex(subject.WorkScreenError,
                                        "clean codec source worktree"):
                subject.run_campaign(
                    self.worker, self.root / "dirty", jobs=1,
                    deadline_seconds=20.0)

    def test_git_source_state_detects_untracked_build_inputs(self) -> None:
        repo = self.root / "source-state"
        (repo / "bench").mkdir(parents=True)
        (repo / "codec").mkdir()
        (repo / "CMakeLists.txt").write_text(
            "cmake_minimum_required(VERSION 3.10)\n")
        (repo / "codec" / "fixture.cpp").write_text("int fixture = 1;\n")
        subprocess.run(["git", "init", "-q"], cwd=str(repo), check=True)
        subprocess.run(
            ["git", "config", "user.email", "test@example.com"],
            cwd=str(repo), check=True)
        subprocess.run(
            ["git", "config", "user.name", "Test"], cwd=str(repo),
            check=True)
        subprocess.run(
            ["git", "add", "CMakeLists.txt", "codec/fixture.cpp"],
            cwd=str(repo), check=True)
        subprocess.run(
            ["git", "commit", "-qm", "fixture"], cwd=str(repo),
            check=True)
        clean = subject._git_source_state(repo)
        self.assertTrue(clean["tracked_worktree_clean"])
        self.assertTrue(clean["untracked_source_clean"])
        (repo / "preferred-route.csv").write_text("user artifact\n")
        self.assertTrue(
            subject._git_source_state(repo)["untracked_source_clean"])
        (repo / "bench" / "untracked.py").write_text("pass\n")
        dirty = subject._git_source_state(repo)
        self.assertFalse(dirty["untracked_source_clean"])
        self.assertNotEqual(
            dirty["untracked_source_sha256"], hashlib.sha256(b"").hexdigest())

    def test_artifact_child_reads_are_pinned_to_open_directory(self) -> None:
        source = self.root / "pinned-artifact"
        source.mkdir()
        (source / "record").write_bytes(b"original")
        directory_fd = os.open(
            str(source), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            source.rename(self.root / "held-artifact")
            source.mkdir()
            (source / "record").write_bytes(b"substitute")
            self.assertEqual(
                subject._read_regular_bytes(
                    Path("record"), 1024, "pinned record", directory_fd),
                b"original")
        finally:
            os.close(directory_fd)

    def test_nonfinite_summary_json_fails_with_domain_error(self) -> None:
        path = self.root / "nonfinite-summary.json"
        path.write_bytes(b'{"value":NaN}\n')
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "nonstandard JSON constant"):
            subject._read_canonical_json(path, 1024)

    def test_partial_artifact_pair_publication_is_rolled_back(self) -> None:
        result_path = self.root / "pair-results.jsonl"
        summary_path = self.root / "pair-summary.json"
        real_atomic_write_at = subject._atomic_write_at
        calls = []

        def fail_second(directory_fd, name, display_path, data, expected):
            calls.append(display_path)
            if len(calls) == 2:
                raise subject.WorkScreenError("injected summary failure")
            return real_atomic_write_at(
                directory_fd, name, display_path, data, expected)

        with mock.patch.object(
                subject, "_atomic_write_at", side_effect=fail_second):
            with self.assertRaisesRegex(subject.WorkScreenError,
                                        "injected summary failure"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertEqual(calls, [result_path, summary_path])
        self.assertFalse(result_path.exists())
        self.assertFalse(summary_path.exists())

    def test_renameat2_capability_fails_before_publication(self) -> None:
        result_path = self.root / "unsupported-results.jsonl"
        summary_path = self.root / "unsupported-summary.json"
        with mock.patch.object(subject, "_RENAMEAT2", None):
            with self.assertRaisesRegex(
                    subject.WorkScreenError, "requires Linux renameat2"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertFalse(result_path.exists())
        self.assertFalse(summary_path.exists())
        self.assertEqual(list(self.root.iterdir()), [self.worker])

    def test_initial_staging_fstat_interrupt_does_not_leak(self) -> None:
        destination = self.root / "initial-fstat.json"
        directory_fd = os.open(
            str(self.root), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            with mock.patch.object(
                    subject.os, "fstat",
                    side_effect=KeyboardInterrupt("injected initial fstat")):
                with self.assertRaisesRegex(
                        KeyboardInterrupt, "injected initial fstat"):
                    subject._atomic_write_at(
                        directory_fd, destination.name, destination,
                        b"artifact\n", {})
        finally:
            os.close(directory_fd)
        self.assertFalse(destination.exists())
        self.assertEqual(
            [path for path in self.root.iterdir() if path != self.worker], [])

    def test_post_publication_base_exception_rolls_back_pair(self) -> None:
        result_path = self.root / "interrupt-results.jsonl"
        summary_path = self.root / "interrupt-summary.json"
        real_atomic_write_at = subject._atomic_write_at

        def interrupt_after_write(
                directory_fd, name, display_path, data, expected):
            identity = real_atomic_write_at(
                directory_fd, name, display_path, data, expected)
            if display_path == summary_path:
                raise KeyboardInterrupt("injected after summary publication")
            return identity

        with mock.patch.object(
                subject, "_atomic_write_at",
                side_effect=interrupt_after_write):
            with self.assertRaisesRegex(
                    KeyboardInterrupt,
                    "injected after summary publication"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertFalse(result_path.exists())
        self.assertFalse(summary_path.exists())

    def test_result_replacement_before_pair_success_is_detected(self) -> None:
        result_path = self.root / "replacement-results.jsonl"
        summary_path = self.root / "replacement-summary.json"
        real_atomic_write_at = subject._atomic_write_at

        def replace_result_after_summary(
                directory_fd, name, display_path, data, expected):
            identity = real_atomic_write_at(
                directory_fd, name, display_path, data, expected)
            if display_path == summary_path:
                os.unlink(result_path.name, dir_fd=directory_fd)
                descriptor = os.open(
                    result_path.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o600, dir_fd=directory_fd)
                try:
                    os.write(descriptor, b"forged\n")
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
            return identity

        with mock.patch.object(
                subject, "_atomic_write_at",
                side_effect=replace_result_after_summary):
            with self.assertRaisesRegex(
                    subject.WorkScreenError,
                    "published artifact identity changed"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertEqual(result_path.read_bytes(), b"forged\n")
        self.assertFalse(summary_path.exists())

    def test_staging_cleanup_interrupt_rolls_back_pair(self) -> None:
        result_path = self.root / "cleanup-results.jsonl"
        summary_path = self.root / "cleanup-summary.json"
        with mock.patch.object(
                subject, "_discard_staged_at",
                side_effect=KeyboardInterrupt("injected staging cleanup")):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "injected staging cleanup"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertFalse(result_path.exists())
        self.assertFalse(summary_path.exists())
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_single_write_cleanup_interrupt_rolls_back(self) -> None:
        destination = self.root / "single-cleanup.json"
        with mock.patch.object(
                subject, "_discard_staged_at",
                side_effect=KeyboardInterrupt("injected staging cleanup")):
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "injected staging cleanup"):
                subject._atomic_write(destination, b"artifact\n")
        self.assertFalse(destination.exists())
        self.assertEqual(list(self.root.glob(".*.tmp")), [])

    def test_parent_swap_during_cleanup_cannot_return_success(self) -> None:
        output_dir = self.root / "cleanup-parent"
        held_dir = self.root / "cleanup-parent-held"
        output_dir.mkdir()
        result_path = output_dir / "result.jsonl"
        summary_path = output_dir / "summary.json"
        real_discard = subject._discard_staged_at
        swapped = []

        def discard_then_swap(directory_fd, expected, display_path):
            real_discard(directory_fd, expected, display_path)
            if not swapped:
                output_dir.rename(held_dir)
                output_dir.mkdir()
                swapped.append(True)

        with mock.patch.object(
                subject, "_discard_staged_at",
                side_effect=discard_then_swap):
            with self.assertRaisesRegex(
                    subject.WorkScreenError, "directory identity changed"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertFalse(result_path.exists())
        self.assertFalse(summary_path.exists())
        self.assertEqual(list(held_dir.iterdir()), [])

    def test_replacement_after_staging_cleanup_survives_rollback(self) -> None:
        result_path = self.root / "aba-results.jsonl"
        summary_path = self.root / "aba-summary.json"
        real_discard = subject._discard_staged_at
        replaced = []

        def discard_then_replace(directory_fd, expected, display_path):
            real_discard(directory_fd, expected, display_path)
            if display_path == result_path and not replaced:
                os.unlink(result_path.name, dir_fd=directory_fd)
                descriptor = os.open(
                    result_path.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o600, dir_fd=directory_fd)
                try:
                    os.write(descriptor, b"foreign replacement\n")
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
                replaced.append(True)

        with mock.patch.object(
                subject, "_discard_staged_at",
                side_effect=discard_then_replace):
            with self.assertRaisesRegex(
                    subject.WorkScreenError,
                    "published artifact identity changed"):
                subject._publish_artifact_pair(
                    result_path, b"result\n", summary_path, b"summary\n")
        self.assertEqual(result_path.read_bytes(), b"foreign replacement\n")
        self.assertFalse(summary_path.exists())

    def test_verifier_opens_raced_fifo_nonblocking(self) -> None:
        destination = self.root / "fifo-race"
        destination.write_bytes(b"x")
        retained_fd = os.open(str(destination), os.O_RDONLY)
        info = os.fstat(retained_fd)
        expected = subject._PublishedArtifact(
            info.st_dev, info.st_ino, retained_fd, "unused-stage",
            1, hashlib.sha256(b"x").hexdigest())
        directory_fd = os.open(
            str(self.root), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        real_open = os.open
        raced = []

        def replace_with_fifo(path, flags, *args, **kwargs):
            if path == destination.name and not raced:
                os.unlink(destination)
                os.mkfifo(destination)
                raced.append(True)
            return real_open(path, flags, *args, **kwargs)

        started = time.monotonic()
        try:
            with mock.patch.object(
                    subject.os, "open", side_effect=replace_with_fifo):
                with self.assertRaisesRegex(
                        subject.WorkScreenError,
                        "identity changed while opening"):
                    subject._verify_identity_at(
                        directory_fd, destination.name, expected,
                        destination, require_staged=False)
        finally:
            os.close(directory_fd)
            os.close(retained_fd)
        self.assertLess(time.monotonic() - started, 0.5)

    def test_atomic_publication_pins_parent_directory(self) -> None:
        output_dir = self.root / "atomic-output"
        held_dir = self.root / "atomic-output-held"
        output_dir.mkdir()
        destination = output_dir / "result.json"
        real_link = os.link
        swapped = []

        def swap_parent_before_link(source, target, **kwargs):
            if not swapped:
                output_dir.rename(held_dir)
                output_dir.mkdir()
                (output_dir / source).write_bytes(b"forged staged bytes")
                swapped.append(True)
            return real_link(source, target, **kwargs)

        with mock.patch.object(
                subject.os, "link", side_effect=swap_parent_before_link):
            with self.assertRaisesRegex(
                    subject.WorkScreenError,
                    "directory identity changed during publication"):
                subject._atomic_write(destination, b"authentic bytes")
        self.assertFalse(destination.exists())
        self.assertEqual(list(held_dir.iterdir()), [])

    def test_full_fake_campaign_publishes_canonical_hashes(self) -> None:
        output_dir = self.root / "output"
        binary_sha256 = hashlib.sha256(self.worker.read_bytes()).hexdigest()
        provenance = {
            "source_git_commit": "1" * 40,
            "tracked_worktree_clean": True,
            "tracked_status_sha256": hashlib.sha256(b"").hexdigest(),
            "untracked_source_clean": True,
            "untracked_source_sha256": hashlib.sha256(b"").hexdigest(),
            "worker_source_git_commit": "1" * 40,
            "worker_binary_sha256": binary_sha256,
            "controller_source_sha256": "4" * 64,
            "wirehair_v2_bench_source_sha256": "2" * 64,
            "wirehair_v2_raw_seed_source_sha256": "5" * 64,
            "wh2_frozen_trace_source_sha256": "3" * 64,
        }
        with mock.patch.object(subject, "_source_provenance",
                               return_value=provenance), \
                mock.patch.object(
                    subject, "build_frozen_domain",
                    return_value=(self.cell_map, self.trace_map, self.proof)):
            summary = subject.run_campaign(
                self.worker, output_dir, jobs=8, deadline_seconds=60.0)
        result_bytes = (output_dir / subject.RESULT_NAME).read_bytes()
        summary_bytes = (output_dir / subject.SUMMARY_NAME).read_bytes()
        self.assertTrue(result_bytes.endswith(b"\n"))
        self.assertEqual(
            hashlib.sha256(result_bytes).hexdigest(),
            summary["result_stream_sha256"])
        self.assertEqual(len(result_bytes.splitlines()), 2880)
        parsed_summary = json.loads(summary_bytes)
        self.assertEqual(parsed_summary, summary)
        unsigned = dict(summary)
        del unsigned["summary_sha256"]
        self.assertEqual(
            summary["summary_sha256"], subject._sha256_json(unsigned))
        self.assertFalse(summary["selection_performed"])
        self.assertFalse(summary["timing_metrics_included"])
        self.assertEqual(summary["construction_seed_basis"],
                         subject.RAW_SEED_BASIS)
        self.assertEqual(summary["seed_schedule_sha256"],
                         subject.RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(summary["source_provenance"], provenance)
        self.assertEqual(
            summary["source_provenance_sha256"],
            subject._sha256_json(provenance))
        self.assertEqual([arm["rows"] for arm in summary["arms"]],
                         [1440, 1440])
        for arm in summary["arms"]:
            self.assertEqual(
                arm["overhead_zero_weak_cell_count"],
                sum(item["failed_strata"] for item in arm["weak_units"]))
            transitions = arm["transitions_vs_control"]
            self.assertEqual(transitions["repair_count"],
                             len(transitions["repairs"]))
            self.assertEqual(transitions["introduction_count"],
                             len(transitions["introductions"]))
        self.assertEqual(
            summary["arms"][1]["transitions_vs_control"]["repair_count"], 0)
        self.assertEqual(
            summary["arms"][1]["transitions_vs_control"]
            ["introduction_count"], 0)
        records = []
        for raw in result_bytes.splitlines():
            value = json.loads(raw)
            self.assertEqual(raw.decode("utf-8"), subject._canonical(value))
            self.assertEqual(
                value["source_provenance_sha256"],
                summary["source_provenance_sha256"])
            self.assertNotIn("solve_ms_mu", value)
            self.assertEqual(value["construction_seed_basis"],
                             subject.RAW_SEED_BASIS)
            self.assertEqual(value["seed_schedule_sha256"],
                             subject.RAW_SEED_SCHEDULE_SHA256)
            self.assertEqual(
                value["effective_precode_seed"],
                subject._effective_precode_seed(value["precode_attempt"]))
            self.assertEqual(
                value["effective_packet_seed"],
                subject._effective_packet_seed(value["packet_attempt"]))
            records.append(value)
        for arm in summary["arms"]:
            transitions = arm["transitions_vs_control"]
            for name, expected_pair in (
                    ("repairs", (0, 1)), ("introductions", (1, 0))):
                for transition in transitions[name]:
                    row = records[transition["record_ordinal"]]
                    control = records[transition["control_record_ordinal"]]
                    self.assertEqual(
                        (row["rank_fail"], control["rank_fail"]),
                        expected_pair)
                    self.assertEqual(row["cell_sha256"],
                                     transition["cell_sha256"])
                    self.assertEqual(row["overhead"], transition["overhead"])

        loaded = subject.load_completed_work_screen(
            self.contract, output_dir)
        self.assertEqual(loaded["summary"], summary)
        self.assertEqual(len(loaded["rows"]), subject.EXPECTED_RECORDS)
        self.assertEqual(loaded["summary_sha256"],
                         summary["summary_sha256"])
        self.assertEqual(loaded["result_stream_sha256"],
                         summary["result_stream_sha256"])
        self.assertEqual(loaded["work_domain_sha256"],
                         summary["work_domain_sha256"])
        self.assertEqual(loaded["source_git_commit"], "1" * 40)
        linked_dir = self.root / "linked-output"
        linked_dir.symlink_to(output_dir, target_is_directory=True)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "real directory"):
            subject.load_completed_work_screen(self.contract, linked_dir)

        def write_forged(name, forged_summary, forged_result):
            directory = self.root / name
            directory.mkdir()
            unsigned_forged = dict(forged_summary)
            unsigned_forged.pop("summary_sha256", None)
            forged_summary = dict(unsigned_forged)
            forged_summary["summary_sha256"] = \
                subject._sha256_json(unsigned_forged)
            (directory / subject.RESULT_NAME).write_bytes(forged_result)
            (directory / subject.SUMMARY_NAME).write_bytes(
                (subject._canonical(forged_summary) + "\n").encode("utf-8"))
            return directory

        numeric_row_aliases = {
            "block-bytes-float": ("block_bytes", 2.0),
            "loss-ppm-float": ("loss_ppm", 100000.0),
            "base-attempt-bool": ("base_seed_attempt", False),
            "base-attempt-float": ("base_seed_attempt", 0.0),
            "packet-count-float": ("packet_count", 2.0),
        }
        for name, (field, alias) in numeric_row_aliases.items():
            aliased_rows = [json.loads(raw)
                            for raw in result_bytes.splitlines()]
            self.assertEqual(aliased_rows[0][field], alias)
            aliased_rows[0][field] = alias
            aliased_result = b"".join(
                (subject._canonical(value) + "\n").encode("utf-8")
                for value in aliased_rows)
            aliased_summary = dict(summary)
            aliased_summary["result_stream_sha256"] = hashlib.sha256(
                aliased_result).hexdigest()
            aliased_dir = write_forged(
                name, aliased_summary, aliased_result)
            with self.subTest(exact_row_alias=name), \
                    self.assertRaises(subject.WorkScreenError):
                subject.load_completed_work_screen(
                    self.contract, aliased_dir)

        summary_aliases = {}
        for name, field in (("record-count-float", "record_count"),
                            ("invocation-count-float", "invocation_count")):
            mutation = copy.deepcopy(summary)
            mutation[field] = float(mutation[field])
            summary_aliases[name] = mutation
        mutation = copy.deepcopy(summary)
        mutation["trace_identity"]["total_attempted_candidates"] = float(
            mutation["trace_identity"]["total_attempted_candidates"])
        summary_aliases["trace-identity-float"] = mutation
        mutation = copy.deepcopy(summary)
        mutation["arms"][0]["rows"] = 1440.0
        summary_aliases["arm-rows-float"] = mutation
        mutation = copy.deepcopy(summary)
        mutation["invocations"][0]["row_count"] = 120.0
        summary_aliases["invocation-rows-float"] = mutation
        for name, mutation in summary_aliases.items():
            alias_dir = write_forged(name, mutation, result_bytes)
            with self.subTest(exact_summary_alias=name), \
                    self.assertRaises(subject.WorkScreenError):
                subject.load_completed_work_screen(
                    self.contract, alias_dir)

        forged_rows = [json.loads(raw) for raw in result_bytes.splitlines()]
        forged_rows[0]["staircase"] += 1
        forged_result = b"".join(
            (subject._canonical(value) + "\n").encode("utf-8")
            for value in forged_rows)
        forged_summary = dict(summary)
        forged_summary["result_stream_sha256"] = \
            hashlib.sha256(forged_result).hexdigest()
        forged_dir = write_forged(
            "forged-realized", forged_summary, forged_result)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "staircase disagrees|realized construction hash"):
            subject.load_completed_work_screen(self.contract, forged_dir)

        malformed_rows = [json.loads(raw) for raw in result_bytes.splitlines()]
        malformed_rows[0]["arm"] = []
        malformed_result = b"".join(
            (subject._canonical(value) + "\n").encode("utf-8")
            for value in malformed_rows)
        malformed_summary = dict(summary)
        malformed_summary["result_stream_sha256"] = \
            hashlib.sha256(malformed_result).hexdigest()
        malformed_dir = write_forged(
            "malformed-arm", malformed_summary, malformed_result)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "unknown raw arm descriptor"):
            subject.load_completed_work_screen(self.contract, malformed_dir)

        malformed_provenance = dict(summary)
        malformed_provenance["source_provenance"] = dict(provenance)
        malformed_provenance["source_provenance"]["source_git_commit"] = []
        malformed_provenance["source_provenance_sha256"] = \
            subject._sha256_json(malformed_provenance["source_provenance"])
        malformed_provenance_dir = write_forged(
            "malformed-provenance", malformed_provenance, result_bytes)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "provenance fields"):
            subject.load_completed_work_screen(
                self.contract, malformed_provenance_dir)

        old_summary = dict(summary)
        old_summary["schema"] = "wirehair.wh2.precodefail-work-summary.v1"
        old_dir = write_forged("old-v1", old_summary, result_bytes)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "summary schema"):
            subject.load_completed_work_screen(self.contract, old_dir)

        provenance_summary = dict(summary)
        provenance_summary["source_provenance_sha256"] = "0" * 64
        provenance_dir = write_forged(
            "forged-provenance", provenance_summary, result_bytes)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "source provenance"):
            subject.load_completed_work_screen(
                self.contract, provenance_dir)

        with mock.patch.object(subject, "_source_provenance",
                               return_value=provenance):
            with self.assertRaisesRegex(subject.WorkScreenError,
                                        "refusing to replace"):
                subject.run_campaign(
                    self.worker, output_dir, jobs=1, deadline_seconds=60.0)


if __name__ == "__main__":
    unittest.main()
