#!/usr/bin/env python3
"""Lightweight tests for the fixed precodefail work/rank sidecar."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
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
    "N,bb,heavy_family,mix_count,overhead,trials,success,rank_fail,error,"
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
    "odd_packet_peel_seed_xor=0x0", "packet_row_seed_multiplier=0x1",
    "packet_row_seed_avalanche=0", "seed_block_bytes_override=0",
    "overhead_stream=paired", "full_payload_solve=1",
    "schedule={}".format(schedule), "exact_attempt_mode=1",
    "exact_precode_attempt={}".format(trial),
    "exact_packet_attempt={}".format(trial),
    "construction_seed_basis={}".format(seed_basis),
    "seed_schedule_sha256={}".format(RAW_SEED_SCHEDULE_SHA256),
]
corruption = os.environ.get("WH2_FAKE_PRECODEFAIL_CORRUPTION", "")
if corruption == "stdout_overflow" and schedule == "iid" and trial == 0 and \
        heavy == 11:
    import time
    sys.stdout.write("x" * (2 * 1024 * 1024 + 1))
    sys.stdout.flush()
    time.sleep(10)
    raise SystemExit(9)
if corruption == "metadata" and schedule == "iid" and trial == 0 and heavy == 11:
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
        row = [
            str(K), "2", "periodic", "3", str(overhead), "1",
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
                heavy == 11 and K == 2 and overhead == 0:
            row[6:10] = ["0", "0", "1", "1.00000000"]
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

    def test_one_invocation_is_exactly_120_bound_rows(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        result = subject.parse_invocation_output(
            invocation, command, stdout, stderr, self.cell_map,
            self.trace_map, hashlib.sha256(self.worker.read_bytes()).hexdigest())
        self.assertEqual(len(result.rows), 120)
        self.assertEqual(len({(row["K"], row["overhead"])
                              for row in result.rows}), 120)
        first = result.rows[0]
        self.assertEqual(
            command[-2:], ["--construction-seed-basis",
                           subject.RAW_SEED_BASIS])
        self.assertEqual(first["heavy_family"], "periodic")
        self.assertEqual(first["arm"],
                         "wirehair2_raw_d12_h11_periodic")
        self.assertEqual(first["precode_attempt"], 0)
        self.assertEqual(first["packet_attempt"], 0)
        self.assertEqual(first["seed_basis"], subject.RAW_SEED_BASIS)
        self.assertEqual(first["seed_schedule_sha256"],
                         subject.RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(first["effective_precode_seed"],
                         "0x487468302aad7105")
        self.assertEqual(first["effective_packet_seed"], "0x4ec72102")
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
                self.trace_map, "0" * 64)

    def test_raw_seed_schedule_tamper_is_rejected(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        stdout = stdout.replace(
            subject.RAW_SEED_SCHEDULE_SHA256.encode("ascii"),
            b"0" * 64, 1)
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "metadata seed_schedule_sha256"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64)

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
                self.trace_map, "0" * 64)

    def test_worker_error_row_is_rejected(self) -> None:
        invocation, command, stdout, stderr = self._one_output()
        lines = stdout.decode("ascii").splitlines()
        fields = lines[2].split(",")
        fields[6:10] = ["0", "0", "1", "1.00000000"]
        lines[2] = ",".join(fields)
        stdout = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaisesRegex(subject.WorkScreenError,
                                    "error or invalid terminal"):
            subject.parse_invocation_output(
                invocation, command, stdout, stderr, self.cell_map,
                self.trace_map, "0" * 64)

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
                    registry, self.cell_map, self.trace_map, "0" * 64)
        self.assertLess(time.monotonic() - started, 5.0)
        self.assertTrue(registry.cancelled.is_set())
        self.assertFalse(registry.processes)

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

    def test_full_fake_campaign_publishes_canonical_hashes(self) -> None:
        output_dir = self.root / "output"
        binary_sha256 = hashlib.sha256(self.worker.read_bytes()).hexdigest()
        provenance = {
            "source_git_commit": "1" * 40,
            "tracked_worktree_clean": True,
            "tracked_status_sha256": hashlib.sha256(b"").hexdigest(),
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
        self.assertEqual(len(result_bytes.splitlines()), 5760)
        parsed_summary = json.loads(summary_bytes)
        self.assertEqual(parsed_summary, summary)
        unsigned = dict(summary)
        del unsigned["summary_sha256"]
        self.assertEqual(
            summary["summary_sha256"], subject._sha256_json(unsigned))
        self.assertFalse(summary["selection_performed"])
        self.assertFalse(summary["timing_metrics_included"])
        self.assertEqual(summary["seed_basis"], subject.RAW_SEED_BASIS)
        self.assertEqual(summary["seed_schedule_sha256"],
                         subject.RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(summary["source_provenance"], provenance)
        self.assertEqual(
            summary["source_provenance_sha256"],
            subject._sha256_json(provenance))
        self.assertEqual([arm["rows"] for arm in summary["arms"]],
                         [1440, 1440, 1440, 1440])
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
            self.assertEqual(value["seed_basis"], subject.RAW_SEED_BASIS)
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

        with mock.patch.object(subject, "_source_provenance",
                               return_value=provenance):
            with self.assertRaisesRegex(subject.WorkScreenError,
                                        "refusing to replace"):
                subject.run_campaign(
                    self.worker, output_dir, jobs=1, deadline_seconds=60.0)


if __name__ == "__main__":
    unittest.main()
