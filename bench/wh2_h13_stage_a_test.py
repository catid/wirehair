#!/usr/bin/env python3
"""Boundary and fail-closed tests for the WH2 H13 Stage-A controller."""

from __future__ import annotations

import csv
from contextlib import redirect_stderr, redirect_stdout
import fcntl
import io
import json
import os
from pathlib import Path
import queue
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
from typing import Any
import unittest
from unittest import mock

import wh2_h13_stage_a as campaign


def wait_for_path(path: Path, *, present: bool, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while path.exists() != present and time.monotonic() < deadline:
        time.sleep(0.01)
    if path.exists() != present:
        raise AssertionError(
            f"path presence did not become {present}: {path}",
        )


def read_pid_when_ready(path: Path, timeout: float = 5.0) -> int:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            text = path.read_text(encoding="ascii").strip()
        except FileNotFoundError:
            text = ""
        if text.isdigit() and int(text) > 1:
            return int(text)
        time.sleep(0.01)
    raise AssertionError(f"PID file did not become complete: {path}")


def task(arm: str, ks: tuple[int, ...] = (2,), *, overhead: int = 0,
         job: int | None = None, construction_index: int = 0) -> campaign.Task:
    return campaign._new_task(
        overhead, (0 if arm == "h12" else 1) if job is None else job,
        0, arm, construction_index, 0, "burst", ks,
    )


def row_for(
    value: campaign.Task, K: int, *, outcome: str = "success",
    overrides: dict[str, str] | None = None,
) -> dict[str, str]:
    row = {field: "0" for field in campaign.BENCH_HEADER}
    row.update({
        "N": str(K), "bb": "64", "heavy_family": "periodic",
        "mix_count": "2", "overhead": str(value.overhead), "trials": "1",
        "success": "1" if outcome == "success" else "0",
        "rank_fail": "1" if outcome == "rank_fail" else "0",
        "error": "1" if outcome == "error" else "0",
        "fail_rate": "0.00000000" if outcome == "success" else "1.00000000",
        "inact_mu": "17.000", "inact_max": "17",
        "binary_def_mu": "12.000", "binary_def_max": "12",
        "heavy_gain_mu": "11.000" if outcome == "rank_fail" else "12.000",
        "heavy_gain_min": "11" if outcome == "rank_fail" else "12",
        "heavy_shortfall": "1" if outcome == "rank_fail" else "0",
        "solve_ms_mu": "0.053",
        "build_ms_mu": "0.001", "peel_ms_mu": "0.009",
        "project_ms_mu": "0.003", "residual_ms_mu": "0.039",
        "backsub_ms_mu": "0.001", "seed_attempt": "0",
        "block_xors_mu": "97.000", "block_muladds_mu": "540.000",
        "first_rank_fail": "0" if outcome == "rank_fail" else "-1",
        "binary_def_hist": "12:1",
        "heavy_gain_hist": "11:1" if outcome == "rank_fail" else "12:1",
        "failure_trials": "" if outcome == "success" else "0",
        "active_packet_peel_seed_xor": "0x0",
        "mixed_joint_source_xors_mu": "4.000",
        "mixed_joint_marginal_xors_mu": "5.000",
        "mixed_joint_marginal_copies_mu": "6.000",
        "mixed_joint_active_deltas_mu": "7.000",
        "mixed_joint_scratch_bytes_mu": "8.000",
        "mixed_dual_source_columns_mu": "9.000",
        "construction_attempt": "0",
        "construction_seed_policy": campaign.CONSTRUCTION_SEED_POLICY,
        "construction_seed": str(value.construction_seed),
        "production_seed_fixups_applied": "0",
        "base_matrix_seed": hex(value.construction_seed),
        "base_peel_seed": hex(campaign.construction_peel_seed(
            value.construction_seed)),
        "matrix_seed": hex(value.construction_seed),
        "peel_seed": hex(campaign.construction_peel_seed(
            value.construction_seed)),
        "actual_staircase_rows": str(K),
        "actual_dense_rows": "12",
        "actual_heavy_rows": "12" if value.arm == "h12" else "13",
        "actual_source_hits": "2" if K < 10000 else "3",
        "actual_dense_two_anchor": "1",
        "actual_dense_two_anchor_phase": "0",
        "systematic_probe_result": "0",
        "precode_count": str(K + (24 if value.arm == "h12" else 25)),
        "packet_trace_seed": hex(campaign.loss_seed(value.seed, K)),
        "packet_trace_sha256": "a" * 64,
    })
    if overrides:
        row.update(overrides)
    return row


def output_for(
    value: campaign.Task, *, outcomes: dict[int, str] | None = None,
    overrides: dict[int, dict[str, str]] | None = None,
) -> str:
    stream = io.StringIO(newline="")
    stream.write(campaign.expected_preamble_line(value) + "\n")
    writer = csv.DictWriter(
        stream, fieldnames=campaign.BENCH_HEADER, lineterminator="\n",
    )
    writer.writeheader()
    for K in value.ks:
        writer.writerow(row_for(
            value, K, outcome=(outcomes or {}).get(K, "success"),
            overrides=(overrides or {}).get(K),
        ))
    return stream.getvalue()


def telemetry_row(
    monotonic: int, busy: str = "99.0", *, cpu_temp: str = "65.0",
    dimm_temp: str = "40.0", dimm_errors: str = "0",
    edac_ce: str = "0", edac_ue: str = "0",
) -> bytes:
    values = [
        f"2026-08-01T00:00:{monotonic % 60:02d}Z", str(monotonic), busy,
        "3000.0", cpu_temp, *(dimm_temp for _ in range(8)), dimm_errors,
        "100.0", "100.0", "100.0", edac_ce, edac_ue,
    ]
    assert len(values) == len(campaign.THERMAL_FIELDS)
    return (",".join(values) + "\n").encode("ascii")


class PlannerTests(unittest.TestCase):
    def test_exact_oh0_cardinality_and_chunking(self) -> None:
        tasks = campaign.build_tasks(0)
        self.assertEqual(len(tasks), 4608)
        self.assertEqual(campaign.PAIRED_CELLS, 575991)
        self.assertEqual(
            sum(len(value.ks) for value in tasks), 2 * 575991,
        )
        self.assertTrue(all(len(value.ks) <= 250 for value in tasks))
        self.assertEqual(tasks[0].ks, tuple(range(2, 252)))
        self.assertEqual(tasks[-1].ks[-1], 64000)
        self.assertEqual(len(tasks[-1].ks), 249)

    def test_arm_argv_has_exactly_one_geometry_delta(self) -> None:
        binary = Path("/frozen/wirehair_v2_bench")
        left = campaign.make_benchmark_argv(binary, task("h12"))
        right = campaign.make_benchmark_argv(binary, task("h13"))
        differences = [
            (index, a, b) for index, (a, b) in enumerate(zip(left, right))
            if a != b
        ]
        self.assertEqual(len(differences), 1)
        self.assertEqual(differences[0][1:], ("10", "11"))
        self.assertIn("--raw-attempt0", left)
        self.assertIn("--paired-overhead-stream", left)
        self.assertIn("--binary-dense-two-anchor", left)
        self.assertNotIn("--seed-block-bytes", left)
        self.assertNotIn("--seed-width", left)
        self.assertEqual(
            left[left.index("--construction-seed") + 1],
            str(campaign.CONSTRUCTION_SEEDS[0]),
        )

    def test_exact_uniform_seed_preamble_line_and_mapping(self) -> None:
        value = task("h12", construction_index=2)
        line = campaign.expected_preamble_line(value)
        self.assertIn(
            " seed=0xd1b54a32d192ed03 seed_role=loss-trace-root "
            "construction_seed_policy=matrix-c-peel-lo32-xor-hi32-v1 "
            f"construction_seed={campaign.CONSTRUCTION_SEEDS[2]} "
            "production_seed_fixups_applied=0 raw_attempt0=1 completion=mixed ",
            line,
        )
        self.assertEqual(
            campaign.construction_peel_seed(campaign.CONSTRUCTION_SEEDS[2]),
            ((campaign.CONSTRUCTION_SEEDS[2] & 0xffffffff) ^
             (campaign.CONSTRUCTION_SEEDS[2] >> 32)),
        )

    def test_construction_roots_are_crossed_without_per_k_changes(self) -> None:
        for construction_index, construction_seed in enumerate(
                campaign.CONSTRUCTION_SEEDS):
            value = campaign._new_task(
                0, 0, 0, "h12", construction_index, 1, "adversarial",
                (2, 64000),
            )
            self.assertEqual(value.construction_seed, construction_seed)
            argv = campaign.make_benchmark_argv(Path("/bench"), value)
            self.assertEqual(
                argv[argv.index("--construction-seed") + 1],
                str(construction_seed),
            )
            rows = campaign.parse_output(output_for(value), value)
            self.assertEqual(
                {row["base_matrix_seed"] for row in rows},
                {hex(construction_seed)},
            )

    def test_nested_planner_runs_both_arms_only_on_union(self) -> None:
        cohort = {(0, "burst", 2), (0, "burst", 64000),
                  (2, "repair-only", 41)}
        tasks = campaign.build_tasks(17, 2, cohort)
        self.assertEqual(len(tasks), 4)
        actual = {
            (value.seed_index, value.schedule, K)
            for value in tasks for K in value.ks
        }
        self.assertEqual(actual, cohort)
        self.assertTrue(all(value.construction_index == 2 for value in tasks))
        self.assertEqual([value.arm for value in tasks], ["h12", "h13"] * 2)
        with self.assertRaises(campaign.CampaignError):
            campaign.build_tasks(18, 2, [
                (0, "burst", 2), (0, "burst", 2),
            ])

    def test_unbounded_job_ids_and_checked_accounting(self) -> None:
        before = task("h12", job=99999)
        after = task("h12", job=100000)
        self.assertNotEqual(before.stem, after.stem)
        self.assertIn("job0099999", before.stem)
        self.assertIn("job0100000", after.stem)
        self.assertEqual(
            campaign.task_from_payload(after.payload()), after,
        )
        for field, bad_value in (("job", True), ("overhead", 1.0),
                                 ("seed_index", "0")):
            malformed = after.payload()
            malformed[field] = bad_value
            with self.subTest(field=field):
                with self.assertRaises(campaign.CampaignError):
                    campaign.task_from_payload(malformed)
        with self.assertRaises(campaign.CampaignError):
            campaign.checked_product(1 << 62, 4)

    def test_exact_loss_seed_formula(self) -> None:
        self.assertEqual(
            hex(campaign.loss_seed(campaign.SEEDS[0], 64)),
            "0x8a7aff2a3a348603",
        )
        self.assertEqual(
            hex(campaign.loss_seed(campaign.SEEDS[0], 2)),
            "0x3bca6207163f7b69",
        )

    def test_nan_timeout_is_rejected_before_dispatch(self) -> None:
        errors = io.StringIO()
        with redirect_stderr(errors):
            result = campaign.main([
                "run", "--result-dir", "/does/not/exist", "--cpus", "0",
                "--timeout", "nan",
            ])
        self.assertEqual(result, 2)
        self.assertIn("finite and positive", errors.getvalue())

    def test_explicit_zero_workers_is_not_treated_as_default(self) -> None:
        self.assertEqual(campaign.resolve_worker_count(None, [1, 2]), 2)
        with self.assertRaises(campaign.CampaignError):
            campaign.resolve_worker_count(0, [1, 2])

    def test_freeze_accepts_trusted_host_tool_hardlink_source_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "taskset-a"
            alias = root / "taskset-b"
            source.write_bytes(b"trusted-host-tool")
            os.link(source, alias)
            with self.assertRaises(campaign.CampaignError):
                campaign._copy_unique(source, root / "strict-copy", True)
            digest = campaign._copy_unique(
                source, root / "host-tool-copy", True,
                allow_source_hardlinks=True,
            )
            self.assertEqual(digest, campaign.sha256_file(root / "host-tool-copy"))
            self.assertEqual((root / "host-tool-copy").stat().st_nlink, 1)

    def test_prepare_freezes_exact_4608_job_plan(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "bench"
            binary.write_bytes(b"fake-test-hook-benchmark")
            binary.chmod(0o755)
            result_dir = root / "campaign"
            output = io.StringIO()
            with redirect_stdout(output):
                campaign.prepare(SimpleNamespace(
                    binary=binary, result_dir=result_dir, telemetry_log=None,
                ))
            contract = campaign.load_contract(
                result_dir.resolve(), require_frozen_controller=False,
            )
            self.assertEqual(contract["r0_paired_cells"], 575991)
            self.assertEqual(contract["r1_paired_cells"], 1151982)
            self.assertEqual(contract["full_paired_cells"], 1727973)
            self.assertEqual(contract["production_seed_fixups_applied"], 0)
            for construction_index in range(3):
                tasks = campaign.load_manifest(
                    campaign.stage_path(result_dir, construction_index, 0),
                    construction_index, 0,
                )
                self.assertEqual(len(tasks), 4608)
            self.assertIn('"full_oh0_jobs":13824', output.getvalue())

    def test_oh0_manifest_rejects_a_resealed_truncated_census(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stage = Path(temporary) / "stage"
            tasks = campaign.build_tasks(0)
            campaign.write_sealed_once(
                stage / "manifest.json", f"{campaign.SCHEMA}.stage_manifest",
                campaign.manifest_payload(0, 0, tasks[:-2]),
            )
            with self.assertRaises(campaign.CampaignError):
                campaign.load_manifest(stage, 0, 0)

    def test_exclusive_campaign_lock_rejects_second_controller(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            descriptor = os.open(root / "controller.lock", os.O_RDWR | os.O_CREAT, 0o600)
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)

            @campaign.exclusive_campaign
            def probe(_args, resolved):
                return resolved

            try:
                with self.assertRaises(campaign.CampaignError):
                    probe(SimpleNamespace(result_dir=root))
            finally:
                os.close(descriptor)
            self.assertEqual(probe(SimpleNamespace(result_dir=root)), root.resolve())

    def test_controller_affinity_is_bound_and_receipted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(campaign.os, "sched_setaffinity") as setter, \
                    mock.patch.object(
                        campaign.os, "sched_getaffinity", return_value={2, 3}):
                campaign.bind_controller_affinity(root, [2, 3])
            setter.assert_called_once_with(0, {2, 3})
            receipt = campaign.load_sealed(
                root / "controller_affinity.json",
                f"{campaign.SCHEMA}.controller_affinity",
            )
            self.assertEqual(receipt["actual_cpus"], [2, 3])
            self.assertTrue(receipt["reserved_cpu_127_excluded"])


class SignalLifecycleTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[subprocess.Popen[bytes], Path, Path, Path]:
        value = task("h12")
        output = root / "output.csv"
        output.write_text(output_for(value), encoding="ascii")
        mode = root / "mode"
        mode.write_text("sleep", encoding="ascii")
        child_pid = root / "child.pid"
        binary = root / "fake-benchmark"
        binary.write_text(
            f"#!{sys.executable}\n"
            "import os\n"
            "from pathlib import Path\n"
            "import sys\n"
            "import time\n"
            f"pid_path = Path({str(child_pid)!r})\n"
            f"mode_path = Path({str(mode)!r})\n"
            f"output_path = Path({str(output)!r})\n"
            "pid_path.write_text(str(os.getpid()), encoding='ascii')\n"
            "if mode_path.read_text(encoding='ascii') == 'sleep':\n"
            "    time.sleep(60)\n"
            "sys.stdout.write(output_path.read_text(encoding='ascii'))\n",
            encoding="utf-8",
        )
        binary.chmod(0o555)
        taskset = root / "taskset"
        taskset.write_text(
            f"#!{sys.executable}\n"
            "import os\n"
            "import sys\n"
            "os.execv(sys.argv[3], sys.argv[3:])\n",
            encoding="utf-8",
        )
        taskset.chmod(0o555)
        harness = root / "harness.py"
        harness.write_text(
            "from pathlib import Path\n"
            "import queue\n"
            "import sys\n"
            f"sys.path.insert(0, {str(Path(campaign.__file__).parent)!r})\n"
            "import wh2_h13_stage_a as campaign\n"
            f"root = Path({str(root)!r})\n"
            f"binary = Path({str(binary)!r})\n"
            f"taskset = Path({str(taskset)!r})\n"
            "value = campaign._new_task(0, 0, 0, 'h12', 0, 0, 'burst', (2,))\n"
            "cpus = queue.Queue()\n"
            "cpus.put(2)\n"
            "try:\n"
            "    campaign.dispatch_tasks(\n"
            "        root, root / 'stage', [value], binary,\n"
            "        campaign.sha256_file(binary), taskset,\n"
            "        campaign.sha256_file(taskset), cpus, 1, 30.0, {2})\n"
            "except campaign.CampaignInterrupted as exc:\n"
            "    raise SystemExit(128 + exc.signum)\n",
            encoding="utf-8",
        )
        process = subprocess.Popen(
            [sys.executable, str(harness)], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        read_pid_when_ready(child_pid)
        return process, child_pid, mode, harness

    def _assert_signal_cleanup(self, signals: tuple[int, ...]) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            process, child_pid_path, mode, harness = self._fixture(root)
            child_pid = int(child_pid_path.read_text(encoding="ascii"))
            for signum in signals:
                try:
                    process.send_signal(signum)
                except ProcessLookupError:
                    break
            stdout, stderr = process.communicate(timeout=10)
            self.assertEqual(stdout, b"")
            self.assertEqual(stderr, b"")
            self.assertEqual(process.returncode, 128 + signals[0])
            wait_for_path(Path(f"/proc/{child_pid}"), present=False)
            self.assertFalse((root / "stage" / "shards").exists())

            # The interrupted task left no trusted shard and restarts once.
            mode.write_text("fast", encoding="ascii")
            restarted = subprocess.run(
                [sys.executable, str(harness)], capture_output=True,
                timeout=10, check=False,
            )
            self.assertEqual(restarted.returncode, 0, restarted.stderr)
            final = campaign.shard_path(root / "stage", task("h12"))
            self.assertTrue(final.is_dir())

    def test_sigint_sigterm_and_double_signal_leave_no_child_or_duplicate(self) -> None:
        for signals in (
                (signal.SIGINT,), (signal.SIGTERM,),
                (signal.SIGINT, signal.SIGINT)):
            with self.subTest(signals=signals):
                self._assert_signal_cleanup(signals)

    def test_sigkill_parent_death_signal_leaves_no_detached_child(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            process, child_pid_path, _mode, _harness = self._fixture(
                Path(temporary),
            )
            child_pid = int(child_pid_path.read_text(encoding="ascii"))
            process.kill()
            process.communicate(timeout=10)
            self.assertEqual(process.returncode, -signal.SIGKILL)
            wait_for_path(Path(f"/proc/{child_pid}"), present=False)

    def test_parent_death_launcher_rejects_wrong_expected_parent(self) -> None:
        process = subprocess.run(
            [str(Path(campaign.__file__).resolve()), "__pdeath_exec",
             str(os.getpid() + 1000000), sys.executable, "-c", "pass"],
            capture_output=True, timeout=10, check=False,
        )
        self.assertEqual(process.returncode, -signal.SIGKILL)


class ParserTests(unittest.TestCase):
    def test_strict_valid_receipt_and_need_more(self) -> None:
        value = task("h12", (2, 3))
        parsed = campaign.parse_output(output_for(
            value, outcomes={3: "rank_fail"},
            overrides={3: {"systematic_probe_result": "1"}},
        ), value)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[1]["systematic_probe_result"], "1")

    def test_source_hits_canonical_cutoff_is_sealed(self) -> None:
        value = task("h12", (9999, 10000))
        parsed = campaign.parse_output(output_for(value), value)
        self.assertEqual(
            [row["actual_source_hits"] for row in parsed], ["2", "3"],
        )
        with self.assertRaises(campaign.CampaignError):
            campaign.parse_output(output_for(
                value, overrides={10000: {"actual_source_hits": "2"}},
            ), value)

    def test_rejects_header_preamble_and_numeric_mutations(self) -> None:
        value = task("h12")
        valid = output_for(value)
        lines = valid.splitlines(keepends=True)
        tokens = lines[0].rstrip("\n").split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        cases = [
            " ".join(tokens) + "\n" + "".join(lines[1:]),
            valid.replace("packet_trace_sha256", "packet_trace_sha257", 1),
            valid.replace("N,bb", '"N",bb', 1),
            valid.replace("\n2,64", '\n"2",64', 1),
            output_for(value, overrides={2: {"success": "01"}}),
            output_for(value, overrides={2: {"base_matrix_seed": "0X1234"}}),
            output_for(value, overrides={2: {"packet_trace_seed": "0x0"}}),
            output_for(value, overrides={2: {"error": "1", "success": "0",
                                             "fail_rate": "1.00000000",
                                             "failure_trials": "0"}}),
            output_for(value, overrides={2: {"actual_dense_rows": "13",
                                             "precode_count": "27"}}),
        ]
        for mutated in cases:
            with self.subTest(prefix=mutated[:50]):
                with self.assertRaises(campaign.CampaignError):
                    campaign.parse_output(mutated, value)

    def test_pairing_rejects_receipt_substitution(self) -> None:
        left = task("h12")
        right = task("h13")
        left_rows = campaign.parse_output(output_for(left), left)
        right_rows = campaign.parse_output(output_for(right), right)
        paired = campaign.pair_rows(((left, left_rows), (right, right_rows)))
        self.assertEqual(len(paired), 1)
        with self.assertRaises(campaign.CampaignError):
                campaign.parse_output(output_for(
                    right, overrides={2: {
                        "base_matrix_seed": "0x999", "matrix_seed": "0x999",
                    }},
                ), right)


@unittest.skipUnless(
    os.environ.get("WIREHAIR_V2_BENCH"),
    "set WIREHAIR_V2_BENCH to exercise the native receipt contract",
)
class NativeIntegrationTests(unittest.TestCase):
    def test_uniform_construction_roots_and_k_boundaries(self) -> None:
        binary = Path(os.environ["WIREHAIR_V2_BENCH"]).resolve(strict=True)
        for construction_index in range(len(campaign.CONSTRUCTION_SEEDS)):
            paired: list[tuple[campaign.Task, list[dict[str, str]]]] = []
            for arm in campaign.ARMS:
                value = campaign._new_task(
                    0, construction_index * len(campaign.ARMS) +
                    campaign.ARMS.index(arm), 0, arm, construction_index,
                    1, "adversarial", (2, 9999, 10000, 64000),
                )
                process = subprocess.run(
                    campaign.make_benchmark_argv(binary, value),
                    capture_output=True, text=True, timeout=60, check=False,
                )
                self.assertEqual(process.returncode, 0, process.stderr)
                self.assertEqual(process.stderr, "")
                rows = campaign.parse_output(process.stdout, value)
                self.assertEqual(
                    {row["construction_seed"] for row in rows},
                    {str(campaign.CONSTRUCTION_SEEDS[construction_index])},
                )
                self.assertEqual(
                    {row["production_seed_fixups_applied"] for row in rows},
                    {"0"},
                )
                paired.append((value, rows))
            self.assertEqual(len(campaign.pair_rows(paired)), 4)


class ReceiptAndReducerTests(unittest.TestCase):
    def setUp(self) -> None:
        # Synthetic telemetry rows use a small deterministic boot clock.
        # Keep the seal-time freshness check on that same clock domain.
        patcher = mock.patch.object(
            campaign, "system_uptime_s", return_value=105.0,
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def _verify_tiny_terminal(
        self, root: Path, contract: dict[str, Any],
    ) -> dict[str, Any]:
        with mock.patch.multiple(
                campaign, OH0_JOBS=2, TOTAL_ARM_OUTCOMES=2):
            return campaign.verify_terminal_campaign(root, contract, {2})

    def _write_shard(
        self, stage: Path, value: campaign.Task, binary: Path, taskset: Path,
        *, receipt_task: campaign.Task | None = None,
        outcome: str = "success",
        outcomes: dict[int, str] | None = None,
    ) -> str:
        output = output_for(
            value,
            outcomes=(outcomes if outcomes is not None else
                      {K: outcome for K in value.ks}),
        )
        stdout = output.encode("ascii")
        stderr = b""
        final = campaign.shard_path(stage, value)
        final.mkdir(parents=True)
        benchmark_argv = campaign.make_benchmark_argv(binary, receipt_task or value)
        command = [str(taskset), "-c", "2", *benchmark_argv]
        campaign.atomic_write(final / "stdout.csv", stdout)
        campaign.atomic_write(final / "stderr.txt", stderr)
        campaign.write_sealed_once(
            final / "receipt.json", f"{campaign.SCHEMA}.job_receipt",
            campaign.receipt_payload(
                receipt_task or value, benchmark_argv, command, 2,
                campaign.sha256_file(binary), stdout, stderr, len(value.ks),
            ),
        )
        return output

    def test_restart_accepts_only_exact_sealed_shard(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            value = task("h12")
            self._write_shard(stage, value, binary, taskset)
            digest = campaign.sha256_file(binary)
            first = campaign.validate_shard(
                stage, value, binary, digest, taskset, {2},
            )
            second = campaign.validate_shard(
                stage, value, binary, digest, taskset, {2},
            )
            self.assertEqual(first, second)
            stdout = campaign.shard_path(stage, value) / "stdout.csv"
            stdout.write_bytes(stdout.read_bytes() + b"substitution\n")
            with self.assertRaises(campaign.CampaignError):
                campaign.validate_shard(
                    stage, value, binary, digest, taskset, {2},
                )

    def test_authenticated_fds_survive_visible_path_substitution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            value = task("h12")
            expected = root / "expected.csv"
            expected.write_text(output_for(value), encoding="ascii")
            binary = root / "wirehair_v2_bench"
            binary.write_text(
                f"#!{sys.executable}\n"
                "from pathlib import Path\n"
                "import sys\n"
                f"sys.stdout.write(Path({str(expected)!r}).read_text(encoding='ascii'))\n",
                encoding="utf-8",
            )
            taskset_source = Path(shutil.which("taskset") or "")
            self.assertTrue(taskset_source.is_file())
            taskset = root / "taskset"
            shutil.copyfile(taskset_source, taskset)
            binary.chmod(0o555)
            taskset.chmod(0o555)
            binary_sha = campaign.sha256_file(binary)
            taskset_sha = campaign.sha256_file(taskset)
            binary_fd = campaign.open_authenticated_executable(
                binary, binary_sha,
            )
            taskset_fd = campaign.open_authenticated_executable(
                taskset, taskset_sha,
            )
            try:
                binary.rename(root / "authenticated-binary")
                taskset.rename(root / "authenticated-taskset")
                shutil.copyfile("/bin/false", binary)
                shutil.copyfile("/bin/false", taskset)
                binary.chmod(0o555)
                taskset.chmod(0o555)
                cpu = next(
                    value for value in sorted(os.sched_getaffinity(0))
                    if value != 127
                )
                cpus: queue.Queue[int] = queue.Queue()
                cpus.put(cpu)
                campaign.execute_task(
                    root, stage, value, binary, binary_sha, taskset,
                    taskset_sha, cpus, 10.0, campaign.ActiveChildren(), {cpu},
                    binary_fd, taskset_fd,
                )
            finally:
                os.close(taskset_fd)
                os.close(binary_fd)
            rows = campaign.validate_shard(
                stage, value, binary, binary_sha, taskset, {cpu},
            )
            self.assertEqual(len(rows), 1)
            receipt = campaign.load_sealed(
                campaign.shard_path(stage, value) / "receipt.json",
                f"{campaign.SCHEMA}.job_receipt",
            )
            self.assertEqual(receipt["execution_mode"], campaign.EXECUTION_MODE)

    def test_cross_arm_receipt_substitution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            left = task("h12")
            right = task("h13")
            self._write_shard(stage, right, binary, taskset, receipt_task=left)
            with self.assertRaises(campaign.CampaignError):
                campaign.validate_shard(
                    stage, right, binary, campaign.sha256_file(binary), taskset,
                    {2},
                )

    def test_job_receipt_rejects_bool_aliases_and_unsealed_cpu(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "stage"
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            value = task("h12")
            self._write_shard(stage, value, binary, taskset)
            receipt_path = campaign.shard_path(stage, value) / "receipt.json"
            original = campaign.load_sealed(
                receipt_path, f"{campaign.SCHEMA}.job_receipt",
            )
            mutations = (
                ("cpu", True), ("cpu", 3), ("returncode", False),
                ("row_count", True),
            )
            for field, replacement in mutations:
                with self.subTest(field=field, replacement=replacement):
                    mutated = json.loads(campaign.canonical_json(original))
                    mutated[field] = replacement
                    receipt_path.unlink()
                    campaign.write_sealed_once(
                        receipt_path, f"{campaign.SCHEMA}.job_receipt",
                        mutated,
                    )
                    with self.assertRaises(campaign.CampaignError):
                        campaign.validate_shard(
                            stage, value, binary,
                            campaign.sha256_file(binary), taskset, {2},
                        )
                    receipt_path.unlink()
                    campaign.write_sealed_once(
                        receipt_path, f"{campaign.SCHEMA}.job_receipt",
                        original,
                    )

    def test_metrics_include_mcnemar_and_joint_common_success_work(self) -> None:
        left = task("h12", (2, 3, 4))
        right = task("h13", (2, 3, 4))
        left_rows = campaign.parse_output(output_for(
            left, outcomes={2: "rank_fail"},
        ), left)
        right_rows = campaign.parse_output(output_for(
            right, outcomes={3: "rank_fail"},
        ), right)
        cells = campaign.pair_rows(((left, left_rows), (right, right_rows)))
        report = campaign.aggregate_cells(cells)
        comparison = report["comparison"]
        self.assertEqual(comparison["repairs"], 1)
        self.assertEqual(comparison["introductions"], 1)
        self.assertEqual(comparison["exact_two_sided_mcnemar"], "1")
        self.assertEqual(comparison["common_success"], 1)
        self.assertEqual(
            comparison["common_success_work"]["h13"]
            ["mixed_joint_source_xors_mu_milli_sum"], 4000,
        )
        self.assertEqual(
            comparison["h13_over_h12_common_success_work_ratios"]
            ["mixed_joint_source_xors_mu"], "1",
        )

    def test_weak_k_construction_and_full_k_multiplicity_domains(self) -> None:
        cells = {}
        for construction_index in range(3):
            for seed_index in range(3):
                for schedule in campaign.SCHEDULES:
                    left_task = campaign._new_task(
                        0, 0, 0, "h12", construction_index, seed_index,
                        schedule, (2,),
                    )
                    right_task = campaign._new_task(
                        0, 1, 0, "h13", construction_index, seed_index,
                        schedule, (2,),
                    )
                    cells[(construction_index, seed_index, schedule, 2)] = {
                        "h12": row_for(left_task, 2, outcome="rank_fail"),
                        "h13": row_for(
                            right_task, 2,
                            outcome=("rank_fail"
                                     if construction_index == 0 else "success"),
                        ),
                    }
        report = campaign.aggregate_cells(cells)
        h12 = report["arms"]["h12"]
        h13 = report["arms"]["h13"]
        self.assertEqual(
            h12["failure_strata_histogram_per_K_construction_0_to_9"]["9"], 3,
        )
        self.assertEqual(
            h12["failure_strata_histogram_per_K_0_to_27"]["27"], 1,
        )
        self.assertEqual(
            h13["failure_strata_histogram_per_K_construction_0_to_9"]["9"], 1,
        )
        self.assertEqual(
            h13["failure_strata_histogram_per_K_0_to_27"]["9"], 1,
        )
        scopes = report["descriptive_scopes"]
        self.assertEqual(len(scopes["by_construction_seed"]), 3)
        self.assertEqual(len(scopes["by_packet_trace_root"]), 3)
        self.assertEqual(len(scopes["by_schedule"]), 3)
        self.assertEqual(set(scopes["by_K_band"]), {"2-63"})
        self.assertEqual(len(scopes["by_construction_loss_schedule"]), 27)
        self.assertEqual(
            report["comparison"]["K_cluster_descriptive"]["clusters"], 1,
        )

    def test_round_gates_are_strict_and_heldout_confirmation_can_reject(self) -> None:
        def metrics(h12: int, h13: int) -> dict[str, Any]:
            return {"arms": {
                "h12": {"rank_fail": h12, "error": 0},
                "h13": {"rank_fail": h13, "error": 0},
            }}
        equal = campaign.failure_gate("R0", [0], {0: metrics(7, 7)})
        self.assertFalse(equal["h13_strictly_fewer_failures"])
        root_records = [
            {"construction_seed_index": index} for index in range(3)
        ]
        terminal = campaign.terminal_payload(
            "confirmation_complete", root_records,
            {0: metrics(7, 6), 1: metrics(1, 2), 2: metrics(1, 2)},
        )
        self.assertTrue(terminal["gates"]["R0"][
            "h13_strictly_fewer_failures"])
        self.assertFalse(terminal["gates"]["R1"][
            "h13_strictly_fewer_failures"])
        self.assertFalse(terminal["architecture_accepted"])

    def test_failure_record_distinguishes_error_from_rank_failure(self) -> None:
        left = row_for(task("h12"), 2, outcome="rank_fail")
        right = row_for(task("h13"), 2, outcome="success")
        left["rank_fail"] = "0"
        left["error"] = "1"
        report = campaign.aggregate_cells({
            (0, 0, "burst", 2): {"h12": left, "h13": right},
        })
        failure = report["arms"]["h12"][
            "weak_K_construction_records"
        ][0]["failures"][0]
        self.assertEqual(failure["cause"], "error")

    def test_sealed_record_never_replaces_different_content(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "record.json"
            campaign.write_sealed_once(path, "test.schema", {"a": 1})
            campaign.write_sealed_once(path, "test.schema", {"a": 1})
            with self.assertRaises(campaign.CampaignError):
                campaign.write_sealed_once(path, "test.schema", {"a": 2})

    def test_sealed_record_create_race_never_overwrites_winner(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "record.json"
            winner = (
                campaign.canonical_json(campaign.sealed_record(
                    "test.schema", {"winner": True},
                )) + "\n"
            ).encode()
            original_link = campaign.os.link

            def install_winner_then_link(
                source: Path, destination: Path, *, follow_symlinks: bool,
            ) -> None:
                campaign.atomic_write(destination, winner)
                original_link(
                    source, destination, follow_symlinks=follow_symlinks,
                )

            with mock.patch.object(
                    campaign.os, "link", side_effect=install_winner_then_link):
                with self.assertRaises(campaign.CampaignError):
                    campaign.write_sealed_once(
                        path, "test.schema", {"winner": False},
                    )
            self.assertEqual(path.read_bytes(), winner)

    def test_orphan_seal_temporaries_are_safely_restartable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = campaign.stage_path(root, 0, 0)
            stage.mkdir(parents=True)
            root_orphan = root / ".terminal.json.create-dead"
            root_orphan.write_bytes(b"pre-link\n")
            stage_orphan = stage / ".complete.json.create-dead"
            stage_orphan.write_bytes(b"post-link\n")
            stage_complete = stage / "complete.json"
            os.link(stage_orphan, stage_complete)
            self.assertEqual(stage_orphan.stat().st_nlink, 2)

            campaign.cleanup_orphan_seal_temporaries(root)
            self.assertFalse(root_orphan.exists())
            self.assertFalse(stage_orphan.exists())
            self.assertEqual(stage_complete.read_bytes(), b"post-link\n")
            self.assertEqual(stage_complete.stat().st_nlink, 1)

            unsafe = root / ".telemetry_end.json.create-unsafe"
            unsafe.mkdir()
            with self.assertRaises(campaign.CampaignError):
                campaign.cleanup_orphan_seal_temporaries(root)

    def _write_predeclared_tiny_plans(self, root: Path) -> None:
        for construction_index in range(3):
            values = [
                task("h12", construction_index=construction_index),
                task("h13", construction_index=construction_index),
            ]
            campaign.write_sealed_once(
                campaign.stage_path(root, construction_index, 0) /
                    "manifest.json",
                f"{campaign.SCHEMA}.stage_manifest",
                campaign.manifest_payload(construction_index, 0, values),
            )

    def test_r0_rejection_is_terminal_and_heldout_plans_stay_unexecuted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, mock.patch.multiple(
                campaign, OH0_JOBS=2, OH0_JOBS_PER_ROOT=2,
                TOTAL_ARM_OUTCOMES=2, PAIRED_CELLS=1,
                PAIRED_CELLS_PER_ROOT=1, R0_PAIRED_CELLS=1):
            root = Path(temporary)
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            self._write_predeclared_tiny_plans(root)
            stage = campaign.stage_path(root, 0, 0)
            def tiny_build(overhead: int, construction_index: int,
                           cohort=None):
                self.assertEqual(overhead, 0)
                self.assertIsNone(cohort)
                return [
                    task("h12", construction_index=construction_index),
                    task("h13", construction_index=construction_index),
                ]
            with mock.patch.object(
                    campaign, "build_tasks", side_effect=tiny_build):
                tasks = campaign.load_manifest(stage, 0, 0)
                for value in tasks:
                    self._write_shard(stage, value, binary, taskset)
                completed = campaign.complete_stage(
                    stage, 0, tasks, binary, campaign.sha256_file(binary),
                    taskset, {2},
                )
                root_terminal = campaign.root_terminal_payload(0, stage, 0, 0)
                terminal = campaign.terminal_payload(
                    "r0_rejected", [root_terminal], {0: completed["metrics"]},
                )
                campaign.write_sealed_once(
                    root / "terminal.json", f"{campaign.SCHEMA}.terminal",
                    terminal,
                )
                contract = {
                    "binary": str(binary),
                    "binary_sha256": campaign.sha256_file(binary),
                    "taskset": str(taskset),
                }
                verified = campaign.verify_terminal_campaign(root, contract, {2})
            self.assertEqual(verified["disposition"], "r0_rejected")
            self.assertFalse(verified["architecture_accepted"])
            for construction_index in (1, 2):
                heldout = campaign.stage_path(root, construction_index, 0)
                self.assertEqual(
                    {entry.name for entry in heldout.iterdir()},
                    {"manifest.json"},
                )
            shutil.rmtree(campaign.shard_path(stage, tasks[0]))
            with mock.patch.object(
                    campaign, "build_tasks", side_effect=tiny_build):
                with self.assertRaises(campaign.CampaignError):
                    campaign.verify_terminal_campaign(root, contract, {2})

    def test_campaign_completion_cross_seal_and_null_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign.write_sealed_once(
                root / "contract.json", "test.contract", {"frozen": True},
            )
            campaign.write_sealed_once(
                root / "controller_affinity.json", "test.affinity",
                {"cpus": [2]},
            )
            for construction_index in range(3):
                stage = campaign.stage_path(root, construction_index, 0)
                campaign.write_sealed_once(
                    stage / "manifest.json", "test.manifest",
                    {"construction_seed_index": construction_index},
                )
                if construction_index == 0:
                    campaign.write_sealed_once(
                        stage / "complete.json", "test.complete", {"overhead": 0},
                    )
            stage = campaign.stage_path(root, 0, 0)
            campaign.write_sealed_once(
                root / "terminal.json", f"{campaign.SCHEMA}.terminal",
                {"root_terminals": [{"construction_seed_index": 0,
                                      "terminal_overhead": 0}]},
            )
            campaign.atomic_write(root / "telemetry_start.json", b"stray\n")
            with self.assertRaises(campaign.CampaignError):
                campaign.seal_campaign_completion(root, None)
            (root / "telemetry_start.json").unlink()
            completed = campaign.seal_campaign_completion(root, None)
            self.assertFalse(completed["external_telemetry_bound"])
            verified = campaign.verify_campaign_completion(
                root, {"external_telemetry": {"path": None}}, None, None,
            )
            self.assertEqual(verified, completed)
            self.assertEqual(completed["stage_entry_count"], 3)
            (stage / "manifest.json").write_bytes(b"changed\n")
            with self.assertRaises(campaign.CampaignError):
                campaign.verify_campaign_completion(
                    root, {"external_telemetry": {"path": None}}, None, None,
                )

    def test_terminal_rejects_per_arm_success_reversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            contract = {
                "binary": str(binary),
                "binary_sha256": campaign.sha256_file(binary),
                "taskset": str(taskset),
            }

            second_tasks = [
                task("h12", overhead=1), task("h13", overhead=1),
            ]
            second = campaign.stage_path(root, 0, 1)
            campaign.write_sealed_once(
                second / "manifest.json", f"{campaign.SCHEMA}.stage_manifest",
                campaign.manifest_payload(0, 1, second_tasks),
            )
            self._write_shard(
                second, second_tasks[0], binary, taskset,
                outcome="rank_fail",
            )
            self._write_shard(
                second, second_tasks[1], binary, taskset, outcome="success",
            )
            campaign.complete_stage(
                second, 1, second_tasks, binary,
                campaign.sha256_file(binary), taskset, {2},
            )
            previous = {(0, 0, "burst", 2): {"h12": False, "h13": True}}
            with self.assertRaisesRegex(
                    campaign.CampaignError, "reverted h12 success"):
                campaign.validate_nested_arm_monotonicity(
                    second, second_tasks, binary, campaign.sha256_file(binary),
                    taskset, {2}, previous,
                )

    def test_legal_nested_reduction_accounts_repairs_and_censoring(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, mock.patch.multiple(
                campaign, PAIRED_CELLS=2, PAIRED_CELLS_PER_ROOT=2,
                R0_PAIRED_CELLS=2, R1_PAIRED_CELLS=4,
                FULL_PAIRED_CELLS=6, TOTAL_ARM_OUTCOMES=4,
                OH0_JOBS=2, OH0_JOBS_PER_ROOT=2, MAX_OVERHEAD=2):
            root = Path(temporary)
            binary = root / "bench"
            taskset = root / "taskset"
            binary.write_bytes(b"binary-v1")
            taskset.write_bytes(b"taskset-v1")
            digest = campaign.sha256_file(binary)

            def make_stage(
                construction_index: int, overhead: int, ks: tuple[int, ...],
                left_outcomes: dict[int, str],
                right_outcomes: dict[int, str],
            ) -> None:
                tasks = [
                    task("h12", ks, overhead=overhead,
                         construction_index=construction_index),
                    task("h13", ks, overhead=overhead,
                         construction_index=construction_index),
                ]
                stage = campaign.stage_path(root, construction_index, overhead)
                campaign.write_sealed_once(
                    stage / "manifest.json",
                    f"{campaign.SCHEMA}.stage_manifest",
                    campaign.manifest_payload(
                        construction_index, overhead, tasks),
                )
                self._write_shard(
                    stage, tasks[0], binary, taskset,
                    outcomes=left_outcomes,
                )
                self._write_shard(
                    stage, tasks[1], binary, taskset,
                    outcomes=right_outcomes,
                )
                campaign.complete_stage(
                    stage, overhead, tasks, binary, digest, taskset, {2},
                )

            make_stage(
                0, 0, (2, 3),
                {2: "rank_fail", 3: "rank_fail"},
                {2: "success", 3: "rank_fail"},
            )
            make_stage(
                0, 1, (2, 3),
                {2: "success", 3: "rank_fail"},
                {2: "success", 3: "success"},
            )
            make_stage(
                0, 2, (3,), {3: "rank_fail"}, {3: "success"},
            )
            for construction_index in (1, 2):
                make_stage(
                    construction_index, 0, (2, 3),
                    {2: "success", 3: "success"},
                    {2: "success", 3: "success"},
                )
            campaign.atomic_write(root / "campaign_complete.json", b"present\n")
            contract = {
                "binary": str(binary), "binary_sha256": digest,
                "taskset": str(taskset),
                "external_telemetry": {
                    "path": None, "policy": campaign.TELEMETRY_POLICY,
                    "sampler_identity": "test sampler",
                    "continuity": campaign.TELEMETRY_CONTINUITY,
                },
            }
            terminal = {
                "disposition": "confirmation_complete",
                "executed_construction_seed_indices": [0, 1, 2],
                "root_terminals": [
                    {"construction_seed_index": 0, "terminal_overhead": 2},
                    {"construction_seed_index": 1, "terminal_overhead": 0},
                    {"construction_seed_index": 2, "terminal_overhead": 0},
                ],
            }
            def tiny_build(overhead: int, construction_index: int,
                           cohort=None):
                if cohort is None:
                    ks = (2, 3)
                else:
                    ks = tuple(sorted({K for _, _, K in cohort}))
                return [
                    task("h12", ks, overhead=overhead,
                         construction_index=construction_index),
                    task("h13", ks, overhead=overhead,
                         construction_index=construction_index),
                ]
            with mock.patch.multiple(
                    campaign,
                    load_contract=mock.Mock(return_value=contract),
                    build_tasks=mock.Mock(side_effect=tiny_build),
                    apply_sealed_controller_affinity=mock.Mock(return_value=[2]),
                    telemetry_start=mock.Mock(return_value=None),
                    verify_terminal_campaign=mock.Mock(return_value=terminal),
                    telemetry_finish=mock.Mock(return_value=None),
                    verify_campaign_completion=mock.Mock(
                        return_value={"verified": True},
                    )):
                with redirect_stdout(io.StringIO()):
                    campaign.reduce_campaign.__wrapped__(
                        SimpleNamespace(), root,
                    )
            summary = campaign.load_sealed(
                root / "analysis.json", f"{campaign.SCHEMA}.analysis_record",
            )
            r0 = summary["rounds"]["R0"]
            self.assertEqual(r0["oh0"]["comparison"]["repairs"], 1)
            self.assertEqual(r0["oh0"]["comparison"]["introductions"], 0)
            self.assertEqual(r0["oh0"]["comparison"]["both_fail"], 1)
            self.assertEqual(
                r0["minimum_success_overhead"]["h12"]["all_cells"],
                {"cells": 2, "p99": None, "censored": 1,
                 "maximum_observed": 1,
                 "histogram": {"1": 1, "censored": 1},
                 "survivors_after_observed_overhead": {"0": 2, "1": 1},
                 "right_censored_values_are_null": True},
            )
            self.assertEqual(
                r0["minimum_success_overhead"]["h13"]["all_cells"]
                ["histogram"], {"0": 1, "1": 1},
            )
            self.assertEqual(
                summary["by_construction_seed"]["0"]
                ["unresolved_union_failures"], 1,
            )
            self.assertEqual(summary["rounds"]["R1"]["paired_cells"], 4)
            self.assertEqual(summary["rounds"]["FULL"]["paired_cells"], 6)

    def test_overhead_reports_distinguish_all_and_failure_cohorts(self) -> None:
        self.assertEqual(campaign.overhead_report([0, 0, 0])["p99"], 0)
        report = campaign.overhead_report([1, 2, None])
        self.assertIsNone(report["p99"])
        self.assertEqual(report["histogram"], {"1": 1, "2": 1, "censored": 1})

    def test_external_telemetry_interval_is_nonempty_and_semantically_sealed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            log.write_bytes(
                (",".join(campaign.THERMAL_FIELDS) + "\n").encode("ascii") +
                telemetry_row(100)
            )
            contract = {"external_telemetry": {
                "path": str(log),
                "prepare_mark": campaign.external_file_mark(log),
            }}
            start = campaign.telemetry_start(root, contract)
            with log.open("ab") as output:
                output.write(telemetry_row(101) + telemetry_row(102))
            end = campaign.telemetry_finish(root, contract, start)
            self.assertEqual(end["semantic_audit"]["samples"], 2)
            self.assertEqual(end["semantic_audit"]["edac_ce_delta"], 0)
            self.assertEqual(end["semantic_audit"]["cpu_busy_min_pct"], 99.0)
            self.assertEqual(end["semantic_audit"]["audit_uptime_s"], 105.0)
            self.assertEqual(end["semantic_audit"]["tail_age_seconds"], 3.0)

    def test_external_telemetry_rejects_stale_final_sample(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101))
            with mock.patch.object(
                    campaign, "system_uptime_s", return_value=106.000001):
                with self.assertRaises(campaign.CampaignError):
                    campaign.telemetry_finish(root, contract, start)

    def _start_telemetry(
        self, root: Path, baseline: bytes | None = None,
    ) -> tuple[Path, dict[str, Any], dict[str, Any]]:
        log = root / "thermal.csv"
        log.write_bytes(
            (",".join(campaign.THERMAL_FIELDS) + "\n").encode("ascii") +
            (telemetry_row(100) if baseline is None else baseline)
        )
        contract = {"external_telemetry": {
            "path": str(log),
            "prepare_mark": campaign.external_file_mark(log),
        }}
        start = campaign.telemetry_start(root, contract)
        assert start is not None
        return log, contract, start

    def test_external_telemetry_policy_boundaries(self) -> None:
        valid_rows = (
            telemetry_row(105),
            telemetry_row(101, "95.0"),
            telemetry_row(101, cpu_temp="89.999"),
            telemetry_row(101, dimm_temp="89.999"),
            telemetry_row(101, edac_ce="42", edac_ue="7"),
        )
        valid_baselines = (
            None, None, None, None,
            telemetry_row(100, edac_ce="42", edac_ue="7"),
        )
        for row, baseline in zip(valid_rows, valid_baselines):
            with self.subTest(valid=row):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    log, contract, start = self._start_telemetry(root, baseline)
                    with log.open("ab") as output:
                        output.write(row)
                    end = campaign.telemetry_finish(root, contract, start)
                    self.assertIsNotNone(end)
                    if baseline is not None:
                        self.assertEqual(end["semantic_audit"]["edac_ce"], 42)
                        self.assertEqual(end["semantic_audit"]["edac_ue"], 7)

        invalid = (
            (None, telemetry_row(100)),
            (None, telemetry_row(99)),
            (None, telemetry_row(106)),
            (None, telemetry_row(101, "94.999")),
            (None, telemetry_row(101, cpu_temp="90.0")),
            (None, telemetry_row(101, dimm_temp="90.0")),
            (None, telemetry_row(101, dimm_errors="1")),
            (telemetry_row(100, edac_ce="1"),
             telemetry_row(101, edac_ce="2")),
            (telemetry_row(100, edac_ce="1"),
             telemetry_row(101, edac_ce="0")),
            (telemetry_row(100, edac_ue="1"),
             telemetry_row(101, edac_ue="2")),
        )
        for baseline, row in invalid:
            with self.subTest(invalid=row):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    log, contract, start = self._start_telemetry(root, baseline)
                    with log.open("ab") as output:
                        output.write(row)
                    with self.assertRaises(campaign.CampaignError):
                        campaign.telemetry_finish(root, contract, start)

    def test_restart_keeps_one_continuous_telemetry_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, first_start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101) + telemetry_row(102))
            restarted = campaign.telemetry_start(root, contract)
            self.assertEqual(restarted, first_start)
            with log.open("ab") as output:
                output.write(telemetry_row(103))
            end = campaign.telemetry_finish(root, contract, restarted)
            self.assertEqual(end["semantic_audit"]["samples"], 3)

        for downtime_row in (
                telemetry_row(101, "94.999"), telemetry_row(106)):
            with self.subTest(downtime_row=downtime_row):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    log, contract, start = self._start_telemetry(root)
                    with log.open("ab") as output:
                        output.write(downtime_row)
                    self.assertEqual(campaign.telemetry_start(root, contract), start)
                    with self.assertRaises(campaign.CampaignError):
                        campaign.telemetry_finish(root, contract, start)

    def test_external_telemetry_rejects_header_and_partial_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            log.write_bytes(b"not,the,thermal,header\n" + telemetry_row(100))
            contract = {"external_telemetry": {
                "path": str(log),
                "prepare_mark": campaign.external_file_mark(log),
            }}
            with self.assertRaises(campaign.CampaignError):
                campaign.telemetry_start(root, contract)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101)[:-1])
            with self.assertRaises(campaign.CampaignError):
                campaign.telemetry_finish(root, contract, start)

    def test_external_telemetry_rejects_rotation_race_and_changed_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101))
            original_mark = campaign.external_file_mark
            rotated = False

            def rotate_after_mark(path: Path) -> dict[str, Any]:
                nonlocal rotated
                mark = original_mark(path)
                if not rotated:
                    replacement = root / "replacement.csv"
                    replacement.write_bytes(path.read_bytes())
                    os.replace(replacement, path)
                    rotated = True
                return mark

            with mock.patch.object(
                    campaign, "external_file_mark",
                    side_effect=rotate_after_mark):
                with self.assertRaises(campaign.CampaignError):
                    campaign.telemetry_finish(root, contract, start)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101))
            campaign.telemetry_finish(root, contract, start)
            with log.open("r+b") as output:
                output.seek(start["offset"])
                original = output.read(128)
                changed = original.replace(b",99.0,", b",98.0,", 1)
                self.assertEqual(len(changed), len(original))
                output.seek(start["offset"])
                output.write(changed)
            with self.assertRaises(campaign.CampaignError):
                campaign.telemetry_finish(root, contract, start)

    def test_external_telemetry_rejects_bool_alias_in_end_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log, contract, start = self._start_telemetry(root)
            with log.open("ab") as output:
                output.write(telemetry_row(101))
            campaign.telemetry_finish(root, contract, start)
            end_path = root / "telemetry_end.json"
            end = campaign.load_sealed(
                end_path, f"{campaign.SCHEMA}.telemetry_end",
            )
            end["semantic_audit"]["samples"] = True
            end_path.unlink()
            campaign.write_sealed_once(
                end_path, f"{campaign.SCHEMA}.telemetry_end", end,
            )
            with self.assertRaises(campaign.CampaignError):
                campaign.telemetry_finish(root, contract, start)

    def test_external_telemetry_rejects_empty_interval(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            log.write_bytes(
                (",".join(campaign.THERMAL_FIELDS) + "\n").encode("ascii") +
                telemetry_row(100)
            )
            contract = {"external_telemetry": {
                "path": str(log),
                "prepare_mark": campaign.external_file_mark(log),
            }}
            start = campaign.telemetry_start(root, contract)
            with self.assertRaises(campaign.CampaignError):
                campaign.telemetry_finish(root, contract, start)


if __name__ == "__main__":
    unittest.main()
