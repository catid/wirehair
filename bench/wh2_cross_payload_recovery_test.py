#!/usr/bin/env python3
"""Bounded invariants and tamper tests for cross-payload recovery evidence."""

from __future__ import annotations

from copy import deepcopy
from decimal import Decimal
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import wh2_cross_payload_recovery as subject


def static_preamble(*, degree_balanced: bool = False) -> dict[str, str]:
    value = {
        "completion": "mixed",
        "mixed_period": "244",
        "mixed_gf256_rows": "10",
        "mixed_gf16_rows": "2",
        "mixed_geometry": "frozen",
        "mixed_residue_skew": "0",
        "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0",
        "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": "0",
        "mixed_grouped_gf256_hash_seed": "0x0",
        "mixed_grouped_final_h_a_columns": "0",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none",
        "binary_dense_rows_override": "0",
        "binary_dense_two_anchor": "1",
        "binary_dense_two_anchor_phase": "0",
        "gf256_heavy_rows_override": "0",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0",
        "full_payload_solve": "0",
    }
    if degree_balanced:
        value["degree_balanced_staircase"] = "1"
    return value


def config_fixture() -> dict[str, object]:
    args = [
        "--completion", "mixed", "--mix-count", "2",
        "--mixed-period", "244", "--mixed-geometry", "frozen",
        "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
        "--mixed-grouped-gf256-rows", "0",
        "--mixed-residue-buckets", "auto",
        "--binary-dense-two-anchor",
    ]
    arms = []
    for index, name in enumerate(subject.ARMS):
        arm_args = list(args)
        if index:
            arm_args.append("--degree-balanced-staircase")
        arms.append({
            "name": name, "label": name,
            "binary": "frozen/" + name, "binary_sha256": "a" * 64,
            "source_commit": str(index + 1) * 40,
            "source_tree": str(index + 3) * 40,
            "architecture_argv": arm_args,
            "mix_count": 2,
            "build_receipt": "frozen/{}.receipt".format(name),
            "build_receipt_sha256": "b" * 64,
            "static_preamble": static_preamble(degree_balanced=bool(index)),
        })
    return {"arms": arms}


def make_task(
    phase: str = "discovery", arm: str = "control",
    Ks: tuple[int, ...] = (2, 3), widths: tuple[int, ...] = (64, 256),
    overheads: tuple[int, ...] = (0,), trials: int = 1,
    loss: str = "0.50", schedule: str = "iid",
) -> dict[str, object]:
    config = config_fixture()
    task = {
        "schema": subject.SCHEMA + ".job", "phase": phase, "job": 0,
        "arm": arm, "root_index": 0,
        "seed": (subject.DISCOVERY_SEEDS[0] if phase == "discovery" else
                 subject.EVALUATION_SEEDS[0]),
        "schedule": schedule, "loss": loss, "trials": trials,
        "widths": list(widths), "overheads": list(overheads), "K": list(Ks),
        "cohort_sha256": None if phase == "discovery" else "c" * 64,
        "stdout_max_bytes": subject.MAX_STDOUT_BYTES[phase],
        "stderr_max_bytes": subject.MAX_STDERR_BYTES,
        "timeout_seconds": subject.TIMEOUT_SECONDS,
    }
    task["argv"] = subject.task_argv(task, config)
    return task


def row_text(
    K: int, width: int, overhead: int, trials: int,
    failures: tuple[int, ...] = (), errors: tuple[int, ...] = (),
) -> str:
    failure_set = tuple(sorted(failures + errors))
    rank_fail = len(failures)
    error = len(errors)
    success = trials - rank_fail - error
    rate = Decimal(rank_fail + error) / Decimal(trials)
    fields = {
        "N": str(K), "bb": str(width), "heavy_family": "periodic",
        "mix_count": "2", "overhead": str(overhead), "trials": str(trials),
        "success": str(success), "rank_fail": str(rank_fail),
        "error": str(error), "fail_rate": "{:.8f}".format(rate),
        "inact_mu": "5.000", "inact_max": "5",
        "binary_def_mu": "0.000", "binary_def_max": "0",
        "heavy_gain_mu": "0.000", "heavy_gain_min": "0",
        "heavy_shortfall": "0", "solve_ms_mu": "0.100",
        "build_ms_mu": "0.010", "peel_ms_mu": "0.020",
        "project_ms_mu": "0.020", "residual_ms_mu": "0.020",
        "backsub_ms_mu": "0.030", "seed_attempt": "0",
        "block_xors_mu": "100.000", "block_muladds_mu": "10.000",
        "first_rank_fail": str(-1 if not failures else failures[0]),
        "binary_def_hist": "0:{}".format(trials),
        "heavy_gain_hist": "0:{}".format(trials),
        "failure_trials": "|".join(str(value) for value in failure_set),
        "active_packet_peel_seed_xor": "0x0",
        "mixed_joint_source_xors_mu": "0.000",
        "mixed_joint_marginal_xors_mu": "0.000",
        "mixed_joint_marginal_copies_mu": "0.000",
        "mixed_joint_active_deltas_mu": "0.000",
        "mixed_joint_scratch_bytes_mu": "0.000",
        "mixed_dual_source_columns_mu": "0.000",
    }
    return ",".join(fields[key] for key in subject.CSV_FIELDS)


def output_bytes(
    task: dict[str, object], *, failure_key: tuple[int, int, int] | None = None,
    failure_trials: tuple[int, ...] = (), error_key: tuple[int, int, int] | None = None,
) -> bytes:
    config = config_fixture()
    static = config["arms"][subject.ARMS.index(str(task["arm"]))]["static_preamble"]
    preamble = dict(static)
    preamble.update(subject.expected_dynamic_preamble(task))
    lines = [
        "# precodefail: " + " ".join(
            "{}={}".format(key, preamble[key]) for key in sorted(preamble)),
        ",".join(subject.CSV_FIELDS),
    ]
    for width in task["widths"]:
        for K in task["K"]:
            for overhead in task["overheads"]:
                key = (K, width, overhead)
                failures = failure_trials if key == failure_key else ()
                errors = (0,) if key == error_key else ()
                lines.append(row_text(
                    K, width, overhead, int(task["trials"]), failures, errors))
    return ("\n".join(lines) + "\n").encode("ascii")


class SeedAndArchitectureTests(unittest.TestCase):
    def test_seed_namespaces_are_stable_and_disjoint(self) -> None:
        self.assertEqual(
            subject.DISCOVERY_SEEDS[0],
            subject.derive_seed(subject.SCHEMA + ".discovery", 0))
        self.assertEqual(
            subject.EVALUATION_SEEDS[-1],
            subject.derive_seed(
                subject.SCHEMA + ".evaluation",
                subject.EVALUATION_ROOT_COUNT - 1))
        self.assertFalse(set(subject.DISCOVERY_SEEDS) & set(subject.EVALUATION_SEEDS))

    def test_controls_cover_payload_regimes_at_every_width(self) -> None:
        controls = set(subject.CONTROL_KS)
        for width in subject.WIDTHS:
            for target in (100 * 1024, 1024 * 1024):
                K = (target + width - 1) // width
                self.assertTrue({K - 1, K, K + 1}.intersection(controls))
        self.assertTrue({2, 64000}.issubset(controls))

    def test_architecture_rejects_axes_and_seed_fixes(self) -> None:
        good = ["--completion", "mixed", "--mix-count", "2"]
        self.assertEqual(subject.parse_architecture_argv(good, "good"), tuple(good))
        for bad in (
            good + ["--N", "10"],
            good + ["--heavy-family", "periodic"],
            good + ["--packet-peel-seed-table", "normalized-h15-v2"],
            good + ["--seed-block-bytes", "2"],
            good + ["--mixed-null-witnesses"],
            ["--completion", "certified", "--mix-count", "2"],
            ["--completion", "mixed", "--mix-count", "2,3"],
        ):
            with self.subTest(bad=bad), self.assertRaises(subject.CampaignError):
                subject.parse_architecture_argv(bad, "bad")

    def test_input_arm_requires_fresh_receipt_bound_to_binary_and_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "bench"
            binary.write_bytes(b"fixture executable\n")
            binary.chmod(0o755)
            receipt = root / "build.json"
            source_commit = "1" * 40
            source_tree = "2" * 40

            def write_receipt(commit: str) -> None:
                receipt.write_bytes(subject.canonical_json({
                    "binary_sha256": subject.sha256_file(binary),
                    "fresh_build": True,
                    "source": {"commit": commit, "tree_oid": source_tree},
                }))

            write_receipt(source_commit)
            record = {
                "name": "control", "label": "raw control",
                "binary": str(binary), "source_commit": source_commit,
                "source_tree": source_tree,
                "architecture_argv": [
                    "--completion", "mixed", "--mix-count", "2"],
                "build_receipt": str(receipt),
                "build_receipt_sha256": subject.sha256_file(receipt),
            }
            arm = subject.ArmInput.from_record(record, "control")
            self.assertEqual(arm.source_tree, source_tree)
            write_receipt("3" * 40)
            record["build_receipt_sha256"] = subject.sha256_file(receipt)
            with self.assertRaises(subject.CampaignError):
                subject.ArmInput.from_record(record, "control")


class LedgerTests(unittest.TestCase):
    def test_evaluation_timeout_has_predeclared_high_k_margin(self) -> None:
        self.assertEqual(subject.TIMEOUT_SECONDS, 1200)

    def test_bounded_discovery_is_cartesian_and_paired(self) -> None:
        config = config_fixture()
        jobs = subject.build_discovery_jobs(config, k_lo=2, k_hi=8, chunk=3)
        subject.validate_jobs(jobs, "discovery", config, exact=False)
        self.assertEqual(
            len(jobs),
            subject.DISCOVERY_ROOT_COUNT * len(subject.SCHEDULES) * 3 * 2)
        for left, right in zip(jobs[::2], jobs[1::2]):
            self.assertEqual((left["arm"], right["arm"]), subject.ARMS)
            for key in ("root_index", "seed", "schedule", "loss", "K", "widths"):
                self.assertEqual(left[key], right[key])
        observed = {
            (job["root_index"], job["schedule"], K, job["arm"])
            for job in jobs for K in job["K"]
        }
        expected = {
            (root, schedule, K, arm)
            for root in range(subject.DISCOVERY_ROOT_COUNT)
            for schedule in subject.SCHEDULES for K in range(2, 9)
            for arm in subject.ARMS
        }
        self.assertEqual(observed, expected)

    def test_evaluation_ledger_binds_cohort_and_common_prefixes(self) -> None:
        config = config_fixture()
        cohort = subject.CONTROL_KS
        jobs = subject.build_evaluation_jobs(config, cohort, "c" * 64)
        subject.validate_jobs(
            jobs, "evaluation", config, cohort, "c" * 64, exact=False)
        self.assertTrue(jobs)
        self.assertTrue(all(job["cohort_sha256"] == "c" * 64 for job in jobs))
        self.assertTrue(all(
            job["overheads"] == list(subject.EVALUATION_OVERHEADS)
            and job["argv"][-1] == "--paired-overhead-stream"
            for job in jobs))
        tampered = deepcopy(jobs[:2])
        tampered[1]["seed"] = subject.EVALUATION_SEEDS[1]
        with self.assertRaises(subject.CampaignError):
            subject.validate_jobs(
                tampered, "evaluation", config, cohort, "c" * 64, exact=False)


class ParserAndReceiptTests(unittest.TestCase):
    def test_parser_accepts_exact_multiwidth_overhead_grid(self) -> None:
        task = make_task(
            "evaluation", Ks=(2, 3), widths=subject.WIDTHS,
            overheads=subject.EVALUATION_OVERHEADS, trials=4,
            loss="0.35", schedule="adversarial")
        data = output_bytes(task, failure_key=(3, 1280, 0), failure_trials=(1, 3))
        rows, semantic = subject.parse_output(data, task, config_fixture())
        self.assertEqual(len(rows), 2 * 4 * 3)
        failed = [row for row in rows if row["rank_fail"]]
        self.assertEqual(len(failed), 1)
        self.assertEqual(failed[0]["failures"], (1, 3))
        self.assertRegex(semantic, r"^[0-9a-f]{64}$")
        self.assertEqual(subject.normalized_loss("0.35"), "0.34999999999999998")

    def test_parser_rejects_preamble_grid_and_failure_tampering(self) -> None:
        task = make_task("evaluation", trials=4, overheads=(0, 1, 2))
        original = output_bytes(task, failure_key=(2, 64, 0), failure_trials=(1,))
        variants = (
            original.replace(b"seed_block_bytes_override=0",
                             b"seed_block_bytes_override=2", 1),
            original.replace(b"schedule=iid", b"schedule=burst", 1),
            original.replace(b",0x0,0.000", b",0x1,0.000", 1),
            original.replace(b"2,64,periodic", b"2,4096,periodic", 1),
            original.replace(
                b",4,3,1,0,0.25000000,", b",4,2,2,0,0.25000000,", 1),
        )
        for value in variants:
            with self.subTest(value=value[:80]), self.assertRaises(subject.CampaignError):
                subject.parse_output(value, task, config_fixture())

    def test_parser_receipts_codec_error_for_explicit_reducer_rejection(self) -> None:
        task = make_task("discovery", Ks=(2,), widths=(64,))
        rows, _semantic = subject.parse_output(
            output_bytes(task, error_key=(2, 64, 0)), task, config_fixture())
        self.assertEqual(rows[0]["error"], 1)
        self.assertEqual(rows[0]["failures"], (0,))

    def test_receipt_binds_task_stdout_stderr_and_semantics(self) -> None:
        config = config_fixture()
        task = make_task("discovery", Ks=(2,), widths=(64,))
        data = output_bytes(task)
        _rows, semantic = subject.parse_output(data, task, config)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = subject.phase_paths(root, "discovery", 0)
            subject.write_once(paths["stdout"], data)
            subject.write_once(paths["stderr"], b"")
            receipt = subject.receipt_record(
                task, data, b"", semantic, 100, 200)
            subject.write_once(paths["receipt"], subject.canonical_json(receipt))
            _receipt, rows = subject.verify_receipt(root, task, config)
            self.assertEqual(len(rows), 1)
            os.chmod(paths["stdout"], 0o644)
            paths["stdout"].write_bytes(data.replace(b"100.000", b"101.000", 1))
            with self.assertRaises(subject.CampaignError):
                subject.verify_receipt(root, task, config)


class CohortTests(unittest.TestCase):
    def records(self) -> dict[str, dict[tuple[int, int, str, int], dict[str, int]]]:
        result: dict[str, dict[tuple[int, int, str, int], dict[str, int]]] = {
            arm: {} for arm in subject.ARMS}
        for arm in subject.ARMS:
            for K in range(2, 5):
                for root in range(subject.DISCOVERY_ROOT_COUNT):
                    for schedule in subject.SCHEDULES:
                        for width in subject.WIDTHS:
                            result[arm][(K, root, schedule, width)] = {
                                "rank_fail": 0, "error": 0}
        return result

    def test_cohort_is_union_across_widths_and_both_arms(self) -> None:
        records = self.records()
        records["control"][(2, 0, "iid", 256)]["rank_fail"] = 1
        records["candidate"][(3, 2, "repair-only", 4096)]["rank_fail"] = 1
        derived = subject.derive_cohort(records, k_lo=2, k_hi=4)
        self.assertEqual(derived["weak_union"], [2, 3])
        self.assertIn(2, derived["weak_K_by_width"]["256"])
        self.assertIn(3, derived["weak_K_by_width"]["4096"])
        self.assertEqual(derived["cohort"], [2, 3, 4])

    def test_codec_error_invalidates_discovery(self) -> None:
        records = self.records()
        records["candidate"][(4, 1, "burst", 1280)]["error"] = 1
        with self.assertRaises(subject.CampaignError):
            subject.derive_cohort(records, k_lo=2, k_hi=4)


class ReductionTests(unittest.TestCase):
    @staticmethod
    def metric_row(failures: tuple[int, ...]) -> dict[str, object]:
        return {
            "trials": 4, "rank_fail": len(failures), "error": 0,
            "failures": failures, "block_xors_mu": Decimal("100"),
            "block_muladds_mu": Decimal("10"),
        }

    def test_paired_reducer_reports_repairs_introductions_and_overhead_gate(self) -> None:
        expected = {
            (2, 0, "iid", "0.35", 64, overhead)
            for overhead in subject.EVALUATION_OVERHEADS}
        control_failures = {0: (0, 1), 1: (2,), 2: ()}
        candidate_failures = {0: (1,), 1: (), 2: (3,)}
        records = {arm: {} for arm in subject.ARMS}
        for key in expected:
            overhead = key[-1]
            records["control"][key] = self.metric_row(control_failures[overhead])
            records["candidate"][key] = self.metric_row(
                candidate_failures[overhead])
        summary = subject.summarize_evaluation_records(records, expected)
        self.assertEqual(summary["paired_outcomes"], {
            "both_success": 8, "repair": 2,
            "introduction": 1, "both_fail": 1})
        self.assertEqual(summary["net_failure_change"], -1)
        self.assertFalse(summary["raw_recovery_gate"]["passed"])
        self.assertEqual(summary["repairs_by_overhead"], {"0": 1, "1": 1, "2": 0})
        self.assertEqual(
            summary["introductions_by_overhead"], {"0": 0, "1": 0, "2": 1})
        self.assertEqual(summary["repairs_by_root"]["0"], 2)
        self.assertEqual(summary["introductions_by_root"]["0"], 1)
        self.assertEqual(
            summary["axis_totals"]["control"]["overhead=0"]["block_xors_sum"],
            "400")

        records["candidate"][(2, 0, "iid", "0.35", 64, 2)] = self.metric_row(())
        improved = subject.summarize_evaluation_records(records, expected)
        self.assertTrue(improved["raw_recovery_gate"]["passed"])


class CanonicalArtifactTests(unittest.TestCase):
    def test_canonical_loader_rejects_duplicate_and_noncanonical_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "x.json"
            for data in (b'{"a":1,"a":2}\n', b'{ "a": 1 }\n', b'{"a":NaN}\n'):
                path.write_bytes(data)
                with self.subTest(data=data), self.assertRaises(subject.CampaignError):
                    subject.load_canonical_object(path, "fixture")
            path.write_bytes(subject.canonical_json({"a": 1}))
            self.assertEqual(subject.load_canonical_object(path, "fixture"), {"a": 1})


class ExecutionBoundTests(unittest.TestCase):
    def test_observed_failure_stops_replenishing_the_worker_window(self) -> None:
        tasks = [
            {"phase": "discovery", "job": job}
            for job in range(20)
        ]
        called: list[int] = []
        hold_nonfailing_worker = threading.Event()

        def execute(
            _root: Path, task: dict[str, object], _config: dict[str, object],
        ) -> tuple[int, bool]:
            job = int(task["job"])
            called.append(job)
            if job == 0:
                raise subject.CampaignError("fixture failure")
            hold_nonfailing_worker.wait(0.2)
            return job, False

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject, "verify_config", return_value={}), \
                mock.patch.object(subject, "load_phase_jobs", return_value=tasks), \
                mock.patch.object(subject, "_execute_job", side_effect=execute):
            with self.assertRaises(subject.CampaignError):
                subject.run_phase(Path(temporary), "discovery", 2, True)
        self.assertLessEqual(len(called), 2)
        self.assertTrue(set(called).issubset({0, 1}))

    def test_preflight_rejects_torn_unreceipted_output(self) -> None:
        tasks = [{"phase": "discovery", "job": 0}]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = subject.phase_paths(root, "discovery", 0)
            paths["stdout"].parent.mkdir(parents=True)
            paths["stdout"].write_bytes(b"torn")
            with mock.patch.object(subject, "verify_config", return_value={}), \
                    mock.patch.object(
                        subject, "load_phase_jobs", return_value=tasks):
                with self.assertRaises(subject.CampaignError):
                    subject.run_phase(root, "discovery", 1, False)


if __name__ == "__main__":
    unittest.main()
