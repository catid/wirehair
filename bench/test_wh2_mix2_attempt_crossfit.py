#!/usr/bin/env python3
"""Focused pure tests for the consumed-root MIX2 attempt cross-fit."""

from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_attempt_crossfit as subject


SOURCE_COMMIT = "1" * 40

METADATA_ORDER = (
    "trials", "threads", "loss", "seed", "source_hits_override",
    "packet_peel_seed_xor", "binary_dense_rows_override",
    "gf256_heavy_rows_override", "dense_anchor_layout",
    "odd_packet_peel_seed_xor", "packet_row_seed_multiplier",
    "packet_row_seed_avalanche", "seed_block_bytes_override",
    "overhead_stream", "full_payload_solve", "schedule",
    "cold_solve_wide_xor", "exact_attempt_mode",
    "exact_precode_attempt", "exact_packet_attempt",
    "construction_seed_basis", "seed_schedule_sha256",
    "source_git_commit",
)


def metadata_line(invocation: subject.Invocation) -> str:
    values = {
        "trials": "1",
        "threads": "1",
        "loss": "0.5",
        "seed": invocation.root,
        "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0",
        "overhead_stream": "paired",
        "full_payload_solve": "1",
        "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy",
        "exact_attempt_mode": "1",
        "exact_precode_attempt": str(invocation.attempt),
        "exact_packet_attempt": str(invocation.attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": SOURCE_COMMIT,
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(field, values[field]) for field in METADATA_ORDER)


def csv_values(invocation: subject.Invocation, K: int,
               weak: bool = False) -> dict:
    deficiency = 1 if weak else 0
    gain = 0
    effective_precode = (
        0x123456789ABCDEF0 +
        invocation.attempt * subject.screen.PRECODE_ATTEMPT_STRIDE) & \
        subject.screen.MASK64
    effective_packet = (
        0x12345678 +
        invocation.attempt * subject.screen.PACKET_ATTEMPT_STRIDE) & \
        subject.screen.MASK32
    return {
        "N": str(K),
        "bb": "2",
        "heavy_family": "periodic",
        "mix_count": "2",
        "staircase": str(subject.STAIRCASE_BY_K[K]),
        "binary_dense_rows": "12",
        "gf256_heavy_rows": "12",
        "source_hits": "3" if K >= 10000 else "2",
        "dense_identity_corner": "0",
        "overhead": "0",
        "trials": "1",
        "success": "0" if weak else "1",
        "rank_fail": "1" if weak else "0",
        "error": "0",
        "fail_rate": "1.00000000" if weak else "0.00000000",
        "inact_mu": "10.000",
        "inact_max": "10",
        "binary_def_mu": "{}.000".format(deficiency),
        "binary_def_max": str(deficiency),
        "heavy_gain_mu": "{}.000".format(gain),
        "heavy_gain_min": str(gain),
        "heavy_shortfall": "1" if weak else "0",
        "solve_ms_mu": "1.000",
        "build_ms_mu": "1.000",
        "peel_ms_mu": "0.100",
        "project_ms_mu": "0.200",
        "residual_ms_mu": "0.300",
        "backsub_ms_mu": "0.400",
        "seed_attempt": "",
        "block_xors_mu": "100.000",
        "block_muladds_mu": "50.000",
        "first_rank_fail": "0" if weak else "-1",
        "binary_def_hist": "{}:1".format(deficiency),
        "heavy_gain_hist": "{}:1".format(gain),
        "failure_trials": "0" if weak else "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(invocation.attempt),
        "packet_attempt": str(invocation.attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed": "0x{:016x}".format(effective_precode),
        "effective_packet_seed": "0x{:08x}".format(effective_packet),
    }


def success_stdout(invocation: subject.Invocation, weak_K=()) -> bytes:
    weak = set(weak_K)
    lines = [metadata_line(invocation), ",".join(subject.screen.CSV_HEADER)]
    for K in invocation.K:
        values = csv_values(invocation, K, K in weak)
        lines.append(",".join(
            values[field] for field in subject.screen.CSV_HEADER))
    return ("\n".join(lines) + "\n").encode("ascii")


def configuration_stdout(invocation: subject.Invocation) -> bytes:
    return (metadata_line(invocation) + "\n" +
            ",".join(subject.screen.CSV_HEADER) + "\n").encode("ascii")


def process_result(invocation: subject.Invocation, stdout: bytes,
                   returncode: int = 0, stderr: bytes = b"") \
        -> subject.ProcessResult:
    command_sha256 = subject.sha256_json({
        "argv": invocation.argv(Path("bench"))[1:],
        "invocation": invocation.identity(),
    })
    return subject.ProcessResult(
        invocation=invocation,
        command_sha256=command_sha256,
        stdout_sha256=hashlib.sha256(stdout).hexdigest(),
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
    )


def mutate_csv(stdout: bytes, row_index: int, field: str, value: str) -> bytes:
    lines = stdout.decode("ascii").splitlines()
    columns = lines[2 + row_index].split(",")
    columns[subject.screen.CSV_HEADER.index(field)] = value
    lines[2 + row_index] = ",".join(columns)
    return ("\n".join(lines) + "\n").encode("ascii")


class SyntheticMatrix:
    def __init__(self, weak=(), block_xor_overrides=None):
        self.weak = set(weak)
        self.block_xor_overrides = block_xor_overrides or {}

    def __getitem__(self, key):
        K, attempt, cell = key
        is_weak = key in self.weak
        effective_precode = (
            0x123456789ABCDEF0 +
            attempt * subject.screen.PRECODE_ATTEMPT_STRIDE) & \
            subject.screen.MASK64
        effective_packet = (
            0x12345678 +
            attempt * subject.screen.PACKET_ATTEMPT_STRIDE) & \
            subject.screen.MASK32
        return {
            "outcome": "rank_fail" if is_weak else "success",
            "success": not is_weak,
            "rank_fail": is_weak,
            "error": False,
            "binary_deficiency": 1 if is_weak else 0,
            "gf256_rank_gain": 0,
            "inactivated_columns": 10,
            "heavy_shortfall": 1 if is_weak else 0,
            "effective_precode_seed": "0x{:016x}".format(effective_precode),
            "effective_packet_seed": "0x{:08x}".format(effective_packet),
            "block_xors": self.block_xor_overrides.get(key, 100),
            "gf256_muladds": 50,
            "K": K,
            "attempt": attempt,
            "cell_ordinal": cell,
        }


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_domain(self):
        subject._validate_constants()
        body = subject._contract_body()
        self.assertEqual(
            subject.sha256_json(body), subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.contract_description()["contract_sha256"],
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(
            subject.ROOTS,
            tuple(subject.screen.SELECTION_ROOTS) +
            tuple(subject.screen.ROOTS))
        self.assertFalse(
            set(subject.ROOTS) & set(subject.screen.FINAL_VALIDATION_ROOTS))
        self.assertEqual(subject.MATRIX_CELL_COUNT, 276480)
        self.assertEqual(subject.OOF_CELL_COUNT, 1080)
        self.assertEqual(subject.MAX_MATRIX_PROCESS_COUNT, 16640)

    def test_invocation_argv_is_exact_and_batched(self):
        invocation = subject.Invocation(
            17, 2, 1, (2, 1000, 64000), "matrix")
        self.assertEqual(invocation.cell_ordinal, 7)
        self.assertEqual(invocation.argv(Path("/pinned/bench")), [
            "/pinned/bench", "precodefail", "--N", "2,1000,64000",
            "--bb-list", "2", "--overhead", "0", "--trials", "1",
            "--threads", "1", "--loss", "0.5", "--seed",
            subject.ROOTS[2], "--schedule", "adversarial",
            "--heavy-family", "periodic", "--mix-count", "2",
            "--binary-dense-rows", "12", "--gf256-heavy-rows", "12",
            "--dense-anchors", "two07", "--paired-overhead-stream",
            "--full-payload-solve", "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", "17", "--exact-packet-attempt",
            "17", "--construction-seed-basis", "production-profile",
        ])
        with self.assertRaises(subject.CrossfitError):
            subject.Invocation(17, 0, 0, (1000, 2), "replay")
        with self.assertRaises(subject.CrossfitError):
            subject.Invocation(256, 0, 0, (2,), "anchor")

    def test_controller_receipt_uses_stable_clean_file(self):
        completed = mock.Mock(returncode=0, stdout=b"", stderr=b"")
        with mock.patch.object(subject.subprocess, "run",
                               return_value=completed) as run_status, \
                mock.patch.object(
                    subject.screen, "_stable_file_bytes",
                    return_value=b"controller") as stable_read:
            receipt = subject._controller_receipt(SOURCE_COMMIT)
        self.assertEqual(receipt, {
            "source_git_commit": SOURCE_COMMIT,
            "path": "bench/wh2_mix2_attempt_crossfit.py",
            "sha256": hashlib.sha256(b"controller").hexdigest(),
        })
        run_status.assert_called_once()
        stable_read.assert_called_once_with(
            subject.CONTROLLER_PATH, 4 * 1024 * 1024)
        with self.assertRaises(subject.CrossfitError):
            subject._controller_receipt("not-a-commit")


class ParserTests(unittest.TestCase):
    def test_uint_parser_rejects_overlong_input_fail_closed(self):
        with self.assertRaises(subject.CrossfitError):
            subject._parse_uint("9" * 5000, subject.UINT64_MAX, "test")

    def test_success_parser_accepts_exact_rows(self):
        invocation = subject.Invocation(7, 1, 2, (2, 10000), "matrix")
        stdout = success_stdout(invocation, weak_K=(10000,))
        rows = subject.parse_success_result(
            process_result(invocation, stdout), SOURCE_COMMIT)
        self.assertEqual([row["K"] for row in rows], [2, 10000])
        self.assertTrue(rows[0]["success"])
        self.assertEqual(rows[1]["outcome"], "rank_fail")
        self.assertEqual(rows[1]["binary_deficiency"], 1)

    def test_success_parser_rejects_noncanonical_or_inconsistent_fields(self):
        invocation = subject.Invocation(7, 1, 2, (2,), "matrix")
        stdout = success_stdout(invocation)
        mutations = (
            ("fail_rate", "0.0000000"),
            ("staircase", "3"),
            ("solve_ms_mu", "nan"),
            ("binary_def_hist", "0:01"),
            ("block_xors_mu", "1e2"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                changed = mutate_csv(stdout, 0, field, value)
                with self.assertRaises(subject.CrossfitError):
                    subject.parse_success_result(
                        process_result(invocation, changed), SOURCE_COMMIT)

        bad_receipt = copy.copy(process_result(invocation, stdout))
        object.__setattr__(bad_receipt, "stdout_sha256", "f" * 64)
        with self.assertRaises(subject.CrossfitError):
            subject.parse_success_result(bad_receipt, SOURCE_COMMIT)

    def test_success_parser_rejects_reordered_K_rows(self):
        invocation = subject.Invocation(7, 1, 2, (2, 3), "matrix")
        lines = success_stdout(invocation).decode("ascii").splitlines()
        lines[2], lines[3] = lines[3], lines[2]
        stdout = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaises(subject.CrossfitError):
            subject.parse_success_result(
                process_result(invocation, stdout), SOURCE_COMMIT)

    def test_configuration_failure_is_exact_and_anchor_only(self):
        invocation = subject.Invocation(9, 0, 0, (2,), "anchor")
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.screen.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=9 result=2\n").encode("ascii")
        proof = subject.parse_configuration_failure(
            process_result(invocation, stdout, 2, stderr), SOURCE_COMMIT)
        self.assertEqual(proof["K"], 2)
        self.assertEqual(proof["stderr_sha256"],
                         hashlib.sha256(stderr).hexdigest())

        with self.assertRaises(subject.CrossfitError):
            subject.parse_configuration_failure(
                process_result(invocation, stdout, 2, stderr + b"extra\n"),
                SOURCE_COMMIT)
        replay = subject.Invocation(9, 0, 0, (2,), "replay")
        replay_stdout = (metadata_line(replay) + "\n" +
                         ",".join(subject.screen.CSV_HEADER) +
                         "\n").encode("ascii")
        with self.assertRaises(subject.CrossfitError):
            subject.parse_configuration_failure(
                process_result(replay, replay_stdout, 2, stderr),
                SOURCE_COMMIT)


class AdjudicationTests(unittest.TestCase):
    def test_crossfit_priority_beats_ascending_without_heldout_weakness(self):
        matrix = SyntheticMatrix(weak={
            (2, 0, 0 * len(subject.SCHEDULES)),
            (3, 0, 2 * len(subject.SCHEDULES)),
        })
        verdict = subject.adjudicate(matrix)
        self.assertEqual(verdict["disposition"], "PASS")
        self.assertEqual(verdict["priority_weak_cells"], 0)
        self.assertEqual(verdict["ascending_weak_cells"], 2)
        self.assertEqual(verdict["priority_selected"][(0, 2)], 1)
        self.assertEqual(verdict["baseline_selected"][(0, 2)], 0)

    def test_no_strict_improvement_rejects(self):
        verdict = subject.adjudicate(SyntheticMatrix())
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(verdict["priority_weak_cells"], 0)
        self.assertEqual(verdict["ascending_weak_cells"], 0)

    def test_unselectable_fold_K_is_a_published_scientific_reject(self):
        weak = {
            (2, attempt, 2 * len(subject.SCHEDULES))
            for attempt in subject.ATTEMPTS
        }
        verdict = subject.adjudicate(SyntheticMatrix(weak=weak))
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(verdict["priority_unselectable_fold_K"], 2)
        self.assertEqual(verdict["ascending_unselectable_fold_K"], 2)
        self.assertNotIn((0, 2), verdict["priority_selected"])

    def test_replay_detects_matrix_mismatch(self):
        selected = {
            (fold, K): 0
            for fold in range(len(subject.FOLDS)) for K in subject.K_VALUES
        }
        invocation = subject.Invocation(
            0, subject.FOLDS[0][0], 0, subject.K_VALUES, "replay")
        result = process_result(invocation, success_stdout(invocation))
        matrix = SyntheticMatrix(block_xor_overrides={(2, 0, 0): 101})
        with mock.patch.object(subject, "_parallel", return_value=[result]):
            with self.assertRaisesRegex(
                    subject.CrossfitError,
                    "selected replay differs from the frozen matrix"):
                subject.replay_selected(
                    matrix, selected, mock.sentinel.pinned,
                    SOURCE_COMMIT, 1)

    def test_replay_rejects_duplicate_and_missing_coordinates(self):
        selected = {
            (fold, K): 0
            for fold in range(len(subject.FOLDS)) for K in subject.K_VALUES
        }
        invocation = subject.Invocation(
            0, subject.FOLDS[0][0], 0, subject.K_VALUES, "replay")
        result = process_result(invocation, success_stdout(invocation))
        with mock.patch.object(subject, "_parallel",
                               return_value=[result] * 36):
            with self.assertRaisesRegex(
                    subject.CrossfitError, "duplicate coordinate"):
                subject.replay_selected(
                    SyntheticMatrix(), selected, mock.sentinel.pinned,
                    SOURCE_COMMIT, 1)

    def test_record_ordinal_is_dense_and_bounded(self):
        self.assertEqual(subject._matrix_record_ordinal(2, 0, 0), 0)
        self.assertEqual(
            subject._matrix_record_ordinal(64000, 255, 35),
            subject.MATRIX_CELL_COUNT - 1)


class MatrixTransportTests(unittest.TestCase):
    def test_raw_process_is_pinned_bounded_and_clean_environment(self):
        invocation = subject.Invocation(0, 0, 0, (2,), "anchor")
        pinned = mock.Mock(
            path=Path("/frozen/wirehair_v2_bench"), descriptor=17)
        completed = mock.Mock(
            returncode=0, stdout=b"stdout\n", stderr=b"")
        with mock.patch.object(subject.subprocess, "run",
                               return_value=completed) as run_process:
            result = subject._run_raw(invocation, pinned)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout_sha256,
                         hashlib.sha256(b"stdout\n").hexdigest())
        arguments, keywords = run_process.call_args
        self.assertEqual(arguments[0], invocation.argv(pinned.path))
        self.assertEqual(keywords["env"], subject.CHILD_ENVIRONMENT)
        self.assertEqual(keywords["timeout"],
                         subject.PROCESS_TIMEOUT_SECONDS)
        self.assertEqual(keywords["executable"], "/proc/self/fd/17")
        self.assertEqual(keywords["pass_fds"], (17,))
        self.assertIs(keywords["stdin"], subject.subprocess.DEVNULL)
        self.assertTrue(keywords["start_new_session"])

    def test_two_phase_matrix_caches_only_exact_anchor_failure(self):
        roots = subject.ROOTS[:2]
        schedules = subject.SCHEDULES[:2]

        def fake_parallel(invocations, jobs, pinned):
            del jobs, pinned
            results = []
            for invocation in invocations:
                if (invocation.phase == "anchor" and
                        invocation.K == (3,) and invocation.attempt == 0):
                    stderr = (
                        "precodefail configuration failed N=3 bb=2 "
                        "heavy_family=periodic mix_count=2 "
                        "attempt_mode=exact precode_attempt=0 "
                        "packet_attempt=0 result=2\n").encode("ascii")
                    results.append(process_result(
                        invocation, configuration_stdout(invocation),
                        2, stderr))
                else:
                    results.append(process_result(
                        invocation, success_stdout(invocation)))
            return results

        with mock.patch.multiple(
                subject,
                K_VALUES=(2, 3),
                ROOTS=roots,
                SCHEDULES=schedules,
                ATTEMPTS=(0, 1),
                STAIRCASE_BY_K={2: 2, 3: 3},
                MATRIX_CELL_COUNT=16), \
                mock.patch.object(subject, "_parallel",
                                  side_effect=fake_parallel):
            matrix, counts = subject.build_matrix(
                mock.sentinel.pinned, SOURCE_COMMIT, 1)

        self.assertEqual(len(matrix), 16)
        self.assertEqual(counts, {
            "anchor_processes": 4,
            "batch_processes": 6,
            "configuration_failed_K_attempts": 1,
        })
        for cell in range(4):
            record = matrix[(3, 0, cell)]
            self.assertEqual(record["outcome"], "construct_failed")
            self.assertFalse(record["trace_executed"])
            self.assertIsNotNone(record["configuration_proof_sha256"])
        self.assertTrue(matrix[(3, 1, 3)]["success"])


if __name__ == "__main__":
    unittest.main()
