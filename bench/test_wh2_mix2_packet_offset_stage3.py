#!/usr/bin/env python3
"""Pure mocked tests for sealed MIX2 packet-offset Stage 3."""

from __future__ import annotations

import copy
from contextlib import redirect_stderr
import hashlib
import io
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_packet_offset_stage3 as subject


CONTROLLER_COMMIT = "1" * 40
SOURCE_SHA = "2" * 64
EMPTY_SHA = hashlib.sha256(b"").hexdigest()


def parser_baseline(coordinate_ordinal: int) -> dict:
    coordinate = subject.stage1.Coordinate(coordinate_ordinal)
    return {
        "record_ordinal": coordinate.baseline_record_ordinal,
        "K": coordinate.K,
        "attempt": coordinate.attempt,
        "root_index": coordinate.root_index,
        "loss_root": coordinate.root,
        "schedule_index": coordinate.schedule_index,
        "schedule": coordinate.schedule,
        "cell_ordinal": coordinate.cell_ordinal,
        "outcome": "success",
        "success": True,
        "rank_fail": False,
        "error": False,
        "binary_deficiency": 10,
        "gf256_rank_gain": 10,
        "inactivated_columns": 20,
        "heavy_shortfall": 0,
        "effective_precode_seed": "0x123456789abcdef0",
        "effective_packet_seed": "0x12345678",
        "block_xors": 100,
        "gf256_muladds": 50,
    }


def metadata_line(invocation: subject.Invocation) -> str:
    coordinate = invocation.coordinate
    values = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": coordinate.root, "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "1", "schedule": coordinate.schedule,
        "cold_solve_wide_xor": "policy", "exact_attempt_mode": "1",
        "exact_precode_attempt": str(coordinate.attempt),
        "exact_packet_attempt": str(coordinate.packet_attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT,
        "mix_pair": subject.CANDIDATE_PAIR,
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(field, values[field])
        for field in subject.METADATA_ORDER)


def csv_values(invocation: subject.Invocation,
               outcome: str = "success") -> dict:
    if outcome == "success":
        success, rank_fail, error, deficiency, gain = 1, 0, 0, 10, 10
    elif outcome == "rank_fail":
        success, rank_fail, error, deficiency, gain = 0, 1, 0, 12, 11
    else:
        success, rank_fail, error, deficiency, gain = 0, 0, 1, 0, 0
    weak = success == 0
    coordinate = invocation.coordinate
    baseline = parser_baseline(invocation.coordinate_ordinal)
    return {
        "N": str(coordinate.K), "bb": "2", "heavy_family": "periodic",
        "mix_count": "2",
        "staircase": str(subject.stage1.STAIRCASE_BY_K[coordinate.K]),
        "binary_dense_rows": "12", "gf256_heavy_rows": "12",
        "source_hits": "3" if coordinate.K >= 10000 else "2",
        "dense_identity_corner": "0", "overhead": "0", "trials": "1",
        "success": str(success), "rank_fail": str(rank_fail),
        "error": str(error),
        "fail_rate": "1.00000000" if weak else "0.00000000",
        "inact_mu": "20.000", "inact_max": "20",
        "binary_def_mu": "{}.000".format(deficiency),
        "binary_def_max": str(deficiency),
        "heavy_gain_mu": "{}.000".format(gain),
        "heavy_gain_min": str(gain),
        "heavy_shortfall": "1" if rank_fail else "0",
        "solve_ms_mu": "1.000", "build_ms_mu": "1.000",
        "peel_ms_mu": "1.000", "project_ms_mu": "1.000",
        "residual_ms_mu": "1.000", "backsub_ms_mu": "1.000",
        "seed_attempt": "", "block_xors_mu": "100.000",
        "block_muladds_mu": "50.000",
        "first_rank_fail": "0" if rank_fail else "-1",
        "binary_def_hist": "{}:1".format(deficiency),
        "heavy_gain_hist": "{}:1".format(gain),
        "failure_trials": "0" if weak else "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(coordinate.attempt),
        "packet_attempt": str(coordinate.packet_attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed": baseline["effective_precode_seed"],
        "effective_packet_seed": subject.stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], coordinate.attempt,
            subject.PACKET_DELTA),
    }


def stdout_for(invocation: subject.Invocation,
               outcome: str = "success") -> bytes:
    values = csv_values(invocation, outcome)
    return ("\n".join((
        metadata_line(invocation), ",".join(subject.CSV_HEADER),
        ",".join(values[field] for field in subject.CSV_HEADER))) +
        "\n").encode("ascii")


def process_result(invocation: subject.Invocation, stdout: bytes,
                   returncode: int = 0,
                   stderr: bytes = b"") -> subject.ProcessResult:
    return subject.ProcessResult(
        invocation,
        subject.sha256_json({
            "argv": invocation.argv(Path("wirehair_v2_bench"))[1:]}),
        subject.sha256_bytes(stdout), subject.sha256_bytes(stderr),
        returncode, stdout, stderr)


def full_baseline() -> tuple[dict, ...]:
    rows = []
    for ordinal in range(subject.COORDINATE_COUNT):
        coordinate = subject.stage1.Coordinate(ordinal)
        success = ordinal not in subject.BASELINE_WEAK_COORDINATES
        row = {
            "record_ordinal": coordinate.baseline_record_ordinal,
            "K": coordinate.K,
            "attempt": coordinate.attempt,
            "root_index": coordinate.root_index,
            "loss_root": coordinate.root,
            "schedule_index": coordinate.schedule_index,
            "schedule": coordinate.schedule,
            "cell_ordinal": coordinate.cell_ordinal,
            "outcome": "success" if success else "rank_fail",
            "success": success,
            "rank_fail": not success,
            "error": False,
            "binary_deficiency": 0 if success else 1,
            "gf256_rank_gain": 0,
            "inactivated_columns": 1,
            "heavy_shortfall": int(not success),
            "effective_precode_seed": "0x123456789abcdef0",
            "effective_packet_seed": "0x12345678",
            "block_xors": 1,
            "gf256_muladds": 1,
        }
        rows.append(row)
    # Put the exact stratum totals on successful cells.  This keeps all frozen
    # baseline aggregate checks meaningful without importing a 197 MB fixture.
    for coordinate, count, totals in (
            (0, len(subject.UNTOUCHED_COORDINATES),
             subject.EXPECTED_UNTOUCHED_BASELINE_SUMS),
            (1, len(subject.CONSULTED_COORDINATES),
             subject.EXPECTED_CONSULTED_BASELINE_SUMS)):
        for field, value in totals.items():
            rows[coordinate][field] = value - (count - 1)
    return tuple(rows)


def synthetic_record(invocation: subject.Invocation,
                     baseline: dict) -> dict:
    coordinate = invocation.coordinate
    record = {
        "schema": subject.RECORD_SCHEMA,
        **invocation.identity(),
        "block_bytes": 2, "loss_ppm": 500000, "overhead": 0,
        "dense_anchor_layout": "two07", "mix_count": 2,
        "binary_dense_rows": 12, "gf256_heavy_rows": 12,
        "heavy_family": "periodic",
        "construction_seed_basis": "production-profile",
        "full_payload_solve": True, "overhead_stream": "paired",
        "cold_solve_wide_xor": "policy",
        "promotion_evidence": False, "fresh_roots_used": False,
        "timing_evidence": False,
        "command_sha256": subject.sha256_json({
            "argv": invocation.argv(Path("wirehair_v2_bench"))[1:]}),
        "stdout_sha256": "3" * 64, "stderr_sha256": EMPTY_SHA,
        "returncode": 0,
        "bench_binary_sha256": subject.EXPECTED_STAGE0_BENCH_SHA256,
        "bench_source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT,
        "controller_source_git_commit": CONTROLLER_COMMIT,
        "source_receipt_sha256": SOURCE_SHA,
        "pair01_baseline_record_sha256": subject.sha256_json(baseline),
        "pair01_baseline": subject._baseline_projection(baseline),
        "baseline_matrix_sha256": subject.stage1.BASELINE_FILES[
            "attempt-crossfit-matrix.jsonl"][1],
        "staircase": subject.stage1.STAIRCASE_BY_K[coordinate.K],
        "source_hits": 3 if coordinate.K >= 10000 else 2,
        "outcome": "success", "success": True, "rank_fail": False,
        "error": False, "configuration_failure": False, "weak": False,
        "inactivated_columns": baseline["inactivated_columns"],
        "binary_deficiency": 0, "gf256_rank_gain": 0,
        "heavy_shortfall": 0, "first_rank_fail": -1,
        "failure_trials": [],
        "effective_precode_seed": baseline["effective_precode_seed"],
        "effective_packet_seed": subject.stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], coordinate.attempt,
            subject.PACKET_DELTA),
        "block_xors": baseline["block_xors"],
        "gf256_muladds": baseline["gf256_muladds"],
        "configuration_stderr_sha256": None,
    }
    record["deterministic_projection_sha256"] = subject.sha256_json(
        subject.deterministic_projection(record))
    if set(record) != subject.RECORD_FIELDS:
        raise AssertionError(set(record) ^ subject.RECORD_FIELDS)
    return record


def rehash(record: dict) -> None:
    record["deterministic_projection_sha256"] = subject.sha256_json(
        subject.deterministic_projection(record))


def passing_roster() -> tuple[list[dict], SimpleNamespace]:
    baseline = full_baseline()
    records = [synthetic_record(invocation,
                                baseline[invocation.coordinate_ordinal])
               for invocation in subject.make_invocations()]
    confirmation = {
        coordinate: subject.scientific_projection(records[coordinate * 2])
        for coordinate in subject.CONSULTED_COORDINATES
    }
    evidence = SimpleNamespace(
        baseline=baseline,
        stage2_confirmation_by_coordinate=confirmation)
    return records, evidence


def make_rank_failure(record: dict) -> None:
    record.update({
        "outcome": "rank_fail", "success": False, "rank_fail": True,
        "weak": True, "binary_deficiency": 1, "gf256_rank_gain": 0,
        "heavy_shortfall": 1, "first_rank_fail": 0,
        "failure_trials": [0],
    })
    if record["inactivated_columns"] < 1:
        record["inactivated_columns"] = 1
    rehash(record)


def make_error(record: dict) -> None:
    record.update({
        "outcome": "error", "success": False, "rank_fail": False,
        "error": True, "weak": True, "binary_deficiency": 0,
        "gf256_rank_gain": 0, "heavy_shortfall": 0,
        "first_rank_fail": -1, "failure_trials": [0],
    })
    rehash(record)


def make_configuration_failure(record: dict) -> None:
    record.update({
        "returncode": 2, "outcome": "configuration_failure",
        "success": False, "rank_fail": False, "error": False,
        "configuration_failure": True, "weak": True,
        "staircase": None, "source_hits": None,
        "inactivated_columns": None, "binary_deficiency": None,
        "gf256_rank_gain": None, "heavy_shortfall": None,
        "first_rank_fail": None, "failure_trials": None,
        "effective_precode_seed": None, "effective_packet_seed": None,
        "block_xors": None, "gf256_muladds": None,
        "configuration_stderr_sha256": "4" * 64,
    })
    rehash(record)


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_partition(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.PACKET_DELTA, 39)
        self.assertEqual(subject.REPETITIONS, (0, 1))
        self.assertEqual(len(subject.CONSULTED_COORDINATES), 47)
        self.assertEqual(len(subject.UNTOUCHED_COORDINATES), 1033)
        self.assertEqual(subject.sha256_json(
            list(subject.UNTOUCHED_COORDINATES)),
            subject.EXPECTED_UNTOUCHED_SHA256)
        candidate = subject._contract_body()["candidate"]
        self.assertTrue(candidate["fallback_delta225_forbidden"])
        self.assertFalse(candidate["production_split_pq_implemented"])
        self.assertIn("no additional offset map",
                      candidate["production_complexity_condition"])
        self.assertIn("63,999-byte",
                      candidate["preexisting_whole_profile_footprint"])
        work = subject._contract_body()["gate"]["work"]
        self.assertEqual(work["accepted_historical_xor_ratio_max"],
                         "0.9829453304")
        self.assertFalse(work["old_five_percent_threshold_applied"])

    def test_full_two_repetition_roster_and_delta39_attempts(self) -> None:
        invocations = subject.make_invocations()
        self.assertEqual(len(invocations), 2160)
        self.assertEqual(
            (invocations[0].coordinate_ordinal, invocations[0].repetition),
            (0, 0))
        self.assertEqual(
            (invocations[-1].coordinate_ordinal,
             invocations[-1].repetition), (1079, 1))
        self.assertEqual(invocations[0].coordinate.attempt, 9)
        self.assertEqual(invocations[0].coordinate.packet_attempt, 48)
        K6 = subject.K_VALUES.index(6) * 36
        self.assertEqual(subject.Coordinate(K6).packet_attempt, 88)
        self.assertEqual(subject.stage0._packet_seed_for_offset(
            "0x12345678", 255, 39), "0xf52e28a7")

    def test_parser_requires_every_sealed_directory(self) -> None:
        parser = subject._parser()
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args(["run", "--bench", "b", "--output-dir", "o"])
        parsed = parser.parse_args([
            "run", "--bench", "b", "--baseline-dir", "base",
            "--v9-dir", "v9", "--stage0-dir", "s0",
            "--stage1-dir", "s1", "--stage2-dir", "s2",
            "--output-dir", "o"])
        self.assertEqual(parsed.stage2_dir, Path("s2"))

    def test_stage2_exact_artifact_pins_are_frozen(self) -> None:
        self.assertEqual(subject.STAGE2_FILES[subject.stage2.SUMMARY_NAME], (
            50540,
            "6940c4c0388fc535c0acb5593b58ce288f8bda20bab1cc47acd92dcde39b7c95"))
        self.assertEqual(subject.EXPECTED_STAGE2_RECORD_COUNT, 9639)
        self.assertEqual(subject.EXPECTED_STAGE2_RAW_DECODED_BYTES, 14955348)


class ParserAndRawTests(unittest.TestCase):
    def test_delta39_parser_and_exact_raw_reparse(self) -> None:
        invocation = subject.Invocation(4, 2, 0)
        stdout = stdout_for(invocation)
        baseline = parser_baseline(invocation.coordinate_ordinal)
        record, raw = subject.parse_process_result(
            process_result(invocation, stdout), CONTROLLER_COMMIT,
            SOURCE_SHA, baseline)
        self.assertEqual(record["delta"], 39)
        self.assertEqual(record["packet_attempt"], 48)
        self.assertEqual(record["stratum"], "consulted")
        self.assertEqual(subject.validate_raw_reparse(
            raw, record, CONTROLLER_COMMIT, SOURCE_SHA, baseline),
            len(stdout))

    def test_tampered_raw_stdout_fails_closed(self) -> None:
        invocation = subject.Invocation(4, 2, 0)
        baseline = parser_baseline(2)
        record, raw = subject.parse_process_result(
            process_result(invocation, stdout_for(invocation)),
            CONTROLLER_COMMIT, SOURCE_SHA, baseline)
        raw = copy.deepcopy(raw)
        decoded = subject._decode_raw(
            raw["stdout_base64"], subject.MAX_STDOUT_BYTES, "stdout")
        raw["stdout_base64"] = subject._canonical_b64(
            decoded.replace(b"packet_attempt,", b"packet_attempz,", 1))
        with self.assertRaises(subject.Stage3Error):
            subject.validate_raw_reparse(
                raw, record, CONTROLLER_COMMIT, SOURCE_SHA, baseline)

    def test_wrong_effective_packet_seed_is_rejected(self) -> None:
        invocation = subject.Invocation(4, 2, 0)
        stdout = stdout_for(invocation)
        expected = csv_values(invocation)["effective_packet_seed"].encode()
        with self.assertRaisesRegex(subject.Stage3Error, "seed"):
            subject.parse_process_result(
                process_result(invocation, stdout.replace(
                    expected, b"0x00000001")), CONTROLLER_COMMIT,
                SOURCE_SHA, parser_baseline(2))

    def test_run_raw_executes_pinned_descriptor_in_frozen_environment(self) \
            -> None:
        invocation = subject.Invocation(0, 0, 0)
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        completed = mock.Mock(stdout=b"stdout\n", stderr=b"", returncode=0)
        with mock.patch.object(subject.subprocess, "run",
                               return_value=completed) as run:
            result = subject._run_raw(invocation, pinned)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(run.call_args.kwargs["executable"],
                         "/proc/self/fd/17")
        self.assertEqual(run.call_args.kwargs["pass_fds"], (17,))
        self.assertEqual(run.call_args.kwargs["env"],
                         subject.CHILD_ENVIRONMENT)

    def test_execution_aborts_on_first_configuration_failure(self) -> None:
        result = mock.Mock()
        historical = SimpleNamespace(baseline=full_baseline())
        with mock.patch.object(subject, "_run_raw", return_value=result) as run, \
                mock.patch.object(subject, "parse_process_result",
                                  return_value=(
                                      {"configuration_failure": True}, {})), \
                self.assertRaisesRegex(subject.Stage3Error,
                                       "before delta39 packet q"):
            subject._execute_roster(
                mock.Mock(), CONTROLLER_COMMIT, SOURCE_SHA, historical)
        run.assert_called_once()


class GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.records, cls.evidence = passing_roster()

    def test_exact_no_regression_roster_passes(self) -> None:
        verdict = subject.adjudicate(self.records, self.evidence)
        self.assertEqual(verdict["disposition"], "PASS")
        self.assertEqual(verdict["candidate_packet_delta"], 39)
        self.assertEqual(verdict["repeat_exact_coordinate_count"], 1080)
        self.assertEqual(
            verdict["stage2_confirmation_projection_exact_observation_count"],
            94)
        self.assertEqual(verdict["selection_leakage"], {
            "consulted_coordinate_count": 47,
            "untouched_coordinate_count": 1033,
            "next_holdout_must_be_disjoint": True,
        })
        self.assertEqual(
            verdict["strata"]["consulted"]["coordinate_count"], 47)
        self.assertEqual(
            verdict["strata"]["untouched"]["coordinate_count"], 1033)
        self.assertFalse(verdict["production_split_pq_implemented"])
        self.assertFalse(verdict["fixed_offset_family_retired"])

    def test_untouched_introduction_rejects_both_repetitions(self) -> None:
        records = list(self.records)
        coordinate = 0
        for repetition in subject.REPETITIONS:
            index = coordinate * 2 + repetition
            records[index] = copy.deepcopy(records[index])
            make_rank_failure(records[index])
        verdict = subject.adjudicate(records, self.evidence)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertFalse(verdict["gates"]["untouched_1033_zero_weak"])
        self.assertFalse(verdict["gates"]["introductions_vs_pair01_zero"])
        self.assertTrue(verdict["fixed_offset_family_retired"])

    def test_consulted_weak_and_unrepaired_pair01_weak_reject(self) -> None:
        for coordinate, gate in (
                (subject.CONSULTED_COORDINATES[0],
                 "consulted_47_zero_weak"),
                (subject.BASELINE_WEAK_COORDINATES[0],
                 "all_three_pair01_weak_coordinates_repaired")):
            with self.subTest(coordinate=coordinate):
                records = list(self.records)
                for repetition in subject.REPETITIONS:
                    index = coordinate * 2 + repetition
                    records[index] = copy.deepcopy(records[index])
                    make_rank_failure(records[index])
                verdict = subject.adjudicate(records, self.evidence)
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate])

    def test_error_and_configuration_failure_are_defensive_rejects(self) \
            -> None:
        for mutator, gate in (
                (make_error, "candidate_errors_zero"),
                (make_configuration_failure,
                 "candidate_configuration_failures_zero")):
            with self.subTest(gate=gate):
                records = list(self.records)
                coordinate = subject.UNTOUCHED_COORDINATES[0]
                for repetition in subject.REPETITIONS:
                    index = coordinate * 2 + repetition
                    records[index] = copy.deepcopy(records[index])
                    mutator(records[index])
                verdict = subject.adjudicate(records, self.evidence)
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate])

    def test_identity_and_precode_seed_tampering_fail_invalid(self) -> None:
        for field, value in (("stratum", "consulted"),
                             ("effective_precode_seed",
                              "0x0000000000000001")):
            with self.subTest(field=field):
                records = list(self.records)
                records[0] = copy.deepcopy(records[0])
                records[0][field] = value
                rehash(records[0])
                with self.assertRaises(subject.Stage3Error):
                    subject.adjudicate(records, self.evidence)

    def test_consulted_projection_mismatch_is_a_reject(self) -> None:
        confirmation = dict(self.evidence.stage2_confirmation_by_coordinate)
        coordinate = subject.CONSULTED_COORDINATES[0]
        confirmation[coordinate] = dict(confirmation[coordinate])
        confirmation[coordinate]["block_xors"] += 1
        evidence = SimpleNamespace(
            baseline=self.evidence.baseline,
            stage2_confirmation_by_coordinate=confirmation)
        verdict = subject.adjudicate(self.records, evidence)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(
            verdict[
                "stage2_confirmation_projection_mismatch_coordinate_ordinals"],
            [coordinate])
        self.assertFalse(verdict["gates"][
            "consulted_94_observations_equal_stage2_confirmation"])

    def test_each_full_aggregate_work_metric_is_a_gate(self) -> None:
        for field, gate in (
                ("block_xors", "aggregate_block_xors_not_above_pair01"),
                ("gf256_muladds",
                 "aggregate_gf256_muladds_not_above_pair01"),
                ("inactivated_columns",
                 "aggregate_inactivated_columns_not_above_pair01")):
            with self.subTest(field=field):
                records = list(self.records)
                for repetition in subject.REPETITIONS:
                    records[repetition] = copy.deepcopy(records[repetition])
                    records[repetition][field] += 1
                    rehash(records[repetition])
                verdict = subject.adjudicate(records, self.evidence)
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate])

    def test_repeat_mismatch_and_roster_changes_fail_closed(self) -> None:
        records = list(self.records)
        records[1] = copy.deepcopy(records[1])
        records[1]["block_xors"] += 1
        rehash(records[1])
        verdict = subject.adjudicate(records, self.evidence)
        self.assertFalse(verdict["gates"][
            "full_repeat_scientific_projection_exact"])
        self.assertFalse(verdict["gates"]["full_repeat_work_exact"])
        for rows in (self.records[:-1], self.records + [self.records[-1]],
                     [self.records[1], self.records[0], *self.records[2:]]):
            with self.assertRaises(subject.Stage3Error):
                subject.adjudicate(rows, self.evidence)


class TransportTests(unittest.TestCase):
    def test_stable_snapshot_detects_same_path_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "input"
            path.write_bytes(b"original")
            snapshots = subject._snapshot_group(
                root, {"input": (8, hashlib.sha256(b"original").hexdigest())},
                "fixture")
            subject._assert_snapshots(snapshots)
            replacement = root / "replacement"
            replacement.write_bytes(b"mutated!")
            os.replace(str(replacement), str(path))
            with self.assertRaisesRegex(subject.Stage3Error,
                                        "historical input changed"):
                subject._assert_snapshots(snapshots)

    def test_exact_consumed_bytes_are_hashed_not_only_fd_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summary_path = root / "summary"
            summary_bytes = b'{"a":1}\n'
            summary_path.write_bytes(summary_bytes)
            with subject.stage1._open_pinned_file(
                    summary_path, len(summary_bytes),
                    hashlib.sha256(summary_bytes).hexdigest(),
                    "summary fixture") as pin:
                self.assertEqual(subject._read_pinned_bytes_exact(
                    pin, len(summary_bytes),
                    hashlib.sha256(summary_bytes).hexdigest(), 100,
                    "summary fixture"), summary_bytes)
                with mock.patch.object(
                        pin, "read_bytes", return_value=b'{"b":1}\n'), \
                        self.assertRaisesRegex(subject.Stage3Error,
                                               "consumed-byte"):
                    subject._read_pinned_bytes_exact(
                        pin, len(summary_bytes),
                        hashlib.sha256(summary_bytes).hexdigest(), 100,
                        "summary fixture")

            stream_path = root / "stream"
            stream_bytes = b'{"ordinal":0,"x":1}\n'
            stream_path.write_bytes(stream_bytes)
            with subject.stage1._open_pinned_file(
                    stream_path, len(stream_bytes),
                    hashlib.sha256(stream_bytes).hexdigest(),
                    "stream fixture") as pin:
                rows = subject._read_pinned_jsonl_exact(
                    pin, len(stream_bytes),
                    hashlib.sha256(stream_bytes).hexdigest(),
                    "stream fixture")
                self.assertEqual(rows, ({"ordinal": 0, "x": 1},))
                transient = b'{"ordinal":0,"x":2}\n'
                with mock.patch.object(
                        pin, "lines", return_value=iter(((0, transient),))), \
                        self.assertRaisesRegex(subject.Stage3Error,
                                               "consumed-byte"):
                    subject._read_pinned_jsonl_exact(
                        pin, len(stream_bytes),
                        hashlib.sha256(stream_bytes).hexdigest(),
                        "stream fixture")

    def test_exact_consumed_reader_rejects_fd_identity_drift(self) -> None:
        data = b'{"a":1}\n'
        pin = SimpleNamespace(descriptor=7, read_bytes=mock.Mock(
            return_value=data))
        before = SimpleNamespace(
            st_dev=1, st_ino=2, st_size=len(data), st_mtime_ns=3,
            st_ctime_ns=4)
        after = SimpleNamespace(
            st_dev=1, st_ino=2, st_size=len(data), st_mtime_ns=3,
            st_ctime_ns=5)
        with mock.patch.object(subject.os, "fstat",
                               side_effect=(before, after)), \
                self.assertRaisesRegex(subject.Stage3Error,
                                       "consumed-byte"):
            subject._read_pinned_bytes_exact(
                pin, len(data), hashlib.sha256(data).hexdigest(), 100,
                "identity fixture")

    def test_historical_preflight_failure_issues_no_work(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(
                    subject, "load_historical_evidence",
                    side_effect=subject.Stage3Error("preflight failed")), \
                mock.patch.object(subject, "_source_receipt") as source, \
                mock.patch.object(subject.screen, "_open_binary") as opened, \
                mock.patch.object(subject, "_run_raw") as run:
            output = Path(temporary) / "fresh-output"
            with self.assertRaisesRegex(subject.Stage3Error,
                                        "preflight failed"):
                subject.run_stage3(
                    Path("bench"), Path("base"), Path("v9"), Path("s0"),
                    Path("s1"), Path("s2"), output)
        source.assert_not_called()
        opened.assert_not_called()
        run.assert_not_called()

    def test_jsonl_writer_is_exclusive_bounded_and_canonical(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rows.jsonl"
            digest, count, size = subject._write_jsonl(
                path, ({"b": 2, "a": 1},), 100)
            self.assertEqual(path.read_bytes(), b'{"a":1,"b":2}\n')
            self.assertEqual(digest,
                             hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual((count, size), (1, len(path.read_bytes())))
            with self.assertRaises(subject.Stage3Error):
                subject._write_jsonl(path, ({"a": 1},), 100)
            with self.assertRaises(subject.Stage3Error):
                subject._write_jsonl(
                    Path(temporary) / "large", ({"x": "y"},), 1)

    def test_output_stream_pair_order_is_reparsed(self) -> None:
        invocation = subject.Invocation(4, 2, 0)
        baseline = parser_baseline(2)
        record, raw = subject.parse_process_result(
            process_result(invocation, stdout_for(invocation)),
            CONTROLLER_COMMIT, SOURCE_SHA, baseline)
        evidence = SimpleNamespace(
            baseline=tuple(parser_baseline(index)
                           for index in range(subject.COORDINATE_COUNT)),
            stage2_confirmation_by_coordinate={})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            records_path = root / "records.jsonl"
            raw_path = root / "raw.jsonl"
            records_sha, _, _ = subject._write_jsonl(records_path, (record,))
            raw_sha, _, _ = subject._write_jsonl(raw_path, (raw,))
            expected = {"disposition": "REJECT"}
            with mock.patch.object(subject, "adjudicate",
                                   return_value=expected):
                audit = subject._audit_output_streams(
                    records_path, raw_path, records_sha, raw_sha, 1,
                    CONTROLLER_COMMIT, SOURCE_SHA, evidence)
            self.assertEqual(audit["verdict"], expected)
            self.assertEqual(audit["decoded_bytes"], len(stdout_for(invocation)))

    def test_output_stream_rejects_reversal_and_cardinality_change(self) \
            -> None:
        invocations = (subject.Invocation(4, 2, 0),
                       subject.Invocation(5, 2, 1))
        baseline = parser_baseline(2)
        pairs = [subject.parse_process_result(
            process_result(invocation, stdout_for(invocation)),
            CONTROLLER_COMMIT, SOURCE_SHA, baseline)
                 for invocation in invocations]
        records = [pair[0] for pair in pairs]
        raw = [pair[1] for pair in pairs]
        evidence = SimpleNamespace(
            baseline=tuple(parser_baseline(index)
                           for index in range(subject.COORDINATE_COUNT)),
            stage2_confirmation_by_coordinate={})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            records_path = root / "records.jsonl"
            raw_path = root / "raw.jsonl"
            reversed_path = root / "raw-reversed.jsonl"
            records_sha, _, _ = subject._write_jsonl(records_path, records)
            raw_sha, _, _ = subject._write_jsonl(raw_path, raw)
            reversed_sha, _, _ = subject._write_jsonl(
                reversed_path, reversed(raw))
            expected = {"disposition": "REJECT"}
            with mock.patch.object(subject, "adjudicate",
                                   return_value=expected), \
                    self.assertRaises(subject.Stage3Error):
                subject._audit_output_streams(
                    records_path, reversed_path, records_sha, reversed_sha, 2,
                    CONTROLLER_COMMIT, SOURCE_SHA, evidence)
            with mock.patch.object(subject, "adjudicate",
                                   return_value=expected), \
                    self.assertRaises(subject.Stage3Error):
                subject._audit_output_streams(
                    records_path, raw_path, records_sha, raw_sha, 3,
                    CONTROLLER_COMMIT, SOURCE_SHA, evidence)

    def test_even_rehashed_raw_timing_change_maps_to_changed_parsed_hash(self) \
            -> None:
        invocation = subject.Invocation(4, 2, 0)
        baseline = parser_baseline(2)
        record, raw = subject.parse_process_result(
            process_result(invocation, stdout_for(invocation)),
            CONTROLLER_COMMIT, SOURCE_SHA, baseline)
        raw = copy.deepcopy(raw)
        decoded = subject._decode_raw(
            raw["stdout_base64"], subject.MAX_STDOUT_BYTES, "stdout")
        changed = decoded.replace(
            b",1.000,1.000,1.000,1.000,1.000,1.000,",
            b",2.000,1.000,1.000,1.000,1.000,1.000,", 1)
        self.assertNotEqual(changed, decoded)
        raw["stdout_base64"] = subject._canonical_b64(changed)
        raw["stdout_sha256"] = subject.sha256_bytes(changed)
        altered = copy.deepcopy(record)
        altered["stdout_sha256"] = raw["stdout_sha256"]
        raw["parsed_record_sha256"] = subject.sha256_json(altered)
        self.assertEqual(subject.validate_raw_reparse(
            raw, altered, CONTROLLER_COMMIT, SOURCE_SHA, baseline),
            len(changed))
        with self.assertRaises(subject.Stage3Error):
            subject.validate_raw_reparse(
                raw, record, CONTROLLER_COMMIT, SOURCE_SHA, baseline)


class Stage2ArtifactLoaderTests(unittest.TestCase):
    def _fixture(self) -> tuple[dict, SimpleNamespace, list[dict], list[dict],
                                dict]:
        historical_receipt = {"frozen": "history"}
        historical = SimpleNamespace(
            receipt=mock.Mock(return_value=historical_receipt),
            baseline_by_coordinate={
                coordinate: parser_baseline(coordinate)
                for coordinate in subject.CONSULTED_COORDINATES
            })
        source = {
            "source_git_commit": subject.EXPECTED_STAGE2_SOURCE_COMMIT,
            "source_receipt_sha256":
                subject.EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256,
        }
        verdict = {
            "disposition": "PASS", "selected_delta": 39,
            "zero_weak_deltas": [39, 225], "gates": {"all": True},
        }
        summary = {
            "schema": subject.stage2.SUMMARY_SCHEMA,
            "contract_sha256": subject.EXPECTED_STAGE2_CONTRACT_SHA256,
            "summary_sha256": subject.EXPECTED_STAGE2_SUMMARY_SELF_SHA256,
            "records_file": subject.stage2.RECORD_NAME,
            "records_sha256": subject.STAGE2_FILES[
                subject.stage2.RECORD_NAME][1],
            "record_count": 47,
            "raw_records_file": subject.stage2.RAW_NAME,
            "raw_records_sha256": subject.STAGE2_FILES[
                subject.stage2.RAW_NAME][1],
            "raw_record_count": 47,
            "raw_decoded_bytes": 47,
            "raw_file_bytes": subject.STAGE2_FILES[
                subject.stage2.RAW_NAME][0],
            "raw_encoding": "canonical-base64-exact-bytes",
            "raw_reparse_equality": True,
            "stage2_only": True, "promotion_evidence": False,
            "fresh_roots_used": False, "timing_evidence": False,
            "source_receipt": source,
            "bench": {
                "sha256": subject.EXPECTED_STAGE0_BENCH_SHA256,
                "size": subject.EXPECTED_STAGE0_BENCH_SIZE,
            },
            "bench_description": {
                "schema": subject.stage2.BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT,
            },
            "historical_evidence": historical_receipt,
            **verdict,
        }
        records = [{
            "ordinal": index,
            "coordinate_ordinal": coordinate,
            "delta": 39,
            "phase": subject.stage2.PHASE_CONFIRMATION,
            "success": True,
        } for index, coordinate in enumerate(subject.CONSULTED_COORDINATES)]
        raw = [{"ordinal": index} for index in range(47)]
        return summary, historical, records, raw, verdict

    def test_full_mocked_stage2_reauth_reparses_and_adjudicates(self) -> None:
        summary, historical, records, raw, verdict = self._fixture()
        pins = {
            subject.stage2.SUMMARY_NAME: object(),
            subject.stage2.RECORD_NAME: object(),
            subject.stage2.RAW_NAME: object(),
        }

        def self_hash(value: dict, field: str) -> str:
            del field
            if value is summary:
                return subject.EXPECTED_STAGE2_SUMMARY_SELF_SHA256
            return subject.EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256

        with mock.patch.object(subject, "EXPECTED_STAGE2_RECORD_COUNT", 47), \
                mock.patch.object(
                    subject, "EXPECTED_STAGE2_RAW_DECODED_BYTES", 47), \
                mock.patch.object(
                    subject.stage2, "PRE_CONFIRMATION_OBSERVATION_COUNT", 0), \
                mock.patch.object(subject, "_read_pinned_bytes_exact",
                                  return_value=b"{}\n"), \
                mock.patch.object(
                    subject.stage2, "_parse_canonical_json_bytes",
                    return_value=summary), \
                mock.patch.object(
                    subject, "_read_pinned_jsonl_exact",
                    side_effect=(tuple(records), tuple(raw))), \
                mock.patch.object(subject.stage2, "self_hash",
                                  side_effect=self_hash), \
                mock.patch.object(
                    subject.stage2, "_invocation_from_record",
                    side_effect=lambda row: SimpleNamespace(
                        coordinate_ordinal=row["coordinate_ordinal"])), \
                mock.patch.object(subject.stage2, "_validate_record") as valid, \
                mock.patch.object(subject.stage2, "validate_raw_reparse",
                                  return_value=1) as reparse, \
                mock.patch.object(subject.stage2, "adjudicate",
                                  return_value=verdict) as adjudicate, \
                mock.patch.object(
                    subject.stage2, "scientific_projection",
                    side_effect=lambda row: {
                        "coordinate_ordinal": row["coordinate_ordinal"]}):
            confirmation, receipt = subject._load_stage2_artifact(
                pins, historical)
        self.assertEqual(tuple(sorted(confirmation)),
                         subject.CONSULTED_COORDINATES)
        self.assertEqual(receipt["selected_delta"], 39)
        self.assertEqual(valid.call_count, 47)
        self.assertEqual(reparse.call_count, 47)
        adjudicate.assert_called_once()
        for index, call in enumerate(reparse.call_args_list):
            self.assertIs(call.args[0], raw[index])
            self.assertIs(call.args[1], records[index])

    def test_stage2_summary_source_or_history_mismatch_stops_pre_stream(self) \
            -> None:
        for field in ("historical_evidence", "source_receipt"):
            with self.subTest(field=field):
                summary, historical, _, _, _ = self._fixture()
                summary[field] = {"tampered": True}
                pins = {subject.stage2.SUMMARY_NAME: object()}

                def self_hash(value: dict, hash_field: str) -> str:
                    del hash_field
                    if value is summary:
                        return subject.EXPECTED_STAGE2_SUMMARY_SELF_SHA256
                    return subject.EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256

                with mock.patch.object(
                        subject, "EXPECTED_STAGE2_RECORD_COUNT", 47), \
                        mock.patch.object(
                            subject, "EXPECTED_STAGE2_RAW_DECODED_BYTES", 47), \
                        mock.patch.object(
                            subject, "_read_pinned_bytes_exact",
                            return_value=b"{}\n"), \
                        mock.patch.object(
                            subject.stage2, "_parse_canonical_json_bytes",
                            return_value=summary), \
                        mock.patch.object(subject.stage2, "self_hash",
                                          side_effect=self_hash), \
                        mock.patch.object(
                            subject, "_read_pinned_jsonl_exact") as streams, \
                        self.assertRaises(subject.Stage3Error):
                    subject._load_stage2_artifact(pins, historical)
                streams.assert_not_called()


if __name__ == "__main__":
    unittest.main()
