#!/usr/bin/env python3
"""Pure mocked tests for the sealed MIX2 pair12 Stage 1 controller."""

from __future__ import annotations

import copy
from contextlib import redirect_stderr
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_pair_stage1 as subject


SOURCE_COMMIT = "1" * 40
BENCH_SHA = "2" * 64
SOURCE_SHA = "3" * 64


def baseline_row(coordinate_ordinal: int, success: bool = True,
                 block_xors: int = 100, gf256_muladds: int = 50,
                 inactivated: int = 20) -> dict:
    coordinate = subject.Coordinate(coordinate_ordinal)
    deficiency = 10 if success else 12
    gain = 10 if success else 11
    return {
        "schema": subject.crossfit.RECORD_SCHEMA,
        "record_ordinal": coordinate.baseline_record_ordinal,
        "K": coordinate.K,
        "attempt": coordinate.attempt,
        "root_index": coordinate.root_index,
        "loss_root": coordinate.root,
        "schedule_index": coordinate.schedule_index,
        "schedule": coordinate.schedule,
        "cell_ordinal": coordinate.cell_ordinal,
        "trace_executed": True,
        "outcome": "success" if success else "rank_fail",
        "success": success,
        "rank_fail": not success,
        "error": False,
        "binary_deficiency": deficiency,
        "gf256_rank_gain": gain,
        "inactivated_columns": inactivated,
        "heavy_shortfall": 1 if not success else 0,
        "effective_precode_seed": "0x123456789abcdef0",
        "effective_packet_seed": "0x12345678",
        "block_xors": block_xors,
        "gf256_muladds": gf256_muladds,
        "command_sha256": "4" * 64,
        "stdout_sha256": "5" * 64,
        "configuration_proof_sha256": None,
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
        "exact_packet_attempt": str(coordinate.attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": SOURCE_COMMIT, "mix_pair": "12",
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(field, values[field])
        for field in subject.METADATA_ORDER)


def csv_values(invocation: subject.Invocation, outcome: str = "success",
               block_xors: int = 100, gf256_muladds: int = 50) -> dict:
    if outcome == "success":
        success, rank_fail, error, deficiency, gain = 1, 0, 0, 10, 10
    elif outcome == "rank_fail":
        success, rank_fail, error, deficiency, gain = 0, 1, 0, 12, 11
    elif outcome == "error":
        success, rank_fail, error, deficiency, gain = 0, 0, 1, 0, 0
    else:
        raise AssertionError("unknown outcome")
    weak = success == 0
    coordinate = invocation.coordinate
    return {
        "N": str(coordinate.K), "bb": "2", "heavy_family": "periodic",
        "mix_count": "2",
        "staircase": str(subject.STAIRCASE_BY_K[coordinate.K]),
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
        "seed_attempt": "",
        "block_xors_mu": "{}.000".format(block_xors),
        "block_muladds_mu": "{}.000".format(gf256_muladds),
        "first_rank_fail": "0" if rank_fail else "-1",
        "binary_def_hist": "{}:1".format(deficiency),
        "heavy_gain_hist": "{}:1".format(gain),
        "failure_trials": "0" if weak else "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(coordinate.attempt),
        "packet_attempt": str(coordinate.attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed": "0x123456789abcdef0",
        "effective_packet_seed": "0x12345678",
    }


def stdout_for(invocation: subject.Invocation, outcome: str = "success",
               block_xors: int = 100, gf256_muladds: int = 50) -> bytes:
    values = csv_values(invocation, outcome, block_xors, gf256_muladds)
    return ("\n".join((
        metadata_line(invocation), ",".join(subject.CSV_HEADER),
        ",".join(values[field] for field in subject.CSV_HEADER))) +
        "\n").encode("ascii")


def process_result(invocation: subject.Invocation, stdout: bytes,
                   returncode: int = 0,
                   stderr: bytes = b"") -> subject.ProcessResult:
    return subject.ProcessResult(
        invocation=invocation,
        command_sha256=subject.sha256_json({
            "argv": invocation.argv(Path("wirehair_v2_bench"))[1:],
        }),
        stdout_sha256=hashlib.sha256(stdout).hexdigest(),
        stderr_sha256=hashlib.sha256(stderr).hexdigest(),
        returncode=returncode, stdout=stdout, stderr=stderr)


def parsed_record(invocation: subject.Invocation,
                  baseline: dict | None = None,
                  outcome: str = "success") -> dict:
    if baseline is None:
        baseline = baseline_row(invocation.coordinate_ordinal)
    return dict(subject.parse_process_result(
        process_result(invocation, stdout_for(invocation, outcome)),
        SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA, baseline))


def synthetic_record(invocation: subject.Invocation, baseline: dict) -> dict:
    return parsed_record(invocation, baseline)


def successful_roster(baseline: list[dict]) -> list[dict]:
    return [synthetic_record(invocation,
                             baseline[invocation.coordinate_ordinal])
            for invocation in subject.make_invocations()]


def rehash(record: dict) -> None:
    record["deterministic_projection_sha256"] = subject.sha256_json(
        subject.deterministic_projection(record))


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_exact_rosters(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.COORDINATE_COUNT, 1080)
        self.assertEqual(len(subject.make_invocations()), 2160)
        self.assertEqual(
            subject.sha256_json([[K, subject.ATTEMPT_BY_K[K]]
                                 for K in subject.K_VALUES]),
            subject.EXPECTED_ATTEMPT_MAP_SHA256)
        self.assertFalse(set(subject.ROOTS) &
                         set(subject.screen.FINAL_VALIDATION_ROOTS))

    def test_exact_pair12_full_payload_command_and_ordinals(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        self.assertEqual(invocation.argv(Path("/pinned/bench")), [
            "/pinned/bench", "precodefail", "--N", "2", "--bb-list",
            "2", "--overhead", "0", "--trials", "1", "--threads", "1",
            "--loss", "0.5", "--seed", "0xd1b54a32d192ed03",
            "--schedule", "burst", "--heavy-family", "periodic",
            "--mix-count", "2", "--mix-pair", "12",
            "--binary-dense-rows", "12", "--gf256-heavy-rows", "12",
            "--dense-anchors", "two07", "--paired-overhead-stream",
            "--full-payload-solve", "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", "9", "--exact-packet-attempt", "9",
            "--construction-seed-basis", "production-profile",
        ])
        last = subject.Coordinate(1079)
        self.assertEqual(last.baseline_record_ordinal, 267299)
        self.assertEqual(subject.Coordinate(0).baseline_record_ordinal, 324)
        with self.assertRaises(subject.Stage1Error):
            subject.Invocation(1, 0, 0)

    def test_run_parser_requires_all_explicit_history_directories(self) -> None:
        parser = subject._parser()
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args([
                "run", "--bench", "b", "--output-dir", "o"])
        parsed = parser.parse_args([
            "run", "--bench", "b", "--baseline-dir", "base",
            "--v9-dir", "v9", "--stage0-dir", "s0", "--output-dir", "o",
        ])
        self.assertEqual(parsed.baseline_dir, Path("base"))


class ParserTests(unittest.TestCase):
    def test_accepts_exact_metadata_csv_and_baseline_receipt(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        baseline = baseline_row(0)
        record = parsed_record(invocation, baseline)
        self.assertTrue(record["success"])
        self.assertEqual(record["mix_pair"], "12")
        self.assertEqual(record["pair01_baseline"],
                         subject._baseline_projection(baseline))
        self.assertEqual(record["block_xors"], 100)
        self.assertEqual(record["gf256_muladds"], 50)

    def test_rejects_missing_reordered_or_wrong_pair_metadata(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        valid = stdout_for(invocation).decode("ascii").splitlines()
        cases = [
            (" ".join(valid[0].split(" ")[:-1]) + "\n" +
             "\n".join(valid[1:]) + "\n").encode("ascii"),
            stdout_for(invocation).replace(b"mix_pair=12", b"mix_pair=02"),
        ]
        tokens = valid[0].split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        cases.append((" ".join(tokens) + "\n" +
                      "\n".join(valid[1:]) + "\n").encode("ascii"))
        for stdout in cases:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.Stage1Error):
                subject.parse_process_result(
                    process_result(invocation, stdout), SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA, baseline_row(0))

    def test_configuration_failure_transport_is_exact_but_weak(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=9 result=2\n").encode("ascii")
        record = subject.parse_process_result(
            process_result(invocation, stdout, 2, stderr), SOURCE_COMMIT,
            BENCH_SHA, SOURCE_SHA, baseline_row(0))
        self.assertTrue(record["configuration_failure"])
        self.assertTrue(record["weak"])
        with self.assertRaises(subject.Stage1Error):
            subject.parse_process_result(
                process_result(invocation, stdout, 2, stderr + b"extra\n"),
                SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA, baseline_row(0))

    def test_rejects_noncanonical_work_and_wrong_baseline_coordinate(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        stdout = stdout_for(invocation).replace(
            b",100.000,50.000,", b",100.00,50.000,")
        with self.assertRaises(subject.Stage1Error):
            subject.parse_process_result(
                process_result(invocation, stdout), SOURCE_COMMIT,
                BENCH_SHA, SOURCE_SHA, baseline_row(0))
        with self.assertRaises(subject.Stage1Error):
            parsed_record(invocation, baseline_row(1))


class HistoricalTransportTests(unittest.TestCase):
    def test_matrix_row_validator_and_projection_field_roster(self) -> None:
        row = baseline_row(0)
        subject._validate_matrix_row(row, 324)
        self.assertEqual(tuple(subject._baseline_projection(row)),
                         subject.BASELINE_PROJECTION_FIELDS)
        damaged = dict(row)
        damaged["record_ordinal"] = 325
        with self.assertRaises(subject.Stage1Error):
            subject._validate_matrix_row(damaged, 324)

    def test_pinned_file_rejects_symlink_and_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "history.json"
            path.write_bytes(b"{}\n")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with subject._open_pinned_file(
                    path, 3, digest, "history") as pinned:
                self.assertEqual(pinned.read_bytes(10), b"{}\n")
                path.write_bytes(b"[]\n")
                with self.assertRaises(subject.Stage1Error):
                    pinned.receipt()
            target = root / "target"
            target.write_bytes(b"{}\n")
            link = root / "link"
            link.symlink_to(target)
            with self.assertRaises(subject.Stage1Error):
                subject._open_pinned_file(link, 3, digest, "link")

    def test_canonical_json_transport_rejects_noncanonical_forms(self) -> None:
        self.assertEqual(subject._parse_canonical_json_bytes(
            b'{"a":1}\n', "object"), {"a": 1})
        for value in (b'{ "a":1}\n', b'{"a":1}', b'{"a":1}\r\n',
                      b'{"a":1}\n\n'):
            with self.subTest(value=value), \
                    self.assertRaises(subject.Stage1Error):
                subject._parse_canonical_json_bytes(value, "object")

    def test_v9_projection_matches_only_exact_deterministic_fields(self) -> None:
        base = baseline_row(0)
        candidate = {
            "success": 1, "rank_fail": 0, "error": 0,
            "inactivated_columns": 20, "block_xors": 100,
            "gf256_muladds": 50,
            "effective_precode_seed": "0x123456789abcdef0",
            "effective_packet_seed": "0x12345678",
        }
        self.assertEqual(subject._v9_candidate_projection(candidate),
                         subject._baseline_v9_projection(base))
        candidate["block_xors"] += 1
        self.assertNotEqual(subject._v9_candidate_projection(candidate),
                            subject._baseline_v9_projection(base))


class GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.baseline = [baseline_row(index) for index in range(1080)]
        for index in (0, 500, 1079):
            cls.baseline[index] = baseline_row(index, success=False)
        cls.records = successful_roster(cls.baseline)

    def fresh(self) -> tuple[list[dict], list[dict]]:
        return copy.deepcopy(self.records), copy.deepcopy(self.baseline)

    def test_pass_repairs_all_weak_and_deduplicates_aggregates(self) -> None:
        records, baseline = self.fresh()
        verdict = subject.adjudicate(records, baseline)
        self.assertEqual(verdict["disposition"], "PASS")
        self.assertEqual(len(verdict["baseline_weak_coordinates"]), 3)
        self.assertEqual(len(verdict["repaired_coordinates"]), 3)
        self.assertEqual(verdict["candidate_weak_coordinates"], [])
        self.assertEqual(
            verdict["aggregates"]["pair12_candidate"]["block_xors"],
            108000)
        self.assertEqual(
            verdict["aggregates"]["candidate_to_baseline_ratios"][
                "block_xors"], "1.000000000000")

    def test_introduction_and_repeat_mismatch_reject(self) -> None:
        records, baseline = self.fresh()
        for repetition in (0, 1):
            row = records[2 * 10 + repetition]
            row.update({
                "outcome": "rank_fail", "success": False,
                "rank_fail": True, "weak": True,
                "binary_deficiency": 12, "gf256_rank_gain": 11,
                "heavy_shortfall": 1, "first_rank_fail": 0,
                "failure_trials": [0],
            })
            rehash(row)
        records[3]["block_xors"] += 1
        rehash(records[3])
        verdict = subject.adjudicate(records, baseline)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(len(verdict["introduced_weak_coordinates"]), 1)
        self.assertFalse(
            verdict["gates"]["full_repeat_deterministic_projection_exact"])

    def test_deduplicated_work_regression_rejects(self) -> None:
        records, baseline = self.fresh()
        for repetition in (0, 1):
            records[repetition]["block_xors"] += 1
            rehash(records[repetition])
        verdict = subject.adjudicate(records, baseline)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertFalse(
            verdict["gates"]["aggregate_block_xors_not_above_pair01"])
        self.assertEqual(
            verdict["aggregates"]["pair12_candidate"]["block_xors"],
            108001)

    def test_record_to_baseline_receipt_mismatch_is_invalid(self) -> None:
        records, baseline = self.fresh()
        records[0]["pair01_baseline"]["block_xors"] += 1
        rehash(records[0])
        with self.assertRaises(subject.Stage1Error):
            subject.adjudicate(records, baseline)

        records, baseline = self.fresh()
        records[0]["K"] = 3
        rehash(records[0])
        with self.assertRaises(subject.Stage1Error):
            subject.adjudicate(records, baseline)


if __name__ == "__main__":
    unittest.main()
