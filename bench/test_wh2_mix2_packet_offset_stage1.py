#!/usr/bin/env python3
"""Pure mocked tests for the sealed MIX2 packet-offset Stage 1 controller."""

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
import wh2_mix2_packet_offset_stage1 as subject


CONTROLLER_SOURCE_COMMIT = "1" * 40
BENCH_SOURCE_COMMIT = subject.EXPECTED_STAGE0_SOURCE_COMMIT
BENCH_SHA = subject.EXPECTED_STAGE0_BENCH_SHA256
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
        "exact_packet_attempt": str(coordinate.packet_attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": BENCH_SOURCE_COMMIT, "mix_pair": "01",
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
        "packet_attempt": str(coordinate.packet_attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed": "0x123456789abcdef0",
        "effective_packet_seed": subject.stage0._packet_seed_for_offset(
            "0x12345678", coordinate.attempt, subject.PACKET_DELTA),
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
        CONTROLLER_SOURCE_COMMIT, BENCH_SOURCE_COMMIT, BENCH_SHA,
        SOURCE_SHA, baseline))


def synthetic_record(invocation: subject.Invocation, baseline: dict) -> dict:
    return parsed_record(invocation, baseline)


def successful_roster(baseline: list[dict]) -> list[dict]:
    return [synthetic_record(invocation,
                             baseline[invocation.coordinate_ordinal])
            for invocation in subject.make_invocations()]


def rehash(record: dict) -> None:
    record["deterministic_projection_sha256"] = subject.sha256_json(
        subject.deterministic_projection(record))


def source_receipt() -> dict:
    return {
        "source_git_commit": CONTROLLER_SOURCE_COMMIT,
        "source_receipt_sha256": SOURCE_SHA,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "working_source_bytes_equal_HEAD_blobs": True,
        "clean_status_scope": [], "source_files": {},
    }


def bench_receipt() -> dict:
    return {
        "path": "/pinned/bench", "sha256": BENCH_SHA,
        "size": subject.EXPECTED_STAGE0_BENCH_SIZE,
    }


class FakeHistorical:
    def __init__(self, baseline: list[dict]) -> None:
        self.baseline = baseline
        self.receipt_value = {"files": {}, "baseline": {}, "v9": {},
                              "stage0": {}}

    def receipt(self) -> dict:
        return self.receipt_value


class MemoryPin:
    def __init__(self, data: bytes) -> None:
        self.data = data

    def read_bytes(self, maximum: int) -> bytes:
        if len(self.data) > maximum:
            raise AssertionError("synthetic pin exceeds test bound")
        return self.data

    def lines(self):
        for ordinal, line in enumerate(self.data.splitlines(keepends=True)):
            yield ordinal, line


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
        self.assertIn(
            "benchmark test hook only",
            subject._contract_body()["candidate"]["implementation_scope"])
        diagonal = subject._contract_body()["domain"]["current_diagonal_arm"]
        self.assertEqual(diagonal["observation_count"], 0)
        self.assertIn(
            "not a current 1080-coordinate diagonal replay",
            diagonal["stage0_control_scope"])

    def test_exact_delta2_full_payload_command_and_ordinals(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        self.assertEqual(invocation.argv(Path("/pinned/bench")), [
            "/pinned/bench", "precodefail", "--N", "2", "--bb-list",
            "2", "--overhead", "0", "--trials", "1", "--threads", "1",
            "--loss", "0.5", "--seed", "0xd1b54a32d192ed03",
            "--schedule", "burst", "--heavy-family", "periodic",
            "--mix-count", "2", "--mix-pair", "01",
            "--binary-dense-rows", "12", "--gf256-heavy-rows", "12",
            "--dense-anchors", "two07", "--paired-overhead-stream",
            "--full-payload-solve", "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", "9", "--exact-packet-attempt", "11",
            "--construction-seed-basis", "production-profile",
        ])
        last = subject.Coordinate(1079)
        self.assertEqual(last.baseline_record_ordinal, 267299)
        self.assertEqual(subject.Coordinate(0).baseline_record_ordinal, 324)
        with self.assertRaises(subject.Stage1Error):
            subject.Invocation(1, 0, 0)

    def test_packet_seed_rule_handles_uint8_wrap(self) -> None:
        self.assertEqual(subject.stage0._packet_seed_for_offset(
            "0x12345678", 255, subject.PACKET_DELTA), "0x172990ea")

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
        self.assertEqual(record["mix_pair"], "01")
        self.assertEqual(record["delta"], 2)
        self.assertEqual(record["precode_attempt"], 9)
        self.assertEqual(record["packet_attempt"], 11)
        self.assertEqual(
            record["controller_source_git_commit"],
            CONTROLLER_SOURCE_COMMIT)
        self.assertEqual(
            record["bench_source_git_commit"], BENCH_SOURCE_COMMIT)
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
            stdout_for(invocation).replace(b"mix_pair=01", b"mix_pair=02"),
        ]
        tokens = valid[0].split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        cases.append((" ".join(tokens) + "\n" +
                      "\n".join(valid[1:]) + "\n").encode("ascii"))
        for stdout in cases:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.Stage1Error):
                subject.parse_process_result(
                    process_result(invocation, stdout),
                    CONTROLLER_SOURCE_COMMIT, BENCH_SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA, baseline_row(0))

    def test_configuration_failure_transport_is_exact_but_invalid(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=11 result=2\n").encode("ascii")
        record = subject.parse_process_result(
            process_result(invocation, stdout, 2, stderr),
            CONTROLLER_SOURCE_COMMIT, BENCH_SOURCE_COMMIT,
            BENCH_SHA, SOURCE_SHA, baseline_row(0))
        self.assertTrue(record["configuration_failure"])
        self.assertTrue(record["weak"])
        records = successful_roster(
            [baseline_row(index) for index in range(1080)])
        records[0] = record
        with self.assertRaisesRegex(subject.Stage1Error, "precode p"):
            subject.adjudicate(
                records, [baseline_row(index) for index in range(1080)])
        with self.assertRaises(subject.Stage1Error):
            subject.parse_process_result(
                process_result(invocation, stdout, 2, stderr + b"extra\n"),
                CONTROLLER_SOURCE_COMMIT, BENCH_SOURCE_COMMIT,
                BENCH_SHA, SOURCE_SHA, baseline_row(0))

    def test_rejects_noncanonical_work_and_wrong_baseline_coordinate(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        stdout = stdout_for(invocation).replace(
            b",100.000,50.000,", b",100.00,50.000,")
        with self.assertRaises(subject.Stage1Error):
            subject.parse_process_result(
                process_result(invocation, stdout),
                CONTROLLER_SOURCE_COMMIT, BENCH_SOURCE_COMMIT,
                BENCH_SHA, SOURCE_SHA, baseline_row(0))
        with self.assertRaises(subject.Stage1Error):
            parsed_record(invocation, baseline_row(1))

    def test_execute_roster_fails_before_second_call_on_config_or_seed(self) \
            -> None:
        invocation = subject.Invocation(0, 0, 0)
        baseline = [baseline_row(index) for index in range(1080)]
        config_stdout = (metadata_line(invocation) + "\n" +
                         ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        config_stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=11 result=2\n").encode("ascii")
        packet_seed = csv_values(invocation)[
            "effective_packet_seed"].encode("ascii")
        cases = (
            (process_result(
                invocation, config_stdout, 2, config_stderr), "precode p"),
            (process_result(
                invocation, stdout_for(invocation).replace(
                    b"0x123456789abcdef0", b"0x0000000000000001")),
             "precode seed"),
            (process_result(
                invocation, stdout_for(invocation).replace(
                    packet_seed, b"0x00000001")),
             "packet seed"),
        )
        for result, message in cases:
            with self.subTest(message=message), mock.patch.object(
                    subject, "_run_raw", return_value=result) as run_raw, \
                    self.assertRaisesRegex(subject.Stage1Error, message):
                subject._execute_roster(
                    mock.Mock(), CONTROLLER_SOURCE_COMMIT,
                    BENCH_SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA, baseline)
            run_raw.assert_called_once()

    def test_run_raw_uses_pinned_descriptor_and_frozen_environment(self) -> None:
        invocation = subject.Invocation(0, 0, 0)
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        completed = mock.Mock(
            stdout=b"stdout\n", stderr=b"", returncode=0)
        with mock.patch.object(
                subject.subprocess, "run", return_value=completed) as run:
            result = subject._run_raw(invocation, pinned)
        self.assertEqual(result.invocation, invocation)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(run.call_args.args[0], invocation.argv(pinned.path))
        self.assertEqual(
            run.call_args.kwargs["executable"], "/proc/self/fd/17")
        self.assertEqual(run.call_args.kwargs["pass_fds"], (17,))
        self.assertEqual(
            run.call_args.kwargs["env"], subject.CHILD_ENVIRONMENT)
        self.assertTrue(run.call_args.kwargs["start_new_session"])
        self.assertTrue(run.call_args.kwargs["close_fds"])


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

    def test_stage0_loader_reaudits_complete_selection_and_v9_link(self) -> None:
        baseline = [baseline_row(index) for index in range(1080)]
        precode_seeds = []
        packet_seeds = []
        controls = []
        for coordinate in subject.stage0.COORDINATES:
            ordinal = (
                subject.K_VALUES.index(coordinate.K) * 36 +
                coordinate.root_index * 3 +
                subject.SCHEDULES.index(coordinate.schedule))
            baseline[ordinal] = baseline_row(ordinal, success=False)
            row = baseline[ordinal]
            precode_seeds.append(row["effective_precode_seed"])
            packet_seeds.append(row["effective_packet_seed"])
            controls.append({
                "staircase": subject.STAIRCASE_BY_K[row["K"]],
                "source_hits": 3 if row["K"] >= 10000 else 2,
                "outcome": row["outcome"], "success": row["success"],
                "rank_fail": row["rank_fail"], "error": row["error"],
                "inactivated_columns": row["inactivated_columns"],
                "effective_precode_seed": row["effective_precode_seed"],
                "effective_packet_seed": row["effective_packet_seed"],
                "block_xors": row["block_xors"],
                "gf256_muladds": row["gf256_muladds"],
            })
        source = {
            "source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT,
            "tracked_and_untracked_linked_sources_clean": True,
            "all_source_paths_tracked_at_HEAD": True,
            "working_source_bytes_equal_HEAD_blobs": True,
            "clean_status_scope": [], "source_files": {},
        }
        source["source_receipt_sha256"] = subject.stage0.self_hash(
            source, "source_receipt_sha256")
        v9 = {
            "contract_sha256": subject.EXPECTED_V9_CONTRACT_SHA256,
            "source_git_commit": subject.EXPECTED_V9_SOURCE_COMMIT,
            "worker_sha256": subject.EXPECTED_V9_WORKER_SHA256,
            "files": {
                name: {"path": "/frozen/" + name, "size": spec[0],
                       "sha256": spec[1]}
                for name, spec in subject.V9_FILES.items()
            },
            "effective_precode_seeds": precode_seeds,
            "effective_packet_seeds": packet_seeds,
            "diagonal_control": controls,
        }
        v9["receipt_sha256"] = subject.stage0.self_hash(
            v9, "receipt_sha256")
        verdict = {
            "diagonal_control": {"disposition": "PASS"},
            "delta_results": [], "survivors": [2, 3],
            "stage1_candidate_delta": 2,
            "effective_precode_seeds": precode_seeds,
            "authenticated_effective_p_packet_seeds": packet_seeds,
            "disposition": "PASS",
        }
        records = [{
            "schema": subject.stage0.RECORD_SCHEMA, "ordinal": ordinal,
            "bench_binary_sha256": subject.EXPECTED_STAGE0_BENCH_SHA256,
            "source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT,
            "source_receipt_sha256": source["source_receipt_sha256"],
            "deterministic_projection_sha256": subject.stage0.sha256_json(
                {"ordinal": ordinal}),
        } for ordinal in range(2)]
        summary = {
            "schema": subject.stage0.SUMMARY_SCHEMA,
            "contract_sha256": subject.EXPECTED_STAGE0_CONTRACT_SHA256,
            "records_file": "mix2-packet-offset-stage0-records.jsonl",
            "records_sha256": subject.STAGE0_FILES[
                "mix2-packet-offset-stage0-records.jsonl"][1],
            "record_count": 2, **verdict,
            "stage0_only": True, "promotion_evidence": False,
            "fresh_roots_used": False, "timing_evidence": False,
            "source_receipt": source,
            "bench": {"sha256": subject.EXPECTED_STAGE0_BENCH_SHA256,
                      "size": subject.EXPECTED_STAGE0_BENCH_SIZE},
            "bench_description": {
                "schema": subject.BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": subject.EXPECTED_STAGE0_SOURCE_COMMIT},
            "v9_prerequisite": v9,
        }
        summary["summary_sha256"] = subject.self_hash(
            summary, "summary_sha256")
        summary_data = (subject.canonical_json(summary) + "\n").encode("ascii")
        record_data = b"".join(
            (subject.canonical_json(row) + "\n").encode("ascii")
            for row in records)
        pins = {
            "mix2-packet-offset-stage0-summary.json": MemoryPin(summary_data),
            "mix2-packet-offset-stage0-records.jsonl": MemoryPin(record_data),
        }
        v9_receipt = {
            "attempt_stream_sha256": subject.V9_FILES[
                "promotion-short-screen-attempts.jsonl"][1],
            "result_stream_sha256": subject.V9_FILES[
                "promotion-short-screen-results.jsonl"][1],
            "source_git_commit": subject.EXPECTED_V9_SOURCE_COMMIT,
        }
        with mock.patch.object(
                subject, "EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256",
                source["source_receipt_sha256"]), \
                mock.patch.object(
                    subject, "EXPECTED_STAGE0_V9_RECEIPT_SHA256",
                    v9["receipt_sha256"]), \
                mock.patch.object(
                    subject, "EXPECTED_STAGE0_SUMMARY_SELF_SHA256",
                    summary["summary_sha256"]), \
                mock.patch.object(
                    subject, "EXPECTED_STAGE0_SURVIVOR_COUNT", 2), \
                mock.patch.object(
                    subject.stage0, "EXPECTED_INVOCATION_COUNT", 2), \
                mock.patch.object(
                    subject.stage0, "deterministic_projection",
                    side_effect=lambda row: {"ordinal": row["ordinal"]}), \
                mock.patch.object(
                    subject.stage0, "adjudicate", return_value=verdict) as gate:
            receipt = subject._load_stage0(pins, v9_receipt, baseline)
        self.assertEqual(receipt["record_count"], 2)
        self.assertEqual(receipt["stage1_candidate_delta"], 2)
        gate.assert_called_once()


class SourceReceiptTests(unittest.TestCase):
    @staticmethod
    def _git_result(stdout: bytes = b"", stderr: bytes = b"",
                    returncode: int = 0):
        return mock.Mock(
            stdout=stdout, stderr=stderr, returncode=returncode)

    def _responses(self, second_blob: bytes = b"beta\n",
                   final_head: str = CONTROLLER_SOURCE_COMMIT):
        return [
            self._git_result(
                (CONTROLLER_SOURCE_COMMIT + "\n").encode("ascii")),
            self._git_result(),
            self._git_result(b"one\0two\0"),
            self._git_result(b"alpha\n"),
            self._git_result(second_blob),
            self._git_result((final_head + "\n").encode("ascii")),
        ]

    def test_source_receipt_binds_working_bytes_to_captured_head(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "one").write_bytes(b"alpha\n")
            (root / "two").write_bytes(b"beta\n")
            with mock.patch.object(subject, "REPO_ROOT", root), \
                    mock.patch.object(subject, "SOURCE_PATHS", ("one", "two")), \
                    mock.patch.object(
                        subject, "SOURCE_STATUS_PATHS", ("one", "two")), \
                    mock.patch.object(
                        subject.subprocess, "run",
                        side_effect=self._responses()) as run:
                receipt = subject._source_receipt()
        self.assertEqual(
            receipt["source_git_commit"], CONTROLLER_SOURCE_COMMIT)
        self.assertTrue(receipt["working_source_bytes_equal_HEAD_blobs"])
        self.assertEqual(
            run.call_args_list[3].args[0],
            ["git", "cat-file", "blob",
             CONTROLLER_SOURCE_COMMIT + ":one"])
        self.assertEqual(
            run.call_args_list[4].args[0],
            ["git", "cat-file", "blob",
             CONTROLLER_SOURCE_COMMIT + ":two"])
        self.assertEqual(
            run.call_args_list[5].args[0], ["git", "rev-parse", "HEAD"])

    def test_source_receipt_rejects_blob_or_final_head_drift(self) -> None:
        cases = (
            (self._responses(second_blob=b"wrong\n"), "byte-for-byte"),
            (self._responses(final_head="9" * 40), "HEAD changed"),
        )
        for responses, message in cases:
            with self.subTest(message=message), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "one").write_bytes(b"alpha\n")
                (root / "two").write_bytes(b"beta\n")
                with mock.patch.object(subject, "REPO_ROOT", root), \
                        mock.patch.object(
                            subject, "SOURCE_PATHS", ("one", "two")), \
                        mock.patch.object(
                            subject, "SOURCE_STATUS_PATHS", ("one", "two")), \
                        mock.patch.object(
                            subject.subprocess, "run", side_effect=responses), \
                        self.assertRaisesRegex(subject.Stage1Error, message):
                    subject._source_receipt()


class PublicationTests(unittest.TestCase):
    def test_rejects_non_stage0_benchmark_before_description(self) -> None:
        baseline = [baseline_row(index) for index in range(1080)]
        historical = FakeHistorical(baseline)
        historical.receipt = mock.Mock(
            return_value=historical.receipt_value)
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = {
            "path": "/pinned/bench", "sha256": "9" * 64,
            "size": subject.EXPECTED_STAGE0_BENCH_SIZE,
        }
        with mock.patch.object(
                subject.screen, "read_bench_description") as description, \
                mock.patch.object(subject, "_execute_roster") as execute, \
                self.assertRaisesRegex(subject.Stage1Error, "exact authenticated"):
            subject._run_stage1_pinned(
                pinned, Path("/unused/out"), source_receipt(), {}, historical,
                historical.receipt_value)
        description.assert_not_called()
        execute.assert_not_called()

    def test_publication_is_no_replace_and_parent_fsynced(self) -> None:
        events = []
        baseline = [baseline_row(index) for index in range(1080)]
        for index in (0, 500, 1079):
            baseline[index] = baseline_row(index, success=False)
        records = successful_roster(baseline)
        historical = FakeHistorical(baseline)
        historical.receipt = mock.Mock(
            return_value=historical.receipt_value)
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.side_effect = lambda: bench_receipt()

        def description(*_args):
            events.append("describe")
            return {"schema": subject.BENCH_DESCRIPTION_SCHEMA,
                    "source_git_commit": BENCH_SOURCE_COMMIT}

        def execute(*_args):
            events.append("work")
            return records

        def publish(source: Path, destination: Path):
            events.append("publish")
            summary = json.loads(
                (source / subject.SUMMARY_NAME).read_text(encoding="ascii"))
            self.assertEqual(summary["record_count"], 2160)
            self.assertEqual(summary["candidate_packet_delta"], 2)
            os.rename(source, destination)

        def fsync(path: Path):
            events.append(
                "fsync-staging" if path.name.startswith(
                    ".wh2-mix2-packet-offset-stage1-") else "fsync-parent")

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  side_effect=description), \
                mock.patch.object(subject, "_execute_roster",
                                  side_effect=execute), \
                mock.patch.object(
                    subject, "_source_receipt",
                    return_value=source_receipt()) as source_attest, \
                mock.patch.object(
                    subject.screen, "_publish_directory_noreplace",
                    side_effect=publish), \
                mock.patch.object(subject, "_fsync_directory",
                                  side_effect=fsync):
            output = Path(temporary) / "out"
            summary = subject._run_stage1_pinned(
                pinned, output, source_receipt(),
                subject.contract_description(), historical,
                historical.receipt_value)
            self.assertTrue(output.is_dir())
            self.assertEqual(summary["disposition"], "PASS")
        self.assertEqual(events, [
            "describe", "work", "fsync-staging", "publish", "fsync-parent"])
        self.assertEqual(pinned.receipt.call_count, 3)
        self.assertEqual(source_attest.call_count, 2)
        self.assertEqual(historical.receipt.call_count, 2)

    def test_publisher_does_not_replace_existing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            destination = root / "destination"
            source.mkdir()
            destination.mkdir()
            (source / "marker").write_bytes(b"source")
            (destination / "marker").write_bytes(b"destination")
            with self.assertRaises(subject.screen.ShortScreenError):
                subject.screen._publish_directory_noreplace(
                    source, destination)
            self.assertEqual((source / "marker").read_bytes(), b"source")
            self.assertEqual(
                (destination / "marker").read_bytes(), b"destination")


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
            verdict["aggregates"]["pair01_delta2_candidate"]["block_xors"],
            108000)
        self.assertEqual(
            verdict["aggregates"]["candidate_to_baseline_ratios"][
                "block_xors"], "1.000000000000")
        self.assertEqual(len(verdict["per_K_work"]), 30)
        self.assertTrue(all(
            item["coordinate_count"] == 36 and
            item["candidate_observation_count"] == 72
            for item in verdict["per_K_work"]))
        self.assertEqual(
            verdict["per_K_work"][0]["metrics"]["block_xors"], {
                "pair01_baseline": 3600,
                "pair01_delta2_candidate": 3600,
                "candidate_minus_baseline": 0,
            })

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

    def test_second_repetition_error_is_reported_not_only_rejected(self) -> None:
        records, baseline = self.fresh()
        row = records[3]
        row.update({
            "outcome": "error", "success": False, "rank_fail": False,
            "error": True, "weak": True, "binary_deficiency": 0,
            "gf256_rank_gain": 0, "heavy_shortfall": 0,
            "first_rank_fail": -1, "failure_trials": [0],
        })
        rehash(row)
        verdict = subject.adjudicate(records, baseline)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(verdict["candidate_errors"], 1)
        self.assertEqual(len(verdict["candidate_weak_coordinates"]), 1)
        self.assertEqual(len(verdict["introduced_weak_coordinates"]), 1)
        self.assertFalse(verdict["gates"]["candidate_errors_zero"])
        self.assertFalse(
            verdict["gates"]["candidate_weak_coordinates_zero"])

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
            verdict["aggregates"]["pair01_delta2_candidate"]["block_xors"],
            108001)

    def test_other_deduplicated_work_regressions_reject(self) -> None:
        for field, gate in (
                ("gf256_muladds",
                 "aggregate_gf256_muladds_not_above_pair01"),
                ("inactivated_columns",
                 "aggregate_inactivated_columns_not_above_pair01")):
            records, baseline = self.fresh()
            for repetition in (0, 1):
                records[repetition][field] += 1
                rehash(records[repetition])
            verdict = subject.adjudicate(records, baseline)
            with self.subTest(field=field):
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate])

    def test_record_to_baseline_receipt_mismatch_is_invalid(self) -> None:
        records, baseline = self.fresh()
        records[0]["pair01_baseline"]["block_xors"] += 1
        rehash(records[0])
        with self.assertRaises(subject.Stage1Error):
            subject.adjudicate(records, baseline)

    def test_candidate_seed_or_packet_identity_drift_is_invalid(self) -> None:
        records, baseline = self.fresh()
        records[0]["effective_precode_seed"] = "0x0000000000000001"
        rehash(records[0])
        with self.assertRaisesRegex(subject.Stage1Error, "precode seed"):
            subject.adjudicate(records, baseline)

        records, baseline = self.fresh()
        records[0]["effective_packet_seed"] = "0x00000001"
        rehash(records[0])
        with self.assertRaisesRegex(subject.Stage1Error, "packet seed"):
            subject.adjudicate(records, baseline)

        records, baseline = self.fresh()
        records[0]["packet_attempt"] = 9
        rehash(records[0])
        with self.assertRaises(subject.Stage1Error):
            subject.adjudicate(records, baseline)

        records, baseline = self.fresh()
        records[0]["K"] = 3
        rehash(records[0])
        with self.assertRaises(subject.Stage1Error):
            subject.adjudicate(records, baseline)

        records, baseline = self.fresh()
        baseline[0]["loss_root"] = "0x0000000000000000"
        with self.assertRaisesRegex(subject.Stage1Error, "baseline roster"):
            subject.adjudicate(records, baseline)

    def test_record_benchmark_or_source_pin_drift_is_invalid(self) -> None:
        cases = (
            ("bench_binary_sha256", "9" * 64),
            ("bench_source_git_commit", "9" * 40),
            ("controller_source_git_commit", "8" * 40),
            ("source_receipt_sha256", "7" * 64),
        )
        for field, value in cases:
            records, baseline = self.fresh()
            records[0][field] = value
            rehash(records[0])
            with self.subTest(field=field), \
                    self.assertRaises(subject.Stage1Error):
                subject.adjudicate(records, baseline)


if __name__ == "__main__":
    unittest.main()
