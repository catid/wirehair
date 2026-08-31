#!/usr/bin/env python3
"""Pure mocked tests for the sealed row-hashed MIX2 Stage 0."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Optional
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_hashed_pair_stage0 as subject


SOURCE_COMMIT = "1" * 40
BENCH_SHA = "2" * 64
SOURCE_SHA = "3" * 64
HISTORY_SHA = "4" * 64

BASELINE_WORK = (
    {"block_xors": 26, "gf256_muladds": 179,
     "inactivated_columns": 14,
     "effective_precode_seed": "0x7cae730f87b736db",
     "effective_packet_seed": "0x278c881b"},
    {"block_xors": 10486, "gf256_muladds": 3787,
     "inactivated_columns": 73,
     "effective_precode_seed": "0xc75d4118cf078e41",
     "effective_packet_seed": "0xe7a48172"},
    {"block_xors": 108130, "gf256_muladds": 4782,
     "inactivated_columns": 157,
     "effective_precode_seed": "0x905ecdf4c1ca0f3e",
     "effective_packet_seed": "0xe2b88058"},
)


def baseline_rows() -> tuple:
    rows = []
    for index, coordinate_ordinal in enumerate(
            subject.BASELINE_COORDINATE_ORDINALS):
        coordinate = subject.Coordinate(coordinate_ordinal)
        work = BASELINE_WORK[index]
        rows.append({
            "schema": "synthetic-baseline",
            "record_ordinal": coordinate.baseline_record_ordinal,
            "K": coordinate.K,
            "attempt": coordinate.attempt,
            "root_index": coordinate.root_index,
            "loss_root": coordinate.root,
            "schedule_index": coordinate.schedule_index,
            "schedule": coordinate.schedule,
            "cell_ordinal": coordinate.cell_ordinal,
            "outcome": "rank_fail",
            "success": False,
            "rank_fail": True,
            "error": False,
            "binary_deficiency": 13,
            "gf256_rank_gain": 12,
            "inactivated_columns": work["inactivated_columns"],
            "heavy_shortfall": 0,
            "effective_precode_seed": work["effective_precode_seed"],
            "effective_packet_seed": work["effective_packet_seed"],
            "block_xors": work["block_xors"],
            "gf256_muladds": work["gf256_muladds"],
        })
    return tuple(rows)


BASELINES = baseline_rows()


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
        "source_git_commit": SOURCE_COMMIT,
        "mix_pair": subject.CANDIDATE_PAIR,
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(field, values[field])
        for field in subject.METADATA_ORDER)


def csv_values(invocation: subject.Invocation,
               outcome: str = "success",
               block_xors: Optional[int] = None,
               gf256_muladds: Optional[int] = None,
               inactivated_columns: Optional[int] = None,
               effective_precode_seed: Optional[str] = None,
               effective_packet_seed: Optional[str] = None) -> dict:
    baseline = BASELINES[invocation.coordinate.selection_index]
    block_xors = (baseline["block_xors"] if block_xors is None
                  else block_xors)
    gf256_muladds = (baseline["gf256_muladds"] if gf256_muladds is None
                     else gf256_muladds)
    inactivated_columns = (
        baseline["inactivated_columns"] if inactivated_columns is None
        else inactivated_columns)
    effective_precode_seed = (
        baseline["effective_precode_seed"]
        if effective_precode_seed is None else effective_precode_seed)
    effective_packet_seed = (
        baseline["effective_packet_seed"]
        if effective_packet_seed is None else effective_packet_seed)
    if outcome == "success":
        success, rank_fail, error, deficiency, gain = 1, 0, 0, 10, 10
    elif outcome == "rank_fail":
        success, rank_fail, error, deficiency, gain = 0, 1, 0, 13, 12
    elif outcome == "error":
        success, rank_fail, error, deficiency, gain = 0, 0, 1, 0, 0
    else:
        raise AssertionError("unknown synthetic outcome")
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
        "inact_mu": "{}.000".format(inactivated_columns),
        "inact_max": str(inactivated_columns),
        "binary_def_mu": "{}.000".format(deficiency),
        "binary_def_max": str(deficiency),
        "heavy_gain_mu": "{}.000".format(gain),
        "heavy_gain_min": str(gain),
        "heavy_shortfall": "0",
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
        "effective_precode_seed": effective_precode_seed,
        "effective_packet_seed": effective_packet_seed,
    }


def stdout_for(invocation: subject.Invocation, **kwargs: object) -> bytes:
    values = csv_values(invocation, **kwargs)
    return ("\n".join((
        metadata_line(invocation),
        ",".join(subject.CSV_HEADER),
        ",".join(values[field] for field in subject.CSV_HEADER),
    )) + "\n").encode("ascii")


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


def parsed_record(invocation: subject.Invocation, **kwargs: object) -> dict:
    baseline = BASELINES[invocation.coordinate.selection_index]
    stdout = stdout_for(invocation, **kwargs)
    return dict(subject.parse_process_result(
        process_result(invocation, stdout), SOURCE_COMMIT, BENCH_SHA,
        SOURCE_SHA, HISTORY_SHA, baseline))


def successful_roster() -> list:
    return [parsed_record(value) for value in subject.make_invocations()]


def source_receipt() -> dict:
    return {
        "source_git_commit": SOURCE_COMMIT,
        "source_receipt_sha256": SOURCE_SHA,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "clean_status_scope": [], "source_files": {},
    }


def bench_receipt() -> dict:
    return {"path": "/pinned/bench", "sha256": BENCH_SHA, "size": 1}


class FakeHistorical:
    def __init__(self, receipt: Optional[dict] = None) -> None:
        self._receipt = receipt or {"authenticated": True}
        self.baseline = [None] * subject.pair_stage1.COORDINATE_COUNT
        for ordinal, row in zip(subject.BASELINE_COORDINATE_ORDINALS,
                                BASELINES):
            self.baseline[ordinal] = row

    def receipt(self) -> dict:
        return dict(self._receipt)

    def __enter__(self) -> "FakeHistorical":
        return self

    def __exit__(self, *_args: object) -> None:
        return None


class ContractAndInvocationTests(unittest.TestCase):
    def test_frozen_contract_and_exact_coordinate_major_roster(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.contract_description()["contract_sha256"],
                         subject.EXPECTED_CONTRACT_SHA256)
        invocations = subject.make_invocations()
        self.assertEqual(len(invocations), 6)
        self.assertEqual(
            [(value.coordinate_ordinal, value.repetition)
             for value in invocations],
            [(coordinate, repetition)
             for coordinate in (30, 573, 820)
             for repetition in (0, 1)])
        self.assertEqual(
            [subject.Coordinate(value).baseline_record_ordinal
             for value in subject.BASELINE_COORDINATE_ORDINALS],
            list(subject.BASELINE_MATRIX_ORDINALS))

    def test_contract_and_imported_helper_drift_fail_closed(self) -> None:
        with mock.patch.object(subject, "EXPECTED_CONTRACT_SHA256", "0" * 64), \
                self.assertRaises(subject.HashedPairStage0Error):
            subject.contract_description()
        with mock.patch.object(
                subject.pair_stage1, "_validate_constants",
                side_effect=subject.pair_stage1.Stage1Error("drift")), \
                self.assertRaisesRegex(subject.HashedPairStage0Error,
                                       "helper changed"):
            subject._validate_constants()

    def test_argv_is_exact_hashed_b2_loss_half_oh0_cell(self) -> None:
        invocation = subject.Invocation(0, 30, 0)
        argv = invocation.argv(Path("/pinned/bench"))
        self.assertEqual(argv[:4],
                         ["/pinned/bench", "precodefail", "--N", "2"])
        self.assertEqual(argv[argv.index("--mix-pair") + 1], "hashed")
        self.assertEqual(argv[argv.index("--exact-precode-attempt") + 1], "9")
        self.assertEqual(argv[argv.index("--exact-packet-attempt") + 1], "9")
        for required in ("--paired-overhead-stream", "--full-payload-solve"):
            self.assertIn(required, argv)
        with self.assertRaises(subject.HashedPairStage0Error):
            subject.Invocation(1, 30, 0)

    def test_run_parser_requires_all_five_paths(self) -> None:
        parser = subject._parser()
        with mock.patch("sys.stderr"), \
                self.assertRaises(SystemExit) as raised:
                parser.parse_args(["run", "--bench", "/b"])
        self.assertEqual(raised.exception.code, 2)


class ParserTests(unittest.TestCase):
    def test_exact_success_binds_source_binary_argv_raw_and_baseline(self) \
            -> None:
        invocation = subject.Invocation(4, 820, 0)
        record = parsed_record(invocation)
        self.assertEqual(record["schema"], subject.RECORD_SCHEMA)
        self.assertEqual(record["pair"], "hashed")
        self.assertEqual(record["coordinate_ordinal"], 820)
        self.assertEqual(record["bench_binary_sha256"], BENCH_SHA)
        self.assertEqual(record["source_receipt_sha256"], SOURCE_SHA)
        self.assertEqual(record["historical_evidence_sha256"], HISTORY_SHA)
        self.assertEqual(record["argv_sha256"], record["command_sha256"])
        self.assertEqual(record["raw_stdout_sha256"], record["stdout_sha256"])
        self.assertEqual(record["raw_stderr_sha256"], record["stderr_sha256"])
        self.assertEqual(record["pair01_baseline"],
                         subject._baseline_projection(BASELINES[2]))
        self.assertEqual(record["deterministic_projection_sha256"],
                         subject.sha256_json(
                             subject.deterministic_projection(record)))

    def test_metadata_drift_and_malformed_framing_are_rejected(self) -> None:
        invocation = subject.Invocation(0, 30, 0)
        valid = stdout_for(invocation)
        cases = (
            valid.replace(b"mix_pair=hashed", b"mix_pair=01"),
            valid[:-1],
            valid.replace(b"\n", b"\r\n", 1),
            valid + b"extra\n",
            valid.replace(b"block_xors_mu", b"block_xors_mx"),
        )
        for stdout in cases:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.HashedPairStage0Error):
                subject.parse_process_result(
                    process_result(invocation, stdout), SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA, HISTORY_SHA, BASELINES[0])

    def test_canonical_configuration_failure_is_scientific_reject(self) \
            -> None:
        invocation = subject.Invocation(0, 30, 0)
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=9 result=2\n").encode("ascii")
        record = subject.parse_process_result(
            process_result(invocation, stdout, 2, stderr), SOURCE_COMMIT,
            BENCH_SHA, SOURCE_SHA, HISTORY_SHA, BASELINES[0])
        self.assertTrue(record["configuration_failure"])
        self.assertEqual(record["outcome"], "configuration_failure")

    def test_wrong_baseline_binding_is_rejected(self) -> None:
        invocation = subject.Invocation(0, 30, 0)
        with self.assertRaisesRegex(subject.HashedPairStage0Error,
                                    "baseline changed"):
            subject.parse_process_result(
                process_result(invocation, stdout_for(invocation)),
                SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA, HISTORY_SHA,
                BASELINES[1])

    def test_baseline_must_be_an_unambiguous_rank_failure(self) -> None:
        invocation = subject.Invocation(0, 30, 0)
        for mutation in (
                {"rank_fail": False}, {"error": True},
                {"outcome": "error"}):
            baseline = dict(BASELINES[0], **mutation)
            with self.subTest(mutation=mutation), \
                    self.assertRaisesRegex(subject.HashedPairStage0Error,
                                           "baseline changed"):
                subject.parse_process_result(
                    process_result(invocation, stdout_for(invocation)),
                    SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA, HISTORY_SHA,
                    baseline)


class GateTests(unittest.TestCase):
    def test_exact_six_successes_pass_and_report_each_work_ledger(self) -> None:
        verdict = subject.adjudicate(successful_roster(), BASELINES)
        self.assertEqual(verdict["disposition"], "PASS")
        self.assertTrue(all(verdict["gates"].values()))
        self.assertEqual(verdict["repeat_exact_coordinate_count"], 3)
        for coordinate in verdict["coordinates"]:
            self.assertTrue(all(
                value["not_above_baseline"]
                for value in coordinate["candidate_vs_pair01"].values()))

    def test_failure_and_repeat_mismatch_reject(self) -> None:
        invocations = subject.make_invocations()
        records = successful_roster()
        records[0] = parsed_record(invocations[0], outcome="rank_fail")
        records[3] = parsed_record(
            invocations[3], block_xors=BASELINES[1]["block_xors"] - 1)
        verdict = subject.adjudicate(records, BASELINES)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertFalse(verdict["gates"]["all_six_observations_succeed"])
        self.assertFalse(verdict["gates"][
            "all_three_repeated_deterministic_projections_exact"])

    def test_each_per_coordinate_work_regression_rejects(self) -> None:
        invocations = subject.make_invocations()
        keyword = {
            "block_xors": "block_xors",
            "gf256_muladds": "gf256_muladds",
            "inactivated_columns": "inactivated_columns",
        }
        gate = {
            "block_xors": "per_coordinate_block_xors_not_above_pair01",
            "gf256_muladds":
                "per_coordinate_gf256_muladds_not_above_pair01",
            "inactivated_columns":
                "per_coordinate_inactivated_columns_not_above_pair01",
        }
        for field in keyword:
            records = successful_roster()
            value = BASELINES[0][field] + 1
            for index in (0, 1):
                records[index] = parsed_record(
                    invocations[index], **{keyword[field]: value})
            verdict = subject.adjudicate(records, BASELINES)
            with self.subTest(field=field):
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate[field]])

    def test_inactivation_cap_is_independent_hard_gate(self) -> None:
        invocations = subject.make_invocations()
        records = successful_roster()
        for index in (4, 5):
            records[index] = parsed_record(
                invocations[index], inactivated_columns=1025)
        verdict = subject.adjudicate(records, BASELINES)
        self.assertFalse(verdict["gates"][
            "all_six_inactivated_columns_at_most_1024"])
        self.assertEqual(verdict["disposition"], "REJECT")

    def test_effective_seed_mismatch_rejects_and_is_reported(self) -> None:
        invocations = subject.make_invocations()
        cases = (
            ("effective_precode_seed", "0x1111111111111111",
             "per_coordinate_effective_precode_seed_equals_pair01"),
            ("effective_packet_seed", "0x11111111",
             "per_coordinate_effective_packet_seed_equals_pair01"),
        )
        for field, value, gate in cases:
            records = successful_roster()
            for index in (0, 1):
                records[index] = parsed_record(
                    invocations[index], **{field: value})
            verdict = subject.adjudicate(records, BASELINES)
            with self.subTest(field=field):
                self.assertEqual(verdict["disposition"], "REJECT")
                self.assertFalse(verdict["gates"][gate])
                self.assertFalse(verdict["coordinates"][0][
                    "effective_seed_identity"][field]["equal"])

    def test_cardinality_and_record_identity_fail_closed(self) -> None:
        records = successful_roster()
        with self.assertRaises(subject.HashedPairStage0Error):
            subject.adjudicate(records[:-1], BASELINES)
        records[0] = dict(records[0], ordinal=1)
        records[0]["deterministic_projection_sha256"] = subject.sha256_json(
            subject.deterministic_projection(records[0]))
        with self.assertRaises(subject.HashedPairStage0Error):
            subject.adjudicate(records, BASELINES)

    def test_terminal_state_outcome_and_weak_must_be_consistent(self) -> None:
        mutations = (
            {"rank_fail": True},
            {"outcome": "rank_fail"},
            {"weak": True},
            {"configuration_failure": True},
        )
        for mutation in mutations:
            records = successful_roster()
            records[0] = dict(records[0], **mutation)
            records[0]["deterministic_projection_sha256"] = \
                subject.sha256_json(
                    subject.deterministic_projection(records[0]))
            with self.subTest(mutation=mutation), \
                    self.assertRaisesRegex(
                        subject.HashedPairStage0Error,
                        "record roster or receipt changed"):
                subject.adjudicate(records, BASELINES)


class StartupAndPublicationTests(unittest.TestCase):
    def test_exact_history_paths_are_required(self) -> None:
        with self.assertRaises(subject.HashedPairStage0Error):
            subject._validate_historical_paths(
                Path("/tmp/baseline"), subject.EXPECTED_V9_DIR,
                subject.EXPECTED_STAGE0_DIR)

    def test_history_authentication_failure_precedes_binary_and_work(self) \
            -> None:
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(
                    subject.pair_stage1, "load_historical_evidence",
                    side_effect=subject.pair_stage1.Stage1Error("bad history")) \
                as load, \
                mock.patch.object(subject.screen, "_open_binary") as open_bin, \
                mock.patch.object(subject, "_execute_roster") as execute, \
                self.assertRaisesRegex(subject.HashedPairStage0Error,
                                       "historical authentication"):
            subject.run_stage0(
                Path("/bench"), subject.EXPECTED_BASELINE_DIR,
                subject.EXPECTED_V9_DIR, subject.EXPECTED_STAGE0_DIR,
                Path(temporary) / "out")
        load.assert_called_once()
        open_bin.assert_not_called()
        execute.assert_not_called()

    def test_describe_is_last_preflight_before_candidate_work(self) -> None:
        events = []
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(
                    subject.screen, "read_bench_description",
                    side_effect=lambda *_: events.append("describe") or {
                        "schema": subject.BENCH_DESCRIPTION_SCHEMA,
                        "source_git_commit": SOURCE_COMMIT}), \
                mock.patch.object(subject, "_execute_roster",
                                  side_effect=lambda *_:
                                  events.append("work") or
                                  successful_roster()), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(subject.screen,
                                  "_publish_directory_noreplace",
                                  side_effect=lambda source, destination:
                                  os.rename(source, destination)):
            summary = subject._run_stage0_pinned(
                pinned, Path(temporary) / "out", source_receipt(),
                subject.contract_description(), history, receipt,
                subject.sha256_json(receipt), BASELINES)
        self.assertEqual(events, ["describe", "work"])
        self.assertEqual(summary["disposition"], "PASS")

    def test_atomic_publication_contains_canonical_records_and_summary(self) \
            -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()
        publications = []

        def publish(source: Path, destination: Path) -> None:
            publications.append((source, destination))
            self.assertTrue((source / subject.RECORD_NAME).is_file())
            self.assertTrue((source / subject.SUMMARY_NAME).is_file())
            summary = json.loads(
                (source / subject.SUMMARY_NAME).read_text(encoding="ascii"))
            self.assertEqual(summary["record_count"], 6)
            self.assertEqual(summary["summary_sha256"],
                             subject.self_hash(summary, "summary_sha256"))
            os.rename(source, destination)

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  return_value={
                                      "schema":
                                          subject.BENCH_DESCRIPTION_SCHEMA,
                                      "source_git_commit": SOURCE_COMMIT}), \
                mock.patch.object(subject, "_execute_roster",
                                  return_value=successful_roster()), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(subject.screen,
                                  "_publish_directory_noreplace",
                                  side_effect=publish):
            output = Path(temporary) / "out"
            summary = subject._run_stage0_pinned(
                pinned, output, source_receipt(),
                subject.contract_description(), history, receipt,
                subject.sha256_json(receipt), BASELINES)
            self.assertTrue(output.is_dir())
            self.assertEqual(summary["records_file"], subject.RECORD_NAME)
        self.assertEqual(len(publications), 1)

    def test_source_or_history_drift_prevents_publication(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()
        drifted_source = dict(source_receipt(), source_receipt_sha256="5" * 64)
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  return_value={}), \
                mock.patch.object(subject, "_execute_roster",
                                  return_value=successful_roster()), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=drifted_source), \
                mock.patch.object(subject.screen,
                                  "_publish_directory_noreplace") as publish, \
                self.assertRaisesRegex(subject.HashedPairStage0Error,
                                       "changed during run"):
            subject._run_stage0_pinned(
                pinned, Path(temporary) / "out", source_receipt(),
                subject.contract_description(), history, receipt,
                subject.sha256_json(receipt), BASELINES)
        publish.assert_not_called()

    def test_history_receipt_drift_prevents_publication(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        original_receipt = history.receipt()
        history._receipt = {"authenticated": False}
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  return_value={}), \
                mock.patch.object(subject, "_execute_roster",
                                  return_value=successful_roster()), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(subject.screen,
                                  "_publish_directory_noreplace") as publish, \
                self.assertRaisesRegex(subject.HashedPairStage0Error,
                                       "changed during run"):
            subject._run_stage0_pinned(
                pinned, Path(temporary) / "out", source_receipt(),
                subject.contract_description(), history, original_receipt,
                subject.sha256_json(original_receipt), BASELINES)
        publish.assert_not_called()

    def test_raw_record_and_summary_write_errors_fail_closed(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()
        cases = (
            (subject.pair_stage0, "_write_jsonl", "stage hashed-pair records"),
            (subject.screen, "_write_exclusive", "stage hashed-pair summary"),
        )
        for owner, name, message in cases:
            with self.subTest(helper=name), \
                    tempfile.TemporaryDirectory() as temporary, \
                    mock.patch.object(
                        subject.screen, "read_bench_description",
                        return_value={}), \
                    mock.patch.object(subject, "_execute_roster",
                                      return_value=successful_roster()), \
                    mock.patch.object(subject, "_source_receipt",
                                      return_value=source_receipt()), \
                    mock.patch.object(owner, name,
                                      side_effect=OSError("disk failure")), \
                    mock.patch.object(
                        subject.screen, "_publish_directory_noreplace") \
                    as publish, \
                    self.assertRaisesRegex(
                        subject.HashedPairStage0Error, message):
                subject._run_stage0_pinned(
                    pinned, Path(temporary) / "out", source_receipt(),
                    subject.contract_description(), history, receipt,
                    subject.sha256_json(receipt), BASELINES)
            publish.assert_not_called()

    def test_parent_sync_error_is_controlled_after_publication(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()

        def publish(source: Path, destination: Path) -> None:
            os.rename(source, destination)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "out"
            with mock.patch.object(
                        subject.screen, "read_bench_description",
                        return_value={}), \
                    mock.patch.object(subject, "_execute_roster",
                                      return_value=successful_roster()), \
                    mock.patch.object(subject, "_source_receipt",
                                      return_value=source_receipt()), \
                    mock.patch.object(
                        subject.screen, "_publish_directory_noreplace",
                        side_effect=publish), \
                    mock.patch.object(
                        subject.pair_stage0, "_fsync_directory",
                        side_effect=(None, OSError("parent fsync failed"))), \
                    self.assertRaisesRegex(
                        subject.HashedPairStage0Error,
                        "was published.*treat the output as invalid"):
                subject._run_stage0_pinned(
                    pinned, output, source_receipt(),
                    subject.contract_description(), history, receipt,
                    subject.sha256_json(receipt), BASELINES)
            self.assertTrue(output.is_dir())

    def test_cleanup_error_does_not_mask_original_failure(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        history = FakeHistorical()
        receipt = history.receipt()
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  return_value={}), \
                mock.patch.object(subject, "_execute_roster",
                                  return_value=successful_roster()), \
                mock.patch.object(
                    subject.pair_stage0, "_write_jsonl",
                    side_effect=OSError("original scientific failure")), \
                mock.patch.object(subject.shutil, "rmtree",
                                  side_effect=OSError("cleanup failed")) \
                as cleanup, \
                self.assertRaisesRegex(
                    subject.HashedPairStage0Error,
                    "original scientific failure") as raised:
            subject._run_stage0_pinned(
                pinned, Path(temporary) / "out", source_receipt(),
                subject.contract_description(), history, receipt,
                subject.sha256_json(receipt), BASELINES)
        cleanup.assert_called_once()
        if hasattr(raised.exception, "add_note"):
            self.assertTrue(any("cleanup failed" in note
                                for note in getattr(raised.exception,
                                                    "__notes__", ())))


if __name__ == "__main__":
    unittest.main()
