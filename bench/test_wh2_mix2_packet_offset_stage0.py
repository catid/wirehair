#!/usr/bin/env python3
"""Pure mocked tests for the sealed MIX2 packet-offset Stage 0."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_packet_offset_stage0 as subject


SOURCE_COMMIT = "1" * 40
BENCH_SHA = "2" * 64
SOURCE_SHA = "3" * 64
SYNTHETIC_PRECODE_SEEDS = {
    0: "0x123456789abcdef0",
    1: "0x123456789abcdef0",
    2: "0x123456789abcdef0",
}
SYNTHETIC_PACKET_SEEDS = {
    0: "0x00000000",
    1: "0x12345678",
    2: "0xffffffff",
}
SYNTHETIC_DIAGONAL_CONTROL = {
    coordinate: {
        "staircase": subject.STAIRCASE_BY_K[value.K],
        "source_hits": 3 if value.K >= 10000 else 2,
        "outcome": "rank_fail",
        "success": False,
        "rank_fail": True,
        "error": False,
        "inactivated_columns": 20,
        "effective_precode_seed": SYNTHETIC_PRECODE_SEEDS[coordinate],
        "effective_packet_seed": SYNTHETIC_PACKET_SEEDS[coordinate],
        "block_xors": 100,
        "gf256_muladds": 50,
    }
    for coordinate, value in enumerate(subject.COORDINATES)
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
        "exact_packet_attempt": str(invocation.packet_attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": SOURCE_COMMIT,
        "mix_pair": subject.MIX_PAIR,
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(field, values[field])
        for field in subject.METADATA_ORDER)


def csv_values(invocation: subject.Invocation,
               outcome: str = "success",
               block_xors: int = 100,
               gf256_muladds: int = 50) -> dict:
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
        "inact_mu": "20.000", "inact_max": "20",
        "binary_def_mu": "{}.000".format(deficiency),
        "binary_def_max": str(deficiency),
        "heavy_gain_mu": "{}.000".format(gain),
        "heavy_gain_min": str(gain),
        "heavy_shortfall": (
            "1" if rank_fail and deficiency <= 12 and gain < deficiency
            else "0"),
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
        "packet_attempt": str(invocation.packet_attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed":
            SYNTHETIC_PRECODE_SEEDS[invocation.coordinate_ordinal],
        "effective_packet_seed": "0x{:08x}".format(
            (int(SYNTHETIC_PACKET_SEEDS[
                invocation.coordinate_ordinal], 16) +
             (invocation.packet_attempt - coordinate.attempt) *
             subject.PACKET_ATTEMPT_STRIDE) & 0xffffffff),
    }


def stdout_for(invocation: subject.Invocation,
               outcome: str = "success",
               block_xors: int = 100,
               gf256_muladds: int = 50) -> bytes:
    values = csv_values(invocation, outcome, block_xors, gf256_muladds)
    lines = [metadata_line(invocation), ",".join(subject.CSV_HEADER),
             ",".join(values[field] for field in subject.CSV_HEADER)]
    return ("\n".join(lines) + "\n").encode("ascii")


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
                  outcome: str = "success",
                  block_xors: int = 100,
                  gf256_muladds: int = 50) -> dict:
    return dict(subject.parse_process_result(
        process_result(
            invocation,
            stdout_for(invocation, outcome, block_xors, gf256_muladds)),
        SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA))


def configuration_failure_record(invocation: subject.Invocation) -> dict:
    stdout = (metadata_line(invocation) + "\n" +
              ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
    coordinate = invocation.coordinate
    stderr = (
        "precodefail configuration failed N={} bb=2 "
        "heavy_family=periodic mix_count=2 attempt_mode=exact "
        "precode_attempt={} packet_attempt={} result=2\n".format(
            coordinate.K, coordinate.attempt,
            invocation.packet_attempt)).encode("ascii")
    return dict(subject.parse_process_result(
        process_result(invocation, stdout, 2, stderr),
        SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA))


def passing_roster() -> list:
    return [parsed_record(
                invocation,
                "rank_fail" if invocation.delta == 0 else "success")
            for invocation in subject.make_invocations()]


def roster_with_survivors(survivors) -> list:
    survivors = set(survivors)
    return [
        parsed_record(
            invocation,
            "success" if invocation.delta in survivors else "rank_fail")
        for invocation in subject.make_invocations()
    ]


def source_receipt() -> dict:
    return {
        "source_git_commit": SOURCE_COMMIT,
        "source_receipt_sha256": SOURCE_SHA,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "working_source_bytes_equal_HEAD_blobs": True,
        "clean_status_scope": [], "source_files": {},
    }


def bench_receipt() -> dict:
    return {"path": "/pinned/bench", "sha256": BENCH_SHA, "size": 1}


class FakeV9:
    def __init__(self) -> None:
        self.expected_precode_seeds = dict(SYNTHETIC_PRECODE_SEEDS)
        self.expected_packet_seeds = dict(SYNTHETIC_PACKET_SEEDS)
        self.expected_diagonal_control = dict(SYNTHETIC_DIAGONAL_CONTROL)
        self.receipt_value = {"files": {}, "receipt_sha256": "4" * 64}

    def receipt(self) -> dict:
        return self.receipt_value


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_exact_cardinality(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.contract_description()["contract_sha256"],
                         subject.EXPECTED_CONTRACT_SHA256)
        invocations = subject.make_invocations()
        self.assertEqual(len(invocations), 1536)
        self.assertEqual(
            (invocations[0].delta, invocations[0].packet_attempt,
             invocations[0].coordinate_ordinal, invocations[0].repetition),
            (0, 9, 0, 0))
        self.assertEqual(
            (invocations[-1].delta, invocations[-1].packet_attempt,
             invocations[-1].coordinate_ordinal,
             invocations[-1].repetition),
            (255, 255, 2, 1))
        self.assertEqual(len({
            (item.delta, item.coordinate_ordinal, item.repetition)
            for item in invocations}), 1536)
        self.assertFalse(
            set(subject.CONSUMED_ROOTS) &
            set(subject.UNCONSUMED_FINAL_ROOTS))

    def test_candidate_packet_attempts_cover_every_off_diagonal_q(self) -> None:
        invocations = subject.make_invocations()
        for coordinate_ordinal, coordinate in enumerate(subject.COORDINATES):
            packet_attempts = {
                item.packet_attempt
                for item in invocations
                if item.delta != 0 and
                item.coordinate_ordinal == coordinate_ordinal
            }
            self.assertEqual(
                packet_attempts, set(range(256)) - {coordinate.attempt})
        wrapped = next(
            item for item in invocations
            if item.delta == 247 and item.coordinate_ordinal == 0 and
            item.repetition == 0)
        self.assertEqual(wrapped.packet_attempt, 0)
        argv = wrapped.argv(Path("/pinned/bench"))
        self.assertEqual(argv[argv.index("--exact-packet-attempt") + 1], "0")

    def test_argv_is_the_frozen_off_diagonal_cell(self) -> None:
        invocation = subject.Invocation(6, 1, 0, 0)
        self.assertEqual(invocation.argv(Path("/pinned/bench")), [
            "/pinned/bench", "precodefail", "--N", "2", "--bb-list",
            "2", "--overhead", "0", "--trials", "1", "--threads", "1",
            "--loss", "0.5", "--seed", "0x7ccd510f122fc160",
            "--schedule", "burst", "--heavy-family", "periodic",
            "--mix-count", "2", "--mix-pair", "01",
            "--binary-dense-rows", "12", "--gf256-heavy-rows", "12",
            "--dense-anchors", "two07", "--paired-overhead-stream",
            "--full-payload-solve", "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", "9", "--exact-packet-attempt", "10",
            "--construction-seed-basis", "production-profile",
        ])
        with self.assertRaises(subject.Stage0Error):
            subject.Invocation(1, 1, 0, 0)

    def test_packet_seed_rule_handles_attempt_wrap_exactly(self) -> None:
        base = "0x278c881b"
        self.assertEqual(subject._packet_seed_for_offset(base, 9, 0), base)
        self.assertEqual(subject._packet_seed_for_offset(base, 9, 246),
                         "0x30db7fe1")
        self.assertEqual(subject._packet_seed_for_offset(base, 9, 247),
                         "0x9799409a")
        self.assertEqual(subject._packet_seed_for_offset(base, 9, 255),
                         "0x89550e62")
        self.assertEqual(subject._packet_seed_for_offset(
            "0xe7a48172", 0, 255), "0x80e6c0b9")


class ParserTests(unittest.TestCase):
    def test_accepts_exact_off_diagonal_metadata_and_success(self) -> None:
        invocation = subject.Invocation(12, 2, 0, 0)
        record = parsed_record(invocation)
        self.assertEqual(record["mix_pair"], "01")
        self.assertEqual(record["delta"], 2)
        self.assertEqual(record["precode_attempt"], 9)
        self.assertEqual(record["packet_attempt"], 11)
        self.assertEqual(record["outcome"], "success")
        self.assertEqual(
            record["deterministic_projection_sha256"],
            subject.sha256_json(subject.deterministic_projection(record)))

    def test_rejects_wrong_packet_attempt_or_pair_metadata(self) -> None:
        invocation = subject.Invocation(6, 1, 0, 0)
        valid = stdout_for(invocation)
        mutations = (
            valid.replace(b"exact_packet_attempt=10",
                          b"exact_packet_attempt=9"),
            valid.replace(b"mix_pair=01", b"mix_pair=12"),
            valid.replace(b",9,10,exact,", b",9,9,exact,"),
        )
        for stdout in mutations:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.Stage0Error):
                subject.parse_process_result(
                    process_result(invocation, stdout), SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA)

    def test_rejects_noncanonical_work_and_inconsistent_outcome(self) -> None:
        invocation = subject.Invocation(6, 1, 0, 0)
        valid = stdout_for(invocation)
        mutations = (
            valid.replace(b",100.000,50.000,", b",100.00,50.000,"),
            valid.replace(b",1,0,0,0.00000000,", b",1,1,0,0.00000000,"),
            valid.replace(b",10.000,10,10.000,10,0,",
                          b",10.000,10,9.000,9,0,"),
        )
        for stdout in mutations:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.Stage0Error):
                subject.parse_process_result(
                    process_result(invocation, stdout), SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA)

    def test_canonical_configuration_failure_records_q(self) -> None:
        invocation = subject.Invocation(6, 1, 0, 0)
        record = configuration_failure_record(invocation)
        self.assertTrue(record["configuration_failure"])
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=10 result=2\n").encode("ascii")
        with self.assertRaises(subject.Stage0Error):
            subject.parse_process_result(
                process_result(invocation, stdout, 2,
                               stderr.replace(b"packet_attempt=10",
                                              b"packet_attempt=9")),
                SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA)


class GateTests(unittest.TestCase):
    def test_all_deltas_survive_and_lowest_is_selected(self) -> None:
        result = subject.adjudicate(
            passing_roster(), SYNTHETIC_PRECODE_SEEDS,
            SYNTHETIC_PACKET_SEEDS, SYNTHETIC_DIAGONAL_CONTROL)
        self.assertEqual(result["disposition"], "PASS")
        self.assertEqual(result["survivors"], list(range(1, 256)))
        self.assertEqual(result["stage1_candidate_delta"], 1)
        self.assertTrue(all(value["disposition"] == "PASS"
                            for value in result["delta_results"]))
        self.assertEqual(result["diagonal_control"]["rank_failures"], 6)

    def test_current_diagonal_must_reproduce_authenticated_v9(self) -> None:
        records = passing_roster()
        records[0] = parsed_record(subject.make_invocations()[0], "success")
        with self.assertRaisesRegex(subject.Stage0Error,
                                    "diagonal control"):
            subject.adjudicate(
                records, SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
                SYNTHETIC_DIAGONAL_CONTROL)

    def test_lowest_failing_delta_falls_through_once(self) -> None:
        invocations = subject.make_invocations()
        records = passing_roster()
        for index in (6, 7):
            records[index] = parsed_record(invocations[index], "rank_fail")
        result = subject.adjudicate(
            records, SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
            SYNTHETIC_DIAGONAL_CONTROL)
        self.assertEqual(result["survivors"], list(range(2, 256)))
        self.assertEqual(result["stage1_candidate_delta"], 2)
        self.assertEqual(result["delta_results"][0]["rank_failures"], 2)

    def test_all_candidate_deltas_reject_without_selection(self) -> None:
        records = roster_with_survivors(set())
        result = subject.adjudicate(
            records, SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
            SYNTHETIC_DIAGONAL_CONTROL)
        self.assertEqual(result["disposition"], "REJECT")
        self.assertEqual(result["survivors"], [])
        self.assertIsNone(result["stage1_candidate_delta"])

    def test_noncontiguous_and_delta255_survivors_are_canonical(self) -> None:
        for survivors, selected in (({2, 255}, 2), ({255}, 255)):
            with self.subTest(survivors=survivors):
                result = subject.adjudicate(
                    roster_with_survivors(survivors),
                    SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
                    SYNTHETIC_DIAGONAL_CONTROL)
                self.assertEqual(result["survivors"], sorted(survivors))
                self.assertEqual(result["stage1_candidate_delta"], selected)
                self.assertEqual(result["disposition"], "PASS")

    def test_repeat_drift_rejects_only_affected_delta(self) -> None:
        records = passing_roster()
        records[13] = parsed_record(subject.make_invocations()[13],
                                   block_xors=101)
        result = subject.adjudicate(
            records, SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
            SYNTHETIC_DIAGONAL_CONTROL)
        self.assertNotIn(2, result["survivors"])
        self.assertEqual(
            result["delta_results"][1]["repeated_work_equal_coordinate_count"],
            2)
        self.assertEqual(
            result["delta_results"][1]["repeated_exact_coordinate_count"], 2)

    def test_configuration_failure_invalidates_screen_in_either_repeat(
            self) -> None:
        invocations = subject.make_invocations()
        for index in (6, 7):
            records = passing_roster()
            records[index] = configuration_failure_record(invocations[index])
            with self.subTest(repetition=index - 6), \
                    self.assertRaisesRegex(subject.Stage0Error, "precode p"):
                subject.adjudicate(
                    records, SYNTHETIC_PRECODE_SEEDS,
                    SYNTHETIC_PACKET_SEEDS, SYNTHETIC_DIAGONAL_CONTROL)

    def test_historical_precode_seed_mismatch_is_invalid(self) -> None:
        wrong = dict(SYNTHETIC_PRECODE_SEEDS)
        wrong[1] = "0x0000000000000001"
        with self.assertRaisesRegex(subject.Stage0Error, "precode seed"):
            subject.adjudicate(
                passing_roster(), wrong, SYNTHETIC_PACKET_SEEDS,
                SYNTHETIC_DIAGONAL_CONTROL)

    def test_historical_packet_seed_or_reported_q_seed_mismatch_is_invalid(
            self) -> None:
        wrong = dict(SYNTHETIC_PACKET_SEEDS)
        wrong[1] = "0x00000001"
        with self.assertRaisesRegex(subject.Stage0Error, "packet seed"):
            subject.adjudicate(
                passing_roster(), SYNTHETIC_PRECODE_SEEDS, wrong,
                SYNTHETIC_DIAGONAL_CONTROL)
        records = passing_roster()
        records[6]["effective_packet_seed"] = "0x00000000"
        records[6]["deterministic_projection_sha256"] = subject.sha256_json(
            subject.deterministic_projection(records[6]))
        with self.assertRaisesRegex(subject.Stage0Error, "packet seed"):
            subject.adjudicate(
                records, SYNTHETIC_PRECODE_SEEDS, SYNTHETIC_PACKET_SEEDS,
                SYNTHETIC_DIAGONAL_CONTROL)


class V9EvidenceTests(unittest.TestCase):
    def test_pins_regular_file_and_detects_mutation_and_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "evidence"
            path.write_bytes(b"frozen\n")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with subject._open_historical_file(
                    path, 7, digest, "test evidence") as pinned:
                self.assertEqual(pinned.read_bytes(8), b"frozen\n")
                path.write_bytes(b"changed\n")
                with self.assertRaises(subject.Stage0Error):
                    pinned.receipt()
            target = root / "target"
            target.write_bytes(b"x")
            link = root / "link"
            link.symlink_to(target)
            with self.assertRaises(subject.Stage0Error):
                subject._open_historical_file(
                    link, 1, hashlib.sha256(b"x").hexdigest(), "link")

    def test_loads_frozen_v9_shape_and_weak_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            attempts_data = b"synthetic attempts\n"
            derivations = [
                {"K": K, "selected_attempt": attempt}
                for K, attempt in zip(subject.V9_K_VALUES,
                                      subject.V9_ATTEMPTS)
            ]
            input_value = {
                "contract_sha256": subject.EXPECTED_V9_CONTRACT_SHA256,
                "attempt_selection_stream_sha256":
                    hashlib.sha256(attempts_data).hexdigest(),
                "selected_attempts": list(subject.V9_ATTEMPTS),
            }
            input_value["input_sha256"] = subject.self_hash(
                input_value, "input_sha256")
            input_data = (subject.canonical_json(input_value) + "\n").encode(
                "ascii")
            weak_ordinals = iter((10, 11, 300, 301, 800, 801))
            weak_by_ordinal = {}
            for coordinate_ordinal, coordinate in enumerate(subject.COORDINATES):
                for repetition in subject.REPETITIONS:
                    ordinal = next(weak_ordinals)
                    weak_by_ordinal[ordinal] = {
                        "ordinal": ordinal,
                        "arm": "candidate_two07_mix2", "weak": True,
                        "K": coordinate.K, "loss_root": coordinate.root,
                        "schedule": coordinate.schedule,
                        "observation_index": repetition,
                        "attempt_mode": "exact",
                        "construction_attempt": coordinate.attempt,
                        "mix_count": 2, "dense_anchor_layout": "two07",
                        "candidate_profile_sha256":
                            subject.EXPECTED_CANDIDATE_PROFILE_SHA256,
                        "success": 0, "rank_fail": 1, "error": 0,
                        "effective_precode_seed":
                            "0x{:016x}".format(coordinate_ordinal + 1),
                        "effective_packet_seed":
                            "0x{:08x}".format(coordinate_ordinal + 1),
                        "staircase":
                            subject.STAIRCASE_BY_K[coordinate.K],
                        "source_hits":
                            3 if coordinate.K >= 10000 else 2,
                        "inactivated_columns": 20 + coordinate_ordinal,
                        "block_xors": 100 + coordinate_ordinal,
                        "gf256_muladds": 50 + coordinate_ordinal,
                    }
            result_rows = []
            for ordinal in range(1080):
                row = weak_by_ordinal.get(
                    ordinal, {"ordinal": ordinal, "arm": "control",
                              "weak": False})
                result_rows.append(
                    (subject.canonical_json(row) + "\n").encode("ascii"))
            results_data = b"".join(result_rows)
            summary = {
                "contract_sha256": subject.EXPECTED_V9_CONTRACT_SHA256,
                "source_git_commit": subject.EXPECTED_V9_SOURCE_COMMIT,
                "repair_worker_binary_sha256":
                    subject.EXPECTED_V9_WORKER_SHA256,
                "candidate_profile_sha256":
                    subject.EXPECTED_CANDIDATE_PROFILE_SHA256,
                "input_sha256": input_value["input_sha256"],
                "input_file_sha256": hashlib.sha256(input_data).hexdigest(),
                "attempt_selection_stream_sha256":
                    hashlib.sha256(attempts_data).hexdigest(),
                "result_stream_sha256":
                    hashlib.sha256(results_data).hexdigest(),
                "record_count": 1080, "candidate_weak_coordinates": 3,
                "candidate_weak_observations": 6, "disposition": "REJECT",
            }
            summary["summary_sha256"] = subject.self_hash(
                summary, "summary_sha256")
            summary_data = (subject.canonical_json(summary) + "\n").encode(
                "ascii")
            values = {
                "promotion-short-screen-attempts.jsonl": attempts_data,
                "promotion-short-screen-input.json": input_data,
                "promotion-short-screen-results.jsonl": results_data,
                "promotion-short-screen-summary.json": summary_data,
            }
            for name, data in values.items():
                (root / name).write_bytes(data)
            file_specs = {
                name: (len(data), hashlib.sha256(data).hexdigest())
                for name, data in values.items()
            }
            with mock.patch.object(subject, "V9_FILES", file_specs), \
                    mock.patch.object(
                        subject, "EXPECTED_V9_INPUT_SELF_SHA256",
                        input_value["input_sha256"]), \
                    mock.patch.object(
                        subject, "EXPECTED_V9_SUMMARY_SELF_SHA256",
                        summary["summary_sha256"]), \
                    mock.patch.object(
                        subject.screen, "parse_derivation_stream",
                        return_value=derivations):
                with subject.load_v9_evidence(root) as evidence:
                    self.assertEqual(evidence.expected_precode_seeds, {
                        0: "0x0000000000000001",
                        1: "0x0000000000000002",
                        2: "0x0000000000000003",
                    })
                    self.assertEqual(evidence.expected_packet_seeds, {
                        0: "0x00000001",
                        1: "0x00000002",
                        2: "0x00000003",
                    })
                    self.assertEqual(
                        evidence.expected_diagonal_control[0]["block_xors"],
                        100)
                    self.assertEqual(
                        evidence.receipt_value["selected_attempts"],
                        list(subject.V9_ATTEMPTS))


class TransportTests(unittest.TestCase):
    def test_raw_process_uses_pinned_descriptor_and_frozen_env(self) -> None:
        invocation = subject.Invocation(6, 1, 0, 0)
        pinned = mock.Mock(path=Path("/frozen/bench"), descriptor=17)
        completed = mock.Mock(returncode=0, stdout=b"stdout\n", stderr=b"")
        with mock.patch.object(subject.subprocess, "run",
                               return_value=completed) as run:
            result = subject._run_raw(invocation, pinned)
        self.assertEqual(result.returncode, 0)
        argv, options = run.call_args
        self.assertEqual(argv[0], invocation.argv(pinned.path))
        self.assertEqual(options["executable"], "/proc/self/fd/17")
        self.assertEqual(options["pass_fds"], (17,))
        self.assertEqual(options["env"], subject.CHILD_ENVIRONMENT)
        self.assertTrue(options["start_new_session"])

    def test_execute_roster_parses_exactly_1536_processes(self) -> None:
        invocations = subject.make_invocations()
        results = [process_result(
                       item, stdout_for(
                           item, "rank_fail" if item.delta == 0 else "success"))
                   for item in invocations]
        with mock.patch.object(subject, "_run_raw",
                               side_effect=results) as run, \
                mock.patch.object(
                    subject, "parse_process_result",
                    wraps=subject.parse_process_result) as parse:
            records = subject._execute_roster(
                mock.sentinel.pinned, SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA,
                SYNTHETIC_DIAGONAL_CONTROL)
        self.assertEqual(len(records), 1536)
        self.assertEqual(run.call_count, 1536)
        self.assertEqual(parse.call_count, 1536)

    def test_execute_roster_stops_after_six_drifting_controls(self) -> None:
        invocations = subject.make_invocations()
        results = [process_result(
                       item, stdout_for(
                           item, "success" if index == 0 else "rank_fail"))
                   for index, item in enumerate(invocations[:6])]
        with mock.patch.object(subject, "_run_raw",
                               side_effect=results) as run, \
                self.assertRaisesRegex(subject.Stage0Error,
                                       "diagonal control"):
            subject._execute_roster(
                mock.sentinel.pinned, SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA,
                SYNTHETIC_DIAGONAL_CONTROL)
        self.assertEqual(run.call_count, 6)


class SourceReceiptTests(unittest.TestCase):
    @staticmethod
    def _git_result(stdout: bytes = b"", stderr: bytes = b"",
                    returncode: int = 0):
        return mock.Mock(
            stdout=stdout, stderr=stderr, returncode=returncode)

    def _responses(self, second_blob: bytes = b"beta\n",
                   final_head: str = SOURCE_COMMIT):
        return [
            self._git_result((SOURCE_COMMIT + "\n").encode("ascii")),
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
        self.assertEqual(receipt["source_git_commit"], SOURCE_COMMIT)
        self.assertEqual(receipt["source_files"], {
            "one": hashlib.sha256(b"alpha\n").hexdigest(),
            "two": hashlib.sha256(b"beta\n").hexdigest(),
        })
        self.assertEqual(
            run.call_args_list[3].args[0],
            ["git", "cat-file", "blob", SOURCE_COMMIT + ":one"])
        self.assertEqual(
            run.call_args_list[4].args[0],
            ["git", "cat-file", "blob", SOURCE_COMMIT + ":two"])
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
                        self.assertRaisesRegex(subject.Stage0Error, message):
                    subject._source_receipt()


class StartupAndPublicationTests(unittest.TestCase):
    def test_stale_description_stops_before_workload(self) -> None:
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.return_value = bench_receipt()
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(
                    subject.screen, "read_bench_description",
                    side_effect=subject.screen.ShortScreenError("stale")) \
                as describe, \
                mock.patch.object(subject, "_execute_roster") as execute, \
                self.assertRaisesRegex(subject.screen.ShortScreenError,
                                       "stale"):
            subject._run_stage0_pinned(
                pinned, Path(temporary) / "out", source_receipt(),
                subject.contract_description(), FakeV9())
        describe.assert_called_once_with(
            Path("/pinned/bench"), SOURCE_COMMIT, 17)
        execute.assert_not_called()

    def test_publication_is_atomic_and_contains_v9_receipt(self) -> None:
        events = []
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.side_effect = lambda: bench_receipt()

        def description(*_args):
            events.append("describe")
            return {"schema": subject.BENCH_DESCRIPTION_SCHEMA,
                    "source_git_commit": SOURCE_COMMIT}

        def execute(*_args):
            events.append("work")
            return passing_roster()

        def publish(source: Path, destination: Path):
            events.append("publish")
            summary = json.loads(
                (source / subject.SUMMARY_NAME).read_text(encoding="ascii"))
            self.assertEqual(summary["record_count"], 1536)
            self.assertEqual(summary["stage1_candidate_delta"], 1)
            self.assertEqual(summary["v9_prerequisite"]["receipt_sha256"],
                             "4" * 64)
            os.rename(source, destination)

        def fsync(path: Path):
            events.append(
                "fsync-staging" if path.name.startswith(
                    ".wh2-mix2-packet-offset-stage0-") else "fsync-parent")

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  side_effect=description), \
                mock.patch.object(subject, "_execute_roster",
                                  side_effect=execute), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(
                    subject.screen, "_publish_directory_noreplace",
                    side_effect=publish), \
                mock.patch.object(subject, "_fsync_directory",
                                  side_effect=fsync):
            output = Path(temporary) / "out"
            summary = subject._run_stage0_pinned(
                pinned, output, source_receipt(),
                subject.contract_description(), FakeV9())
            self.assertTrue(output.is_dir())
            self.assertEqual(summary["record_count"], 1536)
        self.assertEqual(events, [
            "describe", "work", "fsync-staging", "publish", "fsync-parent"])

    def test_publisher_never_replaces_existing_destination(self) -> None:
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

    def test_run_stage0_pins_v9_before_binary(self) -> None:
        events = []
        v9 = mock.MagicMock()
        v9.__enter__.return_value = mock.sentinel.v9
        pinned = mock.MagicMock()
        pinned.__enter__.return_value = mock.sentinel.pinned
        expected = {"disposition": "PASS"}

        def load(*_args):
            events.append("v9")
            return v9

        def open_binary(*_args):
            events.append("binary")
            return pinned

        def attest_source():
            events.append("source")
            return source_receipt()

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject, "load_v9_evidence",
                                  side_effect=load) as load_v9, \
                mock.patch.object(subject, "_source_receipt",
                                  side_effect=attest_source), \
                mock.patch.object(subject.screen, "_open_binary",
                                  side_effect=open_binary) as open_bin, \
                mock.patch.object(subject, "_run_stage0_pinned",
                                  return_value=expected) as run_pinned:
            output = Path(temporary) / "out"
            actual = subject.run_stage0(
                Path("/requested/bench"), Path("/v9"), output)
        self.assertEqual(actual, expected)
        self.assertEqual(events, ["v9", "source", "binary"])
        load_v9.assert_called_once_with(Path("/v9"))
        open_bin.assert_called_once_with(
            Path("/requested/bench"), "wirehair_v2_bench")
        run_pinned.assert_called_once()

    def test_invalid_v9_stops_before_source_binary_or_workload(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(
                    subject, "load_v9_evidence",
                    side_effect=subject.Stage0Error("invalid v9")) as load_v9, \
                mock.patch.object(subject, "_source_receipt") as source, \
                mock.patch.object(subject.screen, "_open_binary") as binary, \
                mock.patch.object(subject, "_execute_roster") as execute, \
                self.assertRaisesRegex(subject.Stage0Error, "invalid v9"):
            subject.run_stage0(
                Path("/requested/bench"), Path("/v9"),
                Path(temporary) / "out")
        load_v9.assert_called_once_with(Path("/v9"))
        source.assert_not_called()
        binary.assert_not_called()
        execute.assert_not_called()


if __name__ == "__main__":
    unittest.main()
