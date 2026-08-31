#!/usr/bin/env python3
"""Pure mocked tests for the sealed MIX2 packet-offset Stage 2 controller."""

from __future__ import annotations

import copy
from contextlib import redirect_stderr
import hashlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_packet_offset_stage2 as subject


CONTROLLER_COMMIT = "1" * 40
SOURCE_SHA = "2" * 64
EMPTY_SHA = hashlib.sha256(b"").hexdigest()


def baseline_projection(coordinate_ordinal: int) -> dict:
    coordinate = subject.stage1.Coordinate(coordinate_ordinal)
    success = coordinate_ordinal not in subject.ORIGINAL_WEAK_COORDINATES
    deficiency = 10 if success else 12
    gain = deficiency if success else deficiency - 1
    return {
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
        "binary_deficiency": deficiency,
        "gf256_rank_gain": gain,
        "inactivated_columns": 20,
        "heavy_shortfall": int(not success),
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


def csv_values(invocation: subject.Invocation, outcome: str = "success",
               block_xors: int = 100, gf256_muladds: int = 50,
               inactivated: int = 20) -> dict:
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
    baseline = baseline_projection(invocation.coordinate_ordinal)
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
        "inact_mu": "{}.000".format(inactivated),
        "inact_max": str(inactivated),
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
        "effective_precode_seed": baseline["effective_precode_seed"],
        "effective_packet_seed": subject.stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], coordinate.attempt,
            invocation.delta),
    }


def stdout_for(invocation: subject.Invocation, outcome: str = "success",
               block_xors: int = 100, gf256_muladds: int = 50,
               inactivated: int = 20) -> bytes:
    values = csv_values(invocation, outcome, block_xors, gf256_muladds,
                        inactivated)
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


def parsed_pair(invocation: subject.Invocation,
                outcome: str = "success") -> tuple[dict, dict]:
    record, raw = subject.parse_process_result(
        process_result(invocation, stdout_for(invocation, outcome)),
        CONTROLLER_COMMIT, SOURCE_SHA,
        baseline_projection(invocation.coordinate_ordinal))
    return dict(record), dict(raw)


def synthetic_record(invocation: subject.Invocation, outcome: str,
                     work: tuple[int, int, int] = (100, 50, 20)) -> dict:
    block_xors, muladds, inactivated = work
    if outcome == "success":
        success, rank_fail, error, deficiency, gain = True, False, False, 10, 10
    elif outcome == "rank_fail":
        success, rank_fail, error, deficiency, gain = False, True, False, 12, 11
    elif outcome == "error":
        success, rank_fail, error, deficiency, gain = False, False, True, 0, 0
    else:
        raise AssertionError("bad outcome")
    coordinate = invocation.coordinate
    baseline = baseline_projection(invocation.coordinate_ordinal)
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
        "baseline_projection_sha256": subject.sha256_json(baseline),
        "staircase": subject.stage1.STAIRCASE_BY_K[coordinate.K],
        "source_hits": 3 if coordinate.K >= 10000 else 2,
        "outcome": outcome, "success": success, "rank_fail": rank_fail,
        "error": error, "configuration_failure": False,
        "weak": not success,
        "inactivated_columns": inactivated,
        "binary_deficiency": deficiency, "gf256_rank_gain": gain,
        "heavy_shortfall": int(rank_fail),
        "first_rank_fail": 0 if rank_fail else -1,
        "failure_trials": [0] if not success else [],
        "effective_precode_seed": baseline["effective_precode_seed"],
        "effective_packet_seed": subject.stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], coordinate.attempt,
            invocation.delta),
        "block_xors": block_xors, "gf256_muladds": muladds,
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


def make_evidence(selected_delta: int = 3,
                  selected_work: tuple[int, int, int] = (100, 50, 20)) \
        -> subject.HistoricalEvidence:
    baseline = {coordinate: baseline_projection(coordinate)
                for coordinate in subject.CONFIRMATION_COORDINATES}
    controls = {}
    for invocation in subject.control_invocations():
        projection = subject.scientific_projection(
            synthetic_record(invocation, "rank_fail"))
        controls[invocation.coordinate_ordinal] = (projection, projection)
    stage0_success = {}
    for delta in subject.SCREEN_DELTAS:
        work = selected_work if delta == selected_delta else (100, 50, 20)
        for coordinate in subject.ORIGINAL_WEAK_COORDINATES:
            invocation = subject.Invocation(
                subject.PRE_CONFIRMATION_OBSERVATION_COUNT +
                subject.CONFIRMATION_COORDINATES.index(coordinate),
                subject.PHASE_CONFIRMATION, delta, coordinate)
            stage0_success[(delta, coordinate)] = \
                subject.scientific_projection(
                    synthetic_record(invocation, "success", work))
    return subject.HistoricalEvidence(
        None, {}, baseline, controls, stage0_success,
        {"stage0": {}, "stage1": {}})


def passing_roster(
        selected_delta: int = 3,
        selected_work: tuple[int, int, int] = (100, 50, 20),
        also_passing: tuple[int, ...] = (5,)) \
        -> tuple[list[dict], subject.HistoricalEvidence]:
    evidence = make_evidence(selected_delta, selected_work)
    records = [synthetic_record(invocation, "rank_fail")
               for invocation in subject.control_invocations()]
    passing = {selected_delta, *also_passing}
    for invocation in subject.screen_invocations():
        outcome = "success" if invocation.delta in passing or \
            invocation.coordinate_ordinal != subject.INTRODUCTION_COORDINATES[0] \
            else "rank_fail"
        work = selected_work if invocation.delta == selected_delta \
            else (100, 50, 20)
        records.append(synthetic_record(invocation, outcome, work))
    for invocation in subject.confirmation_invocations(selected_delta):
        records.append(synthetic_record(invocation, "success", selected_work))
    return records, evidence


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_literal_rosters(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(len(subject.SURVIVOR_DELTAS), 218)
        self.assertEqual(len(subject.SCREEN_DELTAS), 217)
        self.assertEqual(subject.sha256_json(list(subject.SCREEN_DELTAS)),
                         subject.EXPECTED_SCREEN_DELTAS_SHA256)
        self.assertEqual(len(subject.INTRODUCTION_COORDINATES), 44)
        self.assertEqual(subject.ORIGINAL_WEAK_COORDINATES, (30, 573, 820))
        self.assertFalse(set(subject.INTRODUCTION_COORDINATES) &
                         set(subject.ORIGINAL_WEAK_COORDINATES))
        self.assertEqual(subject.MAX_OBSERVATION_COUNT, 9639)
        candidate = subject._contract_body()["candidate"]
        self.assertFalse(candidate["production_split_pq_implemented"])
        self.assertTrue(candidate["test_hook_only"])

    def test_exact_local_to_full_mapping_and_no_fresh_roots(self) -> None:
        self.assertEqual(tuple(
            subject._stage0_global_coordinate(local) for local in range(3)),
            subject.ORIGINAL_WEAK_COORDINATES)
        roots = {subject.stage1.Coordinate(coordinate).root
                 for coordinate in subject.CONFIRMATION_COORDINATES}
        self.assertFalse(roots & set(subject.screen.FINAL_VALIDATION_ROOTS))

    def test_canonical_invocation_order_and_arbitrary_delta_wrap(self) -> None:
        controls = subject.control_invocations()
        screens = subject.screen_invocations()
        self.assertEqual(len(controls), 44)
        self.assertEqual(len(screens), 9548)
        self.assertEqual((screens[0].delta, screens[0].coordinate_ordinal),
                         (3, 1))
        self.assertEqual((screens[-1].delta, screens[-1].coordinate_ordinal),
                         (255, 1021))
        wrap = subject.Invocation(
            subject.CONTROL_OBSERVATION_COUNT +
            subject.SCREEN_DELTAS.index(253) * 44 +
            subject.INTRODUCTION_COORDINATES.index(146),
            subject.PHASE_SCREEN, 253, 146)
        self.assertEqual(wrap.coordinate.attempt, 49)
        self.assertEqual(wrap.coordinate.packet_attempt, 46)
        self.assertEqual(subject.stage0._packet_seed_for_offset(
            "0x12345678", 49, 253), "0x378de94d")
        self.assertIn("--exact-packet-attempt", wrap.argv(Path("bench")))

    def test_parser_requires_explicit_sealed_directories(self) -> None:
        parser = subject._parser()
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args(["run", "--bench", "b", "--output-dir", "o"])
        parsed = parser.parse_args([
            "run", "--bench", "b", "--stage0-dir", "s0",
            "--stage1-dir", "s1", "--output-dir", "o"])
        self.assertEqual(parsed.stage1_dir, Path("s1"))


class ParserAndRawTests(unittest.TestCase):
    def test_arbitrary_delta_parser_and_raw_reparse(self) -> None:
        invocation = subject.screen_invocations()[0]
        record, raw = parsed_pair(invocation)
        self.assertTrue(record["success"])
        self.assertEqual(record["delta"], 3)
        self.assertEqual(record["packet_attempt"], 12)
        decoded = subject.validate_raw_reparse(
            raw, record, CONTROLLER_COMMIT, SOURCE_SHA,
            baseline_projection(invocation.coordinate_ordinal))
        self.assertEqual(decoded, len(stdout_for(invocation)))

    def test_tampered_raw_stdout_fails_reparse_equality(self) -> None:
        invocation = subject.screen_invocations()[0]
        record, raw = parsed_pair(invocation)
        tampered = copy.deepcopy(raw)
        decoded = subject._decode_raw(
            tampered["stdout_base64"], subject.MAX_STDOUT_BYTES, "stdout")
        decoded = decoded.replace(b"solve_ms_mu", b"solve_ms_mx", 1)
        tampered["stdout_base64"] = subject._canonical_b64(decoded)
        with self.assertRaises(subject.Stage2Error):
            subject.validate_raw_reparse(
                tampered, record, CONTROLLER_COMMIT, SOURCE_SHA,
                baseline_projection(invocation.coordinate_ordinal))

    def test_even_rehashed_raw_change_must_equal_parsed_record(self) -> None:
        invocation = subject.screen_invocations()[0]
        record, raw = parsed_pair(invocation)
        tampered = copy.deepcopy(raw)
        decoded = subject._decode_raw(
            tampered["stdout_base64"], subject.MAX_STDOUT_BYTES, "stdout")
        decoded = decoded.replace(b"solve_ms_mu,1.000", b"solve_ms_mu,2.000")
        # The textual replacement above may not match because this is CSV; use
        # a guaranteed canonical timing-token substitution instead.
        if decoded == subject._decode_raw(
                raw["stdout_base64"], subject.MAX_STDOUT_BYTES, "stdout"):
            decoded = decoded.replace(b",1.000,1.000,1.000,1.000,1.000,1.000,",
                                      b",2.000,1.000,1.000,1.000,1.000,1.000,",
                                      1)
        tampered["stdout_base64"] = subject._canonical_b64(decoded)
        tampered["stdout_sha256"] = subject.sha256_bytes(decoded)
        altered_record = copy.deepcopy(record)
        altered_record["stdout_sha256"] = tampered["stdout_sha256"]
        tampered["parsed_record_sha256"] = subject.sha256_json(altered_record)
        # Reparse equality is exact even though timings are not scientific
        # evidence; the persisted parsed record must name the changed bytes.
        self.assertEqual(subject.validate_raw_reparse(
            tampered, altered_record, CONTROLLER_COMMIT, SOURCE_SHA,
            baseline_projection(invocation.coordinate_ordinal)), len(decoded))
        with self.assertRaises(subject.Stage2Error):
            subject.validate_raw_reparse(
                tampered, record, CONTROLLER_COMMIT, SOURCE_SHA,
                baseline_projection(invocation.coordinate_ordinal))

    def test_wrong_effective_packet_seed_and_config_fail_closed(self) -> None:
        invocation = subject.screen_invocations()[0]
        stdout = stdout_for(invocation)
        expected = csv_values(invocation)["effective_packet_seed"].encode()
        with self.assertRaisesRegex(subject.Stage2Error, "packet seed"):
            subject.parse_process_result(
                process_result(invocation, stdout.replace(
                    expected, b"0x00000001")), CONTROLLER_COMMIT, SOURCE_SHA,
                baseline_projection(invocation.coordinate_ordinal))
        config_stdout = (metadata_line(invocation) + "\n" +
                         ",".join(subject.CSV_HEADER) + "\n").encode()
        config_stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=12 result=2\n").encode()
        record, _ = subject.parse_process_result(
            process_result(invocation, config_stdout, 2, config_stderr),
            CONTROLLER_COMMIT, SOURCE_SHA,
            baseline_projection(invocation.coordinate_ordinal))
        self.assertTrue(record["configuration_failure"])

    def test_run_raw_uses_pinned_descriptor_and_frozen_environment(self) -> None:
        invocation = subject.screen_invocations()[0]
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


class GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.records, cls.evidence = passing_roster(
            selected_work=(999, 888, 777))

    def test_pass_screens_every_delta_selects_lowest_without_work(self) -> None:
        verdict = subject.adjudicate(self.records, self.evidence)
        self.assertEqual(verdict["disposition"], "PASS")
        self.assertEqual(verdict["failure_control_rank_failures"], 44)
        self.assertEqual(verdict["screen_delta_count"], 217)
        self.assertEqual(verdict["screen_observation_count"], 9548)
        self.assertEqual(verdict["zero_weak_deltas"][:2], [3, 5])
        self.assertEqual(verdict["selected_delta"], 3)
        self.assertEqual(verdict["confirmation_observation_count"], 47)
        work = verdict["deduplicated_confirmation_work"]
        self.assertFalse(work["selection_used_work"])
        self.assertFalse(work["targeted_work_gate_applied"])
        self.assertGreater(work["selected_confirmation"]["gf256_muladds"],
                           work["pair01_baseline"]["gf256_muladds"])
        self.assertFalse(verdict["production_split_pq_implemented"])
        self.assertEqual(verdict["selection_leakage"], {
            "consulted_coordinate_count": 47,
            "untouched_coordinate_count": 1033,
            "next_replay_is_disjoint": False,
        })

    def test_no_survivor_rejects_without_confirmation(self) -> None:
        records = [copy.deepcopy(row)
                   for row in self.records[
                       :subject.PRE_CONFIRMATION_OBSERVATION_COUNT]]
        for delta in (3, 5):
            index = (subject.CONTROL_OBSERVATION_COUNT +
                     subject.SCREEN_DELTAS.index(delta) * 44)
            row = records[index]
            row.update({
                "outcome": "rank_fail", "success": False,
                "rank_fail": True, "weak": True,
                "binary_deficiency": 12, "gf256_rank_gain": 11,
                "heavy_shortfall": 1, "first_rank_fail": 0,
                "failure_trials": [0],
            })
            rehash(row)
        verdict = subject.adjudicate(records, self.evidence)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertIsNone(verdict["selected_delta"])
        self.assertEqual(verdict["confirmation_observation_count"], 0)
        self.assertTrue(verdict["fixed_offset_family_retired"])

    def test_control_must_match_both_sealed_repetitions(self) -> None:
        controls = dict(self.evidence.stage1_control_by_coordinate)
        coordinate = subject.INTRODUCTION_COORDINATES[0]
        damaged = dict(controls[coordinate][1])
        damaged["block_xors"] += 1
        controls[coordinate] = (controls[coordinate][0], damaged)
        evidence = subject.HistoricalEvidence(
            None, {}, self.evidence.baseline_by_coordinate, controls,
            self.evidence.stage0_success_by_delta_coordinate, {})
        with self.assertRaisesRegex(subject.Stage2Error, "control"):
            subject.adjudicate(self.records, evidence)

    def test_confirmation_mismatch_rejects_without_fallthrough(self) -> None:
        records = list(self.records)
        index = subject.PRE_CONFIRMATION_OBSERVATION_COUNT
        records[index] = copy.deepcopy(records[index])
        records[index]["block_xors"] += 1
        rehash(records[index])
        verdict = subject.adjudicate(records, self.evidence)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(verdict["selected_delta"], 3)
        self.assertEqual(
            verdict["confirmation_projection_mismatch_coordinate_ordinals"],
            [subject.CONFIRMATION_COORDINATES[0]])
        self.assertFalse(
            verdict["gates"]["selected_confirmation_projection_exact"])

    def test_original_confirmation_uses_normalized_stage0_projection(self) \
            -> None:
        records = list(self.records)
        coordinate = subject.ORIGINAL_WEAK_COORDINATES[0]
        index = (subject.PRE_CONFIRMATION_OBSERVATION_COUNT +
                 subject.CONFIRMATION_COORDINATES.index(coordinate))
        records[index] = copy.deepcopy(records[index])
        records[index]["gf256_muladds"] += 1
        rehash(records[index])
        verdict = subject.adjudicate(records, self.evidence)
        self.assertEqual(verdict["disposition"], "REJECT")
        self.assertEqual(verdict["selected_delta"], 3)
        self.assertEqual(
            verdict["confirmation_projection_mismatch_coordinate_ordinals"],
            [coordinate])

    def test_missing_extra_or_reordered_rows_fail_closed(self) -> None:
        cases = (
            self.records[:-1],
            self.records + [self.records[-1]],
            [self.records[1], self.records[0], *self.records[2:]],
        )
        for rows in cases:
            with self.subTest(count=len(rows)), \
                    self.assertRaises(subject.Stage2Error):
                subject.adjudicate(rows, self.evidence)

    def test_delta2_mismatch_stops_before_any_screen_candidate(self) -> None:
        first = subject.control_invocations()[0]
        success_result = process_result(first, stdout_for(first, "success"))
        with mock.patch.object(subject, "_run_raw",
                               return_value=success_result) as run_raw, \
                self.assertRaisesRegex(subject.Stage2Error, "control"):
            subject._execute_stage2(
                mock.Mock(), CONTROLLER_COMMIT, SOURCE_SHA, self.evidence)
        run_raw.assert_called_once()


class TransportTests(unittest.TestCase):
    def test_jsonl_writer_is_exclusive_bounded_and_canonical(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rows.jsonl"
            digest, count, size = subject._write_jsonl(
                path, ({"b": 2, "a": 1},), 100)
            self.assertEqual(path.read_bytes(), b'{"a":1,"b":2}\n')
            self.assertEqual(digest, hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual((count, size), (1, len(path.read_bytes())))
            with self.assertRaises(subject.Stage2Error):
                subject._write_jsonl(path, ({"a": 1},), 100)
            with self.assertRaises(subject.Stage2Error):
                subject._write_jsonl(Path(temporary) / "large", ({"x": "y"},),
                                     1)

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
            self.assertEqual((destination / "marker").read_bytes(),
                             b"destination")

    def test_complete_file_stream_audit_reparses_and_preserves_pair_order(self) \
            -> None:
        invocations = subject.control_invocations()[:2]
        pairs = [parsed_pair(invocation, "rank_fail")
                 for invocation in invocations]
        records = [pair[0] for pair in pairs]
        raw = [pair[1] for pair in pairs]
        evidence = make_evidence()
        expected_verdict = {"disposition": "REJECT"}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            records_path = root / "records.jsonl"
            raw_path = root / "raw.jsonl"
            records_sha, _, _ = subject._write_jsonl(records_path, records)
            raw_sha, _, _ = subject._write_jsonl(raw_path, raw)
            with mock.patch.object(subject, "adjudicate",
                                   return_value=expected_verdict):
                audit = subject._audit_output_streams(
                    records_path, raw_path, records_sha, raw_sha, 2,
                    CONTROLLER_COMMIT, SOURCE_SHA, evidence)
            self.assertEqual(audit["verdict"], expected_verdict)
            self.assertEqual(audit["decoded_bytes"], sum(
                len(stdout_for(invocation, "rank_fail"))
                for invocation in invocations))

            reversed_path = root / "raw-reversed.jsonl"
            reversed_sha, _, _ = subject._write_jsonl(
                reversed_path, reversed(raw))
            with mock.patch.object(subject, "adjudicate",
                                   return_value=expected_verdict), \
                    self.assertRaises(subject.Stage2Error):
                subject._audit_output_streams(
                    records_path, reversed_path, records_sha, reversed_sha, 2,
                    CONTROLLER_COMMIT, SOURCE_SHA, evidence)


if __name__ == "__main__":
    unittest.main()
