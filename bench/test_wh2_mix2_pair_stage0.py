#!/usr/bin/env python3
"""Pure mocked tests for the sealed MIX2 iterator-pair Stage 0."""

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
import wh2_mix2_pair_stage0 as subject


SOURCE_COMMIT = "1" * 40
BENCH_SHA = "2" * 64
SOURCE_SHA = "3" * 64


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
        "source_git_commit": SOURCE_COMMIT, "mix_pair": invocation.pair,
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
        "N": str(coordinate.K),
        "bb": "2",
        "heavy_family": "periodic",
        "mix_count": "2",
        "staircase": str(subject.STAIRCASE_BY_K[coordinate.K]),
        "binary_dense_rows": "12",
        "gf256_heavy_rows": "12",
        "source_hits": "3" if coordinate.K >= 10000 else "2",
        "dense_identity_corner": "0",
        "overhead": "0",
        "trials": "1",
        "success": str(success),
        "rank_fail": str(rank_fail),
        "error": str(error),
        "fail_rate": "1.00000000" if weak else "0.00000000",
        "inact_mu": "20.000",
        "inact_max": "20",
        "binary_def_mu": "{}.000".format(deficiency),
        "binary_def_max": str(deficiency),
        "heavy_gain_mu": "{}.000".format(gain),
        "heavy_gain_min": str(gain),
        "heavy_shortfall": (
            "1" if rank_fail and deficiency <= 12 and gain < deficiency
            else "0"),
        "solve_ms_mu": "1.000",
        "build_ms_mu": "1.000",
        "peel_ms_mu": "1.000",
        "project_ms_mu": "1.000",
        "residual_ms_mu": "1.000",
        "backsub_ms_mu": "1.000",
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


def stdout_for(invocation: subject.Invocation,
               outcome: str = "success",
               block_xors: int = 100,
               gf256_muladds: int = 50) -> bytes:
    values = csv_values(invocation, outcome, block_xors, gf256_muladds)
    lines = [
        metadata_line(invocation),
        ",".join(subject.CSV_HEADER),
        ",".join(values[field] for field in subject.CSV_HEADER),
    ]
    return ("\n".join(lines) + "\n").encode("ascii")


def process_result(invocation: subject.Invocation,
                   stdout: bytes,
                   returncode: int = 0,
                   stderr: bytes = b"") -> subject.ProcessResult:
    return subject.ProcessResult(
        invocation=invocation,
        command_sha256=subject.sha256_json({
            "argv": invocation.argv(Path("wirehair_v2_bench"))[1:],
        }),
        stdout_sha256=hashlib.sha256(stdout).hexdigest(),
        stderr_sha256=hashlib.sha256(stderr).hexdigest(),
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
    )


def parsed_record(invocation: subject.Invocation,
                  outcome: str = "success",
                  block_xors: int = 100,
                  gf256_muladds: int = 50) -> dict:
    return dict(subject.parse_process_result(
        process_result(
            invocation,
            stdout_for(invocation, outcome, block_xors, gf256_muladds)),
        SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA))


def successful_roster() -> list:
    return [parsed_record(invocation)
            for invocation in subject.make_invocations()]


def source_receipt() -> dict:
    return {
        "source_git_commit": SOURCE_COMMIT,
        "source_receipt_sha256": SOURCE_SHA,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "clean_status_scope": [],
        "source_files": {},
    }


def bench_receipt() -> dict:
    return {"path": "/pinned/bench", "sha256": BENCH_SHA, "size": 1}


class ContractTests(unittest.TestCase):
    def test_frozen_contract_and_exact_cardinality(self) -> None:
        subject._validate_constants()
        self.assertEqual(subject.sha256_json(subject._contract_body()),
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(subject.contract_description()["contract_sha256"],
                         subject.EXPECTED_CONTRACT_SHA256)
        self.assertEqual(len(subject._contract_body()["domain"]["coordinates"]),
                         3)
        invocations = subject.make_invocations()
        self.assertEqual(len(invocations), 12)
        self.assertEqual(
            [(item.pair, item.coordinate_ordinal, item.repetition)
             for item in invocations],
            [(pair, coordinate, repetition)
             for pair in ("02", "12")
             for coordinate in range(3)
             for repetition in (0, 1)])
        self.assertFalse(
            set(subject.CONSUMED_ROOTS) &
            set(subject.UNCONSUMED_FINAL_ROOTS))

    def test_argv_is_the_frozen_full_payload_exact_attempt_cell(self) -> None:
        invocation = subject.Invocation(0, "02", 0, 0)
        self.assertEqual(invocation.argv(Path("/pinned/bench")), [
            "/pinned/bench", "precodefail", "--N", "2", "--bb-list",
            "2", "--overhead", "0", "--trials", "1", "--threads", "1",
            "--loss", "0.5", "--seed", "0x7ccd510f122fc160",
            "--schedule", "burst", "--heavy-family", "periodic",
            "--mix-count", "2", "--mix-pair", "02",
            "--binary-dense-rows", "12", "--gf256-heavy-rows", "12",
            "--dense-anchors", "two07", "--paired-overhead-stream",
            "--full-payload-solve", "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", "9", "--exact-packet-attempt", "9",
            "--construction-seed-basis", "production-profile",
        ])
        with self.assertRaises(subject.Stage0Error):
            subject.Invocation(1, "02", 0, 0)


class ParserTests(unittest.TestCase):
    def test_accepts_exact_mix_pair_metadata_and_success_row(self) -> None:
        invocation = subject.Invocation(6, "12", 0, 0)
        record = parsed_record(invocation)
        self.assertEqual(record["pair"], "12")
        self.assertEqual(record["outcome"], "success")
        self.assertEqual(record["block_xors"], 100)
        self.assertEqual(record["gf256_muladds"], 50)
        self.assertEqual(
            record["deterministic_projection_sha256"],
            subject.sha256_json(subject.deterministic_projection(record)))

    def test_rejects_missing_reordered_or_wrong_mix_pair_metadata(self) -> None:
        invocation = subject.Invocation(0, "02", 0, 0)
        valid = stdout_for(invocation).decode("ascii").splitlines()
        cases = []
        cases.append((" ".join(valid[0].split(" ")[:-1]) + "\n" +
                      "\n".join(valid[1:]) + "\n").encode("ascii"))
        tokens = valid[0].split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        cases.append((" ".join(tokens) + "\n" +
                      "\n".join(valid[1:]) + "\n").encode("ascii"))
        cases.append(stdout_for(invocation).replace(b"mix_pair=02",
                                                     b"mix_pair=12"))
        for stdout in cases:
            with self.subTest(stdout=stdout), \
                    self.assertRaises(subject.Stage0Error):
                subject.parse_process_result(
                    process_result(invocation, stdout), SOURCE_COMMIT,
                    BENCH_SHA, SOURCE_SHA)

    def test_rejects_noncanonical_work_and_inconsistent_outcome(self) -> None:
        invocation = subject.Invocation(0, "02", 0, 0)
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

    def test_canonical_configuration_failure_is_recorded_for_gate(self) -> None:
        invocation = subject.Invocation(0, "02", 0, 0)
        stdout = (metadata_line(invocation) + "\n" +
                  ",".join(subject.CSV_HEADER) + "\n").encode("ascii")
        stderr = (
            "precodefail configuration failed N=2 bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt=9 packet_attempt=9 result=2\n").encode("ascii")
        record = subject.parse_process_result(
            process_result(invocation, stdout, 2, stderr), SOURCE_COMMIT,
            BENCH_SHA, SOURCE_SHA)
        self.assertTrue(record["configuration_failure"])
        self.assertEqual(record["outcome"], "configuration_failure")
        self.assertIsNone(record["block_xors"])
        with self.assertRaises(subject.Stage0Error):
            subject.parse_process_result(
                process_result(invocation, stdout, 2, stderr + b"extra\n"),
                SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA)


class GateTests(unittest.TestCase):
    def test_both_pairs_survive_exact_repeated_success(self) -> None:
        result = subject.adjudicate(successful_roster())
        self.assertEqual(result["disposition"], "PASS")
        self.assertEqual(result["survivors"], ["02", "12"])
        self.assertEqual(result["stage1_candidate"], "02")
        self.assertTrue(all(value["disposition"] == "PASS"
                            for value in result["pair_results"]))

    def test_priority_fallback_is_frozen_without_post_result_selection(self) \
            -> None:
        invocations = subject.make_invocations()
        pair02_fails = successful_roster()
        for index in (0, 1):
            pair02_fails[index] = parsed_record(
                invocations[index], "rank_fail")
        result = subject.adjudicate(pair02_fails)
        self.assertEqual(result["survivors"], ["12"])
        self.assertEqual(result["stage1_candidate"], "12")

        pair12_fails = successful_roster()
        for index in (6, 7):
            pair12_fails[index] = parsed_record(
                invocations[index], "rank_fail")
        result = subject.adjudicate(pair12_fails)
        self.assertEqual(result["survivors"], ["02"])
        self.assertEqual(result["stage1_candidate"], "02")

    def test_failure_or_repeated_work_drift_rejects_only_affected_pair(self) -> None:
        records = successful_roster()
        # Pair 02 coordinate zero fails in both repetitions.
        for index in (0, 1):
            records[index] = parsed_record(
                subject.make_invocations()[index], "rank_fail")
        # Pair 12 coordinate zero has nondeterministic work.
        records[7] = parsed_record(subject.make_invocations()[7],
                                   block_xors=101)
        result = subject.adjudicate(records)
        self.assertEqual(result["disposition"], "REJECT")
        self.assertEqual(result["survivors"], [])
        self.assertIsNone(result["stage1_candidate"])
        pair02, pair12 = result["pair_results"]
        self.assertEqual(pair02["rank_failures"], 2)
        self.assertEqual(pair12["repeated_work_equal_coordinate_count"], 2)
        self.assertEqual(pair12["repeated_exact_coordinate_count"], 2)


class TransportTests(unittest.TestCase):
    def test_raw_process_uses_only_the_pinned_descriptor_and_frozen_env(self) -> None:
        invocation = subject.Invocation(0, "02", 0, 0)
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
        self.assertIs(options["stdin"], subject.subprocess.DEVNULL)
        self.assertTrue(options["start_new_session"])

    def test_execute_roster_parses_exactly_twelve_processes(self) -> None:
        invocations = subject.make_invocations()
        results = [process_result(item, stdout_for(item))
                   for item in invocations]
        with mock.patch.object(subject, "_run_raw",
                               side_effect=results) as run, \
                mock.patch.object(
                    subject, "parse_process_result",
                    side_effect=lambda value, *_: {
                        **value.invocation.identity(), "parsed": True,
                    }) as parse:
            records = subject._execute_roster(
                mock.sentinel.pinned, SOURCE_COMMIT, BENCH_SHA, SOURCE_SHA)
        self.assertEqual(len(records), 12)
        self.assertEqual(run.call_count, 12)
        self.assertEqual(parse.call_count, 12)


class StartupAndPublicationTests(unittest.TestCase):
    def test_stale_description_stops_before_any_workload(self) -> None:
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
                subject.contract_description())
        describe.assert_called_once_with(
            Path("/pinned/bench"), SOURCE_COMMIT, 17)
        execute.assert_not_called()

    def test_describe_precedes_work_and_publication_is_atomic_directory(self) -> None:
        events = []
        pinned = mock.Mock(path=Path("/pinned/bench"), descriptor=17)
        pinned.receipt.side_effect = lambda: bench_receipt()

        def description(*_args):
            events.append("describe")
            return {
                "schema": subject.BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": SOURCE_COMMIT,
            }

        def execute(*_args):
            self.assertIn("describe", events)
            events.append("work")
            return successful_roster()

        def publish(source: Path, destination: Path):
            events.append("publish")
            self.assertTrue((source / subject.RECORD_NAME).is_file())
            self.assertTrue((source / subject.SUMMARY_NAME).is_file())
            summary = json.loads(
                (source / subject.SUMMARY_NAME).read_text(encoding="ascii"))
            self.assertEqual(summary["record_count"], 12)
            self.assertEqual(summary["survivors"], ["02", "12"])
            os.rename(source, destination)

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.screen, "read_bench_description",
                                  side_effect=description), \
                mock.patch.object(subject, "_execute_roster",
                                  side_effect=execute), \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(
                    subject.screen, "_publish_directory_noreplace",
                    side_effect=publish):
            output = Path(temporary) / "out"
            summary = subject._run_stage0_pinned(
                pinned, output, source_receipt(),
                subject.contract_description())
            self.assertTrue(output.is_dir())
            self.assertEqual(summary["record_count"], 12)
        self.assertEqual(events, ["describe", "work", "publish"])

    def test_run_stage0_uses_the_binary_pinning_context(self) -> None:
        pinned = mock.MagicMock()
        pinned.__enter__.return_value = mock.sentinel.pinned
        pinned.__exit__.return_value = None
        expected = {"disposition": "PASS"}
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject, "_source_receipt",
                                  return_value=source_receipt()), \
                mock.patch.object(subject.screen, "_open_binary",
                                  return_value=pinned) as open_binary, \
                mock.patch.object(subject, "_run_stage0_pinned",
                                  return_value=expected) as run_pinned:
            output = Path(temporary) / "out"
            actual = subject.run_stage0(Path("/requested/bench"), output)
        self.assertEqual(actual, expected)
        open_binary.assert_called_once_with(
            Path("/requested/bench"), "wirehair_v2_bench")
        run_pinned.assert_called_once()


if __name__ == "__main__":
    unittest.main()
