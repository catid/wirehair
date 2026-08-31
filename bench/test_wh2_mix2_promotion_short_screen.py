#!/usr/bin/env python3
"""Focused pure tests for the MIX2 production-basis confirmation."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_promotion_short_screen as subject


WORKER_SHA = "a" * 64
BENCH_SHA = "b" * 64
SOURCE_SHA = "c" * 64
SOURCE_COMMIT = "1" * 40


def worker_cell(ordinal: int, attempt: int, success: bool = True,
                roots: tuple = subject.SELECTION_ROOTS) -> dict:
    root_index = ordinal // len(subject.SCHEDULES)
    schedule_index = ordinal % len(subject.SCHEDULES)
    return {
        "attempted_candidates": 100 + ordinal,
        "cell_ordinal": ordinal,
        "construction_attempt": attempt,
        "decoded_extra": 0 if success else 1,
        "loss_ppm": subject.LOSS_PPM,
        "loss_root": roots[root_index],
        "outcome": "success" if success else "need_more_at_oh0",
        "result_code": 0,
        "root_index": root_index,
        "schedule": subject.SCHEDULES[schedule_index],
        "trace_sha256": hashlib.sha256(
            "cell:{}".format(ordinal).encode("ascii")).hexdigest(),
    }


def derivation(K: int, attempt: int = 0) -> dict:
    base_precode = "0x{:016x}".format(0x1234000000000000 + K)
    base_packet = "0x{:08x}".format(0x12000000 + K)
    return {
        "K": K,
        "base_packet_seed": base_packet,
        "base_precode_seed": base_precode,
        "candidate_profile_sha256": subject.candidate_profile_sha256(),
        "effective_packet_seed":
            subject._effective_packet_seed(base_packet, attempt),
        "effective_precode_seed":
            subject._effective_precode_seed(base_precode, attempt),
        "lower_attempt_failure_witnesses": [
            worker_cell(index % 18, index, False)
            for index in range(attempt)
        ],
        "mode": "derive",
        "ordinal": K - 2,
        "schema": subject.WORKER_DERIVATION_SCHEMA,
        "selected_attempt": attempt,
        "selected_successes": [
            worker_cell(index, attempt) for index in range(18)
        ],
        "source_sha256": hashlib.sha256(
            "source:{}".format(K).encode("ascii")).hexdigest(),
        "worker_binary_sha256": WORKER_SHA,
    }


def worker_description() -> dict:
    return {
        "binary_sha256": WORKER_SHA,
        "candidate_profile": subject.CANDIDATE_PROFILE,
        "candidate_profile_sha256": subject.candidate_profile_sha256(),
        "contract_schema": subject.WORKER_CONTRACT_SCHEMA,
        "derivation_schema": subject.WORKER_DERIVATION_SCHEMA,
        "protocol": "D ordinal K | V ordinal K attempt | Q",
        "schema": subject.WORKER_DESCRIPTION_SCHEMA,
        "source_git_commit": SOURCE_COMMIT,
        "validation_roster_schema": subject.VALIDATION_ROSTER_SCHEMA,
        "validation_roster_sha256": subject.VALIDATION_ROSTER_SHA256,
        "validation_schema": subject.WORKER_VALIDATION_SCHEMA,
        "worker_schema": subject.WORKER_SCHEMA,
    }


def bench_description() -> dict:
    return {
        "schema": subject.BENCH_DESCRIPTION_SCHEMA,
        "source_git_commit": SOURCE_COMMIT,
    }


def metadata(invocation: subject.Invocation) -> str:
    candidate = invocation.arm == "candidate_two07_mix2"
    values = {
        "trials": "1",
        "threads": "1",
        "loss": "0.5",
        "seed": invocation.root,
        "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07" if candidate else "disabled",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0",
        "overhead_stream": "paired",
        "full_payload_solve": "1",
        "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy",
        "exact_attempt_mode": "1" if candidate else "0",
        "exact_precode_attempt": str(invocation.attempt if candidate else 0),
        "exact_packet_attempt": str(invocation.attempt if candidate else 0),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "source_git_commit": SOURCE_COMMIT,
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(key, values[key])
        for key in (
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
        ))


def bench_stdout(invocation: subject.Invocation, record: dict,
                 *, weak: bool = False) -> bytes:
    candidate = invocation.arm == "candidate_two07_mix2"
    attempt = invocation.attempt if candidate else 0
    effective_precode = (record["effective_precode_seed"] if candidate else
                         "0x0000000000000100")
    effective_packet = (record["effective_packet_seed"] if candidate else
                        "0x00000100")
    values = {
        "N": str(invocation.K),
        "bb": "2",
        "heavy_family": "periodic",
        "mix_count": "2" if candidate else "3",
        "staircase": "2",
        "binary_dense_rows": "12",
        "gf256_heavy_rows": "12",
        "source_hits": "2",
        "dense_identity_corner": "0",
        "overhead": "0",
        "trials": "1",
        "success": "0" if weak else "1",
        "rank_fail": "1" if weak else "0",
        "error": "0",
        "fail_rate": "1.00000000" if weak else "0.00000000",
        "inact_mu": "10.000",
        "inact_max": "10",
        "binary_def_mu": "0.000",
        "binary_def_max": "0",
        "heavy_gain_mu": "0.000",
        "heavy_gain_min": "0",
        "heavy_shortfall": "0",
        "solve_ms_mu": "1.000",
        "build_ms_mu": "1.000",
        "peel_ms_mu": "0.100",
        "project_ms_mu": "0.200",
        "residual_ms_mu": "0.300",
        "backsub_ms_mu": "0.400",
        "seed_attempt": "" if candidate else "0",
        "block_xors_mu": "98.000" if candidate else "100.000",
        "block_muladds_mu": "50.000",
        "first_rank_fail": "0" if weak else "-1",
        "binary_def_hist": "0:1",
        "heavy_gain_hist": "0:1",
        "failure_trials": "0" if weak else "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(attempt),
        "packet_attempt": str(attempt),
        "attempt_mode": "exact" if candidate else "selected",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": subject.ZERO_SHA256,
        "effective_precode_seed": effective_precode,
        "effective_packet_seed": effective_packet,
    }
    lines = [metadata(invocation), ",".join(subject.CSV_HEADER),
             ",".join(values[key] for key in subject.CSV_HEADER)]
    return ("\n".join(lines) + "\n").encode("utf-8")


def synthetic_records(candidate_xors: int = 98) -> list:
    attempts = {K: 0 for K in subject.K_VALUES}
    invocations = subject.make_invocations(attempts)
    records = []
    for invocation in invocations:
        candidate = invocation.arm == "candidate_two07_mix2"
        staircase = 346 if invocation.K == 64000 else 2
        entries = subject.extra_dense_basis_capacity_entries(
            invocation.K, staircase)
        candidate_attempt = invocation.attempt if candidate else 0
        records.append({
            "schema": subject.RESULT_SCHEMA,
            "ordinal": invocation.ordinal,
            "coordinate_ordinal": invocation.coordinate_ordinal,
            "arm": invocation.arm,
            "K": invocation.K,
            "block_bytes": subject.BLOCK_BYTES,
            "loss_ppm": subject.LOSS_PPM,
            "overhead": subject.OVERHEAD,
            "root_index": invocation.root_index,
            "loss_root": invocation.root,
            "schedule": invocation.schedule,
            "cell_ordinal": invocation.cell_ordinal,
            "timing_order": invocation.timing_order,
            "timing_slot": invocation.timing_slot,
            "observation_index": invocation.observation_index,
            "attempt_selection_stream_sha256": "e" * 64,
            "attempt_selection_cell_receipts_used_as_promotion_evidence":
                False,
            "benchmark_loss_trace_hash_recorded": False,
            "full_payload_byte_recovery_verified": True,
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(),
            "construction_attempt": candidate_attempt,
            "attempt_mode": "exact" if candidate else "selected",
            "effective_precode_seed":
                "0x{:016x}".format(0x1000 + invocation.K +
                                     candidate_attempt),
            "effective_packet_seed":
                "0x{:08x}".format(0x1000 + invocation.K +
                                   candidate_attempt),
            "staircase": staircase,
            "binary_dense_rows": 12,
            "gf256_heavy_rows": 12,
            "source_hits": 2,
            "dense_anchor_layout": "two07" if candidate else "disabled",
            "mix_count": 2 if candidate else 3,
            "success": 1,
            "rank_fail": 0,
            "error": 0,
            "weak": False,
            "block_xors": candidate_xors if candidate else 100,
            "gf256_muladds": 50,
            "inactivated_columns": 10,
            "solve_ms": "1.000",
            "build_ms": "1.000",
            "extra_dense_basis_capacity_entries": entries if candidate else 0,
            "extra_dense_basis_capacity_bytes": entries * 4 if candidate else 0,
            "command_sha256": subject.sha256_text(subject.canonical_json(
                invocation.argv(Path("bench")))),
            "bench_binary_sha256": BENCH_SHA,
            "bench_source_git_commit": SOURCE_COMMIT,
            "source_receipt_sha256": SOURCE_SHA,
        })
    return records


class ContractTests(unittest.TestCase):
    def test_preregistered_domain_and_hash_are_stable(self) -> None:
        expected_K = (
            2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512,
            513, 1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000,
            10001, 16384, 20000, 20001, 32768, 49152, 64000,
        )
        expected_selection_roots = (
            "0xd1b54a32d192ed03",
            "0x94d049bb133111eb",
            "0x8538ecb5bd456ea3",
            "0xc0ac29b7c97c50dd",
            "0x3f84d5b5b5470917",
            "0x9216d5d98979fb1b",
        )
        expected_roots = (
            "0xb889883a79549774",
            "0xb5666de0987896af",
            "0x8bfca269b0bc01e0",
        )
        self.assertEqual(subject.K_VALUES, expected_K)
        self.assertEqual(subject.ROOTS, expected_roots)
        self.assertEqual(subject.SELECTION_ROOTS, expected_selection_roots)
        self.assertEqual(subject.SCHEDULES,
                         ("burst", "adversarial", "repair-only"))
        self.assertEqual(len(subject.K_VALUES), 30)
        self.assertEqual(subject.EXPECTED_CELL_COUNT, 270)
        self.assertEqual(subject.EXPECTED_RECORD_COUNT, 1080)
        self.assertEqual(subject.K_VALUES[0], 2)
        self.assertEqual(subject.K_VALUES[-1], 64000)
        self.assertEqual(subject.BLOCK_BYTES, 2)
        for index, (root, digest) in enumerate(zip(
                subject.ROOTS, subject.BENCHMARK_ROOT_FULL_SHA256)):
            preimage = subject.BENCHMARK_ROOT_NAMESPACE_TEMPLATE.format(
                i=index).encode("ascii")
            self.assertEqual(hashlib.sha256(preimage).hexdigest(), digest)
            self.assertEqual(root, "0x" + digest[:16])
        frozen_root_sets = (
            set(subject.SELECTION_ROOTS), set(subject.ROOTS),
            set(subject.DISCARDED_V6_ROOTS),
            set(subject.RETIRED_V7_ROOTS),
            set(subject.FINAL_VALIDATION_ROOTS),
        )
        for left in range(len(frozen_root_sets)):
            for right in range(left + 1, len(frozen_root_sets)):
                self.assertFalse(
                    frozen_root_sets[left] & frozen_root_sets[right])
        contract = subject.contract_description()
        self.assertEqual(
            contract["contract_sha256"],
            subject.self_hash(contract, "contract_sha256"))
        self.assertEqual(
            contract["gates"]["block_xor_ratio_max"], "0.9829453304")
        self.assertEqual(contract["gates"], {
            "candidate_benchmark_weak_observations_max": 0,
            "candidate_benchmark_weak_observations_not_above_control": True,
            "repeated_deterministic_work_equal": True,
            "block_xor_ratio_max": "0.9829453304",
            "gf256_muladd_ratio_max": "1",
            "inactivation_ratio_max": "1",
            "solve_time_ratio_max": "1",
            "build_time_ratio_max": "1.05",
        })
        self.assertEqual(
            contract["arm_definitions"]["current_disabled_mix3"], {
                "configuration": "exact current production selection",
                "dense_anchor_layout": "disabled",
                "mix_count": 3,
                "construction_attempt": "SelectSystematicConfiguration",
            })
        self.assertEqual(
            contract["arm_definitions"]["candidate_two07_mix2"], {
                "configuration": "normalized production candidate",
                "graph_seed_basis": "SelectSeedProfile(K,2)",
                "dense_anchor_layout": "two07",
                "mix_count": 2,
                "construction_attempt":
                    "exact v8 18-cell-selected uint8 attempt indexed by K",
            })
        self.assertEqual(
            contract["complexity"], {
                "allocation_count_delta": 0,
                "asymptotic_construction": "O(K)",
                "asymptotic_storage": "O(K)",
                "production_map_lookups_per_construction": 1,
                "online_retries": 0,
                "extra_dense_basis_capacity_entries":
                    "ceil((K+S+12)/2)-2 uint32 entries",
            })
        semantics = contract["evidence_semantics"]
        self.assertIn("selection provenance only", semantics[
            "attempt_selection_provenance"])
        self.assertIn("no packet IDs or trace hash",
                      semantics["benchmark_loss_trace"])
        self.assertIn("not invoked", semantics["final_validation_holdout"])
        self.assertEqual(semantics["benchmark_identity_preflight"], {
            "command": "--describe",
            "schema": subject.BENCH_DESCRIPTION_SCHEMA,
            "rule":
                "the pinned benchmark must report the exact frozen source "
                "commit before worker description, any D command, or any "
                "holdout-root invocation",
        })
        self.assertEqual(contract["timing_protocol"]["orders"],
                         ["ABBA", "BAAB"])
        self.assertFalse(
            contract["timing_protocol"]["cpu_affinity_required"])
        root_derivation = contract["benchmark_root_derivation"]
        self.assertEqual(root_derivation["full_sha256"],
                         list(subject.BENCHMARK_ROOT_FULL_SHA256))
        self.assertTrue(root_derivation["discarded_v6_roots_excluded"])
        self.assertTrue(root_derivation["retired_v7_roots_excluded"])

    def test_candidate_profile_and_description_are_exact(self) -> None:
        expected = hashlib.sha256(subject.canonical_json(
            subject.CANDIDATE_PROFILE).encode("utf-8")).hexdigest()
        self.assertEqual(subject.candidate_profile_sha256(), expected)
        actual = subject._validate_description(
            worker_description(), WORKER_SHA, SOURCE_COMMIT)
        self.assertEqual(actual["candidate_profile"]["mix_count"], 2)
        self.assertEqual(
            subject.validation_roster_sha256(),
            subject.VALIDATION_ROSTER_SHA256)
        broken = copy.deepcopy(worker_description())
        broken["candidate_profile"]["graph_seed_block_bytes"] = 17
        with self.assertRaises(subject.ShortScreenError):
            subject._validate_description(broken, WORKER_SHA, SOURCE_COMMIT)
        broken = copy.deepcopy(worker_description())
        broken["candidate_profile"]["dense_identity_corner"] = 0
        with self.assertRaises(subject.ShortScreenError):
            subject._validate_description(broken, WORKER_SHA, SOURCE_COMMIT)
        broken = copy.deepcopy(worker_description())
        broken["validation_roster_sha256"] = "f" * 64
        with self.assertRaises(subject.ShortScreenError):
            subject._validate_description(broken, WORKER_SHA, SOURCE_COMMIT)

    def test_invocations_use_selected_control_and_exact_candidate(self) -> None:
        invocations = subject.make_invocations(
            {K: (K % 3) for K in subject.K_VALUES})
        control, candidate = invocations[:2]
        self.assertEqual(control.arm, "current_disabled_mix3")
        self.assertNotIn("--exact-precode-attempt", control.argv(Path("bench")))
        self.assertEqual(candidate.arm, "candidate_two07_mix2")
        command = candidate.argv(Path("bench"))
        self.assertEqual(command[command.index("--mix-count") + 1], "2")
        self.assertEqual(command[command.index("--dense-anchors") + 1],
                         "two07")
        self.assertEqual(
            command[command.index("--exact-precode-attempt") + 1], "2")
        self.assertEqual(
            command[command.index("--construction-seed-basis") + 1],
            "production-profile")
        self.assertEqual(
            [item.arm for item in invocations[:4]],
            ["current_disabled_mix3", "candidate_two07_mix2",
             "candidate_two07_mix2", "current_disabled_mix3"])
        self.assertEqual(
            [item.observation_index for item in invocations[:4]],
            [0, 0, 1, 1])
        self.assertEqual(
            [item.arm for item in invocations[4:8]],
            ["candidate_two07_mix2", "current_disabled_mix3",
             "current_disabled_mix3", "candidate_two07_mix2"])


class ExecutablePinningTests(unittest.TestCase):
    @staticmethod
    def _write_executable(path: Path, marker: str) -> None:
        path.write_text(
            "#!/usr/bin/env python3\nprint({!r})\n".format(marker),
            encoding="ascii")
        path.chmod(0o755)

    def test_path_replacement_cannot_substitute_pinned_executable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "measured"
            replacement = root / "replacement"
            self._write_executable(executable, "pinned-A")
            self._write_executable(replacement, "replacement-B")
            with subject._open_binary(executable, "test") as pinned:
                original_fd = pinned.descriptor
                original_receipt = pinned.receipt()
                os.replace(replacement, executable)
                stdout, stderr = subject._run_process(
                    [str(pinned.path)], None, 5, 1024,
                    "pinned executable test", pinned.descriptor)
                self.assertEqual(stdout, b"pinned-A\n")
                self.assertEqual(stderr, b"")
                self.assertEqual(pinned.receipt(), original_receipt)
                self.assertNotEqual(
                    hashlib.sha256(executable.read_bytes()).hexdigest(),
                    original_receipt["sha256"])
            with self.assertRaises(OSError):
                os.fstat(original_fd)

    def test_pinned_descriptor_closes_on_exception(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "measured"
            self._write_executable(executable, "failure-path")
            pinned = subject._open_binary(executable, "test")
            original_fd = pinned.descriptor
            with self.assertRaisesRegex(RuntimeError, "injected"):
                with pinned:
                    raise RuntimeError("injected")
            self.assertEqual(pinned.descriptor, -1)
            with self.assertRaises(OSError):
                os.fstat(original_fd)


class BenchDescriptionTests(unittest.TestCase):
    def test_exact_canonical_description_uses_pinned_descriptor(self) -> None:
        stream = (subject.canonical_json(bench_description()) + "\n").encode(
            "ascii")
        with mock.patch.object(
                subject, "_run_process", return_value=(stream, b"")) as run:
            actual = subject.read_bench_description(
                Path("bench"), SOURCE_COMMIT, 17)
        self.assertEqual(actual, bench_description())
        self.assertEqual(run.call_args.args[0], ["bench", "--describe"])
        self.assertEqual(run.call_args.args[-1], 17)

    def test_description_rejects_identity_and_framing_drift(self) -> None:
        wrong_commit = bench_description()
        wrong_commit["source_git_commit"] = "2" * 40
        wrong_schema = bench_description()
        wrong_schema["schema"] = "wirehair.wh2.v2-bench-description.v0"
        extra = bench_description()
        extra["extra"] = 1
        valid = subject.canonical_json(bench_description()).encode("ascii")
        cases = (
            (subject.canonical_json(wrong_commit).encode("ascii") + b"\n",
             b""),
            (subject.canonical_json(wrong_schema).encode("ascii") + b"\n",
             b""),
            (subject.canonical_json(extra).encode("ascii") + b"\n", b""),
            (b'{"source_git_commit":"' + SOURCE_COMMIT.encode("ascii") +
             b'"}\n', b""),
            (b'{ "schema": "wirehair.wh2.v2-bench-description.v1", '
             b'"source_git_commit": "' + SOURCE_COMMIT.encode("ascii") +
             b'" }\n', b""),
            (valid, b""),
            (valid + b"\r\n", b""),
            (valid + b"\n" + valid + b"\n", b""),
            (valid + b"\n", b"unexpected"),
        )
        for stdout, stderr in cases:
            with self.subTest(stdout=stdout, stderr=stderr), \
                    mock.patch.object(
                        subject, "_run_process",
                        return_value=(stdout, stderr)), \
                    self.assertRaises(subject.ShortScreenError):
                subject.read_bench_description(
                    Path("bench"), SOURCE_COMMIT, 17)

    def test_stale_benchmark_stops_before_worker_or_holdout_work(self) -> None:
        pinned_bench = mock.Mock()
        pinned_bench.path = Path("/unused/bench")
        pinned_bench.descriptor = 11
        pinned_bench.receipt.return_value = {
            "path": "/unused/bench", "sha256": BENCH_SHA, "size": 1,
        }
        pinned_worker = mock.Mock()
        pinned_worker.path = Path("/unused/worker")
        pinned_worker.descriptor = 12
        pinned_worker.receipt.return_value = {
            "path": "/unused/worker", "sha256": WORKER_SHA, "size": 1,
        }
        source = {
            "source_git_commit": SOURCE_COMMIT,
            "source_receipt_sha256": SOURCE_SHA,
        }
        with mock.patch.object(subject, "_source_receipt",
                               return_value=source), \
                mock.patch.object(
                    subject, "read_bench_description",
                    side_effect=subject.ShortScreenError("stale")) as preflight, \
                mock.patch.object(
                    subject, "read_worker_description") as worker_description_call, \
                mock.patch.object(subject, "derive_attempts") as derive, \
                mock.patch.object(subject, "run_bench_cell") as bench_cell, \
                self.assertRaisesRegex(subject.ShortScreenError, "stale"):
            subject._run_screen_pinned(
                pinned_bench, pinned_worker, Path("/unused/output"),
                Path("/unused"), 0.0)
        preflight.assert_called_once_with(
            Path("/unused/bench"), SOURCE_COMMIT, 11)
        worker_description_call.assert_not_called()
        derive.assert_not_called()
        bench_cell.assert_not_called()


class WorkerProtocolTests(unittest.TestCase):
    def test_derivation_issues_exactly_30_D_commands_and_no_V(self) -> None:
        values = [derivation(K, 0) for K in subject.K_VALUES]
        stream = ("\n".join(subject.canonical_json(value)
                            for value in values) + "\n").encode("utf-8")
        with mock.patch.object(
                subject, "_run_process", return_value=(stream, b"")) as run:
            records, canonical, command_sha = subject.derive_attempts(
                Path("worker"), WORKER_SHA)
        self.assertEqual(len(records), 30)
        self.assertEqual(canonical, stream)
        command_bytes = run.call_args.args[1]
        commands = command_bytes.decode("ascii").splitlines()
        self.assertEqual(commands[-1], "Q")
        self.assertEqual(commands[:-1], [
            "D {} {}".format(K - 2, K) for K in subject.K_VALUES
        ])
        self.assertFalse(any(command.startswith("V ") for command in commands))
        self.assertEqual(command_sha,
                         hashlib.sha256(command_bytes).hexdigest())

    def test_complete_derivation_stream_is_canonical_and_bound(self) -> None:
        records = [derivation(K, K % 3) for K in subject.K_VALUES]
        stream = ("\n".join(subject.canonical_json(value)
                            for value in records) + "\n").encode("utf-8")
        parsed = subject.parse_derivation_stream(stream, WORKER_SHA)
        self.assertEqual(len(parsed), 30)
        self.assertEqual(parsed[-1]["K"], 64000)
        self.assertEqual(parsed[-1]["selected_attempt"], 1)

    def test_derivation_rejects_seed_schedule_drift(self) -> None:
        value = derivation(8, 2)
        value["effective_packet_seed"] = "0x00000000"
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_rejects_bool_identity_aliases(self) -> None:
        value = derivation(2, 0)
        value["K"] = 2.0
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 2, WORKER_SHA)
        value = derivation(3, 0)
        value["ordinal"] = True
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 3, WORKER_SHA)
        value = derivation(8, 0)
        value["selected_successes"][0]["construction_attempt"] = False
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)
        value = derivation(8, 0)
        value["selected_successes"][0]["loss_ppm"] = 500000.0
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_rejects_missing_selection_cell(self) -> None:
        value = derivation(8, 0)
        value["selected_successes"].pop()
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_rejects_lower_witness_trace_drift(self) -> None:
        value = derivation(8, 1)
        value["lower_attempt_failure_witnesses"][0]["trace_sha256"] = \
            "f" * 64
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_rejects_unemittable_failure_receipt(self) -> None:
        for outcome in ("success", "unsupported", "fatal", "invented"):
            value = derivation(8, 1)
            value["lower_attempt_failure_witnesses"][0]["outcome"] = outcome
            with self.subTest(outcome=outcome), self.assertRaises(
                    subject.ShortScreenError):
                subject.validate_derivation_record(value, 8, WORKER_SHA)
        value = derivation(8, 1)
        value["lower_attempt_failure_witnesses"][0]["decoded_extra"] = None
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)
        value = derivation(8, 1)
        value["lower_attempt_failure_witnesses"][0]["attempted_candidates"] = 0
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)
        for outcome, wrong_result in (("need_more_at_cap", 4),
                                      ("construct_failed", 1)):
            value = derivation(8, 1)
            witness = value["lower_attempt_failure_witnesses"][0]
            witness["outcome"] = outcome
            witness["decoded_extra"] = None
            witness["result_code"] = wrong_result
            with self.subTest(outcome=outcome), self.assertRaises(
                    subject.ShortScreenError):
                subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_accepts_worker_emittable_bounded_failures(self) -> None:
        for outcome, result_code in (("need_more_at_cap", 1),
                                     ("construct_failed", 4)):
            value = derivation(8, 1)
            witness = value["lower_attempt_failure_witnesses"][0]
            witness["outcome"] = outcome
            witness["decoded_extra"] = None
            witness["result_code"] = result_code
            with self.subTest(outcome=outcome):
                subject.validate_derivation_record(value, 8, WORKER_SHA)

    def test_derivation_rejects_noncanonical_order(self) -> None:
        value = derivation(8, 0)
        value["selected_successes"][0], value["selected_successes"][1] = \
            value["selected_successes"][1], value["selected_successes"][0]
        with self.assertRaises(subject.ShortScreenError):
            subject.validate_derivation_record(value, 8, WORKER_SHA)


class BenchProtocolTests(unittest.TestCase):
    def setUp(self) -> None:
        self.attempts = {K: 0 for K in subject.K_VALUES}
        self.control, self.candidate = subject.make_invocations(
            self.attempts)[:2]
        self.record = derivation(2, 0)

    def parse(self, invocation: subject.Invocation, payload: bytes) -> dict:
        return subject.parse_bench_stdout(
            payload, invocation, SOURCE_COMMIT, BENCH_SHA, self.record,
            subject.sha256_text(subject.canonical_json(
                invocation.argv(Path("bench")))), SOURCE_SHA, "e" * 64)

    def test_control_and_candidate_rows_parse(self) -> None:
        control = self.parse(
            self.control, bench_stdout(self.control, self.record))
        candidate = self.parse(
            self.candidate, bench_stdout(self.candidate, self.record))
        self.assertEqual(control["attempt_mode"], "selected")
        self.assertEqual(candidate["attempt_mode"], "exact")
        self.assertEqual(candidate["block_xors"], 98)
        self.assertEqual(
            candidate["extra_dense_basis_capacity_entries"], 6)
        self.assertEqual(
            candidate["attempt_selection_stream_sha256"], "e" * 64)
        self.assertFalse(control[
            "attempt_selection_cell_receipts_used_as_promotion_evidence"])
        self.assertFalse(candidate[
            "attempt_selection_cell_receipts_used_as_promotion_evidence"])
        self.assertFalse(control["benchmark_loss_trace_hash_recorded"])
        self.assertFalse(candidate["benchmark_loss_trace_hash_recorded"])
        self.assertTrue(candidate["full_payload_byte_recovery_verified"])

    def test_candidate_rejects_online_selection(self) -> None:
        payload = bench_stdout(self.candidate, self.record).decode("utf-8")
        payload = payload.replace("exact_attempt_mode=1",
                                  "exact_attempt_mode=0", 1)
        with self.assertRaises(subject.ShortScreenError):
            self.parse(self.candidate, payload.encode("utf-8"))

    def test_candidate_rejects_effective_seed_drift(self) -> None:
        payload = bench_stdout(self.candidate, self.record).decode("utf-8")
        payload = payload.replace(
            self.record["effective_packet_seed"], "0x00000000", 1)
        with self.assertRaises(subject.ShortScreenError):
            self.parse(self.candidate, payload.encode("utf-8"))

    def test_row_rejects_non_full_payload_metadata(self) -> None:
        payload = bench_stdout(self.control, self.record).decode("utf-8")
        payload = payload.replace("full_payload_solve=1",
                                  "full_payload_solve=0", 1)
        with self.assertRaises(subject.ShortScreenError):
            self.parse(self.control, payload.encode("utf-8"))

    def test_weak_cell_is_retained(self) -> None:
        parsed = self.parse(
            self.candidate,
            bench_stdout(self.candidate, self.record, weak=True))
        self.assertTrue(parsed["weak"])
        self.assertEqual(parsed["rank_fail"], 1)


class AdjudicationTests(unittest.TestCase):
    def test_all_frozen_gates_pass_at_bound(self) -> None:
        records = synthetic_records(candidate_xors=98)
        summary = subject.adjudicate(records)
        self.assertEqual(summary["disposition"], "PASS")
        self.assertEqual(summary["common_success_cells"], 270)
        self.assertEqual(summary["common_success_timing_pairs"], 540)
        self.assertEqual(summary["ratios"]["block_xors"],
                         "0.980000000000")
        self.assertEqual(summary["aggregates"]["control"]["solve_ms"],
                         "540.000")
        complexity = summary["complexity_receipt"]
        self.assertEqual(complexity["maximum_extra_entries"], 32177)
        self.assertEqual(complexity["maximum_extra_bytes"], 128708)

    def test_xor_gate_rejects_regression(self) -> None:
        summary = subject.adjudicate(synthetic_records(candidate_xors=99))
        self.assertEqual(summary["disposition"], "REJECT")
        self.assertFalse(
            summary["gates"]["block_xor_ratio_at_most_0.9829453304"])

    def test_candidate_weak_cell_rejects_and_is_excluded(self) -> None:
        records = synthetic_records()
        candidate = next(value for value in records
                         if value["arm"] == "candidate_two07_mix2")
        candidate_pair = [value for value in records
                          if value["arm"] == candidate["arm"] and
                          value["coordinate_ordinal"] ==
                          candidate["coordinate_ordinal"]]
        for value in candidate_pair:
            value["weak"] = True
            value["success"] = 0
            value["rank_fail"] = 1
            value["full_payload_byte_recovery_verified"] = False
        summary = subject.adjudicate(records)
        self.assertEqual(summary["disposition"], "REJECT")
        self.assertEqual(summary["candidate_weak_observations"], 2)
        self.assertEqual(summary["candidate_weak_coordinates"], 1)
        self.assertEqual(summary["common_success_cells"], 269)

    def test_control_weak_cell_is_excluded_from_common_sums(self) -> None:
        records = synthetic_records()
        control = next(value for value in records
                       if value["arm"] == "current_disabled_mix3")
        control_pair = [value for value in records
                        if value["arm"] == control["arm"] and
                        value["coordinate_ordinal"] ==
                        control["coordinate_ordinal"]]
        for value in control_pair:
            value["weak"] = True
            value["success"] = 0
            value["rank_fail"] = 1
            value["full_payload_byte_recovery_verified"] = False
            value["block_xors"] = 1
        summary = subject.adjudicate(records)
        self.assertEqual(summary["disposition"], "PASS")
        self.assertEqual(summary["control_weak_observations"], 2)
        self.assertEqual(summary["control_weak_coordinates"], 1)
        self.assertEqual(summary["common_success_cells"], 269)
        self.assertEqual(summary["ratios"]["block_xors"],
                         "0.980000000000")

    def test_capacity_tampering_fails_closed(self) -> None:
        records = synthetic_records()
        candidate = next(value for value in records
                         if value["arm"] == "candidate_two07_mix2")
        candidate["extra_dense_basis_capacity_bytes"] += 4
        with self.assertRaises(subject.ShortScreenError):
            subject.adjudicate(records)

    def test_nonzero_expected_capacity_is_a_comparison_not_a_operand(self) \
            -> None:
        records = synthetic_records()
        candidate = next(value for value in records
                         if value["arm"] == "candidate_two07_mix2" and
                         value["K"] == 2)
        self.assertEqual(candidate["extra_dense_basis_capacity_entries"], 6)
        self.assertEqual(subject.adjudicate(records)["disposition"], "PASS")

    def test_selection_receipts_cannot_be_promotion_evidence(self) -> None:
        records = synthetic_records()
        control = next(value for value in records
                       if value["arm"] == "current_disabled_mix3")
        control[
            "attempt_selection_cell_receipts_used_as_promotion_evidence"] = \
            True
        with self.assertRaises(subject.ShortScreenError):
            subject.adjudicate(records)

        records = synthetic_records()
        candidate = next(value for value in records
                         if value["arm"] == "candidate_two07_mix2")
        candidate["benchmark_loss_trace_hash_recorded"] = True
        with self.assertRaises(subject.ShortScreenError):
            subject.adjudicate(records)

    def test_repeated_work_must_be_equal(self) -> None:
        records = synthetic_records()
        candidate_pair = [value for value in records
                          if value["arm"] == "candidate_two07_mix2" and
                          value["coordinate_ordinal"] == 0]
        candidate_pair[1]["block_xors"] += 1
        with self.assertRaises(subject.ShortScreenError):
            subject.adjudicate(records)

    def test_result_stream_must_bind_one_attempt_stream(self) -> None:
        records = synthetic_records()
        records[-1]["attempt_selection_stream_sha256"] = "f" * 64
        with self.assertRaises(subject.ShortScreenError):
            subject.adjudicate(records, "e" * 64)

    def test_exact_capacity_formula_rounds_up(self) -> None:
        self.assertEqual(subject.extra_dense_basis_capacity_entries(2, 2), 6)
        self.assertEqual(
            subject.extra_dense_basis_capacity_entries(64000, 346), 32177)
        self.assertEqual(32177 * 4, 128708)


if __name__ == "__main__":
    unittest.main()
