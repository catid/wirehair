#!/usr/bin/env python3
"""Bounded synthetic/tamper tests for wh2_h13_stage_b.py."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_h13_stage_b as campaign  # noqa: E402


def native_benchmark_path() -> Path:
    configured = os.environ.get("WIREHAIR_V2_BENCH")
    if configured is not None:
        path = Path(configured)
        if not path.is_file():
            raise AssertionError(
                f"WIREHAIR_V2_BENCH is not a file: {configured}"
            )
        return path
    return (Path(__file__).resolve().parents[1] / "build-h13-stage-a" /
            "codec" / "wirehair_v2_bench")


def freeze_executable_copy(source: Path, destination: Path) -> None:
    destination.write_bytes(source.read_bytes())
    destination.chmod(0o555)


def source_arms(
    key: tuple[int, str, int], h12_fail: bool = False,
    h13_fail: bool = False, h12_probe: int = 0, h13_probe: int = 0,
) -> dict[str, dict[str, str]]:
    seed_index, schedule, K = key
    del schedule
    base_matrix = hex(0x1000000000000000 + K)
    base_peel = hex(0x10000000 + K)
    trace_seed = hex(campaign.loss_seed(campaign.SEEDS[seed_index], K))
    trace_sha = hashlib.sha256(
        f"synthetic-trace|{seed_index}|{key[1]}|{K}".encode("ascii")
    ).hexdigest()

    def row(failed: bool, probe: int) -> dict[str, str]:
        return {
            "success": "0" if failed else "1",
            "rank_fail": "1" if failed else "0", "error": "0",
            "base_matrix_seed": base_matrix, "base_peel_seed": base_peel,
            "matrix_seed": base_matrix, "peel_seed": base_peel,
            "actual_staircase_rows": "3", "actual_dense_rows": "12",
            "actual_source_hits": "2" if K < 10000 else "3",
            "actual_dense_two_anchor": "1",
            "actual_dense_two_anchor_phase": "0",
            "packet_trace_seed": trace_seed,
            "packet_trace_sha256": trace_sha,
            "systematic_probe_result": str(probe),
        }

    return {"h12": row(h12_fail, h12_probe),
            "h13": row(h13_fail, h13_probe)}


def cell_expectation(
    key: tuple[int, str, int], role: str = "union", match_id: int = 0,
    **source_options: object,
) -> dict[str, object]:
    return {
        "role": role, "match_id": match_id,
        "source": campaign.source_cell_record(
            key, source_arms(key, **source_options),
        ),
    }


def grouped_output(
    task: campaign.Task,
    cells: dict[tuple[int, str, int], dict[str, object]],
    results: dict[int, tuple[int, int]],
    probes: tuple[int, int] = (0, 0),
) -> str:
    lines = [campaign.expected_preamble_line(task),
             ",".join(campaign.BENCH_HEADER)]
    for pair_index, K in enumerate(task.ks):
        expected = cells[(task.seed_index, task.schedule, K)]
        source = expected["source"]
        arm_results = dict(zip(campaign.ARMS, results[K]))
        pair_class = (
            "both-success" if arm_results == {"h12": 0, "h13": 0} else
            "h12-only" if arm_results["h12"] == 0 else
            "h13-only" if arm_results["h13"] == 0 else "both-need-more"
        )
        equal = int(arm_results["h12"] == arm_results["h13"])
        trace_seed = int(source["packet_trace_seed"], 16)
        pair_id = campaign.grouped_pair_id(
            task, K, trace_seed, source["packet_trace_sha256"],
        )
        for order, arm in enumerate(campaign.ARMS):
            result = arm_results[arm]
            heavy = campaign.ARM_HEAVY_ROWS[arm]
            precode = source["actual_staircase_rows"] + 12 + heavy
            inactive = 6
            residual_rank = inactive if result == 0 else inactive - 1
            binary_rank = 3
            values = {name: "0" for name in campaign.BENCH_HEADER}
            values.update({
                "pair_id": pair_id, "pair_index": str(pair_index),
                "pair_order_index": str(order), "arm": arm, "N": str(K),
                "bb": "64", "solve_block_bytes": "2", "overhead": "0",
                "schedule": task.schedule,
                "external_seed": str(int(task.seed, 16)), "loss": "0.5",
                "period": "48", "geometry": "shared-x",
                "residue_skew": "0", "residue_schedule": "constant",
                "gf256_rows": str(campaign.ARM_GF256_ROWS[arm]),
                "gf16_rows": "2", "heavy_rows": str(heavy),
                "grouped_rows": "3", "grouped_gf256_row_mask": "0x380",
                "buckets_requested": "separate",
                "grouped_hash_seed": campaign.ARM_GROUP_HASH_SEED[arm],
                "final_h_a_columns": str(heavy), "construction_attempt": "0",
                "systematic_probe_result": str(probes[order]),
                "base_matrix_seed": source["base_matrix_seed"],
                "base_peel_seed": source["base_peel_seed"],
                "matrix_seed": source["matrix_seed"],
                "peel_seed": source["peel_seed"],
                "staircase_rows": str(source["actual_staircase_rows"]),
                "dense_rows": "12",
                "source_hits": str(source["actual_source_hits"]),
                "dense_identity_corner": "0", "dense_two_anchor": "1",
                "dense_two_anchor_phase": "0", "field": "mixed-gf256-gf16",
                "heavy_family": "periodic-cauchy", "mix_count": "2",
                "precode_count": str(precode), "packet_count": str(K),
                "packet_trace_seed": source["packet_trace_seed"],
                "packet_trace_sha256": source["packet_trace_sha256"],
                "pair_class": pair_class, "pair_results_equal": str(equal),
                "result": str(result),
                "result_name": "success" if result == 0 else "need-more",
                "packet_rows": str(K),
                "peeled_columns": str(K + precode - inactive),
                "inactivated_columns": str(inactive),
                "residual_rows": str(inactive),
                "binary_residual_rank": str(binary_rank),
                "residual_rank": str(residual_rank),
                "binary_deficit": str(inactive - binary_rank),
                "heavy_gain": str(residual_rank - binary_rank),
                "binary_row_references": "10",
                "binary_row_storage_bytes": "20",
                "binary_adjacency_storage_bytes": "30",
                "binary_row_storage_allocations": "1",
                "binary_adjacency_storage_allocations": "1",
                "block_xors": "40", "block_muladds": "50",
                "build_ns": "100", "peel_ns": "200", "project_ns": "300",
                "residual_ns": "400", "backsub_ns": "500",
                "packet_seed_attempt": "0", "joint_source_xors": "0",
                "joint_marginal_xors": "0", "joint_marginal_copies": "0",
                "joint_active_deltas": "0", "joint_scratch_bytes": "0",
                "dual_source_columns": "0",
                "intermediate_bytes": str(2 * (K + precode) if result == 0 else 0),
            })
            lines.append(",".join(values[name] for name in campaign.BENCH_HEADER))
    return "\n".join(lines) + "\n"


def mutate_csv(output: str, output_row: int, field: str, value: str) -> str:
    lines = output.splitlines()
    header = next(csv.reader([lines[1]]))
    values = next(csv.reader([lines[2 + output_row]]))
    values[header.index(field)] = value
    lines[2 + output_row] = ",".join(values)
    return "\n".join(lines) + "\n"


def metric_arms(h12_result: int, h13_result: int) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for arm, outcome in zip(campaign.ARMS, (h12_result, h13_result)):
        row = {name: "0" for name in campaign.WORK_RECEIPT_FIELDS}
        row["result"] = str(outcome)
        row["result_name"] = "success" if outcome == 0 else "need-more"
        row["pair_id"] = "a" * 64
        result[arm] = row
    return result


class SealTests(unittest.TestCase):
    def test_seal_canonical_round_trip_and_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "record.json"
            campaign.write_sealed_once(path, "test.schema", {"b": 2, "a": 1})
            self.assertEqual(campaign.load_sealed(path, "test.schema"),
                             {"a": 1, "b": 2})
            data = bytearray(path.read_bytes())
            data[data.index(b"2")] = ord("3")
            path.write_bytes(data)
            with self.assertRaises(campaign.CampaignError):
                campaign.load_sealed(path, "test.schema")

    def test_write_once_refuses_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "record.json"
            campaign.write_sealed_once(path, "test.schema", {"value": 1})
            with self.assertRaises(campaign.CampaignError):
                campaign.write_sealed_once(path, "test.schema", {"value": 2})

    def test_wrong_formal_seal_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "formal.json"
            campaign.write_sealed_once(path, "test.formal", {"ok": True})
            data = path.read_bytes()
            record = json.loads(data)
            file_sha = hashlib.sha256(data).hexdigest()
            self.assertEqual(campaign.require_exact_seal(
                path, "test.formal", record["seal_sha256"], file_sha,
            ), {"ok": True})
            with self.assertRaisesRegex(campaign.CampaignError, "file SHA256"):
                campaign.require_exact_seal(
                    path, "test.formal", record["seal_sha256"], "0" * 64,
                )
            with self.assertRaisesRegex(campaign.CampaignError, "pinned formal"):
                campaign.require_exact_seal(
                    path, "test.formal", "0" * 64, file_sha,
                )


class SourceTrustTests(unittest.TestCase):
    def test_forged_controller_rejected_before_import(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            frozen = root / "frozen"
            frozen.mkdir()
            marker = root / "executed"
            controller = frozen / "evil.py"
            controller.write_text(
                f"from pathlib import Path\nPath({str(marker)!r}).write_text('bad')\n",
                encoding="utf-8",
            )
            binary = frozen / "binary"
            taskset = frozen / "taskset"
            binary.write_bytes(b"binary")
            taskset.write_bytes(b"taskset")
            payload = {
                "result_dir": str(root), "binary": str(binary),
                "binary_sha256": campaign.sha256_file(binary),
                "controller": str(controller),
                "controller_sha256": campaign.sha256_file(controller),
                "taskset": str(taskset),
                "taskset_sha256": campaign.sha256_file(taskset),
            }
            campaign.write_sealed_once(
                root / "contract.json", f"{campaign.SOURCE_SCHEMA}.contract",
                payload,
            )
            with self.assertRaisesRegex(campaign.CampaignError, "pinned formal build"):
                campaign.authenticate_stage_a(root)
            self.assertFalse(marker.exists())

    def test_formal_identity_constants(self) -> None:
        self.assertIs(campaign.FORMAL_STAGE_A_LAUNCHABLE, False)
        self.assertIn("per-K construction-seed fixups",
                      campaign.FORMAL_STAGE_A_BLOCKER)
        with self.assertRaisesRegex(campaign.CampaignError, "launch is blocked"):
            campaign.require_formal_stage_a_launchable()
        self.assertEqual(campaign.FORMAL_STAGE_A_CONTROLLER_SHA256,
                         "f989dd991de3159abf3cb0a8c42deb118c41b2a868c4759c599d38170e6fa823")
        self.assertEqual(campaign.FORMAL_STAGE_A_BINARY_SHA256,
                         "df8da0d908e172f3cfd2ecba367b39d95fa253370f681560d1dd927e39b66a06")
        self.assertEqual(campaign.FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
                         "38c3bcbfa392a96a95a08342a02c46d0404a557aff44401c7a0cba37b2237285")
        self.assertEqual(campaign.FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
                         "3dee34a593203c6009d7f5ea3cbb913900eccbeec5ff40a0236d242f1f542cf0")
        self.assertEqual(campaign.FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
                         "61619d0539f72eea3b5f7045edbd11c2b78c3941bcb629bc5e9b092cc158d1a3")
        self.assertEqual(campaign.FORMAL_STAGE_A_ANALYSIS_SEAL,
                         "ead9727fd0f1684fc9179fc481b5987d2eafe7e896b8eba508d262a4a24c5c34")
        formal = Path("/dev/shm/wh2-h13-stage-a-formal-v1")
        if formal.is_dir():
            campaign.require_exact_seal(
                formal / "campaign_complete.json",
                f"{campaign.SOURCE_SCHEMA}.campaign_complete",
                campaign.FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
                campaign.FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
            )
            campaign.require_exact_seal(
                formal / "analysis.json",
                f"{campaign.SOURCE_SCHEMA}.analysis_record",
                campaign.FORMAL_STAGE_A_ANALYSIS_SEAL,
                campaign.FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
            )

    def test_formal_manifest_maps_every_allk_task_to_adjacent_arm_pair(self) -> None:
        formal = Path("/dev/shm/wh2-h13-stage-a-formal-v1")
        if not formal.is_dir():
            self.skipTest("formal Stage-A campaign is unavailable")
        contract = campaign.load_sealed(
            formal / "contract.json", f"{campaign.SOURCE_SCHEMA}.contract",
        )
        module = campaign._load_stage_a_module(
            Path(contract["controller"]),
            campaign.FORMAL_STAGE_A_CONTROLLER_SHA256,
        )
        source_tasks = module.load_manifest(module.stage_path(formal, 0), 0)
        allk_tasks = campaign.build_allk_tasks()
        self.assertEqual(len(source_tasks), 2 * len(allk_tasks))
        for task, left, right in zip(
                allk_tasks, source_tasks[::2], source_tasks[1::2]):
            self.assertEqual((left.arm, right.arm), campaign.ARMS)
            self.assertEqual((left.pair, right.pair), (task.job, task.job))
            self.assertEqual((left.seed_index, right.seed_index),
                             (task.seed_index, task.seed_index))
            self.assertEqual((left.seed, right.seed), (task.seed, task.seed))
            self.assertEqual((left.schedule, right.schedule),
                             (task.schedule, task.schedule))
            self.assertEqual((tuple(left.ks), tuple(right.ks)),
                             (task.ks, task.ks))

    def test_authenticated_module_executes_exact_hashed_buffer(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            marker = root / "substituted-code-ran"
            controller = root / "stage_a.py"
            controller.write_text(
                f"from pathlib import Path\nPath({str(marker)!r}).write_text('bad')\n",
                encoding="utf-8",
            )
            safe = (
                f"SCHEMA={campaign.SOURCE_SCHEMA!r}\n"
                f"K_MIN={campaign.K_MIN}\nK_MAX={campaign.K_MAX}\n"
                f"SEEDS={campaign.SEEDS!r}\nSCHEDULES={campaign.SCHEDULES!r}\n"
            ).encode("ascii")
            digest = hashlib.sha256(safe).hexdigest()
            with mock.patch.object(campaign, "stable_bytes", return_value=safe):
                module = campaign._load_stage_a_module(controller, digest)
            self.assertEqual(module.SCHEMA, campaign.SOURCE_SCHEMA)
            self.assertFalse(marker.exists())
            sys.modules.pop(f"_wirehair_stage_a_{digest}", None)

    def test_authentication_rejects_control_leaf_changed_before_return(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text).resolve()
            frozen = root / "frozen"
            frozen.mkdir()
            binary = frozen / "binary"
            controller = frozen / "stage_a.py"
            taskset = frozen / "taskset"
            for path in (binary, controller, taskset):
                path.write_bytes(b"fixture")
            independent = {
                "result_dir": str(root), "binary": str(binary),
                "binary_sha256": campaign.FORMAL_STAGE_A_BINARY_SHA256,
                "controller": str(controller),
                "controller_sha256": campaign.FORMAL_STAGE_A_CONTROLLER_SHA256,
                "taskset": str(taskset), "taskset_sha256": "6" * 64,
            }
            stage = root / "stages" / "oh0000"
            control_hashes = {
                root / "contract.json": "1" * 64,
                root / "controller_affinity.json": "2" * 64,
                stage / "manifest.json": "3" * 64,
                stage / "complete.json": "4" * 64,
                root / "terminal.json": "5" * 64,
            }
            completion = {
                "contract_sha256": "1" * 64,
                "controller_affinity_sha256": "2" * 64,
                "terminal_sha256": "5" * 64,
            }
            metrics = {"paired_cells": 0}
            sealed_complete = {
                "overhead": 0, "job_count": 0, "paired_cell_count": 0,
                "metrics": metrics, "union_failures": [],
            }
            analysis = {
                "schema": f"{campaign.SOURCE_SCHEMA}.analysis",
                "campaign_completion": completion,
                "coverage": {
                    "K_min": campaign.K_MIN, "K_max": campaign.K_MAX,
                    "unique_K": campaign.K_COUNT,
                    "packet_trace_roots": len(campaign.SEEDS),
                    "schedules": len(campaign.SCHEDULES),
                    "paired_cells": 0, "outcomes_per_arm": 0,
                    "arm_outcomes_at_OH0": 0, "terminal_overhead": 0,
                },
                "oh0": metrics,
            }
            mutated = {"value": False}

            class FakeDeriver:
                def __init__(self, _union: object) -> None:
                    pass

                def observe(self, _key: object, _arms: object) -> None:
                    pass

                def finish(self) -> dict:
                    mutated["value"] = True
                    return {
                        "source_union_count": 0, "matched_control_count": 0,
                        "paired_cells": 0, "arm_outcomes": 0,
                        "nonempty_strata": 0, "zero_strata": [],
                        "source_identity_stream_sha256": "8" * 64,
                    }

            module = mock.Mock()
            module.load_contract.return_value = independent
            module.stage_path.return_value = stage
            module.verify_terminal_campaign.return_value = {"terminal_overhead": 0}
            module.verify_campaign_completion.return_value = completion
            module.load_manifest.return_value = []
            module.load_sealed.side_effect = lambda path, _schema: (
                analysis if Path(path) == root / "analysis.json" else
                sealed_complete
            )

            def aggregate(stream: object) -> tuple:
                list(stream)
                return metrics, []

            module.aggregate_pair_stream.side_effect = aggregate

            def exact_seal(path: Path, *_args: object) -> dict:
                return completion if Path(path).name == "campaign_complete.json" \
                    else analysis

            def fake_sha(path: Path) -> str:
                path = Path(path)
                if path == binary:
                    return campaign.FORMAL_STAGE_A_BINARY_SHA256
                if path == controller:
                    return campaign.FORMAL_STAGE_A_CONTROLLER_SHA256
                if path == taskset:
                    return "6" * 64
                value = control_hashes[path]
                if path == root / "terminal.json" and mutated["value"]:
                    return "9" * 64
                return value

            with mock.patch.object(campaign, "load_sealed",
                                   return_value=independent), \
                    mock.patch.object(campaign, "sha256_file",
                                      side_effect=fake_sha), \
                    mock.patch.object(campaign, "require_exact_seal",
                                      side_effect=exact_seal), \
                    mock.patch.object(campaign, "_load_stage_a_module",
                                      return_value=module), \
                    mock.patch.object(campaign, "_source_affinity",
                                      return_value=[]), \
                    mock.patch.object(campaign, "_source_telemetry_receipts",
                                      return_value=(None, None)), \
                    mock.patch.object(campaign, "CohortDeriver", FakeDeriver), \
                    mock.patch.object(campaign, "build_screen_tasks",
                                      return_value=[]), \
                    mock.patch.multiple(
                        campaign, SOURCE_PAIRED_CELLS=0, SOURCE_ARM_OUTCOMES=0,
                        EXPECTED_SOURCE_UNION=0, EXPECTED_SCREEN_MATCHES=0,
                        EXPECTED_SCREEN_CELLS=0, EXPECTED_SCREEN_OUTCOMES=0,
                        EXPECTED_NONEMPTY_STRATA=0, EXPECTED_ZERO_STRATA=(),
                        EXPECTED_SCREEN_TASKS=0, ALLK_TASKS=0):
                with self.assertRaisesRegex(campaign.CampaignError,
                                            "changed during authenticated replay"):
                    campaign.authenticate_stage_a(root)


class SelectionTests(unittest.TestCase):
    def derive(self, insufficient: bool = False) -> tuple[dict[str, object], list]:
        keys = [(0, "burst", K) for K in range(2, 8)]
        union_keys = [keys[0], keys[2]]
        sealed_union = [campaign.cell_key_payload(key) for key in union_keys]
        deriver = campaign.CohortDeriver(
            sealed_union, expected_cells=len(keys), expected_union=2,
            require_formal_strata=False,
        )
        for key in keys:
            failed = key in union_keys
            if insufficient and key not in union_keys:
                failed = True
            deriver.observe(key, source_arms(key, failed, False))
        return deriver.finish(), keys

    def test_deterministic_hash_selected_controls(self) -> None:
        cohort, keys = self.derive()
        eligible = [key for key in keys if key[2] not in (2, 4)]
        expected = sorted(
            eligible,
            key=lambda key: campaign._control_digest(0, "burst", 0, key[2]),
        )[:2]
        actual = [
            (match["control"]["seed_index"], match["control"]["schedule"],
             match["control"]["K"]) for match in cohort["matches"]
        ]
        self.assertEqual(actual, expected)
        self.assertEqual(cohort["paired_cells"], 4)
        self.assertEqual(cohort["arm_outcomes"], 8)

    def test_control_is_common_success_and_disjoint(self) -> None:
        cohort, _ = self.derive()
        unions = {match["union"]["K"] for match in cohort["matches"]}
        controls = {match["control"]["K"] for match in cohort["matches"]}
        self.assertFalse(unions & controls)
        self.assertTrue(all(
            not any(match["control"]["raw_failures"].values())
            for match in cohort["matches"]
        ))

    def test_insufficient_controls_fails_closed(self) -> None:
        with self.assertRaisesRegex(campaign.CampaignError, "insufficient"):
            self.derive(insufficient=True)

    def test_noncanonical_source_stream_rejected(self) -> None:
        union = [campaign.cell_key_payload((0, "burst", 2))]
        deriver = campaign.CohortDeriver(
            union, expected_cells=2, expected_union=1,
            require_formal_strata=False,
        )
        deriver.observe((0, "burst", 3), source_arms((0, "burst", 3)))
        with self.assertRaisesRegex(campaign.CampaignError, "canonically ordered"):
            deriver.observe(
                (0, "burst", 2), source_arms((0, "burst", 2), True, False),
            )

    def test_union_outcome_mismatch_rejected(self) -> None:
        union = [campaign.cell_key_payload((0, "burst", 2))]
        deriver = campaign.CohortDeriver(
            union, expected_cells=1, expected_union=1,
            require_formal_strata=False,
        )
        with self.assertRaisesRegex(campaign.CampaignError, "sealed union"):
            deriver.observe((0, "burst", 2), source_arms((0, "burst", 2)))

    def test_formal_stratum_invariants_are_frozen(self) -> None:
        self.assertEqual(campaign.EXPECTED_SOURCE_UNION, 893)
        self.assertEqual(campaign.EXPECTED_NONEMPTY_STRATA, 61)
        self.assertEqual(campaign.EXPECTED_ZERO_STRATA, (
            (0, "repair-only", 2), (2, "burst", 2),
        ))


class PlannerTests(unittest.TestCase):
    def test_small_screen_plan_is_pair_native(self) -> None:
        cohort, _ = SelectionTests().derive()
        tasks = campaign.build_screen_tasks(cohort)
        self.assertEqual(len(tasks), 1)
        self.assertEqual(tasks[0].match_ids, (0, 1))
        self.assertEqual(len(tasks[0].ks), 4)
        self.assertEqual(tasks[0].ks, tuple(sorted(tasks[0].ks)))

    def test_exact_allk_plan_cardinality(self) -> None:
        tasks = campaign.build_allk_tasks()
        self.assertEqual(len(tasks), 2304)
        self.assertEqual(sum(len(task.ks) for task in tasks), 575991)
        self.assertEqual(2 * sum(len(task.ks) for task in tasks), 1151982)
        self.assertTrue(all(not task.match_ids for task in tasks))

    def test_screen_task_requires_two_disjoint_cells_per_match(self) -> None:
        malformed = campaign._new_task("screen", 0, 0, "burst", [0], [3])
        with self.assertRaisesRegex(campaign.CampaignError, "sealed Stage-B geometry"):
            campaign.validate_task(malformed)

    def test_interleaved_control_chunks_do_not_require_global_k_order(self) -> None:
        keys = [(0, "burst", K) for K in range(2, 20)]
        expected = {
            key: {"role": "union" if index % 2 == 0 else "control",
                  "match_id": index // 2, "source": {}}
            for index, key in enumerate(keys)
        }
        order = keys[10:] + keys[:10]
        pairs = [(key, metric_arms(0, 0), expected[key]) for key in order]
        metrics = campaign.aggregate_screen_pairs(pairs, set(keys))
        self.assertEqual(metrics["scopes"]["full_matched"]["paired_cells"], 18)


class ParserTests(unittest.TestCase):
    def setUp(self) -> None:
        self.task = campaign._new_task("allk", 0, 0, "burst", (), [3])
        key = (0, "burst", 3)
        self.cells = {key: cell_expectation(
            key, h12_probe=1, h13_probe=1,
        )}

    def parse(self, output: str) -> list:
        return campaign.parse_output(output, self.task, self.cells)

    def assert_tamper(self, output: str, row: int, field: str, value: str) -> None:
        with self.assertRaises(campaign.CampaignError):
            self.parse(mutate_csv(output, row, field, value))

    def test_live_native_golden_pair_id_and_full_preamble(self) -> None:
        trace = "aacbc866461e70d720c804377f24479e683d0bfd562f78ba9bf377cede6ea45a"
        self.assertEqual(
            campaign.grouped_pair_id(
                self.task, 3, int("dd02fc599574f77c", 16), trace,
            ),
            "acf9cdeb1a5ef14f2a25f5a24fb784bdc55ccc4bcf69442a557e8e510d4a6d90",
        )
        golden = (
            "# groupedrecovery: schema=v1 pair_order=h12,h13 arms_per_N=2 "
            "N_count=1 N=3 overhead=0 seed=15111065706836454659 schedule=burst "
            "policy=h12-h13-q0-grouped-v1 period=48 geometry=shared-x "
            "residue_skew=0 residue_schedule=constant residue_hash_seed=0x0 "
            "extension_residue_seed_xor=0x4e independent_extension_residues=0 "
            "gf256_rows=h12:10|h13:11 gf16_rows=2 grouped_rows=3 "
            "grouped_gf256_row_mask=0x380 buckets=separate "
            "grouped_hash_seed=h12:0xb7e15162|h13:0xb7e15163 "
            "grouped_final_h_a_columns=arm-heavy-rows dense_rows=12 "
            "dense_identity_corner=0 dense_two_anchor=1 dense_two_anchor_phase=0 "
            "source_hits=profile staircase_rows=profile-dense-count "
            "field=mixed-gf256-gf16 heavy_family=periodic-cauchy "
            "construction_attempt=0 systematic_probe=direct-attempt0 "
            "seed_repair=disabled mix=2 bb=64 seed_block_bytes=64 "
            "solve_block_bytes=2 loss=0.5 overhead_stream=paired "
            "packet_row_seed_multiplier=0x1 packet_row_seed_avalanche=0 "
            "odd_packet_peel_seed_xor=0x0 packet_vector=one-shared-per-N "
            "payload=shared-zero-v1 packet_payload_bytes=2 "
            "packet_trace_schema=wirehair-wh2-precodefail-raw-packet-trace-v1 "
            "coefficient_layout=materialized-before-solve "
            "grouped_schedule_prefix=materialized-before-solve"
        )
        self.assertEqual(campaign.expected_preamble_line(self.task), golden)
        self.assertEqual(len(campaign.BENCH_HEADER), 74)

    def test_contract_architecture_covers_every_fixed_preamble_semantic(self) -> None:
        architecture = campaign.architecture_contract()
        rendered = {
            "schema": architecture["native_schema"],
            "pair_order": ",".join(architecture["pair_order"]),
            "arms_per_N": str(architecture["arms_per_N"]),
            "gf256_rows": "|".join(
                f"{arm}:{architecture['gf256_rows'][arm]}"
                for arm in campaign.ARMS
            ),
            "grouped_hash_seed": "|".join(
                f"{arm}:{architecture['grouped_hash_seed'][arm]}"
                for arm in campaign.ARMS
            ),
            "dense_identity_corner": str(int(
                architecture["dense_identity_corner"])),
            "dense_two_anchor": str(int(architecture["dense_two_anchor"])),
            "independent_extension_residues": str(int(
                architecture["independent_extension_residues"])),
            "mix": str(architecture["mix_count"]),
        }
        aliases = {
            "buckets": "buckets", "grouped_final_h_a_columns":
                "grouped_final_h_a_columns",
        }
        for preamble_name, contract_name in aliases.items():
            rendered[preamble_name] = str(architecture[contract_name])
        direct = (
            "policy", "overhead", "period", "geometry", "residue_skew",
            "residue_schedule", "residue_hash_seed",
            "extension_residue_seed_xor", "gf16_rows", "grouped_rows",
            "grouped_gf256_row_mask", "dense_rows", "dense_two_anchor_phase",
            "source_hits", "staircase_rows", "field", "heavy_family",
            "construction_attempt", "systematic_probe", "seed_repair", "bb",
            "seed_block_bytes", "solve_block_bytes", "loss", "overhead_stream",
            "packet_row_seed_multiplier", "packet_row_seed_avalanche",
            "odd_packet_peel_seed_xor", "packet_vector", "payload",
            "packet_payload_bytes", "packet_trace_schema", "coefficient_layout",
            "grouped_schedule_prefix",
        )
        for name in direct:
            rendered[name] = str(architecture[name])
        preamble = dict(campaign.expected_preamble_items(self.task))
        for dynamic in ("N_count", "N", "seed", "schedule"):
            preamble.pop(dynamic)
        self.assertEqual(rendered, preamble)

    def test_live_native_k3_output_parses_against_stage_a_receipt(self) -> None:
        binary = native_benchmark_path()
        if not binary.is_file():
            self.skipTest("bounded native groupedrecovery build is unavailable")
        output = subprocess.check_output([
            str(binary), "groupedrecovery", "--N", "3", "--overhead", "0",
            "--seed", "15111065706836454659", "--schedule", "burst",
        ], text=True, timeout=10)
        expected = self.cells[(0, "burst", 3)]["source"]
        expected.update({
            "base_matrix_seed": "0x679d20a06742baae",
            "base_peel_seed": "0x11f13757",
            "matrix_seed": "0x679d20a06742baae",
            "peel_seed": "0x11f13757", "actual_staircase_rows": 3,
            "actual_dense_rows": 12, "actual_source_hits": 2,
            "actual_dense_two_anchor": 1,
            "actual_dense_two_anchor_phase": 0,
            "packet_trace_seed": "0xdd02fc599574f77c",
            "packet_trace_sha256":
                "aacbc866461e70d720c804377f24479e683d0bfd562f78ba9bf377cede6ea45a",
        })
        pairs = self.parse(output)
        self.assertEqual(pairs[0]["h12"]["pair_id"],
                         "acf9cdeb1a5ef14f2a25f5a24fb784bdc55ccc4bcf69442a557e8e510d4a6d90")

    def test_success_and_need_more_intermediate_rules(self) -> None:
        success = grouped_output(self.task, self.cells, {3: (0, 0)})
        self.assertEqual(self.parse(success)[0]["h13"]["result_name"], "success")
        need_more = grouped_output(self.task, self.cells, {3: (1, 1)})
        self.assertEqual(self.parse(need_more)[0]["h12"]["intermediate_bytes"], "0")
        self.assert_tamper(need_more, 0, "intermediate_bytes", "2")
        self.assert_tamper(success, 1, "intermediate_bytes", "0")

    def test_grouped_probe_need_not_equal_source_probe(self) -> None:
        output = grouped_output(self.task, self.cells, {3: (0, 0)}, probes=(0, 0))
        self.parse(output)
        self.assert_tamper(output, 0, "systematic_probe_result", "2")

    def test_separate_bucket_joint_counter_tamper(self) -> None:
        output = grouped_output(self.task, self.cells, {3: (0, 0)})
        for field in (
            "joint_source_xors", "joint_marginal_xors",
            "joint_marginal_copies", "joint_active_deltas",
            "joint_scratch_bytes", "dual_source_columns",
        ):
            self.assert_tamper(output, 0, field, "1")

    def test_exact_solve_accounting_tampers(self) -> None:
        success = grouped_output(self.task, self.cells, {3: (0, 0)})
        for field, value in (
            ("peeled_columns", "1"),
            ("residual_rows", "5"),
            ("heavy_gain", "13"),
            ("residual_rank", "5"),
        ):
            self.assert_tamper(success, 0, field, value)
        need_more = grouped_output(self.task, self.cells, {3: (1, 1)})
        self.assert_tamper(need_more, 0, "residual_rank", "6")

    def test_architecture_attempt_and_source_tampers(self) -> None:
        output = grouped_output(self.task, self.cells, {3: (0, 1)})
        for field, value in (
            ("grouped_gf256_row_mask", "0x180"),
            ("grouped_hash_seed", "0xb7e15164"),
            ("construction_attempt", "1"),
            ("base_matrix_seed", "0x1"),
            ("packet_trace_sha256", "f" * 64),
            ("packet_seed_attempt", "1"),
        ):
            self.assert_tamper(output, 0, field, value)

    def test_pair_order_class_id_and_header_tampers(self) -> None:
        output = grouped_output(self.task, self.cells, {3: (0, 1)})
        self.assert_tamper(output, 0, "arm", "h13")
        self.assert_tamper(output, 1, "pair_class", "both-success")
        self.assert_tamper(output, 0, "pair_id", "0" * 64)
        lines = output.splitlines()
        lines[1] = lines[1].replace("pair_id", "pair_digest", 1)
        with self.assertRaises(campaign.CampaignError):
            self.parse("\n".join(lines) + "\n")

    def test_preamble_and_line_termination_tamper(self) -> None:
        output = grouped_output(self.task, self.cells, {3: (0, 0)})
        with self.assertRaises(campaign.CampaignError):
            self.parse(output.replace("solve_block_bytes=2", "solve_block_bytes=64", 1))
        with self.assertRaises(campaign.CampaignError):
            self.parse(output.rstrip("\n"))


class DecisionTests(unittest.TestCase):
    def metrics(self, union: list[tuple[int, int]], control: list[tuple[int, int]]) -> dict:
        pairs = []
        expected = set()
        match_id = 0
        K = 2
        for role, outcomes in (("union", union), ("control", control)):
            for h12, h13 in outcomes:
                key = (0, "burst", K)
                expected.add(key)
                pairs.append((key, metric_arms(h12, h13), {
                    "role": role, "match_id": match_id, "source": {},
                }))
                K += 1
                match_id += 1
        return campaign.aggregate_screen_pairs(pairs, expected)

    def test_pass_requires_both_scopes_strictly_better(self) -> None:
        metrics = self.metrics([(1, 0), (1, 1)], [(0, 0), (0, 0)])
        with mock.patch.object(campaign, "EXPECTED_SOURCE_UNION", 2), \
                mock.patch.object(campaign, "EXPECTED_SCREEN_CELLS", 4):
            decision = campaign.screen_decision_payload(metrics)
        self.assertTrue(decision["screen_passed"])

    def test_union_tie_fails_even_if_full_improves(self) -> None:
        metrics = self.metrics([(1, 0), (0, 1)], [(1, 0), (0, 0)])
        with mock.patch.object(campaign, "EXPECTED_SOURCE_UNION", 2), \
                mock.patch.object(campaign, "EXPECTED_SCREEN_CELLS", 4):
            decision = campaign.screen_decision_payload(metrics)
        self.assertFalse(decision["screen_passed"])
        self.assertFalse(decision["union_scope"]["strictly_improved"])

    def test_introductions_are_reported_not_veto(self) -> None:
        metrics = self.metrics(
            [(1, 0), (1, 0), (0, 1)], [(0, 0), (0, 0), (0, 0)],
        )
        with mock.patch.object(campaign, "EXPECTED_SOURCE_UNION", 3), \
                mock.patch.object(campaign, "EXPECTED_SCREEN_CELLS", 6):
            decision = campaign.screen_decision_payload(metrics)
        self.assertTrue(decision["screen_passed"])
        self.assertEqual(decision["union_scope"]["introductions"], 1)

    def test_common_success_work_excludes_discordant_solve_paths(self) -> None:
        common = metric_arms(0, 0)
        discordant = metric_arms(1, 0)
        for name in campaign.WORK_RECEIPT_FIELDS:
            common["h12"][name] = "2"
            common["h13"][name] = "3"
            discordant["h12"][name] = "100"
            discordant["h13"][name] = "200"
        scope = campaign._empty_scope()
        campaign._update_scope(scope, common)
        campaign._update_scope(scope, discordant)
        campaign._finish_scope(scope)
        comparison = scope["comparison"]
        self.assertEqual(comparison["common_success"], 1)
        self.assertEqual(comparison["repairs"], 1)
        for name in campaign.WORK_RECEIPT_FIELDS:
            self.assertEqual(
                comparison["common_success_work"]["h12"][f"{name}_sum"], 2,
            )
            self.assertEqual(
                comparison["common_success_work"]["h13"][f"{name}_sum"], 3,
            )
            self.assertEqual(
                comparison["h13_over_h12_common_success_work_ratios"][name],
                "1.5",
            )


class InventoryTelemetryTests(unittest.TestCase):
    def test_orphan_staging_cleanup_is_manifest_bound(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            stage = Path(root_text) / "screen"
            shards = stage / "shards"
            shards.mkdir(parents=True)
            task = campaign._new_task("screen", 0, 0, "burst", [0], [2])
            orphan = shards / f".staging-{task.task_id}-999"
            orphan.mkdir()
            campaign.cleanup_orphan_screen_staging(stage, [task])
            self.assertFalse(orphan.exists())
            forged = shards / ".staging-000000000000000000000000-999"
            forged.mkdir()
            with self.assertRaisesRegex(campaign.CampaignError, "manifest task"):
                campaign.cleanup_orphan_screen_staging(stage, [task])

    def test_allk_plan_to_completed_inventory_small_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as root_text, \
                mock.patch.object(campaign, "K_MAX", 3), \
                mock.patch.object(campaign, "K_COUNT", 2), \
                mock.patch.object(campaign, "CHUNK_SIZE", 2), \
                mock.patch.object(campaign, "ALLK_TASKS", 9), \
                mock.patch.object(campaign, "ALLK_CELLS", 18):
            root = Path(root_text)
            stage = root / "stages" / "allk"
            tasks = campaign.build_allk_tasks()
            campaign.write_sealed_once(
                stage / "manifest.json", f"{campaign.SCHEMA}.stage_manifest",
                campaign.manifest_payload("allk", tasks),
            )
            self.assertEqual(len(campaign.validate_allk_inventory(
                root, {}, complete_required=False,
            )), 9)
            shards = stage / "shards"
            for task in tasks:
                (shards / task.stem).mkdir(parents=True)
            campaign.write_sealed_once(
                stage / "complete.json", f"{campaign.SCHEMA}.allk_complete",
                {"synthetic": True},
            )
            self.assertEqual(len(campaign.validate_allk_inventory(
                root, {}, complete_required=True,
            )), 9)

    def test_telemetry_requires_all_eight_dram_temperatures(self) -> None:
        values = [
            "2026-08-01T00:00:00Z", "100.0", "99.0", "5000.0", "70.0",
            *(["40.0"] * 8), "0", "100.0", "100.0", "100.0", "0", "0",
        ]
        row = campaign.parse_telemetry_row(
            (",".join(values) + "\n").encode("ascii"), "test",
        )
        self.assertEqual(len(row["dimm_temperatures"]), 8)
        values[5] = ""
        with self.assertRaisesRegex(campaign.CampaignError, "DIMM"):
            campaign.parse_telemetry_row(
                (",".join(values) + "\n").encode("ascii"), "test",
            )


class ShardAndProcessTests(unittest.TestCase):
    def test_authenticated_fd_launch_ignores_frozen_path_replacement(self) -> None:
        binary_source = native_benchmark_path()
        if not binary_source.is_file():
            self.skipTest("bounded native groupedrecovery build is unavailable")
        taskset_text = shutil.which("taskset")
        self.assertIsNotNone(taskset_text)
        false_source = Path(shutil.which("false") or "/bin/false")
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            binary = root / "wirehair_v2_bench"
            taskset = root / "taskset"
            freeze_executable_copy(binary_source, binary)
            freeze_executable_copy(Path(taskset_text), taskset)
            binary_sha = campaign.sha256_file(binary)
            taskset_sha = campaign.sha256_file(taskset)
            binary_fd = campaign.open_authenticated_executable(binary, binary_sha)
            taskset_fd = campaign.open_authenticated_executable(taskset, taskset_sha)
            try:
                replacement_binary = root / "replacement-binary"
                replacement_taskset = root / "replacement-taskset"
                freeze_executable_copy(false_source, replacement_binary)
                freeze_executable_copy(false_source, replacement_taskset)
                os.replace(binary, root / "authenticated-binary-inode")
                os.replace(taskset, root / "authenticated-taskset-inode")
                os.replace(replacement_binary, binary)
                os.replace(replacement_taskset, taskset)
                cpus = sorted(os.sched_getaffinity(0) - {127})
                self.assertTrue(cpus)
                cpu = cpus[0]
                command = [
                    str(taskset), "-c", str(cpu), str(binary),
                    "groupedrecovery", "--N", "3", "--overhead", "0",
                    "--seed", "15111065706836454659", "--schedule", "burst",
                ]
                returncode, stdout, stderr = campaign._run_process(
                    command, 10.0, campaign.ActiveChildren(), cpu=cpu,
                    binary_fd=binary_fd, binary_sha256=binary_sha,
                    taskset_fd=taskset_fd, taskset_sha256=taskset_sha,
                )
            finally:
                os.close(taskset_fd)
                os.close(binary_fd)
            self.assertEqual(returncode, 0)
            self.assertEqual(stderr, b"")
            self.assertIn(
                b"acf9cdeb1a5ef14f2a25f5a24fb784bdc55ccc4bcf69442a557e8e510d4a6d90",
                stdout,
            )

    def test_shard_receipt_and_stdout_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            stage = root / "stage"
            binary = root / "binary"
            taskset = root / "taskset"
            binary.write_bytes(b"binary")
            taskset.write_bytes(b"taskset")
            task = campaign._new_task("allk", 0, 0, "burst", (), [3])
            key = (0, "burst", 3)
            cells = {key: cell_expectation(key)}
            stdout = grouped_output(task, cells, {3: (0, 0)}).encode("ascii")
            stderr = b""
            shard = campaign.shard_path(stage, task)
            shard.mkdir(parents=True)
            (shard / "stdout.csv").write_bytes(stdout)
            (shard / "stderr.txt").write_bytes(stderr)
            argv = campaign.make_benchmark_argv(binary, task)
            command = [str(taskset), "-c", "0", *argv]
            campaign.write_sealed_once(
                shard / "receipt.json", f"{campaign.SCHEMA}.job_receipt",
                campaign.receipt_payload(
                    task, argv, command, 0, campaign.sha256_file(binary),
                    stdout, stderr,
                ),
            )
            receipt = campaign.load_sealed(
                shard / "receipt.json", f"{campaign.SCHEMA}.job_receipt",
            )
            self.assertEqual(receipt["execution_mode"], campaign.EXECUTION_MODE)
            campaign.validate_shard(
                stage, task, binary, campaign.sha256_file(binary), taskset,
                {0}, cells,
            )
            (shard / "stdout.csv").write_bytes(stdout.replace(b"build_ns,", b"build_nx,"))
            with self.assertRaisesRegex(campaign.CampaignError, "receipt"):
                campaign.validate_shard(
                    stage, task, binary, campaign.sha256_file(binary), taskset,
                    {0}, cells,
                )

    def test_unexpected_shard_file_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            stage = root / "stage"
            task = campaign._new_task("screen", 0, 0, "burst", [0], [3])
            shard = campaign.shard_path(stage, task)
            shard.mkdir(parents=True)
            for name in ("stdout.csv", "stderr.txt", "receipt.json", "extra"):
                (shard / name).write_bytes(b"")
            with self.assertRaisesRegex(campaign.CampaignError, "inventory"):
                campaign.validate_shard(
                    stage, task, root / "binary", "0" * 64,
                    root / "taskset", {0}, {},
                )

    def test_active_children_signal_and_late_registration(self) -> None:
        child = campaign.ActiveChildren()
        process = mock.Mock()
        process.pid = 4321
        process.poll.return_value = None
        with mock.patch.object(campaign.os, "killpg") as killpg:
            child.stop(signal.SIGTERM)
            child.register(process)
            killpg.assert_called_with(4321, signal.SIGKILL)
        with self.assertRaises(campaign.CampaignInterrupted) as caught:
            child.check()
        self.assertEqual(caught.exception.signum, signal.SIGTERM)

    def test_timeout_kills_and_reaps_process(self) -> None:
        child = campaign.ActiveChildren()
        process = mock.Mock()
        process.pid = 1234
        process.poll.return_value = None
        process.communicate.side_effect = [
            subprocess.TimeoutExpired(["fake"], 0.01), (b"", b""),
        ]
        with mock.patch.object(campaign.subprocess, "Popen", return_value=process), \
                mock.patch.object(campaign.os, "killpg") as killpg:
            with self.assertRaisesRegex(campaign.CampaignError, "timed out"):
                campaign._run_process(
                    ["/taskset", "-c", "0", "/binary", "argument"],
                    0.01, child, cpu=0, binary_fd=10,
                    binary_sha256="a" * 64, taskset_fd=11,
                    taskset_sha256="b" * 64,
                )
        killpg.assert_called_with(1234, signal.SIGKILL)
        self.assertEqual(process.communicate.call_count, 2)


class AllKOracleTests(unittest.TestCase):
    def oracle_fixture(self) -> tuple:
        task = campaign._new_task("allk", 0, 0, "burst", (), [2, 3])
        left = SimpleNamespace(
            arm="h12", overhead=0, pair=0,
            seed_index=0, seed=campaign.SEEDS[0],
            schedule="burst", ks=(2, 3),
        )
        right = SimpleNamespace(
            arm="h13", overhead=0, pair=0,
            seed_index=0, seed=campaign.SEEDS[0],
            schedule="burst", ks=(2, 3),
        )
        module = mock.Mock()
        module.validate_shard.side_effect = [
            [source_arms((0, "burst", K))["h12"] for K in (2, 3)],
            [source_arms((0, "burst", K))["h13"] for K in (2, 3)],
        ]
        oracle = campaign.StageAAllKOracle.__new__(campaign.StageAAllKOracle)
        oracle.module = module
        oracle.tasks = [left, right]
        oracle.stage = Path("/source/stage")
        oracle.contract = {
            "binary": "/source/binary", "binary_sha256": "a" * 64,
            "taskset": "/source/taskset",
        }
        oracle.allowed_cpus = {0}
        oracle.expected_shard_epoch = ("a" * 64, "b" * 64)
        return task, oracle, module

    def test_bounded_oracle_returns_exact_task_cells(self) -> None:
        task, oracle, module = self.oracle_fixture()
        with mock.patch.object(
                campaign, "source_shard_epoch_sha256",
                side_effect=lambda _module, _stage, source_task:
                    ("a" if source_task.arm == "h12" else "b") * 64):
            cells = oracle.for_task(task)
        self.assertEqual(set(cells), {(0, "burst", 2), (0, "burst", 3)})
        self.assertTrue(all(value["role"] == "allk" for value in cells.values()))
        self.assertEqual(module.validate_shard.call_count, 2)

    def test_oracle_rejects_epoch_mismatch_before_source_parse(self) -> None:
        task, oracle, module = self.oracle_fixture()
        with mock.patch.object(
                campaign, "source_shard_epoch_sha256",
                side_effect=("c" * 64, "b" * 64)):
            with self.assertRaisesRegex(campaign.CampaignError, "replay epoch"):
                oracle.for_task(task)
        module.validate_shard.assert_not_called()

    def test_oracle_rejects_source_mutation_during_parse(self) -> None:
        task, oracle, module = self.oracle_fixture()
        with mock.patch.object(
                campaign, "source_shard_epoch_sha256",
                side_effect=("a" * 64, "b" * 64, "a" * 64, "c" * 64)):
            with self.assertRaisesRegex(campaign.CampaignError, "changed while"):
                oracle.for_task(task)
        self.assertEqual(module.validate_shard.call_count, 2)

    def test_source_shard_epoch_binds_content_and_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            stage = Path(root_text)
            shard = stage / "shard"
            shard.mkdir()
            for name, data in (
                    ("stdout.csv", b"output\n"),
                    ("stderr.txt", b""),
                    ("receipt.json", b"receipt\n")):
                (shard / name).write_bytes(data)
            source_task = SimpleNamespace(payload=lambda: {"task": 1})
            module = SimpleNamespace(
                shard_path=lambda _stage, _task: shard,
            )
            first = campaign.source_shard_epoch_sha256(
                module, stage, source_task,
            )
            (shard / "stdout.csv").write_bytes(b"changed\n")
            second = campaign.source_shard_epoch_sha256(
                module, stage, source_task,
            )
            self.assertNotEqual(first, second)
            (shard / "extra").write_bytes(b"")
            with self.assertRaisesRegex(campaign.CampaignError, "inventory"):
                campaign.source_shard_epoch_sha256(module, stage, source_task)

    def test_oracle_constructor_rejects_control_leaf_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text).resolve()
            frozen = root / "frozen"
            frozen.mkdir()
            controller = frozen / "stage_a.py"
            stage = root / "stages" / "oh0000"
            hashes = {
                "contract_sha256": "1" * 64,
                "controller_affinity_sha256": "2" * 64,
                "oh0_manifest_sha256": "3" * 64,
                "oh0_complete_sha256": "4" * 64,
                "terminal_sha256": "5" * 64,
                "campaign_complete_sha256":
                    campaign.FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
                "analysis_sha256": campaign.FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
            }
            independent = {
                "controller": str(controller),
                "controller_sha256": campaign.FORMAL_STAGE_A_CONTROLLER_SHA256,
                "binary_sha256": campaign.FORMAL_STAGE_A_BINARY_SHA256,
                "taskset_sha256": "6" * 64,
            }
            expected = {
                **hashes, "result_dir": str(root),
                "controller_sha256": campaign.FORMAL_STAGE_A_CONTROLLER_SHA256,
                "binary_sha256": campaign.FORMAL_STAGE_A_BINARY_SHA256,
                "taskset_sha256": "6" * 64,
                "oh0_shard_epoch_sha256s":
                    ["7" * 64] * (2 * campaign.ALLK_TASKS),
            }
            module = mock.Mock()
            module.load_contract.return_value = independent
            module.stage_path.return_value = stage
            loaded = {"manifest": False}

            def load_manifest(_stage: Path, _overhead: int) -> list:
                loaded["manifest"] = True
                return [None] * (2 * campaign.ALLK_TASKS)

            module.load_manifest.side_effect = load_manifest

            def fake_sha(path: Path) -> str:
                path = Path(path)
                if path == controller:
                    return campaign.FORMAL_STAGE_A_CONTROLLER_SHA256
                mapping = {
                    root / "contract.json": hashes["contract_sha256"],
                    root / "controller_affinity.json":
                        hashes["controller_affinity_sha256"],
                    stage / "manifest.json": hashes["oh0_manifest_sha256"],
                    stage / "complete.json": hashes["oh0_complete_sha256"],
                    root / "terminal.json": hashes["terminal_sha256"],
                }
                value = mapping[path]
                if path == stage / "manifest.json" and loaded["manifest"]:
                    return "f" * 64
                return value

            with mock.patch.object(campaign, "sha256_file", side_effect=fake_sha), \
                    mock.patch.object(campaign, "load_sealed",
                                      return_value=independent), \
                    mock.patch.object(campaign, "require_exact_seal"), \
                    mock.patch.object(campaign, "_load_stage_a_module",
                                      return_value=module), \
                    mock.patch.object(campaign, "_source_affinity",
                                      return_value=[0]):
                with self.assertRaisesRegex(campaign.CampaignError,
                                            "control leaf changed"):
                    campaign.StageAAllKOracle(root, expected)

    def test_allk_parser_requires_oracle_cell(self) -> None:
        task = campaign._new_task("allk", 0, 0, "burst", (), [3])
        key = (0, "burst", 3)
        cells = {key: cell_expectation(key, role="allk")}
        output = grouped_output(task, cells, {3: (1, 0)})
        campaign.parse_output(output, task, cells)
        with self.assertRaisesRegex(campaign.CampaignError, "source oracle"):
            campaign.parse_output(output, task, {})

    def test_allk_guard_rejects_nonpassing_screen_before_telemetry(self) -> None:
        with mock.patch.object(campaign, "load_manifest", return_value=[]), \
                mock.patch.object(campaign, "verify_complete_screen",
                                  return_value={"metrics": {}}), \
                mock.patch.object(campaign, "verify_screen_decision",
                                  return_value={"screen_passed": False}), \
                mock.patch.object(campaign, "verify_screen_terminal",
                                  return_value={"allk_escalation_planned": False}), \
                mock.patch.object(campaign, "telemetry_start") as telemetry:
            with self.assertRaisesRegex(campaign.CampaignError, "passing screen"):
                campaign.verify_passing_screen_guard(
                    Path("/synthetic"), {}, {}, {0},
                )
            telemetry.assert_not_called()


if __name__ == "__main__":
    unittest.main()
