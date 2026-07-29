#!/usr/bin/env python3
"""Warning-strict tests for the result-free rv4a repair campaign."""

import copy
import gzip
import hashlib
import importlib
import importlib.util
import io
import json
import math
import multiprocessing
import os
from pathlib import Path
import py_compile
import signal
import stat
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock
from contextlib import redirect_stderr, redirect_stdout


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "wh2_rv4a_campaign_under_test",
    HERE / "wh2_rv4a_campaign.py",
)
campaign = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(campaign)


_PARALLEL_RELEASED = None
_PARALLEL_STARTED_BEFORE_RELEASE = None


def _parallel_slow_first_task(task):
    index, unused_assignment, unused_job = task
    if index == 0:
        time.sleep(1.0)
        _PARALLEL_RELEASED.set()
    elif not _PARALLEL_RELEASED.is_set():
        with _PARALLEL_STARTED_BEFORE_RELEASE.get_lock():
            _PARALLEL_STARTED_BEFORE_RELEASE.value += 1
    return index, {"index": index}


def _parallel_error_task(task):
    index, unused_assignment, unused_job = task
    if index == 0:
        raise RuntimeError("synthetic replay failure")
    time.sleep(0.2)
    return index, {"index": index}


def _ignore_sigterm_forever(ready):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    ready.set()
    while True:
        time.sleep(0.05)


def _synthetic_runtime_bindings():
    source_limits = {
        "native": campaign.NATIVE_SOURCE_BYTE_LIMIT,
        "build": campaign.BUILD_SOURCE_BYTE_LIMIT,
        "python": campaign.PYTHON_SOURCE_BYTE_LIMIT,
        "configuration": campaign.BUILD_CONFIGURATION_BYTE_LIMIT,
    }

    def item(name, byte_limit, index):
        return {
            "path": f"/tmp/rv4a-preflight/{name}",
            "device": 1,
            "inode": index,
            "mode": stat.S_IFREG | 0o644,
            "size": 1,
            "mtime_ns": 1,
            "ctime_ns": 1,
            "sha256": f"{index:064x}",
            "byte_limit": byte_limit,
        }

    return {
        "benchmark": item(
            "benchmark", campaign.BENCHMARK_BINARY_BYTE_LIMIT, 1),
        "sources": {
            name: item(name, source_limits[source_class], index)
            for index, (name, source_class) in enumerate(
                campaign._RUNTIME_PIN_SOURCE_CLASSES.items(), start=2)
        },
        "thermal": {
            "path": "/tmp/rv4a-preflight/thermal",
            "device": 1,
            "inode": 100,
            "mode": stat.S_IFREG | 0o644,
        },
    }


class FrozenPlanTests(unittest.TestCase):
    def test_exact_cardinalities(self):
        validation = campaign.validate_frozen_contract()
        self.assertEqual(validation["cardinalities"], {
            "training": {
                "jobs": 1188,
                "recovery_cells": 912384,
                "unique_selectors": 152064,
                "attempt_rows": 7299072,
                "native_rows_per_job": 92928,
                "native_rows": 110398464,
            },
            "sealed": {
                "jobs": 1188,
                "recovery_cells": 304128,
                "unique_selectors": 25443,
                "attempt_rows": 2433024,
                "native_rows_per_job": 30976,
                "native_rows": 36799488,
            },
        })
        self.assertEqual(len(campaign.build_pre_cpu_jobs("training")), 1188)
        for arm in campaign.ARM_BY_NAME:
            self.assertEqual(
                len(campaign.build_pre_cpu_jobs("sealed", arm)), 1188)

    def test_all_jobs_are_unique_valid_and_deterministic(self):
        for phase, survivor in (
            ("training", None),
            ("sealed", "pure8_s0m1_d3"),
            ("sealed", "pure9_s0m1_d3"),
        ):
            first = campaign.build_pre_cpu_jobs(phase, survivor)
            second = campaign.build_pre_cpu_jobs(phase, survivor)
            self.assertEqual(first, second)
            self.assertEqual(
                len({job["job_id"] for job in first}), len(first))
            self.assertEqual(
                [job["order"] for job in first], list(range(len(first))))
            for job in first:
                self.assertEqual(
                    campaign.canonical_sha256({
                        key: value for key, value in job.items()
                        if key not in ("job_id", "order")
                    }),
                    job["job_id"],
                )
                campaign.validate_job(job)

    def test_training_is_exactly_paired_across_arms(self):
        jobs = campaign.build_pre_cpu_jobs("training")
        coordinates = {}
        for job in jobs:
            key = (job["K"], job["schedule"])
            signature = (
                job["construction_seed_base"],
                job["construction_seed_derivation"],
                job["loss_seed_base"],
                job["loss_seed_derivation"],
                job["replicates"],
            )
            coordinates.setdefault(key, []).append(signature)
        self.assertEqual(
            len(coordinates),
            len(campaign.K_VALUES) * len(campaign.SCHEDULES),
        )
        self.assertEqual({len(items) for items in coordinates.values()}, {2})
        self.assertTrue(all(len(set(items)) == 1
                            for items in coordinates.values()))

    def test_sealed_is_one_arm_and_both_lanes(self):
        for survivor in campaign.ARM_BY_NAME:
            jobs = campaign.build_pre_cpu_jobs("sealed", survivor)
            self.assertEqual({job["arm"] for job in jobs}, {survivor})
            self.assertEqual({job["lane"] for job in jobs},
                             {"random", "production"})
            for block_count in campaign.K_VALUES:
                random_jobs = [
                    job for job in jobs
                    if job["K"] == block_count and job["lane"] == "random"
                ]
                production_jobs = [
                    job for job in jobs
                    if job["K"] == block_count and
                    job["lane"] == "production"
                ]
                self.assertEqual(len(random_jobs), 6)
                self.assertEqual(len(production_jobs), 6)
                self.assertEqual(
                    {job["construction_seed_derivation"]
                     for job in random_jobs},
                    {campaign.CONSTRUCTION_DERIVED},
                )
                self.assertEqual(
                    {job["construction_seed_derivation"]
                     for job in production_jobs},
                    {campaign.CONSTRUCTION_FIXED},
                )
                self.assertEqual(
                    {job["construction_seed_base"]
                     for job in production_jobs},
                    {campaign._SEEDS["construction"][
                        "sealed_production"]["roots_by_K"][
                            str(block_count)]["root"]},
                )

    def test_phase_and_survivor_fail_closed(self):
        with self.assertRaises(campaign.CampaignError):
            campaign.build_pre_cpu_jobs("training", "pure8_s0m1_d3")
        with self.assertRaises(campaign.CampaignError):
            campaign.build_pre_cpu_jobs("sealed")
        with self.assertRaises(campaign.CampaignError):
            campaign.build_pre_cpu_jobs("sealed", "not-an-arm")
        with self.assertRaises(campaign.CampaignError):
            campaign.build_pre_cpu_jobs("other")

    def test_policy_is_repair_specific_and_decoder_solve_first(self):
        policy = campaign.selection_policy()
        self.assertEqual(
            policy["selector"]["policy_sha256"],
            campaign.REPAIR_POLICY_SHA256,
        )
        self.assertEqual(policy["ranking"]["lexicographic"], [
            "decoder_feed_candidate_wh1_log_cost",
            "full_decoder_candidate_wh1_log_cost",
            "selected_direct_candidate_dispatch_log_cost",
            "full_encoder_candidate_wh1_log_cost",
            "candidate_name",
        ])
        self.assertIn(
            "selector_cap_exhaustion",
            policy["eligibility"]["zero"],
        )
        self.assertIn(
            "selector_or-codec_OOM",
            policy["eligibility"]["zero"],
        )
        self.assertIn(
            "candidate_raw_attempt0",
            json.dumps(policy, sort_keys=True),
        )
        self.assertEqual(
            policy["raw_evidence"]["historical_references"],
            {
                "pure8_s0m1_d3": 1208,
                "pure9_s0m1_d3": 1034,
            },
        )

    def test_block_width_is_not_a_roster_axis(self):
        roster = campaign.frozen_roster()
        self.assertEqual(roster["request"]["block_bytes"], 2)
        self.assertEqual(
            roster["correctness_prelaunch"]["block_bytes"],
            [2, 6, 32, 64, 256, 1280, 4096],
        )
        self.assertEqual(
            roster["correctness_prelaunch"]["rule"],
            "not-a-campaign-roster-axis",
        )
        self.assertTrue(all(job["block_bytes"] == 2 for job in
                            campaign.build_pre_cpu_jobs("training")))

    def test_plan_hashes_cover_jobs_cells_policy_and_seeds(self):
        validation = campaign.validate_frozen_contract()
        self.assertEqual(
            {
                key: validation[key]
                for key in (
                    "frozen_roster_sha256",
                    "policy_sha256",
                    "seed_disjointness_sha256",
                    "training_job_list_sha256",
                    "training_cell_set_sha256",
                    "sealed_job_list_sha256",
                    "sealed_cell_set_sha256",
                    "result_free_plan_sha256",
                )
            },
            {
                "frozen_roster_sha256":
                    "c4d29e2cacfe8fdf762dc449e803a8d5bfcdd7817b948c59"
                    "da61580f219d0d7a",
                "policy_sha256":
                    "81d04c117d8fa414971057f8bd69a39d320139a872eac9c78"
                    "5f5ff8c73c02b20",
                "seed_disjointness_sha256":
                    "ca29960e58a231f2d0bb7212843ccdee84e7ff6e6c122e18f"
                    "d29873e2bd8e7f2",
                "training_job_list_sha256":
                    "70e0e4cc184fe478e53d8b877b35b434134d1d38a10b7d24"
                    "c8eb52ecde259c0e",
                "training_cell_set_sha256":
                    "9dc8406a564df1f1d2bf5572e0580688a0dc0032453d60ecc"
                    "31715e576e7d3cb",
                "sealed_job_list_sha256": {
                    "pure8_s0m1_d3":
                        "dc7713a1ad4a4dcd4866c1617fc5b80b7e9b8011fc2e817b"
                        "86f09f919ed31199",
                    "pure9_s0m1_d3":
                        "4a9b72575f7a005cd839ab0c705c56a40b445e11e0122e46d"
                        "bdfebb7f65ce2e9",
                },
                "sealed_cell_set_sha256": {
                    "pure8_s0m1_d3":
                        "b54169a3a539087b96f58934d6ea5009e6d7a6d382bb8bb3"
                        "c0677492d43f9d46",
                    "pure9_s0m1_d3":
                        "c163f8f37a89dc9437304be3e4ea3ff757077049e8c55d45d"
                        "04fe80e901e76b2",
                },
                "result_free_plan_sha256":
                    "568a967393ab1b5a23f9cc518b2e87f91752972c0f5c94d4c"
                    "00099f4d9087291",
            },
        )
        for key in (
            "frozen_roster_sha256",
            "policy_sha256",
            "seed_disjointness_sha256",
            "training_job_list_sha256",
            "training_cell_set_sha256",
            "result_free_plan_sha256",
        ):
            self.assertTrue(campaign._is_sha256(validation[key]), key)
        self.assertEqual(
            set(validation["sealed_job_list_sha256"]),
            set(campaign.ARM_BY_NAME),
        )
        self.assertEqual(
            set(validation["sealed_cell_set_sha256"]),
            set(campaign.ARM_BY_NAME),
        )
        plan = campaign.make_plan("training")
        changed = copy.deepcopy(plan)
        changed["frozen_roster"]["request"]["block_bytes"] = 4
        self.assertNotEqual(
            campaign.canonical_sha256(changed),
            campaign.canonical_sha256(plan),
        )

    def test_result_free_plan_does_not_need_parser(self):
        saved = campaign.repair
        campaign.repair = None
        try:
            plan = campaign.make_plan("training")
        finally:
            campaign.repair = saved
        self.assertEqual(plan["phase"], "training")
        self.assertEqual(len(plan["pre_cpu_jobs"]), 1188)

    def test_source_forced_plan_cli_executes_exact_runner(self):
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                str(HERE / "wh2_rv4a_campaign.py"),
                "plan",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        plan = json.loads(completed.stdout)
        self.assertEqual(
            plan["schema"], "wirehair.wh2.rv4a.result-free-plan.v1")
        self.assertEqual(
            plan["cardinalities"]["training"]["unique_selectors"],
            152064,
        )
        self.assertEqual(
            hashlib.sha256(completed.stdout).hexdigest(),
            "762a3c01984cb57b42f8fb771ba13a14e81f9fedb17e5451c"
            "16d93743a094283",
        )
        self.assertEqual(completed.stderr, b"")

    def test_public_cli_exposes_no_candidate_or_lane_override(self):
        parser = campaign._parser()
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(["plan", "--survivor",
                                   "pure8_s0m1_d3"])
            with self.assertRaises(SystemExit):
                parser.parse_args([
                    "run-sealed",
                    "--directory", "/tmp/x",
                    "--bench", "/tmp/b",
                    "--thermal-source", "/tmp/t",
                    "--training-directory", "/tmp/p",
                    "--survivor", "pure8_s0m1_d3",
                ])
            with self.assertRaises(SystemExit):
                parser.parse_args([
                    "run-sealed",
                    "--directory", "/tmp/x",
                    "--bench", "/tmp/b",
                    "--thermal-source", "/tmp/t",
                    "--training-directory", "/tmp/p",
                    "--lane", "random",
                ])

    def test_run_cli_requires_preflight_pin_and_wires_it(self):
        parser = campaign._parser()
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args([
                    "run-training",
                    "--directory", "/tmp/x",
                    "--bench", "/tmp/b",
                    "--thermal-source", "/tmp/t",
                ])
        args = parser.parse_args([
            "run-training",
            "--directory", "/tmp/x",
            "--bench", "/tmp/b",
            "--thermal-source", "/tmp/t",
            "--pin-record", "/tmp/pins.json",
            "--pin-sha256", "a" * 64,
        ])
        self.assertEqual(args.pin_record, Path("/tmp/pins.json"))
        self.assertEqual(args.pin_sha256, "a" * 64)
        preflight = parser.parse_args([
            "preflight-pin",
            "--bench", "/tmp/b",
            "--thermal-source", "/tmp/t",
        ])
        self.assertEqual(preflight.command, "preflight-pin")

    def test_cli_preflight_to_run_wiring_without_outcomes(self):
        record = {
            "schema": campaign.RUNTIME_PIN_SCHEMA,
            "runtime_files": {},
        }
        output = io.StringIO()
        with (
            mock.patch.object(
                campaign, "make_preflight_pin", return_value=record
            ) as preflight,
            redirect_stdout(output),
        ):
            campaign.main([
                "preflight-pin",
                "--bench", "/tmp/b",
                "--thermal-source", "/tmp/t",
            ])
        preflight.assert_called_once_with(
            Path("/tmp/b"), Path("/tmp/t"))
        self.assertEqual(
            output.getvalue(),
            campaign.canonical_json_bytes(record).decode("ascii"),
        )

        with mock.patch.object(campaign, "run_campaign") as run:
            campaign.main([
                "run-training",
                "--directory", "/tmp/x",
                "--bench", "/tmp/b",
                "--thermal-source", "/tmp/t",
                "--pin-record", "/tmp/pins.json",
                "--pin-sha256", "a" * 64,
                "--workers", "3",
                "--replay-workers", "2",
            ])
        run.assert_called_once_with(
            "training",
            Path("/tmp/x"),
            Path("/tmp/b"),
            Path("/tmp/t"),
            pin_record=Path("/tmp/pins.json"),
            pin_sha256="a" * 64,
            workers=3,
            replay_workers=2,
        )


class SeedDomainTests(unittest.TestCase):
    def test_kdf_framing_is_unambiguous_and_coordinate_sensitive(self):
        domain = campaign.SEED_DOMAINS["training_construction"]
        self.assertNotEqual(
            campaign.domain_value(domain, 32, coordinates=(1, 23)),
            campaign.domain_value(domain, 32, coordinates=(12, 3)),
        )
        self.assertNotEqual(
            campaign.domain_value(domain, 32, coordinates=(1,), counter=0),
            campaign.domain_value(domain, 32, coordinates=(1,), counter=1),
        )
        self.assertNotEqual(
            campaign.domain_value(domain, 32),
            campaign.domain_value(domain, 64) & 0xFFFFFFFF,
        )
        with self.assertRaises(campaign.CampaignError):
            campaign.domain_value(domain, 16)
        with self.assertRaises(campaign.CampaignError):
            campaign.domain_value("non-ascii-\N{SNOWMAN}", 32)

    def test_exact_domains_are_distinct_ascii(self):
        domains = list(campaign.SEED_DOMAINS.values())
        self.assertEqual(len(domains), 8)
        self.assertEqual(len(set(domains)), 8)
        self.assertTrue(all(domain.isascii() for domain in domains))
        self.assertTrue(all(
            domain.startswith("wirehair.wh2.rv4a.") and
            domain.endswith(".v1")
            for domain in domains
        ))

    def test_root_and_loss_populations_are_disjoint(self):
        private = campaign._SEEDS["_private"]
        root_populations = [
            set(private["training_roots"]),
            set(private["sealed_random_roots"]),
            set(private["production_roots"]),
        ]
        for left_index, left in enumerate(root_populations):
            for right in root_populations[left_index + 1:]:
                self.assertFalse(left & right)
        loss_populations = [
            set(private["training_losses"]),
            set(private["sealed_random_losses"]),
            set(private["sealed_production_losses"]),
        ]
        for left_index, left in enumerate(loss_populations):
            for right in loss_populations[left_index + 1:]:
                self.assertFalse(left & right)

    def test_expanded_attempt_spaces_are_disjoint_from_prior(self):
        prior = (
            campaign._root_population(
                campaign.ZA5V_SELECTION_CONSTRUCTION_BASE,
                768,
                campaign.CONSTRUCTION_DERIVED,
            ),
            campaign._root_population(
                campaign.ZA5V_HOLDOUT_CONSTRUCTION_BASE,
                256,
                campaign.CONSTRUCTION_DERIVED,
            ),
            campaign._root_population(
                campaign.ZA5V_SCREEN_CONSTRUCTION_BASE,
                3,
                campaign.CONSTRUCTION_DERIVED,
            ),
        )
        private = campaign._SEEDS["_private"]
        populations = list(prior) + [
            private["training_roots"],
            private["sealed_random_roots"],
            private["production_roots"],
        ]
        expanded = [
            campaign._expanded_sets(population)
            for population in populations
        ]
        for population, (precode, packet) in zip(populations, expanded):
            self.assertEqual(
                len(precode), len(population) * campaign.REPAIR_ATTEMPTS)
            self.assertEqual(
                len(packet), len(population) * campaign.REPAIR_ATTEMPTS)
        for left_index, (left_precode, left_packet) in enumerate(expanded):
            for right_precode, right_packet in expanded[left_index + 1:]:
                self.assertFalse(left_precode & right_precode)
                self.assertFalse(left_packet & right_packet)

    def test_production_roots_and_rejection_counts_are_recorded(self):
        production = campaign.seed_disjointness_proof()[
            "construction"]["sealed_production"]
        self.assertEqual(production["root_count"], 99)
        self.assertEqual(
            set(production["roots_by_K"]),
            {str(block_count) for block_count in campaign.K_VALUES},
        )
        roots = []
        for block_count in campaign.K_VALUES:
            item = production["roots_by_K"][str(block_count)]
            self.assertEqual(item["coordinates"], [block_count])
            self.assertIs(type(item["rejection_count"]), int)
            self.assertGreaterEqual(item["rejection_count"], 0)
            self.assertEqual(
                item["expanded_precode64_count"],
                campaign.REPAIR_ATTEMPTS,
            )
            self.assertEqual(
                item["expanded_packet32_count"],
                campaign.REPAIR_ATTEMPTS,
            )
            precode, packet = campaign._expanded_sets((item["root"],))
            self.assertEqual(len(precode), campaign.REPAIR_ATTEMPTS)
            self.assertEqual(len(packet), campaign.REPAIR_ATTEMPTS)
            roots.append(item["root"])
        self.assertEqual(len(set(roots)), 99)

    def test_proof_does_not_expose_private_lists(self):
        proof = campaign.seed_disjointness_proof()
        self.assertNotIn("_private", proof)
        self.assertNotIn("_roots", json.dumps(proof, sort_keys=True))
        self.assertNotIn("_losses", json.dumps(proof, sort_keys=True))


class MutationTests(unittest.TestCase):
    def setUp(self):
        self.job = dict(campaign.build_pre_cpu_jobs("training")[0])

    def assert_rejected(self, key, value):
        changed = dict(self.job)
        changed[key] = value
        with self.assertRaises(campaign.CampaignError):
            campaign.validate_job(changed)

    def test_job_mutations_are_rejected(self):
        mutations = {
            "phase": "sealed",
            "lane": "random",
            "arm": "nope",
            "arm_provisional_id": "0" * 16,
            "repair_policy_sha256": "0" * 64,
            "K": 101,
            "schedule": "easy",
            "block_bytes": 4,
            "loss": 0.2,
            "construction_seed_base": 0,
            "construction_seed_derivation": campaign.CONSTRUCTION_FIXED,
            "loss_seed_base": 0,
            "replicates": 767,
            "inner_reps": 2,
            "max_overhead": 63,
            "cache_state": "cold",
            "order": -1,
            "job_id": "0" * 64,
        }
        for key, value in mutations.items():
            with self.subTest(key=key):
                self.assert_rejected(key, value)

    def test_typed_json_rejects_boolean_integer_substitution(self):
        self.assert_rejected("K", True)
        self.assert_rejected("replicates", True)
        self.assert_rejected("loss", 1)


class AtomicEvidenceTests(unittest.TestCase):
    def test_atomic_json_refuses_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "value.json"
            campaign.atomic_json(path, {"a": 1})
            with self.assertRaises(campaign.CampaignError):
                campaign.atomic_json(path, {"a": 2})
            self.assertEqual(
                campaign.read_canonical_json(path)[0], {"a": 1})

    def test_gzip_round_trip_and_limits(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "value.json.gz"
            evidence = campaign.atomic_gzip_json(path, {"a": [1, 2, 3]})
            value, replayed = campaign.read_canonical_gzip_json(
                path,
                compressed_limit=evidence["compressed_bytes"],
                uncompressed_limit=evidence["uncompressed_bytes"],
            )
            self.assertEqual(value, {"a": [1, 2, 3]})
            self.assertEqual(evidence, replayed)
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_gzip_json(
                    path,
                    compressed_limit=evidence["compressed_bytes"] - 1,
                    uncompressed_limit=evidence["uncompressed_bytes"],
                )
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_gzip_json(
                    path,
                    compressed_limit=evidence["compressed_bytes"],
                    uncompressed_limit=evidence["uncompressed_bytes"] - 1,
                )

    def test_gzip_rejects_uncompressed_size_before_compression(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "oversized.json.gz"
            with (
                mock.patch.object(
                    campaign.gzip,
                    "compress",
                    side_effect=AssertionError("must not compress"),
                ) as compressor,
                self.assertRaises(campaign.CampaignError),
            ):
                campaign.atomic_gzip_json(
                    path,
                    {"payload": "large"},
                    limits={"compressed": 1024, "uncompressed": 1},
                )
            compressor.assert_not_called()
            self.assertFalse(path.exists())

    def test_gzip_rejects_truncation_trailing_and_concatenation(self):
        payload = campaign.canonical_json_bytes({"x": 1})
        member = gzip.compress(payload, mtime=0)
        variants = {
            "truncated": member[:-1],
            "trailing": member + b"x",
            "concatenated": member + member,
        }
        with tempfile.TemporaryDirectory() as temporary:
            for name, value in variants.items():
                with self.subTest(name=name):
                    path = Path(temporary) / f"{name}.gz"
                    path.write_bytes(value)
                    with self.assertRaises(campaign.CampaignError):
                        campaign.read_canonical_gzip_json(
                            path,
                            compressed_limit=len(value),
                            uncompressed_limit=2 * len(payload),
                        )

    def test_gzip_rejects_noncanonical_and_duplicate_json(self):
        variants = (
            b'{ "x":1}\\n',
            b'{"x":1,"x":2}\n',
            b'{"x":NaN}\n',
        )
        with tempfile.TemporaryDirectory() as temporary:
            for index, payload in enumerate(variants):
                path = Path(temporary) / f"{index}.gz"
                path.write_bytes(gzip.compress(payload, mtime=0))
                with self.assertRaises(campaign.CampaignError):
                    campaign.read_canonical_gzip_json(
                        path,
                        compressed_limit=1024,
                        uncompressed_limit=1024,
                    )

    def test_json_lines_round_trip_caps_and_concatenation(self):
        payload = (
            campaign.canonical_json_bytes({"x": 1}) +
            campaign.canonical_json_bytes({"x": 2})
        )
        compressed = gzip.compress(payload, mtime=0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "ledger.gz"
            path.write_bytes(compressed)
            rows, evidence = campaign.read_canonical_gzip_json_lines(
                path,
                compressed_limit=len(compressed),
                uncompressed_limit=len(payload),
                row_limit=2,
            )
            self.assertEqual(rows, [{"x": 1}, {"x": 2}])
            self.assertEqual(evidence["rows"], 2)
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_gzip_json_lines(
                    path,
                    compressed_limit=len(compressed),
                    uncompressed_limit=len(payload),
                    row_limit=1,
                )
            path.unlink()
            path.write_bytes(compressed + compressed)
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_gzip_json_lines(
                    path,
                    compressed_limit=2 * len(compressed),
                    uncompressed_limit=2 * len(payload),
                    row_limit=4,
                )

    def test_stable_readers_reject_symlinks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target"
            target.write_bytes(campaign.canonical_json_bytes({"x": 1}))
            link = root / "link"
            link.symlink_to(target)
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_json(link)
            with self.assertRaises(campaign.CampaignError):
                campaign._stream_file_sha256(link, byte_limit=1024)
            compressed = root / "target.gz"
            compressed.write_bytes(gzip.compress(
                campaign.canonical_json_bytes({"x": 1}), mtime=0))
            compressed_link = root / "link.gz"
            compressed_link.symlink_to(compressed)
            with self.assertRaises(campaign.CampaignError):
                campaign.read_canonical_gzip_json(
                    compressed_link,
                    compressed_limit=1024,
                    uncompressed_limit=1024,
                )

    def test_stable_reader_rejects_descriptor_mutation(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "value"
            path.write_bytes(b"abc")
            descriptor = path.stat()
            changed = SimpleNamespace(**{
                name: getattr(descriptor, name)
                for name in (
                    "st_dev", "st_ino", "st_mode", "st_size",
                    "st_mtime_ns", "st_ctime_ns",
                )
            })
            changed.st_mtime_ns += 1
            real_fstat = campaign.os.fstat
            calls = 0

            def unstable(fd):
                nonlocal calls
                calls += 1
                return descriptor if calls == 1 else changed

            with mock.patch.object(campaign.os, "fstat", side_effect=unstable):
                with self.assertRaises(campaign.CampaignError):
                    campaign._read_stable_bytes(path, byte_limit=16)
            self.assertIsNotNone(real_fstat)

    def test_file_binding_enforces_source_cap(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.py"
            path.write_bytes(b"x" * 17)
            with self.assertRaises(campaign.CampaignError):
                campaign._stable_file_binding(path, byte_limit=16)
            binding = campaign._stable_file_binding(path, byte_limit=17)
            self.assertEqual(binding["size"], 17)
            self.assertEqual(binding["byte_limit"], 17)

    def test_file_binding_rejects_fifo_before_read(self):
        if not hasattr(os, "mkfifo"):
            self.skipTest("FIFO is unavailable")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fifo"
            os.mkfifo(path)
            started = time.monotonic()
            with self.assertRaises(campaign.CampaignError):
                campaign._stable_file_binding(path, byte_limit=1024)
            self.assertLess(time.monotonic() - started, 1.0)

    def test_thermal_projection_omits_raw_rows(self):
        receipt = {
            "context": {
                "bound": {
                    "cpu_affinity": [3],
                    "thermal_device": 4,
                    "thermal_inode": 5,
                },
                "thermal": {
                    "valid_rows": 2,
                    "rows": [{"large": "raw"}] * 100,
                    "cpu_tctl_max_c": 65.0,
                    "dimm_max_c": 50.0,
                    "edac_ce_max": 0,
                    "edac_ue_max": 0,
                    "first_monotonic_s": 1.0,
                    "last_monotonic_s": 2.0,
                },
            },
        }
        projected = campaign._receipt_context_projection(
            receipt, {"cpu": 3})
        self.assertNotIn("rows", projected)
        self.assertEqual(projected["valid_rows"], 2)

    def test_preflight_pin_exactly_matches_manifest_projection(self):
        bindings = _synthetic_runtime_bindings()
        expected = campaign.make_runtime_pin_record(bindings)
        with mock.patch.object(
                campaign, "runtime_bindings", return_value=bindings):
            preflight = campaign.make_preflight_pin(
                "/unused/bench", "/unused/thermal")
        self.assertEqual(
            campaign.canonical_json_bytes(preflight),
            campaign.canonical_json_bytes(expected),
        )
        manifest_projection = campaign.make_runtime_pin_record(bindings)
        self.assertEqual(preflight, manifest_projection)
        pin_digest = campaign.canonical_sha256(preflight)
        jobs = campaign.build_pre_cpu_jobs("training")
        manifest = campaign._build_manifest(
            "training",
            None,
            jobs,
            [0],
            bindings,
            None,
            pin_digest,
        )
        campaign._validate_manifest(manifest)
        self.assertEqual(
            manifest["preflight_pin_sha256"], pin_digest)

        verifier_spec = importlib.util.spec_from_file_location(
            "_rv4a_pin_projection_test",
            HERE / "wh2_rv4a_parallel_verify.py",
        )
        verifier = importlib.util.module_from_spec(verifier_spec)
        verifier_spec.loader.exec_module(verifier)
        verifier_projection = verifier.make_pin_record(manifest)
        self.assertEqual(
            campaign.canonical_json_bytes(preflight),
            verifier.canonical_json_bytes(verifier_projection),
        )

    def test_wrong_preflight_pin_fails_before_directory_or_workers(self):
        bindings = _synthetic_runtime_bindings()
        pinned = campaign.make_runtime_pin_record(bindings)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pin_path = root / "pins.json"
            pin_path.write_bytes(campaign.canonical_json_bytes(pinned))
            digest = campaign.canonical_sha256(pinned)
            for label in ("runner-content", "thermal-identity"):
                changed = copy.deepcopy(bindings)
                if label == "runner-content":
                    changed["sources"]["runner"]["sha256"] = "f" * 64
                else:
                    changed["thermal"]["inode"] += 1
                output = root / label
                with (
                    self.subTest(label=label),
                    mock.patch.object(
                        campaign, "_require_normal_priority"
                    ),
                    mock.patch.object(
                        campaign, "runtime_bindings", return_value=changed
                    ),
                    mock.patch.object(
                        campaign,
                        "_BOOTSTRAP_RUNNER_SOURCE_SHA256",
                        changed["sources"]["runner"]["sha256"],
                    ),
                    self.assertRaises(campaign.CampaignError),
                ):
                    campaign._run_campaign_impl(
                        "training",
                        output,
                        "/unused/bench",
                        "/unused/thermal",
                        pin_record=pin_path,
                        pin_sha256=digest,
                        workers=1,
                        replay_workers=1,
                    )
                self.assertFalse(output.exists())

    def test_plain_import_cannot_create_outcomes_or_start_workers(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "campaign"
            with (
                mock.patch.object(
                    campaign, "_BOOTSTRAP_RUNNER_SOURCE_SHA256", None
                ),
                mock.patch.object(
                    campaign, "runtime_bindings"
                ) as bindings,
                mock.patch.object(
                    campaign.subprocess, "Popen"
                ) as popen,
                self.assertRaises(campaign.CampaignError),
            ):
                campaign._run_campaign_impl(
                    "training",
                    output,
                    "/unused/bench",
                    "/unused/thermal",
                    pin_record="/unused/pin",
                    pin_sha256="a" * 64,
                    workers=1,
                    replay_workers=1,
                )
            bindings.assert_not_called()
            popen.assert_not_called()
            self.assertFalse(output.exists())

    def test_plain_import_cannot_replay_or_verify_outcomes(self):
        with (
            mock.patch.object(
                campaign, "_BOOTSTRAP_RUNNER_SOURCE_SHA256", None
            ),
            mock.patch.object(
                campaign, "_verify_campaign_impl"
            ) as verify_impl,
            self.assertRaises(campaign.CampaignError),
        ):
            campaign.verify_campaign("/unused", replay_workers=1)
        verify_impl.assert_not_called()

        with (
            mock.patch.object(
                campaign, "_BOOTSTRAP_RUNNER_SOURCE_SHA256", None
            ),
            mock.patch.object(campaign, "_validate_manifest") as validate,
            self.assertRaises(campaign.CampaignError),
        ):
            campaign.build_summaries(
                "/unused",
                {},
                "/unused-output",
                replay_workers=1,
            )
        validate.assert_not_called()

    def test_outcome_gate_requires_manifest_runner_hash(self):
        with mock.patch.object(
            campaign,
            "_BOOTSTRAP_RUNNER_SOURCE_SHA256",
            "a" * 64,
        ):
            campaign._require_source_forced_outcome_runner("a" * 64)
            with self.assertRaises(campaign.CampaignError):
                campaign._require_source_forced_outcome_runner("b" * 64)

    def test_loaded_source_stamp_detects_prebinding_mutation(self):
        bindings = _synthetic_runtime_bindings()
        parser = SimpleNamespace(
            __rv4a_source_sha256__="a" * 64)
        bindings["sources"]["parser"]["sha256"] = "a" * 64
        bindings["sources"]["context_tool"]["sha256"] = \
            campaign.peel_codec.__rv4a_source_sha256__
        campaign._validate_loaded_source_hashes(bindings, parser)
        bindings["sources"]["parser"]["sha256"] = "b" * 64
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_loaded_source_hashes(bindings, parser)
        bindings["sources"]["parser"]["sha256"] = "a" * 64
        bindings["sources"]["runner"]["sha256"] = "b" * 64
        with (
            mock.patch.object(
                campaign,
                "_BOOTSTRAP_RUNNER_SOURCE_SHA256",
                "c" * 64,
            ),
            self.assertRaises(campaign.CampaignError),
        ):
            campaign._validate_loaded_source_hashes(bindings, parser)

    def test_source_loader_ignores_ordinary_and_alternate_pyc(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "poisoned.py"
            forged = "VALUE = 'forged'\n"
            actual = "VALUE = 'source'\n"
            self.assertEqual(len(forged), len(actual))
            source.write_text(forged, encoding="ascii")
            source_stat = source.stat()
            ordinary = Path(importlib.util.cache_from_source(str(source)))
            ordinary.parent.mkdir(parents=True)
            py_compile.compile(
                str(source), cfile=str(ordinary), doraise=True)
            alternate_root = root / ".rv4a-no-bytecode-cache"
            previous_prefix = sys.pycache_prefix
            try:
                sys.pycache_prefix = str(alternate_root)
                alternate = Path(
                    importlib.util.cache_from_source(str(source)))
                alternate.parent.mkdir(parents=True)
                py_compile.compile(
                    str(source), cfile=str(alternate), doraise=True)
            finally:
                sys.pycache_prefix = previous_prefix
            source.write_text(actual, encoding="ascii")
            os.utime(
                source,
                ns=(source_stat.st_atime_ns, source_stat.st_mtime_ns),
            )
            loaded = campaign._bootstrap_source_module(
                "_rv4a_pyc_poison_test", source)
            try:
                self.assertEqual(loaded.VALUE, "source")
                self.assertTrue(ordinary.exists())
                self.assertTrue(alternate.exists())
            finally:
                sys.modules.pop("_rv4a_pyc_poison_test", None)

    def test_v3_consumers_reject_external_context_module(self):
        saved_repair = campaign.repair
        saved_context = campaign.peel_codec
        campaign.repair = SimpleNamespace(
            REPAIRTIMING_PROTOCOL_V3=
                campaign.REQUIRED_REPAIRTIMING_PROTOCOL,
            REPAIRTIMING_SCHEMA_V3=
                campaign.REQUIRED_REPAIRTIMING_SCHEMA,
            REPAIRTIMING_ROWS_PER_CELL=121,
            __file__=str(campaign.TOOLS / "repair_timing_codec.py"),
        )
        campaign.peel_codec = SimpleNamespace(
            __file__="/tmp/malicious/peel_codec.py")
        try:
            with self.assertRaises(campaign.CampaignError):
                campaign._require_v3_parser()
        finally:
            campaign.repair = saved_repair
            campaign.peel_codec = saved_context

    def test_empty_log_verification_rejects_symlink(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target"
            target.write_bytes(b"")
            link = root / "worker.log"
            link.symlink_to(target)
            manifest = {"assignments": [{"log": "worker.log"}]}
            with self.assertRaises(campaign.CampaignError):
                campaign._verify_empty_logs(root, manifest)

    def test_atomic_jsonl_limits_and_error_cleanup(self):
        row = {"x": "abc"}
        payload_size = len(campaign.canonical_json_bytes(row))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            exact = root / "exact.jsonl.gz"
            with campaign.AtomicGzipJsonLines(
                    exact,
                    row_limit=1,
                    compressed_limit=1024,
                    uncompressed_limit=payload_size,
            ) as writer:
                writer.write(row)
            self.assertTrue(exact.is_file())

            too_many = root / "too-many.jsonl.gz"
            with self.assertRaises(campaign.CampaignError):
                with campaign.AtomicGzipJsonLines(
                        too_many,
                        row_limit=1,
                        compressed_limit=1024,
                        uncompressed_limit=2 * payload_size,
                ) as writer:
                    writer.write(row)
                    writer.write(row)
            self.assertFalse(too_many.exists())
            self.assertFalse(
                too_many.with_name(too_many.name + ".tmp").exists())

            too_large = root / "too-large.jsonl.gz"
            with self.assertRaises(campaign.CampaignError):
                with campaign.AtomicGzipJsonLines(
                        too_large,
                        row_limit=1,
                        compressed_limit=1024,
                        uncompressed_limit=payload_size - 1,
                ) as writer:
                    writer.write(row)
            self.assertFalse(too_large.exists())

            compressed = root / "compressed.jsonl.gz"
            with self.assertRaises(campaign.CampaignError):
                with campaign.AtomicGzipJsonLines(
                        compressed,
                        row_limit=1,
                        compressed_limit=1,
                        uncompressed_limit=payload_size,
                ) as writer:
                    writer.write(row)
            self.assertFalse(compressed.exists())

    def test_result_cap_formula_uses_physical_cells(self):
        training = campaign.build_pre_cpu_jobs("training")[0]
        sealed = campaign.build_pre_cpu_jobs(
            "sealed", "pure8_s0m1_d3")[0]
        training_limits = campaign._result_evidence_byte_limits(training)
        sealed_limits = campaign._result_evidence_byte_limits(sealed)
        training_uncompressed = (
            campaign.RESULT_UNCOMPRESSED_FIXED_BYTES +
            campaign.TRAINING_REPLICATES *
            campaign.RESULT_UNCOMPRESSED_BYTES_PER_CELL)
        self.assertEqual(training_limits["uncompressed"],
                         training_uncompressed)
        self.assertEqual(
            training_limits["compressed"],
            training_uncompressed + (training_uncompressed >> 12) +
            (training_uncompressed >> 14) +
            (training_uncompressed >> 25) + 64,
        )
        self.assertEqual(
            sealed_limits["uncompressed"],
            campaign.RESULT_UNCOMPRESSED_FIXED_BYTES +
            campaign.SEALED_REPLICATES *
            campaign.RESULT_UNCOMPRESSED_BYTES_PER_CELL,
        )
        self.assertGreater(
            training_limits["uncompressed"],
            sealed_limits["uncompressed"],
        )


class _FakeRepairParser:
    REPAIRTIMING_PROTOCOL_V3 = campaign.REQUIRED_REPAIRTIMING_PROTOCOL
    REPAIRTIMING_SCHEMA_V3 = campaign.REQUIRED_REPAIRTIMING_SCHEMA
    REPAIRTIMING_ROWS_PER_CELL = 121
    REPAIRTIMING_STATS_FIELDS = ("block_xors", "build_ns")
    REPAIRTIMING_STRUCTURAL_STATS_FIELDS = ("block_xors",)
    REPAIRTIMING_PANELS = (
        "encoder_selector_forced",
        "encoder_selector_aa",
        "encoder_forced_aa",
        "encoder_selector_wh1",
        "encoder_wh1_aa",
        "decoder_feed_wh1",
        "decoder_feed_aa",
        "decoder_feed_wh1_aa",
        "decoder_full_wh1",
        "decoder_full_aa",
        "decoder_full_wh1_aa",
        "direct_dispatch",
        "direct_aa",
        "direct_dispatch_aa",
    )
    REPAIRTIMING_ATTEMPT_FIELDS = (
        "attempt",
        "probe_executed",
        "probe_result",
        "probe_stats_available",
        "probe_block_xors",
        "probe_build_ns",
        "real_executed",
        "real_result",
        "real_stats_available",
        "real_block_xors",
        "real_build_ns",
    )
    REPAIRTIMING_TIMING_FIELDS = (
        "timing_panel",
        "timing_panel_index",
        "timing_slot",
        "timing_pair",
        "timing_label",
        "timing_censor_reason",
        "timing_construct_result",
        "timing_result",
        "timing_recover_result",
        "timing_eligible",
        "timing_executed",
        "timing_stable",
        "timing_saturated",
        "timing_cpu_migrated",
        "timing_fault_contaminated",
        "timing_elapsed_ns",
    )
    REPAIRTIMING_CONTROL_FIELDS = (
        "forced_a_result",
        "forced_b_result",
        "forced_equal",
        "repair_direct_executed",
        "repair_direct_result",
        "repair_direct_witness_equal",
        "dispatch_direct_executed",
        "dispatch_direct_result",
        "dispatch_direct_witness_equal",
    )


def _synthetic_selector(root=7):
    return {
        "schema": "wirehair.wh2.repair-v1.selector-structural.v1",
        "arm_id": int("19cccf775ce0bf09", 16),
        "K": 2,
        "construction_root": root,
        "repair_policy_sha256": campaign.REPAIR_POLICY_SHA256,
        "selector_result": 0,
        "selected_attempt": 0,
        "attempts_executed": 1,
        "calls_executed": 1,
        "real_configuration_calls": 1,
        "structural_probe_calls": 0,
        "cap_exhausted": 0,
        "fatal_attempt_zero_mismatch": 0,
        "oom": 0,
        "committed": 1,
        "descriptor_sha256": "1" * 64,
        "attempts": [],
        "selector_structural_sha256": "2" * 64,
    }


def _synthetic_real(*, cross_left_ns=80, control_ns=100,
                    fault_panel=None, unavailable=False):
    attempt_rows = []
    for attempt in range(campaign.REPAIR_ATTEMPTS):
        attempt_rows.append([
            attempt,
            0, -1, 0, -1, -1,
            int(attempt == 0),
            0 if attempt == 0 else -1,
            0 if unavailable or attempt != 0 else 1,
            -1 if unavailable or attempt != 0 else 7,
            -1 if unavailable or attempt != 0 else 5,
        ])
    cross = {
        "encoder_selector_forced",
        "encoder_selector_wh1",
        "decoder_feed_wh1",
        "decoder_full_wh1",
        "direct_dispatch",
    }
    timing_rows = []
    order = "ABBABAAB"
    for panel_index, panel in enumerate(_FakeRepairParser.REPAIRTIMING_PANELS):
        for slot, label in enumerate(order):
            elapsed = (
                cross_left_ns
                if panel in cross and label == "A"
                else control_ns
            )
            timing_rows.append([
                panel,
                panel_index,
                slot,
                slot // 2,
                label,
                "none",
                0,
                0,
                0,
                1,
                1,
                1,
                0,
                0,
                int(panel == fault_panel and slot == 0),
                elapsed,
            ])
    role = {
        "encode_result": 0,
        "decode_construct_result": 0,
        "feed_result": 0,
        "recover_result": 0,
        "recovery_ok": 1,
        "encoded_symbols": 2,
        "received_symbols": 2,
        "overhead": 0,
        "payload_sha256": "3" * 64,
        "recovered_sha256": "4" * 64,
        "encode_class": "success",
        "feed_class": "success",
        "recover_class": "success",
        "outcome_class": "success",
    }
    return {
        "schema": "wirehair.wh2.repairtiming.real-witness.v3",
        "trace_sha256": "5" * 64,
        "message_sha256": "6" * 64,
        "descriptor_hex": "00" * 17,
        "roles": {
            name: dict(role)
            for name in ("raw", "repaired", "dispatch", "wh1")
        },
        "controls": {
            "forced_a_result": 0,
            "forced_b_result": 0,
            "forced_equal": 1,
            "repair_direct_executed": 1,
            "repair_direct_result": 0,
            "repair_direct_witness_equal": 1,
            "dispatch_direct_executed": 1,
            "dispatch_direct_result": 0,
            "dispatch_direct_witness_equal": 1,
        },
        "attempt_columns":
            list(_FakeRepairParser.REPAIRTIMING_ATTEMPT_FIELDS),
        "attempt_rows": attempt_rows,
        "timing_columns":
            list(_FakeRepairParser.REPAIRTIMING_TIMING_FIELDS),
        "timing_rows": timing_rows,
    }


class CompactAggregationTests(unittest.TestCase):
    def setUp(self):
        self.original = campaign._require_v3_parser
        campaign._require_v3_parser = lambda: _FakeRepairParser
        self.job = dict(campaign.build_pre_cpu_jobs("training")[0])
        self.job["K"] = 2
        self.job["arm"] = "pure8_s0m1_d3"
        self.job["arm_provisional_id"] = "19cccf775ce0bf09"
        self.job["schedule"] = "iid"
        self.job["job_id"] = "a" * 64

    def tearDown(self):
        campaign._require_v3_parser = self.original

    def compact(self, *, cells=4, fault_panel=None, unavailable=False):
        compact = campaign._new_compact_job_aggregate(self.job)
        for index in range(cells):
            campaign._accumulate_cell(
                compact,
                _synthetic_selector(index + 1),
                _synthetic_real(
                    fault_panel=fault_panel,
                    unavailable=unavailable,
                ),
            )
        compact["selector_set_sha256"] = "b" * 64
        compact["unique_selectors"] = cells
        compact["real_stream_sha256"] = "c" * 64
        campaign._finish_compact_job_aggregate(compact)
        return compact

    def test_compact_aggregation_preserves_raw_repaired_and_work(self):
        compact = self.compact()
        self.assertEqual(compact["cells"], 4)
        for role in ("raw", "repaired", "dispatch", "wh1"):
            self.assertEqual(
                compact["recovery"][role]["success_overhead"][0], 4)
            self.assertEqual(
                compact["recovery"][role]["outcome_classes"],
                {"success": 4},
            )
        work = compact["selector"]["work"]
        self.assertEqual(work["block_xors"], {
            "available": 4,
            "unavailable": 0,
            "histogram": {"7": 4},
        })
        self.assertEqual(work["build_ns"]["histogram"], {"5": 4})
        self.assertLess(
            compact["timing"]["encoder_selector_wh1"]["mean"], 0)
        self.assertEqual(
            compact["timing"]["encoder_selector_aa"]["mean"], 0)

    def test_unavailable_work_is_not_folded_into_zero(self):
        compact = self.compact(unavailable=True)
        work = compact["selector"]["work"]["block_xors"]
        self.assertEqual(work["available"], 0)
        self.assertEqual(work["unavailable"], 4)
        self.assertEqual(work["histogram"], {})

    def test_selector_metrics_count_each_key_once(self):
        compact = campaign._new_compact_job_aggregate(self.job)
        selector = _synthetic_selector(7)
        for index in range(5):
            campaign._accumulate_cell(
                compact,
                selector,
                _synthetic_real(),
                accumulate_selector=(index == 0),
            )
        compact["selector_set_sha256"] = "b" * 64
        compact["unique_selectors"] = 1
        compact["real_stream_sha256"] = "c" * 64
        campaign._finish_compact_job_aggregate(compact)
        self.assertEqual(compact["cells"], 5)
        self.assertEqual(compact["selector"]["observations"], 1)
        self.assertEqual(
            sum(compact["selector"]["result_codes"].values()), 1)
        self.assertEqual(
            compact["selector"]["work"]["block_xors"]["available"], 1)
        self.assertEqual(
            compact["selector"]["execution_counts"]["attempts_executed"],
            {"1": 1},
        )
        self.assertEqual(compact["selector"]["weak_constructions"], {
            "raw_attempt0": 0,
            "repaired_final": 0,
        })
        self.assertEqual(
            compact["recovery"]["repaired"]["success_overhead"][0], 5)

    def test_weak_constructions_are_deduplicated_from_physical_cells(self):
        compact = campaign._new_compact_job_aggregate(self.job)
        selector = _synthetic_selector(7)
        real = _synthetic_real()
        for role_name in ("raw", "repaired"):
            role = real["roles"][role_name]
            role.update({
                "encode_result": 3,
                "decode_construct_result": -1,
                "feed_result": -1,
                "recover_result": -1,
                "recovery_ok": 0,
                "encoded_symbols": 0,
                "received_symbols": 0,
                "overhead": -1,
                "payload_sha256": "not_applicable",
                "recovered_sha256": "not_applicable",
                "encode_class": "weak",
                "feed_class": "not_applicable",
                "recover_class": "not_applicable",
                "outcome_class": "weak",
            })
        for index in range(5):
            campaign._accumulate_cell(
                compact,
                selector,
                copy.deepcopy(real),
                accumulate_selector=(index == 0),
            )
        self.assertEqual(compact["selector"]["observations"], 1)
        self.assertEqual(compact["selector"]["weak_constructions"], {
            "raw_attempt0": 1,
            "repaired_final": 1,
        })
        self.assertEqual(
            compact["recovery"]["raw"]["outcome_classes"]["weak"], 5)
        self.assertEqual(compact["candidate_final_weak"], 5)

    def test_selected_attempt_max_ignores_zero_count_bins(self):
        histogram = {
            **{
                str(attempt): 0
                for attempt in range(campaign.REPAIR_ATTEMPTS)
            },
            "none": 0,
        }
        histogram["0"] = 10
        selected = campaign._selected_attempt_statistics(histogram)
        self.assertEqual(selected["selected_statistics"]["mean"], 0.0)
        self.assertEqual(selected["selected_statistics"]["p99"], 0)
        self.assertEqual(selected["selected_statistics"]["max"], 0)
        histogram["3"] = 1
        selected = campaign._selected_attempt_statistics(histogram)
        self.assertEqual(selected["selected_statistics"]["max"], 3)

    def test_selector_execution_count_accounting_is_enforced(self):
        compact = campaign._new_compact_job_aggregate(self.job)
        selector = _synthetic_selector()
        selector["calls_executed"] = 2
        with self.assertRaisesRegex(
                campaign.CampaignError, "accounting is contradictory"):
            campaign._accumulate_cell(
                compact, selector, _synthetic_real())

        compact = self.compact()
        bucket = campaign._new_phase_bucket("invalid-count")
        campaign._merge_compact_into_bucket(
            bucket,
            compact,
            {
                "valid_rows": 2,
                "cpu_tctl_max_c": 65.0,
                "dimm_max_c": 52.0,
                "edac_ce_max": 0,
                "edac_ue_max": 0,
            },
        )
        observations = bucket["selector"]["observations"]
        bucket["selector"]["execution_counts"]["attempts_executed"] = {
            str(campaign.REPAIR_ATTEMPTS + 1): observations,
        }
        with self.assertRaises(campaign.CampaignError):
            campaign._finalize_phase_bucket(bucket)

    def test_cell_aggregate_deduplicates_fixed_production_root(self):
        job = next(
            copy.deepcopy(item)
            for item in campaign.build_pre_cpu_jobs(
                "sealed", "pure8_s0m1_d3")
            if item["lane"] == "production" and
            item["K"] == 2 and item["schedule"] == campaign.SCHEDULES[0]
        )
        job["replicates"] = 3
        root = campaign._expected_root(job, 0)
        selector = _synthetic_selector(root)
        real = _synthetic_real()
        cells = tuple(
            SimpleNamespace(
                replicate=replicate,
                construction_root=root,
                loss_seed=campaign._expected_loss(job, replicate),
                selector_key=(
                    job["arm_provisional_id"], job["K"], root),
                selector=selector,
                real=real,
            )
            for replicate in range(job["replicates"])
        )
        with (
            mock.patch.object(
                _FakeRepairParser,
                "selector_projection",
                side_effect=lambda cell: copy.deepcopy(cell.selector),
                create=True,
            ),
            mock.patch.object(
                _FakeRepairParser,
                "real_projection",
                side_effect=lambda cell: copy.deepcopy(cell.real),
                create=True,
            ),
        ):
            compact = campaign._cell_aggregate(
                _FakeRepairParser,
                SimpleNamespace(cells=cells),
                job,
            )
        self.assertEqual(compact["cells"], 3)
        self.assertEqual(compact["unique_selectors"], 1)
        self.assertEqual(compact["selector"]["observations"], 1)
        self.assertEqual(
            compact["recovery"]["repaired"]["success_overhead"][0], 3)

    def test_selector_roster_totals_use_only_first_schedule(self):
        for phase, survivor, expected in (
            ("training", None, 152064),
            ("sealed", "pure8_s0m1_d3", 25443),
        ):
            total = 0
            production = 0
            for job in campaign.build_pre_cpu_jobs(phase, survivor):
                if job["schedule"] != campaign.SCHEDULES[0]:
                    continue
                unique = (
                    1 if job["lane"] == "production"
                    else job["replicates"]
                )
                total += unique
                if job["lane"] == "production":
                    production += unique
            with self.subTest(phase=phase):
                self.assertEqual(total, expected)
                if phase == "sealed":
                    self.assertEqual(production, len(campaign.K_VALUES))

    def test_uncommitted_forced_control_and_independent_dispatch(self):
        compact = campaign._new_compact_job_aggregate(self.job)
        selector = _synthetic_selector()
        selector["committed"] = 0
        selector["selected_attempt"] = -1
        real = _synthetic_real()
        real["controls"].update({
            "forced_a_result": -1,
            "forced_b_result": -1,
            "forced_equal": 0,
            "repair_direct_executed": 0,
            "repair_direct_result": -1,
            "repair_direct_witness_equal": 0,
        })
        campaign._accumulate_cell(compact, selector, real)
        self.assertEqual(compact["forced_control_failure"], 0)
        self.assertEqual(compact["selector"]["uncommitted"], 1)

        failed = campaign._new_compact_job_aggregate(self.job)
        real["controls"]["dispatch_direct_witness_equal"] = 0
        campaign._accumulate_cell(failed, selector, real)
        self.assertEqual(failed["forced_control_failure"], 1)

    def test_oom_and_result_classes_are_distinct(self):
        for code, runtime, oom in (
            (7, 0, 0),
            (8, 1, 0),
            (9, 0, 1),
        ):
            with self.subTest(code=code):
                compact = campaign._new_compact_job_aggregate(self.job)
                real = _synthetic_real()
                real["roles"]["repaired"]["feed_result"] = code
                campaign._accumulate_cell(
                    compact, _synthetic_selector(), real)
                self.assertEqual(
                    compact["candidate_runtime_error"], runtime)
                self.assertEqual(compact["codec_oom"], oom)

        control = campaign._new_compact_job_aggregate(self.job)
        real = _synthetic_real()
        real["controls"]["dispatch_direct_result"] = 9
        campaign._accumulate_cell(
            control, _synthetic_selector(), real)
        self.assertEqual(control["codec_oom"], 1)
        self.assertEqual(control["candidate_runtime_error"], 0)

        control_error = campaign._new_compact_job_aggregate(self.job)
        real = _synthetic_real()
        real["controls"]["dispatch_direct_result"] = 8
        campaign._accumulate_cell(
            control_error, _synthetic_selector(), real)
        self.assertEqual(control_error["codec_oom"], 0)
        self.assertEqual(
            control_error["candidate_runtime_error"], 1)

    def test_nonrepairable_regression_uses_either_matched_control(self):
        cases = (
            ("dispatch-only", 1, 0, "need_more", 1),
            ("wh1-only", 0, 1, "need_more", 1),
            ("both", 1, 1, "error", 1),
            ("shared", 0, 0, "error", 0),
            ("weak", 1, 1, "weak", 0),
            ("code-7", 1, 0, "need_more", 1),
        )
        for name, dispatch_ok, wh1_ok, outcome, expected in cases:
            with self.subTest(name=name):
                compact = campaign._new_compact_job_aggregate(self.job)
                real = _synthetic_real()
                repaired = real["roles"]["repaired"]
                repaired["recovery_ok"] = 0
                repaired["overhead"] = -1
                repaired["outcome_class"] = outcome
                if name == "code-7":
                    repaired["feed_result"] = 7
                for role_name, recovered in (
                    ("dispatch", dispatch_ok),
                    ("wh1", wh1_ok),
                ):
                    role = real["roles"][role_name]
                    role["recovery_ok"] = recovered
                    role["overhead"] = 0 if recovered else -1
                campaign._accumulate_cell(
                    compact, _synthetic_selector(), real)
                self.assertEqual(
                    compact["candidate_nonrepairable_regression"],
                    expected,
                )

    def test_fault_contamination_censors_whole_cell_panel(self):
        compact = self.compact(fault_panel="direct_dispatch")
        panel = compact["timing"]["direct_dispatch"]
        self.assertEqual(panel["eligible_cells"], 0)
        self.assertEqual(panel["censored_cells"], 4)
        self.assertIsNone(panel["mean"])

    def test_timing_variance_avoids_large_offset_cancellation(self):
        accumulator = {
            "eligible_cells": 0,
            "censored_cells": 0,
            "mean": 0.0,
            "m2": 0.0,
        }
        values = [
            1_000_000_000.0 + delta
            for delta in (-0.003, -0.001, 0.001, 0.003)
        ]
        for value in values:
            campaign._timing_accumulator_add(accumulator, value)
        finished = campaign._finish_timing_accumulator(accumulator)
        centered = math.fsum(values) / len(values)
        expected_m2 = math.fsum(
            (value - centered) ** 2 for value in values)
        self.assertGreater(finished["m2"], 0.0)
        self.assertAlmostEqual(
            finished["m2"], expected_m2, delta=1e-9)
        self.assertGreater(
            finished["ci_high"] - finished["mean"], 0.0)
        self.assertGreater(finished["floor"], 0.0)

    def test_bucket_finalization_and_training_ranking(self):
        pure8_bucket = campaign._new_phase_bucket("pure8")
        pure9_bucket = campaign._new_phase_bucket("pure9")
        compact = self.compact()
        thermal = {
            "valid_rows": 2,
            "cpu_tctl_max_c": 65.0,
            "dimm_max_c": 52.0,
            "edac_ce_max": 0,
            "edac_ue_max": 0,
        }
        for index in range(4):
            item = copy.deepcopy(compact)
            item["job_id"] = f"{index:064x}"
            campaign._merge_compact_into_bucket(
                pure8_bucket, item, thermal)
            campaign._merge_compact_into_bucket(
                pure9_bucket, item, thermal)
        pure8 = campaign._finalize_phase_bucket(pure8_bucket)
        pure9 = campaign._finalize_phase_bucket(pure9_bucket)
        eligibility = campaign._arm_eligibility(pure8)
        self.assertTrue(eligibility["eligible"])
        self.assertIn(
            "full_encoder_selector_forced", pure8["timing"])
        self.assertNotIn(
            "full_encoder_selector_forced",
            eligibility["timing_pass"],
        )
        execution = pure8["selector"]["execution_counts"]
        self.assertEqual(
            execution["attempts_executed"]["histogram"], {"1": 16})
        self.assertEqual(
            execution["attempts_executed"]["statistics"]["mean"], 1.0)
        self.assertEqual(
            execution["attempts_executed"]["statistics"]["max"], 1)
        self.assertEqual(
            execution["structural_probe_calls"]["histogram"], {"0": 16})
        # Decoder solve is the user's primary optimization target.  A faster
        # encoder must not outrank a lower decoder-feed cost once both arms
        # pass every recovery and timing gate.
        pure8["timing"]["decoder_feed_candidate_wh1"]["mean"] = -0.1
        pure9["timing"]["decoder_feed_candidate_wh1"]["mean"] = -0.2
        pure8["timing"]["full_encoder_candidate_wh1"]["mean"] = -0.3
        pure9["timing"]["full_encoder_candidate_wh1"]["mean"] = -0.1
        decision = campaign._derive_training_decision({
            "pure8_s0m1_d3": pure8,
            "pure9_s0m1_d3": pure9,
        })
        self.assertEqual(decision["status"], "winner")
        self.assertEqual(
            decision["selected_survivor"], "pure9_s0m1_d3")
        self.assertEqual(
            decision["candidates"]["pure8_s0m1_d3"][
                "fresh_weak_constructions"
            ],
            {"raw_attempt0": 0, "repaired_final": 0},
        )
        pure8["selector"]["cap_exhausted"] = 1
        pure9["selector"]["cap_exhausted"] = 1
        decision = campaign._derive_training_decision({
            "pure8_s0m1_d3": pure8,
            "pure9_s0m1_d3": pure9,
        })
        self.assertEqual(decision["status"], "no-survivor")
        self.assertIsNone(decision["selected_survivor"])

    def test_sealed_confirmation_requires_each_prebound_lane(self):
        bucket = campaign._new_phase_bucket("sealed")
        compact = self.compact()
        thermal = {
            "valid_rows": 2,
            "cpu_tctl_max_c": 65.0,
            "dimm_max_c": 52.0,
            "edac_ce_max": 0,
            "edac_ue_max": 0,
        }
        for index in range(4):
            item = copy.deepcopy(compact)
            item["job_id"] = f"{index:064x}"
            campaign._merge_compact_into_bucket(bucket, item, thermal)
        eligible = campaign._finalize_phase_bucket(bucket)
        random_lane = copy.deepcopy(eligible)
        production_lane = copy.deepcopy(eligible)
        production_lane["timing"][
            "decoder_feed_candidate_wh1"
        ]["pass"] = False
        confirmation = campaign._derive_sealed_confirmation(
            "pure8_s0m1_d3",
            eligible,
            {
                "random": random_lane,
                "production": production_lane,
            },
        )
        self.assertTrue(confirmation["eligibility"]["eligible"])
        self.assertTrue(
            confirmation["lane_eligibility"]["random"]["eligible"])
        self.assertFalse(
            confirmation["lane_eligibility"]["production"]["eligible"])
        self.assertEqual(confirmation["status"], "failed")
        self.assertEqual(confirmation["fallback"], "forbidden")

    def test_selector_dedup_rejects_schedule_drift(self):
        compact = self.compact()
        groups = {}
        controls = {}
        campaign._validate_cross_job_invariants(
            compact, groups, controls)
        changed = copy.deepcopy(compact)
        changed["schedule"] = "burst"
        changed["selector_set_sha256"] = "d" * 64
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_cross_job_invariants(
                changed, groups, controls)

    def test_recovery_outcome_class_ledger_must_cover_every_cell(self):
        compact = self.compact()
        receipt = copy.deepcopy(compact["recovery"]["repaired"])
        receipt["outcome_classes"]["ghost"] = 1
        with self.assertRaises(campaign.CampaignError):
            campaign._recovery_metrics(receipt, compact["cells"])

    def test_timing_row_reordering_is_rejected(self):
        compact = campaign._new_compact_job_aggregate(self.job)
        real = _synthetic_real()
        real["timing_rows"][0][2] = 1
        with self.assertRaises(campaign.CampaignError):
            campaign._accumulate_cell(
                compact, _synthetic_selector(), real)


class RollingSchedulerTests(unittest.TestCase):
    @staticmethod
    def manifest():
        return {
            "workers": 2,
            "worker_cpus": [0, 1],
            "runtime_bindings": {"synthetic": "binding"},
            "assignments": [
                {
                    "order": order,
                    "cpu": order % 2,
                    "job_file": f"jobs/{order:04d}.json",
                    "output": f"results/{order:04d}.json.gz",
                    "log": f"logs/{order:04d}.log",
                }
                for order in range(4)
            ],
        }

    @staticmethod
    def fake_popen_factory(
            durations, returncodes, launched, handles,
            affinity_pins=None):
        next_pid = [1000]

        class FakeProcess:
            def __init__(self, order):
                self.order = order
                self.remaining = durations[order]
                self.returncode = None
                self.pid = next_pid[0]
                next_pid[0] += 1

            def poll(self):
                if self.returncode is not None:
                    return self.returncode
                self.remaining -= 1
                if self.remaining <= 0:
                    self.returncode = returncodes.get(self.order, 0)
                return self.returncode

            def wait(self, timeout=None):
                del timeout
                return self.returncode

        def popen(command, **kwargs):
            order = int(Path(
                command[command.index("--job-file") + 1]
            ).stem)
            preexec = kwargs.get("preexec_fn")
            if (
                getattr(preexec, "func", None) is not
                    campaign._pin_worker_before_exec or
                getattr(preexec, "args", None) != (order % 2,) or
                kwargs.get("start_new_session") is not True
            ):
                raise AssertionError(
                    "worker was not pinned before source import")
            if affinity_pins is not None:
                affinity_pins.append((order, preexec.args[0]))
            launched.append(order)
            handles.append(kwargs["stdout"])
            if sum(not handle.closed for handle in handles) > 2:
                raise AssertionError("scheduler retained too many logs")
            return FakeProcess(order)

        return popen

    def test_fast_cpu_refills_before_slow_cpu_finishes(self):
        manifest = self.manifest()
        launched = []
        handles = []
        affinity_pins = []
        popen = self.fake_popen_factory(
            {0: 6, 1: 1, 2: 1, 3: 1},
            {},
            launched,
            handles,
            affinity_pins,
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "logs").mkdir()
            with (
                mock.patch.object(
                    campaign.subprocess, "Popen", side_effect=popen
                ),
                mock.patch.object(
                    campaign, "check_runtime_bindings"
                ),
                mock.patch.object(
                    campaign, "_process_group_exists", return_value=False
                ),
                mock.patch.object(campaign.time, "sleep"),
            ):
                monitor = campaign._run_rolling_workers(
                    directory, manifest)
        self.assertEqual(sorted(launched), [0, 1, 2, 3])
        self.assertEqual(sorted(affinity_pins), [
            (0, 0), (1, 1), (2, 0), (3, 1),
        ])
        by_cpu = {
            receipt["cpu"]: receipt
            for receipt in monitor["cpu_queues"]
        }
        self.assertEqual(by_cpu[0]["orders"], [0, 2])
        self.assertEqual(by_cpu[1]["orders"], [1, 3])
        self.assertLess(
            by_cpu[1]["launches"][1]["started_monotonic_ns"],
            by_cpu[0]["launches"][0]["finished_monotonic_ns"],
        )
        self.assertTrue(all(handle.closed for handle in handles))
        self.assertEqual(monitor["pending_idle_checks"], 0)
        campaign._validate_runtime_monitor(monitor, manifest)

        changed = copy.deepcopy(monitor)
        changed["pending_idle_checks"] = 1
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_runtime_monitor(changed, manifest)

    def test_first_failure_stops_refill_and_cleans_all_groups(self):
        manifest = self.manifest()
        launched = []
        handles = []
        cleaned = []
        popen = self.fake_popen_factory(
            {0: 10, 1: 1, 2: 1, 3: 1},
            {1: 17},
            launched,
            handles,
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "logs").mkdir()
            with (
                mock.patch.object(
                    campaign.subprocess, "Popen", side_effect=popen
                ),
                mock.patch.object(
                    campaign, "check_runtime_bindings"
                ),
                mock.patch.object(
                    campaign, "_process_group_exists", return_value=False
                ),
                mock.patch.object(
                    campaign,
                    "_kill_process_groups_with_grace",
                    side_effect=lambda processes, unused_grace:
                        cleaned.extend(processes),
                ),
                self.assertRaises(campaign.CampaignError),
            ):
                campaign._run_rolling_workers(directory, manifest)
        self.assertEqual(launched, [0, 1])
        self.assertEqual([process.order for process in cleaned], [0])
        self.assertTrue(all(handle.closed for handle in handles))

    def test_successful_leader_with_descendants_fails_closed(self):
        manifest = self.manifest()
        launched = []
        handles = []
        cleaned = []
        popen = self.fake_popen_factory(
            {0: 1, 1: 1, 2: 1, 3: 1},
            {},
            launched,
            handles,
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "logs").mkdir()
            with (
                mock.patch.object(
                    campaign.subprocess, "Popen", side_effect=popen
                ),
                mock.patch.object(
                    campaign, "check_runtime_bindings"
                ),
                mock.patch.object(
                    campaign,
                    "_process_group_exists",
                    side_effect=lambda pid: pid == 1000,
                ),
                mock.patch.object(
                    campaign,
                    "_kill_process_groups_with_grace",
                    side_effect=lambda processes, unused_grace:
                        cleaned.extend(processes),
                ),
                self.assertRaisesRegex(
                    campaign.CampaignError, "left descendants"
                ),
            ):
                campaign._run_rolling_workers(directory, manifest)
        self.assertEqual(launched, [0, 1])
        self.assertEqual([process.order for process in cleaned], [0])
        self.assertTrue(all(handle.closed for handle in handles))

    def test_parent_abort_cleans_every_active_group(self):
        manifest = self.manifest()
        launched = []
        handles = []
        cleaned = []
        safe_points = [0]
        popen = self.fake_popen_factory(
            {0: 10, 1: 10, 2: 10, 3: 10},
            {},
            launched,
            handles,
        )

        def safe_point():
            safe_points[0] += 1
            if safe_points[0] == 3:
                raise campaign.CampaignError("synthetic parent signal")

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "logs").mkdir()
            with (
                mock.patch.object(
                    campaign.subprocess, "Popen", side_effect=popen
                ),
                mock.patch.object(
                    campaign, "check_runtime_bindings"
                ),
                mock.patch.object(
                    campaign,
                    "_parent_signal_safe_point",
                    side_effect=safe_point,
                ),
                mock.patch.object(
                    campaign,
                    "_kill_process_groups_with_grace",
                    side_effect=lambda processes, unused_grace:
                        cleaned.extend(processes),
                ),
                self.assertRaises(campaign.CampaignError),
            ):
                campaign._run_rolling_workers(directory, manifest)
        self.assertEqual(launched, [0, 1])
        self.assertEqual(len(cleaned), 2)
        self.assertTrue(all(handle.closed for handle in handles))

    def test_preexec_affinity_helper_pins_exact_cpu(self):
        with (
            mock.patch.object(campaign.os, "sched_setaffinity") as setaffinity,
            mock.patch.object(
                campaign.os, "sched_getaffinity", return_value={3}
            ),
        ):
            campaign._pin_worker_before_exec(3)
        setaffinity.assert_called_once_with(0, {3})

    def test_group_cleanup_returns_promptly_after_sigterm(self):
        process = subprocess.Popen(
            ["/bin/sleep", "30"],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        started = time.monotonic()
        try:
            campaign._kill_process_groups_with_grace([process], 2.0)
        finally:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        self.assertLess(time.monotonic() - started, 1.0)
        self.assertEqual(process.returncode, -signal.SIGTERM)
        self.assertFalse(campaign._process_group_exists(process.pid))

    def test_group_cleanup_rejects_survivors_after_sigkill(self):
        class StubbornProcess:
            pid = 123456789

            @staticmethod
            def poll():
                return None

            @staticmethod
            def wait(timeout=None):
                raise subprocess.TimeoutExpired("stubborn", timeout)

            @staticmethod
            def terminate():
                pass

            @staticmethod
            def kill():
                pass

        with (
            mock.patch.object(
                campaign, "_process_group_exists", return_value=True
            ),
            mock.patch.object(campaign.os, "killpg"),
            mock.patch.object(
                campaign, "WORKER_KILL_GRACE_SECONDS", 0.01
            ),
            self.assertRaisesRegex(
                campaign.CampaignError, "left surviving"
            ),
        ):
            campaign._kill_process_groups_with_grace(
                [StubbornProcess()], 0.01)

    def test_group_cleanup_never_reacquires_a_gone_pgid(self):
        class FinishedProcess:
            pid = 123456788
            returncode = 0

            @staticmethod
            def poll():
                return 0

            @staticmethod
            def wait(timeout=None):
                del timeout
                return 0

        existence = mock.Mock(side_effect=(True, False, True))
        with (
            mock.patch.object(
                campaign, "_process_group_exists", existence
            ),
            mock.patch.object(campaign.os, "killpg") as killpg,
        ):
            campaign._kill_process_groups_with_grace(
                [FinishedProcess()], 0)
        self.assertEqual(existence.call_count, 2)
        killpg.assert_called_once_with(
            FinishedProcess.pid, signal.SIGTERM)


class ParallelReplayTests(unittest.TestCase):
    @staticmethod
    def manifest(count):
        return {
            "assignments": [
                {"order": index} for index in range(count)
            ],
            "pre_cpu_jobs": [
                {"job_id": f"{index:064x}"} for index in range(count)
            ],
        }

    def test_slow_first_replay_is_ordered_and_window_bounded(self):
        global _PARALLEL_RELEASED
        global _PARALLEL_STARTED_BEFORE_RELEASE
        context = multiprocessing.get_context("fork")
        _PARALLEL_RELEASED = context.Event()
        _PARALLEL_STARTED_BEFORE_RELEASE = context.Value("i", 0)
        try:
            with mock.patch.object(
                    campaign,
                    "_replay_process_task",
                    _parallel_slow_first_task,
            ):
                values = list(campaign._ordered_replayed_results(
                    Path("."),
                    self.manifest(100),
                    "a" * 64,
                    4,
                ))
            self.assertEqual(
                [value["index"] for value in values],
                list(range(100)),
            )
            self.assertLessEqual(
                _PARALLEL_STARTED_BEFORE_RELEASE.value + 1,
                campaign.REPLAY_QUEUE_MAX_ITEMS,
            )
        finally:
            _PARALLEL_RELEASED = None
            _PARALLEL_STARTED_BEFORE_RELEASE = None

    def test_worker_error_reaps_every_process(self):
        stopped = []
        original = campaign._stop_replay_processes

        def checked_stop(processes):
            original(processes)
            stopped.extend(processes)
            self.assertTrue(all(
                not process.is_alive() for process in processes))

        with (
            mock.patch.object(
                campaign, "_replay_process_task", _parallel_error_task
            ),
            mock.patch.object(
                campaign, "_stop_replay_processes", checked_stop
            ),
            self.assertRaisesRegex(
                campaign.CampaignError, "synthetic replay failure"
            ),
        ):
            list(campaign._ordered_replayed_results(
                Path("."),
                self.manifest(8),
                "b" * 64,
                4,
            ))
        self.assertEqual(len(stopped), 4)

    def test_stubborn_worker_is_killed_after_term_grace(self):
        context = multiprocessing.get_context("fork")
        ready = context.Event()
        process = context.Process(
            target=_ignore_sigterm_forever, args=(ready,))
        process.start()
        self.assertTrue(ready.wait(2.0))
        campaign._stop_replay_processes([process])
        self.assertFalse(process.is_alive())
        self.assertEqual(process.exitcode, -signal.SIGKILL)


if __name__ == "__main__":
    unittest.main()
