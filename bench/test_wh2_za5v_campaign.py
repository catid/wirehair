#!/usr/bin/env python3
"""Fast contract and aggregation tests for wh2_za5v_campaign.py."""

import gzip
import hashlib
import json
import os
from contextlib import nullcontext
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
from types import SimpleNamespace
import unittest
from unittest import mock

BENCH = Path(__file__).resolve().parent
if str(BENCH) not in sys.path:
    sys.path.insert(0, str(BENCH))

import wh2_za5v_campaign as campaign
import wh2_za5v_parallel_verify as parallel_verify


def adversarial_replay_failure(
        directory, unused_manifest, unused_manifest_sha256,
        assignment, job):
    marker = Path(job["_test_marker"])
    marker.parent.mkdir(parents=True, exist_ok=True)
    (marker.parent / f"worker-{os.getpid()}").touch()
    if assignment["order"] == 0:
        marker.touch()
        time.sleep(30)
        raise AssertionError("sleeping replay unexpectedly completed")
    deadline = time.monotonic() + 5.0
    while not marker.exists():
        if time.monotonic() >= deadline:
            raise AssertionError("slow replay did not start")
        time.sleep(0.005)
    raise campaign.CampaignError("adversarial fast failure")


def adversarial_wrapper_task(task):
    assignment, job = task
    marker = Path(job["_test_marker"])
    marker.parent.mkdir(parents=True, exist_ok=True)
    (marker.parent / f"worker-{os.getpid()}").touch()
    if assignment["order"] == 0:
        marker.touch()
        time.sleep(30)
        raise AssertionError("sleeping replay unexpectedly completed")
    deadline = time.monotonic() + 5.0
    while not marker.exists():
        if time.monotonic() >= deadline:
            raise AssertionError("slow replay did not start")
        time.sleep(0.005)
    return {
        "status": "error",
        "order": assignment["order"],
        "job_id": job["job_id"],
        "exception_kind": "campaign",
        "message": "adversarial fast failure",
    }


def adversarial_wrapper_initializer(*unused_arguments):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGHUP, signal.SIG_DFL)


def synthetic_runtime_bindings():
    file_binding = {
        "path": "/tmp/frozen-runtime",
        "device": 1,
        "inode": 2,
        "mode": 0o100755,
        "size": 3,
        "mtime_ns": 4,
        "ctime_ns": 5,
        "sha256": "0" * 64,
    }
    return {
        name: dict(file_binding)
        for name in ("benchmark", "parser", "context_tool", "runner")
    } | {
        "thermal": {
            "path": "/tmp/frozen-thermal",
            "device": 6,
            "inode": 7,
            "mode": 0o100644,
        },
    }


def parallel_replay_fixture(directory, *, fault=None):
    """Write a small, fully replayable four-candidate selection population."""
    from test_band_timing_tools import (
        band_context,
        bandtiming_stdout,
        parse_kwargs,
    )

    directory = Path(directory)
    for name in ("jobs", "results", "logs"):
        (directory / name).mkdir(parents=True)
    thermal = directory / "thermal.csv"
    thermal.write_text("fixture\n", encoding="ascii")
    bindings = campaign.runtime_bindings(sys.executable, thermal)
    jobs = campaign.build_pre_cpu_jobs("selection")
    assignments = [
        {
            "order": index,
            "job_id": job["job_id"],
            "cpu": 2,
            "job_file": f"jobs/{index:04d}.json",
            "output": f"results/{index:04d}.json.gz",
            "log": f"logs/{index:04d}.log",
        }
        for index, job in enumerate(jobs)
    ]
    manifest = {
        "bandtiming_protocol": campaign.REQUIRED_BANDTIMING_PROTOCOL,
        "phase": "selection",
        "survivor": None,
        "workers": 2,
        "worker_cpus": [2],
        "runtime_bindings": bindings,
        "selection_contract": None,
        "pre_cpu_jobs": list(jobs),
        "assignments": assignments,
        "frozen_roster_sha256": "a" * 64,
        "pre_cpu_job_list_sha256":
            campaign.canonical_sha256(list(jobs)),
    }
    manifest_sha256 = campaign.canonical_sha256(manifest)
    campaign._write_job_files(
        directory, manifest, manifest_sha256)

    for index, (assignment, job) in enumerate(zip(assignments, jobs)):
        context = band_context(cache_state=job["cache_state"])
        context["bound"]["thermal_device"] = \
            bindings["thermal"]["device"]
        context["bound"]["thermal_inode"] = \
            bindings["thermal"]["inode"]
        fixture = {
            "candidate": campaign.descriptor_from_job(job),
            "schedule": job["schedule"],
            "systematic_cache": job["systematic_cache"],
            "cache_state": job["cache_state"],
            "context": context,
            "block_count": job["K"],
            "block_bytes": job["block_bytes"],
            "replicates": job["replicates"],
            "warmups": job["warmup_replicates"],
            "construction_seed": job["construction_seed_base"],
            "loss_seed": job["loss_seed_base"],
            "inner_reps": job["inner_reps"],
            "max_overhead": job["max_overhead"],
            "required_margin": job["required_margin"],
            "loss": job["loss"],
        }
        measurement = campaign.band.parse_bandtiming_output(
            bandtiming_stdout(
                schema=campaign.band.BANDTIMING_SCHEMA_V2,
                **fixture,
            ),
            _protocol=campaign.band.BANDTIMING_PROTOCOL_V2,
            **parse_kwargs(**fixture),
        )
        receipt = measurement.as_dict()
        if fault == "derived" and index == 3:
            receipt["replicates"][0]["loss_seed"] ^= 1
        wall_started = 100 + index * 2
        wall_finished = wall_started + 1
        if fault == "envelope" and index == 3:
            wall_finished = wall_started
        campaign.atomic_gzip_json(
            directory / assignment["output"],
            {
                "schema": campaign.RESULT_SCHEMA,
                "manifest_sha256": manifest_sha256,
                "job": job,
                "assignment": assignment,
                "wall_started_unix_ns": wall_started,
                "wall_finished_unix_ns": wall_finished,
                "runtime_bindings_before": bindings,
                "runtime_bindings_after": bindings,
                "receipt": receipt,
            },
        )
        (directory / assignment["log"]).touch()
    return manifest


def synthetic_clean_selection_blockers(candidate):
    return {
        "schema": campaign.PRODUCTION_BLOCKER_SCHEMA,
        "candidate": candidate,
        "candidate_weak_constructions": 0,
        "candidate_error_outcomes": {
            "constructor": 0,
            "encoder": 0,
            "decoder": 0,
            "direct": 0,
            "total": 0,
        },
        "candidate_control_panels": {
            panel: {
                "recovery_regressions": 0,
                "recovery_improvements": 0,
                "both_failures": 0,
                "candidate_constructor_weak_regressions": 0,
                "candidate_nonrepairable_regressions": 0,
            }
            for panel, unused_scope, unused_control
            in campaign.CANDIDATE_CONTROL_PANELS
        },
        "repairable_seed_blocker": False,
        "nonrepairable_reliability_blocker": False,
    }


def process_group_absent_after_term(unused_process_group, sent_signal):
    if sent_signal == 0:
        raise ProcessLookupError


class FrozenPopulationTests(unittest.TestCase):
    def test_exact_dimensions_and_identical_candidate_cells(self):
        validation = campaign.validate_frozen_contract()
        self.assertEqual(validation["selection_jobs"], 2376)
        self.assertEqual(validation["holdout_jobs_per_survivor"], 594)
        self.assertEqual(
            validation["matched_selection_cells_per_candidate"], 456192)

        jobs = campaign.build_pre_cpu_jobs("selection")
        self.assertEqual(len(jobs), 2376)
        self.assertEqual(len({job["job_id"] for job in jobs}), 2376)
        populations = {}
        for job in jobs:
            population = populations.setdefault(job["candidate"], set())
            population.add((
                job["K"],
                job["schedule"],
                job["construction_seed_base"],
                job["loss_seed_base"],
                job["replicates"],
            ))
        reference = next(iter(populations.values()))
        self.assertTrue(all(value == reference
                            for value in populations.values()))
        self.assertEqual(len(reference), 99 * 6)

        for survivor in campaign.CANDIDATE_BY_NAME:
            holdout = campaign.build_pre_cpu_jobs("holdout", survivor)
            self.assertEqual(len(holdout), 594)
            self.assertEqual(len({job["job_id"] for job in holdout}), 594)
            self.assertEqual({job["candidate"] for job in holdout}, {survivor})

    def test_seed_sets_are_exhaustively_disjoint(self):
        proof = campaign.seed_disjointness_proof()
        self.assertEqual(proof["construction"]["intersection_count"], 0)
        self.assertEqual(proof["loss"]["intersection_count"], 0)
        self.assertEqual(proof["construction"]["selection_count"], 768)
        self.assertEqual(proof["construction"]["holdout_count"], 256)
        self.assertEqual(proof["loss"]["selection_count"], 768)
        self.assertEqual(proof["loss"]["holdout_count"], 256)

        selection_construction = {
            campaign.construction_seed(
                campaign.SELECTION_CONSTRUCTION_BASE, replicate)
            for replicate in range(campaign.SELECTION_REPLICATES)
        }
        holdout_construction = {
            campaign.construction_seed(
                campaign.HOLDOUT_CONSTRUCTION_BASE, replicate)
            for replicate in range(campaign.HOLDOUT_REPLICATES)
        }
        selection_loss = {
            campaign.loss_seed(campaign.SELECTION_LOSS_BASE, replicate)
            for replicate in range(campaign.SELECTION_REPLICATES)
        }
        holdout_loss = {
            campaign.loss_seed(campaign.HOLDOUT_LOSS_BASE, replicate)
            for replicate in range(campaign.HOLDOUT_REPLICATES)
        }
        self.assertFalse(selection_construction & holdout_construction)
        self.assertFalse(selection_loss & holdout_loss)

    def test_deterministic_frozen_hashes(self):
        self.assertEqual(
            campaign.frozen_roster_sha256(),
            "d7f032731d3f5b1a52d37cb27cc4236e92ef8be46efa5445934c8356288421fc",
        )
        self.assertEqual(
            campaign.pre_cpu_job_list_hash("selection"),
            "021dcb0607d83e5ce5da1d4ebbafbab43f1dfcef34209b64326a877fe97e2b63",
        )
        self.assertEqual(
            campaign.validate_frozen_contract()[
                "matched_selection_cell_set_sha256"],
            "6cc049811a0826e1be9fa750497ba03028fd263a2dc3328f9409b8cbaf9a8cdf",
        )
        self.assertEqual(campaign.holdout_hashes(), {
            "pure8_s0m1_d3":
                "621b29e4c897bc47dd4497dc2e990fedbf2236333906d36df77b4ace6c4560b5",
            "pure8_s0_d3":
                "bb15989bcca31a956100cdae64e2ee281d36ca341eca814f0b49b624dd293e4f",
            "pure8_s0m1_d5":
                "04b32e54498252315c74a1d64270d1669cc89fa308ebbbff230353e6bd4466e9",
            "pure9_s0m1_d3":
                "26b0b72a4868aafac50b0f1071e5fe856631d5aa69fe4736bdd41a796511a672",
        })

    def test_dry_plans_are_result_free_and_holdout_is_explicit(self):
        selection = campaign.make_plan("selection")
        self.assertEqual(selection["pre_cpu_job_count"], 2376)
        self.assertEqual(
            selection["bandtiming_protocol"],
            campaign.band.BANDTIMING_PROTOCOL_V2,
        )
        forbidden = (
            "created_unix_ns", "runtime_bindings", "worker_cpus",
            "assignments", "results", "summary", "completion",
        )
        self.assertTrue(all(name not in selection for name in forbidden))
        with self.assertRaises(campaign.CampaignError):
            campaign.make_plan("holdout")
        with self.assertRaises(campaign.CampaignError):
            campaign.make_plan("holdout", "not-frozen")
        holdout = campaign.make_plan("holdout", "pure8_s0m1_d3")
        self.assertEqual(holdout["pre_cpu_job_count"], 594)
        self.assertEqual(holdout["survivor"], "pure8_s0m1_d3")

    def test_campaign_v2_protocol_pin_ignores_mutable_current_alias(self):
        with mock.patch.object(
                campaign.band,
                "BANDTIMING_PROTOCOL",
                campaign.band.BANDTIMING_PROTOCOL_V1):
            selection = campaign.make_plan("selection")
            holdout = campaign.make_plan(
                "holdout", "pure8_s0m1_d3")
        self.assertEqual(
            campaign.REQUIRED_BANDTIMING_PROTOCOL,
            campaign.band.BANDTIMING_PROTOCOL_V2,
        )
        self.assertEqual(
            selection["bandtiming_protocol"],
            campaign.band.BANDTIMING_PROTOCOL_V2,
        )
        self.assertEqual(
            holdout["bandtiming_protocol"],
            campaign.band.BANDTIMING_PROTOCOL_V2,
        )

    def test_result_free_selection_policy_is_bound_into_plan(self):
        policy = campaign.selection_policy()
        plan = campaign.make_plan("selection")
        self.assertEqual(
            plan["frozen_roster"]["selection_policy"], policy)
        self.assertEqual(
            policy["raw_pareto_surface"]["direction"],
            "minimize-every-dimension",
        )
        self.assertEqual(
            policy["selection"]["tie_break"],
            "lexicographically-smallest-candidate-name",
        )
        self.assertIn(
            "both-matched-dispatch-and-matched-WH1",
            policy["selection"]["eligibility"],
        )
        self.assertIn(
            "not-an-eligibility-veto",
            policy["selection"]["weak_seed_rule"],
        )
        self.assertIn(
            "Student-t-95-percent-upper-bound",
            policy["holdout"]["direct_speed_confirmation"]["rule"],
        )
        self.assertEqual(
            policy["raw_evidence"]["timing_completeness"][
                "minimum_eligible_replicates_per_job"],
            {"selection": 384, "holdout": 128},
        )
        self.assertEqual(
            set(policy["production_reliability"][
                "candidate_control_panels"]),
            {
                panel for panel, unused_scope, unused_control
                in campaign.CANDIDATE_CONTROL_PANELS
            },
        )
        self.assertEqual(
            set(policy["production_contract_promotion"]["actions"]),
            {
                "promote-new-contract", "repair-seeds-and-retest",
                "investigate-reliability", "retain-dispatch",
            },
        )

    def test_existing_campaign_directory_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "campaign"
            path.mkdir()
            with self.assertRaisesRegex(
                    campaign.CampaignError, "refusing existing"):
                campaign._create_fresh_directory(path)

    def test_descriptor_seams_at_k2_and_k100(self):
        for block_count in (2, 100):
            dispatch = campaign.band.dispatch_band_descriptor(block_count)
            descriptor = campaign.candidate_descriptor(
                "pure8_s0m1_d3", block_count)
            self.assertEqual(descriptor.staircase, dispatch.staircase - 1)
            self.assertEqual(descriptor.dense_rows, 3)
            self.assertEqual(descriptor.gf256_rows, 8)
            self.assertEqual(descriptor.gf16_rows, 0)
            self.assertEqual(descriptor.period, 244)
            self.assertEqual(descriptor.x_mode, "frozen")

    def test_atomic_no_replace_and_malformed_envelope(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "value.json"
            original_sha = campaign.atomic_json(path, {"value": 1})
            with self.assertRaisesRegex(
                    campaign.CampaignError, "refusing to replace"):
                campaign.atomic_json(path, {"value": 2})
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), original_sha)
            self.assertEqual(json.loads(path.read_text()), {"value": 1})

            jobs = campaign.build_pre_cpu_jobs("selection")
            job = jobs[0]
            assignment = {
                "order": 0,
                "job_id": job["job_id"],
                "cpu": 0,
                "job_file": "jobs/0000.json",
                "output": "results/0000.json.gz",
                "log": "logs/0000.log",
            }
            manifest = {"runtime_bindings": {}, "pre_cpu_jobs": [job]}
            result = directory / assignment["output"]
            result.parent.mkdir()
            campaign.atomic_gzip_json(result, {
                "schema": campaign.RESULT_SCHEMA,
                "job": job,
            })
            with self.assertRaisesRegex(
                    campaign.CampaignError, "envelope is inconsistent"):
                campaign._replay_result(
                    directory, manifest, "0" * 64, assignment, job)

    def test_artifact_replay_rejects_bool_int_aliases(self):
        job = campaign.build_pre_cpu_jobs("selection")[0]
        assignment = {
            "order": job["order"],
            "job_id": job["job_id"],
            "cpu": 0,
            "job_file": "jobs/0000.json",
            "output": "results/0000.json.gz",
            "log": "logs/0000.log",
        }
        bindings = synthetic_runtime_bindings()
        manifest_sha256 = "0" * 64
        manifest = {
            "bandtiming_protocol": campaign.REQUIRED_BANDTIMING_PROTOCOL,
            "runtime_bindings": bindings,
            "pre_cpu_jobs": [job],
            "assignments": [assignment],
        }
        record = campaign._expected_job_record(
            manifest, manifest_sha256, assignment, job)
        envelope = {
            "schema": campaign.RESULT_SCHEMA,
            "manifest_sha256": manifest_sha256,
            "job": job,
            "assignment": assignment,
            "wall_started_unix_ns": 1,
            "wall_finished_unix_ns": 2,
            "runtime_bindings_before": bindings,
            "runtime_bindings_after": bindings,
            "receipt": {},
        }
        mutations = (
            lambda value: value["job"].__setitem__("order", False),
            lambda value: value["assignment"].__setitem__("order", False),
            lambda value: value["runtime_bindings"].get(
                "benchmark").__setitem__("device", True),
        )
        for mutate in mutations:
            with self.subTest(artifact="job", mutation=mutate):
                forged = json.loads(json.dumps(record))
                mutate(forged)
                with mock.patch.object(
                        campaign, "read_canonical_json",
                        return_value=(forged, "1" * 64)):
                    with self.assertRaisesRegex(
                            campaign.CampaignError, "job file changed"):
                        campaign._verify_job_files(
                            Path("/unused"), manifest, manifest_sha256)

        envelope_mutations = (
            lambda value: value["job"].__setitem__("order", False),
            lambda value: value["assignment"].__setitem__("order", False),
            lambda value: value["runtime_bindings_before"][
                "benchmark"].__setitem__("device", True),
            lambda value: value["runtime_bindings_after"][
                "benchmark"].__setitem__("device", True),
        )
        for mutate in envelope_mutations:
            with self.subTest(artifact="result", mutation=mutate):
                forged = json.loads(json.dumps(envelope))
                mutate(forged)
                with (
                        mock.patch.object(
                            campaign, "read_canonical_gzip_json",
                            return_value=(forged, {})),
                        mock.patch.object(
                            campaign.band,
                            "replay_bandtiming_receipt") as replay,
                ):
                    with self.assertRaisesRegex(
                            campaign.CampaignError,
                            "result envelope is inconsistent"):
                        campaign._replay_result(
                            Path("/unused"), manifest, manifest_sha256,
                            assignment, job)
                    replay.assert_not_called()

    def test_result_replay_rejects_legacy_nonpromotable_protocol(self):
        from test_band_timing_tools import (
            band_context,
            bandtiming_stdout,
            parse_kwargs,
        )

        with (
                mock.patch.object(
                    campaign, "SELECTION_REPLICATES", 4),
                tempfile.TemporaryDirectory() as temporary,
        ):
            directory = Path(temporary)
            job = campaign.build_pre_cpu_jobs("selection")[0]
            campaign.validate_job(job)
            candidate = campaign.descriptor_from_job(job)
            context = band_context(cache_state=job["cache_state"])
            fixture = {
                "candidate": candidate,
                "schedule": job["schedule"],
                "systematic_cache": job["systematic_cache"],
                "cache_state": job["cache_state"],
                "block_count": job["K"],
                "block_bytes": job["block_bytes"],
                "replicates": job["replicates"],
                "warmups": job["warmup_replicates"],
                "construction_seed": job["construction_seed_base"],
                "loss_seed": job["loss_seed_base"],
                "inner_reps": job["inner_reps"],
                "max_overhead": job["max_overhead"],
                "required_margin": job["required_margin"],
                "loss": job["loss"],
            }
            legacy_measurement = campaign.band.parse_bandtiming_output(
                bandtiming_stdout(
                    context=context,
                    schema=campaign.band.BANDTIMING_SCHEMA_V1,
                    **fixture,
                ),
                _protocol=campaign.band.BANDTIMING_PROTOCOL_V1,
                **parse_kwargs(context=context, **fixture),
            )
            self.assertFalse(legacy_measurement.valid_for_promotion)
            self.assertEqual(
                campaign.band.replay_bandtiming_receipt(
                    legacy_measurement.as_dict()).protocol,
                campaign.band.BANDTIMING_PROTOCOL_V1,
            )
            assignment = {
                "order": job["order"],
                "job_id": job["job_id"],
                "cpu": 2,
                "job_file": "jobs/0000.json",
                "output": "results/0000.json.gz",
                "log": "logs/0000.log",
            }
            bindings = synthetic_runtime_bindings()
            bindings["thermal"]["device"] = 1
            bindings["thermal"]["inode"] = 2
            manifest_sha256 = "0" * 64
            manifest = {
                "bandtiming_protocol":
                    campaign.REQUIRED_BANDTIMING_PROTOCOL,
                "runtime_bindings": bindings,
                "pre_cpu_jobs": [job],
            }
            envelope = {
                "schema": campaign.RESULT_SCHEMA,
                "manifest_sha256": manifest_sha256,
                "job": job,
                "assignment": assignment,
                "wall_started_unix_ns": 1,
                "wall_finished_unix_ns": 2,
                "runtime_bindings_before": bindings,
                "runtime_bindings_after": bindings,
                "receipt": legacy_measurement.as_dict(),
            }
            with (
                    mock.patch.object(
                        campaign, "read_canonical_gzip_json",
                        return_value=(envelope, {})),
                    mock.patch.object(
                        campaign.band,
                        "BANDTIMING_PROTOCOL",
                        campaign.band.BANDTIMING_PROTOCOL_V1),
            ):
                with self.assertRaisesRegex(
                        campaign.CampaignError, "non-promotable"):
                    campaign._replay_result(
                        directory, manifest, manifest_sha256,
                        assignment, job)

    def test_worker_rejects_probe_and_replay_protocol_downgrades(self):
        bindings = synthetic_runtime_bindings()
        job = campaign.build_pre_cpu_jobs("selection")[0]
        assignment = {
            "order": job["order"],
            "job_id": job["job_id"],
            "cpu": 2,
            "job_file": "jobs/0000.json",
            "output": "results/0000.json.gz",
            "log": "logs/0000.log",
        }
        record = {
            "schema": campaign.JOB_SCHEMA,
            "manifest_sha256": "0" * 64,
            "runtime_bindings": bindings,
            "job": job,
            "assignment": assignment,
        }
        cases = (
            (
                SimpleNamespace(
                    protocol=campaign.band.BANDTIMING_PROTOCOL_V1,
                    valid_for_promotion=False,
                    as_dict=lambda: {}),
                SimpleNamespace(
                    protocol=campaign.band.BANDTIMING_PROTOCOL_V2,
                    valid_for_promotion=True),
                "worker produced a non-promotable receipt",
                False,
            ),
            (
                SimpleNamespace(
                    protocol=campaign.band.BANDTIMING_PROTOCOL_V2,
                    valid_for_promotion=True,
                    as_dict=lambda: {}),
                SimpleNamespace(
                    protocol=campaign.band.BANDTIMING_PROTOCOL_V1,
                    valid_for_promotion=False),
                "worker replay downgraded the receipt protocol",
                True,
            ),
        )
        for probe, replay, diagnostic, replay_called in cases:
            with (
                    self.subTest(diagnostic=diagnostic),
                    tempfile.TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                job_path = directory / assignment["job_file"]
                output_path = directory / assignment["output"]
                with (
                        mock.patch.object(
                            campaign, "_require_normal_priority"),
                        mock.patch.object(
                            campaign, "read_canonical_json",
                            return_value=(record, "unused")),
                        mock.patch.object(
                            campaign.os, "sched_setaffinity"),
                        mock.patch.object(
                            campaign.os, "sched_getaffinity",
                            return_value={assignment["cpu"]}),
                        mock.patch.object(
                            campaign, "runtime_bindings",
                            return_value=bindings),
                        mock.patch.object(
                            campaign.peel_codec,
                            "make_paired_context_config",
                            return_value={}),
                        mock.patch.object(
                            campaign.band, "bandtiming_probe",
                            return_value=probe),
                        mock.patch.object(
                            campaign.band, "replay_bandtiming_receipt",
                            return_value=replay) as replay_mock,
                ):
                    with self.assertRaisesRegex(
                            campaign.CampaignError, diagnostic):
                        campaign.run_worker(job_path, output_path)
                self.assertEqual(replay_mock.called, replay_called)

    def test_malformed_manifest_and_process_group_cleanup(self):
        jobs = campaign.build_pre_cpu_jobs("selection")
        bindings = synthetic_runtime_bindings()
        manifest = campaign._build_manifest(
            "selection", None, jobs, [0], 1, bindings, None)
        campaign._validate_manifest(manifest)
        legacy_protocol = json.loads(json.dumps(manifest))
        legacy_protocol["bandtiming_protocol"] = (
            campaign.band.BANDTIMING_PROTOCOL_V1
        )
        with self.assertRaisesRegex(
                campaign.CampaignError, "runtime contract"):
            campaign._validate_manifest(legacy_protocol)
        aliased_roster_flag = json.loads(json.dumps(manifest))
        aliased_roster_flag["frozen_roster"]["selection_policy"][
            "production_contract_promotion"
        ]["requires_holdout_confirmation"] = 1
        with self.assertRaisesRegex(
                campaign.CampaignError, "runtime contract"):
            campaign._validate_manifest(aliased_roster_flag)
        aliased_population_order = json.loads(json.dumps(manifest))
        aliased_population_order["pre_cpu_jobs"][0]["order"] = False
        aliased_population_order["assignments"][0]["order"] = False
        with self.assertRaisesRegex(
                campaign.CampaignError, "population changed"):
            campaign._validate_manifest(aliased_population_order)
        aliased_assignment_order = json.loads(json.dumps(manifest))
        aliased_assignment_order["assignments"][0]["order"] = False
        with self.assertRaisesRegex(
                campaign.CampaignError, "assignment changed"):
            campaign._validate_manifest(aliased_assignment_order)
        forged = json.loads(json.dumps(manifest))
        forged["assignments"][0]["cpu"] = 1
        with self.assertRaisesRegex(
                campaign.CampaignError, "assignment changed"):
            campaign._validate_manifest(forged)
        malformed_cpu = json.loads(json.dumps(manifest))
        malformed_cpu["worker_cpus"] = [[]]
        with self.assertRaisesRegex(
                campaign.CampaignError, "runtime contract"):
            campaign._validate_manifest(malformed_cpu)
        malformed_job_id = json.loads(json.dumps(manifest))
        malformed_job_id["assignments"][0]["job_id"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "assignments are incomplete"):
            campaign._validate_manifest(malformed_job_id)
        malformed_job = json.loads(json.dumps(jobs[0]))
        malformed_job["candidate"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "outside the frozen roster"):
            campaign.validate_job(malformed_job)
        aliased_job_descriptor = json.loads(json.dumps(jobs[0]))
        aliased_job_descriptor["candidate_descriptor"]["gf16"] = False
        with self.assertRaisesRegex(
                campaign.CampaignError, "differs from its frozen request"):
            campaign.validate_job(aliased_job_descriptor)
        with self.assertRaisesRegex(
                campaign.CampaignError, "explicit survivor"):
            campaign.build_pre_cpu_jobs("holdout", [])
        worker_record = campaign._expected_job_record(
            manifest,
            "0" * 64,
            manifest["assignments"][0],
            jobs[0],
        )
        malformed_worker_candidate = json.loads(
            json.dumps(worker_record))
        malformed_worker_candidate["job"]["candidate"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "outside the frozen roster"):
            campaign._validate_worker_record(
                malformed_worker_candidate)
        malformed_worker_path = json.loads(json.dumps(worker_record))
        malformed_worker_path["assignment"]["output"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "CPU assignment is malformed"):
            campaign._validate_worker_record(malformed_worker_path)
        aliased_worker_order = json.loads(json.dumps(worker_record))
        aliased_worker_order["assignment"]["order"] = False
        with self.assertRaisesRegex(
                campaign.CampaignError, "CPU assignment is malformed"):
            campaign._validate_worker_record(aliased_worker_order)

        survivor = "pure8_s0_d3"
        holdout_jobs = campaign.build_pre_cpu_jobs("holdout", survivor)
        selection_contract = {
            "path": "/tmp/frozen-selection/manifest.json",
            "manifest_sha256": "0" * 64,
            "completion_sha256": "1" * 64,
            "selection_decision_sha256": "2" * 64,
            "selected_survivor": survivor,
            "selected_candidate_production_blockers":
                synthetic_clean_selection_blockers(survivor),
        }
        holdout = campaign._build_manifest(
            "holdout", survivor, holdout_jobs, [0], 1, bindings,
            selection_contract)
        campaign._validate_manifest(holdout)
        forged_holdout = json.loads(json.dumps(holdout))
        forged_holdout["selection_contract"][
            "selected_survivor"] = "pure8_s0m1_d3"
        with self.assertRaisesRegex(
                campaign.CampaignError, "selection contract"):
            campaign._validate_manifest(forged_holdout)
        malformed_survivor = json.loads(json.dumps(holdout))
        malformed_survivor["survivor"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "survivor is malformed"):
            campaign._validate_manifest(malformed_survivor)

        process = mock.Mock()
        process.pid = 12345
        process.poll.return_value = None
        process.wait.return_value = 0
        with mock.patch.object(
                campaign.os, "killpg",
                side_effect=process_group_absent_after_term) as killpg:
            campaign._kill_process_groups(((process, None),))
        self.assertIn(
            mock.call(12345, campaign.signal.SIGTERM),
            killpg.call_args_list,
        )
        process.wait.assert_called_once_with(timeout=5.0)

    def test_selection_contract_rejects_unbound_summary_reread(self):
        jobs = campaign.build_pre_cpu_jobs("selection")
        bindings = synthetic_runtime_bindings()
        survivor = "pure8_s0_d3"
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = campaign._build_manifest(
                "selection", None, jobs, [0], 1, bindings, None)
            manifest_sha256 = campaign.atomic_json(
                directory / "manifest.json", manifest)
            summary_evidence = campaign.atomic_gzip_json(
                directory / "summary.json.gz",
                {
                    "decision": {
                        "schema": campaign.DECISION_SCHEMA,
                        "phase": "selection",
                        "policy_sha256": campaign.canonical_sha256(
                            campaign.selection_policy()),
                        "selected_survivor": survivor,
                    },
                },
            )
            forged_summary_evidence = dict(summary_evidence)
            forged_summary_evidence["compressed_sha256"] = "f" * 64
            completion = {
                "schema": campaign.COMPLETION_SCHEMA,
                "phase": "selection",
                "manifest_sha256": manifest_sha256,
                "jobs": 2376,
                "summary": {
                    "path": "summary.json.gz",
                    **forged_summary_evidence,
                },
            }
            completion_sha256 = campaign.atomic_json(
                directory / "completion.json", completion)
            campaign._write_completion_checksum(
                directory, completion_sha256)
            verified = {
                "manifest_sha256": manifest_sha256,
                "completion_sha256": completion_sha256,
                "jobs": 2376,
            }
            with (
                mock.patch.object(
                    campaign, "verify_campaign", return_value=verified),
                self.assertRaisesRegex(
                    campaign.CampaignError,
                    "not bound by completion summary hash"),
            ):
                campaign._load_selection_contract(
                    directory / "manifest.json", survivor)

    def test_selection_contract_rejects_non_object_artifacts(self):
        jobs = campaign.build_pre_cpu_jobs("selection")
        bindings = synthetic_runtime_bindings()
        survivor = "pure8_s0_d3"
        for artifact in ("completion", "summary"):
            with self.subTest(artifact=artifact), \
                    tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                manifest = campaign._build_manifest(
                    "selection", None, jobs, [0], 1, bindings, None)
                manifest_sha256 = campaign.atomic_json(
                    directory / "manifest.json", manifest)
                if artifact == "completion":
                    completion_sha256 = campaign.atomic_json(
                        directory / "completion.json", [])
                    campaign._write_completion_checksum(
                        directory, completion_sha256)
                    with self.assertRaisesRegex(
                            campaign.CampaignError,
                            "selection completion is not an object"):
                        campaign._load_selection_contract(
                            directory / "manifest.json", survivor)
                    continue

                summary_evidence = campaign.atomic_gzip_json(
                    directory / "summary.json.gz", None)
                completion = {
                    "schema": campaign.COMPLETION_SCHEMA,
                    "phase": "selection",
                    "manifest_sha256": manifest_sha256,
                    "jobs": 2376,
                    "summary": {
                        "path": "summary.json.gz",
                        **summary_evidence,
                    },
                }
                completion_sha256 = campaign.atomic_json(
                    directory / "completion.json", completion)
                campaign._write_completion_checksum(
                    directory, completion_sha256)
                verified = {
                    "manifest_sha256": manifest_sha256,
                    "completion_sha256": completion_sha256,
                    "jobs": 2376,
                }
                with (
                    mock.patch.object(
                        campaign, "verify_campaign",
                        return_value=verified),
                    self.assertRaisesRegex(
                        campaign.CampaignError,
                        "selection summary is not an object"),
                ):
                    campaign._load_selection_contract(
                        directory / "manifest.json", survivor)

    def test_stored_summary_rejects_canonical_non_objects(self):
        for value in ([], None):
            with self.subTest(value=value), \
                    tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                evidence = campaign.atomic_gzip_json(
                    directory / "summary.json.gz", value)
                completion = {
                    "phase": "selection",
                    "pre_cpu_job_list_sha256": "0" * 64,
                    "summary": {
                        "path": "summary.json.gz",
                        **evidence,
                    },
                }
                with self.assertRaisesRegex(
                        campaign.CampaignError,
                        "stored summary does not match completion"):
                    campaign._verify_stored_summary_artifacts(
                        directory, completion)

    def test_capped_gzip_json_preserves_exact_boundary_semantics_and_hashes(
            self):
        value = {
            "schema": campaign.SUMMARY_SCHEMA,
            "phase": "selection",
            "historical_fixture": [None, False, 7, "ascii"],
        }
        payload = campaign.canonical_json_bytes(value)
        compressed = gzip.compress(
            payload, compresslevel=6, mtime=0)
        expected_evidence = {
            "compressed_sha256":
                hashlib.sha256(compressed).hexdigest(),
            "uncompressed_sha256":
                hashlib.sha256(payload).hexdigest(),
            "compressed_bytes": len(compressed),
            "uncompressed_bytes": len(payload),
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "summary.json.gz"
            path.write_bytes(compressed)
            observed, evidence = campaign.read_canonical_gzip_json(
                path,
                compressed_limit=len(compressed),
                uncompressed_limit=len(payload),
            )
            self.assertTrue(
                campaign.peel_codec._same_typed_json(observed, value))
            self.assertEqual(evidence, expected_evidence)

            for dimension, limits, expected in (
                    (
                        "compressed",
                        (len(compressed) - 1, len(payload)),
                        "exceeds compressed byte limit",
                    ),
                    (
                        "uncompressed",
                        (len(compressed), len(payload) - 1),
                        "exceeds uncompressed byte limit",
                    ),
            ):
                with (
                        self.subTest(dimension=dimension),
                        mock.patch.object(
                            campaign, "_strict_json_payload") as parse,
                        self.assertRaisesRegex(
                            campaign.CampaignError, expected),
                ):
                    campaign.read_canonical_gzip_json(
                        path,
                        compressed_limit=limits[0],
                        uncompressed_limit=limits[1],
                    )
                parse.assert_not_called()

    def test_capped_gzip_json_rejects_truncation_and_trailing_data(
            self):
        payload = campaign.canonical_json_bytes({"valid": True})
        compressed = gzip.compress(
            payload, compresslevel=6, mtime=0)
        cases = (
            ("truncated", compressed[:-1], "is truncated"),
            ("zero_padding", compressed + b"\0", "has trailing data"),
            (
                "concatenated_member",
                compressed + gzip.compress(b"", mtime=0),
                "has trailing data",
            ),
        )
        for name, malformed, expected in cases:
            with (
                    self.subTest(name=name),
                    tempfile.TemporaryDirectory() as temporary,
            ):
                path = Path(temporary) / "summary.json.gz"
                path.write_bytes(malformed)
                with (
                        mock.patch.object(
                            campaign, "_strict_json_payload") as parse,
                        self.assertRaisesRegex(
                            campaign.CampaignError, expected),
                ):
                    campaign.read_canonical_gzip_json(
                        path,
                        compressed_limit=len(malformed),
                        uncompressed_limit=len(payload),
                    )
                parse.assert_not_called()

    def test_capped_gzip_json_rejects_compressed_and_inflated_bombs(
            self):
        uncompressed_limit = 1024
        cases = (
            (
                "compressed",
                b"x" * 1025,
                1024,
                uncompressed_limit,
                "exceeds compressed byte limit",
            ),
            (
                "inflated",
                gzip.compress(
                    b"x" * (uncompressed_limit + 1), mtime=0),
                1024,
                uncompressed_limit,
                "exceeds uncompressed byte limit",
            ),
        )
        for name, malformed, compressed_limit, raw_limit, expected in cases:
            with (
                    self.subTest(name=name),
                    tempfile.TemporaryDirectory() as temporary,
            ):
                path = Path(temporary) / "summary.json.gz"
                path.write_bytes(malformed)
                with (
                        mock.patch.object(
                            campaign, "_strict_json_payload") as parse,
                        self.assertRaisesRegex(
                            campaign.CampaignError, expected),
                ):
                    campaign.read_canonical_gzip_json(
                        path,
                        compressed_limit=compressed_limit,
                        uncompressed_limit=raw_limit,
                    )
                parse.assert_not_called()

    def test_stored_evidence_rejects_numeric_type_aliases(self):
        summary_value = {
            "schema": campaign.SUMMARY_SCHEMA,
            "phase": "selection",
            "pre_cpu_job_list_sha256": "0" * 64,
        }
        summary_evidence = {
            "compressed_sha256": "1" * 64,
            "uncompressed_sha256": "2" * 64,
            "compressed_bytes": 10,
            "uncompressed_bytes": 20,
        }
        ledger = {
            "path": "cell-ledger.jsonl.gz",
            "rows": 1,
            "compressed_sha256": "3" * 64,
            "uncompressed_sha256": "4" * 64,
            "compressed_bytes": 30,
            "uncompressed_bytes": 40,
        }
        completion = {
            "phase": "selection",
            "pre_cpu_job_list_sha256": "0" * 64,
            "summary": {
                "path": "summary.json.gz",
                **summary_evidence,
            },
            "cell_ledger": ledger,
        }
        with (
                mock.patch.object(
                    campaign, "read_canonical_gzip_json",
                    return_value=(summary_value, summary_evidence)),
                mock.patch.object(
                    campaign, "_stream_file_sha256",
                    return_value=("3" * 64, 30)),
        ):
            forged_summary = json.loads(json.dumps(completion))
            forged_summary["summary"]["compressed_bytes"] = 10.0
            with self.assertRaisesRegex(
                    campaign.CampaignError,
                    "stored summary does not match completion"):
                campaign._verify_stored_summary_artifacts(
                    Path("/unused"), forged_summary)

            forged_ledger = json.loads(json.dumps(completion))
            forged_ledger["cell_ledger"]["compressed_bytes"] = 30.0
            with self.assertRaisesRegex(
                    campaign.CampaignError,
                    "stored cell ledger does not match completion"):
                campaign._verify_stored_summary_artifacts(
                    Path("/unused"), forged_ledger)

    def test_runtime_binding_nested_schema_fails_closed(self):
        jobs = campaign.build_pre_cpu_jobs("selection")
        bindings = synthetic_runtime_bindings()
        manifest = campaign._build_manifest(
            "selection", None, jobs, [0], 1, bindings, None)
        for name in bindings:
            for malformed in ({}, [], None):
                with self.subTest(name=name, malformed=malformed):
                    forged = json.loads(json.dumps(manifest))
                    forged["runtime_bindings"][name] = malformed
                    with self.assertRaisesRegex(
                            campaign.CampaignError,
                            "runtime .*binding schema"):
                        campaign._validate_manifest(forged)
                    with self.assertRaisesRegex(
                            campaign.CampaignError,
                            "runtime .*binding schema"):
                        campaign.check_runtime_bindings(
                            forged["runtime_bindings"], full_hash=True)
                    record = {
                        "schema": campaign.JOB_SCHEMA,
                        "manifest_sha256": "0" * 64,
                        "runtime_bindings": forged["runtime_bindings"],
                        "job": {},
                        "assignment": {},
                    }
                    with self.assertRaisesRegex(
                            campaign.CampaignError,
                            "runtime .*binding schema"):
                        campaign._validate_worker_record(record)

    def test_post_popen_registration_failure_reaps_process_group(self):
        process = mock.Mock()
        process.pid = 24680
        process.poll.return_value = None
        process.wait.return_value = 0
        assignment = {
            "order": 0,
            "job_file": "jobs/0000.json",
            "output": "results/0000.json.gz",
            "log": "logs/0000.log",
        }
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "logs").mkdir()
            with (
                mock.patch.object(
                    campaign.subprocess, "Popen", return_value=process),
                mock.patch.object(
                    campaign, "_register_process",
                    side_effect=MemoryError("injected registration failure")),
                mock.patch.object(
                    campaign.os, "killpg",
                    side_effect=process_group_absent_after_term) as killpg,
            ):
                with self.assertRaisesRegex(
                        MemoryError, "registration failure"):
                    campaign._run_wave(
                        directory, {"runtime_bindings": {}}, (assignment,))
        self.assertIn(
            mock.call(24680, signal.SIGTERM),
            killpg.call_args_list,
        )
        process.wait.assert_called_once_with(timeout=5.0)

    def test_exited_leader_does_not_orphan_term_ignoring_descendant(self):
        leader_source = r"""
import os
import pathlib
import subprocess
import sys
import time

ready = sys.argv[1]
child_source = r'''
import os
import pathlib
import signal
import sys
import time

signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text(str(os.getpid()), encoding="ascii")
time.sleep(60)
'''
descendant = subprocess.Popen(
    [sys.executable, "-c", child_source, ready],
    stdin=subprocess.DEVNULL,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
deadline = time.monotonic() + 10.0
while not pathlib.Path(ready).exists():
    if time.monotonic() >= deadline:
        raise SystemExit("descendant did not become ready")
    time.sleep(0.01)
os._exit(0)
"""
        with tempfile.TemporaryDirectory() as temporary:
            ready = Path(temporary) / "descendant.pid"
            leader = subprocess.Popen(
                [sys.executable, "-c", leader_source, str(ready)],
                start_new_session=True,
            )
            child_pid = None
            try:
                self.assertEqual(leader.wait(timeout=10), 0)
                child_pid = int(ready.read_text(encoding="ascii"))
                os.kill(child_pid, 0)
                os.killpg(leader.pid, 0)
                campaign._kill_process_groups_with_grace(
                    ((leader, None),), 0.1)
                with self.assertRaises(ProcessLookupError):
                    os.kill(child_pid, 0)
                with self.assertRaises(ProcessLookupError):
                    os.killpg(leader.pid, 0)
            finally:
                try:
                    os.killpg(leader.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                if child_pid is not None:
                    try:
                        os.kill(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                if leader.poll() is None:
                    leader.kill()
                    leader.wait()

    def test_first_signal_during_existing_failure_cannot_abort_cleanup(self):
        worker_source = r"""
import os
import pathlib
import signal
import sys
import time

signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text(str(os.getpid()), encoding="ascii")
time.sleep(60)
"""
        with tempfile.TemporaryDirectory() as temporary:
            ready = Path(temporary) / "worker.pid"
            worker = subprocess.Popen(
                [sys.executable, "-c", worker_source, str(ready)],
                start_new_session=True,
            )
            sender = None
            try:
                deadline = time.monotonic() + 10.0
                while not ready.exists():
                    if worker.poll() is not None:
                        self.fail("TERM-ignoring worker exited before cleanup")
                    if time.monotonic() >= deadline:
                        self.fail("TERM-ignoring worker did not become ready")
                    time.sleep(0.01)

                def send_term_during_grace():
                    time.sleep(0.05)
                    os.kill(os.getpid(), signal.SIGTERM)

                with self.assertRaisesRegex(
                        RuntimeError, "preexisting failure"):
                    with campaign._parent_termination_handlers():
                        try:
                            raise RuntimeError("preexisting failure")
                        except RuntimeError:
                            sender = threading.Thread(
                                target=send_term_during_grace)
                            sender.start()
                            try:
                                campaign._kill_process_groups_with_grace(
                                    ((worker, None),), 0.15)
                            finally:
                                sender.join(timeout=2.0)
                            self.assertFalse(sender.is_alive())
                            self.assertEqual(
                                campaign._PARENT_SIGNAL_STATE["pending"],
                                signal.SIGTERM,
                            )
                            with self.assertRaises(ProcessLookupError):
                                os.killpg(worker.pid, 0)
                            raise
                with self.assertRaises(ProcessLookupError):
                    os.kill(worker.pid, 0)
            finally:
                if sender is not None:
                    sender.join(timeout=2.0)
                try:
                    os.killpg(worker.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                if worker.poll() is None:
                    worker.kill()
                    worker.wait()

    def test_signal_at_exception_cleanup_transition_cannot_orphan_group(self):
        worker_source = r"""
import pathlib
import signal
import sys
import time

signal.signal(signal.SIGTERM, signal.SIG_IGN)
pathlib.Path(sys.argv[1]).write_text(str(os.getpid()), encoding="ascii")
time.sleep(60)
"""
        # Keep the child source self-contained so the tracing regression
        # exercises only the parent signal path.
        worker_source = "import os\n" + worker_source
        with tempfile.TemporaryDirectory() as temporary:
            ready = Path(temporary) / "worker.pid"
            worker = subprocess.Popen(
                [sys.executable, "-c", worker_source, str(ready)],
                start_new_session=True,
            )
            trace_fired = False
            try:
                deadline = time.monotonic() + 10.0
                while not ready.exists():
                    if worker.poll() is not None:
                        self.fail("TERM-ignoring worker exited before cleanup")
                    if time.monotonic() >= deadline:
                        self.fail("TERM-ignoring worker did not become ready")
                    time.sleep(0.01)

                cleanup_code = \
                    campaign._kill_process_groups_with_grace.__code__

                def inject_at_cleanup_call(frame, event, unused_argument):
                    nonlocal trace_fired
                    if event == "call" and frame.f_code is cleanup_code:
                        trace_fired = True
                        sys.settrace(None)
                        os.kill(os.getpid(), signal.SIGTERM)
                    return inject_at_cleanup_call

                with self.assertRaisesRegex(
                        RuntimeError, "preexisting failure"):
                    with campaign._parent_termination_handlers():
                        try:
                            raise RuntimeError("preexisting failure")
                        except RuntimeError:
                            sys.settrace(inject_at_cleanup_call)
                            try:
                                campaign._kill_process_groups_with_grace(
                                    ((worker, None),), 0.1)
                            finally:
                                sys.settrace(None)
                            self.assertTrue(trace_fired)
                            self.assertEqual(
                                campaign._PARENT_SIGNAL_STATE["pending"],
                                signal.SIGTERM,
                            )
                            with self.assertRaises(ProcessLookupError):
                                os.killpg(worker.pid, 0)
                            raise
                with self.assertRaises(ProcessLookupError):
                    os.kill(worker.pid, 0)
            finally:
                sys.settrace(None)
                try:
                    os.killpg(worker.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                if worker.poll() is None:
                    worker.kill()
                    worker.wait()

    def test_run_restores_parent_signal_handlers(self):
        term_before = signal.getsignal(signal.SIGTERM)
        hup_before = signal.getsignal(signal.SIGHUP)
        with mock.patch.object(
                campaign, "_run_campaign_impl", return_value={"ok": True}):
            self.assertEqual(
                campaign.run_campaign(
                    "selection", "unused", "unused", "unused"),
                {"ok": True},
            )
        self.assertIs(signal.getsignal(signal.SIGTERM), term_before)
        self.assertIs(signal.getsignal(signal.SIGHUP), hup_before)

    def test_run_raises_signal_recorded_as_impl_returns(self):
        term_before = signal.getsignal(signal.SIGTERM)
        hup_before = signal.getsignal(signal.SIGHUP)

        def finish_after_term(*unused_args, **unused_kwargs):
            os.kill(os.getpid(), signal.SIGTERM)
            return {"must_not_escape": True}

        with (
            mock.patch.object(
                campaign, "_run_campaign_impl",
                side_effect=finish_after_term),
            self.assertRaisesRegex(
                campaign.CampaignError, "received SIGTERM"),
        ):
            campaign.run_campaign(
                "selection", "unused", "unused", "unused")
        self.assertIs(signal.getsignal(signal.SIGTERM), term_before)
        self.assertIs(signal.getsignal(signal.SIGHUP), hup_before)

    def test_sigterm_reaps_real_worker_process_group(self):
        probe_source = r"""
import os
import pathlib
import subprocess
import sys

sys.path.insert(0, sys.argv[1])
import wh2_za5v_campaign as campaign

directory = pathlib.Path(sys.argv[2])
ready = pathlib.Path(sys.argv[3])
(directory / "logs").mkdir(parents=True)
real_popen = subprocess.Popen
child = {}

def launch_sleeper(unused_command, **kwargs):
    process = real_popen(
        [sys.executable, "-c", "import time; time.sleep(60)"],
        **kwargs)
    child["process"] = process
    ready_temporary = ready.with_name(ready.name + ".tmp")
    ready_temporary.write_text(
        f"{process.pid} {process.pid}\n", encoding="ascii")
    os.replace(ready_temporary, ready)
    # Deliver TERM before Popen returns to _run_wave.  The record-only
    # handler must defer raising until the child group is registered.
    import signal
    os.kill(os.getpid(), signal.SIGTERM)
    return process

def bindings_ready(unused_bindings, full_hash):
    pass

campaign.subprocess.Popen = launch_sleeper
campaign.check_runtime_bindings = bindings_ready
assignment = {
    "order": 0,
    "job_file": "jobs/0000.json",
    "output": "results/0000.json.gz",
    "log": "logs/0000.log",
}
manifest = {"runtime_bindings": {}}

def run_wave(*unused_args, **unused_kwargs):
    return campaign._run_wave(directory, manifest, (assignment,))

campaign._run_campaign_impl = run_wave
try:
    campaign.run_campaign("selection", directory, "unused", "unused")
except campaign.CampaignError as error:
    if "received SIGTERM" not in str(error):
        raise
    raise SystemExit(42)
raise SystemExit("campaign unexpectedly completed")
"""
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            ready = temporary / "ready"
            directory = temporary / "campaign"
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    probe_source,
                    str(BENCH),
                    str(directory),
                    str(ready),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            child_pid = None
            try:
                deadline = time.monotonic() + 10.0
                while time.monotonic() < deadline:
                    if ready.exists():
                        child_pid = int(
                            ready.read_text(encoding="ascii").split()[0])
                        break
                    if process.poll() is not None:
                        break
                    time.sleep(0.01)
                if child_pid is None:
                    stdout, stderr = process.communicate(timeout=2)
                    self.fail(
                        "signal cleanup probe never launched a worker: "
                        f"status={process.returncode}, "
                        f"stdout={stdout!r}, stderr={stderr!r}")
                stdout, stderr = process.communicate(timeout=10)
                self.assertEqual(
                    process.returncode, 42,
                    f"stdout={stdout!r}, stderr={stderr!r}")
                with self.assertRaises(ProcessLookupError):
                    os.kill(child_pid, 0)
                with self.assertRaises(ProcessLookupError):
                    os.killpg(child_pid, 0)
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
                for pipe in (process.stdout, process.stderr):
                    if pipe is not None and not pipe.closed:
                        pipe.close()
                if child_pid is not None:
                    try:
                        os.killpg(child_pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass


def synthetic_replicate(job, replicate_index):
    construction = campaign.construction_seed(
        job["construction_seed_base"], replicate_index)
    loss = campaign.loss_seed(job["loss_seed_base"], replicate_index)
    candidate = "success"
    dispatch = "success"
    candidate_overhead = 1
    dispatch_overhead = 2
    if replicate_index == 0:
        candidate = "need_more"
        candidate_overhead = -1
        dispatch_overhead = 0
    elif replicate_index == 1:
        candidate_overhead = 2
        dispatch = "need_more"
        dispatch_overhead = -1
    elif replicate_index == 2:
        candidate = dispatch = "weak"
        candidate_overhead = dispatch_overhead = -1
    elif replicate_index == 3:
        candidate = "error"
        candidate_overhead = -1
        dispatch_overhead = 0

    candidate_constructor = (
        "weak" if candidate == "weak" else "success")
    dispatch_constructor = (
        "weak" if dispatch == "weak" else "success")
    censored = (
        ["decoder_candidate_dispatch"]
        if candidate != "success" or dispatch != "success" else []
    )
    wh1 = "need_more" if replicate_index in (2, 4) else "success"
    wh1_overhead = -1 if wh1 != "success" else 1
    return {
        "replicate": replicate_index,
        "construction_seed": construction,
        "loss_seed": loss,
        "trace_sha256": hashlib.sha256(
            f"{job['K']}:{replicate_index}".encode("ascii")
        ).hexdigest(),
        "candidate_construct_result":
            4 if candidate_constructor == "weak" else 0,
        "candidate_construct_class": candidate_constructor,
        "dispatch_construct_result":
            4 if dispatch_constructor == "weak" else 0,
        "dispatch_construct_class": dispatch_constructor,
        "wh1_construct_result": 0,
        "wh1_construct_class": "success",
        "encoder_candidate_result_class": candidate_constructor,
        "encoder_candidate_overhead":
            -1 if candidate_constructor == "weak" else 0,
        "encoder_dispatch_result_class": dispatch_constructor,
        "encoder_dispatch_overhead":
            -1 if dispatch_constructor == "weak" else 0,
        "encoder_wh1_result_class": "success",
        "encoder_wh1_overhead": 0,
        "decoder_candidate_result_class": candidate,
        "decoder_candidate_overhead": candidate_overhead,
        "decoder_dispatch_result_class": dispatch,
        "decoder_dispatch_overhead": dispatch_overhead,
        "decoder_wh1_result_class": wh1,
        "decoder_wh1_overhead": wh1_overhead,
        "direct_candidate_result_class": candidate_constructor,
        "direct_candidate_overhead":
            -1 if candidate_constructor == "weak" else 0,
        "direct_dispatch_result_class": dispatch_constructor,
        "direct_dispatch_overhead":
            -1 if dispatch_constructor == "weak" else 0,
        "censored_panels": censored,
        "fault_contaminated_panels": [],
    }


def synthetic_job_evidence(job, cpu, cpu_temp, dimm_temp):
    required_lower, unused_upper = \
        campaign.peel_codec._paired_practical_log_bounds(
            campaign.REQUIRED_MARGIN)
    required_floor = -required_lower
    arms = [
        {
            "scope": scope,
            "arm": arm,
            "eligible": 1,
            "aa_eligible_replicates": job["replicates"],
            "aa_log_cost_mean": 0.0,
            "aa_log_cost_ci_low": -0.001,
            "aa_log_cost_ci_high": 0.001,
            "aa_floor_log": 0.001,
        }
        for scope, arm in campaign.TIMING_ARM_KEYS
    ]
    contrasts = [
        {
            "name": name,
            "eligible_replicates": job["replicates"],
            "log_cost_mean": 0.0,
            "log_cost_ci_low": -0.01,
            "log_cost_ci_high": 0.01,
            "left_aa_floor_log": 0.001,
            "right_aa_floor_log": 0.001,
            "effective_floor_log": required_floor,
            "required_margin": campaign.REQUIRED_MARGIN,
            "recovery_regressions": 0,
            "recovery_improvements": 0,
            "both_failures": 0,
            "timing_ci_available": True,
            "left_faster": False,
        }
        for name in campaign.band.BANDTIMING_CROSS_PANELS
    ]
    digest = hashlib.sha256(job["job_id"].encode("ascii")).hexdigest()

    def replicate(censored, contaminated):
        return {
            "candidate_construct_result": 0,
            "candidate_construct_class": "success",
            "encoder_candidate_result_class": "success",
            "encoder_dispatch_result_class": "success",
            "encoder_wh1_result_class": "success",
            "decoder_candidate_result_class": "success",
            "decoder_dispatch_result_class": "success",
            "decoder_wh1_result_class": "success",
            "direct_candidate_result_class": "success",
            "direct_dispatch_result_class": "success",
            "censored_panels": censored,
            "fault_contaminated_panels": contaminated,
        }

    return {
        "stream_sha256": digest,
        "arms": arms,
        "contrasts": contrasts,
        "replicates": [
            replicate(["decoder_candidate_dispatch"], []),
            replicate([], ["encoder_candidate_dispatch"]),
        ],
        "context": {
            "bound": {
                "cpu_affinity": [cpu],
                "cache_state": campaign.CACHE_STATE,
                "thermal_device": 7,
                "thermal_inode": 11,
                "cpu_model": "test-cpu",
            },
            "thermal": {
                "rows": 3,
                "valid_rows": 3,
                "missing_read_rows": 0,
                "invalid_rows": 0,
                "cpu_tctl_max_c": cpu_temp,
                "dimm_max_c": dimm_temp,
                "edac_ce_max": 0,
                "edac_ue_max": 0,
            },
        },
        "evidence": {
            "context_sha256": digest,
            "final_context_sha256": digest,
        },
    }


def synthetic_decision_inputs(specs, phase="selection"):
    raw_replicates = (
        campaign.SELECTION_REPLICATES
        if phase == "selection"
        else campaign.HOLDOUT_REPLICATES
    )
    required_lower, unused_upper = \
        campaign.peel_codec._paired_practical_log_bounds(
            campaign.REQUIRED_MARGIN)
    required_floor = -required_lower
    groups = []
    timing_jobs = []
    weak = {}
    for candidate, spec in specs.items():
        weak[candidate] = {
            "count": spec.get("weak", 0),
            "cells": [],
        }
        for block_count in campaign.K_VALUES:
            for schedule in campaign.SCHEDULES:
                def arm(level):
                    return {
                        "final_unrecovered": level,
                        "unrecovered_by_overhead_h_0_to_64":
                            [level] * (campaign.MAX_OVERHEAD + 1),
                    }

                groups.append({
                    "candidate_name": candidate,
                    "K": block_count,
                    "schedule": schedule,
                    "candidate": arm(spec["recovery"]),
                    "dispatch": arm(2),
                    "wh1": arm(3),
                })
                arms = [
                    {
                        "scope": scope,
                        "arm": arm_name,
                        "aa_eligible_replicates": spec.get(
                            "aa_eligible_replicates", raw_replicates),
                        "aa_log_cost_mean": 0.0,
                        "aa_log_cost_ci_low": spec.get(
                            "aa_ci_low", -0.001),
                        "aa_log_cost_ci_high": spec.get(
                            "aa_ci_high", 0.001),
                        "aa_floor_log": spec.get("aa_floor", 0.001),
                    }
                    for scope, arm_name in campaign.TIMING_ARM_KEYS
                ]

                def contrast(panel, mean):
                    available = spec.get("timing_complete", True)
                    low = spec.get(
                        "ci_low", mean - spec.get("ci_half", 0.001))
                    high = spec.get(
                        "ci_high", mean + spec.get("ci_half", 0.001))
                    floor = spec.get(
                        "effective_floor", required_floor)
                    regressions = spec.get(
                        f"{panel}_regressions", 0)
                    weak_regressions = spec.get(
                        f"{panel}_weak_regressions", 0)
                    nonrepairable = regressions - weak_regressions
                    return {
                        "name": panel,
                        "eligible_replicates": spec.get(
                            "eligible_replicates", raw_replicates),
                        "log_cost_mean": mean,
                        "log_cost_ci_low": low if available else None,
                        "log_cost_ci_high": high if available else None,
                        "left_aa_floor_log": spec.get("aa_floor", 0.001),
                        "right_aa_floor_log": spec.get("aa_floor", 0.001),
                        "effective_floor_log": floor,
                        "required_margin": campaign.REQUIRED_MARGIN,
                        "recovery_regressions": regressions,
                        "recovery_improvements": spec.get(
                            f"{panel}_improvements", 0),
                        "both_failures": spec.get(
                            f"{panel}_both_failures", 0),
                        "candidate_constructor_weak_regressions":
                            weak_regressions,
                        "candidate_nonrepairable_regressions":
                            nonrepairable,
                        "timing_ci_available": available,
                        "left_faster": (
                            available
                            and regressions == 0
                            and high < -floor
                        ),
                    }

                contrasts = []
                for panel in campaign.band.BANDTIMING_CROSS_PANELS:
                    if panel.startswith("direct_candidate"):
                        mean = spec["direct"]
                    elif panel.startswith("encoder_candidate"):
                        mean = spec["encoder"]
                    elif panel.startswith("decoder_candidate"):
                        mean = spec["decoder"]
                    else:
                        mean = 0.0
                    contrasts.append(contrast(panel, mean))
                timing_jobs.append({
                    "candidate": candidate,
                    "K": block_count,
                    "schedule": schedule,
                    "arms": arms,
                    "contrasts": contrasts,
                    "candidate_error_outcomes": {
                        "constructor": spec.get(
                            "constructor_errors", 0),
                        "encoder": spec.get("encoder_errors", 0),
                        "decoder": spec.get("decoder_errors", 0),
                        "direct": spec.get("direct_errors", 0),
                        "total": (
                            spec.get("constructor_errors", 0)
                            + spec.get("encoder_errors", 0)
                            + spec.get("decoder_errors", 0)
                            + spec.get("direct_errors", 0)
                        ),
                    },
                })
    return (
        {
            "groups": groups,
            "independent_construction_weaknesses": {
                "candidate": weak,
            },
        },
        {"timing_jobs": timing_jobs},
    )


def synthetic_contrast(job, name):
    return next(
        item for item in job["contrasts"]
        if item["name"] == name)


def synthetic_arm(job, scope, arm):
    return next(
        item for item in job["arms"]
        if item["scope"] == scope and item["arm"] == arm)


class ParallelReplayTests(unittest.TestCase):
    def test_serial_parallel_and_import_only_wrapper_are_byte_identical(self):
        with (
                mock.patch.object(campaign, "K_VALUES", (2,)),
                mock.patch.object(campaign, "SELECTION_REPLICATES", 3),
                tempfile.TemporaryDirectory() as temporary,
        ):
            root = Path(temporary)
            evidence = root / "evidence"
            manifest = parallel_replay_fixture(evidence)
            serial_directory = root / "serial"
            parallel_directory = root / "parallel"
            wrapper_directory = root / "wrapper"
            serial = campaign.build_summaries(
                evidence,
                manifest,
                serial_directory,
                replay_workers=1,
            )
            parallel = campaign.build_summaries(
                evidence,
                manifest,
                parallel_directory,
            )
            with mock.patch.object(
                    parallel_verify,
                    "SUPPORTED_RUNNER_SHA256",
                    frozenset({
                        manifest["runtime_bindings"]["runner"]["sha256"],
                    })), mock.patch.object(
                        parallel_verify,
                        "_load_frozen_campaign",
                        return_value=campaign,
                    ), mock.patch.object(
                        parallel_verify,
                        "_assert_source_attestation",
                    ):
                wrapped = parallel_verify.parallel_build_summaries(
                    campaign,
                    Path(campaign.__file__),
                    evidence,
                    manifest,
                    wrapper_directory,
                )
            self.assertTrue(
                campaign.peel_codec._same_typed_json(serial, parallel))
            self.assertTrue(
                campaign.peel_codec._same_typed_json(serial, wrapped))
            for name in ("cell-ledger.jsonl.gz", "summary.json.gz"):
                reference = (serial_directory / name).read_bytes()
                self.assertEqual(
                    (parallel_directory / name).read_bytes(), reference)
                self.assertEqual(
                    (wrapper_directory / name).read_bytes(), reference)

    def test_parallel_replay_preserves_serial_fault_rejection(self):
        for fault in ("envelope", "derived"):
            with (
                    self.subTest(fault=fault),
                    mock.patch.object(campaign, "K_VALUES", (2,)),
                    mock.patch.object(
                        campaign, "SELECTION_REPLICATES", 3),
                    tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                evidence = root / "evidence"
                manifest = parallel_replay_fixture(
                    evidence, fault=fault)
                with self.assertRaises(Exception) as serial_error:
                    campaign.build_summaries(
                        evidence,
                        manifest,
                        root / "serial",
                        replay_workers=1,
                    )
                with self.assertRaises(Exception) as parallel_error:
                    campaign.build_summaries(
                        evidence,
                        manifest,
                        root / "parallel",
                        replay_workers=2,
                    )
                with (
                        mock.patch.object(
                            parallel_verify,
                            "SUPPORTED_RUNNER_SHA256",
                            frozenset({
                                manifest["runtime_bindings"][
                                    "runner"]["sha256"],
                            })),
                        mock.patch.object(
                            parallel_verify,
                            "_load_frozen_campaign",
                            return_value=campaign,
                        ),
                        mock.patch.object(
                            parallel_verify,
                            "_assert_source_attestation",
                        ),
                        self.assertRaises(Exception) as wrapper_error,
                ):
                    parallel_verify.parallel_build_summaries(
                        campaign,
                        Path(campaign.__file__),
                        evidence,
                        manifest,
                        root / "wrapper",
                        replay_workers=2,
                    )
                for observed in (
                        parallel_error.exception,
                        wrapper_error.exception):
                    self.assertIs(
                        type(observed), type(serial_error.exception))
                    self.assertEqual(
                        str(observed), str(serial_error.exception))
                for output in (
                        root / "serial",
                        root / "parallel",
                        root / "wrapper"):
                    self.assertFalse(
                        (output / "cell-ledger.jsonl.gz").exists())
                    self.assertFalse(
                        (output / "summary.json.gz").exists())

    def test_import_only_verifier_patches_and_restores_only_builder(self):
        with (
                mock.patch.object(campaign, "K_VALUES", (2,)),
                mock.patch.object(campaign, "SELECTION_REPLICATES", 3),
                tempfile.TemporaryDirectory() as temporary,
        ):
            root = Path(temporary)
            evidence = root / "evidence"
            manifest = parallel_replay_fixture(evidence)
            serial = campaign.build_summaries(
                evidence,
                manifest,
                root / "serial",
                replay_workers=1,
            )
            manifest_path = evidence / "manifest.json"
            original_read = campaign.read_canonical_json
            original_build = campaign.build_summaries

            def read_fixture(path):
                if Path(path) == manifest_path:
                    return manifest, campaign.canonical_sha256(manifest)
                return original_read(path)

            def invoke_patched_builder(unused_directory):
                self.assertIsNot(
                    campaign.build_summaries, original_build)
                return campaign.build_summaries(
                    evidence, manifest, root / "wrapped")

            with (
                    mock.patch.object(
                        parallel_verify,
                        "_load_frozen_campaign",
                        return_value=campaign),
                    mock.patch.object(
                        parallel_verify,
                        "_assert_source_attestation"),
                    mock.patch.object(
                        parallel_verify,
                        "SUPPORTED_RUNNER_SHA256",
                        frozenset({
                            manifest["runtime_bindings"][
                                "runner"]["sha256"],
                        })),
                    mock.patch.object(
                        campaign,
                        "read_canonical_json",
                        side_effect=read_fixture),
                    mock.patch.object(campaign, "_validate_manifest"),
                    mock.patch.object(
                        campaign,
                        "verify_campaign",
                        side_effect=invoke_patched_builder),
            ):
                imported, wrapped = \
                    parallel_verify.verify_with_parallel_replay(
                        Path(campaign.__file__),
                        evidence,
                        replay_workers=1,
                    )
            self.assertIs(imported, campaign)
            self.assertIs(campaign.build_summaries, original_build)
            self.assertTrue(
                campaign.peel_codec._same_typed_json(serial, wrapped))
            for name in ("cell-ledger.jsonl.gz", "summary.json.gz"):
                self.assertEqual(
                    (root / "wrapped" / name).read_bytes(),
                    (root / "serial" / name).read_bytes(),
                )

    def test_import_only_verifier_rejects_unbound_runner(self):
        bindings = synthetic_runtime_bindings()
        bindings["runner"]["path"] = "/tmp/not-the-imported-runner.py"
        with self.assertRaisesRegex(
                campaign.CampaignError,
                "not manifest-bound"):
            parallel_verify._assert_imported_runner_binding(
                campaign, bindings)

    def test_import_only_verifier_rejects_unsupported_runner_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            thermal = Path(temporary) / "thermal.csv"
            thermal.write_text("fixture\n", encoding="ascii")
            bindings = campaign.runtime_bindings(
                sys.executable, thermal)
            bindings["runner"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(
                    campaign.CampaignError,
                    "does not support"):
                parallel_verify._assert_imported_runner_binding(
                    campaign, bindings)

    def test_replay_worker_counts_are_typed_and_bounded(self):
        manifest = {"workers": 4}
        self.assertEqual(
            campaign._replay_worker_count(manifest, None, 10), 4)
        self.assertEqual(
            parallel_verify._worker_count(
                campaign, manifest, None, 10), 4)
        self.assertEqual(
            campaign._replay_worker_count(manifest, 8, 3), 3)
        self.assertEqual(
            parallel_verify._worker_count(
                campaign, manifest, 8, 3), 3)
        for value in (False, 0, 33, 1.0, "2"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                        campaign.CampaignError, "worker count"):
                    campaign._replay_worker_count(
                        manifest, value, 10)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "worker count"):
                    parallel_verify._worker_count(
                        campaign, manifest, value, 10)

    @staticmethod
    def _adversarial_pool_manifest(marker):
        assignments = [
            {"order": index, "job_id": str(index), "output": "unused"}
            for index in range(2)
        ]
        jobs = [
            {
                "job_id": str(index),
                "_test_marker": str(marker),
            }
            for index in range(2)
        ]
        return {
            "workers": 2,
            "bandtiming_protocol": "test",
            "runtime_bindings": {},
            "assignments": assignments,
            "pre_cpu_jobs": jobs,
        }

    def _assert_worker_markers_reaped(self, root):
        markers = tuple(root.glob("worker-*"))
        self.assertTrue(markers)
        for marker in markers:
            pid = int(marker.name.removeprefix("worker-"))
            self.assertFalse(
                Path("/proc", str(pid)).exists(),
                f"replay child {pid} survived cleanup",
            )

    def test_parallel_error_is_fail_fast_and_reaps_campaign_workers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self._adversarial_pool_manifest(
                root / "slow-started")
            ledger_path = root / "ledger.jsonl.gz"
            started = time.monotonic()
            with (
                    mock.patch.object(
                        campaign,
                        "_replay_result",
                        adversarial_replay_failure,
                    ),
                    self.assertRaisesRegex(
                        campaign.CampaignError,
                        "adversarial fast failure",
                    ),
            ):
                with campaign._ordered_replayed_results(
                        root,
                        manifest,
                        "manifest",
                        replay_workers=2,
                ) as replayed:
                    with campaign.AtomicGzipJsonLines(
                            ledger_path) as ledger:
                        ledger.write({"prefix": True})
                        next(replayed)
            self.assertLess(time.monotonic() - started, 1.0)
            self._assert_worker_markers_reaped(root)
            self.assertFalse(ledger_path.exists())
            self.assertFalse(
                ledger_path.with_name(
                    ledger_path.name + ".tmp").exists())

    def test_parallel_error_is_fail_fast_and_reaps_wrapper_workers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self._adversarial_pool_manifest(
                root / "slow-started")
            ledger_path = root / "ledger.jsonl.gz"
            started = time.monotonic()
            with (
                    mock.patch.object(
                        parallel_verify,
                        "_initialize_worker",
                        adversarial_wrapper_initializer,
                    ),
                    mock.patch.object(
                        parallel_verify,
                        "_replay_task",
                        adversarial_wrapper_task,
                    ),
                    self.assertRaisesRegex(
                        campaign.CampaignError,
                        "adversarial fast failure",
                    ),
            ):
                with parallel_verify._ordered_results(
                        campaign,
                        Path(campaign.__file__),
                        root,
                        manifest,
                        "manifest",
                        replay_workers=2,
                ) as replayed:
                    with campaign.AtomicGzipJsonLines(
                            ledger_path) as ledger:
                        ledger.write({"prefix": True})
                        next(replayed)
            self.assertLess(time.monotonic() - started, 1.0)
            self._assert_worker_markers_reaped(root)
            self.assertFalse(ledger_path.exists())
            self.assertFalse(
                ledger_path.with_name(
                    ledger_path.name + ".tmp").exists())

    def test_result_evidence_caps_reject_compressed_and_inflated_bombs(self):
        for bomb in ("compressed", "uncompressed"):
            with (
                    self.subTest(bomb=bomb),
                    mock.patch.object(campaign, "K_VALUES", (2,)),
                    mock.patch.object(
                        campaign, "SELECTION_REPLICATES", 3),
                    tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                evidence = root / "evidence"
                manifest = parallel_replay_fixture(evidence)
                assignment = manifest["assignments"][0]
                job = manifest["pre_cpu_jobs"][0]
                limits = campaign._result_evidence_byte_limits(job)
                result_path = evidence / assignment["output"]
                if bomb == "compressed":
                    result_path.write_bytes(
                        b"x" * (limits["compressed"] + 1))
                    expected = "exceeds compressed byte limit"
                else:
                    payload = b"x" * (limits["uncompressed"] + 1)
                    compressed = gzip.compress(payload, mtime=0)
                    self.assertLess(
                        len(compressed), limits["compressed"])
                    result_path.write_bytes(compressed)
                    expected = "exceeds uncompressed byte limit"

                invocations = (
                    (
                        "serial",
                        lambda output: campaign.build_summaries(
                            evidence,
                            manifest,
                            output,
                            replay_workers=1,
                        ),
                    ),
                    (
                        "parallel",
                        lambda output: campaign.build_summaries(
                            evidence,
                            manifest,
                            output,
                            replay_workers=2,
                        ),
                    ),
                    (
                        "wrapper",
                        lambda output:
                            parallel_verify.parallel_build_summaries(
                                campaign,
                                Path(campaign.__file__),
                                evidence,
                                manifest,
                                output,
                                replay_workers=2,
                            ),
                    ),
                )
                for name, invoke in invocations:
                    output = root / name
                    patches = (
                        mock.patch.object(
                            parallel_verify,
                            "_load_frozen_campaign",
                            return_value=campaign,
                        ),
                        mock.patch.object(
                            parallel_verify,
                            "_assert_source_attestation",
                        ),
                        mock.patch.object(
                            parallel_verify,
                            "SUPPORTED_RUNNER_SHA256",
                            frozenset({
                                manifest["runtime_bindings"][
                                    "runner"]["sha256"],
                            }),
                        ),
                    ) if name == "wrapper" else (
                        nullcontext(),
                        nullcontext(),
                        nullcontext(),
                    )
                    with (
                            self.subTest(path=name),
                            patches[0],
                            patches[1],
                            patches[2],
                            self.assertRaisesRegex(
                                campaign.CampaignError, expected),
                    ):
                        invoke(output)
                    self.assertFalse(
                        (output / "cell-ledger.jsonl.gz").exists())
                    self.assertFalse(
                        (output / "cell-ledger.jsonl.gz.tmp").exists())
                    self.assertFalse(
                        (output / "summary.json.gz").exists())

    def test_replay_parent_signals_reap_children_and_atomic_temps(self):
        campaign_source = r"""
import os
from pathlib import Path
import signal
import sys
import time
sys.path.insert(0, sys.argv[1])
import wh2_za5v_campaign as campaign
root = Path(sys.argv[2])
def slow(directory, manifest, manifest_sha256, assignment, job):
    (root / ("worker-" + str(os.getpid()))).touch()
    time.sleep(30)
    raise AssertionError("sleep completed")
campaign._replay_result = slow
manifest = {
    "workers": 2,
    "bandtiming_protocol": "test",
    "runtime_bindings": {},
    "assignments": [
        {"order": i, "job_id": str(i), "output": "unused"}
        for i in range(2)
    ],
    "pre_cpu_jobs": [{"job_id": str(i)} for i in range(2)],
}
try:
    with campaign._ordered_replayed_results(
            root, manifest, "manifest", replay_workers=2) as replayed:
        with campaign.AtomicGzipJsonLines(
                root / "ledger.jsonl.gz") as ledger:
            ledger.write({"prefix": True})
            next(replayed)
except campaign.CampaignError:
    raise SystemExit(42)
raise SystemExit(99)
"""
        wrapper_source = r"""
import os
from pathlib import Path
import signal
import sys
import time
sys.path.insert(0, sys.argv[1])
import wh2_za5v_campaign as campaign
import wh2_za5v_parallel_verify as wrapper
root = Path(sys.argv[2])
def initialize(*unused):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGHUP, signal.SIG_DFL)
def slow(task):
    (root / ("worker-" + str(os.getpid()))).touch()
    time.sleep(30)
    raise AssertionError("sleep completed")
wrapper._initialize_worker = initialize
wrapper._replay_task = slow
manifest = {
    "workers": 2,
    "bandtiming_protocol": "test",
    "runtime_bindings": {},
    "assignments": [
        {"order": i, "job_id": str(i), "output": "unused"}
        for i in range(2)
    ],
    "pre_cpu_jobs": [{"job_id": str(i)} for i in range(2)],
}
try:
    with wrapper._ordered_results(
            campaign, Path(campaign.__file__), root, manifest,
            "manifest", replay_workers=2) as replayed:
        with campaign.AtomicGzipJsonLines(
                root / "ledger.jsonl.gz") as ledger:
            ledger.write({"prefix": True})
            next(replayed)
except campaign.CampaignError:
    raise SystemExit(42)
raise SystemExit(99)
"""
        for kind, source in (
                ("campaign", campaign_source),
                ("wrapper", wrapper_source),
        ):
            for signum in (signal.SIGINT, signal.SIGTERM):
                with (
                        self.subTest(kind=kind, signal=signum),
                        tempfile.TemporaryDirectory() as temporary,
                ):
                    root = Path(temporary)
                    process = subprocess.Popen(
                        [
                            sys.executable,
                            "-W",
                            "error",
                            "-c",
                            source,
                            str(BENCH),
                            str(root),
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    try:
                        deadline = time.monotonic() + 10.0
                        while len(tuple(root.glob("worker-*"))) < 2:
                            if process.poll() is not None:
                                stdout, stderr = process.communicate()
                                self.fail(
                                    "signal harness exited before worker "
                                    f"start: {stdout!r} {stderr!r}")
                            if time.monotonic() >= deadline:
                                self.fail(
                                    "signal harness worker did not start")
                            time.sleep(0.01)
                        os.kill(process.pid, signum)
                        stdout, stderr = process.communicate(timeout=5)
                        self.assertEqual(
                            process.returncode,
                            42,
                            (stdout, stderr),
                        )
                        self._assert_worker_markers_reaped(root)
                        self.assertFalse(
                            (root / "ledger.jsonl.gz").exists())
                        self.assertFalse(
                            (root / "ledger.jsonl.gz.tmp").exists())
                    finally:
                        if process.poll() is None:
                            process.kill()
                            process.wait()

    def test_source_loader_ignores_stale_timestamp_pyc(self):
        source = r"""
import os
from pathlib import Path
import py_compile
import sys
import tempfile
sys.path.insert(0, sys.argv[1])
import wh2_za5v_parallel_verify as wrapper
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bench = root / "bench"
    tools = root / "tools"
    bench.mkdir()
    tools.mkdir()
    (tools / "peel_codec.py").write_text(
        "class MeasurementError(Exception):\n    pass\n",
        encoding="ascii")
    (tools / "band_timing_codec.py").write_text(
        "import peel_codec\n", encoding="ascii")
    runner = bench / "frozen_runner.py"
    old = (
        "import band_timing_codec as band\n"
        "import peel_codec\n"
        "MARKER = 'loaded_A'\n"
    )
    new = old.replace("loaded_A", "source_B")
    assert len(old) == len(new)
    runner.write_text(old, encoding="ascii")
    saved = runner.stat().st_mtime_ns
    py_compile.compile(str(runner), doraise=True)
    runner.write_text(new, encoding="ascii")
    os.utime(runner, ns=(saved, saved))
    wrapper.SUPPORTED_RUNNER_SHA256 = frozenset({
        wrapper.hashlib.sha256(runner.read_bytes()).hexdigest()})
    wrapper.SUPPORTED_PARSER_SHA256 = frozenset({
        wrapper.hashlib.sha256(
            (tools / "band_timing_codec.py").read_bytes()).hexdigest()})
    wrapper.SUPPORTED_CONTEXT_TOOL_SHA256 = frozenset({
        wrapper.hashlib.sha256(
            (tools / "peel_codec.py").read_bytes()).hexdigest()})
    loaded = wrapper._load_frozen_campaign(runner)
    assert loaded.MARKER == "source_B", loaded.MARKER
"""
        completed = subprocess.run(
            [
                sys.executable,
                "-W",
                "error",
                "-c",
                source,
                str(BENCH),
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(
            completed.returncode,
            0,
            (completed.stdout, completed.stderr),
        )

    def test_source_snapshot_rejects_oversized_runner_and_dependencies(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for kind, limit in parallel_verify.SOURCE_BYTE_LIMITS.items():
                with self.subTest(kind=kind):
                    path = root / f"{kind}.py"
                    with open(path, "wb") as source:
                        source.truncate(limit + 1)
                    with self.assertRaisesRegex(
                            RuntimeError, "exceeds source byte limit"):
                        parallel_verify._source_snapshot(path, kind)
        source = r"""
from pathlib import Path
import resource
import sys
import tempfile
sys.path.insert(0, sys.argv[1])
import wh2_za5v_parallel_verify as wrapper
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bench = root / "bench"
    bench.mkdir()
    runner = bench / "oversized.py"
    with open(runner, "wb") as output:
        output.truncate(128 * 1024 * 1024)
    resource.setrlimit(
        resource.RLIMIT_AS,
        (96 * 1024 * 1024, 96 * 1024 * 1024))
    try:
        wrapper._load_frozen_campaign(runner)
    except RuntimeError as error:
        assert "exceeds source byte limit" in str(error), error
    except MemoryError:
        raise AssertionError("oversized source was read before rejection")
    else:
        raise AssertionError("oversized source was accepted")
"""
        completed = subprocess.run(
            [
                sys.executable,
                "-W",
                "error",
                "-c",
                source,
                str(BENCH),
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(
            completed.returncode,
            0,
            (completed.stdout, completed.stderr),
        )

    def test_source_loader_rejects_preloaded_and_substituted_modules(self):
        preloaded = r"""
from pathlib import Path
import sys
sys.path.insert(0, sys.argv[1])
import wh2_za5v_parallel_verify as wrapper
runner = Path(sys.argv[1]) / "wh2_za5v_campaign.py"
wrapper.SUPPORTED_RUNNER_SHA256 = frozenset({
    wrapper.hashlib.sha256(runner.read_bytes()).hexdigest()})
sys.path.insert(0, str(runner.parents[1] / "tools"))
import peel_codec
peel_codec.strict_json_loads = lambda unused: {"substituted": True}
try:
    wrapper._load_frozen_campaign(runner)
except RuntimeError as error:
    assert "refusing preloaded" in str(error), error
else:
    raise AssertionError("preloaded same-path module was accepted")
"""
        substituted = r"""
from pathlib import Path
import sys
import types
sys.path.insert(0, sys.argv[1])
import wh2_za5v_parallel_verify as wrapper
runner = Path(sys.argv[1]) / "wh2_za5v_campaign.py"
wrapper.SUPPORTED_RUNNER_SHA256 = frozenset({
    wrapper.hashlib.sha256(runner.read_bytes()).hexdigest()})
campaign = wrapper._load_frozen_campaign(runner)
fake = types.ModuleType("peel_codec")
fake.__file__ = campaign.peel_codec.__file__
sys.modules["peel_codec"] = fake
try:
    wrapper._load_frozen_campaign(runner)
except RuntimeError as error:
    assert "substituted" in str(error), error
else:
    raise AssertionError("same-path cached module substitution was accepted")
"""
        for source in (preloaded, substituted):
            completed = subprocess.run(
                [
                    sys.executable,
                    "-W",
                    "error",
                    "-c",
                    source,
                    str(BENCH),
                ],
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(
                completed.returncode,
                0,
                (completed.stdout, completed.stderr),
            )

    def test_source_loader_cleans_partial_modules_after_compile_failure(self):
        source = r"""
from pathlib import Path
import sys
import tempfile
sys.path.insert(0, sys.argv[1])
import wh2_za5v_parallel_verify as wrapper
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bench = root / "bench"
    tools = root / "tools"
    bench.mkdir()
    tools.mkdir()
    context = tools / "peel_codec.py"
    parser = tools / "band_timing_codec.py"
    runner = bench / "frozen_runner.py"
    context.write_text(
        "class MeasurementError(Exception):\n    pass\n",
        encoding="ascii")
    parser.write_text("import peel_codec\n", encoding="ascii")
    runner.write_text(
        "import band_timing_codec as band\n"
        "import peel_codec\n"
        "this is not valid python !!!\n",
        encoding="ascii")
    digest = wrapper.hashlib.sha256
    wrapper.SUPPORTED_RUNNER_SHA256 = frozenset({
        digest(runner.read_bytes()).hexdigest()})
    wrapper.SUPPORTED_PARSER_SHA256 = frozenset({
        digest(parser.read_bytes()).hexdigest()})
    wrapper.SUPPORTED_CONTEXT_TOOL_SHA256 = frozenset({
        digest(context.read_bytes()).hexdigest()})
    try:
        wrapper._load_frozen_campaign(runner)
    except SyntaxError:
        pass
    else:
        raise AssertionError("invalid frozen runner compiled")
    assert "peel_codec" not in sys.modules
    assert "band_timing_codec" not in sys.modules
    assert not any(
        name.startswith("_wirehair_frozen_wh2_za5v_")
        for name in sys.modules)
"""
        completed = subprocess.run(
            [
                sys.executable,
                "-W",
                "error",
                "-c",
                source,
                str(BENCH),
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(
            completed.returncode,
            0,
            (completed.stdout, completed.stderr),
        )


class RecoveryAggregationTests(unittest.TestCase):
    def setUp(self):
        all_jobs = campaign.build_pre_cpu_jobs("selection")
        self.jobs = tuple(sorted(
            (
                job for job in all_jobs
                if job["K"] == 2
            ),
            key=lambda job: (
                job["candidate"],
                campaign.SCHEDULES.index(job["schedule"]),
            ),
        ))
        self.assertEqual(len(self.jobs), 24)
        self.one_candidate_jobs = tuple(
            job for job in self.jobs
            if job["candidate"] == "pure8_s0m1_d3")
        self.assertEqual(len(self.one_candidate_jobs), 6)
        self.timing_jobs = tuple(
            {**job, "order": index}
            for index, job in enumerate(self.jobs))

    def test_full_class_tail_discordance_and_weak_dedup(self):
        aggregator = campaign.RecoveryAggregator("selection", self.jobs)
        for job in self.jobs:
            for replicate in range(job["replicates"]):
                aggregator.add_replicate(
                    job, synthetic_replicate(job, replicate))
        result = aggregator.finish()
        self.assertEqual(result["raw_cells"], 24 * 768)
        self.assertEqual(
            result["independent_construction_weaknesses"]["candidate"][
                "pure8_s0m1_d3"]["count"],
            1,
        )
        self.assertEqual(
            result["independent_construction_weaknesses"]["dispatch"]["count"],
            1,
        )
        comparison = result["independent_construction_weaknesses"][
            "candidate_dispatch_intersections"]["pure8_s0m1_d3"]
        self.assertEqual(comparison["candidate_only"]["count"], 0)
        self.assertEqual(comparison["shared"]["count"], 1)
        self.assertEqual(comparison["control_only"]["count"], 0)
        wh1_comparison = result["independent_construction_weaknesses"][
            "candidate_wh1_intersections"]["pure8_s0m1_d3"]
        self.assertEqual(wh1_comparison["candidate_only"]["count"], 1)
        self.assertEqual(wh1_comparison["shared"]["count"], 0)
        self.assertEqual(wh1_comparison["control_only"]["count"], 0)
        self.assertEqual(
            result["independent_construction_weaknesses"]["wh1"]["count"],
            0,
        )
        group = next(
            item for item in result["groups"]
            if item["candidate_name"] == "pure8_s0m1_d3"
            and item["schedule"] == "iid")
        self.assertEqual(group["raw_cells"], 768)
        self.assertEqual(group["candidate"]["result_classes"], {
            "success": 765,
            "need_more": 1,
            "weak": 1,
            "error": 1,
        })
        self.assertEqual(group["candidate"]["byte_recovery_failures"], 3)
        self.assertEqual(group["wh1"]["result_classes"], {
            "success": 766,
            "need_more": 2,
            "weak": 0,
            "error": 0,
        })
        self.assertEqual(
            group["wh1"]["unrecovered_by_overhead_h_0_to_64"][0],
            768,
        )
        self.assertEqual(
            group["wh1"]["unrecovered_by_overhead_h_0_to_64"][1],
            2,
        )
        self.assertEqual(
            group["candidate"]["unrecovered_by_overhead_h_0_to_64"][0],
            768,
        )
        self.assertEqual(
            group["candidate"]["unrecovered_by_overhead_h_0_to_64"][1],
            4,
        )
        self.assertEqual(
            group["candidate"]["unrecovered_by_overhead_h_0_to_64"][2],
            3,
        )
        thresholds = group[
            "paired_candidate_dispatch_by_overhead_h_0_to_64"]
        self.assertEqual(thresholds[0], {
            "h": 0,
            "denominator": 768,
            "both_recovered": 0,
            "candidate_only_unrecovered": 2,
            "dispatch_only_unrecovered": 0,
            "both_unrecovered": 766,
        })
        self.assertEqual(thresholds[1], {
            "h": 1,
            "denominator": 768,
            "both_recovered": 0,
            "candidate_only_unrecovered": 2,
            "dispatch_only_unrecovered": 764,
            "both_unrecovered": 2,
        })
        self.assertEqual(thresholds[2], {
            "h": 2,
            "denominator": 768,
            "both_recovered": 764,
            "candidate_only_unrecovered": 2,
            "dispatch_only_unrecovered": 1,
            "both_unrecovered": 1,
        })
        pairs = {
            (item["left_arm"], item["right_arm"]): item
            for item in group["paired_recovery"]
        }
        self.assertEqual(
            set(pairs), set(campaign.RECOVERY_PAIRS))
        candidate_wh1 = pairs[("candidate", "wh1")]
        self.assertEqual(candidate_wh1["final"], {
            "both_recovered": 764,
            "left_only_recovered": 1,
            "right_only_recovered": 2,
            "neither_recovered": 1,
        })
        self.assertEqual(
            candidate_wh1["by_overhead_h_0_to_64"][0],
            {
                "h": 0,
                "denominator": 768,
                "both_recovered": 0,
                "left_only_recovered": 0,
                "right_only_recovered": 0,
                "neither_recovered": 768,
            },
        )
        self.assertEqual(
            candidate_wh1["by_overhead_h_0_to_64"][1],
            {
                "h": 1,
                "denominator": 768,
                "both_recovered": 763,
                "left_only_recovered": 1,
                "right_only_recovered": 3,
                "neither_recovered": 1,
            },
        )
        dispatch_wh1 = pairs[("dispatch", "wh1")]
        self.assertEqual(dispatch_wh1["final"], {
            "both_recovered": 765,
            "left_only_recovered": 1,
            "right_only_recovered": 1,
            "neither_recovered": 1,
        })
        self.assertEqual(
            dispatch_wh1["by_overhead_h_0_to_64"][1],
            {
                "h": 1,
                "denominator": 768,
                "both_recovered": 2,
                "left_only_recovered": 0,
                "right_only_recovered": 764,
                "neither_recovered": 2,
            },
        )
        self.assertEqual(group["censored_replicates"], 4)
        self.assertEqual(
            group["censored_panel_observations"][
                "decoder_candidate_dispatch"],
            4,
        )

    def test_duplicate_and_seed_mismatch_are_rejected(self):
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        receipt = synthetic_replicate(self.one_candidate_jobs[0], 0)
        aggregator.add_replicate(self.one_candidate_jobs[0], receipt)
        with self.assertRaisesRegex(
                campaign.CampaignError, "duplicate matched"):
            aggregator.add_replicate(self.one_candidate_jobs[0], receipt)

        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        forged = synthetic_replicate(self.one_candidate_jobs[0], 0)
        forged["loss_seed"] ^= 1
        with self.assertRaisesRegex(
                campaign.CampaignError, "frozen seeds"):
            aggregator.add_replicate(self.one_candidate_jobs[0], forged)

        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        forged = synthetic_replicate(self.one_candidate_jobs[0], 0)
        forged["censored_panels"] = [[]]
        with self.assertRaisesRegex(
                campaign.CampaignError, "panel censor list is invalid"):
            aggregator.add_replicate(self.one_candidate_jobs[0], forged)

    def test_constructor_class_must_match_across_schedules(self):
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(self.one_candidate_jobs[0], 0)
        second = synthetic_replicate(self.one_candidate_jobs[1], 0)
        aggregator.add_replicate(self.one_candidate_jobs[0], first)
        second["candidate_construct_result"] = 8
        second["candidate_construct_class"] = "error"
        with self.assertRaisesRegex(
                campaign.CampaignError,
                "constructor result/class changed across schedules"):
            aggregator.add_replicate(self.one_candidate_jobs[1], second)

        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(self.one_candidate_jobs[0], 2)
        second = synthetic_replicate(self.one_candidate_jobs[1], 2)
        aggregator.add_replicate(self.one_candidate_jobs[0], first)
        self.assertEqual(second["candidate_construct_class"], "weak")
        second["candidate_construct_result"] = 3
        with self.assertRaisesRegex(
                campaign.CampaignError,
                "constructor result/class changed across schedules"):
            aggregator.add_replicate(self.one_candidate_jobs[1], second)

    def test_shared_postconstruction_error_is_not_constructor_error(self):
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(self.one_candidate_jobs[0], 7)
        second = synthetic_replicate(self.one_candidate_jobs[1], 7)
        aggregator.add_replicate(self.one_candidate_jobs[0], first)
        second["encoder_candidate_result_class"] = "error"
        second["encoder_candidate_overhead"] = -1
        second["decoder_candidate_result_class"] = "error"
        second["decoder_candidate_overhead"] = -1
        second["direct_candidate_result_class"] = "error"
        second["direct_candidate_overhead"] = -1
        cell = aggregator.add_replicate(
            self.one_candidate_jobs[1], second)
        self.assertEqual(cell["candidate_construction_class"], "success")

    def test_wh1_constructor_class_must_match_across_schedules(self):
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(self.one_candidate_jobs[0], 8)
        second = synthetic_replicate(self.one_candidate_jobs[1], 8)
        aggregator.add_replicate(self.one_candidate_jobs[0], first)
        second["wh1_construct_result"] = 8
        second["wh1_construct_class"] = "error"
        with self.assertRaisesRegex(
                campaign.CampaignError,
                "WH1 constructor result/class changed"):
            aggregator.add_replicate(self.one_candidate_jobs[1], second)

    def test_control_drift_across_candidates_is_rejected(self):
        same_cell_jobs = tuple(
            job for job in self.jobs
            if job["schedule"] == "iid")[:2]
        self.assertNotEqual(
            same_cell_jobs[0]["candidate"], same_cell_jobs[1]["candidate"])
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(same_cell_jobs[0], 11)
        second = synthetic_replicate(same_cell_jobs[1], 11)
        aggregator.add_replicate(same_cell_jobs[0], first)
        second["decoder_wh1_overhead"] = 2
        with self.assertRaisesRegex(
                campaign.CampaignError, "control changed"):
            aggregator.add_replicate(same_cell_jobs[1], second)

    def test_full_control_scope_drift_is_rejected(self):
        same_cell_jobs = tuple(
            job for job in self.jobs
            if job["schedule"] == "iid")[:2]
        mutations = (
            ("trace_sha256", "0" * 64),
            ("encoder_dispatch_result_class", "error"),
            ("encoder_dispatch_overhead", 1),
            ("decoder_dispatch_result_class", "error"),
            ("decoder_dispatch_overhead", 1),
            ("direct_dispatch_result_class", "error"),
            ("encoder_wh1_result_class", "error"),
            ("encoder_wh1_overhead", 1),
            ("decoder_wh1_result_class", "error"),
            ("decoder_wh1_overhead", 2),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                aggregator = campaign.RecoveryAggregator(
                    "selection", self.jobs)
                first = synthetic_replicate(same_cell_jobs[0], 11)
                second = synthetic_replicate(same_cell_jobs[1], 11)
                aggregator.add_replicate(same_cell_jobs[0], first)
                second[field] = value
                with self.assertRaisesRegex(
                        campaign.CampaignError, "control changed"):
                    aggregator.add_replicate(same_cell_jobs[1], second)

    def test_candidate_coupled_direct_overhead_is_not_shared_control(self):
        same_cell_jobs = tuple(
            job for job in self.jobs
            if job["schedule"] == "iid")[:2]
        self.assertNotEqual(
            same_cell_jobs[0]["candidate"], same_cell_jobs[1]["candidate"])
        aggregator = campaign.RecoveryAggregator(
            "selection", self.jobs)
        first = synthetic_replicate(same_cell_jobs[0], 11)
        second = synthetic_replicate(same_cell_jobs[1], 11)
        first_prefix = max(
            first["decoder_candidate_overhead"],
            first["decoder_dispatch_overhead"],
        )
        first["direct_candidate_overhead"] = first_prefix
        first["direct_dispatch_overhead"] = first_prefix
        aggregator.add_replicate(same_cell_jobs[0], first)

        # The direct candidate/dispatch panel solves both arms at the
        # max first-success prefix when both recover (otherwise K+64).  A
        # different candidate can therefore move both recorded direct
        # overheads without changing the invariant dispatch result class or
        # shared dispatch/WH1 recovery controls.
        second["decoder_candidate_overhead"] = first_prefix + 1
        second["direct_candidate_overhead"] = first_prefix + 1
        second["direct_dispatch_overhead"] = first_prefix + 1
        cell = aggregator.add_replicate(same_cell_jobs[1], second)
        self.assertEqual(cell["dispatch_overhead"],
                         second["decoder_dispatch_overhead"])

    def test_timing_and_thermal_summary_cardinality(self):
        evidence = campaign.JobEvidenceAggregator(
            "selection", self.timing_jobs)
        for index, job in enumerate(self.timing_jobs):
            evidence.add(
                job,
                synthetic_job_evidence(
                    job, index % 4, 60.0 + index / 10.0, 45.0 + index / 20.0,
                ),
            )
        result = evidence.finish()
        self.assertEqual(len(result["timing_jobs"]), 24)
        self.assertEqual(
            len(result["timing_jobs"][0]["arms"]),
            len(campaign.TIMING_ARM_KEYS),
        )
        self.assertEqual(
            len(result["timing_jobs"][0]["contrasts"]), 7)
        self.assertEqual(
            result["timing_jobs"][0]["censoring"][
                "censored_panel_observations"][
                    "decoder_candidate_dispatch"],
            1,
        )
        thermal = result["thermal_context"]
        self.assertEqual(thermal["jobs"], 24)
        self.assertEqual(thermal["rows"], 72)
        self.assertEqual(thermal["valid_rows"], 72)
        self.assertEqual(thermal["missing_read_rows"], 0)
        self.assertEqual(thermal["invalid_rows"], 0)
        self.assertEqual(thermal["edac_ce_max"], 0)
        self.assertEqual(thermal["edac_ue_max"], 0)
        self.assertEqual(thermal["cpu_assignments"], [0, 1, 2, 3])
        self.assertAlmostEqual(thermal["cpu_tctl_max_c"], 62.3)
        self.assertAlmostEqual(thermal["dimm_max_c"], 46.15)

    def test_job_evidence_rejects_unhashable_timing_keys(self):
        job = self.timing_jobs[0]
        cases = (
            ("arm-scope", "arms", 0, "scope"),
            ("arm-name", "arms", 0, "arm"),
            ("contrast-name", "contrasts", 0, "name"),
        )
        for label, section, index, field in cases:
            with self.subTest(label=label):
                receipt = synthetic_job_evidence(
                    job, 0, 60.0, 45.0)
                receipt[section][index][field] = []
                evidence = campaign.JobEvidenceAggregator(
                    "selection", (job,))
                with self.assertRaisesRegex(
                        campaign.CampaignError,
                        "timing receipt summaries are malformed"):
                    evidence.add(job, receipt)

    def test_job_evidence_classifies_weak_and_nonrepairable_regressions(self):
        job = self.timing_jobs[0]
        cases = (
            ("weak", 4, "weak",
             "candidate_constructor_weak_regressions"),
            ("need_more", 1, "need_more",
             "candidate_nonrepairable_regressions"),
        )
        for label, result, result_class, expected_field in cases:
            with self.subTest(label=label):
                evidence = campaign.JobEvidenceAggregator(
                    "selection", (job,))
                receipt = synthetic_job_evidence(
                    job, 0, 60.0, 45.0)
                replicate = receipt["replicates"][0]
                replicate["candidate_construct_result"] = result
                replicate["candidate_construct_class"] = result_class
                for scope in ("encoder", "decoder", "direct"):
                    replicate[
                        f"{scope}_candidate_result_class"] = result_class
                for panel, unused_scope, unused_control in \
                        campaign.CANDIDATE_CONTROL_PANELS:
                    synthetic_contrast(
                        {"contrasts": receipt["contrasts"]},
                        panel,
                    ).update({
                        "eligible_replicates":
                            job["replicates"] - 1,
                        "recovery_regressions": 1,
                    })
                evidence.add(job, receipt)
                timing_job = evidence.finish()["timing_jobs"][0]
                for panel, unused_scope, unused_control in \
                        campaign.CANDIDATE_CONTROL_PANELS:
                    contrast = synthetic_contrast(timing_job, panel)
                    self.assertEqual(contrast[expected_field], 1)
                    other = (
                        "candidate_nonrepairable_regressions"
                        if expected_field
                            == "candidate_constructor_weak_regressions"
                        else "candidate_constructor_weak_regressions"
                    )
                    self.assertEqual(contrast[other], 0)

    def test_frozen_selection_decision_rejects_cherry_picking(self):
        selected = "pure8_s0_d3"
        specs = {
            "pure8_s0m1_d3": {
                "recovery": 1, "weak": 3,
                "direct": -0.10, "encoder": 0.0, "decoder": 0.0,
            },
            selected: {
                "recovery": 1, "weak": 1,
                "direct": -0.20, "encoder": -0.10, "decoder": -0.10,
            },
            "pure8_s0m1_d5": {
                "recovery": 0, "weak": 4,
                "direct": -0.15, "encoder": -0.20, "decoder": -0.20,
            },
            "pure9_s0m1_d3": {
                "recovery": 4, "weak": 0,
                "direct": -0.30, "encoder": -0.30, "decoder": -0.30,
            },
        }
        aggregates, timing = synthetic_decision_inputs(specs)
        decision = campaign.derive_campaign_decision(
            "selection", aggregates, timing)
        self.assertEqual(decision["selected_survivor"], selected)
        self.assertEqual(
            decision["raw_pareto_candidates"],
            ["pure8_s0_d3", "pure8_s0m1_d5", "pure9_s0m1_d3"],
        )
        self.assertEqual(
            decision["eligible_pareto_candidates"],
            ["pure8_s0_d3", "pure8_s0m1_d5"],
        )
        self.assertEqual(
            campaign._require_selected_survivor(decision, selected),
            selected,
        )
        with self.assertRaisesRegex(
                campaign.CampaignError, "differs from frozen selection"):
            campaign._require_selected_survivor(
                decision, "pure8_s0m1_d5")
        malformed = dict(decision)
        malformed["selected_survivor"] = []
        with self.assertRaisesRegex(
                campaign.CampaignError, "unknown survivor"):
            campaign._require_selected_survivor(
                malformed, selected)

    def test_no_eligible_selection_forbids_holdout(self):
        specs = {
            name: {
                "recovery": 4,
                "weak": index,
                "direct": -0.10 - index / 100.0,
                "encoder": -0.05,
                "decoder": -0.05,
            }
            for index, name in enumerate(campaign.CANDIDATE_BY_NAME)
        }
        aggregates, timing = synthetic_decision_inputs(specs)
        decision = campaign.derive_campaign_decision(
            "selection", aggregates, timing)
        self.assertIsNone(decision["selected_survivor"])
        self.assertEqual(decision["status"], "no-survivor")
        with self.assertRaisesRegex(
                campaign.CampaignError, "no survivor"):
            campaign._require_selected_survivor(
                decision, next(iter(campaign.CANDIDATE_BY_NAME)))

    def test_holdout_decision_separates_confirmation_and_promotion(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.03,
                "encoder": -0.03,
                "decoder": -0.01,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        with self.assertRaisesRegex(
                campaign.CampaignError,
                "production blocker evidence is malformed"):
            campaign.derive_campaign_decision(
                "holdout", aggregates, timing, survivor=survivor)
        decision = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        self.assertTrue(decision["holdout_confirmed"])
        self.assertFalse(decision["production_promotion_ready"])
        self.assertEqual(decision["production_action"], "retain-dispatch")
        self.assertFalse(
            decision["wh1_throughput_objectives"]["decoder_confirmed"])

        spec[survivor]["decoder"] = -0.03
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        promoted = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        self.assertTrue(promoted["production_promotion_ready"])
        self.assertEqual(
            promoted["production_action"], "promote-new-contract")

    def test_wide_per_job_confidence_intervals_reject_promotion(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.03,
                "encoder": -0.03,
                "decoder": -0.03,
                "ci_low": -100.0,
                "ci_high": 100.0,
                "aa_ci_low": -100.0,
                "aa_ci_high": 100.0,
                "aa_floor": 100.0,
                "effective_floor": 100.0,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        decision = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        evidence = decision["direct_speed_evidence"]
        self.assertEqual(evidence["log_cost_mean"], -0.03)
        self.assertEqual(evidence["mean_per_job_ci_high"], 100.0)
        self.assertEqual(evidence["mean_effective_floor_log"], 100.0)
        self.assertEqual(evidence["aggregate_upper_log_cost"], 100.0)
        self.assertFalse(evidence["aggregate_confirmation_pass"])
        self.assertFalse(decision["holdout_confirmed"])
        self.assertFalse(decision["production_promotion_ready"])
        self.assertEqual(decision["production_action"], "retain-dispatch")

    def test_aggregate_confirmation_does_not_require_every_job_to_win(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        first = timing["timing_jobs"][0]
        contrast = synthetic_contrast(
            first, "direct_candidate_dispatch")
        contrast.update({
            "log_cost_mean": 0.01,
            "log_cost_ci_low": -0.01,
            "log_cost_ci_high": 0.03,
            "left_faster": False,
        })
        decision = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        evidence = decision["direct_speed_evidence"]
        self.assertEqual(evidence["left_faster_job_count"], 593)
        self.assertLess(
            evidence["aggregate_upper_log_cost"],
            -evidence["mean_effective_floor_log"],
        )
        self.assertTrue(evidence["aggregate_confirmation_pass"])
        self.assertTrue(decision["production_promotion_ready"])

    def test_cross_and_aa_coverage_floors_reject_sparse_timing(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
                "eligible_replicates": 4,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        sparse_cross = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        cross_evidence = sparse_cross["direct_speed_evidence"]
        self.assertEqual(cross_evidence["coverage_floor_per_job"], 128)
        self.assertEqual(
            cross_evidence["cross_minimum_eligible_replicates"], 4)
        self.assertFalse(cross_evidence["evidence_complete"])
        self.assertFalse(sparse_cross["holdout_confirmed"])

        spec[survivor].pop("eligible_replicates")
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        first = timing["timing_jobs"][0]
        aa = synthetic_arm(first, "direct", "candidate")
        aa.update({
            "aa_eligible_replicates": 0,
            "aa_log_cost_mean": 0.0,
            "aa_log_cost_ci_low": None,
            "aa_log_cost_ci_high": None,
            "aa_floor_log": 0.0,
        })
        contrast = synthetic_contrast(
            first, "direct_candidate_dispatch")
        contrast["left_aa_floor_log"] = 0.0
        contrast["left_faster"] = False
        sparse_aa = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        aa_evidence = sparse_aa["direct_speed_evidence"]
        self.assertEqual(aa_evidence["cross_minimum_eligible_replicates"], 256)
        self.assertEqual(aa_evidence["aa_minimum_eligible_replicates"], 0)
        self.assertFalse(aa_evidence["evidence_complete"])
        self.assertFalse(sparse_aa["holdout_confirmed"])

    def test_timing_confirmation_rejects_contradictory_fields(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        mutations = (
            ("nonboolean-left-faster", {"left_faster": 1}),
            (
                "mean-outside-ci",
                {
                    "log_cost_ci_low": -0.03,
                    "log_cost_ci_high": -0.02,
                },
            ),
            ("negative-floor", {"effective_floor_log": -0.01}),
        )
        for label, mutation in mutations:
            with self.subTest(label=label):
                aggregates, timing = synthetic_decision_inputs(
                    spec, phase="holdout")
                synthetic_contrast(
                    timing["timing_jobs"][0],
                    "direct_candidate_dispatch",
                ).update(mutation)
                with self.assertRaises(campaign.CampaignError):
                    campaign.derive_campaign_decision(
                        "holdout", aggregates, timing,
                        survivor=survivor,
                        selection_production_blockers=
                            synthetic_clean_selection_blockers(survivor))

    def test_malformed_decision_hash_keys_raise_campaign_error(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        selection_blockers = synthetic_clean_selection_blockers(survivor)
        mutations = (
            (
                "group-K",
                aggregates["groups"][0],
                "K",
                "decision recovery group is unexpected",
            ),
            (
                "group-schedule",
                aggregates["groups"][0],
                "schedule",
                "decision recovery group is unexpected",
            ),
            (
                "timing-arm-scope",
                timing["timing_jobs"][0]["arms"][0],
                "scope",
                "decision timing evidence is missing",
            ),
            (
                "timing-contrast-name",
                timing["timing_jobs"][0]["contrasts"][0],
                "name",
                "decision timing evidence is missing",
            ),
        )
        for label, target, field, message in mutations:
            with self.subTest(label=label):
                original = target[field]
                target[field] = []
                try:
                    with self.assertRaisesRegex(
                            campaign.CampaignError, message):
                        campaign.derive_campaign_decision(
                            "holdout", aggregates, timing,
                            survivor=survivor,
                            selection_production_blockers=
                                selection_blockers)
                finally:
                    target[field] = original

    def test_candidate_only_regression_in_any_panel_blocks_production(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        for panel, unused_scope, unused_control in \
                campaign.CANDIDATE_CONTROL_PANELS:
            with self.subTest(panel=panel):
                aggregates, timing = synthetic_decision_inputs(
                    spec, phase="holdout")
                contrast = synthetic_contrast(
                    timing["timing_jobs"][0], panel)
                contrast.update({
                    "eligible_replicates": 255,
                    "recovery_regressions": 1,
                    "candidate_nonrepairable_regressions": 1,
                    "left_faster": False,
                })
                decision = campaign.derive_campaign_decision(
                    "holdout", aggregates, timing, survivor=survivor,
                    selection_production_blockers=
                        synthetic_clean_selection_blockers(survivor))
                self.assertTrue(decision["holdout_confirmed"])
                self.assertFalse(
                    decision["production_promotion_ready"])
                self.assertEqual(
                    decision["production_action"],
                    "investigate-reliability",
                )
                cumulative = decision["production_reliability"][
                    "cumulative"]
                self.assertEqual(
                    cumulative["candidate_control_panels"][panel][
                        "candidate_nonrepairable_regressions"],
                    1,
                )

    def test_shared_failure_is_not_a_candidate_specific_veto(self):
        survivor = "pure8_s0_d3"
        spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            spec, phase="holdout")
        contrast = synthetic_contrast(
            timing["timing_jobs"][0],
            "decoder_candidate_wh1",
        )
        contrast["both_failures"] = 1
        contrast["eligible_replicates"] = 255
        decision = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        self.assertEqual(
            decision["wh1_throughput_objectives"][
                "decoder_evidence"]["both_failures"],
            1,
        )
        self.assertTrue(decision["production_promotion_ready"])
        self.assertEqual(
            decision["production_action"], "promote-new-contract")

    def test_selection_and_holdout_blockers_are_cumulative_and_routed(self):
        survivor = "pure8_s0_d3"
        clean_spec = {
            survivor: {
                "recovery": 1,
                "weak": 0,
                "direct": -0.04,
                "encoder": -0.04,
                "decoder": -0.04,
            },
        }
        aggregates, timing = synthetic_decision_inputs(
            clean_spec, phase="holdout")

        selection_weak = synthetic_clean_selection_blockers(survivor)
        selection_weak["candidate_weak_constructions"] = 1
        selection_weak["repairable_seed_blocker"] = True
        weak = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=selection_weak)
        self.assertTrue(weak["holdout_confirmed"])
        self.assertFalse(weak["production_promotion_ready"])
        self.assertEqual(
            weak["production_action"], "repair-seeds-and-retest")
        self.assertEqual(
            weak["production_reliability"]["cumulative"][
                "candidate_weak_constructions"],
            1,
        )

        selection_bad = synthetic_clean_selection_blockers(survivor)
        panel = selection_bad["candidate_control_panels"][
            "decoder_candidate_dispatch"]
        panel["recovery_regressions"] = 1
        panel["candidate_nonrepairable_regressions"] = 1
        selection_bad["nonrepairable_reliability_blocker"] = True
        bad = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=selection_bad)
        self.assertTrue(bad["holdout_confirmed"])
        self.assertFalse(bad["production_promotion_ready"])
        self.assertEqual(
            bad["production_action"], "investigate-reliability")

        selection_error = synthetic_clean_selection_blockers(survivor)
        selection_error["candidate_error_outcomes"]["decoder"] = 1
        selection_error["candidate_error_outcomes"]["total"] = 1
        selection_error["nonrepairable_reliability_blocker"] = True
        shared_error = campaign.derive_campaign_decision(
            "holdout", aggregates, timing, survivor=survivor,
            selection_production_blockers=selection_error)
        self.assertEqual(
            shared_error["production_action"],
            "investigate-reliability",
        )

        holdout_weak_spec = json.loads(json.dumps(clean_spec))
        holdout_weak_spec[survivor]["weak"] = 1
        weak_aggregates, weak_timing = synthetic_decision_inputs(
            holdout_weak_spec, phase="holdout")
        holdout_weak = campaign.derive_campaign_decision(
            "holdout", weak_aggregates, weak_timing,
            survivor=survivor,
            selection_production_blockers=
                synthetic_clean_selection_blockers(survivor))
        self.assertEqual(
            holdout_weak["production_action"],
            "repair-seeds-and-retest",
        )

        for failed_gate in ("direct", "encoder", "decoder", "recovery"):
            with self.subTest(failed_gate=failed_gate):
                slow_spec = json.loads(json.dumps(clean_spec))
                slow_spec[survivor][failed_gate] = (
                    4 if failed_gate == "recovery" else -0.01
                )
                slow_aggregates, slow_timing = synthetic_decision_inputs(
                    slow_spec, phase="holdout")
                slow = campaign.derive_campaign_decision(
                    "holdout", slow_aggregates, slow_timing,
                    survivor=survivor,
                    selection_production_blockers=selection_weak)
                self.assertFalse(slow["production_promotion_ready"])
                self.assertEqual(
                    slow["production_action"], "retain-dispatch")

    def test_thermal_sensor_or_edac_error_is_rejected(self):
        evidence = campaign.JobEvidenceAggregator(
            "selection", self.timing_jobs)
        receipt = synthetic_job_evidence(
            self.timing_jobs[0], 0, 60.0, 45.0)
        receipt["context"]["thermal"]["edac_ce_max"] = 1
        with self.assertRaisesRegex(
                campaign.CampaignError, "sensor or EDAC"):
            evidence.add(self.timing_jobs[0], receipt)

    def test_holdout_emits_only_the_tested_survivor(self):
        holdout_jobs = tuple(
            job for job in campaign.build_pre_cpu_jobs(
                "holdout", "pure8_s0m1_d3")
            if job["K"] == 2
        )
        self.assertEqual(len(holdout_jobs), 6)
        aggregator = campaign.RecoveryAggregator(
            "holdout", holdout_jobs)
        for job in holdout_jobs:
            for replicate in range(job["replicates"]):
                aggregator.add_replicate(
                    job, synthetic_replicate(job, replicate))
        result = aggregator.finish()
        weakness = result["independent_construction_weaknesses"]
        self.assertEqual(
            set(weakness["candidate"]), {"pure8_s0m1_d3"})
        self.assertEqual(
            set(weakness["candidate_dispatch_intersections"]),
            {"pure8_s0m1_d3"},
        )
        self.assertEqual(
            set(weakness["candidate_wh1_intersections"]),
            {"pure8_s0m1_d3"},
        )
        self.assertEqual(
            result["matched_control_cells"]["copies_per_cell"], 1)


if __name__ == "__main__":
    unittest.main()
