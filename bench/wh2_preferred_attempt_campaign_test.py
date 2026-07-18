#!/usr/bin/env python3
"""Deterministic unit fixtures for the H12 preferred-attempt campaign."""

from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path
import signal
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock
from typing import List, Optional, Set, Tuple


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import wh2_preferred_attempt_campaign as campaign
import wh2_preferred_attempt_search as search
import wh2_rank_floor_two_anchor_screen as screen


CONTEXT = "12" * 32
PARENT = "34" * 32
SEEDS = ("0x0000000000000001", "0x0000000000000002")
TEST_PREPARE_BASELINE_ROW_SHA256 = "ab" * 32


def thermal_row(
    monotonic_s: float,
    *,
    edac_ce: int = 0,
    edac_ue: int = 0,
    cpu_tctl_c: float = 60.0,
) -> bytes:
    fields = [
        "2026-07-18T00:00:00.000Z", f"{monotonic_s:.6f}",
        "100.0", "3600.0", f"{cpu_tctl_c:.1f}",
        *(f"{50.0 + index:.2f}" for index in range(8)),
        "0", "128.0", "128.0", "128.0", str(edac_ce), str(edac_ue),
    ]
    if len(fields) != len(screen.THERMAL_FIELDS):
        raise AssertionError("thermal fixture field count changed")
    return (",".join(fields) + "\n").encode("ascii")


def metrics(
    result: int = 0,
    heavy_shortfall: int = 0,
    inactivated: int = 5,
    binary_deficit: int = 0,
    heavy_gain: int = 0,
    block_xors: int = 100,
    block_muladds: int = 10,
) -> campaign.RecoveryMetrics:
    return campaign.RecoveryMetrics(
        result=result,
        rank_fail=int(result == 1),
        error=int(result not in (0, 1)),
        heavy_shortfall=heavy_shortfall,
        inactivated=inactivated,
        binary_deficit=binary_deficit,
        heavy_gain=heavy_gain,
        block_xors=block_xors,
        block_muladds=block_muladds,
    )


def metric_fields(value: campaign.RecoveryMetrics) -> List[str]:
    return [
        str(value.result), str(value.rank_fail), str(value.error),
        str(value.heavy_shortfall), str(value.inactivated),
        str(value.binary_deficit), str(value.heavy_gain),
        str(value.block_xors), str(value.block_muladds),
    ]


def route_record(
    K: int,
    width: int,
    p: int,
    a0: int = 1,
    canonical_probes: int = 0,
    valid: int = 1,
) -> campaign.RouteRecord:
    noop = int(p == a0)
    if p < a0:
        valid = 0
    fallback = int(not valid)
    direct = int(bool(valid) and not bool(noop))
    actual = p if direct else a0
    return campaign.RouteRecord(
        K, width, "preferred", p, a0, actual, valid, fallback, noop,
        direct, canonical_probes, int(p > a0),
    )


def route_output(job: campaign.JobSpec) -> bytes:
    lines = [
        campaign.expected_route_preamble(str(job.route_context_sha256)),
        ",".join(campaign.ROUTE_HEADER),
    ]
    charged: Set[Tuple[int, int]] = set()
    for K in job.Ks:
        for width in job.widths:
            for p in job.attempt_map()[K]:
                key = (K, width)
                record = route_record(
                    K, width, p, canonical_probes=2 if key not in charged else 0)
                charged.add(key)
                lines.append(",".join(str(value) for value in (
                    record.K, record.width, record.route_status,
                    record.preferred_attempt, record.canonical_attempt,
                    record.actual_attempt, record.preferred_valid,
                    record.fallback, record.no_op, record.direct,
                    record.canonical_probe_solves,
                    record.preferred_probe_solves,
                )))
    return ("\n".join(lines) + "\n").encode("ascii")


def recovery_row(
    K: int,
    width: int,
    arm: str,
    p: int,
    value: campaign.RecoveryMetrics,
    a0: int = 1,
) -> str:
    if arm == "control":
        prefix = [K, width, arm, -1, a0, a0, 0, 1, 0, 0, 0, 1]
    else:
        prefix = [K, width, arm, p, a0, p, 1, 1, 0, 0, 1, 1]
    return ",".join([str(item) for item in prefix] + metric_fields(value))


def control_output(
    job: campaign.JobSpec,
    value: Optional[campaign.RecoveryMetrics] = None,
) -> bytes:
    value = metrics() if value is None else value
    lines = [campaign.expected_recovery_preamble(job),
             ",".join(campaign.RECOVERY_HEADER)]
    for K in job.Ks:
        for width in job.widths:
            lines.append(
                "# preferred_probe_accounting: N={} bb={} "
                "canonical_probe_solves=2 preferred_probe_solves=0".format(
                    K, width))
            lines.append(recovery_row(K, width, "control", -1, value))
    return ("\n".join(lines) + "\n").encode("ascii")


def candidate_output(job: campaign.JobSpec) -> bytes:
    lines = [campaign.expected_recovery_preamble(job),
             ",".join(campaign.RECOVERY_HEADER)]
    for K in job.Ks:
        for width in job.widths:
            for p in job.attempt_map()[K]:
                route = route_record(K, width, p)
                if route.direct:
                    lines.append(recovery_row(K, width, "candidate", p, metrics()))
                else:
                    lines.append(
                        "# preferred_candidate_alias: N={} bb={} p={} a0={} "
                        "actual={} valid={} fallback={} no_op={} direct=0 "
                        "physical_solve=0".format(
                            K, width, p, route.canonical_attempt,
                            route.actual_attempt, route.preferred_valid,
                            route.fallback, route.no_op))
    return ("\n".join(lines) + "\n").encode("ascii")


def make_result_root(temporary: str) -> Tuple[Path, Path]:
    root = Path(temporary).resolve()
    binary = root / "frozen" / "wirehair_v2_bench"
    binary.parent.mkdir(parents=True)
    binary.write_bytes(b"mock frozen binary\n")
    binary.chmod(0o755)
    return root, binary


def trust_patch(contract: Optional[dict] = None) -> mock._patch:
    return mock.patch.object(
        search, "verify_frozen_controller_runtime",
        return_value=({}, {} if contract is None else contract))


def runner_contract(
    cpus: Tuple[int, ...], workers: int, timeout: float, thermal: Path,
) -> dict:
    thermal.write_bytes(b"mock live thermal source\n")
    taskset = thermal.parent / "frozen" / "taskset"
    taskset.parent.mkdir(parents=True, exist_ok=True)
    taskset.write_bytes(b"mock frozen taskset\n")
    taskset.chmod(0o755)
    binary = thermal.parent / "frozen" / "wirehair_v2_bench"
    return {
        "cpu_set": list(cpus), "workers": workers,
        "timeout_seconds": timeout, "thermal": str(thermal.resolve()),
        "binary_sha256": campaign.sha256_bytes(binary.read_bytes()),
        "thermal_baseline": {
            "row_sha256": TEST_PREPARE_BASELINE_ROW_SHA256,
        },
        "thermal_policy": {
            "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
            "timing_cpu_limit_c": 85.0, "timing_dimm_limit_c": 90.0,
            "consecutive_samples": 3, "stale_seconds": 5.0,
            "min_cpu_busy_pct": 95.0,
            "edac_ce_delta": 0, "edac_ue_delta": 0,
        },
        "system_tools": {
            "taskset": {
                "path": str(taskset.resolve()),
                "sha256": campaign.sha256_bytes(taskset.read_bytes()),
            },
        },
    }


class FakeThermalGuard:
    instances = 0

    def __init__(self, path: Path, abort: threading.Event) -> None:
        type(self).instances += 1
        self.path = path
        self.abort = abort
        self.started = False
        phase_offset = type(self).instances * 10
        self.baseline_row = thermal_row(100.0 + phase_offset)
        self.interval_row = thermal_row(101.0 + phase_offset)

    def start(self) -> None:
        self.started = True

    def finish(self, output: Path) -> dict:
        if not self.started:
            raise AssertionError("fake thermal worker did not start")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(
            (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii") +
            self.baseline_row + self.interval_row)
        return {
            "samples": 1, "sealed_samples_including_baseline": 2,
            "cpu_busy_min_pct": 100.0, "cpu_tctl_max_c": 60.0,
            "dimm_max_c": 57.0, "dimm_read_errors_max": 0,
            "edac_ce_delta": 0, "edac_ue_delta": 0,
            "thermal_limit_c": 90.0, "thermal_high_samples": 0,
            "thermal_high_max_consecutive_samples": 0,
            "guard_poll_iterations": 1, "guard_samples": 1,
            "guard_high_samples": 0, "guard_limit_c": 90.0,
            "guard_error": None,
            "baseline_row_sha256": campaign.sha256_bytes(self.baseline_row),
        }


class ThermalSourceTests(unittest.TestCase):
    def test_real_thermal_guard_finish_schema_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, _binary = make_result_root(temporary)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            baseline_row = thermal_row(100.0)
            interval_row = thermal_row(101.0)
            encoded = (
                (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii") +
                baseline_row + interval_row)
            base_summary = {
                "samples": 1, "sealed_samples_including_baseline": 2,
                "cpu_busy_min_pct": 100.0, "cpu_tctl_max_c": 60.0,
                "dimm_max_c": 57.0, "dimm_read_errors_max": 0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
                "thermal_limit_c": 90.0, "thermal_high_samples": 0,
                "thermal_high_max_consecutive_samples": 0,
            }

            def finish_interval(
                _path: Path, _mark: dict, output: Path, **_kwargs: object,
            ) -> dict:
                output.write_bytes(encoded)
                return dict(base_summary)

            guard = object.__new__(campaign.common.ThermalGuard)
            guard.started = True
            guard.stop_event = threading.Event()
            guard.thread = mock.Mock()
            guard.path = thermal
            guard.mark = {"baseline_row": baseline_row}
            guard.limit_c = 90.0
            guard.stale_seconds = 5.0
            guard.poll_iterations = 1
            guard.samples = 1
            guard.high_samples = 0
            guard.error = None
            interval = root / "thermal-interval.csv"
            with mock.patch.object(
                    campaign.common, "thermal_finish",
                    side_effect=finish_interval):
                summary = guard.finish(interval)

            self.assertEqual(
                set(summary), campaign.THERMAL_SUMMARY_FIELDS)
            self.assertEqual(
                campaign.validate_thermal_summary(
                    summary, contract["thermal_policy"],
                    campaign.thermal_artifact_baseline_row_sha256(encoded)),
                summary)

    def test_thermal_summary_binds_exact_sealed_phase_baseline_row(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, _binary = make_result_root(temporary)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            guard = FakeThermalGuard(thermal, threading.Event())
            guard.start()
            interval = root / "thermal-interval.csv"
            summary = guard.finish(interval)
            policy = contract["thermal_policy"]
            phase_baseline = \
                campaign.thermal_artifact_baseline_row_sha256(
                    interval.read_bytes())

            with self.assertRaisesRegex(
                    campaign.CampaignError, "canonical baseline"):
                campaign.thermal_artifact_baseline_row_sha256(
                    b"not a thermal artifact\n")

            self.assertNotEqual(
                phase_baseline, TEST_PREPARE_BASELINE_ROW_SHA256)
            self.assertEqual(
                campaign.validate_thermal_summary(
                    summary, policy, phase_baseline),
                summary)

            missing = dict(summary)
            del missing["baseline_row_sha256"]
            with self.assertRaisesRegex(
                    campaign.CampaignError, "fields are not canonical"):
                campaign.validate_thermal_summary(
                    missing, policy, phase_baseline)

            for digest in ("not-a-digest", "cd" * 32):
                with self.subTest(digest=digest):
                    changed = dict(summary)
                    changed["baseline_row_sha256"] = digest
                    with self.assertRaisesRegex(
                            campaign.CampaignError, "baseline binding"):
                        campaign.validate_thermal_summary(
                            changed, policy, phase_baseline)

    def test_timing_history_uses_stricter_cpu_limit_between_panels(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline_row = thermal_row(100.0, cpu_tctl_c=80.0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(header + baseline_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = campaign.common.thermal_start(
                    path, require_zero_edac=False)
            contract = {
                "thermal_baseline": {
                    "dev": mark["dev"], "ino": mark["ino"],
                    "offset": mark["offset"], "edac_ce": 0, "edac_ue": 0,
                    "monotonic_s": mark["monotonic_s"],
                    "max_temperature_c": mark["max_temperature_c"],
                    "row_sha256": campaign.sha256_bytes(baseline_row),
                },
            }
            policy = {
                "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                "timing_cpu_limit_c": 85.0,
                "timing_dimm_limit_c": 90.0,
                "consecutive_samples": 3, "stale_seconds": 5.0,
            }
            with path.open("ab") as output:
                for monotonic_s in (101.0, 102.0, 103.0):
                    output.write(thermal_row(
                        monotonic_s, cpu_tctl_c=86.0))
            with mock.patch.object(screen, "validate_thermal_current"):
                campaign.validate_frozen_thermal_source(
                    contract, path, policy)
                mark = campaign.common.thermal_start(
                    path, require_zero_edac=False)
            timing_contract = {
                "thermal_baseline": {
                    "dev": mark["dev"], "ino": mark["ino"],
                    "offset": mark["offset"], "edac_ce": mark["edac_ce"],
                    "edac_ue": mark["edac_ue"],
                    "monotonic_s": mark["monotonic_s"],
                    "max_temperature_c": mark["max_temperature_c"],
                    "row_sha256": campaign.sha256_bytes(
                        mark["baseline_row"]),
                },
            }
            with mock.patch.object(screen, "validate_thermal_current"):
                campaign.validate_frozen_thermal_source(
                    timing_contract, path, policy, timing=True)
            with path.open("ab") as output:
                for monotonic_s in (104.0, 105.0, 106.0):
                    output.write(thermal_row(
                        monotonic_s, cpu_tctl_c=86.0))
            with mock.patch.object(screen, "validate_thermal_current"):
                with self.assertRaisesRegex(
                        campaign.CampaignError, "consecutive limit"):
                    campaign.validate_frozen_thermal_source(
                        timing_contract, path, policy, timing=True)

    def test_frozen_source_retains_nonzero_edac_baseline_exactly(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline_row = thermal_row(100.0, edac_ce=5, edac_ue=2)
        current_row = thermal_row(101.0, edac_ce=5, edac_ue=2)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(header + baseline_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = campaign.common.thermal_start(
                    path, require_zero_edac=False)
            contract = {
                "thermal_baseline": {
                    "dev": mark["dev"], "ino": mark["ino"],
                    "offset": mark["offset"], "edac_ce": 5, "edac_ue": 2,
                    "monotonic_s": mark["monotonic_s"],
                    "max_temperature_c": mark["max_temperature_c"],
                    "row_sha256": campaign.sha256_bytes(baseline_row),
                },
            }
            with path.open("ab") as output:
                output.write(current_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                self.assertEqual(
                    campaign.validate_frozen_thermal_source(
                        contract, path, {"stale_seconds": 5.0}),
                    (mark["dev"], mark["ino"], 5, 2))
                with self.assertRaisesRegex(
                        campaign.CampaignError, "identity changed"):
                    campaign.common.ThermalGuard(
                        path, threading.Event(),
                        expected_dev=mark["dev"],
                        expected_ino=mark["ino"] + 1,
                        expected_edac_ce=5, expected_edac_ue=2)
                guard = campaign.common.ThermalGuard(
                    path, threading.Event(),
                    expected_dev=mark["dev"], expected_ino=mark["ino"],
                    expected_edac_ce=5, expected_edac_ue=2)
            guard._validate_sample(
                screen.parse_thermal_sample(
                    current_row, "unchanged EDAC fixture"),
                "unchanged EDAC fixture")
            with self.assertRaises(campaign.CampaignError):
                guard._validate_sample(
                    screen.parse_thermal_sample(
                        thermal_row(102.0, edac_ce=6, edac_ue=2),
                        "changed EDAC fixture"),
                    "changed EDAC fixture")

            changed_baseline = baseline_row.replace(b",60.0,", b",61.0,", 1)
            self.assertEqual(len(changed_baseline), len(baseline_row))
            with path.open("r+b") as output:
                output.seek(len(header))
                output.write(changed_baseline)
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(
                     campaign.CampaignError, "baseline row changed"):
                campaign.validate_frozen_thermal_source(
                    contract, path, {"stale_seconds": 5.0})
            with path.open("r+b") as output:
                output.seek(len(header))
                output.write(baseline_row)

            with path.open("ab") as output:
                output.write(thermal_row(102.0, edac_ce=6, edac_ue=2))
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(
                     campaign.CampaignError, "changed since"):
                campaign.validate_frozen_thermal_source(
                    contract, path, {"stale_seconds": 5.0})
            with path.open("ab") as output:
                output.write(thermal_row(103.0, edac_ce=5, edac_ue=2))
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(
                     campaign.CampaignError, "history changed"):
                campaign.validate_frozen_thermal_source(
                    contract, path, {"stale_seconds": 5.0})

    def test_frozen_source_rejects_prelaunch_sampling_gap(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline_row = thermal_row(100.0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(header + baseline_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = campaign.common.thermal_start(
                    path, require_zero_edac=False)
            contract = {
                "thermal_baseline": {
                    "dev": mark["dev"], "ino": mark["ino"],
                    "offset": mark["offset"], "edac_ce": 0, "edac_ue": 0,
                    "monotonic_s": mark["monotonic_s"],
                    "max_temperature_c": mark["max_temperature_c"],
                    "row_sha256": campaign.sha256_bytes(baseline_row),
                },
            }
            with path.open("ab") as output:
                output.write(thermal_row(106.0))
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(campaign.CampaignError, "contains a gap"):
                campaign.validate_frozen_thermal_source(
                    contract, path, {"stale_seconds": 5.0})


class LowBusyThermalGuard(FakeThermalGuard):
    def finish(self, output: Path) -> dict:
        summary = super().finish(output)
        summary["cpu_busy_min_pct"] = 0.0
        return summary


class ParserAndLedgerTests(unittest.TestCase):
    def test_exact_route_and_candidate_csv_accounting(self) -> None:
        route_job = campaign.JobSpec(
            "r1-route-00000", "route", "r1", 0, 0, (10,),
            ((10, (0, 1, 2, 3)),), (64,), None, None, None, None,
            CONTEXT, expected_logical_rows=4,
        )
        route_bytes = route_output(route_job)
        parsed_routes = campaign.parse_route_output(route_job, route_bytes)
        self.assertEqual([row.direct for row in parsed_routes], [0, 0, 1, 1])
        self.assertEqual(sum(row.canonical_probe_solves for row in parsed_routes), 2)
        self.assertEqual(sum(row.preferred_probe_solves for row in parsed_routes), 2)

        route_hash = campaign.sha256_bytes(route_bytes)
        candidate_job = campaign.JobSpec(
            "r1-candidate-00000", "candidate", "r1", 0, 0, (10,),
            ((10, (0, 1, 2, 3)),), (64,), 0, SEEDS[0], "burst", "0.50",
            CONTEXT, "jobs/r1/route/stdout/r1-route-00000.csv", route_hash,
            route_job.job_id, 4, 2,
        )
        parsed = campaign.parse_recovery_output(
            candidate_job, candidate_output(candidate_job), parsed_routes)
        self.assertEqual(len(parsed.rows), 2)
        self.assertEqual(parsed.alias_rows, 2)
        self.assertEqual(parsed.physical_rows, 2)

        tampered = route_bytes.replace(
            b"10,64,preferred,1,1,1,1,0,1,0,0,0",
            b"10,64,preferred,1,1,1,0,1,1,0,0,0")
        with self.assertRaises(campaign.CampaignError):
            campaign.parse_route_output(route_job, tampered)
        with self.assertRaises(campaign.CampaignError):
            campaign.parse_recovery_output(
                candidate_job, candidate_output(candidate_job) + b"# extra\n",
                parsed_routes)

        noncanonical_alias = candidate_output(candidate_job).replace(
            b"N=10 bb=64 p=0", b"N=010 bb=64 p=0")
        with self.assertRaises(campaign.CampaignError):
            campaign.parse_recovery_output(
                candidate_job, noncanonical_alias, parsed_routes)

        reordered = campaign.JobSpec(
            "r1-route-00001", "route", "r1", 0, 1, (10, 11),
            ((11, (2,)), (10, (2,))), (64,), None, None, None, None,
            CONTEXT, expected_logical_rows=2)
        with self.assertRaises(campaign.CampaignError):
            reordered.validate()

    def test_canonical_ledgers_and_exact_row_arithmetic(self) -> None:
        cohort = (10, 11)
        bins = {0: (10,), 1: (11,)}
        control = campaign.build_control_ledger(
            bins, cohort, SEEDS, PARENT, schedules=("burst",),
            losses=("0.50",), widths=(64, 256))
        self.assertEqual(len(control.jobs), 4)
        self.assertEqual(control.record()["expected_logical_rows"], 8)
        self.assertEqual(control.bytes(), control.bytes())

        spec = campaign.RoundSpec("r1", 4, 2, (0,), (64,))
        survivors = {K: (0, 1, 2, 3) for K in cohort}
        route = campaign.build_route_ledger(
            spec, bins, cohort, survivors, PARENT, CONTEXT, batch_size=2)
        self.assertEqual(len(route.jobs), 4)
        self.assertEqual(route.record()["expected_logical_rows"], 8)
        artifacts = {}
        for job in route.jobs:
            data = route_output(job)
            artifacts[job.job_id] = campaign.RouteArtifact(
                job.job_id,
                "jobs/r1/route/stdout/{}.csv".format(job.job_id),
                campaign.sha256_bytes(data),
                campaign.parse_route_output(job, data),
            )
        candidate = campaign.build_candidate_ledger(
            spec, route, artifacts, SEEDS, route.sha256(),
            schedules=("burst",), losses=("0.50",))
        self.assertEqual(len(candidate.jobs), 4)
        self.assertEqual(candidate.record()["expected_logical_rows"], 8)
        self.assertEqual(candidate.record()["expected_physical_rows"], 4)

        changed = dict(artifacts)
        second_job = route.jobs[1]
        changed_records = tuple(
            route_record(
                row.K, row.width, row.preferred_attempt, a0=2,
                canonical_probes=3 if index == 0 else 0)
            for index, row in enumerate(artifacts[second_job.job_id].records)
        )
        changed[second_job.job_id] = campaign.RouteArtifact(
            second_job.job_id,
            "jobs/r1/route/stdout/{}.csv".format(second_job.job_id),
            "56" * 32, changed_records)
        with self.assertRaises(campaign.CampaignError):
            campaign.build_candidate_ledger(
                spec, route, changed, SEEDS, route.sha256(),
                schedules=("burst",), losses=("0.50",))

        empty = {K: () for K in cohort}
        empty_route = campaign.build_route_ledger(
            campaign.RoundSpec("r2", 2, 1, (1,), (64, 1280)),
            bins, cohort, empty, PARENT, CONTEXT)
        self.assertEqual(empty_route.jobs, ())
        self.assertEqual(empty_route.record()["job_count"], 0)

        duplicate_jobs = (
            campaign.JobSpec(
                "r1-route-00000", "route", "r1", 0, 0, (10,),
                ((10, (2,)),), (64,), None, None, None, None,
                CONTEXT, expected_logical_rows=1),
            campaign.JobSpec(
                "r1-route-00001", "route", "r1", 0, 1, (10,),
                ((10, (2,)),), (64,), None, None, None, None,
                CONTEXT, expected_logical_rows=1),
        )
        with self.assertRaises(campaign.CampaignError):
            campaign.JobLedger("route", "r1", PARENT, duplicate_jobs).record()


class RunnerTests(unittest.TestCase):
    def test_missing_guard_baseline_never_seals_phase(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))

            def execute(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"", 100, 200, cpu)

            class MissingBaselineGuard(FakeThermalGuard):
                def finish(self, output: Path) -> dict:
                    summary = super().finish(output)
                    del summary["baseline_row_sha256"]
                    return summary

            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            runner = campaign.JobRunner(
                root, binary, (3,), 1, 10.0, execute,
                MissingBaselineGuard)
            with trust_patch(contract), self.assertRaisesRegex(
                    campaign.CampaignError, "fields are not canonical"):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertTrue(
                campaign.job_paths(
                    root, ledger.jobs[0])["receipt"].exists())
            self.assertFalse(
                campaign.phase_completion_path(root, ledger).exists())

    def test_later_phase_binds_its_local_thermal_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            cohort = (10,)
            bins = {0: cohort}
            survivors = {10: (2,)}
            ledgers = tuple(
                campaign.build_route_ledger(
                    campaign.RoundSpec(name, 4, 2, (0,), (64,)),
                    bins, cohort, survivors, PARENT, CONTEXT)
                for name in ("r1", "r2"))
            next_job = iter(ledger.jobs[0] for ledger in ledgers)

            def execute(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                job = next(next_job)
                return campaign.ExecutionResult(
                    0, route_output(job), b"", 100, 200, cpu)

            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            runner = campaign.JobRunner(
                root, binary, (3,), 1, 10.0, execute, FakeThermalGuard)
            with trust_patch(contract):
                for ledger in ledgers:
                    runner.run_ledger(
                        ledger, thermal, install_signal_handlers=False)

            phase_digests = []
            for ledger in ledgers:
                completion = campaign.load_canonical_object(
                    campaign.phase_completion_path(root, ledger),
                    "test phase completion")
                summary_digest = completion[
                    "thermal_artifact"]["summary"]["baseline_row_sha256"]
                artifact_digest = \
                    campaign.thermal_artifact_baseline_row_sha256(
                        campaign.phase_thermal_path(
                            root, ledger).read_bytes())
                self.assertEqual(summary_digest, artifact_digest)
                self.assertNotEqual(
                    summary_digest, TEST_PREPARE_BASELINE_ROW_SHA256)
                phase_digests.append(summary_digest)
            self.assertNotEqual(phase_digests[0], phase_digests[1])

    def test_subprocess_executor_bounds_both_streams_and_reaps(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-preferred-bounded-output-") as temporary:
            root = Path(temporary)
            for descriptor in (1, 2):
                with self.subTest(descriptor=descriptor):
                    marker = root / "{}.pid".format(descriptor)
                    writer = (
                        "import os, pathlib, signal; "
                        "pathlib.Path({!r}).write_text(str(os.getpid()), "
                        "encoding='ascii'); "
                        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
                        "descriptor={}; chunk=b'x'*65536; "
                        "exec('while True:\\n os.write(descriptor, chunk)')"
                    ).format(str(marker), descriptor)
                    registry = campaign.common.ProcessRegistry()
                    started = time.monotonic()
                    try:
                        with mock.patch.object(
                                campaign.common, "JOB_OUTPUT_MAX_BYTES", 4096), \
                             self.assertRaisesRegex(
                                 campaign.CampaignError,
                                 "bounded capture limit"):
                            campaign.subprocess_executor(
                                (sys.executable, "-c", writer), 5.0,
                                threading.Event(), registry, 0)
                    finally:
                        if marker.exists():
                            process_group = int(marker.read_text(
                                encoding="ascii"))
                            try:
                                os.killpg(process_group, 0)
                            except ProcessLookupError:
                                pass
                            else:
                                os.killpg(process_group, signal.SIGKILL)
                                self.fail("overflowing executor group survived")
                    self.assertLess(time.monotonic() - started, 2.0)
                    self.assertTrue(marker.exists())
                    self.assertEqual(registry.count(), 0)

    def test_subprocess_executor_abort_and_timeout_prove_group_cleanup(
            self) -> None:
        real_popen = campaign.subprocess.Popen
        for reason in ("campaign abort", "timeout"):
            with self.subTest(reason=reason):
                launched = []

                def launch(*args, **kwargs):
                    process = real_popen(*args, **kwargs)
                    launched.append(process)
                    return process

                abort = threading.Event()
                if reason == "campaign abort":
                    abort.set()
                timeout = 5.0 if abort.is_set() else 0.03
                started = time.monotonic()
                with mock.patch.object(
                        campaign.subprocess, "Popen", side_effect=launch), \
                     self.assertRaisesRegex(
                         campaign.CampaignError, reason):
                    campaign.subprocess_executor(
                        ("/bin/sh", "-c", "trap '' TERM; sleep 20 & wait"),
                        timeout, abort, campaign.common.ProcessRegistry(), 0)
                self.assertLess(time.monotonic() - started, 1.0)
                self.assertEqual(len(launched), 1)
                self.assertFalse(
                    campaign.common.process_group_exists(launched[0]))

    def test_end_history_rescan_rejects_streak_split_across_guard_mark(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            identity = (1, 2, 0, 0)

            def execute(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"", 100, 200, cpu)

            def fake_guard(
                path: Path, abort: threading.Event, **_kwargs: object,
            ) -> FakeThermalGuard:
                return FakeThermalGuard(path, abort)

            runner = campaign.JobRunner(
                root, binary, (3,), 1, 10.0, execute)
            with trust_patch(contract), \
                 mock.patch.object(
                     campaign, "validate_frozen_thermal_source",
                     side_effect=(
                         identity, identity,
                         campaign.CampaignError("split high streak"),
                     ),
                 ) as validate_history, \
                 mock.patch.object(
                     campaign.common, "ThermalGuard",
                     side_effect=fake_guard), \
                 self.assertRaisesRegex(
                     campaign.CampaignError, "split high streak"):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertEqual(validate_history.call_count, 3)
            self.assertFalse(
                campaign.phase_completion_path(root, ledger).exists())

    def test_transient_executable_mutation_is_caught_before_second_launch(
            self) -> None:
        for mutation_target in ("binary", "taskset"):
            with self.subTest(target=mutation_target), \
                 tempfile.TemporaryDirectory() as temporary:
                root, binary = make_result_root(temporary)
                ledger = campaign.build_control_ledger(
                    {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                    schedules=("burst", "adversarial"), losses=("0.50",),
                    widths=(64,))
                thermal = root / "thermal-live.csv"
                contract = runner_contract((3,), 1, 10.0, thermal)
                taskset = Path(contract["system_tools"]["taskset"]["path"])
                original_binary = binary.read_bytes()
                original_taskset = taskset.read_bytes()
                calls = 0

                def execute(
                    _command: Tuple[str, ...], _timeout: float,
                    _abort: threading.Event, _registry: object, cpu: int,
                ) -> campaign.ExecutionResult:
                    nonlocal calls
                    job = ledger.jobs[calls]
                    calls += 1
                    if calls == 1:
                        target = binary if mutation_target == "binary" else taskset
                        target.write_bytes(b"transient executable mutation\n")
                    return campaign.ExecutionResult(
                        0, control_output(job), b"", 100, 200, cpu)

                class RestoringThermal(FakeThermalGuard):
                    def finish(self, output: Path) -> dict:
                        binary.write_bytes(original_binary)
                        binary.chmod(0o755)
                        taskset.write_bytes(original_taskset)
                        taskset.chmod(0o755)
                        return super().finish(output)

                runner = campaign.JobRunner(
                    root, binary, (3,), 1, 10.0, execute,
                    RestoringThermal)
                with trust_patch(contract), self.assertRaisesRegex(
                        campaign.CampaignError, "changed"):
                    runner.run_ledger(
                        ledger, thermal, install_signal_handlers=False)
                self.assertEqual(calls, 1)
                self.assertFalse(
                    campaign.phase_completion_path(root, ledger).exists())

    def test_binary_must_remain_executable_and_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)

            def mutate_binary(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                binary.write_bytes(b"mutated benchmark binary\n")
                binary.chmod(0o755)
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"", 100, 200, cpu)

            runner = campaign.JobRunner(
                root, binary, (3,), 1, 10.0, mutate_binary,
                FakeThermalGuard)
            binary.chmod(0o644)
            with trust_patch(contract), self.assertRaisesRegex(
                    campaign.CampaignError, "nonexecutable"):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertFalse(
                (root / "ledgers/control-control.json").exists())

            binary.write_bytes(b"mock frozen binary\n")
            binary.chmod(0o755)
            with trust_patch(contract), self.assertRaisesRegex(
                    campaign.CampaignError, "binary changed"):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertFalse(
                (root / "completed/control-control.json").exists())

    def test_taskset_atomic_receipt_completion_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            calls = []

            def execute(
                command: Tuple[str, ...], timeout: float,
                abort: threading.Event, registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                calls.append(command)
                self.assertEqual(
                    command[:3],
                    (contract["system_tools"]["taskset"]["path"], "-c", "3"))
                self.assertFalse(abort.is_set())
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"", 100, 200, cpu)

            runner = campaign.JobRunner(
                root, binary, (3,), 1, 10.0, execute, FakeThermalGuard)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((3,), 1, 10.0, thermal)
            taskset_path = Path(
                contract["system_tools"]["taskset"]["path"])
            taskset_alias = taskset_path.with_name("taskset-hardlink")
            os.link(taskset_path, taskset_alias)
            self.assertEqual(taskset_path.stat().st_nlink, 2)
            taskset_loop = taskset_path.with_name("taskset-loop")
            taskset_loop.symlink_to(taskset_loop.name)
            loop_contract = copy.deepcopy(contract)
            loop_contract["system_tools"]["taskset"]["path"] = \
                str(taskset_loop)
            with self.assertRaises(campaign.CampaignError):
                campaign.frozen_taskset_path(loop_contract)
            guards_before = FakeThermalGuard.instances
            with trust_patch(contract) as verifier:
                first = runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)
                second = runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)
            self.assertTrue(taskset_alias.samefile(taskset_path))
            self.assertEqual(taskset_path.stat().st_nlink, 2)
            self.assertEqual(len(calls), 1)
            self.assertEqual(FakeThermalGuard.instances - guards_before, 1)
            self.assertEqual(first[0].receipt, second[0].receipt)
            self.assertGreaterEqual(verifier.call_count, 2)
            completion_path = root / "completed" / "control-control.json"
            completion_sidecar = completion_path.with_suffix(".json.sha256")
            completion_sidecar_bytes = completion_sidecar.read_bytes()
            completion_sidecar.unlink()
            with trust_patch(contract):
                repaired = runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertEqual(repaired[0].receipt, first[0].receipt)
            self.assertEqual(len(calls), 1)
            self.assertEqual(FakeThermalGuard.instances - guards_before, 1)
            self.assertEqual(
                completion_sidecar.read_bytes(), completion_sidecar_bytes)
            completion_sidecar.write_bytes(b"wrong completion hash\n")
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            completion_sidecar.write_bytes(completion_sidecar_bytes)
            completion = campaign.load_canonical_object(
                completion_path, "test completion")
            campaign.verify_sealed_record(
                completion, campaign.SCHEMA_PREFIX + ".phase_completion.v1")
            self.assertEqual(completion["logical_rows"], 1)
            self.assertEqual(completion["physical_rows"], 1)
            interval_path = root / completion["thermal_artifact"]["path"]
            interval_bytes = interval_path.read_bytes()
            self.assertEqual(
                campaign.sha256_bytes(interval_bytes),
                completion["thermal_artifact"]["sha256"])
            self.assertEqual(completion["semantic_rows_sha256"],
                             completion["semantic_rows_sha256"].lower())

            completion_sidecar.unlink()
            interval_path.write_bytes(interval_bytes + b"tamper before repair\n")
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertFalse(completion_sidecar.exists())
            interval_path.write_bytes(interval_bytes)
            with trust_patch(contract):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertEqual(
                completion_sidecar.read_bytes(), completion_sidecar_bytes)

            interval_path.write_bytes(interval_bytes + b"tamper\n")
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)
            interval_path.write_bytes(interval_bytes)

            taskset_bytes = taskset_path.read_bytes()
            taskset_path.write_bytes(taskset_bytes + b"tamper\n")
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)
            taskset_path.write_bytes(taskset_bytes)
            taskset_path.chmod(0o755)

            paths = campaign.job_paths(root, ledger.jobs[0])
            paths["stdout"].write_bytes(paths["stdout"].read_bytes() + b"# tamper\n")
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)

    def test_unreceipted_output_pair_is_discarded_and_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            paths = campaign.job_paths(root, ledger.jobs[0])
            paths["stdout"].parent.mkdir(parents=True)
            paths["stderr"].parent.mkdir(parents=True)
            paths["stdout"].write_bytes(control_output(ledger.jobs[0]))
            paths["stderr"].write_bytes(b"")

            calls = []

            def execute(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                calls.append(cpu)
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"", 300, 400, cpu)

            runner = campaign.JobRunner(
                root, binary, (2,), 1, 10.0, execute,
                FakeThermalGuard)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((2,), 1, 10.0, thermal)
            with trust_patch(contract):
                completed = runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False)
            self.assertEqual(calls, [2])
            self.assertFalse(completed[0].receipt["recovered_after_interrupt"])
            self.assertEqual(completed[0].receipt["cpu"], 2)

    def test_unreceipted_output_cleanup_recovers_only_atomic_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stdout = root / "job.stdout"
            stderr = root / "job.stderr"
            marker = root / ".job.stdout.99999999.partial"
            stdout.write_bytes(b"stdout\n")
            stderr.write_bytes(b"stderr\n")
            os.link(stdout, marker)

            campaign.discard_partial_output_pair(
                {"stdout": stdout, "stderr": stderr}, "job")
            self.assertFalse(stdout.exists())
            self.assertFalse(stderr.exists())
            self.assertFalse(marker.exists())

            malicious = root / "malicious.stdout"
            malicious_alias = root / "malicious-alias"
            malicious.write_bytes(b"untrusted\n")
            os.link(malicious, malicious_alias)
            with self.assertRaisesRegex(campaign.CampaignError, "nonunique"):
                campaign.discard_partial_output_pair(
                    {"stdout": malicious, "stderr": root / "missing"},
                    "malicious")
            self.assertTrue(malicious.exists())
            self.assertTrue(malicious_alias.exists())

    def test_hashed_artifact_resume_distinguishes_sidecar_partial(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "record.json"
            encoded = b"{}\n"
            artifact.write_bytes(encoded)
            sidecar = artifact.with_suffix(".json.sha256")
            stale_sidecar = root / ".record.json.sha256.99999999.partial"
            stale_sidecar.write_bytes(b"interrupted sidecar\n")

            digest = campaign.write_hashed_artifact(artifact, encoded)

            self.assertEqual(digest, hashlib.sha256(encoded).hexdigest())
            self.assertEqual(
                sidecar.read_bytes(),
                (digest + "  record.json\n").encode("ascii"))
            self.assertFalse(stale_sidecar.exists())

    def test_failed_process_never_publishes_resumable_success_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))

            def failed(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                return campaign.ExecutionResult(
                    2, control_output(ledger.jobs[0]), b"", 100, 200, cpu)

            runner = campaign.JobRunner(
                root, binary, (1,), 1, 10.0, failed, LowBusyThermalGuard)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((1,), 1, 10.0, thermal)
            with trust_patch(contract), self.assertRaisesRegex(
                    campaign.CampaignError, "failed with rc=2"):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            paths = campaign.job_paths(root, ledger.jobs[0])
            self.assertFalse(paths["stdout"].exists())
            self.assertFalse(paths["stderr"].exists())
            self.assertFalse(paths["receipt"].exists())

    def test_phase_lock_prevents_concurrent_cleanup_or_execution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            calls = []

            def execute(*_args: object) -> campaign.ExecutionResult:
                calls.append(True)
                raise AssertionError("concurrent phase must not execute")

            runner = campaign.JobRunner(
                root, binary, (1,), 1, 10.0, execute, FakeThermalGuard)
            lock_path = root / "locks" / "control-control.lock"
            with trust_patch({}), campaign.job_lock(lock_path), \
                    self.assertRaisesRegex(
                        campaign.CampaignError, "another worker owns"):
                runner.run_ledger(
                    ledger, root / "unused-live-thermal.csv",
                    install_signal_handlers=False)
            self.assertEqual(calls, [])

    def test_unsealed_low_busy_phase_restarts_all_receipted_jobs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), (SEEDS[0],), PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            calls = []

            def execute(
                _command: Tuple[str, ...], _timeout: float,
                _abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                calls.append(cpu)
                tick = len(calls) * 100
                return campaign.ExecutionResult(
                    0, control_output(ledger.jobs[0]), b"",
                    tick, tick + 10, cpu)

            thermal = root / "thermal-live.csv"
            contract = runner_contract((1,), 1, 10.0, thermal)
            low_runner = campaign.JobRunner(
                root, binary, (1,), 1, 10.0, execute, LowBusyThermalGuard)
            with trust_patch(contract), self.assertRaisesRegex(
                    campaign.CampaignError, "CPU utilization"):
                low_runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertTrue(
                campaign.job_paths(root, ledger.jobs[0])["receipt"].exists())
            good_runner = campaign.JobRunner(
                root, binary, (1,), 1, 10.0, execute, FakeThermalGuard)
            with trust_patch(contract):
                completed = good_runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertEqual(calls, [1, 1])
            self.assertEqual(completed[0].receipt["start_ns"], 200)

    def test_first_failed_future_aborts_other_workers_before_shutdown_wait(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, binary = make_result_root(temporary)
            ledger = campaign.build_control_ledger(
                {0: (10,)}, (10,), SEEDS, PARENT,
                schedules=("burst",), losses=("0.50",), widths=(64,))
            second_started = threading.Event()
            second_saw_abort = []

            def execute(
                command: Tuple[str, ...], _timeout: float,
                abort: threading.Event, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                seed = command[command.index("--seed") + 1]
                if seed == SEEDS[0]:
                    if not second_started.wait(2.0):
                        raise AssertionError("second worker did not start")
                    return campaign.ExecutionResult(2, b"", b"", 100, 200, cpu)
                second_started.set()
                second_saw_abort.append(abort.wait(2.0))
                return campaign.ExecutionResult(2, b"", b"", 100, 200, cpu)

            runner = campaign.JobRunner(
                root, binary, (1, 2), 2, 10.0, execute, FakeThermalGuard)
            thermal = root / "thermal-live.csv"
            contract = runner_contract((1, 2), 2, 10.0, thermal)
            with trust_patch(contract), self.assertRaises(campaign.CampaignError):
                runner.run_ledger(
                    ledger, thermal, install_signal_handlers=False)
            self.assertEqual(second_saw_abort, [True])
            self.assertFalse(
                campaign.phase_completion_path(root, ledger).exists())


class AnalysisTests(unittest.TestCase):
    ROUNDS = (
        campaign.RoundSpec("r1", 4, 2, (0,), (64,)),
        campaign.RoundSpec("r2", 2, 1, (1,), (64, 1280), True),
    )

    @staticmethod
    def install_route(
        evidence: campaign.DevelopmentEvidence,
        round_name: str,
        record: campaign.RouteRecord,
    ) -> None:
        key = (record.K, record.width, record.preferred_attempt)
        old = evidence.routes.get(key)
        if old is not None:
            if campaign.route_semantics(old) != campaign.route_semantics(record):
                raise AssertionError("fixture route semantics changed")
        else:
            evidence.routes[key] = record
        evidence.round_routes[(round_name, *key)] = record

    def test_cumulative_masked_ranking_hard_filter_and_no_backfill(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, _binary = make_result_root(temporary)
            evidence = campaign.DevelopmentEvidence(("burst",), ("0.50",))
            K = 10
            root0 = campaign.CellKey(K, 64, 0, "burst", "0.50")
            root1_64 = campaign.CellKey(K, 64, 1, "burst", "0.50")
            root1_1280 = campaign.CellKey(K, 1280, 1, "burst", "0.50")
            evidence.controls[root0] = metrics(result=1)
            # These future controls and candidates are physically cached before
            # R1.  R1 must nevertheless rank only its root-0 panel.
            evidence.controls[root1_64] = metrics(result=0)
            evidence.controls[root1_1280] = metrics(result=1)

            for p in (0, 1, 2, 3):
                self.install_route(evidence, "r1", route_record(K, 64, p))
            for p in (2, 3):
                self.install_route(evidence, "r2", route_record(K, 64, p))
            self.install_route(
                evidence, "r2", route_record(K, 1280, 2, a0=1, valid=0))
            self.install_route(
                evidence, "r2", route_record(K, 1280, 3, a0=1, valid=1))

            evidence.candidates[(2, root0)] = metrics(result=0)
            evidence.candidates[(3, root0)] = metrics(result=0)
            evidence.candidates[(2, root1_64)] = metrics(result=1)
            evidence.candidates[(3, root1_64)] = metrics(result=0)
            evidence.candidates[(3, root1_1280)] = metrics(result=0)

            incoming1 = {K: (0, 1, 2, 3)}
            with trust_patch():
                first = campaign.analyze_round(
                    root, evidence, self.ROUNDS, self.ROUNDS[0], (K,),
                    incoming1, PARENT)
            first_survivors = campaign.survivors_from_record(
                first, "r1", (K,), 2)
            # p2 and p3 tie on R1, so p is the final tie breaker.  If future
            # root 1 leaked into R1, p2's future introduction would reverse it.
            self.assertEqual(first_survivors[K], (2, 3))
            self.assertEqual(first["current_round_logical_cells"], 4)
            self.assertEqual(first["current_round_physical_candidate_solves"], 2)
            self.assertEqual(first["current_round_alias_cells"], 2)

            derived = campaign.derive_round_input(
                self.ROUNDS, self.ROUNDS[1], (K,), PARENT, first)
            self.assertEqual(derived.survivors, first_survivors)
            self.assertEqual(
                derived.parent_sha256,
                campaign.sha256_bytes(campaign.canonical_json_bytes(first)))

            first_hash = campaign.sha256_bytes(campaign.canonical_json_bytes(first))
            with trust_patch():
                second = campaign.analyze_round(
                    root, evidence, self.ROUNDS, self.ROUNDS[1], (K,),
                    first_survivors, first_hash)
            second_survivors = campaign.survivors_from_record(
                second, "r2", (K,), 1)
            self.assertEqual(second_survivors[K], (3,))
            self.assertEqual(second["current_round_logical_cells"], 4)
            self.assertEqual(second["current_round_physical_candidate_solves"], 3)
            self.assertEqual(second["current_round_alias_cells"], 1)
            row = second["K"][0]
            self.assertEqual(row["input"], [2, 3])
            self.assertNotIn(0, [item["p"] for item in row["candidates"]])
            p2 = next(item for item in row["candidates"] if item["p"] == 2)
            p3 = next(item for item in row["candidates"] if item["p"] == 3)
            self.assertIn("fallback", p2["hard_reasons"])
            self.assertTrue(p3["strictly_better_than_control"])
            self.assertEqual(p3["objective"][5], -2)

            forged = copy.deepcopy(second)
            del forged["schema"]
            del forged["self_sha256_excluding_field"]
            forged["K"][0]["retained"] = [2]
            forged_record = campaign.sealed_record(
                campaign.SCHEMA_PREFIX + ".survivors.v1", forged)
            with self.assertRaises(campaign.CampaignError):
                campaign.survivors_from_record(forged_record, "r2", (K,), 1)

    def test_valid_wide_noop_is_neutral_alias_not_hard_rejection(self) -> None:
        evidence = campaign.DevelopmentEvidence(("burst",), ("0.50",))
        K = 10
        p = 2
        panel64 = campaign.Panel(0, 64)
        panel1280 = campaign.Panel(0, 1280)
        cell64 = campaign.CellKey(K, 64, 0, "burst", "0.50")
        cell1280 = campaign.CellKey(K, 1280, 0, "burst", "0.50")
        evidence.controls[cell64] = metrics(result=0)
        evidence.controls[cell1280] = metrics(result=0)
        evidence.routes[(K, 64, p)] = route_record(K, 64, p, a0=1)
        evidence.routes[(K, 1280, p)] = route_record(K, 1280, p, a0=2)
        evidence.candidates[(p, cell64)] = metrics(result=0)

        score = campaign.score_candidate(
            evidence, K, p, (panel64, panel1280), {panel64, panel1280})
        self.assertTrue(score.eligible)
        self.assertEqual(score.reasons, ())
        self.assertIsNotNone(score.objective)

    def test_inactivation_rank_sums_failure_cells_not_only_common_success(self) -> None:
        evidence = campaign.DevelopmentEvidence(("burst",), ("0.50",))
        K = 10
        p = 2
        panel = campaign.Panel(0, 64)
        cell = campaign.CellKey(K, 64, 0, "burst", "0.50")
        evidence.controls[cell] = metrics(result=1, inactivated=100)
        evidence.routes[(K, 64, p)] = route_record(K, 64, p)
        evidence.candidates[(p, cell)] = metrics(result=1, inactivated=7)
        score = campaign.score_candidate(evidence, K, p, (panel,), {panel})
        self.assertTrue(score.eligible)
        self.assertEqual(score.objective[11], 7)
        self.assertEqual(
            campaign.control_objective(evidence, K, (panel,), {panel})[11],
            100)

    def test_final_r4_global_table_objective_exact_join_and_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, _binary = make_result_root(temporary)
            evidence = campaign.DevelopmentEvidence(("burst",), ("0.50",))
            cohort = (10, 11)
            panels = campaign.declared_panels_through(
                campaign.DEFAULT_ROUNDS, "r4")
            for K in cohort:
                for index, panel in enumerate(panels):
                    cell = campaign.CellKey(
                        K, panel.width, panel.root_index, "burst", "0.50")
                    failed = int(index < 2)
                    if K == 10:
                        evidence.controls[cell] = metrics(result=failed)
                    else:
                        evidence.controls[cell] = metrics(
                            result=failed, heavy_shortfall=failed,
                            inactivated=20, binary_deficit=17)
            for width in campaign.APPROVED_WIDTHS:
                route = route_record(10, width, 2)
                evidence.routes[(10, width, 2)] = route
                evidence.round_routes[("r4", 10, width, 2)] = route
            for panel in panels:
                cell = campaign.CellKey(
                    10, panel.width, panel.root_index, "burst", "0.50")
                evidence.candidates[(2, cell)] = metrics(
                    inactivated=20, binary_deficit=16)

            incoming = {10: (2,), 11: ()}
            with trust_patch():
                final = campaign.analyze_round(
                    root, evidence, campaign.DEFAULT_ROUNDS,
                    campaign.DEFAULT_ROUNDS[-1], cohort, incoming, PARENT)
                report = campaign.global_table_objective_report(
                    root, evidence, campaign.DEFAULT_ROUNDS, cohort, final)
                digest = campaign.write_global_table_objective_report(
                    root, evidence, campaign.DEFAULT_ROUNDS, cohort, final)
            self.assertEqual(
                tuple(report), campaign.GLOBAL_TABLE_OBJECTIVE_FIELDS)
            self.assertEqual(report, {
                "errors": 0,
                "max_K_failure_multiplicity": 2,
                "multi_failure_K_count": 1,
                "total_failures": 2,
                "introductions": 0,
                "worst_adverse_axis_delta": 0,
                "sum_positive_axis_deltas": 0,
                "heavy_shortfall_sum": 2,
                "binary_deficit_gt15_count": 38,
                "binary_deficit_max": 17,
                "binary_deficit_sum": 627,
                "common_success_xors": 3400,
                "common_success_muladds": 340,
            })
            artifact = root / "analysis" / "global_table_objective_report.json"
            self.assertEqual(digest, campaign.sha256_bytes(artifact.read_bytes()))
            loaded_report = campaign.load_canonical_object(
                artifact, "global report")
            self.assertEqual(
                campaign.validate_global_table_objective_report(loaded_report),
                report)

            forged_payload = copy.deepcopy(final)
            del forged_payload["schema"]
            del forged_payload["self_sha256_excluding_field"]
            forged_payload["current_round_logical_cells"] += 1
            forged_payload["current_round_alias_cells"] += 1
            forged = campaign.sealed_record(
                campaign.SCHEMA_PREFIX + ".survivors.v1", forged_payload)
            with trust_patch(), self.assertRaises(campaign.CampaignError):
                campaign.global_table_objective_report(
                    root, evidence, campaign.DEFAULT_ROUNDS, cohort, forged)


if __name__ == "__main__":
    unittest.main()
