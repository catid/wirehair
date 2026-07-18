#!/usr/bin/env python3
"""Deterministic boundary fixtures for sealed preferred-attempt holdouts."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest import mock
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import wh2_preferred_attempt_campaign as campaign
import wh2_preferred_attempt_holdout as subject


CONTEXT = "12" * 32
ROOT_HASH = "34" * 32
PARENT_HASH = "56" * 32


def thermal_row(monotonic_s: float) -> bytes:
    fields = [
        "2026-07-18T00:00:00.000Z", "{:.6f}".format(monotonic_s),
        "100.0", "3600.0", "60.0",
        *("{:.2f}".format(50.0 + index) for index in range(8)),
        "0", "128.0", "128.0", "128.0", "0", "0",
    ]
    if len(fields) != len(campaign.common.THERMAL_FIELDS):
        raise AssertionError("thermal fixture field count changed")
    return (",".join(fields) + "\n").encode("ascii")


PHASE_THERMAL_BASELINE_ROW = thermal_row(100.0)
PHASE_THERMAL_HEADER = (
    ",".join(campaign.common.THERMAL_FIELDS) + "\n").encode("ascii")
PHASE_THERMAL_ARTIFACT_ONE = (
    PHASE_THERMAL_HEADER + PHASE_THERMAL_BASELINE_ROW + thermal_row(101.0))
PHASE_THERMAL_ARTIFACT_TWO = (
    PHASE_THERMAL_ARTIFACT_ONE + thermal_row(102.0))
PHASE_THERMAL_BASELINE_ROW_SHA256 = hashlib.sha256(
    PHASE_THERMAL_BASELINE_ROW).hexdigest()


def metrics(
    result: int = 0,
    heavy_shortfall: int = 0,
    inactivated: int = 5,
    binary_deficit: int = 0,
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
        heavy_gain=0,
        block_xors=block_xors,
        block_muladds=block_muladds,
    )


def preferred_route(
    K: int, width: int, p: int = 2, a0: int = 1,
    canonical_probes: Optional[int] = None,
) -> campaign.RouteRecord:
    valid = int(p >= a0)
    no_op = int(p == a0)
    fallback = int(not valid)
    direct = int(bool(valid) and not bool(no_op))
    return campaign.RouteRecord(
        K, width, "preferred", p, a0, p if direct else a0,
        valid, fallback, no_op, direct,
        a0 + 1 if canonical_probes is None else canonical_probes,
        int(p > a0),
    )


def control_route(K: int, width: int, a0: int = 1) -> campaign.RouteRecord:
    return campaign.RouteRecord(
        K, width, "control", -1, a0, a0, 1, 0, 1, 0, a0 + 1, 0)


def route_bytes(records: Sequence[campaign.RouteRecord]) -> bytes:
    lines = [
        campaign.expected_route_preamble(CONTEXT),
        ",".join(campaign.ROUTE_HEADER),
    ]
    for value in records:
        lines.append(",".join(str(item) for item in (
            value.K, value.width, value.route_status,
            value.preferred_attempt, value.canonical_attempt,
            value.actual_attempt, value.preferred_valid, value.fallback,
            value.no_op, value.direct, value.canonical_probe_solves,
            value.preferred_probe_solves,
        )))
    return ("\n".join(lines) + "\n").encode("ascii")


def metric_fields(value: campaign.RecoveryMetrics) -> List[str]:
    return [str(item) for item in (
        value.result, value.rank_fail, value.error, value.heavy_shortfall,
        value.inactivated, value.binary_deficit, value.heavy_gain,
        value.block_xors, value.block_muladds,
    )]


def recovery_row(
    K: int, width: int, arm: str, p: Optional[int],
    value: campaign.RecoveryMetrics, a0: int = 1,
) -> str:
    if arm == "control":
        prefix = (K, width, "control", -1, a0, a0, 0, 1, 0, 0, 0, 1)
    elif p is None:
        prefix = (K, width, "candidate", -1, a0, a0, 0, 1, 0, 1, 0, 0)
    else:
        prefix = (K, width, "candidate", p, a0, p, 1, 1, 0, 0, 1, 1)
    return ",".join([str(item) for item in prefix] + metric_fields(value))


def paired_bytes(
    job: subject.HoldoutJob,
    control: Optional[campaign.RecoveryMetrics] = None,
    candidate: Optional[campaign.RecoveryMetrics] = None,
) -> bytes:
    control = metrics() if control is None else control
    candidate = control if candidate is None else candidate
    lines = [subject._paired_preamble(job), ",".join(campaign.RECOVERY_HEADER)]
    mapping = job.routed_map()
    for K in job.Ks:
        for width in job.widths:
            lines.append(recovery_row(K, width, "control", None, control))
            lines.append(recovery_row(
                K, width, "candidate", mapping.get(K), candidate))
    return ("\n".join(lines) + "\n").encode("ascii")


def cache_binding(data: bytes, path: str = "cache/routes.csv") \
        -> subject.RouteCacheBinding:
    return subject.RouteCacheBinding(path, hashlib.sha256(data).hexdigest())


def paired_job(
    binding: subject.RouteCacheBinding, routed: bool = True,
) -> subject.HoldoutJob:
    return subject.HoldoutJob(
        "h1-paired-00000", "h1", "paired", 0, (4096,),
        ((4096, 2),) if routed else (), subject.WIDTHS, CONTEXT,
        0, "0x0000000000000001", subject.SCHEDULES[0], subject.LOSSES[0],
        binding,
    )


def bins_for(cohort: Sequence[int]) -> Dict[int, Tuple[int, ...]]:
    return {
        index: tuple(cohort[index::subject.BIN_COUNT])
        for index in range(subject.BIN_COUNT)}


def roots(count: int) -> Tuple[str, ...]:
    return tuple("0x{:016x}".format(index + 1) for index in range(count))


def cell_grid(
    cohort: Sequence[int],
    table: Mapping[int, Optional[int]],
    roots_count: int,
    widths: Sequence[int],
    pairs: Optional[Mapping[Tuple[int, int, int, str, str],
                            Tuple[campaign.RecoveryMetrics,
                                  campaign.RecoveryMetrics]]] = None,
) -> Tuple[subject.PairedCell, ...]:
    pairs = {} if pairs is None else pairs
    result: List[subject.PairedCell] = []
    for K in cohort:
        for width in widths:
            for root in range(roots_count):
                for schedule in subject.SCHEDULES:
                    for loss in subject.LOSSES:
                        key = (K, width, root, schedule, loss)
                        control, candidate = pairs.get(
                            key, (metrics(), metrics()))
                        p = table[K]
                        route = None if p is None else preferred_route(K, width, p)
                        result.append(subject.PairedCell(
                            K, width, root, schedule, loss, p, route,
                            control, candidate))
    return tuple(result)


class HoldoutFixtureTest(unittest.TestCase):
    def test_end_history_rescan_rejects_split_holdout_streak(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            binary = root / "frozen/wirehair_v2_bench"
            binary.parent.mkdir(parents=True)
            binary.write_bytes(b"mock binary\n")
            binary.chmod(0o755)
            manifest = subject.canonical_selected_route_manifest(
                tuple(preferred_route(4096, width)
                      for width in subject.WIDTHS),
                (4096,), subject.WIDTHS, {4096: 2}, CONTEXT)
            cache = root / "cache/routes.csv"
            cache.parent.mkdir(parents=True)
            cache.write_bytes(manifest)
            job = paired_job(cache_binding(manifest))
            ledger = subject.HoldoutLedger(
                "h1", "paired", ROOT_HASH, PARENT_HASH, CONTEXT, (job,))
            thermal = root / "live-thermal.csv"
            thermal.write_bytes(b"mock live thermal\n")
            taskset = root / "taskset"
            taskset.write_bytes(b"mock taskset\n")
            taskset.chmod(0o755)
            policy = {
                "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                "timing_cpu_limit_c": 85.0,
                "timing_dimm_limit_c": 90.0,
                "consecutive_samples": 3, "stale_seconds": 5.0,
                "min_cpu_busy_pct": 95.0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
            }
            contract = {
                "cpu_set": [0], "workers": 1, "timeout_seconds": 10.0,
                "thermal": str(thermal), "thermal_policy": policy,
                "binary_sha256": hashlib.sha256(
                    binary.read_bytes()).hexdigest(),
                "system_tools": {"taskset": {
                    "path": str(taskset),
                    "sha256": hashlib.sha256(
                        taskset.read_bytes()).hexdigest(),
                }},
            }
            summary = {
                "samples": 1, "sealed_samples_including_baseline": 2,
                "cpu_busy_min_pct": 100.0, "cpu_tctl_max_c": 60.0,
                "dimm_max_c": 57.0, "dimm_read_errors_max": 0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
                "thermal_limit_c": 90.0, "thermal_high_samples": 0,
                "thermal_high_max_consecutive_samples": 0,
                "guard_poll_iterations": 1, "guard_samples": 1,
                "guard_high_samples": 0, "guard_limit_c": 90.0,
                "guard_error": None,
                "baseline_row_sha256": PHASE_THERMAL_BASELINE_ROW_SHA256,
            }

            class FakeGuard:
                def __init__(
                    self, _path: Path, _abort: object, **_kwargs: object,
                ) -> None:
                    self.started = False
                    self.error = None

                def start(self) -> None:
                    self.started = True

                def finish(self, path: Path) -> Dict[str, object]:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(PHASE_THERMAL_ARTIFACT_ONE)
                    return dict(summary)

            def executor(
                _command: Tuple[str, ...], _timeout: float,
                _abort: object, _registry: object, cpu: int,
            ) -> campaign.ExecutionResult:
                return campaign.ExecutionResult(
                    0, paired_bytes(job), b"", 10, 20, cpu)

            runner = subject.HoldoutRunner(
                root, binary, (0,), 1, 10.0, executor=executor,
                binding_verifier=lambda _root, _ledger: None,
                runtime_verifier=lambda _root: ({}, contract))
            identity = (1, 2, 0, 0)
            with mock.patch.object(
                    campaign, "validate_frozen_thermal_source",
                    side_effect=(
                        identity, identity,
                        campaign.CampaignError("split holdout streak"),
                    )) as validate_history, \
                 mock.patch.object(
                     campaign.common, "ThermalGuard", FakeGuard), \
                 self.assertRaisesRegex(
                     campaign.CampaignError, "split holdout streak"):
                runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False, exact=False)
            self.assertEqual(validate_history.call_count, 3)
            self.assertFalse(
                (root / "h1" / (subject._phase_stem(ledger) + ".json")).exists())

    def test_transient_executable_mutation_stops_second_holdout_launch(
            self) -> None:
        for mutation_target in ("binary", "taskset"):
            with self.subTest(target=mutation_target), \
                 tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve()
                binary = root / "frozen/wirehair_v2_bench"
                binary.parent.mkdir(parents=True)
                binary.write_bytes(b"mock binary\n")
                binary.chmod(0o755)
                manifest = subject.canonical_selected_route_manifest(
                    tuple(preferred_route(4096, width)
                          for width in subject.WIDTHS),
                    (4096,), subject.WIDTHS, {4096: 2}, CONTEXT)
                cache = root / "cache/routes.csv"
                cache.parent.mkdir(parents=True)
                cache.write_bytes(manifest)
                first = paired_job(cache_binding(manifest))
                second = replace(
                    first, job_id="h1-paired-00001", root_index=1,
                    seed="0x0000000000000002")
                taskset_path = root / "taskset"
                taskset_path.write_bytes(b"mock taskset\n")
                taskset_path.chmod(0o755)
                contract = {
                    "binary_sha256": hashlib.sha256(
                        binary.read_bytes()).hexdigest(),
                    "system_tools": {"taskset": {
                        "path": str(taskset_path),
                        "sha256": hashlib.sha256(
                            taskset_path.read_bytes()).hexdigest(),
                    }},
                }
                calls = 0

                def executor(
                    _command: Tuple[str, ...], _timeout: float,
                    _abort: object, _registry: object, cpu: int,
                ) -> campaign.ExecutionResult:
                    nonlocal calls
                    job = (first, second)[calls]
                    calls += 1
                    if calls == 1:
                        target = (binary if mutation_target == "binary"
                                  else taskset_path)
                        target.write_bytes(b"transient executable mutation\n")
                    return campaign.ExecutionResult(
                        0, paired_bytes(job), b"", 10, 20, cpu)

                runner = subject.HoldoutRunner(
                    root, binary, (0,), 1, 10, executor=executor)
                for job in (first, second):
                    for path in subject.job_paths(root, job).values():
                        path.parent.mkdir(parents=True, exist_ok=True)
                pool = campaign.common.CpuPool((0,), 1)
                abort = threading.Event()
                registry = campaign.common.ProcessRegistry()
                taskset = str(taskset_path)
                runner._run_one(
                    first, pool, abort, registry, taskset, contract)
                with self.assertRaisesRegex(
                        campaign.CampaignError, "changed"):
                    runner._run_one(
                        second, pool, abort, registry, taskset, contract)
                self.assertEqual(calls, 1)

    def test_exact_cardinality_builders_and_malformed_axes(self) -> None:
        h1_cohort = tuple(range(4096, 4096 + subject.H1_K_COUNT))
        h1_table = {K: (2 if K == h1_cohort[0] else None) for K in h1_cohort}
        dummy_caches = {
            index: subject.RouteCacheBinding(
                "cache/h1-{:03d}.csv".format(index), "78" * 32)
            for index in range(subject.BIN_COUNT)}
        h1 = subject.build_h1_ledger(
            bins_for(h1_cohort), h1_cohort, roots(subject.H1_ROOT_COUNT),
            h1_table, dummy_caches, ROOT_HASH, PARENT_HASH, CONTEXT)
        self.assertEqual(len(h1.jobs), 6912)
        self.assertEqual(
            sum(job.expected_logical_rows for job in h1.jobs), 276912)
        self.assertEqual(
            h1.sha256(),
            "2d56d0fd8c0f7a4ff6f4032a56f3694d8d60758c4e86ef3130f6c1f192ec5c75")
        broken = replace(
            h1.jobs[1], schedule=subject.SCHEDULES[1], loss=subject.LOSSES[0])
        with self.assertRaises(subject.CampaignError):
            replace(h1, jobs=(h1.jobs[0], broken, *h1.jobs[2:])).validate()

        all_table = {K: (2 if K == 4096 else None) for K in subject.ALL_K}
        all_bins = bins_for(subject.ALL_K)
        h2_route = subject.build_h2_route_ledger(
            all_bins, roots(subject.H2_ROOT_COUNT), all_table,
            ROOT_HASH, PARENT_HASH, CONTEXT)
        self.assertEqual(len(h2_route.jobs), 128)
        self.assertEqual(
            sum(job.expected_logical_rows for job in h2_route.jobs), 255996)
        self.assertEqual(
            h2_route.sha256(),
            "afc64d6a61a2a3ad6edaf7dda01db3304ad8354c58a5d92f4e941cda1e29c8a7")
        h2_caches = {
            index: subject.RouteCacheBinding(
                "cache/h2-{:03d}.csv".format(index), "9a" * 32)
            for index in range(subject.BIN_COUNT)}
        h2 = subject.build_h2_paired_ledger(
            all_bins, roots(subject.H2_ROOT_COUNT), all_table, h2_caches,
            ROOT_HASH, PARENT_HASH, CONTEXT)
        self.assertEqual(len(h2.jobs), 4608)
        self.assertEqual(
            sum(job.expected_logical_rows for job in h2.jobs), 4607928)
        self.assertEqual(
            h2.sha256(),
            "5cd7a3dc27e529d7499ee884d464048f77e0ee715e56258cef9e4e91a6ab8f68")

    def test_exact_table_hash_and_sign_boundaries(self) -> None:
        digest = subject.table_semantic_sha256((2, 3), {2: None, 3: 7})
        self.assertEqual(
            digest, "cdf6fb846bb7fefa5455cc26d90867600afa1699f91834b281a1f67eb48ffe28")
        self.assertEqual(subject.exact_sign_tail(7, 0), (1, 128, True))
        self.assertEqual(subject.exact_sign_tail(6, 0), (1, 64, False))
        self.assertEqual(subject.exact_sign_tail(0, 0), (1, 1, False))
        for repairs in range(8):
            for introductions in range(8):
                n = repairs + introductions
                expected = sum(
                    math.comb(n, value) for value in range(repairs, n + 1))
                self.assertEqual(
                    subject.exact_sign_tail(repairs, introductions)[0], expected)
        with self.assertRaises(subject.CampaignError):
            subject.exact_sign_tail(True, 0)

    def test_route_derivation_and_strict_paired_parser(self) -> None:
        source = tuple(
            preferred_route(4096, width) for width in subject.WIDTHS)
        selected = {4096: 2}
        exact = subject.canonical_selected_route_manifest(
            source, (4096,), subject.WIDTHS, selected, CONTEXT)
        parsed = subject.parse_route_manifest(
            exact, (4096,), subject.WIDTHS, selected, CONTEXT, True)
        self.assertEqual(len(parsed), 4)
        binding = cache_binding(exact)
        job = paired_job(binding)
        output = paired_bytes(job)
        result = subject.parse_paired_output(job, output, exact)
        self.assertEqual((result.logical_rows, result.physical_rows), (8, 8))
        self.assertEqual(len(result.cells), 4)

        extra = source + (preferred_route(4097, 64),)
        with self.assertRaises(subject.CampaignError):
            subject.canonical_selected_route_manifest(
                extra, (4096,), subject.WIDTHS, selected, CONTEXT)
        repeated = route_bytes((source[0], source[0], *source[1:]))
        with self.assertRaises(subject.CampaignError):
            subject.parse_route_manifest(
                repeated, (4096,), subject.WIDTHS, selected, CONTEXT, False)
        lines = output.decode("ascii").splitlines()
        fields = lines[3].split(",")
        fields[11] = "0"
        lines[3] = ",".join(fields)
        with self.assertRaises(subject.CampaignError):
            subject.parse_paired_output(
                job, ("\n".join(lines) + "\n").encode("ascii"), exact)

    def test_nonrouted_parser_requires_exact_alias(self) -> None:
        exact = subject.canonical_selected_route_manifest(
            tuple(control_route(4096, width) for width in subject.WIDTHS),
            (4096,), subject.WIDTHS, {4096: None}, CONTEXT)
        job = paired_job(cache_binding(exact), routed=False)
        parsed = subject.parse_paired_output(job, paired_bytes(job), exact)
        self.assertEqual(parsed.physical_rows, 4)
        lines = paired_bytes(job).decode("ascii").splitlines()
        candidate = lines[3].split(",")
        candidate[-2] = "101"
        lines[3] = ",".join(candidate)
        with self.assertRaises(subject.CampaignError):
            subject.parse_paired_output(
                job, ("\n".join(lines) + "\n").encode("ascii"), exact)

    def test_h1_keep_neutral_and_prune_boundaries(self) -> None:
        cohort = (4096,)
        table = {4096: 2}
        first = (4096, 64, 0, subject.SCHEDULES[0], subject.LOSSES[0])
        repaired = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {first: (metrics(result=1), metrics())})
        analysis = subject.analyze_h1_cells(
            cohort, table, repaired, strict_domain=False)
        self.assertTrue(analysis["K"][0]["kept"])
        self.assertEqual(
            len(analysis["K"][0]["width_schedule_loss_report"]), 36)

        neutral = tuple(replace(
            cell, candidate=replace(cell.candidate, block_xors=99))
            for cell in cell_grid(cohort, table, 1, subject.WIDTHS))
        neutral_analysis = subject.analyze_h1_cells(
            cohort, table, neutral, strict_domain=False)
        self.assertTrue(neutral_analysis["K"][0]["kept"])

        second = (4096, 64, 0, subject.SCHEDULES[1], subject.LOSSES[0])
        equal = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {
                first: (metrics(result=1), metrics()),
                second: (metrics(), metrics(result=1)),
            })
        equal_analysis = subject.analyze_h1_cells(
            cohort, table, equal, strict_domain=False)
        self.assertFalse(equal_analysis["K"][0]["kept"])
        self.assertIn(
            "equal_nonzero_discordance", equal_analysis["K"][0]["reasons"])

        third = (4096, 256, 0, subject.SCHEDULES[1], subject.LOSSES[1])
        schedule_bad = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {
                first: (metrics(), metrics(result=1)),
                second: (metrics(result=1), metrics()),
                third: (metrics(result=1), metrics()),
            })
        schedule_analysis = subject.analyze_h1_cells(
            cohort, table, schedule_bad, strict_domain=False)
        self.assertFalse(schedule_analysis["K"][0]["kept"])
        self.assertIn("schedule_total", schedule_analysis["K"][0]["reasons"])

        fallback_cells = tuple(
            replace(
                cell,
                route=preferred_route(4096, 256, p=2, a0=3))
            if cell.width == 256 else cell
            for cell in repaired)
        fallback_analysis = subject.analyze_h1_cells(
            cohort, table, fallback_cells, strict_domain=False)
        self.assertFalse(fallback_analysis["K"][0]["kept"])
        self.assertIn("route", fallback_analysis["K"][0]["reasons"])

        no_direct_cells = tuple(
            replace(
                cell, route=preferred_route(4096, 64, p=2, a0=2),
                candidate=cell.control)
            if cell.width == 64 else cell
            for cell in repaired)
        no_direct_analysis = subject.analyze_h1_cells(
            cohort, table, no_direct_cells, strict_domain=False)
        self.assertFalse(no_direct_analysis["K"][0]["kept"])
        self.assertIn("route", no_direct_analysis["K"][0]["reasons"])

        adverse = (4096, 256, 0, subject.SCHEDULES[2], subject.LOSSES[2])
        heavy_bad = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {
                first: (metrics(result=1), metrics()),
                adverse: (
                    metrics(result=1),
                    metrics(result=1, heavy_shortfall=1)),
            })
        heavy_analysis = subject.analyze_h1_cells(
            cohort, table, heavy_bad, strict_domain=False)
        self.assertIn("heavy_shortfall", heavy_analysis["K"][0]["reasons"])

        deficit_bad = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {
                first: (metrics(result=1), metrics()),
                adverse: (
                    metrics(inactivated=20),
                    metrics(inactivated=20, binary_deficit=16)),
            })
        deficit_analysis = subject.analyze_h1_cells(
            cohort, table, deficit_bad, strict_domain=False)
        self.assertIn("deficit_gt15", deficit_analysis["K"][0]["reasons"])
        self.assertIn("max_deficit", deficit_analysis["K"][0]["reasons"])

        totals = subject.ComparisonTotals()
        totals.add(metrics(block_xors=200), metrics(block_xors=201))
        self.assertTrue(totals.work_gate()[0])
        worse = subject.ComparisonTotals()
        worse.add(metrics(block_xors=200), metrics(block_xors=202))
        self.assertFalse(worse.work_gate()[0])

    def test_h1_codec_error_is_hard_and_routes_to_control(self) -> None:
        cohort = (4096,)
        table = {4096: 2}
        key = (4096, 64, 0, subject.SCHEDULES[0], subject.LOSSES[0])
        cells = cell_grid(
            cohort, table, 1, subject.WIDTHS,
            {key: (metrics(result=2), metrics(result=2))})
        analysis = subject.analyze_h1_cells(
            cohort, table, cells, strict_domain=False)
        self.assertFalse(analysis["accepted"])
        self.assertIsNone(analysis["K"][0]["output_preferred_attempt"])
        self.assertIn("codec_error", analysis["K"][0]["reasons"])
        table_record = analysis["frozen_table"]
        payload = {
            name: value for name, value in table_record.items()
            if name not in ("schema", "self_sha256_excluding_field")}
        payload["K_count"] = True
        malformed = campaign.sealed_record(
            subject.SCHEMA_PREFIX + ".h1_table.v1", payload)
        with self.assertRaises(subject.CampaignError):
            subject.parse_h1_table(malformed, strict_domain=False)

    def test_h2_whole_table_exact_sign_and_route_gates(self) -> None:
        cohort = (4096,)
        table = {4096: 2}
        keys = [
            (4096, 64, 0, schedule, loss)
            for schedule in subject.SCHEDULES for loss in subject.LOSSES]
        seven_repairs = {
            key: (metrics(result=1), metrics()) for key in keys[:7]}
        routes = tuple(
            preferred_route(4096, width) for width in subject.WIDTHS)
        accepted = subject.analyze_h2_cells(
            cohort, table, cell_grid(
                cohort, table, 1, (64,), seven_repairs),
            routes, strict_domain=False)
        self.assertTrue(accepted["accepted"])
        self.assertEqual(accepted["table_action"], "accept_unchanged")

        six_repairs = {
            key: (metrics(result=1), metrics()) for key in keys[:6]}
        rejected = subject.analyze_h2_cells(
            cohort, table, cell_grid(
                cohort, table, 1, (64,), six_repairs),
            routes, strict_domain=False)
        self.assertFalse(rejected["accepted"])
        self.assertFalse(rejected["gates"]["exact_sign_p_le_0_01"])
        self.assertEqual(rejected["table_action"], "reject_unchanged")

        duplicated = list(cell_grid(
            cohort, table, 1, (64,), seven_repairs))
        duplicated[-1] = duplicated[0]
        with self.assertRaises(subject.CampaignError):
            subject.analyze_h2_cells(
                cohort, table, duplicated, routes,
                strict_domain=False, streaming=True)

        fallback = list(routes)
        fallback[1] = preferred_route(4096, 256, p=0)
        route_rejected = subject.analyze_h2_cells(
            cohort, table, cell_grid(
                cohort, table, 1, (64,), seven_repairs),
            fallback, strict_domain=False)
        self.assertFalse(route_rejected["gates"]["route_validation"])

    def test_runner_receipt_resume_and_strict_numeric_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            binary = root / "frozen/wirehair_v2_bench"
            binary.parent.mkdir(parents=True)
            binary.write_bytes(b"mock binary\n")
            binary.chmod(0o755)
            exact = subject.canonical_selected_route_manifest(
                tuple(preferred_route(4096, width)
                      for width in subject.WIDTHS),
                (4096,), subject.WIDTHS, {4096: 2}, CONTEXT)
            cache = root / "cache/routes.csv"
            cache.parent.mkdir(parents=True)
            cache.write_bytes(exact)
            binding = cache_binding(exact)
            job = paired_job(binding)
            ledger = subject.HoldoutLedger(
                "h1", "paired", ROOT_HASH, PARENT_HASH, CONTEXT, (job,))
            output = paired_bytes(job)
            calls: List[Tuple[str, ...]] = []
            thermal_starts: List[Path] = []
            live_thermal = root / "live-thermal.csv"
            live_thermal.write_bytes(b"live thermal source\n")
            taskset = root / "taskset"
            taskset.write_bytes(b"mock taskset\n")
            taskset.chmod(0o755)
            policy = {
                "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                "timing_cpu_limit_c": 80.0, "timing_dimm_limit_c": 80.0,
                "consecutive_samples": 3, "stale_seconds": 5.0,
                "min_cpu_busy_pct": 95.0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
            }
            summary = {
                "samples": 2, "sealed_samples_including_baseline": 3,
                "cpu_busy_min_pct": 100.0, "cpu_tctl_max_c": 60.0,
                "dimm_max_c": 57.0, "dimm_read_errors_max": 0,
                "edac_ce_delta": 0, "edac_ue_delta": 0,
                "thermal_limit_c": 90.0, "thermal_high_samples": 0,
                "thermal_high_max_consecutive_samples": 0,
                "guard_poll_iterations": 2, "guard_samples": 2,
                "guard_high_samples": 0, "guard_limit_c": 90.0,
                "guard_error": None,
                "baseline_row_sha256": PHASE_THERMAL_BASELINE_ROW_SHA256,
            }
            contract = {
                "cpu_set": [0], "workers": 1, "timeout_seconds": 10.0,
                "thermal": str(live_thermal), "thermal_policy": policy,
                "binary_sha256": hashlib.sha256(
                    binary.read_bytes()).hexdigest(),
                "system_tools": {"taskset": {
                    "path": str(taskset),
                    "sha256": hashlib.sha256(taskset.read_bytes()).hexdigest(),
                }},
            }

            class FakeThermal:
                def __init__(self, path: Path, _abort: object) -> None:
                    thermal_starts.append(path)
                    self.started = False
                    self.error = None

                def start(self) -> None:
                    self.started = True

                def finish(self, path: Path) -> Dict[str, object]:
                    path.write_bytes(PHASE_THERMAL_ARTIFACT_TWO)
                    return dict(summary)

            def executor(command: Tuple[str, ...], _timeout: float,
                         _abort: object, _registry: object,
                         cpu: int) -> campaign.ExecutionResult:
                calls.append(command)
                return campaign.ExecutionResult(0, output, b"", 10, 20, cpu)

            runner = subject.HoldoutRunner(
                root, binary, (0,), 1, 10, executor=executor,
                binding_verifier=lambda _root, _ledger: None,
                runtime_verifier=lambda _root: ({}, contract),
                thermal_factory=FakeThermal)
            phase_lock = root / "locks/h1-paired.lock"
            with campaign.job_lock(phase_lock), \
                    self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            self.assertEqual(calls, [])
            stale_paths = subject.job_paths(root, job)
            stale_paths["stdout"].parent.mkdir(parents=True, exist_ok=True)
            stale_paths["stderr"].parent.mkdir(parents=True, exist_ok=True)
            stale_paths["stdout"].write_bytes(b"unsealed stale output\n")
            stale_paths["stderr"].write_bytes(b"")
            stale_thermal = subject.phase_thermal_path(root, ledger)
            stale_thermal.parent.mkdir(parents=True, exist_ok=True)
            stale_thermal.write_bytes(b"unsealed stale interval\n")
            summary["cpu_busy_min_pct"] = 1.0
            with self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            self.assertTrue(subject.job_paths(root, job)["receipt"].exists())
            summary["cpu_busy_min_pct"] = 99.0
            runner.run_ledger(
                ledger, live_thermal,
                install_signal_handlers=False, exact=False)
            runner.run_ledger(
                ledger, live_thermal,
                install_signal_handlers=False, exact=False)
            self.assertEqual(len(calls), 2)
            self.assertEqual(calls[0][0], str(taskset))
            self.assertEqual(len(thermal_starts), 2)
            phase_path = root / "h1" / (
                subject._phase_stem(ledger) + ".json")
            phase_sidecar = phase_path.with_suffix(".sha256")
            phase_sidecar_bytes = phase_sidecar.read_bytes()
            phase_sidecar.unlink()
            runner.run_ledger(
                ledger, live_thermal,
                install_signal_handlers=False, exact=False)
            self.assertEqual(len(calls), 2)
            self.assertEqual(len(thermal_starts), 2)
            self.assertEqual(phase_sidecar.read_bytes(), phase_sidecar_bytes)
            phase_sidecar.write_bytes(b"wrong completion hash\n")
            with self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            phase_sidecar.write_bytes(phase_sidecar_bytes)
            original_taskset = taskset.read_bytes()
            taskset.write_bytes(b"tampered taskset\n")
            with self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            taskset.write_bytes(original_taskset)
            interval = subject.phase_thermal_path(root, ledger)
            original_interval = interval.read_bytes()
            phase_sidecar.unlink()
            interval.write_bytes(b"tampered before sidecar repair\n")
            with self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            self.assertFalse(phase_sidecar.exists())
            interval.write_bytes(original_interval)
            runner.run_ledger(
                ledger, live_thermal,
                install_signal_handlers=False, exact=False)
            self.assertEqual(phase_sidecar.read_bytes(), phase_sidecar_bytes)

            interval.write_bytes(b"tampered thermal interval\n")
            with self.assertRaises(subject.CampaignError):
                runner.run_ledger(
                    ledger, live_thermal,
                    install_signal_handlers=False, exact=False)
            interval.write_bytes(original_interval)
            receipt_path = subject.job_paths(root, job)["receipt"]
            receipt = json.loads(receipt_path.read_text("ascii"))
            payload = {
                key: value for key, value in receipt.items()
                if key not in ("schema", "self_sha256_excluding_field")}
            payload["logical_rows"] = True
            malformed = campaign.sealed_record(
                subject.SCHEMA_PREFIX + ".job_receipt.v1", payload)
            receipt_path.write_bytes(campaign.canonical_json_bytes(malformed))
            with self.assertRaises(subject.CampaignError):
                subject.verify_receipt(root, binary, job, {0})

    def test_end_identity_failure_forces_full_holdout_rerun(self) -> None:
        for mutation_target in ("binary", "taskset"):
            with self.subTest(target=mutation_target), \
                 tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve()
                binary = root / "frozen/wirehair_v2_bench"
                binary.parent.mkdir(parents=True)
                binary.write_bytes(b"mock binary\n")
                binary.chmod(0o755)
                manifest = subject.canonical_selected_route_manifest(
                    tuple(preferred_route(4096, width)
                          for width in subject.WIDTHS),
                    (4096,), subject.WIDTHS, {4096: 2}, CONTEXT)
                cache = root / "cache/routes.csv"
                cache.parent.mkdir(parents=True)
                cache.write_bytes(manifest)
                job = paired_job(cache_binding(manifest))
                ledger = subject.HoldoutLedger(
                    "h1", "paired", ROOT_HASH, PARENT_HASH, CONTEXT, (job,))
                output = paired_bytes(job)
                thermal = root / "live-thermal.csv"
                thermal.write_bytes(b"mock live thermal\n")
                taskset = root / "taskset"
                taskset.write_bytes(b"mock taskset\n")
                taskset.chmod(0o755)
                policy = {
                    "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                    "timing_cpu_limit_c": 85.0,
                    "timing_dimm_limit_c": 90.0,
                    "consecutive_samples": 3, "stale_seconds": 5.0,
                    "min_cpu_busy_pct": 95.0,
                    "edac_ce_delta": 0, "edac_ue_delta": 0,
                }
                summary = {
                    "samples": 1, "sealed_samples_including_baseline": 2,
                    "cpu_busy_min_pct": 100.0, "cpu_tctl_max_c": 60.0,
                    "dimm_max_c": 57.0, "dimm_read_errors_max": 0,
                    "edac_ce_delta": 0, "edac_ue_delta": 0,
                    "thermal_limit_c": 90.0, "thermal_high_samples": 0,
                    "thermal_high_max_consecutive_samples": 0,
                    "guard_poll_iterations": 1, "guard_samples": 1,
                    "guard_high_samples": 0, "guard_limit_c": 90.0,
                    "guard_error": None,
                    "baseline_row_sha256": PHASE_THERMAL_BASELINE_ROW_SHA256,
                }
                contract = {
                    "cpu_set": [0], "workers": 1,
                    "timeout_seconds": 10.0, "thermal": str(thermal),
                    "thermal_policy": policy,
                    "binary_sha256": hashlib.sha256(
                        binary.read_bytes()).hexdigest(),
                    "system_tools": {"taskset": {
                        "path": str(taskset),
                        "sha256": hashlib.sha256(
                            taskset.read_bytes()).hexdigest(),
                    }},
                }
                calls = 0
                mutate = True

                class FakeThermal:
                    def __init__(self, _path: Path, _abort: object) -> None:
                        self.started = False
                        self.error = None

                    def start(self) -> None:
                        self.started = True

                    def finish(self, path: Path) -> Dict[str, object]:
                        path.write_bytes(PHASE_THERMAL_ARTIFACT_ONE)
                        return dict(summary)

                def executor(
                    _command: Tuple[str, ...], _timeout: float,
                    _abort: object, _registry: object, cpu: int,
                ) -> campaign.ExecutionResult:
                    nonlocal calls
                    calls += 1
                    if mutate:
                        target = binary if mutation_target == "binary" else taskset
                        target.write_bytes(b"mutated executable\n")
                    return campaign.ExecutionResult(
                        0, output, b"", 10, 20, cpu)

                runner = subject.HoldoutRunner(
                    root, binary, (0,), 1, 10, executor=executor,
                    binding_verifier=lambda _root, _ledger: None,
                    runtime_verifier=lambda _root: ({}, contract),
                    thermal_factory=FakeThermal)
                with self.assertRaises(subject.CampaignError):
                    runner.run_ledger(
                        ledger, thermal,
                        install_signal_handlers=False, exact=False)
                phase = root / "h1" / (subject._phase_stem(ledger) + ".json")
                self.assertFalse(phase.exists())
                binary.write_bytes(b"mock binary\n")
                binary.chmod(0o755)
                taskset.write_bytes(b"mock taskset\n")
                taskset.chmod(0o755)
                mutate = False
                runner.run_ledger(
                    ledger, thermal,
                    install_signal_handlers=False, exact=False)
                self.assertEqual(calls, 2)
                self.assertTrue(phase.exists())


if __name__ == "__main__":
    unittest.main()
