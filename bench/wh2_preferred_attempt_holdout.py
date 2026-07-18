#!/usr/bin/env python3
"""Sealed H1/H2 execution and exact gates for preferred-attempt routing.

The development campaign deliberately stops after producing one optional
preferred attempt per K.  This module consumes that frozen table without ever
looking at a runner-up.  H1 may replace an entry by canonical control; H2 may
only accept or reject the complete H1 table.

The native ``preferredattempt --mode paired`` command is the unit of recovery
work.  Each process emits a matched control/candidate pair for every requested
K/width.  Route manifests are immutable inputs to those processes.  H2 first
re-probes the complete all-K table at all four widths and then derives an exact
bb64-only cache from those validated bytes for the recovery panel.

All production entry points verify the frozen controller/root binding before
touching outcomes.  Tests can exercise the pure builders, parsers, and exact
gate reducers without manufacturing a rooted drand campaign.
"""

from __future__ import annotations

import sys

# Frozen campaign directories have an exact immutable inventory.  Disable
# import-cache writes before loading any local frozen helper module.
sys.dont_write_bytecode = True

from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from contextlib import nullcontext
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import re
import signal
import threading
from typing import (
    Any, Callable, Dict, Iterable, List, Mapping, Optional,
    Sequence, Set, Tuple,
)

import wh2_preferred_attempt_campaign as campaign
import wh2_rank_floor_two_anchor_allk as common


CampaignError = campaign.CampaignError
die = campaign.die

SCHEMA_PREFIX = "wirehair.wh2.h12_preferred_attempt.holdout"
WIDTHS = campaign.APPROVED_WIDTHS
SCHEDULES = campaign.APPROVED_SCHEDULES
LOSSES = campaign.APPROVED_LOSSES
H1_ROOT_COUNT = 6
H2_ROOT_COUNT = 4
H1_K_COUNT = 641
ALL_K = tuple(range(2, 64001))
BIN_COUNT = 128
H1_JOB_COUNT = BIN_COUNT * H1_ROOT_COUNT * len(SCHEDULES) * len(LOSSES)
H2_ROUTE_JOB_COUNT = BIN_COUNT
H2_JOB_COUNT = BIN_COUNT * H2_ROOT_COUNT * len(SCHEDULES) * len(LOSSES)
H1_LOGICAL_ROWS = 2 * H1_K_COUNT * H1_ROOT_COUNT * 9 * len(WIDTHS)
H2_ROUTE_ROWS = len(ALL_K) * len(WIDTHS)
H2_LOGICAL_ROWS = 2 * len(ALL_K) * H2_ROOT_COUNT * 9
U64_MAX = (1 << 64) - 1
U32_MAX = (1 << 32) - 1
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return campaign.canonical_json_bytes(value)


def _sealed(schema: str, payload: Mapping[str, Any]) -> Dict[str, Any]:
    return campaign.sealed_record(schema, payload)


def _bind_analysis(
    record: Mapping[str, Any], schema: str, binding: Mapping[str, Any],
) -> Dict[str, Any]:
    campaign.verify_sealed_record(record, schema)
    payload = {
        key: value for key, value in record.items()
        if key not in ("schema", "self_sha256_excluding_field")}
    if "evidence_binding" in payload:
        die("analysis record was already evidence-bound")
    payload["evidence_binding"] = dict(binding)
    return _sealed(schema, payload)


def _require_hash(value: str, context: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        die("{} is not a lowercase SHA256".format(context))
    return value


def _require_plain_int(
    value: Any, context: str, minimum: int = 0, maximum: int = U64_MAX,
) -> int:
    if (not isinstance(value, int) or isinstance(value, bool) or
            value < minimum or value > maximum):
        die("{} is outside its integer domain".format(context))
    return value


def _safe_relative(value: str, context: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        die("{} is not a safe relative path".format(context))
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != value:
        die("{} is not a safe relative path".format(context))
    return value


def normalize_table(
    cohort: Sequence[int], table: Mapping[int, Optional[int]],
) -> Tuple[Tuple[int, Optional[int]], ...]:
    ordered = tuple(cohort)
    if (ordered != tuple(sorted(ordered)) or len(ordered) != len(set(ordered)) or
            any(not isinstance(K, int) or isinstance(K, bool) or
                K < 2 or K > 64000 for K in ordered) or
            tuple(sorted(table)) != ordered):
        die("preferred table does not exactly cover its canonical cohort")
    result: List[Tuple[int, Optional[int]]] = []
    for K in ordered:
        attempt = table[K]
        if attempt is not None:
            _require_plain_int(attempt, "preferred attempt", 0, 255)
        result.append((K, attempt))
    return tuple(result)


def table_semantic_sha256(
    cohort: Sequence[int], table: Mapping[int, Optional[int]],
) -> str:
    rows = normalize_table(cohort, table)
    return _sha256(_json_bytes([
        {"K": K, "preferred_attempt": p} for K, p in rows]))


@dataclass(frozen=True)
class RouteCacheBinding:
    relative_path: str
    sha256: str

    def validate(self) -> None:
        _safe_relative(self.relative_path, "route-cache path")
        _require_hash(self.sha256, "route-cache hash")

    def record(self) -> Dict[str, str]:
        self.validate()
        return {"path": self.relative_path, "sha256": self.sha256}


@dataclass(frozen=True)
class HoldoutJob:
    job_id: str
    phase: str
    kind: str
    bin_index: int
    Ks: Tuple[int, ...]
    routed: Tuple[Tuple[int, int], ...]
    widths: Tuple[int, ...]
    route_context_sha256: str
    root_index: Optional[int] = None
    seed: Optional[str] = None
    schedule: Optional[str] = None
    loss: Optional[str] = None
    route_cache: Optional[RouteCacheBinding] = None

    def routed_map(self) -> Dict[int, int]:
        return dict(self.routed)

    def validate(self) -> None:
        if self.phase not in ("h1", "h2") or self.kind not in ("route", "paired"):
            die("holdout job phase/kind is invalid")
        if re.fullmatch(
                r"{}-{}-[0-9]{{5}}".format(self.phase, self.kind),
                self.job_id) is None:
            die("holdout job ID is noncanonical")
        _require_plain_int(self.bin_index, "holdout bin", 0, BIN_COUNT - 1)
        if (not self.Ks or self.Ks != tuple(sorted(self.Ks)) or
                len(self.Ks) != len(set(self.Ks)) or
                any(not isinstance(K, int) or isinstance(K, bool) or
                    K < 2 or K > 64000 for K in self.Ks)):
            die("holdout job K list is invalid")
        if self.widths not in (WIDTHS, (64,)):
            die("holdout job width list is not frozen")
        campaign.canonical_context_hash(self.route_context_sha256)
        mapping = self.routed_map()
        if (len(mapping) != len(self.routed) or
                tuple(mapping) != tuple(sorted(mapping)) or
                not set(mapping).issubset(self.Ks) or
                any(not isinstance(p, int) or isinstance(p, bool) or
                    p < 0 or p > 255 for p in mapping.values())):
            die("holdout job preferred map is invalid")
        if self.kind == "route":
            if (self.phase != "h2" or self.widths != WIDTHS or
                    self.root_index is not None or self.seed is not None or
                    self.schedule is not None or self.loss is not None or
                    self.route_cache is not None):
                die("holdout route job has paired-only fields")
            return
        expected_widths = WIDTHS if self.phase == "h1" else (64,)
        expected_roots = H1_ROOT_COUNT if self.phase == "h1" else H2_ROOT_COUNT
        if (self.widths != expected_widths or self.root_index is None or
                not isinstance(self.root_index, int) or
                isinstance(self.root_index, bool) or
                not 0 <= self.root_index < expected_roots or
                self.seed is None or self.schedule not in SCHEDULES or
                self.loss not in LOSSES or self.route_cache is None):
            die("holdout paired job fields are invalid")
        campaign.canonical_seed(self.seed)
        self.route_cache.validate()

    @property
    def expected_logical_rows(self) -> int:
        return (len(self.Ks) * len(self.widths) if self.kind == "route" else
                2 * len(self.Ks) * len(self.widths))

    def record(self) -> Dict[str, Any]:
        self.validate()
        return {
            "job_id": self.job_id, "phase": self.phase, "kind": self.kind,
            "bin": self.bin_index, "K": list(self.Ks),
            "routed": [{"K": K, "p": p} for K, p in self.routed],
            "widths": list(self.widths), "root_index": self.root_index,
            "seed": self.seed, "schedule": self.schedule, "loss": self.loss,
            "route_context_sha256": self.route_context_sha256,
            "route_cache": (
                self.route_cache.record() if self.route_cache is not None else None),
            "expected_logical_rows": self.expected_logical_rows,
        }

    def command(self, binary: Path, result_root: Path) -> Tuple[str, ...]:
        self.validate()
        mapping = self.routed_map()
        if self.kind == "route":
            map_text = "|".join(
                "{}@{}={}".format(
                    K, width, mapping[K] if K in mapping else "none")
                for K in self.Ks for width in self.widths)
        else:
            records = [
                "{}@{}={}".format(K, width, mapping[K])
                for K in self.Ks if K in mapping for width in self.widths
            ]
            map_text = "|".join(records) if records else "none"
        command = [
            str(binary), "preferredattempt", "--mode", self.kind,
            "--N", ",".join(str(K) for K in self.Ks),
            "--bb-list", ",".join(str(width) for width in self.widths),
            "--preferred-map", map_text,
            "--route-context-sha256", self.route_context_sha256,
        ]
        if self.kind == "paired":
            assert self.route_cache is not None
            command.extend((
                "--route-cache", str(result_root / self.route_cache.relative_path),
                "--route-cache-sha256", self.route_cache.sha256,
                "--loss", str(self.loss), "--seed", str(self.seed),
                "--schedule", str(self.schedule),
            ))
        return tuple(command)


@dataclass(frozen=True)
class HoldoutLedger:
    phase: str
    kind: str
    root_file_sha256: str
    parent_table_sha256: str
    route_context_sha256: str
    jobs: Tuple[HoldoutJob, ...]

    def validate(self, exact: bool = True) -> None:
        if self.phase not in ("h1", "h2") or self.kind not in ("route", "paired"):
            die("holdout ledger phase/kind is invalid")
        if self.phase == "h1" and self.kind != "paired":
            die("H1 has no fresh route ledger")
        for value, context in (
                (self.root_file_sha256, "root file"),
                (self.parent_table_sha256, "parent table"),
                (self.route_context_sha256, "route context")):
            _require_hash(value, context + " hash")
        for index, job in enumerate(self.jobs):
            job.validate()
            if (job.phase, job.kind, job.route_context_sha256) != (
                    self.phase, self.kind, self.route_context_sha256) or \
                    job.job_id != "{}-{}-{:05d}".format(
                        self.phase, self.kind, index):
                die("holdout ledger job order/binding changed")
        if not exact:
            return
        expected_jobs = (
            H1_JOB_COUNT if self.phase == "h1" else
            H2_ROUTE_JOB_COUNT if self.kind == "route" else H2_JOB_COUNT)
        expected_rows = (
            H1_LOGICAL_ROWS if self.phase == "h1" else
            H2_ROUTE_ROWS if self.kind == "route" else H2_LOGICAL_ROWS)
        if (len(self.jobs) != expected_jobs or
                sum(job.expected_logical_rows for job in self.jobs) !=
                expected_rows):
            die("holdout ledger exact job/row arithmetic mismatch")
        by_bin: Dict[int, List[HoldoutJob]] = {
            index: [] for index in range(BIN_COUNT)}
        for job in self.jobs:
            by_bin[job.bin_index].append(job)
        if any(not jobs for jobs in by_bin.values()):
            die("holdout ledger has an empty bin")
        representative: Dict[int, HoldoutJob] = {
            index: jobs[0] for index, jobs in by_bin.items()}
        flattened = tuple(
            K for index in range(BIN_COUNT) for K in representative[index].Ks)
        if len(flattened) != len(set(flattened)):
            die("holdout ledger K partition overlaps")
        if self.phase == "h1":
            if len(flattened) != H1_K_COUNT or any(K < 4096 for K in flattened):
                die("H1 ledger does not cover the exact active cardinality")
        elif tuple(sorted(flattened)) != ALL_K:
            die("H2 ledger does not cover exact all-K")
        if self.kind == "route":
            if (any(len(jobs) != 1 for jobs in by_bin.values()) or
                    tuple(job.bin_index for job in self.jobs) !=
                    tuple(range(BIN_COUNT))):
                die("H2 route ledger does not have one job per bin")
            return
        expected_axes = tuple(
            (root, schedule, loss)
            for root in range(
                H1_ROOT_COUNT if self.phase == "h1" else H2_ROOT_COUNT)
            for schedule in SCHEDULES for loss in LOSSES)
        roots: Dict[int, str] = {}
        jobs_per_bin = len(expected_axes)
        if tuple(job.bin_index for job in self.jobs) != tuple(
                index for index in range(BIN_COUNT)
                for _unused in range(jobs_per_bin)):
            die("holdout paired ledger bin-major order changed")
        for bin_index, jobs in by_bin.items():
            first = jobs[0]
            if tuple((job.root_index, job.schedule, job.loss)
                     for job in jobs) != expected_axes:
                die("holdout paired ledger axes/order changed")
            for job in jobs:
                if (job.Ks != first.Ks or job.routed != first.routed or
                        job.route_cache != first.route_cache):
                    die("holdout paired bin bindings changed across axes")
                assert job.root_index is not None and job.seed is not None
                old_seed = roots.setdefault(job.root_index, job.seed)
                if old_seed != job.seed:
                    die("holdout root seed changed across bins/strata")
        if len(roots) != len(expected_axes) // 9 or \
                len(set(roots.values())) != len(roots):
            die("holdout root seed cardinality/uniqueness changed")

    def record(self, exact: bool = True) -> Dict[str, Any]:
        self.validate(exact=exact)
        return _sealed(
            SCHEMA_PREFIX + ".job_ledger.v1",
            {
                "phase": self.phase, "kind": self.kind,
                "root_file_sha256": self.root_file_sha256,
                "parent_table_sha256": self.parent_table_sha256,
                "route_context_sha256": self.route_context_sha256,
                "job_count": len(self.jobs),
                "expected_logical_rows": sum(
                    job.expected_logical_rows for job in self.jobs),
                "jobs": [job.record() for job in self.jobs],
            },
        )

    def sha256(self, exact: bool = True) -> str:
        return _sha256(_json_bytes(self.record(exact=exact)))


def _validate_builder_inputs(
    phase: str,
    bins: Mapping[int, Sequence[int]],
    cohort: Sequence[int],
    roots: Sequence[str],
    table: Mapping[int, Optional[int]],
) -> Tuple[Tuple[str, ...], Tuple[Tuple[int, Optional[int]], ...]]:
    campaign.validate_bins(bins, cohort)
    if len(bins) != BIN_COUNT:
        die("holdout bins do not have the frozen 128-way partition")
    rows = normalize_table(cohort, table)
    expected_roots = H1_ROOT_COUNT if phase == "h1" else H2_ROOT_COUNT
    if len(roots) != expected_roots or len(set(roots)) != len(roots):
        die("holdout root cardinality/uniqueness mismatch")
    root_values = tuple(campaign.canonical_seed(seed) for seed in roots)
    if phase == "h1":
        if len(rows) != H1_K_COUNT or any(K < 4096 for K, _p in rows):
            die("H1 table is not the exact 641-K active cohort")
    elif tuple(K for K, _p in rows) != ALL_K:
        die("H2 table is not the exact all-K domain")
    return root_values, rows


def build_h1_ledger(
    bins: Mapping[int, Sequence[int]],
    cohort: Sequence[int],
    roots: Sequence[str],
    table: Mapping[int, Optional[int]],
    route_caches: Mapping[int, RouteCacheBinding],
    root_file_sha256: str,
    parent_table_sha256: str,
    route_context_sha256: str,
) -> HoldoutLedger:
    roots_tuple, rows = _validate_builder_inputs(
        "h1", bins, cohort, roots, table)
    mapping = dict(rows)
    if set(route_caches) != set(range(BIN_COUNT)):
        die("H1 route-cache bindings do not exactly cover the bins")
    jobs: List[HoldoutJob] = []
    for bin_index in range(BIN_COUNT):
        Ks = tuple(bins[bin_index])
        routed = tuple((K, mapping[K]) for K in Ks if mapping[K] is not None)
        cache = route_caches[bin_index]
        cache.validate()
        for root_index, seed in enumerate(roots_tuple):
            for schedule in SCHEDULES:
                for loss in LOSSES:
                    jobs.append(HoldoutJob(
                        "h1-paired-{:05d}".format(len(jobs)), "h1", "paired",
                        bin_index, Ks, routed, WIDTHS, route_context_sha256,
                        root_index, seed, schedule, loss, cache))
    ledger = HoldoutLedger(
        "h1", "paired", root_file_sha256, parent_table_sha256,
        route_context_sha256, tuple(jobs))
    ledger.validate()
    return ledger


def build_h2_route_ledger(
    bins: Mapping[int, Sequence[int]],
    roots: Sequence[str],
    table: Mapping[int, Optional[int]],
    root_file_sha256: str,
    parent_table_sha256: str,
    route_context_sha256: str,
) -> HoldoutLedger:
    _roots, rows = _validate_builder_inputs(
        "h2", bins, ALL_K, roots, table)
    mapping = dict(rows)
    jobs: List[HoldoutJob] = []
    for bin_index in range(BIN_COUNT):
        Ks = tuple(bins[bin_index])
        routed = tuple((K, mapping[K]) for K in Ks if mapping[K] is not None)
        jobs.append(HoldoutJob(
            "h2-route-{:05d}".format(len(jobs)), "h2", "route",
            bin_index, Ks, routed, WIDTHS, route_context_sha256))
    ledger = HoldoutLedger(
        "h2", "route", root_file_sha256, parent_table_sha256,
        route_context_sha256, tuple(jobs))
    ledger.validate()
    return ledger


def build_h2_paired_ledger(
    bins: Mapping[int, Sequence[int]],
    roots: Sequence[str],
    table: Mapping[int, Optional[int]],
    route_caches: Mapping[int, RouteCacheBinding],
    root_file_sha256: str,
    parent_table_sha256: str,
    route_context_sha256: str,
) -> HoldoutLedger:
    roots_tuple, rows = _validate_builder_inputs(
        "h2", bins, ALL_K, roots, table)
    mapping = dict(rows)
    if set(route_caches) != set(range(BIN_COUNT)):
        die("H2 bb64 route-cache bindings do not exactly cover the bins")
    jobs: List[HoldoutJob] = []
    for bin_index in range(BIN_COUNT):
        Ks = tuple(bins[bin_index])
        routed = tuple((K, mapping[K]) for K in Ks if mapping[K] is not None)
        cache = route_caches[bin_index]
        cache.validate()
        for root_index, seed in enumerate(roots_tuple):
            for schedule in SCHEDULES:
                for loss in LOSSES:
                    jobs.append(HoldoutJob(
                        "h2-paired-{:05d}".format(len(jobs)), "h2", "paired",
                        bin_index, Ks, routed, (64,), route_context_sha256,
                        root_index, seed, schedule, loss, cache))
    ledger = HoldoutLedger(
        "h2", "paired", root_file_sha256, parent_table_sha256,
        route_context_sha256, tuple(jobs))
    ledger.validate()
    return ledger


def _route_preamble(context_sha256: str) -> str:
    return campaign.expected_route_preamble(context_sha256)


def parse_route_manifest(
    data: bytes,
    Ks: Sequence[int],
    widths: Sequence[int],
    table: Mapping[int, Optional[int]],
    route_context_sha256: str,
    exact_one_per_key: bool,
) -> Tuple[campaign.RouteRecord, ...]:
    """Parse a bound manifest and enforce selected/control semantics.

    Development route caches may contain additional preferred attempts for a
    routed key.  Fresh H2 manifests and every derived paired cache contain
    exactly one record per K/width.
    """
    ordered_K = tuple(Ks)
    widths_tuple = tuple(widths)
    if (ordered_K != tuple(sorted(ordered_K)) or len(ordered_K) != len(set(ordered_K)) or
            widths_tuple not in (WIDTHS, (64,))):
        die("route-manifest expected grid is noncanonical")
    selected = dict(normalize_table(ordered_K, table))
    lines = campaign.strict_ascii_lines(data, "holdout route manifest")
    if (len(lines) < 3 or lines[0] != _route_preamble(route_context_sha256) or
            lines[1] != ",".join(campaign.ROUTE_HEADER)):
        die("holdout route manifest preamble/header mismatch")
    records = tuple(
        campaign.parse_route_record(line.split(","), "holdout route row")
        for line in lines[2:])
    by_key: Dict[Tuple[int, int], List[campaign.RouteRecord]] = {}
    previous_grid: Optional[Tuple[int, int]] = None
    grid_order = {(K, width): index for index, (K, width) in enumerate(
        (key for K in ordered_K for key in ((K, width) for width in widths_tuple)))}
    for record in records:
        key = (record.K, record.width)
        if key not in grid_order:
            die("route manifest contains a key outside its exact grid")
        if previous_grid is not None and grid_order[key] < grid_order[previous_grid]:
            die("route manifest key groups are reordered")
        previous_grid = key
        by_key.setdefault(key, []).append(record)
    expected_keys = set(grid_order)
    if set(by_key) != expected_keys:
        die("route manifest does not exactly cover K x widths")
    for K in ordered_K:
        for width in widths_tuple:
            group = by_key[(K, width)]
            if exact_one_per_key and len(group) != 1:
                die("route manifest key does not have exactly one row")
            if len({record.preferred_attempt for record in group}) != len(group):
                die("route manifest repeats a preferred attempt")
            canonical = {record.canonical_attempt for record in group}
            if (len(canonical) != 1 or sum(
                    record.canonical_probe_solves for record in group) !=
                    next(iter(canonical)) + 1):
                die("route manifest canonical/probe accounting changed")
            p = selected[K]
            if p is None:
                if len(group) != 1 or group[0].route_status != "control":
                    die("nonrouted K is not a control route at every width")
            else:
                if any(record.route_status != "preferred" for record in group):
                    die("routed K mixes control and preferred route rows")
                matches = [record for record in group
                           if record.preferred_attempt == p]
                if len(matches) != 1:
                    die("routed K lacks its selected preferred attempt")
    return records


def canonical_selected_route_manifest(
    records: Sequence[campaign.RouteRecord],
    Ks: Sequence[int],
    widths: Sequence[int],
    table: Mapping[int, Optional[int]],
    route_context_sha256: str,
) -> bytes:
    """Produce an exact one-row-per-key cache without doing another probe."""
    ordered_K = tuple(Ks)
    widths_tuple = tuple(widths)
    selected = dict(normalize_table(ordered_K, table))
    expected_keys = {
        (K, width) for K in ordered_K for width in widths_tuple}
    source: Dict[Tuple[int, int], List[campaign.RouteRecord]] = {}
    for record in records:
        campaign.validate_route_record_value(record)
        if (record.K, record.width) not in expected_keys:
            die("route-cache derivation source is outside its exact grid")
        source.setdefault((record.K, record.width), []).append(record)
    if set(source) != expected_keys:
        die("route-cache derivation source does not exactly cover its grid")
    output: List[campaign.RouteRecord] = []
    for K in ordered_K:
        for width in widths_tuple:
            group = source.get((K, width), [])
            if not group:
                die("route-cache derivation lacks a K/width source")
            canonicals = {record.canonical_attempt for record in group}
            if (len(canonicals) != 1 or
                    len({record.preferred_attempt for record in group}) !=
                    len(group)):
                die("route-cache derivation source attempts changed/repeated")
            a0 = next(iter(canonicals))
            if sum(record.canonical_probe_solves for record in group) != a0 + 1:
                die("route-cache derivation source probe accounting changed")
            p = selected[K]
            if p is None:
                chosen = campaign.RouteRecord(
                    K, width, "control", -1, a0, a0, 1, 0, 1, 0,
                    a0 + 1, 0)
            else:
                matches = [record for record in group
                           if record.route_status == "preferred" and
                           record.preferred_attempt == p]
                if len(matches) != 1:
                    die("route-cache derivation lacks selected p")
                old = matches[0]
                chosen = campaign.RouteRecord(
                    old.K, old.width, old.route_status, old.preferred_attempt,
                    old.canonical_attempt, old.actual_attempt,
                    old.preferred_valid, old.fallback, old.no_op, old.direct,
                    old.canonical_attempt + 1, old.preferred_probe_solves)
            campaign.validate_route_record_value(chosen)
            output.append(chosen)
    lines = [
        _route_preamble(route_context_sha256),
        ",".join(campaign.ROUTE_HEADER),
    ]
    for record in output:
        lines.append(",".join(str(value) for value in (
            record.K, record.width, record.route_status,
            record.preferred_attempt, record.canonical_attempt,
            record.actual_attempt, record.preferred_valid, record.fallback,
            record.no_op, record.direct, record.canonical_probe_solves,
            record.preferred_probe_solves)))
    encoded = ("\n".join(lines) + "\n").encode("ascii")
    parse_route_manifest(
        encoded, ordered_K, widths_tuple, selected, route_context_sha256, True)
    return encoded


def _paired_preamble(job: HoldoutJob) -> str:
    if job.kind != "paired" or job.route_cache is None:
        die("paired preamble requested for a nonpaired job")
    return (
        "# preferredattempt: schema=v2 policy=h12-q0-adaptive "
        "canonical=ascending-first-valid-v1 mode=paired loss={} seed={} "
        "schedule={} preferred_batch_max=32 route_cache_sha256={} "
        "route_context_sha256={} probe_route=0 "
        "physical_solve_accounting=explicit "
        "systematic_probe_accounting=explicit"
    ).format(
        campaign.emitted_loss(str(job.loss)), hex(int(str(job.seed), 0)),
        job.schedule, job.route_cache.sha256, job.route_context_sha256)


@dataclass(frozen=True)
class PairedCell:
    K: int
    width: int
    root_index: int
    schedule: str
    loss: str
    preferred_attempt: Optional[int]
    route: Optional[campaign.RouteRecord]
    control: campaign.RecoveryMetrics
    candidate: campaign.RecoveryMetrics


@dataclass(frozen=True)
class ParsedHoldoutOutput:
    route_records: Tuple[campaign.RouteRecord, ...]
    cells: Tuple[PairedCell, ...]
    logical_rows: int
    physical_rows: int
    codec_errors: int
    semantic_sha256: str


def _metrics_tuple(metrics: campaign.RecoveryMetrics) -> Tuple[int, ...]:
    return (
        metrics.result, metrics.rank_fail, metrics.error,
        metrics.heavy_shortfall, metrics.inactivated,
        metrics.binary_deficit, metrics.heavy_gain,
        metrics.block_xors, metrics.block_muladds,
    )


def _validate_metrics(metrics: campaign.RecoveryMetrics, context: str) -> None:
    if not isinstance(metrics, campaign.RecoveryMetrics):
        die(context + " is not a recovery-metrics value")
    values = _metrics_tuple(metrics)
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 0
           for value in values):
        die(context + " contains a noncanonical integer")
    if (metrics.result > U32_MAX or metrics.rank_fail not in (0, 1) or
            metrics.error not in (0, 1) or
            metrics.heavy_shortfall not in (0, 1) or
            metrics.inactivated > U32_MAX or
            metrics.binary_deficit > metrics.inactivated or
            metrics.heavy_gain > U32_MAX or
            metrics.block_xors > U64_MAX or
            metrics.block_muladds > U64_MAX or
            metrics.rank_fail != int(metrics.result == 1) or
            metrics.error != int(metrics.result not in (0, 1)) or
            (metrics.heavy_shortfall and not metrics.rank_fail)):
        die(context + " violates recovery metric invariants")


def _parse_paired_recovery_row(
    fields: Sequence[str], context: str,
) -> campaign.RecoveryRow:
    """Accept the native paired-mode candidate -1 control alias."""
    if len(fields) == len(campaign.RECOVERY_HEADER) and \
            fields[2] == "candidate" and fields[3] == "-1":
        adjusted = list(fields)
        adjusted[3] = "0"
        parsed = campaign.parse_recovery_row(adjusted, context)
        return campaign.RecoveryRow(
            parsed.K, parsed.width, parsed.arm, -1,
            parsed.canonical_attempt, parsed.actual_attempt, parsed.routed,
            parsed.preferred_valid, parsed.fallback, parsed.no_op,
            parsed.direct, parsed.physical_solve, parsed.metrics)
    return campaign.parse_recovery_row(fields, context)


def parse_paired_output(
    job: HoldoutJob, data: bytes, route_data: bytes,
) -> ParsedHoldoutOutput:
    job.validate()
    if job.kind != "paired":
        die("paired parser received a route job")
    mapping = job.routed_map()
    records = parse_route_manifest(
        route_data, job.Ks, job.widths,
        {K: mapping.get(K) for K in job.Ks},
        job.route_context_sha256, True)
    selected_routes: Dict[Tuple[int, int], campaign.RouteRecord] = {}
    by_key: Dict[Tuple[int, int], List[campaign.RouteRecord]] = {}
    for record in records:
        by_key.setdefault((record.K, record.width), []).append(record)
    for K, p in mapping.items():
        for width in job.widths:
            matches = [record for record in by_key[(K, width)]
                       if record.preferred_attempt == p]
            if len(matches) != 1:
                die("paired route cache lost selected p")
            selected_routes[(K, width)] = matches[0]
    lines = campaign.strict_ascii_lines(data, job.job_id + " paired output")
    expected_line_count = 2 + job.expected_logical_rows
    if (len(lines) != expected_line_count or lines[0] != _paired_preamble(job) or
            lines[1] != ",".join(campaign.RECOVERY_HEADER)):
        die("paired output preamble/header/cardinality mismatch")
    cursor = 2
    cells: List[PairedCell] = []
    physical = errors = 0
    digest = hashlib.sha256()
    for K in job.Ks:
        for width in job.widths:
            control_row = _parse_paired_recovery_row(
                lines[cursor].split(","), "{}:row{}".format(job.job_id, cursor + 1))
            cursor += 1
            candidate_row = _parse_paired_recovery_row(
                lines[cursor].split(","), "{}:row{}".format(job.job_id, cursor + 1))
            cursor += 1
            campaign.validate_control_row(control_row, K, width)
            if candidate_row.K != K or candidate_row.width != width or \
                    candidate_row.arm != "candidate":
                die("paired candidate row key/order mismatch")
            p = mapping.get(K)
            route = selected_routes.get((K, width))
            if p is None:
                if ((candidate_row.preferred_attempt, candidate_row.canonical_attempt,
                     candidate_row.actual_attempt, candidate_row.routed,
                     candidate_row.preferred_valid, candidate_row.fallback,
                     candidate_row.no_op, candidate_row.direct,
                     candidate_row.physical_solve) !=
                        (-1, control_row.canonical_attempt,
                         control_row.canonical_attempt, 0, 1, 0, 1, 0, 0)):
                    die("nonrouted paired candidate is not an exact control alias")
            else:
                assert route is not None
                expected_route = (
                    p, route.canonical_attempt, route.actual_attempt, 1,
                    route.preferred_valid, route.fallback, route.no_op,
                    route.direct, route.direct,
                )
                actual_route = (
                    candidate_row.preferred_attempt,
                    candidate_row.canonical_attempt,
                    candidate_row.actual_attempt, candidate_row.routed,
                    candidate_row.preferred_valid, candidate_row.fallback,
                    candidate_row.no_op, candidate_row.direct,
                    candidate_row.physical_solve,
                )
                if actual_route != expected_route or \
                        control_row.canonical_attempt != route.canonical_attempt:
                    die("paired candidate row disagrees with route cache")
            if candidate_row.physical_solve == 0 and \
                    _metrics_tuple(candidate_row.metrics) != _metrics_tuple(
                        control_row.metrics):
                die("logical control alias changed recovery metrics")
            physical += control_row.physical_solve + candidate_row.physical_solve
            errors += control_row.metrics.error + candidate_row.metrics.error
            cell = PairedCell(
                K, width, int(job.root_index), str(job.schedule), str(job.loss),
                p, route, control_row.metrics, candidate_row.metrics)
            cells.append(cell)
            digest.update(_json_bytes({
                "K": K, "width": width, "root_index": job.root_index,
                "schedule": job.schedule, "loss": job.loss, "p": p,
                "route": None if route is None else list(
                    campaign.route_semantics(route)),
                "control": list(_metrics_tuple(control_row.metrics)),
                "candidate": list(_metrics_tuple(candidate_row.metrics)),
            }))
    return ParsedHoldoutOutput(
        (), tuple(cells), job.expected_logical_rows, physical, errors,
        digest.hexdigest())


def parse_holdout_output(
    job: HoldoutJob, data: bytes, route_data: Optional[bytes] = None,
) -> ParsedHoldoutOutput:
    if job.kind == "paired":
        if route_data is None:
            die("paired output lacks bound route bytes")
        return parse_paired_output(job, data, route_data)
    records = parse_route_manifest(
        data, job.Ks, job.widths,
        {K: job.routed_map().get(K) for K in job.Ks},
        job.route_context_sha256, True)
    digest = hashlib.sha256()
    for record in records:
        digest.update(_json_bytes(list(campaign.route_semantics(record))))
    return ParsedHoldoutOutput(
        records, (), len(records), 0, 0, digest.hexdigest())


@dataclass
class ComparisonTotals:
    cells: int = 0
    control_failures: int = 0
    candidate_failures: int = 0
    repairs: int = 0
    introductions: int = 0
    control_heavy_shortfall: int = 0
    candidate_heavy_shortfall: int = 0
    control_deficit_gt15: int = 0
    candidate_deficit_gt15: int = 0
    control_max_deficit: int = 0
    candidate_max_deficit: int = 0
    control_xors: int = 0
    candidate_xors: int = 0
    control_muladds: int = 0
    candidate_muladds: int = 0
    control_inactivated: int = 0
    candidate_inactivated: int = 0
    errors: int = 0

    def add(self, control: campaign.RecoveryMetrics,
            candidate: campaign.RecoveryMetrics) -> None:
        cf = campaign.failure_bit(control) if not control.error else 0
        tf = campaign.failure_bit(candidate) if not candidate.error else 0
        self.cells += 1
        self.errors += control.error + candidate.error
        self.control_failures += cf
        self.candidate_failures += tf
        self.repairs += int(cf == 1 and tf == 0)
        self.introductions += int(cf == 0 and tf == 1)
        self.control_heavy_shortfall += control.heavy_shortfall
        self.candidate_heavy_shortfall += candidate.heavy_shortfall
        self.control_deficit_gt15 += int(control.binary_deficit > 15)
        self.candidate_deficit_gt15 += int(candidate.binary_deficit > 15)
        self.control_max_deficit = max(
            self.control_max_deficit, control.binary_deficit)
        self.candidate_max_deficit = max(
            self.candidate_max_deficit, candidate.binary_deficit)
        if control.success and candidate.success:
            self.control_xors += control.block_xors
            self.candidate_xors += candidate.block_xors
            self.control_muladds += control.block_muladds
            self.candidate_muladds += candidate.block_muladds
            self.control_inactivated += control.inactivated
            self.candidate_inactivated += candidate.inactivated

    def merge(self, other: "ComparisonTotals") -> None:
        for field in (
                "cells", "control_failures", "candidate_failures", "repairs",
                "introductions", "control_heavy_shortfall",
                "candidate_heavy_shortfall", "control_deficit_gt15",
                "candidate_deficit_gt15", "control_xors", "candidate_xors",
                "control_muladds", "candidate_muladds",
                "control_inactivated", "candidate_inactivated", "errors"):
            setattr(self, field, getattr(self, field) + getattr(other, field))
        self.control_max_deficit = max(
            self.control_max_deficit, other.control_max_deficit)
        self.candidate_max_deficit = max(
            self.candidate_max_deficit, other.candidate_max_deficit)

    def work_gate(self, exact_nonworse: bool = False) -> Tuple[bool, bool]:
        pairs = (
            (self.candidate_xors, self.control_xors),
            (self.candidate_muladds, self.control_muladds),
            (self.candidate_inactivated, self.control_inactivated),
        )
        if exact_nonworse:
            return (all(candidate <= control for candidate, control in pairs),
                    any(candidate < control for candidate, control in pairs))
        return (all(200 * candidate <= 201 * control
                    for candidate, control in pairs), False)

    def record(self) -> Dict[str, int]:
        return {
            "cells": self.cells,
            "control_failures": self.control_failures,
            "candidate_failures": self.candidate_failures,
            "repairs": self.repairs, "introductions": self.introductions,
            "control_heavy_shortfall": self.control_heavy_shortfall,
            "candidate_heavy_shortfall": self.candidate_heavy_shortfall,
            "control_deficit_gt15": self.control_deficit_gt15,
            "candidate_deficit_gt15": self.candidate_deficit_gt15,
            "control_max_deficit": self.control_max_deficit,
            "candidate_max_deficit": self.candidate_max_deficit,
            "control_xors": self.control_xors,
            "candidate_xors": self.candidate_xors,
            "control_muladds": self.control_muladds,
            "candidate_muladds": self.candidate_muladds,
            "control_inactivated": self.control_inactivated,
            "candidate_inactivated": self.candidate_inactivated,
            "errors": self.errors,
        }


class HoldoutReducer:
    def __init__(
        self,
        cohort: Sequence[int],
        table: Mapping[int, Optional[int]],
        per_K_axes: bool = False,
        track_seen: bool = True,
    ) -> None:
        self.cohort = tuple(cohort)
        self.table = dict(normalize_table(self.cohort, table))
        self.per_K_axes_enabled = per_K_axes
        # Kept as an argument for callers, but exact coverage is always a
        # compact per-K bit mask.  This avoids a 2.3-million-tuple H2 set.
        self.track_seen = track_seen
        self.per_K: Dict[int, ComparisonTotals] = {
            K: ComparisonTotals() for K in self.cohort}
        self.schedule: Dict[str, ComparisonTotals] = {
            value: ComparisonTotals() for value in SCHEDULES}
        self.loss: Dict[str, ComparisonTotals] = {
            value: ComparisonTotals() for value in LOSSES}
        self.strata: Dict[Tuple[str, str], ComparisonTotals] = {
            (schedule, loss): ComparisonTotals()
            for schedule in SCHEDULES for loss in LOSSES}
        self.slices: Dict[Tuple[int, str, str], ComparisonTotals] = {}
        self.per_K_schedule: Dict[int, Dict[str, ComparisonTotals]] = {}
        self.per_K_loss: Dict[int, Dict[str, ComparisonTotals]] = {}
        self.per_K_slices: Dict[
            int, Dict[Tuple[int, str, str], ComparisonTotals]] = {}
        if per_K_axes:
            self.per_K_schedule = {
                K: {value: ComparisonTotals() for value in SCHEDULES}
                for K in self.cohort}
            self.per_K_loss = {
                K: {value: ComparisonTotals() for value in LOSSES}
                for K in self.cohort}
            self.per_K_slices = {K: {} for K in self.cohort}
        self.routes: Dict[Tuple[int, int], Tuple[Any, ...]] = {}
        self.coverage_masks: Dict[int, int] = {K: 0 for K in self.cohort}
        self.cell_count = 0

    def add(self, cell: PairedCell) -> None:
        if (not isinstance(cell, PairedCell) or
                not isinstance(cell.K, int) or isinstance(cell.K, bool) or
                not isinstance(cell.width, int) or isinstance(cell.width, bool) or
                cell.width not in WIDTHS or
                not isinstance(cell.root_index, int) or
                isinstance(cell.root_index, bool) or
                not 0 <= cell.root_index < H1_ROOT_COUNT):
            die("holdout semantic cell key is malformed")
        _validate_metrics(cell.control, "holdout control metrics")
        _validate_metrics(cell.candidate, "holdout candidate metrics")
        if (cell.K not in self.per_K or cell.schedule not in self.schedule or
                cell.loss not in self.loss or
                cell.preferred_attempt != self.table[cell.K]):
            die("holdout semantic cell is duplicate or outside its frozen grid")
        bit_index = (((WIDTHS.index(cell.width) * H1_ROOT_COUNT +
                       cell.root_index) * len(SCHEDULES) +
                      SCHEDULES.index(cell.schedule)) * len(LOSSES) +
                     LOSSES.index(cell.loss))
        bit = 1 << bit_index
        if self.coverage_masks[cell.K] & bit:
            die("holdout semantic cell is duplicate or outside its frozen grid")
        if self.table[cell.K] is None:
            if cell.route is not None or \
                    _metrics_tuple(cell.control) != _metrics_tuple(cell.candidate):
                die("nonrouted holdout cell is not an exact control alias")
        elif not isinstance(cell.route, campaign.RouteRecord) or \
                cell.route.preferred_attempt != self.table[cell.K]:
            die("routed holdout cell lacks its frozen route")
        elif not cell.route.direct and \
                _metrics_tuple(cell.control) != _metrics_tuple(cell.candidate):
            die("nonphysical routed cell changed recovery metrics")
        self.coverage_masks[cell.K] |= bit
        self.cell_count += 1
        for totals in (
                self.per_K[cell.K], self.schedule[cell.schedule],
                self.loss[cell.loss], self.strata[(cell.schedule, cell.loss)],
                self.slices.setdefault(
                    (cell.width, cell.schedule, cell.loss), ComparisonTotals())):
            totals.add(cell.control, cell.candidate)
        if self.per_K_axes_enabled:
            self.per_K_schedule[cell.K][cell.schedule].add(
                cell.control, cell.candidate)
            self.per_K_loss[cell.K][cell.loss].add(
                cell.control, cell.candidate)
            self.per_K_slices[cell.K].setdefault(
                (cell.width, cell.schedule, cell.loss),
                ComparisonTotals()).add(cell.control, cell.candidate)
        signature: Tuple[Any, ...]
        if cell.route is None:
            signature = (cell.K, cell.width, None, "control")
        else:
            if (cell.route.K, cell.route.width) != (cell.K, cell.width):
                die("holdout cell route key disagrees with the cell")
            signature = campaign.route_semantics(cell.route)
        route_key = (cell.K, cell.width)
        old = self.routes.get(route_key)
        if old is None:
            if cell.route is not None:
                campaign.validate_route_record_value(cell.route)
            self.routes[route_key] = signature
        elif old != signature:
            die("route semantics changed across holdout roots/strata")

    def validate_coverage(
        self, root_count: int, widths: Sequence[int],
    ) -> None:
        expected_count = (
            len(self.cohort) * len(widths) * root_count *
            len(SCHEDULES) * len(LOSSES))
        expected_mask = 0
        for width in widths:
            for root in range(root_count):
                for schedule in SCHEDULES:
                    for loss in LOSSES:
                        bit_index = (((WIDTHS.index(width) * H1_ROOT_COUNT +
                                       root) * len(SCHEDULES) +
                                      SCHEDULES.index(schedule)) * len(LOSSES) +
                                     LOSSES.index(loss))
                        expected_mask |= 1 << bit_index
        if (self.cell_count != expected_count or any(
                value.cells != len(widths) * root_count * 9
                for value in self.per_K.values()) or any(
                mask != expected_mask for mask in self.coverage_masks.values())):
            die("holdout streaming evidence cardinality mismatch")
        if set(self.routes) != {
                (K, width) for K in self.cohort for width in widths}:
            die("holdout route semantics do not exactly cover K x widths")


def _axis_nonworse(values: Mapping[Any, ComparisonTotals]) -> bool:
    return all(value.candidate_failures <= value.control_failures
               for value in values.values())


def _h1_route_gate(reducer: HoldoutReducer, K: int) -> bool:
    p = reducer.table[K]
    if p is None:
        return True
    for width in WIDTHS:
        signature = reducer.routes[(K, width)]
        # Route semantics tuple fields: K,bb,status,p,a0,actual,valid,
        # fallback,no_op,direct.
        if signature[2] != "preferred" or signature[3] != p or signature[7]:
            return False
        if width == 64 and not signature[9]:
            return False
    return True


def _h1_decision(
    reducer: HoldoutReducer, K: int,
) -> Tuple[bool, Tuple[str, ...], Dict[str, Any]]:
    totals = reducer.per_K[K]
    reasons: List[str] = []
    if totals.errors:
        reasons.append("codec_error")
    if not _h1_route_gate(reducer, K):
        reasons.append("route")
    if totals.repairs > totals.introductions:
        repair_gate = True
    elif totals.repairs == totals.introductions == 0:
        exact_nonworse, strict = totals.work_gate(exact_nonworse=True)
        repair_gate = exact_nonworse and strict
        if not exact_nonworse:
            reasons.append("neutral_work_worse")
        elif not strict:
            reasons.append("neutral_without_strict_work_gain")
    else:
        repair_gate = False
        reasons.append(
            "equal_nonzero_discordance" if
            totals.repairs == totals.introductions else
            "repairs_not_greater")
    if not repair_gate and not any(reason in reasons for reason in (
            "neutral_work_worse", "neutral_without_strict_work_gain",
            "equal_nonzero_discordance", "repairs_not_greater")):
        reasons.append("repair_gate")
    work_ok, _unused = totals.work_gate()
    if not work_ok:
        reasons.append("work_1_005")
    if totals.candidate_heavy_shortfall > totals.control_heavy_shortfall:
        reasons.append("heavy_shortfall")
    if totals.candidate_deficit_gt15 > totals.control_deficit_gt15:
        reasons.append("deficit_gt15")
    if totals.candidate_max_deficit > totals.control_max_deficit:
        reasons.append("max_deficit")
    return not reasons, tuple(sorted(set(reasons))), {
        "totals": totals.record(), "work_1_005": work_ok,
        "repair_gate": repair_gate,
    }


def analyze_h1_cells(
    cohort: Sequence[int],
    table: Mapping[int, Optional[int]],
    cells: Iterable[PairedCell],
    strict_domain: bool = True,
) -> Dict[str, Any]:
    if strict_domain and (len(cohort) != H1_K_COUNT or
                          any(K < 4096 for K in cohort)):
        die("H1 analysis cohort is not frozen")
    if strict_domain:
        cell_source = cells
        root_count = H1_ROOT_COUNT
    else:
        materialized = tuple(cells)
        cell_source = materialized
        root_count = 1 + max(
            (cell.root_index for cell in materialized), default=-1)
    reducer = HoldoutReducer(cohort, table, per_K_axes=True)
    for cell in cell_source:
        reducer.add(cell)
    reducer.validate_coverage(root_count, WIDTHS)
    errors_zero = all(value.errors == 0 for value in reducer.per_K.values())
    rows: List[Dict[str, Any]] = []
    output: Dict[int, Optional[int]] = {}
    for K in reducer.cohort:
        p = reducer.table[K]
        if p is None:
            output[K] = None
            reasons = ["already_control"]
            if reducer.per_K[K].errors:
                reasons.append("codec_error")
            rows.append({
                "K": K, "input_preferred_attempt": None,
                "output_preferred_attempt": None, "kept": False,
                "reasons": reasons,
                "totals": reducer.per_K[K].record(),
                "schedule_failures": {
                    name: value.record() for name, value in
                    reducer.per_K_schedule[K].items()},
                "loss_failures": {
                    name: value.record() for name, value in
                    reducer.per_K_loss[K].items()},
                "width_schedule_loss_report": {
                    "{}|{}|{}".format(*key): value.record()
                    for key, value in sorted(reducer.per_K_slices[K].items())},
            })
            continue
        keep, reasons, detail = _h1_decision(reducer, K)
        schedule_ok = _axis_nonworse(reducer.per_K_schedule[K])
        loss_ok = _axis_nonworse(reducer.per_K_loss[K])
        reason_list = list(reasons)
        if not schedule_ok:
            reason_list.append("schedule_total")
        if not loss_ok:
            reason_list.append("loss_total")
        keep = keep and schedule_ok and loss_ok
        reason_list = sorted(set(reason_list))
        output[K] = p if keep else None
        rows.append({
            "K": K, "input_preferred_attempt": p,
            "output_preferred_attempt": output[K], "kept": keep,
            "reasons": reason_list, "totals": detail["totals"],
            "schedule_failures": {
                name: value.record() for name, value in
                reducer.per_K_schedule[K].items()},
            "loss_failures": {
                name: value.record() for name, value in
                reducer.per_K_loss[K].items()},
            "width_schedule_loss_report": {
                "{}|{}|{}".format(*key): value.record()
                for key, value in sorted(reducer.per_K_slices[K].items())},
        })
    table_record = _sealed(
        SCHEMA_PREFIX + ".h1_table.v1",
        {
            "policy": "KPreferredAttemptV1", "K_count": len(rows),
            "rows": [
                {"K": row["K"], "preferred_attempt":
                 row["output_preferred_attempt"]} for row in rows],
            "routed_K_count": sum(value is not None for value in output.values()),
        },
    )
    return _sealed(
        SCHEMA_PREFIX + ".h1_analysis.v1",
        {
            "accepted": errors_zero,
            "gates": {"errors_zero": errors_zero},
            "input_table_sha256": table_semantic_sha256(cohort, table),
            "output_table_sha256": table_semantic_sha256(cohort, output),
            "frozen_table_file_sha256": _sha256(_json_bytes(table_record)),
            "K_count": len(rows), "K": rows,
            "report_only_width_schedule_loss": {
                "{}|{}|{}".format(*key): value.record()
                for key, value in sorted(reducer.slices.items())},
            "frozen_table": table_record,
        },
    )


def parse_h1_table(
    record: Mapping[str, Any], strict_domain: bool = True,
) -> Dict[int, Optional[int]]:
    campaign.verify_sealed_record(record, SCHEMA_PREFIX + ".h1_table.v1")
    if set(record) != {
            "schema", "policy", "K_count", "rows", "routed_K_count",
            "self_sha256_excluding_field"} or \
            record.get("policy") != "KPreferredAttemptV1" or \
            not isinstance(record.get("rows"), list):
        die("H1 frozen table schema/fields mismatch")
    K_count = _require_plain_int(record.get("K_count"), "H1 table K count")
    routed_count = _require_plain_int(
        record.get("routed_K_count"), "H1 table routed K count")
    result: Dict[int, Optional[int]] = {}
    for row in record["rows"]:
        if not isinstance(row, dict) or set(row) != {"K", "preferred_attempt"}:
            die("H1 frozen table row is malformed")
        K = _require_plain_int(row.get("K"), "H1 table K", 2, 64000)
        p = row.get("preferred_attempt")
        if p is not None:
            _require_plain_int(p, "H1 table p", 0, 255)
        if K in result:
            die("H1 frozen table repeats K")
        result[K] = p
    if (tuple(result) != tuple(sorted(result)) or
            K_count != len(result) or
            routed_count !=
                sum(value is not None for value in result.values())):
        die("H1 frozen table arithmetic/order mismatch")
    if strict_domain and (len(result) != H1_K_COUNT or
                          any(K < 4096 for K in result)):
        die("H1 frozen table is not the exact active-K cardinality/domain")
    return result


def expand_h1_table_all_k(
    record: Mapping[str, Any],
) -> Dict[int, Optional[int]]:
    active = parse_h1_table(record, strict_domain=True)
    return {K: active.get(K) for K in ALL_K}


def exact_sign_tail(repairs: int, introductions: int) -> Tuple[int, int, bool]:
    max_cells = H2_LOGICAL_ROWS // 2
    repairs = _require_plain_int(repairs, "repairs", 0, max_cells)
    introductions = _require_plain_int(
        introductions, "introductions", 0, max_cells)
    n = repairs + introductions
    if n > max_cells:
        die("discordant sign-test count exceeds the H2 cell domain")
    # By binomial symmetry the upper tail j=repairs..n is the lower tail
    # j=0..introductions.  Build adjacent coefficients once instead of asking
    # ``math.comb`` to independently reconstruct every enormous integer.
    term = 1
    tail = 1
    for value in range(1, introductions + 1):
        term = term * (n - value + 1) // value
        tail += term
    denominator = 1 << n
    return tail, denominator, 100 * tail <= denominator


def analyze_h2_cells(
    cohort: Sequence[int],
    table: Mapping[int, Optional[int]],
    cells: Iterable[PairedCell],
    route_records: Sequence[campaign.RouteRecord],
    strict_domain: bool = True,
    streaming: bool = False,
) -> Dict[str, Any]:
    if strict_domain and tuple(cohort) != ALL_K:
        die("H2 analysis cohort is not exact all-K")
    if strict_domain:
        root_count = H2_ROOT_COUNT
        cell_source = cells
    else:
        materialized = tuple(cells)
        root_count = 1 + max(
            (cell.root_index for cell in materialized), default=-1)
        cell_source = materialized
    reducer = HoldoutReducer(cohort, table, track_seen=not streaming)
    for cell in cell_source:
        reducer.add(cell)
    reducer.validate_coverage(root_count, (64,))
    expected_route_keys = [
        (K, width) for K in cohort for width in WIDTHS]
    if len(route_records) != len(expected_route_keys):
        die("H2 route evidence cardinality mismatch")
    parsed_routes: List[campaign.RouteRecord] = []
    for expected_key, record in zip(expected_route_keys, route_records):
        campaign.validate_route_record_value(record)
        if (record.K, record.width) != expected_key or \
                record.canonical_probe_solves != record.canonical_attempt + 1:
            die("H2 route evidence order/probe accounting mismatch")
        parsed_routes.append(record)
    route_ok = True
    by_key = {(record.K, record.width): record for record in parsed_routes}
    if len(by_key) != len(parsed_routes):
        die("H2 route evidence repeats a K/width key")
    for K in cohort:
        p = reducer.table[K]
        for width in WIDTHS:
            record = by_key[(K, width)]
            if p is None:
                route_ok = route_ok and record.route_status == "control"
            else:
                route_ok = route_ok and record.route_status == "preferred" and \
                    record.preferred_attempt == p and not bool(record.fallback)
                if width == 64:
                    route_ok = route_ok and bool(record.direct)
            if width == 64:
                paired_signature = reducer.routes[(K, width)]
                expected_signature = (
                    (K, width, None, "control") if p is None else
                    campaign.route_semantics(record))
                route_ok = route_ok and paired_signature == expected_signature
    total = ComparisonTotals()
    routed = ComparisonTotals()
    weak_control = weak_candidate = 0
    multi_control = multi_candidate = 0
    max_control = max_candidate = 0
    for K in cohort:
        value = reducer.per_K[K]
        total.merge(value)
        if reducer.table[K] is not None:
            routed.merge(value)
        cf = value.control_failures
        tf = value.candidate_failures
        weak_control += int(cf >= 1)
        weak_candidate += int(tf >= 1)
        multi_control += int(cf >= 2)
        multi_candidate += int(tf >= 2)
        max_control = max(max_control, cf)
        max_candidate = max(max_candidate, tf)
    tail, _denominator, sign_ok = exact_sign_tail(
        total.repairs, total.introductions)
    all_work, _unused = total.work_gate()
    routed_work, _unused = routed.work_gate()
    gates = {
        "route_validation": route_ok,
        "errors_zero": total.errors == 0,
        "failures_strictly_lower":
            total.candidate_failures < total.control_failures,
        "repairs_gt_introductions": total.repairs > total.introductions,
        "exact_sign_p_le_0_01": sign_ok,
        "schedule_totals_nonworse": _axis_nonworse(reducer.schedule),
        "loss_totals_nonworse": _axis_nonworse(reducer.loss),
        "weak_K_nonworse": weak_candidate <= weak_control,
        "multi_K_nonworse": multi_candidate <= multi_control,
        "max_multiplicity_nonworse": max_candidate <= max_control,
        "heavy_shortfall_nonworse":
            total.candidate_heavy_shortfall <= total.control_heavy_shortfall,
        "deficit_gt15_nonworse":
            total.candidate_deficit_gt15 <= total.control_deficit_gt15,
        "max_deficit_nonworse":
            total.candidate_max_deficit <= total.control_max_deficit,
        "all_K_work_1_005": all_work,
        "routed_work_1_005": routed_work,
    }
    return _sealed(
        SCHEMA_PREFIX + ".h2_analysis.v1",
        {
            "accepted": all(gates.values()), "gates": gates,
            "table_action": "accept_unchanged" if all(gates.values()) else
                "reject_unchanged",
            "input_table_sha256": table_semantic_sha256(cohort, table),
            "output_table_sha256": table_semantic_sha256(cohort, table),
            "K_count": len(cohort), "routed_K_count": sum(
                reducer.table[K] is not None for K in cohort),
            "totals": total.record(), "routed_totals": routed.record(),
            "failure_multiplicity": {
                "control_weak_K": weak_control,
                "candidate_weak_K": weak_candidate,
                "control_multi_K": multi_control,
                "candidate_multi_K": multi_candidate,
                "control_max": max_control, "candidate_max": max_candidate,
            },
            "exact_sign": {
                "repairs": total.repairs,
                "introductions": total.introductions,
                "n": total.repairs + total.introductions,
                "tail_numerator_hex": "0x{:x}".format(tail),
                "denominator": "2^{}".format(total.repairs +
                                                total.introductions),
                "hundred_tail_le_denominator": sign_ok,
            },
            "schedule": {
                key: value.record() for key, value in reducer.schedule.items()},
            "loss": {
                key: value.record() for key, value in reducer.loss.items()},
            "report_only_schedule_loss": {
                "{}|{}".format(*key): value.record()
                for key, value in reducer.strata.items()},
        },
    )


# Execution and persistent sealing -------------------------------------------------


def job_paths(result_root: Path, job: HoldoutJob) -> Dict[str, Path]:
    base = result_root / job.phase / "jobs" / job.kind
    return {
        "stdout": base / "stdout" / (job.job_id + ".csv"),
        "stderr": base / "stderr" / (job.job_id + ".txt"),
        "receipt": base / "receipts" / (job.job_id + ".json"),
        "lock": base / "locks" / (job.job_id + ".lock"),
    }


def _route_bytes_for_job(
    result_root: Path, job: HoldoutJob,
) -> Optional[bytes]:
    if job.kind != "paired":
        return None
    assert job.route_cache is not None
    raw = campaign.stable_result_file(
        result_root, job.route_cache.relative_path, job.job_id + " route cache")
    if _sha256(raw) != job.route_cache.sha256:
        die("paired route cache changed")
    return raw


def _receipt_record(
    result_root: Path,
    binary: Path,
    job: HoldoutJob,
    command: Sequence[str],
    parsed: ParsedHoldoutOutput,
    execution: campaign.ExecutionResult,
) -> Dict[str, Any]:
    paths = job_paths(result_root, job)
    stdout = common.stable_bytes(paths["stdout"])
    stderr = common.stable_bytes(paths["stderr"])
    if stderr:
        die("holdout job stderr is nonempty")
    return _sealed(
        SCHEMA_PREFIX + ".job_receipt.v1",
        {
            "job_id": job.job_id,
            "job_sha256": _sha256(_json_bytes(job.record())),
            "command": list(command), "returncode": 0,
            "cpu": execution.cpu, "start_ns": execution.start_ns,
            "end_ns": execution.end_ns,
            "recovered_after_interrupt": False,
            "stdout": paths["stdout"].relative_to(result_root).as_posix(),
            "stdout_sha256": _sha256(stdout),
            "stderr": paths["stderr"].relative_to(result_root).as_posix(),
            "stderr_sha256": _sha256(stderr),
            "logical_rows": parsed.logical_rows,
            "physical_rows": parsed.physical_rows,
            "codec_errors": parsed.codec_errors,
            "semantic_sha256": parsed.semantic_sha256,
        },
    )


def verify_receipt(
    result_root: Path,
    binary: Path,
    job: HoldoutJob,
    allowed_cpus: Optional[Set[int]] = None,
    taskset_path: Optional[str] = None,
) -> Tuple[Dict[str, Any], ParsedHoldoutOutput]:
    paths = job_paths(result_root, job)
    receipt = campaign.load_canonical_object(
        paths["receipt"], job.job_id + " holdout receipt")
    campaign.verify_sealed_record(receipt, SCHEMA_PREFIX + ".job_receipt.v1")
    expected_fields = {
        "schema", "job_id", "job_sha256", "command", "returncode", "cpu",
        "start_ns", "end_ns", "recovered_after_interrupt", "stdout",
        "stdout_sha256", "stderr", "stderr_sha256", "logical_rows",
        "physical_rows", "codec_errors", "semantic_sha256",
        "self_sha256_excluding_field",
    }
    if set(receipt) != expected_fields:
        die("holdout receipt fields changed")
    for field in (
            "job_sha256", "stdout_sha256", "stderr_sha256",
            "semantic_sha256"):
        _require_hash(receipt.get(field), "receipt " + field)
    returncode = _require_plain_int(
        receipt.get("returncode"), "receipt returncode", 0, 255)
    logical_rows = _require_plain_int(
        receipt.get("logical_rows"), "receipt logical rows")
    physical_rows = _require_plain_int(
        receipt.get("physical_rows"), "receipt physical rows")
    codec_errors = _require_plain_int(
        receipt.get("codec_errors"), "receipt codec errors")
    if (returncode != 0 or physical_rows > logical_rows or
            codec_errors > logical_rows or
            not isinstance(receipt.get("command"), list) or
            any(not isinstance(value, str) or not value
                for value in receipt["command"])):
        die("holdout receipt numeric/command domain changed")
    stdout = common.stable_bytes(paths["stdout"])
    stderr = common.stable_bytes(paths["stderr"])
    route_data = _route_bytes_for_job(result_root, job)
    parsed = parse_holdout_output(job, stdout, route_data)
    recovered = receipt.get("recovered_after_interrupt")
    base = job.command(binary, result_root)
    if taskset_path is None:
        _prepare, contract = campaign.verify_frozen_controller_runtime(
            result_root)
        taskset_path = frozen_taskset_path(contract)
    elif not isinstance(taskset_path, str) or not Path(taskset_path).is_absolute():
        die("receipt verifier taskset path is not absolute")
    if recovered is False:
        cpu = receipt.get("cpu")
        start = receipt.get("start_ns")
        end = receipt.get("end_ns")
        live_fields_ok = (
            isinstance(cpu, int) and not isinstance(cpu, bool) and cpu >= 0 and
            (allowed_cpus is None or cpu in allowed_cpus) and
            isinstance(start, int) and not isinstance(start, bool) and start >= 0 and
            isinstance(end, int) and not isinstance(end, bool) and end >= start)
        expected_command = (taskset_path, "-c", str(cpu), *base)
    else:
        die("holdout receipt cannot claim unmonitored interrupted recovery")
    if (not live_fields_ok or receipt.get("job_id") != job.job_id or
            receipt.get("job_sha256") != _sha256(_json_bytes(job.record())) or
            receipt.get("command") != list(expected_command) or
            returncode != 0 or stderr or
            receipt.get("stdout_sha256") != _sha256(stdout) or
            receipt.get("stderr_sha256") != _sha256(stderr) or
            receipt.get("logical_rows") != parsed.logical_rows or
            receipt.get("physical_rows") != parsed.physical_rows or
            receipt.get("codec_errors") != parsed.codec_errors or
            receipt.get("semantic_sha256") != parsed.semantic_sha256 or
            receipt.get("stdout") !=
                paths["stdout"].relative_to(result_root).as_posix() or
            receipt.get("stderr") !=
                paths["stderr"].relative_to(result_root).as_posix()):
        die("holdout receipt/output semantic binding mismatch")
    return receipt, parsed


def verify_root_binding(result_root: Path, ledger: HoldoutLedger) -> None:
    """Revalidate the frozen controller, rooted drand record, and parent table."""
    import wh2_preferred_attempt_search as search

    result_root = result_root.resolve(strict=True)
    _prepare, contract = campaign.verify_frozen_controller_runtime(result_root)
    spec = search.HOLDOUT_PHASES[ledger.phase]
    state_path = search.seal_state_path(result_root, ledger.phase)
    state = search.load_fixed_json(state_path, ledger.phase + " controller state")
    if state.get("status") != "ROOTED":
        die("holdout phase is not rooted")
    attempt, manifest, seal_record, seal_sha256 = search.load_seal_attempt(
        result_root, spec, state)
    search.validate_loaded_seal(spec, manifest, seal_record, contract)
    search.validate_state_binding(
        state, spec, manifest, seal_record, seal_sha256)
    search.verify_manifest_bytes(result_root, spec, manifest)
    search.verify_seal_publication(
        result_root, state, attempt, seal_record, seal_sha256, contract)
    search.verify_rooted_record(
        result_root, spec, contract, seal_sha256, seal_record)
    search.validate_rooted_state_binding(result_root, spec, state)
    root_path = result_root / spec.root_file
    if common.sha256_file(root_path) != ledger.root_file_sha256:
        die("holdout root file changed after ledger construction")
    parent = result_root / (
        "development/frozen_table.json" if ledger.phase == "h1" else
        "h1/frozen_table.json")
    if common.sha256_file(parent) != ledger.parent_table_sha256:
        die("holdout parent table changed")
    if (common.sha256_file(
            result_root / "frozen/route_context.json") !=
            ledger.route_context_sha256):
        die("holdout route context changed")


Executor = Callable[
    [Tuple[str, ...], float, threading.Event, common.ProcessRegistry, int],
    campaign.ExecutionResult,
]
ThermalFactory = Callable[[Path, threading.Event], Any]
RuntimeVerifier = Callable[[Path], Tuple[Dict[str, Any], Dict[str, Any]]]


def frozen_taskset_path(contract: Mapping[str, Any]) -> str:
    return str(campaign.frozen_taskset_path(contract))


def phase_thermal_path(result_root: Path, ledger: HoldoutLedger) -> Path:
    return result_root / "thermal" / (
        ledger.phase + "-" + ledger.kind + ".csv")


def _phase_stem(ledger: HoldoutLedger) -> str:
    return "phase_complete" if ledger.phase == "h1" else \
        "phase_" + ledger.kind


def _discard_unsealed_phase_artifacts(
    result_root: Path, ledger: HoldoutLedger,
) -> None:
    """Do not splice old receipts into a newly monitored thermal interval."""
    paths: List[Path] = []
    for job in ledger.jobs:
        job_files = job_paths(result_root, job)
        for name in ("receipt", "stdout", "stderr"):
            if campaign.path_present(job_files[name]):
                paths.append(job_files[name])
    thermal = phase_thermal_path(result_root, ledger)
    if campaign.path_present(thermal):
        paths.append(thermal)
    parents: Set[Path] = set()
    for path in paths:
        common.stable_bytes(path)
    for path in paths:
        path.unlink()
        parents.add(path.parent)
    for parent in parents:
        descriptor = os.open(
            str(parent), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


class HoldoutRunner:
    """Parallel runner with immutable completed-phase resume and telemetry."""

    def __init__(
        self,
        result_root: Path,
        binary: Path,
        cpus: Sequence[int],
        workers: int,
        timeout: float,
        executor: Optional[Executor] = None,
        binding_verifier: Callable[[Path, HoldoutLedger], None] =
            verify_root_binding,
        runtime_verifier: RuntimeVerifier =
            campaign.verify_frozen_controller_runtime,
        thermal_factory: Optional[ThermalFactory] = None,
    ) -> None:
        if (not isinstance(workers, int) or isinstance(workers, bool) or
                workers <= 0 or workers != len(cpus) or
                len(cpus) != len(set(cpus)) or
                any(not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
                    for cpu in cpus) or
                not isinstance(timeout, (int, float)) or isinstance(timeout, bool) or
                not math.isfinite(float(timeout)) or timeout <= 0):
            die("holdout runner CPU/worker/timeout policy is invalid")
        self.result_root = result_root.resolve()
        self.binary = binary.absolute()
        self.cpus = tuple(cpus)
        self.workers = workers
        self.timeout = float(timeout)
        self.executor = executor or campaign.subprocess_executor
        self.binding_verifier = binding_verifier
        self.runtime_verifier = runtime_verifier
        self.thermal_factory = thermal_factory

    def _run_one(
        self,
        job: HoldoutJob,
        pool: common.CpuPool,
        abort: threading.Event,
        registry: common.ProcessRegistry,
        taskset: str,
        contract: Mapping[str, Any],
    ) -> Tuple[Dict[str, Any], ParsedHoldoutOutput]:
        paths = job_paths(self.result_root, job)
        with campaign.job_lock(paths["lock"]):
            if campaign.path_present(paths["receipt"]):
                return verify_receipt(
                    self.result_root, self.binary, job, set(self.cpus), taskset)
            if (campaign.path_present(paths["stdout"]) or
                    campaign.path_present(paths["stderr"])):
                campaign.discard_partial_output_pair(paths, job.job_id)
            base = job.command(self.binary, self.result_root)
            if abort.is_set():
                die("holdout abort set before launch")
            _route_bytes_for_job(self.result_root, job)
            cpu = pool.acquire()
            try:
                if abort.is_set():
                    die("holdout abort set while waiting for a CPU")
                campaign.verify_frozen_binary(self.binary, contract)
                if frozen_taskset_path(contract) != taskset:
                    die("frozen taskset changed before holdout job launch")
                command = (taskset, "-c", str(cpu), *base)
                execution = self.executor(
                    tuple(command), self.timeout, abort, registry, cpu)
            finally:
                pool.release(cpu)
            if (not isinstance(execution, campaign.ExecutionResult) or
                    not isinstance(execution.stdout, bytes) or
                    not isinstance(execution.stderr, bytes) or
                    not isinstance(execution.returncode, int) or
                    isinstance(execution.returncode, bool) or
                    not isinstance(execution.cpu, int) or
                    isinstance(execution.cpu, bool) or execution.cpu != cpu or
                    not isinstance(execution.start_ns, int) or
                    isinstance(execution.start_ns, bool) or
                    not isinstance(execution.end_ns, int) or
                    isinstance(execution.end_ns, bool) or
                    execution.start_ns < 0 or
                    execution.end_ns < execution.start_ns):
                die("{} executor result is malformed".format(job.job_id))
            if execution.returncode != 0 or execution.stderr:
                die("{} failed with rc={}".format(
                    job.job_id, execution.returncode))
            campaign.durable_atomic_write(paths["stderr"], execution.stderr)
            campaign.durable_atomic_write(paths["stdout"], execution.stdout)
            parsed = parse_holdout_output(
                job, execution.stdout, _route_bytes_for_job(self.result_root, job))
            receipt = _receipt_record(
                self.result_root, self.binary, job, command, parsed,
                execution)
            campaign.durable_atomic_write(paths["receipt"], _json_bytes(receipt))
            return verify_receipt(
                self.result_root, self.binary, job, set(self.cpus), taskset)

    def run_ledger(
        self,
        ledger: HoldoutLedger,
        thermal_path: Optional[Path] = None,
        install_signal_handlers: bool = True,
        exact: bool = True,
    ) -> Tuple[Tuple[Dict[str, Any], ParsedHoldoutOutput], ...]:
        ledger.validate(exact=exact)
        # Cross the frozen/root trust boundary before creating even a lock.
        self.binding_verifier(self.result_root, ledger)
        self.runtime_verifier(self.result_root)
        phase_lock = self.result_root / "locks" / (
            ledger.phase + "-" + ledger.kind + ".lock")
        with campaign.job_lock(phase_lock):
            return self._run_ledger_locked(
                ledger, thermal_path, install_signal_handlers, exact)

    def _run_ledger_locked(
        self,
        ledger: HoldoutLedger,
        thermal_path: Optional[Path],
        install_signal_handlers: bool,
        exact: bool,
    ) -> Tuple[Tuple[Dict[str, Any], ParsedHoldoutOutput], ...]:
        ledger.validate(exact=exact)
        self.binding_verifier(self.result_root, ledger)
        _prepare, contract = self.runtime_verifier(self.result_root)
        if not isinstance(contract, dict) or not isinstance(thermal_path, Path):
            die("holdout runner requires the frozen runtime and thermal log")
        expected_cpus = contract.get("cpu_set")
        expected_workers = contract.get("workers")
        expected_timeout = contract.get("timeout_seconds")
        expected_thermal = contract.get("thermal")
        thermal_policy = campaign.validate_thermal_policy(contract)
        taskset = frozen_taskset_path(contract)
        if (not isinstance(expected_cpus, list) or
                tuple(expected_cpus) != self.cpus or
                any(not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
                    for cpu in expected_cpus) or
                not isinstance(expected_workers, int) or
                isinstance(expected_workers, bool) or
                expected_workers != self.workers or
                not isinstance(expected_timeout, (int, float)) or
                isinstance(expected_timeout, bool) or
                float(expected_timeout) != self.timeout or
                not isinstance(expected_thermal, str)):
            die("holdout CPU/worker/timeout policy differs from frozen runtime")
        try:
            actual_thermal = thermal_path.resolve(strict=True)
            frozen_thermal = Path(expected_thermal).resolve(strict=True)
        except OSError as error:
            die("frozen holdout thermal log is unavailable: {}".format(error))
        if actual_thermal != frozen_thermal:
            die("holdout thermal log differs from frozen runtime")
        frozen_thermal_baseline: Optional[Tuple[int, int, int, int]] = None
        if self.thermal_factory is None:
            frozen_thermal_baseline = campaign.validate_frozen_thermal_source(
                contract, actual_thermal, thermal_policy)
        expected_binary = (
            self.result_root / "frozen/wirehair_v2_bench").resolve(strict=True)
        if (self.binary.resolve(strict=True) != expected_binary or
                self.binary != expected_binary):
            die("holdout runner binary is not the frozen binary")
        campaign.verify_frozen_binary(self.binary, contract)
        ledger_path = self.result_root / ledger.phase / "ledgers" / (
            ledger.kind + ".json")
        campaign.write_hashed_artifact(
            ledger_path, _json_bytes(ledger.record(exact=exact)))
        phase_path = self.result_root / ledger.phase / (
            _phase_stem(ledger) + ".json")
        phase_sidecar = phase_path.with_suffix(".sha256")
        if campaign.path_present(phase_path):
            _load_phase_completion(
                self.result_root, ledger, self.binary, set(self.cpus),
                exact=exact, thermal_policy=thermal_policy,
                taskset_path=taskset, repair_missing_sidecar=True)
            return tuple(
                verify_receipt(
                    self.result_root, self.binary, job, set(self.cpus), taskset)
                for job in ledger.jobs)
        if campaign.path_present(phase_sidecar):
            die("holdout phase sidecar exists without its sealed record")

        _discard_unsealed_phase_artifacts(self.result_root, ledger)
        if not ledger.jobs:
            completion = phase_completion_record(
                self.result_root, self.binary, ledger, set(self.cpus),
                taskset_path=taskset, exact=exact)
            _publish_named_record(
                self.result_root / ledger.phase, _phase_stem(ledger), completion)
            return ()
        abort = threading.Event()
        registry = common.ProcessRegistry()
        pool = common.CpuPool(self.cpus, self.workers)
        guard = common.CampaignSignalGuard(abort, registry) \
            if install_signal_handlers else nullcontext()
        completed: Dict[str, Tuple[Dict[str, Any], ParsedHoldoutOutput]] = {}
        thermal = None
        thermal_summary: Optional[Dict[str, Any]] = None
        run_error: Optional[BaseException] = None
        cleanup_error: Optional[BaseException] = None
        with guard:
            executor: Optional[ThreadPoolExecutor] = None
            futures: Dict[Future[Any], HoldoutJob] = {}
            try:
                if self.thermal_factory is None:
                    if frozen_thermal_baseline is None:
                        die("frozen holdout thermal baseline is unavailable")
                    thermal = common.ThermalGuard(
                        actual_thermal, abort,
                        stale_seconds=float(thermal_policy["stale_seconds"]),
                        limit_c=float(thermal_policy["cpu_limit_c"]),
                        consecutive_limit=thermal_policy["consecutive_samples"],
                        expected_dev=frozen_thermal_baseline[0],
                        expected_ino=frozen_thermal_baseline[1],
                        expected_edac_ce=frozen_thermal_baseline[2],
                        expected_edac_ue=frozen_thermal_baseline[3])
                    if campaign.validate_frozen_thermal_source(
                            contract, actual_thermal,
                            thermal_policy) != frozen_thermal_baseline:
                        die(
                            "frozen holdout thermal baseline changed before "
                            "guard start")
                else:
                    thermal = self.thermal_factory(actual_thermal, abort)
                thermal.start()
                executor = ThreadPoolExecutor(max_workers=self.workers)
                futures = {
                    executor.submit(
                        self._run_one, job, pool, abort, registry, taskset,
                        contract): job
                    for job in ledger.jobs}
                for future in as_completed(futures):
                    completed[futures[future].job_id] = future.result()
                    if getattr(thermal, "error", None) is not None:
                        die("holdout thermal guard failed: {}".format(
                            thermal.error))
            except BaseException as error:
                run_error = error
                abort.set()
                registry.signal_all(signal.SIGTERM)
                for future in futures:
                    future.cancel()
            finally:
                try:
                    if executor is not None:
                        executor.shutdown(wait=True)
                    if registry.count() != 0:
                        registry.drain()
                except BaseException as error:
                    cleanup_error = error
                finally:
                    if thermal is not None and getattr(thermal, "started", False):
                        try:
                            output = phase_thermal_path(self.result_root, ledger)
                            parent_fd = common.open_durable_directory(
                                output.parent, create=True)
                            os.close(parent_fd)
                            thermal_summary = thermal.finish(output)
                            campaign.fsync_existing_regular_file(output)
                        except BaseException as error:
                            cleanup_error = error
        if cleanup_error is not None:
            raise cleanup_error from run_error
        if thermal_summary is None:
            if run_error is not None:
                raise run_error
            die("holdout thermal guard did not produce a summary")
        validated_summary = campaign.validate_thermal_summary(
            thermal_summary, thermal_policy, require_busy=run_error is None)
        if run_error is not None:
            raise run_error
        if registry.count() != 0:
            die("holdout phase left registered child processes")
        self.binding_verifier(self.result_root, ledger)
        _end_prepare, end_contract = self.runtime_verifier(self.result_root)
        campaign.verify_frozen_binary(self.binary, end_contract)
        if frozen_taskset_path(end_contract) != taskset:
            die("frozen holdout taskset identity changed during the phase")
        if self.thermal_factory is None:
            end_policy = campaign.validate_thermal_policy(end_contract)
            if (end_policy != thermal_policy or
                    frozen_thermal_baseline is None or
                    campaign.validate_frozen_thermal_source(
                        end_contract, actual_thermal, end_policy) !=
                    frozen_thermal_baseline):
                die("frozen holdout thermal history changed during the phase")
        ordered = tuple(completed[job.job_id] for job in ledger.jobs)
        if len(ordered) != len(ledger.jobs):
            die("holdout phase completion cardinality mismatch")
        completion = phase_completion_record(
            self.result_root, self.binary, ledger, set(self.cpus),
            phase_thermal_path(self.result_root, ledger), validated_summary,
            thermal_policy, taskset, exact=exact)
        _publish_named_record(
            self.result_root / ledger.phase, _phase_stem(ledger), completion)
        return ordered


def phase_completion_record(
    result_root: Path,
    binary: Path,
    ledger: HoldoutLedger,
    allowed_cpus: Set[int],
    thermal_path: Optional[Path] = None,
    thermal_summary: Optional[Mapping[str, Any]] = None,
    thermal_policy: Optional[Mapping[str, Any]] = None,
    taskset_path: Optional[str] = None,
    exact: bool = True,
) -> Dict[str, Any]:
    ledger.validate(exact=exact)
    if taskset_path is None and ledger.jobs:
        _prepare, contract = campaign.verify_frozen_controller_runtime(
            result_root)
        taskset_path = frozen_taskset_path(contract)
    elif taskset_path is not None and (
            not isinstance(taskset_path, str) or
            not Path(taskset_path).is_absolute()):
        die("holdout completion taskset path is not absolute")
    jobs: List[Dict[str, Any]] = []
    logical = physical = errors = 0
    semantic = hashlib.sha256()
    for job in ledger.jobs:
        receipt, parsed = verify_receipt(
            result_root, binary, job, allowed_cpus, taskset_path)
        receipt_bytes = _json_bytes(receipt)
        jobs.append({
            "job_id": job.job_id, "receipt_sha256": _sha256(receipt_bytes),
            "stdout_sha256": receipt["stdout_sha256"],
            "logical_rows": parsed.logical_rows,
            "physical_rows": parsed.physical_rows,
            "codec_errors": parsed.codec_errors,
            "semantic_sha256": parsed.semantic_sha256,
        })
        logical += parsed.logical_rows
        physical += parsed.physical_rows
        errors += parsed.codec_errors
        semantic.update(bytes.fromhex(parsed.semantic_sha256))
    if logical != sum(job.expected_logical_rows for job in ledger.jobs):
        die("holdout phase logical row arithmetic mismatch")
    if ledger.jobs:
        expected_thermal = phase_thermal_path(result_root, ledger)
        if (thermal_path is None or thermal_summary is None or
                thermal_path.absolute() != expected_thermal.absolute()):
            die("nonempty holdout completion lacks exact thermal provenance")
        if thermal_policy is None:
            die("nonempty holdout completion lacks frozen thermal policy")
        summary = campaign.validate_thermal_summary(
            thermal_summary, thermal_policy)
        thermal_bytes = common.stable_bytes(expected_thermal)
        thermal_artifact: Optional[Dict[str, Any]] = {
            "path": expected_thermal.relative_to(result_root).as_posix(),
            "sha256": _sha256(thermal_bytes), "summary": summary,
        }
    else:
        if thermal_path is not None or thermal_summary is not None:
            die("empty holdout completion unexpectedly claims telemetry")
        thermal_artifact = None
    return _sealed(
        SCHEMA_PREFIX + ".phase_completion.v1",
        {
            "phase": ledger.phase, "kind": ledger.kind,
            "job_ledger_sha256": ledger.sha256(exact=exact),
            "job_count": len(jobs), "jobs": jobs,
            "logical_rows": logical, "physical_rows": physical,
            "codec_errors": errors,
            "taskset_path": taskset_path,
            "thermal_artifact": thermal_artifact,
            "semantic_stream_sha256": semantic.hexdigest(),
        },
    )


def _publish_named_record(directory: Path, stem: str,
                          record: Mapping[str, Any]) -> str:
    encoded = _json_bytes(record)
    digest = _sha256(encoded)
    json_path = directory / (stem + ".json")
    campaign.durable_write_once_or_same(json_path, encoded)
    sidecar = (digest + "  " + json_path.name + "\n").encode("ascii")
    campaign.durable_write_once_or_same(directory / (stem + ".sha256"), sidecar)
    return digest


def _load_phase_completion(
    result_root: Path, ledger: HoldoutLedger, binary: Path,
    allowed_cpus: Set[int], exact: bool = True,
    thermal_policy: Optional[Mapping[str, Any]] = None,
    taskset_path: Optional[str] = None,
    repair_missing_sidecar: bool = False,
) -> Dict[str, Any]:
    stem = _phase_stem(ledger)
    path = result_root / ledger.phase / (stem + ".json")
    stored = campaign.load_canonical_object(path, "holdout phase completion")
    campaign.verify_sealed_record(
        stored, SCHEMA_PREFIX + ".phase_completion.v1")
    stored_bytes = _json_bytes(stored)
    expected_sidecar = (
        _sha256(stored_bytes) + "  " + path.name + "\n").encode("ascii")
    sidecar = path.with_suffix(".sha256")
    sidecar_present = campaign.path_present(sidecar)
    if (sidecar_present and
            common.stable_bytes(sidecar) != expected_sidecar):
        die("holdout phase completion sidecar changed")
    if not sidecar_present and not repair_missing_sidecar:
        die("holdout phase completion sidecar is missing")
    if thermal_policy is None or taskset_path is None:
        _prepare, contract = campaign.verify_frozen_controller_runtime(
            result_root)
        if thermal_policy is None:
            thermal_policy = campaign.validate_thermal_policy(contract)
        if taskset_path is None:
            taskset_path = frozen_taskset_path(contract)
    artifact = stored.get("thermal_artifact")
    if ledger.jobs:
        expected_thermal = phase_thermal_path(result_root, ledger)
        expected_relative = expected_thermal.relative_to(result_root).as_posix()
        if (not isinstance(artifact, dict) or
                set(artifact) != {"path", "sha256", "summary"} or
                artifact.get("path") != expected_relative):
            die("holdout phase thermal binding is malformed")
        _require_hash(artifact.get("sha256"), "phase thermal artifact hash")
        if _sha256(common.stable_bytes(expected_thermal)) != artifact["sha256"]:
            die("completed holdout thermal artifact changed")
        summary = campaign.validate_thermal_summary(
            artifact.get("summary"), thermal_policy)
        thermal_path: Optional[Path] = expected_thermal
    else:
        if artifact is not None:
            die("empty completed holdout phase claims telemetry")
        thermal_path = None
        summary = None
    recomputed = phase_completion_record(
        result_root, binary, ledger, allowed_cpus,
        thermal_path, summary, thermal_policy, taskset_path, exact=exact)
    if stored != recomputed:
        die("holdout phase completion changed after publication")
    if not sidecar_present:
        # Publication orders the canonical record before its checksum.  Heal
        # only the missing checksum after the complete semantic record,
        # receipts, route bindings, and telemetry have been recomputed.  A
        # present but incorrect checksum remains a permanent error.
        campaign.durable_write_once_or_same(sidecar, expected_sidecar)
    return stored


def _all_cells_from_ledger(
    result_root: Path, binary: Path, ledger: HoldoutLedger,
    allowed_cpus: Set[int], exact: bool = True,
    taskset_path: Optional[str] = None,
) -> Tuple[PairedCell, ...]:
    if taskset_path is None:
        _prepare, contract = campaign.verify_frozen_controller_runtime(result_root)
        taskset_path = frozen_taskset_path(contract)
    _load_phase_completion(
        result_root, ledger, binary, allowed_cpus, exact=exact,
        taskset_path=taskset_path)
    cells: List[PairedCell] = []
    for job in ledger.jobs:
        _receipt, parsed = verify_receipt(
            result_root, binary, job, allowed_cpus, taskset_path)
        cells.extend(parsed.cells)
    return tuple(cells)


def _iter_cells_from_ledger(
    result_root: Path, binary: Path, ledger: HoldoutLedger,
    allowed_cpus: Set[int], exact: bool = True,
    verify_completion: bool = True,
    taskset_path: Optional[str] = None,
) -> Iterable[PairedCell]:
    if taskset_path is None:
        _prepare, contract = campaign.verify_frozen_controller_runtime(result_root)
        taskset_path = frozen_taskset_path(contract)
    if verify_completion:
        _load_phase_completion(
            result_root, ledger, binary, allowed_cpus, exact=exact,
            taskset_path=taskset_path)
    for job in ledger.jobs:
        _receipt, parsed = verify_receipt(
            result_root, binary, job, allowed_cpus, taskset_path)
        for cell in parsed.cells:
            yield cell


def analyze_h1(
    result_root: Path,
    binary: Path,
    ledger: HoldoutLedger,
    cohort: Sequence[int],
    table: Mapping[int, Optional[int]],
    allowed_cpus: Set[int],
) -> Dict[str, Any]:
    ledger.validate()
    if (ledger.phase, ledger.kind) != ("h1", "paired"):
        die("H1 analyzer received another phase")
    verify_root_binding(result_root, ledger)
    cells = _all_cells_from_ledger(
        result_root, binary, ledger, allowed_cpus)
    directory = result_root / "h1"
    phase_bytes = common.stable_bytes(directory / "phase_complete.json")
    analysis = _bind_analysis(
        analyze_h1_cells(cohort, table, cells),
        SCHEMA_PREFIX + ".h1_analysis.v1",
        {
            "phase_completion_sha256": _sha256(phase_bytes),
            "root_file_sha256": ledger.root_file_sha256,
            "parent_table_sha256": ledger.parent_table_sha256,
            "route_context_sha256": ledger.route_context_sha256,
            "input_table_sha256": table_semantic_sha256(cohort, table),
        },
    )
    _publish_named_record(directory, "analysis_complete", analysis)
    frozen_table = analysis["frozen_table"]
    campaign.durable_write_once_or_same(
        directory / "frozen_table.json", _json_bytes(frozen_table))
    return analysis


def analyze_h2(
    result_root: Path,
    binary: Path,
    route_ledger: HoldoutLedger,
    paired_ledger: HoldoutLedger,
    table: Mapping[int, Optional[int]],
    allowed_cpus: Set[int],
) -> Dict[str, Any]:
    route_ledger.validate()
    paired_ledger.validate()
    if ((route_ledger.phase, route_ledger.kind) != ("h2", "route") or
            (paired_ledger.phase, paired_ledger.kind) != ("h2", "paired") or
            route_ledger.root_file_sha256 != paired_ledger.root_file_sha256 or
            route_ledger.parent_table_sha256 != paired_ledger.parent_table_sha256):
        die("H2 analyzer ledger bindings disagree")
    verify_root_binding(result_root, route_ledger)
    verify_root_binding(result_root, paired_ledger)
    _prepare, contract = campaign.verify_frozen_controller_runtime(result_root)
    thermal_policy = campaign.validate_thermal_policy(contract)
    taskset = frozen_taskset_path(contract)
    route_completion = _load_phase_completion(
        result_root, route_ledger, binary, allowed_cpus,
        thermal_policy=thermal_policy, taskset_path=taskset)
    paired_completion = _load_phase_completion(
        result_root, paired_ledger, binary, allowed_cpus,
        thermal_policy=thermal_policy, taskset_path=taskset)
    combined_completion = _sealed(
        SCHEMA_PREFIX + ".h2_phase_completion.v1",
        {
            "phase": "h2",
            "root_file_sha256": route_ledger.root_file_sha256,
            "parent_table_sha256": route_ledger.parent_table_sha256,
            "route_completion_sha256": _sha256(_json_bytes(route_completion)),
            "paired_completion_sha256": _sha256(_json_bytes(paired_completion)),
            "route_rows": route_completion["logical_rows"],
            "recovery_logical_rows": paired_completion["logical_rows"],
            "codec_errors": (
                route_completion["codec_errors"] +
                paired_completion["codec_errors"]),
        },
    )
    _publish_named_record(
        result_root / "h2", "phase_complete", combined_completion)
    cells = _iter_cells_from_ledger(
        result_root, binary, paired_ledger, allowed_cpus,
        verify_completion=False, taskset_path=taskset)
    routes: List[campaign.RouteRecord] = []
    for job in route_ledger.jobs:
        _receipt, parsed = verify_receipt(
            result_root, binary, job, allowed_cpus, taskset)
        routes.extend(parsed.route_records)
    routes.sort(key=lambda record: (record.K, WIDTHS.index(record.width)))
    analysis = _bind_analysis(
        analyze_h2_cells(
            ALL_K, table, cells, routes, streaming=True),
        SCHEMA_PREFIX + ".h2_analysis.v1",
        {
            "phase_completion_sha256": _sha256(
                _json_bytes(combined_completion)),
            "route_completion_sha256": _sha256(
                _json_bytes(route_completion)),
            "paired_completion_sha256": _sha256(
                _json_bytes(paired_completion)),
            "root_file_sha256": route_ledger.root_file_sha256,
            "parent_table_sha256": route_ledger.parent_table_sha256,
            "route_context_sha256": route_ledger.route_context_sha256,
            "input_table_sha256": table_semantic_sha256(ALL_K, table),
        },
    )
    _publish_named_record(result_root / "h2", "analysis_complete", analysis)
    return analysis


def derive_h2_bb64_route_caches(
    result_root: Path,
    binary: Path,
    route_ledger: HoldoutLedger,
    table: Mapping[int, Optional[int]],
    allowed_cpus: Set[int],
) -> Dict[int, RouteCacheBinding]:
    """Verify the fresh four-width phase, then publish exact bb64 caches."""
    route_ledger.validate()
    if (route_ledger.phase, route_ledger.kind) != ("h2", "route"):
        die("H2 cache derivation received another phase")
    normalize_table(ALL_K, table)
    verify_root_binding(result_root, route_ledger)
    _prepare, contract = campaign.verify_frozen_controller_runtime(result_root)
    thermal_policy = campaign.validate_thermal_policy(contract)
    taskset = frozen_taskset_path(contract)
    completion = _load_phase_completion(
        result_root, route_ledger, binary, allowed_cpus,
        thermal_policy=thermal_policy, taskset_path=taskset)
    bindings: Dict[int, RouteCacheBinding] = {}
    index_rows: List[Dict[str, Any]] = []
    for job in route_ledger.jobs:
        _receipt, parsed = verify_receipt(
            result_root, binary, job, allowed_cpus, taskset)
        bb64 = tuple(
            record for record in parsed.route_records if record.width == 64)
        selected = {K: table[K] for K in job.Ks}
        relative = "h2/route_cache_bb64/bin-{:03d}.csv".format(
            job.bin_index)
        binding = publish_derived_route_cache(
            result_root, relative, bb64, job.Ks, (64,), selected,
            route_ledger.route_context_sha256)
        bindings[job.bin_index] = binding
        index_rows.append({
            "bin": job.bin_index, "K_count": len(job.Ks),
            "path": binding.relative_path, "sha256": binding.sha256,
        })
    if set(bindings) != set(range(BIN_COUNT)):
        die("H2 derived bb64 cache index is incomplete")
    index_record = _sealed(
        SCHEMA_PREFIX + ".h2_bb64_route_cache_index.v1",
        {
            "route_ledger_sha256": route_ledger.sha256(),
            "route_phase_completion_sha256": _sha256(
                _json_bytes(completion)),
            "input_table_sha256": table_semantic_sha256(ALL_K, table),
            "cache_count": len(index_rows), "caches": index_rows,
        },
    )
    _publish_named_record(
        result_root / "h2/route_cache_bb64", "index", index_record)
    return bindings


def publish_derived_route_cache(
    result_root: Path,
    relative_path: str,
    records: Sequence[campaign.RouteRecord],
    Ks: Sequence[int],
    widths: Sequence[int],
    table: Mapping[int, Optional[int]],
    route_context_sha256: str,
) -> RouteCacheBinding:
    relative = _safe_relative(relative_path, "derived route-cache path")
    encoded = canonical_selected_route_manifest(
        records, Ks, widths, table, route_context_sha256)
    path = result_root / relative
    digest = campaign.write_hashed_artifact(path, encoded)
    return RouteCacheBinding(relative, digest)


__all__ = [
    "ALL_K", "BIN_COUNT", "H1_JOB_COUNT", "H1_LOGICAL_ROWS",
    "H2_JOB_COUNT", "H2_LOGICAL_ROWS", "H2_ROUTE_JOB_COUNT",
    "H2_ROUTE_ROWS", "CampaignError", "ComparisonTotals", "HoldoutJob",
    "HoldoutLedger", "HoldoutReducer", "HoldoutRunner", "PairedCell",
    "ParsedHoldoutOutput", "RouteCacheBinding", "analyze_h1",
    "analyze_h1_cells", "analyze_h2", "analyze_h2_cells",
    "build_h1_ledger", "build_h2_paired_ledger", "build_h2_route_ledger",
    "canonical_selected_route_manifest", "derive_h2_bb64_route_caches",
    "exact_sign_tail", "expand_h1_table_all_k", "frozen_taskset_path",
    "parse_h1_table",
    "parse_holdout_output", "parse_paired_output",
    "parse_route_manifest", "phase_completion_record", "phase_thermal_path",
    "publish_derived_route_cache", "table_semantic_sha256",
    "verify_receipt", "verify_root_binding",
]
