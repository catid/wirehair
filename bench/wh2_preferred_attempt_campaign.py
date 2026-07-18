#!/usr/bin/env python3
"""Deterministic execution and analysis primitives for the WH2 H12 search.

This module deliberately contains no policy-discovery shortcuts.  The caller
supplies a frozen cohort/bin assignment, roots, and the survivor seal from the
previous round.  It returns canonical job ledgers and a canonical next-round
survivor seal.  The preferred-attempt controller can import these primitives
after it has verified its own frozen trust anchors.

The implementation is Python 3.8 compatible.  Recovery jobs use the existing
``wirehair_v2_bench preferredattempt`` v2 CSV contract and the process, signal,
thermal, and atomic-file safety primitives from the all-K campaign runner.
"""

from __future__ import annotations

import sys

# Frozen campaign directories have an exact immutable inventory.  Disable
# import-cache writes before loading any local frozen helper module.
sys.dont_write_bytecode = True

from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import stat
import subprocess
import threading
import time
from typing import (
    Any, Callable, Dict, Iterator, List, Mapping, Optional, Sequence,
    Set, Tuple,
)

import wh2_rank_floor_two_anchor_allk as common


CampaignError = common.CampaignError
die = common.die

SCHEMA_PREFIX = "wirehair.wh2.h12_preferred_attempt"
APPROVED_WIDTHS = (64, 256, 1280, 4096)
APPROVED_SCHEDULES = ("burst", "adversarial", "repair-only")
APPROVED_LOSSES = ("0.35", "0.50", "0.65")
MAX_ATTEMPT = 255
PREFERRED_BATCH = 32
U32_MAX = (1 << 32) - 1
U64_MAX = (1 << 64) - 1

ROUTE_HEADER = (
    "N", "bb", "route_status", "preferred_attempt", "canonical_attempt",
    "actual_attempt", "preferred_valid", "fallback", "no_op", "direct",
    "canonical_probe_solves", "preferred_probe_solves",
)
RECOVERY_HEADER = (
    "N", "bb", "arm", "preferred_attempt", "canonical_attempt",
    "actual_attempt", "routed", "preferred_valid", "fallback", "no_op",
    "direct", "physical_solve", "result", "rank_fail", "error",
    "heavy_shortfall", "inactivated", "binary_def", "heavy_gain",
    "block_xors", "block_muladds",
)


@dataclass(frozen=True)
class RoundSpec:
    name: str
    input_max: int
    retain_max: int
    root_indexes: Tuple[int, ...]
    widths: Tuple[int, ...]
    require_control_improvement: bool = False


DEFAULT_ROUNDS = (
    RoundSpec("r1", 256, 64, (0,), (64,)),
    RoundSpec("r2", 64, 16, (1,), (64, 1280)),
    RoundSpec("r3", 16, 8, (2, 3), APPROVED_WIDTHS),
    RoundSpec("r4", 8, 1, (4, 5), APPROVED_WIDTHS, True),
)


@dataclass(frozen=True, order=True)
class Panel:
    root_index: int
    width: int


@dataclass(frozen=True, order=True)
class CellKey:
    K: int
    width: int
    root_index: int
    schedule: str
    loss: str


@dataclass(frozen=True)
class RecoveryMetrics:
    result: int
    rank_fail: int
    error: int
    heavy_shortfall: int
    inactivated: int
    binary_deficit: int
    heavy_gain: int
    block_xors: int
    block_muladds: int

    @property
    def success(self) -> bool:
        return self.result == 0 and self.rank_fail == 0 and self.error == 0

    @property
    def failure(self) -> bool:
        return self.result == 1 and self.rank_fail == 1 and self.error == 0


@dataclass(frozen=True)
class RecoveryRow:
    K: int
    width: int
    arm: str
    preferred_attempt: int
    canonical_attempt: int
    actual_attempt: int
    routed: int
    preferred_valid: int
    fallback: int
    no_op: int
    direct: int
    physical_solve: int
    metrics: RecoveryMetrics


@dataclass(frozen=True)
class RouteRecord:
    K: int
    width: int
    route_status: str
    preferred_attempt: int
    canonical_attempt: int
    actual_attempt: int
    preferred_valid: int
    fallback: int
    no_op: int
    direct: int
    canonical_probe_solves: int
    preferred_probe_solves: int


@dataclass(frozen=True)
class ExecutionResult:
    returncode: int
    stdout: bytes
    stderr: bytes
    start_ns: int
    end_ns: int
    cpu: int


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False,
            ensure_ascii=True,
        ) + "\n"
    ).encode("ascii")


def sealed_record(schema: str, payload: Mapping[str, Any]) -> Dict[str, Any]:
    if (not isinstance(schema, str) or not schema or
            "schema" in payload or "self_sha256_excluding_field" in payload):
        die("sealed record payload uses a reserved field")
    record = {"schema": schema, **dict(payload)}
    record["self_sha256_excluding_field"] = sha256_bytes(
        canonical_json_bytes(record))
    return record


def verify_sealed_record(record: Mapping[str, Any], schema: str) -> None:
    if not isinstance(record, dict) or record.get("schema") != schema:
        die("sealed record schema mismatch")
    digest = record.get("self_sha256_excluding_field")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        die("sealed record has an invalid self hash")
    unsigned = dict(record)
    del unsigned["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json_bytes(unsigned)) != digest:
        die("sealed record self hash mismatch")


def durable_atomic_write(path: Path, data: bytes) -> None:
    common.atomic_write_once_or_same(path, data)
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(str(path.parent), flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_existing_regular_file(path: Path) -> None:
    """Durably order an already-written regular artifact before its seal."""
    descriptor = os.open(
        str(path), os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            die("refusing to seal a nonregular artifact: {}".format(path))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    parent = os.open(
        str(path.parent), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(parent)
    finally:
        os.close(parent)


def strict_ascii_lines(data: bytes, context: str) -> List[str]:
    if not data or not data.endswith(b"\n") or b"\r" in data or b"\0" in data:
        die("{} is not canonical LF text".format(context))
    try:
        return data.decode("ascii").splitlines()
    except UnicodeDecodeError:
        die("{} is not ASCII".format(context))
    return []


def strict_uint(text: str, context: str, maximum: int = U64_MAX) -> int:
    if (not text or not text.isdigit() or
            (len(text) > 1 and text.startswith("0"))):
        die("{} is not a canonical unsigned integer".format(context))
    value = int(text, 10)
    if value > maximum:
        die("{} is outside its integer domain".format(context))
    return value


def strict_attempt(text: str, context: str, allow_control: bool) -> int:
    if allow_control and text == "-1":
        return -1
    value = strict_uint(text, context, MAX_ATTEMPT)
    return value


def strict_bit(text: str, context: str) -> int:
    value = strict_uint(text, context, 1)
    if value not in (0, 1):
        die("{} is not a bit".format(context))
    return value


def is_plain_int(value: Any, minimum: int = 0) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= minimum


def canonical_seed(seed: str) -> str:
    if not isinstance(seed, str) or not re.fullmatch(r"0x[0-9a-f]{16}", seed):
        die("root seed is not a canonical u64")
    return seed


def canonical_context_hash(digest: str) -> str:
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        die("route context SHA256 is not canonical")
    return digest


def canonical_loss(loss: str) -> str:
    if loss not in APPROVED_LOSSES:
        die("loss is outside the frozen development grid")
    return loss


def emitted_loss(loss: str) -> str:
    return format(float(canonical_loss(loss)), ".17g")


def verify_frozen_controller_runtime(result_root: Path) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Delegate the trust decision to the frozen controller implementation.

    The lazy import avoids a module-import cycle when the frozen controller
    imports this helper.  There is intentionally no local fallback or reduced
    verification path: production run/resume/analyze entries must pass the
    controller's public-freeze verification before touching outcomes.
    """
    import wh2_preferred_attempt_search as search

    return search.verify_frozen_controller_runtime(result_root)


def validate_grid_subset(
    values: Sequence[Any], approved: Sequence[Any], context: str,
) -> Tuple[Any, ...]:
    result = tuple(values)
    if (not result or any(value not in approved for value in result) or
            result != tuple(value for value in approved if value in result) or
            len(result) != len(set(result))):
        die("{} is not a nonempty canonical approved subset".format(context))
    return result


def validate_round_spec_value(spec: RoundSpec) -> None:
    if (not isinstance(spec, RoundSpec) or
            not isinstance(spec.name, str) or
            not re.fullmatch(r"r[1-9][0-9]*", spec.name) or
            not isinstance(spec.input_max, int) or
            isinstance(spec.input_max, bool) or spec.input_max <= 0 or
            spec.input_max > MAX_ATTEMPT + 1 or
            not isinstance(spec.retain_max, int) or
            isinstance(spec.retain_max, bool) or spec.retain_max <= 0 or
            spec.retain_max > spec.input_max or
            not isinstance(spec.require_control_improvement, bool) or
            not isinstance(spec.root_indexes, tuple) or
            not spec.root_indexes or
            any(not isinstance(root, int) or isinstance(root, bool) or root < 0
                for root in spec.root_indexes) or
            tuple(sorted(spec.root_indexes)) != spec.root_indexes or
            len(set(spec.root_indexes)) != len(spec.root_indexes) or
            not isinstance(spec.widths, tuple) or
            not spec.widths or
            tuple(width for width in APPROVED_WIDTHS if width in spec.widths) !=
                spec.widths):
        die("development round object is malformed")


def validate_rounds(rounds: Sequence[RoundSpec]) -> None:
    if not rounds:
        die("development round list is empty")
    seen_roots: Set[int] = set()
    previous_retain: Optional[int] = None
    for index, spec in enumerate(rounds):
        validate_round_spec_value(spec)
        if spec.name != "r{}".format(index + 1):
            die("development rounds are not canonically named")
        if previous_retain is not None and spec.input_max != previous_retain:
            die("round input quota does not equal the prior retain quota")
        if seen_roots.intersection(spec.root_indexes):
            die("development roots are repeated or noncanonical")
        seen_roots.update(spec.root_indexes)
        previous_retain = spec.retain_max
    if sum(spec.require_control_improvement for spec in rounds) != 1 or \
            not rounds[-1].require_control_improvement:
        die("only the final round may require strict control improvement")


def declared_panels_through(
    rounds: Sequence[RoundSpec], round_name: str,
) -> Tuple[Panel, ...]:
    validate_rounds(rounds)
    result: List[Panel] = []
    found = False
    for spec in rounds:
        for root_index in spec.root_indexes:
            for width in spec.widths:
                result.append(Panel(root_index, width))
        if spec.name == round_name:
            found = True
            break
    if not found:
        die("unknown development round {}".format(round_name))
    if len(result) != len(set(result)):
        die("declared development panels overlap")
    return tuple(result)


def validate_bins(bins: Mapping[int, Sequence[int]], cohort: Sequence[int]) -> None:
    supplied = tuple(cohort)
    if any(not isinstance(K, int) or isinstance(K, bool) or K < 2
           for K in supplied):
        die("cohort is unsorted, duplicated, or outside the K domain")
    expected = tuple(sorted(supplied))
    if supplied != expected or len(expected) != len(set(expected)):
        die("cohort is unsorted, duplicated, or outside the K domain")
    bin_keys = tuple(bins)
    if (any(not isinstance(index, int) or isinstance(index, bool)
            for index in bin_keys) or
            tuple(sorted(bin_keys)) != tuple(range(len(bins)))):
        die("bin indexes are not contiguous from zero")
    flattened: List[int] = []
    for index in range(len(bins)):
        values = tuple(bins[index])
        if (not values or
                any(not isinstance(K, int) or isinstance(K, bool) or K < 2
                    for K in values) or
                values != tuple(sorted(values)) or len(values) != len(set(values))):
            die("bin {} is empty, unsorted, or duplicated".format(index))
        flattened.extend(values)
    if tuple(sorted(flattened)) != expected or len(flattened) != len(expected):
        die("bins do not exactly partition the frozen cohort")


def validate_survivors(
    survivors: Mapping[int, Sequence[int]], cohort: Sequence[int], input_max: int,
) -> None:
    ordered_cohort = tuple(cohort)
    survivor_keys = tuple(survivors)
    if (not isinstance(input_max, int) or isinstance(input_max, bool) or
            input_max < 0 or input_max > MAX_ATTEMPT + 1 or
            any(not isinstance(K, int) or isinstance(K, bool)
                for K in survivor_keys) or
            ordered_cohort != tuple(sorted(ordered_cohort)) or
            len(ordered_cohort) != len(set(ordered_cohort)) or
            any(not isinstance(K, int) or isinstance(K, bool) or K < 2
                for K in ordered_cohort) or
            tuple(sorted(survivor_keys)) != ordered_cohort):
        die("survivor map does not exactly cover the cohort")
    for K in sorted(survivors):
        attempts = tuple(survivors[K])
        if (len(attempts) > input_max or len(attempts) != len(set(attempts)) or
                any(not isinstance(p, int) or isinstance(p, bool) or
                    p < 0 or p > MAX_ATTEMPT for p in attempts)):
            die("K={} survivor sequence is invalid".format(K))


def initial_survivors(cohort: Sequence[int]) -> Dict[int, Tuple[int, ...]]:
    """Return the frozen R1 input: all byte-valued preferred attempts per K."""
    ordered = tuple(cohort)
    if (ordered != tuple(sorted(ordered)) or len(ordered) != len(set(ordered)) or
            any(not isinstance(K, int) or isinstance(K, bool) or K < 2
                for K in ordered)):
        die("initial survivor cohort is noncanonical")
    attempts = tuple(range(MAX_ATTEMPT + 1))
    return {K: attempts for K in ordered}


@dataclass(frozen=True)
class JobSpec:
    job_id: str
    kind: str
    round_name: str
    bin_index: int
    batch_index: int
    Ks: Tuple[int, ...]
    attempts: Tuple[Tuple[int, Tuple[int, ...]], ...]
    widths: Tuple[int, ...]
    root_index: Optional[int]
    seed: Optional[str]
    schedule: Optional[str]
    loss: Optional[str]
    route_context_sha256: Optional[str]
    route_cache_path: Optional[str] = None
    route_cache_sha256: Optional[str] = None
    route_job_id: Optional[str] = None
    expected_logical_rows: int = 0
    expected_physical_rows: int = 0

    def attempt_map(self) -> Dict[int, Tuple[int, ...]]:
        return {K: values for K, values in self.attempts}

    def validate(self) -> None:
        if (not isinstance(self.job_id, str) or
                not isinstance(self.kind, str) or
                not isinstance(self.round_name, str) or
                not isinstance(self.Ks, tuple) or
                not isinstance(self.attempts, tuple) or
                any(not isinstance(item, tuple) or len(item) != 2 or
                    not isinstance(item[1], tuple) for item in self.attempts) or
                not isinstance(self.widths, tuple)):
            die("job record container fields are malformed")
        if self.kind not in ("control", "route", "candidate"):
            die("job kind is invalid")
        if not re.fullmatch(
                r"(?:control|r[1-9][0-9]*-(?:route|candidate))-[0-9]{5}",
                self.job_id):
            die("job ID is noncanonical")
        if (not isinstance(self.bin_index, int) or
                isinstance(self.bin_index, bool) or self.bin_index < 0 or
                not isinstance(self.batch_index, int) or
                isinstance(self.batch_index, bool) or self.batch_index < 0 or
                not self.Ks or
                any(not isinstance(K, int) or isinstance(K, bool) or K < 2
                    for K in self.Ks) or
                self.Ks != tuple(sorted(self.Ks)) or
                len(self.Ks) != len(set(self.Ks))):
            die("job K/bin/batch identity is invalid")
        if (not self.widths or
                any(not isinstance(width, int) or isinstance(width, bool)
                    for width in self.widths) or
                tuple(width for width in APPROVED_WIDTHS
                                     if width in self.widths) != self.widths):
            die("job width list is invalid")
        if (not isinstance(self.expected_logical_rows, int) or
                isinstance(self.expected_logical_rows, bool) or
                not isinstance(self.expected_physical_rows, int) or
                isinstance(self.expected_physical_rows, bool)):
            die("job row counts are not canonical integers")
        mapping = self.attempt_map()
        if (len(mapping) != len(self.attempts) or
                (self.kind != "control" and
                 tuple(K for K, _values in self.attempts) != self.Ks)):
            die("job attempt map is repeated, reordered, or incomplete")
        if self.kind == "control":
            if (self.job_id.split("-")[0] != "control" or
                    self.round_name != "control" or self.attempts or
                    self.batch_index != 0 or self.root_index is None or
                    not isinstance(self.root_index, int) or
                    isinstance(self.root_index, bool) or self.root_index < 0 or
                    self.seed is None or self.schedule not in APPROVED_SCHEDULES or
                    self.loss not in APPROVED_LOSSES or
                    self.route_context_sha256 is not None or
                    self.route_cache_path is not None or
                    self.route_cache_sha256 is not None or
                    self.route_job_id is not None):
                die("control job has candidate-only fields")
            canonical_seed(self.seed)
            expected = len(self.Ks) * len(self.widths)
            if (self.expected_logical_rows != expected or
                    self.expected_physical_rows != expected):
                die("control job row arithmetic mismatch")
            return
        if (not re.fullmatch(r"r[1-9][0-9]*", self.round_name) or
                not self.job_id.startswith(
                    "{}-{}-".format(self.round_name, self.kind))):
            die("preferred job ID/round/kind binding is invalid")
        if set(mapping) != set(self.Ks):
            die("preferred job map does not exactly cover its K list")
        for K, values in self.attempts:
            if not values or len(values) > PREFERRED_BATCH or len(values) != len(set(values)):
                die("preferred job attempt batch is empty, duplicated, or oversized")
            if any(not isinstance(p, int) or isinstance(p, bool) or
                   p < 0 or p > MAX_ATTEMPT for p in values):
                die("preferred job attempt is outside the byte domain")
        if self.route_context_sha256 is None:
            die("preferred job lacks its route context")
        canonical_context_hash(self.route_context_sha256)
        logical = sum(len(values) for values in mapping.values()) * len(self.widths)
        if self.expected_logical_rows != logical:
            die("preferred job logical-row arithmetic mismatch")
        if self.kind == "route":
            if (self.root_index is not None or self.seed is not None or
                    self.schedule is not None or self.loss is not None or
                    self.route_cache_path is not None or
                    self.route_cache_sha256 is not None or
                    self.route_job_id is not None or
                    self.expected_physical_rows != 0):
                die("route job has recovery-only fields")
        else:
            if (self.root_index is None or self.seed is None or
                    not isinstance(self.root_index, int) or
                    isinstance(self.root_index, bool) or self.root_index < 0 or
                    self.schedule not in APPROVED_SCHEDULES or
                    self.loss not in APPROVED_LOSSES or
                    not isinstance(self.route_cache_path, str) or
                    self.route_cache_sha256 is None or
                    not isinstance(self.route_job_id, str) or
                    not re.fullmatch(
                        r"{}-route-[0-9]{{5}}".format(self.round_name),
                        self.route_job_id) or
                    not isinstance(self.route_cache_sha256, str) or
                    not re.fullmatch(r"[0-9a-f]{64}", self.route_cache_sha256) or
                    self.expected_physical_rows < 0 or
                    self.expected_physical_rows > logical):
                die("candidate job has invalid recovery/cache fields")
            canonical_seed(self.seed)
            route_path = Path(self.route_cache_path)
            # Candidate job ordinals expand each route batch over strata, so
            # the route job ID is not derived from the candidate ordinal.
            expected_cache = (
                "jobs/{}/route/stdout/{}.csv".format(
                    self.round_name, self.route_job_id))
            if (not self.route_cache_path or "\\" in self.route_cache_path or
                    route_path.is_absolute() or ".." in route_path.parts or
                    route_path.as_posix() != self.route_cache_path or
                    self.route_cache_path != expected_cache):
                die("candidate route-cache path is unsafe")

    def to_record(self) -> Dict[str, Any]:
        self.validate()
        return {
            "job_id": self.job_id,
            "kind": self.kind,
            "round": self.round_name,
            "bin": self.bin_index,
            "batch": self.batch_index,
            "K": list(self.Ks),
            "attempts": [
                {"K": K, "p": list(values)} for K, values in self.attempts
            ],
            "widths": list(self.widths),
            "root_index": self.root_index,
            "seed": self.seed,
            "schedule": self.schedule,
            "loss": self.loss,
            "route_context_sha256": self.route_context_sha256,
            "route_cache_path": self.route_cache_path,
            "route_cache_sha256": self.route_cache_sha256,
            "route_job_id": self.route_job_id,
            "expected_logical_rows": self.expected_logical_rows,
            "expected_physical_rows": self.expected_physical_rows,
        }

    def command(self, binary: Path, result_root: Path) -> Tuple[str, ...]:
        self.validate()
        command = [
            str(binary), "preferredattempt", "--mode", self.kind,
            "--N", ",".join(str(K) for K in self.Ks),
            "--bb-list", ",".join(str(width) for width in self.widths),
        ]
        if self.kind != "control":
            mapping = self.attempt_map()
            command.extend((
                "--preferred-map",
                "|".join(
                    "{}@{}={}".format(
                        K, width, ",".join(str(p) for p in mapping[K]))
                    for K in self.Ks for width in self.widths
                ),
                "--route-context-sha256", str(self.route_context_sha256),
            ))
        if self.kind == "candidate":
            command.extend((
                "--route-cache", str(result_root / str(self.route_cache_path)),
                "--route-cache-sha256", str(self.route_cache_sha256),
            ))
        if self.kind in ("control", "candidate"):
            command.extend((
                "--loss", str(self.loss), "--seed", str(self.seed),
                "--schedule", str(self.schedule),
            ))
        return tuple(command)


@dataclass(frozen=True)
class JobLedger:
    kind: str
    phase: str
    parent_sha256: str
    jobs: Tuple[JobSpec, ...]

    def record(self) -> Dict[str, Any]:
        if (not isinstance(self.kind, str) or not isinstance(self.phase, str) or
                not isinstance(self.jobs, tuple)):
            die("job ledger container fields are malformed")
        if self.kind not in ("control", "route", "candidate"):
            die("job ledger kind is invalid")
        if ((self.kind == "control") != (self.phase == "control") or
                (self.kind != "control" and
                 not re.fullmatch(r"r[1-9][0-9]*", self.phase))):
            die("job ledger phase/kind identity is invalid")
        if self.kind == "control" and not self.jobs:
            die("control ledger may not be empty")
        if (not isinstance(self.parent_sha256, str) or
                not re.fullmatch(r"[0-9a-f]{64}", self.parent_sha256)):
            die("job ledger parent hash is invalid")
        logical_keys: Set[Tuple[Any, ...]] = set()
        candidate_route_sources: Dict[Tuple[int, int, int], Tuple[Any, ...]] = {}
        for index, job in enumerate(self.jobs):
            job.validate()
            if job.kind != self.kind or job.round_name != self.phase:
                die("job ledger mixes job kinds")
            expected_id = (
                "control-{:05d}".format(index) if self.kind == "control" else
                "{}-{}-{:05d}".format(self.phase, self.kind, index)
            )
            if job.job_id != expected_id:
                die("job ledger IDs are not contiguous and canonical")
            if job.kind == "control":
                keys = (
                    (K, width, job.root_index, job.schedule, job.loss)
                    for K in job.Ks for width in job.widths
                )
            else:
                mapping = job.attempt_map()
                keys = (
                    (K, width, p) if job.kind == "route" else
                    (K, width, p, job.root_index, job.schedule, job.loss)
                    for K in job.Ks for width in job.widths
                    for p in mapping[K]
                )
                if job.kind == "candidate":
                    source = (
                        job.route_job_id, job.route_cache_path,
                        job.route_cache_sha256, job.route_context_sha256)
                    for K in job.Ks:
                        for width in job.widths:
                            for p in mapping[K]:
                                route_key = (K, width, p)
                                old_source = candidate_route_sources.setdefault(
                                    route_key, source)
                                if old_source != source:
                                    die("candidate ledger changes a route-cache source")
            for key in keys:
                if key in logical_keys:
                    die("job ledger repeats a logical campaign cell")
                logical_keys.add(key)
        logical = sum(job.expected_logical_rows for job in self.jobs)
        physical = sum(job.expected_physical_rows for job in self.jobs)
        return sealed_record(
            SCHEMA_PREFIX + ".job_ledger.v1",
            {
                "kind": self.kind,
                "phase": self.phase,
                "parent_sha256": self.parent_sha256,
                "jobs": [job.to_record() for job in self.jobs],
                "job_count": len(self.jobs),
                "expected_logical_rows": logical,
                "expected_physical_rows": physical,
            },
        )

    def bytes(self) -> bytes:
        return canonical_json_bytes(self.record())

    def sha256(self) -> str:
        return sha256_bytes(self.bytes())


def durable_write_once_or_same(path: Path, data: bytes) -> None:
    """Publish immutable campaign bytes, accepting only an identical resume."""
    try:
        common.atomic_write_once_or_same(path, data)
    except common.CampaignError as exc:
        die("immutable campaign artifact publication failed: {}: {}".format(
            path, exc))


def write_hashed_artifact(path: Path, data: bytes) -> str:
    digest = sha256_bytes(data)
    durable_write_once_or_same(path, data)
    sidecar = (digest + "  " + path.name + "\n").encode("ascii")
    durable_write_once_or_same(path.with_suffix(path.suffix + ".sha256"), sidecar)
    return digest


def _write_job_ledger(result_root: Path, ledger: JobLedger) -> str:
    path = result_root / "ledgers" / (
        "{}-{}.json".format(ledger.phase, ledger.kind))
    return write_hashed_artifact(path, ledger.bytes())


def write_job_ledger(result_root: Path, ledger: JobLedger) -> str:
    """Verify the frozen runtime and publish a canonical immutable job ledger."""
    result_root = result_root.resolve()
    verify_frozen_controller_runtime(result_root)
    return _write_job_ledger(result_root, ledger)


@dataclass(frozen=True)
class RouteArtifact:
    job_id: str
    relative_path: str
    sha256: str
    records: Tuple[RouteRecord, ...]

    @property
    def physical_rows(self) -> int:
        return sum(record.direct for record in self.records)


def normalized_attempts(
    survivors: Mapping[int, Sequence[int]], Ks: Sequence[int], begin: int,
    batch_size: int,
) -> Tuple[Tuple[int, Tuple[int, ...]], ...]:
    result = []
    for K in Ks:
        values = tuple(survivors[K][begin:begin + batch_size])
        if values:
            result.append((K, values))
    return tuple(result)


def build_control_ledger(
    bins: Mapping[int, Sequence[int]],
    cohort: Sequence[int],
    roots: Sequence[str],
    parent_sha256: str,
    schedules: Sequence[str] = APPROVED_SCHEDULES,
    losses: Sequence[str] = APPROVED_LOSSES,
    widths: Sequence[int] = APPROVED_WIDTHS,
) -> JobLedger:
    validate_bins(bins, cohort)
    if not roots:
        die("control root list is empty")
    for seed in roots:
        canonical_seed(seed)
    schedules_tuple = validate_grid_subset(
        schedules, APPROVED_SCHEDULES, "control schedules")
    losses_tuple = validate_grid_subset(
        losses, APPROVED_LOSSES, "control losses")
    widths_tuple = validate_grid_subset(
        widths, APPROVED_WIDTHS, "control widths")
    jobs: List[JobSpec] = []
    for bin_index in range(len(bins)):
        Ks = tuple(bins[bin_index])
        for root_index, seed in enumerate(roots):
            for schedule in schedules_tuple:
                for loss in losses_tuple:
                    row_count = len(Ks) * len(widths_tuple)
                    jobs.append(JobSpec(
                        "control-{:05d}".format(len(jobs)), "control", "control",
                        bin_index, 0, Ks, (), widths_tuple, root_index, seed,
                        schedule, loss, None,
                        expected_logical_rows=row_count,
                        expected_physical_rows=row_count,
                    ))
    return JobLedger("control", "control", parent_sha256, tuple(jobs))


def build_route_ledger(
    spec: RoundSpec,
    bins: Mapping[int, Sequence[int]],
    cohort: Sequence[int],
    survivors: Mapping[int, Sequence[int]],
    parent_sha256: str,
    route_context_sha256: str,
    batch_size: int = PREFERRED_BATCH,
) -> JobLedger:
    validate_round_spec_value(spec)
    validate_bins(bins, cohort)
    validate_survivors(survivors, cohort, spec.input_max)
    canonical_context_hash(route_context_sha256)
    if (not isinstance(batch_size, int) or isinstance(batch_size, bool) or
            batch_size <= 0 or batch_size > PREFERRED_BATCH):
        die("preferred batch size is outside 1..32")
    jobs: List[JobSpec] = []
    for bin_index in range(len(bins)):
        bin_Ks = tuple(bins[bin_index])
        maximum = max(len(survivors[K]) for K in bin_Ks)
        for begin in range(0, maximum, batch_size):
            mapping = normalized_attempts(
                survivors, bin_Ks, begin, batch_size)
            if not mapping:
                continue
            Ks = tuple(K for K, _values in mapping)
            logical = sum(len(values) for _K, values in mapping) * len(spec.widths)
            jobs.append(JobSpec(
                "{}-route-{:05d}".format(spec.name, len(jobs)),
                "route", spec.name, bin_index, begin // batch_size, Ks,
                mapping, spec.widths, None, None, None, None,
                route_context_sha256, expected_logical_rows=logical,
            ))
    return JobLedger("route", spec.name, parent_sha256, tuple(jobs))


def build_candidate_ledger(
    spec: RoundSpec,
    route_ledger: JobLedger,
    route_artifacts: Mapping[str, RouteArtifact],
    roots: Sequence[str],
    parent_sha256: str,
    schedules: Sequence[str] = APPROVED_SCHEDULES,
    losses: Sequence[str] = APPROVED_LOSSES,
) -> JobLedger:
    validate_round_spec_value(spec)
    if route_ledger.kind != "route" or route_ledger.phase != spec.name:
        die("candidate ledger received another round's route ledger")
    route_ledger.record()
    if (not isinstance(route_artifacts, Mapping) or
            set(route_artifacts) != {job.job_id for job in route_ledger.jobs}):
        die("route artifact set does not exactly cover the route ledger")
    if any(index < 0 or index >= len(roots) for index in spec.root_indexes):
        die("round root index is outside the frozen root ledger")
    schedules_tuple = validate_grid_subset(
        schedules, APPROVED_SCHEDULES, "candidate schedules")
    losses_tuple = validate_grid_subset(
        losses, APPROVED_LOSSES, "candidate losses")
    jobs: List[JobSpec] = []
    canonical_by_key: Dict[Tuple[int, int], int] = {}
    for route_job in route_ledger.jobs:
        route_job.validate()
        if route_job.widths != spec.widths:
            die("route ledger widths do not match the declared round")
        artifact = route_artifacts[route_job.job_id]
        if (not isinstance(artifact, RouteArtifact) or
                not isinstance(artifact.sha256, str) or
                not isinstance(artifact.relative_path, str) or
                not isinstance(artifact.records, tuple) or
                artifact.job_id != route_job.job_id or
                not re.fullmatch(r"[0-9a-f]{64}", artifact.sha256)):
            die("route artifact identity is invalid")
        expected_relative = (
            "jobs/{}/route/stdout/{}.csv".format(spec.name, route_job.job_id))
        if artifact.relative_path != expected_relative:
            die("route artifact path is not the canonical job output")
        expected_sequence = [
            (K, width, p)
            for K, values in route_job.attempts
            for width in route_job.widths for p in values
        ]
        actual_sequence = [
            (record.K, record.width, record.preferred_attempt)
            for record in artifact.records if record.route_status == "preferred"
        ]
        if (actual_sequence != expected_sequence or
                len(artifact.records) != len(expected_sequence)):
            die("route artifact does not exactly cover its job")
        charge_seen: Set[Tuple[int, int]] = set()
        for record in artifact.records:
            validate_route_record_value(record)
            key = (record.K, record.width)
            previous = canonical_by_key.setdefault(key, record.canonical_attempt)
            if previous != record.canonical_attempt:
                die("canonical attempt changed between route batches")
            expected_charge = (
                record.canonical_attempt + 1 if key not in charge_seen else 0)
            if record.canonical_probe_solves != expected_charge:
                die("route artifact canonical probe charge is noncanonical")
            charge_seen.add(key)
        for root_index in spec.root_indexes:
            seed = canonical_seed(roots[root_index])
            for schedule in schedules_tuple:
                for loss in losses_tuple:
                    jobs.append(JobSpec(
                        "{}-candidate-{:05d}".format(spec.name, len(jobs)),
                        "candidate", spec.name, route_job.bin_index,
                        route_job.batch_index, route_job.Ks, route_job.attempts,
                        route_job.widths, root_index, seed, schedule, loss,
                        route_job.route_context_sha256,
                        artifact.relative_path, artifact.sha256,
                        route_job.job_id, route_job.expected_logical_rows,
                        artifact.physical_rows,
                    ))
    return JobLedger("candidate", spec.name, parent_sha256, tuple(jobs))


def expected_route_preamble(context_sha256: str) -> str:
    return (
        "# preferredattempt-route-manifest: schema=v1 "
        "policy=h12-q0-adaptive canonical=ascending-first-valid-v1 "
        "max_attempt=255 context_sha256={}".format(
            canonical_context_hash(context_sha256))
    )


def parse_route_record(fields: Sequence[str], context: str) -> RouteRecord:
    if len(fields) != len(ROUTE_HEADER) or any(field == "" for field in fields):
        die("{} has a malformed route row".format(context))
    K = strict_uint(fields[0], context + ":N", U32_MAX)
    width = strict_uint(fields[1], context + ":bb", U32_MAX)
    status = fields[2]
    if status not in ("preferred", "control"):
        die("{} has an invalid route status".format(context))
    preferred = strict_attempt(fields[3], context + ":p", status == "control")
    canonical = strict_uint(fields[4], context + ":a0", MAX_ATTEMPT)
    actual = strict_uint(fields[5], context + ":actual", MAX_ATTEMPT)
    valid = strict_bit(fields[6], context + ":valid")
    fallback = strict_bit(fields[7], context + ":fallback")
    no_op = strict_bit(fields[8], context + ":no_op")
    direct = strict_bit(fields[9], context + ":direct")
    canonical_probes = strict_uint(
        fields[10], context + ":canonical_probe_solves", 256)
    preferred_probes = strict_uint(
        fields[11], context + ":preferred_probe_solves", 1)
    if status == "control":
        if (preferred != -1 or actual != canonical or valid != 1 or
                fallback != 0 or no_op != 1 or direct != 0 or
                preferred_probes != 0):
            die("{} control route invariants failed".format(context))
    else:
        expected_noop = int(preferred == canonical)
        expected_fallback = int(not valid)
        expected_direct = int(bool(valid) and not bool(expected_noop))
        expected_actual = preferred if expected_direct else canonical
        expected_probe = int(preferred > canonical)
        if (no_op != expected_noop or fallback != expected_fallback or
                direct != expected_direct or actual != expected_actual or
                preferred_probes != expected_probe or
                (preferred < canonical and valid != 0) or
                (preferred == canonical and valid != 1)):
            die("{} preferred route invariants failed".format(context))
    return RouteRecord(
        K, width, status, preferred, canonical, actual, valid, fallback,
        no_op, direct, canonical_probes, preferred_probes)


def validate_route_record_value(record: RouteRecord) -> None:
    fields = (
        record.K, record.width, record.route_status,
        record.preferred_attempt, record.canonical_attempt,
        record.actual_attempt, record.preferred_valid, record.fallback,
        record.no_op, record.direct, record.canonical_probe_solves,
        record.preferred_probe_solves,
    )
    reparsed = parse_route_record(
        tuple(str(value) for value in fields), "route artifact")
    if reparsed != record:
        die("route artifact record is not canonical")


def parse_route_output(job: JobSpec, data: bytes) -> Tuple[RouteRecord, ...]:
    job.validate()
    if job.kind != "route":
        die("route parser received a non-route job")
    lines = strict_ascii_lines(data, job.job_id + " route output")
    if len(lines) < 2 or lines[0] != expected_route_preamble(
            str(job.route_context_sha256)) or lines[1] != ",".join(ROUTE_HEADER):
        die("{} route preamble/header mismatch".format(job.job_id))
    mapping = job.attempt_map()
    expected = [
        (K, width, p)
        for K in job.Ks for width in job.widths for p in mapping[K]
    ]
    if len(lines) != len(expected) + 2:
        die("{} route row cardinality mismatch".format(job.job_id))
    records: List[RouteRecord] = []
    charge_seen: Set[Tuple[int, int]] = set()
    canonical_by_key: Dict[Tuple[int, int], int] = {}
    for index, ((K, width, p), line) in enumerate(zip(expected, lines[2:]), 1):
        record = parse_route_record(
            line.split(","), "{}:row{}".format(job.job_id, index))
        if (record.K, record.width, record.preferred_attempt) != (K, width, p) or \
                record.route_status != "preferred":
            die("{} route row ordering/key mismatch".format(job.job_id))
        key = (K, width)
        old = canonical_by_key.setdefault(key, record.canonical_attempt)
        if old != record.canonical_attempt:
            die("{} route canonical attempt changed within a key".format(job.job_id))
        expected_charge = record.canonical_attempt + 1 if key not in charge_seen else 0
        if record.canonical_probe_solves != expected_charge:
            die("{} canonical probe accounting mismatch".format(job.job_id))
        charge_seen.add(key)
        records.append(record)
    if len(records) != job.expected_logical_rows:
        die("{} route ledger/output arithmetic mismatch".format(job.job_id))
    return tuple(records)


def parse_recovery_metrics(fields: Sequence[str], context: str) -> RecoveryMetrics:
    result = strict_uint(fields[12], context + ":result", U32_MAX)
    rank_fail = strict_bit(fields[13], context + ":rank_fail")
    error = strict_bit(fields[14], context + ":error")
    heavy_shortfall = strict_bit(fields[15], context + ":heavy_shortfall")
    inactivated = strict_uint(fields[16], context + ":inactivated", U32_MAX)
    binary_deficit = strict_uint(fields[17], context + ":binary_def", U32_MAX)
    heavy_gain = strict_uint(fields[18], context + ":heavy_gain", U32_MAX)
    block_xors = strict_uint(fields[19], context + ":block_xors", U64_MAX)
    block_muladds = strict_uint(fields[20], context + ":block_muladds", U64_MAX)
    expected_rank_fail = int(result == 1)
    expected_error = int(result not in (0, 1))
    if (rank_fail != expected_rank_fail or error != expected_error or
            binary_deficit > inactivated or
            (heavy_shortfall and not rank_fail)):
        die("{} recovery metric invariants failed".format(context))
    return RecoveryMetrics(
        result, rank_fail, error, heavy_shortfall, inactivated,
        binary_deficit, heavy_gain, block_xors, block_muladds)


def validate_recovery_metrics_value(
    value: RecoveryMetrics, context: str,
) -> None:
    if not isinstance(value, RecoveryMetrics):
        die("{} is not a recovery-metric record".format(context))
    raw = (
        value.result, value.rank_fail, value.error, value.heavy_shortfall,
        value.inactivated, value.binary_deficit, value.heavy_gain,
        value.block_xors, value.block_muladds,
    )
    if any(not isinstance(item, int) or isinstance(item, bool) for item in raw):
        die("{} has noncanonical recovery-metric integers".format(context))
    fields = ["0"] * 12 + [str(item) for item in raw]
    if parse_recovery_metrics(fields, context) != value:
        die("{} recovery metrics do not round-trip".format(context))


def parse_recovery_row(fields: Sequence[str], context: str) -> RecoveryRow:
    if len(fields) != len(RECOVERY_HEADER) or any(field == "" for field in fields):
        die("{} has a malformed recovery row".format(context))
    K = strict_uint(fields[0], context + ":N", U32_MAX)
    width = strict_uint(fields[1], context + ":bb", U32_MAX)
    arm = fields[2]
    if arm not in ("control", "candidate"):
        die("{} has an invalid recovery arm".format(context))
    preferred = strict_attempt(fields[3], context + ":p", arm == "control")
    canonical = strict_uint(fields[4], context + ":a0", MAX_ATTEMPT)
    actual = strict_uint(fields[5], context + ":actual", MAX_ATTEMPT)
    routed = strict_bit(fields[6], context + ":routed")
    valid = strict_bit(fields[7], context + ":valid")
    fallback = strict_bit(fields[8], context + ":fallback")
    no_op = strict_bit(fields[9], context + ":no_op")
    direct = strict_bit(fields[10], context + ":direct")
    physical = strict_bit(fields[11], context + ":physical")
    metrics = parse_recovery_metrics(fields, context)
    return RecoveryRow(
        K, width, arm, preferred, canonical, actual, routed, valid, fallback,
        no_op, direct, physical, metrics)


def expected_recovery_preamble(job: JobSpec) -> str:
    if job.kind not in ("control", "candidate"):
        die("recovery preamble requested for a route job")
    return (
        "# preferredattempt: schema=v2 policy=h12-q0-adaptive "
        "canonical=ascending-first-valid-v1 mode={} loss={} seed={} "
        "schedule={} preferred_batch_max=32 route_cache_sha256={} "
        "route_context_sha256={} probe_route=0 "
        "physical_solve_accounting=explicit "
        "systematic_probe_accounting=explicit"
    ).format(
        job.kind, emitted_loss(str(job.loss)), hex(int(str(job.seed), 0)),
        job.schedule,
        job.route_cache_sha256 if job.kind == "candidate" else "none",
        job.route_context_sha256 if job.kind == "candidate" else "none",
    )


ALIAS_PATTERN = re.compile(
    r"^# preferred_candidate_alias: N=(0|[1-9][0-9]*) "
    r"bb=(0|[1-9][0-9]*) p=(0|[1-9][0-9]*) "
    r"a0=(0|[1-9][0-9]*) actual=(0|[1-9][0-9]*) "
    r"valid=([01]) fallback=([01]) "
    r"no_op=([01]) direct=0 physical_solve=0$"
)
PROBE_PATTERN = re.compile(
    r"^# preferred_probe_accounting: N=(0|[1-9][0-9]*) "
    r"bb=(0|[1-9][0-9]*) "
    r"canonical_probe_solves=(0|[1-9][0-9]*) "
    r"preferred_probe_solves=0$"
)


@dataclass(frozen=True)
class ParsedRecoveryOutput:
    rows: Tuple[RecoveryRow, ...]
    alias_rows: int
    physical_rows: int
    canonical_probe_solves: int


def validate_control_row(row: RecoveryRow, K: int, width: int) -> None:
    if ((row.K, row.width, row.arm, row.preferred_attempt) !=
            (K, width, "control", -1) or
            row.actual_attempt != row.canonical_attempt or
            (row.routed, row.preferred_valid, row.fallback, row.no_op,
             row.direct, row.physical_solve) != (0, 1, 0, 0, 0, 1)):
        die("control recovery row invariant/key mismatch")


def validate_candidate_row(row: RecoveryRow, route: RouteRecord) -> None:
    if ((row.K, row.width, row.arm, row.preferred_attempt) !=
            (route.K, route.width, "candidate", route.preferred_attempt) or
            (row.canonical_attempt, row.actual_attempt, row.routed,
             row.preferred_valid, row.fallback, row.no_op, row.direct,
             row.physical_solve) !=
            (route.canonical_attempt, route.actual_attempt, 1,
             route.preferred_valid, route.fallback, route.no_op,
             route.direct, 1) or not route.direct):
        die("candidate recovery row disagrees with its route cache")


def parse_alias_comment(line: str, route: RouteRecord) -> None:
    match = ALIAS_PATTERN.fullmatch(line)
    if match is None:
        die("candidate output has a malformed alias record")
    values = tuple(int(value, 10) for value in match.groups())
    expected = (
        route.K, route.width, route.preferred_attempt,
        route.canonical_attempt, route.actual_attempt,
        route.preferred_valid, route.fallback, route.no_op,
    )
    if values != expected or route.direct:
        die("candidate alias record disagrees with its route cache")


def parse_recovery_output(
    job: JobSpec,
    data: bytes,
    route_records: Sequence[RouteRecord] = (),
) -> ParsedRecoveryOutput:
    job.validate()
    if job.kind not in ("control", "candidate"):
        die("recovery parser received a route job")
    lines = strict_ascii_lines(data, job.job_id + " recovery output")
    if len(lines) < 2 or lines[0] != expected_recovery_preamble(job) or \
            lines[1] != ",".join(RECOVERY_HEADER):
        die("{} recovery preamble/header mismatch".format(job.job_id))
    cursor = 2
    rows: List[RecoveryRow] = []
    alias_rows = 0
    canonical_probes = 0
    if job.kind == "control":
        for K in job.Ks:
            for width in job.widths:
                if cursor >= len(lines):
                    die("{} control output ended early".format(job.job_id))
                match = PROBE_PATTERN.fullmatch(lines[cursor])
                if match is None:
                    die("{} lacks canonical probe accounting".format(job.job_id))
                values = tuple(int(value, 10) for value in match.groups())
                if values[0:2] != (K, width) or values[2] <= 0 or values[2] > 256:
                    die("{} control probe accounting mismatch".format(job.job_id))
                canonical_probes += values[2]
                cursor += 1
                if cursor >= len(lines):
                    die("{} control row is missing".format(job.job_id))
                row = parse_recovery_row(
                    lines[cursor].split(","), "{}:{}".format(job.job_id, cursor + 1))
                validate_control_row(row, K, width)
                if row.canonical_attempt + 1 != values[2]:
                    die("{} control canonical probes do not match a0".format(job.job_id))
                rows.append(row)
                cursor += 1
    else:
        expected_routes = [
            record for record in route_records
            if record.route_status == "preferred"
        ]
        mapping = job.attempt_map()
        expected_keys = [
            (K, width, p)
            for K in job.Ks for width in job.widths for p in mapping[K]
        ]
        route_by_key = {
            (record.K, record.width, record.preferred_attempt): record
            for record in expected_routes
        }
        if len(route_by_key) != len(expected_routes) or set(route_by_key) != set(expected_keys):
            die("candidate route cache coverage mismatch")
        for key in expected_keys:
            if cursor >= len(lines):
                die("{} candidate output ended early".format(job.job_id))
            route = route_by_key[key]
            if route.direct:
                row = parse_recovery_row(
                    lines[cursor].split(","), "{}:{}".format(job.job_id, cursor + 1))
                validate_candidate_row(row, route)
                rows.append(row)
            else:
                parse_alias_comment(lines[cursor], route)
                alias_rows += 1
            cursor += 1
    if cursor != len(lines):
        die("{} output contains unexpected trailing records".format(job.job_id))
    physical = sum(row.physical_solve for row in rows)
    if (job.kind == "control" and
            (len(rows) != job.expected_logical_rows or
             physical != job.expected_physical_rows)) or \
            (job.kind == "candidate" and
             (physical != job.expected_physical_rows or
              alias_rows + physical != job.expected_logical_rows)):
        die("{} logical/physical output accounting mismatch".format(job.job_id))
    return ParsedRecoveryOutput(tuple(rows), alias_rows, physical, canonical_probes)


@dataclass(frozen=True)
class CompletedJob:
    job: JobSpec
    receipt: Dict[str, Any]
    route_records: Tuple[RouteRecord, ...] = ()
    recovery: Optional[ParsedRecoveryOutput] = None


def job_paths(result_root: Path, job: JobSpec) -> Dict[str, Path]:
    base = result_root / "jobs" / job.round_name / job.kind
    return {
        "stdout": base / "stdout" / (job.job_id + ".csv"),
        "stderr": base / "stderr" / (job.job_id + ".txt"),
        "receipt": base / "receipts" / (job.job_id + ".json"),
        "lock": base / "locks" / (job.job_id + ".lock"),
    }


def path_present(path: Path) -> bool:
    """Return true for every directory entry, including a broken symlink."""
    return os.path.lexists(str(path))


def discard_partial_output_pair(paths: Mapping[str, Path], job_id: str) -> None:
    """Remove any unreceipted outputs; bytes alone prove no execution origin."""
    existing = [
        paths[name] for name in ("stdout", "stderr")
        if path_present(paths[name])
    ]
    if not existing:
        die("{} uncommitted-output cleanup found no outputs".format(job_id))
    for path in existing:
        common.stable_bytes(path)  # Reject symlinks/nonregular/changing inputs.
    for path in existing:
        descriptor = common.open_durable_directory(path.parent)
        try:
            os.unlink(path.name, dir_fd=descriptor)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


@contextmanager
def job_lock(path: Path) -> Iterator[None]:
    directory_fd = common.open_durable_directory(path.parent, create=True)
    try:
        descriptor = os.open(
            path.name,
            os.O_CREAT | os.O_RDWR | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0),
            0o600, dir_fd=directory_fd)
    except BaseException:
        os.close(directory_fd)
        raise
    try:
        metadata = os.fstat(descriptor)
        named = os.stat(
            path.name, dir_fd=directory_fd, follow_symlinks=False)
        if (not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1 or
                (metadata.st_dev, metadata.st_ino) !=
                (named.st_dev, named.st_ino)):
            die("job lock is not a unique regular file")
        os.fsync(directory_fd)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            die("another worker owns {}".format(path.name))
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)
            os.close(directory_fd)


def route_job_from_candidate(job: JobSpec) -> JobSpec:
    if job.kind != "candidate" or job.route_job_id is None:
        die("cannot reconstruct a route job from this job")
    return JobSpec(
        job.route_job_id, "route", job.round_name, job.bin_index,
        job.batch_index, job.Ks, job.attempts, job.widths,
        None, None, None, None, job.route_context_sha256,
        expected_logical_rows=job.expected_logical_rows,
    )


def stable_result_file(result_root: Path, relative: str, context: str) -> bytes:
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        die("{} path is unsafe".format(context))
    try:
        root = result_root.resolve(strict=True)
    except OSError as error:
        die("{} result root is unavailable: {}".format(context, error))
    lexical = root / relative_path
    current = root
    for part in relative_path.parts:
        current = current / part
        if current.is_symlink():
            die("{} path traverses a symlink".format(context))
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as error:
        die("{} file is unavailable: {}".format(context, error))
    if resolved == root or root not in resolved.parents:
        die("{} path escapes the result root".format(context))
    return common.stable_bytes(resolved)


def candidate_route_bytes(job: JobSpec, result_root: Path) -> bytes:
    if job.kind != "candidate" or job.route_cache_path is None:
        die("route-cache reader received a noncandidate job")
    route_bytes = stable_result_file(
        result_root, job.route_cache_path, job.job_id + " route cache")
    if sha256_bytes(route_bytes) != job.route_cache_sha256:
        die("candidate route cache changed before output validation")
    return route_bytes


def parse_job_bytes(
    job: JobSpec, stdout: bytes, result_root: Path,
) -> Tuple[Tuple[RouteRecord, ...], Optional[ParsedRecoveryOutput]]:
    if job.kind == "route":
        return parse_route_output(job, stdout), None
    if job.kind == "control":
        return (), parse_recovery_output(job, stdout)
    route_bytes = candidate_route_bytes(job, result_root)
    route_records = parse_route_output(route_job_from_candidate(job), route_bytes)
    return (), parse_recovery_output(job, stdout, route_records)


def validate_execution_result(execution: Any) -> ExecutionResult:
    if (not isinstance(execution, ExecutionResult) or
            not isinstance(execution.returncode, int) or
            isinstance(execution.returncode, bool) or
            not isinstance(execution.stdout, bytes) or
            not isinstance(execution.stderr, bytes) or
            not isinstance(execution.cpu, int) or
            isinstance(execution.cpu, bool) or execution.cpu < 0 or
            not isinstance(execution.start_ns, int) or
            isinstance(execution.start_ns, bool) or execution.start_ns < 0 or
            not isinstance(execution.end_ns, int) or
            isinstance(execution.end_ns, bool) or execution.end_ns < 0):
        die("job executor returned a malformed result")
    return execution


def receipt_for(
    job: JobSpec,
    command: Sequence[str],
    result_root: Path,
    paths: Mapping[str, Path],
    execution: ExecutionResult,
    route_records: Sequence[RouteRecord],
    recovery: Optional[ParsedRecoveryOutput],
) -> Dict[str, Any]:
    validate_execution_result(execution)
    if execution.returncode != 0 or execution.stderr:
        die("cannot seal an unsuccessful job execution")
    stdout = common.stable_bytes(paths["stdout"])
    stderr = common.stable_bytes(paths["stderr"])
    if stdout != execution.stdout or stderr != execution.stderr:
        die("{} published output changed before receipt sealing".format(
            job.job_id))
    if stderr:
        die("{} has nonempty stderr".format(job.job_id))
    logical = job.expected_logical_rows
    if job.kind == "route":
        physical = 0
        aliases = sum(1 - record.direct for record in route_records)
        route_rows = len(route_records)
        canonical_probes = sum(
            record.canonical_probe_solves for record in route_records)
        preferred_probes = sum(
            record.preferred_probe_solves for record in route_records)
    else:
        if recovery is None:
            die("internal recovery receipt lacks parsed output")
        physical = recovery.physical_rows
        aliases = recovery.alias_rows
        route_rows = 0
        canonical_probes = recovery.canonical_probe_solves
        preferred_probes = 0
    payload = {
        "job_id": job.job_id,
        "job_sha256": sha256_bytes(canonical_json_bytes(job.to_record())),
        "command": list(command),
        "returncode": 0,
        "cpu": execution.cpu,
        "start_ns": execution.start_ns,
        "end_ns": execution.end_ns,
        # Kept as an explicit schema assertion so an older, provenance-free
        # recovered receipt is rejected rather than silently upgraded.
        "recovered_after_interrupt": False,
        "stdout": paths["stdout"].relative_to(result_root).as_posix(),
        "stdout_sha256": sha256_bytes(stdout),
        "stderr": paths["stderr"].relative_to(result_root).as_posix(),
        "stderr_sha256": sha256_bytes(stderr),
        "logical_rows": logical,
        "physical_rows": physical,
        "alias_rows": aliases,
        "route_rows": route_rows,
        "canonical_probe_solves": canonical_probes,
        "preferred_probe_solves": preferred_probes,
    }
    if payload["end_ns"] < payload["start_ns"]:
        die("job wall-clock timestamps are reversed")
    return sealed_record(SCHEMA_PREFIX + ".job_receipt.v1", payload)


RECEIPT_FIELDS = {
    "schema", "job_id", "job_sha256", "command", "returncode", "cpu",
    "start_ns", "end_ns", "recovered_after_interrupt", "stdout",
    "stdout_sha256", "stderr", "stderr_sha256", "logical_rows",
    "physical_rows", "alias_rows", "route_rows", "canonical_probe_solves",
    "preferred_probe_solves", "self_sha256_excluding_field",
}


def load_canonical_object(path: Path, context: str) -> Dict[str, Any]:
    raw = common.stable_bytes(path)
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        die("{} is not canonical JSON".format(context))
    if not isinstance(value, dict) or canonical_json_bytes(value) != raw:
        die("{} is not a canonical JSON object".format(context))
    return value


def verify_job_receipt(
    job: JobSpec,
    result_root: Path,
    paths: Mapping[str, Path],
    binary: Path,
    taskset_path: Path,
    allowed_cpus: Optional[Set[int]] = None,
) -> CompletedJob:
    receipt = load_canonical_object(paths["receipt"], job.job_id + " receipt")
    verify_sealed_record(receipt, SCHEMA_PREFIX + ".job_receipt.v1")
    if set(receipt) != RECEIPT_FIELDS:
        die("{} receipt fields do not exactly match the schema".format(job.job_id))
    numeric_fields = (
        "returncode", "logical_rows", "physical_rows", "alias_rows",
        "route_rows", "canonical_probe_solves", "preferred_probe_solves",
    )
    if any(not is_plain_int(receipt.get(field)) for field in numeric_fields):
        die("{} receipt count fields are not canonical integers".format(job.job_id))
    expected_job_hash = sha256_bytes(canonical_json_bytes(job.to_record()))
    stdout = common.stable_bytes(paths["stdout"])
    stderr = common.stable_bytes(paths["stderr"])
    if (receipt.get("job_id") != job.job_id or
            receipt.get("job_sha256") != expected_job_hash or
            receipt.get("stdout") !=
                paths["stdout"].relative_to(result_root).as_posix() or
            receipt.get("stderr") !=
                paths["stderr"].relative_to(result_root).as_posix() or
            receipt.get("stdout_sha256") != sha256_bytes(stdout) or
            receipt.get("stderr_sha256") != sha256_bytes(stderr) or stderr or
            receipt.get("logical_rows") != job.expected_logical_rows or
            receipt.get("returncode") != 0):
        die("{} receipt/output binding mismatch".format(job.job_id))
    recovered = receipt.get("recovered_after_interrupt")
    if recovered is not False:
        die("{} receipt has provenance-free interrupted recovery".format(
            job.job_id))
    command = receipt.get("command")
    base_command = job.command(binary, result_root)
    cpu = receipt.get("cpu")
    start_ns = receipt.get("start_ns")
    end_ns = receipt.get("end_ns")
    if (not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0 or
            (allowed_cpus is not None and cpu not in allowed_cpus) or
            not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            start_ns < 0 or end_ns < start_ns):
        die("{} receipt execution fields are invalid".format(job.job_id))
    expected_command = (str(taskset_path), "-c", str(cpu), *base_command)
    if (not isinstance(command, list) or
            any(not isinstance(value, str) for value in command) or
            tuple(command) != expected_command):
        die("{} receipt command binding mismatch".format(job.job_id))
    routes, recovery = parse_job_bytes(job, stdout, result_root)
    expected_physical = (
        0 if job.kind == "route" else
        recovery.physical_rows if recovery is not None else -1)
    expected_alias = (
        sum(1 - record.direct for record in routes) if job.kind == "route" else
        recovery.alias_rows if recovery is not None else -1)
    expected_canonical_probes = (
        sum(record.canonical_probe_solves for record in routes)
        if job.kind == "route" else
        recovery.canonical_probe_solves if recovery is not None else -1)
    expected_preferred_probes = (
        sum(record.preferred_probe_solves for record in routes)
        if job.kind == "route" else 0)
    if (receipt.get("physical_rows") != expected_physical or
            receipt.get("alias_rows") != expected_alias or
            receipt.get("route_rows") != len(routes) or
            receipt.get("canonical_probe_solves") != expected_canonical_probes or
            receipt.get("preferred_probe_solves") != expected_preferred_probes):
        die("{} receipt parsed-count mismatch".format(job.job_id))
    return CompletedJob(job, receipt, routes, recovery)


Executor = Callable[
    [Tuple[str, ...], float, threading.Event, common.ProcessRegistry, int],
    ExecutionResult,
]
ThermalFactory = Callable[[Path, threading.Event], Any]


THERMAL_SUMMARY_FIELDS = {
    "samples", "sealed_samples_including_baseline", "cpu_busy_min_pct",
    "cpu_tctl_max_c", "dimm_max_c", "dimm_read_errors_max",
    "edac_ce_delta", "edac_ue_delta", "thermal_limit_c",
    "thermal_high_samples", "thermal_high_max_consecutive_samples",
    "guard_poll_iterations", "guard_samples", "guard_high_samples",
    "guard_limit_c", "guard_error",
}
THERMAL_POLICY_FIELDS = {
    "cpu_limit_c", "dimm_limit_c", "timing_cpu_limit_c",
    "timing_dimm_limit_c", "consecutive_samples", "stale_seconds",
    "min_cpu_busy_pct", "edac_ce_delta", "edac_ue_delta",
}


def validate_thermal_policy(contract: Mapping[str, Any]) -> Dict[str, Any]:
    policy = contract.get("thermal_policy")
    if not isinstance(policy, dict) or set(policy) != THERMAL_POLICY_FIELDS:
        die("frozen campaign thermal policy fields changed")
    float_fields = (
        "cpu_limit_c", "dimm_limit_c", "timing_cpu_limit_c",
        "timing_dimm_limit_c", "stale_seconds", "min_cpu_busy_pct",
    )
    if any(not isinstance(policy.get(field), (int, float)) or
           isinstance(policy.get(field), bool) or
           not math.isfinite(float(policy[field])) for field in float_fields):
        die("frozen campaign thermal policy has a nonfinite value")
    if (float(policy["cpu_limit_c"]) != float(policy["dimm_limit_c"]) or
            float(policy["cpu_limit_c"]) != 90.0 or
            float(policy["stale_seconds"]) != 5.0 or
            float(policy["min_cpu_busy_pct"]) != 95.0 or
            not is_plain_int(policy.get("consecutive_samples"), 1) or
            policy["consecutive_samples"] != 3 or
            not is_plain_int(policy.get("edac_ce_delta")) or
            not is_plain_int(policy.get("edac_ue_delta")) or
            policy["edac_ce_delta"] != 0 or policy["edac_ue_delta"] != 0):
        die("frozen campaign thermal safety policy changed")
    return dict(policy)


def frozen_taskset_path(contract: Mapping[str, Any]) -> Path:
    tools = contract.get("system_tools")
    entry = tools.get("taskset") if isinstance(tools, dict) else None
    if (not isinstance(entry, dict) or set(entry) != {"path", "sha256"} or
            not isinstance(entry.get("path"), str) or
            not isinstance(entry.get("sha256"), str) or
            not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"])):
        die("frozen contract lacks exact taskset path/SHA256 provenance")
    declared = Path(entry["path"])
    if not declared.is_absolute():
        die("frozen taskset path is not absolute")
    try:
        resolved = declared.resolve(strict=True)
    except OSError as error:
        die("frozen taskset is unavailable: {}".format(error))
    if resolved != declared or not os.access(str(declared), os.X_OK):
        die("frozen taskset is indirect or nonexecutable")
    tool_bytes = common.stable_bytes(declared)
    if sha256_bytes(tool_bytes) != entry["sha256"]:
        die("frozen taskset SHA256 changed")
    return declared


def validate_thermal_summary(
    summary: Mapping[str, Any], policy: Mapping[str, Any],
    require_busy: bool = True,
) -> Dict[str, Any]:
    if not isinstance(summary, dict) or set(summary) != THERMAL_SUMMARY_FIELDS:
        die("thermal guard summary fields are not canonical")
    integer_fields = (
        "samples", "sealed_samples_including_baseline",
        "dimm_read_errors_max", "edac_ce_delta", "edac_ue_delta",
        "thermal_high_samples", "thermal_high_max_consecutive_samples",
        "guard_poll_iterations", "guard_samples", "guard_high_samples",
    )
    if any(not is_plain_int(summary.get(field)) for field in integer_fields):
        die("thermal guard summary counters are not canonical integers")
    float_fields = (
        "cpu_busy_min_pct", "cpu_tctl_max_c", "dimm_max_c",
        "thermal_limit_c", "guard_limit_c",
    )
    if any(not isinstance(summary.get(field), (int, float)) or
           isinstance(summary.get(field), bool) or
           not math.isfinite(float(summary[field])) for field in float_fields):
        die("thermal guard summary contains a nonfinite measurement")
    samples = int(summary["samples"])
    if (samples <= 0 or
            summary["sealed_samples_including_baseline"] != samples + 1 or
            not 0.0 <= float(summary["cpu_busy_min_pct"]) <= 100.0 or
            not 0.0 <= float(summary["cpu_tctl_max_c"]) <= 120.0 or
            not 0.0 <= float(summary["dimm_max_c"]) <= 120.0 or
            summary["thermal_high_samples"] > samples + 1 or
            summary["thermal_high_max_consecutive_samples"] >
                summary["thermal_high_samples"] or
            summary["guard_samples"] > samples or
            summary["guard_high_samples"] > summary["guard_samples"] or
            summary["guard_high_samples"] >
                summary["thermal_high_samples"] or
            float(summary["thermal_limit_c"]) !=
                float(policy["cpu_limit_c"]) or
            float(summary["guard_limit_c"]) !=
                float(policy["cpu_limit_c"])):
        die("thermal guard summary arithmetic/policy binding failed")
    if summary["guard_error"] is not None:
        die("thermal guard failed: {}".format(summary["guard_error"]))
    if (require_busy and float(summary["cpu_busy_min_pct"]) <
            float(policy["min_cpu_busy_pct"])):
        die("CPU utilization fell below the frozen campaign gate")
    if (summary["thermal_high_max_consecutive_samples"] >=
            policy["consecutive_samples"] or
            summary["dimm_read_errors_max"] != 0 or
            summary["edac_ce_delta"] != policy["edac_ce_delta"] or
            summary["edac_ue_delta"] != policy["edac_ue_delta"]):
        die("memory telemetry changed during the campaign")
    return dict(summary)


def subprocess_executor(
    command: Tuple[str, ...],
    timeout: float,
    abort: threading.Event,
    registry: common.ProcessRegistry,
    cpu: int,
) -> ExecutionResult:
    start_ns = time.time_ns()
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True)
    registered = False
    try:
        registry.add(process)
        registered = True
        deadline = time.monotonic() + timeout
        while True:
            remaining = deadline - time.monotonic()
            try:
                stdout, stderr = process.communicate(
                    timeout=max(0.001, min(0.25, remaining)))
                break
            except subprocess.TimeoutExpired:
                if not abort.is_set() and remaining > 0:
                    continue
                common.stop_and_reap_process_group(process)
                reason = "campaign abort" if abort.is_set() else \
                    "timeout after {:g}s".format(timeout)
                die("job process terminated on " + reason)
    finally:
        leader_live = process.poll() is None
        descendants_live = common.process_group_exists(process)
        if leader_live or descendants_live:
            common.stop_and_reap_process_group(process)
        else:
            common.close_process_streams(process)
        if registered:
            registry.remove(process)
        if process.poll() is None or common.process_group_exists(process):
            die("job child process group cleanup was not proven")
        if not leader_live and descendants_live:
            die("job left an unexpected child process")
    return ExecutionResult(
        int(process.returncode), stdout, stderr, start_ns, time.time_ns(), cpu)


def phase_completion_path(result_root: Path, ledger: JobLedger) -> Path:
    return result_root / "completed" / (
        "{}-{}.json".format(ledger.phase, ledger.kind))


def phase_thermal_path(result_root: Path, ledger: JobLedger) -> Path:
    return result_root / "thermal" / (
        ledger.phase + "-" + ledger.kind + ".csv")


def verify_hash_sidecar(path: Path, data: bytes, context: str) -> None:
    digest = sha256_bytes(data)
    expected = (digest + "  " + path.name + "\n").encode("ascii")
    if common.stable_bytes(path.with_suffix(path.suffix + ".sha256")) != expected:
        die("{} SHA256 sidecar mismatch".format(context))


def discard_incomplete_phase_artifacts(
    result_root: Path, ledger: JobLedger,
) -> None:
    """Restart a phase whose receipts were never bound to sealed telemetry.

    A job receipt proves its process/CPU/output, but the phase completion is
    what binds that receipt set to a particular thermal interval.  Keeping
    receipts after a missing completion and overwriting the interval would
    splice evidence from different executions.
    """
    paths: List[Path] = []
    for job in ledger.jobs:
        job_files = job_paths(result_root, job)
        for name in ("receipt", "stdout", "stderr"):
            if path_present(job_files[name]):
                paths.append(job_files[name])
    thermal = phase_thermal_path(result_root, ledger)
    if path_present(thermal):
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


class JobRunner:
    def __init__(
        self,
        result_root: Path,
        binary: Path,
        cpus: Sequence[int],
        workers: int,
        timeout: float,
        executor: Optional[Executor] = None,
        thermal_factory: Optional[ThermalFactory] = None,
    ) -> None:
        if (not isinstance(workers, int) or isinstance(workers, bool) or
                workers <= 0 or workers != len(cpus) or
                not isinstance(timeout, (int, float)) or isinstance(timeout, bool) or
                not math.isfinite(timeout) or
                timeout <= 0 or len(set(cpus)) != len(cpus) or
                any(not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
                    for cpu in cpus)):
            die("job runner worker/CPU/timeout policy is invalid")
        self.result_root = result_root.resolve()
        self.binary = binary.absolute()
        self.cpus = tuple(cpus)
        self.workers = workers
        self.timeout = timeout
        self.executor = executor if executor is not None else subprocess_executor
        self.thermal_factory = (
            thermal_factory if thermal_factory is not None else common.ThermalGuard)

    def _run_one(
        self,
        job: JobSpec,
        pool: common.CpuPool,
        abort: threading.Event,
        registry: common.ProcessRegistry,
        taskset_path: Path,
    ) -> CompletedJob:
        paths = job_paths(self.result_root, job)
        with job_lock(paths["lock"]):
            if path_present(paths["receipt"]):
                return verify_job_receipt(
                    job, self.result_root, paths, self.binary, taskset_path,
                    set(self.cpus))
            stdout_exists = path_present(paths["stdout"])
            stderr_exists = path_present(paths["stderr"])
            if stdout_exists or stderr_exists:
                discard_partial_output_pair(paths, job.job_id)
            command_without_taskset = job.command(self.binary, self.result_root)
            if abort.is_set():
                die("campaign abort set before {}".format(job.job_id))
            if job.kind == "candidate":
                # Validate the exact immutable route bytes before allowing the
                # native process to consume them.
                parse_route_output(
                    route_job_from_candidate(job),
                    candidate_route_bytes(job, self.result_root))
            cpu = pool.acquire()
            try:
                if abort.is_set():
                    die("campaign abort set before {} launch".format(job.job_id))
                command = (
                    str(taskset_path), "-c", str(cpu), *command_without_taskset)
                execution = self.executor(
                    tuple(command), self.timeout, abort, registry, cpu)
            finally:
                pool.release(cpu)
            validate_execution_result(execution)
            if execution.returncode != 0 or execution.stderr:
                die("{} failed with rc={}".format(
                    job.job_id, execution.returncode))
            durable_atomic_write(paths["stderr"], execution.stderr)
            durable_atomic_write(paths["stdout"], execution.stdout)
            routes, recovery = parse_job_bytes(
                job, execution.stdout, self.result_root)
            receipt = receipt_for(
                job, command, self.result_root, paths, execution,
                routes, recovery)
            durable_atomic_write(paths["receipt"], canonical_json_bytes(receipt))
            return verify_job_receipt(
                job, self.result_root, paths, self.binary, taskset_path,
                set(self.cpus))

    def run_ledger(
        self,
        ledger: JobLedger,
        thermal_path: Path,
        install_signal_handlers: bool = True,
    ) -> Tuple[CompletedJob, ...]:
        ledger.record()
        # Do not create even a lock artifact until the public-freeze trust
        # boundary has passed.  The locked implementation repeats this check
        # so the exact contract is fresh at the point of use.
        verify_frozen_controller_runtime(self.result_root)
        phase_lock_path = self.result_root / "locks" / (
            ledger.phase + "-" + ledger.kind + ".lock")
        with job_lock(phase_lock_path):
            return self._run_ledger_locked(
                ledger, thermal_path, install_signal_handlers)

    def _run_ledger_locked(
        self,
        ledger: JobLedger,
        thermal_path: Path,
        install_signal_handlers: bool,
    ) -> Tuple[CompletedJob, ...]:
        _prepare_record, contract = verify_frozen_controller_runtime(
            self.result_root)
        if not isinstance(thermal_path, Path):
            die("job runner requires the live CPU/DIMM thermal log path")
        if not isinstance(contract, dict):
            die("frozen controller contract is not an object")
        expected_cpus = contract.get("cpu_set")
        expected_workers = contract.get("workers")
        expected_timeout = contract.get("timeout_seconds")
        expected_thermal = contract.get("thermal")
        thermal_policy = validate_thermal_policy(contract)
        taskset_path = frozen_taskset_path(contract)
        if (not isinstance(expected_cpus, list) or
                any(not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
                    for cpu in expected_cpus) or
                tuple(expected_cpus) != self.cpus or
                not isinstance(expected_workers, int) or
                isinstance(expected_workers, bool) or
                expected_workers != self.workers or
                not isinstance(expected_timeout, (int, float)) or
                isinstance(expected_timeout, bool) or
                float(expected_timeout) != float(self.timeout) or
                not isinstance(expected_thermal, str)):
            die("job runner CPU/worker/timeout policy differs from the frozen contract")
        try:
            actual_thermal = thermal_path.resolve(strict=True)
            frozen_thermal = Path(expected_thermal).resolve(strict=True)
        except OSError as error:
            die("frozen thermal log is unavailable: {}".format(error))
        if actual_thermal != frozen_thermal:
            die("job runner thermal log differs from the frozen contract")
        expected_binary = (
            self.result_root / "frozen" / "wirehair_v2_bench").resolve(strict=True)
        if (self.binary.resolve(strict=True) != expected_binary or
                self.binary != expected_binary):
            die("job runner binary is not the frozen controller binary")
        # Building the record validates contiguous IDs and all job invariants
        # before any subprocess is launched.
        ledger.record()
        _write_job_ledger(self.result_root, ledger)
        completion_path = phase_completion_path(self.result_root, ledger)
        completion_sidecar = completion_path.with_suffix(
            completion_path.suffix + ".sha256")
        if path_present(completion_path):
            return _verify_phase_completion(
                self.result_root, self.binary, taskset_path, ledger,
                set(self.cpus), thermal_policy)
        if path_present(completion_sidecar):
            die("phase completion sidecar exists without its sealed record")

        # A missing phase seal means prior receipts were never bound to an
        # immutable thermal interval.  Restart the exact ledger instead of
        # mixing old process evidence with a newly overwritten interval.
        discard_incomplete_phase_artifacts(self.result_root, ledger)
        if not ledger.jobs:
            ordered: Tuple[CompletedJob, ...] = ()
            _write_phase_completion(
                self.result_root, self.binary, taskset_path, ledger, ordered,
                set(self.cpus), None, None)
            return ordered

        abort = threading.Event()
        registry = common.ProcessRegistry()
        pool = common.CpuPool(self.cpus, self.workers)
        guard = (
            common.CampaignSignalGuard(abort, registry)
            if install_signal_handlers else nullcontext())
        thermal = None
        thermal_summary: Optional[Dict[str, Any]] = None
        campaign_error: Optional[BaseException] = None
        cleanup_error: Optional[BaseException] = None
        completed: Dict[str, CompletedJob] = {}
        with guard:
            executor: Optional[ThreadPoolExecutor] = None
            futures: Dict[Future[CompletedJob], JobSpec] = {}
            try:
                thermal = self.thermal_factory(thermal_path, abort)
                thermal.start()
                executor = ThreadPoolExecutor(max_workers=self.workers)
                futures = {
                    executor.submit(
                        self._run_one, job, pool, abort, registry,
                        taskset_path): job
                    for job in ledger.jobs
                }
                for future in as_completed(futures):
                    result = future.result()
                    completed[result.job.job_id] = result
                    if getattr(thermal, "error", None) is not None:
                        die("thermal guard failed: {}".format(thermal.error))
            except BaseException as error:
                campaign_error = error
                abort.set()
                registry.signal_all(signal.SIGTERM)
                for future in futures:
                    future.cancel()
            finally:
                # This ordering is deliberate: on the first failed future,
                # abort and SIGTERM happen before shutdown waits for workers.
                try:
                    if executor is not None:
                        executor.shutdown(wait=True)
                    if registry.count() != 0:
                        registry.drain()
                        die("phase workers left registered child processes")
                except BaseException as error:
                    cleanup_error = error
                finally:
                    if thermal is not None and getattr(thermal, "started", False):
                        try:
                            thermal_output = phase_thermal_path(
                                self.result_root, ledger)
                            parent_fd = common.open_durable_directory(
                                thermal_output.parent, create=True)
                            os.close(parent_fd)
                            thermal_summary = thermal.finish(thermal_output)
                            fsync_existing_regular_file(thermal_output)
                        except BaseException as error:
                            cleanup_error = error
        if cleanup_error is not None:
            raise cleanup_error from campaign_error
        if thermal_summary is None:
            if campaign_error is not None:
                raise campaign_error
            die("thermal guard did not produce a summary")
        validated_summary = validate_thermal_summary(
            thermal_summary, thermal_policy,
            require_busy=campaign_error is None)
        if campaign_error is not None:
            raise campaign_error
        if registry.count() != 0:
            die("phase completion left registered child processes")
        if frozen_taskset_path(contract) != taskset_path:
            die("frozen taskset identity changed during the phase")
        ordered = tuple(completed[job.job_id] for job in ledger.jobs)
        if len(ordered) != len(ledger.jobs):
            die("job ledger completion cardinality mismatch")
        _write_phase_completion(
            self.result_root, self.binary, taskset_path, ledger, ordered,
            set(self.cpus), phase_thermal_path(self.result_root, ledger),
            validated_summary)
        return ordered


def _phase_completion_record(
    result_root: Path,
    binary: Path,
    taskset_path: Path,
    ledger: JobLedger,
    completed: Sequence[CompletedJob],
    allowed_cpus: Set[int],
    thermal_path: Optional[Path],
    thermal_summary: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    """Re-parse every output and seal exact physical and semantic totals."""
    if tuple(item.job.job_id for item in completed) != tuple(
            job.job_id for job in ledger.jobs):
        die("phase completion jobs are missing, duplicated, or reordered")
    verified = tuple(
        verify_job_receipt(
            job, result_root, job_paths(result_root, job), binary,
            taskset_path,
            allowed_cpus)
        for job in ledger.jobs
    )
    receipt_rows: List[Dict[str, Any]] = []
    logical = physical = aliases = route_rows = 0
    canonical_probes = preferred_probes = 0
    successes = rank_failures = codec_errors = 0
    direct_routes = fallback_routes = noop_routes = 0
    semantic_records: List[Dict[str, Any]] = []
    for item in verified:
        receipt = item.receipt
        receipt_bytes = canonical_json_bytes(receipt)
        receipt_rows.append({
            "job_id": item.job.job_id,
            "receipt_sha256": sha256_bytes(receipt_bytes),
            "receipt_self_sha256": receipt["self_sha256_excluding_field"],
            "stdout_sha256": receipt["stdout_sha256"],
            "logical_rows": receipt["logical_rows"],
            "physical_rows": receipt["physical_rows"],
            "alias_rows": receipt["alias_rows"],
            "route_rows": receipt["route_rows"],
        })
        logical += int(receipt["logical_rows"])
        physical += int(receipt["physical_rows"])
        aliases += int(receipt["alias_rows"])
        route_rows += int(receipt["route_rows"])
        canonical_probes += int(receipt["canonical_probe_solves"])
        preferred_probes += int(receipt["preferred_probe_solves"])
        for route in item.route_records:
            direct_routes += route.direct
            fallback_routes += route.fallback
            noop_routes += route.no_op
            semantic_records.append({
                "K": route.K, "width": route.width,
                "p": route.preferred_attempt,
                "a0": route.canonical_attempt,
                "actual": route.actual_attempt,
                "valid": route.preferred_valid,
                "fallback": route.fallback,
                "no_op": route.no_op,
                "direct": route.direct,
            })
        if item.recovery is not None:
            for row in item.recovery.rows:
                successes += int(row.metrics.success)
                rank_failures += int(row.metrics.failure)
                codec_errors += row.metrics.error
                semantic_records.append({
                    "K": row.K, "width": row.width,
                    "root_index": item.job.root_index,
                    "schedule": item.job.schedule, "loss": item.job.loss,
                    "p": row.preferred_attempt,
                    "result": row.metrics.result,
                    "rank_fail": row.metrics.rank_fail,
                    "error": row.metrics.error,
                    "heavy_shortfall": row.metrics.heavy_shortfall,
                    "inactivated": row.metrics.inactivated,
                    "binary_deficit": row.metrics.binary_deficit,
                    "heavy_gain": row.metrics.heavy_gain,
                    "block_xors": row.metrics.block_xors,
                    "block_muladds": row.metrics.block_muladds,
                })
    ledger_record = ledger.record()
    expected_logical = ledger_record["expected_logical_rows"]
    expected_physical = ledger_record["expected_physical_rows"]
    if logical != expected_logical or physical != expected_physical:
        die("phase completion totals disagree with the job ledger")
    if ledger.kind == "candidate" and logical != physical + aliases:
        die("candidate phase logical/physical/alias identity failed")
    if ledger.kind == "route" and (
            route_rows != logical or direct_routes + aliases != logical or
            fallback_routes + noop_routes != aliases):
        die("route phase row/classification identity failed")
    if ledger.kind == "control" and (physical != logical or aliases or route_rows):
        die("control phase row identity failed")
    if (ledger.kind in ("control", "candidate") and
            successes + rank_failures + codec_errors != physical):
        die("recovery phase outcome classes do not cover physical solves")
    if ledger.jobs:
        if thermal_path is None or not isinstance(thermal_summary, dict):
            die("nonempty phase completion lacks thermal provenance")
        expected_thermal = phase_thermal_path(result_root, ledger)
        if thermal_path.absolute() != expected_thermal.absolute():
            die("phase completion thermal artifact path changed")
        thermal_bytes = common.stable_bytes(expected_thermal)
        thermal_artifact: Optional[Dict[str, Any]] = {
            "path": expected_thermal.relative_to(result_root).as_posix(),
            "sha256": sha256_bytes(thermal_bytes),
            "summary": dict(thermal_summary),
        }
    else:
        if thermal_path is not None or thermal_summary is not None:
            die("empty phase completion unexpectedly claims a thermal interval")
        thermal_artifact = None
    return sealed_record(
        SCHEMA_PREFIX + ".phase_completion.v1",
        {
            "phase": ledger.phase,
            "kind": ledger.kind,
            "job_ledger_sha256": sha256_bytes(canonical_json_bytes(ledger_record)),
            "job_count": len(verified),
            "jobs": receipt_rows,
            "logical_rows": logical,
            "physical_rows": physical,
            "alias_rows": aliases,
            "route_rows": route_rows,
            "canonical_probe_solves": canonical_probes,
            "preferred_probe_solves": preferred_probes,
            "successes": successes,
            "rank_failures": rank_failures,
            "codec_errors": codec_errors,
            "direct_routes": direct_routes,
            "fallback_routes": fallback_routes,
            "noop_routes": noop_routes,
            "thermal_artifact": thermal_artifact,
            "semantic_rows_sha256": sha256_bytes(
                canonical_json_bytes(semantic_records)),
        },
    )


def _write_phase_completion(
    result_root: Path,
    binary: Path,
    taskset_path: Path,
    ledger: JobLedger,
    completed: Sequence[CompletedJob],
    allowed_cpus: Set[int],
    thermal_path: Optional[Path],
    thermal_summary: Optional[Mapping[str, Any]],
) -> str:
    record = _phase_completion_record(
        result_root, binary, taskset_path, ledger, completed, allowed_cpus,
        thermal_path, thermal_summary)
    path = phase_completion_path(result_root, ledger)
    return write_hashed_artifact(path, canonical_json_bytes(record))


def _verify_phase_completion(
    result_root: Path,
    binary: Path,
    taskset_path: Path,
    ledger: JobLedger,
    allowed_cpus: Set[int],
    thermal_policy: Mapping[str, Any],
) -> Tuple[CompletedJob, ...]:
    """Verify a fully sealed phase without starting or replacing telemetry."""
    path = phase_completion_path(result_root, ledger)
    encoded = common.stable_bytes(path)
    sidecar = path.with_suffix(path.suffix + ".sha256")
    sidecar_present = path_present(sidecar)
    if sidecar_present:
        verify_hash_sidecar(path, encoded, "phase completion")
    record = load_canonical_object(path, "phase completion")
    verify_sealed_record(record, SCHEMA_PREFIX + ".phase_completion.v1")
    completed = tuple(
        verify_job_receipt(
            job, result_root, job_paths(result_root, job), binary,
            taskset_path,
            allowed_cpus)
        for job in ledger.jobs
    )
    if ledger.jobs:
        artifact = record.get("thermal_artifact")
        expected_path = phase_thermal_path(result_root, ledger)
        expected_relative = expected_path.relative_to(result_root).as_posix()
        if (not isinstance(artifact, dict) or
                set(artifact) != {"path", "sha256", "summary"} or
                artifact.get("path") != expected_relative or
                not isinstance(artifact.get("sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", artifact["sha256"])):
            die("phase completion thermal binding is malformed")
        thermal_bytes = common.stable_bytes(expected_path)
        if sha256_bytes(thermal_bytes) != artifact["sha256"]:
            die("completed phase thermal artifact changed")
        summary = validate_thermal_summary(
            artifact.get("summary"), thermal_policy)
        expected = _phase_completion_record(
            result_root, binary, taskset_path, ledger, completed, allowed_cpus,
            expected_path, summary)
    else:
        if record.get("thermal_artifact") is not None:
            die("empty completed phase claims a thermal artifact")
        expected = _phase_completion_record(
            result_root, binary, taskset_path, ledger, completed,
            allowed_cpus, None, None)
    if canonical_json_bytes(expected) != encoded:
        die("phase completion is not the exact recomputed record")
    if not sidecar_present:
        # The completion record is written and directory-fsynced before its
        # checksum sidecar.  A hard interruption in that narrow window must
        # not strand an otherwise exact, fully recomputable phase.  Repair
        # only after receipts, telemetry, schema, and canonical bytes have all
        # been independently revalidated above; an existing bad sidecar is
        # never replaced.
        digest = sha256_bytes(encoded)
        durable_write_once_or_same(
            sidecar, (digest + "  " + path.name + "\n").encode("ascii"))
    return completed


def route_artifacts_from_completed(
    result_root: Path, completed: Sequence[CompletedJob],
) -> Dict[str, RouteArtifact]:
    artifacts: Dict[str, RouteArtifact] = {}
    for result in completed:
        if result.job.kind != "route" or not result.route_records:
            die("route artifact source is not a completed route job")
        path = job_paths(result_root, result.job)["stdout"]
        relative = path.relative_to(result_root).as_posix()
        stdout_bytes = common.stable_bytes(path)
        digest = sha256_bytes(stdout_bytes)
        if digest != result.receipt.get("stdout_sha256"):
            die("completed route output changed before artifact collection")
        artifact = RouteArtifact(
            result.job.job_id, relative, digest,
            result.route_records)
        if artifact.job_id in artifacts:
            die("duplicate completed route job")
        artifacts[artifact.job_id] = artifact
    return artifacts


def route_semantics(record: RouteRecord) -> Tuple[Any, ...]:
    return (
        record.K, record.width, record.route_status,
        record.preferred_attempt, record.canonical_attempt,
        record.actual_attempt, record.preferred_valid, record.fallback,
        record.no_op, record.direct,
    )


class DevelopmentEvidence:
    """Validated physical evidence, with future controls kept masked by users."""

    def __init__(
        self,
        schedules: Sequence[str] = APPROVED_SCHEDULES,
        losses: Sequence[str] = APPROVED_LOSSES,
    ) -> None:
        self.schedules = validate_grid_subset(
            schedules, APPROVED_SCHEDULES, "evidence schedules")
        self.losses = validate_grid_subset(
            losses, APPROVED_LOSSES, "evidence losses")
        self.controls: Dict[CellKey, RecoveryMetrics] = {}
        self.control_canonical: Dict[Tuple[int, int], int] = {}
        self.routes: Dict[Tuple[int, int, int], RouteRecord] = {}
        self.route_canonical: Dict[Tuple[int, int], int] = {}
        self.round_routes: Dict[Tuple[str, int, int, int], RouteRecord] = {}
        self.candidates: Dict[Tuple[int, CellKey], RecoveryMetrics] = {}

    def add_completed(self, results: Sequence[CompletedJob]) -> None:
        for completed in results:
            job = completed.job
            if job.kind == "route":
                for route in completed.route_records:
                    validate_route_record_value(route)
                    canonical_key = (route.K, route.width)
                    canonical = self.route_canonical.setdefault(
                        canonical_key, route.canonical_attempt)
                    if canonical != route.canonical_attempt:
                        die("canonical attempt changed between route batches")
                    control_canonical = self.control_canonical.get(canonical_key)
                    if (control_canonical is not None and
                            control_canonical != route.canonical_attempt):
                        die("route/control canonical attempts disagree")
                    key = (route.K, route.width, route.preferred_attempt)
                    old = self.routes.get(key)
                    if old is not None and route_semantics(old) != route_semantics(route):
                        die("route semantics changed between rounds")
                    self.routes[key] = route
                    round_key = (job.round_name, *key)
                    if round_key in self.round_routes:
                        die("duplicate route evidence in one round")
                    self.round_routes[round_key] = route
                continue
            if completed.recovery is None or job.root_index is None or \
                    job.schedule is None or job.loss is None:
                die("completed recovery job lacks its stratum identity")
            for row in completed.recovery.rows:
                validate_recovery_metrics_value(
                    row.metrics, "completed recovery evidence")
                cell = CellKey(
                    row.K, row.width, job.root_index, job.schedule, job.loss)
                if job.kind == "control":
                    canonical_key = (row.K, row.width)
                    canonical = self.control_canonical.setdefault(
                        canonical_key, row.canonical_attempt)
                    if canonical != row.canonical_attempt:
                        die("control canonical attempt changed between strata")
                    route_canonical = self.route_canonical.get(canonical_key)
                    if (route_canonical is not None and
                            route_canonical != row.canonical_attempt):
                        die("control/route canonical attempts disagree")
                    if cell in self.controls:
                        die("duplicate physical control cell")
                    self.controls[cell] = row.metrics
                else:
                    key = (row.preferred_attempt, cell)
                    if key in self.candidates:
                        die("duplicate physical candidate cell")
                    self.candidates[key] = row.metrics

    def control(self, cell: CellKey, allowed: Set[Panel]) -> RecoveryMetrics:
        # The allowlist is mandatory even though all six roots may already be
        # physically cached.  This is the future-root masking boundary.
        if Panel(cell.root_index, cell.width) not in allowed:
            die("analysis attempted to consume a masked future control panel")
        try:
            value = self.controls[cell]
        except KeyError:
            die("control cube is missing {}".format(cell))
        validate_recovery_metrics_value(value, "control cube cell")
        if value.error:
            die("control cube contains a codec error")
        return value

    def route(self, round_name: str, K: int, width: int, p: int) -> RouteRecord:
        try:
            return self.round_routes[(round_name, K, width, p)]
        except KeyError:
            die("round route evidence is missing K={} bb={} p={}".format(
                K, width, p))
        raise AssertionError("unreachable")

    def logical_candidate(
        self, p: int, cell: CellKey, allowed: Set[Panel],
    ) -> RecoveryMetrics:
        control = self.control(cell, allowed)
        try:
            route = self.routes[(cell.K, cell.width, p)]
        except KeyError:
            die("cumulative route evidence is missing")
        physical = self.candidates.get((p, cell))
        if route.direct:
            if physical is None:
                die("direct logical cell lacks its physical candidate solve")
            validate_recovery_metrics_value(
                physical, "physical candidate cell")
            return physical
        if physical is not None:
            die("alias/fallback cell consumed a physical candidate solve")
        return control


@dataclass(frozen=True)
class CandidateScore:
    p: int
    eligible: bool
    reasons: Tuple[str, ...]
    objective: Optional[Tuple[int, ...]]
    full_rank: Optional[Tuple[int, ...]]


GLOBAL_TABLE_OBJECTIVE_FIELDS = (
    "errors",
    "max_K_failure_multiplicity",
    "multi_failure_K_count",
    "total_failures",
    "introductions",
    "worst_adverse_axis_delta",
    "sum_positive_axis_deltas",
    "heavy_shortfall_sum",
    "binary_deficit_gt15_count",
    "binary_deficit_max",
    "binary_deficit_sum",
    "common_success_xors",
    "common_success_muladds",
)


def objective_is_canonical(values: Any, control: bool = False) -> bool:
    if (not isinstance(values, list) or len(values) != 12 or
            any(not isinstance(value, int) or isinstance(value, bool)
                for value in values)):
        return False
    if (any(value < 0 for index, value in enumerate(values) if index != 5) or
            values[5] > 0 or values[0] > values[4] or
            values[2] > values[3] or values[7] > values[8]):
        return False
    adverse_count, adverse_max, adverse_sum = values[1:4]
    if ((adverse_count == 0) != (adverse_max == 0 and adverse_sum == 0) or
            (adverse_count > 0 and
             (adverse_max <= 0 or adverse_sum < adverse_count))):
        return False
    if control and (values[0:4] != [0, 0, 0, 0] or values[5] != 0):
        return False
    return True


def failure_bit(metrics: RecoveryMetrics) -> int:
    if metrics.error:
        die("codec error cannot be interpreted as a rank failure")
    if metrics.success:
        return 0
    if metrics.failure:
        return 1
    die("recovery result is neither success nor rank failure")
    return 1


def score_candidate(
    evidence: DevelopmentEvidence,
    K: int,
    p: int,
    panels: Sequence[Panel],
    allowed: Set[Panel],
) -> CandidateScore:
    reasons: Set[str] = set()
    widths = sorted({panel.width for panel in panels})
    any_direct = False
    for width in widths:
        try:
            route = evidence.routes[(K, width, p)]
        except KeyError:
            reasons.add("missing_route")
            continue
        try:
            validate_route_record_value(route)
        except CampaignError:
            reasons.add("route_invariant")
            continue
        if ((route.K, route.width, route.preferred_attempt,
             route.route_status) != (K, width, p, "preferred")):
            reasons.add("route_invariant")
            continue
        if route.fallback:
            reasons.add("fallback")
        if not route.preferred_valid:
            reasons.add("ineligible_route")
        any_direct = any_direct or bool(route.direct)
        if width == 64 and not route.direct:
            reasons.add("not_direct_bb64")
    if not any_direct:
        reasons.add("no_direct_width")
    if reasons:
        return CandidateScore(p, False, tuple(sorted(reasons)), None, None)
    cells: List[Tuple[CellKey, RecoveryMetrics, RecoveryMetrics]] = []
    for panel in panels:
        for schedule in evidence.schedules:
            for loss in evidence.losses:
                cell = CellKey(K, panel.width, panel.root_index, schedule, loss)
                control = evidence.control(cell, allowed)
                candidate = evidence.logical_candidate(p, cell, allowed)
                if candidate.error:
                    reasons.add("candidate_error")
                cells.append((cell, control, candidate))
    if reasons:
        return CandidateScore(p, False, tuple(sorted(reasons)), None, None)
    introductions = 0
    repairs = 0
    candidate_failures = 0
    slice_deltas: Dict[Tuple[int, str, str], int] = {}
    heavy_shortfall_sum = 0
    binary_deficit_max = 0
    binary_deficit_sum = 0
    common_muladds = 0
    common_xors = 0
    common_inactivated = 0
    for cell, control, candidate in cells:
        control_failed = failure_bit(control)
        candidate_failed = failure_bit(candidate)
        introductions += int(not control_failed and candidate_failed)
        repairs += int(control_failed and not candidate_failed)
        candidate_failures += candidate_failed
        slice_key = (cell.width, cell.schedule, cell.loss)
        slice_deltas[slice_key] = (
            slice_deltas.get(slice_key, 0) +
            candidate_failed - control_failed)
        heavy_shortfall_sum += candidate.heavy_shortfall
        binary_deficit_max = max(binary_deficit_max, candidate.binary_deficit)
        binary_deficit_sum += candidate.binary_deficit
        common_inactivated += candidate.inactivated
        if control.success and candidate.success:
            common_muladds += candidate.block_muladds
            common_xors += candidate.block_xors
    positives = [delta for delta in slice_deltas.values() if delta > 0]
    objective = (
        introductions,
        len(positives),
        max(positives) if positives else 0,
        sum(positives),
        candidate_failures,
        -repairs,
        heavy_shortfall_sum,
        binary_deficit_max,
        binary_deficit_sum,
        common_muladds,
        common_xors,
        common_inactivated,
    )
    return CandidateScore(p, True, (), objective, (*objective, p))


def control_objective(
    evidence: DevelopmentEvidence,
    K: int,
    panels: Sequence[Panel],
    allowed: Set[Panel],
) -> Tuple[int, ...]:
    failures = 0
    heavy_shortfall = 0
    deficit_max = 0
    deficit_sum = 0
    muladds = 0
    xors = 0
    inactivated = 0
    for panel in panels:
        for schedule in evidence.schedules:
            for loss in evidence.losses:
                metrics = evidence.control(
                    CellKey(K, panel.width, panel.root_index, schedule, loss),
                    allowed)
                failures += failure_bit(metrics)
                heavy_shortfall += metrics.heavy_shortfall
                deficit_max = max(deficit_max, metrics.binary_deficit)
                deficit_sum += metrics.binary_deficit
                inactivated += metrics.inactivated
                if metrics.success:
                    muladds += metrics.block_muladds
                    xors += metrics.block_xors
    return (
        0, 0, 0, 0, failures, 0, heavy_shortfall, deficit_max,
        deficit_sum, muladds, xors, inactivated,
    )


def validate_round_evidence(
    evidence: DevelopmentEvidence,
    spec: RoundSpec,
    cohort: Sequence[int],
    incoming: Mapping[int, Sequence[int]],
) -> Tuple[int, int, int]:
    expected_route = {
        (spec.name, K, width, p)
        for K in cohort for p in incoming[K] for width in spec.widths
    }
    actual_route = {
        key for key in evidence.round_routes
        if key[0] == spec.name
    }
    if actual_route != expected_route:
        die("{} route evidence does not exactly cover its input".format(spec.name))
    expected_physical: Set[Tuple[int, CellKey]] = set()
    logical = 0
    for K in cohort:
        for p in incoming[K]:
            for width in spec.widths:
                route = evidence.route(spec.name, K, width, p)
                for root_index in spec.root_indexes:
                    for schedule in evidence.schedules:
                        for loss in evidence.losses:
                            logical += 1
                            if route.direct:
                                expected_physical.add((p, CellKey(
                                    K, width, root_index, schedule, loss)))
    current_roots = set(spec.root_indexes)
    actual_physical = {
        key for key in evidence.candidates
        if key[1].root_index in current_roots
    }
    if actual_physical != expected_physical:
        die("{} physical candidate evidence coverage mismatch".format(spec.name))
    physical = len(actual_physical)
    return logical, physical, logical - physical


def analyze_round(
    result_root: Path,
    evidence: DevelopmentEvidence,
    rounds: Sequence[RoundSpec],
    spec: RoundSpec,
    cohort: Sequence[int],
    incoming: Mapping[int, Sequence[int]],
    parent_survivor_sha256: str,
) -> Dict[str, Any]:
    result_root = result_root.resolve()
    verify_frozen_controller_runtime(result_root)
    validate_rounds(rounds)
    matches = tuple(item for item in rounds if item.name == spec.name)
    if len(matches) != 1 or matches[0] != spec:
        die("analyzed round spec is not the declared round object")
    validate_survivors(incoming, cohort, spec.input_max)
    if (not isinstance(parent_survivor_sha256, str) or
            not re.fullmatch(r"[0-9a-f]{64}", parent_survivor_sha256)):
        die("round parent survivor hash is invalid")
    panels = declared_panels_through(rounds, spec.name)
    allowed = set(panels)
    logical, physical, aliases = validate_round_evidence(
        evidence, spec, cohort, incoming)
    K_records: List[Dict[str, Any]] = []
    retained_total = 0
    for K in sorted(cohort):
        scores = [
            score_candidate(evidence, K, p, panels, allowed)
            for p in incoming[K]
        ]
        eligible = sorted(
            (score for score in scores if score.eligible),
            key=lambda score: score.full_rank)
        retained_scores = eligible[:spec.retain_max]
        control = control_objective(evidence, K, panels, allowed)
        if spec.require_control_improvement:
            retained_scores = [
                score for score in retained_scores
                if score.objective is not None and score.objective < control
            ]
        retained = tuple(score.p for score in retained_scores)
        retained_total += len(retained)
        score_records = []
        retained_set = set(retained)
        for score in sorted(scores, key=lambda item: item.p):
            score_records.append({
                "p": score.p,
                "eligible": score.eligible,
                "hard_reasons": list(score.reasons),
                "objective": list(score.objective)
                    if score.objective is not None else None,
                "rank": list(score.full_rank)
                    if score.full_rank is not None else None,
                "retained": score.p in retained_set,
                "strictly_better_than_control": (
                    score.objective < control
                    if score.objective is not None else False),
            })
        K_records.append({
            "K": K,
            "input": list(incoming[K]),
            "control_objective": list(control),
            "candidates": score_records,
            "eligible_ranked": [score.p for score in eligible],
            "retained": list(retained),
        })
    record = sealed_record(
        SCHEMA_PREFIX + ".survivors.v1",
        {
            "round": spec.name,
            "parent_survivor_sha256": parent_survivor_sha256,
            "input_max": spec.input_max,
            "retain_max": spec.retain_max,
            "strict_control_improvement": spec.require_control_improvement,
            "declared_panels_cumulative": [
                {"root_index": panel.root_index, "width": panel.width}
                for panel in panels
            ],
            "K": K_records,
            "K_count": len(K_records),
            "input_candidates": sum(len(incoming[K]) for K in cohort),
            "retained_candidates": retained_total,
            "current_round_logical_cells": logical,
            "current_round_physical_candidate_solves": physical,
            "current_round_alias_cells": aliases,
        },
    )
    survivors_from_record(record, spec.name, cohort, spec.retain_max)
    return record


def survivors_from_record(
    record: Mapping[str, Any],
    expected_round: str,
    cohort: Sequence[int],
    retain_max: int,
) -> Dict[int, Tuple[int, ...]]:
    verify_sealed_record(record, SCHEMA_PREFIX + ".survivors.v1")
    cohort_tuple = tuple(cohort)
    if (any(not isinstance(K, int) or isinstance(K, bool) or K < 2
            for K in cohort_tuple) or cohort_tuple != tuple(sorted(cohort_tuple)) or
            len(cohort_tuple) != len(set(cohort_tuple))):
        die("survivor ledger expected cohort is noncanonical")
    expected_fields = {
        "schema", "round", "parent_survivor_sha256", "input_max",
        "retain_max", "strict_control_improvement",
        "declared_panels_cumulative", "K", "K_count", "input_candidates",
        "retained_candidates", "current_round_logical_cells",
        "current_round_physical_candidate_solves",
        "current_round_alias_cells", "self_sha256_excluding_field",
    }
    if set(record) != expected_fields:
        die("survivor ledger fields do not exactly match the schema")
    rows = record.get("K")
    strict = record.get("strict_control_improvement")
    input_max = record.get("input_max")
    parent_hash = record.get("parent_survivor_sha256")
    if (record.get("round") != expected_round or not isinstance(rows, list) or
            not isinstance(strict, bool) or
            not isinstance(input_max, int) or isinstance(input_max, bool) or
            input_max <= 0 or input_max > MAX_ATTEMPT + 1 or
            record.get("retain_max") != retain_max or
            not isinstance(retain_max, int) or isinstance(retain_max, bool) or
            retain_max <= 0 or retain_max > input_max or
            not isinstance(parent_hash, str) or
            not re.fullmatch(r"[0-9a-f]{64}", parent_hash)):
        die("survivor ledger round/rows mismatch")
    panels = record.get("declared_panels_cumulative")
    panel_values: List[Panel] = []
    if not isinstance(panels, list) or not panels:
        die("survivor ledger has no cumulative panels")
    for panel in panels:
        if (not isinstance(panel, dict) or set(panel) != {"root_index", "width"} or
                not isinstance(panel.get("root_index"), int) or
                isinstance(panel.get("root_index"), bool) or
                panel["root_index"] < 0 or panel.get("width") not in APPROVED_WIDTHS):
            die("survivor ledger panel is malformed")
        panel_values.append(Panel(panel["root_index"], panel["width"]))
    if len(panel_values) != len(set(panel_values)):
        die("survivor ledger repeats a cumulative panel")
    result: Dict[int, Tuple[int, ...]] = {}
    total_input = total_retained = 0
    hard_reason_names = {
        "missing_route", "route_invariant", "fallback", "ineligible_route",
        "not_direct_bb64", "no_direct_width", "candidate_error",
    }
    K_row_fields = {
        "K", "input", "control_objective", "candidates",
        "eligible_ranked", "retained",
    }
    candidate_fields = {
        "p", "eligible", "hard_reasons", "objective", "rank", "retained",
        "strictly_better_than_control",
    }
    for row in rows:
        if not isinstance(row, dict) or set(row) != K_row_fields:
            die("survivor K row is not an object")
        K = row.get("K")
        incoming = row.get("input")
        control = row.get("control_objective")
        candidates = row.get("candidates")
        retained = row.get("retained")
        ranked = row.get("eligible_ranked")
        if (not isinstance(K, int) or isinstance(K, bool) or
                not isinstance(incoming, list) or
                len(incoming) > input_max or
                any(not isinstance(p, int) or isinstance(p, bool) or
                    p < 0 or p > MAX_ATTEMPT for p in incoming) or
                len(incoming) != len(set(incoming)) or
                not objective_is_canonical(control, control=True) or
                not isinstance(candidates, list) or
                not isinstance(retained, list) or not isinstance(ranked, list) or
                len(retained) > retain_max or
                any(not isinstance(p, int) or isinstance(p, bool) or
                    p < 0 or p > MAX_ATTEMPT for p in retained) or
                len(retained) != len(set(retained)) or K in result):
            die("survivor K row violates quota/order")
        score_by_p: Dict[int, Dict[str, Any]] = {}
        eligible_scores: List[Dict[str, Any]] = []
        for score in candidates:
            if not isinstance(score, dict) or set(score) != candidate_fields:
                die("survivor candidate score fields are malformed")
            p = score.get("p")
            eligible_flag = score.get("eligible")
            reasons = score.get("hard_reasons")
            objective = score.get("objective")
            rank = score.get("rank")
            retained_flag = score.get("retained")
            better = score.get("strictly_better_than_control")
            if (not isinstance(p, int) or isinstance(p, bool) or
                    p < 0 or p > MAX_ATTEMPT or p in score_by_p or
                    not isinstance(eligible_flag, bool) or
                    not isinstance(reasons, list) or
                    any(not isinstance(reason, str) for reason in reasons) or
                    reasons != sorted(set(reasons)) or
                    any(reason not in hard_reason_names for reason in reasons) or
                    not isinstance(retained_flag, bool) or
                    not isinstance(better, bool)):
                die("survivor candidate score identity is malformed")
            if eligible_flag:
                if (reasons or not objective_is_canonical(objective) or
                        not isinstance(rank, list) or
                        rank != objective + [p] or
                        better != (tuple(objective) < tuple(control))):
                    die("eligible survivor candidate rank is malformed")
                eligible_scores.append(score)
            elif (not reasons or objective is not None or rank is not None or better):
                die("ineligible survivor candidate is not hard-filtered")
            score_by_p[p] = score
        if (list(score_by_p) != sorted(score_by_p) or
                set(score_by_p) != set(incoming)):
            die("survivor candidate scores do not exactly cover the input")
        eligible_scores.sort(key=lambda score: tuple(score["rank"]))
        expected_ranked = [score["p"] for score in eligible_scores]
        if (any(not isinstance(p, int) or isinstance(p, bool) or
                p < 0 or p > MAX_ATTEMPT for p in ranked) or
                ranked != expected_ranked or len(ranked) != len(set(ranked))):
            die("survivor eligible ranking is not canonical")
        expected_retained = [
            score["p"] for score in eligible_scores[:retain_max]
            if not strict or score["strictly_better_than_control"]
        ]
        if retained != expected_retained:
            die("survivor retained set violates quota/strict-control policy")
        retained_set = set(retained)
        if any(score["retained"] != (score["p"] in retained_set)
               for score in candidates):
            die("survivor retained flags disagree with the retained ledger")
        result[K] = tuple(retained)
        total_input += len(incoming)
        total_retained += len(retained)
    logical = record.get("current_round_logical_cells")
    physical = record.get("current_round_physical_candidate_solves")
    aliases = record.get("current_round_alias_cells")
    if (list(result) != list(cohort_tuple) or
            not is_plain_int(record.get("K_count")) or
            record.get("K_count") != len(rows) or
            not is_plain_int(record.get("input_candidates")) or
            record.get("input_candidates") != total_input or
            not is_plain_int(record.get("retained_candidates")) or
            record.get("retained_candidates") != total_retained or
            not isinstance(logical, int) or isinstance(logical, bool) or logical < 0 or
            not isinstance(physical, int) or isinstance(physical, bool) or physical < 0 or
            not isinstance(aliases, int) or isinstance(aliases, bool) or aliases < 0 or
            logical != physical + aliases):
        die("survivor ledger does not exactly cover the cohort")
    return result


def write_survivor_ledger(
    result_root: Path, path: Path, record: Mapping[str, Any],
) -> str:
    result_root = result_root.resolve()
    verify_frozen_controller_runtime(result_root)
    verify_sealed_record(record, SCHEMA_PREFIX + ".survivors.v1")
    rows = record.get("K")
    retain_max = record.get("retain_max")
    round_name = record.get("round")
    if (not isinstance(rows, list) or not isinstance(round_name, str) or
            not isinstance(retain_max, int) or isinstance(retain_max, bool)):
        die("survivor ledger publication metadata is malformed")
    cohort = tuple(
        row.get("K") if isinstance(row, dict) else None for row in rows)
    survivors_from_record(record, round_name, cohort, retain_max)
    encoded = canonical_json_bytes(record)
    target = path if path.is_absolute() else result_root / path
    resolved_parent = target.absolute().parent.resolve()
    if resolved_parent != result_root and result_root not in resolved_parent.parents:
        die("survivor ledger path escapes the result root")
    return write_hashed_artifact(target, encoded)


@dataclass(frozen=True)
class DerivedRoundInput:
    survivors: Dict[int, Tuple[int, ...]]
    parent_sha256: str


def derive_round_input(
    rounds: Sequence[RoundSpec],
    spec: RoundSpec,
    cohort: Sequence[int],
    r1_parent_sha256: str,
    previous_survivor_record: Optional[Mapping[str, Any]] = None,
) -> DerivedRoundInput:
    """Derive, rather than repopulate, a round's immutable candidate input.

    R1 receives exactly p=0..255.  Every later round is reconstructed solely
    from the immediately preceding canonical survivor record, which is the
    module-level no-backfill boundary.
    """
    validate_rounds(rounds)
    matches = [index for index, value in enumerate(rounds) if value == spec]
    if len(matches) != 1:
        die("round input spec is not uniquely declared")
    index = matches[0]
    if index == 0:
        if previous_survivor_record is not None:
            die("R1 may not consume a previous survivor record")
        canonical_context_hash(r1_parent_sha256)
        return DerivedRoundInput(
            initial_survivors(cohort), r1_parent_sha256)
    if previous_survivor_record is None:
        die("later round lacks its previous survivor record")
    previous = rounds[index - 1]
    survivors = survivors_from_record(
        previous_survivor_record, previous.name, cohort, previous.retain_max)
    expected_panels = [
        {"root_index": panel.root_index, "width": panel.width}
        for panel in declared_panels_through(rounds, previous.name)
    ]
    if (previous_survivor_record.get("input_max") != previous.input_max or
            previous_survivor_record.get("retain_max") != previous.retain_max or
            previous_survivor_record.get("strict_control_improvement") !=
                previous.require_control_improvement or
            previous_survivor_record.get("declared_panels_cumulative") !=
                expected_panels):
        die("previous survivor record does not match its declared round policy")
    parent_bytes = canonical_json_bytes(previous_survivor_record)
    return DerivedRoundInput(survivors, sha256_bytes(parent_bytes))


def validate_global_table_objective_report(
    report: Mapping[str, Any],
) -> Dict[str, int]:
    """Validate the exact frozen whole-development-table objective fields."""
    if (not isinstance(report, dict) or
            set(report) != set(GLOBAL_TABLE_OBJECTIVE_FIELDS) or
            any(not is_plain_int(report.get(field))
                for field in GLOBAL_TABLE_OBJECTIVE_FIELDS)):
        die("global table objective report is not canonical")
    if (report["worst_adverse_axis_delta"] >
            report["sum_positive_axis_deltas"] or
            report["max_K_failure_multiplicity"] >
            report["total_failures"] or
            report["introductions"] > report["total_failures"] or
            2 * report["multi_failure_K_count"] >
                report["total_failures"] or
            report["binary_deficit_max"] > report["binary_deficit_sum"]):
        die("global table objective report arithmetic is impossible")
    return {field: int(report[field]) for field in GLOBAL_TABLE_OBJECTIVE_FIELDS}


def global_table_objective_report(
    result_root: Path,
    evidence: DevelopmentEvidence,
    rounds: Sequence[RoundSpec],
    cohort: Sequence[int],
    final_survivor_record: Mapping[str, Any],
) -> Dict[str, int]:
    """Join the immutable final selection to every declared dev cell.

    Selection is never repeated or substituted here.  Instead, the supplied
    final survivor seal is first reproduced byte-for-byte from the same
    evidence.  A K with no retained R4 preferred attempt is joined to its
    matched control.  A retained attempt uses its direct physical observation
    or its exact no-op control alias at each width.  Adverse axes are the exact
    ``(BlockBytes, schedule, loss)`` slices, aggregated across all K and every
    declared root, matching the frozen per-K ranking definition.
    """
    result_root = result_root.resolve()
    verify_frozen_controller_runtime(result_root)
    validate_rounds(rounds)
    if not isinstance(final_survivor_record, dict):
        die("global table objective requires a canonical final survivor object")
    final_spec = rounds[-1]
    if (final_spec.name != "r4" or final_spec.retain_max != 1 or
            not final_spec.require_control_improvement):
        die("global table objective requires the frozen final R4 policy")
    selected = survivors_from_record(
        final_survivor_record, final_spec.name, cohort, final_spec.retain_max)
    rows = final_survivor_record.get("K")
    if not isinstance(rows, list):  # Also proven by survivors_from_record.
        die("final R4 survivor rows are unavailable")
    incoming = {
        int(row["K"]): tuple(row["input"])
        for row in rows
    }
    reproduced = analyze_round(
        result_root, evidence, rounds, final_spec, cohort, incoming,
        str(final_survivor_record["parent_survivor_sha256"]))
    if canonical_json_bytes(reproduced) != canonical_json_bytes(
            final_survivor_record):
        die("final R4 survivor seal does not reproduce from development evidence")

    panels = declared_panels_through(rounds, final_spec.name)
    allowed = set(panels)
    failure_multiplicity: Dict[int, int] = {K: 0 for K in cohort}
    axis_deltas: Dict[Tuple[int, str, str], int] = {}
    errors = total_failures = introductions = 0
    heavy_shortfall_sum = 0
    binary_deficit_gt15_count = 0
    binary_deficit_max = 0
    binary_deficit_sum = 0
    common_success_xors = 0
    common_success_muladds = 0
    for K in cohort:
        retained = selected[K]
        if len(retained) > 1:
            die("final R4 selection contains more than one preferred attempt")
        preferred = retained[0] if retained else None
        for panel in panels:
            for schedule in evidence.schedules:
                for loss in evidence.losses:
                    cell = CellKey(
                        K, panel.width, panel.root_index, schedule, loss)
                    control = evidence.control(cell, allowed)
                    chosen = control if preferred is None else \
                        evidence.logical_candidate(preferred, cell, allowed)
                    errors += chosen.error
                    # Reproducing the final ledger proves chosen errors are
                    # zero, but retain explicit report accounting and reject
                    # any impossible result classification defensively.
                    chosen_failed = failure_bit(chosen)
                    control_failed = failure_bit(control)
                    failure_multiplicity[K] += chosen_failed
                    total_failures += chosen_failed
                    introductions += int(not control_failed and chosen_failed)
                    axis = (panel.width, schedule, loss)
                    axis_deltas[axis] = (
                        axis_deltas.get(axis, 0) +
                        chosen_failed - control_failed)
                    heavy_shortfall_sum += chosen.heavy_shortfall
                    binary_deficit_gt15_count += int(
                        chosen.binary_deficit > 15)
                    binary_deficit_max = max(
                        binary_deficit_max, chosen.binary_deficit)
                    binary_deficit_sum += chosen.binary_deficit
                    if control.success and chosen.success:
                        common_success_xors += chosen.block_xors
                        common_success_muladds += chosen.block_muladds
    positive_axes = [delta for delta in axis_deltas.values() if delta > 0]
    report = {
        "errors": errors,
        "max_K_failure_multiplicity": max(
            failure_multiplicity.values(), default=0),
        "multi_failure_K_count": sum(
            count >= 2 for count in failure_multiplicity.values()),
        "total_failures": total_failures,
        "introductions": introductions,
        "worst_adverse_axis_delta": max(positive_axes, default=0),
        "sum_positive_axis_deltas": sum(positive_axes),
        "heavy_shortfall_sum": heavy_shortfall_sum,
        "binary_deficit_gt15_count": binary_deficit_gt15_count,
        "binary_deficit_max": binary_deficit_max,
        "binary_deficit_sum": binary_deficit_sum,
        "common_success_xors": common_success_xors,
        "common_success_muladds": common_success_muladds,
    }
    return validate_global_table_objective_report(report)


def write_global_table_objective_report(
    result_root: Path,
    evidence: DevelopmentEvidence,
    rounds: Sequence[RoundSpec],
    cohort: Sequence[int],
    final_survivor_record: Mapping[str, Any],
) -> str:
    """Publish the exact canonical R4 whole-table report and SHA256 sidecar."""
    result_root = result_root.resolve()
    report = global_table_objective_report(
        result_root, evidence, rounds, cohort, final_survivor_record)
    path = result_root / "analysis" / "global_table_objective_report.json"
    return write_hashed_artifact(path, canonical_json_bytes(report))
