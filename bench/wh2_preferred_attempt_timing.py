#!/usr/bin/env python3
"""Frozen WH2 preferred-attempt timing runner and exact analyzer.

The controller owns campaign sealing and supplies the routed ``K -> attempt``
map.  This module owns deterministic timing sampling, strict parsing of the
``preferredtiming`` CLI, isolated paired-cycle execution with bounded retries,
and exact-rational acceptance.  It intentionally has no production-profile
write path.

Integration entry points are ``select_timing_sample``,
``inspect_linux_isolation``, ``run_or_resume_timing_panel``, and
``analyze_timing``.
The runner requires an evidence probe with
``start(spec, process_index, cycle_index)`` and ``finish(token)`` methods; the
probe must retain the raw thermal and performance-counter intervals and return
their SHA256 bindings in :class:`EnvironmentEvidence`.
"""

from __future__ import annotations

import sys

# Frozen campaign directories have an exact immutable inventory.  Disable
# import-cache writes before loading any local frozen helper module.
sys.dont_write_bytecode = True

from dataclasses import asdict, dataclass, replace
import errno
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import stat
import subprocess
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import wh2_rank_floor_two_anchor_allk as common


SAMPLE_DOMAIN = b"wirehair.wh2.h12-preferred-attempt.timing-sample.v1|"
TRACE_SEED_DOMAIN = b"wirehair.wh2.h12-preferred-attempt.timing-trace-seed.v1"
BANDS = ((4096, 16383), (16384, 32767), (32768, 49151), (49152, 64000))
WIDTHS = (64, 1280, 4096)
SCHEDULES = ("burst", "adversarial", "repair-only")
ORDER = "CTTCTCCT"
TIMED_CYCLES = (1, 2, 3)
MAX_RETRIES_PER_CYCLE = 2
MAX_PROCESS_ATTEMPTS = 1 + len(TIMED_CYCLES) * MAX_RETRIES_PER_CYCLE
MIN_SAMPLE = 5
MAX_SAMPLE = 32
MIN_EVICTION_BYTES = 256 * 1024 * 1024
CPU_LIMIT_MILLIC = 85_000
DIMM_LIMIT_MILLIC = 90_000
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")

CSV_FIELDS = (
    "N", "bb", "metric", "cycle", "slot", "arm", "attempt", "result",
    "elapsed_ns", "saturated", "cpu_before", "cpu_after", "minflt_delta",
    "majflt_delta", "inactivated", "binary_def", "heavy_gain", "block_xors",
    "block_muladds", "source_bytes", "intermediate_bytes",
)
PREAMBLE_KEYS = (
    "schema", "policy", "metric", "cycles", "order", "discard_cycle",
    "cycle_mode", "cycle_index", "overhead", "payload",
    "payload_alignment", "payload_prefaulted", "route_cache_sha256",
    "route_context_sha256", "trace_sha256",
)
PREDECLARED_CONTAMINATION = frozenset((
    "migration", "major-fault", "throttle", "aperf-mperf-excursion",
    "temperature-excursion", "telemetry-gap",
))
PANEL_RESULT_SCHEMA = (
    "wirehair.wh2.h12_preferred_attempt.timing_panel_result.v1")
PANEL_RESULT_NAME = "panel.json"
PANEL_RESULT_SIDECAR_NAME = "panel.json.sha256"


class TimingError(RuntimeError):
    """Raised when a frozen timing invariant is violated."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_json(value: Any) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ensure_ascii=True,
    ) + "\n").encode("ascii")


def _sha256_file(path: Path, *, require_unique: bool = True) -> str:
    try:
        return common.sha256_file(path, require_unique=require_unique)
    except common.CampaignError as exc:
        raise TimingError(str(exc)) from exc


def _stable_file_bytes(path: Path) -> bytes:
    try:
        return common.stable_bytes(path)
    except common.CampaignError as exc:
        raise TimingError(str(exc)) from exc


def _canonical_object_file(path: Path, name: str) -> Tuple[Dict[str, Any], bytes]:
    try:
        raw = _stable_file_bytes(path)
        value = json.loads(raw.decode("ascii"))
        canonical = _canonical_json(value)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TimingError("%s is not readable canonical JSON" % name) from exc
    if not isinstance(value, dict) or canonical != raw:
        raise TimingError("%s is not a canonical JSON object" % name)
    return value, raw


def _strict_uint(text: str, name: str, maximum: int = (1 << 64) - 1) -> int:
    if UINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical unsigned integer" % name)
    value = int(text, 10)
    if value > maximum:
        raise TimingError("%s exceeds its integer domain" % name)
    return value


def _strict_sint(
    text: str, name: str, minimum: int = -(1 << 63), maximum: int = (1 << 63) - 1,
) -> int:
    if SINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical signed integer" % name)
    value = int(text, 10)
    if not minimum <= value <= maximum:
        raise TimingError("%s exceeds its integer domain" % name)
    return value


def _require_sha256(value: str, name: str) -> None:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise TimingError("%s is not a lowercase SHA256" % name)


def _upper_median(values: Sequence[Fraction]) -> Fraction:
    if not values:
        raise TimingError("upper median requires at least one value")
    return sorted(values)[len(values) // 2]


def _fraction_record(value: Fraction) -> Dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def sample_priority(K: int) -> bytes:
    if not isinstance(K, int) or isinstance(K, bool) or not 4096 <= K <= 64000:
        raise TimingError("timing K is outside the routed domain")
    return hashlib.sha256(SAMPLE_DOMAIN + str(K).encode("ascii")).digest()


def select_timing_sample(routed_ks: Iterable[int]) -> Tuple[int, ...]:
    """Select the immutable S=min(32,routed count) set, returned K-ascending."""
    values = list(routed_ks)
    if any(not isinstance(K, int) or isinstance(K, bool) for K in values):
        raise TimingError("routed K values must be integers")
    if len(values) != len(set(values)):
        raise TimingError("routed K values contain a duplicate")
    if any(not 4096 <= K <= 64000 for K in values):
        raise TimingError("routed K is outside all timing bands")
    target = min(MAX_SAMPLE, len(values))
    if target < MIN_SAMPLE:
        raise TimingError("fewer than five routed K values cannot power timing")
    priority = lambda K: (sample_priority(K), K)
    chosen = set()
    for low, high in BANDS:
        band = sorted((K for K in values if low <= K <= high), key=priority)
        chosen.update(band[:8])
    if len(chosen) > target:
        raise TimingError("band selection exceeded the timing target")
    remaining = sorted((K for K in values if K not in chosen), key=priority)
    chosen.update(remaining[:target - len(chosen)])
    if len(chosen) != target:
        raise TimingError("global timing fill did not reach S")
    return tuple(sorted(chosen))


def timing_panel_seed(
    metric: str, K: int, block_bytes: int, schedule: str,
) -> int:
    if (not isinstance(metric, str) or metric not in ("solve", "setup") or
            not isinstance(K, int) or isinstance(K, bool) or
            not 4096 <= K <= 64000 or
            not isinstance(block_bytes, int) or isinstance(block_bytes, bool) or
            block_bytes not in WIDTHS or
            not isinstance(schedule, str) or schedule not in SCHEDULES):
        raise TimingError("timing trace-seed coordinates are outside the grid")
    payload = b"\0".join((
        TRACE_SEED_DOMAIN, metric.encode("ascii"), str(K).encode("ascii"),
        str(block_bytes).encode("ascii"), schedule.encode("ascii"),
    ))
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


@dataclass(frozen=True)
class TimingPanelSpec:
    K: int
    block_bytes: int
    preferred_attempt: int
    metric: str
    seed: int
    schedule: str

    def __post_init__(self) -> None:
        if (not isinstance(self.K, int) or isinstance(self.K, bool) or
                not 4096 <= self.K <= 64000):
            raise TimingError("panel K is outside the timing domain")
        if (not isinstance(self.block_bytes, int) or
                isinstance(self.block_bytes, bool) or
                self.block_bytes not in WIDTHS):
            raise TimingError("panel width is not frozen")
        if (not isinstance(self.preferred_attempt, int) or
                isinstance(self.preferred_attempt, bool) or
                not 0 <= self.preferred_attempt <= 255 or
                (self.block_bytes == 64 and self.preferred_attempt == 0)):
            raise TimingError("preferred attempt is outside its eligible domain")
        if (not isinstance(self.metric, str) or
                self.metric not in ("solve", "setup")):
            raise TimingError("timing metric must be solve or setup")
        if (not isinstance(self.seed, int) or isinstance(self.seed, bool) or
                not 0 <= self.seed < 1 << 64):
            raise TimingError("timing seed is outside u64")
        if (not isinstance(self.schedule, str) or
                self.schedule not in SCHEDULES):
            raise TimingError("timing schedule is not frozen")
        if self.seed != timing_panel_seed(
                self.metric, self.K, self.block_bytes, self.schedule):
            raise TimingError("timing panel seed does not match frozen derivation")

    def key(self) -> Tuple[str, int, int, str]:
        return (self.metric, self.K, self.block_bytes, self.schedule)


@dataclass(frozen=True)
class TimingRow:
    K: int
    block_bytes: int
    metric: str
    cycle: int
    slot: int
    arm: str
    attempt: int
    result: int
    elapsed_ns: int
    saturated: int
    cpu_before: int
    cpu_after: int
    minflt_delta: int
    majflt_delta: int
    inactivated: int
    binary_def: int
    heavy_gain: int
    block_xors: int
    block_muladds: int
    source_bytes: int
    intermediate_bytes: int

    def work_signature(self) -> Tuple[int, ...]:
        return (
            self.result, self.inactivated, self.binary_def, self.heavy_gain,
            self.block_xors, self.block_muladds, self.source_bytes,
            self.intermediate_bytes,
        )


@dataclass(frozen=True)
class ParsedTimingOutput:
    spec: TimingPanelSpec
    cycle_index: Optional[int]
    canonical_attempt: int
    trace_sha256: str
    stdout_sha256: str
    rows: Tuple[TimingRow, ...]

    def cycle_rows(self, cycle: int) -> Tuple[TimingRow, ...]:
        result = tuple(row for row in self.rows if row.cycle == cycle)
        if len(result) != 8:
            raise TimingError("timing cycle does not contain eight rows")
        return result


def _parse_preamble(line: str) -> Tuple[Tuple[str, str], ...]:
    prefix = "# preferredtiming: "
    if not line.startswith(prefix):
        raise TimingError("preferredtiming preamble prefix mismatch")
    tokens = line[len(prefix):].split(" ")
    pairs: List[Tuple[str, str]] = []
    for token in tokens:
        if token.count("=") != 1:
            raise TimingError("preferredtiming preamble token is malformed")
        key, value = token.split("=", 1)
        if not key or not value:
            raise TimingError("preferredtiming preamble token is empty")
        pairs.append((key, value))
    if tuple(key for key, _value in pairs) != PREAMBLE_KEYS:
        raise TimingError("preferredtiming preamble key order/schema mismatch")
    return tuple(pairs)


def parse_timing_output(
    raw: bytes,
    spec: TimingPanelSpec,
    route_cache_sha256: str,
    route_context_sha256: str,
    cycle_index: Optional[int] = None,
) -> ParsedTimingOutput:
    """Strictly bind a full 32-row run or one indexed 8-row replacement."""
    _require_sha256(route_cache_sha256, "route-cache hash")
    _require_sha256(route_context_sha256, "route-context hash")
    if (cycle_index is not None and
            (not isinstance(cycle_index, int) or isinstance(cycle_index, bool) or
             not 0 <= cycle_index <= 3)):
        raise TimingError("replacement cycle index is outside 0..3")
    if not raw or not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw:
        raise TimingError("preferredtiming output is not canonical LF text")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TimingError("preferredtiming output is not ASCII") from exc
    lines = text.splitlines()
    expected_line_count = 34 if cycle_index is None else 10
    if len(lines) != expected_line_count:
        raise TimingError(
            "preferredtiming output has the wrong full/replacement line count")
    preamble_pairs = _parse_preamble(lines[0])
    preamble = dict(preamble_pairs)
    expected_metric = {
        "schema": "v1", "policy": "h12-q0-adaptive", "metric": spec.metric,
        "cycles": "4" if cycle_index is None else "1",
        "order": ORDER, "discard_cycle": "0",
        "cycle_mode": "full" if cycle_index is None else "replacement",
        "cycle_index": "all" if cycle_index is None else str(cycle_index),
        "overhead": "4" if spec.metric == "solve" else "none",
        "payload": "distinct-zero-v1" if spec.metric == "solve" else "none",
        "payload_alignment": "64" if spec.metric == "solve" else "none",
        "payload_prefaulted": "1" if spec.metric == "solve" else "none",
        "route_cache_sha256": route_cache_sha256,
        "route_context_sha256": route_context_sha256,
    }
    if any(preamble.get(key) != value for key, value in expected_metric.items()):
        raise TimingError("preferredtiming preamble binding mismatch")
    trace_sha256 = preamble["trace_sha256"]
    _require_sha256(trace_sha256, "trace hash")
    if tuple(lines[1].split(",")) != CSV_FIELDS:
        raise TimingError("preferredtiming CSV header mismatch")

    rows: List[TimingRow] = []
    control_attempt: Optional[int] = None
    for index, line in enumerate(lines[2:]):
        fields = line.split(",")
        if len(fields) != len(CSV_FIELDS) or any(value == "" for value in fields):
            raise TimingError("preferredtiming row field-count mismatch")
        row = dict(zip(CSV_FIELDS, fields))
        emitted_cycle, slot = divmod(index, 8)
        cycle = emitted_cycle if cycle_index is None else cycle_index
        arm = "control" if ORDER[slot] == "C" else "candidate"
        parsed = TimingRow(
            K=_strict_uint(row["N"], "N", (1 << 32) - 1),
            block_bytes=_strict_uint(row["bb"], "bb", (1 << 32) - 1),
            metric=row["metric"],
            cycle=_strict_uint(row["cycle"], "cycle", 3),
            slot=_strict_uint(row["slot"], "slot", 7),
            arm=row["arm"],
            attempt=_strict_uint(row["attempt"], "attempt", 255),
            result=_strict_sint(row["result"], "result", -(1 << 31), (1 << 31) - 1),
            elapsed_ns=_strict_uint(row["elapsed_ns"], "elapsed_ns"),
            saturated=_strict_uint(row["saturated"], "saturated", 1),
            cpu_before=_strict_sint(row["cpu_before"], "cpu_before", -1, (1 << 31) - 1),
            cpu_after=_strict_sint(row["cpu_after"], "cpu_after", -1, (1 << 31) - 1),
            minflt_delta=_strict_sint(row["minflt_delta"], "minflt_delta"),
            majflt_delta=_strict_sint(row["majflt_delta"], "majflt_delta"),
            inactivated=_strict_uint(row["inactivated"], "inactivated", (1 << 32) - 1),
            binary_def=_strict_uint(row["binary_def"], "binary_def", (1 << 32) - 1),
            heavy_gain=_strict_uint(row["heavy_gain"], "heavy_gain", (1 << 32) - 1),
            block_xors=_strict_uint(row["block_xors"], "block_xors"),
            block_muladds=_strict_uint(row["block_muladds"], "block_muladds"),
            source_bytes=_strict_uint(row["source_bytes"], "source_bytes"),
            intermediate_bytes=_strict_uint(
                row["intermediate_bytes"], "intermediate_bytes"),
        )
        if (parsed.K != spec.K or parsed.block_bytes != spec.block_bytes or
                parsed.metric != spec.metric or parsed.cycle != cycle or
                parsed.slot != slot or parsed.arm != arm):
            raise TimingError("preferredtiming row order/binding mismatch")
        if parsed.result != 0 or parsed.elapsed_ns == 0 or parsed.saturated != 0:
            raise TimingError("preferredtiming row is not a clean successful timing")
        if arm == "candidate":
            if parsed.attempt != spec.preferred_attempt:
                raise TimingError("candidate attempt changed")
        elif control_attempt is None:
            control_attempt = parsed.attempt
        elif parsed.attempt != control_attempt:
            raise TimingError("control attempt changed")
        if spec.metric == "setup":
            if parsed.work_signature() != (0,) * 8:
                raise TimingError("setup row contains solve work")
        elif (parsed.source_bytes != spec.K * spec.block_bytes or
              parsed.intermediate_bytes <= parsed.source_bytes or
              parsed.intermediate_bytes % spec.block_bytes != 0):
            raise TimingError("solve source/intermediate byte count mismatch")
        rows.append(parsed)
    if (control_attempt is None or control_attempt > spec.preferred_attempt or
            (spec.block_bytes == 64 and
             control_attempt == spec.preferred_attempt)):
        raise TimingError("control/preferred route is not eligible and ordered")
    for arm in ("control", "candidate"):
        signatures = {row.work_signature() for row in rows if row.arm == arm}
        if len(signatures) != 1:
            raise TimingError("deterministic work counters changed within an arm")
    if control_attempt == spec.preferred_attempt and len({
            row.work_signature() for row in rows}) != 1:
        raise TimingError("neutral alias changed work between timing arms")
    if (len({row.source_bytes for row in rows}) != 1 or
            len({row.intermediate_bytes for row in rows}) != 1):
        raise TimingError("source/intermediate byte counts changed across arms")
    return ParsedTimingOutput(
        spec=spec, cycle_index=cycle_index, canonical_attempt=control_attempt,
        trace_sha256=trace_sha256, stdout_sha256=_sha256(raw), rows=tuple(rows),
    )


@dataclass(frozen=True)
class EnvironmentEvidence:
    core: int
    numa_node: int
    exclusive_core: bool
    load_workers_stopped: bool
    sibling_busy_ticks: int
    cpu_temperature_millic: int
    dimm_temperature_millic: int
    dimm_read_errors: int
    edac_ce_delta: int
    edac_ue_delta: int
    throttled: bool
    aperf_mperf_excursion: bool
    telemetry_gap: bool
    thermal_interval_sha256: str
    performance_interval_sha256: str

    def __post_init__(self) -> None:
        integer_fields = (
            self.core, self.numa_node, self.sibling_busy_ticks,
            self.cpu_temperature_millic, self.dimm_temperature_millic,
            self.dimm_read_errors, self.edac_ce_delta, self.edac_ue_delta,
        )
        if any(not isinstance(value, int) or isinstance(value, bool)
               for value in integer_fields):
            raise TimingError("environment evidence integer type mismatch")
        if any(value < 0 for value in integer_fields):
            raise TimingError("environment evidence contains a negative counter")
        for value in (
                self.exclusive_core, self.load_workers_stopped, self.throttled,
                self.aperf_mperf_excursion, self.telemetry_gap):
            if not isinstance(value, bool):
                raise TimingError("environment evidence boolean type mismatch")
        if (not self.telemetry_gap and
                (self.cpu_temperature_millic == 0 or
                 self.dimm_temperature_millic == 0)):
            raise TimingError("complete telemetry contains a zero temperature")
        _require_sha256(self.thermal_interval_sha256, "thermal interval hash")
        _require_sha256(
            self.performance_interval_sha256, "performance interval hash")

    def record(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class LinuxIsolation:
    core: int
    sibling_cpus: Tuple[int, ...]
    numa_node: int
    llc_bytes: int
    eviction_bytes: int

    def __post_init__(self) -> None:
        if (not isinstance(self.core, int) or isinstance(self.core, bool) or
                self.core < 0 or
                not isinstance(self.numa_node, int) or
                isinstance(self.numa_node, bool) or self.numa_node < 0 or
                not isinstance(self.llc_bytes, int) or
                isinstance(self.llc_bytes, bool) or self.llc_bytes <= 0 or
                not isinstance(self.eviction_bytes, int) or
                isinstance(self.eviction_bytes, bool) or
                self.eviction_bytes <= 0):
            raise TimingError("Linux timing isolation has an invalid integer")
        if (not isinstance(self.sibling_cpus, tuple) or
                any(not isinstance(cpu, int) or isinstance(cpu, bool) or cpu < 0
                    for cpu in self.sibling_cpus) or
                tuple(sorted(set(self.sibling_cpus))) != self.sibling_cpus or
                self.core in self.sibling_cpus):
            raise TimingError("Linux timing sibling topology is not canonical")
        if self.eviction_bytes != max(2 * self.llc_bytes, MIN_EVICTION_BYTES):
            raise TimingError("Linux timing eviction size is not frozen")


def _parse_cpu_list(text: str) -> Tuple[int, ...]:
    if not text or text.strip() != text:
        raise TimingError("CPU-list text is not canonical")
    values: List[int] = []
    for token in text.split(","):
        if "-" in token:
            fields = token.split("-")
            if len(fields) != 2:
                raise TimingError("CPU-list range is malformed")
            low = _strict_uint(fields[0], "CPU-list low", (1 << 31) - 1)
            high = _strict_uint(fields[1], "CPU-list high", (1 << 31) - 1)
            if low >= high:
                raise TimingError("CPU-list range is not increasing")
            values.extend(range(low, high + 1))
        else:
            values.append(_strict_uint(token, "CPU-list value", (1 << 31) - 1))
    if values != sorted(set(values)):
        raise TimingError("CPU-list is duplicate or nonascending")
    return tuple(values)


def _parse_cache_size(text: str) -> int:
    match = re.fullmatch(r"([1-9][0-9]*)([KMG])", text)
    if match is None:
        raise TimingError("Linux cache size is not canonical")
    scale = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3}[match.group(2)]
    return int(match.group(1)) * scale


def inspect_linux_isolation(
    core: int, sysfs_root: Path = Path("/sys/devices/system/cpu"),
) -> LinuxIsolation:
    if not isinstance(core, int) or isinstance(core, bool) or core < 0:
        raise TimingError("timing core is invalid")
    if not isinstance(sysfs_root, Path):
        raise TimingError("CPU sysfs root must be pathlib.Path")
    cpu = sysfs_root / ("cpu%d" % core)
    if cpu.is_symlink() or not cpu.is_dir():
        raise TimingError("timing CPU sysfs directory is unavailable")
    online = cpu / "online"
    if online.exists() and online.read_text(encoding="ascii").strip() != "1":
        raise TimingError("timing CPU is not online")
    siblings = _parse_cpu_list(
        (cpu / "topology/thread_siblings_list").read_text(encoding="ascii").strip())
    if core not in siblings:
        raise TimingError("timing core is absent from its sibling topology")
    nodes = []
    for path in cpu.glob("node[0-9]*"):
        match = re.fullmatch(r"node([0-9]+)", path.name)
        if match is not None:
            nodes.append(int(match.group(1)))
    if len(set(nodes)) != 1:
        raise TimingError("timing core does not have one NUMA node")
    cache_records: List[Tuple[int, int]] = []
    for index in (cpu / "cache").glob("index[0-9]*"):
        cache_type = (index / "type").read_text(encoding="ascii").strip()
        if cache_type not in ("Data", "Unified"):
            continue
        level = _strict_uint(
            (index / "level").read_text(encoding="ascii").strip(),
            "cache level", 32)
        shared = _parse_cpu_list(
            (index / "shared_cpu_list").read_text(encoding="ascii").strip())
        if core not in shared:
            raise TimingError("cache record does not contain timing core")
        size = _parse_cache_size(
            (index / "size").read_text(encoding="ascii").strip())
        cache_records.append((level, size))
    if not cache_records:
        raise TimingError("timing core has no data/unified cache record")
    maximum_level = max(level for level, _size in cache_records)
    llc_bytes = max(size for level, size in cache_records if level == maximum_level)
    eviction = max(2 * llc_bytes, MIN_EVICTION_BYTES)
    return LinuxIsolation(
        core=core, sibling_cpus=tuple(value for value in siblings if value != core),
        numa_node=nodes[0], llc_bytes=llc_bytes, eviction_bytes=eviction,
    )


def _parse_launcher_schema(
    launcher: Tuple[str, ...],
) -> Tuple[int, int, Tuple[Tuple[int, str], ...]]:
    if (len(launcher) == 6 and
            Path(launcher[0]).name == "numactl" and
            Path(launcher[0]).is_absolute() and
            launcher[1].startswith("--physcpubind=") and
            launcher[2].startswith("--membind=") and
            Path(launcher[3]).name == "taskset" and
            Path(launcher[3]).is_absolute() and
            launcher[4] == "-c"):
        core_text = launcher[1][len("--physcpubind="):]
        node_text = launcher[2][len("--membind="):]
        taskset_core = launcher[5]
        executables = ((0, "numactl"), (3, "taskset"))
    elif (len(launcher) == 15 and
          Path(launcher[0]).name == "sudo" and
          Path(launcher[0]).is_absolute() and launcher[1] == "-n" and
          Path(launcher[2]).name == "systemd-run" and
          Path(launcher[2]).is_absolute() and
          launcher[3:7] == (
              "--scope", "--quiet", "--collect",
              "--slice=wirehair-wh2-timing.slice") and
          launcher[7].startswith("--property=AllowedCPUs=") and
          launcher[8] == "--" and
          Path(launcher[9]).name == "numactl" and
          Path(launcher[9]).is_absolute() and
          launcher[10].startswith("--physcpubind=") and
          launcher[11].startswith("--membind=") and
          Path(launcher[12]).name == "taskset" and
          Path(launcher[12]).is_absolute() and launcher[13] == "-c"):
        core_text = launcher[10][len("--physcpubind="):]
        node_text = launcher[11][len("--membind="):]
        taskset_core = launcher[14]
        allowed_core = launcher[7][len("--property=AllowedCPUs="):]
        if allowed_core != core_text:
            raise TimingError("systemd timing AllowedCPUs binding changed")
        executables = (
            (0, "sudo"), (2, "systemd-run"),
            (9, "numactl"), (12, "taskset"),
        )
    else:
        raise TimingError("timing launcher schema is not frozen")
    core = _strict_uint(core_text, "launcher physical CPU", (1 << 31) - 1)
    node = _strict_uint(node_text, "launcher NUMA node", (1 << 31) - 1)
    if taskset_core != str(core):
        raise TimingError("timing launcher taskset CPU changed")
    return core, node, executables


@dataclass(frozen=True)
class TimingRunnerConfig:
    binary: Path
    binary_sha256: str
    route_cache: Path
    route_cache_sha256: str
    route_context_sha256: str
    isolation: LinuxIsolation
    launcher: Tuple[str, ...]
    timeout_seconds: float = 1800.0
    launcher_sha256: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.binary, Path) or not isinstance(self.route_cache, Path):
            raise TimingError("timing binary/cache paths must be pathlib.Path")
        if not isinstance(self.isolation, LinuxIsolation):
            raise TimingError("timing isolation is not a LinuxIsolation record")
        _require_sha256(self.binary_sha256, "binary hash")
        _require_sha256(self.route_cache_sha256, "route-cache hash")
        _require_sha256(self.route_context_sha256, "route-context hash")
        if (not isinstance(self.timeout_seconds, (int, float)) or
                isinstance(self.timeout_seconds, bool) or
                not math.isfinite(float(self.timeout_seconds)) or
                self.timeout_seconds <= 0):
            raise TimingError("timing timeout must be finite and positive")
        if self.isolation.eviction_bytes != max(
                2 * self.isolation.llc_bytes, MIN_EVICTION_BYTES):
            raise TimingError("eviction bytes do not equal max(2*LLC,256MiB)")
        if (not isinstance(self.launcher, tuple) or
                any(not isinstance(value, str) or not value
                    for value in self.launcher)):
            raise TimingError("timing launcher is not a canonical argument tuple")
        core, node, executables = _parse_launcher_schema(self.launcher)
        if core != self.isolation.core or node != self.isolation.numa_node:
            raise TimingError(
                "timing launcher does not bind NUMA and taskset isolation")
        if not isinstance(self.launcher_sha256, tuple):
            raise TimingError("timing launcher hashes are not a tuple")
        if not self.launcher_sha256:
            hashes = []
            for index, name in executables:
                path = Path(self.launcher[index])
                try:
                    resolved = path.resolve(strict=True)
                    digest = _sha256_file(path, require_unique=False)
                except (OSError, TimingError) as exc:
                    raise TimingError(
                        "timing launcher executable is unavailable: %s" %
                        name) from exc
                if (path.is_symlink() or resolved != path or
                        not path.is_file() or
                        not os.access(str(path), os.X_OK)):
                    raise TimingError(
                        "timing launcher executable is not frozen: %s" % name)
                hashes.append(digest)
            object.__setattr__(self, "launcher_sha256", tuple(hashes))
        elif len(self.launcher_sha256) != len(executables):
            raise TimingError("timing launcher hash count changed")
        for value in self.launcher_sha256:
            _require_sha256(value, "timing launcher executable hash")


def _verify_launcher_files(config: TimingRunnerConfig) -> None:
    core, node, executables = _parse_launcher_schema(config.launcher)
    if core != config.isolation.core or node != config.isolation.numa_node:
        raise TimingError("timing launcher isolation binding changed")
    if len(config.launcher_sha256) != len(executables):
        raise TimingError("timing launcher hash count changed")
    for (index, name), expected_sha256 in zip(
            executables, config.launcher_sha256):
        path = Path(config.launcher[index])
        try:
            resolved = path.resolve(strict=True)
            digest = _sha256_file(path, require_unique=False)
        except (OSError, TimingError) as exc:
            raise TimingError(
                "timing launcher executable is unavailable: %s" % name) from exc
        if (path.is_symlink() or resolved != path or
                not path.is_file() or not os.access(str(path), os.X_OK) or
                path.name != name or digest != expected_sha256):
            raise TimingError("timing launcher executable is not frozen: %s" % name)


def _verified_timing_inputs(
    config: TimingRunnerConfig,
) -> Tuple[Path, Path]:
    _verify_launcher_files(config)
    if config.binary.is_symlink() or config.route_cache.is_symlink():
        raise TimingError("timing binary/cache must not be symlinks")
    try:
        binary = config.binary.resolve(strict=True)
        route_cache = config.route_cache.resolve(strict=True)
    except OSError as exc:
        raise TimingError("timing binary/cache is unavailable") from exc
    if not os.access(str(binary), os.X_OK):
        raise TimingError("timing binary is not executable")
    try:
        if _sha256_file(config.binary) != config.binary_sha256:
            raise TimingError("timing binary hash mismatch")
        if _sha256_file(config.route_cache) != config.route_cache_sha256:
            raise TimingError("timing route-cache hash mismatch")
    except OSError as exc:
        raise TimingError("timing binary/cache cannot be hashed") from exc
    return binary, route_cache


def linux_runner_config(
    binary: Path,
    binary_sha256: str,
    route_cache: Path,
    route_cache_sha256: str,
    route_context_sha256: str,
    core: int,
    timeout_seconds: float = 1800.0,
) -> TimingRunnerConfig:
    isolation = inspect_linux_isolation(core)
    numactl = shutil.which("numactl")
    taskset = shutil.which("taskset")
    if numactl is None or taskset is None:
        raise TimingError("timing isolation requires numactl and taskset")
    launcher = (
        str(Path(numactl).resolve()), "--physcpubind=%d" % core,
        "--membind=%d" % isolation.numa_node,
        str(Path(taskset).resolve()), "-c", str(core),
    )
    return TimingRunnerConfig(
        binary=binary, binary_sha256=binary_sha256,
        route_cache=route_cache, route_cache_sha256=route_cache_sha256,
        route_context_sha256=route_context_sha256, isolation=isolation,
        launcher=launcher, timeout_seconds=timeout_seconds,
    )


@dataclass(frozen=True)
class TimingCycle:
    cycle: int
    attempt_index: int
    rows: Tuple[TimingRow, ...]

    def ratio(self) -> Fraction:
        if (not isinstance(self.cycle, int) or isinstance(self.cycle, bool) or
                not 0 <= self.cycle <= 3 or
                not isinstance(self.attempt_index, int) or
                isinstance(self.attempt_index, bool) or
                not 0 <= self.attempt_index < MAX_PROCESS_ATTEMPTS or
                not isinstance(self.rows, tuple) or
                any(not isinstance(row, TimingRow) for row in self.rows) or
                len(self.rows) != 8 or
                any(row.cycle != self.cycle for row in self.rows) or
                tuple(row.slot for row in self.rows) != tuple(range(8)) or
                tuple(row.arm for row in self.rows) != tuple(
                    "control" if arm == "C" else "candidate"
                    for arm in ORDER)):
            raise TimingError("timing cycle row layout is invalid")
        candidate = sum(row.elapsed_ns for row in self.rows if row.arm == "candidate")
        control = sum(row.elapsed_ns for row in self.rows if row.arm == "control")
        if candidate <= 0 or control <= 0:
            raise TimingError("timing cycle has a nonpositive arm sum")
        return Fraction(candidate, control)


COMMAND_OPTIONS = (
    "--N", "--bb", "--preferred-attempt", "--evict-bytes", "--metric",
    "--route-cache", "--route-cache-sha256", "--route-context-sha256",
    "--loss", "--seed", "--schedule",
)


def _attempt_command_options(
    command: Tuple[str, ...],
) -> Tuple[Tuple[str, ...], Dict[str, str]]:
    if command.count("preferredtiming") != 1:
        raise TimingError("attempt command has no unique preferredtiming verb")
    marker = command.index("preferredtiming")
    if marker < 1:
        raise TimingError("attempt command omits its executable")
    arguments = command[marker + 1:]
    if len(arguments) % 2 != 0:
        raise TimingError("attempt command option/value layout is malformed")
    names = tuple(arguments[::2])
    expected_names = COMMAND_OPTIONS + (
        ("--cycle-index",) if "--cycle-index" in names else ())
    if names != expected_names:
        raise TimingError("attempt command option order/schema mismatch")
    values = tuple(arguments[1::2])
    if any(not value for value in values):
        raise TimingError("attempt command contains an empty option value")
    return command[:marker], dict(zip(names, values))


def _validate_attempt_isolation_prefix(
    prefix: Tuple[str, ...], evidence: EnvironmentEvidence,
) -> int:
    if len(prefix) < 2 or not Path(prefix[-1]).is_absolute():
        raise TimingError(
            "attempt command does not bind its measured NUMA/core isolation")
    core, node, _executables = _parse_launcher_schema(prefix[:-1])
    if node != evidence.numa_node:
        raise TimingError(
            "attempt command does not bind its measured NUMA/core isolation")
    return core


@dataclass(frozen=True)
class AttemptRecord:
    attempt_index: int
    cycle_index: Optional[int]
    command: Tuple[str, ...]
    returncode: int
    stdout_sha256: str
    stderr_sha256: str
    trace_sha256: str
    environment: EnvironmentEvidence
    cycle_contamination: Mapping[int, Tuple[str, ...]]
    accepted_cycles: Tuple[int, ...]
    record_sha256: str

    def core_record(self) -> Dict[str, Any]:
        return {
            "attempt_index": self.attempt_index,
            "cycle_index": self.cycle_index,
            "command": list(self.command), "returncode": self.returncode,
            "stdout_sha256": self.stdout_sha256,
            "stderr_sha256": self.stderr_sha256,
            "trace_sha256": self.trace_sha256,
            "environment": self.environment.record(),
            "environment_sha256": _sha256(_canonical_json(self.environment.record())),
            "cycle_contamination": {
                str(key): list(value)
                for key, value in sorted(self.cycle_contamination.items())
            },
            "accepted_cycles": list(self.accepted_cycles),
        }

    def validate(self) -> None:
        if (not isinstance(self.attempt_index, int) or
                isinstance(self.attempt_index, bool) or
                not 0 <= self.attempt_index < MAX_PROCESS_ATTEMPTS):
            raise TimingError("attempt index is invalid")
        if (self.cycle_index is not None and
                (not isinstance(self.cycle_index, int) or
                 isinstance(self.cycle_index, bool) or
                 self.cycle_index not in TIMED_CYCLES)):
            raise TimingError("attempt replacement cycle index is invalid")
        if (not isinstance(self.command, tuple) or not self.command or
                any(not isinstance(value, str) or not value
                    for value in self.command)):
            raise TimingError("attempt command is invalid")
        if (not isinstance(self.cycle_contamination, Mapping) or
                not isinstance(self.accepted_cycles, tuple)):
            raise TimingError("attempt cycle evidence containers are invalid")
        prefix, options = _attempt_command_options(self.command)
        if not isinstance(self.environment, EnvironmentEvidence):
            raise TimingError("attempt environment evidence is invalid")
        commanded_core = _validate_attempt_isolation_prefix(
            prefix, self.environment)
        replacement_text = options.get("--cycle-index")
        if ((self.cycle_index is None and replacement_text is not None) or
                (self.cycle_index is not None and
                 replacement_text != str(self.cycle_index))):
            raise TimingError("attempt command replacement index mismatch")
        _strict_uint(options["--evict-bytes"], "command eviction bytes")
        if int(options["--evict-bytes"]) < MIN_EVICTION_BYTES:
            raise TimingError("attempt command eviction is below 256MiB")
        _require_sha256(options["--route-cache-sha256"], "command route hash")
        _require_sha256(
            options["--route-context-sha256"], "command route-context hash")
        if options["--loss"] != "0.50":
            raise TimingError("attempt command loss is not frozen")
        if (not isinstance(self.returncode, int) or
                isinstance(self.returncode, bool) or self.returncode != 0):
            raise TimingError("attempt return code is invalid")
        if self.stderr_sha256 != _sha256(b""):
            raise TimingError("successful attempt retained nonempty stderr")
        if (not self.environment.exclusive_core or
                not self.environment.load_workers_stopped or
                self.environment.sibling_busy_ticks != 0):
            raise TimingError("successful attempt retained invalid isolation")
        expected_cycles = (
            set(range(4)) if self.cycle_index is None else {self.cycle_index})
        if set(self.cycle_contamination) != expected_cycles:
            raise TimingError("attempt contamination coverage is invalid")
        for cycle, reasons in self.cycle_contamination.items():
            if (not isinstance(cycle, int) or isinstance(cycle, bool) or
                    not isinstance(reasons, tuple) or
                    tuple(sorted(set(reasons))) != reasons or
                    not set(reasons).issubset(PREDECLARED_CONTAMINATION)):
                raise TimingError("attempt contamination is not canonical")
        if (self.environment.core != commanded_core and
                any("migration" not in reasons
                    for reasons in self.cycle_contamination.values())):
            raise TimingError("attempt omitted commanded-core migration")
        environment_reasons = set(
            _environment_contamination(self.environment))
        if any(not environment_reasons.issubset(reasons)
               for reasons in self.cycle_contamination.values()):
            raise TimingError("attempt omitted environment contamination")
        if (tuple(sorted(set(self.accepted_cycles))) != self.accepted_cycles or
                any(not isinstance(cycle, int) or isinstance(cycle, bool)
                    for cycle in self.accepted_cycles) or
                not set(self.accepted_cycles).issubset(TIMED_CYCLES) or
                any(self.cycle_contamination[cycle]
                    for cycle in self.accepted_cycles)):
            raise TimingError("attempt accepted-cycle set is invalid")
        for value, name in (
                (self.stdout_sha256, "attempt stdout hash"),
                (self.stderr_sha256, "attempt stderr hash"),
                (self.trace_sha256, "attempt trace hash"),
                (self.record_sha256, "attempt record hash")):
            _require_sha256(value, name)
        if _sha256(_canonical_json(self.core_record())) != self.record_sha256:
            raise TimingError("attempt record hash mismatch")


@dataclass(frozen=True)
class TimingPanelResult:
    spec: TimingPanelSpec
    trace_sha256: str
    canonical_attempt: int
    cycles: Tuple[TimingCycle, ...]
    attempts: Tuple[AttemptRecord, ...]

    def validate(self) -> None:
        if (not isinstance(self.spec, TimingPanelSpec) or
                not isinstance(self.cycles, tuple) or
                any(not isinstance(cycle, TimingCycle) for cycle in self.cycles) or
                not isinstance(self.attempts, tuple) or
                any(not isinstance(record, AttemptRecord)
                    for record in self.attempts)):
            raise TimingError("timing panel result container types are invalid")
        _require_sha256(self.trace_sha256, "panel trace hash")
        if (not isinstance(self.canonical_attempt, int) or
                isinstance(self.canonical_attempt, bool) or
                not 0 <= self.canonical_attempt <= self.spec.preferred_attempt or
                (self.spec.block_bytes == 64 and
                 self.canonical_attempt == self.spec.preferred_attempt)):
            raise TimingError("panel control attempt is not ordered")
        if not 1 <= len(self.attempts) <= MAX_PROCESS_ATTEMPTS:
            raise TimingError("panel retained too many process attempts")
        if tuple(record.attempt_index for record in self.attempts) != tuple(
                range(len(self.attempts))):
            raise TimingError("panel attempt records are not contiguous")
        runtime_binding: Optional[Tuple[Any, ...]] = None
        for record in self.attempts:
            record.validate()
            if record.trace_sha256 != self.trace_sha256:
                raise TimingError("panel attempt trace hash changed")
            if (record.environment.dimm_read_errors != 0 or
                    record.environment.edac_ce_delta != 0 or
                    record.environment.edac_ue_delta != 0):
                raise TimingError("panel contains EDAC or DIMM read errors")
            prefix, options = _attempt_command_options(record.command)
            if (options["--N"] != str(self.spec.K) or
                    options["--bb"] != str(self.spec.block_bytes) or
                    options["--preferred-attempt"] !=
                    str(self.spec.preferred_attempt) or
                    options["--metric"] != self.spec.metric or
                    options["--seed"] != str(self.spec.seed) or
                    options["--schedule"] != self.spec.schedule):
                raise TimingError("panel attempt command does not bind its spec")
            binding = (
                prefix, options["--evict-bytes"], options["--route-cache"],
                options["--route-cache-sha256"],
                options["--route-context-sha256"],
            )
            if runtime_binding is None:
                runtime_binding = binding
            elif binding != runtime_binding:
                raise TimingError("panel runtime binding changed across retries")
        if self.attempts[0].cycle_index is not None:
            raise TimingError("panel first attempt is not a full four-cycle run")
        if any(record.cycle_index is None for record in self.attempts[1:]):
            raise TimingError("panel repeated a full run instead of one dirty cycle")
        pending = set(TIMED_CYCLES)
        full = self.attempts[0]
        full_accepted = tuple(
            cycle for cycle in TIMED_CYCLES
            if not full.cycle_contamination[cycle])
        if full.accepted_cycles != full_accepted:
            raise TimingError("panel full attempt accepted the wrong cycle set")
        pending.difference_update(full_accepted)
        attempt_index = 1
        for _retry_round in range(MAX_RETRIES_PER_CYCLE):
            targets = tuple(sorted(pending))
            for target in targets:
                if attempt_index >= len(self.attempts):
                    break
                record = self.attempts[attempt_index]
                if record.cycle_index != target:
                    raise TimingError("panel replacement retry order changed")
                expected_accepted = (
                    (target,) if not record.cycle_contamination[target] else ())
                if record.accepted_cycles != expected_accepted:
                    raise TimingError(
                        "panel replacement accepted the wrong cycle set")
                pending.difference_update(expected_accepted)
                attempt_index += 1
            if not pending:
                break
        if attempt_index != len(self.attempts):
            raise TimingError("panel retained an extra cycle retry")
        if pending:
            raise TimingError("panel attempts leave contaminated timing cycles")
        if tuple(cycle.cycle for cycle in self.cycles) != TIMED_CYCLES:
            raise TimingError("panel does not contain timed cycles 1,2,3")
        accepted = {
            (record.attempt_index, cycle)
            for record in self.attempts for cycle in record.accepted_cycles
        }
        if len(accepted) != len(TIMED_CYCLES):
            raise TimingError("panel accepted cycles are duplicate or incomplete")
        for cycle in self.cycles:
            if (cycle.attempt_index, cycle.cycle) not in accepted:
                raise TimingError("panel cycle is not backed by its attempt record")
            if any(row.K != self.spec.K or
                   row.block_bytes != self.spec.block_bytes or
                   row.metric != self.spec.metric
                   for row in cycle.rows):
                raise TimingError("panel cycle rows do not bind to the panel spec")
            record = self.attempts[cycle.attempt_index]
            for row in cycle.rows:
                if (row.result != 0 or row.saturated != 0 or
                        row.elapsed_ns <= 0 or row.minflt_delta < 0 or
                        row.majflt_delta != 0 or
                        row.cpu_before != record.environment.core or
                        row.cpu_after != record.environment.core):
                    raise TimingError("accepted panel cycle is contaminated")
                expected_attempt = (
                    self.canonical_attempt if row.arm == "control" else
                    self.spec.preferred_attempt)
                if row.attempt != expected_attempt:
                    raise TimingError("accepted panel cycle attempt changed")
                if self.spec.metric == "setup":
                    if row.work_signature() != (0,) * 8:
                        raise TimingError("accepted setup cycle contains solve work")
                elif (row.source_bytes != self.spec.K * self.spec.block_bytes or
                      row.intermediate_bytes <= row.source_bytes or
                      row.intermediate_bytes % self.spec.block_bytes != 0):
                    raise TimingError("accepted solve cycle byte counts are invalid")
            cycle.ratio()
        for arm in ("control", "candidate"):
            signatures = {
                row.work_signature() for cycle in self.cycles
                for row in cycle.rows if row.arm == arm
            }
            if len(signatures) != 1:
                raise TimingError("accepted cycle work changed across retries")
        if self.canonical_attempt == self.spec.preferred_attempt and len({
                row.work_signature() for cycle in self.cycles
                for row in cycle.rows}) != 1:
            raise TimingError("accepted neutral alias work differs between arms")
        all_rows = [row for cycle in self.cycles for row in cycle.rows]
        if (len({row.source_bytes for row in all_rows}) != 1 or
                len({row.intermediate_bytes for row in all_rows}) != 1):
            raise TimingError("accepted cycle bytes changed across arms")

    def panel_ratio(self) -> Fraction:
        self.validate()
        return _upper_median([cycle.ratio() for cycle in self.cycles])


def _environment_contamination(
    evidence: EnvironmentEvidence,
) -> Tuple[str, ...]:
    reasons = set()
    if (evidence.cpu_temperature_millic >= CPU_LIMIT_MILLIC or
            evidence.dimm_temperature_millic >= DIMM_LIMIT_MILLIC):
        reasons.add("temperature-excursion")
    if evidence.throttled:
        reasons.add("throttle")
    if evidence.aperf_mperf_excursion:
        reasons.add("aperf-mperf-excursion")
    if evidence.telemetry_gap:
        reasons.add("telemetry-gap")
    if not reasons.issubset(PREDECLARED_CONTAMINATION):
        raise TimingError("unrecognized timing contamination")
    return tuple(sorted(reasons))


def _global_contamination(
    evidence: EnvironmentEvidence, config: TimingRunnerConfig,
) -> Tuple[str, ...]:
    reasons = set(_environment_contamination(evidence))
    if evidence.core != config.isolation.core:
        reasons.add("migration")
    return tuple(sorted(reasons))


def _isolation_failure(
    evidence: EnvironmentEvidence, config: TimingRunnerConfig,
) -> Optional[str]:
    failures = []
    if not evidence.exclusive_core:
        failures.append("exclusive physical core was lost")
    if not evidence.load_workers_stopped:
        failures.append("load workers were not stopped")
    if evidence.sibling_busy_ticks != 0:
        failures.append("SMT sibling was busy")
    if evidence.numa_node != config.isolation.numa_node:
        failures.append("NUMA placement changed")
    return "; ".join(failures) if failures else None


def _cycle_contamination(
    parsed: ParsedTimingOutput,
    cycle: int,
    evidence: EnvironmentEvidence,
    config: TimingRunnerConfig,
) -> Tuple[str, ...]:
    reasons = set(_global_contamination(evidence, config))
    for row in parsed.cycle_rows(cycle):
        if (row.cpu_before != config.isolation.core or
                row.cpu_after != config.isolation.core or
                row.cpu_before != row.cpu_after):
            reasons.add("migration")
        if row.majflt_delta < 0 or row.minflt_delta < 0:
            reasons.add("telemetry-gap")
        elif row.majflt_delta != 0:
            reasons.add("major-fault")
    return tuple(sorted(reasons))


def _discard_stale_atomic_partials(path: Path) -> bool:
    try:
        return common._discard_stale_atomic_partials(path)
    except common.CampaignError as exc:
        raise TimingError(str(exc)) from exc


def _fsync_directory(path: Path) -> None:
    directory_fd = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _atomic_write(path: Path, data: bytes) -> None:
    """Durably publish immutable timing bytes without target replacement."""
    try:
        common.atomic_write_once_or_same(path, data)
    except common.CampaignError as exc:
        raise TimingError(str(exc)) from exc


def _build_attempt_record(
    index: int,
    cycle_index: Optional[int],
    command: Tuple[str, ...],
    completed: subprocess.CompletedProcess,
    parsed: ParsedTimingOutput,
    evidence: EnvironmentEvidence,
    contamination: Mapping[int, Tuple[str, ...]],
    accepted: Tuple[int, ...],
) -> AttemptRecord:
    provisional = AttemptRecord(
        attempt_index=index, cycle_index=cycle_index,
        command=command, returncode=completed.returncode,
        stdout_sha256=_sha256(completed.stdout),
        stderr_sha256=_sha256(completed.stderr),
        trace_sha256=parsed.trace_sha256, environment=evidence,
        cycle_contamination=dict(contamination), accepted_cycles=accepted,
        record_sha256="0" * 64,
    )
    return replace(
        provisional,
        record_sha256=_sha256(_canonical_json(provisional.core_record())),
    )


def _persist_attempt(
    directory: Path,
    record: AttemptRecord,
    stdout: bytes,
    stderr: bytes,
) -> None:
    stem = "attempt-%02d" % record.attempt_index
    _atomic_write(directory / (stem + ".stdout"), stdout)
    _atomic_write(directory / (stem + ".stderr"), stderr)
    payload = record.core_record()
    payload["record_sha256"] = record.record_sha256
    _atomic_write(directory / (stem + ".json"), _canonical_json(payload))


def _persist_failed_attempt(
    directory: Path,
    index: int,
    cycle_index: Optional[int],
    command: Tuple[str, ...],
    stdout: bytes,
    stderr: bytes,
    returncode: Optional[int],
    evidence: Optional[EnvironmentEvidence],
    failure: str,
) -> None:
    """Retain and hash an attempted invocation that cannot enter analysis."""
    if not isinstance(failure, str) or not failure:
        raise TimingError("failed-attempt reason is empty")
    payload: Dict[str, Any] = {
        "schema": "wirehair.wh2.h12_preferred_attempt.failed_attempt.v1",
        "attempt_index": index,
        "cycle_index": cycle_index,
        "command": list(command),
        "returncode": returncode,
        "stdout_sha256": _sha256(stdout),
        "stderr_sha256": _sha256(stderr),
        "failure": failure,
        "environment": evidence.record() if evidence is not None else None,
        "environment_sha256": (
            _sha256(_canonical_json(evidence.record()))
            if evidence is not None else None),
    }
    payload["record_sha256"] = _sha256(_canonical_json(payload))
    stem = "attempt-%02d" % index
    _atomic_write(directory / (stem + ".stdout"), stdout)
    _atomic_write(directory / (stem + ".stderr"), stderr)
    _atomic_write(directory / (stem + ".json"), _canonical_json(payload))


def _timing_command(
    spec: TimingPanelSpec,
    config: TimingRunnerConfig,
    binary: Optional[Path] = None,
    route_cache: Optional[Path] = None,
    cycle_index: Optional[int] = None,
) -> Tuple[str, ...]:
    binary = config.binary if binary is None else binary
    route_cache = config.route_cache if route_cache is None else route_cache
    command = config.launcher + (
        str(binary), "preferredtiming", "--N", str(spec.K),
        "--bb", str(spec.block_bytes), "--preferred-attempt",
        str(spec.preferred_attempt), "--evict-bytes",
        str(config.isolation.eviction_bytes), "--metric", spec.metric,
        "--route-cache", str(route_cache), "--route-cache-sha256",
        config.route_cache_sha256, "--route-context-sha256",
        config.route_context_sha256, "--loss", "0.50", "--seed",
        str(spec.seed), "--schedule", spec.schedule,
    )
    if cycle_index is not None:
        command += ("--cycle-index", str(cycle_index))
    return command


ATTEMPT_RECORD_FIELDS = frozenset((
    "attempt_index", "cycle_index", "command", "returncode",
    "stdout_sha256", "stderr_sha256", "trace_sha256", "environment",
    "environment_sha256", "cycle_contamination", "accepted_cycles",
    "record_sha256",
))
PANEL_RESULT_FIELDS = frozenset((
    "schema", "spec", "trace_sha256", "canonical_attempt",
    "attempt_count", "attempt_record_sha256", "accepted_cycles",
    "runtime_binding", "self_sha256_excluding_field",
))
PANEL_SPEC_FIELDS = frozenset((
    "K", "block_bytes", "preferred_attempt", "metric", "seed", "schedule",
))
PANEL_RUNTIME_FIELDS = frozenset((
    "binary", "binary_sha256", "route_cache", "route_cache_sha256",
    "route_context_sha256", "launcher", "launcher_sha256", "core",
    "numa_node", "llc_bytes", "eviction_bytes", "sibling_cpus",
    "timeout_seconds",
))


def _attempt_file_paths(
    directory: Path, index: int,
) -> Tuple[Path, Path, Path]:
    stem = "attempt-%02d" % index
    return (
        directory / (stem + ".json"),
        directory / (stem + ".stdout"),
        directory / (stem + ".stderr"),
    )


def _attempt_inventory(count: int) -> set[str]:
    names = set()
    for index in range(count):
        names.update(path.name for path in _attempt_file_paths(Path("."), index))
    return names


def _directory_inventory(directory: Path) -> set[str]:
    try:
        return {path.name for path in directory.iterdir()}
    except OSError as exc:
        raise TimingError("timing panel directory cannot be inventoried") from exc


def _decode_attempt_record(path: Path) -> AttemptRecord:
    payload, _raw = _canonical_object_file(path, "timing attempt record")
    if set(payload) != ATTEMPT_RECORD_FIELDS:
        raise TimingError("timing attempt record fields do not match schema")
    environment_value = payload.get("environment")
    environment_fields = frozenset(EnvironmentEvidence.__dataclass_fields__)
    if (not isinstance(environment_value, dict) or
            set(environment_value) != environment_fields):
        raise TimingError("timing attempt environment fields do not match schema")
    try:
        environment = EnvironmentEvidence(**environment_value)
    except TypeError as exc:
        raise TimingError("timing attempt environment is malformed") from exc
    environment_sha256 = payload.get("environment_sha256")
    _require_sha256(environment_sha256, "attempt environment hash")
    if environment_sha256 != _sha256(_canonical_json(environment.record())):
        raise TimingError("timing attempt environment hash mismatch")

    command_value = payload.get("command")
    if (not isinstance(command_value, list) or
            any(not isinstance(value, str) for value in command_value)):
        raise TimingError("timing attempt command is malformed")
    contamination_value = payload.get("cycle_contamination")
    if not isinstance(contamination_value, dict):
        raise TimingError("timing attempt contamination is malformed")
    contamination: Dict[int, Tuple[str, ...]] = {}
    for key, reasons in contamination_value.items():
        if not isinstance(key, str):
            raise TimingError("timing attempt contamination key is malformed")
        cycle = _strict_uint(key, "attempt contamination cycle", 3)
        if (not isinstance(reasons, list) or
                any(not isinstance(reason, str) for reason in reasons)):
            raise TimingError("timing attempt contamination reasons are malformed")
        if cycle in contamination:
            raise TimingError("timing attempt contamination cycle is duplicated")
        contamination[cycle] = tuple(reasons)
    accepted_value = payload.get("accepted_cycles")
    if not isinstance(accepted_value, list):
        raise TimingError("timing attempt accepted cycles are malformed")
    record = AttemptRecord(
        attempt_index=payload.get("attempt_index"),
        cycle_index=payload.get("cycle_index"),
        command=tuple(command_value), returncode=payload.get("returncode"),
        stdout_sha256=payload.get("stdout_sha256"),
        stderr_sha256=payload.get("stderr_sha256"),
        trace_sha256=payload.get("trace_sha256"), environment=environment,
        cycle_contamination=contamination,
        accepted_cycles=tuple(accepted_value),
        record_sha256=payload.get("record_sha256"),
    )
    record.validate()
    expected = record.core_record()
    expected["record_sha256"] = record.record_sha256
    if payload != expected:
        raise TimingError("timing attempt record is not its canonical value")
    return record


def _load_retained_attempt(
    directory: Path,
    index: int,
    spec: TimingPanelSpec,
    config: TimingRunnerConfig,
    binary: Path,
    route_cache: Path,
) -> Tuple[AttemptRecord, ParsedTimingOutput]:
    record_path, stdout_path, stderr_path = _attempt_file_paths(directory, index)
    record = _decode_attempt_record(record_path)
    if record.attempt_index != index:
        raise TimingError("retained timing attempt index changed")
    try:
        stdout = _stable_file_bytes(stdout_path)
        stderr = _stable_file_bytes(stderr_path)
    except OSError as exc:
        raise TimingError("retained timing attempt stream is unavailable") from exc
    if (_sha256(stdout) != record.stdout_sha256 or
            _sha256(stderr) != record.stderr_sha256):
        raise TimingError("retained timing attempt stream hash mismatch")
    expected_command = _timing_command(
        spec, config, binary, route_cache, record.cycle_index)
    if record.command != expected_command:
        raise TimingError("retained timing attempt command binding changed")
    parsed = parse_timing_output(
        stdout, spec, config.route_cache_sha256,
        config.route_context_sha256, record.cycle_index)
    if (parsed.stdout_sha256 != record.stdout_sha256 or
            parsed.trace_sha256 != record.trace_sha256):
        raise TimingError("retained timing attempt output binding changed")
    emitted_cycles = (
        tuple(range(4)) if record.cycle_index is None else
        (record.cycle_index,))
    contamination = {
        cycle: _cycle_contamination(parsed, cycle, record.environment, config)
        for cycle in emitted_cycles
    }
    if contamination != dict(record.cycle_contamination):
        raise TimingError("retained timing attempt contamination changed")
    if (record.environment.edac_ce_delta != 0 or
            record.environment.edac_ue_delta != 0 or
            record.environment.dimm_read_errors != 0):
        raise TimingError("retained timing attempt contains EDAC/DIMM errors")
    isolation_failure = _isolation_failure(record.environment, config)
    if isolation_failure is not None:
        raise TimingError(
            "retained timing attempt isolation failed: %s" % isolation_failure)
    return record, parsed


def _panel_spec_record(spec: TimingPanelSpec) -> Dict[str, Any]:
    return {
        "K": spec.K, "block_bytes": spec.block_bytes,
        "preferred_attempt": spec.preferred_attempt, "metric": spec.metric,
        "seed": spec.seed, "schedule": spec.schedule,
    }


def _runtime_binding_record(
    config: TimingRunnerConfig, binary: Path, route_cache: Path,
) -> Dict[str, Any]:
    return {
        "binary": str(binary), "binary_sha256": config.binary_sha256,
        "route_cache": str(route_cache),
        "route_cache_sha256": config.route_cache_sha256,
        "route_context_sha256": config.route_context_sha256,
        "launcher": list(config.launcher),
        "launcher_sha256": list(config.launcher_sha256),
        "core": config.isolation.core,
        "numa_node": config.isolation.numa_node,
        "sibling_cpus": list(config.isolation.sibling_cpus),
        "llc_bytes": config.isolation.llc_bytes,
        "eviction_bytes": config.isolation.eviction_bytes,
        "timeout_seconds": config.timeout_seconds,
    }


def _panel_result_record(
    result: TimingPanelResult,
    config: TimingRunnerConfig,
    binary: Path,
    route_cache: Path,
) -> Dict[str, Any]:
    core: Dict[str, Any] = {
        "schema": PANEL_RESULT_SCHEMA,
        "spec": _panel_spec_record(result.spec),
        "trace_sha256": result.trace_sha256,
        "canonical_attempt": result.canonical_attempt,
        "attempt_count": len(result.attempts),
        "attempt_record_sha256": [
            record.record_sha256 for record in result.attempts],
        "runtime_binding": _runtime_binding_record(
            config, binary, route_cache),
        "accepted_cycles": [
            {"cycle": cycle.cycle, "attempt_index": cycle.attempt_index}
            for cycle in result.cycles
        ],
    }
    record = dict(core)
    record["self_sha256_excluding_field"] = _sha256(_canonical_json(core))
    return record


def write_timing_panel_result(
    attempt_directory: Path,
    result: TimingPanelResult,
    config: TimingRunnerConfig,
) -> str:
    """Seal one successful panel after verifying every retained attempt file."""
    if not isinstance(attempt_directory, Path):
        raise TimingError("timing panel directory must be pathlib.Path")
    if (attempt_directory.is_symlink() or
            not attempt_directory.is_dir()):
        raise TimingError("timing panel directory is unavailable or a symlink")
    result.validate()
    binary, route_cache = _verified_timing_inputs(config)
    expected_inventory = _attempt_inventory(len(result.attempts))
    if _directory_inventory(attempt_directory) != expected_inventory:
        raise TimingError("successful timing panel attempt inventory is not exact")
    parsed_by_attempt: Dict[int, ParsedTimingOutput] = {}
    for index, expected_record in enumerate(result.attempts):
        record, parsed = _load_retained_attempt(
            attempt_directory, index, result.spec, config, binary, route_cache)
        if record != expected_record:
            raise TimingError("retained timing attempt differs from panel result")
        parsed_by_attempt[index] = parsed
    for cycle in result.cycles:
        if cycle.rows != parsed_by_attempt[cycle.attempt_index].cycle_rows(
                cycle.cycle):
            raise TimingError("panel cycle differs from retained timing stdout")

    record = _panel_result_record(result, config, binary, route_cache)
    encoded = _canonical_json(record)
    digest = _sha256(encoded)
    manifest_path = attempt_directory / PANEL_RESULT_NAME
    sidecar_path = attempt_directory / PANEL_RESULT_SIDECAR_NAME
    if (manifest_path.exists() or manifest_path.is_symlink() or
            sidecar_path.exists() or sidecar_path.is_symlink()):
        raise TimingError("timing panel result seal already exists")
    _atomic_write(manifest_path, encoded)
    _atomic_write(
        sidecar_path,
        (digest + "  " + PANEL_RESULT_NAME + "\n").encode("ascii"))
    loaded = load_timing_panel_result(attempt_directory, result.spec, config)
    if loaded != result:
        raise TimingError("published timing panel result did not round-trip")
    return digest


def _load_timing_panel_result(
    attempt_directory: Path,
    expected_spec: TimingPanelSpec,
    config: TimingRunnerConfig,
    *,
    require_sidecar: bool,
) -> Tuple[TimingPanelResult, bytes]:
    """Validate and recompute one panel, optionally before sidecar publication."""
    if not isinstance(attempt_directory, Path):
        raise TimingError("timing panel directory must be pathlib.Path")
    if (attempt_directory.is_symlink() or
            not attempt_directory.is_dir()):
        raise TimingError("timing panel directory is unavailable or a symlink")
    binary, route_cache = _verified_timing_inputs(config)
    manifest_path = attempt_directory / PANEL_RESULT_NAME
    sidecar_path = attempt_directory / PANEL_RESULT_SIDECAR_NAME
    manifest_present = manifest_path.exists() or manifest_path.is_symlink()
    sidecar_present = sidecar_path.exists() or sidecar_path.is_symlink()
    if not manifest_present or (require_sidecar and not sidecar_present):
        raise TimingError(
            "timing panel is partial or failed; controller-level rerun is forbidden")
    manifest, encoded = _canonical_object_file(
        manifest_path, "timing panel result")
    digest = _sha256(encoded)
    if sidecar_present:
        try:
            sidecar = _stable_file_bytes(sidecar_path)
        except OSError as exc:
            raise TimingError("timing panel result sidecar is unavailable") from exc
        if sidecar != (
                digest + "  " + PANEL_RESULT_NAME + "\n").encode("ascii"):
            raise TimingError("timing panel result file hash mismatch")
    if set(manifest) != PANEL_RESULT_FIELDS:
        raise TimingError("timing panel result fields do not match schema")
    if manifest.get("schema") != PANEL_RESULT_SCHEMA:
        raise TimingError("timing panel result schema mismatch")
    self_sha256 = manifest.get("self_sha256_excluding_field")
    _require_sha256(self_sha256, "timing panel result self hash")
    unsigned = dict(manifest)
    del unsigned["self_sha256_excluding_field"]
    if _sha256(_canonical_json(unsigned)) != self_sha256:
        raise TimingError("timing panel result self hash mismatch")

    spec_value = manifest.get("spec")
    if not isinstance(spec_value, dict) or set(spec_value) != PANEL_SPEC_FIELDS:
        raise TimingError("timing panel result spec fields do not match schema")
    try:
        spec = TimingPanelSpec(**spec_value)
    except TypeError as exc:
        raise TimingError("timing panel result spec is malformed") from exc
    if spec != expected_spec:
        raise TimingError("timing panel result does not match expected spec")
    runtime_binding = manifest.get("runtime_binding")
    expected_runtime_binding = _runtime_binding_record(
        config, binary, route_cache)
    if (not isinstance(runtime_binding, dict) or
            set(runtime_binding) != PANEL_RUNTIME_FIELDS or
            _canonical_json(runtime_binding) !=
            _canonical_json(expected_runtime_binding)):
        raise TimingError("timing panel result runtime binding changed")
    attempt_count = manifest.get("attempt_count")
    if (not isinstance(attempt_count, int) or isinstance(attempt_count, bool) or
            not 1 <= attempt_count <= MAX_PROCESS_ATTEMPTS):
        raise TimingError("timing panel result attempt count is invalid")
    attempt_hashes = manifest.get("attempt_record_sha256")
    if (not isinstance(attempt_hashes, list) or
            len(attempt_hashes) != attempt_count):
        raise TimingError("timing panel result attempt hashes are malformed")
    for value in attempt_hashes:
        _require_sha256(value, "timing panel attempt record hash")
    expected_inventory = _attempt_inventory(attempt_count)
    expected_inventory.add(PANEL_RESULT_NAME)
    if sidecar_present:
        expected_inventory.add(PANEL_RESULT_SIDECAR_NAME)
    if _directory_inventory(attempt_directory) != expected_inventory:
        raise TimingError("sealed timing panel file inventory is not exact")

    records: List[AttemptRecord] = []
    parsed_by_attempt: Dict[int, ParsedTimingOutput] = {}
    for index in range(attempt_count):
        record, parsed = _load_retained_attempt(
            attempt_directory, index, spec, config, binary, route_cache)
        if record.record_sha256 != attempt_hashes[index]:
            raise TimingError("timing panel attempt record hash changed")
        records.append(record)
        parsed_by_attempt[index] = parsed
    trace_sha256 = manifest.get("trace_sha256")
    _require_sha256(trace_sha256, "timing panel trace hash")
    canonical_attempt = manifest.get("canonical_attempt")
    if any(parsed.trace_sha256 != trace_sha256 or
           parsed.canonical_attempt != canonical_attempt
           for parsed in parsed_by_attempt.values()):
        raise TimingError("timing panel output trace or route changed")

    accepted: Dict[int, TimingCycle] = {}
    for record in records:
        parsed = parsed_by_attempt[record.attempt_index]
        for cycle in record.accepted_cycles:
            if cycle in accepted:
                raise TimingError("timing panel accepted a cycle twice")
            accepted[cycle] = TimingCycle(
                cycle=cycle, attempt_index=record.attempt_index,
                rows=parsed.cycle_rows(cycle))
    if set(accepted) != set(TIMED_CYCLES):
        raise TimingError("timing panel accepted-cycle coverage is incomplete")
    cycles = tuple(accepted[cycle] for cycle in TIMED_CYCLES)
    accepted_value = manifest.get("accepted_cycles")
    expected_accepted = [
        {"cycle": cycle.cycle, "attempt_index": cycle.attempt_index}
        for cycle in cycles
    ]
    if (not isinstance(accepted_value, list) or
            _canonical_json(accepted_value) != _canonical_json(expected_accepted)):
        raise TimingError("timing panel accepted-cycle manifest changed")
    result = TimingPanelResult(
        spec=spec, trace_sha256=trace_sha256,
        canonical_attempt=canonical_attempt, cycles=cycles,
        attempts=tuple(records),
    )
    result.validate()
    if _canonical_json(_panel_result_record(
            result, config, binary, route_cache)) != encoded:
        raise TimingError("timing panel result is not the exact recomputed record")
    return result, encoded


def load_timing_panel_result(
    attempt_directory: Path,
    expected_spec: TimingPanelSpec,
    config: TimingRunnerConfig,
) -> TimingPanelResult:
    """Load only a completely sealed panel; failed/partial panels always reject."""
    result, _encoded = _load_timing_panel_result(
        attempt_directory, expected_spec, config, require_sidecar=True)
    return result


def _execute_timing_process(
    command: Tuple[str, ...],
    *,
    stdout: Any,
    stderr: Any,
    check: bool,
    timeout: float,
) -> subprocess.CompletedProcess:
    """Execute one timing scope and reap its entire local process group."""
    if check:
        raise TimingError("timing executor must use explicit return-code handling")
    process = subprocess.Popen(
        command, stdout=stdout, stderr=stderr, start_new_session=True)
    try:
        output, errors = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        output = exc.output if exc.output is not None else b""
        errors = exc.stderr if exc.stderr is not None else b""
        try:
            # Never perform an unbounded drain here: an escaped session could
            # retain our pipe descriptors after the original group is dead.
            common.stop_and_reap_process_group(process, 1.0)
        except common.CampaignError as cleanup_error:
            raise TimingError(
                "timed-out process group cleanup could not be proven") \
                from cleanup_error
        exc.stdout = output
        exc.output = output
        exc.stderr = errors
        raise
    except BaseException as error:
        try:
            common.stop_and_reap_process_group(process, 1.0)
        except common.CampaignError as cleanup_error:
            raise TimingError(
                "interrupted process group cleanup could not be proven") \
                from cleanup_error
        raise
    descendants_live = common.process_group_exists(process)
    if descendants_live:
        try:
            common.stop_and_reap_process_group(process)
        except common.CampaignError as exc:
            raise TimingError(
                "timing process group cleanup could not be proven") from exc
        raise TimingError("timing process left an unexpected child process")
    common.close_process_streams(process)
    return subprocess.CompletedProcess(
        command, process.returncode, stdout=output, stderr=errors)


def run_timing_panel(
    spec: TimingPanelSpec,
    config: TimingRunnerConfig,
    evidence_probe: Any,
    attempt_directory: Path,
    executor: Callable[..., subprocess.CompletedProcess] = _execute_timing_process,
) -> TimingPanelResult:
    """Run a full panel, then retry each dirty paired cycle at most twice."""
    if not isinstance(attempt_directory, Path):
        raise TimingError("timing panel directory must be pathlib.Path")
    binary, route_cache = _verified_timing_inputs(config)
    try:
        directory_fd = common.create_durable_directory(attempt_directory)
    except common.CampaignError as exc:
        raise TimingError(str(exc)) from exc
    os.close(directory_fd)
    pending = set(TIMED_CYCLES)
    accepted: Dict[int, TimingCycle] = {}
    records: List[AttemptRecord] = []
    expected_trace: Optional[str] = None
    expected_control: Optional[int] = None
    expected_work: Optional[Dict[str, Tuple[int, ...]]] = None
    process_index = 0

    for retry_round in range(MAX_RETRIES_PER_CYCLE + 1):
        if not pending:
            break
        cycle_targets: Tuple[Optional[int], ...] = (
            (None,) if retry_round == 0 else
            tuple(sorted(pending)))
        for cycle_index in cycle_targets:
            if cycle_index is not None and cycle_index not in pending:
                continue
            try:
                _verify_launcher_files(config)
                inputs_changed = (
                    _sha256_file(config.binary) != config.binary_sha256 or
                    _sha256_file(config.route_cache) !=
                    config.route_cache_sha256)
            except (OSError, TimingError) as exc:
                raise TimingError(
                    "timing launcher, executable, or route cache became unavailable "
                    "between retries") from exc
            if inputs_changed:
                raise TimingError(
                    "timing executable or route cache changed between retries")
            command = _timing_command(
                spec, config, binary, route_cache, cycle_index)
            try:
                token = evidence_probe.start(spec, process_index, cycle_index)
            except Exception as exc:
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    b"", b"", None, None,
                    "environment probe start failed: %s: %s" % (
                        type(exc).__name__, exc))
                raise TimingError("timing environment probe failed to start") from exc
            completed: Optional[subprocess.CompletedProcess] = None
            execution_failure: Optional[BaseException] = None
            evidence: Optional[EnvironmentEvidence] = None
            evidence_failure: Optional[BaseException] = None
            try:
                completed = executor(
                    command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=False, timeout=config.timeout_seconds,
                )
            except Exception as exc:
                execution_failure = exc
            finally:
                try:
                    evidence = evidence_probe.finish(token)
                except Exception as exc:
                    evidence_failure = exc

            output_owner = (
                completed if completed is not None else execution_failure)
            stdout_value = getattr(output_owner, "stdout", b"") or b""
            stderr_value = getattr(output_owner, "stderr", b"") or b""
            stdout = stdout_value if isinstance(stdout_value, bytes) else b""
            stderr = stderr_value if isinstance(stderr_value, bytes) else b""
            returncode = getattr(completed, "returncode", None)
            if (evidence_failure is not None or
                    not isinstance(evidence, EnvironmentEvidence)):
                reason = (
                    "%s: %s" % (type(evidence_failure).__name__, evidence_failure)
                    if evidence_failure is not None else
                    "probe returned non-EnvironmentEvidence")
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    stdout, stderr, returncode, None,
                    "environment evidence failed: %s" % reason)
                raise TimingError(
                    "timing environment evidence failed") from evidence_failure
            input_failure: Optional[str] = None
            try:
                _verify_launcher_files(config)
                if _sha256_file(config.binary) != config.binary_sha256:
                    input_failure = "timing binary changed during invocation"
                elif _sha256_file(config.route_cache) != config.route_cache_sha256:
                    input_failure = "timing route cache changed during invocation"
            except (OSError, TimingError) as exc:
                input_failure = (
                    "timing launcher/input became unreadable or changed: %s" % exc)
            if input_failure is not None:
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    stdout, stderr, returncode, evidence, input_failure)
                raise TimingError(input_failure)
            if execution_failure is not None:
                failure = "%s: %s" % (
                    type(execution_failure).__name__, str(execution_failure))
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    stdout, stderr, None, evidence, failure)
                raise TimingError(
                    "preferredtiming process execution failed") from execution_failure
            if completed is None:
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    b"", b"", None, evidence,
                    "executor returned no process result")
                raise TimingError("timing executor returned no process result")
            if not isinstance(completed, subprocess.CompletedProcess):
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    stdout, stderr, returncode, evidence,
                    "executor returned a non-CompletedProcess result")
                raise TimingError("timing executor returned an invalid process result")
            if (not isinstance(completed.stdout, bytes) or
                    not isinstance(completed.stderr, bytes)):
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    b"", b"", completed.returncode, evidence,
                    "executor returned non-byte stdout or stderr")
                raise TimingError("timing executor must return byte streams")
            if (not isinstance(completed.returncode, int) or
                    isinstance(completed.returncode, bool) or
                    completed.returncode != 0 or completed.stderr):
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    completed.stdout, completed.stderr, completed.returncode,
                    evidence, "nonzero return code or nonempty stderr")
                raise TimingError("preferredtiming process failed or wrote stderr")
            try:
                parsed = parse_timing_output(
                    completed.stdout, spec, config.route_cache_sha256,
                    config.route_context_sha256, cycle_index)
            except TimingError as exc:
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    completed.stdout, completed.stderr, completed.returncode,
                    evidence, "strict output parse failed: %s" % exc)
                raise

            work = {
                arm: next(iter({
                    row.work_signature() for row in parsed.rows
                    if row.arm == arm
                }))
                for arm in ("control", "candidate")
            }
            if expected_trace is None:
                if cycle_index is not None:
                    raise TimingError("first timing attempt was not full")
                expected_trace = parsed.trace_sha256
                expected_control = parsed.canonical_attempt
                expected_work = work
            elif (parsed.trace_sha256 != expected_trace or
                  parsed.canonical_attempt != expected_control or
                  work != expected_work):
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    completed.stdout, completed.stderr, completed.returncode,
                    evidence, "timing trace, route, or work changed across retries")
                raise TimingError("timing trace/route/work changed across retries")

            emitted_cycles = (
                tuple(range(4)) if cycle_index is None else (cycle_index,))
            contamination = {
                cycle: _cycle_contamination(parsed, cycle, evidence, config)
                for cycle in emitted_cycles
            }
            if (evidence.edac_ce_delta != 0 or evidence.edac_ue_delta != 0 or
                    evidence.dimm_read_errors != 0):
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    completed.stdout, completed.stderr, completed.returncode,
                    evidence, "EDAC or DIMM telemetry changed during timing")
                raise TimingError("EDAC or DIMM telemetry changed during timing")

            isolation_failure = _isolation_failure(evidence, config)
            if isolation_failure is not None:
                _persist_failed_attempt(
                    attempt_directory, process_index, cycle_index, command,
                    completed.stdout, completed.stderr, completed.returncode,
                    evidence, "timing isolation failed: %s" % isolation_failure)
                raise TimingError("timing isolation failed: %s" % isolation_failure)

            newly_accepted = tuple(
                cycle for cycle in emitted_cycles
                if cycle in pending and not contamination[cycle])
            record = _build_attempt_record(
                process_index, cycle_index, command, completed, parsed,
                evidence, contamination, newly_accepted)
            record.validate()
            _persist_attempt(
                attempt_directory, record, completed.stdout, completed.stderr)
            records.append(record)
            for cycle in newly_accepted:
                accepted[cycle] = TimingCycle(
                    cycle=cycle, attempt_index=process_index,
                    rows=parsed.cycle_rows(cycle))
                pending.remove(cycle)
            process_index += 1
    if pending:
        raise TimingError(
            "timing contamination survived one original plus two retries: %s" %
            sorted(pending))
    result = TimingPanelResult(
        spec=spec, trace_sha256=expected_trace or "",
        canonical_attempt=expected_control if expected_control is not None else -1,
        cycles=tuple(accepted[cycle] for cycle in TIMED_CYCLES),
        attempts=tuple(records),
    )
    result.panel_ratio()
    write_timing_panel_result(attempt_directory, result, config)
    return result


def run_or_resume_timing_panel(
    spec: TimingPanelSpec,
    config: TimingRunnerConfig,
    evidence_probe: Any,
    attempt_directory: Path,
    executor: Callable[..., subprocess.CompletedProcess] = _execute_timing_process,
) -> TimingPanelResult:
    """Resume a clean seal, or run only when no panel directory exists.

    The frozen protocol has no controller-level retry after a failed panel.
    The sole recoverable partial state is a fully published and independently
    recomputable ``panel.json`` whose checksum sidecar was not yet published
    when the controller stopped.  The caller holds the timing-host lock while
    that sidecar is repaired; no process attempt is invoked or repeated.
    """
    if not isinstance(attempt_directory, Path):
        raise TimingError("timing panel directory must be pathlib.Path")
    if attempt_directory.exists() or attempt_directory.is_symlink():
        flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_DIRECTORY", 0) |
                 getattr(os, "O_NOFOLLOW", 0))
        try:
            directory_fd = os.open(attempt_directory, flags)
        except OSError as exc:
            raise TimingError(
                "timing panel must be a plain directory") from exc
        try:
            metadata = os.fstat(directory_fd)
            pathname_metadata = os.stat(
                attempt_directory, follow_symlinks=False)
        except OSError as exc:
            raise TimingError(
                "timing panel directory identity is unstable") from exc
        finally:
            os.close(directory_fd)
        if (not stat.S_ISDIR(metadata.st_mode) or
                (metadata.st_dev, metadata.st_ino) !=
                (pathname_metadata.st_dev, pathname_metadata.st_ino)):
            raise TimingError("timing panel must be a plain directory")
        manifest_path = attempt_directory / PANEL_RESULT_NAME
        sidecar_path = attempt_directory / PANEL_RESULT_SIDECAR_NAME
        manifest_present = (
            manifest_path.exists() or manifest_path.is_symlink())
        sidecar_present = sidecar_path.exists() or sidecar_path.is_symlink()
        if manifest_present and not sidecar_present:
            # SIGKILL can leave the exact O_EXCL temporary for this missing
            # sidecar.  Remove only a dead writer's regular, unique partial
            # before exact inventory validation; live or malformed debris is
            # fatal and is never hidden.
            if _discard_stale_atomic_partials(sidecar_path):
                _fsync_directory(attempt_directory)
            result, encoded = _load_timing_panel_result(
                attempt_directory, spec, config, require_sidecar=False)
            # An existing sidecar, including a dangling symlink, is never
            # replaced.  Under the caller's timing-host lock this recheck is
            # expected to remain false; it also fails closed if that contract
            # is violated.
            if sidecar_path.exists() or sidecar_path.is_symlink():
                raise TimingError(
                    "timing panel result sidecar appeared during validation")
            digest = _sha256(encoded)
            _atomic_write(
                sidecar_path,
                (digest + "  " + PANEL_RESULT_NAME + "\n").encode("ascii"))
            loaded = load_timing_panel_result(attempt_directory, spec, config)
            if loaded != result:
                raise TimingError(
                    "repaired timing panel result did not round-trip")
            return loaded
        return load_timing_panel_result(attempt_directory, spec, config)
    return run_timing_panel(
        spec, config, evidence_probe, attempt_directory, executor=executor)


def _validate_panel_coverage(
    sample: Sequence[int], panels: Sequence[TimingPanelResult],
    setup_schedule: str,
) -> Tuple[Dict[Tuple[int, int, str], TimingPanelResult], Dict[Tuple[int, int], TimingPanelResult]]:
    solve: Dict[Tuple[int, int, str], TimingPanelResult] = {}
    setup: Dict[Tuple[int, int], TimingPanelResult] = {}
    sample_set = set(sample)
    runtime_binding: Optional[Tuple[Any, ...]] = None
    for panel in panels:
        if not isinstance(panel, TimingPanelResult):
            raise TimingError("timing campaign contains a non-panel result")
        if panel.spec.K not in sample_set:
            raise TimingError("timing panel K is outside the sealed sample")
        panel.panel_ratio()
        prefix, options = _attempt_command_options(panel.attempts[0].command)
        binding = (
            prefix, options["--evict-bytes"], options["--route-cache"],
            options["--route-cache-sha256"],
            options["--route-context-sha256"],
        )
        if runtime_binding is None:
            runtime_binding = binding
        elif binding != runtime_binding:
            raise TimingError("timing runtime binding changed across panels")
        if panel.spec.metric == "solve":
            key = (panel.spec.K, panel.spec.block_bytes, panel.spec.schedule)
            if key in solve:
                raise TimingError("duplicate solve timing panel")
            solve[key] = panel
        else:
            if panel.spec.schedule != setup_schedule:
                raise TimingError("setup panel used the wrong frozen schedule")
            key2 = (panel.spec.K, panel.spec.block_bytes)
            if key2 in setup:
                raise TimingError("duplicate setup timing panel")
            setup[key2] = panel
    expected_solve = {
        (K, width, schedule)
        for K in sample for width in WIDTHS for schedule in SCHEDULES
    }
    expected_setup = {(K, width) for K in sample for width in WIDTHS}
    if set(solve) != expected_solve or set(setup) != expected_setup:
        raise TimingError("timing panel coverage is not exact")
    for K in sample:
        K_attempts = set()
        for width in WIDTHS:
            attempts = {
                solve[(K, width, schedule)].spec.preferred_attempt
                for schedule in SCHEDULES
            }
            attempts.add(setup[(K, width)].spec.preferred_attempt)
            if len(attempts) != 1:
                raise TimingError("preferred route changed across panel strata")
            K_attempts.update(attempts)
            canonical_attempts = {
                solve[(K, width, schedule)].canonical_attempt
                for schedule in SCHEDULES
            }
            canonical_attempts.add(setup[(K, width)].canonical_attempt)
            if len(canonical_attempts) != 1:
                raise TimingError("canonical route changed across panel strata")
        if len(K_attempts) != 1:
            raise TimingError("preferred attempt changed across widths for one K")
    return solve, setup


def analyze_timing(
    sample: Sequence[int],
    panels: Sequence[TimingPanelResult],
    setup_schedule: str = "burst",
) -> Dict[str, Any]:
    """Apply all frozen solve/setup gates with exact ``Fraction`` arithmetic."""
    sample_tuple = tuple(sample)
    if (len(sample_tuple) < MIN_SAMPLE or len(sample_tuple) > MAX_SAMPLE or
            any(not isinstance(K, int) or isinstance(K, bool) or
                not 4096 <= K <= 64000 for K in sample_tuple) or
            tuple(sorted(sample_tuple)) != sample_tuple or
            len(set(sample_tuple)) != len(sample_tuple)):
        raise TimingError("timing sample is not canonical K-ascending S=5..32")
    if setup_schedule not in SCHEDULES:
        raise TimingError("setup schedule is invalid")
    panels_tuple = tuple(panels)
    solve, setup = _validate_panel_coverage(
        sample_tuple, panels_tuple, setup_schedule)

    solve_k: Dict[int, Fraction] = {}
    solve_width_k: Dict[int, Dict[int, Fraction]] = {width: {} for width in WIDTHS}
    for K in sample_tuple:
        nine = []
        for width in WIDTHS:
            schedules = [
                solve[(K, width, schedule)].panel_ratio()
                for schedule in SCHEDULES
            ]
            solve_width_k[width][K] = _upper_median(schedules)
            nine.extend(schedules)
        solve_k[K] = _upper_median(nine)
    solve_global = _upper_median(list(solve_k.values()))
    width_solve = {
        width: _upper_median(list(solve_width_k[width].values()))
        for width in WIDTHS
    }
    wins = sum(value < 1 for value in solve_k.values())
    losses = sum(value > 1 for value in solve_k.values())
    ties = len(sample_tuple) - wins - losses
    n = wins + losses
    sign_tail = sum(math.comb(n, j) for j in range(wins, n + 1)) if n else 0
    sign_pass = n > 0 and 20 * sign_tail <= (1 << n)
    solve_width_pass = all(
        200 * value.numerator <= 201 * value.denominator
        for value in width_solve.values())
    solve_pass = solve_global < 1 and sign_pass and solve_width_pass

    setup_k: Dict[int, Fraction] = {}
    setup_width_k: Dict[int, Dict[int, Fraction]] = {width: {} for width in WIDTHS}
    for K in sample_tuple:
        widths = []
        for width in WIDTHS:
            value = setup[(K, width)].panel_ratio()
            setup_width_k[width][K] = value
            widths.append(value)
        setup_k[K] = _upper_median(widths)
    setup_global = _upper_median(list(setup_k.values()))
    width_setup = {
        width: _upper_median(list(setup_width_k[width].values()))
        for width in WIDTHS
    }
    setup_global_pass = (
        200 * setup_global.numerator <= 201 * setup_global.denominator)
    setup_width_pass = all(
        200 * value.numerator <= 201 * value.denominator
        for value in width_setup.values())
    setup_pass = setup_global_pass and setup_width_pass

    S = len(sample_tuple)
    attempted_commands = sum(len(panel.attempts) for panel in panels_tuple)
    full_commands = sum(
        attempt.cycle_index is None
        for panel in panels_tuple for attempt in panel.attempts)
    replacement_commands = attempted_commands - full_commands
    physical_invocations = 32 * full_commands + 8 * replacement_commands
    return {
        "schema": "wirehair.wh2.h12_preferred_attempt.timing_analysis.v1",
        "sample": list(sample_tuple), "S": S,
        "accepted": solve_pass and setup_pass,
        "realized": {
            "solve_panels": 9 * S, "setup_panels": 3 * S,
            "solve_total_invocations": 288 * S,
            "solve_timed_invocations": 216 * S,
            "setup_total_invocations": 96 * S,
            "setup_timed_invocations": 72 * S,
            "accepted_total_invocations": 384 * S,
            "accepted_timed_invocations": 288 * S,
            "physical_process_attempts_including_retries": attempted_commands,
            "physical_full_process_attempts": full_commands,
            "physical_replacement_process_attempts": replacement_commands,
            "physical_invocations_including_retries": physical_invocations,
        },
        "solve": {
            "accepted": solve_pass,
            "global_ratio": _fraction_record(solve_global),
            "global_strictly_below_one": solve_global < 1,
            "per_width_ratio": {
                str(width): _fraction_record(value)
                for width, value in width_solve.items()
            },
            "per_width_at_most_1_005": solve_width_pass,
            "sign": {
                "wins": wins, "losses": losses, "ties": ties, "n": n,
                "tail_numerator": sign_tail,
                "denominator": 1 << n if n else 1,
                "twenty_tail_le_denominator": sign_pass,
            },
            "per_K_ratio": {
                str(K): _fraction_record(solve_k[K]) for K in sample_tuple
            },
        },
        "setup": {
            "accepted": setup_pass,
            "global_ratio": _fraction_record(setup_global),
            "global_at_most_1_005": setup_global_pass,
            "per_width_ratio": {
                str(width): _fraction_record(value)
                for width, value in width_setup.items()
            },
            "per_width_at_most_1_005": setup_width_pass,
            "per_K_ratio": {
                str(K): _fraction_record(setup_k[K]) for K in sample_tuple
            },
        },
    }
