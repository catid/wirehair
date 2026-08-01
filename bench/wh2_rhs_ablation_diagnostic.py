#!/usr/bin/env python3
"""Sealed initial diagnostic for the WH2 unused-row RHS ablation.

This is a narrow specialization of the authenticated grouped-commit timing
controller from commit ``5cb8d46``.  It deliberately does not contain a
second process runner.  Instead it verifies and imports that exact controller,
then replaces only the experiment contract, ``rhstiming`` parser, command,
task grid, and measurement overlay.  The inherited engine continues to own
bounded process groups, ABBABAAB outer ordering, CPU/NUMA placement, IRQ and
sibling-idle gates, sole-reader thermal evidence, immutable receipts, replay,
and publication.

The initial grid contains exactly 30 cells::

    (K, bb) = (3200,4096), (20000,4096), (32000,4096),
              (3200,1280), (20000,1280)
    schedule = burst, adversarial, repair-only
    cache = cold, warm

Each cell executes the untouched ``be0bc94`` baseline and ``d6439d6``
ablation in outer ABBABAAB order.  Each process uses the authenticated
``rhstiming`` v1 inner ABBABAAB fixture and discards inner cycle zero.  The
controller requires exact cross-binary and cold/warm graph, packet trace,
output, recovery, and non-time work identity before retaining timing.

``rhstiming`` v1 intentionally rejects ``R == 0`` and does not emit the
compile-time ``WH_COUNT`` counters.  The reduction therefore reports observed
``R == 0`` and ``R > 0`` populations and clears ``timing_promotional`` when
both are not represented.  It also records that low-level GF256 calls/bytes
are unavailable from the timing binary.  A future ``WH_COUNT`` diagnostic
must be a separate non-timing build; this controller never enables it in a
timed binary.  Native nonzero exits and stderr are substantive campaign-fatal
failures, including pre-emission ``rhstiming`` rejections; only a successful,
parseable execution with post-parse environmental contamination is eligible
for the inherited bounded retry loop.  Every emitted cycle must receipt
exactly zero minor and major faults to be accepted.

Typical use (only after competing formal workloads have stopped)::

    python3 bench/wh2_rhs_ablation_diagnostic.py selftest
    python3 bench/wh2_rhs_ablation_diagnostic.py \
      --engine /path/to/5cb8d46/bench/wh2_grouped_commit_timing.py \
      prepare --repo /path/to/wirehair --result-dir /dev/shm/result \
      --thermal-sampler /path/to/wirehair_expo_thermal_sampler.py \
      --core 6 --controller-core 80 --numa-node 0

Run the printed frozen controller path for ``run``, ``reduce``, and
``verify``.  Supplying the prepare-receipt SHA256 is mandatory.  No command
changes a production profile or the caller's worktree.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path
import py_compile
import re
import shutil
import subprocess
import sys
import tempfile
import types
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

BASE_COMMIT = "be0bc94b97d03d6ddbc23db3b7058aa7f575b5cd"
CANDIDATE_COMMIT = "d6439d6a10fee690bf9a39c172f50a06b0a1ab88"
RHS_OVERLAY_COMMIT = "eaf07929cb0e4c6a5895840bdeb89d466aace739"
RHS_OVERLAY_PARENT_COMMIT = "243f8ed86b7bf102fa1cb7156a481c170935e57b"
RHS_OVERLAY_TREE = "d5a1b4ccebfa44e8d87c2d6ae25c872bb72f0a39"
FORMAL_ENGINE_COMMIT = "5cb8d46eb0d9ffd8225586ce5ab64255020d939b"
FORMAL_ENGINE_SHA256 = \
    "772175233ae356ea89942be4d2a980a815a14f3b5975cd8bb8a566fbd0284dde"
FORMAL_HELPER_SHA256 = \
    "387929f39ec696ada3be32e47af9e97e1fa4444163af3cbf93af223c38edd47c"
RHS_TIMING_ONLY_BENCH_SHA256 = \
    "a9de1ab0df0d16750dc0b20f16f8caddb79fed5573d65f6a1b8c82520b825d3f"
DESIGN_SCHEMA = "wirehair.wh2.rhs-ablation-diagnostic.design.v1"
PREPARE_SCHEMA = "wirehair.wh2.rhs-ablation-diagnostic.prepare-receipt.v1"
SUMMARY_SCHEMA = "wirehair.wh2.rhs-ablation-diagnostic.validated-summary.v1"

RHS_OVERLAY_FILES = (
    "codec/WirehairV2Bench.cpp",
    "codec/V2BenchCliTest.cmake",
)
RHS_OVERLAY_BLOBS = {
    "codec/WirehairV2Bench.cpp": (
        "e36290e64b538ed0e7753ffb8ae7aab178891129",
        "1cb16c862f47ae2ada262cbe8987dcce0b9e13b6a178d24a6860a2c83756ed40",
    ),
    "codec/V2BenchCliTest.cmake": (
        "b428e146ec5f2aceee8121b9a5f8959a767e53bd",
        "061d2f5624b5292423b15505bb81c59b930c4a5eabf1cc5dfdf8ed7185ee34f0",
    ),
}

CELL_PAIRS = (
    (3200, 4096),
    (20000, 4096),
    (32000, 4096),
    (3200, 1280),
    (20000, 1280),
)
SCHEDULE_SEEDS = (
    ("burst", 15111065706836454659),
    ("adversarial", 10723151780598845931),
    ("repair-only", 9599682871048892067),
)
CACHE_STATES = ("cold", "warm")
OUTER_ORDER = "ABBABAAB"
INNER_ORDER = "ABBABAAB"
OVERHEAD = 4
EXPECT_Q = 8
LOSS_TEXT = "0.5"
OVERHEAD_STREAM = "salted"
TASK_ORDER_DOMAIN = b"wirehair.wh2.rhs-ablation.initial-grid.order.v1\0"

TIMING_C_FLAGS_RELEASE = "-O3 -DNDEBUG -ffunction-sections -fdata-sections"
TIMING_CXX_FLAGS_RELEASE = (
    TIMING_C_FLAGS_RELEASE +
    " -DWIREHAIR_V2_RHS_TIMING_ONLY=1"
    " -DWIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION=1"
    " -Wno-unused-function"
)
TIMING_EXE_LINKER_FLAGS = "-Wl,--gc-sections"
FORBIDDEN_TIMED_BINARY_SYMBOLS = (
    "CheckMixedRhsFusionOracleForTesting",
    "CheckPackedBinaryResidualOracleForTesting",
    "CmdSelfTest()", "CmdCompare(", "CmdPrecodeCheck(",
    "CmdSeedTable(", "CmdPeelCost(", "CmdDenseCheck(",
    "CmdDenseTune(", "CmdDenseCount(", "CmdDenseGrid(",
    "CmdPrecodeFail(",
    "CmdGroupedTiming(", "CmdPreferredTiming(", "CmdPreferredAttempt(",
    "gf256_count_bytes", "gf256_count_calls", "gf256_count_reset",
    " U fork@", " U vfork@", " U clone@", " U clone3@",
    " U posix_spawn@", " U popen@", " U pthread_create@",
    " U system@", " U syscall@",
)

RHS_FIELDS = (
    "N", "bb", "overhead", "expect_q", "overhead_stream", "schedule",
    "seed", "loss", "cache_state", "cycle", "slot", "arm", "field",
    "completion", "H", "D2", "S", "source_hits", "mix",
    "dense_two_anchor", "heavy_family", "seed_attempt",
    "packet_seed_attempt", "matrix_seed", "peel_seed", "graph_sha256",
    "trace_sha256", "packet_prefix_sha256", "payload_contract",
    "payload_sha256", "preflight_result", "output_sha256",
    "output_validation", "all_zero", "cell_class", "common_success",
    "result", "outcome_stable", "elapsed_ns", "saturated", "cpu_before",
    "cpu_after", "cpu_migrated", "minflt_delta", "majflt_delta",
    "fault_contaminated", "packet_rows", "peeled_columns", "inactivated",
    "residual_rows", "residual_rank", "binary_residual_rank", "q",
    "heavy_gain", "unused_binary_rows", "rhs_init_destination_bytes",
    "binary_row_references", "binary_row_storage_bytes",
    "binary_adjacency_storage_bytes", "binary_row_storage_allocations",
    "binary_adjacency_storage_allocations", "block_xors", "block_muladds",
    "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns",
    "phase_sum_ns", "mixed_joint_source_xors", "mixed_joint_marginal_xors",
    "mixed_joint_marginal_copies", "mixed_joint_active_deltas",
    "mixed_joint_scratch_bytes", "mixed_dual_source_columns", "source_bytes",
    "packet_payload_bytes", "intermediate_bytes",
)

RHS_PREAMBLE_FIELDS = (
    "schema", "policy", "timing_scope", "native_pair", "cycles", "order",
    "discard_cycle", "cycle_mode", "cycle_index", "N", "bb", "overhead",
    "expect_q", "overhead_stream", "loss", "seed", "schedule",
    "cache_state", "evict_bytes", "eviction_prefaulted", "completion",
    "field", "H", "D2", "S", "source_hits", "mix", "dense_two_anchor",
    "heavy_family", "base_matrix_seed", "base_peel_seed", "selected_attempt",
    "preflight_packet_seed_attempt", "matrix_seed", "peel_seed",
    "graph_sha256", "trace_sha256", "packet_prefix_sha256", "payload",
    "payload_count", "payload_bytes", "payload_alignment",
    "payload_prefaulted", "payload_sha256", "system_build",
    "output_validation", "allocator_state", "preflight_control_result",
    "preflight_candidate_result", "cell_class", "common_success",
    "preflight_output_all_zero", "preflight_output_sha256",
    "preflight_intermediate_bytes", "preflight_packet_rows",
    "preflight_peeled_columns", "preflight_inactivated",
    "preflight_residual_rows", "preflight_residual_rank",
    "preflight_binary_residual_rank", "preflight_q", "preflight_heavy_gain",
    "preflight_unused_binary_rows", "preflight_rhs_init_destination_bytes",
)

PHASE_FIELDS = (
    "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns",
)
WORK_FIELDS = (
    "preflight_result", "output_sha256", "all_zero", "result",
    "outcome_stable", "packet_rows", "peeled_columns", "inactivated",
    "residual_rows", "residual_rank", "binary_residual_rank", "q",
    "heavy_gain", "unused_binary_rows", "rhs_init_destination_bytes",
    "binary_row_references", "binary_row_storage_bytes",
    "binary_adjacency_storage_bytes", "binary_row_storage_allocations",
    "binary_adjacency_storage_allocations", "block_xors", "block_muladds",
    "mixed_joint_source_xors", "mixed_joint_marginal_xors",
    "mixed_joint_marginal_copies", "mixed_joint_active_deltas",
    "mixed_joint_scratch_bytes", "mixed_dual_source_columns", "source_bytes",
    "packet_payload_bytes", "intermediate_bytes",
)

SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
HEX_RE = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")
UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")
MAX_U64 = (1 << 64) - 1
MAX_I64 = (1 << 63) - 1
MAX_CPU_ID = (1 << 31) - 1


class DiagnosticError(RuntimeError):
    """The specialization or RHS evidence violated its frozen contract."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ensure_ascii=True) + "\n").encode("ascii")


def parse_uint(text: object, name: str, maximum: int = MAX_U64) -> int:
    if not isinstance(text, str) or UINT_RE.fullmatch(text) is None:
        raise DiagnosticError("%s is not a canonical unsigned integer" % name)
    value = int(text, 10)
    if value > maximum:
        raise DiagnosticError("%s exceeds its integer domain" % name)
    return value


def parse_sint(
    text: object, name: str, minimum: int = -MAX_I64 - 1,
    maximum: int = MAX_I64,
) -> int:
    if not isinstance(text, str) or SINT_RE.fullmatch(text) is None:
        raise DiagnosticError("%s is not a canonical signed integer" % name)
    value = int(text, 10)
    if not minimum <= value <= maximum:
        raise DiagnosticError("%s exceeds its integer domain" % name)
    return value


def require_sha256(value: object, name: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise DiagnosticError("%s is not a lowercase SHA256" % name)
    return value


def _read_plain_file(path: Path, maximum: int) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise DiagnosticError("not a plain file: %s" % path)
    metadata_before = path.stat()
    if metadata_before.st_size > maximum:
        raise DiagnosticError("file exceeds its bound: %s" % path)
    raw = path.read_bytes()
    metadata_after = path.stat()
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns)
    if len(raw) != metadata_before.st_size or \
            identity(metadata_before) != identity(metadata_after):
        raise DiagnosticError("file identity changed while reading: %s" % path)
    return raw


_FORMAL_ENGINE_MODULE_NAME = "_wirehair_rhs_formal_engine"
_FORMAL_HELPER_MODULE_NAME = "wh2_timing_evidence_io"
_MISSING_MODULE = object()


def _restore_module(name: str, previous: object) -> None:
    if previous is _MISSING_MODULE:
        sys.modules.pop(name, None)
    else:
        sys.modules[name] = previous


def _exec_authenticated_module(
    name: str, path: Path, source: bytes,
) -> types.ModuleType:
    """Compile and execute exactly ``source`` without an import loader."""
    module = types.ModuleType(name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    module.__cached__ = None
    sys.modules[name] = module
    try:
        code = compile(source, str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
    except BaseException:
        if sys.modules.get(name) is module:
            del sys.modules[name]
        raise
    return module


def load_formal_engine(path: Path):
    """Execute only the exact authenticated 5cb8d46 source bytes."""
    engine_path = path.resolve()
    helper_path = engine_path.with_name("wh2_timing_evidence_io.py")
    engine_raw = _read_plain_file(engine_path, 1024 * 1024)
    helper_raw = _read_plain_file(helper_path, 1024 * 1024)
    if sha256_bytes(engine_raw) != FORMAL_ENGINE_SHA256:
        raise DiagnosticError("formal timing engine SHA256 mismatch")
    # This value is assigned once below; the separate local makes it difficult
    # for a mistaken engine-only update to bypass the helper binding.
    helper_sha256 = sha256_bytes(helper_raw)
    if helper_sha256 != _formal_helper_sha256():
        raise DiagnosticError("formal timing helper SHA256 mismatch")

    # Both modules are compiled directly from the authenticated byte strings.
    # In particular, the mutable source directory is never put on sys.path,
    # and neither a sibling shadow module nor a timestamp-valid .pyc can be
    # selected by an import loader.  The exact helper is present before the
    # engine executes its ordinary ``import wh2_timing_evidence_io`` line.
    previous_helper = sys.modules.get(
        _FORMAL_HELPER_MODULE_NAME, _MISSING_MODULE)
    previous_engine = sys.modules.get(
        _FORMAL_ENGINE_MODULE_NAME, _MISSING_MODULE)
    try:
        helper = _exec_authenticated_module(
            _FORMAL_HELPER_MODULE_NAME, helper_path, helper_raw)
        module = _exec_authenticated_module(
            _FORMAL_ENGINE_MODULE_NAME, engine_path, engine_raw)
        if (getattr(module, "evidence_io", None) is not helper or
                sys.modules.get(_FORMAL_ENGINE_MODULE_NAME) is not module or
                sys.modules.get(_FORMAL_HELPER_MODULE_NAME) is not helper or
                module.__file__ != str(engine_path) or
                helper.__file__ != str(helper_path) or
                module.__spec__ is not None or helper.__spec__ is not None or
                getattr(module, "BASE_COMMIT", None) !=
                    "48d14bc77e3f9e98605fca4d226aa218d7d03a0d" or
                getattr(module, "CANDIDATE_COMMIT", None) != BASE_COMMIT or
                getattr(module, "OUTER_ORDER", None) != OUTER_ORDER):
            raise DiagnosticError(
                "formal timing engine semantic anchors changed")
        if (_read_plain_file(engine_path, 1024 * 1024) != engine_raw or
                _read_plain_file(helper_path, 1024 * 1024) != helper_raw):
            raise DiagnosticError(
                "formal timing engine/helper changed during loading")
    except BaseException:
        _restore_module(_FORMAL_ENGINE_MODULE_NAME, previous_engine)
        _restore_module(_FORMAL_HELPER_MODULE_NAME, previous_helper)
        raise
    return module


def _formal_helper_sha256() -> str:
    return FORMAL_HELPER_SHA256


def generate_tasks() -> Tuple[Dict[str, object], ...]:
    pending: List[Tuple[str, Dict[str, object]]] = []
    for K, width in CELL_PAIRS:
        for seed_index, (schedule, seed) in enumerate(SCHEDULE_SEEDS):
            for cache_state in CACHE_STATES:
                task: Dict[str, object] = {
                    "K": K, "bb": width, "seed_index": seed_index,
                    "seed": seed, "schedule": schedule,
                    "cache_state": cache_state,
                }
                priority = sha256_bytes(TASK_ORDER_DOMAIN + canonical_json(task))
                pending.append((priority, task))
    pending.sort(key=lambda item: item[0])
    result: List[Dict[str, object]] = []
    for job, (_priority, task) in enumerate(pending):
        value = dict(task)
        value["job"] = job
        value["task_id"] = (
            "%03d.K%d.bb%d.seed%d.%s.%s" %
            (job, task["K"], task["bb"], task["seed_index"],
             task["schedule"], task["cache_state"]))
        result.append(value)
    if len(result) != 30 or len({row["task_id"] for row in result}) != 30:
        raise DiagnosticError("internal RHS task grid changed")
    return tuple(result)


@dataclass(frozen=True)
class ParsedOutput:
    label: str
    schema: str
    preamble: Mapping[str, str]
    rows: Tuple[Mapping[str, str], ...]
    work_signature: Tuple[str, ...]
    semantic_sha256: str
    stdout_sha256: str
    timed_elapsed_ns: int
    all_elapsed_ns: int
    timed_phase_ns: Mapping[str, int]
    all_phase_ns: Mapping[str, int]
    timed_minor_faults: int
    discard_minor_faults: int
    max_minor_faults: int
    contaminations: Tuple[str, ...]


def _parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# rhstiming: "
    if not line.startswith(prefix):
        raise DiagnosticError("missing rhstiming preamble")
    tokens = line[len(prefix):].split(" ")
    if any(token.count("=") != 1 for token in tokens):
        raise DiagnosticError("rhstiming preamble token is malformed")
    pairs = tuple(tuple(token.split("=", 1)) for token in tokens)
    if tuple(pair[0] for pair in pairs) != RHS_PREAMBLE_FIELDS:
        raise DiagnosticError("rhstiming preamble order/schema mismatch")
    if len({pair[0] for pair in pairs}) != len(pairs):
        raise DiagnosticError("rhstiming preamble contains a duplicate")
    return dict(pairs)


def _expected_preamble(
    task: Mapping[str, object], evict_bytes: int,
) -> Dict[str, str]:
    K = int(task["K"])
    width = int(task["bb"])
    return {
        "schema": "v1", "policy": "certified-gf256-rhs",
        "timing_scope": "solve",
        "native_pair": "identical-negative-control", "cycles": "4",
        "order": INNER_ORDER, "discard_cycle": "0", "cycle_mode": "full",
        "cycle_index": "all", "N": str(K), "bb": str(width),
        "overhead": str(OVERHEAD), "expect_q": str(EXPECT_Q),
        "overhead_stream": OVERHEAD_STREAM, "loss": LOSS_TEXT,
        "seed": str(task["seed"]), "schedule": str(task["schedule"]),
        "cache_state": str(task["cache_state"]),
        "evict_bytes": str(evict_bytes), "eviction_prefaulted": "1",
        "completion": "gf256", "field": "gf256", "H": "12", "D2": "12",
        "source_hits": "2", "mix": "3", "dense_two_anchor": "0",
        "heavy_family": "periodic", "payload": "distinct-packet-zero-v1",
        "payload_count": str(K + OVERHEAD),
        "payload_bytes": str((K + OVERHEAD) * width),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "system_build": "outside-timer",
        "output_validation": "exact-all-zero-full-scan-v1",
        "allocator_state": "preflight-warmed",
        "preflight_control_result": "0",
        "preflight_candidate_result": "0", "cell_class": "common-success",
        "common_success": "1", "preflight_output_all_zero": "1",
        "preflight_packet_rows": str(K + OVERHEAD),
        "preflight_q": str(EXPECT_Q),
        "preflight_heavy_gain": str(EXPECT_Q),
    }


def _hex_value(text: object, name: str, maximum: int) -> int:
    if not isinstance(text, str) or HEX_RE.fullmatch(text) is None:
        raise DiagnosticError("%s is not canonical hexadecimal" % name)
    value = int(text, 16)
    if value > maximum:
        raise DiagnosticError("%s exceeds its integer domain" % name)
    return value


def parse_rhs_output(
    raw: bytes, label: str, task: Mapping[str, object], evict_bytes: int,
    expected_core: int,
) -> ParsedOutput:
    if label not in ("base", "candidate"):
        raise DiagnosticError("unknown binary label")
    if (not raw or len(raw) > 2 * 1024 * 1024 or not raw.endswith(b"\n") or
            b"\r" in raw or b"\0" in raw or b'"' in raw):
        raise DiagnosticError("rhstiming output is not canonical LF text")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise DiagnosticError("rhstiming output is not ASCII") from exc
    if len(lines) != 34:
        raise DiagnosticError("rhstiming output does not have 34 lines")
    if any(line.count(",") != len(RHS_FIELDS) - 1 for line in lines[1:]):
        raise DiagnosticError("rhstiming CSV row arity changed")
    preamble = _parse_preamble(lines[0])
    for key, expected in _expected_preamble(task, evict_bytes).items():
        if preamble.get(key) != expected:
            raise DiagnosticError(
                "rhstiming preamble mismatch %s: %r != %r" %
                (key, preamble.get(key), expected))
    for field in (
        "graph_sha256", "trace_sha256", "packet_prefix_sha256",
        "payload_sha256", "preflight_output_sha256",
    ):
        require_sha256(preamble.get(field), field)
    S = parse_uint(preamble.get("S"), "S", (1 << 32) - 1)
    parse_uint(preamble.get("selected_attempt"), "selected_attempt", 255)
    parse_uint(
        preamble.get("preflight_packet_seed_attempt"),
        "preflight_packet_seed_attempt", 255)
    for field, maximum in (
        ("base_matrix_seed", MAX_U64), ("matrix_seed", MAX_U64),
        ("base_peel_seed", (1 << 32) - 1),
        ("peel_seed", (1 << 32) - 1),
    ):
        _hex_value(preamble.get(field), field, maximum)
    numeric_preamble = (
        "preflight_intermediate_bytes", "preflight_peeled_columns",
        "preflight_inactivated", "preflight_residual_rows",
        "preflight_residual_rank", "preflight_binary_residual_rank",
        "preflight_unused_binary_rows",
        "preflight_rhs_init_destination_bytes",
    )
    for field in numeric_preamble:
        parse_uint(preamble.get(field), field)
    preflight_R = parse_uint(
        preamble["preflight_inactivated"], "preflight R")
    preflight_peeled = parse_uint(
        preamble["preflight_peeled_columns"], "preflight peeled")
    preflight_intermediate = parse_uint(
        preamble["preflight_intermediate_bytes"],
        "preflight intermediate")
    preflight_residual_rows = parse_uint(
        preamble["preflight_residual_rows"], "preflight residual rows")
    preflight_unused = parse_uint(
        preamble["preflight_unused_binary_rows"],
        "preflight unused binary rows")
    if preflight_R == 0:
        raise DiagnosticError(
            "authenticated rhstiming v1 unexpectedly accepted R == 0")
    if (preflight_intermediate % int(task["bb"]) or
            preflight_intermediate // int(task["bb"]) !=
            int(task["K"]) + S + 12 + 12 or
            preflight_peeled + preflight_R !=
            preflight_intermediate // int(task["bb"]) or
            preflight_residual_rows < 12 or
            preflight_unused == 0 or
            preflight_unused != preflight_residual_rows - 12 or
            parse_uint(
                preamble["preflight_rhs_init_destination_bytes"],
                "preflight RHS destination bytes") !=
            preflight_unused * int(task["bb"])):
        raise DiagnosticError("preflight dimensional receipt changed")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != RHS_FIELDS:
        raise DiagnosticError("rhstiming CSV schema mismatch")
    rows = tuple(dict(row) for row in reader)
    if len(rows) != 32 or any(
            set(row) != set(RHS_FIELDS) or
            any(value is None for value in row.values()) for row in rows):
        raise DiagnosticError("rhstiming CSV contains malformed rows")
    K = int(task["K"])
    width = int(task["bb"])
    preflight_receipts = {
        "packet_rows": preamble["preflight_packet_rows"],
        "peeled_columns": preamble["preflight_peeled_columns"],
        "inactivated": preamble["preflight_inactivated"],
        "residual_rows": preamble["preflight_residual_rows"],
        "residual_rank": preamble["preflight_residual_rank"],
        "binary_residual_rank": preamble["preflight_binary_residual_rank"],
        "q": preamble["preflight_q"],
        "heavy_gain": preamble["preflight_heavy_gain"],
        "unused_binary_rows": preamble["preflight_unused_binary_rows"],
        "rhs_init_destination_bytes":
            preamble["preflight_rhs_init_destination_bytes"],
    }
    signatures = set()
    contaminations: List[str] = []
    timed_elapsed = 0
    all_elapsed = 0
    timed_phase = {field: 0 for field in PHASE_FIELDS}
    all_phase = {field: 0 for field in PHASE_FIELDS}
    timed_minor = 0
    discard_minor = 0
    max_minor = 0
    for index, row in enumerate(rows):
        cycle, slot = divmod(index, 8)
        arm = "control" if INNER_ORDER[slot] == "A" else "candidate"
        exact = {
            "N": str(K), "bb": str(width), "overhead": str(OVERHEAD),
            "expect_q": str(EXPECT_Q), "overhead_stream": OVERHEAD_STREAM,
            "schedule": str(task["schedule"]), "seed": str(task["seed"]),
            "loss": LOSS_TEXT, "cache_state": str(task["cache_state"]),
            "cycle": str(cycle), "slot": str(slot), "arm": arm,
            "field": "gf256", "completion": "gf256", "H": "12",
            "D2": "12", "S": preamble["S"], "source_hits": "2",
            "mix": "3", "dense_two_anchor": "0",
            "heavy_family": "periodic",
            "seed_attempt": preamble["selected_attempt"],
            "packet_seed_attempt":
                preamble["preflight_packet_seed_attempt"],
            "matrix_seed": preamble["matrix_seed"],
            "peel_seed": preamble["peel_seed"],
            "graph_sha256": preamble["graph_sha256"],
            "trace_sha256": preamble["trace_sha256"],
            "packet_prefix_sha256": preamble["packet_prefix_sha256"],
            "payload_contract": "distinct-packet-zero-v1",
            "payload_sha256": preamble["payload_sha256"],
            "preflight_result": "0",
            "output_sha256": preamble["preflight_output_sha256"],
            "output_validation": "exact-all-zero-full-scan-v1",
            "all_zero": "1", "cell_class": "common-success",
            "common_success": "1", "result": "0", "outcome_stable": "1",
            "packet_rows": str(K + OVERHEAD), "q": str(EXPECT_Q),
            "heavy_gain": str(EXPECT_Q),
            "source_bytes": str(K * width),
            "packet_payload_bytes": str((K + OVERHEAD) * width),
            **preflight_receipts,
        }
        for key, expected in exact.items():
            if row.get(key) != expected:
                raise DiagnosticError("row %d mismatch %s" % (index, key))
        unsigned = (
            "elapsed_ns", "saturated", "packet_rows", "peeled_columns",
            "inactivated", "residual_rows", "residual_rank",
            "binary_residual_rank", "q", "heavy_gain", "unused_binary_rows",
            "rhs_init_destination_bytes", "binary_row_references",
            "binary_row_storage_bytes", "binary_adjacency_storage_bytes",
            "binary_row_storage_allocations",
            "binary_adjacency_storage_allocations", "block_xors",
            "block_muladds", "build_ns", "peel_ns", "project_ns",
            "residual_ns", "backsub_ns", "phase_sum_ns",
            "mixed_joint_source_xors", "mixed_joint_marginal_xors",
            "mixed_joint_marginal_copies", "mixed_joint_active_deltas",
            "mixed_joint_scratch_bytes", "mixed_dual_source_columns",
            "source_bytes", "packet_payload_bytes", "intermediate_bytes",
        )
        for field in unsigned:
            parse_uint(row.get(field), "row %d %s" % (index, field))
        elapsed = parse_uint(row["elapsed_ns"], "elapsed_ns")
        if elapsed == 0:
            raise DiagnosticError("timing row has zero elapsed time")
        phases = [parse_uint(row[field], field) for field in PHASE_FIELDS]
        if sum(phases) != parse_uint(row["phase_sum_ns"], "phase_sum_ns") or \
                sum(phases) > elapsed:
            raise DiagnosticError("timing phase receipt is inconsistent")
        inactivated = parse_uint(row["inactivated"], "inactivated")
        binary_rank = parse_uint(
            row["binary_residual_rank"], "binary_residual_rank")
        residual_rank = parse_uint(row["residual_rank"], "residual_rank")
        if (inactivated == 0 or binary_rank > inactivated or
                residual_rank < binary_rank or
                inactivated - binary_rank != EXPECT_Q or
                residual_rank - binary_rank != EXPECT_Q):
            raise DiagnosticError("residual rank relationship changed")
        if (parse_uint(row["peeled_columns"], "peeled columns") +
                inactivated !=
                parse_uint(row["intermediate_bytes"], "intermediate") //
                width or
                parse_uint(row["unused_binary_rows"], "unused rows") !=
                parse_uint(row["residual_rows"], "residual rows") - 12):
            raise DiagnosticError("row dimensional receipt changed")
        if parse_uint(
                row["rhs_init_destination_bytes"],
                "rhs_init_destination_bytes") != \
                parse_uint(row["unused_binary_rows"], "unused rows") * width:
            raise DiagnosticError("RHS destination-byte receipt changed")
        intermediate = parse_uint(row["intermediate_bytes"], "intermediate")
        if (intermediate !=
                parse_uint(preamble["preflight_intermediate_bytes"],
                           "preflight intermediate") or
                intermediate < K * width or intermediate % width):
            raise DiagnosticError("intermediate storage receipt changed")
        for field in (
            "mixed_joint_source_xors", "mixed_joint_marginal_xors",
            "mixed_joint_marginal_copies", "mixed_joint_active_deltas",
            "mixed_joint_scratch_bytes", "mixed_dual_source_columns",
        ):
            if row[field] != "0":
                raise DiagnosticError("generic GF256 row emitted mixed work")
        cpu_before = parse_sint(row["cpu_before"], "cpu_before", -1, MAX_CPU_ID)
        cpu_after = parse_sint(row["cpu_after"], "cpu_after", -1, MAX_CPU_ID)
        migrated = parse_sint(row["cpu_migrated"], "cpu_migrated", -1, 1)
        minflt = parse_sint(row["minflt_delta"], "minflt_delta", -1, MAX_I64)
        majflt = parse_sint(row["majflt_delta"], "majflt_delta", -1, MAX_I64)
        fault = parse_sint(row["fault_contaminated"], "fault", -1, 1)
        saturated = parse_uint(row["saturated"], "saturated", 1)
        if cpu_before != expected_core or cpu_after != expected_core or migrated != 0:
            contaminations.append(
                "row%d:migration:%d->%d:%d" %
                (index, cpu_before, cpu_after, migrated))
        if saturated:
            contaminations.append("row%d:saturated" % index)
        expected_fault = -1 if minflt < 0 or majflt < 0 else \
            (1 if minflt or majflt else 0)
        if fault != expected_fault:
            raise DiagnosticError("fault receipt disagrees with counters")
        if minflt != 0:
            contaminations.append("row%d:minor-fault:%d" % (index, minflt))
        if majflt != 0:
            contaminations.append("row%d:major-fault:%d" % (index, majflt))
        all_elapsed += elapsed
        for field, value in zip(PHASE_FIELDS, phases):
            all_phase[field] += value
        max_minor = max(max_minor, minflt)
        if cycle == 0:
            discard_minor += max(minflt, 0)
        else:
            timed_elapsed += elapsed
            timed_minor += max(minflt, 0)
            for field, value in zip(PHASE_FIELDS, phases):
                timed_phase[field] += value
        signatures.add(tuple(row[field] for field in WORK_FIELDS))
    if len(signatures) != 1:
        raise DiagnosticError("native A/B arms changed deterministic work")
    signature = next(iter(signatures))
    semantic = {
        "preamble": {
            key: value for key, value in preamble.items() if key != "schema"},
        "work_fields": list(WORK_FIELDS),
        "work_signature": list(signature), "row_count": 32,
        "timed_cycles": [1, 2, 3],
    }
    return ParsedOutput(
        label=label, schema="v1", preamble=preamble, rows=rows,
        work_signature=signature,
        semantic_sha256=sha256_bytes(canonical_json(semantic)),
        stdout_sha256=sha256_bytes(raw), timed_elapsed_ns=timed_elapsed,
        all_elapsed_ns=all_elapsed, timed_phase_ns=timed_phase,
        all_phase_ns=all_phase, timed_minor_faults=timed_minor,
        discard_minor_faults=discard_minor, max_minor_faults=max_minor,
        contaminations=tuple(contaminations))


def execution_name(task: Mapping[str, object], slot: int, label: str) -> str:
    if label not in ("base", "candidate") or not 0 <= slot < len(OUTER_ORDER):
        raise DiagnosticError("execution name coordinates are invalid")
    return "%s.outer%d.%s.csv" % (task["task_id"], slot, label)


def command_for(
    design: Mapping[str, object], task: Mapping[str, object], label: str,
) -> List[str]:
    if label not in ("base", "candidate"):
        raise DiagnosticError("unknown binary label")
    root = Path(str(design["root"]))
    tools = design["tools"]
    binary_name = "wirehair_v2_bench.A" if label == "base" else \
        "wirehair_v2_bench.B"
    core = int(design["core"])
    node = int(design["numa_node"])
    return [
        str(tools["env"]["path"]), "-i", "LC_ALL=C", "TZ=UTC",
        "PATH=/usr/bin:/bin", "MALLOC_MMAP_THRESHOLD_=1073741824",
        "MALLOC_TRIM_THRESHOLD_=-1", str(tools["taskset"]["path"]),
        "-c", str(core), str(tools["numactl"]["path"]),
        "--physcpubind=" + str(core), "--membind=" + str(node),
        str(root / "frozen" / binary_name), "rhstiming",
        "--N", str(task["K"]), "--bb", str(task["bb"]),
        "--overhead", str(OVERHEAD), "--expect-q", str(EXPECT_Q),
        "--overhead-stream", OVERHEAD_STREAM,
        "--evict-bytes", str(design["evict_bytes"]),
        "--cache-state", str(task["cache_state"]), "--loss", LOSS_TEXT,
        "--seed", str(task["seed"]), "--schedule", str(task["schedule"]),
    ]


def _register_cross_cache_identity(
    ledger: Dict[Tuple[object, ...], Dict[str, Dict[str, object]]],
    task: Mapping[str, object], parsed: ParsedOutput,
) -> None:
    key = (
        task["K"], task["bb"], task["seed_index"], task["seed"],
        task["schedule"])
    cache = str(task["cache_state"])
    if cache not in CACHE_STATES:
        raise DiagnosticError("cross-cache ledger has unknown state")
    identity: Dict[str, object] = {
        "graph_sha256": parsed.preamble["graph_sha256"],
        "trace_sha256": parsed.preamble["trace_sha256"],
        "packet_prefix_sha256": parsed.preamble["packet_prefix_sha256"],
        "payload_sha256": parsed.preamble["payload_sha256"],
        "selected_attempt": parsed.preamble["selected_attempt"],
        "matrix_seed": parsed.preamble["matrix_seed"],
        "peel_seed": parsed.preamble["peel_seed"],
        "output_sha256": parsed.preamble["preflight_output_sha256"],
        "work_signature": list(parsed.work_signature),
    }
    states = ledger.setdefault(key, {})
    if cache in states:
        raise DiagnosticError("cross-cache ledger contains duplicate cell")
    states[cache] = identity
    if len(states) == len(CACHE_STATES) and len({
            sha256_bytes(canonical_json(value)) for value in states.values()
    }) != 1:
        raise DiagnosticError("cold/warm graph, trace, output, or work changed")


def _validate_cross_cache_ledger(
    ledger: Mapping[Tuple[object, ...], Mapping[str, Mapping[str, object]]],
) -> None:
    expected = len(CELL_PAIRS) * len(SCHEDULE_SEEDS)
    if len(ledger) != expected or any(
            set(states) != set(CACHE_STATES) for states in ledger.values()):
        raise DiagnosticError("cross-cache ledger is incomplete")


_USAGE_NEEDLE = b"""    if (argc < 2) {\n#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \\\n    !defined(WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT)\n"""
_USAGE_REPLACEMENT = b"""    if (argc < 2) {\n#if defined(WIREHAIR_V2_RHS_TIMING_ONLY)\n        std::fprintf(stderr,\n            \"usage: wirehair_v2_bench rhstiming [opts]\\n\");\n#elif defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \\\n    !defined(WIREHAIR_V2_BENCH_DISABLE_PREFERRED_ATTEMPT)\n"""
_DISPATCH_NEEDLE = b"""    try\n    {\n        if (!std::strcmp(argv[1], \"compare\")) {\n"""
_DISPATCH_REPLACEMENT = b"""    try\n    {\n#if defined(WIREHAIR_V2_RHS_TIMING_ONLY)\n        if (!std::strcmp(argv[1], \"rhstiming\")) {\n            return CmdRhsTiming(argc - 2, argv + 2);\n        }\n#else\n        if (!std::strcmp(argv[1], \"compare\")) {\n"""
_CATCH_NEEDLE = b"""    }\n    catch (const std::bad_alloc&)\n    {\n        std::fprintf(stderr,\n            \"wirehair_v2_bench: allocation failed for requested workload\\n\");\n"""
_CATCH_REPLACEMENT = b"""#endif\n    }\n    catch (const std::bad_alloc&)\n    {\n        std::fprintf(stderr,\n            \"wirehair_v2_bench: allocation failed for requested workload\\n\");\n"""


def rhs_timing_only_overlay(raw: bytes) -> bytes:
    """Apply the minimal d6ab-style timing-only main dispatch to eaf0792."""
    result = raw
    for needle, replacement, name in (
        (_USAGE_NEEDLE, _USAGE_REPLACEMENT, "usage"),
        (_DISPATCH_NEEDLE, _DISPATCH_REPLACEMENT, "dispatch"),
        (_CATCH_NEEDLE, _CATCH_REPLACEMENT, "catch"),
    ):
        if result.count(needle) != 1:
            raise DiagnosticError(
                "authenticated rhstiming source has %d %s anchors" %
                (result.count(needle), name))
        result = result.replace(needle, replacement, 1)
    if (result.count(b"WIREHAIR_V2_RHS_TIMING_ONLY") != 2 or
            result.count(b"return CmdRhsTiming(argc - 2, argv + 2);") != 2):
        raise DiagnosticError("RHS timing-only transformation is malformed")
    return result


def _formal_helper_digest_selftest() -> None:
    # Keep this assertion separate from engine loading so ``selftest`` catches
    # source-level digest typos without needing a repository or engine path.
    digests = (
        FORMAL_ENGINE_SHA256, FORMAL_HELPER_SHA256,
        RHS_TIMING_ONLY_BENCH_SHA256,
        *(digest for _blob, digest in RHS_OVERLAY_BLOBS.values()),
    )
    if (FORMAL_HELPER_SHA256 != _formal_helper_sha256() or
            any(SHA256_RE.fullmatch(digest) is None for digest in digests)):
        raise DiagnosticError("authenticated digest literal is malformed")


def _fault_policy() -> Dict[str, object]:
    return {
        "minimum": 0, "maximum": 0, "major_faults": 0,
        "all_cycles_receipted": True,
        "unavailable_or_nonzero_is_contamination": True,
    }


def _execution_failure_policy(engine) -> Dict[str, object]:
    return {
        "post_parse_contamination_max_attempts":
            engine.MAX_ENVIRONMENTAL_ATTEMPTS,
        "native_nonzero_exit_or_stderr": "substantive-campaign-fatal",
        "native_pre_emission_rejection":
            "substantive-campaign-fatal-no-retry",
    }


def _configure_engine(engine) -> None:
    """Install the narrow RHS contract into the authenticated engine."""
    if getattr(engine, "_rhs_ablation_configured", False):
        return
    engine._rhs_original_verify_immutable = engine._verify_immutable
    engine._rhs_original_replay_campaign = engine.replay_campaign
    engine.BASE_COMMIT = BASE_COMMIT
    engine.CANDIDATE_COMMIT = CANDIDATE_COMMIT
    engine.MEASUREMENT_OVERLAY_COMMIT = RHS_OVERLAY_COMMIT
    engine.MEASUREMENT_OVERLAY_PARENT_COMMIT = RHS_OVERLAY_PARENT_COMMIT
    engine.MEASUREMENT_OVERLAY_FILES = RHS_OVERLAY_FILES
    engine.KS = tuple(sorted({K for K, _width in CELL_PAIRS}))
    engine.WIDTHS = tuple(sorted({width for _K, width in CELL_PAIRS}))
    engine.SCHEDULE_SEEDS = SCHEDULE_SEEDS
    engine.CACHE_STATES = CACHE_STATES
    engine.OUTER_ORDER = OUTER_ORDER
    engine.INNER_ORDER = INNER_ORDER
    engine.OVERHEAD = OVERHEAD
    engine.LOSS_TEXT = LOSS_TEXT
    engine.ARCHITECTURE = {
        "policy": "certified-gf256-rhs", "H": 12, "D2": 12,
        "mix": 3, "expect_q": EXPECT_Q,
        "overhead_stream": OVERHEAD_STREAM,
    }
    engine.TIMING_C_FLAGS_RELEASE = TIMING_C_FLAGS_RELEASE
    engine.TIMING_CXX_FLAGS_RELEASE = TIMING_CXX_FLAGS_RELEASE
    engine.TIMING_EXE_LINKER_FLAGS = TIMING_EXE_LINKER_FLAGS
    engine.FORBIDDEN_TIMED_BINARY_SYMBOLS = FORBIDDEN_TIMED_BINARY_SYMBOLS
    engine.GROUPED_FIELDS = RHS_FIELDS
    engine.GROUPED_PREAMBLE_FIELDS = RHS_PREAMBLE_FIELDS
    engine.WORK_FIELDS = WORK_FIELDS
    engine.PHASE_FIELDS = PHASE_FIELDS
    engine.PRIMARY_PHASE_FIELD = "residual_ns"
    engine.NEGATIVE_CONTROL_PHASE_FIELDS = (
        "build_ns", "peel_ns", "project_ns", "backsub_ns")
    # The specialized parser classifies every unavailable or nonzero per-row
    # fault counter as contamination, so the inherited acceptance ceiling must
    # say zero as well.  Keeping the engine's generic value of 64 would make the
    # sealed design less strict than the execution/replay path.
    engine.MAX_MINOR_FAULTS = 0
    engine.MAX_GROUPED_OUTPUT_BYTES = 2 * 1024 * 1024
    engine.generate_tasks = generate_tasks

    def parse_for_engine(*args, **kwargs):
        try:
            return parse_rhs_output(*args, **kwargs)
        except DiagnosticError as exc:
            raise engine.TimingError(str(exc)) from exc

    engine.parse_grouped_output = parse_for_engine
    engine.execution_name = execution_name
    engine.command_for = command_for
    engine._register_cross_cache_identity = _register_cross_cache_identity
    engine._validate_cross_cache_ledger = _validate_cross_cache_ledger
    engine._overlay_identity = lambda record: _overlay_identity(engine, record)
    engine._verify_overlay_status = lambda git, source, files=None: \
        _verify_overlay_status(engine, git, source, files)
    engine._apply_measurement_overlay = \
        lambda label, repo, source, provenance, git: \
        _apply_measurement_overlay(engine, label, repo, source, provenance, git)
    engine._prepare_cross_binary_smoke = \
        lambda staging, provenance, tools, core, node: \
        _prepare_smoke(engine, staging, provenance, tools, core, node)
    engine._load_design = lambda root: _load_design(engine, root)
    engine._validate_prepare_receipt = \
        lambda root, design: _validate_prepare_receipt(engine, root, design)
    engine._verify_immutable = \
        lambda root, design: _verify_immutable(engine, root, design)
    engine.replay_campaign = lambda root: _replay_campaign(engine, root)
    engine._rhs_ablation_configured = True


def _overlay_identity(engine, record: object) -> Dict[str, object]:
    identity_fields = (
        "source_commit", "source_parent_commit", "source_tree", "files",
        "diff_options", "diff_sha256", "stable_patch_id",
        "timing_only_transform",
    )
    if not isinstance(record, dict) or set(record) != \
            set(identity_fields) | {"diff_evidence_name"}:
        raise engine.TimingError("RHS overlay provenance fields changed")
    if (record.get("source_commit") != RHS_OVERLAY_COMMIT or
            record.get("source_parent_commit") != RHS_OVERLAY_PARENT_COMMIT or
            record.get("source_tree") != RHS_OVERLAY_TREE or
            record.get("diff_options") != list(engine.MEASUREMENT_DIFF_OPTIONS)):
        raise engine.TimingError("RHS overlay identity changed")
    engine.require_sha256(record.get("diff_sha256"), "RHS overlay diff")
    patch_id = record.get("stable_patch_id")
    if not isinstance(patch_id, str) or engine.GIT_OID_RE.fullmatch(patch_id) is None:
        raise engine.TimingError("RHS overlay patch-id is malformed")
    transform = record.get("timing_only_transform")
    if (not isinstance(transform, dict) or set(transform) != {
            "schema", "macro", "authenticated_bench_sha256",
            "timing_only_bench_sha256"} or
            transform.get("schema") !=
                "wirehair.wh2.rhs-timing-only-transform.v1" or
            transform.get("macro") != "WIREHAIR_V2_RHS_TIMING_ONLY"):
        raise engine.TimingError("RHS timing-only transformation changed")
    for field in ("authenticated_bench_sha256", "timing_only_bench_sha256"):
        engine.require_sha256(transform.get(field), field)
    files = record.get("files")
    expected_fields = {
        "path", "parent_blob_oid", "parent_sha256",
        "authenticated_blob_oid", "authenticated_sha256", "overlay_sha256"}
    if (not isinstance(files, list) or len(files) != len(RHS_OVERLAY_FILES) or
            [item.get("path") if isinstance(item, dict) else None
             for item in files] != list(RHS_OVERLAY_FILES)):
        raise engine.TimingError("RHS overlay file ledger changed")
    for item in files:
        if not isinstance(item, dict) or set(item) != expected_fields:
            raise engine.TimingError("RHS overlay file receipt changed")
        path = str(item["path"])
        expected_blob, expected_sha = RHS_OVERLAY_BLOBS[path]
        if (item["authenticated_blob_oid"] != expected_blob or
                item["authenticated_sha256"] != expected_sha):
            raise engine.TimingError("authenticated RHS overlay blob changed")
        if (not isinstance(item.get("parent_blob_oid"), str) or
                engine.GIT_OID_RE.fullmatch(item["parent_blob_oid"]) is None):
            raise engine.TimingError("RHS parent blob ID is malformed")
        engine.require_sha256(item.get("parent_sha256"), "parent file")
        engine.require_sha256(item.get("overlay_sha256"), "overlay file")
    if transform["authenticated_bench_sha256"] != \
            RHS_OVERLAY_BLOBS["codec/WirehairV2Bench.cpp"][1]:
        raise engine.TimingError("RHS timing transform input changed")
    overlays = {str(item["path"]): item["overlay_sha256"] for item in files}
    if (transform["timing_only_bench_sha256"] !=
            RHS_TIMING_ONLY_BENCH_SHA256 or
            overlays.get("codec/WirehairV2Bench.cpp") !=
            RHS_TIMING_ONLY_BENCH_SHA256 or
            overlays.get("codec/V2BenchCliTest.cmake") !=
            RHS_OVERLAY_BLOBS["codec/V2BenchCliTest.cmake"][1]):
        raise engine.TimingError("RHS transformed overlay content changed")
    return {field: record[field] for field in identity_fields}


def _verify_overlay_status(engine, git: Path, source: Path, files=None) -> None:
    status = engine._git_status_lines(git, source)
    expected = {" M " + path for path in RHS_OVERLAY_FILES}
    if set(status) != expected or len(status) != len(expected):
        raise engine.TimingError("build worktree differs from exact RHS overlay")
    if files is not None:
        if [item.get("path") for item in files] != list(RHS_OVERLAY_FILES):
            raise engine.TimingError("RHS overlay paths changed")
        for item in files:
            path = source / str(item["path"])
            if (path.is_symlink() or not path.is_file() or
                    engine.sha256_file(path) != item.get("overlay_sha256")):
                raise engine.TimingError("RHS overlay content changed")


def _apply_measurement_overlay(
    engine, label: str, repo: Path, source: Path, provenance: Path, git: Path,
) -> Dict[str, object]:
    if engine._git_status_lines(git, source):
        raise engine.TimingError("detached codec worktree is not clean")
    files: List[Dict[str, object]] = []
    authenticated_bench_sha = ""
    timing_bench_sha = ""
    for relative in RHS_OVERLAY_FILES:
        destination = source / relative
        if destination.is_symlink() or not destination.is_file():
            raise engine.TimingError("RHS overlay target is not regular")
        parent_revision = "HEAD:" + relative
        overlay_revision = RHS_OVERLAY_COMMIT + ":" + relative
        parent = engine._git_blob(git, source, parent_revision)
        authenticated = engine._git_blob(git, repo, overlay_revision)
        expected_blob, expected_sha = RHS_OVERLAY_BLOBS[relative]
        if (engine._git_value(git, repo, "rev-parse", overlay_revision) !=
                expected_blob or sha256_bytes(authenticated) != expected_sha):
            raise engine.TimingError("authenticated RHS overlay blob mismatch")
        if engine.stable_bytes(destination, attempts=3) != parent:
            raise engine.TimingError("codec-parent file disagrees with Git")
        final = rhs_timing_only_overlay(authenticated) if \
            relative == "codec/WirehairV2Bench.cpp" else authenticated
        destination.write_bytes(final)
        if relative == "codec/WirehairV2Bench.cpp":
            authenticated_bench_sha = expected_sha
            timing_bench_sha = sha256_bytes(final)
        files.append({
            "path": relative,
            "parent_blob_oid": engine._git_value(
                git, source, "rev-parse", parent_revision),
            "parent_sha256": sha256_bytes(parent),
            "authenticated_blob_oid": expected_blob,
            "authenticated_sha256": expected_sha,
            "overlay_sha256": sha256_bytes(final),
        })
    _verify_overlay_status(engine, git, source, files)
    diff_result = subprocess.run(
        (str(git), "-C", str(source), "diff", *engine.MEASUREMENT_DIFF_OPTIONS,
         "--", *RHS_OVERLAY_FILES), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if diff_result.returncode != 0 or diff_result.stderr or not diff_result.stdout:
        raise engine.TimingError("could not capture RHS measurement overlay")
    name = label + ".measurement-overlay.diff"
    engine.write_new(provenance / name, diff_result.stdout)
    record = {
        "source_commit": RHS_OVERLAY_COMMIT,
        "source_parent_commit": RHS_OVERLAY_PARENT_COMMIT,
        "source_tree": RHS_OVERLAY_TREE, "files": files,
        "diff_options": list(engine.MEASUREMENT_DIFF_OPTIONS),
        "diff_sha256": sha256_bytes(diff_result.stdout),
        "stable_patch_id": engine._stable_patch_id(git, diff_result.stdout),
        "timing_only_transform": {
            "schema": "wirehair.wh2.rhs-timing-only-transform.v1",
            "macro": "WIREHAIR_V2_RHS_TIMING_ONLY",
            "authenticated_bench_sha256": authenticated_bench_sha,
            "timing_only_bench_sha256": timing_bench_sha,
        },
        "diff_evidence_name": name,
    }
    _overlay_identity(engine, record)
    return record


def _prepare_smoke(
    engine, staging: Path, provenance: Path, tools: Mapping[str, Path],
    core: int, node: int,
) -> Dict[str, object]:
    """Run one bounded compatibility check; never treat it as timing."""
    task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 1280 and
        value["schedule"] == "burst" and value["cache_state"] == "warm")
    design = {
        "root": str(staging),
        "tools": {
            name: {"path": str(path), "sha256": engine.sha256_file(path)}
            for name, path in tools.items()
        },
        "core": core, "numa_node": node, "evict_bytes": 4096,
    }
    observations: List[ParsedOutput] = []
    records: Dict[str, Dict[str, object]] = {}
    for label in ("base", "candidate"):
        argv = command_for(design, task, label)
        result = subprocess.run(
            argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=300, check=False, start_new_session=True)
        if result.returncode != 0 or result.stderr:
            raise engine.TimingError(
                "%s RHS prepare smoke failed exit=%d stderr=%r" %
                (label, result.returncode, result.stderr[:1000]))
        try:
            observation = parse_rhs_output(
                result.stdout, label, task, 4096, core)
        except DiagnosticError as exc:
            raise engine.TimingError(str(exc)) from exc
        name = "smoke.%s.csv" % label
        engine.write_new(provenance / name, result.stdout)
        records[label] = {
            "argv": argv, "stdout_name": name,
            "stdout_sha256": observation.stdout_sha256,
            "schema_version": observation.schema,
            "semantic_sha256": observation.semantic_sha256,
            "graph_sha256": observation.preamble["graph_sha256"],
            "trace_sha256": observation.preamble["trace_sha256"],
            "packet_prefix_sha256":
                observation.preamble["packet_prefix_sha256"],
            "output_sha256":
                observation.preamble["preflight_output_sha256"],
            "work_signature_sha256": sha256_bytes(canonical_json({
                "fields": list(WORK_FIELDS),
                "values": list(observation.work_signature),
            })),
            "inactivated": int(observation.preamble["preflight_inactivated"]),
            "nonpromotional_contaminations":
                list(observation.contaminations),
        }
        observations.append(observation)
    identity = {
        (item.semantic_sha256, item.preamble["graph_sha256"],
         item.preamble["trace_sha256"],
         item.preamble["packet_prefix_sha256"],
         item.preamble["preflight_output_sha256"], item.work_signature)
        for item in observations
    }
    if len(identity) != 1:
        raise engine.TimingError(
            "RHS prepare smoke found cross-binary identity drift")
    return {
        "scope": "bounded compatibility smoke under ordinary machine load",
        "timing_evidence": False, "task": task, "binaries": records,
        "graph_trace_packet_output_work_identity_exact": True,
        "recovery_identity_exact": True,
    }


_SOLVE_SYMBOL = "wirehair_v2::SolvePrecodeSystemImpl("
_GF256_SYMBOL_PREFIXES = (
    "gf256_add_mem", "gf256_add2_mem", "gf256_addset_multi_mem",
    "gf256_add_multi_mem", "gf256_mul_mem", "gf256_muladd_mem",
)


def _layout_command(
    engine, argv: Sequence[str], description: str,
) -> bytes:
    result = subprocess.run(
        list(argv), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "TZ": "UTC"})
    if (result.returncode != 0 or result.stderr or not result.stdout or
            not result.stdout.endswith(b"\n") or b"\0" in result.stdout or
            len(result.stdout) > 16 * 1024 * 1024):
        raise engine.TimingError(
            "%s failed exit=%d stderr=%r" %
            (description, result.returncode, result.stderr[:1000]))
    return result.stdout


def _parse_layout(
    engine, nm_raw: bytes, elf_raw: bytes, binary_sha256: str,
) -> Dict[str, object]:
    try:
        nm_text = nm_raw.decode("utf-8", "strict")
        elf_text = elf_raw.decode("utf-8", "strict")
    except UnicodeDecodeError as exc:
        raise engine.TimingError("binary-layout evidence is not UTF-8") from exc
    header: Dict[str, str] = {}
    for key in ("Class", "Data", "Type", "Machine"):
        matches = re.findall(r"^\s*%s:\s*(.*?)\s*$" % key, elf_text, re.M)
        if len(matches) != 1 or not matches[0]:
            raise engine.TimingError("ELF header field is missing: %s" % key)
        header[key.lower()] = matches[0]
    text_section = None
    for line in elf_text.splitlines():
        fields = line.split()
        try:
            index = fields.index(".text")
        except ValueError:
            continue
        if len(fields) < index + 10 or text_section is not None:
            raise engine.TimingError("ELF .text section row is malformed")
        address, offset, size = fields[index + 2:index + 5]
        alignment = fields[index + 9]
        if (re.fullmatch(r"[0-9a-fA-F]+", address) is None or
                re.fullmatch(r"[0-9a-fA-F]+", offset) is None or
                re.fullmatch(r"[0-9a-fA-F]+", size) is None or
                UINT_RE.fullmatch(alignment) is None):
            raise engine.TimingError("ELF .text coordinates are malformed")
        text_section = {
            "address": int(address, 16), "offset": int(offset, 16),
            "size": int(size, 16), "alignment": int(alignment, 10),
        }
    if text_section is None or text_section["size"] <= 0:
        raise engine.TimingError("ELF .text section is unavailable")
    symbols: List[Dict[str, object]] = []
    for line in nm_text.splitlines():
        fields = line.split(None, 3)
        if (len(fields) != 4 or
                re.fullmatch(r"[0-9a-fA-F]+", fields[0]) is None or
                re.fullmatch(r"[0-9a-fA-F]+", fields[1]) is None or
                len(fields[2]) != 1):
            continue
        symbols.append({
            "address": int(fields[0], 16), "size": int(fields[1], 16),
            "type": fields[2], "name": fields[3],
        })
    solve = [
        row for row in symbols
        if str(row["name"]).startswith(_SOLVE_SYMBOL) and
        str(row["name"]).endswith(", bool)")
    ]
    if len(solve) != 1 or int(solve[0]["size"]) <= 0:
        raise engine.TimingError(
            "SolvePrecodeSystemImpl layout symbol is not unique")
    solve_address = int(solve[0]["address"])
    prefix = [row for row in symbols if int(row["address"]) < solve_address]
    gf256 = [
        row for row in symbols
        if any(str(row["name"]).startswith(prefix_name)
               for prefix_name in _GF256_SYMBOL_PREFIXES)
    ]
    if not gf256:
        raise engine.TimingError("GF256 layout symbols are unavailable")
    return {
        "schema": "wirehair.wh2.rhs-ablation.binary-layout.v1",
        "binary_sha256": binary_sha256,
        "nm_sha256": sha256_bytes(nm_raw),
        "readelf_sha256": sha256_bytes(elf_raw),
        "elf_header": header, ".text": text_section,
        "solve_symbol": solve[0],
        "solve_section_offset": solve_address - int(text_section["address"]),
        "solve_alignment_mod_64": solve_address % 64,
        "prefix_symbol_ledger_sha256": sha256_bytes(canonical_json(prefix)),
        "prefix_symbol_count": len(prefix),
        "gf256_symbols": gf256,
    }


def _capture_layout(
    engine, binary: Path, label: str, provenance: Path,
    tools: Mapping[str, Path],
) -> Dict[str, object]:
    nm_raw = _layout_command(
        engine, (str(tools["nm"]), "-S", "-n", "-C", str(binary)),
        label + " layout nm")
    elf_raw = _layout_command(
        engine, (str(tools["readelf"]), "-W", "-h", "-S", str(binary)),
        label + " layout readelf")
    nm_name = label + ".layout.nm.txt"
    elf_name = label + ".layout.readelf.txt"
    engine.write_new(provenance / nm_name, nm_raw)
    engine.write_new(provenance / elf_name, elf_raw)
    parsed = _parse_layout(
        engine, nm_raw, elf_raw, engine.sha256_file(binary))
    parsed["nm_evidence_name"] = nm_name
    parsed["readelf_evidence_name"] = elf_name
    return parsed


def _validate_layout_record(
    engine, value: object, label: str,
) -> Mapping[str, object]:
    required = {
        "schema", "binary_sha256", "nm_sha256", "readelf_sha256",
        "elf_header", ".text", "solve_symbol", "solve_section_offset",
        "solve_alignment_mod_64", "prefix_symbol_ledger_sha256",
        "prefix_symbol_count", "gf256_symbols",
    }
    if (not isinstance(value, dict) or
            set(value) not in (required, required | {
                "nm_evidence_name", "readelf_evidence_name"}) or
            value.get("schema") !=
                "wirehair.wh2.rhs-ablation.binary-layout.v1"):
        raise engine.TimingError("%s binary-layout fields changed" % label)
    for field in (
            "binary_sha256", "nm_sha256", "readelf_sha256",
            "prefix_symbol_ledger_sha256"):
        try:
            require_sha256(value.get(field), "%s %s" % (label, field))
        except DiagnosticError as exc:
            raise engine.TimingError(str(exc)) from exc
    header = value.get("elf_header")
    if (not isinstance(header, dict) or
            set(header) != {"class", "data", "type", "machine"} or
            any(not isinstance(item, str) or not item
                for item in header.values())):
        raise engine.TimingError("%s ELF header receipt changed" % label)
    text = value.get(".text")
    if not isinstance(text, dict) or set(text) != {
            "address", "offset", "size", "alignment"}:
        raise engine.TimingError("%s .text receipt changed" % label)
    for field in ("address", "offset", "size", "alignment"):
        item = text.get(field)
        if (not isinstance(item, int) or isinstance(item, bool) or item < 0 or
                item > MAX_U64):
            raise engine.TimingError(
                "%s .text coordinate changed: %s" % (label, field))
    if text["size"] == 0 or text["alignment"] == 0:
        raise engine.TimingError("%s .text has an empty domain" % label)
    solve = value.get("solve_symbol")
    if (not isinstance(solve, dict) or
            set(solve) != {"address", "size", "type", "name"} or
            not isinstance(solve.get("address"), int) or
            isinstance(solve.get("address"), bool) or
            not isinstance(solve.get("size"), int) or
            isinstance(solve.get("size"), bool) or solve.get("size", 0) <= 0 or
            not isinstance(solve.get("type"), str) or
            len(solve.get("type", "")) != 1 or
            not isinstance(solve.get("name"), str) or
            not solve["name"].startswith(_SOLVE_SYMBOL) or
            not solve["name"].endswith(", bool)")):
        raise engine.TimingError("%s solve-symbol receipt changed" % label)
    solve_address = int(solve["address"])
    section_offset = value.get("solve_section_offset")
    alignment_mod = value.get("solve_alignment_mod_64")
    if (not isinstance(section_offset, int) or isinstance(section_offset, bool) or
            section_offset != solve_address - int(text["address"]) or
            section_offset < 0 or section_offset >= int(text["size"]) or
            not isinstance(alignment_mod, int) or isinstance(alignment_mod, bool) or
            alignment_mod != solve_address % 64):
        raise engine.TimingError("%s solve placement receipt changed" % label)
    prefix_count = value.get("prefix_symbol_count")
    if (not isinstance(prefix_count, int) or isinstance(prefix_count, bool) or
            prefix_count < 0):
        raise engine.TimingError("%s prefix-symbol count changed" % label)
    gf256 = value.get("gf256_symbols")
    if not isinstance(gf256, list) or not gf256:
        raise engine.TimingError("%s GF256 layout ledger changed" % label)
    previous = -1
    for row in gf256:
        if (not isinstance(row, dict) or
                set(row) != {"address", "size", "type", "name"} or
                not isinstance(row.get("address"), int) or
                isinstance(row.get("address"), bool) or
                not isinstance(row.get("size"), int) or
                isinstance(row.get("size"), bool) or row.get("size", 0) <= 0 or
                not isinstance(row.get("type"), str) or
                len(row.get("type", "")) != 1 or
                not isinstance(row.get("name"), str) or
                not any(row["name"].startswith(prefix)
                        for prefix in _GF256_SYMBOL_PREFIXES) or
                row["address"] <= previous):
            raise engine.TimingError("%s GF256 symbol receipt changed" % label)
        previous = int(row["address"])
    if not set(_GF256_SYMBOL_PREFIXES) <= {
            str(row["name"]) for row in gf256}:
        raise engine.TimingError("%s GF256 entry symbols are incomplete" % label)
    for field in ("nm_evidence_name", "readelf_evidence_name"):
        if field in value and (
                not isinstance(value[field], str) or not value[field] or
                "/" in value[field] or "\\" in value[field]):
            raise engine.TimingError("%s layout evidence name changed" % label)
    return value


def _layout_gate(engine, layouts: Mapping[str, Mapping[str, object]]) \
        -> Dict[str, object]:
    if set(layouts) != {"base", "candidate"}:
        raise engine.TimingError("binary-layout pair is incomplete")
    base = _validate_layout_record(engine, layouts["base"], "base")
    candidate = _validate_layout_record(
        engine, layouts["candidate"], "candidate")
    base_text = base[".text"]
    candidate_text = candidate[".text"]
    for key in ("elf_header", "solve_section_offset",
                "solve_alignment_mod_64", "prefix_symbol_ledger_sha256",
                "prefix_symbol_count"):
        if base.get(key) != candidate.get(key):
            raise engine.TimingError("same-layout gate failed: %s" % key)
    for key in ("address", "offset", "alignment"):
        if base_text.get(key) != candidate_text.get(key):
            raise engine.TimingError("same-layout .text gate failed: %s" % key)
    base_solve = base["solve_symbol"]
    candidate_solve = candidate["solve_symbol"]
    for key in ("address", "type", "name"):
        if base_solve.get(key) != candidate_solve.get(key):
            raise engine.TimingError("same-layout solve gate failed: %s" % key)
    def gf_contract(layout: Mapping[str, object]) -> List[Dict[str, object]]:
        return [{
            "name": row["name"], "type": row["type"], "size": row["size"],
            "alignment_mod_64": int(row["address"]) % 64,
        } for row in layout["gf256_symbols"]]
    if gf_contract(base) != gf_contract(candidate):
        raise engine.TimingError("same-layout GF256 symbol contract changed")
    return {
        "schema": "wirehair.wh2.rhs-ablation.same-layout-gate.v1",
        "passed": True,
        "contract": (
            "same ELF ABI, .text origin/file offset/alignment, exact symbol "
            "prefix, solve entry/type/alignment, and GF256 name/type/size/"
            "alignment; changed solve/.text sizes are measured, not hidden"),
        "solve_address": base_solve["address"],
        "base_solve_size": base_solve["size"],
        "candidate_solve_size": candidate_solve["size"],
        "solve_size_delta":
            int(candidate_solve["size"]) - int(base_solve["size"]),
        "base_text_size": base_text["size"],
        "candidate_text_size": candidate_text["size"],
        "text_size_delta":
            int(candidate_text["size"]) - int(base_text["size"]),
        "gf256_contract_sha256":
            sha256_bytes(canonical_json(gf_contract(base))),
    }


def _cpu_context(engine, core: int) -> Dict[str, object]:
    try:
        raw = Path("/proc/cpuinfo").read_bytes()
    except OSError as exc:
        raise engine.TimingError("cannot read CPU context") from exc
    if not raw or len(raw) > 16 * 1024 * 1024 or b"\0" in raw:
        raise engine.TimingError("CPU context input is malformed")
    try:
        text = raw.decode("ascii", "strict")
    except UnicodeDecodeError as exc:
        raise engine.TimingError("CPU context is not ASCII") from exc
    selected: Optional[Dict[str, str]] = None
    for stanza in text.strip().split("\n\n"):
        fields: Dict[str, str] = {}
        for line in stanza.splitlines():
            if "\t:" not in line:
                continue
            key, value = line.split("\t:", 1)
            key = key.strip()
            value = value.strip()
            if not key or key in fields:
                raise engine.TimingError("CPU context contains duplicate fields")
            fields[key] = value
        if fields.get("processor") == str(core):
            if selected is not None:
                raise engine.TimingError("CPU context core is duplicated")
            selected = fields
    required = (
        "vendor_id", "cpu family", "model", "model name", "stepping",
        "microcode", "flags",
    )
    if selected is None or any(not selected.get(field) for field in required):
        raise engine.TimingError("CPU context core is unavailable")
    flags = selected["flags"].split()
    if (not flags or len(flags) != len(set(flags)) or
            any(re.fullmatch(r"[a-z0-9_]+", flag) is None for flag in flags) or
            "avx2" not in flags):
        raise engine.TimingError("CPU context lacks a canonical AVX2 domain")
    return {
        "schema": "wirehair.wh2.rhs-ablation.cpu-context.v1",
        "processor": core,
        "vendor_id": selected["vendor_id"],
        "cpu_family": selected["cpu family"],
        "model": selected["model"],
        "model_name": selected["model name"],
        "stepping": selected["stepping"],
        "microcode": selected["microcode"],
        "flags": flags,
        "dispatch_capabilities": {
            feature: feature in flags
            for feature in ("ssse3", "avx2", "gfni", "avx512f")
        },
        "measurement_claim": (
            "host capability and -march=native context only; rhstiming v1 "
            "does not emit the selected low-level call path"),
    }


def prepare_campaign(engine, args: argparse.Namespace) -> None:
    """Prepare a sealed campaign without running the 30-cell timing grid."""
    result = Path(args.result_dir).resolve()
    repo = Path(args.repo).resolve()
    if result.exists():
        raise engine.TimingError("result directory already exists")
    if not repo.is_dir():
        raise engine.TimingError("repository directory does not exist")
    sampler_source = Path(args.thermal_sampler).resolve()
    sampler_raw = engine.stable_bytes(
        sampler_source, attempts=3, max_bytes=1024 * 1024,
        require_unique=False)
    if not sampler_raw.startswith(b"#!/usr/bin/env python3\n"):
        raise engine.TimingError("thermal sampler has an unexpected header")
    if (args.core < 0 or args.controller_core < 0 or args.numa_node < 0 or
            args.evict_bytes < 4096 or args.build_jobs <= 0):
        raise engine.TimingError("prepare integer is outside its domain")
    if (engine.CLOCK_TICKS_PER_SECOND <= 0 or
            engine.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS !=
            (1_000_000_000 + engine.CLOCK_TICKS_PER_SECOND - 1) //
            engine.CLOCK_TICKS_PER_SECOND):
        raise engine.TimingError("SC_CLK_TCK is outside its domain")
    engine.require_pidfd_runtime()
    tool_names = (
        "git", "cmake", "env", "taskset", "numactl", "ldd", "nm",
        "readelf", "sudo", "fuser", "true",
    )
    tools = {name: engine.resolve_tool(name) for name in tool_names}
    python = Path(sys.executable).resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise engine.TimingError("active Python is not executable")
    tools["python"] = python
    sudo_probe = subprocess.run(
        (str(tools["sudo"]), "-n", str(tools["true"])),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if sudo_probe.returncode != 0 or sudo_probe.stdout or sudo_probe.stderr:
        raise engine.TimingError(
            "passwordless sudo is unavailable for thermal verification")
    ninja = shutil.which("ninja")
    if ninja is not None:
        tools["ninja"] = Path(ninja).resolve()
    c_compiler = Path(args.c_compiler).resolve() if args.c_compiler else \
        engine.resolve_tool("cc")
    cxx_compiler = Path(args.cxx_compiler).resolve() if args.cxx_compiler else \
        engine.resolve_tool("c++")
    for compiler, name in ((c_compiler, "C compiler"),
                           (cxx_compiler, "C++ compiler")):
        if not compiler.is_file() or not os.access(compiler, os.X_OK):
            raise engine.TimingError("%s is not executable" % name)
    for commit, name in (
            (BASE_COMMIT, "base"), (CANDIDATE_COMMIT, "candidate"),
            (RHS_OVERLAY_PARENT_COMMIT, "RHS overlay parent"),
            (RHS_OVERLAY_COMMIT, "RHS overlay")):
        if engine._git_value(
                tools["git"], repo, "rev-parse", commit + "^{commit}") != commit:
            raise engine.TimingError("%s commit is unavailable" % name)
    if engine._git_value(
            tools["git"], repo, "rev-parse", CANDIDATE_COMMIT + "^") != \
            BASE_COMMIT:
        raise engine.TimingError("candidate is not the direct child of base")
    if (engine._git_value(
            tools["git"], repo, "rev-parse", RHS_OVERLAY_COMMIT + "^") !=
            RHS_OVERLAY_PARENT_COMMIT or
            engine._git_value(
                tools["git"], repo, "rev-parse",
                RHS_OVERLAY_COMMIT + "^{tree}") != RHS_OVERLAY_TREE):
        raise engine.TimingError("authenticated RHS overlay ancestry changed")
    topology = engine.topology_record(args.core, args.numa_node)
    cpu_context = _cpu_context(engine, args.core)
    controller_topology = engine.topology_record(
        args.controller_core, args.numa_node)
    if args.controller_core in topology["llc_shared_cpus"]:
        raise engine.TimingError("controller CPU shares the timing LLC")
    if args.evict_bytes < max(
            2 * int(topology["llc_bytes"]), 256 * 1024 * 1024):
        raise engine.TimingError(
            "eviction allocation is smaller than the frozen LLC gate")
    staging = result.with_name(result.name + ".prepare.%d" % os.getpid())
    workspace = result.with_name(result.name + ".build.%d" % os.getpid())
    if staging.exists() or workspace.exists():
        raise engine.TimingError("stale prepare workspace exists")
    staging.mkdir(mode=0o700)
    workspace.mkdir(mode=0o700)
    os.chmod(staging, 0o700)
    os.chmod(workspace, 0o700)
    staging_sealed = False
    try:
        frozen = staging / "frozen"
        provenance = staging / "provenance"
        frozen.mkdir(mode=0o700)
        provenance.mkdir(mode=0o700)
        for directory in engine.PREPARED_CAMPAIGN_DIRECTORIES[2:]:
            (staging / directory).mkdir(mode=0o700)
        wrapper_source = Path(__file__).resolve()
        formal_source = Path(str(engine.__file__)).resolve()
        helper_source = Path(str(engine.evidence_io.__file__)).resolve()
        if (formal_source.parent != helper_source.parent or
                helper_source.name != "wh2_timing_evidence_io.py" or
                engine.sha256_file(formal_source) != FORMAL_ENGINE_SHA256 or
                engine.sha256_file(helper_source) != FORMAL_HELPER_SHA256):
            raise engine.TimingError("formal engine/helper binding changed")
        wrapper_frozen = frozen / "wh2_rhs_ablation_diagnostic.py"
        formal_frozen = frozen / "wh2_grouped_commit_timing.py"
        helper_frozen = frozen / "wh2_timing_evidence_io.py"
        sampler_frozen = frozen / engine.THERMAL_SAMPLER_NAME
        engine.write_new(
            wrapper_frozen, engine.stable_bytes(wrapper_source, attempts=3),
            mode=0o555)
        engine.write_new(
            formal_frozen, engine.stable_bytes(formal_source, attempts=3),
            mode=0o555)
        engine.write_new(
            helper_frozen, engine.stable_bytes(helper_source, attempts=3))
        engine.write_new(sampler_frozen, sampler_raw, mode=0o555)
        builds: Dict[str, Dict[str, object]] = {}
        for label, commit in (("base", BASE_COMMIT),
                              ("candidate", CANDIDATE_COMMIT)):
            builds[label] = engine._build_one(
                label, commit, repo, workspace, frozen, provenance, tools,
                args.build_jobs, c_compiler, cxx_compiler)
        base_overlay = _overlay_identity(
            engine, builds["base"]["record"].get("measurement_overlay"))
        candidate_overlay = _overlay_identity(
            engine,
            builds["candidate"]["record"].get("measurement_overlay"))
        if base_overlay != candidate_overlay:
            raise engine.TimingError(
                "binaries did not receive an identical RHS overlay")
        layouts = {
            label: _capture_layout(
                engine, frozen / engine.BINARY_NAMES[label], label,
                provenance, tools)
            for label in ("base", "candidate")
        }
        layout_gate = _layout_gate(engine, layouts)
        smoke = _prepare_smoke(
            engine, staging, provenance, tools, args.core, args.numa_node)
        if (engine.sha256_file(wrapper_source) !=
                engine.sha256_file(wrapper_frozen) or
                engine.sha256_file(formal_source) !=
                engine.sha256_file(formal_frozen) or
                engine.sha256_file(helper_source) !=
                engine.sha256_file(helper_frozen)):
            raise engine.TimingError(
                "controller source changed during preparation")
        tasks = generate_tasks()
        manifest = b"".join(canonical_json(task) for task in tasks)
        engine.write_new(staging / "tasks_manifest.jsonl", manifest)
        immutable: Dict[str, str] = {}
        for directory in (frozen, provenance):
            for path in sorted(directory.iterdir()):
                if path.is_file():
                    immutable[str(path.relative_to(staging))] = \
                        engine.sha256_file(path)
        tool_records = {
            name: {"path": str(path), "sha256": engine.sha256_file(path)}
            for name, path in sorted(tools.items())
        }
        design_payload: Dict[str, object] = {
            "root": str(result), "base_commit": BASE_COMMIT,
            "candidate_commit": CANDIDATE_COMMIT,
            "codec_lineage": {
                "candidate_parent_commit": BASE_COMMIT,
                "rhs_overlay_commit": RHS_OVERLAY_COMMIT,
                "rhs_overlay_parent_commit": RHS_OVERLAY_PARENT_COMMIT,
                "rhs_overlay_tree": RHS_OVERLAY_TREE,
                "overlay_is_measurement_only_and_branch_independent": True,
            },
            "measurement_overlay": base_overlay,
            "formal_engine": {
                "commit": FORMAL_ENGINE_COMMIT,
                "engine_sha256": FORMAL_ENGINE_SHA256,
                "helper_sha256": FORMAL_HELPER_SHA256,
                "delegated_controls": [
                    "bounded-process-groups", "outer-ABBABAAB",
                    "CPU-NUMA-affinity", "IRQ-and-sibling-idle-gates",
                    "sole-reader-CPU-DIMM-thermal-evidence",
                    "tmpfs-identity",
                    "bounded-post-parse-contamination-retries",
                    "exact-raw-receipt-replay",
                ],
            },
            "controller": {
                "name": wrapper_frozen.name,
                "sha256": immutable["frozen/" + wrapper_frozen.name],
            },
            "timing_scope": "full-payload decoder precode solve",
            "architecture": engine.ARCHITECTURE,
            "cell_pairs": [list(item) for item in CELL_PAIRS],
            "K": list(engine.KS), "bb": list(engine.WIDTHS),
            "schedule_seeds": [list(item) for item in SCHEDULE_SEEDS],
            "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
            "expect_q": EXPECT_Q, "overhead_stream": OVERHEAD_STREAM,
            "loss": LOSS_TEXT, "outer_order": OUTER_ORDER,
            "inner_order": INNER_ORDER, "inner_cycles": 4,
            "discard_inner_cycle": 0, "task_count": len(tasks),
            "processes_per_task": len(OUTER_ORDER),
            "rows_per_process": 32, "timed_rows_per_process": 24,
            "timed_rows_per_binary_per_task": 96,
            "allocator_environment": {
                "MALLOC_MMAP_THRESHOLD_": engine.MALLOC_MMAP_THRESHOLD,
                "MALLOC_TRIM_THRESHOLD_": engine.MALLOC_TRIM_THRESHOLD,
                "purpose": (
                    "retain and prefault solve-output pages in each process "
                    "while preserving all timed decoder work"),
            },
            "phase_timing": {
                "primary": engine.PRIMARY_PHASE_FIELD,
                "negative_controls":
                    list(engine.NEGATIVE_CONTROL_PHASE_FIELDS),
                "fields": list(engine.PHASE_FIELDS),
            },
            "timed_binary_build": {
                "rhs_timing_only": True,
                "mode": "rhstiming",
                "wh_count": False,
                "c_release_flags": engine.TIMING_C_FLAGS_RELEASE,
                "cxx_release_flags": engine.TIMING_CXX_FLAGS_RELEASE,
                "exe_linker_flags": engine.TIMING_EXE_LINKER_FLAGS,
                "forbidden_symbols_absent":
                    list(engine.FORBIDDEN_TIMED_BINARY_SYMBOLS),
                "unrelated_modes_rejected": True,
            },
            "r_population_policy": {
                "report_r_zero_and_positive_separately": True,
                "both_required_for_promotion": True,
                "authenticated_rhstiming_v1_rejects_r_zero": True,
            },
            "low_level_gf256_policy": {
                "timing_binary_wh_count": False,
                "rhstiming_v1_exposes_call_or_byte_counters": False,
                "separate_non_timing_diagnostic_required": True,
            },
            "binary_layout": {"binaries": layouts, "gate": layout_gate},
            "minor_fault_policy": _fault_policy(),
            "execution_failure_policy": _execution_failure_policy(engine),
            "sibling_idle_policy": engine.SIBLING_IDLE_POLICY,
            "core": args.core, "numa_node": args.numa_node,
            "topology": topology, "cpu_context": cpu_context,
            "controller_core": args.controller_core,
            "controller_topology": controller_topology,
            "evict_bytes": args.evict_bytes,
            "thermal_limits_c": {
                "cpu": engine.MAX_CPU_TEMP_C,
                "dimm": engine.MAX_DIMM_TEMP_C,
                "max_gap_s": engine.MAX_THERMAL_GAP_S,
                "max_coverage_margin_s": engine.MAX_THERMAL_MARGIN_S,
            },
            "thermal_sampler": {
                "name": engine.THERMAL_SAMPLER_NAME,
                "sha256": sha256_bytes(sampler_raw),
            },
            "fresh_only": True,
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "immutable_files": immutable, "tools": tool_records,
            "build_provenance_sha256": {
                label: engine.sha256_file(value["path"])
                for label, value in builds.items()
            },
            "prepare_smoke": smoke,
        }
        design = engine.sealed_record(DESIGN_SCHEMA, design_payload)
        design_path = staging / "design.json"
        engine.write_new(design_path, engine.canonical_json(design))
        receipt = engine.sealed_record(PREPARE_SCHEMA, {
            "prepared_utc": engine.utc_now(),
            "design_sha256": engine.sha256_file(design_path),
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "immutable_files": immutable,
            "base_binary_sha256":
                immutable["frozen/" + engine.BINARY_NAMES["base"]],
            "candidate_binary_sha256":
                immutable["frozen/" + engine.BINARY_NAMES["candidate"]],
        })
        engine.write_new(
            staging / "prepare_receipt.json", engine.canonical_json(receipt))
        for directory in (frozen, provenance):
            os.chmod(directory, 0o500)
            with engine.evidence_io.held_private_directory(
                    directory, require_writable=False,
                    error_type=engine.TimingError) as directory_fd:
                os.fsync(directory_fd)
        os.chmod(staging, 0o500)
        staging_sealed = True
        engine.evidence_io.publish_directory_noreplace(
            staging, result, error_type=engine.TimingError)
        staging_sealed = False
        os.chmod(result, 0o700)
        with engine.evidence_io.held_private_directory(
                result, require_writable=True,
                error_type=engine.TimingError) as result_fd:
            os.fsync(result_fd)
        with engine.evidence_io.held_private_directory(
                result.parent, require_writable=True,
                error_type=engine.TimingError) as parent_fd:
            os.fsync(parent_fd)
        print(json.dumps({
            "result_dir": str(result), "task_count": len(tasks),
            "design_sha256": engine.sha256_file(result / "design.json"),
            "manifest_sha256":
                engine.sha256_file(result / "tasks_manifest.jsonl"),
            "prepare_receipt_sha256":
                engine.sha256_file(result / "prepare_receipt.json"),
            "base_binary_sha256": receipt["base_binary_sha256"],
            "candidate_binary_sha256": receipt["candidate_binary_sha256"],
            "same_layout_gate": True,
            "frozen_controller":
                str(result / "frozen" / wrapper_frozen.name),
            "frozen_engine": str(result / "frozen" / formal_frozen.name),
        }, sort_keys=True))
    finally:
        if staging.exists() and not staging_sealed:
            shutil.rmtree(staging)
        if workspace.exists():
            shutil.rmtree(workspace)


_DESIGN_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "root", "base_commit",
    "candidate_commit", "codec_lineage", "measurement_overlay",
    "formal_engine", "controller", "timing_scope", "architecture",
    "cell_pairs", "K", "bb", "schedule_seeds", "cache_states",
    "overhead", "expect_q", "overhead_stream", "loss", "outer_order",
    "inner_order", "inner_cycles", "discard_inner_cycle", "task_count",
    "processes_per_task", "rows_per_process", "timed_rows_per_process",
    "timed_rows_per_binary_per_task", "allocator_environment",
    "phase_timing", "timed_binary_build", "r_population_policy",
    "low_level_gf256_policy", "binary_layout", "minor_fault_policy",
    "execution_failure_policy", "sibling_idle_policy", "core",
    "numa_node", "topology", "cpu_context", "controller_core",
    "controller_topology",
    "evict_bytes", "thermal_limits_c", "thermal_sampler", "fresh_only",
    "tasks_manifest_sha256", "immutable_files", "tools",
    "build_provenance_sha256", "prepare_smoke",
))


def _load_design(engine, root: Path) -> Dict[str, object]:
    design = engine.load_canonical(root / "design.json", "RHS timing design")
    engine.verify_sealed_record(design, DESIGN_SCHEMA, "RHS timing design")
    if set(design) != _DESIGN_FIELDS:
        raise engine.TimingError("RHS timing design fields changed")
    if (design.get("root") != str(root.resolve()) or
            design.get("base_commit") != BASE_COMMIT or
            design.get("candidate_commit") != CANDIDATE_COMMIT):
        raise engine.TimingError("RHS timing design identity changed")
    expected: Dict[str, object] = {
        "codec_lineage": {
            "candidate_parent_commit": BASE_COMMIT,
            "rhs_overlay_commit": RHS_OVERLAY_COMMIT,
            "rhs_overlay_parent_commit": RHS_OVERLAY_PARENT_COMMIT,
            "rhs_overlay_tree": RHS_OVERLAY_TREE,
            "overlay_is_measurement_only_and_branch_independent": True,
        },
        "formal_engine": {
            "commit": FORMAL_ENGINE_COMMIT,
            "engine_sha256": FORMAL_ENGINE_SHA256,
            "helper_sha256": FORMAL_HELPER_SHA256,
            "delegated_controls": [
                "bounded-process-groups", "outer-ABBABAAB",
                "CPU-NUMA-affinity", "IRQ-and-sibling-idle-gates",
                "sole-reader-CPU-DIMM-thermal-evidence", "tmpfs-identity",
                "bounded-post-parse-contamination-retries",
                "exact-raw-receipt-replay",
            ],
        },
        "timing_scope": "full-payload decoder precode solve",
        "architecture": engine.ARCHITECTURE,
        "cell_pairs": [list(item) for item in CELL_PAIRS],
        "K": list(engine.KS), "bb": list(engine.WIDTHS),
        "schedule_seeds": [list(item) for item in SCHEDULE_SEEDS],
        "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
        "expect_q": EXPECT_Q, "overhead_stream": OVERHEAD_STREAM,
        "loss": LOSS_TEXT, "outer_order": OUTER_ORDER,
        "inner_order": INNER_ORDER, "inner_cycles": 4,
        "discard_inner_cycle": 0, "task_count": 30,
        "processes_per_task": 8, "rows_per_process": 32,
        "timed_rows_per_process": 24,
        "timed_rows_per_binary_per_task": 96,
        "allocator_environment": {
            "MALLOC_MMAP_THRESHOLD_": engine.MALLOC_MMAP_THRESHOLD,
            "MALLOC_TRIM_THRESHOLD_": engine.MALLOC_TRIM_THRESHOLD,
            "purpose": (
                "retain and prefault solve-output pages in each process "
                "while preserving all timed decoder work"),
        },
        "phase_timing": {
            "primary": "residual_ns",
            "negative_controls": list(engine.NEGATIVE_CONTROL_PHASE_FIELDS),
            "fields": list(PHASE_FIELDS),
        },
        "timed_binary_build": {
            "rhs_timing_only": True, "mode": "rhstiming",
            "wh_count": False,
            "c_release_flags": TIMING_C_FLAGS_RELEASE,
            "cxx_release_flags": TIMING_CXX_FLAGS_RELEASE,
            "exe_linker_flags": TIMING_EXE_LINKER_FLAGS,
            "forbidden_symbols_absent": list(FORBIDDEN_TIMED_BINARY_SYMBOLS),
            "unrelated_modes_rejected": True,
        },
        "r_population_policy": {
            "report_r_zero_and_positive_separately": True,
            "both_required_for_promotion": True,
            "authenticated_rhstiming_v1_rejects_r_zero": True,
        },
        "low_level_gf256_policy": {
            "timing_binary_wh_count": False,
            "rhstiming_v1_exposes_call_or_byte_counters": False,
            "separate_non_timing_diagnostic_required": True,
        },
        "minor_fault_policy": _fault_policy(),
        "execution_failure_policy": _execution_failure_policy(engine),
        "sibling_idle_policy": engine.SIBLING_IDLE_POLICY,
        "thermal_limits_c": {
            "cpu": engine.MAX_CPU_TEMP_C, "dimm": engine.MAX_DIMM_TEMP_C,
            "max_gap_s": engine.MAX_THERMAL_GAP_S,
            "max_coverage_margin_s": engine.MAX_THERMAL_MARGIN_S,
        },
        "fresh_only": True,
    }
    for key, value in expected.items():
        if key == "sibling_idle_policy":
            engine.require_frozen_sibling_idle_policy(design.get(key))
        elif design.get(key) != value:
            raise engine.TimingError("RHS timing policy changed: %s" % key)
    overlay = design.get("measurement_overlay")
    if _overlay_identity(
            engine, {**overlay, "diff_evidence_name": "design.diff"}
            if isinstance(overlay, dict) else overlay) != overlay:
        raise engine.TimingError("RHS design overlay is malformed")
    for key in ("core", "controller_core", "numa_node", "evict_bytes"):
        value = design.get(key)
        if (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise engine.TimingError("RHS design integer is malformed: %s" % key)
    topology = design.get("topology")
    controller_topology = design.get("controller_topology")
    if (not isinstance(topology, dict) or
            not isinstance(controller_topology, dict) or
            design["evict_bytes"] < 4096 or
            design["controller_core"] in topology.get("llc_shared_cpus", [])):
        raise engine.TimingError("RHS timing isolation domain is malformed")
    if design.get("cpu_context") != _cpu_context(engine, int(design["core"])):
        raise engine.TimingError("RHS timing CPU/ISA context changed")
    immutable = design.get("immutable_files")
    controller = design.get("controller")
    if (not isinstance(immutable, dict) or not isinstance(controller, dict) or
            set(controller) != {"name", "sha256"} or
            controller.get("name") != "wh2_rhs_ablation_diagnostic.py" or
            controller.get("sha256") != immutable.get(
                "frozen/wh2_rhs_ablation_diagnostic.py")):
        raise engine.TimingError("RHS controller binding is malformed")
    sampler = design.get("thermal_sampler")
    if sampler != {
            "name": engine.THERMAL_SAMPLER_NAME,
            "sha256": immutable.get(
                "frozen/" + engine.THERMAL_SAMPLER_NAME)}:
        raise engine.TimingError("RHS thermal sampler binding changed")
    layout = design.get("binary_layout")
    if (not isinstance(layout, dict) or
            set(layout) != {"binaries", "gate"} or
            not isinstance(layout.get("binaries"), dict) or
            _layout_gate(engine, layout["binaries"]) != layout.get("gate")):
        raise engine.TimingError("RHS binary-layout receipt is malformed")
    smoke = design.get("prepare_smoke")
    smoke_task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 1280 and
        value["schedule"] == "burst" and value["cache_state"] == "warm")
    if (not isinstance(smoke, dict) or smoke.get("timing_evidence") is not False or
            smoke.get("task") != smoke_task or
            smoke.get("graph_trace_packet_output_work_identity_exact") is not True or
            smoke.get("recovery_identity_exact") is not True):
        raise engine.TimingError("RHS compatibility-smoke receipt is malformed")
    engine.require_sha256(
        design.get("tasks_manifest_sha256"), "task manifest hash")
    return design


def _validate_prepare_receipt(
    engine, root: Path, design: Mapping[str, object],
) -> Dict[str, object]:
    prepare = engine.load_canonical(
        root / "prepare_receipt.json", "RHS prepare receipt")
    engine.verify_sealed_record(prepare, PREPARE_SCHEMA, "RHS prepare receipt")
    if set(prepare) != engine.PREPARE_RECEIPT_FIELDS:
        raise engine.TimingError("RHS prepare receipt fields changed")
    immutable = design.get("immutable_files")
    if (prepare.get("design_sha256") !=
            engine.sha256_file(root / "design.json") or
            prepare.get("tasks_manifest_sha256") !=
            engine.sha256_file(root / "tasks_manifest.jsonl") or
            prepare.get("immutable_files") != immutable or
            prepare.get("base_binary_sha256") != immutable.get(
                "frozen/" + engine.BINARY_NAMES["base"]) or
            prepare.get("candidate_binary_sha256") != immutable.get(
                "frozen/" + engine.BINARY_NAMES["candidate"]) or
            not isinstance(prepare.get("prepared_utc"), str) or
            engine.UTC_RE.fullmatch(prepare["prepared_utc"]) is None):
        raise engine.TimingError("RHS prepare receipt binding changed")
    return prepare


def _verify_immutable(
    engine, root: Path, design: Mapping[str, object],
) -> None:
    engine._rhs_original_verify_immutable(root, design)
    wrapper = Path(__file__).resolve()
    formal = Path(str(engine.__file__)).resolve()
    helper = Path(str(engine.evidence_io.__file__)).resolve()
    expected_wrapper = root / "frozen" / "wh2_rhs_ablation_diagnostic.py"
    expected_formal = root / "frozen" / "wh2_grouped_commit_timing.py"
    expected_helper = root / "frozen" / "wh2_timing_evidence_io.py"
    if (wrapper != expected_wrapper or formal != expected_formal or
            helper != expected_helper or
            engine.sha256_file(wrapper) != design["controller"]["sha256"] or
            engine.sha256_file(formal) != FORMAL_ENGINE_SHA256 or
            engine.sha256_file(helper) != FORMAL_HELPER_SHA256):
        raise engine.TimingError(
            "campaign must run through its exact frozen wrapper and engine")
    tools = {
        name: Path(str(record["path"]))
        for name, record in design["tools"].items()
    }
    replayed: Dict[str, Mapping[str, object]] = {}
    for label in ("base", "candidate"):
        record = design["binary_layout"]["binaries"].get(label)
        if not isinstance(record, dict):
            raise engine.TimingError("binary-layout label is missing")
        nm_name = record.get("nm_evidence_name")
        elf_name = record.get("readelf_evidence_name")
        if (not isinstance(nm_name, str) or not isinstance(elf_name, str) or
                "/" in nm_name or "/" in elf_name):
            raise engine.TimingError("binary-layout evidence name changed")
        nm_raw = engine.stable_bytes(
            root / "provenance" / nm_name,
            max_bytes=16 * 1024 * 1024)
        elf_raw = engine.stable_bytes(
            root / "provenance" / elf_name,
            max_bytes=16 * 1024 * 1024)
        binary = root / "frozen" / engine.BINARY_NAMES[label]
        if (_layout_command(
                engine,
                (str(tools["nm"]), "-S", "-n", "-C", str(binary)),
                label + " replay nm") != nm_raw or
                _layout_command(
                    engine,
                    (str(tools["readelf"]), "-W", "-h", "-S", str(binary)),
                    label + " replay readelf") != elf_raw):
            raise engine.TimingError("live binary-layout evidence changed")
        parsed = _parse_layout(
            engine, nm_raw, elf_raw, engine.sha256_file(binary))
        expected = {
            key: value for key, value in record.items()
            if key not in {"nm_evidence_name", "readelf_evidence_name"}
        }
        if parsed != expected:
            raise engine.TimingError("binary-layout receipt does not replay")
        replayed[label] = record
    if _layout_gate(engine, replayed) != design["binary_layout"]["gate"]:
        raise engine.TimingError("same-layout gate does not replay")


def _replay_campaign(engine, root: Path) \
        -> Tuple[Dict[str, object], set[str]]:
    payload, expected = engine._rhs_original_replay_campaign(root)
    design = _load_design(engine, root)
    tasks = engine._load_tasks(root, design)
    population = {
        label: {
            "processes": {"R_eq_0": 0, "R_gt_0": 0},
            "all_rows": {"R_eq_0": 0, "R_gt_0": 0},
            "timed_rows": {"R_eq_0": 0, "R_gt_0": 0},
            "distinct_R": set(),
        }
        for label in ("base", "candidate")
    }
    cell_population: List[Dict[str, object]] = []
    for task in tasks:
        task_R: Dict[str, List[int]] = {"base": [], "candidate": []}
        for slot, marker in enumerate(OUTER_ORDER):
            label = "base" if marker == "A" else "candidate"
            name = execution_name(task, slot, label)
            raw = engine.stable_bytes(
                root / "raw" / name, max_bytes=2 * 1024 * 1024)
            try:
                parsed = parse_rhs_output(
                    raw, label, task, int(design["evict_bytes"]),
                    int(design["core"]))
            except DiagnosticError as exc:
                raise engine.TimingError(str(exc)) from exc
            process_R = int(parsed.preamble["preflight_inactivated"])
            process_class = "R_eq_0" if process_R == 0 else "R_gt_0"
            population[label]["processes"][process_class] += 1
            population[label]["distinct_R"].add(process_R)
            task_R[label].append(process_R)
            for index, row in enumerate(parsed.rows):
                row_R = int(row["inactivated"])
                row_class = "R_eq_0" if row_R == 0 else "R_gt_0"
                population[label]["all_rows"][row_class] += 1
                if index >= 8:
                    population[label]["timed_rows"][row_class] += 1
        if (any(len(values) != 4 or len(set(values)) != 1
                for values in task_R.values()) or
                task_R["base"][0] != task_R["candidate"][0]):
            raise engine.TimingError("per-cell R population changed by arm")
        cell_population.append({
            "task_id": task["task_id"], "K": task["K"], "bb": task["bb"],
            "schedule": task["schedule"],
            "cache_state": task["cache_state"],
            "R": task_R["base"][0],
            "processes_per_binary": 4,
            "timed_rows_per_binary": 96,
        })
    population_json: Dict[str, object] = {}
    both_classes = True
    for label, counts in population.items():
        row = dict(counts)
        row["distinct_R"] = sorted(counts["distinct_R"])
        population_json[label] = row
        both_classes = both_classes and all(
            int(counts[domain][kind]) > 0
            for domain in ("processes", "all_rows", "timed_rows")
            for kind in ("R_eq_0", "R_gt_0"))
    payload["r_population"] = {
        "by_binary": population_json,
        "by_cell": cell_population,
        "both_R_zero_and_positive_observed": both_classes,
        "promotion_coverage_complete": both_classes,
        "note": (
            "authenticated rhstiming v1 rejects R=0 before emission; this "
            "initial large-K timing grid therefore cannot establish R=0 "
            "performance without a separately authenticated fixture"),
    }
    payload["low_level_gf256_receipts"] = {
        "availability": "not-exposed-by-authenticated-rhstiming-v1",
        "timing_binary_wh_count": False,
        "gf256_calls": None, "gf256_bytes": None,
        "separate_non_timing_build_required": True,
        "timing_not_mixed_with_diagnostic_dispatch": True,
    }
    payload["binary_layout"] = design["binary_layout"]
    payload["timing_promotional"] = bool(
        payload.get("timing_promotional") is True and both_classes)
    return payload, expected


def reduce_campaign(engine, args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    engine._require_external_prepare_anchor(
        root, args.expected_prepare_sha256)
    payload, expected = _replay_campaign(engine, root)
    summary = engine.sealed_record(SUMMARY_SCHEMA, payload)
    summary_path = root / "validated_summary.json"
    publication = {
        "validated_summary.json": engine._publish_or_verify_reduction_artifact(
            summary_path, engine.canonical_json(summary)),
    }
    expected.add("validated_summary.json")
    manifest = engine._data_manifest(root, expected)
    manifest_path = root / "data_manifest.json"
    publication["data_manifest.json"] = \
        engine._publish_or_verify_reduction_artifact(
            manifest_path, engine.canonical_json({"files": manifest}))
    publication["data_manifest.sha256"] = \
        engine._publish_or_verify_reduction_artifact(
            root / "data_manifest.sha256",
            (engine.sha256_file(manifest_path) + "\n").encode("ascii"))
    print(json.dumps({
        "summary_sha256": engine.sha256_file(summary_path),
        "data_manifest_sha256": engine.sha256_file(manifest_path),
        "overall_ratio": summary["overall"]["ratio_of_sums"]["decimal"],
        "residual_ratio": summary["phase_timing"]["metrics"][
            "residual_ns"]["overall"]["ratio_of_sums"]["decimal"],
        "work_and_recovery_exact": True,
        "timing_promotional": summary["timing_promotional"],
        "r_population": summary["r_population"],
        "publication": publication,
    }, sort_keys=True))


def verify_campaign(engine, args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    engine._require_external_prepare_anchor(
        root, args.expected_prepare_sha256)
    payload, expected = _replay_campaign(engine, root)
    summary = engine.load_canonical(
        root / "validated_summary.json", "RHS validated summary")
    engine.verify_sealed_record(
        summary, SUMMARY_SCHEMA, "RHS validated summary")
    if summary != engine.sealed_record(SUMMARY_SCHEMA, payload):
        raise engine.TimingError("RHS validated summary does not reproduce")
    expected.add("validated_summary.json")
    manifest = engine.load_canonical(
        root / "data_manifest.json", "RHS data manifest")
    if (set(manifest) != {"files"} or
            manifest["files"] != engine._data_manifest(root, expected)):
        raise engine.TimingError("RHS data manifest does not reproduce")
    try:
        sidecar = engine.stable_bytes(
            root / "data_manifest.sha256", max_bytes=128).decode(
                "ascii", "strict")
    except UnicodeDecodeError as exc:
        raise engine.TimingError("RHS manifest sidecar is not ASCII") from exc
    if sidecar != engine.sha256_file(root / "data_manifest.json") + "\n":
        raise engine.TimingError("RHS manifest sidecar changed")
    actual = {
        str(path.relative_to(root))
        for path in root.rglob("*") if path.is_file()
    }
    allowed = set(expected) | {"data_manifest.json", "data_manifest.sha256"}
    if actual != allowed:
        raise engine.TimingError("RHS campaign file inventory changed")
    print(json.dumps({
        "verified": True,
        "summary_sha256":
            engine.sha256_file(root / "validated_summary.json"),
        "data_manifest_sha256":
            engine.sha256_file(root / "data_manifest.json"),
        "file_count": len(actual),
        "timing_promotional": summary["timing_promotional"],
    }, sort_keys=True))


def _synthetic_rhs_output(
    task: Mapping[str, object], core: int,
) -> bytes:
    K = int(task["K"])
    width = int(task["bb"])
    R = 20
    peeled = K + 16
    intermediate = (K + 36) * width
    residual_rows = 32
    unused = residual_rows - 12
    hashes = {
        "graph_sha256": "1" * 64,
        "trace_sha256": "2" * 64,
        "packet_prefix_sha256": "3" * 64,
        "payload_sha256": "4" * 64,
        "preflight_output_sha256": "5" * 64,
    }
    preamble = {
        **_expected_preamble(task, 4096),
        "S": "12", "base_matrix_seed": "0x1",
        "base_peel_seed": "0x2", "selected_attempt": "0",
        "preflight_packet_seed_attempt": "0", "matrix_seed": "0x1",
        "peel_seed": "0x2", **hashes,
        "preflight_intermediate_bytes": str(intermediate),
        "preflight_peeled_columns": str(peeled),
        "preflight_inactivated": str(R),
        "preflight_residual_rows": str(residual_rows),
        "preflight_residual_rank": "20",
        "preflight_binary_residual_rank": "12",
        "preflight_unused_binary_rows": str(unused),
        "preflight_rhs_init_destination_bytes": str(unused * width),
    }
    if set(preamble) != set(RHS_PREAMBLE_FIELDS):
        raise DiagnosticError("synthetic preamble fields changed")
    lines = [
        "# rhstiming: " + " ".join(
            "%s=%s" % (field, preamble[field])
            for field in RHS_PREAMBLE_FIELDS),
        ",".join(RHS_FIELDS),
    ]
    for index in range(32):
        cycle, slot = divmod(index, 8)
        values = {
            "N": str(K), "bb": str(width), "overhead": str(OVERHEAD),
            "expect_q": str(EXPECT_Q),
            "overhead_stream": OVERHEAD_STREAM,
            "schedule": str(task["schedule"]), "seed": str(task["seed"]),
            "loss": LOSS_TEXT, "cache_state": str(task["cache_state"]),
            "cycle": str(cycle), "slot": str(slot),
            "arm": "control" if INNER_ORDER[slot] == "A" else "candidate",
            "field": "gf256", "completion": "gf256", "H": "12",
            "D2": "12", "S": "12", "source_hits": "2", "mix": "3",
            "dense_two_anchor": "0", "heavy_family": "periodic",
            "seed_attempt": "0", "packet_seed_attempt": "0",
            "matrix_seed": "0x1", "peel_seed": "0x2",
            "graph_sha256": hashes["graph_sha256"],
            "trace_sha256": hashes["trace_sha256"],
            "packet_prefix_sha256": hashes["packet_prefix_sha256"],
            "payload_contract": "distinct-packet-zero-v1",
            "payload_sha256": hashes["payload_sha256"],
            "preflight_result": "0",
            "output_sha256": hashes["preflight_output_sha256"],
            "output_validation": "exact-all-zero-full-scan-v1",
            "all_zero": "1", "cell_class": "common-success",
            "common_success": "1", "result": "0",
            "outcome_stable": "1", "elapsed_ns": "1000",
            "saturated": "0", "cpu_before": str(core),
            "cpu_after": str(core), "cpu_migrated": "0",
            "minflt_delta": "0", "majflt_delta": "0",
            "fault_contaminated": "0", "packet_rows": str(K + OVERHEAD),
            "peeled_columns": str(peeled), "inactivated": str(R),
            "residual_rows": str(residual_rows), "residual_rank": "20",
            "binary_residual_rank": "12", "q": "8", "heavy_gain": "8",
            "unused_binary_rows": str(unused),
            "rhs_init_destination_bytes": str(unused * width),
            "binary_row_references": "100", "binary_row_storage_bytes": "200",
            "binary_adjacency_storage_bytes": "300",
            "binary_row_storage_allocations": "1",
            "binary_adjacency_storage_allocations": "1",
            "block_xors": "123", "block_muladds": "456",
            "build_ns": "100", "peel_ns": "100", "project_ns": "100",
            "residual_ns": "500", "backsub_ns": "100",
            "phase_sum_ns": "900", "mixed_joint_source_xors": "0",
            "mixed_joint_marginal_xors": "0",
            "mixed_joint_marginal_copies": "0",
            "mixed_joint_active_deltas": "0",
            "mixed_joint_scratch_bytes": "0",
            "mixed_dual_source_columns": "0",
            "source_bytes": str(K * width),
            "packet_payload_bytes": str((K + OVERHEAD) * width),
            "intermediate_bytes": str(intermediate),
        }
        if set(values) != set(RHS_FIELDS):
            raise DiagnosticError("synthetic row fields changed")
        lines.append(",".join(values[field] for field in RHS_FIELDS))
    return ("\n".join(lines) + "\n").encode("ascii")


def _expect_diagnostic_failure(function, message: str) -> None:
    try:
        function()
    except DiagnosticError:
        return
    raise DiagnosticError("selftest accepted %s" % message)


def _poison_source(size: int, marker: Path, *, engine: bool) -> bytes:
    lines = [
        "from pathlib import Path",
        "Path(%r).write_text('poisoned')" % str(marker),
    ]
    if engine:
        lines.extend((
            "import wh2_timing_evidence_io as evidence_io",
            "BASE_COMMIT='48d14bc77e3f9e98605fca4d226aa218d7d03a0d'",
            "CANDIDATE_COMMIT=%r" % BASE_COMMIT,
            "OUTER_ORDER=%r" % OUTER_ORDER,
        ))
    prefix = ("\n".join(lines) + "\n#").encode("ascii")
    if len(prefix) > size:
        raise DiagnosticError("poisoned loader fixture exceeds source size")
    return prefix + b"x" * (size - len(prefix))


def _install_poisoned_pyc(
    source_path: Path, authenticated: bytes, marker: Path, *, engine: bool,
) -> None:
    poisoned = _poison_source(len(authenticated), marker, engine=engine)
    source_path.write_bytes(poisoned)
    timestamp = 1700000000
    os.utime(str(source_path), (timestamp, timestamp))
    cached = py_compile.compile(str(source_path), doraise=True)
    if cached is None or not Path(cached).is_file():
        raise DiagnosticError("could not create poisoned-pyc fixture")
    source_path.write_bytes(authenticated)
    os.utime(str(source_path), (timestamp, timestamp))


def _authenticated_loader_selftest(engine_path: Path) -> None:
    """Prove byte execution ignores pyc/shadows and rejects source swaps."""
    engine_path = engine_path.resolve()
    helper_path = engine_path.with_name("wh2_timing_evidence_io.py")
    engine_raw = _read_plain_file(engine_path, 1024 * 1024)
    helper_raw = _read_plain_file(helper_path, 1024 * 1024)
    if (sha256_bytes(engine_raw) != FORMAL_ENGINE_SHA256 or
            sha256_bytes(helper_raw) != FORMAL_HELPER_SHA256):
        raise DiagnosticError("loader selftest source identity mismatch")
    previous_engine = sys.modules.get(
        _FORMAL_ENGINE_MODULE_NAME, _MISSING_MODULE)
    previous_helper = sys.modules.get(
        _FORMAL_HELPER_MODULE_NAME, _MISSING_MODULE)
    try:
        with tempfile.TemporaryDirectory(
                prefix="wh2-rhs-loader-selftest-") as temporary:
            root = Path(temporary)

            pyc_root = root / "poisoned-pyc"
            pyc_root.mkdir()
            pyc_engine = pyc_root / engine_path.name
            pyc_helper = pyc_root / helper_path.name
            engine_marker = root / "engine-pyc-executed"
            helper_marker = root / "helper-pyc-executed"
            _install_poisoned_pyc(
                pyc_engine, engine_raw, engine_marker, engine=True)
            _install_poisoned_pyc(
                pyc_helper, helper_raw, helper_marker, engine=False)
            loaded = load_formal_engine(pyc_engine)
            if (engine_marker.exists() or helper_marker.exists() or
                    loaded.__file__ != str(pyc_engine.resolve()) or
                    loaded.evidence_io.__file__ !=
                        str(pyc_helper.resolve()) or
                    loaded.BoundCommandResult.__module__ !=
                        _FORMAL_ENGINE_MODULE_NAME):
                raise DiagnosticError("authenticated loader used poisoned pyc")

            swap_root = root / "source-swap"
            swap_root.mkdir()
            swap_engine = swap_root / engine_path.name
            swap_helper = swap_root / helper_path.name
            swap_engine.write_bytes(engine_raw)
            swap_helper.write_bytes(helper_raw)
            swap_marker = root / "swapped-source-executed"
            original_reader = globals()["_read_plain_file"]
            swapped = [False]

            def swapping_reader(path: Path, maximum: int) -> bytes:
                raw = original_reader(path, maximum)
                if Path(path).resolve() == swap_engine.resolve() and \
                        not swapped[0]:
                    swap_engine.write_bytes(
                        _poison_source(512, swap_marker, engine=True))
                    swapped[0] = True
                return raw

            sentinel_engine = types.ModuleType("prior_engine_sentinel")
            sentinel_helper = types.ModuleType("prior_helper_sentinel")
            sys.modules[_FORMAL_ENGINE_MODULE_NAME] = sentinel_engine
            sys.modules[_FORMAL_HELPER_MODULE_NAME] = sentinel_helper
            globals()["_read_plain_file"] = swapping_reader
            try:
                _expect_diagnostic_failure(
                    lambda: load_formal_engine(swap_engine),
                    "authenticated source swap")
            finally:
                globals()["_read_plain_file"] = original_reader
            if (not swapped[0] or swap_marker.exists() or
                    sys.modules.get(_FORMAL_ENGINE_MODULE_NAME) is not
                        sentinel_engine or
                    sys.modules.get(_FORMAL_HELPER_MODULE_NAME) is not
                        sentinel_helper):
                raise DiagnosticError(
                    "source-swap failure left partial authenticated modules")

            shadow_root = root / "shadow"
            exact_root = root / "shadow-target"
            shadow_root.mkdir()
            exact_root.mkdir()
            shadow_marker = root / "shadow-helper-executed"
            (shadow_root / helper_path.name).write_bytes(
                _poison_source(512, shadow_marker, engine=False))
            exact_engine = exact_root / engine_path.name
            exact_helper = exact_root / helper_path.name
            exact_engine.write_bytes(engine_raw)
            exact_helper.write_bytes(helper_raw)
            original_path = list(sys.path)
            sys.path.insert(0, str(shadow_root))
            loader_path = list(sys.path)
            try:
                loaded = load_formal_engine(exact_engine)
                if sys.path != loader_path:
                    raise DiagnosticError(
                        "authenticated loader changed the import path")
            finally:
                sys.path[:] = original_path
            if (shadow_marker.exists() or
                    loaded.evidence_io.__file__ !=
                        str(exact_helper.resolve()) or
                    sha256_bytes(_read_plain_file(
                        Path(loaded.evidence_io.__file__), 1024 * 1024)) !=
                        FORMAL_HELPER_SHA256):
                raise DiagnosticError("authenticated loader used shadow helper")
    finally:
        _restore_module(_FORMAL_ENGINE_MODULE_NAME, previous_engine)
        _restore_module(_FORMAL_HELPER_MODULE_NAME, previous_helper)


def selftest(engine_path: Optional[str] = None) -> None:
    _formal_helper_digest_selftest()
    if engine_path is not None:
        _authenticated_loader_selftest(Path(engine_path))
        engine = load_formal_engine(Path(engine_path))
        _configure_engine(engine)
        if (engine.BASE_COMMIT != BASE_COMMIT or
                engine.CANDIDATE_COMMIT != CANDIDATE_COMMIT or
                engine.MAX_MINOR_FAULTS != 0 or
                _execution_failure_policy(engine) != {
                    "post_parse_contamination_max_attempts": 10,
                    "native_nonzero_exit_or_stderr":
                        "substantive-campaign-fatal",
                    "native_pre_emission_rejection":
                        "substantive-campaign-fatal-no-retry",
                } or
                len(engine.generate_tasks()) != 30):
            raise DiagnosticError("formal engine specialization failed")
    tasks = generate_tasks()
    if (len(tasks) != 30 or tuple(task["job"] for task in tasks) !=
            tuple(range(30)) or
            {(task["K"], task["bb"]) for task in tasks} != set(CELL_PAIRS)):
        raise DiagnosticError("task-grid selftest failed")
    task = next(
        value for value in tasks
        if value["K"] == 3200 and value["bb"] == 1280 and
        value["schedule"] == "burst" and value["cache_state"] == "warm")
    raw = _synthetic_rhs_output(task, 6)
    parsed = parse_rhs_output(raw, "base", task, 4096, 6)
    if (parsed.timed_elapsed_ns != 24000 or
            parsed.all_elapsed_ns != 32000 or parsed.contaminations or
            _fault_policy() != {
                "minimum": 0, "maximum": 0, "major_faults": 0,
                "all_cycles_receipted": True,
                "unavailable_or_nonzero_is_contamination": True,
            }):
        raise DiagnosticError("RHS parser aggregation selftest failed")
    lines = raw.decode("ascii").splitlines()
    fault_header = lines[1].split(",")
    fault_row = lines[2].split(",")
    fault_row[fault_header.index("minflt_delta")] = "1"
    fault_row[fault_header.index("fault_contaminated")] = "1"
    lines[2] = ",".join(fault_row)
    faulted = parse_rhs_output(
        ("\n".join(lines) + "\n").encode("ascii"),
        "base", task, 4096, 6)
    if (faulted.max_minor_faults != 1 or
            faulted.contaminations != ("row0:minor-fault:1",)):
        raise DiagnosticError("exact-zero fault policy selftest failed")
    lines = raw.decode("ascii").splitlines()
    fault_header = lines[1].split(",")
    fault_row = lines[2].split(",")
    fault_row[fault_header.index("majflt_delta")] = "1"
    fault_row[fault_header.index("fault_contaminated")] = "1"
    lines[2] = ",".join(fault_row)
    major_faulted = parse_rhs_output(
        ("\n".join(lines) + "\n").encode("ascii"),
        "base", task, 4096, 6)
    if major_faulted.contaminations != ("row0:major-fault:1",):
        raise DiagnosticError("exact-zero major-fault selftest failed")
    lines = raw.decode("ascii").splitlines()
    header = lines[1].split(",")
    trace_index = header.index("trace_sha256")
    changed = lines[2].split(",")
    changed[trace_index] = "6" * 64
    lines[2] = ",".join(changed)
    bad_trace = ("\n".join(lines) + "\n").encode("ascii")
    _expect_diagnostic_failure(
        lambda: parse_rhs_output(bad_trace, "base", task, 4096, 6),
        "one-row trace mutation")
    lines = raw.decode("ascii").splitlines()
    block_index = lines[1].split(",").index("block_xors")
    changed = lines[-1].split(",")
    changed[block_index] = "124"
    lines[-1] = ",".join(changed)
    bad_work = ("\n".join(lines) + "\n").encode("ascii")
    _expect_diagnostic_failure(
        lambda: parse_rhs_output(bad_work, "base", task, 4096, 6),
        "one-row work mutation")
    r0 = raw.replace(b"preflight_inactivated=20", b"preflight_inactivated=0", 1)
    _expect_diagnostic_failure(
        lambda: parse_rhs_output(r0, "base", task, 4096, 6),
        "authenticated-fixture R=0")
    fake_design = {
        "root": "/sealed", "core": 6, "numa_node": 0,
        "evict_bytes": 4096,
        "tools": {
            name: {"path": "/usr/bin/" + name}
            for name in ("env", "taskset", "numactl")
        },
    }
    argv = command_for(fake_design, task, "candidate")
    required = {
        "rhstiming", "--N", "3200", "--bb", "1280", "--expect-q", "8",
        "--overhead-stream", "salted", "--cache-state", "warm",
        "--schedule", "burst",
    }
    if not required <= set(argv) or "--cycle-index" in argv:
        raise DiagnosticError("RHS argv selftest failed")
    ledger: Dict[Tuple[object, ...], Dict[str, Dict[str, object]]] = {}
    _register_cross_cache_identity(ledger, task, parsed)
    cold_task = dict(task)
    cold_task["cache_state"] = "cold"
    cold = parse_rhs_output(
        _synthetic_rhs_output(cold_task, 6), "candidate", cold_task, 4096, 6)
    _register_cross_cache_identity(ledger, cold_task, cold)
    if len(ledger) != 1:
        raise DiagnosticError("cross-cache identity selftest failed")
    synthetic_source = (
        _USAGE_NEEDLE + b"usage-body\n#endif\n" +
        _DISPATCH_NEEDLE + b"compare-body\n" +
        b"        if (!std::strcmp(argv[1], \"rhstiming\")) {\n"
        b"            return CmdRhsTiming(argc - 2, argv + 2);\n"
        b"        }\n" + _CATCH_NEEDLE)
    transformed = rhs_timing_only_overlay(synthetic_source)
    if transformed == synthetic_source or b"#else" not in transformed:
        raise DiagnosticError("timing-only overlay selftest failed")
    _expect_diagnostic_failure(
        lambda: rhs_timing_only_overlay(
            synthetic_source.replace(_DISPATCH_NEEDLE, b"changed\n", 1)),
        "timing-only overlay anchor mutation")
    class LayoutEngine:
        TimingError = DiagnosticError
    layout_engine = LayoutEngine()
    nm = ("".join(
        "%016x 0000000000000020 T %s\n" % (0x1000 + index * 0x40, name)
        for index, name in enumerate(_GF256_SYMBOL_PREFIXES)) +
        "0000000000001180 0000000000000020 t prior_symbol\n" +
        "0000000000001200 0000000000000100 t " + _SOLVE_SYMBOL +
        "args, bool)\n").encode("utf-8")
    elf = (
        "ELF Header:\n  Class: ELF64\n"
        "  Data: 2's complement, little endian\n"
        "  Type: DYN (Position-Independent Executable file)\n"
        "  Machine: Advanced Micro Devices X86-64\n"
        "Section Headers:\n"
        "  [16] .text PROGBITS 0000000000001000 001000 003000 00 AX 0 0 64\n"
    ).encode("utf-8")
    base_layout = _parse_layout(layout_engine, nm, elf, "a" * 64)
    candidate_layout = json.loads(json.dumps(base_layout))
    candidate_layout["binary_sha256"] = "b" * 64
    candidate_layout["solve_symbol"]["size"] = 240
    candidate_layout[".text"]["size"] = 12016
    gate = _layout_gate(
        layout_engine, {"base": base_layout, "candidate": candidate_layout})
    if gate["solve_size_delta"] != -16 or gate["passed"] is not True:
        raise DiagnosticError("same-layout gate selftest failed")
    broken_layout = json.loads(json.dumps(candidate_layout))
    broken_layout["solve_section_offset"] += 16
    broken_layout["solve_symbol"]["address"] += 16
    broken_layout["solve_alignment_mod_64"] = \
        broken_layout["solve_symbol"]["address"] % 64
    _expect_diagnostic_failure(
        lambda: _layout_gate(
            layout_engine,
            {"base": base_layout, "candidate": broken_layout}),
        "solve-layout mutation")
    print(json.dumps({
        "selftest": "PASS", "task_count": len(tasks),
        "parser_rows": len(parsed.rows), "same_layout_gate": True,
        "formal_engine_imported": engine_path is not None,
        "authenticated_loader_adversarial": engine_path is not None,
        "formal_engine_sha256": FORMAL_ENGINE_SHA256,
        "formal_helper_sha256": FORMAL_HELPER_SHA256,
    }, sort_keys=True))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--engine",
        help=(
            "exact 5cb8d46 wh2_grouped_commit_timing.py; required except "
            "for selftest"))
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "selftest", help="run bounded synthetic fail-closed tests")
    prepare = subparsers.add_parser(
        "prepare", help="build and freeze the exact paired diagnostic")
    prepare.add_argument("--repo", required=True)
    prepare.add_argument("--result-dir", required=True)
    prepare.add_argument("--thermal-sampler", required=True)
    prepare.add_argument("--core", type=int, required=True)
    prepare.add_argument("--controller-core", type=int, required=True)
    prepare.add_argument("--numa-node", type=int, required=True)
    prepare.add_argument("--evict-bytes", type=int, default=256 * 1024 * 1024)
    prepare.add_argument(
        "--build-jobs", type=int, default=max(1, os.cpu_count() or 1))
    prepare.add_argument("--c-compiler")
    prepare.add_argument("--cxx-compiler")
    run = subparsers.add_parser(
        "run", help="run one fresh campaign after competing loads stop")
    run.add_argument("--result-dir", required=True)
    run.add_argument("--thermal-csv", required=True)
    run.add_argument("--thermal-pid-file", required=True)
    run.add_argument("--expected-prepare-sha256", required=True)
    run.add_argument("--timeout-seconds", type=int, default=1800)
    reduce = subparsers.add_parser(
        "reduce", help="strictly replay and publish a summary")
    reduce.add_argument("--result-dir", required=True)
    reduce.add_argument("--expected-prepare-sha256", required=True)
    verify = subparsers.add_parser(
        "verify", help="replay a reduced campaign")
    verify.add_argument("--result-dir", required=True)
    verify.add_argument("--expected-prepare-sha256", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.command == "selftest":
        try:
            selftest(args.engine)
        except (OSError, ValueError, DiagnosticError) as exc:
            print("RHS diagnostic selftest failed: %s" % exc,
                  file=sys.stderr, flush=True)
            return 2
        return 0
    if not args.engine:
        parser.error("--engine is required for prepare/run/reduce/verify")
    try:
        engine = load_formal_engine(Path(args.engine))
        _configure_engine(engine)
        if getattr(args, "timeout_seconds", 1) <= 0:
            raise engine.TimingError("timeout must be positive")
        result = Path(args.result_dir).resolve()
        lock_name = ".wirehair-wh2-rhs-ablation-%s.lock" % sha256_bytes(
            str(result).encode("utf-8"))[:24]
        with engine.evidence_io.nonblocking_global_flock(
                result.parent / lock_name, error_type=engine.TimingError):
            if args.command == "prepare":
                prepare_campaign(engine, args)
            elif args.command == "run":
                engine.run_campaign(args)
            elif args.command == "reduce":
                reduce_campaign(engine, args)
            elif args.command == "verify":
                verify_campaign(engine, args)
            else:
                raise engine.TimingError("unknown RHS diagnostic command")
    except (OSError, ValueError, DiagnosticError) as exc:
        print("RHS ablation diagnostic failed: %s" % exc,
              file=sys.stderr, flush=True)
        return 2
    except Exception as exc:
        # The authenticated engine's TimingError type is available only after
        # its exact file is loaded.  Keep argparse/import failures above narrow,
        # then preserve the engine's concise fail-closed CLI behavior here.
        if "engine" in locals() and isinstance(exc, engine.TimingError):
            print("RHS ablation diagnostic failed: %s" % exc,
                  file=sys.stderr, flush=True)
            return 2
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
