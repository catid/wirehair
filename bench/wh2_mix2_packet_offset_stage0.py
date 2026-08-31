#!/usr/bin/env python3
"""Sealed fixed off-diagonal MIX2 packet-attempt Stage 0 screen.

This controller keeps the v9 precode attempt on the three preregistered weak
coordinates and scans one global nonzero packet-attempt offset.  It consumes
no new root, collects no timing evidence, and cannot promote a profile.  The
lowest surviving offset only earns a separately frozen consumed-root Stage 1
replay.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage0-contract.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage0-record.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage0-summary.v1"
BENCH_DESCRIPTION_SCHEMA = "wirehair.wh2.v2-bench-description.v1"
RECORD_NAME = "mix2-packet-offset-stage0-records.jsonl"
SUMMARY_NAME = "mix2-packet-offset-stage0-summary.json"

MIX_PAIR = "01"
DELTAS = tuple(range(1, 256))
OFFSETS = (0,) + DELTAS
PACKET_ATTEMPT_STRIDE = 0x9e3779b9
REPETITIONS = (0, 1)
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
BINARY_DENSE_ROWS = 12
GF256_HEAVY_ROWS = 12
ZERO_SHA256 = "0" * 64
EXPECTED_CANDIDATE_PROFILE_SHA256 = (
    "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0")
EXPECTED_V9_CONTRACT_SHA256 = (
    "28e6affc80680377df2e8099825b1b475b65e51d04116a8c2fd6465566041cdb")
EXPECTED_V9_SOURCE_COMMIT = "52177d5d798935e13e83f01a436bd1b7198e80ed"
EXPECTED_V9_WORKER_SHA256 = (
    "23e2548804fe9cc7f30b2864429f37f48ac6f6e92c306f61533be51c0454481b")
EXPECTED_V9_INPUT_SELF_SHA256 = (
    "5c7f3309d3a05a458adbcbe1a50ff25e3b70c829ae9806b7592c9d188a1840c4")
EXPECTED_V9_SUMMARY_SELF_SHA256 = (
    "bbbf7e357e9d5cd84b9c79d6d91d0109f83609763452588e264e6cf38a200f4e")
V9_FILES = {
    "promotion-short-screen-attempts.jsonl": (
        281208,
        "fc4efa00fb152e2bd442eca31d3f381af0c0b1d4468fa2225ac606a33162d315"),
    "promotion-short-screen-input.json": (
        897281,
        "d75bd30efe62eed134dd95097759c8d84b47c5e5889226b236079f956e1e85da"),
    "promotion-short-screen-results.jsonl": (
        1805192,
        "2d8b0269a8b103daaf8d4736c3c73c2e4bee2dd33d0f69eeeda6938e7070b371"),
    "promotion-short-screen-summary.json": (
        5503,
        "d90cf1912fb250a4a00ca63acbdb637f7df475ca856db68d74285b4032f5c2f8"),
}
V9_K_VALUES = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512, 513,
    1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000, 10001, 16384,
    20000, 20001, 32768, 49152, 64000,
)
V9_ATTEMPTS = (
    9, 8, 0, 2, 49, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
)

# Selection roots 0..8 followed by the already-consumed v9 roots 9..11.
CONSUMED_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
    "0xc0ac29b7c97c50dd",
    "0x3f84d5b5b5470917",
    "0x9216d5d98979fb1b",
    "0xb889883a79549774",
    "0xb5666de0987896af",
    "0x8bfca269b0bc01e0",
    "0xc4695292d9835286",
    "0x7ccd510f122fc160",
    "0x7001a960b7d9c0a4",
)
UNCONSUMED_FINAL_ROOTS = (
    "0xefd20c982041a46b",
    "0x8827bc36ed906555",
    "0x86029f23d6132efa",
)
STAIRCASE_BY_K = {2: 2, 1000: 50, 10000: 86}


@dataclass(frozen=True)
class Coordinate:
    K: int
    attempt: int
    root_index: int
    schedule: str

    @property
    def root(self) -> str:
        return CONSUMED_ROOTS[self.root_index]

    def description(self, ordinal: int) -> Mapping[str, Any]:
        return {
            "coordinate_ordinal": ordinal,
            "K": self.K,
            "attempt": self.attempt,
            "root_index": self.root_index,
            "loss_root": self.root,
            "schedule": self.schedule,
        }


COORDINATES = (
    Coordinate(2, 9, 10, "burst"),
    Coordinate(1000, 0, 11, "burst"),
    Coordinate(10000, 0, 9, "adversarial"),
)
EXPECTED_INVOCATION_COUNT = (
    len(OFFSETS) * len(COORDINATES) * len(REPETITIONS))

CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_identity_corner", "overhead", "trials", "success",
    "rank_fail", "error", "fail_rate", "inact_mu", "inact_max",
    "binary_def_mu", "binary_def_max", "heavy_gain_mu", "heavy_gain_min",
    "heavy_shortfall", "solve_ms_mu", "build_ms_mu", "peel_ms_mu",
    "project_ms_mu", "residual_ms_mu", "backsub_ms_mu", "seed_attempt",
    "block_xors_mu", "block_muladds_mu", "first_rank_fail",
    "binary_def_hist", "heavy_gain_hist", "failure_trials",
    "active_packet_peel_seed_xor", "precode_attempt", "packet_attempt",
    "attempt_mode", "construction_seed_basis", "seed_schedule_sha256",
    "effective_precode_seed", "effective_packet_seed",
)
METADATA_ORDER = (
    "trials", "threads", "loss", "seed", "source_hits_override",
    "packet_peel_seed_xor", "binary_dense_rows_override",
    "gf256_heavy_rows_override", "dense_anchor_layout",
    "odd_packet_peel_seed_xor", "packet_row_seed_multiplier",
    "packet_row_seed_avalanche", "seed_block_bytes_override",
    "overhead_stream", "full_payload_solve", "schedule",
    "cold_solve_wide_xor", "exact_attempt_mode",
    "exact_precode_attempt", "exact_packet_attempt",
    "construction_seed_basis", "seed_schedule_sha256",
    "source_git_commit", "mix_pair",
)

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent

# Files that directly define, build, or configure wirehair_v2_bench and its
# linked libraries, plus this controller.  The binary hash is separately pinned.
SOURCE_PATHS = (
    "CMakeLists.txt",
    "codec/CMakeLists.txt",
    "bench/test_wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_promotion_short_screen.py",
    "wirehair.cpp",
    "include/wirehair/wirehair.h",
    "include/wirehair/wirehair.hpp",
    "gf256.cpp",
    "gf256.h",
    "WirehairCodec.cpp",
    "WirehairCodec.h",
    "WirehairCompiler.h",
    "WirehairEnvironment.h",
    "WirehairHeavy.h",
    "WirehairTools.cpp",
    "WirehairTools.h",
    "WirehairDenseFixups.inc",
    "WirehairPeelFixups.inc",
    "codec/WirehairV2Bench.cpp",
    "codec/WirehairV2Codec.cpp",
    "codec/WirehairV2Codec.h",
    "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Peel.h",
    "codec/WirehairV2Plan.cpp",
    "codec/WirehairV2Plan.h",
    "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Policy.h",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h",
    "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeDecode.h",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2PrecodeEncode.h",
    "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2RawSeed.h",
    "codec/WirehairV2Seeds.cpp",
    "codec/WirehairV2Seeds.h",
    "codec/WirehairV2Solve.cpp",
    "codec/WirehairV2Solve.h",
)
SOURCE_STATUS_PATHS = (
    "CMakeLists.txt", "codec", "include", "wirehair.cpp", "gf256.cpp",
    "gf256.h", "WirehairCodec.cpp", "WirehairCodec.h",
    "WirehairCompiler.h", "WirehairEnvironment.h", "WirehairHeavy.h",
    "WirehairTools.cpp", "WirehairTools.h", "WirehairDenseFixups.inc",
    "WirehairPeelFixups.inc",
    "bench/test_wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_promotion_short_screen.py",
)

MAX_STDOUT_BYTES = 256 * 1024
MAX_STDERR_BYTES = 64 * 1024
PROCESS_TIMEOUT_SECONDS = 600
CHILD_ENVIRONMENT = {"LANG": "C", "LC_ALL": "C", "TZ": "UTC"}
UINT64_MAX = (1 << 64) - 1

COMMIT = re.compile(r"[0-9a-f]{40}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
UINT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
THREE_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")
EIGHT_DECIMAL = re.compile(r"(?:0|1)\.[0-9]{8}\Z")
HEX64 = re.compile(r"0x[0-9a-f]{16}\Z")
HEX32 = re.compile(r"0x[0-9a-f]{8}\Z")

DIAGONAL_CONTROL_FIELDS = (
    "staircase", "source_hits", "outcome", "success", "rank_fail", "error",
    "inactivated_columns", "effective_precode_seed", "effective_packet_seed",
    "block_xors", "gf256_muladds",
)

# Replaced once after the contract body is finalized.  A literal digest makes
# a later roster or gate edit fail closed instead of blessing itself.
EXPECTED_CONTRACT_SHA256 = (
    "aa3ce118502c362ad3ef6f232e350b2d7a0e0ecae4dc270f25ee678bd61b7557")


class Stage0Error(RuntimeError):
    """The frozen Stage 0 screen cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise Stage0Error(message)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("ascii"))


def self_hash(value: Mapping[str, Any], field: str) -> str:
    body = dict(value)
    body.pop(field, None)
    return sha256_json(body)


def _contract_body() -> Mapping[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "bead": "wirehair-sxvz.16.1.20.38.2",
        "scope": "fixed off-diagonal MIX2 packet-attempt Stage 0 only",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "candidate": {
            "base_profile_sha256": EXPECTED_CANDIDATE_PROFILE_SHA256,
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "mix_pair": MIX_PAIR,
            "pair_semantics": [0, 1],
            "binary_dense_rows": BINARY_DENSE_ROWS,
            "gf256_heavy_rows": GF256_HEAVY_ROWS,
            "heavy_family": "periodic",
            "construction_seed_basis": "production-profile",
            "precode_attempt_rule": "exact frozen v9 p(K)",
            "packet_attempt_rule": "q=(p+delta) mod 256",
            "packet_seed_rule": (
                "effective_q=(effective_p+(((p+delta) mod 256)-p)*"
                "0x9e3779b9) mod 2^32"),
            "packet_attempt_stride": "0x9e3779b9",
            "global_packet_offset_domain": list(DELTAS),
            "profile_representation": (
                "one serialized p byte plus one profile-ID-bound delta; "
                "no second map or descriptor field"),
        },
        "domain": {
            "coordinates": [
                coordinate.description(index)
                for index, coordinate in enumerate(COORDINATES)
            ],
            "consumed_root_roster": list(CONSUMED_ROOTS),
            "unconsumed_final_roots_excluded": list(UNCONSUMED_FINAL_ROOTS),
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "trials": 1,
            "threads": 1,
            "overhead_stream": "paired",
            "full_payload_solve": True,
            "solve_semantics": (
                "deterministic two-byte-wide all-zero-RHS rank/work solve; "
                "the separate payload-e2e byte reconstruction is not run"),
            "cold_solve_wide_xor": "policy",
            "repetitions": list(REPETITIONS),
            "invocation_count": EXPECTED_INVOCATION_COUNT,
            "current_build_diagonal_control_observations": 6,
        },
        "transport": {
            "v9_prerequisite": {
                "path": "explicit --v9-dir",
                "contract_sha256": EXPECTED_V9_CONTRACT_SHA256,
                "source_commit": EXPECTED_V9_SOURCE_COMMIT,
                "worker_sha256": EXPECTED_V9_WORKER_SHA256,
                "input_self_sha256": EXPECTED_V9_INPUT_SELF_SHA256,
                "summary_self_sha256": EXPECTED_V9_SUMMARY_SELF_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in V9_FILES.items()
                },
                "authentication_order": (
                    "pin and semantically validate before benchmark preflight "
                    "or candidate work"),
            },
            "benchmark_pin": "open regular executable descriptor and rehash",
            "benchmark_preflight": {
                "command": "--describe",
                "schema": BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": "must equal clean source HEAD",
                "ordering": "must complete before any precodefail workload",
            },
            "source_receipt": {
                "clean_at_HEAD": True,
                "working_source_bytes_equal_HEAD_blobs": True,
                "paths": list(SOURCE_PATHS),
            },
            "child_environment": dict(CHILD_ENVIRONMENT),
            "process_timeout_seconds": PROCESS_TIMEOUT_SECONDS,
            "publication": "canonical JSONL and summary, atomic no-replace dir",
        },
        "gate": {
            "any_configuration_failure_invalidates_screen": True,
            "diagonal_control": {
                "offset": 0,
                "packet_attempt_rule": "q=p",
                "all_three_historical_weak_coordinates_rank_fail_twice": True,
                "errors": 0,
                "repeated_non_timing_fields_exact": True,
                "repeated_block_xors_equal": True,
                "repeated_gf256_muladds_equal": True,
                "authenticated_v9_fields_equal": list(
                    DIAGONAL_CONTROL_FIELDS),
                "execution_order": (
                    "six controls first; any drift stops before delta 1"),
                "rule": (
                    "candidate offsets are ineligible unless the current "
                    "pinned build reproduces all six diagonal weaknesses"),
            },
            "per_delta": {
                "all_three_coordinates_succeed_in_both_repetitions": True,
                "errors": 0,
                "configuration_failures": 0,
                "repeated_non_timing_fields_exact": True,
                "repeated_block_xors_equal": True,
                "repeated_gf256_muladds_equal": True,
            },
            "selection_order": list(DELTAS),
            "survivors": (
                "ascending subset of deltas 1..255 passing every per-delta "
                "gate"),
            "stage1_candidate_delta": (
                "lowest numeric survivor, or null; only this fixed delta is "
                "licensed for the 1080-cell consumed-root Stage 1 replay"),
            "overall": "PASS iff survivors is nonempty",
        },
        "stop_rule": (
            "no timing, broad consumed-root replay, all-K campaign, fresh-root "
            "screen, or promotion is licensed by this controller"),
    }


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("MIX2 packet-offset Stage 0 contract differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _validate_constants() -> None:
    if (MIX_PAIR != "01" or DELTAS != tuple(range(1, 256)) or
            OFFSETS != tuple(range(256)) or
            PACKET_ATTEMPT_STRIDE != 0x9e3779b9 or
            REPETITIONS != (0, 1)):
        fail("pair, delta, or repetition roster changed")
    expected_coordinates = (
        (2, 9, 10, "burst", "0x7ccd510f122fc160"),
        (1000, 0, 11, "burst", "0x7001a960b7d9c0a4"),
        (10000, 0, 9, "adversarial", "0xc4695292d9835286"),
    )
    realized = tuple(
        (value.K, value.attempt, value.root_index,
         value.schedule, value.root) for value in COORDINATES)
    if realized != expected_coordinates:
        fail("Stage 0 coordinate roster changed")
    if (len(CONSUMED_ROOTS) != 12 or
            len(set(CONSUMED_ROOTS)) != len(CONSUMED_ROOTS) or
            any(re.fullmatch(r"0x[0-9a-f]{16}", root) is None
                for root in CONSUMED_ROOTS)):
        fail("consumed root roster is malformed")
    if set(CONSUMED_ROOTS) & set(UNCONSUMED_FINAL_ROOTS):
        fail("Stage 0 root roster overlaps unconsumed final roots")
    if CONSUMED_ROOTS != tuple(screen.SELECTION_ROOTS) + tuple(screen.ROOTS):
        fail("Stage 0 consumed roots differ from the sealed v9 provenance")
    if UNCONSUMED_FINAL_ROOTS != tuple(screen.FINAL_VALIDATION_ROOTS):
        fail("Stage 0 final-root exclusion differs from the v9 holdout")
    if screen.candidate_profile_sha256() != \
            EXPECTED_CANDIDATE_PROFILE_SHA256:
        fail("Stage 0 base candidate profile identity changed")
    if any(coordinate.root not in CONSUMED_ROOTS for coordinate in COORDINATES):
        fail("Stage 0 names a root outside the consumed roster")
    if EXPECTED_INVOCATION_COUNT != 1536:
        fail("Stage 0 invocation cardinality changed")
    if set(STAIRCASE_BY_K) != {value.K for value in COORDINATES}:
        fail("Stage 0 staircase roster changed")
    if tuple(CSV_HEADER) != tuple(screen.CSV_HEADER):
        fail("precodefail CSV protocol changed")


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    delta: int
    coordinate_ordinal: int
    repetition: int

    def __post_init__(self) -> None:
        if (type(self.ordinal) is not int or
                not 0 <= self.ordinal < EXPECTED_INVOCATION_COUNT):
            fail("invocation ordinal is outside Stage 0")
        if type(self.delta) is not int or self.delta not in OFFSETS:
            fail("invocation packet-attempt delta is outside Stage 0")
        if type(self.coordinate_ordinal) is not int or not \
                0 <= self.coordinate_ordinal < len(COORDINATES):
            fail("invocation coordinate is outside Stage 0")
        if type(self.repetition) is not int or self.repetition not in REPETITIONS:
            fail("invocation repetition is outside Stage 0")
        expected_ordinal = (
            (OFFSETS.index(self.delta) * len(COORDINATES) +
             self.coordinate_ordinal) * len(REPETITIONS) + self.repetition)
        if self.ordinal != expected_ordinal:
            fail("invocation ordinal differs from the canonical Stage 0 order")

    @property
    def coordinate(self) -> Coordinate:
        return COORDINATES[self.coordinate_ordinal]

    @property
    def packet_attempt(self) -> int:
        return (self.coordinate.attempt + self.delta) & 255

    def argv(self, bench: Path) -> List[str]:
        coordinate = self.coordinate
        return [
            str(bench), "precodefail", "--N", str(coordinate.K),
            "--bb-list", "2", "--overhead", "0", "--trials", "1",
            "--threads", "1", "--loss", "0.5", "--seed",
            coordinate.root, "--schedule", coordinate.schedule,
            "--heavy-family", "periodic", "--mix-count", "2",
            "--mix-pair", MIX_PAIR, "--binary-dense-rows", "12",
            "--gf256-heavy-rows", "12", "--dense-anchors", "two07",
            "--paired-overhead-stream", "--full-payload-solve",
            "--cold-solve-wide-xor", "policy", "--exact-precode-attempt",
            str(coordinate.attempt), "--exact-packet-attempt",
            str(self.packet_attempt), "--construction-seed-basis",
            "production-profile",
        ]

    def identity(self) -> Mapping[str, Any]:
        coordinate = self.coordinate
        return {
            "ordinal": self.ordinal,
            "mix_pair": MIX_PAIR,
            "delta": self.delta,
            "precode_attempt": coordinate.attempt,
            "packet_attempt": self.packet_attempt,
            "coordinate_ordinal": self.coordinate_ordinal,
            "repetition": self.repetition,
            "K": coordinate.K,
            "root_index": coordinate.root_index,
            "loss_root": coordinate.root,
            "schedule": coordinate.schedule,
        }


def make_invocations() -> Tuple[Invocation, ...]:
    invocations: List[Invocation] = []
    for delta in OFFSETS:
        for coordinate_ordinal in range(len(COORDINATES)):
            for repetition in REPETITIONS:
                invocations.append(Invocation(
                    len(invocations), delta, coordinate_ordinal, repetition))
    if len(invocations) != EXPECTED_INVOCATION_COUNT:
        fail("Stage 0 invocation generation changed cardinality")
    identities = {
        (item.delta, item.coordinate_ordinal, item.repetition)
        for item in invocations
    }
    if len(identities) != EXPECTED_INVOCATION_COUNT:
        fail("Stage 0 invocation generation duplicated a cell")
    return tuple(invocations)


@dataclass(frozen=True)
class ProcessResult:
    invocation: Invocation
    command_sha256: str
    stdout_sha256: str
    stderr_sha256: str
    returncode: int
    stdout: bytes
    stderr: bytes


def _fd_identity(value: os.stat_result) -> Tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)


def _hash_descriptor(descriptor: int, context: str) -> Tuple[str, int]:
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} is not a regular file".format(context))
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            chunk = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset), offset)
            if not chunk:
                fail("{} was truncated while hashing".format(context))
            digest.update(chunk)
            offset += len(chunk)
        if os.pread(descriptor, 1, offset):
            fail("{} grew while hashing".format(context))
        after = os.fstat(descriptor)
        if _fd_identity(before) != _fd_identity(after):
            fail("{} changed while hashing".format(context))
        return digest.hexdigest(), before.st_size
    except OSError as exc:
        fail("cannot hash {}: {}".format(context, exc))
    raise AssertionError("unreachable")


@dataclass
class PinnedHistoricalFile:
    path: Path
    descriptor: int
    sha256: str
    size: int
    label: str

    def receipt(self) -> Mapping[str, Any]:
        if self.descriptor < 0:
            fail("{} is already closed".format(self.label))
        digest, size = _hash_descriptor(self.descriptor, self.label)
        if digest != self.sha256 or size != self.size:
            fail("{} changed after it was pinned".format(self.label))
        return {"path": str(self.path), "sha256": digest, "size": size}

    def read_bytes(self, maximum: int) -> bytes:
        if self.size > maximum:
            fail("{} exceeds its read bound".format(self.label))
        chunks: List[bytes] = []
        offset = 0
        while offset < self.size:
            try:
                chunk = os.pread(
                    self.descriptor,
                    min(1024 * 1024, self.size - offset), offset)
            except OSError as exc:
                fail("cannot read {}: {}".format(self.label, exc))
            if not chunk:
                fail("{} was truncated while reading".format(self.label))
            chunks.append(chunk)
            offset += len(chunk)
        value = b"".join(chunks)
        if self.receipt()["size"] != len(value):
            fail("{} read length changed".format(self.label))
        return value

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __enter__(self) -> "PinnedHistoricalFile":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def _open_historical_file(path: Path, expected_size: int,
                          expected_sha256: str,
                          label: str) -> PinnedHistoricalFile:
    descriptor = -1
    try:
        original = path.lstat()
        if stat.S_ISLNK(original.st_mode) or not stat.S_ISREG(original.st_mode):
            fail("{} must be a non-symlink regular file".format(label))
        resolved = path.resolve(strict=True)
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed".format(label))
        descriptor = os.open(
            str(resolved), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow)
        opened = os.fstat(descriptor)
        current = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISREG(opened.st_mode) or
                (original.st_dev, original.st_ino) !=
                    (opened.st_dev, opened.st_ino) or
                (current.st_dev, current.st_ino) !=
                    (opened.st_dev, opened.st_ino)):
            fail("{} is not one stable regular file".format(label))
        digest, size = _hash_descriptor(descriptor, label)
        if digest != expected_sha256 or size != expected_size:
            fail("{} differs from its frozen size or SHA-256".format(label))
        result = PinnedHistoricalFile(
            resolved, descriptor, digest, size, label)
        descriptor = -1
        return result
    except OSError as exc:
        fail("cannot inspect {}: {}".format(label, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    raise AssertionError("unreachable")


def _canonical_json_object(data: bytes, context: str) -> Mapping[str, Any]:
    if not data.endswith(b"\n") or b"\r" in data:
        fail("{} framing changed".format(context))
    try:
        value = json.loads(data.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail("{} is not canonical JSON: {}".format(context, exc))
    if type(value) is not dict or data != \
            (canonical_json(value) + "\n").encode("ascii"):
        fail("{} is not one canonical JSON object".format(context))
    return value


def _self_hash_matches(value: Mapping[str, Any], field: str,
                       expected: str, context: str) -> None:
    if value.get(field) != expected or self_hash(value, field) != expected:
        fail("{} self-hash changed".format(context))


@dataclass
class V9Evidence:
    stack: ExitStack
    pins: Mapping[str, PinnedHistoricalFile]
    receipt_value: Mapping[str, Any]
    expected_precode_seeds: Mapping[int, str]
    expected_packet_seeds: Mapping[int, str]
    expected_diagonal_control: Mapping[int, Mapping[str, Any]]

    def receipt(self) -> Mapping[str, Any]:
        files = {name: pin.receipt() for name, pin in self.pins.items()}
        if files != self.receipt_value["files"]:
            fail("V9 pinned files changed")
        return self.receipt_value

    def close(self) -> None:
        self.stack.close()

    def __enter__(self) -> "V9Evidence":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def load_v9_evidence(directory: Path) -> V9Evidence:
    stack = ExitStack()
    try:
        original = directory.lstat()
        if stat.S_ISLNK(original.st_mode) or not stat.S_ISDIR(original.st_mode):
            fail("V9 directory must be a non-symlink directory")
        root = directory.resolve(strict=True)
        pins = {
            name: stack.enter_context(_open_historical_file(
                root / name, size, digest, "V9 {}".format(name)))
            for name, (size, digest) in V9_FILES.items()
        }

        attempts_data = pins[
            "promotion-short-screen-attempts.jsonl"].read_bytes(1024 * 1024)
        derivations = screen.parse_derivation_stream(
            attempts_data, EXPECTED_V9_WORKER_SHA256)
        if (tuple(record["K"] for record in derivations) != V9_K_VALUES or
                tuple(record["selected_attempt"] for record in derivations) !=
                    V9_ATTEMPTS):
            fail("V9 attempt stream differs from the frozen p(K) map")

        input_value = _canonical_json_object(
            pins["promotion-short-screen-input.json"].read_bytes(
                2 * 1024 * 1024), "V9 input")
        summary = _canonical_json_object(
            pins["promotion-short-screen-summary.json"].read_bytes(
                1024 * 1024), "V9 summary")
        _self_hash_matches(
            input_value, "input_sha256", EXPECTED_V9_INPUT_SELF_SHA256,
            "V9 input")
        _self_hash_matches(
            summary, "summary_sha256", EXPECTED_V9_SUMMARY_SELF_SHA256,
            "V9 summary")
        if (input_value.get("contract_sha256") !=
                EXPECTED_V9_CONTRACT_SHA256 or
                input_value.get("attempt_selection_stream_sha256") !=
                    V9_FILES["promotion-short-screen-attempts.jsonl"][1] or
                input_value.get("selected_attempts") != list(V9_ATTEMPTS) or
                summary.get("contract_sha256") !=
                    EXPECTED_V9_CONTRACT_SHA256 or
                summary.get("source_git_commit") !=
                    EXPECTED_V9_SOURCE_COMMIT or
                summary.get("repair_worker_binary_sha256") !=
                    EXPECTED_V9_WORKER_SHA256 or
                summary.get("candidate_profile_sha256") !=
                    EXPECTED_CANDIDATE_PROFILE_SHA256 or
                summary.get("input_sha256") !=
                    EXPECTED_V9_INPUT_SELF_SHA256 or
                summary.get("input_file_sha256") !=
                    V9_FILES["promotion-short-screen-input.json"][1] or
                summary.get("attempt_selection_stream_sha256") !=
                    V9_FILES["promotion-short-screen-attempts.jsonl"][1] or
                summary.get("result_stream_sha256") !=
                    V9_FILES["promotion-short-screen-results.jsonl"][1] or
                summary.get("record_count") != 1080 or
                summary.get("candidate_weak_coordinates") != 3 or
                summary.get("candidate_weak_observations") != 6 or
                summary.get("disposition") != "REJECT"):
            fail("V9 input or summary provenance changed")

        results_data = pins[
            "promotion-short-screen-results.jsonl"].read_bytes(2 * 1024 * 1024)
        weak: Dict[Tuple[int, str, str, int], Mapping[str, Any]] = {}
        lines = results_data.splitlines(keepends=True)
        if len(lines) != 1080:
            fail("V9 result stream has the wrong cardinality")
        for ordinal, line in enumerate(lines):
            if not line.endswith(b"\n") or b"\r" in line:
                fail("V9 result framing changed")
            try:
                row = json.loads(line.decode("ascii"))
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                fail("V9 result {} is malformed: {}".format(ordinal, exc))
            if (type(row) is not dict or row.get("ordinal") != ordinal or
                    line != (canonical_json(row) + "\n").encode("ascii")):
                fail("V9 result {} is not canonical".format(ordinal))
            if row.get("arm") != "candidate_two07_mix2" or \
                    row.get("weak") is not True:
                continue
            key = (row.get("K"), row.get("loss_root"), row.get("schedule"),
                   row.get("observation_index"))
            if (key in weak or row.get("attempt_mode") != "exact" or
                    row.get("construction_attempt") !=
                        dict(zip(V9_K_VALUES, V9_ATTEMPTS)).get(row.get("K")) or
                    row.get("mix_count") != 2 or
                    row.get("dense_anchor_layout") != "two07" or
                    row.get("candidate_profile_sha256") !=
                        EXPECTED_CANDIDATE_PROFILE_SHA256 or
                    row.get("success") != 0 or row.get("rank_fail") != 1 or
                    row.get("error") != 0 or
                    HEX64.fullmatch(row.get("effective_precode_seed", "")) is
                        None):
                fail("V9 weak result differs from its frozen protocol")
            weak[key] = row
        expected_weak = {
            (coordinate.K, coordinate.root, coordinate.schedule, repetition)
            for coordinate in COORDINATES for repetition in REPETITIONS
        }
        if set(weak) != expected_weak:
            fail("V9 weak coordinates differ from the frozen Stage 0 domain")
        precode_seeds = {}
        packet_seeds = {}
        diagonal_control = {}
        for coordinate_ordinal, coordinate in enumerate(COORDINATES):
            seeds = {
                weak[(coordinate.K, coordinate.root, coordinate.schedule,
                      repetition)]["effective_precode_seed"]
                for repetition in REPETITIONS
            }
            if len(seeds) != 1:
                fail("V9 weak coordinate changed effective precode seed")
            precode_seeds[coordinate_ordinal] = next(iter(seeds))
            packet_values = {
                weak[(coordinate.K, coordinate.root, coordinate.schedule,
                      repetition)].get("effective_packet_seed")
                for repetition in REPETITIONS
            }
            if (len(packet_values) != 1 or
                    HEX32.fullmatch(next(iter(packet_values)) or "") is None):
                fail("V9 weak coordinate changed effective packet seed")
            packet_seeds[coordinate_ordinal] = next(iter(packet_values))
            control_values = []
            for repetition in REPETITIONS:
                row = weak[(coordinate.K, coordinate.root,
                            coordinate.schedule, repetition)]
                value = {
                    "staircase": row.get("staircase"),
                    "source_hits": row.get("source_hits"),
                    "outcome": "rank_fail",
                    "success": False,
                    "rank_fail": True,
                    "error": False,
                    "inactivated_columns": row.get("inactivated_columns"),
                    "effective_precode_seed":
                        row.get("effective_precode_seed"),
                    "effective_packet_seed":
                        row.get("effective_packet_seed"),
                    "block_xors": row.get("block_xors"),
                    "gf256_muladds": row.get("gf256_muladds"),
                }
                if (type(value["staircase"]) is not int or
                        value["staircase"] != STAIRCASE_BY_K[coordinate.K] or
                        type(value["source_hits"]) is not int or
                        value["source_hits"] !=
                            (3 if coordinate.K >= 10000 else 2) or
                        type(value["inactivated_columns"]) is not int or
                        value["inactivated_columns"] < 0 or
                        type(value["block_xors"]) is not int or
                        value["block_xors"] < 0 or
                        type(value["gf256_muladds"]) is not int or
                        value["gf256_muladds"] < 0):
                    fail("V9 diagonal control counters are malformed")
                control_values.append(value)
            if control_values[0] != control_values[1]:
                fail("V9 diagonal control repetitions differ")
            diagonal_control[coordinate_ordinal] = control_values[0]

        file_receipts = {name: pin.receipt() for name, pin in pins.items()}
        receipt: Dict[str, Any] = {
            "contract_sha256": EXPECTED_V9_CONTRACT_SHA256,
            "source_git_commit": EXPECTED_V9_SOURCE_COMMIT,
            "worker_sha256": EXPECTED_V9_WORKER_SHA256,
            "input_self_sha256": EXPECTED_V9_INPUT_SELF_SHA256,
            "summary_self_sha256": EXPECTED_V9_SUMMARY_SELF_SHA256,
            "files": file_receipts,
            "selected_attempts": list(V9_ATTEMPTS),
            "weak_coordinates": [
                coordinate.description(index)
                for index, coordinate in enumerate(COORDINATES)
            ],
            "effective_precode_seeds": [
                precode_seeds[index] for index in range(len(COORDINATES))
            ],
            "effective_packet_seeds": [
                packet_seeds[index] for index in range(len(COORDINATES))
            ],
            "diagonal_control": [
                diagonal_control[index]
                for index in range(len(COORDINATES))
            ],
        }
        receipt["receipt_sha256"] = self_hash(receipt, "receipt_sha256")
        evidence = V9Evidence(
            stack, pins, receipt, precode_seeds, packet_seeds,
            diagonal_control)
        evidence.receipt()
        return evidence
    except BaseException:
        stack.close()
        raise


def _parse_uint(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or UINT.fullmatch(text) is None:
        fail("{} is not canonical unsigned decimal".format(context))
    if len(text) > len(str(maximum)):
        fail("{} exceeds its bound".format(context))
    value = int(text)
    if value > maximum:
        fail("{} exceeds its bound".format(context))
    return value


def _packet_seed_for_offset(effective_p: str, p: int, delta: int) -> str:
    if (type(effective_p) is not str or
            HEX32.fullmatch(effective_p) is None or
            type(p) is not int or not 0 <= p <= 255 or
            type(delta) is not int or delta not in OFFSETS):
        fail("packet-offset seed input is malformed")
    q = (p + delta) & 255
    value = (int(effective_p, 16) +
             (q - p) * PACKET_ATTEMPT_STRIDE) & 0xffffffff
    return "0x{:08x}".format(value)


def _parse_integral_decimal(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical three-place decimal".format(context))
    whole, fraction = text.split(".")
    if fraction != "000":
        fail("{} is not integral for a one-trial cell".format(context))
    return _parse_uint(whole, maximum, context)


def _parse_timing_decimal(text: str, context: str) -> None:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical nonnegative timing".format(context))
    _parse_uint(text.split(".")[0], UINT64_MAX, context)


def _parse_histogram(text: str, expected_key: int, context: str) -> None:
    if text != "{}:1".format(expected_key):
        fail("{} differs from its one-trial counter".format(context))


def _parse_metadata(line: str, invocation: Invocation,
                    source_commit: str) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("precodefail metadata is missing")
    tokens = line[len(prefix):].split(" ")
    values: Dict[str, str] = {}
    order: List[str] = []
    for token in tokens:
        if token.count("=") != 1:
            fail("precodefail metadata token is malformed")
        key, value = token.split("=", 1)
        if not key or not value or key in values:
            fail("precodefail metadata is ambiguous")
        order.append(key)
        values[key] = value
    if tuple(order) != METADATA_ORDER:
        fail("precodefail metadata fields or canonical order changed")
    coordinate = invocation.coordinate
    expected = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": coordinate.root, "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "1", "schedule": coordinate.schedule,
        "cold_solve_wide_xor": "policy", "exact_attempt_mode": "1",
        "exact_precode_attempt": str(coordinate.attempt),
        "exact_packet_attempt": str(invocation.packet_attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
        "source_git_commit": source_commit, "mix_pair": MIX_PAIR,
    }
    if values != expected:
        fail("precodefail metadata differs from the frozen invocation")
    return values


def _parse_csv_line(line: str, expected_width: int) -> List[str]:
    try:
        rows = list(csv.reader([line], strict=True))
    except csv.Error as exc:
        fail("precodefail CSV is malformed: {}".format(exc))
    if len(rows) != 1 or len(rows[0]) != expected_width:
        fail("precodefail CSV has the wrong field count")
    if line != ",".join(rows[0]):
        fail("precodefail CSV is not canonical unquoted CSV")
    return rows[0]


def _framed_lines(result: ProcessResult,
                  source_commit: str) -> Tuple[str, ...]:
    expected_command = sha256_json({
        "argv": result.invocation.argv(Path("wirehair_v2_bench"))[1:],
    })
    if (SHA256.fullmatch(result.command_sha256) is None or
            result.command_sha256 != expected_command or
            SHA256.fullmatch(result.stdout_sha256) is None or
            result.stdout_sha256 != sha256_bytes(result.stdout) or
            SHA256.fullmatch(result.stderr_sha256) is None or
            result.stderr_sha256 != sha256_bytes(result.stderr)):
        fail("precodefail process receipt is inconsistent")
    if (len(result.stdout) > MAX_STDOUT_BYTES or
            len(result.stderr) > MAX_STDERR_BYTES):
        fail("precodefail output exceeds its bound")
    if not result.stdout.endswith(b"\n") or b"\r" in result.stdout:
        fail("precodefail stdout framing changed")
    try:
        text = result.stdout.decode("ascii")
    except UnicodeDecodeError:
        fail("precodefail stdout is not ASCII")
    lines = tuple(text.splitlines())
    if len(lines) < 2 or not lines[0] or not lines[1]:
        fail("precodefail omitted metadata or CSV header")
    _parse_metadata(lines[0], result.invocation, source_commit)
    if tuple(_parse_csv_line(lines[1], len(CSV_HEADER))) != CSV_HEADER:
        fail("precodefail CSV header changed")
    return lines


def _empty_terminal() -> Dict[str, Any]:
    return {
        "staircase": None,
        "source_hits": None,
        "outcome": "configuration_failure",
        "success": False,
        "rank_fail": False,
        "error": False,
        "configuration_failure": True,
        "weak": True,
        "inactivated_columns": None,
        "binary_deficiency": None,
        "gf256_rank_gain": None,
        "heavy_shortfall": None,
        "first_rank_fail": None,
        "failure_trials": None,
        "effective_precode_seed": None,
        "effective_packet_seed": None,
        "block_xors": None,
        "gf256_muladds": None,
        "configuration_stderr_sha256": None,
    }


def _parse_terminal_row(row: Mapping[str, str], invocation: Invocation) \
        -> Mapping[str, Any]:
    coordinate = invocation.coordinate
    expected = {
        "N": str(coordinate.K), "bb": "2", "heavy_family": "periodic",
        "mix_count": "2", "binary_dense_rows": "12",
        "gf256_heavy_rows": "12", "dense_identity_corner": "0",
        "overhead": "0", "trials": "1", "seed_attempt": "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(coordinate.attempt),
        "packet_attempt": str(invocation.packet_attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
    }
    for field, value in expected.items():
        if row[field] != value:
            fail("precodefail field {} differs from Stage 0".format(field))
    staircase = _parse_uint(row["staircase"], 4096, "staircase")
    source_hits = _parse_uint(row["source_hits"], 8, "source hits")
    if staircase != STAIRCASE_BY_K[coordinate.K] or \
            source_hits != (3 if coordinate.K >= 10000 else 2):
        fail("precodefail realized profile geometry changed")

    success = _parse_uint(row["success"], 1, "success")
    rank_fail = _parse_uint(row["rank_fail"], 1, "rank failure")
    error = _parse_uint(row["error"], 1, "error")
    if success + rank_fail + error != 1:
        fail("precodefail terminal counts do not sum to one")
    weak = success != 1
    expected_rate = "0.00000000" if success else "1.00000000"
    if EIGHT_DECIMAL.fullmatch(row["fail_rate"]) is None or \
            row["fail_rate"] != expected_rate:
        fail("precodefail failure rate disagrees with its outcome")

    inactivated = _parse_integral_decimal(
        row["inact_mu"], 65535, "inactivation")
    if _parse_uint(row["inact_max"], 65535, "maximum inactivation") != \
            inactivated:
        fail("precodefail inactivation counters disagree")
    deficiency = _parse_integral_decimal(
        row["binary_def_mu"], 65535, "binary deficiency")
    gain = _parse_integral_decimal(
        row["heavy_gain_mu"], 65535, "GF256 rank gain")
    if (_parse_uint(row["binary_def_max"], 65535,
                    "maximum binary deficiency") != deficiency or
            _parse_uint(row["heavy_gain_min"], 65535,
                        "minimum GF256 rank gain") != gain or
            deficiency > inactivated or gain > deficiency or gain > 12):
        fail("precodefail rank counters are inconsistent")
    if success and gain != deficiency:
        fail("precodefail success is rank-deficient")
    if rank_fail and gain >= deficiency:
        fail("precodefail rank failure has full rank")
    shortfall = _parse_uint(row["heavy_shortfall"], 1, "heavy shortfall")
    expected_shortfall = int(
        bool(rank_fail) and deficiency <= 12 and gain < deficiency)
    if shortfall != expected_shortfall:
        fail("precodefail heavy shortfall is inconsistent")
    _parse_histogram(row["binary_def_hist"], deficiency,
                     "binary deficiency histogram")
    _parse_histogram(row["heavy_gain_hist"], gain,
                     "GF256 rank gain histogram")
    expected_first = "0" if rank_fail else "-1"
    if row["first_rank_fail"] != expected_first or \
            row["failure_trials"] != ("0" if weak else ""):
        fail("precodefail failure diagnostics are inconsistent")
    for field in (
            "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
            "residual_ms_mu", "backsub_ms_mu"):
        _parse_timing_decimal(row[field], field)
    if HEX64.fullmatch(row["effective_precode_seed"]) is None or \
            HEX32.fullmatch(row["effective_packet_seed"]) is None:
        fail("precodefail effective seed is malformed")

    outcome = "success" if success else ("rank_fail" if rank_fail else "error")
    return {
        "staircase": staircase,
        "source_hits": source_hits,
        "outcome": outcome,
        "success": bool(success),
        "rank_fail": bool(rank_fail),
        "error": bool(error),
        "configuration_failure": False,
        "weak": weak,
        "inactivated_columns": inactivated,
        "binary_deficiency": deficiency,
        "gf256_rank_gain": gain,
        "heavy_shortfall": shortfall,
        "first_rank_fail": 0 if rank_fail else -1,
        "failure_trials": [0] if weak else [],
        "effective_precode_seed": row["effective_precode_seed"],
        "effective_packet_seed": row["effective_packet_seed"],
        "block_xors": _parse_integral_decimal(
            row["block_xors_mu"], UINT64_MAX, "block XORs"),
        "gf256_muladds": _parse_integral_decimal(
            row["block_muladds_mu"], UINT64_MAX, "GF256 muladds"),
        "configuration_stderr_sha256": None,
    }


def _base_record(result: ProcessResult,
                 source_commit: str,
                 bench_sha256: str,
                 source_receipt_sha256: str) -> Dict[str, Any]:
    identity = result.invocation.identity()
    return {
        "schema": RECORD_SCHEMA,
        **identity,
        "block_bytes": BLOCK_BYTES,
        "loss_ppm": LOSS_PPM,
        "overhead": OVERHEAD,
        "dense_anchor_layout": "two07",
        "mix_count": 2,
        "binary_dense_rows": BINARY_DENSE_ROWS,
        "gf256_heavy_rows": GF256_HEAVY_ROWS,
        "heavy_family": "periodic",
        "construction_seed_basis": "production-profile",
        "full_payload_solve": True,
        "overhead_stream": "paired",
        "cold_solve_wide_xor": "policy",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "timing_evidence": False,
        "command_sha256": result.command_sha256,
        "stdout_sha256": result.stdout_sha256,
        "stderr_sha256": result.stderr_sha256,
        "returncode": result.returncode,
        "bench_binary_sha256": bench_sha256,
        "source_git_commit": source_commit,
        "source_receipt_sha256": source_receipt_sha256,
    }


def parse_process_result(result: ProcessResult,
                         source_commit: str,
                         bench_sha256: str,
                         source_receipt_sha256: str) -> Mapping[str, Any]:
    if COMMIT.fullmatch(source_commit) is None or \
            SHA256.fullmatch(bench_sha256) is None or \
            SHA256.fullmatch(source_receipt_sha256) is None:
        fail("Stage 0 receipt identity is malformed")
    lines = _framed_lines(result, source_commit)
    record = _base_record(
        result, source_commit, bench_sha256, source_receipt_sha256)
    if result.returncode == 0:
        if result.stderr or len(lines) != 3:
            fail("successful precodefail process framing changed")
        values = _parse_csv_line(lines[2], len(CSV_HEADER))
        terminal = _parse_terminal_row(dict(zip(CSV_HEADER, values)),
                                       result.invocation)
    elif result.returncode == 2:
        coordinate = result.invocation.coordinate
        expected_stderr = (
            "precodefail configuration failed N={} bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt={} packet_attempt={} result=2\n".format(
                coordinate.K, coordinate.attempt,
                result.invocation.packet_attempt)).encode("ascii")
        if len(lines) != 2 or result.stderr != expected_stderr:
            fail("precodefail configuration failure is noncanonical")
        terminal = _empty_terminal()
        terminal["configuration_stderr_sha256"] = result.stderr_sha256
    else:
        fail("precodefail process exited with an inadmissible status")
    record.update(terminal)
    record["deterministic_projection_sha256"] = sha256_json(
        deterministic_projection(record))
    return record


DETERMINISTIC_FIELDS = (
    "mix_pair", "delta", "precode_attempt", "packet_attempt",
    "coordinate_ordinal", "K", "root_index",
    "loss_root", "schedule", "block_bytes", "loss_ppm", "overhead",
    "dense_anchor_layout", "mix_count", "binary_dense_rows",
    "gf256_heavy_rows", "heavy_family", "construction_seed_basis",
    "full_payload_solve", "overhead_stream", "cold_solve_wide_xor",
    "returncode", "staircase", "source_hits", "outcome", "success",
    "rank_fail", "error", "configuration_failure", "weak",
    "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
    "heavy_shortfall", "first_rank_fail", "failure_trials",
    "effective_precode_seed", "effective_packet_seed", "block_xors",
    "gf256_muladds", "configuration_stderr_sha256", "command_sha256",
    "bench_binary_sha256", "source_git_commit", "source_receipt_sha256",
)


def deterministic_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in DETERMINISTIC_FIELDS if field not in record]
    if missing:
        fail("Stage 0 record omits deterministic fields: {}".format(missing))
    return {field: record[field] for field in DETERMINISTIC_FIELDS}


def _run_raw(invocation: Invocation,
             pinned: screen.PinnedExecutable) -> ProcessResult:
    argv = invocation.argv(pinned.path)
    try:
        completed = subprocess.run(
            argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=PROCESS_TIMEOUT_SECONDS,
            check=False, shell=False, close_fds=True,
            start_new_session=True, env=CHILD_ENVIRONMENT,
            executable="/proc/self/fd/{}".format(pinned.descriptor),
            pass_fds=(pinned.descriptor,))
    except (OSError, subprocess.SubprocessError) as exc:
        fail("precodefail workload failed to execute: {}".format(exc))
    return ProcessResult(
        invocation=invocation,
        command_sha256=sha256_json({"argv": argv[1:]}),
        stdout_sha256=sha256_bytes(completed.stdout),
        stderr_sha256=sha256_bytes(completed.stderr),
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )


def _execute_roster(pinned: screen.PinnedExecutable,
                    source_commit: str,
                    bench_sha256: str,
                    source_receipt_sha256: str,
                    expected_diagonal_control:
                        Mapping[int, Mapping[str, Any]]) \
        -> Tuple[Mapping[str, Any], ...]:
    records = []
    for invocation in make_invocations():
        result = _run_raw(invocation, pinned)
        records.append(parse_process_result(
            result, source_commit, bench_sha256, source_receipt_sha256))
        if len(records) == len(COORDINATES) * len(REPETITIONS):
            _validate_diagonal_control(records, expected_diagonal_control)
    if len(records) != EXPECTED_INVOCATION_COUNT:
        fail("Stage 0 executed the wrong number of invocations")
    return tuple(records)


def _validate_diagonal_control(
        records: Sequence[Mapping[str, Any]],
        expected: Mapping[int, Mapping[str, Any]]) -> Mapping[str, Any]:
    control_count = len(COORDINATES) * len(REPETITIONS)
    if (len(records) != control_count or
            set(expected) != set(range(len(COORDINATES)))):
        fail("diagonal control roster or historical projection is malformed")
    exact_coordinates = 0
    work_equal_coordinates = 0
    coordinate_results = []
    for coordinate in range(len(COORDINATES)):
        first = records[coordinate * len(REPETITIONS)]
        second = records[coordinate * len(REPETITIONS) + 1]
        for repetition, record in enumerate((first, second)):
            invocation = Invocation(
                coordinate * len(REPETITIONS) + repetition,
                0, coordinate, repetition)
            identity = invocation.identity()
            if (type(record) is not dict or
                    any(record.get(field) != value
                        for field, value in identity.items()) or
                    record.get("configuration_failure") is not False):
                fail("current diagonal control identity or construction failed")
            historical = expected[coordinate]
            if (type(historical) is not dict or
                    set(historical) != set(DIAGONAL_CONTROL_FIELDS) or
                    any(record.get(field) != historical[field]
                        for field in DIAGONAL_CONTROL_FIELDS)):
                fail("current diagonal control differs from authenticated v9")
        exact = deterministic_projection(first) == \
            deterministic_projection(second)
        work_equal = (
            type(first["block_xors"]) is int and
            type(second["block_xors"]) is int and
            type(first["gf256_muladds"]) is int and
            type(second["gf256_muladds"]) is int and
            first["block_xors"] == second["block_xors"] and
            first["gf256_muladds"] == second["gf256_muladds"])
        if not exact or not work_equal:
            fail("current diagonal control repetitions differ")
        exact_coordinates += 1
        work_equal_coordinates += 1
        coordinate_results.append({
            "coordinate_ordinal": coordinate,
            "both_rank_fail": True,
            "repeated_non_timing_fields_exact": True,
            "repeated_work_equal": True,
        })
    return {
        "offset": 0,
        "record_count": control_count,
        "successes": 0,
        "rank_failures": control_count,
        "errors": 0,
        "repeated_exact_coordinate_count": exact_coordinates,
        "repeated_work_equal_coordinate_count": work_equal_coordinates,
        "coordinates": coordinate_results,
        "disposition": "PASS",
    }


def adjudicate(records: Sequence[Mapping[str, Any]],
               expected_precode_seeds: Mapping[int, str],
               expected_packet_seeds: Mapping[int, str],
               expected_diagonal_control:
                   Mapping[int, Mapping[str, Any]]) \
        -> Mapping[str, Any]:
    if len(records) != EXPECTED_INVOCATION_COUNT:
        fail("Stage 0 result roster has the wrong cardinality")
    expected_coordinate_ordinals = set(range(len(COORDINATES)))
    if (set(expected_precode_seeds) != expected_coordinate_ordinals or
            any(type(value) is not str or HEX64.fullmatch(value) is None
                for value in expected_precode_seeds.values())):
        fail("authenticated v9 precode seed map is malformed")
    if (set(expected_packet_seeds) != expected_coordinate_ordinals or
            any(type(value) is not str or HEX32.fullmatch(value) is None
                for value in expected_packet_seeds.values())):
        fail("authenticated v9 packet seed map is malformed")
    expected = {
        (delta, coordinate, repetition)
        for delta in OFFSETS
        for coordinate in range(len(COORDINATES))
        for repetition in REPETITIONS
    }
    by_key: Dict[Tuple[int, int, int], Mapping[str, Any]] = {}
    ordinals = set()
    invocations = make_invocations()
    for expected_ordinal, record in enumerate(records):
        if type(record) is not dict:
            fail("Stage 0 record is not an object")
        expected_identity = invocations[expected_ordinal].identity()
        if any(record.get(field) != value
               for field, value in expected_identity.items()):
            fail("Stage 0 record identity differs from the frozen roster")
        try:
            key = (record["delta"], record["coordinate_ordinal"],
                   record["repetition"])
            ordinal = record["ordinal"]
        except KeyError as exc:
            fail("Stage 0 record omits identity: {}".format(exc))
        if (type(ordinal) is not int or ordinal != expected_ordinal or
                key in by_key or key not in expected or ordinal in ordinals):
            fail("Stage 0 result roster has duplicate or unknown identity")
        projection_digest = record.get("deterministic_projection_sha256")
        if (record.get("schema") != RECORD_SCHEMA or
                record.get("promotion_evidence") is not False or
                record.get("fresh_roots_used") is not False or
                record.get("timing_evidence") is not False or
                type(record.get("success")) is not bool or
                type(record.get("rank_fail")) is not bool or
                type(record.get("error")) is not bool or
                type(record.get("configuration_failure")) is not bool or
                type(projection_digest) is not str or
                SHA256.fullmatch(projection_digest) is None or
                projection_digest !=
                    sha256_json(deterministic_projection(record))):
            fail("Stage 0 record schema or deterministic receipt changed")
        coordinate = record["coordinate_ordinal"]
        delta = record["delta"]
        if record["configuration_failure"]:
            fail("fixed precode p failed before packet offset q was used")
        if record.get("effective_precode_seed") != \
                expected_precode_seeds[coordinate]:
            fail("candidate precode seed differs from authenticated v9 p")
        if record.get("effective_packet_seed") != \
                _packet_seed_for_offset(
                    expected_packet_seeds[coordinate],
                    COORDINATES[coordinate].attempt, delta):
            fail("candidate packet seed differs from the authenticated q rule")
        by_key[key] = record
        ordinals.add(ordinal)
    if set(by_key) != expected or \
            ordinals != set(range(EXPECTED_INVOCATION_COUNT)):
        fail("Stage 0 result roster is incomplete or noncanonical")

    baseline_records = [
        by_key[(0, coordinate, repetition)]
        for coordinate in range(len(COORDINATES))
        for repetition in REPETITIONS
    ]
    baseline = _validate_diagonal_control(
        baseline_records, expected_diagonal_control)

    delta_results = []
    survivors: List[int] = []
    for delta in DELTAS:
        delta_records = [
            by_key[(delta, coordinate, repetition)]
            for coordinate in range(len(COORDINATES))
            for repetition in REPETITIONS
        ]
        exact_coordinates = 0
        work_equal_coordinates = 0
        coordinate_results = []
        for coordinate in range(len(COORDINATES)):
            first = by_key[(delta, coordinate, 0)]
            second = by_key[(delta, coordinate, 1)]
            exact = deterministic_projection(first) == \
                deterministic_projection(second)
            work_values = (
                first["block_xors"], second["block_xors"],
                first["gf256_muladds"], second["gf256_muladds"])
            work_equal = (
                all(type(value) is int for value in work_values) and
                first["block_xors"] == second["block_xors"] and
                first["gf256_muladds"] == second["gf256_muladds"])
            exact_coordinates += int(exact)
            work_equal_coordinates += int(work_equal)
            coordinate_results.append({
                "coordinate_ordinal": coordinate,
                "both_success": first["success"] and second["success"],
                "repeated_non_timing_fields_exact": exact,
                "repeated_work_equal": work_equal,
            })
        successes = sum(int(record["success"]) for record in delta_records)
        rank_failures = sum(
            int(record["rank_fail"]) for record in delta_records)
        errors = sum(int(record["error"]) for record in delta_records)
        configuration_failures = sum(
            int(record["configuration_failure"]) for record in delta_records)
        passed = (
            successes == 6 and rank_failures == 0 and errors == 0 and
            configuration_failures == 0 and exact_coordinates == 3 and
            work_equal_coordinates == 3)
        if passed:
            survivors.append(delta)
        delta_results.append({
            "delta": delta,
            "record_count": len(delta_records),
            "successes": successes,
            "rank_failures": rank_failures,
            "errors": errors,
            "configuration_failures": configuration_failures,
            "repeated_exact_coordinate_count": exact_coordinates,
            "repeated_work_equal_coordinate_count": work_equal_coordinates,
            "coordinates": coordinate_results,
            "disposition": "PASS" if passed else "REJECT",
        })
    return {
        "diagonal_control": baseline,
        "delta_results": delta_results,
        "survivors": survivors,
        "stage1_candidate_delta": survivors[0] if survivors else None,
        "effective_precode_seeds": [
            expected_precode_seeds[coordinate]
            for coordinate in range(len(COORDINATES))
        ],
        "authenticated_effective_p_packet_seeds": [
            expected_packet_seeds[coordinate]
            for coordinate in range(len(COORDINATES))
        ],
        "disposition": "PASS" if survivors else "REJECT",
    }


def _source_receipt() -> Mapping[str, Any]:
    for item in SOURCE_PATHS:
        if not (REPO_ROOT / item).is_file():
            fail("required Stage 0 source is missing: {}".format(item))
    try:
        head_process = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
        status_process = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all",
             "--"] + list(SOURCE_STATUS_PATHS), cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
        tracked_process = subprocess.run(
            ["git", "ls-files", "-z", "--error-unmatch", "--"] +
            list(SOURCE_PATHS), cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot attest Stage 0 source checkout: {}".format(exc))
    if (head_process.returncode != 0 or head_process.stderr or
            status_process.returncode != 0 or status_process.stderr or
            tracked_process.returncode != 0 or tracked_process.stderr):
        fail("cannot attest Stage 0 source checkout")
    try:
        head = head_process.stdout.decode("ascii").strip()
    except UnicodeDecodeError:
        fail("Stage 0 source HEAD is not ASCII")
    if COMMIT.fullmatch(head) is None or status_process.stdout:
        fail("Stage 0 linked sources are not clean at a full source HEAD")
    try:
        tracked = tracked_process.stdout.decode("ascii").split("\0")
    except UnicodeDecodeError:
        fail("Stage 0 tracked source roster is not ASCII")
    if tracked and tracked[-1] == "":
        tracked.pop()
    if len(tracked) != len(SOURCE_PATHS) or set(tracked) != set(SOURCE_PATHS):
        fail("Stage 0 source receipt includes a path not tracked at HEAD")
    hashes = {
        item: sha256_bytes(screen._stable_file_bytes(REPO_ROOT / item))
        for item in SOURCE_PATHS
    }
    for item in SOURCE_PATHS:
        try:
            blob = subprocess.run(
                ["git", "cat-file", "blob", "{}:{}".format(head, item)],
                cwd=str(REPO_ROOT), stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False, timeout=30)
        except (OSError, subprocess.SubprocessError) as exc:
            fail("cannot read HEAD blob {}: {}".format(item, exc))
        if (blob.returncode != 0 or blob.stderr or
                sha256_bytes(blob.stdout) != hashes[item]):
            fail("Stage 0 source {} differs byte-for-byte from HEAD".format(
                item))
    try:
        final_head_process = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot re-attest Stage 0 source HEAD: {}".format(exc))
    if (final_head_process.returncode != 0 or final_head_process.stderr or
            final_head_process.stdout != (head + "\n").encode("ascii")):
        fail("Stage 0 source HEAD changed while the receipt was built")
    receipt: Dict[str, Any] = {
        "source_git_commit": head,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "working_source_bytes_equal_HEAD_blobs": True,
        "clean_status_scope": list(SOURCE_STATUS_PATHS),
        "source_files": hashes,
    }
    receipt["source_receipt_sha256"] = self_hash(
        receipt, "source_receipt_sha256")
    return receipt


def _write_jsonl(path: Path,
                 rows: Iterable[Mapping[str, Any]]) -> Tuple[str, int]:
    digest = hashlib.sha256()
    count = 0
    try:
        descriptor = os.open(
            str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except OSError as exc:
        fail("cannot create {}: {}".format(path, exc))
    try:
        for row in rows:
            data = (canonical_json(row) + "\n").encode("ascii")
            digest.update(data)
            count += 1
            offset = 0
            while offset < len(data):
                written = os.write(descriptor, data[offset:])
                if written <= 0:
                    fail("short write creating {}".format(path))
                offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return digest.hexdigest(), count


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(str(path), flags)
    except OSError as exc:
        fail("cannot open output staging directory: {}".format(exc))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def run_stage0(bench: Path, v9_dir: Path,
               output_dir: Path) -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path must be one fresh path below an existing directory")
    with load_v9_evidence(v9_dir) as v9:
        source = _source_receipt()
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            return _run_stage0_pinned(
                pinned, output_dir, source, contract, v9)


def _run_stage0_pinned(
        pinned: screen.PinnedExecutable,
        output_dir: Path,
        source: Mapping[str, Any],
        contract: Mapping[str, Any],
        v9: V9Evidence) -> Mapping[str, Any]:
    source_commit = source["source_git_commit"]
    source_receipt_sha256 = source["source_receipt_sha256"]
    bench_receipt = pinned.receipt()
    # This is deliberately the last action before any precodefail invocation.
    description = screen.read_bench_description(
        pinned.path, source_commit, pinned.descriptor)
    records = _execute_roster(
        pinned, source_commit, bench_receipt["sha256"],
        source_receipt_sha256, v9.expected_diagonal_control)
    verdict = adjudicate(
        records, v9.expected_precode_seeds, v9.expected_packet_seeds,
        v9.expected_diagonal_control)

    temporary = Path(tempfile.mkdtemp(
        prefix=".wh2-mix2-packet-offset-stage0-", dir=str(output_dir.parent)))
    published = False
    try:
        records_sha256, record_count = _write_jsonl(
            temporary / RECORD_NAME, records)
        if record_count != EXPECTED_INVOCATION_COUNT:
            fail("Stage 0 JSONL cardinality changed")
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                v9.receipt() != v9.receipt_value):
            fail("Stage 0 benchmark or linked source changed during execution")
        summary: Dict[str, Any] = {
            "schema": SUMMARY_SCHEMA,
            "contract_sha256": contract["contract_sha256"],
            "source_receipt": source,
            "bench": bench_receipt,
            "bench_description": description,
            "v9_prerequisite": v9.receipt(),
            "records_file": RECORD_NAME,
            "records_sha256": records_sha256,
            "record_count": record_count,
            **verdict,
            "stage0_only": True,
            "promotion_evidence": False,
            "fresh_roots_used": False,
            "timing_evidence": False,
        }
        summary["summary_sha256"] = self_hash(summary, "summary_sha256")
        screen._write_exclusive(
            temporary / SUMMARY_NAME,
            (canonical_json(summary) + "\n").encode("ascii"))
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                v9.receipt() != v9.receipt_value):
            fail("Stage 0 benchmark or linked source changed before publication")
        _fsync_directory(temporary)
        screen._publish_directory_noreplace(temporary, output_dir)
        published = True
        _fsync_directory(output_dir.parent)
        return summary
    finally:
        if not published and temporary.exists():
            shutil.rmtree(temporary)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("describe", help="print the frozen Stage 0 contract")
    run_parser = subparsers.add_parser(
        "run", help="run the 1,536-observation Stage 0")
    run_parser.add_argument("--bench", type=Path, required=True)
    run_parser.add_argument("--v9-dir", type=Path, required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        _validate_constants()
        if arguments.command == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run_stage0(
            arguments.bench, arguments.v9_dir, arguments.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (Stage0Error, screen.ShortScreenError) as exc:
        print("wh2 MIX2 packet-offset Stage 0: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
