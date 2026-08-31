#!/usr/bin/env python3
"""Consumed-root filter for the remaining fixed MIX2 packet offsets.

Stage 2 authenticates the sealed Stage-0 survivor roster and the sealed
Stage-1 delta-2 introductions before executing anything.  It replays delta 2
as a failure control, screens every other survivor on exactly those 44
coordinates, and confirms the lowest zero-weak delta on the 47-coordinate
union.  This is recovery evidence only: it uses no fresh roots or timing and
cannot modify or license the production codec.

Unlike the earlier stages, every child stdout/stderr byte is durably published
in a bounded canonical base64 JSONL stream.  Publication reparses that stream
and requires byte-for-byte agreement with every parsed result record.
"""

from __future__ import annotations

import argparse
import base64
import binascii
from contextlib import ExitStack
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import (Any, Dict, Iterable, List, Mapping, Optional,
                    Sequence, Tuple)

try:
    from bench import wh2_mix2_packet_offset_stage0 as stage0
    from bench import wh2_mix2_packet_offset_stage1 as stage1
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_packet_offset_stage0 as stage0
    import wh2_mix2_packet_offset_stage1 as stage1
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage2-contract.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage2-record.v1"
RAW_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage2-raw.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage2-summary.v1"
BENCH_DESCRIPTION_SCHEMA = "wirehair.wh2.v2-bench-description.v1"
RECORD_NAME = "mix2-packet-offset-stage2-records.jsonl"
RAW_NAME = "mix2-packet-offset-stage2-raw.jsonl"
SUMMARY_NAME = "mix2-packet-offset-stage2-summary.json"

CANDIDATE_PAIR = "01"
CONTROL_DELTA = 2
PACKET_ATTEMPT_STRIDE = 0x9e3779b9
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
BINARY_DENSE_ROWS = 12
GF256_HEAVY_ROWS = 12
ZERO_SHA256 = "0" * 64

# Exact Stage-0 list.  Do not infer this as "every delta except failures": the
# literal roster and its canonical digest are part of the Stage-2 contract.
SURVIVOR_DELTAS = (
    2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22, 23,
    24, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60,
    63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 76, 77, 81, 82, 83, 84,
    85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97, 98, 99, 100, 102, 103,
    107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
    120, 121, 122, 124, 125, 126, 127, 128, 130, 131, 132, 133, 134,
    135, 137, 138, 139, 140, 141, 142, 144, 145, 147, 148, 149, 150,
    151, 152, 153, 154, 155, 156, 157, 159, 160, 161, 162, 163, 164,
    166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
    193, 195, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 218, 219, 220, 221,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236,
    237, 238, 239, 240, 241, 242, 243, 244, 246, 247, 248, 249, 250,
    253, 254, 255,
)
EXPECTED_SURVIVOR_DELTAS_SHA256 = (
    "4144abe36b10f00654c6bc1c3d9cae8b2095aeb2771daed9de937e78cd1b86a0")
SCREEN_DELTAS = tuple(delta for delta in SURVIVOR_DELTAS
                      if delta != CONTROL_DELTA)
EXPECTED_SCREEN_DELTAS_SHA256 = (
    "1cd7148ff6d45d4155955d6e86dda9d61a9d50146d860e3fce013a4409169c9c")

# Full Stage-1 coordinate ordinals.  Stage-0's local ordinals 0,1,2 must never
# enter this set: local 1 and 2 collide with genuine Stage-1 introductions.
INTRODUCTION_COORDINATES = (
    1, 2, 5, 8, 10, 14, 16, 20, 23, 26, 28, 31, 32, 35, 38, 44, 50,
    56, 59, 109, 112, 115, 118, 121, 127, 129, 130, 134, 136, 140, 142,
    146, 152, 161, 167, 175, 177, 181, 218, 224, 751, 851, 894, 1021,
)
ORIGINAL_WEAK_COORDINATES = (30, 573, 820)
CONFIRMATION_COORDINATES = tuple(sorted(
    set(INTRODUCTION_COORDINATES) | set(ORIGINAL_WEAK_COORDINATES)))
EXPECTED_INTRODUCTION_ORDINALS_SHA256 = (
    "86f4cd7b65c5461bf23dd8627c8f5fd8856c5412a2b0e574fdc6fd6b25b36b66")
EXPECTED_ORIGINAL_ORDINALS_SHA256 = (
    "af5600c0777e28d808933b5efd225bed8848f884ce2662b6eca7a6849d37329d")
EXPECTED_CONFIRMATION_ORDINALS_SHA256 = (
    "1c71a0ec8bf9721cd017e787d93eaaf762e4877096454346d7bc29f2b89bd553")
EXPECTED_STAGE1_INTRODUCTION_PROJECTION_SHA256 = (
    "6161035afed5425439d8e572607173378ae3987398815629f2db638bc434ca48")
EXPECTED_STAGE1_BASELINE_WEAK_PROJECTION_SHA256 = (
    "dfe88904adbcfd32032588f64ccf866c0369e1f09f50a59b44bf3cd0e9da0685")

CONTROL_OBSERVATION_COUNT = len(INTRODUCTION_COORDINATES)
SCREEN_OBSERVATION_COUNT = len(SCREEN_DELTAS) * len(INTRODUCTION_COORDINATES)
CONFIRMATION_OBSERVATION_COUNT = len(CONFIRMATION_COORDINATES)
PRE_CONFIRMATION_OBSERVATION_COUNT = (
    CONTROL_OBSERVATION_COUNT + SCREEN_OBSERVATION_COUNT)
MAX_OBSERVATION_COUNT = (
    PRE_CONFIRMATION_OBSERVATION_COUNT + CONFIRMATION_OBSERVATION_COUNT)

EXPECTED_CANDIDATE_PROFILE_SHA256 = stage1.EXPECTED_CANDIDATE_PROFILE_SHA256
EXPECTED_STAGE0_CONTRACT_SHA256 = stage1.EXPECTED_STAGE0_CONTRACT_SHA256
EXPECTED_STAGE1_CONTRACT_SHA256 = stage1.EXPECTED_CONTRACT_SHA256
EXPECTED_STAGE0_HELPER_SHA256 = stage1.EXPECTED_STAGE0_HELPER_SHA256
EXPECTED_STAGE1_HELPER_SHA256 = (
    "8d0eef6a9d6d073c1fdd84ec18f86ca1c3fe03ccff1f474b7cc46108ef8b4f01")
EXPECTED_STAGE0_SOURCE_COMMIT = stage1.EXPECTED_STAGE0_SOURCE_COMMIT
EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 = (
    stage1.EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256)
EXPECTED_STAGE0_BENCH_SHA256 = stage1.EXPECTED_STAGE0_BENCH_SHA256
EXPECTED_STAGE0_BENCH_SIZE = stage1.EXPECTED_STAGE0_BENCH_SIZE
EXPECTED_STAGE1_SOURCE_COMMIT = "16a20c0ed88fe8670e8a64df3847ae97092849bc"
EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256 = (
    "4be1b178dffa298d30225473153087ffdfa4244975171a0dc4a506d0127c4049")

STAGE0_FILES = dict(stage1.STAGE0_FILES)
STAGE1_FILES = {
    "mix2-packet-offset-stage1-records.jsonl": (
        5284872,
        "dbfa86fe1bf2457d3649e0ebe067be40d823be4edecaed72961b39f1f39d3b0a"),
    "mix2-packet-offset-stage1-summary.json": (
        61836,
        "fca9f5aea567684d8acf1e5e7b22dc0481aba7fc610b8a30771afcb50b0dbd80"),
}
EXPECTED_STAGE0_SUMMARY_SELF_SHA256 = (
    stage1.EXPECTED_STAGE0_SUMMARY_SELF_SHA256)
EXPECTED_STAGE1_SUMMARY_SELF_SHA256 = (
    "7d00f0d1afaf4240c1ddbaa1bbeeb8b23c55d998160453a61bfd0d10b9c802b4")

CSV_HEADER = tuple(stage1.CSV_HEADER)
METADATA_ORDER = tuple(stage1.METADATA_ORDER)
MAX_STDOUT_BYTES = 16 * 1024
MAX_STDERR_BYTES = 64 * 1024
MAX_RAW_DECODED_BYTES = 64 * 1024 * 1024
MAX_RAW_FILE_BYTES = 96 * 1024 * 1024
MAX_RAW_JSONL_LINE_BYTES = 32 * 1024
PROCESS_TIMEOUT_SECONDS = 600
CHILD_ENVIRONMENT = {"LANG": "C", "LC_ALL": "C", "TZ": "UTC"}
UINT64_MAX = (1 << 64) - 1

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent
SOURCE_PATHS = (
    "CMakeLists.txt", "codec/CMakeLists.txt",
    "bench/test_wh2_mix2_packet_offset_stage2.py",
    "bench/wh2_mix2_packet_offset_stage2.py",
    "bench/test_wh2_mix2_packet_offset_stage1.py",
    "bench/wh2_mix2_packet_offset_stage1.py",
    "bench/test_wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_attempt_crossfit.py",
    "bench/wh2_mix2_promotion_short_screen.py",
    "wirehair.cpp", "include/wirehair/wirehair.h",
    "include/wirehair/wirehair.hpp", "gf256.cpp", "gf256.h",
    "WirehairCodec.cpp", "WirehairCodec.h", "WirehairCompiler.h",
    "WirehairEnvironment.h", "WirehairHeavy.h", "WirehairTools.cpp",
    "WirehairTools.h", "WirehairDenseFixups.inc", "WirehairPeelFixups.inc",
    "codec/WirehairV2Bench.cpp", "codec/WirehairV2Codec.cpp",
    "codec/WirehairV2Codec.h", "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Peel.h", "codec/WirehairV2Plan.cpp",
    "codec/WirehairV2Plan.h", "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Policy.h", "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h", "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeDecode.h", "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2PrecodeEncode.h", "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2RawSeed.h", "codec/WirehairV2Seeds.cpp",
    "codec/WirehairV2Seeds.h", "codec/WirehairV2Solve.cpp",
    "codec/WirehairV2Solve.h",
)
SOURCE_STATUS_PATHS = (
    "CMakeLists.txt", "codec", "include", "wirehair.cpp", "gf256.cpp",
    "gf256.h", "WirehairCodec.cpp", "WirehairCodec.h",
    "WirehairCompiler.h", "WirehairEnvironment.h", "WirehairHeavy.h",
    "WirehairTools.cpp", "WirehairTools.h", "WirehairDenseFixups.inc",
    "WirehairPeelFixups.inc", "bench/wh2_mix2_packet_offset_stage2.py",
    "bench/test_wh2_mix2_packet_offset_stage2.py",
    "bench/wh2_mix2_packet_offset_stage1.py",
    "bench/test_wh2_mix2_packet_offset_stage1.py",
    "bench/wh2_mix2_packet_offset_stage0.py",
    "bench/test_wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_attempt_crossfit.py",
    "bench/wh2_mix2_promotion_short_screen.py",
)

COMMIT = re.compile(r"[0-9a-f]{40}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
HEX64 = re.compile(r"0x[0-9a-f]{16}\Z")
HEX32 = re.compile(r"0x[0-9a-f]{8}\Z")

# Replaced after the literal body is final.
EXPECTED_CONTRACT_SHA256 = (
    "91c1a0f1a9f97ac6c27baf12e6e1d112798e74d26246734e98c8712572cd449a")


class Stage2Error(RuntimeError):
    """Stage 2 cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise Stage2Error(message)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("ascii"))


def self_hash(value: Mapping[str, Any], field: str) -> str:
    body = dict(value)
    body.pop(field, None)
    return sha256_json(body)


def _exact_fields(value: Mapping[str, Any], fields: Iterable[str],
                  context: str) -> None:
    if type(value) is not dict or set(value) != set(fields):
        fail("{} fields differ from the frozen schema".format(context))


def _uint(value: Any, maximum: int, context: str) -> int:
    if type(value) is not int or value < 0 or value > maximum:
        fail("{} is not a bounded unsigned integer".format(context))
    return value


def _sha(value: Any, context: str) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        fail("{} is not a lowercase SHA-256".format(context))
    return value


def _contract_body() -> Mapping[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "bead": "wirehair-sxvz.16.1.20.38.3",
        "scope": "remaining fixed global packet offsets on sealed failures",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "candidate": {
            "base_profile_sha256": EXPECTED_CANDIDATE_PROFILE_SHA256,
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "mix_pair": CANDIDATE_PAIR,
            "binary_dense_rows": BINARY_DENSE_ROWS,
            "gf256_heavy_rows": GF256_HEAVY_ROWS,
            "construction_seed_basis": "production-profile",
            "precode_attempt_rule": "exact frozen v9 p(K)",
            "packet_attempt_rule": "q=(p+delta) mod 256",
            "packet_seed_rule": (
                "effective_q=(effective_p+(q-p)*0x9e3779b9) mod 2^32"),
            "proposed_implementation_complexity": (
                "hypothesis requiring a later production implementation "
                "audit: one immutable profile constant and one "
                "configuration-time add/mask, with no map, row, allocation, "
                "retry, per-packet branch, or asymptotic work; Stage2 is not "
                "proof that production meets this condition"),
            "test_hook_only": True,
            "production_split_pq_implemented": False,
        },
        "domain": {
            "stage0_survivor_deltas": list(SURVIVOR_DELTAS),
            "stage0_survivor_deltas_sha256":
                EXPECTED_SURVIVOR_DELTAS_SHA256,
            "failure_control_delta": CONTROL_DELTA,
            "screen_deltas": list(SCREEN_DELTAS),
            "screen_deltas_sha256": EXPECTED_SCREEN_DELTAS_SHA256,
            "screen_delta_count": len(SCREEN_DELTAS),
            "introduction_coordinate_ordinals":
                list(INTRODUCTION_COORDINATES),
            "introduction_ordinals_sha256":
                EXPECTED_INTRODUCTION_ORDINALS_SHA256,
            "original_weak_coordinate_ordinals":
                list(ORIGINAL_WEAK_COORDINATES),
            "original_ordinals_sha256": EXPECTED_ORIGINAL_ORDINALS_SHA256,
            "confirmation_coordinate_ordinals":
                list(CONFIRMATION_COORDINATES),
            "confirmation_ordinals_sha256":
                EXPECTED_CONFIRMATION_ORDINALS_SHA256,
            "control_observations": CONTROL_OBSERVATION_COUNT,
            "screen_observations": SCREEN_OBSERVATION_COUNT,
            "confirmation_observations_if_selected":
                CONFIRMATION_OBSERVATION_COUNT,
            "maximum_observations": MAX_OBSERVATION_COUNT,
            "order": (
                "delta2 control over introductions; every other survivor in "
                "ascending delta then introduction ordinal; selected "
                "confirmation in ascending full coordinate ordinal"),
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "trials": 1,
            "threads": 1,
            "full_payload_solve": True,
            "timing_used": False,
        },
        "prerequisites": {
            "stage0": {
                "contract_sha256": EXPECTED_STAGE0_CONTRACT_SHA256,
                "helper_sha256": EXPECTED_STAGE0_HELPER_SHA256,
                "source_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
                "source_receipt_sha256":
                    EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256,
                "summary_self_sha256":
                    EXPECTED_STAGE0_SUMMARY_SELF_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in STAGE0_FILES.items()
                },
                "required_disposition": "PASS",
                "required_survivors": list(SURVIVOR_DELTAS),
                "required_stage1_candidate_delta": CONTROL_DELTA,
            },
            "stage1": {
                "contract_sha256": EXPECTED_STAGE1_CONTRACT_SHA256,
                "helper_sha256": EXPECTED_STAGE1_HELPER_SHA256,
                "source_commit": EXPECTED_STAGE1_SOURCE_COMMIT,
                "source_receipt_sha256":
                    EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256,
                "summary_self_sha256":
                    EXPECTED_STAGE1_SUMMARY_SELF_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in STAGE1_FILES.items()
                },
                "required_disposition": "REJECT",
                "introduction_projection_sha256":
                    EXPECTED_STAGE1_INTRODUCTION_PROJECTION_SHA256,
                "baseline_weak_projection_sha256":
                    EXPECTED_STAGE1_BASELINE_WEAK_PROJECTION_SHA256,
                "required_introductions": 44,
                "required_original_weak": 3,
            },
        },
        "selection": {
            "failure_control": (
                "delta2 must reproduce all 44 sealed deterministic failures "
                "before survivor screening"),
            "screen": (
                "all 217 remaining survivors run once on all 44 sealed "
                "introductions; no early stop"),
            "rule": "lowest numeric delta with 44/44 success",
            "work_blind": True,
            "confirmation": (
                "selected delta runs once on the full 47-cell union; the 44 "
                "screen rows repeat exactly and the three original rows equal "
                "their authenticated Stage0 success projections"),
        },
        "gate": {
            "recovery": (
                "selected delta exists and confirmation has zero weak, errors, "
                "or configuration failures with exact projections"),
            "work": (
                "targeted 47-cell work is descriptive only because the three "
                "failed pair01 baselines contaminate comparison; the later "
                "1080-coordinate replay must require the selected offset no "
                "worse than diagonal pair01 for block XORs, GF256 muladds, "
                "and inactivation, preserving the already accepted "
                "1.705467-percent MIX2-versus-current XOR gain"),
            "overall": "PASS iff every Stage2 recovery gate passes",
        },
        "transport": {
            "benchmark": {
                "sha256": EXPECTED_STAGE0_BENCH_SHA256,
                "size": EXPECTED_STAGE0_BENCH_SIZE,
                "source_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
            },
            "historical_preflight": (
                "pin exact non-symlink files and fully authenticate Stage0 and "
                "Stage1 before the current source receipt, benchmark "
                "description, or workload"),
            "source_receipt": {
                "clean_at_HEAD": True,
                "working_source_bytes_equal_HEAD_blobs": True,
                "final_HEAD_recheck": True,
                "paths": list(SOURCE_PATHS),
            },
            "raw_output": {
                "schema": RAW_SCHEMA,
                "encoding": "canonical base64 of exact stdout/stderr bytes",
                "per_stdout_bound": MAX_STDOUT_BYTES,
                "per_stderr_bound": MAX_STDERR_BYTES,
                "aggregate_decoded_bound": MAX_RAW_DECODED_BYTES,
                "file_bound": MAX_RAW_FILE_BYTES,
                "publication_reparse_equality": True,
            },
            "child_environment": dict(CHILD_ENVIRONMENT),
            "process_timeout_seconds": PROCESS_TIMEOUT_SECONDS,
            "publication": (
                "canonical fsynced parsed and raw JSONL plus summary; reparse "
                "audit; fsynced staging; atomic no-replace rename; parent fsync"),
        },
        "stop_rule": (
            "PASS licenses only a separately frozen 1080-coordinate "
            "consumed-root replay (47 consulted coordinates and 1033 untouched "
            "coordinates, explicitly selection-leaked); that replay alone may "
            "license a separately "
            "frozen disjoint holdout.  REJECT retires the fixed-offset family. "
            "No timing, fresh roots, all-K work, production change, or "
            "promotion is licensed here."),
    }


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("MIX2 packet-offset Stage 2 contract differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _helper_hash(path: Path) -> str:
    return sha256_bytes(screen._stable_file_bytes(path, 4 * 1024 * 1024))


def _validate_constants() -> None:
    if (CANDIDATE_PAIR != "01" or CONTROL_DELTA != 2 or
            PACKET_ATTEMPT_STRIDE != 0x9e3779b9 or
            len(SURVIVOR_DELTAS) != 218 or
            tuple(sorted(set(SURVIVOR_DELTAS))) != SURVIVOR_DELTAS or
            SURVIVOR_DELTAS[0] != CONTROL_DELTA or
            len(SCREEN_DELTAS) != 217 or CONTROL_DELTA in SCREEN_DELTAS or
            sha256_json(list(SURVIVOR_DELTAS)) !=
                EXPECTED_SURVIVOR_DELTAS_SHA256 or
            sha256_json(list(SCREEN_DELTAS)) !=
                EXPECTED_SCREEN_DELTAS_SHA256):
        fail("Stage 2 survivor roster changed")
    if (len(INTRODUCTION_COORDINATES) != 44 or
            tuple(sorted(set(INTRODUCTION_COORDINATES))) !=
                INTRODUCTION_COORDINATES or
            len(ORIGINAL_WEAK_COORDINATES) != 3 or
            tuple(sorted(set(ORIGINAL_WEAK_COORDINATES))) !=
                ORIGINAL_WEAK_COORDINATES or
            set(INTRODUCTION_COORDINATES) &
                set(ORIGINAL_WEAK_COORDINATES) or
            len(CONFIRMATION_COORDINATES) != 47 or
            sha256_json(list(INTRODUCTION_COORDINATES)) !=
                EXPECTED_INTRODUCTION_ORDINALS_SHA256 or
            sha256_json(list(ORIGINAL_WEAK_COORDINATES)) !=
                EXPECTED_ORIGINAL_ORDINALS_SHA256 or
            sha256_json(list(CONFIRMATION_COORDINATES)) !=
                EXPECTED_CONFIRMATION_ORDINALS_SHA256):
        fail("Stage 2 full-coordinate roster changed")
    if (CONTROL_OBSERVATION_COUNT != 44 or
            SCREEN_OBSERVATION_COUNT != 9548 or
            PRE_CONFIRMATION_OBSERVATION_COUNT != 9592 or
            CONFIRMATION_OBSERVATION_COUNT != 47 or
            MAX_OBSERVATION_COUNT != 9639):
        fail("Stage 2 observation cardinality changed")
    if (tuple(CSV_HEADER) != tuple(stage0.CSV_HEADER) or
            tuple(CSV_HEADER) != tuple(screen.CSV_HEADER)):
        fail("precodefail CSV protocol changed")
    if screen.candidate_profile_sha256() != \
            EXPECTED_CANDIDATE_PROFILE_SHA256:
        fail("Stage 2 base candidate profile changed")
    if stage0._packet_seed_for_offset("0x12345678", 255, 3) != "0xb5610aa3":
        fail("wrapped arbitrary-delta packet seed rule changed")

    # Authenticate imported prerequisite code and contracts.  Stage1's own
    # validation recursively authenticates V9/cross-fit dependencies.
    stage0._validate_constants()
    stage1._validate_constants()
    if (stage0.contract_description()["contract_sha256"] !=
            EXPECTED_STAGE0_CONTRACT_SHA256 or
            stage1.contract_description()["contract_sha256"] !=
            EXPECTED_STAGE1_CONTRACT_SHA256):
        fail("an imported prerequisite contract changed")
    for path, expected in (
            (Path(stage0.__file__).resolve(), EXPECTED_STAGE0_HELPER_SHA256),
            (Path(stage1.__file__).resolve(), EXPECTED_STAGE1_HELPER_SHA256)):
        if _helper_hash(path) != expected:
            fail("imported helper differs from frozen source: {}".format(path))

    expected_local_to_full = ((0, 30), (1, 573), (2, 820))
    realized = []
    for local, old in enumerate(stage0.COORDINATES):
        matches = []
        for full in range(stage1.COORDINATE_COUNT):
            current = stage1.Coordinate(full)
            if (current.K, current.attempt, current.root_index, current.root,
                    current.schedule) == (
                    old.K, old.attempt, old.root_index, old.root,
                    old.schedule):
                matches.append(full)
        if len(matches) != 1:
            fail("Stage0 original coordinate has ambiguous Stage1 identity")
        realized.append((local, matches[0]))
    if tuple(realized) != expected_local_to_full:
        fail("Stage0 local-to-full coordinate normalization changed")


@dataclass(frozen=True)
class RunCoordinate:
    ordinal: int
    delta: int

    def __post_init__(self) -> None:
        if (type(self.ordinal) is not int or
                not 0 <= self.ordinal < stage1.COORDINATE_COUNT or
                type(self.delta) is not int or not 0 <= self.delta <= 255):
            fail("Stage 2 run coordinate is out of range")

    @property
    def base(self) -> stage1.Coordinate:
        return stage1.Coordinate(self.ordinal)

    @property
    def K(self) -> int:
        return self.base.K

    @property
    def attempt(self) -> int:
        return self.base.attempt

    @property
    def root_index(self) -> int:
        return self.base.root_index

    @property
    def root(self) -> str:
        return self.base.root

    @property
    def schedule_index(self) -> int:
        return self.base.schedule_index

    @property
    def schedule(self) -> str:
        return self.base.schedule

    @property
    def cell_ordinal(self) -> int:
        return self.base.cell_ordinal

    @property
    def packet_attempt(self) -> int:
        return (self.attempt + self.delta) & 255

    def identity(self) -> Mapping[str, Any]:
        return {
            "coordinate_ordinal": self.ordinal,
            "K": self.K,
            "attempt": self.attempt,
            "root_index": self.root_index,
            "loss_root": self.root,
            "schedule_index": self.schedule_index,
            "schedule": self.schedule,
            "cell_ordinal": self.cell_ordinal,
            "baseline_record_ordinal": self.base.baseline_record_ordinal,
        }


PHASE_CONTROL = "delta2-failure-control"
PHASE_SCREEN = "survivor-screen"
PHASE_CONFIRMATION = "selected-confirmation"


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    phase: str
    delta: int
    coordinate_ordinal: int

    def __post_init__(self) -> None:
        if type(self.ordinal) is not int or self.ordinal < 0:
            fail("Stage 2 invocation ordinal is invalid")
        if self.phase == PHASE_CONTROL:
            if (self.delta != CONTROL_DELTA or
                    self.coordinate_ordinal not in INTRODUCTION_COORDINATES or
                    self.ordinal != INTRODUCTION_COORDINATES.index(
                        self.coordinate_ordinal)):
                fail("delta2 control invocation is noncanonical")
        elif self.phase == PHASE_SCREEN:
            if (self.delta not in SCREEN_DELTAS or
                    self.coordinate_ordinal not in INTRODUCTION_COORDINATES):
                fail("survivor-screen invocation is outside the roster")
            expected = (CONTROL_OBSERVATION_COUNT +
                        SCREEN_DELTAS.index(self.delta) *
                        len(INTRODUCTION_COORDINATES) +
                        INTRODUCTION_COORDINATES.index(
                            self.coordinate_ordinal))
            if self.ordinal != expected:
                fail("survivor-screen invocation order changed")
        elif self.phase == PHASE_CONFIRMATION:
            if (self.delta not in SCREEN_DELTAS or
                    self.coordinate_ordinal not in CONFIRMATION_COORDINATES or
                    self.ordinal != PRE_CONFIRMATION_OBSERVATION_COUNT +
                        CONFIRMATION_COORDINATES.index(
                            self.coordinate_ordinal)):
                fail("selected-confirmation invocation is noncanonical")
        else:
            fail("Stage 2 invocation phase is unknown")

    @property
    def coordinate(self) -> RunCoordinate:
        return RunCoordinate(self.coordinate_ordinal, self.delta)

    def argv(self, bench: Path) -> List[str]:
        coordinate = self.coordinate
        return [
            str(bench), "precodefail", "--N", str(coordinate.K),
            "--bb-list", "2", "--overhead", "0", "--trials", "1",
            "--threads", "1", "--loss", "0.5", "--seed",
            coordinate.root, "--schedule", coordinate.schedule,
            "--heavy-family", "periodic", "--mix-count", "2",
            "--mix-pair", CANDIDATE_PAIR, "--binary-dense-rows", "12",
            "--gf256-heavy-rows", "12", "--dense-anchors", "two07",
            "--paired-overhead-stream", "--full-payload-solve",
            "--cold-solve-wide-xor", "policy", "--exact-precode-attempt",
            str(coordinate.attempt), "--exact-packet-attempt",
            str(coordinate.packet_attempt), "--construction-seed-basis",
            "production-profile",
        ]

    def identity(self) -> Mapping[str, Any]:
        return {
            "ordinal": self.ordinal,
            "phase": self.phase,
            "mix_pair": CANDIDATE_PAIR,
            "delta": self.delta,
            "precode_attempt": self.coordinate.attempt,
            "packet_attempt": self.coordinate.packet_attempt,
            **self.coordinate.identity(),
        }


def control_invocations() -> Tuple[Invocation, ...]:
    return tuple(Invocation(index, PHASE_CONTROL, CONTROL_DELTA, coordinate)
                 for index, coordinate in enumerate(INTRODUCTION_COORDINATES))


def screen_invocations() -> Tuple[Invocation, ...]:
    values = []
    for delta in SCREEN_DELTAS:
        for coordinate in INTRODUCTION_COORDINATES:
            values.append(Invocation(
                CONTROL_OBSERVATION_COUNT + len(values), PHASE_SCREEN,
                delta, coordinate))
    if len(values) != SCREEN_OBSERVATION_COUNT:
        fail("survivor-screen invocation roster changed cardinality")
    return tuple(values)


def confirmation_invocations(delta: int) -> Tuple[Invocation, ...]:
    return tuple(Invocation(
        PRE_CONFIRMATION_OBSERVATION_COUNT + index, PHASE_CONFIRMATION,
        delta, coordinate)
        for index, coordinate in enumerate(CONFIRMATION_COORDINATES))


STAGE1_RECORD_FIELDS = frozenset((
    "K", "attempt", "baseline_matrix_sha256", "baseline_record_ordinal",
    "bench_binary_sha256", "bench_source_git_commit", "binary_deficiency",
    "binary_dense_rows", "block_bytes", "block_xors", "cell_ordinal",
    "cold_solve_wide_xor", "command_sha256", "configuration_failure",
    "configuration_stderr_sha256", "construction_seed_basis",
    "controller_source_git_commit", "coordinate_ordinal", "delta",
    "dense_anchor_layout", "deterministic_projection_sha256",
    "effective_packet_seed", "effective_precode_seed", "error",
    "failure_trials", "first_rank_fail", "fresh_roots_used",
    "full_payload_solve", "gf256_heavy_rows", "gf256_muladds",
    "gf256_rank_gain", "heavy_family", "heavy_shortfall",
    "inactivated_columns", "loss_ppm", "loss_root", "mix_count",
    "mix_pair", "ordinal", "outcome", "overhead", "overhead_stream",
    "packet_attempt", "pair01_baseline", "pair01_baseline_record_sha256",
    "precode_attempt", "promotion_evidence", "rank_fail", "repetition",
    "returncode", "root_index", "schedule", "schedule_index", "schema",
    "source_hits", "source_receipt_sha256", "staircase", "stderr_sha256",
    "stdout_sha256", "success", "timing_evidence", "weak",
))

BASELINE_FIELDS = tuple(stage1.BASELINE_PROJECTION_FIELDS)
SCIENTIFIC_FIELDS = (
    "mix_pair", "delta", "precode_attempt", "packet_attempt",
    "coordinate_ordinal", "K", "root_index", "loss_root", "schedule",
    "block_bytes", "loss_ppm", "overhead", "dense_anchor_layout",
    "mix_count", "binary_dense_rows", "gf256_heavy_rows", "heavy_family",
    "construction_seed_basis", "full_payload_solve", "overhead_stream",
    "cold_solve_wide_xor", "returncode", "staircase", "source_hits",
    "outcome", "success", "rank_fail", "error", "configuration_failure",
    "weak", "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
    "heavy_shortfall", "first_rank_fail", "failure_trials",
    "effective_precode_seed", "effective_packet_seed", "block_xors",
    "gf256_muladds", "configuration_stderr_sha256",
)
DETERMINISTIC_FIELDS = (
    "phase", *SCIENTIFIC_FIELDS, "attempt", "schedule_index",
    "cell_ordinal", "baseline_record_ordinal", "command_sha256",
    "bench_binary_sha256", "bench_source_git_commit",
    "controller_source_git_commit", "source_receipt_sha256",
    "baseline_projection_sha256",
)
RECORD_FIELDS = frozenset((
    "schema", *DETERMINISTIC_FIELDS, "ordinal", "promotion_evidence",
    "fresh_roots_used", "timing_evidence", "stdout_sha256",
    "stderr_sha256", "deterministic_projection_sha256",
))
RAW_FIELDS = frozenset((
    "schema", "ordinal", "invocation", "command_sha256", "stdout_sha256",
    "stderr_sha256", "returncode", "stdout_base64", "stderr_base64",
    "parsed_record_sha256",
))


def _baseline_projection(value: Mapping[str, Any]) -> Mapping[str, Any]:
    if type(value) is not dict or set(value) != set(BASELINE_FIELDS):
        fail("pair01 baseline projection fields changed")
    return {field: value[field] for field in BASELINE_FIELDS}


def scientific_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in SCIENTIFIC_FIELDS if field not in record]
    if missing:
        fail("scientific record omits fields: {}".format(missing))
    return {field: record[field] for field in SCIENTIFIC_FIELDS}


def deterministic_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in DETERMINISTIC_FIELDS if field not in record]
    if missing:
        fail("Stage 2 record omits deterministic fields: {}".format(missing))
    return {field: record[field] for field in DETERMINISTIC_FIELDS}


def _terminal_record_is_valid(record: Mapping[str, Any], context: str) -> None:
    for field in ("success", "rank_fail", "error", "configuration_failure",
                  "weak"):
        if type(record.get(field)) is not bool:
            fail("{} terminal Boolean is malformed".format(context))
    if record["configuration_failure"]:
        if (record["outcome"] != "configuration_failure" or
                any(record[field] for field in
                    ("success", "rank_fail", "error")) or
                record["weak"] is not True):
            fail("{} configuration outcome is inconsistent".format(context))
        return
    if sum(int(record[field]) for field in
           ("success", "rank_fail", "error")) != 1:
        fail("{} terminal state is ambiguous".format(context))
    expected = ("success" if record["success"] else
                ("rank_fail" if record["rank_fail"] else "error"))
    if record["outcome"] != expected or record["weak"] is record["success"]:
        fail("{} outcome and weakness disagree".format(context))
    for field in ("inactivated_columns", "binary_deficiency",
                  "gf256_rank_gain", "heavy_shortfall", "block_xors",
                  "gf256_muladds", "staircase", "source_hits"):
        _uint(record.get(field), UINT64_MAX, "{} {}".format(context, field))
    if (record["gf256_rank_gain"] > record["binary_deficiency"] or
            record["binary_deficiency"] > record["inactivated_columns"] or
            record["gf256_rank_gain"] > GF256_HEAVY_ROWS or
            (record["success"] and record["gf256_rank_gain"] !=
             record["binary_deficiency"]) or
            (record["rank_fail"] and record["gf256_rank_gain"] >=
             record["binary_deficiency"])):
        fail("{} rank counters disagree".format(context))
    if (type(record.get("effective_precode_seed")) is not str or
            HEX64.fullmatch(record["effective_precode_seed"]) is None or
            type(record.get("effective_packet_seed")) is not str or
            HEX32.fullmatch(record["effective_packet_seed"]) is None):
        fail("{} effective seed is malformed".format(context))


def _stage0_global_coordinate(local: int) -> int:
    old = stage0.COORDINATES[local]
    matches = []
    for full in range(stage1.COORDINATE_COUNT):
        current = stage1.Coordinate(full)
        if (current.K, current.attempt, current.root_index, current.root,
                current.schedule) == (
                old.K, old.attempt, old.root_index, old.root, old.schedule):
            matches.append(full)
    if len(matches) != 1:
        fail("cannot normalize Stage0 local coordinate {}".format(local))
    return matches[0]


def _stage0_scientific_projection(
        record: Mapping[str, Any], full_coordinate: int) -> Mapping[str, Any]:
    coordinate = RunCoordinate(full_coordinate, record["delta"])
    augmented = dict(record)
    augmented.update({
        "coordinate_ordinal": full_coordinate,
        "attempt": coordinate.attempt,
        "schedule_index": coordinate.schedule_index,
        "cell_ordinal": coordinate.cell_ordinal,
        "baseline_record_ordinal": coordinate.base.baseline_record_ordinal,
    })
    return scientific_projection(augmented)


def _parse_canonical_json_bytes(data: bytes, context: str) \
        -> Mapping[str, Any]:
    if not data.endswith(b"\n") or b"\r" in data or data.count(b"\n") != 1:
        fail("{} is not one newline-canonical JSON object".format(context))
    try:
        text = data[:-1].decode("ascii")
        value = json.loads(text)
    except (UnicodeDecodeError, TypeError, ValueError) as exc:
        fail("{} is invalid JSON: {}".format(context, exc))
    if type(value) is not dict or text != canonical_json(value):
        fail("{} is not a canonical JSON object".format(context))
    return value


def _parse_canonical_json_line(line: bytes, context: str) \
        -> Mapping[str, Any]:
    if not line.endswith(b"\n") or b"\r" in line:
        fail("{} has invalid JSONL framing".format(context))
    try:
        text = line[:-1].decode("ascii")
        value = json.loads(text)
    except (UnicodeDecodeError, TypeError, ValueError) as exc:
        fail("{} is invalid JSON: {}".format(context, exc))
    if type(value) is not dict or text != canonical_json(value):
        fail("{} is not canonical JSON".format(context))
    return value


@dataclass
class HistoricalEvidence:
    stack: Optional[ExitStack]
    pins: Mapping[str, Mapping[str, stage1.PinnedFile]]
    baseline_by_coordinate: Mapping[int, Mapping[str, Any]]
    stage1_control_by_coordinate: Mapping[
        int, Tuple[Mapping[str, Any], Mapping[str, Any]]]
    stage0_success_by_delta_coordinate: Mapping[
        Tuple[int, int], Mapping[str, Any]]
    receipt_value: Mapping[str, Any]

    def file_receipts(self) -> Mapping[str, Any]:
        return {
            group: {
                name: pin.receipt() for name, pin in sorted(values.items())
            }
            for group, values in sorted(self.pins.items())
        }

    def receipt(self) -> Mapping[str, Any]:
        value = dict(self.receipt_value)
        if self.pins:
            value["files"] = self.file_receipts()
        return value

    def close(self) -> None:
        if self.stack is not None:
            self.stack.close()
            self.stack = None

    def __enter__(self) -> "HistoricalEvidence":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def _pin_group(stack: ExitStack, directory: Path,
               specs: Mapping[str, Tuple[int, str]], label: str) \
        -> Mapping[str, stage1.PinnedFile]:
    root = stage1._validate_historical_directory(directory, label)
    return {
        name: stack.enter_context(stage1._open_pinned_file(
            root / name, size, digest, "{} {}".format(label, name)))
        for name, (size, digest) in specs.items()
    }


def _load_stage0(pins: Mapping[str, stage1.PinnedFile]) \
        -> Tuple[Mapping[int, Mapping[str, Any]],
                 Mapping[Tuple[int, int], Mapping[str, Any]],
                 Mapping[str, Any]]:
    summary = _parse_canonical_json_bytes(
        pins["mix2-packet-offset-stage0-summary.json"].read_bytes(
            1024 * 1024), "Stage0 summary")
    if (summary.get("summary_sha256") !=
            EXPECTED_STAGE0_SUMMARY_SELF_SHA256 or
            self_hash(summary, "summary_sha256") !=
                EXPECTED_STAGE0_SUMMARY_SELF_SHA256 or
            summary.get("schema") != stage0.SUMMARY_SCHEMA or
            summary.get("contract_sha256") !=
                EXPECTED_STAGE0_CONTRACT_SHA256 or
            summary.get("records_sha256") != STAGE0_FILES[
                "mix2-packet-offset-stage0-records.jsonl"][1] or
            summary.get("record_count") != stage0.EXPECTED_INVOCATION_COUNT or
            summary.get("disposition") != "PASS" or
            summary.get("stage1_candidate_delta") != CONTROL_DELTA or
            summary.get("stage0_only") is not True or
            summary.get("promotion_evidence") is not False or
            summary.get("fresh_roots_used") is not False or
            summary.get("timing_evidence") is not False or
            tuple(summary.get("survivors", ())) != SURVIVOR_DELTAS or
            sha256_json(summary.get("survivors")) !=
                EXPECTED_SURVIVOR_DELTAS_SHA256):
        fail("Stage0 summary provenance or survivor roster changed")
    source = summary.get("source_receipt")
    bench = summary.get("bench")
    description = summary.get("bench_description")
    if (type(source) is not dict or
            source.get("source_git_commit") != EXPECTED_STAGE0_SOURCE_COMMIT or
            source.get("source_receipt_sha256") !=
                EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 or
            self_hash(source, "source_receipt_sha256") !=
                EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 or
            type(bench) is not dict or
            bench.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench.get("size") != EXPECTED_STAGE0_BENCH_SIZE or
            type(description) is not dict or
            description.get("schema") != BENCH_DESCRIPTION_SCHEMA or
            description.get("source_git_commit") !=
                EXPECTED_STAGE0_SOURCE_COMMIT):
        fail("Stage0 source or benchmark receipt changed")

    records = []
    for ordinal, line in pins[
            "mix2-packet-offset-stage0-records.jsonl"].lines():
        row = _parse_canonical_json_line(
            line, "Stage0 record {}".format(ordinal))
        if (row.get("schema") != stage0.RECORD_SCHEMA or
                row.get("ordinal") != ordinal or
                row.get("bench_binary_sha256") !=
                    EXPECTED_STAGE0_BENCH_SHA256 or
                row.get("source_git_commit") !=
                    EXPECTED_STAGE0_SOURCE_COMMIT or
                row.get("source_receipt_sha256") !=
                    EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 or
                row.get("deterministic_projection_sha256") !=
                    stage0.sha256_json(stage0.deterministic_projection(row))):
            fail("Stage0 record {} receipt changed".format(ordinal))
        records.append(row)
    if len(records) != stage0.EXPECTED_INVOCATION_COUNT:
        fail("Stage0 record roster has the wrong cardinality")

    precode_seeds = summary.get("effective_precode_seeds")
    packet_seeds = summary.get("authenticated_effective_p_packet_seeds")
    diagonal = summary.get("v9_prerequisite", {}).get("diagonal_control")
    if (type(precode_seeds) is not list or len(precode_seeds) != 3 or
            type(packet_seeds) is not list or len(packet_seeds) != 3 or
            type(diagonal) is not list or len(diagonal) != 3):
        fail("Stage0 authenticated seed/control maps changed")
    verdict = stage0.adjudicate(
        records, dict(enumerate(precode_seeds)),
        dict(enumerate(packet_seeds)), dict(enumerate(diagonal)))
    for field in ("diagonal_control", "delta_results", "survivors",
                  "stage1_candidate_delta", "effective_precode_seeds",
                  "authenticated_effective_p_packet_seeds", "disposition"):
        if summary.get(field) != verdict[field]:
            fail("Stage0 independent adjudication changed {}".format(field))

    by_key = {(row["delta"], row["coordinate_ordinal"], row["repetition"]):
              row for row in records}
    baseline: Dict[int, Mapping[str, Any]] = {}
    success: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    for local in range(3):
        full = _stage0_global_coordinate(local)
        diagonal_row = by_key[(0, local, 0)]
        projection = {
            "K": diagonal_row["K"],
            "attempt": diagonal_row["precode_attempt"],
            "root_index": diagonal_row["root_index"],
            "loss_root": diagonal_row["loss_root"],
            "schedule_index": stage1.Coordinate(full).schedule_index,
            "schedule": diagonal_row["schedule"],
            "cell_ordinal": stage1.Coordinate(full).cell_ordinal,
            "outcome": diagonal_row["outcome"],
            "success": diagonal_row["success"],
            "rank_fail": diagonal_row["rank_fail"],
            "error": diagonal_row["error"],
            "binary_deficiency": diagonal_row["binary_deficiency"],
            "gf256_rank_gain": diagonal_row["gf256_rank_gain"],
            "inactivated_columns": diagonal_row["inactivated_columns"],
            "heavy_shortfall": diagonal_row["heavy_shortfall"],
            "effective_precode_seed": diagonal_row["effective_precode_seed"],
            "effective_packet_seed": diagonal_row["effective_packet_seed"],
            "block_xors": diagonal_row["block_xors"],
            "gf256_muladds": diagonal_row["gf256_muladds"],
        }
        baseline[full] = _baseline_projection(projection)
        for delta in SURVIVOR_DELTAS:
            first = by_key[(delta, local, 0)]
            second = by_key[(delta, local, 1)]
            if (scientific_projection({
                    **_stage0_scientific_projection(first, full)}) !=
                    scientific_projection({
                        **_stage0_scientific_projection(second, full)}) or
                    first["success"] is not True or second["success"] is not True):
                fail("Stage0 survivor projection is not an exact success")
            success[(delta, full)] = _stage0_scientific_projection(first, full)
    if tuple(sorted(baseline)) != ORIGINAL_WEAK_COORDINATES:
        fail("Stage0 original coordinates were not normalized to full ordinals")
    receipt = {
        "summary_self_sha256": EXPECTED_STAGE0_SUMMARY_SELF_SHA256,
        "records_sha256": STAGE0_FILES[
            "mix2-packet-offset-stage0-records.jsonl"][1],
        "record_count": len(records),
        "source_git_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
        "source_receipt_sha256": EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256,
        "survivor_deltas_sha256": EXPECTED_SURVIVOR_DELTAS_SHA256,
        "survivor_count": len(SURVIVOR_DELTAS),
        "local_to_full_coordinate_ordinals": [[0, 30], [1, 573], [2, 820]],
    }
    return baseline, success, receipt


def _validate_stage1_record(row: Mapping[str, Any], ordinal: int) -> None:
    _exact_fields(row, STAGE1_RECORD_FIELDS,
                  "Stage1 record {}".format(ordinal))
    invocation = stage1.Invocation(
        ordinal, ordinal // len(stage1.REPETITIONS),
        ordinal % len(stage1.REPETITIONS))
    coordinate = invocation.coordinate
    baseline = _baseline_projection(row["pair01_baseline"])
    if (row["schema"] != stage1.RECORD_SCHEMA or
            row["ordinal"] != ordinal or
            row["coordinate_ordinal"] != invocation.coordinate_ordinal or
            row["repetition"] != invocation.repetition or
            any(row.get(field) != value for field, value in
                coordinate.identity().items()) or
            row["mix_pair"] != CANDIDATE_PAIR or
            row["delta"] != CONTROL_DELTA or
            row["precode_attempt"] != coordinate.attempt or
            row["packet_attempt"] != (coordinate.attempt + 2) & 255 or
            row["block_bytes"] != BLOCK_BYTES or
            row["loss_ppm"] != LOSS_PPM or row["overhead"] != OVERHEAD or
            row["dense_anchor_layout"] != "two07" or
            row["mix_count"] != 2 or
            row["binary_dense_rows"] != BINARY_DENSE_ROWS or
            row["gf256_heavy_rows"] != GF256_HEAVY_ROWS or
            row["heavy_family"] != "periodic" or
            row["construction_seed_basis"] != "production-profile" or
            row["full_payload_solve"] is not True or
            row["overhead_stream"] != "paired" or
            row["cold_solve_wide_xor"] != "policy" or
            row["promotion_evidence"] is not False or
            row["fresh_roots_used"] is not False or
            row["timing_evidence"] is not False or
            row["bench_binary_sha256"] != EXPECTED_STAGE0_BENCH_SHA256 or
            row["bench_source_git_commit"] != EXPECTED_STAGE0_SOURCE_COMMIT or
            row["controller_source_git_commit"] !=
                EXPECTED_STAGE1_SOURCE_COMMIT or
            row["source_receipt_sha256"] !=
                EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256 or
            row["deterministic_projection_sha256"] !=
                stage1.sha256_json(stage1.deterministic_projection(row))):
        fail("Stage1 record {} identity or receipt changed".format(ordinal))
    reconstructed = dict(baseline)
    reconstructed["record_ordinal"] = coordinate.baseline_record_ordinal
    stage1._validate_candidate_seed_receipt(
        row, reconstructed, invocation.coordinate_ordinal)
    _terminal_record_is_valid(row, "Stage1 record {}".format(ordinal))
    for field in ("command_sha256", "stdout_sha256", "stderr_sha256",
                  "pair01_baseline_record_sha256"):
        _sha(row[field], "Stage1 record {} {}".format(ordinal, field))


def _load_stage1(pins: Mapping[str, stage1.PinnedFile]) \
        -> Tuple[Mapping[int, Mapping[str, Any]],
                 Mapping[int, Tuple[Mapping[str, Any], Mapping[str, Any]]],
                 Mapping[str, Any]]:
    summary = _parse_canonical_json_bytes(
        pins["mix2-packet-offset-stage1-summary.json"].read_bytes(
            1024 * 1024), "Stage1 summary")
    if (summary.get("summary_sha256") !=
            EXPECTED_STAGE1_SUMMARY_SELF_SHA256 or
            self_hash(summary, "summary_sha256") !=
                EXPECTED_STAGE1_SUMMARY_SELF_SHA256 or
            summary.get("schema") != stage1.SUMMARY_SCHEMA or
            summary.get("contract_sha256") !=
                EXPECTED_STAGE1_CONTRACT_SHA256 or
            summary.get("records_sha256") != STAGE1_FILES[
                "mix2-packet-offset-stage1-records.jsonl"][1] or
            summary.get("record_count") != stage1.OBSERVATION_COUNT or
            summary.get("disposition") != "REJECT" or
            summary.get("candidate_packet_delta") != CONTROL_DELTA or
            summary.get("candidate_pair") != CANDIDATE_PAIR or
            summary.get("candidate_errors") != 0 or
            summary.get("candidate_configuration_failures") != 0 or
            summary.get("stage1_only") is not True or
            summary.get("promotion_evidence") is not False or
            summary.get("fresh_roots_used") is not False or
            summary.get("timing_evidence") is not False):
        fail("Stage1 summary provenance changed")
    source = summary.get("source_receipt")
    bench = summary.get("bench")
    description = summary.get("bench_description")
    if (type(source) is not dict or
            source.get("source_git_commit") != EXPECTED_STAGE1_SOURCE_COMMIT or
            source.get("source_receipt_sha256") !=
                EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256 or
            self_hash(source, "source_receipt_sha256") !=
                EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256 or
            type(bench) is not dict or
            bench.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench.get("size") != EXPECTED_STAGE0_BENCH_SIZE or
            type(description) is not dict or
            description.get("schema") != BENCH_DESCRIPTION_SCHEMA or
            description.get("source_git_commit") !=
                EXPECTED_STAGE0_SOURCE_COMMIT):
        fail("Stage1 source or benchmark receipt changed")
    introductions = summary.get("introduced_weak_coordinates")
    baseline_weak = summary.get("baseline_weak_coordinates")
    if (type(introductions) is not list or len(introductions) != 44 or
            sha256_json(introductions) !=
                EXPECTED_STAGE1_INTRODUCTION_PROJECTION_SHA256 or
            tuple(item.get("coordinate_ordinal") for item in introductions) !=
                INTRODUCTION_COORDINATES or
            type(baseline_weak) is not list or len(baseline_weak) != 3 or
            sha256_json(baseline_weak) !=
                EXPECTED_STAGE1_BASELINE_WEAK_PROJECTION_SHA256 or
            tuple(item.get("coordinate_ordinal") for item in baseline_weak) !=
                ORIGINAL_WEAK_COORDINATES):
        fail("Stage1 weak-coordinate projections changed")

    by_coordinate: Dict[int, List[Mapping[str, Any]]] = {}
    baseline_by_coordinate: Dict[int, Mapping[str, Any]] = {}
    count = 0
    for ordinal, line in pins[
            "mix2-packet-offset-stage1-records.jsonl"].lines():
        row = _parse_canonical_json_line(
            line, "Stage1 record {}".format(ordinal))
        _validate_stage1_record(row, ordinal)
        coordinate = row["coordinate_ordinal"]
        baseline = _baseline_projection(row["pair01_baseline"])
        previous = baseline_by_coordinate.setdefault(coordinate, baseline)
        if previous != baseline:
            fail("Stage1 repeated baseline projection changed")
        by_coordinate.setdefault(coordinate, []).append(row)
        count += 1
    if count != stage1.OBSERVATION_COUNT or len(by_coordinate) != 1080:
        fail("Stage1 record roster has the wrong cardinality")

    candidate_weak = []
    computed_baseline_weak = []
    repairs = []
    computed_introductions = []
    controls: Dict[
        int, Tuple[Mapping[str, Any], Mapping[str, Any]]] = {}
    for coordinate in range(stage1.COORDINATE_COUNT):
        rows = by_coordinate.get(coordinate, [])
        if len(rows) != 2 or stage1.deterministic_projection(rows[0]) != \
                stage1.deterministic_projection(rows[1]):
            fail("Stage1 repeated deterministic projection changed")
        baseline = dict(baseline_by_coordinate[coordinate])
        baseline["record_ordinal"] = \
            stage1.Coordinate(coordinate).baseline_record_ordinal
        receipt = stage1._coordinate_result(coordinate, rows, baseline)
        candidate_success = all(row["success"] for row in rows)
        if not candidate_success:
            candidate_weak.append(receipt)
        if not baseline["success"]:
            computed_baseline_weak.append(receipt)
        if not baseline["success"] and candidate_success:
            repairs.append(receipt)
        if baseline["success"] and not candidate_success:
            computed_introductions.append(receipt)
        if coordinate in INTRODUCTION_COORDINATES:
            if (rows[0]["rank_fail"] is not True or
                    rows[1]["rank_fail"] is not True or
                    baseline["success"] is not True):
                fail("Stage1 introduction is not a repeated rank failure")
            controls[coordinate] = (
                scientific_projection(rows[0]),
                scientific_projection(rows[1]))
    if (candidate_weak != summary.get("candidate_weak_coordinates") or
            computed_introductions != introductions or
            computed_baseline_weak != baseline_weak or
            repairs != summary.get("repaired_coordinates") or
            tuple(sorted(controls)) != INTRODUCTION_COORDINATES):
        fail("Stage1 independent weak-coordinate adjudication changed")
    receipt = {
        "summary_self_sha256": EXPECTED_STAGE1_SUMMARY_SELF_SHA256,
        "records_sha256": STAGE1_FILES[
            "mix2-packet-offset-stage1-records.jsonl"][1],
        "record_count": count,
        "source_git_commit": EXPECTED_STAGE1_SOURCE_COMMIT,
        "source_receipt_sha256": EXPECTED_STAGE1_SOURCE_RECEIPT_SHA256,
        "introduction_projection_sha256":
            EXPECTED_STAGE1_INTRODUCTION_PROJECTION_SHA256,
        "baseline_weak_projection_sha256":
            EXPECTED_STAGE1_BASELINE_WEAK_PROJECTION_SHA256,
        "introduction_count": len(controls),
    }
    return ({coordinate: baseline_by_coordinate[coordinate]
             for coordinate in INTRODUCTION_COORDINATES}, controls, receipt)


def load_historical_evidence(stage0_dir: Path,
                             stage1_dir: Path) -> HistoricalEvidence:
    stack = ExitStack()
    try:
        pins = {
            "stage0": _pin_group(stack, stage0_dir, STAGE0_FILES, "Stage0"),
            "stage1": _pin_group(stack, stage1_dir, STAGE1_FILES, "Stage1"),
        }
        original_baseline, stage0_success, stage0_receipt = \
            _load_stage0(pins["stage0"])
        intro_baseline, controls, stage1_receipt = \
            _load_stage1(pins["stage1"])
        baseline = {**original_baseline, **intro_baseline}
        if (tuple(sorted(baseline)) != CONFIRMATION_COORDINATES or
                len(baseline) != 47):
            fail("historical 47-cell union is incomplete or colliding")
        receipt = {
            "stage0": stage0_receipt,
            "stage1": stage1_receipt,
            "full_coordinate_union": list(CONFIRMATION_COORDINATES),
            "full_coordinate_union_sha256":
                sha256_json(list(CONFIRMATION_COORDINATES)),
        }
        result = HistoricalEvidence(
            stack, pins, baseline, controls, stage0_success, receipt)
        result.receipt()
        return result
    except BaseException:
        stack.close()
        raise


@dataclass(frozen=True)
class ProcessResult:
    invocation: Invocation
    command_sha256: str
    stdout_sha256: str
    stderr_sha256: str
    returncode: int
    stdout: bytes
    stderr: bytes


def _parse_result_record(
        result: ProcessResult, controller_source_commit: str,
        source_receipt_sha256: str, baseline: Mapping[str, Any]) \
        -> Mapping[str, Any]:
    if (COMMIT.fullmatch(controller_source_commit) is None or
            SHA256.fullmatch(source_receipt_sha256) is None):
        fail("Stage 2 source identity is malformed")
    baseline = _baseline_projection(baseline)
    coordinate = result.invocation.coordinate
    if (baseline["K"] != coordinate.K or
            baseline["attempt"] != coordinate.attempt or
            baseline["root_index"] != coordinate.root_index or
            baseline["loss_root"] != coordinate.root or
            baseline["schedule_index"] != coordinate.schedule_index or
            baseline["schedule"] != coordinate.schedule or
            baseline["cell_ordinal"] != coordinate.cell_ordinal):
        fail("pair01 baseline does not match the Stage 2 coordinate")
    try:
        lines = stage1._framed_lines(result, EXPECTED_STAGE0_SOURCE_COMMIT)
    except stage1.Stage1Error as exc:
        fail(str(exc))
    record: Dict[str, Any] = {
        "schema": RECORD_SCHEMA,
        **result.invocation.identity(),
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
        "bench_binary_sha256": EXPECTED_STAGE0_BENCH_SHA256,
        "bench_source_git_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
        "controller_source_git_commit": controller_source_commit,
        "source_receipt_sha256": source_receipt_sha256,
        "baseline_projection_sha256": sha256_json(baseline),
    }
    if result.returncode == 0:
        if result.stderr or len(lines) != 3:
            fail("successful precodefail process framing changed")
        try:
            values = stage1._parse_csv_line(lines[2], len(CSV_HEADER))
            terminal = stage1._parse_terminal_row(
                dict(zip(CSV_HEADER, values)), result.invocation)
        except stage1.Stage1Error as exc:
            fail(str(exc))
    elif result.returncode == 2:
        expected_stderr = (
            "precodefail configuration failed N={} bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt={} packet_attempt={} result=2\n".format(
                coordinate.K, coordinate.attempt,
                coordinate.packet_attempt)).encode("ascii")
        if len(lines) != 2 or result.stderr != expected_stderr:
            fail("precodefail configuration failure is noncanonical")
        terminal = stage1._empty_terminal()
        terminal["configuration_stderr_sha256"] = result.stderr_sha256
    else:
        fail("precodefail process exited with an inadmissible status")
    record.update(terminal)
    expected_packet_seed = stage0._packet_seed_for_offset(
        baseline["effective_packet_seed"], coordinate.attempt,
        result.invocation.delta)
    if record["configuration_failure"] is False:
        if record["effective_precode_seed"] != \
                baseline["effective_precode_seed"]:
            fail("candidate precode seed differs from frozen p")
        if record["effective_packet_seed"] != expected_packet_seed:
            fail("candidate packet seed differs from arbitrary-delta q")
    record["deterministic_projection_sha256"] = sha256_json(
        deterministic_projection(record))
    _exact_fields(record, RECORD_FIELDS, "Stage 2 parsed record")
    return record


def _canonical_b64(value: bytes) -> str:
    return base64.b64encode(value).decode("ascii")


def _raw_record(result: ProcessResult,
                parsed: Mapping[str, Any]) -> Mapping[str, Any]:
    return {
        "schema": RAW_SCHEMA,
        "ordinal": result.invocation.ordinal,
        "invocation": result.invocation.identity(),
        "command_sha256": result.command_sha256,
        "stdout_sha256": result.stdout_sha256,
        "stderr_sha256": result.stderr_sha256,
        "returncode": result.returncode,
        "stdout_base64": _canonical_b64(result.stdout),
        "stderr_base64": _canonical_b64(result.stderr),
        "parsed_record_sha256": sha256_json(parsed),
    }


def parse_process_result(
        result: ProcessResult, controller_source_commit: str,
        source_receipt_sha256: str, baseline: Mapping[str, Any]) \
        -> Tuple[Mapping[str, Any], Mapping[str, Any]]:
    if (result.stdout_sha256 != sha256_bytes(result.stdout) or
            result.stderr_sha256 != sha256_bytes(result.stderr) or
            len(result.stdout) > MAX_STDOUT_BYTES or
            len(result.stderr) > MAX_STDERR_BYTES):
        fail("precodefail byte receipt or bound is inconsistent")
    record = _parse_result_record(
        result, controller_source_commit, source_receipt_sha256, baseline)
    raw = _raw_record(result, record)
    _exact_fields(raw, RAW_FIELDS, "Stage 2 raw record")
    return record, raw


def _decode_raw(value: Any, maximum: int, context: str) -> bytes:
    if type(value) is not str or len(value) > 4 * maximum + 16:
        fail("{} base64 text exceeds its bound".format(context))
    try:
        decoded = base64.b64decode(value.encode("ascii"), validate=True)
    except (UnicodeEncodeError, binascii.Error) as exc:
        fail("{} is not canonical base64: {}".format(context, exc))
    if _canonical_b64(decoded) != value:
        fail("{} is not canonical base64".format(context))
    return decoded


def _invocation_from_record(record: Mapping[str, Any]) -> Invocation:
    return Invocation(record["ordinal"], record["phase"], record["delta"],
                      record["coordinate_ordinal"])


def validate_raw_reparse(
        raw: Mapping[str, Any], record: Mapping[str, Any],
        controller_source_commit: str, source_receipt_sha256: str,
        baseline: Mapping[str, Any]) -> int:
    _exact_fields(raw, RAW_FIELDS, "raw reparse record")
    invocation = _invocation_from_record(record)
    if (raw["schema"] != RAW_SCHEMA or raw["ordinal"] != record["ordinal"] or
            raw["invocation"] != invocation.identity() or
            raw["command_sha256"] != record["command_sha256"] or
            raw["stdout_sha256"] != record["stdout_sha256"] or
            raw["stderr_sha256"] != record["stderr_sha256"] or
            raw["returncode"] != record["returncode"] or
            raw["parsed_record_sha256"] != sha256_json(record)):
        fail("raw stream does not map exactly to its parsed record")
    stdout = _decode_raw(raw["stdout_base64"], MAX_STDOUT_BYTES, "stdout")
    stderr = _decode_raw(raw["stderr_base64"], MAX_STDERR_BYTES, "stderr")
    if (len(stdout) > MAX_STDOUT_BYTES or len(stderr) > MAX_STDERR_BYTES or
            sha256_bytes(stdout) != raw["stdout_sha256"] or
            sha256_bytes(stderr) != raw["stderr_sha256"]):
        fail("raw stream byte hashes or bounds changed")
    result = ProcessResult(
        invocation, raw["command_sha256"], raw["stdout_sha256"],
        raw["stderr_sha256"], raw["returncode"], stdout, stderr)
    reparsed = _parse_result_record(
        result, controller_source_commit, source_receipt_sha256, baseline)
    if reparsed != record:
        fail("raw stdout/stderr reparsing differs from parsed record")
    return len(stdout) + len(stderr)


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


def _append_invocation(
        invocation: Invocation, pinned: screen.PinnedExecutable,
        controller_source_commit: str, source_receipt_sha256: str,
        historical: HistoricalEvidence,
        records: List[Mapping[str, Any]], raw_records: List[Mapping[str, Any]],
        decoded_size: List[int]) -> None:
    result = _run_raw(invocation, pinned)
    baseline = historical.baseline_by_coordinate[
        invocation.coordinate_ordinal]
    record, raw = parse_process_result(
        result, controller_source_commit, source_receipt_sha256, baseline)
    if record["configuration_failure"]:
        fail("fixed precode p failed before packet offset q was used")
    decoded_size[0] += len(result.stdout) + len(result.stderr)
    if decoded_size[0] > MAX_RAW_DECODED_BYTES:
        fail("aggregate raw child output exceeds its bound")
    records.append(record)
    raw_records.append(raw)


def _select_lowest_zero_weak(
        screen_records: Sequence[Mapping[str, Any]]) -> Optional[int]:
    if len(screen_records) != SCREEN_OBSERVATION_COUNT:
        fail("selection received an incomplete survivor screen")
    by_delta: Dict[int, List[Mapping[str, Any]]] = {
        delta: [] for delta in SCREEN_DELTAS}
    for index, record in enumerate(screen_records):
        delta_index, coordinate_index = divmod(
            index, len(INTRODUCTION_COORDINATES))
        expected = Invocation(
            CONTROL_OBSERVATION_COUNT + index, PHASE_SCREEN,
            SCREEN_DELTAS[delta_index],
            INTRODUCTION_COORDINATES[coordinate_index])
        if (record.get("delta") != expected.delta or
                record.get("coordinate_ordinal") !=
                    expected.coordinate_ordinal or
                record.get("phase") != PHASE_SCREEN):
            fail("selection screen order changed")
        by_delta[expected.delta].append(record)
    passing = [delta for delta in SCREEN_DELTAS
               if len(by_delta[delta]) == 44 and
               all(record.get("success") is True
                   for record in by_delta[delta])]
    return passing[0] if passing else None


def _execute_stage2(
        pinned: screen.PinnedExecutable, controller_source_commit: str,
        source_receipt_sha256: str, historical: HistoricalEvidence) \
        -> Tuple[Tuple[Mapping[str, Any], ...],
                 Tuple[Mapping[str, Any], ...], int]:
    records: List[Mapping[str, Any]] = []
    raw_records: List[Mapping[str, Any]] = []
    decoded_size = [0]
    for invocation in control_invocations():
        _append_invocation(
            invocation, pinned, controller_source_commit,
            source_receipt_sha256, historical, records, raw_records,
            decoded_size)
        expected_pair = historical.stage1_control_by_coordinate[
            invocation.coordinate_ordinal]
        observed = scientific_projection(records[-1])
        if any(observed != expected for expected in expected_pair):
            fail("delta2 failure control differs from sealed Stage1")
    if not all(record["rank_fail"] for record in records):
        fail("delta2 failure control did not reproduce all 44 failures")

    for invocation in screen_invocations():
        _append_invocation(
            invocation, pinned, controller_source_commit,
            source_receipt_sha256, historical, records, raw_records,
            decoded_size)
    selected = _select_lowest_zero_weak(
        records[CONTROL_OBSERVATION_COUNT:])
    if selected is not None:
        for invocation in confirmation_invocations(selected):
            _append_invocation(
                invocation, pinned, controller_source_commit,
                source_receipt_sha256, historical, records, raw_records,
                decoded_size)
        # A confirmation mismatch is a canonical scientific REJECT, not an
        # invalid prerequisite.  Adjudication records it and never falls
        # through to a higher delta.
    return tuple(records), tuple(raw_records), decoded_size[0]


def _validate_record(
        record: Mapping[str, Any], invocation: Invocation,
        controller_source_commit: str, source_receipt_sha256: str,
        historical: HistoricalEvidence) -> None:
    _exact_fields(record, RECORD_FIELDS,
                  "Stage 2 record {}".format(invocation.ordinal))
    baseline = historical.baseline_by_coordinate[
        invocation.coordinate_ordinal]
    if (record["schema"] != RECORD_SCHEMA or
            any(record.get(field) != value
                for field, value in invocation.identity().items()) or
            record["block_bytes"] != BLOCK_BYTES or
            record["loss_ppm"] != LOSS_PPM or record["overhead"] != OVERHEAD or
            record["dense_anchor_layout"] != "two07" or
            record["mix_count"] != 2 or
            record["binary_dense_rows"] != BINARY_DENSE_ROWS or
            record["gf256_heavy_rows"] != GF256_HEAVY_ROWS or
            record["heavy_family"] != "periodic" or
            record["construction_seed_basis"] != "production-profile" or
            record["full_payload_solve"] is not True or
            record["overhead_stream"] != "paired" or
            record["cold_solve_wide_xor"] != "policy" or
            record["promotion_evidence"] is not False or
            record["fresh_roots_used"] is not False or
            record["timing_evidence"] is not False or
            record["bench_binary_sha256"] != EXPECTED_STAGE0_BENCH_SHA256 or
            record["bench_source_git_commit"] != EXPECTED_STAGE0_SOURCE_COMMIT or
            record["controller_source_git_commit"] !=
                controller_source_commit or
            record["source_receipt_sha256"] != source_receipt_sha256 or
            record["baseline_projection_sha256"] != sha256_json(baseline) or
            record["deterministic_projection_sha256"] !=
                sha256_json(deterministic_projection(record))):
        fail("Stage 2 record identity or deterministic receipt changed")
    _terminal_record_is_valid(record, "Stage 2 record")
    if record["configuration_failure"]:
        fail("Stage 2 record contains a configuration failure")
    expected_packet = stage0._packet_seed_for_offset(
        baseline["effective_packet_seed"], invocation.coordinate.attempt,
        invocation.delta)
    if (record["effective_precode_seed"] !=
            baseline["effective_precode_seed"] or
            record["effective_packet_seed"] != expected_packet):
        fail("Stage 2 record effective seed differs from frozen p/q")
    for field in ("command_sha256", "stdout_sha256", "stderr_sha256"):
        _sha(record[field], "Stage 2 {}".format(field))


def _ratio(numerator: int, denominator: int) -> str:
    return stage1._ratio_text(numerator, denominator)


def adjudicate(records: Sequence[Mapping[str, Any]],
               historical: HistoricalEvidence) -> Mapping[str, Any]:
    if len(records) not in (
            PRE_CONFIRMATION_OBSERVATION_COUNT, MAX_OBSERVATION_COUNT):
        fail("Stage 2 record roster has the wrong cardinality")
    pre = tuple(control_invocations()) + tuple(screen_invocations())
    if len(records) < len(pre):
        fail("Stage 2 record roster is incomplete")
    controller_commit = records[0].get("controller_source_git_commit")
    source_sha = records[0].get("source_receipt_sha256")
    if (type(controller_commit) is not str or
            COMMIT.fullmatch(controller_commit) is None or
            type(source_sha) is not str or SHA256.fullmatch(source_sha) is None):
        fail("Stage 2 record source identity is malformed")
    for invocation, record in zip(pre, records[:len(pre)]):
        _validate_record(record, invocation, controller_commit, source_sha,
                         historical)

    controls = records[:CONTROL_OBSERVATION_COUNT]
    for record in controls:
        expected_pair = historical.stage1_control_by_coordinate[
            record["coordinate_ordinal"]]
        observed = scientific_projection(record)
        if any(observed != expected for expected in expected_pair):
            fail("delta2 failure control projection mismatch")
    if not all(record["rank_fail"] for record in controls):
        fail("delta2 failure control did not reproduce 44 rank failures")

    screen_rows = records[
        CONTROL_OBSERVATION_COUNT:PRE_CONFIRMATION_OBSERVATION_COUNT]
    selected = _select_lowest_zero_weak(screen_rows)
    expected_count = (MAX_OBSERVATION_COUNT if selected is not None else
                      PRE_CONFIRMATION_OBSERVATION_COUNT)
    if len(records) != expected_count:
        fail("selected confirmation presence differs from selection")

    screen_results = []
    zero_weak = []
    offset = 0
    for delta in SCREEN_DELTAS:
        rows = screen_rows[offset:offset + len(INTRODUCTION_COORDINATES)]
        offset += len(rows)
        successes = sum(int(row["success"]) for row in rows)
        rank_failures = sum(int(row["rank_fail"]) for row in rows)
        errors = sum(int(row["error"]) for row in rows)
        weak = [row["coordinate_ordinal"] for row in rows if row["weak"]]
        passed = successes == 44 and not weak and errors == 0
        if passed:
            zero_weak.append(delta)
        screen_results.append({
            "delta": delta,
            "observation_count": len(rows),
            "successes": successes,
            "rank_failures": rank_failures,
            "errors": errors,
            "configuration_failures": 0,
            "weak_coordinate_ordinals": weak,
            "disposition": "PASS" if passed else "REJECT",
        })
    if selected != (zero_weak[0] if zero_weak else None):
        fail("lowest-numeric selection rule changed")

    confirmation = records[PRE_CONFIRMATION_OBSERVATION_COUNT:]
    work = None
    confirmation_weak = []
    confirmation_errors = 0
    confirmation_configurations = 0
    projection_exact = selected is not None
    projection_mismatches: List[int] = []
    gates = {
        "delta2_control_exact_44_failures": True,
        "all_217_survivors_screened_once": True,
        "selected_zero_weak_delta_exists": selected is not None,
        "selected_confirmation_zero_weak": False,
        "selected_confirmation_errors_zero": False,
        "selected_confirmation_configuration_failures_zero": False,
        "selected_confirmation_projection_exact": False,
    }
    if selected is not None:
        invocations = confirmation_invocations(selected)
        screen_selected = {
            row["coordinate_ordinal"]: row for row in screen_rows
            if row["delta"] == selected
        }
        for invocation, record in zip(invocations, confirmation):
            _validate_record(record, invocation, controller_commit, source_sha,
                             historical)
            coordinate = invocation.coordinate_ordinal
            observed = scientific_projection(record)
            if coordinate in INTRODUCTION_COORDINATES:
                expected = scientific_projection(screen_selected[coordinate])
            else:
                expected = historical.stage0_success_by_delta_coordinate[
                    (selected, coordinate)]
            if observed != expected:
                projection_exact = False
                projection_mismatches.append(coordinate)
        confirmation_weak = [row["coordinate_ordinal"]
                             for row in confirmation if row["weak"]]
        confirmation_errors = sum(int(row["error"]) for row in confirmation)
        confirmation_configurations = sum(
            int(row["configuration_failure"]) for row in confirmation)
        baseline_sums = {
            field: sum(historical.baseline_by_coordinate[coordinate][field]
                       for coordinate in CONFIRMATION_COORDINATES)
            for field in ("block_xors", "gf256_muladds",
                          "inactivated_columns")
        }
        candidate_sums = {
            field: sum(row[field] for row in confirmation)
            for field in baseline_sums
        }
        work = {
            "coordinate_count": len(CONFIRMATION_COORDINATES),
            "selection_used_work": False,
            "pair01_baseline": baseline_sums,
            "selected_confirmation": candidate_sums,
            "selected_to_pair01_ratios": {
                field: _ratio(candidate_sums[field], baseline_sums[field])
                for field in baseline_sums
            },
            "selected_minus_pair01": {
                field: candidate_sums[field] - baseline_sums[field]
                for field in baseline_sums
            },
            "block_xor_threshold_gate_applied": False,
            "targeted_work_gate_applied": False,
            "comparison_caveat": (
                "descriptive only: three pair01 baselines are failed solves"),
        }
        gates.update({
            "selected_confirmation_zero_weak": not confirmation_weak,
            "selected_confirmation_errors_zero": confirmation_errors == 0,
            "selected_confirmation_configuration_failures_zero":
                confirmation_configurations == 0,
            "selected_confirmation_projection_exact": projection_exact,
        })
    return {
        "failure_control_delta": CONTROL_DELTA,
        "failure_control_observation_count": len(controls),
        "failure_control_rank_failures": sum(
            int(row["rank_fail"]) for row in controls),
        "screen_delta_count": len(SCREEN_DELTAS),
        "screen_observation_count": len(screen_rows),
        "screen_results": screen_results,
        "zero_weak_deltas": zero_weak,
        "selected_delta": selected,
        "selection_rule": "lowest numeric zero-weak delta; work blind",
        "confirmation_observation_count": len(confirmation),
        "confirmation_weak_coordinate_ordinals": confirmation_weak,
        "confirmation_errors": confirmation_errors,
        "confirmation_configuration_failures": confirmation_configurations,
        "confirmation_projection_mismatch_coordinate_ordinals":
            projection_mismatches,
        "deduplicated_confirmation_work": work,
        "production_split_pq_implemented": False,
        "selection_leakage": {
            "consulted_coordinate_count": 47,
            "untouched_coordinate_count": 1033,
            "next_replay_is_disjoint": False,
        },
        "gates": gates,
        "fixed_offset_family_retired": not all(gates.values()),
        "disposition": "PASS" if all(gates.values()) else "REJECT",
    }


def _source_receipt() -> Mapping[str, Any]:
    for item in SOURCE_PATHS:
        if not (REPO_ROOT / item).is_file():
            fail("required Stage 2 source is missing: {}".format(item))
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
            list(SOURCE_PATHS), cwd=str(REPO_ROOT), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot attest Stage 2 source checkout: {}".format(exc))
    if (head_process.returncode != 0 or head_process.stderr or
            status_process.returncode != 0 or status_process.stderr or
            tracked_process.returncode != 0 or tracked_process.stderr):
        fail("cannot attest Stage 2 source checkout")
    try:
        head = head_process.stdout.decode("ascii").strip()
        tracked = tracked_process.stdout.decode("ascii").split("\0")
    except UnicodeDecodeError:
        fail("Stage 2 source receipt is not ASCII")
    if COMMIT.fullmatch(head) is None or status_process.stdout:
        fail("Stage 2 linked sources are not clean at a full source HEAD")
    if tracked and tracked[-1] == "":
        tracked.pop()
    if len(tracked) != len(SOURCE_PATHS) or set(tracked) != set(SOURCE_PATHS):
        fail("Stage 2 source receipt includes a path not tracked at HEAD")
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
            fail("Stage 2 source {} differs byte-for-byte from HEAD".format(
                item))
    try:
        final_head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot re-attest Stage 2 source HEAD: {}".format(exc))
    if (final_head.returncode != 0 or final_head.stderr or
            final_head.stdout != (head + "\n").encode("ascii")):
        fail("Stage 2 source HEAD changed while the receipt was built")
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


def _write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]],
                 maximum_bytes: Optional[int] = None) -> Tuple[str, int, int]:
    digest = hashlib.sha256()
    count = 0
    total = 0
    try:
        descriptor = os.open(
            str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except OSError as exc:
        fail("cannot create {}: {}".format(path, exc))
    try:
        for row in rows:
            data = (canonical_json(row) + "\n").encode("ascii")
            total += len(data)
            if maximum_bytes is not None and total > maximum_bytes:
                fail("{} exceeds its publication bound".format(path))
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
    return digest.hexdigest(), count, total


def _read_jsonl(path: Path, maximum_line: int) \
        -> Tuple[Tuple[Mapping[str, Any], ...], str, int]:
    digest = hashlib.sha256()
    rows = []
    total = 0
    try:
        with path.open("rb") as stream:
            for ordinal, line in enumerate(stream):
                total += len(line)
                if len(line) > maximum_line or not line.endswith(b"\n"):
                    fail("published JSONL line {} has invalid framing".format(
                        ordinal))
                digest.update(line)
                rows.append(_parse_canonical_json_line(
                    line, "published JSONL row {}".format(ordinal)))
    except OSError as exc:
        fail("cannot audit published JSONL: {}".format(exc))
    return tuple(rows), digest.hexdigest(), total


def _audit_output_streams(
        records_path: Path, raw_path: Path, expected_records_sha: str,
        expected_raw_sha: str, expected_count: int,
        controller_source_commit: str, source_receipt_sha256: str,
        historical: HistoricalEvidence) -> Mapping[str, Any]:
    records, records_sha, _ = _read_jsonl(records_path, 32 * 1024)
    raw_records, raw_sha, raw_size = _read_jsonl(
        raw_path, MAX_RAW_JSONL_LINE_BYTES)
    if (records_sha != expected_records_sha or raw_sha != expected_raw_sha or
            len(records) != expected_count or
            len(raw_records) != expected_count or
            raw_size > MAX_RAW_FILE_BYTES):
        fail("published parsed/raw stream receipt changed")
    decoded = 0
    for record, raw in zip(records, raw_records):
        invocation = _invocation_from_record(record)
        _validate_record(record, invocation, controller_source_commit,
                         source_receipt_sha256, historical)
        decoded += validate_raw_reparse(
            raw, record, controller_source_commit, source_receipt_sha256,
            historical.baseline_by_coordinate[
                invocation.coordinate_ordinal])
        if decoded > MAX_RAW_DECODED_BYTES:
            fail("published decoded raw stream exceeds its bound")
    verdict = adjudicate(records, historical)
    return {"decoded_bytes": decoded, "raw_file_bytes": raw_size,
            "verdict": verdict}


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(str(path), flags)
    except OSError as exc:
        fail("cannot open output directory for fsync: {}".format(exc))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def run_stage2(bench: Path, stage0_dir: Path, stage1_dir: Path,
               output_dir: Path) -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path must be one fresh path below an existing directory")
    # Historical evidence comes first: no current source/bench/workload action
    # can occur unless both immutable prior stages fully authenticate.
    with load_historical_evidence(stage0_dir, stage1_dir) as historical:
        historical_receipt = historical.receipt()
        source = _source_receipt()
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            return _run_stage2_pinned(
                pinned, output_dir, source, contract, historical,
                historical_receipt)


def _run_stage2_pinned(
        pinned: screen.PinnedExecutable, output_dir: Path,
        source: Mapping[str, Any], contract: Mapping[str, Any],
        historical: HistoricalEvidence,
        historical_receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    controller_commit = source["source_git_commit"]
    source_sha = source["source_receipt_sha256"]
    bench_receipt = pinned.receipt()
    if (bench_receipt.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench_receipt.get("size") != EXPECTED_STAGE0_BENCH_SIZE):
        fail("Stage 2 requires the exact authenticated Stage0 benchmark")
    # Final pre-work action.
    description = screen.read_bench_description(
        pinned.path, EXPECTED_STAGE0_SOURCE_COMMIT, pinned.descriptor)
    records, raw_records, decoded_size = _execute_stage2(
        pinned, controller_commit, source_sha, historical)
    verdict = adjudicate(records, historical)

    temporary = Path(tempfile.mkdtemp(
        prefix=".wh2-mix2-packet-offset-stage2-",
        dir=str(output_dir.parent)))
    published = False
    try:
        records_sha, record_count, _ = _write_jsonl(
            temporary / RECORD_NAME, records)
        raw_sha, raw_count, raw_file_size = _write_jsonl(
            temporary / RAW_NAME, raw_records, MAX_RAW_FILE_BYTES)
        if record_count != len(records) or raw_count != record_count:
            fail("Stage 2 parsed/raw publication cardinality changed")
        audit = _audit_output_streams(
            temporary / RECORD_NAME, temporary / RAW_NAME, records_sha,
            raw_sha, record_count, controller_commit, source_sha, historical)
        if (audit["decoded_bytes"] != decoded_size or
                audit["raw_file_bytes"] != raw_file_size or
                audit["verdict"] != verdict):
            fail("Stage 2 durable raw reparse audit changed the verdict")
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                historical.receipt() != historical_receipt):
            fail("Stage 2 executable, source, or history changed during run")
        summary: Dict[str, Any] = {
            "schema": SUMMARY_SCHEMA,
            "contract_sha256": contract["contract_sha256"],
            "source_receipt": source,
            "bench": bench_receipt,
            "bench_description": description,
            "historical_evidence": historical_receipt,
            "records_file": RECORD_NAME,
            "records_sha256": records_sha,
            "record_count": record_count,
            "raw_records_file": RAW_NAME,
            "raw_records_sha256": raw_sha,
            "raw_record_count": raw_count,
            "raw_decoded_bytes": decoded_size,
            "raw_file_bytes": raw_file_size,
            "raw_encoding": "canonical-base64-exact-bytes",
            "raw_reparse_equality": True,
            **verdict,
            "stage2_only": True,
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
                historical.receipt() != historical_receipt):
            fail("Stage 2 inputs changed before publication")
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
    subparsers.add_parser("describe", help="print the frozen Stage 2 contract")
    run_parser = subparsers.add_parser(
        "run", help="run the bounded remaining-global-offset filter")
    run_parser.add_argument("--bench", type=Path, required=True)
    run_parser.add_argument("--stage0-dir", type=Path, required=True)
    run_parser.add_argument("--stage1-dir", type=Path, required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        _validate_constants()
        if arguments.command == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run_stage2(
            arguments.bench, arguments.stage0_dir, arguments.stage1_dir,
            arguments.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (Stage2Error, stage1.Stage1Error, stage0.Stage0Error,
            screen.ShortScreenError) as exc:
        print("wh2 MIX2 packet-offset Stage 2: {}".format(exc),
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
