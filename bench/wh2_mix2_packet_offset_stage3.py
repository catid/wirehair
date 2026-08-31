#!/usr/bin/env python3
"""Sealed 1,080-coordinate consumed-root replay for MIX2 delta 39.

Stage 3 authenticates the complete pair01 baseline and every durable Stage-2
stdout/stderr byte before executing anything.  It then runs fixed packet
offset 39 twice on all 1,080 already-consumed coordinates.  The 47 coordinates
consulted while selecting the offset and the untouched 1,033-coordinate
complement are reported separately.  This is selection-leaked consumed-root
evidence only: it uses no fresh roots or timing and cannot promote or modify
the production codec.
"""

from __future__ import annotations

import argparse
import base64
import binascii
from contextlib import ExitStack
from dataclasses import dataclass
from decimal import Decimal, localcontext
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
from typing import (Any, Dict, Iterable, List, Mapping, Optional, Sequence,
                    Tuple)

try:
    from bench import wh2_mix2_packet_offset_stage0 as stage0
    from bench import wh2_mix2_packet_offset_stage1 as stage1
    from bench import wh2_mix2_packet_offset_stage2 as stage2
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_packet_offset_stage0 as stage0
    import wh2_mix2_packet_offset_stage1 as stage1
    import wh2_mix2_packet_offset_stage2 as stage2
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage3-contract.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage3-record.v1"
RAW_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage3-raw.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage3-summary.v1"
RECORD_NAME = "mix2-packet-offset-stage3-records.jsonl"
RAW_NAME = "mix2-packet-offset-stage3-raw.jsonl"
SUMMARY_NAME = "mix2-packet-offset-stage3-summary.json"

CANDIDATE_PAIR = "01"
PACKET_DELTA = 39
REPETITIONS = (0, 1)
PACKET_ATTEMPT_STRIDE = 0x9e3779b9
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
BINARY_DENSE_ROWS = 12
GF256_HEAVY_ROWS = 12
ZERO_SHA256 = "0" * 64

K_VALUES = tuple(stage1.K_VALUES)
V9_ATTEMPTS = tuple(stage1.V9_ATTEMPTS)
ROOTS = tuple(stage1.ROOTS)
SCHEDULES = tuple(stage1.SCHEDULES)
COORDINATE_COUNT = 1080
OBSERVATION_COUNT = 2160

# Literal full-coordinate Stage-2 selection footprint.  Do not reconstruct it
# from weak/failure status: its identity and complement are frozen evidence.
CONSULTED_COORDINATES = (
    1, 2, 5, 8, 10, 14, 16, 20, 23, 26, 28, 30, 31, 32, 35, 38, 44,
    50, 56, 59, 109, 112, 115, 118, 121, 127, 129, 130, 134, 136, 140,
    142, 146, 152, 161, 167, 175, 177, 181, 218, 224, 573, 751, 820,
    851, 894, 1021,
)
UNTOUCHED_COORDINATES = tuple(
    value for value in range(COORDINATE_COUNT)
    if value not in frozenset(CONSULTED_COORDINATES))
EXPECTED_CONSULTED_SHA256 = (
    "1c71a0ec8bf9721cd017e787d93eaaf762e4877096454346d7bc29f2b89bd553")
EXPECTED_UNTOUCHED_SHA256 = (
    "8d968ad4c43498265b79bfb433f122cd262731955a4ac35d9a5d6e040f4b3ccd")
EXPECTED_FULL_ORDINALS_SHA256 = (
    "e599f73bec6735d9785f466c71fb8012f7cd8f228821f080b496fbdab90ef485")

EXPECTED_BASELINE_SUMS = dict(stage1.EXPECTED_BASELINE_SUMS)
EXPECTED_CONSULTED_BASELINE_SUMS = {
    "block_xors": 1819147,
    "gf256_muladds": 41842,
    "inactivated_columns": 1875,
}
EXPECTED_UNTOUCHED_BASELINE_SUMS = {
    "block_xors": 184527320,
    "gf256_muladds": 3637374,
    "inactivated_columns": 133200,
}
BASELINE_WEAK_COORDINATES = tuple(stage2.ORIGINAL_WEAK_COORDINATES)

EXPECTED_CANDIDATE_PROFILE_SHA256 = stage1.EXPECTED_CANDIDATE_PROFILE_SHA256
EXPECTED_STAGE0_BENCH_SHA256 = stage1.EXPECTED_STAGE0_BENCH_SHA256
EXPECTED_STAGE0_BENCH_SIZE = stage1.EXPECTED_STAGE0_BENCH_SIZE
EXPECTED_STAGE0_SOURCE_COMMIT = stage1.EXPECTED_STAGE0_SOURCE_COMMIT
EXPECTED_STAGE1_HELPER_SHA256 = stage2.EXPECTED_STAGE1_HELPER_SHA256
EXPECTED_STAGE2_HELPER_SHA256 = (
    "53161487ce1bae0c4596e0ab3ee306fcdec878bb5bec11f22c53728eb442af8c")
EXPECTED_STAGE2_CONTRACT_SHA256 = stage2.EXPECTED_CONTRACT_SHA256
EXPECTED_STAGE2_SOURCE_COMMIT = "2a89408e4fb18ce51bdd5007052167bbd5336067"
EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256 = (
    "6c3f760d95fddf533e50fbc27b3352013f0deace6b42836960d7a3bf66f75739")
EXPECTED_STAGE2_SUMMARY_SELF_SHA256 = (
    "7b9ffdbd4ffee32794d687f3e441bfb367a617a328e071152e2d473fa915be27")
STAGE2_FILES = {
    "mix2-packet-offset-stage2-records.jsonl": (
        18644503,
        "beaf82f1e2bf994d0af3bf55b8851e564fbcfe62c99d7e8a3fe82b3befbb7cee"),
    "mix2-packet-offset-stage2-raw.jsonl": (
        27421985,
        "ba273aff7ebc5271db261a5ecea788b67f887106b3de89a4d5766e10ff52a0e3"),
    "mix2-packet-offset-stage2-summary.json": (
        50540,
        "6940c4c0388fc535c0acb5593b58ce288f8bda20bab1cc47acd92dcde39b7c95"),
}
EXPECTED_STAGE2_RECORD_COUNT = 9639
EXPECTED_STAGE2_RAW_DECODED_BYTES = 14955348

# The user explicitly accepts the historical 1.70546696-percent XOR margin.
# Stage 3 preserves it by requiring delta39 to be no worse than pair01; the
# retired five-percent threshold must not be reintroduced here.
ACCEPTED_HISTORICAL_XOR_RATIO_MAX = "0.9829453304"

CSV_HEADER = tuple(stage1.CSV_HEADER)
METADATA_ORDER = tuple(stage1.METADATA_ORDER)
MAX_STDOUT_BYTES = stage2.MAX_STDOUT_BYTES
MAX_STDERR_BYTES = stage2.MAX_STDERR_BYTES
MAX_RAW_DECODED_BYTES = 16 * 1024 * 1024
MAX_RAW_FILE_BYTES = 32 * 1024 * 1024
MAX_JSONL_LINE_BYTES = 32 * 1024
PROCESS_TIMEOUT_SECONDS = 600
CHILD_ENVIRONMENT = {"LANG": "C", "LC_ALL": "C", "TZ": "UTC"}
UINT64_MAX = (1 << 64) - 1

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent
SOURCE_PATHS = (
    "CMakeLists.txt", "codec/CMakeLists.txt",
    "bench/test_wh2_mix2_packet_offset_stage3.py",
    "bench/wh2_mix2_packet_offset_stage3.py",
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
    "WirehairPeelFixups.inc", "bench/wh2_mix2_packet_offset_stage3.py",
    "bench/test_wh2_mix2_packet_offset_stage3.py",
    "bench/wh2_mix2_packet_offset_stage2.py",
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

# Replaced after the literal contract body is finalized.
EXPECTED_CONTRACT_SHA256 = (
    "e99d5a09b51b4d691113295a4ba004445e35f69f020bde1a04c19e4a8fc31216")


class Stage3Error(RuntimeError):
    """Stage 3 cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise Stage3Error(message)


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


def _sha(value: Any, context: str) -> str:
    if type(value) is not str or SHA256.fullmatch(value) is None:
        fail("{} is not a lowercase SHA-256".format(context))
    return value


def _contract_body() -> Mapping[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "bead": "wirehair-sxvz.16.1.20.38.4",
        "scope": "fixed delta39 1080-coordinate consumed-root replay only",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "candidate": {
            "base_profile_sha256": EXPECTED_CANDIDATE_PROFILE_SHA256,
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "mix_pair": CANDIDATE_PAIR,
            "packet_delta": PACKET_DELTA,
            "binary_dense_rows": BINARY_DENSE_ROWS,
            "gf256_heavy_rows": GF256_HEAVY_ROWS,
            "construction_seed_basis": "production-profile",
            "precode_attempt_rule": "exact frozen v9 p(K)",
            "packet_attempt_rule": "q=(p+39) mod 256",
            "packet_seed_rule": (
                "effective_q=(effective_p+(q-p)*0x9e3779b9) mod 2^32"),
            "selection_rule": "Stage2 lowest numeric zero-weak delta",
            "fallback_delta225_forbidden": True,
            "test_hook_only": True,
            "production_split_pq_implemented": False,
            "production_complexity_condition": (
                "later production review must prove one immutable profile "
                "constant and bounded construction-time arithmetic with no "
                "additional offset map, allocation, retry, new row, per-packet "
                "branch, hot-path work, or asymptotic change"),
            "preexisting_whole_profile_footprint": (
                "the authenticated v9 architecture already requires its "
                "63,999-byte K-indexed repair map and other whole-profile "
                "state; those costs are not caused by packet offset 39 and "
                "must remain visible in any whole-profile promotion review"),
            "offset_only_complexity_scope": (
                "delta39 may add only its immutable constant and bounded "
                "configuration-time arithmetic; it may not add another map "
                "or any online search"),
        },
        "domain": {
            "K": list(K_VALUES),
            "v9_attempts": list(V9_ATTEMPTS),
            "consumed_roots": list(ROOTS),
            "unconsumed_final_roots_excluded":
                list(screen.FINAL_VALIDATION_ROOTS),
            "schedules": list(SCHEDULES),
            "coordinate_order": "K-major, root-major, schedule-major",
            "coordinate_count": COORDINATE_COUNT,
            "repetitions": list(REPETITIONS),
            "invocation_order": "coordinate-major, repetition-minor",
            "observation_count": OBSERVATION_COUNT,
            "consulted_coordinate_ordinals": list(CONSULTED_COORDINATES),
            "consulted_ordinals_sha256": EXPECTED_CONSULTED_SHA256,
            "consulted_coordinate_count": len(CONSULTED_COORDINATES),
            "untouched_coordinate_ordinals": list(UNTOUCHED_COORDINATES),
            "untouched_ordinals_sha256": EXPECTED_UNTOUCHED_SHA256,
            "untouched_coordinate_count": len(UNTOUCHED_COORDINATES),
            "full_ordinals_sha256": EXPECTED_FULL_ORDINALS_SHA256,
            "selection_leakage": True,
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "trials": 1,
            "threads": 1,
            "full_payload_solve": True,
            "timing_used": False,
        },
        "prerequisites": {
            "stage1_full_history": {
                "loader": "sealed Stage1 full baseline/V9/Stage0 loader",
                "baseline_files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in stage1.BASELINE_FILES.items()
                },
                "v9_files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in stage1.V9_FILES.items()
                },
                "stage0_files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in stage1.STAGE0_FILES.items()
                },
                "selected_record_count": COORDINATE_COUNT,
                "baseline_sums": dict(EXPECTED_BASELINE_SUMS),
                "baseline_weak_coordinates":
                    list(BASELINE_WEAK_COORDINATES),
            },
            "stage2": {
                "contract_sha256": EXPECTED_STAGE2_CONTRACT_SHA256,
                "source_commit": EXPECTED_STAGE2_SOURCE_COMMIT,
                "source_receipt_sha256":
                    EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256,
                "helper_sha256": EXPECTED_STAGE2_HELPER_SHA256,
                "summary_self_sha256":
                    EXPECTED_STAGE2_SUMMARY_SELF_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in STAGE2_FILES.items()
                },
                "record_count": EXPECTED_STAGE2_RECORD_COUNT,
                "raw_record_count": EXPECTED_STAGE2_RECORD_COUNT,
                "raw_decoded_bytes": EXPECTED_STAGE2_RAW_DECODED_BYTES,
                "required_disposition": "PASS",
                "required_selected_delta": PACKET_DELTA,
                "required_zero_weak_deltas": [39, 225],
                "authentication": (
                    "descriptor-pin all files; hash the exact summary, parsed, "
                    "and raw bytes consumed under unchanged fd identity; "
                    "canonical-parse both streams; reparse every stdout/stderr "
                    "byte; independently adjudicate; and compare the entire "
                    "published verdict"),
                "baseline_overlap_rule": (
                    "all 47 Stage2 baseline projections equal the independently "
                    "loaded full Stage1 pair01 baseline"),
            },
        },
        "gate": {
            "recovery": {
                "candidate_weak_coordinates": 0,
                "errors": 0,
                "configuration_failures": 0,
                "introductions_vs_pair01": 0,
                "repair_all_three_pair01_weak_coordinates": True,
                "repeat_scientific_projection_exact": True,
                "repeat_work_exact": True,
                "consulted_47_zero_weak": True,
                "untouched_1033_zero_weak": True,
                "both_repetitions_of_consulted_47_equal_stage2_confirmation":
                    True,
            },
            "work": {
                "deduplication": "repetition zero after exact repeat equality",
                "full_1080_block_xors": "candidate <= pair01",
                "full_1080_gf256_muladds": "candidate <= pair01",
                "full_1080_inactivated_columns": "candidate <= pair01",
                "accepted_historical_xor_ratio_max":
                    ACCEPTED_HISTORICAL_XOR_RATIO_MAX,
                "accepted_historical_xor_gain_percent": "1.70546696",
                "old_five_percent_threshold_applied": False,
                "consulted_untouched_and_per_K_work": "reported, not gated",
            },
            "overall": "PASS iff every frozen recovery and work gate passes",
        },
        "transport": {
            "historical_preflight": (
                "all Stage1/Stage2 history authenticates before current source, "
                "benchmark description, or workload; pre/post and retained "
                "device/inode/size/mtime_ns/ctime_ns snapshots guard every "
                "imported loader input"),
            "benchmark": {
                "sha256": EXPECTED_STAGE0_BENCH_SHA256,
                "size": EXPECTED_STAGE0_BENCH_SIZE,
                "source_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
            },
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
                "canonical fsynced parsed/raw JSONL and summary; full raw "
                "reparse and independent adjudication; atomic no-replace "
                "directory publication; parent fsync"),
        },
        "stop_rule": (
            "PASS licenses only a separately frozen disjoint fresh-root "
            "holdout. REJECT retires the fixed-offset family without trying "
            "delta225. No timing, all-K campaign, production change, profile "
            "identifier, or promotion is licensed here."),
    }


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("MIX2 packet-offset Stage 3 contract differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _helper_hash(path: Path) -> str:
    return sha256_bytes(screen._stable_file_bytes(path, 4 * 1024 * 1024))


def _validate_constants() -> None:
    if (CANDIDATE_PAIR != "01" or PACKET_DELTA != 39 or
            REPETITIONS != (0, 1) or PACKET_ATTEMPT_STRIDE != 0x9e3779b9 or
            COORDINATE_COUNT != 1080 or OBSERVATION_COUNT != 2160 or
            K_VALUES != tuple(stage1.K_VALUES) or
            V9_ATTEMPTS != tuple(stage1.V9_ATTEMPTS) or
            ROOTS != tuple(stage1.ROOTS) or
            SCHEDULES != tuple(stage1.SCHEDULES)):
        fail("Stage 3 candidate or full consumed roster changed")
    if (len(CONSULTED_COORDINATES) != 47 or
            tuple(sorted(set(CONSULTED_COORDINATES))) !=
                CONSULTED_COORDINATES or
            CONSULTED_COORDINATES != tuple(stage2.CONFIRMATION_COORDINATES) or
            len(UNTOUCHED_COORDINATES) != 1033 or
            set(CONSULTED_COORDINATES) & set(UNTOUCHED_COORDINATES) or
            tuple(sorted(CONSULTED_COORDINATES + UNTOUCHED_COORDINATES)) !=
                tuple(range(COORDINATE_COUNT)) or
            sha256_json(list(CONSULTED_COORDINATES)) !=
                EXPECTED_CONSULTED_SHA256 or
            sha256_json(list(UNTOUCHED_COORDINATES)) !=
                EXPECTED_UNTOUCHED_SHA256 or
            sha256_json(list(range(COORDINATE_COUNT))) !=
                EXPECTED_FULL_ORDINALS_SHA256):
        fail("Stage 3 consulted/untouched partition changed")
    if (EXPECTED_BASELINE_SUMS != {
            "block_xors": 186346467, "gf256_muladds": 3679216,
            "inactivated_columns": 135075} or
            any(EXPECTED_CONSULTED_BASELINE_SUMS[field] +
                EXPECTED_UNTOUCHED_BASELINE_SUMS[field] !=
                EXPECTED_BASELINE_SUMS[field]
                for field in EXPECTED_BASELINE_SUMS)):
        fail("Stage 3 frozen pair01 work receipt changed")
    if set(ROOTS) & set(screen.FINAL_VALIDATION_ROOTS):
        fail("Stage 3 includes a fresh root")
    if (tuple(CSV_HEADER) != tuple(stage2.CSV_HEADER) or
            screen.candidate_profile_sha256() !=
                EXPECTED_CANDIDATE_PROFILE_SHA256 or
            stage0._packet_seed_for_offset("0x12345678", 255, 39) !=
                "0xf52e28a7"):
        fail("Stage 3 benchmark protocol or delta39 seed rule changed")
    stage2._validate_constants()
    if stage2.contract_description()["contract_sha256"] != \
            EXPECTED_STAGE2_CONTRACT_SHA256:
        fail("Stage 2 imported contract changed")
    for path, expected in (
            (Path(stage1.__file__).resolve(), EXPECTED_STAGE1_HELPER_SHA256),
            (Path(stage2.__file__).resolve(), EXPECTED_STAGE2_HELPER_SHA256)):
        if _helper_hash(path) != expected:
            fail("imported helper differs from frozen source: {}".format(path))


@dataclass(frozen=True)
class Coordinate:
    ordinal: int

    def __post_init__(self) -> None:
        if type(self.ordinal) is not int or not 0 <= self.ordinal < 1080:
            fail("Stage 3 coordinate ordinal is out of range")

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
    def packet_attempt(self) -> int:
        return (self.attempt + PACKET_DELTA) & 255

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
    def baseline_record_ordinal(self) -> int:
        return self.base.baseline_record_ordinal

    @property
    def stratum(self) -> str:
        return "consulted" if self.ordinal in CONSULTED_COORDINATES \
            else "untouched"

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
            "baseline_record_ordinal": self.baseline_record_ordinal,
            "stratum": self.stratum,
        }


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    coordinate_ordinal: int
    repetition: int

    def __post_init__(self) -> None:
        if (type(self.ordinal) is not int or
                type(self.coordinate_ordinal) is not int or
                type(self.repetition) is not int or
                not 0 <= self.coordinate_ordinal < COORDINATE_COUNT or
                self.repetition not in REPETITIONS or
                self.ordinal != self.coordinate_ordinal * 2 + self.repetition):
            fail("Stage 3 invocation differs from canonical order")

    @property
    def coordinate(self) -> Coordinate:
        return Coordinate(self.coordinate_ordinal)

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
            "repetition": self.repetition,
            "mix_pair": CANDIDATE_PAIR,
            "delta": PACKET_DELTA,
            "precode_attempt": self.coordinate.attempt,
            "packet_attempt": self.coordinate.packet_attempt,
            **self.coordinate.identity(),
        }


def make_invocations() -> Tuple[Invocation, ...]:
    values = tuple(
        Invocation(coordinate * 2 + repetition, coordinate, repetition)
        for coordinate in range(COORDINATE_COUNT)
        for repetition in REPETITIONS)
    if len(values) != OBSERVATION_COUNT or len({
            (value.coordinate_ordinal, value.repetition)
            for value in values}) != OBSERVATION_COUNT:
        fail("Stage 3 invocation roster is incomplete or duplicated")
    return values


@dataclass
class HistoricalEvidence:
    stack: Optional[ExitStack]
    stage1_history: stage1.HistoricalEvidence
    stage2_history: stage2.HistoricalEvidence
    stage2_pins: Mapping[str, stage1.PinnedFile]
    baseline: Sequence[Mapping[str, Any]]
    stage2_confirmation_by_coordinate: Mapping[int, Mapping[str, Any]]
    stage2_receipt: Mapping[str, Any]
    historical_snapshots: Mapping[str, Tuple[int, int, int, int, int]]

    def receipt(self) -> Mapping[str, Any]:
        _assert_snapshots(self.historical_snapshots)
        return {
            "stage1_full_history": self.stage1_history.receipt(),
            "stage2_recovery_history": self.stage2_history.receipt(),
            "stage2": {
                **self.stage2_receipt,
                "files": {
                    name: pin.receipt()
                    for name, pin in sorted(self.stage2_pins.items())
                },
            },
            "stable_path_snapshots": {
                path: {
                    "device": identity[0], "inode": identity[1],
                    "size": identity[2], "mtime_ns": identity[3],
                    "ctime_ns": identity[4],
                }
                for path, identity in sorted(
                    self.historical_snapshots.items())
            },
        }

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


def _snapshot_group(
        directory: Path, specs: Mapping[str, Tuple[int, str]], label: str) \
        -> Mapping[str, Tuple[int, int, int, int, int]]:
    try:
        root = stage1._validate_historical_directory(directory, label)
    except stage1.Stage1Error as exc:
        fail(str(exc))
    snapshots = {}
    for name, (expected_size, _) in specs.items():
        path = root / name
        try:
            state = path.lstat()
        except OSError as exc:
            fail("cannot snapshot {} {}: {}".format(label, name, exc))
        if stat.S_ISLNK(state.st_mode) or not stat.S_ISREG(state.st_mode):
            fail("{} {} must be a non-symlink regular file".format(
                label, name))
        if state.st_size != expected_size:
            fail("{} {} changed size before parsing".format(label, name))
        snapshots[str(path.resolve(strict=True))] = stage1._fd_identity(state)
    return snapshots


def _merge_snapshots(
        destination: Dict[str, Tuple[int, int, int, int, int]],
        values: Mapping[str, Tuple[int, int, int, int, int]]) -> None:
    for path, identity in values.items():
        if path in destination and destination[path] != identity:
            fail("historical path identity changed between prerequisite loads")
        destination[path] = identity


def _assert_snapshots(
        snapshots: Mapping[str, Tuple[int, int, int, int, int]]) -> None:
    for text, expected in snapshots.items():
        path = Path(text)
        try:
            state = path.lstat()
        except OSError as exc:
            fail("cannot re-snapshot historical input {}: {}".format(
                path, exc))
        if (stat.S_ISLNK(state.st_mode) or not stat.S_ISREG(state.st_mode) or
                stage1._fd_identity(state) != expected):
            fail("historical input changed while its parsed state was used: {}"
                 .format(path))


def _read_pinned_bytes_exact(
        pin: stage1.PinnedFile, expected_size: int, expected_sha: str,
        maximum: int, context: str) -> bytes:
    try:
        before = stage1._fd_identity(os.fstat(pin.descriptor))
    except OSError as exc:
        fail("cannot snapshot {} descriptor: {}".format(context, exc))
    try:
        data = pin.read_bytes(maximum)
    except stage1.Stage1Error as exc:
        fail(str(exc))
    try:
        after = stage1._fd_identity(os.fstat(pin.descriptor))
    except OSError as exc:
        fail("cannot re-snapshot {} descriptor: {}".format(context, exc))
    if (before != after or len(data) != expected_size or
            sha256_bytes(data) != expected_sha):
        fail("{} exact consumed-byte receipt changed".format(context))
    return data


def _read_pinned_jsonl_exact(
        pin: stage1.PinnedFile, expected_size: int, expected_sha: str,
        context: str) -> Tuple[Mapping[str, Any], ...]:
    try:
        before = stage1._fd_identity(os.fstat(pin.descriptor))
    except OSError as exc:
        fail("cannot snapshot {} descriptor: {}".format(context, exc))
    digest = hashlib.sha256()
    total = 0
    rows = []
    try:
        iterator = pin.lines(MAX_JSONL_LINE_BYTES)
        for ordinal, line in iterator:
            digest.update(line)
            total += len(line)
            row = stage2._parse_canonical_json_line(
                line, "{} row {}".format(context, ordinal))
            if row.get("ordinal") != ordinal:
                fail("{} ordinal changed".format(context))
            rows.append(row)
    except (stage1.Stage1Error, stage2.Stage2Error) as exc:
        fail(str(exc))
    try:
        after = stage1._fd_identity(os.fstat(pin.descriptor))
    except OSError as exc:
        fail("cannot re-snapshot {} descriptor: {}".format(context, exc))
    if (before != after or total != expected_size or
            digest.hexdigest() != expected_sha):
        fail("{} exact consumed-byte stream changed".format(context))
    return tuple(rows)


def _load_stage2_artifact(
        pins: Mapping[str, stage1.PinnedFile],
        historical: stage2.HistoricalEvidence) \
        -> Tuple[Mapping[int, Mapping[str, Any]], Mapping[str, Any]]:
    summary_bytes = _read_pinned_bytes_exact(
        pins[stage2.SUMMARY_NAME], STAGE2_FILES[stage2.SUMMARY_NAME][0],
        STAGE2_FILES[stage2.SUMMARY_NAME][1], 1024 * 1024,
        "Stage2 summary")
    summary = stage2._parse_canonical_json_bytes(
        summary_bytes, "Stage2 summary")
    if (summary.get("schema") != stage2.SUMMARY_SCHEMA or
            summary.get("contract_sha256") !=
                EXPECTED_STAGE2_CONTRACT_SHA256 or
            summary.get("summary_sha256") !=
                EXPECTED_STAGE2_SUMMARY_SELF_SHA256 or
            stage2.self_hash(summary, "summary_sha256") !=
                EXPECTED_STAGE2_SUMMARY_SELF_SHA256 or
            summary.get("records_file") != stage2.RECORD_NAME or
            summary.get("records_sha256") !=
                STAGE2_FILES[stage2.RECORD_NAME][1] or
            summary.get("record_count") != EXPECTED_STAGE2_RECORD_COUNT or
            summary.get("raw_records_file") != stage2.RAW_NAME or
            summary.get("raw_records_sha256") !=
                STAGE2_FILES[stage2.RAW_NAME][1] or
            summary.get("raw_record_count") != EXPECTED_STAGE2_RECORD_COUNT or
            summary.get("raw_decoded_bytes") !=
                EXPECTED_STAGE2_RAW_DECODED_BYTES or
            summary.get("raw_file_bytes") !=
                STAGE2_FILES[stage2.RAW_NAME][0] or
            summary.get("raw_encoding") != "canonical-base64-exact-bytes" or
            summary.get("raw_reparse_equality") is not True or
            summary.get("stage2_only") is not True or
            summary.get("promotion_evidence") is not False or
            summary.get("fresh_roots_used") is not False or
            summary.get("timing_evidence") is not False):
        fail("Stage2 summary identity, stream receipt, or scope changed")
    source = summary.get("source_receipt")
    bench = summary.get("bench")
    description = summary.get("bench_description")
    if (type(source) is not dict or
            source.get("source_git_commit") != EXPECTED_STAGE2_SOURCE_COMMIT or
            source.get("source_receipt_sha256") !=
                EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256 or
            stage2.self_hash(source, "source_receipt_sha256") !=
                EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256 or
            type(bench) is not dict or
            bench.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench.get("size") != EXPECTED_STAGE0_BENCH_SIZE or
            type(description) is not dict or
            description.get("schema") != stage2.BENCH_DESCRIPTION_SCHEMA or
            description.get("source_git_commit") !=
                EXPECTED_STAGE0_SOURCE_COMMIT or
            summary.get("historical_evidence") != historical.receipt()):
        fail("Stage2 source, benchmark, or historical receipt changed")

    records = _read_pinned_jsonl_exact(
        pins[stage2.RECORD_NAME], STAGE2_FILES[stage2.RECORD_NAME][0],
        STAGE2_FILES[stage2.RECORD_NAME][1], "Stage2 parsed stream")
    raw_records = _read_pinned_jsonl_exact(
        pins[stage2.RAW_NAME], STAGE2_FILES[stage2.RAW_NAME][0],
        STAGE2_FILES[stage2.RAW_NAME][1], "Stage2 raw stream")
    if (len(records) != EXPECTED_STAGE2_RECORD_COUNT or
            len(raw_records) != EXPECTED_STAGE2_RECORD_COUNT):
        fail("Stage2 parsed/raw roster changed")
    decoded = 0
    for record, raw in zip(records, raw_records):
        invocation = stage2._invocation_from_record(record)
        stage2._validate_record(
            record, invocation, EXPECTED_STAGE2_SOURCE_COMMIT,
            EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256, historical)
        decoded += stage2.validate_raw_reparse(
            raw, record, EXPECTED_STAGE2_SOURCE_COMMIT,
            EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256,
            historical.baseline_by_coordinate[
                invocation.coordinate_ordinal])
        if decoded > stage2.MAX_RAW_DECODED_BYTES:
            fail("Stage2 decoded raw bytes exceed their frozen bound")
    if decoded != EXPECTED_STAGE2_RAW_DECODED_BYTES:
        fail("Stage2 decoded raw-byte receipt changed")
    verdict = stage2.adjudicate(records, historical)
    if any(summary.get(field) != value for field, value in verdict.items()):
        fail("Stage2 independent adjudication differs from its summary")
    if (verdict.get("disposition") != "PASS" or
            verdict.get("selected_delta") != PACKET_DELTA or
            verdict.get("zero_weak_deltas") != [39, 225] or
            not all(verdict.get("gates", {}).values())):
        fail("Stage2 did not admit fixed delta39 exactly")
    confirmation = records[stage2.PRE_CONFIRMATION_OBSERVATION_COUNT:]
    by_coordinate = {
        row["coordinate_ordinal"]: stage2.scientific_projection(row)
        for row in confirmation
    }
    if (len(confirmation) != 47 or
            tuple(sorted(by_coordinate)) != CONSULTED_COORDINATES or
            any(row.get("delta") != PACKET_DELTA or
                row.get("phase") != stage2.PHASE_CONFIRMATION or
                row.get("success") is not True for row in confirmation)):
        fail("Stage2 selected-confirmation projection changed")
    receipt = {
        "contract_sha256": EXPECTED_STAGE2_CONTRACT_SHA256,
        "source_git_commit": EXPECTED_STAGE2_SOURCE_COMMIT,
        "source_receipt_sha256": EXPECTED_STAGE2_SOURCE_RECEIPT_SHA256,
        "summary_self_sha256": EXPECTED_STAGE2_SUMMARY_SELF_SHA256,
        "records_sha256": STAGE2_FILES[stage2.RECORD_NAME][1],
        "raw_records_sha256": STAGE2_FILES[stage2.RAW_NAME][1],
        "record_count": len(records),
        "raw_decoded_bytes": decoded,
        "selected_delta": PACKET_DELTA,
        "zero_weak_deltas": [39, 225],
        "disposition": "PASS",
    }
    return by_coordinate, receipt


def load_historical_evidence(
        baseline_dir: Path, v9_dir: Path, stage0_dir: Path,
        stage1_dir: Path, stage2_dir: Path) -> HistoricalEvidence:
    stack = ExitStack()
    try:
        snapshots: Dict[str, Tuple[int, int, int, int, int]] = {}
        full_before = {}
        for directory, specs, label in (
                (baseline_dir, stage1.BASELINE_FILES, "pair01 baseline"),
                (v9_dir, stage1.V9_FILES, "V9"),
                (stage0_dir, stage1.STAGE0_FILES, "Stage0")):
            values = _snapshot_group(directory, specs, label)
            _merge_snapshots(full_before, values)
            _merge_snapshots(snapshots, values)
        full = stack.enter_context(stage1.load_historical_evidence(
            baseline_dir, v9_dir, stage0_dir))
        _assert_snapshots(full_before)
        recovery_before = {}
        for directory, specs, label in (
                (stage0_dir, stage2.STAGE0_FILES, "Stage0"),
                (stage1_dir, stage2.STAGE1_FILES, "Stage1")):
            values = _snapshot_group(directory, specs, label)
            _merge_snapshots(recovery_before, values)
            _merge_snapshots(snapshots, values)
        recovery = stack.enter_context(stage2.load_historical_evidence(
            stage0_dir, stage1_dir))
        _assert_snapshots(recovery_before)
        stage2_before = _snapshot_group(stage2_dir, STAGE2_FILES, "Stage2")
        _merge_snapshots(snapshots, stage2_before)
        pins = _pin_group(stack, stage2_dir, STAGE2_FILES, "Stage2")
        confirmation, stage2_receipt = _load_stage2_artifact(pins, recovery)
        _assert_snapshots(stage2_before)
        if len(full.baseline) != COORDINATE_COUNT:
            fail("Stage1 full pair01 baseline has changed cardinality")
        for coordinate in CONSULTED_COORDINATES:
            expected = stage1._baseline_projection(full.baseline[coordinate])
            observed = recovery.baseline_by_coordinate[coordinate]
            if expected != observed:
                fail("Stage1/Stage2 baseline overlap differs at {}".format(
                    coordinate))
        result = HistoricalEvidence(
            stack, full, recovery, pins, full.baseline, confirmation,
            stage2_receipt, snapshots)
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


SCIENTIFIC_FIELDS = tuple(stage2.SCIENTIFIC_FIELDS)
DETERMINISTIC_FIELDS = (
    "stratum", *SCIENTIFIC_FIELDS, "attempt", "schedule_index",
    "cell_ordinal", "baseline_record_ordinal", "command_sha256",
    "bench_binary_sha256", "bench_source_git_commit",
    "controller_source_git_commit", "source_receipt_sha256",
    "pair01_baseline_record_sha256", "pair01_baseline",
    "baseline_matrix_sha256",
)
RECORD_FIELDS = frozenset((
    "schema", *DETERMINISTIC_FIELDS, "ordinal", "repetition",
    "promotion_evidence", "fresh_roots_used", "timing_evidence",
    "stdout_sha256", "stderr_sha256", "deterministic_projection_sha256",
))
RAW_FIELDS = frozenset((
    "schema", "ordinal", "invocation", "command_sha256", "stdout_sha256",
    "stderr_sha256", "returncode", "stdout_base64", "stderr_base64",
    "parsed_record_sha256",
))


def scientific_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    try:
        return stage2.scientific_projection(record)
    except stage2.Stage2Error as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def deterministic_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in DETERMINISTIC_FIELDS if field not in record]
    if missing:
        fail("Stage 3 record omits deterministic fields: {}".format(missing))
    return {field: record[field] for field in DETERMINISTIC_FIELDS}


def _baseline_projection(row: Mapping[str, Any]) -> Mapping[str, Any]:
    try:
        return stage1._baseline_projection(row)
    except stage1.Stage1Error as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _parse_result_record(
        result: ProcessResult, controller_source_commit: str,
        source_receipt_sha256: str, baseline: Mapping[str, Any]) \
        -> Mapping[str, Any]:
    if (COMMIT.fullmatch(controller_source_commit) is None or
            SHA256.fullmatch(source_receipt_sha256) is None or
            len(result.stdout) > MAX_STDOUT_BYTES or
            len(result.stderr) > MAX_STDERR_BYTES or
            result.stdout_sha256 != sha256_bytes(result.stdout) or
            result.stderr_sha256 != sha256_bytes(result.stderr)):
        fail("Stage 3 process or source receipt is malformed")
    coordinate = result.invocation.coordinate
    if (baseline.get("record_ordinal") !=
            coordinate.baseline_record_ordinal or
            baseline.get("K") != coordinate.K or
            baseline.get("attempt") != coordinate.attempt or
            baseline.get("root_index") != coordinate.root_index or
            baseline.get("loss_root") != coordinate.root or
            baseline.get("schedule_index") != coordinate.schedule_index or
            baseline.get("schedule") != coordinate.schedule):
        fail("pair01 baseline row does not match Stage 3 coordinate")
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
        "pair01_baseline_record_sha256": sha256_json(baseline),
        "pair01_baseline": _baseline_projection(baseline),
        "baseline_matrix_sha256": stage1.BASELINE_FILES[
            "attempt-crossfit-matrix.jsonl"][1],
    }
    if result.returncode == 0:
        if result.stderr or len(lines) != 3:
            fail("successful precodefail framing changed")
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
        try:
            terminal = stage1._empty_terminal()
        except stage1.Stage1Error as exc:
            fail(str(exc))
        terminal["configuration_stderr_sha256"] = result.stderr_sha256
    else:
        fail("precodefail exited with an inadmissible status")
    record.update(terminal)
    if not record["configuration_failure"]:
        expected_packet = stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], coordinate.attempt,
            PACKET_DELTA)
        if (record["effective_precode_seed"] !=
                baseline["effective_precode_seed"] or
                record["effective_packet_seed"] != expected_packet):
            fail("candidate effective seed differs from frozen delta39 p/q")
    record["deterministic_projection_sha256"] = sha256_json(
        deterministic_projection(record))
    _exact_fields(record, RECORD_FIELDS, "Stage 3 parsed record")
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
    record = _parse_result_record(
        result, controller_source_commit, source_receipt_sha256, baseline)
    raw = _raw_record(result, record)
    _exact_fields(raw, RAW_FIELDS, "Stage 3 raw record")
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


def validate_raw_reparse(
        raw: Mapping[str, Any], record: Mapping[str, Any],
        controller_source_commit: str, source_receipt_sha256: str,
        baseline: Mapping[str, Any]) -> int:
    _exact_fields(raw, RAW_FIELDS, "Stage 3 raw reparse record")
    invocation = Invocation(
        record["ordinal"], record["coordinate_ordinal"],
        record["repetition"])
    if (raw["schema"] != RAW_SCHEMA or raw["ordinal"] != record["ordinal"] or
            raw["invocation"] != invocation.identity() or
            raw["command_sha256"] != record["command_sha256"] or
            raw["stdout_sha256"] != record["stdout_sha256"] or
            raw["stderr_sha256"] != record["stderr_sha256"] or
            raw["returncode"] != record["returncode"] or
            raw["parsed_record_sha256"] != sha256_json(record)):
        fail("Stage 3 raw stream does not map exactly to parsed record")
    stdout = _decode_raw(raw["stdout_base64"], MAX_STDOUT_BYTES, "stdout")
    stderr = _decode_raw(raw["stderr_base64"], MAX_STDERR_BYTES, "stderr")
    if (sha256_bytes(stdout) != raw["stdout_sha256"] or
            sha256_bytes(stderr) != raw["stderr_sha256"]):
        fail("Stage 3 raw stream byte hash changed")
    reparsed = _parse_result_record(
        ProcessResult(
            invocation, raw["command_sha256"], raw["stdout_sha256"],
            raw["stderr_sha256"], raw["returncode"], stdout, stderr),
        controller_source_commit, source_receipt_sha256, baseline)
    if reparsed != record:
        fail("Stage 3 raw stdout/stderr reparse differs from parsed record")
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
        invocation, sha256_json({"argv": argv[1:]}),
        sha256_bytes(completed.stdout), sha256_bytes(completed.stderr),
        completed.returncode, completed.stdout, completed.stderr)


def _execute_roster(
        pinned: screen.PinnedExecutable, controller_source_commit: str,
        source_receipt_sha256: str, historical: HistoricalEvidence) \
        -> Tuple[Tuple[Mapping[str, Any], ...],
                 Tuple[Mapping[str, Any], ...], int]:
    records = []
    raw_records = []
    decoded = 0
    for invocation in make_invocations():
        result = _run_raw(invocation, pinned)
        record, raw = parse_process_result(
            result, controller_source_commit, source_receipt_sha256,
            historical.baseline[invocation.coordinate_ordinal])
        if record["configuration_failure"]:
            fail("fixed precode p failed before delta39 packet q was used")
        decoded += len(result.stdout) + len(result.stderr)
        if decoded > MAX_RAW_DECODED_BYTES:
            fail("Stage 3 aggregate raw child output exceeds its bound")
        records.append(record)
        raw_records.append(raw)
    return tuple(records), tuple(raw_records), decoded


def _validate_record(
        record: Mapping[str, Any], invocation: Invocation,
        controller_source_commit: str, source_receipt_sha256: str,
        baseline: Mapping[str, Any]) -> None:
    _exact_fields(record, RECORD_FIELDS,
                  "Stage 3 record {}".format(invocation.ordinal))
    projection = _baseline_projection(baseline)
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
            record["bench_source_git_commit"] !=
                EXPECTED_STAGE0_SOURCE_COMMIT or
            record["controller_source_git_commit"] !=
                controller_source_commit or
            record["source_receipt_sha256"] != source_receipt_sha256 or
            record["pair01_baseline_record_sha256"] != sha256_json(baseline) or
            record["pair01_baseline"] != projection or
            record["baseline_matrix_sha256"] != stage1.BASELINE_FILES[
                "attempt-crossfit-matrix.jsonl"][1] or
            record["deterministic_projection_sha256"] !=
                sha256_json(deterministic_projection(record))):
        fail("Stage 3 record identity or deterministic receipt changed")
    try:
        stage2._terminal_record_is_valid(record, "Stage 3 record")
    except stage2.Stage2Error as exc:
        fail(str(exc))
    if not record["configuration_failure"]:
        expected_packet = stage0._packet_seed_for_offset(
            baseline["effective_packet_seed"], invocation.coordinate.attempt,
            PACKET_DELTA)
        if (record["effective_precode_seed"] !=
                baseline["effective_precode_seed"] or
                record["effective_packet_seed"] != expected_packet):
            fail("Stage 3 effective seed differs from fixed p/q")
    for field in ("command_sha256", "stdout_sha256", "stderr_sha256",
                  "pair01_baseline_record_sha256"):
        _sha(record[field], "Stage 3 {}".format(field))


def _ratio_text(numerator: int, denominator: int) -> str:
    if type(numerator) is not int or type(denominator) is not int or \
            numerator < 0 or denominator <= 0:
        fail("aggregate ratio has an invalid operand")
    with localcontext() as context:
        context.prec = 60
        return format(Decimal(numerator) / Decimal(denominator), ".12f")


WORK_FIELDS = ("block_xors", "gf256_muladds", "inactivated_columns")


def _work_report(
        coordinates: Sequence[int], deduplicated: Sequence[Mapping[str, Any]],
        baseline: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    baseline_sums = {
        field: sum(baseline[index][field] for index in coordinates)
        for field in WORK_FIELDS
    }
    candidate_sums = {
        field: (sum(deduplicated[index][field] for index in coordinates)
                if all(type(deduplicated[index][field]) is int
                       for index in coordinates) else None)
        for field in WORK_FIELDS
    }
    return {
        "coordinate_count": len(coordinates),
        "pair01_baseline": baseline_sums,
        "delta39_candidate": candidate_sums,
        "candidate_minus_pair01": {
            field: (candidate_sums[field] - baseline_sums[field]
                    if type(candidate_sums[field]) is int else None)
            for field in WORK_FIELDS
        },
        "candidate_to_pair01_ratios": {
            field: (_ratio_text(candidate_sums[field], baseline_sums[field])
                    if type(candidate_sums[field]) is int else None)
            for field in WORK_FIELDS
        },
    }


def adjudicate(records: Sequence[Mapping[str, Any]],
               historical: HistoricalEvidence) -> Mapping[str, Any]:
    if (len(records) != OBSERVATION_COUNT or
            len(historical.baseline) != COORDINATE_COUNT):
        fail("Stage 3 result or baseline roster has wrong cardinality")
    controller_commit = records[0].get("controller_source_git_commit")
    source_sha = records[0].get("source_receipt_sha256")
    if (type(controller_commit) is not str or
            COMMIT.fullmatch(controller_commit) is None or
            type(source_sha) is not str or SHA256.fullmatch(source_sha) is None):
        fail("Stage 3 result source receipt is malformed")
    by_key = {}
    for invocation, record in zip(make_invocations(), records):
        _validate_record(
            record, invocation, controller_commit, source_sha,
            historical.baseline[invocation.coordinate_ordinal])
        key = (record["coordinate_ordinal"], record["repetition"])
        if key in by_key:
            fail("Stage 3 result roster is duplicated")
        by_key[key] = record
    if len(by_key) != OBSERVATION_COUNT:
        fail("Stage 3 result roster is incomplete")

    deduplicated = []
    repeat_exact = 0
    repeat_work_exact = 0
    for coordinate in range(COORDINATE_COUNT):
        first = by_key[(coordinate, 0)]
        second = by_key[(coordinate, 1)]
        exact = scientific_projection(first) == scientific_projection(second)
        work_exact = all(
            type(first[field]) is int and
            first[field] == second[field] for field in WORK_FIELDS)
        repeat_exact += int(exact)
        repeat_work_exact += int(work_exact)
        deduplicated.append(first)

    candidate_weak = []
    baseline_weak = []
    repairs = []
    introductions = []
    for coordinate in range(COORDINATE_COUNT):
        pair = (by_key[(coordinate, 0)], by_key[(coordinate, 1)])
        candidate_success = all(row["success"] for row in pair)
        baseline_success = historical.baseline[coordinate]["success"]
        receipt = {
            **Coordinate(coordinate).identity(),
            "candidate_repetition_outcomes": [row["outcome"] for row in pair],
            "candidate_success_all_repetitions": candidate_success,
            "baseline_outcome": historical.baseline[coordinate]["outcome"],
            "baseline_success": baseline_success,
        }
        if not candidate_success:
            candidate_weak.append(receipt)
        if not baseline_success:
            baseline_weak.append(receipt)
        if not baseline_success and candidate_success:
            repairs.append(receipt)
        if baseline_success and not candidate_success:
            introductions.append(receipt)

    errors = sum(int(row["error"]) for row in records)
    configurations = sum(int(row["configuration_failure"]) for row in records)
    projection_mismatches = []
    projection_exact_observations = 0
    for coordinate in CONSULTED_COORDINATES:
        expected = historical.stage2_confirmation_by_coordinate[coordinate]
        exact = True
        for repetition in REPETITIONS:
            observed = scientific_projection(by_key[(coordinate, repetition)])
            if observed == expected:
                projection_exact_observations += 1
            else:
                exact = False
        if not exact:
            projection_mismatches.append(coordinate)

    def stratum_report(name: str, coordinates: Sequence[int]) \
            -> Mapping[str, Any]:
        weak = [coordinate for coordinate in coordinates
                if not all(by_key[(coordinate, repetition)]["success"]
                           for repetition in REPETITIONS)]
        baseline_weak_values = [coordinate for coordinate in coordinates
                                if not historical.baseline[coordinate]["success"]]
        introduced = [coordinate for coordinate in weak
                      if historical.baseline[coordinate]["success"]]
        repaired = [coordinate for coordinate in baseline_weak_values
                    if coordinate not in weak]
        observations = [by_key[(coordinate, repetition)]
                        for coordinate in coordinates
                        for repetition in REPETITIONS]
        return {
            "name": name,
            "coordinate_count": len(coordinates),
            "observation_count": len(observations),
            "weak_coordinate_ordinals": weak,
            "errors": sum(int(row["error"]) for row in observations),
            "configuration_failures": sum(
                int(row["configuration_failure"]) for row in observations),
            "baseline_weak_coordinate_ordinals": baseline_weak_values,
            "repaired_coordinate_ordinals": repaired,
            "introduced_weak_coordinate_ordinals": introduced,
            "work": _work_report(
                coordinates, deduplicated, historical.baseline),
        }

    consulted = stratum_report("consulted", CONSULTED_COORDINATES)
    untouched = stratum_report("untouched", UNTOUCHED_COORDINATES)
    full_work = _work_report(
        tuple(range(COORDINATE_COUNT)), deduplicated, historical.baseline)
    if (full_work["pair01_baseline"] != EXPECTED_BASELINE_SUMS or
            consulted["work"]["pair01_baseline"] !=
                EXPECTED_CONSULTED_BASELINE_SUMS or
            untouched["work"]["pair01_baseline"] !=
                EXPECTED_UNTOUCHED_BASELINE_SUMS or
            tuple(item["coordinate_ordinal"] for item in baseline_weak) !=
                BASELINE_WEAK_COORDINATES):
        fail("Stage 3 pair01 baseline aggregate or weak roster changed")

    per_K = []
    for K_index, K in enumerate(K_VALUES):
        coordinates = tuple(range(K_index * 36, K_index * 36 + 36))
        per_K.append({
            "K": K,
            "coordinate_count": 36,
            "candidate_observation_count": 72,
            "work": _work_report(
                coordinates, deduplicated, historical.baseline),
        })
    candidate_sums = full_work["delta39_candidate"]
    gates = {
        "candidate_weak_coordinates_zero": not candidate_weak,
        "candidate_errors_zero": errors == 0,
        "candidate_configuration_failures_zero": configurations == 0,
        "introductions_vs_pair01_zero": not introductions,
        "all_three_pair01_weak_coordinates_repaired":
            len(repairs) == len(BASELINE_WEAK_COORDINATES),
        "full_repeat_scientific_projection_exact":
            repeat_exact == COORDINATE_COUNT,
        "full_repeat_work_exact": repeat_work_exact == COORDINATE_COUNT,
        "consulted_47_zero_weak": not consulted["weak_coordinate_ordinals"],
        "untouched_1033_zero_weak": not untouched["weak_coordinate_ordinals"],
        "consulted_94_observations_equal_stage2_confirmation":
            (projection_exact_observations == 94 and
             not projection_mismatches),
        "aggregate_block_xors_not_above_pair01":
            (type(candidate_sums["block_xors"]) is int and
             candidate_sums["block_xors"] <=
                EXPECTED_BASELINE_SUMS["block_xors"]),
        "aggregate_gf256_muladds_not_above_pair01":
            (type(candidate_sums["gf256_muladds"]) is int and
             candidate_sums["gf256_muladds"] <=
                EXPECTED_BASELINE_SUMS["gf256_muladds"]),
        "aggregate_inactivated_columns_not_above_pair01":
            (type(candidate_sums["inactivated_columns"]) is int and
             candidate_sums["inactivated_columns"] <=
                EXPECTED_BASELINE_SUMS["inactivated_columns"]),
    }
    return {
        "candidate_pair": CANDIDATE_PAIR,
        "candidate_packet_delta": PACKET_DELTA,
        "candidate_observation_count": len(records),
        "deduplicated_candidate_coordinate_count": len(deduplicated),
        "repeat_exact_coordinate_count": repeat_exact,
        "repeat_work_exact_coordinate_count": repeat_work_exact,
        "candidate_errors": errors,
        "candidate_configuration_failures": configurations,
        "baseline_weak_coordinates": baseline_weak,
        "candidate_weak_coordinates": candidate_weak,
        "repaired_coordinates": repairs,
        "introduced_weak_coordinates": introductions,
        "stage2_confirmation_projection_exact_observation_count":
            projection_exact_observations,
        "stage2_confirmation_projection_mismatch_coordinate_ordinals":
            projection_mismatches,
        "selection_leakage": {
            "consulted_coordinate_count": 47,
            "untouched_coordinate_count": 1033,
            "next_holdout_must_be_disjoint": True,
        },
        "strata": {"consulted": consulted, "untouched": untouched},
        "aggregates": {
            **full_work,
            "deduplication": "repetition zero after exact repeat equality",
            "accepted_historical_xor_ratio_max":
                ACCEPTED_HISTORICAL_XOR_RATIO_MAX,
            "old_five_percent_threshold_applied": False,
        },
        "per_K_work": per_K,
        "production_split_pq_implemented": False,
        "promotion_evidence": False,
        "gates": gates,
        "fixed_offset_family_retired": not all(gates.values()),
        "disposition": "PASS" if all(gates.values()) else "REJECT",
    }


def _source_receipt() -> Mapping[str, Any]:
    for item in SOURCE_PATHS:
        if not (REPO_ROOT / item).is_file():
            fail("required Stage 3 source is missing: {}".format(item))
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
        fail("cannot attest Stage 3 source checkout: {}".format(exc))
    if (head_process.returncode != 0 or head_process.stderr or
            status_process.returncode != 0 or status_process.stderr or
            tracked_process.returncode != 0 or tracked_process.stderr):
        fail("cannot attest Stage 3 source checkout")
    try:
        head = head_process.stdout.decode("ascii").strip()
        tracked = tracked_process.stdout.decode("ascii").split("\0")
    except UnicodeDecodeError:
        fail("Stage 3 source receipt is not ASCII")
    if COMMIT.fullmatch(head) is None or status_process.stdout:
        fail("Stage 3 linked sources are not clean at full source HEAD")
    if tracked and tracked[-1] == "":
        tracked.pop()
    if len(tracked) != len(SOURCE_PATHS) or set(tracked) != set(SOURCE_PATHS):
        fail("Stage 3 source receipt includes a path not tracked at HEAD")
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
            fail("Stage 3 source {} differs from HEAD".format(item))
    try:
        final_head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot re-attest Stage 3 source HEAD: {}".format(exc))
    if (final_head.returncode != 0 or final_head.stderr or
            final_head.stdout != (head + "\n").encode("ascii")):
        fail("Stage 3 source HEAD changed while receipt was built")
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
                fail("{} exceeds publication bound".format(path))
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


def _read_jsonl(path: Path) \
        -> Tuple[Tuple[Mapping[str, Any], ...], str, int]:
    digest = hashlib.sha256()
    rows = []
    total = 0
    try:
        with path.open("rb") as stream:
            for ordinal, line in enumerate(stream):
                total += len(line)
                if len(line) > MAX_JSONL_LINE_BYTES or not line.endswith(b"\n"):
                    fail("published JSONL row {} has invalid framing".format(
                        ordinal))
                digest.update(line)
                try:
                    rows.append(stage2._parse_canonical_json_line(
                        line, "published row {}".format(ordinal)))
                except stage2.Stage2Error as exc:
                    fail(str(exc))
    except OSError as exc:
        fail("cannot audit published JSONL: {}".format(exc))
    return tuple(rows), digest.hexdigest(), total


def _audit_output_streams(
        records_path: Path, raw_path: Path, expected_records_sha: str,
        expected_raw_sha: str, expected_count: int,
        controller_source_commit: str, source_receipt_sha256: str,
        historical: HistoricalEvidence) -> Mapping[str, Any]:
    records, records_sha, _ = _read_jsonl(records_path)
    raw_records, raw_sha, raw_size = _read_jsonl(raw_path)
    if (records_sha != expected_records_sha or raw_sha != expected_raw_sha or
            len(records) != expected_count or
            len(raw_records) != expected_count or
            raw_size > MAX_RAW_FILE_BYTES):
        fail("published Stage 3 parsed/raw stream receipt changed")
    decoded = 0
    for record, raw in zip(records, raw_records):
        invocation = Invocation(
            record["ordinal"], record["coordinate_ordinal"],
            record["repetition"])
        baseline = historical.baseline[invocation.coordinate_ordinal]
        _validate_record(
            record, invocation, controller_source_commit,
            source_receipt_sha256, baseline)
        decoded += validate_raw_reparse(
            raw, record, controller_source_commit, source_receipt_sha256,
            baseline)
        if decoded > MAX_RAW_DECODED_BYTES:
            fail("published decoded raw stream exceeds Stage 3 bound")
    return {
        "decoded_bytes": decoded,
        "raw_file_bytes": raw_size,
        "verdict": adjudicate(records, historical),
    }


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


def run_stage3(
        bench: Path, baseline_dir: Path, v9_dir: Path, stage0_dir: Path,
        stage1_dir: Path, stage2_dir: Path, output_dir: Path) \
        -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path must be one fresh path below existing directory")
    # Historical authentication deliberately precedes current source, binary
    # description, and every scientific child invocation.
    with load_historical_evidence(
            baseline_dir, v9_dir, stage0_dir, stage1_dir,
            stage2_dir) as historical:
        historical_receipt = historical.receipt()
        source = _source_receipt()
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            return _run_stage3_pinned(
                pinned, output_dir, source, contract, historical,
                historical_receipt)


def _run_stage3_pinned(
        pinned: screen.PinnedExecutable, output_dir: Path,
        source: Mapping[str, Any], contract: Mapping[str, Any],
        historical: HistoricalEvidence,
        historical_receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    controller_commit = source["source_git_commit"]
    source_sha = source["source_receipt_sha256"]
    bench_receipt = pinned.receipt()
    if (bench_receipt.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench_receipt.get("size") != EXPECTED_STAGE0_BENCH_SIZE):
        fail("Stage 3 requires exact authenticated Stage0 benchmark")
    description = screen.read_bench_description(
        pinned.path, EXPECTED_STAGE0_SOURCE_COMMIT, pinned.descriptor)
    records, raw_records, decoded = _execute_roster(
        pinned, controller_commit, source_sha, historical)
    verdict = adjudicate(records, historical)

    temporary = Path(tempfile.mkdtemp(
        prefix=".wh2-mix2-packet-offset-stage3-",
        dir=str(output_dir.parent)))
    published = False
    try:
        records_sha, record_count, _ = _write_jsonl(
            temporary / RECORD_NAME, records)
        raw_sha, raw_count, raw_size = _write_jsonl(
            temporary / RAW_NAME, raw_records, MAX_RAW_FILE_BYTES)
        if record_count != OBSERVATION_COUNT or raw_count != record_count:
            fail("Stage 3 publication cardinality changed")
        audit = _audit_output_streams(
            temporary / RECORD_NAME, temporary / RAW_NAME, records_sha,
            raw_sha, record_count, controller_commit, source_sha, historical)
        if (audit["decoded_bytes"] != decoded or
                audit["raw_file_bytes"] != raw_size or
                audit["verdict"] != verdict):
            fail("Stage 3 durable raw reparse changed verdict")
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                historical.receipt() != historical_receipt):
            fail("Stage 3 executable, source, or history changed during run")
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
            "raw_decoded_bytes": decoded,
            "raw_file_bytes": raw_size,
            "raw_encoding": "canonical-base64-exact-bytes",
            "raw_reparse_equality": True,
            **verdict,
            "stage3_only": True,
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
            fail("Stage 3 inputs changed before publication")
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
    subparsers.add_parser("describe", help="print frozen Stage 3 contract")
    run_parser = subparsers.add_parser(
        "run", help="run fixed delta39 consumed-root replay")
    run_parser.add_argument("--bench", type=Path, required=True)
    run_parser.add_argument("--baseline-dir", type=Path, required=True)
    run_parser.add_argument("--v9-dir", type=Path, required=True)
    run_parser.add_argument("--stage0-dir", type=Path, required=True)
    run_parser.add_argument("--stage1-dir", type=Path, required=True)
    run_parser.add_argument("--stage2-dir", type=Path, required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        _validate_constants()
        if arguments.command == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run_stage3(
            arguments.bench, arguments.baseline_dir, arguments.v9_dir,
            arguments.stage0_dir, arguments.stage1_dir,
            arguments.stage2_dir, arguments.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (Stage3Error, stage2.Stage2Error, stage1.Stage1Error,
            stage0.Stage0Error, screen.ShortScreenError) as exc:
        print("wh2 MIX2 packet-offset Stage 3: {}".format(exc),
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
