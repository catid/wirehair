#!/usr/bin/env python3
"""Sealed consumed-root Stage 1 replay for MIX2 packet offset delta 2.

The controller validates the complete immutable pair01 attempt matrix and its
V9 and Stage-0 provenance before it starts any candidate process.  It then
runs pair01 with q=(p+2) mod 256 twice on each of the 1,080 already-consumed
coordinates.  Timing is parsed only as transport framing and is never used as
evidence.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import csv
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
from typing import (Any, Dict, Iterable, Iterator, List, Mapping, Optional,
                    Sequence, Tuple)

try:
    from bench import wh2_mix2_attempt_crossfit as crossfit
    from bench import wh2_mix2_packet_offset_stage0 as stage0
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_attempt_crossfit as crossfit
    import wh2_mix2_packet_offset_stage0 as stage0
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage1-contract.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage1-record.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-packet-offset-stage1-summary.v1"
BENCH_DESCRIPTION_SCHEMA = "wirehair.wh2.v2-bench-description.v1"
RECORD_NAME = "mix2-packet-offset-stage1-records.jsonl"
SUMMARY_NAME = "mix2-packet-offset-stage1-summary.json"

CANDIDATE_PAIR = "01"
PACKET_DELTA = 2
PACKET_ATTEMPT_STRIDE = 0x9e3779b9
REPETITIONS = (0, 1)
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
BINARY_DENSE_ROWS = 12
GF256_HEAVY_ROWS = 12
ZERO_SHA256 = "0" * 64

K_VALUES = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512, 513,
    1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000, 10001, 16384,
    20000, 20001, 32768, 49152, 64000,
)
V9_ATTEMPTS = (
    9, 8, 0, 2, 49, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
)
ATTEMPT_BY_K = dict(zip(K_VALUES, V9_ATTEMPTS))
ROOTS = (
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
SCHEDULES = ("burst", "adversarial", "repair-only")
STAIRCASE_BY_K = dict(crossfit.STAIRCASE_BY_K)

COORDINATE_COUNT = len(K_VALUES) * len(ROOTS) * len(SCHEDULES)
OBSERVATION_COUNT = COORDINATE_COUNT * len(REPETITIONS)
BASELINE_MATRIX_RECORD_COUNT = (
    len(K_VALUES) * 256 * len(ROOTS) * len(SCHEDULES))
V9_RESULT_RECORD_COUNT = len(K_VALUES) * 3 * len(SCHEDULES) * 4
V9_CANDIDATE_OBSERVATION_COUNT = len(K_VALUES) * 3 * len(SCHEDULES) * 2

EXPECTED_CANDIDATE_PROFILE_SHA256 = (
    "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0")
EXPECTED_STAGE0_CONTRACT_SHA256 = (
    "aa3ce118502c362ad3ef6f232e350b2d7a0e0ecae4dc270f25ee678bd61b7557")
EXPECTED_V9_CONTRACT_SHA256 = (
    "28e6affc80680377df2e8099825b1b475b65e51d04116a8c2fd6465566041cdb")
EXPECTED_CROSSFIT_CONTRACT_SHA256 = (
    "cb15c941b3f303a544945934c6a31181c2109f3ae61f348153515d037dae119c")

EXPECTED_STAGE0_HELPER_SHA256 = (
    "b20c618924ebacd189ff06f0f7c2d5912450825a6967cc1442dc0a46ac82268f")
EXPECTED_V9_HELPER_SHA256 = (
    "b36a256277c215a10a54cac5e81bc0546391f36989abf26ca33016cada730af1")
EXPECTED_CROSSFIT_HELPER_SHA256 = (
    "c94bef05c103cc0c663975de17f0b5556880abe9a0c14e545acf1858e4cfddf2")

EXPECTED_STAGE0_SOURCE_COMMIT = "56e6a2ee4d567a7077a80e50bc6c3fc0705c0a60"
EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 = (
    "9264cab112c9f3d950a72ee07d54beb155a04b9f1f76b3e2c9c07f2e56d8a456")
EXPECTED_STAGE0_BENCH_SHA256 = (
    "9094ae526c811d5697e2324fdaa340f1b8d7b706c092fdced2ff88d50235d8c2")
EXPECTED_STAGE0_BENCH_SIZE = 665712
EXPECTED_STAGE0_V9_RECEIPT_SHA256 = (
    "3e7ec94b80640105033968d723ec5ce0a9f688b4fd29012f2b61f18b085c3c06")
EXPECTED_STAGE0_SURVIVOR_COUNT = 218
EXPECTED_V9_SOURCE_COMMIT = "52177d5d798935e13e83f01a436bd1b7198e80ed"
EXPECTED_CROSSFIT_SOURCE_COMMIT = "c9cfbc09503b2a68f3833fa6e691f524d207ae7c"
EXPECTED_V9_WORKER_SHA256 = (
    "23e2548804fe9cc7f30b2864429f37f48ac6f6e92c306f61533be51c0454481b")

# name -> (exact byte count, exact file SHA-256)
BASELINE_FILES = {
    "attempt-crossfit-input.json": (
        6291,
        "5899eab5d615d0fabc40d60ef8b66f4f491da7a778993bedcd733efc8a8ae2ea"),
    "attempt-crossfit-matrix.jsonl": (
        197150043,
        "91ceef3e4082e6f9854a4db8b15437397f8c50bf22df6ba0403fedf5de43ab2b"),
    "attempt-crossfit-summary.json": (
        2688,
        "4a11280b33a557d58b76af3de41237815649b9cec02dd27d3c9ca138ae380540"),
}
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
STAGE0_FILES = {
    "mix2-packet-offset-stage0-records.jsonl": (
        2559698,
        "362d45c43e6e209ae39bc91a54feedbe12abed96b3b8b884d635a0a23eb92bfe"),
    "mix2-packet-offset-stage0-summary.json": (
        150131,
        "4807530fde16ecf097d9f389bc06a2ec699c2f2f116d66addfa0b0a5212d8cbf"),
}
EXPECTED_BASELINE_INPUT_SELF_SHA256 = (
    "f664cdf7e9ee130415fc7de88f120cb095d1e6c8ebcca6ddb79e5fc321cfa523")
EXPECTED_BASELINE_SUMMARY_SELF_SHA256 = (
    "75146cd0e84f7e19180964384e8042b5805e6f066c526d23bdb2d2f09678f46a")
EXPECTED_V9_INPUT_SELF_SHA256 = (
    "5c7f3309d3a05a458adbcbe1a50ff25e3b70c829ae9806b7592c9d188a1840c4")
EXPECTED_V9_SUMMARY_SELF_SHA256 = (
    "bbbf7e357e9d5cd84b9c79d6d91d0109f83609763452588e264e6cf38a200f4e")
EXPECTED_STAGE0_SUMMARY_SELF_SHA256 = (
    "b4ebb3819ff16621c4a1ed7666ddcfd0f005e5847b3a5cf33721fa36124cc0f5")
EXPECTED_ATTEMPT_MAP_SHA256 = (
    "229c93bf4d47013f96900d7496f82acc1dbc9913ef7995116c23aeaae30d816a")
EXPECTED_SELECTED_RAW_SHA256 = (
    "b210e2f4b319a818c6a71ec71e89e855322792490437fc2c79b5580865306c4f")
EXPECTED_SELECTED_PROJECTION_SHA256 = (
    "384d8e91d6a578b3754e831a7e1fbc4448f50e9fa7c58d100b696c571a678a15")
EXPECTED_BASELINE_SUMS = {
    "block_xors": 186346467,
    "gf256_muladds": 3679216,
    "inactivated_columns": 135075,
}
EXPECTED_BASELINE_WEAK_ORDINALS = (354, 138273, 202780)

BASELINE_PROJECTION_FIELDS = (
    "K", "attempt", "root_index", "loss_root", "schedule_index",
    "schedule", "cell_ordinal", "outcome", "success", "rank_fail",
    "error", "binary_deficiency", "gf256_rank_gain",
    "inactivated_columns", "heavy_shortfall", "effective_precode_seed",
    "effective_packet_seed", "block_xors", "gf256_muladds",
)
BASELINE_MATRIX_FIELDS = frozenset((
    "schema", "record_ordinal", "K", "attempt", "root_index",
    "loss_root", "schedule_index", "schedule", "cell_ordinal",
    "trace_executed", "outcome", "success", "rank_fail", "error",
    "binary_deficiency", "gf256_rank_gain", "inactivated_columns",
    "heavy_shortfall", "effective_precode_seed", "effective_packet_seed",
    "block_xors", "gf256_muladds", "command_sha256", "stdout_sha256",
    "configuration_proof_sha256",
))

CSV_HEADER = tuple(stage0.CSV_HEADER)
METADATA_ORDER = tuple(stage0.METADATA_ORDER)
MAX_STDOUT_BYTES = 256 * 1024
MAX_STDERR_BYTES = 64 * 1024
PROCESS_TIMEOUT_SECONDS = 600
CHILD_ENVIRONMENT = {"LANG": "C", "LC_ALL": "C", "TZ": "UTC"}
UINT64_MAX = (1 << 64) - 1

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent
SOURCE_PATHS = (
    "CMakeLists.txt", "codec/CMakeLists.txt",
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
    "WirehairPeelFixups.inc", "bench/wh2_mix2_packet_offset_stage1.py",
    "bench/test_wh2_mix2_packet_offset_stage1.py",
    "bench/wh2_mix2_packet_offset_stage0.py",
    "bench/test_wh2_mix2_packet_offset_stage0.py",
    "bench/wh2_mix2_attempt_crossfit.py",
    "bench/wh2_mix2_promotion_short_screen.py",
)

COMMIT = re.compile(r"[0-9a-f]{40}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
UINT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
THREE_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")
EIGHT_DECIMAL = re.compile(r"(?:0|1)\.[0-9]{8}\Z")
HEX64 = re.compile(r"0x[0-9a-f]{16}\Z")
HEX32 = re.compile(r"0x[0-9a-f]{8}\Z")

# Replaced once after the literal contract body is finalized.
EXPECTED_CONTRACT_SHA256 = (
    "e8ebc4f4ce6f9b67051e39b522272bc37d1731ca691d57a299de49159f62d9e2")


class Stage1Error(RuntimeError):
    """The frozen Stage 1 replay cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise Stage1Error(message)


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


def _exact_fields(value: Mapping[str, Any], expected: Iterable[str],
                  context: str) -> None:
    if set(value) != set(expected):
        fail("{} fields differ from the frozen schema".format(context))


def _uint_value(value: Any, maximum: int, context: str) -> int:
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
        "bead": "wirehair-sxvz.16.1.20.38.2.1",
        "scope": "fixed packet-offset delta2 consumed-root Stage 1 replay only",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "candidate": {
            "base_profile_sha256": EXPECTED_CANDIDATE_PROFILE_SHA256,
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "mix_pair": CANDIDATE_PAIR,
            "pair_semantics": [0, 1],
            "packet_delta": PACKET_DELTA,
            "binary_dense_rows": BINARY_DENSE_ROWS,
            "gf256_heavy_rows": GF256_HEAVY_ROWS,
            "heavy_family": "periodic",
            "construction_seed_basis": "production-profile",
            "precode_attempt_rule": "exact frozen v9 p(K)",
            "packet_attempt_rule": "q=(p+2) mod 256",
            "packet_seed_rule": (
                "effective_q=(effective_p+(q-p)*0x9e3779b9) mod 2^32"),
            "profile_identity_rule": (
                "delta changes packet equations and requires a new immutable "
                "profile ID if later promoted; the base profile hash alone is "
                "not the candidate equation-profile identity"),
            "implementation_scope": (
                "exact split precode/packet attempts are exercised through a "
                "benchmark test hook only; PASS cannot license production "
                "routing, which requires a separately versioned implementation "
                "and review"),
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
            "candidate_observation_count": OBSERVATION_COUNT,
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "trials": 1,
            "threads": 1,
            "overhead_stream": "paired",
            "full_payload_solve": True,
            "solve_semantics": (
                "deterministic two-byte-wide all-zero-RHS rank/work solve; "
                "full-payload-solve does not make this payload-e2e byte "
                "reconstruction evidence"),
            "cold_solve_wide_xor": "policy",
            "current_diagonal_arm": {
                "observation_count": 0,
                "omission_basis": (
                    "the exact Stage-0 benchmark is reused; audited c9cfbc0-to-"
                    "56e6a2e mode-01 code maps pair01 to the unchanged ordinal "
                    "and equation identity, while authenticated Stage-0 controls "
                    "reproduced the three selected weak coordinates twice"),
                "stage0_control_scope": (
                    "three selected weak coordinates times two repetitions, not "
                    "a current 1080-coordinate diagonal replay"),
            },
        },
        "prerequisites": {
            "stage0": {
                "contract_sha256": EXPECTED_STAGE0_CONTRACT_SHA256,
                "source_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
                "helper_sha256": EXPECTED_STAGE0_HELPER_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in STAGE0_FILES.items()
                },
                "summary_self_sha256":
                    EXPECTED_STAGE0_SUMMARY_SELF_SHA256,
                "source_receipt_sha256":
                    EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256,
                "benchmark_sha256": EXPECTED_STAGE0_BENCH_SHA256,
                "benchmark_size": EXPECTED_STAGE0_BENCH_SIZE,
                "v9_receipt_sha256": EXPECTED_STAGE0_V9_RECEIPT_SHA256,
                "required_disposition": "PASS",
                "required_survivor_count": EXPECTED_STAGE0_SURVIVOR_COUNT,
                "required_lowest_survivor": PACKET_DELTA,
                "required_stage1_candidate_delta": PACKET_DELTA,
            },
            "v9": {
                "contract_sha256": EXPECTED_V9_CONTRACT_SHA256,
                "source_commit": EXPECTED_V9_SOURCE_COMMIT,
                "helper_sha256": EXPECTED_V9_HELPER_SHA256,
                "worker_sha256": EXPECTED_V9_WORKER_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in V9_FILES.items()
                },
                "input_self_sha256": EXPECTED_V9_INPUT_SELF_SHA256,
                "summary_self_sha256": EXPECTED_V9_SUMMARY_SELF_SHA256,
                "attempt_map_sha256": EXPECTED_ATTEMPT_MAP_SHA256,
                "candidate_matrix_crosscheck_observations":
                    V9_CANDIDATE_OBSERVATION_COUNT,
                "local_to_global_root_index": "+9",
            },
            "pair01_baseline": {
                "contract_sha256": EXPECTED_CROSSFIT_CONTRACT_SHA256,
                "source_commit": EXPECTED_CROSSFIT_SOURCE_COMMIT,
                "helper_sha256": EXPECTED_CROSSFIT_HELPER_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in BASELINE_FILES.items()
                },
                "input_self_sha256": EXPECTED_BASELINE_INPUT_SELF_SHA256,
                "summary_self_sha256":
                    EXPECTED_BASELINE_SUMMARY_SELF_SHA256,
                "matrix_record_count": BASELINE_MATRIX_RECORD_COUNT,
                "matrix_order": (
                    "K-major, then attempts 0..255, then "
                    "cell=root_index*3+schedule_index"),
                "record_ordinal_formula": (
                    "((K_index*256+attempt)*36+cell_ordinal)"),
                "selected_record_count": COORDINATE_COUNT,
                "selected_raw_stream_sha256":
                    EXPECTED_SELECTED_RAW_SHA256,
                "selected_projection_stream_sha256":
                    EXPECTED_SELECTED_PROJECTION_SHA256,
                "projection_fields": list(BASELINE_PROJECTION_FIELDS),
                "baseline_sums": dict(EXPECTED_BASELINE_SUMS),
                "weak_matrix_ordinals":
                    list(EXPECTED_BASELINE_WEAK_ORDINALS),
                "replay_file_forbidden": "attempt-crossfit-replay.jsonl",
            },
        },
        "transport": {
            "historical_paths": (
                "explicit --baseline-dir, --v9-dir, and --stage0-dir"),
            "historical_preflight": (
                "open non-symlink regular descriptors; exact size/hash; "
                "canonical parse; validate full matrix roster, V9 replay "
                "cross-check, and Stage-0 receipts before candidate work"),
            "benchmark_pin": {
                "rule": "open regular executable descriptor and rehash",
                "sha256": EXPECTED_STAGE0_BENCH_SHA256,
                "size": EXPECTED_STAGE0_BENCH_SIZE,
                "source_git_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
                "reuse_rule": (
                    "must be the exact benchmark that produced the authenticated "
                    "Stage-0 controls; no rebuilt current-HEAD binary is allowed"),
            },
            "benchmark_preflight": {
                "command": "--describe",
                "schema": BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
                "ordering": "last action before any precodefail workload",
            },
            "source_receipt": {
                "scope": "current Stage-1 controller and linked source checkout",
                "clean_at_HEAD": True,
                "working_source_bytes_equal_HEAD_blobs": True,
                "final_HEAD_recheck": True,
                "paths": list(SOURCE_PATHS),
            },
            "child_environment": dict(CHILD_ENVIRONMENT),
            "process_timeout_seconds": PROCESS_TIMEOUT_SECONDS,
            "publication": (
                "canonical fsynced JSONL and summary, fsynced staging, atomic "
                "no-replace directory rename, then parent-directory fsync"),
        },
        "gate": {
            "candidate_pair01_delta2": {
                "weak_coordinates": 0,
                "errors": 0,
                "configuration_failures": 0,
                "introductions_vs_pair01": 0,
                "repair_every_pair01_weak_coordinate": True,
                "full_repeat_deterministic_outcome_seed_rank_work_exact": True,
            },
            "deduplicated_aggregate_sums": {
                "rule": (
                    "one candidate value per coordinate after exact repeats "
                    "versus one extracted pair01 baseline value"),
                "block_xors": "candidate <= pair01",
                "gf256_muladds": "candidate <= pair01",
                "inactivated_columns": "candidate <= pair01",
            },
            "timing_used": False,
            "per_K_work_reported_without_extra_gate": True,
            "overall": "PASS iff every recovery, repeat, and work gate passes",
        },
        "stop_rule": (
            "PASS licenses only a separately preregistered next recovery or "
            "validation step; no fresh roots, timing claim, all-K campaign, "
            "production change, or promotion is licensed here"),
    }


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("MIX2 packet-offset Stage 1 contract differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _helper_hash(path: Path) -> str:
    return sha256_bytes(screen._stable_file_bytes(path, 4 * 1024 * 1024))


def _validate_constants() -> None:
    if (CANDIDATE_PAIR != "01" or PACKET_DELTA != 2 or
            PACKET_ATTEMPT_STRIDE != 0x9e3779b9 or REPETITIONS != (0, 1) or
            len(K_VALUES) != 30 or tuple(sorted(set(K_VALUES))) != K_VALUES or
            len(V9_ATTEMPTS) != len(K_VALUES) or
            any(type(value) is not int or not 0 <= value <= 255
                for value in V9_ATTEMPTS)):
        fail("Stage 1 candidate, K, attempt, or repetition roster changed")
    if (K_VALUES != tuple(screen.K_VALUES) or
            K_VALUES != tuple(crossfit.K_VALUES) or
            ROOTS != tuple(stage0.CONSUMED_ROOTS) or
            ROOTS != tuple(crossfit.ROOTS) or
            ROOTS != tuple(screen.SELECTION_ROOTS) + tuple(screen.ROOTS) or
            SCHEDULES != tuple(screen.SCHEDULES) or
            SCHEDULES != tuple(crossfit.SCHEDULES)):
        fail("Stage 1 coordinate provenance changed")
    if (set(ROOTS) & set(screen.FINAL_VALIDATION_ROOTS) or
            set(STAIRCASE_BY_K) != set(K_VALUES) or
            STAIRCASE_BY_K != dict(crossfit.STAIRCASE_BY_K)):
        fail("Stage 1 includes a fresh root or changed geometry")
    if (COORDINATE_COUNT != 1080 or OBSERVATION_COUNT != 2160 or
            BASELINE_MATRIX_RECORD_COUNT != 276480 or
            V9_RESULT_RECORD_COUNT != 1080 or
            V9_CANDIDATE_OBSERVATION_COUNT != 540):
        fail("Stage 1 domain cardinality changed")
    if (tuple(CSV_HEADER) != tuple(screen.CSV_HEADER) or
            tuple(CSV_HEADER) != tuple(stage0.CSV_HEADER)):
        fail("precodefail CSV protocol changed")
    if (stage0._packet_seed_for_offset("0x12345678", 255, PACKET_DELTA) !=
            "0x172990ea"):
        fail("wrapped packet-offset seed rule changed")
    if (screen.candidate_profile_sha256() !=
            EXPECTED_CANDIDATE_PROFILE_SHA256):
        fail("Stage 1 base candidate profile changed")

    stage0._validate_constants()
    crossfit._validate_constants()
    if (stage0.contract_description()["contract_sha256"] !=
            EXPECTED_STAGE0_CONTRACT_SHA256 or
            screen.contract_description()["contract_sha256"] !=
            EXPECTED_V9_CONTRACT_SHA256 or
            crossfit.contract_description()["contract_sha256"] !=
            EXPECTED_CROSSFIT_CONTRACT_SHA256):
        fail("an imported prerequisite contract changed")
    helper_receipts = {
        Path(stage0.__file__).resolve(): EXPECTED_STAGE0_HELPER_SHA256,
        Path(screen.__file__).resolve(): EXPECTED_V9_HELPER_SHA256,
        Path(crossfit.__file__).resolve(): EXPECTED_CROSSFIT_HELPER_SHA256,
    }
    for path, expected in helper_receipts.items():
        if _helper_hash(path) != expected:
            fail("imported helper differs from frozen source: {}".format(path))
    for K_index, K in enumerate(K_VALUES):
        for attempt in (0, ATTEMPT_BY_K[K], 255):
            for cell in (0, 35):
                expected = ((K_index * 256 + attempt) * 36 + cell)
                if crossfit._matrix_record_ordinal(K, attempt, cell) != expected:
                    fail("cross-fit matrix ordinal formula changed")


@dataclass(frozen=True)
class Coordinate:
    ordinal: int

    def __post_init__(self) -> None:
        if type(self.ordinal) is not int or not 0 <= self.ordinal < 1080:
            fail("coordinate ordinal is outside Stage 1")

    @property
    def K_index(self) -> int:
        return self.ordinal // (len(ROOTS) * len(SCHEDULES))

    @property
    def K(self) -> int:
        return K_VALUES[self.K_index]

    @property
    def cell_ordinal(self) -> int:
        return self.ordinal % (len(ROOTS) * len(SCHEDULES))

    @property
    def root_index(self) -> int:
        return self.cell_ordinal // len(SCHEDULES)

    @property
    def schedule_index(self) -> int:
        return self.cell_ordinal % len(SCHEDULES)

    @property
    def root(self) -> str:
        return ROOTS[self.root_index]

    @property
    def schedule(self) -> str:
        return SCHEDULES[self.schedule_index]

    @property
    def attempt(self) -> int:
        return ATTEMPT_BY_K[self.K]

    @property
    def packet_attempt(self) -> int:
        return (self.attempt + PACKET_DELTA) & 255

    @property
    def baseline_record_ordinal(self) -> int:
        return ((self.K_index * 256 + self.attempt) * 36 +
                self.cell_ordinal)

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
                self.ordinal !=
                    self.coordinate_ordinal * len(REPETITIONS) +
                    self.repetition):
            fail("invocation identity differs from canonical Stage 1 order")

    @property
    def coordinate(self) -> Coordinate:
        return Coordinate(self.coordinate_ordinal)

    @property
    def pair(self) -> str:
        return CANDIDATE_PAIR

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
    invocations = tuple(
        Invocation(coordinate * 2 + repetition, coordinate, repetition)
        for coordinate in range(COORDINATE_COUNT)
        for repetition in REPETITIONS)
    if len(invocations) != OBSERVATION_COUNT or len({
            (item.coordinate_ordinal, item.repetition)
            for item in invocations}) != OBSERVATION_COUNT:
        fail("Stage 1 invocation roster is incomplete or duplicated")
    return invocations


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
class PinnedFile:
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
                    self.descriptor, min(1024 * 1024, self.size - offset),
                    offset)
            except OSError as exc:
                fail("cannot read {}: {}".format(self.label, exc))
            if not chunk:
                fail("{} was truncated while reading".format(self.label))
            chunks.append(chunk)
            offset += len(chunk)
        if self.receipt()["size"] != len(b"".join(chunks)):
            fail("{} read length changed".format(self.label))
        return b"".join(chunks)

    def lines(self, maximum_line: int = 16 * 1024) \
            -> Iterator[Tuple[int, bytes]]:
        try:
            duplicated = os.dup(self.descriptor)
        except OSError as exc:
            fail("cannot duplicate {}: {}".format(self.label, exc))
        with os.fdopen(duplicated, "rb", closefd=True) as stream:
            stream.seek(0)
            for ordinal, line in enumerate(stream):
                if (not line.endswith(b"\n") or b"\r" in line or
                        len(line) > maximum_line or line == b"\n"):
                    fail("{} line {} has invalid framing".format(
                        self.label, ordinal))
                yield ordinal, line
        self.receipt()

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __enter__(self) -> "PinnedFile":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def _open_pinned_file(path: Path, expected_size: int,
                      expected_sha256: str, label: str) -> PinnedFile:
    descriptor = -1
    try:
        original = path.lstat()
        if stat.S_ISLNK(original.st_mode) or not stat.S_ISREG(original.st_mode):
            fail("{} must be a non-symlink regular file".format(label))
        resolved = path.resolve(strict=True)
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        cloexec = getattr(os, "O_CLOEXEC", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed".format(label))
        descriptor = os.open(str(resolved), os.O_RDONLY | cloexec | nofollow)
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
        result = PinnedFile(resolved, descriptor, digest, size, label)
        descriptor = -1
        return result
    except OSError as exc:
        fail("cannot inspect {}: {}".format(label, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    raise AssertionError("unreachable")


def _validate_historical_directory(path: Path, label: str) -> Path:
    try:
        state = path.lstat()
    except OSError as exc:
        fail("cannot inspect {} directory: {}".format(label, exc))
    if stat.S_ISLNK(state.st_mode) or not stat.S_ISDIR(state.st_mode):
        fail("{} must be an explicit non-symlink directory".format(label))
    return path.resolve(strict=True)


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


def _validate_self_hash(value: Mapping[str, Any], field: str,
                        expected: str, context: str) -> None:
    if value.get(field) != expected or self_hash(value, field) != expected:
        fail("{} self hash is inconsistent".format(context))


def _baseline_projection(row: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in BASELINE_PROJECTION_FIELDS
               if field not in row]
    if missing:
        fail("baseline row omits projection fields: {}".format(missing))
    return {field: row[field] for field in BASELINE_PROJECTION_FIELDS}


def _validate_matrix_row(row: Mapping[str, Any], ordinal: int) -> None:
    if set(row) != BASELINE_MATRIX_FIELDS:
        fail("baseline matrix row {} fields changed".format(ordinal))
    K_index, remainder = divmod(ordinal, 256 * 36)
    attempt, cell_ordinal = divmod(remainder, 36)
    root_index, schedule_index = divmod(cell_ordinal, len(SCHEDULES))
    if not 0 <= K_index < len(K_VALUES):
        fail("baseline matrix has excess rows")
    K = K_VALUES[K_index]
    expected_identity = {
        "schema": crossfit.RECORD_SCHEMA,
        "record_ordinal": ordinal,
        "K": K,
        "attempt": attempt,
        "root_index": root_index,
        "loss_root": ROOTS[root_index],
        "schedule_index": schedule_index,
        "schedule": SCHEDULES[schedule_index],
        "cell_ordinal": cell_ordinal,
    }
    if any(row.get(field) != expected
           for field, expected in expected_identity.items()):
        fail("baseline matrix row {} has noncanonical identity".format(
            ordinal))
    if (crossfit._matrix_record_ordinal(K, attempt, cell_ordinal) != ordinal or
            row["trace_executed"] is not True or
            row["configuration_proof_sha256"] is not None):
        fail("baseline matrix row {} violates execution provenance".format(
            ordinal))
    for field in ("success", "rank_fail", "error"):
        if type(row[field]) is not bool:
            fail("baseline matrix row {} terminal type changed".format(
                ordinal))
    if sum(int(row[field]) for field in ("success", "rank_fail", "error")) != 1:
        fail("baseline matrix row {} terminal state is ambiguous".format(
            ordinal))
    expected_outcome = (
        "success" if row["success"] else
        ("rank_fail" if row["rank_fail"] else "error"))
    if row["outcome"] != expected_outcome:
        fail("baseline matrix row {} outcome disagrees".format(ordinal))
    for field in (
            "binary_deficiency", "gf256_rank_gain", "inactivated_columns",
            "heavy_shortfall", "block_xors", "gf256_muladds"):
        _uint_value(row[field], UINT64_MAX,
                    "baseline row {} {}".format(ordinal, field))
    if (row["gf256_rank_gain"] > row["binary_deficiency"] or
            row["binary_deficiency"] > row["inactivated_columns"] or
            row["gf256_rank_gain"] > GF256_HEAVY_ROWS or
            (row["success"] and
             row["gf256_rank_gain"] != row["binary_deficiency"]) or
            (row["rank_fail"] and
             row["gf256_rank_gain"] >= row["binary_deficiency"])):
        fail("baseline matrix row {} rank counters disagree".format(ordinal))
    expected_shortfall = int(
        row["rank_fail"] and row["binary_deficiency"] <= GF256_HEAVY_ROWS and
        row["gf256_rank_gain"] < row["binary_deficiency"])
    if row["heavy_shortfall"] != expected_shortfall:
        fail("baseline matrix row {} shortfall disagrees".format(ordinal))
    if (HEX64.fullmatch(row["effective_precode_seed"]) is None or
            HEX32.fullmatch(row["effective_packet_seed"]) is None):
        fail("baseline matrix row {} seed is malformed".format(ordinal))
    _sha(row["command_sha256"], "baseline command receipt")
    _sha(row["stdout_sha256"], "baseline stdout receipt")


def _load_baseline(pins: Mapping[str, PinnedFile]) \
        -> Tuple[Tuple[Mapping[str, Any], ...], Mapping[str, Any]]:
    input_value = _parse_canonical_json_bytes(
        pins["attempt-crossfit-input.json"].read_bytes(1024 * 1024),
        "cross-fit input")
    summary = _parse_canonical_json_bytes(
        pins["attempt-crossfit-summary.json"].read_bytes(1024 * 1024),
        "cross-fit summary")
    _validate_self_hash(
        input_value, "input_sha256", EXPECTED_BASELINE_INPUT_SELF_SHA256,
        "cross-fit input")
    _validate_self_hash(
        summary, "summary_sha256", EXPECTED_BASELINE_SUMMARY_SELF_SHA256,
        "cross-fit summary")
    if (input_value.get("schema") != crossfit.INPUT_SCHEMA or
            input_value.get("contract") != crossfit.contract_description() or
            input_value.get("controller", {}).get("sha256") !=
                EXPECTED_CROSSFIT_HELPER_SHA256 or
            input_value.get("controller", {}).get("source_git_commit") !=
                EXPECTED_CROSSFIT_SOURCE_COMMIT):
        fail("cross-fit input provenance changed")
    if (summary.get("schema") != crossfit.SUMMARY_SCHEMA or
            summary.get("contract_sha256") !=
                EXPECTED_CROSSFIT_CONTRACT_SHA256 or
            summary.get("input_sha256") !=
                EXPECTED_BASELINE_INPUT_SELF_SHA256 or
            summary.get("input_file_sha256") !=
                BASELINE_FILES["attempt-crossfit-input.json"][1] or
            summary.get("matrix_sha256") !=
                BASELINE_FILES["attempt-crossfit-matrix.jsonl"][1] or
            summary.get("matrix_record_count") !=
                BASELINE_MATRIX_RECORD_COUNT or
            summary.get("controller_sha256") !=
                EXPECTED_CROSSFIT_HELPER_SHA256 or
            summary.get("source_git_commit") !=
                EXPECTED_CROSSFIT_SOURCE_COMMIT or
            summary.get("candidate_profile_sha256") !=
                EXPECTED_CANDIDATE_PROFILE_SHA256 or
            summary.get("configuration_failed_K_attempts") != 0):
        fail("cross-fit summary provenance changed")

    selected: List[Mapping[str, Any]] = []
    selected_raw = hashlib.sha256()
    selected_projection = hashlib.sha256()
    weak_ordinals: List[int] = []
    sums = {field: 0 for field in EXPECTED_BASELINE_SUMS}
    count = 0
    for ordinal, line in pins["attempt-crossfit-matrix.jsonl"].lines():
        row = _parse_canonical_json_line(
            line, "baseline matrix row {}".format(ordinal))
        _validate_matrix_row(row, ordinal)
        count += 1
        if row["attempt"] == ATTEMPT_BY_K[row["K"]]:
            selected_raw.update(line)
            projection_line = (
                canonical_json(_baseline_projection(row)) + "\n").encode(
                    "ascii")
            selected_projection.update(projection_line)
            selected.append(row)
            for field in sums:
                sums[field] += row[field]
            if not row["success"]:
                weak_ordinals.append(row["record_ordinal"])
    if count != BASELINE_MATRIX_RECORD_COUNT:
        fail("baseline matrix roster has the wrong cardinality")
    if len(selected) != COORDINATE_COUNT:
        fail("baseline selected roster has the wrong cardinality")
    for coordinate_ordinal, row in enumerate(selected):
        coordinate = Coordinate(coordinate_ordinal)
        if (row["record_ordinal"] != coordinate.baseline_record_ordinal or
                any(row[field] != value for field, value in
                    coordinate.identity().items()
                    if field not in ("coordinate_ordinal",
                                     "baseline_record_ordinal"))):
            fail("baseline selected roster order changed")
    if (selected_raw.hexdigest() != EXPECTED_SELECTED_RAW_SHA256 or
            selected_projection.hexdigest() !=
                EXPECTED_SELECTED_PROJECTION_SHA256 or
            sums != EXPECTED_BASELINE_SUMS or
            tuple(weak_ordinals) != EXPECTED_BASELINE_WEAK_ORDINALS):
        fail("baseline extraction receipt differs from the frozen audit")
    receipt = {
        "input_self_sha256": EXPECTED_BASELINE_INPUT_SELF_SHA256,
        "summary_self_sha256": EXPECTED_BASELINE_SUMMARY_SELF_SHA256,
        "matrix_record_count": count,
        "selected_record_count": len(selected),
        "selected_raw_stream_sha256": selected_raw.hexdigest(),
        "selected_projection_stream_sha256":
            selected_projection.hexdigest(),
        "selected_sums": sums,
        "selected_weak_matrix_ordinals": weak_ordinals,
        "excluded_replay_file": "attempt-crossfit-replay.jsonl",
        "excluded_replay_opened": False,
    }
    return tuple(selected), receipt


def _load_stage0(pins: Mapping[str, PinnedFile],
                 v9_receipt: Mapping[str, Any],
                 baseline: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    summary = _parse_canonical_json_bytes(
        pins["mix2-packet-offset-stage0-summary.json"].read_bytes(
            1024 * 1024), "packet-offset Stage-0 summary")
    _validate_self_hash(
        summary, "summary_sha256", EXPECTED_STAGE0_SUMMARY_SELF_SHA256,
        "packet-offset Stage-0 summary")
    if (summary.get("schema") != stage0.SUMMARY_SCHEMA or
            summary.get("contract_sha256") !=
                EXPECTED_STAGE0_CONTRACT_SHA256 or
            summary.get("records_file") !=
                "mix2-packet-offset-stage0-records.jsonl" or
            summary.get("records_sha256") != STAGE0_FILES[
                "mix2-packet-offset-stage0-records.jsonl"][1] or
            summary.get("record_count") != stage0.EXPECTED_INVOCATION_COUNT or
            summary.get("disposition") != "PASS" or
            summary.get("stage1_candidate_delta") != PACKET_DELTA or
            summary.get("stage0_only") is not True or
            summary.get("promotion_evidence") is not False or
            summary.get("fresh_roots_used") is not False or
            summary.get("timing_evidence") is not False):
        fail("Stage-0 summary does not license packet offset delta 2")
    source = summary.get("source_receipt")
    bench = summary.get("bench")
    description = summary.get("bench_description")
    v9 = summary.get("v9_prerequisite")
    if (type(source) is not dict or
            source.get("source_git_commit") != EXPECTED_STAGE0_SOURCE_COMMIT or
            source.get("source_receipt_sha256") !=
                EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 or
            source.get("source_receipt_sha256") !=
                stage0.self_hash(source, "source_receipt_sha256") or
            type(bench) is not dict or
            bench.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench.get("size") != EXPECTED_STAGE0_BENCH_SIZE or
            type(description) is not dict or
            description.get("schema") != BENCH_DESCRIPTION_SCHEMA or
            description.get("source_git_commit") !=
                EXPECTED_STAGE0_SOURCE_COMMIT or
            type(v9) is not dict or
            v9.get("receipt_sha256") != EXPECTED_STAGE0_V9_RECEIPT_SHA256 or
            v9.get("receipt_sha256") !=
                stage0.self_hash(v9, "receipt_sha256") or
            v9.get("contract_sha256") != EXPECTED_V9_CONTRACT_SHA256 or
            v9.get("source_git_commit") != EXPECTED_V9_SOURCE_COMMIT or
            v9.get("worker_sha256") != EXPECTED_V9_WORKER_SHA256 or
            v9.get("files") is None):
        fail("Stage-0 source, benchmark, or V9 receipt is inconsistent")
    v9_files = v9["files"]
    if (type(v9_files) is not dict or set(v9_files) != set(V9_FILES)):
        fail("Stage-0 embedded V9 file roster changed")
    for name, (size, digest) in V9_FILES.items():
        value = v9_files[name]
        if (type(value) is not dict or value.get("size") != size or
                value.get("sha256") != digest):
            fail("Stage-0 embedded V9 file receipt changed")
    if (v9_receipt.get("attempt_stream_sha256") !=
            V9_FILES["promotion-short-screen-attempts.jsonl"][1] or
            v9_receipt.get("result_stream_sha256") !=
            V9_FILES["promotion-short-screen-results.jsonl"][1] or
            v9_receipt.get("source_git_commit") != EXPECTED_V9_SOURCE_COMMIT):
        fail("Stage-0 embedded V9 receipt differs from independent V9 audit")

    if len(baseline) != COORDINATE_COUNT:
        fail("Stage-0 audit received an incomplete pair01 baseline")
    expected_precode_seeds = []
    expected_packet_seeds = []
    expected_diagonal_control = []
    for stage0_coordinate in stage0.COORDINATES:
        K_index = K_VALUES.index(stage0_coordinate.K)
        schedule_index = SCHEDULES.index(stage0_coordinate.schedule)
        coordinate_ordinal = (
            K_index * len(ROOTS) * len(SCHEDULES) +
            stage0_coordinate.root_index * len(SCHEDULES) + schedule_index)
        row = baseline[coordinate_ordinal]
        if (row["K"] != stage0_coordinate.K or
                row["attempt"] != stage0_coordinate.attempt or
                row["root_index"] != stage0_coordinate.root_index or
                row["loss_root"] != stage0_coordinate.root or
                row["schedule"] != stage0_coordinate.schedule):
            fail("Stage-0 weak coordinate differs from pair01 baseline")
        expected_precode_seeds.append(row["effective_precode_seed"])
        expected_packet_seeds.append(row["effective_packet_seed"])
        expected_diagonal_control.append({
            "staircase": STAIRCASE_BY_K[row["K"]],
            "source_hits": 3 if row["K"] >= 10000 else 2,
            "outcome": row["outcome"],
            "success": row["success"],
            "rank_fail": row["rank_fail"],
            "error": row["error"],
            "inactivated_columns": row["inactivated_columns"],
            "effective_precode_seed": row["effective_precode_seed"],
            "effective_packet_seed": row["effective_packet_seed"],
            "block_xors": row["block_xors"],
            "gf256_muladds": row["gf256_muladds"],
        })
    if (v9.get("effective_precode_seeds") != expected_precode_seeds or
            v9.get("effective_packet_seeds") != expected_packet_seeds or
            v9.get("diagonal_control") != expected_diagonal_control):
        fail("Stage-0 embedded V9 seed or control roster changed")

    records: List[Mapping[str, Any]] = []
    for ordinal, line in pins[
            "mix2-packet-offset-stage0-records.jsonl"].lines():
        row = _parse_canonical_json_line(
            line, "packet-offset Stage-0 record {}".format(ordinal))
        if (row.get("schema") != stage0.RECORD_SCHEMA or
                row.get("ordinal") != ordinal or
                row.get("bench_binary_sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
                row.get("source_git_commit") != EXPECTED_STAGE0_SOURCE_COMMIT or
                row.get("source_receipt_sha256") !=
                    EXPECTED_STAGE0_SOURCE_RECEIPT_SHA256 or
                row.get("deterministic_projection_sha256") !=
                    stage0.sha256_json(stage0.deterministic_projection(row))):
            fail("packet-offset Stage-0 record {} receipt changed".format(
                ordinal))
        records.append(row)
    if len(records) != stage0.EXPECTED_INVOCATION_COUNT:
        fail("packet-offset Stage-0 record roster has the wrong cardinality")
    verdict = stage0.adjudicate(
        records,
        dict(enumerate(expected_precode_seeds)),
        dict(enumerate(expected_packet_seeds)),
        dict(enumerate(expected_diagonal_control)))
    for field in (
            "diagonal_control", "delta_results", "survivors",
            "stage1_candidate_delta", "effective_precode_seeds",
            "authenticated_effective_p_packet_seeds", "disposition"):
        if summary.get(field) != verdict[field]:
            fail("packet-offset Stage-0 adjudication receipt changed")
    if (len(verdict["survivors"]) != EXPECTED_STAGE0_SURVIVOR_COUNT or
            not verdict["survivors"] or
            verdict["survivors"][0] != PACKET_DELTA):
        fail("packet-offset Stage-0 lowest-survivor rule changed")
    return {
        "summary_self_sha256": EXPECTED_STAGE0_SUMMARY_SELF_SHA256,
        "record_count": len(records),
        "records_sha256": STAGE0_FILES[
            "mix2-packet-offset-stage0-records.jsonl"][1],
        "source_git_commit": EXPECTED_STAGE0_SOURCE_COMMIT,
        "bench_binary_sha256": bench["sha256"],
        "source_receipt_sha256": source["source_receipt_sha256"],
        "v9_receipt_sha256": v9["receipt_sha256"],
        "survivor_count": len(verdict["survivors"]),
        "stage1_candidate_delta": verdict["stage1_candidate_delta"],
    }


def _v9_candidate_projection(row: Mapping[str, Any]) -> Mapping[str, Any]:
    return {
        "success": bool(row["success"]),
        "rank_fail": bool(row["rank_fail"]),
        "error": bool(row["error"]),
        "inactivated_columns": row["inactivated_columns"],
        "block_xors": row["block_xors"],
        "gf256_muladds": row["gf256_muladds"],
        "effective_precode_seed": row["effective_precode_seed"],
        "effective_packet_seed": row["effective_packet_seed"],
    }


def _baseline_v9_projection(row: Mapping[str, Any]) -> Mapping[str, Any]:
    return {
        "success": row["success"],
        "rank_fail": row["rank_fail"],
        "error": row["error"],
        "inactivated_columns": row["inactivated_columns"],
        "block_xors": row["block_xors"],
        "gf256_muladds": row["gf256_muladds"],
        "effective_precode_seed": row["effective_precode_seed"],
        "effective_packet_seed": row["effective_packet_seed"],
    }


def _load_v9(pins: Mapping[str, PinnedFile],
             baseline: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    attempts_bytes = pins[
        "promotion-short-screen-attempts.jsonl"].read_bytes(1024 * 1024)
    derivations = screen.parse_derivation_stream(
        attempts_bytes, EXPECTED_V9_WORKER_SHA256)
    realized_attempts = tuple(record["selected_attempt"]
                              for record in derivations)
    attempt_pairs = [[K, ATTEMPT_BY_K[K]] for K in K_VALUES]
    if (realized_attempts != V9_ATTEMPTS or
            sha256_json(attempt_pairs) != EXPECTED_ATTEMPT_MAP_SHA256):
        fail("V9 attempt stream differs from the frozen attempt map")

    input_value = _parse_canonical_json_bytes(
        pins["promotion-short-screen-input.json"].read_bytes(2 * 1024 * 1024),
        "V9 input")
    summary = _parse_canonical_json_bytes(
        pins["promotion-short-screen-summary.json"].read_bytes(1024 * 1024),
        "V9 summary")
    _validate_self_hash(
        input_value, "input_sha256", EXPECTED_V9_INPUT_SELF_SHA256,
        "V9 input")
    _validate_self_hash(
        summary, "summary_sha256", EXPECTED_V9_SUMMARY_SELF_SHA256,
        "V9 summary")
    if (input_value.get("contract") != screen.contract_description() or
            input_value.get("contract_sha256") != EXPECTED_V9_CONTRACT_SHA256 or
            input_value.get("attempt_selection_stream_sha256") !=
                V9_FILES["promotion-short-screen-attempts.jsonl"][1] or
            input_value.get("attempt_selection_record_count") != 30 or
            input_value.get("repair_worker_binary", {}).get("sha256") !=
                EXPECTED_V9_WORKER_SHA256):
        fail("V9 input provenance changed")
    invocations = input_value.get("invocations")
    if type(invocations) is not list or len(invocations) != V9_RESULT_RECORD_COUNT:
        fail("V9 input invocation roster has the wrong cardinality")
    if (summary.get("schema") != screen.SUMMARY_SCHEMA or
            summary.get("contract_sha256") != EXPECTED_V9_CONTRACT_SHA256 or
            summary.get("input_sha256") != EXPECTED_V9_INPUT_SELF_SHA256 or
            summary.get("input_file_sha256") !=
                V9_FILES["promotion-short-screen-input.json"][1] or
            summary.get("result_stream_sha256") !=
                V9_FILES["promotion-short-screen-results.jsonl"][1] or
            summary.get("attempt_selection_stream_sha256") !=
                V9_FILES["promotion-short-screen-attempts.jsonl"][1] or
            summary.get("attempt_selection_record_count") != 30 or
            summary.get("record_count") != V9_RESULT_RECORD_COUNT or
            summary.get("source_git_commit") != EXPECTED_V9_SOURCE_COMMIT or
            summary.get("repair_worker_binary_sha256") !=
                EXPECTED_V9_WORKER_SHA256 or
            summary.get("candidate_profile_sha256") !=
                EXPECTED_CANDIDATE_PROFILE_SHA256):
        fail("V9 summary provenance changed")

    seen: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    total = 0
    for ordinal, line in pins["promotion-short-screen-results.jsonl"].lines():
        row = _parse_canonical_json_line(line, "V9 result {}".format(ordinal))
        if ordinal >= len(invocations):
            fail("V9 result roster has excess rows")
        invocation = invocations[ordinal]
        if (type(invocation) is not dict or
                row.get("schema") != screen.RESULT_SCHEMA or
                row.get("ordinal") != ordinal or
                invocation.get("ordinal") != ordinal):
            fail("V9 result {} has noncanonical identity".format(ordinal))
        for field in (
                "K", "root_index", "loss_root", "schedule", "cell_ordinal",
                "coordinate_ordinal", "arm", "construction_attempt",
                "timing_order", "timing_slot", "observation_index",
                "command_sha256"):
            if row.get(field) != invocation.get(field):
                fail("V9 result {} differs from its input receipt".format(
                    ordinal))
        total += 1
        if row["arm"] != "candidate_two07_mix2":
            continue
        K = row["K"]
        if K not in ATTEMPT_BY_K:
            fail("V9 candidate result names an unknown K")
        local_root = _uint_value(row["root_index"], 2, "V9 local root")
        schedule_index = SCHEDULES.index(row["schedule"])
        local_cell = local_root * len(SCHEDULES) + schedule_index
        K_index = K_VALUES.index(K)
        expected_local_coordinate = K_index * 9 + local_cell
        global_coordinate = K_index * 36 + (local_root + 9) * 3 + schedule_index
        observation = _uint_value(
            row["observation_index"], 1, "V9 observation index")
        if (row["coordinate_ordinal"] != expected_local_coordinate or
                row["cell_ordinal"] != local_cell or
                row["loss_root"] != ROOTS[local_root + 9] or
                row["construction_attempt"] != ATTEMPT_BY_K[K] or
                row.get("attempt_mode") != "exact" or
                row.get("mix_count") != 2 or
                row.get("dense_anchor_layout") != "two07" or
                row.get("block_bytes") != BLOCK_BYTES or
                row.get("loss_ppm") != LOSS_PPM or
                row.get("overhead") != OVERHEAD or
                row.get("full_payload_byte_recovery_verified") is not
                    bool(row["success"]) or
                row.get("weak") is not (row["success"] == 0)):
            fail("V9 candidate result differs from pair01 baseline protocol")
        for field in ("success", "rank_fail", "error"):
            if type(row[field]) is not int or row[field] not in (0, 1):
                fail("V9 candidate terminal count is not Boolean")
        if row["success"] + row["rank_fail"] + row["error"] != 1:
            fail("V9 candidate terminal state is ambiguous")
        key = (global_coordinate, observation)
        if key in seen:
            fail("V9 candidate cross-check duplicated an observation")
        baseline_row = baseline[global_coordinate]
        if _v9_candidate_projection(row) != \
                _baseline_v9_projection(baseline_row):
            fail("V9 candidate result does not reproduce cross-fit matrix")
        seen[key] = row
    expected_seen = {
        (K_index * 36 + root * 3 + schedule, observation)
        for K_index in range(len(K_VALUES))
        for root in range(9, 12)
        for schedule in range(3)
        for observation in REPETITIONS
    }
    if (total != V9_RESULT_RECORD_COUNT or set(seen) != expected_seen or
            len(seen) != V9_CANDIDATE_OBSERVATION_COUNT):
        fail("V9 result cross-check roster is incomplete")
    return {
        "input_self_sha256": EXPECTED_V9_INPUT_SELF_SHA256,
        "summary_self_sha256": EXPECTED_V9_SUMMARY_SELF_SHA256,
        "attempt_stream_sha256":
            V9_FILES["promotion-short-screen-attempts.jsonl"][1],
        "attempt_map_sha256": EXPECTED_ATTEMPT_MAP_SHA256,
        "result_stream_sha256":
            V9_FILES["promotion-short-screen-results.jsonl"][1],
        "result_record_count": total,
        "candidate_crosscheck_observations": len(seen),
        "source_git_commit": EXPECTED_V9_SOURCE_COMMIT,
    }


@dataclass
class HistoricalEvidence:
    stack: ExitStack
    pins: Mapping[str, Mapping[str, PinnedFile]]
    baseline: Tuple[Mapping[str, Any], ...]
    baseline_receipt: Mapping[str, Any]
    v9_receipt: Mapping[str, Any]
    stage0_receipt: Mapping[str, Any]

    def file_receipts(self) -> Mapping[str, Any]:
        return {
            group: {
                name: pin.receipt() for name, pin in sorted(values.items())
            }
            for group, values in sorted(self.pins.items())
        }

    def receipt(self) -> Mapping[str, Any]:
        return {
            "files": self.file_receipts(),
            "baseline": self.baseline_receipt,
            "v9": self.v9_receipt,
            "stage0": self.stage0_receipt,
        }

    def close(self) -> None:
        self.stack.close()

    def __enter__(self) -> "HistoricalEvidence":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def _pin_group(stack: ExitStack, directory: Path,
               specs: Mapping[str, Tuple[int, str]], label: str) \
        -> Mapping[str, PinnedFile]:
    root = _validate_historical_directory(directory, label)
    return {
        name: stack.enter_context(_open_pinned_file(
            root / name, size, digest, "{} {}".format(label, name)))
        for name, (size, digest) in specs.items()
    }


def load_historical_evidence(baseline_dir: Path, v9_dir: Path,
                             stage0_dir: Path) -> HistoricalEvidence:
    stack = ExitStack()
    try:
        pins = {
            "baseline": _pin_group(
                stack, baseline_dir, BASELINE_FILES, "pair01 baseline"),
            "v9": _pin_group(stack, v9_dir, V9_FILES, "V9"),
            "stage0": _pin_group(
                stack, stage0_dir, STAGE0_FILES, "Stage-0"),
        }
        baseline, baseline_receipt = _load_baseline(pins["baseline"])
        v9_receipt = _load_v9(pins["v9"], baseline)
        stage0_receipt = _load_stage0(
            pins["stage0"], v9_receipt, baseline)
        evidence = HistoricalEvidence(
            stack, pins, baseline, baseline_receipt, v9_receipt,
            stage0_receipt)
        evidence.receipt()
        return evidence
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


def _parse_uint_text(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or UINT.fullmatch(text) is None:
        fail("{} is not canonical unsigned decimal".format(context))
    if len(text) > len(str(maximum)):
        fail("{} exceeds its bound".format(context))
    value = int(text)
    if value > maximum:
        fail("{} exceeds its bound".format(context))
    return value


def _parse_integral_decimal(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical three-place decimal".format(context))
    whole, fraction = text.split(".")
    if fraction != "000":
        fail("{} is not integral for a one-trial cell".format(context))
    return _parse_uint_text(whole, maximum, context)


def _parse_timing_decimal(text: str, context: str) -> None:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical nonnegative timing".format(context))
    _parse_uint_text(text.split(".")[0], UINT64_MAX, context)


def _parse_histogram(text: str, expected_key: int, context: str) -> None:
    if text != "{}:1".format(expected_key):
        fail("{} differs from its one-trial counter".format(context))


def _parse_metadata(line: str, invocation: Invocation,
                    bench_source_commit: str) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("precodefail metadata is missing")
    values: Dict[str, str] = {}
    order: List[str] = []
    for token in line[len(prefix):].split(" "):
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
        "exact_packet_attempt": str(coordinate.packet_attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
        "source_git_commit": bench_source_commit, "mix_pair": CANDIDATE_PAIR,
    }
    if values != expected:
        fail("precodefail metadata differs from the frozen Stage 1 cell")
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


def _framed_lines(result: ProcessResult, bench_source_commit: str) \
        -> Tuple[str, ...]:
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
    _parse_metadata(lines[0], result.invocation, bench_source_commit)
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
        "packet_attempt": str(coordinate.packet_attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
    }
    for field, value in expected.items():
        if row[field] != value:
            fail("precodefail field {} differs from Stage 1".format(field))
    staircase = _parse_uint_text(row["staircase"], 4096, "staircase")
    source_hits = _parse_uint_text(row["source_hits"], 8, "source hits")
    if (staircase != STAIRCASE_BY_K[coordinate.K] or
            source_hits != (3 if coordinate.K >= 10000 else 2)):
        fail("precodefail realized profile geometry changed")

    success = _parse_uint_text(row["success"], 1, "success")
    rank_fail = _parse_uint_text(row["rank_fail"], 1, "rank failure")
    error = _parse_uint_text(row["error"], 1, "error")
    if success + rank_fail + error != 1:
        fail("precodefail terminal counts do not sum to one")
    weak = success != 1
    expected_rate = "0.00000000" if success else "1.00000000"
    if (EIGHT_DECIMAL.fullmatch(row["fail_rate"]) is None or
            row["fail_rate"] != expected_rate):
        fail("precodefail failure rate disagrees with its outcome")

    inactivated = _parse_integral_decimal(
        row["inact_mu"], 65535, "inactivation")
    if _parse_uint_text(row["inact_max"], 65535, "maximum inactivation") != \
            inactivated:
        fail("precodefail inactivation counters disagree")
    deficiency = _parse_integral_decimal(
        row["binary_def_mu"], 65535, "binary deficiency")
    gain = _parse_integral_decimal(
        row["heavy_gain_mu"], 65535, "GF256 rank gain")
    if (_parse_uint_text(row["binary_def_max"], 65535,
                         "maximum binary deficiency") != deficiency or
            _parse_uint_text(row["heavy_gain_min"], 65535,
                             "minimum GF256 rank gain") != gain or
            deficiency > inactivated or gain > deficiency or
            gain > GF256_HEAVY_ROWS):
        fail("precodefail rank counters are inconsistent")
    if success and gain != deficiency:
        fail("precodefail success is rank-deficient")
    if rank_fail and gain >= deficiency:
        fail("precodefail rank failure has full rank")
    shortfall = _parse_uint_text(row["heavy_shortfall"], 1,
                                 "heavy shortfall")
    expected_shortfall = int(
        bool(rank_fail) and deficiency <= GF256_HEAVY_ROWS and
        gain < deficiency)
    if shortfall != expected_shortfall:
        fail("precodefail heavy shortfall is inconsistent")
    _parse_histogram(row["binary_def_hist"], deficiency,
                     "binary deficiency histogram")
    _parse_histogram(row["heavy_gain_hist"], gain,
                     "GF256 rank gain histogram")
    if (row["first_rank_fail"] != ("0" if rank_fail else "-1") or
            row["failure_trials"] != ("0" if weak else "")):
        fail("precodefail failure diagnostics are inconsistent")
    for field in (
            "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
            "residual_ms_mu", "backsub_ms_mu"):
        _parse_timing_decimal(row[field], field)
    if (HEX64.fullmatch(row["effective_precode_seed"]) is None or
            HEX32.fullmatch(row["effective_packet_seed"]) is None):
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


DETERMINISTIC_FIELDS = (
    "mix_pair", "delta", "precode_attempt", "packet_attempt",
    "coordinate_ordinal", "K", "attempt", "root_index",
    "loss_root", "schedule_index", "schedule", "cell_ordinal",
    "baseline_record_ordinal", "block_bytes", "loss_ppm", "overhead",
    "dense_anchor_layout", "mix_count", "binary_dense_rows",
    "gf256_heavy_rows", "heavy_family", "construction_seed_basis",
    "full_payload_solve", "overhead_stream", "cold_solve_wide_xor",
    "returncode", "staircase", "source_hits", "outcome", "success",
    "rank_fail", "error", "configuration_failure", "weak",
    "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
    "heavy_shortfall", "first_rank_fail", "failure_trials",
    "effective_precode_seed", "effective_packet_seed", "block_xors",
    "gf256_muladds", "configuration_stderr_sha256", "command_sha256",
    "bench_binary_sha256", "bench_source_git_commit",
    "controller_source_git_commit", "source_receipt_sha256",
    "pair01_baseline_record_sha256", "pair01_baseline",
    "baseline_matrix_sha256",
)


def deterministic_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in DETERMINISTIC_FIELDS if field not in record]
    if missing:
        fail("Stage 1 record omits deterministic fields: {}".format(missing))
    return {field: record[field] for field in DETERMINISTIC_FIELDS}


def parse_process_result(result: ProcessResult, controller_source_commit: str,
                         bench_source_commit: str, bench_sha256: str,
                         source_receipt_sha256: str,
                         baseline: Mapping[str, Any]) -> Mapping[str, Any]:
    if (COMMIT.fullmatch(controller_source_commit) is None or
            bench_source_commit != EXPECTED_STAGE0_SOURCE_COMMIT or
            bench_sha256 != EXPECTED_STAGE0_BENCH_SHA256 or
            SHA256.fullmatch(source_receipt_sha256) is None):
        fail("Stage 1 receipt identity is malformed")
    coordinate = result.invocation.coordinate
    if (baseline.get("record_ordinal") !=
            coordinate.baseline_record_ordinal or
            baseline.get("K") != coordinate.K or
            baseline.get("attempt") != coordinate.attempt or
            baseline.get("root_index") != coordinate.root_index or
            baseline.get("schedule") != coordinate.schedule):
        fail("pair01 baseline row does not match the candidate coordinate")
    lines = _framed_lines(result, bench_source_commit)
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
        "bench_binary_sha256": bench_sha256,
        "bench_source_git_commit": bench_source_commit,
        "controller_source_git_commit": controller_source_commit,
        "source_receipt_sha256": source_receipt_sha256,
        "pair01_baseline_record_sha256": sha256_json(baseline),
        "pair01_baseline": _baseline_projection(baseline),
        "baseline_matrix_sha256": BASELINE_FILES[
            "attempt-crossfit-matrix.jsonl"][1],
    }
    if result.returncode == 0:
        if result.stderr or len(lines) != 3:
            fail("successful precodefail process framing changed")
        values = _parse_csv_line(lines[2], len(CSV_HEADER))
        terminal = _parse_terminal_row(
            dict(zip(CSV_HEADER, values)), result.invocation)
    elif result.returncode == 2:
        expected_stderr = (
            "precodefail configuration failed N={} bb=2 "
            "heavy_family=periodic mix_count=2 attempt_mode=exact "
            "precode_attempt={} packet_attempt={} result=2\n".format(
                coordinate.K, coordinate.attempt,
                coordinate.packet_attempt)).encode("ascii")
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


def _validate_candidate_seed_receipt(
        record: Mapping[str, Any], baseline: Mapping[str, Any],
        coordinate_ordinal: int) -> None:
    coordinate = Coordinate(coordinate_ordinal)
    if record.get("configuration_failure") is not False:
        fail("fixed precode p failed before packet offset q was used")
    if record.get("effective_precode_seed") != \
            baseline.get("effective_precode_seed"):
        fail("candidate precode seed differs from diagonal pair01 p")
    expected_packet_seed = stage0._packet_seed_for_offset(
        baseline.get("effective_packet_seed"), coordinate.attempt, PACKET_DELTA)
    if record.get("effective_packet_seed") != expected_packet_seed:
        fail("candidate packet seed differs from wrapped delta2 q")


def _execute_roster(pinned: screen.PinnedExecutable,
                    controller_source_commit: str,
                    bench_source_commit: str, bench_sha256: str,
                    source_receipt_sha256: str,
                    baseline: Sequence[Mapping[str, Any]]) \
        -> Tuple[Mapping[str, Any], ...]:
    if len(baseline) != COORDINATE_COUNT:
        fail("candidate execution received an incomplete pair01 baseline")
    records: List[Mapping[str, Any]] = []
    for invocation in make_invocations():
        result = _run_raw(invocation, pinned)
        baseline_row = baseline[invocation.coordinate_ordinal]
        record = parse_process_result(
            result, controller_source_commit, bench_source_commit,
            bench_sha256, source_receipt_sha256, baseline_row)
        _validate_candidate_seed_receipt(
            record, baseline_row, invocation.coordinate_ordinal)
        records.append(record)
    if len(records) != OBSERVATION_COUNT:
        fail("Stage 1 executed the wrong number of observations")
    return tuple(records)


def _ratio_text(numerator: int, denominator: int) -> str:
    if type(numerator) is not int or type(denominator) is not int or \
            numerator < 0 or denominator <= 0:
        fail("aggregate ratio has an invalid integer operand")
    with localcontext() as context:
        context.prec = 60
        return format(Decimal(numerator) / Decimal(denominator), ".12f")


def _coordinate_result(coordinate_ordinal: int,
                       candidates: Sequence[Mapping[str, Any]],
                       baseline: Mapping[str, Any]) -> Mapping[str, Any]:
    if len(candidates) != len(REPETITIONS):
        fail("candidate coordinate receipt has the wrong repetition count")
    coordinate = Coordinate(coordinate_ordinal)
    return {
        **coordinate.identity(),
        "packet_delta": PACKET_DELTA,
        "candidate_precode_attempt": coordinate.attempt,
        "candidate_packet_attempt": coordinate.packet_attempt,
        "candidate_repetition_outcomes": [
            candidate["outcome"] for candidate in candidates],
        "baseline_outcome": baseline["outcome"],
        "candidate_success_all_repetitions": all(
            candidate["success"] for candidate in candidates),
        "baseline_success": baseline["success"],
    }


def adjudicate(records: Sequence[Mapping[str, Any]],
               baseline: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    if (len(records) != OBSERVATION_COUNT or
            len(baseline) != COORDINATE_COUNT):
        fail("Stage 1 result or baseline roster has the wrong cardinality")
    for coordinate_ordinal, baseline_row in enumerate(baseline):
        coordinate = Coordinate(coordinate_ordinal)
        if (baseline_row.get("record_ordinal") !=
                coordinate.baseline_record_ordinal or
                any(baseline_row.get(field) != value
                    for field, value in coordinate.identity().items()
                    if field not in (
                        "coordinate_ordinal", "baseline_record_ordinal"))):
            fail("Stage 1 baseline roster is noncanonical")
    by_key: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    controller_source_commit: Optional[str] = None
    source_receipt_sha256: Optional[str] = None
    for expected_ordinal, record in enumerate(records):
        if type(record) is not dict:
            fail("Stage 1 record is not an object")
        if expected_ordinal == 0:
            controller_source_commit = record.get(
                "controller_source_git_commit")
            source_receipt_sha256 = record.get("source_receipt_sha256")
            if (type(controller_source_commit) is not str or
                    COMMIT.fullmatch(controller_source_commit) is None or
                    type(source_receipt_sha256) is not str or
                    SHA256.fullmatch(source_receipt_sha256) is None):
                fail("Stage 1 controller source receipt is malformed")
        key = (record.get("coordinate_ordinal"), record.get("repetition"))
        coordinate = expected_ordinal // 2
        coordinate_identity = Coordinate(coordinate).identity()
        baseline_row = baseline[coordinate]
        if (record.get("ordinal") != expected_ordinal or key in by_key or
                key != (expected_ordinal // 2, expected_ordinal % 2) or
                any(record.get(field) != value for field, value in
                    coordinate_identity.items()) or
                record.get("schema") != RECORD_SCHEMA or
                record.get("mix_pair") != CANDIDATE_PAIR or
                record.get("delta") != PACKET_DELTA or
                record.get("precode_attempt") !=
                    coordinate_identity["attempt"] or
                record.get("packet_attempt") !=
                    ((coordinate_identity["attempt"] + PACKET_DELTA) & 255) or
                record.get("block_bytes") != BLOCK_BYTES or
                record.get("loss_ppm") != LOSS_PPM or
                record.get("overhead") != OVERHEAD or
                record.get("dense_anchor_layout") != "two07" or
                record.get("mix_count") != 2 or
                record.get("binary_dense_rows") != BINARY_DENSE_ROWS or
                record.get("gf256_heavy_rows") != GF256_HEAVY_ROWS or
                record.get("heavy_family") != "periodic" or
                record.get("construction_seed_basis") !=
                    "production-profile" or
                record.get("full_payload_solve") is not True or
                record.get("overhead_stream") != "paired" or
                record.get("cold_solve_wide_xor") != "policy" or
                record.get("promotion_evidence") is not False or
                record.get("fresh_roots_used") is not False or
                record.get("timing_evidence") is not False or
                type(record.get("success")) is not bool or
                type(record.get("rank_fail")) is not bool or
                type(record.get("error")) is not bool or
                type(record.get("configuration_failure")) is not bool or
                record.get("bench_binary_sha256") !=
                    EXPECTED_STAGE0_BENCH_SHA256 or
                record.get("bench_source_git_commit") !=
                    EXPECTED_STAGE0_SOURCE_COMMIT or
                record.get("controller_source_git_commit") !=
                    controller_source_commit or
                record.get("source_receipt_sha256") !=
                    source_receipt_sha256 or
                record.get("baseline_matrix_sha256") != BASELINE_FILES[
                    "attempt-crossfit-matrix.jsonl"][1] or
                record.get("baseline_record_ordinal") !=
                    baseline_row["record_ordinal"] or
                record.get("pair01_baseline_record_sha256") !=
                    sha256_json(baseline_row) or
                record.get("pair01_baseline") !=
                    _baseline_projection(baseline_row) or
                record.get("deterministic_projection_sha256") !=
                    sha256_json(deterministic_projection(record))):
            fail("Stage 1 result roster or deterministic receipt changed")
        _validate_candidate_seed_receipt(record, baseline_row, coordinate)
        by_key[key] = record
    if len(by_key) != OBSERVATION_COUNT:
        fail("Stage 1 result roster is incomplete")

    deduplicated: List[Mapping[str, Any]] = []
    repeat_exact_count = 0
    repeat_work_exact_count = 0
    for coordinate in range(COORDINATE_COUNT):
        first = by_key[(coordinate, 0)]
        second = by_key[(coordinate, 1)]
        exact = deterministic_projection(first) == \
            deterministic_projection(second)
        work_exact = (
            type(first["block_xors"]) is int and
            type(first["gf256_muladds"]) is int and
            type(first["inactivated_columns"]) is int and
            first["block_xors"] == second["block_xors"] and
            first["gf256_muladds"] == second["gf256_muladds"] and
            first["inactivated_columns"] == second["inactivated_columns"])
        repeat_exact_count += int(exact)
        repeat_work_exact_count += int(work_exact)
        deduplicated.append(first)

    candidate_weak = []
    baseline_weak = []
    repairs = []
    introductions = []
    for coordinate, (candidate_row, baseline_row) in enumerate(
            zip(deduplicated, baseline)):
        candidate_rows = (
            by_key[(coordinate, 0)], by_key[(coordinate, 1)])
        candidate_success = all(row["success"] for row in candidate_rows)
        receipt = _coordinate_result(
            coordinate, candidate_rows, baseline_row)
        if not candidate_success:
            candidate_weak.append(receipt)
        if not baseline_row["success"]:
            baseline_weak.append(receipt)
        if not baseline_row["success"] and candidate_success:
            repairs.append(receipt)
        if baseline_row["success"] and not candidate_success:
            introductions.append(receipt)

    errors = sum(int(row["error"]) for row in records)
    configuration_failures = sum(
        int(row["configuration_failure"]) for row in records)
    baseline_sums = {
        field: sum(row[field] for row in baseline)
        for field in EXPECTED_BASELINE_SUMS
    }
    candidate_sums = {
        field: (sum(row[field] for row in deduplicated)
                if all(type(row[field]) is int for row in deduplicated)
                else None)
        for field in EXPECTED_BASELINE_SUMS
    }
    ratios = {
        field: (_ratio_text(candidate_sums[field], baseline_sums[field])
                if type(candidate_sums[field]) is int else None)
        for field in EXPECTED_BASELINE_SUMS
    }
    per_K = []
    coordinates_per_K = len(ROOTS) * len(SCHEDULES)
    for K_index, K in enumerate(K_VALUES):
        begin = K_index * coordinates_per_K
        end = begin + coordinates_per_K
        metrics = {}
        for field in EXPECTED_BASELINE_SUMS:
            baseline_value = sum(row[field] for row in baseline[begin:end])
            candidate_value = sum(
                row[field] for row in deduplicated[begin:end])
            metrics[field] = {
                "pair01_baseline": baseline_value,
                "pair01_delta2_candidate": candidate_value,
                "candidate_minus_baseline": candidate_value - baseline_value,
            }
        per_K.append({
            "K": K,
            "coordinate_count": coordinates_per_K,
            "candidate_observation_count":
                coordinates_per_K * len(REPETITIONS),
            "metrics": metrics,
        })
    gates = {
        "candidate_weak_coordinates_zero": len(candidate_weak) == 0,
        "candidate_errors_zero": errors == 0,
        "candidate_configuration_failures_zero":
            configuration_failures == 0,
        "introductions_vs_pair01_zero": len(introductions) == 0,
        "all_pair01_weak_coordinates_repaired":
            len(repairs) == len(baseline_weak),
        "full_repeat_deterministic_projection_exact":
            repeat_exact_count == COORDINATE_COUNT,
        "full_repeat_work_exact":
            repeat_work_exact_count == COORDINATE_COUNT,
        "aggregate_block_xors_not_above_pair01":
            (type(candidate_sums["block_xors"]) is int and
             candidate_sums["block_xors"] <= baseline_sums["block_xors"]),
        "aggregate_gf256_muladds_not_above_pair01":
            (type(candidate_sums["gf256_muladds"]) is int and
             candidate_sums["gf256_muladds"] <=
                baseline_sums["gf256_muladds"]),
        "aggregate_inactivated_columns_not_above_pair01":
            (type(candidate_sums["inactivated_columns"]) is int and
             candidate_sums["inactivated_columns"] <=
                baseline_sums["inactivated_columns"]),
    }
    return {
        "candidate_pair": CANDIDATE_PAIR,
        "candidate_packet_delta": PACKET_DELTA,
        "candidate_observation_count": len(records),
        "deduplicated_candidate_coordinate_count": len(deduplicated),
        "repeat_exact_coordinate_count": repeat_exact_count,
        "repeat_work_exact_coordinate_count": repeat_work_exact_count,
        "candidate_errors": errors,
        "candidate_configuration_failures": configuration_failures,
        "baseline_weak_coordinates": baseline_weak,
        "candidate_weak_coordinates": candidate_weak,
        "repaired_coordinates": repairs,
        "introduced_weak_coordinates": introductions,
        "aggregates": {
            "pair01_baseline": baseline_sums,
            "pair01_delta2_candidate": candidate_sums,
            "candidate_to_baseline_ratios": ratios,
            "deduplication": (
                "one candidate value per exact repeated coordinate"),
        },
        "per_K_work": per_K,
        "gates": gates,
        "disposition": "PASS" if all(gates.values()) else "REJECT",
    }


def _source_receipt() -> Mapping[str, Any]:
    for item in SOURCE_PATHS:
        if not (REPO_ROOT / item).is_file():
            fail("required Stage 1 source is missing: {}".format(item))
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
        fail("cannot attest Stage 1 source checkout: {}".format(exc))
    if (head_process.returncode != 0 or head_process.stderr or
            status_process.returncode != 0 or status_process.stderr or
            tracked_process.returncode != 0 or tracked_process.stderr):
        fail("cannot attest Stage 1 source checkout")
    try:
        head = head_process.stdout.decode("ascii").strip()
    except UnicodeDecodeError:
        fail("Stage 1 source HEAD is not ASCII")
    if COMMIT.fullmatch(head) is None or status_process.stdout:
        fail("Stage 1 linked sources are not clean at a full source HEAD")
    try:
        tracked = tracked_process.stdout.decode("ascii").split("\0")
    except UnicodeDecodeError:
        fail("Stage 1 tracked source roster is not ASCII")
    if tracked and tracked[-1] == "":
        tracked.pop()
    if len(tracked) != len(SOURCE_PATHS) or set(tracked) != set(SOURCE_PATHS):
        fail("Stage 1 source receipt includes a path not tracked at HEAD")
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
            fail("Stage 1 source {} differs byte-for-byte from HEAD".format(
                item))
    try:
        final_head_process = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot re-attest Stage 1 source HEAD: {}".format(exc))
    if (final_head_process.returncode != 0 or final_head_process.stderr or
            final_head_process.stdout != (head + "\n").encode("ascii")):
        fail("Stage 1 source HEAD changed while the receipt was built")
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


def _write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) \
        -> Tuple[str, int]:
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


def run_stage1(bench: Path, baseline_dir: Path, v9_dir: Path,
               stage0_dir: Path, output_dir: Path) -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path must be one fresh path below an existing directory")
    source = _source_receipt()
    # All historical validation, including the complete 276,480-row matrix
    # parse and the 540-observation V9 cross-check, precedes candidate work.
    with load_historical_evidence(
            baseline_dir, v9_dir, stage0_dir) as historical:
        historical_receipt = historical.receipt()
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            return _run_stage1_pinned(
                pinned, output_dir, source, contract, historical,
                historical_receipt)


def _run_stage1_pinned(
        pinned: screen.PinnedExecutable, output_dir: Path,
        source: Mapping[str, Any], contract: Mapping[str, Any],
        historical: HistoricalEvidence,
        historical_receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    controller_source_commit = source["source_git_commit"]
    source_receipt_sha256 = source["source_receipt_sha256"]
    bench_receipt = pinned.receipt()
    if (bench_receipt.get("sha256") != EXPECTED_STAGE0_BENCH_SHA256 or
            bench_receipt.get("size") != EXPECTED_STAGE0_BENCH_SIZE):
        fail("Stage 1 requires the exact authenticated Stage-0 benchmark")
    # Deliberately the final action before the first candidate invocation.
    description = screen.read_bench_description(
        pinned.path, EXPECTED_STAGE0_SOURCE_COMMIT, pinned.descriptor)
    records = _execute_roster(
        pinned, controller_source_commit, EXPECTED_STAGE0_SOURCE_COMMIT,
        bench_receipt["sha256"], source_receipt_sha256,
        historical.baseline)
    verdict = adjudicate(records, historical.baseline)

    temporary = Path(tempfile.mkdtemp(
        prefix=".wh2-mix2-packet-offset-stage1-", dir=str(output_dir.parent)))
    published = False
    try:
        records_sha256, record_count = _write_jsonl(
            temporary / RECORD_NAME, records)
        if record_count != OBSERVATION_COUNT:
            fail("Stage 1 JSONL cardinality changed")
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                historical.receipt() != historical_receipt):
            fail("Stage 1 executable, source, or history changed during run")
        summary: Dict[str, Any] = {
            "schema": SUMMARY_SCHEMA,
            "contract_sha256": contract["contract_sha256"],
            "source_receipt": source,
            "bench": bench_receipt,
            "bench_description": description,
            "historical_evidence": historical_receipt,
            "records_file": RECORD_NAME,
            "records_sha256": records_sha256,
            "record_count": record_count,
            **verdict,
            "stage1_only": True,
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
            fail("Stage 1 inputs changed before publication")
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
    subparsers.add_parser("describe", help="print the frozen Stage 1 contract")
    run_parser = subparsers.add_parser(
        "run", help="run the 2,160-observation consumed-root replay")
    run_parser.add_argument("--bench", type=Path, required=True)
    run_parser.add_argument("--baseline-dir", type=Path, required=True)
    run_parser.add_argument("--v9-dir", type=Path, required=True)
    run_parser.add_argument("--stage0-dir", type=Path, required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        _validate_constants()
        if arguments.command == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run_stage1(
            arguments.bench, arguments.baseline_dir, arguments.v9_dir,
            arguments.stage0_dir, arguments.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (Stage1Error, stage0.Stage0Error, crossfit.CrossfitError,
            screen.ShortScreenError) as exc:
        print("wh2 MIX2 packet-offset Stage 1: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
