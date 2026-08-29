#!/usr/bin/env python3
"""Fail-closed one-shot controller and adjudicator for the WH2 direct-systematic complement screen."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
import os
from pathlib import Path
import re
import select
import selectors
import signal
import stat
import subprocess
import sys
import tempfile
import time
import types
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


CONFIG_SCHEMA = "wirehair.wh2.direct-systematic-complement-screen.config.v2"
PANEL_SCHEMA = "wirehair.wh2.direct-systematic-complement-screen.panel.v2"
TERMINAL_SCHEMA = "wirehair.wh2.direct-systematic-complement-screen.terminal.v2"
SUMMARY_SCHEMA = "wirehair.wh2.direct-systematic-complement-screen.summary.v2"
DESCRIPTOR_SCHEMA = "wirehair.wh2.direct-systematic-complement-arm.v2"

CELLS: Tuple[Tuple[int, int], ...] = (
    (32, 1280),
    (100, 64),
    (100, 1280),
    (1000, 1280),
    (2048, 64),
    (8192, 64),
    (8192, 1280),
    (20000, 64),
    (20000, 1280),
    (32768, 64),
    (32768, 1280),
)
PANEL_SCOPES = ("fresh-systematic-total-v2",)
METRICS: Tuple[Tuple[str, str, str], ...] = (
    ("fresh-systematic-total-v2", "fresh-systematic-total-v2", "elapsed_ns"),
)
COMPARISONS = ("baseline-aa", "candidate-aa", "candidate-over-baseline")
REPLICATES = 12
INVOCATION_BUDGET = 65536
MINIMUM_INVOCATIONS = 24
INTERNAL_DEADLINE_SECONDS = 105
CONTROLLER_DEADLINE_SECONDS = 110.0
OUTER_DEADLINE_SECONDS = 120.0
EXPECTED_PANEL_COUNT = (
    REPLICATES * len(CELLS) * len(PANEL_SCOPES) * len(COMPARISONS)
)
EXPECTED_RECORD_COUNT = EXPECTED_PANEL_COUNT + 2
EXPECTED_TIMED_INVOCATIONS = 43224
EXPECTED_WARMUP_INVOCATIONS = 792
EXPECTED_TOTAL_INVOCATIONS = 44016
T11_975 = 2.200985160082949
MAX_INT63 = (1 << 63) - 1
MAX_UINT32 = (1 << 32) - 1
MAX_UINT64 = (1 << 64) - 1
MAX_STDOUT_BYTES = 16 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
MAX_FINAL_ARTIFACT_BYTES = MAX_STDOUT_BYTES + 1
MAX_BINARY_BYTES = 64 * 1024 * 1024
MAX_SAMPLER_CSV_BYTES = 8 * 1024 * 1024
MAX_SAMPLER_VALIDATION_BYTES = 4 * 1024 * 1024
MAX_SAMPLER_RECEIPT_BYTES = 1024 * 1024
MAX_THERMAL_WINDOW_BYTES = 2 * 1024 * 1024
MAX_SOURCE_FILE_BYTES = 1024 * 1024
MAX_PROC_VECTOR_BYTES = 64 * 1024
MAX_PREFLIGHT_RECEIPT_BYTES = 4 * 1024 * 1024
SAMPLER_HEARTBEAT_MAX_GAP_NS = 5_000_000_000
MAX_FAILURE_TEXT_CHARS = 64 * 1024
MAX_FAILED_GATES = 64
MAX_PUBLICATION_FAILURES = 64
EXPECTED_TARGET_CPU = 120
EXPECTED_CONTROLLER_CPU = 121
EXPECTED_SAMPLER_CPU = 122
EXPECTED_CAMPAIGN_UID = 1000
EXPECTED_CAMPAIGN_GID = 1000
EXPECTED_SAMPLER_I2C_GID = 113
EXPECTED_TARGET_CORE = (0, 56)
EXPECTED_TARGET_THREADS = (56, 120)
EXPECTED_CONTROLLER_CORE = (0, 57)
EXPECTED_SAMPLER_CORE = (0, 58)
SIBLING_NON_IDLE_TICK_CAP = 5
THERMAL_MAX_MILLIC = 85000
HEALTH_SCHEMA = "wirehair.wh2.direct-systematic-complement-health.v3"
HEALTH_LOADER_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-health-source-loader.v2"
)
CONTROLLER_SCHEMA = "wirehair.wh2.direct-systematic-complement-controller.v2"
CONTROLLER_PROVENANCE_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-controller-provenance.v3"
)
CLAIM_SCHEMA = "wirehair.wh2.direct-systematic-complement-claim.v3"
COMPLETE_SCHEMA = "wirehair.wh2.direct-systematic-complement-complete.v2"
VERIFY_RESULT_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-retained-verification.v2"
)
PREFLIGHT_SEAL_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-preflight-seal.v1"
)
PREFLIGHT_CONTROLLER_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-preflight-controller.v1"
)
FIXED_OUTPUT_DIR = Path(
    "/var/tmp/wh2-retained-direct-systematic-complement-v2-r0"
)
FIXED_CLAIM_PATH = Path(str(FIXED_OUTPUT_DIR) + ".claim")
GIT_EXECUTABLE = Path("/usr/bin/git")
SEALED_LAUNCH_ENVIRONMENT = {
    "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
    "PATH": "/usr/bin:/bin", "TZ": "UTC",
}
SAMPLER_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-sampler-attestation.v3"
)
SAMPLER_TERMINAL_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-sampler-terminal.v1"
)
PROCESS_SECURITY_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-process-security.v1"
)
THERMAL_VALIDATION_STREAM_SCHEMA = "wirehair.wh2.thermal_validation_stream.v1"
THERMAL_VALIDATION_SAMPLE_SCHEMA = "wirehair.wh2.thermal_validation_sample.v1"
THERMAL_SAMPLER_SCHEMA = "wirehair.wh2.thermal_sampler.v2"
EXPECTED_SAMPLER_PYTHON = "/usr/bin/python3"
EXPECTED_SAMPLER_PYTHON_FLAGS = ("-I", "-S", "-B")
EXPECTED_SAMPLER_INTERVAL_TEXT = "1.0"
EXPECTED_SAMPLER_DIMM_ATTEMPTS_TEXT = "5"
EXPECTED_SAMPLER_DIMM_RETRY_DELAY_TEXT = "0.01"
THERMAL_SAMPLER_THRESHOLDS = {
    "dimm_safety_c_inclusive": 90.0,
    "hot_confirm_samples": 3,
    "max_dimm_jump_c": 12.0,
    "max_dimm_rate_c_per_s": 6.0,
    "max_plausible_dimm_c_exclusive": 130.0,
    "min_plausible_dimm_c_exclusive": 0.0,
    "telemetry_fault_abort_samples": 8,
}
CPU_TICK_RECEIPT_KEYS = {
    "cpu", "non_idle_ticks", "read_monotonic_ns", "tick_fields",
}
CPU_TICK_FIELD_KEYS = {
    "idle", "iowait", "irq", "nice", "softirq", "steal", "system", "user",
}
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
LOWER40 = re.compile(r"^[0-9a-f]{40}$")

# Frozen from the reviewed native complement receipt table.  Exact roster and
# field checks below keep any later accidental omission fail closed.
EXPECTED_CELL_RECEIPTS: Mapping[Tuple[int, int], Mapping[str, str]] = {
    (32, 1280): {
        "direct_hit_oracle_sha256": "57fe8cc94797c30786446d7d21feefd78523fe0a42f55106532530a54956ab68",
        "equation_configuration_sha256": "b283539ee89370374952a920b7f4605de7930fa22a428d8edac7d40f9f305971",
        "first_repair_sha256": "e4ec3da26bc910aab4a7f1d5143496526e7b040a019aef5bb30a98e6cbde7e11",
        "full_repair_sha256": "da9ba77df030059b0e68f78d544590b5b0cbf5faaa673356d63c6d391cefa062",
        "high_id_oracle_sha256": "006ba6b059f6b86dc5a92d7b5c53b8dcdf29e4d6e9a83276ea4e6397606b4fcb",
        "solved_state_equivalence_receipt_sha256": "4cf3c08f22a4fd9055c894fd8c7d253c746c1089d1483d5a84f112c69b7b47da",
        "source_seed": "0x547000f8f0f19096",
        "source_sha256": "1e296fa19681d28489ea4c245f975e8f61c63085c9ef19897d4ff7f1b2315d71",
    },
    (100, 64): {
        "direct_hit_oracle_sha256": "c49f184f939634d268b9ac30f452e172a6acf6cfcd5952ee51fc9c58b0b32dfc",
        "equation_configuration_sha256": "58b35a698480884fb8632a5be4ecd2b12252fed394eaa410d439258c2db61ba4",
        "first_repair_sha256": "91ab0bb101199ff7143cc8ac3e6ee460fce3cd3264d1678b8824c2cf45fc3a3a",
        "full_repair_sha256": "d5ce6769b52e553aabf28f17e65628250c0c6d7ead7a83e46ec2879131060325",
        "high_id_oracle_sha256": "072925662f69d22a448aeffb7d9ec48d05ec3a8eaf3e96b85e3f80944ccb5ad5",
        "solved_state_equivalence_receipt_sha256": "dbb815ef7b18adb456c552044eb1e3843830dddb310667e1694de2792faf1ab4",
        "source_seed": "0x547000f8f07995d6",
        "source_sha256": "a2e57558b1c8027ae2a5c4a0634903a105a0f89935bb72748c503e55a65bd655",
    },
    (100, 1280): {
        "direct_hit_oracle_sha256": "c49f184f939634d268b9ac30f452e172a6acf6cfcd5952ee51fc9c58b0b32dfc",
        "equation_configuration_sha256": "07cbf4890e7b9e9650c6d16db922550f562fef1ad7d85cc91ce105edd546e215",
        "first_repair_sha256": "2032d39d529b666847e477d34e3727606f3fa1261b2e12d5ffe0cd9caac72b97",
        "full_repair_sha256": "04c2233417e2797c59bdd22d05344d0128b3c7120d9c43f20e694e00dbf5888f",
        "high_id_oracle_sha256": "b8af8ef80d2f64a62503afcfd4aa1a390f2f8270f663a5a1a7955a001f1a4083",
        "solved_state_equivalence_receipt_sha256": "ee95c5f56c1ac3bc5e749884e8e0639dd274a32691d48328d07cf6e81e19f2ee",
        "source_seed": "0x547000f8f0799096",
        "source_sha256": "17f6f7bc90fb70cd21e0968809b57aa627991388537105ce5749df3a7d11736c",
    },
    (1000, 1280): {
        "direct_hit_oracle_sha256": "a6f9cb15de11c8cb55c30f134843006627d11bfa4a196886392c8a9616801139",
        "equation_configuration_sha256": "92cf189ef939a4aba39c68a701aa0378ead1a04c6714c613c366e67847bec4ed",
        "first_repair_sha256": "695f23f5df60b3bd7f04225f53c4963fe2af815aede65128701791a8661680ea",
        "full_repair_sha256": "857cbb1384675e86abd40a7fca1560892e3f191b40c1d03f377d9b68e85f929d",
        "high_id_oracle_sha256": "48d051490a1702234fba2216ed74c249e834c60e37f5db60e800e2b334f1f2a8",
        "solved_state_equivalence_receipt_sha256": "4f4a4af87c304a13810f251582c3cbceeab0cc704ac07c29da9f1cd31af2c7ac",
        "source_seed": "0x547000f8f7619096",
        "source_sha256": "ff645fe4936b314fd4c471c626e7b2e020e533ae1206e1875d2d46924f16d7ca",
    },
    (2048, 64): {
        "direct_hit_oracle_sha256": "88364e26249a7f897da55030719c17543675f4ddf040dc5239cd75fc82389795",
        "equation_configuration_sha256": "53b2410721b558e8705c5512c74268233cff1353220a8981ba0f89f594dec147",
        "first_repair_sha256": "96db9bc74ffe6ef938c84fa678361f4a5c7afa10b799e74a3ec312a386f056c5",
        "full_repair_sha256": "c039209bc65ab8557abace8e1a4931dcde135fb19e54c7f5f532a6bb95752be1",
        "high_id_oracle_sha256": "817a3c100489f2e8a838faedb71facbb4d23558c70c090b1d5676e1db3ec3e83",
        "solved_state_equivalence_receipt_sha256": "69c492fa7e3b5c3332d7b92273f06ed7928fbd0ba4cc354fc799e734f9dbb12f",
        "source_seed": "0x547000f8e0b195d6",
        "source_sha256": "9c250263d5b2e75f8f82d28a3abe72cc644e0347cc014caed4164f1e6ee0724a",
    },
    (8192, 64): {
        "direct_hit_oracle_sha256": "61ed57e6b1b924da3ed34a053a840beca82faf648b958f85eb849001d2ad7181",
        "equation_configuration_sha256": "14ce58ca623b635a98c9aa802edfb0a4f24095d1cedfc9069e5fdca88b16daf2",
        "first_repair_sha256": "7d6337d1b55f1e59715df8c825a82288178f1432361cfaec49241f15da257179",
        "full_repair_sha256": "9aa31b1e1feb78010cf81fbb3c2127b7fdaa469612a78bb9ed4f8bc347517800",
        "high_id_oracle_sha256": "b566c35130def561985c04e42607753600035d69c462e88b1bcfd5058df63192",
        "solved_state_equivalence_receipt_sha256": "34d482bff01ee09cb458bb912b490161cf1ffa723bbf3fc9cd4784de3de32862",
        "source_seed": "0x547000f8b0b195d6",
        "source_sha256": "4c1a75ed90a9df9aee57aec7649656147ce8d2ea4b8cd75f0bd6bfc46c062952",
    },
    (8192, 1280): {
        "direct_hit_oracle_sha256": "61ed57e6b1b924da3ed34a053a840beca82faf648b958f85eb849001d2ad7181",
        "equation_configuration_sha256": "c3c69374cec038a6f38127707e78cce35e9293f9f8e607625a00bda608b07511",
        "first_repair_sha256": "71a876791a4e443544e2730e2a9a1eb361f709b7c591eff980b437e9f0ec25b1",
        "full_repair_sha256": "d4917f59db9e3a57baa26edd9f0f14f7aeead9372e6d5678749efb565c5bf307",
        "high_id_oracle_sha256": "4e3e79768942bfa306d5198840c44bd90e7a923dbdf2575bd4bc48d4ff0d6ae8",
        "solved_state_equivalence_receipt_sha256": "29448aec3574f98177cddbc61d4a8a3bb754359ffb0ea7faf351268fd34387c8",
        "source_seed": "0x547000f8b0b19096",
        "source_sha256": "4ca30eb14e3e1b516390933fc18a8e6773ee52b5265082d23ef42e9566ed673f",
    },
    (20000, 64): {
        "direct_hit_oracle_sha256": "4486b30e27d851608c266d583e55ebf9ccbabc706638a2a214c12499fb3074e7",
        "equation_configuration_sha256": "f877ec4bb14d2c19023620c58d7438c0b8c6ffeebe1a71671eef4e980fc93a79",
        "first_repair_sha256": "e64980ab7053b4b2acd93b2f33b5cba875429144c58b30c3e7642ce1e6a4cb64",
        "full_repair_sha256": "d701db59cab3dae48e8bc5ea33ec11b9d3afcf62b4266ef2e5b29e370a79c4d3",
        "high_id_oracle_sha256": "656b06de7ef753f0999d89f176b298510592a9574f6c0969029fb02a6675d0a6",
        "solved_state_equivalence_receipt_sha256": "92622eb136713cc4330ac7304c71e02783bcc0cf7079dcba8e66050c481d605e",
        "source_seed": "0x547000f86cf195d6",
        "source_sha256": "3db8c777ae42ad27e188d074cd30ce5f993f1a9e682f33d08d549f39d6b76d5b",
    },
    (20000, 1280): {
        "direct_hit_oracle_sha256": "4486b30e27d851608c266d583e55ebf9ccbabc706638a2a214c12499fb3074e7",
        "equation_configuration_sha256": "ccbfc46c15380255317dda067e6a8cb814c43fe2f28a6b347a804482b37cfb68",
        "first_repair_sha256": "dc8c2ed8ed9a8947ef73376222c23eba0af5da221748218bddbe22292898d169",
        "full_repair_sha256": "502588e72e1ac409690d4897a063ec81e5e0a45bcd3a32da47432e0cd27a4d7e",
        "high_id_oracle_sha256": "64d6afc42f2c61acee4384143191e5d432f5f7f3ee5967bebcee4aa06bc49923",
        "solved_state_equivalence_receipt_sha256": "1c32565b969f3e443d5ac1943fcf3f3dffbd1b1a24398f1462be244227bafd81",
        "source_seed": "0x547000f86cf19096",
        "source_sha256": "0b7004dd9c291f91a1aef3219c77862ca35d792ec2a6fc276695f3e47044d3a4",
    },
    (32768, 64): {
        "direct_hit_oracle_sha256": "4338dea9a6928832fbfb1e45e23f51313a13b4ad4bde7ccf87b9a9e0c3737565",
        "equation_configuration_sha256": "e361af19658f217f663436822d1df28b42fc7644622b2b384e5081934b710626",
        "first_repair_sha256": "bb73eca6a8929f39be371badaa3485fd45efbdb0303c2a94c9bdf16a1a3493ac",
        "full_repair_sha256": "80720243cc518ea1b559c9f23ea46b6fc87dd9a4f276350bcfca72f1d929d982",
        "high_id_oracle_sha256": "c7989e346ba0a3bb693ec5f909bcdc93545a3c43893c17af740d4923ee8fca45",
        "solved_state_equivalence_receipt_sha256": "5ed759119604192b2f8f038825f09783113b362218f6e59ac6c9389923af53af",
        "source_seed": "0x547000f9f0b195d6",
        "source_sha256": "3628fea7aada74f0dedf60f34ffcc7e76cb634b23d9f28738a1ab9a72461a15e",
    },
    (32768, 1280): {
        "direct_hit_oracle_sha256": "4338dea9a6928832fbfb1e45e23f51313a13b4ad4bde7ccf87b9a9e0c3737565",
        "equation_configuration_sha256": "4d9978d65eacfd17cb7b17d5dbc8cae28c0bfd989e4aed3f5f4ec9c24ee85fed",
        "first_repair_sha256": "40d3fe7e775ad2993747cafe0fcd8663f5ba91cf92db129355305533790f1903",
        "full_repair_sha256": "9cfec479624a860ecc6d3feffdb4f73f47dfac77a84ffce7287032b84cfa3769",
        "high_id_oracle_sha256": "9e6b0375d9b4c8d9ffa9c48e4d957962f5239a4916f718754add9a39d2c7594e",
        "solved_state_equivalence_receipt_sha256": "5e00fc2f973da725376fb6d606e4f04c0e467e692483f245cf1f7f6618788f49",
        "source_seed": "0x547000f9f0b19096",
        "source_sha256": "13aaa0e514bcdf73cc24645cfa67155df23fe645d8ba4161601779668f93f955",
    },
}

SAMPLER_SOURCE_PATH = "bench/wirehair_expo_thermal_sampler.py"
SOURCE_PATHS = (
    "CMakeLists.txt",
    "bench/Wh2DirectSystematicComplementLaunch.py",
    "bench/Wh2DirectSystematicComplementScreen.cpp",
    "bench/Wh2DirectSystematicComplementScreen.py",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "bench/Wh2NativePanel.cpp",
    "bench/Wh2NativePanel.h",
    "bench/wh2_benchmark_contract.py",
    "bench/wh2_benchmark_contract_v4.json",
    "bench/wh2_native_short_screen.py",
    "bench/wh2_run_native_short_screen.py",
    SAMPLER_SOURCE_PATH,
    "cmake/Wh2DirectSystematicComplementSymbolAudit.cmake",
    "cmake/Wh2TimingPolicySymbolAudit.cmake",
)
BUILD_PATHS = (
    "CMakeCache.txt",
    "compile_commands.json",
    "build.ninja",
    "CMakeFiles/rules.ninja",
)
FINAL_OUTPUT_NAMES = (
    "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
    "controller.json", "COMPLETE",
)
HEALTH_COLLECTION_DEADLINE_SECONDS = 115.0
SAMPLER_STALE_RECOVERY_WORK_SECONDS = 1.0
SAMPLER_HEALTH_FINALIZE_RESERVE_SECONDS = 0.5
FINAL_COMMIT_START_DEADLINE_SECONDS = 117.0
RUN_ONCE_OPTION_ORDER = (
    "--binary", "--build-dir", "--cpu", "--controller-cpu",
    "--sampler-pid", "--sampler-cpu", "--sampler-script", "--sampler-csv",
    "--sampler-pid-file", "--sampler-validation-jsonl", "--sampler-receipt",
    "--expected-source-commit", "--expected-binary-sha256",
    "--expected-binary-uid", "--expected-build-manifest-sha256",
    "--expected-controller-sha256", "--expected-controller-uid",
    "--expected-controller-gid",
    "--expected-git-sha256", "--expected-python-sha256",
    "--expected-sampler-process-start-ticks",
    "--expected-sampler-script-sha256",
    "--expected-sampler-csv-device", "--expected-sampler-csv-inode",
    "--expected-sampler-pid-file-device",
    "--expected-sampler-pid-file-inode",
    "--expected-sampler-validation-device",
    "--expected-sampler-validation-inode",
    "--expected-sampler-receipt-device",
    "--expected-sampler-receipt-inode",
    "--expected-sampler-cmdline-sha256", "--expected-sampler-environ-sha256",
    "--expected-sampler-executable-sha256", "--expected-sampler-uid",
    "--expected-sampler-gid", "--expected-sampler-i2c-gid",
    "--expected-source-manifest-sha256",
)
PREFLIGHT_SEAL_OPTION_ORDER = (
    "--binary", "--build-dir", "--sampler-pid", "--sampler-script",
    "--sampler-csv", "--sampler-pid-file", "--sampler-validation-jsonl",
    "--sampler-receipt", "--expected-source-commit",
)
PREFLIGHT_SEAL_DEADLINE_SECONDS = 30.0
PREFLIGHT_SIDE_EFFECT_GUARD = False


class ValidationError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise ValidationError(message)


def require_preflight_side_effects_allowed(where: str) -> None:
    if PREFLIGHT_SIDE_EFFECT_GUARD:
        fail("preflight side-effect tripwire: " + where)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob_oid(data: bytes) -> str:
    preimage = b"blob " + str(len(data)).encode("ascii") + b"\0" + data
    return hashlib.new(
        "sha1", preimage, usedforsecurity=False
    ).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def bounded_failure_text(
    parts: Iterable[str], limit: int = MAX_FAILURE_TEXT_CHARS,
) -> str:
    if type(limit) is not int or limit < 128 or limit > MAX_FAILURE_TEXT_CHARS:
        fail("failure-text limit differs")
    values = [part for part in parts if type(part) is str and part]
    text = "; ".join(values) or "none"
    if len(text) <= limit:
        return text
    trailer = "...[sha256:" + sha256_bytes(text.encode("utf-8")) + "]"
    return text[:limit - len(trailer)] + trailer


def exact_keys(value: Mapping[str, Any], keys: Iterable[str], where: str) -> None:
    expected = set(keys)
    actual = set(value.keys())
    if actual != expected:
        fail(f"{where} keys differ: expected {sorted(expected)}, got {sorted(actual)}")


def exact_int(value: Any, lower: int, upper: int, where: str) -> int:
    if type(value) is not int or value < lower or value > upper:
        fail(f"{where} is not an integer in [{lower},{upper}]")
    return value


def exact_string(value: Any, expected: str, where: str) -> None:
    if type(value) is not str or value != expected:
        fail(f"{where} differs")


def lower_hash(value: Any, where: str) -> str:
    if type(value) is not str or LOWER64.fullmatch(value) is None:
        fail(f"{where} is not a lowercase SHA-256")
    return value


def unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    fail(f"non-finite JSON number: {value}")


def parse_transcript(raw: bytes) -> List[Dict[str, Any]]:
    if not raw or not raw.endswith(b"\n") or b"\r" in raw:
        fail("raw transcript must be nonempty LF-terminated JSONL without CR")
    lines = raw.splitlines(keepends=True)
    if len(lines) != EXPECTED_PANEL_COUNT + 2:
        fail(f"raw record count is {len(lines)}, expected {EXPECTED_PANEL_COUNT + 2}")
    records: List[Dict[str, Any]] = []
    for index, line in enumerate(lines):
        if len(line) > 1024 * 1024:
            fail(f"raw line {index} exceeds 1 MiB")
        try:
            record = json.loads(
                line[:-1].decode("ascii"),
                object_pairs_hook=unique_object,
                parse_constant=reject_constant,
            )
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            fail(f"raw line {index} is malformed: {exc}")
        if type(record) is not dict:
            fail(f"raw line {index} is not an object")
        if canonical_bytes(record) + b"\n" != line:
            fail(f"raw line {index} is not canonical JSON")
        records.append(record)
    return records


def descriptor(metric: str, emission: str) -> Dict[str, Any]:
    if metric != "fresh-systematic-total-v2" or emission not in (
        "equation-eval-v1", "retained-source-direct-v1"
    ):
        fail("invalid descriptor request")
    return {
        "arm": "wirehair2_head",
        "codec": "wirehair2_certified",
        "equation_transform": "none",
        "metric": metric,
        "schema": DESCRIPTOR_SCHEMA,
        "source_acquisition": "fixture-copy-before-clock-move-v1",
        "source_storage": "native-arm-owned-kxb-v1",
        "systematic_emission": emission,
        "timed_work": (
            "fresh-eager-init-plus-systematic-ids-0-through-k-minus-1-v2"
        ),
    }


def descriptor_hash(value: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_bytes(value))


CONFIG_KEYS = {
    "campaign", "cells", "comparisons", "descriptors",
    "expected_invocations", "expected_measured_invocations",
    "expected_panels", "expected_records", "expected_warmup_invocations",
    "internal_deadline_seconds", "invocation_budget", "minimum_invocations",
    "panel_key_schema", "panel_replicates", "schema", "scope",
    "source_git_commit", "source_seed_base", "target_cpu",
}
CELL_KEYS = {
    "K", "block_bytes", "construction_equivalent",
    "direct_hit_oracle_sha256", "direct_hit_oracle_verified",
    "equation_configuration", "equation_configuration_sha256",
    "first_repair_sha256", "full_repair_oracle_verified",
    "full_repair_sha256", "high_id_oracle_sha256",
    "high_id_oracle_verified", "invocations_by_replicate",
    "solved_state_equivalence_receipt_sha256",
    "solved_state_equivalent", "solved_state_scope", "source_seed",
    "source_sha256", "systematic_oracle_sha256",
    "systematic_oracle_verified",
}
CONFIGURATION_KEYS = {
    "K", "block_bytes", "dense_anchor_layout", "dense_identity_corner",
    "dense_rows", "heavy_family", "heavy_rows", "mix_count", "packet_attempt",
    "packet_peel_seed", "precode_attempt", "precode_seed", "source_hits",
    "staircase",
}
COMPARISON_KEYS = {"left_mode", "name", "right_mode"}
DESCRIPTOR_ENTRY_KEYS = {"descriptor", "descriptor_sha256", "mode"}
PANEL_KEYS = {
    "K", "block_bytes", "comparison", "invocations_per_slot",
    "left_descriptor_sha256", "left_direct_systematic_packets",
    "left_outcome_code", "order", "panel_key_sha256", "primary_metric",
    "replicate", "right_descriptor_sha256",
    "right_direct_systematic_packets", "right_outcome_code", "schema",
    "scope", "slots", "target_cpu",
}
SLOT_KEYS = {"elapsed_ns", "invocation_count", "side"}
PANEL_KEY_SCHEMA = "wirehair.wh2.direct-systematic-complement-screen.panel-key.v2"
CAMPAIGN = "wh2-retained-direct-systematic-complement-v2-r0"
SCOPE = "fresh-systematic-total-v2"
GOLDEN_PANEL_KEY_SHA256 = (
    "a0a0933ae824bacc6d11eb00d14d63efa6a55ec4e2bfe5f2f1e0d590830e1b4e"
)
SOURCE_SEED_BASE = "0x547000f8f0b19596"
SOLVED_STATE_SCOPE = (
    "identity-system-rows-packet-runtime-configuration-intermediate-bytes-"
    "non-timing-solve-statistics-v1"
)
EXPECTED_COMPARISON_RECEIPTS = (
    {"left_mode": "equation", "name": "baseline-aa", "right_mode": "equation"},
    {"left_mode": "retained", "name": "candidate-aa", "right_mode": "retained"},
    {
        "left_mode": "retained", "name": "candidate-over-baseline",
        "right_mode": "equation",
    },
)
EXPECTED_RECEIPT_FIELDS = {
    "direct_hit_oracle_sha256", "equation_configuration_sha256",
    "first_repair_sha256", "full_repair_sha256", "high_id_oracle_sha256",
    "solved_state_equivalence_receipt_sha256", "source_seed", "source_sha256",
}


def expected_source_seed(k: int, block_bytes: int) -> str:
    value = 0x547000F8F0B19596 ^ (k << 17) ^ block_bytes
    return "0x" + format(value & MAX_UINT64, "x")


def direct_hit_oracle_preimage(cell: Mapping[str, Any]) -> Dict[str, Any]:
    k = cell["K"]
    return {
        "equation_direct_systematic_packets": 0,
        "retained_direct_systematic_packets": k,
        "verified_ids": [0, k - 1, k, MAX_UINT32],
    }


def solved_state_preimage(cell: Mapping[str, Any]) -> Dict[str, Any]:
    return {
        "configuration_sha256": cell["equation_configuration_sha256"],
        "direct_hit_oracle_sha256": cell["direct_hit_oracle_sha256"],
        "first_repair_sha256": cell["first_repair_sha256"],
        "full_repair_sha256": cell["full_repair_sha256"],
        "high_id_oracle_sha256": cell["high_id_oracle_sha256"],
        "source_sha256": cell["source_sha256"],
        "systematic_oracle_sha256": cell["systematic_oracle_sha256"],
        "verified_state_scope": SOLVED_STATE_SCOPE,
    }


def validate_config(
    config: Dict[str, Any], cpu: int, expected_commit: str
) -> Dict[Tuple[str, str], str]:
    exact_keys(config, CONFIG_KEYS, "config")
    exact_string(config["campaign"], CAMPAIGN, "config.campaign")
    exact_string(config["schema"], CONFIG_SCHEMA, "config.schema")
    exact_string(config["scope"], SCOPE, "config.scope")
    exact_string(
        config["panel_key_schema"], PANEL_KEY_SCHEMA,
        "config.panel_key_schema",
    )
    exact_string(config["source_seed_base"], SOURCE_SEED_BASE, "source seed base")
    exact_int(config["target_cpu"], cpu, cpu, "config.target_cpu")
    for field, expected in (
        ("internal_deadline_seconds", INTERNAL_DEADLINE_SECONDS),
        ("invocation_budget", INVOCATION_BUDGET),
        ("minimum_invocations", MINIMUM_INVOCATIONS),
        ("panel_replicates", REPLICATES),
        ("expected_panels", EXPECTED_PANEL_COUNT),
        ("expected_records", EXPECTED_RECORD_COUNT),
        ("expected_invocations", EXPECTED_TOTAL_INVOCATIONS),
        ("expected_measured_invocations", EXPECTED_TIMED_INVOCATIONS),
        ("expected_warmup_invocations", EXPECTED_WARMUP_INVOCATIONS),
    ):
        exact_int(config[field], expected, expected, f"config.{field}")
    exact_string(config["source_git_commit"], expected_commit, "source_git_commit")

    comparisons = config["comparisons"]
    if type(comparisons) is not list or len(comparisons) != len(COMPARISONS):
        fail("config.comparisons roster differs")
    for index, (actual, expected) in enumerate(
        zip(comparisons, EXPECTED_COMPARISON_RECEIPTS)
    ):
        if type(actual) is not dict:
            fail(f"config.comparisons[{index}] is not an object")
        exact_keys(actual, COMPARISON_KEYS, f"config.comparisons[{index}]")
        if actual != expected:
            fail(f"config.comparisons[{index}] differs")

    if set(EXPECTED_CELL_RECEIPTS) != set(CELLS):
        fail("frozen expected cell receipt roster is incomplete")
    cells = config["cells"]
    if type(cells) is not list or len(cells) != len(CELLS):
        fail("config.cells roster differs")
    for index, ((expected_k, expected_b), cell) in enumerate(zip(CELLS, cells)):
        where = f"config.cells[{index}]"
        if type(cell) is not dict:
            fail(f"{where} is not an object")
        exact_keys(cell, CELL_KEYS, where)
        exact_int(cell["K"], expected_k, expected_k, f"{where}.K")
        exact_int(cell["block_bytes"], expected_b, expected_b, f"{where}.B")
        for field in (
            "construction_equivalent", "direct_hit_oracle_verified",
            "full_repair_oracle_verified", "high_id_oracle_verified",
            "solved_state_equivalent", "systematic_oracle_verified",
        ):
            if cell[field] is not True:
                fail(f"{where}.{field} is not true")
        exact_string(cell["solved_state_scope"], SOLVED_STATE_SCOPE, where + ".scope")
        invocations = cell["invocations_by_replicate"]
        expected_invocations = [
            invocations_for(expected_k, replicate)
            for replicate in range(REPLICATES)
        ]
        if (
            type(invocations) is not list
            or len(invocations) != REPLICATES
            or any(type(value) is not int for value in invocations)
            or invocations != expected_invocations
            or sum(invocations) != total_invocations(expected_k)
        ):
            fail(f"{where}.invocations_by_replicate differs")
        configuration = cell["equation_configuration"]
        if type(configuration) is not dict:
            fail(f"{where}.equation_configuration is not an object")
        exact_keys(configuration, CONFIGURATION_KEYS, f"{where}.configuration")
        exact_int(configuration["K"], expected_k, expected_k, "configuration.K")
        exact_int(
            configuration["block_bytes"], expected_b, expected_b,
            "configuration.block_bytes",
        )
        for field in CONFIGURATION_KEYS - {"dense_identity_corner", "precode_seed"}:
            exact_int(configuration[field], 0, MAX_UINT32, f"configuration.{field}")
        exact_int(
            configuration["precode_seed"], 0, MAX_UINT64,
            "configuration.precode_seed",
        )
        if type(configuration["dense_identity_corner"]) is not bool:
            fail("configuration.dense_identity_corner is not bool")
        configuration_hash = lower_hash(
            cell["equation_configuration_sha256"],
            f"{where}.equation_configuration_sha256",
        )
        if configuration_hash != sha256_bytes(canonical_bytes(configuration)):
            fail(f"{where} configuration hash differs from its object")
        expected_receipts = EXPECTED_CELL_RECEIPTS[(expected_k, expected_b)]
        exact_keys(expected_receipts, EXPECTED_RECEIPT_FIELDS, where + " frozen")
        for field in EXPECTED_RECEIPT_FIELDS:
            value = cell[field]
            expected_value = expected_receipts[field]
            if field == "source_seed":
                if (
                    type(value) is not str or not value.startswith("0x")
                    or re.fullmatch(r"0x[0-9a-f]+", value) is None
                    or value != expected_source_seed(expected_k, expected_b)
                ):
                    fail(f"{where}.source_seed is not canonical uint64 hex")
            else:
                lower_hash(value, f"{where}.{field}")
            exact_string(value, expected_value, f"{where}.{field}")
        exact_string(
            cell["systematic_oracle_sha256"], cell["source_sha256"],
            f"{where}.systematic_oracle_sha256",
        )
        if cell["direct_hit_oracle_sha256"] != sha256_bytes(
            canonical_bytes(direct_hit_oracle_preimage(cell))
        ):
            fail(f"{where}.direct_hit_oracle_sha256 preimage differs")
        if cell["solved_state_equivalence_receipt_sha256"] != sha256_bytes(
            canonical_bytes(solved_state_preimage(cell))
        ):
            fail(f"{where}.solved-state receipt preimage differs")

    entries = config["descriptors"]
    if type(entries) is not list or len(entries) != 2:
        fail("config.descriptors roster differs")
    hashes: Dict[Tuple[str, str], str] = {}
    values: Dict[str, Dict[str, Any]] = {}
    for index, (mode, emission) in enumerate((
        ("equation", "equation-eval-v1"),
        ("retained", "retained-source-direct-v1"),
    )):
        entry = entries[index]
        where = f"config.descriptors[{index}]"
        if type(entry) is not dict:
            fail(f"{where} is not an object")
        exact_keys(entry, DESCRIPTOR_ENTRY_KEYS, where)
        exact_string(entry["mode"], mode, f"{where}.mode")
        expected = descriptor(SCOPE, emission)
        if entry["descriptor"] != expected:
            fail(f"{where}.descriptor differs")
        digest = lower_hash(entry["descriptor_sha256"], f"{where}.sha256")
        if digest != descriptor_hash(expected):
            fail(f"{where}.descriptor hash differs")
        hashes[(SCOPE, mode)] = digest
        values[mode] = expected
    normalized = dict(values["retained"])
    normalized["systematic_emission"] = "equation-eval-v1"
    if normalized != values["equation"]:
        fail("descriptors differ outside systematic emission")
    if hashes[(SCOPE, "equation")] == hashes[(SCOPE, "retained")]:
        fail("descriptor hashes collide")
    return hashes


def total_invocations(k: int) -> int:
    return max(MINIMUM_INVOCATIONS, (INVOCATION_BUDGET + k - 1) // k)


def invocations_for(k: int, replicate: int) -> int:
    quotient, remainder = divmod(total_invocations(k), REPLICATES)
    return quotient + (1 if replicate < remainder else 0)


def comparison_modes(comparison: str) -> Tuple[str, str]:
    if comparison == "baseline-aa":
        return "equation", "equation"
    if comparison == "candidate-aa":
        return "retained", "retained"
    if comparison == "candidate-over-baseline":
        return "retained", "equation"
    fail("unknown comparison")
    raise AssertionError


def expected_sides(order: str) -> Tuple[str, ...]:
    if order == "ABBA":
        return ("left", "right", "right", "left",
                "right", "left", "left", "right")
    if order == "BAAB":
        return ("right", "left", "left", "right",
                "left", "right", "right", "left")
    fail("unknown panel order")
    raise AssertionError


def panel_key_sha256(k: int, block_bytes: int, comparison: str, scope: str) -> str:
    return sha256_bytes(canonical_bytes({
        "K": k,
        "block_bytes": block_bytes,
        "campaign": CAMPAIGN,
        "comparison": comparison,
        "schema": PANEL_KEY_SCHEMA,
        "scope": scope,
    }))


def panel_order(
    k: int, block_bytes: int, comparison: str, scope: str, replicate: int
) -> str:
    phase_bit = bytes.fromhex(
        panel_key_sha256(k, block_bytes, comparison, scope)
    )[-1] & 1
    return "ABBA" if (replicate & 1) == phase_bit else "BAAB"


def lane_contrast(elapsed: Sequence[int], order: str) -> float:
    logs = [math.log(value) for value in elapsed]

    def contrast(first: int, block_order: str) -> float:
        if block_order == "ABBA":
            return ((logs[first] - logs[first + 1])
                    + (logs[first + 3] - logs[first + 2])) / 2.0
        return ((logs[first + 1] - logs[first])
                + (logs[first + 2] - logs[first + 3])) / 2.0

    opposite = "BAAB" if order == "ABBA" else "ABBA"
    return 0.5 * (contrast(0, order) + contrast(4, opposite))


def validate_panel(
    panel: Dict[str, Any], cpu: int, replicate: int, k: int, block_bytes: int,
    scope: str, comparison: str, hashes: Mapping[Tuple[str, str], str],
) -> Dict[str, float]:
    where = f"panel r{replicate} K{k} B{block_bytes} {scope} {comparison}"
    exact_keys(panel, PANEL_KEYS, where)
    exact_int(panel["K"], k, k, f"{where}.K")
    exact_int(panel["block_bytes"], block_bytes, block_bytes, f"{where}.B")
    exact_string(panel["schema"], PANEL_SCHEMA, f"{where}.schema")
    exact_string(panel["scope"], scope, f"{where}.scope")
    exact_string(panel["comparison"], comparison, f"{where}.comparison")
    exact_int(panel["replicate"], replicate, replicate, f"{where}.replicate")
    exact_int(panel["target_cpu"], cpu, cpu, f"{where}.target_cpu")
    exact_string(
        panel["panel_key_sha256"],
        panel_key_sha256(k, block_bytes, comparison, scope),
        f"{where}.panel_key_sha256",
    )
    order = panel_order(k, block_bytes, comparison, scope, replicate)
    exact_string(panel["order"], order, f"{where}.order")
    expected_n = invocations_for(k, replicate)
    exact_int(
        panel["invocations_per_slot"], expected_n, expected_n,
        f"{where}.invocations_per_slot",
    )
    primary = next(
        metric for metric, metric_scope_name, _ in METRICS
        if metric_scope_name == scope
    )
    exact_string(panel["primary_metric"], primary, f"{where}.primary_metric")
    left_mode, right_mode = comparison_modes(comparison)
    exact_string(
        panel["left_descriptor_sha256"], hashes[(primary, left_mode)],
        f"{where}.left_descriptor_sha256",
    )
    exact_string(
        panel["right_descriptor_sha256"], hashes[(primary, right_mode)],
        f"{where}.right_descriptor_sha256",
    )
    expected_left_direct = (
        k if scope == SCOPE and left_mode == "retained" else 0
    )
    expected_right_direct = (
        k if scope == SCOPE and right_mode == "retained" else 0
    )
    exact_int(
        panel["left_direct_systematic_packets"], expected_left_direct,
        expected_left_direct, f"{where}.left_direct_systematic_packets",
    )
    exact_int(
        panel["right_direct_systematic_packets"], expected_right_direct,
        expected_right_direct, f"{where}.right_direct_systematic_packets",
    )
    exact_int(panel["left_outcome_code"], 0, 0, f"{where}.left_outcome_code")
    exact_int(panel["right_outcome_code"], 0, 0, f"{where}.right_outcome_code")

    slots = panel["slots"]
    sides = expected_sides(order)
    if type(slots) is not list or len(slots) != 8:
        fail(f"{where}.slots does not contain eight entries")
    upper_half = (expected_n + 1) // 2
    lower_half = expected_n // 2
    expected_slot_counts = (upper_half,) * 4 + (lower_half,) * 4
    primary_elapsed: List[int] = []
    for index, (slot, side, invocation_count) in enumerate(
        zip(slots, sides, expected_slot_counts)
    ):
        if type(slot) is not dict:
            fail(f"{where}.slots[{index}] is not an object")
        exact_keys(slot, SLOT_KEYS, f"{where}.slots[{index}]")
        exact_string(slot["side"], side, f"{where}.slots[{index}].side")
        exact_int(
            slot["invocation_count"], invocation_count, invocation_count,
            f"{where}.slots[{index}].invocation_count",
        )
        primary_elapsed.append(exact_int(
            slot["elapsed_ns"], 1, MAX_INT63, f"{where}.slots[{index}].elapsed_ns"
        ))
    if sum(expected_slot_counts) != 4 * expected_n:
        fail(f"{where} measured invocation split differs")
    return {primary: lane_contrast(primary_elapsed, order)}


def sample_summary(values: Sequence[float]) -> Dict[str, Any]:
    if len(values) != REPLICATES or any(not math.isfinite(value) for value in values):
        fail("statistical group is incomplete or non-finite")
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
    if variance < 0.0 or not math.isfinite(variance):
        fail("statistical variance is invalid")
    standard_error = math.sqrt(variance / len(values))
    lower_log = mean - T11_975 * standard_error
    upper_log = mean + T11_975 * standard_error
    return {
        "geometric_mean": math.exp(mean),
        "lower95": math.exp(lower_log),
        "lower95_log": lower_log,
        "log_mean": mean,
        "log_standard_error": standard_error,
        "n": len(values),
        "upper95": math.exp(upper_log),
        "upper95_log": upper_log,
    }


def band_for(k: int) -> str:
    if k in (32, 100):
        return "small"
    if k in (1000, 2048, 8192):
        return "medium"
    if k in (20000, 32768):
        return "large"
    raise AssertionError("unsealed K")


def aggregate_replicates(
    logs: Mapping[Tuple[str, str, int, int], List[float]], metric: str,
    cells: Sequence[Tuple[int, int]],
) -> List[float]:
    return [
        math.fsum(logs[(metric, "candidate-over-baseline", k, b)][replicate]
                  for k, b in cells) / len(cells)
        for replicate in range(REPLICATES)
    ]


def validate_transcript(
    records: Sequence[Dict[str, Any]], cpu: int, expected_commit: str
) -> Tuple[Dict[str, Any], Dict[str, Any], List[str]]:
    hashes = validate_config(records[0], cpu, expected_commit)
    logs: Dict[Tuple[str, str, int, int], List[float]] = {}
    measured_invocations = 0
    warmup_invocations = 0
    index = 1
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in PANEL_SCOPES:
                for comparison in COMPARISONS:
                    lanes = validate_panel(
                        records[index], cpu, replicate, k, block_bytes, scope,
                        comparison, hashes,
                    )
                    for metric, lane in lanes.items():
                        logs.setdefault((metric, comparison, k, block_bytes), []).append(lane)
                    measured_invocations += sum(
                        slot["invocation_count"] for slot in records[index]["slots"]
                    )
                    warmup_invocations += 2
                    index += 1
    terminal = records[index]
    exact_keys(
        terminal,
        {
            "invocation_count", "measured_invocation_count", "panel_count",
            "record_count", "schema", "status", "warmup_invocation_count",
        },
        "terminal",
    )
    for field, expected in (
        ("invocation_count", EXPECTED_TOTAL_INVOCATIONS),
        ("measured_invocation_count", EXPECTED_TIMED_INVOCATIONS),
        ("warmup_invocation_count", EXPECTED_WARMUP_INVOCATIONS),
        ("panel_count", EXPECTED_PANEL_COUNT),
        ("record_count", EXPECTED_RECORD_COUNT),
    ):
        exact_int(terminal[field], expected, expected, f"terminal.{field}")
    if (
        measured_invocations != EXPECTED_TIMED_INVOCATIONS
        or warmup_invocations != EXPECTED_WARMUP_INVOCATIONS
        or measured_invocations + warmup_invocations
        != EXPECTED_TOTAL_INVOCATIONS
    ):
        fail("independently reconstructed invocation counts differ")
    exact_string(terminal["schema"], TERMINAL_SCHEMA, "terminal.schema")
    exact_string(terminal["status"], "complete", "terminal.status")

    failures: List[str] = []
    group_stats: Dict[str, Any] = {}
    practical_log = math.log1p(0.02)
    for metric, _, _ in METRICS:
        for comparison in COMPARISONS:
            for k, block_bytes in CELLS:
                label = f"{metric}:{comparison}:K{k}:B{block_bytes}"
                summary = sample_summary(logs[(metric, comparison, k, block_bytes)])
                group_stats[label] = summary
                if comparison != "candidate-over-baseline":
                    if not (
                        summary["lower95_log"] > -practical_log
                        and summary["upper95_log"] < practical_log
                    ):
                        failures.append(f"aa_ci:{label}")
                else:
                    if summary["log_mean"] > practical_log:
                        failures.append(f"cell_point:{label}")

    aggregates: Dict[str, Any] = {}
    for metric, _, _ in METRICS:
        limit_log = math.log(0.99)
        groups: List[Tuple[str, Sequence[Tuple[int, int]]]] = [
            ("all", CELLS),
            ("B64", tuple(cell for cell in CELLS if cell[1] == 64)),
            ("B1280", tuple(cell for cell in CELLS if cell[1] == 1280)),
        ]
        groups.extend(
            (band, tuple(cell for cell in CELLS if band_for(cell[0]) == band))
            for band in ("small", "medium", "large")
        )
        metric_stats: Dict[str, Any] = {}
        for name, cells in groups:
            summary = sample_summary(aggregate_replicates(logs, metric, cells))
            metric_stats[name] = summary
            if not summary["upper95_log"] < limit_log:
                failures.append(f"aggregate_upper95:{metric}:{name}")
        aggregates[metric] = metric_stats

    statistics = {
        "aggregates": aggregates,
        "failed_gates": failures,
        "groups": group_stats,
        "student_t_975_df11": T11_975,
    }
    return records[0], statistics, failures


def read_sealed_source(
    root: Path, relative: str,
    expected: Optional[Mapping[str, Any]] = None,
    deadline: Optional[float] = None,
) -> Tuple[bytes, os.stat_result]:
    if deadline is not None and time.monotonic() >= deadline:
        fail("sealed source read reached its global deadline")
    path = root / relative
    fd = os.open(str(path), nonblocking_read_flags("source sealing"))
    try:
        info = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_nlink != 1
            or info.st_size < 0
            or info.st_size > MAX_SOURCE_FILE_BYTES
            or stat.S_IMODE(info.st_mode) != 0o444
            or info.st_uid != EXPECTED_CAMPAIGN_UID
            or info.st_gid != EXPECTED_CAMPAIGN_GID
            or not same_file_receipt(info, named)
        ):
            fail(f"sealed source file policy differs: {relative}")
        data = bytearray()
        offset = 0
        while offset < info.st_size:
            if deadline is not None and time.monotonic() >= deadline:
                fail("sealed source read reached its global deadline")
            block = os.pread(fd, min(1024 * 1024, info.st_size - offset), offset)
            if not block:
                fail(f"sealed source file is short: {relative}")
            data.extend(block)
            offset += len(block)
        info_after = os.fstat(fd)
        named_after = os.stat(str(path), follow_symlinks=False)
        if (
            not same_file_receipt(info, info_after)
            or not same_file_receipt(info_after, named_after)
        ):
            fail(f"sealed source changed during held-FD read: {relative}")
        if expected is not None and (
            len(data) != expected["bytes"]
            or sha256_bytes(bytes(data)) != expected["sha256"]
            or git_blob_oid(bytes(data)) != expected["git_blob_oid"]
            or not canonical_equal(stat_receipt(info_after), expected["stat"])
        ):
            fail(f"sealed source receipt changed: {relative}")
        return bytes(data), info_after
    finally:
        os.close(fd)


def source_manifest(
    root: Path, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    entries: List[Dict[str, Any]] = []
    preimage = bytearray()
    for relative in SOURCE_PATHS:
        data, info = read_sealed_source(root, relative, deadline=deadline)
        digest = sha256_bytes(data)
        entries.append({
            "bytes": len(data), "git_blob_oid": git_blob_oid(data),
            "path": relative, "sha256": digest, "stat": stat_receipt(info),
        })
        preimage.extend(f"{digest}  {relative}\n".encode("ascii"))
    return {"entries": entries, "sha256": sha256_bytes(bytes(preimage))}


def validate_sampler_source_authority(
    args: argparse.Namespace, manifest: Mapping[str, Any], source_root: Path,
) -> Mapping[str, Any]:
    validate_source_manifest_receipt(manifest, "sampler source manifest")
    expected_path = source_root / SAMPLER_SOURCE_PATH
    matches = [
        entry for entry in manifest["entries"]
        if entry["path"] == SAMPLER_SOURCE_PATH
    ]
    if len(matches) != 1:
        fail("sampler source manifest entry differs")
    entry = matches[0]
    if (
        args.sampler_script != expected_path
        or args.expected_sampler_script_sha256 != entry["sha256"]
        or entry["stat"]["mode"] != 0o444
        or entry["stat"]["uid"] != args.expected_sampler_uid
        or entry["stat"]["gid"] != args.expected_sampler_gid
    ):
        fail("sampler run authority differs from its sealed source entry")
    return entry


@dataclass(frozen=True)
class SealedHealthModules:
    contract: Any
    native: Any
    runner: Any
    receipt: Dict[str, Any]


def load_sealed_health_modules(
    root: Path, manifest: Mapping[str, Any],
) -> SealedHealthModules:
    """Compile the manifest-bound helper source bytes without consulting pyc."""
    if (
        not sys.flags.isolated
        or not sys.dont_write_bytecode
        or sys.flags.optimize != 0
    ):
        fail("sealed helper loading requires unoptimized Python -I -B")
    validate_source_manifest_receipt(manifest, "health loader source manifest")
    by_path = {entry["path"]: entry for entry in manifest["entries"]}
    modules: Dict[str, Any] = {}
    installed: List[str] = []
    original_path = list(sys.path)
    try:
        for module_name, relative in HEALTH_MODULE_SOURCES:
            if module_name in sys.modules:
                fail(f"health helper module was imported before sealing: {module_name}")
            module = types.ModuleType(module_name)
            module.__file__ = str(root / relative)
            module.__package__ = ""
            modules[module_name] = module
            sys.modules[module_name] = module
            installed.append(module_name)
        module_receipts: List[Dict[str, Any]] = []
        for module_name, relative in HEALTH_MODULE_SOURCES:
            expected = by_path[relative]
            data, _ = read_sealed_source(root, relative, expected)
            digest = sha256_bytes(data)
            try:
                code = compile(
                    data, str(root / relative), "exec",
                    dont_inherit=True, optimize=0,
                )
                exec(code, modules[module_name].__dict__)
            finally:
                # The helper sources insert their directory for standalone use.
                # Restore the isolated interpreter path before executing the
                # next helper so untracked repo files cannot shadow stdlib.
                sys.path[:] = original_path
            module_receipts.append({
                "bytes": len(data), "module": module_name,
                "path": relative, "sha256": digest,
            })
        receipt: Dict[str, Any] = {
            "dont_write_bytecode": True,
            "isolated": True,
            "modules": module_receipts,
            "optimize": 0,
            "receipt_sha256": None,
            "schema": HEALTH_LOADER_SCHEMA,
        }
        receipt["receipt_sha256"] = sha256_bytes(canonical_bytes(receipt))
        return SealedHealthModules(
            contract=modules["wh2_benchmark_contract"],
            native=modules["wh2_native_short_screen"],
            runner=modules["wh2_run_native_short_screen"],
            receipt=receipt,
        )
    except BaseException:
        for name in installed:
            if sys.modules.get(name) is modules.get(name):
                del sys.modules[name]
        raise
    finally:
        sys.path[:] = original_path


def file_sha256_fd(
    fd: int, size: int, deadline: Optional[float] = None,
) -> str:
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        if deadline is not None and time.monotonic() >= deadline:
            fail("held-file hashing reached its global deadline")
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("short read while hashing binary")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def stat_receipt(st: os.stat_result) -> Dict[str, Any]:
    return {
        "device": st.st_dev, "gid": st.st_gid, "inode": st.st_ino,
        "mode": stat.S_IMODE(st.st_mode), "mtime_ns": st.st_mtime_ns,
        "nlink": st.st_nlink, "size": st.st_size, "uid": st.st_uid,
    }


def same_file_receipt(left: os.stat_result, right: os.stat_result) -> bool:
    return stat_receipt(left) == stat_receipt(right)


def exception_text(exc: BaseException) -> str:
    text = f"{type(exc).__name__}: {exc}"
    return text if len(text) <= 4000 else text[:3997] + "..."


def terminate_process_group(process: subprocess.Popen) -> Optional[str]:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return None
    except OSError as exc:
        return exception_text(exc)
    return None


class SignalGuard:
    """Latch termination signals without allowing them to strand artifacts."""

    SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)

    def __init__(self) -> None:
        self.first_signal: Optional[str] = None
        self.observed_signals: List[str] = []
        self.kill_error: Optional[str] = None
        self.process: Optional[subprocess.Popen] = None
        self.old_handlers: Dict[int, Any] = {}
        self.old_mask: Optional[set] = None
        self.output_bundle: Optional[Any] = None
        self.entered = False
        self.seal_blocked = False
        self.signal_decision_cutoff = False
        self.evidence_committed = False

    def _handler(self, signum: int, _frame: Any) -> None:
        if self.signal_decision_cutoff:
            return
        name = signal.Signals(signum).name
        if self.first_signal is None:
            self.first_signal = name
        if name not in self.observed_signals:
            self.observed_signals.append(name)
        if self.process is not None:
            error = terminate_process_group(self.process)
            if error is not None and self.kill_error is None:
                self.kill_error = error

    def __enter__(self) -> "SignalGuard":
        if self.entered:
            fail("signal guard cannot be reused")
        if not all(
            hasattr(signal, name)
            for name in ("pthread_sigmask", "sigpending", "sigtimedwait")
        ):
            fail(
                "signal-safe publication requires pthread_sigmask, "
                "sigpending, and sigtimedwait"
            )
        # Block first so no guarded signal can land between observing and
        # replacing the old dispositions.  Reservation starts only after this
        # method has installed the handlers and restored the incoming mask.
        self.old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, self.SIGNALS)
        if any(signum in self.old_mask for signum in self.SIGNALS):
            signal.pthread_sigmask(signal.SIG_SETMASK, self.old_mask)
            fail("guarded termination signals are already blocked")
        try:
            for signum in self.SIGNALS:
                self.old_handlers[signum] = signal.getsignal(signum)
                signal.signal(signum, self._handler)
        except BaseException:
            for signum, handler in self.old_handlers.items():
                signal.signal(signum, handler)
            signal.pthread_sigmask(signal.SIG_SETMASK, self.old_mask)
            raise
        signal.pthread_sigmask(signal.SIG_SETMASK, self.old_mask)
        self.entered = True
        return self

    def attach(self, process: subprocess.Popen) -> None:
        self.process = process
        if self.first_signal is not None:
            error = terminate_process_group(process)
            if error is not None and self.kill_error is None:
                self.kill_error = error

    def detach(self, process: subprocess.Popen) -> None:
        if self.process is process:
            self.process = None

    def own_output_bundle(self, bundle: Any) -> None:
        if self.output_bundle is not None:
            fail("signal guard already owns an output bundle")
        self.output_bundle = bundle

    def block_for_commit(self) -> None:
        if not self.entered or self.signal_decision_cutoff:
            fail("signal guard commit transition is invalid")
        if not self.seal_blocked:
            signal.pthread_sigmask(signal.SIG_BLOCK, self.SIGNALS)
            self.seal_blocked = True

    def _collect_pending_unchecked(self) -> bool:
        changed = False
        pending = signal.sigpending()
        for signum in self.SIGNALS:
            if signum not in pending:
                continue
            name = signal.Signals(signum).name
            if self.first_signal is None:
                self.first_signal = name
            if name not in self.observed_signals:
                self.observed_signals.append(name)
                changed = True
            # Consume the blocked instance after recording it.  Otherwise it
            # would be delivered under the restored pre-run disposition on
            # context exit, defeating the guarded invalid-evidence return.
            consumed = signal.sigtimedwait({signum}, 0.0)
            if consumed is None and signum in signal.sigpending():
                fail(f"could not consume pending {name}")
        return changed

    def collect_pending(self) -> bool:
        if not self.seal_blocked or self.signal_decision_cutoff:
            fail("pending signals require an uncommitted blocked guard")
        return self._collect_pending_unchecked()

    def commit_logical_seal(self) -> bool:
        if not self.seal_blocked or self.signal_decision_cutoff:
            fail("signal guard logical commit transition is invalid")
        # The kernel set is still blocked, so freezing this in-memory cutoff
        # before the single conservative drain cannot lose a pre-cutoff
        # signal: any such instance remains pending until the drain consumes
        # and records it.  Instances racing just after the cutoff may also be
        # included, which is deliberately fail-closed.
        self.signal_decision_cutoff = True
        return self._collect_pending_unchecked()

    def __exit__(self, _kind: Any, _value: Any, _traceback: Any) -> None:
        if not self.entered:
            return
        # A durable COMPLETE is authoritative.  Keep the guarded set blocked
        # until immediate process exit so a signal arriving after that sole
        # commit cannot contradict the already-frozen artifact/exit class.
        if self.evidence_committed:
            self.entered = False
            return
        bundle = self.output_bundle
        if bundle is not None and not bundle.closed:
            bundle.signal_name = self.first_signal
            bundle.signal_names = list(self.observed_signals)
            try:
                emergency = bundle.emergency_summary(
                    "signal-guard outer ownership drain", []
                )
            except BaseException:
                emergency = b""
            try:
                bundle.finish(emergency, self)
            except BaseException:
                # OutputBundle.finish owns its own idempotent descriptor drain.
                pass
        if self.signal_decision_cutoff:
            # A post-cutoff failure has already poisoned/closed its output in
            # the final owner.  Keep signals blocked through immediate invalid
            # process exit, but do not mislabel evidence as committed.
            self.entered = False
            return
        if not self.seal_blocked:
            signal.pthread_sigmask(signal.SIG_BLOCK, self.SIGNALS)
        try:
            for signum, handler in self.old_handlers.items():
                signal.signal(signum, handler)
        finally:
            if self.old_mask is None:
                raise RuntimeError("signal guard lost its incoming signal mask")
            signal.pthread_sigmask(signal.SIG_SETMASK, self.old_mask)
            self.entered = False


@dataclass
class CaptureResult:
    stdout: bytes = b""
    stderr: bytes = b""
    timed_out: bool = False
    output_overflow: bool = False
    sampler_event: str = "none"
    error: str = "none"


def bounded_capture(
    process: subprocess.Popen, deadline: float,
    signal_guard: Optional[SignalGuard] = None,
    sampler_monitor_handles: Optional[Mapping[str, Any]] = None,
    *, stdout_limit: int = MAX_STDOUT_BYTES,
    stderr_limit: int = MAX_STDERR_BYTES,
    kill_grace_seconds: float = 5.0,
) -> CaptureResult:
    result = CaptureResult()
    exact_int(stdout_limit, 1, MAX_STDOUT_BYTES, "stdout capture limit")
    exact_int(stderr_limit, 1, MAX_STDERR_BYTES, "stderr capture limit")
    if (
        type(kill_grace_seconds) is not float
        or not math.isfinite(kill_grace_seconds)
        or not 0.0 < kill_grace_seconds <= 5.0
    ):
        fail("capture kill grace differs")
    if process.stdout is None or process.stderr is None:
        result.error = "ValidationError: screen process pipes are absent"
        error = terminate_process_group(process)
        if error is not None:
            result.error += f"; process-group kill: {error}"
        try:
            process.wait(timeout=5.0)
        except BaseException as exc:
            result.error += f"; reap: {exception_text(exc)}"
        return result
    selector = selectors.DefaultSelector()
    streams = ((process.stdout, "stdout"), (process.stderr, "stderr"))
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    limits = {"stdout": stdout_limit, "stderr": stderr_limit}
    killed_at: Optional[float] = None
    completed = False
    errors: List[str] = []
    if signal_guard is not None:
        signal_guard.attach(process)
    try:
        for stream, name in streams:
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, name)
        # A child can close both capture pipes and remain alive.  Keep the
        # deadline, signal, and sampler interlocks active until both the
        # streams are drained and the child has actually exited.
        while selector.get_map() or process.poll() is None:
            now = time.monotonic()
            if (
                sampler_monitor_handles is not None
                and result.sampler_event == "none"
            ):
                sampler_event = poll_sampler_supervision(
                    sampler_monitor_handles
                )
                if sampler_event != "none":
                    result.sampler_event = sampler_event
                    if killed_at is None:
                        error = terminate_process_group(process)
                        if error is not None:
                            errors.append(
                                "sampler-event process-group kill: " + error
                            )
                        killed_at = now
            if (
                signal_guard is not None
                and signal_guard.first_signal is not None
                and killed_at is None
            ):
                error = terminate_process_group(process)
                if error is not None:
                    errors.append(f"signal process-group kill: {error}")
                killed_at = now
            if killed_at is None and now >= deadline:
                result.timed_out = True
                error = terminate_process_group(process)
                if error is not None:
                    errors.append(f"timeout process-group kill: {error}")
                killed_at = now
            if killed_at is not None and now - killed_at >= kill_grace_seconds:
                errors.append("screen process pipes did not close after SIGKILL")
                break
            wait_until = (
                killed_at + kill_grace_seconds
                if killed_at is not None else deadline
            )
            wait_seconds = max(0.0, min(0.05, wait_until - now))
            if selector.get_map():
                events = selector.select(wait_seconds)
            else:
                time.sleep(wait_seconds)
                events = ()
            for key, _ in events:
                stream = key.fileobj
                name = key.data
                try:
                    data = os.read(stream.fileno(), 65536)
                except BlockingIOError:
                    continue
                if not data:
                    selector.unregister(stream)
                    stream.close()
                    continue
                # Retain exactly a bounded prefix.  Continue draining after a
                # kill, but never append beyond limit+1.
                available = max(0, limits[name] + 1 - len(buffers[name]))
                if available:
                    buffers[name].extend(data[:available])
                if len(data) > available or len(buffers[name]) > limits[name]:
                    result.output_overflow = True
                    if killed_at is None:
                        error = terminate_process_group(process)
                        if error is not None:
                            errors.append(f"overflow process-group kill: {error}")
                        killed_at = time.monotonic()
        remaining = (
            max(0.0, killed_at + kill_grace_seconds - time.monotonic())
            if killed_at is not None else max(0.0, deadline - time.monotonic())
        )
        process.wait(timeout=remaining)
        completed = True
    except BaseException as exc:
        errors.append(f"capture: {exception_text(exc)}")
    finally:
        if not completed:
            error = terminate_process_group(process)
            if error is not None:
                errors.append(f"cleanup process-group kill: {error}")
            try:
                process.wait(timeout=kill_grace_seconds)
                completed = True
            except BaseException as exc:
                errors.append(f"cleanup reap: {exception_text(exc)}")
        if signal_guard is not None:
            signal_guard.detach(process)
        try:
            selector.close()
        except BaseException as exc:
            errors.append(f"selector close: {exception_text(exc)}")
        for stream, _ in streams:
            if not stream.closed:
                try:
                    stream.close()
                except BaseException as exc:
                    errors.append(f"stream close: {exception_text(exc)}")
    result.stdout = bytes(buffers["stdout"])
    result.stderr = bytes(buffers["stderr"])
    if errors:
        result.error = "; ".join(errors)
    return result


OUTPUT_NAMES = ("raw.jsonl", "stderr.txt", "summary.json")


@dataclass(frozen=True)
class BundleFinishResult:
    """Separate the immutable decision from post-commit durability status.

    Pre-commit failures are reflected in a full-schema invalid summary before
    ``logical_commit_succeeded`` can become true.  Post-commit failures happen
    after that immutable summary decision: they make publication fail closed,
    but can neither rewrite nor reclassify the committed summary.
    """

    logical_commit_succeeded: bool = False
    precommit_failures: Tuple[str, ...] = tuple()
    postcommit_failures: Tuple[str, ...] = tuple()

SUMMARY_KEYS = {
    "binary", "config", "elapsed_seconds", "expected_source_git_commit",
    "expected_source_manifest_sha256", "failure", "git_after", "git_before",
    "health", "health_module_loader", "outcome",
    "output_bundle", "process_exit_code", "publication_failures",
    "raw_bytes", "raw_complete", "raw_record_count", "raw_sha256",
    "schema", "signal", "signal_names", "source_manifest_after",
    "source_manifest_before",
    "statistics", "stderr_bytes", "stderr_sha256",
    "summary_preimage_sha256", "target_cpu",
}
BINARY_RECEIPT_KEYS = {
    "capture_error", "child_affinity_after_spawn", "child_pid",
    "child_process_start_ticks", "child_reap_monotonic_ns",
    "child_start_monotonic_ns",
    "execution_finished_monotonic_ns", "execution_started_monotonic_ns",
    "expected_sha256", "output_overflow", "path", "path_stat_after",
    "process_started", "sha256_after", "sha256_before", "stat_after",
    "stat_before", "timed_out",
}
STAT_RECEIPT_KEYS = {
    "device", "gid", "inode", "mode", "mtime_ns", "nlink", "size", "uid",
}
OUTPUT_BUNDLE_KEYS = {"directory", "files", "parent", "path"}
OUTPUT_DIRECTORY_KEYS = {
    "device", "inode", "nlink", "reserved_mode", "sealed_mode",
}
OUTPUT_FILE_KEYS = OUTPUT_DIRECTORY_KEYS
OUTPUT_PARENT_KEYS = {"device", "initial_nlink", "inode", "mode"}
HEALTH_KEYS = {
    "admission_sibling_ticks", "child_reap_monotonic_ns",
    "child_start_monotonic_ns", "collection_failures",
    "controller_core", "controller_cpu", "controller_initial_affinity",
    "controller_singleton_affinity_end", "edac_policy", "evidence_status",
    "receipt_sha256", "violations",
    "sampler", "sampler_admission", "sampler_admission_receipt_sha256",
    "sampler_core", "sampler_cpu", "sampler_receipt_sha256",
    "sampler_terminal", "sampler_terminal_receipt_sha256",
    "schema", "sibling_non_idle_tick_cap", "sibling_tick_policy",
    "sibling_ticks", "target_core", "target_cpu", "target_threads",
    "thermal", "thermal_max_millic", "terminal_status",
}
SAMPLER_RECEIPT_KEYS = {
    "cmdline_argv", "cmdline_sha256", "cpu", "csv_device", "csv_inode", "csv_path",
    "csv_bytes", "csv_sha256", "csv_stat",
    "evidence_parent",
    "environ_sha256", "environment", "environment_sha256",
    "executable_device", "executable_inode", "executable_path",
    "executable_sha256", "executable_stat", "pid", "pid_file",
    "process_affinity", "process_gid", "process_security",
    "process_start_ticks", "process_uid",
    "receipt_file",
    "schema", "script_device", "script_inode", "script_path",
    "script_sha256", "script_stat", "terminal_status", "window_end_monotonic_ns",
    "window_start_monotonic_ns", "validation_header_ascii", "validation_jsonl",
}
SAMPLER_PARENT_KEYS = {"path", "stat"}
HEALTH_LOADER_KEYS = {
    "dont_write_bytecode", "isolated", "modules", "optimize",
    "receipt_sha256", "schema",
}
HEALTH_LOADER_MODULE_KEYS = {"bytes", "module", "path", "sha256"}
HEALTH_MODULE_SOURCES = (
    ("wh2_benchmark_contract", "bench/wh2_benchmark_contract.py"),
    ("wh2_native_short_screen", "bench/wh2_native_short_screen.py"),
    ("wh2_run_native_short_screen", "bench/wh2_run_native_short_screen.py"),
)
SOURCE_MANIFEST_KEYS = {"entries", "sha256"}
SOURCE_MANIFEST_ENTRY_KEYS = {
    "bytes", "git_blob_oid", "path", "sha256", "stat",
}
GIT_RECEIPT_KEYS = {
    "executable", "head", "source_blob_oids", "source_blob_roster_sha256",
    "source_roster_sha256", "tracked_index_flags_sha256",
    "tracked_status_sha256", "worktree_root",
}
GIT_EXECUTABLE_KEYS = {"path", "sha256", "stat"}
GIT_SOURCE_BLOB_KEYS = {"head_oid", "path", "worktree_oid"}
THERMAL_RECEIPT_KEYS = {
    "cpu", "cpu_tctl_max_millic", "csv_device", "csv_inode", "csv_path",
    "dimm_max_millic", "dimm_read_errors", "edac_ce_max", "edac_ue_max",
    "invalid_sample_count", "parse_failures", "pid", "process_start_ticks",
    "sample_count", "script_path", "script_sha256", "terminal_status",
    "valid_sample_count", "window_csv_ascii", "window_csv_bytes",
    "window_csv_sha256",
    "validation_attempt_errors_total", "validation_device",
    "validation_failures", "validation_inode", "validation_jsonl_ascii",
    "validation_jsonl_bytes", "validation_jsonl_sha256", "validation_path",
    "validation_sample_count",
    "window_end_monotonic_ns", "window_start_monotonic_ns",
}
STATISTICS_KEYS = {
    "aggregates", "failed_gates", "groups", "student_t_975_df11",
}
SAMPLE_SUMMARY_KEYS = {
    "geometric_mean", "lower95", "lower95_log", "log_mean",
    "log_standard_error", "n", "upper95", "upper95_log",
}
AGGREGATE_GROUPS = (
    "all", "B64", "B1280", "small", "medium", "large",
)


def make_summary_preimage(summary: Dict[str, Any]) -> Dict[str, Any]:
    exact_keys(summary, SUMMARY_KEYS, "summary")
    preimage = dict(summary)
    preimage["summary_preimage_sha256"] = None
    summary["summary_preimage_sha256"] = sha256_bytes(canonical_bytes(preimage))
    return summary


def canonical_equal(left: Any, right: Any) -> bool:
    """Type-exact equality for JSON values (notably, bool never equals int)."""
    return canonical_bytes(left) == canonical_bytes(right)


def exact_absolute_path(value: Any, where: str) -> str:
    if (
        type(value) is not str
        or not value
        or not Path(value).is_absolute()
        or os.path.normpath(value) != value
    ):
        fail(f"{where} is not a canonical absolute path")
    return value


def nonblocking_read_flags(
    where: str, *, nofollow: bool = True, noatime: bool = False,
) -> int:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    nonblock = getattr(os, "O_NONBLOCK", 0)
    if not nonblock:
        fail(f"{where} requires Linux O_NONBLOCK")
    flags |= nonblock
    if nofollow:
        nofollow_flag = getattr(os, "O_NOFOLLOW", 0)
        if not nofollow_flag:
            fail(f"{where} requires Linux O_NOFOLLOW")
        flags |= nofollow_flag
    if noatime:
        noatime_flag = getattr(os, "O_NOATIME", 0)
        if not noatime_flag:
            fail(f"{where} requires Linux O_NOATIME")
        flags |= noatime_flag
    return flags


def validate_stat_receipt(
    receipt: Any, where: str, *, minimum_nlink: int = 1,
) -> None:
    if type(receipt) is not dict:
        fail(f"{where} is not an object")
    exact_keys(receipt, STAT_RECEIPT_KEYS, where)
    exact_int(receipt["device"], 0, MAX_UINT64, f"{where}.device")
    exact_int(receipt["gid"], 0, MAX_UINT32, f"{where}.gid")
    exact_int(receipt["inode"], 1, MAX_UINT64, f"{where}.inode")
    exact_int(receipt["mode"], 0, 0o7777, f"{where}.mode")
    exact_int(
        receipt["mtime_ns"], -(1 << 63), MAX_INT63, f"{where}.mtime_ns"
    )
    exact_int(
        receipt["nlink"], minimum_nlink, MAX_UINT64,
        f"{where}.nlink",
    )
    exact_int(receipt["size"], 0, MAX_INT63, f"{where}.size")
    exact_int(receipt["uid"], 0, MAX_UINT32, f"{where}.uid")


def validate_output_bundle_receipt(output: Any) -> None:
    if type(output) is not dict:
        fail("summary output bundle is not an object")
    exact_keys(output, OUTPUT_BUNDLE_KEYS, "summary output bundle")
    exact_absolute_path(output["path"], "summary output path")
    directory = output["directory"]
    parent = output["parent"]
    files = output["files"]
    if type(directory) is not dict or type(parent) is not dict or type(files) is not dict:
        fail("summary output nested receipt is not an object")
    exact_keys(directory, OUTPUT_DIRECTORY_KEYS, "output directory")
    exact_keys(parent, OUTPUT_PARENT_KEYS, "output parent")
    exact_keys(files, OUTPUT_NAMES, "output files")
    exact_int(directory["device"], 0, MAX_UINT64, "output directory device")
    exact_int(directory["inode"], 1, MAX_UINT64, "output directory inode")
    exact_int(directory["nlink"], 1, MAX_UINT64, "output directory nlink")
    exact_int(directory["reserved_mode"], 0o700, 0o700, "output reserved mode")
    exact_int(directory["sealed_mode"], 0o500, 0o500, "output sealed mode")
    exact_int(parent["device"], 0, MAX_UINT64, "output parent device")
    exact_int(parent["inode"], 1, MAX_UINT64, "output parent inode")
    exact_int(parent["initial_nlink"], 1, MAX_UINT64, "output parent nlink")
    exact_int(parent["mode"], 0, 0o7777, "output parent mode")
    for name in OUTPUT_NAMES:
        receipt = files[name]
        if type(receipt) is not dict:
            fail(f"output file receipt {name} is not an object")
        exact_keys(receipt, OUTPUT_FILE_KEYS, f"output file receipt {name}")
        exact_int(receipt["device"], 0, MAX_UINT64, f"output {name} device")
        exact_int(receipt["inode"], 1, MAX_UINT64, f"output {name} inode")
        exact_int(receipt["nlink"], 1, 1, f"output {name} nlink")
        exact_int(receipt["reserved_mode"], 0o600, 0o600, f"output {name} reserved mode")
        exact_int(receipt["sealed_mode"], 0o400, 0o400, f"output {name} sealed mode")


def validate_source_manifest_receipt(manifest: Any, where: str) -> None:
    if type(manifest) is not dict:
        fail(f"{where} is not an object")
    exact_keys(manifest, SOURCE_MANIFEST_KEYS, where)
    entries = manifest["entries"]
    if type(entries) is not list or len(entries) != len(SOURCE_PATHS):
        fail(f"{where} source roster differs")
    preimage = bytearray()
    for index, (entry, expected_path) in enumerate(zip(entries, SOURCE_PATHS)):
        entry_where = f"{where}.entries[{index}]"
        if type(entry) is not dict:
            fail(f"{entry_where} is not an object")
        exact_keys(entry, SOURCE_MANIFEST_ENTRY_KEYS, entry_where)
        exact_string(entry["path"], expected_path, f"{entry_where}.path")
        exact_int(
            entry["bytes"], 0, MAX_SOURCE_FILE_BYTES,
            f"{entry_where}.bytes",
        )
        digest = lower_hash(entry["sha256"], f"{entry_where}.sha256")
        if (
            type(entry["git_blob_oid"]) is not str
            or LOWER40.fullmatch(entry["git_blob_oid"]) is None
        ):
            fail(f"{entry_where}.git_blob_oid differs")
        validate_stat_receipt(entry["stat"], f"{entry_where}.stat")
        if (
            entry["stat"]["size"] != entry["bytes"]
            or entry["stat"]["size"] > MAX_SOURCE_FILE_BYTES
            or entry["stat"]["nlink"] != 1
            or entry["stat"]["mode"] != 0o444
            or entry["stat"]["uid"] != EXPECTED_CAMPAIGN_UID
            or entry["stat"]["gid"] != EXPECTED_CAMPAIGN_GID
        ):
            fail(f"{entry_where} source stat policy differs")
        preimage.extend(f"{digest}  {expected_path}\n".encode("ascii"))
    digest = lower_hash(manifest["sha256"], f"{where}.sha256")
    if sha256_bytes(bytes(preimage)) != digest:
        fail(f"{where} preimage SHA-256 differs")


def validate_health_loader_receipt(receipt: Any) -> None:
    if type(receipt) is not dict:
        fail("health module loader receipt is not an object")
    exact_keys(receipt, HEALTH_LOADER_KEYS, "health module loader")
    exact_string(receipt["schema"], HEALTH_LOADER_SCHEMA, "health loader schema")
    if receipt["isolated"] is not True or receipt["dont_write_bytecode"] is not True:
        fail("health loader interpreter policy differs")
    exact_int(receipt["optimize"], 0, 0, "health loader optimization level")
    modules = receipt["modules"]
    if type(modules) is not list or len(modules) != len(HEALTH_MODULE_SOURCES):
        fail("health loader module roster differs")
    for index, (module_receipt, expected) in enumerate(
        zip(modules, HEALTH_MODULE_SOURCES)
    ):
        where = f"health loader modules[{index}]"
        if type(module_receipt) is not dict:
            fail(f"{where} is not an object")
        exact_keys(module_receipt, HEALTH_LOADER_MODULE_KEYS, where)
        exact_string(module_receipt["module"], expected[0], f"{where}.module")
        exact_string(module_receipt["path"], expected[1], f"{where}.path")
        exact_int(
            module_receipt["bytes"], 0, MAX_SOURCE_FILE_BYTES,
            f"{where}.bytes",
        )
        lower_hash(module_receipt["sha256"], f"{where}.sha256")
    digest = lower_hash(receipt["receipt_sha256"], "health loader receipt")
    preimage = dict(receipt)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail("health loader receipt SHA-256 differs")


def validate_git_receipt(receipt: Any, expected_commit: str, where: str) -> None:
    if type(receipt) is not dict:
        fail(f"{where} is not an object")
    exact_keys(receipt, GIT_RECEIPT_KEYS, where)
    executable = receipt["executable"]
    if type(executable) is not dict:
        fail(f"{where}.executable is not an object")
    exact_keys(executable, GIT_EXECUTABLE_KEYS, where + ".executable")
    exact_string(
        executable["path"], str(GIT_EXECUTABLE), where + ".executable.path"
    )
    lower_hash(executable["sha256"], where + ".executable.sha256")
    validate_stat_receipt(executable["stat"], where + ".executable.stat")
    if (
        executable["stat"]["mode"] != 0o755
        or executable["stat"]["uid"] != 0
        or executable["stat"]["gid"] != 0
        or executable["stat"]["nlink"] != 1
        or not 1 <= executable["stat"]["size"] <= MAX_BINARY_BYTES
    ):
        fail(f"{where} Git executable policy differs")
    exact_string(receipt["head"], expected_commit, f"{where}.head")
    exact_string(
        receipt["worktree_root"],
        str(Path(__file__).resolve(strict=True).parents[1]),
        f"{where}.worktree_root",
    )
    roster_hash = lower_hash(
        receipt["source_roster_sha256"], f"{where}.source roster"
    )
    status_hash = lower_hash(
        receipt["tracked_status_sha256"], f"{where}.status"
    )
    lower_hash(
        receipt["tracked_index_flags_sha256"],
        f"{where}.tracked index flags",
    )
    blobs = receipt["source_blob_oids"]
    if type(blobs) is not list or len(blobs) != len(SOURCE_PATHS):
        fail(f"{where} source blob roster differs")
    for index, (blob, expected_path) in enumerate(zip(blobs, SOURCE_PATHS)):
        blob_where = f"{where}.source_blob_oids[{index}]"
        if type(blob) is not dict:
            fail(f"{blob_where} is not an object")
        exact_keys(blob, GIT_SOURCE_BLOB_KEYS, blob_where)
        exact_string(blob["path"], expected_path, blob_where + ".path")
        head_oid = blob["head_oid"]
        worktree_oid = blob["worktree_oid"]
        if (
            type(head_oid) is not str or LOWER40.fullmatch(head_oid) is None
            or type(worktree_oid) is not str
            or LOWER40.fullmatch(worktree_oid) is None
            or head_oid != worktree_oid
        ):
            fail(f"{blob_where} raw worktree/HEAD blob binding differs")
    blob_roster_hash = lower_hash(
        receipt["source_blob_roster_sha256"],
        f"{where}.source blob roster",
    )
    expected_roster = b"".join(
        (relative + "\n").encode("ascii") for relative in SOURCE_PATHS
    )
    if (
        roster_hash != sha256_bytes(expected_roster)
        or status_hash != sha256_bytes(b"")
        or blob_roster_hash != sha256_bytes(canonical_bytes(blobs))
    ):
        fail(f"{where} clean tracked-source binding differs")


def validate_sample_summary(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, SAMPLE_SUMMARY_KEYS, where)
    exact_int(value["n"], REPLICATES, REPLICATES, f"{where}.n")
    for name in SAMPLE_SUMMARY_KEYS - {"n"}:
        number = value[name]
        if type(number) is not float or not math.isfinite(number):
            fail(f"{where}.{name} is not finite")
    if (
        value["geometric_mean"] <= 0.0
        or value["lower95"] <= 0.0
        or value["upper95"] <= 0.0
        or value["log_standard_error"] < 0.0
        or value["lower95"] > value["geometric_mean"]
        or value["geometric_mean"] > value["upper95"]
        or value["lower95_log"] > value["log_mean"]
        or value["log_mean"] > value["upper95_log"]
    ):
        fail(f"{where} statistical interval is invalid")


def validate_statistics_receipt(value: Any) -> None:
    if type(value) is not dict:
        fail("summary statistics is not an object")
    exact_keys(value, STATISTICS_KEYS, "summary statistics")
    student_t = value["student_t_975_df11"]
    if type(student_t) is not float or student_t != T11_975:
        fail("summary Student-t constant differs")
    failures = value["failed_gates"]
    if (
        type(failures) is not list
        or len(failures) > MAX_FAILED_GATES
        or any(type(item) is not str or not item or len(item) > 1024 for item in failures)
        or len(set(failures)) != len(failures)
    ):
        fail("summary failed-gate roster differs")
    groups = value["groups"]
    aggregates = value["aggregates"]
    if type(groups) is not dict or type(aggregates) is not dict:
        fail("summary statistical groups are not objects")
    expected_group_names = {
        f"{metric}:{comparison}:K{k}:B{block_bytes}"
        for metric, _, _ in METRICS
        for comparison in COMPARISONS
        for k, block_bytes in CELLS
    }
    exact_keys(groups, expected_group_names, "summary cell statistics")
    for name, sample in groups.items():
        validate_sample_summary(sample, f"summary cell statistic {name}")
    exact_keys(
        aggregates, (metric for metric, _, _ in METRICS),
        "summary aggregate metrics",
    )
    for metric, _, _ in METRICS:
        metric_groups = aggregates[metric]
        if type(metric_groups) is not dict:
            fail(f"summary aggregate {metric} is not an object")
        exact_keys(metric_groups, AGGREGATE_GROUPS, f"summary aggregate {metric}")
        for name, sample in metric_groups.items():
            validate_sample_summary(sample, f"summary aggregate {metric}:{name}")


def validate_binary_receipt(binary: Any, require_complete: bool) -> None:
    if type(binary) is not dict:
        fail("summary binary is not an object")
    exact_keys(binary, BINARY_RECEIPT_KEYS, "summary binary")
    if (
        type(binary["capture_error"]) is not str
        or not binary["capture_error"]
        or len(binary["capture_error"]) > 1024 * 1024
    ):
        fail("summary binary capture error differs")
    exact_absolute_path(binary["path"], "summary binary path")
    for name in ("process_started", "timed_out", "output_overflow"):
        if type(binary[name]) is not bool:
            fail(f"summary binary {name} is not bool")
    expected_hash = lower_hash(binary["expected_sha256"], "binary expected SHA-256")
    for name in ("sha256_before", "sha256_after"):
        if binary[name] is not None:
            lower_hash(binary[name], f"binary {name}")
    start = binary["child_start_monotonic_ns"]
    reap = binary["child_reap_monotonic_ns"]
    execution_start = binary["execution_started_monotonic_ns"]
    execution_finish = binary["execution_finished_monotonic_ns"]
    if (execution_start is None) != (execution_finish is None):
        fail("binary execution interval is incomplete")
    if execution_start is not None and execution_finish is not None:
        exact_int(
            execution_start, 0, MAX_INT63,
            "binary execution start timestamp",
        )
        exact_int(
            execution_finish, execution_start, MAX_INT63,
            "binary execution finish timestamp",
        )
    if binary["process_started"]:
        if binary["child_pid"] is None or start is None:
            fail("started binary lacks its PID/start timestamp")
        exact_int(binary["child_pid"], 1, MAX_INT63, "binary child PID")
        if binary["child_process_start_ticks"] is not None:
            exact_int(
                binary["child_process_start_ticks"], 1, MAX_UINT64,
                "binary child process start ticks",
            )
        affinity = binary["child_affinity_after_spawn"]
        if affinity:
            exact_int_list(
                affinity, "binary child affinity after spawn",
                sorted_unique=True,
            )
            if binary["child_process_start_ticks"] is None:
                fail("binary affinity lacks prior process-start identity")
        exact_int(start, 0, MAX_INT63, "binary child start timestamp")
        if reap is not None:
            exact_int(
                reap, start if start is not None else 0, MAX_INT63,
                "binary child reap timestamp",
            )
        if (
            execution_start is None or execution_finish is None
            or start is not None and start < execution_start
            or reap is not None and reap > execution_finish
        ):
            fail("binary child interval escapes its execution interval")
    elif (
        start is not None or reap is not None or binary["child_pid"] is not None
        or binary["child_process_start_ticks"] is not None
        or binary["child_affinity_after_spawn"]
        or binary["timed_out"] or binary["output_overflow"]
        or binary["capture_error"] != "none"
    ):
        fail("unstarted binary has child/capture process evidence")
    for name in ("stat_before", "stat_after", "path_stat_after"):
        receipt = binary[name]
        if type(receipt) is not dict:
            fail(f"binary {name} is not an object")
        if receipt:
            allow_unlinked_held_file = (
                name == "stat_after"
                or name == "stat_before" and not binary["process_started"]
            )
            validate_stat_receipt(
                receipt, f"binary {name}",
                minimum_nlink=0 if allow_unlinked_held_file else 1,
            )
    if binary["sha256_before"] is not None and (
        not binary["stat_before"]
        or binary["stat_before"]["nlink"] != 1
        or binary["stat_before"]["mode"] != 0o555
        or not 1 <= binary["stat_before"]["size"] <= MAX_BINARY_BYTES
    ):
        fail("binary pre-exec hash lacks its prior path/stat policy")
    if binary["sha256_after"] is not None and (
        not binary["stat_after"]
        or not 1 <= binary["stat_after"]["size"] <= MAX_BINARY_BYTES
    ):
        fail("binary post-exec hash lacks its prior held-file stat policy")
    if binary["path_stat_after"] and binary["sha256_after"] is None:
        fail("binary post-exec path stat lacks its prior held-file hash")
    if binary["process_started"] and (
        not binary["stat_before"]
        or binary["stat_before"]["nlink"] != 1
        or binary["stat_before"]["mode"] != 0o555
        or binary["stat_before"]["size"] < 1
        or binary["stat_before"]["size"] > MAX_BINARY_BYTES
        or binary["sha256_before"] != expected_hash
    ):
        fail("started binary lacks its exact pre-exec image seal")
    if require_complete:
        if (
            not binary["process_started"]
            or binary["child_pid"] is None
            or binary["child_process_start_ticks"] is None
            or binary["child_affinity_after_spawn"] != [EXPECTED_TARGET_CPU]
            or execution_start is None or execution_finish is None
            or start is None
            or binary["timed_out"]
            or binary["output_overflow"]
            or binary["capture_error"] != "none"
            or reap is None
            or binary["sha256_before"] != expected_hash
            or binary["sha256_after"] != expected_hash
            or not binary["stat_before"]
            or not canonical_equal(binary["stat_before"], binary["stat_after"])
            or not canonical_equal(binary["stat_after"], binary["path_stat_after"])
        ):
            fail("complete binary receipt differs")
        info = binary["stat_before"]
        if info["nlink"] != 1 or info["mode"] != 0o555:
            fail("complete binary executable policy requires exact mode 0555")


def binary_execution_elapsed(binary: Mapping[str, Any]) -> float:
    if not binary:
        return 0.0
    start = binary.get("execution_started_monotonic_ns")
    finish = binary.get("execution_finished_monotonic_ns")
    if start is None and finish is None:
        return 0.0
    if type(start) is not int or type(finish) is not int or finish < start:
        fail("binary execution interval cannot derive elapsed time")
    return (finish - start) / 1_000_000_000.0


def validate_execution_reap_binding(
    binary: Mapping[str, Any], returncode: Optional[int], where: str,
) -> None:
    if binary.get("process_started"):
        if (
            returncode is None
            or binary.get("child_reap_monotonic_ns") is None
        ):
            fail(f"{where} lacks confirmed worker reap evidence")
    elif returncode is not None:
        fail(f"{where} has an exit code for an unstarted worker")


def exact_int_list(
    value: Any, where: str, *, length: Optional[int] = None,
    sorted_unique: bool = False,
) -> List[int]:
    if type(value) is not list or (length is not None and len(value) != length):
        fail(f"{where} list shape differs")
    for index, item in enumerate(value):
        exact_int(item, 0, MAX_INT63, f"{where}[{index}]")
    if sorted_unique and value != sorted(set(value)):
        fail(f"{where} is not sorted and unique")
    return value


def validate_bounded_messages(value: Any, where: str) -> List[str]:
    if (
        type(value) is not list
        or len(value) > MAX_PUBLICATION_FAILURES
        or any(type(item) is not str or not item or len(item) > 4000 for item in value)
        or len(set(value)) != len(value)
    ):
        fail(f"{where} message roster differs")
    return value


def validate_sampler_receipt(sampler: Any) -> None:
    if type(sampler) is not dict:
        fail("health sampler receipt is not an object")
    exact_keys(sampler, SAMPLER_RECEIPT_KEYS, "health sampler")
    exact_string(sampler["schema"], SAMPLER_SCHEMA, "sampler schema")
    exact_string(sampler["terminal_status"], "ok", "sampler status")
    exact_int(sampler["pid"], 1, MAX_INT63, "sampler pid")
    exact_int(
        sampler["cpu"], EXPECTED_SAMPLER_CPU, EXPECTED_SAMPLER_CPU,
        "sampler cpu",
    )
    for field in (
        "process_start_ticks", "csv_device", "csv_inode", "executable_device",
        "executable_inode", "script_device", "script_inode",
    ):
        exact_int(sampler[field], 0 if field.endswith("device") else 1,
                  MAX_UINT64, f"sampler {field}")
    exact_int(sampler["process_uid"], 0, MAX_UINT32, "sampler process UID")
    exact_int(sampler["process_gid"], 0, MAX_UINT32, "sampler process GID")
    if exact_int_list(
        sampler["process_affinity"], "sampler process affinity",
        length=1, sorted_unique=True,
    ) != [EXPECTED_SAMPLER_CPU]:
        fail("sampler process affinity differs")
    process_security = sampler["process_security"]
    if type(process_security) is not dict or type(process_security.get("groups")) is not list:
        fail("sampler process security is not an object")
    validate_process_security(
        process_security, sampler["process_uid"], sampler["process_gid"],
        process_security["groups"], "sampler process security",
    )
    for name in ("csv_stat", "executable_stat", "script_stat"):
        validate_stat_receipt(sampler[name], "sampler " + name)
    if type(sampler["evidence_parent"]) is not dict:
        fail("sampler evidence parent is not an object")
    exact_keys(sampler["evidence_parent"], SAMPLER_PARENT_KEYS, "sampler evidence parent")
    exact_absolute_path(
        sampler["evidence_parent"]["path"], "sampler evidence parent path"
    )
    validate_stat_receipt(
        sampler["evidence_parent"]["stat"], "sampler evidence parent stat"
    )
    if (
        sampler["evidence_parent"]["stat"]["mode"] != 0o700
        or sampler["evidence_parent"]["stat"]["uid"] != sampler["process_uid"]
        or sampler["evidence_parent"]["stat"]["gid"] != sampler["process_gid"]
    ):
        fail("sampler evidence parent mode/owner differs")
    for name, mode, maximum in (
        ("pid_file", 0o444, 64),
        ("validation_jsonl", 0o444, MAX_SAMPLER_VALIDATION_BYTES),
        ("receipt_file", 0o444, MAX_SAMPLER_RECEIPT_BYTES),
    ):
        artifact = sampler[name]
        if type(artifact) is not dict or type(artifact.get("path")) is not str:
            fail("sampler " + name + " is not an artifact object")
        validate_sampler_artifact_receipt(
            artifact, Path(artifact["path"]), mode,
            sampler["process_uid"], sampler["process_gid"], maximum,
            "sampler " + name,
        )
    if (
        sampler["pid_file"]["bytes"] != len(str(sampler["pid"])) + 1
        or sampler["pid_file"]["sha256"]
        != sha256_bytes((str(sampler["pid"]) + "\n").encode("ascii"))
        or sampler["receipt_file"]["bytes"] != 0
        or sampler["receipt_file"]["sha256"] != sha256_bytes(b"")
    ):
        fail("sampler live PID/receipt artifact state differs")
    cmdline_argv = sampler["cmdline_argv"]
    if (
        type(cmdline_argv) is not list or not cmdline_argv
        or len(cmdline_argv) > 64
        or any(type(item) is not str or not item or len(item) > 4096
               for item in cmdline_argv)
    ):
        fail("sampler command-line argv differs")
    if sha256_bytes(
        b"".join(os.fsencode(item) + b"\0" for item in cmdline_argv)
    ) != sampler["cmdline_sha256"]:
        fail("sampler command-line argv/hash binding differs")
    validation_header_ascii = sampler["validation_header_ascii"]
    if type(validation_header_ascii) is not str:
        fail("sampler validation header is not ASCII text")
    try:
        validation_header_raw = validation_header_ascii.encode("ascii")
    except UnicodeEncodeError:
        fail("sampler validation header is not ASCII")
    validation_header = parse_canonical_json_line(
        validation_header_raw, "sampler validation header"
    )
    if (
        set(validation_header) != {
            "expected_output_owner_uid", "raw_columns",
            "sampler_source_expected_sha256", "sampling", "schema", "thresholds",
        }
        or validation_header.get("schema") != THERMAL_VALIDATION_STREAM_SCHEMA
        or validation_header.get("raw_columns") != list(THERMAL_HEADER)
        or validation_header.get("expected_output_owner_uid")
        != sampler["process_uid"]
        or validation_header.get("sampler_source_expected_sha256")
        != sampler["script_sha256"]
        or validation_header.get("sampling") != {
            "dimm_attempts": 5, "dimm_retry_delay_s": 0.01,
            "interval_s": 1.0,
        }
        or validation_header.get("thresholds") != THERMAL_SAMPLER_THRESHOLDS
    ):
        fail("sampler validation header contract differs")
    if (
        sampler["csv_stat"]["device"] != sampler["csv_device"]
        or sampler["csv_stat"]["inode"] != sampler["csv_inode"]
        or sampler["csv_stat"]["mode"] != 0o600
        or sampler["csv_stat"]["nlink"] != 1
        or sampler["csv_stat"]["uid"] != sampler["process_uid"]
        or sampler["csv_stat"]["gid"] != sampler["process_gid"]
        or not 1 <= sampler["csv_stat"]["size"] <= MAX_SAMPLER_CSV_BYTES
        or sampler["script_stat"]["device"] != sampler["script_device"]
        or sampler["script_stat"]["inode"] != sampler["script_inode"]
        or sampler["script_stat"]["mode"] != 0o444
        or sampler["script_stat"]["nlink"] != 1
        or sampler["script_stat"]["uid"] != sampler["process_uid"]
        or sampler["script_stat"]["gid"] != sampler["process_gid"]
        or sampler["script_stat"]["size"] > MAX_SOURCE_FILE_BYTES
        or sampler["executable_stat"]["device"]
        != sampler["executable_device"]
        or sampler["executable_stat"]["inode"]
        != sampler["executable_inode"]
        or sampler["executable_stat"]["mode"] != 0o755
        or sampler["executable_stat"]["uid"] != 0
        or sampler["executable_stat"]["gid"] != 0
        or not 1 <= sampler["executable_stat"]["size"] <= MAX_BINARY_BYTES
    ):
        fail("sampler file stat/owner/mode binding differs")
    start = exact_int(
        sampler["window_start_monotonic_ns"], 0, MAX_INT63,
        "sampler window start",
    )
    exact_int(
        sampler["window_end_monotonic_ns"], start, MAX_INT63,
        "sampler window end",
    )
    for field in ("script_path", "csv_path", "executable_path"):
        exact_absolute_path(sampler[field], f"sampler {field}")
    artifact_paths = [
        sampler["csv_path"], sampler["pid_file"]["path"],
        sampler["validation_jsonl"]["path"], sampler["receipt_file"]["path"],
    ]
    if (
        len(set(artifact_paths)) != 4
        or any(str(Path(path).parent) != sampler["evidence_parent"]["path"]
               for path in artifact_paths)
    ):
        fail("sampler evidence path/parent roster differs")
    for field in (
        "script_sha256", "cmdline_sha256", "executable_sha256",
        "environ_sha256", "environment_sha256", "csv_sha256",
    ):
        lower_hash(sampler[field], f"sampler {field}")
    exact_int(
        sampler["csv_bytes"], 1, sampler["csv_stat"]["size"],
        "sampler CSV complete-prefix bytes",
    )
    if (
        not canonical_equal(sampler["environment"], SEALED_LAUNCH_ENVIRONMENT)
        or sampler["environment_sha256"]
        != sha256_bytes(canonical_bytes(SEALED_LAUNCH_ENVIRONMENT))
    ):
        fail("sampler launch environment binding differs")


def validate_sampler_growth_binding(
    bound: Mapping[str, Any], current: Mapping[str, Any], where: str,
) -> None:
    validate_sampler_receipt(bound)
    validate_sampler_receipt(current)
    bound_outer = dict(bound)
    current_outer = dict(current)
    bound_start = bound_outer.pop("window_start_monotonic_ns")
    bound_end = bound_outer.pop("window_end_monotonic_ns")
    current_start = current_outer.pop("window_start_monotonic_ns")
    current_end = current_outer.pop("window_end_monotonic_ns")
    bound_outer.pop("csv_sha256")
    current_outer.pop("csv_sha256")
    bound_prefix_size = bound_outer.pop("csv_bytes")
    current_prefix_size = current_outer.pop("csv_bytes")
    bound_csv = dict(bound_outer.pop("csv_stat"))
    current_csv = dict(current_outer.pop("csv_stat"))
    bound_validation = dict(bound_outer.pop("validation_jsonl"))
    current_validation = dict(current_outer.pop("validation_jsonl"))
    if (
        current_start > current_end
        or not (
            bound_start == 0 and bound_end == 0
            or bound_start == current_start and bound_end == current_end
        )
        or not canonical_equal(bound_outer, current_outer)
    ):
        fail(f"{where} immutable sampler attestation changed")
    bound_size = bound_csv.pop("size")
    current_size = current_csv.pop("size")
    bound_mtime = bound_csv.pop("mtime_ns")
    current_mtime = current_csv.pop("mtime_ns")
    if (
        not canonical_equal(bound_csv, current_csv)
        or current_size < bound_size
        or current_mtime < bound_mtime
        or current_prefix_size < bound_prefix_size
        or current_prefix_size > current_size
    ):
        fail(f"{where} sampler CSV identity regressed or changed")
    bound_validation_stat = dict(bound_validation.pop("stat"))
    current_validation_stat = dict(current_validation.pop("stat"))
    bound_validation_bytes = bound_validation.pop("bytes")
    current_validation_bytes = current_validation.pop("bytes")
    bound_validation.pop("sha256")
    current_validation.pop("sha256")
    bound_validation_size = bound_validation_stat.pop("size")
    current_validation_size = current_validation_stat.pop("size")
    bound_validation_mtime = bound_validation_stat.pop("mtime_ns")
    current_validation_mtime = current_validation_stat.pop("mtime_ns")
    if (
        not canonical_equal(bound_validation, current_validation)
        or not canonical_equal(bound_validation_stat, current_validation_stat)
        or current_validation_size < bound_validation_size
        or current_validation_mtime < bound_validation_mtime
        or bound_validation_bytes != bound_validation_size
        or current_validation_bytes != current_validation_size
    ):
        fail(f"{where} sampler validation stream identity regressed or changed")


def validate_thermal_receipt(
    thermal: Any, sampler: Mapping[str, Any],
    child_start: Optional[int], child_reap: Optional[int],
    *, allow_child_shortfall: bool = False,
) -> None:
    if type(thermal) is not dict:
        fail("health thermal receipt is not an object")
    exact_keys(thermal, THERMAL_RECEIPT_KEYS, "health thermal receipt")
    exact_int(
        thermal["cpu"], EXPECTED_SAMPLER_CPU, EXPECTED_SAMPLER_CPU,
        "thermal CPU",
    )
    for name, lower in (
        ("pid", 1), ("process_start_ticks", 1), ("csv_device", 0),
        ("csv_inode", 1), ("sample_count", 0),
        ("valid_sample_count", 0), ("invalid_sample_count", 0),
        ("window_csv_bytes", 1), ("window_start_monotonic_ns", 0),
        ("window_end_monotonic_ns", 0),
        ("validation_attempt_errors_total", 0),
        ("validation_device", 0), ("validation_inode", 1),
        ("validation_jsonl_bytes", 1), ("validation_sample_count", 1),
    ):
        exact_int(thermal[name], lower, MAX_UINT64, f"thermal {name}")
    if (
        thermal["sample_count"]
        != thermal["valid_sample_count"] + thermal["invalid_sample_count"]
    ):
        fail("thermal valid/invalid sample accounting differs")
    for name in ("cpu_tctl_max_millic", "dimm_max_millic"):
        if thermal[name] is not None:
            exact_int(thermal[name], -(1 << 63), MAX_INT63, f"thermal {name}")
    for name in ("dimm_read_errors", "edac_ce_max", "edac_ue_max"):
        if thermal[name] is not None:
            exact_int(thermal[name], 0, MAX_UINT64, f"thermal {name}")
    validate_bounded_messages(thermal["parse_failures"], "thermal parse failures")
    validate_bounded_messages(
        thermal["validation_failures"], "thermal validation failures"
    )
    exact_string(thermal["terminal_status"], "complete", "thermal status")
    exact_absolute_path(thermal["script_path"], "thermal script path")
    exact_absolute_path(thermal["csv_path"], "thermal CSV path")
    lower_hash(thermal["script_sha256"], "thermal script SHA-256")
    digest = lower_hash(thermal["window_csv_sha256"], "thermal window SHA-256")
    window_ascii = thermal["window_csv_ascii"]
    if type(window_ascii) is not str:
        fail("thermal raw window is not text")
    try:
        window_bytes = window_ascii.encode("ascii")
    except UnicodeEncodeError:
        fail("thermal raw window is not ASCII")
    if (
        len(window_bytes) != thermal["window_csv_bytes"]
        or len(window_bytes) > MAX_THERMAL_WINDOW_BYTES
        or sha256_bytes(window_bytes) != digest
    ):
        fail("thermal raw window binding differs")
    validation_ascii = thermal["validation_jsonl_ascii"]
    if type(validation_ascii) is not str:
        fail("thermal validation window is not text")
    try:
        validation_bytes = validation_ascii.encode("ascii")
    except UnicodeEncodeError:
        fail("thermal validation window is not ASCII")
    if (
        len(validation_bytes) != thermal["validation_jsonl_bytes"]
        or len(validation_bytes) > MAX_SAMPLER_VALIDATION_BYTES
        or sha256_bytes(validation_bytes)
        != lower_hash(
            thermal["validation_jsonl_sha256"],
            "thermal validation window SHA-256",
        )
    ):
        fail("thermal validation window binding differs")
    exact_absolute_path(thermal["validation_path"], "thermal validation path")
    physical = window_bytes.splitlines(keepends=True)
    if not physical or tuple(
        parse_csv_physical_line(physical[0], "thermal receipt header")
    ) != THERMAL_HEADER:
        fail("thermal receipt header differs")
    recomputed = summarize_thermal_window_bytes(
        window_bytes, thermal["window_start_monotonic_ns"],
        thermal["window_end_monotonic_ns"],
    )
    for name in (
        "cpu_tctl_max_millic", "dimm_max_millic", "dimm_read_errors",
        "edac_ce_max", "edac_ue_max", "invalid_sample_count",
        "parse_failures", "sample_count", "valid_sample_count",
    ):
        if not canonical_equal(thermal[name], recomputed[name]):
            fail(f"thermal raw-window recomputation differs at {name}")
    validation_recomputed = summarize_validation_window_bytes(
        validation_bytes, window_bytes, sampler
    )
    for name, expected_name in (
        ("validation_attempt_errors_total", "attempt_errors_total"),
        ("validation_failures", "failures"),
        ("validation_sample_count", "sample_count"),
    ):
        if not canonical_equal(thermal[name], validation_recomputed[expected_name]):
            fail(f"thermal validation-window recomputation differs at {name}")
    validation_artifact = sampler["validation_jsonl"]
    if (
        thermal["validation_path"] != validation_artifact["path"]
        or thermal["validation_device"] != validation_artifact["stat"]["device"]
        or thermal["validation_inode"] != validation_artifact["stat"]["inode"]
        or thermal["validation_sample_count"] != thermal["sample_count"]
    ):
        fail("thermal validation/sampler identity differs")
    for name in (
        "pid", "cpu", "process_start_ticks", "script_path", "script_sha256",
        "csv_path", "csv_device", "csv_inode", "window_start_monotonic_ns",
        "window_end_monotonic_ns",
    ):
        if not canonical_equal(thermal[name], sampler[name]):
            fail(f"thermal/sampler identity differs at {name}")
    if (child_start is None) != (child_reap is None):
        fail("thermal child interval is incomplete")
    if child_start is not None and child_reap is not None:
        covers_child = bool(
            thermal["window_start_monotonic_ns"] <= child_start
            <= child_reap <= thermal["window_end_monotonic_ns"]
        )
        authentic_shortfall = bool(
            allow_child_shortfall
            and thermal["window_start_monotonic_ns"] <= child_start
            and thermal["window_end_monotonic_ns"] < child_reap
        )
        if not covers_child and not authentic_shortfall:
            fail("thermal window does not cover the complete child interval")


def validate_sibling_receipt(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(
        value,
        {
            "cpu", "delta_non_idle_ticks", "end", "interval_end_monotonic_ns",
            "interval_start_monotonic_ns", "start",
        },
        where,
    )
    exact_int(value["cpu"], 56, 56, f"{where}.cpu")
    start_ns = exact_int(
        value["interval_start_monotonic_ns"], 0, MAX_INT63,
        f"{where}.interval_start",
    )
    end_ns = exact_int(
        value["interval_end_monotonic_ns"], start_ns, MAX_INT63,
        f"{where}.interval_end",
    )
    validate_cpu_tick_receipt(value["start"], where + " start")
    validate_cpu_tick_receipt(value["end"], where + " end")
    recomputed = unbounded_sibling_tick_receipt(
        value["start"], value["end"], start_ns, end_ns
    )
    exact_int(
        value["delta_non_idle_ticks"], recomputed["delta_non_idle_ticks"],
        recomputed["delta_non_idle_ticks"], f"{where}.delta",
    )


def validate_health_receipt(
    health: Any, binary: Mapping[str, Any], require_complete: bool = False,
) -> None:
    if type(health) is not dict:
        fail("summary health is not an object")
    exact_keys(health, HEALTH_KEYS, "summary health")
    exact_string(health["schema"], HEALTH_SCHEMA, "health schema")
    collection = validate_bounded_messages(
        health["collection_failures"], "health collection failures"
    )
    violations = validate_bounded_messages(
        health["violations"], "health violations"
    )
    if health["evidence_status"] not in ("complete", "partial"):
        fail("health evidence status differs")
    if health["terminal_status"] not in ("ok", "invalid"):
        fail("health terminal status differs")
    exact_int(
        health["target_cpu"], EXPECTED_TARGET_CPU, EXPECTED_TARGET_CPU,
        "health target CPU",
    )
    exact_int(
        health["controller_cpu"], EXPECTED_CONTROLLER_CPU,
        EXPECTED_CONTROLLER_CPU, "health controller CPU",
    )
    exact_int(
        health["sampler_cpu"], EXPECTED_SAMPLER_CPU,
        EXPECTED_SAMPLER_CPU, "health sampler CPU",
    )
    exact_int(
        health["sibling_non_idle_tick_cap"], SIBLING_NON_IDLE_TICK_CAP,
        SIBLING_NON_IDLE_TICK_CAP, "health sibling cap",
    )
    exact_int(
        health["thermal_max_millic"], THERMAL_MAX_MILLIC,
        THERMAL_MAX_MILLIC, "health thermal cap",
    )
    exact_string(
        health["edac_policy"], "ce-and-ue-every-sample-zero-v1",
        "health EDAC policy",
    )
    exact_string(
        health["sibling_tick_policy"],
        "linux-proc-stat-user-nice-system-irq-softirq-steal-v1",
        "health sibling tick policy",
    )
    child_start = health["child_start_monotonic_ns"]
    child_reap = health["child_reap_monotonic_ns"]
    if child_start is not None:
        exact_int(child_start, 0, MAX_INT63, "health child start")
    if child_reap is not None:
        exact_int(
            child_reap, child_start if child_start is not None else 0,
            MAX_INT63, "health child reap",
        )
    if (
        binary.get("child_start_monotonic_ns") != child_start
        or binary.get("child_reap_monotonic_ns") != child_reap
    ):
        fail("health child timestamps differ from binary")

    topology_fields = (
        ("target_core", EXPECTED_TARGET_CORE),
        ("controller_core", EXPECTED_CONTROLLER_CORE),
        ("sampler_core", EXPECTED_SAMPLER_CORE),
        ("target_threads", EXPECTED_TARGET_THREADS),
    )
    for name, expected in topology_fields:
        value = health[name]
        if value:
            exact_int_list(value, "health " + name, length=len(expected))
            if tuple(value) != tuple(expected):
                fail(f"health {name} differs")
        elif require_complete:
            fail(f"complete health lacks {name}")
    initial_affinity = health["controller_initial_affinity"]
    if initial_affinity:
        exact_int_list(
            initial_affinity, "health initial affinity", sorted_unique=True
        )
        if (
            EXPECTED_TARGET_CPU not in initial_affinity
            or EXPECTED_CONTROLLER_CPU not in initial_affinity
            or EXPECTED_SAMPLER_CPU not in initial_affinity
        ):
            fail("health initial affinity lacks frozen CPUs")
    elif require_complete:
        fail("complete health lacks initial affinity")
    end_affinity = health["controller_singleton_affinity_end"]
    drift_marker = "controller CPU121 affinity changed"
    affinity_marker_interrupted = bool(end_affinity) and any(
        message.startswith("controller affinity end:")
        for message in collection
    )
    if end_affinity:
        exact_int_list(
            end_affinity, "health controller affinity end",
            sorted_unique=True,
        )
        affinity_drift = end_affinity != [EXPECTED_CONTROLLER_CPU]
        drift_marker_present = drift_marker in violations
        if (
            not affinity_drift and drift_marker_present
            or affinity_drift
            and not affinity_marker_interrupted
            and not drift_marker_present
        ):
            fail("health controller affinity drift/violation binding differs")
        if require_complete and affinity_drift:
            fail("complete health controller affinity end differs")
    else:
        if drift_marker in violations:
            fail("missing controller affinity has a false drift violation")
        if require_complete:
            fail("complete health lacks controller affinity end")

    admission = health["admission_sibling_ticks"]
    if admission:
        validate_sibling_receipt(admission, "health admission sibling")
    elif require_complete:
        fail("complete health lacks admission sibling receipt")
    siblings = health["sibling_ticks"]
    run_cap_marker = "CPU56 run non-idle tick cap exceeded"
    run_marker_interrupted = bool(siblings) and any(
        message.startswith("run sibling ticks:") for message in collection
    )
    if type(siblings) is not list or len(siblings) > 1:
        fail("health run sibling roster differs")
    if siblings:
        validate_sibling_receipt(siblings[0], "health run sibling")
        if siblings[0]["interval_start_monotonic_ns"] != child_start:
            fail("run sibling start differs from child start")
        if siblings[0]["interval_end_monotonic_ns"] != child_reap:
            fail("run sibling end differs from child reap")
        execution_start = binary.get("execution_started_monotonic_ns")
        execution_finish = binary.get("execution_finished_monotonic_ns")
        if (
            execution_start is not None
            and siblings[0]["start"]["read_monotonic_ns"] < execution_start
            or execution_finish is not None
            and siblings[0]["end"]["read_monotonic_ns"] > execution_finish
        ):
            fail("run sibling tick reads escape the execution interval")
        run_cap_exceeded = (
            siblings[0]["delta_non_idle_ticks"]
            > SIBLING_NON_IDLE_TICK_CAP
        )
        run_cap_marker_present = run_cap_marker in violations
        if (
            not run_cap_exceeded and run_cap_marker_present
            or run_cap_exceeded
            and not run_marker_interrupted
            and not run_cap_marker_present
        ):
            fail("health run sibling cap/violation binding differs")
    elif require_complete:
        fail("complete health lacks run sibling receipt")
    elif run_cap_marker in violations:
        fail("missing run sibling receipt has a false cap violation")

    sampler_admission = health["sampler_admission"]
    if sampler_admission:
        validate_sampler_receipt(sampler_admission)
        if (
            sampler_admission["window_start_monotonic_ns"] != 0
            or sampler_admission["window_end_monotonic_ns"] != 0
        ):
            fail("health sampler admission window differs")
        admission_digest = lower_hash(
            health["sampler_admission_receipt_sha256"],
            "sampler admission receipt SHA-256",
        )
        if sha256_bytes(canonical_bytes(sampler_admission)) != admission_digest:
            fail("sampler admission receipt binding differs")
    elif health["sampler_admission_receipt_sha256"] is not None:
        fail("partial health has an unbound sampler admission digest")
    sampler = health["sampler"]
    sampler_terminal = health["sampler_terminal"]
    thermal = health["thermal"]
    if sampler:
        validate_sampler_receipt(sampler)
        execution_start = binary.get("execution_started_monotonic_ns")
        execution_finish = binary.get("execution_finished_monotonic_ns")
        if (
            execution_start is not None
            and sampler["window_start_monotonic_ns"] > execution_start
        ):
            fail("health sampler start follows binary execution start")
        if (
            execution_finish is not None
            and sampler["window_end_monotonic_ns"] > execution_finish
        ):
            fail("health sampler endpoint follows binary execution finish")
        digest = lower_hash(
            health["sampler_receipt_sha256"], "sampler receipt SHA-256"
        )
        if sha256_bytes(canonical_bytes(sampler)) != digest:
            fail("sampler receipt binding differs")
        if not sampler_admission:
            fail("final sampler receipt lacks its admission receipt")
        validate_sampler_growth_binding(
            sampler_admission, sampler, "health sampler admission/final"
        )
    elif health["sampler_receipt_sha256"] is not None:
        fail("partial health has an unbound sampler digest")
    if sampler_terminal:
        if sampler or thermal or not sampler_admission:
            fail("terminal sampler evidence overlaps live sampler evidence")
        validate_sampler_terminal_receipt(
            sampler_terminal, sampler_admission, child_start, child_reap
        )
        terminal_digest = lower_hash(
            health["sampler_terminal_receipt_sha256"],
            "sampler terminal receipt SHA-256",
        )
        if sha256_bytes(canonical_bytes(sampler_terminal)) != terminal_digest:
            fail("sampler terminal receipt binding differs")
        execution_finish = binary.get("execution_finished_monotonic_ns")
        if (
            execution_finish is not None
            and sampler_terminal["process_exit_observed_monotonic_ns"]
            > execution_finish
        ):
            fail("sampler terminal exit observation follows execution finish")
        if not any(message.startswith("sampler endpoint:") for message in collection):
            fail("terminal sampler evidence lacks its endpoint diagnostic")
    elif health["sampler_terminal_receipt_sha256"] is not None:
        fail("partial health has an unbound sampler terminal digest")
    if thermal:
        if not sampler:
            fail("thermal receipt lacks sampler identity")
        child_shortfall_markers = {
            (
                "sampler endpoint: sampler-monitor-invalid:"
                "validation-heartbeat-stale"
            ),
            (
                "sampler endpoint: sampler-monitor-invalid:"
                "coverage-recovery-cutoff"
            ),
        }
        validate_thermal_receipt(
            thermal, sampler, child_start, child_reap,
            allow_child_shortfall=bool(
                child_shortfall_markers.intersection(collection)
            ),
        )
    thermal_collection_markers = sum(
        message.startswith("thermal collection:") for message in collection
    )
    expected_thermal_collection_markers = 1 if sampler and not thermal else 0
    if thermal_collection_markers != expected_thermal_collection_markers:
        fail("health thermal collection/receipt binding differs")
    # finish_health starts with a clean preclaim violation roster and is the
    # only postclaim producer.  Affinity/run markers are already checked for
    # causal legality above; sampler/thermal publication is transactional.
    # The complete ordered roster is therefore exact even when thermal
    # collection fails before producing a receipt.
    expected_violations: List[str] = []
    if drift_marker in violations:
        append_unique_failure(expected_violations, drift_marker)
    if run_cap_marker in violations:
        append_unique_failure(expected_violations, run_cap_marker)
    if thermal:
        for message in thermal_policy_violation_messages(thermal):
            append_unique_failure(expected_violations, message)
    if violations != expected_violations:
        fail("health policy/violation roster differs")
    complete = bool(
        not collection and child_start is not None and child_reap is not None
        and admission and len(siblings) == 1 and sampler_admission
        and sampler and not sampler_terminal and thermal
        and all(health[name] for name, _ in topology_fields)
        and initial_affinity and end_affinity
    )
    if health["evidence_status"] != ("complete" if complete else "partial"):
        fail("health computed evidence status differs")
    expected_terminal = "ok" if complete and not violations else "invalid"
    if health["terminal_status"] != expected_terminal:
        fail("health computed terminal status differs")
    receipt_hash = lower_hash(health["receipt_sha256"], "health receipt SHA-256")
    preimage = dict(health)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != receipt_hash:
        fail("health receipt SHA-256 differs")
    if require_complete:
        if (
            not complete or collection or violations
            or health["terminal_status"] != "ok"
            or siblings[0]["delta_non_idle_ticks"] > SIBLING_NON_IDLE_TICK_CAP
            or thermal["sample_count"] < 2
            or thermal["invalid_sample_count"] != 0
            or thermal["parse_failures"]
            or thermal["validation_failures"]
            or thermal["cpu_tctl_max_millic"] is None
            or thermal["cpu_tctl_max_millic"] > THERMAL_MAX_MILLIC
            or thermal["dimm_max_millic"] is None
            or thermal["dimm_max_millic"] > THERMAL_MAX_MILLIC
            or thermal["dimm_read_errors"] != 0
            or thermal["edac_ce_max"] != 0
            or thermal["edac_ue_max"] != 0
        ):
            fail("valid outcome lacks clean complete health evidence")


def parse_summary_bytes(data: bytes) -> Dict[str, Any]:
    if not data or len(data) > 8 * 1024 * 1024 or not data.endswith(b"\n"):
        fail("summary is not one bounded LF-terminated object")
    if b"\r" in data or data.count(b"\n") != 1:
        fail("summary contains a noncanonical line ending")
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(f"summary is malformed: {exc}")
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail("summary is not canonical JSON")
    exact_keys(value, SUMMARY_KEYS, "summary")
    exact_string(value["schema"], SUMMARY_SCHEMA, "summary schema")
    digest = lower_hash(value["summary_preimage_sha256"], "summary preimage")
    preimage = dict(value)
    preimage["summary_preimage_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail("summary self-preimage SHA-256 differs")
    if value["outcome"] not in ("pass", "reject", "invalid"):
        fail("summary outcome differs")
    publication_failures = value["publication_failures"]
    if (
        type(publication_failures) is not list
        or len(publication_failures) > MAX_PUBLICATION_FAILURES
        or any(
            type(item) is not str or not item or len(item) > 4096
            for item in publication_failures
        )
        or len(set(publication_failures)) != len(publication_failures)
    ):
        fail("summary publication failures differ")
    if value["signal"] is not None and value["signal"] not in (
        "SIGINT", "SIGTERM", "SIGHUP"
    ):
        fail("summary signal differs")
    if (
        type(value["signal_names"]) is not list
        or len(value["signal_names"]) > len(SignalGuard.SIGNALS)
        or any(type(name) is not str or name not in (
                   "SIGINT", "SIGTERM", "SIGHUP")
               for name in value["signal_names"])
        or len(set(value["signal_names"])) != len(value["signal_names"])
        or value["signal"] != (
            value["signal_names"][0] if value["signal_names"] else None
        )
    ):
        fail("summary observed-signal roster differs")
    if type(value["raw_complete"]) is not bool:
        fail("summary raw-complete flag differs")
    for name in ("raw_bytes", "raw_record_count", "stderr_bytes", "target_cpu"):
        exact_int(value[name], 0, MAX_INT63, f"summary {name}")
    for name in ("raw_sha256", "stderr_sha256"):
        lower_hash(value[name], f"summary {name}")
    expected_commit = value["expected_source_git_commit"]
    if type(expected_commit) is not str or LOWER40.fullmatch(expected_commit) is None:
        fail("summary expected source commit differs")
    expected_manifest_sha256 = lower_hash(
        value["expected_source_manifest_sha256"],
        "summary expected source manifest",
    )
    if (
        type(value["elapsed_seconds"]) is not float
        or not math.isfinite(value["elapsed_seconds"])
        or value["elapsed_seconds"] < 0.0
        or type(value["failure"]) is not str
        or not value["failure"]
        or len(value["failure"]) > MAX_FAILURE_TEXT_CHARS
    ):
        fail("summary scalar receipt differs")
    if value["process_exit_code"] is not None:
        exact_int(
            value["process_exit_code"], -(1 << 31), (1 << 31) - 1,
            "summary process exit code",
        )
    for name in (
        "binary", "config", "git_after", "git_before", "health",
        "health_module_loader", "output_bundle", "source_manifest_after",
        "source_manifest_before", "statistics",
    ):
        if type(value[name]) is not dict:
            fail(f"summary {name} is not an object")
    validate_output_bundle_receipt(value["output_bundle"])
    for name in ("source_manifest_before", "source_manifest_after"):
        if value[name]:
            validate_source_manifest_receipt(value[name], f"summary {name}")
    for name in ("git_before", "git_after"):
        if value[name]:
            validate_git_receipt(value[name], expected_commit, f"summary {name}")
    if value["config"]:
        validate_config(value["config"], value["target_cpu"], expected_commit)
    if value["statistics"]:
        validate_statistics_receipt(value["statistics"])
    binary = value["binary"]
    if binary:
        validate_binary_receipt(binary, value["outcome"] in ("pass", "reject"))
        validate_execution_reap_binding(
            binary, value["process_exit_code"], "summary execution"
        )
    if value["elapsed_seconds"] != binary_execution_elapsed(binary):
        fail("summary elapsed time differs from its binary execution interval")
    if value["health"]:
        if not binary:
            fail("health receipt lacks its binary receipt")
        validate_health_receipt(
            value["health"], binary,
            require_complete=value["outcome"] in ("pass", "reject"),
        )
    if value["health_module_loader"]:
        validate_health_loader_receipt(value["health_module_loader"])
    if value["outcome"] in ("pass", "reject"):
        if not value["statistics"]:
            fail("valid summary lacks statistics")
        expected_failure = (
            "statistical-gates"
            if value["statistics"]["failed_gates"] else "none"
        )
        expected_outcome = (
            "reject" if value["statistics"]["failed_gates"] else "pass"
        )
        if (
            publication_failures
            or value["signal"] is not None
            or value["signal_names"]
            or not value["raw_complete"]
            or value["raw_record_count"] != EXPECTED_PANEL_COUNT + 2
            or value["process_exit_code"] != 0
            or value["stderr_bytes"] != 0
            or value["stderr_sha256"] != sha256_bytes(b"")
            or value["elapsed_seconds"] >= OUTER_DEADLINE_SECONDS
            or value["target_cpu"] != EXPECTED_TARGET_CPU
            or not binary
            or not value["health"]
            or not value["health_module_loader"]
            or not value["config"]
            or not value["statistics"]
            or not value["source_manifest_before"]
            or not value["source_manifest_after"]
            or not value["git_before"]
            or not value["git_after"]
            or value["failure"] != expected_failure
            or value["outcome"] != expected_outcome
            or value["source_manifest_before"]["sha256"]
            != expected_manifest_sha256
            or not canonical_equal(
                value["source_manifest_before"], value["source_manifest_after"]
            )
            or not canonical_equal(value["git_before"], value["git_after"])
        ):
            fail("valid summary outcome lacks complete exact execution evidence")
        manifest_entries = {
            entry["path"]: entry
            for entry in value["source_manifest_before"]["entries"]
        }
        for module in value["health_module_loader"]["modules"]:
            source = manifest_entries[module["path"]]
            if (
                module["bytes"] != source["bytes"]
                or module["sha256"] != source["sha256"]
            ):
                fail("health loader differs from the source manifest")
    return value


class OutputBundle:
    """Owner for one exact, crash-durable, immutable evidence directory."""

    def __init__(
        self, output_dir: Path, target_cpu: int, expected_commit: str,
        expected_source_manifest_sha256: str,
        expected_git_sha256: str,
        source_root: Path,
    ) -> None:
        self.output_dir = output_dir
        self.target_cpu = target_cpu
        self.expected_commit = expected_commit
        self.expected_source_manifest_sha256 = expected_source_manifest_sha256
        self.expected_git_sha256 = expected_git_sha256
        self.source_root = source_root
        self.parent_fd = -1
        self.directory_fd = -1
        self.file_fds: Dict[str, int] = {}
        self.parent_identity: Tuple[int, int, int] = (0, 0, 0)
        self.parent_initial_nlink = 0
        self.directory_identity: Tuple[int, int, int] = (0, 0, 0)
        self.file_identities: Dict[str, Tuple[int, int, int]] = {}
        self.staged_raw = b""
        self.staged_stderr = b""
        self.raw_validation_config: Dict[str, Any] = {}
        self.raw_validation_statistics: Dict[str, Any] = {}
        self.raw_validation_failures: List[str] = []
        self.expected_summary_components: Dict[str, bytes] = {}
        self.fallback_summary_components: Dict[str, bytes] = {}
        self.signal_name: Optional[str] = None
        self.signal_names: List[str] = []
        self.closed = False
        self.logically_sealed = False
        self.summary_mode_poisoned = False
        self.summary_unlinked = False

    @classmethod
    def reserve(
        cls, output_dir: Path, target_cpu: int, expected_commit: str,
        expected_source_manifest_sha256: str,
        expected_git_sha256: str,
        source_root: Path,
        signal_guard: Optional[SignalGuard] = None,
    ) -> "OutputBundle":
        require_preflight_side_effects_allowed("output reservation")
        if (
            not output_dir.is_absolute()
            or output_dir.name in ("", ".", "..")
            or os.path.normpath(str(output_dir)) != str(output_dir)
        ):
            fail("--output-dir must be one canonical absolute new directory")
        lower_hash(
            expected_source_manifest_sha256,
            "expected source manifest SHA-256",
        )
        lower_hash(expected_git_sha256, "expected Git executable SHA-256")
        if (
            not source_root.is_absolute()
            or os.path.normpath(str(source_root)) != str(source_root)
        ):
            fail("source root must be one canonical absolute directory")
        try:
            resolved_source_root = source_root.resolve(strict=True)
        except OSError as exc:
            fail(f"cannot resolve source root: {exc}")
        if resolved_source_root != source_root or not source_root.is_dir():
            fail("source root path is not its canonical directory")
        bundle = cls(
            output_dir, target_cpu, expected_commit,
            expected_source_manifest_sha256, expected_git_sha256, source_root,
        )
        created = False
        try:
            directory_flag = getattr(os, "O_DIRECTORY", 0)
            nofollow = getattr(os, "O_NOFOLLOW", 0)
            if not directory_flag or not nofollow:
                fail("output reservation requires O_DIRECTORY and O_NOFOLLOW")
            common = getattr(os, "O_CLOEXEC", 0) | nofollow
            bundle.parent_fd = os.open(
                str(output_dir.parent), os.O_RDONLY | directory_flag | common
            )
            bundle._verify_parent_path()
            os.mkdir(output_dir.name, 0o700, dir_fd=bundle.parent_fd)
            created = True
            os.fsync(bundle.parent_fd)
            bundle.directory_fd = os.open(
                output_dir.name, os.O_RDONLY | directory_flag | common,
                dir_fd=bundle.parent_fd,
            )
            os.fchmod(bundle.directory_fd, 0o700)
            directory_info = os.fstat(bundle.directory_fd)
            if not stat.S_ISDIR(directory_info.st_mode):
                fail("reserved output is not a directory")
            bundle.directory_identity = (
                directory_info.st_dev, directory_info.st_ino,
                directory_info.st_nlink,
            )
            flags = (
                os.O_RDWR | os.O_CREAT | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0) | nofollow
            )
            for name in OUTPUT_NAMES:
                fd = os.open(name, flags, 0o600, dir_fd=bundle.directory_fd)
                bundle.file_fds[name] = fd
                os.fchmod(fd, 0o600)
                info = os.fstat(fd)
                if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                    fail(f"reserved {name} is not a single-link regular file")
                bundle.file_identities[name] = (
                    info.st_dev, info.st_ino, info.st_nlink,
                )
            parent_info = os.fstat(bundle.parent_fd)
            bundle.parent_identity = (
                parent_info.st_dev, parent_info.st_ino,
                stat.S_IMODE(parent_info.st_mode),
            )
            bundle.parent_initial_nlink = parent_info.st_nlink
            bundle._verify_identities(0o700, 0o600)
            bundle._verify_roster()
            os.fsync(bundle.directory_fd)
            if signal_guard is not None:
                signal_guard.own_output_bundle(bundle)
            return bundle
        except BaseException as exc:
            bundle._cleanup_incomplete(created)
            if isinstance(exc, (ValidationError, KeyboardInterrupt, SystemExit)):
                raise
            fail(f"cannot reserve output bundle: {exception_text(exc)}")

    def _cleanup_incomplete(self, created: bool) -> None:
        for fd in self.file_fds.values():
            try:
                os.close(fd)
            except OSError:
                pass
        self.file_fds.clear()
        if self.directory_fd >= 0:
            for name in OUTPUT_NAMES:
                try:
                    os.unlink(name, dir_fd=self.directory_fd)
                except OSError:
                    pass
            try:
                os.close(self.directory_fd)
            except OSError:
                pass
            self.directory_fd = -1
        if created and self.parent_fd >= 0:
            try:
                os.rmdir(self.output_dir.name, dir_fd=self.parent_fd)
                os.fsync(self.parent_fd)
            except OSError:
                pass
        if self.parent_fd >= 0:
            try:
                os.close(self.parent_fd)
            except OSError:
                pass
            self.parent_fd = -1

    def _verify_parent_path(self) -> None:
        retained = os.fstat(self.parent_fd)
        named = os.stat(str(self.output_dir.parent), follow_symlinks=False)
        if (
            not stat.S_ISDIR(retained.st_mode)
            or not stat.S_ISDIR(named.st_mode)
            or (retained.st_dev, retained.st_ino) != (named.st_dev, named.st_ino)
        ):
            fail("output parent path/FD identity differs")

    def _verify_roster(
        self, expected_names: Sequence[str] = OUTPUT_NAMES,
    ) -> None:
        names: List[str] = []
        scan_fd = os.open(
            ".", os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=self.directory_fd,
        )
        try:
            with os.scandir(scan_fd) as entries:
                for entry in entries:
                    if len(names) >= len(expected_names):
                        fail("output directory has too many entries")
                    if not entry.is_file(follow_symlinks=False):
                        fail("output directory contains a non-regular entry")
                    names.append(entry.name)
        finally:
            os.close(scan_fd)
        if set(names) != set(expected_names) or len(names) != len(expected_names):
            fail("output directory roster differs")

    def _verify_identities(self, directory_mode: int, file_mode: int) -> None:
        self._verify_parent_path()
        parent = os.fstat(self.parent_fd)
        if (
            self.parent_identity != (0, 0, 0)
            and self.parent_identity != (
                parent.st_dev, parent.st_ino, stat.S_IMODE(parent.st_mode),
            )
        ):
            fail("output parent receipt changed")
        retained_dir = os.fstat(self.directory_fd)
        named_dir = os.stat(
            self.output_dir.name, dir_fd=self.parent_fd, follow_symlinks=False
        )
        expected_dir = self.directory_identity
        for info in (retained_dir, named_dir):
            if (
                not stat.S_ISDIR(info.st_mode)
                or (info.st_dev, info.st_ino, info.st_nlink) != expected_dir
                or stat.S_IMODE(info.st_mode) != directory_mode
            ):
                fail("output directory path/FD receipt differs")
        for name in OUTPUT_NAMES:
            retained = os.fstat(self.file_fds[name])
            named = os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
            expected = self.file_identities[name]
            for info in (retained, named):
                if (
                    not stat.S_ISREG(info.st_mode)
                    or (info.st_dev, info.st_ino, info.st_nlink) != expected
                    or stat.S_IMODE(info.st_mode) != file_mode
                ):
                    fail(f"output {name} path/FD receipt differs")

    def receipt(self) -> Dict[str, Any]:
        return {
            "directory": {
                "device": self.directory_identity[0],
                "inode": self.directory_identity[1],
                "nlink": self.directory_identity[2],
                "reserved_mode": 0o700,
                "sealed_mode": 0o500,
            },
            "files": {
                name: {
                    "device": self.file_identities[name][0],
                    "inode": self.file_identities[name][1],
                    "nlink": self.file_identities[name][2],
                    "reserved_mode": 0o600,
                    "sealed_mode": 0o400,
                }
                for name in OUTPUT_NAMES
            },
            "parent": {
                "device": self.parent_identity[0],
                "initial_nlink": self.parent_initial_nlink,
                "inode": self.parent_identity[1],
                "mode": self.parent_identity[2],
            },
            "path": str(self.output_dir),
        }

    @staticmethod
    def _write_exact(fd: int, data: bytes) -> bytes:
        os.ftruncate(fd, 0)
        os.lseek(fd, 0, os.SEEK_SET)
        view = memoryview(data)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                fail("short output write")
            view = view[written:]
        os.fsync(fd)
        info = os.fstat(fd)
        if info.st_size != len(data):
            fail("output size differs after write")
        chunks: List[bytes] = []
        offset = 0
        while offset < len(data):
            block = os.pread(fd, min(1024 * 1024, len(data) - offset), offset)
            if not block:
                fail("short output readback")
            chunks.append(block)
            offset += len(block)
        if os.pread(fd, 1, len(data)):
            fail("output has trailing bytes after readback")
        readback = b"".join(chunks)
        if readback != data:
            fail("output readback differs from requested bytes")
        return readback

    @staticmethod
    def _read_current(fd: int, limit: int) -> bytes:
        info = os.fstat(fd)
        if info.st_size < 0 or info.st_size > limit:
            fail("partial output exceeds its bounded size")
        chunks: List[bytes] = []
        offset = 0
        while offset < info.st_size:
            block = os.pread(fd, min(1024 * 1024, info.st_size - offset), offset)
            if not block:
                fail("partial output readback is short")
            chunks.append(block)
            offset += len(block)
        return b"".join(chunks)

    @staticmethod
    def _close_owned_fd(fd: int) -> None:
        os.close(fd)

    def stage_evidence(self, raw: bytes, stderr: bytes) -> List[str]:
        failures: List[str] = []
        self.raw_validation_config = {}
        self.raw_validation_statistics = {}
        self.raw_validation_failures = []
        for name, payload, limit in (
            ("raw.jsonl", raw, MAX_STDOUT_BYTES + 1),
            ("stderr.txt", stderr, MAX_STDERR_BYTES + 1),
        ):
            try:
                readback = self._write_exact(self.file_fds[name], payload)
            except BaseException as exc:
                failures.append(f"stage {name}: {exception_text(exc)}")
                try:
                    readback = self._read_current(self.file_fds[name], limit)
                except BaseException as read_exc:
                    failures.append(
                        f"read partial {name}: {exception_text(read_exc)}"
                    )
                    readback = b""
            if name == "raw.jsonl":
                self.staged_raw = readback
            else:
                self.staged_stderr = readback
        try:
            records = parse_transcript(self.staged_raw)
            (
                self.raw_validation_config,
                self.raw_validation_statistics,
                self.raw_validation_failures,
            ) = validate_transcript(
                records, self.target_cpu, self.expected_commit
            )
        except ValidationError:
            # Transcript incompleteness is scientific evidence, not a storage
            # failure.  Execution/process diagnostics adjudicate it invalid.
            pass
        except BaseException as exc:
            failures.append(f"canonical raw reparse: {exception_text(exc)}")
        return failures

    def bind_expected_summary_components(
        self, components: Mapping[str, Any],
    ) -> None:
        expected_names = {
            "binary", "config", "elapsed_seconds", "failure", "git_after",
            "git_before", "health", "health_module_loader", "outcome",
            "process_exit_code",
            "source_manifest_after", "source_manifest_before", "statistics",
            "target_cpu",
        }
        exact_keys(components, expected_names, "expected summary components")
        if self.expected_summary_components:
            fail("expected summary components were already bound")
        self.expected_summary_components = {
            name: canonical_bytes(components[name]) for name in expected_names
        }

    def update_fallback_summary_components(
        self, components: Mapping[str, Any],
    ) -> None:
        allowed = {
            "binary", "config", "elapsed_seconds", "failure", "git_after",
            "git_before", "health", "health_module_loader", "outcome",
            "process_exit_code", "source_manifest_after",
            "source_manifest_before", "statistics", "target_cpu",
        }
        if type(components) is not dict or not set(components).issubset(allowed):
            fail("fallback summary component roster differs")
        encoded = {
            name: canonical_bytes(value) for name, value in components.items()
        }
        updated = dict(self.fallback_summary_components)
        updated.update(encoded)
        self.fallback_summary_components = updated

    def _bound_component(self, name: str, default: Any) -> Any:
        encoded = self.expected_summary_components.get(name)
        if encoded is None:
            encoded = self.fallback_summary_components.get(name)
        if encoded is None:
            return default
        return json.loads(encoded.decode("ascii"))

    @staticmethod
    def _bounded_publication_failures(
        failures: Sequence[str], emergency_reason: Optional[str] = None,
    ) -> List[str]:
        bounded: List[str] = []
        values = list(failures)
        if emergency_reason is not None:
            values.append("emergency summary: " + emergency_reason)
        for value in values:
            text = value if type(value) is str else repr(value)
            if not text:
                text = "empty publication diagnostic"
            if len(text) > 4096:
                text = text[:4093] + "..."
            if text not in bounded:
                bounded.append(text)
            if len(bounded) >= MAX_PUBLICATION_FAILURES:
                break
        return bounded

    def emergency_summary(self, reason: str, failures: Sequence[str]) -> bytes:
        bounded_reason = bounded_failure_text([reason])
        binary = self._bound_component("binary", {})
        config = self._bound_component(
            "config", self.raw_validation_config
        )
        statistics = self._bound_component(
            "statistics", self.raw_validation_statistics
        )
        raw_complete = bool(
            self.raw_validation_config and self.raw_validation_statistics
        )
        publication_failures = self._bounded_publication_failures(
            failures,
            bounded_reason if (
                self.expected_summary_components
                or self.fallback_summary_components
            ) else None,
        )
        summary = make_summary_preimage({
            "binary": binary,
            "config": config,
            "elapsed_seconds": binary_execution_elapsed(binary),
            "expected_source_git_commit": self.expected_commit,
            "expected_source_manifest_sha256": (
                self.expected_source_manifest_sha256
            ),
            "failure": bounded_reason,
            "git_after": self._bound_component("git_after", {}),
            "git_before": self._bound_component("git_before", {}),
            "health": self._bound_component("health", {}),
            "health_module_loader": self._bound_component(
                "health_module_loader", {}
            ),
            "outcome": "invalid", "output_bundle": self.receipt(),
            "process_exit_code": self._bound_component(
                "process_exit_code", None
            ),
            "publication_failures": publication_failures,
            "raw_bytes": len(self.staged_raw), "raw_complete": raw_complete,
            "raw_record_count": self.staged_raw.count(b"\n"),
            "raw_sha256": sha256_bytes(self.staged_raw),
            "schema": SUMMARY_SCHEMA, "signal": self.signal_name,
            "signal_names": list(self.signal_names),
            "source_manifest_after": self._bound_component(
                "source_manifest_after", {}
            ),
            "source_manifest_before": self._bound_component(
                "source_manifest_before", {}
            ),
            "statistics": statistics,
            "stderr_bytes": len(self.staged_stderr),
            "stderr_sha256": sha256_bytes(self.staged_stderr),
            "summary_preimage_sha256": None,
            "target_cpu": self._bound_component("target_cpu", self.target_cpu),
        })
        return canonical_bytes(summary) + b"\n"

    def _validate_summary_cross(self, data: bytes) -> Dict[str, Any]:
        summary = parse_summary_bytes(data)
        if (
            summary["raw_bytes"] != len(self.staged_raw)
            or summary["raw_record_count"] != self.staged_raw.count(b"\n")
            or summary["raw_sha256"] != sha256_bytes(self.staged_raw)
            or summary["stderr_bytes"] != len(self.staged_stderr)
            or summary["stderr_sha256"] != sha256_bytes(self.staged_stderr)
            or not canonical_equal(summary["output_bundle"], self.receipt())
            or summary["expected_source_git_commit"] != self.expected_commit
            or summary["expected_source_manifest_sha256"]
            != self.expected_source_manifest_sha256
            or summary["target_cpu"] != self.target_cpu
        ):
            fail("summary cross-artifact binding differs")
        if summary["raw_complete"]:
            if (
                not self.raw_validation_config
                or not self.raw_validation_statistics
                or not canonical_equal(
                    summary["config"], self.raw_validation_config
                )
                or not canonical_equal(
                    summary["statistics"], self.raw_validation_statistics
                )
            ):
                fail("summary differs from the independently validated raw data")
            if summary["outcome"] in ("pass", "reject") and (
                summary["outcome"] != (
                    "reject" if self.raw_validation_failures else "pass"
                )
                or summary["failure"] != (
                    "statistical-gates"
                    if self.raw_validation_failures else "none"
                )
            ):
                fail("valid summary decision differs from the raw data")
        if summary["outcome"] in ("pass", "reject") and not (
            self.expected_summary_components
        ):
            fail("valid summary lacks independently bound components")
        if self.expected_summary_components:
            evidence_names = set(self.expected_summary_components) - {
                "failure", "outcome",
            }
            for name in evidence_names:
                expected = self.expected_summary_components[name]
                if canonical_bytes(summary[name]) != expected:
                    fail(f"summary independently bound component differs: {name}")
            if not summary["publication_failures"] and not summary["signal_names"]:
                for name in ("failure", "outcome"):
                    expected = self.expected_summary_components[name]
                    if canonical_bytes(summary[name]) != expected:
                        fail(
                            "summary independently bound decision differs: "
                            + name
                        )
        return summary

    @staticmethod
    def _append_unique(failures: List[str], message: str) -> bool:
        if message in failures:
            return False
        failures.append(message)
        return True

    def _invalid_summary(
        self, summary: Mapping[str, Any], failures: Sequence[str],
        signal_guard: Optional[SignalGuard],
    ) -> bytes:
        invalid = dict(summary)
        signal_names = (
            list(signal_guard.observed_signals)
            if signal_guard is not None else list(self.signal_names)
        )
        signal_name = (
            signal_guard.first_signal
            if signal_guard is not None else self.signal_name
        )
        existing_failure = str(invalid.get("failure", ""))
        generated_markers = (
            " | publication-failures:", " | guarded-signals:",
        )
        marker_offsets = [
            offset for marker in generated_markers
            if (offset := existing_failure.find(marker)) >= 0
        ]
        if marker_offsets:
            existing_failure = existing_failure[:min(marker_offsets)]
        existing_failure = existing_failure.rstrip()
        if not existing_failure or existing_failure == "none":
            existing_failure = "invalid-publication"
        publication_failures = self._bounded_publication_failures(
            list(invalid["publication_failures"]) + list(failures)
        )
        compact_suffixes: List[str] = []
        if publication_failures:
            compact_suffixes.append(
                f"publication-failures:{len(publication_failures)}"
            )
        if signal_names:
            compact_suffixes.append(
                "guarded-signals:" + ",".join(signal_names)
            )
        suffix = " | ".join(compact_suffixes)
        base_limit = MAX_FAILURE_TEXT_CHARS - (
            len(suffix) + 3 if suffix else 0
        )
        base_failure = bounded_failure_text(
            [existing_failure], limit=base_limit
        )
        invalid["outcome"] = "invalid"
        invalid["failure"] = (
            base_failure + " | " + suffix if suffix else base_failure
        )
        invalid["publication_failures"] = publication_failures
        invalid["signal"] = signal_name
        invalid["signal_names"] = signal_names
        invalid["summary_preimage_sha256"] = None
        return canonical_bytes(make_summary_preimage(invalid)) + b"\n"

    def _live_source_manifest(self) -> Dict[str, Any]:
        return source_manifest(self.source_root, time.monotonic() + 1.0)

    def _live_git_receipt(self) -> Dict[str, Any]:
        return git_receipt(
            self.source_root, self.expected_commit, self.expected_git_sha256,
            time.monotonic() + 1.0,
        )

    def _recheck_positive_source_git(
        self, summary: Mapping[str, Any],
    ) -> None:
        # Sandwich the final Git receipt between two complete source reads.
        # This bounds a source mutation concurrent with Git inspection and
        # makes the second manifest the last live source read before the
        # pending-signal snapshot and logical commit.
        source_first = self._live_source_manifest()
        git_current = self._live_git_receipt()
        source_second = self._live_source_manifest()
        validate_source_manifest_receipt(
            source_first, "final source manifest before Git"
        )
        validate_git_receipt(
            git_current, self.expected_commit, "final Git receipt"
        )
        validate_source_manifest_receipt(
            source_second, "final source manifest after Git"
        )
        if not canonical_equal(source_first, source_second):
            fail("source manifest changed during final source/Git sandwich")
        if source_second["sha256"] != self.expected_source_manifest_sha256:
            fail("final source manifest differs from the presealed value")
        for name in ("source_manifest_before", "source_manifest_after"):
            if not canonical_equal(summary[name], source_second):
                fail(f"final source manifest differs from bound {name}")
        for name in ("git_before", "git_after"):
            if not canonical_equal(summary[name], git_current):
                fail(f"final Git receipt differs from bound {name}")

    def _mutable_precommit_pass(self) -> List[str]:
        failures: List[str] = []
        for name in OUTPUT_NAMES:
            try:
                os.fsync(self.file_fds[name])
            except BaseException as exc:
                failures.append(f"precommit fsync {name}: {exception_text(exc)}")
        try:
            if (
                self._read_current(
                    self.file_fds["raw.jsonl"], MAX_STDOUT_BYTES + 1
                ) != self.staged_raw
                or self._read_current(
                    self.file_fds["stderr.txt"], MAX_STDERR_BYTES + 1
                ) != self.staged_stderr
            ):
                fail("sealed evidence bytes changed after readback")
            summary = self._validate_summary_cross(self._read_current(
                self.file_fds["summary.json"], 8 * 1024 * 1024
            ))
            if summary["outcome"] in ("pass", "reject"):
                self._recheck_positive_source_git(summary)
            self._verify_roster()
            os.fsync(self.directory_fd)
            self._verify_identities(0o700, 0o600)
            os.fsync(self.parent_fd)
            self._verify_identities(0o700, 0o600)
        except BaseException as exc:
            failures.append(f"mutable precommit check: {exception_text(exc)}")
        return failures

    def _harden_after_commit(self) -> List[str]:
        failures: List[str] = []
        for name in OUTPUT_NAMES:
            try:
                os.fchmod(self.file_fds[name], 0o400)
                os.fsync(self.file_fds[name])
            except BaseException as exc:
                failures.append(f"postcommit seal {name}: {exception_text(exc)}")
        try:
            if (
                self._read_current(
                    self.file_fds["raw.jsonl"], MAX_STDOUT_BYTES + 1
                ) != self.staged_raw
                or self._read_current(
                    self.file_fds["stderr.txt"], MAX_STDERR_BYTES + 1
                ) != self.staged_stderr
            ):
                fail("postcommit evidence bytes changed")
            self._validate_summary_cross(self._read_current(
                self.file_fds["summary.json"], 8 * 1024 * 1024
            ))
            self._verify_roster()
            self._verify_identities(0o700, 0o400)
            os.fsync(self.directory_fd)
            os.fchmod(self.directory_fd, 0o500)
            os.fsync(self.directory_fd)
            self._verify_identities(0o500, 0o400)
            self._verify_roster()
            os.fsync(self.parent_fd)
            self._verify_identities(0o500, 0o400)
            self._verify_roster()
        except BaseException as exc:
            failures.append(f"postcommit directory seal: {exception_text(exc)}")
        return failures

    def _truncate_summary_for_poison(self) -> None:
        fd = self.file_fds["summary.json"]
        os.ftruncate(fd, 0)
        os.fsync(fd)
        if self._read_current(fd, 0) != b"":
            fail("poisoned summary did not read back empty")

    def _poison_noncommitted_summary(self) -> Tuple[bool, List[str]]:
        """Ensure a stale positive candidate cannot survive failed rewrites."""
        diagnostics: List[str] = []
        try:
            self._verify_identities(0o700, 0o600)
            self._truncate_summary_for_poison()
            self.summary_mode_poisoned = False
            diagnostics.append("summary publication poisoned by empty truncation")
            return True, diagnostics
        except BaseException as exc:
            diagnostics.append(
                "summary poison truncation: " + exception_text(exc)
            )
        try:
            self._verify_identities(0o700, 0o600)
            fd = self.file_fds["summary.json"]
            retained = os.fstat(fd)
            named = os.stat(
                "summary.json", dir_fd=self.directory_fd,
                follow_symlinks=False,
            )
            expected = self.file_identities["summary.json"]
            for info in (retained, named):
                if (
                    not stat.S_ISREG(info.st_mode)
                    or (info.st_dev, info.st_ino, info.st_nlink) != expected
                ):
                    fail("summary poison path/FD identity differs")
            os.unlink("summary.json", dir_fd=self.directory_fd)
            os.fsync(self.directory_fd)
            os.fsync(self.parent_fd)
            if os.fstat(fd).st_nlink != 0:
                fail("unlinked poisoned summary retained a directory link")
            try:
                os.stat(
                    "summary.json", dir_fd=self.directory_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                fail("poisoned summary path still exists after unlink")
            self.summary_unlinked = True
            self.summary_mode_poisoned = False
            diagnostics.append("summary publication poisoned by unlink")
            return True, diagnostics
        except BaseException as exc:
            diagnostics.append("summary poison unlink: " + exception_text(exc))
        # Last-resort fail-closed state: preserve the bytes only behind a mode
        # that can never satisfy the sealed-summary contract.  Cleanup retains
        # mode 000 instead of accidentally promoting stale positive bytes to
        # the normal 0400 publication mode.
        self.summary_mode_poisoned = True
        try:
            fd = self.file_fds["summary.json"]
            os.fchmod(fd, 0o000)
            os.fsync(fd)
            if stat.S_IMODE(os.fstat(fd).st_mode) != 0o000:
                fail("mode-poisoned summary did not retain mode 000")
            diagnostics.append("summary publication poisoned by mode 000")
            return True, diagnostics
        except BaseException as exc:
            diagnostics.append("summary poison mode: " + exception_text(exc))
            return False, diagnostics

    def _verify_cleanup_state(self) -> None:
        self._verify_parent_path()
        parent = os.fstat(self.parent_fd)
        if self.parent_identity != (
            parent.st_dev, parent.st_ino, stat.S_IMODE(parent.st_mode),
        ):
            fail("cleanup output-parent receipt changed")
        retained_dir = os.fstat(self.directory_fd)
        named_dir = os.stat(
            self.output_dir.name, dir_fd=self.parent_fd,
            follow_symlinks=False,
        )
        for info in (retained_dir, named_dir):
            if (
                not stat.S_ISDIR(info.st_mode)
                or (info.st_dev, info.st_ino, info.st_nlink)
                != self.directory_identity
                or stat.S_IMODE(info.st_mode) != 0o500
            ):
                fail("cleanup output-directory receipt differs")
        for name in OUTPUT_NAMES:
            retained = os.fstat(self.file_fds[name])
            expected = self.file_identities[name]
            if name == "summary.json" and self.summary_unlinked:
                if (
                    not stat.S_ISREG(retained.st_mode)
                    or (retained.st_dev, retained.st_ino) != expected[:2]
                    or retained.st_nlink != 0
                    or stat.S_IMODE(retained.st_mode) != 0o400
                ):
                    fail("cleanup unlinked-summary receipt differs")
                try:
                    os.stat(
                        name, dir_fd=self.directory_fd,
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    continue
                fail("cleanup unlinked-summary path reappeared")
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False
            )
            expected_mode = (
                0o000
                if name == "summary.json" and self.summary_mode_poisoned
                else 0o400
            )
            for info in (retained, named):
                if (
                    not stat.S_ISREG(info.st_mode)
                    or (info.st_dev, info.st_ino, info.st_nlink) != expected
                    or stat.S_IMODE(info.st_mode) != expected_mode
                ):
                    fail(f"cleanup output {name} receipt differs")
        expected_roster = (
            ("raw.jsonl", "stderr.txt")
            if self.summary_unlinked else OUTPUT_NAMES
        )
        self._verify_roster(expected_roster)

    def finish(
        self, summary_bytes: bytes,
        signal_guard: Optional[SignalGuard] = None,
    ) -> BundleFinishResult:
        if self.closed:
            failure = "output bundle was already closed"
            return BundleFinishResult(
                logical_commit_succeeded=self.logically_sealed,
                precommit_failures=(
                    tuple() if self.logically_sealed else (failure,)
                ),
                postcommit_failures=(
                    (failure,) if self.logically_sealed else tuple()
                ),
            )
        failures: List[str] = []
        postcommit_failures: List[str] = []
        published_summary: Optional[Dict[str, Any]] = None
        committed = self.logically_sealed
        hardening_completed = False
        if committed:
            postcommit_failures.append(
                "output finalization resumed after logical commit"
            )
        try:
            if not committed:
                try:
                    summary = self._validate_summary_cross(summary_bytes)
                    readback = self._write_exact(
                        self.file_fds["summary.json"], summary_bytes
                    )
                    if self._validate_summary_cross(readback) != summary:
                        fail("summary readback semantic object differs")
                    published_summary = summary
                except BaseException as exc:
                    self._append_unique(
                        failures, f"summary publication: {exception_text(exc)}"
                    )
                if published_summary is None:
                    try:
                        fallback = self.emergency_summary(
                            "summary-publication-failure", failures
                        )
                        published_summary = self._validate_summary_cross(
                            self._write_exact(
                                self.file_fds["summary.json"], fallback
                            )
                        )
                    except BaseException as exc:
                        self._append_unique(
                            failures,
                            f"fallback summary publication: {exception_text(exc)}",
                        )

                if signal_guard is not None:
                    try:
                        signal_guard.block_for_commit()
                    except BaseException as exc:
                        self._append_unique(
                            failures,
                            f"block commit signals: {exception_text(exc)}",
                        )

                # With the guarded set blocked, grow the finite observed-name
                # set, rewrite the same full schema as invalid whenever it grows
                # (or a mutable precommit failure appears), and commit only
                # after one complete pending/readback snapshot adds nothing.
                # Files remain 0600 and the directory 0700 throughout this loop,
                # so no external observer can mistake a rewritable candidate for
                # the final immutable artifact.
                reflected_failures: Tuple[str, ...] = tuple()
                reflected_signals: Tuple[str, ...] = tuple()
                for _ in range(len(SignalGuard.SIGNALS) + 6):
                    if signal_guard is not None and signal_guard.seal_blocked:
                        try:
                            signal_guard.collect_pending()
                        except BaseException as exc:
                            self._append_unique(
                                failures,
                                f"collect commit signals: {exception_text(exc)}",
                            )
                    current_signals = tuple(
                        signal_guard.observed_signals
                        if signal_guard is not None else self.signal_names
                    )
                    if signal_guard is not None:
                        self.signal_name = signal_guard.first_signal
                        self.signal_names = list(signal_guard.observed_signals)
                    need_rewrite = (
                        published_summary is None
                        or tuple(failures) != reflected_failures
                        or current_signals != reflected_signals
                        or (
                            (failures or current_signals)
                            and published_summary.get("outcome") != "invalid"
                        )
                    )
                    if need_rewrite:
                        try:
                            base = published_summary or parse_summary_bytes(
                                self.emergency_summary(
                                    "unavailable-summary-preimage", failures
                                )
                            )
                            replacement = self._invalid_summary(
                                base, failures, signal_guard
                            )
                            published_summary = self._validate_summary_cross(
                                self._write_exact(
                                    self.file_fds["summary.json"], replacement
                                )
                            )
                            reflected_failures = tuple(failures)
                            reflected_signals = current_signals
                        except BaseException as exc:
                            self._append_unique(
                                failures,
                                f"invalid summary rewrite: {exception_text(exc)}",
                            )
                            continue

                    try:
                        precommit_failures = self._mutable_precommit_pass()
                    except BaseException as exc:
                        precommit_failures = [
                            f"mutable precommit call: {exception_text(exc)}"
                        ]
                    new_failure = False
                    for message in precommit_failures:
                        new_failure = (
                            self._append_unique(failures, message) or new_failure
                        )
                    new_signal = False
                    final_signal_check_failed = False
                    if signal_guard is not None and signal_guard.seal_blocked:
                        try:
                            new_signal = signal_guard.collect_pending()
                        except BaseException as exc:
                            final_signal_check_failed = True
                            new_failure = self._append_unique(
                                failures,
                                f"collect final commit signals: {exception_text(exc)}",
                            ) or new_failure
                    if new_failure or new_signal:
                        continue
                    if precommit_failures or final_signal_check_failed:
                        # A persistent mutable-precommit or final signal-snapshot
                        # failure is represented by the invalid summary, but
                        # cannot produce a logical commit.
                        break
                    if signal_guard is not None:
                        try:
                            signal_guard.commit_logical_seal()
                        except BaseException as exc:
                            new_commit_failure = self._append_unique(
                                failures,
                                f"logical signal commit: {exception_text(exc)}",
                            )
                            if new_commit_failure:
                                continue
                            break
                    self.logically_sealed = True
                    committed = True
                    try:
                        hardening_failures = self._harden_after_commit()
                        hardening_completed = True
                        for message in hardening_failures:
                            self._append_unique(postcommit_failures, message)
                    except BaseException as exc:
                        self._append_unique(
                            postcommit_failures,
                            f"postcommit hardening call: {exception_text(exc)}",
                        )
                    break
            if not committed:
                self._append_unique(
                    failures, "logical output commit did not complete"
                )
                # A persistent precommit/checkpoint failure cannot earn a
                # logical commit.  Before any best-effort immutable modes are
                # exposed, nevertheless rewrite and read back one full-schema
                # invalid object containing the complete terminal diagnosis.
                # This path intentionally bypasses an injected/broken
                # `_mutable_precommit_pass`: the prior successful exact-write
                # primitive remains the narrowest available evidence sink.
                for _ in range(3):
                    if signal_guard is not None and signal_guard.seal_blocked:
                        try:
                            signal_guard.collect_pending()
                        except BaseException as exc:
                            self._append_unique(
                                failures,
                                "collect abort signals: "
                                + exception_text(exc),
                            )
                    try:
                        base = published_summary or parse_summary_bytes(
                            self.emergency_summary(
                                "uncommitted-output-publication", failures
                            )
                        )
                        replacement = self._invalid_summary(
                            base, failures, signal_guard
                        )
                        published_summary = self._validate_summary_cross(
                            self._write_exact(
                                self.file_fds["summary.json"], replacement
                            )
                        )
                        break
                    except BaseException as exc:
                        if not self._append_unique(
                            failures,
                            "terminal invalid summary: "
                            + exception_text(exc),
                        ):
                            break
                invalid_summary_secured = bool(
                    published_summary is not None
                    and published_summary.get("outcome") == "invalid"
                )
                if not invalid_summary_secured:
                    # Set the fail-closed mode before entering the whole poison
                    # call so even a BaseException at its first instruction
                    # cannot promote stale positive bytes to mode 0400.
                    self.summary_mode_poisoned = True
                    try:
                        _, poison_diagnostics = (
                            self._poison_noncommitted_summary()
                        )
                    except BaseException as exc:
                        poison_diagnostics = [
                            "summary poison call: " + exception_text(exc)
                        ]
                    for message in poison_diagnostics:
                        self._append_unique(failures, message)
        except BaseException as exc:
            # Catch failures in the commit-loop bookkeeping itself, not only
            # failures raised by the publication primitives it calls.  Before
            # a logical commit, establish the fail-closed mode first so no
            # secondary exception below can promote an earlier pass/reject
            # candidate to the normal 0400 summary contract.
            if committed:
                try:
                    self._append_unique(
                        postcommit_failures,
                        "postcommit finalization: " + exception_text(exc),
                    )
                except BaseException:
                    pass
            else:
                self.summary_mode_poisoned = True
                try:
                    self._append_unique(
                        failures,
                        "uncommitted finalization: " + exception_text(exc),
                    )
                except BaseException:
                    pass
                invalid_summary_secured = False
                try:
                    if signal_guard is not None:
                        if not signal_guard.seal_blocked:
                            signal_guard.block_for_commit()
                        signal_guard.collect_pending()
                        self.signal_name = signal_guard.first_signal
                        self.signal_names = list(
                            signal_guard.observed_signals
                        )
                    base = published_summary or parse_summary_bytes(
                        self.emergency_summary(
                            "uncommitted-finalization-exception", failures
                        )
                    )
                    replacement = self._invalid_summary(
                        base, failures, signal_guard
                    )
                    published_summary = self._validate_summary_cross(
                        self._write_exact(
                            self.file_fds["summary.json"], replacement
                        )
                    )
                    invalid_summary_secured = (
                        published_summary.get("outcome") == "invalid"
                    )
                except BaseException as publish_exc:
                    try:
                        self._append_unique(
                            failures,
                            "exception invalid summary: "
                            + exception_text(publish_exc),
                        )
                    except BaseException:
                        pass
                if invalid_summary_secured:
                    self.summary_mode_poisoned = False
                else:
                    try:
                        _, poison_diagnostics = (
                            self._poison_noncommitted_summary()
                        )
                    except BaseException as poison_exc:
                        poison_diagnostics = [
                            "summary poison call: "
                            + exception_text(poison_exc)
                        ]
                    for message in poison_diagnostics:
                        try:
                            self._append_unique(failures, message)
                        except BaseException:
                            break
        finally:
            # Before a successful logical commit, cleanup diagnostics are part
            # of the invalid outcome.  After commit, close-only diagnostics do
            # not retroactively contradict the immutable decision.
            cleanup_failures: List[str] = []
            close_diagnostics: List[str] = []
            # If the authoritative hardening pass failed after logical commit,
            # make one last best-effort mode/durability sweep.  The original
            # post-commit diagnostic remains authoritative even if this sweep
            # repairs the on-disk modes; it is never folded into or written over
            # the committed summary decision.
            run_cleanup_hardening = (
                not committed
                or not hardening_completed
                or bool(postcommit_failures)
            )
            if run_cleanup_hardening:
                for name, fd in tuple(self.file_fds.items()):
                    try:
                        final_mode = (
                            0o000
                            if name == "summary.json"
                            and self.summary_mode_poisoned
                            else 0o400
                        )
                        os.fchmod(fd, final_mode)
                        os.fsync(fd)
                    except BaseException as exc:
                        cleanup_failures.append(
                            f"finalize {name}: {exception_text(exc)}"
                        )
                if self.directory_fd >= 0:
                    try:
                        os.fchmod(self.directory_fd, 0o500)
                        os.fsync(self.directory_fd)
                        if self.parent_fd >= 0:
                            os.fsync(self.parent_fd)
                        # Verify through the still-open retained file, directory,
                        # and parent descriptors before any close can erase the
                        # evidence needed to distinguish a path substitution.
                        self._verify_cleanup_state()
                    except BaseException as exc:
                        cleanup_failures.append(
                            f"finalize directory: {exception_text(exc)}"
                        )
            for name, fd in tuple(self.file_fds.items()):
                try:
                    self._close_owned_fd(fd)
                except BaseException as exc:
                    close_diagnostics.append(
                        f"close {name}: {exception_text(exc)}"
                    )
            self.file_fds.clear()
            if self.directory_fd >= 0:
                try:
                    self._close_owned_fd(self.directory_fd)
                except BaseException as exc:
                    close_diagnostics.append(
                        f"close directory: {exception_text(exc)}"
                    )
                self.directory_fd = -1
            if self.parent_fd >= 0:
                try:
                    self._close_owned_fd(self.parent_fd)
                except BaseException as exc:
                    close_diagnostics.append(
                        f"close parent: {exception_text(exc)}"
                    )
                self.parent_fd = -1
            if not committed:
                for message in cleanup_failures + close_diagnostics:
                    self._append_unique(failures, message)
            else:
                # Once authoritative hardening and verification completed,
                # close-only noise cannot reclassify the already committed
                # summary or controller status.  Durability/mode/identity
                # cleanup failures remain exit-affecting postcommit failures.
                for message in cleanup_failures:
                    self._append_unique(postcommit_failures, message)
            self.closed = True
        return BundleFinishResult(
            logical_commit_succeeded=committed,
            precommit_failures=tuple(failures),
            postcommit_failures=tuple(postcommit_failures),
        )


def parse_cpu_tick_receipt(
    data: bytes, cpu: int, read_monotonic_ns: int
) -> Dict[str, Any]:
    if type(cpu) is not int or cpu < 0:
        fail("CPU tick receipt requires a nonnegative CPU")
    if type(read_monotonic_ns) is not int or read_monotonic_ns < 0:
        fail("CPU tick receipt requires a nonnegative timestamp")
    label = f"cpu{cpu}"
    matches = []
    for physical in data.splitlines():
        fields = physical.split()
        if fields and fields[0] == label.encode("ascii"):
            matches.append(fields)
    if len(matches) != 1 or len(matches[0]) < 9:
        fail(f"/proc/stat lacks one complete {label} row")
    values: List[int] = []
    for index, raw in enumerate(matches[0][1:9]):
        try:
            text = raw.decode("ascii")
        except UnicodeDecodeError:
            fail(f"{label} tick {index} is not ASCII")
        if not text.isdecimal() or (len(text) > 1 and text.startswith("0")):
            fail(f"{label} tick {index} is not canonical")
        value = int(text)
        if value > MAX_UINT64:
            fail(f"{label} tick {index} exceeds uint64")
        values.append(value)
    # Linux accounts guest time inside user/nice, so exclude the guest fields
    # and define non-idle from the first eight canonical counters only.
    non_idle = sum(values[index] for index in (0, 1, 2, 5, 6, 7))
    receipt = {
        "cpu": cpu,
        "non_idle_ticks": non_idle,
        "read_monotonic_ns": read_monotonic_ns,
        "tick_fields": {
            "idle": values[3], "iowait": values[4], "irq": values[5],
            "nice": values[1], "softirq": values[6], "steal": values[7],
            "system": values[2], "user": values[0],
        },
    }
    validate_cpu_tick_receipt(receipt, label)
    return receipt


def validate_cpu_tick_receipt(
    receipt: Mapping[str, Any], where: str,
) -> None:
    if type(receipt) is not dict:
        fail(f"{where} tick receipt is not an object")
    exact_keys(receipt, CPU_TICK_RECEIPT_KEYS, f"{where} tick receipt")
    exact_int(receipt["cpu"], 0, MAX_INT63, f"{where} tick CPU")
    exact_int(
        receipt["read_monotonic_ns"], 0, MAX_INT63,
        f"{where} tick timestamp",
    )
    fields = receipt["tick_fields"]
    if type(fields) is not dict:
        fail(f"{where} tick fields are not an object")
    exact_keys(fields, CPU_TICK_FIELD_KEYS, f"{where} tick fields")
    for name in CPU_TICK_FIELD_KEYS:
        exact_int(fields[name], 0, MAX_UINT64, f"{where} tick field {name}")
    recomputed = sum(
        fields[name]
        for name in ("user", "nice", "system", "irq", "softirq", "steal")
    )
    exact_int(recomputed, 0, MAX_UINT64, f"{where} non-idle tick sum")
    exact_int(
        receipt["non_idle_ticks"], recomputed, recomputed,
        f"{where} non-idle ticks",
    )


def read_cpu_tick_receipt(cpu: int) -> Dict[str, Any]:
    try:
        data = Path("/proc/stat").read_bytes()
    except OSError as exc:
        fail(f"cannot read /proc/stat: {exc}")
    return parse_cpu_tick_receipt(data, cpu, time.monotonic_ns())


THERMAL_HEADER = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)


def append_unique_failure(values: List[str], message: str) -> None:
    if type(message) is not str or not message:
        message = "unspecified health failure"
    message = message if len(message) <= 4000 else message[:3997] + "..."
    if message not in values and len(values) < MAX_PUBLICATION_FAILURES:
        values.append(message)


def empty_health_receipt(args: argparse.Namespace) -> Dict[str, Any]:
    receipt: Dict[str, Any] = {
        "admission_sibling_ticks": {},
        "child_reap_monotonic_ns": None,
        "child_start_monotonic_ns": None,
        "collection_failures": [],
        "controller_core": [],
        "controller_cpu": EXPECTED_CONTROLLER_CPU,
        "controller_initial_affinity": [],
        "controller_singleton_affinity_end": [],
        "edac_policy": "ce-and-ue-every-sample-zero-v1",
        "evidence_status": "partial",
        "receipt_sha256": None,
        "sampler": {},
        "sampler_admission": {},
        "sampler_admission_receipt_sha256": None,
        "sampler_core": [],
        "sampler_cpu": EXPECTED_SAMPLER_CPU,
        "sampler_receipt_sha256": None,
        "sampler_terminal": {},
        "sampler_terminal_receipt_sha256": None,
        "schema": HEALTH_SCHEMA,
        "sibling_non_idle_tick_cap": SIBLING_NON_IDLE_TICK_CAP,
        "sibling_tick_policy": (
            "linux-proc-stat-user-nice-system-irq-softirq-steal-v1"
        ),
        "sibling_ticks": [],
        "target_core": [],
        "target_cpu": EXPECTED_TARGET_CPU,
        "target_threads": [],
        "thermal": {},
        "thermal_max_millic": THERMAL_MAX_MILLIC,
        "terminal_status": "invalid",
        "violations": [],
    }
    receipt["receipt_sha256"] = sha256_bytes(canonical_bytes(receipt))
    return receipt


def finalize_health_receipt(
    receipt: Dict[str, Any], collection_failures: Sequence[str],
    violations: Sequence[str],
) -> Dict[str, Any]:
    receipt["collection_failures"] = list(dict.fromkeys(collection_failures))
    receipt["violations"] = list(dict.fromkeys(violations))
    complete = (
        not receipt["collection_failures"]
        and receipt["child_start_monotonic_ns"] is not None
        and receipt["child_reap_monotonic_ns"] is not None
        and bool(receipt["admission_sibling_ticks"])
        and len(receipt["sibling_ticks"]) == 1
        and bool(receipt["sampler_admission"])
        and bool(receipt["sampler"])
        and not receipt["sampler_terminal"]
        and bool(receipt["thermal"])
        and bool(receipt["target_core"])
        and bool(receipt["controller_core"])
        and bool(receipt["sampler_core"])
        and bool(receipt["target_threads"])
        and bool(receipt["controller_initial_affinity"])
        and bool(receipt["controller_singleton_affinity_end"])
    )
    receipt["evidence_status"] = "complete" if complete else "partial"
    receipt["terminal_status"] = (
        "ok" if complete and not receipt["violations"] else "invalid"
    )
    receipt["receipt_sha256"] = None
    receipt["receipt_sha256"] = sha256_bytes(canonical_bytes(receipt))
    return receipt


def sampler_snapshot(path: Path) -> Tuple[bytes, os.stat_result]:
    fd = os.open(str(path), nonblocking_read_flags("sampler snapshot"))
    try:
        info = os.fstat(fd)
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_nlink != 1
            or info.st_size < 1
            or info.st_size > MAX_SAMPLER_CSV_BYTES
        ):
            fail("sampler CSV file policy differs")
        data = bytearray()
        offset = 0
        while offset < info.st_size:
            block = os.pread(fd, min(1024 * 1024, info.st_size - offset), offset)
            if not block:
                fail("sampler CSV snapshot is short")
            data.extend(block)
            offset += len(block)
        if not data or b"\n" not in data:
            fail("sampler CSV has no complete physical line")
        # The live sampler may have one incomplete append.  Bind only the
        # complete prefix; the held inode and endpoint sample bind the window.
        complete = bytes(data[: data.rfind(b"\n") + 1])
        return complete, info
    finally:
        os.close(fd)


def stable_regular_snapshot(
    path: Path, max_bytes: int, where: str, *, attempts: int = 4,
) -> Tuple[bytes, os.stat_result]:
    exact_int(max_bytes, 0, MAX_FINAL_ARTIFACT_BYTES, where + " size bound")
    exact_int(attempts, 1, 8, where + " snapshot attempts")
    last_error = "changed during bounded snapshot"
    for _ in range(attempts):
        fd = os.open(str(path), nonblocking_read_flags(where))
        try:
            before = os.fstat(fd)
            if (
                not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
                or before.st_size < 0 or before.st_size > max_bytes
            ):
                fail(f"{where} regular-file policy differs")
            data = bytearray()
            offset = 0
            while offset < before.st_size:
                block = os.pread(
                    fd, min(1024 * 1024, before.st_size - offset), offset
                )
                if not block:
                    fail(f"{where} snapshot is short")
                data.extend(block)
                offset += len(block)
            after = os.fstat(fd)
            named = os.stat(str(path), follow_symlinks=False)
            if same_file_receipt(before, after) and same_file_receipt(after, named):
                return bytes(data), after
            last_error = "changed during bounded snapshot"
        except FileNotFoundError:
            last_error = "disappeared during bounded snapshot"
        finally:
            os.close(fd)
    fail(f"{where} {last_error}")


def sampler_artifact_receipt(
    path: Path, max_bytes: int, expected_mode: int, expected_uid: int,
    expected_gid: int, where: str,
) -> Tuple[Dict[str, Any], bytes, os.stat_result]:
    data, info = stable_regular_snapshot(path, max_bytes, where)
    if (
        stat.S_IMODE(info.st_mode) != expected_mode
        or info.st_uid != expected_uid or info.st_gid != expected_gid
        or info.st_nlink != 1
    ):
        fail(f"{where} mode/owner/link policy differs")
    return {
        "bytes": len(data), "path": str(path),
        "sha256": sha256_bytes(data), "stat": stat_receipt(info),
    }, data, info


SAMPLER_ARTIFACT_KEYS = {"bytes", "path", "sha256", "stat"}


def validate_sampler_artifact_receipt(
    value: Any, expected_path: Path, expected_mode: int,
    expected_uid: int, expected_gid: int, max_bytes: int, where: str,
) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, SAMPLER_ARTIFACT_KEYS, where)
    exact_string(value["path"], str(expected_path), where + ".path")
    exact_int(value["bytes"], 0, max_bytes, where + ".bytes")
    lower_hash(value["sha256"], where + ".sha256")
    validate_stat_receipt(value["stat"], where + ".stat")
    if (
        value["stat"]["size"] != value["bytes"]
        or value["stat"]["mode"] != expected_mode
        or value["stat"]["uid"] != expected_uid
        or value["stat"]["gid"] != expected_gid
        or value["stat"]["nlink"] != 1
    ):
        fail(f"{where} stat/content policy differs")


def parse_canonical_json_line(raw: bytes, where: str) -> Dict[str, Any]:
    if not raw.endswith(b"\n") or b"\r" in raw or raw.count(b"\n") != 1:
        fail(f"{where} is not one LF-terminated line")
    try:
        value = json.loads(
            raw[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        fail(f"{where} is not canonical ASCII JSON: {exc}")
    if type(value) is not dict or canonical_bytes(value) + b"\n" != raw:
        fail(f"{where} canonical JSON encoding differs")
    return value


def expected_validation_header(args: argparse.Namespace) -> Dict[str, Any]:
    return {
        "expected_output_owner_uid": args.expected_sampler_uid,
        "raw_columns": list(THERMAL_HEADER),
        "sampler_source_expected_sha256": args.expected_sampler_script_sha256,
        "sampling": {
            "dimm_attempts": int(EXPECTED_SAMPLER_DIMM_ATTEMPTS_TEXT),
            "dimm_retry_delay_s": float(EXPECTED_SAMPLER_DIMM_RETRY_DELAY_TEXT),
            "interval_s": float(EXPECTED_SAMPLER_INTERVAL_TEXT),
        },
        "schema": THERMAL_VALIDATION_STREAM_SCHEMA,
        "thresholds": dict(THERMAL_SAMPLER_THRESHOLDS),
    }


def validate_validation_stream_header(
    data: bytes, args: argparse.Namespace, where: str,
) -> None:
    if not data or b"\n" not in data:
        fail(f"{where} has no complete header")
    header_raw = data[:data.find(b"\n") + 1]
    header = parse_canonical_json_line(header_raw, where + " header")
    if not canonical_equal(header, expected_validation_header(args)):
        fail(f"{where} header authority differs")


VALIDATION_SAMPLE_KEYS = {
    "consecutive_fault_rows", "decision", "edac_ce_delta", "edac_ue_delta",
    "fault_count", "hot_sensors", "monotonic_s", "read_error_count",
    "sample_index", "schema", "sensors",
}
VALIDATION_SENSOR_KEYS = {
    "attempt_errors", "hot", "hot_streak", "jump_c", "rate_c_per_s",
    "raw_c", "reason", "valid",
}
VALIDATION_DIMM_FIELDS = tuple(THERMAL_HEADER[5:13])


def parse_validation_sample(value: Any, where: str) -> Dict[str, Any]:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, VALIDATION_SAMPLE_KEYS, where)
    exact_string(value["schema"], THERMAL_VALIDATION_SAMPLE_SCHEMA, where + ".schema")
    for name in (
        "consecutive_fault_rows", "edac_ce_delta", "edac_ue_delta",
        "fault_count", "read_error_count", "sample_index",
    ):
        exact_int(value[name], 0, MAX_UINT64, where + "." + name)
    if value["decision"] not in (
        "continue", "thermal_abort", "telemetry_abort", "edac_abort",
    ):
        fail(f"{where} decision differs")
    monotonic_s = value["monotonic_s"]
    if type(monotonic_s) not in (int, float) or not math.isfinite(monotonic_s) or monotonic_s < 0:
        fail(f"{where} monotonic timestamp differs")
    hot_sensors = value["hot_sensors"]
    if (
        type(hot_sensors) is not list
        or any(type(item) is not str or item not in VALIDATION_DIMM_FIELDS
               for item in hot_sensors)
        or hot_sensors != sorted(set(hot_sensors))
    ):
        fail(f"{where} hot sensor roster differs")
    sensors = value["sensors"]
    if type(sensors) is not dict or set(sensors) != set(VALIDATION_DIMM_FIELDS):
        fail(f"{where} sensor roster differs")
    for name in VALIDATION_DIMM_FIELDS:
        sensor = sensors[name]
        sensor_where = where + ".sensors." + name
        if type(sensor) is not dict:
            fail(f"{sensor_where} is not an object")
        exact_keys(sensor, VALIDATION_SENSOR_KEYS, sensor_where)
        exact_int(sensor["attempt_errors"], 0, MAX_UINT64, sensor_where + ".attempt_errors")
        exact_int(sensor["hot_streak"], 0, MAX_UINT64, sensor_where + ".hot_streak")
        if type(sensor["hot"]) is not bool or type(sensor["valid"]) is not bool:
            fail(f"{sensor_where} boolean fields differ")
        if type(sensor["reason"]) is not str or not sensor["reason"] or len(sensor["reason"]) > 128:
            fail(f"{sensor_where} reason differs")
        for scalar_name in ("jump_c", "rate_c_per_s", "raw_c"):
            scalar = sensor[scalar_name]
            if scalar is not None and (
                type(scalar) not in (int, float) or not math.isfinite(scalar)
            ):
                fail(f"{sensor_where}.{scalar_name} differs")
    return value


def validate_live_sampler_stream_pair(
    raw_csv: bytes, validation_jsonl: bytes,
    args: argparse.Namespace, where: str,
) -> bool:
    """Validate one live raw/validation snapshot; report whether it settled.

    The sampler deliberately flushes each raw row before its matching
    validation record.  A raw-ahead snapshot is therefore retryable, but a
    complete terminal validation decision is never retryable or ignorable.
    """
    raw_lines = raw_csv.splitlines(keepends=True)
    if (
        not raw_csv.endswith(b"\n") or len(raw_lines) < 1
        or tuple(parse_csv_physical_line(raw_lines[0], where + " raw header"))
        != THERMAL_HEADER
    ):
        fail(f"{where} raw stream framing differs")
    complete_end = validation_jsonl.rfind(b"\n") + 1
    validation_complete = validation_jsonl[:complete_end]
    validation_partial = validation_jsonl[complete_end:]
    validate_validation_stream_header(
        validation_complete, args, where + " validation stream"
    )
    validation_lines = validation_complete.splitlines(keepends=True)
    previous_monotonic_ns: Optional[int] = None
    previous_csv_timestamp_ns: Optional[int] = None
    for index, raw in enumerate(validation_lines[1:]):
        sample = parse_validation_sample(
            parse_canonical_json_line(
                raw, f"{where} validation record {index}"
            ),
            f"{where} validation record {index}",
        )
        if sample["sample_index"] != index:
            fail(f"{where} validation sample indices differ")
        monotonic_ns = exact_int(
            int(round(float(sample["monotonic_s"]) * 1_000_000_000.0)),
            0, MAX_INT63, f"{where} validation timestamp",
        )
        csv_timestamp_ns = exact_int(
            int(round(
                float(f"{float(sample['monotonic_s']):.6f}")
                * 1_000_000_000.0
            )),
            0, MAX_INT63, f"{where} validation CSV timestamp",
        )
        if (
            previous_monotonic_ns is not None
            and monotonic_ns <= previous_monotonic_ns
            or previous_csv_timestamp_ns is not None
            and csv_timestamp_ns <= previous_csv_timestamp_ns
        ):
            fail(f"{where} validation timestamps are not strictly increasing")
        previous_monotonic_ns = monotonic_ns
        previous_csv_timestamp_ns = csv_timestamp_ns
        if sample["decision"] != "continue":
            fail(
                f"{where} live validation contains terminal decision "
                + sample["decision"]
            )
    if validation_partial:
        return False
    raw_sample_count = len(raw_lines) - 1
    validation_sample_count = len(validation_lines) - 1
    if validation_sample_count > raw_sample_count:
        fail(f"{where} validation stream advanced ahead of raw evidence")
    if raw_sample_count != validation_sample_count:
        return False
    expected_header = canonical_bytes(expected_validation_header(args)) + b"\n"
    replay_terminal_validation(
        raw_csv, validation_jsonl,
        {"validation_header_ascii": expected_header.decode("ascii")},
    )
    return True


def validate_stale_sampler_stream_prefix(
    raw_csv: bytes, validation_jsonl: bytes,
    args: argparse.Namespace, where: str, observed_monotonic_ns: int,
) -> Dict[str, Any]:
    """Validate the maximal producer-authentic paired prefix of a stall."""
    exact_int(
        observed_monotonic_ns, 0, MAX_INT63,
        where + " snapshot observation",
    )
    try:
        raw_csv.decode("ascii")
        validation_jsonl.decode("ascii")
    except UnicodeDecodeError:
        fail(f"{where} stream suffix is not ASCII")
    stream_parts = split_terminal_streams(
        raw_csv, validation_jsonl, allow_unpaired=True
    )
    suffix_shape = stream_parts["suffix_shape"]
    if (
        suffix_shape[2]
        and not stream_parts["validation_partial"].startswith(b"{")
        or suffix_shape[1]
        and (
            not stream_parts["raw_suffix"]
            or b"\r" in stream_parts["raw_suffix"]
            or b"\n" in stream_parts["raw_suffix"]
        )
    ):
        fail(f"{where} partial stream suffix framing differs")
    if not validate_live_sampler_stream_pair(
        stream_parts["paired_raw"], stream_parts["validation_complete"],
        args, where + " paired prefix",
    ):
        fail(f"{where} paired prefix did not settle")
    expected_header = canonical_bytes(expected_validation_header(args)) + b"\n"
    _, final_decision = replay_terminal_validation(
        stream_parts["paired_raw"], stream_parts["validation_complete"],
        {"validation_header_ascii": expected_header.decode("ascii")},
    )
    validation_lines = stream_parts["validation_complete"].splitlines(
        keepends=True
    )
    last_monotonic_ns: Optional[int] = None
    last_csv_timestamp_ns: Optional[int] = None
    for index, raw in enumerate(validation_lines[1:]):
        sample = parse_validation_sample(
            parse_canonical_json_line(
                raw, f"{where} validation record {index}"
            ),
            f"{where} validation record {index}",
        )
        if sample["decision"] != "continue":
            fail(
                f"{where} live validation contains terminal decision "
                + sample["decision"]
            )
        last_monotonic_ns = exact_int(
            int(round(float(sample["monotonic_s"]) * 1_000_000_000.0)),
            0, MAX_INT63, f"{where} validation timestamp",
        )
        last_csv_timestamp_ns = exact_int(
            int(round(
                float(f"{float(sample['monotonic_s']):.6f}")
                * 1_000_000_000.0
            )),
            0, MAX_INT63, f"{where} validation CSV timestamp",
        )
    if (
        final_decision != "continue"
        or last_monotonic_ns is None
        or last_csv_timestamp_ns is None
    ):
        fail(f"{where} paired prefix lacks a continuing heartbeat")
    if stream_parts["suffix_shape"][0] == 1:
        unpaired_row = parse_csv_physical_line(
            stream_parts["raw_suffix"], where + " unpaired raw row"
        )
        if len(unpaired_row) != len(THERMAL_HEADER):
            fail(f"{where} unpaired raw row width differs")
        unpaired_value = finite_scalar(
            unpaired_row[1], where + " unpaired raw timestamp"
        )
        if f"{unpaired_value:.6f}" != unpaired_row[1]:
            fail(f"{where} unpaired raw timestamp is not canonical")
        unpaired_ns = exact_int(
            int(round(unpaired_value * 1_000_000_000.0)),
            0, MAX_INT63, where + " unpaired raw timestamp",
        )
        if not (
            last_csv_timestamp_ns < unpaired_ns <= observed_monotonic_ns
        ):
            fail(f"{where} unpaired raw chronology differs")
    stream_parts["last_monotonic_ns"] = last_monotonic_ns
    stream_parts["last_csv_timestamp_ns"] = last_csv_timestamp_ns
    return stream_parts


def settle_live_sampler_streams(
    args: argparse.Namespace, health_api: Any,
    deadline: Optional[float], where: str, *, require_fresh: bool = True,
    allow_stale_suffixes: bool = False,
) -> Dict[str, Any]:
    """Return a stable, fully paired, all-continue live sampler snapshot."""
    settle_deadline = min(
        deadline if deadline is not None else time.monotonic() + 1.0,
        time.monotonic() + 1.0,
    )

    def read_once() -> Dict[str, Any]:
        pid_file, pid_data, pid_info = sampler_artifact_receipt(
            args.sampler_pid_file, 64, 0o444,
            args.expected_sampler_uid, args.expected_sampler_gid,
            where + " PID file",
        )
        receipt_file, receipt_data, receipt_info = sampler_artifact_receipt(
            args.sampler_receipt, MAX_SAMPLER_RECEIPT_BYTES, 0o444,
            args.expected_sampler_uid, args.expected_sampler_gid,
            where + " terminal receipt reservation",
        )
        validation_jsonl, validation_data, validation_info = (
            sampler_artifact_receipt(
                args.sampler_validation_jsonl,
                MAX_SAMPLER_VALIDATION_BYTES, 0o444,
                args.expected_sampler_uid, args.expected_sampler_gid,
                where + " validation JSONL",
            )
        )
        # Read validation before raw.  Since the producer flushes raw first,
        # a raw-ahead result can be an honest in-flight sample; validation
        # ahead of this later raw snapshot cannot be producer-authentic.
        csv_full_data, csv_info = stable_regular_snapshot(
            args.sampler_csv, MAX_SAMPLER_CSV_BYTES, where + " raw CSV"
        )
        complete_end = csv_full_data.rfind(b"\n") + 1
        if complete_end == 0:
            fail(f"{where} raw CSV has no complete physical line")
        csv_data = csv_full_data[:complete_end]
        if pid_data != (str(args.sampler_pid) + "\n").encode("ascii"):
            fail(f"{where} PID file content differs")
        if receipt_data:
            fail(f"{where} terminal intent is already published")
        for info, expected_device, expected_inode, label in (
            (
                pid_info, args.expected_sampler_pid_file_device,
                args.expected_sampler_pid_file_inode, "PID file",
            ),
            (
                validation_info, args.expected_sampler_validation_device,
                args.expected_sampler_validation_inode, "validation JSONL",
            ),
            (
                receipt_info, args.expected_sampler_receipt_device,
                args.expected_sampler_receipt_inode, "terminal receipt",
            ),
            (
                csv_info, args.expected_sampler_csv_device,
                args.expected_sampler_csv_inode, "raw CSV",
            ),
        ):
            if info.st_dev != expected_device or info.st_ino != expected_inode:
                fail(f"{where} {label} identity differs")
        if (
            stat.S_IMODE(csv_info.st_mode) != 0o600
            or csv_info.st_uid != args.expected_sampler_uid
            or csv_info.st_gid != args.expected_sampler_gid
            or csv_info.st_nlink != 1
        ):
            fail(f"{where} raw CSV mode/owner/link policy differs")
        return {
            "csv_data": csv_data,
            "csv_full_data": csv_full_data,
            "csv_info": csv_info,
            "pid_file": pid_file,
            "receipt_file": receipt_file,
            "validation_data": validation_data,
            "validation_info": validation_info,
            "validation_jsonl": validation_jsonl,
            "snapshot_observed_monotonic_ns": time.monotonic_ns(),
        }

    def classify(snapshot: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
        paired = validate_live_sampler_stream_pair(
            snapshot["csv_data"], snapshot["validation_data"], args, where
        )
        raw_complete = (
            len(snapshot["csv_data"]) == len(snapshot["csv_full_data"])
        )
        if paired and raw_complete:
            return validate_stale_sampler_stream_prefix(
                snapshot["csv_full_data"], snapshot["validation_data"],
                args, where, snapshot["snapshot_observed_monotonic_ns"],
            )
        if not allow_stale_suffixes:
            return None
        return validate_stale_sampler_stream_prefix(
            snapshot["csv_full_data"], snapshot["validation_data"],
            args, where, snapshot["snapshot_observed_monotonic_ns"],
        )

    while True:
        if time.monotonic() >= settle_deadline:
            fail(f"{where} did not reach a paired live-stream fixed point")
        health_api._verify_live_sampler_process(
            args.sampler_pid, EXPECTED_SAMPLER_CPU,
            args.expected_sampler_process_start_ticks,
            args.sampler_script, args.sampler_csv,
        )
        if (
            process_start_ticks(args.sampler_pid)
            != args.expected_sampler_process_start_ticks
        ):
            fail(f"{where} process start identity differs")
        first = read_once()
        first_shape = classify(first)
        if first_shape is None:
            time.sleep(min(
                0.05, max(0.0, settle_deadline - time.monotonic())
            ))
            continue
        security = read_process_security(args.sampler_pid, where + " security")
        validate_process_security(
            security, args.expected_sampler_uid, args.expected_sampler_gid,
            [args.expected_sampler_i2c_gid], where + " security",
        )
        if sorted(os.sched_getaffinity(args.sampler_pid)) != [EXPECTED_SAMPLER_CPU]:
            fail(f"{where} process affinity differs")
        second = read_once()
        second_shape = classify(second)
        if (
            second_shape is None
            or first["csv_full_data"] != second["csv_full_data"]
            or stat_receipt(first["csv_info"])
            != stat_receipt(second["csv_info"])
            or not canonical_equal(first["pid_file"], second["pid_file"])
            or not canonical_equal(
                first["validation_jsonl"], second["validation_jsonl"]
            )
            or not canonical_equal(
                first["receipt_file"], second["receipt_file"]
            )
        ):
            time.sleep(min(
                0.05, max(0.0, settle_deadline - time.monotonic())
            ))
            continue
        health_api._verify_live_sampler_process(
            args.sampler_pid, EXPECTED_SAMPLER_CPU,
            args.expected_sampler_process_start_ticks,
            args.sampler_script, args.sampler_csv,
        )
        if (
            process_start_ticks(args.sampler_pid)
            != args.expected_sampler_process_start_ticks
        ):
            fail(f"{where} process identity changed at the fixed point")
        third_validation, third_validation_data, _ = (
            sampler_artifact_receipt(
                args.sampler_validation_jsonl,
                MAX_SAMPLER_VALIDATION_BYTES, 0o444,
                args.expected_sampler_uid, args.expected_sampler_gid,
                where + " confirming validation JSONL",
            )
        )
        if (
            third_validation_data != second["validation_data"]
            or not canonical_equal(
                third_validation, second["validation_jsonl"]
            )
        ):
            time.sleep(min(
                0.05, max(0.0, settle_deadline - time.monotonic())
            ))
            continue
        health_api._verify_live_sampler_process(
            args.sampler_pid, EXPECTED_SAMPLER_CPU,
            args.expected_sampler_process_start_ticks,
            args.sampler_script, args.sampler_csv,
        )
        if (
            process_start_ticks(args.sampler_pid)
            != args.expected_sampler_process_start_ticks
        ):
            fail(f"{where} process identity changed after confirmation")
        final_receipt, final_receipt_data, final_receipt_info = (
            sampler_artifact_receipt(
                args.sampler_receipt, MAX_SAMPLER_RECEIPT_BYTES, 0o444,
                args.expected_sampler_uid, args.expected_sampler_gid,
                where + " final terminal receipt reservation",
            )
        )
        if (
            final_receipt_data
            or final_receipt_info.st_dev
            != args.expected_sampler_receipt_device
            or final_receipt_info.st_ino
            != args.expected_sampler_receipt_inode
            or not canonical_equal(second["receipt_file"], final_receipt)
        ):
            fail(f"{where} terminal intent appeared at the fixed point")
        last_monotonic_ns = second_shape["last_monotonic_ns"]
        if require_fresh and (
            time.monotonic_ns()
            > last_monotonic_ns + SAMPLER_HEARTBEAT_MAX_GAP_NS
        ):
            fail(f"{where} validation heartbeat is stale")
        second["paired_end_monotonic_ns"] = second_shape[
            "last_csv_timestamp_ns"
        ]
        second["paired_raw"] = second_shape["paired_raw"]
        second["validation_complete"] = second_shape["validation_complete"]
        return second


def select_validation_window(
    validation_data: bytes, raw_window: bytes, where: str,
) -> bytes:
    validation_lines = validation_data.splitlines(keepends=True)
    raw_lines = raw_window.splitlines(keepends=True)
    if len(validation_lines) < 2 or len(raw_lines) < 2:
        fail(f"{where} has insufficient records")
    header = parse_canonical_json_line(validation_lines[0], where + " header")
    if header.get("schema") != THERMAL_VALIDATION_STREAM_SCHEMA:
        fail(f"{where} header schema differs")
    by_timestamp: Dict[str, bytes] = {}
    previous_index: Optional[int] = None
    for index, raw in enumerate(validation_lines[1:], start=1):
        sample = parse_validation_sample(
            parse_canonical_json_line(raw, f"{where} record {index}"),
            f"{where} record {index}",
        )
        if previous_index is not None and sample["sample_index"] != previous_index + 1:
            fail(f"{where} sample indices are not contiguous")
        previous_index = sample["sample_index"]
        timestamp_text = f"{float(sample['monotonic_s']):.6f}"
        if timestamp_text in by_timestamp:
            fail(f"{where} has duplicate formatted timestamps")
        by_timestamp[timestamp_text] = raw
    selected = [validation_lines[0]]
    for index, raw in enumerate(raw_lines[1:], start=2):
        row = parse_csv_physical_line(raw, f"{where} raw row {index}")
        if len(row) != len(THERMAL_HEADER) or row[1] not in by_timestamp:
            fail(f"{where} lacks a validation record for raw row {index}")
        selected.append(by_timestamp[row[1]])
    return b"".join(selected)


def summarize_validation_window_bytes(
    validation_window: bytes, raw_window: bytes,
    sampler: Mapping[str, Any],
) -> Dict[str, Any]:
    validation_lines = validation_window.splitlines(keepends=True)
    raw_lines = raw_window.splitlines(keepends=True)
    if len(validation_lines) != len(raw_lines) or len(validation_lines) < 2:
        fail("thermal validation/raw window row count differs")
    header = parse_canonical_json_line(
        validation_lines[0], "thermal validation header"
    )
    expected_header = {
        "expected_output_owner_uid": sampler["process_uid"],
        "raw_columns": list(THERMAL_HEADER),
        "sampler_source_expected_sha256": sampler["script_sha256"],
        "sampling": {
            "dimm_attempts": 5, "dimm_retry_delay_s": 0.01,
            "interval_s": 1.0,
        },
        "schema": THERMAL_VALIDATION_STREAM_SCHEMA,
        "thresholds": THERMAL_SAMPLER_THRESHOLDS,
    }
    if not canonical_equal(header, expected_header):
        fail("thermal validation header differs from sampler authority")
    failures: List[str] = []
    previous_sample_index: Optional[int] = None
    previous_valid: Dict[str, Tuple[float, float]] = {}
    attempt_errors_total = 0
    for offset, (validation_raw, csv_raw) in enumerate(
        zip(validation_lines[1:], raw_lines[1:]), start=0
    ):
        where = f"thermal validation record {offset}"
        sample = parse_validation_sample(
            parse_canonical_json_line(validation_raw, where), where
        )
        if previous_sample_index is not None and sample["sample_index"] != previous_sample_index + 1:
            fail("thermal validation sample indices are not contiguous")
        previous_sample_index = sample["sample_index"]
        row = parse_csv_physical_line(csv_raw, f"thermal validation raw row {offset}")
        if len(row) != len(THERMAL_HEADER):
            fail("thermal validation raw row width differs")
        sample_time = float(sample["monotonic_s"])
        if f"{sample_time:.6f}" != row[1]:
            fail("thermal validation/raw timestamp differs")
        if canonical_counter(row[13], where + " raw DIMM errors") != sample["read_error_count"]:
            fail("thermal validation/raw DIMM error count differs")
        if canonical_counter(row[17], where + " raw EDAC CE") < sample["edac_ce_delta"]:
            fail("thermal validation EDAC CE delta exceeds raw count")
        if canonical_counter(row[18], where + " raw EDAC UE") < sample["edac_ue_delta"]:
            fail("thermal validation EDAC UE delta exceeds raw count")
        computed_faults = 0
        computed_reads = 0
        computed_hot: List[str] = []
        for sensor_index, name in enumerate(VALIDATION_DIMM_FIELDS, start=5):
            sensor = sample["sensors"][name]
            attempt_errors_total += sensor["attempt_errors"]
            raw_text = row[sensor_index]
            raw_value = None if raw_text == "" else finite_scalar(
                raw_text, where + " raw " + name
            )
            if not canonical_equal(sensor["raw_c"], raw_value):
                fail("thermal validation/raw DIMM value differs")
            if raw_value is None:
                computed_faults += 1
                computed_reads += 1
                if (
                    sensor["valid"] or sensor["reason"] != "read_error"
                    or sensor["hot"] or sensor["hot_streak"] != 0
                    or sensor["jump_c"] is not None or sensor["rate_c_per_s"] is not None
                ):
                    fail("thermal validation read-error state differs")
                continue
            hot = raw_value >= THERMAL_SAMPLER_THRESHOLDS["dimm_safety_c_inclusive"]
            if sensor["hot"] != hot:
                fail("thermal validation hot classification differs")
            if hot and sensor["hot_streak"] >= THERMAL_SAMPLER_THRESHOLDS["hot_confirm_samples"]:
                computed_hot.append(name)
            if sensor["valid"]:
                if sensor["reason"] != "ok":
                    fail("thermal validation valid sensor reason differs")
                if name in previous_valid:
                    previous_value, previous_time = previous_valid[name]
                    elapsed = sample_time - previous_time
                    if elapsed <= 0:
                        fail("thermal validation sensor time is not increasing")
                    jump = abs(raw_value - previous_value)
                    rate = jump / elapsed
                    if (
                        not canonical_equal(sensor["jump_c"], jump)
                        or not canonical_equal(sensor["rate_c_per_s"], rate)
                    ):
                        fail("thermal validation jump/rate replay differs")
                previous_valid[name] = (raw_value, sample_time)
            else:
                computed_faults += 1
        if sample["fault_count"] != computed_faults or sample["read_error_count"] != computed_reads:
            fail("thermal validation fault/read accounting differs")
        if sample["hot_sensors"] != computed_hot:
            fail("thermal validation confirmed-hot roster differs")
        if (
            sample["decision"] != "continue" or sample["fault_count"] != 0
            or sample["read_error_count"] != 0 or sample["hot_sensors"]
            or sample["edac_ce_delta"] != 0 or sample["edac_ue_delta"] != 0
            or any(
                not sensor["valid"] or sensor["reason"] != "ok"
                for sensor in sample["sensors"].values()
            )
        ):
            append_unique_failure(
                failures, f"thermal validation record {offset} is not clean"
            )
    return {
        "attempt_errors_total": attempt_errors_total,
        "failures": failures,
        "sample_count": len(validation_lines) - 1,
    }


SAMPLER_BINDING_KEYS = {
    "device", "gid", "inode", "mode", "nlink", "sha256", "size", "uid",
}
SAMPLER_DESTINATION_KEYS = {
    "basename", "expected_owner_uid", "parent", "path",
}
SAMPLER_DESTINATION_PARENT_KEYS = {"binding", "path"}
SAMPLER_PARENT_BINDING_KEYS = {
    "device", "gid", "inode", "mode", "nlink", "uid",
}
SAMPLER_TERMINAL_KEYS = {
    "admission_receipt_sha256", "coverage", "pid_file_held_binding",
    "process_exit_observation", "process_exit_observed_monotonic_ns",
    "producer_receipt_ascii", "producer_receipt_binding", "raw_csv_ascii",
    "raw_csv_binding", "schema", "terminal_status", "validation_jsonl_ascii",
    "validation_jsonl_binding", "stream_suffixes", "window_csv_ascii",
    "window_validation_jsonl_ascii",
}
SAMPLER_TERMINAL_STREAM_SUFFIX_KEYS = {
    "paired_sample_count", "raw_unpaired_suffix_bytes",
    "raw_unpaired_suffix_sha256", "validation_complete_prefix_bytes",
    "validation_complete_prefix_sha256", "validation_partial_suffix_bytes",
    "validation_partial_suffix_sha256",
}
SAMPLER_TERMINAL_COVERAGE_KEYS = {
    "child_reap_monotonic_ns", "child_start_monotonic_ns",
    "coverage_shortfall_ns", "covers_child_interval",
    "window_end_monotonic_ns", "window_start_monotonic_ns",
}
SAMPLER_PRODUCER_RECEIPT_KEYS = {
    "argv", "cpu_tctl_max_c", "edac_ce_paths", "edac_ue_paths", "errors",
    "exit_code", "finished_monotonic_ns", "finished_utc",
    "expected_output_owner_uid", "outcome", "pid", "pid_file", "raw_csv",
    "raw_samples_preserved", "receipt_file", "sampler_source", "sampling",
    "schema", "self_sha256_excluding_field", "signal",
    "started_monotonic_ns", "started_utc", "summary", "thresholds",
    "validation_jsonl",
}


def validate_sampler_binding(value: Any, where: str) -> Dict[str, Any]:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, SAMPLER_BINDING_KEYS, where)
    for name, lower in (
        ("device", 0), ("gid", 0), ("inode", 1), ("mode", 0),
        ("nlink", 0), ("size", 0), ("uid", 0),
    ):
        exact_int(value[name], lower, MAX_UINT64, where + "." + name)
    lower_hash(value["sha256"], where + ".sha256")
    return value


def binding_from_snapshot(data: bytes, info: os.stat_result) -> Dict[str, Any]:
    return {
        "device": info.st_dev, "gid": info.st_gid, "inode": info.st_ino,
        "mode": stat.S_IMODE(info.st_mode), "nlink": info.st_nlink,
        "sha256": sha256_bytes(data), "size": len(data), "uid": info.st_uid,
    }


def validate_sampler_destination(
    value: Any, expected_path: str, expected_uid: Optional[int],
    expected_parent: Optional[Mapping[str, Any]], where: str,
) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, SAMPLER_DESTINATION_KEYS, where)
    exact_string(value["path"], expected_path, where + ".path")
    exact_string(value["basename"], Path(expected_path).name, where + ".basename")
    if value["expected_owner_uid"] != expected_uid:
        fail(f"{where} expected owner differs")
    parent = value["parent"]
    if type(parent) is not dict:
        fail(f"{where}.parent is not an object")
    exact_keys(parent, SAMPLER_DESTINATION_PARENT_KEYS, where + ".parent")
    exact_string(parent["path"], str(Path(expected_path).parent), where + ".parent.path")
    binding = parent["binding"]
    if type(binding) is not dict:
        fail(f"{where}.parent.binding is not an object")
    exact_keys(binding, SAMPLER_PARENT_BINDING_KEYS, where + ".parent.binding")
    for name in SAMPLER_PARENT_BINDING_KEYS:
        exact_int(binding[name], 0, MAX_UINT64, where + ".parent.binding." + name)
    if expected_parent is not None:
        stat_value = expected_parent["stat"]
        expected_binding = {
            name: stat_value[name] for name in SAMPLER_PARENT_BINDING_KEYS
        }
        if not canonical_equal(binding, expected_binding):
            fail(f"{where} parent binding differs from admission")


def split_terminal_streams(
    raw_csv: bytes, validation_jsonl: bytes, *, allow_unpaired: bool,
) -> Dict[str, Any]:
    raw_complete_end = raw_csv.rfind(b"\n") + 1
    raw_complete = raw_csv[:raw_complete_end]
    raw_partial = raw_csv[raw_complete_end:]
    raw_lines = raw_complete.splitlines(keepends=True)
    if (
        len(raw_lines) < 2
        or tuple(parse_csv_physical_line(raw_lines[0], "terminal raw header"))
        != THERMAL_HEADER
    ):
        fail("terminal raw stream framing differs")
    complete_end = validation_jsonl.rfind(b"\n") + 1
    validation_complete = validation_jsonl[:complete_end]
    validation_partial = validation_jsonl[complete_end:]
    validation_lines = validation_complete.splitlines(keepends=True)
    if not validation_lines:
        fail("terminal validation stream lacks its complete header")
    paired_count = len(validation_lines) - 1
    raw_count = len(raw_lines) - 1
    raw_delta = raw_count - paired_count
    suffix_shape = (
        raw_delta, bool(raw_partial), bool(validation_partial)
    )
    allowed_unpaired_shapes = {
        (0, False, False),
        (1, False, False),
        (1, False, True),
        (0, True, False),
    }
    if (
        not allow_unpaired and suffix_shape != (0, False, False)
        or allow_unpaired and suffix_shape not in allowed_unpaired_shapes
    ):
        fail("terminal raw/validation paired-row shape differs")
    paired_raw = b"".join(raw_lines[:paired_count + 1])
    raw_suffix = raw_csv[len(paired_raw):]
    return {
        "paired_raw": paired_raw,
        "raw_complete": raw_complete,
        "raw_suffix": raw_suffix,
        "suffix_shape": suffix_shape,
        "suffix_receipt": {
            "paired_sample_count": paired_count,
            "raw_unpaired_suffix_bytes": len(raw_suffix),
            "raw_unpaired_suffix_sha256": sha256_bytes(raw_suffix),
            "validation_complete_prefix_bytes": len(validation_complete),
            "validation_complete_prefix_sha256": sha256_bytes(
                validation_complete
            ),
            "validation_partial_suffix_bytes": len(validation_partial),
            "validation_partial_suffix_sha256": sha256_bytes(
                validation_partial
            ),
        },
        "validation_complete": validation_complete,
        "validation_partial": validation_partial,
    }


def replay_terminal_validation(
    raw_csv: bytes, validation_jsonl: bytes,
    sampler_admission: Mapping[str, Any],
) -> Tuple[Dict[str, Any], str]:
    raw_lines = raw_csv.splitlines(keepends=True)
    validation_lines = validation_jsonl.splitlines(keepends=True)
    if (
        len(raw_lines) < 2 or len(validation_lines) != len(raw_lines)
        or tuple(parse_csv_physical_line(raw_lines[0], "terminal raw header"))
        != THERMAL_HEADER
    ):
        fail("terminal raw/validation row roster differs")
    header = parse_canonical_json_line(
        validation_lines[0], "terminal validation header"
    )
    expected_header_raw = sampler_admission["validation_header_ascii"].encode("ascii")
    if canonical_bytes(header) + b"\n" != expected_header_raw:
        fail("terminal validation header differs from admission")
    sensor_totals: Dict[str, Dict[str, Any]] = {
        name: {
            "attempt_errors": 0, "invalid_samples": 0, "max_hot_streak": 0,
            "max_raw_c": None, "max_valid_c": None, "raw_samples": 0,
            "read_error_samples": 0, "valid_samples": 0,
        }
        for name in VALIDATION_DIMM_FIELDS
    }
    ce_baseline: Optional[int] = None
    ue_baseline: Optional[int] = None
    last_ce = 0
    last_ue = 0
    max_consecutive = 0
    cpu_values: List[float] = []
    last_sample: Optional[Dict[str, Any]] = None
    for index, (raw_line, validation_line) in enumerate(
        zip(raw_lines[1:], validation_lines[1:])
    ):
        row = parse_csv_physical_line(raw_line, f"terminal raw row {index}")
        if len(row) != len(THERMAL_HEADER):
            fail("terminal raw row width differs")
        sample = parse_validation_sample(
            parse_canonical_json_line(
                validation_line, f"terminal validation record {index}"
            ),
            f"terminal validation record {index}",
        )
        if sample["sample_index"] != index:
            fail("terminal validation sample indices differ")
        if f"{float(sample['monotonic_s']):.6f}" != row[1]:
            fail("terminal validation/raw timestamp differs")
        cpu_values.append(finite_scalar(row[4], f"terminal CPU row {index}"))
        raw_ce = canonical_counter(row[17], f"terminal EDAC CE row {index}")
        raw_ue = canonical_counter(row[18], f"terminal EDAC UE row {index}")
        current_ce_baseline = raw_ce - sample["edac_ce_delta"]
        current_ue_baseline = raw_ue - sample["edac_ue_delta"]
        if current_ce_baseline < 0 or current_ue_baseline < 0:
            fail("terminal validation EDAC baseline is negative")
        if ce_baseline is None:
            ce_baseline = current_ce_baseline
            ue_baseline = current_ue_baseline
        elif ce_baseline != current_ce_baseline or ue_baseline != current_ue_baseline:
            fail("terminal validation EDAC baseline changed")
        last_ce = raw_ce
        last_ue = raw_ue
        raw_read_errors = canonical_counter(
            row[13], f"terminal DIMM errors row {index}"
        )
        if raw_read_errors != sample["read_error_count"]:
            fail("terminal raw/validation read-error count differs")
        for sensor_offset, name in enumerate(VALIDATION_DIMM_FIELDS, start=5):
            sensor = sample["sensors"][name]
            totals = sensor_totals[name]
            raw_value = None if row[sensor_offset] == "" else finite_scalar(
                row[sensor_offset], f"terminal {name} row {index}"
            )
            if not canonical_equal(raw_value, sensor["raw_c"]):
                fail("terminal raw/validation DIMM value differs")
            totals["attempt_errors"] += sensor["attempt_errors"]
            totals["max_hot_streak"] = max(
                totals["max_hot_streak"], sensor["hot_streak"]
            )
            if raw_value is None:
                totals["read_error_samples"] += 1
            else:
                totals["raw_samples"] += 1
                totals["max_raw_c"] = (
                    raw_value if totals["max_raw_c"] is None
                    else max(totals["max_raw_c"], raw_value)
                )
            if sensor["valid"]:
                totals["valid_samples"] += 1
                totals["max_valid_c"] = (
                    raw_value if totals["max_valid_c"] is None
                    else max(totals["max_valid_c"], raw_value)
                )
            elif sensor["reason"] != "read_error":
                totals["invalid_samples"] += 1
        max_consecutive = max(max_consecutive, sample["consecutive_fault_rows"])
        last_sample = sample
    if last_sample is None or ce_baseline is None or ue_baseline is None:
        fail("terminal validation has no samples")
    summary = {
        "consecutive_fault_rows": last_sample["consecutive_fault_rows"],
        "decision": last_sample["decision"],
        "dimm_attempt_errors_total": sum(
            value["attempt_errors"] for value in sensor_totals.values()
        ),
        "dimm_invalid_samples_total": sum(
            value["invalid_samples"] for value in sensor_totals.values()
        ),
        "dimm_read_error_samples_total": sum(
            value["read_error_samples"] for value in sensor_totals.values()
        ),
        "edac_ce_baseline": ce_baseline,
        "edac_ce_delta": last_ce - ce_baseline,
        "edac_ce_last": last_ce,
        "edac_ue_baseline": ue_baseline,
        "edac_ue_delta": last_ue - ue_baseline,
        "edac_ue_last": last_ue,
        "max_consecutive_fault_rows": max_consecutive,
        "sample_count": len(raw_lines) - 1,
        "sensors": sensor_totals,
    }
    return summary, last_sample["decision"]


def validate_sampler_producer_receipt(
    receipt: Any, receipt_ascii: bytes, sampler_admission: Mapping[str, Any],
    raw_csv: bytes, raw_binding: Mapping[str, Any],
    validation_jsonl: bytes, validation_binding: Mapping[str, Any],
    producer_receipt_binding: Mapping[str, Any],
    pid_file_held_binding: Mapping[str, Any],
) -> Tuple[str, Dict[str, Any]]:
    if type(receipt) is not dict:
        fail("sampler producer receipt is not an object")
    exact_keys(receipt, SAMPLER_PRODUCER_RECEIPT_KEYS, "sampler producer receipt")
    exact_string(receipt["schema"], THERMAL_SAMPLER_SCHEMA, "sampler producer schema")
    self_hash = lower_hash(
        receipt["self_sha256_excluding_field"], "sampler producer self-hash"
    )
    preimage = dict(receipt)
    del preimage["self_sha256_excluding_field"]
    if sha256_bytes(canonical_bytes(preimage) + b"\n") != self_hash:
        fail("sampler producer receipt self-hash differs")
    if canonical_bytes(receipt) + b"\n" != receipt_ascii:
        fail("sampler producer receipt canonical bytes differ")
    if receipt["argv"] != sampler_admission["cmdline_argv"][5:]:
        fail("sampler producer receipt argv differs from admission")
    if (
        receipt["expected_output_owner_uid"] != sampler_admission["process_uid"]
        or receipt["pid"] != sampler_admission["pid"]
        or receipt["sampling"] != {
            "dimm_attempts": 5, "dimm_retry_delay_s": 0.01,
            "interval_s": 1.0,
        }
        or receipt["thresholds"] != THERMAL_SAMPLER_THRESHOLDS
    ):
        fail("sampler producer run authority differs")
    outcome = receipt["outcome"]
    outcome_laws = {
        "stopped": (0, ("SIGINT", "SIGTERM"), False),
        "thermal_abort": (3, (None, "SIGINT", "SIGTERM"), False),
        "telemetry_abort": (4, (None, "SIGINT", "SIGTERM"), False),
        "sampler_error": (5, (None, "SIGINT", "SIGTERM"), True),
        "edac_abort": (6, (None, "SIGINT", "SIGTERM"), False),
    }
    if outcome not in outcome_laws:
        fail("sampler producer outcome differs")
    expected_exit, expected_signal, require_errors = outcome_laws[outcome]
    errors = receipt["errors"]
    if (
        receipt["exit_code"] != expected_exit
        or (
            receipt["signal"] not in expected_signal
            if type(expected_signal) is tuple
            else receipt["signal"] != expected_signal
        )
        or type(errors) is not list or len(errors) > MAX_PUBLICATION_FAILURES
        or any(
            type(item) is not dict or set(item) != {"message", "phase", "type"}
            or any(type(item[key]) is not str or not item[key] or len(item[key]) > 4000
                   for key in ("message", "phase", "type"))
            for item in errors
        )
        or bool(errors) != require_errors
    ):
        fail("sampler producer outcome/exit/error binding differs")
    started_monotonic_ns = exact_int(
        receipt["started_monotonic_ns"], 0, MAX_INT63,
        "sampler producer started_monotonic_ns",
    )
    if (
        type(receipt["started_utc"]) is not str
        or not receipt["started_utc"]
        or len(receipt["started_utc"]) > 64
    ):
        fail("sampler producer start UTC timestamp differs")
    timestamp_errors = [
        item for item in errors
        if item["phase"] == "build_terminal_timestamps"
    ]
    finish_missing = (
        receipt["finished_monotonic_ns"] is None
        and receipt["finished_utc"] is None
    )
    if finish_missing:
        if outcome != "sampler_error" or len(timestamp_errors) != 1:
            fail("sampler producer missing finish timestamps are unlicensed")
    else:
        if (
            receipt["finished_monotonic_ns"] is None
            or receipt["finished_utc"] is None
            or timestamp_errors
        ):
            fail("sampler producer finish timestamp pair differs")
        finished_monotonic_ns = exact_int(
            receipt["finished_monotonic_ns"], 0, MAX_INT63,
            "sampler producer finished_monotonic_ns",
        )
        if finished_monotonic_ns < started_monotonic_ns:
            fail("sampler producer timestamps are reversed")
        if (
            type(receipt["finished_utc"]) is not str
            or not receipt["finished_utc"]
            or len(receipt["finished_utc"]) > 64
        ):
            fail("sampler producer finish UTC timestamp differs")
    for name in ("edac_ce_paths", "edac_ue_paths"):
        paths = receipt[name]
        if (
            type(paths) is not list or not paths
            or any(type(path) is not str or not Path(path).is_absolute()
                   for path in paths)
            or paths != sorted(set(paths))
        ):
            fail("sampler producer EDAC path roster differs")

    parent = sampler_admission["evidence_parent"]
    artifact_specs = (
        ("raw_csv", sampler_admission["csv_path"], raw_binding),
        (
            "validation_jsonl", sampler_admission["validation_jsonl"]["path"],
            validation_binding,
        ),
    )
    producer_artifact_bindings: Dict[str, Optional[Dict[str, Any]]] = {}
    for name, expected_path, actual_binding in artifact_specs:
        value = receipt[name]
        if type(value) is not dict or set(value) != {"binding", "destination", "path"}:
            fail("sampler producer " + name + " shape differs")
        exact_string(value["path"], expected_path, "sampler producer " + name + " path")
        seal_error = any(item["phase"] == "seal_" + name for item in errors)
        binding_value = value["binding"]
        if binding_value is None:
            if outcome != "sampler_error" or not seal_error:
                fail("sampler producer " + name + " binding is absent")
            binding = None
        else:
            binding = validate_sampler_binding(
                binding_value, "sampler producer " + name + " binding"
            )
            if seal_error or not canonical_equal(binding, actual_binding):
                fail("sampler producer " + name + " binding differs")
        producer_artifact_bindings[name] = binding
        validate_sampler_destination(
            value["destination"], expected_path, sampler_admission["process_uid"],
            parent, "sampler producer " + name + " destination",
        )
    raw_samples_preserved = receipt["raw_samples_preserved"]
    if (
        type(raw_samples_preserved) is not bool
        or raw_samples_preserved
        != (producer_artifact_bindings["raw_csv"] is not None)
    ):
        fail("sampler producer raw preservation flag differs")
    if (
        raw_binding["sha256"] != sha256_bytes(raw_csv)
        or raw_binding["size"] != len(raw_csv)
        or validation_binding["sha256"] != sha256_bytes(validation_jsonl)
        or validation_binding["size"] != len(validation_jsonl)
    ):
        fail("sampler terminal embedded artifact binding differs")
    admission_raw_stat = sampler_admission["csv_stat"]
    admission_validation = sampler_admission["validation_jsonl"]
    admission_validation_stat = admission_validation["stat"]
    if (
        raw_binding["mode"] != 0o444
        or raw_binding["nlink"] != 1
        or any(
            raw_binding[name] != admission_raw_stat[name]
            for name in ("device", "gid", "inode", "uid")
        )
        or raw_binding["size"] < sampler_admission["csv_bytes"]
        or sha256_bytes(raw_csv[:sampler_admission["csv_bytes"]])
        != sampler_admission["csv_sha256"]
        or validation_binding["mode"] != 0o444
        or validation_binding["nlink"] != 1
        or any(
            validation_binding[name] != admission_validation_stat[name]
            for name in ("device", "gid", "inode", "uid")
        )
        or validation_binding["size"] < admission_validation["bytes"]
        or sha256_bytes(validation_jsonl[:admission_validation["bytes"]])
        != admission_validation["sha256"]
    ):
        fail("sampler terminal artifacts differ from admission inodes/prefixes")

    pid_value = receipt["pid_file"]
    if type(pid_value) is not dict or set(pid_value) != {
        "binding", "destination", "path", "removed",
    }:
        fail("sampler producer PID receipt shape differs")
    exact_string(
        pid_value["path"], sampler_admission["pid_file"]["path"],
        "sampler producer PID path",
    )
    pid_binding = validate_sampler_binding(
        pid_value["binding"], "sampler producer PID binding"
    )
    admission_pid = sampler_admission["pid_file"]
    expected_pid_binding = {
        "device": admission_pid["stat"]["device"],
        "gid": admission_pid["stat"]["gid"],
        "inode": admission_pid["stat"]["inode"],
        "mode": admission_pid["stat"]["mode"],
        "nlink": admission_pid["stat"]["nlink"],
        "sha256": admission_pid["sha256"],
        "size": admission_pid["bytes"],
        "uid": admission_pid["stat"]["uid"],
    }
    removed = pid_value["removed"]
    unlink_error = any(
        item["phase"] == "unlink_pid_file" for item in errors
    )
    if (
        type(removed) is not bool
        or (
            not removed
            and (outcome != "sampler_error" or not unlink_error)
        )
        or not canonical_equal(pid_binding, expected_pid_binding)
        or pid_file_held_binding["device"] != pid_binding["device"]
        or pid_file_held_binding["inode"] != pid_binding["inode"]
        or pid_file_held_binding["nlink"] != 0
        or any(
            pid_file_held_binding[name] != pid_binding[name]
            for name in ("gid", "mode", "sha256", "size", "uid")
        )
    ):
        fail("sampler producer PID unlink/held binding differs")
    validate_sampler_destination(
        pid_value["destination"], admission_pid["path"],
        sampler_admission["process_uid"], parent,
        "sampler producer PID destination",
    )

    receipt_file = receipt["receipt_file"]
    if type(receipt_file) is not dict or set(receipt_file) != {
        "destination", "path", "reservation_binding",
    }:
        fail("sampler producer receipt-file shape differs")
    exact_string(
        receipt_file["path"], sampler_admission["receipt_file"]["path"],
        "sampler producer receipt-file path",
    )
    reservation = validate_sampler_binding(
        receipt_file["reservation_binding"],
        "sampler producer receipt reservation",
    )
    admission_receipt_file = sampler_admission["receipt_file"]
    expected_reservation = {
        "device": admission_receipt_file["stat"]["device"],
        "gid": admission_receipt_file["stat"]["gid"],
        "inode": admission_receipt_file["stat"]["inode"],
        "mode": admission_receipt_file["stat"]["mode"],
        "nlink": admission_receipt_file["stat"]["nlink"],
        "sha256": admission_receipt_file["sha256"],
        "size": admission_receipt_file["bytes"],
        "uid": admission_receipt_file["stat"]["uid"],
    }
    if (
        not canonical_equal(reservation, expected_reservation)
        or producer_receipt_binding["device"] != reservation["device"]
        or producer_receipt_binding["inode"] != reservation["inode"]
        or producer_receipt_binding["mode"] != 0o444
        or producer_receipt_binding["nlink"] != 1
        or producer_receipt_binding["uid"] != sampler_admission["process_uid"]
        or producer_receipt_binding["gid"] != sampler_admission["process_gid"]
        or producer_receipt_binding["sha256"] != sha256_bytes(receipt_ascii)
        or producer_receipt_binding["size"] != len(receipt_ascii)
    ):
        fail("sampler producer receipt reservation/final binding differs")
    validate_sampler_destination(
        receipt_file["destination"], admission_receipt_file["path"],
        sampler_admission["process_uid"], parent,
        "sampler producer receipt destination",
    )

    source = receipt["sampler_source"]
    if type(source) is not dict or set(source) != {
        "binding", "binding_finished", "destination", "expected_sha256",
        "path", "sha256",
    }:
        fail("sampler producer source receipt shape differs")
    source_binding = validate_sampler_binding(
        source["binding"], "sampler producer source binding"
    )
    source_finished = validate_sampler_binding(
        source["binding_finished"], "sampler producer source final binding"
    )
    admission_script = sampler_admission["script_stat"]
    if (
        source["path"] != sampler_admission["script_path"]
        or source["expected_sha256"] != sampler_admission["script_sha256"]
        or source["sha256"] != sampler_admission["script_sha256"]
        or not canonical_equal(source_binding, source_finished)
        or source_binding["sha256"] != sampler_admission["script_sha256"]
        or any(
            source_binding[name] != admission_script[name]
            for name in ("device", "gid", "inode", "mode", "nlink", "size", "uid")
        )
    ):
        fail("sampler producer source binding differs from admission")
    validate_sampler_destination(
        source["destination"], sampler_admission["script_path"], None, None,
        "sampler producer source destination",
    )

    stream_parts = split_terminal_streams(
        raw_csv, validation_jsonl, allow_unpaired=outcome == "sampler_error"
    )
    summary, final_decision = replay_terminal_validation(
        stream_parts["paired_raw"], stream_parts["validation_complete"],
        sampler_admission,
    )
    suffix_shape = stream_parts["suffix_shape"]
    if final_decision != "continue" and suffix_shape != (0, False, False):
        fail("sampler terminal evidence follows a terminal validation decision")
    if suffix_shape[0] == 1:
        unpaired_row = parse_csv_physical_line(
            stream_parts["raw_suffix"],
            "sampler terminal unpaired raw row",
        )
        if len(unpaired_row) != len(THERMAL_HEADER):
            fail("sampler terminal unpaired raw row width differs")
        unpaired_ns = exact_int(
            int(round(
                finite_scalar(
                    unpaired_row[1],
                    "sampler terminal unpaired raw timestamp",
                ) * 1_000_000_000.0
            )),
            0, MAX_INT63, "sampler terminal unpaired raw timestamp",
        )
        _, paired_timestamps = sampler_rows(stream_parts["paired_raw"])
        paired_end_ns = max(paired_timestamps.values())
        if (
            unpaired_ns <= paired_end_ns
            or receipt["finished_monotonic_ns"] is not None
            and unpaired_ns > receipt["finished_monotonic_ns"]
        ):
            fail("sampler terminal unpaired raw chronology differs")
        stream_parts["unpaired_raw_monotonic_ns"] = unpaired_ns
    else:
        stream_parts["unpaired_raw_monotonic_ns"] = None
    if outcome != "sampler_error" and not canonical_equal(
        receipt["summary"], summary
    ):
        fail("sampler producer summary does not replay")
    if outcome == "sampler_error" and receipt["summary"] is not None and type(
        receipt["summary"]
    ) is not dict:
        fail("sampler-error producer summary shape differs")
    raw_lines = stream_parts["raw_complete"].splitlines(keepends=True)
    observed_complete_cpu_max = max(
        finite_scalar(
            parse_csv_physical_line(raw, "terminal CPU maximum row")[4],
            "terminal CPU maximum",
        )
        for raw in raw_lines[1:]
    )
    if outcome == "sampler_error":
        reported_cpu_max = receipt["cpu_tctl_max_c"]
        if (
            type(reported_cpu_max) not in (int, float)
            or not math.isfinite(reported_cpu_max)
            or reported_cpu_max < observed_complete_cpu_max
        ):
            fail("sampler-error producer CPU maximum contradicts retained rows")
    elif not canonical_equal(
        receipt["cpu_tctl_max_c"], observed_complete_cpu_max
    ):
        fail("sampler producer CPU maximum does not replay")
    expected_final_decision = "continue" if outcome == "stopped" else outcome
    if outcome != "sampler_error" and final_decision != expected_final_decision:
        fail("sampler producer final validation decision differs from outcome")
    return outcome, stream_parts


def validate_sampler_terminal_receipt(
    terminal: Any, sampler_admission: Mapping[str, Any],
    child_start: Optional[int], child_reap: Optional[int],
) -> None:
    if type(terminal) is not dict:
        fail("sampler terminal evidence is not an object")
    exact_keys(terminal, SAMPLER_TERMINAL_KEYS, "sampler terminal evidence")
    exact_string(
        terminal["schema"], SAMPLER_TERMINAL_SCHEMA,
        "sampler terminal evidence schema",
    )
    exact_string(
        terminal["terminal_status"], "invalid",
        "sampler terminal evidence status",
    )
    exact_string(
        terminal["process_exit_observation"],
        "linux-pidfd-readable-nonchild-v1",
        "sampler terminal process-exit observation",
    )
    exit_observed = exact_int(
        terminal["process_exit_observed_monotonic_ns"], 0, MAX_INT63,
        "sampler terminal process-exit timestamp",
    )
    admission_digest = lower_hash(
        terminal["admission_receipt_sha256"],
        "sampler terminal admission digest",
    )
    if admission_digest != sha256_bytes(canonical_bytes(sampler_admission)):
        fail("sampler terminal admission binding differs")
    ascii_fields = (
        ("producer_receipt_ascii", MAX_SAMPLER_RECEIPT_BYTES),
        ("raw_csv_ascii", MAX_SAMPLER_CSV_BYTES),
        ("validation_jsonl_ascii", MAX_SAMPLER_VALIDATION_BYTES),
        ("window_csv_ascii", MAX_THERMAL_WINDOW_BYTES),
        ("window_validation_jsonl_ascii", MAX_SAMPLER_VALIDATION_BYTES),
    )
    decoded: Dict[str, bytes] = {}
    for name, maximum in ascii_fields:
        value = terminal[name]
        if type(value) is not str:
            fail("sampler terminal " + name + " is not text")
        try:
            raw = value.encode("ascii")
        except UnicodeEncodeError:
            fail("sampler terminal " + name + " is not ASCII")
        if not raw or len(raw) > maximum:
            fail("sampler terminal " + name + " exceeds its evidence bound")
        decoded[name] = raw
    producer_receipt = parse_canonical_json_line(
        decoded["producer_receipt_ascii"], "sampler terminal producer receipt"
    )
    producer_binding = validate_sampler_binding(
        terminal["producer_receipt_binding"], "sampler terminal receipt binding"
    )
    raw_binding = validate_sampler_binding(
        terminal["raw_csv_binding"], "sampler terminal raw binding"
    )
    validation_binding = validate_sampler_binding(
        terminal["validation_jsonl_binding"],
        "sampler terminal validation binding",
    )
    pid_binding = validate_sampler_binding(
        terminal["pid_file_held_binding"], "sampler terminal held PID binding"
    )
    _, stream_parts = validate_sampler_producer_receipt(
        producer_receipt, decoded["producer_receipt_ascii"], sampler_admission,
        decoded["raw_csv_ascii"], raw_binding,
        decoded["validation_jsonl_ascii"], validation_binding,
        producer_binding, pid_binding,
    )
    suffixes = terminal["stream_suffixes"]
    if type(suffixes) is not dict:
        fail("sampler terminal stream-suffix receipt is not an object")
    exact_keys(
        suffixes, SAMPLER_TERMINAL_STREAM_SUFFIX_KEYS,
        "sampler terminal stream-suffix receipt",
    )
    for name in (
        "paired_sample_count", "raw_unpaired_suffix_bytes",
        "validation_complete_prefix_bytes", "validation_partial_suffix_bytes",
    ):
        exact_int(
            suffixes[name], 0, MAX_SAMPLER_VALIDATION_BYTES,
            "sampler terminal stream-suffix " + name,
        )
    for name in (
        "raw_unpaired_suffix_sha256", "validation_complete_prefix_sha256",
        "validation_partial_suffix_sha256",
    ):
        lower_hash(suffixes[name], "sampler terminal stream-suffix " + name)
    if not canonical_equal(suffixes, stream_parts["suffix_receipt"]):
        fail("sampler terminal stream-suffix receipt does not replay")
    coverage = terminal["coverage"]
    if type(coverage) is not dict:
        fail("sampler terminal coverage is not an object")
    exact_keys(coverage, SAMPLER_TERMINAL_COVERAGE_KEYS, "sampler terminal coverage")
    start = exact_int(
        coverage["window_start_monotonic_ns"], 0, MAX_INT63,
        "sampler terminal window start",
    )
    end = exact_int(
        coverage["window_end_monotonic_ns"], start, MAX_INT63,
        "sampler terminal window end",
    )
    producer_finish = producer_receipt["finished_monotonic_ns"]
    if producer_finish is None:
        chronology_valid = (
            producer_receipt["started_monotonic_ns"]
            <= start <= end <= exit_observed
        )
    else:
        chronology_valid = (
            producer_receipt["started_monotonic_ns"] <= start <= end
            <= producer_finish <= exit_observed
        )
    if (
        not chronology_valid
        or (
            stream_parts["unpaired_raw_monotonic_ns"] is not None
            and stream_parts["unpaired_raw_monotonic_ns"] > exit_observed
        )
    ):
        fail("sampler terminal producer/window/exit chronology differs")
    if (
        coverage["child_start_monotonic_ns"] != child_start
        or coverage["child_reap_monotonic_ns"] != child_reap
        or type(coverage["covers_child_interval"]) is not bool
    ):
        fail("sampler terminal child interval binding differs")
    expected_covers = bool(
        child_start is not None and child_reap is not None
        and start <= child_start <= child_reap <= end
    )
    expected_shortfall = (
        None if child_reap is None else max(0, child_reap - end)
    )
    if (
        coverage["covers_child_interval"] != expected_covers
        or coverage["coverage_shortfall_ns"] != expected_shortfall
    ):
        fail("sampler terminal coverage computation differs")
    paired_raw = stream_parts["paired_raw"]
    physical, timestamps = sampler_rows(paired_raw)
    start_lines = [index for index, value in timestamps.items() if value == start]
    end_lines = [index for index, value in timestamps.items() if value == end]
    if len(start_lines) != 1 or len(end_lines) != 1 or start_lines[0] > end_lines[0]:
        fail("sampler terminal window endpoints differ")
    expected_window = physical[0] + b"".join(
        physical[start_lines[0]:end_lines[0] + 1]
    )
    if decoded["window_csv_ascii"] != expected_window:
        fail("sampler terminal retained CSV window differs")
    expected_validation_window = select_validation_window(
        stream_parts["validation_complete"], expected_window,
        "sampler terminal validation selection",
    )
    if decoded["window_validation_jsonl_ascii"] != expected_validation_window:
        fail("sampler terminal retained validation window differs")


def parse_csv_physical_line(raw: bytes, where: str) -> List[str]:
    if not raw.endswith(b"\n") or b"\r" in raw or len(raw) > 16384:
        fail(f"{where} physical line differs")
    try:
        text = raw[:-1].decode("ascii")
        rows = list(csv.reader([text], strict=True))
    except (UnicodeDecodeError, csv.Error) as exc:
        fail(f"{where} is not canonical ASCII CSV: {exc}")
    if len(rows) != 1:
        fail(f"{where} CSV row count differs")
    if (",".join(rows[0]) + "\n").encode("ascii") != raw:
        fail(f"{where} is not canonical unquoted CSV")
    return rows[0]


def sampler_rows(data: bytes) -> Tuple[List[bytes], Dict[int, int]]:
    physical = data.splitlines(keepends=True)
    if not physical:
        fail("sampler CSV is empty")
    if tuple(parse_csv_physical_line(physical[0], "sampler header")) != THERMAL_HEADER:
        fail("sampler CSV header differs")
    timestamps: Dict[int, int] = {}
    for index, raw in enumerate(physical[1:], start=1):
        try:
            row = parse_csv_physical_line(raw, f"sampler row {index + 1}")
            if len(row) != len(THERMAL_HEADER):
                continue
            value = float(row[1])
            if not math.isfinite(value) or value < 0.0:
                continue
            timestamp = int(round(value * 1_000_000_000.0))
            exact_int(
                timestamp, 0, MAX_INT63,
                f"sampler row {index + 1} monotonic timestamp",
            )
            timestamps[index] = timestamp
        except (ValidationError, ValueError, OverflowError):
            continue
    return physical, timestamps


def wait_for_sampler_timestamp(
    path: Path, deadline: float, *, greater_than: Optional[int] = None,
    at_or_after: Optional[int] = None,
) -> int:
    if (greater_than is None) == (at_or_after is None):
        fail("sampler wait requires exactly one timestamp predicate")
    while True:
        data, _ = sampler_snapshot(path)
        _, timestamps = sampler_rows(data)
        candidates = sorted(timestamps.values())
        if greater_than is not None:
            matches = [value for value in candidates if value > greater_than]
        else:
            matches = [value for value in candidates if value >= at_or_after]
        if matches:
            return matches[0]
        if time.monotonic() >= deadline:
            fail("sampler did not append a required endpoint before deadline")
        time.sleep(0.05)


def wait_for_sampler_timestamp_supervised(
    path: Path, deadline: float, handles: Mapping[str, Any],
    *, at_or_after: int,
) -> Tuple[Optional[int], str]:
    while True:
        event = poll_sampler_supervision(handles)
        if event != "none":
            return None, event
        data, _ = sampler_snapshot(path)
        _, timestamps = sampler_rows(data)
        matches = sorted(
            value for value in timestamps.values() if value >= at_or_after
        )
        event = poll_sampler_supervision(handles)
        if event != "none":
            return None, event
        if matches:
            return matches[0], "none"
        if time.monotonic() >= deadline:
            fail("sampler did not append a supervised endpoint before deadline")
        time.sleep(min(0.05, max(0.0, deadline - time.monotonic())))


def wait_for_sampler_validation_timestamp_supervised(
    handles: Mapping[str, Any], timestamp_ns: int, deadline: float,
) -> str:
    exact_int(
        timestamp_ns, 0, MAX_INT63,
        "sampler admission validation timestamp",
    )
    while True:
        event = poll_sampler_supervision(handles)
        if event != "none":
            return event
        observed = handles.get("validation_csv_timestamp_ns", [])
        if timestamp_ns in observed:
            # A terminal latch or later terminal validation record must win
            # over a just-observed continue record before admission returns.
            event = poll_sampler_supervision(handles)
            if event != "none":
                return event
            if timestamp_ns not in handles.get(
                "validation_csv_timestamp_ns", []
            ):
                fail("sampler validation state changed during admission")
            return "none"
        last_timestamp = handles.get("validation_last_csv_timestamp_ns")
        if last_timestamp is not None and last_timestamp > timestamp_ns:
            fail("sampler validation skipped the admission CSV endpoint")
        if time.monotonic() >= deadline:
            fail("sampler did not validate the admission CSV endpoint")
        time.sleep(min(0.05, max(0.0, deadline - time.monotonic())))


def process_start_ticks(pid: int) -> int:
    data = Path(f"/proc/{pid}/stat").read_text(encoding="ascii")
    close = data.rfind(")")
    fields = data[close + 2:].split() if close >= 0 else []
    if len(fields) <= 19 or not fields[19].isdecimal():
        fail("sampler /proc stat is malformed")
    return exact_int(int(fields[19]), 1, MAX_UINT64, "sampler start ticks")


def read_bounded_proc_vector(path: Path, where: str) -> bytes:
    fd = os.open(
        str(path), nonblocking_read_flags(where)
    )
    try:
        data = bytearray()
        while len(data) <= MAX_PROC_VECTOR_BYTES:
            block = os.read(fd, min(16384, MAX_PROC_VECTOR_BYTES + 1 - len(data)))
            if not block:
                break
            data.extend(block)
        if len(data) > MAX_PROC_VECTOR_BYTES:
            fail(f"{where} exceeds its 64-KiB bound")
        return bytes(data)
    finally:
        os.close(fd)


PROCESS_SECURITY_KEYS = {
    "cap_ambient", "cap_bounding", "cap_effective", "cap_inheritable",
    "cap_permitted", "gids", "groups", "no_new_privs", "schema", "uids",
}


def parse_process_security(data: bytes, where: str) -> Dict[str, Any]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        fail(f"{where} is not ASCII")
    wanted = {
        "Uid": "uids", "Gid": "gids", "Groups": "groups",
        "CapInh": "cap_inheritable", "CapPrm": "cap_permitted",
        "CapEff": "cap_effective", "CapBnd": "cap_bounding",
        "CapAmb": "cap_ambient", "NoNewPrivs": "no_new_privs",
    }
    observed: Dict[str, str] = {}
    for raw_line in text.splitlines():
        if ":" not in raw_line:
            continue
        raw_name, raw_value = raw_line.split(":", 1)
        if raw_name not in wanted:
            continue
        name = wanted[raw_name]
        if name in observed:
            fail(f"{where} contains duplicate {raw_name}")
        observed[name] = raw_value.strip()
    if set(observed) != set(wanted.values()):
        fail(f"{where} security field roster differs")

    def decimal_vector(raw: str, count: Optional[int], name: str) -> List[int]:
        fields = raw.split() if raw else []
        if count is not None and len(fields) != count:
            fail(f"{where} {name} field count differs")
        values: List[int] = []
        for field in fields:
            if (
                not field.isascii() or not field.isdecimal()
                or (len(field) > 1 and field.startswith("0"))
            ):
                fail(f"{where} {name} is not canonical decimal")
            values.append(exact_int(int(field), 0, MAX_UINT32, f"{where} {name}"))
        return values

    uids = decimal_vector(observed["uids"], 4, "Uid")
    gids = decimal_vector(observed["gids"], 4, "Gid")
    groups = decimal_vector(observed["groups"], None, "Groups")
    if groups != sorted(set(groups)):
        fail(f"{where} supplementary group roster is not sorted and unique")
    capabilities: Dict[str, str] = {}
    for name in (
        "cap_inheritable", "cap_permitted", "cap_effective",
        "cap_bounding", "cap_ambient",
    ):
        value = observed[name]
        if re.fullmatch(r"[0-9a-f]{16}", value) is None:
            fail(f"{where} {name} is not canonical 64-bit lowercase hex")
        capabilities[name] = value
    nnp_text = observed["no_new_privs"]
    if nnp_text not in ("0", "1"):
        fail(f"{where} NoNewPrivs differs")
    return {
        **capabilities,
        "gids": gids,
        "groups": groups,
        "no_new_privs": int(nnp_text),
        "schema": PROCESS_SECURITY_SCHEMA,
        "uids": uids,
    }


def validate_process_security(
    value: Any, expected_uid: int, expected_gid: int,
    expected_groups: Sequence[int], where: str,
) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, PROCESS_SECURITY_KEYS, where)
    exact_string(value["schema"], PROCESS_SECURITY_SCHEMA, where + ".schema")
    exact_int(expected_uid, 0, MAX_UINT32, where + " expected UID")
    exact_int(expected_gid, 0, MAX_UINT32, where + " expected GID")
    groups = exact_int_list(
        list(expected_groups), where + " expected groups", sorted_unique=True,
    )
    if (
        value["uids"] != [expected_uid] * 4
        or value["gids"] != [expected_gid] * 4
        or value["groups"] != groups
        or value["no_new_privs"] != 1
        or any(
            value[name] != "0000000000000000"
            for name in (
                "cap_inheritable", "cap_permitted", "cap_effective",
                "cap_bounding", "cap_ambient",
            )
        )
    ):
        fail(f"{where} credentials/capabilities differ from authorization")


def read_process_security(pid: int, where: str) -> Dict[str, Any]:
    exact_int(pid, 1, MAX_INT63, where + " PID")
    return parse_process_security(
        read_bounded_proc_vector(Path(f"/proc/{pid}/status"), where + " status"),
        where,
    )


def parse_sealed_environment(data: bytes, where: str) -> Dict[str, str]:
    if not data or not data.endswith(b"\0"):
        fail(f"{where} is not terminal-NUL framed")
    environment: Dict[str, str] = {}
    for index, raw in enumerate(data[:-1].split(b"\0")):
        if not raw or b"=" not in raw:
            fail(f"{where} entry {index} framing differs")
        raw_name, raw_value = raw.split(b"=", 1)
        try:
            name = raw_name.decode("ascii")
            value = raw_value.decode("ascii")
        except UnicodeDecodeError:
            fail(f"{where} entry {index} is not ASCII")
        if (
            re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name) is None
            or name in environment
        ):
            fail(f"{where} entry {index} name differs")
        environment[name] = value
    if not canonical_equal(environment, SEALED_LAUNCH_ENVIRONMENT):
        fail(f"{where} exact environment differs")
    return environment


def expected_sampler_argv(args: argparse.Namespace) -> List[str]:
    return [
        EXPECTED_SAMPLER_PYTHON,
        *EXPECTED_SAMPLER_PYTHON_FLAGS,
        str(args.sampler_script),
        "--csv", str(args.sampler_csv),
        "--pid-file", str(args.sampler_pid_file),
        "--validation-jsonl", str(args.sampler_validation_jsonl),
        "--receipt", str(args.sampler_receipt),
        "--expected-source-sha256", args.expected_sampler_script_sha256,
        "--expected-output-owner-uid", str(args.expected_sampler_uid),
        "--interval", EXPECTED_SAMPLER_INTERVAL_TEXT,
        "--dimm-attempts", EXPECTED_SAMPLER_DIMM_ATTEMPTS_TEXT,
        "--dimm-retry-delay", EXPECTED_SAMPLER_DIMM_RETRY_DELAY_TEXT,
    ]


def expected_sampler_cmdline_bytes(args: argparse.Namespace) -> bytes:
    return b"".join(os.fsencode(value) + b"\0" for value in expected_sampler_argv(args))


def validate_sampler_cmdline(
    data: bytes, args: argparse.Namespace, where: str,
) -> List[str]:
    if not data or not data.endswith(b"\0") or b"\0\0" in data:
        fail(f"{where} is not an exact terminal-NUL argv vector")
    expected = expected_sampler_cmdline_bytes(args)
    if data != expected:
        fail(f"{where} differs from the frozen ordered sampler launch")
    exact_string(
        sha256_bytes(data), args.expected_sampler_cmdline_sha256,
        where + " SHA-256",
    )
    return expected_sampler_argv(args)


def hash_regular_path(
    path: Path, where: str, *, single_link: bool = True,
    max_size: int = MAX_INT63, deadline: Optional[float] = None,
) -> Tuple[str, os.stat_result]:
    fd = os.open(str(path), nonblocking_read_flags(where))
    try:
        info = os.fstat(fd)
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_nlink < 1
            or (single_link and info.st_nlink != 1)
            or info.st_size < 0 or info.st_size > max_size
        ):
            fail(f"{where} regular-file link policy differs")
        digest = file_sha256_fd(fd, info.st_size, deadline)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            not same_file_receipt(info, after)
            or not same_file_receipt(after, named)
        ):
            fail(f"{where} changed during held-FD hashing")
        return digest, after
    finally:
        os.close(fd)


def make_sampler_attestation(
    args: argparse.Namespace, start_ns: int, end_ns: int,
    health_modules: SealedHealthModules, deadline: Optional[float] = None,
    *, allow_stale_evidence: bool = False,
    stale_stream_capture: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    if deadline is not None and time.monotonic() >= deadline:
        fail("sampler attestation reached its global deadline")
    health_api = health_modules.native
    evidence_paths = (
        args.sampler_csv, args.sampler_pid_file,
        args.sampler_validation_jsonl, args.sampler_receipt,
    )
    for path, where in (
        (args.sampler_script, "sampler script"),
        (args.sampler_csv, "sampler CSV"),
        (args.sampler_pid_file, "sampler PID file"),
        (args.sampler_validation_jsonl, "sampler validation JSONL"),
        (args.sampler_receipt, "sampler terminal receipt"),
    ):
        if path.resolve(strict=True) != path:
            fail(f"{where} path is not canonical")
    evidence_parent = args.sampler_csv.parent
    if (
        any(path.parent != evidence_parent for path in evidence_paths)
        or len(set(evidence_paths)) != len(evidence_paths)
    ):
        fail("sampler evidence files do not share one distinct canonical parent")
    parent_before = os.stat(str(evidence_parent), follow_symlinks=False)
    if (
        not stat.S_ISDIR(parent_before.st_mode)
        or stat.S_IMODE(parent_before.st_mode) != 0o700
        or parent_before.st_uid != args.expected_sampler_uid
        or parent_before.st_gid != args.expected_sampler_gid
    ):
        fail("sampler evidence parent mode/owner differs")
    health_api._verify_live_sampler_process(
        args.sampler_pid, EXPECTED_SAMPLER_CPU,
        args.expected_sampler_process_start_ticks,
        args.sampler_script, args.sampler_csv,
    )
    actual_start_ticks = process_start_ticks(args.sampler_pid)
    exact_int(
        actual_start_ticks, args.expected_sampler_process_start_ticks,
        args.expected_sampler_process_start_ticks, "sampler process start ticks",
    )
    security_before = read_process_security(args.sampler_pid, "sampler pre-attestation")
    validate_process_security(
        security_before, args.expected_sampler_uid, args.expected_sampler_gid,
        [args.expected_sampler_i2c_gid], "sampler pre-attestation security",
    )
    affinity_before = sorted(os.sched_getaffinity(args.sampler_pid))
    if affinity_before != [EXPECTED_SAMPLER_CPU]:
        fail("sampler pre-attestation affinity differs")
    script_hash, script_info = hash_regular_path(
        args.sampler_script, "sampler script",
        max_size=MAX_SOURCE_FILE_BYTES, deadline=deadline,
    )
    exact_string(
        script_hash, args.expected_sampler_script_sha256,
        "sampler script SHA-256",
    )
    csv_data, csv_info = sampler_snapshot(args.sampler_csv)
    pid_file, pid_data, pid_info = sampler_artifact_receipt(
        args.sampler_pid_file, 64, 0o444, args.expected_sampler_uid,
        args.expected_sampler_gid, "sampler PID file",
    )
    validation_jsonl, validation_data, validation_info = sampler_artifact_receipt(
        args.sampler_validation_jsonl, MAX_SAMPLER_VALIDATION_BYTES, 0o444,
        args.expected_sampler_uid, args.expected_sampler_gid,
        "sampler validation JSONL",
    )
    receipt_file, receipt_data, receipt_info = sampler_artifact_receipt(
        args.sampler_receipt, MAX_SAMPLER_RECEIPT_BYTES, 0o444,
        args.expected_sampler_uid, args.expected_sampler_gid,
        "sampler terminal receipt reservation",
    )
    if pid_data != (str(args.sampler_pid) + "\n").encode("ascii"):
        fail("sampler PID file content differs")
    if receipt_data:
        fail("live sampler terminal receipt reservation is not empty")
    validate_validation_stream_header(
        validation_data, args, "sampler validation JSONL"
    )
    validation_header_raw = validation_data[:validation_data.find(b"\n") + 1]
    if (
        stat.S_IMODE(script_info.st_mode) != 0o444
        or script_info.st_uid != args.expected_sampler_uid
        or script_info.st_gid != args.expected_sampler_gid
        or stat.S_IMODE(csv_info.st_mode) != 0o600
        or csv_info.st_uid != args.expected_sampler_uid
        or csv_info.st_gid != args.expected_sampler_gid
    ):
        fail("sampler script/CSV mode or owner differs")
    exact_int(
        csv_info.st_dev, args.expected_sampler_csv_device,
        args.expected_sampler_csv_device, "sampler CSV device",
    )
    exact_int(
        csv_info.st_ino, args.expected_sampler_csv_inode,
        args.expected_sampler_csv_inode, "sampler CSV inode",
    )
    for info, expected_device, expected_inode, where in (
        (
            pid_info, args.expected_sampler_pid_file_device,
            args.expected_sampler_pid_file_inode, "sampler PID file",
        ),
        (
            validation_info, args.expected_sampler_validation_device,
            args.expected_sampler_validation_inode, "sampler validation JSONL",
        ),
        (
            receipt_info, args.expected_sampler_receipt_device,
            args.expected_sampler_receipt_inode, "sampler terminal receipt",
        ),
    ):
        exact_int(info.st_dev, expected_device, expected_device, where + " device")
        exact_int(info.st_ino, expected_inode, expected_inode, where + " inode")
    cmdline = read_bounded_proc_vector(
        Path(f"/proc/{args.sampler_pid}/cmdline"), "sampler command line"
    )
    cmdline_argv = validate_sampler_cmdline(cmdline, args, "sampler command line")
    environ = read_bounded_proc_vector(
        Path(f"/proc/{args.sampler_pid}/environ"), "sampler environment"
    )
    environment = parse_sealed_environment(environ, "sampler environment")
    exact_string(
        sha256_bytes(environ), args.expected_sampler_environ_sha256,
        "sampler environment SHA-256",
    )
    proc_executable = Path(f"/proc/{args.sampler_pid}/exe")
    executable = proc_executable.resolve(strict=True)
    executable_hash_named, executable_info_named = hash_regular_path(
        executable, "sampler named executable", single_link=False,
        max_size=MAX_BINARY_BYTES, deadline=deadline,
    )
    executable_fd = os.open(
        str(proc_executable),
        nonblocking_read_flags("sampler live executable", nofollow=False),
    )
    try:
        executable_info = os.fstat(executable_fd)
        if (
            executable_info.st_size < 1
            or executable_info.st_size > MAX_BINARY_BYTES
        ):
            fail("sampler live executable size differs")
        executable_hash = file_sha256_fd(
            executable_fd, executable_info.st_size, deadline
        )
        executable_info_after = os.fstat(executable_fd)
        if (
            not stat.S_ISREG(executable_info.st_mode)
            or not same_file_receipt(executable_info, executable_info_named)
            or not same_file_receipt(executable_info, executable_info_after)
            or executable_hash != executable_hash_named
            or stat.S_IMODE(executable_info.st_mode) != 0o755
            or executable_info.st_uid != 0
            or executable_info.st_gid != 0
        ):
            fail("sampler live executable differs from its canonical path")
    finally:
        os.close(executable_fd)
    exact_string(
        executable_hash, args.expected_sampler_executable_sha256,
        "sampler executable SHA-256",
    )
    proc_info = Path(f"/proc/{args.sampler_pid}").stat()
    exact_int(
        proc_info.st_uid, args.expected_sampler_uid,
        args.expected_sampler_uid, "sampler process UID",
    )
    exact_int(
        proc_info.st_gid, args.expected_sampler_gid,
        args.expected_sampler_gid, "sampler process GID",
    )
    security_after = read_process_security(args.sampler_pid, "sampler post-attestation")
    validate_process_security(
        security_after, args.expected_sampler_uid, args.expected_sampler_gid,
        [args.expected_sampler_i2c_gid], "sampler post-attestation security",
    )
    affinity_after = sorted(os.sched_getaffinity(args.sampler_pid))
    parent_after = os.stat(str(evidence_parent), follow_symlinks=False)
    if (
        process_start_ticks(args.sampler_pid) != actual_start_ticks
        or not canonical_equal(security_before, security_after)
        or affinity_after != affinity_before
        or not same_file_receipt(parent_before, parent_after)
    ):
        fail("sampler process identity changed during attestation")
    settled_streams = settle_live_sampler_streams(
        args, health_api, deadline, "sampler final live streams",
        require_fresh=not allow_stale_evidence,
        allow_stale_suffixes=allow_stale_evidence,
    )
    csv_data = (
        settled_streams["paired_raw"]
        if allow_stale_evidence else settled_streams["csv_data"]
    )
    csv_info = settled_streams["csv_info"]
    pid_file = settled_streams["pid_file"]
    receipt_file = settled_streams["receipt_file"]
    validation_data = settled_streams["validation_data"]
    validation_jsonl = settled_streams["validation_jsonl"]
    validation_header_raw = validation_data[:validation_data.find(b"\n") + 1]
    if deadline is not None and time.monotonic() >= deadline:
        fail("sampler attestation reached its global deadline")
    if allow_stale_evidence:
        end_ns = settled_streams["paired_end_monotonic_ns"]
    if start_ns > end_ns or not csv_data:
        fail("sampler attestation window differs")
    attestation = {
        "cmdline_argv": cmdline_argv,
        "cmdline_sha256": args.expected_sampler_cmdline_sha256,
        "cpu": EXPECTED_SAMPLER_CPU,
        "csv_bytes": len(csv_data),
        "csv_device": csv_info.st_dev,
        "csv_inode": csv_info.st_ino,
        "csv_path": str(args.sampler_csv),
        "csv_sha256": sha256_bytes(csv_data),
        "csv_stat": stat_receipt(csv_info),
        "evidence_parent": {
            "path": str(evidence_parent), "stat": stat_receipt(parent_after),
        },
        "environ_sha256": sha256_bytes(environ),
        "environment": environment,
        "environment_sha256": sha256_bytes(canonical_bytes(environment)),
        "executable_device": executable_info.st_dev,
        "executable_inode": executable_info.st_ino,
        "executable_path": str(executable),
        "executable_sha256": executable_hash,
        "executable_stat": stat_receipt(executable_info),
        "pid": args.sampler_pid,
        "pid_file": pid_file,
        "process_affinity": affinity_after,
        "process_gid": proc_info.st_gid,
        "process_security": security_after,
        "process_start_ticks": actual_start_ticks,
        "process_uid": proc_info.st_uid,
        "receipt_file": receipt_file,
        "schema": SAMPLER_SCHEMA,
        "script_device": script_info.st_dev,
        "script_inode": script_info.st_ino,
        "script_path": str(args.sampler_script),
        "script_sha256": script_hash,
        "script_stat": stat_receipt(script_info),
        "terminal_status": "ok",
        "window_end_monotonic_ns": end_ns,
        "window_start_monotonic_ns": start_ns,
        "validation_header_ascii": validation_header_raw.decode("ascii"),
        "validation_jsonl": validation_jsonl,
    }
    if stale_stream_capture is not None:
        if not allow_stale_evidence:
            fail("sampler stream capture requires stale evidence mode")
        stale_stream_capture.clear()
        stale_stream_capture.update({
            "csv_full_data": settled_streams["csv_full_data"],
            "csv_info": settled_streams["csv_info"],
            "snapshot_observed_monotonic_ns": settled_streams[
                "snapshot_observed_monotonic_ns"
            ],
            "validation_data": settled_streams["validation_data"],
            "validation_info": settled_streams["validation_info"],
        })
    return attestation


def sampler_admission_stat(
    sampler: Mapping[str, Any], name: str,
) -> Mapping[str, Any]:
    if name == "raw_csv":
        return sampler["csv_stat"]
    return sampler[name]["stat"]


def consume_sampler_validation_bytes(
    handles: Mapping[str, Any], appended: bytes, where: str,
) -> str:
    buffer = handles.get("validation_buffer", b"") + appended
    lines = buffer.splitlines(keepends=True)
    tail = b""
    if lines and not lines[-1].endswith(b"\n"):
        tail = lines.pop()
    previous_index = handles.get("validation_sample_index")
    previous_monotonic_ns = handles.get("validation_last_monotonic_ns")
    observed_monotonic_ns = list(
        handles.get("validation_monotonic_ns", [])
    )
    previous_csv_timestamp_ns = handles.get(
        "validation_last_csv_timestamp_ns"
    )
    observed_csv_timestamp_ns = list(
        handles.get("validation_csv_timestamp_ns", [])
    )
    for offset, raw in enumerate(lines):
        sample = parse_validation_sample(
            parse_canonical_json_line(raw, f"{where} record {offset}"),
            f"{where} record {offset}",
        )
        if (
            previous_index is not None
            and sample["sample_index"] != previous_index + 1
        ):
            fail(where + " sample indices are not contiguous")
        sample_monotonic_ns = exact_int(
            int(round(float(sample["monotonic_s"]) * 1_000_000_000.0)),
            0, MAX_INT63, where + " sample monotonic timestamp",
        )
        sample_csv_timestamp_ns = exact_int(
            int(round(
                float(f"{float(sample['monotonic_s']):.6f}")
                * 1_000_000_000.0
            )),
            0, MAX_INT63, where + " sample CSV timestamp",
        )
        if (
            previous_monotonic_ns is not None
            and sample_monotonic_ns <= previous_monotonic_ns
        ):
            fail(where + " sample timestamps are not strictly increasing")
        if (
            previous_csv_timestamp_ns is not None
            and sample_csv_timestamp_ns <= previous_csv_timestamp_ns
        ):
            fail(where + " CSV timestamps are not strictly increasing")
        previous_index = sample["sample_index"]
        previous_monotonic_ns = sample_monotonic_ns
        previous_csv_timestamp_ns = sample_csv_timestamp_ns
        observed_monotonic_ns.append(sample_monotonic_ns)
        observed_csv_timestamp_ns.append(sample_csv_timestamp_ns)
        if sample["decision"] != "continue":
            if type(handles) is dict:
                handles["validation_buffer"] = tail
                handles["validation_last_monotonic_ns"] = previous_monotonic_ns
                handles["validation_monotonic_ns"] = observed_monotonic_ns
                handles["validation_last_csv_timestamp_ns"] = (
                    previous_csv_timestamp_ns
                )
                handles["validation_csv_timestamp_ns"] = (
                    observed_csv_timestamp_ns
                )
                handles["validation_sample_index"] = previous_index
            return "sampler-terminal-decision:" + sample["decision"]
    if type(handles) is dict:
        handles["validation_buffer"] = tail
        handles["validation_last_monotonic_ns"] = previous_monotonic_ns
        handles["validation_monotonic_ns"] = observed_monotonic_ns
        handles["validation_last_csv_timestamp_ns"] = (
            previous_csv_timestamp_ns
        )
        handles["validation_csv_timestamp_ns"] = observed_csv_timestamp_ns
        handles["validation_sample_index"] = previous_index
    return "none"


def open_sampler_admission_handles(
    args: argparse.Namespace, sampler: Mapping[str, Any],
) -> Dict[str, Any]:
    if not hasattr(os, "pidfd_open"):
        fail("sampler supervision requires Linux pidfd_open")
    parent_fd = -1
    pidfd = -1
    file_fds: Dict[str, int] = {}
    try:
        flags = (
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        )
        parent_fd = os.open(sampler["evidence_parent"]["path"], flags)
        parent_info = os.fstat(parent_fd)
        if stat_receipt(parent_info) != sampler["evidence_parent"]["stat"]:
            fail("sampler held evidence parent differs from admission")
        path_roster = {
            "raw_csv": args.sampler_csv,
            "pid_file": args.sampler_pid_file,
            "validation_jsonl": args.sampler_validation_jsonl,
            "receipt_file": args.sampler_receipt,
        }
        read_flags = nonblocking_read_flags("sampler held evidence")
        for name, path in path_roster.items():
            fd = os.open(path.name, read_flags, dir_fd=parent_fd)
            file_fds[name] = fd
            info = os.fstat(fd)
            expected = sampler_admission_stat(sampler, name)
            if (
                not stat.S_ISREG(info.st_mode)
                or info.st_dev != expected["device"]
                or info.st_ino != expected["inode"]
                or info.st_uid != expected["uid"]
                or info.st_gid != expected["gid"]
                or stat.S_IMODE(info.st_mode) != expected["mode"]
                or info.st_nlink != 1
                or name not in ("raw_csv", "validation_jsonl")
                and stat_receipt(info) != expected
                or name in ("raw_csv", "validation_jsonl")
                and (info.st_size < expected["size"] or info.st_mtime_ns < expected["mtime_ns"])
            ):
                fail("sampler held " + name + " differs from admission")
        raw_prefix = os.pread(
            file_fds["raw_csv"], sampler["csv_bytes"], 0
        )
        if (
            len(raw_prefix) != sampler["csv_bytes"]
            or sha256_bytes(raw_prefix) != sampler["csv_sha256"]
            or not raw_prefix.endswith(b"\n")
        ):
            fail("sampler held raw CSV admission prefix differs")
        validation_info = os.fstat(file_fds["validation_jsonl"])
        validation_prefix = os.pread(
            file_fds["validation_jsonl"], validation_info.st_size, 0
        )
        if len(validation_prefix) != validation_info.st_size:
            fail("sampler held validation prefix is incomplete")
        admission_validation = sampler["validation_jsonl"]
        admission_size = admission_validation["bytes"]
        if (
            validation_info.st_size < admission_size
            or sha256_bytes(validation_prefix[:admission_size])
            != admission_validation["sha256"]
            or not validation_prefix[:admission_size].endswith(b"\n")
        ):
            fail("sampler held validation admission prefix differs")
        validate_validation_stream_header(
            validation_prefix[:admission_size], args,
            "sampler held validation prefix",
        )
        header_size = len(sampler["validation_header_ascii"].encode("ascii"))
        handles: Dict[str, Any] = {
            "file_fds": file_fds,
            "parent_fd": parent_fd,
            "pidfd": -1,
            "pidfd_exit_observed_ns": None,
            "receipt_admission_stat": dict(
                sampler["receipt_file"]["stat"]
            ),
            "validation_buffer": b"",
            "validation_csv_timestamp_ns": [],
            "validation_last_csv_timestamp_ns": None,
            "validation_last_monotonic_ns": None,
            "validation_monotonic_ns": [],
            "validation_offset": validation_info.st_size,
            "validation_sample_index": None,
        }
        initial_event = consume_sampler_validation_bytes(
            handles, validation_prefix[header_size:],
            "sampler held validation admission",
        )
        if initial_event != "none":
            fail("sampler terminal decision preceded worker admission")
        start_before = process_start_ticks(args.sampler_pid)
        security_before = read_process_security(
            args.sampler_pid, "sampler pidfd before"
        )
        pidfd = os.pidfd_open(args.sampler_pid, 0)
        start_after = process_start_ticks(args.sampler_pid)
        security_after = read_process_security(
            args.sampler_pid, "sampler pidfd after"
        )
        if (
            start_before != args.expected_sampler_process_start_ticks
            or start_after != start_before
            or not canonical_equal(security_before, security_after)
            or select.select([pidfd], [], [], 0.0)[0]
        ):
            fail("sampler identity/liveness changed while acquiring pidfd")
        validate_process_security(
            security_after, args.expected_sampler_uid,
            args.expected_sampler_gid, [args.expected_sampler_i2c_gid],
            "sampler pidfd security",
        )
        handles["pidfd"] = pidfd
        return handles
    except BaseException:
        for fd in file_fds.values():
            try:
                os.close(fd)
            except OSError:
                pass
        if pidfd >= 0:
            os.close(pidfd)
        if parent_fd >= 0:
            os.close(parent_fd)
        raise


def close_sampler_admission_handles(handles: Optional[Mapping[str, Any]]) -> None:
    if not handles or handles.get("closed"):
        return
    errors: List[BaseException] = []
    for fd in list(handles.get("file_fds", {}).values()):
        try:
            os.close(fd)
        except BaseException as exc:
            errors.append(exc)
    for name in ("pidfd", "parent_fd"):
        fd = handles.get(name, -1)
        if type(fd) is int and fd >= 0:
            try:
                os.close(fd)
            except BaseException as exc:
                errors.append(exc)
    if type(handles) is dict:
        handles["closed"] = True
        handles["file_fds"] = {}
        handles["pidfd"] = -1
        handles["parent_fd"] = -1
    if errors:
        fail("sampler admission descriptor close failed: " + exception_text(errors[0]))


def poll_sampler_supervision(handles: Mapping[str, Any]) -> str:
    if not handles or handles.get("closed"):
        return "sampler-monitor-invalid:admission-handles-closed"
    receipt_info = os.fstat(handles["file_fds"]["receipt_file"])
    receipt_admission = handles["receipt_admission_stat"]
    if (
        not stat.S_ISREG(receipt_info.st_mode)
        or receipt_info.st_dev != receipt_admission["device"]
        or receipt_info.st_ino != receipt_admission["inode"]
        or receipt_info.st_uid != receipt_admission["uid"]
        or receipt_info.st_gid != receipt_admission["gid"]
        or stat.S_IMODE(receipt_info.st_mode) != 0o444
        or receipt_info.st_nlink != 1
        or receipt_info.st_size < 0
        or receipt_info.st_size > MAX_SAMPLER_RECEIPT_BYTES
    ):
        return "sampler-monitor-invalid:receipt-file-policy"
    if receipt_info.st_size > 0:
        return "sampler-terminal-intent"
    pidfd = handles["pidfd"]
    if select.select([pidfd], [], [], 0.0)[0]:
        if type(handles) is dict and handles.get("pidfd_exit_observed_ns") is None:
            handles["pidfd_exit_observed_ns"] = time.monotonic_ns()
        return "sampler-process-exited"
    fd = handles["file_fds"]["validation_jsonl"]
    info = os.fstat(fd)
    if (
        not stat.S_ISREG(info.st_mode) or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != 0o444
        or info.st_size < handles["validation_offset"]
        or info.st_size > MAX_SAMPLER_VALIDATION_BYTES
    ):
        return "sampler-monitor-invalid:validation-file-policy"
    offset = handles["validation_offset"]
    appended = bytearray()
    while offset < info.st_size:
        block = os.pread(fd, min(65536, info.st_size - offset), offset)
        if not block:
            return "sampler-monitor-invalid:validation-short-read"
        appended.extend(block)
        offset += len(block)
    if type(handles) is dict:
        handles["validation_offset"] = offset
    try:
        event = consume_sampler_validation_bytes(
            handles, bytes(appended), "live sampler validation event"
        )
    except BaseException as exc:
        return "sampler-monitor-invalid:" + exception_text(exc)
    if event != "none":
        return event
    last_monotonic_ns = handles.get("validation_last_monotonic_ns")
    if type(last_monotonic_ns) is not int:
        return "sampler-monitor-invalid:validation-heartbeat-absent"
    if (
        time.monotonic_ns()
        > last_monotonic_ns + SAMPLER_HEARTBEAT_MAX_GAP_NS
    ):
        return "sampler-monitor-invalid:validation-heartbeat-stale"
    return "none"


def held_sampler_snapshot(
    fd: int, max_bytes: int, where: str, *, allowed_nlinks: Sequence[int] = (1,),
) -> Tuple[bytes, os.stat_result]:
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode) or before.st_nlink not in allowed_nlinks
        or before.st_size < 0 or before.st_size > max_bytes
    ):
        fail(where + " held-file policy differs")
    data = bytearray()
    offset = 0
    while offset < before.st_size:
        block = os.pread(fd, min(1024 * 1024, before.st_size - offset), offset)
        if not block:
            fail(where + " held snapshot is short")
        data.extend(block)
        offset += len(block)
    after = os.fstat(fd)
    if not same_file_receipt(before, after):
        fail(where + " changed during held snapshot")
    return bytes(data), after


def collect_sampler_terminal_evidence(
    args: argparse.Namespace, sampler_admission: Mapping[str, Any],
    handles: Mapping[str, Any], start_ns: int,
    child_start_ns: Optional[int], child_reap_ns: Optional[int],
    deadline: float,
) -> Dict[str, Any]:
    exact_int(start_ns, 0, MAX_INT63, "sampler terminal start endpoint")
    if not handles or handles.get("closed"):
        fail("sampler terminal evidence lacks held admission descriptors")
    pidfd = handles["pidfd"]
    while not select.select([pidfd], [], [], 0.0)[0]:
        monitor = poll_sampler_supervision(handles)
        if monitor.startswith("sampler-monitor-invalid:"):
            fail("sampler terminal monitor failed: " + monitor)
        if time.monotonic() >= deadline:
            fail("sampler terminal process did not exit before deadline")
        select.select([pidfd], [], [], min(0.05, deadline - time.monotonic()))
    exit_observed_ns = handles.get("pidfd_exit_observed_ns")
    if exit_observed_ns is None:
        exit_observed_ns = time.monotonic_ns()
        if type(handles) is dict:
            handles["pidfd_exit_observed_ns"] = exit_observed_ns

    file_fds = handles["file_fds"]
    raw_csv, raw_info = held_sampler_snapshot(
        file_fds["raw_csv"], MAX_SAMPLER_CSV_BYTES,
        "sampler terminal raw CSV",
    )
    validation_jsonl, validation_info = held_sampler_snapshot(
        file_fds["validation_jsonl"], MAX_SAMPLER_VALIDATION_BYTES,
        "sampler terminal validation JSONL",
    )
    producer_receipt_ascii, receipt_info = held_sampler_snapshot(
        file_fds["receipt_file"], MAX_SAMPLER_RECEIPT_BYTES,
        "sampler terminal producer receipt",
    )
    pid_data, pid_info = held_sampler_snapshot(
        file_fds["pid_file"], 64, "sampler terminal held PID file",
        allowed_nlinks=(0,),
    )
    if (
        not producer_receipt_ascii.endswith(b"\n")
        or pid_data != (str(args.sampler_pid) + "\n").encode("ascii")
        or stat.S_IMODE(raw_info.st_mode) != 0o444
        or stat.S_IMODE(validation_info.st_mode) != 0o444
        or stat.S_IMODE(receipt_info.st_mode) != 0o444
        or stat.S_IMODE(pid_info.st_mode) != 0o444
        or any(
            info.st_uid != args.expected_sampler_uid
            or info.st_gid != args.expected_sampler_gid
            for info in (raw_info, validation_info, receipt_info, pid_info)
        )
    ):
        fail("sampler terminal held artifact mode/owner/content differs")

    parent_fd = handles["parent_fd"]
    parent_info = os.fstat(parent_fd)
    admission_parent = sampler_admission["evidence_parent"]["stat"]
    if (
        not stat.S_ISDIR(parent_info.st_mode)
        or parent_info.st_dev != admission_parent["device"]
        or parent_info.st_ino != admission_parent["inode"]
        or parent_info.st_uid != admission_parent["uid"]
        or parent_info.st_gid != admission_parent["gid"]
        or stat.S_IMODE(parent_info.st_mode) != 0o700
    ):
        fail("sampler terminal held evidence parent differs")
    expected_named = {
        args.sampler_csv.name: raw_info,
        args.sampler_validation_jsonl.name: validation_info,
        args.sampler_receipt.name: receipt_info,
    }
    observed_names: Dict[str, os.stat_result] = {}
    with os.scandir(parent_fd) as entries:
        for entry in entries:
            if entry.name in observed_names or entry.name not in expected_named:
                fail("sampler terminal evidence parent roster differs")
            if not entry.is_file(follow_symlinks=False):
                fail("sampler terminal evidence contains a non-regular name")
            observed_names[entry.name] = entry.stat(follow_symlinks=False)
    if set(observed_names) != set(expected_named):
        fail("sampler terminal evidence parent roster differs")
    for name, expected_info in expected_named.items():
        if not same_file_receipt(observed_names[name], expected_info):
            fail("sampler terminal named artifact differs from held inode")
    try:
        os.stat(args.sampler_pid_file.name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        fail("sampler terminal PID pathname still exists")

    producer_receipt_preview = parse_canonical_json_line(
        producer_receipt_ascii, "sampler terminal producer receipt preview"
    )
    stream_parts = split_terminal_streams(
        raw_csv, validation_jsonl,
        allow_unpaired=producer_receipt_preview.get("outcome") == "sampler_error",
    )
    physical, timestamps = sampler_rows(stream_parts["paired_raw"])
    start_lines = [index for index, value in timestamps.items() if value == start_ns]
    if len(start_lines) != 1:
        fail("sampler terminal start endpoint is not unique")
    end_ns = max(timestamps.values()) if timestamps else None
    if end_ns is None or end_ns < start_ns:
        fail("sampler terminal lacks an end endpoint")
    end_lines = [index for index, value in timestamps.items() if value == end_ns]
    if len(end_lines) != 1 or end_lines[0] < start_lines[0]:
        fail("sampler terminal end endpoint is not unique")
    window_csv = physical[0] + b"".join(
        physical[start_lines[0]:end_lines[0] + 1]
    )
    window_validation = select_validation_window(
        stream_parts["validation_complete"], window_csv,
        "sampler terminal retained window",
    )
    covers = bool(
        child_start_ns is not None and child_reap_ns is not None
        and start_ns <= child_start_ns <= child_reap_ns <= end_ns
    )
    terminal = {
        "admission_receipt_sha256": sha256_bytes(
            canonical_bytes(sampler_admission)
        ),
        "coverage": {
            "child_reap_monotonic_ns": child_reap_ns,
            "child_start_monotonic_ns": child_start_ns,
            "coverage_shortfall_ns": (
                None if child_reap_ns is None
                else max(0, child_reap_ns - end_ns)
            ),
            "covers_child_interval": covers,
            "window_end_monotonic_ns": end_ns,
            "window_start_monotonic_ns": start_ns,
        },
        "pid_file_held_binding": binding_from_snapshot(pid_data, pid_info),
        "process_exit_observation": "linux-pidfd-readable-nonchild-v1",
        "process_exit_observed_monotonic_ns": exit_observed_ns,
        "producer_receipt_ascii": producer_receipt_ascii.decode("ascii"),
        "producer_receipt_binding": binding_from_snapshot(
            producer_receipt_ascii, receipt_info
        ),
        "raw_csv_ascii": raw_csv.decode("ascii"),
        "raw_csv_binding": binding_from_snapshot(raw_csv, raw_info),
        "schema": SAMPLER_TERMINAL_SCHEMA,
        "stream_suffixes": stream_parts["suffix_receipt"],
        "terminal_status": "invalid",
        "validation_jsonl_ascii": validation_jsonl.decode("ascii"),
        "validation_jsonl_binding": binding_from_snapshot(
            validation_jsonl, validation_info
        ),
        "window_csv_ascii": window_csv.decode("ascii"),
        "window_validation_jsonl_ascii": window_validation.decode("ascii"),
    }
    validate_sampler_terminal_receipt(
        terminal, sampler_admission, child_start_ns, child_reap_ns
    )
    return terminal


def canonical_counter(text: str, where: str) -> int:
    if not text.isdecimal() or (len(text) > 1 and text.startswith("0")):
        fail(f"{where} is not a canonical nonnegative integer")
    return exact_int(int(text), 0, MAX_UINT64, where)


def finite_scalar(text: str, where: str) -> float:
    try:
        value = float(text)
    except ValueError:
        fail(f"{where} is not numeric")
    if not math.isfinite(value):
        fail(f"{where} is not finite")
    return value


def summarize_thermal_window_bytes(
    window_bytes: bytes, start_ns: int, end_ns: int,
) -> Dict[str, Any]:
    physical = window_bytes.splitlines(keepends=True)
    if (
        not physical
        or tuple(parse_csv_physical_line(physical[0], "thermal window header"))
        != THERMAL_HEADER
        or len(physical) < 2
    ):
        fail("thermal window framing/header differs")
    parse_failures: List[str] = []
    valid_count = 0
    invalid_count = 0
    prior_ns: Optional[int] = None
    observed_timestamps: List[int] = []
    cpu_values: List[float] = []
    dimm_values: List[float] = []
    dimm_error_values: List[int] = []
    edac_ce_values: List[int] = []
    edac_ue_values: List[int] = []
    for offset, raw in enumerate(physical[1:]):
        label = f"thermal window row {offset + 2}"
        row_valid = True
        try:
            row = parse_csv_physical_line(raw, label)
            if len(row) != len(THERMAL_HEADER):
                fail(f"{label} width differs")
            mono = finite_scalar(row[1], label + " monotonic")
            mono_ns = int(round(mono * 1_000_000_000.0))
            if not start_ns <= mono_ns <= end_ns:
                fail(f"{label} timestamp is outside the window")
            if prior_ns is not None and (
                mono_ns <= prior_ns
                or mono_ns - prior_ns > SAMPLER_HEARTBEAT_MAX_GAP_NS
            ):
                fail(f"{label} timestamp order/gap differs")
            prior_ns = mono_ns
            observed_timestamps.append(mono_ns)
            busy = finite_scalar(row[2], label + " busy")
            mhz = finite_scalar(row[3], label + " MHz")
            cpu = finite_scalar(row[4], label + " CPU temperature")
            dimms = [
                finite_scalar(value, label + f" DIMM {index}")
                for index, value in enumerate(row[5:13])
            ]
            dimm_errors = canonical_counter(row[13], label + " DIMM errors")
            loads = [
                finite_scalar(value, label + " load")
                for value in row[14:17]
            ]
            edac_ce = canonical_counter(row[17], label + " EDAC CE")
            edac_ue = canonical_counter(row[18], label + " EDAC UE")
            cpu_values.append(cpu)
            dimm_values.extend(dimms)
            dimm_error_values.append(dimm_errors)
            edac_ce_values.append(edac_ce)
            edac_ue_values.append(edac_ue)
            if (
                not 0.0 <= busy <= 100.0 or mhz <= 0.0
                or not 0.0 <= cpu <= 120.0
                or any(not 0.0 < value <= 100.0 for value in dimms)
                or any(value < 0.0 for value in loads)
            ):
                row_valid = False
                append_unique_failure(
                    parse_failures, label + " range violation"
                )
        except (ValidationError, ValueError, OverflowError) as exc:
            row_valid = False
            append_unique_failure(parse_failures, exception_text(exc))
        if row_valid:
            valid_count += 1
        else:
            invalid_count += 1
    sample_count = len(physical) - 1
    if sample_count < 2:
        append_unique_failure(
            parse_failures, "thermal window has fewer than two samples"
        )
    if (
        not observed_timestamps
        or observed_timestamps[0] != start_ns
        or observed_timestamps[-1] != end_ns
    ):
        append_unique_failure(
            parse_failures, "thermal window endpoint timestamps differ"
        )

    def millic(values: Sequence[float], name: str) -> Optional[int]:
        if not values:
            return None
        converted = int(round(max(values) * 1000.0))
        if not -(1 << 63) <= converted <= MAX_INT63:
            append_unique_failure(
                parse_failures, name + " maximum exceeds signed int63"
            )
            return None
        return converted

    return {
        "cpu_tctl_max_millic": millic(cpu_values, "CPU temperature"),
        "dimm_max_millic": millic(dimm_values, "DIMM temperature"),
        "dimm_read_errors": (
            max(dimm_error_values) if dimm_error_values else None
        ),
        "edac_ce_max": max(edac_ce_values) if edac_ce_values else None,
        "edac_ue_max": max(edac_ue_values) if edac_ue_values else None,
        "invalid_sample_count": invalid_count,
        "parse_failures": parse_failures,
        "sample_count": sample_count,
        "valid_sample_count": valid_count,
    }


def thermal_policy_violation_messages(
    thermal: Mapping[str, Any],
) -> List[str]:
    """Reconstruct the producer's exact ordered thermal-policy markers."""
    violations: List[str] = []
    for message in thermal["parse_failures"]:
        append_unique_failure(violations, "thermal parse: " + message)
    for message in thermal["validation_failures"]:
        append_unique_failure(violations, "thermal validation: " + message)
    if (
        thermal["cpu_tctl_max_millic"] is None
        or thermal["cpu_tctl_max_millic"] > THERMAL_MAX_MILLIC
    ):
        append_unique_failure(
            violations, "CPU thermal cap violated or unavailable"
        )
    if (
        thermal["dimm_max_millic"] is None
        or thermal["dimm_max_millic"] > THERMAL_MAX_MILLIC
    ):
        append_unique_failure(
            violations, "DIMM thermal cap violated or unavailable"
        )
    if (
        thermal["dimm_read_errors"] is None
        or thermal["dimm_read_errors"] != 0
    ):
        append_unique_failure(
            violations,
            "DIMM read-error policy violated or unavailable",
        )
    if thermal["edac_ce_max"] is None or thermal["edac_ce_max"] != 0:
        append_unique_failure(
            violations, "EDAC CE policy violated or unavailable"
        )
    if thermal["edac_ue_max"] is None or thermal["edac_ue_max"] != 0:
        append_unique_failure(
            violations, "EDAC UE policy violated or unavailable"
        )
    return violations


def tolerant_thermal_window(
    args: argparse.Namespace, sampler: Mapping[str, Any], start_ns: int,
    end_ns: int, *, allow_stale_suffixes: bool = False,
    deadline: Optional[float] = None,
    stale_stream_capture: Optional[Mapping[str, Any]] = None,
) -> Tuple[Dict[str, Any], List[str], List[str]]:
    collection: List[str] = []
    violations: List[str] = []
    try:
        if deadline is not None and time.monotonic() >= deadline:
            fail("thermal collection reached its recovery deadline")
        if allow_stale_suffixes:
            if not stale_stream_capture:
                fail("stale thermal evidence lacks its cutoff snapshot")
            raw_data = stale_stream_capture["csv_full_data"]
            info = stale_stream_capture["csv_info"]
            if stat_receipt(info) != sampler["csv_stat"]:
                fail("stale thermal raw CSV differs from attestation")
            complete_end = raw_data.rfind(b"\n") + 1
            if complete_end == 0:
                fail("stale thermal raw CSV has no complete line")
            data = raw_data[:complete_end]
        else:
            data, info = sampler_snapshot(args.sampler_csv)
        if deadline is not None and time.monotonic() >= deadline:
            fail("thermal raw snapshot reached its recovery deadline")
        if (
            info.st_dev != sampler["csv_device"]
            or info.st_ino != sampler["csv_inode"]
        ):
            fail("thermal CSV identity changed after sampler attestation")
        physical, timestamps = sampler_rows(data)
        start_lines = [index for index, value in timestamps.items() if value == start_ns]
        end_lines = [index for index, value in timestamps.items() if value == end_ns]
        if len(start_lines) != 1 or len(end_lines) != 1:
            fail("thermal window endpoints are not unique exact CSV samples")
        first = start_lines[0]
        last = end_lines[0]
        if first > last:
            fail("thermal window endpoints are reversed")
        window_bytes = physical[0] + b"".join(physical[first:last + 1])
        if len(window_bytes) > MAX_THERMAL_WINDOW_BYTES:
            fail("thermal window exceeds its evidence bound")
        try:
            window_ascii = window_bytes.decode("ascii")
        except UnicodeDecodeError:
            fail("thermal window is not ASCII")
        if allow_stale_suffixes:
            validation_data = stale_stream_capture["validation_data"]
            validation_info = stale_stream_capture["validation_info"]
        else:
            validation_data, validation_info = stable_regular_snapshot(
                args.sampler_validation_jsonl, MAX_SAMPLER_VALIDATION_BYTES,
                "thermal validation JSONL",
            )
        if deadline is not None and time.monotonic() >= deadline:
            fail("thermal validation snapshot reached its recovery deadline")
        validation_artifact = sampler["validation_jsonl"]
        if (
            validation_info.st_dev != validation_artifact["stat"]["device"]
            or validation_info.st_ino != validation_artifact["stat"]["inode"]
            or validation_info.st_uid != sampler["process_uid"]
            or validation_info.st_gid != sampler["process_gid"]
            or validation_info.st_nlink != 1
            or stat.S_IMODE(validation_info.st_mode) != 0o444
        ):
            fail("thermal validation JSONL identity changed after attestation")
        if allow_stale_suffixes:
            if (
                len(validation_data) != validation_artifact["bytes"]
                or sha256_bytes(validation_data)
                != validation_artifact["sha256"]
                or stat_receipt(validation_info)
                != validation_artifact["stat"]
            ):
                fail("stale thermal validation differs from attestation")
            stale_parts = validate_stale_sampler_stream_prefix(
                raw_data, validation_data, args,
                "stale thermal stream prefix",
                stale_stream_capture["snapshot_observed_monotonic_ns"],
            )
            if (
                sampler["csv_bytes"] != len(stale_parts["paired_raw"])
                or sampler["csv_sha256"]
                != sha256_bytes(stale_parts["paired_raw"])
                or stale_parts["last_csv_timestamp_ns"] != end_ns
            ):
                fail("stale thermal paired prefix differs from attestation")
            validation_data = stale_parts["validation_complete"]
        validation_window = select_validation_window(
            validation_data, window_bytes, "thermal validation JSONL"
        )
        try:
            validation_ascii = validation_window.decode("ascii")
        except UnicodeDecodeError:
            fail("thermal validation window is not ASCII")
        summary = summarize_thermal_window_bytes(window_bytes, start_ns, end_ns)
        validation_summary = summarize_validation_window_bytes(
            validation_window, window_bytes, sampler
        )
        if deadline is not None and time.monotonic() >= deadline:
            fail("thermal replay reached its recovery deadline")
        parse_failures = summary["parse_failures"]
        sample_count = summary["sample_count"]
        valid_count = summary["valid_sample_count"]
        invalid_count = summary["invalid_sample_count"]
        cpu_max = summary["cpu_tctl_max_millic"]
        dimm_max = summary["dimm_max_millic"]
        dimm_errors_max = summary["dimm_read_errors"]
        edac_ce_max = summary["edac_ce_max"]
        edac_ue_max = summary["edac_ue_max"]
        thermal = {
            "cpu": EXPECTED_SAMPLER_CPU,
            "cpu_tctl_max_millic": cpu_max,
            "csv_device": info.st_dev,
            "csv_inode": info.st_ino,
            "csv_path": str(args.sampler_csv),
            "dimm_max_millic": dimm_max,
            "dimm_read_errors": dimm_errors_max,
            "edac_ce_max": edac_ce_max,
            "edac_ue_max": edac_ue_max,
            "invalid_sample_count": invalid_count,
            "parse_failures": parse_failures,
            "pid": args.sampler_pid,
            "process_start_ticks": sampler["process_start_ticks"],
            "sample_count": sample_count,
            "script_path": str(args.sampler_script),
            "script_sha256": sampler["script_sha256"],
            "terminal_status": "complete",
            "valid_sample_count": valid_count,
            "validation_attempt_errors_total": validation_summary[
                "attempt_errors_total"
            ],
            "validation_device": validation_info.st_dev,
            "validation_failures": validation_summary["failures"],
            "validation_inode": validation_info.st_ino,
            "validation_jsonl_ascii": validation_ascii,
            "validation_jsonl_bytes": len(validation_window),
            "validation_jsonl_sha256": sha256_bytes(validation_window),
            "validation_path": str(args.sampler_validation_jsonl),
            "validation_sample_count": validation_summary["sample_count"],
            "window_csv_ascii": window_ascii,
            "window_csv_bytes": len(window_bytes),
            "window_csv_sha256": sha256_bytes(window_bytes),
            "window_end_monotonic_ns": end_ns,
            "window_start_monotonic_ns": start_ns,
        }
        violations = thermal_policy_violation_messages(thermal)
        return thermal, collection, violations
    except BaseException as exc:
        append_unique_failure(collection, "thermal collection: " + exception_text(exc))
        return {}, collection, violations


def unbounded_sibling_tick_receipt(
    start: Mapping[str, Any], end: Mapping[str, Any], interval_start_ns: int,
    interval_end_ns: int,
) -> Dict[str, Any]:
    validate_cpu_tick_receipt(start, "sibling start")
    validate_cpu_tick_receipt(end, "sibling end")
    if start["cpu"] != 56 or end["cpu"] != 56:
        fail("sibling endpoint CPU differs")
    if not (
        start["read_monotonic_ns"] <= interval_start_ns
        <= interval_end_ns <= end["read_monotonic_ns"]
    ):
        fail("sibling endpoints do not cover the interval")
    start_fields = start["tick_fields"]
    end_fields = end["tick_fields"]
    if any(end_fields[name] < start_fields[name] for name in CPU_TICK_FIELD_KEYS):
        fail("sibling tick fields are not monotone")
    delta = end["non_idle_ticks"] - start["non_idle_ticks"]
    exact_int(delta, 0, MAX_UINT64, "sibling non-idle delta")
    return {
        "cpu": 56,
        "delta_non_idle_ticks": delta,
        "end": end,
        "interval_end_monotonic_ns": interval_end_ns,
        "interval_start_monotonic_ns": interval_start_ns,
        "start": start,
    }


def prepare_health(
    args: argparse.Namespace, target_cpu: int, deadline: float,
    health_modules: SealedHealthModules,
    sampler_admission: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    receipt = empty_health_receipt(args)
    collection: List[str] = []
    violations: List[str] = []
    state: Dict[str, Any] = {
        "collection_failures": collection,
        "deadline": deadline,
        "health_modules": health_modules,
        "receipt": receipt,
        "ready": False,
        "violations": violations,
    }
    health_api = health_modules.native
    try:
        exact_int(target_cpu, EXPECTED_TARGET_CPU, EXPECTED_TARGET_CPU, "target CPU")
        exact_int(
            args.controller_cpu, EXPECTED_CONTROLLER_CPU,
            EXPECTED_CONTROLLER_CPU, "controller CPU",
        )
        exact_int(
            args.sampler_cpu, EXPECTED_SAMPLER_CPU,
            EXPECTED_SAMPLER_CPU, "sampler CPU",
        )
        allowed = list(getattr(args, "controller_initial_affinity", []))
        receipt["controller_initial_affinity"] = allowed
        if any(
            cpu not in allowed for cpu in (
                EXPECTED_TARGET_CPU, EXPECTED_CONTROLLER_CPU,
                EXPECTED_SAMPLER_CPU,
            )
        ):
            fail("frozen CPUs 120/121/122 are outside initial affinity")
        if os.sched_getaffinity(0) != {EXPECTED_CONTROLLER_CPU}:
            fail("controller did not retain singleton CPU121 affinity")
        target_core = health_api._cpu_physical_core(EXPECTED_TARGET_CPU)
        controller_core = health_api._cpu_physical_core(EXPECTED_CONTROLLER_CPU)
        sampler_core = health_api._cpu_physical_core(EXPECTED_SAMPLER_CPU)
        siblings = tuple(health_api.timing_sibling_cpus([EXPECTED_TARGET_CPU]))
        threads = tuple(sorted((EXPECTED_TARGET_CPU, *siblings)))
        if (
            target_core != EXPECTED_TARGET_CORE
            or controller_core != EXPECTED_CONTROLLER_CORE
            or sampler_core != EXPECTED_SAMPLER_CORE
            or threads != EXPECTED_TARGET_THREADS
            or siblings != (56,)
        ):
            fail("frozen 120/121/122 topology differs")
        receipt["target_core"] = list(target_core)
        receipt["controller_core"] = list(controller_core)
        receipt["sampler_core"] = list(sampler_core)
        receipt["target_threads"] = list(threads)
        health_modules.runner._preflight_sampler(
            args.sampler_pid, EXPECTED_SAMPLER_CPU,
            args.sampler_script, args.sampler_csv,
        )
        if process_start_ticks(args.sampler_pid) != args.expected_sampler_process_start_ticks:
            fail("sampler process start identity differs at admission")
        admission_receipt = dict(
            sampler_admission
            if sampler_admission is not None
            else make_sampler_attestation(
                args, 0, 0, health_modules, deadline
            )
        )
        validate_sampler_receipt(admission_receipt)
        receipt["sampler_admission"] = admission_receipt
        receipt["sampler_admission_receipt_sha256"] = sha256_bytes(
            canonical_bytes(admission_receipt)
        )
        state["sampler_admission"] = admission_receipt
        state["sampler_handles"] = open_sampler_admission_handles(
            args, admission_receipt
        )
        initial_data, initial_info = sampler_snapshot(args.sampler_csv)
        exact_int(
            initial_info.st_dev, args.expected_sampler_csv_device,
            args.expected_sampler_csv_device, "sampler CSV device",
        )
        exact_int(
            initial_info.st_ino, args.expected_sampler_csv_inode,
            args.expected_sampler_csv_inode, "sampler CSV inode",
        )
        _, initial_timestamps = sampler_rows(initial_data)
        if not initial_timestamps:
            fail("sampler CSV has no timestamped baseline")
        admission_start = read_cpu_tick_receipt(56)
        admission_start_ns = admission_start["read_monotonic_ns"]
        admission_sample_deadline = min(deadline, time.monotonic() + 6.0)
        sample_start_ns = wait_for_sampler_timestamp(
            args.sampler_csv, admission_sample_deadline,
            greater_than=max(initial_timestamps.values()),
        )
        admission_monitor = wait_for_sampler_validation_timestamp_supervised(
            state["sampler_handles"], sample_start_ns,
            admission_sample_deadline,
        )
        if admission_monitor != "none":
            fail(
                "sampler terminal state preceded claim: "
                + admission_monitor
            )
        admission_end = read_cpu_tick_receipt(56)
        admission_end_ns = admission_end["read_monotonic_ns"]
        admission = unbounded_sibling_tick_receipt(
            admission_start, admission_end, admission_start_ns, admission_end_ns
        )
        receipt["admission_sibling_ticks"] = admission
        # Admission is preserved as a diagnostic receipt.  The frozen <=5
        # acceptance cap applies only to the worker interval below.
        state["sample_start_ns"] = sample_start_ns
        state["sibling_start"] = read_cpu_tick_receipt(56)
        state["ready"] = not collection and not violations
    except BaseException as exc:
        append_unique_failure(collection, "health admission: " + exception_text(exc))
        try:
            close_sampler_admission_handles(state.get("sampler_handles"))
        except BaseException as close_exc:
            append_unique_failure(
                collection, "health admission close: " + exception_text(close_exc)
            )
    finalize_health_receipt(receipt, collection, violations)
    return state


def finish_health(
    args: argparse.Namespace, target_cpu: int, child_start_ns: Optional[int],
    child_reap_ns: Optional[int], state: Mapping[str, Any], deadline: float,
    health_modules: SealedHealthModules,
) -> Tuple[Dict[str, Any], List[str]]:
    receipt = dict(state.get("receipt", empty_health_receipt(args)))
    collection = list(state.get("collection_failures", []))
    violations = list(state.get("violations", []))
    if child_start_ns is not None:
        receipt["child_start_monotonic_ns"] = child_start_ns
    if child_reap_ns is not None:
        receipt["child_reap_monotonic_ns"] = child_reap_ns
    try:
        affinity_end = sorted(os.sched_getaffinity(0))
        receipt["controller_singleton_affinity_end"] = affinity_end
        if affinity_end != [EXPECTED_CONTROLLER_CPU]:
            append_unique_failure(violations, "controller CPU121 affinity changed")
    except BaseException as exc:
        append_unique_failure(collection, "controller affinity end: " + exception_text(exc))
    if child_start_ns is not None and child_reap_ns is not None and state.get("sibling_start"):
        try:
            sibling_end = read_cpu_tick_receipt(56)
            during = unbounded_sibling_tick_receipt(
                state["sibling_start"], sibling_end, child_start_ns, child_reap_ns
            )
            receipt["sibling_ticks"] = [during]
            if during["delta_non_idle_ticks"] > SIBLING_NON_IDLE_TICK_CAP:
                append_unique_failure(violations, "CPU56 run non-idle tick cap exceeded")
        except BaseException as exc:
            append_unique_failure(collection, "run sibling ticks: " + exception_text(exc))
    else:
        append_unique_failure(collection, "run sibling ticks unavailable")

    start_ns = state.get("sample_start_ns")
    endpoint_target = child_reap_ns if child_reap_ns is not None else time.monotonic_ns()
    sampler_handles = state.get("sampler_handles")
    supervision_event = state.get("sampler_supervision_event", "none")

    def retain_terminal_event(event: str) -> None:
        append_unique_failure(collection, "sampler endpoint: " + event)
        if event in {
            "sampler-monitor-invalid:validation-heartbeat-stale",
            "sampler-monitor-invalid:coverage-recovery-cutoff",
        }:
            # The worker has already been killed by supervision.  Preserve
            # the complete paired prefix from the still-bound live sampler,
            # but do not wait for terminal artifacts from a process that may
            # remain stuck indefinitely.
            try:
                stale_stream_capture: Dict[str, Any] = {}
                recovery_deadline = min(
                    time.monotonic()
                    + SAMPLER_STALE_RECOVERY_WORK_SECONDS,
                    deadline - SAMPLER_HEALTH_FINALIZE_RESERVE_SECONDS,
                )
                if time.monotonic() >= recovery_deadline:
                    fail("sampler stale recovery lacks its reserved margin")
                sampler = make_sampler_attestation(
                    args, start_ns, start_ns, health_modules,
                    recovery_deadline,
                    allow_stale_evidence=True,
                    stale_stream_capture=stale_stream_capture,
                )
                stale_end_ns = sampler["window_end_monotonic_ns"]
                thermal, thermal_collection, thermal_violations = (
                    tolerant_thermal_window(
                        args, sampler, start_ns, stale_end_ns,
                        allow_stale_suffixes=True,
                        deadline=recovery_deadline,
                        stale_stream_capture=stale_stream_capture,
                    )
                )
                if time.monotonic() >= recovery_deadline:
                    fail("sampler stale recovery reached its strict deadline")
                next_collection = list(collection)
                for message in thermal_collection:
                    append_unique_failure(next_collection, message)
                next_violations = list(violations)
                for message in thermal_violations:
                    append_unique_failure(next_violations, message)
                receipt["sampler"] = sampler
                receipt["sampler_receipt_sha256"] = sha256_bytes(
                    canonical_bytes(sampler)
                )
                receipt["sampler_terminal"] = {}
                receipt["sampler_terminal_receipt_sha256"] = None
                receipt["thermal"] = thermal
                collection[:] = next_collection
                violations[:] = next_violations
            except BaseException as exc:
                append_unique_failure(
                    collection,
                    "sampler stale evidence: " + exception_text(exc),
                )
            return
        try:
            terminal = collect_sampler_terminal_evidence(
                args, receipt["sampler_admission"], sampler_handles,
                start_ns, child_start_ns, child_reap_ns, deadline,
            )
            receipt["sampler"] = {}
            receipt["sampler_receipt_sha256"] = None
            receipt["sampler_terminal"] = terminal
            receipt["sampler_terminal_receipt_sha256"] = sha256_bytes(
                canonical_bytes(terminal)
            )
            receipt["thermal"] = {}
        except BaseException as exc:
            append_unique_failure(
                collection,
                "sampler terminal evidence: " + exception_text(exc),
            )

    if type(start_ns) is int:
        if sampler_handles and supervision_event == "none":
            supervision_event = poll_sampler_supervision(sampler_handles)
        if supervision_event != "none":
            retain_terminal_event(supervision_event)
        else:
            coverage_recovery_cutoff = deadline - (
                SAMPLER_STALE_RECOVERY_WORK_SECONDS
                + SAMPLER_HEALTH_FINALIZE_RESERVE_SECONDS
            )
            if time.monotonic() >= coverage_recovery_cutoff:
                retain_terminal_event(
                    "sampler-monitor-invalid:coverage-recovery-cutoff"
                )
                supervision_event = (
                    "sampler-monitor-invalid:coverage-recovery-cutoff"
                )
        if supervision_event == "none":
            try:
                endpoint_deadline = min(
                    coverage_recovery_cutoff, time.monotonic() + 6.0
                )
                sample_end_ns, wait_event = wait_for_sampler_timestamp_supervised(
                    args.sampler_csv,
                    endpoint_deadline,
                    sampler_handles,
                    at_or_after=max(endpoint_target, start_ns + 1),
                )
                if wait_event != "none":
                    retain_terminal_event(wait_event)
                elif sample_end_ns is None:
                    fail("supervised sampler endpoint lacks a result")
                else:
                    validation_event = (
                        wait_for_sampler_validation_timestamp_supervised(
                            sampler_handles, sample_end_ns,
                            endpoint_deadline,
                        )
                    )
                    if validation_event != "none":
                        retain_terminal_event(validation_event)
                    else:
                        sampler = make_sampler_attestation(
                            args, start_ns, sample_end_ns,
                            health_modules, coverage_recovery_cutoff,
                        )
                        sampler_digest = sha256_bytes(canonical_bytes(sampler))
                        thermal, thermal_collection, thermal_violations = (
                            tolerant_thermal_window(
                                args, sampler, start_ns, sample_end_ns,
                                deadline=coverage_recovery_cutoff,
                            )
                        )
                        if time.monotonic() >= coverage_recovery_cutoff:
                            fail(
                                "sampler endpoint collection reached its "
                                "recovery cutoff"
                            )
                        next_collection = list(collection)
                        for message in thermal_collection:
                            append_unique_failure(next_collection, message)
                        next_violations = list(violations)
                        for message in thermal_violations:
                            append_unique_failure(next_violations, message)
                        next_receipt = dict(receipt)
                        next_receipt["sampler"] = sampler
                        next_receipt["sampler_receipt_sha256"] = sampler_digest
                        next_receipt["thermal"] = thermal
                        receipt, collection, violations = (
                            next_receipt, next_collection, next_violations
                        )
            except BaseException as exc:
                terminal_event = (
                    poll_sampler_supervision(sampler_handles)
                    if sampler_handles else "none"
                )
                if terminal_event != "none":
                    retain_terminal_event(terminal_event)
                elif time.monotonic() >= coverage_recovery_cutoff:
                    retain_terminal_event(
                        "sampler-monitor-invalid:coverage-recovery-cutoff"
                    )
                else:
                    append_unique_failure(
                        collection,
                        "sampler endpoint: " + exception_text(exc),
                    )
    else:
        append_unique_failure(collection, "sampler start endpoint unavailable")
    if time.monotonic() >= deadline:
        append_unique_failure(
            collection,
            "health collection reached the strict 115-second deadline",
        )
    try:
        close_sampler_admission_handles(sampler_handles)
    except BaseException as exc:
        append_unique_failure(
            collection,
            "health terminal collection: " + exception_text(exc),
        )
    finalized = finalize_health_receipt(receipt, collection, violations)
    errors = ["health collection: " + item for item in collection]
    errors.extend("health violation: " + item for item in violations)
    return finalized, errors


@dataclass
class ExecutionResult:
    raw: bytes = b""
    stderr: bytes = b""
    returncode: Optional[int] = None
    elapsed_seconds: float = 0.0
    binary: Dict[str, Any] = field(default_factory=dict)
    health: Dict[str, Any] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)

    def raw_complete(self) -> bool:
        return bool(
            self.binary.get("process_started")
            and self.returncode == 0
            and not self.binary.get("timed_out")
            and not self.binary.get("output_overflow")
            and self.binary.get("capture_error") == "none"
        )


def empty_binary_receipt(binary: Path, expected_sha256: str) -> Dict[str, Any]:
    return {
        "capture_error": "none",
        "child_affinity_after_spawn": [],
        "child_pid": None,
        "child_process_start_ticks": None,
        "child_reap_monotonic_ns": None,
        "child_start_monotonic_ns": None,
        "execution_finished_monotonic_ns": None,
        "execution_started_monotonic_ns": None,
        "expected_sha256": expected_sha256,
        "output_overflow": False,
        "path": str(binary),
        "path_stat_after": {},
        "process_started": False,
        "sha256_after": None,
        "sha256_before": None,
        "stat_after": {},
        "stat_before": {},
        "timed_out": False,
    }


def pin_worker_before_exec() -> None:
    os.sched_setaffinity(0, {EXPECTED_TARGET_CPU})


def reap_reached_controller_deadline(
    child_reap_monotonic_ns: Optional[int], controller_started_monotonic_ns: int,
) -> bool:
    """Use one inclusive integer-nanosecond law for the child 110 s bound."""
    if child_reap_monotonic_ns is None:
        return True
    return monotonic_deadline_reached(
        child_reap_monotonic_ns, controller_started_monotonic_ns,
        CONTROLLER_DEADLINE_SECONDS,
    )


def monotonic_deadline_reached(
    observed_ns: int, origin_ns: int, seconds: float,
) -> bool:
    return observed_ns >= origin_ns + int(seconds * 1_000_000_000)


def run_binary(
    binary: Path, cpu: int, expected_sha256: str,
    controller_started: Optional[float] = None,
    controller_started_ns: Optional[int] = None,
    health_args: Optional[argparse.Namespace] = None,
    health_modules: Optional[SealedHealthModules] = None,
    prepared_health_state: Optional[Dict[str, Any]] = None,
    signal_guard: Optional[SignalGuard] = None,
) -> ExecutionResult:
    require_preflight_side_effects_allowed("worker launch")
    if (controller_started is None) != (controller_started_ns is None):
        fail("controller float/integer start origins must be supplied together")
    if controller_started is None and controller_started_ns is None:
        execution_started_ns = time.monotonic_ns()
        controller_started_ns = execution_started_ns
        controller_started = controller_started_ns / 1_000_000_000.0
    else:
        exact_int(
            controller_started_ns, 0, MAX_INT63,
            "controller integer start origin",
        )
        if (
            type(controller_started) is not float
            or not math.isfinite(controller_started)
            or controller_started
            != controller_started_ns / 1_000_000_000.0
        ):
            fail("controller float/integer start origins differ")
        execution_started_ns = time.monotonic_ns()
    initial_health = (
        dict(prepared_health_state["receipt"])
        if prepared_health_state is not None
        else (empty_health_receipt(health_args) if health_args is not None else {})
    )
    result = ExecutionResult(
        binary=empty_binary_receipt(binary, expected_sha256),
        health=initial_health,
    )
    result.binary["execution_started_monotonic_ns"] = execution_started_ns
    child_deadline = controller_started + CONTROLLER_DEADLINE_SECONDS
    health_deadline = controller_started + HEALTH_COLLECTION_DEADLINE_SECONDS
    finalization_deadline = controller_started + 116.0
    fd = -1
    process: Optional[subprocess.Popen] = None
    fd_before: Optional[os.stat_result] = None
    digest_before: Optional[str] = None
    health_state: Dict[str, Any] = prepared_health_state or {}
    child_start_ns: Optional[int] = None
    capture_started = False
    try:
        path_before = os.stat(str(binary), follow_symlinks=False)
        if (
            not stat.S_ISREG(path_before.st_mode)
            or path_before.st_nlink != 1
            or path_before.st_size < 1
            or path_before.st_size > MAX_BINARY_BYTES
            or path_before.st_uid != (
                health_args.expected_binary_uid
                if health_args is not None else os.getuid()
            )
            or path_before.st_gid != (
                health_args.expected_controller_gid
                if health_args is not None else os.getgid()
            )
        ):
            fail("screen binary must be a regular single-link file")
        if stat.S_IMODE(path_before.st_mode) != 0o555:
            fail("screen binary must have exact mode 0555")
        fd = os.open(
            str(binary), nonblocking_read_flags("screen binary execution")
        )
        fd_before = os.fstat(fd)
        result.binary["stat_before"] = stat_receipt(fd_before)
        if not same_file_receipt(path_before, fd_before):
            fail("screen binary path/FD identity differs")
        digest_before = file_sha256_fd(
            fd, fd_before.st_size, child_deadline
        )
        result.binary["sha256_before"] = digest_before
        if digest_before != expected_sha256:
            fail("screen binary SHA-256 differs from the presealed value")
        if health_args is not None:
            if health_modules is None:
                fail("health execution lacks sealed helper modules")
            if prepared_health_state is None:
                health_state = prepare_health(
                    health_args, cpu, child_deadline, health_modules
                )
            result.health = health_state["receipt"]
            if not health_state.get("ready"):
                fail("health admission did not license worker launch")
        if signal_guard is not None and signal_guard.first_signal is not None:
            fail(f"runner received {signal_guard.first_signal} before child start")
        fd_path = f"/proc/self/fd/{fd}"
        if time.monotonic() >= child_deadline:
            fail("controller reached the worker deadline before launch")
        if health_args is not None:
            health_state["sibling_start"] = read_cpu_tick_receipt(56)
            prelaunch_sampler_event = poll_sampler_supervision(
                health_state["sampler_handles"]
            )
            if prelaunch_sampler_event != "none":
                health_state["sampler_supervision_event"] = (
                    prelaunch_sampler_event
                )
                fail(
                    "sampler terminal state preceded worker launch: "
                    + prelaunch_sampler_event
                )
        child_start_ns = time.monotonic_ns()
        process = subprocess.Popen(
            [fd_path, "--run", "--cpu", str(cpu)], executable=fd_path,
            pass_fds=(fd,), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env={"LANG": "C", "LC_ALL": "C"},
            start_new_session=True,
            preexec_fn=pin_worker_before_exec,
        )
        result.binary["process_started"] = True
        result.binary["child_pid"] = process.pid
        result.binary["child_start_monotonic_ns"] = child_start_ns
        result.binary["child_process_start_ticks"] = process_start_ticks(
            process.pid
        )
        result.binary["child_affinity_after_spawn"] = sorted(
            os.sched_getaffinity(process.pid)
        )
        if result.binary["child_affinity_after_spawn"] != [EXPECTED_TARGET_CPU]:
            fail("worker did not inherit exact singleton CPU120 affinity")
        capture = bounded_capture(
            process, child_deadline, signal_guard,
            sampler_monitor_handles=health_state.get("sampler_handles"),
        )
        capture_started = True
        result.raw = capture.stdout
        result.stderr = capture.stderr
        result.binary["timed_out"] = capture.timed_out
        result.binary["output_overflow"] = capture.output_overflow
        result.binary["capture_error"] = capture.error
        health_state["sampler_supervision_event"] = capture.sampler_event
        result.returncode = process.returncode
        child_reap_ns: Optional[int] = (
            time.monotonic_ns() if process.returncode is not None else None
        )
        result.binary["child_reap_monotonic_ns"] = child_reap_ns
        if reap_reached_controller_deadline(child_reap_ns, controller_started_ns):
            result.binary["timed_out"] = True
            result.errors.append(
                "execution: child reap reached the inclusive 110-second deadline"
            )
        if capture.error != "none":
            result.errors.append(capture.error)
        if capture.sampler_event != "none":
            result.errors.append(
                "sampler supervision: " + capture.sampler_event
            )
    except BaseException as exc:
        result.errors.append(f"execution: {exception_text(exc)}")
    finally:
        if process is not None and not capture_started:
            error = terminate_process_group(process)
            if error is not None:
                result.errors.append(f"execution cleanup kill: {error}")
            try:
                process.wait(timeout=5.0)
            except BaseException as exc:
                result.errors.append(f"execution cleanup reap: {exception_text(exc)}")
            result.returncode = process.returncode
            if (
                process.returncode is not None
                and result.binary["process_started"]
            ):
                result.binary["child_reap_monotonic_ns"] = time.monotonic_ns()
                if reap_reached_controller_deadline(
                    result.binary["child_reap_monotonic_ns"],
                    controller_started_ns,
                ):
                    result.binary["timed_out"] = True
                    result.errors.append(
                        "execution cleanup: child reap reached the inclusive "
                        "110-second deadline"
                    )
            if signal_guard is not None:
                signal_guard.detach(process)
        if health_args is not None:
            if health_modules is None:
                health_finish_error = "sealed helper modules absent"
                result.errors.append("health finish: " + health_finish_error)
                partial_health = dict(result.health)
                partial_health["child_start_monotonic_ns"] = (
                    result.binary.get("child_start_monotonic_ns")
                )
                partial_health["child_reap_monotonic_ns"] = (
                    result.binary.get("child_reap_monotonic_ns")
                )
                collection = list(
                    partial_health.get("collection_failures", [])
                )
                append_unique_failure(
                    collection,
                    "health terminal collection: " + health_finish_error,
                )
                result.health = finalize_health_receipt(
                    partial_health, collection,
                    partial_health.get("violations", []),
                )
            else:
                try:
                    finished_health, health_errors = finish_health(
                        health_args, cpu,
                        result.binary.get("child_start_monotonic_ns"),
                        result.binary.get("child_reap_monotonic_ns"),
                        health_state, health_deadline, health_modules,
                    )
                    result.health = finished_health
                    result.errors.extend(health_errors)
                except BaseException as exc:
                    health_finish_error = exception_text(exc)
                    result.errors.append("health finish: " + health_finish_error)
                    partial_health = dict(result.health)
                    partial_health["child_start_monotonic_ns"] = (
                        result.binary.get("child_start_monotonic_ns")
                    )
                    partial_health["child_reap_monotonic_ns"] = (
                        result.binary.get("child_reap_monotonic_ns")
                    )
                    collection = list(
                        partial_health.get("collection_failures", [])
                    )
                    append_unique_failure(
                        collection,
                        "health terminal collection: " + health_finish_error,
                    )
                    result.health = finalize_health_receipt(
                        partial_health, collection,
                        partial_health.get("violations", []),
                    )
            try:
                close_sampler_admission_handles(
                    health_state.get("sampler_handles")
                )
            except BaseException as exc:
                close_error = exception_text(exc)
                result.errors.append("health descriptor close: " + close_error)
                partial_health = dict(result.health)
                collection = list(
                    partial_health.get("collection_failures", [])
                )
                append_unique_failure(
                    collection,
                    "health terminal collection: " + close_error,
                )
                result.health = finalize_health_receipt(
                    partial_health, collection,
                    partial_health.get("violations", []),
                )
        if fd >= 0:
            try:
                fd_after = os.fstat(fd)
                result.binary["stat_after"] = stat_receipt(fd_after)
                if fd_after.st_size < 1 or fd_after.st_size > MAX_BINARY_BYTES:
                    fail("binary postcheck image size differs")
                digest_after = file_sha256_fd(
                    fd, fd_after.st_size, finalization_deadline
                )
                result.binary["sha256_after"] = digest_after
                path_after = os.stat(str(binary), follow_symlinks=False)
                result.binary["path_stat_after"] = stat_receipt(path_after)
                if (
                    fd_before is None
                    or digest_before is None
                    or not same_file_receipt(fd_before, fd_after)
                    or not same_file_receipt(fd_after, path_after)
                    or digest_before != digest_after
                ):
                    result.errors.append("binary postcheck: image changed")
            except BaseException as exc:
                result.errors.append(f"binary postcheck: {exception_text(exc)}")
            try:
                os.close(fd)
            except BaseException as exc:
                result.errors.append(f"binary close: {exception_text(exc)}")
        execution_finished_ns = time.monotonic_ns()
        result.binary["execution_finished_monotonic_ns"] = execution_finished_ns
        result.elapsed_seconds = (
            execution_finished_ns - execution_started_ns
        ) / 1_000_000_000.0
    return result


def validate_git_index_flag_roster(data: bytes) -> None:
    expected = b"".join(
        ("H " + relative + "\n").encode("ascii")
        for relative in SOURCE_PATHS
    )
    if data != expected:
        fail(
            "screen source index flags differ: every path must have exact H status"
        )


def git_receipt(
    root: Path, expected_commit: str, expected_git_sha256: str,
    deadline: Optional[float] = None,
) -> Dict[str, Any]:
    if (
        not root.is_absolute() or os.path.normpath(str(root)) != str(root)
        or root.resolve(strict=True) != root or not root.is_dir()
    ):
        fail("Git receipt root is not one canonical absolute directory")
    environment = {
        "GIT_ASKPASS": "/bin/false", "GIT_NO_LAZY_FETCH": "1",
        "GIT_NO_REPLACE_OBJECTS": "1", "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0", "LANG": "C", "LC_ALL": "C",
        "PATH": "/usr/bin:/bin", "SSH_ASKPASS": "/bin/false",
    }
    lower_hash(expected_git_sha256, "expected Git executable SHA-256")
    receipt_deadline = (
        time.monotonic() + 10.0 if deadline is None else deadline
    )
    if (
        type(receipt_deadline) is not float
        or not math.isfinite(receipt_deadline)
    ):
        fail("Git receipt deadline differs")
    flags = nonblocking_read_flags("sealed Git executable")
    git_fd = os.open(str(GIT_EXECUTABLE), flags)
    git_before = os.fstat(git_fd)
    git_named_before = os.stat(str(GIT_EXECUTABLE), follow_symlinks=False)
    if (
        not stat.S_ISREG(git_before.st_mode)
        or not same_file_receipt(git_before, git_named_before)
        or not 1 <= git_before.st_size <= MAX_BINARY_BYTES
        or git_before.st_uid != 0 or git_before.st_gid != 0
        or git_before.st_nlink != 1
        or stat.S_IMODE(git_before.st_mode) != 0o755
    ):
        os.close(git_fd)
        fail("sealed /usr/bin/git executable policy differs")
    git_hash_before = file_sha256_fd(
        git_fd, git_before.st_size, receipt_deadline
    )
    if git_hash_before != expected_git_sha256:
        os.close(git_fd)
        fail("sealed /usr/bin/git executable hash differs")
    executable_receipt = {
        "path": str(GIT_EXECUTABLE),
        "sha256": git_hash_before,
        "stat": stat_receipt(git_before),
    }
    git_fd_path = f"/proc/self/fd/{git_fd}"
    def git(*arguments: str) -> bytes:
        if time.monotonic() >= receipt_deadline:
            fail("static Git receipt reached its global deadline")
        process = subprocess.Popen(
            [
                git_fd_path, "-c", "core.fsmonitor=false",
                "-c", "core.filemode=false",
                "-c", "safe.directory=" + str(root), *arguments,
            ],
            executable=git_fd_path, pass_fds=(git_fd,), cwd=root,
            env=environment, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True,
        )
        captured = bounded_capture(
            process, receipt_deadline,
            stdout_limit=1024 * 1024, stderr_limit=64 * 1024,
            kill_grace_seconds=0.5,
        )
        if (
            captured.timed_out or captured.output_overflow
            or captured.error != "none" or process.returncode != 0
            or captured.stderr
        ):
            fail(f"static git command failed: {' '.join(arguments)}")
        return captured.stdout

    try:
        worktree = git("rev-parse", "--show-toplevel")
        if worktree != (str(root) + "\n").encode("ascii"):
            fail("Git worktree root differs from the sealed source root")
        commit = git("rev-parse", "--verify", "HEAD^{commit}")
        if commit != (expected_commit + "\n").encode("ascii"):
            fail("Git HEAD differs from the expected source commit")
        tracked_status = git("status", "--porcelain=v1", "--untracked-files=no")
        if tracked_status:
            fail("tracked worktree is not clean")
        tracked_index_flags = git("ls-files", "-v")
        if (
            not tracked_index_flags
            or not tracked_index_flags.endswith(b"\n")
            or any(
                not line.startswith(b"H ")
                for line in tracked_index_flags.splitlines()
            )
        ):
            fail("tracked worktree contains non-normal index flags")
        tracked_paths = git("ls-files", "--error-unmatch", "--", *SOURCE_PATHS)
        expected_paths = b"".join(
            (relative + "\n").encode("ascii") for relative in SOURCE_PATHS
        )
        if tracked_paths != expected_paths:
            fail("screen source roster is not tracked in the expected order")
        index_flags = git(
            "ls-files", "-v", "--error-unmatch", "--", *SOURCE_PATHS
        )
        validate_git_index_flag_roster(index_flags)
        tree_rows = git(
            "ls-tree", "--full-tree", "-r", "HEAD", "--", *SOURCE_PATHS
        )
        tree_oids: Dict[str, str] = {}
        try:
            for raw_line in tree_rows.decode("ascii").splitlines():
                metadata, path = raw_line.split("\t", 1)
                mode, object_type, oid = metadata.split(" ")
                if (
                    mode != "100644" or object_type != "blob"
                    or LOWER40.fullmatch(oid) is None or path in tree_oids
                ):
                    fail("Git HEAD source blob row differs")
                tree_oids[path] = oid
        except (UnicodeDecodeError, ValueError) as exc:
            fail("Git HEAD source blob roster is malformed: " + str(exc))
        if set(tree_oids) != set(SOURCE_PATHS):
            fail("Git HEAD source blob roster differs")
        worktree_oids_raw = git(
            "hash-object", "--no-filters", "--", *SOURCE_PATHS
        )
        try:
            worktree_oids = worktree_oids_raw.decode("ascii").splitlines()
        except UnicodeDecodeError as exc:
            fail("Git raw worktree blob roster is malformed: " + str(exc))
        if (
            len(worktree_oids) != len(SOURCE_PATHS)
            or any(LOWER40.fullmatch(oid) is None for oid in worktree_oids)
        ):
            fail("Git raw worktree blob roster differs")
        source_blob_oids = [
            {
                "head_oid": tree_oids[path], "path": path,
                "worktree_oid": worktree_oid,
            }
            for path, worktree_oid in zip(SOURCE_PATHS, worktree_oids)
        ]
        if any(
            entry["head_oid"] != entry["worktree_oid"]
            for entry in source_blob_oids
        ):
            fail("raw source bytes differ from their exact HEAD blobs")
        git(
            "diff", "--quiet", "--no-ext-diff", "--no-textconv", "HEAD",
            "--", *SOURCE_PATHS,
        )
        git_after = os.fstat(git_fd)
        git_named_after = os.stat(str(GIT_EXECUTABLE), follow_symlinks=False)
        if (
            not same_file_receipt(git_before, git_after)
            or not same_file_receipt(git_after, git_named_after)
            or file_sha256_fd(
                git_fd, git_after.st_size, receipt_deadline
            ) != git_hash_before
        ):
            fail("Git executable changed during receipt collection")
        return {
            "executable": executable_receipt,
            "head": expected_commit,
            "source_blob_oids": source_blob_oids,
            "source_blob_roster_sha256": sha256_bytes(
                canonical_bytes(source_blob_oids)
            ),
            "source_roster_sha256": sha256_bytes(tracked_paths),
            "tracked_index_flags_sha256": sha256_bytes(tracked_index_flags),
            "tracked_status_sha256": sha256_bytes(tracked_status),
            "worktree_root": str(root),
        }
    finally:
        os.close(git_fd)


def build_manifest(
    root: Path, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    if (
        not root.is_absolute() or os.path.normpath(str(root)) != str(root)
        or root.resolve(strict=True) != root or not root.is_dir()
    ):
        fail("build directory is not one canonical absolute directory")
    entries: List[Dict[str, Any]] = []
    for relative in BUILD_PATHS:
        if deadline is not None and time.monotonic() >= deadline:
            fail("build manifest reached its global deadline")
        path = root / relative
        digest, info = hash_regular_path(
            path, "build provenance " + relative,
            max_size=64 * 1024 * 1024, deadline=deadline,
        )
        if (
            stat.S_IMODE(info.st_mode) != 0o444
            or info.st_uid != EXPECTED_CAMPAIGN_UID
            or info.st_gid != EXPECTED_CAMPAIGN_GID
        ):
            fail(
                f"sealed build provenance {relative} must be "
                "UID/GID 1000 mode 0444"
            )
        entries.append({
            "bytes": info.st_size,
            "device": info.st_dev,
            "gid": info.st_gid,
            "inode": info.st_ino,
            "mode": stat.S_IMODE(info.st_mode),
            "nlink": info.st_nlink,
            "path": relative,
            "sha256": digest,
            "uid": info.st_uid,
        })
    value: Dict[str, Any] = {
        "entries": entries,
        "root": str(root),
        "schema": "wirehair.wh2.direct-systematic-complement-build.v2",
        "sha256": None,
    }
    value["sha256"] = sha256_bytes(canonical_bytes(value))
    return value


def controller_provenance(
    args: argparse.Namespace, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    process_start_before = process_start_ticks(os.getpid())
    process_security_before = read_process_security(
        os.getpid(), "controller provenance before"
    )
    validate_process_security(
        process_security_before, args.expected_controller_uid,
        args.expected_controller_gid, [], "controller provenance before",
    )
    script_path = Path(__file__).resolve(strict=True)
    source_root = script_path.parents[1]
    expected_script = source_root / "bench/Wh2DirectSystematicComplementScreen.py"
    if script_path != expected_script:
        fail("controller script path differs from its source-manifest path")
    script_hash, script_info = hash_regular_path(
        script_path, "controller script", max_size=MAX_SOURCE_FILE_BYTES,
        deadline=deadline,
    )
    exact_string(
        script_hash, args.expected_controller_sha256,
        "controller script SHA-256",
    )
    if (
        stat.S_IMODE(script_info.st_mode) != 0o444
        or script_info.st_uid != args.expected_controller_uid
        or script_info.st_gid != args.expected_controller_gid
    ):
        fail("controller script mode/owner seal differs")
    declared_python = Path(sys.executable)
    if not declared_python.is_absolute():
        fail("Python executable path is not absolute")
    python_path = declared_python.resolve(strict=True)
    python_hash, python_info = hash_regular_path(
        python_path, "Python executable", single_link=False,
        max_size=MAX_BINARY_BYTES, deadline=deadline,
    )
    live_python_path = Path("/proc/self/exe").resolve(strict=True)
    live_fd = os.open(
        "/proc/self/exe",
        nonblocking_read_flags("live Python executable", nofollow=False),
    )
    try:
        live_python_info = os.fstat(live_fd)
        if (
            live_python_info.st_size < 1
            or live_python_info.st_size > MAX_BINARY_BYTES
        ):
            fail("live Python executable size differs")
        live_python_hash = file_sha256_fd(
            live_fd, live_python_info.st_size, deadline
        )
        live_python_after = os.fstat(live_fd)
        if (
            live_python_path != python_path
            or live_python_hash != python_hash
            or not same_file_receipt(live_python_info, python_info)
            or not same_file_receipt(live_python_info, live_python_after)
        ):
            fail("declared Python executable differs from /proc/self/exe")
    finally:
        os.close(live_fd)
    exact_string(
        python_hash, args.expected_python_sha256, "Python executable SHA-256"
    )
    if (
        stat.S_IMODE(python_info.st_mode) != 0o755
        or python_info.st_uid != 0 or python_info.st_gid != 0
    ):
        fail("Python executable mode/owner differs")
    git_hash, git_info = hash_regular_path(
        GIT_EXECUTABLE, "Git executable", max_size=MAX_BINARY_BYTES,
        deadline=deadline,
    )
    if (
        git_hash != args.expected_git_sha256
        or stat.S_IMODE(git_info.st_mode) != 0o755
        or git_info.st_uid != 0 or git_info.st_gid != 0
        or git_info.st_nlink != 1
    ):
        fail("Git executable provenance differs")
    process_security_after = read_process_security(
        os.getpid(), "controller provenance after"
    )
    process_info = os.stat(f"/proc/{os.getpid()}", follow_symlinks=False)
    if (
        process_start_ticks(os.getpid()) != process_start_before
        or not canonical_equal(process_security_before, process_security_after)
        or process_info.st_uid != args.expected_controller_uid
        or process_info.st_gid != args.expected_controller_gid
    ):
        fail("controller process identity changed during provenance collection")
    argv_bytes = b"".join(os.fsencode(value) + b"\0" for value in sys.argv)
    value: Dict[str, Any] = {
        "argv": list(sys.argv),
        "argv_sha256": sha256_bytes(argv_bytes),
        "controller_cpu": EXPECTED_CONTROLLER_CPU,
        "dont_write_bytecode": bool(sys.dont_write_bytecode),
        "environment": dict(os.environ),
        "git_path": str(GIT_EXECUTABLE),
        "git_sha256": git_hash,
        "git_stat": stat_receipt(git_info),
        "isolated": bool(sys.flags.isolated),
        "optimize": sys.flags.optimize,
        "pid": os.getpid(),
        "process_gid": process_info.st_gid,
        "process_security": process_security_after,
        "process_start_ticks": process_start_before,
        "process_uid": process_info.st_uid,
        "python_declared_path": str(declared_python),
        "python_path": str(python_path),
        "python_sha256": python_hash,
        "python_stat": stat_receipt(python_info),
        "receipt_sha256": None,
        "schema": CONTROLLER_PROVENANCE_SCHEMA,
        "script_path": str(script_path),
        "script_sha256": script_hash,
        "script_stat": stat_receipt(script_info),
        "singleton_affinity": sorted(os.sched_getaffinity(0)),
    }
    if (
        not value["isolated"] or not value["dont_write_bytecode"]
        or value["optimize"] != 0
        or not canonical_equal(value["environment"], SEALED_LAUNCH_ENVIRONMENT)
        or value["singleton_affinity"] != [EXPECTED_CONTROLLER_CPU]
    ):
        fail("controller interpreter/affinity seal differs")
    value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
    return value


PREFLIGHT_CONTROLLER_IMAGE_KEYS = {
    "git_path", "git_sha256", "git_stat", "python_declared_path",
    "python_path", "python_sha256", "python_stat", "script_path",
    "script_sha256", "script_stat",
}
PREFLIGHT_CONTROLLER_KEYS = {
    "argv", "argv_sha256", "controller_cpu", "dont_write_bytecode",
    "environment", "image", "isolated", "optimize", "pid", "process_gid",
    "process_security", "process_start_ticks", "process_uid",
    "receipt_sha256", "schema", "singleton_affinity",
}


def preflight_controller_receipt(
    run_args: argparse.Namespace, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    observed = controller_provenance(run_args, deadline)
    image = {
        name: observed[name] for name in PREFLIGHT_CONTROLLER_IMAGE_KEYS
    }
    value: Dict[str, Any] = {
        "argv": observed["argv"],
        "argv_sha256": observed["argv_sha256"],
        "controller_cpu": observed["controller_cpu"],
        "dont_write_bytecode": observed["dont_write_bytecode"],
        "environment": observed["environment"],
        "image": image,
        "isolated": observed["isolated"],
        "optimize": observed["optimize"],
        "pid": observed["pid"],
        "process_gid": observed["process_gid"],
        "process_security": observed["process_security"],
        "process_start_ticks": observed["process_start_ticks"],
        "process_uid": observed["process_uid"],
        "receipt_sha256": None,
        "schema": PREFLIGHT_CONTROLLER_SCHEMA,
        "singleton_affinity": observed["singleton_affinity"],
    }
    value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
    validate_preflight_controller_receipt(value, "preflight controller")
    return value


def validate_preflight_controller_receipt(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, PREFLIGHT_CONTROLLER_KEYS, where)
    exact_string(value["schema"], PREFLIGHT_CONTROLLER_SCHEMA, where + ".schema")
    exact_int(
        value["controller_cpu"], EXPECTED_CONTROLLER_CPU,
        EXPECTED_CONTROLLER_CPU, where + ".controller_cpu",
    )
    exact_int(value["pid"], 1, MAX_INT63, where + ".pid")
    exact_int(
        value["process_uid"], EXPECTED_CAMPAIGN_UID,
        EXPECTED_CAMPAIGN_UID, where + ".process_uid",
    )
    exact_int(
        value["process_gid"], EXPECTED_CAMPAIGN_GID,
        EXPECTED_CAMPAIGN_GID, where + ".process_gid",
    )
    exact_int(
        value["process_start_ticks"], 1, MAX_UINT64,
        where + ".process_start_ticks",
    )
    validate_process_security(
        value["process_security"], EXPECTED_CAMPAIGN_UID,
        EXPECTED_CAMPAIGN_GID, [], where + ".process_security",
    )
    exact_int(value["optimize"], 0, 0, where + ".optimize")
    affinity = exact_int_list(
        value["singleton_affinity"], where + ".singleton_affinity",
        length=1, sorted_unique=True,
    )
    if (
        value["isolated"] is not True
        or value["dont_write_bytecode"] is not True
        or affinity != [EXPECTED_CONTROLLER_CPU]
        or not canonical_equal(value["environment"], SEALED_LAUNCH_ENVIRONMENT)
    ):
        fail(f"{where} interpreter/environment/affinity policy differs")
    argv = value["argv"]
    if (
        type(argv) is not list
        or len(argv) != 2 + 2 * len(PREFLIGHT_SEAL_OPTION_ORDER)
        or argv[0] != str(Path(__file__).resolve(strict=True))
        or argv[1] != "--preflight-seal"
        or any(
            type(item) is not str or not item or len(item) > 4096
            or "\0" in item or not item.isascii()
            for item in argv
        )
        or any(
            argv[2 + 2 * index] != option
            for index, option in enumerate(PREFLIGHT_SEAL_OPTION_ORDER)
        )
    ):
        fail(f"{where} argv differs")
    argv_digest = sha256_bytes(
        b"".join(os.fsencode(item) + b"\0" for item in argv)
    )
    if argv_digest != lower_hash(value["argv_sha256"], where + ".argv_sha256"):
        fail(f"{where} argv hash differs")
    image = value["image"]
    if type(image) is not dict:
        fail(f"{where}.image is not an object")
    exact_keys(image, PREFLIGHT_CONTROLLER_IMAGE_KEYS, where + ".image")
    for name in ("git_path", "python_declared_path", "python_path", "script_path"):
        exact_absolute_path(image[name], f"{where}.image.{name}")
    exact_string(image["git_path"], str(GIT_EXECUTABLE), where + ".image.git_path")
    exact_string(
        image["python_declared_path"], EXPECTED_SAMPLER_PYTHON,
        where + ".image.python_declared_path",
    )
    exact_string(
        image["python_path"],
        str(Path(EXPECTED_SAMPLER_PYTHON).resolve(strict=True)),
        where + ".image.python_path",
    )
    exact_string(
        image["script_path"], str(Path(__file__).resolve(strict=True)),
        where + ".image.script_path",
    )
    for name in ("git_sha256", "python_sha256", "script_sha256"):
        lower_hash(image[name], f"{where}.image.{name}")
    for name in ("git_stat", "python_stat", "script_stat"):
        validate_stat_receipt(image[name], f"{where}.image.{name}")
    if (
        image["git_stat"]["mode"] != 0o755
        or image["git_stat"]["uid"] != 0
        or image["git_stat"]["gid"] != 0
        or image["git_stat"]["nlink"] != 1
        or not 1 <= image["git_stat"]["size"] <= MAX_BINARY_BYTES
        or image["python_stat"]["mode"] != 0o755
        or image["python_stat"]["uid"] != 0
        or image["python_stat"]["gid"] != 0
        or not 1 <= image["python_stat"]["size"] <= MAX_BINARY_BYTES
        or image["script_stat"]["mode"] != 0o444
        or image["script_stat"]["uid"] != EXPECTED_CAMPAIGN_UID
        or image["script_stat"]["gid"] != EXPECTED_CAMPAIGN_GID
        or image["script_stat"]["nlink"] != 1
        or not 1 <= image["script_stat"]["size"] <= MAX_SOURCE_FILE_BYTES
    ):
        fail(f"{where} image policy differs")
    digest = lower_hash(value["receipt_sha256"], where + ".receipt_sha256")
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail(f"{where} self-hash differs")


@dataclass
class ClaimReservation:
    fd: int
    parent_fd: int
    identity: Tuple[int, int, int]
    parent_identity: Tuple[int, int, int, int, int]
    document: Dict[str, Any]
    raw: bytes

    @classmethod
    def reserve(
        cls, args: argparse.Namespace, controller: Mapping[str, Any],
        controller_started_ns: int,
    ) -> "ClaimReservation":
        require_preflight_side_effects_allowed("claim reservation")
        require_authorized_process_identity(
            args.expected_controller_uid, args.expected_controller_gid
        )
        parent = FIXED_CLAIM_PATH.parent
        if parent != Path("/var/tmp") or not parent.is_dir():
            fail("fixed claim parent differs")
        flags = (
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        )
        parent_fd = os.open(str(parent), flags)
        fd = -1
        try:
            parent_info = os.fstat(parent_fd)
            parent_named = os.stat(str(parent), follow_symlinks=False)
            parent_identity = (
                parent_info.st_dev, parent_info.st_ino,
                stat.S_IMODE(parent_info.st_mode), parent_info.st_uid,
                parent_info.st_gid,
            )
            if (
                not stat.S_ISDIR(parent_info.st_mode)
                or not stat.S_ISDIR(parent_named.st_mode)
                or (parent_info.st_dev, parent_info.st_ino)
                != (parent_named.st_dev, parent_named.st_ino)
                or parent_identity[2:] != (0o1777, 0, 0)
            ):
                fail("fixed claim parent must be held root:root mode 01777 /var/tmp")
            fd = os.open(
                FIXED_CLAIM_PATH.name,
                os.O_RDWR | os.O_CREAT | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
                0o600, dir_fd=parent_fd,
            )
            document: Dict[str, Any] = {
                "binary_sha256": args.expected_binary_sha256,
                "build_manifest_sha256": args.expected_build_manifest_sha256,
                "campaign": CAMPAIGN,
                "controller_receipt_sha256": controller["receipt_sha256"],
                "controller_started_monotonic_ns": controller_started_ns,
                "git_sha256": args.expected_git_sha256,
                "gid": os.getgid(),
                "output_path": str(FIXED_OUTPUT_DIR),
                "parent_device": parent_identity[0],
                "parent_inode": parent_identity[1],
                "pid": os.getpid(),
                "process_start_ticks": controller["process_start_ticks"],
                "schema": CLAIM_SCHEMA,
                "source_commit": args.expected_source_commit,
                "source_manifest_sha256": args.expected_source_manifest_sha256,
                "uid": os.getuid(),
            }
            raw = canonical_bytes(document) + b"\n"
            OutputBundle._write_exact(fd, raw)
            os.fchmod(fd, 0o400)
            os.fsync(fd)
            os.fsync(parent_fd)
            info = os.fstat(fd)
            named = os.stat(
                FIXED_CLAIM_PATH.name, dir_fd=parent_fd,
                follow_symlinks=False,
            )
            if (
                not same_file_receipt(info, named)
                or stat.S_IMODE(info.st_mode) != 0o400
                or info.st_nlink != 1
                or info.st_uid != document["uid"]
                or info.st_gid != document["gid"]
            ):
                fail("fixed claim identity/mode differs")
            return cls(
                fd=fd, parent_fd=parent_fd,
                identity=(info.st_dev, info.st_ino, info.st_nlink),
                parent_identity=parent_identity,
                document=document, raw=raw,
            )
        except BaseException:
            if fd >= 0:
                os.close(fd)
            os.close(parent_fd)
            raise

    def receipt(self) -> Dict[str, Any]:
        info = os.fstat(self.fd)
        return {
            "bytes": len(self.raw),
            "device": self.identity[0],
            "inode": self.identity[1],
            "nlink": self.identity[2],
            "parent": {
                "device": self.parent_identity[0],
                "gid": self.parent_identity[4],
                "inode": self.parent_identity[1],
                "mode": self.parent_identity[2],
                "path": str(FIXED_CLAIM_PATH.parent),
                "uid": self.parent_identity[3],
            },
            "path": str(FIXED_CLAIM_PATH),
            "sha256": sha256_bytes(self.raw),
            "stat": stat_receipt(info),
        }

    def verify(self) -> None:
        os.fsync(self.fd)
        os.fsync(self.parent_fd)
        parent = os.fstat(self.parent_fd)
        parent_named = os.stat(str(FIXED_CLAIM_PATH.parent), follow_symlinks=False)
        info = os.fstat(self.fd)
        named = os.stat(
            FIXED_CLAIM_PATH.name, dir_fd=self.parent_fd,
            follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(parent.st_mode)
            or not stat.S_ISDIR(parent_named.st_mode)
            or (
                parent.st_dev, parent.st_ino, stat.S_IMODE(parent.st_mode),
                parent.st_uid, parent.st_gid,
            ) != self.parent_identity
            or (parent.st_dev, parent.st_ino)
            != (parent_named.st_dev, parent_named.st_ino)
            or not same_file_receipt(info, named)
            or (info.st_dev, info.st_ino, info.st_nlink) != self.identity
            or stat.S_IMODE(info.st_mode) != 0o400
            or info.st_uid != self.document["uid"]
            or info.st_gid != self.document["gid"]
            or info.st_size != len(self.raw)
            or file_sha256_fd(self.fd, info.st_size) != sha256_bytes(self.raw)
        ):
            fail("fixed claim changed")


def reserve_output_file(directory_fd: int, name: str, mode: int) -> int:
    return os.open(
        name,
        os.O_RDWR | os.O_CREAT | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        mode, dir_fd=directory_fd,
    )


def held_artifact(fd: int, expected_mode: int) -> Dict[str, Any]:
    info = os.fstat(fd)
    if (
        not stat.S_ISREG(info.st_mode) or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != expected_mode
        or not 0 <= info.st_size <= MAX_FINAL_ARTIFACT_BYTES
    ):
        fail("held artifact policy differs")
    digest = file_sha256_fd(fd, info.st_size)
    after = os.fstat(fd)
    if not same_file_receipt(info, after):
        fail("held artifact changed during hashing")
    return {
        "bytes": info.st_size,
        "device": info.st_dev,
        "gid": info.st_gid,
        "inode": info.st_ino,
        "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink,
        "sha256": digest,
        "uid": info.st_uid,
    }


def controller_document_bytes(value: Dict[str, Any]) -> bytes:
    value["receipt_sha256"] = None
    value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
    return canonical_bytes(value) + b"\n"


ARTIFACT_RECEIPT_KEYS = {
    "bytes", "device", "gid", "inode", "mode", "nlink", "sha256", "uid",
}
BUILD_ENTRY_KEYS = ARTIFACT_RECEIPT_KEYS | {"path"}
CONTROLLER_PROVENANCE_KEYS = {
    "argv", "argv_sha256", "controller_cpu", "dont_write_bytecode",
    "environment", "git_path", "git_sha256", "git_stat", "isolated",
    "optimize", "pid", "process_gid", "process_security",
    "process_start_ticks", "process_uid",
    "python_declared_path", "python_path", "python_sha256", "python_stat",
    "receipt_sha256", "schema", "script_path", "script_sha256",
    "script_stat", "singleton_affinity",
}


def validate_artifact_receipt(value: Any, where: str, mode: int = 0o400) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, ARTIFACT_RECEIPT_KEYS, where)
    for name, lower in (
        ("bytes", 0), ("device", 0), ("gid", 0), ("inode", 1),
        ("nlink", 1), ("uid", 0),
    ):
        exact_int(value[name], lower, MAX_UINT64, f"{where}.{name}")
    exact_int(value["mode"], mode, mode, f"{where}.mode")
    exact_int(value["nlink"], 1, 1, f"{where}.nlink")
    lower_hash(value["sha256"], f"{where}.sha256")


def validate_build_manifest(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, {"entries", "root", "schema", "sha256"}, where)
    exact_string(
        value["schema"],
        "wirehair.wh2.direct-systematic-complement-build.v2",
        where + ".schema",
    )
    exact_absolute_path(value["root"], where + ".root")
    entries = value["entries"]
    if type(entries) is not list or len(entries) != len(BUILD_PATHS):
        fail(f"{where} entry roster differs")
    for index, (entry, expected_path) in enumerate(zip(entries, BUILD_PATHS)):
        entry_where = f"{where}.entries[{index}]"
        if type(entry) is not dict:
            fail(f"{entry_where} is not an object")
        exact_keys(entry, BUILD_ENTRY_KEYS, entry_where)
        exact_string(entry["path"], expected_path, entry_where + ".path")
        validate_artifact_receipt(
            {name: entry[name] for name in ARTIFACT_RECEIPT_KEYS},
            entry_where + ".artifact", mode=0o444,
        )
        if (
            entry["bytes"] > 64 * 1024 * 1024
            or entry["uid"] != EXPECTED_CAMPAIGN_UID
            or entry["gid"] != EXPECTED_CAMPAIGN_GID
        ):
            fail(f"{entry_where} exceeds the build provenance size bound")
    digest = lower_hash(value["sha256"], where + ".sha256")
    preimage = dict(value)
    preimage["sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail(f"{where} self-hash differs")


def validate_controller_provenance(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, CONTROLLER_PROVENANCE_KEYS, where)
    exact_string(
        value["schema"], CONTROLLER_PROVENANCE_SCHEMA, where + ".schema"
    )
    exact_int(
        value["controller_cpu"], EXPECTED_CONTROLLER_CPU,
        EXPECTED_CONTROLLER_CPU, where + ".controller_cpu",
    )
    exact_int(value["pid"], 1, MAX_INT63, where + ".pid")
    exact_int(value["process_uid"], 0, MAX_UINT32, where + ".process_uid")
    exact_int(value["process_gid"], 0, MAX_UINT32, where + ".process_gid")
    validate_process_security(
        value["process_security"], value["process_uid"], value["process_gid"],
        [], where + ".process_security",
    )
    exact_int(
        value["process_start_ticks"], 1, MAX_UINT64,
        where + ".process_start_ticks",
    )
    exact_int(value["optimize"], 0, 0, where + ".optimize")
    if not canonical_equal(value["environment"], SEALED_LAUNCH_ENVIRONMENT):
        fail(f"{where} startup environment differs")
    affinity = exact_int_list(
        value["singleton_affinity"], where + ".singleton_affinity",
        length=1, sorted_unique=True,
    )
    if (
        value["isolated"] is not True
        or value["dont_write_bytecode"] is not True
        or affinity[0] != EXPECTED_CONTROLLER_CPU
    ):
        fail(f"{where} interpreter/affinity policy differs")
    argv = value["argv"]
    expected_run_argv_length = 2 + 2 * len(RUN_ONCE_OPTION_ORDER)
    if (
        type(argv) is not list or len(argv) != expected_run_argv_length
        or argv[1] != "--run-once"
        or any(type(item) is not str or len(item) > 4096 for item in argv)
    ):
        fail(f"{where} argv differs")
    argv_bytes = b"".join(os.fsencode(item) + b"\0" for item in argv)
    if sha256_bytes(argv_bytes) != lower_hash(
        value["argv_sha256"], where + ".argv_sha256"
    ):
        fail(f"{where} argv hash differs")
    for name in ("git_path", "python_declared_path", "python_path", "script_path"):
        exact_absolute_path(value[name], f"{where}.{name}")
    exact_string(value["git_path"], str(GIT_EXECUTABLE), where + ".git_path")
    lower_hash(value["git_sha256"], where + ".git_sha256")
    lower_hash(value["python_sha256"], where + ".python_sha256")
    lower_hash(value["script_sha256"], where + ".script_sha256")
    validate_stat_receipt(value["git_stat"], where + ".git_stat")
    validate_stat_receipt(value["python_stat"], where + ".python_stat")
    validate_stat_receipt(value["script_stat"], where + ".script_stat")
    if (
        value["python_stat"]["mode"] != 0o755
        or value["python_stat"]["uid"] != 0
        or value["python_stat"]["gid"] != 0
        or not 1 <= value["python_stat"]["size"] <= MAX_BINARY_BYTES
        or value["git_stat"]["mode"] != 0o755
        or value["git_stat"]["uid"] != 0
        or value["git_stat"]["gid"] != 0
        or value["git_stat"]["nlink"] != 1
        or not 1 <= value["git_stat"]["size"] <= MAX_BINARY_BYTES
        or value["script_stat"]["mode"] != 0o444
        or value["script_stat"]["size"] > MAX_SOURCE_FILE_BYTES
    ):
        fail(f"{where} executable/script policy differs")
    digest = lower_hash(value["receipt_sha256"], where + ".receipt_sha256")
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail(f"{where} self-hash differs")


def validate_claim_receipt(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, {
        "bytes", "device", "inode", "nlink", "parent", "path", "sha256",
        "stat",
    }, where)
    exact_absolute_path(value["path"], where + ".path")
    exact_string(value["path"], str(FIXED_CLAIM_PATH), where + ".path")
    exact_int(value["bytes"], 1, 1024 * 1024, where + ".bytes")
    exact_int(value["device"], 0, MAX_UINT64, where + ".device")
    exact_int(value["inode"], 1, MAX_UINT64, where + ".inode")
    exact_int(value["nlink"], 1, 1, where + ".nlink")
    lower_hash(value["sha256"], where + ".sha256")
    validate_stat_receipt(value["stat"], where + ".stat")
    if (
        value["stat"]["device"] != value["device"]
        or value["stat"]["inode"] != value["inode"]
        or value["stat"]["nlink"] != 1
        or value["stat"]["mode"] != 0o400
        or value["stat"]["size"] != value["bytes"]
    ):
        fail(f"{where} stat binding differs")
    parent = value["parent"]
    if type(parent) is not dict:
        fail(f"{where}.parent is not an object")
    exact_keys(parent, {"device", "gid", "inode", "mode", "path", "uid"}, where + ".parent")
    exact_string(parent["path"], "/var/tmp", where + ".parent.path")
    for name, expected in (("mode", 0o1777), ("uid", 0), ("gid", 0)):
        exact_int(parent[name], expected, expected, f"{where}.parent.{name}")
    exact_int(parent["device"], 0, MAX_UINT64, where + ".parent.device")
    exact_int(parent["inode"], 1, MAX_UINT64, where + ".parent.inode")


def parse_claim_document(data: bytes) -> Dict[str, Any]:
    if (
        not data or len(data) > 1024 * 1024 or not data.endswith(b"\n")
        or data.count(b"\n") != 1 or b"\r" in data
    ):
        fail("claim document framing differs")
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(f"claim document is malformed: {exc}")
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail("claim document is not canonical JSON")
    exact_keys(value, {
        "binary_sha256", "build_manifest_sha256", "campaign",
        "controller_receipt_sha256", "controller_started_monotonic_ns",
        "gid", "git_sha256", "output_path", "parent_device", "parent_inode",
        "pid", "process_start_ticks", "schema", "source_commit",
        "source_manifest_sha256", "uid",
    }, "claim document")
    exact_string(value["schema"], CLAIM_SCHEMA, "claim schema")
    exact_string(value["campaign"], CAMPAIGN, "claim campaign")
    exact_string(value["output_path"], str(FIXED_OUTPUT_DIR), "claim output")
    for name in (
        "binary_sha256", "build_manifest_sha256",
        "controller_receipt_sha256", "git_sha256",
        "source_manifest_sha256",
    ):
        lower_hash(value[name], "claim " + name)
    if (
        type(value["source_commit"]) is not str
        or LOWER40.fullmatch(value["source_commit"]) is None
    ):
        fail("claim source commit differs")
    for name, lower, upper in (
        ("controller_started_monotonic_ns", 0, MAX_INT63),
        ("gid", 0, MAX_UINT32),
        ("parent_device", 0, MAX_UINT64),
        ("parent_inode", 1, MAX_UINT64),
        ("pid", 1, MAX_INT63),
        ("process_start_ticks", 1, MAX_UINT64),
        ("uid", 0, MAX_UINT32),
    ):
        exact_int(value[name], lower, upper, "claim " + name)
    return value


def validate_claim_runtime_binding(
    document: Mapping[str, Any], receipt: Mapping[str, Any],
    expected_document: Mapping[str, Any], raw: bytes,
) -> None:
    validate_claim_receipt(receipt, "bound claim receipt")
    if (
        not canonical_equal(document, expected_document)
        or receipt["sha256"] != sha256_bytes(raw)
        or receipt["bytes"] != len(raw)
        or receipt["parent"]["device"] != document["parent_device"]
        or receipt["parent"]["inode"] != document["parent_inode"]
        or receipt["stat"]["uid"] != document["uid"]
        or receipt["stat"]["gid"] != document["gid"]
    ):
        fail("claim document/runtime provenance binding differs")


def validate_controller_source_git_binding(
    controller: Mapping[str, Any], git: Mapping[str, Any],
    source: Mapping[str, Any], where: str,
) -> None:
    script_relative = "bench/Wh2DirectSystematicComplementScreen.py"
    script_entry = source["entries"][SOURCE_PATHS.index(script_relative)]
    executable = git["executable"]
    if (
        controller["script_path"]
        != str(Path(__file__).resolve(strict=True))
        or script_entry["path"] != script_relative
        or script_entry["sha256"] != controller["script_sha256"]
        or script_entry["bytes"] != controller["script_stat"]["size"]
        or not canonical_equal(
            script_entry["stat"], controller["script_stat"]
        )
        or executable["path"] != controller["git_path"]
        or executable["sha256"] != controller["git_sha256"]
        or not canonical_equal(executable["stat"], controller["git_stat"])
    ):
        fail(f"{where} controller/source/Git provenance binding differs")
    for source_entry, blob_entry in zip(
        source["entries"], git["source_blob_oids"]
    ):
        if (
            source_entry["path"] != blob_entry["path"]
            or source_entry["git_blob_oid"] != blob_entry["head_oid"]
            or source_entry["git_blob_oid"] != blob_entry["worktree_oid"]
        ):
            fail(f"{where} source bytes/HEAD blob binding differs")


def parse_controller_document(data: bytes) -> Dict[str, Any]:
    if (
        not data or len(data) > 8 * 1024 * 1024 or not data.endswith(b"\n")
        or data.count(b"\n") != 1 or b"\r" in data
    ):
        fail("controller document framing differs")
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(f"controller document is malformed: {exc}")
    expected_keys = {
        "artifacts", "binary", "build_after", "build_before", "campaign",
        "claim", "controller_after", "controller_before",
        "controller_deadline_seconds", "controller_elapsed_seconds",
        "controller_observed_monotonic_ns",
        "controller_started_monotonic_ns", "failure", "git_after",
        "git_before", "health_receipt_sha256", "outcome",
        "outer_deadline_seconds", "output_path", "receipt_sha256", "schema",
        "sampler_after", "signals", "source_manifest_after", "source_manifest_before",
        "summary", "target_cpu",
    }
    if type(value) is not dict:
        fail("controller document is not an object")
    exact_keys(value, expected_keys, "controller document")
    if canonical_bytes(value) + b"\n" != data:
        fail("controller document is not canonical")
    exact_string(value["schema"], CONTROLLER_SCHEMA, "controller schema")
    exact_string(value["campaign"], CAMPAIGN, "controller campaign")
    exact_string(value["output_path"], str(FIXED_OUTPUT_DIR), "controller output")
    if value["outcome"] not in ("pass", "reject", "invalid"):
        fail("controller outcome differs")
    if (
        type(value["failure"]) is not str or not value["failure"]
        or len(value["failure"]) > MAX_FAILURE_TEXT_CHARS
    ):
        fail("controller failure text differs")
    exact_int(
        value["target_cpu"], EXPECTED_TARGET_CPU, EXPECTED_TARGET_CPU,
        "controller target CPU",
    )
    for field, expected in (
        ("controller_deadline_seconds", CONTROLLER_DEADLINE_SECONDS),
        ("outer_deadline_seconds", OUTER_DEADLINE_SECONDS),
    ):
        if type(value[field]) is not float or value[field] != expected:
            fail(f"controller {field} differs")
    controller_started_ns = exact_int(
        value["controller_started_monotonic_ns"], 0, MAX_INT63,
        "controller start timestamp",
    )
    controller_observed_ns = exact_int(
        value["controller_observed_monotonic_ns"], controller_started_ns,
        MAX_INT63, "controller observation timestamp",
    )
    expected_elapsed = (
        controller_observed_ns - controller_started_ns
    ) / 1_000_000_000.0
    if (
        type(value["controller_elapsed_seconds"]) is not float
        or not math.isfinite(value["controller_elapsed_seconds"])
        or value["controller_elapsed_seconds"] != expected_elapsed
        or controller_observed_ns >= (
            controller_started_ns
            + int(FINAL_COMMIT_START_DEADLINE_SECONDS * 1_000_000_000)
        )
    ):
        fail("controller elapsed/observation time differs")
    lower_hash(value["health_receipt_sha256"], "controller health receipt")
    signals = validate_bounded_messages(value["signals"], "controller signals")
    if any(name not in ("SIGINT", "SIGTERM", "SIGHUP") for name in signals):
        fail("controller signal roster differs")
    for name in (
        "artifacts", "binary", "build_after", "build_before", "claim",
        "controller_after", "controller_before", "git_after", "git_before",
        "sampler_after", "source_manifest_after", "source_manifest_before", "summary",
    ):
        if type(value[name]) is not dict:
            fail(f"controller {name} is not an object")
    artifacts = value["artifacts"]
    exact_keys(
        artifacts,
        {"raw.jsonl", "stderr.txt", "summary.json", "thermal.csv"},
        "controller artifacts",
    )
    for name, receipt in artifacts.items():
        validate_artifact_receipt(receipt, "controller artifact " + name)
    summary = value["summary"]
    exact_keys(
        summary, {"outcome", "preimage_sha256", "sha256"},
        "controller summary binding",
    )
    if summary["outcome"] not in ("pass", "reject", "invalid"):
        fail("controller summary outcome differs")
    lower_hash(summary["preimage_sha256"], "controller summary preimage")
    lower_hash(summary["sha256"], "controller summary SHA-256")
    if summary["sha256"] != artifacts["summary.json"]["sha256"]:
        fail("controller summary artifact hash differs")
    validate_claim_receipt(value["claim"], "controller claim")
    validate_binary_receipt(
        value["binary"], require_complete=value["outcome"] in ("pass", "reject")
    )
    execution_start = value["binary"]["execution_started_monotonic_ns"]
    execution_finish = value["binary"]["execution_finished_monotonic_ns"]
    if execution_start is not None and (
        execution_start < controller_started_ns
        or execution_finish > controller_observed_ns
    ):
        fail("binary execution interval escapes the controller interval")
    validate_build_manifest(value["build_before"], "controller build before")
    validate_controller_provenance(
        value["controller_before"], "controller provenance before"
    )
    validate_source_manifest_receipt(
        value["source_manifest_before"], "controller source before"
    )
    if (
        type(value["git_before"].get("head")) is not str
        or LOWER40.fullmatch(value["git_before"]["head"]) is None
    ):
        fail("controller Git before commit differs")
    validate_git_receipt(
        value["git_before"], value["git_before"].get("head", ""),
        "controller Git before",
    )
    validate_controller_source_git_binding(
        value["controller_before"], value["git_before"],
        value["source_manifest_before"], "controller before",
    )
    for name, validator in (
        ("build_after", validate_build_manifest),
        ("controller_after", validate_controller_provenance),
        ("source_manifest_after", validate_source_manifest_receipt),
    ):
        if value[name]:
            validator(value[name], "controller " + name.replace("_", " "))
    if value["git_after"]:
        validate_git_receipt(
            value["git_after"], value["git_before"]["head"],
            "controller Git after",
        )
    for before_name, after_name in (
        ("build_before", "build_after"),
        ("controller_before", "controller_after"),
        ("git_before", "git_after"),
        ("source_manifest_before", "source_manifest_after"),
    ):
        if value[after_name] and not canonical_equal(
            value[before_name], value[after_name]
        ):
            fail(
                "controller nonempty final provenance changed: "
                + after_name
            )
    if value["sampler_after"]:
        validate_sampler_receipt(value["sampler_after"])
        if (
            value["sampler_after"]["window_start_monotonic_ns"]
            < controller_started_ns
            or execution_start is not None
            and value["sampler_after"]["window_start_monotonic_ns"]
            > execution_start
        ):
            fail("controller sampler start escapes admission/execution order")
        if (
            value["sampler_after"]["window_end_monotonic_ns"]
            > controller_observed_ns
        ):
            fail("controller observation predates its sampler endpoint")
        if (
            execution_finish is not None
            and value["sampler_after"]["window_end_monotonic_ns"]
            > execution_finish
        ):
            fail("controller sampler endpoint follows binary execution finish")
    child_start_observed = value["binary"]["child_start_monotonic_ns"]
    child_reap_observed = value["binary"]["child_reap_monotonic_ns"]
    if (
        child_start_observed is not None
        and child_start_observed < controller_started_ns
    ):
        fail("controller child start predates its origin")
    if (
        child_reap_observed is not None
        and child_reap_observed > controller_observed_ns
    ):
        fail("controller observation predates its child reap")
    if (
        child_reap_observed is not None
        and reap_reached_controller_deadline(
            child_reap_observed, controller_started_ns
        )
        and value["binary"]["timed_out"] is not True
    ):
        fail("late child reap lacks its inclusive timeout classification")
    if (
        value["controller_after"] and value["git_after"]
        and value["source_manifest_after"]
    ):
        validate_controller_source_git_binding(
            value["controller_after"], value["git_after"],
            value["source_manifest_after"], "controller after",
        )
    if summary["outcome"] != value["outcome"]:
        fail("controller/summary outcomes differ")
    if value["outcome"] in ("pass", "reject"):
        expected_failure = (
            "statistical-gates" if value["outcome"] == "reject" else "none"
        )
        if (
            signals
            or value["failure"] != expected_failure
            or not value["build_after"]
            or not value["controller_after"]
            or not value["git_after"]
            or not value["sampler_after"]
            or not value["source_manifest_after"]
            or not canonical_equal(value["build_before"], value["build_after"])
            or not canonical_equal(
                value["controller_before"], value["controller_after"]
            )
            or not canonical_equal(value["git_before"], value["git_after"])
            or not canonical_equal(
                value["source_manifest_before"], value["source_manifest_after"]
            )
        ):
            fail("positive controller decision lacks exact final provenance")
        child_start = value["binary"]["child_start_monotonic_ns"]
        child_reap = value["binary"]["child_reap_monotonic_ns"]
        if (
            child_start < value["controller_started_monotonic_ns"]
            or reap_reached_controller_deadline(
                child_reap, value["controller_started_monotonic_ns"]
            )
            or value["sampler_after"]["window_end_monotonic_ns"] > (
                value["controller_started_monotonic_ns"]
                + int(
                    HEALTH_COLLECTION_DEADLINE_SECONDS * 1_000_000_000
                )
            )
        ):
            fail("positive controller decision violates a global deadline")
    digest = lower_hash(value["receipt_sha256"], "controller receipt SHA-256")
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail("controller receipt digest differs")
    return value


def parse_complete_document(data: bytes) -> Dict[str, Any]:
    if (
        not data or len(data) > 4096 or not data.endswith(b"\n")
        or data.count(b"\n") != 1 or b"\r" in data
    ):
        fail("COMPLETE framing differs")
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(f"COMPLETE is malformed: {exc}")
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail("COMPLETE is not canonical JSON")
    exact_keys(value, {
        "campaign", "claim_sha256", "controller_outcome",
        "controller_sha256", "complete_observed_monotonic_ns",
        "elapsed_seconds_before_commit", "schema", "status",
    }, "COMPLETE")
    exact_string(value["campaign"], CAMPAIGN, "COMPLETE campaign")
    exact_string(value["schema"], COMPLETE_SCHEMA, "COMPLETE schema")
    exact_string(value["status"], "complete", "COMPLETE status")
    if value["controller_outcome"] not in ("pass", "reject", "invalid"):
        fail("COMPLETE controller outcome differs")
    lower_hash(value["claim_sha256"], "COMPLETE claim SHA-256")
    lower_hash(value["controller_sha256"], "COMPLETE controller SHA-256")
    exact_int(
        value["complete_observed_monotonic_ns"], 0, MAX_INT63,
        "COMPLETE observation timestamp",
    )
    elapsed = value["elapsed_seconds_before_commit"]
    if (
        type(elapsed) is not float or not math.isfinite(elapsed)
        or not 0.0 <= elapsed < OUTER_DEADLINE_SECONDS
    ):
        fail("COMPLETE elapsed time differs")
    return value


def health_thermal_artifact_bytes(health: Mapping[str, Any]) -> bytes:
    thermal = health.get("thermal", {})
    terminal = health.get("sampler_terminal", {})
    if thermal and terminal:
        fail("health cannot bind both live and terminal thermal artifacts")
    text = (
        thermal.get("window_csv_ascii") if thermal
        else terminal.get("window_csv_ascii") if terminal
        else ""
    )
    if type(text) is not str:
        fail("health thermal artifact is not text")
    try:
        data = text.encode("ascii")
    except UnicodeEncodeError:
        fail("health thermal artifact is not ASCII")
    if len(data) > MAX_THERMAL_WINDOW_BYTES:
        fail("health thermal artifact exceeds its bound")
    return data


def validate_controller_bundle_binding(
    controller: Mapping[str, Any], summary: Mapping[str, Any],
    summary_bytes: bytes, claim_receipt: Mapping[str, Any],
    binary: Mapping[str, Any],
) -> None:
    sampler = summary["health"].get("sampler", {})
    admission = summary["health"].get("admission_sibling_ticks", {})
    thermal_receipt = summary["health"].get("thermal", {})
    claimed_collection = summary["health"]["collection_failures"]
    claimed_collection_exact = {
        "run sibling ticks unavailable",
        "sampler start endpoint unavailable",
        "health collection reached the strict 115-second deadline",
    }
    claimed_collection_prefixes = (
        "controller affinity end: ",
        "run sibling ticks: ",
        "thermal collection: ",
        "sampler endpoint: ",
        "sampler stale evidence: ",
        "sampler terminal evidence: ",
        "health terminal collection: ",
    )
    for message in claimed_collection:
        if message in claimed_collection_exact or any(
            message.startswith(prefix) and len(message) > len(prefix)
            for prefix in claimed_collection_prefixes
        ):
            continue
        fail("claimed bundle has impossible health collection evidence")
    for field in (
        "admission_sibling_ticks", "target_core", "controller_core",
        "sampler_core", "target_threads", "controller_initial_affinity",
    ):
        if not summary["health"].get(field):
            fail(f"claimed bundle lacks ready preclaim health field {field}")
    loader = summary["health_module_loader"]
    source_before = summary["source_manifest_before"]
    if not loader or not source_before:
        fail("claimed bundle lacks its preclaim sealed source loader")
    source_by_path = {
        entry["path"]: entry for entry in source_before["entries"]
    }
    for module in loader["modules"]:
        source_entry = source_by_path[module["path"]]
        if (
            module["bytes"] != source_entry["bytes"]
            or module["sha256"] != source_entry["sha256"]
        ):
            fail("claimed bundle health loader/source binding differs")
    thermal_bytes = health_thermal_artifact_bytes(summary["health"])
    artifacts = controller["artifacts"]
    execution_start = binary.get("execution_started_monotonic_ns")
    execution_finish = binary.get("execution_finished_monotonic_ns")
    if (
        controller["outcome"] != summary["outcome"]
        or controller["summary"]["outcome"] != summary["outcome"]
        or controller["target_cpu"] != summary["target_cpu"]
        or not canonical_equal(controller["signals"], summary["signal_names"])
        or controller["summary"]["sha256"] != sha256_bytes(summary_bytes)
        or controller["summary"]["preimage_sha256"]
        != summary["summary_preimage_sha256"]
        or controller["health_receipt_sha256"]
        != summary["health"]["receipt_sha256"]
        or not canonical_equal(controller["binary"], binary)
        or not canonical_equal(controller["claim"], claim_receipt)
        or controller["controller_elapsed_seconds"] < summary["elapsed_seconds"]
        or artifacts["raw.jsonl"]["bytes"] != summary["raw_bytes"]
        or artifacts["raw.jsonl"]["sha256"] != summary["raw_sha256"]
        or artifacts["stderr.txt"]["bytes"] != summary["stderr_bytes"]
        or artifacts["stderr.txt"]["sha256"] != summary["stderr_sha256"]
        or artifacts["summary.json"]["bytes"] != len(summary_bytes)
        or artifacts["summary.json"]["sha256"] != sha256_bytes(summary_bytes)
        or artifacts["thermal.csv"]["bytes"] != len(thermal_bytes)
        or artifacts["thermal.csv"]["sha256"] != sha256_bytes(thermal_bytes)
    ):
        fail("controller/summary/claim evidence binding differs")
    if controller["outcome"] == "invalid" and (
        not summary["publication_failures"]
        or summary["failure"] in ("none", "statistical-gates")
        or controller["failure"] in ("none", "statistical-gates")
    ):
        fail("invalid controller bundle lacks its causal failure evidence")
    validate_execution_reap_binding(
        binary, summary["process_exit_code"], "controller bundle execution"
    )
    child_pid = binary.get("child_pid")
    if child_pid is not None and (
        child_pid == controller["controller_before"]["pid"]
        or sampler and child_pid == sampler["pid"]
    ):
        fail("controller/worker/sampler process identities overlap")
    if controller["sampler_after"]:
        if not sampler:
            fail("controller sampler lacks its health sampler binding")
        validate_sampler_growth_binding(
            sampler, controller["sampler_after"],
            "controller/summary sampler",
        )
    if (
        sampler
        and sampler["window_end_monotonic_ns"]
        > controller["controller_observed_monotonic_ns"]
    ):
        fail("controller observation predates its health sampler endpoint")
    if sampler and (
        sampler["window_start_monotonic_ns"]
        < controller["controller_started_monotonic_ns"]
        or execution_start is not None
        and sampler["window_start_monotonic_ns"] > execution_start
    ):
        fail("health sampler start escapes controller/execution order")
    if admission:
        admission_start = admission["interval_start_monotonic_ns"]
        admission_end = admission["interval_end_monotonic_ns"]
        if (
            admission_start != admission["start"]["read_monotonic_ns"]
            or admission_end != admission["end"]["read_monotonic_ns"]
            or admission_start < controller["controller_started_monotonic_ns"]
            or admission_end >= (
                controller["controller_started_monotonic_ns"]
                + 10_000_000_000
            )
            or admission_end > controller["controller_observed_monotonic_ns"]
            or execution_start is not None and admission_end > execution_start
            or sampler
            and sampler["window_start_monotonic_ns"] > admission_end
        ):
            fail("health admission chronology escapes its preclaim interval")
    if (
        sampler and execution_finish is not None
        and sampler["window_end_monotonic_ns"] > execution_finish
    ):
        fail("health sampler endpoint follows binary execution finish")
    if controller["outcome"] in ("pass", "reject") and (
        not sampler
        or sampler["window_end_monotonic_ns"] > (
            controller["controller_started_monotonic_ns"]
            + int(HEALTH_COLLECTION_DEADLINE_SECONDS * 1_000_000_000)
        )
    ):
        fail("positive controller health reached the collection deadline")


def validate_complete_bundle_binding(
    complete: Mapping[str, Any], controller: Mapping[str, Any],
    controller_bytes: bytes, claim_receipt: Mapping[str, Any],
) -> None:
    controller_started_ns = controller["controller_started_monotonic_ns"]
    complete_observed_ns = complete["complete_observed_monotonic_ns"]
    expected_elapsed = (
        complete_observed_ns - controller_started_ns
    ) / 1_000_000_000.0
    if (
        complete["claim_sha256"] != claim_receipt["sha256"]
        or complete["controller_sha256"] != sha256_bytes(controller_bytes)
        or complete["controller_outcome"] != controller["outcome"]
        or complete_observed_ns
        < controller["controller_observed_monotonic_ns"]
        or complete_observed_ns >= (
            controller_started_ns
            + int(OUTER_DEADLINE_SECONDS * 1_000_000_000)
        )
        or complete["elapsed_seconds_before_commit"] != expected_elapsed
    ):
        fail("COMPLETE/controller/claim binding differs")


def collect_final_attestations(
    args: argparse.Namespace, claim: ClaimReservation,
    controller_before: Mapping[str, Any], build_before: Mapping[str, Any],
    source_before: Mapping[str, Any], git_before: Mapping[str, Any],
    execution: ExecutionResult, source_root: Path,
    health_modules: SealedHealthModules, errors: List[str],
    deadline: float,
) -> Dict[str, Any]:
    values: Dict[str, Any] = {
        "binary": {}, "build": {}, "controller": {}, "git": {},
        "sampler": {}, "source": {},
    }
    try:
        claim.verify()
    except BaseException as exc:
        append_unique_failure(errors, "claim final check: " + exception_text(exc))
    try:
        source_first = source_manifest(source_root, deadline)
        git_current = git_receipt(
            source_root, args.expected_source_commit,
            args.expected_git_sha256,
            deadline,
        )
        source_second = source_manifest(source_root, deadline)
        if (
            not canonical_equal(source_first, source_second)
            or not canonical_equal(source_second, source_before)
            or source_second["sha256"] != args.expected_source_manifest_sha256
            or not canonical_equal(git_current, git_before)
        ):
            fail("final source/Git/source sandwich differs")
        values["source"] = source_second
        values["git"] = git_current
    except BaseException as exc:
        append_unique_failure(errors, "source/Git final check: " + exception_text(exc))
    try:
        current_build = build_manifest(args.build_dir, deadline)
        if (
            current_build["sha256"] != args.expected_build_manifest_sha256
            or not canonical_equal(current_build, build_before)
        ):
            fail("final build manifest differs")
        values["build"] = current_build
    except BaseException as exc:
        append_unique_failure(errors, "build final check: " + exception_text(exc))
    try:
        current_controller = controller_provenance(args, deadline)
        if not canonical_equal(current_controller, controller_before):
            fail("final controller provenance differs")
        values["controller"] = current_controller
    except BaseException as exc:
        append_unique_failure(errors, "controller final check: " + exception_text(exc))
    try:
        binary_hash, binary_info = hash_regular_path(
            Path(args.binary), "final binary", max_size=MAX_BINARY_BYTES,
            deadline=deadline,
        )
        current_binary = stat_receipt(binary_info)
        if (
            binary_hash != args.expected_binary_sha256
            or binary_info.st_uid != args.expected_binary_uid
            or binary_info.st_gid != args.expected_controller_gid
            or binary_info.st_nlink != 1
            or stat.S_IMODE(binary_info.st_mode) != 0o555
            or not canonical_equal(
                current_binary, execution.binary.get("path_stat_after", {})
            )
        ):
            fail("final binary image differs")
        values["binary"] = {
            "sha256": binary_hash, "stat": current_binary,
        }
    except BaseException as exc:
        append_unique_failure(errors, "binary final check: " + exception_text(exc))
    sampler = execution.health.get("sampler", {})
    if sampler:
        try:
            current_sampler = make_sampler_attestation(
                args, sampler["window_start_monotonic_ns"],
                sampler["window_end_monotonic_ns"], health_modules,
                deadline,
            )
            validate_sampler_growth_binding(
                sampler, current_sampler, "final sampler attestation"
            )
            thermal, thermal_collection, _ = tolerant_thermal_window(
                args, current_sampler,
                sampler["window_start_monotonic_ns"],
                sampler["window_end_monotonic_ns"],
            )
            if thermal_collection or not canonical_equal(
                thermal, execution.health.get("thermal", {})
            ):
                fail("final sampler thermal window changed in place")
            values["sampler"] = current_sampler
        except BaseException as exc:
            append_unique_failure(
                errors, "sampler final check: " + exception_text(exc)
            )
    return values


def validate_final_attestation_growth_binding(
    bound: Mapping[str, Any], current: Mapping[str, Any],
) -> None:
    expected = {"binary", "build", "controller", "git", "sampler", "source"}
    exact_keys(bound, expected, "bound final attestations")
    exact_keys(current, expected, "current final attestations")
    for name in expected - {"sampler"}:
        if not canonical_equal(bound[name], current[name]):
            fail(f"final {name} attestation changed after controller binding")
    if bound["sampler"] or current["sampler"]:
        if not bound["sampler"] or not current["sampler"]:
            fail("final sampler attestation disappeared")
        validate_sampler_growth_binding(
            bound["sampler"], current["sampler"],
            "final fixed-point sampler",
        )


def verify_final_artifacts(
    bundle: OutputBundle, held: Mapping[str, int],
    expected: Mapping[str, Mapping[str, Any]], complete_mode: int,
) -> Dict[str, Dict[str, Any]]:
    bundle._verify_identities(0o500, 0o400)
    bundle._verify_roster(FINAL_OUTPUT_NAMES)
    current: Dict[str, Dict[str, Any]] = {}
    for name in FINAL_OUTPUT_NAMES:
        mode = complete_mode if name == "COMPLETE" else 0o400
        receipt = held_artifact(held[name], mode)
        named = os.stat(name, dir_fd=bundle.directory_fd, follow_symlinks=False)
        retained = os.fstat(held[name])
        if not same_file_receipt(retained, named):
            fail(f"final artifact {name} named/held identity differs")
        prior = expected.get(name)
        if prior is not None and not canonical_equal(receipt, prior):
            fail(f"final artifact {name} changed after binding")
        current[name] = receipt
    return current


def finalize_complete_marker(
    fd: int, directory_fd: int, complete_bytes: bytes,
    controller_started_ns: int,
) -> None:
    """Expose COMPLETE, or force it back to an unambiguous non-marker."""
    try:
        current = OutputBundle._read_current(fd, 4096)
        if current != complete_bytes or stat.S_IMODE(os.fstat(fd).st_mode) != 0:
            fail("prepared COMPLETE marker differs")
        if monotonic_deadline_reached(
            time.monotonic_ns(), controller_started_ns,
            OUTER_DEADLINE_SECONDS,
        ):
            fail("COMPLETE semantic chmod reached the strict outer deadline")
        # This is the last enforceable userspace observation before exposure.
        # Kernel scheduling between this check and the one semantic fchmod is
        # part of the trusted syscall-linearization boundary; no post-exposure
        # timer or status check may create a contradictory result.
        # The bytes/dentry were already fsynced and read back while mode000.
        # This chmod is the sole visible semantic commit.  A later fsync could
        # fail after exposing a valid-looking marker, so none is permitted.
        os.fchmod(fd, 0o400)
        return
    except BaseException as original:
        reset_error: Optional[BaseException] = None
        try:
            os.fchmod(fd, 0o000)
            os.ftruncate(fd, 0)
            os.fsync(fd)
            if (
                OutputBundle._read_current(fd, 1) != b""
                or stat.S_IMODE(os.fstat(fd).st_mode) != 0
            ):
                fail("COMPLETE reset readback differs")
        except BaseException as exc:
            reset_error = exc
            try:
                retained = os.fstat(fd)
                named = os.stat(
                    "COMPLETE", dir_fd=directory_fd, follow_symlinks=False
                )
                if not same_file_receipt(retained, named):
                    fail("COMPLETE fallback unlink identity differs")
                os.fchmod(directory_fd, 0o700)
                os.unlink("COMPLETE", dir_fd=directory_fd)
                os.fsync(directory_fd)
                os.fchmod(directory_fd, 0o500)
                os.fsync(directory_fd)
                reset_error = None
            except BaseException as unlink_exc:
                reset_error = unlink_exc
        if reset_error is not None:
            fail(
                "COMPLETE exposure failed and reset failed: "
                + exception_text(original) + "; " + exception_text(reset_error)
            )
        raise original


def open_readonly_final_evidence(
    bundle: OutputBundle, held: Mapping[str, int], claim: ClaimReservation,
    expected: Mapping[str, Mapping[str, Any]],
) -> Dict[str, int]:
    mirrors: Dict[str, int] = {}
    flags = nonblocking_read_flags(
        "read-only final evidence", noatime=True
    )
    try:
        for name in (
            "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
            "controller.json",
        ):
            fd = os.open(name, flags, dir_fd=bundle.directory_fd)
            mirrors[name] = fd
            retained = os.fstat(held[name])
            readonly = os.fstat(fd)
            receipt = held_artifact(fd, 0o400)
            if (
                not same_file_receipt(retained, readonly)
                or not canonical_equal(receipt, expected[name])
            ):
                fail(f"read-only final evidence mirror differs for {name}")
        claim_fd = os.open(
            FIXED_CLAIM_PATH.name, flags, dir_fd=claim.parent_fd
        )
        mirrors["claim"] = claim_fd
        claim_retained = os.fstat(claim.fd)
        claim_readonly = os.fstat(claim_fd)
        if (
            not same_file_receipt(claim_retained, claim_readonly)
            or claim_readonly.st_size != len(claim.raw)
            or file_sha256_fd(claim_fd, claim_readonly.st_size)
            != sha256_bytes(claim.raw)
            or stat.S_IMODE(claim_readonly.st_mode) != 0o400
        ):
            fail("read-only final claim mirror differs")
        return mirrors
    except BaseException:
        for fd in mirrors.values():
            try:
                os.close(fd)
            except OSError:
                pass
        raise


def verify_surviving_named_artifact(
    directory_fd: int, fd: int, name: str,
    expected: Mapping[str, Any], expected_mode: int,
    expected_uid: int, expected_gid: int,
) -> Dict[str, Any]:
    before = os.fstat(fd)
    receipt = held_artifact(fd, expected_mode)
    after = os.fstat(fd)
    named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        not same_file_receipt(before, after)
        or not same_file_receipt(after, named)
        or not canonical_equal(receipt, expected)
        or receipt["uid"] != expected_uid
        or receipt["gid"] != expected_gid
    ):
        fail(f"final surviving artifact authority differs for {name}")
    return receipt


def verify_surviving_final_authorities(
    bundle: OutputBundle, mirrors: Mapping[str, int], complete_fd: int,
    expected: Mapping[str, Mapping[str, Any]], claim: ClaimReservation,
    claim_receipt: Mapping[str, Any],
) -> None:
    """Recheck only the descriptors that survive the writable-FD drain."""
    evidence_names = (
        "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
        "controller.json",
    )
    exact_keys(mirrors, (*evidence_names, "claim"), "final read-only mirrors")
    exact_keys(expected, FINAL_OUTPUT_NAMES, "final artifact snapshot")
    validate_claim_receipt(claim_receipt, "final bound claim receipt")

    bundle._verify_parent_path()
    output_parent = os.fstat(bundle.parent_fd)
    claim_parent = os.fstat(claim.parent_fd)
    parent_named = os.stat("/var/tmp", follow_symlinks=False)
    expected_parent = claim_receipt["parent"]
    for info in (output_parent, claim_parent, parent_named):
        if (
            not stat.S_ISDIR(info.st_mode)
            or info.st_dev != expected_parent["device"]
            or info.st_ino != expected_parent["inode"]
            or stat.S_IMODE(info.st_mode) != expected_parent["mode"]
            or info.st_uid != expected_parent["uid"]
            or info.st_gid != expected_parent["gid"]
        ):
            fail("final surviving /var/tmp authority differs")

    directory_info = os.fstat(bundle.directory_fd)
    directory_named = os.stat(
        FIXED_OUTPUT_DIR.name, dir_fd=bundle.parent_fd,
        follow_symlinks=False,
    )
    if (
        not stat.S_ISDIR(directory_info.st_mode)
        or not same_file_receipt(directory_info, directory_named)
        or (
            directory_info.st_dev, directory_info.st_ino,
            directory_info.st_nlink,
        ) != bundle.directory_identity
        or stat.S_IMODE(directory_info.st_mode) != 0o500
        or directory_info.st_uid != claim_receipt["stat"]["uid"]
        or directory_info.st_gid != claim_receipt["stat"]["gid"]
    ):
        fail("final surviving output-directory authority differs")
    bundle._verify_roster(FINAL_OUTPUT_NAMES)

    for name in evidence_names:
        verify_surviving_named_artifact(
            bundle.directory_fd, mirrors[name], name, expected[name], 0o400,
            claim_receipt["stat"]["uid"], claim_receipt["stat"]["gid"],
        )
    verify_surviving_named_artifact(
        bundle.directory_fd, complete_fd, "COMPLETE", expected["COMPLETE"],
        0o000, claim_receipt["stat"]["uid"],
        claim_receipt["stat"]["gid"],
    )

    claim_fd = mirrors["claim"]
    claim_before = os.fstat(claim_fd)
    if claim_before.st_size != claim_receipt["bytes"]:
        fail("final surviving claim size differs")
    claim_hash = file_sha256_fd(claim_fd, claim_before.st_size)
    claim_after = os.fstat(claim_fd)
    claim_named = os.stat(
        FIXED_CLAIM_PATH.name, dir_fd=claim.parent_fd,
        follow_symlinks=False,
    )
    if (
        not same_file_receipt(claim_before, claim_after)
        or not same_file_receipt(claim_after, claim_named)
        or not canonical_equal(stat_receipt(claim_after), claim_receipt["stat"])
        or claim_hash != claim_receipt["sha256"]
    ):
        fail("final surviving claim authority differs")


def poison_noncommitted_publication(
    bundle: OutputBundle, held: Mapping[str, int], errors: Sequence[str],
    guarded_signals: SignalGuard,
    base_summary: Optional[Mapping[str, Any]],
    summary_identity_fd: int = -1,
) -> None:
    complete_fd = held.get("COMPLETE", -1)
    if complete_fd >= 0:
        try:
            os.fchmod(complete_fd, 0o000)
            os.ftruncate(complete_fd, 0)
            os.fsync(complete_fd)
            if (
                OutputBundle._read_current(complete_fd, 1) != b""
                or stat.S_IMODE(os.fstat(complete_fd).st_mode) != 0
            ):
                fail("noncommitted COMPLETE poison differs")
        except BaseException:
            try:
                retained = os.fstat(complete_fd)
                named = os.stat(
                    "COMPLETE", dir_fd=bundle.directory_fd,
                    follow_symlinks=False,
                )
                if not same_file_receipt(retained, named):
                    fail("noncommitted COMPLETE unlink identity differs")
                os.fchmod(bundle.directory_fd, 0o700)
                os.unlink("COMPLETE", dir_fd=bundle.directory_fd)
                os.fsync(bundle.directory_fd)
                os.fchmod(bundle.directory_fd, 0o500)
                os.fsync(bundle.directory_fd)
            except BaseException:
                pass

    summary_fd = held.get("summary.json", -1)
    writable_owned = False
    identity_fd = summary_identity_fd
    identity_owned = False
    try:
        expected = bundle.file_identities["summary.json"]
        if identity_fd < 0:
            identity_fd = os.open(
                "summary.json",
                nonblocking_read_flags("summary poison identity"),
                dir_fd=bundle.directory_fd,
            )
            identity_owned = True
        identity_info = os.fstat(identity_fd)
        named_path = os.stat(
            "summary.json", dir_fd=bundle.directory_fd,
            follow_symlinks=False,
        )
        if (
            identity_info.st_dev, identity_info.st_ino, identity_info.st_nlink
        ) != expected or not same_file_receipt(identity_info, named_path):
            fail("summary poison named identity differs")
        try:
            retained_info = os.fstat(summary_fd)
            if not same_file_receipt(retained_info, identity_info):
                fail("summary poison retained identity differs")
        except OSError:
            # Keep this bound read-only descriptor live through the writable
            # reopen.  If either chmod or reopen faults, it remains an
            # identity-capable handle for mode000 poisoning or safe unlink.
            os.fchmod(identity_fd, 0o600)
            summary_fd = os.open(
                "summary.json",
                os.O_RDWR | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0)
                | getattr(os, "O_NONBLOCK", 0),
                dir_fd=bundle.directory_fd,
            )
            writable_owned = True
            if not same_file_receipt(os.fstat(summary_fd), identity_info):
                fail("summary poison writable reopen identity differs")
        os.fchmod(summary_fd, 0o600)
        if base_summary is None:
            fail("final abort lacks a summary preimage")
        abort_errors = list(errors)
        append_unique_failure(abort_errors, "final publication aborted")
        invalid_bytes = bundle._invalid_summary(
            base_summary, abort_errors, guarded_signals
        )
        parsed = parse_summary_bytes(
            OutputBundle._write_exact(summary_fd, invalid_bytes)
        )
        if parsed["outcome"] != "invalid":
            fail("final abort summary is not invalid")
        os.fchmod(summary_fd, 0o400)
        os.fsync(summary_fd)
    except BaseException:
        poison_fd = -1
        for candidate in (summary_fd, identity_fd):
            if candidate < 0:
                continue
            try:
                candidate_info = os.fstat(candidate)
            except OSError:
                continue
            if (
                candidate_info.st_dev, candidate_info.st_ino,
                candidate_info.st_nlink,
            ) == bundle.file_identities["summary.json"]:
                poison_fd = candidate
                break
        try:
            if poison_fd < 0:
                fail("summary poison lost every bound identity descriptor")
            os.fchmod(poison_fd, 0o000)
            try:
                os.ftruncate(poison_fd, 0)
                os.fsync(poison_fd)
            except OSError:
                # O_RDONLY identity handles cannot truncate, but mode000 is an
                # authoritative non-summary state.  Prefer truncation when a
                # writable handle survived; otherwise retain the safe mode.
                if stat.S_IMODE(os.fstat(poison_fd).st_mode) != 0:
                    raise
        except BaseException:
            try:
                if poison_fd < 0:
                    fail("summary poison unlink lacks a bound descriptor")
                retained = os.fstat(poison_fd)
                named = os.stat(
                    "summary.json", dir_fd=bundle.directory_fd,
                    follow_symlinks=False,
                )
                if not same_file_receipt(retained, named):
                    fail("summary poison unlink identity differs")
                os.fchmod(bundle.directory_fd, 0o700)
                os.unlink("summary.json", dir_fd=bundle.directory_fd)
                os.fsync(bundle.directory_fd)
                os.fchmod(bundle.directory_fd, 0o500)
                os.fsync(bundle.directory_fd)
            except BaseException:
                pass
    finally:
        if writable_owned:
            try:
                os.close(summary_fd)
            except OSError:
                pass
        if identity_owned and identity_fd >= 0:
            try:
                os.close(identity_fd)
            except OSError:
                pass


def final_publish(
    args: argparse.Namespace, bundle: OutputBundle, claim: ClaimReservation,
    guarded_signals: SignalGuard, controller_started: float,
    controller_started_ns: int, proposed_outcome: str, errors: List[str],
    execution: ExecutionResult, source_before: Mapping[str, Any],
    git_before: Mapping[str, Any], build_before: Mapping[str, Any],
    controller_before: Mapping[str, Any], source_root: Path,
    health_modules: SealedHealthModules, base_summary_bytes: bytes,
) -> Tuple[bool, str, List[str]]:
    if proposed_outcome not in ("pass", "reject", "invalid"):
        fail("proposed controller outcome differs")
    held: Dict[str, int] = dict(bundle.file_fds)
    complete_created = False
    committed = False
    base_summary: Optional[Dict[str, Any]] = None
    readonly_mirrors: Dict[str, int] = {}
    try:
        validate_execution_reap_binding(
            execution.binary, execution.returncode,
            "final publication execution",
        )
        if monotonic_deadline_reached(
            time.monotonic_ns(), controller_started_ns,
            FINAL_COMMIT_START_DEADLINE_SECONDS,
        ):
            append_unique_failure(
                errors,
                "final publication did not start before 117 seconds",
            )
        claim_document = parse_claim_document(claim.raw)
        expected_claim_document = {
            "binary_sha256": args.expected_binary_sha256,
            "build_manifest_sha256": args.expected_build_manifest_sha256,
            "campaign": CAMPAIGN,
            "controller_receipt_sha256": controller_before["receipt_sha256"],
            "controller_started_monotonic_ns": controller_started_ns,
            "gid": os.getgid(),
            "git_sha256": args.expected_git_sha256,
            "output_path": str(FIXED_OUTPUT_DIR),
            "parent_device": claim.parent_identity[0],
            "parent_inode": claim.parent_identity[1],
            "pid": os.getpid(),
            "process_start_ticks": controller_before["process_start_ticks"],
            "schema": CLAIM_SCHEMA,
            "source_commit": args.expected_source_commit,
            "source_manifest_sha256": args.expected_source_manifest_sha256,
            "uid": os.getuid(),
        }
        claim_receipt = claim.receipt()
        validate_claim_runtime_binding(
            claim_document, claim_receipt, expected_claim_document, claim.raw
        )
        output_receipt = bundle.receipt()
        if (
            output_receipt["path"] != str(FIXED_OUTPUT_DIR)
            or output_receipt["parent"]["device"]
            != claim_receipt["parent"]["device"]
            or output_receipt["parent"]["inode"]
            != claim_receipt["parent"]["inode"]
        ):
            fail("fixed output/claim parent identity binding differs")
        bundle._verify_identities(0o700, 0o600)
        bundle._verify_roster(OUTPUT_NAMES)
        held["thermal.csv"] = reserve_output_file(
            bundle.directory_fd, "thermal.csv", 0o600
        )
        held["controller.json"] = reserve_output_file(
            bundle.directory_fd, "controller.json", 0o600
        )
        bundle._verify_roster((
            "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
            "controller.json",
        ))
        thermal_bytes = health_thermal_artifact_bytes(execution.health)
        OutputBundle._write_exact(held["thermal.csv"], thermal_bytes)
        base_summary = parse_summary_bytes(base_summary_bytes)
        bundle._validate_summary_cross(base_summary_bytes)
        OutputBundle._write_exact(held["summary.json"], base_summary_bytes)
        guarded_signals.block_for_commit()

        # Harden the five evidence objects before the COMPLETE dentry exists.
        for name in (
            "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
            "controller.json",
        ):
            os.fchmod(held[name], 0o400)
            os.fsync(held[name])

        frozen_result: Optional[Tuple[bool, str, List[str]]] = None
        signal_cutoff_frozen = False
        for _ in range(len(SignalGuard.SIGNALS) + 8):
            if not signal_cutoff_frozen:
                guarded_signals.collect_pending()
            for signal_name in guarded_signals.observed_signals:
                append_unique_failure(errors, "controller signal: " + signal_name)
            if guarded_signals.kill_error is not None:
                append_unique_failure(
                    errors,
                    "signal process-group kill: " + guarded_signals.kill_error,
                )
            attestations = collect_final_attestations(
                args, claim, controller_before, build_before, source_before,
                git_before, execution, source_root, health_modules, errors,
                controller_started + 116.0,
            )
            outcome = "invalid" if errors else proposed_outcome
            if errors or guarded_signals.observed_signals:
                summary_bytes = bundle._invalid_summary(
                    base_summary, errors, guarded_signals
                )
            else:
                summary_bytes = base_summary_bytes
            summary = bundle._validate_summary_cross(
                OutputBundle._write_exact(held["summary.json"], summary_bytes)
            )
            raw_now = OutputBundle._read_current(
                held["raw.jsonl"], MAX_STDOUT_BYTES + 1
            )
            stderr_now = OutputBundle._read_current(
                held["stderr.txt"], MAX_STDERR_BYTES + 1
            )
            if (
                summary["raw_bytes"] != len(raw_now)
                or summary["raw_record_count"] != raw_now.count(b"\n")
                or summary["raw_sha256"] != sha256_bytes(raw_now)
                or summary["stderr_bytes"] != len(stderr_now)
                or summary["stderr_sha256"] != sha256_bytes(stderr_now)
                or thermal_bytes != OutputBundle._read_current(
                    held["thermal.csv"], MAX_THERMAL_WINDOW_BYTES
                )
            ):
                fail("final summary/evidence cross-binding differs")
            artifacts = {
                name: held_artifact(held[name], 0o400)
                for name in (
                    "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv"
                )
            }
            controller_observed_ns = time.monotonic_ns()
            elapsed = (
                controller_observed_ns - controller_started_ns
            ) / 1_000_000_000.0
            document = {
                "artifacts": artifacts,
                "binary": execution.binary,
                "build_after": attestations["build"],
                "build_before": dict(build_before),
                "campaign": CAMPAIGN,
                "claim": claim.receipt(),
                "controller_after": attestations["controller"],
                "controller_before": dict(controller_before),
                "controller_deadline_seconds": CONTROLLER_DEADLINE_SECONDS,
                "controller_elapsed_seconds": elapsed,
                "controller_observed_monotonic_ns": controller_observed_ns,
                "controller_started_monotonic_ns": controller_started_ns,
                "failure": (
                    bounded_failure_text(errors) if errors
                    else ("statistical-gates" if outcome == "reject" else "none")
                ),
                "git_after": attestations["git"],
                "git_before": dict(git_before),
                "health_receipt_sha256": execution.health["receipt_sha256"],
                "outcome": outcome,
                "outer_deadline_seconds": OUTER_DEADLINE_SECONDS,
                "output_path": str(FIXED_OUTPUT_DIR),
                "receipt_sha256": None,
                "sampler_after": attestations["sampler"],
                "schema": CONTROLLER_SCHEMA,
                "signals": list(guarded_signals.observed_signals),
                "source_manifest_after": attestations["source"],
                "source_manifest_before": dict(source_before),
                "summary": {
                    "outcome": summary["outcome"],
                    "preimage_sha256": summary["summary_preimage_sha256"],
                    "sha256": sha256_bytes(summary_bytes),
                },
                "target_cpu": EXPECTED_TARGET_CPU,
            }
            controller_bytes = controller_document_bytes(document)
            parsed_controller = parse_controller_document(
                OutputBundle._write_exact(
                    held["controller.json"], controller_bytes
                )
            )
            validate_controller_bundle_binding(
                parsed_controller, summary, summary_bytes, claim_receipt,
                execution.binary,
            )
            if (
                parsed_controller["summary"]["sha256"]
                != artifacts["summary.json"]["sha256"]
            ):
                fail("controller summary/artifact binding differs")

            if not complete_created:
                bundle._verify_roster((
                    "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
                    "controller.json",
                ))
                held["COMPLETE"] = reserve_output_file(
                    bundle.directory_fd, "COMPLETE", 0o000
                )
                complete_created = True
                os.fsync(bundle.directory_fd)
                os.fchmod(bundle.directory_fd, 0o500)
                os.fsync(bundle.directory_fd)
                os.fsync(bundle.parent_fd)

            artifact_snapshot = dict(artifacts)
            artifact_snapshot["controller.json"] = held_artifact(
                held["controller.json"], 0o400
            )
            artifact_snapshot["COMPLETE"] = held_artifact(
                held["COMPLETE"], 0o000
            )
            verify_final_artifacts(
                bundle, held, artifact_snapshot, complete_mode=0o000
            )

            before_final_errors = tuple(errors)
            before_final_signals = tuple(guarded_signals.observed_signals)
            final_attestations = collect_final_attestations(
                args, claim, controller_before, build_before, source_before,
                git_before, execution, source_root, health_modules, errors,
                controller_started + 116.0,
            )
            try:
                validate_final_attestation_growth_binding(
                    attestations, final_attestations
                )
            except BaseException as exc:
                append_unique_failure(
                    errors,
                    "final attestation changed after controller binding: "
                    + exception_text(exc),
                )
            if not signal_cutoff_frozen:
                guarded_signals.collect_pending()
            for signal_name in guarded_signals.observed_signals:
                append_unique_failure(errors, "controller signal: " + signal_name)
            if (
                tuple(errors) != before_final_errors
                or tuple(guarded_signals.observed_signals) != before_final_signals
            ):
                continue
            complete_observed_ns = time.monotonic_ns()
            if monotonic_deadline_reached(
                complete_observed_ns, controller_started_ns,
                OUTER_DEADLINE_SECONDS,
            ):
                fail("controller exhausted the outer publication deadline")
            complete = {
                "campaign": CAMPAIGN,
                "claim_sha256": claim.receipt()["sha256"],
                "complete_observed_monotonic_ns": complete_observed_ns,
                "controller_outcome": outcome,
                "controller_sha256": sha256_bytes(controller_bytes),
                "elapsed_seconds_before_commit": (
                    (complete_observed_ns - controller_started_ns)
                    / 1_000_000_000.0
                ),
                "schema": COMPLETE_SCHEMA,
                "status": "complete",
            }
            complete_bytes = canonical_bytes(complete) + b"\n"
            parsed_complete = parse_complete_document(
                OutputBundle._write_exact(held["COMPLETE"], complete_bytes)
            )
            validate_complete_bundle_binding(
                parsed_complete, parsed_controller, controller_bytes,
                claim_receipt,
            )
            final_artifact_snapshot = {
                **artifact_snapshot,
                "COMPLETE": held_artifact(held["COMPLETE"], 0o000),
            }
            verify_final_artifacts(
                bundle, held, final_artifact_snapshot,
                complete_mode=0o000,
            )
            readonly_mirrors = open_readonly_final_evidence(
                bundle, held, claim, final_artifact_snapshot,
            )
            # One last fixed-point signal snapshot.  A newly observed signal
            # rewrites summary/controller/COMPLETE while COMPLETE is mode000.
            if (
                not signal_cutoff_frozen
                and guarded_signals.collect_pending()
            ):
                for fd in readonly_mirrors.values():
                    os.close(fd)
                readonly_mirrors = {}
                continue
            if tuple(guarded_signals.observed_signals) != before_final_signals:
                for fd in readonly_mirrors.values():
                    os.close(fd)
                readonly_mirrors = {}
                continue
            if monotonic_deadline_reached(
                time.monotonic_ns(), controller_started_ns,
                OUTER_DEADLINE_SECONDS,
            ):
                fail("controller reached the strict outer semantic cutoff")
            if not signal_cutoff_frozen:
                # Freeze the in-memory decision while the guarded set remains
                # blocked, then conservatively drain once across that cutoff.
                # A signal found by the drain is folded into a newly written
                # invalid summary/controller/COMPLETE on the next iteration.
                # Signals arriving later are post-cutoff and remain blocked
                # through the controller's immediate exit.
                cutoff_changed = guarded_signals.commit_logical_seal()
                signal_cutoff_frozen = True
                if cutoff_changed:
                    for fd in readonly_mirrors.values():
                        os.close(fd)
                    readonly_mirrors = {}
                    for signal_name in guarded_signals.observed_signals:
                        append_unique_failure(
                            errors,
                            "controller signal at cutoff: " + signal_name,
                        )
                    continue
            # No writable descriptor for evidence or the claim survives the
            # semantic commit.  COMPLETE alone remains writable for its final
            # mode transition; all mirrors were already rehashed O_RDONLY.
            for name in (
                "raw.jsonl", "stderr.txt", "summary.json", "thermal.csv",
                "controller.json",
            ):
                os.close(held.pop(name))
            os.close(claim.fd)
            verify_surviving_final_authorities(
                bundle, readonly_mirrors, held["COMPLETE"],
                final_artifact_snapshot, claim, claim_receipt,
            )
            frozen_result = (True, outcome, errors)
            # Sole visible semantic transition.  On success the guard remains
            # blocked until immediate process exit; no status work follows.
            finalize_complete_marker(
                held["COMPLETE"], bundle.directory_fd, complete_bytes,
                controller_started_ns,
            )
            bundle.closed = True
            guarded_signals.evidence_committed = True
            committed = True
            return frozen_result
        fail("final signal/attestation fixed point did not converge")
    except BaseException as exc:
        append_unique_failure(errors, "final publication: " + exception_text(exc))
        return False, "invalid", errors
    finally:
        if not committed:
            poison_noncommitted_publication(
                bundle, held, errors, guarded_signals, base_summary,
                readonly_mirrors.get("summary.json", -1),
            )
            for fd in readonly_mirrors.values():
                try:
                    os.close(fd)
                except OSError:
                    pass
            bundle.closed = True
            for name, fd in tuple(held.items()):
                try:
                    if name not in ("COMPLETE", "summary.json"):
                        os.fchmod(fd, 0o400)
                        os.fsync(fd)
                except OSError:
                    pass
                try:
                    os.close(fd)
                except OSError:
                    pass
            try:
                os.fchmod(bundle.directory_fd, 0o500)
                os.fsync(bundle.directory_fd)
                os.fsync(bundle.parent_fd)
            except OSError:
                pass
            for fd in (bundle.directory_fd, bundle.parent_fd):
                try:
                    os.close(fd)
                except OSError:
                    pass
            bundle.file_fds.clear()


def require_fixed_namespace_absent() -> None:
    for path, where in (
        (FIXED_CLAIM_PATH, "fixed campaign claim"),
        (FIXED_OUTPUT_DIR, "fixed campaign output"),
    ):
        try:
            os.stat(str(path), follow_symlinks=False)
        except FileNotFoundError:
            continue
        fail(f"{where} already exists")


def preflight_binary(
    args: argparse.Namespace, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    binary = Path(args.binary)
    if binary.resolve(strict=True) != binary:
        fail("screen binary path is not canonical")
    digest, info = hash_regular_path(
        binary, "screen binary preflight", max_size=MAX_BINARY_BYTES,
        deadline=deadline,
    )
    if (
        digest != args.expected_binary_sha256
        or info.st_uid != args.expected_binary_uid
        or info.st_gid != args.expected_controller_gid
        or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != 0o555
    ):
        fail("screen binary preflight seal differs")
    return {"sha256": digest, "stat": stat_receipt(info)}


PREFLIGHT_BINARY_KEYS = {"path", "sha256", "stat"}


def observe_preflight_binary(
    path: Path, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    if path.resolve(strict=True) != path:
        fail("preflight screen binary path is not canonical")
    digest, info = hash_regular_path(
        path, "preflight screen binary", max_size=MAX_BINARY_BYTES,
        deadline=deadline,
    )
    if (
        not 1 <= info.st_size <= MAX_BINARY_BYTES
        or info.st_uid != EXPECTED_CAMPAIGN_UID
        or info.st_gid != EXPECTED_CAMPAIGN_GID
        or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != 0o555
    ):
        fail("preflight screen binary seal differs")
    return {
        "path": str(path), "sha256": digest, "stat": stat_receipt(info),
    }


def validate_preflight_binary_observation(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, PREFLIGHT_BINARY_KEYS, where)
    exact_absolute_path(value["path"], where + ".path")
    lower_hash(value["sha256"], where + ".sha256")
    validate_stat_receipt(value["stat"], where + ".stat")
    if (
        value["stat"]["mode"] != 0o555
        or value["stat"]["uid"] != EXPECTED_CAMPAIGN_UID
        or value["stat"]["gid"] != EXPECTED_CAMPAIGN_GID
        or value["stat"]["nlink"] != 1
        or not 1 <= value["stat"]["size"] <= MAX_BINARY_BYTES
    ):
        fail(f"{where} binary policy differs")


def build_run_authority_argv(
    values: Mapping[str, str], controller_pid: int,
) -> List[str]:
    if type(values) is not dict:
        fail("preflight run authority values are not an object")
    exact_keys(values, RUN_ONCE_OPTION_ORDER, "preflight run authority values")
    controller_path = str(Path(__file__).resolve(strict=True))
    argv = [controller_path, "--run-once"]
    for option in RUN_ONCE_OPTION_ORDER:
        value = values[option]
        if type(value) is not str or not value or len(value) > 4096 or "\0" in value:
            fail(f"preflight run authority value differs: {option}")
        try:
            value.encode("ascii")
        except UnicodeEncodeError:
            fail(f"preflight run authority value is not ASCII: {option}")
        argv.extend((option, value))
    if len(argv) != 2 + 2 * len(RUN_ONCE_OPTION_ORDER):
        fail("preflight run authority argv length differs")
    parse_recorded_run_argv(argv, controller_pid)
    return argv


def run_argv_sha256(argv: Sequence[str]) -> str:
    if (
        type(argv) is not list
        or len(argv) != 2 + 2 * len(RUN_ONCE_OPTION_ORDER)
        or argv[0] != str(Path(__file__).resolve(strict=True))
        or argv[1] != "--run-once"
        or any(
            type(item) is not str or not item or len(item) > 4096
            or "\0" in item or not item.isascii()
            for item in argv
        )
        or any(
            argv[2 + 2 * index] != option
            for index, option in enumerate(RUN_ONCE_OPTION_ORDER)
        )
    ):
        fail("run authority argv is not a string array")
    return sha256_bytes(
        b"".join(os.fsencode(item) + b"\0" for item in argv)
    )


def manifest_entry(manifest: Mapping[str, Any], relative: str) -> Mapping[str, Any]:
    matches = [
        entry for entry in manifest["entries"] if entry["path"] == relative
    ]
    if len(matches) != 1:
        fail(f"source manifest entry differs: {relative}")
    return matches[0]


def observe_preflight_python(deadline: Optional[float]) -> Dict[str, Any]:
    declared = Path(sys.executable)
    if not declared.is_absolute():
        fail("preflight Python executable path is not absolute")
    path = declared.resolve(strict=True)
    digest, info = hash_regular_path(
        path, "preflight Python executable", single_link=False,
        max_size=MAX_BINARY_BYTES, deadline=deadline,
    )
    if (
        not 1 <= info.st_size <= MAX_BINARY_BYTES
        or stat.S_IMODE(info.st_mode) != 0o755
        or info.st_uid != 0 or info.st_gid != 0
    ):
        fail("preflight Python executable policy differs")
    return {
        "declared_path": str(declared), "path": str(path),
        "sha256": digest, "stat": stat_receipt(info),
    }


def bootstrap_sampler_run_authority(
    args: argparse.Namespace, deadline: Optional[float],
) -> Dict[str, Any]:
    start_before = process_start_ticks(args.sampler_pid)
    security_before = read_process_security(
        args.sampler_pid, "preflight sampler bootstrap before"
    )
    validate_process_security(
        security_before, EXPECTED_CAMPAIGN_UID, EXPECTED_CAMPAIGN_GID,
        [EXPECTED_SAMPLER_I2C_GID], "preflight sampler bootstrap before",
    )
    affinity_before = sorted(os.sched_getaffinity(args.sampler_pid))
    if affinity_before != [EXPECTED_SAMPLER_CPU]:
        fail("preflight sampler bootstrap affinity differs")
    artifact_infos: Dict[str, os.stat_result] = {}
    for name, path in (
        ("csv", args.sampler_csv),
        ("pid_file", args.sampler_pid_file),
        ("validation", args.sampler_validation_jsonl),
        ("receipt", args.sampler_receipt),
    ):
        if path.resolve(strict=True) != path:
            fail(f"preflight sampler {name} path is not canonical")
        artifact_infos[name] = os.stat(str(path), follow_symlinks=False)
    cmdline = read_bounded_proc_vector(
        Path(f"/proc/{args.sampler_pid}/cmdline"),
        "preflight sampler bootstrap command line",
    )
    environ = read_bounded_proc_vector(
        Path(f"/proc/{args.sampler_pid}/environ"),
        "preflight sampler bootstrap environment",
    )
    parse_sealed_environment(environ, "preflight sampler bootstrap environment")
    executable = Path(f"/proc/{args.sampler_pid}/exe").resolve(strict=True)
    executable_sha256, executable_info = hash_regular_path(
        executable, "preflight sampler bootstrap executable",
        single_link=False, max_size=MAX_BINARY_BYTES, deadline=deadline,
    )
    if (
        not 1 <= executable_info.st_size <= MAX_BINARY_BYTES
        or stat.S_IMODE(executable_info.st_mode) != 0o755
        or executable_info.st_uid != 0 or executable_info.st_gid != 0
    ):
        fail("preflight sampler bootstrap executable policy differs")
    start_after = process_start_ticks(args.sampler_pid)
    security_after = read_process_security(
        args.sampler_pid, "preflight sampler bootstrap after"
    )
    affinity_after = sorted(os.sched_getaffinity(args.sampler_pid))
    if (
        start_before != start_after
        or not canonical_equal(security_before, security_after)
        or affinity_after != affinity_before
    ):
        fail("preflight sampler bootstrap identity changed")
    return {
        "cmdline_sha256": sha256_bytes(cmdline),
        "csv_device": artifact_infos["csv"].st_dev,
        "csv_inode": artifact_infos["csv"].st_ino,
        "environ_sha256": sha256_bytes(environ),
        "executable_sha256": executable_sha256,
        "pid_file_device": artifact_infos["pid_file"].st_dev,
        "pid_file_inode": artifact_infos["pid_file"].st_ino,
        "process_start_ticks": start_after,
        "receipt_device": artifact_infos["receipt"].st_dev,
        "receipt_inode": artifact_infos["receipt"].st_ino,
        "validation_device": artifact_infos["validation"].st_dev,
        "validation_inode": artifact_infos["validation"].st_ino,
    }


def preflight_run_values(
    preflight_args: argparse.Namespace, source: Mapping[str, Any],
    git: Mapping[str, Any], build: Mapping[str, Any],
    binary: Mapping[str, Any], python: Mapping[str, Any],
    sampler: Mapping[str, Any],
) -> Dict[str, str]:
    controller_entry = manifest_entry(
        source, "bench/Wh2DirectSystematicComplementScreen.py"
    )
    sampler_entry = manifest_entry(source, SAMPLER_SOURCE_PATH)
    if "pid_file" in sampler:
        sampler_values = {
            "cmdline_sha256": sampler["cmdline_sha256"],
            "csv_device": sampler["csv_device"],
            "csv_inode": sampler["csv_inode"],
            "environ_sha256": sampler["environ_sha256"],
            "executable_sha256": sampler["executable_sha256"],
            "pid_file_device": sampler["pid_file"]["stat"]["device"],
            "pid_file_inode": sampler["pid_file"]["stat"]["inode"],
            "process_start_ticks": sampler["process_start_ticks"],
            "receipt_device": sampler["receipt_file"]["stat"]["device"],
            "receipt_inode": sampler["receipt_file"]["stat"]["inode"],
            "validation_device": sampler["validation_jsonl"]["stat"]["device"],
            "validation_inode": sampler["validation_jsonl"]["stat"]["inode"],
        }
    else:
        sampler_values = dict(sampler)
    values = {
        "--binary": binary["path"],
        "--build-dir": build["root"],
        "--cpu": str(EXPECTED_TARGET_CPU),
        "--controller-cpu": str(EXPECTED_CONTROLLER_CPU),
        "--sampler-pid": str(preflight_args.sampler_pid),
        "--sampler-cpu": str(EXPECTED_SAMPLER_CPU),
        "--sampler-script": str(preflight_args.sampler_script),
        "--sampler-csv": str(preflight_args.sampler_csv),
        "--sampler-pid-file": str(preflight_args.sampler_pid_file),
        "--sampler-validation-jsonl": str(
            preflight_args.sampler_validation_jsonl
        ),
        "--sampler-receipt": str(preflight_args.sampler_receipt),
        "--expected-source-commit": preflight_args.expected_source_commit,
        "--expected-binary-sha256": binary["sha256"],
        "--expected-binary-uid": str(EXPECTED_CAMPAIGN_UID),
        "--expected-build-manifest-sha256": build["sha256"],
        "--expected-controller-sha256": controller_entry["sha256"],
        "--expected-controller-uid": str(EXPECTED_CAMPAIGN_UID),
        "--expected-controller-gid": str(EXPECTED_CAMPAIGN_GID),
        "--expected-git-sha256": git["executable"]["sha256"],
        "--expected-python-sha256": python["sha256"],
        "--expected-sampler-process-start-ticks": str(
            sampler_values["process_start_ticks"]
        ),
        "--expected-sampler-script-sha256": sampler_entry["sha256"],
        "--expected-sampler-csv-device": str(sampler_values["csv_device"]),
        "--expected-sampler-csv-inode": str(sampler_values["csv_inode"]),
        "--expected-sampler-pid-file-device": str(
            sampler_values["pid_file_device"]
        ),
        "--expected-sampler-pid-file-inode": str(
            sampler_values["pid_file_inode"]
        ),
        "--expected-sampler-validation-device": str(
            sampler_values["validation_device"]
        ),
        "--expected-sampler-validation-inode": str(
            sampler_values["validation_inode"]
        ),
        "--expected-sampler-receipt-device": str(
            sampler_values["receipt_device"]
        ),
        "--expected-sampler-receipt-inode": str(
            sampler_values["receipt_inode"]
        ),
        "--expected-sampler-cmdline-sha256": sampler_values["cmdline_sha256"],
        "--expected-sampler-environ-sha256": sampler_values["environ_sha256"],
        "--expected-sampler-executable-sha256": sampler_values[
            "executable_sha256"
        ],
        "--expected-sampler-uid": str(EXPECTED_CAMPAIGN_UID),
        "--expected-sampler-gid": str(EXPECTED_CAMPAIGN_GID),
        "--expected-sampler-i2c-gid": str(EXPECTED_SAMPLER_I2C_GID),
        "--expected-source-manifest-sha256": source["sha256"],
    }
    exact_keys(values, RUN_ONCE_OPTION_ORDER, "preflight run values")
    return values


def validate_preflight_sampler_run_binding(
    run: argparse.Namespace, sampler: Mapping[str, Any],
    source_entry: Mapping[str, Any], controller_image: Mapping[str, Any],
    where: str,
) -> None:
    """Cross-bind one live sampler receipt to the future run authority."""
    validate_sampler_receipt(sampler)
    if (
        sampler["pid"] != run.sampler_pid
        or sampler["cpu"] != run.sampler_cpu
        or sampler["process_start_ticks"]
        != run.expected_sampler_process_start_ticks
        or sampler["process_uid"] != run.expected_sampler_uid
        or sampler["process_gid"] != run.expected_sampler_gid
        or sampler["script_path"] != str(run.sampler_script)
        or sampler["script_sha256"]
        != run.expected_sampler_script_sha256
        or sampler["script_sha256"] != source_entry["sha256"]
        or not canonical_equal(sampler["script_stat"], source_entry["stat"])
        or sampler["csv_path"] != str(run.sampler_csv)
        or sampler["csv_device"] != run.expected_sampler_csv_device
        or sampler["csv_inode"] != run.expected_sampler_csv_inode
        or sampler["pid_file"]["path"] != str(run.sampler_pid_file)
        or sampler["pid_file"]["stat"]["device"]
        != run.expected_sampler_pid_file_device
        or sampler["pid_file"]["stat"]["inode"]
        != run.expected_sampler_pid_file_inode
        or sampler["validation_jsonl"]["path"]
        != str(run.sampler_validation_jsonl)
        or sampler["validation_jsonl"]["stat"]["device"]
        != run.expected_sampler_validation_device
        or sampler["validation_jsonl"]["stat"]["inode"]
        != run.expected_sampler_validation_inode
        or sampler["receipt_file"]["path"] != str(run.sampler_receipt)
        or sampler["receipt_file"]["stat"]["device"]
        != run.expected_sampler_receipt_device
        or sampler["receipt_file"]["stat"]["inode"]
        != run.expected_sampler_receipt_inode
        or sampler["cmdline_argv"] != expected_sampler_argv(run)
        or sampler["cmdline_sha256"]
        != sha256_bytes(expected_sampler_cmdline_bytes(run))
        or sampler["cmdline_sha256"]
        != run.expected_sampler_cmdline_sha256
        or sampler["environ_sha256"]
        != run.expected_sampler_environ_sha256
        or sampler["executable_path"] != controller_image["python_path"]
        or sampler["executable_sha256"]
        != run.expected_sampler_executable_sha256
        or sampler["executable_sha256"]
        != controller_image["python_sha256"]
        or not canonical_equal(
            sampler["executable_stat"], controller_image["python_stat"]
        )
    ):
        fail(where + " differs from the future run authority")
    validate_process_security(
        sampler["process_security"], run.expected_sampler_uid,
        run.expected_sampler_gid, [run.expected_sampler_i2c_gid],
        where + " process security",
    )


PREFLIGHT_SAMPLER_PREFIX_KEYS = {
    "final_pid_file_stat", "final_raw_stat", "final_receipt_stat",
    "final_validation_stat", "monitor_event", "pid_file_sha256",
    "raw_after_bytes", "raw_after_sha256", "raw_before_bytes",
    "raw_before_sha256", "receipt_bytes", "sampler_after_sha256",
    "sampler_before_sha256", "validation_after_bytes",
    "validation_after_sha256", "validation_before_bytes",
    "validation_before_sha256",
}


def preflight_sampler_prefix_binding(
    args: argparse.Namespace, before: Mapping[str, Any],
    after: Mapping[str, Any], handles: Mapping[str, Any],
) -> Dict[str, Any]:
    validate_sampler_growth_binding(before, after, "preflight sampler A/B")
    if not handles or handles.get("closed"):
        fail("preflight sampler prefix proof lacks held descriptors")
    fds = handles["file_fds"]

    def exact_prefix(
        fd: int, size: int, digest: str, where: str,
    ) -> bytes:
        data = os.pread(fd, size, 0)
        if len(data) != size or sha256_bytes(data) != digest:
            fail(where + " differs from its held prefix")
        return data

    raw_before = exact_prefix(
        fds["raw_csv"], before["csv_bytes"], before["csv_sha256"],
        "preflight sampler raw A",
    )
    raw_after = exact_prefix(
        fds["raw_csv"], after["csv_bytes"], after["csv_sha256"],
        "preflight sampler raw B",
    )
    validation_before = exact_prefix(
        fds["validation_jsonl"], before["validation_jsonl"]["bytes"],
        before["validation_jsonl"]["sha256"],
        "preflight sampler validation A",
    )
    validation_after = exact_prefix(
        fds["validation_jsonl"], after["validation_jsonl"]["bytes"],
        after["validation_jsonl"]["sha256"],
        "preflight sampler validation B",
    )
    if (
        not raw_before.endswith(b"\n") or not raw_after.endswith(b"\n")
        or not validation_before.endswith(b"\n")
        or not validation_after.endswith(b"\n")
        or not raw_after.startswith(raw_before)
        or not validation_after.startswith(validation_before)
    ):
        fail("preflight sampler A/B stream prefix law differs")
    final_infos = {
        name: os.fstat(fds[name])
        for name in ("raw_csv", "pid_file", "validation_jsonl", "receipt_file")
    }
    expected = {
        "raw_csv": after["csv_stat"],
        "pid_file": after["pid_file"]["stat"],
        "validation_jsonl": after["validation_jsonl"]["stat"],
        "receipt_file": after["receipt_file"]["stat"],
    }
    size_limits = {
        "raw_csv": MAX_SAMPLER_CSV_BYTES,
        "pid_file": 64,
        "validation_jsonl": MAX_SAMPLER_VALIDATION_BYTES,
        "receipt_file": MAX_SAMPLER_RECEIPT_BYTES,
    }
    for name, info in final_infos.items():
        receipt = expected[name]
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_dev != receipt["device"]
            or info.st_ino != receipt["inode"]
            or info.st_uid != receipt["uid"]
            or info.st_gid != receipt["gid"]
            or stat.S_IMODE(info.st_mode) != receipt["mode"]
            or info.st_nlink != 1
            or info.st_size < receipt["size"]
            or info.st_mtime_ns < receipt["mtime_ns"]
            or info.st_size > size_limits[name]
        ):
            fail("preflight sampler held identity changed: " + name)
    pid_info = final_infos["pid_file"]
    pid_data = os.pread(fds["pid_file"], pid_info.st_size, 0)
    receipt_info = final_infos["receipt_file"]
    receipt_data = os.pread(fds["receipt_file"], receipt_info.st_size, 0)
    if (
        pid_data != (str(args.sampler_pid) + "\n").encode("ascii")
        or receipt_data
    ):
        fail("preflight sampler PID/terminal reservation changed")
    event = poll_sampler_supervision(handles)
    if event != "none":
        fail("preflight sampler supervision event: " + event)
    value = {
        "final_pid_file_stat": stat_receipt(pid_info),
        "final_raw_stat": stat_receipt(final_infos["raw_csv"]),
        "final_receipt_stat": stat_receipt(receipt_info),
        "final_validation_stat": stat_receipt(final_infos["validation_jsonl"]),
        "monitor_event": event,
        "pid_file_sha256": sha256_bytes(pid_data),
        "raw_after_bytes": len(raw_after),
        "raw_after_sha256": sha256_bytes(raw_after),
        "raw_before_bytes": len(raw_before),
        "raw_before_sha256": sha256_bytes(raw_before),
        "receipt_bytes": len(receipt_data),
        "sampler_after_sha256": sha256_bytes(canonical_bytes(after)),
        "sampler_before_sha256": sha256_bytes(canonical_bytes(before)),
        "validation_after_bytes": len(validation_after),
        "validation_after_sha256": sha256_bytes(validation_after),
        "validation_before_bytes": len(validation_before),
        "validation_before_sha256": sha256_bytes(validation_before),
    }
    validate_preflight_sampler_prefix_binding(value, before, after)
    return value


def validate_preflight_sampler_prefix_binding(
    value: Any, before: Mapping[str, Any], after: Mapping[str, Any],
) -> None:
    if type(value) is not dict:
        fail("preflight sampler prefix binding is not an object")
    exact_keys(
        value, PREFLIGHT_SAMPLER_PREFIX_KEYS,
        "preflight sampler prefix binding",
    )
    validate_sampler_growth_binding(before, after, "preflight sampler A/B")
    for name in (
        "raw_before_sha256", "raw_after_sha256",
        "validation_before_sha256", "validation_after_sha256",
        "pid_file_sha256", "sampler_before_sha256", "sampler_after_sha256",
    ):
        lower_hash(value[name], "preflight sampler prefix " + name)
    for name in (
        "raw_before_bytes", "raw_after_bytes", "validation_before_bytes",
        "validation_after_bytes", "receipt_bytes",
    ):
        exact_int(
            value[name], 0, MAX_SAMPLER_CSV_BYTES + MAX_SAMPLER_VALIDATION_BYTES,
            "preflight sampler prefix " + name,
        )
    if (
        value["monitor_event"] != "none"
        or value["receipt_bytes"] != 0
        or value["raw_before_bytes"] != before["csv_bytes"]
        or value["raw_before_sha256"] != before["csv_sha256"]
        or value["raw_after_bytes"] != after["csv_bytes"]
        or value["raw_after_sha256"] != after["csv_sha256"]
        or value["validation_before_bytes"]
        != before["validation_jsonl"]["bytes"]
        or value["validation_before_sha256"]
        != before["validation_jsonl"]["sha256"]
        or value["validation_after_bytes"]
        != after["validation_jsonl"]["bytes"]
        or value["validation_after_sha256"]
        != after["validation_jsonl"]["sha256"]
        or value["sampler_before_sha256"]
        != sha256_bytes(canonical_bytes(before))
        or value["sampler_after_sha256"]
        != sha256_bytes(canonical_bytes(after))
    ):
        fail("preflight sampler prefix digest binding differs")
    for field, expected in (
        ("final_raw_stat", after["csv_stat"]),
        ("final_pid_file_stat", after["pid_file"]["stat"]),
        ("final_validation_stat", after["validation_jsonl"]["stat"]),
        ("final_receipt_stat", after["receipt_file"]["stat"]),
    ):
        validate_stat_receipt(value[field], "preflight sampler prefix " + field)
        if (
            value[field]["device"] != expected["device"]
            or value[field]["inode"] != expected["inode"]
            or value[field]["uid"] != expected["uid"]
            or value[field]["gid"] != expected["gid"]
            or value[field]["mode"] != expected["mode"]
            or value[field]["nlink"] != 1
            or value[field]["size"] < expected["size"]
            or value[field]["mtime_ns"] < expected["mtime_ns"]
        ):
            fail("preflight sampler prefix final stat differs")
    if (
        value["final_raw_stat"]["size"] > MAX_SAMPLER_CSV_BYTES
        or value["final_pid_file_stat"]["size"] > 64
        or value["final_validation_stat"]["size"]
        > MAX_SAMPLER_VALIDATION_BYTES
        or value["final_receipt_stat"]["size"]
        > MAX_SAMPLER_RECEIPT_BYTES
    ):
        fail("preflight sampler prefix final stat exceeds its bound")
    if (
        value["final_receipt_stat"]["size"] != 0
        or value["final_pid_file_stat"]["size"]
        != after["pid_file"]["bytes"]
        or value["final_raw_stat"]["size"] < value["raw_after_bytes"]
        or value["final_validation_stat"]["size"]
        < value["validation_after_bytes"]
    ):
        fail("preflight sampler prefix final size binding differs")


PREFLIGHT_SEAL_KEYS = {
    "binary_after", "binary_before", "build_manifest_after",
    "build_manifest_before", "expected_source_commit", "git_after",
    "git_before", "health_module_loader", "preflight_controller_after",
    "preflight_controller_before", "receipt_sha256", "run_argv",
    "run_argv_sha256", "sampler_after", "sampler_before",
    "sampler_prefix_binding", "schema", "source_manifest_after",
    "source_manifest_before", "source_root",
}


def validate_preflight_seal_receipt(value: Any) -> None:
    if type(value) is not dict:
        fail("preflight seal receipt is not an object")
    exact_keys(value, PREFLIGHT_SEAL_KEYS, "preflight seal receipt")
    exact_string(value["schema"], PREFLIGHT_SEAL_SCHEMA, "preflight schema")
    source_root = Path(__file__).resolve(strict=True).parents[1]
    exact_string(value["source_root"], str(source_root), "preflight source root")
    commit = value["expected_source_commit"]
    if type(commit) is not str or LOWER40.fullmatch(commit) is None:
        fail("preflight expected source commit differs")
    source_before = value["source_manifest_before"]
    source_after = value["source_manifest_after"]
    validate_source_manifest_receipt(source_before, "preflight source before")
    validate_source_manifest_receipt(source_after, "preflight source after")
    if not canonical_equal(source_before, source_after):
        fail("preflight source manifest changed")
    git_before = value["git_before"]
    git_after = value["git_after"]
    validate_git_receipt(git_before, commit, "preflight Git before")
    validate_git_receipt(git_after, commit, "preflight Git after")
    if not canonical_equal(git_before, git_after):
        fail("preflight Git receipt changed")
    for source_entry, blob_entry in zip(
        source_before["entries"], git_before["source_blob_oids"]
    ):
        if (
            source_entry["path"] != blob_entry["path"]
            or source_entry["git_blob_oid"] != blob_entry["head_oid"]
            or source_entry["git_blob_oid"] != blob_entry["worktree_oid"]
        ):
            fail("preflight source bytes/HEAD blob binding differs")
    build_before = value["build_manifest_before"]
    build_after = value["build_manifest_after"]
    validate_build_manifest(build_before, "preflight build before")
    validate_build_manifest(build_after, "preflight build after")
    if not canonical_equal(build_before, build_after):
        fail("preflight build manifest changed")
    binary_before = value["binary_before"]
    binary_after = value["binary_after"]
    validate_preflight_binary_observation(binary_before, "preflight binary before")
    validate_preflight_binary_observation(binary_after, "preflight binary after")
    if not canonical_equal(binary_before, binary_after):
        fail("preflight binary changed")
    controller_before = value["preflight_controller_before"]
    controller_after = value["preflight_controller_after"]
    validate_preflight_controller_receipt(
        controller_before, "preflight controller before"
    )
    validate_preflight_controller_receipt(
        controller_after, "preflight controller after"
    )
    if not canonical_equal(controller_before, controller_after):
        fail("preflight controller process/image changed")
    image = controller_after["image"]
    preflight_args = parse_cli_tokens(
        controller_after["argv"][1:], live_process=False
    )
    if commit != preflight_args.expected_source_commit:
        fail("preflight receipt/input source commit binding differs")
    controller_entry = manifest_entry(
        source_after, "bench/Wh2DirectSystematicComplementScreen.py"
    )
    if (
        controller_entry["sha256"] != image["script_sha256"]
        or controller_entry["bytes"] != image["script_stat"]["size"]
        or not canonical_equal(controller_entry["stat"], image["script_stat"])
        or git_after["executable"]["sha256"] != image["git_sha256"]
        or not canonical_equal(
            git_after["executable"]["stat"], image["git_stat"]
        )
    ):
        fail("preflight controller source/Git image binding differs")
    loader = value["health_module_loader"]
    validate_health_loader_receipt(loader)
    source_by_path = {
        entry["path"]: entry for entry in source_after["entries"]
    }
    for module in loader["modules"]:
        source_entry = source_by_path[module["path"]]
        if (
            module["bytes"] != source_entry["bytes"]
            or module["sha256"] != source_entry["sha256"]
        ):
            fail("preflight health loader/source binding differs")
    sampler_before = value["sampler_before"]
    sampler_after = value["sampler_after"]
    validate_sampler_receipt(sampler_before)
    validate_sampler_receipt(sampler_after)
    validate_preflight_sampler_prefix_binding(
        value["sampler_prefix_binding"], sampler_before, sampler_after
    )
    run_argv = value["run_argv"]
    run_args = parse_recorded_run_argv(run_argv, controller_after["pid"])
    digest = run_argv_sha256(run_argv)
    if digest != lower_hash(value["run_argv_sha256"], "preflight run argv"):
        fail("preflight run argv SHA-256 differs")
    sampler_source_entry = validate_sampler_source_authority(
        run_args, source_after, source_root
    )
    validate_preflight_sampler_run_binding(
        run_args, sampler_before, sampler_source_entry, image,
        "preflight sampler before",
    )
    validate_preflight_sampler_run_binding(
        run_args, sampler_after, sampler_source_entry, image,
        "preflight sampler after",
    )
    for name in (
        "binary", "build_dir", "sampler_pid", "sampler_script", "sampler_csv",
        "sampler_pid_file", "sampler_validation_jsonl", "sampler_receipt",
        "expected_source_commit",
    ):
        if getattr(preflight_args, name) != getattr(run_args, name):
            fail("preflight input/future run binding differs: " + name)
    expected_values = preflight_run_values(
        preflight_args, source_after, git_after, build_after, binary_after,
        {"sha256": image["python_sha256"]}, sampler_after,
    )
    expected_argv = build_run_authority_argv(
        expected_values, controller_after["pid"]
    )
    if not canonical_equal(run_argv, expected_argv):
        fail("preflight run argv differs from final observations")
    prefix = value["sampler_prefix_binding"]
    if prefix["pid_file_sha256"] != sampler_after["pid_file"]["sha256"]:
        fail("preflight sampler PID prefix binding differs")
    receipt_digest = lower_hash(
        value["receipt_sha256"], "preflight receipt SHA-256"
    )
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != receipt_digest:
        fail("preflight receipt self-hash differs")


def parse_preflight_seal_bytes(data: bytes) -> Dict[str, Any]:
    if (
        not data or len(data) > MAX_PREFLIGHT_RECEIPT_BYTES
        or not data.endswith(b"\n") or data.count(b"\n") != 1
        or b"\r" in data
    ):
        fail("preflight seal receipt framing differs")
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
        canonical = canonical_bytes(value) + b"\n"
    except (
        UnicodeDecodeError, json.JSONDecodeError, RecursionError, ValueError,
    ) as exc:
        fail("preflight seal receipt is malformed: " + str(exc))
    if canonical != data:
        fail("preflight seal receipt is not canonical JSON")
    validate_preflight_seal_receipt(value)
    return value


def _collect_preflight_seal(
    args: argparse.Namespace, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    if deadline is None:
        deadline = time.monotonic() + PREFLIGHT_SEAL_DEADLINE_SECONDS
    if not math.isfinite(deadline) or time.monotonic() >= deadline:
        fail("preflight seal reached its global deadline")
    launch_affinity = sorted(os.sched_getaffinity(0))
    if launch_affinity != [
        EXPECTED_TARGET_CPU, EXPECTED_CONTROLLER_CPU, EXPECTED_SAMPLER_CPU
    ]:
        fail("preflight launch affinity is not exact CPUs 120/121/122")
    args.controller_initial_affinity = launch_affinity
    os.sched_setaffinity(0, {EXPECTED_CONTROLLER_CPU})
    if os.sched_getaffinity(0) != {EXPECTED_CONTROLLER_CPU}:
        fail("preflight controller did not pin to singleton CPU121")
    require_fixed_namespace_absent()
    source_root = Path(__file__).resolve(strict=True).parents[1]
    expected_sampler_script = source_root / SAMPLER_SOURCE_PATH
    if args.sampler_script != expected_sampler_script:
        fail("preflight sampler script is not its source-manifest path")
    source_before = source_manifest(source_root, deadline)
    git_sha256, git_info = hash_regular_path(
        GIT_EXECUTABLE, "preflight Git bootstrap",
        max_size=MAX_BINARY_BYTES, deadline=deadline,
    )
    if (
        not 1 <= git_info.st_size <= MAX_BINARY_BYTES
        or stat.S_IMODE(git_info.st_mode) != 0o755
        or git_info.st_uid != 0 or git_info.st_gid != 0
        or git_info.st_nlink != 1
    ):
        fail("preflight Git bootstrap policy differs")
    git_before = git_receipt(
        source_root, args.expected_source_commit, git_sha256, deadline
    )
    build_before = build_manifest(args.build_dir, deadline)
    binary_before = observe_preflight_binary(Path(args.binary), deadline)
    python_before = observe_preflight_python(deadline)
    sampler_bootstrap = bootstrap_sampler_run_authority(args, deadline)
    initial_values = preflight_run_values(
        args, source_before, git_before, build_before, binary_before,
        python_before, sampler_bootstrap,
    )
    initial_argv = build_run_authority_argv(initial_values, os.getpid())
    run_args = parse_recorded_run_argv(initial_argv, os.getpid())
    run_args.controller_initial_affinity = launch_affinity
    validate_sampler_source_authority(run_args, source_before, source_root)
    health_modules = load_sealed_health_modules(source_root, source_before)
    controller_before = preflight_controller_receipt(run_args, deadline)
    observed_binary = preflight_binary(run_args, deadline)
    if (
        observed_binary["sha256"] != binary_before["sha256"]
        or not canonical_equal(observed_binary["stat"], binary_before["stat"])
        or controller_before["image"]["python_sha256"]
        != python_before["sha256"]
    ):
        fail("preflight bootstrap/static authority changed")
    sampler_before = make_sampler_attestation(
        run_args, 0, 0, health_modules, deadline
    )
    handles: Optional[Dict[str, Any]] = None
    try:
        handles = open_sampler_admission_handles(run_args, sampler_before)
        sampler_after = make_sampler_attestation(
            run_args, 0, 0, health_modules, deadline
        )
        binary_after = observe_preflight_binary(Path(args.binary), deadline)
        controller_after = preflight_controller_receipt(run_args, deadline)
        build_after = build_manifest(args.build_dir, deadline)
        git_after = git_receipt(
            source_root, args.expected_source_commit, git_sha256, deadline
        )
        source_after = source_manifest(source_root, deadline)
        if (
            not canonical_equal(source_before, source_after)
            or not canonical_equal(git_before, git_after)
            or not canonical_equal(build_before, build_after)
            or not canonical_equal(binary_before, binary_after)
            or not canonical_equal(controller_before, controller_after)
        ):
            fail("preflight static authority changed between passes")
        sampler_prefix = preflight_sampler_prefix_binding(
            run_args, sampler_before, sampler_after, handles
        )
        final_values = preflight_run_values(
            args, source_after, git_after, build_after, binary_after,
            {"sha256": controller_after["image"]["python_sha256"]},
            sampler_after,
        )
        run_argv = build_run_authority_argv(
            final_values, controller_after["pid"]
        )
        if not canonical_equal(run_argv, initial_argv):
            fail("preflight run authority changed between passes")
        if poll_sampler_supervision(handles) != "none":
            fail("preflight sampler changed after its prefix proof")
        require_fixed_namespace_absent()
        value: Dict[str, Any] = {
            "binary_after": binary_after,
            "binary_before": binary_before,
            "build_manifest_after": build_after,
            "build_manifest_before": build_before,
            "expected_source_commit": args.expected_source_commit,
            "git_after": git_after,
            "git_before": git_before,
            "health_module_loader": health_modules.receipt,
            "preflight_controller_after": controller_after,
            "preflight_controller_before": controller_before,
            "receipt_sha256": None,
            "run_argv": run_argv,
            "run_argv_sha256": run_argv_sha256(run_argv),
            "sampler_after": sampler_after,
            "sampler_before": sampler_before,
            "sampler_prefix_binding": sampler_prefix,
            "schema": PREFLIGHT_SEAL_SCHEMA,
            "source_manifest_after": source_after,
            "source_manifest_before": source_before,
            "source_root": str(source_root),
        }
        value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
        validate_preflight_seal_receipt(value)
        if time.monotonic() >= deadline:
            fail("preflight seal reached its global deadline")
        if poll_sampler_supervision(handles) != "none":
            fail("preflight sampler changed before receipt emission")
        require_fixed_namespace_absent()
        return value
    finally:
        close_sampler_admission_handles(handles)


def collect_preflight_seal(
    args: argparse.Namespace, deadline: Optional[float] = None,
) -> Dict[str, Any]:
    global PREFLIGHT_SIDE_EFFECT_GUARD
    if PREFLIGHT_SIDE_EFFECT_GUARD:
        fail("preflight seal collection is already active")
    PREFLIGHT_SIDE_EFFECT_GUARD = True
    try:
        return _collect_preflight_seal(args, deadline)
    finally:
        PREFLIGHT_SIDE_EFFECT_GUARD = False


def write_preflight_payload(fd: int, payload: bytes, deadline: float) -> None:
    if (
        type(fd) is not int or fd < 0 or type(payload) is not bytes
        or not payload or len(payload) > MAX_PREFLIGHT_RECEIPT_BYTES
        or not math.isfinite(deadline)
    ):
        fail("preflight seal output parameters differ")
    try:
        was_blocking = os.get_blocking(fd)
    except OSError as exc:
        fail("preflight seal output descriptor failed: " + str(exc))
    try:
        try:
            os.set_blocking(fd, False)
        except (OSError, ValueError) as exc:
            fail("preflight seal output nonblocking setup failed: " + str(exc))
        view = memoryview(payload)
        while view:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                fail("preflight seal output reached its global deadline")
            try:
                _, writable, _ = select.select([], [fd], [], remaining)
            except InterruptedError:
                continue
            except (OSError, ValueError) as exc:
                fail("preflight seal output poll failed: " + str(exc))
            if not writable:
                fail("preflight seal output reached its global deadline")
            try:
                written = os.write(fd, view)
            except (BlockingIOError, InterruptedError):
                continue
            except OSError as exc:
                fail("preflight seal receipt write failed: " + str(exc))
            if written <= 0:
                fail("preflight seal receipt write was short")
            view = view[written:]
    finally:
        try:
            os.set_blocking(fd, was_blocking)
        except (OSError, ValueError) as exc:
            fail("preflight seal output blocking restore failed: " + str(exc))


def emit_preflight_seal(args: argparse.Namespace) -> int:
    deadline = time.monotonic() + PREFLIGHT_SEAL_DEADLINE_SECONDS
    receipt = collect_preflight_seal(args, deadline)
    payload = canonical_bytes(receipt) + b"\n"
    if len(payload) > MAX_PREFLIGHT_RECEIPT_BYTES:
        fail("preflight seal receipt exceeds its output bound")
    write_preflight_payload(1, payload, deadline)
    return 0


def run_once(args: argparse.Namespace) -> int:
    require_preflight_side_effects_allowed("run-once dispatch")
    controller_started_ns = time.monotonic_ns()
    controller_started = controller_started_ns / 1_000_000_000.0
    if not canonical_equal(dict(os.environ), SEALED_LAUNCH_ENVIRONMENT):
        fail("controller startup environment differs")
    source_root = Path(__file__).resolve(strict=True).parents[1]
    binary_path = Path(args.binary)
    guarded_signals = SignalGuard()
    with guarded_signals:
        require_authorized_process_identity(
            args.expected_controller_uid, args.expected_controller_gid
        )
        launch_affinity = sorted(os.sched_getaffinity(0))
        if (
            EXPECTED_TARGET_CPU not in launch_affinity
            or EXPECTED_CONTROLLER_CPU not in launch_affinity
            or EXPECTED_SAMPLER_CPU not in launch_affinity
        ):
            fail("controller launch affinity lacks frozen CPUs 120/121/122")
        args.controller_initial_affinity = launch_affinity
        os.sched_setaffinity(0, {EXPECTED_CONTROLLER_CPU})
        if os.sched_getaffinity(0) != {EXPECTED_CONTROLLER_CPU}:
            fail("controller did not pin to exact singleton CPU121")

        # Every fallible identity/topology/sampler admission check precedes the
        # durable one-shot claim.  No campaign namespace is consumed on a
        # failed preflight.
        require_fixed_namespace_absent()
        source_before = source_manifest(
            source_root, controller_started + 9.0
        )
        if source_before["sha256"] != args.expected_source_manifest_sha256:
            fail("source manifest differs from the presealed value")
        validate_sampler_source_authority(args, source_before, source_root)
        git_before = git_receipt(
            source_root, args.expected_source_commit,
            args.expected_git_sha256,
            controller_started + 9.0,
        )
        health_modules = load_sealed_health_modules(source_root, source_before)
        health_module_loader = health_modules.receipt
        controller_before = controller_provenance(
            args, controller_started + 9.0
        )
        controller_entry = next(
            entry for entry in source_before["entries"]
            if entry["path"] == "bench/Wh2DirectSystematicComplementScreen.py"
        )
        if (
            controller_entry["sha256"] != args.expected_controller_sha256
            or controller_entry["bytes"] != controller_before["script_stat"]["size"]
            or not canonical_equal(
                controller_entry["stat"], controller_before["script_stat"]
            )
        ):
            fail("controller/source-manifest binding differs")
        build_before = build_manifest(
            args.build_dir, controller_started + 9.0
        )
        if build_before["sha256"] != args.expected_build_manifest_sha256:
            fail("build manifest differs from the presealed value")
        preflight_binary(args, controller_started + 9.0)
        sampler_admission = make_sampler_attestation(
            args, 0, 0, health_modules, controller_started + 9.0
        )
        health_state = prepare_health(
            args, EXPECTED_TARGET_CPU,
            controller_started + CONTROLLER_DEADLINE_SECONDS,
            health_modules, sampler_admission,
        )
        if not health_state.get("ready"):
            failures = health_state.get("collection_failures", [])
            failures += health_state.get("violations", [])
            fail("health admission failed before claim: " + bounded_failure_text(failures))
        if guarded_signals.first_signal is not None:
            return 1
        if monotonic_deadline_reached(
            time.monotonic_ns(), controller_started_ns, 10.0,
        ):
            fail("preflight consumed the bounded ten-second admission allowance")

        execution = ExecutionResult(
            binary=empty_binary_receipt(
                binary_path, args.expected_binary_sha256
            ),
            health=dict(health_state["receipt"]),
        )
        source_after: Dict[str, Any] = {}
        git_after: Dict[str, Any] = {}
        config: Dict[str, Any] = {}
        statistics: Dict[str, Any] = {}
        gate_failures: List[str] = []
        errors: List[str] = []
        stage_failures: List[str] = []
        worker_returned = False
        claim = ClaimReservation.reserve(
            args, controller_before, controller_started_ns
        )
        try:
            bundle = OutputBundle.reserve(
                FIXED_OUTPUT_DIR, EXPECTED_TARGET_CPU,
                args.expected_source_commit,
                args.expected_source_manifest_sha256,
                args.expected_git_sha256, source_root,
                guarded_signals,
            )
        except BaseException:
            # The claim is intentionally durable and never removed.  Its
            # existence records consumption even if output reservation fails.
            raise

        try:
            # This is deliberately the first operation after ownership succeeds:
            # SignalGuard's emergency path can never fall back to health={}.  The
            # provisional receipt is already canonical and self-hashed.
            bundle.update_fallback_summary_components({
                "binary": execution.binary,
                "config": config,
                "elapsed_seconds": 0.0,
                "failure": "post-claim controller initialization",
                "git_after": git_after,
                "git_before": git_before,
                "health": execution.health,
                "health_module_loader": health_module_loader,
                "outcome": "invalid",
                "process_exit_code": execution.returncode,
                "source_manifest_after": source_after,
                "source_manifest_before": source_before,
                "statistics": statistics,
                "target_cpu": EXPECTED_TARGET_CPU,
            })
            try:
                execution = run_binary(
                    binary_path, EXPECTED_TARGET_CPU,
                    args.expected_binary_sha256,
                    controller_started=controller_started,
                    controller_started_ns=controller_started_ns,
                    health_args=args, health_modules=health_modules,
                    prepared_health_state=health_state,
                    signal_guard=guarded_signals,
                )
                worker_returned = True
                bundle.update_fallback_summary_components({
                    "binary": execution.binary,
                    "elapsed_seconds": execution.elapsed_seconds,
                    "failure": "post-claim controller execution",
                    "health": execution.health,
                    "process_exit_code": execution.returncode,
                })
                for message in execution.errors:
                    append_unique_failure(errors, message)
            except BaseException as exc:
                append_unique_failure(errors, "runner body: " + exception_text(exc))
            finally:
                if not worker_returned:
                    try:
                        execution.health, health_errors = finish_health(
                            args, EXPECTED_TARGET_CPU, None, None, health_state,
                            controller_started + HEALTH_COLLECTION_DEADLINE_SECONDS,
                            health_modules,
                        )
                        for message in health_errors:
                            append_unique_failure(errors, message)
                    except BaseException as exc:
                        partial = dict(health_state["receipt"])
                        collection = list(partial.get("collection_failures", []))
                        append_unique_failure(
                            collection, "health terminal collection: " + exception_text(exc)
                        )
                        execution.health = finalize_health_receipt(
                            partial, collection, partial.get("violations", [])
                        )

            bundle.update_fallback_summary_components({
                "binary": execution.binary,
                "elapsed_seconds": execution.elapsed_seconds,
                "failure": "post-claim controller execution",
                "health": execution.health,
                "process_exit_code": execution.returncode,
            })

            try:
                validate_binary_receipt(execution.binary, require_complete=False)
            except BaseException as exc:
                append_unique_failure(errors, "binary receipt: " + exception_text(exc))
            try:
                validate_health_receipt(
                    execution.health, execution.binary, require_complete=False
                )
            except BaseException as exc:
                append_unique_failure(errors, "health receipt: " + exception_text(exc))
            if execution.binary.get("timed_out"):
                append_unique_failure(
                    errors, "screen reached the hard 110-second child deadline"
                )
            if execution.binary.get("output_overflow"):
                append_unique_failure(errors, "screen exceeded its output allowance")
            if execution.binary.get("process_started") and execution.returncode != 0:
                append_unique_failure(errors, f"screen exited {execution.returncode}")
            if execution.stderr:
                append_unique_failure(errors, "screen emitted stderr")
            if guarded_signals.first_signal is not None:
                append_unique_failure(
                    errors, "runner received " + guarded_signals.first_signal
                )
            try:
                source_after = source_manifest(
                    source_root, controller_started + 116.0
                )
                git_after = git_receipt(
                    source_root, args.expected_source_commit,
                    args.expected_git_sha256,
                    controller_started + 116.0,
                )
                if (
                    not canonical_equal(source_after, source_before)
                    or not canonical_equal(git_after, git_before)
                ):
                    fail("source/Git receipt changed during worker execution")
            except BaseException as exc:
                append_unique_failure(errors, "source postcheck: " + exception_text(exc))
            bundle.update_fallback_summary_components({
                "git_after": git_after,
                "source_manifest_after": source_after,
            })

            if execution.raw_complete():
                try:
                    records = parse_transcript(execution.raw)
                    config, statistics, gate_failures = validate_transcript(
                        records, EXPECTED_TARGET_CPU, args.expected_source_commit
                    )
                except BaseException as exc:
                    append_unique_failure(
                        errors, "transcript validation: " + exception_text(exc)
                    )
            try:
                stage_failures = bundle.stage_evidence(
                    execution.raw, execution.stderr
                )
            except BaseException as exc:
                stage_failures = ["evidence staging: " + exception_text(exc)]
            for message in stage_failures:
                append_unique_failure(errors, message)
            config = bundle.raw_validation_config
            statistics = bundle.raw_validation_statistics
            gate_failures = list(bundle.raw_validation_failures)

            outcome = "invalid"
            if not errors and config and statistics:
                outcome = "reject" if gate_failures else "pass"
            elif not errors:
                append_unique_failure(errors, "runner did not reach validation")
            try:
                validate_health_receipt(
                    execution.health, execution.binary,
                    require_complete=outcome in ("pass", "reject"),
                )
            except BaseException as exc:
                append_unique_failure(
                    errors, "health adjudication: " + exception_text(exc)
                )
                outcome = "invalid"
            if errors:
                outcome = "invalid"
            snapshot_signal = guarded_signals.first_signal
            snapshot_signals = list(guarded_signals.observed_signals)
            bundle.signal_name = snapshot_signal
            bundle.signal_names = snapshot_signals
            failure = (
                bounded_failure_text(errors) if errors
                else ("statistical-gates" if gate_failures else "none")
            )
            raw_complete = bool(config and statistics)
            components = {
                "binary": execution.binary,
                "config": config,
                "elapsed_seconds": execution.elapsed_seconds,
                "failure": failure,
                "git_after": git_after,
                "git_before": git_before,
                "health": execution.health,
                "health_module_loader": health_module_loader,
                "outcome": outcome,
                "process_exit_code": execution.returncode,
                "source_manifest_after": source_after,
                "source_manifest_before": source_before,
                "statistics": statistics,
                "target_cpu": EXPECTED_TARGET_CPU,
            }
            bundle.bind_expected_summary_components(components)
            summary = make_summary_preimage({
                **components,
                "expected_source_git_commit": args.expected_source_commit,
                "expected_source_manifest_sha256": (
                    args.expected_source_manifest_sha256
                ),
                "output_bundle": bundle.receipt(),
                "publication_failures": list(stage_failures),
                "raw_bytes": len(bundle.staged_raw),
                "raw_complete": raw_complete,
                "raw_record_count": bundle.staged_raw.count(b"\n"),
                "raw_sha256": sha256_bytes(bundle.staged_raw),
                "schema": SUMMARY_SCHEMA,
                "signal": snapshot_signal,
                "signal_names": snapshot_signals,
                "stderr_bytes": len(bundle.staged_stderr),
                "stderr_sha256": sha256_bytes(bundle.staged_stderr),
                "summary_preimage_sha256": None,
            })
            summary_bytes = canonical_bytes(summary) + b"\n"
            committed, final_outcome, _ = final_publish(
                args, bundle, claim, guarded_signals, controller_started,
                controller_started_ns, outcome, errors, execution, source_before,
                git_before, build_before, controller_before, source_root,
                health_modules, summary_bytes,
            )
            if not committed or final_outcome == "invalid":
                return 1
            return 0 if final_outcome == "pass" else 2
        except BaseException as exc:
            append_unique_failure(
                errors, "post-claim owner: " + exception_text(exc)
            )
            if bundle.closed:
                return 1
            # Refresh the most complete health/binary evidence before any
            # fallback summary is constructed.  This owner is entered for
            # every representable exception after output reservation.
            try:
                bundle.update_fallback_summary_components({
                    "binary": execution.binary,
                    "config": config,
                    "elapsed_seconds": execution.elapsed_seconds,
                    "failure": bounded_failure_text(errors),
                    "git_after": git_after,
                    "git_before": git_before,
                    "health": execution.health,
                    "health_module_loader": health_module_loader,
                    "outcome": "invalid",
                    "process_exit_code": execution.returncode,
                    "source_manifest_after": source_after,
                    "source_manifest_before": source_before,
                    "statistics": statistics,
                    "target_cpu": EXPECTED_TARGET_CPU,
                })
            except BaseException as fallback_exc:
                append_unique_failure(
                    errors,
                    "fallback component refresh: "
                    + exception_text(fallback_exc),
                )
            try:
                recovery_stage = bundle.stage_evidence(
                    execution.raw, execution.stderr
                )
                for message in recovery_stage:
                    append_unique_failure(errors, message)
                bundle.update_fallback_summary_components({
                    "config": bundle.raw_validation_config,
                    "statistics": bundle.raw_validation_statistics,
                })
            except BaseException as stage_exc:
                append_unique_failure(
                    errors, "fallback evidence staging: " + exception_text(stage_exc)
                )
            bundle.signal_name = guarded_signals.first_signal
            bundle.signal_names = list(guarded_signals.observed_signals)
            try:
                emergency_bytes = bundle.emergency_summary(
                    "post-claim v2 owner", errors
                )
                final_publish(
                    args, bundle, claim, guarded_signals, controller_started,
                    controller_started_ns, "invalid", errors, execution,
                    source_before, git_before, build_before, controller_before,
                    source_root, health_modules, emergency_bytes,
                )
                return 1
            except BaseException as final_exc:
                append_unique_failure(
                    errors, "post-claim final owner: " + exception_text(final_exc)
                )
                try:
                    poison_noncommitted_publication(
                        bundle, bundle.file_fds, errors, guarded_signals, None
                    )
                except BaseException:
                    pass
                for name, fd in bundle.file_fds.items():
                    try:
                        if name != "summary.json":
                            os.fchmod(fd, 0o400)
                            os.fsync(fd)
                    except OSError:
                        pass
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                bundle.file_fds.clear()
                try:
                    os.fchmod(bundle.directory_fd, 0o500)
                    os.fsync(bundle.directory_fd)
                    os.fsync(bundle.parent_fd)
                except OSError:
                    pass
                for fd in (bundle.directory_fd, bundle.parent_fd):
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                bundle.closed = True
                return 1


def retained_read_exact(
    fd: int, limit: int, where: str,
) -> Tuple[bytes, os.stat_result]:
    exact_int(limit, 1, MAX_INT63, where + " limit")
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
        or stat.S_IMODE(before.st_mode) != 0o400
        or before.st_size < 0 or before.st_size > limit
    ):
        fail(f"{where} retained-file policy differs")
    data = bytearray()
    offset = 0
    while offset < before.st_size:
        block = os.pread(fd, min(1024 * 1024, before.st_size - offset), offset)
        if not block:
            fail(f"{where} retained read is short")
        data.extend(block)
        offset += len(block)
    after = os.fstat(fd)
    if not same_file_receipt(before, after):
        fail(f"{where} changed during retained read")
    return bytes(data), after


def retained_artifact_receipt(
    info: os.stat_result, data: bytes,
) -> Dict[str, Any]:
    return {
        "bytes": len(data), "device": info.st_dev, "gid": info.st_gid,
        "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink, "sha256": sha256_bytes(data),
        "uid": info.st_uid,
    }


def retained_optional_stat(parent_fd: int, name: str) -> Optional[os.stat_result]:
    try:
        return os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        return None


def retained_roster(directory_fd: int) -> None:
    names: List[str] = []
    with os.scandir(directory_fd) as entries:
        for entry in entries:
            if len(names) >= len(FINAL_OUTPUT_NAMES):
                fail("retained output has too many directory entries")
            if not entry.is_file(follow_symlinks=False):
                fail("retained output contains a non-regular entry")
            names.append(entry.name)
    if len(names) != len(FINAL_OUTPUT_NAMES) or set(names) != set(FINAL_OUTPUT_NAMES):
        fail("retained output roster differs")


def validate_recorded_run_authority(
    expected_argv_sha256: str, controller: Mapping[str, Any],
    summary: Mapping[str, Any], claim_document: Mapping[str, Any],
) -> argparse.Namespace:
    before = controller["controller_before"]
    if before["argv_sha256"] != expected_argv_sha256:
        fail("retained run argv differs from the external authority seal")
    run = parse_recorded_run_argv(before["argv"], before["pid"])
    binary = controller["binary"]
    build = controller["build_before"]
    source = controller["source_manifest_before"]
    git = controller["git_before"]
    health_sampler = summary["health"].get("sampler", {})
    sampler_admission = summary["health"].get("sampler_admission", {})
    child_pid = binary.get("child_pid")
    sampler_identity_collision = (
        bool(health_sampler) and child_pid == run.sampler_pid
    )
    script_entry = source["entries"][SOURCE_PATHS.index(
        "bench/Wh2DirectSystematicComplementScreen.py"
    )]
    sampler_source_entry = validate_sampler_source_authority(
        run, source, Path(__file__).resolve(strict=True).parents[1]
    )
    if (
        run.binary != binary["path"]
        or run.expected_binary_sha256 != binary["expected_sha256"]
        or run.expected_binary_sha256 != claim_document["binary_sha256"]
        or str(run.build_dir) != build["root"]
        or run.expected_build_manifest_sha256 != build["sha256"]
        or run.expected_build_manifest_sha256
        != claim_document["build_manifest_sha256"]
        or run.expected_controller_sha256 != before["script_sha256"]
        or run.expected_controller_sha256 != script_entry["sha256"]
        or run.expected_controller_uid != claim_document["uid"]
        or run.expected_controller_gid != claim_document["gid"]
        or run.expected_controller_uid != before["script_stat"]["uid"]
        or run.expected_controller_gid != before["script_stat"]["gid"]
        or run.expected_controller_uid != before["process_uid"]
        or run.expected_controller_gid != before["process_gid"]
        or run.expected_git_sha256 != before["git_sha256"]
        or run.expected_git_sha256 != git["executable"]["sha256"]
        or run.expected_git_sha256 != claim_document["git_sha256"]
        or run.expected_python_sha256 != before["python_sha256"]
        or run.expected_source_commit != claim_document["source_commit"]
        or run.expected_source_commit != git["head"]
        or run.expected_source_commit != summary["expected_source_git_commit"]
        or run.expected_source_manifest_sha256 != source["sha256"]
        or run.expected_source_manifest_sha256
        != claim_document["source_manifest_sha256"]
        or run.expected_source_manifest_sha256
        != summary["expected_source_manifest_sha256"]
        or run.cpu != controller["target_cpu"]
        or run.controller_cpu != before["controller_cpu"]
        or child_pid == before["pid"]
        or sampler_identity_collision
    ):
        fail("recorded run argv authority binding differs")
    validate_process_security(
        before["process_security"], run.expected_controller_uid,
        run.expected_controller_gid, [],
        "recorded controller process security authority",
    )
    # The external owner seal authorizes the exact image checked before exec.
    # Invalid evidence may truthfully preserve a later chown/replacement in
    # either post-exec receipt.  Positive receipts already require both later
    # stats to equal stat_before in validate_binary_receipt(require_complete).
    stat_before = binary["stat_before"]
    if (
        binary["sha256_before"] is not None and stat_before
        and (
            stat_before["uid"] != run.expected_binary_uid
            or stat_before["gid"] != run.expected_controller_gid
        )
    ):
        fail("recorded pre-exec binary owner authority binding differs")
    samplers = [
        sampler_admission, health_sampler, controller["sampler_after"],
    ]
    require_frozen_sampler_stat = controller["outcome"] in ("pass", "reject")
    for sampler in samplers:
        if not sampler:
            continue
        if (
            sampler["pid"] != run.sampler_pid
            or sampler["cpu"] != run.sampler_cpu
            or sampler["script_path"] != str(run.sampler_script)
            or sampler["csv_path"] != str(run.sampler_csv)
            or sampler["process_start_ticks"]
            != run.expected_sampler_process_start_ticks
            or sampler["script_sha256"]
            != run.expected_sampler_script_sha256
            or sampler["csv_device"] != run.expected_sampler_csv_device
            or sampler["csv_inode"] != run.expected_sampler_csv_inode
            or sampler["cmdline_sha256"]
            != run.expected_sampler_cmdline_sha256
            or sampler["cmdline_argv"] != expected_sampler_argv(run)
            or sampler["cmdline_sha256"]
            != sha256_bytes(expected_sampler_cmdline_bytes(run))
            or sampler["environ_sha256"]
            != run.expected_sampler_environ_sha256
            or sampler["executable_sha256"]
            != run.expected_sampler_executable_sha256
            or sampler["process_uid"] != run.expected_sampler_uid
            or sampler["process_gid"] != run.expected_sampler_gid
            or sampler["pid_file"]["path"] != str(run.sampler_pid_file)
            or sampler["pid_file"]["stat"]["device"]
            != run.expected_sampler_pid_file_device
            or sampler["pid_file"]["stat"]["inode"]
            != run.expected_sampler_pid_file_inode
            or sampler["validation_jsonl"]["path"]
            != str(run.sampler_validation_jsonl)
            or sampler["validation_jsonl"]["stat"]["device"]
            != run.expected_sampler_validation_device
            or sampler["validation_jsonl"]["stat"]["inode"]
            != run.expected_sampler_validation_inode
            or sampler["receipt_file"]["path"] != str(run.sampler_receipt)
            or sampler["receipt_file"]["stat"]["device"]
            != run.expected_sampler_receipt_device
            or sampler["receipt_file"]["stat"]["inode"]
            != run.expected_sampler_receipt_inode
            or sampler["script_stat"]["size"]
            != sampler_source_entry["bytes"]
            or require_frozen_sampler_stat and not canonical_equal(
                sampler["script_stat"], sampler_source_entry["stat"]
            )
        ):
            fail("recorded sampler authority binding differs")
        validate_process_security(
            sampler["process_security"], run.expected_sampler_uid,
            run.expected_sampler_gid, [run.expected_sampler_i2c_gid],
            "recorded sampler process security authority",
        )
    return run


def validate_retained_claim_controller_binding(
    claim_document: Mapping[str, Any], claim_receipt: Mapping[str, Any],
    controller: Mapping[str, Any], summary: Mapping[str, Any],
) -> None:
    before = controller["controller_before"]
    if (
        claim_document["controller_receipt_sha256"]
        != before["receipt_sha256"]
        or claim_document["controller_started_monotonic_ns"]
        != controller["controller_started_monotonic_ns"]
        or claim_document["pid"] != before["pid"]
        or claim_document["process_start_ticks"] != before["process_start_ticks"]
        or claim_document["uid"] != claim_receipt["stat"]["uid"]
        or claim_document["gid"] != claim_receipt["stat"]["gid"]
        or claim_document["output_path"] != str(FIXED_OUTPUT_DIR)
        or claim_document["parent_device"]
        != claim_receipt["parent"]["device"]
        or claim_document["parent_inode"] != claim_receipt["parent"]["inode"]
        or not canonical_equal(
            controller["source_manifest_before"],
            summary["source_manifest_before"],
        )
        or not canonical_equal(controller["git_before"], summary["git_before"])
    ):
        fail("retained claim/controller/summary provenance binding differs")
    if controller["outcome"] in ("pass", "reject"):
        for name in ("source_manifest_after", "git_after"):
            if controller[name] and summary[name] and not canonical_equal(
                controller[name], summary[name]
            ):
                fail(f"retained controller/summary {name} differs")


def validate_live_retained_verifier(
    controller: Mapping[str, Any], source: Mapping[str, Any],
) -> None:
    before = controller["controller_before"]
    script_relative = "bench/Wh2DirectSystematicComplementScreen.py"
    script_entry = source["entries"][SOURCE_PATHS.index(script_relative)]
    source_root = Path(__file__).resolve(strict=True).parents[1]
    data, info = read_sealed_source(source_root, script_relative, script_entry)
    if (
        sha256_bytes(data) != before["script_sha256"]
        or not canonical_equal(stat_receipt(info), before["script_stat"])
    ):
        fail("live retained verifier script differs from the recorded controller")
    declared = Path(sys.executable)
    if str(declared) != before["python_declared_path"]:
        fail("live verifier declared Python path differs")
    resolved = declared.resolve(strict=True)
    digest, named_info = hash_regular_path(
        resolved, "retained verifier Python", single_link=False,
        max_size=MAX_BINARY_BYTES,
    )
    live_fd = os.open(
        "/proc/self/exe",
        nonblocking_read_flags(
            "retained verifier live Python", nofollow=False
        ),
    )
    try:
        live_before = os.fstat(live_fd)
        if (
            not stat.S_ISREG(live_before.st_mode)
            or not 1 <= live_before.st_size <= MAX_BINARY_BYTES
        ):
            fail("retained verifier live Python size/type differs")
        live_digest = file_sha256_fd(live_fd, live_before.st_size)
        live_after = os.fstat(live_fd)
    finally:
        os.close(live_fd)
    if (
        str(resolved) != before["python_path"]
        or digest != before["python_sha256"] or live_digest != digest
        or not canonical_equal(stat_receipt(named_info), before["python_stat"])
        or not same_file_receipt(named_info, live_before)
        or not same_file_receipt(live_before, live_after)
    ):
        fail("live retained verifier Python differs from the recorded runtime")


def validate_retained_science(
    raw: bytes, stderr: bytes, thermal: bytes, summary: Mapping[str, Any],
) -> None:
    if (
        summary["raw_bytes"] != len(raw)
        or summary["raw_record_count"] != raw.count(b"\n")
        or summary["raw_sha256"] != sha256_bytes(raw)
        or summary["stderr_bytes"] != len(stderr)
        or summary["stderr_sha256"] != sha256_bytes(stderr)
    ):
        fail("retained raw/stderr summary binding differs")
    health = summary["health"]
    if not health:
        fail("retained claimed outcome lacks a health receipt")
    expected_thermal = health_thermal_artifact_bytes(health)
    if thermal != expected_thermal:
        fail("retained thermal artifact differs from the health receipt")
    raw_config: Dict[str, Any] = {}
    raw_statistics: Dict[str, Any] = {}
    raw_gates: List[str] = []
    raw_complete = False
    try:
        records = parse_transcript(raw)
        raw_config, raw_statistics, raw_gates = validate_transcript(
            records, summary["target_cpu"], summary["expected_source_git_commit"]
        )
        raw_complete = True
    except Exception:
        # Producer staging treats every ordinary parser/validator exception as
        # bounded invalid scientific evidence.  Mirror that total classifier
        # here so ValueError (for example Python's integer-digit guard) and
        # RecursionError cannot turn an authentic invalid into malformed.
        raw_complete = False
    if summary["raw_complete"] != raw_complete:
        fail("retained raw completeness classification differs")
    if raw_complete:
        if not canonical_equal(raw_config, summary["config"]) or not canonical_equal(
            raw_statistics, summary["statistics"]
        ):
            fail("retained scientific recomputation differs")
        if summary["outcome"] in ("pass", "reject") and (
            summary["outcome"] != ("reject" if raw_gates else "pass")
        ):
            fail("retained scientific decision differs")
    elif summary["config"] or summary["statistics"]:
        fail("incomplete retained raw data carries scientific receipts")


def validate_retained_output_binding(
    parent_info: os.stat_result, directory_info: os.stat_result,
    file_infos: Mapping[str, os.stat_result], data: Mapping[str, bytes],
    claim_document: Mapping[str, Any], claim_receipt: Mapping[str, Any],
    summary: Mapping[str, Any], controller: Mapping[str, Any],
) -> None:
    output = summary["output_bundle"]
    if output["path"] != str(FIXED_OUTPUT_DIR):
        fail("retained summary output path differs")
    directory = output["directory"]
    if (
        directory["device"] != directory_info.st_dev
        or directory["inode"] != directory_info.st_ino
        or directory["nlink"] != directory_info.st_nlink
        or output["parent"]["device"] != parent_info.st_dev
        or output["parent"]["inode"] != parent_info.st_ino
        or output["parent"]["mode"] != stat.S_IMODE(parent_info.st_mode)
    ):
        fail("retained output directory/parent receipt differs")
    for name in OUTPUT_NAMES:
        expected = output["files"][name]
        info = file_infos[name]
        if (
            expected["device"] != info.st_dev
            or expected["inode"] != info.st_ino
            or expected["nlink"] != info.st_nlink
        ):
            fail(f"retained summary output identity differs for {name}")
    actual_artifacts = {
        name: retained_artifact_receipt(file_infos[name], data[name])
        for name in ("raw.jsonl", "stderr.txt", "summary.json", "thermal.csv")
    }
    if not canonical_equal(actual_artifacts, controller["artifacts"]):
        fail("retained controller artifact receipts differ")
    owner_uid = claim_document["uid"]
    owner_gid = claim_document["gid"]
    if (
        directory_info.st_uid != owner_uid or directory_info.st_gid != owner_gid
        or any(
            info.st_uid != owner_uid or info.st_gid != owner_gid
            for info in file_infos.values()
        )
    ):
        fail("retained output ownership differs from the claim owner")


def emit_retained_verification(outcome: str) -> int:
    codes = {"pass": 0, "reject": 2, "invalid": 3, "absent": 4}
    if outcome not in codes:
        fail("retained verification classification differs")
    document = {
        "outcome": outcome, "schema": VERIFY_RESULT_SCHEMA,
        "status": "absent" if outcome == "absent" else "verified",
    }
    payload = canonical_bytes(document) + b"\n"
    sys.stdout.flush()
    view = memoryview(payload)
    while view:
        written = os.write(1, view)
        if written <= 0:
            fail("retained verification classification write was short")
        view = view[written:]
    return codes[outcome]


def inspect_retained(expected_argv_sha256: str) -> str:
    lower_hash(expected_argv_sha256, "expected retained run argv SHA-256")
    opath = getattr(os, "O_PATH", 0)
    directory_flag = getattr(os, "O_DIRECTORY", 0)
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    noatime = getattr(os, "O_NOATIME", 0)
    nonblock = getattr(os, "O_NONBLOCK", 0)
    cloexec = getattr(os, "O_CLOEXEC", 0)
    if (
        not opath or not directory_flag or not nofollow or not noatime
        or not nonblock
    ):
        fail(
            "retained verification requires Linux O_PATH/O_NOATIME/"
            "O_NOFOLLOW/O_NONBLOCK"
        )
    parent_fd = -1
    output_fd = -1
    claim_fd = -1
    file_fds: Dict[str, int] = {}
    close_errors: List[str] = []
    classification: Optional[str] = None
    try:
        parent_fd = os.open(
            "/var/tmp", opath | directory_flag | nofollow | cloexec
        )
        parent_info = os.fstat(parent_fd)
        parent_named = os.stat("/var/tmp", follow_symlinks=False)
        if (
            not stat.S_ISDIR(parent_info.st_mode)
            or (parent_info.st_dev, parent_info.st_ino)
            != (parent_named.st_dev, parent_named.st_ino)
            or stat.S_IMODE(parent_info.st_mode) != 0o1777
            or parent_info.st_uid != 0 or parent_info.st_gid != 0
        ):
            fail("retained /var/tmp parent policy differs")
        claim_named = retained_optional_stat(parent_fd, FIXED_CLAIM_PATH.name)
        output_named = retained_optional_stat(parent_fd, FIXED_OUTPUT_DIR.name)
        if claim_named is None and output_named is None:
            if (
                retained_optional_stat(parent_fd, FIXED_CLAIM_PATH.name) is not None
                or retained_optional_stat(parent_fd, FIXED_OUTPUT_DIR.name) is not None
            ):
                fail("retained namespace appeared during absence classification")
            classification = "absent"
        elif claim_named is None or output_named is None:
            fail("retained namespace is one-sided or partial")
        else:
            if (
                not stat.S_ISREG(claim_named.st_mode)
                or claim_named.st_nlink != 1
                or stat.S_IMODE(claim_named.st_mode) != 0o400
                or not 0 < claim_named.st_size <= 1024 * 1024
                or not stat.S_ISDIR(output_named.st_mode)
                or stat.S_IMODE(output_named.st_mode) != 0o500
            ):
                fail("retained namespace entry policy differs")
            read_flags = (
                os.O_RDONLY | noatime | nofollow | nonblock | cloexec
            )
            claim_fd = os.open(
                FIXED_CLAIM_PATH.name, read_flags, dir_fd=parent_fd
            )
            output_fd = os.open(
                FIXED_OUTPUT_DIR.name, read_flags | directory_flag,
                dir_fd=parent_fd,
            )
            claim_raw, claim_info = retained_read_exact(
                claim_fd, 1024 * 1024, "retained claim"
            )
            if not same_file_receipt(claim_info, claim_named):
                fail("retained claim named/held identity differs")
            directory_info = os.fstat(output_fd)
            output_named_after = os.stat(
                FIXED_OUTPUT_DIR.name, dir_fd=parent_fd,
                follow_symlinks=False,
            )
            if (
                not stat.S_ISDIR(directory_info.st_mode)
                or stat.S_IMODE(directory_info.st_mode) != 0o500
                or not same_file_receipt(directory_info, output_named_after)
            ):
                fail("retained output directory policy differs")
            retained_roster(output_fd)
            limits = {
                "raw.jsonl": MAX_STDOUT_BYTES + 1,
                "stderr.txt": MAX_STDERR_BYTES + 1,
                "summary.json": 8 * 1024 * 1024,
                "thermal.csv": MAX_THERMAL_WINDOW_BYTES,
                "controller.json": 8 * 1024 * 1024,
                "COMPLETE": 4096,
            }
            file_data: Dict[str, bytes] = {}
            file_infos: Dict[str, os.stat_result] = {}
            for name in FINAL_OUTPUT_NAMES:
                fd = os.open(name, read_flags, dir_fd=output_fd)
                file_fds[name] = fd
                payload, info = retained_read_exact(fd, limits[name], name)
                named = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
                if not same_file_receipt(info, named):
                    fail(f"retained named/held identity differs for {name}")
                file_data[name] = payload
                file_infos[name] = info

            claim_document = parse_claim_document(claim_raw)
            claim_receipt = {
                "bytes": len(claim_raw), "device": claim_info.st_dev,
                "inode": claim_info.st_ino, "nlink": claim_info.st_nlink,
                "parent": {
                    "device": parent_info.st_dev, "gid": parent_info.st_gid,
                    "inode": parent_info.st_ino,
                    "mode": stat.S_IMODE(parent_info.st_mode),
                    "path": "/var/tmp", "uid": parent_info.st_uid,
                },
                "path": str(FIXED_CLAIM_PATH),
                "sha256": sha256_bytes(claim_raw),
                "stat": stat_receipt(claim_info),
            }
            validate_claim_runtime_binding(
                claim_document, claim_receipt, claim_document, claim_raw
            )
            summary = parse_summary_bytes(file_data["summary.json"])
            controller = parse_controller_document(file_data["controller.json"])
            complete = parse_complete_document(file_data["COMPLETE"])
            run = validate_recorded_run_authority(
                expected_argv_sha256, controller, summary, claim_document
            )
            if claim_info.st_uid != run.expected_controller_uid:
                fail("retained claim owner differs from run authorization")
            validate_controller_bundle_binding(
                controller, summary, file_data["summary.json"],
                claim_receipt, summary["binary"],
            )
            validate_complete_bundle_binding(
                complete, controller, file_data["controller.json"], claim_receipt
            )
            if complete["elapsed_seconds_before_commit"] < controller[
                "controller_elapsed_seconds"
            ]:
                fail("retained COMPLETE predates its controller receipt")
            validate_retained_claim_controller_binding(
                claim_document, claim_receipt, controller, summary
            )
            validate_retained_science(
                file_data["raw.jsonl"], file_data["stderr.txt"],
                file_data["thermal.csv"], summary,
            )
            validate_retained_output_binding(
                parent_info, directory_info, file_infos, file_data,
                claim_document, claim_receipt, summary, controller,
            )
            if summary["health_module_loader"]:
                by_path = {
                    entry["path"]: entry
                    for entry in summary["source_manifest_before"]["entries"]
                }
                for module in summary["health_module_loader"]["modules"]:
                    entry = by_path[module["path"]]
                    if (
                        module["bytes"] != entry["bytes"]
                        or module["sha256"] != entry["sha256"]
                    ):
                        fail("retained health loader/source binding differs")
            else:
                fail("retained claimed outcome lacks its sealed health loader")
            validate_live_retained_verifier(
                controller, controller["source_manifest_before"]
            )

            # Final name/FD/data sandwich immediately before classification.
            current_parent = os.fstat(parent_fd)
            if (
                current_parent.st_dev, current_parent.st_ino,
                stat.S_IMODE(current_parent.st_mode), current_parent.st_uid,
                current_parent.st_gid,
            ) != (
                parent_info.st_dev, parent_info.st_ino,
                stat.S_IMODE(parent_info.st_mode), parent_info.st_uid,
                parent_info.st_gid,
            ):
                fail("retained parent changed during verification")
            retained_roster(output_fd)
            output_final_named = os.stat(
                FIXED_OUTPUT_DIR.name, dir_fd=parent_fd,
                follow_symlinks=False,
            )
            if (
                not same_file_receipt(os.fstat(output_fd), directory_info)
                or not same_file_receipt(directory_info, output_final_named)
            ):
                fail("retained directory changed during verification")
            claim_again, claim_again_info = retained_read_exact(
                claim_fd, 1024 * 1024, "retained claim final"
            )
            if (
                claim_again != claim_raw
                or not same_file_receipt(claim_again_info, claim_info)
                or not same_file_receipt(
                    claim_again_info, os.stat(
                    FIXED_CLAIM_PATH.name, dir_fd=parent_fd,
                    follow_symlinks=False,
                    ),
                )
            ):
                fail("retained claim changed before classification")
            for name in FINAL_OUTPUT_NAMES:
                payload, info = retained_read_exact(
                    file_fds[name], limits[name], "retained final " + name
                )
                if (
                    payload != file_data[name]
                    or not same_file_receipt(info, file_infos[name])
                    or not same_file_receipt(
                        info,
                        os.stat(name, dir_fd=output_fd, follow_symlinks=False),
                    )
                ):
                    fail(f"retained {name} changed before classification")
            classification = controller["outcome"]
    finally:
        for fd in list(file_fds.values()) + [claim_fd, output_fd, parent_fd]:
            if fd < 0:
                continue
            try:
                os.close(fd)
            except OSError as exc:
                close_errors.append(exception_text(exc))
    if close_errors:
        fail("retained verification descriptor close failed: " + close_errors[0])
    if classification is None:
        fail("retained verification did not classify the namespace")
    return classification


def verify_retained(expected_argv_sha256: str) -> int:
    return emit_retained_verification(inspect_retained(expected_argv_sha256))


def synthetic_config_v2(
    cpu: int, commit: str,
) -> Tuple[Dict[str, Any], Dict[Tuple[int, int], Dict[str, str]]]:
    cells: List[Dict[str, Any]] = []
    receipts: Dict[Tuple[int, int], Dict[str, str]] = {}
    for index, (k, block_bytes) in enumerate(CELLS):
        configuration = {
            "K": k,
            "block_bytes": block_bytes,
            "dense_anchor_layout": index,
            "dense_identity_corner": bool(index & 1),
            "dense_rows": 12 + index,
            "heavy_family": index % 3,
            "heavy_rows": 6 + index,
            "mix_count": 3 + index,
            "packet_attempt": index,
            "packet_peel_seed": 100 + index,
            "precode_attempt": index + 1,
            "precode_seed": (1 << 63) + k + block_bytes,
            "source_hits": 2 + index,
            "staircase": 1,
        }
        source_hash = sha256_bytes(
            f"synthetic-source:{k}:{block_bytes}".encode("ascii")
        )
        cell: Dict[str, Any] = {
            "K": k,
            "block_bytes": block_bytes,
            "construction_equivalent": True,
            "direct_hit_oracle_sha256": "0" * 64,
            "direct_hit_oracle_verified": True,
            "equation_configuration": configuration,
            "equation_configuration_sha256": sha256_bytes(
                canonical_bytes(configuration)
            ),
            "first_repair_sha256": sha256_bytes(
                f"synthetic-first:{k}:{block_bytes}".encode("ascii")
            ),
            "full_repair_oracle_verified": True,
            "full_repair_sha256": sha256_bytes(
                f"synthetic-full:{k}:{block_bytes}".encode("ascii")
            ),
            "high_id_oracle_sha256": sha256_bytes(
                f"synthetic-high:{k}:{block_bytes}".encode("ascii")
            ),
            "high_id_oracle_verified": True,
            "invocations_by_replicate": [
                invocations_for(k, replicate)
                for replicate in range(REPLICATES)
            ],
            "solved_state_equivalence_receipt_sha256": "0" * 64,
            "solved_state_equivalent": True,
            "solved_state_scope": SOLVED_STATE_SCOPE,
            "source_seed": expected_source_seed(k, block_bytes),
            "source_sha256": source_hash,
            "systematic_oracle_sha256": source_hash,
            "systematic_oracle_verified": True,
        }
        cell["direct_hit_oracle_sha256"] = sha256_bytes(
            canonical_bytes(direct_hit_oracle_preimage(cell))
        )
        cell["solved_state_equivalence_receipt_sha256"] = sha256_bytes(
            canonical_bytes(solved_state_preimage(cell))
        )
        receipts[(k, block_bytes)] = {
            field: cell[field] for field in EXPECTED_RECEIPT_FIELDS
        }
        cells.append(cell)
    descriptors = []
    for mode, emission in (
        ("equation", "equation-eval-v1"),
        ("retained", "retained-source-direct-v1"),
    ):
        value = descriptor(SCOPE, emission)
        descriptors.append({
            "descriptor": value,
            "descriptor_sha256": descriptor_hash(value),
            "mode": mode,
        })
    config = {
        "campaign": CAMPAIGN,
        "cells": cells,
        "comparisons": [dict(value) for value in EXPECTED_COMPARISON_RECEIPTS],
        "descriptors": descriptors,
        "expected_invocations": EXPECTED_TOTAL_INVOCATIONS,
        "expected_measured_invocations": EXPECTED_TIMED_INVOCATIONS,
        "expected_panels": EXPECTED_PANEL_COUNT,
        "expected_records": EXPECTED_RECORD_COUNT,
        "expected_warmup_invocations": EXPECTED_WARMUP_INVOCATIONS,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "minimum_invocations": MINIMUM_INVOCATIONS,
        "panel_key_schema": PANEL_KEY_SCHEMA,
        "panel_replicates": REPLICATES,
        "schema": CONFIG_SCHEMA,
        "scope": SCOPE,
        "source_git_commit": commit,
        "source_seed_base": SOURCE_SEED_BASE,
        "target_cpu": cpu,
    }
    return config, receipts


def synthetic_transcript_v2(
    cpu: int, commit: str,
) -> Tuple[bytes, Dict[Tuple[int, int], Dict[str, str]]]:
    config, receipts = synthetic_config_v2(cpu, commit)
    hashes = {
        (SCOPE, entry["mode"]): entry["descriptor_sha256"]
        for entry in config["descriptors"]
    }
    records: List[Dict[str, Any]] = [config]
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for comparison in COMPARISONS:
                order = panel_order(k, block_bytes, comparison, SCOPE, replicate)
                sides = expected_sides(order)
                left_mode, right_mode = comparison_modes(comparison)
                left_elapsed = (
                    900 if comparison == "candidate-over-baseline" else 1000
                )
                right_elapsed = 1000
                n = invocations_for(k, replicate)
                counts = ((n + 1) // 2,) * 4 + (n // 2,) * 4
                records.append({
                    "K": k,
                    "block_bytes": block_bytes,
                    "comparison": comparison,
                    "invocations_per_slot": n,
                    "left_descriptor_sha256": hashes[(SCOPE, left_mode)],
                    "left_direct_systematic_packets": (
                        k if left_mode == "retained" else 0
                    ),
                    "left_outcome_code": 0,
                    "order": order,
                    "panel_key_sha256": panel_key_sha256(
                        k, block_bytes, comparison, SCOPE
                    ),
                    "primary_metric": SCOPE,
                    "replicate": replicate,
                    "right_descriptor_sha256": hashes[(SCOPE, right_mode)],
                    "right_direct_systematic_packets": (
                        k if right_mode == "retained" else 0
                    ),
                    "right_outcome_code": 0,
                    "schema": PANEL_SCHEMA,
                    "scope": SCOPE,
                    "slots": [
                        {
                            "elapsed_ns": (
                                left_elapsed if side == "left" else right_elapsed
                            ),
                            "invocation_count": count,
                            "side": side,
                        }
                        for side, count in zip(sides, counts)
                    ],
                    "target_cpu": cpu,
                })
    records.append({
        "invocation_count": EXPECTED_TOTAL_INVOCATIONS,
        "measured_invocation_count": EXPECTED_TIMED_INVOCATIONS,
        "panel_count": EXPECTED_PANEL_COUNT,
        "record_count": EXPECTED_RECORD_COUNT,
        "schema": TERMINAL_SCHEMA,
        "status": "complete",
        "warmup_invocation_count": EXPECTED_WARMUP_INVOCATIONS,
    })
    return (
        b"".join(canonical_bytes(record) + b"\n" for record in records),
        receipts,
    )


def synthetic_health_fixture() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    child_start = 1_200_000_000
    child_reap = 1_800_000_000
    binary_hash = "a" * 64
    fake_stat = {
        "device": 1, "gid": 1000, "inode": 2, "mode": 0o555,
        "mtime_ns": 3, "nlink": 1, "size": 123, "uid": 1000,
    }
    binary = {
        "capture_error": "none",
        "child_affinity_after_spawn": [EXPECTED_TARGET_CPU],
        "child_pid": 123,
        "child_process_start_ticks": 456,
        "child_reap_monotonic_ns": child_reap,
        "child_start_monotonic_ns": child_start,
        "execution_finished_monotonic_ns": 2_100_000_000,
        "execution_started_monotonic_ns": 1_100_000_000,
        "expected_sha256": binary_hash,
        "output_overflow": False,
        "path": "/tmp/synthetic-worker",
        "path_stat_after": dict(fake_stat),
        "process_started": True,
        "sha256_after": binary_hash,
        "sha256_before": binary_hash,
        "stat_after": dict(fake_stat),
        "stat_before": dict(fake_stat),
        "timed_out": False,
    }
    evidence_parent = "/tmp/synthetic-sampler-evidence"
    csv_path = evidence_parent + "/thermal.csv"
    pid_path = evidence_parent + "/sampler.pid"
    validation_path = evidence_parent + "/validation.jsonl"
    receipt_path = evidence_parent + "/receipt.json"
    script_path = str(
        Path(__file__).resolve(strict=True).parents[1] / SAMPLER_SOURCE_PATH
    )
    header = ",".join(THERMAL_HEADER) + "\n"
    rows = []
    for utc, monotonic_value in (
        ("2026-01-01T00:00:01Z", "1.000000"),
        ("2026-01-01T00:00:02Z", "2.000000"),
    ):
        row = [
            utc, monotonic_value, "50.0", "3000.0", "50.0",
            *("40.0" for _ in range(8)), "0", "0.1", "0.1", "0.1",
            "0", "0",
        ]
        rows.append(",".join(row) + "\n")
    window_ascii = header + "".join(rows)
    window_bytes = window_ascii.encode("ascii")
    admission_csv_bytes = (header + rows[0]).encode("ascii")
    validation_header = {
        "expected_output_owner_uid": EXPECTED_CAMPAIGN_UID,
        "raw_columns": list(THERMAL_HEADER),
        "sampler_source_expected_sha256": "d" * 64,
        "sampling": {
            "dimm_attempts": 5, "dimm_retry_delay_s": 0.01,
            "interval_s": 1.0,
        },
        "schema": THERMAL_VALIDATION_STREAM_SCHEMA,
        "thresholds": dict(THERMAL_SAMPLER_THRESHOLDS),
    }
    validation_header_raw = canonical_bytes(validation_header) + b"\n"

    def validation_sample(index: int, timestamp: float) -> Dict[str, Any]:
        has_previous = index > 0
        return {
            "consecutive_fault_rows": 0,
            "decision": "continue",
            "edac_ce_delta": 0,
            "edac_ue_delta": 0,
            "fault_count": 0,
            "hot_sensors": [],
            "monotonic_s": timestamp,
            "read_error_count": 0,
            "sample_index": index,
            "schema": THERMAL_VALIDATION_SAMPLE_SCHEMA,
            "sensors": {
                name: {
                    "attempt_errors": 0,
                    "hot": False,
                    "hot_streak": 0,
                    "jump_c": 0.0 if has_previous else None,
                    "rate_c_per_s": 0.0 if has_previous else None,
                    "raw_c": 40.0,
                    "reason": "ok",
                    "valid": True,
                }
                for name in VALIDATION_DIMM_FIELDS
            },
        }

    validation_records = [
        canonical_bytes(validation_sample(index, timestamp)) + b"\n"
        for index, timestamp in enumerate((1.0, 2.0))
    ]
    validation_bytes = validation_header_raw + b"".join(validation_records)
    admission_validation_bytes = validation_header_raw + validation_records[0]

    def synthetic_stat(
        device: int, inode: int, mode: int, size: int,
        mtime_ns: int, *, uid: int = 1000, gid: int = 1000,
        nlink: int = 1,
    ) -> Dict[str, int]:
        return {
            "device": device, "gid": gid, "inode": inode, "mode": mode,
            "mtime_ns": mtime_ns, "nlink": nlink, "size": size, "uid": uid,
        }

    def artifact(
        path: str, device: int, inode: int, mode: int, data: bytes,
        mtime_ns: int,
    ) -> Dict[str, Any]:
        return {
            "bytes": len(data), "path": path, "sha256": sha256_bytes(data),
            "stat": synthetic_stat(
                device, inode, mode, len(data), mtime_ns,
            ),
        }

    process_security = {
        "cap_ambient": "0000000000000000",
        "cap_bounding": "0000000000000000",
        "cap_effective": "0000000000000000",
        "cap_inheritable": "0000000000000000",
        "cap_permitted": "0000000000000000",
        "gids": [EXPECTED_CAMPAIGN_GID] * 4,
        "groups": [EXPECTED_SAMPLER_I2C_GID],
        "no_new_privs": 1,
        "schema": PROCESS_SECURITY_SCHEMA,
        "uids": [EXPECTED_CAMPAIGN_UID] * 4,
    }
    cmdline_argv = [
        EXPECTED_SAMPLER_PYTHON, *EXPECTED_SAMPLER_PYTHON_FLAGS, script_path,
        "--csv", csv_path, "--pid-file", pid_path,
        "--validation-jsonl", validation_path, "--receipt", receipt_path,
        "--expected-source-sha256", "d" * 64,
        "--expected-output-owner-uid", str(EXPECTED_CAMPAIGN_UID),
        "--interval", EXPECTED_SAMPLER_INTERVAL_TEXT,
        "--dimm-attempts", EXPECTED_SAMPLER_DIMM_ATTEMPTS_TEXT,
        "--dimm-retry-delay", EXPECTED_SAMPLER_DIMM_RETRY_DELAY_TEXT,
    ]
    pid_bytes = b"124\n"
    sampler = {
        "cmdline_argv": cmdline_argv,
        "cmdline_sha256": sha256_bytes(
            b"".join(os.fsencode(item) + b"\0" for item in cmdline_argv)
        ),
        "cpu": EXPECTED_SAMPLER_CPU,
        "csv_bytes": len(window_bytes),
        "csv_device": 3,
        "csv_inode": 4,
        "csv_path": csv_path,
        "csv_sha256": sha256_bytes(window_bytes),
        "csv_stat": synthetic_stat(3, 4, 0o600, len(window_bytes), 6),
        "environ_sha256": "e" * 64,
        "environment": dict(SEALED_LAUNCH_ENVIRONMENT),
        "environment_sha256": sha256_bytes(
            canonical_bytes(SEALED_LAUNCH_ENVIRONMENT)
        ),
        "evidence_parent": {
            "path": evidence_parent,
            "stat": synthetic_stat(
                3, 30, 0o700, 4096, 5, nlink=2,
            ),
        },
        "executable_device": 22,
        "executable_inode": 23,
        "executable_path": "/usr/bin/python3.12",
        "executable_sha256": "6" * 64,
        "executable_stat": synthetic_stat(
            22, 23, 0o755, 8192, 1, uid=0, gid=0, nlink=2,
        ),
        "pid": 124,
        "pid_file": artifact(pid_path, 3, 31, 0o444, pid_bytes, 5),
        "process_affinity": [EXPECTED_SAMPLER_CPU],
        "process_gid": EXPECTED_CAMPAIGN_GID,
        "process_security": process_security,
        "process_start_ticks": 457,
        "process_uid": EXPECTED_CAMPAIGN_UID,
        "receipt_file": artifact(receipt_path, 3, 33, 0o444, b"", 5),
        "schema": SAMPLER_SCHEMA,
        "script_device": 7,
        "script_inode": 8,
        "script_path": script_path,
        "script_sha256": "d" * 64,
        "script_stat": synthetic_stat(7, 8, 0o444, 321, 9),
        "terminal_status": "ok",
        "validation_header_ascii": validation_header_raw.decode("ascii"),
        "validation_jsonl": artifact(
            validation_path, 3, 32, 0o444, validation_bytes, 6,
        ),
        "window_end_monotonic_ns": 2_000_000_000,
        "window_start_monotonic_ns": 1_000_000_000,
    }
    sampler_admission = json.loads(canonical_bytes(sampler).decode("ascii"))
    sampler_admission.update({
        "csv_bytes": len(admission_csv_bytes),
        "csv_sha256": sha256_bytes(admission_csv_bytes),
        "window_end_monotonic_ns": 0,
        "window_start_monotonic_ns": 0,
    })
    sampler_admission["csv_stat"].update({
        "mtime_ns": 5, "size": len(admission_csv_bytes),
    })
    sampler_admission["validation_jsonl"] = artifact(
        validation_path, 3, 32, 0o444, admission_validation_bytes, 5,
    )
    thermal_summary = summarize_thermal_window_bytes(
        window_bytes, sampler["window_start_monotonic_ns"],
        sampler["window_end_monotonic_ns"],
    )
    validation_summary = summarize_validation_window_bytes(
        validation_bytes, window_bytes, sampler,
    )
    thermal = {
        "cpu": EXPECTED_SAMPLER_CPU,
        "csv_device": sampler["csv_device"],
        "csv_inode": sampler["csv_inode"],
        "csv_path": sampler["csv_path"],
        "pid": sampler["pid"],
        "process_start_ticks": sampler["process_start_ticks"],
        "script_path": sampler["script_path"],
        "script_sha256": sampler["script_sha256"],
        "terminal_status": "complete",
        "validation_attempt_errors_total": validation_summary[
            "attempt_errors_total"
        ],
        "validation_device": sampler["validation_jsonl"]["stat"]["device"],
        "validation_failures": validation_summary["failures"],
        "validation_inode": sampler["validation_jsonl"]["stat"]["inode"],
        "validation_jsonl_ascii": validation_bytes.decode("ascii"),
        "validation_jsonl_bytes": len(validation_bytes),
        "validation_jsonl_sha256": sha256_bytes(validation_bytes),
        "validation_path": validation_path,
        "validation_sample_count": validation_summary["sample_count"],
        "window_csv_ascii": window_ascii,
        "window_csv_bytes": len(window_bytes),
        "window_csv_sha256": sha256_bytes(window_bytes),
        "window_end_monotonic_ns": sampler["window_end_monotonic_ns"],
        "window_start_monotonic_ns": sampler["window_start_monotonic_ns"],
        **thermal_summary,
    }

    def tick(non_idle: int, read_ns: int, idle: int) -> Dict[str, Any]:
        return {
            "cpu": 56,
            "non_idle_ticks": non_idle,
            "read_monotonic_ns": read_ns,
            "tick_fields": {
                "idle": idle, "iowait": 0, "irq": 0, "nice": 0,
                "softirq": 0, "steal": 0, "system": 0, "user": non_idle,
            },
        }

    admission = unbounded_sibling_tick_receipt(
        tick(10, 1_050_000_000, 100), tick(11, 1_100_000_000, 101),
        1_050_000_000, 1_100_000_000,
    )
    during = unbounded_sibling_tick_receipt(
        tick(20, 1_100_000_000, 110), tick(21, 1_900_000_000, 111),
        child_start, child_reap,
    )
    health = {
        "admission_sibling_ticks": admission,
        "child_reap_monotonic_ns": child_reap,
        "child_start_monotonic_ns": child_start,
        "collection_failures": [],
        "controller_core": list(EXPECTED_CONTROLLER_CORE),
        "controller_cpu": EXPECTED_CONTROLLER_CPU,
        "controller_initial_affinity": [
            EXPECTED_TARGET_CPU, EXPECTED_CONTROLLER_CPU, EXPECTED_SAMPLER_CPU,
        ],
        "controller_singleton_affinity_end": [EXPECTED_CONTROLLER_CPU],
        "edac_policy": "ce-and-ue-every-sample-zero-v1",
        "evidence_status": "complete",
        "receipt_sha256": None,
        "sampler": sampler,
        "sampler_admission": sampler_admission,
        "sampler_admission_receipt_sha256": sha256_bytes(
            canonical_bytes(sampler_admission)
        ),
        "sampler_core": list(EXPECTED_SAMPLER_CORE),
        "sampler_cpu": EXPECTED_SAMPLER_CPU,
        "sampler_receipt_sha256": sha256_bytes(canonical_bytes(sampler)),
        "sampler_terminal": {},
        "sampler_terminal_receipt_sha256": None,
        "schema": HEALTH_SCHEMA,
        "sibling_non_idle_tick_cap": SIBLING_NON_IDLE_TICK_CAP,
        "sibling_tick_policy": (
            "linux-proc-stat-user-nice-system-irq-softirq-steal-v1"
        ),
        "sibling_ticks": [during],
        "target_core": list(EXPECTED_TARGET_CORE),
        "target_cpu": EXPECTED_TARGET_CPU,
        "target_threads": list(EXPECTED_TARGET_THREADS),
        "thermal": thermal,
        "thermal_max_millic": THERMAL_MAX_MILLIC,
        "terminal_status": "ok",
        "violations": [],
    }
    return finalize_health_receipt(health, [], []), binary


def synthetic_terminal_health_fixture(
    health: Mapping[str, Any], binary: Mapping[str, Any],
    *, terminal_signal: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a fully bound invalid EDAC-abort receipt for selftest coverage."""
    terminal_health = json.loads(canonical_bytes(health).decode("ascii"))
    admission = terminal_health["sampler_admission"]
    raw_lines = terminal_health["thermal"][
        "window_csv_ascii"
    ].encode("ascii").splitlines(keepends=True)
    validation_lines = terminal_health["thermal"][
        "validation_jsonl_ascii"
    ].encode("ascii").splitlines(keepends=True)
    final_raw_row = parse_csv_physical_line(
        raw_lines[-1], "synthetic terminal raw row"
    )
    final_raw_row[17] = "1"
    raw_lines[-1] = (",".join(final_raw_row) + "\n").encode("ascii")
    raw_csv = b"".join(raw_lines)
    final_validation = parse_canonical_json_line(
        validation_lines[-1], "synthetic terminal validation record"
    )
    final_validation["decision"] = "edac_abort"
    final_validation["edac_ce_delta"] = 1
    validation_lines[-1] = canonical_bytes(final_validation) + b"\n"
    validation_jsonl = b"".join(validation_lines)
    summary, final_decision = replay_terminal_validation(
        raw_csv, validation_jsonl, admission
    )
    if final_decision != "edac_abort":
        raise AssertionError("synthetic terminal decision did not replay")

    def binding_from_artifact(
        artifact: Mapping[str, Any], data: bytes, mode: int,
        *, nlink: int = 1,
    ) -> Dict[str, Any]:
        info = artifact["stat"]
        return {
            "device": info["device"], "gid": info["gid"],
            "inode": info["inode"], "mode": mode, "nlink": nlink,
            "sha256": sha256_bytes(data), "size": len(data),
            "uid": info["uid"],
        }

    parent_stat = admission["evidence_parent"]["stat"]
    parent_binding = {
        name: parent_stat[name] for name in SAMPLER_PARENT_BINDING_KEYS
    }

    def destination(path: str, *, source: bool = False) -> Dict[str, Any]:
        return {
            "basename": Path(path).name,
            "expected_owner_uid": (
                None if source else admission["process_uid"]
            ),
            "parent": {
                "binding": (
                    {
                        "device": 7, "gid": admission["process_gid"],
                        "inode": 70, "mode": 0o755, "nlink": 2,
                        "uid": admission["process_uid"],
                    }
                    if source else dict(parent_binding)
                ),
                "path": str(Path(path).parent),
            },
            "path": path,
        }

    raw_binding = binding_from_artifact(
        {"stat": admission["csv_stat"]}, raw_csv, 0o444
    )
    validation_binding = binding_from_artifact(
        admission["validation_jsonl"], validation_jsonl, 0o444
    )
    pid_data = (str(admission["pid"]) + "\n").encode("ascii")
    pid_binding = binding_from_artifact(
        admission["pid_file"], pid_data, 0o444
    )
    pid_held_binding = dict(pid_binding)
    pid_held_binding["nlink"] = 0
    receipt_reservation_binding = binding_from_artifact(
        admission["receipt_file"], b"", 0o444
    )
    script_stat = admission["script_stat"]
    source_binding = {
        "device": script_stat["device"], "gid": script_stat["gid"],
        "inode": script_stat["inode"], "mode": script_stat["mode"],
        "nlink": script_stat["nlink"],
        "sha256": admission["script_sha256"],
        "size": script_stat["size"], "uid": script_stat["uid"],
    }
    producer_receipt: Dict[str, Any] = {
        "argv": admission["cmdline_argv"][5:],
        "cpu_tctl_max_c": 50.0,
        "edac_ce_paths": ["/sys/devices/system/edac/ce_count"],
        "edac_ue_paths": ["/sys/devices/system/edac/ue_count"],
        "errors": [],
        "exit_code": 6,
        "finished_monotonic_ns": 2_070_000_000,
        "finished_utc": "2026-01-01T00:00:03Z",
        "expected_output_owner_uid": admission["process_uid"],
        "outcome": "edac_abort",
        "pid": admission["pid"],
        "pid_file": {
            "binding": pid_binding,
            "destination": destination(admission["pid_file"]["path"]),
            "path": admission["pid_file"]["path"], "removed": True,
        },
        "raw_csv": {
            "binding": raw_binding,
            "destination": destination(admission["csv_path"]),
            "path": admission["csv_path"],
        },
        "raw_samples_preserved": True,
        "receipt_file": {
            "destination": destination(admission["receipt_file"]["path"]),
            "path": admission["receipt_file"]["path"],
            "reservation_binding": receipt_reservation_binding,
        },
        "sampler_source": {
            "binding": source_binding,
            "binding_finished": dict(source_binding),
            "destination": destination(admission["script_path"], source=True),
            "expected_sha256": admission["script_sha256"],
            "path": admission["script_path"],
            "sha256": admission["script_sha256"],
        },
        "sampling": {
            "dimm_attempts": 5, "dimm_retry_delay_s": 0.01,
            "interval_s": 1.0,
        },
        "schema": THERMAL_SAMPLER_SCHEMA,
        "signal": terminal_signal,
        "started_monotonic_ns": 500_000_000,
        "started_utc": "2026-01-01T00:00:00Z",
        "summary": summary,
        "thresholds": dict(THERMAL_SAMPLER_THRESHOLDS),
        "validation_jsonl": {
            "binding": validation_binding,
            "destination": destination(
                admission["validation_jsonl"]["path"]
            ),
            "path": admission["validation_jsonl"]["path"],
        },
    }
    producer_receipt["self_sha256_excluding_field"] = sha256_bytes(
        canonical_bytes(producer_receipt) + b"\n"
    )
    producer_receipt_ascii = canonical_bytes(producer_receipt) + b"\n"
    producer_binding = binding_from_artifact(
        admission["receipt_file"], producer_receipt_ascii, 0o444
    )
    stream_parts = split_terminal_streams(
        raw_csv, validation_jsonl, allow_unpaired=False
    )
    child_start = binary["child_start_monotonic_ns"]
    child_reap = binary["child_reap_monotonic_ns"]
    terminal = {
        "admission_receipt_sha256": sha256_bytes(canonical_bytes(admission)),
        "coverage": {
            "child_reap_monotonic_ns": child_reap,
            "child_start_monotonic_ns": child_start,
            "coverage_shortfall_ns": 0,
            "covers_child_interval": True,
            "window_end_monotonic_ns": 2_000_000_000,
            "window_start_monotonic_ns": 1_000_000_000,
        },
        "pid_file_held_binding": pid_held_binding,
        "process_exit_observation": "linux-pidfd-readable-nonchild-v1",
        "process_exit_observed_monotonic_ns": 2_080_000_000,
        "producer_receipt_ascii": producer_receipt_ascii.decode("ascii"),
        "producer_receipt_binding": producer_binding,
        "raw_csv_ascii": raw_csv.decode("ascii"),
        "raw_csv_binding": raw_binding,
        "schema": SAMPLER_TERMINAL_SCHEMA,
        "stream_suffixes": stream_parts["suffix_receipt"],
        "terminal_status": "invalid",
        "validation_jsonl_ascii": validation_jsonl.decode("ascii"),
        "validation_jsonl_binding": validation_binding,
        "window_csv_ascii": raw_csv.decode("ascii"),
        "window_validation_jsonl_ascii": validation_jsonl.decode("ascii"),
    }
    validate_sampler_terminal_receipt(
        terminal, admission, child_start, child_reap
    )
    terminal_health["sampler"] = {}
    terminal_health["sampler_receipt_sha256"] = None
    terminal_health["sampler_terminal"] = terminal
    terminal_health["sampler_terminal_receipt_sha256"] = sha256_bytes(
        canonical_bytes(terminal)
    )
    terminal_health["thermal"] = {}
    terminal_health = finalize_health_receipt(
        terminal_health,
        ["sampler endpoint: sampler-terminal-decision:edac_abort"], [],
    )
    validate_health_receipt(terminal_health, binary, require_complete=False)
    return terminal_health


def synthetic_final_bundle_fixture(
    health: Mapping[str, Any], binary: Mapping[str, Any], commit: str,
) -> Dict[str, Any]:
    def receipt_stat(
        device: int, inode: int, mode: int, size: int,
        uid: int = 1000, gid: int = 1000, nlink: int = 1,
    ) -> Dict[str, int]:
        return {
            "device": device, "gid": gid, "inode": inode, "mode": mode,
            "mtime_ns": 1, "nlink": nlink, "size": size, "uid": uid,
        }

    sampler_authority = health["sampler_admission"]
    git_sha = "9" * 64
    git_stat = receipt_stat(20, 21, 0o755, 4096, 0, 0)
    source_blob_oids = []
    for path in SOURCE_PATHS:
        oid = git_blob_oid(("synthetic-source:" + path).encode("ascii"))
        source_blob_oids.append({
            "head_oid": oid, "path": path, "worktree_oid": oid,
        })
    git = {
        "executable": {
            "path": str(GIT_EXECUTABLE), "sha256": git_sha,
            "stat": git_stat,
        },
        "head": commit,
        "source_blob_oids": source_blob_oids,
        "source_blob_roster_sha256": sha256_bytes(
            canonical_bytes(source_blob_oids)
        ),
        "source_roster_sha256": sha256_bytes(b"".join(
            (path + "\n").encode("ascii") for path in SOURCE_PATHS
        )),
        "tracked_index_flags_sha256": "8" * 64,
        "tracked_status_sha256": sha256_bytes(b""),
        "worktree_root": str(Path(__file__).resolve(strict=True).parents[1]),
    }
    script_relative = "bench/Wh2DirectSystematicComplementScreen.py"
    script_sha = "7" * 64
    script_stat = receipt_stat(24, 25, 0o444, 321)
    entries = []
    manifest_preimage = bytearray()
    for index, path in enumerate(SOURCE_PATHS):
        if path == script_relative:
            digest = script_sha
            size = 321
            source_stat = script_stat
        elif path == SAMPLER_SOURCE_PATH:
            digest = sampler_authority["script_sha256"]
            size = sampler_authority["script_stat"]["size"]
            source_stat = dict(sampler_authority["script_stat"])
        else:
            digest = sha256_bytes(("synthetic-source:" + path).encode("ascii"))
            size = index + 1
            source_stat = receipt_stat(24, 100 + index, 0o444, size)
        entries.append({
            "bytes": size,
            "git_blob_oid": source_blob_oids[index]["head_oid"],
            "path": path, "sha256": digest, "stat": source_stat,
        })
        manifest_preimage.extend(f"{digest}  {path}\n".encode("ascii"))
    source = {
        "entries": entries,
        "sha256": sha256_bytes(bytes(manifest_preimage)),
    }
    argv = [str(Path(__file__).resolve(strict=True)), "--run-once"]
    controller = {
        "argv": argv,
        "argv_sha256": sha256_bytes(
            b"".join(os.fsencode(item) + b"\0" for item in argv)
        ),
        "controller_cpu": EXPECTED_CONTROLLER_CPU,
        "dont_write_bytecode": True,
        "environment": dict(SEALED_LAUNCH_ENVIRONMENT),
        "git_path": str(GIT_EXECUTABLE),
        "git_sha256": git_sha,
        "git_stat": git_stat,
        "isolated": True,
        "optimize": 0,
        "pid": 111,
        "process_gid": EXPECTED_CAMPAIGN_GID,
        "process_security": {
            "cap_ambient": "0000000000000000",
            "cap_bounding": "0000000000000000",
            "cap_effective": "0000000000000000",
            "cap_inheritable": "0000000000000000",
            "cap_permitted": "0000000000000000",
            "gids": [EXPECTED_CAMPAIGN_GID] * 4,
            "groups": [],
            "no_new_privs": 1,
            "schema": PROCESS_SECURITY_SCHEMA,
            "uids": [EXPECTED_CAMPAIGN_UID] * 4,
        },
        "process_start_ticks": 222,
        "process_uid": EXPECTED_CAMPAIGN_UID,
        "python_declared_path": "/usr/bin/python3",
        "python_path": "/usr/bin/python3.12",
        "python_sha256": "6" * 64,
        "python_stat": receipt_stat(
            22, 23, 0o755, 8192, 0, 0, nlink=2
        ),
        "receipt_sha256": None,
        "schema": CONTROLLER_PROVENANCE_SCHEMA,
        "script_path": str(Path(__file__).resolve(strict=True)),
        "script_sha256": script_sha,
        "script_stat": script_stat,
        "singleton_affinity": [EXPECTED_CONTROLLER_CPU],
    }
    build_entries = []
    for index, path in enumerate(BUILD_PATHS):
        build_entries.append({
            **{
                "bytes": index + 10, "device": 30, "gid": 1000,
                "inode": 31 + index, "mode": 0o444, "nlink": 1,
                "sha256": sha256_bytes(path.encode("ascii")), "uid": 1000,
            },
            "path": path,
        })
    build = {
        "entries": build_entries,
        "root": "/tmp/synthetic-build",
        "schema": "wirehair.wh2.direct-systematic-complement-build.v2",
        "sha256": None,
    }
    build["sha256"] = sha256_bytes(canonical_bytes(build))
    run_values = (
        binary["path"], build["root"], str(EXPECTED_TARGET_CPU),
        str(EXPECTED_CONTROLLER_CPU), str(sampler_authority["pid"]),
        str(EXPECTED_SAMPLER_CPU), sampler_authority["script_path"],
        sampler_authority["csv_path"], sampler_authority["pid_file"]["path"],
        sampler_authority["validation_jsonl"]["path"],
        sampler_authority["receipt_file"]["path"],
        commit, binary["expected_sha256"],
        str(binary["stat_before"]["uid"]), build["sha256"], script_sha,
        str(script_stat["uid"]), str(script_stat["gid"]),
        git_sha, controller["python_sha256"],
        str(sampler_authority["process_start_ticks"]),
        sampler_authority["script_sha256"],
        str(sampler_authority["csv_device"]),
        str(sampler_authority["csv_inode"]),
        str(sampler_authority["pid_file"]["stat"]["device"]),
        str(sampler_authority["pid_file"]["stat"]["inode"]),
        str(sampler_authority["validation_jsonl"]["stat"]["device"]),
        str(sampler_authority["validation_jsonl"]["stat"]["inode"]),
        str(sampler_authority["receipt_file"]["stat"]["device"]),
        str(sampler_authority["receipt_file"]["stat"]["inode"]),
        sampler_authority["cmdline_sha256"],
        sampler_authority["environ_sha256"],
        sampler_authority["executable_sha256"],
        str(sampler_authority["process_uid"]),
        str(sampler_authority["process_gid"]),
        str(EXPECTED_SAMPLER_I2C_GID), source["sha256"],
    )
    if len(run_values) != len(RUN_ONCE_OPTION_ORDER):
        raise AssertionError("synthetic run authority roster differs")
    argv = [str(Path(__file__).resolve(strict=True)), "--run-once"]
    for option, value in zip(RUN_ONCE_OPTION_ORDER, run_values):
        argv.extend((option, value))
    controller["argv"] = argv
    controller["argv_sha256"] = sha256_bytes(
        b"".join(os.fsencode(item) + b"\0" for item in argv)
    )
    controller["receipt_sha256"] = sha256_bytes(canonical_bytes(controller))
    claim_document = {
        "binary_sha256": binary["expected_sha256"],
        "build_manifest_sha256": build["sha256"],
        "campaign": CAMPAIGN,
        "controller_receipt_sha256": controller["receipt_sha256"],
        "controller_started_monotonic_ns": 1_000_000_000,
        "gid": 1000,
        "git_sha256": git_sha,
        "output_path": str(FIXED_OUTPUT_DIR),
        "parent_device": 40,
        "parent_inode": 41,
        "pid": controller["pid"],
        "process_start_ticks": controller["process_start_ticks"],
        "schema": CLAIM_SCHEMA,
        "source_commit": commit,
        "source_manifest_sha256": source["sha256"],
        "uid": 1000,
    }
    claim_raw = canonical_bytes(claim_document) + b"\n"
    claim_stat = receipt_stat(42, 43, 0o400, len(claim_raw))
    claim_receipt = {
        "bytes": len(claim_raw), "device": 42, "inode": 43, "nlink": 1,
        "parent": {
            "device": 40, "gid": 0, "inode": 41, "mode": 0o1777,
            "path": "/var/tmp", "uid": 0,
        },
        "path": str(FIXED_CLAIM_PATH),
        "sha256": sha256_bytes(claim_raw), "stat": claim_stat,
    }
    output = {
        "directory": {
            "device": 50, "inode": 51, "nlink": 2,
            "reserved_mode": 0o700, "sealed_mode": 0o500,
        },
        "files": {
            name: {
                "device": 50, "inode": 52 + index, "nlink": 1,
                "reserved_mode": 0o600, "sealed_mode": 0o400,
            }
            for index, name in enumerate(OUTPUT_NAMES)
        },
        "parent": {
            "device": 40, "initial_nlink": 2, "inode": 41,
            "mode": 0o1777,
        },
        "path": str(FIXED_OUTPUT_DIR),
    }
    source_by_path = {entry["path"]: entry for entry in source["entries"]}
    health_loader = {
        "dont_write_bytecode": True,
        "isolated": True,
        "modules": [
            {
                "bytes": source_by_path[path]["bytes"],
                "module": module,
                "path": path,
                "sha256": source_by_path[path]["sha256"],
            }
            for module, path in HEALTH_MODULE_SOURCES
        ],
        "optimize": 0,
        "receipt_sha256": None,
        "schema": HEALTH_LOADER_SCHEMA,
    }
    health_loader["receipt_sha256"] = sha256_bytes(
        canonical_bytes(health_loader)
    )
    summary = make_summary_preimage({
        "binary": dict(binary), "config": {}, "elapsed_seconds": 1.0,
        "expected_source_git_commit": commit,
        "expected_source_manifest_sha256": source["sha256"],
        "failure": "synthetic-invalid", "git_after": git,
        "git_before": git, "health": dict(health),
        "health_module_loader": health_loader, "outcome": "invalid",
        "output_bundle": output, "process_exit_code": 0,
        "publication_failures": ["synthetic invalid evidence"], "raw_bytes": 0,
        "raw_complete": False, "raw_record_count": 0,
        "raw_sha256": sha256_bytes(b""), "schema": SUMMARY_SCHEMA,
        "signal": None, "signal_names": [],
        "source_manifest_after": source, "source_manifest_before": source,
        "statistics": {}, "stderr_bytes": 0,
        "stderr_sha256": sha256_bytes(b""),
        "summary_preimage_sha256": None, "target_cpu": EXPECTED_TARGET_CPU,
    })
    summary_bytes = canonical_bytes(summary) + b"\n"
    thermal_bytes = health_thermal_artifact_bytes(health)
    payloads = {
        "raw.jsonl": b"", "stderr.txt": b"",
        "summary.json": summary_bytes, "thermal.csv": thermal_bytes,
    }
    artifacts = {
        name: {
            "bytes": len(payload), "device": 50, "gid": 1000,
            "inode": 52 + index,
            "mode": 0o400, "nlink": 1,
            "sha256": sha256_bytes(payload), "uid": 1000,
        }
        for index, (name, payload) in enumerate(payloads.items())
    }
    controller_document = {
        "artifacts": artifacts, "binary": dict(binary),
        "build_after": build, "build_before": build, "campaign": CAMPAIGN,
        "claim": claim_receipt, "controller_after": controller,
        "controller_before": controller,
        "controller_deadline_seconds": CONTROLLER_DEADLINE_SECONDS,
        "controller_elapsed_seconds": 2.0,
        "controller_observed_monotonic_ns": 3_000_000_000,
        "controller_started_monotonic_ns": 1_000_000_000,
        "failure": "synthetic-invalid", "git_after": git,
        "git_before": git,
        "health_receipt_sha256": health["receipt_sha256"],
        "outcome": "invalid", "outer_deadline_seconds": OUTER_DEADLINE_SECONDS,
        "output_path": str(FIXED_OUTPUT_DIR), "receipt_sha256": None,
        "sampler_after": health["sampler"], "schema": CONTROLLER_SCHEMA,
        "signals": [], "source_manifest_after": source,
        "source_manifest_before": source,
        "summary": {
            "outcome": "invalid",
            "preimage_sha256": summary["summary_preimage_sha256"],
            "sha256": sha256_bytes(summary_bytes),
        },
        "target_cpu": EXPECTED_TARGET_CPU,
    }
    controller_bytes = controller_document_bytes(controller_document)
    complete = {
        "campaign": CAMPAIGN, "claim_sha256": claim_receipt["sha256"],
        "complete_observed_monotonic_ns": 4_000_000_000,
        "controller_outcome": "invalid",
        "controller_sha256": sha256_bytes(controller_bytes),
        "elapsed_seconds_before_commit": 3.0, "schema": COMPLETE_SCHEMA,
        "status": "complete",
    }
    return {
        "binary": dict(binary), "build": build, "claim_document": claim_document,
        "claim_raw": claim_raw, "claim_receipt": claim_receipt,
        "complete": complete, "controller": controller_document,
        "controller_bytes": controller_bytes, "git": git, "source": source,
        "run_argv_sha256": controller["argv_sha256"],
        "summary": summary, "summary_bytes": summary_bytes,
    }


def synthetic_preflight_seal_fixture(
    final_fixture: Mapping[str, Any], health: Mapping[str, Any], commit: str,
) -> Dict[str, Any]:
    source = final_fixture["source"]
    git = final_fixture["git"]
    build = final_fixture["build"]
    binary = final_fixture["binary"]
    run_controller = final_fixture["controller"]["controller_before"]
    sampler = health["sampler_admission"]
    preflight_values = (
        binary["path"], build["root"], str(sampler["pid"]),
        sampler["script_path"], sampler["csv_path"],
        sampler["pid_file"]["path"], sampler["validation_jsonl"]["path"],
        sampler["receipt_file"]["path"], commit,
    )
    preflight_argv = [
        str(Path(__file__).resolve(strict=True)), "--preflight-seal",
    ]
    for option, item in zip(PREFLIGHT_SEAL_OPTION_ORDER, preflight_values):
        preflight_argv.extend((option, item))
    image = {
        name: run_controller[name]
        for name in PREFLIGHT_CONTROLLER_IMAGE_KEYS
    }
    preflight_controller: Dict[str, Any] = {
        "argv": preflight_argv,
        "argv_sha256": sha256_bytes(
            b"".join(os.fsencode(item) + b"\0" for item in preflight_argv)
        ),
        "controller_cpu": EXPECTED_CONTROLLER_CPU,
        "dont_write_bytecode": True,
        "environment": dict(SEALED_LAUNCH_ENVIRONMENT),
        "image": image,
        "isolated": True,
        "optimize": 0,
        "pid": run_controller["pid"],
        "process_gid": EXPECTED_CAMPAIGN_GID,
        "process_security": run_controller["process_security"],
        "process_start_ticks": run_controller["process_start_ticks"],
        "process_uid": EXPECTED_CAMPAIGN_UID,
        "receipt_sha256": None,
        "schema": PREFLIGHT_CONTROLLER_SCHEMA,
        "singleton_affinity": [EXPECTED_CONTROLLER_CPU],
    }
    preflight_controller["receipt_sha256"] = sha256_bytes(
        canonical_bytes(preflight_controller)
    )
    binary_observation = {
        "path": binary["path"], "sha256": binary["expected_sha256"],
        "stat": binary["stat_before"],
    }
    prefix = {
        "final_pid_file_stat": sampler["pid_file"]["stat"],
        "final_raw_stat": sampler["csv_stat"],
        "final_receipt_stat": sampler["receipt_file"]["stat"],
        "final_validation_stat": sampler["validation_jsonl"]["stat"],
        "monitor_event": "none",
        "pid_file_sha256": sampler["pid_file"]["sha256"],
        "raw_after_bytes": sampler["csv_bytes"],
        "raw_after_sha256": sampler["csv_sha256"],
        "raw_before_bytes": sampler["csv_bytes"],
        "raw_before_sha256": sampler["csv_sha256"],
        "receipt_bytes": 0,
        "sampler_after_sha256": sha256_bytes(canonical_bytes(sampler)),
        "sampler_before_sha256": sha256_bytes(canonical_bytes(sampler)),
        "validation_after_bytes": sampler["validation_jsonl"]["bytes"],
        "validation_after_sha256": sampler["validation_jsonl"]["sha256"],
        "validation_before_bytes": sampler["validation_jsonl"]["bytes"],
        "validation_before_sha256": sampler["validation_jsonl"]["sha256"],
    }
    value: Dict[str, Any] = {
        "binary_after": binary_observation,
        "binary_before": binary_observation,
        "build_manifest_after": build,
        "build_manifest_before": build,
        "expected_source_commit": commit,
        "git_after": git,
        "git_before": git,
        "health_module_loader": final_fixture["summary"]["health_module_loader"],
        "preflight_controller_after": preflight_controller,
        "preflight_controller_before": preflight_controller,
        "receipt_sha256": None,
        "run_argv": run_controller["argv"],
        "run_argv_sha256": run_controller["argv_sha256"],
        "sampler_after": sampler,
        "sampler_before": sampler,
        "sampler_prefix_binding": prefix,
        "schema": PREFLIGHT_SEAL_SCHEMA,
        "source_manifest_after": source,
        "source_manifest_before": source,
        "source_root": str(Path(__file__).resolve(strict=True).parents[1]),
    }
    value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
    return value


def selftest_v2() -> int:
    cpu = EXPECTED_TARGET_CPU
    commit = "1" * 40
    if panel_key_sha256(32, 1280, "baseline-aa", SCOPE) != (
        GOLDEN_PANEL_KEY_SHA256
    ):
        raise AssertionError("cross-language panel-key golden digest differs")
    raw, synthetic_receipts = synthetic_transcript_v2(cpu, commit)
    global EXPECTED_CELL_RECEIPTS
    frozen_receipts = EXPECTED_CELL_RECEIPTS
    EXPECTED_CELL_RECEIPTS = synthetic_receipts
    try:
        records = parse_transcript(raw)
        _, statistics, failures = validate_transcript(records, cpu, commit)
        if failures:
            raise AssertionError(f"passing transcript failed: {failures}")
        if statistics["student_t_975_df11"] != T11_975:
            raise AssertionError("Student-t constant differs")
        for group in AGGREGATE_GROUPS:
            if not (
                statistics["aggregates"][SCOPE][group]["upper95_log"]
                < math.log(0.99)
            ):
                raise AssertionError("synthetic aggregate did not strictly pass")

        mutations: List[Tuple[str, List[Dict[str, Any]]]] = []
        wrong_count = json.loads(canonical_bytes(records).decode("ascii"))
        wrong_count[-1]["measured_invocation_count"] -= 1
        mutations.append(("terminal measured count", wrong_count))
        wrong_split = json.loads(canonical_bytes(records).decode("ascii"))
        wrong_split[1]["slots"][0]["invocation_count"] += 1
        mutations.append(("slot invocation split", wrong_split))
        wrong_seed = json.loads(canonical_bytes(records).decode("ascii"))
        wrong_seed[0]["cells"][0]["source_seed"] = "0x0"
        mutations.append(("source seed", wrong_seed))
        wrong_direct = json.loads(canonical_bytes(records).decode("ascii"))
        wrong_direct[0]["cells"][0]["direct_hit_oracle_sha256"] = "0" * 64
        mutations.append(("direct-hit receipt", wrong_direct))
        wrong_order = json.loads(canonical_bytes(records).decode("ascii"))
        wrong_order[1]["order"] = (
            "BAAB" if wrong_order[1]["order"] == "ABBA" else "ABBA"
        )
        mutations.append(("panel order", wrong_order))
        for label, mutated in mutations:
            try:
                validate_transcript(mutated, cpu, commit)
            except ValidationError:
                pass
            else:
                raise AssertionError(f"tampered {label} was accepted")

        nonimproving = json.loads(canonical_bytes(records).decode("ascii"))
        for panel in nonimproving[1:-1]:
            if panel["comparison"] == "candidate-over-baseline":
                for slot in panel["slots"]:
                    slot["elapsed_ns"] = 1000
        _, _, failed = validate_transcript(nonimproving, cpu, commit)
        if set(failed) != {
            f"aggregate_upper95:{SCOPE}:{group}" for group in AGGREGATE_GROUPS
        }:
            raise AssertionError("strict 0.99 aggregate gates differ")

        cell_slow = json.loads(canonical_bytes(records).decode("ascii"))
        for panel in cell_slow[1:-1]:
            if (
                panel["comparison"] == "candidate-over-baseline"
                and panel["K"] == CELLS[0][0]
                and panel["block_bytes"] == CELLS[0][1]
            ):
                for slot in panel["slots"]:
                    slot["elapsed_ns"] = 1030 if slot["side"] == "left" else 1000
        _, _, failed = validate_transcript(cell_slow, cpu, commit)
        if f"cell_point:{SCOPE}:candidate-over-baseline:K32:B1280" not in failed:
            raise AssertionError("2% candidate cell point gate differs")

        aa_bias = json.loads(canonical_bytes(records).decode("ascii"))
        for panel in aa_bias[1:-1]:
            if (
                panel["comparison"] == "baseline-aa"
                and panel["K"] == CELLS[0][0]
                and panel["block_bytes"] == CELLS[0][1]
            ):
                for slot in panel["slots"]:
                    slot["elapsed_ns"] = 1030 if slot["side"] == "left" else 1000
        _, _, failed = validate_transcript(aa_bias, cpu, commit)
        if f"aa_ci:{SCOPE}:baseline-aa:K32:B1280" not in failed:
            raise AssertionError("A/A two-sided CI gate differs")

        duplicate = b'{"campaign":"x","campaign":"y"}\n' + raw.split(b"\n", 1)[1]
        try:
            parse_transcript(duplicate)
        except ValidationError:
            pass
        else:
            raise AssertionError("duplicate JSON key was accepted")
    finally:
        EXPECTED_CELL_RECEIPTS = frozen_receipts

    expected_seeds = (
        "0x547000f8f0f19096", "0x547000f8f07995d6",
        "0x547000f8f0799096", "0x547000f8f7619096",
        "0x547000f8e0b195d6", "0x547000f8b0b195d6",
        "0x547000f8b0b19096", "0x547000f86cf195d6",
        "0x547000f86cf19096", "0x547000f9f0b195d6",
        "0x547000f9f0b19096",
    )
    if tuple(expected_source_seed(*cell) for cell in CELLS) != expected_seeds:
        raise AssertionError("frozen source-seed formula differs")
    if (
        sum(4 * invocations_for(k, replicate) * len(COMPARISONS)
            for k, _ in CELLS for replicate in range(REPLICATES))
        != EXPECTED_TIMED_INVOCATIONS
        or EXPECTED_TIMED_INVOCATIONS + EXPECTED_WARMUP_INVOCATIONS
        != EXPECTED_TOTAL_INVOCATIONS
    ):
        raise AssertionError("invocation arithmetic differs")

    health, binary = synthetic_health_fixture()
    validate_binary_receipt(binary, require_complete=True)
    validate_health_receipt(health, binary, require_complete=True)
    stale_health = json.loads(canonical_bytes(health).decode("ascii"))
    stale_sampler = stale_health["sampler"]
    stale_sampler["window_end_monotonic_ns"] = (
        stale_sampler["window_start_monotonic_ns"]
    )
    full_raw_lines = stale_health["thermal"][
        "window_csv_ascii"
    ].encode("ascii").splitlines(keepends=True)
    full_validation_lines = stale_health["thermal"][
        "validation_jsonl_ascii"
    ].encode("ascii").splitlines(keepends=True)
    stale_raw = b"".join(full_raw_lines[:2])
    stale_validation = b"".join(full_validation_lines[:2])
    stale_raw_full = stale_raw + full_raw_lines[2]
    stale_validation_full = (
        stale_validation
        + full_validation_lines[2][:len(full_validation_lines[2]) // 2]
    )
    stale_sampler["csv_bytes"] = len(stale_raw)
    stale_sampler["csv_sha256"] = sha256_bytes(stale_raw)
    stale_sampler["csv_stat"]["size"] = len(stale_raw_full)
    stale_sampler["validation_jsonl"]["bytes"] = len(
        stale_validation_full
    )
    stale_sampler["validation_jsonl"]["sha256"] = sha256_bytes(
        stale_validation_full
    )
    stale_sampler["validation_jsonl"]["stat"]["size"] = len(
        stale_validation_full
    )
    stale_raw_summary = summarize_thermal_window_bytes(
        stale_raw, stale_sampler["window_start_monotonic_ns"],
        stale_sampler["window_end_monotonic_ns"],
    )
    stale_validation_summary = summarize_validation_window_bytes(
        stale_validation, stale_raw, stale_sampler,
    )
    stale_thermal = stale_health["thermal"]
    stale_thermal.update({
        **stale_raw_summary,
        "validation_attempt_errors_total": stale_validation_summary[
            "attempt_errors_total"
        ],
        "validation_failures": stale_validation_summary["failures"],
        "validation_jsonl_ascii": stale_validation.decode("ascii"),
        "validation_jsonl_bytes": len(stale_validation),
        "validation_jsonl_sha256": sha256_bytes(stale_validation),
        "validation_sample_count": stale_validation_summary[
            "sample_count"
        ],
        "window_csv_ascii": stale_raw.decode("ascii"),
        "window_csv_bytes": len(stale_raw),
        "window_csv_sha256": sha256_bytes(stale_raw),
        "window_end_monotonic_ns": stale_sampler[
            "window_end_monotonic_ns"
        ],
        "window_start_monotonic_ns": stale_sampler[
            "window_start_monotonic_ns"
        ],
    })
    stale_health["sampler_receipt_sha256"] = sha256_bytes(
        canonical_bytes(stale_sampler)
    )
    stale_marker = (
        "sampler endpoint: sampler-monitor-invalid:"
        "validation-heartbeat-stale"
    )
    stale_health = finalize_health_receipt(
        stale_health, [stale_marker],
        thermal_policy_violation_messages(stale_thermal),
    )
    validate_health_receipt(stale_health, binary, require_complete=False)
    coverage_cutoff_health = json.loads(
        canonical_bytes(stale_health).decode("ascii")
    )
    coverage_cutoff_marker = (
        "sampler endpoint: sampler-monitor-invalid:"
        "coverage-recovery-cutoff"
    )
    coverage_cutoff_health = finalize_health_receipt(
        coverage_cutoff_health, [coverage_cutoff_marker],
        coverage_cutoff_health["violations"],
    )
    validate_health_receipt(
        coverage_cutoff_health, binary, require_complete=False
    )
    stale_without_marker = json.loads(
        canonical_bytes(stale_health).decode("ascii")
    )
    stale_without_marker = finalize_health_receipt(
        stale_without_marker, [],
        thermal_policy_violation_messages(stale_without_marker["thermal"]),
    )
    try:
        validate_health_receipt(
            stale_without_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("unlicensed short thermal window was accepted")
    appended_sampler = json.loads(
        canonical_bytes(health["sampler"]).decode("ascii")
    )
    appended_sampler["csv_stat"]["size"] += 4096
    appended_sampler["csv_stat"]["mtime_ns"] += 1
    validate_sampler_growth_binding(
        health["sampler_admission"], appended_sampler, "synthetic append"
    )
    changed_sampler_mode = json.loads(
        canonical_bytes(appended_sampler).decode("ascii")
    )
    changed_sampler_mode["csv_stat"]["mode"] = 0o644
    try:
        validate_sampler_growth_binding(
            health["sampler_admission"], changed_sampler_mode,
            "synthetic mode mutation",
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("sampler CSV mode mutation was accepted")
    invalid_health = json.loads(canonical_bytes(health).decode("ascii"))
    header = ",".join(THERMAL_HEADER) + "\n"
    bad_rows = []
    for utc, monotonic_value in (
        ("2026-01-01T00:00:01Z", "1.000000"),
        ("2026-01-01T00:00:02Z", "2.000000"),
    ):
        row = [
            utc, monotonic_value, "50.0", "3000.0", "-1.0",
            *("40.0" for _ in range(8)), "0", "0.1", "0.1", "0.1",
            "1", "0",
        ]
        bad_rows.append(",".join(row) + "\n")
    bad_ascii = header + "".join(bad_rows)
    bad_bytes = bad_ascii.encode("ascii")
    bad_summary = summarize_thermal_window_bytes(
        bad_bytes, 1_000_000_000, 2_000_000_000
    )
    invalid_health["thermal"].update({
        **bad_summary,
        "window_csv_ascii": bad_ascii,
        "window_csv_bytes": len(bad_bytes),
        "window_csv_sha256": sha256_bytes(bad_bytes),
    })
    invalid_thermal_violations = thermal_policy_violation_messages(
        invalid_health["thermal"]
    )
    invalid_health = finalize_health_receipt(
        invalid_health, [], invalid_thermal_violations
    )
    validate_health_receipt(invalid_health, binary, require_complete=False)
    try:
        validate_health_receipt(invalid_health, binary, require_complete=True)
    except ValidationError:
        pass
    else:
        raise AssertionError("invalid thermal/DIMM/EDAC receipt licensed a pass")
    if invalid_health["thermal"]["cpu_tctl_max_millic"] != -1000:
        raise AssertionError("negative invalid temperature evidence was lost")

    missing_thermal_marker = json.loads(
        canonical_bytes(invalid_health).decode("ascii")
    )
    missing_thermal_marker = finalize_health_receipt(
        missing_thermal_marker, [], invalid_thermal_violations[1:]
    )
    try:
        validate_health_receipt(
            missing_thermal_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("missing deterministic thermal marker was accepted")
    interrupted_thermal_merge = json.loads(
        canonical_bytes(invalid_health).decode("ascii")
    )
    interrupted_thermal_merge = finalize_health_receipt(
        interrupted_thermal_merge,
        ["sampler endpoint: injected marker merge interruption"],
        invalid_thermal_violations[:1],
    )
    try:
        validate_health_receipt(
            interrupted_thermal_merge, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "transactionally impossible thermal marker prefix was accepted"
        )
    partial_thermal_collection = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    partial_thermal_collection["thermal"] = {}
    partial_thermal_collection = finalize_health_receipt(
        partial_thermal_collection,
        ["thermal collection: injected pre-publication failure"], [],
    )
    validate_health_receipt(
        partial_thermal_collection, binary, require_complete=False
    )
    missing_thermal_collection_marker = json.loads(
        canonical_bytes(partial_thermal_collection).decode("ascii")
    )
    missing_thermal_collection_marker = finalize_health_receipt(
        missing_thermal_collection_marker, [], []
    )
    try:
        validate_health_receipt(
            missing_thermal_collection_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "empty thermal receipt without its collection marker was accepted"
        )
    false_thermal_collection_marker = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    false_thermal_collection_marker = finalize_health_receipt(
        false_thermal_collection_marker,
        ["thermal collection: false successful-collection marker"], [],
    )
    try:
        validate_health_receipt(
            false_thermal_collection_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "nonempty thermal receipt carried a false collection marker"
        )
    false_empty_thermal_marker = json.loads(
        canonical_bytes(partial_thermal_collection).decode("ascii")
    )
    false_empty_thermal_marker = finalize_health_receipt(
        false_empty_thermal_marker,
        false_empty_thermal_marker["collection_failures"],
        ["CPU thermal cap violated or unavailable"],
    )
    try:
        validate_health_receipt(
            false_empty_thermal_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("empty thermal receipt carried a false marker")
    half_bound_sampler = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    half_bound_sampler["sampler_receipt_sha256"] = None
    half_bound_sampler["receipt_sha256"] = None
    half_bound_sampler["receipt_sha256"] = sha256_bytes(
        canonical_bytes(half_bound_sampler)
    )
    try:
        validate_health_receipt(
            half_bound_sampler, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("half-bound sampler identity was accepted")
    false_thermal_marker = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    false_thermal_marker = finalize_health_receipt(
        false_thermal_marker, [],
        ["CPU thermal cap violated or unavailable"],
    )
    try:
        validate_health_receipt(
            false_thermal_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("false deterministic thermal marker was accepted")

    run_cap_marker = "CPU56 run non-idle tick cap exceeded"
    run_cap_health = json.loads(canonical_bytes(health).decode("ascii"))
    run_sibling = run_cap_health["sibling_ticks"][0]
    run_sibling["end"]["tick_fields"]["user"] = (
        run_sibling["start"]["tick_fields"]["user"]
        + SIBLING_NON_IDLE_TICK_CAP + 1
    )
    run_sibling["end"]["non_idle_ticks"] = sum(
        run_sibling["end"]["tick_fields"][name]
        for name in ("user", "nice", "system", "irq", "softirq", "steal")
    )
    run_sibling["delta_non_idle_ticks"] = (
        run_sibling["end"]["non_idle_ticks"]
        - run_sibling["start"]["non_idle_ticks"]
    )
    run_cap_health = finalize_health_receipt(
        run_cap_health, [], [run_cap_marker]
    )
    validate_health_receipt(run_cap_health, binary, require_complete=False)
    missing_run_cap_marker = json.loads(
        canonical_bytes(run_cap_health).decode("ascii")
    )
    missing_run_cap_marker = finalize_health_receipt(
        missing_run_cap_marker, [], []
    )
    try:
        validate_health_receipt(
            missing_run_cap_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("missing run sibling cap marker was accepted")
    interrupted_run_marker = json.loads(
        canonical_bytes(run_cap_health).decode("ascii")
    )
    interrupted_run_marker = finalize_health_receipt(
        interrupted_run_marker,
        ["run sibling ticks: injected marker interruption"], [],
    )
    validate_health_receipt(
        interrupted_run_marker, binary, require_complete=False
    )
    interrupted_run_after_marker = json.loads(
        canonical_bytes(run_cap_health).decode("ascii")
    )
    interrupted_run_after_marker = finalize_health_receipt(
        interrupted_run_after_marker,
        ["run sibling ticks: injected post-marker interruption"],
        [run_cap_marker],
    )
    validate_health_receipt(
        interrupted_run_after_marker, binary, require_complete=False
    )
    false_run_cap_marker = json.loads(canonical_bytes(health).decode("ascii"))
    false_run_cap_marker = finalize_health_receipt(
        false_run_cap_marker, [], [run_cap_marker]
    )
    try:
        validate_health_receipt(
            false_run_cap_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("false run sibling cap marker was accepted")
    false_interrupted_run_marker = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    false_interrupted_run_marker = finalize_health_receipt(
        false_interrupted_run_marker,
        ["run sibling ticks: injected low-delta interruption"],
        [run_cap_marker],
    )
    try:
        validate_health_receipt(
            false_interrupted_run_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "interrupted low-delta run carried a false cap marker"
        )
    missing_run_false_marker = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    missing_run_false_marker["sibling_ticks"] = []
    missing_run_false_marker = finalize_health_receipt(
        missing_run_false_marker,
        ["run sibling ticks unavailable"], [run_cap_marker],
    )
    try:
        validate_health_receipt(
            missing_run_false_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "missing run sibling carried a false cap marker"
        )

    affinity_drift_marker = "controller CPU121 affinity changed"
    for observed_affinity in ([120], [120, 121]):
        drift_health = json.loads(canonical_bytes(health).decode("ascii"))
        drift_health["controller_singleton_affinity_end"] = observed_affinity
        drift_health = finalize_health_receipt(
            drift_health, [], [affinity_drift_marker]
        )
        validate_health_receipt(drift_health, binary, require_complete=False)
        expect_drift = lambda drift_health=drift_health: (
            validate_health_receipt(
                drift_health, binary, require_complete=True
            )
        )
        try:
            expect_drift()
        except ValidationError:
            pass
        else:
            raise AssertionError(
                "controller affinity drift licensed a positive outcome"
            )

    unmarked_drift = json.loads(canonical_bytes(health).decode("ascii"))
    unmarked_drift["controller_singleton_affinity_end"] = [120, 121]
    unmarked_drift["receipt_sha256"] = None
    unmarked_drift["receipt_sha256"] = sha256_bytes(
        canonical_bytes(unmarked_drift)
    )
    try:
        validate_health_receipt(
            unmarked_drift, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("unmarked controller affinity drift was accepted")
    for markers in ([], [affinity_drift_marker]):
        interrupted_drift = json.loads(
            canonical_bytes(health).decode("ascii")
        )
        interrupted_drift["controller_singleton_affinity_end"] = [120]
        interrupted_drift = finalize_health_receipt(
            interrupted_drift,
            ["controller affinity end: injected marker interruption"],
            markers,
        )
        validate_health_receipt(
            interrupted_drift, binary, require_complete=False
        )

    false_drift_marker = finalize_health_receipt(
        json.loads(canonical_bytes(health).decode("ascii")), [],
        [affinity_drift_marker],
    )
    try:
        validate_health_receipt(
            false_drift_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError("false controller affinity drift marker was accepted")

    missing_affinity_false_marker = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    missing_affinity_false_marker[
        "controller_singleton_affinity_end"
    ] = []
    missing_affinity_false_marker = finalize_health_receipt(
        missing_affinity_false_marker, [], [affinity_drift_marker]
    )
    try:
        validate_health_receipt(
            missing_affinity_false_marker, binary, require_complete=False
        )
    except ValidationError:
        pass
    else:
        raise AssertionError(
            "missing controller affinity carried a false drift marker"
        )
    missing_affinity_partial = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    missing_affinity_partial[
        "controller_singleton_affinity_end"
    ] = []
    missing_affinity_partial = finalize_health_receipt(
        missing_affinity_partial,
        ["controller affinity end: injected collection failure"], [],
    )
    validate_health_receipt(
        missing_affinity_partial, binary, require_complete=False
    )

    partial_binary = empty_binary_receipt(Path("/tmp/partial-worker"), "e" * 64)
    partial_health = empty_health_receipt(argparse.Namespace())
    partial_health = finalize_health_receipt(
        partial_health, ["binary preflight failed"], []
    )
    validate_binary_receipt(partial_binary, require_complete=False)
    validate_health_receipt(partial_health, partial_binary, require_complete=False)

    def expect_validation_error(label: str, operation: Any) -> None:
        try:
            operation()
        except ValidationError:
            return
        raise AssertionError(f"{label} mutation was accepted")

    def expect_exact_validation_error(
        label: str, expected: str, operation: Any,
    ) -> None:
        try:
            operation()
        except ValidationError as exc:
            if str(exc) == expected:
                return
            raise AssertionError(
                f"{label} raised the wrong validation error: {exc}"
            ) from exc
        raise AssertionError(f"{label} mutation was accepted")

    for field, value in (
        ("timed_out", True),
        ("output_overflow", True),
        ("capture_error", "injected capture diagnostic"),
    ):
        unstarted_capture_evidence = json.loads(
            canonical_bytes(partial_binary).decode("ascii")
        )
        unstarted_capture_evidence[field] = value
        expect_validation_error(
            "unstarted binary " + field,
            lambda unstarted_capture_evidence=unstarted_capture_evidence:
            validate_binary_receipt(
                unstarted_capture_evidence, require_complete=False
            ),
        )

    incomplete_execution_binary = json.loads(
        canonical_bytes(binary).decode("ascii")
    )
    incomplete_execution_binary["execution_finished_monotonic_ns"] = None
    expect_validation_error(
        "binary execution interval pair",
        lambda: validate_binary_receipt(
            incomplete_execution_binary, require_complete=False
        ),
    )
    for label, field in (
        ("started binary missing PID", "child_pid"),
        ("started binary missing start timestamp", "child_start_monotonic_ns"),
    ):
        missing_started_field = json.loads(
            canonical_bytes(binary).decode("ascii")
        )
        missing_started_field[field] = None
        expect_validation_error(
            label,
            lambda missing_started_field=missing_started_field:
            validate_binary_receipt(
                missing_started_field, require_complete=False
            ),
        )
    affinity_without_ticks = json.loads(
        canonical_bytes(binary).decode("ascii")
    )
    affinity_without_ticks["child_process_start_ticks"] = None
    expect_validation_error(
        "binary affinity without process-start identity",
        lambda: validate_binary_receipt(
            affinity_without_ticks, require_complete=False
        ),
    )
    ticks_without_affinity = json.loads(
        canonical_bytes(binary).decode("ascii")
    )
    ticks_without_affinity["child_affinity_after_spawn"] = []
    validate_binary_receipt(ticks_without_affinity, require_complete=False)
    for label, mutate in (
        (
            "started binary missing pre-exec stat",
            lambda receipt: receipt.__setitem__("stat_before", {}),
        ),
        (
            "started binary wrong pre-exec mode",
            lambda receipt: receipt["stat_before"].__setitem__("mode", 0o444),
        ),
        (
            "started binary wrong pre-exec link count",
            lambda receipt: receipt["stat_before"].__setitem__("nlink", 2),
        ),
        (
            "started binary wrong pre-exec hash",
            lambda receipt: receipt.__setitem__("sha256_before", "0" * 64),
        ),
        (
            "started binary oversized pre-exec image",
            lambda receipt: receipt["stat_before"].__setitem__(
                "size", MAX_BINARY_BYTES + 1
            ),
        ),
    ):
        changed_preexec = json.loads(canonical_bytes(binary).decode("ascii"))
        mutate(changed_preexec)
        expect_validation_error(
            label,
            lambda changed_preexec=changed_preexec: validate_binary_receipt(
                changed_preexec, require_complete=False
            ),
        )

    late_sampler_health = json.loads(canonical_bytes(health).decode("ascii"))
    late_sampler_health["sampler"]["window_end_monotonic_ns"] = (
        binary["execution_finished_monotonic_ns"] + 1
    )
    late_sampler_health["sampler_receipt_sha256"] = sha256_bytes(
        canonical_bytes(late_sampler_health["sampler"])
    )
    late_sampler_health["receipt_sha256"] = None
    late_sampler_health["receipt_sha256"] = sha256_bytes(
        canonical_bytes(late_sampler_health)
    )
    expect_validation_error(
        "health sampler after execution finish",
        lambda: validate_health_receipt(
            late_sampler_health, binary, require_complete=False
        ),
    )
    late_sampler_start_health = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    late_sampler_start_health["sampler"]["window_start_monotonic_ns"] = (
        binary["execution_started_monotonic_ns"] + 1
    )
    late_sampler_start_health["sampler_receipt_sha256"] = sha256_bytes(
        canonical_bytes(late_sampler_start_health["sampler"])
    )
    late_sampler_start_health["receipt_sha256"] = None
    late_sampler_start_health["receipt_sha256"] = sha256_bytes(
        canonical_bytes(late_sampler_start_health)
    )
    expect_validation_error(
        "health sampler start after execution start",
        lambda: validate_health_receipt(
            late_sampler_start_health, binary, require_complete=False
        ),
    )
    for label, endpoint, timestamp in (
        (
            "run sibling start read before execution",
            "start", binary["execution_started_monotonic_ns"] - 1,
        ),
        (
            "run sibling end read after execution",
            "end", binary["execution_finished_monotonic_ns"] + 1,
        ),
    ):
        escaped_sibling_health = json.loads(
            canonical_bytes(health).decode("ascii")
        )
        escaped_sibling_health["sibling_ticks"][0][endpoint][
            "read_monotonic_ns"
        ] = timestamp
        escaped_sibling_health["receipt_sha256"] = None
        escaped_sibling_health["receipt_sha256"] = sha256_bytes(
            canonical_bytes(escaped_sibling_health)
        )
        expect_validation_error(
            label,
            lambda escaped_sibling_health=escaped_sibling_health:
            validate_health_receipt(
                escaped_sibling_health, binary, require_complete=False
            ),
        )

    terminal_health = synthetic_terminal_health_fixture(health, binary)
    validate_health_receipt(
        synthetic_terminal_health_fixture(
            health, binary, terminal_signal="SIGTERM"
        ),
        binary, require_complete=False,
    )

    def rebind_synthetic_terminal_producer(
        terminal: Dict[str, Any], mutate: Any,
    ) -> None:
        producer = parse_canonical_json_line(
            terminal["producer_receipt_ascii"].encode("ascii"),
            "synthetic rebound sampler receipt",
        )
        mutate(producer)
        preimage = dict(producer)
        del preimage["self_sha256_excluding_field"]
        producer["self_sha256_excluding_field"] = sha256_bytes(
            canonical_bytes(preimage) + b"\n"
        )
        producer_ascii = canonical_bytes(producer) + b"\n"
        terminal["producer_receipt_ascii"] = producer_ascii.decode("ascii")
        terminal["producer_receipt_binding"]["size"] = len(producer_ascii)
        terminal["producer_receipt_binding"]["sha256"] = sha256_bytes(
            producer_ascii
        )

    terminal_timestamp_failure_health = json.loads(
        canonical_bytes(terminal_health).decode("ascii")
    )
    timestamp_failure_terminal = terminal_timestamp_failure_health[
        "sampler_terminal"
    ]

    def make_terminal_timestamp_failure(producer: Dict[str, Any]) -> None:
        producer["errors"] = [{
            "message": "synthetic final timestamp failure",
            "phase": "build_terminal_timestamps", "type": "OSError",
        }]
        producer["exit_code"] = 5
        producer["finished_monotonic_ns"] = None
        producer["finished_utc"] = None
        producer["outcome"] = "sampler_error"

    rebind_synthetic_terminal_producer(
        timestamp_failure_terminal, make_terminal_timestamp_failure
    )
    terminal_timestamp_failure_health[
        "sampler_terminal_receipt_sha256"
    ] = sha256_bytes(canonical_bytes(timestamp_failure_terminal))
    terminal_timestamp_failure_health = finalize_health_receipt(
        terminal_timestamp_failure_health,
        terminal_timestamp_failure_health["collection_failures"],
        terminal_timestamp_failure_health["violations"],
    )
    validate_health_receipt(
        terminal_timestamp_failure_health, binary, require_complete=False
    )
    timestamp_failure_without_phase = json.loads(
        canonical_bytes(timestamp_failure_terminal).decode("ascii")
    )

    def replace_timestamp_failure_phase(producer: Dict[str, Any]) -> None:
        producer["errors"][0]["phase"] = "sampling"

    rebind_synthetic_terminal_producer(
        timestamp_failure_without_phase, replace_timestamp_failure_phase
    )
    expect_validation_error(
        "missing finish timestamps without exact producer error",
        lambda: validate_sampler_terminal_receipt(
            timestamp_failure_without_phase,
            terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )
    timestamp_failure_half_pair = json.loads(
        canonical_bytes(timestamp_failure_terminal).decode("ascii")
    )
    rebind_synthetic_terminal_producer(
        timestamp_failure_half_pair,
        lambda producer: producer.__setitem__(
            "finished_utc", "2026-01-01T00:00:03Z"
        ),
    )
    expect_validation_error(
        "half-missing producer finish timestamp pair",
        lambda: validate_sampler_terminal_receipt(
            timestamp_failure_half_pair,
            terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )

    sampler_error_suffix_health = json.loads(
        canonical_bytes(terminal_health).decode("ascii")
    )
    suffix_terminal = sampler_error_suffix_health["sampler_terminal"]
    paired_raw = health["thermal"]["window_csv_ascii"].encode("ascii")
    paired_validation = health["thermal"][
        "validation_jsonl_ascii"
    ].encode("ascii")
    unpaired_row = parse_csv_physical_line(
        paired_raw.splitlines(keepends=True)[-1],
        "synthetic sampler-error unpaired raw row",
    )
    unpaired_row[0] = "2026-01-01T00:00:02.060000Z"
    unpaired_row[1] = "2.060000"
    unpaired_raw_line = (",".join(unpaired_row) + "\n").encode("ascii")
    next_validation = parse_canonical_json_line(
        paired_validation.splitlines(keepends=True)[-1],
        "synthetic sampler-error partial validation source",
    )
    next_validation["sample_index"] += 1
    next_validation["monotonic_s"] = 2.06
    sampler_error_summary, _ = replay_terminal_validation(
        paired_raw + unpaired_raw_line,
        paired_validation + canonical_bytes(next_validation) + b"\n",
        terminal_health["sampler_admission"],
    )
    partial_validation = canonical_bytes(next_validation)
    partial_validation = partial_validation[:len(partial_validation) // 2]
    suffix_raw = paired_raw + unpaired_raw_line
    suffix_validation = paired_validation + partial_validation
    suffix_terminal["raw_csv_ascii"] = suffix_raw.decode("ascii")
    suffix_terminal["raw_csv_binding"]["size"] = len(suffix_raw)
    suffix_terminal["raw_csv_binding"]["sha256"] = sha256_bytes(suffix_raw)
    suffix_terminal["validation_jsonl_ascii"] = suffix_validation.decode(
        "ascii"
    )
    suffix_terminal["validation_jsonl_binding"]["size"] = len(
        suffix_validation
    )
    suffix_terminal["validation_jsonl_binding"]["sha256"] = sha256_bytes(
        suffix_validation
    )
    suffix_terminal["window_csv_ascii"] = paired_raw.decode("ascii")
    suffix_terminal["window_validation_jsonl_ascii"] = (
        paired_validation.decode("ascii")
    )
    suffix_terminal["stream_suffixes"] = split_terminal_streams(
        suffix_raw, suffix_validation, allow_unpaired=True
    )["suffix_receipt"]

    def make_sampler_error_suffix(producer: Dict[str, Any]) -> None:
        producer["errors"] = [{
            "message": "synthetic validation write stalled",
            "phase": "sampling", "type": "OSError",
        }]
        producer["exit_code"] = 5
        producer["outcome"] = "sampler_error"
        producer["summary"] = sampler_error_summary
        producer["raw_csv"]["binding"] = dict(
            suffix_terminal["raw_csv_binding"]
        )
        producer["validation_jsonl"]["binding"] = dict(
            suffix_terminal["validation_jsonl_binding"]
        )

    rebind_synthetic_terminal_producer(
        suffix_terminal, make_sampler_error_suffix
    )
    sampler_error_suffix_health["sampler_terminal_receipt_sha256"] = (
        sha256_bytes(canonical_bytes(suffix_terminal))
    )
    sampler_error_suffix_health = finalize_health_receipt(
        sampler_error_suffix_health,
        ["sampler endpoint: sampler-terminal-intent"],
        sampler_error_suffix_health["violations"],
    )
    validate_health_receipt(
        sampler_error_suffix_health, binary, require_complete=False
    )

    validation_only_partial = json.loads(
        canonical_bytes(suffix_terminal).decode("ascii")
    )
    validation_only_partial["raw_csv_ascii"] = paired_raw.decode("ascii")
    validation_only_partial["raw_csv_binding"]["size"] = len(paired_raw)
    validation_only_partial["raw_csv_binding"]["sha256"] = sha256_bytes(
        paired_raw
    )
    validation_only_partial["stream_suffixes"][
        "raw_unpaired_suffix_bytes"
    ] = 0
    validation_only_partial["stream_suffixes"][
        "raw_unpaired_suffix_sha256"
    ] = sha256_bytes(b"")

    def bind_validation_only_partial(producer: Dict[str, Any]) -> None:
        producer["raw_csv"]["binding"] = dict(
            validation_only_partial["raw_csv_binding"]
        )

    rebind_synthetic_terminal_producer(
        validation_only_partial, bind_validation_only_partial
    )
    expect_validation_error(
        "sampler validation partial without preceding raw row",
        lambda: validate_sampler_terminal_receipt(
            validation_only_partial, terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )

    future_unpaired_raw = json.loads(
        canonical_bytes(suffix_terminal).decode("ascii")
    )
    future_raw_lines = future_unpaired_raw["raw_csv_ascii"].encode(
        "ascii"
    ).splitlines(keepends=True)
    future_row = parse_csv_physical_line(
        future_raw_lines[-1], "synthetic future unpaired raw row"
    )
    future_row[1] = "9999.000000"
    future_raw_lines[-1] = (",".join(future_row) + "\n").encode("ascii")
    future_raw = b"".join(future_raw_lines)
    future_unpaired_raw["raw_csv_ascii"] = future_raw.decode("ascii")
    future_unpaired_raw["raw_csv_binding"]["size"] = len(future_raw)
    future_unpaired_raw["raw_csv_binding"]["sha256"] = sha256_bytes(
        future_raw
    )
    future_unpaired_raw["stream_suffixes"] = split_terminal_streams(
        future_raw, suffix_validation, allow_unpaired=True
    )["suffix_receipt"]

    def bind_future_unpaired_raw(producer: Dict[str, Any]) -> None:
        producer["raw_csv"]["binding"] = dict(
            future_unpaired_raw["raw_csv_binding"]
        )

    rebind_synthetic_terminal_producer(
        future_unpaired_raw, bind_future_unpaired_raw
    )
    expect_validation_error(
        "sampler unpaired raw row after producer finish",
        lambda: validate_sampler_terminal_receipt(
            future_unpaired_raw, terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )

    pid_unlink_after_removal_health = json.loads(
        canonical_bytes(terminal_health).decode("ascii")
    )

    def make_unlink_after_removal_error(producer: Dict[str, Any]) -> None:
        producer["errors"] = [{
            "message": "synthetic directory fsync after unlink failed",
            "phase": "unlink_pid_file", "type": "OSError",
        }]
        producer["exit_code"] = 5
        producer["outcome"] = "sampler_error"
        producer["pid_file"]["removed"] = False

    rebind_synthetic_terminal_producer(
        pid_unlink_after_removal_health["sampler_terminal"],
        make_unlink_after_removal_error,
    )
    pid_unlink_after_removal_health["sampler_terminal_receipt_sha256"] = (
        sha256_bytes(
            canonical_bytes(
                pid_unlink_after_removal_health["sampler_terminal"]
            )
        )
    )
    pid_unlink_after_removal_health = finalize_health_receipt(
        pid_unlink_after_removal_health,
        pid_unlink_after_removal_health["collection_failures"],
        pid_unlink_after_removal_health["violations"],
    )
    validate_health_receipt(
        pid_unlink_after_removal_health, binary, require_complete=False
    )
    seal_recovery_healths: List[Tuple[str, Dict[str, Any]]] = []
    for artifact_name, error_phase in (
        ("raw_csv", "seal_raw_csv"),
        ("validation_jsonl", "seal_validation_jsonl"),
    ):
        recovered_health = json.loads(
            canonical_bytes(terminal_health).decode("ascii")
        )

        def make_recoverable_seal_error(
            producer: Dict[str, Any], *, artifact_name: str = artifact_name,
            error_phase: str = error_phase,
        ) -> None:
            producer["errors"] = [{
                "message": "synthetic post-fchmod fsync failed",
                "phase": error_phase, "type": "OSError",
            }]
            producer["exit_code"] = 5
            producer["outcome"] = "sampler_error"
            producer[artifact_name]["binding"] = None
            if artifact_name == "raw_csv":
                producer["raw_samples_preserved"] = False

        rebind_synthetic_terminal_producer(
            recovered_health["sampler_terminal"],
            make_recoverable_seal_error,
        )
        recovered_health["sampler_terminal_receipt_sha256"] = sha256_bytes(
            canonical_bytes(recovered_health["sampler_terminal"])
        )
        recovered_health = finalize_health_receipt(
            recovered_health, recovered_health["collection_failures"],
            recovered_health["violations"],
        )
        validate_health_receipt(
            recovered_health, binary, require_complete=False
        )
        seal_recovery_healths.append((error_phase, recovered_health))

    unlicensed_missing_binding = json.loads(
        canonical_bytes(terminal_health["sampler_terminal"]).decode("ascii")
    )

    def make_unlicensed_missing_binding(producer: Dict[str, Any]) -> None:
        producer["errors"] = [{
            "message": "synthetic close failure",
            "phase": "close_raw_csv", "type": "OSError",
        }]
        producer["exit_code"] = 5
        producer["outcome"] = "sampler_error"
        producer["raw_csv"]["binding"] = None
        producer["raw_samples_preserved"] = False

    rebind_synthetic_terminal_producer(
        unlicensed_missing_binding, make_unlicensed_missing_binding
    )
    expect_validation_error(
        "sampler missing raw binding without seal error",
        lambda: validate_sampler_terminal_receipt(
            unlicensed_missing_binding, terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )
    unlicensed_false_removed = json.loads(
        canonical_bytes(terminal_health["sampler_terminal"]).decode("ascii")
    )
    rebind_synthetic_terminal_producer(
        unlicensed_false_removed,
        lambda producer: producer["pid_file"].__setitem__("removed", False),
    )
    expect_validation_error(
        "sampler false PID removed flag without unlink error",
        lambda: validate_sampler_terminal_receipt(
            unlicensed_false_removed, terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )
    terminal_started_after_window = json.loads(
        canonical_bytes(terminal_health["sampler_terminal"]).decode("ascii")
    )
    late_start_producer = parse_canonical_json_line(
        terminal_started_after_window["producer_receipt_ascii"].encode(
            "ascii"
        ),
        "synthetic late-start sampler receipt",
    )
    late_start_producer["started_monotonic_ns"] = (
        terminal_started_after_window["coverage"][
            "window_start_monotonic_ns"
        ] + 1
    )
    late_start_preimage = dict(late_start_producer)
    del late_start_preimage["self_sha256_excluding_field"]
    late_start_producer["self_sha256_excluding_field"] = sha256_bytes(
        canonical_bytes(late_start_preimage) + b"\n"
    )
    late_start_ascii = canonical_bytes(late_start_producer) + b"\n"
    terminal_started_after_window["producer_receipt_ascii"] = (
        late_start_ascii.decode("ascii")
    )
    terminal_started_after_window["producer_receipt_binding"]["size"] = (
        len(late_start_ascii)
    )
    terminal_started_after_window["producer_receipt_binding"]["sha256"] = (
        sha256_bytes(late_start_ascii)
    )
    expect_validation_error(
        "sampler terminal producer start after retained window",
        lambda: validate_sampler_terminal_receipt(
            terminal_started_after_window,
            terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )
    terminal_exit_before_finish = json.loads(
        canonical_bytes(terminal_health["sampler_terminal"]).decode("ascii")
    )
    terminal_producer = parse_canonical_json_line(
        terminal_exit_before_finish["producer_receipt_ascii"].encode(
            "ascii"
        ),
        "synthetic early-exit sampler receipt",
    )
    terminal_exit_before_finish["process_exit_observed_monotonic_ns"] = (
        terminal_producer["finished_monotonic_ns"] - 1
    )
    expect_validation_error(
        "sampler terminal process exit before producer finish",
        lambda: validate_sampler_terminal_receipt(
            terminal_exit_before_finish,
            terminal_health["sampler_admission"],
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
        ),
    )
    terminal_exit_after_execution = json.loads(
        canonical_bytes(terminal_health).decode("ascii")
    )
    terminal_exit_after_execution["sampler_terminal"][
        "process_exit_observed_monotonic_ns"
    ] = binary["execution_finished_monotonic_ns"] + 1
    terminal_exit_after_execution["sampler_terminal_receipt_sha256"] = (
        sha256_bytes(
            canonical_bytes(terminal_exit_after_execution["sampler_terminal"])
        )
    )
    terminal_exit_after_execution = finalize_health_receipt(
        terminal_exit_after_execution,
        terminal_exit_after_execution["collection_failures"],
        terminal_exit_after_execution["violations"],
    )
    expect_validation_error(
        "sampler terminal exit after binary execution finish",
        lambda: validate_health_receipt(
            terminal_exit_after_execution, binary, require_complete=False
        ),
    )
    stale_recovery_failure_health = json.loads(
        canonical_bytes(health).decode("ascii")
    )
    stale_recovery_failure_health["sampler"] = {}
    stale_recovery_failure_health["sampler_receipt_sha256"] = None
    stale_recovery_failure_health["sampler_terminal"] = {}
    stale_recovery_failure_health["sampler_terminal_receipt_sha256"] = None
    stale_recovery_failure_health["thermal"] = {}
    stale_recovery_failure_health = finalize_health_receipt(
        stale_recovery_failure_health,
        [
            stale_marker,
            "sampler stale evidence: injected fixed-point failure",
        ],
        [],
    )
    validate_health_receipt(
        stale_recovery_failure_health, binary, require_complete=False
    )

    final_fixture = synthetic_final_bundle_fixture(health, binary, commit)
    preflight_fixture = synthetic_preflight_seal_fixture(
        final_fixture, health, commit
    )
    validate_preflight_seal_receipt(preflight_fixture)
    preflight_payload = canonical_bytes(preflight_fixture) + b"\n"
    if not canonical_equal(
        parse_preflight_seal_bytes(preflight_payload), preflight_fixture
    ):
        raise AssertionError("canonical preflight receipt parsing differs")
    expect_validation_error(
        "preflight noncanonical JSON",
        lambda: parse_preflight_seal_bytes(b" " + preflight_payload),
    )
    duplicate_preflight_payload = preflight_payload.replace(
        b'{"binary_after":', b'{"binary_after":{},"binary_after":', 1
    )
    expect_validation_error(
        "preflight duplicate JSON key",
        lambda: parse_preflight_seal_bytes(duplicate_preflight_payload),
    )
    preflight_tokens = preflight_fixture[
        "preflight_controller_after"
    ]["argv"][1:]
    parsed_preflight = parse_cli_tokens(preflight_tokens, live_process=False)
    if (
        not parsed_preflight.preflight_seal
        or parsed_preflight.expected_source_commit != commit
        or parsed_preflight.sampler_pid != health["sampler_admission"]["pid"]
    ):
        raise AssertionError("preflight public argument parsing differs")
    reconstructed_values = preflight_run_values(
        parsed_preflight, final_fixture["source"], final_fixture["git"],
        final_fixture["build"], preflight_fixture["binary_after"],
        {
            "sha256": preflight_fixture["preflight_controller_after"][
                "image"
            ]["python_sha256"]
        },
        preflight_fixture["sampler_after"],
    )
    reconstructed_argv = build_run_authority_argv(
        reconstructed_values,
        preflight_fixture["preflight_controller_after"]["pid"],
    )
    if (
        reconstructed_argv != preflight_fixture["run_argv"]
        or run_argv_sha256(reconstructed_argv)
        != preflight_fixture["run_argv_sha256"]
        or len(reconstructed_argv) != 76
    ):
        raise AssertionError("preflight future run authority differs")
    for label, tokens in (
        (
            "preflight reordered locators",
            preflight_tokens[:2] + preflight_tokens[4:6]
            + preflight_tokens[2:4] + preflight_tokens[6:],
        ),
        (
            "preflight equals-form locator",
            ["--preflight-seal", "--binary=" + binary["path"]]
            + preflight_tokens[3:],
        ),
        ("preflight missing locator", preflight_tokens[:-2]),
        ("preflight extra locator", preflight_tokens + ["--cpu", "120"]),
        (
            "preflight noncanonical PID",
            [
                "+" + item if preflight_tokens[index - 1] == "--sampler-pid"
                else item
                for index, item in enumerate(preflight_tokens)
            ],
        ),
    ):
        expect_validation_error(
            label,
            lambda tokens=tokens: parse_cli_tokens(tokens, live_process=False),
        )

    mutated_preflight_run = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    binary_hash_value_index = (
        2 + 2 * RUN_ONCE_OPTION_ORDER.index("--expected-binary-sha256") + 1
    )
    mutated_preflight_run["run_argv"][binary_hash_value_index] = "1" * 64
    mutated_preflight_run["run_argv_sha256"] = run_argv_sha256(
        mutated_preflight_run["run_argv"]
    )
    mutated_preflight_run["receipt_sha256"] = None
    mutated_preflight_run["receipt_sha256"] = sha256_bytes(
        canonical_bytes(mutated_preflight_run)
    )
    expect_exact_validation_error(
        "preflight run/final observation mismatch",
        "preflight run argv differs from final observations",
        lambda: validate_preflight_seal_receipt(mutated_preflight_run),
    )
    surrogate_run_argv = list(reconstructed_argv)
    surrogate_run_argv[3] = "/tmp/\ud800"
    expect_exact_validation_error(
        "preflight future argv surrogate",
        "run authority argv is not a string array",
        lambda: run_argv_sha256(surrogate_run_argv),
    )
    surrogate_preflight = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    for key in ("preflight_controller_before", "preflight_controller_after"):
        surrogate_preflight[key]["argv"][3] = "/tmp/\ud800"
    surrogate_payload = canonical_bytes(surrogate_preflight) + b"\n"
    expect_validation_error(
        "preflight receipt surrogate argv",
        lambda: parse_preflight_seal_bytes(surrogate_payload),
    )
    deeply_nested_preflight = (
        b"[" * 2000 + b"0" + b"]" * 2000 + b"\n"
    )
    expect_validation_error(
        "preflight deeply nested JSON",
        lambda: parse_preflight_seal_bytes(deeply_nested_preflight),
    )
    mutated_preflight_prefix = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    mutated_preflight_prefix["sampler_prefix_binding"][
        "raw_after_sha256"
    ] = "0" * 64
    mutated_preflight_prefix["receipt_sha256"] = None
    mutated_preflight_prefix["receipt_sha256"] = sha256_bytes(
        canonical_bytes(mutated_preflight_prefix)
    )
    expect_validation_error(
        "preflight held sampler prefix mutation",
        lambda: validate_preflight_seal_receipt(mutated_preflight_prefix),
    )
    oversized_preflight_raw = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    oversized_preflight_raw["sampler_prefix_binding"]["final_raw_stat"][
        "size"
    ] = MAX_SAMPLER_CSV_BYTES + 1
    oversized_preflight_raw["receipt_sha256"] = None
    oversized_preflight_raw["receipt_sha256"] = sha256_bytes(
        canonical_bytes(oversized_preflight_raw)
    )
    expect_validation_error(
        "preflight oversized held raw stream",
        lambda: validate_preflight_seal_receipt(oversized_preflight_raw),
    )
    regressed_preflight_raw_mtime = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    regressed_preflight_raw_mtime["sampler_prefix_binding"][
        "final_raw_stat"
    ]["mtime_ns"] = (
        preflight_fixture["sampler_after"]["csv_stat"]["mtime_ns"] - 1
    )
    regressed_preflight_raw_mtime["receipt_sha256"] = None
    regressed_preflight_raw_mtime["receipt_sha256"] = sha256_bytes(
        canonical_bytes(regressed_preflight_raw_mtime)
    )
    expect_validation_error(
        "preflight regressed held raw mtime",
        lambda: validate_preflight_seal_receipt(regressed_preflight_raw_mtime),
    )
    for label, field, replacement in (
        (
            "preflight sampler CSV locator mismatch", "csv_path",
            str(Path(health["sampler_admission"]["csv_path"]).with_name(
                "other.csv"
            )),
        ),
        (
            "preflight sampler executable locator mismatch", "executable_path",
            "/usr/bin/not-python",
        ),
    ):
        mutated_sampler_binding = json.loads(
            canonical_bytes(preflight_fixture).decode("ascii")
        )
        for sampler_name in ("sampler_before", "sampler_after"):
            mutated_sampler_binding[sampler_name][field] = replacement
        prefix_binding = mutated_sampler_binding["sampler_prefix_binding"]
        prefix_binding["sampler_before_sha256"] = sha256_bytes(
            canonical_bytes(mutated_sampler_binding["sampler_before"])
        )
        prefix_binding["sampler_after_sha256"] = sha256_bytes(
            canonical_bytes(mutated_sampler_binding["sampler_after"])
        )
        mutated_sampler_binding["receipt_sha256"] = None
        mutated_sampler_binding["receipt_sha256"] = sha256_bytes(
            canonical_bytes(mutated_sampler_binding)
        )
        expect_validation_error(
            label,
            lambda value=mutated_sampler_binding:
                validate_preflight_seal_receipt(value),
        )
    mutated_sampler_pid = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    replacement_pid = health["sampler_admission"]["pid"] + 1
    replacement_pid_bytes = (str(replacement_pid) + "\n").encode("ascii")
    for sampler_name in ("sampler_before", "sampler_after"):
        sampler_value = mutated_sampler_pid[sampler_name]
        sampler_value["pid"] = replacement_pid
        sampler_value["pid_file"]["bytes"] = len(replacement_pid_bytes)
        sampler_value["pid_file"]["sha256"] = sha256_bytes(
            replacement_pid_bytes
        )
    pid_prefix = mutated_sampler_pid["sampler_prefix_binding"]
    pid_prefix["pid_file_sha256"] = sha256_bytes(replacement_pid_bytes)
    pid_prefix["sampler_before_sha256"] = sha256_bytes(
        canonical_bytes(mutated_sampler_pid["sampler_before"])
    )
    pid_prefix["sampler_after_sha256"] = sha256_bytes(
        canonical_bytes(mutated_sampler_pid["sampler_after"])
    )
    mutated_sampler_pid["receipt_sha256"] = None
    mutated_sampler_pid["receipt_sha256"] = sha256_bytes(
        canonical_bytes(mutated_sampler_pid)
    )
    expect_validation_error(
        "preflight sampler PID mismatch",
        lambda: validate_preflight_seal_receipt(mutated_sampler_pid),
    )
    mutated_preflight_input = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    for field in ("preflight_controller_before", "preflight_controller_after"):
        controller_receipt = mutated_preflight_input[field]
        controller_receipt["argv"][3] = "/tmp/different-worker"
        controller_receipt["argv_sha256"] = sha256_bytes(
            b"".join(
                os.fsencode(item) + b"\0" for item in controller_receipt["argv"]
            )
        )
        controller_receipt["receipt_sha256"] = None
        controller_receipt["receipt_sha256"] = sha256_bytes(
            canonical_bytes(controller_receipt)
        )
    mutated_preflight_input["receipt_sha256"] = None
    mutated_preflight_input["receipt_sha256"] = sha256_bytes(
        canonical_bytes(mutated_preflight_input)
    )
    expect_validation_error(
        "preflight input/future run mismatch",
        lambda: validate_preflight_seal_receipt(mutated_preflight_input),
    )
    mutated_preflight_commit = json.loads(
        canonical_bytes(preflight_fixture).decode("ascii")
    )
    mutated_preflight_commit["expected_source_commit"] = "2" * 40
    mutated_preflight_commit["git_before"]["head"] = "2" * 40
    mutated_preflight_commit["git_after"]["head"] = "2" * 40
    mutated_preflight_commit["receipt_sha256"] = None
    mutated_preflight_commit["receipt_sha256"] = sha256_bytes(
        canonical_bytes(mutated_preflight_commit)
    )
    expect_validation_error(
        "preflight receipt/input commit mismatch",
        lambda: validate_preflight_seal_receipt(mutated_preflight_commit),
    )

    original_preflight_collector = globals()["_collect_preflight_seal"]
    try:
        def injected_preflight_collector(
            _args: argparse.Namespace, _deadline: Optional[float] = None,
        ) -> Dict[str, Any]:
            if not PREFLIGHT_SIDE_EFFECT_GUARD:
                raise AssertionError("preflight side-effect guard was not active")
            expect_exact_validation_error(
                "preflight guarded run-once",
                "preflight side-effect tripwire: run-once dispatch",
                lambda: run_once(argparse.Namespace()),
            )
            expect_exact_validation_error(
                "preflight guarded worker launch",
                "preflight side-effect tripwire: worker launch",
                lambda: run_binary(Path("/tmp/worker"), 120, "0" * 64),
            )
            expect_exact_validation_error(
                "preflight guarded claim reservation",
                "preflight side-effect tripwire: claim reservation",
                lambda: ClaimReservation.reserve(
                    argparse.Namespace(), {}, 0
                ),
            )
            expect_exact_validation_error(
                "preflight guarded output reservation",
                "preflight side-effect tripwire: output reservation",
                lambda: OutputBundle.reserve(
                    Path("/tmp/output"), 120, "0" * 40, "0" * 64,
                    "0" * 64, Path("/tmp"),
                ),
            )
            return preflight_fixture

        globals()["_collect_preflight_seal"] = injected_preflight_collector
        guarded_receipt = collect_preflight_seal(argparse.Namespace())
        if (
            not canonical_equal(guarded_receipt, preflight_fixture)
            or PREFLIGHT_SIDE_EFFECT_GUARD
        ):
            raise AssertionError("preflight side-effect guard lifecycle differs")
    finally:
        globals()["_collect_preflight_seal"] = original_preflight_collector
        globals()["PREFLIGHT_SIDE_EFFECT_GUARD"] = False

    output_read_fd, output_write_fd = os.pipe()
    try:
        output_payload = b'{"preflight":"bounded"}\n'
        write_preflight_payload(
            output_write_fd, output_payload, time.monotonic() + 1.0
        )
        if os.read(output_read_fd, len(output_payload) + 1) != output_payload:
            raise AssertionError("preflight bounded output bytes differ")
    finally:
        os.close(output_read_fd)
        os.close(output_write_fd)

    full_read_fd, full_write_fd = os.pipe()
    try:
        os.set_blocking(full_write_fd, False)
        while True:
            try:
                os.write(full_write_fd, b"x" * 65536)
            except BlockingIOError:
                break
        expect_exact_validation_error(
            "preflight blocked output deadline",
            "preflight seal output reached its global deadline",
            lambda: write_preflight_payload(
                full_write_fd, b"x", time.monotonic() + 0.02
            ),
        )
    finally:
        os.close(full_read_fd)
        os.close(full_write_fd)

    parsed_claim = parse_claim_document(final_fixture["claim_raw"])
    validate_claim_runtime_binding(
        parsed_claim, final_fixture["claim_receipt"],
        final_fixture["claim_document"], final_fixture["claim_raw"],
    )
    parsed_summary = parse_summary_bytes(final_fixture["summary_bytes"])
    parsed_controller = parse_controller_document(
        final_fixture["controller_bytes"]
    )
    oversized_loader_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    oversized_loader_summary["health_module_loader"]["modules"][0][
        "bytes"
    ] = MAX_SOURCE_FILE_BYTES + 1
    oversized_loader_summary["summary_preimage_sha256"] = None
    make_summary_preimage(oversized_loader_summary)
    expect_validation_error(
        "oversized sealed health module receipt",
        lambda: parse_summary_bytes(
            canonical_bytes(oversized_loader_summary) + b"\n"
        ),
    )
    oversized_sampler_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    oversized_sampler_summary["health"]["sampler"]["csv_stat"]["size"] = (
        MAX_SAMPLER_CSV_BYTES + 1
    )
    oversized_sampler_summary["summary_preimage_sha256"] = None
    make_summary_preimage(oversized_sampler_summary)
    expect_validation_error(
        "oversized sampler CSV receipt",
        lambda: parse_summary_bytes(
            canonical_bytes(oversized_sampler_summary) + b"\n"
        ),
    )
    for label, mutate in (
        (
            "oversized controller Python receipt",
            lambda document: document["controller_before"]["python_stat"].__setitem__(
                "size", MAX_BINARY_BYTES + 1
            ),
        ),
        (
            "oversized Git executable receipt",
            lambda document: document["git_before"]["executable"][
                "stat"
            ].__setitem__("size", MAX_BINARY_BYTES + 1),
        ),
        (
            "oversized build provenance receipt",
            lambda document: document["build_before"]["entries"][0].__setitem__(
                "bytes", 64 * 1024 * 1024 + 1
            ),
        ),
    ):
        oversized_controller_document = json.loads(
            canonical_bytes(final_fixture["controller"]).decode("ascii")
        )
        mutate(oversized_controller_document)
        expect_validation_error(
            label,
            lambda oversized_controller_document=oversized_controller_document:
            parse_controller_document(
                controller_document_bytes(oversized_controller_document)
            ),
        )
    parsed_run = validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], parsed_controller, parsed_summary,
        parsed_claim,
    )
    if (
        not parsed_run.run_once
        or parsed_run.expected_sampler_environ_sha256
        != health["sampler"]["environ_sha256"]
    ):
        raise AssertionError("synthetic retained run authority differs")

    def drift_sampler_source_stat(sampler: Dict[str, Any]) -> None:
        sampler["script_inode"] += 1000
        sampler["script_stat"]["inode"] = sampler["script_inode"]
        sampler["script_stat"]["mtime_ns"] += 1

    drifted_sampler_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    drift_sampler_source_stat(drifted_sampler_summary["health"]["sampler"])
    validate_sampler_receipt(drifted_sampler_summary["health"]["sampler"])
    drifted_sampler_controller = json.loads(
        canonical_bytes(parsed_controller).decode("ascii")
    )
    drift_sampler_source_stat(drifted_sampler_controller["sampler_after"])
    validate_sampler_receipt(drifted_sampler_controller["sampler_after"])
    validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], drifted_sampler_controller,
        drifted_sampler_summary, parsed_claim,
    )
    positive_drifted_sampler_controller = dict(drifted_sampler_controller)
    positive_drifted_sampler_controller["outcome"] = "pass"
    expect_validation_error(
        "positive retained sampler source stat drift",
        lambda: validate_recorded_run_authority(
            final_fixture["run_argv_sha256"],
            positive_drifted_sampler_controller, drifted_sampler_summary,
            parsed_claim,
        ),
    )
    for label, field, value in (
        (
            "sampler source path authority", "sampler_script",
            Path("/tmp/different-sampler.py"),
        ),
        (
            "sampler source hash authority", "expected_sampler_script_sha256",
            "0" * 64,
        ),
    ):
        mutated_run = argparse.Namespace(**vars(parsed_run))
        setattr(mutated_run, field, value)
        expect_validation_error(
            label,
            lambda mutated_run=mutated_run: validate_sampler_source_authority(
                mutated_run, parsed_controller["source_manifest_before"],
                Path(__file__).resolve(strict=True).parents[1],
            ),
        )
    expect_validation_error(
        "external run argv authority seal",
        lambda: validate_recorded_run_authority(
            "0" * 64, parsed_controller, parsed_summary, parsed_claim
        ),
    )

    def parsed_binary_evidence_mutation(
        mutate: Any, *, positive: bool = False,
    ) -> Dict[str, Any]:
        document = json.loads(
            canonical_bytes(final_fixture["controller"]).decode("ascii")
        )
        mutate(document["binary"])
        if positive:
            document["outcome"] = "pass"
            document["failure"] = "none"
            document["summary"]["outcome"] = "pass"
        return parse_controller_document(controller_document_bytes(document))

    def make_unstarted(binary_document: Dict[str, Any]) -> None:
        binary_document["process_started"] = False
        binary_document["child_affinity_after_spawn"] = []
        binary_document["child_pid"] = None
        binary_document["child_process_start_ticks"] = None
        binary_document["child_reap_monotonic_ns"] = None
        binary_document["child_start_monotonic_ns"] = None

    sampler_pid_overlap = parsed_binary_evidence_mutation(
        lambda binary_document: binary_document.__setitem__(
            "child_pid", parsed_run.sampler_pid
        )
    )
    expect_validation_error(
        "live sampler/worker PID overlap",
        lambda: validate_recorded_run_authority(
            final_fixture["run_argv_sha256"], sampler_pid_overlap,
            parsed_summary, parsed_claim,
        ),
    )

    dead_sampler_health = json.loads(canonical_bytes(health).decode("ascii"))
    dead_sampler_health["sampler"] = {}
    dead_sampler_health["sampler_receipt_sha256"] = None
    dead_sampler_health["thermal"] = {}
    dead_sampler_health = finalize_health_receipt(
        dead_sampler_health, ["sampler endpoint: process exited before receipt"],
        dead_sampler_health["violations"],
    )
    reused_pid_summary_document = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    reused_pid_summary_document["binary"]["child_pid"] = parsed_run.sampler_pid
    reused_pid_summary_document["health"] = dead_sampler_health
    reused_pid_summary_document["summary_preimage_sha256"] = None
    make_summary_preimage(reused_pid_summary_document)
    reused_pid_summary_bytes = (
        canonical_bytes(reused_pid_summary_document) + b"\n"
    )
    reused_pid_summary = parse_summary_bytes(reused_pid_summary_bytes)
    reused_pid_controller_document = json.loads(
        canonical_bytes(final_fixture["controller"]).decode("ascii")
    )
    reused_pid_controller_document["binary"]["child_pid"] = (
        parsed_run.sampler_pid
    )
    reused_pid_controller_document["health_receipt_sha256"] = (
        dead_sampler_health["receipt_sha256"]
    )
    reused_pid_controller_document["sampler_after"] = {}
    reused_pid_controller_document["summary"]["preimage_sha256"] = (
        reused_pid_summary["summary_preimage_sha256"]
    )
    reused_pid_controller_document["summary"]["sha256"] = sha256_bytes(
        reused_pid_summary_bytes
    )
    reused_pid_controller_document["artifacts"]["summary.json"]["bytes"] = (
        len(reused_pid_summary_bytes)
    )
    reused_pid_controller_document["artifacts"]["summary.json"]["sha256"] = (
        sha256_bytes(reused_pid_summary_bytes)
    )
    reused_pid_controller_document["artifacts"]["thermal.csv"]["bytes"] = 0
    reused_pid_controller_document["artifacts"]["thermal.csv"]["sha256"] = (
        sha256_bytes(b"")
    )
    reused_pid_controller = parse_controller_document(
        controller_document_bytes(reused_pid_controller_document)
    )
    validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], reused_pid_controller,
        reused_pid_summary, parsed_claim,
    )
    validate_controller_bundle_binding(
        reused_pid_controller, reused_pid_summary, reused_pid_summary_bytes,
        final_fixture["claim_receipt"], reused_pid_controller["binary"],
    )
    same_tick_pid_reuse = json.loads(
        canonical_bytes(reused_pid_controller).decode("ascii")
    )
    same_tick_pid_reuse["binary"]["child_process_start_ticks"] = (
        parsed_run.expected_sampler_process_start_ticks
    )
    same_tick_pid_reuse = parse_controller_document(
        controller_document_bytes(same_tick_pid_reuse)
    )
    # Linux start ticks have jiffy resolution; a dead sampler PID can be
    # reused by a worker in the same tick.  Only a retained live sampler
    # attestation makes numeric PID overlap contradictory.
    validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], same_tick_pid_reuse,
        reused_pid_summary, parsed_claim,
    )

    def unstarted_pre_hash_uid_drift(
        binary_document: Dict[str, Any],
    ) -> None:
        make_unstarted(binary_document)
        binary_document["sha256_before"] = None
        binary_document["stat_before"]["uid"] += 1

    unstarted_uid_drift = parsed_binary_evidence_mutation(
        unstarted_pre_hash_uid_drift
    )
    validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], unstarted_uid_drift,
        parsed_summary, parsed_claim,
    )

    started_preexec_uid_drift = parsed_binary_evidence_mutation(
        lambda binary_document: binary_document["stat_before"].__setitem__(
            "uid", binary_document["stat_before"]["uid"] + 1
        )
    )
    expect_validation_error(
        "started pre-exec binary UID authority drift",
        lambda: validate_recorded_run_authority(
            final_fixture["run_argv_sha256"], started_preexec_uid_drift,
            parsed_summary, parsed_claim,
        ),
    )

    for receipt_name in ("stat_after", "path_stat_after"):
        def drift_postexec_uid(
            binary_document: Dict[str, Any],
            receipt_name: str = receipt_name,
        ) -> None:
            receipt = binary_document[receipt_name]
            receipt["uid"] += 1

        invalid_postexec_uid_drift = parsed_binary_evidence_mutation(
            drift_postexec_uid
        )
        validate_recorded_run_authority(
            final_fixture["run_argv_sha256"], invalid_postexec_uid_drift,
            parsed_summary, parsed_claim,
        )
        expect_validation_error(
            "positive post-exec binary UID drift " + receipt_name,
            lambda drift_postexec_uid=drift_postexec_uid:
            parsed_binary_evidence_mutation(
                drift_postexec_uid, positive=True
            ),
        )

    def unlinked_unstarted_preexec(binary_document: Dict[str, Any]) -> None:
        make_unstarted(binary_document)
        binary_document["sha256_before"] = None
        binary_document["stat_before"]["nlink"] = 0
        binary_document["stat_after"]["nlink"] = 0
        binary_document["path_stat_after"] = {}

    parsed_binary_evidence_mutation(unlinked_unstarted_preexec)

    def unlinked_held_postexec(binary_document: Dict[str, Any]) -> None:
        binary_document["stat_after"]["nlink"] = 0

    parsed_binary_evidence_mutation(unlinked_held_postexec)
    expect_validation_error(
        "positive unlinked held binary after exec",
        lambda: parsed_binary_evidence_mutation(
            unlinked_held_postexec, positive=True
        ),
    )

    def oversized_unstarted_preexec(binary_document: Dict[str, Any]) -> None:
        make_unstarted(binary_document)
        binary_document["sha256_before"] = None
        binary_document["stat_before"]["size"] = MAX_BINARY_BYTES + 1
        binary_document["stat_after"]["size"] = MAX_BINARY_BYTES + 1
        binary_document["sha256_after"] = None
        binary_document["path_stat_after"] = {}

    oversized_unstarted_controller = parsed_binary_evidence_mutation(
        oversized_unstarted_preexec
    )
    validate_recorded_run_authority(
        final_fixture["run_argv_sha256"], oversized_unstarted_controller,
        parsed_summary, parsed_claim,
    )
    for receipt_name in ("stat_after", "path_stat_after"):
        def oversize_postexec(
            binary_document: Dict[str, Any],
            receipt_name: str = receipt_name,
        ) -> None:
            binary_document[receipt_name]["size"] = MAX_BINARY_BYTES + 1
            if receipt_name == "stat_after":
                binary_document["sha256_after"] = None
                binary_document["path_stat_after"] = {}

        oversized_invalid_controller = parsed_binary_evidence_mutation(
            oversize_postexec
        )
        validate_recorded_run_authority(
            final_fixture["run_argv_sha256"],
            oversized_invalid_controller, parsed_summary, parsed_claim,
        )
        expect_validation_error(
            "positive oversized post-exec binary " + receipt_name,
            lambda oversize_postexec=oversize_postexec:
            parsed_binary_evidence_mutation(
                oversize_postexec, positive=True
            ),
        )
    expect_validation_error(
        "started empty pre-exec binary image",
        lambda: parsed_binary_evidence_mutation(
            lambda binary_document: binary_document["stat_before"].__setitem__(
                "size", 0
            )
        ),
    )
    expect_validation_error(
        "started oversized pre-exec binary image",
        lambda: parsed_binary_evidence_mutation(
            lambda binary_document: binary_document["stat_before"].__setitem__(
                "size", MAX_BINARY_BYTES + 1
            )
        ),
    )
    expect_validation_error(
        "post-exec path stat without held-file hash",
        lambda: parsed_binary_evidence_mutation(
            lambda binary_document: binary_document.__setitem__(
                "sha256_after", None
            )
        ),
    )

    validate_retained_claim_controller_binding(
        parsed_claim, final_fixture["claim_receipt"], parsed_controller,
        parsed_summary,
    )
    for label, mutate in (
        (
            "source changed then restored",
            lambda summary_document: summary_document[
                "source_manifest_after"
            ]["entries"][0]["stat"].__setitem__(
                "mtime_ns",
                summary_document["source_manifest_after"]["entries"][0][
                    "stat"
                ]["mtime_ns"] + 1,
            ),
        ),
        (
            "Git executable changed then restored",
            lambda summary_document: summary_document["git_after"][
                "executable"
            ]["stat"].__setitem__(
                "mtime_ns",
                summary_document["git_after"]["executable"]["stat"][
                    "mtime_ns"
                ] + 1,
            ),
        ),
    ):
        changed_worker_summary_document = json.loads(
            canonical_bytes(parsed_summary).decode("ascii")
        )
        mutate(changed_worker_summary_document)
        changed_worker_summary_document["summary_preimage_sha256"] = None
        make_summary_preimage(changed_worker_summary_document)
        changed_worker_summary = parse_summary_bytes(
            canonical_bytes(changed_worker_summary_document) + b"\n"
        )
        validate_retained_claim_controller_binding(
            parsed_claim, final_fixture["claim_receipt"], parsed_controller,
            changed_worker_summary,
        )
        positive_controller = dict(parsed_controller)
        positive_controller["outcome"] = "pass"
        expect_validation_error(
            "positive retained " + label,
            lambda positive_controller=positive_controller,
            changed_worker_summary=changed_worker_summary:
            validate_retained_claim_controller_binding(
                parsed_claim, final_fixture["claim_receipt"],
                positive_controller, changed_worker_summary,
            ),
        )
    validate_retained_science(
        b"", b"", health["thermal"]["window_csv_ascii"].encode("ascii"),
        parsed_summary,
    )
    huge_integer_raw = (
        b'{"value":' + b"1" * 5000 + b"}\n"
        + b"{}\n" * (EXPECTED_PANEL_COUNT + 1)
    )
    huge_integer_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    huge_integer_summary.update({
        "config": {}, "raw_bytes": len(huge_integer_raw),
        "raw_complete": False,
        "raw_record_count": huge_integer_raw.count(b"\n"),
        "raw_sha256": sha256_bytes(huge_integer_raw), "statistics": {},
        "summary_preimage_sha256": None,
    })
    huge_integer_summary["summary_preimage_sha256"] = sha256_bytes(
        canonical_bytes(huge_integer_summary)
    )
    parsed_huge_integer_summary = parse_summary_bytes(
        canonical_bytes(huge_integer_summary) + b"\n"
    )
    validate_retained_science(
        huge_integer_raw, b"",
        health["thermal"]["window_csv_ascii"].encode("ascii"),
        parsed_huge_integer_summary,
    )
    validate_controller_bundle_binding(
        parsed_controller, parsed_summary, final_fixture["summary_bytes"],
        final_fixture["claim_receipt"], binary,
    )

    impossible_collection_summary_document = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    impossible_collection_health = impossible_collection_summary_document[
        "health"
    ]
    impossible_collection_health = finalize_health_receipt(
        impossible_collection_health,
        ["health admission: impossible postclaim collection"],
        impossible_collection_health["violations"],
    )
    impossible_collection_summary_document["health"] = (
        impossible_collection_health
    )
    impossible_collection_summary_document["summary_preimage_sha256"] = None
    make_summary_preimage(impossible_collection_summary_document)
    impossible_collection_summary_bytes = (
        canonical_bytes(impossible_collection_summary_document) + b"\n"
    )
    impossible_collection_summary = parse_summary_bytes(
        impossible_collection_summary_bytes
    )
    impossible_collection_controller_document = json.loads(
        canonical_bytes(parsed_controller).decode("ascii")
    )
    impossible_collection_controller_document["health_receipt_sha256"] = (
        impossible_collection_health["receipt_sha256"]
    )
    impossible_collection_controller_document["summary"] = {
        "outcome": impossible_collection_summary["outcome"],
        "preimage_sha256": impossible_collection_summary[
            "summary_preimage_sha256"
        ],
        "sha256": sha256_bytes(impossible_collection_summary_bytes),
    }
    impossible_collection_controller_document["artifacts"]["summary.json"][
        "bytes"
    ] = len(impossible_collection_summary_bytes)
    impossible_collection_controller_document["artifacts"]["summary.json"][
        "sha256"
    ] = sha256_bytes(impossible_collection_summary_bytes)
    impossible_collection_controller = parse_controller_document(
        controller_document_bytes(impossible_collection_controller_document)
    )
    expect_validation_error(
        "claimed bundle impossible health collection prefix",
        lambda: validate_controller_bundle_binding(
            impossible_collection_controller, impossible_collection_summary,
            impossible_collection_summary_bytes,
            final_fixture["claim_receipt"], binary,
        ),
    )

    def invalid_causality_mutation(
        mutate_summary: Any, mutate_controller: Any,
    ) -> None:
        summary_document = json.loads(
            canonical_bytes(parsed_summary).decode("ascii")
        )
        mutate_summary(summary_document)
        summary_document["summary_preimage_sha256"] = None
        make_summary_preimage(summary_document)
        mutated_summary_bytes = canonical_bytes(summary_document) + b"\n"
        mutated_summary = parse_summary_bytes(mutated_summary_bytes)
        controller_document = json.loads(
            canonical_bytes(parsed_controller).decode("ascii")
        )
        mutate_controller(controller_document)
        controller_document["summary"] = {
            "outcome": mutated_summary["outcome"],
            "preimage_sha256": mutated_summary[
                "summary_preimage_sha256"
            ],
            "sha256": sha256_bytes(mutated_summary_bytes),
        }
        controller_document["artifacts"]["summary.json"]["bytes"] = (
            len(mutated_summary_bytes)
        )
        controller_document["artifacts"]["summary.json"]["sha256"] = (
            sha256_bytes(mutated_summary_bytes)
        )
        mutated_controller = parse_controller_document(
            controller_document_bytes(controller_document)
        )
        validate_controller_bundle_binding(
            mutated_controller, mutated_summary, mutated_summary_bytes,
            final_fixture["claim_receipt"], binary,
        )

    expect_validation_error(
        "invalid bundle empty publication-failure roster",
        lambda: invalid_causality_mutation(
            lambda document: document.__setitem__(
                "publication_failures", []
            ),
            lambda _document: None,
        ),
    )
    for field, reserved in (
        ("summary", "none"), ("controller", "statistical-gates"),
    ):
        expect_validation_error(
            "invalid bundle reserved " + field + " failure",
            lambda field=field, reserved=reserved:
            invalid_causality_mutation(
                (
                    lambda document: document.__setitem__(
                        "failure", reserved
                    )
                    if field == "summary" else None
                ),
                (
                    lambda document: document.__setitem__(
                        "failure", reserved
                    )
                    if field == "controller" else None
                ),
            ),
        )
    for field in (
        "admission_sibling_ticks", "target_core", "controller_core",
        "sampler_core", "target_threads", "controller_initial_affinity",
    ):
        missing_preclaim_summary = json.loads(
            canonical_bytes(parsed_summary).decode("ascii")
        )
        missing_preclaim_summary["health"][field] = (
            {} if field == "admission_sibling_ticks" else []
        )
        expect_validation_error(
            "claimed bundle missing preclaim health " + field,
            lambda missing_preclaim_summary=missing_preclaim_summary:
            validate_controller_bundle_binding(
                parsed_controller, missing_preclaim_summary,
                final_fixture["summary_bytes"],
                final_fixture["claim_receipt"], binary,
            ),
        )
    missing_loader_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    missing_loader_summary["health_module_loader"] = {}
    expect_validation_error(
        "claimed bundle missing sealed health loader",
        lambda: validate_controller_bundle_binding(
            parsed_controller, missing_loader_summary,
            final_fixture["summary_bytes"],
            final_fixture["claim_receipt"], binary,
        ),
    )
    different_loader_summary = json.loads(
        canonical_bytes(parsed_summary).decode("ascii")
    )
    different_loader_summary["health_module_loader"]["modules"][0][
        "sha256"
    ] = "0" * 64
    expect_validation_error(
        "claimed bundle loader/source mismatch",
        lambda: validate_controller_bundle_binding(
            parsed_controller, different_loader_summary,
            final_fixture["summary_bytes"],
            final_fixture["claim_receipt"], binary,
        ),
    )
    for label, artifact_name, artifact_field, replacement in (
        ("raw artifact hash", "raw.jsonl", "sha256", "0" * 64),
        (
            "stderr artifact size", "stderr.txt", "bytes",
            parsed_summary["stderr_bytes"] + 1,
        ),
        ("thermal artifact hash", "thermal.csv", "sha256", "0" * 64),
        (
            "summary artifact size", "summary.json", "bytes",
            len(final_fixture["summary_bytes"]) + 1,
        ),
    ):
        mismatched_artifact_controller_document = json.loads(
            canonical_bytes(parsed_controller).decode("ascii")
        )
        mismatched_artifact_controller_document["artifacts"][artifact_name][
            artifact_field
        ] = replacement
        mismatched_artifact_controller = parse_controller_document(
            controller_document_bytes(mismatched_artifact_controller_document)
        )
        expect_validation_error(
            "controller/summary " + label,
            lambda mismatched_artifact_controller=(
                mismatched_artifact_controller
            ): validate_controller_bundle_binding(
                mismatched_artifact_controller, parsed_summary,
                final_fixture["summary_bytes"],
                final_fixture["claim_receipt"], binary,
            ),
        )
    def set_admission_start(
        summary_document: Dict[str, Any], timestamp: int,
    ) -> None:
        admission_receipt = summary_document["health"][
            "admission_sibling_ticks"
        ]
        admission_receipt["interval_start_monotonic_ns"] = timestamp
        admission_receipt["start"]["read_monotonic_ns"] = timestamp

    def set_admission_end(
        summary_document: Dict[str, Any], timestamp: int,
    ) -> None:
        admission_receipt = summary_document["health"][
            "admission_sibling_ticks"
        ]
        admission_receipt["interval_end_monotonic_ns"] = timestamp
        admission_receipt["end"]["read_monotonic_ns"] = timestamp

    for label, mutate in (
        (
            "admission endpoint receipt mismatch",
            lambda summary_document: summary_document["health"][
                "admission_sibling_ticks"
            ].__setitem__(
                "interval_start_monotonic_ns",
                summary_document["health"]["admission_sibling_ticks"][
                    "start"
                ]["read_monotonic_ns"] + 1,
            ),
        ),
        (
            "admission before controller origin",
            lambda summary_document: set_admission_start(
                summary_document,
                parsed_controller["controller_started_monotonic_ns"] - 1,
            ),
        ),
        (
            "admission after execution start",
            lambda summary_document: set_admission_end(
                summary_document,
                binary["execution_started_monotonic_ns"] + 1,
            ),
        ),
        (
            "admission after controller observation",
            lambda summary_document: set_admission_end(
                summary_document,
                parsed_controller["controller_observed_monotonic_ns"] + 1,
            ),
        ),
        (
            "admission exact ten-second boundary",
            lambda summary_document: set_admission_end(
                summary_document,
                parsed_controller["controller_started_monotonic_ns"]
                + 10_000_000_000,
            ),
        ),
    ):
        impossible_admission_summary = json.loads(
            canonical_bytes(parsed_summary).decode("ascii")
        )
        mutate(impossible_admission_summary)
        expect_validation_error(
            label,
            lambda impossible_admission_summary=impossible_admission_summary:
            validate_controller_bundle_binding(
                parsed_controller, impossible_admission_summary,
                final_fixture["summary_bytes"],
                final_fixture["claim_receipt"], binary,
            ),
        )
    for label, mutate in (
        (
            "controller/summary target CPU mismatch",
            lambda document: document.__setitem__(
                "target_cpu", EXPECTED_CONTROLLER_CPU
            ),
        ),
        (
            "controller/summary signal roster mismatch",
            lambda document: document.__setitem__("signals", ["SIGTERM"]),
        ),
    ):
        mismatched_controller = json.loads(
            canonical_bytes(parsed_controller).decode("ascii")
        )
        mutate(mismatched_controller)
        expect_validation_error(
            label,
            lambda mismatched_controller=mismatched_controller:
            validate_controller_bundle_binding(
                mismatched_controller, parsed_summary,
                final_fixture["summary_bytes"],
                final_fixture["claim_receipt"], binary,
            ),
        )
    complete_bytes = canonical_bytes(final_fixture["complete"]) + b"\n"
    parsed_complete = parse_complete_document(complete_bytes)
    validate_complete_bundle_binding(
        parsed_complete, parsed_controller,
        final_fixture["controller_bytes"], final_fixture["claim_receipt"],
    )

    for label, evidence_health in (
        ("terminal EDAC abort", terminal_health),
        (
            "terminal final timestamp acquisition failure",
            terminal_timestamp_failure_health,
        ),
        (
            "terminal PID unlink post-removal failure",
            pid_unlink_after_removal_health,
        ),
        (
            "terminal sampler-error partial validation",
            sampler_error_suffix_health,
        ),
        *seal_recovery_healths,
        ("stale raw-ahead partial-validation evidence", stale_health),
        ("coverage-cutoff paired-prefix evidence", coverage_cutoff_health),
        ("stale evidence recovery failure", stale_recovery_failure_health),
    ):
        evidence_fixture = synthetic_final_bundle_fixture(
            evidence_health, binary, commit
        )
        evidence_claim = parse_claim_document(evidence_fixture["claim_raw"])
        evidence_summary = parse_summary_bytes(
            evidence_fixture["summary_bytes"]
        )
        evidence_controller = parse_controller_document(
            evidence_fixture["controller_bytes"]
        )
        validate_recorded_run_authority(
            evidence_fixture["run_argv_sha256"], evidence_controller,
            evidence_summary, evidence_claim,
        )
        validate_controller_bundle_binding(
            evidence_controller, evidence_summary,
            evidence_fixture["summary_bytes"],
            evidence_fixture["claim_receipt"], binary,
        )
        evidence_complete = parse_complete_document(
            canonical_bytes(evidence_fixture["complete"]) + b"\n"
        )
        validate_complete_bundle_binding(
            evidence_complete, evidence_controller,
            evidence_fixture["controller_bytes"],
            evidence_fixture["claim_receipt"],
        )
        validate_retained_claim_controller_binding(
            evidence_claim, evidence_fixture["claim_receipt"],
            evidence_controller, evidence_summary,
        )
        validate_retained_science(
            b"", b"", health_thermal_artifact_bytes(evidence_health),
            evidence_summary,
        )
        if evidence_controller["outcome"] != "invalid":
            raise AssertionError(label + " did not remain invalid")

    mutated_claim = json.loads(canonical_bytes(parsed_claim).decode("ascii"))
    mutated_claim["binary_sha256"] = "0" * 64
    mutated_claim_raw = canonical_bytes(mutated_claim) + b"\n"
    mutated_claim_parsed = parse_claim_document(mutated_claim_raw)
    expect_validation_error(
        "claim expected hash binding",
        lambda: validate_claim_runtime_binding(
            mutated_claim_parsed, final_fixture["claim_receipt"],
            final_fixture["claim_document"], mutated_claim_raw,
        ),
    )
    mutated_claim_id = json.loads(canonical_bytes(parsed_claim).decode("ascii"))
    mutated_claim_id["parent_inode"] += 1
    mutated_claim_id_raw = canonical_bytes(mutated_claim_id) + b"\n"
    expect_validation_error(
        "claim runtime identity binding",
        lambda: validate_claim_runtime_binding(
            parse_claim_document(mutated_claim_id_raw),
            final_fixture["claim_receipt"], final_fixture["claim_document"],
            mutated_claim_id_raw,
        ),
    )
    mutated_claim_owner = json.loads(
        canonical_bytes(final_fixture["claim_receipt"]).decode("ascii")
    )
    mutated_claim_owner["stat"]["uid"] += 1
    expect_validation_error(
        "claim owner binding",
        lambda: validate_claim_runtime_binding(
            parsed_claim, mutated_claim_owner,
            final_fixture["claim_document"], final_fixture["claim_raw"],
        ),
    )

    def mutated_controller_bytes(mutate: Any) -> bytes:
        document = json.loads(
            canonical_bytes(final_fixture["controller"]).decode("ascii")
        )
        mutate(document)
        return controller_document_bytes(document)

    def changed_final_build(document: Dict[str, Any]) -> None:
        final = document["build_after"]
        final["root"] = "/tmp/changed-synthetic-build"
        final["sha256"] = None
        final["sha256"] = sha256_bytes(canonical_bytes(final))

    def changed_final_controller(document: Dict[str, Any]) -> None:
        final = document["controller_after"]
        final["pid"] += 1
        final["receipt_sha256"] = None
        final["receipt_sha256"] = sha256_bytes(canonical_bytes(final))

    def changed_final_git(document: Dict[str, Any]) -> None:
        document["git_after"]["tracked_index_flags_sha256"] = "0" * 64

    def changed_final_source(document: Dict[str, Any]) -> None:
        final = document["source_manifest_after"]
        final["entries"][0]["sha256"] = "0" * 64
        preimage = b"".join(
            (
                entry["sha256"] + "  " + entry["path"] + "\n"
            ).encode("ascii")
            for entry in final["entries"]
        )
        final["sha256"] = sha256_bytes(preimage)

    for label, mutate in (
        ("changed nonempty final build", changed_final_build),
        ("changed nonempty final controller", changed_final_controller),
        ("changed nonempty final Git", changed_final_git),
        ("changed nonempty final source", changed_final_source),
    ):
        expect_validation_error(
            label,
            lambda mutate=mutate: parse_controller_document(
                mutated_controller_bytes(mutate)
            ),
        )

    missing_final_sampler = parse_controller_document(
        mutated_controller_bytes(
            lambda document: document.__setitem__("sampler_after", {})
        )
    )
    validate_controller_bundle_binding(
        missing_final_sampler, parsed_summary, final_fixture["summary_bytes"],
        final_fixture["claim_receipt"], binary,
    )

    health_mismatch_bytes = mutated_controller_bytes(
        lambda document: document.__setitem__(
            "health_receipt_sha256", "0" * 64
        )
    )
    health_mismatch = parse_controller_document(health_mismatch_bytes)
    expect_validation_error(
        "controller health binding",
        lambda: validate_controller_bundle_binding(
            health_mismatch, parsed_summary, final_fixture["summary_bytes"],
            final_fixture["claim_receipt"], binary,
        ),
    )
    summary_mismatch_bytes = mutated_controller_bytes(
        lambda document: document["summary"].__setitem__(
            "preimage_sha256", "0" * 64
        )
    )
    summary_mismatch = parse_controller_document(summary_mismatch_bytes)
    expect_validation_error(
        "controller summary preimage binding",
        lambda: validate_controller_bundle_binding(
            summary_mismatch, parsed_summary, final_fixture["summary_bytes"],
            final_fixture["claim_receipt"], binary,
        ),
    )

    def contradict_git_source(document: Dict[str, Any]) -> None:
        document["controller_before"]["git_sha256"] = "5" * 64
        document["controller_before"]["receipt_sha256"] = None
        document["controller_before"]["receipt_sha256"] = sha256_bytes(
            canonical_bytes(document["controller_before"])
        )

    expect_validation_error(
        "controller source/Git contradiction",
        lambda: parse_controller_document(
            mutated_controller_bytes(contradict_git_source)
        ),
    )

    def redirected_worktree(document: Dict[str, Any]) -> None:
        document["git_before"]["worktree_root"] = "/tmp/redirected-worktree"

    expect_validation_error(
        "redirected Git worktree",
        lambda: parse_controller_document(
            mutated_controller_bytes(redirected_worktree)
        ),
    )

    def contradict_source_blob(document: Dict[str, Any]) -> None:
        blob = document["git_before"]["source_blob_oids"][0]
        blob["head_oid"] = "5" * 40
        blob["worktree_oid"] = "5" * 40
        document["git_before"]["source_blob_roster_sha256"] = sha256_bytes(
            canonical_bytes(document["git_before"]["source_blob_oids"])
        )

    expect_validation_error(
        "controller raw source/HEAD blob contradiction",
        lambda: parse_controller_document(
            mutated_controller_bytes(contradict_source_blob)
        ),
    )

    def writable_source(document: Dict[str, Any]) -> None:
        document["source_manifest_before"]["entries"][0]["stat"][
            "mode"
        ] = 0o664

    expect_validation_error(
        "writable source provenance",
        lambda: parse_controller_document(
            mutated_controller_bytes(writable_source)
        ),
    )
    oversized_source = json.loads(
        canonical_bytes(final_fixture["source"]).decode("ascii")
    )
    oversized_source["entries"][0]["bytes"] = MAX_SOURCE_FILE_BYTES + 1
    oversized_source["entries"][0]["stat"]["size"] = (
        MAX_SOURCE_FILE_BYTES + 1
    )
    expect_validation_error(
        "oversized source receipt",
        lambda: validate_source_manifest_receipt(
            oversized_source, "oversized synthetic source"
        ),
    )
    for label, oid in (
        ("short source Git object ID", "a" * 39),
        ("long source Git object ID", "a" * 41),
        ("uppercase source Git object ID", "A" * 40),
    ):
        malformed_oid_source = json.loads(
            canonical_bytes(final_fixture["source"]).decode("ascii")
        )
        malformed_oid_source["entries"][0]["git_blob_oid"] = oid
        expect_validation_error(
            label,
            lambda malformed_oid_source=malformed_oid_source:
            validate_source_manifest_receipt(
                malformed_oid_source, "malformed synthetic source OID"
            ),
        )

    for label, field, value in (
        ("controller optimize bool alias", "optimize", False),
        ("controller affinity float alias", "singleton_affinity", [121.0]),
    ):
        def mutate_provenance(
            document: Dict[str, Any], field: str = field, value: Any = value,
        ) -> None:
            document["controller_before"][field] = value
            document["controller_before"]["receipt_sha256"] = None
            document["controller_before"]["receipt_sha256"] = sha256_bytes(
                canonical_bytes(document["controller_before"])
            )

        expect_validation_error(
            label,
            lambda mutate_provenance=mutate_provenance: parse_controller_document(
                mutated_controller_bytes(mutate_provenance)
            ),
        )

    def extra_controller_environment(document: Dict[str, Any]) -> None:
        document["controller_before"]["environment"]["LD_PRELOAD"] = "/tmp/x"
        document["controller_before"]["receipt_sha256"] = None
        document["controller_before"]["receipt_sha256"] = sha256_bytes(
            canonical_bytes(document["controller_before"])
        )

    expect_validation_error(
        "controller startup environment extension",
        lambda: parse_controller_document(
            mutated_controller_bytes(extra_controller_environment)
        ),
    )

    mutated_sampler_environment = json.loads(
        canonical_bytes(health["sampler"]).decode("ascii")
    )
    mutated_sampler_environment["environment"].pop("TZ")
    mutated_sampler_environment["environment_sha256"] = sha256_bytes(
        canonical_bytes(mutated_sampler_environment["environment"])
    )
    expect_validation_error(
        "sampler startup environment omission",
        lambda: validate_sampler_receipt(mutated_sampler_environment),
    )

    for label, field, integer_value in (
        ("controller deadline integer alias", "controller_deadline_seconds", 110),
        ("outer deadline integer alias", "outer_deadline_seconds", 120),
        ("controller elapsed integer alias", "controller_elapsed_seconds", 2),
    ):
        expect_validation_error(
            label,
            lambda field=field, integer_value=integer_value:
            parse_controller_document(
                mutated_controller_bytes(
                    lambda document: document.__setitem__(field, integer_value)
                )
            ),
        )

    def set_controller_observation(
        document: Dict[str, Any], observed_ns: int,
    ) -> None:
        started_ns = document["controller_started_monotonic_ns"]
        document["controller_observed_monotonic_ns"] = observed_ns
        document["controller_elapsed_seconds"] = (
            observed_ns - started_ns
        ) / 1_000_000_000.0

    def set_invalid_reap_timeout(
        document: Dict[str, Any], reap_ns: int, timed_out: bool,
    ) -> None:
        document["binary"]["child_reap_monotonic_ns"] = reap_ns
        document["binary"]["execution_finished_monotonic_ns"] = reap_ns
        document["binary"]["timed_out"] = timed_out
        set_controller_observation(document, reap_ns)

    timeout_boundary = (
        final_fixture["controller"]["controller_started_monotonic_ns"]
        + int(CONTROLLER_DEADLINE_SECONDS * 1_000_000_000)
    )
    parse_controller_document(
        mutated_controller_bytes(
            lambda document: set_invalid_reap_timeout(
                document, timeout_boundary - 1, True
            )
        )
    )
    expect_validation_error(
        "late invalid child reap without timeout classification",
        lambda: parse_controller_document(
            mutated_controller_bytes(
                lambda document: set_invalid_reap_timeout(
                    document, timeout_boundary, False
                )
            )
        ),
    )

    expect_validation_error(
        "binary execution start before controller origin",
        lambda: parse_controller_document(
            mutated_controller_bytes(
                lambda document: document["binary"].__setitem__(
                    "execution_started_monotonic_ns",
                    document["controller_started_monotonic_ns"] - 1,
                )
            )
        ),
    )
    expect_validation_error(
        "binary execution finish after controller observation",
        lambda: parse_controller_document(
            mutated_controller_bytes(
                lambda document: document["binary"].__setitem__(
                    "execution_finished_monotonic_ns",
                    document["controller_observed_monotonic_ns"] + 1,
                )
            )
        ),
    )
    for label, sampler_start in (
        (
            "sampler start before controller origin",
            final_fixture["controller"]["controller_started_monotonic_ns"] - 1,
        ),
        (
            "sampler start after execution start",
            binary["execution_started_monotonic_ns"] + 1,
        ),
    ):
        expect_validation_error(
            label,
            lambda sampler_start=sampler_start: parse_controller_document(
                mutated_controller_bytes(
                    lambda document: document["sampler_after"].__setitem__(
                        "window_start_monotonic_ns", sampler_start
                    )
                )
            ),
        )

    for label, elapsed_value in (
        ("invalid controller at final-start deadline", 117.0),
        ("invalid controller after final-start deadline", 118.0),
    ):
        expect_validation_error(
            label,
            lambda elapsed_value=elapsed_value: parse_controller_document(
                mutated_controller_bytes(
                    lambda document: set_controller_observation(
                        document,
                        document["controller_started_monotonic_ns"]
                        + int(elapsed_value * 1_000_000_000),
                    )
                )
            ),
        )

    expect_validation_error(
        "controller derived elapsed mismatch",
        lambda: parse_controller_document(
            mutated_controller_bytes(
                lambda document: document.__setitem__(
                    "controller_elapsed_seconds", 1.0
                )
            )
        ),
    )

    for label, observed_ns in (
        ("controller observation before child reap", 1_700_000_000),
        ("controller observation before sampler endpoint", 1_900_000_000),
    ):
        expect_validation_error(
            label,
            lambda observed_ns=observed_ns: parse_controller_document(
                mutated_controller_bytes(
                    lambda document: set_controller_observation(
                        document, observed_ns
                    )
                )
            ),
        )

    def impossible_short_elapsed(document: Dict[str, Any]) -> None:
        started_ns = document["controller_started_monotonic_ns"]
        document["outcome"] = "pass"
        document["failure"] = "none"
        document["summary"]["outcome"] = "pass"
        document["binary"]["child_reap_monotonic_ns"] = (
            started_ns + 109_000_000_000
        )
        document["sampler_after"]["window_end_monotonic_ns"] = (
            started_ns + 114_000_000_000
        )
        document["controller_observed_monotonic_ns"] = (
            started_ns + 114_000_000_000
        )
        document["controller_elapsed_seconds"] = 1.0

    expect_validation_error(
        "impossible short controller elapsed for late evidence",
        lambda: parse_controller_document(
            mutated_controller_bytes(impossible_short_elapsed)
        ),
    )

    def positive_controller_at_reap(reap_ns: int) -> bytes:
        def mutate(document: Dict[str, Any]) -> None:
            document["outcome"] = "pass"
            document["failure"] = "none"
            document["summary"]["outcome"] = "pass"
            document["binary"]["child_reap_monotonic_ns"] = reap_ns
            document["binary"]["execution_finished_monotonic_ns"] = reap_ns
            set_controller_observation(document, reap_ns)

        return mutated_controller_bytes(mutate)

    controller_start_ns = final_fixture["controller"][
        "controller_started_monotonic_ns"
    ]
    parse_controller_document(
        positive_controller_at_reap(
            controller_start_ns + 109_999_999_999
        )
    )
    expect_validation_error(
        "positive controller inclusive 110-second reap",
        lambda: parse_controller_document(
            positive_controller_at_reap(
                controller_start_ns + 110_000_000_000
            )
        ),
    )

    def positive_controller_at_health_endpoint(endpoint_ns: int) -> bytes:
        def mutate(document: Dict[str, Any]) -> None:
            document["outcome"] = "pass"
            document["failure"] = "none"
            document["summary"]["outcome"] = "pass"
            document["sampler_after"]["window_end_monotonic_ns"] = endpoint_ns
            document["binary"]["execution_finished_monotonic_ns"] = endpoint_ns
            set_controller_observation(document, endpoint_ns)

        return mutated_controller_bytes(mutate)

    quantized_health_boundary = (
        controller_start_ns
        + int(HEALTH_COLLECTION_DEADLINE_SECONDS * 1_000_000_000)
    )
    # The CSV seconds field is rounded to integer nanoseconds.  Equality can
    # represent a sample taken just before the live strict float deadline.
    parse_controller_document(
        positive_controller_at_health_endpoint(quantized_health_boundary)
    )
    expect_validation_error(
        "positive controller post-health-deadline endpoint",
        lambda: parse_controller_document(
            positive_controller_at_health_endpoint(
                quantized_health_boundary + 1
            )
        ),
    )

    def positive_late_controller(document: Dict[str, Any]) -> None:
        document["outcome"] = "pass"
        document["failure"] = "none"
        document["summary"]["outcome"] = "pass"
        set_controller_observation(
            document,
            document["controller_started_monotonic_ns"]
            + int(FINAL_COMMIT_START_DEADLINE_SECONDS * 1_000_000_000),
        )

    expect_validation_error(
        "positive controller final-start deadline",
        lambda: parse_controller_document(
            mutated_controller_bytes(positive_late_controller)
        ),
    )

    late_summary_document = json.loads(
        canonical_bytes(final_fixture["summary"]).decode("ascii")
    )
    late_summary_document["binary"][
        "execution_finished_monotonic_ns"
    ] = 3_600_000_000
    late_summary_document["elapsed_seconds"] = 2.5
    late_summary_document["summary_preimage_sha256"] = None
    late_summary_document["summary_preimage_sha256"] = sha256_bytes(
        canonical_bytes(late_summary_document)
    )
    late_summary_bytes = canonical_bytes(late_summary_document) + b"\n"
    late_summary = parse_summary_bytes(late_summary_bytes)
    early_controller = json.loads(
        canonical_bytes(parsed_controller).decode("ascii")
    )
    early_controller["binary"] = late_summary["binary"]
    early_controller["summary"] = {
        "outcome": late_summary["outcome"],
        "preimage_sha256": late_summary["summary_preimage_sha256"],
        "sha256": sha256_bytes(late_summary_bytes),
    }
    expect_validation_error(
        "controller predates summary elapsed time",
        lambda: validate_controller_bundle_binding(
            early_controller, late_summary, late_summary_bytes,
            final_fixture["claim_receipt"], late_summary["binary"],
        ),
    )

    def positive_without_final_sampler(document: Dict[str, Any]) -> None:
        document["outcome"] = "pass"
        document["failure"] = "none"
        document["summary"]["outcome"] = "pass"
        document["sampler_after"] = {}

    expect_validation_error(
        "positive controller missing final sampler",
        lambda: parse_controller_document(
            mutated_controller_bytes(positive_without_final_sampler)
        ),
    )

    invalid_late_health_summary = parsed_summary
    invalid_late_health_controller = json.loads(
        canonical_bytes(parsed_controller).decode("ascii")
    )
    invalid_late_health_bytes = final_fixture["summary_bytes"]
    invalid_late_health_controller["controller_observed_monotonic_ns"] = (
        invalid_late_health_summary["health"]["sampler"][
            "window_end_monotonic_ns"
        ] - 1
    )
    invalid_late_health_controller["sampler_after"] = {}
    expect_validation_error(
        "invalid health endpoint after controller observation",
        lambda: validate_controller_bundle_binding(
            invalid_late_health_controller, invalid_late_health_summary,
            invalid_late_health_bytes, final_fixture["claim_receipt"], binary,
        ),
    )

    for field, replacement in (
        ("claim_sha256", "0" * 64),
        ("controller_sha256", "0" * 64),
        ("controller_outcome", "pass"),
    ):
        mutated_complete = dict(final_fixture["complete"])
        mutated_complete[field] = replacement
        mutated_complete_parsed = parse_complete_document(
            canonical_bytes(mutated_complete) + b"\n"
        )
        expect_validation_error(
            "COMPLETE " + field + " binding",
            lambda mutated_complete_parsed=mutated_complete_parsed:
            validate_complete_bundle_binding(
                mutated_complete_parsed, parsed_controller,
                final_fixture["controller_bytes"],
                final_fixture["claim_receipt"],
            ),
        )
    integer_complete_elapsed = dict(final_fixture["complete"])
    integer_complete_elapsed["elapsed_seconds_before_commit"] = 3
    expect_validation_error(
        "COMPLETE elapsed integer alias",
        lambda: parse_complete_document(
            canonical_bytes(integer_complete_elapsed) + b"\n"
        ),
    )
    early_complete = dict(final_fixture["complete"])
    early_complete["complete_observed_monotonic_ns"] = 2_500_000_000
    early_complete["elapsed_seconds_before_commit"] = 1.5
    expect_validation_error(
        "COMPLETE predates controller elapsed time",
        lambda: validate_complete_bundle_binding(
            parse_complete_document(canonical_bytes(early_complete) + b"\n"),
            parsed_controller, final_fixture["controller_bytes"],
            final_fixture["claim_receipt"],
        ),
    )
    mismatched_complete_elapsed = dict(final_fixture["complete"])
    mismatched_complete_elapsed["elapsed_seconds_before_commit"] = 2.5
    expect_validation_error(
        "COMPLETE derived elapsed mismatch",
        lambda: validate_complete_bundle_binding(
            parse_complete_document(
                canonical_bytes(mismatched_complete_elapsed) + b"\n"
            ),
            parsed_controller, final_fixture["controller_bytes"],
            final_fixture["claim_receipt"],
        ),
    )
    commit_margin_boundary = dict(final_fixture["complete"])
    commit_margin_boundary["complete_observed_monotonic_ns"] = (
        final_fixture["controller"]["controller_started_monotonic_ns"]
        + int(FINAL_COMMIT_START_DEADLINE_SECONDS * 1_000_000_000)
    )
    commit_margin_boundary["elapsed_seconds_before_commit"] = (
        FINAL_COMMIT_START_DEADLINE_SECONDS
    )
    validate_complete_bundle_binding(
        parse_complete_document(
            canonical_bytes(commit_margin_boundary) + b"\n"
        ),
        parsed_controller, final_fixture["controller_bytes"],
        final_fixture["claim_receipt"],
    )
    before_commit_margin = dict(final_fixture["complete"])
    before_commit_margin["complete_observed_monotonic_ns"] = (
        commit_margin_boundary["complete_observed_monotonic_ns"] - 1
    )
    before_commit_margin["elapsed_seconds_before_commit"] = (
        (before_commit_margin["complete_observed_monotonic_ns"]
         - final_fixture["controller"]["controller_started_monotonic_ns"])
        / 1_000_000_000.0
    )
    validate_complete_bundle_binding(
        parse_complete_document(canonical_bytes(before_commit_margin) + b"\n"),
        parsed_controller, final_fixture["controller_bytes"],
        final_fixture["claim_receipt"],
    )
    before_outer_boundary = dict(final_fixture["complete"])
    before_outer_boundary["complete_observed_monotonic_ns"] = (
        final_fixture["controller"]["controller_started_monotonic_ns"]
        + int(OUTER_DEADLINE_SECONDS * 1_000_000_000) - 1
    )
    before_outer_boundary["elapsed_seconds_before_commit"] = (
        (before_outer_boundary["complete_observed_monotonic_ns"]
         - final_fixture["controller"]["controller_started_monotonic_ns"])
        / 1_000_000_000.0
    )
    validate_complete_bundle_binding(
        parse_complete_document(canonical_bytes(before_outer_boundary) + b"\n"),
        parsed_controller, final_fixture["controller_bytes"],
        final_fixture["claim_receipt"],
    )
    outer_boundary_complete = dict(final_fixture["complete"])
    outer_boundary_complete["complete_observed_monotonic_ns"] = (
        final_fixture["controller"]["controller_started_monotonic_ns"]
        + int(OUTER_DEADLINE_SECONDS * 1_000_000_000)
    )
    outer_boundary_complete["elapsed_seconds_before_commit"] = (
        OUTER_DEADLINE_SECONDS
    )
    expect_validation_error(
        "COMPLETE exact outer boundary",
        lambda: parse_complete_document(
            canonical_bytes(outer_boundary_complete) + b"\n"
        ),
    )
    integer_elapsed = dict(final_fixture["summary"])
    integer_elapsed["elapsed_seconds"] = 1
    integer_elapsed["summary_preimage_sha256"] = None
    integer_elapsed["summary_preimage_sha256"] = sha256_bytes(
        canonical_bytes(integer_elapsed)
    )
    expect_validation_error(
        "summary elapsed integer alias",
        lambda: parse_summary_bytes(canonical_bytes(integer_elapsed) + b"\n"),
    )
    first_sample = next(iter(statistics["groups"].values()))
    integer_sample = dict(first_sample)
    integer_sample["log_standard_error"] = int(
        integer_sample["log_standard_error"]
    )
    expect_validation_error(
        "sample-summary integer alias",
        lambda: validate_sample_summary(integer_sample, "mutated sample"),
    )
    integer_student_t = json.loads(canonical_bytes(statistics).decode("ascii"))
    integer_student_t["student_t_975_df11"] = 2
    expect_validation_error(
        "Student-t integer alias",
        lambda: validate_statistics_receipt(integer_student_t),
    )

    missing_sampler_cpu = json.loads(canonical_bytes(health).decode("ascii"))
    missing_sampler_cpu["controller_initial_affinity"].remove(
        EXPECTED_SAMPLER_CPU
    )
    missing_sampler_cpu["receipt_sha256"] = None
    missing_sampler_cpu["receipt_sha256"] = sha256_bytes(
        canonical_bytes(missing_sampler_cpu)
    )
    expect_validation_error(
        "health initial sampler CPU",
        lambda: validate_health_receipt(
            missing_sampler_cpu, binary, require_complete=False
        ),
    )
    try:
        parse_csv_physical_line(b'"utc",monotonic_s\n', "quoted row")
    except ValidationError:
        pass
    else:
        raise AssertionError("quoted noncanonical CSV was accepted")

    def synthetic_sampler_row(monotonic_text: str) -> bytes:
        fields = ["0"] * len(THERMAL_HEADER)
        fields[0] = "2026-01-01T00:00:00Z"
        fields[1] = monotonic_text
        return (",".join(fields) + "\n").encode("ascii")

    sampler_timestamp_fixture = (
        (",".join(THERMAL_HEADER) + "\n").encode("ascii")
        + synthetic_sampler_row("1.25")
        + synthetic_sampler_row("9223372036.0")
        + synthetic_sampler_row("9223372037.0")
        + synthetic_sampler_row("1e100")
    )
    _, bounded_sampler_timestamps = sampler_rows(
        sampler_timestamp_fixture
    )
    if (
        bounded_sampler_timestamps != {
            1: 1_250_000_000,
            2: 9_223_372_036_000_000_000,
        }
    ):
        raise AssertionError(
            "sampler timestamps escaped their exact int63 bound"
        )

    live_pair_args = argparse.Namespace(
        expected_sampler_script_sha256=health["sampler"]["script_sha256"],
        expected_sampler_uid=health["sampler"]["process_uid"],
    )
    live_raw = health["thermal"]["window_csv_ascii"].encode("ascii")
    live_validation = health["thermal"][
        "validation_jsonl_ascii"
    ].encode("ascii")
    if not validate_live_sampler_stream_pair(
        live_raw, live_validation, live_pair_args,
        "synthetic settled live sampler",
    ):
        raise AssertionError("paired live sampler streams did not settle")
    live_raw_lines = live_raw.splitlines(keepends=True)
    appended_row = parse_csv_physical_line(
        live_raw_lines[-1], "synthetic appended live raw row"
    )
    appended_row[0] = "2026-01-01T00:00:03Z"
    appended_row[1] = "3.000000"
    appended_raw_line = (",".join(appended_row) + "\n").encode("ascii")
    raw_ahead = live_raw + appended_raw_line
    if validate_live_sampler_stream_pair(
        raw_ahead, live_validation, live_pair_args,
        "synthetic raw-ahead live sampler",
    ):
        raise AssertionError("raw-ahead live sampler was treated as settled")
    appended_validation = parse_canonical_json_line(
        live_validation.splitlines(keepends=True)[-1],
        "synthetic appended live validation",
    )
    appended_validation["sample_index"] = 2
    appended_validation["monotonic_s"] = 3.0
    paired_continue_validation = (
        live_validation + canonical_bytes(appended_validation) + b"\n"
    )
    if not validate_live_sampler_stream_pair(
        raw_ahead, paired_continue_validation, live_pair_args,
        "synthetic delayed paired live sampler",
    ):
        raise AssertionError("delayed continue validation did not settle")
    appended_validation["decision"] = "edac_abort"
    terminal_live_validation = (
        live_validation + canonical_bytes(appended_validation) + b"\n"
    )
    expect_validation_error(
        "post-window live terminal validation",
        lambda: validate_live_sampler_stream_pair(
            raw_ahead, terminal_live_validation, live_pair_args,
            "synthetic post-window terminal live sampler",
        ),
    )

    with tempfile.TemporaryDirectory() as directory:
        settle_root = Path(directory)
        settle_raw = settle_root / "thermal.csv"
        settle_pid = settle_root / "sampler.pid"
        settle_validation = settle_root / "validation.jsonl"
        settle_receipt = settle_root / "receipt.json"
        settle_raw.write_bytes(raw_ahead)
        settle_pid.write_bytes((str(os.getpid()) + "\n").encode("ascii"))
        settle_validation_fd = os.open(
            settle_validation,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o600,
        )
        settle_receipt_fd = os.open(
            settle_receipt,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o600,
        )
        os.write(settle_validation_fd, live_validation)
        os.fchmod(settle_validation_fd, 0o444)
        os.fchmod(settle_receipt_fd, 0o444)
        os.chmod(settle_raw, 0o600)
        os.chmod(settle_pid, 0o444)
        settle_paths = (
            settle_raw, settle_pid, settle_validation, settle_receipt,
        )
        settle_infos = [
            os.stat(path, follow_symlinks=False) for path in settle_paths
        ]
        settle_args = argparse.Namespace(
            expected_sampler_csv_device=settle_infos[0].st_dev,
            expected_sampler_csv_inode=settle_infos[0].st_ino,
            expected_sampler_gid=os.getgid(),
            expected_sampler_i2c_gid=EXPECTED_SAMPLER_I2C_GID,
            expected_sampler_pid_file_device=settle_infos[1].st_dev,
            expected_sampler_pid_file_inode=settle_infos[1].st_ino,
            expected_sampler_process_start_ticks=process_start_ticks(
                os.getpid()
            ),
            expected_sampler_receipt_device=settle_infos[3].st_dev,
            expected_sampler_receipt_inode=settle_infos[3].st_ino,
            expected_sampler_script_sha256=health["sampler"][
                "script_sha256"
            ],
            expected_sampler_uid=os.getuid(),
            expected_sampler_validation_device=settle_infos[2].st_dev,
            expected_sampler_validation_inode=settle_infos[2].st_ino,
            sampler_csv=settle_raw,
            sampler_pid=os.getpid(),
            sampler_pid_file=settle_pid,
            sampler_receipt=settle_receipt,
            sampler_script=Path(health["sampler"]["script_path"]),
            sampler_validation_jsonl=settle_validation,
        )

        class SyntheticHealthApi:
            def __init__(self) -> None:
                self.calls = 0
                self.publish_final_intent = False

            def _verify_live_sampler_process(self, *_args: Any) -> None:
                self.calls += 1
                if self.publish_final_intent and self.calls == 2:
                    os.lseek(settle_receipt_fd, 0, os.SEEK_SET)
                    os.write(settle_receipt_fd, b"!")
                    os.fsync(settle_receipt_fd)

        settle_health_api = SyntheticHealthApi()
        original_read_process_security = globals()["read_process_security"]
        original_sched_getaffinity = os.sched_getaffinity
        original_monotonic_ns = time.monotonic_ns
        original_sleep = time.sleep
        delayed_validation_written = False
        appended_validation_line = paired_continue_validation[
            len(live_validation):
        ]

        def synthetic_process_security(
            _pid: int, _where: str,
        ) -> Dict[str, Any]:
            security = dict(health["sampler"]["process_security"])
            security["uids"] = [os.getuid()] * 4
            security["gids"] = [os.getgid()] * 4
            return security

        def publish_delayed_validation(_seconds: float) -> None:
            nonlocal delayed_validation_written
            if not delayed_validation_written:
                os.lseek(settle_validation_fd, 0, os.SEEK_END)
                os.write(settle_validation_fd, appended_validation_line)
                os.fsync(settle_validation_fd)
                delayed_validation_written = True

        try:
            globals()["read_process_security"] = synthetic_process_security
            os.sched_getaffinity = (  # type: ignore[assignment]
                lambda _pid: {EXPECTED_SAMPLER_CPU}
            )
            time.monotonic_ns = (  # type: ignore[assignment]
                lambda: 3_000_000_000 + SAMPLER_HEARTBEAT_MAX_GAP_NS
            )
            time.sleep = publish_delayed_validation  # type: ignore[assignment]
            settled = settle_live_sampler_streams(
                settle_args, settle_health_api, time.monotonic() + 1.0,
                "synthetic delayed live streams",
            )
            if (
                not delayed_validation_written
                or settled["csv_data"] != raw_ahead
                or settled["validation_data"] != paired_continue_validation
            ):
                raise AssertionError(
                    "live stream settler did not wait for delayed validation"
                )
            settle_health_api.calls = 0
            settle_health_api.publish_final_intent = True
            expect_validation_error(
                "terminal intent after paired stream snapshots",
                lambda: settle_live_sampler_streams(
                    settle_args, settle_health_api,
                    time.monotonic() + 1.0,
                    "synthetic final-intent live streams",
                ),
            )
            os.lseek(settle_receipt_fd, 0, os.SEEK_SET)
            os.ftruncate(settle_receipt_fd, 0)
            os.fsync(settle_receipt_fd)
            settle_health_api.calls = 0
            settle_health_api.publish_final_intent = False
            time.monotonic_ns = (  # type: ignore[assignment]
                lambda: (
                    3_000_000_000 + SAMPLER_HEARTBEAT_MAX_GAP_NS + 1
                )
            )
            expect_validation_error(
                "final live sampler stale heartbeat",
                lambda: settle_live_sampler_streams(
                    settle_args, settle_health_api,
                    time.monotonic() + 1.0,
                    "synthetic stale live streams",
                ),
            )
            os.lseek(settle_validation_fd, 0, os.SEEK_SET)
            os.ftruncate(settle_validation_fd, 0)
            os.write(settle_validation_fd, live_validation)
            os.fsync(settle_validation_fd)
            stale_raw_ahead = settle_live_sampler_streams(
                settle_args, settle_health_api, time.monotonic() + 1.0,
                "synthetic stale raw-ahead streams", require_fresh=False,
                allow_stale_suffixes=True,
            )
            if (
                stale_raw_ahead["csv_full_data"] != raw_ahead
                or stale_raw_ahead["validation_data"] != live_validation
                or stale_raw_ahead["paired_end_monotonic_ns"]
                != 2_000_000_000
            ):
                raise AssertionError(
                    "stale raw-ahead sampler prefix was not preserved"
                )
            partial_validation = (
                live_validation
                + appended_validation_line[:len(appended_validation_line) // 2]
            )
            os.lseek(settle_validation_fd, 0, os.SEEK_SET)
            os.ftruncate(settle_validation_fd, 0)
            os.write(settle_validation_fd, partial_validation)
            os.fsync(settle_validation_fd)
            stale_validation_partial = settle_live_sampler_streams(
                settle_args, settle_health_api, time.monotonic() + 1.0,
                "synthetic stale partial-validation streams",
                require_fresh=False, allow_stale_suffixes=True,
            )
            if (
                stale_validation_partial["validation_data"]
                != partial_validation
                or stale_validation_partial["paired_end_monotonic_ns"]
                != 2_000_000_000
            ):
                raise AssertionError(
                    "stale partial validation suffix was not preserved"
                )
            cutoff_sampler = json.loads(
                canonical_bytes(health["sampler"]).decode("ascii")
            )
            cutoff_sampler.update({
                "csv_bytes": len(stale_validation_partial["paired_raw"]),
                "csv_device": stale_validation_partial["csv_info"].st_dev,
                "csv_inode": stale_validation_partial["csv_info"].st_ino,
                "csv_path": str(settle_raw),
                "csv_sha256": sha256_bytes(
                    stale_validation_partial["paired_raw"]
                ),
                "csv_stat": stat_receipt(
                    stale_validation_partial["csv_info"]
                ),
                "pid": os.getpid(),
                "process_gid": os.getgid(),
                "process_start_ticks": (
                    settle_args.expected_sampler_process_start_ticks
                ),
                "process_uid": os.getuid(),
                "validation_jsonl": stale_validation_partial[
                    "validation_jsonl"
                ],
                "window_end_monotonic_ns": 2_000_000_000,
                "window_start_monotonic_ns": 1_000_000_000,
            })
            cutoff_capture = {
                "csv_full_data": stale_validation_partial["csv_full_data"],
                "csv_info": stale_validation_partial["csv_info"],
                "snapshot_observed_monotonic_ns": (
                    stale_validation_partial[
                        "snapshot_observed_monotonic_ns"
                    ]
                ),
                "validation_data": stale_validation_partial[
                    "validation_data"
                ],
                "validation_info": stale_validation_partial[
                    "validation_info"
                ],
            }
            with settle_raw.open("ab", buffering=0) as moving_raw:
                moving_raw.write(appended_raw_line)
                os.fsync(moving_raw.fileno())
            os.lseek(settle_validation_fd, 0, os.SEEK_END)
            os.write(settle_validation_fd, b"post-cutoff-growth")
            os.fsync(settle_validation_fd)
            cutoff_thermal, cutoff_collection, cutoff_violations = (
                tolerant_thermal_window(
                    settle_args, cutoff_sampler,
                    cutoff_sampler["window_start_monotonic_ns"],
                    cutoff_sampler["window_end_monotonic_ns"],
                    allow_stale_suffixes=True,
                    deadline=time.monotonic() + 1.0,
                    stale_stream_capture=cutoff_capture,
                )
            )
            if (
                cutoff_collection or cutoff_violations
                or cutoff_thermal["window_csv_ascii"]
                != health["thermal"]["window_csv_ascii"]
                or cutoff_thermal["validation_jsonl_ascii"]
                != health["thermal"]["validation_jsonl_ascii"]
            ):
                raise AssertionError(
                    "post-cutoff sampler growth discarded captured evidence"
                )
            expect_validation_error(
                "validation partial without preceding raw row",
                lambda: validate_stale_sampler_stream_prefix(
                    live_raw, partial_validation, live_pair_args,
                    "synthetic impossible validation-only partial",
                    time.monotonic_ns(),
                ),
            )
            partial_raw = live_raw + appended_raw_line[
                :len(appended_raw_line) // 2
            ]
            settle_raw.write_bytes(partial_raw)
            os.lseek(settle_validation_fd, 0, os.SEEK_SET)
            os.ftruncate(settle_validation_fd, 0)
            os.write(settle_validation_fd, live_validation)
            os.fsync(settle_validation_fd)
            stale_raw_partial = settle_live_sampler_streams(
                settle_args, settle_health_api, time.monotonic() + 1.0,
                "synthetic stale partial-raw streams", require_fresh=False,
                allow_stale_suffixes=True,
            )
            if stale_raw_partial["csv_full_data"] != partial_raw:
                raise AssertionError(
                    "stale partial raw suffix was not preserved"
                )
        finally:
            time.sleep = original_sleep  # type: ignore[assignment]
            time.monotonic_ns = original_monotonic_ns  # type: ignore[assignment]
            os.sched_getaffinity = original_sched_getaffinity  # type: ignore[assignment]
            globals()["read_process_security"] = original_read_process_security
            os.close(settle_receipt_fd)
            os.close(settle_validation_fd)

    validation_lines = health["thermal"][
        "validation_jsonl_ascii"
    ].encode("ascii").splitlines(keepends=True)
    quantized_handles: Dict[str, Any] = {
        "validation_buffer": b"",
        "validation_csv_timestamp_ns": [],
        "validation_last_csv_timestamp_ns": None,
        "validation_last_monotonic_ns": None,
        "validation_monotonic_ns": [],
        "validation_sample_index": None,
    }
    quantized_first = parse_canonical_json_line(
        validation_lines[1], "synthetic quantized validation first"
    )
    quantized_first["monotonic_s"] = 1.2345674
    quantized_second = parse_canonical_json_line(
        validation_lines[2], "synthetic quantized validation second"
    )
    quantized_second["monotonic_s"] = 2.0000006
    if (
        consume_sampler_validation_bytes(
            quantized_handles,
            canonical_bytes(quantized_first) + b"\n",
            "synthetic quantized validation first",
        ) != "none"
        or consume_sampler_validation_bytes(
            quantized_handles,
            canonical_bytes(quantized_second) + b"\n",
            "synthetic quantized validation second",
        ) != "none"
        or quantized_handles["validation_monotonic_ns"]
        != [1_234_567_400, 2_000_000_600]
        or quantized_handles["validation_csv_timestamp_ns"]
        != [1_234_567_000, 2_000_001_000]
    ):
        raise AssertionError("sampler validation/CSV timestamp quantization differs")

    original_poll_sampler_supervision = globals()["poll_sampler_supervision"]
    try:
        paired_target = 1_234_567_000
        delayed_pair_handles: Dict[str, Any] = {
            "validation_csv_timestamp_ns": [],
            "validation_last_csv_timestamp_ns": None,
        }
        delayed_pair_polls = 0

        def delayed_pair_poll(_handles: Mapping[str, Any]) -> str:
            nonlocal delayed_pair_polls
            delayed_pair_polls += 1
            if delayed_pair_polls == 2:
                delayed_pair_handles["validation_csv_timestamp_ns"] = [
                    paired_target
                ]
                delayed_pair_handles["validation_last_csv_timestamp_ns"] = (
                    paired_target
                )
            return "none"

        globals()["poll_sampler_supervision"] = delayed_pair_poll
        if (
            wait_for_sampler_validation_timestamp_supervised(
                delayed_pair_handles, paired_target, time.monotonic() + 1.0
            ) != "none"
            or delayed_pair_polls < 3
        ):
            raise AssertionError("raw admission endpoint bypassed validation pairing")

        terminal_pair_handles: Dict[str, Any] = {
            "validation_csv_timestamp_ns": [],
            "validation_last_csv_timestamp_ns": None,
        }

        def terminal_pair_poll(_handles: Mapping[str, Any]) -> str:
            terminal_pair_handles["validation_csv_timestamp_ns"] = [
                paired_target
            ]
            terminal_pair_handles["validation_last_csv_timestamp_ns"] = (
                paired_target
            )
            return "sampler-terminal-decision:thermal_abort"

        globals()["poll_sampler_supervision"] = terminal_pair_poll
        if wait_for_sampler_validation_timestamp_supervised(
            terminal_pair_handles, paired_target, time.monotonic() + 1.0
        ) != "sampler-terminal-decision:thermal_abort":
            raise AssertionError("terminal validation decision lost to paired endpoint")

        confirmed_pair_handles: Dict[str, Any] = {
            "validation_csv_timestamp_ns": [paired_target],
            "validation_last_csv_timestamp_ns": paired_target,
        }
        confirmation_polls = 0

        def confirmation_poll(_handles: Mapping[str, Any]) -> str:
            nonlocal confirmation_polls
            confirmation_polls += 1
            return (
                "sampler-terminal-intent"
                if confirmation_polls == 2 else "none"
            )

        globals()["poll_sampler_supervision"] = confirmation_poll
        if wait_for_sampler_validation_timestamp_supervised(
            confirmed_pair_handles, paired_target, time.monotonic() + 1.0
        ) != "sampler-terminal-intent":
            raise AssertionError("post-pair sampler intent did not win admission")

        skipped_pair_handles: Dict[str, Any] = {
            "validation_csv_timestamp_ns": [paired_target + 1000],
            "validation_last_csv_timestamp_ns": paired_target + 1000,
        }
        globals()["poll_sampler_supervision"] = lambda _handles: "none"
        expect_validation_error(
            "skipped sampler validation endpoint",
            lambda: wait_for_sampler_validation_timestamp_supervised(
                skipped_pair_handles, paired_target, time.monotonic() + 1.0
            ),
        )
        missing_pair_handles: Dict[str, Any] = {
            "validation_csv_timestamp_ns": [],
            "validation_last_csv_timestamp_ns": None,
        }
        expect_validation_error(
            "sampler validation endpoint deadline",
            lambda: wait_for_sampler_validation_timestamp_supervised(
                missing_pair_handles, paired_target, time.monotonic()
            ),
        )
    finally:
        globals()["poll_sampler_supervision"] = original_poll_sampler_supervision

    cutoff_receipt = json.loads(canonical_bytes(health).decode("ascii"))
    cutoff_receipt.update({
        "child_reap_monotonic_ns": None,
        "child_start_monotonic_ns": None,
        "controller_singleton_affinity_end": [],
        "sampler": {},
        "sampler_receipt_sha256": None,
        "sampler_terminal": {},
        "sampler_terminal_receipt_sha256": None,
        "sibling_ticks": [],
        "thermal": {},
    })
    cutoff_receipt = finalize_health_receipt(cutoff_receipt, [], [])
    cutoff_state = {
        "collection_failures": [],
        "ready": True,
        "receipt": cutoff_receipt,
        "sample_start_ns": health["sampler"][
            "window_start_monotonic_ns"
        ],
        "sampler_handles": {"synthetic": True},
        "sampler_supervision_event": "none",
        "sibling_start": health["sibling_ticks"][0]["start"],
        "violations": [],
    }
    cutoff_args = argparse.Namespace()
    cutoff_make_calls = 0
    cutoff_thermal_calls = 0
    cutoff_clock = 98.5
    cutoff_deadline = 100.0
    original_make_sampler_attestation = globals()["make_sampler_attestation"]
    original_tolerant_thermal_window = globals()["tolerant_thermal_window"]
    original_wait_sampler_timestamp = globals()[
        "wait_for_sampler_timestamp_supervised"
    ]
    original_wait_sampler_validation = globals()[
        "wait_for_sampler_validation_timestamp_supervised"
    ]
    original_cutoff_poll_sampler = globals()["poll_sampler_supervision"]
    original_close_sampler_handles = globals()[
        "close_sampler_admission_handles"
    ]
    original_read_tick = globals()["read_cpu_tick_receipt"]
    original_sched_getaffinity = os.sched_getaffinity
    original_monotonic = time.monotonic

    def cutoff_make_sampler(
        _args: argparse.Namespace, start_ns: int, end_ns: int,
        _modules: Any, deadline: Optional[float] = None, *,
        allow_stale_evidence: bool = False,
        stale_stream_capture: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        nonlocal cutoff_make_calls
        cutoff_make_calls += 1
        if (
            start_ns != health["sampler"]["window_start_monotonic_ns"]
            or end_ns != start_ns or deadline != 99.5
            or not allow_stale_evidence
            or stale_stream_capture is None
        ):
            raise AssertionError("coverage cutoff stale-attestation call differs")
        stale_stream_capture["cutoff-snapshot"] = True
        return json.loads(canonical_bytes(health["sampler"]).decode("ascii"))

    def cutoff_thermal_window(
        _args: argparse.Namespace, sampler: Mapping[str, Any],
        start_ns: int, end_ns: int, *,
        allow_stale_suffixes: bool = False,
        deadline: Optional[float] = None,
        stale_stream_capture: Optional[Mapping[str, Any]] = None,
    ) -> Tuple[Dict[str, Any], List[str], List[str]]:
        nonlocal cutoff_thermal_calls
        cutoff_thermal_calls += 1
        if (
            not canonical_equal(sampler, health["sampler"])
            or start_ns != health["sampler"]["window_start_monotonic_ns"]
            or end_ns != health["sampler"]["window_end_monotonic_ns"]
            or not allow_stale_suffixes or deadline != 99.5
            or stale_stream_capture != {"cutoff-snapshot": True}
        ):
            raise AssertionError("coverage cutoff thermal call differs")
        return (
            json.loads(canonical_bytes(health["thermal"]).decode("ascii")),
            [], [],
        )

    def forbidden_cutoff_wait(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("coverage cutoff attempted an endpoint wait")

    try:
        globals()["make_sampler_attestation"] = cutoff_make_sampler
        globals()["tolerant_thermal_window"] = cutoff_thermal_window
        globals()["wait_for_sampler_timestamp_supervised"] = (
            forbidden_cutoff_wait
        )
        globals()["wait_for_sampler_validation_timestamp_supervised"] = (
            forbidden_cutoff_wait
        )
        globals()["poll_sampler_supervision"] = lambda _handles: "none"
        globals()["close_sampler_admission_handles"] = lambda _handles: None
        globals()["read_cpu_tick_receipt"] = (
            lambda _cpu: health["sibling_ticks"][0]["end"]
        )
        os.sched_getaffinity = (  # type: ignore[assignment]
            lambda _pid: {EXPECTED_CONTROLLER_CPU}
        )
        time.monotonic = lambda: cutoff_clock  # type: ignore[assignment]
        cutoff_health, cutoff_errors = finish_health(
            cutoff_args, EXPECTED_TARGET_CPU,
            binary["child_start_monotonic_ns"],
            binary["child_reap_monotonic_ns"],
            cutoff_state, cutoff_deadline, object(),
        )
    finally:
        time.monotonic = original_monotonic  # type: ignore[assignment]
        os.sched_getaffinity = original_sched_getaffinity  # type: ignore[assignment]
        globals()["read_cpu_tick_receipt"] = original_read_tick
        globals()["close_sampler_admission_handles"] = (
            original_close_sampler_handles
        )
        globals()["wait_for_sampler_validation_timestamp_supervised"] = (
            original_wait_sampler_validation
        )
        globals()["wait_for_sampler_timestamp_supervised"] = (
            original_wait_sampler_timestamp
        )
        globals()["poll_sampler_supervision"] = original_cutoff_poll_sampler
        globals()["tolerant_thermal_window"] = original_tolerant_thermal_window
        globals()["make_sampler_attestation"] = original_make_sampler_attestation
    cutoff_marker = (
        "sampler endpoint: sampler-monitor-invalid:"
        "coverage-recovery-cutoff"
    )
    if (
        cutoff_make_calls != 1 or cutoff_thermal_calls != 1
        or cutoff_health["collection_failures"] != [cutoff_marker]
        or cutoff_errors != ["health collection: " + cutoff_marker]
        or not canonical_equal(cutoff_health["sampler"], health["sampler"])
        or not canonical_equal(cutoff_health["thermal"], health["thermal"])
    ):
        raise AssertionError("coverage recovery cutoff branch differs")
    validate_health_receipt(cutoff_health, binary, require_complete=False)

    with tempfile.TemporaryDirectory() as directory:
        replaced_csv = Path(directory) / "sampler.csv"
        replaced_csv.write_bytes(
            health["thermal"]["window_csv_ascii"].encode("ascii")
        )
        os.chmod(replaced_csv, 0o600)
        replacement_info = os.stat(replaced_csv, follow_symlinks=False)
        stale_sampler = json.loads(
            canonical_bytes(health["sampler"]).decode("ascii")
        )
        stale_sampler["csv_device"] = replacement_info.st_dev
        stale_sampler["csv_inode"] = replacement_info.st_ino + 1
        thermal_args = argparse.Namespace(
            sampler_csv=replaced_csv,
            sampler_pid=stale_sampler["pid"],
            sampler_script=Path(stale_sampler["script_path"]),
        )
        replacement_thermal, replacement_collection, _ = (
            tolerant_thermal_window(
                thermal_args, stale_sampler,
                stale_sampler["window_start_monotonic_ns"],
                stale_sampler["window_end_monotonic_ns"],
            )
        )
        if replacement_thermal or not any(
            "thermal CSV identity changed" in message
            for message in replacement_collection
        ):
            raise AssertionError(
                "replacement sampler CSV did not fail as partial evidence"
            )

    if (
        reap_reached_controller_deadline(110_000_000_000, 0) is not True
        or reap_reached_controller_deadline(109_999_999_999, 0) is not False
        or monotonic_deadline_reached(
            119_999_999_999, 0, OUTER_DEADLINE_SECONDS
        ) is not False
        or monotonic_deadline_reached(
            120_000_000_000, 0, OUTER_DEADLINE_SECONDS
        ) is not True
    ):
        raise AssertionError("inclusive integer deadline law differs")
    expect_validation_error(
        "controller origin pair",
        lambda: run_binary(
            Path("/tmp/never-run-complement-worker"), EXPECTED_TARGET_CPU,
            "e" * 64, controller_started=0.0,
        ),
    )
    expect_validation_error(
        "controller origin value mismatch",
        lambda: run_binary(
            Path("/tmp/never-run-complement-worker"), EXPECTED_TARGET_CPU,
            "e" * 64, controller_started=1.0,
            controller_started_ns=2_000_000_000,
        ),
    )

    with tempfile.TemporaryDirectory() as directory:
        fake_binary = Path(directory) / "synthetic-worker"
        fake_binary_bytes = b"synthetic worker image; never executed\n"
        fake_binary.write_bytes(fake_binary_bytes)
        os.chmod(fake_binary, 0o555)
        health_args = argparse.Namespace(
            expected_binary_uid=os.getuid(),
            expected_controller_gid=os.getgid(),
        )
        origin_ns = time.monotonic_ns()
        prepared_receipt = empty_health_receipt(health_args)
        prepared_admission = json.loads(canonical_bytes(
            health["admission_sibling_ticks"]
        ).decode("ascii"))
        prepared_admission["start"]["read_monotonic_ns"] = origin_ns
        prepared_admission["end"]["read_monotonic_ns"] = origin_ns
        prepared_admission["interval_start_monotonic_ns"] = origin_ns
        prepared_admission["interval_end_monotonic_ns"] = origin_ns
        prepared_receipt["admission_sibling_ticks"] = prepared_admission
        for field in (
            "target_core", "controller_core", "sampler_core",
            "target_threads", "controller_initial_affinity",
        ):
            prepared_receipt[field] = list(health[field])
        prepared_receipt = finalize_health_receipt(prepared_receipt, [], [])
        prepared_state = {
            "collection_failures": [], "ready": True,
            "receipt": prepared_receipt,
            "sampler_handles": {"synthetic": True}, "violations": [],
        }
        observed_finish_intervals: List[
            Tuple[Optional[int], Optional[int]]
        ] = []
        original_popen = subprocess.Popen
        original_finish_health = globals()["finish_health"]
        original_read_tick = globals()["read_cpu_tick_receipt"]
        original_process_start_ticks = globals()["process_start_ticks"]
        original_poll_sampler_supervision = globals()[
            "poll_sampler_supervision"
        ]
        original_terminate = globals()["terminate_process_group"]

        def injected_finish_health(
            _args: argparse.Namespace, _cpu: int,
            child_start_ns: Optional[int], child_reap_ns: Optional[int],
            state: Mapping[str, Any], _deadline: float, _modules: Any,
        ) -> Tuple[Dict[str, Any], List[str]]:
            observed_finish_intervals.append((child_start_ns, child_reap_ns))
            receipt = dict(state["receipt"])
            receipt["child_start_monotonic_ns"] = child_start_ns
            receipt["child_reap_monotonic_ns"] = child_reap_ns
            receipt = finalize_health_receipt(
                receipt, ["injected launch/attestation failure"], []
            )
            return receipt, ["health collection: injected failure"]

        def injected_finish_health_exception(
            _args: argparse.Namespace, _cpu: int,
            _child_start_ns: Optional[int], _child_reap_ns: Optional[int],
            _state: Mapping[str, Any], _deadline: float, _modules: Any,
        ) -> Tuple[Dict[str, Any], List[str]]:
            raise OSError("injected terminal health failure after reap")

        def injected_tick(_cpu: int) -> Dict[str, Any]:
            return {"synthetic": True}

        def injected_popen_failure(*_args: Any, **_kwargs: Any) -> Any:
            raise OSError("injected Popen failure")

        class SyntheticProcess:
            pid = 987654321
            returncode: Optional[int] = None

            def wait(self, timeout: float) -> int:
                del timeout
                self.returncode = -9
                return self.returncode

        class SyntheticUnreapedProcess(SyntheticProcess):
            def wait(self, timeout: float) -> int:
                raise subprocess.TimeoutExpired("synthetic-worker", timeout)

        try:
            globals()["finish_health"] = injected_finish_health
            globals()["read_cpu_tick_receipt"] = injected_tick
            globals()["poll_sampler_supervision"] = lambda _handles: "none"
            subprocess.Popen = injected_popen_failure  # type: ignore[assignment]
            launch_failure = run_binary(
                fake_binary, EXPECTED_TARGET_CPU,
                sha256_bytes(fake_binary_bytes),
                controller_started=origin_ns / 1_000_000_000.0,
                controller_started_ns=origin_ns,
                health_args=health_args, health_modules=object(),
                prepared_health_state=prepared_state,
            )
            if (
                launch_failure.binary["process_started"]
                or launch_failure.binary["child_start_monotonic_ns"] is not None
                or launch_failure.binary["child_reap_monotonic_ns"] is not None
                or observed_finish_intervals[-1] != (None, None)
            ):
                raise AssertionError("Popen failure invented a child interval")
            validate_binary_receipt(
                launch_failure.binary, require_complete=False
            )
            validate_health_receipt(
                launch_failure.health, launch_failure.binary,
                require_complete=False,
            )

            observed_finish_intervals.clear()
            synthetic_process = SyntheticProcess()
            subprocess.Popen = (  # type: ignore[assignment]
                lambda *_args, **_kwargs: synthetic_process
            )
            globals()["process_start_ticks"] = (
                lambda _pid: (_ for _ in ()).throw(
                    OSError("injected live-attestation failure")
                )
            )
            globals()["terminate_process_group"] = lambda _process: None
            post_launch_failure = run_binary(
                fake_binary, EXPECTED_TARGET_CPU,
                sha256_bytes(fake_binary_bytes),
                controller_started=origin_ns / 1_000_000_000.0,
                controller_started_ns=origin_ns,
                health_args=health_args, health_modules=object(),
                prepared_health_state=prepared_state,
            )
            if (
                not post_launch_failure.binary["process_started"]
                or post_launch_failure.binary["child_start_monotonic_ns"] is None
                or post_launch_failure.binary["child_reap_monotonic_ns"] is None
                or observed_finish_intervals[-1] != (
                    post_launch_failure.binary["child_start_monotonic_ns"],
                    post_launch_failure.binary["child_reap_monotonic_ns"],
                )
            ):
                raise AssertionError(
                    "post-launch failure lost its partial child interval"
                )
            validate_binary_receipt(
                post_launch_failure.binary, require_complete=False
            )
            validate_health_receipt(
                post_launch_failure.health, post_launch_failure.binary,
                require_complete=False,
            )

            health_finish_process = SyntheticProcess()
            subprocess.Popen = (  # type: ignore[assignment]
                lambda *_args, **_kwargs: health_finish_process
            )
            globals()["finish_health"] = injected_finish_health_exception
            health_finish_failure = run_binary(
                fake_binary, EXPECTED_TARGET_CPU,
                sha256_bytes(fake_binary_bytes),
                controller_started=origin_ns / 1_000_000_000.0,
                controller_started_ns=origin_ns,
                health_args=health_args, health_modules=object(),
                prepared_health_state=prepared_state,
            )
            if (
                not health_finish_failure.binary["process_started"]
                or health_finish_failure.returncode is None
                or health_finish_failure.binary[
                    "child_start_monotonic_ns"
                ] is None
                or health_finish_failure.binary[
                    "child_reap_monotonic_ns"
                ] is None
                or health_finish_failure.health[
                    "child_start_monotonic_ns"
                ] != health_finish_failure.binary[
                    "child_start_monotonic_ns"
                ]
                or health_finish_failure.health[
                    "child_reap_monotonic_ns"
                ] != health_finish_failure.binary[
                    "child_reap_monotonic_ns"
                ]
                or not any(
                    "injected terminal health failure after reap" in message
                    for message in health_finish_failure.health[
                        "collection_failures"
                    ]
                )
            ):
                raise AssertionError(
                    "terminal health failure lost confirmed child evidence"
                )
            validate_binary_receipt(
                health_finish_failure.binary, require_complete=False
            )
            validate_health_receipt(
                health_finish_failure.health,
                health_finish_failure.binary,
                require_complete=False,
            )
            validate_execution_reap_binding(
                health_finish_failure.binary,
                health_finish_failure.returncode,
                "synthetic terminal health failure",
            )

            failure_summary_document = json.loads(
                canonical_bytes(final_fixture["summary"]).decode("ascii")
            )
            failure_summary_document["binary"] = (
                health_finish_failure.binary
            )
            failure_summary_document["elapsed_seconds"] = (
                binary_execution_elapsed(health_finish_failure.binary)
            )
            failure_summary_document["failure"] = (
                "synthetic terminal health failure"
            )
            failure_summary_document["health"] = (
                health_finish_failure.health
            )
            failure_summary_document["process_exit_code"] = (
                health_finish_failure.returncode
            )
            failure_summary_document["summary_preimage_sha256"] = None
            make_summary_preimage(failure_summary_document)
            failure_summary_bytes = (
                canonical_bytes(failure_summary_document) + b"\n"
            )
            failure_summary = parse_summary_bytes(failure_summary_bytes)
            validate_retained_science(b"", b"", b"", failure_summary)

            failure_controller_document = json.loads(
                canonical_bytes(final_fixture["controller"]).decode("ascii")
            )
            failure_controller_document["binary"] = (
                health_finish_failure.binary
            )
            failure_controller_document["controller_started_monotonic_ns"] = (
                origin_ns
            )
            failure_controller_document["controller_observed_monotonic_ns"] = (
                health_finish_failure.binary[
                    "execution_finished_monotonic_ns"
                ]
            )
            failure_controller_document["controller_elapsed_seconds"] = (
                (
                    failure_controller_document[
                        "controller_observed_monotonic_ns"
                    ] - origin_ns
                ) / 1_000_000_000.0
            )
            failure_controller_document["failure"] = (
                "synthetic terminal health failure"
            )
            failure_controller_document["health_receipt_sha256"] = (
                health_finish_failure.health["receipt_sha256"]
            )
            failure_controller_document["sampler_after"] = {}
            failure_controller_document["summary"] = {
                "outcome": "invalid",
                "preimage_sha256": failure_summary[
                    "summary_preimage_sha256"
                ],
                "sha256": sha256_bytes(failure_summary_bytes),
            }
            failure_controller_document["artifacts"]["summary.json"][
                "bytes"
            ] = len(failure_summary_bytes)
            failure_controller_document["artifacts"]["summary.json"][
                "sha256"
            ] = sha256_bytes(failure_summary_bytes)
            failure_controller_document["artifacts"]["thermal.csv"][
                "bytes"
            ] = 0
            failure_controller_document["artifacts"]["thermal.csv"][
                "sha256"
            ] = sha256_bytes(b"")
            failure_controller_bytes = controller_document_bytes(
                failure_controller_document
            )
            failure_controller = parse_controller_document(
                failure_controller_bytes
            )
            validate_controller_bundle_binding(
                failure_controller, failure_summary, failure_summary_bytes,
                final_fixture["claim_receipt"],
                health_finish_failure.binary,
            )

            failure_complete_document = json.loads(
                canonical_bytes(final_fixture["complete"]).decode("ascii")
            )
            failure_complete_document["controller_sha256"] = sha256_bytes(
                failure_controller_bytes
            )
            failure_complete_document["complete_observed_monotonic_ns"] = (
                failure_controller_document[
                    "controller_observed_monotonic_ns"
                ] + 1
            )
            failure_complete_document["elapsed_seconds_before_commit"] = (
                (
                    failure_complete_document[
                        "complete_observed_monotonic_ns"
                    ] - origin_ns
                ) / 1_000_000_000.0
            )
            failure_complete = parse_complete_document(
                canonical_bytes(failure_complete_document) + b"\n"
            )
            validate_complete_bundle_binding(
                failure_complete, failure_controller,
                failure_controller_bytes, final_fixture["claim_receipt"],
            )

            observed_finish_intervals.clear()
            unreaped_process = SyntheticUnreapedProcess()
            subprocess.Popen = (  # type: ignore[assignment]
                lambda *_args, **_kwargs: unreaped_process
            )
            globals()["finish_health"] = injected_finish_health
            unreaped_failure = run_binary(
                fake_binary, EXPECTED_TARGET_CPU,
                sha256_bytes(fake_binary_bytes),
                controller_started=origin_ns / 1_000_000_000.0,
                controller_started_ns=origin_ns,
                health_args=health_args, health_modules=object(),
                prepared_health_state=prepared_state,
            )
            if (
                not unreaped_failure.binary["process_started"]
                or unreaped_failure.returncode is not None
                or unreaped_failure.binary["child_reap_monotonic_ns"] is not None
            ):
                raise AssertionError("unreaped-child injection did not hold")
            validate_binary_receipt(
                unreaped_failure.binary, require_complete=False
            )
            validate_health_receipt(
                unreaped_failure.health, unreaped_failure.binary,
                require_complete=False,
            )
            expect_validation_error(
                "unconfirmed worker reap publication",
                lambda: validate_execution_reap_binding(
                    unreaped_failure.binary, unreaped_failure.returncode,
                    "synthetic final publication",
                ),
            )
        finally:
            subprocess.Popen = original_popen  # type: ignore[assignment]
            globals()["finish_health"] = original_finish_health
            globals()["read_cpu_tick_receipt"] = original_read_tick
            globals()["process_start_ticks"] = original_process_start_ticks
            globals()["poll_sampler_supervision"] = (
                original_poll_sampler_supervision
            )
            globals()["terminate_process_group"] = original_terminate

    previous_receipts = EXPECTED_CELL_RECEIPTS
    EXPECTED_CELL_RECEIPTS = synthetic_receipts
    try:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            transcript_bundle = OutputBundle.reserve(
                root / "transcript-invalid", EXPECTED_TARGET_CPU, commit,
                "2" * 64, "3" * 64, root,
            )
            original_write = OutputBundle._write_exact
            stderr_fd = transcript_bundle.file_fds["stderr.txt"]

            def stderr_write_then_raise(fd: int, data: bytes) -> bytes:
                readback = original_write(fd, data)
                if fd == stderr_fd:
                    raise OSError("injected stderr post-fsync failure")
                return readback

            OutputBundle._write_exact = staticmethod(stderr_write_then_raise)
            try:
                transcript_failures = transcript_bundle.stage_evidence(
                    raw, b"post-terminal deadline\n"
                )
            finally:
                OutputBundle._write_exact = staticmethod(original_write)
            if (
                not transcript_failures
                or not transcript_bundle.raw_validation_config
                or not transcript_bundle.raw_validation_statistics
            ):
                raise AssertionError(
                    "complete raw was coupled to stderr/publication status"
                )
            transcript_bundle.update_fallback_summary_components({
                "binary": binary, "elapsed_seconds": 106.0,
                "failure": "post-terminal deadline",
                "git_after": {}, "git_before": {}, "health": health,
                "health_module_loader": {}, "outcome": "invalid",
                "process_exit_code": 1, "source_manifest_after": {},
                "source_manifest_before": {}, "target_cpu": EXPECTED_TARGET_CPU,
            })
            transcript_invalid = parse_summary_bytes(
                transcript_bundle.emergency_summary(
                    "post-terminal deadline", transcript_failures
                )
            )
            if (
                transcript_invalid["outcome"] != "invalid"
                or transcript_invalid["raw_complete"] is not True
                or not transcript_invalid["config"]
                or not transcript_invalid["statistics"]
            ):
                raise AssertionError(
                    "complete transcript invalid evidence was not preserved"
                )
            validate_retained_science(
                transcript_bundle.staged_raw,
                transcript_bundle.staged_stderr,
                health["thermal"]["window_csv_ascii"].encode("ascii"),
                transcript_invalid,
            )
            transcript_bundle._cleanup_incomplete(True)
    finally:
        EXPECTED_CELL_RECEIPTS = previous_receipts

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        directory_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
        try:
            complete_fd = reserve_output_file(directory_fd, "COMPLETE", 0o000)
            complete_value = {
                "campaign": CAMPAIGN,
                "claim_sha256": "1" * 64,
                "complete_observed_monotonic_ns": 1_000_000_000,
                "controller_outcome": "invalid",
                "controller_sha256": "2" * 64,
                "elapsed_seconds_before_commit": 1.0,
                "schema": COMPLETE_SCHEMA,
                "status": "complete",
            }
            complete_bytes = canonical_bytes(complete_value) + b"\n"
            OutputBundle._write_exact(complete_fd, complete_bytes)
            parse_complete_document(complete_bytes)
            finalize_complete_marker(
                complete_fd, directory_fd, complete_bytes, time.monotonic_ns()
            )
            if stat.S_IMODE(os.fstat(complete_fd).st_mode) != 0o400:
                raise AssertionError("scratch COMPLETE did not become visible")
            os.close(complete_fd)

            failed_fd = reserve_output_file(directory_fd, "FAILED", 0o000)
            OutputBundle._write_exact(failed_fd, complete_bytes)
            original_read = OutputBundle._read_current

            def injected_read_failure(_fd: int, _limit: int) -> bytes:
                raise RuntimeError("injected precommit read failure")

            OutputBundle._read_current = staticmethod(injected_read_failure)
            try:
                try:
                    finalize_complete_marker(
                        failed_fd, directory_fd, complete_bytes,
                        time.monotonic_ns(),
                    )
                except (RuntimeError, ValidationError):
                    pass
                else:
                    raise AssertionError("injected COMPLETE failure committed")
            finally:
                OutputBundle._read_current = staticmethod(original_read)
            failed_info = os.fstat(failed_fd)
            if (
                stat.S_IMODE(failed_info.st_mode) == 0o400
                and failed_info.st_size > 0
            ):
                raise AssertionError("failed COMPLETE remained valid-looking")
            os.close(failed_fd)

            exposure = root / "exposure"
            exposure.mkdir(mode=0o700)
            exposure_fd = os.open(exposure, os.O_RDONLY | os.O_DIRECTORY)
            exposed_fd = reserve_output_file(
                exposure_fd, "COMPLETE", 0o000
            )
            OutputBundle._write_exact(exposed_fd, complete_bytes)
            real_fchmod = os.fchmod
            injected = []

            def expose_then_raise(fd: int, mode: int) -> None:
                real_fchmod(fd, mode)
                if fd == exposed_fd and mode == 0o400 and not injected:
                    injected.append(True)
                    raise RuntimeError("injected post-exposure failure")

            os.fchmod = expose_then_raise  # type: ignore[assignment]
            try:
                try:
                    finalize_complete_marker(
                        exposed_fd, exposure_fd, complete_bytes,
                        time.monotonic_ns(),
                    )
                except RuntimeError:
                    pass
                else:
                    raise AssertionError("post-exposure injection committed")
            finally:
                os.fchmod = real_fchmod  # type: ignore[assignment]
            exposed_info = os.fstat(exposed_fd)
            if (
                not injected
                or exposed_info.st_size != 0
                or stat.S_IMODE(exposed_info.st_mode) != 0o000
            ):
                raise AssertionError(
                    "post-exposure failure did not poison COMPLETE"
                )
            os.close(exposed_fd)
            os.close(exposure_fd)
        finally:
            os.close(directory_fd)

    for poison_case in ("identity-fchmod", "writable-reopen"):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "poison"
            output.mkdir(mode=0o700)
            directory_fd = os.open(output, os.O_RDONLY | os.O_DIRECTORY)
            summary_fd = reserve_output_file(
                directory_fd, "summary.json", 0o600
            )
            OutputBundle._write_exact(summary_fd, b"stale-positive\n")
            os.fchmod(summary_fd, 0o400)
            summary_info = os.fstat(summary_fd)
            stale_fd = os.dup(summary_fd)
            os.close(summary_fd)
            os.close(stale_fd)
            poison_bundle = OutputBundle(
                output, EXPECTED_TARGET_CPU, "1" * 40, "2" * 64,
                "3" * 64, root,
            )
            poison_bundle.directory_fd = directory_fd
            poison_bundle.file_identities["summary.json"] = (
                summary_info.st_dev, summary_info.st_ino,
                summary_info.st_nlink,
            )
            original_fchmod = os.fchmod
            original_open = os.open
            injected: List[bool] = []

            def poison_fchmod(fd: int, mode: int) -> None:
                if (
                    poison_case == "identity-fchmod" and mode == 0o600
                    and not injected
                ):
                    injected.append(True)
                    raise OSError("injected identity chmod failure")
                original_fchmod(fd, mode)

            def poison_open(
                path: Any, flags: int, mode: int = 0o777,
                *, dir_fd: Optional[int] = None,
            ) -> int:
                if (
                    poison_case == "writable-reopen"
                    and path == "summary.json"
                    and flags & os.O_ACCMODE == os.O_RDWR
                    and not injected
                ):
                    injected.append(True)
                    raise OSError("injected writable reopen failure")
                return original_open(path, flags, mode, dir_fd=dir_fd)

            os.fchmod = poison_fchmod  # type: ignore[assignment]
            os.open = poison_open  # type: ignore[assignment]
            try:
                poison_noncommitted_publication(
                    poison_bundle, {"summary.json": stale_fd},
                    ["injected abort"], SignalGuard(), None,
                )
            finally:
                os.fchmod = original_fchmod  # type: ignore[assignment]
                os.open = original_open  # type: ignore[assignment]
            if not injected:
                raise AssertionError(f"{poison_case} injection did not fire")
            try:
                poisoned_info = os.stat(
                    "summary.json", dir_fd=directory_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                poisoned_info = None
            if (
                poisoned_info is not None
                and stat.S_IMODE(poisoned_info.st_mode) != 0o000
            ):
                raise AssertionError(
                    f"{poison_case} left a readable stale summary"
                )
            os.fchmod(directory_fd, 0o700)
            if poisoned_info is not None:
                os.unlink("summary.json", dir_fd=directory_fd)
            os.close(directory_fd)

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        fallback_bundle = OutputBundle.reserve(
            root / "fallback", EXPECTED_TARGET_CPU, "1" * 40,
            "2" * 64, "3" * 64, root,
        )
        try:
            fallback_bundle.update_fallback_summary_components({
                "binary": partial_binary, "config": {},
                "elapsed_seconds": 0.0,
                "failure": "immediate post-reserve failure",
                "git_after": {}, "git_before": {},
                "health": partial_health, "health_module_loader": {},
                "outcome": "invalid", "process_exit_code": None,
                "source_manifest_after": {},
                "source_manifest_before": {}, "statistics": {},
                "target_cpu": EXPECTED_TARGET_CPU,
            })
            fallback_bundle.stage_evidence(b"", b"")
            immediate = parse_summary_bytes(
                fallback_bundle.emergency_summary(
                    "immediate post-reserve injection", []
                )
            )
            if not canonical_equal(immediate["health"], partial_health):
                raise AssertionError("post-reserve fallback lost partial health")
            fallback_snapshot = dict(
                fallback_bundle.fallback_summary_components
            )
            try:
                fallback_bundle.update_fallback_summary_components({
                    "binary": binary, "health": object(),
                })
            except TypeError:
                pass
            else:
                raise AssertionError("fallback encoding injection did not fail")
            if fallback_bundle.fallback_summary_components != fallback_snapshot:
                raise AssertionError("fallback update was not transactional")
            coherent_after_fault = parse_summary_bytes(
                fallback_bundle.emergency_summary(
                    "transactional fallback encoding injection", []
                )
            )
            if (
                not canonical_equal(
                    coherent_after_fault["binary"], partial_binary
                )
                or not canonical_equal(
                    coherent_after_fault["health"], partial_health
                )
            ):
                raise AssertionError(
                    "failed fallback update corrupted its prior snapshot"
                )
            fallback_bundle.update_fallback_summary_components({
                "binary": binary, "health": invalid_health,
                "elapsed_seconds": 1.0, "process_exit_code": 0,
            })
            returned = parse_summary_bytes(
                fallback_bundle.emergency_summary(
                    "immediate post-run-return injection", []
                )
            )
            if not canonical_equal(returned["health"], invalid_health):
                raise AssertionError("post-return fallback lost current health")
        finally:
            fallback_bundle._cleanup_incomplete(True)

    cutoff_guard = SignalGuard()
    cutoff_guard.entered = True
    cutoff_guard.seal_blocked = True

    def injected_cutoff_drain() -> bool:
        cutoff_guard.first_signal = "SIGTERM"
        cutoff_guard.observed_signals = ["SIGTERM"]
        return True

    cutoff_guard._collect_pending_unchecked = (  # type: ignore[method-assign]
        injected_cutoff_drain
    )
    if not cutoff_guard.commit_logical_seal():
        raise AssertionError("cutoff drain did not report its injected signal")
    cutoff_bundle = OutputBundle(
        Path("/tmp/cutoff"), EXPECTED_TARGET_CPU, "1" * 40,
        "2" * 64, "3" * 64, Path("/tmp"),
    )
    cutoff_summary = parse_summary_bytes(
        cutoff_bundle._invalid_summary(
            parsed_summary, ["signal cutoff rewrite"], cutoff_guard
        )
    )
    if (
        cutoff_summary["outcome"] != "invalid"
        or cutoff_summary["signal_names"] != ["SIGTERM"]
    ):
        raise AssertionError("cutoff signal did not rewrite invalid evidence")
    cutoff_guard._handler(signal.SIGINT, None)
    if cutoff_guard.observed_signals != ["SIGTERM"]:
        raise AssertionError("post-cutoff signal changed the frozen decision")

    verifier_tokens = [
        "--verify-retained", "--expected-run-argv-sha256",
        final_fixture["run_argv_sha256"],
    ]
    parsed_verifier = parse_cli_tokens(verifier_tokens, live_process=False)
    if (
        not parsed_verifier.verify_retained
        or parsed_verifier.expected_run_argv_sha256
        != final_fixture["run_argv_sha256"]
    ):
        raise AssertionError("retained verifier argument parsing differs")
    for label, tokens in (
        ("help exit collision", ["--help"]),
        ("missing verifier authority", ["--verify-retained"]),
        (
            "equals-form verifier authority",
            [
                "--verify-retained",
                "--expected-run-argv-sha256="
                + final_fixture["run_argv_sha256"],
            ],
        ),
        ("abbreviated verifier mode", ["--verify-ret"]),
        (
            "extra verifier argument",
            verifier_tokens + ["--expected-run-argv-sha256", "0" * 64],
        ),
    ):
        expect_validation_error(
            label,
            lambda tokens=tokens: parse_cli_tokens(tokens, live_process=False),
        )

    original_environment = dict(os.environ)
    try:
        for label, environment in (
            (
                "extra authority environment",
                {**SEALED_LAUNCH_ENVIRONMENT, "LD_PRELOAD": "/tmp/x"},
            ),
            (
                "missing authority environment",
                {
                    name: value for name, value in SEALED_LAUNCH_ENVIRONMENT.items()
                    if name != "TZ"
                },
            ),
            (
                "wrong authority environment",
                {**SEALED_LAUNCH_ENVIRONMENT, "TZ": "Etc/UTC"},
            ),
        ):
            os.environ.clear()
            os.environ.update(environment)
            expect_validation_error(
                label,
                lambda: parse_cli_tokens(verifier_tokens, live_process=True),
            )
    finally:
        os.environ.clear()
        os.environ.update(original_environment)

    original_argv0 = sys.argv[0]
    original_environment = dict(os.environ)
    try:
        os.environ.clear()
        os.environ.update(SEALED_LAUNCH_ENVIRONMENT)
        canonical_argv0 = str(Path(__file__).resolve(strict=True))
        sys.argv[0] = canonical_argv0
        if not parse_cli_tokens(
            verifier_tokens, live_process=True
        ).verify_retained:
            raise AssertionError("canonical live verifier argv was rejected")
        for label, argv0 in (
            (
                "relative authority controller argv",
                "bench/Wh2DirectSystematicComplementScreen.py",
            ),
            (
                "symlink-spelled authority controller argv",
                "/tmp/synthetic-controller-link.py",
            ),
        ):
            sys.argv[0] = argv0
            expect_validation_error(
                label,
                lambda: parse_cli_tokens(verifier_tokens, live_process=True),
            )
    finally:
        sys.argv[0] = original_argv0
        os.environ.clear()
        os.environ.update(original_environment)

    with tempfile.TemporaryDirectory() as directory:
        retained_root = Path(directory) / "retained"
        retained_root.mkdir(mode=0o700)
        retained_directory_fd = os.open(
            retained_root, os.O_RDONLY | os.O_DIRECTORY
        )
        retained_fds: List[int] = []
        try:
            for index, name in enumerate(FINAL_OUTPUT_NAMES):
                fd = reserve_output_file(retained_directory_fd, name, 0o400)
                retained_fds.append(fd)
                payload = (name + ":" + str(index) + "\n").encode("ascii")
                OutputBundle._write_exact(fd, payload)
                os.fchmod(fd, 0o400)
            retained_roster(retained_directory_fd)
            read_flags = (
                os.O_RDONLY | os.O_NOFOLLOW | os.O_NOATIME
                | os.O_NONBLOCK | getattr(os, "O_CLOEXEC", 0)
            )
            for index, name in enumerate(FINAL_OUTPUT_NAMES):
                read_fd = os.open(
                    name, read_flags, dir_fd=retained_directory_fd
                )
                try:
                    payload, info = retained_read_exact(
                        read_fd, 4096, "scratch retained " + name
                    )
                    expected = (name + ":" + str(index) + "\n").encode(
                        "ascii"
                    )
                    if payload != expected or info.st_uid != os.geteuid():
                        raise AssertionError("scratch retained read differs")
                finally:
                    os.close(read_fd)
            surviving_fd = os.open(
                "raw.jsonl", read_flags, dir_fd=retained_directory_fd
            )
            try:
                surviving_expected = held_artifact(
                    retained_fds[0], 0o400
                )
                verify_surviving_named_artifact(
                    retained_directory_fd, surviving_fd, "raw.jsonl",
                    surviving_expected, 0o400,
                    surviving_expected["uid"], surviving_expected["gid"],
                )
                os.rename(
                    "raw.jsonl", "raw.pre-replacement",
                    src_dir_fd=retained_directory_fd,
                    dst_dir_fd=retained_directory_fd,
                )
                replacement_fd = reserve_output_file(
                    retained_directory_fd, "raw.jsonl", 0o400
                )
                try:
                    OutputBundle._write_exact(
                        replacement_fd, b"raw.jsonl:0\n"
                    )
                    os.fchmod(replacement_fd, 0o400)
                    expect_validation_error(
                        "surviving mirror named replacement",
                        lambda: verify_surviving_named_artifact(
                            retained_directory_fd, surviving_fd,
                            "raw.jsonl", surviving_expected, 0o400,
                            surviving_expected["uid"],
                            surviving_expected["gid"],
                        ),
                    )
                finally:
                    os.close(replacement_fd)
                    os.unlink("raw.jsonl", dir_fd=retained_directory_fd)
                    os.rename(
                        "raw.pre-replacement", "raw.jsonl",
                        src_dir_fd=retained_directory_fd,
                        dst_dir_fd=retained_directory_fd,
                    )
            finally:
                os.close(surviving_fd)
            extra_fd = reserve_output_file(
                retained_directory_fd, "unexpected", 0o400
            )
            retained_fds.append(extra_fd)
            expect_validation_error(
                "retained exact file roster",
                lambda: retained_roster(retained_directory_fd),
            )
            os.mkfifo("claim-fifo", 0o400, dir_fd=retained_directory_fd)
            fifo_fd = os.open(
                "claim-fifo", read_flags | os.O_NONBLOCK,
                dir_fd=retained_directory_fd,
            )
            try:
                expect_validation_error(
                    "retained FIFO bounded rejection",
                    lambda: retained_read_exact(
                        fifo_fd, 4096, "scratch retained FIFO"
                    ),
                )
            finally:
                os.close(fifo_fd)
        finally:
            for fd in retained_fds:
                try:
                    os.close(fd)
                except OSError:
                    pass
            os.close(retained_directory_fd)

    sleeper = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time; time.sleep(10)"],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, start_new_session=True,
    )
    captured = bounded_capture(sleeper, time.monotonic() + 0.05)
    if not captured.timed_out or sleeper.returncode is None:
        raise AssertionError("bounded scratch timeout did not kill/reap")

    with tempfile.TemporaryDirectory() as directory:
        heartbeat_root = Path(directory)
        heartbeat_receipt = heartbeat_root / "receipt.json"
        heartbeat_validation = heartbeat_root / "validation.jsonl"
        heartbeat_receipt.write_bytes(b"")
        heartbeat_validation.write_bytes(b"")
        os.chmod(heartbeat_receipt, 0o444)
        os.chmod(heartbeat_validation, 0o444)
        heartbeat_receipt_fd = os.open(
            heartbeat_receipt, nonblocking_read_flags("heartbeat receipt")
        )
        heartbeat_validation_fd = os.open(
            heartbeat_validation,
            nonblocking_read_flags("heartbeat validation"),
        )
        heartbeat_pidfd_writer: Optional[int] = None
        if hasattr(os, "pidfd_open"):
            heartbeat_pidfd = os.pidfd_open(os.getpid(), 0)
        else:
            # Python 3.8 predates os.pidfd_open.  A quiet pipe read end has the
            # same not-readable selector state needed by this scratch-only
            # heartbeat-boundary test; authority modes still require a real
            # Linux pidfd and fail closed in open_sampler_admission_handles.
            heartbeat_pidfd, heartbeat_pidfd_writer = os.pipe()
        heartbeat_base = 10_000_000_000
        heartbeat_handles: Dict[str, Any] = {
            "file_fds": {
                "receipt_file": heartbeat_receipt_fd,
                "validation_jsonl": heartbeat_validation_fd,
            },
            "pidfd": heartbeat_pidfd,
            "pidfd_exit_observed_ns": None,
            "receipt_admission_stat": stat_receipt(
                os.fstat(heartbeat_receipt_fd)
            ),
            "validation_buffer": b"",
            "validation_csv_timestamp_ns": [heartbeat_base],
            "validation_last_csv_timestamp_ns": heartbeat_base,
            "validation_last_monotonic_ns": heartbeat_base,
            "validation_monotonic_ns": [heartbeat_base],
            "validation_offset": 0,
            "validation_sample_index": 0,
        }
        original_monotonic_ns = time.monotonic_ns
        try:
            time.monotonic_ns = (  # type: ignore[assignment]
                lambda: heartbeat_base + SAMPLER_HEARTBEAT_MAX_GAP_NS
            )
            if poll_sampler_supervision(heartbeat_handles) != "none":
                raise AssertionError("exact heartbeat boundary was rejected")
            time.monotonic_ns = (  # type: ignore[assignment]
                lambda: heartbeat_base + SAMPLER_HEARTBEAT_MAX_GAP_NS + 1
            )
            if poll_sampler_supervision(heartbeat_handles) != (
                "sampler-monitor-invalid:validation-heartbeat-stale"
            ):
                raise AssertionError("stale sampler heartbeat was not detected")
        finally:
            time.monotonic_ns = original_monotonic_ns  # type: ignore[assignment]
            os.close(heartbeat_pidfd)
            if heartbeat_pidfd_writer is not None:
                os.close(heartbeat_pidfd_writer)
            os.close(heartbeat_validation_fd)
            os.close(heartbeat_receipt_fd)

    closed_pipe_process = subprocess.Popen(
        [
            sys.executable, "-I", "-B", "-c",
            "import os,time;os.close(1);os.close(2);time.sleep(10)",
        ],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, start_new_session=True,
    )
    original_poll_sampler_supervision = globals()["poll_sampler_supervision"]
    try:
        def closed_pipe_sampler_event(_handles: Mapping[str, Any]) -> str:
            return (
                "sampler-monitor-invalid:validation-heartbeat-stale"
                if closed_pipe_process.stdout is not None
                and closed_pipe_process.stderr is not None
                and closed_pipe_process.stdout.closed
                and closed_pipe_process.stderr.closed
                else "none"
            )

        globals()["poll_sampler_supervision"] = closed_pipe_sampler_event
        closed_pipe_capture = bounded_capture(
            closed_pipe_process, time.monotonic() + 1.0,
            sampler_monitor_handles={"synthetic": True},
        )
        if (
            closed_pipe_capture.sampler_event
            != "sampler-monitor-invalid:validation-heartbeat-stale"
            or closed_pipe_capture.timed_out
            or closed_pipe_process.returncode is None
        ):
            raise AssertionError(
                "closed worker pipes disabled sampler supervision"
            )
    finally:
        globals()["poll_sampler_supervision"] = (
            original_poll_sampler_supervision
        )
        if closed_pipe_process.returncode is None:
            terminate_process_group(closed_pipe_process)
            try:
                closed_pipe_process.wait(timeout=5.0)
            except BaseException:
                pass
    print("WH2 direct systematic complement analyzer selftest passed")
    return 0


class AuthorityArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        fail("argument validation: " + message)


def argument_parser() -> AuthorityArgumentParser:
    parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest", action="store_true")
    mode.add_argument("--preflight-seal", action="store_true")
    mode.add_argument("--run-once", action="store_true")
    mode.add_argument("--verify-retained", action="store_true")
    parser.add_argument("--binary")
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--controller-cpu", type=int)
    parser.add_argument("--sampler-pid", type=int)
    parser.add_argument("--sampler-cpu", type=int)
    parser.add_argument("--sampler-script", type=Path)
    parser.add_argument("--sampler-csv", type=Path)
    parser.add_argument("--sampler-pid-file", type=Path)
    parser.add_argument("--sampler-validation-jsonl", type=Path)
    parser.add_argument("--sampler-receipt", type=Path)
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--expected-binary-sha256")
    parser.add_argument("--expected-binary-uid", type=int)
    parser.add_argument("--expected-build-manifest-sha256")
    parser.add_argument("--expected-controller-sha256")
    parser.add_argument("--expected-controller-uid", type=int)
    parser.add_argument("--expected-controller-gid", type=int)
    parser.add_argument("--expected-git-sha256")
    parser.add_argument("--expected-python-sha256")
    parser.add_argument("--expected-sampler-process-start-ticks", type=int)
    parser.add_argument("--expected-sampler-script-sha256")
    parser.add_argument("--expected-sampler-csv-device", type=int)
    parser.add_argument("--expected-sampler-csv-inode", type=int)
    parser.add_argument("--expected-sampler-pid-file-device", type=int)
    parser.add_argument("--expected-sampler-pid-file-inode", type=int)
    parser.add_argument("--expected-sampler-validation-device", type=int)
    parser.add_argument("--expected-sampler-validation-inode", type=int)
    parser.add_argument("--expected-sampler-receipt-device", type=int)
    parser.add_argument("--expected-sampler-receipt-inode", type=int)
    parser.add_argument("--expected-sampler-cmdline-sha256")
    parser.add_argument("--expected-sampler-environ-sha256")
    parser.add_argument("--expected-sampler-executable-sha256")
    parser.add_argument("--expected-sampler-uid", type=int)
    parser.add_argument("--expected-sampler-gid", type=int)
    parser.add_argument("--expected-sampler-i2c-gid", type=int)
    parser.add_argument("--expected-source-manifest-sha256")
    parser.add_argument("--expected-run-argv-sha256")
    return parser


def validate_cli_token_shape(argv: Sequence[str]) -> None:
    tokens = list(argv)
    if tokens == ["--selftest"]:
        return
    if (
        len(tokens) == 1 + 2 * len(PREFLIGHT_SEAL_OPTION_ORDER)
        and tokens[0] == "--preflight-seal"
    ):
        for index, expected in enumerate(PREFLIGHT_SEAL_OPTION_ORDER):
            if tokens[1 + 2 * index] != expected:
                fail("preflight-seal option roster/order differs")
        return
    if (
        len(tokens) == 3
        and tokens[0] == "--verify-retained"
        and tokens[1] == "--expected-run-argv-sha256"
    ):
        return
    if len(tokens) == 1 + 2 * len(RUN_ONCE_OPTION_ORDER) and tokens[0] == "--run-once":
        for index, expected in enumerate(RUN_ONCE_OPTION_ORDER):
            if tokens[1 + 2 * index] != expected:
                fail("run-once option roster/order differs")
        return
    fail("argument vector does not match one exact public mode")


def require_authorized_process_identity(expected_uid: int, expected_gid: int) -> None:
    if not hasattr(os, "getresuid") or not hasattr(os, "getresgid"):
        fail("controller identity policy requires Linux getresuid/getresgid")
    if os.getresuid() != (expected_uid, expected_uid, expected_uid):
        fail("controller real/effective/saved UIDs differ from authorization")
    if os.getresgid() != (expected_gid, expected_gid, expected_gid):
        fail("controller real/effective/saved GIDs differ from authorization")
    if os.getgroups() != []:
        fail("controller supplementary group roster is not empty")
    security = read_process_security(os.getpid(), "controller live identity")
    validate_process_security(
        security, expected_uid, expected_gid, [], "controller live identity"
    )
    proc_info = os.stat(f"/proc/{os.getpid()}", follow_symlinks=False)
    if proc_info.st_uid != expected_uid or proc_info.st_gid != expected_gid:
        fail("controller /proc directory ownership differs")


def parse_cli_tokens(
    argv: Sequence[str], *, live_process: bool,
) -> argparse.Namespace:
    validate_cli_token_shape(argv)
    parser = argument_parser()
    args = parser.parse_args(argv)
    if live_process and (
        not sys.flags.isolated
        or not sys.dont_write_bytecode
        or sys.flags.optimize != 0
    ):
        parser.error(
            "this sealed runner requires unoptimized Python -I -B"
        )
    if args.selftest:
        # The selftest is a non-authoritative scratch-only mode registered by
        # CTest; launch-environment sealing applies to retained authority modes.
        return args
    if live_process and sys.argv[0] != str(Path(__file__).resolve(strict=True)):
        parser.error("authority mode controller argv[0] is not canonical")
    if live_process and not canonical_equal(
        dict(os.environ), SEALED_LAUNCH_ENVIRONMENT
    ):
        parser.error("authority mode startup environment differs")
    if args.preflight_seal:
        source_root = Path(__file__).resolve(strict=True).parents[1]
        token_values = {
            argv[1 + 2 * index]: argv[2 + 2 * index]
            for index in range(len(PREFLIGHT_SEAL_OPTION_ORDER))
        }
        paths = (
            args.build_dir, args.sampler_script, args.sampler_csv,
            args.sampler_pid_file, args.sampler_validation_jsonl,
            args.sampler_receipt,
        )
        if (
            not args.binary
            or not Path(args.binary).is_absolute()
            or os.path.normpath(args.binary) != args.binary
            or args.sampler_pid is None
            or not 1 <= args.sampler_pid <= MAX_INT63
            or token_values["--sampler-pid"] != str(args.sampler_pid)
            or (live_process and args.sampler_pid == os.getpid())
            or type(args.expected_source_commit) is not str
            or LOWER40.fullmatch(args.expected_source_commit) is None
            or any(
                path is None or not path.is_absolute()
                or os.path.normpath(str(path)) != str(path)
                for path in paths
            )
            or args.sampler_script != source_root / SAMPLER_SOURCE_PATH
            or len({args.binary, *(str(path) for path in paths)})
            != len(paths) + 1
            or any(
                value is not None
                for value in (
                    args.cpu, args.controller_cpu, args.sampler_cpu,
                    args.expected_binary_sha256, args.expected_binary_uid,
                    args.expected_build_manifest_sha256,
                    args.expected_controller_sha256,
                    args.expected_controller_uid, args.expected_controller_gid,
                    args.expected_git_sha256, args.expected_python_sha256,
                    args.expected_sampler_process_start_ticks,
                    args.expected_sampler_script_sha256,
                    args.expected_sampler_csv_device,
                    args.expected_sampler_csv_inode,
                    args.expected_sampler_pid_file_device,
                    args.expected_sampler_pid_file_inode,
                    args.expected_sampler_validation_device,
                    args.expected_sampler_validation_inode,
                    args.expected_sampler_receipt_device,
                    args.expected_sampler_receipt_inode,
                    args.expected_sampler_cmdline_sha256,
                    args.expected_sampler_environ_sha256,
                    args.expected_sampler_executable_sha256,
                    args.expected_sampler_uid, args.expected_sampler_gid,
                    args.expected_sampler_i2c_gid,
                    args.expected_source_manifest_sha256,
                    args.expected_run_argv_sha256,
                )
            )
        ):
            parser.error(
                "--preflight-seal requires one canonical binary/build/sampler "
                "locator roster and one lowercase source commit"
            )
        if live_process:
            require_authorized_process_identity(
                EXPECTED_CAMPAIGN_UID, EXPECTED_CAMPAIGN_GID
            )
        return args
    if args.verify_retained:
        if (
            type(args.expected_run_argv_sha256) is not str
            or LOWER64.fullmatch(args.expected_run_argv_sha256) is None
        ):
            parser.error("retained verification requires one run-argv SHA-256")
        return args
    hashes = (
        args.expected_binary_sha256,
        args.expected_build_manifest_sha256,
        args.expected_controller_sha256,
        args.expected_git_sha256,
        args.expected_python_sha256,
        args.expected_sampler_script_sha256,
        args.expected_sampler_cmdline_sha256,
        args.expected_sampler_environ_sha256,
        args.expected_sampler_executable_sha256,
        args.expected_source_manifest_sha256,
    )
    paths = (
        args.build_dir, args.sampler_script, args.sampler_csv,
        args.sampler_pid_file, args.sampler_validation_jsonl,
        args.sampler_receipt,
    )
    identities = (
        args.expected_binary_uid,
        args.expected_controller_uid,
        args.expected_controller_gid,
        args.expected_sampler_process_start_ticks,
        args.expected_sampler_csv_device,
        args.expected_sampler_csv_inode,
        args.expected_sampler_pid_file_device,
        args.expected_sampler_pid_file_inode,
        args.expected_sampler_validation_device,
        args.expected_sampler_validation_inode,
        args.expected_sampler_receipt_device,
        args.expected_sampler_receipt_inode,
        args.expected_sampler_uid,
        args.expected_sampler_gid,
        args.expected_sampler_i2c_gid,
    )
    if (
        not args.binary
        or not Path(args.binary).is_absolute()
        or os.path.normpath(args.binary) != args.binary
        or args.cpu != EXPECTED_TARGET_CPU
        or args.controller_cpu != EXPECTED_CONTROLLER_CPU
        or args.sampler_pid is None or args.sampler_pid <= 0
        or (live_process and args.sampler_pid == os.getpid())
        or args.sampler_cpu != EXPECTED_SAMPLER_CPU
        or args.expected_binary_uid != EXPECTED_CAMPAIGN_UID
        or args.expected_controller_uid != EXPECTED_CAMPAIGN_UID
        or args.expected_controller_gid != EXPECTED_CAMPAIGN_GID
        or args.expected_sampler_uid != EXPECTED_CAMPAIGN_UID
        or args.expected_sampler_gid != EXPECTED_CAMPAIGN_GID
        or args.expected_sampler_i2c_gid != EXPECTED_SAMPLER_I2C_GID
        or any(
            path is None or not path.is_absolute()
            or os.path.normpath(str(path)) != str(path)
            for path in paths
        )
        or len({str(Path(args.binary)), *(str(path) for path in paths)})
        != len(paths) + 1
        or type(args.expected_source_commit) is not str
        or LOWER40.fullmatch(args.expected_source_commit) is None
        or any(type(value) is not str or LOWER64.fullmatch(value) is None
               for value in hashes)
        or any(type(value) is not int or value < 0 for value in identities)
        or args.expected_sampler_process_start_ticks == 0
        or args.expected_sampler_csv_inode == 0
        or args.expected_sampler_pid_file_inode == 0
        or args.expected_sampler_validation_inode == 0
        or args.expected_sampler_receipt_inode == 0
        or args.expected_run_argv_sha256 is not None
    ):
        parser.error(
            "--run-once requires exact CPUs 120/121/122, canonical distinct "
            "binary/build/sampler paths, and every source/build/controller/"
            "binary/sampler identity seal"
        )
    args.controller_initial_affinity = []
    if live_process:
        require_authorized_process_identity(
            args.expected_controller_uid, args.expected_controller_gid
        )
    return args


def parse_recorded_run_argv(
    argv: Any, controller_pid: int,
) -> argparse.Namespace:
    if (
        type(argv) is not list or not argv
        or argv[0] != str(Path(__file__).resolve(strict=True))
    ):
        fail("recorded controller argv path differs")
    args = parse_cli_tokens(argv[1:], live_process=False)
    if not args.run_once or args.sampler_pid == controller_pid:
        fail("recorded run authority process identities differ")
    return args


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    return parse_cli_tokens(argv, live_process=True)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest_v2()
    if args.preflight_seal:
        return emit_preflight_seal(args)
    if args.verify_retained:
        return verify_retained(args.expected_run_argv_sha256)
    return run_once(args)


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Exception as exc:
        print(
            "direct systematic complement analyzer failed: "
            + exception_text(exc),
            file=sys.stderr,
        )
        sys.exit(1)
