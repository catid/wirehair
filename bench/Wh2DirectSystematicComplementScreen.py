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
MAX_THERMAL_WINDOW_BYTES = 2 * 1024 * 1024
MAX_SOURCE_FILE_BYTES = 1024 * 1024
MAX_PROC_VECTOR_BYTES = 64 * 1024
MAX_FAILURE_TEXT_CHARS = 64 * 1024
MAX_FAILED_GATES = 64
MAX_PUBLICATION_FAILURES = 64
EXPECTED_TARGET_CPU = 120
EXPECTED_CONTROLLER_CPU = 121
EXPECTED_SAMPLER_CPU = 122
EXPECTED_TARGET_CORE = (0, 56)
EXPECTED_TARGET_THREADS = (56, 120)
EXPECTED_CONTROLLER_CORE = (0, 57)
EXPECTED_SAMPLER_CORE = (0, 58)
SIBLING_NON_IDLE_TICK_CAP = 5
THERMAL_MAX_MILLIC = 85000
HEALTH_SCHEMA = "wirehair.wh2.direct-systematic-complement-health.v2"
HEALTH_LOADER_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-health-source-loader.v2"
)
CONTROLLER_SCHEMA = "wirehair.wh2.direct-systematic-complement-controller.v2"
CONTROLLER_PROVENANCE_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-controller-provenance.v2"
)
CLAIM_SCHEMA = "wirehair.wh2.direct-systematic-complement-claim.v2"
COMPLETE_SCHEMA = "wirehair.wh2.direct-systematic-complement-complete.v2"
VERIFY_RESULT_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-retained-verification.v2"
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
    "wirehair.wh2.direct-systematic-complement-sampler-attestation.v2"
)
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

SOURCE_PATHS = (
    "CMakeLists.txt",
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
FINAL_COMMIT_START_DEADLINE_SECONDS = 117.0
RUN_ONCE_OPTION_ORDER = (
    "--binary", "--build-dir", "--cpu", "--controller-cpu",
    "--sampler-pid", "--sampler-cpu", "--sampler-script", "--sampler-csv",
    "--expected-source-commit", "--expected-binary-sha256",
    "--expected-binary-uid", "--expected-build-manifest-sha256",
    "--expected-controller-sha256", "--expected-controller-uid",
    "--expected-git-sha256", "--expected-python-sha256",
    "--expected-sampler-process-start-ticks",
    "--expected-sampler-script-sha256",
    "--expected-sampler-csv-device", "--expected-sampler-csv-inode",
    "--expected-sampler-cmdline-sha256", "--expected-sampler-environ-sha256",
    "--expected-sampler-executable-sha256", "--expected-sampler-uid",
    "--expected-source-manifest-sha256",
)


class ValidationError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise ValidationError(message)


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
            or stat.S_IMODE(info.st_mode) & 0o022
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
    error: str = "none"


def bounded_capture(
    process: subprocess.Popen, deadline: float,
    signal_guard: Optional[SignalGuard] = None,
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
        while selector.get_map():
            now = time.monotonic()
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
            events = selector.select(max(0.0, min(0.25, wait_until - now)))
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
    "sampler", "sampler_core", "sampler_cpu", "sampler_receipt_sha256",
    "schema", "sibling_non_idle_tick_cap", "sibling_tick_policy",
    "sibling_ticks", "target_core", "target_cpu", "target_threads",
    "thermal", "thermal_max_millic", "terminal_status",
}
SAMPLER_RECEIPT_KEYS = {
    "cmdline_sha256", "cpu", "csv_device", "csv_inode", "csv_path",
    "csv_stat",
    "environ_sha256", "environment", "environment_sha256",
    "executable_device", "executable_inode", "executable_path",
    "executable_sha256", "executable_stat", "pid", "process_start_ticks", "process_uid",
    "schema", "script_device", "script_inode", "script_path",
    "script_sha256", "script_stat", "terminal_status", "window_end_monotonic_ns",
    "window_start_monotonic_ns",
}
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
            or entry["stat"]["mode"] & 0o022
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
    for name in ("csv_stat", "executable_stat", "script_stat"):
        validate_stat_receipt(sampler[name], "sampler " + name)
    if (
        sampler["csv_stat"]["device"] != sampler["csv_device"]
        or sampler["csv_stat"]["inode"] != sampler["csv_inode"]
        or sampler["csv_stat"]["mode"] != 0o600
        or sampler["csv_stat"]["nlink"] != 1
        or sampler["csv_stat"]["uid"] != sampler["process_uid"]
        or not 1 <= sampler["csv_stat"]["size"] <= MAX_SAMPLER_CSV_BYTES
        or sampler["script_stat"]["device"] != sampler["script_device"]
        or sampler["script_stat"]["inode"] != sampler["script_inode"]
        or sampler["script_stat"]["mode"] != 0o444
        or sampler["script_stat"]["nlink"] != 1
        or sampler["script_stat"]["uid"] != sampler["process_uid"]
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
    for field in (
        "script_sha256", "cmdline_sha256", "executable_sha256",
        "environ_sha256", "environment_sha256",
    ):
        lower_hash(sampler[field], f"sampler {field}")
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
    bound_csv = dict(bound_outer.pop("csv_stat"))
    current_csv = dict(current_outer.pop("csv_stat"))
    if not canonical_equal(bound_outer, current_outer):
        fail(f"{where} immutable sampler attestation changed")
    bound_size = bound_csv.pop("size")
    current_size = current_csv.pop("size")
    bound_mtime = bound_csv.pop("mtime_ns")
    current_mtime = current_csv.pop("mtime_ns")
    if (
        not canonical_equal(bound_csv, current_csv)
        or current_size < bound_size
        or current_mtime < bound_mtime
    ):
        fail(f"{where} sampler CSV identity regressed or changed")


def validate_thermal_receipt(
    thermal: Any, sampler: Mapping[str, Any],
    child_start: Optional[int], child_reap: Optional[int],
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
        if not (
            thermal["window_start_monotonic_ns"] <= child_start
            <= child_reap <= thermal["window_end_monotonic_ns"]
        ):
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

    sampler = health["sampler"]
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
    elif health["sampler_receipt_sha256"] is not None:
        fail("partial health has an unbound sampler digest")
    if thermal:
        if not sampler:
            fail("thermal receipt lacks sampler identity")
        validate_thermal_receipt(thermal, sampler, child_start, child_reap)
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
        and admission and len(siblings) == 1 and sampler and thermal
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
        "sampler_core": [],
        "sampler_cpu": EXPECTED_SAMPLER_CPU,
        "sampler_receipt_sha256": None,
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
        and bool(receipt["sampler"])
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
) -> Dict[str, Any]:
    if deadline is not None and time.monotonic() >= deadline:
        fail("sampler attestation reached its global deadline")
    health_api = health_modules.native
    for path, where in (
        (args.sampler_script, "sampler script"),
        (args.sampler_csv, "sampler CSV"),
    ):
        if path.resolve(strict=True) != path:
            fail(f"{where} path is not canonical")
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
    script_hash, script_info = hash_regular_path(
        args.sampler_script, "sampler script",
        max_size=MAX_SOURCE_FILE_BYTES, deadline=deadline,
    )
    exact_string(
        script_hash, args.expected_sampler_script_sha256,
        "sampler script SHA-256",
    )
    csv_data, csv_info = sampler_snapshot(args.sampler_csv)
    if (
        stat.S_IMODE(script_info.st_mode) != 0o444
        or script_info.st_uid != args.expected_sampler_uid
        or stat.S_IMODE(csv_info.st_mode) != 0o600
        or csv_info.st_uid != args.expected_sampler_uid
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
    cmdline = read_bounded_proc_vector(
        Path(f"/proc/{args.sampler_pid}/cmdline"), "sampler command line"
    )
    exact_string(
        sha256_bytes(cmdline), args.expected_sampler_cmdline_sha256,
        "sampler command-line SHA-256",
    )
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
    if process_start_ticks(args.sampler_pid) != actual_start_ticks:
        fail("sampler process identity changed during attestation")
    if deadline is not None and time.monotonic() >= deadline:
        fail("sampler attestation reached its global deadline")
    if start_ns > end_ns or not csv_data:
        fail("sampler attestation window differs")
    return {
        "cmdline_sha256": args.expected_sampler_cmdline_sha256,
        "cpu": EXPECTED_SAMPLER_CPU,
        "csv_device": csv_info.st_dev,
        "csv_inode": csv_info.st_ino,
        "csv_path": str(args.sampler_csv),
        "csv_stat": stat_receipt(csv_info),
        "environ_sha256": sha256_bytes(environ),
        "environment": environment,
        "environment_sha256": sha256_bytes(canonical_bytes(environment)),
        "executable_device": executable_info.st_dev,
        "executable_inode": executable_info.st_ino,
        "executable_path": str(executable),
        "executable_sha256": executable_hash,
        "executable_stat": stat_receipt(executable_info),
        "pid": args.sampler_pid,
        "process_start_ticks": actual_start_ticks,
        "process_uid": proc_info.st_uid,
        "schema": SAMPLER_SCHEMA,
        "script_device": script_info.st_dev,
        "script_inode": script_info.st_ino,
        "script_path": str(args.sampler_script),
        "script_sha256": script_hash,
        "script_stat": stat_receipt(script_info),
        "terminal_status": "ok",
        "window_end_monotonic_ns": end_ns,
        "window_start_monotonic_ns": start_ns,
    }


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
                mono_ns <= prior_ns or mono_ns - prior_ns > 5_000_000_000
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
    end_ns: int,
) -> Tuple[Dict[str, Any], List[str], List[str]]:
    collection: List[str] = []
    violations: List[str] = []
    try:
        data, info = sampler_snapshot(args.sampler_csv)
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
        summary = summarize_thermal_window_bytes(window_bytes, start_ns, end_ns)
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
        sample_start_ns = wait_for_sampler_timestamp(
            args.sampler_csv, min(deadline, time.monotonic() + 6.0),
            greater_than=max(initial_timestamps.values()),
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
    if type(start_ns) is int:
        try:
            sample_end_ns = wait_for_sampler_timestamp(
                args.sampler_csv, min(deadline, time.monotonic() + 6.0),
                at_or_after=max(endpoint_target, start_ns + 1),
            )
            sampler = make_sampler_attestation(
                args, start_ns, sample_end_ns, health_modules, deadline
            )
            sampler_digest = sha256_bytes(canonical_bytes(sampler))
            thermal, thermal_collection, thermal_violations = (
                tolerant_thermal_window(
                    args, sampler, start_ns, sample_end_ns
                )
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
            append_unique_failure(collection, "sampler endpoint: " + exception_text(exc))
    else:
        append_unique_failure(collection, "sampler start endpoint unavailable")
    if time.monotonic() >= deadline:
        append_unique_failure(
            collection,
            "health collection reached the strict 115-second deadline",
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
        capture = bounded_capture(process, child_deadline, signal_guard)
        capture_started = True
        result.raw = capture.stdout
        result.stderr = capture.stderr
        result.binary["timed_out"] = capture.timed_out
        result.binary["output_overflow"] = capture.output_overflow
        result.binary["capture_error"] = capture.error
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
            [git_fd_path, "-c", "core.fsmonitor=false", *arguments],
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
        if stat.S_IMODE(info.st_mode) != 0o444:
            fail(f"sealed build provenance {relative} must have mode 0444")
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
        "process_start_ticks": process_start_ticks(os.getpid()),
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
        require_authorized_process_uid(args.expected_controller_uid)
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
    "optimize", "pid",
    "process_start_ticks",
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
        if entry["bytes"] > 64 * 1024 * 1024:
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
    if (
        type(argv) is not list or not argv or len(argv) > 64
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
        "git_sha256", "output_path", "parent_device", "parent_inode",
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
    thermal_bytes = (
        thermal_receipt["window_csv_ascii"].encode("ascii")
        if thermal_receipt else b""
    )
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
        thermal_ascii = execution.health.get("thermal", {}).get(
            "window_csv_ascii", ""
        )
        if type(thermal_ascii) is not str:
            fail("health thermal window is not text")
        thermal_bytes = thermal_ascii.encode("ascii")
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
        or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != 0o555
    ):
        fail("screen binary preflight seal differs")
    return {"sha256": digest, "stat": stat_receipt(info)}


def run_once(args: argparse.Namespace) -> int:
    controller_started_ns = time.monotonic_ns()
    controller_started = controller_started_ns / 1_000_000_000.0
    if not canonical_equal(dict(os.environ), SEALED_LAUNCH_ENVIRONMENT):
        fail("controller startup environment differs")
    source_root = Path(__file__).resolve(strict=True).parents[1]
    binary_path = Path(args.binary)
    guarded_signals = SignalGuard()
    with guarded_signals:
        require_authorized_process_uid(args.expected_controller_uid)
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
        make_sampler_attestation(
            args, 0, 0, health_modules, controller_started + 9.0
        )
        health_state = prepare_health(
            args, EXPECTED_TARGET_CPU,
            controller_started + CONTROLLER_DEADLINE_SECONDS,
            health_modules,
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
    child_pid = binary.get("child_pid")
    sampler_identity_collision = (
        bool(health_sampler) and child_pid == run.sampler_pid
    )
    script_entry = source["entries"][SOURCE_PATHS.index(
        "bench/Wh2DirectSystematicComplementScreen.py"
    )]
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
        or run.expected_controller_uid != before["script_stat"]["uid"]
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
    # The external owner seal authorizes the exact image checked before exec.
    # Invalid evidence may truthfully preserve a later chown/replacement in
    # either post-exec receipt.  Positive receipts already require both later
    # stats to equal stat_before in validate_binary_receipt(require_complete).
    stat_before = binary["stat_before"]
    if (
        binary["sha256_before"] is not None and stat_before
        and stat_before["uid"] != run.expected_binary_uid
    ):
        fail("recorded pre-exec binary UID authority binding differs")
    samplers = [
        health_sampler, controller["sampler_after"],
    ]
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
            or sampler["environ_sha256"]
            != run.expected_sampler_environ_sha256
            or sampler["executable_sha256"]
            != run.expected_sampler_executable_sha256
            or sampler["process_uid"] != run.expected_sampler_uid
        ):
            fail("recorded sampler authority binding differs")
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
    thermal_receipt = health.get("thermal", {})
    expected_thermal = (
        thermal_receipt["window_csv_ascii"].encode("ascii")
        if thermal_receipt else b""
    )
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
    owner_gid = claim_receipt["stat"]["gid"]
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
    sampler = {
        "cmdline_sha256": "b" * 64,
        "cpu": EXPECTED_SAMPLER_CPU,
        "csv_device": 3,
        "csv_inode": 4,
        "csv_path": "/tmp/synthetic-sampler.csv",
        "csv_stat": {
            "device": 3, "gid": 1000, "inode": 4, "mode": 0o600,
            "mtime_ns": 5, "nlink": 1, "size": 456, "uid": 1000,
        },
        "environ_sha256": "e" * 64,
        "environment": dict(SEALED_LAUNCH_ENVIRONMENT),
        "environment_sha256": sha256_bytes(
            canonical_bytes(SEALED_LAUNCH_ENVIRONMENT)
        ),
        "executable_device": 5,
        "executable_inode": 6,
        "executable_path": "/usr/bin/python3",
        "executable_sha256": "c" * 64,
        "executable_stat": {
            "device": 5, "gid": 0, "inode": 6, "mode": 0o755,
            "mtime_ns": 7, "nlink": 2, "size": 789, "uid": 0,
        },
        "pid": 124,
        "process_start_ticks": 457,
        "process_uid": 1000,
        "schema": SAMPLER_SCHEMA,
        "script_device": 7,
        "script_inode": 8,
        "script_path": "/tmp/synthetic-sampler.py",
        "script_sha256": "d" * 64,
        "script_stat": {
            "device": 7, "gid": 1000, "inode": 8, "mode": 0o444,
            "mtime_ns": 9, "nlink": 1, "size": 321, "uid": 1000,
        },
        "terminal_status": "ok",
        "window_end_monotonic_ns": 2_000_000_000,
        "window_start_monotonic_ns": 1_000_000_000,
    }
    header = ",".join(THERMAL_HEADER) + "\n"
    rows = []
    for utc, monotonic_value in (
        ("2026-01-01T00:00:01Z", "1.000000000"),
        ("2026-01-01T00:00:02Z", "2.000000000"),
    ):
        row = [
            utc, monotonic_value, "50.0", "3000.0", "50.0",
            *("40.0" for _ in range(8)), "0", "0.1", "0.1", "0.1",
            "0", "0",
        ]
        rows.append(",".join(row) + "\n")
    window_ascii = header + "".join(rows)
    window_bytes = window_ascii.encode("ascii")
    thermal_summary = summarize_thermal_window_bytes(
        window_bytes, sampler["window_start_monotonic_ns"],
        sampler["window_end_monotonic_ns"],
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
        "sampler_core": list(EXPECTED_SAMPLER_CORE),
        "sampler_cpu": EXPECTED_SAMPLER_CPU,
        "sampler_receipt_sha256": sha256_bytes(canonical_bytes(sampler)),
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
        digest = script_sha if path == script_relative else sha256_bytes(
            ("synthetic-source:" + path).encode("ascii")
        )
        size = 321 if path == script_relative else index + 1
        source_stat = (
            script_stat if path == script_relative
            else receipt_stat(24, 100 + index, 0o444, size)
        )
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
        "process_start_ticks": 222,
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
        str(EXPECTED_CONTROLLER_CPU), str(health["sampler"]["pid"]),
        str(EXPECTED_SAMPLER_CPU), health["sampler"]["script_path"],
        health["sampler"]["csv_path"], commit, binary["expected_sha256"],
        str(binary["stat_before"]["uid"]), build["sha256"], script_sha,
        str(script_stat["uid"]), git_sha, controller["python_sha256"],
        str(health["sampler"]["process_start_ticks"]),
        health["sampler"]["script_sha256"],
        str(health["sampler"]["csv_device"]),
        str(health["sampler"]["csv_inode"]),
        health["sampler"]["cmdline_sha256"],
        health["sampler"]["environ_sha256"],
        health["sampler"]["executable_sha256"],
        str(health["sampler"]["process_uid"]), source["sha256"],
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
    thermal_bytes = health["thermal"]["window_csv_ascii"].encode("ascii")
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
    appended_sampler = json.loads(
        canonical_bytes(health["sampler"]).decode("ascii")
    )
    appended_sampler["csv_stat"]["size"] += 4096
    appended_sampler["csv_stat"]["mtime_ns"] += 1
    validate_sampler_growth_binding(
        health["sampler"], appended_sampler, "synthetic append"
    )
    changed_sampler_mode = json.loads(
        canonical_bytes(appended_sampler).decode("ascii")
    )
    changed_sampler_mode["csv_stat"]["mode"] = 0o644
    try:
        validate_sampler_growth_binding(
            health["sampler"], changed_sampler_mode,
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
        ("2026-01-01T00:00:01Z", "1.000000000"),
        ("2026-01-01T00:00:02Z", "2.000000000"),
    ):
        row = [
            utc, monotonic_value, "50.0", "3000.0", "-1.0",
            *("40.0" for _ in range(8)), "1", "0.1", "0.1", "0.1",
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

    final_fixture = synthetic_final_bundle_fixture(health, binary, commit)
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
        health_args = argparse.Namespace(expected_binary_uid=os.getuid())
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
            "receipt": prepared_receipt, "violations": [],
        }
        observed_finish_intervals: List[
            Tuple[Optional[int], Optional[int]]
        ] = []
        original_popen = subprocess.Popen
        original_finish_health = globals()["finish_health"]
        original_read_tick = globals()["read_cpu_tick_receipt"]
        original_process_start_ticks = globals()["process_start_ticks"]
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
    print("WH2 direct systematic complement analyzer selftest passed")
    return 0


class AuthorityArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        fail("argument validation: " + message)


def argument_parser() -> AuthorityArgumentParser:
    parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest", action="store_true")
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
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--expected-binary-sha256")
    parser.add_argument("--expected-binary-uid", type=int)
    parser.add_argument("--expected-build-manifest-sha256")
    parser.add_argument("--expected-controller-sha256")
    parser.add_argument("--expected-controller-uid", type=int)
    parser.add_argument("--expected-git-sha256")
    parser.add_argument("--expected-python-sha256")
    parser.add_argument("--expected-sampler-process-start-ticks", type=int)
    parser.add_argument("--expected-sampler-script-sha256")
    parser.add_argument("--expected-sampler-csv-device", type=int)
    parser.add_argument("--expected-sampler-csv-inode", type=int)
    parser.add_argument("--expected-sampler-cmdline-sha256")
    parser.add_argument("--expected-sampler-environ-sha256")
    parser.add_argument("--expected-sampler-executable-sha256")
    parser.add_argument("--expected-sampler-uid", type=int)
    parser.add_argument("--expected-source-manifest-sha256")
    parser.add_argument("--expected-run-argv-sha256")
    return parser


def validate_cli_token_shape(argv: Sequence[str]) -> None:
    tokens = list(argv)
    if tokens == ["--selftest"]:
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


def require_authorized_process_uid(expected_uid: int) -> None:
    if not hasattr(os, "getresuid"):
        fail("controller UID policy requires Linux getresuid")
    if os.getresuid() != (expected_uid, expected_uid, expected_uid):
        fail("controller real/effective/saved UIDs differ from authorization")


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
    paths = (args.build_dir, args.sampler_script, args.sampler_csv)
    identities = (
        args.expected_binary_uid,
        args.expected_controller_uid,
        args.expected_sampler_process_start_ticks,
        args.expected_sampler_csv_device,
        args.expected_sampler_csv_inode,
        args.expected_sampler_uid,
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
        or args.expected_run_argv_sha256 is not None
    ):
        parser.error(
            "--run-once requires exact CPUs 120/121/122, canonical distinct "
            "binary/build/sampler paths, and every source/build/controller/"
            "binary/sampler identity seal"
        )
    args.controller_initial_affinity = []
    if live_process:
        require_authorized_process_uid(args.expected_controller_uid)
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
