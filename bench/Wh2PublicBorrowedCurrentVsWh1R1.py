#!/usr/bin/env python3
"""One-attempt supervisor/reducer for the public borrowed-source r1 screen.

The native worker owns all timing and emits one canonical JSONL transcript.
This controller only freezes provenance, validates the exact 1,442-record
roster, reduces matched log contrasts, and publishes the fixed one-attempt
bundle.  ``--selftest`` is synthetic: it does not launch a worker, create the
fixed namespace, or perform timing work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import shlex
import signal
import stat
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


CAMPAIGN = "wh2-public-borrowed-current-vs-wh1-r1"
SCHEMA_PREFIX = "wirehair.wh2.public-borrowed-current-vs-wh1-r1"
CONFIG_SCHEMA = SCHEMA_PREFIX + ".config.v1"
PANEL_SCHEMA = SCHEMA_PREFIX + ".panel.v1"
TERMINAL_SCHEMA = SCHEMA_PREFIX + ".terminal.v1"
SUMMARY_SCHEMA = SCHEMA_PREFIX + ".summary.v1"
PROVENANCE_SCHEMA = SCHEMA_PREFIX + ".provenance.v1"
CLAIM_SCHEMA = SCHEMA_PREFIX + ".claim.v1"
COMPLETE_SCHEMA = SCHEMA_PREFIX + ".complete.v1"
WORKER_STARTED_SCHEMA = SCHEMA_PREFIX + ".worker-started.v1"
PANEL_KEY_SCHEMA = SCHEMA_PREFIX + ".panel-key.v1"

FIXED_OUTPUT_DIR = Path("/var/tmp/wh2-public-borrowed-current-vs-wh1-r1")
R0_ROOTS = (
    Path("/var/tmp/wh2-v2-facade-default-parent-falsifier-r0.controller"),
    Path("/var/tmp/wh2-v2-facade-default-parent-falsifier-r0.root"),
    Path("/var/tmp/wh2-v2-facade-default-parent-falsifier-r0.sampler"),
    Path("/usr/local/libexec/wirehair/Wh2V2FacadeTimingLaunch.py"),
    Path("/var/lib/wirehair/wh2-v2-facade-default-parent-falsifier-r0-build-authority.json"),
)
CELLS: Tuple[Tuple[int, int], ...] = tuple(
    (k, block_bytes)
    for k in (8, 128, 512, 8192, 64000)
    for block_bytes in (64, 1280)
)
SCOPES = (
    "prebuilt-systematic",
    "fresh-systematic",
    "fresh-repair",
    "prebuilt-repair",
)
SCOPE_SPECS = {
    "prebuilt-systematic": ("0", False),
    "fresh-systematic": ("0", True),
    "fresh-repair": ("K", True),
    "prebuilt-repair": ("K", False),
}
SYSTEMATIC_SCOPES = frozenset((
    "prebuilt-systematic", "fresh-systematic"))
REPAIR_SCOPES = frozenset(("fresh-repair", "prebuilt-repair"))
COMPARISONS: Tuple[Tuple[str, str, str], ...] = (
    ("C/C", "C", "C"),
    ("L/L", "L", "L"),
    ("C/L", "C", "L"),
)
COMPARISON_ARMS = {
    name: (left, right) for name, left, right in COMPARISONS
}
AA_COMPARISONS = frozenset(("C/C", "L/L"))
REPLICATES = 12
INVOCATION_BUDGET = 65536
MINIMUM_INVOCATIONS = 24
INTERNAL_DEADLINE_SECONDS = 330
OUTER_DEADLINE_SECONDS = 360.0
TARGET_CPU = 120
TARGET_SIBLING = 56
CONTROLLER_CPU = 0
SIBLING_NON_IDLE_TICK_CAP = 5
THERMAL_ENDPOINT_MAX_MILLIC = 85000
T11_975 = 2.200985160082949
EXPECTED_PANEL_COUNT = 1440
EXPECTED_RECORD_COUNT = 1442
EXPECTED_MEASURED_INVOCATIONS = 852480
EXPECTED_WARMUP_INVOCATIONS = 2880
EXPECTED_MEASURED_ENCODE_CALLS = 185204736
EXPECTED_WARMUP_ENCODE_CALLS = 41955840
EXPECTED_ENCODE_CALLS = 227160576
SOURCE_SEED_BASE = 0xC199F24886210F53
SOURCE_SEED_BASE_TEXT = "0xc199f24886210f53"
SOURCE_SHA256 = {
    (8, 64, 64): "cba6894442c1997e328ad878e0e770b03a68b8eab70f2b03c630ae3c1c5eb5b4",
    (8, 64, 1): "e36a9993e92d7e9106d0e3fbe42d0ad259cdcd080044f0a599585fd78555a641",
    (8, 64, 63): "fea47c9ecca49a46e96266744a845e799599b64fceec41f31a96f2a27d665349",
    (8, 1280, 1280): "3f51720a962bd111409ae7495298431f80da83cc8cb9abee887d7c17f668a5d7",
    (8, 1280, 1): "fa8fa5499dc664577138d36ba4a3eb0e2b08fa09373c77b28cf51c604d84e95a",
    (8, 1280, 1279): "e8ce97e0b634eda74c154171b1fceaac92f541d7c5bac78c1a45831082059c00",
    (128, 64, 64): "658c12b5fdbcb59c4127570e8c365afd10e03b44926b4d914201ee0b69706fad",
    (128, 64, 1): "50663383e31016e6ae60797be76ac28aaf6a1e84449332d3d3462fd8f12b3cf1",
    (128, 64, 63): "4ff6ba4e2e7e98a41e24fc5a2887d4e0f8826d7e6ac7a4f5497f2dc4046fbf14",
    (128, 1280, 1280): "e62dfcb4ed8024009c497a53c5a0c34bf89e185ba1f720bc6b894c1ab48a0edf",
    (128, 1280, 1): "74aec531aad140deea34612d11d763802c227450891860f1fd61f8bdc266fc98",
    (128, 1280, 1279): "de1e21f35434cbad9daccc972177eff146b5b5e6c41ee967d78b6e0ef848c442",
    (512, 64, 64): "f7ddd00bddcab771aafb93e77ed5628e80151a1b220d8648e9a8d783c52e4d24",
    (512, 64, 1): "58f381a5f86cfa9c41f69c0308c881eed30566724040667a2722fc60bab9ab82",
    (512, 64, 63): "21f0175b7afd9c71201d5e30041f81974ccf1040b701278406b4b7d92838cae8",
    (512, 1280, 1280): "94790fe07ece812c776bde61701c718f7b33bb63d447a3d4f330a27657c9270d",
    (512, 1280, 1): "926e4f461b8da287a5d5e32eeded6e56ede8cbe1e80cfef00d3be1d911e20fd9",
    (512, 1280, 1279): "5b3970a1fb51642e6f3992a25cfad60d62e2be132c4f7fda78911d6046efe353",
    (8192, 64, 64): "ebbc9afef82223cd3bb9d9ab217104d6042d6fd3091416be0e86361bcb4c96c5",
    (8192, 64, 1): "d14d0042f8c9d0a0dc952e88f9c8661c5668710a681ed695d43ca787d43a4377",
    (8192, 64, 63): "ab6fcf51ff6976da64d5de0f55191c0ffd520b8761aeb1f83e14bccb4803151d",
    (8192, 1280, 1280): "79b809d539b6ea395db0c143e5e012d1400ea221a17b87f0bb971fbcb6d585bc",
    (8192, 1280, 1): "1599a99a970f9d411011b4eddef2bee6dde6b907502add9b24da1f75840b9db1",
    (8192, 1280, 1279): "ee5fce72e3a064a7a5107e65ce7007c2f7e3d22a5980853050f12784b61306a2",
    (64000, 64, 64): "8a14f4fc0252e4391665ab08fcd60a44b709e00894fd8c025c0e177feb35d931",
    (64000, 64, 1): "3633074091164c1c037829add6a76f6ac111ebdbac1dd784709248a8acefd4cf",
    (64000, 64, 63): "244416cd99df3c06d11ae3ed67a5d1412d74ea6569d85cd19539160614c9e5ec",
    (64000, 1280, 1280): "273aeab909363cfd4f1a09177ba5623f53520f1c170e1e75dcf6b208aa92b529",
    (64000, 1280, 1): "5a52b383b61ed658240aaf491c9f51392ccf1c8e492e010a19a0be85e9ad8e7c",
    (64000, 1280, 1279): "53b0f2ad551409d76bd7765f96cbfc0fa4f3795ac1532f2167b8e91a5516ea48",
}
MAX_INT63 = (1 << 63) - 1
MAX_LINE_BYTES = 1024 * 1024
MAX_RAW_BYTES = 64 * 1024 * 1024
MAX_STDERR_BYTES = 2 * 1024 * 1024
MAX_PROBE_BYTES = 8 * 1024 * 1024
LOWER40 = re.compile(r"^[0-9a-f]{40}$")
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
HEX_U64 = re.compile(r"^0x[0-9a-f]{16}$")
REQUIRED_DYNAMIC_SYMBOLS = frozenset((
    "wirehair_init_",
    "wirehair_encoder_create_ex",
    "wirehair_encode",
    "wirehair_decoder_create",
    "wirehair_decode",
    "wirehair_recover",
    "wirehair_free",
    "wirehair_v2_encoder_create_with_options",
    "wirehair_v2_profile_deserialize",
    "wirehair_v2_profile_validate",
    "wirehair_v2_encode",
    "wirehair_v2_decoder_create",
    "wirehair_v2_decode",
    "wirehair_v2_recover",
    "wirehair_v2_free",
))
FORBIDDEN_SYMBOL_FRAGMENTS = (
    "ForTesting",
    "TestMode",
    "wirehair_wh2_benchmark",
    "wirehair_wh2_test",
)
WORKER_FORBIDDEN_SYMBOL_FRAGMENTS = (
    "ForTesting",
    "ForBenchmark",
    "TestHook",
    "TestMode",
    "wirehair_v2::",
    "wirehair::Codec",
    "WirehairCodec::",
    "InitializeForValidatedSystem",
)
WIREHAIR_SYMBOL_VERSION = "WIREHAIR_2.0"
R1_TRACKED_INPUTS = (
    "CMakeLists.txt",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativePanel.cpp",
    "bench/Wh2NativePanel.h",
    "bench/Wh2PublicBorrowedCurrentVsWh1R1.cpp",
    "bench/Wh2PublicBorrowedCurrentVsWh1R1.py",
    "bench/Wh2PublicBorrowedTargetIdentity.cpp",
    "bench/Wh2PublicBorrowedTargetIdentity.h",
    "bench/Wh2RdpruTargetIdentityV2.cpp",
    "bench/Wh2RdpruTargetIdentityV2.h",
    "bench/test_Wh2PublicBorrowedCurrentVsWh1R1.py",
    "cmake/Wh2PublicBorrowedCurrentVsWh1R1SymbolAudit.cmake",
)
R1_TARGET_SOURCES = (
    "bench/Wh2PublicBorrowedCurrentVsWh1R1.cpp",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2NativePanel.cpp",
    "bench/Wh2PublicBorrowedTargetIdentity.cpp",
    "bench/Wh2RdpruTargetIdentityV2.cpp",
)
R1_LIBRARY_SOURCES = (
    "wirehair.cpp",
    "gf256.cpp",
    "WirehairCodec.cpp",
    "WirehairTools.cpp",
    "codec/WirehairV2Codec.cpp",
    "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Plan.cpp",
    "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.cpp",
    "codec/WirehairV2Solve.cpp",
)
OUTPUT_CORE_FILES = (
    "CLAIM", "raw.jsonl", "worker.stderr", "summary.json",
    "provenance.json", "COMPLETE",
)
WORKER_STARTED_FILE = "WORKER_STARTED"


class ValidationError(RuntimeError):
    """Evidence invalidity, distinct from a valid scientific rejection."""


def fail(message: str) -> None:
    raise ValidationError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_sha256(path: Path, maximum: Optional[int] = None) -> str:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        while True:
            data = source.read(1024 * 1024)
            if not data:
                break
            size += len(data)
            if maximum is not None and size > maximum:
                fail("{} exceeds its byte cap".format(path))
            digest.update(data)
    return digest.hexdigest()


def unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("duplicate JSON key: {}".format(key))
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    fail("non-finite JSON number: {}".format(value))


def parse_canonical_line(data: bytes, where: str) -> Dict[str, Any]:
    if (
        not data or len(data) > MAX_LINE_BYTES or
        not data.endswith(b"\n") or data.count(b"\n") != 1 or b"\r" in data
    ):
        fail("{} is not one bounded LF-terminated line".format(where))
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail("{} is malformed: {}".format(where, exc))
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail("{} is not a canonical JSON object".format(where))
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], where: str) -> None:
    actual = set(value.keys())
    wanted = set(expected)
    if actual != wanted:
        fail("{} keys differ: expected {}, got {}".format(
            where, sorted(wanted), sorted(actual)))


def exact_int(value: Any, lower: int, upper: int, where: str) -> int:
    if type(value) is not int or value < lower or value > upper:
        fail("{} is not an integer in [{},{}]".format(where, lower, upper))
    return value


def exact_bool(value: Any, expected: bool, where: str) -> bool:
    if type(value) is not bool or value is not expected:
        fail("{} differs".format(where))
    return value


def exact_string(value: Any, expected: str, where: str) -> str:
    if type(value) is not str or value != expected:
        fail("{} differs".format(where))
    return value


def bounded_string(value: Any, where: str, maximum: int = 16384) -> str:
    if type(value) is not str or not value or len(value) > maximum:
        fail("{} is not a bounded nonempty string".format(where))
    return value


def lower_hash(value: Any, where: str) -> str:
    if type(value) is not str or LOWER64.fullmatch(value) is None:
        fail("{} is not a lowercase SHA-256".format(where))
    return value


def lower_commit(value: Any, where: str) -> str:
    if type(value) is not str or LOWER40.fullmatch(value) is None:
        fail("{} is not a lowercase Git object ID".format(where))
    return value


def total_invocations(k: int) -> int:
    return max(MINIMUM_INVOCATIONS, (INVOCATION_BUDGET + k - 1) // k)


def invocations_by_replicate(k: int) -> List[int]:
    quotient, remainder = divmod(total_invocations(k), REPLICATES)
    return [
        quotient + (1 if replicate < remainder else 0)
        for replicate in range(REPLICATES)
    ]


def message_bytes(k: int, block_bytes: int) -> int:
    return k * block_bytes


def source_seed_text(k: int, block_bytes: int, tail_bytes: int) -> str:
    value = SOURCE_SEED_BASE ^ (k << 33) ^ (block_bytes << 1) ^ tail_bytes
    return "0x{:016x}".format(value & ((1 << 64) - 1))


def band_for(k: int) -> str:
    if k in (8, 128):
        return "small"
    if k in (512, 8192):
        return "medium"
    if k == 64000:
        return "large"
    raise AssertionError("unfrozen K")


def panel_key_sha256(
    k: int, block_bytes: int, scope: str, comparison: str,
) -> str:
    return sha256_bytes(canonical_bytes({
        "K": k,
        "block_bytes": block_bytes,
        "campaign": CAMPAIGN,
        "comparison": comparison,
        "schema": PANEL_KEY_SCHEMA,
        "scope": scope,
    }))


def panel_order(
    replicate: int, k: int, block_bytes: int, scope: str, comparison: str,
) -> str:
    digest = panel_key_sha256(k, block_bytes, scope, comparison)
    phase = bytes.fromhex(digest)[-1] & 1
    return "ABBA" if (replicate & 1) == phase else "BAAB"


def sides_for(order: str) -> Tuple[str, ...]:
    if order == "ABBA":
        return (
            "left", "right", "right", "left",
            "right", "left", "left", "right",
        )
    if order == "BAAB":
        return (
            "right", "left", "left", "right",
            "left", "right", "right", "left",
        )
    fail("unknown panel order")
    raise AssertionError


def lane_contrast(values: Sequence[int], order: str) -> float:
    if len(values) != 8 or any(
        type(value) is not int or value <= 0 for value in values
    ):
        fail("lane contrast requires eight positive integer sums")
    logs = [math.log(value) for value in values]

    def contrast(offset: int, block_order: str) -> float:
        if block_order == "ABBA":
            return ((logs[offset] - logs[offset + 1]) +
                    (logs[offset + 3] - logs[offset + 2])) / 2.0
        return ((logs[offset + 1] - logs[offset]) +
                (logs[offset + 2] - logs[offset + 3])) / 2.0

    opposite = "BAAB" if order == "ABBA" else "ABBA"
    return 0.5 * (contrast(0, order) + contrast(4, opposite))


def sample_summary(values: Sequence[float]) -> Dict[str, Any]:
    if len(values) != REPLICATES or any(
        type(item) not in (int, float) or not math.isfinite(item)
        for item in values
    ):
        fail("statistical sample is incomplete or non-finite")
    mean = math.fsum(values) / REPLICATES
    variance = math.fsum(
        (item - mean) ** 2 for item in values) / (REPLICATES - 1)
    if not math.isfinite(variance) or variance < 0.0:
        fail("statistical variance is invalid")
    standard_error = math.sqrt(variance / REPLICATES)
    lower_log = mean - T11_975 * standard_error
    upper_log = mean + T11_975 * standard_error
    return {
        "geometric_mean": math.exp(mean),
        "log_mean": mean,
        "log_standard_error": standard_error,
        "lower95": math.exp(lower_log),
        "lower95_log": lower_log,
        "n": REPLICATES,
        "upper95": math.exp(upper_log),
        "upper95_log": upper_log,
    }


CONFIG_KEYS = {
    "arm_descriptors", "campaign", "cells", "comparisons", "compile",
    "expected_encode_calls",
    "expected_measured_encode_calls", "expected_measured_invocations",
    "expected_panels", "expected_records", "expected_warmup_encode_calls",
    "expected_warmup_invocations", "internal_deadline_seconds",
    "invocation_budget", "minimum_invocations",
    "native_config_identity_sha256", "panel_key_schema",
    "panel_replicates", "roster_order", "runtime_library_maps_sha256",
    "runtime_library_path", "schema", "scopes", "sibling_cpu",
    "source_seed_base", "target_cpu", "target_identity",
}
ARM_DESCRIPTOR_KEYS = {"arm", "descriptor", "descriptor_sha256"}
COMPILE_KEYS = {
    "compiler_path", "compiler_sha256", "compiler_version",
    "harness_git_commit", "implementation_git_commit",
}
CELL_KEYS = {
    "K", "block_bytes", "final_bytes", "invocations_by_replicate",
    "message_bytes", "oracles", "partial_final_oracles", "source_seed",
    "source_sha256",
}
ORACLE_KEYS = {
    "arm", "attached_overlap_no_write_verified",
    "borrowed_eligible_systematic_ids", "descriptor_sha256",
    "equation_configuration_sha256", "first_repair_sha256",
    "full_repair_sha256", "high_id_sha256", "profile_bytes", "profile_hex",
    "profile_sha256", "roundtrip_first_id", "roundtrip_packet_count",
    "roundtrip_repair_only_verified", "roundtrip_sha256",
    "roundtrip_verified", "systematic_sha256",
}
PARTIAL_KEYS = {"arms", "source_sha256", "tail_bytes"}
PARTIAL_ARM_KEYS = {
    "arm", "attached_overlap_no_write_verified", "first_repair_sha256",
    "high_id_sha256", "profile_sha256", "roundtrip_first_id",
    "roundtrip_packet_count", "roundtrip_repair_only_verified",
    "roundtrip_sha256", "source_immutable_verified", "systematic_sha256",
    "systematic_verified",
}
TARGET_IDENTITY_KEYS = {
    "after_cpu", "before_cpu", "canonical_bytes", "canonical_hex",
    "canonical_sha256", "capabilities", "derived", "raw_capture_complete",
    "requested_cpu", "semantic_validation_passed",
    "singleton_affinity_verified",
}
TARGET_CAPABILITY_KEYS = {
    "leaf1_ecx", "leaf1_edx", "leaf6_eax", "leaf6_ecx",
    "leaf80000001_ecx", "leaf80000001_edx", "leaf80000008_ebx",
    "leaf80000021_eax", "max_basic_leaf", "max_extended_leaf",
}
TARGET_DERIVED_KEYS = {
    "ccd_id", "complex_id", "core_id", "family", "full_apic_id",
    "initial_apic_id_8", "logical_processors_per_package", "model",
    "package_id", "stepping", "thread_id", "threads_per_core",
}
WORKER_STARTED_KEYS = {
    "campaign", "claim_sha256", "native_config_identity_sha256", "schema",
    "source_commit", "status", "worker_sha256",
}
COMPARISON_KEYS = {"left_arm", "name", "right_arm"}
SCOPE_KEYS = {"first_id", "name", "timed_construction"}
PANEL_KEYS = {
    "K", "affinity_verified", "block_bytes", "comparison",
    "cpu_migration_verified", "invocations_per_slot", "left_arm",
    "left_descriptor_sha256", "left_outcome_code", "left_output_sha256",
    "left_public_encode_calls", "order", "panel_key_sha256", "replicate",
    "right_arm", "right_descriptor_sha256", "right_outcome_code",
    "right_output_sha256", "right_public_encode_calls",
    "runtime_library_maps_sha256", "schema", "scope", "sequence", "slots",
    "target_cpu",
}
SLOT_KEYS = {"elapsed_ns", "invocation_count", "side"}
TERMINAL_KEYS = {
    "encode_call_count", "measured_invocation_count", "panel_count",
    "record_count", "schema", "status", "warmup_invocation_count",
}
CLAIM_KEYS = {
    "campaign", "controller_interpreter_sha256", "controller_sha256",
    "created_unix_ns", "library_sha256", "native_config_identity_sha256",
    "schema", "source_commit", "worker_sha256",
}
SUMMARY_KEYS = {
    "campaign", "elapsed_seconds", "infrastructure_failures", "outcome",
    "raw_sha256", "reject_gates", "schema", "statistics",
    "transcript_invalid_gates",
}
COMPLETE_KEYS = {
    "campaign", "claim_sha256", "outcome", "provenance_sha256",
    "raw_sha256", "schema", "status", "summary_sha256",
}
STAT_RECEIPT_KEYS = {
    "device", "inode", "mode", "mtime_ns", "path", "size",
}
ARTIFACT_KEYS = STAT_RECEIPT_KEYS | {"sha256"}
COMPILER_RECEIPT_KEYS = {
    "path", "sha256", "version", "version_sha256", "version_text",
    "version_text_sha256",
}
INTERPRETER_RECEIPT_KEYS = {
    "artifact", "dont_write_bytecode", "ignore_environment",
    "implementation", "invoked_path", "isolated", "no_site",
    "no_user_site", "path", "version", "version_info",
}
GIT_RECEIPT_KEYS = {
    "commit", "git_artifact", "git_invoked_path", "git_version_sha256",
    "git_version_text", "required_tracked_inputs", "tracked_clean", "tree",
    "upstream_commit",
}
BUILD_RECEIPT_KEYS = {
    "cache", "cmake_build_root", "cmake_compiler",
    "cmake_compiler_resolved", "cmake_interpreter",
    "cmake_interpreter_resolved", "cmake_source_root",
    "compile_command_sha256", "compile_commands", "compile_definitions",
    "compiler_invocation",
    "dry_run_sha256", "dry_run_text", "library_compile_command_sha256",
    "library_compile_commands", "library_compile_definitions",
    "library_link_command", "library_link_command_sha256",
    "library_link_object_roster", "library_sources", "link_command",
    "link_command_sha256", "link_object_roster", "ninja_commands_sha256",
    "ninja_commands_text", "target_sources",
}
DYNAMIC_RECEIPT_KEYS = {
    "forbidden_symbol_fragments", "ldd_library_path",
    "ldd_library_resolved_path", "ldd_library_stat",
    "ldd_normalized_sha256", "ldd_normalized_text", "library_nm_argv",
    "library_symbol_table_text", "nm_artifact", "nm_invoked_path",
    "nm_version_sha256", "nm_version_text", "readelf_dynamic_sha256",
    "readelf_dynamic_text", "required_dynamic_symbols",
    "required_dynamic_symbols_verified", "symbol_table_sha256",
    "worker_forbidden_symbol_fragments", "worker_full_nm_argv",
    "worker_full_symbol_table_sha256", "worker_full_symbol_table_text",
    "worker_import_symbol_version", "worker_import_symbols",
    "worker_import_symbols_verified", "worker_nm_argv",
    "worker_symbol_table_sha256", "worker_symbol_table_text",
}
NATIVE_CONFIG_RECEIPT_KEYS = {
    "config_sha256", "native_config_identity_sha256",
    "runtime_library_maps_sha256", "runtime_library_path", "stderr_sha256",
    "target_identity_sha256",
}
EXPECTED_KEYS = {
    "compiler_path", "compiler_sha256", "compiler_version", "library_path",
    "native_config_identity_sha256", "source_commit",
    "target_identity_sha256",
}
R0_KEYS = {"entries", "snapshot_sha256"}
HEALTH_ADJUDICATION_KEYS = {
    "edac_no_increase_verified", "governor_stable_verified",
    "sibling_non_idle_tick_cap", "sibling_non_idle_tick_delta",
    "thermal_endpoint_max_millic",
    "thermal_endpoint_samples_within_cap_verified",
    "thermal_roster_stable_verified", "topology_stable_verified",
}
PROVENANCE_POSTFLIGHT_FIELDS = (
    "health_after", "health_adjudication", "r0_after", "artifacts_after",
    "interpreter_after", "git_after", "build_after", "dynamic_after",
)
PROVENANCE_KEYS = {
    "artifacts", "artifacts_after", "build", "build_after", "build_root",
    "campaign", "claim_sha256", "compiler", "controller_affinity",
    "controller_path", "dynamic",
    "dynamic_after", "expected", "git", "git_after",
    "health_adjudication", "health_after", "health_before", "interpreter",
    "interpreter_after", "library_path", "native_config",
    "outer_deadline_seconds", "r0_after", "r0_before", "raw_sha256",
    "schema", "source_root", "stderr_sha256", "worker_exit", "worker_path",
    "worker_started_sha256", "worker_timed_out",
}


def validate_target_identity(value: Any, where: str) -> Dict[str, Any]:
    if type(value) is not dict:
        fail("{} is not an object".format(where))
    exact_keys(value, TARGET_IDENTITY_KEYS, where)
    for field in ("requested_cpu", "before_cpu", "after_cpu"):
        exact_int(value[field], TARGET_CPU, TARGET_CPU, where + "." + field)
    exact_bool(value["raw_capture_complete"], True,
               where + ".raw_capture_complete")
    exact_bool(value["semantic_validation_passed"], True,
               where + ".semantic_validation_passed")
    exact_bool(value["singleton_affinity_verified"], True,
               where + ".singleton_affinity_verified")
    byte_count = exact_int(
        value["canonical_bytes"], 1, 4096, where + ".canonical_bytes")
    canonical_hex = bounded_string(
        value["canonical_hex"], where + ".canonical_hex",
        maximum=byte_count * 2)
    if re.fullmatch(r"[0-9a-f]{%d}" % (byte_count * 2), canonical_hex) is None:
        fail("{} canonical bytes are not exact lowercase hex".format(where))
    canonical = bytes.fromhex(canonical_hex)
    exact_string(value["canonical_sha256"], sha256_bytes(canonical),
                 where + ".canonical_sha256")
    capabilities = value["capabilities"]
    derived = value["derived"]
    if type(capabilities) is not dict or type(derived) is not dict:
        fail("{} capability/derived receipts are malformed".format(where))
    exact_keys(capabilities, TARGET_CAPABILITY_KEYS, where + ".capabilities")
    exact_keys(derived, TARGET_DERIVED_KEYS, where + ".derived")
    for field in TARGET_CAPABILITY_KEYS:
        exact_int(capabilities[field], 0, (1 << 32) - 1,
                  where + ".capabilities." + field)
    for field in TARGET_DERIVED_KEYS:
        exact_int(derived[field], 0, (1 << 32) - 1,
                  where + ".derived." + field)
    exact_int(derived["initial_apic_id_8"], 0, 255,
              where + ".derived.initial_apic_id_8")
    for field in ("threads_per_core", "logical_processors_per_package"):
        exact_int(derived[field], 1, (1 << 32) - 1,
                  where + ".derived." + field)
    return value


def native_config_identity_sha256(value: Mapping[str, Any]) -> str:
    projection = dict(value)
    projection["runtime_library_maps_sha256"] = "0" * 64
    projection["native_config_identity_sha256"] = "0" * 64
    return sha256_bytes(canonical_bytes(projection) + b"\n")


def validate_oracle(
    value: Any, arm: str, k: int, source_sha256: str, where: str,
) -> Dict[str, Any]:
    if type(value) is not dict:
        fail("{} is not an object".format(where))
    exact_keys(value, ORACLE_KEYS, where)
    exact_string(value["arm"], arm, where + ".arm")
    if arm == "C":
        exact_bool(value["attached_overlap_no_write_verified"], True,
                   where + ".attached_overlap_no_write_verified")
    elif value["attached_overlap_no_write_verified"] is not None:
        fail("{} WH1 overlap field is not null".format(where))
    exact_bool(value["roundtrip_verified"], True, where + ".roundtrip_verified")
    exact_bool(value["roundtrip_repair_only_verified"], True,
               where + ".roundtrip_repair_only_verified")
    exact_int(value["roundtrip_first_id"], k, k,
              where + ".roundtrip_first_id")
    exact_int(value["roundtrip_packet_count"], k, k + 64,
              where + ".roundtrip_packet_count")
    exact_int(value["borrowed_eligible_systematic_ids"],
              k if arm == "C" else 0, k if arm == "C" else 0,
              where + ".borrowed_eligible_systematic_ids")
    for field in (
        "descriptor_sha256", "equation_configuration_sha256",
        "first_repair_sha256", "full_repair_sha256", "high_id_sha256",
        "roundtrip_sha256", "systematic_sha256",
    ):
        lower_hash(value[field], where + "." + field)
    exact_string(value["systematic_sha256"], source_sha256,
                 where + ".systematic_sha256")
    exact_string(value["roundtrip_sha256"], source_sha256,
                 where + ".roundtrip_sha256")
    profile_bytes = exact_int(value["profile_bytes"], 0, 4096,
                              where + ".profile_bytes")
    if arm == "C":
        if profile_bytes <= 0:
            fail("{} lacks the serialized V2 profile".format(where))
        profile_hex = bounded_string(
            value["profile_hex"], where + ".profile_hex",
            maximum=profile_bytes * 2)
        if re.fullmatch(r"[0-9a-f]{%d}" % (profile_bytes * 2), profile_hex) is None:
            fail("{} profile is not exact lowercase hexadecimal".format(where))
        try:
            serialized = bytes.fromhex(profile_hex)
        except ValueError:
            fail("{} profile is not lowercase hexadecimal".format(where))
        if len(serialized) != profile_bytes:
            fail("{} decoded profile length differs".format(where))
        exact_string(value["profile_sha256"], sha256_bytes(serialized),
                     where + ".profile_sha256")
    elif value["profile_hex"] is not None or value["profile_sha256"] is not None:
        fail("{} WH1 profile fields are not null".format(where))
    return value


def validate_partial_oracles(
    value: Any, k: int, block_bytes: int, where: str,
) -> None:
    if type(value) is not list or len(value) != 2:
        fail("{} roster differs".format(where))
    for index, (partial, tail) in enumerate(zip(value, (1, block_bytes - 1))):
        item_where = "{}[{}]".format(where, index)
        if type(partial) is not dict:
            fail("{} is not an object".format(item_where))
        exact_keys(partial, PARTIAL_KEYS, item_where)
        exact_int(partial["tail_bytes"], tail, tail, item_where + ".tail_bytes")
        exact_string(
            partial["source_sha256"], SOURCE_SHA256[(k, block_bytes, tail)],
            item_where + ".source_sha256")
        arms = partial["arms"]
        if type(arms) is not list or len(arms) != 2:
            fail("{} arm roster differs".format(item_where))
        for arm_index, (receipt, arm) in enumerate(zip(arms, ("C", "L"))):
            arm_where = "{}.arms[{}]".format(item_where, arm_index)
            if type(receipt) is not dict:
                fail("{} is not an object".format(arm_where))
            exact_keys(receipt, PARTIAL_ARM_KEYS, arm_where)
            exact_string(receipt["arm"], arm, arm_where + ".arm")
            if arm == "C":
                exact_bool(receipt["attached_overlap_no_write_verified"], True,
                           arm_where + ".attached_overlap_no_write_verified")
            elif receipt["attached_overlap_no_write_verified"] is not None:
                fail("{} WH1 overlap field is not null".format(arm_where))
            exact_bool(receipt["source_immutable_verified"], True,
                       arm_where + ".source_immutable_verified")
            exact_bool(receipt["systematic_verified"], True,
                       arm_where + ".systematic_verified")
            exact_bool(receipt["roundtrip_repair_only_verified"], True,
                       arm_where + ".roundtrip_repair_only_verified")
            exact_int(receipt["roundtrip_first_id"], k, k,
                      arm_where + ".roundtrip_first_id")
            exact_int(receipt["roundtrip_packet_count"], k, k + 64,
                      arm_where + ".roundtrip_packet_count")
            for field in (
                "first_repair_sha256", "high_id_sha256", "roundtrip_sha256",
                "systematic_sha256",
            ):
                lower_hash(receipt[field], arm_where + "." + field)
            exact_string(receipt["roundtrip_sha256"], partial["source_sha256"],
                         arm_where + ".roundtrip_sha256")
            if arm == "C":
                lower_hash(receipt["profile_sha256"],
                           arm_where + ".profile_sha256")
            elif receipt["profile_sha256"] is not None:
                fail("{} WH1 profile SHA-256 is not null".format(arm_where))
        if arms[0]["first_repair_sha256"] == arms[1]["first_repair_sha256"]:
            fail("{} C/L first repair receipts match".format(item_where))
        if arms[0]["high_id_sha256"] == arms[1]["high_id_sha256"]:
            fail("{} C/L high-ID receipts match".format(item_where))
        exact_string(
            arms[0]["systematic_sha256"], arms[1]["systematic_sha256"],
            item_where + ".C/L systematic_sha256")


def validate_config(
    value: Dict[str, Any], expected: Mapping[str, str],
) -> Dict[Tuple[int, int, str], Dict[str, Any]]:
    exact_keys(value, CONFIG_KEYS, "config")
    exact_string(value["campaign"], CAMPAIGN, "config.campaign")
    exact_string(value["schema"], CONFIG_SCHEMA, "config.schema")
    exact_string(value["panel_key_schema"], PANEL_KEY_SCHEMA,
                 "config.panel_key_schema")
    lower_hash(value["runtime_library_maps_sha256"],
               "config.runtime_library_maps_sha256")
    exact_string(value["runtime_library_path"], expected["library_path"],
                 "config.runtime_library_path")
    exact_int(value["target_cpu"], TARGET_CPU, TARGET_CPU, "config.target_cpu")
    exact_int(value["sibling_cpu"], TARGET_SIBLING, TARGET_SIBLING,
              "config.sibling_cpu")
    target_identity = validate_target_identity(
        value["target_identity"], "config.target_identity")
    target_identity_sha256 = target_identity["canonical_sha256"]
    if "target_identity_sha256" in expected:
        exact_string(
            target_identity_sha256, expected["target_identity_sha256"],
            "config.target_identity.canonical_sha256")
    config_identity = lower_hash(
        value["native_config_identity_sha256"],
        "config.native_config_identity_sha256")
    if "native_config_identity_sha256" in expected:
        exact_string(
            config_identity, expected["native_config_identity_sha256"],
            "config.native_config_identity_sha256")
    if value["roster_order"] != ["replicate", "cell", "scope", "comparison"]:
        fail("config roster order differs")
    for field, number in (
        ("internal_deadline_seconds", INTERNAL_DEADLINE_SECONDS),
        ("invocation_budget", INVOCATION_BUDGET),
        ("minimum_invocations", MINIMUM_INVOCATIONS),
        ("panel_replicates", REPLICATES),
        ("expected_panels", EXPECTED_PANEL_COUNT),
        ("expected_records", EXPECTED_RECORD_COUNT),
        ("expected_measured_invocations", EXPECTED_MEASURED_INVOCATIONS),
        ("expected_warmup_invocations", EXPECTED_WARMUP_INVOCATIONS),
        ("expected_measured_encode_calls", EXPECTED_MEASURED_ENCODE_CALLS),
        ("expected_warmup_encode_calls", EXPECTED_WARMUP_ENCODE_CALLS),
        ("expected_encode_calls", EXPECTED_ENCODE_CALLS),
    ):
        exact_int(value[field], number, number, "config." + field)
    exact_string(value["source_seed_base"], SOURCE_SEED_BASE_TEXT,
                 "config.source_seed_base")

    compile_value = value["compile"]
    if type(compile_value) is not dict:
        fail("config.compile is not an object")
    exact_keys(compile_value, COMPILE_KEYS, "config.compile")
    exact_string(compile_value["compiler_path"], expected["compiler_path"],
                 "config.compile.compiler_path")
    exact_string(compile_value["compiler_sha256"], expected["compiler_sha256"],
                 "config.compile.compiler_sha256")
    exact_string(compile_value["compiler_version"], expected["compiler_version"],
                 "config.compile.compiler_version")
    exact_string(compile_value["harness_git_commit"], expected["source_commit"],
                 "config.compile.harness_git_commit")
    exact_string(compile_value["implementation_git_commit"],
                 expected["source_commit"],
                 "config.compile.implementation_git_commit")

    expected_descriptors = {
        "C": {
            "api": "wirehair_v2_encoder_create_with_options",
            "arm": "C", "codec": "wirehair2",
            "equation_profile": "certified-2026-07",
            "schema": SCHEMA_PREFIX + ".arm.v1",
            "source_policy": "borrowed-immutable",
        },
        "L": {
            "api": "wirehair_encoder_create_ex", "arm": "L",
            "codec": "wirehair1", "equation_profile": "legacy-current",
            "schema": SCHEMA_PREFIX + ".arm.v1",
            "source_policy": "borrowed",
        },
    }
    descriptors = value["arm_descriptors"]
    if type(descriptors) is not list or len(descriptors) != 2:
        fail("config arm descriptor roster differs")
    for index, (entry, arm) in enumerate(zip(descriptors, ("C", "L"))):
        where = "config.arm_descriptors[{}]".format(index)
        if type(entry) is not dict:
            fail("{} is not an object".format(where))
        exact_keys(entry, ARM_DESCRIPTOR_KEYS, where)
        exact_string(entry["arm"], arm, where + ".arm")
        if entry["descriptor"] != expected_descriptors[arm]:
            fail("{} descriptor differs".format(where))
        exact_string(entry["descriptor_sha256"],
                     sha256_bytes(canonical_bytes(expected_descriptors[arm])),
                     where + ".descriptor_sha256")

    expected_comparisons = [
        {"left_arm": left, "name": name, "right_arm": right}
        for name, left, right in COMPARISONS
    ]
    if value["comparisons"] != expected_comparisons:
        fail("config comparison roster differs")
    expected_scopes = []
    for name in SCOPES:
        first_id, timed = SCOPE_SPECS[name]
        expected_scopes.append({
            "first_id": first_id,
            "name": name,
            "timed_construction": timed,
        })
    if value["scopes"] != expected_scopes:
        fail("config scope roster differs")

    cells = value["cells"]
    if type(cells) is not list or len(cells) != len(CELLS):
        fail("config cell roster length differs")
    oracle_map: Dict[Tuple[int, int, str], Dict[str, Any]] = {}
    for index, ((k, block_bytes), cell) in enumerate(zip(CELLS, cells)):
        where = "config.cells[{}]".format(index)
        if type(cell) is not dict:
            fail("{} is not an object".format(where))
        exact_keys(cell, CELL_KEYS, where)
        exact_int(cell["K"], k, k, where + ".K")
        exact_int(cell["block_bytes"], block_bytes, block_bytes,
                  where + ".block_bytes")
        exact_int(cell["final_bytes"], block_bytes, block_bytes,
                  where + ".final_bytes")
        exact_int(cell["message_bytes"], message_bytes(k, block_bytes),
                  message_bytes(k, block_bytes), where + ".message_bytes")
        if cell["invocations_by_replicate"] != invocations_by_replicate(k):
            fail("{} invocation allocation differs".format(where))
        exact_string(
            cell["source_seed"], source_seed_text(k, block_bytes, block_bytes),
            where + ".source_seed")
        source_hash = SOURCE_SHA256[(k, block_bytes, block_bytes)]
        exact_string(cell["source_sha256"], source_hash,
                     where + ".source_sha256")
        validate_partial_oracles(
            cell["partial_final_oracles"], k, block_bytes,
            where + ".partial_final_oracles")
        oracles = cell["oracles"]
        if type(oracles) is not list or len(oracles) != 2:
            fail("{} oracle roster differs".format(where))
        for arm, oracle in zip(("C", "L"), oracles):
            validated = validate_oracle(
                oracle, arm, k, source_hash, where + ".oracles[{}]".format(arm))
            oracle_map[(k, block_bytes, arm)] = validated
        if (oracles[0]["full_repair_sha256"] ==
                oracles[1]["full_repair_sha256"]):
            fail("{} C and L repair streams unexpectedly match".format(where))
    exact_string(
        config_identity, native_config_identity_sha256(value),
        "config recomputed native-config identity")
    return oracle_map


def validate_panel(
    value: Dict[str, Any], replicate: int, k: int, block_bytes: int,
    scope: str, comparison: str, sequence: int,
    oracles: Mapping[Tuple[int, int, str], Mapping[str, Any]],
    runtime_maps_sha256: str,
) -> Tuple[float, int, int]:
    where = "panel {}".format(sequence)
    exact_keys(value, PANEL_KEYS, where)
    exact_string(value["schema"], PANEL_SCHEMA, where + ".schema")
    exact_int(value["sequence"], sequence, sequence, where + ".sequence")
    exact_int(value["K"], k, k, where + ".K")
    exact_int(value["block_bytes"], block_bytes, block_bytes,
              where + ".block_bytes")
    exact_int(value["replicate"], replicate, replicate, where + ".replicate")
    exact_string(value["scope"], scope, where + ".scope")
    exact_string(value["comparison"], comparison, where + ".comparison")
    left, right = COMPARISON_ARMS[comparison]
    exact_string(value["left_arm"], left, where + ".left_arm")
    exact_string(value["right_arm"], right, where + ".right_arm")
    exact_string(value["left_descriptor_sha256"],
                 oracles[(k, block_bytes, left)]["descriptor_sha256"],
                 where + ".left_descriptor_sha256")
    exact_string(value["right_descriptor_sha256"],
                 oracles[(k, block_bytes, right)]["descriptor_sha256"],
                 where + ".right_descriptor_sha256")
    for side in ("left", "right"):
        exact_int(value[side + "_outcome_code"], 0, 0,
                  where + "." + side + "_outcome_code")
        arm = left if side == "left" else right
        expected_output = oracles[(k, block_bytes, arm)][
            "systematic_sha256" if scope in SYSTEMATIC_SCOPES
            else "full_repair_sha256"]
        exact_string(value[side + "_output_sha256"], expected_output,
                     where + "." + side + "_output_sha256")
    exact_bool(value["affinity_verified"], True, where + ".affinity_verified")
    exact_bool(value["cpu_migration_verified"], True,
               where + ".cpu_migration_verified")
    exact_int(value["target_cpu"], TARGET_CPU, TARGET_CPU, where + ".target_cpu")
    exact_string(value["runtime_library_maps_sha256"], runtime_maps_sha256,
                 where + ".runtime_library_maps_sha256")
    expected_key = panel_key_sha256(k, block_bytes, scope, comparison)
    exact_string(value["panel_key_sha256"], expected_key,
                 where + ".panel_key_sha256")
    order = panel_order(replicate, k, block_bytes, scope, comparison)
    exact_string(value["order"], order, where + ".order")
    n = invocations_by_replicate(k)[replicate]
    exact_int(value["invocations_per_slot"], n, n,
              where + ".invocations_per_slot")
    slots = value["slots"]
    if type(slots) is not list or len(slots) != 8:
        fail("{} slot roster differs".format(where))
    elapsed: List[int] = []
    left_count = 0
    right_count = 0
    for index, (slot, side) in enumerate(zip(slots, sides_for(order))):
        slot_where = "{}.slots[{}]".format(where, index)
        if type(slot) is not dict:
            fail("{} is not an object".format(slot_where))
        exact_keys(slot, SLOT_KEYS, slot_where)
        exact_string(slot["side"], side, slot_where + ".side")
        count = (n + 1) // 2 if index < 4 else n // 2
        exact_int(slot["invocation_count"], count, count,
                  slot_where + ".invocation_count")
        elapsed.append(exact_int(
            slot["elapsed_ns"], count, MAX_INT63, slot_where + ".elapsed_ns"))
        if side == "left":
            left_count += count
        else:
            right_count += count
    if left_count != 2 * n or right_count != 2 * n:
        fail("{} side invocation accounting differs".format(where))
    left_calls = (left_count + 1) * k
    right_calls = (right_count + 1) * k
    exact_int(value["left_public_encode_calls"], left_calls, left_calls,
              where + ".left_public_encode_calls")
    exact_int(value["right_public_encode_calls"], right_calls, right_calls,
              where + ".right_public_encode_calls")
    return lane_contrast(elapsed, order), left_count + right_count, \
        left_calls + right_calls


def aggregate_logs(
    logs: Mapping[Tuple[str, int, int, str], List[float]],
    scope: str, cells: Sequence[Tuple[int, int]],
) -> List[float]:
    return [
        math.fsum(logs[("C/L", k, block_bytes, scope)][replicate]
                  for k, block_bytes in cells) / len(cells)
        for replicate in range(REPLICATES)
    ]


def adjudicate(
    logs: Mapping[Tuple[str, int, int, str], List[float]],
) -> Tuple[Dict[str, Any], List[str], List[str]]:
    invalid: List[str] = []
    reject: List[str] = []
    cells_out: Dict[str, Any] = {}
    aa_limit = math.log1p(0.02)
    point_limit = math.log(1.02)
    systematic_limit = math.log(0.99)
    repair_limit = math.log(1.02)
    for comparison, _, _ in COMPARISONS:
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                key = (comparison, k, block_bytes, scope)
                label = "{}:{}:K{}:B{}".format(
                    comparison, scope, k, block_bytes)
                summary = sample_summary(logs[key])
                cells_out[label] = summary
                if comparison in AA_COMPARISONS:
                    if not (
                        summary["lower95_log"] > -aa_limit and
                        summary["upper95_log"] < aa_limit
                    ):
                        invalid.append("aa_ci:" + label)
                elif summary["log_mean"] > point_limit:
                    reject.append("cell_point:" + label)
                if (
                    comparison == "C/L" and scope in REPAIR_SCOPES and
                    summary["upper95_log"] > repair_limit
                ):
                    reject.append("repair_cell_upper95:" + label)

    groups: List[Tuple[str, Tuple[Tuple[int, int], ...]]] = [
        ("equal-cell", CELLS),
        ("width-B64", tuple(cell for cell in CELLS if cell[1] == 64)),
        ("width-B1280", tuple(cell for cell in CELLS if cell[1] == 1280)),
    ]
    groups.extend((
        "band-" + band,
        tuple(cell for cell in CELLS if band_for(cell[0]) == band),
    ) for band in ("small", "medium", "large"))
    aggregates: Dict[str, Any] = {}
    for scope in SCOPES:
        scope_out: Dict[str, Any] = {}
        for name, selected in groups:
            summary = sample_summary(aggregate_logs(logs, scope, selected))
            scope_out[name] = summary
            if scope in SYSTEMATIC_SCOPES:
                if not summary["upper95_log"] < systematic_limit:
                    reject.append("systematic_group_upper95:C/L:{}:{}".format(
                        scope, name))
            elif summary["upper95_log"] > repair_limit:
                reject.append("repair_group_upper95:C/L:{}:{}".format(
                    scope, name))
        aggregates[scope] = scope_out
    return ({
        "aggregates": {"C/L": aggregates},
        "cells": cells_out,
        "gate_policy": {
            "aa_two_sided_ratio_strict": 1.02,
            "repair_cell_and_group_upper_ratio_inclusive": 1.02,
            "systematic_group_upper_ratio_strict": 0.99,
            "systematic_point_upper_ratio_inclusive": 1.02,
        },
        "student_t_975_df11": T11_975,
    }, invalid, reject)


def parse_transcript(
    raw: bytes, expected: Mapping[str, str], runtime_maps_sha256: str,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], List[str], List[str]]:
    if (
        not raw or len(raw) > MAX_RAW_BYTES or not raw.endswith(b"\n") or
        b"\r" in raw
    ):
        fail("raw transcript is not bounded LF-terminated JSONL")
    lines = raw.splitlines(keepends=True)
    if len(lines) != EXPECTED_RECORD_COUNT:
        fail("raw transcript has {} records, expected {}".format(
            len(lines), EXPECTED_RECORD_COUNT))
    records = [
        parse_canonical_line(line, "raw line {}".format(index))
        for index, line in enumerate(lines)
    ]
    config = records[0]
    oracles = validate_config(config, expected)
    logs: Dict[Tuple[str, int, int, str], List[float]] = {}
    sequence = 0
    measured = 0
    encode_calls = 0
    record_index = 1
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison, _, _ in COMPARISONS:
                    contrast, count, calls = validate_panel(
                        records[record_index], replicate, k, block_bytes,
                        scope, comparison, sequence, oracles,
                        runtime_maps_sha256)
                    logs.setdefault(
                        (comparison, k, block_bytes, scope), []).append(contrast)
                    measured += count
                    encode_calls += calls
                    sequence += 1
                    record_index += 1
    terminal = records[-1]
    exact_keys(terminal, TERMINAL_KEYS, "terminal")
    exact_string(terminal["schema"], TERMINAL_SCHEMA, "terminal.schema")
    exact_string(terminal["status"], "complete", "terminal.status")
    exact_int(terminal["panel_count"], EXPECTED_PANEL_COUNT,
              EXPECTED_PANEL_COUNT, "terminal.panel_count")
    exact_int(terminal["record_count"], EXPECTED_RECORD_COUNT,
              EXPECTED_RECORD_COUNT, "terminal.record_count")
    exact_int(terminal["measured_invocation_count"], measured, measured,
              "terminal.measured_invocation_count")
    exact_int(terminal["warmup_invocation_count"], EXPECTED_WARMUP_INVOCATIONS,
              EXPECTED_WARMUP_INVOCATIONS, "terminal.warmup_invocation_count")
    exact_int(terminal["encode_call_count"], encode_calls, encode_calls,
              "terminal.encode_call_count")
    if measured != EXPECTED_MEASURED_INVOCATIONS:
        fail("derived measured invocation total differs")
    if encode_calls != EXPECTED_ENCODE_CALLS:
        fail("derived encode-call total differs")
    statistics, invalid, reject = adjudicate(logs)
    return config, terminal, statistics, invalid, reject


def run_checked(
    argv: Sequence[str], maximum: int = MAX_PROBE_BYTES,
    env: Optional[Mapping[str, str]] = None,
    require_empty_stderr: bool = False,
) -> bytes:
    try:
        completed = subprocess.run(
            list(argv), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=30.0, check=False,
            env=dict(env) if env is not None else {
                "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
                "TZ": "UTC",
            },
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("probe failed for {}: {}".format(argv[0], exc))
    if (
        completed.returncode != 0 or len(completed.stdout) > maximum or
        len(completed.stderr) > maximum or
        (require_empty_stderr and completed.stderr)
    ):
        fail("probe failed for {} with exit {}: {}".format(
            argv[0], completed.returncode,
            completed.stderr[:4096].decode("utf-8", "replace")))
    return completed.stdout


def exact_absolute_file(text: str, where: str) -> Path:
    path = Path(text)
    if not path.is_absolute() or path.resolve() != path or not path.is_file():
        fail("{} is not an exact absolute regular file".format(where))
    info = path.stat()
    if stat.S_ISLNK(info.st_mode):
        fail("{} is a symbolic link".format(where))
    return path


def exact_absolute_directory(text: str, where: str) -> Path:
    path = Path(text)
    if not path.is_absolute() or path.resolve() != path or not path.is_dir():
        fail("{} is not an exact absolute directory".format(where))
    return path


def stat_receipt(path: Path) -> Dict[str, Any]:
    info = os.stat(str(path), follow_symlinks=False)
    return {
        "device": info.st_dev,
        "inode": info.st_ino,
        "mode": stat.S_IMODE(info.st_mode),
        "mtime_ns": info.st_mtime_ns,
        "path": str(path),
        "size": info.st_size,
    }


def artifact_receipt(path: Path) -> Dict[str, Any]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(str(path), flags)
    except OSError as exc:
        fail("cannot open artifact {}: {}".format(path, exc))
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink < 1 or
            before.st_size < 0 or before.st_size > MAX_RAW_BYTES
        ):
            fail("{} is not a bounded regular artifact".format(path))
        digest = hashlib.sha256()
        total = 0
        while True:
            chunk = os.read(fd, 1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_RAW_BYTES:
                fail("artifact {} exceeds its byte cap".format(path))
            digest.update(chunk)
        after = os.fstat(fd)
        stable_fields = (
            "st_dev", "st_ino", "st_mode", "st_nlink", "st_size",
            "st_mtime_ns",
        )
        if total != before.st_size or any(
            getattr(before, field) != getattr(after, field)
            for field in stable_fields
        ):
            fail("artifact {} changed while it was hashed".format(path))
        return {
            "device": before.st_dev,
            "inode": before.st_ino,
            "mode": stat.S_IMODE(before.st_mode),
            "mtime_ns": before.st_mtime_ns,
            "path": str(path),
            "sha256": digest.hexdigest(),
            "size": before.st_size,
        }
    finally:
        os.close(fd)


def snapshot_r0() -> Dict[str, Any]:
    entries: List[Dict[str, Any]] = []
    for root in R0_ROOTS:
        try:
            root_info = os.stat(str(root), follow_symlinks=False)
        except FileNotFoundError:
            entries.append({"exists": False, "path": str(root)})
            continue
        except OSError as exc:
            entries.append({
                "accessible": False,
                "errno": exc.errno,
                "exists": None,
                "path": str(root),
            })
            continue
        paths = [root]
        if stat.S_ISDIR(root_info.st_mode):
            try:
                paths.extend(sorted(root.rglob("*"), key=lambda item: str(item)))
            except OSError as exc:
                entries.append({
                    "accessible": False,
                    "errno": exc.errno,
                    "exists": True,
                    "path": str(root) + "/<children>",
                })
        for path in paths:
            try:
                info = os.stat(str(path), follow_symlinks=False)
            except OSError as exc:
                fail("cannot snapshot r0 {}: {}".format(path, exc))
            entry = {
                "device": info.st_dev,
                "exists": True,
                "inode": info.st_ino,
                "mode": stat.S_IMODE(info.st_mode),
                "mtime_ns": info.st_mtime_ns,
                "path": str(path),
                "size": info.st_size,
                "type": (
                    "file" if stat.S_ISREG(info.st_mode) else
                    "directory" if stat.S_ISDIR(info.st_mode) else
                    "symlink" if stat.S_ISLNK(info.st_mode) else "other"
                ),
            }
            if stat.S_ISREG(info.st_mode) and os.access(str(path), os.R_OK):
                entry["sha256"] = file_sha256(path, maximum=MAX_RAW_BYTES)
            entries.append(entry)
    value = {"entries": entries}
    value["snapshot_sha256"] = sha256_bytes(canonical_bytes(entries))
    return value


def parse_cache(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    try:
        data = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        fail("cannot read CMakeCache.txt: {}".format(exc))
    if len(data.encode("utf-8")) > MAX_PROBE_BYTES:
        fail("CMakeCache.txt exceeds its byte cap")
    for line in data.splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        result[key] = value
    return result


def compiler_receipt(path: Path) -> Dict[str, Any]:
    version_bytes = run_checked([str(path), "--version"], maximum=65536)
    short_bytes = run_checked(
        [str(path), "-dumpfullversion", "-dumpversion"], maximum=1024)
    try:
        version_text = version_bytes.decode("utf-8")
        version = short_bytes.decode("ascii").strip()
    except UnicodeDecodeError:
        fail("compiler version is not UTF-8")
    if (
        not version_text or not version_text.endswith("\n") or
        "gcc" not in version_text.lower() and "g++" not in version_text.lower()
    ):
        fail("compiler is not GCC")
    if re.fullmatch(r"13(?:\.[0-9]+){1,2}", version) is None:
        fail("compiler is not GCC 13")
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "version": version,
        "version_text": version_text,
        "version_text_sha256": sha256_bytes(version_bytes),
        "version_sha256": sha256_bytes(short_bytes),
    }


def interpreter_receipt() -> Dict[str, Any]:
    invoked_path = Path(sys.executable)
    if not invoked_path.is_absolute():
        fail("controller Python interpreter path is not absolute")
    try:
        invoked_resolved = invoked_path.resolve(strict=True)
        proc_resolved = Path("/proc/self/exe").resolve(strict=True)
    except OSError as exc:
        fail("cannot resolve the running Python image: {}".format(exc))
    if invoked_resolved != proc_resolved:
        fail("sys.executable differs from the running Python image")
    path = exact_absolute_file(str(proc_resolved),
                               "running controller Python image")
    version_info = [
        sys.version_info.major,
        sys.version_info.minor,
        sys.version_info.micro,
        sys.version_info.releaselevel,
        sys.version_info.serial,
    ]
    return {
        "artifact": artifact_receipt(path),
        "dont_write_bytecode": sys.flags.dont_write_bytecode == 1,
        "ignore_environment": sys.flags.ignore_environment == 1,
        "implementation": sys.implementation.name,
        "invoked_path": str(invoked_path),
        "isolated": sys.flags.isolated == 1,
        "no_site": sys.flags.no_site == 1,
        "no_user_site": sys.flags.no_user_site == 1,
        "path": str(path),
        "version": sys.version,
        "version_info": version_info,
    }


def parse_posix_nm(data: bytes, where: str) -> List[Tuple[str, str]]:
    try:
        text = data.decode("utf-8", "strict")
    except UnicodeDecodeError:
        fail("{} is not UTF-8".format(where))
    result: List[Tuple[str, str]] = []
    for line_number, line in enumerate(text.splitlines(), 1):
        fields = line.split()
        if len(fields) < 2 or len(fields[1]) != 1:
            fail("{} line {} is malformed".format(where, line_number))
        result.append((fields[0], fields[1]))
    return result


def dynamic_receipt(worker: Path, library: Path) -> Dict[str, Any]:
    nm_path = exact_absolute_file("/usr/bin/x86_64-linux-gnu-nm",
                                  "system nm")
    # /usr/bin/nm is the fixed invocation path, while the resolved artifact is
    # recorded separately so alternatives/symlink changes cannot be hidden.
    if Path("/usr/bin/nm").resolve() != nm_path:
        fail("/usr/bin/nm does not resolve to the frozen nm artifact")
    nm_version = run_checked(
        ["/usr/bin/nm", "--version"], maximum=65536,
        require_empty_stderr=True)
    readelf = run_checked(
        ["/usr/bin/readelf", "-d", str(worker)],
        require_empty_stderr=True)
    ldd_env = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        "LD_LIBRARY_PATH": str(library.parent),
    }
    ldd = run_checked(
        ["/usr/bin/ldd", str(worker)], env=ldd_env,
        require_empty_stderr=True)
    library_nm = run_checked([
        "/usr/bin/nm", "-D", "--format=posix", "--defined-only",
        str(library),
    ], require_empty_stderr=True)
    library_nm_text = library_nm.decode("utf-8", "strict")
    library_rows = parse_posix_nm(library_nm, "library dynamic nm")
    dynamic_symbols = {name.split("@", 1)[0] for name, _ in library_rows}
    missing = sorted(REQUIRED_DYNAMIC_SYMBOLS - dynamic_symbols)
    if missing:
        fail("library lacks required public dynamic symbols: {}".format(missing))
    forbidden = sorted(
        symbol for symbol in dynamic_symbols
        if any(fragment in symbol for fragment in FORBIDDEN_SYMBOL_FRAGMENTS)
    )
    if forbidden:
        fail("library exports forbidden internal/test symbols: {}".format(
            forbidden[:32]))

    worker_dynamic_nm = run_checked([
        "/usr/bin/nm", "-D", "--format=posix", str(worker),
    ], require_empty_stderr=True)
    worker_dynamic_text = worker_dynamic_nm.decode("utf-8", "strict")
    worker_dynamic_rows = parse_posix_nm(
        worker_dynamic_nm, "worker dynamic nm")
    wirehair_imports = [
        (name, symbol_type) for name, symbol_type in worker_dynamic_rows
        if name.split("@", 1)[0].startswith("wirehair")
    ]
    expected_imports = sorted(
        symbol + "@" + WIREHAIR_SYMBOL_VERSION
        for symbol in REQUIRED_DYNAMIC_SYMBOLS)
    actual_imports = sorted(name for name, _ in wirehair_imports)
    if (
        len(wirehair_imports) != len(REQUIRED_DYNAMIC_SYMBOLS) or
        len(set(actual_imports)) != len(actual_imports) or
        actual_imports != expected_imports or
        any(symbol_type != "U" for _, symbol_type in wirehair_imports) or
        any("@@" in name for name in actual_imports)
    ):
        fail("worker public Wirehair import roster/version/type differs")

    worker_full_nm = run_checked([
        "/usr/bin/nm", "-C", "--format=posix", str(worker),
    ], require_empty_stderr=True)
    try:
        worker_full_text = worker_full_nm.decode("utf-8", "strict")
    except UnicodeDecodeError:
        fail("worker full nm is not UTF-8")
    forbidden_worker = [
        fragment for fragment in WORKER_FORBIDDEN_SYMBOL_FRAGMENTS
        if fragment in worker_full_text
    ]
    if forbidden_worker:
        fail("worker contains forbidden codec/test-hook symbols: {}".format(
            forbidden_worker))

    ldd_text = ldd.decode("utf-8", "strict")
    matches = [
        line for line in ldd_text.splitlines()
        if "libwirehair" in line
    ]
    if len(matches) != 1:
        fail("ldd does not bind exactly to the intended libwirehair")
    binding_match = re.fullmatch(
        r"\s*libwirehair\.so[^\s]*\s+=>\s+(\S+)\s+\(0x[0-9a-fA-F]+\)\s*",
        matches[0])
    if binding_match is None:
        fail("ldd libwirehair binding is not in the required resolved form")
    loaded_path = Path(binding_match.group(1))
    if not loaded_path.is_absolute() or not loaded_path.is_file():
        fail("ldd libwirehair target is not an absolute regular file")
    try:
        loaded_info = loaded_path.stat()
        intended_info = library.stat()
        loaded_resolved = loaded_path.resolve()
    except OSError as exc:
        fail("cannot stat ldd libwirehair target: {}".format(exc))
    if (
        loaded_resolved != library or loaded_info.st_dev != intended_info.st_dev or
        loaded_info.st_ino != intended_info.st_ino
    ):
        fail("ldd does not bind to the intended libwirehair file identity")
    readelf_text = readelf.decode("utf-8", "strict")
    if "Shared library: [libwirehair" not in readelf_text:
        fail("worker does not declare a libwirehair DT_NEEDED entry")
    normalized_ldd_text = re.sub(
        r"\(0x[0-9a-fA-F]+\)", "(0xADDR)", ldd_text)
    return {
        "forbidden_symbol_fragments": list(FORBIDDEN_SYMBOL_FRAGMENTS),
        "ldd_normalized_sha256": sha256_bytes(
            normalized_ldd_text.encode("utf-8")),
        "ldd_library_path": str(loaded_path),
        "ldd_library_resolved_path": str(loaded_resolved),
        "ldd_library_stat": stat_receipt(loaded_path),
        "ldd_normalized_text": normalized_ldd_text,
        "library_nm_argv": [
            "/usr/bin/nm", "-D", "--format=posix", "--defined-only",
            str(library),
        ],
        "library_symbol_table_text": library_nm_text,
        "nm_artifact": artifact_receipt(nm_path),
        "nm_invoked_path": "/usr/bin/nm",
        "nm_version_sha256": sha256_bytes(nm_version),
        "nm_version_text": nm_version.decode("utf-8", "strict"),
        "readelf_dynamic_sha256": sha256_bytes(readelf),
        "readelf_dynamic_text": readelf_text,
        "required_dynamic_symbols": sorted(REQUIRED_DYNAMIC_SYMBOLS),
        "required_dynamic_symbols_verified": True,
        "symbol_table_sha256": sha256_bytes(library_nm),
        "worker_forbidden_symbol_fragments":
            list(WORKER_FORBIDDEN_SYMBOL_FRAGMENTS),
        "worker_full_symbol_table_sha256": sha256_bytes(worker_full_nm),
        "worker_full_symbol_table_text": worker_full_text,
        "worker_full_nm_argv": [
            "/usr/bin/nm", "-C", "--format=posix", str(worker),
        ],
        "worker_import_symbol_version": WIREHAIR_SYMBOL_VERSION,
        "worker_import_symbols": actual_imports,
        "worker_import_symbols_verified": True,
        "worker_symbol_table_sha256": sha256_bytes(worker_dynamic_nm),
        "worker_symbol_table_text": worker_dynamic_text,
        "worker_nm_argv": [
            "/usr/bin/nm", "-D", "--format=posix", str(worker),
        ],
    }


def git_receipt(source_root: Path, expected_commit: str) -> Dict[str, Any]:
    git_invoked = Path("/usr/bin/git")
    try:
        git_path = exact_absolute_file(
            str(git_invoked.resolve(strict=True)), "system Git")
    except OSError as exc:
        fail("cannot resolve system Git: {}".format(exc))
    git_env = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_OPTIONAL_LOCKS": "0", "GIT_TERMINAL_PROMPT": "0",
    }
    prefix = [
        "/usr/bin/git", "--no-replace-objects",
        "-c", "core.attributesFile=/dev/null",
        "-c", "core.excludesFile=/dev/null",
        "-C", str(source_root),
    ]
    version_data = run_checked(
        ["/usr/bin/git", "--version"], maximum=65536, env=git_env,
        require_empty_stderr=True)
    head = run_checked(
        prefix + ["rev-parse", "HEAD"], maximum=1024, env=git_env,
        require_empty_stderr=True).decode(
        "ascii").strip()
    tree = run_checked(
        prefix + ["rev-parse", "HEAD^{tree}"], maximum=1024, env=git_env,
        require_empty_stderr=True).decode(
        "ascii").strip()
    upstream = run_checked(
        prefix + ["rev-parse", "@{upstream}"], maximum=1024,
        env=git_env, require_empty_stderr=True).decode("ascii").strip()
    # The frozen contract requires a clean tracked tree.  Deliberately ignore
    # unrelated protected untracked artifacts in the shared worktree.
    status_data = run_checked(
        prefix + ["status", "--porcelain=v1", "--untracked-files=no"],
        env=git_env, require_empty_stderr=True)
    lower_commit(head, "git HEAD")
    lower_commit(tree, "git tree")
    lower_commit(upstream, "git upstream")
    if head != expected_commit or upstream != head or status_data:
        fail("source is not the exact clean pushed commit")
    tracked_input_blobs: Dict[str, str] = {}
    for relative in R1_TRACKED_INPUTS:
        committed_blob = run_checked(
            prefix + ["rev-parse", "{}:{}".format(head, relative)],
            maximum=1024, env=git_env,
            require_empty_stderr=True).decode("ascii").strip()
        worktree_blob = run_checked(
            prefix + ["hash-object", "--no-filters", "--", relative],
            maximum=1024, env=git_env,
            require_empty_stderr=True).decode("ascii").strip()
        lower_commit(committed_blob, "committed input blob " + relative)
        lower_commit(worktree_blob, "worktree input blob " + relative)
        if committed_blob != worktree_blob:
            fail("required r1 input is not exact at HEAD: " + relative)
        tracked_input_blobs[relative] = committed_blob
    return {
        "commit": head,
        "git_artifact": artifact_receipt(git_path),
        "git_invoked_path": str(git_invoked),
        "git_version_sha256": sha256_bytes(version_data),
        "git_version_text": version_data.decode("utf-8", "strict"),
        "required_tracked_inputs": tracked_input_blobs,
        "tracked_clean": True,
        "tree": tree,
        "upstream_commit": upstream,
    }


def path_is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def build_receipt(
    source_root: Path, build_root: Path, compiler: Mapping[str, Any],
    interpreter: Mapping[str, Any], worker: Path, library: Path,
    source_commit: str,
) -> Dict[str, Any]:
    cache = parse_cache(build_root / "CMakeCache.txt")
    required = {
        "BUILD_SHARED_LIBS": "ON",
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_GENERATOR": "Ninja",
        "MARCH_NATIVE": "ON",
        "WH_LTO": "OFF",
        "WH_PGO_MODE": "OFF",
    }
    actual = {key: cache.get(key) for key in required}
    if actual != required:
        fail("build configuration differs: {}".format(actual))
    configured_source = cache.get("CMAKE_HOME_DIRECTORY")
    configured_build = cache.get("CMAKE_CACHEFILE_DIR")
    if (
        configured_source is None or Path(configured_source).resolve() != source_root or
        configured_build is None or Path(configured_build).resolve() != build_root
    ):
        fail("CMake source/build roots differ from the exact roots")
    if not path_is_within(worker, build_root) or not path_is_within(library, build_root):
        fail("worker or library is outside the exact build root")
    configured_compiler = cache.get("CMAKE_CXX_COMPILER")
    if configured_compiler is None:
        fail("CMake compiler is absent")
    configured_compiler_path = Path(configured_compiler)
    if configured_compiler_path.is_absolute():
        configured_compiler_resolved = configured_compiler_path.resolve()
    elif "/" not in configured_compiler:
        candidates = {
            path.resolve() for path in (
                Path("/usr/bin") / configured_compiler,
                Path("/bin") / configured_compiler,
            ) if path.is_file()
        }
        if len(candidates) != 1:
            fail("relative CMake compiler does not resolve uniquely")
        configured_compiler_resolved = candidates.pop()
    else:
        fail("relative CMake compiler contains a path component")
    if configured_compiler_resolved != Path(compiler["path"]):
        fail("CMake compiler differs from the exact compiler")
    configured_interpreter = cache.get("_Python3_EXECUTABLE")
    if configured_interpreter is None:
        fail("CMake Python interpreter is absent")
    configured_interpreter_path = Path(configured_interpreter)
    if not configured_interpreter_path.is_absolute():
        fail("CMake Python interpreter is not absolute")
    try:
        configured_interpreter_resolved = configured_interpreter_path.resolve(
            strict=True)
    except OSError as exc:
        fail("cannot resolve CMake Python interpreter: {}".format(exc))
    if configured_interpreter_resolved != Path(interpreter["path"]):
        fail("CMake Python interpreter differs from the running controller")
    probe_env = {"LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin"}
    dry_run = run_checked([
        "/usr/bin/ninja", "-C", str(build_root), "-n",
        "wirehair_wh2_public_borrowed_current_vs_wh1_r1",
    ], env=probe_env)
    dry_run_text = dry_run.decode("utf-8", "strict")
    dry_run_lines = dry_run_text.splitlines()
    if not (
        dry_run_lines and dry_run_lines[-1] == "ninja: no work to do." and
        len(dry_run_lines) in (1, 2) and
        (len(dry_run_lines) == 1 or
         dry_run_lines[0].startswith("ninja: Entering directory `"))
    ):
        fail("exact worker target is not fully up to date")
    commands = run_checked([
        "/usr/bin/ninja", "-C", str(build_root), "-t", "commands",
        "wirehair_wh2_public_borrowed_current_vs_wh1_r1",
    ], env=probe_env)
    command_text = commands.decode("utf-8", "strict")
    parsed_rows: List[Tuple[str, List[str]]] = []
    try:
        for line in command_text.splitlines():
            parsed_rows.append((line, shlex.split(line)))
    except ValueError as exc:
        fail("cannot parse exact Ninja commands: {}".format(exc))

    expected_definitions = [
        "-DWIREHAIR_DLL=1",
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_ID="GNU"',
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION="{}"'.format(
            compiler["version"]),
        '-DWIREHAIR_WH2_CXX_COMPILER_PATH="{}"'.format(compiler["path"]),
        '-DWIREHAIR_WH2_CXX_COMPILER_SHA256="{}"'.format(
            compiler["sha256"]),
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(source_commit),
    ]
    target_prefix = (
        "CMakeFiles/wirehair_wh2_public_borrowed_current_vs_wh1_r1.dir/")
    expected_objects = {
        relative: target_prefix + relative + ".o"
        for relative in R1_TARGET_SOURCES
    }
    compile_commands: Dict[str, str] = {}
    compile_command_hashes: Dict[str, str] = {}
    compiler_invocation: Optional[str] = None
    for relative, expected_object in expected_objects.items():
        source_path = str(source_root / relative)
        matches: List[Tuple[str, List[str]]] = []
        for line, tokens in parsed_rows:
            if "-c" not in tokens:
                continue
            indices = [index for index, token in enumerate(tokens)
                       if token == "-c"]
            if len(indices) == 1 and indices[0] + 1 < len(tokens) and \
                    tokens[indices[0] + 1] == source_path:
                matches.append((line, tokens))
        if len(matches) != 1:
            fail("cannot identify the exact compile command for {}".format(
                relative))
        line, tokens = matches[0]
        try:
            compiler_token = Path(tokens[0]).resolve(strict=True)
        except (OSError, IndexError) as exc:
            fail("cannot resolve worker compile command: {}".format(exc))
        if compiler_token != Path(compiler["path"]):
            fail("worker compile command uses a different compiler")
        if compiler_invocation is None:
            compiler_invocation = tokens[0]
        elif tokens[0] != compiler_invocation:
            fail("worker compile compiler invocation differs")
        expected_tail = expected_definitions + [
            "-I" + str(source_root / "include"),
            "-O3", "-DNDEBUG", "-std=gnu++11", "-Wall", "-Wextra",
            "-march=native", "-MD", "-MT", expected_object,
            "-MF", expected_object + ".d", "-o", expected_object,
            "-c", source_path,
        ]
        if tokens[1:] != expected_tail:
            fail("worker compile command differs for {}".format(relative))
        compile_commands[relative] = line
        compile_command_hashes[relative] = sha256_bytes(line.encode("utf-8"))

    library_definitions = [
        "-DWIREHAIR_BUILDING=1", "-DWIREHAIR_DLL=1",
        "-Dwirehair_EXPORTS",
    ]
    library_prefix = "CMakeFiles/wirehair.dir/"
    expected_library_objects = {
        relative: library_prefix + relative + ".o"
        for relative in R1_LIBRARY_SOURCES
    }
    library_compile_commands: Dict[str, str] = {}
    library_compile_hashes: Dict[str, str] = {}
    for relative, expected_object in expected_library_objects.items():
        source_path = str(source_root / relative)
        matches = []
        for line, tokens in parsed_rows:
            indices = [index for index, token in enumerate(tokens)
                       if token == "-c"]
            if len(indices) == 1 and indices[0] + 1 < len(tokens) and \
                    tokens[indices[0] + 1] == source_path:
                matches.append((line, tokens))
        if len(matches) != 1:
            fail("cannot identify exact library compile command for {}".format(
                relative))
        line, tokens = matches[0]
        try:
            compiler_token = Path(tokens[0]).resolve(strict=True)
        except (OSError, IndexError) as exc:
            fail("cannot resolve library compile command: {}".format(exc))
        expected_tail = library_definitions + [
            "-I" + str(source_root / "include"),
            "-O3", "-DNDEBUG", "-std=gnu++11", "-fPIC", "-Wall",
            "-Wextra", "-march=native", "-MD", "-MT", expected_object,
            "-MF", expected_object + ".d", "-o", expected_object,
            "-c", source_path,
        ]
        if (
            compiler_token != Path(compiler["path"]) or
            tokens[0] != compiler_invocation or tokens[1:] != expected_tail
        ):
            fail("library compile command differs for {}".format(relative))
        library_compile_commands[relative] = line
        library_compile_hashes[relative] = sha256_bytes(line.encode("utf-8"))

    if library.name != "libwirehair.so.2.0.0":
        fail("library filename differs from the exact shared object")
    library_link_matches: List[Tuple[str, List[str]]] = []
    for line, tokens in parsed_rows:
        if "-c" in tokens:
            continue
        output_indices = [index for index, token in enumerate(tokens)
                          if token == "-o"]
        if len(output_indices) == 1 and output_indices[0] + 1 < len(tokens) and \
                tokens[output_indices[0] + 1] == library.name:
            library_link_matches.append((line, tokens))
    if len(library_link_matches) != 1:
        fail("cannot identify the exact library link command")
    library_link_command, library_link_tokens = library_link_matches[0]
    try:
        library_link_compiler = Path(library_link_tokens[2]).resolve(strict=True)
    except (OSError, IndexError) as exc:
        fail("cannot resolve library link command: {}".format(exc))
    library_object_roster = [
        expected_library_objects[item] for item in R1_LIBRARY_SOURCES
    ]
    expected_library_link_tokens = [
        ":", "&&", library_link_tokens[2] if len(library_link_tokens) > 2 else "",
        "-fPIC", "-O3", "-DNDEBUG",
        "-Wl,--version-script={}".format(source_root / "abi" / "wirehair.map"),
        "-shared", "-Wl,-soname,libwirehair.so.2", "-o", library.name,
    ] + library_object_roster + ["-lm", "&&", ":"]
    if library_link_compiler != Path(compiler["path"]) or \
            library_link_tokens[2] != compiler_invocation or \
            library_link_tokens != expected_library_link_tokens:
        fail("library link command/input roster differs")

    link_matches: List[Tuple[str, List[str]]] = []
    for line, tokens in parsed_rows:
        if "-c" in tokens:
            continue
        output_indices = [index for index, token in enumerate(tokens)
                          if token == "-o"]
        if len(output_indices) == 1 and output_indices[0] + 1 < len(tokens) and \
                tokens[output_indices[0] + 1] == worker.name:
            link_matches.append((line, tokens))
    if len(link_matches) != 1:
        fail("cannot identify the exact worker link command")
    link_command, link_tokens = link_matches[0]
    compiler_positions = []
    for index, token in enumerate(link_tokens):
        if not token.startswith("/"):
            continue
        try:
            if Path(token).resolve(strict=True) == Path(compiler["path"]):
                compiler_positions.append(index)
        except OSError:
            continue
    expected_object_roster = [expected_objects[item] for item in R1_TARGET_SOURCES]
    expected_link_tokens = [
        ":", "&&", link_tokens[2] if len(link_tokens) > 2 else "",
        "-O3", "-DNDEBUG",
    ] + expected_object_roster + [
        "-o", worker.name, "-Wl,-rpath,{}".format(build_root),
        library.name, "-ldl", "&&", ":",
    ]
    if len(compiler_positions) != 1 or compiler_positions != [2] or \
            link_tokens[2] != compiler_invocation or \
            link_tokens != expected_link_tokens:
        fail("worker link command/input roster differs")
    return {
        "cache": actual,
        "cmake_build_root": configured_build,
        "cmake_compiler": configured_compiler,
        "cmake_compiler_resolved": str(configured_compiler_resolved),
        "cmake_interpreter": configured_interpreter,
        "cmake_interpreter_resolved": str(configured_interpreter_resolved),
        "cmake_source_root": configured_source,
        "compiler_invocation": compiler_invocation,
        "compile_command_sha256": compile_command_hashes,
        "compile_commands": compile_commands,
        "compile_definitions": expected_definitions,
        "dry_run_sha256": sha256_bytes(dry_run),
        "dry_run_text": dry_run_text,
        "library_compile_command_sha256": library_compile_hashes,
        "library_compile_commands": library_compile_commands,
        "library_compile_definitions": library_definitions,
        "library_link_command": library_link_command,
        "library_link_command_sha256": sha256_bytes(
            library_link_command.encode("utf-8")),
        "library_link_object_roster": library_object_roster,
        "library_sources": list(R1_LIBRARY_SOURCES),
        "link_command": link_command,
        "link_command_sha256": sha256_bytes(link_command.encode("utf-8")),
        "link_object_roster": expected_object_roster,
        "ninja_commands_sha256": sha256_bytes(commands),
        "ninja_commands_text": command_text,
        "target_sources": list(R1_TARGET_SOURCES),
    }


def pin_controller() -> Dict[str, Any]:
    try:
        before = sorted(os.sched_getaffinity(0))
    except OSError as exc:
        fail("cannot read inherited controller affinity: {}".format(exc))
    required = {CONTROLLER_CPU, TARGET_CPU, TARGET_SIBLING}
    if not required.issubset(set(before)):
        fail("controller inherited affinity lacks a required CPU")
    if CONTROLLER_CPU in (TARGET_CPU, TARGET_SIBLING):
        fail("controller CPU collides with the timed core or sibling")
    try:
        os.sched_setaffinity(0, {CONTROLLER_CPU})
        after = sorted(os.sched_getaffinity(0))
    except OSError as exc:
        fail("cannot pin controller to its fixed CPU: {}".format(exc))
    if after != [CONTROLLER_CPU]:
        fail("controller singleton affinity was not realized")
    return {
        "affinity_after": after,
        "affinity_before": before,
        "controller_cpu": CONTROLLER_CPU,
        "sibling_available_before": True,
        "singleton_verified": True,
        "target_available_before": True,
    }


def parse_cpu_list(text: str, where: str) -> List[int]:
    sibling_values: List[int] = []
    try:
        for part in text.split(","):
            if not part:
                fail("{} contains an empty component".format(where))
            if "-" in part:
                bounds = part.split("-", 1)
                if len(bounds) != 2:
                    fail("{} contains a malformed range".format(where))
                first, last = (int(item) for item in bounds)
                if first < 0 or last < first or last > 1048575:
                    fail("{} contains an invalid range".format(where))
                sibling_values.extend(range(first, last + 1))
            else:
                value = int(part)
                if value < 0 or value > 1048575:
                    fail("{} contains an invalid CPU".format(where))
                sibling_values.append(value)
    except ValueError:
        fail("{} is not a valid CPU list".format(where))
    if len(sibling_values) != len(set(sibling_values)):
        fail("{} contains duplicate CPUs".format(where))
    return sibling_values


def cpu_receipt() -> Dict[str, Any]:
    try:
        available = sorted(os.sched_getaffinity(0))
    except OSError as exc:
        fail("cannot read controller affinity: {}".format(exc))
    if available != [CONTROLLER_CPU]:
        fail("controller escaped its fixed singleton affinity")
    sibling_path = Path(
        "/sys/devices/system/cpu/cpu{}/topology/thread_siblings_list".format(
            TARGET_CPU))
    try:
        sibling_text = sibling_path.read_text(encoding="ascii").strip()
    except (OSError, UnicodeDecodeError) as exc:
        fail("cannot read target sibling topology: {}".format(exc))
    sibling_values = parse_cpu_list(sibling_text, "target sibling topology")
    if sorted(sibling_values) != sorted((TARGET_CPU, TARGET_SIBLING)):
        fail("target CPU sibling topology differs")
    governor_path = Path(
        "/sys/devices/system/cpu/cpu{}/cpufreq/scaling_governor".format(TARGET_CPU))
    try:
        governor = (
            governor_path.read_text(encoding="ascii").strip()
            if governor_path.is_file() else "unavailable"
        )
        load = Path("/proc/loadavg").read_text(encoding="ascii").strip()
        stat_text = Path("/proc/stat").read_text(encoding="ascii")
    except (OSError, UnicodeDecodeError) as exc:
        fail("cannot read CPU health inputs: {}".format(exc))
    matches = [line for line in stat_text.splitlines()
               if line.startswith("cpu{} ".format(TARGET_SIBLING))]
    if len(matches) != 1:
        fail("sibling CPU activity row is unavailable")
    tick_parts = matches[0].split()
    if len(tick_parts) < 9:
        fail("sibling CPU activity row is incomplete")
    try:
        tick_values = [int(value) for value in tick_parts[1:9]]
    except ValueError:
        fail("sibling CPU activity row is not numeric")
    if any(value < 0 for value in tick_values):
        fail("sibling CPU activity row contains a negative value")
    non_idle_ticks = sum(tick_values[index] for index in (0, 1, 2, 5, 6, 7))
    thermal_paths = set(Path("/sys/class/thermal").glob("thermal_zone*/temp"))
    thermal_paths.update(
        Path("/sys/class/hwmon").glob("hwmon*/temp*_input"))
    thermal: List[Dict[str, Any]] = []
    for path in sorted(thermal_paths, key=lambda item: str(item)):
        try:
            raw_value = path.read_text(encoding="ascii").strip()
        except (OSError, UnicodeDecodeError):
            continue
        try:
            value = int(raw_value)
        except ValueError:
            fail("thermal path {} is not numeric".format(path))
        if value < 0 or value > THERMAL_ENDPOINT_MAX_MILLIC:
            fail("thermal path {} exceeds the frozen cap".format(path))
        thermal.append({"path": str(path), "value_millic": value})
    if not thermal:
        fail("no readable thermal endpoint is available")
    edac: List[Dict[str, Any]] = []
    for counter in ("ce_count", "ue_count"):
        for path in sorted(
            Path("/sys/devices/system/edac/mc").glob("mc*/" + counter)
        ):
            try:
                raw_value = path.read_text(encoding="ascii").strip()
            except OSError:
                continue
            try:
                value = int(raw_value)
            except ValueError:
                fail("EDAC counter {} is not numeric".format(path))
            if value < 0:
                fail("EDAC counter {} is negative".format(path))
            edac.append({"counter": counter, "path": str(path), "value": value})
    return {
        "available_cpu_count": len(available),
        "controller_affinity": available,
        "controller_cpu": CONTROLLER_CPU,
        "edac": edac,
        "governor": governor,
        "loadavg": load,
        "sibling_cpu": TARGET_SIBLING,
        "sibling_non_idle_ticks": non_idle_ticks,
        "sibling_proc_stat": matches[0],
        "sibling_tick_fields": tick_values,
        "target_cpu": TARGET_CPU,
        "thermal": thermal,
        "thermal_endpoint_max_millic": THERMAL_ENDPOINT_MAX_MILLIC,
        "thread_siblings_list": sibling_text,
    }


def validate_health_receipt(value: Any, where: str) -> Dict[str, Any]:
    if type(value) is not dict:
        fail("{} is not an object".format(where))
    exact_keys(value, {
        "available_cpu_count", "controller_affinity", "controller_cpu",
        "edac", "governor", "loadavg",
        "sibling_cpu", "sibling_non_idle_ticks", "sibling_proc_stat",
        "sibling_tick_fields", "target_cpu", "thermal",
        "thermal_endpoint_max_millic", "thread_siblings_list",
    }, where)
    exact_int(value["available_cpu_count"], 1, 1,
              where + ".available_cpu_count")
    exact_int(value["controller_cpu"], CONTROLLER_CPU, CONTROLLER_CPU,
              where + ".controller_cpu")
    if value["controller_affinity"] != [CONTROLLER_CPU]:
        fail("{}.controller_affinity differs".format(where))
    exact_int(value["target_cpu"], TARGET_CPU, TARGET_CPU,
              where + ".target_cpu")
    exact_int(value["sibling_cpu"], TARGET_SIBLING, TARGET_SIBLING,
              where + ".sibling_cpu")
    exact_int(
        value["thermal_endpoint_max_millic"],
        THERMAL_ENDPOINT_MAX_MILLIC, THERMAL_ENDPOINT_MAX_MILLIC,
        where + ".thermal_endpoint_max_millic")
    governor = bounded_string(value["governor"], where + ".governor", 256)
    if not governor:
        fail("{}.governor is empty".format(where))
    bounded_string(value["loadavg"], where + ".loadavg", 1024)
    sibling_text = bounded_string(
        value["thread_siblings_list"], where + ".thread_siblings_list", 256)
    if sorted(parse_cpu_list(sibling_text, where + ".thread_siblings_list")) != \
            sorted((TARGET_CPU, TARGET_SIBLING)):
        fail("{} sibling topology differs".format(where))
    tick_fields = value["sibling_tick_fields"]
    if type(tick_fields) is not list or len(tick_fields) != 8:
        fail("{}.sibling_tick_fields differs".format(where))
    parsed_ticks = [
        exact_int(item, 0, MAX_INT63,
                  "{}.sibling_tick_fields[{}]".format(where, index))
        for index, item in enumerate(tick_fields)
    ]
    non_idle = sum(parsed_ticks[index] for index in (0, 1, 2, 5, 6, 7))
    exact_int(value["sibling_non_idle_ticks"], non_idle, non_idle,
              where + ".sibling_non_idle_ticks")
    proc_stat = bounded_string(
        value["sibling_proc_stat"], where + ".sibling_proc_stat", 4096)
    proc_fields = proc_stat.split()
    if len(proc_fields) < 9 or proc_fields[0] != "cpu{}".format(TARGET_SIBLING):
        fail("{}.sibling_proc_stat differs".format(where))
    try:
        proc_ticks = [int(item) for item in proc_fields[1:9]]
    except ValueError:
        fail("{}.sibling_proc_stat is not numeric".format(where))
    if proc_ticks != parsed_ticks:
        fail("{}.sibling_proc_stat does not match tick fields".format(where))

    thermal = value["thermal"]
    if type(thermal) is not list or not thermal or len(thermal) > 4096:
        fail("{}.thermal roster differs".format(where))
    thermal_paths: List[str] = []
    for index, item in enumerate(thermal):
        item_where = "{}.thermal[{}]".format(where, index)
        if type(item) is not dict:
            fail("{} is not an object".format(item_where))
        exact_keys(item, {"path", "value_millic"}, item_where)
        path = bounded_string(item["path"], item_where + ".path", 4096)
        if not Path(path).is_absolute() or not path.startswith("/sys/class/"):
            fail("{} is not an absolute sysfs endpoint".format(item_where))
        exact_int(item["value_millic"], 0, THERMAL_ENDPOINT_MAX_MILLIC,
                  item_where + ".value_millic")
        thermal_paths.append(path)
    if thermal_paths != sorted(thermal_paths) or \
            len(thermal_paths) != len(set(thermal_paths)):
        fail("{}.thermal roster is not unique and sorted".format(where))

    edac = value["edac"]
    if type(edac) is not list or len(edac) > 4096:
        fail("{}.edac roster differs".format(where))
    edac_keys: List[Tuple[str, str]] = []
    for index, item in enumerate(edac):
        item_where = "{}.edac[{}]".format(where, index)
        if type(item) is not dict:
            fail("{} is not an object".format(item_where))
        exact_keys(item, {"counter", "path", "value"}, item_where)
        counter = bounded_string(item["counter"], item_where + ".counter", 32)
        if counter not in ("ce_count", "ue_count"):
            fail("{} counter differs".format(item_where))
        path = bounded_string(item["path"], item_where + ".path", 4096)
        if not Path(path).is_absolute():
            fail("{} path is not absolute".format(item_where))
        exact_int(item["value"], 0, MAX_INT63, item_where + ".value")
        edac_keys.append((counter, path))
    if edac_keys != sorted(edac_keys) or len(edac_keys) != len(set(edac_keys)):
        fail("{}.edac roster is not unique and sorted".format(where))
    return value


def validate_health_pair(
    before: Mapping[str, Any], after: Mapping[str, Any],
) -> Dict[str, Any]:
    validate_health_receipt(before, "health_before")
    validate_health_receipt(after, "health_after")
    for field in (
        "available_cpu_count", "controller_affinity", "controller_cpu",
        "governor", "sibling_cpu", "target_cpu", "thread_siblings_list",
    ):
        if before.get(field) != after.get(field):
            fail("health field {} changed during the worker".format(field))
    start_ticks = before.get("sibling_non_idle_ticks")
    end_ticks = after.get("sibling_non_idle_ticks")
    if type(start_ticks) is not int or type(end_ticks) is not int:
        fail("sibling non-idle tick receipt is malformed")
    tick_delta = end_ticks - start_ticks
    if tick_delta < 0 or tick_delta > SIBLING_NON_IDLE_TICK_CAP:
        fail("SMT sibling exceeded the frozen non-idle tick cap")
    before_edac = {
        (item["counter"], item["path"]): item["value"]
        for item in before.get("edac", [])
    }
    after_edac = {
        (item["counter"], item["path"]): item["value"]
        for item in after.get("edac", [])
    }
    if set(before_edac) != set(after_edac):
        fail("EDAC counter roster changed during the worker")
    if any(after_edac[key] != before_edac[key] for key in before_edac):
        fail("EDAC counter increased during the worker")
    before_thermal = {item["path"] for item in before.get("thermal", [])}
    after_thermal = {item["path"] for item in after.get("thermal", [])}
    if before_thermal != after_thermal:
        fail("thermal sensor roster changed during the worker")
    return {
        "edac_no_increase_verified": True,
        "governor_stable_verified": True,
        "sibling_non_idle_tick_cap": SIBLING_NON_IDLE_TICK_CAP,
        "sibling_non_idle_tick_delta": tick_delta,
        "thermal_endpoint_max_millic": THERMAL_ENDPOINT_MAX_MILLIC,
        "thermal_endpoint_samples_within_cap_verified": True,
        "thermal_roster_stable_verified": True,
        "topology_stable_verified": True,
    }


def expected_from_args(
    args: argparse.Namespace, compiler: Mapping[str, Any], library: Path,
) -> Dict[str, str]:
    expected = {
        "compiler_path": str(compiler["path"]),
        "compiler_sha256": str(compiler["sha256"]),
        "compiler_version": str(compiler["version"]),
        "library_path": str(library),
        "source_commit": args.expected_source_commit,
    }
    lower_commit(
        args.expected_source_commit, "argument --expected-source-commit")
    for field in (
        "expected_controller_sha256", "expected_worker_sha256",
        "expected_library_sha256",
    ):
        lower_hash(getattr(args, field), "argument --" + field.replace("_", "-"))
    return expected


def native_config_receipt(
    worker: Path, library: Path, expected: Mapping[str, str],
) -> Dict[str, Any]:
    namespace_before = fixed_namespace_snapshot()
    env = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
        "LD_LIBRARY_PATH": str(library.parent),
    }
    try:
        completed = subprocess.run(
            [str(worker), "--emit-config-selftest"], stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300.0,
            check=False, env=env,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("native config selftest failed: {}".format(exc))
    if (
        completed.returncode != 0 or not completed.stdout or completed.stderr or
        len(completed.stdout) > MAX_LINE_BYTES or
        len(completed.stderr) > MAX_STDERR_BYTES
    ):
        fail("native config selftest returned an invalid process result")
    config = parse_canonical_line(completed.stdout, "native config selftest")
    validate_config(config, expected)
    if fixed_namespace_snapshot() != namespace_before:
        fail("native config selftest touched the fixed one-attempt namespace")
    return {
        "config_sha256": sha256_bytes(completed.stdout),
        "native_config_identity_sha256":
            config["native_config_identity_sha256"],
        "runtime_library_maps_sha256": config["runtime_library_maps_sha256"],
        "runtime_library_path": config["runtime_library_path"],
        "stderr_sha256": sha256_bytes(completed.stderr),
        "target_identity_sha256":
            config["target_identity"]["canonical_sha256"],
    }


def preflight(
    args: argparse.Namespace, controller_affinity: Mapping[str, Any],
) -> Dict[str, Any]:
    source_root = exact_absolute_directory(args.source_root, "--source-root")
    build_root = exact_absolute_directory(args.build_root, "--build-root")
    worker = exact_absolute_file(args.worker, "--worker")
    library = exact_absolute_file(args.library, "--library")
    compiler_path = exact_absolute_file(args.compiler, "--compiler")
    controller = Path(__file__).resolve()
    if file_sha256(controller) != args.expected_controller_sha256:
        fail("controller SHA-256 differs from its exact argument")
    if file_sha256(worker) != args.expected_worker_sha256:
        fail("worker SHA-256 differs from its exact argument")
    if file_sha256(library) != args.expected_library_sha256:
        fail("library SHA-256 differs from its exact argument")
    compiler = compiler_receipt(compiler_path)
    interpreter = interpreter_receipt()
    if (
        interpreter["implementation"] != "cpython" or
        not interpreter["isolated"] or
        not interpreter["dont_write_bytecode"] or
        not interpreter["ignore_environment"] or
        not interpreter["no_site"] or
        not interpreter["no_user_site"]
    ):
        fail("real controller must run under CPython with exact -I -B -S isolation")
    expected = expected_from_args(args, compiler, library)
    git = git_receipt(source_root, args.expected_source_commit)
    build = build_receipt(
        source_root, build_root, compiler, interpreter, worker, library,
        args.expected_source_commit)
    dynamic = dynamic_receipt(worker, library)
    native_config = native_config_receipt(worker, library, expected)
    expected["native_config_identity_sha256"] = (
        native_config["native_config_identity_sha256"])
    expected["target_identity_sha256"] = (
        native_config["target_identity_sha256"])
    r0 = snapshot_r0()
    health = cpu_receipt()
    artifacts = {
        "compiler": artifact_receipt(compiler_path),
        "controller": artifact_receipt(controller),
        "controller_interpreter": interpreter["artifact"],
        "library": artifact_receipt(library),
        "worker": artifact_receipt(worker),
    }
    for name, wanted in (
        ("controller", args.expected_controller_sha256),
        ("library", args.expected_library_sha256),
        ("worker", args.expected_worker_sha256),
        ("compiler", compiler["sha256"]),
    ):
        exact_string(
            artifacts[name]["sha256"], wanted,
            "immediate pre-CLAIM {} SHA-256".format(name))
    return {
        "artifacts": artifacts,
        "build": build,
        "build_root": str(build_root),
        "compiler": compiler,
        "controller_affinity": dict(controller_affinity),
        "controller_path": str(controller),
        "dynamic": dynamic,
        "dynamic_after": None,
        "expected": expected,
        "git": git,
        "git_after": None,
        "health_adjudication": None,
        "health_after": None,
        "health_before": health,
        "interpreter": interpreter,
        "interpreter_after": None,
        "library_path": str(library),
        "native_config": native_config,
        "r0_after": None,
        "r0_before": r0,
        "source_root": str(source_root),
        "worker_path": str(worker),
        "artifacts_after": None,
        "build_after": None,
    }


def write_exclusive(directory_fd: int, name: str, data: bytes, mode: int) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    flags |= getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(name, flags, mode, dir_fd=directory_fd)
    try:
        os.fchmod(fd, mode)
        info = os.fstat(fd)
        if (
            not stat.S_ISREG(info.st_mode) or
            stat.S_IMODE(info.st_mode) != mode or info.st_nlink != 1
        ):
            fail("new evidence file metadata differs for {}".format(name))
        offset = 0
        while offset < len(data):
            written = os.write(fd, data[offset:])
            if written <= 0:
                fail("short write for {}".format(name))
            offset += written
        os.fsync(fd)
    finally:
        os.close(fd)


def publish_json(directory_fd: int, name: str, value: Any) -> None:
    write_exclusive(directory_fd, name, canonical_bytes(value) + b"\n", 0o400)


def claim_namespace(claim: Mapping[str, Any]) -> int:
    parent_fd = os.open(
        str(FIXED_OUTPUT_DIR.parent),
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_CLOEXEC", 0))
    try:
        try:
            os.mkdir(str(FIXED_OUTPUT_DIR), 0o700)
        except FileExistsError:
            fail("fixed r1 namespace already exists; a real attempt is already spent")
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    directory_fd = -1
    try:
        flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        directory_fd = os.open(str(FIXED_OUTPUT_DIR), flags)
        info = os.fstat(directory_fd)
        if stat.S_IMODE(info.st_mode) != 0o700:
            fail("fixed namespace mode differs")
        publish_json(directory_fd, "CLAIM", claim)
        os.fsync(directory_fd)
        return directory_fd
    except BaseException:
        if directory_fd >= 0:
            os.close(directory_fd)
        raise


def read_worker_started(
    directory_fd: int, claim: Mapping[str, Any], claim_sha256: str,
) -> Optional[Tuple[Dict[str, Any], str]]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(WORKER_STARTED_FILE, flags, dir_fd=directory_fd)
    except FileNotFoundError:
        return None
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or
            stat.S_IMODE(before.st_mode) != 0o400 or before.st_nlink != 1 or
            before.st_size < 0 or before.st_size > MAX_PROBE_BYTES
        ):
            fail("WORKER_STARTED metadata differs")
        chunks: List[bytes] = []
        total = 0
        while True:
            chunk = os.read(fd, min(65536, MAX_PROBE_BYTES + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > MAX_PROBE_BYTES:
                fail("WORKER_STARTED exceeds its byte cap")
        after = os.fstat(fd)
        stable_fields = (
            "st_dev", "st_ino", "st_mode", "st_nlink", "st_size",
            "st_mtime_ns",
        )
        if any(getattr(before, field) != getattr(after, field)
               for field in stable_fields):
            fail("WORKER_STARTED changed while it was read")
    finally:
        os.close(fd)
    data = b"".join(chunks)
    marker = parse_canonical_line(data, "WORKER_STARTED")
    exact_keys(marker, WORKER_STARTED_KEYS, "WORKER_STARTED")
    for value, expected, where in (
        (marker["campaign"], CAMPAIGN, "WORKER_STARTED.campaign"),
        (marker["claim_sha256"], claim_sha256,
         "WORKER_STARTED.claim_sha256"),
        (marker["native_config_identity_sha256"],
         claim["native_config_identity_sha256"],
         "WORKER_STARTED.native_config_identity_sha256"),
        (marker["schema"], WORKER_STARTED_SCHEMA, "WORKER_STARTED.schema"),
        (marker["source_commit"], claim["source_commit"],
         "WORKER_STARTED.source_commit"),
        (marker["status"], "started", "WORKER_STARTED.status"),
        (marker["worker_sha256"], claim["worker_sha256"],
         "WORKER_STARTED.worker_sha256"),
    ):
        exact_string(value, expected, where)
    return marker, sha256_bytes(data)


def run_worker(
    worker: str, library: str, deadline: float, claim_sha256: str,
    worker_sha256: str, native_config_identity_sha256: str,
) -> Tuple[bytes, bytes, int, float, bool]:
    env = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
        "LD_LIBRARY_PATH": str(Path(library).parent),
    }
    started = time.monotonic()
    argv = [
        worker, "--run", "--cpu", str(TARGET_CPU),
        "--evidence-dir", str(FIXED_OUTPUT_DIR),
        "--claim-sha256", claim_sha256,
        "--worker-sha256", worker_sha256,
        "--config-identity-sha256", native_config_identity_sha256,
    ]
    process: Optional[subprocess.Popen] = None
    selector: Optional[selectors.BaseSelector] = None
    streams: Dict[int, Tuple[str, Any]] = {}
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    caps = {"stdout": MAX_RAW_BYTES, "stderr": MAX_STDERR_BYTES}
    timed_out = False
    shutdown_stage: Optional[str] = None
    shutdown_deadline = 0.0
    overflow: Optional[str] = None
    leaked_pipe = False
    exit_seen_at: Optional[float] = None

    def signal_group(sig: int) -> None:
        if process is None:
            return
        try:
            os.killpg(process.pid, sig)
        except ProcessLookupError:
            pass
        except OSError as exc:
            fail("cannot signal worker process group: {}".format(exc))

    try:
        process = subprocess.Popen(
            argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env=env, start_new_session=True,
            bufsize=0,
        )
        if process.stdout is None or process.stderr is None:
            fail("worker pipes are unavailable")
        selector = selectors.DefaultSelector()
        for name, stream in (
            ("stdout", process.stdout), ("stderr", process.stderr)
        ):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, name)
            streams[stream.fileno()] = (name, stream)

        while streams or process.poll() is None:
            now = time.monotonic()
            if shutdown_stage is None and (now >= deadline or overflow is not None):
                if now >= deadline:
                    timed_out = True
                signal_group(signal.SIGTERM)
                shutdown_stage = "term"
                shutdown_deadline = now + 2.0
            elif shutdown_stage == "term" and now >= shutdown_deadline:
                signal_group(signal.SIGKILL)
                shutdown_stage = "kill"
                shutdown_deadline = now + 5.0
            elif shutdown_stage == "kill" and now >= shutdown_deadline:
                fail("worker process group did not quiesce after SIGKILL")

            returncode = process.poll()
            if returncode is not None and exit_seen_at is None:
                exit_seen_at = now
            if (
                returncode is not None and streams and
                exit_seen_at is not None and now - exit_seen_at >= 2.0 and
                shutdown_stage is None
            ):
                leaked_pipe = True
                signal_group(signal.SIGKILL)
                shutdown_stage = "kill"
                shutdown_deadline = now + 5.0

            next_deadline = deadline if shutdown_stage is None else shutdown_deadline
            timeout = max(0.0, min(0.25, next_deadline - now))
            events = selector.select(timeout)
            for key, _ in events:
                stream = key.fileobj
                name = key.data
                try:
                    chunk = os.read(stream.fileno(), 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    fail("cannot read worker {}: {}".format(name, exc))
                if not chunk:
                    selector.unregister(stream)
                    streams.pop(stream.fileno(), None)
                    stream.close()
                    continue
                room = caps[name] + 1 - len(buffers[name])
                if room > 0:
                    buffers[name].extend(chunk[:room])
                if len(chunk) > max(0, room) or len(buffers[name]) > caps[name]:
                    overflow = name

        try:
            returncode = process.wait(timeout=0.0)
        except subprocess.TimeoutExpired:
            fail("worker was not reaped after pipe closure")
        elapsed = time.monotonic() - started
        if overflow is not None:
            fail("worker {} exceeds its byte cap".format(overflow))
        if leaked_pipe:
            fail("worker descendant retained an output pipe")
        return (
            bytes(buffers["stdout"]), bytes(buffers["stderr"]),
            returncode, elapsed, timed_out,
        )
    except BaseException as exc:
        if process is not None:
            try:
                signal_group(signal.SIGTERM)
                process.wait(timeout=0.5)
            except (OSError, subprocess.TimeoutExpired, ValidationError):
                pass
            finally:
                # The leader may already have exited while descendants still
                # hold pipes.  Always finish cleanup of the owned process group.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except OSError:
                    pass
                try:
                    process.wait(timeout=5.0)
                except (OSError, subprocess.TimeoutExpired):
                    pass
        if isinstance(exc, ValidationError):
            raise
        raise ValidationError("worker supervision failed: {}".format(exc)) from None
    finally:
        if selector is not None:
            selector.close()
        # Include pipes whose selector registration failed part way through
        # setup, so exception cleanup does not leak their file descriptors.
        for stream in (() if process is None else
                       (process.stdout, process.stderr)):
            if stream is None:
                continue
            try:
                stream.close()
            except OSError:
                pass


def make_summary(
    outcome: str, statistics: Optional[Mapping[str, Any]],
    infrastructure: Sequence[str], transcript_invalid: Sequence[str],
    reject: Sequence[str], raw_sha256: str, elapsed: float,
) -> Dict[str, Any]:
    return {
        "campaign": CAMPAIGN,
        "elapsed_seconds": elapsed,
        "infrastructure_failures": list(infrastructure),
        "outcome": outcome,
        "raw_sha256": raw_sha256,
        "reject_gates": list(reject),
        "schema": SUMMARY_SCHEMA,
        "statistics": statistics,
        "transcript_invalid_gates": list(transcript_invalid),
    }


def postflight(
    args: argparse.Namespace, provenance: Dict[str, Any],
) -> None:
    provenance["health_after"] = cpu_receipt()
    provenance["health_adjudication"] = validate_health_pair(
        provenance["health_before"], provenance["health_after"])
    provenance["r0_after"] = snapshot_r0()
    if canonical_bytes(provenance["r0_before"]) != canonical_bytes(
        provenance["r0_after"]
    ):
        fail("r0 artifact metadata changed")

    compiler_path = Path(provenance["compiler"]["path"])
    controller = Path(provenance["controller_path"])
    library = Path(provenance["library_path"])
    worker = Path(provenance["worker_path"])
    source_root = Path(provenance["source_root"])
    build_root = Path(provenance["build_root"])
    provenance["artifacts_after"] = {
        "compiler": artifact_receipt(compiler_path),
        "controller": artifact_receipt(controller),
        "controller_interpreter": artifact_receipt(
            Path(provenance["interpreter"]["path"])),
        "library": artifact_receipt(library),
        "worker": artifact_receipt(worker),
    }
    if canonical_bytes(provenance["artifacts_after"]) != canonical_bytes(
        provenance["artifacts"]
    ):
        fail("artifact identity changed during the worker")
    provenance["interpreter_after"] = interpreter_receipt()
    if canonical_bytes(provenance["interpreter_after"]) != canonical_bytes(
        provenance["interpreter"]
    ):
        fail("controller interpreter identity changed during the worker")
    provenance["git_after"] = git_receipt(
        source_root, args.expected_source_commit)
    if canonical_bytes(provenance["git_after"]) != canonical_bytes(
        provenance["git"]
    ):
        fail("Git receipt changed during the worker")
    provenance["build_after"] = build_receipt(
        source_root, build_root, provenance["compiler"],
        provenance["interpreter"], worker, library,
        args.expected_source_commit)
    if canonical_bytes(provenance["build_after"]) != canonical_bytes(
        provenance["build"]
    ):
        fail("build receipt changed during the worker")
    provenance["dynamic_after"] = dynamic_receipt(worker, library)
    if canonical_bytes(provenance["dynamic_after"]) != canonical_bytes(
        provenance["dynamic"]
    ):
        fail("dynamic-link receipt changed during the worker")


def run_once(args: argparse.Namespace) -> int:
    # All fallible read-only provenance checks precede consumption of the fixed
    # one-attempt namespace.  The worker cannot run until CLAIM is durable.
    controller_affinity = pin_controller()
    provenance = preflight(args, controller_affinity)
    claim = {
        "campaign": CAMPAIGN,
        "controller_interpreter_sha256":
            provenance["interpreter"]["artifact"]["sha256"],
        "controller_sha256": args.expected_controller_sha256,
        "created_unix_ns": time.time_ns(),
        "library_sha256": args.expected_library_sha256,
        "native_config_identity_sha256":
            provenance["native_config"]["native_config_identity_sha256"],
        "schema": CLAIM_SCHEMA,
        "source_commit": args.expected_source_commit,
        "worker_sha256": args.expected_worker_sha256,
    }
    directory_fd = claim_namespace(claim)
    claim_sha256 = sha256_bytes(canonical_bytes(claim) + b"\n")
    deadline = time.monotonic() + OUTER_DEADLINE_SECONDS
    raw = b""
    stderr = b""
    elapsed = 0.0
    infrastructure: List[str] = []
    transcript_invalid: List[str] = []
    reject: List[str] = []
    statistics: Optional[Dict[str, Any]] = None
    worker_exit: Optional[int] = None
    timed_out = False
    worker_started_sha256: Optional[str] = None
    try:
        try:
            raw, stderr, worker_exit, elapsed, timed_out = run_worker(
                provenance["worker_path"], provenance["library_path"], deadline,
                claim_sha256, args.expected_worker_sha256,
                provenance["native_config"][
                    "native_config_identity_sha256"])
            if timed_out:
                infrastructure.append("worker_outer_deadline")
            if worker_exit != 0:
                infrastructure.append("worker_exit:{}".format(worker_exit))
            if stderr:
                infrastructure.append("worker_stderr_nonempty")
        except Exception as exc:
            infrastructure.append("execution:" + str(exc))

        try:
            marker_receipt = read_worker_started(
                directory_fd, claim, claim_sha256)
            if marker_receipt is None:
                infrastructure.append("worker_started_missing")
            else:
                _, worker_started_sha256 = marker_receipt
        except Exception as exc:
            infrastructure.append("worker_started:" + str(exc))

        try:
            postflight(args, provenance)
        except Exception as exc:
            infrastructure.append("postflight:" + str(exc))

        if worker_exit == 0 and not timed_out:
            try:
                config_record = parse_canonical_line(
                    raw.splitlines(keepends=True)[0], "raw config")
                runtime_maps_sha256 = lower_hash(
                    config_record["runtime_library_maps_sha256"],
                    "config.runtime_library_maps_sha256")
                _, _, statistics, transcript_invalid, reject = parse_transcript(
                    raw, provenance["expected"], runtime_maps_sha256)
            except (IndexError, KeyError, ValidationError) as exc:
                infrastructure.append("transcript:" + str(exc))
        if time.monotonic() > deadline:
            infrastructure.append("outer_deadline")
        outcome = (
            "invalid" if infrastructure or transcript_invalid else
            "reject" if reject else "pass")
        raw_sha256 = sha256_bytes(raw)
        summary = make_summary(
            outcome, statistics, infrastructure, transcript_invalid, reject,
            raw_sha256, elapsed)
        provenance.update({
            "campaign": CAMPAIGN,
            "claim_sha256": claim_sha256,
            "outer_deadline_seconds": int(OUTER_DEADLINE_SECONDS),
            "raw_sha256": raw_sha256,
            "schema": PROVENANCE_SCHEMA,
            "stderr_sha256": sha256_bytes(stderr),
            "worker_exit": worker_exit,
            "worker_started_sha256": worker_started_sha256,
            "worker_timed_out": timed_out,
        })
        summary_sha256 = sha256_bytes(canonical_bytes(summary) + b"\n")
        provenance_sha256 = sha256_bytes(canonical_bytes(provenance) + b"\n")
        complete = {
            "campaign": CAMPAIGN,
            "claim_sha256": claim_sha256,
            "outcome": outcome,
            "provenance_sha256": provenance_sha256,
            "raw_sha256": raw_sha256,
            "schema": COMPLETE_SCHEMA,
            "status": "complete",
            "summary_sha256": summary_sha256,
        }
        # CLAIM already exists.  COMPLETE is deliberately the final file and
        # the final durable publication operation.
        write_exclusive(directory_fd, "raw.jsonl", raw, 0o400)
        write_exclusive(directory_fd, "worker.stderr", stderr, 0o400)
        publish_json(directory_fd, "summary.json", summary)
        publish_json(directory_fd, "provenance.json", provenance)
        publish_json(directory_fd, "COMPLETE", complete)
        os.fsync(directory_fd)
        return 0 if outcome == "pass" else 2 if outcome == "reject" else 1
    finally:
        os.close(directory_fd)


def synthetic_hash(label: str) -> str:
    return sha256_bytes(label.encode("ascii"))


def synthetic_oracle(arm: str, k: int, source_hash: str, label: str) -> Dict[str, Any]:
    profile = synthetic_hash("profile-" + label) if arm == "C" else None
    return {
        "arm": arm,
        "attached_overlap_no_write_verified": True if arm == "C" else None,
        "borrowed_eligible_systematic_ids": k if arm == "C" else 0,
        "descriptor_sha256": synthetic_hash("descriptor-" + arm),
        "equation_configuration_sha256": synthetic_hash("equation-" + label),
        "first_repair_sha256": synthetic_hash("first-repair-" + label),
        "full_repair_sha256": synthetic_hash("full-repair-" + label),
        "high_id_sha256": synthetic_hash("high-id-" + label),
        "profile_bytes": 32 if arm == "C" else 0,
        "profile_hex": profile if arm == "C" else None,
        "profile_sha256": (
            sha256_bytes(bytes.fromhex(profile)) if profile is not None else None
        ),
        "roundtrip_first_id": k,
        "roundtrip_packet_count": k,
        "roundtrip_repair_only_verified": True,
        "roundtrip_sha256": source_hash,
        "roundtrip_verified": True,
        "systematic_sha256": source_hash,
    }


def synthetic_partial(
    k: int, block_bytes: int, tail: int, cell_label: str,
) -> Dict[str, Any]:
    label = "{}-tail{}".format(cell_label, tail)
    source_hash = SOURCE_SHA256[(k, block_bytes, tail)]
    systematic_hash = synthetic_hash("partial-systematic-" + label)
    arms = []
    for arm in ("C", "L"):
        arm_label = label + "-" + arm
        arms.append({
            "arm": arm,
            "attached_overlap_no_write_verified": True if arm == "C" else None,
            "first_repair_sha256": synthetic_hash("partial-first-" + arm_label),
            "high_id_sha256": synthetic_hash("partial-high-" + arm_label),
            "profile_sha256": (
                synthetic_hash("partial-profile-" + arm_label)
                if arm == "C" else None
            ),
            "roundtrip_first_id": k,
            "roundtrip_packet_count": k,
            "roundtrip_repair_only_verified": True,
            "roundtrip_sha256": source_hash,
            "source_immutable_verified": True,
            "systematic_sha256": systematic_hash,
            "systematic_verified": True,
        })
    return {"arms": arms, "source_sha256": source_hash, "tail_bytes": tail}


def synthetic_config(expected: Mapping[str, str], maps_hash: str) -> Dict[str, Any]:
    cells = []
    for k, block_bytes in CELLS:
        label = "K{}-B{}".format(k, block_bytes)
        source_hash = SOURCE_SHA256[(k, block_bytes, block_bytes)]
        cells.append({
            "K": k,
            "block_bytes": block_bytes,
            "final_bytes": block_bytes,
            "invocations_by_replicate": invocations_by_replicate(k),
            "message_bytes": message_bytes(k, block_bytes),
            "oracles": [
                synthetic_oracle("C", k, source_hash, label + "-C"),
                synthetic_oracle("L", k, source_hash, label + "-L"),
            ],
            "partial_final_oracles": [
                synthetic_partial(k, block_bytes, 1, label),
                synthetic_partial(k, block_bytes, block_bytes - 1, label),
            ],
            "source_seed": source_seed_text(k, block_bytes, block_bytes),
            "source_sha256": source_hash,
        })
    descriptors = [
        {
            "api": "wirehair_v2_encoder_create_with_options",
            "arm": "C", "codec": "wirehair2",
            "equation_profile": "certified-2026-07",
            "schema": SCHEMA_PREFIX + ".arm.v1",
            "source_policy": "borrowed-immutable",
        },
        {
            "api": "wirehair_encoder_create_ex", "arm": "L",
            "codec": "wirehair1", "equation_profile": "legacy-current",
            "schema": SCHEMA_PREFIX + ".arm.v1",
            "source_policy": "borrowed",
        },
    ]
    identity_bytes = b"synthetic-target-identity"
    target_identity = {
        "after_cpu": TARGET_CPU,
        "before_cpu": TARGET_CPU,
        "canonical_bytes": len(identity_bytes),
        "canonical_hex": identity_bytes.hex(),
        "canonical_sha256": sha256_bytes(identity_bytes),
        "capabilities": {key: 0 for key in TARGET_CAPABILITY_KEYS},
        "derived": {key: 0 for key in TARGET_DERIVED_KEYS},
        "raw_capture_complete": True,
        "requested_cpu": TARGET_CPU,
        "semantic_validation_passed": True,
        "singleton_affinity_verified": True,
    }
    target_identity["derived"]["threads_per_core"] = 1
    target_identity["derived"]["logical_processors_per_package"] = 1
    config = {
        "arm_descriptors": [
            {"arm": value["arm"], "descriptor": value,
             "descriptor_sha256": sha256_bytes(canonical_bytes(value))}
            for value in descriptors
        ],
        "campaign": CAMPAIGN,
        "cells": cells,
        "comparisons": [
            {"left_arm": left, "name": name, "right_arm": right}
            for name, left, right in COMPARISONS
        ],
        "compile": {
            "compiler_path": expected["compiler_path"],
            "compiler_sha256": expected["compiler_sha256"],
            "compiler_version": expected["compiler_version"],
            "harness_git_commit": expected["source_commit"],
            "implementation_git_commit": expected["source_commit"],
        },
        "expected_encode_calls": EXPECTED_ENCODE_CALLS,
        "expected_measured_encode_calls": EXPECTED_MEASURED_ENCODE_CALLS,
        "expected_measured_invocations": EXPECTED_MEASURED_INVOCATIONS,
        "expected_panels": EXPECTED_PANEL_COUNT,
        "expected_records": EXPECTED_RECORD_COUNT,
        "expected_warmup_encode_calls": EXPECTED_WARMUP_ENCODE_CALLS,
        "expected_warmup_invocations": EXPECTED_WARMUP_INVOCATIONS,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "minimum_invocations": MINIMUM_INVOCATIONS,
        "native_config_identity_sha256": "0" * 64,
        "panel_key_schema": PANEL_KEY_SCHEMA,
        "panel_replicates": REPLICATES,
        "roster_order": ["replicate", "cell", "scope", "comparison"],
        "runtime_library_maps_sha256": maps_hash,
        "runtime_library_path": expected["library_path"],
        "schema": CONFIG_SCHEMA,
        "scopes": [
            {"first_id": SCOPE_SPECS[name][0], "name": name,
             "timed_construction": SCOPE_SPECS[name][1]}
            for name in SCOPES
        ],
        "sibling_cpu": TARGET_SIBLING,
        "source_seed_base": SOURCE_SEED_BASE_TEXT,
        "target_cpu": TARGET_CPU,
        "target_identity": target_identity,
    }
    config["native_config_identity_sha256"] = native_config_identity_sha256(
        config)
    return config


def synthetic_panel(
    config: Mapping[str, Any], replicate: int, k: int, block_bytes: int,
    scope: str, comparison: str, sequence: int, ratio: float,
) -> Dict[str, Any]:
    order = panel_order(replicate, k, block_bytes, scope, comparison)
    left, right = COMPARISON_ARMS[comparison]
    n = invocations_by_replicate(k)[replicate]
    slots = []
    left_count = 0
    right_count = 0
    for index, side in enumerate(sides_for(order)):
        count = (n + 1) // 2 if index < 4 else n // 2
        # Integer nanoseconds remain exact while preserving the intended ratio.
        unit = 980000 if side == "left" and ratio == 0.98 else 1000000
        if ratio not in (0.98, 1.0):
            unit = int(round(1000000 * ratio)) if side == "left" else 1000000
        slots.append({
            "elapsed_ns": count * unit,
            "invocation_count": count,
            "side": side,
        })
        if side == "left":
            left_count += count
        else:
            right_count += count
    oracle_map = {
        (cell["K"], cell["block_bytes"], oracle["arm"]): oracle
        for cell in config["cells"] for oracle in cell["oracles"]
    }
    return {
        "K": k,
        "affinity_verified": True,
        "block_bytes": block_bytes,
        "comparison": comparison,
        "cpu_migration_verified": True,
        "invocations_per_slot": n,
        "left_arm": left,
        "left_descriptor_sha256": oracle_map[
            (k, block_bytes, left)]["descriptor_sha256"],
        "left_outcome_code": 0,
        "left_output_sha256": oracle_map[(k, block_bytes, left)][
            "systematic_sha256" if scope in SYSTEMATIC_SCOPES
            else "full_repair_sha256"],
        "left_public_encode_calls": (left_count + 1) * k,
        "order": order,
        "panel_key_sha256": panel_key_sha256(
            k, block_bytes, scope, comparison),
        "replicate": replicate,
        "right_arm": right,
        "right_descriptor_sha256": oracle_map[
            (k, block_bytes, right)]["descriptor_sha256"],
        "right_outcome_code": 0,
        "right_output_sha256": oracle_map[(k, block_bytes, right)][
            "systematic_sha256" if scope in SYSTEMATIC_SCOPES
            else "full_repair_sha256"],
        "right_public_encode_calls": (right_count + 1) * k,
        "runtime_library_maps_sha256": config["runtime_library_maps_sha256"],
        "schema": PANEL_SCHEMA,
        "scope": scope,
        "sequence": sequence,
        "slots": slots,
        "target_cpu": TARGET_CPU,
    }


def synthetic_expected() -> Dict[str, str]:
    return {
        "compiler_path": "/usr/bin/g++-13",
        "compiler_sha256": "1" * 64,
        "compiler_version": "13.3.0",
        "library_path": "/tmp/r1/libwirehair.so.2.0.0",
        "source_commit": "2" * 40,
    }


def synthetic_transcript(
    cl_systematic_ratio: float = 0.98, cl_repair_ratio: float = 1.0,
    aa_ratio: float = 1.0,
) -> Tuple[bytes, Dict[str, str], str]:
    expected = synthetic_expected()
    maps_hash = synthetic_hash("runtime-maps")
    config = synthetic_config(expected, maps_hash)
    records: List[Dict[str, Any]] = [config]
    sequence = 0
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison, _, _ in COMPARISONS:
                    ratio = aa_ratio if comparison in AA_COMPARISONS else (
                        cl_systematic_ratio if scope in SYSTEMATIC_SCOPES
                        else cl_repair_ratio)
                    records.append(synthetic_panel(
                        config, replicate, k, block_bytes, scope, comparison,
                        sequence, ratio))
                    sequence += 1
    records.append({
        "encode_call_count": EXPECTED_ENCODE_CALLS,
        "measured_invocation_count": EXPECTED_MEASURED_INVOCATIONS,
        "panel_count": EXPECTED_PANEL_COUNT,
        "record_count": EXPECTED_RECORD_COUNT,
        "schema": TERMINAL_SCHEMA,
        "status": "complete",
        "warmup_invocation_count": EXPECTED_WARMUP_INVOCATIONS,
    })
    raw = b"".join(canonical_bytes(record) + b"\n" for record in records)
    return raw, expected, maps_hash


def expect_invalid(callback: Any, label: str) -> None:
    try:
        callback()
    except ValidationError:
        return
    fail("selftest mutation was accepted: {}".format(label))


def fixed_namespace_snapshot() -> Optional[Tuple[int, int, int, int]]:
    try:
        info = os.stat(str(FIXED_OUTPUT_DIR), follow_symlinks=False)
    except FileNotFoundError:
        return None
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)


def selftest() -> int:
    namespace_before = fixed_namespace_snapshot()
    if (
        len(CELLS) != 10 or len(SCOPES) != 4 or len(COMPARISONS) != 3 or
        EXPECTED_PANEL_COUNT != len(CELLS) * len(SCOPES) * len(
            COMPARISONS) * REPLICATES or EXPECTED_RECORD_COUNT != 1442
    ):
        fail("frozen roster constants differ")
    if invocations_by_replicate(8) != [683] * 8 + [682] * 4:
        fail("K=8 invocation allocation differs")
    if invocations_by_replicate(64000) != [2] * REPLICATES:
        fail("K=64000 invocation allocation differs")
    if sum(
        2 * len(SCOPES) * len(COMPARISONS) * 4 * total_invocations(k)
        for k in (8, 128, 512, 8192, 64000)
    ) != EXPECTED_MEASURED_INVOCATIONS:
        fail("measured invocation accounting differs")
    for order in ("ABBA", "BAAB"):
        values = [950000 if side == "left" else 1000000
                  for side in sides_for(order)]
        if abs(math.exp(lane_contrast(values, order)) - 0.95) > 1e-12:
            fail("lane contrast differs")

    raw, expected, maps_hash = synthetic_transcript()
    if len(raw.splitlines()) != EXPECTED_RECORD_COUNT:
        fail("synthetic transcript roster differs")
    config, terminal, statistics, invalid, reject = parse_transcript(
        raw, expected, maps_hash)
    if invalid or reject or terminal["record_count"] != EXPECTED_RECORD_COUNT:
        fail("synthetic passing transcript did not pass")
    replay = parse_transcript(raw, expected, maps_hash)
    if canonical_bytes(replay[2]) != canonical_bytes(statistics):
        fail("synthetic reduction is not deterministic")
    if config["cells"][0]["partial_final_oracles"][1]["tail_bytes"] != 63:
        fail("synthetic partial-final roster differs")

    expect_invalid(
        lambda: parse_canonical_line(b'{"a":1,"a":2}\n', "duplicate"),
        "duplicate JSON key")
    expect_invalid(
        lambda: parse_canonical_line(b'{"a": 1}\n', "spacing"),
        "noncanonical JSON spacing")
    truncated = b"\n".join(raw.splitlines()[:-1]) + b"\n"
    expect_invalid(
        lambda: parse_transcript(truncated, expected, maps_hash),
        "truncated roster")
    records = [json.loads(line) for line in raw.splitlines()]
    records[1]["sequence"] = 1
    mutated = b"".join(canonical_bytes(record) + b"\n" for record in records)
    expect_invalid(
        lambda: parse_transcript(mutated, expected, maps_hash),
        "sequence mutation")
    records = [json.loads(line) for line in raw.splitlines()]
    records[1]["slots"][0]["invocation_count"] += 1
    mutated = b"".join(canonical_bytes(record) + b"\n" for record in records)
    expect_invalid(
        lambda: parse_transcript(mutated, expected, maps_hash),
        "slot-count mutation")

    aa_raw, expected, maps_hash = synthetic_transcript(aa_ratio=1.03)
    _, _, _, invalid, _ = parse_transcript(aa_raw, expected, maps_hash)
    if not invalid or not all(item.startswith("aa_ci:") for item in invalid):
        fail("synthetic A/A invalidity was not classified")
    reject_raw, expected, maps_hash = synthetic_transcript(
        cl_systematic_ratio=1.0)
    _, _, _, invalid, reject = parse_transcript(
        reject_raw, expected, maps_hash)
    if invalid or not reject:
        fail("synthetic performance rejection was not classified")
    if fixed_namespace_snapshot() != namespace_before:
        fail("selftest mutated the fixed one-attempt namespace")
    print("WH2 public borrowed current-vs-WH1 r1 controller selftest passed")
    return 0


def native_config_selftest(args: argparse.Namespace) -> int:
    namespace_before = fixed_namespace_snapshot()
    worker = exact_absolute_file(args.worker, "--worker")
    library = exact_absolute_file(args.library, "--library")
    compiler_path = exact_absolute_file(args.compiler, "--compiler")
    source_commit = lower_commit(
        args.expected_source_commit, "--expected-source-commit")
    compiler = compiler_receipt(compiler_path)
    expected = {
        "compiler_path": str(compiler_path),
        "compiler_sha256": compiler["sha256"],
        "compiler_version": compiler["version"],
        "library_path": str(library),
        "source_commit": source_commit,
    }
    dynamic_receipt(worker, library)
    receipt = native_config_receipt(worker, library, expected)
    if fixed_namespace_snapshot() != namespace_before:
        fail("native config controller selftest touched the fixed namespace")
    print(canonical_bytes({
        "campaign": CAMPAIGN,
        "native_config": receipt,
        "status": "passed",
    }).decode("ascii"))
    return 0


def parse_canonical_document(path: Path, maximum: int) -> Dict[str, Any]:
    with path.open("rb") as source:
        data = source.read(maximum + 1)
    if len(data) > maximum:
        fail("{} exceeds its byte cap".format(path))
    return parse_canonical_line(data, str(path))


def replay_object(value: Any, keys: Iterable[str], where: str) -> Dict[str, Any]:
    if type(value) is not dict:
        fail("{} is not an object".format(where))
    exact_keys(value, keys, where)
    return value


def replay_equal(value: Any, expected: Any, where: str) -> None:
    if canonical_bytes(value) != canonical_bytes(expected):
        fail("{} differs".format(where))


def replay_path(value: Any, where: str) -> str:
    text = bounded_string(value, where, 4096)
    path = Path(text)
    # Replay validates recorded paths without requiring the original build or
    # machine to remain available.  It never executes any recorded command.
    if not path.is_absolute() or str(path) != text or ".." in path.parts:
        fail("{} is not a normalized absolute path".format(where))
    return text


def replay_text_hash(
    value: Mapping[str, Any], text_key: str, hash_key: str, where: str,
) -> str:
    text = bounded_string(value[text_key], where + "." + text_key,
                          MAX_PROBE_BYTES)
    exact_string(value[hash_key], sha256_bytes(text.encode("utf-8")),
                 where + "." + hash_key)
    return text


def validate_stat_receipt(
    value: Any, where: str, expected_path: Optional[str] = None,
    artifact: bool = False,
) -> Dict[str, Any]:
    replay_object(value, ARTIFACT_KEYS if artifact else STAT_RECEIPT_KEYS, where)
    replay_path(value["path"], where + ".path")
    if expected_path is not None:
        exact_string(value["path"], expected_path, where + ".path")
    for key in ("device", "inode", "mtime_ns", "size"):
        exact_int(value[key], 0, MAX_INT63, where + "." + key)
    exact_int(value["mode"], 0, 0o7777, where + ".mode")
    if artifact:
        exact_int(value["size"], 0, MAX_RAW_BYTES, where + ".size")
        lower_hash(value["sha256"], where + ".sha256")
    return value


def validate_interpreter_receipt(value: Any, where: str) -> None:
    replay_object(value, INTERPRETER_RECEIPT_KEYS, where)
    replay_path(value["path"], where + ".path")
    replay_path(value["invoked_path"], where + ".invoked_path")
    validate_stat_receipt(value["artifact"], where + ".artifact",
                          value["path"], artifact=True)
    exact_string(value["implementation"], "cpython", where + ".implementation")
    for key in ("dont_write_bytecode", "ignore_environment", "isolated",
                "no_site", "no_user_site"):
        exact_bool(value[key], True, where + "." + key)
    bounded_string(value["version"], where + ".version", 4096)
    version = value["version_info"]
    if type(version) is not list or len(version) != 5:
        fail("{}.version_info differs".format(where))
    exact_int(version[0], 3, 3, where + ".version_info.major")
    for index in (1, 2, 4):
        exact_int(version[index], 0, 65535, where + ".version_info")
    if version[3] not in ("alpha", "beta", "candidate", "final"):
        fail("{}.version_info releaselevel differs".format(where))
    if not value["version"].startswith("{}.{}.{}".format(*version[:3])):
        fail("{}.version does not match version_info".format(where))


def validate_git_receipt(value: Any, commit: str, where: str) -> None:
    replay_object(value, GIT_RECEIPT_KEYS, where)
    exact_string(value["commit"], commit, where + ".commit")
    exact_string(value["upstream_commit"], commit, where + ".upstream_commit")
    lower_commit(value["tree"], where + ".tree")
    exact_bool(value["tracked_clean"], True, where + ".tracked_clean")
    exact_string(value["git_invoked_path"], "/usr/bin/git",
                 where + ".git_invoked_path")
    validate_stat_receipt(value["git_artifact"], where + ".git_artifact",
                          artifact=True)
    version = replay_text_hash(value, "git_version_text", "git_version_sha256",
                               where)
    if not version.startswith("git version ") or not version.endswith("\n"):
        fail("{}.git_version_text differs".format(where))
    blobs = replay_object(value["required_tracked_inputs"], R1_TRACKED_INPUTS,
                           where + ".required_tracked_inputs")
    for name, digest in blobs.items():
        lower_commit(digest, where + ".required_tracked_inputs." + name)


def validate_build_receipt(
    value: Any, provenance: Mapping[str, Any], where: str,
) -> None:
    replay_object(value, BUILD_RECEIPT_KEYS, where)
    compiler = provenance["compiler"]
    source_root = Path(provenance["source_root"])
    build_root = Path(provenance["build_root"])
    worker = Path(provenance["worker_path"])
    library = Path(provenance["library_path"])
    replay_equal(value["cache"], {
        "BUILD_SHARED_LIBS": "ON", "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_GENERATOR": "Ninja", "MARCH_NATIVE": "ON", "WH_LTO": "OFF",
        "WH_PGO_MODE": "OFF",
    }, where + ".cache")
    for key, expected in (
        ("cmake_source_root", str(source_root)),
        ("cmake_build_root", str(build_root)),
        ("cmake_compiler_resolved", compiler["path"]),
        ("cmake_interpreter_resolved", provenance["interpreter"]["path"]),
    ):
        exact_string(value[key], expected, where + "." + key)
    configured_compiler = bounded_string(value["cmake_compiler"],
                                         where + ".cmake_compiler", 4096)
    if "/" in configured_compiler and not configured_compiler.startswith("/"):
        fail("{} has a relative compiler path".format(where))
    replay_path(value["cmake_interpreter"], where + ".cmake_interpreter")
    invocation = replay_path(value["compiler_invocation"],
                             where + ".compiler_invocation")
    dry_run = replay_text_hash(value, "dry_run_text", "dry_run_sha256", where)
    dry_lines = dry_run.splitlines()
    if not (len(dry_lines) in (1, 2) and
            dry_lines[-1] == "ninja: no work to do." and
            (len(dry_lines) == 1 or
             dry_lines[0].startswith("ninja: Entering directory `"))):
        fail("{}.dry_run_text does not show an up-to-date target".format(where))
    commands_text = replay_text_hash(value, "ninja_commands_text",
                                     "ninja_commands_sha256", where)
    try:
        command_rows = [shlex.split(line) for line in commands_text.splitlines()]
    except ValueError as exc:
        fail("{}.ninja_commands_text is malformed: {}".format(where, exc))
    definitions = [
        "-DWIREHAIR_DLL=1",
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_ID="GNU"',
        '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION="{}"'.format(
            compiler["version"]),
        '-DWIREHAIR_WH2_CXX_COMPILER_PATH="{}"'.format(compiler["path"]),
        '-DWIREHAIR_WH2_CXX_COMPILER_SHA256="{}"'.format(compiler["sha256"]),
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(
            provenance["expected"]["source_commit"]),
    ]
    library_definitions = [
        "-DWIREHAIR_BUILDING=1", "-DWIREHAIR_DLL=1", "-Dwirehair_EXPORTS",
    ]
    replay_equal(value["compile_definitions"], definitions,
                 where + ".compile_definitions")
    replay_equal(value["library_compile_definitions"], library_definitions,
                 where + ".library_compile_definitions")
    all_objects: List[List[str]] = []
    for prefix, sources, object_prefix, defs in (
        ("", R1_TARGET_SOURCES,
         "CMakeFiles/wirehair_wh2_public_borrowed_current_vs_wh1_r1.dir/",
         definitions),
        ("library_", R1_LIBRARY_SOURCES, "CMakeFiles/wirehair.dir/",
         library_definitions),
    ):
        source_key = "library_sources" if prefix else "target_sources"
        replay_equal(value[source_key], list(sources), where + "." + source_key)
        commands = replay_object(value[prefix + "compile_commands"], sources,
                                  where + "." + prefix + "compile_commands")
        hashes = replay_object(value[prefix + "compile_command_sha256"], sources,
                                where + "." + prefix + "compile_command_sha256")
        objects = [object_prefix + source + ".o" for source in sources]
        all_objects.append(objects)
        replay_equal(value[prefix + "link_object_roster"], objects,
                     where + "." + prefix + "link_object_roster")
        for source, obj in zip(sources, objects):
            line = bounded_string(commands[source], where + ".compile_command",
                                  MAX_PROBE_BYTES)
            exact_string(hashes[source], sha256_bytes(line.encode("utf-8")),
                         where + ".compile_command_sha256." + source)
            expected_tokens = [invocation] + defs + [
                "-I" + str(source_root / "include"), "-O3", "-DNDEBUG",
                "-std=gnu++11",
            ] + (["-fPIC"] if prefix else []) + [
                "-Wall", "-Wextra", "-march=native", "-MD", "-MT", obj,
                "-MF", obj + ".d", "-o", obj, "-c", str(source_root / source),
            ]
            try:
                tokens = shlex.split(line)
            except ValueError as exc:
                fail("{} compile command is malformed: {}".format(where, exc))
            replay_equal(tokens, expected_tokens, where + ".compile_command")
            if commands_text.splitlines().count(line) != 1:
                fail("{} compile command is not unique in Ninja receipt".format(where))
            matching = [row for row in command_rows if "-c" in row and
                        row[row.index("-c") + 1:] == [str(source_root / source)]]
            if matching != [tokens]:
                fail("{} Ninja compile source is not unique".format(where))
    target_objects, library_objects = all_objects
    expected_links = (
        ("", [":", "&&", invocation, "-O3", "-DNDEBUG"] + target_objects + [
            "-o", worker.name, "-Wl,-rpath,{}".format(build_root),
            library.name, "-ldl", "&&", ":",
        ]),
        ("library_", [":", "&&", invocation, "-fPIC", "-O3", "-DNDEBUG",
                      "-Wl,--version-script={}".format(
                          source_root / "abi" / "wirehair.map"),
                      "-shared", "-Wl,-soname,libwirehair.so.2", "-o",
                      library.name] + library_objects + ["-lm", "&&", ":"]),
    )
    for prefix, tokens in expected_links:
        line = replay_text_hash(value, prefix + "link_command",
                                 prefix + "link_command_sha256", where)
        try:
            actual = shlex.split(line)
        except ValueError as exc:
            fail("{} link command is malformed: {}".format(where, exc))
        replay_equal(actual, tokens, where + "." + prefix + "link_command")
        if commands_text.splitlines().count(line) != 1:
            fail("{} link command is not unique in Ninja receipt".format(where))


def validate_dynamic_receipt(
    value: Any, provenance: Mapping[str, Any], where: str,
) -> None:
    replay_object(value, DYNAMIC_RECEIPT_KEYS, where)
    worker = provenance["worker_path"]
    library = provenance["library_path"]
    for key, expected in (
        ("forbidden_symbol_fragments", list(FORBIDDEN_SYMBOL_FRAGMENTS)),
        ("worker_forbidden_symbol_fragments",
         list(WORKER_FORBIDDEN_SYMBOL_FRAGMENTS)),
        ("required_dynamic_symbols", sorted(REQUIRED_DYNAMIC_SYMBOLS)),
        ("worker_import_symbol_version", WIREHAIR_SYMBOL_VERSION),
        ("nm_invoked_path", "/usr/bin/nm"),
        ("library_nm_argv", ["/usr/bin/nm", "-D", "--format=posix",
                             "--defined-only", library]),
        ("worker_nm_argv", ["/usr/bin/nm", "-D", "--format=posix", worker]),
        ("worker_full_nm_argv", ["/usr/bin/nm", "-C", "--format=posix", worker]),
        ("ldd_library_resolved_path", library),
    ):
        replay_equal(value[key], expected, where + "." + key)
    for key in ("required_dynamic_symbols_verified", "worker_import_symbols_verified"):
        exact_bool(value[key], True, where + "." + key)
    validate_stat_receipt(value["nm_artifact"], where + ".nm_artifact",
                          "/usr/bin/x86_64-linux-gnu-nm", artifact=True)
    replay_text_hash(value, "nm_version_text", "nm_version_sha256", where)
    library_text = replay_text_hash(value, "library_symbol_table_text",
                                    "symbol_table_sha256", where)
    library_rows = parse_posix_nm(library_text.encode("utf-8"), where + ".library_nm")
    symbols = {name.split("@", 1)[0] for name, _ in library_rows}
    if REQUIRED_DYNAMIC_SYMBOLS - symbols or any(
        fragment in name for name in symbols for fragment in FORBIDDEN_SYMBOL_FRAGMENTS
    ):
        fail("{} library public symbol roster differs".format(where))
    worker_text = replay_text_hash(value, "worker_symbol_table_text",
                                   "worker_symbol_table_sha256", where)
    rows = parse_posix_nm(worker_text.encode("utf-8"), where + ".worker_nm")
    imports = [(name, kind) for name, kind in rows
               if name.split("@", 1)[0].startswith("wirehair")]
    expected_imports = sorted(name + "@" + WIREHAIR_SYMBOL_VERSION
                              for name in REQUIRED_DYNAMIC_SYMBOLS)
    replay_equal(sorted(imports), [(name, "U") for name in expected_imports],
                 where + ".worker imports")
    replay_equal(value["worker_import_symbols"], expected_imports,
                 where + ".worker_import_symbols")
    full_text = replay_text_hash(value, "worker_full_symbol_table_text",
                                 "worker_full_symbol_table_sha256", where)
    if any(fragment in full_text for fragment in WORKER_FORBIDDEN_SYMBOL_FRAGMENTS):
        fail("{} worker contains internal/test-hook symbols".format(where))
    readelf = replay_text_hash(value, "readelf_dynamic_text",
                               "readelf_dynamic_sha256", where)
    if "Shared library: [libwirehair" not in readelf:
        fail("{} lacks a libwirehair DT_NEEDED entry".format(where))
    ldd = replay_text_hash(value, "ldd_normalized_text",
                           "ldd_normalized_sha256", where)
    matches = [line for line in ldd.splitlines() if "libwirehair" in line]
    if len(matches) != 1:
        fail("{} ldd binding count differs".format(where))
    match = re.fullmatch(
        r"\s*libwirehair\.so[^\s]*\s+=>\s+(\S+)\s+\(0xADDR\)\s*", matches[0])
    if match is None:
        fail("{} ldd binding is malformed".format(where))
    replay_path(value["ldd_library_path"], where + ".ldd_library_path")
    exact_string(value["ldd_library_path"], match.group(1), where + ".ldd binding")
    # The loader spelling is commonly a SONAME symlink.  Its recorded lstat
    # identity is distinct from the resolved library artifact's identity.
    validate_stat_receipt(value["ldd_library_stat"], where + ".ldd_library_stat",
                          value["ldd_library_path"])


def validate_r0_receipt(value: Any, where: str) -> None:
    replay_object(value, R0_KEYS, where)
    entries = value["entries"]
    if type(entries) is not list or not entries or len(entries) > 65536:
        fail("{}.entries differs".format(where))
    exact_string(value["snapshot_sha256"], sha256_bytes(canonical_bytes(entries)),
                 where + ".snapshot_sha256")
    seen: set = set()
    roots = {str(path) for path in R0_ROOTS}
    for index, entry in enumerate(entries):
        item_where = "{}.entries[{}]".format(where, index)
        if type(entry) is not dict:
            fail("{} is not an object".format(item_where))
        path = replay_path(entry.get("path"), item_where + ".path")
        if path in seen or not any(path == root or path.startswith(root + "/")
                                   for root in roots):
            fail("{} path is duplicate or outside the r0 roots".format(item_where))
        seen.add(path)
        if entry.get("accessible") is False:
            replay_object(entry, {"accessible", "errno", "exists", "path"}, item_where)
            exact_int(entry["errno"], 1, 4095, item_where + ".errno")
            if entry["exists"] is not None and entry["exists"] is not True:
                fail("{} inaccessible existence receipt differs".format(item_where))
        elif entry.get("exists") is False:
            replay_object(entry, {"exists", "path"}, item_where)
            if path not in roots:
                fail("{} missing path is not an r0 root".format(item_where))
        else:
            replay_object(entry, STAT_RECEIPT_KEYS | {"exists", "type"} |
                          ({"sha256"} if "sha256" in entry else set()), item_where)
            exact_bool(entry["exists"], True, item_where + ".exists")
            if entry["type"] not in ("file", "directory", "symlink", "other"):
                fail("{} type differs".format(item_where))
            if "sha256" in entry:
                if entry["type"] != "file":
                    fail("{} non-file contains a file hash".format(item_where))
                lower_hash(entry["sha256"], item_where + ".sha256")
            validate_stat_receipt({key: entry[key] for key in STAT_RECEIPT_KEYS},
                                  item_where)
    if not roots.issubset(seen):
        fail("{} omits an r0 root".format(where))


def validate_provenance_receipt(
    value: Any, claim: Mapping[str, Any], infrastructure: Sequence[str],
) -> Dict[str, Any]:
    replay_object(value, PROVENANCE_KEYS, "provenance")
    exact_string(value["campaign"], CAMPAIGN, "provenance.campaign")
    exact_string(value["schema"], PROVENANCE_SCHEMA, "provenance.schema")
    exact_int(value["outer_deadline_seconds"], int(OUTER_DEADLINE_SECONDS),
              int(OUTER_DEADLINE_SECONDS), "provenance.outer_deadline_seconds")
    for key in ("claim_sha256", "raw_sha256", "stderr_sha256"):
        lower_hash(value[key], "provenance." + key)
    for key in ("source_root", "build_root", "controller_path", "worker_path",
                "library_path"):
        replay_path(value[key], "provenance." + key)
    source_root, build_root = Path(value["source_root"]), Path(value["build_root"])
    exact_string(value["controller_path"],
                 str(source_root / "bench/Wh2PublicBorrowedCurrentVsWh1R1.py"),
                 "provenance.controller_path")
    if not all(path_is_within(Path(value[key]), build_root)
               for key in ("worker_path", "library_path")) or \
            Path(value["library_path"]).name != "libwirehair.so.2.0.0":
        fail("provenance build artifact paths differ")
    compiler = replay_object(value["compiler"], COMPILER_RECEIPT_KEYS,
                              "provenance.compiler")
    replay_path(compiler["path"], "provenance.compiler.path")
    lower_hash(compiler["sha256"], "provenance.compiler.sha256")
    version = bounded_string(compiler["version"], "provenance.compiler.version", 128)
    if re.fullmatch(r"13(?:\.[0-9]+){1,2}", version) is None:
        fail("provenance compiler is not GCC 13")
    version_text = replay_text_hash(compiler, "version_text", "version_text_sha256",
                                    "provenance.compiler")
    if not version_text.endswith("\n") or not any(
        name in version_text.lower() for name in ("gcc", "g++")
    ):
        fail("provenance compiler version text differs")
    exact_string(compiler["version_sha256"],
                 sha256_bytes((version + "\n").encode("ascii")),
                 "provenance.compiler.version_sha256")
    validate_interpreter_receipt(value["interpreter"], "provenance.interpreter")
    artifact_paths = {
        "compiler": compiler["path"], "controller": value["controller_path"],
        "controller_interpreter": value["interpreter"]["path"],
        "library": value["library_path"], "worker": value["worker_path"],
    }
    artifacts = replay_object(value["artifacts"], artifact_paths, "provenance.artifacts")
    for name, path in artifact_paths.items():
        validate_stat_receipt(artifacts[name], "provenance.artifacts." + name,
                              path, artifact=True)
        digest = compiler["sha256"] if name == "compiler" else claim[name + "_sha256"]
        exact_string(artifacts[name]["sha256"], digest,
                     "provenance.artifacts." + name + ".sha256")
    replay_equal(artifacts["controller_interpreter"], value["interpreter"]["artifact"],
                 "provenance.interpreter.artifact")
    expected = replay_object(value["expected"], EXPECTED_KEYS, "provenance.expected")
    replay_equal(expected, {
        "compiler_path": compiler["path"], "compiler_sha256": compiler["sha256"],
        "compiler_version": version, "library_path": value["library_path"],
        "native_config_identity_sha256": claim["native_config_identity_sha256"],
        "source_commit": claim["source_commit"],
        "target_identity_sha256": expected["target_identity_sha256"],
    }, "provenance.expected")
    lower_hash(expected["target_identity_sha256"], "provenance.expected.target_identity_sha256")
    native = replay_object(value["native_config"], NATIVE_CONFIG_RECEIPT_KEYS,
                            "provenance.native_config")
    for key in NATIVE_CONFIG_RECEIPT_KEYS - {"runtime_library_path"}:
        lower_hash(native[key], "provenance.native_config." + key)
    for key, wanted in (
        ("runtime_library_path", value["library_path"]),
        ("native_config_identity_sha256", claim["native_config_identity_sha256"]),
        ("target_identity_sha256", expected["target_identity_sha256"]),
        ("stderr_sha256", sha256_bytes(b"")),
    ):
        exact_string(native[key], wanted, "provenance.native_config." + key)
    affinity = replay_object(value["controller_affinity"], {
        "affinity_after", "affinity_before", "controller_cpu",
        "sibling_available_before", "singleton_verified", "target_available_before",
    }, "provenance.controller_affinity")
    exact_int(affinity["controller_cpu"], CONTROLLER_CPU, CONTROLLER_CPU,
              "provenance.controller_affinity.controller_cpu")
    replay_equal(affinity["affinity_after"], [CONTROLLER_CPU],
                 "provenance.controller_affinity.affinity_after")
    before = affinity["affinity_before"]
    if type(before) is not list or not before or len(before) > 1048576:
        fail("provenance inherited controller affinity differs")
    for cpu in before:
        exact_int(cpu, 0, 1048575, "provenance.controller_affinity.affinity_before")
    if before != sorted(set(before)) or not {
        CONTROLLER_CPU, TARGET_CPU, TARGET_SIBLING
    }.issubset(before):
        fail("provenance inherited controller affinity lacks required CPUs")
    for key in ("sibling_available_before", "singleton_verified", "target_available_before"):
        exact_bool(affinity[key], True, "provenance.controller_affinity." + key)
    validate_git_receipt(value["git"], claim["source_commit"], "provenance.git")
    validate_build_receipt(value["build"], value, "provenance.build")
    validate_dynamic_receipt(value["dynamic"], value, "provenance.dynamic")
    validate_r0_receipt(value["r0_before"], "provenance.r0_before")
    validate_health_receipt(value["health_before"], "provenance.health_before")
    if type(value["worker_timed_out"]) is not bool:
        fail("provenance.worker_timed_out is not boolean")
    if value["worker_exit"] is not None:
        exact_int(value["worker_exit"], -255, 255, "provenance.worker_exit")
    if value["worker_started_sha256"] is not None:
        lower_hash(value["worker_started_sha256"], "provenance.worker_started_sha256")

    # postflight assigns receipts in order and stops at its first failure.
    # An invalid run may therefore retain a prefix, including the receipt
    # whose comparison failed.  Successful postflight requires all receipts.
    missing = False
    postflight_failed = any(item.startswith("postflight:") for item in infrastructure)
    for key in PROVENANCE_POSTFLIGHT_FIELDS:
        receipt = value[key]
        if receipt is None:
            missing = True
        elif missing:
            fail("provenance postflight receipts are not a prefix")
        elif type(receipt) is not dict:
            fail("provenance.{} is not an object".format(key))
    if missing and not postflight_failed:
        fail("provenance lacks receipts without a postflight failure")
    if value["health_after"] is not None:
        validate_health_receipt(value["health_after"], "provenance.health_after")
    if value["health_adjudication"] is not None:
        replay_object(value["health_adjudication"], HEALTH_ADJUDICATION_KEYS,
                      "provenance.health_adjudication")
        replay_equal(value["health_adjudication"], validate_health_pair(
            value["health_before"], value["health_after"]),
            "provenance.health_adjudication")
    if value["r0_after"] is not None:
        validate_r0_receipt(value["r0_after"], "provenance.r0_after")
    if value["artifacts_after"] is not None:
        replay_object(value["artifacts_after"], artifact_paths, "provenance.artifacts_after")
        for name, path in artifact_paths.items():
            validate_stat_receipt(value["artifacts_after"][name],
                                  "provenance.artifacts_after." + name,
                                  path, artifact=True)
    if value["interpreter_after"] is not None:
        validate_interpreter_receipt(value["interpreter_after"], "provenance.interpreter_after")
    if value["git_after"] is not None:
        validate_git_receipt(value["git_after"], claim["source_commit"], "provenance.git_after")
    if value["build_after"] is not None:
        validate_build_receipt(value["build_after"], value, "provenance.build_after")
    if value["dynamic_after"] is not None:
        validate_dynamic_receipt(value["dynamic_after"], value, "provenance.dynamic_after")
    if not postflight_failed:
        for before_key, after_key in (
            ("r0_before", "r0_after"), ("artifacts", "artifacts_after"),
            ("interpreter", "interpreter_after"), ("git", "git_after"),
            ("build", "build_after"), ("dynamic", "dynamic_after"),
        ):
            replay_equal(value[before_key], value[after_key], "provenance." + after_key)
    return expected


def replay(args: argparse.Namespace) -> int:
    bundle = exact_absolute_directory(args.bundle, "--bundle")
    if stat.S_IMODE(bundle.stat().st_mode) != 0o700:
        fail("bundle directory mode differs")
    names = sorted(path.name for path in bundle.iterdir())
    if names not in (sorted(OUTPUT_CORE_FILES),
                     sorted(OUTPUT_CORE_FILES + (WORKER_STARTED_FILE,))):
        fail("bundle file roster differs")
    paths = {name: exact_absolute_file(str(bundle / name), name)
             for name in names}
    for name, path in paths.items():
        info = path.stat()
        if stat.S_IMODE(info.st_mode) != 0o400 or info.st_nlink != 1:
            fail("bundle file mode differs for {}".format(name))
    claim = parse_canonical_document(paths["CLAIM"], MAX_PROBE_BYTES)
    summary = parse_canonical_document(paths["summary.json"], MAX_PROBE_BYTES)
    provenance = parse_canonical_document(
        paths["provenance.json"], MAX_PROBE_BYTES)
    complete = parse_canonical_document(paths["COMPLETE"], MAX_PROBE_BYTES)
    with paths["raw.jsonl"].open("rb") as source:
        raw = source.read(MAX_RAW_BYTES + 1)
    with paths["worker.stderr"].open("rb") as source:
        stderr = source.read(MAX_STDERR_BYTES + 1)
    if len(raw) > MAX_RAW_BYTES or len(stderr) > MAX_STDERR_BYTES:
        fail("bundle raw or stderr exceeds its byte cap")
    exact_keys(claim, CLAIM_KEYS, "CLAIM")
    exact_keys(summary, SUMMARY_KEYS, "summary")
    exact_keys(complete, COMPLETE_KEYS, "COMPLETE")
    exact_string(claim["campaign"], CAMPAIGN, "CLAIM.campaign")
    exact_string(claim["schema"], CLAIM_SCHEMA, "CLAIM.schema")
    exact_int(claim["created_unix_ns"], 1, MAX_INT63, "CLAIM.created_unix_ns")
    lower_commit(claim["source_commit"], "CLAIM.source_commit")
    for key in ("controller_interpreter_sha256", "controller_sha256",
                "library_sha256", "native_config_identity_sha256", "worker_sha256"):
        lower_hash(claim[key], "CLAIM." + key)
    exact_string(summary["campaign"], CAMPAIGN, "summary.campaign")
    exact_string(summary["schema"], SUMMARY_SCHEMA, "summary.schema")
    exact_string(provenance.get("schema"), PROVENANCE_SCHEMA,
                 "provenance.schema")
    exact_string(complete["campaign"], CAMPAIGN, "COMPLETE.campaign")
    exact_string(complete["schema"], COMPLETE_SCHEMA, "COMPLETE.schema")
    exact_string(complete["status"], "complete", "COMPLETE.status")
    claim_hash = file_sha256(paths["CLAIM"], MAX_PROBE_BYTES)
    raw_hash = sha256_bytes(raw)
    stderr_hash = sha256_bytes(stderr)
    summary_hash = file_sha256(paths["summary.json"], MAX_PROBE_BYTES)
    provenance_hash = file_sha256(paths["provenance.json"], MAX_PROBE_BYTES)
    for value, expected_value, where in (
        (complete["claim_sha256"], claim_hash, "COMPLETE.claim_sha256"),
        (complete["raw_sha256"], raw_hash, "COMPLETE.raw_sha256"),
        (complete["summary_sha256"], summary_hash, "COMPLETE.summary_sha256"),
        (complete["provenance_sha256"], provenance_hash,
         "COMPLETE.provenance_sha256"),
        (provenance.get("claim_sha256"), claim_hash,
         "provenance.claim_sha256"),
        (provenance.get("raw_sha256"), raw_hash, "provenance.raw_sha256"),
        (provenance.get("stderr_sha256"), stderr_hash,
         "provenance.stderr_sha256"),
        (summary["raw_sha256"], raw_hash, "summary.raw_sha256"),
    ):
        exact_string(value, expected_value, where)
    infrastructure = summary["infrastructure_failures"]
    recorded_reject = summary["reject_gates"]
    recorded_invalid = summary["transcript_invalid_gates"]
    for key in ("infrastructure_failures", "reject_gates", "transcript_invalid_gates"):
        items = summary[key]
        if type(items) is not list or len(items) > EXPECTED_PANEL_COUNT:
            fail("summary.{} is not a bounded list".format(key))
        for item in items:
            bounded_string(item, "summary." + key, MAX_PROBE_BYTES)
    elapsed = summary["elapsed_seconds"]
    if type(elapsed) not in (float, int) or not math.isfinite(elapsed) or elapsed < 0:
        fail("summary.elapsed_seconds is not finite and nonnegative")
    expected = validate_provenance_receipt(provenance, claim, infrastructure)
    marker_hash = None
    if WORKER_STARTED_FILE in paths:
        marker = parse_canonical_document(paths[WORKER_STARTED_FILE], MAX_PROBE_BYTES)
        replay_object(marker, WORKER_STARTED_KEYS, "WORKER_STARTED")
        replay_equal(marker, {
            "campaign": CAMPAIGN, "claim_sha256": claim_hash,
            "native_config_identity_sha256": claim["native_config_identity_sha256"],
            "schema": WORKER_STARTED_SCHEMA, "source_commit": claim["source_commit"],
            "status": "started", "worker_sha256": claim["worker_sha256"],
        }, "WORKER_STARTED")
        marker_hash = file_sha256(paths[WORKER_STARTED_FILE], MAX_PROBE_BYTES)
    replay_equal(provenance["worker_started_sha256"], marker_hash,
                 "provenance.worker_started_sha256")
    if marker_hash is None and "worker_started_missing" not in infrastructure:
        fail("bundle lacks the worker marker without a missing-marker failure")
    if marker_hash is not None and "worker_started_missing" in infrastructure:
        fail("bundle records a missing worker marker that is present")
    timed_out = provenance["worker_timed_out"]
    worker_exit = provenance["worker_exit"]
    if ("worker_outer_deadline" in infrastructure) != timed_out:
        fail("worker timeout receipt and infrastructure disposition differ")
    if ("worker_stderr_nonempty" in infrastructure) != bool(stderr):
        fail("worker stderr receipt and infrastructure disposition differ")
    if worker_exit is None:
        if not any(item.startswith("execution:") for item in infrastructure):
            fail("missing worker exit lacks an execution failure")
    elif ("worker_exit:{}".format(worker_exit) in infrastructure) != (worker_exit != 0):
        fail("worker exit receipt and infrastructure disposition differ")
    statistics = None
    invalid: List[str] = []
    reject: List[str] = []
    transcript_error: Optional[str] = None
    if worker_exit == 0 and not timed_out:
        try:
            first = parse_canonical_line(
                raw.splitlines(keepends=True)[0], "raw config")
            maps_hash = lower_hash(first["runtime_library_maps_sha256"],
                                   "config.runtime_library_maps_sha256")
            config, _, statistics, invalid, reject = parse_transcript(raw, expected, maps_hash)
        except (IndexError, KeyError, ValidationError) as exc:
            transcript_error = str(exc)
        if transcript_error is None:
            # ASLR changes only the runtime map hash between the untimed
            # native preflight and worker capture.  Reconstruct the former
            # config to verify its recorded full-document hash as well.
            native_config = dict(config)
            native_config["runtime_library_maps_sha256"] = provenance[
                "native_config"]["runtime_library_maps_sha256"]
            exact_string(provenance["native_config"]["config_sha256"],
                         sha256_bytes(canonical_bytes(native_config) + b"\n"),
                         "provenance.native_config.config_sha256")
    recorded_errors = [item for item in infrastructure if item.startswith("transcript:")]
    replay_equal(recorded_errors, [] if transcript_error is None else
                 ["transcript:" + transcript_error], "summary transcript failures")
    replay_equal(statistics, summary["statistics"], "recomputed statistics")
    replay_equal(invalid, recorded_invalid, "recomputed transcript invalid gates")
    replay_equal(reject, recorded_reject, "recomputed reject gates")
    outcome = "invalid" if infrastructure or invalid else "reject" if reject else "pass"
    exact_string(summary["outcome"], outcome, "summary.outcome")
    exact_string(complete["outcome"], outcome, "COMPLETE.outcome")
    print(canonical_bytes({
        "infrastructure_failures": infrastructure,
        "outcome": outcome,
        "reject_gates": recorded_reject,
        "statistics": statistics,
        "transcript_invalid_gates": invalid,
    }).decode("ascii"))
    return 0 if outcome == "pass" else 2 if outcome == "reject" else 1


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--selftest", action="store_true")
    modes.add_argument("--native-config-selftest", action="store_true")
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--replay", action="store_true")
    parser.add_argument("--worker")
    parser.add_argument("--library")
    parser.add_argument("--source-root")
    parser.add_argument("--build-root")
    parser.add_argument("--compiler")
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--expected-controller-sha256")
    parser.add_argument("--expected-worker-sha256")
    parser.add_argument("--expected-library-sha256")
    parser.add_argument("--bundle")
    args = parser.parse_args(argv)
    run_fields = (
        "worker", "library", "source_root", "build_root", "compiler",
        "expected_source_commit", "expected_controller_sha256",
        "expected_worker_sha256", "expected_library_sha256",
    )
    native_fields = (
        "worker", "library", "compiler", "expected_source_commit",
    )
    if args.run:
        missing = [field for field in run_fields if getattr(args, field) is None]
        if missing:
            parser.error("--run requires {}".format(
                ", ".join("--" + item.replace("_", "-") for item in missing)))
        if args.bundle is not None:
            parser.error("--run does not accept replay inputs")
    elif args.native_config_selftest:
        missing = [
            field for field in native_fields if getattr(args, field) is None
        ]
        if missing:
            parser.error("--native-config-selftest requires {}".format(
                ", ".join("--" + item.replace("_", "-") for item in missing)))
        forbidden = set(run_fields) - set(native_fields)
        if (
            any(getattr(args, field) is not None for field in forbidden) or
            args.bundle is not None
        ):
            parser.error("--native-config-selftest accepts only its four inputs")
    elif args.replay:
        if args.bundle is None:
            parser.error("--replay requires --bundle")
        if any(getattr(args, field) is not None for field in run_fields):
            parser.error("--replay does not accept run authority arguments")
    elif any(getattr(args, field) is not None for field in run_fields) or \
            args.bundle is not None:
        parser.error("--selftest accepts no workload or artifact arguments")
    return args


def main(argv: Sequence[str]) -> int:
    try:
        args = parse_args(argv)
        if args.selftest:
            return selftest()
        if args.native_config_selftest:
            return native_config_selftest(args)
        if args.replay:
            return replay(args)
        return run_once(args)
    except ValidationError as exc:
        print("invalid: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
