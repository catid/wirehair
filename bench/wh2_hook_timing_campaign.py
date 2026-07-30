#!/usr/bin/env python3
"""Run the source-authenticated paired WH2 hook-overhead campaign.

This protocol is deliberately measurement-only.  It records hook-on/off
effects and A/A floors, but it cannot make a production-promotion decision.
"""

import argparse
import base64
import ctypes
import dataclasses
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import selectors
import signal
import stat
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


SCHEMA = "wirehair.wh2.hook-timing-campaign.v1"
RAW_SCHEMA = "wirehair.wh2.hook-timing-campaign.raw.v1"
BUILD_SCHEMA = "wirehair.wh2.hook-build-matrix.v1"
TIMING_SCHEMA = "wirehair.wh2.hook-timing.v1"
BUILD_TARGET = "wirehair_v2_hook_timing"
HOOK_MACRO = "WIREHAIR_V2_ENABLE_TEST_HOOKS"

TREATMENTS = ("hook-on-a", "hook-on-b", "hook-off-a", "hook-off-b")
TREATMENT_STATE = {
    "hook-on-a": 1,
    "hook-on-b": 1,
    "hook-off-a": 0,
    "hook-off-b": 0,
}
WILLIAMS_ROWS = (
    ("hook-on-a", "hook-on-b", "hook-off-b", "hook-off-a"),
    ("hook-on-b", "hook-off-a", "hook-on-a", "hook-off-b"),
    ("hook-off-a", "hook-off-b", "hook-on-b", "hook-on-a"),
    ("hook-off-b", "hook-on-a", "hook-off-a", "hook-on-b"),
)

PRIMARY_CELLS = (
    (64, 2),
    (100, 1280),
    (101, 1280),
    (512, 2),
    (513, 2),
    (1000, 1280),
    (9999, 1280),
    (10000, 1280),
    (32000, 1280),
    (64000, 1280),
    (64, 4096),
    (1000, 4096),
    (10000, 4096),
    (32000, 4096),
)
SCOPES = ("row", "encoder", "decoder-feed", "decoder-full", "direct")
SCHEDULES = (
    "iid",
    "burst",
    "permutation",
    "systematic-first",
    "repair-only",
    "adversarial",
)
PRIMARY_SCHEDULE = "repair-only"
PRIMARY_LOSS_PPM = 350000
TIMED_ANCHORS = (
    {"schedule": "iid", "K": 64, "block_bytes": 4096},
    {"schedule": "burst", "K": 10000, "block_bytes": 1280},
    {"schedule": "adversarial", "K": 64000, "block_bytes": 1280},
)

SEMANTIC_SCOPE = "semantic"
TARGET_BATCH_NS = 100 * 1000 * 1000
WARMUP_BLOCKS = 2
MEASURED_BLOCKS = 16
WARMUPS_PER_TREATMENT = 4
OBSERVATIONS_PER_TREATMENT = 32
MAX_INNER_REPS = 1000000
MAX_CALIBRATION_STEPS = 24
MAX_BUILD_RECEIPT_BYTES = 256 * 1024 * 1024
MAX_BINARY_BYTES = 512 * 1024 * 1024
MAX_STREAM_BYTES = 8 * 1024 * 1024
MAX_RAW_RECORD_BYTES = 24 * 1024 * 1024
MAX_JSON_NESTING = 64
PROCESS_REAP_GRACE_SECONDS = 1.0
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1
F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
F_SEAL_EXEC = getattr(fcntl, "F_SEAL_EXEC", 0x0020)
EXECUTION_SEALS = (
    F_SEAL_SEAL | F_SEAL_SHRINK | F_SEAL_GROW | F_SEAL_WRITE | F_SEAL_EXEC)
EXECUTION_MODE = 0o500
MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
MANAGED_SIGNALS = (signal.SIGINT, signal.SIGTERM)

HEX64 = re.compile(r"[0-9a-f]{64}")
BUILD_ID = re.compile(r"[0-9a-f]+")

HEADER_KEYS = (
    "record", "schema", "profile", "hooks_compiled", "K", "block_bytes",
    "construction_seed", "loss_seed", "loss_ppm", "schedule",
    "scope_request", "warmup_reps", "inner_reps", "max_working_mib",
    "context_sha256", "cpu", "clock", "start_barrier",
)
SEMANTIC_KEYS = (
    "record", "schema", "profile", "status", "encoder_result",
    "decoder_result", "recover_result", "direct_result", "trace_packets",
    "delivered_packets", "overhead_packets", "decoder_solve_attempts",
    "profile_sha256", "system_sha256", "coefficients_sha256",
    "trace_sha256", "row_sha256", "message_sha256", "payload_sha256",
    "intermediate_sha256", "recovered_sha256", "encoder_stats_sha256",
    "decoder_stats_sha256", "direct_stats_sha256", "semantic_sha256",
)
TIMING_KEYS = (
    "record", "schema", "profile", "hooks_compiled", "scope", "lifecycle",
    "semantic_sha256", "unit", "clock", "warmup_reps", "measured_reps",
    "inner_reps", "work_items_per_rep", "elapsed_ns", "min_ns", "max_ns",
    "minor_faults", "major_faults", "voluntary_context_switches",
    "involuntary_context_switches", "sink", "result",
)
DONE_KEYS = (
    "record", "schema", "status", "records_before_done", "stream_sha256",
)
LIFECYCLE = {
    "row": "caller-row-buffer-reuse-v1",
    "encoder":
        "fresh-first-then-transactional-reinitialize-including-prior-release-v1",
    "decoder-feed": "distinct-preinitialized-endpoints-v1",
    "decoder-full":
        "fresh-first-then-transactional-reinitialize-including-prior-release-v1",
    "direct": "fresh-first-then-transactional-output-reuse-v1",
}

# Two-sided Student-t 0.975 quantiles.  The protocol uses n=16 (df=15);
# the complete small table keeps the statistics helper independently testable.
T_CRITICAL_975 = (
    None,
    12.706204736, 4.302652730, 3.182446305, 2.776445105,
    2.570581836, 2.446911851, 2.364624252, 2.306004135,
    2.262157163, 2.228138852, 2.200985160, 2.178812830,
    2.160368656, 2.144786688, 2.131449546, 2.119905299,
    2.109815578, 2.100922040, 2.093024054, 2.085963447,
    2.079613845, 2.073873068, 2.068657610, 2.063898562,
    2.059538553, 2.055529439, 2.051830516, 2.048407142,
    2.045229642, 2.042272456,
)


class CampaignError(RuntimeError):
    """A protocol or execution invariant failed."""


class CampaignInterrupted(CampaignError):
    """The campaign was interrupted after its children were reaped."""


def _int(value: Any, name: str, minimum: int = 0,
         maximum: Optional[int] = None) -> int:
    if type(value) is not int or value < minimum:
        raise CampaignError("{} must be an integer >= {}".format(
            name, minimum))
    if maximum is not None and value > maximum:
        raise CampaignError("{} exceeds {}".format(name, maximum))
    return value


def _string(value: Any, name: str) -> str:
    if not isinstance(value, str) or "\0" in value:
        raise CampaignError("{} must be a NUL-free string".format(name))
    return value


def _sha256(value: Any, name: str) -> str:
    value = _string(value, name)
    if HEX64.fullmatch(value) is None:
        raise CampaignError("{} must be lower-case SHA-256".format(name))
    return value


def _build_id(value: Any, name: str) -> str:
    value = _string(value, name)
    if (BUILD_ID.fullmatch(value) is None or len(value) % 2 != 0 or
            len(value) < 8):
        raise CampaignError("{} is not a lower-case build ID".format(name))
    return value


def canonical_json_bytes(value: Any) -> bytes:
    try:
        return (
            json.dumps(
                value, sort_keys=True, separators=(",", ":"),
                ensure_ascii=True, allow_nan=False,
            ) + "\n"
        ).encode("ascii")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise CampaignError("value is not canonical JSON: {}".format(error))


def _reject_float(unused: str) -> Any:
    raise CampaignError("floating-point JSON input is forbidden")


def _reject_constant(value: str) -> Any:
    raise CampaignError("non-finite JSON token is forbidden: {}".format(value))


def _bounded_json_int(text: str) -> int:
    digits = text[1:] if text.startswith("-") else text
    if len(digits) > 20:
        raise CampaignError("JSON integer exceeds uint64 magnitude")
    value = int(text)
    if value < -UINT64_MAX or value > UINT64_MAX:
        raise CampaignError("JSON integer exceeds uint64 magnitude")
    return value


def _validate_json_nesting(text: str) -> None:
    depth = 0
    in_string = False
    escaped = False
    for character in text:
        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False
            continue
        if character == '"':
            in_string = True
        elif character in "[{":
            depth += 1
            if depth > MAX_JSON_NESTING:
                raise CampaignError("JSON nesting exceeds protocol limit")
        elif character in "]}":
            depth -= 1
            if depth < 0:
                raise CampaignError("JSON closing delimiter is unmatched")


def _unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    output = {}
    for key, value in pairs:
        if key in output:
            raise CampaignError("duplicate JSON key: {}".format(key))
        output[key] = value
    return output


def strict_json_document(payload: bytes, description: str,
                         require_sorted_canonical: bool = True) -> Any:
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as error:
        raise CampaignError("{} is not ASCII JSON: {}".format(
            description, error))
    _validate_json_nesting(text)
    try:
        value = json.loads(
            text,
            object_pairs_hook=_unique_object,
            parse_int=_bounded_json_int,
            parse_float=_reject_float,
            parse_constant=_reject_constant,
        )
    except (
            json.JSONDecodeError, CampaignError, ValueError,
            RecursionError) as error:
        raise CampaignError("invalid {}: {}".format(description, error))
    if require_sorted_canonical:
        if canonical_json_bytes(value) != payload:
            raise CampaignError("{} is not exact canonical JSON".format(
                description))
    return value


def strict_json_line(line: bytes, description: str) -> Dict[str, Any]:
    if not line.endswith(b"\n") or line in (b"\n", b""):
        raise CampaignError("{} is not one terminated JSON line".format(
            description))
    value = strict_json_document(
        line[:-1], description, require_sorted_canonical=False)
    if not isinstance(value, dict):
        raise CampaignError("{} is not a JSON object".format(description))
    try:
        encoded = (
            json.dumps(
                value, sort_keys=False, separators=(",", ":"),
                ensure_ascii=True, allow_nan=False,
            ) + "\n"
        ).encode("ascii")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise CampaignError("{} cannot be canonicalized: {}".format(
            description, error))
    if encoded != line:
        raise CampaignError("{} is not exact compact JSON".format(description))
    return value


def _stable_file_bytes(path: Path, limit: int,
                       description: str) -> Tuple[bytes, Dict[str, Any]]:
    path = Path(os.path.abspath(os.fspath(path)))
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
    except OSError as error:
        raise CampaignError("cannot open {}: {}".format(description, error))
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise CampaignError("{} is not a regular file".format(description))
        if before.st_size < 0 or before.st_size > limit:
            raise CampaignError("{} exceeds byte limit".format(description))
        chunks = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                raise CampaignError("{} was truncated while read".format(
                    description))
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise CampaignError("{} grew while read".format(description))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError("{} disappeared: {}".format(description, error))
    stable = ("st_dev", "st_ino", "st_mode", "st_size",
              "st_mtime_ns", "st_ctime_ns")
    if (any(getattr(before, key) != getattr(after, key) for key in stable) or
            any(getattr(after, key) != getattr(current, key)
                for key in stable)):
        raise CampaignError("{} changed while read".format(description))
    payload = b"".join(chunks)
    return payload, {
        "path": str(path),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": stat.S_IMODE(after.st_mode),
    }


def _fsync_directory(path: Path, description: str) -> None:
    directory_flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_DIRECTORY"):
        directory_flags |= os.O_DIRECTORY
    try:
        descriptor = os.open(str(path), directory_flags)
    except OSError as error:
        raise CampaignError(
            "cannot open {} directory: {}".format(description, error))
    sync_error = None
    try:
        os.fsync(descriptor)
    except OSError as error:
        sync_error = error
    try:
        os.close(descriptor)
    except OSError as error:
        if sync_error is None:
            sync_error = error
    if sync_error is not None:
        raise CampaignError(
            "cannot sync {} directory: {}".format(
                description, sync_error))


@dataclasses.dataclass(frozen=True)
class BinaryBinding:
    treatment: str
    hook_state: int
    replicate: str
    path: str
    bytes: int
    sha256: str
    build_id: str
    device: int
    inode: int
    mode: int
    execution_fd: int


def _sealed_executable_memfd(name: str, payload: bytes) -> int:
    try:
        encoded_name = name.encode("ascii")
    except UnicodeEncodeError:
        raise CampaignError("executable memfd name is not ASCII")
    if not encoded_name or len(encoded_name) > 249 or b"\0" in encoded_name:
        raise CampaignError("executable memfd name is invalid")
    flags = MFD_CLOEXEC | MFD_ALLOW_SEALING
    try:
        if hasattr(os, "memfd_create"):
            descriptor = os.memfd_create(name, flags)
        else:
            libc = ctypes.CDLL(None, use_errno=True)
            function = getattr(libc, "memfd_create", None)
            if function is None:
                raise CampaignError("sealed executable memfd is unavailable")
            function.argtypes = (ctypes.c_char_p, ctypes.c_uint)
            function.restype = ctypes.c_int
            descriptor = function(encoded_name, flags)
            if descriptor < 0:
                error_number = ctypes.get_errno()
                raise OSError(error_number, os.strerror(error_number))
    except OSError as error:
        raise CampaignError("cannot create executable memfd: {}".format(error))
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("short memfd write")
            view = view[written:]
        os.fchmod(descriptor, EXECUTION_MODE)
        fcntl.fcntl(descriptor, F_ADD_SEALS, EXECUTION_SEALS)
        actual_seals = fcntl.fcntl(descriptor, F_GET_SEALS)
        if actual_seals & EXECUTION_SEALS != EXECUTION_SEALS:
            raise CampaignError("executable memfd seal set is incomplete")
        os.lseek(descriptor, 0, os.SEEK_SET)
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def close_execution_bindings(
        bindings: Dict[str, BinaryBinding]) -> None:
    closed = set()
    for binding in list(bindings.values()):
        if binding.execution_fd in closed:
            continue
        try:
            os.close(binding.execution_fd)
        except OSError:
            pass
        closed.add(binding.execution_fd)
    bindings.clear()


def _require_mapping(value: Any, name: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise CampaignError("{} must be an object".format(name))
    return value


def authenticate_build_receipt(
        path: Path, expected_sha256: str,
        binding_sink: Optional[Dict[str, BinaryBinding]] = None) -> Tuple[
            Dict[str, Any], Dict[str, BinaryBinding], Dict[str, Any]]:
    expected_sha256 = _sha256(expected_sha256, "build receipt SHA-256")
    payload, receipt_binding = _stable_file_bytes(
        path, MAX_BUILD_RECEIPT_BYTES, "build receipt")
    if receipt_binding["sha256"] != expected_sha256:
        raise CampaignError("build receipt SHA-256 mismatch")
    receipt = strict_json_document(payload, "build receipt")
    receipt = _require_mapping(receipt, "build receipt")
    if receipt.get("schema") != BUILD_SCHEMA:
        raise CampaignError("wrong build receipt schema")
    if receipt.get("timing_executed") is not False:
        raise CampaignError("build receipt already claims timing execution")
    macro = _require_mapping(receipt.get("hook_macro"), "hook_macro")
    if macro != {
            "name": HOOK_MACRO,
            "off_contract": "macro-absent",
            "on_token": "-D{}=1".format(HOOK_MACRO)}:
        raise CampaignError("build receipt hook macro contract changed")
    configuration = _require_mapping(
        receipt.get("configuration"), "configuration")
    expected_configuration = {
        "target": BUILD_TARGET,
        "binary_relative": "codec/wirehair_v2_hook_timing",
        "build_type": "Release",
        "build_order": [
            "hook-on-a", "hook-off-a", "hook-on-b", "hook-off-b"],
        "generator": "cmake-default",
        "definitions": [],
        "flags": {
            "c_flags": "",
            "cxx_flags": "",
            "exe_linker_flags": "",
            "module_linker_flags": "",
            "shared_linker_flags": "",
        },
        "build_shared": "OFF",
        "static_pic": "ON",
        "march_native": "OFF",
        "lto": "OFF",
        "pgo": "OFF",
    }
    for key, value in expected_configuration.items():
        if configuration.get(key) != value:
            raise CampaignError(
                "build configuration changed at {}".format(key))
    _int(configuration.get("jobs"), "build jobs", 1, 256)
    comparisons = _require_mapping(
        receipt.get("comparability"), "comparability")
    required_true = (
        "cache_normalized_sha256_identical",
        "compile_commands_normalized_sha256_identical",
        "link_normalized_sha256_identical",
        "toolchains_identical",
        "hook_on_ab_binary_identical",
        "hook_off_ab_binary_identical",
    )
    if any(comparisons.get(key) is not True for key in required_true):
        raise CampaignError("build comparability proof is incomplete")
    if comparisons.get("hook_on_off_binary_identical") is not False:
        raise CampaignError("hook-on/off binaries are not distinct")
    if comparisons.get("only_allowed_compile_difference") != \
            "-D{}=1".format(HOOK_MACRO):
        raise CampaignError("build compile-difference proof changed")
    builds = _require_mapping(receipt.get("builds"), "builds")
    if set(builds) != set(TREATMENTS):
        raise CampaignError("build receipt does not contain four treatments")

    tools = _require_mapping(receipt.get("tools"), "tools")
    if set(tools) != {"python", "git", "cmake", "readelf"}:
        raise CampaignError("build receipt tool set changed")
    authenticated_tools = {}
    for tool_name in ("python", "readelf"):
        tool = _require_mapping(tools[tool_name], tool_name)
        tool_path = Path(_string(
            tool.get("resolved_path"),
            "{} resolved path".format(tool_name)))
        expected_tool_bytes = _int(
            tool.get("bytes"), "{} bytes".format(tool_name), 1,
            MAX_BINARY_BYTES)
        expected_tool_hash = _sha256(
            tool.get("sha256"), "{} SHA-256".format(tool_name))
        unused_payload, actual_tool = _stable_file_bytes(
            tool_path, MAX_BINARY_BYTES, "{} tool".format(tool_name))
        if (actual_tool["bytes"] != expected_tool_bytes or
                actual_tool["sha256"] != expected_tool_hash):
            raise CampaignError(
                "{} tool changed after build".format(tool_name))
        authenticated_tools[tool_name] = {
            "path": str(tool_path),
            "bytes": expected_tool_bytes,
            "sha256": expected_tool_hash,
        }
    current_python = Path(sys.executable).resolve()
    unused_python, current_python_binding = _stable_file_bytes(
        current_python, MAX_BINARY_BYTES, "running Python")
    if (
            str(current_python) != authenticated_tools["python"]["path"] or
            current_python_binding["bytes"] !=
                authenticated_tools["python"]["bytes"] or
            current_python_binding["sha256"] !=
                authenticated_tools["python"]["sha256"]):
        raise CampaignError(
            "campaign Python does not match authenticated build Python")

    bindings = {} if binding_sink is None else binding_sink
    if bindings:
        raise CampaignError("binary binding sink must start empty")
    identities = set()
    try:
        for treatment in TREATMENTS:
            build = _require_mapping(builds[treatment], treatment)
            state_text = "on" if TREATMENT_STATE[treatment] else "off"
            replicate = treatment[-1]
            if (build.get("name") != treatment or
                    build.get("hook_state") != state_text or
                    build.get("replicate") != replicate):
                raise CampaignError("build identity mismatch for {}".format(
                    treatment))
            binary = _require_mapping(
                build.get("binary"), "{} binary".format(treatment))
            expected_bytes = _int(
                binary.get("bytes"), "{} binary bytes".format(treatment), 1,
                MAX_BINARY_BYTES)
            expected_hash = _sha256(
                binary.get("sha256"), "{} binary SHA-256".format(treatment))
            expected_mode = _int(
                binary.get("mode"), "{} binary mode".format(treatment),
                0, 0o7777)
            build_id_object = _require_mapping(
                binary.get("build_id"), "{} build ID".format(treatment))
            expected_build_id = _build_id(
                build_id_object.get("value"),
                "{} build ID value".format(treatment))
            binary_path = Path(_string(
                binary.get("path"), "{} binary path".format(treatment)))
            if not binary_path.is_absolute():
                raise CampaignError("binary path is not absolute")
            binary_payload, actual = _stable_file_bytes(
                binary_path, MAX_BINARY_BYTES,
                "{} binary".format(treatment))
            if (actual["bytes"] != expected_bytes or
                    actual["sha256"] != expected_hash or
                    actual["mode"] != expected_mode):
                raise CampaignError("binary binding mismatch for {}".format(
                    treatment))
            if expected_mode & 0o7000 or expected_mode & stat.S_IXUSR == 0:
                raise CampaignError(
                    "binary mode must be ordinary and owner-executable")
            identity = (actual["device"], actual["inode"])
            if identity in identities:
                raise CampaignError(
                    "treatments must use distinct binary inodes")
            identities.add(identity)
            execution_fd = _sealed_executable_memfd(
                "wh2-{}".format(treatment), binary_payload)
            bindings[treatment] = BinaryBinding(
                treatment=treatment,
                hook_state=TREATMENT_STATE[treatment],
                replicate=replicate,
                path=str(binary_path),
                bytes=expected_bytes,
                sha256=expected_hash,
                build_id=expected_build_id,
                device=actual["device"],
                inode=actual["inode"],
                mode=expected_mode,
                execution_fd=execution_fd,
            )
    except BaseException:
        close_execution_bindings(bindings)
        raise

    try:
        for state in ("on", "off"):
            left = bindings["hook-{}-a".format(state)]
            right = bindings["hook-{}-b".format(state)]
            if (left.bytes, left.sha256, left.build_id, left.mode) != (
                    right.bytes, right.sha256, right.build_id, right.mode):
                raise CampaignError(
                    "same-state A/B binaries are not identical")
        if len({binding.mode for binding in bindings.values()}) != 1:
            raise CampaignError(
                "all treatment binary modes must be identical")
        on = bindings["hook-on-a"]
        off = bindings["hook-off-a"]
        if (on.bytes, on.sha256) == (off.bytes, off.sha256):
            raise CampaignError("actual hook-on/off binaries are identical")
    except BaseException:
        close_execution_bindings(bindings)
        raise

    receipt_binding["authenticated_tools"] = authenticated_tools
    return receipt, bindings, receipt_binding


@dataclasses.dataclass(frozen=True)
class InvocationExpected:
    treatment: str
    K: int
    block_bytes: int
    construction_seed: int
    loss_seed: int
    loss_ppm: int
    schedule: str
    scope: str
    inner_reps: int
    max_working_mib: int
    context_sha256: str
    cpu: int


@dataclasses.dataclass(frozen=True)
class ParsedTimingStream:
    header: Dict[str, Any]
    semantic: Dict[str, Any]
    timing: Optional[Dict[str, Any]]
    done: Dict[str, Any]
    status: str


def _exact_keys(value: Dict[str, Any], keys: Tuple[str, ...],
                description: str) -> None:
    if tuple(value.keys()) != keys:
        raise CampaignError("{} keys/order changed".format(description))


def _validate_hash_or_na(value: Any, name: str, allow_na: bool) -> str:
    if allow_na and value == "not_applicable":
        return value
    return _sha256(value, name)


def parse_timing_stream(
        payload: bytes, expected: InvocationExpected) -> ParsedTimingStream:
    if not payload or len(payload) > MAX_STREAM_BYTES:
        raise CampaignError("timing stdout is empty or oversized")
    if not payload.endswith(b"\n") or b"\r" in payload or b"\0" in payload:
        raise CampaignError("timing stdout framing is not canonical")
    raw_lines = payload.splitlines(keepends=True)
    if len(raw_lines) < 3:
        raise CampaignError("timing stream is truncated")
    records = [
        strict_json_line(line, "timing record {}".format(index))
        for index, line in enumerate(raw_lines)
    ]
    header, semantic, done = records[0], records[1], records[-1]
    _exact_keys(header, HEADER_KEYS, "header")
    _exact_keys(semantic, SEMANTIC_KEYS, "semantic")
    _exact_keys(done, DONE_KEYS, "done")
    for name, minimum, maximum in (
            ("hooks_compiled", 0, 1),
            ("K", 2, 64000),
            ("block_bytes", 1, 0x7fffffff),
            ("construction_seed", 0, UINT64_MAX),
            ("loss_seed", 0, UINT64_MAX),
            ("loss_ppm", 0, 999999),
            ("warmup_reps", 0, 1000),
            ("inner_reps", 1, MAX_INNER_REPS),
            ("max_working_mib", 1, 1048576),
            ("cpu", 0, 1048575)):
        _int(header.get(name), "header {}".format(name), minimum, maximum)
    if header != {
            "record": "header",
            "schema": TIMING_SCHEMA,
            "profile": "dispatch-v1",
            "hooks_compiled": TREATMENT_STATE[expected.treatment],
            "K": expected.K,
            "block_bytes": expected.block_bytes,
            "construction_seed": expected.construction_seed,
            "loss_seed": expected.loss_seed,
            "loss_ppm": expected.loss_ppm,
            "schedule": expected.schedule,
            "scope_request": expected.scope,
            "warmup_reps": 0,
            "inner_reps": expected.inner_reps,
            "max_working_mib": expected.max_working_mib,
            "context_sha256": expected.context_sha256,
            "cpu": expected.cpu,
            "clock": "CLOCK_MONOTONIC",
            "start_barrier":
                "none" if expected.scope == SEMANTIC_SCOPE
                else "ready-go-pipe-v1"}:
        raise CampaignError("timing header does not match invocation")
    if (semantic.get("record") != "semantic" or
            semantic.get("schema") != TIMING_SCHEMA or
            semantic.get("profile") != "dispatch-v1"):
        raise CampaignError("semantic record identity changed")
    status = semantic.get("status")
    if status not in ("success", "weak-root"):
        raise CampaignError("unknown semantic status")
    _int(semantic.get("encoder_result"), "encoder_result", 0, 11)
    _int(semantic.get("trace_packets"), "trace_packets", 0, 200000)
    _int(semantic.get("delivered_packets"), "delivered_packets", 0, 200000)
    _int(semantic.get("decoder_solve_attempts"),
         "decoder_solve_attempts", 0, 1000000)
    common_hashes = (
        "profile_sha256", "system_sha256", "coefficients_sha256",
        "trace_sha256", "row_sha256", "message_sha256", "semantic_sha256",
    )
    for name in common_hashes:
        _sha256(semantic.get(name), name)
    late_hashes = (
        "payload_sha256", "intermediate_sha256", "recovered_sha256",
        "encoder_stats_sha256", "decoder_stats_sha256",
        "direct_stats_sha256",
    )
    if status == "success":
        for name in (
                "decoder_result", "recover_result", "direct_result"):
            _int(semantic.get(name), name, 0, 11)
        if (semantic.get("encoder_result") != 0 or
                semantic.get("decoder_result") != 0 or
                semantic.get("recover_result") != 0 or
                semantic.get("direct_result") != 0):
            raise CampaignError("successful semantic result is nonzero")
        _int(
            semantic.get("overhead_packets"), "overhead_packets",
            0, 200000)
        if (
                semantic["delivered_packets"] !=
                expected.K + semantic["overhead_packets"]):
            raise CampaignError("successful overhead accounting changed")
        for name in late_hashes:
            _validate_hash_or_na(semantic.get(name), name, False)
    else:
        if (semantic.get("encoder_result") not in (1, 4, 8) or
                semantic.get("decoder_result") is not None or
                semantic.get("recover_result") is not None or
                semantic.get("direct_result") is not None or
                semantic.get("delivered_packets") != 0 or
                semantic.get("overhead_packets") is not None or
                semantic.get("decoder_solve_attempts") != 0):
            raise CampaignError("weak-root structure changed")
        for name in late_hashes:
            if semantic.get(name) != "not_applicable":
                raise CampaignError("weak-root {} is applicable".format(name))

    timing_records = records[2:-1]
    timing = None
    if expected.scope == SEMANTIC_SCOPE or status == "weak-root":
        if timing_records:
            raise CampaignError("semantic/weak invocation emitted timing")
    else:
        if len(timing_records) != 1:
            raise CampaignError("single-scope invocation must emit one timing")
        timing = timing_records[0]
        _exact_keys(timing, TIMING_KEYS, "timing")
        for name, minimum, maximum in (
                ("hooks_compiled", 0, 1),
                ("warmup_reps", 0, 1000),
                ("measured_reps", 1, 1000000),
                ("inner_reps", 1, MAX_INNER_REPS),
                ("work_items_per_rep", 1, UINT64_MAX),
                ("elapsed_ns", 1, UINT64_MAX),
                ("min_ns", 1, UINT64_MAX),
                ("max_ns", 1, UINT64_MAX),
                ("minor_faults", 0, UINT64_MAX),
                ("major_faults", 0, UINT64_MAX),
                ("voluntary_context_switches", 0, UINT64_MAX),
                ("involuntary_context_switches", 0, UINT64_MAX),
                ("sink", 0, UINT64_MAX),
                ("result", 0, 11)):
            _int(timing.get(name), "timing {}".format(name), minimum, maximum)
        if (timing.get("record") != "timing" or
                timing.get("schema") != TIMING_SCHEMA or
                timing.get("profile") != "dispatch-v1" or
                timing.get("hooks_compiled") !=
                TREATMENT_STATE[expected.treatment] or
                timing.get("scope") != expected.scope or
                timing.get("lifecycle") != LIFECYCLE[expected.scope] or
                timing.get("semantic_sha256") !=
                semantic["semantic_sha256"] or
                timing.get("unit") != "ns" or
                timing.get("clock") != "CLOCK_MONOTONIC" or
                timing.get("warmup_reps") != 0 or
                timing.get("measured_reps") != 1 or
                timing.get("inner_reps") != expected.inner_reps or
                timing.get("result") != 0):
            raise CampaignError("timing metadata changed")
        work = timing["work_items_per_rep"]
        expected_work = (
            expected.K if expected.scope == "encoder"
            else semantic["delivered_packets"])
        if work != expected_work:
            raise CampaignError("timing work count does not match semantic")
        elapsed = timing["elapsed_ns"]
        if (timing.get("min_ns") != elapsed or
                timing.get("max_ns") != elapsed):
            raise CampaignError("single-batch min/max changed")
    expected_before_done = len(records) - 1
    _int(
        done.get("records_before_done"), "records_before_done",
        2, len(records))
    prefix = b"".join(raw_lines[:-1])
    if done != {
            "record": "done",
            "schema": TIMING_SCHEMA,
            "status": status,
            "records_before_done": expected_before_done,
            "stream_sha256": hashlib.sha256(prefix).hexdigest()}:
        raise CampaignError("done record/hash does not authenticate stream")
    if expected_before_done != 2 + len(timing_records):
        raise CampaignError("done record count is inconsistent")
    return ParsedTimingStream(
        header=header, semantic=semantic, timing=timing,
        done=done, status=status)


def williams_block(block_index: int) -> List[Dict[str, Any]]:
    """Return four simultaneous rounds from one adjacent Williams-row pair."""
    _int(block_index, "block index", 0)
    first_row = block_index % 4
    rows = (first_row, (first_row + 1) % 4)
    swap = (block_index // 4) % 2 != 0
    cpu_rows = (rows[1], rows[0]) if swap else rows
    rounds = []
    for round_index in range(4):
        rounds.append({
            "round": round_index,
            "cpu_a_treatment": WILLIAMS_ROWS[cpu_rows[0]][round_index],
            "cpu_b_treatment": WILLIAMS_ROWS[cpu_rows[1]][round_index],
            "cpu_a_row": cpu_rows[0],
            "cpu_b_row": cpu_rows[1],
        })
    return rounds


def validate_frozen_design() -> None:
    if (len(PRIMARY_CELLS) != 14 or len(set(PRIMARY_CELLS)) != 14 or
            len(SCOPES) != 5 or len(set(SCOPES)) != 5 or
            len(SCHEDULES) != 6 or len(set(SCHEDULES)) != 6 or
            len(TIMED_ANCHORS) != 3 or
            WARMUP_BLOCKS != 2 or WARMUPS_PER_TREATMENT != 4 or
            MEASURED_BLOCKS != 16 or
            OBSERVATIONS_PER_TREATMENT != 32):
        raise CampaignError("frozen campaign dimensions changed")
    expected_cells = {
        (64, 2), (100, 1280), (101, 1280), (512, 2), (513, 2),
        (1000, 1280), (9999, 1280), (10000, 1280), (32000, 1280),
        (64000, 1280), (64, 4096), (1000, 4096),
        (10000, 4096), (32000, 4096),
    }
    if set(PRIMARY_CELLS) != expected_cells:
        raise CampaignError("frozen primary cells changed")
    if WILLIAMS_ROWS != (
            ("hook-on-a", "hook-on-b", "hook-off-b", "hook-off-a"),
            ("hook-on-b", "hook-off-a", "hook-on-a", "hook-off-b"),
            ("hook-off-a", "hook-off-b", "hook-on-b", "hook-on-a"),
            ("hook-off-b", "hook-on-a", "hook-off-a", "hook-on-b")):
        raise CampaignError("frozen Williams rows changed")
    if tuple(
            (row["schedule"], row["K"], row["block_bytes"])
            for row in TIMED_ANCHORS) != (
                ("iid", 64, 4096),
                ("burst", 10000, 1280),
                ("adversarial", 64000, 1280)):
        raise CampaignError("frozen timed anchors changed")
    for row in WILLIAMS_ROWS:
        if len(row) != 4 or set(row) != set(TREATMENTS):
            raise CampaignError("Williams row is not a treatment permutation")
    counts = {
        treatment: {"cpu_a": 0, "cpu_b": 0, "total": 0}
        for treatment in TREATMENTS
    }
    pair_kinds = {"on-aa": 0, "off-aa": 0, "cross": 0}
    pair_kinds_by_period = [
        {"on-aa": 0, "off-aa": 0, "cross": 0}
        for unused_period in range(4)
    ]
    exact_pair_orientations: Dict[
        Tuple[str, str], Dict[Tuple[str, str], int]] = {}
    cpu_row_counts = {
        cpu: {row: 0 for row in range(4)}
        for cpu in ("cpu_a", "cpu_b")}
    cpu_period_treatments = {
        cpu: [
            {treatment: 0 for treatment in TREATMENTS}
            for unused_period in range(4)
        ]
        for cpu in ("cpu_a", "cpu_b")
    }
    carryovers = {
        cpu: {} for cpu in ("cpu_a", "cpu_b")}
    for block in range(MEASURED_BLOCKS):
        block_counts = {treatment: 0 for treatment in TREATMENTS}
        round_rows = williams_block(block)
        cpu_sequences = {"cpu_a": [], "cpu_b": []}
        for cpu in ("cpu_a", "cpu_b"):
            row_ids = {
                round_row["{}_row".format(cpu)]
                for round_row in round_rows}
            if len(row_ids) != 1:
                raise CampaignError("Williams CPU row changes within block")
            cpu_row_counts[cpu][row_ids.pop()] += 1
        for round_row in round_rows:
            pair = (
                round_row["cpu_a_treatment"],
                round_row["cpu_b_treatment"])
            states = tuple(TREATMENT_STATE[item] for item in pair)
            exact_pair = tuple(sorted(pair))
            orientations = exact_pair_orientations.setdefault(
                exact_pair, {exact_pair: 0, exact_pair[::-1]: 0})
            orientations[pair] += 1
            if states == (1, 1):
                pair_kinds["on-aa"] += 1
                pair_kinds_by_period[round_row["round"]]["on-aa"] += 1
            elif states == (0, 0):
                pair_kinds["off-aa"] += 1
                pair_kinds_by_period[round_row["round"]]["off-aa"] += 1
            elif states[0] != states[1]:
                pair_kinds["cross"] += 1
                pair_kinds_by_period[round_row["round"]]["cross"] += 1
            else:
                raise CampaignError("unknown Williams pair")
            for cpu_name, treatment in (
                    ("cpu_a", pair[0]), ("cpu_b", pair[1])):
                counts[treatment][cpu_name] += 1
                counts[treatment]["total"] += 1
                block_counts[treatment] += 1
                cpu_sequences[cpu_name].append(treatment)
                cpu_period_treatments[
                    cpu_name][round_row["round"]][treatment] += 1
        if any(value != 2 for value in block_counts.values()):
            raise CampaignError("Williams block is treatment-imbalanced")
        for cpu, sequence in cpu_sequences.items():
            for transition in zip(sequence, sequence[1:]):
                carryovers[cpu][transition] = (
                    carryovers[cpu].get(transition, 0) + 1)
    for treatment in TREATMENTS:
        if counts[treatment] != {
                "cpu_a": 16, "cpu_b": 16,
                "total": OBSERVATIONS_PER_TREATMENT}:
            raise CampaignError("Williams CPU balance changed")
    if pair_kinds != {"on-aa": 16, "off-aa": 16, "cross": 32}:
        raise CampaignError("Williams pair balance changed")
    if any(row != {"on-aa": 4, "off-aa": 4, "cross": 8}
           for row in pair_kinds_by_period):
        raise CampaignError("Williams period balance changed")
    expected_exact_pairs = {
        tuple(sorted(pair))
        for pair in (
            ("hook-on-a", "hook-on-b"),
            ("hook-off-a", "hook-off-b"),
            ("hook-on-b", "hook-off-a"),
            ("hook-on-a", "hook-off-b"),
        )
    }
    if (set(exact_pair_orientations) != expected_exact_pairs or
            any(sorted(orientations.values()) != [8, 8]
                for orientations in exact_pair_orientations.values())):
        raise CampaignError("Williams exact-pair CPU orientation changed")
    if any(
            any(count != 4 for count in rows.values())
            for rows in cpu_row_counts.values()):
        raise CampaignError("Williams CPU row balance changed")
    if any(
            any(count != 4 for count in period.values())
            for cpu_periods in cpu_period_treatments.values()
            for period in cpu_periods):
        raise CampaignError("Williams CPU-period balance changed")
    directed_nonself = {
        (left, right)
        for left in TREATMENTS for right in TREATMENTS
        if left != right
    }
    if any(
            set(rows) != directed_nonself or
            any(count != 4 for count in rows.values())
            for rows in carryovers.values()):
        raise CampaignError("Williams directed carryover balance changed")
    warm_counts = {treatment: 0 for treatment in TREATMENTS}
    for block in range(WARMUP_BLOCKS):
        for round_row in williams_block(block):
            warm_counts[round_row["cpu_a_treatment"]] += 1
            warm_counts[round_row["cpu_b_treatment"]] += 1
    if any(value != WARMUPS_PER_TREATMENT
           for value in warm_counts.values()):
        raise CampaignError("warmup Williams blocks changed")


def frozen_timed_cases() -> List[Dict[str, Any]]:
    cases = []
    for K, block_bytes in PRIMARY_CELLS:
        for scope in SCOPES:
            cases.append({
                "case_id": "primary-k{}-b{}-{}".format(
                    K, block_bytes, scope),
                "kind": "primary",
                "K": K,
                "block_bytes": block_bytes,
                "schedule": PRIMARY_SCHEDULE,
                "loss_ppm": PRIMARY_LOSS_PPM,
                "scope": scope,
            })
    for anchor in TIMED_ANCHORS:
        for scope in SCOPES:
            cases.append({
                "case_id": "anchor-{}-k{}-b{}-{}".format(
                    anchor["schedule"], anchor["K"],
                    anchor["block_bytes"], scope),
                "kind": "anchor",
                "K": anchor["K"],
                "block_bytes": anchor["block_bytes"],
                "schedule": anchor["schedule"],
                "loss_ppm": PRIMARY_LOSS_PPM,
                "scope": scope,
            })
    if (len(cases) != 14 * 5 + 3 * 5 or
            len({row["case_id"] for row in cases}) != len(cases)):
        raise CampaignError("frozen timed-case expansion changed")
    return cases


def frozen_preflight_cases() -> List[Dict[str, Any]]:
    cases = []
    for K, block_bytes in PRIMARY_CELLS:
        for schedule in SCHEDULES:
            cases.append({
                "case_id": "semantic-{}-k{}-b{}".format(
                    schedule, K, block_bytes),
                "K": K,
                "block_bytes": block_bytes,
                "schedule": schedule,
                "loss_ppm": PRIMARY_LOSS_PPM,
                "scope": SEMANTIC_SCOPE,
            })
    if len(cases) != 14 * 6:
        raise CampaignError("frozen preflight expansion changed")
    return cases


def pre_results_contract(
        build_receipt_sha256: str,
        raw_root: int,
        loss_root: int,
        max_working_mib: int,
        cpu_a: int,
        cpu_b: int,
        helper_sha256: str,
        timeout_seconds: int = 1800,
        calibration_max_inner_reps: int = MAX_INNER_REPS,
        calibration_max_steps: int = MAX_CALIBRATION_STEPS,
        runner_python_sha256: Optional[str] = None) -> Dict[str, Any]:
    contract = {
        "schema": SCHEMA,
        "build_receipt_sha256": build_receipt_sha256,
        "campaign_helper_sha256": helper_sha256,
        "roots": {
            "construction_raw_root": raw_root,
            "loss_root": loss_root,
            "derivation": "passed-verbatim-to-every-cell-v1",
        },
        "cpus": {
            "cpu_a": cpu_a,
            "cpu_b": cpu_b,
            "contract": "validated-mutual-smt-siblings-v1",
        },
        "max_working_mib": max_working_mib,
        "primary": {
            "cells": [list(cell) for cell in PRIMARY_CELLS],
            "scopes": list(SCOPES),
            "schedule": PRIMARY_SCHEDULE,
            "loss_ppm": PRIMARY_LOSS_PPM,
        },
        "semantic_preflight": {
            "cells": [list(cell) for cell in PRIMARY_CELLS],
            "schedules": list(SCHEDULES),
            "scope": SEMANTIC_SCOPE,
            "treatments": list(TREATMENTS),
        },
        "timed_anchors": [
            {
                **anchor,
                "scopes": list(SCOPES),
                "loss_ppm": PRIMARY_LOSS_PPM,
            }
            for anchor in TIMED_ANCHORS
        ],
        "calibration": {
            "initial_inner_reps": 1,
            "shared_across_treatments": True,
            "minimum_elapsed_ns": TARGET_BATCH_NS,
            "maximum_inner_reps": calibration_max_inner_reps,
            "maximum_steps": calibration_max_steps,
            "scale_formula":
                "ceil(current*target/min_elapsed); at least current+1",
        },
        "warmup": {
            "blocks": WARMUP_BLOCKS,
            "invocations_per_treatment": WARMUPS_PER_TREATMENT,
            "native_warmup_reps": 0,
            "included_in_estimates": False,
        },
        "measurement": {
            "blocks": MEASURED_BLOCKS,
            "scheduled_invocations_per_treatment":
                OBSERVATIONS_PER_TREATMENT,
            "successful_timing_observations_per_treatment":
                OBSERVATIONS_PER_TREATMENT,
            "weak_root_timing_observations_per_treatment": 0,
            "weak_root_summary": None,
            "native_warmup_reps": 0,
            "williams_rows": [list(row) for row in WILLIAMS_ROWS],
            "block_rule":
                "start rows cycle 0,1,2,3 with the next row modulo 4; "
                "swap CPU rows on alternate four-block cycles",
        },
        "statistics": {
            "raw_cross_log_ratio":
                "ln(elapsed_ns_hook_on/elapsed_ns_hook_off)",
            "effect_block_log_ratio":
                "within each Williams block, arithmetic mean of its "
                "two raw cross log ratios",
            "effect_sampling_unit":
                "16 Williams-block inferential sampling units (df=15); "
                "32 within-block cross pairs are not treated as independent",
            "effect_components": [
                "hook-on-a/hook-off-b: one pair per block",
                "hook-on-b/hook-off-a: one pair per block",
            ],
            "aa_contrasts": [
                "hook-on-a/hook-on-b: one pair per block",
                "hook-off-a/hook-off-b: one pair per block",
            ],
            "ci_assumption":
                "Student-t CI assumes independence across the 16 "
                "Williams-block inferential sampling units",
            "mean": "math.fsum(log_ratios)/n",
            "ci":
                "mean +/- t_0.975_(n-1)*sqrt("
                "fsum((x-mean)^2)/(n-1))/sqrt(n)",
            "bootstrap": False,
            "aa_floor":
                "max(abs(expm1(ci_low)),abs(expm1(ci_high))) "
                "for hook-on A/A and hook-off A/A",
            "promotion_conclusion": None,
        },
        "unsupported_integrations": {
            "filler_stop_resume": "fail-closed-not-implemented",
            "thermal_window": "fail-closed-not-implemented",
        },
        "execution_limits": {
            "pair_timeout_seconds": timeout_seconds,
            "stdout_stderr_bytes_per_stream": MAX_STREAM_BYTES,
        },
        "binary_execution": {
            "contract": "sealed-authenticated-memfd-v1",
            "argv0": "/proc/self/fd/<execution-fd>",
            "source": "exact-authenticated-binary-bytes",
            "execution_mode": EXECUTION_MODE,
            "required_seal_mask": EXECUTION_SEALS,
            "required_seals": [
                "F_SEAL_SEAL", "F_SEAL_SHRINK",
                "F_SEAL_GROW", "F_SEAL_WRITE", "F_SEAL_EXEC"],
            "dynamic_runtime_closure":
                "ambient-host; bind and compare immediately before/after run",
        },
        "runner_python_sha256": runner_python_sha256,
    }
    return contract


def _parse_cpu_list(value: str) -> set:
    cpus = set()
    if not value or value.strip() != value:
        raise CampaignError("malformed CPU sibling list")
    for part in value.split(","):
        match = re.fullmatch(r"([0-9]+)(?:-([0-9]+))?", part)
        if match is None:
            raise CampaignError("malformed CPU sibling list")
        endpoints = (match.group(1), match.group(2) or match.group(1))
        if any(
                len(endpoint) > 7 or
                re.fullmatch(r"0|[1-9][0-9]*", endpoint) is None
                for endpoint in endpoints):
            raise CampaignError("malformed CPU sibling list")
        first, last = (int(endpoint) for endpoint in endpoints)
        if first > 1048575 or last > 1048575:
            raise CampaignError("CPU sibling exceeds protocol limit")
        if last < first or last - first > 4096:
            raise CampaignError("invalid CPU sibling range")
        cpus.update(range(first, last + 1))
    return cpus


def _parse_topology_uint(value: str, name: str) -> int:
    if re.fullmatch(r"0|[1-9][0-9]*", value) is None or len(value) > 10:
        raise CampaignError("{} is not canonical unsigned decimal".format(
            name))
    return _int(int(value), name, 0, UINT32_MAX)


def validate_smt_siblings(cpu_a: int, cpu_b: int,
                          topology_root: Path = Path(
                              "/sys/devices/system/cpu")) -> Dict[str, Any]:
    _int(cpu_a, "cpu-a", 0, 1048575)
    _int(cpu_b, "cpu-b", 0, 1048575)
    if cpu_a == cpu_b:
        raise CampaignError("two distinct CPUs are required")
    if not hasattr(os, "sched_getaffinity"):
        raise CampaignError("CPU affinity is unavailable")
    allowed = os.sched_getaffinity(0)
    if cpu_a not in allowed or cpu_b not in allowed:
        raise CampaignError("supplied CPUs are outside campaign affinity")

    rows = {}
    for cpu in (cpu_a, cpu_b):
        base = topology_root / "cpu{}".format(cpu) / "topology"
        values = {}
        for name in (
                "thread_siblings_list", "core_id", "physical_package_id"):
            path = base / name
            try:
                payload = path.read_bytes()
            except OSError as error:
                raise CampaignError("cannot read CPU topology: {}".format(
                    error))
            if len(payload) > 4096:
                raise CampaignError("CPU topology value is oversized")
            try:
                text = payload.decode("ascii").strip()
            except UnicodeDecodeError:
                raise CampaignError("CPU topology is not ASCII")
            if name == "thread_siblings_list":
                values[name] = text
            else:
                values[name] = _parse_topology_uint(
                    text, "CPU {} {}".format(cpu, name))
        rows[cpu] = values
    siblings_a = _parse_cpu_list(rows[cpu_a]["thread_siblings_list"])
    siblings_b = _parse_cpu_list(rows[cpu_b]["thread_siblings_list"])
    if (cpu_b not in siblings_a or cpu_a not in siblings_b or
            rows[cpu_a]["core_id"] != rows[cpu_b]["core_id"] or
            rows[cpu_a]["physical_package_id"] !=
            rows[cpu_b]["physical_package_id"]):
        raise CampaignError("supplied CPUs are not mutual SMT siblings")
    return {
        "cpu_a": cpu_a,
        "cpu_b": cpu_b,
        "cpu_a_thread_siblings": sorted(siblings_a),
        "cpu_b_thread_siblings": sorted(siblings_b),
        "core_id": rows[cpu_a]["core_id"],
        "physical_package_id": rows[cpu_a]["physical_package_id"],
    }


@dataclasses.dataclass(frozen=True)
class InvocationSpec:
    ordinal: int
    pair_ordinal: int
    phase: str
    case_id: str
    block: Optional[int]
    round_index: Optional[int]
    treatment: str
    cpu: int
    executable_fd: int
    expected: InvocationExpected
    command: Tuple[str, ...]


@dataclasses.dataclass(frozen=True)
class ProcessCapture:
    returncode: int
    stdout: bytes
    stderr: bytes
    launched_command: Tuple[str, ...]
    timed_out: bool = False


@dataclasses.dataclass
class InvocationOutcome:
    spec: InvocationSpec
    capture: ProcessCapture
    parsed: Optional[ParsedTimingStream]
    raw_ordinal: Optional[int] = None


def invocation_command(
        binary: BinaryBinding, expected: InvocationExpected) -> Tuple[str, ...]:
    return (
        "/proc/self/fd/{}".format(binary.execution_fd),
        "--K", str(expected.K),
        "--bb", str(expected.block_bytes),
        "--construction-seed", str(expected.construction_seed),
        "--loss-seed", str(expected.loss_seed),
        "--loss-ppm", str(expected.loss_ppm),
        "--schedule", expected.schedule,
        "--scope", expected.scope,
        "--warmup-reps", "0",
        "--inner-reps", str(expected.inner_reps),
        "--max-working-mib", str(expected.max_working_mib),
        "--context-sha256", expected.context_sha256,
    )


def _child_setup(cpu: int, inherited_signal_mask: set) -> None:
    os.sched_setaffinity(0, {cpu})
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    signal.pthread_sigmask(signal.SIG_SETMASK, inherited_signal_mask)


def _group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False


def _terminate_children(
        children: Sequence[subprocess.Popen],
        pgids: Sequence[int]) -> None:
    for pgid in pgids:
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + PROCESS_REAP_GRACE_SECONDS
    while time.monotonic() < deadline:
        for child in children:
            if child.poll() is not None:
                child.wait()
        if not any(_group_exists(pgid) for pgid in pgids):
            break
        if pgids:
            time.sleep(0.01)
    for pgid in pgids:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    for child in children:
        try:
            child.wait(timeout=PROCESS_REAP_GRACE_SECONDS)
        except subprocess.TimeoutExpired:
            raise CampaignError("child could not be reaped")


def _consume_pending_signals(signals_to_consume: set) -> List[int]:
    consumed = []
    while signals_to_consume:
        pending = set(signal.sigpending()).intersection(signals_to_consume)
        if not pending:
            break
        for signum in sorted(pending):
            received = signal.sigwait({signum})
            consumed.append(int(received))
    return consumed


def _block_managed_signals_safely() -> set:
    # CPython calls pthread_sigmask before it allocates the returned set and
    # checks pending handlers.  Snapshot first so an exception after the real
    # mask change cannot strand SIGINT/SIGTERM blocked.
    previous = signal.pthread_sigmask(signal.SIG_BLOCK, set())
    try:
        signal.pthread_sigmask(signal.SIG_BLOCK, set(MANAGED_SIGNALS))
    except BaseException as block_error:
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous)
        except BaseException as restore_error:
            raise CampaignError(
                "managed signal block failed and the prior mask could not "
                "be restored: {}: {}".format(
                    type(restore_error).__name__, restore_error)
            ) from block_error
        raise
    return previous


def run_simultaneous_pair(
        specs: Tuple[InvocationSpec, InvocationSpec],
        timeout_seconds: float,
        environment: Optional[Dict[str, str]] = None,
        stream_limit: int = MAX_STREAM_BYTES,
) -> Tuple[ProcessCapture, ProcessCapture]:
    if len(specs) != 2 or specs[0].cpu == specs[1].cpu:
        raise CampaignError("simultaneous pair needs two distinct CPUs")
    if timeout_seconds <= 0 or timeout_seconds > 86400:
        raise CampaignError("invalid invocation timeout")
    if stream_limit <= 0 or stream_limit > MAX_STREAM_BYTES:
        raise CampaignError("invalid stream limit")
    if environment is None:
        environment = {"LC_ALL": "C", "LANG": "C", "TZ": "UTC"}
    forbidden = (
        "WIREHAIR_V2_PEEL_DEGREES",
        "WIREHAIR_V2_STAIRCASE_DEGREES",
        "WIREHAIR_V2_STAIRCASE_ROW_DEGREES",
        "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE",
        "WIREHAIR_V2_BAND_TRACKING_X",
    )
    if any(name in environment for name in forbidden):
        raise CampaignError("timing environment contains equation override")
    if not all(hasattr(signal, name) for name in (
            "pthread_sigmask", "sigpending", "sigwait")):
        raise CampaignError("managed POSIX signal masking is unavailable")

    barrier_enabled = specs[0].expected.scope != SEMANTIC_SCOPE
    if (specs[1].expected.scope != SEMANTIC_SCOPE) != barrier_enabled:
        raise CampaignError("pair mixes barrier and non-barrier scopes")
    children = []
    pgids = []
    selector = None
    buffers: Dict[Tuple[int, str], bytearray] = {}
    ready_buffers: Dict[int, bytearray] = {}
    ready_eof = set()
    parent_pipe_fds = set()
    pending_child_pipe_fds = set()
    go_writers: Dict[int, int] = {}
    launched_commands: Dict[int, Tuple[str, ...]] = {}
    original_signal_mask = _block_managed_signals_safely()
    try:
        consumable_signals = (
            set(MANAGED_SIGNALS) - set(original_signal_mask))
        selector = selectors.DefaultSelector()
        for index, spec in enumerate(specs):
            child_fds: Tuple[int, ...] = ()
            launched_command = tuple(spec.command)
            if barrier_enabled:
                ready_read, ready_write = os.pipe2(os.O_CLOEXEC)
                parent_pipe_fds.add(ready_read)
                pending_child_pipe_fds.add(ready_write)
                go_read, go_write = os.pipe2(os.O_CLOEXEC)
                pending_child_pipe_fds.add(go_read)
                parent_pipe_fds.add(go_write)
                go_writers[index] = go_write
                child_fds = (ready_write, go_read)
                launched_command += (
                    "--ready-fd", str(ready_write),
                    "--go-fd", str(go_read),
                )
            child = subprocess.Popen(
                list(launched_command),
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd="/",
                env=environment,
                close_fds=True,
                pass_fds=(spec.executable_fd,) + child_fds,
                start_new_session=True,
                preexec_fn=lambda cpu=spec.cpu, mask=original_signal_mask:
                    _child_setup(cpu, mask),
            )
            children.append(child)
            pgids.append(child.pid)
            launched_commands[index] = launched_command
            for descriptor in child_fds:
                os.close(descriptor)
                pending_child_pipe_fds.discard(descriptor)
            assert child.stdout is not None and child.stderr is not None
            for name, stream in (
                    ("stdout", child.stdout), ("stderr", child.stderr)):
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ, (index, name))
                buffers[(index, name)] = bytearray()
            if barrier_enabled:
                selector.register(
                    ready_read, selectors.EVENT_READ, (index, "ready"))
                ready_buffers[index] = bytearray()
            received = _consume_pending_signals(consumable_signals)
            if received:
                raise CampaignInterrupted(
                    "campaign interrupted by signal {}".format(received[0]))
        deadline = time.monotonic() + timeout_seconds
        barrier_released = not barrier_enabled
        while selector.get_map() or any(
                child.poll() is None for child in children):
            received = _consume_pending_signals(consumable_signals)
            if received:
                raise CampaignInterrupted(
                    "campaign interrupted by signal {}".format(received[0]))
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("simultaneous invocation timed out")
            events = selector.select(min(0.05, remaining))
            for key, unused_mask in events:
                index, name = key.data
                try:
                    descriptor = (
                        key.fileobj if isinstance(key.fileobj, int)
                        else key.fileobj.fileno())
                    chunk = os.read(
                        descriptor, 2 if name == "ready" else 65536)
                except BlockingIOError:
                    continue
                if not chunk:
                    selector.unregister(key.fileobj)
                    if isinstance(key.fileobj, int):
                        os.close(key.fileobj)
                        parent_pipe_fds.discard(key.fileobj)
                        ready_eof.add(index)
                    else:
                        key.fileobj.close()
                    continue
                if name == "ready":
                    ready_buffers[index].extend(chunk)
                    if ready_buffers[index] != bytearray(b"R"):
                        raise CampaignError(
                            "child emitted invalid ready barrier bytes")
                    continue
                buffer = buffers[(index, name)]
                buffer.extend(chunk)
                if len(buffer) > stream_limit:
                    raise CampaignError(
                        "{} exceeded stream limit".format(name))
            if (barrier_enabled and not barrier_released and
                    any(child.poll() is not None for child in children)):
                raise CampaignError("child exited before start barrier")
            if (barrier_enabled and not barrier_released and
                    ready_eof == {0, 1}):
                if any(ready_buffers[index] != bytearray(b"R")
                       for index in (0, 1)):
                    raise CampaignError("child closed ready pipe without R")
                for index in (0, 1):
                    descriptor = go_writers[index]
                    if os.write(descriptor, b"G") != 1:
                        raise CampaignError("short go-barrier write")
                    os.close(descriptor)
                    parent_pipe_fds.discard(descriptor)
                barrier_released = True
        captures = []
        for index, child in enumerate(children):
            returncode = child.wait()
            captures.append(ProcessCapture(
                returncode=returncode,
                stdout=bytes(buffers[(index, "stdout")]),
                stderr=bytes(buffers[(index, "stderr")]),
                launched_command=launched_commands[index],
            ))
        leaked_groups = [pgid for pgid in pgids if _group_exists(pgid)]
        if leaked_groups:
            _terminate_children(children, leaked_groups)
            raise CampaignError("timing child left descendant processes")
        received = _consume_pending_signals(consumable_signals)
        if received:
            raise CampaignInterrupted(
                "campaign interrupted by signal {}".format(received[0]))
        return captures[0], captures[1]
    except BaseException as error:
        for descriptor in list(parent_pipe_fds):
            try:
                os.close(descriptor)
            except OSError:
                pass
            parent_pipe_fds.discard(descriptor)
        for descriptor in list(pending_child_pipe_fds):
            try:
                os.close(descriptor)
            except OSError:
                pass
            pending_child_pipe_fds.discard(descriptor)
        _terminate_children(children, pgids)
        if isinstance(error, TimeoutError):
            raise CampaignError(str(error))
        raise
    finally:
        try:
            if selector is not None:
                selector.close()
            for descriptor in list(parent_pipe_fds):
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            for descriptor in list(pending_child_pipe_fds):
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            for child in children:
                for stream in (child.stdout, child.stderr):
                    if stream is not None and not stream.closed:
                        stream.close()
        finally:
            signal.pthread_sigmask(
                signal.SIG_SETMASK, original_signal_mask)


class AppendOnlyRawLog:
    def __init__(self, path: Path):
        # Preserve the final path component so O_EXCL rejects symlinks instead
        # of resolving a dangling link to an unintended destination.
        self.path = Path(os.path.abspath(os.fspath(path)))
        if not self.path.parent.is_dir():
            raise CampaignError("raw-log parent directory does not exist")
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_APPEND | os.O_CLOEXEC
        try:
            self.descriptor = os.open(str(self.path), flags, 0o644)
        except OSError as error:
            raise CampaignError("cannot create raw log: {}".format(error))
        try:
            _fsync_directory(self.path.parent, "raw-log parent")
        except BaseException:
            try:
                os.close(self.descriptor)
            except OSError:
                pass
            try:
                os.unlink(self.path)
            except OSError:
                pass
            raise
        self.sha256 = hashlib.sha256()
        self.bytes = 0
        self.records = 0
        self.closed = False

    def append_many(self, rows: Sequence[Dict[str, Any]]) -> List[int]:
        if self.closed:
            raise CampaignError("raw log is closed")
        payloads = []
        ordinals = []
        for index, row in enumerate(rows):
            expected_ordinal = self.records + index
            if row.get("raw_ordinal") != expected_ordinal:
                raise CampaignError("raw record ordinal is not append-only")
            payload = canonical_json_bytes(row)
            if len(payload) > MAX_RAW_RECORD_BYTES:
                raise CampaignError("raw invocation record is oversized")
            payloads.append(payload)
            ordinals.append(expected_ordinal)
        combined = b"".join(payloads)
        view = memoryview(combined)
        try:
            while view:
                written = os.write(self.descriptor, view)
                if written <= 0:
                    raise CampaignError("short raw-log write")
                view = view[written:]
            os.fsync(self.descriptor)
        except OSError as error:
            raise CampaignError("raw-log append failed: {}".format(error))
        self.sha256.update(combined)
        self.bytes += len(combined)
        self.records += len(rows)
        return ordinals

    def close(self) -> Dict[str, Any]:
        if not self.closed:
            close_error = None
            descriptor = self.descriptor
            try:
                os.fsync(descriptor)
            except OSError as error:
                close_error = error
            try:
                os.close(descriptor)
            except OSError as error:
                if close_error is None:
                    close_error = error
            finally:
                self.descriptor = -1
                self.closed = True
            if close_error is not None:
                raise CampaignError(
                    "raw-log close failed: {}".format(close_error))
            _fsync_directory(self.path.parent, "raw-log parent")
        payload, binding = _stable_file_bytes(
            self.path, max(self.bytes, 1), "raw invocation log")
        digest = hashlib.sha256(payload).hexdigest()
        if (len(payload) != self.bytes or digest != self.sha256.hexdigest()):
            raise CampaignError("raw invocation log changed")
        return {
            "path": str(self.path),
            "bytes": self.bytes,
            "records": self.records,
            "sha256": digest,
            "device": binding["device"],
            "inode": binding["inode"],
        }


def _raw_record(
        raw_ordinal: int,
        outcome: InvocationOutcome,
        validation_error: Optional[str] = None) -> Dict[str, Any]:
    capture = outcome.capture
    parsed = outcome.parsed
    return {
        "schema": RAW_SCHEMA,
        "raw_ordinal": raw_ordinal,
        "invocation_ordinal": outcome.spec.ordinal,
        "pair_ordinal": outcome.spec.pair_ordinal,
        "phase": outcome.spec.phase,
        "case_id": outcome.spec.case_id,
        "block": outcome.spec.block,
        "round": outcome.spec.round_index,
        "treatment": outcome.spec.treatment,
        "cpu": outcome.spec.cpu,
        "command": list(capture.launched_command),
        "returncode": capture.returncode,
        "timed_out": capture.timed_out,
        "stdout": {
            "bytes": len(capture.stdout),
            "sha256": hashlib.sha256(capture.stdout).hexdigest(),
            "base64": base64.b64encode(capture.stdout).decode("ascii"),
        },
        "stderr": {
            "bytes": len(capture.stderr),
            "sha256": hashlib.sha256(capture.stderr).hexdigest(),
            "base64": base64.b64encode(capture.stderr).decode("ascii"),
        },
        "validation": {
            "ok": validation_error is None,
            "error": validation_error,
            "status": parsed.status if parsed is not None else None,
            "semantic_sha256":
                parsed.semantic["semantic_sha256"]
                if parsed is not None else None,
            "done_stream_sha256":
                parsed.done["stream_sha256"]
                if parsed is not None else None,
        },
    }


def _mean_ci(samples: Sequence[float]) -> Dict[str, Any]:
    if len(samples) < 2:
        raise CampaignError("at least two paired samples are required")
    if len(samples) - 1 >= len(T_CRITICAL_975):
        raise CampaignError("Student-t table does not cover sample count")
    if any(not math.isfinite(value) for value in samples):
        raise CampaignError("non-finite statistical sample")
    count = len(samples)
    mean = math.fsum(samples) / count
    sum_squares = math.fsum((value - mean) ** 2 for value in samples)
    sample_sd = math.sqrt(sum_squares / (count - 1))
    standard_error = sample_sd / math.sqrt(count)
    critical = T_CRITICAL_975[count - 1]
    assert critical is not None
    half_width = critical * standard_error
    low = mean - half_width
    high = mean + half_width
    return {
        "n": count,
        "degrees_of_freedom": count - 1,
        "mean_log_ratio": mean,
        "sample_sd_log_ratio": sample_sd,
        "standard_error_log_ratio": standard_error,
        "t_critical_975": critical,
        "ci95_low_log_ratio": low,
        "ci95_high_log_ratio": high,
        "geometric_ratio": math.exp(mean),
        "geometric_percent": math.expm1(mean) * 100.0,
        "ci95_low_ratio": math.exp(low),
        "ci95_high_ratio": math.exp(high),
    }


def summarize_measurements(
        pairs: Sequence[Tuple[InvocationOutcome, InvocationOutcome]],
) -> Dict[str, Any]:
    cross_by_block = {block: [] for block in range(MEASURED_BLOCKS)}
    on_aa_by_block = {block: [] for block in range(MEASURED_BLOCKS)}
    off_aa_by_block = {block: [] for block in range(MEASURED_BLOCKS)}
    by_treatment = {treatment: [] for treatment in TREATMENTS}
    seen_slots = set()
    seen_pair_ordinals = set()
    cross_components = {
        ("hook-on-a", "hook-off-b"): {
            "pairs": 0, "numerator_cpu_a": 0, "numerator_cpu_b": 0},
        ("hook-on-b", "hook-off-a"): {
            "pairs": 0, "numerator_cpu_a": 0, "numerator_cpu_b": 0},
    }
    aa_numerator_roles = {
        "hook-on-a_over_hook-on-b": {
            "pairs": 0, "numerator_cpu_a": 0, "numerator_cpu_b": 0},
        "hook-off-a_over_hook-off-b": {
            "pairs": 0, "numerator_cpu_a": 0, "numerator_cpu_b": 0},
    }
    for left, right in pairs:
        if left.parsed is None or right.parsed is None:
            raise CampaignError("unparsed measured outcome")
        if left.parsed.timing is None or right.parsed.timing is None:
            raise CampaignError("successful measurement lacks timing")
        pair = {left.spec.treatment: left, right.spec.treatment: right}
        if len(pair) != 2:
            raise CampaignError("measured pair repeats a treatment")
        block = left.spec.block
        round_index = left.spec.round_index
        if (type(block) is not int or block not in cross_by_block or
                type(round_index) is not int or
                round_index < 0 or round_index >= 4 or
                right.spec.block != block or
                right.spec.round_index != round_index or
                left.spec.phase != "measured" or
                right.spec.phase != "measured" or
                left.spec.pair_ordinal != right.spec.pair_ordinal or
                left.spec.ordinal + 1 != right.spec.ordinal or
                left.spec.cpu == right.spec.cpu):
            raise CampaignError("measured block/round metadata changed")
        slot = (block, round_index)
        if slot in seen_slots:
            raise CampaignError("duplicate measured block/round slot")
        if left.spec.pair_ordinal in seen_pair_ordinals:
            raise CampaignError("duplicate measured pair ordinal")
        seen_slots.add(slot)
        seen_pair_ordinals.add(left.spec.pair_ordinal)
        designed = williams_block(block)[round_index]
        if (
                left.spec.treatment != designed["cpu_a_treatment"] or
                right.spec.treatment != designed["cpu_b_treatment"]):
            raise CampaignError("measured Williams treatment order changed")
        for outcome in (left, right):
            by_treatment[outcome.spec.treatment].append(outcome)
        states = {
            TREATMENT_STATE[left.spec.treatment],
            TREATMENT_STATE[right.spec.treatment],
        }
        if states == {0, 1}:
            on = left if TREATMENT_STATE[left.spec.treatment] else right
            off = right if on is left else left
            component = cross_components.get(
                (on.spec.treatment, off.spec.treatment))
            if component is None:
                raise CampaignError("unknown measured cross component")
            component["pairs"] += 1
            component[
                "numerator_cpu_a" if on is left
                else "numerator_cpu_b"] += 1
            cross_by_block[block].append(math.log(
                on.parsed.timing["elapsed_ns"] /
                off.parsed.timing["elapsed_ns"]))
        elif states == {1}:
            numerator = pair.get("hook-on-a")
            denominator = pair.get("hook-on-b")
            if numerator is None or denominator is None:
                raise CampaignError("hook-on A/A pair is malformed")
            role = aa_numerator_roles["hook-on-a_over_hook-on-b"]
            role["pairs"] += 1
            role[
                "numerator_cpu_a" if numerator is left
                else "numerator_cpu_b"] += 1
            on_aa_by_block[block].append(math.log(
                numerator.parsed.timing["elapsed_ns"] /
                denominator.parsed.timing["elapsed_ns"]))
        elif states == {0}:
            numerator = pair.get("hook-off-a")
            denominator = pair.get("hook-off-b")
            if numerator is None or denominator is None:
                raise CampaignError("hook-off A/A pair is malformed")
            role = aa_numerator_roles["hook-off-a_over_hook-off-b"]
            role["pairs"] += 1
            role[
                "numerator_cpu_a" if numerator is left
                else "numerator_cpu_b"] += 1
            off_aa_by_block[block].append(math.log(
                numerator.parsed.timing["elapsed_ns"] /
                denominator.parsed.timing["elapsed_ns"]))
        else:
            raise CampaignError("unknown measured pair state")
    if (
            seen_slots != {
                (block, round_index)
                for block in range(MEASURED_BLOCKS)
                for round_index in range(4)} or
            any(len(cross_by_block[block]) != 2
                for block in range(MEASURED_BLOCKS)) or
            any(len(on_aa_by_block[block]) != 1
                for block in range(MEASURED_BLOCKS)) or
            any(len(off_aa_by_block[block]) != 1
                for block in range(MEASURED_BLOCKS)) or
            any(len(rows) != OBSERVATIONS_PER_TREATMENT
                for rows in by_treatment.values())):
        raise CampaignError("measured Williams sample counts changed")
    expected_role_counts = {
        "pairs": MEASURED_BLOCKS,
        "numerator_cpu_a": MEASURED_BLOCKS // 2,
        "numerator_cpu_b": MEASURED_BLOCKS // 2,
    }
    if (
            any(row != expected_role_counts
                for row in cross_components.values()) or
            any(row != expected_role_counts
                for row in aa_numerator_roles.values())):
        raise CampaignError("measured contrast CPU-role balance changed")
    cross_block_effects = [
        math.fsum(cross_by_block[block]) / 2.0
        for block in range(MEASURED_BLOCKS)
    ]
    on_aa = [
        on_aa_by_block[block][0]
        for block in range(MEASURED_BLOCKS)
    ]
    off_aa = [
        off_aa_by_block[block][0]
        for block in range(MEASURED_BLOCKS)
    ]
    effect = _mean_ci(cross_block_effects)
    effect["raw_cross_pairs"] = 2 * MEASURED_BLOCKS
    effect["effect_blocks"] = MEASURED_BLOCKS
    effect["contrast_components"] = [
        {
            "numerator": numerator,
            "denominator": denominator,
            **cross_components[(numerator, denominator)],
        }
        for numerator, denominator in (
            ("hook-on-a", "hook-off-b"),
            ("hook-on-b", "hook-off-a"),
        )
    ]
    hook_on_aa = _mean_ci(on_aa)
    hook_off_aa = _mean_ci(off_aa)
    hook_on_aa["contrast"] = {
        "numerator": "hook-on-a",
        "denominator": "hook-on-b",
        **aa_numerator_roles["hook-on-a_over_hook-on-b"],
    }
    hook_off_aa["contrast"] = {
        "numerator": "hook-off-a",
        "denominator": "hook-off-b",
        **aa_numerator_roles["hook-off-a_over_hook-off-b"],
    }
    for summary in (hook_on_aa, hook_off_aa):
        summary["conservative_floor_fraction"] = max(
            abs(math.expm1(summary["ci95_low_log_ratio"])),
            abs(math.expm1(summary["ci95_high_log_ratio"])),
        )
    metrics = (
        "elapsed_ns", "minor_faults", "major_faults",
        "voluntary_context_switches", "involuntary_context_switches",
    )
    resources = {}
    for treatment in TREATMENTS:
        rows = by_treatment[treatment]
        resources[treatment] = {}
        for metric in metrics:
            values = [row.parsed.timing[metric] for row in rows]
            total = sum(values)
            resources[treatment][metric] = {
                "sum": total,
                "mean": total / len(values),
                "min": min(values),
                "max": max(values),
            }
    return {
        "hook_on_over_hook_off": effect,
        "hook_on_a_over_b": hook_on_aa,
        "hook_off_a_over_b": hook_off_aa,
        "aa_floor_fraction": max(
            hook_on_aa["conservative_floor_fraction"],
            hook_off_aa["conservative_floor_fraction"],
        ),
        "resource_usage": resources,
    }


def _binary_stat_matches(binding: BinaryBinding) -> bool:
    try:
        current = os.stat(binding.path, follow_symlinks=False)
    except OSError:
        return False
    return (
        stat.S_ISREG(current.st_mode) and
        current.st_dev == binding.device and
        current.st_ino == binding.inode and
        current.st_size == binding.bytes and
        stat.S_IMODE(current.st_mode) == binding.mode
    )


class CampaignRunner:
    def __init__(
            self,
            bindings: Dict[str, BinaryBinding],
            raw_log: AppendOnlyRawLog,
            context_sha256: str,
            raw_root: int,
            loss_root: int,
            max_working_mib: int,
            cpu_a: int,
            cpu_b: int,
            timeout_seconds: float,
            calibration_max_inner_reps: int,
            calibration_max_steps: int,
            process_environment: Optional[Dict[str, str]] = None):
        self.bindings = bindings
        self.raw_log = raw_log
        self.context_sha256 = context_sha256
        self.raw_root = raw_root
        self.loss_root = loss_root
        self.max_working_mib = max_working_mib
        self.cpu_a = cpu_a
        self.cpu_b = cpu_b
        self.timeout_seconds = timeout_seconds
        self.calibration_max_inner_reps = calibration_max_inner_reps
        self.calibration_max_steps = calibration_max_steps
        self.process_environment = process_environment
        self.next_invocation = 0
        self.next_pair = 0
        self.execution_order: List[Dict[str, Any]] = []
        self.semantic_registry: Dict[Tuple[int, int, str], Dict[str, Any]] = {}

    def _spec(
            self, treatment: str, cpu: int, phase: str,
            case: Dict[str, Any], inner_reps: int,
            block: Optional[int], round_index: Optional[int],
            pair_ordinal: int) -> InvocationSpec:
        expected = InvocationExpected(
            treatment=treatment,
            K=case["K"],
            block_bytes=case["block_bytes"],
            construction_seed=self.raw_root,
            loss_seed=self.loss_root,
            loss_ppm=case["loss_ppm"],
            schedule=case["schedule"],
            scope=case["scope"],
            inner_reps=inner_reps,
            max_working_mib=self.max_working_mib,
            context_sha256=self.context_sha256,
            cpu=cpu,
        )
        command = invocation_command(self.bindings[treatment], expected)
        spec = InvocationSpec(
            ordinal=self.next_invocation,
            pair_ordinal=pair_ordinal,
            phase=phase,
            case_id=case["case_id"],
            block=block,
            round_index=round_index,
            treatment=treatment,
            cpu=cpu,
            executable_fd=self.bindings[treatment].execution_fd,
            expected=expected,
            command=command,
        )
        self.next_invocation += 1
        self.execution_order.append({
            "invocation_ordinal": spec.ordinal,
            "pair_ordinal": pair_ordinal,
            "phase": phase,
            "case_id": case["case_id"],
            "block": block,
            "round": round_index,
            "treatment": treatment,
            "cpu": cpu,
            "inner_reps": inner_reps,
        })
        return spec

    def execute_pair(
            self, cpu_a_treatment: str, cpu_b_treatment: str,
            phase: str, case: Dict[str, Any], inner_reps: int,
            block: Optional[int] = None,
            round_index: Optional[int] = None,
    ) -> Tuple[InvocationOutcome, InvocationOutcome]:
        pair_ordinal = self.next_pair
        self.next_pair += 1
        specs = (
            self._spec(
                cpu_a_treatment, self.cpu_a, phase, case, inner_reps,
                block, round_index, pair_ordinal),
            self._spec(
                cpu_b_treatment, self.cpu_b, phase, case, inner_reps,
                block, round_index, pair_ordinal),
        )
        for spec in specs:
            if not _binary_stat_matches(self.bindings[spec.treatment]):
                raise CampaignError("binary changed before invocation")
        captures = run_simultaneous_pair(
            specs, self.timeout_seconds,
            environment=self.process_environment)
        outcomes = []
        validation_errors: List[Optional[str]] = []
        for spec, capture in zip(specs, captures):
            parsed = None
            validation_error = None
            try:
                if capture.returncode != 0:
                    raise CampaignError(
                        "timing child exited {}".format(capture.returncode))
                if capture.stderr:
                    raise CampaignError("successful timing child wrote stderr")
                parsed = parse_timing_stream(capture.stdout, spec.expected)
            except CampaignError as error:
                validation_error = str(error)
            outcome = InvocationOutcome(
                spec=spec, capture=capture, parsed=parsed)
            outcomes.append(outcome)
            validation_errors.append(validation_error)
        left, right = outcomes
        pair_error = None
        if any(error is not None for error in validation_errors):
            for index, error in enumerate(validation_errors):
                if error is None:
                    validation_errors[index] = (
                        "simultaneous peer failed stream validation")
        else:
            assert left.parsed is not None and right.parsed is not None
            if left.parsed.semantic != right.parsed.semantic:
                pair_error = "simultaneous pair has semantic drift"
            else:
                key = (
                    case["K"], case["block_bytes"], case["schedule"])
                previous = self.semantic_registry.get(key)
                if previous is not None and previous != left.parsed.semantic:
                    pair_error = (
                        "semantic result changed across invocations")
                elif previous is None:
                    self.semantic_registry[key] = left.parsed.semantic
            if pair_error is not None:
                validation_errors = [pair_error, pair_error]
        raw_start = self.raw_log.records
        raw_ordinals = self.raw_log.append_many([
            _raw_record(
                raw_start + index, outcome, validation_errors[index])
            for index, outcome in enumerate(outcomes)
        ])
        for outcome, raw_ordinal in zip(outcomes, raw_ordinals):
            outcome.raw_ordinal = raw_ordinal
        if any(error is not None for error in validation_errors):
            errors = [
                "{}: {}".format(outcome.spec.treatment, error)
                for outcome, error in zip(outcomes, validation_errors)
            ]
            raise CampaignError("; ".join(errors))
        return left, right

    def run_semantic_preflight(self) -> List[Dict[str, Any]]:
        summaries = []
        for case_index, case in enumerate(frozen_preflight_cases()):
            outcomes = []
            pairs = (
                ("hook-on-a", "hook-off-a"),
                ("hook-on-b", "hook-off-b"),
            )
            raw_ordinals = []
            for pair_index, pair in enumerate(pairs):
                if (case_index + pair_index) % 2:
                    pair = (pair[1], pair[0])
                result = self.execute_pair(
                    pair[0], pair[1], "semantic-preflight", case, 1,
                    block=None, round_index=pair_index)
                outcomes.extend(result)
                raw_ordinals.extend(
                    outcome.raw_ordinal for outcome in result)
            statuses = {
                outcome.parsed.status for outcome in outcomes
                if outcome.parsed is not None}
            semantics = {
                canonical_json_bytes(outcome.parsed.semantic)
                for outcome in outcomes if outcome.parsed is not None}
            if len(statuses) != 1 or len(semantics) != 1:
                raise CampaignError("four-build semantic preflight disagrees")
            semantic = outcomes[0].parsed.semantic
            summaries.append({
                **case,
                "status": outcomes[0].parsed.status,
                "semantic_sha256": semantic["semantic_sha256"],
                "raw_ordinals": raw_ordinals,
            })
        return summaries

    def calibrate(self, case: Dict[str, Any]) -> Dict[str, Any]:
        inner_reps = 1
        iterations = []
        for step in range(self.calibration_max_steps):
            on_pair = ("hook-on-a", "hook-on-b")
            off_pair = ("hook-off-a", "hook-off-b")
            if step % 2:
                on_pair = (on_pair[1], on_pair[0])
                off_pair = (off_pair[1], off_pair[0])
            outcomes = []
            outcomes.extend(self.execute_pair(
                on_pair[0], on_pair[1], "calibration", case, inner_reps,
                block=step, round_index=0))
            outcomes.extend(self.execute_pair(
                off_pair[1], off_pair[0], "calibration", case, inner_reps,
                block=step, round_index=1))
            statuses = {outcome.parsed.status for outcome in outcomes}
            raw_ordinals = [outcome.raw_ordinal for outcome in outcomes]
            if len(statuses) != 1:
                raise CampaignError("calibration weak-root disagreement")
            if statuses == {"weak-root"}:
                iterations.append({
                    "step": step,
                    "inner_reps": inner_reps,
                    "status": "weak-root",
                    "raw_ordinals": raw_ordinals,
                })
                return {
                    "status": "weak-root",
                    "inner_reps": inner_reps,
                    "iterations": iterations,
                }
            elapsed = {
                outcome.spec.treatment:
                    outcome.parsed.timing["elapsed_ns"]
                for outcome in outcomes
            }
            minimum = min(elapsed.values())
            iterations.append({
                "step": step,
                "inner_reps": inner_reps,
                "status": "success",
                "elapsed_ns": elapsed,
                "minimum_elapsed_ns": minimum,
                "raw_ordinals": raw_ordinals,
            })
            if minimum >= TARGET_BATCH_NS:
                return {
                    "status": "success",
                    "inner_reps": inner_reps,
                    "iterations": iterations,
                }
            numerator = inner_reps * TARGET_BATCH_NS
            next_inner = max(
                inner_reps + 1,
                (numerator + minimum - 1) // minimum)
            if next_inner > self.calibration_max_inner_reps:
                raise CampaignError(
                    "calibration needs inner_reps above cap")
            inner_reps = next_inner
        raise CampaignError("calibration did not reach 100 ms before step cap")

    def _run_blocks(
            self, phase: str, case: Dict[str, Any], inner_reps: int,
            block_count: int, expected_status: str,
    ) -> Tuple[
            List[Tuple[InvocationOutcome, InvocationOutcome]],
            List[Dict[str, Any]]]:
        pairs = []
        receipts = []
        counts = {treatment: 0 for treatment in TREATMENTS}
        cpu_counts = {
            treatment: {self.cpu_a: 0, self.cpu_b: 0}
            for treatment in TREATMENTS}
        for block in range(block_count):
            for round_row in williams_block(block):
                pair = self.execute_pair(
                    round_row["cpu_a_treatment"],
                    round_row["cpu_b_treatment"],
                    phase, case, inner_reps,
                    block=block, round_index=round_row["round"])
                if any(outcome.parsed.status != expected_status
                       for outcome in pair):
                    raise CampaignError(
                        "{} status changed after calibration".format(phase))
                for outcome in pair:
                    counts[outcome.spec.treatment] += 1
                    cpu_counts[outcome.spec.treatment][outcome.spec.cpu] += 1
                pairs.append(pair)
                receipts.append({
                    "block": block,
                    "round": round_row["round"],
                    "cpu_a_treatment": pair[0].spec.treatment,
                    "cpu_b_treatment": pair[1].spec.treatment,
                    "raw_ordinals": [
                        pair[0].raw_ordinal, pair[1].raw_ordinal],
                })
        expected_count = block_count * 2
        if any(value != expected_count for value in counts.values()):
            raise CampaignError("{} treatment count is imbalanced".format(
                phase))
        if any(
                cpu_counts[treatment][self.cpu_a] != block_count or
                cpu_counts[treatment][self.cpu_b] != block_count
                for treatment in TREATMENTS):
            raise CampaignError("{} CPU assignment is imbalanced".format(
                phase))
        return pairs, receipts

    def run_timed_case(self, case: Dict[str, Any]) -> Dict[str, Any]:
        calibration = self.calibrate(case)
        status = calibration["status"]
        inner_reps = calibration["inner_reps"]
        unused_warmup_pairs, warmup_receipts = self._run_blocks(
            "warmup", case, inner_reps, WARMUP_BLOCKS, status)
        measured_pairs, measured_receipts = self._run_blocks(
            "measured", case, inner_reps, MEASURED_BLOCKS, status)
        if status == "success":
            summary = summarize_measurements(measured_pairs)
        else:
            if any(outcome.parsed.timing is not None
                   for pair in measured_pairs for outcome in pair):
                raise CampaignError("weak-root measurement emitted timing")
            summary = None
        return {
            **case,
            "status": status,
            "semantic_sha256":
                measured_pairs[0][0].parsed.semantic["semantic_sha256"],
            "calibration": calibration,
            "warmup_pairs": warmup_receipts,
            "measured_pairs": measured_receipts,
            "summary": summary,
            "promotion_conclusion": None,
        }


def _reverify_binaries(bindings: Dict[str, BinaryBinding]) -> None:
    for treatment in TREATMENTS:
        binding = bindings[treatment]
        unused_payload, actual = _stable_file_bytes(
            Path(binding.path), MAX_BINARY_BYTES,
            "{} binary final verification".format(treatment))
        if (actual["bytes"] != binding.bytes or
                actual["sha256"] != binding.sha256 or
                actual["device"] != binding.device or
                actual["inode"] != binding.inode or
                actual["mode"] != binding.mode):
            raise CampaignError("binary changed during campaign: {}".format(
                treatment))


def _publish_no_replace(path: Path, payload: bytes) -> Dict[str, Any]:
    # Keep the final component unresolved so a dangling symlink cannot redirect
    # publication.  The output is linked from an anonymous O_TMPFILE while both
    # that file and its parent directory remain open.  No mutable temporary
    # pathname participates in the commit.
    path = Path(os.path.abspath(os.fspath(path)))
    if not path.name or not path.parent.is_dir():
        raise CampaignError("output parent directory does not exist")
    if not hasattr(os, "O_TMPFILE"):
        raise CampaignError("anonymous O_TMPFILE publication is unavailable")

    directory_flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_DIRECTORY"):
        directory_flags |= os.O_DIRECTORY
    anonymous_flags = os.O_WRONLY | os.O_CLOEXEC | os.O_TMPFILE
    directory = -1
    descriptor = -1
    committed = False
    post_commit_warnings = []
    post_commit_interruptions = []

    try:
        directory = os.open(str(path.parent), directory_flags)
        directory_stat = os.fstat(directory)
        if not stat.S_ISDIR(directory_stat.st_mode):
            raise CampaignError("output parent is not a directory")
        try:
            os.stat(path.name, dir_fd=directory, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise CampaignError("output already exists")

        descriptor = os.open(
            ".", anonymous_flags, 0o644, dir_fd=directory)
        os.fchmod(descriptor, 0o644)
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("short receipt write")
            view = view[written:]
        os.fsync(descriptor)
        anonymous_stat = os.fstat(descriptor)

        previous_mask = _block_managed_signals_safely()
        try:
            link_error = None
            try:
                os.link(
                    "/proc/self/fd/{}".format(descriptor),
                    path.name,
                    dst_dir_fd=directory,
                    follow_symlinks=True,
                )
            except BaseException as error:
                link_error = error
            else:
                # A successful linkat is the no-replace commit point.  Later
                # validation is useful evidence but cannot turn this back into
                # a precommit failure.
                committed = True

            verification_error = None
            try:
                published_stat = os.stat(
                    path.name, dir_fd=directory, follow_symlinks=False)
                published_matches = (
                    stat.S_ISREG(published_stat.st_mode) and
                    anonymous_stat.st_dev == published_stat.st_dev and
                    anonymous_stat.st_ino == published_stat.st_ino)
            except BaseException as error:
                verification_error = error
                published_matches = False

            if not committed:
                if verification_error is not None:
                    raise CampaignError(
                        "receipt publication outcome is ambiguous after "
                        "{}: {}; verification failed with {}: {}".format(
                            type(link_error).__name__, link_error,
                            type(verification_error).__name__,
                            verification_error)
                    ) from link_error
                if not published_matches:
                    raise link_error
                committed = True
                post_commit_warnings.append(
                    "exception deferred after receipt link: {}: {}".format(
                        type(link_error).__name__, link_error))
            elif verification_error is not None:
                message = (
                    "receipt link verification failed after commit: {}: "
                    "{}".format(
                        type(verification_error).__name__,
                        verification_error))
                post_commit_warnings.append(message)
                if isinstance(verification_error, CampaignInterrupted):
                    post_commit_interruptions.append(message)
            elif not published_matches:
                post_commit_warnings.append(
                    "published receipt path changed after verified link")

            # Durability of the directory reached by the held descriptor does
            # not depend on whether its reported pathname has been renamed.
            try:
                os.fsync(directory)
            except BaseException as error:
                message = (
                    "directory sync failed after receipt link: {}: {}".format(
                        type(error).__name__, error))
                post_commit_warnings.append(message)
                if isinstance(error, CampaignInterrupted):
                    post_commit_interruptions.append(message)
            try:
                current_parent = os.stat(path.parent)
                if (current_parent.st_dev != directory_stat.st_dev or
                        current_parent.st_ino != directory_stat.st_ino):
                    raise CampaignError(
                        "reported output parent changed during publication")
            except BaseException as error:
                message = (
                    "reported output parent validation failed after receipt "
                    "link: {}: {}".format(type(error).__name__, error))
                post_commit_warnings.append(message)
                if isinstance(error, CampaignInterrupted):
                    post_commit_interruptions.append(message)
        finally:
            restore_preserves_primary = sys.exc_info()[0] is not None
            try:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
            except CampaignInterrupted as error:
                if committed:
                    message = (
                        "campaign interrupted after receipt link: {}".format(
                            error))
                    post_commit_warnings.append(message)
                    post_commit_interruptions.append(message)
                else:
                    raise
            except BaseException as error:
                if committed:
                    post_commit_warnings.append(
                        "signal-mask restoration failed after receipt link: "
                        "{}: {}".format(type(error).__name__, error))
                elif not restore_preserves_primary:
                    raise
    except OSError as error:
        if committed:
            post_commit_warnings.append(
                "post-commit publication error: {}: {}".format(
                    type(error).__name__, error))
        else:
            raise CampaignError(
                "could not publish receipt: {}".format(error))
    except BaseException as error:
        if committed:
            message = "post-commit publication error: {}: {}".format(
                type(error).__name__, error)
            post_commit_warnings.append(message)
            if isinstance(error, CampaignInterrupted):
                post_commit_interruptions.append(message)
        else:
            raise
    finally:
        close_preserves_primary = sys.exc_info()[0] is not None
        close_errors = []
        for descriptor_to_close, description in (
                (descriptor, "anonymous receipt"),
                (directory, "output parent")):
            if descriptor_to_close < 0:
                continue
            try:
                os.close(descriptor_to_close)
            except BaseException as error:
                close_errors.append((
                    "{} close failed: {}: {}".format(
                        description, type(error).__name__, error),
                    error,
                ))
        if close_errors:
            if committed:
                post_commit_warnings.extend(
                    message for message, unused in close_errors)
                post_commit_interruptions.extend(
                    message for message, error in close_errors
                    if isinstance(error, CampaignInterrupted))
            elif not close_preserves_primary:
                raise CampaignError(
                    "; ".join(message for message, unused in close_errors))

    return {
        "path": str(path),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "committed": True,
        "post_commit_warnings": post_commit_warnings,
        "post_commit_interruptions": post_commit_interruptions,
    }


def _canonical_uint(text: str) -> int:
    if (len(text or "") > 20 or
            re.fullmatch(r"0|[1-9][0-9]*", text or "") is None):
        raise argparse.ArgumentTypeError("must be canonical unsigned decimal")
    value = int(text)
    if value > (1 << 64) - 1:
        raise argparse.ArgumentTypeError("unsigned decimal exceeds uint64")
    return value


def _canonical_sha_arg(text: str) -> str:
    if HEX64.fullmatch(text or "") is None:
        raise argparse.ArgumentTypeError("must be lower-case SHA-256")
    return text


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-receipt", type=Path, required=True)
    parser.add_argument(
        "--build-receipt-sha256", type=_canonical_sha_arg, required=True)
    parser.add_argument("--raw-log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu-a", type=_canonical_uint, required=True)
    parser.add_argument("--cpu-b", type=_canonical_uint, required=True)
    parser.add_argument("--raw-root", type=_canonical_uint, required=True)
    parser.add_argument("--loss-root", type=_canonical_uint, required=True)
    parser.add_argument(
        "--max-working-mib", type=_canonical_uint, required=True)
    parser.add_argument(
        "--timeout-seconds", type=_canonical_uint, default=1800)
    parser.add_argument(
        "--calibration-max-inner-reps",
        type=_canonical_uint, default=MAX_INNER_REPS)
    parser.add_argument(
        "--calibration-max-steps",
        type=_canonical_uint, default=MAX_CALIBRATION_STEPS)
    parser.add_argument(
        "--filler-controller",
        help="reserved interface; any value fails closed in protocol v1")
    parser.add_argument(
        "--thermal-window-controller",
        help="reserved interface; any value fails closed in protocol v1")
    return parser


def _binary_receipts(
        bindings: Dict[str, BinaryBinding]) -> Dict[str, Any]:
    return {
        treatment: {
            "treatment": binding.treatment,
            "hooks_compiled": binding.hook_state,
            "replicate": binding.replicate,
            "path": binding.path,
            "bytes": binding.bytes,
            "sha256": binding.sha256,
            "build_id": binding.build_id,
            "device": binding.device,
            "inode": binding.inode,
            "mode": binding.mode,
            "execution": {
                "contract": "sealed-authenticated-memfd-v1",
                "bytes": binding.bytes,
                "sha256": binding.sha256,
                "mode": EXECUTION_MODE,
                "seal_mask": EXECUTION_SEALS,
            },
        }
        for treatment, binding in sorted(bindings.items())
    }


def _close_raw_log(
        raw_log: AppendOnlyRawLog,
        preserve_primary_exception: bool) -> None:
    if raw_log.closed:
        return
    try:
        raw_log.close()
    except (CampaignError, OSError) as error:
        if not preserve_primary_exception:
            if isinstance(error, CampaignError):
                raise
            raise CampaignError("raw-log close failed: {}".format(error))


def _run_campaign_unclosed(
        args: argparse.Namespace,
        binding_sink: Dict[str, BinaryBinding],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    if sys.version_info < (3, 8) or sys.version_info >= (4, 0):
        raise CampaignError("Python 3.8+ and below 4.0 is required")
    validate_frozen_design()
    if args.filler_controller is not None:
        raise CampaignError(
            "filler stop/resume integration is not implemented in v1")
    if args.thermal_window_controller is not None:
        raise CampaignError(
            "thermal-window integration is not implemented in v1")
    cpu_a = _int(args.cpu_a, "cpu-a", 0, 1048575)
    cpu_b = _int(args.cpu_b, "cpu-b", 0, 1048575)
    raw_root = _int(args.raw_root, "raw-root", 0, (1 << 64) - 1)
    loss_root = _int(args.loss_root, "loss-root", 0, (1 << 64) - 1)
    max_working_mib = _int(
        args.max_working_mib, "max-working-mib", 1, 1048576)
    timeout_seconds = _int(
        args.timeout_seconds, "timeout-seconds", 1, 86400)
    calibration_max_inner_reps = _int(
        args.calibration_max_inner_reps,
        "calibration-max-inner-reps", 1, MAX_INNER_REPS)
    calibration_max_steps = _int(
        args.calibration_max_steps,
        "calibration-max-steps", 1, MAX_CALIBRATION_STEPS)

    # Preserve final components.  The descriptor-based readers and O_EXCL
    # creators must see and reject final-component symlinks themselves.
    receipt_path = Path(os.path.abspath(os.fspath(args.build_receipt)))
    raw_path = Path(os.path.abspath(os.fspath(args.raw_log)))
    output_path = Path(os.path.abspath(os.fspath(args.output)))
    if len({receipt_path, raw_path, output_path}) != 3:
        raise CampaignError("build receipt, raw log, and output must differ")
    if not raw_path.parent.is_dir() or not output_path.parent.is_dir():
        raise CampaignError(
            "raw-log and output parent directories must already exist")
    if raw_path.exists() or output_path.exists():
        raise CampaignError("raw log and output must be fresh paths")

    build_receipt, bindings, build_binding = authenticate_build_receipt(
        receipt_path, args.build_receipt_sha256, binding_sink)
    del build_receipt
    topology = validate_smt_siblings(cpu_a, cpu_b)
    helper_payload, helper_binding = _stable_file_bytes(
        Path(__file__).resolve(), 8 * 1024 * 1024, "campaign helper")
    del helper_payload
    unused_python, runner_python_binding = _stable_file_bytes(
        Path(sys.executable).resolve(), MAX_BINARY_BYTES, "runner Python")
    contract = pre_results_contract(
        args.build_receipt_sha256,
        raw_root,
        loss_root,
        max_working_mib,
        cpu_a,
        cpu_b,
        helper_binding["sha256"],
        timeout_seconds,
        calibration_max_inner_reps,
        calibration_max_steps,
        runner_python_binding["sha256"],
    )
    context_sha256 = hashlib.sha256(
        canonical_json_bytes(contract)).hexdigest()

    raw_log = AppendOnlyRawLog(raw_path)
    raw_binding = None
    try:
        runner = CampaignRunner(
            bindings=bindings,
            raw_log=raw_log,
            context_sha256=context_sha256,
            raw_root=raw_root,
            loss_root=loss_root,
            max_working_mib=max_working_mib,
            cpu_a=cpu_a,
            cpu_b=cpu_b,
            timeout_seconds=timeout_seconds,
            calibration_max_inner_reps=calibration_max_inner_reps,
            calibration_max_steps=calibration_max_steps,
        )
        preflight = runner.run_semantic_preflight()
        timed_results = [
            runner.run_timed_case(case)
            for case in frozen_timed_cases()
        ]
        if (len(preflight) != 84 or len(timed_results) != 85 or
                sum(row["kind"] == "primary" for row in timed_results) != 70 or
                sum(row["kind"] == "anchor" for row in timed_results) != 15):
            raise CampaignError("campaign result cardinality changed")
        _reverify_binaries(bindings)
        raw_binding = raw_log.close()
        if raw_binding["records"] != runner.next_invocation:
            raise CampaignError("raw log does not cover every invocation")
        execution_order_sha256 = hashlib.sha256(
            canonical_json_bytes(runner.execution_order)).hexdigest()
        receipt = {
            "schema": SCHEMA,
            "status": "complete",
            "measurement_only": True,
            "promotion_conclusion": None,
            "context_sha256": context_sha256,
            "pre_results_contract": contract,
            "campaign_helper": helper_binding,
            "python": {
                **runner_python_binding,
                "executable": str(Path(sys.executable).resolve()),
                "version": list(sys.version_info[:3]),
            },
            "build_receipt": build_binding,
            "binaries": _binary_receipts(bindings),
            "cpu_topology": topology,
            "semantic_preflight": preflight,
            "timed_results": timed_results,
            "execution": {
                "invocations": runner.next_invocation,
                "simultaneous_pairs": runner.next_pair,
                "order_sha256": execution_order_sha256,
                "order": runner.execution_order,
            },
            "raw_invocations": raw_binding,
            "statistics_contract": contract["statistics"],
        }
        payload = canonical_json_bytes(receipt)
        publication = _publish_no_replace(output_path, payload)
        return receipt, publication
    finally:
        outer_exception_active = sys.exc_info()[0] is not None
        _close_raw_log(raw_log, outer_exception_active)


def run_campaign(args: argparse.Namespace) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    bindings: Dict[str, BinaryBinding] = {}
    try:
        return _run_campaign_unclosed(args, bindings)
    finally:
        close_execution_bindings(bindings)


def _install_signal_handlers() -> Dict[int, Any]:
    previous = {}

    def interrupted(signum: int, unused_frame: Any) -> None:
        for managed in MANAGED_SIGNALS:
            signal.signal(managed, signal.SIG_IGN)
        raise CampaignInterrupted(
            "campaign interrupted by signal {}".format(signum))

    for signum in MANAGED_SIGNALS:
        previous[signum] = signal.getsignal(signum)
        signal.signal(signum, interrupted)
    return previous


def _restore_signal_handlers(previous: Dict[int, Any]) -> None:
    for signum, handler in previous.items():
        signal.signal(signum, handler)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = _parser().parse_args(argv)
    previous = _install_signal_handlers()
    try:
        unused_receipt, publication = run_campaign(args)
    finally:
        _restore_signal_handlers(previous)
    print(canonical_json_bytes({
        "committed": publication["committed"],
        "output": publication["path"],
        "post_commit_interruptions":
            publication["post_commit_interruptions"],
        "post_commit_warnings": publication["post_commit_warnings"],
        "schema": SCHEMA,
        "sha256": publication["sha256"],
    }).decode("ascii"), end="", flush=True)
    if publication["post_commit_interruptions"]:
        raise CampaignInterrupted(
            "campaign was interrupted after the receipt commit")


if __name__ == "__main__":
    try:
        main()
    except CampaignError as error:
        print("error: {}".format(error), file=sys.stderr)
        raise SystemExit(2)
