#!/usr/bin/env python3
"""Strict one-shot runner and adjudicator for the WH2 emission screen."""

from __future__ import annotations

import argparse
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
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


CONFIG_SCHEMA = "wirehair.wh2.systematic-emission-screen.config.v1"
PANEL_SCHEMA = "wirehair.wh2.systematic-emission-screen.panel.v1"
TERMINAL_SCHEMA = "wirehair.wh2.systematic-emission-screen.terminal.v1"
SUMMARY_SCHEMA = "wirehair.wh2.systematic-emission-screen.summary.v1"
DESCRIPTOR_SCHEMA = "wirehair.wh2.systematic-emission-arm.v1"
CELLS: Tuple[Tuple[int, int], ...] = ((32, 64), (1000, 64), (2048, 1280))
EXPECTED_CONFIGURATIONS: Mapping[Tuple[int, int], Mapping[str, Any]] = {
    (32, 64): {
        "K": 32,
        "block_bytes": 64,
        "dense_anchor_layout": 0,
        "dense_identity_corner": False,
        "dense_rows": 12,
        "heavy_family": 0,
        "heavy_rows": 12,
        "mix_count": 3,
        "packet_attempt": 0,
        "packet_peel_seed": 3308310999,
        "precode_attempt": 0,
        "precode_seed": 17232033509700759258,
        "source_hits": 2,
        "staircase": 13,
    },
    (1000, 64): {
        "K": 1000,
        "block_bytes": 64,
        "dense_anchor_layout": 0,
        "dense_identity_corner": False,
        "dense_rows": 12,
        "heavy_family": 0,
        "heavy_rows": 12,
        "mix_count": 3,
        "packet_attempt": 0,
        "packet_peel_seed": 1379224174,
        "precode_attempt": 0,
        "precode_seed": 10797265746045855562,
        "source_hits": 2,
        "staircase": 50,
    },
    (2048, 1280): {
        "K": 2048,
        "block_bytes": 1280,
        "dense_anchor_layout": 0,
        "dense_identity_corner": False,
        "dense_rows": 12,
        "heavy_family": 0,
        "heavy_rows": 12,
        "mix_count": 3,
        "packet_attempt": 0,
        "packet_peel_seed": 2756916472,
        "precode_attempt": 0,
        "precode_seed": 17633020099416924731,
        "source_hits": 2,
        "staircase": 54,
    },
}
EXPECTED_SOURCE_SHA256 = {
    (32, 64): "c837f2e189651cb420f6d2de1af12d8a5808cc2a9221f7ac1c2a5a12e3864e08",
    (1000, 64): "931232d4fbe6172bd35eb9ba8a3090a39403ce98028fc630c3da833c5e82e492",
    (2048, 1280): "5dddcad4ab034d18b0fa11fd096935dd3f57c2b10d4d0dac1e533c2efdeb3402",
}
EXPECTED_REPAIR_ORACLE_SHA256 = {
    (32, 64): "323a6d8e0e9462ce71f470a0a2e00994b0ea38652912c7440c777dcebed504ba",
    (1000, 64): "f7262b5a3c47b416306df10ce56a7933450098587a8980b97a99d79ebc57a6f1",
    (2048, 1280): "28b358a4dbabea124cd61f239673f078f5cbcdf5c6d92f82b0eaac56f527c154",
}
SCOPES = ("encoder-systematic", "repair-evaluate")
COMPARISONS = (
    "baseline-aa",
    "candidate-aa",
    "candidate-over-baseline",
)
REPLICATES = 12
INVOCATION_BUDGET = 65536
INTERNAL_DEADLINE_SECONDS = 115
OUTER_DEADLINE_SECONDS = 120.0
EXPECTED_PANEL_COUNT = REPLICATES * len(CELLS) * len(SCOPES) * len(COMPARISONS)
T11_975 = 2.200985160082949
MAX_INT63 = (1 << 63) - 1
MAX_STDOUT_BYTES = 8 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
LOWER40 = re.compile(r"^[0-9a-f]{40}$")

SOURCE_PATHS = (
    "CMakeLists.txt",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "bench/Wh2NativeCodecTest.cpp",
    "bench/Wh2NativePanel.cpp",
    "bench/Wh2NativePanel.h",
    "bench/Wh2SystematicEmissionScreen.cpp",
    "bench/Wh2SystematicEmissionScreen.py",
    "cmake/Wh2TimingPolicySymbolAudit.cmake",
)


class ValidationError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise ValidationError(message)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


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


def descriptor(scope: str, emission: str) -> Dict[str, Any]:
    acquisition = (
        "fixture-copy-before-clock-move-v1"
        if scope == "encoder-systematic"
        else "native-arm-copy-before-repair-clock-v1"
    )
    return {
        "arm": "wirehair2_head",
        "codec": "wirehair2_certified",
        "equation_transform": "none",
        "schema": DESCRIPTOR_SCHEMA,
        "source_acquisition": acquisition,
        "source_storage": "native-arm-owned-kxb-v1",
        "systematic_emission": emission,
    }


def descriptor_hash(value: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_bytes(value))


def validate_descriptor(
    value: Any, scope: str, emission: str, where: str
) -> Dict[str, Any]:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    expected = descriptor(scope, emission)
    if value != expected:
        fail(f"{where} differs from the frozen descriptor")
    return value


CONFIG_KEYS = {
    "cells",
    "encoder_equation_descriptor",
    "encoder_equation_descriptor_sha256",
    "encoder_retained_descriptor",
    "encoder_retained_descriptor_sha256",
    "internal_deadline_seconds",
    "invocation_budget",
    "panel_replicates",
    "repair_equation_descriptor",
    "repair_equation_descriptor_sha256",
    "repair_retained_descriptor",
    "repair_retained_descriptor_sha256",
    "schema",
    "source_git_commit",
    "target_cpu",
}
CELL_KEYS = {
    "K",
    "block_bytes",
    "construction_equivalent",
    "equation_configuration",
    "equation_configuration_sha256",
    "high_repair_controls_verified",
    "repair_oracle_sha256",
    "source_sha256",
}
CONFIGURATION_KEYS = {
    "K",
    "block_bytes",
    "dense_anchor_layout",
    "dense_identity_corner",
    "dense_rows",
    "heavy_family",
    "heavy_rows",
    "mix_count",
    "packet_attempt",
    "packet_peel_seed",
    "precode_attempt",
    "precode_seed",
    "source_hits",
    "staircase",
}
PANEL_KEYS = {
    "K",
    "block_bytes",
    "comparison",
    "invocations_per_slot",
    "left_descriptor_sha256",
    "left_direct_systematic_packets",
    "left_outcome_code",
    "order",
    "panel_key_sha256",
    "replicate",
    "right_descriptor_sha256",
    "right_direct_systematic_packets",
    "right_outcome_code",
    "schema",
    "scope",
    "slots",
    "target_cpu",
}
SLOT_KEYS = {"elapsed_ns", "side"}


def validate_config(
    config: Dict[str, Any], cpu: int, expected_commit: str
) -> Dict[Tuple[str, str], str]:
    exact_keys(config, CONFIG_KEYS, "config")
    exact_string(config["schema"], CONFIG_SCHEMA, "config.schema")
    if config["target_cpu"] != cpu or type(config["target_cpu"]) is not int:
        fail("config.target_cpu differs")
    exact_int(
        config["internal_deadline_seconds"],
        INTERNAL_DEADLINE_SECONDS,
        INTERNAL_DEADLINE_SECONDS,
        "config.internal_deadline_seconds",
    )
    exact_int(
        config["invocation_budget"],
        INVOCATION_BUDGET,
        INVOCATION_BUDGET,
        "config.invocation_budget",
    )
    exact_int(
        config["panel_replicates"],
        REPLICATES,
        REPLICATES,
        "config.panel_replicates",
    )
    exact_string(config["source_git_commit"], expected_commit, "source_git_commit")

    cells = config["cells"]
    if type(cells) is not list or len(cells) != len(CELLS):
        fail("config.cells roster differs")
    seen_configuration_hashes = set()
    for index, ((expected_k, expected_b), cell) in enumerate(zip(CELLS, cells)):
        if type(cell) is not dict:
            fail(f"config.cells[{index}] is not an object")
        exact_keys(cell, CELL_KEYS, f"config.cells[{index}]")
        exact_int(cell["K"], expected_k, expected_k, f"cell[{index}].K")
        exact_int(
            cell["block_bytes"],
            expected_b,
            expected_b,
            f"cell[{index}].block_bytes",
        )
        if cell["construction_equivalent"] is not True:
            fail(f"cell[{index}] lacks solved-state equivalence")
        if cell["high_repair_controls_verified"] is not True:
            fail(f"cell[{index}] lacks high-ID repair controls")
        configuration = cell["equation_configuration"]
        if type(configuration) is not dict:
            fail(f"cell[{index}].equation_configuration is not an object")
        exact_keys(
            configuration,
            CONFIGURATION_KEYS,
            f"cell[{index}].equation_configuration",
        )
        exact_int(configuration["K"], expected_k, expected_k, "configuration.K")
        exact_int(
            configuration["block_bytes"],
            expected_b,
            expected_b,
            "configuration.block_bytes",
        )
        expected_configuration = EXPECTED_CONFIGURATIONS[(expected_k, expected_b)]
        if canonical_bytes(configuration) != canonical_bytes(expected_configuration):
            fail("equation configuration differs from the presealed cell law")
        digest = lower_hash(
            cell["equation_configuration_sha256"],
            f"cell[{index}].equation_configuration_sha256",
        )
        if digest != sha256_bytes(canonical_bytes(configuration)):
            fail("equation configuration hash differs from its object")
        exact_string(
            cell["repair_oracle_sha256"],
            EXPECTED_REPAIR_ORACLE_SHA256[(expected_k, expected_b)],
            "cell.repair_oracle_sha256",
        )
        exact_string(
            cell["source_sha256"],
            EXPECTED_SOURCE_SHA256[(expected_k, expected_b)],
            "cell.source_sha256",
        )
        if digest in seen_configuration_hashes:
            fail("equation configuration hashes are not cell-distinct")
        seen_configuration_hashes.add(digest)

    hashes: Dict[Tuple[str, str], str] = {}
    for scope, prefix in (
        ("encoder-systematic", "encoder"),
        ("repair-evaluate", "repair"),
    ):
        equation = validate_descriptor(
            config[f"{prefix}_equation_descriptor"],
            scope,
            "equation-eval-v1",
            f"config.{prefix}_equation_descriptor",
        )
        retained = validate_descriptor(
            config[f"{prefix}_retained_descriptor"],
            scope,
            "retained-source-direct-v1",
            f"config.{prefix}_retained_descriptor",
        )
        normalized = dict(retained)
        normalized["systematic_emission"] = "equation-eval-v1"
        if normalized != equation:
            fail(f"{scope} descriptors differ outside systematic_emission")
        for emission, value in (
            ("equation", equation),
            ("retained", retained),
        ):
            actual = lower_hash(
                config[f"{prefix}_{emission}_descriptor_sha256"],
                f"config.{prefix}_{emission}_descriptor_sha256",
            )
            if actual != descriptor_hash(value):
                fail(f"{scope} {emission} descriptor hash differs")
            hashes[(scope, emission)] = actual
        if hashes[(scope, "equation")] == hashes[(scope, "retained")]:
            fail(f"{scope} descriptor hashes collide")
    return hashes


def invocations_for(k: int, replicate: int) -> int:
    total = (INVOCATION_BUDGET + k - 1) // k
    quotient, remainder = divmod(total, REPLICATES)
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
    return sha256_bytes(
        canonical_bytes(
            {
                "K": k,
                "block_bytes": block_bytes,
                "comparison": comparison,
                "scope": scope,
            }
        )
    )


def panel_order(
    k: int, block_bytes: int, comparison: str, scope: str, replicate: int
) -> str:
    phase_bit = bytes.fromhex(
        panel_key_sha256(k, block_bytes, comparison, scope)
    )[-1] & 1
    return "ABBA" if (replicate & 1) == phase_bit else "BAAB"


def validate_panel(
    panel: Dict[str, Any],
    cpu: int,
    replicate: int,
    k: int,
    block_bytes: int,
    scope: str,
    comparison: str,
    hashes: Mapping[Tuple[str, str], str],
) -> float:
    where = f"panel r{replicate} K{k} B{block_bytes} {scope} {comparison}"
    exact_keys(panel, PANEL_KEYS, where)
    exact_int(panel["K"], k, k, f"{where}.K")
    exact_int(panel["block_bytes"], block_bytes, block_bytes, f"{where}.B")
    exact_string(panel["schema"], PANEL_SCHEMA, f"{where}.schema")
    exact_string(panel["scope"], scope, f"{where}.scope")
    exact_string(panel["comparison"], comparison, f"{where}.comparison")
    exact_int(panel["replicate"], replicate, replicate, f"{where}.replicate")
    if panel["target_cpu"] != cpu or type(panel["target_cpu"]) is not int:
        fail(f"{where}.target_cpu differs")
    expected_key = panel_key_sha256(k, block_bytes, comparison, scope)
    exact_string(
        panel["panel_key_sha256"], expected_key, f"{where}.panel_key_sha256"
    )
    expected_order = panel_order(k, block_bytes, comparison, scope, replicate)
    exact_string(panel["order"], expected_order, f"{where}.order")
    expected_invocations = invocations_for(k, replicate)
    exact_int(
        panel["invocations_per_slot"],
        expected_invocations,
        expected_invocations,
        f"{where}.invocations_per_slot",
    )

    left_mode, right_mode = comparison_modes(comparison)
    exact_string(
        panel["left_descriptor_sha256"],
        hashes[(scope, left_mode)],
        f"{where}.left_descriptor_sha256",
    )
    exact_string(
        panel["right_descriptor_sha256"],
        hashes[(scope, right_mode)],
        f"{where}.right_descriptor_sha256",
    )
    expected_left_direct = k if scope == "encoder-systematic" and left_mode == "retained" else 0
    expected_right_direct = k if scope == "encoder-systematic" and right_mode == "retained" else 0
    exact_int(
        panel["left_direct_systematic_packets"],
        expected_left_direct,
        expected_left_direct,
        f"{where}.left_direct_systematic_packets",
    )
    exact_int(
        panel["right_direct_systematic_packets"],
        expected_right_direct,
        expected_right_direct,
        f"{where}.right_direct_systematic_packets",
    )
    exact_int(panel["left_outcome_code"], 0, 0, f"{where}.left_outcome_code")
    exact_int(panel["right_outcome_code"], 0, 0, f"{where}.right_outcome_code")

    slots = panel["slots"]
    sides = expected_sides(expected_order)
    if type(slots) is not list or len(slots) != 8:
        fail(f"{where}.slots does not contain eight entries")
    elapsed: List[int] = []
    for index, (slot, side) in enumerate(zip(slots, sides)):
        if type(slot) is not dict:
            fail(f"{where}.slots[{index}] is not an object")
        exact_keys(slot, SLOT_KEYS, f"{where}.slots[{index}]")
        exact_string(slot["side"], side, f"{where}.slots[{index}].side")
        elapsed.append(
            exact_int(
                slot["elapsed_ns"],
                1,
                MAX_INT63,
                f"{where}.slots[{index}].elapsed_ns",
            )
        )

    logs = [math.log(value) for value in elapsed]

    def contrast(first: int, order: str) -> float:
        if order == "ABBA":
            return (
                (logs[first] - logs[first + 1])
                + (logs[first + 3] - logs[first + 2])
            ) / 2.0
        return (
            (logs[first + 1] - logs[first])
            + (logs[first + 2] - logs[first + 3])
        ) / 2.0

    opposite_order = "BAAB" if expected_order == "ABBA" else "ABBA"
    return 0.5 * (contrast(0, expected_order) + contrast(4, opposite_order))


def sample_summary(values: Sequence[float]) -> Dict[str, Any]:
    if len(values) != REPLICATES or any(not math.isfinite(value) for value in values):
        fail("statistical group is incomplete or non-finite")
    mean = math.fsum(values) / len(values)
    variance = math.fsum(
        (value - mean) ** 2 for value in values
    ) / (len(values) - 1)
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


def validate_transcript(
    records: Sequence[Dict[str, Any]], cpu: int, expected_commit: str
) -> Tuple[Dict[str, Any], Dict[str, Any], List[str]]:
    hashes = validate_config(records[0], cpu, expected_commit)
    logs: Dict[Tuple[str, str, int], List[float]] = {}
    index = 1
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison in COMPARISONS:
                    lane = validate_panel(
                        records[index],
                        cpu,
                        replicate,
                        k,
                        block_bytes,
                        scope,
                        comparison,
                        hashes,
                    )
                    logs.setdefault((scope, comparison, k), []).append(lane)
                    index += 1
    terminal = records[index]
    exact_keys(terminal, {"panel_count", "schema", "status"}, "terminal")
    exact_int(
        terminal["panel_count"],
        EXPECTED_PANEL_COUNT,
        EXPECTED_PANEL_COUNT,
        "terminal.panel_count",
    )
    exact_string(terminal["schema"], TERMINAL_SCHEMA, "terminal.schema")
    exact_string(terminal["status"], "complete", "terminal.status")

    group_stats: Dict[str, Any] = {}
    failures: List[str] = []
    practical_log = math.log1p(0.02)
    for scope in SCOPES:
        for comparison in COMPARISONS:
            for k, block_bytes in CELLS:
                key = (scope, comparison, k)
                summary = sample_summary(logs[key])
                label = f"{scope}:{comparison}:K{k}:B{block_bytes}"
                group_stats[label] = summary
                if comparison != "candidate-over-baseline":
                    if not (
                        summary["lower95_log"] > -practical_log
                        and summary["upper95_log"] < practical_log
                    ):
                        failures.append(f"aa_ci:{label}")
                elif summary["log_mean"] > practical_log:
                    failures.append(f"cell_point:{label}")

    pooled_stats: Dict[str, Any] = {}
    for scope in SCOPES:
        per_replicate = [
            math.fsum(
                logs[(scope, "candidate-over-baseline", k)][replicate]
                for k, _ in CELLS
            ) / len(CELLS)
            for replicate in range(REPLICATES)
        ]
        summary = sample_summary(per_replicate)
        pooled_stats[scope] = summary
        if scope == "encoder-systematic":
            if not summary["upper95_log"] < math.log(0.99):
                failures.append("pooled_upper95:encoder-systematic")
        elif summary["upper95_log"] > practical_log:
            failures.append("pooled_upper95:repair-evaluate")

    statistics = {
        "failed_gates": failures,
        "groups": group_stats,
        "pooled": pooled_stats,
        "student_t_975_df11": T11_975,
    }
    return records[0], statistics, failures


def source_manifest(root: Path) -> Dict[str, Any]:
    entries: List[Dict[str, Any]] = []
    preimage = bytearray()
    for relative in SOURCE_PATHS:
        path = root / relative
        data = path.read_bytes()
        digest = sha256_bytes(data)
        entries.append({"bytes": len(data), "path": relative, "sha256": digest})
        preimage.extend(f"{digest}  {relative}\n".encode("ascii"))
    return {"entries": entries, "sha256": sha256_bytes(bytes(preimage))}


def file_sha256_fd(fd: int, size: int) -> str:
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("short read while hashing binary")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def stat_receipt(st: os.stat_result) -> Dict[str, Any]:
    return {
        "device": st.st_dev,
        "inode": st.st_ino,
        "mode": stat.S_IMODE(st.st_mode),
        "mtime_ns": st.st_mtime_ns,
        "nlink": st.st_nlink,
        "size": st.st_size,
    }


def same_file_receipt(left: os.stat_result, right: os.stat_result) -> bool:
    return stat_receipt(left) == stat_receipt(right)


def terminate_process_group(process: subprocess.Popen) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def bounded_capture(
    process: subprocess.Popen, deadline: float
) -> Tuple[bytes, bytes, bool, bool]:
    if process.stdout is None or process.stderr is None:
        fail("screen process pipes are absent")
    selector = selectors.DefaultSelector()
    streams = ((process.stdout, "stdout"), (process.stderr, "stderr"))
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    limits = {"stdout": MAX_STDOUT_BYTES, "stderr": MAX_STDERR_BYTES}
    for stream, name in streams:
        os.set_blocking(stream.fileno(), False)
        selector.register(stream, selectors.EVENT_READ, name)
    timed_out = False
    output_overflow = False
    killed_at: Optional[float] = None
    completed = False
    try:
        while selector.get_map():
            now = time.monotonic()
            if killed_at is None and now >= deadline:
                timed_out = True
                terminate_process_group(process)
                killed_at = now
            if killed_at is not None and now - killed_at >= 5.0:
                fail("screen process pipes did not close after SIGKILL")
            wait_until = killed_at + 5.0 if killed_at is not None else deadline
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
                buffers[name].extend(data)
                if len(buffers[name]) > limits[name] and killed_at is None:
                    output_overflow = True
                    terminate_process_group(process)
                    killed_at = time.monotonic()
        remaining = (
            max(0.0, killed_at + 5.0 - time.monotonic())
            if killed_at is not None
            else max(0.0, deadline - time.monotonic())
        )
        process.wait(timeout=remaining)
        completed = True
    except subprocess.TimeoutExpired:
        terminate_process_group(process)
        fail("screen process did not reap within its deadline")
    finally:
        if not completed and process.poll() is None:
            terminate_process_group(process)
            try:
                process.wait(timeout=5.0)
            except subprocess.TimeoutExpired:
                pass
        selector.close()
        for stream, _ in streams:
            if not stream.closed:
                stream.close()
    return (
        bytes(buffers["stdout"]),
        bytes(buffers["stderr"]),
        timed_out,
        output_overflow,
    )


def reserve_outputs(paths: Sequence[Path]) -> List[int]:
    if len(set(paths)) != len(paths):
        fail("output paths must be distinct")
    fds: List[int] = []
    created: List[Path] = []
    try:
        for path in paths:
            if not path.is_absolute() or not path.parent.is_dir():
                fail(f"output path is not absolute or parent is absent: {path}")
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            fds.append(fd)
            created.append(path)
        return fds
    except Exception:
        for fd in fds:
            os.close(fd)
        for path in created:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        raise


def write_all(fd: int, data: bytes) -> None:
    os.ftruncate(fd, 0)
    os.lseek(fd, 0, os.SEEK_SET)
    view = memoryview(data)
    while view:
        written = os.write(fd, view)
        if written <= 0:
            fail("short output write")
        view = view[written:]
    os.fsync(fd)
    os.fchmod(fd, 0o400)


def run_binary(
    binary: Path, cpu: int, expected_sha256: str
) -> Tuple[bytes, bytes, int, float, Dict[str, Any]]:
    path_before = binary.stat()
    if not stat.S_ISREG(path_before.st_mode) or path_before.st_nlink != 1:
        fail("screen binary must be a regular single-link file")
    if stat.S_IMODE(path_before.st_mode) & 0o022:
        fail("screen binary must not be group/other writable")
    if stat.S_IMODE(path_before.st_mode) & 0o111 == 0:
        fail("screen binary is not executable")
    fd = os.open(binary, os.O_RDONLY | os.O_CLOEXEC)
    try:
        fd_before = os.fstat(fd)
        if not same_file_receipt(path_before, fd_before):
            fail("screen binary path/FD identity differs")
        digest_before = file_sha256_fd(fd, fd_before.st_size)
        if digest_before != expected_sha256:
            fail("screen binary SHA-256 differs from the presealed value")
        fd_path = f"/proc/self/fd/{fd}"
        started = time.monotonic()
        deadline = started + OUTER_DEADLINE_SECONDS
        process = subprocess.Popen(
            [fd_path, "--run", "--cpu", str(cpu)],
            executable=fd_path,
            pass_fds=(fd,),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env={"LANG": "C", "LC_ALL": "C"},
            start_new_session=True,
        )
        stdout, stderr, timed_out, output_overflow = bounded_capture(
            process, deadline
        )
        elapsed = time.monotonic() - started
        fd_after = os.fstat(fd)
        path_after = binary.stat()
        digest_after = file_sha256_fd(fd, fd_after.st_size)
        if (
            not same_file_receipt(fd_before, fd_after)
            or not same_file_receipt(fd_after, path_after)
            or digest_before != digest_after
        ):
            fail("screen binary changed during execution")
        receipt = {
            "path": str(binary),
            "sha256": digest_before,
            "stat": stat_receipt(fd_before),
            "timed_out": timed_out,
            "output_overflow": output_overflow,
        }
        return stdout, stderr, process.returncode, elapsed, receipt
    finally:
        os.close(fd)


def git_receipt(root: Path, expected_commit: str) -> Dict[str, Any]:
    environment = {
        "GIT_OPTIONAL_LOCKS": "0",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
    }

    def git(*arguments: str) -> bytes:
        completed = subprocess.run(
            ["git", "-c", "core.fsmonitor=false", *arguments],
            cwd=root,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10.0,
            check=False,
        )
        if completed.returncode != 0 or completed.stderr:
            fail(f"static git command failed: {' '.join(arguments)}")
        return completed.stdout

    commit = git("rev-parse", "--verify", "HEAD^{commit}")
    if commit != (expected_commit + "\n").encode("ascii"):
        fail("Git HEAD differs from the expected source commit")
    tracked_status = git("status", "--porcelain=v1", "--untracked-files=no")
    if tracked_status:
        fail("tracked worktree is not clean")
    tracked_paths = git("ls-files", "--error-unmatch", "--", *SOURCE_PATHS)
    expected_paths = b"".join(
        (relative + "\n").encode("ascii") for relative in SOURCE_PATHS
    )
    if tracked_paths != expected_paths:
        fail("screen source roster is not tracked in the expected order")
    git("diff", "--quiet", "HEAD", "--", *SOURCE_PATHS)
    return {
        "head": expected_commit,
        "source_roster_sha256": sha256_bytes(tracked_paths),
        "tracked_status_sha256": sha256_bytes(tracked_status),
    }


def make_summary_preimage(summary: Dict[str, Any]) -> Dict[str, Any]:
    preimage = dict(summary)
    preimage["summary_preimage_sha256"] = None
    summary["summary_preimage_sha256"] = sha256_bytes(canonical_bytes(preimage))
    return summary


def run_once(args: argparse.Namespace) -> int:
    raw_path = Path(args.raw)
    stderr_path = Path(args.stderr)
    summary_path = Path(args.summary)
    raw_fd, stderr_fd, summary_fd = reserve_outputs(
        (raw_path, stderr_path, summary_path)
    )
    raw = b""
    stderr = b""
    exit_code = 1
    elapsed = 0.0
    binary_receipt: Dict[str, Any] = {}
    failure = "runner did not reach validation"
    source_before: Dict[str, Any] = {}
    source_after: Dict[str, Any] = {}
    git_before: Dict[str, Any] = {}
    git_after: Dict[str, Any] = {}
    statistics: Dict[str, Any] = {}
    config: Dict[str, Any] = {}
    outcome = "invalid"
    try:
        root = Path(__file__).resolve().parents[1]
        git_before = git_receipt(root, args.expected_source_commit)
        source_before = source_manifest(root)
        if source_before["sha256"] != args.expected_source_manifest_sha256:
            fail("source manifest differs from the presealed value")
        raw, stderr, exit_code, elapsed, binary_receipt = run_binary(
            Path(args.binary).resolve(), args.cpu, args.expected_binary_sha256
        )
        source_after = source_manifest(root)
        git_after = git_receipt(root, args.expected_source_commit)
        if source_before != source_after:
            fail("screen source changed during execution")
        if binary_receipt.get("timed_out"):
            fail("screen exceeded the hard 120-second deadline")
        if binary_receipt.get("output_overflow"):
            fail("screen exceeded its bounded output allowance")
        if elapsed >= OUTER_DEADLINE_SECONDS:
            fail("screen returned at or after the hard 120-second deadline")
        if exit_code != 0:
            fail(f"screen exited {exit_code}")
        if stderr != b"":
            fail("screen emitted stderr")
        records = parse_transcript(raw)
        config, statistics, failures = validate_transcript(
            records, args.cpu, args.expected_source_commit
        )
        outcome = "pass" if not failures else "reject"
        failure = "none" if not failures else "statistical-gates"
    except Exception as exc:
        failure = f"{type(exc).__name__}: {exc}"
        outcome = "invalid"
    finally:
        raw_sha = sha256_bytes(raw)
        stderr_sha = sha256_bytes(stderr)
        summary = make_summary_preimage(
            {
                "binary": binary_receipt,
                "config": config,
                "elapsed_seconds": elapsed,
                "expected_source_git_commit": args.expected_source_commit,
                "failure": failure,
                "git_after": git_after,
                "git_before": git_before,
                "outcome": outcome,
                "process_exit_code": exit_code,
                "raw_bytes": len(raw),
                "raw_record_count": raw.count(b"\n"),
                "raw_sha256": raw_sha,
                "schema": SUMMARY_SCHEMA,
                "source_manifest_after": source_after,
                "source_manifest_before": source_before,
                "statistics": statistics,
                "stderr_bytes": len(stderr),
                "stderr_sha256": stderr_sha,
                "summary_preimage_sha256": None,
                "target_cpu": args.cpu,
            }
        )
        try:
            write_all(raw_fd, raw)
            write_all(stderr_fd, stderr)
            write_all(summary_fd, canonical_bytes(summary) + b"\n")
        finally:
            os.close(raw_fd)
            os.close(stderr_fd)
            os.close(summary_fd)
    if outcome == "pass":
        return 0
    if outcome == "reject":
        return 2
    return 1


def synthetic_transcript(cpu: int, commit: str) -> bytes:
    descriptors: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for scope in SCOPES:
        descriptors[(scope, "equation")] = descriptor(scope, "equation-eval-v1")
        descriptors[(scope, "retained")] = descriptor(
            scope, "retained-source-direct-v1"
        )
    config: Dict[str, Any] = {
        "cells": [
            {
                "K": k,
                "block_bytes": block_bytes,
                "construction_equivalent": True,
                "equation_configuration": dict(EXPECTED_CONFIGURATIONS[(k, block_bytes)]),
                "equation_configuration_sha256": sha256_bytes(
                    canonical_bytes(EXPECTED_CONFIGURATIONS[(k, block_bytes)])
                ),
                "high_repair_controls_verified": True,
                "repair_oracle_sha256": EXPECTED_REPAIR_ORACLE_SHA256[(k, block_bytes)],
                "source_sha256": EXPECTED_SOURCE_SHA256[(k, block_bytes)],
            }
            for k, block_bytes in CELLS
        ],
        "encoder_equation_descriptor": descriptors[("encoder-systematic", "equation")],
        "encoder_equation_descriptor_sha256": descriptor_hash(
            descriptors[("encoder-systematic", "equation")]
        ),
        "encoder_retained_descriptor": descriptors[("encoder-systematic", "retained")],
        "encoder_retained_descriptor_sha256": descriptor_hash(
            descriptors[("encoder-systematic", "retained")]
        ),
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "panel_replicates": REPLICATES,
        "repair_equation_descriptor": descriptors[("repair-evaluate", "equation")],
        "repair_equation_descriptor_sha256": descriptor_hash(
            descriptors[("repair-evaluate", "equation")]
        ),
        "repair_retained_descriptor": descriptors[("repair-evaluate", "retained")],
        "repair_retained_descriptor_sha256": descriptor_hash(
            descriptors[("repair-evaluate", "retained")]
        ),
        "schema": CONFIG_SCHEMA,
        "source_git_commit": commit,
        "target_cpu": cpu,
    }
    records: List[Dict[str, Any]] = [config]
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison in COMPARISONS:
                    order = panel_order(k, block_bytes, comparison, scope, replicate)
                    sides = expected_sides(order)
                    left_mode, right_mode = comparison_modes(comparison)
                    left_time = 1000
                    right_time = 1000
                    if (
                        scope == "encoder-systematic"
                        and comparison == "candidate-over-baseline"
                    ):
                        left_time = 900
                    record = {
                        "K": k,
                        "block_bytes": block_bytes,
                        "comparison": comparison,
                        "invocations_per_slot": invocations_for(k, replicate),
                        "left_descriptor_sha256": descriptor_hash(
                            descriptors[(scope, left_mode)]
                        ),
                        "left_direct_systematic_packets": (
                            k
                            if scope == "encoder-systematic" and left_mode == "retained"
                            else 0
                        ),
                        "left_outcome_code": 0,
                        "order": order,
                        "panel_key_sha256": panel_key_sha256(
                            k, block_bytes, comparison, scope
                        ),
                        "replicate": replicate,
                        "right_descriptor_sha256": descriptor_hash(
                            descriptors[(scope, right_mode)]
                        ),
                        "right_direct_systematic_packets": (
                            k
                            if scope == "encoder-systematic" and right_mode == "retained"
                            else 0
                        ),
                        "right_outcome_code": 0,
                        "schema": PANEL_SCHEMA,
                        "scope": scope,
                        "slots": [
                            {
                                "elapsed_ns": left_time if side == "left" else right_time,
                                "side": side,
                            }
                            for side in sides
                        ],
                        "target_cpu": cpu,
                    }
                    records.append(record)
    records.append(
        {
            "panel_count": EXPECTED_PANEL_COUNT,
            "schema": TERMINAL_SCHEMA,
            "status": "complete",
        }
    )
    return b"".join(canonical_bytes(record) + b"\n" for record in records)


def selftest() -> int:
    cpu = 7
    commit = "1" * 40
    raw = synthetic_transcript(cpu, commit)
    records = parse_transcript(raw)
    _, statistics, failures = validate_transcript(records, cpu, commit)
    if failures or not statistics["pooled"]["encoder-systematic"]["upper95"] < 0.99:
        raise AssertionError("passing synthetic transcript did not pass")
    hashes = validate_config(records[0], cpu, commit)
    probe = json.loads(json.dumps(records[1]))
    elapsed = [1, 2, 3, 4, 5, 6, 7, 8]
    for slot, value in zip(probe["slots"], elapsed):
        slot["elapsed_ns"] = value
    lane = validate_panel(
        probe,
        cpu,
        probe["replicate"],
        probe["K"],
        probe["block_bytes"],
        probe["scope"],
        probe["comparison"],
        hashes,
    )
    logs = [math.log(value) for value in elapsed]
    if probe["order"] == "ABBA":
        primary = ((logs[0] - logs[1]) + (logs[3] - logs[2])) / 2.0
        opposite = ((logs[5] - logs[4]) + (logs[6] - logs[7])) / 2.0
    else:
        primary = ((logs[1] - logs[0]) + (logs[2] - logs[3])) / 2.0
        opposite = ((logs[4] - logs[5]) + (logs[7] - logs[6])) / 2.0
    expected_lane = (primary + opposite) / 2.0
    sides = [slot["side"] for slot in probe["slots"]]
    wrong_lane = math.fsum(
        math.log(
            math.fsum(
                elapsed[index]
                for index in range(first, first + 4)
                if sides[index] == "left"
            )
        )
        - math.log(
            math.fsum(
                elapsed[index]
                for index in range(first, first + 4)
                if sides[index] == "right"
            )
        )
        for first in (0, 4)
    ) / 2.0
    if not math.isclose(lane, expected_lane, rel_tol=0.0, abs_tol=1e-15):
        raise AssertionError("nonuniform lane contrast differs")
    if math.isclose(lane, wrong_lane, rel_tol=0.0, abs_tol=1e-12):
        raise AssertionError("lane contrast regressed to ratio-of-sums")
    bad = bytearray(raw)
    terminal = canonical_bytes(
        {"panel_count": EXPECTED_PANEL_COUNT - 1, "schema": TERMINAL_SCHEMA, "status": "complete"}
    ) + b"\n"
    bad[-len(canonical_bytes(records[-1])) - 1 :] = terminal
    try:
        validate_transcript(parse_transcript(bytes(bad)), cpu, commit)
    except ValidationError:
        pass
    else:
        raise AssertionError("bad terminal count was accepted")
    missing = raw.rsplit(b"\n", 2)[0] + b"\n"
    try:
        parse_transcript(missing)
    except ValidationError:
        pass
    else:
        raise AssertionError("missing terminal record was accepted")
    extra_record = dict(records[1])
    extra_record["unexpected"] = 1
    extra = (
        canonical_bytes(records[0])
        + b"\n"
        + canonical_bytes(extra_record)
        + b"\n"
        + b"".join(canonical_bytes(record) + b"\n" for record in records[2:])
    )
    try:
        validate_transcript(parse_transcript(extra), cpu, commit)
    except ValidationError:
        pass
    else:
        raise AssertionError("extra panel key was accepted")
    for field, value in (
        ("panel_key_sha256", "0" * 64),
        ("order", "BAAB" if records[1]["order"] == "ABBA" else "ABBA"),
        ("left_direct_systematic_packets", True),
    ):
        mutated = json.loads(json.dumps(records))
        mutated[1][field] = value
        try:
            validate_transcript(mutated, cpu, commit)
        except ValidationError:
            pass
        else:
            raise AssertionError(f"tampered panel {field} was accepted")
    typed_config = json.loads(json.dumps(records))
    typed_config[0]["cells"][0]["equation_configuration"][
        "dense_anchor_layout"
    ] = False
    typed_config[0]["cells"][0]["equation_configuration_sha256"] = sha256_bytes(
        canonical_bytes(typed_config[0]["cells"][0]["equation_configuration"])
    )
    try:
        validate_transcript(typed_config, cpu, commit)
    except ValidationError:
        pass
    else:
        raise AssertionError("bool-for-integer configuration was accepted")
    failing = [dict(record) for record in records]
    for record in failing[1:-1]:
        if (
            record["scope"] == "encoder-systematic"
            and record["comparison"] == "candidate-over-baseline"
        ):
            record["slots"] = [
                {"elapsed_ns": 1000, "side": slot["side"]}
                for slot in record["slots"]
            ]
    _, _, failed_gates = validate_transcript(failing, cpu, commit)
    if "pooled_upper95:encoder-systematic" not in failed_gates:
        raise AssertionError("non-improving candidate passed the pooled gate")
    cell_failure = json.loads(json.dumps(records))
    for record in cell_failure[1:-1]:
        if (
            record["scope"] == "encoder-systematic"
            and record["comparison"] == "candidate-over-baseline"
            and record["K"] == 32
        ):
            record["slots"] = [
                {
                    "elapsed_ns": 1030 if slot["side"] == "left" else 1000,
                    "side": slot["side"],
                }
                for slot in record["slots"]
            ]
    _, _, cell_failed_gates = validate_transcript(cell_failure, cpu, commit)
    if "cell_point:encoder-systematic:candidate-over-baseline:K32:B64" not in cell_failed_gates:
        raise AssertionError("slow encoder cell passed its point gate")
    aa_failure = json.loads(json.dumps(records))
    for record in aa_failure[1:-1]:
        if (
            record["scope"] == "encoder-systematic"
            and record["comparison"] == "baseline-aa"
            and record["K"] == 32
        ):
            record["slots"] = [
                {
                    "elapsed_ns": 1030 if slot["side"] == "left" else 1000,
                    "side": slot["side"],
                }
                for slot in record["slots"]
            ]
    _, _, aa_failed_gates = validate_transcript(aa_failure, cpu, commit)
    if "aa_ci:encoder-systematic:baseline-aa:K32:B64" not in aa_failed_gates:
        raise AssertionError("biased A/A cell passed")
    repair_failure = json.loads(json.dumps(records))
    for record in repair_failure[1:-1]:
        if (
            record["scope"] == "repair-evaluate"
            and record["comparison"] == "candidate-over-baseline"
        ):
            record["slots"] = [
                {
                    "elapsed_ns": 1030 if slot["side"] == "left" else 1000,
                    "side": slot["side"],
                }
                for slot in record["slots"]
            ]
    _, _, repair_failed_gates = validate_transcript(repair_failure, cpu, commit)
    if "pooled_upper95:repair-evaluate" not in repair_failed_gates:
        raise AssertionError("slow repair control passed")
    duplicate = b'{"cells":[],"cells":[]}\n' + raw.split(b"\n", 1)[1]
    try:
        parse_transcript(duplicate)
    except ValidationError:
        pass
    else:
        raise AssertionError("duplicate JSON key was accepted")
    # The one-shot experiment runner is intentionally Linux-only: it relies on
    # process groups, pread, and /proc/self/fd to bind the executed image.  Keep
    # the protocol/statistics selftests portable while exercising those runner
    # mechanisms only on the platform where --run is supported.
    if not sys.platform.startswith("linux"):
        print("WH2 systematic emission analyzer selftest passed")
        return 0
    sleeper = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(10)"],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    _, _, timed_out, overflow = bounded_capture(
        sleeper, time.monotonic() + 0.05
    )
    if not timed_out or overflow or sleeper.returncode is None:
        raise AssertionError("bounded timeout did not kill and reap")
    flooder = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import os; os.write(2, b'x' * 1114112)",
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    _, _, timed_out, overflow = bounded_capture(
        flooder, time.monotonic() + 5.0
    )
    if timed_out or not overflow or flooder.returncode is None:
        raise AssertionError("bounded output overflow did not kill and reap")
    with tempfile.TemporaryDirectory(prefix="wh2-systematic-selftest-") as tmp:
        root = Path(tmp)
        reserved_paths = (
            root / "raw.jsonl",
            root / "stderr.txt",
            root / "summary.json",
        )
        reserved = reserve_outputs(reserved_paths)
        try:
            try:
                reserve_outputs(reserved_paths)
            except FileExistsError:
                pass
            else:
                raise AssertionError("output reservation allowed reuse")
            write_all(reserved[0], b"sealed\n")
            if stat.S_IMODE(reserved_paths[0].stat().st_mode) != 0o400:
                raise AssertionError("sealed output mode differs")
        finally:
            for fd in reserved:
                os.close(fd)
        helper = root / "helper.py"
        helper.write_bytes(
            b"#!" + os.fsencode(sys.executable) + b"\n"
            b"import sys\n"
            b"sys.stdout.write('helper-ok\\n')\n"
        )
        helper.chmod(0o555)
        helper_digest = sha256_bytes(helper.read_bytes())
        stdout, stderr, returncode, _, receipt = run_binary(
            helper, 0, helper_digest
        )
        if (
            stdout != b"helper-ok\n"
            or stderr
            or returncode != 0
            or receipt["sha256"] != helper_digest
        ):
            raise AssertionError("sealed helper execution receipt differs")
        try:
            run_binary(helper, 0, "0" * 64)
        except ValidationError:
            pass
        else:
            raise AssertionError("wrong expected binary hash was accepted")
    print("WH2 systematic emission analyzer selftest passed")
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--binary")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--raw")
    parser.add_argument("--stderr")
    parser.add_argument("--summary")
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--expected-binary-sha256")
    parser.add_argument("--expected-source-manifest-sha256")
    args = parser.parse_args(argv)
    if args.selftest:
        return args
    if (
        not args.binary
        or args.cpu is None
        or args.cpu < 0
        or not args.raw
        or not args.stderr
        or not args.summary
        or type(args.expected_source_commit) is not str
        or LOWER40.fullmatch(args.expected_source_commit) is None
        or type(args.expected_binary_sha256) is not str
        or LOWER64.fullmatch(args.expected_binary_sha256) is None
        or type(args.expected_source_manifest_sha256) is not str
        or LOWER64.fullmatch(args.expected_source_manifest_sha256) is None
    ):
        parser.error(
            "--run requires binary, CPU, three outputs, and commit/binary/source seals"
        )
    return args


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    return run_once(args)


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except (ValidationError, AssertionError, OSError) as exc:
        print(f"systematic emission analyzer failed: {exc}", file=sys.stderr)
        sys.exit(1)
