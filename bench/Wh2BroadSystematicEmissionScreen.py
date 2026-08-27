#!/usr/bin/env python3
"""Strict one-shot runner and adjudicator for the broad WH2 emission screen."""

from __future__ import annotations

import argparse
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


CONFIG_SCHEMA = "wirehair.wh2.broad-systematic-emission-screen.config.v1"
PANEL_SCHEMA = "wirehair.wh2.broad-systematic-emission-screen.panel.v1"
TERMINAL_SCHEMA = "wirehair.wh2.broad-systematic-emission-screen.terminal.v1"
SUMMARY_SCHEMA = "wirehair.wh2.broad-systematic-emission-screen.summary.v1"
DESCRIPTOR_SCHEMA = "wirehair.wh2.broad-systematic-emission-arm.v1"

CELLS: Tuple[Tuple[int, int], ...] = tuple(
    (k, block_bytes)
    for k in (8, 128, 512, 5000, 64000)
    for block_bytes in (64, 1280)
)
PANEL_SCOPES = ("fresh-systematic", "fresh-repair", "steady-repair")
METRICS: Tuple[Tuple[str, str, str], ...] = (
    ("fresh-systematic-total", "fresh-systematic", "elapsed_ns"),
    ("fresh-repair-total", "fresh-repair", "elapsed_ns"),
    ("fresh-repair-init-first", "fresh-repair", "nested_elapsed_ns"),
    ("steady-repair-total", "steady-repair", "elapsed_ns"),
)
COMPARISONS = ("baseline-aa", "candidate-aa", "candidate-over-baseline")
REPLICATES = 12
INVOCATION_BUDGET = 65536
MINIMUM_INVOCATIONS = 24
INTERNAL_DEADLINE_SECONDS = 115
OUTER_DEADLINE_SECONDS = 120.0
EXPECTED_PANEL_COUNT = (
    REPLICATES * len(CELLS) * len(PANEL_SCOPES) * len(COMPARISONS)
)
T11_975 = 2.200985160082949
MAX_INT63 = (1 << 63) - 1
MAX_UINT32 = (1 << 32) - 1
MAX_UINT64 = (1 << 64) - 1
MAX_STDOUT_BYTES = 16 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
MAX_FAILURE_TEXT_CHARS = 64 * 1024
MAX_FAILED_GATES = 174
MAX_PUBLICATION_FAILURES = 64
EXPECTED_TARGET_CPU = 120
EXPECTED_TARGET_CORE = (0, 56)
EXPECTED_TARGET_THREADS = (56, 120)
SIBLING_NON_IDLE_TICK_CAP = 5
THERMAL_MAX_MILLIC = 85000
HEALTH_SCHEMA = "wirehair.wh2.broad-systematic-emission-health.v1"
HEALTH_LOADER_SCHEMA = (
    "wirehair.wh2.broad-systematic-emission-health-source-loader.v1"
)
SAMPLER_SCHEMA = "wirehair.wh2.sampler-attestation.v1"
CPU_TICK_RECEIPT_KEYS = {
    "cpu", "non_idle_ticks", "read_monotonic_ns", "tick_fields",
}
CPU_TICK_FIELD_KEYS = {
    "idle", "iowait", "irq", "nice", "softirq", "steal", "system", "user",
}
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
LOWER40 = re.compile(r"^[0-9a-f]{40}$")

# Frozen from the native untimed selftest receipt.  These maps deliberately
# cover the exact Cartesian roster.
EXPECTED_CONFIGURATION_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "f696d3f825d9d8b56d4f667d423f1fca0d76db5cea156b2bae18ed27e5e3369e",
    (8, 1280): "39b3f23adeabdb6de6f63144e74008b526d653ff06fd01748fe0b1f07a7eed24",
    (128, 64): "1cecbe75d7c35c5f4413ba602e8919f6bdc90a6c8b1dfdc467894ba5c9270572",
    (128, 1280): "9e1ee0fa094d14e1ced2226af42fa762405df1e4f59a4f212bdd53a522617951",
    (512, 64): "5f201020c4e541b37c31e2cecb20088e1c77794ebb266a5d8f9f42e412bfac71",
    (512, 1280): "f3772040f91e7a58eb7fab8b4be44a726afed8b606ab0d21af3b007faeea3af8",
    (5000, 64): "96b443b6fe64b4a6d29ca4f307f3372ed13951356e3de9a568952ce7f9d5829b",
    (5000, 1280): "ad0769817c845f6ed6eda1ada8a5c7d048e7ee31f9303302485831983eea2696",
    (64000, 64): "112530ebf8646e224936bfb32552591dbecb0cb133df02e0f841b61eccc267bc",
    (64000, 1280): "6ef35f5bf51ad0688504b1590c179c476c9d6cd31a73424ea8be364f692d0455",
}
EXPECTED_SOURCE_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "da7d0d04aaeb241cc351650a839fd974f99be26766b7697de0af6b23a72c1276",
    (8, 1280): "b609fa6e4de08903bd33716155bf3936f53fe9cb14fe4ad9c44a1cbb55604adb",
    (128, 64): "24389e04d91e44f451a9d8f801be054aa6c78485819eb25676838ae2e9765f4d",
    (128, 1280): "d8216125397beca9146bb3dddc87a904adca60fcee659506b8546f29a2a1ebb9",
    (512, 64): "5bf428d7f837c9e3983f7116d6dd1516e9abfdc8a0fde2ede78b837653195a26",
    (512, 1280): "4fcecbc851e2ab882f98c9c831eae9045cfc481d69232081964f76c704f7a1ac",
    (5000, 64): "df192ef37e7f2a187108d7659714570d070b2dcbae74e11f437bf53b64a0c840",
    (5000, 1280): "a61f8833c1deff31737a411237b47cb2494469dba9ffbe85e3dab7d102beb5f9",
    (64000, 64): "a37eee40a067926c7979a159ddaea717cc96a21090d67b712cc8db2d7df16f18",
    (64000, 1280): "8a0df20ca26df770afc859c858877acf95c435926a7a57d8e238f3ecf076c80d",
}
EXPECTED_REPAIR_ORACLE_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "75112775332232f5c4f2e70ef9dcfbabd41afc1aa89217374d2fbf025c42ddbf",
    (8, 1280): "893f50a4b50f6136d515ec44229a0808adf006dd41f25db406dd92c0233787e3",
    (128, 64): "04bc95892a4eee3afea70fc33030735a76c5a69af2738dd12c626bb6fe721b47",
    (128, 1280): "5353ff791ae97b657d7d0fde48e6ea537bca80097bb6de65253ecb66cf76d6a4",
    (512, 64): "17f5b2ca5321774bb9740b0cfcbb4e06053a1af88800928b067f389d037c0dbb",
    (512, 1280): "e75a45a5d8af915c4f24dd61e2348c16474d7452d8e238555b6b83bd287000b8",
    (5000, 64): "da5a781cabfcacac7b765c2cd9d87f00e840a33df943619b81a37d6d314a6ee8",
    (5000, 1280): "cc0ab889b23c9f2a440e9edf46a1f85e4c09e1a1a9620952cdd0791058f50acd",
    (64000, 64): "38f7f201fac458b70d3d73e6a8d2c8703119c7bf91f3b8bc949352564c09d18e",
    (64000, 1280): "a6083d31b6ba6b2dc0e4e5214317b2a16e8b3ce0a7b0643fedd8fc7a5807e4cd",
}
EXPECTED_FIRST_REPAIR_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "422a576a15228c4fe99ae3d72b220c4b9c3eff12981778733ffcda41ff86c4ab",
    (8, 1280): "5f5246fb387706ca99faa6e958c916f8882cf3b70339fa249984e5690012b94d",
    (128, 64): "aee24d5a978e4618e4dd9b86f9f8a82273e1200a9ad817fe89c954565c9b2438",
    (128, 1280): "e7f4ec656551ae770f10245da4c187b7a16c61d70ffcfc8e7877898c47e404e3",
    (512, 64): "4b78dcec571a957965b27e22fc213b87afb3e5c52059a9fa128fd7946bff67e5",
    (512, 1280): "f063214b55a1be96d89bf6f49d7c212d99c30ae6f835d1bd8c2b50cde9aca4b7",
    (5000, 64): "061d006b12581421ae9f1cf0f50de6284e138c716b5b01cd29862266f48c05e0",
    (5000, 1280): "951d348686c4c5b249e009e7e22584757ebfe123b42d0a2260de6bb0453d97b3",
    (64000, 64): "f667a2be87515e63e70f9597314fc16c98e0c33e8c0f3870a6a73c62a9134219",
    (64000, 1280): "fb85f157c691fda71203e13a2ce775a12ef240aadcfd9156efbbe0b115fb4db4",
}
EXPECTED_HIGH_ID_ORACLE_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "cea33c827c0102fffd20f8d6ddefb39474f8a63ce70a6fae17d5f0eba51eb608",
    (8, 1280): "cbd3ea46f46890873ba75bf0cef18b0981ac3f32a526550eff2d311b0ae2ba42",
    (128, 64): "6ab0519fc415a5c6e3832e6e30df31ca7d23b4be3340acd4611ecdaf87aca051",
    (128, 1280): "ccb3d70885fca5b305a733c90d9c2b086168b1d689323eb7f10929dc63ce5385",
    (512, 64): "0bd4ee7de5e9d967adf465d0b1795eec26f0dc10a5764fe5f2b6489e3b2c3b5a",
    (512, 1280): "93f9a3ebf6c1020620c7384e1cc1db78a81f9ec9353e1240f609489e731638a2",
    (5000, 64): "79b7f9d4c4cf2232e750da5661d8cf43d4b38eb9fdb645030d832e77e72dd2ef",
    (5000, 1280): "6a89ce732b75af1e513c404f6ede3808d020f9fb3fcbabea8b838a48589fb42e",
    (64000, 64): "3d2c21f2ddfb32bb088193cf9daac9426ba26bb2d0afbd66f3e41d4c8ffa4e8b",
    (64000, 1280): "ef170b0bcf6ba809a25a47369018782eac71c2c89d33052268c574d397fd9662",
}
EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256: Mapping[Tuple[int, int], str] = {
    (8, 64): "c765f5bfe152e6eb9376bbfcc82eda06164077650e28aa12c92170d57dc3f534",
    (8, 1280): "123f643f8c36098aaa2b2a47c4408bf120f0d31a5a66b12bc85defc97a57859b",
    (128, 64): "398551197d004d987ed033e3f2a30e7289de0e6db990128cfdae2eff1519c28a",
    (128, 1280): "3e4b62fb4a5b4dfb35045eb0c8c1130c9321fe2bb6fde029ff7da498a02cceb8",
    (512, 64): "8f3dc098fc821babc939505b0d7961dac12f799be6f24604965904055ce92983",
    (512, 1280): "e1661d99ce94c5c007fbad170e5c63424c258f647e5318d7ab153a84e94f5165",
    (5000, 64): "fce0ca96ed238b58804f7126625150b18cdc9f5c2b6bf1b252ff856b2177c0d2",
    (5000, 1280): "942cf4574936e99f0af206d3eb0be6e23aefc4995528ff749d40e19c3343635a",
    (64000, 64): "e5d820cd7a0fc7ba09dd4d7aa79332f40e5be077dcbf978ae8eae2ec958fc76c",
    (64000, 1280): "220c9abba2e111e616f9e729eff286aa2d8bdb2a9a8b4c1dcfdc83c321b68b13",
}

SOURCE_PATHS = (
    "CMakeLists.txt",
    "bench/Wh2BroadSystematicEmissionScreen.cpp",
    "bench/Wh2BroadSystematicEmissionScreen.py",
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
    "cmake/Wh2BroadSystematicEmissionSymbolAudit.cmake",
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


def metric_scope(metric: str) -> str:
    for expected_metric, scope, _ in METRICS:
        if metric == expected_metric:
            return scope
    fail("unknown metric")
    raise AssertionError


def metric_field(metric: str) -> str:
    for expected_metric, _, field in METRICS:
        if metric == expected_metric:
            return field
    fail("unknown metric")
    raise AssertionError


def descriptor(metric: str, emission: str) -> Dict[str, Any]:
    timed_work = {
        "fresh-systematic-total": (
            "fresh-eager-init-plus-systematic-ids-0-through-k-minus-1-v1"
        ),
        "fresh-repair-total": (
            "fresh-eager-init-plus-repair-ids-k-through-2k-minus-1-v1"
        ),
        "fresh-repair-init-first": (
            "nested-fresh-eager-init-plus-first-repair-id-k-v1"
        ),
        "steady-repair-total": (
            "prebuilt-repair-ids-k-through-2k-minus-1-v1"
        ),
    }
    acquisition = {
        "fresh-systematic-total": "fixture-copy-before-clock-move-v1",
        "fresh-repair-total": "source-copy-in-fresh-init-clock-v1",
        "fresh-repair-init-first": "source-copy-in-fresh-init-clock-v1",
        "steady-repair-total": "native-arm-copy-before-steady-clock-v1",
    }
    if metric not in timed_work or emission not in (
        "equation-eval-v1", "retained-source-direct-v1"
    ):
        fail("invalid descriptor request")
    return {
        "arm": "wirehair2_head",
        "codec": "wirehair2_certified",
        "equation_transform": "none",
        "metric": metric,
        "schema": DESCRIPTOR_SCHEMA,
        "source_acquisition": acquisition[metric],
        "source_storage": "native-arm-owned-kxb-v1",
        "systematic_emission": emission,
        "timed_work": timed_work[metric],
    }


def descriptor_hash(value: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_bytes(value))


CONFIG_KEYS = {
    "cells", "descriptors", "internal_deadline_seconds", "invocation_budget",
    "minimum_invocations", "panel_replicates", "schema", "source_git_commit",
    "target_cpu",
}
CELL_KEYS = {
    "K", "block_bytes", "construction_equivalent", "equation_configuration",
    "equation_configuration_sha256", "first_repair_sha256",
    "high_id_oracle_sha256", "high_repair_controls_verified",
    "repair_oracle_sha256", "solved_state_equivalence_receipt_sha256",
    "source_sha256",
}
CONFIGURATION_KEYS = {
    "K", "block_bytes", "dense_anchor_layout", "dense_identity_corner",
    "dense_rows", "heavy_family", "heavy_rows", "mix_count", "packet_attempt",
    "packet_peel_seed", "precode_attempt", "precode_seed", "source_hits",
    "staircase",
}
DESCRIPTOR_ENTRY_KEYS = {"descriptor", "descriptor_sha256", "mode"}
PANEL_KEYS = {
    "K", "block_bytes", "comparison", "invocations_per_slot",
    "left_descriptor_sha256", "left_direct_systematic_packets",
    "left_nested_descriptor_sha256", "left_outcome_code", "nested_metric",
    "order", "panel_key_sha256",
    "primary_metric", "replicate", "right_descriptor_sha256",
    "right_direct_systematic_packets", "right_nested_descriptor_sha256",
    "right_outcome_code", "schema", "scope", "slots", "target_cpu",
}
SLOT_KEYS = {"elapsed_ns", "nested_elapsed_ns", "side"}


def solved_state_equivalence_receipt(
    cell: Mapping[str, Any]
) -> Dict[str, str]:
    return {
        "configuration_sha256": cell["equation_configuration_sha256"],
        "first_repair_sha256": cell["first_repair_sha256"],
        "high_id_oracle_sha256": cell["high_id_oracle_sha256"],
        "repair_oracle_sha256": cell["repair_oracle_sha256"],
        "source_sha256": cell["source_sha256"],
        "systematic_oracle_sha256": cell["source_sha256"],
    }


def validate_config(
    config: Dict[str, Any], cpu: int, expected_commit: str
) -> Dict[Tuple[str, str], str]:
    exact_keys(config, CONFIG_KEYS, "config")
    exact_string(config["schema"], CONFIG_SCHEMA, "config.schema")
    exact_int(config["target_cpu"], cpu, cpu, "config.target_cpu")
    exact_int(
        config["internal_deadline_seconds"],
        INTERNAL_DEADLINE_SECONDS,
        INTERNAL_DEADLINE_SECONDS,
        "config.internal_deadline_seconds",
    )
    exact_int(
        config["invocation_budget"], INVOCATION_BUDGET, INVOCATION_BUDGET,
        "config.invocation_budget",
    )
    exact_int(
        config["minimum_invocations"], MINIMUM_INVOCATIONS, MINIMUM_INVOCATIONS,
        "config.minimum_invocations",
    )
    exact_int(
        config["panel_replicates"], REPLICATES, REPLICATES,
        "config.panel_replicates",
    )
    exact_string(config["source_git_commit"], expected_commit, "source_git_commit")

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
        if cell["construction_equivalent"] is not True:
            fail(f"{where} lacks solved-state equivalence")
        if cell["high_repair_controls_verified"] is not True:
            fail(f"{where} lacks high-ID repair controls")
        configuration = cell["equation_configuration"]
        if type(configuration) is not dict:
            fail(f"{where}.equation_configuration is not an object")
        exact_keys(configuration, CONFIGURATION_KEYS, f"{where}.configuration")
        exact_int(configuration["K"], expected_k, expected_k, "configuration.K")
        exact_int(
            configuration["block_bytes"], expected_b, expected_b,
            "configuration.block_bytes",
        )
        # Reject Python bools in all integer-valued configuration fields.
        for field in CONFIGURATION_KEYS - {
            "dense_identity_corner", "precode_seed"
        }:
            exact_int(
                configuration[field], 0, MAX_UINT32,
                f"configuration.{field}",
            )
        exact_int(
            configuration["precode_seed"], 0, MAX_UINT64,
            "configuration.precode_seed",
        )
        if type(configuration["dense_identity_corner"]) is not bool:
            fail("configuration.dense_identity_corner is not bool")
        key = (expected_k, expected_b)
        config_hash = lower_hash(
            cell["equation_configuration_sha256"], f"{where}.configuration hash"
        )
        if config_hash != sha256_bytes(canonical_bytes(configuration)):
            fail(f"{where} configuration hash differs from its object")
        exact_string(
            config_hash, EXPECTED_CONFIGURATION_SHA256[key],
            f"{where}.equation_configuration_sha256",
        )
        for field, expected_map in (
            ("source_sha256", EXPECTED_SOURCE_SHA256),
            ("repair_oracle_sha256", EXPECTED_REPAIR_ORACLE_SHA256),
            ("first_repair_sha256", EXPECTED_FIRST_REPAIR_SHA256),
            ("high_id_oracle_sha256", EXPECTED_HIGH_ID_ORACLE_SHA256),
            (
                "solved_state_equivalence_receipt_sha256",
                EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256,
            ),
        ):
            exact_string(cell[field], expected_map[key], f"{where}.{field}")
            lower_hash(cell[field], f"{where}.{field}")
        if cell["solved_state_equivalence_receipt_sha256"] != sha256_bytes(
            canonical_bytes(solved_state_equivalence_receipt(cell))
        ):
            fail(f"{where} solved-state receipt differs")

    entries = config["descriptors"]
    if type(entries) is not list or len(entries) != 2 * len(METRICS):
        fail("config.descriptors roster differs")
    hashes: Dict[Tuple[str, str], str] = {}
    index = 0
    for metric, _, _ in METRICS:
        values: Dict[str, Dict[str, Any]] = {}
        for mode, emission in (
            ("equation", "equation-eval-v1"),
            ("retained", "retained-source-direct-v1"),
        ):
            entry = entries[index]
            index += 1
            where = f"config.descriptors[{index - 1}]"
            if type(entry) is not dict:
                fail(f"{where} is not an object")
            exact_keys(entry, DESCRIPTOR_ENTRY_KEYS, where)
            exact_string(entry["mode"], mode, f"{where}.mode")
            expected = descriptor(metric, emission)
            if entry["descriptor"] != expected:
                fail(f"{where}.descriptor differs")
            digest = lower_hash(entry["descriptor_sha256"], f"{where}.sha256")
            if digest != descriptor_hash(expected):
                fail(f"{where}.descriptor hash differs")
            hashes[(metric, mode)] = digest
            values[mode] = expected
        normalized = dict(values["retained"])
        normalized["systematic_emission"] = "equation-eval-v1"
        if normalized != values["equation"]:
            fail(f"{metric} descriptors differ outside emission")
        if hashes[(metric, "equation")] == hashes[(metric, "retained")]:
            fail(f"{metric} descriptor hashes collide")
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
        "K": k, "block_bytes": block_bytes, "comparison": comparison,
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
    primary = next(metric for metric, metric_scope_name, _ in METRICS
                   if metric_scope_name == scope and metric != "fresh-repair-init-first")
    nested = "fresh-repair-init-first" if scope == "fresh-repair" else "none"
    exact_string(panel["primary_metric"], primary, f"{where}.primary_metric")
    exact_string(panel["nested_metric"], nested, f"{where}.nested_metric")
    left_mode, right_mode = comparison_modes(comparison)
    exact_string(
        panel["left_descriptor_sha256"], hashes[(primary, left_mode)],
        f"{where}.left_descriptor_sha256",
    )
    exact_string(
        panel["right_descriptor_sha256"], hashes[(primary, right_mode)],
        f"{where}.right_descriptor_sha256",
    )
    if scope == "fresh-repair":
        exact_string(
            panel["left_nested_descriptor_sha256"],
            hashes[("fresh-repair-init-first", left_mode)],
            f"{where}.left_nested_descriptor_sha256",
        )
        exact_string(
            panel["right_nested_descriptor_sha256"],
            hashes[("fresh-repair-init-first", right_mode)],
            f"{where}.right_nested_descriptor_sha256",
        )
    elif (
        panel["left_nested_descriptor_sha256"] is not None
        or panel["right_nested_descriptor_sha256"] is not None
    ):
        fail(f"{where} has a nested descriptor outside fresh repair")
    expected_left_direct = (
        k if scope == "fresh-systematic" and left_mode == "retained" else 0
    )
    expected_right_direct = (
        k if scope == "fresh-systematic" and right_mode == "retained" else 0
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
    primary_elapsed: List[int] = []
    nested_elapsed: List[int] = []
    for index, (slot, side) in enumerate(zip(slots, sides)):
        if type(slot) is not dict:
            fail(f"{where}.slots[{index}] is not an object")
        exact_keys(slot, SLOT_KEYS, f"{where}.slots[{index}]")
        exact_string(slot["side"], side, f"{where}.slots[{index}].side")
        total = exact_int(
            slot["elapsed_ns"], 1, MAX_INT63, f"{where}.slots[{index}].elapsed_ns"
        )
        primary_elapsed.append(total)
        if scope == "fresh-repair":
            nested_value = exact_int(
                slot["nested_elapsed_ns"], 1, total,
                f"{where}.slots[{index}].nested_elapsed_ns",
            )
            nested_elapsed.append(nested_value)
        else:
            exact_int(
                slot["nested_elapsed_ns"], 0, 0,
                f"{where}.slots[{index}].nested_elapsed_ns",
            )
    result = {primary: lane_contrast(primary_elapsed, order)}
    if scope == "fresh-repair":
        result["fresh-repair-init-first"] = lane_contrast(nested_elapsed, order)
    return result


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
    if k in (8, 128):
        return "small"
    if k in (512, 5000):
        return "medium"
    if k == 64000:
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
                    index += 1
    terminal = records[index]
    exact_keys(terminal, {"panel_count", "schema", "status"}, "terminal")
    exact_int(
        terminal["panel_count"], EXPECTED_PANEL_COUNT, EXPECTED_PANEL_COUNT,
        "terminal.panel_count",
    )
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
                    if metric != "fresh-systematic-total" and not (
                        summary["upper95_log"] < practical_log
                    ):
                        failures.append(f"cell_upper95:{label}")

    aggregates: Dict[str, Any] = {}
    for metric, _, _ in METRICS:
        limit = 0.99 if metric == "fresh-systematic-total" else 1.02
        limit_log = math.log(limit)
        groups: List[Tuple[str, Sequence[Tuple[int, int]]]] = [
            ("pooled", CELLS),
            ("width-B64", tuple(cell for cell in CELLS if cell[1] == 64)),
            ("width-B1280", tuple(cell for cell in CELLS if cell[1] == 1280)),
        ]
        groups.extend(
            (f"band-{band}", tuple(cell for cell in CELLS if band_for(cell[0]) == band))
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


def source_manifest(root: Path) -> Dict[str, Any]:
    entries: List[Dict[str, Any]] = []
    preimage = bytearray()
    for relative in SOURCE_PATHS:
        data = (root / relative).read_bytes()
        digest = sha256_bytes(data)
        entries.append({"bytes": len(data), "path": relative, "sha256": digest})
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
            data = (root / relative).read_bytes()
            expected = by_path[relative]
            digest = sha256_bytes(data)
            if len(data) != expected["bytes"] or digest != expected["sha256"]:
                fail(f"health helper source bytes changed: {relative}")
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
        "device": st.st_dev, "inode": st.st_ino, "mode": stat.S_IMODE(st.st_mode),
        "mtime_ns": st.st_mtime_ns, "nlink": st.st_nlink, "size": st.st_size,
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
        self.logical_sealed = False

    def _handler(self, signum: int, _frame: Any) -> None:
        if self.logical_sealed:
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
        if not self.entered or self.logical_sealed:
            fail("signal guard commit transition is invalid")
        if not self.seal_blocked:
            signal.pthread_sigmask(signal.SIG_BLOCK, self.SIGNALS)
            self.seal_blocked = True

    def collect_pending(self) -> bool:
        if not self.seal_blocked or self.logical_sealed:
            fail("pending signals require an uncommitted blocked guard")
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

    def commit_logical_seal(self) -> None:
        if not self.seal_blocked or self.logical_sealed:
            fail("signal guard logical commit transition is invalid")
        self.logical_sealed = True

    def __exit__(self, _kind: Any, _value: Any, _traceback: Any) -> None:
        if not self.entered:
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
) -> CaptureResult:
    result = CaptureResult()
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
    limits = {"stdout": MAX_STDOUT_BYTES, "stderr": MAX_STDERR_BYTES}
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
            if killed_at is not None and now - killed_at >= 5.0:
                errors.append("screen process pipes did not close after SIGKILL")
                break
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
            max(0.0, killed_at + 5.0 - time.monotonic())
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
                process.wait(timeout=5.0)
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

    def all_failures(self) -> List[str]:
        return list(self.precommit_failures + self.postcommit_failures)

    def merged(self, other: "BundleFinishResult") -> "BundleFinishResult":
        return BundleFinishResult(
            logical_commit_succeeded=(
                self.logical_commit_succeeded
                or other.logical_commit_succeeded
            ),
            precommit_failures=tuple(dict.fromkeys(
                self.precommit_failures + other.precommit_failures
            )),
            postcommit_failures=tuple(dict.fromkeys(
                self.postcommit_failures + other.postcommit_failures
            )),
        )


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
    "capture_error", "child_reap_monotonic_ns", "child_start_monotonic_ns",
    "expected_sha256", "output_overflow", "path", "path_stat_after",
    "process_started", "sha256_after", "sha256_before", "stat_after",
    "stat_before", "timed_out",
}
STAT_RECEIPT_KEYS = {
    "device", "inode", "mode", "mtime_ns", "nlink", "size",
}
OUTPUT_BUNDLE_KEYS = {"directory", "files", "parent", "path"}
OUTPUT_DIRECTORY_KEYS = {
    "device", "inode", "nlink", "reserved_mode", "sealed_mode",
}
OUTPUT_FILE_KEYS = OUTPUT_DIRECTORY_KEYS
OUTPUT_PARENT_KEYS = {"device", "initial_nlink", "inode", "mode"}
HEALTH_KEYS = {
    "child_reap_monotonic_ns", "child_start_monotonic_ns",
    "controller_core", "controller_cpu", "controller_initial_affinity",
    "controller_singleton_affinity_end", "edac_policy", "receipt_sha256",
    "sampler", "sampler_core", "sampler_cpu", "sampler_receipt_sha256",
    "schema", "sibling_non_idle_tick_cap", "sibling_tick_policy",
    "sibling_ticks", "target_core", "target_cpu", "target_threads",
    "thermal", "thermal_max_millic", "terminal_status",
}
SAMPLER_RECEIPT_KEYS = {
    "cpu", "csv_device", "csv_inode", "csv_path", "pid",
    "process_start_ticks", "schema", "script_path", "script_sha256",
    "terminal_status", "window_end_monotonic_ns", "window_start_monotonic_ns",
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
SOURCE_MANIFEST_ENTRY_KEYS = {"bytes", "path", "sha256"}
GIT_RECEIPT_KEYS = {
    "head", "source_roster_sha256", "tracked_status_sha256",
}
THERMAL_RECEIPT_KEYS = {
    "cpu", "cpu_tctl_max_millic", "csv_device", "csv_inode", "csv_path",
    "dimm_max_millic", "dimm_read_errors", "edac_ce_max", "edac_ue_max",
    "pid", "process_start_ticks", "sample_count", "script_path",
    "script_sha256", "terminal_status", "window_csv_sha256",
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
    "pooled", "width-B64", "width-B1280", "band-small", "band-medium",
    "band-large",
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


def validate_stat_receipt(receipt: Any, where: str) -> None:
    if type(receipt) is not dict:
        fail(f"{where} is not an object")
    exact_keys(receipt, STAT_RECEIPT_KEYS, where)
    exact_int(receipt["device"], 0, MAX_UINT64, f"{where}.device")
    exact_int(receipt["inode"], 1, MAX_UINT64, f"{where}.inode")
    exact_int(receipt["mode"], 0, 0o7777, f"{where}.mode")
    exact_int(
        receipt["mtime_ns"], -(1 << 63), MAX_INT63, f"{where}.mtime_ns"
    )
    exact_int(receipt["nlink"], 1, MAX_UINT64, f"{where}.nlink")
    exact_int(receipt["size"], 0, MAX_INT63, f"{where}.size")


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
        exact_int(entry["bytes"], 0, MAX_INT63, f"{entry_where}.bytes")
        digest = lower_hash(entry["sha256"], f"{entry_where}.sha256")
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
        exact_int(module_receipt["bytes"], 0, MAX_INT63, f"{where}.bytes")
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
    exact_string(receipt["head"], expected_commit, f"{where}.head")
    roster_hash = lower_hash(
        receipt["source_roster_sha256"], f"{where}.source roster"
    )
    status_hash = lower_hash(
        receipt["tracked_status_sha256"], f"{where}.status"
    )
    expected_roster = b"".join(
        (relative + "\n").encode("ascii") for relative in SOURCE_PATHS
    )
    if (
        roster_hash != sha256_bytes(expected_roster)
        or status_hash != sha256_bytes(b"")
    ):
        fail(f"{where} clean tracked-source binding differs")


def validate_sample_summary(value: Any, where: str) -> None:
    if type(value) is not dict:
        fail(f"{where} is not an object")
    exact_keys(value, SAMPLE_SUMMARY_KEYS, where)
    exact_int(value["n"], REPLICATES, REPLICATES, f"{where}.n")
    for name in SAMPLE_SUMMARY_KEYS - {"n"}:
        number = value[name]
        if type(number) not in (int, float) or not math.isfinite(number):
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
    if type(student_t) not in (int, float) or student_t != T11_975:
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
    if binary["process_started"]:
        exact_int(start, 0, MAX_INT63, "binary child start timestamp")
        if reap is not None:
            exact_int(reap, start, MAX_INT63, "binary child reap timestamp")
    elif start is not None or reap is not None:
        fail("unstarted binary has child timestamps")
    for name in ("stat_before", "stat_after", "path_stat_after"):
        receipt = binary[name]
        if type(receipt) is not dict:
            fail(f"binary {name} is not an object")
        if receipt:
            validate_stat_receipt(receipt, f"binary {name}")
    if require_complete:
        if (
            not binary["process_started"]
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


def validate_sampler_receipt(sampler: Any) -> None:
    if type(sampler) is not dict:
        fail("health sampler receipt is not an object")
    exact_keys(sampler, SAMPLER_RECEIPT_KEYS, "health sampler")
    exact_string(sampler["schema"], SAMPLER_SCHEMA, "sampler schema")
    exact_string(sampler["terminal_status"], "ok", "sampler status")
    exact_int(sampler["pid"], 1, MAX_INT63, "sampler pid")
    exact_int(sampler["cpu"], 0, MAX_INT63, "sampler cpu")
    exact_int(
        sampler["process_start_ticks"], 1, MAX_UINT64,
        "sampler process start ticks",
    )
    exact_int(sampler["csv_device"], 0, MAX_UINT64, "sampler CSV device")
    exact_int(sampler["csv_inode"], 1, MAX_UINT64, "sampler CSV inode")
    start = exact_int(
        sampler["window_start_monotonic_ns"], 0, MAX_INT63,
        "sampler window start",
    )
    exact_int(
        sampler["window_end_monotonic_ns"], start, MAX_INT63,
        "sampler window end",
    )
    exact_absolute_path(sampler["script_path"], "sampler script path")
    exact_absolute_path(sampler["csv_path"], "sampler CSV path")
    lower_hash(sampler["script_sha256"], "sampler script SHA-256")


def validate_thermal_receipt(
    thermal: Any, sampler: Mapping[str, Any], child_start: int, child_reap: int,
) -> None:
    if type(thermal) is not dict:
        fail("health thermal receipt is not an object")
    exact_keys(thermal, THERMAL_RECEIPT_KEYS, "health thermal receipt")
    for name, lower, upper in (
        ("pid", 1, MAX_INT63), ("cpu", 0, MAX_INT63),
        ("process_start_ticks", 1, MAX_UINT64),
        ("csv_device", 0, MAX_UINT64), ("csv_inode", 1, MAX_UINT64),
        ("sample_count", 2, MAX_INT63),
        ("cpu_tctl_max_millic", 0, MAX_INT63),
        ("dimm_max_millic", 0, MAX_INT63),
        ("dimm_read_errors", 0, MAX_UINT64),
        ("edac_ce_max", 0, MAX_UINT64), ("edac_ue_max", 0, MAX_UINT64),
        ("window_start_monotonic_ns", 0, MAX_INT63),
        ("window_end_monotonic_ns", 0, MAX_INT63),
    ):
        exact_int(thermal[name], lower, upper, f"thermal {name}")
    exact_string(thermal["terminal_status"], "ok", "thermal status")
    exact_absolute_path(thermal["script_path"], "thermal script path")
    exact_absolute_path(thermal["csv_path"], "thermal CSV path")
    lower_hash(thermal["script_sha256"], "thermal script SHA-256")
    lower_hash(thermal["window_csv_sha256"], "thermal window SHA-256")
    for name in (
        "pid", "cpu", "process_start_ticks", "script_path", "script_sha256",
        "csv_path", "csv_device", "csv_inode", "terminal_status",
        "window_start_monotonic_ns", "window_end_monotonic_ns",
    ):
        if not canonical_equal(thermal[name], sampler[name]):
            fail(f"thermal/sampler identity differs at {name}")
    if not (
        thermal["window_start_monotonic_ns"] <= child_start
        <= child_reap <= thermal["window_end_monotonic_ns"]
    ):
        fail("thermal window does not cover the complete child interval")
    if (
        thermal["cpu_tctl_max_millic"] > THERMAL_MAX_MILLIC
        or thermal["dimm_max_millic"] > THERMAL_MAX_MILLIC
        or thermal["dimm_read_errors"] != 0
        or thermal["edac_ce_max"] != 0
        or thermal["edac_ue_max"] != 0
    ):
        fail("health thermal receipt violates the frozen policy")


def validate_health_receipt(health: Any, binary: Mapping[str, Any]) -> None:
    if type(health) is not dict:
        fail("summary health is not an object")
    exact_keys(health, HEALTH_KEYS, "summary health")
    exact_string(health["schema"], HEALTH_SCHEMA, "health schema")
    exact_string(health["terminal_status"], "ok", "health status")
    child_start = exact_int(
        health["child_start_monotonic_ns"], 0, MAX_INT63,
        "health child start",
    )
    child_reap = exact_int(
        health["child_reap_monotonic_ns"], child_start, MAX_INT63,
        "health child reap",
    )
    if (
        binary.get("child_start_monotonic_ns") != child_start
        or binary.get("child_reap_monotonic_ns") != child_reap
    ):
        fail("health child timestamps differ from the binary receipt")
    target_cpu = exact_int(
        health["target_cpu"], EXPECTED_TARGET_CPU, EXPECTED_TARGET_CPU,
        "health target CPU",
    )
    controller_cpu = exact_int(
        health["controller_cpu"], 0, MAX_INT63, "health controller CPU"
    )
    sampler_cpu = exact_int(
        health["sampler_cpu"], 0, MAX_INT63, "health sampler CPU"
    )
    target_core = exact_int_list(health["target_core"], "health target core", length=2)
    controller_core = exact_int_list(
        health["controller_core"], "health controller core", length=2
    )
    sampler_core = exact_int_list(
        health["sampler_core"], "health sampler core", length=2
    )
    target_threads = exact_int_list(
        health["target_threads"], "health target threads", sorted_unique=True
    )
    initial_affinity = exact_int_list(
        health["controller_initial_affinity"],
        "health controller initial affinity", sorted_unique=True,
    )
    singleton = exact_int_list(
        health["controller_singleton_affinity_end"],
        "health controller singleton affinity", length=1,
        sorted_unique=True,
    )
    if (
        tuple(target_core) != EXPECTED_TARGET_CORE
        or tuple(target_threads) != EXPECTED_TARGET_THREADS
        or singleton != [controller_cpu]
        or target_cpu not in initial_affinity
        or controller_cpu not in initial_affinity
        or len({tuple(target_core), tuple(controller_core), tuple(sampler_core)}) != 3
        or controller_cpu in target_threads
        or sampler_cpu in target_threads
    ):
        fail("health frozen topology/affinity differs")
    exact_string(
        health["edac_policy"], "ce-and-ue-every-sample-zero-v1",
        "health EDAC policy",
    )
    exact_string(
        health["sibling_tick_policy"],
        "linux-proc-stat-user-nice-system-irq-softirq-steal-v1",
        "health sibling tick policy",
    )
    exact_int(
        health["sibling_non_idle_tick_cap"], SIBLING_NON_IDLE_TICK_CAP,
        SIBLING_NON_IDLE_TICK_CAP, "health sibling cap",
    )
    exact_int(
        health["thermal_max_millic"], THERMAL_MAX_MILLIC,
        THERMAL_MAX_MILLIC, "health thermal cap",
    )
    sampler = health["sampler"]
    validate_sampler_receipt(sampler)
    if sampler["cpu"] != sampler_cpu:
        fail("sampler CPU differs from the health receipt")
    sampler_digest = lower_hash(
        health["sampler_receipt_sha256"], "sampler receipt SHA-256"
    )
    if sha256_bytes(canonical_bytes(sampler)) != sampler_digest:
        fail("health sampler receipt binding differs")
    if not (
        sampler["window_start_monotonic_ns"] <= child_start
        <= child_reap <= sampler["window_end_monotonic_ns"]
    ):
        fail("sampler endpoints do not cover the complete child interval")
    siblings = health["sibling_ticks"]
    if type(siblings) is not list or len(siblings) != 1:
        fail("health sibling tick roster differs")
    sibling = siblings[0]
    if type(sibling) is not dict:
        fail("health sibling tick receipt is not an object")
    exact_keys(
        sibling, {"cpu", "delta_non_idle_ticks", "end", "start"},
        "health sibling tick receipt",
    )
    exact_int(sibling["cpu"], 56, 56, "health sibling CPU")
    validate_cpu_tick_receipt(sibling["start"], "health sibling start")
    validate_cpu_tick_receipt(sibling["end"], "health sibling end")
    if sibling["start"]["cpu"] != 56 or sibling["end"]["cpu"] != 56:
        fail("health sibling tick endpoint CPU differs")
    recomputed_delta = sibling_tick_delta(
        sibling["start"], sibling["end"], child_start, child_reap
    )
    exact_int(
        sibling["delta_non_idle_ticks"], recomputed_delta, recomputed_delta,
        "health sibling delta",
    )
    validate_thermal_receipt(health["thermal"], sampler, child_start, child_reap)
    receipt_hash = lower_hash(health["receipt_sha256"], "health receipt SHA-256")
    health_preimage = dict(health)
    health_preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(health_preimage)) != receipt_hash:
        fail("health receipt SHA-256 differs")


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
        type(value["elapsed_seconds"]) not in (int, float)
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
    if value["health"]:
        if not binary:
            fail("health receipt lacks its binary receipt")
        validate_health_receipt(value["health"], binary)
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
        source_root: Path,
    ) -> None:
        self.output_dir = output_dir
        self.target_cpu = target_cpu
        self.expected_commit = expected_commit
        self.expected_source_manifest_sha256 = expected_source_manifest_sha256
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
            expected_source_manifest_sha256, source_root,
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

    def stage_evidence(
        self, raw: bytes, stderr: bytes, raw_complete: bool,
    ) -> List[str]:
        failures: List[str] = []
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
        if raw_complete and not failures:
            try:
                records = parse_transcript(self.staged_raw)
                (
                    self.raw_validation_config,
                    self.raw_validation_statistics,
                    self.raw_validation_failures,
                ) = validate_transcript(
                    records, self.target_cpu, self.expected_commit
                )
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

    def _bound_component(self, name: str, default: Any) -> Any:
        encoded = self.expected_summary_components.get(name)
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
            bounded_reason if self.expected_summary_components else None,
        )
        summary = make_summary_preimage({
            "binary": self._bound_component("binary", {}),
            "config": config,
            "elapsed_seconds": self._bound_component("elapsed_seconds", 0.0),
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
        return source_manifest(self.source_root)

    def _live_git_receipt(self) -> Dict[str, Any]:
        return git_receipt(self.source_root, self.expected_commit)

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


def sibling_tick_delta(
    start: Mapping[str, Any], end: Mapping[str, Any],
    child_start_ns: int, child_reap_ns: int,
) -> int:
    validate_cpu_tick_receipt(start, "sibling start")
    validate_cpu_tick_receipt(end, "sibling end")
    exact_int(child_start_ns, 0, MAX_INT63, "child start timestamp")
    exact_int(child_reap_ns, child_start_ns, MAX_INT63, "child reap timestamp")
    if not (
        start["read_monotonic_ns"] <= child_start_ns
        <= child_reap_ns <= end["read_monotonic_ns"]
    ):
        fail("sibling tick endpoints do not cover the complete child interval")
    if start["cpu"] != end["cpu"]:
        fail("SMT sibling tick CPU changed")
    start_fields = start["tick_fields"]
    end_fields = end["tick_fields"]
    if set(start_fields) != set(end_fields) or any(
        end_fields[field] < start_fields[field] for field in start_fields
    ):
        fail("SMT sibling /proc/stat fields are not monotone")
    delta = end["non_idle_ticks"] - start["non_idle_ticks"]
    if delta < 0 or delta > SIBLING_NON_IDLE_TICK_CAP:
        fail("SMT sibling exceeded the frozen non-idle tick cap")
    return delta


def make_sampler_attestation(
    args: argparse.Namespace, start_ns: int, end_ns: int,
    health_modules: SealedHealthModules,
) -> Dict[str, Any]:
    health_api = health_modules.native
    process_start_ticks = health_api._parse_proc_start_ticks(args.sampler_pid)
    health_api._require_process_predates_window(
        process_start_ticks, start_ns, "broad-screen sampler"
    )
    health_api._verify_live_sampler_process(
        args.sampler_pid, args.sampler_cpu, process_start_ticks,
        args.sampler_script, args.sampler_csv,
    )
    try:
        resolved_script = str(args.sampler_script.resolve(strict=True))
        resolved_csv = str(args.sampler_csv.resolve(strict=True))
    except OSError as exc:
        fail(f"cannot resolve sampler artifacts: {exc}")
    script_fd = -1
    csv_fd = -1
    try:
        script_fd, _ = health_api._open_regular_nofollow(
            args.sampler_script, "broad-screen sampler script"
        )
        csv_fd, csv_info = health_api._open_regular_nofollow(
            args.sampler_csv, "broad-screen sampler CSV"
        )
        result = {
            "cpu": args.sampler_cpu,
            "csv_device": csv_info.st_dev,
            "csv_inode": csv_info.st_ino,
            "csv_path": resolved_csv,
            "pid": args.sampler_pid,
            "process_start_ticks": process_start_ticks,
            "schema": health_api.SAMPLER_SCHEMA,
            "script_path": resolved_script,
            "script_sha256": health_api._sha256_descriptor(script_fd),
            "terminal_status": "ok",
            "window_end_monotonic_ns": end_ns,
            "window_start_monotonic_ns": start_ns,
        }
    finally:
        if csv_fd >= 0:
            os.close(csv_fd)
        if script_fd >= 0:
            os.close(script_fd)
    health_api._verify_live_sampler_process(
        args.sampler_pid, args.sampler_cpu, process_start_ticks,
        args.sampler_script, args.sampler_csv,
    )
    return result


def prepare_health(
    args: argparse.Namespace, target_cpu: int, deadline: float,
    health_modules: SealedHealthModules,
) -> Dict[str, Any]:
    health_api = health_modules.native
    health_runner = health_modules.runner
    if target_cpu != EXPECTED_TARGET_CPU:
        fail(f"frozen health policy requires target CPU {EXPECTED_TARGET_CPU}")
    target_core = health_api._cpu_physical_core(target_cpu)
    siblings = tuple(health_api.timing_sibling_cpus([target_cpu]))
    threads = tuple(sorted((target_cpu, *siblings)))
    if target_core != EXPECTED_TARGET_CORE or threads != EXPECTED_TARGET_THREADS:
        fail("target CPU physical-core topology differs from the frozen pair")
    controller_core = health_api._cpu_physical_core(args.controller_cpu)
    sampler_core = health_api._cpu_physical_core(args.sampler_cpu)
    if (
        args.controller_cpu in threads or args.sampler_cpu in threads
        or controller_core == target_core or sampler_core == target_core
        or controller_core == sampler_core
    ):
        fail("target, controller, and sampler physical cores are not distinct")
    health_runner._preflight_sampler(
        args.sampler_pid, args.sampler_cpu,
        args.sampler_script, args.sampler_csv,
    )
    try:
        allowed = os.sched_getaffinity(0)
    except (AttributeError, OSError) as exc:
        fail(f"cannot inspect controller affinity: {exc}")
    if args.controller_cpu not in allowed or target_cpu not in allowed:
        fail("target/controller CPUs are outside the runner's initial affinity")
    try:
        os.sched_setaffinity(0, {args.controller_cpu})
    except (AttributeError, OSError) as exc:
        fail(f"cannot pin the controller CPU: {exc}")
    if os.sched_getaffinity(0) != {args.controller_cpu}:
        fail("controller did not retain singleton affinity")

    sample_start_ns = health_runner.choose_new_sampler_start(
        args.sampler_csv, deadline
    )
    if sample_start_ns > time.monotonic_ns():
        fail("sampler start timestamp is in the future")
    tick_start = read_cpu_tick_receipt(siblings[0])
    return {
        "controller_core": controller_core,
        "initial_affinity": sorted(allowed),
        "sample_start_ns": sample_start_ns,
        "sampler_core": sampler_core,
        "sibling_cpus": list(siblings),
        "target_core": target_core,
        "target_threads": list(threads),
        "tick_start": tick_start,
    }


def finish_health(
    args: argparse.Namespace, target_cpu: int, child_start_ns: int,
    child_reap_ns: int, state: Mapping[str, Any], deadline: float,
    health_modules: SealedHealthModules,
) -> Dict[str, Any]:
    health_api = health_modules.native
    health_runner = health_modules.runner
    try:
        controller_affinity_end = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError) as exc:
        fail(f"cannot recheck controller affinity: {exc}")
    if controller_affinity_end != [args.controller_cpu]:
        fail("controller singleton affinity changed during the child interval")
    target_core_end = health_api._cpu_physical_core(target_cpu)
    sibling_cpus_end = tuple(health_api.timing_sibling_cpus([target_cpu]))
    target_threads_end = tuple(sorted((target_cpu, *sibling_cpus_end)))
    controller_core_end = health_api._cpu_physical_core(args.controller_cpu)
    sampler_core_end = health_api._cpu_physical_core(args.sampler_cpu)
    if (
        target_core_end != state["target_core"]
        or sibling_cpus_end != tuple(state["sibling_cpus"])
        or target_threads_end != tuple(state["target_threads"])
        or controller_core_end != state["controller_core"]
        or sampler_core_end != state["sampler_core"]
    ):
        fail("target/controller/sampler topology changed during the child interval")
    tick_end = read_cpu_tick_receipt(state["sibling_cpus"][0])
    tick_delta = sibling_tick_delta(
        state["tick_start"], tick_end, child_start_ns, child_reap_ns
    )
    sample_end_ns = health_runner._wait_for_sampler_sample(
        args.sampler_csv, deadline, at_or_after_ns=child_reap_ns
    )
    sampler = make_sampler_attestation(
        args, state["sample_start_ns"], sample_end_ns, health_modules
    )
    try:
        thermal = health_api._thermal_window(
            sampler, child_start_ns, child_reap_ns, [target_cpu], True,
            controller_cpu=args.controller_cpu,
        )
    except health_api.NativeEvidenceError as exc:
        fail(str(exc))
    if (
        thermal["sample_count"] < 2
        or thermal["cpu_tctl_max_millic"] > THERMAL_MAX_MILLIC
        or thermal["dimm_max_millic"] > THERMAL_MAX_MILLIC
        or thermal["dimm_read_errors"] != 0
        or thermal["edac_ce_max"] != 0
        or thermal["edac_ue_max"] != 0
    ):
        fail("thermal/EDAC receipt violates the frozen health cap")
    if not (
        sampler["window_start_monotonic_ns"] <= child_start_ns
        and sampler["window_end_monotonic_ns"] >= child_reap_ns
    ):
        fail("sampler endpoints do not cover the complete child interval")
    receipt: Dict[str, Any] = {
        "child_reap_monotonic_ns": child_reap_ns,
        "child_start_monotonic_ns": child_start_ns,
        "controller_core": list(state["controller_core"]),
        "controller_cpu": args.controller_cpu,
        "controller_initial_affinity": list(state["initial_affinity"]),
        "controller_singleton_affinity_end": controller_affinity_end,
        "edac_policy": "ce-and-ue-every-sample-zero-v1",
        "receipt_sha256": None,
        "sampler": sampler,
        "sampler_core": list(state["sampler_core"]),
        "sampler_cpu": args.sampler_cpu,
        "sampler_receipt_sha256": sha256_bytes(canonical_bytes(sampler)),
        "schema": HEALTH_SCHEMA,
        "sibling_non_idle_tick_cap": SIBLING_NON_IDLE_TICK_CAP,
        "sibling_tick_policy": (
            "linux-proc-stat-user-nice-system-irq-softirq-steal-v1"
        ),
        "sibling_ticks": [{
            "cpu": state["sibling_cpus"][0], "delta_non_idle_ticks": tick_delta,
            "end": tick_end, "start": state["tick_start"],
        }],
        "target_core": list(state["target_core"]),
        "target_cpu": target_cpu,
        "target_threads": list(state["target_threads"]),
        "thermal": dict(thermal),
        "thermal_max_millic": THERMAL_MAX_MILLIC,
        "terminal_status": "ok",
    }
    receipt["receipt_sha256"] = sha256_bytes(canonical_bytes(receipt))
    return receipt


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
        "child_reap_monotonic_ns": None,
        "child_start_monotonic_ns": None,
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


def run_binary(
    binary: Path, cpu: int, expected_sha256: str,
    health_args: Optional[argparse.Namespace] = None,
    health_modules: Optional[SealedHealthModules] = None,
    signal_guard: Optional[SignalGuard] = None,
) -> ExecutionResult:
    result = ExecutionResult(binary=empty_binary_receipt(binary, expected_sha256))
    started = time.monotonic()
    deadline = started + OUTER_DEADLINE_SECONDS
    fd = -1
    process: Optional[subprocess.Popen] = None
    fd_before: Optional[os.stat_result] = None
    digest_before: Optional[str] = None
    health_state: Dict[str, Any] = {}
    capture_started = False
    try:
        path_before = os.stat(str(binary), follow_symlinks=False)
        if not stat.S_ISREG(path_before.st_mode) or path_before.st_nlink != 1:
            fail("screen binary must be a regular single-link file")
        if stat.S_IMODE(path_before.st_mode) != 0o555:
            fail("screen binary must have exact mode 0555")
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if not nofollow:
            fail("screen binary execution requires O_NOFOLLOW")
        fd = os.open(
            str(binary), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | nofollow
        )
        fd_before = os.fstat(fd)
        result.binary["stat_before"] = stat_receipt(fd_before)
        if not same_file_receipt(path_before, fd_before):
            fail("screen binary path/FD identity differs")
        digest_before = file_sha256_fd(fd, fd_before.st_size)
        result.binary["sha256_before"] = digest_before
        if digest_before != expected_sha256:
            fail("screen binary SHA-256 differs from the presealed value")
        if health_args is not None:
            if health_modules is None:
                fail("health execution lacks sealed helper modules")
            health_state = prepare_health(
                health_args, cpu, deadline, health_modules
            )
        if signal_guard is not None and signal_guard.first_signal is not None:
            fail(f"runner received {signal_guard.first_signal} before child start")
        fd_path = f"/proc/self/fd/{fd}"
        child_start_ns = time.monotonic_ns()
        process = subprocess.Popen(
            [fd_path, "--run", "--cpu", str(cpu)], executable=fd_path,
            pass_fds=(fd,), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env={"LANG": "C", "LC_ALL": "C"},
            start_new_session=True,
        )
        result.binary["process_started"] = True
        result.binary["child_start_monotonic_ns"] = child_start_ns
        capture = bounded_capture(process, deadline, signal_guard)
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
        if capture.error != "none":
            result.errors.append(capture.error)
        if health_args is not None and child_reap_ns is not None:
            if health_modules is None:
                fail("health finish lacks sealed helper modules")
            try:
                result.health = finish_health(
                    health_args, cpu, child_start_ns, child_reap_ns,
                    health_state, deadline, health_modules,
                )
            except BaseException as exc:
                result.errors.append(f"health finish: {exception_text(exc)}")
        elif health_args is not None:
            result.errors.append("health finish: child was not reaped")
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
            if process.returncode is not None:
                result.binary["child_reap_monotonic_ns"] = time.monotonic_ns()
            if signal_guard is not None:
                signal_guard.detach(process)
        if fd >= 0:
            try:
                fd_after = os.fstat(fd)
                result.binary["stat_after"] = stat_receipt(fd_after)
                digest_after = file_sha256_fd(fd, fd_after.st_size)
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
        result.elapsed_seconds = time.monotonic() - started
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


def git_receipt(root: Path, expected_commit: str) -> Dict[str, Any]:
    environment = {
        "GIT_OPTIONAL_LOCKS": "0", "LANG": "C", "LC_ALL": "C",
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
    }

    def git(*arguments: str) -> bytes:
        completed = subprocess.run(
            ["git", "-c", "core.fsmonitor=false", *arguments], cwd=root,
            env=environment, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=10.0, check=False,
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
    expected_paths = b"".join((relative + "\n").encode("ascii") for relative in SOURCE_PATHS)
    if tracked_paths != expected_paths:
        fail("screen source roster is not tracked in the expected order")
    index_flags = git(
        "ls-files", "-v", "--error-unmatch", "--", *SOURCE_PATHS
    )
    validate_git_index_flag_roster(index_flags)
    git("diff", "--quiet", "HEAD", "--", *SOURCE_PATHS)
    return {
        "head": expected_commit, "source_roster_sha256": sha256_bytes(tracked_paths),
        "tracked_status_sha256": sha256_bytes(tracked_status),
    }


def run_once(args: argparse.Namespace) -> int:
    guarded_signals = SignalGuard()
    with guarded_signals:
        output_dir = Path(args.output_dir)
        binary_path = Path(args.binary)
        execution = ExecutionResult(
            binary=empty_binary_receipt(binary_path, args.expected_binary_sha256)
        )
        source_before: Dict[str, Any] = {}
        source_after: Dict[str, Any] = {}
        git_before: Dict[str, Any] = {}
        git_after: Dict[str, Any] = {}
        statistics: Dict[str, Any] = {}
        config: Dict[str, Any] = {}
        health_modules: Optional[SealedHealthModules] = None
        health_module_loader: Dict[str, Any] = {}
        infrastructure_errors: List[str] = []
        gate_failures: List[str] = []
        outcome = "invalid"
        bundle: Optional[OutputBundle] = None
        stage_failures: List[str] = []
        finish_result = BundleFinishResult()
        source_root = Path(__file__).resolve().parents[1]

        # A signal that predates any reservation mutation consumes no output
        # path.  A signal after this check is latched into the owned bundle.
        if guarded_signals.first_signal is not None:
            return 1
        bundle = OutputBundle.reserve(
            output_dir, args.cpu, args.expected_source_commit,
            args.expected_source_manifest_sha256, source_root,
            guarded_signals,
        )
        try:
            try:
                git_before = git_receipt(
                    source_root, args.expected_source_commit
                )
                source_before = source_manifest(source_root)
                if source_before["sha256"] != args.expected_source_manifest_sha256:
                    fail("source manifest differs from the presealed value")
                health_modules = load_sealed_health_modules(
                    source_root, source_before
                )
                health_module_loader = health_modules.receipt
            except BaseException as exc:
                infrastructure_errors.append(
                    f"source preflight: {exception_text(exc)}"
                )

            if not infrastructure_errors and guarded_signals.first_signal is None:
                execution = run_binary(
                    binary_path, args.cpu, args.expected_binary_sha256,
                    health_args=args, health_modules=health_modules,
                    signal_guard=guarded_signals,
                )
                infrastructure_errors.extend(execution.errors)

            # Preserve post-run source evidence even when health or binary
            # adjudication failed after the child was already captured.
            try:
                source_after = source_manifest(source_root)
                git_after = git_receipt(
                    source_root, args.expected_source_commit
                )
                if source_before and source_before != source_after:
                    fail("screen source changed during execution")
            except BaseException as exc:
                infrastructure_errors.append(
                    f"source postcheck: {exception_text(exc)}"
                )

            if execution.binary.get("timed_out"):
                infrastructure_errors.append(
                    "screen exceeded the hard 120-second deadline"
                )
            if execution.binary.get("output_overflow"):
                infrastructure_errors.append(
                    "screen exceeded its bounded output allowance"
                )
            if execution.elapsed_seconds >= OUTER_DEADLINE_SECONDS:
                infrastructure_errors.append(
                    "screen returned at or after the hard 120-second deadline"
                )
            if execution.binary.get("process_started") and execution.returncode != 0:
                infrastructure_errors.append(
                    f"screen exited {execution.returncode}"
                )
            if execution.stderr:
                infrastructure_errors.append("screen emitted stderr")
            if guarded_signals.kill_error is not None:
                infrastructure_errors.append(
                    f"signal process-group kill: {guarded_signals.kill_error}"
                )
            if guarded_signals.first_signal is not None:
                infrastructure_errors.append(
                    f"runner received {guarded_signals.first_signal}"
                )

            if execution.raw_complete():
                try:
                    records = parse_transcript(execution.raw)
                    config, statistics, gate_failures = validate_transcript(
                        records, args.cpu, args.expected_source_commit
                    )
                except BaseException as exc:
                    infrastructure_errors.append(
                        f"transcript validation: {exception_text(exc)}"
                    )
            if infrastructure_errors:
                outcome = "invalid"
            elif config and statistics:
                outcome = "pass" if not gate_failures else "reject"
        except BaseException as exc:
            infrastructure_errors.append(
                f"runner body: {exception_text(exc)}"
            )
            outcome = "invalid"
        finally:
            # This nested owner guard spans every path after reservation,
            # including BaseException and summary-construction failures.
            try:
                stage_failures = bundle.stage_evidence(
                    execution.raw, execution.stderr, execution.raw_complete()
                )
            except BaseException as exc:
                stage_failures = [f"evidence staging: {exception_text(exc)}"]
            if stage_failures:
                infrastructure_errors.extend(stage_failures)
                outcome = "invalid"
            if (
                bundle.raw_validation_config
                and bundle.raw_validation_statistics
            ):
                config = bundle.raw_validation_config
                statistics = bundle.raw_validation_statistics
                gate_failures = list(bundle.raw_validation_failures)
                if not infrastructure_errors:
                    outcome = "pass" if not gate_failures else "reject"
            if guarded_signals.kill_error is not None and not any(
                "signal process-group kill" in item
                for item in infrastructure_errors
            ):
                infrastructure_errors.append(
                    f"signal process-group kill: {guarded_signals.kill_error}"
                )
            if guarded_signals.first_signal is not None and not any(
                guarded_signals.first_signal in item
                for item in infrastructure_errors
            ):
                infrastructure_errors.append(
                    f"runner received {guarded_signals.first_signal}"
                )
            if infrastructure_errors:
                outcome = "invalid"
            elif outcome not in ("pass", "reject"):
                infrastructure_errors.append("runner did not reach validation")
                outcome = "invalid"

            try:
                snapshot_signal = guarded_signals.first_signal
                snapshot_signal_names = list(guarded_signals.observed_signals)
                if snapshot_signal is not None and not any(
                    snapshot_signal in item for item in infrastructure_errors
                ):
                    infrastructure_errors.append(
                        f"runner received {snapshot_signal}"
                    )
                    outcome = "invalid"
                bundle.signal_name = snapshot_signal
                bundle.signal_names = snapshot_signal_names
                failure = (
                    bounded_failure_text(infrastructure_errors)
                    if infrastructure_errors
                    else ("statistical-gates" if gate_failures else "none")
                )
                raw_complete = execution.raw_complete() and not stage_failures
                bundle.bind_expected_summary_components({
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
                    "target_cpu": args.cpu,
                })
                summary = make_summary_preimage({
                    "binary": execution.binary,
                    "config": config,
                    "elapsed_seconds": execution.elapsed_seconds,
                    "expected_source_git_commit": args.expected_source_commit,
                    "expected_source_manifest_sha256": (
                        args.expected_source_manifest_sha256
                    ),
                    "failure": failure,
                    "git_after": git_after,
                    "git_before": git_before,
                    "health": execution.health,
                    "health_module_loader": health_module_loader,
                    "outcome": outcome,
                    "output_bundle": bundle.receipt(),
                    "process_exit_code": execution.returncode,
                    "publication_failures": list(stage_failures),
                    "raw_bytes": len(bundle.staged_raw),
                    "raw_complete": raw_complete,
                    "raw_record_count": bundle.staged_raw.count(b"\n"),
                    "raw_sha256": sha256_bytes(bundle.staged_raw),
                    "schema": SUMMARY_SCHEMA,
                    "signal": snapshot_signal,
                    "signal_names": snapshot_signal_names,
                    "source_manifest_after": source_after,
                    "source_manifest_before": source_before,
                    "statistics": statistics,
                    "stderr_bytes": len(bundle.staged_stderr),
                    "stderr_sha256": sha256_bytes(bundle.staged_stderr),
                    "summary_preimage_sha256": None,
                    "target_cpu": args.cpu,
                })
                summary_bytes = canonical_bytes(summary) + b"\n"
            except BaseException as exc:
                outcome = "invalid"
                bundle.signal_name = guarded_signals.first_signal
                infrastructure_errors.append(
                    f"summary construction: {exception_text(exc)}"
                )
                bundle.signal_names = list(guarded_signals.observed_signals)
                summary_bytes = bundle.emergency_summary(
                    bounded_failure_text(infrastructure_errors),
                    stage_failures,
                )
            try:
                finish_result = bundle.finish(summary_bytes, guarded_signals)
            except BaseException as exc:
                failure = f"bundle finish: {exception_text(exc)}"
                finish_result = BundleFinishResult(
                    logical_commit_succeeded=bundle.logically_sealed,
                    precommit_failures=(
                        tuple() if bundle.logically_sealed else (failure,)
                    ),
                    postcommit_failures=(
                        (failure,) if bundle.logically_sealed else tuple()
                    ),
                )
            if (
                not finish_result.logical_commit_succeeded
                or finish_result.precommit_failures
                or guarded_signals.first_signal is not None
            ):
                outcome = "invalid"
        if not bundle.closed:
            # `finish` is designed not to raise, but this second idempotent
            # ownership drain protects future refactors and injected failures.
            bundle.signal_name = guarded_signals.first_signal
            bundle.signal_names = list(guarded_signals.observed_signals)
            emergency = bundle.emergency_summary(
                "outer output ownership guard", finish_result.all_failures()
            )
            try:
                drain_result = bundle.finish(emergency, guarded_signals)
            except BaseException as exc:
                failure = f"outer bundle finish: {exception_text(exc)}"
                drain_result = BundleFinishResult(
                    logical_commit_succeeded=bundle.logically_sealed,
                    precommit_failures=(
                        tuple() if bundle.logically_sealed else (failure,)
                    ),
                    postcommit_failures=(
                        (failure,) if bundle.logically_sealed else tuple()
                    ),
                )
            finish_result = finish_result.merged(drain_result)

    # A post-commit durability failure is a distinct publication failure.  It
    # fails the controller exit status without pretending the immutable summary
    # decision was rewritten from pass/reject to invalid after logical commit.
    if (
        not finish_result.logical_commit_succeeded
        or finish_result.precommit_failures
        or finish_result.postcommit_failures
        or outcome == "invalid"
    ):
        return 1
    return 0 if outcome == "pass" else 2


def synthetic_config(cpu: int, commit: str) -> Dict[str, Any]:
    cells = []
    for k, block_bytes in CELLS:
        key = (k, block_bytes)
        configuration = {
            "K": k, "block_bytes": block_bytes, "dense_anchor_layout": 0,
            "dense_identity_corner": False, "dense_rows": 12,
            "heavy_family": 0, "heavy_rows": 12, "mix_count": 3,
            "packet_attempt": 0, "packet_peel_seed": 1,
            "precode_attempt": 0, "precode_seed": (1 << 63) + k,
            "source_hits": 2,
            "staircase": 1,
        }
        # Selftests run after real seals have replaced the placeholders.  Use
        # the exact digest map while retaining a syntactically valid object by
        # allowing a local override below in synthetic_transcript().
        cell = {
            "K": k, "block_bytes": block_bytes,
            "construction_equivalent": True,
            "equation_configuration": configuration,
            "equation_configuration_sha256": sha256_bytes(canonical_bytes(configuration)),
            "first_repair_sha256": EXPECTED_FIRST_REPAIR_SHA256[key],
            "high_id_oracle_sha256": EXPECTED_HIGH_ID_ORACLE_SHA256[key],
            "high_repair_controls_verified": True,
            "repair_oracle_sha256": EXPECTED_REPAIR_ORACLE_SHA256[key],
            "solved_state_equivalence_receipt_sha256": (
                EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256[key]
            ),
            "source_sha256": EXPECTED_SOURCE_SHA256[key],
        }
        cells.append(cell)
    descriptors = []
    for metric, _, _ in METRICS:
        for mode, emission in (
            ("equation", "equation-eval-v1"),
            ("retained", "retained-source-direct-v1"),
        ):
            value = descriptor(metric, emission)
            descriptors.append({
                "descriptor": value, "descriptor_sha256": descriptor_hash(value),
                "mode": mode,
            })
    return {
        "cells": cells, "descriptors": descriptors,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "minimum_invocations": MINIMUM_INVOCATIONS,
        "panel_replicates": REPLICATES, "schema": CONFIG_SCHEMA,
        "source_git_commit": commit, "target_cpu": cpu,
    }


def synthetic_transcript(cpu: int, commit: str) -> bytes:
    config = synthetic_config(cpu, commit)
    # Synthetic configuration objects need hashes that pass the exact sealed
    # maps.  The selftest substitutes those maps temporarily.
    records: List[Dict[str, Any]] = [config]
    hashes: Dict[Tuple[str, str], str] = {}
    for entry in config["descriptors"]:
        metric = entry["descriptor"]["metric"]
        hashes[(metric, entry["mode"])] = entry["descriptor_sha256"]
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in PANEL_SCOPES:
                for comparison in COMPARISONS:
                    order = panel_order(k, block_bytes, comparison, scope, replicate)
                    sides = expected_sides(order)
                    left_mode, right_mode = comparison_modes(comparison)
                    primary = next(
                        metric for metric, metric_scope_name, _ in METRICS
                        if metric_scope_name == scope
                        and metric != "fresh-repair-init-first"
                    )
                    left_time = 1000
                    right_time = 1000
                    if (
                        primary == "fresh-systematic-total"
                        and comparison == "candidate-over-baseline"
                    ):
                        left_time = 900
                    records.append({
                        "K": k, "block_bytes": block_bytes,
                        "comparison": comparison,
                        "invocations_per_slot": invocations_for(k, replicate),
                        "left_descriptor_sha256": hashes[(primary, left_mode)],
                        "left_direct_systematic_packets": (
                            k if scope == "fresh-systematic" and left_mode == "retained" else 0
                        ),
                        "left_nested_descriptor_sha256": (
                            hashes[("fresh-repair-init-first", left_mode)]
                            if scope == "fresh-repair" else None
                        ),
                        "left_outcome_code": 0,
                        "nested_metric": (
                            "fresh-repair-init-first" if scope == "fresh-repair" else "none"
                        ),
                        "order": order,
                        "panel_key_sha256": panel_key_sha256(
                            k, block_bytes, comparison, scope
                        ),
                        "primary_metric": primary, "replicate": replicate,
                        "right_descriptor_sha256": hashes[(primary, right_mode)],
                        "right_direct_systematic_packets": (
                            k if scope == "fresh-systematic" and right_mode == "retained" else 0
                        ),
                        "right_nested_descriptor_sha256": (
                            hashes[("fresh-repair-init-first", right_mode)]
                            if scope == "fresh-repair" else None
                        ),
                        "right_outcome_code": 0, "schema": PANEL_SCHEMA,
                        "scope": scope,
                        "slots": [{
                            "elapsed_ns": left_time if side == "left" else right_time,
                            "nested_elapsed_ns": 800 if scope == "fresh-repair" else 0,
                            "side": side,
                        } for side in sides],
                        "target_cpu": cpu,
                    })
    records.append({
        "panel_count": EXPECTED_PANEL_COUNT, "schema": TERMINAL_SCHEMA,
        "status": "complete",
    })
    return b"".join(canonical_bytes(record) + b"\n" for record in records)


def selftest() -> int:
    cpu = 7
    commit = "1" * 40
    project_root = Path(__file__).resolve().parents[1]
    loader_manifest = source_manifest(project_root)
    loaded_health = load_sealed_health_modules(project_root, loader_manifest)
    validate_health_loader_receipt(loaded_health.receipt)
    loaded_entries = {
        entry["path"]: entry for entry in loader_manifest["entries"]
    }
    if any(
        module["bytes"] != loaded_entries[module["path"]]["bytes"]
        or module["sha256"] != loaded_entries[module["path"]]["sha256"]
        for module in loaded_health.receipt["modules"]
    ) or (
        loaded_health.native.contract_api is not loaded_health.contract
        or loaded_health.runner.contract_api is not loaded_health.contract
        or loaded_health.runner.native_api is not loaded_health.native
    ):
        raise AssertionError("sealed health source-loader receipt differs")
    # The native selftest owns the real frozen hashes.  Build internally
    # coherent synthetic seals to exercise every protocol/statistical law.
    global EXPECTED_CONFIGURATION_SHA256, EXPECTED_SOURCE_SHA256
    global EXPECTED_REPAIR_ORACLE_SHA256, EXPECTED_FIRST_REPAIR_SHA256
    global EXPECTED_HIGH_ID_ORACLE_SHA256
    global EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256
    synthetic = synthetic_config(cpu, commit)
    EXPECTED_CONFIGURATION_SHA256 = {
        (cell["K"], cell["block_bytes"]): cell["equation_configuration_sha256"]
        for cell in synthetic["cells"]
    }
    EXPECTED_SOURCE_SHA256 = {
        cell: sha256_bytes(f"source:{cell[0]}:{cell[1]}".encode("ascii"))
        for cell in CELLS
    }
    EXPECTED_REPAIR_ORACLE_SHA256 = {
        cell: sha256_bytes(f"repair:{cell[0]}:{cell[1]}".encode("ascii"))
        for cell in CELLS
    }
    EXPECTED_FIRST_REPAIR_SHA256 = {
        cell: sha256_bytes(f"first:{cell[0]}:{cell[1]}".encode("ascii"))
        for cell in CELLS
    }
    EXPECTED_HIGH_ID_ORACLE_SHA256 = {
        cell: sha256_bytes(f"high:{cell[0]}:{cell[1]}".encode("ascii"))
        for cell in CELLS
    }
    EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256 = {}
    for cell in CELLS:
        receipt = {
            "configuration_sha256": EXPECTED_CONFIGURATION_SHA256[cell],
            "first_repair_sha256": EXPECTED_FIRST_REPAIR_SHA256[cell],
            "high_id_oracle_sha256": EXPECTED_HIGH_ID_ORACLE_SHA256[cell],
            "repair_oracle_sha256": EXPECTED_REPAIR_ORACLE_SHA256[cell],
            "source_sha256": EXPECTED_SOURCE_SHA256[cell],
            "systematic_oracle_sha256": EXPECTED_SOURCE_SHA256[cell],
        }
        EXPECTED_SOLVED_STATE_EQUIVALENCE_RECEIPT_SHA256[cell] = sha256_bytes(
            canonical_bytes(receipt)
        )

    raw = synthetic_transcript(cpu, commit)
    records = parse_transcript(raw)
    _, statistics, failures = validate_transcript(records, cpu, commit)
    if failures:
        raise AssertionError(f"passing synthetic transcript failed: {failures}")
    if not statistics["aggregates"]["fresh-systematic-total"]["pooled"]["upper95"] < 0.99:
        raise AssertionError("synthetic systematic improvement was not measured")
    if not records[0]["cells"][0]["equation_configuration"]["precode_seed"] > MAX_INT63:
        raise AssertionError("uint64 precode-seed selftest lacks >int63 coverage")

    hashes = validate_config(records[0], cpu, commit)
    probe = json.loads(json.dumps(records[1]))
    elapsed = [1, 2, 3, 4, 5, 6, 7, 8]
    for slot, value in zip(probe["slots"], elapsed):
        slot["elapsed_ns"] = value
    lanes = validate_panel(
        probe, cpu, probe["replicate"], probe["K"], probe["block_bytes"],
        probe["scope"], probe["comparison"], hashes,
    )
    if not math.isclose(
        lanes[probe["primary_metric"]], lane_contrast(elapsed, probe["order"]),
        rel_tol=0.0, abs_tol=1e-15,
    ):
        raise AssertionError("nonuniform lane contrast differs")

    mutations: List[Tuple[str, List[Dict[str, Any]]]] = []
    bad_count = json.loads(json.dumps(records))
    bad_count[-1]["panel_count"] -= 1
    mutations.append(("terminal count", bad_count))
    bad_order = json.loads(json.dumps(records))
    bad_order[1]["order"] = "BAAB" if records[1]["order"] == "ABBA" else "ABBA"
    mutations.append(("panel order", bad_order))
    bad_nested = json.loads(json.dumps(records))
    fresh_index = next(
        index for index, record in enumerate(records)
        if index > 0 and record.get("scope") == "fresh-repair"
    )
    bad_nested[fresh_index]["slots"][0]["nested_elapsed_ns"] = 0
    mutations.append(("nested zero", bad_nested))
    bad_hash = json.loads(json.dumps(records))
    bad_hash[0]["cells"][0][
        "solved_state_equivalence_receipt_sha256"
    ] = "0" * 64
    mutations.append(("solved-state hash", bad_hash))
    for label, mutated in mutations:
        try:
            validate_transcript(mutated, cpu, commit)
        except ValidationError:
            pass
        else:
            raise AssertionError(f"tampered {label} was accepted")

    nonimproving = json.loads(json.dumps(records))
    for record in nonimproving[1:-1]:
        if (
            record["scope"] == "fresh-systematic"
            and record["comparison"] == "candidate-over-baseline"
        ):
            for slot in record["slots"]:
                slot["elapsed_ns"] = 1000
    _, _, failed = validate_transcript(nonimproving, cpu, commit)
    if "aggregate_upper95:fresh-systematic-total:pooled" not in failed:
        raise AssertionError("non-improving systematic candidate passed")

    slow_repair = json.loads(json.dumps(records))
    for record in slow_repair[1:-1]:
        if (
            record["scope"] == "steady-repair"
            and record["comparison"] == "candidate-over-baseline"
            and record["K"] == 8 and record["block_bytes"] == 64
        ):
            for slot in record["slots"]:
                slot["elapsed_ns"] = 1030 if slot["side"] == "left" else 1000
    _, _, failed = validate_transcript(slow_repair, cpu, commit)
    expected_failure = (
        "cell_point:steady-repair-total:candidate-over-baseline:K8:B64"
    )
    if expected_failure not in failed:
        raise AssertionError("slow repair cell passed its point gate")

    aa_bias = json.loads(json.dumps(records))
    for record in aa_bias[1:-1]:
        if (
            record["scope"] == "fresh-repair"
            and record["comparison"] == "candidate-aa"
            and record["K"] == 128 and record["block_bytes"] == 1280
        ):
            for slot in record["slots"]:
                slot["nested_elapsed_ns"] = (
                    824 if slot["side"] == "left" else 800
                )
    _, _, failed = validate_transcript(aa_bias, cpu, commit)
    expected_failure = (
        "aa_ci:fresh-repair-init-first:candidate-aa:K128:B1280"
    )
    if expected_failure not in failed:
        raise AssertionError("biased nested A/A cell passed")

    duplicate = b'{"cells":[],"cells":[]}\n' + raw.split(b"\n", 1)[1]
    try:
        parse_transcript(duplicate)
    except ValidationError:
        pass
    else:
        raise AssertionError("duplicate JSON key was accepted")
    missing = raw.rsplit(b"\n", 2)[0] + b"\n"
    try:
        parse_transcript(missing)
    except ValidationError:
        pass
    else:
        raise AssertionError("missing terminal record was accepted")

    tick_start = parse_cpu_tick_receipt(
        b"cpu55 1 1 1 1 1 1 1 1 1 1\n"
        b"cpu56 10 2 3 100 4 5 6 7 999 888\n",
        56, 10,
    )
    tick_end = parse_cpu_tick_receipt(
        b"cpu56 15 2 3 101 4 5 6 7 5000 4000\n",
        56, 40,
    )
    if tick_start["non_idle_ticks"] != 33:
        raise AssertionError("guest ticks entered the sibling non-idle sum")
    if sibling_tick_delta(tick_start, tick_end, 20, 30) != 5:
        raise AssertionError("sibling tick boundary calculation differs")
    too_busy = parse_cpu_tick_receipt(
        b"cpu56 16 2 3 101 4 5 6 7 5000 4000\n", 56, 40
    )
    try:
        sibling_tick_delta(tick_start, too_busy, 20, 30)
    except ValidationError:
        pass
    else:
        raise AssertionError("sibling non-idle delta above five was accepted")
    regressed = parse_cpu_tick_receipt(
        b"cpu56 15 2 3 99 4 5 6 7 5000 4000\n", 56, 40
    )
    try:
        sibling_tick_delta(tick_start, regressed, 20, 30)
    except ValidationError:
        pass
    else:
        raise AssertionError("regressed sibling tick component was accepted")

    normal_index_flags = b"".join(
        ("H " + relative + "\n").encode("ascii")
        for relative in SOURCE_PATHS
    )
    validate_git_index_flag_roster(normal_index_flags)
    for marker in (b"h", b"S"):
        mutated_index_flags = marker + normal_index_flags[1:]
        try:
            validate_git_index_flag_roster(mutated_index_flags)
        except ValidationError:
            pass
        else:
            raise AssertionError(
                "hidden Git index flag was accepted: "
                + marker.decode("ascii")
            )

    if not sys.platform.startswith("linux"):
        print("WH2 broad systematic emission analyzer selftest passed")
        return 0
    sleeper = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(10)"],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
    )
    capture = bounded_capture(sleeper, time.monotonic() + 0.05)
    if not capture.timed_out or capture.output_overflow or sleeper.returncode is None:
        raise AssertionError("bounded timeout did not kill and reap")
    flooder = subprocess.Popen(
        [sys.executable, "-c", "import os; os.write(2, b'x' * 1114112)"],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
    )
    capture = bounded_capture(flooder, time.monotonic() + 5.0)
    if (
        capture.timed_out or not capture.output_overflow
        or flooder.returncode is None
        or len(capture.stderr) != MAX_STDERR_BYTES + 1
    ):
        raise AssertionError("bounded output overflow did not kill and reap")
    with tempfile.TemporaryDirectory(prefix="wh2-broad-systematic-selftest-") as tmp:
        root = Path(tmp)
        output = root / "one-shot"
        selftest_source_manifest = "2" * 64
        bundle = OutputBundle.reserve(
            output, cpu, commit, selftest_source_manifest, project_root
        )
        try:
            try:
                OutputBundle.reserve(
                    output, cpu, commit, selftest_source_manifest,
                    project_root,
                )
            except ValidationError:
                pass
            else:
                raise AssertionError("output directory reuse was allowed")
        finally:
            failures = bundle.stage_evidence(b"sealed\n", b"", False)
            candidate_summary = bundle.emergency_summary(
                "selftest-invalid-evidence", failures
            )
            parsed_candidate = parse_summary_bytes(candidate_summary)
            for label, path, replacement in (
                (
                    "true output nlink",
                    ("output_bundle", "files", "raw.jsonl", "nlink"),
                    True,
                ),
                (
                    "string output device",
                    ("output_bundle", "directory", "device"),
                    "0",
                ),
            ):
                mutated_summary = json.loads(
                    canonical_bytes(parsed_candidate).decode("ascii")
                )
                target = mutated_summary
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = replacement
                mutated_summary["summary_preimage_sha256"] = None
                mutated_bytes = canonical_bytes(
                    make_summary_preimage(mutated_summary)
                ) + b"\n"
                try:
                    parse_summary_bytes(mutated_bytes)
                except ValidationError:
                    pass
                else:
                    raise AssertionError(
                        f"nested summary mutation was accepted: {label}"
                    )
            rewrite_diagnostics = [
                f"rewrite-diagnostic-{index:02d}:" + "x" * 1000
                for index in range(MAX_PUBLICATION_FAILURES + 6)
            ]
            rewrite_base = parsed_candidate
            stable_failure: Optional[str] = None
            rewrite_attempts = len(SignalGuard.SIGNALS) + 8
            if rewrite_attempts <= len(SignalGuard.SIGNALS) + 6:
                raise AssertionError("rewrite stability selftest is too short")
            for _ in range(rewrite_attempts):
                rewritten = bundle._invalid_summary(
                    rewrite_base, rewrite_diagnostics, None
                )
                rewrite_base = parse_summary_bytes(rewritten)
                current_failure = rewrite_base["failure"]
                if (
                    len(current_failure) > MAX_FAILURE_TEXT_CHARS
                    or len(rewrite_base["publication_failures"])
                    != MAX_PUBLICATION_FAILURES
                    or (
                        stable_failure is not None
                        and current_failure != stable_failure
                    )
                ):
                    raise AssertionError(
                        "repeated invalid-summary failure text is not stable"
                    )
                stable_failure = current_failure
            finish_result = bundle.finish(candidate_summary)
            failures.extend(finish_result.all_failures())
            if not finish_result.logical_commit_succeeded:
                failures.append("selftest output lacked a logical commit")
            if failures:
                raise AssertionError(f"output bundle seal failed: {failures}")
        if stat.S_IMODE(output.stat().st_mode) != 0o500:
            raise AssertionError("sealed output directory mode differs")
        if any(
            stat.S_IMODE((output / name).stat().st_mode) != 0o400
            for name in OUTPUT_NAMES
        ):
            raise AssertionError("sealed output mode differs")
        preserved_output = root / "preserved-invalid-one-shot"
        preserved_bundle = OutputBundle.reserve(
            preserved_output, cpu, commit, selftest_source_manifest,
            project_root,
        )
        failures = preserved_bundle.stage_evidence(raw, b"", True)
        preserved_binary = empty_binary_receipt(
            root / "captured-screen", "3" * 64
        )
        preserved_binary["capture_error"] = "preserved health/postcheck failure"
        preserved_bundle.bind_expected_summary_components({
            "binary": preserved_binary,
            "config": records[0],
            "elapsed_seconds": 1.0,
            "failure": "health finish: injected failure",
            "git_after": {}, "git_before": {}, "health": {},
            "health_module_loader": {}, "outcome": "invalid",
            "process_exit_code": 0, "source_manifest_after": {},
            "source_manifest_before": {}, "statistics": statistics,
            "target_cpu": cpu,
        })
        preserved_summary = preserved_bundle.emergency_summary(
            "post-capture adjudication failure", failures
        )
        finish_result = preserved_bundle.finish(preserved_summary)
        parsed_preserved = parse_summary_bytes(
            (preserved_output / "summary.json").read_bytes()
        )
        if (
            failures
            or not finish_result.logical_commit_succeeded
            or finish_result.all_failures()
            or not parsed_preserved["raw_complete"]
            or not canonical_equal(parsed_preserved["binary"], preserved_binary)
            or not canonical_equal(parsed_preserved["config"], records[0])
            or not canonical_equal(parsed_preserved["statistics"], statistics)
        ):
            raise AssertionError(
                "post-capture invalid summary discarded sealed evidence"
            )
        signal_output = root / "signal-one-shot"
        with SignalGuard() as signal_guard:
            signal_bundle = OutputBundle.reserve(
                signal_output, cpu, commit, selftest_source_manifest,
                project_root, signal_guard,
            )
            failures = signal_bundle.stage_evidence(b"prefix\n", b"", False)
            candidate_summary = signal_bundle.emergency_summary(
                "pre-signal-candidate", failures
            )
            signal_guard._handler(signal.SIGTERM, None)
            signal_guard._handler(signal.SIGINT, None)
            finish_result = signal_bundle.finish(candidate_summary, signal_guard)
            failures.extend(finish_result.all_failures())
            if not finish_result.logical_commit_succeeded:
                failures.append("signal output lacked a logical commit")
            if failures:
                raise AssertionError(
                    f"signal-invalid output seal failed: {failures}"
                )
        signal_summary = parse_summary_bytes(
            (signal_output / "summary.json").read_bytes()
        )
        if (
            signal_summary["outcome"] != "invalid"
            or signal_summary["signal"] != "SIGTERM"
            or signal_summary["signal_names"] != ["SIGTERM", "SIGINT"]
        ):
            raise AssertionError("guarded signal commit receipt differs")
        postcommit_output = root / "postcommit-one-shot"
        postcommit_bundle = OutputBundle.reserve(
            postcommit_output, cpu, commit, selftest_source_manifest,
            project_root,
        )
        failures = postcommit_bundle.stage_evidence(b"postcommit\n", b"", False)
        postcommit_summary = postcommit_bundle.emergency_summary(
            "postcommit-selftest-candidate", failures
        )
        def injected_raising_hardener() -> List[str]:
            raise RuntimeError("injected postcommit durability failure")

        postcommit_bundle._harden_after_commit = (  # type: ignore[method-assign]
            injected_raising_hardener
        )
        finish_result = postcommit_bundle.finish(postcommit_summary)
        if (
            failures
            or not finish_result.logical_commit_succeeded
            or finish_result.precommit_failures
            or finish_result.postcommit_failures
            != (
                "postcommit hardening call: RuntimeError: "
                "injected postcommit durability failure",
            )
            or (postcommit_output / "summary.json").read_bytes()
            != postcommit_summary
            or stat.S_IMODE(postcommit_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((postcommit_output / name).stat().st_mode) != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "postcommit failure changed the immutable summary decision"
            )
        positive_raw = synthetic_transcript(EXPECTED_TARGET_CPU, commit)
        positive_records = parse_transcript(positive_raw)
        positive_config, positive_statistics, positive_gates = (
            validate_transcript(
                positive_records, EXPECTED_TARGET_CPU, commit
            )
        )
        if positive_gates:
            raise AssertionError("positive publication fixture failed its gates")

        def make_positive_summary(
            positive_bundle: OutputBundle,
        ) -> Tuple[bytes, List[str]]:
            stage_errors = positive_bundle.stage_evidence(
                positive_raw, b"", True
            )
            fake_stat = {
                "device": 1, "inode": 1, "mode": 0o555,
                "mtime_ns": 1, "nlink": 1, "size": 1,
            }
            binary_digest = "4" * 64
            binary_receipt = {
                "capture_error": "none",
                "child_reap_monotonic_ns": 200,
                "child_start_monotonic_ns": 100,
                "expected_sha256": binary_digest,
                "output_overflow": False,
                "path": str(root / "positive-binary"),
                "path_stat_after": dict(fake_stat),
                "process_started": True,
                "sha256_after": binary_digest,
                "sha256_before": binary_digest,
                "stat_after": dict(fake_stat),
                "stat_before": dict(fake_stat),
                "timed_out": False,
            }
            sampler = {
                "cpu": 1, "csv_device": 1, "csv_inode": 2,
                "csv_path": str(root / "sampler.csv"), "pid": 1,
                "process_start_ticks": 1, "schema": SAMPLER_SCHEMA,
                "script_path": str(root / "sampler.py"),
                "script_sha256": "5" * 64, "terminal_status": "ok",
                "window_end_monotonic_ns": 210,
                "window_start_monotonic_ns": 90,
            }
            thermal = {
                "cpu": 1, "cpu_tctl_max_millic": 50000,
                "csv_device": 1, "csv_inode": 2,
                "csv_path": str(root / "sampler.csv"),
                "dimm_max_millic": 50000, "dimm_read_errors": 0,
                "edac_ce_max": 0, "edac_ue_max": 0, "pid": 1,
                "process_start_ticks": 1, "sample_count": 2,
                "script_path": str(root / "sampler.py"),
                "script_sha256": "5" * 64, "terminal_status": "ok",
                "window_csv_sha256": "6" * 64,
                "window_end_monotonic_ns": 210,
                "window_start_monotonic_ns": 90,
            }
            tick_start_receipt = {
                "cpu": 56, "non_idle_ticks": 10,
                "read_monotonic_ns": 90,
                "tick_fields": {
                    "idle": 100, "iowait": 0, "irq": 0, "nice": 0,
                    "softirq": 0, "steal": 0, "system": 0, "user": 10,
                },
            }
            tick_end_receipt = {
                "cpu": 56, "non_idle_ticks": 11,
                "read_monotonic_ns": 210,
                "tick_fields": {
                    "idle": 101, "iowait": 0, "irq": 0, "nice": 0,
                    "softirq": 0, "steal": 0, "system": 0, "user": 11,
                },
            }
            health_receipt: Dict[str, Any] = {
                "child_reap_monotonic_ns": 200,
                "child_start_monotonic_ns": 100,
                "controller_core": [0, 0], "controller_cpu": 0,
                "controller_initial_affinity": [0, 120],
                "controller_singleton_affinity_end": [0],
                "edac_policy": "ce-and-ue-every-sample-zero-v1",
                "receipt_sha256": None, "sampler": sampler,
                "sampler_core": [0, 1], "sampler_cpu": 1,
                "sampler_receipt_sha256": sha256_bytes(
                    canonical_bytes(sampler)
                ),
                "schema": HEALTH_SCHEMA,
                "sibling_non_idle_tick_cap": SIBLING_NON_IDLE_TICK_CAP,
                "sibling_tick_policy": (
                    "linux-proc-stat-user-nice-system-irq-softirq-steal-v1"
                ),
                "sibling_ticks": [{
                    "cpu": 56, "delta_non_idle_ticks": 1,
                    "end": tick_end_receipt, "start": tick_start_receipt,
                }],
                "target_core": list(EXPECTED_TARGET_CORE),
                "target_cpu": EXPECTED_TARGET_CPU,
                "target_threads": list(EXPECTED_TARGET_THREADS),
                "thermal": thermal,
                "thermal_max_millic": THERMAL_MAX_MILLIC,
                "terminal_status": "ok",
            }
            health_receipt["receipt_sha256"] = sha256_bytes(
                canonical_bytes(health_receipt)
            )
            git_value = {
                "head": commit,
                "source_roster_sha256": sha256_bytes(b"".join(
                    (relative + "\n").encode("ascii")
                    for relative in SOURCE_PATHS
                )),
                "tracked_status_sha256": sha256_bytes(b""),
            }
            components = {
                "binary": binary_receipt,
                "config": positive_config,
                "elapsed_seconds": 1.0,
                "failure": "none", "git_after": git_value,
                "git_before": git_value, "health": health_receipt,
                "health_module_loader": loaded_health.receipt,
                "outcome": "pass", "process_exit_code": 0,
                "source_manifest_after": loader_manifest,
                "source_manifest_before": loader_manifest,
                "statistics": positive_statistics,
                "target_cpu": EXPECTED_TARGET_CPU,
            }
            positive_bundle.bind_expected_summary_components(components)
            summary = make_summary_preimage({
                **components,
                "expected_source_git_commit": commit,
                "expected_source_manifest_sha256": loader_manifest["sha256"],
                "output_bundle": positive_bundle.receipt(),
                "publication_failures": [],
                "raw_bytes": len(positive_raw), "raw_complete": True,
                "raw_record_count": positive_raw.count(b"\n"),
                "raw_sha256": sha256_bytes(positive_raw),
                "schema": SUMMARY_SCHEMA, "signal": None,
                "signal_names": [], "stderr_bytes": 0,
                "stderr_sha256": sha256_bytes(b""),
                "summary_preimage_sha256": None,
            })
            summary_bytes = canonical_bytes(summary) + b"\n"
            if parse_summary_bytes(summary_bytes)["outcome"] != "pass":
                raise AssertionError("positive publication fixture is not valid")
            return summary_bytes, stage_errors

        source_recheck_output = root / "source-recheck-one-shot"
        source_recheck_bundle = OutputBundle.reserve(
            source_recheck_output, EXPECTED_TARGET_CPU, commit,
            loader_manifest["sha256"], project_root,
        )
        positive_summary, failures = make_positive_summary(
            source_recheck_bundle
        )
        positive_object = parse_summary_bytes(positive_summary)
        mutated_manifest = json.loads(
            canonical_bytes(loader_manifest).decode("ascii")
        )
        mutated_manifest["entries"][0]["sha256"] = "0" * 64
        mutated_preimage = b"".join(
            (
                entry["sha256"] + "  " + entry["path"] + "\n"
            ).encode("ascii")
            for entry in mutated_manifest["entries"]
        )
        mutated_manifest["sha256"] = sha256_bytes(mutated_preimage)
        validate_source_manifest_receipt(
            mutated_manifest, "injected final source mutation"
        )
        source_reads = [loader_manifest, mutated_manifest]

        def injected_live_source_manifest() -> Dict[str, Any]:
            return source_reads.pop(0) if source_reads else mutated_manifest

        def injected_live_git_receipt() -> Dict[str, Any]:
            return positive_object["git_before"]

        source_recheck_bundle._live_source_manifest = (  # type: ignore[method-assign]
            injected_live_source_manifest
        )
        source_recheck_bundle._live_git_receipt = (  # type: ignore[method-assign]
            injected_live_git_receipt
        )
        finish_result = source_recheck_bundle.finish(positive_summary)
        parsed_source_recheck = parse_summary_bytes(
            (source_recheck_output / "summary.json").read_bytes()
        )
        source_recheck_failure = (
            "mutable precommit check: ValidationError: source manifest "
            "changed during final source/Git sandwich"
        )
        if (
            failures
            or not finish_result.logical_commit_succeeded
            or source_recheck_failure not in finish_result.precommit_failures
            or parsed_source_recheck["outcome"] != "invalid"
            or source_recheck_failure
            not in parsed_source_recheck["publication_failures"]
            or source_reads
            or stat.S_IMODE(source_recheck_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((source_recheck_output / name).stat().st_mode)
                != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "final live source mutation retained a positive summary"
            )

        close_output = root / "close-only-one-shot"
        close_bundle = OutputBundle.reserve(
            close_output, EXPECTED_TARGET_CPU, commit,
            loader_manifest["sha256"], project_root,
        )
        positive_summary, failures = make_positive_summary(close_bundle)
        positive_object = parse_summary_bytes(positive_summary)
        close_bundle._live_source_manifest = (  # type: ignore[method-assign]
            lambda: loader_manifest
        )
        close_bundle._live_git_receipt = (  # type: ignore[method-assign]
            lambda: positive_object["git_before"]
        )
        summary_fd = close_bundle.file_fds["summary.json"]
        original_close_owned_fd = close_bundle._close_owned_fd
        close_diagnostics_injected: List[int] = []

        def injected_close_only_diagnostic(fd: int) -> None:
            original_close_owned_fd(fd)
            if fd == summary_fd:
                close_diagnostics_injected.append(fd)
                raise RuntimeError("injected close-only diagnostic")

        close_bundle._close_owned_fd = (  # type: ignore[method-assign]
            injected_close_only_diagnostic
        )
        finish_result = close_bundle.finish(positive_summary)
        if (
            failures
            or not finish_result.logical_commit_succeeded
            or finish_result.all_failures()
            or close_diagnostics_injected != [summary_fd]
            or (close_output / "summary.json").read_bytes()
            != positive_summary
            or parse_summary_bytes(positive_summary)["outcome"] != "pass"
            or stat.S_IMODE(close_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((close_output / name).stat().st_mode) != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "close-only diagnostic changed a committed positive decision"
            )

        poison_output = root / "poison-one-shot"
        poison_bundle = OutputBundle.reserve(
            poison_output, EXPECTED_TARGET_CPU, commit,
            loader_manifest["sha256"], project_root,
        )
        positive_summary, failures = make_positive_summary(poison_bundle)
        original_write_exact = poison_bundle._write_exact

        def reject_invalid_replacement(fd: int, data: bytes) -> bytes:
            if (
                fd == poison_bundle.file_fds["summary.json"]
                and b'"outcome":"invalid"' in data
            ):
                raise RuntimeError("injected persistent replacement failure")
            return original_write_exact(fd, data)

        poison_bundle._write_exact = (  # type: ignore[method-assign]
            reject_invalid_replacement
        )
        poison_bundle._mutable_precommit_pass = (  # type: ignore[method-assign]
            lambda: ["injected positive-candidate precommit failure"]
        )
        finish_result = poison_bundle.finish(positive_summary)
        if (
            failures
            or finish_result.logical_commit_succeeded
            or not finish_result.precommit_failures
            or (poison_output / "summary.json").read_bytes() != b""
            or stat.S_IMODE(poison_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((poison_output / name).stat().st_mode) != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "persistent invalid rewrite retained a positive summary"
            )
        unlink_output = root / "poison-unlink-one-shot"
        unlink_bundle = OutputBundle.reserve(
            unlink_output, EXPECTED_TARGET_CPU, commit,
            loader_manifest["sha256"], project_root,
        )
        positive_summary, failures = make_positive_summary(unlink_bundle)
        original_write_exact = unlink_bundle._write_exact

        def reject_unlink_invalid_replacement(fd: int, data: bytes) -> bytes:
            if (
                fd == unlink_bundle.file_fds["summary.json"]
                and b'"outcome":"invalid"' in data
            ):
                raise RuntimeError("injected persistent replacement failure")
            return original_write_exact(fd, data)

        unlink_bundle._write_exact = (  # type: ignore[method-assign]
            reject_unlink_invalid_replacement
        )
        unlink_bundle._mutable_precommit_pass = (  # type: ignore[method-assign]
            lambda: ["injected positive-candidate precommit failure"]
        )

        def injected_poison_truncate_failure() -> None:
            raise RuntimeError("injected poison truncation failure")

        unlink_bundle._truncate_summary_for_poison = (  # type: ignore[method-assign]
            injected_poison_truncate_failure
        )
        finish_result = unlink_bundle.finish(positive_summary)
        if (
            failures
            or finish_result.logical_commit_succeeded
            or not finish_result.precommit_failures
            or (unlink_output / "summary.json").exists()
            or stat.S_IMODE(unlink_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((unlink_output / name).stat().st_mode) != 0o400
                for name in ("raw.jsonl", "stderr.txt")
            )
        ):
            raise AssertionError(
                "failed summary truncation did not use identity-checked unlink"
            )
        bookkeeping_output = root / "bookkeeping-exception-one-shot"
        bookkeeping_bundle = OutputBundle.reserve(
            bookkeeping_output, EXPECTED_TARGET_CPU, commit,
            loader_manifest["sha256"], project_root,
        )
        positive_summary, failures = make_positive_summary(
            bookkeeping_bundle
        )

        def injected_invalid_precommit_result() -> Any:
            return None

        bookkeeping_bundle._mutable_precommit_pass = (  # type: ignore[method-assign]
            injected_invalid_precommit_result
        )
        finish_result = bookkeeping_bundle.finish(positive_summary)
        parsed_bookkeeping = parse_summary_bytes(
            (bookkeeping_output / "summary.json").read_bytes()
        )
        if (
            failures
            or finish_result.logical_commit_succeeded
            or not finish_result.precommit_failures
            or parsed_bookkeeping["outcome"] != "invalid"
            or stat.S_IMODE(bookkeeping_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((bookkeeping_output / name).stat().st_mode)
                != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "commit-loop bookkeeping exception retained a positive summary"
            )
        precommit_output = root / "precommit-one-shot"
        precommit_bundle = OutputBundle.reserve(
            precommit_output, cpu, commit, selftest_source_manifest,
            project_root,
        )
        failures = precommit_bundle.stage_evidence(
            b"precommit\n", b"", False
        )
        precommit_summary = precommit_bundle.emergency_summary(
            "precommit-selftest-candidate", failures
        )

        def injected_raising_precommit() -> List[str]:
            raise RuntimeError("injected mutable precommit failure")

        precommit_bundle._mutable_precommit_pass = (  # type: ignore[method-assign]
            injected_raising_precommit
        )
        finish_result = precommit_bundle.finish(precommit_summary)
        parsed_precommit = parse_summary_bytes(
            (precommit_output / "summary.json").read_bytes()
        )
        precommit_failure = (
            "mutable precommit call: RuntimeError: "
            "injected mutable precommit failure"
        )
        if (
            failures
            or finish_result.logical_commit_succeeded
            or precommit_failure not in finish_result.precommit_failures
            or parsed_precommit["outcome"] != "invalid"
            or precommit_failure not in parsed_precommit["publication_failures"]
            or stat.S_IMODE(precommit_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((precommit_output / name).stat().st_mode) != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "raising mutable-precommit pass escaped invalid publication"
            )
        commit_output = root / "commit-transition-one-shot"
        with SignalGuard() as commit_guard:
            commit_bundle = OutputBundle.reserve(
                commit_output, cpu, commit, selftest_source_manifest,
                project_root, commit_guard,
            )
            failures = commit_bundle.stage_evidence(
                b"commit-transition\n", b"", False
            )
            commit_summary = commit_bundle.emergency_summary(
                "commit-transition-selftest-candidate", failures
            )

            def injected_raising_commit() -> None:
                raise RuntimeError("injected logical commit failure")

            commit_guard.commit_logical_seal = (  # type: ignore[method-assign]
                injected_raising_commit
            )
            finish_result = commit_bundle.finish(
                commit_summary, commit_guard
            )
        parsed_commit = parse_summary_bytes(
            (commit_output / "summary.json").read_bytes()
        )
        commit_failure = (
            "logical signal commit: RuntimeError: "
            "injected logical commit failure"
        )
        if (
            failures
            or finish_result.logical_commit_succeeded
            or commit_failure not in finish_result.precommit_failures
            or parsed_commit["outcome"] != "invalid"
            or commit_failure not in parsed_commit["publication_failures"]
            or stat.S_IMODE(commit_output.stat().st_mode) != 0o500
            or any(
                stat.S_IMODE((commit_output / name).stat().st_mode) != 0o400
                for name in OUTPUT_NAMES
            )
        ):
            raise AssertionError(
                "raising commit transition escaped invalid publication"
            )
        helper = root / "helper.py"
        helper.write_bytes(
            b"#!" + os.fsencode(sys.executable) + b"\n"
            b"import sys\nsys.stdout.write('helper-ok\\n')\n"
        )
        helper.chmod(0o555)
        helper_digest = sha256_bytes(helper.read_bytes())
        helper.chmod(0o755)
        try:
            writable_result = run_binary(helper, 0, helper_digest)
        finally:
            helper.chmod(0o555)
        if (
            not writable_result.errors
            or writable_result.binary["process_started"]
        ):
            raise AssertionError(
                "owner-writable same-hash helper was accepted"
            )
        result = run_binary(helper, 0, helper_digest)
        if (
            result.raw != b"helper-ok\n" or result.stderr
            or result.returncode != 0 or result.errors
            or result.binary["sha256_before"] != helper_digest or result.health
        ):
            raise AssertionError("sealed helper execution receipt differs")
        result = run_binary(helper, 0, "0" * 64)
        if not result.errors or result.binary["process_started"]:
            raise AssertionError("wrong expected binary hash was accepted")
    print("WH2 broad systematic emission analyzer selftest passed")
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--binary")
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--controller-cpu", type=int)
    parser.add_argument("--output-dir")
    parser.add_argument("--sampler-pid", type=int)
    parser.add_argument("--sampler-cpu", type=int)
    parser.add_argument("--sampler-script", type=Path)
    parser.add_argument("--sampler-csv", type=Path)
    parser.add_argument("--expected-source-commit")
    parser.add_argument("--expected-binary-sha256")
    parser.add_argument("--expected-source-manifest-sha256")
    args = parser.parse_args(argv)
    if (
        not sys.flags.isolated
        or not sys.dont_write_bytecode
        or sys.flags.optimize != 0
    ):
        parser.error(
            "this sealed runner requires unoptimized Python -I -B"
        )
    if args.selftest:
        return args
    if (
        not args.binary
        or not Path(args.binary).is_absolute()
        or os.path.normpath(args.binary) != args.binary
        or args.cpu != EXPECTED_TARGET_CPU
        or args.controller_cpu is None or args.controller_cpu < 0
        or not args.output_dir
        or args.sampler_pid is None or args.sampler_pid <= 0
        or args.sampler_cpu is None or args.sampler_cpu < 0
        or args.sampler_script is None or not args.sampler_script.is_absolute()
        or args.sampler_csv is None or not args.sampler_csv.is_absolute()
        or type(args.expected_source_commit) is not str
        or LOWER40.fullmatch(args.expected_source_commit) is None
        or type(args.expected_binary_sha256) is not str
        or LOWER64.fullmatch(args.expected_binary_sha256) is None
        or type(args.expected_source_manifest_sha256) is not str
        or LOWER64.fullmatch(args.expected_source_manifest_sha256) is None
    ):
        parser.error(
            "--run requires the frozen target/controller, one live absolute "
            "sampler, one new output directory, and commit/binary/source seals"
        )
    return args


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    return selftest() if args.selftest else run_once(args)


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except (ValidationError, AssertionError, OSError) as exc:
        print(f"broad systematic emission analyzer failed: {exc}", file=sys.stderr)
        sys.exit(1)
