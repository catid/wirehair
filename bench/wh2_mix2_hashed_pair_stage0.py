#!/usr/bin/env python3
"""Sealed consumed-coordinate Stage 0 for row-hashed MIX2 pairs.

This controller authenticates the complete pair01/V9/fixed-pair history before
running the row-hashed candidate on the three immutable pair01 weak cells.  It
uses no fresh root, records no timing evidence, and cannot promote a profile.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from bench import wh2_mix2_pair_stage0 as pair_stage0
    from bench import wh2_mix2_pair_stage1 as pair_stage1
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_pair_stage0 as pair_stage0
    import wh2_mix2_pair_stage1 as pair_stage1
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-hashed-pair-stage0-contract.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-hashed-pair-stage0-record.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-hashed-pair-stage0-summary.v1"
RECORD_NAME = "mix2-hashed-pair-stage0-records.jsonl"
SUMMARY_NAME = "mix2-hashed-pair-stage0-summary.json"

CANDIDATE_PAIR = "hashed"
PAIR_SELECTOR_SALT = "0x4d495832"
REPETITIONS = (0, 1)
BASELINE_COORDINATE_ORDINALS = (30, 573, 820)
BASELINE_MATRIX_ORDINALS = (354, 138273, 202780)
EXPECTED_INVOCATION_COUNT = 6
INACTIVATION_CAP = 1024

BLOCK_BYTES = pair_stage1.BLOCK_BYTES
LOSS_PPM = pair_stage1.LOSS_PPM
OVERHEAD = pair_stage1.OVERHEAD
BINARY_DENSE_ROWS = pair_stage1.BINARY_DENSE_ROWS
GF256_HEAVY_ROWS = pair_stage1.GF256_HEAVY_ROWS
ZERO_SHA256 = pair_stage1.ZERO_SHA256
CSV_HEADER = tuple(pair_stage0.CSV_HEADER)
METADATA_ORDER = tuple(pair_stage0.METADATA_ORDER)
STAIRCASE_BY_K = {
    K: pair_stage1.STAIRCASE_BY_K[K] for K in (2, 1000, 10000)
}
BENCH_DESCRIPTION_SCHEMA = pair_stage0.BENCH_DESCRIPTION_SCHEMA
CHILD_ENVIRONMENT = dict(pair_stage0.CHILD_ENVIRONMENT)

EXPECTED_BASELINE_DIR = Path(
    "/dev/shm/wh2-mix2-attempt-crossfit-c9cfbc0-run")
EXPECTED_V9_DIR = Path("/dev/shm/wh2-mix2-v9-short-52177d5d7989-run")
EXPECTED_STAGE0_DIR = Path(
    "/dev/shm/wh2-mix2-pair-stage0-99acef9-run")

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent

# pair_stage1's source roster already contains both imported pair controllers,
# the cross-fit/V9 authenticators, the benchmark, and every linked codec file.
SOURCE_PATHS = (
    "bench/wh2_mix2_hashed_pair_stage0.py",
    "bench/test_wh2_mix2_hashed_pair_stage0.py",
) + tuple(pair_stage1.SOURCE_PATHS)
SOURCE_STATUS_PATHS = (
    "bench/wh2_mix2_hashed_pair_stage0.py",
    "bench/test_wh2_mix2_hashed_pair_stage0.py",
) + tuple(pair_stage1.SOURCE_STATUS_PATHS)

SHA256 = pair_stage0.SHA256
COMMIT = pair_stage0.COMMIT
canonical_json = pair_stage0.canonical_json
sha256_bytes = pair_stage0.sha256_bytes
sha256_json = pair_stage0.sha256_json
self_hash = pair_stage0.self_hash
ProcessResult = pair_stage0.ProcessResult


class HashedPairStage0Error(RuntimeError):
    """The sealed hashed-pair screen cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise HashedPairStage0Error(message)


@dataclass(frozen=True)
class Coordinate:
    ordinal: int

    def __post_init__(self) -> None:
        if type(self.ordinal) is not int or \
                self.ordinal not in BASELINE_COORDINATE_ORDINALS:
            fail("coordinate ordinal is outside hashed-pair Stage 0")

    @property
    def stage1(self) -> pair_stage1.Coordinate:
        return pair_stage1.Coordinate(self.ordinal)

    @property
    def selection_index(self) -> int:
        return BASELINE_COORDINATE_ORDINALS.index(self.ordinal)

    @property
    def K(self) -> int:
        return self.stage1.K

    @property
    def attempt(self) -> int:
        return self.stage1.attempt

    @property
    def root_index(self) -> int:
        return self.stage1.root_index

    @property
    def root(self) -> str:
        return self.stage1.root

    @property
    def schedule_index(self) -> int:
        return self.stage1.schedule_index

    @property
    def schedule(self) -> str:
        return self.stage1.schedule

    @property
    def cell_ordinal(self) -> int:
        return self.stage1.cell_ordinal

    @property
    def baseline_record_ordinal(self) -> int:
        return self.stage1.baseline_record_ordinal

    def identity(self) -> Mapping[str, Any]:
        return {
            "selection_index": self.selection_index,
            **self.stage1.identity(),
        }

    def description(self) -> Mapping[str, Any]:
        return dict(self.identity())


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    coordinate_ordinal: int
    repetition: int

    def __post_init__(self) -> None:
        if (type(self.ordinal) is not int or
                type(self.coordinate_ordinal) is not int or
                type(self.repetition) is not int or
                self.coordinate_ordinal not in
                    BASELINE_COORDINATE_ORDINALS or
                self.repetition not in REPETITIONS):
            fail("invocation identity is outside hashed-pair Stage 0")
        selection = BASELINE_COORDINATE_ORDINALS.index(
            self.coordinate_ordinal)
        if self.ordinal != selection * len(REPETITIONS) + self.repetition:
            fail("invocation differs from coordinate-major repetition order")

    @property
    def pair(self) -> str:
        return CANDIDATE_PAIR

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
            str(coordinate.attempt), "--construction-seed-basis",
            "production-profile",
        ]

    def identity(self) -> Mapping[str, Any]:
        return {
            "ordinal": self.ordinal,
            "pair": CANDIDATE_PAIR,
            "repetition": self.repetition,
            **self.coordinate.identity(),
        }


def make_invocations() -> Tuple[Invocation, ...]:
    invocations = tuple(
        Invocation(len(REPETITIONS) * selection + repetition,
                   coordinate_ordinal, repetition)
        for selection, coordinate_ordinal in enumerate(
            BASELINE_COORDINATE_ORDINALS)
        for repetition in REPETITIONS
    )
    if len(invocations) != EXPECTED_INVOCATION_COUNT or len({
            (item.coordinate_ordinal, item.repetition)
            for item in invocations}) != EXPECTED_INVOCATION_COUNT:
        fail("hashed-pair invocation roster is incomplete or duplicated")
    return invocations


def _contract_body() -> Mapping[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "bead": "wirehair-sxvz.16.1.20.39",
        "scope": "row-hashed MIX2 consumed-coordinate Stage 0 only",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "timing_evidence": False,
        "candidate": {
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "mix_pair_token": CANDIDATE_PAIR,
            "test_mode": 3,
            "selector": (
                "x=Avalanche32(block_id^PacketRowConfig.PeelSeed^0x4d495832);"
                "pair=high32(uint64(x)*3)"),
            "pair_map": {"0": [0, 1], "1": [0, 2], "2": [1, 2]},
            "selector_full_uint32_counts": [
                1431655766, 1431655765, 1431655765],
            "binary_dense_rows": BINARY_DENSE_ROWS,
            "gf256_heavy_rows": GF256_HEAVY_ROWS,
            "heavy_family": "periodic",
            "construction_seed_basis": "production-profile",
            "attempt_rule": "exact V9 precode attempt equals packet attempt",
        },
        "domain": {
            "coordinates": [
                Coordinate(value).description()
                for value in BASELINE_COORDINATE_ORDINALS
            ],
            "baseline_coordinate_ordinals":
                list(BASELINE_COORDINATE_ORDINALS),
            "baseline_matrix_ordinals": list(BASELINE_MATRIX_ORDINALS),
            "coordinate_order": "coordinate-major, repetition 0 then 1",
            "repetitions": list(REPETITIONS),
            "invocation_count": EXPECTED_INVOCATION_COUNT,
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "trials": 1,
            "threads": 1,
            "overhead_stream": "paired",
            "full_payload_solve": True,
            "cold_solve_wide_xor": "policy",
        },
        "historical_prerequisites": {
            "directories": {
                "baseline": str(EXPECTED_BASELINE_DIR),
                "v9": str(EXPECTED_V9_DIR),
                "fixed_pair_stage0": str(EXPECTED_STAGE0_DIR),
            },
            "authentication": (
                "pair_stage1.load_historical_evidence completes before any "
                "candidate process"),
            "pair01_baseline": {
                "contract_sha256":
                    pair_stage1.EXPECTED_CROSSFIT_CONTRACT_SHA256,
                "source_commit": pair_stage1.EXPECTED_CROSSFIT_SOURCE_COMMIT,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in pair_stage1.BASELINE_FILES.items()
                },
                "selected_coordinate_ordinals":
                    list(BASELINE_COORDINATE_ORDINALS),
                "selected_matrix_ordinals": list(BASELINE_MATRIX_ORDINALS),
            },
            "v9": {
                "contract_sha256": pair_stage1.EXPECTED_V9_CONTRACT_SHA256,
                "source_commit": pair_stage1.EXPECTED_V9_SOURCE_COMMIT,
                "attempt_map_sha256": pair_stage1.EXPECTED_ATTEMPT_MAP_SHA256,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in pair_stage1.V9_FILES.items()
                },
            },
            "fixed_pair_stage0": {
                "contract_sha256":
                    pair_stage1.EXPECTED_STAGE0_CONTRACT_SHA256,
                "source_commit": pair_stage1.EXPECTED_STAGE0_SOURCE_COMMIT,
                "files": {
                    name: {"size": spec[0], "sha256": spec[1]}
                    for name, spec in pair_stage1.STAGE0_FILES.items()
                },
            },
        },
        "transport": {
            "benchmark_pin": "open regular executable descriptor and rehash",
            "benchmark_preflight": {
                "command": "--describe",
                "schema": BENCH_DESCRIPTION_SCHEMA,
                "source_git_commit": "must equal clean source HEAD",
                "ordering": "last action before candidate precodefail work",
            },
            "source_receipt": {
                "clean_at_HEAD": True,
                "paths": list(SOURCE_PATHS),
                "includes_all_imported_linked_source": True,
            },
            "process_receipts": (
                "SHA256 of exact argv excluding executable plus raw stdout "
                "and stderr bytes for every invocation"),
            "child_environment": dict(CHILD_ENVIRONMENT),
            "publication": "canonical JSONL and summary, atomic no-replace dir",
        },
        "gate": {
            "all_six_observations_succeed": True,
            "all_three_repeated_deterministic_projections_exact": True,
            "every_observation_inactivated_columns_at_most":
                INACTIVATION_CAP,
            "per_coordinate_vs_immutable_pair01": {
                "block_xors": "candidate <= baseline",
                "gf256_muladds": "candidate <= baseline",
                "inactivated_columns": "candidate <= baseline",
                "effective_precode_seed": "candidate == baseline",
                "effective_packet_seed": "candidate == baseline",
            },
            "overall": "PASS iff every listed recovery, repeat, cap, and work gate passes",
        },
        "stop_rule": (
            "no fresh roots, timing claim, production change, all-K work, "
            "or promotion is licensed by this controller"),
    }


# Replaced once after the complete literal contract body is finalized.
EXPECTED_CONTRACT_SHA256 = (
    "599475efd2bce98f41e3ad13fddcc7e9edf6d8b47b729bb2d56b39eff95feede")


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("hashed-pair Stage 0 contract differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _validate_constants() -> None:
    try:
        pair_stage0._validate_constants()
        pair_stage1._validate_constants()
    except (pair_stage0.Stage0Error, pair_stage1.Stage1Error,
            pair_stage1.crossfit.CrossfitError,
            screen.ShortScreenError) as exc:
        fail("imported historical helper changed: {}".format(exc))
    expected = (
        (30, 2, 9, 10, "burst", 354),
        (573, 1000, 0, 11, "burst", 138273),
        (820, 10000, 0, 9, "adversarial", 202780),
    )
    realized = tuple(
        (value, Coordinate(value).K, Coordinate(value).attempt,
         Coordinate(value).root_index, Coordinate(value).schedule,
         Coordinate(value).baseline_record_ordinal)
        for value in BASELINE_COORDINATE_ORDINALS)
    if (realized != expected or
            BASELINE_MATRIX_ORDINALS != tuple(item[-1] for item in expected) or
            pair_stage1.EXPECTED_BASELINE_WEAK_ORDINALS !=
                BASELINE_MATRIX_ORDINALS or
            REPETITIONS != (0, 1) or EXPECTED_INVOCATION_COUNT != 6 or
            CANDIDATE_PAIR != "hashed" or PAIR_SELECTOR_SALT != "0x4d495832"):
        fail("hashed-pair Stage 0 roster or candidate changed")
    if (tuple(CSV_HEADER) != tuple(pair_stage1.CSV_HEADER) or
            tuple(METADATA_ORDER) != tuple(pair_stage1.METADATA_ORDER)):
        fail("imported precodefail protocol changed")
    expected_sources = (
        "bench/wh2_mix2_hashed_pair_stage0.py",
        "bench/test_wh2_mix2_hashed_pair_stage0.py",
    ) + tuple(pair_stage1.SOURCE_PATHS)
    if (SOURCE_PATHS != expected_sources or
            len(SOURCE_PATHS) != len(set(SOURCE_PATHS))):
        fail("hashed-pair linked source roster changed")
    make_invocations()


def _validate_historical_paths(baseline_dir: Path, v9_dir: Path,
                               stage0_dir: Path) -> None:
    if (baseline_dir != EXPECTED_BASELINE_DIR or
            v9_dir != EXPECTED_V9_DIR or
            stage0_dir != EXPECTED_STAGE0_DIR):
        fail("historical directories differ from the three frozen /dev/shm inputs")


def _baseline_projection(row: Mapping[str, Any]) -> Mapping[str, Any]:
    try:
        return pair_stage1._baseline_projection(row)
    except pair_stage1.Stage1Error as exc:
        fail("selected pair01 baseline is malformed: {}".format(exc))
    raise AssertionError("unreachable")


def _validate_selected_baseline(row: Mapping[str, Any],
                                coordinate: Coordinate) -> None:
    expected = coordinate.stage1.identity()
    if (type(row) is not dict or
            row.get("record_ordinal") != coordinate.baseline_record_ordinal or
            row.get("success") is not False or
            row.get("rank_fail") is not True or
            row.get("error") is not False or
            row.get("outcome") != "rank_fail" or
            any(row.get(field) != value for field, value in expected.items()
                if field not in ("coordinate_ordinal",
                                 "baseline_record_ordinal"))):
        fail("selected pair01 weak baseline changed")
    _baseline_projection(row)


def select_baselines(full_baseline: Sequence[Mapping[str, Any]]) \
        -> Tuple[Mapping[str, Any], ...]:
    if len(full_baseline) != pair_stage1.COORDINATE_COUNT:
        fail("authenticated pair01 baseline has the wrong cardinality")
    selected = tuple(full_baseline[value]
                     for value in BASELINE_COORDINATE_ORDINALS)
    for coordinate_ordinal, row in zip(
            BASELINE_COORDINATE_ORDINALS, selected):
        _validate_selected_baseline(row, Coordinate(coordinate_ordinal))
    if tuple(row["record_ordinal"] for row in selected) != \
            BASELINE_MATRIX_ORDINALS:
        fail("selected pair01 matrix ordinals changed")
    return selected


DETERMINISTIC_FIELDS = tuple(pair_stage0.DETERMINISTIC_FIELDS) + (
    "selection_index", "schedule_index", "cell_ordinal",
    "baseline_record_ordinal", "argv_sha256",
    "historical_evidence_sha256", "pair01_baseline_record_sha256",
    "pair01_baseline",
)


def deterministic_projection(record: Mapping[str, Any]) -> Mapping[str, Any]:
    missing = [field for field in DETERMINISTIC_FIELDS if field not in record]
    if missing:
        fail("hashed-pair record omits deterministic fields: {}".format(
            missing))
    return {field: record[field] for field in DETERMINISTIC_FIELDS}


def parse_process_result(result: ProcessResult,
                         source_commit: str,
                         bench_sha256: str,
                         source_receipt_sha256: str,
                         historical_evidence_sha256: str,
                         baseline: Mapping[str, Any]) -> Mapping[str, Any]:
    if (type(result.invocation) is not Invocation or
            SHA256.fullmatch(historical_evidence_sha256) is None):
        fail("hashed-pair process or history identity is malformed")
    _validate_selected_baseline(baseline, result.invocation.coordinate)
    try:
        record = dict(pair_stage0.parse_process_result(
            result, source_commit, bench_sha256, source_receipt_sha256))
    except pair_stage0.Stage0Error as exc:
        fail("candidate process is inadmissible: {}".format(exc))
    record["schema"] = RECORD_SCHEMA
    record.update(result.invocation.identity())
    record["argv_sha256"] = result.command_sha256
    record["raw_stdout_sha256"] = result.stdout_sha256
    record["raw_stderr_sha256"] = result.stderr_sha256
    record["historical_evidence_sha256"] = historical_evidence_sha256
    record["pair01_baseline_record_sha256"] = sha256_json(baseline)
    record["pair01_baseline"] = _baseline_projection(baseline)
    record["deterministic_projection_sha256"] = sha256_json(
        deterministic_projection(record))
    return record


def _run_raw(invocation: Invocation,
             pinned: screen.PinnedExecutable) -> ProcessResult:
    try:
        return pair_stage0._run_raw(invocation, pinned)
    except pair_stage0.Stage0Error as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _execute_roster(pinned: screen.PinnedExecutable,
                    source_commit: str,
                    bench_sha256: str,
                    source_receipt_sha256: str,
                    historical_evidence_sha256: str,
                    baselines: Sequence[Mapping[str, Any]]) \
        -> Tuple[Mapping[str, Any], ...]:
    if len(baselines) != len(BASELINE_COORDINATE_ORDINALS):
        fail("selected baseline roster has the wrong cardinality")
    records: List[Mapping[str, Any]] = []
    for invocation in make_invocations():
        result = _run_raw(invocation, pinned)
        records.append(parse_process_result(
            result, source_commit, bench_sha256, source_receipt_sha256,
            historical_evidence_sha256,
            baselines[invocation.coordinate.selection_index]))
    if len(records) != EXPECTED_INVOCATION_COUNT:
        fail("hashed-pair Stage 0 executed the wrong invocation count")
    return tuple(records)


def adjudicate(records: Sequence[Mapping[str, Any]],
               baselines: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    if (len(records) != EXPECTED_INVOCATION_COUNT or
            len(baselines) != len(BASELINE_COORDINATE_ORDINALS)):
        fail("hashed-pair result or baseline roster has wrong cardinality")
    for coordinate_ordinal, baseline in zip(
            BASELINE_COORDINATE_ORDINALS, baselines):
        _validate_selected_baseline(baseline, Coordinate(coordinate_ordinal))

    by_key: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    for invocation, record in zip(make_invocations(), records):
        if type(record) is not dict:
            fail("hashed-pair record roster or receipt changed")
        identity = invocation.identity()
        key = (invocation.coordinate_ordinal, invocation.repetition)
        terminal_fields = (
            "success", "rank_fail", "error", "configuration_failure")
        terminal_types_valid = all(
            type(record.get(field)) is bool for field in terminal_fields)
        terminal_state_valid = terminal_types_valid and sum(
            int(record[field]) for field in terminal_fields) == 1
        expected_outcome = None
        if terminal_state_valid:
            expected_outcome = next(
                field for field in terminal_fields if record[field])
        if (key in by_key or
                record.get("schema") != RECORD_SCHEMA or
                any(record.get(field) != value for field, value in
                    identity.items()) or
                record.get("promotion_evidence") is not False or
                record.get("fresh_roots_used") is not False or
                record.get("timing_evidence") is not False or
                not terminal_state_valid or
                record.get("outcome") != expected_outcome or
                type(record.get("weak")) is not bool or
                record.get("weak") is not (not record.get("success")) or
                record.get("raw_stdout_sha256") != record.get("stdout_sha256") or
                record.get("raw_stderr_sha256") != record.get("stderr_sha256") or
                record.get("argv_sha256") != record.get("command_sha256") or
                record.get("deterministic_projection_sha256") !=
                    sha256_json(deterministic_projection(record))):
            fail("hashed-pair record roster or receipt changed")
        baseline = baselines[invocation.coordinate.selection_index]
        if (record.get("pair01_baseline_record_sha256") !=
                sha256_json(baseline) or
                record.get("pair01_baseline") !=
                    _baseline_projection(baseline)):
            fail("hashed-pair record is not bound to its pair01 baseline")
        by_key[key] = record
    if len(by_key) != EXPECTED_INVOCATION_COUNT:
        fail("hashed-pair result roster is incomplete")

    successes = 0
    errors = 0
    configuration_failures = 0
    repeat_exact_count = 0
    cap_observations = 0
    coordinate_results: List[Mapping[str, Any]] = []
    work_gates = {
        "block_xors": True,
        "gf256_muladds": True,
        "inactivated_columns": True,
    }
    seed_gates = {
        "effective_precode_seed": True,
        "effective_packet_seed": True,
    }
    for selection, coordinate_ordinal in enumerate(
            BASELINE_COORDINATE_ORDINALS):
        first = by_key[(coordinate_ordinal, 0)]
        second = by_key[(coordinate_ordinal, 1)]
        exact = deterministic_projection(first) == \
            deterministic_projection(second)
        repeat_exact_count += int(exact)
        pair_records = (first, second)
        successes += sum(int(row["success"]) for row in pair_records)
        errors += sum(int(row["error"]) for row in pair_records)
        configuration_failures += sum(
            int(row["configuration_failure"]) for row in pair_records)
        cap_observations += sum(int(
            type(row["inactivated_columns"]) is int and
            row["inactivated_columns"] <= INACTIVATION_CAP)
            for row in pair_records)

        baseline = baselines[selection]
        comparisons: Dict[str, Any] = {}
        for field in work_gates:
            candidate_value = first[field]
            baseline_value = baseline[field]
            passed = (type(candidate_value) is int and
                      type(baseline_value) is int and
                      candidate_value <= baseline_value)
            work_gates[field] = work_gates[field] and passed
            comparisons[field] = {
                "candidate": candidate_value,
                "pair01_baseline": baseline_value,
                "candidate_minus_baseline": (
                    candidate_value - baseline_value
                    if type(candidate_value) is int and
                    type(baseline_value) is int else None),
                "not_above_baseline": passed,
            }
        seed_comparisons: Dict[str, Any] = {}
        for field in seed_gates:
            candidate_value = first[field]
            baseline_value = baseline[field]
            passed = candidate_value == baseline_value
            seed_gates[field] = seed_gates[field] and passed
            seed_comparisons[field] = {
                "candidate": candidate_value,
                "pair01_baseline": baseline_value,
                "equal": passed,
            }
        coordinate_results.append({
            **Coordinate(coordinate_ordinal).identity(),
            "both_success": first["success"] and second["success"],
            "repeated_deterministic_projection_exact": exact,
            "inactivation_cap_both": all(
                type(row["inactivated_columns"]) is int and
                row["inactivated_columns"] <= INACTIVATION_CAP
                for row in pair_records),
            "candidate_vs_pair01": comparisons,
            "effective_seed_identity": seed_comparisons,
        })

    gates = {
        "all_six_observations_succeed":
            successes == EXPECTED_INVOCATION_COUNT,
        "candidate_errors_zero": errors == 0,
        "candidate_configuration_failures_zero":
            configuration_failures == 0,
        "all_three_repeated_deterministic_projections_exact":
            repeat_exact_count == len(BASELINE_COORDINATE_ORDINALS),
        "all_six_inactivated_columns_at_most_1024":
            cap_observations == EXPECTED_INVOCATION_COUNT,
        "per_coordinate_block_xors_not_above_pair01":
            work_gates["block_xors"],
        "per_coordinate_gf256_muladds_not_above_pair01":
            work_gates["gf256_muladds"],
        "per_coordinate_inactivated_columns_not_above_pair01":
            work_gates["inactivated_columns"],
        "per_coordinate_effective_precode_seed_equals_pair01":
            seed_gates["effective_precode_seed"],
        "per_coordinate_effective_packet_seed_equals_pair01":
            seed_gates["effective_packet_seed"],
    }
    return {
        "candidate_pair": CANDIDATE_PAIR,
        "candidate_observation_count": len(records),
        "candidate_successes": successes,
        "candidate_errors": errors,
        "candidate_configuration_failures": configuration_failures,
        "repeat_exact_coordinate_count": repeat_exact_count,
        "inactivation_cap_observation_count": cap_observations,
        "coordinates": coordinate_results,
        "gates": gates,
        "disposition": "PASS" if all(gates.values()) else "REJECT",
    }


def _source_receipt() -> Mapping[str, Any]:
    for item in SOURCE_PATHS:
        if not (REPO_ROOT / item).is_file():
            fail("required hashed-pair source is missing: {}".format(item))
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
        fail("cannot attest hashed-pair source checkout: {}".format(exc))
    if (head_process.returncode != 0 or head_process.stderr or
            status_process.returncode != 0 or status_process.stderr or
            tracked_process.returncode != 0 or tracked_process.stderr):
        fail("cannot attest hashed-pair source checkout")
    try:
        head = head_process.stdout.decode("ascii").strip()
        tracked = tracked_process.stdout.decode("ascii").split("\0")
    except UnicodeDecodeError:
        fail("hashed-pair source receipt is not ASCII")
    if COMMIT.fullmatch(head) is None or status_process.stdout:
        fail("hashed-pair linked sources are not clean at a full source HEAD")
    if tracked and tracked[-1] == "":
        tracked.pop()
    if len(tracked) != len(SOURCE_PATHS) or set(tracked) != set(SOURCE_PATHS):
        fail("hashed-pair source receipt includes a path not tracked at HEAD")
    hashes: Dict[str, str] = {}
    for item in SOURCE_PATHS:
        working = screen._stable_file_bytes(REPO_ROOT / item)
        try:
            blob_process = subprocess.run(
                ["git", "cat-file", "blob", "{}:{}".format(head, item)],
                cwd=str(REPO_ROOT), stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False, timeout=30)
        except (OSError, subprocess.SubprocessError) as exc:
            fail("cannot read HEAD blob for {}: {}".format(item, exc))
        if (blob_process.returncode != 0 or blob_process.stderr or
                blob_process.stdout != working):
            fail("hashed-pair working source differs from HEAD: {}".format(
                item))
        hashes[item] = sha256_bytes(working)
    try:
        final_head_process = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot recheck hashed-pair source HEAD: {}".format(exc))
    if (final_head_process.returncode != 0 or final_head_process.stderr or
            final_head_process.stdout != (head + "\n").encode("ascii")):
        fail("hashed-pair source HEAD changed during attestation")
    receipt: Dict[str, Any] = {
        "source_git_commit": head,
        "tracked_and_untracked_linked_sources_clean": True,
        "all_source_paths_tracked_at_HEAD": True,
        "working_bytes_equal_HEAD_blobs": True,
        "final_HEAD_rechecked": True,
        "clean_status_scope": list(SOURCE_STATUS_PATHS),
        "source_files": hashes,
    }
    receipt["source_receipt_sha256"] = self_hash(
        receipt, "source_receipt_sha256")
    return receipt


def run_stage0(bench: Path, baseline_dir: Path, v9_dir: Path,
               stage0_dir: Path, output_dir: Path) -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    _validate_historical_paths(baseline_dir, v9_dir, stage0_dir)
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path must be one fresh path below an existing directory")
    source = _source_receipt()
    # Complete historical authentication and selection before opening a path to
    # any candidate precodefail invocation.
    try:
        historical_context = pair_stage1.load_historical_evidence(
            baseline_dir, v9_dir, stage0_dir)
    except (pair_stage1.Stage1Error, pair_stage0.Stage0Error,
            pair_stage1.crossfit.CrossfitError,
            screen.ShortScreenError) as exc:
        fail("historical authentication failed: {}".format(exc))
    with historical_context as historical:
        baselines = select_baselines(historical.baseline)
        historical_receipt = historical.receipt()
        historical_evidence_sha256 = sha256_json(historical_receipt)
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            return _run_stage0_pinned(
                pinned, output_dir, source, contract, historical,
                historical_receipt, historical_evidence_sha256, baselines)


def _run_stage0_pinned(
        pinned: screen.PinnedExecutable,
        output_dir: Path,
        source: Mapping[str, Any],
        contract: Mapping[str, Any],
        historical: pair_stage1.HistoricalEvidence,
        historical_receipt: Mapping[str, Any],
        historical_evidence_sha256: str,
        baselines: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    source_commit = source["source_git_commit"]
    source_receipt_sha256 = source["source_receipt_sha256"]
    bench_receipt = pinned.receipt()
    # Deliberately the final preflight before the first candidate process.
    description = screen.read_bench_description(
        pinned.path, source_commit, pinned.descriptor)
    records = _execute_roster(
        pinned, source_commit, bench_receipt["sha256"],
        source_receipt_sha256, historical_evidence_sha256, baselines)
    verdict = adjudicate(records, baselines)

    try:
        temporary = Path(tempfile.mkdtemp(
            prefix=".wh2-mix2-hashed-pair-stage0-",
            dir=str(output_dir.parent)))
    except OSError as exc:
        fail("cannot create hashed-pair staging directory: {}".format(exc))
    published = False
    try:
        try:
            records_sha256, record_count = pair_stage0._write_jsonl(
                temporary / RECORD_NAME, records)
        except (pair_stage0.Stage0Error, OSError) as exc:
            fail("cannot durably stage hashed-pair records: {}".format(exc))
        if record_count != EXPECTED_INVOCATION_COUNT:
            fail("hashed-pair JSONL cardinality changed")
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                historical.receipt() != historical_receipt or
                sha256_json(historical.receipt()) !=
                    historical_evidence_sha256):
            fail("hashed-pair executable, source, or history changed during run")
        summary: Dict[str, Any] = {
            "schema": SUMMARY_SCHEMA,
            "contract_sha256": contract["contract_sha256"],
            "source_receipt": source,
            "bench": bench_receipt,
            "bench_description": description,
            "historical_evidence": historical_receipt,
            "historical_evidence_sha256": historical_evidence_sha256,
            "selected_pair01_baselines": [
                {
                    "coordinate": Coordinate(value).identity(),
                    "record_sha256": sha256_json(baseline),
                    "projection": _baseline_projection(baseline),
                }
                for value, baseline in zip(
                    BASELINE_COORDINATE_ORDINALS, baselines)
            ],
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
        try:
            screen._write_exclusive(
                temporary / SUMMARY_NAME,
                (canonical_json(summary) + "\n").encode("ascii"))
        except (screen.ShortScreenError, OSError) as exc:
            fail("cannot durably stage hashed-pair summary: {}".format(exc))
        if (pinned.receipt() != bench_receipt or
                _source_receipt() != source or
                historical.receipt() != historical_receipt or
                sha256_json(historical.receipt()) !=
                    historical_evidence_sha256):
            fail("hashed-pair inputs changed before publication")
        try:
            pair_stage0._fsync_directory(temporary)
        except (pair_stage0.Stage0Error, OSError) as exc:
            fail("cannot sync hashed-pair staging directory: {}".format(exc))
        try:
            screen._publish_directory_noreplace(temporary, output_dir)
        except (screen.ShortScreenError, OSError) as exc:
            fail("cannot publish hashed-pair output: {}".format(exc))
        published = True
        try:
            pair_stage0._fsync_directory(output_dir.parent)
        except (pair_stage0.Stage0Error, OSError) as exc:
            fail(
                "hashed-pair output was published, but its parent directory "
                "could not be synced; treat the output as invalid: {}".format(
                    exc))
        return summary
    finally:
        active_error = sys.exc_info()[1]
        if not published and temporary.exists():
            try:
                shutil.rmtree(temporary)
            except OSError as exc:
                message = "cannot remove hashed-pair staging directory: {}" \
                    .format(exc)
                if active_error is None:
                    fail(message)
                if hasattr(active_error, "add_note"):
                    active_error.add_note(message)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("describe", help="print the frozen Stage 0 contract")
    run_parser = subparsers.add_parser(
        "run", help="run the six-observation hashed-pair Stage 0")
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
        summary = run_stage0(
            arguments.bench, arguments.baseline_dir, arguments.v9_dir,
            arguments.stage0_dir, arguments.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (HashedPairStage0Error, pair_stage0.Stage0Error,
            pair_stage1.Stage1Error, pair_stage1.crossfit.CrossfitError,
            screen.ShortScreenError) as exc:
        print("wh2 MIX2 hashed-pair Stage 0: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
