#!/usr/bin/env python3
"""Fail-closed tests for the minimal Wirehair2 benchmark contract."""

from __future__ import annotations

import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as subject


EXPECTED_RECOVERY_HASHES = {
    "development": "f97f28c211428cd77aed97160073b192d93014cb4a61a844bc7d76375ac61b77",
    "final_raw": "c33160c0379675dac9de3162463e3100041bc8eeb55c002dcf6c3473973e7749",
    "final_repaired": "e247b0ff51f47383f415b1747a783e0c4b4a91e97862e9e72aa660241ff0aeee",
    "final_validation": "a63611ad300396a1b8977281fcd387fde4bd207a3ca82043c5e234e8fc39d4de",
    "cross_width_validation": "80fb4c19b31f88101ef4b3c480e6d5ff4a3bba7d45c96f02dc3518cd7c1c0399",
}
EXPECTED_TIMING_HASHES = {
    "development": "1a15ee48e893280013b74f448d81607e45630fbf631d4e73960b3a2194e204c3",
    "final": "90bd3ea12b21b5a04290d7e55aac95287f9d0fe3150ff27047027dfe84ba63ec",
}


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


Outcome = Tuple[str, Any]
OutcomeFunction = Callable[[str, int, Mapping[str, Any]], Outcome]


def successful(_arm: str, _ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
    return "success", 0


def ledger_rows(
        contract: Mapping[str, Any],
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        outcome: OutcomeFunction = successful,
        trace_hashes: Sequence[str] = (),
        phase: str = "development",
        repair_attempts: Optional[Mapping[str, Mapping[int, int]]] = None,
        repair_map_hashes: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, Any]]:
    if repair_attempts is None:
        repair_attempts = {}
    if repair_map_hashes is None:
        repair_map_hashes = {}
    rows: List[Dict[str, Any]] = []
    for ordinal, cell in enumerate(subject.iter_recovery_cells(contract, phase)):
        trace = trace_hashes[ordinal] if trace_hashes else digest(
            "trace:" + subject.sha256_json(cell))
        for arm in arms:
            result, extra = outcome(arm, ordinal, cell)
            descriptor = digest("descriptor:" + arm)
            codec = "wirehair2_certified" if arm == "wirehair2_head" else \
                "wirehair1" if arm == "wirehair1" else "wirehair2_experiment"
            attempt = repair_attempts.get(arm, {}).get(
                cell["K"],
                0 if codec == "wirehair1" else cell["base_seed_attempt"])
            rows.append({
                "arm": arm,
                **cell,
                "outcome": result,
                "decoded_extra": extra,
                "cell_sha256": subject.sha256_json(cell),
                "trace_sha256": trace,
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": descriptor,
                "construction_attempt": attempt,
                "realized_construction_sha256":
                    subject.generic_realized_construction_sha256(
                        codec, descriptor, cell["K"], cell["block_bytes"],
                        attempt),
                "repair_map_sha256": repair_map_hashes.get(arm, "0" * 64),
            })
    return rows


def write_ledger(root: Path, rows: Sequence[Mapping[str, Any]],
                 *, terminal_newline: bool = True) -> Path:
    path = root / "ledger.jsonl"
    data = "\n".join(subject.canonical_json(row) for row in rows)
    if terminal_newline:
        data += "\n"
    path.write_bytes(data.encode("utf-8"))
    return path


def trace_rows(contract: Mapping[str, Any],
               trace_hashes: Sequence[str] = (),
               phase: str = "development") -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for ordinal, cell in enumerate(subject.iter_recovery_cells(
            contract, phase)):
        trace = trace_hashes[ordinal] if trace_hashes else digest(
            "trace:" + subject.sha256_json(cell))
        rows.append({
            "ordinal": ordinal,
            "cell_sha256": subject.sha256_json(cell),
            "trace_sha256": trace,
        })
    return rows


def write_trace_manifest(
        root: Path, contract: Mapping[str, Any],
        trace_hashes: Sequence[str] = (),
        phase: str = "development") -> Path:
    path = root / "traces.jsonl"
    data = "".join(
        subject.canonical_json(row) + "\n"
        for row in trace_rows(contract, trace_hashes, phase))
    path.write_bytes(data.encode("utf-8"))
    return path


def write_freeze_manifest(
        root: Path, contract: Mapping[str, Any], trace_path: Path,
        arms: Sequence[str], phase: str = "development",
        repair_map_hashes: Optional[Mapping[str, str]] = None,
        training_trace_sha256: str = "0" * 64) -> Path:
    if repair_map_hashes is None:
        repair_map_hashes = {}
    arm_records = []
    for arm in arms:
        codec = "wirehair2_certified" if arm == "wirehair2_head" else \
            "wirehair1" if arm == "wirehair1" else "wirehair2_experiment"
        mapped = arm in repair_map_hashes
        arm_records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256": digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else
                "repair_map" if mapped else "raw_base",
            "repair_map_sha256": repair_map_hashes.get(arm, "0" * 64),
        })
    manifest = {
        "schema": subject.FREEZE_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": phase,
        "domain_sha256": contract["recovery"]["domains"][phase][
            "domain_sha256"],
        "source_git_commit": "1" * 40,
        "arm_roster": list(arms),
        "arm_roster_sha256": subject.arm_roster_sha256(arms),
        "trace_manifest_sha256": subject.trace_manifest_sha256(
            contract, "recovery", phase, trace_path),
        "repair_training_trace_manifest_sha256": training_trace_sha256,
        "commands": [["wirehair_v2_bench", "raw-cell", "--frozen"]],
        "cpu_affinity": [0, 1],
        "host_identity": {"name": "unit-test-host"},
        "arms": arm_records,
    }
    path = root / "freeze.json"
    path.write_bytes((subject.canonical_json(manifest) + "\n").encode("utf-8"))
    return path


def timing_trace_hashes(contract: Mapping[str, Any]) -> List[str]:
    return [
        digest("timing-trace:" + subject.sha256_json(cell))
        for cell in subject.iter_timing_cells(contract, "development")
    ]


def write_timing_trace_manifest(
        root: Path, contract: Mapping[str, Any],
        hashes: Sequence[str]) -> Path:
    rows = []
    for ordinal, cell in enumerate(subject.iter_timing_cells(
            contract, "development")):
        rows.append({
            "ordinal": ordinal,
            "cell_sha256": subject.sha256_json(cell),
            "trace_sha256": hashes[ordinal],
        })
    path = root / "timing-traces.jsonl"
    path.write_bytes("".join(
        subject.canonical_json(row) + "\n" for row in rows).encode("utf-8"))
    return path


def arm_codec(arm: str) -> str:
    if arm == "wirehair2_head":
        return "wirehair2_certified"
    if arm == "wirehair1":
        return "wirehair1"
    return "wirehair2_experiment"


def write_timing_freeze_manifest(
        root: Path, contract: Mapping[str, Any], trace_path: Path,
        arms: Sequence[str]) -> Path:
    arm_records = []
    for arm in arms:
        codec = arm_codec(arm)
        arm_records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256": digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else "raw_base",
            "repair_map_sha256": "0" * 64,
        })
    manifest = {
        "schema": subject.FREEZE_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "evidence_kind": "timing",
        "phase": "development",
        "domain_sha256": contract["timing"]["domains"]["development"][
            "domain_sha256"],
        "source_git_commit": "1" * 40,
        "arm_roster": list(arms),
        "arm_roster_sha256": subject.arm_roster_sha256(arms),
        "trace_manifest_sha256": subject.trace_manifest_sha256(
            contract, "timing", "development", trace_path),
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_v2_bench", "timing-cell", "--frozen"]],
        "cpu_affinity": [0, 1],
        "host_identity": {"name": "unit-test-host"},
        "arms": arm_records,
    }
    path = root / "timing-freeze.json"
    path.write_bytes((subject.canonical_json(manifest) + "\n").encode("utf-8"))
    return path


def write_final_freeze_manifests(
        root: Path, contract: Mapping[str, Any],
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        ) -> Dict[Tuple[str, str], Path]:
    training_trace = digest("final repaired training traces")
    map_hashes = {
        arm: digest("production map:" + arm)
        for arm in arms if arm != "wirehair1"
    }
    paths: Dict[Tuple[str, str], Path] = {}
    keys = (
        ("recovery", "final_raw"),
        ("recovery", "final_repaired"),
        ("recovery", "final_validation"),
        ("recovery", "cross_width_validation"),
        ("timing", "final"),
    )
    for evidence_kind, phase in keys:
        raw = phase == "final_raw"
        records = []
        for arm in arms:
            codec = arm_codec(arm)
            records.append({
                "arm": arm,
                "codec": codec,
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": digest("descriptor:" + arm),
                "construction_policy":
                    "not_applicable" if codec == "wirehair1" else
                    "raw_base" if raw else "repair_map",
                "repair_map_sha256":
                    "0" * 64 if raw or codec == "wirehair1" else
                    map_hashes[arm],
            })
        trace_sha256 = training_trace if phase == "final_repaired" else \
            digest("trace:{}:{}".format(evidence_kind, phase))
        manifest = {
            "schema": subject.FREEZE_SCHEMA,
            "contract_sha256": subject.contract_sha256(contract),
            "evidence_kind": evidence_kind,
            "phase": phase,
            "domain_sha256": contract[evidence_kind]["domains"][phase][
                "domain_sha256"],
            "source_git_commit": "1" * 40,
            "arm_roster": list(arms),
            "arm_roster_sha256": subject.arm_roster_sha256(arms),
            "trace_manifest_sha256": trace_sha256,
            "repair_training_trace_manifest_sha256":
                "0" * 64 if raw else training_trace,
            "commands": [["wirehair_v2_bench", evidence_kind, phase]],
            "cpu_affinity": [0, 1],
            "host_identity": {"name": "unit-test-host"},
            "arms": records,
        }
        path = root / (evidence_kind + "-" + phase + "-freeze.json")
        path.write_bytes(
            (subject.canonical_json(manifest) + "\n").encode("utf-8"))
        paths[(evidence_kind, phase)] = path
    return paths


def architecture_selection_receipt(
        contract: Mapping[str, Any], selected: str = "candidate",
        ) -> Dict[str, Any]:
    selected_artifact = {
        "codec": "wirehair2_experiment",
        "arm_descriptor_sha256": digest("descriptor:" + selected),
    }
    result = {
        "schema": subject.SCHEMA + ".architecture-selection.v1",
        "contract_sha256": subject.contract_sha256(contract),
        "recovery_domain_sha256": contract["recovery"]["domains"]
            ["development"]["domain_sha256"],
        "timing_domain_sha256": contract["timing"]["domains"]
            ["development"]["domain_sha256"],
        "recovery_freeze_manifest_sha256": digest("recovery freeze"),
        "timing_freeze_manifest_sha256": digest("timing freeze"),
        "architecture_artifact_sha256": digest("architecture artifacts"),
        "recovery_cells_per_arm": 360,
        "timing_rows": 2112,
        "candidate_roster": [selected],
        "eligible_candidates": [selected],
        "eligible_overhead0_failures": {selected: 0},
        "minimum_overhead0_failures": 0,
        "recovery_equivalence_allowance": 1,
        "recovery_equivalent_candidates": [selected],
        "ranking": [{
            "arm": selected,
            "decoder_solve_mean_log_ratio": -0.1,
            "failures_overhead0": 0,
            "failures_overhead1": 0,
            "failures_overhead2": 0,
            "failures_overhead4": 0,
        }],
        "selected_arm": selected,
        "selected_codec": selected_artifact["codec"],
        "selected_arm_descriptor_sha256":
            selected_artifact["arm_descriptor_sha256"],
        "selected_architecture_sha256":
            subject.selected_architecture_sha256(selected, selected_artifact),
    }
    result["selection_sha256"] = subject.sha256_json(result)
    return result


def timing_receipt_rows(
        contract: Mapping[str, Any], hashes: Sequence[str],
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        candidate_scale: float = 0.8) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    panels = subject.timing_panels(contract, arms)
    for ordinal, cell in enumerate(subject.iter_timing_cells(
            contract, "development")):
        for panel in panels:
            left = panel["left_arm"]
            right = panel["right_arm"]
            order = subject.timing_order(panel, cell["replicate"])
            if panel["panel_kind"] == "AA":
                left_ns = right_ns = 100000
            else:
                left_ns = int(100000 * candidate_scale)
                right_ns = 100000
            elapsed = [left_ns, right_ns, right_ns, left_ns] \
                if order == "ABBA" else \
                [right_ns, left_ns, left_ns, right_ns]
            row: Dict[str, Any] = {
                **cell,
                **panel,
                "order": order,
                "left_outcome": "success",
                "right_outcome": "success",
                "left_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else 4,
                "right_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else 4,
                "elapsed_ns": elapsed,
                "cell_sha256": subject.sha256_json(cell),
                "trace_sha256": hashes[ordinal],
            }
            for side, arm in (("left", left), ("right", right)):
                codec = arm_codec(arm)
                descriptor = digest("descriptor:" + arm)
                attempt = 0 if codec == "wirehair1" else \
                    cell["base_seed_attempt"]
                row[side + "_binary_sha256"] = digest("binary:" + arm)
                row[side + "_arm_descriptor_sha256"] = descriptor
                row[side + "_construction_attempt"] = attempt
                row[side + "_realized_construction_sha256"] = \
                    subject.generic_realized_construction_sha256(
                        codec, descriptor, cell["K"], cell["block_bytes"],
                        attempt)
                row[side + "_repair_map_sha256"] = "0" * 64
            rows.append(row)
    return rows


def write_timing_receipt(root: Path,
                         rows: Sequence[Mapping[str, Any]]) -> Path:
    path = root / "timing.jsonl"
    path.write_bytes("".join(
        subject.canonical_json(row) + "\n" for row in rows).encode("utf-8"))
    return path


def write_repair_map(
        root: Path, contract: Mapping[str, Any], arm: str,
        training_trace_sha256: str,
        retry_overrides: Optional[Mapping[int, int]] = None,
        ) -> Tuple[Path, str, Dict[int, int]]:
    if retry_overrides is None:
        retry_overrides = {}
    offsets = [retry_overrides.get(K, 0) for K in range(2, 64001)]
    value = {
        "schema": subject.REPAIR_MAP_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "training_domain_sha256": contract["recovery"]["domains"][
            "final_repaired"]["domain_sha256"],
        "training_trace_manifest_sha256": training_trace_sha256,
        "arm": arm,
        "source_git_commit": "1" * 40,
        "binary_sha256": digest("binary:" + arm),
        "arm_descriptor_sha256": digest("descriptor:" + arm),
        "production_base_seed_attempt": 0,
        "entry_kind": "retry_offset_indexed_by_K_minus_2",
        "attempt_derivation": subject.EXPECTED_ATTEMPT_DERIVATION,
        "repair_rule": subject.EXPECTED_REPAIR_RULE,
        "retry_offsets": offsets,
    }
    path = root / (arm + "-repair-map.json")
    path.write_bytes((subject.canonical_json(value) + "\n").encode("utf-8"))
    return path, subject.repair_map_sha256(value), {
        K: (contract["seeds"]["production_base_seed_attempt"] + offset) & 0xff
        for K, offset in zip(range(2, 64001), offsets)
    }


class ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()

    def test_exact_recovery_domains_and_hashes(self) -> None:
        expected_counts = {
            "development": 360,
            "final_raw": 575991,
            "final_repaired": 575991,
            "final_validation": 575991,
            "cross_width_validation": 1440,
        }
        for phase, expected_hash in EXPECTED_RECOVERY_HASHES.items():
            with self.subTest(phase=phase):
                domain = self.contract["recovery"]["domains"][phase]
                self.assertEqual(domain["expected_cells_per_arm"], expected_counts[phase])
                self.assertEqual(domain["domain_sha256"], expected_hash)
                self.assertEqual(subject.recovery_domain_sha256(
                    self.contract, phase), expected_hash)

    def test_exact_timing_domains_and_hashes(self) -> None:
        for phase, expected_hash in EXPECTED_TIMING_HASHES.items():
            with self.subTest(phase=phase):
                domain = self.contract["timing"]["domains"][phase]
                self.assertEqual(subject.timing_domain_sha256(
                    self.contract, phase), expected_hash)
                self.assertEqual(domain["domain_sha256"], expected_hash)
        self.assertEqual(
            self.contract["timing"]["domains"]["development"]["expected_cells"],
            192,
        )
        self.assertEqual(
            self.contract["timing"]["domains"]["final"]["expected_cells"],
            3600,
        )

    def test_every_k_hard_domain_is_575991_not_a_cross_product(self) -> None:
        cells = self.contract["recovery"]["domains"]["final_raw"][
            "expected_cells_per_arm"]
        self.assertEqual(cells, 63999 * 3 * 3)
        first = next(subject.iter_recovery_cells(self.contract, "final_raw"))
        self.assertEqual(first["K"], 2)
        self.assertEqual(first["base_seed_attempt"],
                         self.contract["seeds"]["raw_base_seed_attempts"][0])
        self.assertEqual(first["loss_seed"],
                         self.contract["seeds"]["training_loss_roots"][0])

    def test_v1_bands_and_short_cohorts_are_immutable(self) -> None:
        for mutation in ("band", "short", "timing"):
            changed = copy.deepcopy(self.contract)
            if mutation == "band":
                changed["k_bands"][0]["last"] = 99
                changed["k_bands"][1]["first"] = 100
            elif mutation == "short":
                changed["k_sets"]["short"].insert(9, 99)
            else:
                changed["k_sets"]["timing_short"][0] = 7
            with self.subTest(mutation=mutation), \
                    self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_mutated_domain_hash_rejects(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["recovery"]["domains"]["development"]["domain_sha256"] = "0" * 64
        with self.assertRaises(subject.ContractError):
            subject._validate_structure(changed, check_domain_hashes=True)

    def test_v1_semantics_cannot_be_resealed_to_an_easier_contract(self) -> None:
        changes = []
        easy_loss = copy.deepcopy(self.contract)
        easy_loss["strata_sets"]["development"][1] = {
            "schedule": "easy", "loss_ppm": 1}
        easy_loss["strata_sets"]["hard"][0] = {
            "schedule": "easy", "loss_ppm": 1}
        changes.append(("loss", easy_loss))
        short_final = copy.deepcopy(self.contract)
        short_final["recovery"]["domains"]["final_raw"]["k_set"] = "short"
        short_final["recovery"]["domains"]["final_raw"][
            "expected_cells_per_arm"] = 270
        changes.append(("domain", short_final))
        easy_timing = copy.deepcopy(self.contract)
        easy_timing["timing"]["schedule"] = "made-up"
        easy_timing["timing"]["loss_ppm"] = -7
        changes.append(("timing", easy_timing))
        for name, changed in changes:
            with self.subTest(name=name), self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_bool_float_numeric_aliases_reject_as_contract_errors(self) -> None:
        mutations = (
            ("goal", ("goal", "primary_overhead"), False),
            ("range", ("k_sets", "all", "first"), 2.0),
            ("count", ("recovery", "domains", "development",
                       "expected_cells_per_arm"), 360.0),
        )
        for name, path, replacement in mutations:
            changed = copy.deepcopy(self.contract)
            target = changed
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = replacement
            with self.subTest(name=name), self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_malformed_roots_and_extreme_json_numbers_fail_closed(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["seeds"]["training_loss_roots"][0] = {}
        with self.assertRaises(subject.ContractError):
            subject._validate_structure(changed, check_domain_hashes=False)
        for payload in (b"1e400", b"1" + b"0" * 5000):
            with self.subTest(payload_bytes=len(payload)), \
                    self.assertRaises(subject.ContractError):
                subject._load_json_bytes(payload, "mutant")

    def test_loss_roots_are_disjoint_and_raw_attempts_are_explicit(self) -> None:
        seeds = self.contract["seeds"]
        self.assertIn(seeds["production_base_seed_attempt"],
                      seeds["raw_base_seed_attempts"])
        self.assertTrue(set(seeds["training_loss_roots"]).isdisjoint(
            seeds["validation_loss_roots"]))
        self.assertEqual(self.contract["recovery"]["raw_construction_attempts"], 1)
        self.assertEqual(
            self.contract["recovery"]["production_seed_fixups_in_raw_phase"], 0)

    def test_certified_realized_profile_receipt_is_canonical_little_endian(
            self) -> None:
        profile = subject.serialized_wh2_profile(2, 64, 7)
        self.assertEqual(len(profile), 32)
        self.assertEqual(profile[:8], b"WHV2\x01\x00\x20\x00")
        self.assertEqual(profile[8:16], bytes.fromhex("c9f9f447bb5b294b"))
        self.assertEqual(profile[16:24], (128).to_bytes(8, "little"))
        self.assertEqual(profile[24:29], b"\x40\x00\x00\x00\x07")
        self.assertEqual(profile[29:], b"\0\0\0")

    def test_final_freezes_bind_one_architecture_and_production_map(self) -> None:
        mutations = ("binary", "map", "training_trace")
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                paths = write_final_freeze_manifests(root, self.contract)
                summary = subject.validate_final_freeze_continuity(
                    self.contract, paths,
                    architecture_selection_receipt(self.contract))
                self.assertEqual(summary["selected_candidate"], "candidate")
                path = paths[("recovery", "final_validation")]
                value = json.loads(path.read_text(encoding="utf-8"))
                candidate = next(
                    arm for arm in value["arms"]
                    if arm["arm"] == "candidate")
                if mutation == "binary":
                    candidate["binary_sha256"] = digest("substitute binary")
                elif mutation == "map":
                    candidate["repair_map_sha256"] = digest("substitute map")
                else:
                    value["repair_training_trace_manifest_sha256"] = \
                        digest("substitute training trace")
                path.write_bytes(
                    (subject.canonical_json(value) + "\n").encode("utf-8"))
                with self.assertRaises(subject.ContractError):
                    subject.validate_final_freeze_continuity(
                        self.contract, paths,
                        architecture_selection_receipt(self.contract))

        with tempfile.TemporaryDirectory() as temporary:
            paths = write_final_freeze_manifests(Path(temporary), self.contract)
            with self.assertRaises(subject.ContractError):
                subject.validate_final_freeze_continuity(
                    self.contract, paths,
                    architecture_selection_receipt(self.contract, "substitute"))
            substituted = architecture_selection_receipt(self.contract)
            substituted["selected_arm_descriptor_sha256"] = \
                digest("different selected equations")
            substituted["selected_architecture_sha256"] = \
                subject.selected_architecture_sha256("candidate", {
                    "codec": substituted["selected_codec"],
                    "arm_descriptor_sha256":
                        substituted["selected_arm_descriptor_sha256"],
                })
            unsigned = {key: value for key, value in substituted.items()
                        if key != "selection_sha256"}
            substituted["selection_sha256"] = subject.sha256_json(unsigned)
            with self.assertRaises(subject.ContractError):
                subject.validate_final_freeze_continuity(
                    self.contract, paths, substituted)

    def test_reserved_codec_kinds_are_bijective_in_final_freezes(self) -> None:
        for codec, policy, map_hash in (
                ("wirehair1", "not_applicable", "0" * 64),
                ("wirehair2_certified", "repair_map",
                 digest("production map:candidate"))):
            with self.subTest(codec=codec), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                paths = write_final_freeze_manifests(root, self.contract)
                path = paths[("timing", "final")]
                value = json.loads(path.read_text(encoding="utf-8"))
                candidate = next(
                    arm for arm in value["arms"]
                    if arm["arm"] == "candidate")
                candidate["codec"] = codec
                candidate["construction_policy"] = policy
                candidate["repair_map_sha256"] = map_hash
                path.write_bytes(
                    (subject.canonical_json(value) + "\n").encode("utf-8"))
                with self.assertRaises(subject.ContractError):
                    subject.load_freeze_manifest(
                        self.contract, "final", path, "timing")


class LedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.trace_hashes = [
            digest("trace:" + subject.sha256_json(cell))
            for cell in subject.iter_recovery_cells(cls.contract, "development")
        ]
        cls.valid_rows = ledger_rows(
            cls.contract, trace_hashes=cls.trace_hashes)

    def validate(self, rows: Sequence[Mapping[str, Any]],
                 arms: Sequence[str] = (
                     "wirehair2_head", "wirehair1", "candidate"),
                 trace_hashes: Sequence[str] = ()) -> Mapping[str, Any]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = write_ledger(root, rows)
            trace_path = write_trace_manifest(
                root, self.contract, trace_hashes or self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path, arms)
            return subject.validate_ledger(
                self.contract, "development", path, freeze_path, trace_path)

    def assert_rejects(self, rows: Sequence[Mapping[str, Any]],
                       arms: Sequence[str] = (
                           "wirehair2_head", "wirehair1", "candidate")) -> None:
        with self.assertRaises(subject.ContractError):
            self.validate(rows, arms)

    def test_valid_common_ledger_reports_repairs_introductions_and_tails(self) -> None:
        def outcomes(arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair2_head" and ordinal in (0, 1):
                return "need_more_at_cap", None
            if arm == "candidate" and ordinal in (1, 2):
                return "construct_failed", None
            if ordinal == 3:
                return "success", 2
            return "success", 0

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        self.assertEqual(summary["excluded_cells"], 0)
        self.assertEqual(summary["arms"]["wirehair2_head"]["failure_by_overhead"],
                         {"0": 3, "1": 3, "2": 2, "4": 2})
        self.assertEqual(summary["arms"]["candidate"]["failure_by_overhead"],
                         {"0": 3, "1": 3, "2": 2, "4": 2})
        comparison = summary["comparisons"]["candidate"]["controls"][
            "wirehair2_head"]
        self.assertEqual(comparison["overhead0_repairs"], 1)
        self.assertEqual(comparison["overhead0_introductions"], 1)
        self.assertEqual(comparison["overhead0_shared_failures"], 2)
        self.assertTrue(comparison["all_noninferiority_gates_pass"])
        self.assertIn("wirehair1",
                      summary["comparisons"]["candidate"]["controls"])

    def test_explicit_unsupported_rows_stay_in_denominator_and_charge_cap(self) -> None:
        def outcomes(_arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            return ("unsupported", None) if ordinal < 40 else ("success", 0)

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        for arm in ("wirehair2_head", "wirehair1", "candidate"):
            self.assertEqual(summary["arms"][arm]["cells"], 360)
            self.assertEqual(summary["arms"][arm]["status_counts"]["unsupported"], 40)
            self.assertEqual(summary["arms"][arm]["capped_overhead_sum"], 160)
        self.assertFalse(
            summary["comparisons"]["candidate"]["all_controls_noninferior"])

    def test_mandatory_control_unsupported_makes_candidates_ineligible(self) -> None:
        def outcomes(arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair1" and ordinal == 0:
                return "unsupported", None
            return "success", 0

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        self.assertFalse(summary["mandatory_controls_supported"])
        self.assertFalse(summary["comparisons"]["candidate"][
            "architecture_eligible"])

        def construct_failure(
                arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair1" and ordinal == 0:
                return "construct_failed", None
            return "success", 0

        scored = self.validate(ledger_rows(
            self.contract, outcome=construct_failure))
        self.assertTrue(scored["mandatory_controls_supported"])
        self.assertTrue(scored["comparisons"]["candidate"][
            "architecture_eligible"])

    def test_missing_duplicate_extra_seed_swap_and_partial_roster_reject(self) -> None:
        cases: List[Tuple[str, List[Dict[str, Any]], Sequence[str]]] = []
        cases.append(("missing", copy.deepcopy(self.valid_rows[1:]),
                      ("wirehair2_head", "wirehair1", "candidate")))
        cases.append(("duplicate", copy.deepcopy(self.valid_rows + [self.valid_rows[0]]),
                      ("wirehair2_head", "wirehair1", "candidate")))
        extra = copy.deepcopy(self.valid_rows)
        extra[0]["K"] = 7
        cases.append(("extra K", extra,
                      ("wirehair2_head", "wirehair1", "candidate")))
        swapped = copy.deepcopy(self.valid_rows)
        swapped[0]["loss_seed"] = \
            self.contract["seeds"]["training_loss_roots"][1]
        cases.append(("seed swap", swapped,
                      ("wirehair2_head", "wirehair1", "candidate")))
        cases.append(("partial roster", copy.deepcopy(self.valid_rows),
                      ("wirehair2_head", "wirehair1", "candidate", "third")))
        for name, rows, arms in cases:
            with self.subTest(name=name):
                self.assert_rejects(rows, arms)

    def test_wrong_band_trial_loss_and_cell_hash_reject(self) -> None:
        mutations = (
            ("band", "101-1000"),
            ("trial", 2),
            ("loss_ppm", 350000),
            ("cell_sha256", "0" * 64),
        )
        for field, value in mutations:
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_trace_arm_binary_and_descriptor_drift_reject(self) -> None:
        for field in ("trace_sha256", "binary_sha256", "arm_descriptor_sha256",
                      "realized_construction_sha256", "repair_map_sha256"):
            rows = copy.deepcopy(self.valid_rows)
            # Row 3 is the Wirehair2 control's next cell.  This makes each field drift
            # either between common-cell arms or within one frozen arm.
            rows[3][field] = digest("mutated:" + field)
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_duplicate_frozen_trace_contents_are_allowed(self) -> None:
        duplicate_traces = [digest("same legitimate trace")] * 360
        rows = ledger_rows(self.contract, trace_hashes=duplicate_traces)
        self.validate(rows, trace_hashes=duplicate_traces)

    def test_forged_trace_receipt_rejects_against_frozen_manifest(self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        rows[0]["trace_sha256"] = digest("forged trace")
        self.assert_rejects(rows)

    def test_trace_manifest_is_hashed_and_loaded_from_one_open(self) -> None:
        first_rows = trace_rows(self.contract, self.trace_hashes)
        second_rows = copy.deepcopy(first_rows)
        second_rows[0]["trace_sha256"] = digest("unfrozen replacement")
        first_data = "".join(
            subject.canonical_json(row) + "\n" for row in first_rows
        ).encode("utf-8")
        second_data = "".join(
            subject.canonical_json(row) + "\n" for row in second_rows
        ).encode("utf-8")

        class ReplacedPath:
            def __init__(self) -> None:
                self.opens = 0

            def open(self, _mode: str) -> io.BytesIO:
                self.opens += 1
                return io.BytesIO(first_data if self.opens == 1 else second_data)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stable = root / "stable-traces.jsonl"
            stable.write_bytes(first_data)
            expected = subject.trace_manifest_sha256(
                self.contract, "recovery", "development", stable)
        replaced = ReplacedPath()
        loaded = subject.load_trace_manifest(
            self.contract, "development", replaced, expected)  # type: ignore[arg-type]
        self.assertEqual(replaced.opens, 1)
        self.assertEqual(loaded[0], first_rows[0]["trace_sha256"])

    def test_raw_attempt_or_repair_map_rejects(self) -> None:
        for field, value in (
                ("construction_attempt", 1),
                ("construction_attempt", True),
                ("repair_map_sha256", digest("unexpected repair map"))):
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field, value=value):
                self.assert_rejects(rows)

    def test_repaired_validation_replays_exact_authenticated_map(self) -> None:
        phase = "cross_width_validation"
        trace_hashes = [
            digest("cross-trace:" + subject.sha256_json(cell))
            for cell in subject.iter_recovery_cells(self.contract, phase)
        ]
        training_trace_sha = digest("frozen final_repaired trace manifest")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            map_paths: Dict[str, Path] = {}
            map_hashes: Dict[str, str] = {}
            attempts: Dict[str, Mapping[int, int]] = {}
            for arm, offset in (("wirehair2_head", 7), ("candidate", 8)):
                map_path, map_hash, arm_attempts = write_repair_map(
                    root, self.contract, arm, training_trace_sha, {2: offset})
                map_paths[arm] = map_path
                map_hashes[arm] = map_hash
                attempts[arm] = arm_attempts
            trace_path = write_trace_manifest(
                root, self.contract, trace_hashes, phase)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"), phase,
                map_hashes, training_trace_sha)
            rows = ledger_rows(
                self.contract, trace_hashes=trace_hashes, phase=phase,
                repair_attempts=attempts, repair_map_hashes=map_hashes)
            ledger_path = write_ledger(root, rows)
            summary = subject.validate_ledger(
                self.contract, phase, ledger_path, freeze_path, trace_path,
                map_paths)
            self.assertTrue(summary["arms"]["candidate"][
                "phase_recovery_gate_pass"])

            mutated = copy.deepcopy(rows)
            candidate_row = next(
                row for row in mutated
                if row["arm"] == "candidate" and row["K"] == 2)
            candidate_row["construction_attempt"] = 9
            candidate_row["realized_construction_sha256"] = \
                subject.generic_realized_construction_sha256(
                    "wirehair2_experiment",
                    candidate_row["arm_descriptor_sha256"], 2,
                    candidate_row["block_bytes"], 9)
            ledger_path = write_ledger(root, mutated)
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, phase, ledger_path, freeze_path, trace_path,
                    map_paths)

    def test_missing_declared_control_rejects_before_scoring(self) -> None:
        rows = ledger_rows(
            self.contract, arms=("wirehair2_head", "candidate"))
        self.assert_rejects(rows, ("wirehair2_head", "candidate"))

    def test_fatal_unknown_and_contradictory_outcomes_reject(self) -> None:
        mutations = (
            ("fatal", None),
            ("success", None),
            ("unsupported", 4),
            ("success", 5),
        )
        for outcome, extra in mutations:
            rows = copy.deepcopy(self.valid_rows)
            rows[0]["outcome"] = outcome
            rows[0]["decoded_extra"] = extra
            with self.subTest(outcome=outcome, extra=extra):
                self.assert_rejects(rows)

    def test_noncanonical_types_reject(self) -> None:
        for field, value in (("trial", True), ("K", 2.0), ("arm", ["control"])):
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_duplicate_json_key_and_partial_line_reject(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            duplicate = root / "duplicate.jsonl"
            duplicate.write_bytes(b'{"arm":"control","arm":"candidate"}\n')
            trace_path = write_trace_manifest(root, self.contract, self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"))
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, "development", duplicate,
                    freeze_path, trace_path)
            partial = write_ledger(root, self.valid_rows, terminal_newline=False)
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, "development", partial,
                    freeze_path, trace_path)

    def test_canonical_crlf_evidence_is_portable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_trace_manifest(
                root, self.contract, self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"))
            ledger_path = write_ledger(root, self.valid_rows)
            for path in (trace_path, freeze_path, ledger_path):
                path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
            subject.validate_ledger(
                self.contract, "development", ledger_path,
                freeze_path, trace_path)


class TimingReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.arms = ("wirehair2_head", "wirehair1", "candidate")
        cls.hashes = timing_trace_hashes(cls.contract)
        cls.valid_rows = timing_receipt_rows(cls.contract, cls.hashes)

    def validate(self, rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_timing_trace_manifest(
                root, self.contract, self.hashes)
            freeze_path = write_timing_freeze_manifest(
                root, self.contract, trace_path, self.arms)
            receipt_path = write_timing_receipt(root, rows)
            return subject.validate_timing_receipt(
                self.contract, "development", receipt_path,
                freeze_path, trace_path)

    def assert_rejects(self, rows: Sequence[Mapping[str, Any]]) -> None:
        with self.assertRaises(subject.ContractError):
            self.validate(rows)

    def test_exact_panels_counterbalancing_and_student_t_gate(self) -> None:
        summary = self.validate(self.valid_rows)
        self.assertEqual(summary["rows"], 192 * 11)
        self.assertTrue(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])
        decisions = summary["candidates"]["candidate"]["panels"]
        for panel in decisions.values():
            for decision in panel["by_band_width"].values():
                self.assertEqual(decision["status"], "faster")
                self.assertLess(decision["upper_95"],
                                -decision["effective_floor"])

    def test_resolved_slowdown_fails_development_gate(self) -> None:
        rows = timing_receipt_rows(
            self.contract, self.hashes, candidate_scale=1.2)
        summary = self.validate(rows)
        self.assertFalse(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])
        statuses = {
            decision["status"]
            for panel in summary["candidates"]["candidate"]["panels"].values()
            for decision in panel["by_band_width"].values()
        }
        self.assertEqual(statuses, {"resolved_slower"})

    def test_aa_noise_cannot_expand_the_allowed_regression(self) -> None:
        rows = timing_receipt_rows(
            self.contract, self.hashes, candidate_scale=1.2)
        for row in rows:
            if row["panel_kind"] != "AA":
                continue
            row["elapsed_ns"] = \
                [200000, 100000, 100000, 200000] \
                if row["order"] == "ABBA" else \
                [100000, 200000, 200000, 100000]
        summary = self.validate(rows)
        candidate = summary["candidates"]["candidate"]
        self.assertFalse(candidate["phase_speed_gate_pass"])
        reasons = {
            decision.get("reason")
            for panel in candidate["panels"].values()
            for decision in panel["by_band_width"].values()
        }
        self.assertEqual(reasons, {"aa_repeatability_threshold"})

    def test_final_speed_gate_requires_all_four_declared_wins(self) -> None:
        faster = {"by_band_width": {"band|64": {"status": "faster"}},
                  "aggregate": {"status": "faster"}}
        decisions = {
            "decoder_solve|wirehair2_head": copy.deepcopy(faster),
            "receive_to_success|wirehair1": copy.deepcopy(faster),
            "encoder_init_plus_first_K_symbols|wirehair2_head":
                copy.deepcopy(faster),
            "encoder_init_plus_first_K_symbols|wirehair1":
                copy.deepcopy(faster),
        }
        self.assertTrue(subject._timing_speed_gate("final", decisions))
        mutations = (
            ("decoder_solve|wirehair2_head", "band", "noninferior"),
            ("receive_to_success|wirehair1", "band", "noninferior"),
            ("encoder_init_plus_first_K_symbols|wirehair2_head",
             "band", "unresolved"),
            ("encoder_init_plus_first_K_symbols|wirehair2_head",
             "aggregate", "noninferior"),
            ("encoder_init_plus_first_K_symbols|wirehair1",
             "aggregate", "noninferior"),
        )
        for panel, scope, status in mutations:
            changed = copy.deepcopy(decisions)
            if scope == "band":
                changed[panel]["by_band_width"]["band|64"]["status"] = status
            else:
                changed[panel]["aggregate"]["status"] = status
            with self.subTest(panel=panel, scope=scope):
                self.assertFalse(subject._timing_speed_gate("final", changed))

    def test_missing_order_trace_timer_and_identity_mutants_reject(self) -> None:
        cases: List[Tuple[str, List[Dict[str, Any]]]] = []
        cases.append(("missing", copy.deepcopy(self.valid_rows[1:])))
        wrong_order = copy.deepcopy(self.valid_rows)
        wrong_order[0]["order"] = \
            "BAAB" if wrong_order[0]["order"] == "ABBA" else "ABBA"
        cases.append(("order", wrong_order))
        wrong_trace = copy.deepcopy(self.valid_rows)
        wrong_trace[0]["trace_sha256"] = digest("wrong timing trace")
        cases.append(("trace", wrong_trace))
        zero_timer = copy.deepcopy(self.valid_rows)
        zero_timer[0]["elapsed_ns"][0] = 0
        cases.append(("timer", zero_timer))
        descriptor_drift = copy.deepcopy(self.valid_rows)
        descriptor_drift[0]["left_arm_descriptor_sha256"] = digest("drift")
        cases.append(("identity", descriptor_drift))
        for name, rows in cases:
            with self.subTest(name=name):
                self.assert_rejects(rows)

    def test_explicit_failure_is_not_dropped_and_makes_panel_nonselectable(
            self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        target = rows[0]
        target_cell = {key: target[key] for key in self.contract["timing"][
            "cell_key"]}
        target_scope = target["scope"]
        target_arm = target["left_arm"]
        for row in rows:
            if (all(row[key] == value for key, value in target_cell.items()) and
                    row["scope"] == target_scope):
                changed = False
                for side in ("left", "right"):
                    if row[side + "_arm"] == target_arm:
                        row[side + "_outcome"] = "need_more_at_cap"
                        row[side + "_decoded_extra"] = None
                        changed = True
                if changed:
                    row["elapsed_ns"] = [None, None, None, None]
        summary = self.validate(rows)
        self.assertFalse(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])


class SelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()

    def test_recovery_equivalent_set_precedes_speed_ranking(self) -> None:
        artifact = digest("one development architecture roster")
        candidates = {
            "fast_but_worse": (7, -0.40),
            "reliable": (5, -0.05),
            "selected": (6, -0.20),
        }
        controls = ("wirehair2_head", "wirehair1")
        arm_failures = {
            **{arm: 0 for arm in controls},
            **{arm: value[0] for arm, value in candidates.items()},
        }
        arm_artifacts = {
            arm: {
                "codec": arm_codec(arm),
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": digest("descriptor:" + arm),
            }
            for arm in arm_failures
        }
        recovery: Dict[str, Any] = {
            "schema": subject.SCHEMA + ".ledger-summary.v1",
            "phase": "development",
            "contract_sha256": subject.contract_sha256(self.contract),
            "domain_sha256": self.contract["recovery"]["domains"]
                ["development"]["domain_sha256"],
            "freeze_manifest_sha256": digest("recovery freeze"),
            "architecture_artifact_sha256": artifact,
            "arm_artifacts": copy.deepcopy(arm_artifacts),
            "excluded_cells": 0,
            "mandatory_controls_supported": True,
            "comparisons": {
                arm: {
                    "architecture_eligible": True,
                    "controls": {control: {} for control in controls},
                }
                for arm in candidates
            },
            "arms": {
                arm: {
                    "cells": 360,
                    "failure_by_overhead": {
                        "0": failure0, "1": max(0, failure0 - 2),
                        "2": max(0, failure0 - 4), "4": 0,
                    },
                }
                for arm, failure0 in arm_failures.items()
            },
        }
        timing: Dict[str, Any] = {
            "schema": subject.SCHEMA + ".timing-summary.v1",
            "phase": "development",
            "contract_sha256": subject.contract_sha256(self.contract),
            "domain_sha256": self.contract["timing"]["domains"]
                ["development"]["domain_sha256"],
            "freeze_manifest_sha256": digest("timing freeze"),
            "architecture_artifact_sha256": artifact,
            "arm_artifacts": copy.deepcopy(arm_artifacts),
            "excluded_cells": 0,
            "rows": 192 * (4 + 7 * len(candidates)),
            "candidates": {
                arm: {
                    "phase_speed_gate_pass": True,
                    "panels": {
                        "decoder_solve|wirehair2_head": {
                            "aggregate": {
                                "selectable": True,
                                "mean_log_ratio": mean,
                            },
                        },
                    },
                }
                for arm, (_failure0, mean) in candidates.items()
            },
        }
        result = subject.select_development_architecture(
            self.contract, recovery, timing)
        subject.validate_selection_receipt(self.contract, result)
        self.assertEqual(result["recovery_equivalence_allowance"], 1)
        self.assertEqual(
            result["recovery_equivalent_candidates"],
            ["reliable", "selected"])
        self.assertEqual(result["selected_arm"], "selected")
        self.assertEqual(result["candidate_roster"], sorted(candidates))
        tampered_receipt = copy.deepcopy(result)
        tampered_receipt["recovery_freeze_manifest_sha256"] = \
            digest("substituted after selection")
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract, tampered_receipt)
        forged_minimum = copy.deepcopy(result)
        forged_minimum["minimum_overhead0_failures"] = 0
        unsigned = {key: value for key, value in forged_minimum.items()
                    if key != "selection_sha256"}
        forged_minimum["selection_sha256"] = subject.sha256_json(unsigned)
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract, forged_minimum)
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract,
                architecture_selection_receipt(self.contract, "wirehair1"))

        mutations = []
        changed = copy.deepcopy(recovery)
        changed["domain_sha256"] = digest("other recovery domain")
        mutations.append(("recovery domain", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["excluded_cells"] = 1
        mutations.append(("recovery exclusion", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["arms"]["selected"]["cells"] = 1
        mutations.append(("recovery cardinality", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["arms"]["selected"]["failure_by_overhead"]["1"] = 7
        mutations.append(("nonnested tail", changed, timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["domain_sha256"] = digest("other timing domain")
        mutations.append(("timing domain", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["rows"] = 1
        mutations.append(("timing cardinality", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["freeze_manifest_sha256"] = "bad"
        mutations.append(("timing freeze", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["architecture_artifact_sha256"] = digest("other binary")
        mutations.append(("architecture", recovery, changed_timing))
        for name, recovery_value, timing_value in mutations:
            with self.subTest(mutation=name), \
                    self.assertRaises(subject.ContractError):
                subject.select_development_architecture(
                    self.contract, recovery_value, timing_value)


if __name__ == "__main__":
    unittest.main()
