#!/usr/bin/env python3
"""Regression tests for the additive repair-v1 v3 classification erratum."""

import copy
import hashlib
import importlib
import os
from pathlib import Path
import subprocess
import sys
import unittest


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent
TOOLS = REPOSITORY / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(HERE))

erratum = importlib.import_module("wh2_rv4a_v3_erratum")
classification = erratum.classification
parser_module = erratum.campaign._require_v3_parser()


ARM_ID = 0x19cccf775ce0bf09
CONSTRUCTION_ROOT = 0x12345678
SCHEDULES = (
    "iid",
    "burst",
    "permutation",
    "systematic-first",
    "repair-only",
    "adversarial",
)


def _descriptor(root, selected_attempt, arm_id=ARM_ID):
    payload = (
        b"W2R3" +
        arm_id.to_bytes(8, "little") +
        root.to_bytes(4, "little") +
        bytes((selected_attempt,))
    )
    return payload.hex(), hashlib.sha256(payload).hexdigest()


def _attempt_row(attempt):
    row = {
        field: -1
        for field in parser_module.REPAIRTIMING_ATTEMPT_FIELDS
    }
    row.update({
        "attempt": attempt,
        "probe_executed": 0,
        "probe_result": -1,
        "probe_stats_available": 0,
        "real_executed": 0,
        "real_result": -1,
        "real_stats_available": 0,
    })
    return row


def _execute_attempt_call(row, kind, result):
    row[f"{kind}_executed"] = 1
    row[f"{kind}_result"] = result
    row[f"{kind}_stats_available"] = 1
    for field in parser_module.REPAIRTIMING_STATS_FIELDS:
        row[f"{kind}_{field}"] = 1


def _selector_attempt(row):
    return {
        "attempt": row["attempt"],
        "probe_executed": row["probe_executed"],
        "probe_result": row["probe_result"],
        "probe_stats_available": row["probe_stats_available"],
        "probe_stats": {
            field: row[f"probe_{field}"]
            for field in parser_module.REPAIRTIMING_STRUCTURAL_STATS_FIELDS
        },
    }


def _timing_row(panel, panel_index, slot):
    row = {
        field: 0
        for field in parser_module.REPAIRTIMING_TIMING_FIELDS
    }
    row.update({
        "timing_panel": panel,
        "timing_panel_index": panel_index,
        "timing_slot": slot,
        "timing_pair": slot // 2,
        "timing_label": "ABBABAAB"[slot],
        "timing_role": "repair_forced_encoder",
        "timing_scope": "encoder_forced",
        "timing_censor_reason": "none",
        "timing_construct_result": 0,
        "timing_result": 0,
        "timing_recover_result": 0,
        "timing_eligible": 1,
        "timing_executed": 1,
        "timing_recovery_ok": 1,
        "timing_stable": 1,
        "timing_elapsed_ns": 100,
        "timing_inner_reps": 1,
        "timing_saturated": 0,
        "timing_cpu_before": 0,
        "timing_cpu_after": 0,
        "timing_cpu_migrated": 0,
        "timing_minflt": 0,
        "timing_majflt": 0,
        "timing_fault_contaminated": 0,
        "timing_stats_available": 1,
        "timing_fixed_prefix_symbols": 2,
        "timing_intermediate_bytes": 2,
        "timing_intermediate_sha256": "a" * 64,
    })
    return row


def _role(*, encode=0, construct=0, feed=0, recover=0,
          recovery_ok=1, outcome="success"):
    return {
        "encode_result": encode,
        "decode_construct_result": construct,
        "feed_result": feed,
        "recover_result": recover,
        "recovery_ok": recovery_ok,
        "encoded_symbols": 2 if encode == 0 else 0,
        "received_symbols": 2 if recovery_ok else 0,
        "overhead": 0 if recovery_ok else -1,
        "payload_sha256": "b" * 64 if encode == 0 else "not_applicable",
        "recovered_sha256":
            "c" * 64 if recovery_ok else "not_applicable",
        "encode_class": (
            "success" if encode == 0
            else "need_more" if encode == 1
            else "weak" if encode in (3, 4)
            else "error"
        ),
        "feed_class": "success" if feed == 0 else "not_applicable",
        "recover_class": "success" if recover == 0 else "not_applicable",
        "outcome_class": outcome,
    }


def _controls():
    controls = {
        field: 0
        for field in parser_module.REPAIRTIMING_CONTROL_FIELDS
    }
    controls.update({
        "forced_a_result": 0,
        "forced_b_result": 0,
        "forced_a_payload_sha256": "b" * 64,
        "forced_b_payload_sha256": "b" * 64,
        "forced_equal": 1,
        "direct_fixed_prefix_symbols": 2,
        "repair_intermediate_bytes": 2,
        "repair_intermediate_sha256": "d" * 64,
        "repair_direct_executed": 1,
        "repair_direct_result": 0,
        "repair_direct_intermediate_bytes": 2,
        "repair_direct_intermediate_sha256": "d" * 64,
        "repair_direct_witness_equal": 1,
        "dispatch_intermediate_bytes": 2,
        "dispatch_intermediate_sha256": "e" * 64,
        "dispatch_direct_executed": 1,
        "dispatch_direct_result": 0,
        "dispatch_direct_intermediate_bytes": 2,
        "dispatch_direct_intermediate_sha256": "e" * 64,
        "dispatch_direct_witness_equal": 1,
    })
    return controls


def _make_fixture(kind="corroborated_error", *, root=CONSTRUCTION_ROOT):
    attempts = [_attempt_row(index) for index in range(8)]
    if kind == "success":
        _execute_attempt_call(attempts[0], "real", 0)
        selected = 0
        selector_result = 0
        committed = 1
        raw = _role()
        repaired = _role()
    elif kind in ("explicit_weak_3", "explicit_weak_4"):
        code = int(kind[-1])
        _execute_attempt_call(attempts[0], "real", code)
        selected = -1
        selector_result = code
        committed = 0
        raw = _role(
            encode=code, construct=-1, feed=-1, recover=-1,
            recovery_ok=0, outcome="weak")
        repaired = copy.deepcopy(raw)
    else:
        initial = 1 if kind == "need_more" else 8
        _execute_attempt_call(attempts[0], "probe", 1)
        _execute_attempt_call(attempts[0], "real", initial)
        _execute_attempt_call(attempts[1], "probe", 0)
        _execute_attempt_call(attempts[1], "real", 0)
        selected = 1
        selector_result = 0
        committed = 1
        raw = _role(
            encode=initial, construct=-1, feed=-1, recover=-1,
            recovery_ok=0,
            outcome="need_more" if initial == 1 else "error")
        repaired = _role()

    if committed:
        descriptor_hex, descriptor_sha256 = _descriptor(root, selected)
    else:
        descriptor_hex = descriptor_sha256 = "not_applicable"
    selector = {
        "schema": "wirehair.wh2.repair-v1.selector-structural.v1",
        "arm_id": ARM_ID,
        "K": 2,
        "construction_root": root,
        "repair_policy_sha256":
            erratum.campaign.REPAIR_POLICY_SHA256,
        "selector_result": selector_result,
        "selected_attempt": selected,
        "attempts_executed": 1 if selected in (-1, 0) else selected + 1,
        "calls_executed": sum(
            row[f"{call}_executed"]
            for row in attempts
            for call in ("probe", "real")
        ),
        "real_configuration_calls": sum(
            row["real_executed"] for row in attempts),
        "structural_probe_calls": sum(
            row["probe_executed"] for row in attempts),
        "cap_exhausted": 0,
        "fatal_attempt_zero_mismatch": 0,
        "oom": 0,
        "committed": committed,
        "descriptor_sha256": descriptor_sha256,
        "attempts": [_selector_attempt(row) for row in attempts],
        "selector_structural_sha256": "f" * 64,
    }
    timing = [
        _timing_row(panel, panel_index, slot)
        for panel_index, panel in enumerate(
            parser_module.REPAIRTIMING_PANELS)
        for slot in range(8)
    ]
    real = {
        "schema": "wirehair.wh2.repairtiming.real-witness.v3",
        "trace_sha256": "1" * 64,
        "message_sha256": "2" * 64,
        "descriptor_hex": descriptor_hex,
        "roles": {
            "raw": raw,
            "repaired": repaired,
            "dispatch": _role(),
            "wh1": _role(),
        },
        "controls": _controls(),
        "attempt_columns":
            list(parser_module.REPAIRTIMING_ATTEMPT_FIELDS),
        "attempt_rows": [
            [row[field] for field in parser_module.REPAIRTIMING_ATTEMPT_FIELDS]
            for row in attempts
        ],
        "timing_columns":
            list(parser_module.REPAIRTIMING_TIMING_FIELDS),
        "timing_rows": [
            [row[field] for field in parser_module.REPAIRTIMING_TIMING_FIELDS]
            for row in timing
        ],
    }
    return selector, real


def _set_attempt(real, attempt, field, value):
    column = real["attempt_columns"].index(field)
    real["attempt_rows"][attempt][column] = value


def _set_timing(real, row, field, value):
    column = real["timing_columns"].index(field)
    real["timing_rows"][row][column] = value


def _classify(selector, real):
    return classification.classify_repair_v1_cell(
        selector, real, parser_module=parser_module)


class RepairV1CellClassificationTests(unittest.TestCase):
    def assert_fails_closed(self, selector, real):
        try:
            result = _classify(selector, real)
        except classification.ClassificationError:
            return
        self.assertEqual(result["raw_attempt0_kind"], "none")
        self.assertEqual(result["raw_attempt0_structural_weak"], 0)
        self.assertEqual(result["candidate_runtime_error"], 1)

    def test_golden_classes(self):
        expected = {
            "success": ("none", 0, 0, 0),
            "explicit_weak_3": ("explicit_weak", 1, 0, 1),
            "explicit_weak_4": ("explicit_weak", 1, 0, 1),
            "need_more": ("need_more", 1, 0, 0),
            "corroborated_error": ("corroborated_error", 1, 0, 0),
        }
        for fixture, wanted in expected.items():
            with self.subTest(fixture=fixture):
                result = _classify(*_make_fixture(fixture))
                self.assertEqual((
                    result["raw_attempt0_kind"],
                    result["raw_attempt0_structural_weak"],
                    result["candidate_runtime_error"],
                    result["repaired_final_weak"],
                ), wanted)
                self.assertEqual(result["runtime_error_observations"], {})

    def test_every_error_corroboration_field_fails_closed(self):
        mutations = {
            "raw-result": lambda s, r:
                r["roles"]["raw"].__setitem__("encode_result", 0),
            "attempt0-real-executed": lambda s, r:
                _set_attempt(r, 0, "real_executed", 0),
            "attempt0-real-result": lambda s, r:
                _set_attempt(r, 0, "real_result", 0),
            "selector-attempt0-probe-executed": lambda s, r:
                s["attempts"][0].__setitem__("probe_executed", 0),
            "selector-attempt0-probe-result": lambda s, r:
                s["attempts"][0].__setitem__("probe_result", 0),
            "real-attempt0-probe-executed": lambda s, r:
                _set_attempt(r, 0, "probe_executed", 0),
            "real-attempt0-probe-result": lambda s, r:
                _set_attempt(r, 0, "probe_result", 0),
            "selected-attempt": lambda s, r:
                s.__setitem__("selected_attempt", 0),
            "selector-selected-probe-executed": lambda s, r:
                s["attempts"][1].__setitem__("probe_executed", 0),
            "selector-selected-probe-result": lambda s, r:
                s["attempts"][1].__setitem__("probe_result", 1),
            "real-selected-probe-executed": lambda s, r:
                _set_attempt(r, 1, "probe_executed", 0),
            "real-selected-probe-result": lambda s, r:
                _set_attempt(r, 1, "probe_result", 1),
            "selected-real-executed": lambda s, r:
                _set_attempt(r, 1, "real_executed", 0),
            "selected-real-result": lambda s, r:
                _set_attempt(r, 1, "real_result", 1),
            "selector-result": lambda s, r:
                s.__setitem__("selector_result", 1),
            "committed": lambda s, r: s.__setitem__("committed", 0),
            "cap-exhausted": lambda s, r:
                s.__setitem__("cap_exhausted", 1),
            "fatal-mismatch": lambda s, r:
                s.__setitem__("fatal_attempt_zero_mismatch", 1),
            "oom": lambda s, r: s.__setitem__("oom", 1),
            "descriptor-sha": lambda s, r:
                s.__setitem__("descriptor_sha256", "0" * 64),
            "descriptor-magic": lambda s, r:
                r.__setitem__("descriptor_hex", "00" * 17),
            "descriptor-attempt": lambda s, r:
                r.__setitem__(
                    "descriptor_hex", _descriptor(
                        CONSTRUCTION_ROOT, 2)[0]),
            "descriptor-arm": lambda s, r:
                s.__setitem__("arm_id", ARM_ID + 1),
            "descriptor-root": lambda s, r:
                s.__setitem__("construction_root", CONSTRUCTION_ROOT + 1),
            "repaired-encode": lambda s, r:
                r["roles"]["repaired"].__setitem__("encode_result", 1),
            "repaired-construct": lambda s, r:
                r["roles"]["repaired"].__setitem__(
                    "decode_construct_result", 1),
            "repaired-feed": lambda s, r:
                r["roles"]["repaired"].__setitem__("feed_result", 1),
            "repaired-recover": lambda s, r:
                r["roles"]["repaired"].__setitem__("recover_result", 1),
            "repaired-recovery": lambda s, r:
                r["roles"]["repaired"].__setitem__("recovery_ok", 0),
            "repaired-outcome": lambda s, r:
                r["roles"]["repaired"].__setitem__(
                    "outcome_class", "need_more"),
        }
        for name, mutate in mutations.items():
            with self.subTest(field=name):
                selector, real = _make_fixture()
                mutate(selector, real)
                self.assert_fails_closed(selector, real)

    def test_every_nonexempt_role_and_control_error_remains_fatal(self):
        locations = [(
            "selector.selector_result",
            lambda s, r: s.__setitem__("selector_result", 8),
        )]
        for role in ("raw", "repaired", "dispatch", "wh1"):
            for field in (
                    "encode_result", "decode_construct_result",
                    "feed_result", "recover_result"):
                if role == "raw" and field == "encode_result":
                    continue
                locations.append((
                    f"roles.{role}.{field}",
                    lambda s, r, role=role, field=field:
                        r["roles"][role].__setitem__(field, 8),
                ))
        for field in (
                "forced_a_result", "forced_b_result",
                "repair_direct_result", "dispatch_direct_result"):
            locations.append((
                f"controls.{field}",
                lambda s, r, field=field:
                    r["controls"].__setitem__(field, 8),
            ))
        for location, mutate in locations:
            with self.subTest(location=location):
                selector, real = _make_fixture()
                mutate(selector, real)
                result = _classify(selector, real)
                self.assertEqual(result["candidate_runtime_error"], 1)
                self.assertIn(
                    f"8:{location}", result["runtime_error_observations"])

    def test_every_nonexempt_attempt_error_remains_fatal(self):
        cases = (
            ("attempts.0.probe_result", 0, "probe_result"),
            ("attempts.1.probe_result", 1, "probe_result"),
            ("attempts.1.real_result", 1, "real_result"),
            ("attempts.2.probe_result", 2, "probe_result"),
            ("attempts.2.real_result", 2, "real_result"),
        )
        for location, attempt, result_field in cases:
            with self.subTest(location=location):
                selector, real = _make_fixture()
                self.assertTrue(result_field.endswith("_result"))
                kind = result_field[:-len("_result")]
                _set_attempt(real, attempt, f"{kind}_executed", 1)
                _set_attempt(real, attempt, result_field, 8)
                if attempt == 0 and kind == "probe":
                    selector["attempts"][0]["probe_result"] = 8
                result = _classify(selector, real)
                self.assertEqual(result["candidate_runtime_error"], 1)
                self.assertIn(
                    f"8:{location}", result["runtime_error_observations"])

    def test_every_timing_result_position_remains_fatal(self):
        rows = len(parser_module.REPAIRTIMING_PANELS) * 8
        for row in range(rows):
            for field in (
                    "timing_construct_result",
                    "timing_result",
                    "timing_recover_result"):
                with self.subTest(row=row, field=field):
                    selector, real = _make_fixture()
                    _set_timing(real, row, field, 8)
                    result = _classify(selector, real)
                    self.assertEqual(
                        result["candidate_runtime_error"], 1)
                    self.assertIn(
                        f"8:timing.{field}",
                        result["runtime_error_observations"],
                    )

    def test_other_runtime_codes_and_legacy_nonruntime_codes_are_preserved(self):
        for code in (2, 5, 6, 10):
            with self.subTest(runtime_code=code):
                selector, real = _make_fixture()
                real["controls"]["dispatch_direct_result"] = code
                result = _classify(selector, real)
                self.assertEqual(result["candidate_runtime_error"], 1)
                self.assertIn(
                    f"{code}:controls.dispatch_direct_result",
                    result["runtime_error_observations"],
                )
        for code in (-1, 0, 1, 3, 4, 7, 9):
            with self.subTest(nonruntime_code=code):
                selector, real = _make_fixture()
                real["controls"]["dispatch_direct_result"] = code
                result = _classify(selector, real)
                self.assertEqual(result["candidate_runtime_error"], 0)


class SelectorScheduleAggregationTests(unittest.TestCase):
    def records(self, roots=(CONSTRUCTION_ROOT,), fixture="need_more"):
        records = []
        for replicate, root in enumerate(roots):
            classification_result = _classify(
                *_make_fixture(fixture, root=root))
            for schedule in SCHEDULES:
                records.append({
                    "selector_key": [ARM_ID, 2, root],
                    "schedule": schedule,
                    "replicate": replicate,
                    "classification": copy.deepcopy(
                        classification_result),
                })
        return records

    def aggregate(self, records):
        return classification.aggregate_selector_classifications(
            records,
            schedules=SCHEDULES,
            first_schedule=SCHEDULES[0],
        )

    def test_six_schedules_count_one_selector(self):
        result = self.aggregate(self.records())
        self.assertEqual(result, {
            "unique_selectors": 1,
            "raw_attempt0_structural_weak": {
                "total": 1,
                "by_kind": {
                    "explicit_weak": 0,
                    "need_more": 1,
                    "corroborated_error": 0,
                },
            },
            "repaired_final_weak": {
                "unique_selectors": 0,
            },
        })

    def test_two_roots_count_two_not_twelve(self):
        result = self.aggregate(self.records(
            roots=(CONSTRUCTION_ROOT, CONSTRUCTION_ROOT + 1),
            fixture="corroborated_error",
        ))
        self.assertEqual(result["unique_selectors"], 2)
        self.assertEqual(
            result["raw_attempt0_structural_weak"]["total"], 2)
        self.assertEqual(
            result["raw_attempt0_structural_weak"]["by_kind"]
            ["corroborated_error"],
            2,
        )

    def test_missing_and_duplicate_schedule_are_rejected(self):
        missing = self.records()[:-1]
        with self.assertRaisesRegex(
                classification.ClassificationError, "lacks"):
            self.aggregate(missing)

        duplicate = self.records()
        duplicate.append(copy.deepcopy(duplicate[0]))
        with self.assertRaisesRegex(
                classification.ClassificationError, "repeats"):
            self.aggregate(duplicate)

    def test_schedule_dependent_classification_is_rejected(self):
        mutations = (
            {
                "raw_attempt0_kind": "none",
                "raw_attempt0_structural_weak": 0,
            },
            {
                "repaired_final_weak": 1,
            },
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                records = self.records()
                records[-1]["classification"].update(mutation)
                with self.assertRaisesRegex(
                        classification.ClassificationError,
                        "changed by schedule"):
                    self.aggregate(records)

    def test_schedule_dependent_replicate_is_rejected(self):
        records = self.records()
        records[-1]["replicate"] += 1
        with self.assertRaisesRegex(
                classification.ClassificationError, "replicate"):
            self.aggregate(records)

    def test_malformed_schedule_contract_and_record_are_rejected(self):
        records = self.records()
        with self.assertRaisesRegex(
                classification.ClassificationError, "contract"):
            classification.aggregate_selector_classifications(
                records,
                schedules=SCHEDULES,
                first_schedule=SCHEDULES[1],
            )
        malformed = self.records()
        malformed[0]["selector_key"] = [ARM_ID, 2, []]
        with self.assertRaisesRegex(
                classification.ClassificationError, "unhashable"):
            self.aggregate(malformed)


class SourceAndArtifactBindingTests(unittest.TestCase):
    def test_policy_scopes_authorization_to_this_no_survivor_artifact(self):
        self.assertEqual(
            erratum.classification_policy()["outcome_use"],
            "reporting-erratum-only-this-historical-v1-no-survivor-"
            "decision-and-this-erratum-grant-no-sealed-or-public-"
            "promotion-authorization",
        )

    def test_classifier_is_the_auditor_source_loaded_module(self):
        self.assertIs(classification, erratum.classification)
        source = erratum._source_bytes(
            TOOLS / "repair_v1_classification.py")
        self.assertEqual(
            classification.__rv4a_source_sha256__,
            hashlib.sha256(source).hexdigest(),
        )

    def test_git_source_receipt_binds_exact_committed_bytes(self):
        commit = subprocess.run(
            ["git", "-C", str(REPOSITORY), "rev-parse", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
            timeout=30,
        ).stdout.strip()
        relative = "bench/wh2_rv4a_campaign.py"
        receipt = erratum._git_source_receipt(commit, [relative])
        payload = subprocess.run(
            ["git", "-C", str(REPOSITORY), "show", f"{commit}:{relative}"],
            check=True,
            stdout=subprocess.PIPE,
            timeout=30,
        ).stdout
        self.assertEqual(receipt["commit"], commit)
        self.assertEqual(receipt["files"][relative], {
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        })

    def test_git_source_receipt_rejects_path_escape(self):
        with self.assertRaises(erratum.ErratumError):
            erratum._git_source_receipt("HEAD", ["../outside"])

    @unittest.skipUnless(
        os.environ.get("WH2_RV4A3_TRAINING_DIR"),
        "set WH2_RV4A3_TRAINING_DIR for authenticated artifact binding",
    )
    def test_authenticated_v3_artifact_hashes_are_immutable(self):
        directory = Path(os.environ["WH2_RV4A3_TRAINING_DIR"])
        source = erratum._read_completion_source(directory)
        expected = {
            "manifest.json":
                "c11d15c1c0c17e226ea45930cf11fb9d8cd9b43c05a4a3ff8d95c7b9f23a534a",
            "completion.json":
                "ab2a269488e64a5fba03716df9c284084cc3a71f6cdf69cb07b61d969456059c",
            "summary.json.gz":
                "21a23a5d940350dfd60e72276c82b929ce3738ce878d8adcd4fa510c3ed3b816",
            "cell-ledger.jsonl.gz":
                "df659fbde0965ea38f03bc44f20ccebd32368240da456df4184d0943b56adae1",
        }
        for relative, digest in expected.items():
            with self.subTest(path=relative):
                self.assertEqual(
                    erratum.campaign._stable_file_binding(
                        directory / relative,
                        byte_limit=4 * 1024 * 1024,
                    )["sha256"],
                    digest,
                )
        self.assertEqual(
            source["completion_sha256"], expected["completion.json"])
        self.assertEqual(
            source["manifest_sha256"], expected["manifest.json"])
        corrected = {
            "pure8_s0m1_d3": {
                "error": 902,
                "need_more": 384,
                "total": 1286,
            },
            "pure9_s0m1_d3": {
                "error": 737,
                "need_more": 390,
                "total": 1127,
            },
        }
        for arm, expected_counts in corrected.items():
            with self.subTest(arm=arm):
                legacy = source["summary"]["arms"][arm]
                outcomes = legacy["recovery"]["raw"]["outcome_classes"]
                counts = {
                    "error": outcomes["error"] // len(SCHEDULES),
                    "need_more": outcomes["need_more"] // len(SCHEDULES),
                }
                self.assertEqual(
                    outcomes["error"] % len(SCHEDULES), 0)
                self.assertEqual(
                    outcomes["need_more"] % len(SCHEDULES), 0)
                counts["total"] = counts["error"] + counts["need_more"]
                self.assertEqual(counts, expected_counts)
                self.assertEqual(
                    legacy["selector"]["weak_constructions"]
                    ["repaired_final"],
                    0,
                )
                self.assertEqual(legacy["candidate_final_weak"], 0)

    @unittest.skipUnless(
        os.environ.get("WH2_RV4A3_ERRATUM_RECEIPT"),
        "set WH2_RV4A3_ERRATUM_RECEIPT for corrected receipt binding",
    )
    def test_authenticated_erratum_has_exact_corrected_counts(self):
        path = Path(os.environ["WH2_RV4A3_ERRATUM_RECEIPT"])
        expected_receipt_sha256 = os.environ.get(
            "WH2_RV4A3_ERRATUM_RECEIPT_SHA256")
        expected_classifier_commit = os.environ.get(
            "WH2_RV4A3_ERRATUM_CLASSIFIER_COMMIT")
        self.assertRegex(
            expected_receipt_sha256 or "", r"^[0-9a-f]{64}$")
        self.assertRegex(
            expected_classifier_commit or "", r"^[0-9a-f]{40}$")
        receipt, receipt_sha256 = \
            erratum.campaign.read_canonical_json(path)
        self.assertEqual(
            receipt_sha256, expected_receipt_sha256,
        )
        self.assertEqual(receipt["schema"], erratum.ERRATUM_SCHEMA)
        expected_policy = erratum.classification_policy()
        self.assertEqual(receipt["policy"], expected_policy)
        self.assertEqual(
            receipt["policy_sha256"],
            erratum.campaign.canonical_sha256(expected_policy),
        )
        classifier_paths = [
            "bench/wh2_rv4a_campaign.py",
            "bench/wh2_rv4a_v3_erratum.py",
            "tools/repair_v1_classification.py",
        ]
        self.assertEqual(
            receipt["source_evidence"]["classification_source"],
            erratum._git_source_receipt(
                expected_classifier_commit, classifier_paths),
        )
        for relative in classifier_paths:
            with self.subTest(classifier_source=relative):
                binding = erratum.campaign._stable_file_binding(
                    REPOSITORY / relative,
                    byte_limit=erratum.SOURCE_BYTE_LIMIT,
                )
                self.assertEqual(
                    receipt["source_evidence"]["classification_source"]
                    ["files"][relative],
                    {
                        "bytes": binding["size"],
                        "sha256": binding["sha256"],
                    },
                )
        self.assertEqual(
            receipt["source_evidence"]["manifest_sha256"],
            "c11d15c1c0c17e226ea45930cf11fb9d8cd9b43c05a4a3ff8d95c7b9f23a534a",
        )
        self.assertEqual(
            receipt["source_evidence"]["completion_sha256"],
            "ab2a269488e64a5fba03716df9c284084cc3a71f6cdf69cb07b61d969456059c",
        )
        self.assertEqual(
            receipt["outcome_authorization"],
            "forbidden-v1-decision-retained-only-as-historical-evidence",
        )
        self.assertEqual(
            receipt["promotion_authorization"],
            "forbidden-reporting-erratum-is-not-a-promotion-decision",
        )
        expected = {
            "pure8_s0m1_d3": {
                "total": 1286,
                "by_kind": {
                    "explicit_weak": 0,
                    "need_more": 384,
                    "corroborated_error": 902,
                },
            },
            "pure9_s0m1_d3": {
                "total": 1127,
                "by_kind": {
                    "explicit_weak": 0,
                    "need_more": 390,
                    "corroborated_error": 737,
                },
            },
        }
        for arm, expected_weak in expected.items():
            with self.subTest(arm=arm):
                arm_receipt = receipt["arms"][arm]
                self.assertEqual(
                    arm_receipt["unique_selectors"], 76032)
                self.assertEqual(
                    arm_receipt["raw_attempt0_structural_weak"],
                    expected_weak,
                )
                self.assertEqual(
                    arm_receipt["candidate_runtime_error"], 0)
                self.assertEqual(
                    arm_receipt["repaired_final_weak"], {
                        "unique_selectors": 0,
                        "physical_cells": 0,
                    })

    @unittest.skipUnless(
        os.environ.get("WH2_RV4A3_VERIFICATION_RECEIPT"),
        "set WH2_RV4A3_VERIFICATION_RECEIPT for receipt binding",
    )
    def test_authenticated_verification_receipt_hash_is_immutable(self):
        path = Path(os.environ["WH2_RV4A3_VERIFICATION_RECEIPT"])
        self.assertEqual(
            erratum.campaign._stable_file_binding(
                path, byte_limit=4096)["sha256"],
            "32ff54529883e5c10d53cdd1340c22ff61dbc01e7fa366ba673595b7d16bffcd",
        )
        receipt, unused_digest = \
            erratum.campaign.read_canonical_json(path)
        self.assertEqual(
            unused_digest,
            erratum.campaign.canonical_sha256(receipt),
        )
        self.assertEqual(
            receipt["verified"]["training_decision_sha256"],
            "ee676d07ce82f71e9252fa2a26b031630470cebc650d9910b098fe83655a7e52",
        )


if __name__ == "__main__":
    unittest.main()
