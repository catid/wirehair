#!/usr/bin/env python3
"""Adversarial synthetic tests for the independent repairtiming v3 parser."""

import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import peel_codec
import repair_timing_codec as repair


def _thermal_row(monotonic, cpu="60.0", dimm="40.0"):
    fields = [
        "2026-07-29T00:00:00Z",
        f"{monotonic:.9f}",
        "100.0",
        "4000.0",
        cpu,
        *(dimm for unused in range(8)),
        "0",
        "1.0",
        "1.0",
        "1.0",
        "0",
        "0",
    ]
    return ",".join(fields) + "\n"


def _thermal_csv(*rows):
    return (
        ",".join(peel_codec._THERMAL_CSV_COLUMNS) + "\n" +
        "".join(rows)
    ).encode("ascii")


def synthetic_context(start=0.5, end=3.0):
    middle = (start + end) / 2.0
    thermal = peel_codec._thermal_window_from_csv(
        _thermal_csv(
            _thermal_row(start),
            _thermal_row(middle),
            _thermal_row(end),
        ),
        int(round(start * 1e9)),
        int(round(end * 1e9)),
        sampling_interval_ms=1000,
    )
    return {
        "schema": peel_codec.PAIRED_CONTEXT_SCHEMA,
        "bound": {
            "cpu_model": "synthetic-test-cpu",
            "cpu_affinity": [2],
            "cpu_governors": {"2": "performance"},
            "cache_state": "warm",
            "evict_bytes": 4096,
            "thermal_source": "/tmp/synthetic-repairtiming-thermal.csv",
            "thermal_sampling_interval_ms": 1000,
            "clock_domain": peel_codec.PAIRED_CLOCK_DOMAIN,
            "thermal_device": 11,
            "thermal_inode": 22,
            "thermal_prelaunch_monotonic_s": float(start),
            "thermal_prelaunch_row_sha256": hashlib.sha256(
                thermal["prelaunch_row"].encode("ascii")
            ).hexdigest(),
        },
        "thermal": thermal,
    }


def _digest(label):
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def _stats(available, *, block_count=17):
    if not available:
        return {name: -1 for name in repair.REPAIRTIMING_STATS_FIELDS}
    values = {name: 0 for name in repair.REPAIRTIMING_STATS_FIELDS}
    values["packet_rows"] = block_count
    values["peeled_columns"] = block_count
    values["block_xors"] = 1
    values["block_copies"] = 1
    values["build_ns"] = 1
    return values


def _role(block_count, max_overhead, payload, message):
    return {
        "encode_result": 0,
        "decode_construct_result": 0,
        "feed_result": 0,
        "recover_result": 0,
        "recovery_ok": 1,
        "encoded_symbols": block_count + max_overhead,
        "received_symbols": block_count,
        "overhead": 0,
        "payload_sha256": payload,
        "recovered_sha256": message,
    }


def synthetic_stream(
        *, block_count=17, block_bytes=2, arm=None,
        construction_seed=0xfedcba98,
        construction_seed_derivation=repair.REPAIRTIMING_CONSTRUCTION_SEED_FIXED,
        loss=0.1, loss_seed=0x0123456789abcdef,
        schedule="iid", warmups=0, replicates=3, inner_reps=1,
        max_overhead=64, context=None,
        selector_result=0,
        started_ns=1_000_000_000, finished_ns=2_000_000_000):
    """Build valid non-roster evidence without executing a campaign outcome."""
    arm = arm or repair.REPAIR_V1_ARMS[0]
    context = context or synthetic_context()
    cells_count = warmups + replicates
    context_sha = peel_codec._canonical_json_sha256(context["bound"])
    manifest = {
        "protocol": repair.REPAIRTIMING_PROTOCOL_V3,
        "schema": repair.REPAIRTIMING_SCHEMA_V3,
        "repair_arm": arm.name,
        "repair_contract_id": f"{arm.provisional_id:016x}",
        "repair_contract_sha256": arm.contract_sha256,
        "repair_policy_sha256": repair.REPAIR_V1_POLICY_SHA256,
        "repair_attempt_cap": 8,
        "descriptor_magic": "W2R3",
        "descriptor_encoding": repair.REPAIR_V1_DESCRIPTOR_ENCODING,
        "descriptor_bytes": 17,
        "descriptor_policy": repair.REPAIR_V1_DESCRIPTOR_POLICY,
        "K": block_count,
        "bb": block_bytes,
        "message_bytes": (block_count - 1) * block_bytes + 1,
        "message_tail_bytes": 1,
        "dispatch_profile": peel_codec.TARGET_PROFILE,
        "dispatch_contract_id": peel_codec.TARGET_CONTRACT["contract_id"],
        "construction_seed_base": construction_seed,
        "construction_seed_derivation": construction_seed_derivation,
        "loss": f"{loss:.17g}",
        "loss_seed_base": loss_seed,
        "loss_seed_derivation": repair.REPAIRTIMING_LOSS_SEED_DERIVATION,
        "schedule": schedule,
        "loss_model": "packet-schedule-v1",
        "trace_encoding": "wirehair-wh2-repairtiming-loss-trace-v1",
        "message_seed_policy": "replicate-loss-seed-partial-final-v1",
        "warmup_replicates": warmups,
        "replicates": replicates,
        "cells": cells_count,
        "summary_rows_per_cell": 1,
        "attempt_rows_per_cell": 8,
        "timing_panels": "+".join(repair.REPAIRTIMING_PANELS),
        "timing_panels_per_cell": len(repair.REPAIRTIMING_PANELS),
        "timing_slots_per_panel": len(repair.REPAIRTIMING_ORDER),
        "timing_rows_per_cell":
            len(repair.REPAIRTIMING_PANELS) *
            len(repair.REPAIRTIMING_ORDER),
        "timing_order": repair.REPAIRTIMING_ORDER,
        "timing_label_policy": "fixed-logical-role-v1",
        "inner_reps": inner_reps,
        "max_overhead": max_overhead,
        "cache_state": "warm",
        "systematic_cache": "off",
        "evict_bytes": 4096,
        "payload_alignment": 64,
        "prefault": 1,
        "cpu_affinity_policy": "first-allowed-affinity-v1",
        "encoder_selector_scope":
            "lazy-selector-init-first-k-encode-w2r3-serialize-v1",
        "encoder_forced_scope":
            "forced-selected-init-first-k-encode-w2r3-serialize-v1",
        "decoder_feed_scope":
            "exact-selected-descriptor-bound-feed-through-first-success-v1",
        "decoder_full_scope":
            "w2r3-parse-create-feed-through-first-success-recover-v1",
        "direct_scope":
            "selected-attempt-pair-local-fixed-prefix-encoder-"
            "intermediate-witnessed-solve-v1",
        "codec_reuse": "none-fresh-object-every-inner-v1",
        "scratch_reporting": "complete-precode-solve-stats-per-call-v1",
        "timing_stats_policy":
            "first-inner-structural-counters-summed-phase-ns-v1",
        "selector_timing_stats_policy":
            "all-executed-selector-subcalls-summed-per-inner-v1",
        "structural_digest_encoding":
            "wirehair-wh2-repair-v1-selector-structural-v1",
        "context_sha256": context_sha,
        "uncertainty": repair.REPAIRTIMING_UNCERTAINTY,
        "required_margin": "0",
        "margin_rule":
            "upper-log-cost-lt-negative-required-margin-and-aa-floors-v1",
        "clock_domain": repair.REPAIRTIMING_CLOCK_DOMAIN,
        "stream_hash_scope": "body-plus-done-prefix-v1",
        "started_monotonic_ns": started_ns,
        "expected_summary_rows": cells_count,
        "expected_attempt_rows": cells_count * 8,
        "expected_timing_rows": cells_count * 112,
        "expected_rows": cells_count * 121,
    }
    lines = [
        "# repairtiming," + ",".join(
            f"{name}={manifest[name]}"
            for name in repair.REPAIRTIMING_MANIFEST_FIELDS
        ) + "\n",
        "# repairtiming_summary_columns," +
        ",".join(repair.REPAIRTIMING_SUMMARY_COLUMNS) + "\n",
        "# repairtiming_attempt_columns," +
        ",".join(repair.REPAIRTIMING_ATTEMPT_COLUMNS) + "\n",
        "# repairtiming_timing_columns," +
        ",".join(repair.REPAIRTIMING_TIMING_COLUMNS) + "\n",
    ]

    for cell_index in range(cells_count):
        measured = int(cell_index >= warmups)
        root = repair._expected_construction_root(
            construction_seed, cell_index, construction_seed_derivation)
        cell_loss = repair._expected_loss_seed(loss_seed, cell_index)
        trace = _digest(f"synthetic-trace-{cell_index}")
        message = _digest(f"synthetic-message-{cell_index}")
        payload = _digest(f"synthetic-payload-{cell_index}")
        repair_intermediate = _digest(
            f"synthetic-repair-intermediate-{cell_index}")
        dispatch_intermediate = _digest(
            f"synthetic-dispatch-intermediate-{cell_index}")
        common = {
            "row_kind": "summary",
            "cell_index": cell_index,
            "measured": measured,
            "construction_root": root,
            "loss_seed": cell_loss,
            "trace_sha256": trace,
            "message_sha256": message,
        }
        committed = selector_result == 0
        descriptor = (
            repair._descriptor_bytes(arm, root, 0)
            if committed else None)
        summary = dict(common)
        summary.update({
            "selector_result": selector_result,
            "selected_attempt": 0 if committed else -1,
            "attempts_executed": 1,
            "calls_executed": 1,
            "real_configuration_calls": 1,
            "structural_probe_calls": 0,
            "cap_exhausted": 0,
            "fatal_attempt_zero_mismatch": 0,
            "oom": int(selector_result == 9),
            "committed": int(committed),
            "descriptor_hex":
                descriptor.hex() if committed else "not_applicable",
            "descriptor_sha256": (
                hashlib.sha256(descriptor).hexdigest()
                if committed else "not_applicable"
            ),
            "selector_structural_sha256": "0" * 64,
        })
        for role in repair.REPAIRTIMING_ROLES:
            role_value = _role(block_count, max_overhead, payload, message)
            if not committed and role in ("raw", "repaired"):
                role_value = {
                    "encode_result": selector_result,
                    "decode_construct_result": -1,
                    "feed_result": -1,
                    "recover_result": -1,
                    "recovery_ok": 0,
                    "encoded_symbols": 0,
                    "received_symbols": 0,
                    "overhead": -1,
                    "payload_sha256": "not_applicable",
                    "recovered_sha256": "not_applicable",
                }
            for suffix, value in role_value.items():
                summary[f"{role}_{suffix}"] = value
        summary.update({
            "forced_a_result": 0 if committed else -1,
            "forced_b_result": 0 if committed else -1,
            "forced_a_payload_sha256":
                payload if committed else "not_applicable",
            "forced_b_payload_sha256":
                payload if committed else "not_applicable",
            "forced_equal": int(committed),
            "direct_fixed_prefix_symbols":
                block_count if committed else block_count + max_overhead,
            "repair_intermediate_bytes": (
                repair._candidate_intermediate_bytes(
                    arm, block_count, block_bytes)
                if committed else -1
            ),
            "repair_intermediate_sha256":
                repair_intermediate if committed else "not_applicable",
            "repair_direct_executed": int(committed),
            "repair_direct_result": 0 if committed else -1,
            "repair_direct_intermediate_bytes": (
                repair._candidate_intermediate_bytes(
                    arm, block_count, block_bytes)
                if committed else -1
            ),
            "repair_direct_intermediate_sha256":
                repair_intermediate if committed else "not_applicable",
            "repair_direct_witness_equal": int(committed),
            "dispatch_intermediate_bytes":
                repair._dispatch_intermediate_bytes(
                    block_count, block_bytes),
            "dispatch_intermediate_sha256": dispatch_intermediate,
            "dispatch_direct_executed": 1,
            "dispatch_direct_result": 0,
            "dispatch_direct_intermediate_bytes":
                repair._dispatch_intermediate_bytes(
                    block_count, block_bytes),
            "dispatch_direct_intermediate_sha256": dispatch_intermediate,
            "dispatch_direct_witness_equal": 1,
        })

        attempt_rows = []
        normalized_attempts = []
        for attempt in range(8):
            row = dict(common)
            row["row_kind"] = "attempt"
            row.update({
                "attempt": attempt,
                "probe_executed": 0,
                "probe_result": -1,
                "probe_stats_available": 0,
                "real_executed": int(attempt == 0),
                "real_result": selector_result if attempt == 0 else -1,
                "real_stats_available": int(attempt == 0 and committed),
            })
            for name, value in _stats(False).items():
                row[f"probe_{name}"] = value
            for name, value in _stats(
                    attempt == 0 and committed,
                    block_count=block_count).items():
                row[f"real_{name}"] = value
            attempt_rows.append(row)
            normalized_attempts.append(repair._parse_attempt_row(
                copy.deepcopy(row),
                expected_attempt=attempt,
                label=f"synthetic attempt {attempt}",
            ))
        summary["selector_structural_sha256"] = hashlib.sha256(
            repair._structural_digest_text(
                arm, block_count, root, summary,
                normalized_attempts).encode("ascii")
        ).hexdigest()

        lines.append(",".join(
            str(summary[name])
            for name in repair.REPAIRTIMING_SUMMARY_COLUMNS
        ) + "\n")
        for row in attempt_rows:
            lines.append(",".join(
                str(row[name])
                for name in repair.REPAIRTIMING_ATTEMPT_COLUMNS
            ) + "\n")

        parsed_roles = repair._validate_summary_row(
            copy.deepcopy(summary),
            arm=arm,
            block_count=block_count,
            block_bytes=block_bytes,
            max_overhead=max_overhead,
            label="synthetic summary",
        )
        for panel_index, (
                panel, unused_scope, left, right
        ) in enumerate(repair.REPAIRTIMING_PANEL_SPECS):
            for slot, timing_label in enumerate(repair.REPAIRTIMING_ORDER):
                role = left if timing_label == "A" else right
                scope = repair._timing_scope_for_role(role)
                expected = repair._timing_role_summary(
                    role, parsed_roles, summary)
                requires_repair = (
                    left.startswith("repair_") or
                    right.startswith("repair_"))
                panel_success = all(
                    repair._timing_role_summary(
                        panel_role, parsed_roles, summary)["success"]
                    for panel_role in (left, right)
                )
                execute = (
                    (committed or not requires_repair) and panel_success)
                available = int(not role.startswith("wh1_"))
                intermediate_bytes, intermediate_sha = \
                    repair._timing_intermediate(
                        role, arm, block_count, block_bytes, summary)
                row = dict(common)
                row["row_kind"] = "timing"
                row.update({
                    "timing_panel": panel,
                    "timing_panel_index": panel_index,
                    "timing_slot": slot,
                    "timing_pair": slot // 2,
                    "timing_label": timing_label,
                    "timing_role": role,
                    "timing_scope": scope,
                    "timing_eligible": int(execute),
                    "timing_executed": int(execute),
                    "timing_censor_reason": (
                        "none" if execute else
                        (
                            "selector_uncommitted"
                            if requires_repair and not committed
                            else "role_preflight_failed"
                        )
                    ),
                    "timing_construct_result":
                        expected["construct_result"] if execute else -1,
                    "timing_result": expected["result"] if execute else -1,
                    "timing_recover_result":
                        expected["recover_result"] if execute else -1,
                    "timing_recovery_ok":
                        expected["recovery_ok"] if execute else -1,
                    "timing_stable": 1 if execute else -1,
                    "timing_elapsed_ns": 1000 if execute else -1,
                    "timing_inner_reps": inner_reps if execute else -1,
                    "timing_saturated": 0 if execute else -1,
                    "timing_cpu_before": 2 if execute else -1,
                    "timing_cpu_after": 2 if execute else -1,
                    "timing_cpu_migrated": 0 if execute else -1,
                    "timing_minflt": 0 if execute else -1,
                    "timing_majflt": 0 if execute else -1,
                    "timing_fault_contaminated": 0 if execute else -1,
                    "timing_stats_available": available if execute else -1,
                    "timing_fixed_prefix_symbols":
                        (
                            summary["direct_fixed_prefix_symbols"]
                            if execute and scope == "direct" else -1
                        ),
                    "timing_intermediate_bytes":
                        intermediate_bytes if execute else -1,
                    "timing_intermediate_sha256":
                        intermediate_sha if execute else "not_applicable",
                })
                for name, value in _stats(
                        execute and available,
                        block_count=block_count).items():
                    row[f"timing_{name}"] = value
                lines.append(",".join(
                    str(row[name])
                    for name in repair.REPAIRTIMING_TIMING_COLUMNS
                ) + "\n")

    done_prefix = (
        "# repairtiming_done,"
        f"complete=1,cells={cells_count},"
        f"summary_rows={cells_count},"
        f"attempt_rows={cells_count * 8},"
        f"timing_rows={cells_count * 112},"
        f"rows={cells_count * 121},"
        f"finished_monotonic_ns={finished_ns},"
        "stream_sha256="
    )
    stream_sha = hashlib.sha256(
        ("".join(lines) + done_prefix).encode("ascii")
    ).hexdigest()
    lines.append(done_prefix + stream_sha + "\n")
    return "".join(lines), context


def _rehash(lines):
    done_prefix = lines[-1].rsplit("stream_sha256=", 1)[0]
    lines[-1] = done_prefix + "stream_sha256=" + hashlib.sha256(
        ("".join(lines[:-1]) + done_prefix + "stream_sha256=")
        .encode("ascii")
    ).hexdigest() + "\n"
    return "".join(lines)


def _mutate_csv_row(stdout, line_index, column, value, columns):
    lines = stdout.splitlines(keepends=True)
    fields = next(csv.reader([lines[line_index][:-1]], strict=True))
    fields[columns.index(column)] = str(value)
    lines[line_index] = ",".join(fields) + "\n"
    return _rehash(lines)


def parse_kwargs(context):
    return {
        "block_count": 17,
        "block_bytes": 2,
        "repair_arm": repair.REPAIR_V1_ARMS[0],
        "construction_seed": 0xfedcba98,
        "construction_seed_derivation":
            repair.REPAIRTIMING_CONSTRUCTION_SEED_FIXED,
        "loss": 0.1,
        "loss_seed": 0x0123456789abcdef,
        "schedule": "iid",
        "warmup_replicates": 0,
        "replicates": 3,
        "inner_reps": 1,
        "max_overhead": 64,
        "cache_state": "warm",
        "systematic_cache": "off",
        "evict_bytes": 4096,
        "context": context,
        "required_margin": 0.0,
    }


class RepairTimingRoundTripTests(unittest.TestCase):
    def setUp(self):
        self.stdout, self.context = synthetic_stream()
        self.kwargs = parse_kwargs(self.context)

    def parse(self, stdout=None):
        return repair.parse_repairtiming_output(
            self.stdout if stdout is None else stdout, **self.kwargs)

    def test_roundtrip_and_compact_persisted_receipt(self):
        evidence = self.parse()
        self.assertEqual(len(evidence.cells), 3)
        cell = evidence.cells[0]
        self.assertEqual(
            cell.selector_key,
            ("19cccf775ce0bf09", 17, 0xfedcba98))
        self.assertEqual(len(cell.selector_witness["attempts"]), 8)
        self.assertNotIn(
            "build_ns",
            cell.selector_witness["attempts"][0]["probe_stats"])
        self.assertEqual(
            len(cell.real_witness["timing_rows"]), 14 * 8)
        receipt = evidence.as_dict()
        self.assertNotIn("cells", receipt)
        replayed = repair.replay_repairtiming_receipt(
            receipt, expected_request=self.kwargs)
        self.assertEqual(
            repair.selector_projection(replayed.cells[0]),
            repair.selector_projection(cell))
        self.assertEqual(
            repair.real_projection(replayed.cells[0]),
            repair.real_projection(cell))

    def test_typed_header_hashes_match_native_contract(self):
        expected = (
            (
                "summary",
                repair.REPAIRTIMING_SUMMARY_COLUMNS,
                "0f7ba631537c4f4fd69c81966e9a684b"
                "ed2731286bdb74d40fe0d37291549667",
            ),
            (
                "attempt",
                repair.REPAIRTIMING_ATTEMPT_COLUMNS,
                "2ea38a4a483b96bf591ec974f402863a"
                "25226886c6ab3564ac255b8121d3c842",
            ),
            (
                "timing",
                repair.REPAIRTIMING_TIMING_COLUMNS,
                "631e51b68dd213652f24f15f17ef6840"
                "68943d8acedf7333d197878e0c0d700c",
            ),
        )
        for kind, columns, digest in expected:
            with self.subTest(kind=kind):
                header = (
                    f"# repairtiming_{kind}_columns," +
                    ",".join(columns) + "\n"
                )
                self.assertEqual(
                    hashlib.sha256(header.encode("ascii")).hexdigest(),
                    digest)

    def test_uncommitted_selector_keeps_independent_control_aa(self):
        stdout, context = synthetic_stream(selector_result=4)
        parsed = repair.parse_repairtiming_output(
            stdout, **parse_kwargs(context))
        cell = parsed.cells[0]
        self.assertEqual(cell.selector_witness["committed"], 0)
        self.assertEqual(
            cell.real_witness["controls"]["forced_a_result"], -1)
        columns = cell.real_witness["timing_columns"]
        rows = cell.real_witness["timing_rows"]
        panel_index = columns.index("timing_panel")
        executed_index = columns.index("timing_executed")
        reason_index = columns.index("timing_censor_reason")
        by_panel = {}
        for row in rows:
            by_panel.setdefault(row[panel_index], []).append(row)
        for panel in (
                "encoder_wh1_aa", "decoder_feed_wh1_aa",
                "decoder_full_wh1_aa", "direct_dispatch_aa"):
            self.assertEqual(
                {row[executed_index] for row in by_panel[panel]}, {1})
            self.assertEqual(
                {row[reason_index] for row in by_panel[panel]}, {"none"})
        self.assertEqual(
            {
                row[reason_index]
                for row in by_panel["encoder_selector_aa"]
            },
            {"selector_uncommitted"},
        )

    def test_schema_derived_caps_cover_fixture_without_duplicate_cells(self):
        evidence = self.parse()
        receipt = evidence.as_dict()
        canonical = json.dumps(
            receipt, sort_keys=True, separators=(",", ":")).encode("ascii")
        self.assertLess(
            len(self.stdout),
            repair.REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT)
        self.assertLess(
            len(canonical),
            repair.REPAIRTIMING_RECEIPT_CANONICAL_BYTE_LIMIT)
        self.assertLessEqual(
            repair.REPAIRTIMING_RECEIPT_CELL_CANONICAL_BYTE_LIMIT,
            256 * 1024)

    def test_dimension_boundaries_and_type_aliases(self):
        base = dict(self.kwargs)
        del base["context"]
        for name, value in (
            ("block_count", 1),
            ("block_count", 101),
            ("block_bytes", 1),
            ("block_bytes", 4098),
            ("max_overhead", 63),
            ("max_overhead", 65),
            ("replicates", 0),
            ("replicates", 2),
            ("replicates", 769),
            ("warmup_replicates", 766),
            ("inner_reps", 0),
            ("inner_reps", 1025),
            ("construction_seed", True),
            ("construction_seed_derivation", []),
            ("loss", -0.0),
            ("loss", -0.001),
            ("loss", 0.9900000000001),
            ("loss", 1.0),
            ("loss", float("nan")),
            ("loss", float("inf")),
            ("schedule", "unknown"),
            ("evict_bytes", 4095),
            (
                "evict_bytes",
                peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1,
            ),
            ("required_margin", -0.0),
            ("required_margin", 1.001),
            ("required_margin", float("nan")),
        ):
            with self.subTest(name=name, value=value):
                changed = dict(base)
                changed[name] = value
                with self.assertRaises(ValueError):
                    repair.validate_repairtiming_dimensions(**changed)
        valid = dict(base)
        valid.update({
            "block_count": 2,
            "block_bytes": 4096,
            "construction_seed": 0xffffffff,
            "loss": 0.99,
            "loss_seed": 0xffffffffffffffff,
            "replicates": 3,
            "inner_reps": 1024,
            "evict_bytes": peel_codec.PEELTIMING_EVICT_BYTES_MAX,
            "required_margin": 1.0,
        })
        repair.validate_repairtiming_dimensions(**valid)
        valid.update({
            "block_count": 100,
            "block_bytes": 2,
            "cache_state": "cold",
            "inner_reps": 1,
        })
        repair.validate_repairtiming_dimensions(**valid)
        boundary_stdout, boundary_context = synthetic_stream(loss=0.99)
        boundary_kwargs = parse_kwargs(boundary_context)
        boundary_kwargs["loss"] = 0.99
        self.assertEqual(
            repair.parse_repairtiming_output(
                boundary_stdout, **boundary_kwargs).manifest["loss"],
            0.99,
        )
        over_boundary_lines = boundary_stdout.splitlines(keepends=True)
        over_boundary_lines[0] = over_boundary_lines[0].replace(
            f"loss={0.99:.17g},",
            f"loss={0.9900000000001:.17g},",
            1,
        )
        with self.assertRaises(repair.MeasurementError):
            repair.parse_repairtiming_output(
                _rehash(over_boundary_lines), **boundary_kwargs)
        self.assertEqual(
            repair._dispatch_intermediate_bytes(16, 2), 74)
        self.assertEqual(
            repair._candidate_intermediate_bytes(
                repair.REPAIR_V1_ARMS[0], 2, 2),
            28)


class RepairTimingMutationTests(RepairTimingRoundTripTests):
    def assertRejected(self, stdout):
        with self.assertRaises(repair.MeasurementError):
            self.parse(stdout)

    def test_wrong_manifest_protocol_schema_contract_and_policy_rejected(self):
        for name, value in (
            ("protocol", "wirehair-v2-bench:repairtiming:repair-v1:v2"),
            ("schema", "wirehair.wh2.repairtiming.v2"),
            ("repair_contract_id", "0000000000000000"),
            ("repair_contract_sha256", "0" * 64),
            ("repair_policy_sha256", "0" * 64),
            ("selector_timing_stats_policy", "selected-real-only-v1"),
            ("timing_panels", "encoder_selector_forced"),
        ):
            with self.subTest(name=name):
                lines = self.stdout.splitlines(keepends=True)
                fields = lines[0][len("# repairtiming,"):-1].split(",")
                index = repair.REPAIRTIMING_MANIFEST_FIELDS.index(name)
                fields[index] = f"{name}={value}"
                lines[0] = "# repairtiming," + ",".join(fields) + "\n"
                self.assertRejected(_rehash(lines))

    def test_typed_headers_missing_reordered_or_unionized_rejected(self):
        for mutate in ("missing", "reordered", "union"):
            with self.subTest(mutate=mutate):
                lines = self.stdout.splitlines(keepends=True)
                if mutate == "missing":
                    del lines[2]
                elif mutate == "reordered":
                    fields = lines[2].split(",")
                    fields[-1], fields[-2] = fields[-2], fields[-1]
                    lines[2] = ",".join(fields)
                else:
                    lines[1] = (
                        "# repairtiming_summary_columns," +
                        ",".join(repair.REPAIRTIMING_COLUMNS) + "\n")
                self.assertRejected(_rehash(lines))

    def test_missing_duplicate_and_reordered_attempts_rejected(self):
        for mutate in ("missing", "duplicate", "reordered"):
            with self.subTest(mutate=mutate):
                lines = self.stdout.splitlines(keepends=True)
                if mutate == "missing":
                    del lines[6]
                elif mutate == "duplicate":
                    lines.insert(6, lines[5])
                else:
                    lines[5], lines[6] = lines[6], lines[5]
                self.assertRejected(_rehash(lines))

    def test_impossible_policy_accounting_and_stats_rejected(self):
        # Data starts after four preamble lines. Summary=4, attempts=5..12.
        mutations = (
            (4, "calls_executed", 2, repair.REPAIRTIMING_SUMMARY_COLUMNS),
            (4, "selected_attempt", 1, repair.REPAIRTIMING_SUMMARY_COLUMNS),
            (4, "selected_attempt", -1, repair.REPAIRTIMING_SUMMARY_COLUMNS),
            (5, "real_executed", 0, repair.REPAIRTIMING_ATTEMPT_COLUMNS),
            (
                6, "probe_stats_available", 1,
                repair.REPAIRTIMING_ATTEMPT_COLUMNS,
            ),
            (
                5, "real_residual_rank", 2,
                repair.REPAIRTIMING_ATTEMPT_COLUMNS,
            ),
            (
                5, "real_stats_available", 0,
                repair.REPAIRTIMING_ATTEMPT_COLUMNS,
            ),
        )
        for line, field, value, columns in mutations:
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    self.stdout, line, field, value, columns))

    def test_descriptor_structural_digest_and_common_witness_rejected(self):
        for line, field, value, columns in (
            (
                4, "descriptor_hex", "00" * 17,
                repair.REPAIRTIMING_SUMMARY_COLUMNS,
            ),
            (
                4, "selector_structural_sha256", "0" * 64,
                repair.REPAIRTIMING_SUMMARY_COLUMNS,
            ),
            (
                5, "trace_sha256", _digest("changed"),
                repair.REPAIRTIMING_ATTEMPT_COLUMNS,
            ),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    self.stdout, line, field, value, columns))

    def test_stage_direct_prefix_and_intermediate_drift_rejected(self):
        for field, value in (
            ("raw_encoded_symbols", 17),
            ("raw_decode_construct_result", -1),
            ("raw_recovery_ok", 0),
            ("raw_received_symbols", 18),
            ("direct_fixed_prefix_symbols", 18),
            ("repair_direct_executed", 0),
            ("repair_direct_intermediate_sha256", _digest("drift")),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    self.stdout, 4, field, value,
                    repair.REPAIRTIMING_SUMMARY_COLUMNS))

    def test_failed_direct_uses_executed_failure_sentinels(self):
        lines = self.stdout.splitlines(keepends=True)
        fields = next(csv.reader([lines[4][:-1]], strict=True))
        for name, value in (
            ("repair_direct_result", 4),
            ("repair_direct_intermediate_bytes", -1),
            ("repair_direct_intermediate_sha256", "not_applicable"),
            ("repair_direct_witness_equal", 0),
        ):
            fields[repair.REPAIRTIMING_SUMMARY_COLUMNS.index(name)] = \
                str(value)
        lines[4] = ",".join(fields) + "\n"

        # A failed summary direct preflight censors every repair-direct panel
        # row without claiming any measured work.
        first_timing = 13
        for direct_panel in (11, 12):
            for offset in range(8):
                row_index = first_timing + direct_panel * 8 + offset
                fields = next(csv.reader(
                    [lines[row_index][:-1]], strict=True))
                for name, value in (
                    ("timing_eligible", 0),
                    ("timing_executed", 0),
                    ("timing_censor_reason", "role_preflight_failed"),
                    ("timing_construct_result", -1),
                    ("timing_result", -1),
                    ("timing_recover_result", -1),
                    ("timing_recovery_ok", -1),
                    ("timing_stable", -1),
                    ("timing_elapsed_ns", -1),
                    ("timing_inner_reps", -1),
                    ("timing_saturated", -1),
                    ("timing_cpu_before", -1),
                    ("timing_cpu_after", -1),
                    ("timing_cpu_migrated", -1),
                    ("timing_minflt", -1),
                    ("timing_majflt", -1),
                    ("timing_fault_contaminated", -1),
                    ("timing_stats_available", -1),
                    ("timing_fixed_prefix_symbols", -1),
                    ("timing_intermediate_bytes", -1),
                    ("timing_intermediate_sha256", "not_applicable"),
                ):
                    fields[
                        repair.REPAIRTIMING_TIMING_COLUMNS.index(name)
                    ] = str(value)
                for name in repair.REPAIRTIMING_STATS_FIELDS:
                    fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                        f"timing_{name}")] = "-1"
                lines[row_index] = ",".join(fields) + "\n"
        failed = _rehash(lines)
        parsed = self.parse(failed)
        controls = parsed.cells[0].real_witness["controls"]
        self.assertEqual(controls["repair_direct_executed"], 1)
        self.assertEqual(controls["repair_direct_result"], 4)
        self.assertEqual(
            controls["repair_direct_intermediate_bytes"], -1)

        for field, value in (
            ("repair_direct_intermediate_bytes", 1),
            ("repair_direct_intermediate_sha256", _digest("partial")),
            ("repair_direct_witness_equal", 1),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    failed, 4, field, value,
                    repair.REPAIRTIMING_SUMMARY_COLUMNS))

    def test_timing_order_fast_failure_and_fake_eligibility_rejected(self):
        first_timing = 13
        for field, value in (
            ("timing_panel_index", 1),
            ("timing_slot", 1),
            ("timing_role", "wh1_encoder"),
            ("timing_elapsed_ns", 0),
            ("timing_stats_available", 0),
            ("timing_block_xors", 2),
            ("timing_fixed_prefix_symbols", 17),
            ("timing_intermediate_sha256", _digest("drift")),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    self.stdout, first_timing, field, value,
                    repair.REPAIRTIMING_TIMING_COLUMNS))

    def test_selector_timing_uses_call_aware_not_inner_aware_bounds(self):
        row = {
            f"timing_{name}": value
            for name, value in _stats(True).items()
        }
        calls = repair.REPAIRTIMING_SELECTOR_CALL_CAP
        row["timing_packet_seed_attempt"] = 255 * calls
        row["timing_inactivated_columns"] = 0xffff * calls
        row["timing_mixed_joint_scratch_bytes"] = (
            peel_codec.PEELTIMING_WORKING_SET_BYTES_MAX * calls
        )
        repair._validate_stats(
            row, "timing", 1, "selector aggregate",
            acceptance_max=calls, call_limit=calls)
        for field in (
                "timing_packet_seed_attempt",
                "timing_inactivated_columns",
                "timing_mixed_joint_scratch_bytes"):
            with self.subTest(field=field):
                invalid = dict(row)
                invalid[field] += 1
                with self.assertRaises(repair.MeasurementError):
                    repair._validate_stats(
                        invalid, "timing", 1, "selector aggregate",
                        acceptance_max=calls, call_limit=calls)

    def test_failed_timed_call_retains_telemetry_not_partial_witness(self):
        lines = self.stdout.splitlines(keepends=True)
        first_timing = 13
        for offset in range(8):
            fields = next(csv.reader(
                [lines[first_timing + offset][:-1]], strict=True))
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_eligible")] = "0"
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_censor_reason")] = "outcome_unstable"
            if offset == 0:
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_result")] = "4"
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_stable")] = "0"
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_intermediate_bytes")] = "-1"
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_intermediate_sha256")] = "not_applicable"
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_stats_available")] = "0"
                for name in repair.REPAIRTIMING_STATS_FIELDS:
                    fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                        f"timing_{name}")] = "-1"
            lines[first_timing + offset] = ",".join(fields) + "\n"
        failed = _rehash(lines)
        cell = self.parse(failed).cells[0]
        columns = cell.real_witness["timing_columns"]
        row = cell.real_witness["timing_rows"][0]
        self.assertEqual(row[columns.index("timing_executed")], 1)
        self.assertEqual(row[columns.index("timing_result")], 4)
        self.assertEqual(row[columns.index("timing_elapsed_ns")], 1000)
        self.assertEqual(
            row[columns.index("timing_intermediate_bytes")], -1)

        for field, value in (
            (
                "timing_intermediate_bytes",
                repair._candidate_intermediate_bytes(
                    repair.REPAIR_V1_ARMS[0], 17, 2),
            ),
            ("timing_intermediate_sha256", _digest("partial")),
            ("timing_stable", 1),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    failed, first_timing, field, value,
                    repair.REPAIRTIMING_TIMING_COLUMNS))

    def test_failed_full_decode_uses_recover_and_stats_sentinels(self):
        lines = self.stdout.splitlines(keepends=True)
        first_timing = 13
        full_panel = 8
        panel_start = first_timing + full_panel * 8
        for offset in range(8):
            fields = next(csv.reader(
                [lines[panel_start + offset][:-1]], strict=True))
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_eligible")] = "0"
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_censor_reason")] = "outcome_unstable"
            if offset == 0:
                for name, value in (
                    ("timing_result", 4),
                    ("timing_recover_result", -1),
                    ("timing_recovery_ok", 0),
                    ("timing_stable", 0),
                    ("timing_stats_available", 0),
                ):
                    fields[
                        repair.REPAIRTIMING_TIMING_COLUMNS.index(name)
                    ] = str(value)
                for name in repair.REPAIRTIMING_STATS_FIELDS:
                    fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                        f"timing_{name}")] = "-1"
            lines[panel_start + offset] = ",".join(fields) + "\n"
        failed = _rehash(lines)
        cell = self.parse(failed).cells[0]
        columns = cell.real_witness["timing_columns"]
        row = cell.real_witness["timing_rows"][full_panel * 8]
        self.assertEqual(row[columns.index("timing_result")], 4)
        self.assertEqual(row[columns.index("timing_recover_result")], -1)
        self.assertEqual(row[columns.index("timing_stats_available")], 0)

        for field, value in (
            ("timing_recover_result", 0),
            ("timing_recovery_ok", 1),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    failed, panel_start, field, value,
                    repair.REPAIRTIMING_TIMING_COLUMNS))

    def test_failed_timed_direct_keeps_stats_but_not_partial_witness(self):
        lines = self.stdout.splitlines(keepends=True)
        first_timing = 13
        direct_panel = 11
        panel_start = first_timing + direct_panel * 8
        for offset in range(8):
            fields = next(csv.reader(
                [lines[panel_start + offset][:-1]], strict=True))
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_eligible")] = "0"
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_censor_reason")] = "outcome_unstable"
            if offset == 0:
                for name, value in (
                    ("timing_result", 4),
                    ("timing_stable", 0),
                    ("timing_intermediate_bytes", -1),
                    (
                        "timing_intermediate_sha256",
                        "not_applicable",
                    ),
                ):
                    fields[
                        repair.REPAIRTIMING_TIMING_COLUMNS.index(name)
                    ] = str(value)
            lines[panel_start + offset] = ",".join(fields) + "\n"
        failed = _rehash(lines)
        cell = self.parse(failed).cells[0]
        columns = cell.real_witness["timing_columns"]
        row = cell.real_witness["timing_rows"][direct_panel * 8]
        self.assertEqual(row[columns.index("timing_result")], 4)
        self.assertEqual(row[columns.index("timing_stats_available")], 1)
        self.assertEqual(
            row[columns.index("timing_intermediate_bytes")], -1)

        self.assertRejected(_mutate_csv_row(
            failed, panel_start, "timing_intermediate_bytes",
            repair._candidate_intermediate_bytes(
                repair.REPAIR_V1_ARMS[0], 17, 2),
            repair.REPAIRTIMING_TIMING_COLUMNS))

    def test_panel_wide_contamination_retains_work_and_rejects_wrong_reason(self):
        lines = self.stdout.splitlines(keepends=True)
        first_timing = 13
        # Contaminate one slot; every slot in the panel must carry the
        # panel-wide highest-priority reason and eligible=0.
        for offset in range(8):
            fields = next(csv.reader(
                [lines[first_timing + offset][:-1]], strict=True))
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_eligible")] = "0"
            fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                "timing_censor_reason")] = "cpu_migrated"
            if offset == 0:
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_cpu_after")] = "3"
                fields[repair.REPAIRTIMING_TIMING_COLUMNS.index(
                    "timing_cpu_migrated")] = "1"
            lines[first_timing + offset] = ",".join(fields) + "\n"
        contaminated = _rehash(lines)
        parsed = self.parse(contaminated)
        timing_columns = parsed.cells[0].real_witness["timing_columns"]
        row = parsed.cells[0].real_witness["timing_rows"][0]
        self.assertEqual(
            row[timing_columns.index("timing_elapsed_ns")], 1000)
        self.assertEqual(
            row[timing_columns.index("timing_censor_reason")],
            "cpu_migrated")

        wrong = _mutate_csv_row(
            contaminated, first_timing + 1,
            "timing_censor_reason", "fault_contaminated",
            repair.REPAIRTIMING_TIMING_COLUMNS)
        self.assertRejected(wrong)

    def test_noncanonical_numeric_csv_and_stream_hash_rejected(self):
        for field, value in (
            ("cell_index", "+0"),
            ("cell_index", "00"),
            ("measured", "2"),
        ):
            with self.subTest(field=field):
                self.assertRejected(_mutate_csv_row(
                    self.stdout, 4, field, value,
                    repair.REPAIRTIMING_SUMMARY_COLUMNS))
        quoted = self.stdout.replace(
            "summary,0,1,", '"summary",0,1,', 1)
        self.assertRejected(quoted)
        lines = self.stdout.splitlines(keepends=True)
        lines[-1] = lines[-1].replace(
            "stream_sha256=", "stream_sha256=0", 1)
        self.assertRejected("".join(lines))

    def test_receipt_replay_rejects_derived_field_and_request_forgery(self):
        receipt = self.parse().as_dict()
        forged = copy.deepcopy(receipt)
        forged["evidence"]["rows"] += 1
        with self.assertRaises(repair.MeasurementError):
            repair.replay_repairtiming_receipt(forged)
        with self.assertRaises(repair.MeasurementError):
            repair.replay_repairtiming_receipt(
                receipt, expected_request={"max_overhead": 63})

    def test_stdout_byte_cap_checked_before_line_parsing(self):
        with mock.patch.object(
                repair, "REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT",
                len(self.stdout)):
            self.parse()
        with mock.patch.object(
                repair, "REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT",
                len(self.stdout) - 1):
            self.assertRejected(self.stdout)

    def test_context_canonical_byte_cap_exact_boundary(self):
        size = len(json.dumps(
            self.context, sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False).encode("ascii"))
        with mock.patch.object(
                repair, "REPAIRTIMING_CONTEXT_CANONICAL_BYTE_LIMIT",
                size):
            self.parse()
        with mock.patch.object(
                repair, "REPAIRTIMING_CONTEXT_CANONICAL_BYTE_LIMIT",
                size - 1):
            with self.assertRaises(repair.MeasurementError):
                self.parse()


class RepairTimingProcessCaptureTests(unittest.TestCase):
    @staticmethod
    def run_script(
            source, *, stdout_limit=4096, stderr_limit=4096,
            diagnostic_limit=256):
        return repair._run_repairtiming_checked(
            [sys.executable, "-c", source],
            peel_codec.isolated_codec_env(),
            stdout_limit=stdout_limit,
            stderr_limit=stderr_limit,
            diagnostic_limit=diagnostic_limit,
        )

    def assert_process_absent(self, process_id, *, process_group=False):
        deadline = time.monotonic() + 5.0
        while True:
            try:
                if process_group:
                    os.killpg(process_id, 0)
                else:
                    os.kill(process_id, 0)
            except ProcessLookupError:
                return
            if time.monotonic() >= deadline:
                self.fail(f"process {process_id} survived cleanup")
            time.sleep(0.01)

    def wait_for_path(self, path, process=None):
        deadline = time.monotonic() + 5.0
        while not path.exists():
            if process is not None and process.poll() is not None:
                self.fail(
                    f"process exited {process.returncode} before {path} existed")
            if time.monotonic() >= deadline:
                self.fail(f"timed out waiting for {path}")
            time.sleep(0.01)

    def test_success_has_devnull_stdin_isolated_env_and_private_group(self):
        source = "\n".join((
            "import os",
            "import sys",
            "assert sys.stdin.buffer.read() == b''",
            "assert 'WIREHAIR_V2_CAPTURE_SENTINEL' not in os.environ",
            "assert not [",
            "    name for name in os.listdir('/proc/self/fd')",
            "    if name.isdecimal() and int(name) > 2 and",
            "    os.path.exists('/proc/self/fd/' + name)",
            "]",
            "sys.stdout.write(",
            "    f'{os.getpid()} {os.getpgrp()} {os.getsid(0)}')",
        ))
        with mock.patch.dict(
                os.environ,
                {"WIREHAIR_V2_CAPTURE_SENTINEL": "must-not-leak"}):
            stdout = self.run_script(source)
        process_id, process_group, session_id = map(int, stdout.split())
        self.assertNotEqual(process_id, process_group)
        self.assertEqual(session_id, os.getsid(0))
        self.assertNotEqual(process_group, os.getpgrp())
        self.assert_process_absent(process_group, process_group=True)

    def test_exact_stdout_cap_and_success_stderr_rejection(self):
        self.assertEqual(
            self.run_script(
                "import os; os.write(1, b'x' * 32)",
                stdout_limit=32),
            "x" * 32,
        )
        with self.assertRaisesRegex(
                repair.MeasurementError, "stdout exceeds"):
            self.run_script(
                "import os; os.write(1, b'x' * 33)",
                stdout_limit=32)
        with self.assertRaisesRegex(
                repair.MeasurementError, "unexpected stderr"):
            self.run_script(
                "import os; os.write(2, b'unexpected warning')")

    def test_exec_failure_releases_and_reaps_ready_watchdog(self):
        watchdogs = []
        original = repair._start_repairtiming_watchdog

        def traced_start(env):
            watchdog = original(env)
            watchdogs.append(watchdog)
            return watchdog

        with mock.patch.object(
                repair, "_start_repairtiming_watchdog",
                side_effect=traced_start):
            with self.assertRaisesRegex(
                    repair.MeasurementError, "could not execute"):
                repair._run_repairtiming_checked(
                    ["/definitely/missing/repairtiming-benchmark"],
                    peel_codec.isolated_codec_env(),
                    stdout_limit=4096,
                    stderr_limit=4096,
                )
        self.assertEqual(len(watchdogs), 1)
        watchdog = watchdogs[0]
        self.assertEqual(watchdog.liveness_fd, -1)
        self.assertIsNotNone(watchdog.process.returncode)
        self.assert_process_absent(
            watchdog.process_group, process_group=True)

    def test_launch_interrupt_still_reaps_ready_watchdog(self):
        watchdog_processes = []
        original = subprocess.Popen

        def interrupted_launch(*args, **kwargs):
            if kwargs.get("process_group") == 0:
                process = original(*args, **kwargs)
                watchdog_processes.append(process)
                return process
            raise KeyboardInterrupt

        with mock.patch.object(
                repair.subprocess, "Popen", side_effect=interrupted_launch):
            with self.assertRaises(KeyboardInterrupt):
                self.run_script("raise AssertionError('not launched')")
        self.assertEqual(len(watchdog_processes), 1)
        watchdog = watchdog_processes[0]
        self.assertIsNotNone(watchdog.returncode)
        self.assert_process_absent(watchdog.pid, process_group=True)

    def test_stdout_flood_is_killed_reaped_and_bounded_during_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            pid_path = Path(temporary) / "child.pid"
            source = "\n".join((
                "import os",
                "import signal",
                "from pathlib import Path",
                f"Path({str(pid_path)!r}).write_text(str(os.getpid()))",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "while True:",
                "    os.write(1, b'x' * 4096)",
            ))
            started = time.monotonic()
            with self.assertRaisesRegex(
                    repair.MeasurementError, "stdout exceeds"):
                self.run_script(source, stdout_limit=31)
            elapsed = time.monotonic() - started
            self.assertLess(elapsed, 5.0)
            child_pid = int(pid_path.read_text())
            self.assert_process_absent(child_pid, process_group=True)

    def test_post_kill_group_survivor_fails_closed_after_bounded_confirmation(
            self):
        with tempfile.TemporaryDirectory() as temporary:
            identity_path = Path(temporary) / "survivor-target.identity"
            source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                f"Path({str(identity_path)!r}).write_text(",
                "    f'{os.getpid()} {os.getpgrp()}')",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "while True:",
                "    os.write(1, b'x' * 4096)",
            ))
            original = repair._repairtiming_process_group_members

            def injected_survivor(process_group):
                members = original(process_group)
                members[999999999] = "D"
                return members

            with mock.patch.object(
                    repair, "_repairtiming_process_group_members",
                    side_effect=injected_survivor), mock.patch.object(
                        repair, "REPAIRTIMING_GROUP_CONFIRM_SECONDS", 0.05):
                with self.assertRaisesRegex(
                        repair.MeasurementError,
                        "process group survived SIGKILL"):
                    self.run_script(source, stdout_limit=31)
            process_id, process_group = map(
                int, identity_path.read_text().split())
            self.assert_process_absent(process_id)
            self.assert_process_absent(process_group, process_group=True)

    def test_stderr_flood_is_killed_promptly(self):
        source = "\n".join((
            "import os",
            "import time",
            "os.write(2, b'e' * 4096)",
            "time.sleep(60)",
        ))
        started = time.monotonic()
        with self.assertRaisesRegex(
                repair.MeasurementError, "stderr exceeds"):
            self.run_script(source, stderr_limit=31)
        self.assertLess(time.monotonic() - started, 5.0)

    def test_descendant_writer_is_killed_with_private_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            child_pid_path = Path(temporary) / "descendant.pid"
            child_source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                f"Path({str(child_pid_path)!r}).write_text(str(os.getpid()))",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "while True:",
                "    os.write(1, b'd' * 4096)",
            ))
            leader_source = "\n".join((
                "import subprocess",
                "import sys",
                "import time",
                f"subprocess.Popen([sys.executable, '-c', {child_source!r}])",
                "time.sleep(60)",
            ))
            with self.assertRaisesRegex(
                    repair.MeasurementError, "stdout exceeds"):
                self.run_script(leader_source, stdout_limit=31)
            child_pid = int(child_pid_path.read_text())
            self.assert_process_absent(child_pid)

    def test_exited_leader_with_silent_pipe_holder_does_not_deadlock(self):
        with tempfile.TemporaryDirectory() as temporary:
            child_pid_path = Path(temporary) / "holder.pid"
            child_source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                "import time",
                f"Path({str(child_pid_path)!r}).write_text(str(os.getpid()))",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "time.sleep(60)",
            ))
            leader_source = "\n".join((
                "from pathlib import Path",
                "import subprocess",
                "import sys",
                "import time",
                f"path = Path({str(child_pid_path)!r})",
                f"subprocess.Popen([sys.executable, '-c', {child_source!r}])",
                "while not path.exists():",
                "    time.sleep(0.001)",
            ))
            started = time.monotonic()
            with self.assertRaisesRegex(
                    repair.MeasurementError, "held an output pipe open"):
                self.run_script(leader_source)
            self.assertLess(time.monotonic() - started, 5.0)
            child_pid = int(child_pid_path.read_text())
            self.assert_process_absent(child_pid)

    def test_exited_leader_with_detached_pipes_reports_and_kills_descendant(self):
        with tempfile.TemporaryDirectory() as temporary:
            child_pid_path = Path(temporary) / "detached-holder.pid"
            child_source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                "import time",
                f"Path({str(child_pid_path)!r}).write_text(str(os.getpid()))",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "time.sleep(60)",
            ))
            leader_source = "\n".join((
                "from pathlib import Path",
                "import subprocess",
                "import sys",
                "import time",
                f"path = Path({str(child_pid_path)!r})",
                "with open('/dev/null', 'wb') as sink:",
                "    subprocess.Popen(",
                f"        [sys.executable, '-c', {child_source!r}],",
                "        stdin=subprocess.DEVNULL, stdout=sink, stderr=sink)",
                "while not path.exists():",
                "    time.sleep(0.001)",
            ))
            with self.assertRaisesRegex(
                    repair.MeasurementError, "unexpected descendant"):
                self.run_script(leader_source)
            child_pid = int(child_pid_path.read_text())
            self.assert_process_absent(child_pid)

    def test_closed_pipes_before_exit_times_out_and_reaps_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            identity_path = Path(temporary) / "identity"
            source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                "import time",
                f"Path({str(identity_path)!r}).write_text(",
                "    f'{os.getpid()} {os.getpgrp()}')",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "os.close(1)",
                "os.close(2)",
                "time.sleep(60)",
            ))
            started = time.monotonic()
            with self.assertRaisesRegex(
                    repair.MeasurementError,
                    "closed its output pipes before exit"):
                self.run_script(source)
            self.assertLess(time.monotonic() - started, 5.0)
            process_id, process_group = map(
                int, identity_path.read_text().split())
            self.assert_process_absent(process_id)
            self.assert_process_absent(process_group, process_group=True)

    def test_parent_sigkill_kills_leader_and_descendant_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            target_identity_path = Path(temporary) / "target.identity"
            descendant_pid_path = Path(temporary) / "descendant.pid"
            descendant_source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                "import time",
                f"Path({str(descendant_pid_path)!r}).write_text(",
                "    str(os.getpid()))",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                "time.sleep(60)",
            ))
            target_source = "\n".join((
                "import os",
                "from pathlib import Path",
                "import signal",
                "import subprocess",
                "import sys",
                "import time",
                f"Path({str(target_identity_path)!r}).write_text(",
                "    f'{os.getpid()} {os.getpgrp()}')",
                "signal.signal(signal.SIGTERM, signal.SIG_IGN)",
                f"subprocess.Popen([sys.executable, '-c', "
                f"{descendant_source!r}])",
                "time.sleep(60)",
            ))
            helper_source = "\n".join((
                "import sys",
                f"sys.path.insert(0, {HERE!r})",
                "import peel_codec",
                "import repair_timing_codec as repair",
                "repair._run_repairtiming_checked(",
                f"    [{sys.executable!r}, '-c', {target_source!r}],",
                "    peel_codec.isolated_codec_env(),",
                "    stdout_limit=4096, stderr_limit=4096)",
            ))
            helper = subprocess.Popen(
                [sys.executable, "-c", helper_source],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            target_pid = None
            target_group = None
            descendant_pid = None
            try:
                self.wait_for_path(target_identity_path, helper)
                self.wait_for_path(descendant_pid_path, helper)
                target_pid, target_group = map(
                    int, target_identity_path.read_text().split())
                descendant_pid = int(descendant_pid_path.read_text())
                self.assertNotEqual(target_pid, target_group)
                helper.kill()
                helper.wait(timeout=5.0)
                self.assert_process_absent(target_pid)
                self.assert_process_absent(descendant_pid)
                self.assert_process_absent(target_group, process_group=True)
            finally:
                if helper.poll() is None:
                    os.killpg(helper.pid, signal.SIGKILL)
                    helper.wait(timeout=5.0)
                if (
                    target_group is not None and
                    repair._repairtiming_process_group_exists(target_group)
                ):
                    os.killpg(target_group, signal.SIGKILL)

    def test_nonzero_diagnostic_is_stderr_first_and_bounded(self):
        source = "\n".join((
            "import os",
            "os.write(1, b'ignored stdout')",
            "os.write(2, b'failure-detail:' + b'z' * 256)",
            "raise SystemExit(7)",
        ))
        with self.assertRaises(repair.MeasurementError) as caught:
            self.run_script(source, diagnostic_limit=32)
        message = str(caught.exception)
        self.assertIn("exited 7", message)
        self.assertIn("stderr='failure-detail:", message)
        self.assertIn("[truncated]", message)
        self.assertNotIn("ignored stdout", message)
        self.assertLess(len(message), 512)

    def test_concurrent_pipe_drain_and_strict_ascii_decode(self):
        source = "\n".join((
            "import os",
            "for unused in range(16):",
            "    os.write(1, b'o' * 4096)",
            "    os.write(2, b'e' * 4096)",
            "raise SystemExit(3)",
        ))
        with self.assertRaisesRegex(
                repair.MeasurementError, "exited 3"):
            self.run_script(
                source, stdout_limit=65536, stderr_limit=65536)
        with self.assertRaisesRegex(
                repair.MeasurementError, "non-ASCII stdout"):
            self.run_script("import os; os.write(1, b'\\xff')")


class RepairTimingGzipTests(unittest.TestCase):
    def test_exact_member_roundtrip_and_hashes(self):
        payload = b"synthetic repairtiming evidence\n"
        compressed = gzip.compress(payload, mtime=0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            path.write_bytes(compressed)
            text, evidence = repair.read_repairtiming_gzip_text(path)
            exact_text, exact_evidence = repair.read_repairtiming_gzip_text(
                path,
                compressed_limit=len(compressed),
                decompressed_limit=len(payload),
            )
        self.assertEqual(text.encode("ascii"), payload)
        self.assertEqual(exact_text, text)
        self.assertEqual(exact_evidence, evidence)
        self.assertEqual(
            evidence["compressed_sha256"],
            hashlib.sha256(compressed).hexdigest())
        self.assertEqual(
            evidence["decompressed_sha256"],
            hashlib.sha256(payload).hexdigest())

    def test_truncation_trailing_concatenation_and_bombs_rejected(self):
        payload = b"x" * 1024
        compressed = gzip.compress(payload, mtime=0)
        variants = {
            "truncated": compressed[:-1],
            "trailing": compressed + b"x",
            "concatenated": compressed + gzip.compress(b"", mtime=0),
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, value in variants.items():
                with self.subTest(name=name):
                    path = root / f"{name}.gz"
                    path.write_bytes(value)
                    with self.assertRaises(repair.MeasurementError):
                        repair.read_repairtiming_gzip_text(path)
            path = root / "bomb.gz"
            path.write_bytes(compressed)
            with self.assertRaises(repair.MeasurementError):
                repair.read_repairtiming_gzip_text(
                    path, decompressed_limit=len(payload) - 1)
            with self.assertRaises(repair.MeasurementError):
                repair.read_repairtiming_gzip_text(
                    path, compressed_limit=len(compressed) - 1)

    def test_symlink_and_toctou_identity_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target.gz"
            target.write_bytes(gzip.compress(b"x\n", mtime=0))
            link = root / "link.gz"
            link.symlink_to(target)
            with self.assertRaises(repair.MeasurementError):
                repair.read_repairtiming_gzip_text(link)
            actual = os.stat(target, follow_symlinks=False)
            with mock.patch.object(
                    repair.os, "stat",
                    return_value=mock.Mock(
                        **{
                            name: (
                                getattr(actual, name) + 1
                                if name == "st_mtime_ns"
                                else getattr(actual, name)
                            )
                            for name in (
                                "st_dev", "st_ino", "st_mode", "st_size",
                                "st_mtime_ns", "st_ctime_ns",
                            )
                        })):
                with self.assertRaises(repair.MeasurementError):
                    repair.read_repairtiming_gzip_text(target)


if __name__ == "__main__":
    unittest.main()
