#!/usr/bin/env python3
"""Adversarial tests for the independent bandtiming v1 consumer."""

import copy
import gc
import hashlib
import math
import os
import sys
import time
import unittest
from types import SimpleNamespace
from unittest import mock


HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import band_timing_codec as band
import peel_codec


def thermal_row(monotonic, cpu="60.0", dimm="40.0"):
    fields = [
        "2026-07-28T00:00:00Z", f"{monotonic:.9f}",
        "100.0", "4000.0", cpu,
        *(dimm for unused in range(8)),
        "0", "1.0", "1.0", "1.0", "0", "0",
    ]
    return ",".join(fields) + "\n"


def thermal_csv(*rows):
    return (
        ",".join(peel_codec._THERMAL_CSV_COLUMNS) + "\n" +
        "".join(rows)
    ).encode("ascii")


def band_context(start=0.5, end=3.0, *, evict_bytes=4096,
                 cache_state="warm"):
    middle = (start + end) / 2.0
    thermal = peel_codec._thermal_window_from_csv(
        thermal_csv(
            thermal_row(start), thermal_row(middle), thermal_row(end)),
        int(round(start * 1e9)), int(round(end * 1e9)),
        sampling_interval_ms=1000)
    return {
        "schema": peel_codec.PAIRED_CONTEXT_SCHEMA,
        "bound": {
            "cpu_model": "test-cpu",
            "cpu_affinity": [2],
            "cpu_governors": {"2": "performance"},
            "cache_state": cache_state,
            "evict_bytes": evict_bytes,
            "thermal_source": "/tmp/test-band-sampler.csv",
            "thermal_sampling_interval_ms": 1000,
            "clock_domain": peel_codec.PAIRED_CLOCK_DOMAIN,
            "thermal_device": 1,
            "thermal_inode": 2,
            "thermal_prelaunch_monotonic_s": float(start),
            "thermal_prelaunch_row_sha256":
                hashlib.sha256(
                    thermal["prelaunch_row"].encode("ascii")).hexdigest(),
        },
        "thermal": thermal,
    }


def default_candidate(block_count=16):
    unused = block_count
    return band.BandDescriptor(
        staircase=6, dense_rows=4, gf256_rows=8, gf16_rows=2,
        period=97, x_mode="tracking-x")


def parse_kwargs(*, context=None, candidate=None, schedule="iid",
                 systematic_cache="off", cache_state="warm",
                 block_count=16, block_bytes=2, replicates=4,
                 warmups=0, construction_seed=7, loss_seed=9,
                 inner_reps=1, max_overhead=16, required_margin=0.0,
                 loss=0.1):
    return {
        "block_count": block_count,
        "block_bytes": block_bytes,
        "candidate": candidate or default_candidate(block_count),
        "construction_seed": construction_seed,
        "loss": loss,
        "loss_seed": loss_seed,
        "schedule": schedule,
        "warmup_replicates": warmups,
        "replicates": replicates,
        "inner_reps": inner_reps,
        "max_overhead": max_overhead,
        "cache_state": cache_state,
        "systematic_cache": systematic_cache,
        "evict_bytes": 4096,
        "context": context or band_context(cache_state=cache_state),
        "required_margin": required_margin,
    }


def _class_for(code):
    return band._generic_result_class(code)


def _weak_config(weak, replicate, scope, base):
    return (weak or {}).get((replicate, scope, base))


def _fixture_failure_class(config, scope):
    if config is None:
        return None
    construct_result = config.get("construct_result", 3)
    result = config.get(
        "result", construct_result if construct_result != 0 else 4)
    result_class = config["class"]
    if construct_result != 0 or result != 0:
        return result_class
    if scope == "decoder":
        recovery_class = config.get("recovery_class", "success")
        if recovery_class != "success":
            return recovery_class
    return None


def bandtiming_stdout(
        *, candidate=None, schedule="iid", systematic_cache="off",
        cache_state="warm", context=None, block_count=16, block_bytes=2,
        warmups=0, replicates=4, construction_seed=7, loss_seed=9,
        inner_reps=1, max_overhead=16, required_margin=0.0,
        weak=None, mutate_row=None, started_ns=1_000_000_000,
        finished_ns=2_000_000_000, loss=0.1):
    candidate = candidate or default_candidate(block_count)
    dispatch = band.dispatch_band_descriptor(block_count)
    context = context or band_context(cache_state=cache_state)
    total = warmups + replicates
    expected_rows = total * 120
    message_bytes = (block_count - 1) * block_bytes + 1
    manifest = {
        "schema": band.BANDTIMING_SCHEMA,
        "dispatch_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "contract_id": peel_codec.TARGET_CONTRACT["contract_id"],
        "K": block_count,
        "bb": block_bytes,
        "message_bytes": message_bytes,
        "completion": "mixed",
        "candidate_S": candidate.staircase,
        "candidate_D2": candidate.dense_rows,
        "candidate_gf256": candidate.gf256_rows,
        "candidate_gf16": candidate.gf16_rows,
        "candidate_P": candidate.period,
        "candidate_x_geometry": candidate.x_mode,
        "candidate_descriptor_sha256":
            band.descriptor_sha256(candidate, block_count),
        "dispatch_S": dispatch.staircase,
        "dispatch_D2": dispatch.dense_rows,
        "dispatch_gf256": dispatch.gf256_rows,
        "dispatch_gf16": dispatch.gf16_rows,
        "dispatch_P": dispatch.period,
        "dispatch_x_geometry": dispatch.x_mode,
        "dispatch_descriptor_sha256":
            band.descriptor_sha256(dispatch, block_count),
        "descriptor_encoding": "S-D2-gf256-gf16-P-x-newline-v1",
        "recovery_mix_count": 3,
        "construction_seed_base": construction_seed,
        "construction_seed_derivation":
            band.BANDTIMING_CONSTRUCTION_SEED_DERIVATION,
        "wh2_seed_mapping": "zero-extend-u32-precode_xor-fold-packet-v1",
        "wh1_seed_mapping": "low16-peel-high16-dense-v1",
        "wh1_dense_count": peel_codec._legacy_dense_count(block_count),
        "semantic_seed_derivation":
            band.BANDTIMING_SEMANTIC_SEED_DERIVATION,
        "loss": f"{loss:.17g}",
        "loss_seed_base": loss_seed,
        "loss_seed_derivation": band.BANDTIMING_LOSS_SEED_DERIVATION,
        "message_seed_policy": "replicate-loss-seed-partial-final-v1",
        "schedule": schedule,
        "loss_model": "packet-schedule-v1",
        "trace_encoding": "wirehair-wh2-bandtiming-loss-trace-v1",
        "panels": "+".join(band.BANDTIMING_PANELS),
        "scope_order": "+".join(band.BANDTIMING_SCOPE_ORDER),
        "warmup_replicates": warmups,
        "replicates": replicates,
        "slots_per_panel": 8,
        "panels_per_replicate": 15,
        "order": band.BANDTIMING_ORDER,
        "label_swap": band.BANDTIMING_LABEL_SWAP,
        "inner_reps": inner_reps,
        "max_overhead": max_overhead,
        "cache_state": cache_state,
        "systematic_cache": systematic_cache,
        "wh2_source_cache": int(systematic_cache == "on"),
        "wh2_received_systematic_cache": int(systematic_cache == "on"),
        "wh1_encoder_source_policy":
            "borrow-when-off-copy-when-on-v1",
        "wh1_decoder_systematic_policy":
            "native-staging-uncontrolled-v1",
        "wh1_intermediate_policy": "unavailable-zero-v1",
        "evict_bytes": 4096,
        "payload_alignment": 64,
        "prefault": 1,
        "cpu_affinity_policy": "first-allowed-affinity-v1",
        "encoder_scope":
            "fresh-object-init-through-first-K-symbols-v1",
        "decoder_scope":
            "fresh-init-outside-timer-first-feed-through-own-success-v1",
        "direct_scope":
            "candidate-dispatch-pair-local-fixed-prefix-solve-v1",
        "weak_seed_policy": "panel-local-balanced-censor-v1",
        "hook_path":
            "caller-pinned-explicit-transaction-attempt-zero-v2",
        "codec_reuse": "none-fresh-object-every-inner-v1",
        "context_sha256":
            peel_codec._canonical_json_sha256(context["bound"]),
        "uncertainty": band.BANDTIMING_UNCERTAINTY,
        "required_margin": f"{required_margin:.17g}",
        "margin_rule":
            "upper-log-cost-lt-negative-required-margin-and-arm-aa-floors-v1",
        "clock_domain": band.BANDTIMING_CLOCK_DOMAIN,
        "stream_hash_scope": "body-plus-done-prefix-v1",
        "started_monotonic_ns": started_ns,
        "expected_rows": expected_rows,
    }
    digest = hashlib.sha256(b"same").hexdigest()
    message_digest = hashlib.sha256(b"message").hexdigest()
    semantic_bytes = (
        block_count + dispatch.staircase + dispatch.dense_rows +
        dispatch.heavy_rows
    ) * block_bytes
    measured_seeds = {
        band._expected_seed(construction_seed, replicate)
        for replicate in range(total)
    }
    semantic_attempt = next(
        attempt for attempt in range(256)
        if ((construction_seed + attempt) & 0xffffffff) not in measured_seeds)
    semantic_seed = (construction_seed + semantic_attempt) & 0xffffffff
    semantic = {
        "timed": 0,
        "construction_seed": semantic_seed,
        "seed_attempt": semantic_attempt,
        "seed_attempt_cap": 256,
        "canonical_S": dispatch.staircase,
        "canonical_D2": dispatch.dense_rows,
        "canonical_gf256": dispatch.gf256_rows,
        "canonical_gf16": dispatch.gf16_rows,
        "canonical_P": dispatch.period,
        "canonical_x": dispatch.x_mode,
        "trace_sha256": hashlib.sha256(b"semantic-trace").hexdigest(),
        "message_sha256": message_digest,
        "explicit_construct_result": 0,
        "dispatch_construct_result": 0,
        "explicit_params_sha256": digest,
        "dispatch_params_sha256": digest,
        "params_equal": 1,
        "explicit_coefficients_sha256": digest,
        "dispatch_coefficients_sha256": digest,
        "coefficients_equal": 1,
        "explicit_packet_rows_sha256": digest,
        "dispatch_packet_rows_sha256": digest,
        "packet_rows_equal": 1,
        "explicit_intermediate_sha256": digest,
        "dispatch_intermediate_sha256": digest,
        "explicit_intermediate_bytes": semantic_bytes,
        "dispatch_intermediate_bytes": semantic_bytes,
        "intermediate_equal": 1,
        "explicit_payload_sha256": digest,
        "dispatch_payload_sha256": digest,
        "payload_equal": 1,
        "explicit_direct_result": 0,
        "dispatch_direct_result": 0,
        "explicit_solve_sha256": digest,
        "dispatch_solve_sha256": digest,
        "direct_equal": 1,
        "explicit_decode_result": 0,
        "dispatch_decode_result": 0,
        "explicit_overhead": 1,
        "dispatch_overhead": 1,
        "explicit_recovered_sha256": message_digest,
        "dispatch_recovered_sha256": message_digest,
        "recovery_equal": 1,
        "message_equal": 1,
        "pass": 1,
    }
    lines = [
        "# bandtiming," + ",".join(
            f"{name}={manifest[name]}"
            for name in band.BANDTIMING_MANIFEST_FIELDS),
        "# band_semantic," + ",".join(
            f"{name}={semantic[name]}"
            for name in band.BANDTIMING_SEMANTIC_FIELDS),
        ",".join(band.BANDTIMING_COLUMNS),
    ]
    base_overhead = {"candidate": 2, "dispatch": 3, "wh1": 4}
    base_elapsed = {"candidate": 80, "dispatch": 100, "wh1": 120}

    def decoder_receipt(replicate, base):
        config = _weak_config(weak, replicate, "decoder", base)
        if config is None:
            return True, block_count + base_overhead[base]
        construct_result = config.get("construct_result", 3)
        result = config.get(
            "result", construct_result if construct_result != 0 else 4)
        result_class = config["class"]
        recovery_result = (
            config.get("recovery_result", 0)
            if construct_result == 0 and result_class == "success" else -1
        )
        recovery_class = (
            config.get(
                "recovery_class", _class_for(recovery_result))
            if recovery_result != -1 else "not_applicable"
        )
        recovery_ok = config.get(
            "recovery_ok", int(recovery_class == "success"))
        received = config.get(
            "received_symbols", block_count + base_overhead[base])
        return (
            construct_result == 0 and result == 0 and
            result_class == "success" and recovery_result == 0 and
            recovery_class == "success" and recovery_ok == 1,
            received,
        )

    row_index = 0
    for replicate in range(total):
        construction = band._expected_seed(construction_seed, replicate)
        replicate_loss_seed = band._expected_loss_seed(
            loss_seed, replicate)
        trace = hashlib.sha256(
            f"trace-{replicate}".encode("ascii")).hexdigest()
        for panel_index, (scope, panel, left, right) in enumerate(
                band.BANDTIMING_PANEL_SPECS):
            panel_bases = []
            for role in (left, right):
                base = band._base_arm(role)
                if base not in panel_bases:
                    panel_bases.append(base)
            failures = []
            for base in panel_bases:
                config = _weak_config(weak, replicate, scope, base)
                failure_class = _fixture_failure_class(config, scope)
                if failure_class is not None:
                    failures.append(f"{base}_{failure_class}")
            reason = "none" if not failures else "+".join(failures)
            eligible = int(reason == "none")
            direct_prefix = -1
            if scope == "direct":
                decoder_receipts = [
                    decoder_receipt(replicate, base)
                    for base in panel_bases
                ]
                direct_prefix = (
                    max(received for unused_ok, received
                        in decoder_receipts)
                    if all(ok for ok, unused_received
                           in decoder_receipts) else
                    block_count + max_overhead
                )
            for slot, label in enumerate(band.BANDTIMING_ORDER):
                role = band._expected_role(panel, replicate, label)
                base = band._base_arm(role)
                config = _weak_config(weak, replicate, scope, base)
                if config is None:
                    construct_result = 0
                    result = 0
                    result_class = "success"
                else:
                    construct_result = config.get("construct_result", 3)
                    result = config.get(
                        "result",
                        construct_result if construct_result != 0 else 4)
                    result_class = config["class"]
                if (scope == "decoder" and
                        construct_result == 0 and result_class == "success"):
                    recovery_result = (
                        config.get("recovery_result", 0)
                        if config is not None else 0
                    )
                    recovery_class = (
                        config.get(
                            "recovery_class",
                            _class_for(recovery_result))
                        if config is not None else "success"
                    )
                    recovery_ok = (
                        config.get(
                            "recovery_ok",
                            int(recovery_class == "success"))
                        if config is not None else 1
                    )
                else:
                    recovery_result = -1
                    recovery_class = "not_applicable"
                    recovery_ok = 0
                zero_prefix = construct_result != 0
                if zero_prefix:
                    encoded = received = 0
                    overhead = fixed = -1
                    packet_bytes = 0
                elif scope == "encoder":
                    encoded = block_count
                    received = 0
                    overhead = 0
                    fixed = -1
                    packet_bytes = message_bytes
                elif scope == "decoder":
                    encoded = 0
                    received = (
                        config.get(
                            "received_symbols",
                            block_count + base_overhead[base])
                        if config is not None else
                        block_count + base_overhead[base]
                    )
                    overhead = (
                        received - block_count if recovery_ok else -1)
                    fixed = -1
                    packet_bytes = (
                        received * block_bytes
                        if schedule in ("repair-only", "adversarial") else
                        max(
                            0,
                            received * block_bytes - (block_bytes - 1))
                    )
                else:
                    encoded = 0
                    received = (
                        config.get(
                            "received_symbols", direct_prefix)
                        if config is not None else direct_prefix
                    )
                    overhead = received - block_count
                    fixed = received
                    packet_bytes = (
                        received * block_bytes
                        if schedule in ("repair-only", "adversarial") else
                        max(
                            0,
                            received * block_bytes - (block_bytes - 1))
                    )
                if config is not None and "packet_payload_bytes" in config:
                    packet_bytes = config["packet_payload_bytes"]
                stats_available = int(
                    eligible and scope in ("decoder", "direct") and
                    base != "wh1")
                descriptor = candidate if base == "candidate" else dispatch
                if stats_available:
                    stats = {
                        "inactivated": descriptor.heavy_rows + 3,
                        "binary_def": descriptor.heavy_rows,
                        "heavy_gain": descriptor.heavy_rows,
                        "block_xors": 20,
                        "block_muladds": 10,
                        "build_ns_sum": 5 + slot,
                        "peel_ns_sum": 6 + slot,
                        "project_ns_sum": 7 + slot,
                        "residual_ns_sum": 8 + slot,
                        "backsub_ns_sum": 9 + slot,
                    }
                else:
                    stats = {
                        name: -1 for name in band._SOLVE_STATS_FIELDS
                    }
                intermediate = (
                    0 if base == "wh1" else
                    (
                        block_count + descriptor.staircase +
                        descriptor.dense_rows + descriptor.heavy_rows
                    ) * block_bytes
                )
                elapsed = base_elapsed[base] if eligible else 0
                row = {
                    "replicate": replicate,
                    "measured": int(replicate >= warmups),
                    "scope": scope,
                    "panel": panel,
                    "panel_index": panel_index,
                    "slot": slot,
                    "pair": slot // 2,
                    "label": label,
                    "role": role,
                    "label_swap": replicate & 1,
                    "construction_seed": construction,
                    "loss_seed": replicate_loss_seed,
                    "trace_sha256": trace,
                    "construct_result": construct_result,
                    "result": result,
                    "result_class": result_class,
                    "recovery_result": recovery_result,
                    "recovery_class": recovery_class,
                    "recovery_ok": recovery_ok,
                    "encoded_symbols": encoded,
                    "received_symbols": received,
                    "arm_overhead": overhead,
                    "fixed_prefix_symbols": fixed,
                    "timing_eligible": eligible,
                    "panel_censored": 1 - eligible,
                    "censor_reason": reason,
                    "preflight_result": result,
                    "timing_result": result if eligible else -1,
                    "outcome_stable": eligible,
                    "elapsed_ns": elapsed,
                    "inner_reps": inner_reps if eligible else 0,
                    "saturated": 0,
                    "cpu_before": 2 if eligible else -1,
                    "cpu_after": 2 if eligible else -1,
                    "cpu_migrated": 0,
                    "minflt_delta": 0,
                    "majflt_delta": 0,
                    "fault_contaminated": 0,
                    "stats_available": stats_available,
                    **stats,
                    "source_bytes": message_bytes,
                    "packet_payload_bytes": packet_bytes,
                    "intermediate_bytes": intermediate,
                }
                if mutate_row is not None:
                    mutate_row(row_index, row)
                lines.append(",".join(
                    str(row[name]) for name in band.BANDTIMING_COLUMNS))
                row_index += 1
    prefix = "\n".join(lines) + "\n"
    done_prefix = (
        f"# bandtiming_done,complete=1,rows={expected_rows},"
        f"finished_monotonic_ns={finished_ns},stream_sha256="
    )
    stream_sha = hashlib.sha256(
        (prefix + done_prefix).encode("ascii")).hexdigest()
    return prefix + done_prefix + stream_sha + "\n"


def authenticated_edit(stdout, edit):
    lines = stdout.splitlines()
    edit(lines)
    done_values = dict(
        field.split("=", 1)
        for field in lines[-1][len("# bandtiming_done,"):].split(","))
    prefix = "\n".join(lines[:-1]) + "\n"
    done_prefix = (
        f"# bandtiming_done,complete={done_values['complete']},"
        f"rows={done_values['rows']},"
        f"finished_monotonic_ns={done_values['finished_monotonic_ns']},"
        "stream_sha256="
    )
    lines[-1] = done_prefix + hashlib.sha256(
        (prefix + done_prefix).encode("ascii")).hexdigest()
    return "\n".join(lines) + "\n"


def edit_row(stdout, row_index, **updates):
    def edit(lines):
        values = lines[3 + row_index].split(",")
        row = dict(zip(band.BANDTIMING_COLUMNS, values))
        row.update({name: str(value) for name, value in updates.items()})
        lines[3 + row_index] = ",".join(
            row[name] for name in band.BANDTIMING_COLUMNS)
    return authenticated_edit(stdout, edit)


class BandTimingParserTests(unittest.TestCase):

    def test_golden_stream_has_exact_480_rows_and_independent_summaries(self):
        context = band_context()
        kwargs = parse_kwargs(context=context)
        stdout = bandtiming_stdout(context=context)
        result = band.parse_bandtiming_output(stdout, **kwargs)
        self.assertEqual(len(stdout.splitlines()), 484)
        self.assertEqual(len(result.arms), 8)
        self.assertEqual(len(result.contrasts), 7)
        self.assertEqual(len(result.replicate_receipts), 4)
        self.assertEqual(result.weak_cells, ())
        self.assertEqual(
            {
                name: result.replicate_receipts[0][name]
                for name in (
                    "candidate_construct_result",
                    "candidate_construct_class",
                    "dispatch_construct_result",
                    "dispatch_construct_class",
                    "wh1_construct_result",
                    "wh1_construct_class",
                )
            },
            {
                "candidate_construct_result": 0,
                "candidate_construct_class": "success",
                "dispatch_construct_result": 0,
                "dispatch_construct_class": "success",
                "wh1_construct_result": 0,
                "wh1_construct_class": "success",
            },
        )
        self.assertEqual(
            [contrast.name for contrast in result.contrasts],
            list(band.BANDTIMING_CROSS_PANELS))
        self.assertEqual(
            [(arm.scope, arm.arm) for arm in result.arms],
            list(band._ARM_SCOPE_ORDER))
        self.assertTrue(all(arm.trials == 4 for arm in result.arms))
        self.assertTrue(all(arm.successes == 4 for arm in result.arms))
        self.assertTrue(all(arm.aa_eligible_replicates == 4
                            for arm in result.arms))
        self.assertTrue(all(arm.aa_floor_log == 0.0 for arm in result.arms))
        self.assertEqual(
            result.semantic["explicit_intermediate_bytes"],
            (
                kwargs["block_count"] +
                result.dispatch_descriptor.staircase +
                result.dispatch_descriptor.dense_rows +
                result.dispatch_descriptor.heavy_rows
            ) * kwargs["block_bytes"])
        context["bound"]["cpu_affinity"].append(99)
        self.assertEqual(
            result.context["bound"]["cpu_affinity"], [2])

    def test_parser_scaling_is_linear_through_512_replicates(self):
        # Generate outside the timed region, use process CPU time to exclude
        # scheduler stalls, and suppress nondeterministic cyclic-GC pauses.
        # The 2x-linear endpoint allowance is intentionally loose for slow
        # CI hosts while still separating the former quadratic path, whose
        # measured 512:32 ratio exceeded 140x.
        timings = []
        gc_was_enabled = gc.isenabled()
        for replicates in (32, 64, 128, 256, 512):
            stdout = bandtiming_stdout(replicates=replicates)
            gc.collect()
            if gc_was_enabled:
                gc.disable()
            try:
                started = time.process_time()
                parsed = band.parse_bandtiming_output(
                    stdout, **parse_kwargs(replicates=replicates))
                elapsed = time.process_time() - started
            finally:
                if gc_was_enabled:
                    gc.enable()
            self.assertEqual(len(parsed.replicate_receipts), replicates)
            timings.append((replicates, elapsed))
            del parsed, stdout
            gc.collect()

        endpoint_ratio = timings[-1][1] / timings[0][1]
        normalized = [
            elapsed / replicates for replicates, elapsed in timings
        ]
        details = ", ".join(
            f"R={replicates}:{elapsed:.6f}s"
            for replicates, elapsed in timings)
        self.assertLessEqual(
            endpoint_ratio, 32.0,
            f"parser endpoint scaling is superlinear: {details}")
        self.assertLessEqual(
            max(normalized) / min(normalized), 2.5,
            f"parser per-replicate cost is unstable: {details}")

    def test_every_schedule_label_swap_and_cache_regime_round_trip(self):
        for schedule in sorted(band.BANDTIMING_SCHEDULES):
            for systematic_cache in ("off", "on"):
                with self.subTest(
                        schedule=schedule,
                        systematic_cache=systematic_cache):
                    context = band_context()
                    result = band.parse_bandtiming_output(
                        bandtiming_stdout(
                            schedule=schedule,
                            systematic_cache=systematic_cache,
                            context=context),
                        **parse_kwargs(
                            schedule=schedule,
                            systematic_cache=systematic_cache,
                            context=context))
                    self.assertEqual(result.manifest["schedule"], schedule)
                    self.assertEqual(
                        result.manifest["systematic_cache"],
                        systematic_cache)
                    self.assertEqual(
                        result.manifest["wh2_source_cache"],
                        int(systematic_cache == "on"))
        repair_stdout = bandtiming_stdout(
            schedule="repair-only", block_bytes=64)
        with self.assertRaisesRegex(
                band.MeasurementError, "payload byte count"):
            band.parse_bandtiming_output(
                edit_row(
                    repair_stdout, 6 * 8,
                    packet_payload_bytes=(16 + 2) * 64 - 63),
                **parse_kwargs(schedule="repair-only", block_bytes=64))
        repair_direct_failure = {
            (0, "direct", "candidate"): {
                "construct_result": 0, "result": 8, "class": "error",
            },
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "payload byte count"):
            band.parse_bandtiming_output(
                edit_row(
                    bandtiming_stdout(
                        schedule="repair-only", block_bytes=64,
                        weak=repair_direct_failure),
                    12 * 8, packet_payload_bytes=1),
                **parse_kwargs(schedule="repair-only", block_bytes=64))
        def wrong_odd_label(unused_index, row):
            if (row["replicate"] == 1 and row["panel_index"] == 0 and
                    row["slot"] == 0):
                row["role"] = "candidate"
        with self.assertRaisesRegex(
                band.MeasurementError, "reordered or mislabeled"):
            band.parse_bandtiming_output(
                bandtiming_stdout(mutate_row=wrong_odd_label),
                **parse_kwargs())

    def test_warmup_outcome_cannot_change_across_panels(self):
        def forge_one_warmup_panel(lines):
            for offset in range(8):
                index = 3 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if band._base_arm(row["role"]) == "candidate":
                    row.update({
                        "result": "8",
                        "result_class": "error",
                        "preflight_result": "8",
                    })
                row.update({
                    "timing_eligible": "0",
                    "panel_censored": "1",
                    "censor_reason": "candidate_error",
                    "timing_result": "-1",
                    "outcome_stable": "0",
                    "elapsed_ns": "0",
                    "inner_reps": "0",
                    "cpu_before": "-1",
                    "cpu_after": "-1",
                })
                lines[index] = ",".join(
                    row[name] for name in band.BANDTIMING_COLUMNS)

        stdout = authenticated_edit(
            bandtiming_stdout(warmups=1, replicates=3),
            forge_one_warmup_panel)
        with self.assertRaisesRegex(
                band.MeasurementError, "outcome changed across panels"):
            band.parse_bandtiming_output(
                stdout, **parse_kwargs(warmups=1, replicates=3))

    def test_three_seed_screen_keeps_confidence_intervals_unavailable(self):
        result = band.parse_bandtiming_output(
            bandtiming_stdout(replicates=3),
            **parse_kwargs(replicates=3))
        self.assertTrue(all(
            contrast.eligible_replicates == 3 and
            not contrast.timing_ci_available
            for contrast in result.contrasts))
        contrast = next(
            item for item in result.contrasts
            if item.name == "encoder_candidate_dispatch")
        self.assertAlmostEqual(
            contrast.log_cost_mean, math.log(0.8), places=14)
        self.assertIsNone(contrast.log_cost_ci_low)
        self.assertIsNone(contrast.log_cost_ci_high)
        self.assertTrue(all(
            arm.aa_eligible_replicates == 3 and
            arm.aa_log_cost_mean == 0.0 and
            arm.aa_log_cost_ci_low is None and
            arm.aa_log_cost_ci_high is None and
            arm.aa_floor_log == 0.0
            for arm in result.arms))

    def test_fault_contamination_excludes_the_whole_paired_panel(self):
        stdout = edit_row(
            bandtiming_stdout(), 0,
            minflt_delta=1, fault_contaminated=1)
        result = band.parse_bandtiming_output(stdout, **parse_kwargs())
        contrast = next(
            item for item in result.contrasts
            if item.name == "encoder_candidate_dispatch")
        self.assertEqual(contrast.eligible_replicates, 3)
        candidate = next(
            arm for arm in result.arms
            if (arm.scope, arm.arm) == ("encoder", "candidate"))
        dispatch = next(
            arm for arm in result.arms
            if (arm.scope, arm.arm) == ("encoder", "dispatch"))
        self.assertEqual(candidate.timed_slots, 60)
        self.assertEqual(dispatch.timed_slots, 60)
        self.assertEqual(
            result.replicate_receipts[0]["fault_contaminated_panels"],
            ["encoder_candidate_dispatch"])
        # A cross-panel CI alone cannot claim a winner when either arm has
        # fewer than four clean replicates in its own A/A control.
        stdout = edit_row(
            bandtiming_stdout(), 3 * 8,
            minflt_delta=1, fault_contaminated=1)
        result = band.parse_bandtiming_output(stdout, **parse_kwargs())
        contrast = next(
            item for item in result.contrasts
            if item.name == "encoder_candidate_dispatch")
        self.assertTrue(contrast.timing_ci_available)
        self.assertFalse(contrast.left_faster)
        candidate = next(
            arm for arm in result.arms
            if (arm.scope, arm.arm) == ("encoder", "candidate"))
        self.assertEqual(candidate.aa_eligible_replicates, 3)

    def test_candidate_weak_construction_censors_only_its_panels(self):
        weak = {
            (0, scope, "candidate"): {
                "construct_result": 3, "class": "weak",
            }
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        result = band.parse_bandtiming_output(
            bandtiming_stdout(weak=weak), **parse_kwargs())
        self.assertEqual(
            [(item["scope"], item["category"])
             for item in result.weak_cells],
            [
                ("encoder", "candidate-only"),
                ("decoder", "candidate-only"),
                ("direct", "candidate-only"),
            ])
        first = result.replicate_receipts[0]
        self.assertEqual(first["candidate_construct_result"], 3)
        self.assertEqual(first["candidate_construct_class"], "weak")
        self.assertIn("encoder_candidate_dispatch",
                      first["censored_panels"])
        self.assertIn("encoder_candidate_wh1", first["censored_panels"])
        self.assertIn("encoder_candidate_aa", first["censored_panels"])
        self.assertNotIn("encoder_dispatch_wh1",
                         first["censored_panels"])
        self.assertNotIn("encoder_dispatch_aa",
                         first["censored_panels"])
        dispatch = next(
            arm for arm in result.arms
            if arm.scope == "encoder" and arm.arm == "dispatch")
        self.assertGreater(dispatch.timed_slots, 0)

    def test_three_seed_probe_certified_singularity_is_normalized_weak(self):
        candidate = band.BandDescriptor(
            staircase=1, dense_rows=1, gf256_rows=1, gf16_rows=0,
            period=1, x_mode="tracking-x")
        weak = {
            (replicate, scope, "candidate"): {
                "construct_result": 4,
                "class": "weak",
            }
            for replicate in range(3)
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        result = band.parse_bandtiming_output(
            bandtiming_stdout(
                block_count=10, candidate=candidate,
                replicates=3, weak=weak),
            **parse_kwargs(
                block_count=10, candidate=candidate, replicates=3))
        self.assertEqual(
            [(item["replicate"], item["scope"], item["category"],
              item["weak_arms"])
             for item in result.weak_cells],
            [
                (replicate, scope, "candidate-only", ["candidate"])
                for replicate in range(3)
                for scope in band.BANDTIMING_SCOPE_ORDER
            ])
        for receipt in result.replicate_receipts:
            for scope in band.BANDTIMING_SCOPE_ORDER:
                self.assertEqual(
                    [panel for panel in receipt["censored_panels"]
                     if panel.startswith(f"{scope}_") and
                     "candidate" in panel],
                    [
                        f"{scope}_candidate_dispatch"
                        if scope != "direct" else
                        "direct_candidate_dispatch",
                        f"{scope}_candidate_wh1",
                        f"{scope}_candidate_aa",
                    ] if scope != "direct" else [
                        "direct_candidate_dispatch",
                        "direct_candidate_aa",
                    ])
                self.assertNotIn(
                    f"{scope}_dispatch_wh1",
                    receipt["censored_panels"])

    def test_postconstruction_and_recovery_failures_keep_native_counters(self):
        failures = {
            (0, "encoder", "candidate"): {
                "construct_result": 0, "result": 8, "class": "error",
            },
            (0, "decoder", "candidate"): {
                "construct_result": 0, "result": 1,
                "class": "need_more", "received_symbols": 32,
            },
            (1, "decoder", "candidate"): {
                "construct_result": 0, "result": 0, "class": "success",
                "recovery_result": 8, "recovery_class": "error",
                "recovery_ok": 0,
            },
            (0, "direct", "candidate"): {
                "construct_result": 0, "result": 1,
                "class": "need_more",
            },
        }
        result = band.parse_bandtiming_output(
            bandtiming_stdout(weak=failures), **parse_kwargs())
        candidate = {
            arm.scope: arm for arm in result.arms
            if arm.arm == "candidate"
        }
        self.assertEqual(candidate["encoder"].errors, 1)
        self.assertEqual(candidate["decoder"].need_more, 1)
        self.assertEqual(candidate["decoder"].errors, 1)
        self.assertEqual(candidate["direct"].need_more, 1)
        self.assertEqual(result.weak_cells, ())
        first = result.replicate_receipts[0]
        self.assertEqual(first["encoder_candidate_overhead"], 0)
        self.assertEqual(first["decoder_candidate_overhead"], -1)
        self.assertEqual(first["direct_candidate_overhead"], 16)
        self.assertEqual(first["candidate_construct_result"], 0)
        self.assertEqual(first["candidate_construct_class"], "success")

        stdout = bandtiming_stdout()
        with self.assertRaisesRegex(
                band.MeasurementError, "result class contradicts"):
            band.parse_bandtiming_output(
                edit_row(stdout, 6 * 8, recovery_ok=2),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "output witness"):
            band.parse_bandtiming_output(
                edit_row(
                    stdout, 6 * 8, recovery_result=-1,
                    recovery_class="not_applicable", recovery_ok=0),
                **parse_kwargs())

    def test_all_scope_payload_error_is_not_construction_failure(self):
        failures = {
            (0, scope, "candidate"): {
                "construct_result": 0,
                "result": 8,
                "class": "error",
            }
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        result = band.parse_bandtiming_output(
            bandtiming_stdout(weak=failures), **parse_kwargs())
        first = result.replicate_receipts[0]
        self.assertEqual(first["candidate_construct_result"], 0)
        self.assertEqual(first["candidate_construct_class"], "success")
        self.assertEqual(
            [
                first[f"{scope}_candidate_result_class"]
                for scope in band.BANDTIMING_SCOPE_ORDER
            ],
            ["error", "error", "error"],
        )

    def test_hash_order_cardinality_trace_and_semantic_tampering_fail(self):
        stdout = bandtiming_stdout()
        with self.assertRaisesRegex(band.MeasurementError, "SHA-256"):
            band.parse_bandtiming_output(
                stdout[:-2] + ("0" if stdout[-2] != "0" else "1") + "\n",
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "reordered or mislabeled"):
            def swap_rows(lines):
                lines[3], lines[4] = lines[4], lines[3]
            band.parse_bandtiming_output(
                authenticated_edit(stdout, swap_rows), **parse_kwargs())
        with self.assertRaisesRegex(band.MeasurementError, "emitted"):
            def drop_row(lines):
                del lines[3]
            band.parse_bandtiming_output(
                authenticated_edit(stdout, drop_row), **parse_kwargs())
        changed_trace = edit_row(
            stdout, 1,
            trace_sha256=hashlib.sha256(b"foreign").hexdigest())
        with self.assertRaisesRegex(band.MeasurementError, "changed trace"):
            band.parse_bandtiming_output(changed_trace, **parse_kwargs())
        def fail_semantic(lines):
            lines[1] = lines[1].replace("params_equal=1", "params_equal=0")
        with self.assertRaisesRegex(
                band.MeasurementError, "semantic bridge"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, fail_semantic), **parse_kwargs())
        def exceed_semantic_overhead(lines):
            lines[1] = lines[1].replace(
                "explicit_overhead=1,dispatch_overhead=1",
                "explicit_overhead=17,dispatch_overhead=17")
        with self.assertRaisesRegex(
                band.MeasurementError, "semantic bridge"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, exceed_semantic_overhead),
                **parse_kwargs())
        def forge_descriptor(lines):
            lines[0] = lines[0].replace(
                "candidate_descriptor_sha256=" +
                band.descriptor_sha256(default_candidate(), 16),
                "candidate_descriptor_sha256=" + "0" * 64)
        with self.assertRaisesRegex(
                band.MeasurementError, "candidate_descriptor_sha256"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, forge_descriptor),
                **parse_kwargs())

    def test_scope_counters_prefixes_and_byte_receipts_are_honest(self):
        stdout = bandtiming_stdout()
        cases = (
            (
                {"stats_available": 1, "inactivated": 0},
                "stats availability",
            ),
            (
                {"encoded_symbols": 15},
                "scope-specific symbol",
            ),
            (
                {"source_bytes": 32},
                "byte provenance",
            ),
        )
        # Row zero is an encoder candidate: every solver counter is sentinel.
        for updates, message in cases:
            with self.subTest(updates=updates):
                with self.assertRaisesRegex(band.MeasurementError, message):
                    band.parse_bandtiming_output(
                        edit_row(stdout, 0, **updates), **parse_kwargs())
        # First WH1 decoder row: panel 7, and its solver stats are unavailable.
        wh1_decoder = 7 * 8 + 1
        with self.assertRaisesRegex(
                band.MeasurementError, "stats availability"):
            band.parse_bandtiming_output(
                edit_row(
                    stdout, wh1_decoder,
                    stats_available=1, inactivated=0),
                **parse_kwargs())
        # First direct panel is index 12.  Preserve candidate stability while
        # making it disagree with dispatch inside every pair.
        with self.assertRaisesRegex(
                band.MeasurementError, "fixed prefix"):
            def change_candidate_prefix(lines):
                for offset in range(8):
                    line_index = 3 + 12 * 8 + offset
                    values = lines[line_index].split(",")
                    row = dict(zip(band.BANDTIMING_COLUMNS, values))
                    if row["role"] == "candidate":
                        row.update({
                            "received_symbols": "20",
                            "arm_overhead": "4",
                            "fixed_prefix_symbols": "20",
                            "packet_payload_bytes": "39",
                        })
                        lines[line_index] = ",".join(
                            row[name] for name in band.BANDTIMING_COLUMNS)
            band.parse_bandtiming_output(
                authenticated_edit(stdout, change_candidate_prefix),
                **parse_kwargs())
        # Decoder never advertises a common fixed prefix.
        with self.assertRaisesRegex(
                band.MeasurementError, "scope-specific symbol"):
            band.parse_bandtiming_output(
                edit_row(stdout, 6 * 8, fixed_prefix_symbols=18),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "contaminated"):
            band.parse_bandtiming_output(
                edit_row(stdout, 0, minflt_delta=-1),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "contaminated"):
            band.parse_bandtiming_output(
                edit_row(stdout, 0, cpu_before=99, cpu_after=99),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "solver statistics"):
            band.parse_bandtiming_output(
                edit_row(stdout, 6 * 8, build_ns_sum=60),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "solver statistics"):
            band.parse_bandtiming_output(
                edit_row(stdout, 6 * 8, inactivated=999999),
                **parse_kwargs())
        with self.assertRaisesRegex(
                band.MeasurementError, "scope-specific symbol"):
            band.parse_bandtiming_output(
                edit_row(
                    stdout, 6 * 8, received_symbols=0,
                    arm_overhead=-16, packet_payload_bytes=0),
                **parse_kwargs())
        # First-success is reached before recovery.  A later recovery failure
        # cannot make a fresh decoder's authenticated receive prefix disappear.
        with self.assertRaisesRegex(
                band.MeasurementError, "scope-specific symbol"):
            band.parse_bandtiming_output(
                edit_row(
                    stdout, 6 * 8,
                    recovery_result=2, recovery_class="error", recovery_ok=0,
                    received_symbols=0, arm_overhead=-1,
                    packet_payload_bytes=0),
                **parse_kwargs())
        # A complete prefix contains full blocks and, at most, the one-byte
        # final systematic packet.  Intermediate byte counts are impossible.
        bb64_stdout = bandtiming_stdout(block_bytes=64)
        with self.assertRaisesRegex(
                band.MeasurementError, "payload byte count"):
            band.parse_bandtiming_output(
                edit_row(
                    bb64_stdout, 6 * 8,
                    packet_payload_bytes=(16 + 2) * 64 - 1),
                **parse_kwargs(block_bytes=64))
        # Native preflights one decoder outcome per arm and replicate, so a
        # secondary panel cannot advertise a different first-success prefix.
        def change_secondary_decoder_prefix(lines):
            for offset in range(8):
                index = 3 + 7 * 8 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if band._base_arm(row["role"]) == "candidate":
                    row.update({
                        "received_symbols": "19",
                        "arm_overhead": "3",
                        "packet_payload_bytes": "37",
                    })
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "outcome changed across panels"):
            band.parse_bandtiming_output(
                authenticated_edit(
                    stdout, change_secondary_decoder_prefix),
                **parse_kwargs())

        def change_direct_prefix(lines):
            for offset in range(8):
                index = 3 + 12 * 8 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                row.update({
                    "received_symbols": "16",
                    "arm_overhead": "0",
                    "fixed_prefix_symbols": "16",
                    "packet_payload_bytes": "31",
                })
                lines[index] = ",".join(
                    row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "decoder-bound fixed prefix"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, change_direct_prefix),
                **parse_kwargs())

        def change_direct_candidate_payload(lines):
            for offset in range(8):
                index = 3 + 12 * 8 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if band._base_arm(row["role"]) == "candidate":
                    row["packet_payload_bytes"] = str(
                        int(row["packet_payload_bytes"]) + 1)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "decoder-bound fixed prefix"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, change_direct_candidate_payload),
                **parse_kwargs())

        bb64_stdout = bandtiming_stdout(block_bytes=64)
        with self.assertRaisesRegex(
                band.MeasurementError, "payload byte count"):
            band.parse_bandtiming_output(
                edit_row(
                    bb64_stdout, 12 * 8,
                    result=8, result_class="error", preflight_result=8,
                    packet_payload_bytes=2),
                **parse_kwargs(block_bytes=64))

        early_payload_failures = {
            (0, "decoder", "candidate"): {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            },
            (0, "direct", "candidate"): {
                "construct_result": 0, "result": 1, "class": "need_more",
            },
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "shared payload failure"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=early_payload_failures),
                **parse_kwargs())

        stopped_need_more = {
            (0, "decoder", "candidate"): {
                "construct_result": 0, "result": 1,
                "class": "need_more", "received_symbols": 16,
            },
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "stopped early"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=stopped_need_more),
                **parse_kwargs())

        matching_early_payload_failure = {
            (0, "decoder", "candidate"): {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            },
            (0, "direct", "candidate"): {
                "construct_result": 0, "result": 8, "class": "error",
                "packet_payload_bytes": 0,
            },
        }

        def claim_full_failed_payload(lines):
            for index in range(3, 3 + 120):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (row["scope"] == "direct" and
                        band._base_arm(row["role"]) == "candidate"):
                    row["packet_payload_bytes"] = str(
                        int(row["received_symbols"]) * 2)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "complete direct prefix"):
            band.parse_bandtiming_output(
                authenticated_edit(
                    bandtiming_stdout(
                        weak=matching_early_payload_failure),
                    claim_full_failed_payload),
                **parse_kwargs())

        impossible_partial_history = {}
        for arm, packet_bytes in (("candidate", 9), ("dispatch", 12)):
            impossible_partial_history[(0, "decoder", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            }
            impossible_partial_history[(0, "direct", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "packet_payload_bytes": packet_bytes,
            }
        with self.assertRaisesRegex(
                band.MeasurementError, "packet-byte history"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=impossible_partial_history),
                **parse_kwargs())

        lossless_missing_partial_tail = {}
        for arm in ("candidate", "dispatch"):
            lossless_missing_partial_tail[(0, "decoder", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            }
            lossless_missing_partial_tail[(0, "direct", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "packet_payload_bytes": 16 * 2,
            }
        with self.assertRaisesRegex(
                band.MeasurementError, "deterministic systematic tail"):
            band.parse_bandtiming_output(
                bandtiming_stdout(
                    loss=0.0, weak=lossless_missing_partial_tail),
                **parse_kwargs(loss=0.0))

        lossless_premature_partial_tail = {
            (0, "decoder", arm): {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            }
            for arm in ("candidate", "dispatch", "wh1")
        }
        for arm in ("candidate", "dispatch"):
            lossless_premature_partial_tail[(0, "direct", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "packet_payload_bytes": 9,
            }
        with self.assertRaisesRegex(
                band.MeasurementError, "deterministic systematic tail"):
            band.parse_bandtiming_output(
                bandtiming_stdout(
                    loss=0.0, weak=lossless_premature_partial_tail),
                **parse_kwargs(loss=0.0))

        lossless_permutation_missing_tail = {
            (0, "decoder", arm): {
                "construct_result": 0, "result": 8, "class": "error",
                "received_symbols": 0,
            }
            for arm in ("candidate", "dispatch", "wh1")
        }
        for arm in ("candidate", "dispatch"):
            lossless_permutation_missing_tail[(0, "direct", arm)] = {
                "construct_result": 0, "result": 8, "class": "error",
                "packet_payload_bytes": (16 + 512) * 2,
            }
        with self.assertRaisesRegex(
                band.MeasurementError, "deterministic systematic tail"):
            band.parse_bandtiming_output(
                bandtiming_stdout(
                    loss=0.0, schedule="permutation", max_overhead=528,
                    weak=lossless_permutation_missing_tail),
                **parse_kwargs(
                    loss=0.0, schedule="permutation", max_overhead=528))

        encoder_only_constructor_failure = {
            (0, "encoder", "candidate"): {
                "construct_result": 4, "class": "weak",
            },
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "constructor code changed"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=encoder_only_constructor_failure),
                **parse_kwargs())

        def make_dispatch_decoder_prefix_full(lines):
            for index in range(3, 3 + 120):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (row["scope"] == "decoder" and
                        band._base_arm(row["role"]) == "dispatch"):
                    row["packet_payload_bytes"] = str(
                        int(row["received_symbols"]) * 64)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "shared payload prefix"):
            band.parse_bandtiming_output(
                authenticated_edit(
                    bb64_stdout, make_dispatch_decoder_prefix_full),
                **parse_kwargs(block_bytes=64))

        def make_every_prefix_19_full(lines):
            make_dispatch_decoder_prefix_full(lines)
            for index in range(3, 3 + 120):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (row["scope"] == "direct" and
                        row["received_symbols"] == "19"):
                    row["packet_payload_bytes"] = str(19 * 64)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "shared payload prefix"):
            band.parse_bandtiming_output(
                authenticated_edit(bb64_stdout, make_every_prefix_19_full),
                **parse_kwargs(block_bytes=64))

        def deliver_tail_after_its_candidate_deadline(lines):
            for index in range(3, 3 + 120):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (band._base_arm(row["role"]) == "candidate" and
                        row["scope"] in ("decoder", "direct") and
                        row["received_symbols"] == "18"):
                    row["packet_payload_bytes"] = str(18 * 64)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "packet-byte history"):
            band.parse_bandtiming_output(
                authenticated_edit(
                    bb64_stdout, deliver_tail_after_its_candidate_deadline),
                **parse_kwargs(block_bytes=64))

        def make_every_complete_prefix_full(lines):
            for index in range(3, len(lines) - 1):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (row["scope"] in ("decoder", "direct") and
                        row["received_symbols"] != "0"):
                    row["packet_payload_bytes"] = str(
                        int(row["received_symbols"]) * 64)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "payload byte count"):
            band.parse_bandtiming_output(
                authenticated_edit(
                    bandtiming_stdout(
                        block_bytes=64, loss=0.0, schedule="iid"),
                    make_every_complete_prefix_full),
                **parse_kwargs(
                    block_bytes=64, loss=0.0, schedule="iid"))

        def change_decoder_counts(lines):
            for offset in range(8):
                index = 3 + 7 * 8 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if (band._base_arm(row["role"]) == "candidate" and
                        row["stats_available"] == "1"):
                    row["inactivated"] = str(int(row["inactivated"]) + 1)
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "solver counts changed"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, change_decoder_counts),
                **parse_kwargs())

        def change_direct_dispatch_counts(lines):
            for offset in range(8):
                index = 3 + 14 * 8 + offset
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                row["inactivated"] = str(int(row["inactivated"]) + 1)
                lines[index] = ",".join(
                    row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "direct/dispatch solver counts"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, change_direct_dispatch_counts),
                **parse_kwargs())

        with self.assertRaisesRegex(band.MeasurementError, "reports OOM"):
            band.parse_bandtiming_output(
                edit_row(
                    stdout, 0, construct_result=9, result=9,
                    result_class="error", preflight_result=9,
                    timing_result=9),
                **parse_kwargs())

    def test_result_codes_cannot_be_relabelled_as_weak(self):
        stdout = bandtiming_stdout()
        cases = (
            (
                {"result": 2, "result_class": "weak"},
                "class contradicts",
            ),
            (
                {"result": 11, "result_class": "error"},
                "class contradicts",
            ),
        )
        for updates, message in cases:
            with self.subTest(updates=updates):
                with self.assertRaisesRegex(band.MeasurementError, message):
                    band.parse_bandtiming_output(
                        edit_row(stdout, 0, **updates), **parse_kwargs())
        oom_as_weak = {
            (0, "encoder", "candidate"): {
                "construct_result": 2, "class": "weak",
            },
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "construction failure class"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=oom_as_weak), **parse_kwargs())
        # A probe-certified raw singularity is normalized to an explicit
        # weak-seed code before it reaches the parser.
        weak = {
            (0, scope, "candidate"): {
                "construct_result": 4,
                "class": "weak",
            }
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        parsed = band.parse_bandtiming_output(
            bandtiming_stdout(weak=weak), **parse_kwargs())
        self.assertEqual(parsed.weak_cells[0]["category"], "candidate-only")
        with self.assertRaisesRegex(
                band.MeasurementError, "construction failure class"):
            # A genuine invariant Error cannot be textually relabeled weak,
            # including in the explicit encoder scope.
            encoder_error_as_weak = {
                (0, "encoder", "candidate"): {
                    "construct_result": 8, "class": "weak",
                },
            }
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=encoder_error_as_weak),
                **parse_kwargs())
        unnormalized_need_more = {
            (0, scope, "candidate"): {
                "construct_result": 1, "class": "need_more",
            }
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        with self.assertRaisesRegex(
                band.MeasurementError, "construction normalization"):
            band.parse_bandtiming_output(
                bandtiming_stdout(weak=unnormalized_need_more),
                **parse_kwargs())
        invariant_error = {
            (0, scope, "candidate"): {
                "construct_result": 8, "class": "error",
            }
            for scope in band.BANDTIMING_SCOPE_ORDER
        }
        parsed = band.parse_bandtiming_output(
            bandtiming_stdout(weak=invariant_error), **parse_kwargs())
        self.assertEqual(parsed.weak_cells, ())

    def test_thermal_context_and_numeric_overflow_fail_closed(self):
        stdout = bandtiming_stdout()
        context = band_context()
        bad = copy.deepcopy(context)
        bad["thermal"]["last_monotonic_s"] = 1.5
        with self.assertRaises(band.MeasurementError):
            band.parse_bandtiming_output(
                stdout, **parse_kwargs(context=bad))
        bad = copy.deepcopy(context)
        bad["bound"]["thermal_inode"] += 1
        with self.assertRaisesRegex(band.MeasurementError, "manifest"):
            band.parse_bandtiming_output(
                stdout, **parse_kwargs(context=bad))
        with self.assertRaisesRegex(band.MeasurementError, "integer domain"):
            band.parse_bandtiming_output(
                edit_row(stdout, 0, elapsed_ns=1 << 64),
                **parse_kwargs())

        def inflate_all_timers(lines):
            for index in range(3, len(lines) - 1):
                row = dict(zip(
                    band.BANDTIMING_COLUMNS, lines[index].split(",")))
                if row["timing_eligible"] == "1":
                    row["elapsed_ns"] = "1000000000"
                    lines[index] = ",".join(
                        row[name] for name in band.BANDTIMING_COLUMNS)
        with self.assertRaisesRegex(
                band.MeasurementError, "aggregate timed interval"):
            band.parse_bandtiming_output(
                authenticated_edit(stdout, inflate_all_timers),
                **parse_kwargs())

    def test_descriptor_and_request_boundaries_match_native(self):
        accept = band.BandDescriptor(
            staircase=64000, dense_rows=0, gf256_rows=1, gf16_rows=0,
            period=1, x_mode="frozen")
        band.validate_bandtiming_dimensions(
            **{
                name: value for name, value in parse_kwargs(
                    block_count=2, candidate=accept).items()
                if name != "context"
            },
            dispatch_profile=peel_codec.TARGET_PROFILE,
            seed_policy=peel_codec.TARGET_SEED_POLICY)
        rejects = (
            band.BandDescriptor(64001, 0, 1, 0, 1, "frozen"),
            band.BandDescriptor(1, 65, 1, 0, 1, "frozen"),
            band.BandDescriptor(1, 0, 0, 1, 1, "frozen"),
            band.BandDescriptor(1, 0, 11, 0, 11, "frozen"),
            band.BandDescriptor(1, 0, 1, 5, 6, "frozen"),
            band.BandDescriptor(1, 0, 4, 1, 4, "frozen"),
            band.BandDescriptor(1, 0, 1, 0, 245, "frozen"),
            band.BandDescriptor(1, 0, 1, 0, 1, "tracking"),
        )
        for descriptor in rejects:
            with self.subTest(descriptor=descriptor):
                with self.assertRaisesRegex(ValueError, "descriptor"):
                    band.validate_bandtiming_dimensions(
                        **{
                            name: value
                            for name, value in parse_kwargs(
                                block_count=2,
                                candidate=descriptor).items()
                            if name != "context"
                        },
                        dispatch_profile=peel_codec.TARGET_PROFILE,
                        seed_policy=peel_codec.TARGET_SEED_POLICY)
        for updates in (
                {"block_count": 101},
                {"block_count": True},
                {"block_count": "16"},
                {"construction_seed": 1 << 32},
                {"systematic_cache": "pool"},
        ):
            values = parse_kwargs()
            values.update(updates)
            values.pop("context")
            if "block_count" in updates:
                values["candidate"] = default_candidate(
                    updates["block_count"])
            with self.assertRaises(ValueError):
                band.validate_bandtiming_dimensions(
                    **values,
                    dispatch_profile=peel_codec.TARGET_PROFILE,
                    seed_policy=peel_codec.TARGET_SEED_POLICY)

    def test_working_set_limit_matches_native_formula(self):
        values = parse_kwargs(block_count=100, candidate=default_candidate(100))
        values.pop("context")
        candidate = values["candidate"]
        dispatch = band.dispatch_band_descriptor(100)
        packet_slots = 100 + values["max_overhead"]
        blocks = (
            6 * packet_slots +
            4 * (100 + candidate.staircase + candidate.dense_rows +
                 candidate.heavy_rows) +
            4 * (100 + dispatch.staircase + dispatch.dense_rows +
                 dispatch.heavy_rows) +
            8 * 100 + 256
        )
        available = (
            band.BANDTIMING_WORKING_SET_BYTES_MAX -
            100 * 16384 - values["evict_bytes"])
        boundary = (available // blocks) & ~1
        values["block_bytes"] = boundary
        band.validate_bandtiming_dimensions(
            **values,
            dispatch_profile=peel_codec.TARGET_PROFILE,
            seed_policy=peel_codec.TARGET_SEED_POLICY)
        values["block_bytes"] = boundary + 2
        with self.assertRaisesRegex(ValueError, "invalid bandtiming"):
            band.validate_bandtiming_dimensions(
                **values,
                dispatch_profile=peel_codec.TARGET_PROFILE,
                seed_policy=peel_codec.TARGET_SEED_POLICY)

    def test_replay_recomputes_every_derived_summary(self):
        measurement = band.parse_bandtiming_output(
            bandtiming_stdout(), **parse_kwargs())
        receipt = measurement.as_dict()
        receipt["context"]["bound"]["cpu_affinity"].append(99)
        receipt["replicates"][0]["censored_panels"].append("forged")
        fresh = measurement.as_dict()
        self.assertEqual(fresh["context"]["bound"]["cpu_affinity"], [2])
        self.assertEqual(fresh["replicates"][0]["censored_panels"], [])
        receipt = fresh
        replayed = band.replay_bandtiming_receipt(receipt)
        self.assertEqual(replayed.stream_sha256, measurement.stream_sha256)
        legacy = copy.deepcopy(receipt)
        for replicate in legacy["replicates"]:
            for field in band._REPLICATE_CONSTRUCTION_FIELDS:
                del replicate[field]
            self.assertEqual(
                set(replicate), band._LEGACY_REPLICATE_FIELDS)
        replayed_legacy = band.replay_bandtiming_receipt(legacy)
        self.assertEqual(
            replayed_legacy.stream_sha256, measurement.stream_sha256)

        partially_enriched = copy.deepcopy(legacy)
        partially_enriched["replicates"][0][
            "candidate_construct_result"
        ] = 0
        with self.assertRaisesRegex(
                band.MeasurementError, "stale or forged"):
            band.replay_bandtiming_receipt(partially_enriched)
        incomplete_legacy = copy.deepcopy(legacy)
        del incomplete_legacy["replicates"][0][
            "encoder_candidate_overhead"
        ]
        with self.assertRaisesRegex(
                band.MeasurementError, "stale or forged"):
            band.replay_bandtiming_receipt(incomplete_legacy)
        forged_legacy = copy.deepcopy(legacy)
        forged_legacy["replicates"][0][
            "encoder_candidate_overhead"
        ] = 1
        with self.assertRaisesRegex(
                band.MeasurementError, "stale or forged"):
            band.replay_bandtiming_receipt(forged_legacy)
        forged_enriched = copy.deepcopy(receipt)
        forged_enriched["replicates"][0][
            "candidate_construct_result"
        ] = 8
        forged_enriched["replicates"][0][
            "candidate_construct_class"
        ] = "error"
        with self.assertRaisesRegex(
                band.MeasurementError, "stale or forged"):
            band.replay_bandtiming_receipt(forged_enriched)

        stale = copy.deepcopy(receipt)
        stale["arms"][0]["throughput_mbps"] += 1.0
        with self.assertRaisesRegex(
                band.MeasurementError, "stale or forged"):
            band.replay_bandtiming_receipt(stale)
        with self.assertRaisesRegex(
                band.MeasurementError, "request changed block_count"):
            band.replay_bandtiming_receipt(
                receipt, expected_request={"block_count": 17})
        malformed = copy.deepcopy(receipt)
        malformed["manifest"]["K"] = "16"
        with self.assertRaisesRegex(
                band.MeasurementError, "manifest is invalid"):
            band.replay_bandtiming_receipt(malformed)


class BandTimingProbeTests(unittest.TestCase):

    @mock.patch("band_timing_codec.os.close")
    @mock.patch("peel_codec._finalize_paired_context")
    @mock.patch("peel_codec._run_checked")
    @mock.patch("peel_codec._validate_paired_bound_context")
    @mock.patch("peel_codec._prepare_paired_context")
    def test_probe_uses_exact_long_cli_and_run_bounds_callback(
            self, prepare, validate, run, finalize, close):
        context = band_context()
        stdout = bandtiming_stdout(context=context)
        capture = SimpleNamespace(fd=91)
        prepare.return_value = (context, capture)
        validate.return_value = peel_codec._canonical_json_sha256(
            context["bound"])
        run.return_value = stdout
        finalize.return_value = context
        result = band.bandtiming_probe(
            "bench", 16, 2, default_candidate(),
            construction_seed=7, loss=0.1, loss_seed=9, schedule="iid",
            warmup_replicates=0, replicates=4, inner_reps=1,
            max_overhead=16, cache_state="warm", systematic_cache="off",
            evict_bytes=4096, context=context, required_margin=0.0)
        self.assertEqual(len(result.contrasts), 7)
        command = run.call_args.args[0]
        self.assertEqual(command[:2], ["bench", "bandtiming"])
        for option in (
                "--candidate-staircase", "--candidate-dense-rows",
                "--candidate-gf256-rows", "--candidate-gf16-rows",
                "--candidate-period", "--candidate-x-geometry",
                "--systematic-cache"):
            self.assertIn(option, command)
        self.assertIs(
            finalize.call_args.kwargs["run_bounds"],
            band._bandtiming_run_bounds)
        close.assert_called_once_with(91)

    @mock.patch("band_timing_codec.os.close")
    @mock.patch("peel_codec._run_checked")
    @mock.patch("peel_codec._validate_paired_bound_context")
    @mock.patch("peel_codec._prepare_paired_context")
    def test_probe_closes_retained_sampler_on_subprocess_failure(
            self, prepare, validate, run, close):
        context = band_context()
        prepare.return_value = (context, SimpleNamespace(fd=92))
        validate.return_value = "0" * 64
        run.side_effect = band.MeasurementError("native failed")
        with self.assertRaisesRegex(band.MeasurementError, "native failed"):
            band.bandtiming_probe(
                "bench", 16, 2, default_candidate(),
                construction_seed=7, loss=0.1, loss_seed=9, schedule="iid",
                warmup_replicates=0, replicates=4, inner_reps=1,
                max_overhead=16, cache_state="warm",
                systematic_cache="off", evict_bytes=4096,
                context=context, required_margin=0.0)
        close.assert_called_once_with(92)

    @mock.patch("band_timing_codec.os.close")
    @mock.patch("peel_codec._validate_paired_bound_context")
    @mock.patch("peel_codec._prepare_paired_context")
    def test_probe_closes_retained_sampler_on_context_validation_failure(
            self, prepare, validate, close):
        context = band_context()
        prepare.return_value = (context, SimpleNamespace(fd=93))
        validate.side_effect = band.MeasurementError("context changed")
        with self.assertRaisesRegex(band.MeasurementError, "context changed"):
            band.bandtiming_probe(
                "bench", 16, 2, default_candidate(),
                construction_seed=7, loss=0.1, loss_seed=9, schedule="iid",
                warmup_replicates=0, replicates=4, inner_reps=1,
                max_overhead=16, cache_state="warm",
                systematic_cache="off", evict_bytes=4096,
                context=context, required_margin=0.0)
        close.assert_called_once_with(93)

    def test_run_bounds_rejects_non_ascii_before_hashing(self):
        stdout = bandtiming_stdout()
        with self.assertRaisesRegex(band.MeasurementError, "not ASCII"):
            band._bandtiming_run_bounds(stdout.replace("mixed", "mixéd", 1))


if __name__ == "__main__":
    unittest.main()
