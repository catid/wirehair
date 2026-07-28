#!/usr/bin/env python3
"""Focused correctness tests for the peel-training Python tools."""

import hashlib
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from types import SimpleNamespace
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import peel_codec
import peel_direct
import peel_funnel
import peel_params
import peel_sweep
import peel_validate


def compare_stdout(
        fail=0, oh_mean=0.25, oh_sd=0.5, oh50=0.0, oh95=1.0,
        oh99=2.0, oh_max=3.0, decode_mbps=1234.5,
        construction_seed=0x1122334455667788,
        loss_seed=0x8877665544332211, pmf_digest=None,
        staircase_scale="unset", schedule="iid", loss=0.1):
    profile = fake_profile()
    if pmf_digest is None:
        pmf_digest = peel_codec.STOCK_PMF_DIGEST
    target = target_receipt(
        profile, construction_seed, loss_seed, 64, pmf_digest,
        staircase_scale, loss=loss, schedule=schedule)
    target_line = "# wh2_target," + ",".join(
        f"{name}={value}" for name, value in target.items())
    completion = " ".join(
        f"{name}={value}"
        for name, value in peel_codec.PRODUCTION_COMPLETION_REGIME.items())
    compare_arm = " ".join(
        f"{name}={value}"
        for name, value in peel_codec._COMPARE_BANNER_ARM.items())
    return (
        f"{target_line}\n"
        f"# compare: N=[64,64] trials/bb=10 loss={loss:.17g} "
        f"seed={loss_seed:#x} "
        f"max_message_mib=0 schedule={schedule} "
        f"schedule_seed={loss_seed:#x} loss_trace=common-id-v2 "
        f"{completion} {compare_arm}\n"
        "codec bb trials fail N_mean OH_mean OH_sd OH50 OH95 OH99 OH_max "
        "create_MBps encode_MBps decode_MBps recover_MBps\n"
        f"baseline 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} "
        f"{oh95} {oh99} {oh_max} 1 2 {decode_mbps} 4\n"
        f"v2 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} "
        f"{oh95} {oh99} {oh_max} 1 2 {decode_mbps} 4\n"
        f"v2_target 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} {oh95} "
        f"{oh99} {oh_max} 1 2 {decode_mbps} 4\n"
    )


def target_receipt(
        profile, construction_seed, loss_seed, block_bytes=64,
        pmf_digest=None, staircase_scale="unset", loss=0.1, schedule="iid"):
    if pmf_digest is None:
        pmf_digest = peel_codec.STOCK_PMF_DIGEST
    return {
        "profile": peel_codec.TARGET_PROFILE,
        "contract_id": peel_codec.TARGET_CONTRACT["contract_id"],
        "contract_sha256": peel_codec.TARGET_CONTRACT["contract_sha256"],
        "precode_contract":
            str(peel_codec.TARGET_CONTRACT["precode_contract"]),
        "packet_contract":
            str(peel_codec.TARGET_CONTRACT["packet_contract"]),
        "architecture": peel_codec.TARGET_CONTRACT["architecture"],
        **peel_codec.TARGET_MEASUREMENT_POLICY,
        "N": str(profile.block_count),
        "bb": str(block_bytes),
        "staircase": str(profile.staircase),
        "dense_rows": str(profile.dense_rows),
        "heavy_rows": str(profile.heavy_rows),
        "completion": profile.completion,
        "gf256_rows": str(profile.gf256_rows),
        "gf16_rows": str(profile.gf16_rows),
        "period": str(profile.period),
        "geometry": profile.geometry,
        "residue_schedule": profile.residue_schedule,
        "residue_skew": str(profile.residue_skew),
        "source_hits": str(profile.source_hits),
        "target_mean": peel_codec._target_mean_spec(
            profile.block_count, profile.staircase, profile.source_hits,
            staircase_scale),
        "mix_count": str(profile.mix_count),
        "packet_peel_seed": str(
            ((construction_seed & 0xffffffff) ^
             ((construction_seed >> 32) & 0xffffffff)) & 0xffffffff),
        "construction_seed": f"0x{construction_seed:x}",
        "loss_rate": f"{loss:.17g}",
        "loss_seed": f"0x{loss_seed:x}",
        "schedule": schedule,
        "pmf_sha256": pmf_digest,
        "pmf_encoding": peel_codec.TARGET_CONTRACT["pmf_encoding"],
        "staircase_scale": staircase_scale,
    }


def metrics(
        construction_seed, loss_seed=None, fail=0, decode_mbps=1000.0,
        profile=None, block_bytes=64, pmf_digest=None,
        staircase_scale="unset"):
    if loss_seed is None:
        loss_seed = construction_seed
    if profile is None:
        profile = fake_profile()
    return peel_codec.RecoveryMetrics(
        construction_seed=construction_seed,
        loss_seed=loss_seed,
        target_receipt=target_receipt(
            profile, construction_seed, loss_seed, block_bytes,
            pmf_digest, staircase_scale),
        fail=fail,
        oh_mean=0.25,
        oh_sd=0.5,
        oh50=0.0,
        oh95=1.0,
        oh99=2.0,
        oh_max=3.0,
        decode_mbps=decode_mbps,
    )


def metrics_for_probe(args, kwargs, *, fail=0, decode_mbps=1000.0):
    peel_weights = kwargs.get("peel_weights")
    pmf_digest = (
        peel_codec.STOCK_PMF_DIGEST if peel_weights is None else
        peel_codec.pmf_sha256(peel_weights))
    scale = kwargs.get("degree_scale")
    return metrics(
        kwargs["construction_seed"],
        kwargs["loss_seed"],
        fail=fail,
        decode_mbps=decode_mbps,
        profile=kwargs["native_profile"],
        block_bytes=args[3],
        pmf_digest=pmf_digest,
        staircase_scale=(
            "unset" if scale is None else
            peel_codec._canonical_staircase_scale_spec(scale)),
    )


def fake_profile(block_count=64, target_mean=None):
    staircase = peel_codec._dispatch_staircase_count(block_count)
    source_hits = peel_codec._dispatch_source_hits(block_count)
    if target_mean is None:
        target_mean = (
            block_count * min(source_hits, staircase) / staircase
        )
    return peel_codec.StockProfile(
        block_count=block_count,
        target_profile=peel_codec.TARGET_PROFILE,
        contract_id=peel_codec.TARGET_CONTRACT["contract_id"],
        contract_sha256=peel_codec.TARGET_CONTRACT["contract_sha256"],
        precode_contract=5,
        packet_contract=4,
        architecture="smallband100-d4",
        staircase=staircase,
        dense_rows=4,
        heavy_rows=12,
        source_hits=source_hits,
        completion="mixed",
        gf256_rows=10,
        gf16_rows=2,
        period=244,
        geometry="frozen",
        residue_schedule="constant",
        residue_skew=0,
        mix_count=3,
        target_mean=target_mean,
        native_pmf_sha256=peel_codec.STOCK_PMF_DIGEST,
        pmf_encoding="wirehair-v2-peel-spec-v1",
        pmf=(1.0 / 64.0,) * 64,
    )


def probe_kwargs(
        profile=None, construction_seed=0x1122334455667788,
        loss_seed=0x8877665544332211, loss=0.1, schedule="iid"):
    return {
        "native_profile": profile or fake_profile(),
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "loss": loss,
        "schedule": schedule,
        "construction_seed": construction_seed,
        "loss_seed": loss_seed,
    }


def paired_parse_kwargs(
        context=None, schedule="iid", profile=None,
        candidate_pmf=None, block_bytes=64, degree_scale=None,
        construction_seed=9, loss_seed=10, warmups=0, replicates=4,
        max_overhead=16, required_margin=0.0):
    profile = profile or fake_profile()
    return {
        "block_count": profile.block_count,
        "block_bytes": block_bytes,
        "candidate_pmf": (
            [0.25, 0.75] if candidate_pmf is None
            else list(candidate_pmf)),
        "degree_scale": degree_scale,
        "native_profile": profile,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "construction_seed": construction_seed,
        "loss": 0.1,
        "loss_seed": loss_seed,
        "schedule": schedule,
        "warmup_replicates": warmups,
        "replicates": replicates,
        "inner_reps": 1,
        "max_overhead": max_overhead,
        "cache_state": "warm",
        "evict_bytes": 4096,
        "context": context or paired_context(),
        "required_margin": required_margin,
    }


def table_kwargs():
    return {
        "construction_seed_base": 1,
        "loss_seed_base": 2,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "loss": 0.1,
        "schedule": "iid",
    }


def direct_settings(kmin=64, kmax=64):
    return {
        "kmin": kmin,
        "kmax": kmax,
        "screen": 1,
        "refine": 0,
        "gate_trials": 1,
        "gate_block_bytes": 64,
        "screen_paired_replicates": 4,
        "paired_replicates": 4,
        "rank_block_bytes": 64,
        "rank_top": 1,
        "paired_context": "/tmp/paired-context.json",
        "paired_warmups": 0,
        "paired_inner_reps": 1,
        "max_overhead": 16,
        "cache_state": "warm",
        "evict_bytes": 4096,
        "rank_margin": 0.0,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "loss": 0.1,
        "schedule": "iid",
    }


def sweep_settings(block_counts=(2,), paired_replicates=4):
    return {
        "proxy_k_ladder": list(block_counts),
        "paired_replicates": paired_replicates,
        "paired_context": "/tmp/paired-context.json",
        "paired_warmups": 0,
        "paired_inner_reps": 1,
        "max_overhead": 16,
        "cache_state": "warm",
        "evict_bytes": 4096,
        "rank_margin": 0.0,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "loss": 0.1,
        "schedule": "iid",
        "proxy_cost_model": peel_codec.PROXY_COST_MODEL,
        "proxy_measure_regime": dict(peel_codec.PROXY_MEASURE_REGIME),
        "proxy_ordering": peel_codec.PROXY_ORDERING_PROTOCOL,
        "allow_unverified_cost_model": True,
        "search_box": peel_codec.SEARCH_BOX_PROTOCOL,
    }


def exact_cli_args(construction_seed=9, loss_seed=10):
    return [
        "--target-profile", "dispatch-v1",
        "--seed-policy", "raw",
        "--schedule", "iid",
        "--loss", "0.1",
        "--construction-seed", str(construction_seed),
        "--loss-seed", str(loss_seed),
    ]


def peelpmf_stdout(
        probabilities, *, staircase=10, source_hits=2, target_mean=12.8):
    profile = fake_profile(target_mean=target_mean)
    values = {
        "N": 64,
        "target_profile": profile.target_profile,
        "contract_id": profile.contract_id,
        "contract_sha256": profile.contract_sha256,
        "precode_contract": profile.precode_contract,
        "packet_contract": profile.packet_contract,
        "architecture": profile.architecture,
        "degrees": 64,
        "staircase": staircase,
        "dense_rows": profile.dense_rows,
        "heavy_rows": profile.heavy_rows,
        "source_hits": source_hits,
        "completion": profile.completion,
        "gf256_rows": profile.gf256_rows,
        "gf16_rows": profile.gf16_rows,
        "period": profile.period,
        "geometry": profile.geometry,
        "residue_schedule": profile.residue_schedule,
        "residue_skew": profile.residue_skew,
        "mix_count": profile.mix_count,
        "target_mean": target_mean,
        "pmf_sha256": peel_codec.STOCK_PMF_DIGEST,
        "pmf_encoding": profile.pmf_encoding,
    }
    lines = [
        "# peelpmf," + ",".join(
            f"{name}={value}" for name, value in values.items()),
        "degree,probability",
    ]
    lines.extend(
        f"{degree},{probability:.17g}"
        for degree, probability in enumerate(probabilities, 1))
    return "\n".join(lines) + "\n"


def paired_context(start=0.5, end=3.0, *, evict_bytes=4096):
    middle = (start + end) / 2.0
    payload = thermal_csv(
        thermal_row(start), thermal_row(middle), thermal_row(end))
    thermal = peel_codec._thermal_window_from_csv(
        payload,
        int(round(start * 1e9)),
        int(round(end * 1e9)),
        sampling_interval_ms=1000)
    return {
        "schema": peel_codec.PAIRED_CONTEXT_SCHEMA,
        "bound": {
            "cpu_model": "test-cpu",
            "cpu_affinity": [2],
            "cpu_governors": {"2": "performance"},
            "cache_state": "warm",
            "evict_bytes": evict_bytes,
            "thermal_source": "/tmp/test-sampler.csv",
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


def write_unfinalized_paired_context(path, thermal_source):
    """Write a strict pre-launch context fixture with one protected sampler."""
    context = paired_context()
    for name in peel_codec._PAIRED_BOUND_CAPTURE_FIELDS:
        context["bound"].pop(name)
    context["bound"]["thermal_source"] = os.path.realpath(thermal_source)
    context["thermal"] = None
    peel_codec.write_json_atomic(path, context)
    return context


def peeltiming_stdout(
        *, profile=None, candidate_pmf=(0.25, 0.75), candidate_ns=80,
        identity_ns=100, aa_left_ns=100, aa_right_ns=100,
        warmups=0, replicates=4, construction_seed=9, loss_seed=10,
        schedule="iid", context=None, required_margin=0.0, mutate_row=None,
        block_bytes=64, degree_scale=None, max_overhead=16,
        evict_bytes=4096, started_ns=1_000_000_000,
        finished_ns=2_000_000_000):
    if profile is None:
        profile = fake_profile()
    if context is None:
        context = paired_context(
            max(0.0, started_ns / 1e9 - 0.5),
            finished_ns / 1e9 + 1.0,
            evict_bytes=evict_bytes)
    total = warmups + replicates
    expected_rows = 16 * total
    candidate_digest = peel_codec._paired_pmf_sha256(candidate_pmf)
    identity_digest = peel_codec._paired_pmf_sha256(profile.pmf)
    values = {
        "schema": peel_codec._PEELTIMING_SCHEMA,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "contract_id": peel_codec.TARGET_CONTRACT["contract_id"],
        "K": profile.block_count,
        "bb": block_bytes,
        "S": profile.staircase,
        "D2": profile.dense_rows,
        "H": profile.heavy_rows,
        "construction_seed_base": construction_seed,
        "construction_seed_derivation":
            peel_codec.PAIRED_CONSTRUCTION_SEED_DERIVATION,
        "semantic_seed_derivation":
            "base-plus-attempt-mod2^64-skip-measured-alias-v1",
        "loss": "0.10000000000000001",
        "loss_seed_base": loss_seed,
        "loss_seed_derivation": peel_codec.PAIRED_LOSS_SEED_DERIVATION,
        "message_seed_policy": "replicate-loss-seed-v1",
        "schedule": schedule,
        "loss_model": peel_codec.PAIRED_LOSS_MODEL,
        "trace_encoding": "wirehair-wh2-peeltiming-loss-trace-v1",
        "panels": "candidate_control+identity_aa",
        "candidate_pmf_sha256": candidate_digest,
        "candidate_pmf_encoding": "binary64-17g-comma-v1",
        "candidate_scale_requested": (
            "identity" if degree_scale is None else f"{degree_scale:.17g}"),
        "candidate_scale_effective": (
            f"{profile.target_mean:.17g}"
            if degree_scale is None else f"{degree_scale:.17g}"),
        "identity_pmf_sha256": identity_digest,
        "identity_pmf_encoding": "stock-recovered-explicit-17g-v1",
        "identity_scale_effective": f"{profile.target_mean:.17g}",
        "warmup_replicates": warmups,
        "replicates": replicates,
        "slots_per_panel": 8,
        "panels_per_replicate": 2,
        "order": peel_codec.PAIRED_ORDER,
        "label_swap": "alternating",
        "inner_reps": 1,
        "max_overhead": max_overhead,
        "cache_state": "warm",
        "evict_bytes": evict_bytes,
        "payload_alignment": 64,
        "prefault": 1,
        "cpu_affinity_policy": "first-allowed-affinity-v1",
        "payload": "common-source-role-encoded-consistent-rhs-v1",
        "timing_scope": "decoder_solve_only",
        "timing_prefix": "common-max-recovery-overhead-v1",
        "recovery_scope": "full_encode_decode_recover_memcmp",
        "weak_seed_policy": "balanced-timing-ineligible-rows-v1",
        "hook_path": "explicit-tls-peel-pmf-and-scale-v1",
        "no_override_scope": "untimed-semantic-only",
        "system_build": "per_replicate_per_scale_outside_timer",
        "startup_amortization":
            "per-replicate-build-and-recovery-preflight-excluded-v1",
        "slot_prewarm":
            "validated-plus-conditioning-matching-role-solves-"
            "same-cpu-before-cache-v1",
        "context_sha256": peel_codec._canonical_json_sha256(context["bound"]),
        "uncertainty": peel_codec.PAIRED_UNCERTAINTY_PROTOCOL,
        "required_margin": f"{required_margin:.17g}",
        "margin_rule":
            "upper-log-cost-lt-negative-required-margin-and-aa-floor-v1",
        "clock_domain": "posix-clock-monotonic-ns-v1",
        "stream_hash_scope": "body-plus-done-prefix-v1",
        "started_monotonic_ns": started_ns,
        "expected_rows": expected_rows,
    }
    semantic_construction = construction_seed
    semantic_loss = loss_seed
    semantic_trace = hashlib.sha256(b"semantic-trace").hexdigest()
    message_digest = hashlib.sha256(b"message").hexdigest()
    system_digest = hashlib.sha256(b"system").hexdigest()
    packet_digest = hashlib.sha256(b"packets").hexdigest()
    payload_digest = hashlib.sha256(b"payload").hexdigest()
    solve_digest = hashlib.sha256(b"solve").hexdigest()
    intermediate_bytes = (
        profile.block_count + profile.staircase +
        profile.dense_rows + profile.heavy_rows
    ) * block_bytes
    semantic = {
        "timed": 0,
        "construction_seed": semantic_construction,
        "seed_attempt": 0,
        "seed_attempt_cap": 256,
        "loss_seed": semantic_loss,
        "trace_sha256": semantic_trace,
        "message_sha256": message_digest,
        "identity_pmf_sha256": identity_digest,
        "identity_pmf_encoding": "stock-recovered-explicit-17g-v1",
        "identity_scale_effective": f"{profile.target_mean:.17g}",
        "nohook_result": 0,
        "nohook_recovery_ok": 1,
        "nohook_overhead": 1,
        "identity_result": 0,
        "identity_recovery_ok": 1,
        "identity_overhead": 1,
        "nohook_system_sha256": system_digest,
        "identity_system_sha256": system_digest,
        "system_equal": 1,
        "nohook_packet_rows_sha256": packet_digest,
        "identity_packet_rows_sha256": packet_digest,
        "packet_rows_equal": 1,
        "nohook_payload_sha256": payload_digest,
        "identity_payload_sha256": payload_digest,
        "payload_equal": 1,
        "nohook_direct_result": 0,
        "identity_direct_result": 0,
        "nohook_intermediate_bytes": intermediate_bytes,
        "identity_intermediate_bytes": intermediate_bytes,
        "nohook_solve_sha256": solve_digest,
        "identity_solve_sha256": solve_digest,
        "solve_equal": 1,
        "full_recovery_equal": 1,
        "pass": 1,
    }
    lines = [
        "# peeltiming," + ",".join(
            f"{name}={values[name]}"
            for name in peel_codec.PEELTIMING_MANIFEST_FIELDS),
        "# peel_semantic," + ",".join(
            f"{name}={semantic[name]}"
            for name in peel_codec.PEELTIMING_SEMANTIC_FIELDS),
        ",".join(peel_codec.PEELTIMING_COLUMNS),
    ]
    labels = peel_codec.PAIRED_ORDER
    row_index = 0
    for replicate in range(total):
        construction = (
            construction_seed ^
            (((replicate + 1) * 0xd1b54a32d192ed03) & 0xffffffffffffffff)
        )
        loss = (
            loss_seed ^
            (((replicate + 1) * 0x9e3779b97f4a7c15) & 0xffffffffffffffff)
        )
        trace = hashlib.sha256(f"trace-{replicate}".encode()).hexdigest()
        for panel in ("candidate_control", "identity_aa"):
            for slot, label in enumerate(labels):
                role = peel_codec._expected_peeltiming_roles(
                    panel, replicate, label)
                elapsed = {
                    "candidate": candidate_ns,
                    "identity": identity_ns,
                    "identity_a": aa_left_ns,
                    "identity_b": aa_right_ns,
                }[role]
                row = {
                    "replicate": replicate,
                    "measured": int(replicate >= warmups),
                    "panel": panel,
                    "slot": slot,
                    "pair": slot // 2,
                    "label": label,
                    "role": role,
                    "label_swap": replicate & 1,
                    "construction_seed": construction,
                    "loss_seed": loss,
                    "trace_sha256": trace,
                    "common_overhead": 1,
                    "arm_overhead": 1,
                    "recovery_result": 0,
                    "recovery_ok": 1,
                    "timing_eligible": 1,
                    "preflight_result": 0,
                    "timing_result": 0,
                    "outcome_stable": 1,
                    "elapsed_ns": elapsed,
                    "inner_reps": 1,
                    "saturated": 0,
                    "cpu_before": 2,
                    "cpu_after": 2,
                    "cpu_migrated": 0,
                    "minflt_delta": 0,
                    "majflt_delta": 0,
                    "fault_contaminated": 0,
                    "inactivated": 3,
                    "binary_def": 1,
                    "heavy_gain": 1,
                    "block_xors": 10,
                    "block_muladds": 5,
                    "build_ns_sum": 10,
                    "peel_ns_sum": 10,
                    "project_ns_sum": 10,
                    "residual_ns_sum": 10,
                    "backsub_ns_sum": 10,
                    "source_bytes": profile.block_count * block_bytes,
                    "packet_payload_bytes":
                        (profile.block_count + 1) * block_bytes,
                    "intermediate_bytes": (
                        profile.block_count + profile.staircase +
                        profile.dense_rows + profile.heavy_rows
                    ) * block_bytes,
                }
                if mutate_row is not None:
                    mutate_row(row_index, row)
                lines.append(",".join(
                    str(row[name]) for name in peel_codec.PEELTIMING_COLUMNS))
                row_index += 1
    prefix = "\n".join(lines) + "\n"
    done_prefix = (
        f"# peeltiming_done,complete=1,rows={expected_rows},"
        f"finished_monotonic_ns={finished_ns},stream_sha256="
    )
    stream_sha = hashlib.sha256(
        (prefix + done_prefix).encode("ascii")).hexdigest()
    return prefix + done_prefix + stream_sha + "\n"


def authenticated_peeltiming_edit(stdout, edit):
    """Edit a native fixture while preserving its stream authentication."""
    lines = stdout.splitlines()
    edit(lines)
    completion = dict(
        field.split("=", 1)
        for field in lines[-1][len("# peeltiming_done,"):].split(",")
    )
    prefix = "\n".join(lines[:-1]) + "\n"
    done_prefix = (
        f"# peeltiming_done,complete={completion['complete']},"
        f"rows={completion['rows']},"
        f"finished_monotonic_ns={completion['finished_monotonic_ns']},"
        "stream_sha256="
    )
    lines[-1] = (
        done_prefix +
        hashlib.sha256((prefix + done_prefix).encode("ascii")).hexdigest()
    )
    return "\n".join(lines) + "\n"


def skip_peeltiming_row(row, *, common_overhead):
    """Turn one fixture row into an authenticated recovery-only placeholder."""
    row.update({
        "common_overhead": common_overhead,
        "timing_eligible": 0,
        "preflight_result": -1,
        "timing_result": -1,
        "outcome_stable": 0,
        "elapsed_ns": 0,
        "inner_reps": 0,
        "saturated": 0,
        "cpu_before": -1,
        "cpu_after": -1,
        "cpu_migrated": -1,
        "minflt_delta": 0,
        "majflt_delta": 0,
        "fault_contaminated": 0,
        "inactivated": 0,
        "binary_def": 0,
        "heavy_gain": 0,
        "block_xors": 0,
        "block_muladds": 0,
        "build_ns_sum": 0,
        "peel_ns_sum": 0,
        "project_ns_sum": 0,
        "residual_ns_sum": 0,
        "backsub_ns_sum": 0,
        "packet_payload_bytes": 0,
        "intermediate_bytes": 0,
    })


def thermal_row(
        monotonic, *, dimm_read_errors="0", missing_dimm=False,
        edac_ce="0", edac_ue="0", cpu="60.0", dimm="40.0"):
    fields = [
        "2026-07-28T00:00:00Z",
        monotonic if isinstance(monotonic, str) else f"{monotonic:.9f}",
        "100.0",
        "4000.0",
        cpu,
        *(dimm for unused in range(8)),
        dimm_read_errors,
        "1.0",
        "1.0",
        "1.0",
        edac_ce,
        edac_ue,
    ]
    if missing_dimm:
        fields[5] = ""
    return ",".join(fields) + "\n"


def thermal_csv(*rows):
    return (
        ",".join(peel_codec._THERMAL_CSV_COLUMNS) + "\n" +
        "".join(rows)
    ).encode("ascii")


def _legacy_complete_entry(block_count=64, mode="shipped"):
    profile = fake_profile(block_count)
    shipped = mode == "shipped"
    coordinates = {
        "K": block_count,
        # Direct search never supplies a staircase-scale override.
        "scale": -1.0,
        "p1": 100 if shipped else 200,
        "tilt": 0,
        "dmax": 64,
        "absorb": 100,
        **({"reverted_to_shipped": True} if shipped else {}),
    }
    selected_pmf = (
        list(profile.pmf) if shipped else
        peel_codec.family(
            profile.pmf, coordinates["p1"], coordinates["tilt"],
            coordinates["dmax"], coordinates["absorb"])
    )
    construction_seed = peel_codec.derive_seed(
        1, "direct-search", block_count, "rank", 10, 64, "construction")
    loss_seed = peel_codec.derive_seed(
        2, "direct-search", block_count, "rank", 10, 64, "loss")
    selected_digest = (
        peel_codec.STOCK_PMF_DIGEST if shipped else
        peel_codec.pmf_sha256(selected_pmf))
    measurement = metrics(
        construction_seed, loss_seed, profile=profile,
        pmf_digest=selected_digest)
    shipped_measurement = metrics(
        construction_seed, loss_seed,
        decode_mbps=900.0 if not shipped else measurement.decode_mbps,
        profile=profile)
    return {
        **coordinates,
        **measurement.as_dict(),
        "goodput": measurement.goodput(block_count),
        "native": profile.as_dict(),
        "peel_pmf": selected_pmf,
        "seconds": 0.0,
        "probes": 1,
        "search_receipt": peel_codec.make_search_receipt(
            measurement,
            mode=mode,
            goodput=measurement.goodput(block_count),
            trials=10,
            block_bytes=64,
            search_kind="direct-real-codec",
            construction_seed_base=1,
            loss_seed_base=2,
            seed_domain="direct-search",
            coordinates={
                name: coordinates[name]
                for name in ("scale", "p1", "tilt", "dmax", "absorb")
            },
            peel_pmf=selected_pmf,
            shipped_control=shipped_measurement,
            shipped_goodput=shipped_measurement.goodput(block_count),
            context={
                "warm_start": None,
                "sampling_seed": peel_codec.derive_seed(
                    1, "direct-search", block_count, "sampling"),
                "screen": 1,
                "refine": 0,
                "gate_trials": 1,
                "gate_block_bytes": 64,
                "rank_top": 1,
            },
        ),
    }


def _legacy_funnel_entry(block_count=64, mode="shipped", native_shape=False):
    """Build one schema-valid funnel receipt, including shipped controls."""
    entry = complete_entry(block_count, mode)
    profile = fake_profile(block_count)
    del entry["probes"]
    entry.update({
        "S": profile.staircase,
        "source_hits": profile.source_hits,
        "target_mean": profile.target_mean,
        "screen": 3000,
        "screen_cells": 1000,
        "finals": 16,
        "rejected": 0,
        "real_trials": 10,
    })
    construction_seed = peel_codec.derive_seed(
        1, "funnel-search", block_count, "rank", 10, 4096, "construction")
    loss_seed = peel_codec.derive_seed(
        2, "funnel-search", block_count, "rank", 10, 4096, "loss")
    receipt = entry["search_receipt"]
    for target in (entry, receipt, receipt["shipped_control"]):
        target["construction_seed"] = construction_seed
        target["loss_seed"] = loss_seed
    candidate_digest = (
        peel_codec.STOCK_PMF_DIGEST if mode == "shipped" else
        peel_codec.pmf_sha256(receipt["peel_pmf"]))
    entry["target_receipt"] = receipt["target_receipt"] = target_receipt(
        profile, construction_seed, loss_seed, 4096, candidate_digest)
    receipt["shipped_control"]["target_receipt"] = target_receipt(
        profile, construction_seed, loss_seed, 4096)
    receipt["block_bytes"] = 4096
    receipt["search_kind"] = "unverified-proxy-funnel"
    receipt["seed_domain"] = "funnel-search"
    scale_max_centi = peel_funnel.search_box(profile)[0][2]
    receipt["context"] = {
        "proxy_cost_model": peel_funnel.PROXY_COST_MODEL,
        "proxy_measure_regime": dict(peel_funnel.PROXY_MEASURE_REGIME),
        "proxy_ordering": peel_funnel.PROXY_ORDERING_PROTOCOL,
        "search_box": peel_funnel.SEARCH_BOX_PROTOCOL,
        "scale_centi": [0, scale_max_centi],
        "warm_start": None,
        "sampling_seed": peel_codec.derive_seed(
            1, "funnel-search", block_count, "sampling"),
        "screen": 3000,
        "refine": 400,
        "finals": 16,
        "screen_cells": 1000,
        "gate_trials": 25,
        "gate_block_bytes": 64,
        "rank_top": 3,
        "threads": 64,
        "batch": 60,
        "cell_base": 900_000_000,
    }
    if mode == "trained":
        scale = min(12.0, scale_max_centi / 100.0)
        scale_spec = f"{scale:.17g}"
        entry["scale"] = receipt["coordinates"]["scale"] = scale
        if native_shape:
            entry["p1"] = receipt["coordinates"]["p1"] = 100
            entry["peel_pmf"] = list(profile.pmf)
            receipt["peel_pmf"] = list(profile.pmf)
            receipt["peel_pmf_sha256"] = peel_codec.pmf_sha256(profile.pmf)
            receipt["mode"] = "scale-only"
        candidate_digest = (
            peel_codec.STOCK_PMF_DIGEST if native_shape else
            peel_codec.pmf_sha256(receipt["peel_pmf"]))
        entry["target_receipt"] = receipt["target_receipt"] = target_receipt(
            profile, construction_seed, loss_seed, 4096, candidate_digest,
            scale_spec)
        receipt["shipped_control"]["target_receipt"] = target_receipt(
            profile, construction_seed, loss_seed, 4096,
            peel_codec.STOCK_PMF_DIGEST)
    return entry


def paired_fixture_measurement(
        profile, candidate_pmf, construction_seed, loss_seed, *,
        block_bytes=64, degree_scale=None, candidate_ns=80, identity_ns=100,
        required_margin=0.0, started_ns=1_000_000_000,
        finished_ns=2_000_000_000, replicates=4, mutate_row=None):
    context = paired_context(
        max(0.0, started_ns / 1e9 - 0.5),
        finished_ns / 1e9 + 1.0)
    stdout = peeltiming_stdout(
        profile=profile,
        candidate_pmf=candidate_pmf,
        candidate_ns=candidate_ns,
        identity_ns=identity_ns,
        construction_seed=construction_seed,
        loss_seed=loss_seed,
        context=context,
        block_bytes=block_bytes,
        degree_scale=degree_scale,
        required_margin=required_margin,
        replicates=replicates,
        mutate_row=mutate_row,
        started_ns=started_ns,
        finished_ns=finished_ns,
    )
    return peel_codec.parse_peeltiming_output(
        stdout,
        **paired_parse_kwargs(
            context=context,
            profile=profile,
            candidate_pmf=candidate_pmf,
            block_bytes=block_bytes,
            degree_scale=degree_scale,
            construction_seed=construction_seed,
            loss_seed=loss_seed,
            replicates=replicates,
            required_margin=required_margin,
        ),
    )


def complete_entry(block_count=64, mode="shipped"):
    """Build one exact v4 direct-search entry from replayable native rows."""
    profile = fake_profile(block_count)
    stock_coordinates = {
        "scale": -1.0, "p1": 100, "tilt": 0,
        "dmax": 64, "absorb": 100,
    }
    coordinates = (
        dict(stock_coordinates) if mode == "shipped" else
        {
            "scale": -1.0, "p1": 200, "tilt": 0,
            "dmax": 64, "absorb": 100,
        }
    )
    evaluated_coordinates = dict(coordinates)
    evaluated_pmf = (
        list(profile.pmf) if mode == "shipped" else
        peel_codec.family(
            profile.pmf, coordinates["p1"], coordinates["tilt"],
            coordinates["dmax"], coordinates["absorb"])
    )
    construction_seed = peel_codec.derive_seed(
        1, "direct-search", block_count, "rank", 4, 64, "construction")
    loss_seed = peel_codec.derive_seed(
        2, "direct-search", block_count, "rank", 4, 64, "loss")
    measurement = paired_fixture_measurement(
        profile, evaluated_pmf, construction_seed, loss_seed,
        candidate_ns=100 if mode == "shipped" else 80)
    stock_control = paired_fixture_measurement(
        profile, profile.pmf, construction_seed, loss_seed,
        candidate_ns=100, started_ns=3_000_000_000,
        finished_ns=4_000_000_000)
    selected = (
        measurement.identity if mode == "shipped"
        else measurement.candidate)
    search_context = {
        "warm_start": None,
        "sampling_seed": peel_codec.derive_seed(
            1, "direct-search", block_count, "sampling"),
        "screen": 1,
        "refine": 0,
        "gate_trials": 1,
        "gate_block_bytes": 64,
        "screen_paired_replicates": 4,
        "rank_top": 1,
    }
    return {
        "K": block_count,
        **coordinates,
        **({"reverted_to_shipped": True} if mode == "shipped" else {}),
        **selected.as_dict(),
        "goodput": selected.goodput(block_count),
        "native": profile.as_dict(),
        "peel_pmf": (
            list(profile.pmf) if mode == "shipped"
            else list(evaluated_pmf)),
        "seconds": 0.0,
        "probes": 1,
        "search_receipt": peel_codec.make_paired_search_receipt(
            measurement,
            mode=mode,
            block_count=block_count,
            block_bytes=64,
            search_kind="direct-real-codec",
            construction_seed_base=1,
            loss_seed_base=2,
            seed_domain="direct-search",
            coordinates=coordinates,
            peel_pmf=(
                list(profile.pmf) if mode == "shipped"
                else evaluated_pmf),
            evaluated_coordinates=evaluated_coordinates,
            evaluated_pmf=evaluated_pmf,
            stock_control_measurement=(
                None if mode == "shipped" else stock_control),
            context=search_context,
        ),
    }


def funnel_entry(
        block_count=64, mode="shipped", native_shape=False,
        scale_override=None, trained_p1=200, trained_tilt=0,
        trained_dmax=64, trained_absorb=100):
    """Build one exact v4 funnel entry from replayable native rows."""
    profile = fake_profile(block_count)
    scale_max_centi = peel_funnel.search_box(profile)[0][2]
    if scale_override is None:
        scale_centi = min(1200, scale_max_centi)
        while (scale_centi > 0 and
               peel_codec._scale_realizes_legacy_edges(
                   scale_centi / 100.0, block_count, profile.staircase,
                   profile.source_hits)):
            scale_centi -= 1
        scale = scale_centi / 100.0
    else:
        scale = float(scale_override)
    actual_mode = (
        "scale-only" if native_shape or mode == "scale-only" else mode)
    if actual_mode == "shipped":
        coordinates = {
            "scale": -1.0, "p1": 100, "tilt": 0,
            "dmax": 64, "absorb": 100,
        }
    elif actual_mode == "scale-only":
        coordinates = {
            "scale": float(scale), "p1": 100, "tilt": 0,
            "dmax": 64, "absorb": 100,
        }
    else:
        coordinates = {
            "scale": float(scale), "p1": trained_p1,
            "tilt": trained_tilt,
            "dmax": trained_dmax, "absorb": trained_absorb,
        }
    evaluated_coordinates = dict(coordinates)
    evaluated_pmf = (
        list(profile.pmf)
        if actual_mode in ("shipped", "scale-only") else
        peel_codec.family(
            profile.pmf, coordinates["p1"], coordinates["tilt"],
            coordinates["dmax"], coordinates["absorb"])
    )
    if (actual_mode == "trained" and
            peel_codec._native_peel_cdf_signature(
                evaluated_pmf, block_count) ==
            peel_codec._native_peel_cdf_signature(
                profile.pmf, block_count)):
        actual_mode = "scale-only"
        coordinates = {
            "scale": float(scale), "p1": 100, "tilt": 0,
            "dmax": 64, "absorb": 100,
        }
        evaluated_coordinates = dict(coordinates)
        evaluated_pmf = list(profile.pmf)
    construction_seed = peel_codec.derive_seed(
        1, "funnel-search", block_count, "rank", 4, 4096, "construction")
    loss_seed = peel_codec.derive_seed(
        2, "funnel-search", block_count, "rank", 4, 4096, "loss")
    measurement = paired_fixture_measurement(
        profile, evaluated_pmf, construction_seed, loss_seed,
        block_bytes=4096,
        degree_scale=(
            None if coordinates["scale"] == -1.0 else coordinates["scale"]),
        candidate_ns=100 if actual_mode == "shipped" else 80,
    )
    stock_control = paired_fixture_measurement(
        profile, profile.pmf, construction_seed, loss_seed,
        block_bytes=4096, candidate_ns=100,
        started_ns=3_000_000_000, finished_ns=4_000_000_000)
    selected = (
        measurement.identity if actual_mode == "shipped"
        else measurement.candidate)
    screen, unused_refine, finals, screen_cells = peel_codec._sweep_budget(
        block_count)
    search_context = {
        "proxy_cost_model": peel_funnel.PROXY_COST_MODEL,
        "proxy_measure_regime": dict(peel_funnel.PROXY_MEASURE_REGIME),
        "proxy_ordering": peel_funnel.PROXY_ORDERING_PROTOCOL,
        "search_box": peel_funnel.SEARCH_BOX_PROTOCOL,
        "scale_centi": [0, scale_max_centi],
        "warm_start": None,
        "sampling_seed": peel_codec.derive_seed(
            1, "funnel-search", block_count, "sampling"),
        "screen": screen,
        "refine": unused_refine,
        "finals": finals,
        "screen_cells": screen_cells,
        "gate_trials": 25,
        "gate_block_bytes": 64,
        "rank_top": 3,
        "threads": 64,
        "batch": 60,
        "cell_base": 900_000_000,
    }
    return {
        "K": block_count,
        **coordinates,
        **({"reverted_to_shipped": True}
           if actual_mode == "shipped" else {}),
        **selected.as_dict(),
        "goodput": selected.goodput(block_count),
        "native": profile.as_dict(),
        "peel_pmf": (
            list(profile.pmf) if actual_mode == "shipped"
            else list(evaluated_pmf)),
        "S": profile.staircase,
        "source_hits": profile.source_hits,
        "target_mean": profile.target_mean,
        "seconds": 0.0,
        "screen": screen,
        "screen_cells": screen_cells,
        "finals": finals,
        "rejected": 0,
        "paired_replicates": 4,
        "search_receipt": peel_codec.make_paired_search_receipt(
            measurement,
            mode=actual_mode,
            block_count=block_count,
            block_bytes=4096,
            search_kind="unverified-proxy-funnel",
            construction_seed_base=1,
            loss_seed_base=2,
            seed_domain="funnel-search",
            coordinates=coordinates,
            peel_pmf=(
                list(profile.pmf) if actual_mode == "shipped"
                else evaluated_pmf),
            evaluated_coordinates=evaluated_coordinates,
            evaluated_pmf=evaluated_pmf,
            stock_control_measurement=(
                None if actual_mode == "shipped" else stock_control),
            context=search_context,
        ),
    }


class PeelCodecTests(unittest.TestCase):
    def setUp(self):
        peel_codec._STOCK_CACHE.clear()

    def test_peeltiming_golden_stream_is_independently_aggregated(self):
        result = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(), **paired_parse_kwargs())

        self.assertTrue(result.valid_for_promotion)
        self.assertEqual(result.candidate.fail, 0)
        self.assertEqual(result.identity.fail, 0)
        self.assertEqual(result.candidate.oh_mean, 1.0)
        self.assertEqual(result.identity.oh_mean, 1.0)
        self.assertAlmostEqual(
            result.candidate_log_cost_mean, math.log(0.8), places=14)
        self.assertAlmostEqual(
            result.candidate_log_cost_ci_low, math.log(0.8), places=14)
        self.assertAlmostEqual(
            result.candidate_log_cost_ci_high, math.log(0.8), places=14)
        self.assertEqual(result.aa_log_cost_mean, 0.0)
        self.assertEqual(result.aa_floor_log, 0.0)
        self.assertEqual(result.candidate.solve_mbps, 51200.0)
        self.assertEqual(result.identity.solve_mbps, 40960.0)
        self.assertEqual(len(result.replicate_receipts), 4)
        self.assertEqual(
            [receipt["replicate"]
             for receipt in result.replicate_receipts],
            list(range(4)))
        self.assertEqual(
            result.as_dict()["protocol"], peel_codec.NATIVE_PAIRED_PROTOCOL)
        self.assertEqual(
            result.as_dict()["evidence"]["final_context_sha256"],
            peel_codec._canonical_json_sha256(paired_context()))

    def test_exact_identity_candidates_never_promote_and_controls_are_sane(
            self):
        below_half = math.nextafter(0.25, 0.0) * 2
        at_half = 0.25 * 2
        above_half = math.nextafter(0.25, math.inf) * 2
        self.assertEqual(
            [
                peel_codec._nonnegative_llround(value)
                for value in (below_half, at_half, above_half)
            ],
            [0, 1, 1])
        profile = fake_profile()
        for scale in (None, profile.target_mean, profile.target_mean + 0.01):
            with self.subTest(scale=scale):
                measurement = paired_fixture_measurement(
                    profile, profile.pmf, 9, 10,
                    degree_scale=scale, candidate_ns=80)
                self.assertLess(
                    measurement.candidate_log_cost_ci_high, 0.0)
                self.assertFalse(measurement.valid_for_promotion)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "biased"):
                    peel_codec.require_paired_stock_control(measurement)
        saturated = fake_profile(2)
        saturated_alias = paired_fixture_measurement(
            saturated, saturated.pmf, 9, 10,
            degree_scale=64000.0, candidate_ns=80)
        self.assertFalse(saturated_alias.valid_for_promotion)

        proportional = [2.0 * value for value in profile.pmf]
        proportional_alias = paired_fixture_measurement(
            profile, proportional, 9, 10, candidate_ns=80)
        self.assertNotEqual(
            peel_codec._paired_pmf_sha256(proportional),
            peel_codec._paired_pmf_sha256(profile.pmf))
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(proportional, 64),
            peel_codec._native_peel_cdf_signature(profile.pmf, 64))
        self.assertFalse(proportional_alias.valid_for_promotion)

        short_profile = replace(profile, pmf=(0.25, 0.75))
        zero_padded = [0.25, 0.75] + [0.0] * 62
        zero_padding_alias = paired_fixture_measurement(
            short_profile, zero_padded, 9, 10, candidate_ns=80)
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(zero_padded, 64),
            peel_codec._native_peel_cdf_signature(short_profile.pmf, 64))
        self.assertFalse(zero_padding_alias.valid_for_promotion)

        subgrid_alias = list(profile.pmf)
        subgrid_alias[-1] += 1e-12
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(subgrid_alias, 64),
            peel_codec._native_peel_cdf_signature(profile.pmf, 64))
        self.assertFalse(paired_fixture_measurement(
            profile, subgrid_alias, 9, 10,
            candidate_ns=80).valid_for_promotion)

        reachable_family_alias = peel_codec.family(
            profile.pmf, 100, 0, 32, 100)
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(
                reachable_family_alias, 64),
            peel_codec._native_peel_cdf_signature(profile.pmf, 64))
        self.assertFalse(paired_fixture_measurement(
            profile, reachable_family_alias, 9, 10,
            candidate_ns=80).valid_for_promotion)

        k2 = fake_profile(2)
        self.assertEqual(
            peel_codec._native_peel_cdf_signature([1.0], 2),
            peel_codec._native_peel_cdf_signature([0.0, 1.0], 2))
        self.assertFalse(paired_fixture_measurement(
            k2, [0.0, 1.0], 9, 10,
            candidate_ns=80).valid_for_promotion)

        nearby_nonalias = list(profile.pmf)
        nearby_nonalias[-1] += 1e-8
        self.assertNotEqual(
            peel_codec._native_peel_cdf_signature(nearby_nonalias, 64),
            peel_codec._native_peel_cdf_signature(profile.pmf, 64))
        self.assertTrue(paired_fixture_measurement(
            profile, nearby_nonalias, 9, 10,
            candidate_ns=80).valid_for_promotion)

        biased_aa = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(
                profile=profile, candidate_pmf=profile.pmf,
                candidate_ns=100, aa_left_ns=125, aa_right_ns=100),
            **paired_parse_kwargs(
                profile=profile, candidate_pmf=profile.pmf))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "biased"):
            peel_codec.require_paired_stock_control(biased_aa)

    def test_stock_control_cross_calibrates_practical_bias_symmetrically(
            self):
        profile = fake_profile()

        def control(
                *, candidate_ns=100, identity_ns=100, aa_left_ns=100,
                aa_right_ns=100, required_margin=0.0, mutate_row=None):
            return peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    profile=profile,
                    candidate_pmf=profile.pmf,
                    candidate_ns=candidate_ns,
                    identity_ns=identity_ns,
                    aa_left_ns=aa_left_ns,
                    aa_right_ns=aa_right_ns,
                    required_margin=required_margin,
                    mutate_row=mutate_row,
                ),
                **paired_parse_kwargs(
                    profile=profile,
                    candidate_pmf=profile.pmf,
                    required_margin=required_margin,
                ),
            )

        # Location is not uncertainty.  Two precise panels with the same large
        # bias, or with opposite large biases, must never legitimize each
        # other.
        for label, kwargs in (
                ("same-positive", {
                    "candidate_ns": 200, "aa_left_ns": 200,
                }),
                ("opposite-sign", {
                    "candidate_ns": 200, "aa_left_ns": 100,
                    "aa_right_ns": 200,
                }),
                ("same-negative", {
                    "candidate_ns": 100, "identity_ns": 200,
                    "aa_left_ns": 100, "aa_right_ns": 200,
                })):
            with self.subTest(precise_confounded=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "biased"):
                    peel_codec.require_paired_stock_control(control(**kwargs))

        # Comparable measured noise in either direction can cover a small
        # off-zero interval without requiring both 95% intervals to hit zero.
        def noisy_control(candidate_ratios, aa_ratios):
            denominator = 100000

            def vary_replicates(unused_index, row):
                replicate = row["replicate"]
                if row["role"] == "candidate":
                    row["elapsed_ns"] = round(
                        denominator * candidate_ratios[replicate])
                elif row["role"] == "identity":
                    row["elapsed_ns"] = denominator
                elif row["role"] == "identity_a":
                    row["elapsed_ns"] = round(
                        denominator * aa_ratios[replicate])
                elif row["role"] == "identity_b":
                    row["elapsed_ns"] = denominator

            return control(
                candidate_ns=denominator,
                identity_ns=denominator,
                aa_left_ns=denominator,
                aa_right_ns=denominator,
                mutate_row=vary_replicates,
            )

        positive_ratios = [1.01, 1.015, 1.025, 1.03]
        negative_ratios = [0.99, 0.985, 0.975, 0.97]
        for label, measurement in (
                ("same-positive", noisy_control(
                    positive_ratios, positive_ratios)),
                ("opposite-sign", noisy_control(
                    positive_ratios, negative_ratios))):
            with self.subTest(noisy_cross_calibration=label):
                self.assertFalse(
                    measurement.candidate_log_cost_ci_low <= 0.0 <=
                    measurement.candidate_log_cost_ci_high)
                self.assertFalse(
                    measurement.aa_log_cost_ci_low <= 0.0 <=
                    measurement.aa_log_cost_ci_high)
                peel_codec.require_paired_stock_control(measurement)

        # Exact equality with the independent half-width is acceptable,
        # matching the strict inequality used by the promotion decision.
        equality = replace(
            control(),
            candidate_log_cost_mean=0.01,
            candidate_log_cost_ci_low=0.01,
            candidate_log_cost_ci_high=0.01,
            aa_log_cost_mean=0.0,
            aa_log_cost_ci_low=-0.01,
            aa_log_cost_ci_high=0.01,
        )
        self.assertEqual(
            equality.candidate_log_cost_ci_low,
            (equality.aa_log_cost_ci_high -
             equality.aa_log_cost_ci_low) / 2.0,
        )
        peel_codec.require_paired_stock_control(equality)

        # With no declared practical margin and a perfectly quiet independent
        # panel, even a small wholly off-zero interval remains relevant.
        for label, kwargs in (
                ("candidate-positive", {"candidate_ns": 101}),
                ("candidate-negative", {"candidate_ns": 99}),
                ("aa-positive", {"aa_left_ns": 101}),
                ("aa-negative", {"aa_left_ns": 99})):
            with self.subTest(zero_floor=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "biased"):
                    peel_codec.require_paired_stock_control(control(**kwargs))

        # A predeclared practical margin accepts sub-margin stock drift even
        # when the other panel is exactly zero.
        for label, kwargs in (
                ("candidate", {"candidate_ns": 101}),
                ("aa", {"aa_left_ns": 101})):
            with self.subTest(practical_margin=label):
                peel_codec.require_paired_stock_control(control(
                    required_margin=0.02, **kwargs))

        # The exact decimal margin boundaries use the same floating-point
        # ratios as the measured log costs.  Equality is accepted for stock
        # controls in both directions and in both panels.
        for label, kwargs in (
                ("candidate-positive", {
                    "candidate_ns": 101, "identity_ns": 100,
                }),
                ("candidate-negative", {
                    "candidate_ns": 100, "identity_ns": 101,
                }),
                ("aa-positive", {
                    "aa_left_ns": 101, "aa_right_ns": 100,
                }),
                ("aa-negative", {
                    "aa_left_ns": 100, "aa_right_ns": 101,
                })):
            with self.subTest(practical_margin_equality=label):
                peel_codec.require_paired_stock_control(control(
                    required_margin=0.01, **kwargs))

        exact_margin = paired_fixture_measurement(
            profile, (0.25, 0.75), 9, 10,
            candidate_ns=100, identity_ns=101, required_margin=0.01)
        required_lower, unused_required_upper = \
            peel_codec._paired_practical_log_bounds(0.01)
        self.assertEqual(
            exact_margin.candidate_log_cost_ci_high, required_lower)
        self.assertFalse(exact_margin.valid_for_promotion)
        self.assertTrue(paired_fixture_measurement(
            profile, (0.25, 0.75), 9, 10,
            candidate_ns=99, identity_ns=100,
            required_margin=0.01).valid_for_promotion)

        # Public timing validation accepts integer margins as well as floats;
        # parsing must preserve that contract all the way through aggregation.
        for integer_margin in (0, 1):
            with self.subTest(integer_margin=integer_margin):
                measurement = peel_codec.parse_peeltiming_output(
                    peeltiming_stdout(
                        profile=profile,
                        candidate_pmf=profile.pmf,
                        required_margin=integer_margin,
                    ),
                    **paired_parse_kwargs(
                        profile=profile,
                        candidate_pmf=profile.pmf,
                        required_margin=integer_margin,
                    ),
                )
                self.assertEqual(
                    measurement.required_margin, float(integer_margin))

        # Either panel, in either direction, must still veto when it is
        # materially farther from zero than its independent envelope.
        for label, kwargs in (
                ("candidate-positive", {
                    "candidate_ns": 125, "aa_left_ns": 101,
                }),
                ("candidate-negative", {
                    "candidate_ns": 80, "aa_left_ns": 101,
                }),
                ("aa-positive", {
                    "candidate_ns": 101, "aa_left_ns": 125,
                }),
                ("aa-negative", {
                    "candidate_ns": 101, "aa_left_ns": 80,
                })):
            with self.subTest(material_bias=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "biased"):
                    peel_codec.require_paired_stock_control(control(**kwargs))

    def test_peeltiming_supports_every_bound_packet_schedule(self):
        for schedule in sorted(peel_codec.TARGET_SCHEDULES):
            with self.subTest(schedule=schedule):
                result = peel_codec.parse_peeltiming_output(
                    peeltiming_stdout(schedule=schedule),
                    **paired_parse_kwargs(schedule=schedule))
                self.assertEqual(result.manifest["schedule"], schedule)
                self.assertEqual(
                    result.manifest["loss_model"],
                    peel_codec.PAIRED_LOSS_MODEL)
                self.assertTrue(result.valid_for_promotion)

    def test_compare_probe_supports_every_bound_packet_schedule(self):
        with mock.patch("peel_codec.subprocess.run") as run:
            for schedule in sorted(peel_codec.TARGET_SCHEDULES):
                with self.subTest(schedule=schedule):
                    run.return_value = subprocess.CompletedProcess(
                        ["bench"], 0,
                        stdout=compare_stdout(schedule=schedule), stderr="")
                    result = peel_codec.compare_probe(
                        "unused", 64, 10, 64,
                        **probe_kwargs(schedule=schedule))
                    self.assertEqual(
                        result.target_receipt["schedule"], schedule)
                    self.assertEqual(
                        result.target_receipt["loss"],
                        "packet-schedule-v1")
                    self.assertEqual(
                        run.call_args.args[0][
                            run.call_args.args[0].index("--schedule") + 1],
                        schedule)

    def test_peeltiming_rejects_authenticated_balance_corruption(self):
        trace_zero = hashlib.sha256(b"trace-0").hexdigest()
        first_construction = (
            9 ^ (0xd1b54a32d192ed03 & 0xffffffffffffffff)
        )

        def changed_trace(index, row):
            if index == 1:
                row["trace_sha256"] = hashlib.sha256(
                    b"changed-trace").hexdigest()

        def reused_trace(index, row):
            if 16 <= index < 32:
                row["trace_sha256"] = trace_zero

        def reused_construction_seed(index, row):
            if 16 <= index < 32:
                row["construction_seed"] = first_construction

        def wrong_loss_seed(index, row):
            if index == 0:
                row["loss_seed"] ^= 1

        def wrong_common_overhead(unused_index, row):
            if row["replicate"] == 0:
                row["common_overhead"] = 2
                row["packet_payload_bytes"] = (64 + 2) * 64

        def negative_counter(index, row):
            if index == 0:
                row["block_xors"] = -1

        cases = {
            "slot": lambda index, row: (
                row.__setitem__("slot", 7) if index == 0 else None),
            "label": lambda index, row: (
                row.__setitem__("label", "B") if index == 0 else None),
            "role": lambda index, row: (
                row.__setitem__("role", "candidate") if index == 0 else None),
            "changed trace": changed_trace,
            "reused trace": reused_trace,
            "reused construction seed": reused_construction_seed,
            "wrong loss seed": wrong_loss_seed,
            "wrong common overhead": wrong_common_overhead,
            "negative counter": negative_counter,
        }
        for name, mutation in cases.items():
            with self.subTest(corruption=name):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(mutate_row=mutation),
                        **paired_parse_kwargs())

    def test_peeltiming_rejects_contaminated_timing_rows(self):
        def migrated(index, row):
            if index == 0:
                row["cpu_after"] = 3
                row["cpu_migrated"] = 1

        def faulted(index, row):
            if index == 0:
                row["minflt_delta"] = 1
                row["fault_contaminated"] = 1

        def saturated(index, row):
            if index == 0:
                row["saturated"] = 1

        for name, mutation in {
                "migration": migrated,
                "fault": faulted,
                "saturation": saturated,
        }.items():
            with self.subTest(contamination=name):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "balance contract"):
                    peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(mutate_row=mutation),
                        **paired_parse_kwargs())

    def test_peeltiming_rejects_truncation_hash_and_record_reordering(self):
        stdout = peeltiming_stdout()
        with self.assertRaises(peel_codec.MeasurementError):
            peel_codec.parse_peeltiming_output(
                stdout[:-1], **paired_parse_kwargs())
        prefix, marker, unused_digest = stdout.rpartition("stream_sha256=")
        self.assertEqual(marker, "stream_sha256=")
        bad_hash = prefix + marker + ("0" * 64) + "\n"
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "SHA-256"):
            peel_codec.parse_peeltiming_output(
                bad_hash, **paired_parse_kwargs())
        changed_finish = stdout.replace(
            "finished_monotonic_ns=2000000000",
            "finished_monotonic_ns=2000000001")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "SHA-256"):
            peel_codec.parse_peeltiming_output(
                changed_finish, **paired_parse_kwargs())

        def reorder(lines):
            lines[3], lines[4] = lines[4], lines[3]

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "balance contract"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, reorder),
                **paired_parse_kwargs())

    def test_peeltiming_semantic_and_header_receipts_are_exact(self):
        stdout = peeltiming_stdout()

        def fail_semantic(lines):
            lines[1] = lines[1].replace(
                "solve_equal=1", "solve_equal=0")

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "semantic identity"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, fail_semantic),
                **paired_parse_kwargs())

        def forge_need_more(lines):
            lines[1] = (
                lines[1]
                .replace("nohook_result=0", "nohook_result=1")
                .replace(
                    "nohook_recovery_ok=1", "nohook_recovery_ok=0")
                .replace("identity_result=0", "identity_result=1")
                .replace(
                    "identity_recovery_ok=1", "identity_recovery_ok=0")
            )

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "semantic identity"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, forge_need_more),
                **paired_parse_kwargs())

        def forge_semantic_trace(lines):
            original = hashlib.sha256(b"semantic-trace").hexdigest()
            forged = hashlib.sha256(b"trace-0").hexdigest()
            lines[1] = lines[1].replace(
                f"trace_sha256={original}", f"trace_sha256={forged}")

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "semantic trace"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, forge_semantic_trace),
                **paired_parse_kwargs())

        def forge_identity_scale(lines):
            lines[0] = lines[0].replace(
                "candidate_scale_effective=12.8",
                "candidate_scale_effective=12.9").replace(
                    "identity_scale_effective=12.8",
                    "identity_scale_effective=12.9")
            lines[1] = lines[1].replace(
                "identity_scale_effective=12.8",
                "identity_scale_effective=12.9")

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "identity scale"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, forge_identity_scale),
                **paired_parse_kwargs())

        def reorder_header(lines):
            columns = lines[2].split(",")
            columns[0], columns[1] = columns[1], columns[0]
            lines[2] = ",".join(columns)

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "header"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, reorder_header),
                **paired_parse_kwargs())

    def test_peeltiming_aa_recovery_and_overhead_gates_fail_closed(self):
        aa_biased = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(aa_left_ns=125, aa_right_ns=100),
            **paired_parse_kwargs())
        self.assertFalse(aa_biased.valid_for_promotion)
        self.assertGreater(aa_biased.aa_log_cost_ci_low, 0.0)

        def recovery_regression(unused_index, row):
            if row["replicate"] != 5:
                return
            skip_peeltiming_row(row, common_overhead=16)
            if row["role"] == "candidate":
                row["recovery_result"] = 1
                row["recovery_ok"] = 0
                row["arm_overhead"] = -1

        recovery = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(
                replicates=6, mutate_row=recovery_regression),
            **paired_parse_kwargs(replicates=6))
        self.assertTrue(recovery.timing_ci_available)
        self.assertFalse(recovery.valid_for_promotion)
        self.assertEqual(recovery.recovery_regressions, 1)
        self.assertEqual(recovery.candidate.fail, 1)

        def overhead_regression(unused_index, row):
            if row["replicate"] == 0:
                row["common_overhead"] = 2
                row["packet_payload_bytes"] = (64 + 2) * 64
                if row["role"] == "candidate":
                    row["arm_overhead"] = 2

        overhead = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(mutate_row=overhead_regression),
            **paired_parse_kwargs())
        self.assertTrue(overhead.valid_for_promotion)
        self.assertEqual(overhead.overhead_regressions, 1)

    def test_peeltiming_records_weak_seed_without_fabricating_timing(self):
        def weak_final_replicate(unused_index, row):
            if row["replicate"] != 3:
                return
            skip_peeltiming_row(row, common_overhead=16)
            row.update({
                "arm_overhead": -1,
                "recovery_result": 4,
                "recovery_ok": 0,
            })

        measurement = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(mutate_row=weak_final_replicate),
            **paired_parse_kwargs())
        self.assertEqual(measurement.timing_eligible_replicates, 3)
        self.assertFalse(measurement.timing_ci_available)
        self.assertFalse(measurement.valid_for_promotion)
        self.assertEqual(
            measurement.candidate_recovery_result_counts["bad_peel_seed"], 1)
        self.assertEqual(
            measurement.identity_recovery_result_counts["bad_peel_seed"], 1)
        self.assertEqual(measurement.candidate_weak_seed_failures, 1)
        self.assertEqual(measurement.identity_weak_seed_failures, 1)
        self.assertEqual(measurement.both_failures, 1)
        skipped = measurement.replicate_receipts[-1]
        self.assertFalse(skipped["timing_eligible"])
        self.assertIsNone(skipped["candidate_log_cost"])
        self.assertIsNone(skipped["aa_log_cost"])

    def test_shared_raw_weak_seed_is_descriptive_not_a_promotion_veto(self):
        profile = fake_profile()

        def shared_weak(unused_index, row):
            if row["replicate"] != 5:
                return
            skip_peeltiming_row(row, common_overhead=16)
            row.update({
                "arm_overhead": -1,
                "recovery_result": 4,
                "recovery_ok": 0,
            })

        measurement = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(
                profile=profile, replicates=6, mutate_row=shared_weak),
            **paired_parse_kwargs(profile=profile, replicates=6))
        self.assertEqual(measurement.timing_eligible_replicates, 5)
        self.assertTrue(measurement.timing_ci_available)
        self.assertTrue(measurement.valid_for_promotion)
        self.assertEqual(measurement.candidate.fail, 1)
        self.assertEqual(measurement.identity.fail, 1)
        self.assertEqual(measurement.both_failures, 1)
        self.assertGreater(measurement.candidate.goodput(64), 0.0)
        self.assertGreater(measurement.identity.goodput(64), 0.0)
        self.assertEqual(
            peel_codec.paired_arm_goodput(measurement.candidate, 64),
            peel_codec.paired_arm_goodput(
                measurement.candidate.as_dict(), 64))

        stock_control = peel_codec.parse_peeltiming_output(
            peeltiming_stdout(
                profile=profile, candidate_pmf=profile.pmf,
                candidate_ns=100, replicates=6, mutate_row=shared_weak),
            **paired_parse_kwargs(
                profile=profile, candidate_pmf=profile.pmf,
                replicates=6))
        peel_codec.require_paired_stock_control(stock_control)

    def test_peeltiming_rejects_fatal_recovery_result(self):
        def fatal_result(index, row):
            if index == 0:
                row["recovery_result"] = 8
                row["recovery_ok"] = 0
                row["arm_overhead"] = -1

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "balance contract"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(mutate_row=fatal_result),
                **paired_parse_kwargs())

    def test_peeltiming_thermal_context_must_exist_and_cover_run(self):
        invalid_device = paired_context()
        invalid_device["bound"]["thermal_device"] = -1
        invalid_contexts = [
            {
                "schema": peel_codec.PAIRED_CONTEXT_SCHEMA,
                "bound": paired_context()["bound"],
                "thermal": None,
            },
            paired_context(start=1.000000001),
            paired_context(end=1.999999999),
            invalid_device,
        ]
        for context in invalid_contexts:
            with self.subTest(context=context):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(context=context),
                        **paired_parse_kwargs(context=context))

        original = paired_context()
        tampered = paired_context()
        tampered["bound"]["thermal_inode"] += 1
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "manifest"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(context=original),
                **paired_parse_kwargs(context=tampered))

    def test_thermal_window_validates_only_selected_raw_bracket(self):
        payload = thermal_csv(
            # A sampler restart/gap long before this run is irrelevant.
            thermal_row(1.0),
            thermal_row(10.0),
            thermal_row(100.0),
            thermal_row(
                101.0, dimm_read_errors="1", missing_dimm=True),
            thermal_row(102.0),
            thermal_row(200.0),
        )
        summary = peel_codec._thermal_window_from_csv(
            payload, 100_200_000_000, 101_800_000_000)
        expected_artifact = thermal_csv(
            thermal_row(100.0),
            thermal_row(
                101.0, dimm_read_errors="1", missing_dimm=True),
            thermal_row(102.0),
        )
        self.assertEqual(
            summary["csv_sha256"],
            hashlib.sha256(expected_artifact).hexdigest())
        self.assertEqual(summary["rows"], 3)
        self.assertEqual(summary["valid_rows"], 2)
        self.assertEqual(summary["missing_read_rows"], 1)
        self.assertEqual(summary["invalid_rows"], 1)
        self.assertEqual(summary["first_monotonic_s"], 100.0)
        self.assertEqual(summary["last_monotonic_s"], 102.0)
        self.assertEqual(
            summary["summary_sha256"],
            peel_codec._canonical_json_sha256({
                name: value for name, value in summary.items()
                if name != "summary_sha256"
            }))

    def test_thermal_window_rejects_noncanonical_edac_and_partial_rows(self):
        cases = {
            "leading-zero counter": thermal_csv(
                thermal_row(100.0),
                thermal_row(101.0, dimm_read_errors="01"),
                thermal_row(102.0)),
            "scientific number": thermal_csv(
                thermal_row(100.0),
                thermal_row("1.01e2"),
                thermal_row(102.0)),
            "EDAC": thermal_csv(
                thermal_row(100.0),
                thermal_row(101.0, edac_ce="1"),
                thermal_row(102.0)),
            "partial append": (
                thermal_csv(thermal_row(100.0)) +
                thermal_row(102.0).encode("ascii")[:-1]
            ),
        }
        for name, payload in cases.items():
            with self.subTest(corruption=name):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec._thermal_window_from_csv(
                        payload, 100_200_000_000, 101_800_000_000)

    def test_paired_context_summary_hash_is_recomputed(self):
        context = paired_context()
        context["thermal"]["cpu_tctl_max_c"] += 1.0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "thermal context"):
            peel_codec._validate_paired_context(
                context, cache_state="warm", evict_bytes=4096,
                started_ns=1_000_000_000, finished_ns=2_000_000_000)

    @mock.patch("peel_codec._validate_runtime_identity")
    def test_thermal_capture_rejects_prefill_stale_interval_and_replacement(
            self, unused_identity):
        now = time.clock_gettime(time.CLOCK_MONOTONIC)

        def config(path, interval=1000):
            base = paired_context()["bound"]
            for name in peel_codec._PAIRED_BOUND_CAPTURE_FIELDS:
                del base[name]
            base["thermal_source"] = path
            base["thermal_sampling_interval_ms"] = interval
            return {
                "schema": peel_codec.PAIRED_CONTEXT_SCHEMA,
                "bound": base,
                "thermal": None,
            }

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "thermal.csv")
            with open(path, "wb") as output:
                output.write(thermal_csv(
                    thermal_row(now - 1.0),
                    thermal_row(now)))

            prefilled = config(path)
            prefilled["thermal"] = paired_context()["thermal"]
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "unfinalized"):
                peel_codec._prepare_paired_context(
                    prefilled, cache_state="warm", evict_bytes=4096)

            stale_path = os.path.join(directory, "stale.csv")
            with open(stale_path, "wb") as output:
                output.write(thermal_csv(
                    thermal_row(now - 11.0),
                    thermal_row(now - 10.0)))
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "stale"):
                peel_codec._prepare_paired_context(
                    config(stale_path),
                    cache_state="warm", evict_bytes=4096)

            mismatch_path = os.path.join(directory, "mismatch.csv")
            with open(mismatch_path, "wb") as output:
                output.write(thermal_csv(
                    thermal_row(now - 1.0),
                    thermal_row(now)))
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "interval"):
                peel_codec._prepare_paired_context(
                    config(mismatch_path, interval=10),
                    cache_state="warm", evict_bytes=4096)

            frozen, capture = peel_codec._prepare_paired_context(
                config(path), cache_state="warm", evict_bytes=4096)
            try:
                self.assertIsNone(frozen["thermal"])
                self.assertEqual(
                    frozen["bound"]["thermal_prelaunch_row_sha256"],
                    hashlib.sha256(thermal_row(now).encode("ascii")).hexdigest())
                replacement = os.path.join(directory, "replacement.csv")
                with open(replacement, "wb") as output:
                    output.write(thermal_csv(
                        thermal_row(now), thermal_row(now + 1.0)))
                real_pread = os.pread
                replaced_path = False

                def replace_during_read(fd, count, offset):
                    nonlocal replaced_path
                    payload = real_pread(fd, count, offset)
                    if not replaced_path:
                        os.replace(replacement, path)
                        replaced_path = True
                    return payload

                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "identity changed"):
                    with mock.patch(
                            "peel_codec.os.pread",
                            side_effect=replace_during_read):
                        peel_codec._snapshot_capture(capture)
            finally:
                os.close(capture.fd)

    @mock.patch("peel_codec._validate_runtime_identity")
    def test_thermal_capture_retries_a_partial_sampler_append(
            self, unused_identity):
        now = time.clock_gettime(time.CLOCK_MONOTONIC)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "thermal.csv")
            next_row = thermal_row(now).encode("ascii")
            split = len(next_row) // 2
            with open(path, "wb") as output:
                output.write(
                    thermal_csv(thermal_row(now - 1.0)) +
                    next_row[:split])

            context = paired_context()
            for name in peel_codec._PAIRED_BOUND_CAPTURE_FIELDS:
                del context["bound"][name]
            context["bound"]["thermal_source"] = path
            context["bound"]["thermal_sampling_interval_ms"] = 1000
            context["thermal"] = None

            def finish_append(unused_seconds):
                with open(path, "ab") as output:
                    output.write(next_row[split:])

            with mock.patch(
                    "peel_codec.time.sleep",
                    side_effect=finish_append) as sleep:
                frozen, capture = peel_codec._prepare_paired_context(
                    context, cache_state="warm", evict_bytes=4096)
            try:
                sleep.assert_called_once_with(0.1)
                self.assertEqual(
                    frozen["bound"]["thermal_prelaunch_row_sha256"],
                    hashlib.sha256(next_row).hexdigest())
                self.assertEqual(capture.prefix, thermal_csv(
                    thermal_row(now - 1.0), thermal_row(now)))
            finally:
                os.close(capture.fd)

    @mock.patch("peel_codec._validate_runtime_identity")
    def test_thermal_capture_is_bounded_and_authenticates_retained_tail(
            self, unused_identity):
        now = time.clock_gettime(time.CLOCK_MONOTONIC)
        header = (
            ",".join(peel_codec._THERMAL_CSV_COLUMNS) + "\n"
        ).encode("ascii")
        old_row = thermal_row(now - 100.0).encode("ascii")
        repeat_count = (4 * 1024 * 1024) // len(old_row) + 1
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "thermal.csv")
            with open(path, "wb") as output:
                output.write(header)
                output.write(old_row * repeat_count)
                output.write(thermal_row(now - 1.0).encode("ascii"))
                output.write(thermal_row(now).encode("ascii"))

            context = paired_context()
            for name in peel_codec._PAIRED_BOUND_CAPTURE_FIELDS:
                del context["bound"][name]
            context["bound"]["thermal_source"] = path
            context["bound"]["thermal_sampling_interval_ms"] = 1000
            context["thermal"] = None
            real_pread = os.pread
            requests = []

            def record_pread(fd, count, offset):
                requests.append((count, offset))
                return real_pread(fd, count, offset)

            with mock.patch(
                    "peel_codec.os.pread",
                    side_effect=record_pread):
                unused_frozen, capture = peel_codec._prepare_paired_context(
                    context, cache_state="warm", evict_bytes=4096)
                try:
                    self.assertGreater(capture.tail_offset, len(header))
                    self.assertLessEqual(
                        len(capture.prefix),
                        len(header) +
                        peel_codec._THERMAL_CAPTURE_TAIL_BYTES)
                    appended = thermal_row(now + 1.0).encode("ascii")
                    with open(path, "ab") as output:
                        output.write(appended)
                    snapshot = peel_codec._snapshot_capture(capture)
                    summary = peel_codec._thermal_window_from_csv(
                        snapshot,
                        int((now + 0.1) * 1e9),
                        int((now + 0.9) * 1e9),
                        capture,
                        sampling_interval_ms=1000)
                    self.assertEqual(summary["rows"], 2)
                    self.assertTrue(requests)
                    self.assertLessEqual(
                        max(count for count, unused_offset in requests),
                        peel_codec._THERMAL_CAPTURE_TAIL_BYTES +
                        len(appended))

                    relative = capture.prefix.find(capture.prelaunch_line)
                    self.assertGreaterEqual(relative, len(header))
                    absolute = (
                        capture.tail_offset + relative - len(header))
                    tampered = capture.prelaunch_line.replace(
                        b",60.0,", b",61.0,", 1)
                    self.assertEqual(len(tampered), len(capture.prelaunch_line))
                    with open(path, "r+b") as output:
                        output.seek(absolute)
                        output.write(tampered)
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError, "modified"):
                        peel_codec._snapshot_capture(capture)

                    with open(path, "r+b") as output:
                        output.seek(absolute)
                        output.write(capture.prelaunch_line)
                    os.truncate(path, capture.initial_size - 1)
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError, "truncated"):
                        peel_codec._snapshot_capture(capture)
                finally:
                    os.close(capture.fd)

    @mock.patch("peel_codec._validate_runtime_identity")
    def test_prepare_context_partial_retry_is_narrow_and_bounded(
            self, unused_identity):
        now = time.clock_gettime(time.CLOCK_MONOTONIC)
        complete = thermal_csv(
            thermal_row(now - 1.0), thermal_row(now))
        header = (
            ",".join(peel_codec._THERMAL_CSV_COLUMNS) + "\n"
        ).encode("ascii")

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "thermal.csv")
            with open(path, "wb") as output:
                output.write(complete)

            def config():
                context = paired_context()
                for name in peel_codec._PAIRED_BOUND_CAPTURE_FIELDS:
                    del context["bound"][name]
                context["bound"]["thermal_source"] = path
                context["bound"]["thermal_sampling_interval_ms"] = 1000
                context["thermal"] = None
                return context

            for label, partial in (
                    ("tail", complete[:-1]),
                    ("header", header[:len(header) // 2])):
                with self.subTest(success=label):
                    with mock.patch(
                            "peel_codec._snapshot_capture",
                            side_effect=[partial, complete]) as snapshot, \
                            mock.patch("peel_codec.time.sleep") as sleep:
                        frozen, capture = peel_codec._prepare_paired_context(
                            config(), cache_state="warm", evict_bytes=4096)
                    try:
                        self.assertIsNone(frozen["thermal"])
                        self.assertEqual(snapshot.call_count, 2)
                        self.assertIs(
                            snapshot.call_args_list[0].args[0],
                            snapshot.call_args_list[1].args[0])
                        sleep.assert_called_once_with(0.1)
                    finally:
                        os.close(capture.fd)

            malformed = (
                thermal_csv(thermal_row(now - 1.0)) +
                b"malformed,row\n")
            for label, payload in (
                    ("divergent header", b"not-a-header"),
                    ("complete malformed row", malformed)):
                with self.subTest(immediate_failure=label):
                    with mock.patch(
                            "peel_codec._snapshot_capture",
                            return_value=payload) as snapshot, \
                            mock.patch("peel_codec.time.sleep") as sleep:
                        with self.assertRaises(peel_codec.MeasurementError):
                            peel_codec._prepare_paired_context(
                                config(), cache_state="warm",
                                evict_bytes=4096)
                    self.assertEqual(snapshot.call_count, 1)
                    sleep.assert_not_called()

            with mock.patch(
                    "peel_codec._snapshot_capture",
                    return_value=complete[:-1]) as snapshot, \
                    mock.patch(
                        "peel_codec.time.monotonic",
                        side_effect=[0.0, 0.0, 10.0]), \
                    mock.patch("peel_codec.time.sleep") as sleep:
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "tail remained incomplete"):
                    peel_codec._prepare_paired_context(
                        config(), cache_state="warm", evict_bytes=4096)
            self.assertEqual(snapshot.call_count, 2)
            sleep.assert_called_once_with(0.1)

    def test_peeltiming_rejects_cpu_opstat_and_stock_control_drift(self):
        def off_affinity(unused_index, row):
            row["cpu_before"] = row["cpu_after"] = 999

        def impossible_rank(unused_index, row):
            row["heavy_gain"] = row["binary_def"] + 1

        for label, mutation in (
                ("CPU affinity", off_affinity),
                ("rank chain", impossible_rank)):
            with self.subTest(corruption=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "balance contract"):
                    peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(mutate_row=mutation),
                        **paired_parse_kwargs())

        multi_cpu_context = paired_context()
        multi_cpu_context["bound"]["cpu_affinity"] = [2, 3]
        multi_cpu_context["bound"]["cpu_governors"]["3"] = "performance"

        def cross_cpu_pair(unused_index, row):
            if row["role"] == "candidate" and row["pair"] == 0:
                row["cpu_before"] = row["cpu_after"] = 3

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "balance contract"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    context=multi_cpu_context,
                    mutate_row=cross_cpu_pair),
                **paired_parse_kwargs(context=multi_cpu_context))

        def second_allowed_cpu(unused_index, row):
            row["cpu_before"] = row["cpu_after"] = 3

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "balance contract"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    context=multi_cpu_context,
                    mutate_row=second_allowed_cpu),
                **paired_parse_kwargs(context=multi_cpu_context))

        profile = fake_profile()

        def drift_stock(unused_index, row):
            if row["role"] == "candidate":
                row["block_xors"] += 1

        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "stock candidate/control metadata"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    profile=profile,
                    candidate_pmf=profile.pmf,
                    candidate_ns=100,
                    mutate_row=drift_stock),
                **paired_parse_kwargs(
                    profile=profile, candidate_pmf=profile.pmf))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "stock candidate/control metadata"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    profile=profile,
                    candidate_pmf=profile.pmf,
                    candidate_ns=100,
                    degree_scale=12.81,
                    mutate_row=drift_stock),
                **paired_parse_kwargs(
                    profile=profile, candidate_pmf=profile.pmf,
                    degree_scale=12.81))
        proportional = [2.0 * value for value in profile.pmf]
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "stock candidate/control metadata"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    profile=profile,
                    candidate_pmf=proportional,
                    candidate_ns=100,
                    mutate_row=drift_stock),
                **paired_parse_kwargs(
                    profile=profile, candidate_pmf=proportional))
        k2 = fake_profile(2)
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "stock candidate/control metadata"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(
                    profile=k2,
                    candidate_pmf=[0.0, 1.0],
                    candidate_ns=100,
                    mutate_row=drift_stock),
                **paired_parse_kwargs(
                    profile=k2, candidate_pmf=[0.0, 1.0]))

    def test_peeltiming_rejects_native_integer_overflow_and_width_drift(self):
        def overflow(name, value):
            def mutate(unused_index, row):
                row[name] = value
            return mutate

        cases = (
            ("elapsed uint64", overflow("elapsed_ns", 1 << 64),
             "integer domain"),
            ("op uint64", overflow("block_xors", 1 << 64),
             "integer domain"),
            ("CPU int", overflow("cpu_before", 1 << 31),
             "integer domain"),
            (
                "system width",
                overflow(
                    "inactivated",
                    64 + fake_profile().staircase +
                    fake_profile().dense_rows +
                    fake_profile().heavy_rows + 1),
                "balance contract",
            ),
            ("success rank gap", overflow("heavy_gain", 0),
             "balance contract"),
        )
        for label, mutation, message in cases:
            with self.subTest(corruption=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, message):
                    peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(mutate_row=mutation),
                        **paired_parse_kwargs())

        def too_many_heavy_rows(unused_index, row):
            value = fake_profile().heavy_rows + 1
            row["heavy_gain"] = value
            row["binary_def"] = value
            row["inactivated"] = value

        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "balance contract"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(mutate_row=too_many_heavy_rows),
                **paired_parse_kwargs())

    def test_peeltiming_rejects_signed_zero_policy_and_capture_values(self):
        with self.assertRaisesRegex(ValueError, "dimensions"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(required_margin=-0.0),
                **paired_parse_kwargs(required_margin=-0.0))
        context = paired_context()
        context["bound"]["thermal_prelaunch_monotonic_s"] = -0.0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "capture context"):
            peel_codec._validate_paired_bound_context(
                context, cache_state="warm", evict_bytes=4096,
                require_capture=True)
        for label, mutate in (
                (
                    "CPU domain",
                    lambda bound: (
                        bound.__setitem__("cpu_affinity", [1 << 100]),
                        bound.__setitem__(
                            "cpu_governors", {str(1 << 100): "performance"}),
                    ),
                ),
                (
                    "governor whitespace",
                    lambda bound: bound["cpu_governors"].__setitem__(
                        "2", "performance bad"),
                ),
                (
                    "relative thermal path",
                    lambda bound: bound.__setitem__(
                        "thermal_source", "relative.csv"),
                ),
                (
                    "double slash thermal path",
                    lambda bound: bound.__setitem__(
                        "thermal_source", "//tmp/thermal.csv"),
                ),
                (
                    "NUL thermal path",
                    lambda bound: bound.__setitem__(
                        "thermal_source", "/tmp/test\0-sampler.csv"),
                ),
                (
                    "device overflow",
                    lambda bound: bound.__setitem__(
                        "thermal_device", 1 << 64),
                ),
                (
                    "inode overflow",
                    lambda bound: bound.__setitem__(
                        "thermal_inode", 1 << 64),
                ),
        ):
            forged = paired_context()
            mutate(forged["bound"])
            with self.subTest(context_alias=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "(?:bound|capture) context"):
                    peel_codec._validate_paired_bound_context(
                        forged, cache_state="warm", evict_bytes=4096,
                        require_capture=True)
        boundary = paired_context()
        boundary["bound"]["thermal_device"] = 0xffffffffffffffff
        boundary["bound"]["thermal_inode"] = 0xffffffffffffffff
        peel_codec._validate_paired_bound_context(
            boundary, cache_state="warm", evict_bytes=4096,
            require_capture=True)

    def test_peeltiming_rejects_clock_overflow_and_scale_type_aliases(self):
        stdout = peeltiming_stdout(degree_scale=1.0)

        def started_overflow(lines):
            fields = lines[0].split(",")
            index = next(
                i for i, field in enumerate(fields)
                if field.startswith("started_monotonic_ns="))
            fields[index] = f"started_monotonic_ns={1 << 64}"
            lines[0] = ",".join(fields)

        def finished_overflow(lines):
            lines[-1] = lines[-1].replace(
                "finished_monotonic_ns=2000000000",
                f"finished_monotonic_ns={1 << 64}",
            )

        def equal_finish(lines):
            lines[-1] = lines[-1].replace(
                "finished_monotonic_ns=2000000000",
                "finished_monotonic_ns=1000000000",
            )

        for label, edit in (
                ("started", started_overflow),
                ("finished", finished_overflow)):
            with self.subTest(clock=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "integer domain"):
                    peel_codec.parse_peeltiming_output(
                        authenticated_peeltiming_edit(stdout, edit),
                        **paired_parse_kwargs(degree_scale=1.0))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "clock did not advance"):
            peel_codec.parse_peeltiming_output(
                authenticated_peeltiming_edit(stdout, equal_finish),
                **paired_parse_kwargs(degree_scale=1.0))

        for alias in (True, 1, "1"):
            with self.subTest(scale_alias=alias):
                with self.assertRaisesRegex(ValueError, "canonical float"):
                    peel_codec.parse_peeltiming_output(
                        stdout,
                        **paired_parse_kwargs(degree_scale=alias))
                with mock.patch("peel_codec.subprocess.run") as run:
                    with self.assertRaisesRegex(
                            ValueError, "canonical float"):
                        peel_codec.paired_probe(
                            "unused", 64, 64, [0.25, 0.75],
                            degree_scale=alias,
                            warmup_replicates=0,
                            replicates=4,
                            inner_reps=1,
                            max_overhead=16,
                            cache_state="warm",
                            evict_bytes=4096,
                            context=paired_context(),
                            required_margin=0.0,
                            **probe_kwargs())
                    run.assert_not_called()

    def test_thermal_interval_domain_and_counter_overflow_are_exact(self):
        bound_context = paired_context()
        bound_context["bound"]["thermal_sampling_interval_ms"] = 2000
        peel_codec._validate_paired_bound_context(
            bound_context, cache_state="warm", evict_bytes=4096,
            require_capture=True)
        bound_context["bound"]["thermal_sampling_interval_ms"] = 2001
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "bound context"):
            peel_codec._validate_paired_bound_context(
                bound_context, cache_state="warm", evict_bytes=4096,
                require_capture=True)

        jittered = thermal_csv(
            thermal_row(0.0),
            thermal_row(2.001),
            thermal_row(4.002),
        )
        summary = peel_codec._thermal_window_from_csv(
            jittered, 100_000_000, 3_900_000_000,
            sampling_interval_ms=2000)
        self.assertEqual(summary["rows"], 3)
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "sampling interval"):
            peel_codec._thermal_window_from_csv(
                jittered, 100_000_000, 3_900_000_000,
                sampling_interval_ms=2001)

        overflow = peel_codec._parse_thermal_record(
            thermal_row(1.0, edac_ce=str(1 << 64)).encode("ascii"), 1)
        self.assertTrue(overflow["hard_invalid"])
        huge = peel_codec._parse_thermal_record(
            thermal_row(1.0, edac_ce="9" * 5000).encode("ascii"), 1)
        self.assertTrue(huge["hard_invalid"])

    def test_peeltiming_request_and_csv_spellings_are_exact(self):
        with self.assertRaisesRegex(ValueError, "dimensions"):
            peel_codec.parse_peeltiming_output(
                peeltiming_stdout(warmups=1),
                **paired_parse_kwargs(warmups=1))
        stdout = peeltiming_stdout()

        def float_alias(lines):
            lines[0] = lines[0].replace(
                "loss=0.10000000000000001", "loss=1e-1")

        def decimal_seed_alias(lines):
            fields = lines[3].split(",")
            seed_index = peel_codec.PEELTIMING_COLUMNS.index(
                "construction_seed")
            fields[seed_index] = hex(int(fields[seed_index]))
            lines[3] = ",".join(fields)

        def quoted_row(lines):
            fields = lines[3].split(",")
            fields[0] = f'"{fields[0]}"'
            lines[3] = ",".join(fields)

        def quoted_header(lines):
            lines[2] = lines[2].replace("replicate", '"replicate"', 1)

        cases = (
            ("manifest", float_alias),
            ("integer", decimal_seed_alias),
            ("canonical CSV", quoted_row),
            ("header", quoted_header),
        )
        for message, mutation in cases:
            with self.subTest(alias=message):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, message):
                    peel_codec.parse_peeltiming_output(
                        authenticated_peeltiming_edit(stdout, mutation),
                        **paired_parse_kwargs())

    def test_embedded_thermal_replay_rejects_tamper_counts_and_gaps(self):
        context = paired_context()
        bad_counts = json.loads(json.dumps(context))
        bad_counts["thermal"]["valid_rows"] -= 1
        payload = {
            name: value for name, value in bad_counts["thermal"].items()
            if name != "summary_sha256"
        }
        bad_counts["thermal"]["summary_sha256"] = (
            peel_codec._canonical_json_sha256(payload))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "thermal context"):
            peel_codec._validate_paired_context(
                bad_counts, cache_state="warm", evict_bytes=4096,
                started_ns=1_000_000_000, finished_ns=2_000_000_000)

        tampered = json.loads(json.dumps(context))
        tampered["thermal"]["csv_window"] = tampered["thermal"][
            "csv_window"].replace(",60.0,", ",61.0,", 1)
        tampered["thermal"]["csv_sha256"] = hashlib.sha256(
            tampered["thermal"]["csv_window"].encode("ascii")).hexdigest()
        payload = {
            name: value for name, value in tampered["thermal"].items()
            if name != "summary_sha256"
        }
        tampered["thermal"]["summary_sha256"] = (
            peel_codec._canonical_json_sha256(payload))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "embedded CSV"):
            peel_codec._validate_paired_context(
                tampered, cache_state="warm", evict_bytes=4096,
                started_ns=1_000_000_000, finished_ns=2_000_000_000)

        gapped = thermal_csv(
            thermal_row(0.0), thermal_row(0.3), thermal_row(0.6))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "gapped"):
            peel_codec._thermal_window_from_csv(
                gapped, 100_000_000, 500_000_000,
                sampling_interval_ms=100)

    def test_thermal_missing_row_spikes_are_included_in_maxima(self):
        payload = thermal_csv(
            thermal_row(0.0),
            thermal_row(
                1.0, dimm_read_errors="1", missing_dimm=True,
                cpu="110.0", dimm="99.0"),
            thermal_row(2.0))
        summary = peel_codec._thermal_window_from_csv(
            payload, 500_000_000, 1_500_000_000,
            sampling_interval_ms=1000)
        self.assertEqual(summary["cpu_tctl_max_c"], 110.0)
        self.assertEqual(summary["dimm_max_c"], 99.0)

    @mock.patch("peel_codec._validate_runtime_identity")
    @mock.patch("peel_codec._snapshot_capture")
    def test_thermal_integrity_failure_is_not_retried(
            self, snapshot, unused_identity):
        snapshot.side_effect = peel_codec.MeasurementError(
            "path identity changed")
        frozen = paired_context()
        frozen["thermal"] = None
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "identity changed"):
            peel_codec._finalize_paired_context(
                frozen, peeltiming_stdout(), SimpleNamespace())
        self.assertEqual(snapshot.call_count, 1)

    def test_student_t_critical_matches_reference_values(self):
        references = {
            1: 12.7062047364,
            2: 4.30265272975,
            3: 3.18244630528,
            10: 2.22813885196,
            30: 2.04227245630,
        }
        for degrees_freedom, expected in references.items():
            with self.subTest(degrees_freedom=degrees_freedom):
                self.assertAlmostEqual(
                    peel_codec._student_t_critical_95(degrees_freedom),
                    expected, places=9)

    @mock.patch.dict(
        os.environ,
        {
            "WIREHAIR_V2_PEEL_DEGREES": "ambient-bad",
            "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE": "ambient-bad",
            "SAFE_PARENT_VALUE": "kept",
        },
        clear=True)
    @mock.patch("peel_codec.os.close")
    @mock.patch("peel_codec._finalize_paired_context")
    @mock.patch("peel_codec._prepare_paired_context")
    @mock.patch("peel_codec.subprocess.run")
    def test_paired_probe_runs_one_native_process_with_frozen_context(
            self, run, prepare, finalize, close):
        context = paired_context()
        capture = SimpleNamespace(fd=99)
        prepare.return_value = (context, capture)
        finalize.return_value = context
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=peeltiming_stdout(
                candidate_pmf=(0.0, 1.0), schedule="burst",
                context=context),
            stderr="")
        result = peel_codec.paired_probe(
            "bench", 64, 64, [-0.0, 1.0],
            native_profile=fake_profile(),
            target_profile=peel_codec.TARGET_PROFILE,
            seed_policy=peel_codec.TARGET_SEED_POLICY,
            loss=0.1,
            schedule="burst",
            construction_seed=9,
            loss_seed=10,
            warmup_replicates=0,
            replicates=4,
            inner_reps=1,
            max_overhead=16,
            cache_state="warm",
            evict_bytes=4096,
            context={"untrusted": "input"},
            required_margin=0.0,
        )
        self.assertTrue(result.valid_for_promotion)
        command = run.call_args.args[0]
        self.assertEqual(command[:2], ["bench", "peeltiming"])
        self.assertEqual(
            command[command.index("--candidate-pmf") + 1], "0,1")
        self.assertEqual(
            command[command.index("--schedule") + 1], "burst")
        self.assertEqual(
            command[command.index("--context-sha256") + 1],
            peel_codec._canonical_json_sha256(context["bound"]))
        environment = run.call_args.kwargs["env"]
        self.assertNotIn("WIREHAIR_V2_PEEL_DEGREES", environment)
        self.assertNotIn(
            "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE", environment)
        self.assertEqual(environment["SAFE_PARENT_VALUE"], "kept")
        finalize.assert_called_once_with(
            context, run.return_value.stdout, capture)
        close.assert_called_once_with(99)

    def test_dispatch_geometry_derivation_covers_every_rule_transition(self):
        bands = (
            (2, 5, 2),
            (6, 10, 3),
            (11, 15, 4),
            (16, 23, 5),
            (24, 31, 6),
            (32, 40, 7),
            (41, 51, 8),
            (52, 63, 9),
            (64, 77, 10),
            (78, 92, 11),
            (93, 100, 12),
        )
        covered = []
        for first, last, expected in bands:
            for block_count in range(first, last + 1):
                with self.subTest(block_count=block_count):
                    self.assertEqual(
                        peel_codec._dispatch_staircase_count(block_count),
                        expected)
                covered.append(block_count)
        self.assertEqual(covered, list(range(2, 101)))
        self.assertEqual(peel_codec._dispatch_staircase_count(101), 26)
        self.assertEqual(peel_codec._dispatch_staircase_count(1), 0)
        self.assertEqual(peel_codec._dispatch_staircase_count(64001), 0)

        # Exercise every GetDenseCount branch handoff, one descending segment
        # where Python floor division would disagree with C++, and the native
        # graph's intentional duplicate point.
        dense_seams = {
            64: 19,
            65: 26,
            500: 38,
            501: 38,
            1000: 50,
            1001: 50,
            2047: 62,
            2048: 54,
            3962: 70,
            4119: 70,
            4120: 66,
            4277: 66,
            54810: 386,
            54811: 382,
            54812: 382,
            55667: 374,
            64000: 346,
        }
        for block_count, expected in dense_seams.items():
            with self.subTest(dense_count=block_count):
                self.assertEqual(
                    peel_codec._legacy_dense_count(block_count), expected)
                if block_count > 100:
                    self.assertEqual(
                        peel_codec._dispatch_staircase_count(block_count),
                        expected)
        for block_count, raw_count in peel_codec._LEGACY_DENSE_POINTS:
            expected = raw_count + (2 - raw_count % 4) % 4
            with self.subTest(dense_point=block_count):
                self.assertEqual(
                    peel_codec._legacy_dense_count(block_count), expected)

        self.assertEqual(peel_codec._dispatch_source_hits(9999), 2)
        self.assertEqual(peel_codec._dispatch_source_hits(10000), 3)
        self.assertEqual(peel_codec._dispatch_source_hits(64000), 3)

    @mock.patch("peel_codec.subprocess.run")
    def test_invalid_schedule_and_loss_fail_before_spawn(self, run):
        for schedule in ([], {}):
            with self.subTest(schedule=schedule):
                with self.assertRaisesRegex(
                        ValueError, "invalid target schedule"):
                    peel_codec.compare_probe(
                        "unused", 64, 10, 64,
                        **probe_kwargs(schedule=schedule))
        for loss in (-0.0, 0.9900001, 1 << 4096):
            with self.subTest(loss=loss):
                with self.assertRaisesRegex(
                        ValueError, "loss must be finite"):
                    peel_codec.compare_probe(
                        "unused", 64, 10, 64,
                        **probe_kwargs(loss=loss))
        run.assert_not_called()

    def test_all_training_clis_reject_non_native_loss_values(self):
        for loss in ("-0.0", "0.9900001"):
            exact = exact_cli_args()
            exact[exact.index("--loss") + 1] = loss
            commands = (
                (peel_direct.main, [
                    "--gate-bb", "64", "--rank-bb", "64", *exact]),
                (peel_funnel.main, [
                    "--K", "64", "--gate-bb", "64", "--rank-bb", "64",
                    "--allow-unverified-cost-model", *exact]),
                (peel_sweep.main, [
                    "--allow-unverified-cost-model", *exact]),
                (peel_validate.main, [
                    "--table", "unused", "--bb", "64", *exact]),
            )
            for main, arguments in commands:
                with self.subTest(main=main.__module__, loss=loss):
                    with self.assertRaises(SystemExit):
                        main(arguments)

    def test_training_clis_reject_invalid_native_timing_dimensions_early(self):
        direct_base = [
            "--gate-bb", "64", "--rank-bb", "64",
            "--paired-context", "/tmp/context.json", *exact_cli_args(),
        ]
        funnel_base = [
            "--K", "64", "--gate-bb", "64", "--rank-bb", "64",
            "--paired-context", "/tmp/context.json",
            "--allow-unverified-cost-model", *exact_cli_args(),
        ]
        mutations = (
            ("rank payload", "--rank-bb", "1"),
            ("rank payload", "--rank-bb", "3"),
            ("inner repetitions", "--paired-inner-reps", "1025"),
            ("cold inner repetitions", "--cache-state", "cold",
             "--paired-inner-reps", "2"),
            ("replicate total", "--paired-warmups", "9998"),
            ("K-specific overhead", "--max-overhead", "577"),
            ("gate trials", "--gate-trials", str(1 << 32)),
            ("gate payload", "--gate-bb", str(1 << 31)),
            (
                "eviction cap", "--evict-bytes",
                str(peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1),
            ),
        )
        for main, base in (
                (peel_direct.main, direct_base),
                (peel_funnel.main, funnel_base)):
            for mutation in mutations:
                with self.subTest(
                        main=main.__module__, corruption=mutation[0]):
                    with self.assertRaises(SystemExit):
                        main([*base, *mutation[1:]])
        with self.assertRaises(SystemExit):
            peel_funnel.main([*funnel_base, "--threads", "1025"])

        sweep_base = [
            "--paired-context", "/tmp/context.json",
            "--allow-unverified-cost-model", *exact_cli_args(),
        ]
        for label, extra in (
                ("inner repetitions",
                 ["--paired-inner-reps", "1025"]),
                ("cold inner repetitions",
                 ["--cache-state", "cold", "--paired-inner-reps", "2"]),
                ("replicate total", ["--paired-warmups", "9998"]),
                ("K-specific overhead", ["--max-overhead", "515"]),
                (
                    "eviction cap",
                    [
                        "--evict-bytes",
                        str(peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1),
                    ],
                )):
            with self.subTest(main="peel_sweep", corruption=label):
                with self.assertRaises(SystemExit):
                    peel_sweep.main([*sweep_base, *extra])
        with self.assertRaises(SystemExit):
            peel_validate.main([
                "--table", "unused", "--bb", "64",
                "--evict-bytes",
                str(peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1),
                *exact_cli_args(),
            ])

    def test_peeltiming_memory_caps_match_native_boundaries(self):
        dimensions = {
            "block_count": 64,
            "block_bytes": 64,
            "target_profile": "dispatch-v1",
            "seed_policy": "raw",
            "construction_seed": 1,
            "loss": 0.1,
            "loss_seed": 2,
            "schedule": "iid",
            "warmup_replicates": 0,
            "replicates": 4,
            "inner_reps": 1,
            "max_overhead": 16,
            "cache_state": "warm",
            "evict_bytes": peel_codec.PEELTIMING_EVICT_BYTES_MAX,
            "required_margin": 0.0,
        }
        peel_codec.validate_peeltiming_dimensions(**dimensions)
        with self.assertRaisesRegex(ValueError, "invalid peeltiming"):
            peel_codec.validate_peeltiming_dimensions(
                **dict(
                    dimensions,
                    evict_bytes=
                        peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1))

        # Exact 4-GiB combined working-set boundary for dispatch-v1 at K=64000.
        boundary = dict(
            dimensions,
            block_count=64000,
            block_bytes=8192,
            evict_bytes=812875776,
        )
        peel_codec.validate_peeltiming_dimensions(**boundary)
        with self.assertRaisesRegex(ValueError, "invalid peeltiming"):
            peel_codec.validate_peeltiming_dimensions(
                **dict(boundary, evict_bytes=812875777))
        with self.assertRaisesRegex(ValueError, "invalid peeltiming"):
            peel_codec.validate_peeltiming_dimensions(
                **dict(
                    boundary, block_bytes=16384, evict_bytes=4096))
        peel_codec.validate_peeltiming_dimensions(
            **dict(
                boundary, block_bytes=4096,
                evict_bytes=peel_codec.PEELTIMING_EVICT_BYTES_MAX))

        # Dispatch staircase drops make the native working-set function
        # non-monotone in K.  Both endpoints fit here while the interior K
        # exceeds 4 GiB, so endpoint-only settings replay is unsound.
        nonmonotone = dict(
            dimensions,
            block_bytes=255500,
            max_overhead=0,
            evict_bytes=4096,
        )
        peel_codec.validate_peeltiming_dimensions(
            **dict(nonmonotone, block_count=2))
        peel_codec.validate_peeltiming_dimensions(
            **dict(nonmonotone, block_count=2048))
        with self.assertRaisesRegex(ValueError, "invalid peeltiming"):
            peel_codec.validate_peeltiming_dimensions(
                **dict(nonmonotone, block_count=2047))

        context = paired_context(
            evict_bytes=peel_codec.PEELTIMING_EVICT_BYTES_MAX)
        peel_codec._validate_paired_bound_context(
            context, cache_state="warm",
            evict_bytes=peel_codec.PEELTIMING_EVICT_BYTES_MAX)
        context["bound"]["evict_bytes"] += 1
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "bound context"):
            peel_codec._validate_paired_bound_context(
                context, cache_state="warm",
                evict_bytes=peel_codec.PEELTIMING_EVICT_BYTES_MAX + 1)

    def test_training_outputs_cannot_alias_executable_context_or_each_other(
            self):
        with tempfile.TemporaryDirectory() as directory:
            bench = os.path.join(directory, "bench")
            context = os.path.join(directory, "context.json")
            output = os.path.join(directory, "out.json")
            for path, payload in (
                    (bench, "executable-sentinel"),
                    (context, "context-sentinel")):
                with open(path, "w", encoding="ascii") as handle:
                    handle.write(payload)

            direct_args = [
                "--bench", bench, "--out", output,
                "--gate-bb", "64", "--rank-bb", "64",
                "--paired-context", context, *exact_cli_args(),
            ]
            self.assertEqual(
                peel_direct.main([
                    *direct_args,
                    "--out", bench,
                ]),
                2)
            hardlink = os.path.join(directory, "context-link")
            os.link(context, hardlink)
            self.assertEqual(
                peel_direct.main([
                    *direct_args,
                    "--out", hardlink,
                ]),
                2)

            sweep_args = [
                "--bench", bench, "--out", output,
                "--paired-context", context,
                "--allow-unverified-cost-model", *exact_cli_args(),
            ]
            self.assertEqual(
                peel_sweep.main([*sweep_args, "--log", bench]), 2)
            self.assertEqual(
                peel_sweep.main([*sweep_args, "--log", output]), 2)
            with open(bench, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "executable-sentinel")

            source = os.path.join(directory, "source.json")
            with open(source, "w", encoding="ascii") as handle:
                handle.write("{}")
            validate_args = [
                "--bench", bench, "--table", source, "--out", bench,
                "--bb", "64", "--paired-context", context,
                *exact_cli_args(),
            ]
            self.assertEqual(peel_validate.main(validate_args), 2)

            protected_source = os.path.join(directory, "protected.cpp")
            protected_hardlink = os.path.join(
                directory, "protected-hardlink.cpp")
            source_thermal = os.path.join(directory, "source-thermal.csv")
            source_context = os.path.join(directory, "source-context.json")
            with open(protected_source, "w", encoding="ascii") as handle:
                handle.write("source-sentinel")
            with open(source_thermal, "w", encoding="ascii") as handle:
                handle.write("thermal-sentinel")
            write_unfinalized_paired_context(source_context, source_thermal)
            os.link(protected_source, protected_hardlink)
            source_files = (("codec/protected.cpp", protected_source),)
            source_direct_args = [
                "--bench", bench, "--out", output,
                "--gate-bb", "64", "--rank-bb", "64",
                "--paired-context", source_context, *exact_cli_args(),
            ]
            source_sweep_args = [
                "--bench", bench, "--out", output,
                "--paired-context", source_context,
                "--allow-unverified-cost-model", *exact_cli_args(),
            ]
            source_validate_args = [
                "--bench", bench, "--table", source, "--out", output,
                "--bb", "64", "--paired-context", source_context,
                *exact_cli_args(),
            ]
            with mock.patch.object(
                    peel_codec, "_source_provenance_files",
                    return_value=source_files), mock.patch.object(
                        peel_direct, "capture_artifact_identity",
                        side_effect=AssertionError("source check bypassed")
                    ), mock.patch.object(
                        peel_sweep, "capture_artifact_identity",
                        side_effect=AssertionError("source check bypassed")
                    ), mock.patch.object(
                        peel_validate, "capture_artifact_identity",
                        side_effect=AssertionError("source check bypassed")):
                self.assertEqual(
                    peel_direct.main([
                        *source_direct_args,
                        "--out", protected_source,
                    ]),
                    2)
                self.assertEqual(
                    peel_sweep.main([
                        *source_sweep_args,
                        "--log", protected_hardlink,
                    ]),
                    2)
                self.assertEqual(
                    peel_validate.main([
                        *source_validate_args,
                        "--out",
                        os.path.join(directory, ".", "protected.cpp"),
                    ]),
                    2)
            with open(protected_source, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "source-sentinel")

            thermal = os.path.join(directory, "thermal.csv")
            strict_context = os.path.join(directory, "strict-context.json")
            with open(thermal, "w", encoding="ascii") as handle:
                handle.write("thermal-sentinel")
            write_unfinalized_paired_context(strict_context, thermal)
            strict_direct = [
                "--bench", bench, "--out", thermal,
                "--gate-bb", "64", "--rank-bb", "64",
                "--paired-context", strict_context, *exact_cli_args(),
            ]
            self.assertEqual(peel_direct.main(strict_direct), 2)
            strict_sweep = [
                "--bench", bench, "--out", output, "--log", thermal,
                "--paired-context", strict_context,
                "--allow-unverified-cost-model", *exact_cli_args(),
            ]
            self.assertEqual(peel_sweep.main(strict_sweep), 2)
            self.assertEqual(
                peel_validate.main([
                    "--bench", bench, "--table", source, "--out", thermal,
                    "--bb", "64", "--paired-context", strict_context,
                    *exact_cli_args(),
                ]),
                2)
            with open(thermal, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "thermal-sentinel")

    def test_result_paths_reject_source_provenance_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            generator = os.path.join(directory, "generator.py")
            source = os.path.join(directory, "source.cpp")
            generator_hardlink = os.path.join(
                directory, "generator-hardlink.py")
            source_hardlink = os.path.join(directory, "source-hardlink.cpp")
            source_symlink = os.path.join(directory, "source-symlink.cpp")
            distinct = os.path.join(directory, "result.json")
            alias_parent = os.path.join(directory, "alias-parent")
            git_directory = os.path.join(directory, "git-worktree")
            common_directory = os.path.join(directory, "git-common")
            os.mkdir(alias_parent)
            os.mkdir(git_directory)
            os.mkdir(common_directory)
            for path, payload in (
                    (generator, "generator"),
                    (source, "source"),
                    (distinct, "result")):
                with open(path, "w", encoding="ascii") as handle:
                    handle.write(payload)
            os.link(generator, generator_hardlink)
            os.link(source, source_hardlink)
            os.symlink(source, source_symlink)

            source_files = (("codec/source.cpp", source),)
            head = os.path.join(git_directory, "HEAD")
            index = os.path.join(git_directory, "index")
            packed_refs = os.path.join(common_directory, "packed-refs")
            git_locator = os.path.join(directory, "worktree-dot-git")
            symbolic_ref = os.path.join(
                common_directory, "refs", "heads", "topic")
            os.makedirs(os.path.dirname(symbolic_ref))
            for path in (
                    head, index, packed_refs, symbolic_ref, git_locator):
                with open(path, "w", encoding="ascii") as handle:
                    handle.write(path)
            head_hardlink = os.path.join(directory, "head-hardlink")
            index_hardlink = os.path.join(directory, "index-hardlink")
            locator_hardlink = os.path.join(directory, "locator-hardlink")
            ref_symlink = os.path.join(directory, "ref-symlink")
            os.link(head, head_hardlink)
            os.link(index, index_hardlink)
            os.link(git_locator, locator_hardlink)
            os.symlink(symbolic_ref, ref_symlink)
            git_paths = (
                (
                    ("the worktree Git directory", git_directory),
                    ("the common Git directory", common_directory),
                ),
                (
                    ("the worktree Git locator", git_locator),
                    ("Git HEAD", head),
                    ("the Git index", index),
                    ("Git packed refs", packed_refs),
                    ("the Git symbolic HEAD ref", symbolic_ref),
                ),
            )
            aliases = (
                generator,
                generator_hardlink,
                source,
                os.path.join(alias_parent, "..", "source.cpp"),
                source_hardlink,
                source_symlink,
                git_directory,
                os.path.join(git_directory, "future", "metadata"),
                os.path.join(common_directory, "refs", "future"),
                head_hardlink,
                index_hardlink,
                locator_hardlink,
                ref_symlink,
            )
            with mock.patch.object(
                    peel_codec, "_git_provenance_paths",
                    return_value=git_paths), mock.patch.object(
                    peel_codec, "_source_provenance_files",
                    return_value=source_files):
                for alias in aliases:
                    with self.subTest(alias=alias):
                        with self.assertRaisesRegex(
                                peel_codec.MeasurementError,
                                "must not replace"):
                            peel_codec.require_distinct_from_source_provenance(
                                alias, "--out", generator)
                peel_codec.require_distinct_from_source_provenance(
                    distinct, "--out", generator)

    def test_source_provenance_includes_direct_build_and_ignore_inputs(self):
        source_files = dict(peel_codec._source_provenance_files())
        required_inputs = (
            "abi/wirehair.map",
            "cmake/wirehair.pc.in",
            "cmake/wirehairConfig.cmake.in",
            ".gitignore",
            ".beads/.gitignore",
            "@git/info-exclude",
        )
        for relative in required_inputs:
            with self.subTest(relative=relative):
                self.assertIn(relative, source_files)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "must not (?:create or )?replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        source_files[relative], "--out",
                        "tools/peel_direct.py")
        with tempfile.TemporaryDirectory() as directory:
            worktree_alias = os.path.join(directory, "worktree-alias")
            os.symlink(peel_codec.REPO_ROOT, worktree_alias)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "source-provenance input"):
                peel_codec.require_distinct_from_source_provenance(
                    os.path.join(
                        worktree_alias, "future-source", ".gitignore"),
                    "--out", "tools/peel_direct.py")

    def test_source_provenance_includes_self_ignored_gitignore(self):
        with tempfile.TemporaryDirectory() as directory:
            repository = os.path.join(directory, "repository")
            subprocess.run(
                ["git", "init", "-q", repository],
                check=True, capture_output=True)
            repository_ignore = os.path.join(repository, ".gitignore")
            with open(repository_ignore, "w", encoding="ascii") as handle:
                handle.write("*\n")
            ignore_hardlink = os.path.join(directory, "ignore-hardlink")
            os.link(repository_ignore, ignore_hardlink)

            with mock.patch.object(peel_codec, "REPO_ROOT", repository):
                source_files = dict(peel_codec._source_provenance_files())
                self.assertIn(".gitignore", source_files)
                self.assertEqual(
                    source_files[".gitignore"], repository_ignore)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        ignore_hardlink, "--log", "tools/peel_sweep.py")

    def test_git_provenance_protects_relative_xdg_and_missing_configs(self):
        with tempfile.TemporaryDirectory() as directory:
            missing_global = os.path.join(directory, "missing-global-config")
            with mock.patch.dict(
                    os.environ,
                    {
                        "GIT_CONFIG_GLOBAL": missing_global,
                        "GIT_CONFIG_NOSYSTEM": "1",
                        "XDG_CONFIG_HOME": "relative-xdg",
                    },
                    clear=False):
                expected_excludes = os.path.join(
                    peel_codec.REPO_ROOT,
                    "relative-xdg", "git", "ignore")
                self.assertEqual(
                    peel_codec._effective_global_git_exclude_path(),
                    expected_excludes)
                self.assertIn(
                    missing_global,
                    peel_codec._git_configuration_provenance_paths())
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        missing_global, "--out", "tools/peel_direct.py")

            previous_directory = os.getcwd()
            try:
                os.chdir(directory)
                with mock.patch.dict(
                        os.environ,
                        {
                            "HOME": "relative-home",
                            "GIT_CONFIG_GLOBAL": "",
                            "GIT_CONFIG_NOSYSTEM": "1",
                        },
                        clear=False):
                    os.environ.pop("XDG_CONFIG_HOME", None)
                    self.assertEqual(
                        peel_codec._effective_global_git_exclude_path(),
                        os.path.join(
                            peel_codec.REPO_ROOT, "relative-home", ".config",
                            "git", "ignore"))
            finally:
                os.chdir(previous_directory)

            newline_config = os.path.join(directory, "global\nconfig")
            newline_alias = os.path.join(directory, "newline-config-alias")
            with open(newline_config, "w", encoding="ascii") as handle:
                handle.write("")
            os.link(newline_config, newline_alias)
            with mock.patch.dict(
                    os.environ,
                    {
                        "GIT_CONFIG_GLOBAL": newline_config,
                        "GIT_CONFIG_NOSYSTEM": "1",
                    },
                    clear=False):
                configuration_paths = (
                    peel_codec._git_configuration_provenance_paths())
                self.assertIn(newline_config, configuration_paths)
                self.assertNotIn(
                    newline_config.partition("\n")[0], configuration_paths)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        newline_alias, "--log", "tools/peel_sweep.py")

            literal_tilde = os.path.join(
                peel_codec.REPO_ROOT, "~", "literal-global-config")
            with mock.patch.dict(
                    os.environ,
                    {
                        "GIT_CONFIG_GLOBAL": "~/literal-global-config",
                        "GIT_CONFIG_NOSYSTEM": "1",
                    },
                    clear=False):
                configuration_paths = (
                    peel_codec._git_configuration_provenance_paths())
                self.assertIn(literal_tilde, configuration_paths)
                self.assertNotIn(
                    os.path.expanduser("~/literal-global-config"),
                    configuration_paths)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        literal_tilde, "--out", "tools/peel_direct.py")

            for label, home, xdg, expected_relatives in (
                    (
                        "newline-home",
                        "relative\nhome",
                        None,
                        (
                            os.path.join(
                                "relative\nhome", ".config", "git", "config"),
                            os.path.join("relative\nhome", ".gitconfig"),
                        ),
                    ),
                    (
                        "newline-xdg",
                        "relative-home",
                        "relative\nxdg",
                        (
                            os.path.join("relative\nxdg", "git", "config"),
                            os.path.join("relative-home", ".gitconfig"),
                        ),
                    )):
                with self.subTest(global_default=label), mock.patch.dict(
                        os.environ,
                        {
                            "HOME": home,
                            "GIT_CONFIG_NOSYSTEM": "1",
                        },
                        clear=False):
                    os.environ.pop("GIT_CONFIG_GLOBAL", None)
                    if xdg is None:
                        os.environ.pop("XDG_CONFIG_HOME", None)
                    else:
                        os.environ["XDG_CONFIG_HOME"] = xdg
                    configuration_paths = (
                        peel_codec._git_configuration_provenance_paths())
                    for relative in expected_relatives:
                        self.assertIn(
                            os.path.join(peel_codec.REPO_ROOT, relative),
                            configuration_paths)

            including_config = os.path.join(directory, "including-config")
            missing_include = os.path.join(directory, "missing-include")
            with open(including_config, "w", encoding="ascii") as handle:
                handle.write("[include]\n\tpath = missing-include\n")
            with mock.patch.dict(
                    os.environ,
                    {
                        "GIT_CONFIG_GLOBAL": including_config,
                        "GIT_CONFIG_NOSYSTEM": "1",
                    },
                    clear=False):
                self.assertIn(
                    missing_include,
                    peel_codec._git_configuration_provenance_paths())
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        missing_include, "--log", "tools/peel_sweep.py")

    def test_git_provenance_protects_linked_worktree_metadata_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            repository = os.path.join(directory, "repository")
            worktree = os.path.join(directory, "worktree")
            subprocess.run(
                ["git", "init", "-q", repository],
                check=True, capture_output=True)
            tracked_source = os.path.join(repository, "tracked.cpp")
            with open(tracked_source, "w", encoding="ascii") as handle:
                handle.write("int tracked;\n")
            subprocess.run(
                ["git", "-C", repository, "add", "tracked.cpp"],
                check=True, capture_output=True)
            subprocess.run(
                [
                    "git", "-C", repository,
                    "-c", "user.name=Peel Test",
                    "-c", "user.email=peel-test@example.invalid",
                    "commit", "-q", "-m", "initial",
                ],
                check=True, capture_output=True)
            subprocess.run(
                [
                    "git", "-C", repository, "worktree", "add", "-q",
                    "-b", "linked-test", worktree,
                ],
                check=True, capture_output=True)
            subprocess.run(
                ["git", "-C", worktree, "update-index", "--split-index"],
                check=True, capture_output=True)
            subprocess.run(
                [
                    "git", "-C", worktree, "config",
                    "core.excludesFile", "local-exclude",
                ],
                check=True, capture_output=True)
            local_config_result = subprocess.run(
                [
                    "git", "-C", worktree, "rev-parse",
                    "--path-format=absolute", "--git-path", "config",
                ],
                check=True, capture_output=True, text=True)
            local_config = local_config_result.stdout.rstrip("\n")
            symbolic_paths = {}
            for name in ("HEAD", "refs/heads/alias", "refs/heads/linked-test"):
                result = subprocess.run(
                    [
                        "git", "-C", worktree, "rev-parse",
                        "--path-format=absolute", "--git-path", name,
                    ],
                    check=True, capture_output=True, text=True)
                symbolic_paths[name] = result.stdout.rstrip("\n")
            os.makedirs(
                os.path.dirname(symbolic_paths["refs/heads/alias"]),
                exist_ok=True)
            with open(
                    symbolic_paths["refs/heads/alias"], "w",
                    encoding="ascii") as handle:
                handle.write("ref: refs/heads/linked-test\n")
            with open(
                    symbolic_paths["HEAD"], "w",
                    encoding="ascii") as handle:
                handle.write("ref: refs/heads/alias\n")
            override_config = os.path.join(directory, "override-config")
            with open(override_config, "w", encoding="ascii") as handle:
                handle.write("[core]\n\texcludesFile = override-exclude\n")
            commondir_alias = os.path.join(directory, "commondir-alias")
            shared_index_alias = os.path.join(
                directory, "shared-index-alias")
            local_config_alias = os.path.join(
                directory, "local-config-alias")
            symbolic_ref_alias = os.path.join(
                directory, "symbolic-ref-alias")
            os.link(local_config, local_config_alias)
            os.link(
                symbolic_paths["refs/heads/alias"], symbolic_ref_alias)

            with mock.patch.object(
                    peel_codec, "REPO_ROOT", worktree), mock.patch.dict(
                        os.environ, {"GIT_CONFIG": override_config},
                        clear=False):
                unused_directories, git_files = (
                    peel_codec._git_provenance_paths())
                metadata = dict(git_files)
                commondir = metadata["the Git common-directory locator"]
                shared_index = metadata["the Git shared index"]
                self.assertEqual(
                    metadata["the Git symbolic HEAD ref"],
                    symbolic_paths["refs/heads/alias"])
                self.assertIn(
                    symbolic_paths["refs/heads/linked-test"],
                    metadata.values())
                self.assertTrue(os.path.isfile(commondir))
                self.assertTrue(os.path.isfile(shared_index))
                self.assertIn(
                    local_config,
                    peel_codec._git_configuration_provenance_paths())
                self.assertEqual(
                    peel_codec._effective_global_git_exclude_path(),
                    os.path.join(worktree, "local-exclude"))
                os.link(commondir, commondir_alias)
                os.link(shared_index, shared_index_alias)
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        commondir_alias, "--log", "tools/peel_sweep.py")
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        shared_index_alias, "--log", "tools/peel_sweep.py")
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        local_config_alias, "--log", "tools/peel_sweep.py")
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "must not replace"):
                    peel_codec.require_distinct_from_source_provenance(
                        symbolic_ref_alias, "--log", "tools/peel_sweep.py")

    def test_training_outputs_reject_actual_git_provenance_paths(self):
        git_directories, git_files = peel_codec._git_provenance_paths()
        directories = dict(git_directories)
        metadata = dict(git_files)
        git_directory = directories["the worktree Git directory"]
        common_directory = directories.get(
            "the common Git directory", git_directory)
        git_locator = metadata["the worktree Git locator"]
        head = metadata["Git HEAD"]
        index = metadata["the Git index"]
        info_exclude = metadata["Git exclude input @git/info-exclude"]
        global_exclude = metadata["Git exclude input @git/global-excludes"]
        configuration = next(
            path for label, path in git_files
            if label.startswith("Git configuration input "))
        reference = metadata.get(
            "the Git symbolic HEAD ref",
            os.path.join(common_directory, "refs", "heads", "future"))
        repository_ignore = os.path.join(peel_codec.REPO_ROOT, ".gitignore")
        future_ignore = os.path.join(
            peel_codec.REPO_ROOT, "future-source", ".gitignore")
        future_source = os.path.join(
            peel_codec.REPO_ROOT, "future-source", "receipt.cpp")
        future_cmake = os.path.join(
            peel_codec.REPO_ROOT, "future-source", "CMakeLists.txt")

        # Use the worktree filesystem so the hard-link alias of Git's
        # per-repository excludes file is representable.
        with tempfile.TemporaryDirectory(
                dir=peel_codec.REPO_ROOT) as directory:
            bench = os.path.join(directory, "bench")
            table = os.path.join(directory, "table.json")
            context = os.path.join(directory, "context.json")
            thermal = os.path.join(directory, "thermal.csv")
            output = os.path.join(directory, "result.json")
            log = os.path.join(directory, "result.log")
            for path, payload in (
                    (bench, "benchmark"),
                    (table, "{}"),
                    (thermal, "thermal")):
                with open(path, "w", encoding="ascii") as handle:
                    handle.write(payload)
            write_unfinalized_paired_context(context, thermal)
            info_exclude_hardlink = os.path.join(
                directory, "info-exclude-hardlink")
            os.link(info_exclude, info_exclude_hardlink)
            dangling_source_log = os.path.join(
                directory, "dangling-source-log")
            os.symlink(future_source, dangling_source_log)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "source-provenance input"):
                peel_codec.require_distinct_from_source_provenance(
                    dangling_source_log, "--log", "tools/peel_sweep.py")
            peel_codec.require_distinct_from_source_provenance(
                output, "--out", "tools/peel_direct.py")
            peel_codec.require_distinct_from_source_provenance(
                log, "--log", "tools/peel_sweep.py")

            with mock.patch.object(
                    peel_direct, "capture_artifact_identity",
                    side_effect=AssertionError("Git check bypassed")
            ), mock.patch.object(
                    peel_sweep, "capture_artifact_identity",
                    side_effect=AssertionError("Git check bypassed")
            ), mock.patch.object(
                    peel_validate, "capture_artifact_identity",
                    side_effect=AssertionError("Git check bypassed")):
                direct_args = [
                    "--bench", bench, "--out", head,
                    "--gate-bb", "64", "--rank-bb", "64",
                    "--paired-context", context, *exact_cli_args(),
                ]
                self.assertEqual(peel_direct.main(direct_args), 2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", git_locator,
                    ]),
                    2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", future_source,
                    ]),
                    2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", os.path.join(
                            git_directory, "future", "result.json"),
                    ]),
                    2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", repository_ignore,
                    ]),
                    2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", future_ignore,
                    ]),
                    2)
                self.assertEqual(
                    peel_direct.main([
                        *direct_args,
                        "--out", configuration,
                    ]),
                    2)

                sweep_args = [
                    "--bench", bench, "--out", output, "--log", index,
                    "--paired-context", context,
                    "--allow-unverified-cost-model", *exact_cli_args(),
                ]
                self.assertEqual(peel_sweep.main(sweep_args), 2)
                self.assertEqual(
                    peel_sweep.main([
                        *sweep_args,
                        "--log", info_exclude_hardlink,
                    ]),
                    2)
                self.assertEqual(
                    peel_sweep.main([
                        *sweep_args,
                        "--log", dangling_source_log,
                    ]),
                    2)

                validate_args = [
                    "--bench", bench, "--table", table, "--out", reference,
                    "--bb", "64", "--paired-context", context,
                    *exact_cli_args(),
                ]
                self.assertEqual(peel_validate.main(validate_args), 2)
                self.assertEqual(
                    peel_validate.main([
                        *validate_args,
                        "--out", global_exclude,
                    ]),
                    2)
                self.assertEqual(
                    peel_validate.main([
                        *validate_args,
                        "--out", future_cmake,
                    ]),
                    2)

    @mock.patch("peel_codec.subprocess.run")
    def test_nonzero_compare_exit_fails_closed(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 7, stdout=compare_stdout(), stderr="fatal")
        with self.assertRaisesRegex(
            peel_codec.MeasurementError, "exited 7.*fatal"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, **probe_kwargs())

    @mock.patch.dict(
        os.environ,
        {
            "WIREHAIR_V2_PEEL_DEGREES": "inherited-bad-value",
            "WIREHAIR_V2_OTHER_TEST_HOOK": "also-bad",
            "WIREHAIR_V2_STAIRCASE_DEGREES": "forbidden",
            "WIREHAIR_V2_STAIRCASE_ROW_DEGREES": "forbidden",
            "WIREHAIR_V2_BAND_TRACKING_X": "forbidden",
            "SAFE_PARENT_VALUE": "kept",
        },
        clear=True)
    @mock.patch("peel_codec.subprocess.run")
    def test_compare_isolates_environment_and_preserves_metrics(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=compare_stdout(
                pmf_digest=peel_codec.pmf_sha256([0.25, 0.75]),
                staircase_scale="2.5"),
            stderr="")
        result = peel_codec.compare_probe(
            "bench", 64, 10, 64,
            peel_weights=[0.25, 0.75],
            degree_scale=2.5,
            **probe_kwargs())

        command = run.call_args.args[0]
        environment = run.call_args.kwargs["env"]
        self.assertEqual(
            command[command.index("--construction-seed") + 1],
            str(0x1122334455667788))
        self.assertEqual(
            command[command.index("--loss-seed") + 1],
            str(0x8877665544332211))
        self.assertNotIn("--mixed-gf16-rows", command)
        self.assertEqual(
            command[command.index("--max-message-mib") + 1], "0")
        self.assertNotIn("WIREHAIR_V2_OTHER_TEST_HOOK", environment)
        self.assertNotIn("WIREHAIR_V2_STAIRCASE_DEGREES", environment)
        self.assertNotIn("WIREHAIR_V2_STAIRCASE_ROW_DEGREES", environment)
        self.assertNotIn("WIREHAIR_V2_BAND_TRACKING_X", environment)
        self.assertEqual(environment["SAFE_PARENT_VALUE"], "kept")
        self.assertEqual(
            environment["WIREHAIR_V2_PEEL_DEGREES"], "0.25,0.75")
        self.assertEqual(
            environment["WIREHAIR_V2_STAIRCASE_DEGREE_SCALE"], "2.5")
        self.assertEqual(result.oh_sd, 0.5)
        self.assertEqual(result.oh50, 0.0)
        self.assertEqual(result.oh95, 1.0)
        self.assertEqual(result.oh99, 2.0)
        self.assertEqual(result.oh_max, 3.0)
        self.assertEqual(result.decode_mbps, 1234.5)
        self.assertEqual(result.construction_seed, 0x1122334455667788)
        self.assertEqual(result.loss_seed, 0x8877665544332211)

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_canonicalizes_negative_zero_scale(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=compare_stdout(staircase_scale="0"),
            stderr="")
        result = peel_codec.compare_probe(
            "bench", 64, 10, 64, degree_scale=-0.0, **probe_kwargs())

        environment = run.call_args.kwargs["env"]
        self.assertEqual(
            environment["WIREHAIR_V2_STAIRCASE_DEGREE_SCALE"], "0")
        self.assertEqual(result.target_receipt["staircase_scale"], "0")
        self.assertEqual(
            peel_codec._canonical_staircase_scale_spec(-0.0), "0")

    @mock.patch("peel_codec.subprocess.run")
    def test_exact_native_pmf_and_metadata_are_parsed(self, run):
        probabilities = [0.125] + [0.875 / 63.0] * 63
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=peelpmf_stdout(probabilities), stderr="")

        profile = peel_codec.stock_profile(
            sys.executable, 64, target_profile="dispatch-v1")
        self.assertEqual(profile.block_count, 64)
        self.assertEqual(profile.staircase, 10)
        self.assertEqual(profile.source_hits, 2)
        self.assertEqual(profile.target_mean, 12.8)
        self.assertAlmostEqual(profile.pmf[0], 0.125)
        self.assertNotAlmostEqual(profile.pmf[0], 1.0 / 64.0)
        self.assertEqual(
            run.call_args.args[0],
            [
                sys.executable, "peelpmf", "--N", "64",
                "--target-profile", "dispatch-v1",
            ])
        self.assertFalse(any(
            key.startswith("WIREHAIR_V2_")
            for key in run.call_args.kwargs["env"]))

        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=peelpmf_stdout(
                probabilities,
                target_mean=math.nextafter(12.8, math.inf)),
            stderr="")
        peel_codec._STOCK_CACHE.clear()
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid peelpmf metadata"):
            peel_codec.stock_profile(
                sys.executable, 64, target_profile="dispatch-v1")

    @mock.patch("peel_codec.subprocess.run")
    def test_native_pmf_rejects_self_consistent_wrong_dispatch_geometry(
            self, run):
        probabilities = [1.0 / 64.0] * 64
        for label, staircase, source_hits in (
                ("staircase", 11, 2),
                ("source hits", 10, 3)):
            with self.subTest(label=label):
                target_mean = (
                    64.0 * min(source_hits, staircase) / staircase
                )
                run.return_value = subprocess.CompletedProcess(
                    ["bench"], 0,
                    stdout=peelpmf_stdout(
                        probabilities,
                        staircase=staircase,
                        source_hits=source_hits,
                        target_mean=target_mean),
                    stderr="")
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "invalid peelpmf metadata"):
                    peel_codec.stock_profile(
                        sys.executable, 64, target_profile="dispatch-v1")

    def test_pair_seed_is_shared_but_domains_are_separate(self):
        direct_rank = peel_codec.derive_seed(
            7, "direct-search", 64, "rank", 20, 4096, "construction")
        loss_rank = peel_codec.derive_seed(
            8, "direct-search", 64, "rank", 20, 4096, "loss")
        self.assertEqual(direct_rank, 0x9e8154799396c85f)
        self.assertEqual(loss_rank, 0x67db35c40e3e3531)
        self.assertEqual(
            direct_rank,
            peel_codec.derive_seed(
                7, "direct-search", 64, "rank", 20, 4096,
                "construction"))
        self.assertNotEqual(
            direct_rank,
            peel_codec.derive_seed(
                7, "validation", 64, 20, 4096, "construction"))

    def test_family_uses_validated_numeric_conversion(self):
        stock = [0.5, 0.5]
        self.assertEqual(
            peel_codec.family(stock, "200", "0", "2", "100"),
            peel_codec.family(stock, 200, 0, 2, 100))
        huge = 1 << 4096
        self.assertIsNone(
            peel_codec.family(stock, huge, 0, 2, 100))
        with self.assertRaisesRegex(ValueError, "invalid PMF"):
            peel_codec.pmf_sha256([huge, 1.0])
        with self.assertRaisesRegex(ValueError, "invalid peel weight vector"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64,
                peel_weights=[huge, 1.0], **probe_kwargs())
        with self.assertRaisesRegex(
                ValueError, "invalid staircase degree scale"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64, degree_scale=huge,
                **probe_kwargs())
        for trials, block_bytes in (
                (peel_codec._COMPARE_TRIALS_MAX + 1, 64),
                (10, peel_codec._COMPARE_BLOCK_BYTES_MAX + 1),
                (10, 63)):
            with self.subTest(trials=trials, block_bytes=block_bytes):
                with self.assertRaisesRegex(
                        ValueError, "invalid compare K, trial count, or block"):
                    peel_codec.compare_probe(
                        "unused", 64, trials, block_bytes, **probe_kwargs())

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_boolean_numeric_aliases(self, run):
        with self.assertRaisesRegex(ValueError, "invalid peel weight vector"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64,
                peel_weights=(value for value in (True, 1.0)),
                **probe_kwargs())
        with self.assertRaisesRegex(
                ValueError, "invalid staircase degree scale"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64, degree_scale=True, **probe_kwargs())
        run.assert_not_called()

    @mock.patch("peel_codec.subprocess.run")
    def test_native_pmf_rejects_staircase_outside_production_span(self, run):
        staircase = peel_codec._production_staircase_max(64) + 1
        probabilities = [1.0 / 64.0] * 64
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=peelpmf_stdout(
                probabilities, staircase=staircase,
                target_mean=64.0 * 2.0 / staircase),
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid peelpmf metadata"):
            peel_codec.stock_profile(
                sys.executable, 64, target_profile="dispatch-v1")

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_semantically_wrong_n_mean(self, run):
        stdout = compare_stdout().replace(
            "v2_target 64 10 0 64 ", "v2_target 64 10 0 65 ")
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=stdout, stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "semantically wrong compare row"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, **probe_kwargs())

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_unstructured_trailing_output(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=compare_stdout() + "unexpected junk\n",
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "unexpected compare output"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, **probe_kwargs())

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_mutated_exact_target_receipt(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=compare_stdout().replace(
                ",dense_rows=4,", ",dense_rows=12,", 1),
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "wh2_target receipt does not match"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, **probe_kwargs())

        scaled = compare_stdout(staircase_scale="12.345")
        self.assertIn(",target_mean=12.300000000000001,", scaled)
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=scaled.replace(
                ",target_mean=12.300000000000001,",
                ",target_mean=12.345,", 1),
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "wh2_target receipt does not match"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, degree_scale=12.345,
                **probe_kwargs())

    @mock.patch("peel_codec.subprocess.run")
    def test_native_pmf_rejects_partial_or_trailing_output(self, run):
        probabilities = [1.0 / 64.0] * 64
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=peelpmf_stdout(probabilities[:-1]), stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "63 of 64 degrees"):
            peel_codec.stock_profile(
                sys.executable, 64, target_profile="dispatch-v1")

        peel_codec._STOCK_CACHE.clear()
        output = peelpmf_stdout(probabilities) + "unexpected junk\n"
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=output, stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "unexpected peelpmf output"):
            peel_codec.stock_profile(
                sys.executable, 64, target_profile="dispatch-v1")

    def test_source_identity_includes_compiled_root_sources(self):
        original = peel_codec._sha256_file
        with mock.patch(
                "peel_codec._sha256_file", wraps=original) as digest:
            peel_codec.capture_artifact_identity(
                sys.executable, "tools/peel_direct.py")
        paths = [call.args[0] for call in digest.call_args_list]
        self.assertTrue(any(
            path.endswith("/WirehairTools.cpp") for path in paths))

    @mock.patch("peel_codec._run_checked")
    @mock.patch("peel_codec._benchmark_identity")
    def test_stock_cache_is_keyed_by_executable_digest(self, identity, run):
        identity.side_effect = [
            {"path": "/bench", "sha256": "a" * 64, "size": 1},
            {"path": "/bench", "sha256": "a" * 64, "size": 1},
            {"path": "/bench", "sha256": "b" * 64, "size": 1},
            {"path": "/bench", "sha256": "b" * 64, "size": 1},
        ]
        probabilities = [1.0 / 64.0] * 64
        run.return_value = peelpmf_stdout(probabilities)
        peel_codec.stock_profile(
            "/bench", 64, target_profile="dispatch-v1")
        peel_codec.stock_profile(
            "/bench", 64, target_profile="dispatch-v1")
        self.assertEqual(run.call_count, 2)
        self.assertTrue(all(
            key[2:5] == (
                "dispatch-v1",
                peel_codec.TARGET_CONTRACT["contract_id"],
                peel_codec.TARGET_CONTRACT["contract_sha256"])
            for key in peel_codec._STOCK_CACHE))

    @mock.patch("peel_codec._artifact_identity")
    def test_publication_rejects_artifact_drift(self, identity):
        before = {
            "benchmark": {
                "path": "/bench", "sha256": "a" * 64, "size": 1,
            },
            "source": {
                "git_commit": "a" * 40,
                "state_sha256": "b" * 64,
                "file_count": 1,
                "generator_sha256": "c" * 64,
            },
        }
        after = json.loads(json.dumps(before))
        after["benchmark"]["sha256"] = "d" * 64
        identity.return_value = after
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "changed during measurement"):
            peel_codec.make_table_document(
                {64: complete_entry(64)},
                generator="tools/peel_direct.py",
                bench="/bench",
                settings=direct_settings(),
                artifact_identity=before,
                **table_kwargs(),
            )

    def test_benchmark_identity_mismatch_is_refused(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            settings=direct_settings(),
            **table_kwargs(),
        )
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "benchmark identity mismatch"):
            peel_codec.verify_benchmark_identity(document, "/bin/true")


class LegacyPeelTableV3Examples:
    def test_table_constructor_rejects_unhashable_schedule_and_bad_loss(self):
        for field, value in (
                ("schedule", []),
                ("schedule", {}),
                ("loss", -0.0),
                ("loss", 0.9900001),
                ("loss", 1 << 4096)):
            with self.subTest(field=field, value=value):
                kwargs = table_kwargs()
                kwargs[field] = value
                with self.assertRaises(ValueError):
                    peel_codec.make_table_document(
                        {64: complete_entry(64)},
                        generator="tools/peel_direct.py",
                        bench=sys.executable,
                        settings=direct_settings(),
                        **kwargs,
                    )

    def test_schema_requires_artifact_native_and_complete_search_receipts(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            settings=direct_settings(),
            **table_kwargs(),
        )
        provenance = document["provenance"]
        self.assertEqual(len(provenance["benchmark"]["sha256"]), 64)
        self.assertEqual(len(provenance["source"]["state_sha256"]), 64)
        self.assertEqual(
            peel_codec.validate_table_document(document)[64]["K"], 64)
        removals = (
            ("benchmark digest", ("provenance", "benchmark", "sha256")),
            ("metric scope", ("provenance", "recovery_metric_scope")),
            ("native metadata", ("entries", "64", "native")),
            ("search tail metric", (
                "entries", "64", "search_receipt", "OH99")),
        )
        for label, path in removals:
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                target = damaged
                for component in path[:-1]:
                    target = target[component]
                del target[path[-1]]
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        for mutation in ("wrong", "missing"):
            with self.subTest(target_mean=mutation):
                damaged = json.loads(json.dumps(document))
                target = damaged["entries"]["64"]["search_receipt"][
                    "target_receipt"]
                if mutation == "wrong":
                    target["target_mean"] = "12.7"
                else:
                    del target["target_mean"]
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        contradictions = (
            ("forged goodput", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__(
                 "goodput", 123456.0)),
            ("forged construction seed", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__(
                 "construction_seed", 5)),
            ("wrong target contract", lambda d:
             d["entries"]["64"]["search_receipt"]["target_receipt"].__setitem__(
                 "contract_id", "0000000000000000")),
            ("wrong target D2", lambda d:
             d["entries"]["64"]["search_receipt"]["target_receipt"].__setitem__(
                 "dense_rows", "12")),
            ("wrong target seed policy", lambda d:
             d["entries"]["64"]["search_receipt"]["target_receipt"].__setitem__(
                 "seed_policy", "table")),
            ("candidate/control scale mismatch", lambda d:
             d["entries"]["64"]["search_receipt"]["shipped_control"][
                 "target_receipt"].__setitem__("staircase_scale", "12")),
            ("oversized scale", lambda d:
             d["entries"]["64"].__setitem__("scale", 64000.01)),
            ("wrong seed domain", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__(
                 "seed_domain", "funnel-search")),
            ("wrong sampling seed", lambda d:
             d["entries"]["64"]["search_receipt"]["context"].__setitem__(
                 "sampling_seed", 0)),
            ("unexpected direct context", lambda d:
             d["entries"]["64"]["search_receipt"]["context"].__setitem__(
                 "extra", 1)),
            ("PMF-coordinate contradiction", lambda d:
             d["entries"]["64"]["search_receipt"].update({
                 "peel_pmf": [0.5, 0.5],
                 "peel_pmf_sha256": peel_codec.pmf_sha256([0.5, 0.5]),
             })),
            ("integer coordinate encoded as float", lambda d: (
             d["entries"]["64"].__setitem__("p1", 100.0),
             d["entries"]["64"]["search_receipt"]["coordinates"].__setitem__(
                 "p1", 100.0))),
            ("wrong native shipped mean", lambda d:
             d["entries"]["64"]["native"].__setitem__(
                 "target_mean", 99.0)),
            ("boolean top-level fail alias", lambda d:
             d["entries"]["64"].__setitem__("fail", False)),
            ("boolean top-level quantile alias", lambda d:
             d["entries"]["64"].__setitem__("OH95", True)),
            ("obsolete scaled-control marker", lambda d:
             d["entries"]["64"].__setitem__("reverted_to_control", True)),
        )
        for label, damage in contradictions:
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                damage(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        for field, value in (
                ("schedule", []),
                ("schedule", {}),
                ("loss", -0.0),
                ("loss", 0.9900001)):
            with self.subTest(measurement_policy=field, value=value):
                damaged = json.loads(json.dumps(document))
                damaged["provenance"]["measurement_policy"][field] = value
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        damaged = json.loads(json.dumps(document))
        native = damaged["entries"]["64"]["native"]
        native["pmf"] = [1.0 / 63.0] * 63
        native["pmf_sha256"] = peel_codec.pmf_sha256(native["pmf"])
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid native metadata"):
            peel_codec.validate_table_document(damaged)

        trained_document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        forged_tie = json.loads(json.dumps(trained_document))
        forged_entry = forged_tie["entries"]["64"]
        forged_receipt = forged_entry["search_receipt"]
        forged_control = forged_receipt["shipped_control"]
        for name in (
                "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
                "OH_max", "decode_mbps"):
            forged_control[name] = forged_receipt[name]
        canonical_tie = peel_codec._canonical_goodput(forged_receipt, 64)
        forged_receipt["shipped_goodput"] = canonical_tie
        forged_receipt["goodput"] = math.nextafter(
            canonical_tie, math.inf)
        forged_entry["goodput"] = forged_receipt["goodput"]
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "inconsistent goodput"):
            peel_codec.validate_table_document(forged_tie)

        damaged = json.loads(json.dumps(trained_document))
        receipt = damaged["entries"]["64"]["search_receipt"]
        receipt["shipped_control"]["decode_mbps"] = 1100.0
        receipt["shipped_goodput"] = (
            1100.0 * 64.0 /
            (64.0 + receipt["shipped_control"]["oh_mean"]))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "did not beat shipped"):
            peel_codec.validate_table_document(damaged)

        damaged = json.loads(json.dumps(trained_document))
        damaged["entries"]["64"]["scale"] = 12.0
        damaged["entries"]["64"]["search_receipt"]["coordinates"][
            "scale"] = 12.0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid search coordinates"):
            peel_codec.validate_table_document(damaged)

    def test_entry_top_level_schema_is_generator_and_state_exact(self):
        trained = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        for field, value in (
                ("unknown_contract_field", "forged"),
                ("target_profile", "not-the-target"),
                ("construction_seed_base", 0),
                ("verified_mbps", 1e300),
                ("reverted_to_shipped", False)):
            with self.subTest(unexpected_field=field):
                damaged = json.loads(json.dumps(trained))
                damaged["entries"]["64"][field] = value
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "unexpected/missing top-level fields"):
                    peel_codec.validate_table_document(damaged)
        for field in ("seconds", "probes"):
            with self.subTest(missing_direct_diagnostic=field):
                damaged = json.loads(json.dumps(trained))
                del damaged["entries"]["64"][field]
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "unexpected/missing top-level fields"):
                    peel_codec.validate_table_document(damaged)
        for field, value in (
                ("seconds", -0.1),
                ("seconds", 0),
                ("probes", False),
                ("probes", -1)):
            with self.subTest(invalid_direct_diagnostic=field, value=value):
                damaged = json.loads(json.dumps(trained))
                damaged["entries"]["64"][field] = value
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        shipped = peel_codec.make_table_document(
            {64: complete_entry(64, "shipped")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        for mutation in ("missing", "false"):
            with self.subTest(shipped_marker=mutation):
                damaged = json.loads(json.dumps(shipped))
                if mutation == "missing":
                    del damaged["entries"]["64"]["reverted_to_shipped"]
                else:
                    damaged["entries"]["64"]["reverted_to_shipped"] = False
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        sweep = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        sweep_diagnostics = (
            "S", "source_hits", "target_mean", "seconds", "screen",
            "screen_cells", "finals", "rejected", "real_trials",
        )
        for field in sweep_diagnostics:
            with self.subTest(missing_sweep_diagnostic=field):
                damaged = json.loads(json.dumps(sweep))
                del damaged["entries"]["2"][field]
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "unexpected/missing top-level fields"):
                    peel_codec.validate_table_document(damaged)
        for field, value in (
                ("S", False),
                ("source_hits", 99),
                ("target_mean", int(
                    sweep["entries"]["2"]["target_mean"])),
                ("screen", 2999),
                ("rejected", 33),
                ("real_trials", 11)):
            with self.subTest(invalid_sweep_diagnostic=field):
                damaged = json.loads(json.dumps(sweep))
                damaged["entries"]["2"][field] = value
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "invalid proxy-sweep diagnostics"):
                    peel_codec.validate_table_document(damaged)

    def test_generator_settings_are_exact_and_match_entry_receipts(self):
        direct = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        direct_mutations = (
            ("rank trials", lambda document:
             document["provenance"]["settings"].__setitem__(
                 "rank_trials", 11)),
            ("gate context", lambda document:
             document["entries"]["64"]["search_receipt"]["context"].__setitem__(
                 "screen", 2)),
            ("K coverage", lambda document:
             document["provenance"]["settings"].__setitem__("kmax", 65)),
            ("measurement policy", lambda document:
             document["provenance"]["settings"].__setitem__("loss", 0.2)),
            ("unexpected setting", lambda document:
             document["provenance"]["settings"].__setitem__("extra", 1)),
        )
        for label, mutate in direct_mutations:
            with self.subTest(direct_settings=label):
                damaged = json.loads(json.dumps(direct))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        sweep = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        sweep_mutations = (
            ("real trials", lambda document:
             document["provenance"]["settings"].__setitem__(
                 "real_trials_override", 11)),
            ("K coverage", lambda document:
             document["provenance"]["settings"].__setitem__(
                 "proxy_k_ladder", [2, 3])),
            ("fixed budget", lambda document:
             document["entries"]["2"]["search_receipt"]["context"].__setitem__(
                 "finals", 15)),
            ("proxy setting", lambda document:
             document["provenance"]["settings"].__setitem__(
                 "proxy_ordering", "forged")),
        )
        for label, mutate in sweep_mutations:
            with self.subTest(sweep_settings=label):
                damaged = json.loads(json.dumps(sweep))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        prefix = peel_codec.make_table_document(
            {
                2: funnel_entry(2, "shipped"),
                3: funnel_entry(3, "shipped"),
                4: funnel_entry(4, "trained", native_shape=True),
            },
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings((2, 3, 4)),
        )
        damaged = json.loads(json.dumps(prefix))
        del damaged["entries"]["3"]
        damaged["provenance"]["settings"]["proxy_k_ladder"] = [2, 4]
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "exact proxy-sweep settings receipt"):
            peel_codec.validate_table_document(damaged)

    def test_contract_native_and_proxy_integer_receipts_reject_alias_types(self):
        direct = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        damaged = json.loads(json.dumps(direct))
        damaged["provenance"]["target_contract"]["residue_skew"] = False
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "exact dispatch-v1 contract"):
            peel_codec.validate_table_document(damaged)

        damaged = json.loads(json.dumps(direct))
        entry = damaged["entries"]["64"]
        entry["native"]["residue_skew"] = False
        for target in (
                entry["target_receipt"],
                entry["search_receipt"]["target_receipt"],
                entry["search_receipt"]["shipped_control"]["target_receipt"]):
            target["residue_skew"] = "False"
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "non-canonical integer metadata"):
            peel_codec.validate_table_document(damaged)

        sweep = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        proxy_mutations = (
            ("settings bool", lambda document:
             document["provenance"]["settings"]["proxy_measure_regime"].
             __setitem__("solve_block_bytes", False)),
            ("context bool", lambda document:
             document["entries"]["2"]["search_receipt"]["context"][
                 "proxy_measure_regime"].__setitem__("band_tracking_x", True)),
            ("scale bound bool", lambda document:
             document["entries"]["2"]["search_receipt"]["context"][
                 "scale_centi"].__setitem__(0, False)),
        )
        for label, mutate in proxy_mutations:
            with self.subTest(proxy_alias=label):
                damaged = json.loads(json.dumps(sweep))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

    def test_float_schema_fields_reject_integer_aliases(self):
        direct = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        mutations = (
            ("top scale", lambda document:
             document["entries"]["64"].__setitem__("scale", -1)),
            ("search scale", lambda document:
             document["entries"]["64"]["search_receipt"]["coordinates"].
             __setitem__("scale", -1)),
            ("top recovery metric", lambda document:
             document["entries"]["64"].__setitem__("OH50", 0)),
            ("search recovery metric", lambda document:
             document["entries"]["64"]["search_receipt"].
             __setitem__("OH50", 0)),
            ("control recovery metric", lambda document:
             document["entries"]["64"]["search_receipt"]["shipped_control"].
             __setitem__("OH50", 0)),
        )
        for label, mutate in mutations:
            with self.subTest(field=label):
                damaged = json.loads(json.dumps(direct))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        sweep = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        damaged = json.loads(json.dumps(sweep))
        entry = damaged["entries"]["2"]
        off_lattice = math.nextafter(entry["scale"], math.inf)
        entry["scale"] = off_lattice
        entry["search_receipt"]["coordinates"]["scale"] = off_lattice
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "invalid search coordinates"):
            peel_codec.validate_table_document(damaged)

        zero_entry = complete_entry(64)
        zero_entry["decode_mbps"] = 0.0
        zero_entry["goodput"] = 0.0
        zero_entry["search_receipt"]["decode_mbps"] = 0.0
        zero_entry["search_receipt"]["goodput"] = 0.0
        zero_entry["search_receipt"]["shipped_control"]["decode_mbps"] = 0.0
        zero_entry["search_receipt"]["shipped_goodput"] = 0.0
        zero_document = peel_codec.make_table_document(
            {64: zero_entry},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        goodput_mutations = (
            ("top goodput", lambda document:
             document["entries"]["64"].__setitem__("goodput", 0)),
            ("search goodput", lambda document:
             document["entries"]["64"]["search_receipt"].
             __setitem__("goodput", 0)),
            ("control goodput", lambda document:
             document["entries"]["64"]["search_receipt"].
             __setitem__("shipped_goodput", 0)),
        )
        for label, mutate in goodput_mutations:
            with self.subTest(field=label):
                damaged = json.loads(json.dumps(zero_document))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "non-canonical probabilities"):
            peel_codec._validate_pmf([0, 1.0], "test PMF")

        profile = fake_profile()
        trials = 10
        block_bytes = 64
        construction_seed = peel_codec.derive_seed(
            1, "validation", 64, trials, block_bytes, "construction")
        loss_seed = peel_codec.derive_seed(
            2, "validation", 64, trials, block_bytes, "loss")
        candidate_digest = peel_codec.pmf_sha256([0.5, 0.5])
        trained = replace(
            metrics(
                construction_seed, loss_seed, profile=profile,
                pmf_digest=candidate_digest),
            oh_mean=0.0,
        )
        shipped = replace(
            metrics(construction_seed, loss_seed, profile=profile),
            oh_mean=0.0,
        )
        validation_receipt = {
            "verdict": "control",
            "margin_percent": 0.0,
            "trials": trials,
            "block_bytes": block_bytes,
            "scale": -1.0,
            "trained_pmf_sha256": candidate_digest,
            "trained_goodput": trained.goodput(64),
            "shipped_goodput": shipped.goodput(64),
            "trained": trained.as_dict(),
            "shipped": shipped.as_dict(),
        }
        peel_codec._validate_validation_receipt(
            validation_receipt, 64, 1, 2, profile.as_dict(),
            {
                "target_profile": "dispatch-v1",
                "seed_policy": "raw",
                "loss": 0.1,
                "schedule": "iid",
            },
            "test validation receipt",
        )
        validation_aliases = (
            ("margin", "margin_percent", 0),
            ("scale", "scale", -1),
            ("trained goodput", "trained_goodput", 1000),
            ("shipped goodput", "shipped_goodput", 1000),
        )
        for label, field, alias in validation_aliases:
            with self.subTest(validation_field=label):
                damaged = json.loads(json.dumps(validation_receipt))
                damaged[field] = alias
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "non-canonical"):
                    peel_codec._validate_validation_receipt(
                        damaged, 64, 1, 2, profile.as_dict(),
                        {
                            "target_profile": "dispatch-v1",
                            "seed_policy": "raw",
                            "loss": 0.1,
                            "schedule": "iid",
                        },
                        "test validation receipt",
                    )

    def test_integer_schema_fields_reject_float_aliases(self):
        direct = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        mutations = (
            ("provenance seed", lambda document:
             document["provenance"].__setitem__(
                 "construction_seed_base", 1.0)),
            ("generator setting", lambda document:
             document["provenance"]["settings"].__setitem__(
                 "screen", 1.0)),
            ("top K", lambda document:
             document["entries"]["64"].__setitem__("K", 64.0)),
            ("top coordinate", lambda document:
             document["entries"]["64"].__setitem__("p1", 100.0)),
            ("top fail", lambda document:
             document["entries"]["64"].__setitem__("fail", 0.0)),
            ("top probes", lambda document:
             document["entries"]["64"].__setitem__("probes", 1.0)),
            ("native block count", lambda document:
             document["entries"]["64"]["native"].__setitem__(
                 "block_count", 64.0)),
            ("search trials", lambda document:
             document["entries"]["64"]["search_receipt"].__setitem__(
                 "trials", 10.0)),
            ("search coordinate", lambda document:
             document["entries"]["64"]["search_receipt"]["coordinates"].
             __setitem__("p1", 100.0)),
            ("search fail", lambda document:
             document["entries"]["64"]["search_receipt"].__setitem__(
                 "fail", 0.0)),
            ("control fail", lambda document:
             document["entries"]["64"]["search_receipt"][
                 "shipped_control"].__setitem__("fail", 0.0)),
            ("search context", lambda document:
             document["entries"]["64"]["search_receipt"]["context"].
             __setitem__("screen", 1.0)),
        )
        for label, mutate in mutations:
            with self.subTest(field=label):
                damaged = json.loads(json.dumps(direct))
                mutate(damaged)
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(damaged)

    def test_snapshot_hashes_exact_parsed_bytes_and_rejects_duplicate_keys(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "table.json")
            payload = (
                json.dumps(document, indent=2, allow_nan=False) +
                "\n\n").encode("utf-8")
            with open(path, "wb") as handle:
                handle.write(payload)
            parsed, entries, digest = (
                peel_codec.read_table_document_snapshot(path))
            self.assertEqual(parsed, document)
            self.assertEqual(entries[64]["K"], 64)
            self.assertEqual(
                digest, __import__("hashlib").sha256(payload).hexdigest())

            with open(path, "w") as handle:
                handle.write('{"schema":"first","schema":"second"}\n')
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "duplicate JSON object key"):
                peel_codec.read_table_document_snapshot(path)

            with open(path, "w") as handle:
                handle.write('{"schema":1e999}\n')
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "non-finite JSON number"):
                peel_codec.read_table_document_snapshot(path)

    def test_huge_integer_fields_fail_as_measurement_errors(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        huge = 1 << 4096
        huge_decimal = 1 << 20000
        overflow_mutations = (
            ("search scale", lambda entry:
             entry["search_receipt"]["coordinates"].__setitem__(
                 "scale", huge)),
            ("recovery metric", lambda entry:
             entry["search_receipt"].__setitem__("decode_mbps", huge)),
            ("PMF probability", lambda entry:
             entry["search_receipt"]["peel_pmf"].__setitem__(0, huge)),
            ("native mean", lambda entry:
             entry["native"].__setitem__("target_mean", huge)),
        )
        domain_mutations = (
            ("trial count", lambda entry:
             entry["search_receipt"].__setitem__("trials", huge_decimal)),
            ("block bytes", lambda entry:
             entry["search_receipt"].__setitem__(
                 "block_bytes", huge_decimal)),
            ("odd block bytes", lambda entry:
             entry["search_receipt"].__setitem__("block_bytes", 63)),
        )
        for label, mutate in overflow_mutations:
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                mutate(damaged["entries"]["64"])
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "out-of-range numeric value|invalid target mean|"
                        "non-numeric search coordinates|"
                        "non-canonical recovery metrics|"
                        "non-canonical probabilities"):
                    peel_codec.validate_table_document(damaged)
        for label, mutate in domain_mutations:
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                mutate(damaged["entries"]["64"])
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "invalid search metadata"):
                    peel_codec.validate_table_document(damaged)

    def test_stored_native_metadata_rejects_impossible_staircase(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        native = document["entries"]["64"]["native"]
        native["staircase"] = 10 ** 400
        # Without a native-domain bound the expected mean underflows to zero,
        # and the old absolute tolerance accepted this positive subnormal.
        native["target_mean"] = 5e-324
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid native metadata"):
            peel_codec.validate_table_document(document)

    def test_direct_table_rejects_nextafter_native_target_mean(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        native = document["entries"]["64"]["native"]
        native["target_mean"] = math.nextafter(
            native["target_mean"], math.inf)
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid native metadata"):
            peel_codec.validate_table_document(document)

    def test_stored_native_metadata_rejects_forged_dispatch_geometry(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        for label, field, forged_value in (
                ("staircase", "staircase", 11),
                ("source hits", "source_hits", 3)):
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                native = damaged["entries"]["64"]["native"]
                native[field] = forged_value
                forged_mean = (
                    64.0 * min(
                        native["source_hits"], native["staircase"]) /
                    native["staircase"]
                )
                native["target_mean"] = forged_mean
                for target in (
                        damaged["entries"]["64"]["target_receipt"],
                        damaged["entries"]["64"]["search_receipt"][
                            "target_receipt"],
                        damaged["entries"]["64"]["search_receipt"][
                            "shipped_control"]["target_receipt"]):
                    target[field] = str(forged_value)
                    target["target_mean"] = peel_codec._target_mean_spec(
                        64, native["staircase"], native["source_hits"],
                        target["staircase_scale"])
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "invalid native metadata"):
                    peel_codec.validate_table_document(damaged)

    def test_snapshot_rejects_nonregular_input_before_reading_it(self):
        source = mock.MagicMock()
        source.__enter__.return_value = source
        source.fileno.return_value = 7
        fifo_stat = SimpleNamespace(st_mode=__import__("stat").S_IFIFO)
        with mock.patch("builtins.open", return_value=source), mock.patch(
                "peel_codec.os.fstat", return_value=fifo_stat):
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "not a regular file"):
                peel_codec.read_table_document_snapshot("fifo")
        source.read.assert_not_called()

    def test_funnel_schema_accepts_shipped_and_scaled_native_shape_arms(self):
        for label, entry in (
                ("shipped", funnel_entry(2, mode="shipped")),
                ("scaled native shape",
                 funnel_entry(2, mode="trained", native_shape=True))):
            with self.subTest(label=label):
                document = peel_codec.make_table_document(
                    {2: entry},
                    generator="tools/peel_sweep.py",
                    bench=sys.executable,
                    **table_kwargs(),
                    settings=sweep_settings(),
                )
                self.assertEqual(
                    peel_codec.validate_table_document(document)[2]["K"], 2)
        document = peel_codec.make_table_document(
            {2: funnel_entry(2, mode="trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        for field, value in (
                ("proxy_cost_model", "forged"),
                ("scale_centi", [0, 5119]),
                ("gate_block_bytes", 63),
                ("screen_cells", 0x100000001)):
            with self.subTest(funnel_context=field):
                damaged = json.loads(json.dumps(document))
                damaged["entries"]["2"]["search_receipt"]["context"][
                    field] = value
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "funnel context"):
                    peel_codec.validate_table_document(damaged)
        damaged = json.loads(json.dumps(document))
        entry = damaged["entries"]["2"]
        entry["scale"] = 1.2345
        entry["search_receipt"]["coordinates"]["scale"] = 1.2345
        entry["target_receipt"]["staircase_scale"] = "1.2344999999999999"
        entry["search_receipt"]["target_receipt"][
            "staircase_scale"] = "1.2344999999999999"
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid search coordinates"):
            peel_codec.validate_table_document(damaged)

    def test_funnel_schema_rejects_native_pmf_labeled_trained(self):
        entry = funnel_entry(2, mode="trained", native_shape=True)
        entry["search_receipt"]["mode"] = "trained"
        explicit_digest = peel_codec.pmf_sha256(entry["peel_pmf"])
        entry["target_receipt"]["pmf_sha256"] = explicit_digest
        entry["search_receipt"]["target_receipt"][
            "pmf_sha256"] = explicit_digest
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "labels the stock PMF as a trained arm"):
            peel_codec.make_table_document(
                {2: entry},
                generator="tools/peel_sweep.py",
                bench=sys.executable,
                **table_kwargs(),
                settings=sweep_settings(),
            )

    def test_validator_refuses_legacy_table_without_touching_output(self):
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "legacy.json")
            output = os.path.join(directory, "validated.json")
            with open(table, "w") as handle:
                json.dump({"64": {"K": 64, "p1": 100}}, handle)
            with open(output, "w") as handle:
                handle.write("sentinel\n")

            status = peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", "unused",
                "--bb", "64",
                *exact_cli_args(),
            ])
            self.assertEqual(status, 2)
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    def test_params_refuses_legacy_table(self):
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "legacy.json")
            with open(table, "w") as handle:
                json.dump({"64": {"K": 64, "p1": 100}}, handle)
            self.assertEqual(
                peel_params.main([
                    "--table", table, "--K", "64",
                    "--target-profile", "dispatch-v1",
                ]), 2)

    def test_validator_refuses_destructive_in_place_output(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            peel_codec.write_json_atomic(table, document)
            before = peel_codec.file_sha256(table)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", table,
                "--bench", sys.executable,
                "--bb", "64",
                *exact_cli_args(),
            ]), 2)
            self.assertEqual(peel_codec.file_sha256(table), before)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_failed_validation_does_not_publish_partial_table(
            self, probe, unused_profile):
        entries = {
            64: complete_entry(64),
            65: complete_entry(65),
        }
        document = peel_codec.make_table_document(
            entries,
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(64, 65),
        )
        probe.side_effect = [
            metrics(1),
            metrics(2),
            peel_codec.MeasurementError("injected failure"),
        ]
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            output = os.path.join(directory, "output.json")
            peel_codec.write_json_atomic(table, document)
            with open(output, "w") as handle:
                handle.write("sentinel\n")

            status = peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                *exact_cli_args(),
            ])
            self.assertEqual(status, 1)
            first = probe.call_args_list[0].kwargs
            second = probe.call_args_list[1].kwargs
            self.assertEqual(
                first["construction_seed"], second["construction_seed"])
            self.assertEqual(first["loss_seed"], second["loss_seed"])
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count, target_mean=99.0))
    @mock.patch("peel_validate.compare_probe")
    def test_native_profile_drift_does_not_publish(
            self, probe, unused_profile):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            output = os.path.join(directory, "output.json")
            peel_codec.write_json_atomic(table, document)
            with open(output, "w") as handle:
                handle.write("sentinel\n")

            status = peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                *exact_cli_args(),
            ])
            self.assertEqual(status, 1)
            probe.assert_not_called()
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_publishes_complete_paired_receipt(
            self, probe, unused_profile):
        call_count = 0

        def measured(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            return metrics_for_probe(
                args, kwargs,
                decode_mbps=1100.0 if call_count == 1 else 1000.0)

        probe.side_effect = measured
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, document)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                "--margin", "0",
                *exact_cli_args(),
            ]), 0)
            validated, entries = peel_codec.read_table_document(output)
            receipt = entries[64]["validation_receipt"]
            self.assertEqual(receipt["verdict"], "keep")
            self.assertEqual(
                receipt["trained"]["construction_seed"],
                receipt["shipped"]["construction_seed"])
            self.assertEqual(
                receipt["trained"]["loss_seed"],
                receipt["shipped"]["loss_seed"])
            self.assertIn("source_provenance", validated)
            self.assertEqual(entries[64]["decode_mbps"], 1100.0)
            self.assertEqual(entries[64]["goodput"], receipt["trained_goodput"])
            self.assertEqual(
                probe.call_args_list[0].kwargs["peel_weights"],
                document["entries"]["64"]["peel_pmf"])
            self.assertEqual(
                validated["source_provenance"]["document_sha256"],
                peel_codec.file_sha256(table))
            self.assertEqual(
                validated["source_provenance"]["entry_count"], 1)
            self.assertEqual(
                validated["source_provenance"]["selected_entry_count"], 1)
            self.assertEqual(
                validated["source_provenance"]["selected_K"], [64])
            self.assertEqual(
                peel_params.load(output, bench=sys.executable)[64]["K"], 64)
            self.assertEqual(entries[64]["seconds"], 0.0)
            self.assertEqual(entries[64]["probes"], 1)
            damaged = json.loads(json.dumps(validated))
            native = damaged["entries"]["64"]["native"]
            native["target_mean"] = math.nextafter(
                native["target_mean"], math.inf)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "invalid native metadata"):
                peel_codec.validate_table_document(damaged)
            forged_tie = json.loads(json.dumps(validated))
            forged_entry = forged_tie["entries"]["64"]
            forged_validation = forged_entry["validation_receipt"]
            trained = forged_validation["trained"]
            shipped = forged_validation["shipped"]
            for name in (
                    "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
                    "OH_max", "decode_mbps"):
                trained[name] = shipped[name]
            canonical_tie = peel_codec._canonical_goodput(shipped, 64)
            forged_validation["shipped_goodput"] = canonical_tie
            forged_validation["trained_goodput"] = math.nextafter(
                canonical_tie, math.inf)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "inconsistent trained goodput"):
                peel_codec.validate_table_document(forged_tie)
            damaged = json.loads(json.dumps(validated))
            damaged["entries"]["64"]["gain_pct"] = 10
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "stale validation summary"):
                peel_codec.validate_table_document(damaged)
            validation_integer_aliases = (
                ("selected K", lambda document:
                 document["provenance"]["settings"]["selected_K"].
                 __setitem__(0, 64.0)),
                ("validation trials", lambda document:
                 document["entries"]["64"]["validation_receipt"].
                 __setitem__("trials", 10.0)),
                ("trained fail", lambda document:
                 document["entries"]["64"]["validation_receipt"]["trained"].
                 __setitem__("fail", 0.0)),
            )
            for label, mutate in validation_integer_aliases:
                with self.subTest(validation_integer_alias=label):
                    damaged = json.loads(json.dumps(validated))
                    mutate(damaged)
                    with self.assertRaises(peel_codec.MeasurementError):
                        peel_codec.validate_table_document(damaged)
            for field in (
                    "validation_receipt", "verified_mbps", "verified_oh",
                    "shipped_mbps", "gain_pct"):
                with self.subTest(missing_validation_summary=field):
                    damaged = json.loads(json.dumps(validated))
                    del damaged["entries"]["64"][field]
                    with self.assertRaises(peel_codec.MeasurementError):
                        peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            damaged["source_provenance"] = {}
            with self.assertRaises(peel_codec.MeasurementError):
                peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            damaged["source_provenance"]["selected_K"] = [65]
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "settings contradict|source selection contradicts"):
                peel_codec.validate_table_document(damaged)
            for field, value in (
                    ("schedule", []),
                    ("schedule", {}),
                    ("loss", -0.0),
                    ("loss", 0.9900001)):
                with self.subTest(source_measurement_policy=field, value=value):
                    damaged = json.loads(json.dumps(validated))
                    damaged["source_provenance"]["provenance"][
                        "measurement_policy"][field] = value
                    with self.assertRaises(peel_codec.MeasurementError):
                        peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            damaged["source_provenance"]["provenance"]["settings"][
                "rank_trials"] = 11
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "contradicts its direct-search settings"):
                peel_codec.validate_table_document(damaged)
            for name, replacement in (
                    ("source_entry_count", True),
                    ("selected_entry_count", True),
                    ("source_table", None),
                    ("trials", "10"),
                    ("block_bytes", 63),
                    ("margin_percent", False),
                    ("loss", False)):
                with self.subTest(validation_setting=name):
                    damaged = json.loads(json.dumps(validated))
                    damaged["provenance"]["settings"][name] = replacement
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError,
                            "settings contradict"):
                        peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            damaged["provenance"]["settings"]["trials"] = 11
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "receipt contradicts"):
                peel_codec.validate_table_document(damaged)
            for name, value in (
                    ("trials", peel_codec._COMPARE_TRIALS_MAX + 1),
                    ("block_bytes", 63)):
                with self.subTest(validation_field=name):
                    damaged = json.loads(json.dumps(validated))
                    damaged["entries"]["64"]["validation_receipt"][
                        name] = value
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError,
                            "invalid validation metadata"):
                        peel_codec.validate_table_document(damaged)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validated_sweep_retains_nested_unverified_opt_in(
            self, probe, unused_profile):
        call_count = 0

        def measured(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            return metrics_for_probe(
                args, kwargs,
                decode_mbps=1100.0 if call_count == 1 else 1000.0)

        probe.side_effect = measured
        source = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained")},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "sweep.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, source)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                "--margin", "0",
                *exact_cli_args(),
            ]), 0)
            validated, _ = peel_codec.read_table_document(output)
            nested_settings = validated[
                "source_provenance"]["provenance"]["settings"]
            self.assertIs(
                nested_settings["allow_unverified_cost_model"], True)
            del nested_settings["allow_unverified_cost_model"]
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "exact proxy-sweep settings receipt"):
                peel_codec.validate_table_document(validated)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_compares_scaled_native_shape_to_true_shipped(
            self, probe, unused_profile):
        calls = 0

        def measured(*args, **kwargs):
            nonlocal calls
            calls += 1
            return metrics_for_probe(
                args, kwargs,
                decode_mbps=1100.0 if calls == 1 else 1000.0)

        probe.side_effect = measured
        source = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "sweep.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, source)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                "--margin", "0",
                *exact_cli_args(),
            ]), 0)
            _, entries = peel_codec.read_table_document(output)
            entry = entries[2]
            self.assertNotIn("reverted_to_control", entry)
            self.assertNotIn("reverted_to_shipped", entry)
            self.assertEqual(
                entry["validation_receipt"]["verdict"], "keep")
            self.assertEqual(
                entry["validation_receipt"]["trained_pmf_sha256"],
                peel_codec.STOCK_PMF_DIGEST)
            self.assertEqual(len(probe.call_args_list), 2)
            candidate, shipped = probe.call_args_list
            self.assertNotIn("peel_weights", candidate.kwargs)
            self.assertEqual(candidate.kwargs["degree_scale"], 2.0)
            self.assertNotIn("peel_weights", shipped.kwargs)
            self.assertNotIn("degree_scale", shipped.kwargs)
            with open(output, encoding="utf-8") as handle:
                damaged = peel_codec.strict_json_loads(handle.read())
            damaged["entries"]["2"]["scale"] = 1.0
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "coordinates contradict"):
                peel_codec.validate_table_document(damaged)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_scaled_native_shape_loss_reverts_to_true_shipped(
            self, probe, unused_profile):
        calls = 0

        def measured(*args, **kwargs):
            nonlocal calls
            calls += 1
            return metrics_for_probe(
                args, kwargs,
                decode_mbps=900.0 if calls == 1 else 1000.0)

        probe.side_effect = measured
        source = peel_codec.make_table_document(
            {2: funnel_entry(2, "trained", native_shape=True)},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "sweep.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, source)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                "--margin", "0",
                *exact_cli_args(),
            ]), 0)
            _, entries = peel_codec.read_table_document(output)
            entry = entries[2]
            self.assertTrue(entry["reverted_to_shipped"])
            self.assertNotIn("reverted_to_control", entry)
            self.assertEqual(entry["scale"], -1.0)
            receipt = entry["validation_receipt"]
            self.assertEqual(receipt["verdict"], "control")
            self.assertEqual(receipt["scale"], 2.0)
            self.assertEqual(
                receipt["trained"]["target_receipt"]["staircase_scale"], "2")
            self.assertEqual(
                receipt["shipped"]["target_receipt"]["staircase_scale"],
                "unset")
            candidate, shipped = probe.call_args_list
            self.assertEqual(candidate.kwargs["degree_scale"], 2.0)
            self.assertNotIn("degree_scale", shipped.kwargs)

    def test_output_identity_check_fails_closed_on_samefile_error(self):
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "output.json")
            for path in (source, output):
                with open(path, "w") as handle:
                    handle.write("{}\n")
            with mock.patch(
                    "peel_validate.os.path.samefile",
                    side_effect=OSError("injected stat failure")):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "could not establish"):
                    peel_validate.require_distinct_output(source, output)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_revert_uses_fresh_shipped_metrics_and_keeps_search(
            self, probe, unused_profile):
        speeds = iter((900.0, 1000.0))
        probe.side_effect = lambda *args, **kwargs: metrics_for_probe(
            args, kwargs, decode_mbps=next(speeds))
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        original_search = document["entries"]["64"]["search_receipt"]
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, document)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                "--margin", "0",
                *exact_cli_args(),
            ]), 0)
            validated, entries = peel_codec.read_table_document(output)
            entry = entries[64]
            receipt = entry["validation_receipt"]
            self.assertTrue(entry["reverted_to_shipped"])
            self.assertEqual(entry["decode_mbps"], 1000.0)
            self.assertEqual(
                entry["construction_seed"],
                receipt["shipped"]["construction_seed"])
            self.assertEqual(entry["loss_seed"], receipt["shipped"]["loss_seed"])
            self.assertEqual(entry["goodput"], receipt["shipped_goodput"])
            self.assertEqual(entry["peel_pmf"], entry["native"]["pmf"])
            self.assertEqual(entry["search_receipt"], original_search)
            self.assertEqual(
                entry["search_receipt"]["coordinates"]["scale"], -1.0)
            self.assertNotIn("seconds", entry)
            self.assertNotIn("probes", entry)
            damaged = json.loads(json.dumps(validated))
            damaged["entries"]["64"]["probes"] = 1
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "unexpected/missing top-level fields"):
                peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            del damaged["entries"]["64"]["search_would_have_lost_pct"]
            with self.assertRaises(peel_codec.MeasurementError):
                peel_codec.validate_table_document(damaged)
            damaged = json.loads(json.dumps(validated))
            damaged["entries"]["64"]["gain_pct"] = False
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "stale validation summary"):
                peel_codec.validate_table_document(damaged)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count, **unused:
        fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_refuses_a_nondecoding_selected_shipped_arm(
            self, probe, unused_profile):
        calls = 0

        def measured(*args, **kwargs):
            nonlocal calls
            calls += 1
            return metrics_for_probe(
                args, kwargs, fail=1 if calls == 2 else 0)

        probe.side_effect = measured
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "shipped")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            output = os.path.join(directory, "validated.json")
            peel_codec.write_json_atomic(table, document)
            with open(output, "w") as handle:
                handle.write("sentinel\n")
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", output,
                "--bench", sys.executable,
                "--trials", "10",
                "--bb", "64",
                *exact_cli_args(),
            ]), 1)
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    def test_shipped_anchor_prevents_hybrid_interpolation(self):
        table = {
            64: {
                "K": 64, "scale": 12.0, "p1": 150, "tilt": 20,
                "dmax": 32, "absorb": 80,
            },
            128: {
                "K": 128, "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100, "reverted_to_shipped": True,
            },
        }
        result, how = peel_params.params_for_k(96, table)
        self.assertEqual(result, peel_params.shipped_params())
        self.assertEqual(how, "shipped-unmeasured")
        result, how = peel_params.params_for_k(
            96, table, allow_unverified_interpolation=True)
        self.assertEqual(result, peel_params.shipped_params())
        self.assertEqual(how, "shipped-control 64..128")

    def test_unmeasured_trained_parameters_require_experimental_opt_in(self):
        table = {
            64: {
                "K": 64, "scale": 10.0, "p1": 120, "tilt": 10,
                "dmax": 32, "absorb": 80,
            },
            128: {
                "K": 128, "scale": 20.0, "p1": 220, "tilt": 30,
                "dmax": 48, "absorb": 60,
            },
        }
        result, how = peel_params.params_for_k(96, table)
        self.assertEqual(result, peel_params.shipped_params())
        self.assertEqual(how, "shipped-unmeasured")
        self.assertTrue(peel_params.uses_shipped_arm(96, table))

        result, how = peel_params.params_for_k(
            96, table, allow_unverified_interpolation=True)
        self.assertEqual(result, {
            "scale": 15.0, "p1": 170, "tilt": 20,
            "dmax": 40, "absorb": 70,
        })
        self.assertEqual(how, "EXPERIMENTAL interp 64..128")
        self.assertFalse(peel_params.uses_shipped_arm(
            96, table, allow_unverified_interpolation=True))

        result, how = peel_params.params_for_k(32, table)
        self.assertEqual(result, peel_params.shipped_params())
        self.assertEqual(how, "shipped-unmeasured")
        result, how = peel_params.params_for_k(
            32, table, allow_unverified_interpolation=True)
        self.assertEqual(result["scale"], 10.0)
        self.assertEqual(how, "EXPERIMENTAL clamped-low from 64")

    @mock.patch("peel_params.stock_pmf", return_value=[0.5, 0.5])
    def test_trained_unset_scale_still_applies_peel_family(self, unused_stock):
        table = {
            64: {
                "K": 64, "scale": -1.0, "p1": 200, "tilt": 0,
                "dmax": 2, "absorb": 100,
                "peel_pmf": [2.0 / 3.0, 1.0 / 3.0],
            },
        }
        self.assertEqual(
            peel_params.pmf_for_k(
                64, table, "bench", target_profile="dispatch-v1"),
            [2.0 / 3.0, 1.0 / 3.0])

    def test_params_requires_opt_in_for_unverified_proxy_table(self):
        entries = {2: funnel_entry(2, "trained")}
        document = peel_codec.make_table_document(
            entries,
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=sweep_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "proxy.json")
            peel_codec.write_json_atomic(table, document)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "2",
                    "--target-profile", "dispatch-v1",
                ]), 2)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "2",
                    "--target-profile", "dispatch-v1",
                    "--allow-unverified-cost-model",
                ]),
                2)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "2",
                    "--target-profile", "dispatch-v1",
                    "--allow-unvalidated-search",
                ]),
                2)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "2",
                    "--target-profile", "dispatch-v1",
                    "--allow-unverified-cost-model",
                    "--allow-unvalidated-search",
                ]),
                0)

    def test_params_requires_opt_in_for_unvalidated_direct_table(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            **table_kwargs(),
            settings=direct_settings(),
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "direct.json")
            peel_codec.write_json_atomic(table, document)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "unvalidated search result"):
                peel_params.load(table, bench=sys.executable)
            entries = peel_params.load(
                table, bench=sys.executable,
                allow_unvalidated_search=True)
            self.assertEqual(entries[64]["K"], 64)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "64",
                    "--target-profile", "dispatch-v1",
                ]),
                2)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "64",
                    "--target-profile", "dispatch-v1",
                    "--allow-unvalidated-search",
                ]),
                0)


class PeelTableV4Tests(unittest.TestCase):
    def make_document(self, entry, *, generator="tools/peel_direct.py"):
        settings = (
            direct_settings(entry["K"], entry["K"])
            if generator == "tools/peel_direct.py" else
            sweep_settings((entry["K"],)))
        return peel_codec.make_table_document(
            {entry["K"]: entry},
            generator=generator,
            bench=sys.executable,
            settings=settings,
            **table_kwargs())

    def test_v4_rejects_stale_native_compare_protocol(self):
        document = self.make_document(complete_entry(64, "trained"))
        self.assertEqual(
            document["provenance"]["native_compare"],
            "wirehair-v2-bench:compare:wh2-target:v4")
        document["provenance"]["native_compare"] = (
            "wirehair-v2-bench:compare:wh2-target:v3")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "unsupported compare protocol"):
            peel_codec.validate_table_document(document)

    def test_direct_and_funnel_paired_documents_round_trip(self):
        direct = self.make_document(complete_entry(64, "trained"))
        funnel = self.make_document(
            funnel_entry(2, "trained"),
            generator="tools/peel_sweep.py")
        self.assertEqual(sorted(peel_codec.validate_table_document(direct)), [64])
        self.assertEqual(sorted(peel_codec.validate_table_document(funnel)), [2])
        self.assertEqual(
            direct["entries"]["64"]["search_receipt"]["protocol"],
            peel_codec.NATIVE_PAIRED_PROTOCOL)

    def test_params_lists_v4_search_and_validated_solve_rates(self):
        base = {
            "K": 64, "scale": -1.0, "p1": 100, "tilt": 0,
            "dmax": 64, "absorb": 100, "oh_mean": 1.0,
        }
        cases = (
            ("search", {**base, "solve_mbps": 123.5}, "123.5"),
            (
                "validated",
                {**base, "verified_solve_mbps": 456.5},
                "456.5",
            ),
        )
        for label, entry, expected in cases:
            with self.subTest(document=label):
                output = io.StringIO()
                with mock.patch(
                        "peel_params.load", return_value={64: entry}), \
                        redirect_stdout(output):
                    self.assertEqual(
                        peel_params.main([
                            "--table", "unused.json",
                            "--bench", sys.executable,
                            "--target-profile", "dispatch-v1",
                        ]),
                        0)
                self.assertIn(expected, output.getvalue())

    def test_paired_receipt_rejects_bool_float_and_signed_zero_aliases(self):
        document = self.make_document(complete_entry(64, "trained"))
        mutations = [
            ("bool fail", lambda r: r["candidate"].__setitem__("fail", False)),
            ("int mean", lambda r: r["candidate"].__setitem__("oh_mean", 1)),
            ("signed zero", lambda r: r.__setitem__("aa_log_cost_mean", -0.0)),
        ]
        for label, mutate in mutations:
            forged = json.loads(json.dumps(document))
            receipt = forged["entries"]["64"]["search_receipt"][
                "paired_measurement"]
            mutate(receipt)
            with self.subTest(alias=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "aggregate does not match"):
                    peel_codec.validate_table_document(forged)

    def test_v4_rejects_signed_zero_in_outer_receipts_settings_and_pmfs(self):
        scale_zero = self.make_document(
            funnel_entry(2, "scale-only", scale_override=0.0),
            generator="tools/peel_sweep.py")
        mutations = (
            (
                "selected scale",
                lambda entry, document:
                entry.__setitem__("scale", -0.0),
            ),
            (
                "search scale",
                lambda entry, document:
                entry["search_receipt"]["coordinates"].__setitem__(
                    "scale", -0.0),
            ),
            (
                "evaluated scale",
                lambda entry, document:
                entry["search_receipt"][
                    "evaluated_coordinates"].__setitem__("scale", -0.0),
            ),
            (
                "top OH",
                lambda entry, document:
                entry.__setitem__("OH_sd", -0.0),
            ),
            (
                "seconds",
                lambda entry, document:
                entry.__setitem__("seconds", -0.0),
            ),
            (
                "rank margin",
                lambda entry, document:
                document["provenance"]["settings"].__setitem__(
                    "rank_margin", -0.0),
            ),
        )
        for label, mutate in mutations:
            forged = json.loads(json.dumps(scale_zero))
            mutate(forged["entries"]["2"], forged)
            with self.subTest(alias=label):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(forged)

        with_zero = self.make_document(
            funnel_entry(2, "trained", trained_dmax=2),
            generator="tools/peel_sweep.py")
        entry = with_zero["entries"]["2"]
        zero_index = len(entry["peel_pmf"]) - 1
        for target in (
                entry["peel_pmf"],
                entry["search_receipt"]["peel_pmf"],
                entry["search_receipt"]["evaluated_pmf"]):
            target[zero_index] = -0.0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "probabil"):
            peel_codec.validate_table_document(with_zero)

    def test_v4_rejects_noncanonical_provenance_paths(self):
        document = self.make_document(complete_entry(64, "trained"))
        for label, path in (
                ("relative", "relative/bench"),
                ("non-normal", "/tmp/a/../bench"),
                ("double slash", "//tmp/bench"),
                ("NUL", "/tmp/bench\0alias"),
                ("unencodable surrogate", "/tmp/\ud800"),
        ):
            forged = json.loads(json.dumps(document))
            forged["provenance"]["benchmark"]["path"] = path
            with self.subTest(path=label):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "provenance"):
                    peel_codec.validate_table_document(forged)

    def test_direct_settings_replay_screen_peeltiming_native_limits(self):
        document = self.make_document(complete_entry(64, "trained"))
        policy = document["provenance"]["measurement_policy"]

        boundary = direct_settings()
        boundary["paired_warmups"] = 2
        boundary["screen_paired_replicates"] = 9998
        peel_codec._validate_search_generator_settings(
            "tools/peel_direct.py", boundary, policy, "settings")

        over_sum = dict(boundary)
        over_sum["screen_paired_replicates"] = 10000
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "impossible direct peeltiming dimensions"):
            peel_codec._validate_search_generator_settings(
                "tools/peel_direct.py", over_sum, policy, "settings")

        over_work = direct_settings(64000, 64000)
        over_work["screen_paired_replicates"] = 9998
        over_work["paired_inner_reps"] = 1024
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "impossible direct peeltiming dimensions"):
            peel_codec._validate_search_generator_settings(
                "tools/peel_direct.py", over_work, policy, "settings")

        nonmonotone = direct_settings(2, 2048)
        nonmonotone["rank_block_bytes"] = 338484
        nonmonotone["max_overhead"] = 0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "impossible direct peeltiming dimensions"):
            peel_codec._validate_search_generator_settings(
                "tools/peel_direct.py", nonmonotone, policy, "settings")

    def test_direct_rejects_scale_only_and_stock_pmf_trained_aliases(self):
        profile = fake_profile()

        def alias_entry(mode, absorb):
            coordinates = {
                "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": absorb,
            }
            construction_seed = peel_codec.derive_seed(
                1, "direct-search", 64, "rank", 4, 64, "construction")
            loss_seed = peel_codec.derive_seed(
                2, "direct-search", 64, "rank", 4, 64, "loss")
            measurement = paired_fixture_measurement(
                profile, profile.pmf, construction_seed, loss_seed,
                candidate_ns=80)
            stock_control = paired_fixture_measurement(
                profile, profile.pmf, construction_seed, loss_seed,
                candidate_ns=100, started_ns=3_000_000_000,
                finished_ns=4_000_000_000)
            selected = measurement.candidate
            return {
                "K": 64,
                **coordinates,
                **selected.as_dict(),
                "goodput": selected.goodput(64),
                "native": profile.as_dict(),
                "peel_pmf": list(profile.pmf),
                "seconds": 0.0,
                "probes": 1,
                "search_receipt": peel_codec.make_paired_search_receipt(
                    measurement,
                    mode=mode,
                    block_count=64,
                    block_bytes=64,
                    search_kind="direct-real-codec",
                    construction_seed_base=1,
                    loss_seed_base=2,
                    seed_domain="direct-search",
                    coordinates=coordinates,
                    peel_pmf=profile.pmf,
                    evaluated_coordinates=coordinates,
                    evaluated_pmf=profile.pmf,
                    stock_control_measurement=stock_control,
                    context={
                        "warm_start": None,
                        "sampling_seed": peel_codec.derive_seed(
                            1, "direct-search", 64, "sampling"),
                        "screen": 1,
                        "refine": 0,
                        "gate_trials": 1,
                        "gate_block_bytes": 64,
                        "screen_paired_replicates": 4,
                        "rank_top": 1,
                    },
                ),
            }

        for mode, absorb, message in (
                ("scale-only", 100, "scale-only"),
                ("trained", 0, "stock-PMF alias")):
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, message):
                    self.make_document(alias_entry(mode, absorb))

        k2_alias = funnel_entry(2, "scale-only")
        k2_alias["search_receipt"]["mode"] = "trained"
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "realized stock-PMF alias"):
            self.make_document(
                k2_alias,
                generator="tools/peel_sweep.py")

    def test_direct_warm_start_chain_round_trips_and_rejects_tamper(self):
        first = complete_entry(64, "trained")
        second = complete_entry(65, "trained")
        second["search_receipt"]["context"]["warm_start"] = [
            first["search_receipt"]["coordinates"][name]
            for name in ("p1", "tilt", "dmax", "absorb")
        ]
        document = peel_codec.make_table_document(
            {64: first, 65: second},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            settings=direct_settings(64, 65),
            **table_kwargs())
        self.assertEqual(
            sorted(peel_codec.validate_table_document(document)),
            [64, 65])
        forged = json.loads(json.dumps(document))
        forged["entries"]["65"]["search_receipt"]["context"][
            "warm_start"][0] += 1
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "warm-start chain"):
            peel_codec.validate_table_document(forged)

    def test_sweep_warm_start_chain_clamps_to_child_scale_box(self):
        def lightweight(scale, scale_max, warm_start):
            return {
                "search_receipt": {
                    "mode": "trained",
                    "coordinates": {
                        "scale": scale, "p1": 200, "tilt": 10,
                        "dmax": 32, "absorb": 75,
                    },
                    "context": {
                        "scale_centi": [0, scale_max],
                        "warm_start": warm_start,
                    },
                },
            }

        complete = list(peel_codec.PROXY_K_LADDER)
        entries = {
            128: lightweight(70.0, 10240, None),
            96: lightweight(
                60.0, 7680, [7000, 200, 110, 32, 75]),
            # Raw 96->64 scale 6000 is clamped to the K=64 ceiling 5120.
            64: lightweight(
                50.0, 5120, [5120, 200, 110, 32, 75]),
        }
        peel_codec._validate_warm_start_chain(
            entries, "tools/peel_sweep.py", complete)
        entries[64]["search_receipt"]["context"]["warm_start"][0] = 5119
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "warm-start chain"):
            peel_codec._validate_warm_start_chain(
                entries, "tools/peel_sweep.py", complete)

    def test_authenticated_native_stdout_is_replayed_not_trusted(self):
        document = self.make_document(complete_entry(64, "trained"))
        forged = json.loads(json.dumps(document))
        receipt = forged["entries"]["64"]["search_receipt"][
            "paired_measurement"]
        receipt["candidate_log_cost_ci_high"] = -99.0
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "aggregate does not match"):
            peel_codec.validate_table_document(forged)

    def test_trained_receipts_replay_independent_stock_control(self):
        cases = (
            (
                "Direct", "tools/peel_direct.py",
                complete_entry(64, "trained"), "direct-search", 64, 64,
            ),
            (
                "Funnel", "tools/peel_sweep.py",
                funnel_entry(2, "trained"), "funnel-search", 4096, 2,
            ),
        )
        for label, generator, entry, domain, block_bytes, block_count in cases:
            with self.subTest(generator=label):
                document = self.make_document(
                    entry, generator=generator)
                search = document["entries"][str(block_count)][
                    "search_receipt"]
                self.assertEqual(
                    search["stock_control_source"],
                    peel_codec.STOCK_CONTROL_EMBEDDED)
                self.assertIsInstance(
                    search["stock_control_measurement"], dict)

                construction_seed = peel_codec.derive_seed(
                    1, domain, block_count, "rank", 4, block_bytes,
                    "construction")
                loss_seed = peel_codec.derive_seed(
                    2, domain, block_count, "rank", 4, block_bytes, "loss")
                biased = paired_fixture_measurement(
                    fake_profile(block_count), fake_profile(block_count).pmf,
                    construction_seed, loss_seed,
                    block_bytes=block_bytes, candidate_ns=80,
                    started_ns=3_000_000_000,
                    finished_ns=4_000_000_000)
                forged = json.loads(json.dumps(document))
                forged["entries"][str(block_count)]["search_receipt"][
                    "stock_control_measurement"] = biased.as_dict()
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError, "biased"):
                    peel_codec.validate_table_document(forged)

                missing = json.loads(json.dumps(document))
                missing["entries"][str(block_count)]["search_receipt"][
                    "stock_control_measurement"] = None
                with self.assertRaisesRegex(
                        peel_codec.MeasurementError,
                        "independent stock control"):
                    peel_codec.validate_table_document(missing)

                if label == "Direct":
                    old_control = paired_fixture_measurement(
                        fake_profile(), fake_profile().pmf,
                        construction_seed, loss_seed,
                        block_bytes=block_bytes, candidate_ns=100)
                    old_splice = json.loads(json.dumps(document))
                    old_splice["entries"]["64"]["search_receipt"][
                        "stock_control_measurement"] = old_control.as_dict()
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError,
                            "did not run after"):
                        peel_codec.validate_table_document(old_splice)

                    foreign_context = paired_context(2.5, 5.0)
                    foreign_context["bound"]["cpu_model"] = "other-test-cpu"
                    foreign_control = peel_codec.parse_peeltiming_output(
                        peeltiming_stdout(
                            profile=fake_profile(),
                            candidate_pmf=fake_profile().pmf,
                            candidate_ns=100,
                            construction_seed=construction_seed,
                            loss_seed=loss_seed,
                            context=foreign_context,
                            block_bytes=block_bytes,
                            started_ns=3_000_000_000,
                            finished_ns=4_000_000_000),
                        **paired_parse_kwargs(
                            profile=fake_profile(),
                            candidate_pmf=fake_profile().pmf,
                            construction_seed=construction_seed,
                            loss_seed=loss_seed,
                            context=foreign_context,
                            block_bytes=block_bytes))
                    spliced = json.loads(json.dumps(document))
                    spliced["entries"]["64"]["search_receipt"][
                        "stock_control_measurement"] = (
                            foreign_control.as_dict())
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError,
                            "host and sampler context"):
                        peel_codec.validate_table_document(spliced)

                    def replace_traces(lines):
                        semantic = lines[1].split(",")
                        for index, field in enumerate(semantic):
                            if field.startswith("trace_sha256="):
                                unused_name, old = field.split("=", 1)
                                semantic[index] = (
                                    "trace_sha256=" +
                                    hashlib.sha256(
                                        f"splice-semantic:{old}".encode(
                                            "ascii")).hexdigest())
                        lines[1] = ",".join(semantic)
                        trace_column = peel_codec.PEELTIMING_COLUMNS.index(
                            "trace_sha256")
                        mapping = {}
                        for line_index in range(3, len(lines) - 1):
                            fields = lines[line_index].split(",")
                            old = fields[trace_column]
                            fields[trace_column] = mapping.setdefault(
                                old,
                                hashlib.sha256(
                                    f"splice:{old}".encode(
                                        "ascii")).hexdigest())
                            lines[line_index] = ",".join(fields)

                    changed_stdout = authenticated_peeltiming_edit(
                        search["stock_control_measurement"]["native_stdout"],
                        replace_traces)
                    changed_control = peel_codec.parse_peeltiming_output(
                        changed_stdout,
                        **paired_parse_kwargs(
                            profile=fake_profile(),
                            candidate_pmf=fake_profile().pmf,
                            construction_seed=construction_seed,
                            loss_seed=loss_seed,
                            block_bytes=block_bytes,
                            context=search[
                                "stock_control_measurement"]["context"]))
                    trace_splice = json.loads(json.dumps(document))
                    trace_splice["entries"]["64"]["search_receipt"][
                        "stock_control_measurement"] = (
                            changed_control.as_dict())
                    with self.assertRaisesRegex(
                            peel_codec.MeasurementError,
                            "loss trace or identity recovery"):
                        peel_codec.validate_table_document(trace_splice)

    def test_stock_control_replay_binds_exact_shared_failure_result(self):
        profile = fake_profile()
        construction_seed = 123
        loss_seed = 456

        def weak_result(result):
            def mutate(unused_index, row):
                if row["replicate"] != 5:
                    return
                skip_peeltiming_row(row, common_overhead=16)
                row.update({
                    "arm_overhead": -1,
                    "recovery_result": result,
                    "recovery_ok": 0,
                })
            return mutate

        selected = paired_fixture_measurement(
            profile, [0.25, 0.75], construction_seed, loss_seed,
            replicates=6, mutate_row=weak_result(4))
        control = paired_fixture_measurement(
            profile, profile.pmf, construction_seed, loss_seed,
            candidate_ns=100, started_ns=3_000_000_000,
            finished_ns=4_000_000_000, replicates=6,
            mutate_row=weak_result(3))
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "loss trace or identity recovery"):
            peel_codec.validate_paired_stock_control_receipt(
                mode="trained",
                source=peel_codec.STOCK_CONTROL_EMBEDDED,
                receipt=control.as_dict(),
                selected_measurement=selected,
                block_count=64,
                block_bytes=64,
                native=profile.as_dict(),
                measurement_policy={
                    "target_profile": peel_codec.TARGET_PROFILE,
                    "seed_policy": peel_codec.TARGET_SEED_POLICY,
                    "loss": 0.1,
                    "schedule": "iid",
                },
                construction_seed=construction_seed,
                loss_seed=loss_seed,
                label="test stock control")

        def candidate_recovers_identity_does_not(unused_index, row):
            if row["replicate"] != 5:
                return
            skip_peeltiming_row(row, common_overhead=16)
            if row["role"] != "candidate":
                row.update({
                    "arm_overhead": -1,
                    "recovery_result": 4,
                    "recovery_ok": 0,
                })

        improved = paired_fixture_measurement(
            profile, [0.25, 0.75], construction_seed, loss_seed,
            candidate_ns=70, replicates=6,
            mutate_row=candidate_recovers_identity_does_not)
        matching_control = paired_fixture_measurement(
            profile, profile.pmf, construction_seed, loss_seed,
            candidate_ns=100, started_ns=3_000_000_000,
            finished_ns=4_000_000_000, replicates=6,
            mutate_row=weak_result(4))
        self.assertTrue(improved.valid_for_promotion)
        self.assertEqual(improved.recovery_improvements, 1)
        accepted_control = peel_codec.validate_paired_stock_control_receipt(
                mode="trained",
                source=peel_codec.STOCK_CONTROL_EMBEDDED,
                receipt=matching_control.as_dict(),
                selected_measurement=improved,
                block_count=64,
                block_bytes=64,
                native=profile.as_dict(),
                measurement_policy={
                    "target_profile": peel_codec.TARGET_PROFILE,
                    "seed_policy": peel_codec.TARGET_SEED_POLICY,
                    "loss": 0.1,
                    "schedule": "iid",
                },
                construction_seed=construction_seed,
                loss_seed=loss_seed,
                label="test stock control")
        self.assertEqual(accepted_control.as_dict(), matching_control.as_dict())

    def test_shipped_validation_rejects_biased_stock_control(self):
        search = complete_entry(64, "shipped")["search_receipt"]
        profile = fake_profile()
        construction_seed = peel_codec.derive_seed(
            3, "validation", 64, 4, 64, "construction")
        loss_seed = peel_codec.derive_seed(
            4, "validation", 64, 4, 64, "loss")
        biased = paired_fixture_measurement(
            profile, profile.pmf, construction_seed, loss_seed,
            candidate_ns=80)
        receipt = {
            "protocol": peel_codec.NATIVE_PAIRED_PROTOCOL,
            "verdict": "control",
            "selected_arm": "identity",
            "source_mode": "shipped",
            "trials": 4,
            "block_bytes": 64,
            "scale": -1.0,
            "trained_pmf_sha256": search["peel_pmf_sha256"],
            "trained_goodput": biased.candidate.goodput(64),
            "shipped_goodput": biased.identity.goodput(64),
            "paired_measurement": biased.as_dict(),
        }
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "biased"):
            peel_codec._validate_validation_receipt(
                receipt, 64, 3, 4, profile.as_dict(),
                {
                    "target_profile": "dispatch-v1",
                    "seed_policy": "raw",
                    "loss": 0.1,
                    "schedule": "iid",
                },
                search, "validation")

    def test_validated_source_requires_exact_native_paired_protocol(self):
        source_entry = complete_entry(64, "trained")
        source = self.make_document(source_entry)
        profile = fake_profile()
        trials = 4
        block_bytes = 64
        construction_seed = peel_codec.derive_seed(
            3, "validation", 64, trials, block_bytes, "construction")
        loss_seed = peel_codec.derive_seed(
            4, "validation", 64, trials, block_bytes, "loss")
        measurement = paired_fixture_measurement(
            profile, source_entry["peel_pmf"],
            construction_seed, loss_seed,
            block_bytes=block_bytes, candidate_ns=70,
            required_margin=0.02)
        validation_receipt = {
            "protocol": peel_codec.NATIVE_PAIRED_PROTOCOL,
            "verdict": "keep",
            "selected_arm": "candidate",
            "source_mode": "trained",
            "trials": trials,
            "block_bytes": block_bytes,
            "scale": -1.0,
            "trained_pmf_sha256":
                source_entry["search_receipt"]["peel_pmf_sha256"],
            "trained_goodput": measurement.candidate.goodput(64),
            "shipped_goodput": measurement.identity.goodput(64),
            "paired_measurement": measurement.as_dict(),
        }
        solve_gain = round(
            100.0 * (
                measurement.candidate.solve_mbps -
                measurement.identity.solve_mbps) /
            measurement.identity.solve_mbps,
            2)
        validated_entry = {
            **source_entry,
            **measurement.candidate.as_dict(),
            "goodput": measurement.candidate.goodput(64),
            "verified_solve_mbps": measurement.candidate.solve_mbps,
            "verified_oh": measurement.candidate.oh_mean,
            "shipped_solve_mbps": measurement.identity.solve_mbps,
            "solve_gain_pct": solve_gain,
            "validation_receipt": validation_receipt,
        }
        source_hash = "a" * 64
        source_provenance = {
            "schema": source["schema"],
            "document_sha256": source_hash,
            "provenance": source["provenance"],
            "entry_count": 1,
            "selected_entry_count": 1,
            "selected_K": [64],
        }
        settings = {
            "source_table": "/tmp/source.json",
            "source_table_sha256": source_hash,
            "source_entry_count": 1,
            "selected_entry_count": 1,
            "selected_K": [64],
            "paired_replicates": trials,
            "block_bytes": block_bytes,
            "kmax": 64,
            "margin_percent": 2.0,
            "paired_context": "/tmp/paired-context.json",
            "paired_warmups": 0,
            "paired_inner_reps": 1,
            "max_overhead": 16,
            "cache_state": "warm",
            "evict_bytes": 4096,
            "target_profile": "dispatch-v1",
            "seed_policy": "raw",
            "loss": 0.1,
            "schedule": "iid",
        }
        validated = peel_codec.make_table_document(
            {64: validated_entry},
            generator="tools/peel_validate.py",
            bench=sys.executable,
            construction_seed_base=3,
            loss_seed_base=4,
            target_profile="dispatch-v1",
            seed_policy="raw",
            loss=0.1,
            schedule="iid",
            settings=settings,
            source_provenance=source_provenance)
        forged = json.loads(json.dumps(validated))
        forged["source_provenance"]["provenance"]["native_paired"] = "forged"
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "source provenance"):
            peel_codec.validate_table_document(forged)

    def validation_cli(self, source, output, *, replicates=4):
        context_path = source + ".paired-context.json"
        thermal_path = source + ".thermal.csv"
        with open(thermal_path, "w", encoding="ascii") as handle:
            handle.write(thermal_csv(
                thermal_row(0.5), thermal_row(1.5), thermal_row(3.0)
            ).decode("ascii"))
        write_unfinalized_paired_context(context_path, thermal_path)
        return [
            "--bench", sys.executable,
            "--table", source,
            "--out", output,
            "--paired-replicates", str(replicates),
            "--bb", "64",
            "--margin", "2",
            "--paired-context", context_path,
            "--paired-warmups", "0",
            "--max-overhead", "16",
            "--evict-bytes", "4096",
            "--target-profile", "dispatch-v1",
            "--seed-policy", "raw",
            "--schedule", "iid",
            "--loss", "0.1",
            "--construction-seed", "3",
            "--loss-seed", "4",
        ]

    def write_source(self, path, entries):
        block_counts = sorted(entries)
        document = peel_codec.make_table_document(
            entries,
            generator="tools/peel_direct.py",
            bench=sys.executable,
            settings=direct_settings(block_counts[0], block_counts[-1]),
            **table_kwargs())
        peel_codec.write_json_atomic(path, document)
        return document

    def write_sweep_source(self, path, entry):
        document = peel_codec.make_table_document(
            {entry["K"]: entry},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            settings=sweep_settings((entry["K"],)),
            **table_kwargs())
        peel_codec.write_json_atomic(path, document)
        return document

    @staticmethod
    def validation_measurement(args, kwargs, *, candidate_ns=70):
        return paired_fixture_measurement(
            kwargs["native_profile"], args[3],
            kwargs["construction_seed"], kwargs["loss_seed"],
            block_bytes=args[2],
            degree_scale=kwargs.get("degree_scale"),
            candidate_ns=candidate_ns,
            required_margin=kwargs["required_margin"])

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_validation_fallback_is_paired_and_replayable(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))
        paired.side_effect = (
            lambda *args, **kwargs:
            self.validation_measurement(args, kwargs, candidate_ns=110))
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 0)
            document, entries, unused_digest = (
                peel_codec.read_table_document_snapshot(output))
            entry = entries[64]
            self.assertTrue(entry["reverted_to_shipped"])
            self.assertEqual(
                entry["validation_receipt"]["protocol"],
                peel_codec.NATIVE_PAIRED_PROTOCOL)
            self.assertEqual(
                entry["validation_receipt"]["selected_arm"], "identity")
            self.assertEqual(document["provenance"]["generator"],
                             "tools/peel_validate.py")
        self.assertEqual(paired.call_count, 1)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_shipped_source_biased_control_never_publishes(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))
        paired.side_effect = (
            lambda *args, **kwargs:
            self.validation_measurement(args, kwargs, candidate_ns=80))
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "shipped")})
            with open(output, "w", encoding="ascii") as handle:
                handle.write("sentinel")
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 1)
            with open(output, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "sentinel")
        self.assertEqual(paired.call_count, 1)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_near_tie_revert_publishes_canonical_positive_zero(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))

        def near_tie(*args, **kwargs):
            return paired_fixture_measurement(
                kwargs["native_profile"], args[3],
                kwargs["construction_seed"], kwargs["loss_seed"],
                block_bytes=args[2],
                degree_scale=kwargs.get("degree_scale"),
                candidate_ns=20001,
                identity_ns=20000,
                required_margin=kwargs["required_margin"])

        paired.side_effect = near_tie
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 0)
            unused_document, entries, unused_digest = (
                peel_codec.read_table_document_snapshot(output))
            delta = entries[64]["search_would_have_solve_delta_pct"]
            self.assertEqual(delta, 0.0)
            self.assertEqual(math.copysign(1.0, delta), 1.0)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_validation_keeps_promotable_trained_candidate(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))
        paired.side_effect = (
            lambda *args, **kwargs:
            self.validation_measurement(args, kwargs, candidate_ns=70))
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 0)
            unused_document, entries, unused_digest = (
                peel_codec.read_table_document_snapshot(output))
            entry = entries[64]
            self.assertNotIn("reverted_to_shipped", entry)
            self.assertEqual(
                entry["validation_receipt"]["verdict"], "keep")
            self.assertEqual(
                entry["validation_receipt"]["selected_arm"], "candidate")
            self.assertEqual(
                entry["solve_mbps"],
                entry["validation_receipt"]["paired_measurement"][
                    "candidate"]["solve_mbps"])

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_validation_round_trip_keeps_candidate_with_shared_weak_seed(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))

        def shared_weak(unused_index, row):
            if row["replicate"] != 5:
                return
            skip_peeltiming_row(row, common_overhead=16)
            row.update({
                "arm_overhead": -1,
                "recovery_result": 4,
                "recovery_ok": 0,
            })

        def measured(*args, **kwargs):
            return paired_fixture_measurement(
                kwargs["native_profile"], args[3],
                kwargs["construction_seed"], kwargs["loss_seed"],
                block_bytes=args[2],
                degree_scale=kwargs.get("degree_scale"),
                candidate_ns=70,
                required_margin=kwargs["required_margin"],
                replicates=6,
                mutate_row=shared_weak)

        paired.side_effect = measured
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(
                    self.validation_cli(source, output, replicates=6)),
                0)
            document, entries, unused_digest = (
                peel_codec.read_table_document_snapshot(output))
            entry = entries[64]
            measurement = entry["validation_receipt"]["paired_measurement"]
            self.assertEqual(entry["validation_receipt"]["verdict"], "keep")
            self.assertEqual(entry["validation_receipt"]["selected_arm"],
                             "candidate")
            self.assertEqual(measurement["candidate"]["fail"], 1)
            self.assertEqual(measurement["identity"]["fail"], 1)
            self.assertEqual(measurement["both_failures"], 1)
            self.assertEqual(
                peel_codec.validate_table_document(document)[64]["K"], 64)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_scale_only_validation_keeps_or_reverts_from_paired_evidence(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))
        for candidate_ns, expected_verdict in ((70, "keep"), (110, "control")):
            with self.subTest(verdict=expected_verdict), \
                    tempfile.TemporaryDirectory() as directory:
                source = os.path.join(directory, "source.json")
                output = os.path.join(directory, "validated.json")
                self.write_sweep_source(
                    source, funnel_entry(2, "scale-only"))
                paired.side_effect = (
                    lambda *args, _ns=candidate_ns, **kwargs:
                    self.validation_measurement(
                        args, kwargs, candidate_ns=_ns))
                self.assertEqual(
                    peel_validate.main(
                        self.validation_cli(source, output)),
                    0)
                unused_document, entries, unused_digest = (
                    peel_codec.read_table_document_snapshot(output))
                entry = entries[2]
                self.assertEqual(
                    entry["validation_receipt"]["source_mode"],
                    "scale-only")
                self.assertEqual(
                    entry["validation_receipt"]["verdict"],
                    expected_verdict)
                self.assertEqual(
                    entry.get("reverted_to_shipped", False),
                    expected_verdict == "control")

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_validated_scale_and_provenance_paths_reject_aliases(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))
        paired.side_effect = (
            lambda *args, **kwargs:
            self.validation_measurement(args, kwargs, candidate_ns=70))
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_sweep_source(
                source,
                funnel_entry(2, "scale-only", scale_override=0.0))
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 0)
            with open(output, encoding="utf-8") as handle:
                document = peel_codec.strict_json_loads(handle.read())

        mutations = (
            (
                "validation scale",
                lambda value: value["entries"]["2"][
                    "validation_receipt"].__setitem__("scale", -0.0),
            ),
            (
                "source table path",
                lambda value: value["provenance"]["settings"].__setitem__(
                    "source_table", "relative/source.json"),
            ),
            (
                "validation margin",
                lambda value: value["provenance"]["settings"].__setitem__(
                    "margin_percent", -0.0),
            ),
            (
                "benchmark path",
                lambda value: value["provenance"]["benchmark"].__setitem__(
                    "path", "relative/bench"),
            ),
            (
                "source benchmark path",
                lambda value: value["source_provenance"]["provenance"][
                    "benchmark"].__setitem__("path", "relative/bench"),
            ),
        )
        for label, mutate in mutations:
            forged = json.loads(json.dumps(document))
            mutate(forged)
            with self.subTest(alias=label):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(forged)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_failed_identity_control_preserves_existing_output(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))

        def failed_identity(*args, **kwargs):
            measurement = self.validation_measurement(
                args, kwargs, candidate_ns=110)
            return replace(
                measurement,
                identity=replace(measurement.identity, fail=1),
                valid_for_promotion=False)

        paired.side_effect = failed_identity
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            with open(output, "w", encoding="ascii") as handle:
                handle.write("sentinel")
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 1)
            with open(output, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "sentinel")

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_validation_failure_never_publishes_partial_output(
            self, paired, profile):
        profile.side_effect = (
            lambda unused_bench, block_count, **unused:
            fake_profile(block_count))

        def measure(*args, **kwargs):
            if args[1] == 65:
                raise peel_codec.MeasurementError("injected failure")
            return self.validation_measurement(args, kwargs)

        paired.side_effect = measure
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            first = complete_entry(64, "trained")
            second = complete_entry(65, "trained")
            second["search_receipt"]["context"]["warm_start"] = [
                first["search_receipt"]["coordinates"][name]
                for name in ("p1", "tilt", "dmax", "absorb")
            ]
            self.write_source(
                source,
                {
                    64: first,
                    65: second,
                })
            with open(output, "w", encoding="ascii") as handle:
                handle.write("sentinel")
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 1)
            with open(output, encoding="ascii") as handle:
                self.assertEqual(handle.read(), "sentinel")

    def test_validation_refuses_destructive_output_and_params_need_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            original = self.write_source(
                source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, source)), 2)
            current, unused_entries, unused_hash = (
                peel_codec.read_table_document_snapshot(source))
            self.assertEqual(current, original)
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError, "unvalidated"):
                peel_params.load(source, bench=sys.executable)
            loaded = peel_params.load(
                source, bench=sys.executable,
                allow_unvalidated_search=True)
            self.assertEqual(loaded[64]["K"], 64)

    def test_v4_exact_top_level_and_native_types_are_fail_closed(self):
        document = self.make_document(complete_entry(64, "trained"))
        for label, mutation in (
                (
                    "top fail bool",
                    lambda value: value.__setitem__("fail", False),
                ),
                (
                    "native staircase float",
                    lambda value: value["native"].__setitem__(
                        "staircase", float(value["native"]["staircase"])),
                ),
                (
                    "unexpected field",
                    lambda value: value.__setitem__("decode_mbps", 1.0),
                )):
            forged = json.loads(json.dumps(document))
            mutation(forged["entries"]["64"])
            with self.subTest(corruption=label):
                with self.assertRaises(peel_codec.MeasurementError):
                    peel_codec.validate_table_document(forged)

    @mock.patch("peel_validate.stock_profile")
    @mock.patch("peel_validate.paired_probe")
    def test_native_profile_drift_refuses_validation_publication(
            self, paired, profile):
        profile.return_value = replace(fake_profile(), source_hits=3)
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.json")
            output = os.path.join(directory, "validated.json")
            self.write_source(source, {64: complete_entry(64, "trained")})
            self.assertEqual(
                peel_validate.main(self.validation_cli(source, output)), 1)
            self.assertFalse(os.path.exists(output))
        paired.assert_not_called()


class FunnelTests(unittest.TestCase):
    @mock.patch("peel_funnel.stock_profile", return_value=fake_profile())
    def test_proxy_structure_uses_receipted_gf256_only_regime(
            self, unused_profile):
        funnel = peel_funnel.Funnel(SimpleNamespace(
            allow_unverified_cost_model=True,
            bench="bench",
            K=64,
            target_profile="dispatch-v1",
        ))
        self.assertEqual(funnel.struct, "10:10:4")
        self.assertEqual(peel_funnel.PROXY_MEASURE_REGIME["gf16_rows"], 0)
        self.assertEqual(funnel.native_profile.heavy_rows, 12)

    def test_unverified_proxy_requires_explicit_opt_in(self):
        self.assertEqual(peel_funnel.main([
            "--K", "64",
            "--paired-context", "/tmp/context.json",
            *exact_cli_args(),
        ]), 2)

    def legacy_v3_shipped_control_is_always_gated_and_ranked(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5,
            gate_bb=64,
            real_trials=20,
            rank_bb=4096,
            rank_top=1,
        )
        trained_a = [100, 200, 100, 64, 100]
        trained_b = [200, 300, 100, 64, 100]
        calls = []

        def probe(vector, trials, block_bytes, tier, *, stock=False):
            calls.append((vector, trials, block_bytes, tier))
            return metrics(len(calls), decode_mbps=900.0 if vector else 1000.0)

        funnel.real_probe = probe
        ranked, dead = funnel.real_select([trained_a, trained_b, None])
        self.assertFalse(dead)
        shipped_calls = [call for call in calls if call[0] is None]
        self.assertEqual(
            shipped_calls,
            [(None, 5, 64, "gate"), (None, 20, 4096, "rank")])
        self.assertTrue(any(vector is None for _, _, vector in ranked))

    def test_proxy_failure_orders_before_predicted_work(self):
        rows = [
            (1.0, 10.0, "fast-but-fragile"),
            (100.0, 0.0, "slower-but-reliable"),
        ]
        rows.sort(key=peel_funnel.proxy_order)
        self.assertEqual(rows[0][2], "slower-but-reliable")

    def test_funnel_thread_domain_matches_native_essearch(self):
        entry = funnel_entry(2, "trained")
        context = entry["search_receipt"]["context"]
        sampling_seed = peel_codec.derive_seed(
            1, "funnel-search", 2, "sampling")
        for threads in (256, 257, 1024):
            with self.subTest(threads=threads):
                candidate = json.loads(json.dumps(context))
                candidate["threads"] = threads
                peel_codec._validate_search_context(
                    candidate, 2, "unverified-proxy-funnel",
                    entry["native"], sampling_seed, "test context")
        candidate = json.loads(json.dumps(context))
        candidate["threads"] = 1025
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "funnel context"):
            peel_codec._validate_search_context(
                candidate, 2, "unverified-proxy-funnel",
                entry["native"], sampling_seed, "test context")

    def test_empty_proxy_screen_still_measures_and_returns_shipped_control(
            self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.struct = "10:10:4"
        funnel.native_profile = fake_profile()
        funnel.box = peel_funnel.search_box(funnel.native_profile)
        funnel.calls = 0
        funnel.a = SimpleNamespace(
            construction_seed=7,
            screen=2,
            init=None,
            init_frac=0.2,
            screen_cells=1,
            refine=0,
            finals=2,
            gate_bb=64,
            gate_trials=1,
            rank_top=1,
            rank_bb=64,
            paired_replicates=4,
            show=1,
        )
        funnel.lhs = mock.Mock(return_value=[
            [100, 100, 100, 64, 100],
            [200, 200, 100, 64, 100],
        ])
        funnel.measure = mock.Mock(
            return_value=[(None, None), (None, None)])
        control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            candidate_ns=100)

        def select(pool):
            self.assertEqual(pool, [None])
            funnel.gate_arm_count = 1
            funnel.control_evaluated_vector = None
            return [(
                control.identity.goodput(64), control, None,
            )], []

        funnel.real_select = mock.Mock(side_effect=select)
        ranked = funnel.run()
        self.assertEqual(len(ranked), 1)
        self.assertIsNone(ranked[0][2])
        funnel.real_select.assert_called_once_with([None])

    def test_trained_winner_receipt_embeds_fixed_later_stock_control(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.struct = "10:10:4"
        funnel.native_profile = fake_profile()
        funnel.box = peel_funnel.search_box(funnel.native_profile)
        funnel.calls = 0
        funnel.a = SimpleNamespace(
            construction_seed=7,
            screen=1,
            init=None,
            init_frac=0.2,
            screen_cells=1,
            refine=0,
            finals=2,
            gate_bb=64,
            gate_trials=1,
            rank_top=1,
            rank_bb=4096,
            paired_replicates=4,
            show=1,
        )
        trained = [1200, 200, 100, 64, 100]
        trained_pmf = peel_codec.family(
            funnel.native_profile.pmf, 200, 0, 64, 100)
        candidate = paired_fixture_measurement(
            funnel.native_profile, trained_pmf, 1, 2,
            block_bytes=4096)
        control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            block_bytes=4096, candidate_ns=100,
            started_ns=3_000_000_000, finished_ns=4_000_000_000)
        funnel.lhs = mock.Mock(return_value=[trained])
        funnel.measure = mock.Mock(return_value=[(1.0, 0.0)])

        def select(pool):
            self.assertEqual(pool, [trained, None])
            funnel.gate_arm_count = 2
            funnel.control_evaluated_vector = None
            funnel.stock_control_measurement = control
            return [
                (candidate.candidate.goodput(64), candidate, trained),
                (control.identity.goodput(64), control, None),
            ], []

        funnel.real_select = mock.Mock(side_effect=select)
        output = io.StringIO()
        with redirect_stdout(output):
            ranked = funnel.run()
        self.assertEqual(ranked[0][2], trained)
        line = next(
            line for line in output.getvalue().splitlines()
            if line.startswith(peel_funnel.FUNNEL_RESULT_PREFIX))
        receipt = json.loads(
            line[len(peel_funnel.FUNNEL_RESULT_PREFIX):])
        self.assertEqual(
            receipt["stock_control_source"],
            peel_codec.STOCK_CONTROL_EMBEDDED)
        self.assertEqual(
            receipt["stock_control_measurement"], control.as_dict())
        self.assertGreater(
            receipt["stock_control_measurement"]["manifest"][
                "started_monotonic_ns"],
            receipt["paired_measurement"]["evidence"][
                "finished_monotonic_ns"])

    def test_realized_stock_family_alias_canonicalizes_to_scale_only(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5, gate_bb=64, paired_replicates=4,
            rank_bb=4096, rank_top=1)
        alias = [1200, 100, 100, 32, 100]
        canonical = [1200, *peel_funnel.SCALE_ONLY_SHAPE]
        calls = []
        scaled = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            block_bytes=4096, degree_scale=12.0)
        control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            block_bytes=4096, candidate_ns=100,
            started_ns=3_000_000_000, finished_ns=4_000_000_000)

        def probe(vector, unused_trials, unused_bb, tier):
            calls.append((vector, tier))
            if tier == "gate":
                return metrics(1)
            return control if vector is None else scaled

        funnel.real_probe = probe
        ranked, unused_dead = funnel.real_select([alias, None])
        self.assertFalse(any(vector == alias for vector, unused in calls))
        self.assertTrue(any(vector == canonical for vector, unused in calls))
        self.assertEqual(ranked[0][2], canonical)

    def legacy_v3_dedupes_complete_config(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5,
            gate_bb=64,
            real_trials=20,
            rank_bb=4096,
            rank_top=2,
        )
        scaled_native = [1200, 100, 100, 64, 100]
        calls = []

        def probe(vector, trials, block_bytes, tier, *, stock=False):
            calls.append((vector, tier, stock))
            return metrics(
                len(calls),
                decode_mbps=900.0 if vector is not None else 1000.0)

        funnel.real_probe = probe
        ranked, dead = funnel.real_select(
            [scaled_native, list(scaled_native), None])
        self.assertFalse(dead)
        self.assertEqual(
            [call for call in calls if call[0] is not None],
            [
                (scaled_native, "gate", False),
                (scaled_native, "rank", False),
            ])
        self.assertEqual(
            [call for call in calls if call[0] is None],
            [(None, "gate", False), (None, "rank", False)])
        self.assertEqual(
            sum(vector == scaled_native for _, _, vector in ranked), 1)

    def legacy_v3_finalists_share_control_map(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5,
            gate_bb=64,
            real_trials=20,
            rank_bb=4096,
            rank_top=2,
        )
        trained_a = [1200, 200, 100, 64, 100]
        trained_b = [1200, 300, 100, 64, 100]
        calls = []

        def probe(vector, trials, block_bytes, tier, *, stock=False):
            result = metrics(
                len(calls) + 1,
                decode_mbps=1000.0 if stock else 1100.0)
            calls.append((vector, tier, stock, result))
            return result

        funnel.real_probe = probe
        ranked, dead = funnel.real_select([trained_a, trained_b, None])
        self.assertFalse(dead)
        shipped_calls = [
            call for call in calls
            if call[1] == "rank" and call[0] is None
        ]
        self.assertEqual(len(shipped_calls), 1)
        control = shipped_calls[0][3]
        self.assertIs(funnel.rank_controls[tuple(trained_a)], control)
        self.assertIs(funnel.rank_controls[tuple(trained_b)], control)
        self.assertIs(funnel.rank_controls[None], control)
        self.assertFalse(any(
            stock and vector is not None
            for vector, tier, stock, unused in calls
            if tier == "rank"))
        self.assertEqual(
            sum(vector is None for unused_g, unused_r, vector in ranked), 1)

    def legacy_v3_goodput_tie_ordering(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5,
            gate_bb=64,
            real_trials=20,
            rank_bb=4096,
            rank_top=1,
        )
        trained = [1200, 200, 100, 64, 100]

        def probe(vector, trials, block_bytes, tier, *, stock=False):
            return metrics(1, decode_mbps=1000.0)

        funnel.real_probe = probe
        ranked, dead = funnel.real_select([trained, None])
        self.assertFalse(dead)
        self.assertIsNone(ranked[0][2])
        self.assertEqual(ranked[1][2], [1200, 100, 100, 64, 100])
        self.assertEqual(ranked[2][2], trained)
        self.assertEqual(ranked[0][0], ranked[1][0])
        self.assertEqual(ranked[1][0], ranked[2][0])

    def legacy_v3_scale_only_goodput_ordering(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5,
            gate_bb=64,
            real_trials=20,
            rank_bb=4096,
            rank_top=1,
        )
        trained = [1200, 200, 100, 64, 100]

        def probe(vector, trials, block_bytes, tier, *, stock=False):
            if vector is None:
                speed = 900.0
            elif peel_funnel.is_scale_only(vector):
                speed = 1100.0
            else:
                speed = 1000.0
            return metrics(1, decode_mbps=speed)

        funnel.real_probe = probe
        ranked, dead = funnel.real_select([trained, None])
        self.assertFalse(dead)
        self.assertEqual(ranked[0][2], [1200, 100, 100, 64, 100])
        self.assertEqual(ranked[1][2], trained)
        self.assertIsNone(ranked[2][2])

    def test_funnel_orders_only_by_paired_log_ci_and_uses_fixed_control(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5, gate_bb=64, paired_replicates=4,
            rank_bb=4096, rank_top=3)
        a = [1200, 200, 100, 64, 100]
        b = [1200, 300, 100, 64, 100]
        base_a = paired_fixture_measurement(
            funnel.native_profile, [0.25, 0.75], 1, 2)
        base_b = paired_fixture_measurement(
            funnel.native_profile, [0.25, 0.75], 1, 2)
        # A has the better paired CI but deliberately terrible displayed
        # goodput; B has enormous displayed goodput.  CI must still choose A.
        ma = replace(
            base_a,
            candidate=replace(base_a.candidate, solve_mbps=1.0),
            candidate_log_cost_ci_high=-0.30,
            valid_for_promotion=True)
        mb = replace(
            base_b,
            candidate=replace(base_b.candidate, solve_mbps=1e9),
            candidate_log_cost_ci_high=-0.20,
            valid_for_promotion=True)
        control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            candidate_ns=100)

        def probe(vector, unused_trials, unused_bb, tier):
            if tier == "gate":
                return metrics(1)
            if vector is None:
                return control
            if peel_funnel.is_scale_only(vector):
                return replace(control, valid_for_promotion=False)
            return ma if vector == a else mb

        funnel.real_probe = probe
        ranked, unused_dead = funnel.real_select([a, b, None])
        self.assertEqual(ranked[0][2], a)
        self.assertEqual(ranked[1][2], b)
        self.assertIsNone(ranked[-1][2])
        self.assertIs(ranked[-1][1], control)

    def test_failed_predeclared_stock_control_aborts_funnel(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5, gate_bb=64, paired_replicates=4,
            rank_bb=4096, rank_top=1)
        trained = [1200, 200, 100, 64, 100]
        candidate = paired_fixture_measurement(
            funnel.native_profile, [0.25, 0.75], 1, 2)
        control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            candidate_ns=100)
        failed = replace(
            control,
            candidate=replace(control.candidate, fail=1),
            identity=replace(control.identity, fail=1),
            timing_ci_available=False,
            valid_for_promotion=False)

        def probe(vector, unused_trials, unused_bb, tier):
            if tier == "gate":
                return metrics(1)
            return failed if vector is None else candidate

        funnel.real_probe = probe
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "lacks powered timing evidence"):
            funnel.real_select([trained, None])

    def test_biased_predeclared_stock_control_aborts_funnel(self):
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.native_profile = fake_profile()
        funnel.a = SimpleNamespace(
            gate_trials=5, gate_bb=64, paired_replicates=4,
            rank_bb=4096, rank_top=1)
        trained = [1200, 200, 100, 64, 100]
        candidate = paired_fixture_measurement(
            funnel.native_profile, [0.25, 0.75], 1, 2)
        biased_control = paired_fixture_measurement(
            funnel.native_profile, funnel.native_profile.pmf, 1, 2,
            candidate_ns=80)

        def probe(vector, unused_trials, unused_bb, tier):
            if tier == "gate":
                return metrics(1)
            return biased_control if vector is None else candidate

        funnel.real_probe = probe
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "biased"):
            funnel.real_select([trained, None])

    def test_large_k_box_contains_native_density_neighborhood(self):
        box = peel_funnel.search_box(
            fake_profile(64000, target_mean=554.91329479768785))
        self.assertLessEqual(box[0][1], 55491)
        self.assertGreaterEqual(box[0][2], 55492)

    @staticmethod
    def _measure_funnel():
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.struct = "10:10:12"
        funnel.native_profile = replace(
            fake_profile(), heavy_rows=10, dense_rows=12)
        funnel.calls = 0
        funnel.a = SimpleNamespace(
            bench="bench", batch=60, threads=2, cell_base=100)
        return funnel

    @staticmethod
    def _measure_output(funnel, vectors, statuses):
        header = (
            "S,H,D2,scale,shape,p1,p2,p3,c1,c2,c3,c4,c5,dmax,"
            "peel,peel_p1,peel_tilt,peel_dmax,peel_absorb,"
            "cells,failures,fail_rate,solved,pred_ns,xors,muladds,copies,"
            "zerofills,status"
        )
        rows = []
        for vector, status in zip(vectors, statuses):
            rejected = status == "rejected"
            tail = (
                ["2", "2" if rejected else "0",
                 "1.00000000" if rejected else "0.00000000",
                 "0" if rejected else "2",
                 "0.000000" if rejected else "10.500000",
                 "1", "2", "3", "4", status]
            )
            rows.append(",".join(funnel.token(vector).split(",") + tail))
        banner = (
            "# essearch measure,N=64,cells=[100,102),solve_bb=0,"
            "cost_model_bb=1280,cost_model_verified=0,band_tracking_x=1,"
            "loss=0.100000,seed_base=55,completion=mixed,geometry=frozen,"
            "period=244,gf16_rows=0,threads=2"
        )
        return banner + "\n" + header + "\n" + "\n".join(rows) + "\n"

    @mock.patch("peel_funnel.subprocess.run")
    def test_measure_forwards_opt_in_and_filters_rejected_rows(self, run):
        funnel = self._measure_funnel()
        vectors = [
            [100, 100, 100, 64, 100],
            [200, 200, 110, 32, 75],
        ]
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=self._measure_output(
                funnel, vectors, ["ok", "rejected"]),
            stderr="")
        measured = funnel.measure(vectors, 2)
        self.assertEqual(measured[0], (10.5, 0.0))
        self.assertEqual(measured[1], (None, None))
        command = run.call_args.args[0]
        self.assertIn("--allow-unverified-cost-model", command)
        self.assertEqual(command[command.index("--h-hi") + 1], "10")
        self.assertEqual(command[command.index("--d2-hi") + 1], "12")

    @mock.patch("peel_funnel.subprocess.run")
    def test_measure_rejects_reordered_config_receipts(self, run):
        funnel = self._measure_funnel()
        vectors = [
            [100, 100, 100, 64, 100],
            [200, 200, 110, 32, 75],
        ]
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0,
            stdout=self._measure_output(
                funnel, vectors[::-1], ["ok", "ok"]),
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError,
                "echoed configuration mismatch"):
            funnel.measure(vectors, 2)

    @mock.patch("peel_funnel.paired_probe")
    @mock.patch("peel_funnel.compare_probe")
    def test_rank_arms_use_one_paired_seed(self, compare, paired):
        compare.side_effect = lambda *args, **kwargs: metrics_for_probe(
            args, kwargs)
        paired.side_effect = lambda *args, **kwargs: paired_fixture_measurement(
            kwargs["native_profile"], args[3],
            kwargs["construction_seed"], kwargs["loss_seed"],
            block_bytes=args[2],
            degree_scale=kwargs.get("degree_scale"))
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.a = SimpleNamespace(
            construction_seed=7,
            loss_seed=8,
            bench="bench",
            target_profile="dispatch-v1",
            seed_policy="raw",
            loss=0.1,
            schedule="iid",
            paired_warmups=0,
            paired_inner_reps=1,
            max_overhead=16,
            cache_state="warm",
            evict_bytes=4096,
            paired_context="/tmp/context.json",
            rank_margin=0.0,
        )
        funnel.native_profile = fake_profile()
        trained = [1200, 100, 100, 2, 100]
        trained_result = funnel.real_probe(
            trained, 20, 4096, "rank")
        shipped_result = funnel.real_probe(
            None, 20, 4096, "rank")
        gate_result = funnel.real_probe(
            None, 5, 64, "gate")
        self.assertEqual(
            trained_result.manifest["construction_seed_base"],
            shipped_result.manifest["construction_seed_base"])
        self.assertEqual(
            trained_result.manifest["loss_seed_base"],
            shipped_result.manifest["loss_seed_base"])
        self.assertNotEqual(
            trained_result.manifest["construction_seed_base"],
            gate_result.construction_seed)
        self.assertNotEqual(
            trained_result.manifest["loss_seed_base"],
            gate_result.loss_seed)

    @mock.patch("peel_funnel.paired_probe")
    @mock.patch("peel_funnel.compare_probe")
    def test_scale_only_uses_stock_hook_even_when_native_sum_is_inexact(
            self, compare, paired):
        compare.side_effect = lambda *args, **kwargs: metrics_for_probe(
            args, kwargs)
        paired.return_value = paired_fixture_measurement(
            fake_profile(), fake_profile().pmf, 1, 2,
            block_bytes=4096, degree_scale=12.0)
        pmf = [0.2, 0.3] + [0.0] * 61 + [0.4999999999999999]
        self.assertEqual(len(pmf), 64)
        self.assertNotEqual(sum(pmf), 1.0)
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.a = SimpleNamespace(
            construction_seed=7,
            loss_seed=8,
            bench="bench",
            target_profile="dispatch-v1",
            seed_policy="raw",
            loss=0.1,
            schedule="iid",
            paired_warmups=0,
            paired_inner_reps=1,
            max_overhead=16,
            cache_state="warm",
            evict_bytes=4096,
            paired_context="/tmp/context.json",
            rank_margin=0.0,
        )
        funnel.native_profile = replace(fake_profile(), pmf=tuple(pmf))
        funnel.real_probe(
            [1200, 100, 100, 64, 100], 20, 4096, "rank")
        call = paired.call_args
        self.assertNotIn("peel_weights", call.kwargs)
        self.assertEqual(call.kwargs["degree_scale"], 12.0)


class DirectTests(unittest.TestCase):
    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_realized_stock_family_alias_is_skipped_before_probing(
            self, unused_profile):
        args = SimpleNamespace(
            bench="bench", construction_seed=7, loss_seed=8,
            target_profile="dispatch-v1", screen=1, refine=0,
            gate_trials=5, gate_bb=64, screen_paired_replicates=4,
            paired_replicates=4, rank_bb=64, rank_top=1,
        )
        direct = peel_direct.Direct(args)
        alias = [100, 0, 32, 100]
        alias_pmf = peel_codec.family(fake_profile().pmf, *alias)
        self.assertNotEqual(alias_pmf, list(fake_profile().pmf))
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(alias_pmf, 64),
            peel_codec._native_peel_cdf_signature(
                fake_profile().pmf, 64))
        direct.lhs = mock.Mock(return_value=[alias])
        control = paired_fixture_measurement(
            fake_profile(), fake_profile().pmf, 1, 2,
            candidate_ns=100)
        calls = []

        def probe(unused_k, vector, unused_trials, unused_bb, tier):
            calls.append((vector, tier))
            if vector is not None:
                self.fail("realized stock alias reached a native probe")
            return control

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(calls, [(None, "rank")])

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def legacy_v3_shipped_control_recovery_metrics(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            construction_seed=7,
            loss_seed=8,
            target_profile="dispatch-v1",
            screen=2,
            refine=0,
            gate_trials=5,
            gate_bb=64,
            rank_trials=20,
            rank_bb=4096,
            rank_top=2,
        ))
        calls = []

        def probe(block_count, vector, trials, block_bytes, tier):
            calls.append((vector, trials, block_bytes, tier))
            return metrics(
                len(calls), fail=1 if vector is not None else 0)

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(calls[-1], (None, 20, 4096, "rank"))

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def legacy_v3_rank_top_goodput(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            construction_seed=7,
            loss_seed=8,
            target_profile="dispatch-v1",
            screen=1,
            refine=0,
            gate_trials=5,
            gate_bb=64,
            rank_trials=20,
            rank_bb=4096,
            rank_top=1,
        ))
        calls = []

        def probe(block_count, vector, trials, block_bytes, tier):
            calls.append((vector, trials, block_bytes, tier))
            return metrics(
                len(calls),
                decode_mbps=900.0 if vector is not None else 1000.0)

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertIsNotNone(calls[-2][0])
        self.assertIsNone(calls[-1][0])

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def legacy_v3_native_pmf_alias_mock(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            construction_seed=7,
            loss_seed=8,
            target_profile="dispatch-v1",
            screen=2,
            refine=0,
            gate_trials=5,
            gate_bb=64,
            rank_trials=20,
            rank_bb=4096,
            rank_top=2,
        ))
        identity = [100, 0, 64, 100]
        direct.lhs = mock.Mock(
            return_value=[identity, list(identity)])
        calls = []

        def probe(block_count, vector, trials, block_bytes, tier):
            calls.append((vector, tier))
            return metrics(len(calls))

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(calls, [(None, "rank")])

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_compare_goodput_cannot_change_paired_ci_selection(
            self, unused_profile):
        args = SimpleNamespace(
            bench="bench", construction_seed=7, loss_seed=8,
            target_profile="dispatch-v1", screen=2, refine=0,
            gate_trials=5, gate_bb=64, screen_paired_replicates=4,
            paired_replicates=4, rank_bb=64, rank_top=2,
        )
        direct = peel_direct.Direct(args)
        slow = [150, 0, 64, 100]
        fast = [250, 0, 64, 100]
        direct.lhs = mock.Mock(return_value=[slow, fast])
        profile = fake_profile()
        paired = {
            tuple(slow): paired_fixture_measurement(
                profile, peel_codec.family(profile.pmf, *slow), 1, 2,
                candidate_ns=90),
            tuple(fast): paired_fixture_measurement(
                profile, peel_codec.family(profile.pmf, *fast), 1, 2,
                candidate_ns=70),
            None: paired_fixture_measurement(
                profile, profile.pmf, 1, 2, candidate_ns=100),
        }

        def probe(unused_k, vector, unused_trials, unused_bb, tier):
            if tier == "gate":
                # The slower arm has absurdly better compare goodput.  The
                # gate may observe only recovery, so this cannot order it.
                return metrics(
                    1, decode_mbps=1e12 if vector == slow else 1.0)
            return paired[None if vector is None else tuple(vector)]

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertEqual(
            [result[name] for name in ("p1", "tilt", "dmax", "absorb")],
            fast)
        self.assertEqual(
            result["search_receipt"]["stock_control_source"],
            peel_codec.STOCK_CONTROL_EMBEDDED)
        self.assertEqual(
            result["search_receipt"]["stock_control_measurement"],
            paired[None].as_dict())

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_failed_predeclared_stock_control_aborts_direct(
            self, unused_profile):
        args = SimpleNamespace(
            bench="bench", construction_seed=7, loss_seed=8,
            target_profile="dispatch-v1", screen=1, refine=0,
            gate_trials=5, gate_bb=64, screen_paired_replicates=4,
            paired_replicates=4, rank_bb=64, rank_top=1,
        )
        direct = peel_direct.Direct(args)
        vector = [200, 0, 64, 100]
        direct.lhs = mock.Mock(return_value=[vector])
        profile = fake_profile()
        candidate = paired_fixture_measurement(
            profile, peel_codec.family(profile.pmf, *vector), 1, 2)
        control = paired_fixture_measurement(
            profile, profile.pmf, 1, 2, candidate_ns=100)
        failed_control = replace(
            control,
            candidate=replace(control.candidate, fail=1),
            identity=replace(control.identity, fail=1),
            timing_ci_available=False,
            valid_for_promotion=False)

        def probe(unused_k, selected, unused_trials, unused_bb, tier):
            if tier == "gate":
                return metrics(1)
            return failed_control if selected is None else candidate

        direct.probe = probe
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "lacks powered timing evidence"):
            direct.solve(64, None)

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_biased_predeclared_stock_control_aborts_direct(
            self, unused_profile):
        args = SimpleNamespace(
            bench="bench", construction_seed=7, loss_seed=8,
            target_profile="dispatch-v1", screen=1, refine=0,
            gate_trials=5, gate_bb=64, screen_paired_replicates=4,
            paired_replicates=4, rank_bb=64, rank_top=1,
        )
        direct = peel_direct.Direct(args)
        vector = [200, 0, 64, 100]
        direct.lhs = mock.Mock(return_value=[vector])
        profile = fake_profile()
        candidate = paired_fixture_measurement(
            profile, peel_codec.family(profile.pmf, *vector), 1, 2)
        biased_control = paired_fixture_measurement(
            profile, profile.pmf, 1, 2, candidate_ns=80)

        def probe(unused_k, selected, unused_trials, unused_bb, tier):
            if tier == "gate":
                return metrics(1)
            return biased_control if selected is None else candidate

        direct.probe = probe
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "biased"):
            direct.solve(64, None)

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_all_gate_recovery_regressions_use_one_predeclared_stock_rank(
            self, unused_profile):
        args = SimpleNamespace(
            bench="bench", construction_seed=7, loss_seed=8,
            target_profile="dispatch-v1", screen=2, refine=0,
            gate_trials=5, gate_bb=64, screen_paired_replicates=4,
            paired_replicates=4, rank_bb=64, rank_top=2,
        )
        direct = peel_direct.Direct(args)
        direct.lhs = mock.Mock(return_value=[
            [150, 0, 64, 100],
            [250, 0, 64, 100],
        ])
        profile = fake_profile()
        construction_seed = peel_codec.derive_seed(
            7, "direct-search", 64, "rank", 4, 64, "construction")
        loss_seed = peel_codec.derive_seed(
            8, "direct-search", 64, "rank", 4, 64, "loss")
        control = paired_fixture_measurement(
            profile, profile.pmf, construction_seed, loss_seed,
            candidate_ns=100)
        calls = []

        def probe(unused_k, vector, unused_trials, unused_bb, tier):
            calls.append((vector, tier))
            if tier == "gate":
                return metrics(1, fail=1 if vector is None else 2)
            self.assertIsNone(vector)
            return control

        direct.probe = probe
        result = direct.solve(64, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(result["peel_pmf"], list(profile.pmf))
        self.assertEqual(
            [call for call in calls if call[1] == "rank"],
            [(None, "rank")])


class SweepTests(unittest.TestCase):
    def test_optional_log_failures_do_not_mask_the_primary_result(self):
        class BrokenLog:
            def write(self, unused_text):
                raise OSError("write failed")

            def close(self):
                raise OSError("close failed")

        self.assertEqual(
            peel_sweep._close_optional_log(
                BrokenLog(), peel_codec.MeasurementError("measurement failed")),
            [
                "could not append failure note: write failed",
                "could not close log: close failed",
            ])

    def test_shipped_pivot_does_not_create_negative_scale_warm_start(self):
        self.assertIsNone(peel_sweep.warm_start({
            "scale": -1.0,
            "p1": 100,
            "tilt": 0,
            "dmax": 64,
            "absorb": 100,
            "reverted_to_shipped": True,
        }))
        self.assertEqual(peel_sweep.warm_start({
            "scale": 12.5,
            "p1": 200,
            "tilt": 10,
            "dmax": 32,
            "absorb": 75,
        }), [1250, 200, 110, 32, 75])

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.budget", return_value=(1, 0, 2, 1))
    @mock.patch("peel_sweep.subprocess.run")
    def legacy_v2_run_one_recovery_receipt(
            self, run, unused_budget, profile):
        construction_seed = peel_codec.derive_seed(
            7, "funnel-search", 64, "rank", 1, 4096, "construction")
        loss_seed = peel_codec.derive_seed(
            8, "funnel-search", 64, "rank", 1, 4096, "loss")
        native = fake_profile()
        winner_pmf = peel_codec.family(
            native.pmf, 200, 10, 32, 75)
        winner = metrics(
            construction_seed, loss_seed, decode_mbps=100.0,
            profile=native, block_bytes=4096,
            pmf_digest=peel_codec.pmf_sha256(winner_pmf),
            staircase_scale="12.5")
        control = metrics(
            construction_seed, loss_seed, decode_mbps=90.0,
            profile=native, block_bytes=4096)
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "trained",
            "coordinates": {
                "scale": 12.5, "p1": 200, "tilt": 10,
                "dmax": 32, "absorb": 75,
            },
            "peel_pmf": winner_pmf,
            "goodput": 100.0 * 64.0 / 64.25,
            "trials": 1,
            "block_bytes": 4096,
            "rejected": 0,
            "shipped_control": control.as_dict(),
            "shipped_goodput": control.goodput(64),
            **winner.as_dict(),
        }
        run.return_value = subprocess.CompletedProcess(
            ["funnel"], 0,
            stdout=(
                "  finals 2 gated @bb64/1tr, top 1 timed @bb64/1tr "
                "0.0s  0 rejected as non-decoding\n"
                f"{peel_funnel.FUNNEL_RESULT_PREFIX}"
                f"{json.dumps(receipt)}\n"
            ),
            stderr="")
        profile.return_value = native
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 7, 8, 1, True,
            target_profile="dispatch-v1", seed_policy="raw",
            loss=0.1, schedule="iid")
        self.assertIn(
            "--allow-unverified-cost-model", run.call_args.args[0])
        self.assertEqual(result["S"], 10)
        self.assertEqual(result["OH99"], 2.0)
        self.assertEqual(result["construction_seed"], construction_seed)
        self.assertEqual(result["loss_seed"], loss_seed)

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.budget", return_value=(1, 0, 2, 1))
    @mock.patch("peel_sweep.subprocess.run")
    def legacy_v2_run_one_shipped_receipt(
            self, run, unused_budget, profile):
        construction_seed = peel_codec.derive_seed(
            7, "funnel-search", 64, "rank", 1, 4096, "construction")
        loss_seed = peel_codec.derive_seed(
            8, "funnel-search", 64, "rank", 1, 4096, "loss")
        native = fake_profile()
        shipped = metrics(
            construction_seed, loss_seed, profile=native,
            block_bytes=4096)
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "shipped",
            "coordinates": {
                "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100,
            },
            "peel_pmf": list(native.pmf),
            "goodput": shipped.goodput(64),
            "trials": 1,
            "block_bytes": 4096,
            "rejected": 0,
            "shipped_control": shipped.as_dict(),
            "shipped_goodput": shipped.goodput(64),
            **shipped.as_dict(),
        }
        run.return_value = subprocess.CompletedProcess(
            ["funnel"], 0,
            stdout=(
                f"{peel_funnel.FUNNEL_RESULT_PREFIX}"
                f"{json.dumps(receipt)}\n"
            ),
            stderr="")
        profile.return_value = native
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 7, 8, 1, True,
            target_profile="dispatch-v1", seed_policy="raw",
            loss=0.1, schedule="iid")
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(result["scale"], -1.0)

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.budget", return_value=(1, 0, 2, 1))
    @mock.patch("peel_sweep.subprocess.run")
    def test_run_one_replays_complete_paired_funnel_receipt(
            self, run, unused_budget, profile_mock):
        native = fake_profile()
        coordinates = {
            "scale": 12.5, "p1": 200, "tilt": 10,
            "dmax": 32, "absorb": 75,
        }
        candidate_pmf = peel_codec.family(
            native.pmf, 200, 10, 32, 75)
        construction_seed = peel_codec.derive_seed(
            7, "funnel-search", 64, "rank", 4, 4096, "construction")
        loss_seed = peel_codec.derive_seed(
            8, "funnel-search", 64, "rank", 4, 4096, "loss")
        measurement = paired_fixture_measurement(
            native, candidate_pmf, construction_seed, loss_seed,
            block_bytes=4096, degree_scale=12.5)
        stock_control = paired_fixture_measurement(
            native, native.pmf, construction_seed, loss_seed,
            block_bytes=4096, candidate_ns=100,
            started_ns=3_000_000_000, finished_ns=4_000_000_000)
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "trained",
            "coordinates": coordinates,
            "peel_pmf": candidate_pmf,
            "goodput": measurement.candidate.goodput(64),
            "trials": 4,
            "block_bytes": 4096,
            "rejected": 0,
            "evaluated_coordinates": dict(coordinates),
            "evaluated_pmf": candidate_pmf,
            "paired_measurement": measurement.as_dict(),
            "stock_control_source": peel_codec.STOCK_CONTROL_EMBEDDED,
            "stock_control_measurement": stock_control.as_dict(),
            "shipped_goodput": measurement.identity.goodput(64),
            **measurement.candidate.as_dict(),
        }
        run.return_value = subprocess.CompletedProcess(
            ["funnel"], 0,
            stdout=(
                f"{peel_funnel.FUNNEL_RESULT_PREFIX}"
                f"{json.dumps(receipt)}\n"),
            stderr="")
        profile_mock.return_value = native
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 7, 8, 4, True,
            target_profile="dispatch-v1", seed_policy="raw",
            loss=0.1, schedule="iid",
            paired_context="/tmp/context.json",
            paired_warmups=0, paired_inner_reps=1,
            max_overhead=16, cache_state="warm",
            evict_bytes=4096, rank_margin=0.0)
        self.assertEqual(result["solve_mbps"],
                         measurement.candidate.solve_mbps)
        self.assertEqual(
            result["search_receipt"]["paired_measurement"],
            measurement.as_dict())
        self.assertIn("--paired-replicates", run.call_args.args[0])

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.subprocess.run")
    def test_run_one_rejects_trained_coordinates_that_alias_stock(
            self, run, profile_mock):
        native = fake_profile()
        alias_coordinates = {
            "scale": 12.5, "p1": 100, "tilt": 0,
            "dmax": 32, "absorb": 100,
        }
        alias_pmf = peel_codec.family(
            native.pmf,
            alias_coordinates["p1"],
            alias_coordinates["tilt"],
            alias_coordinates["dmax"],
            alias_coordinates["absorb"])
        self.assertNotEqual(alias_pmf, list(native.pmf))
        self.assertEqual(
            peel_codec._native_peel_cdf_signature(alias_pmf, 64),
            peel_codec._native_peel_cdf_signature(native.pmf, 64))
        construction_seed = peel_codec.derive_seed(
            1, "funnel-search", 64, "rank", 4, 4096, "construction")
        loss_seed = peel_codec.derive_seed(
            2, "funnel-search", 64, "rank", 4, 4096, "loss")
        measurement = paired_fixture_measurement(
            native, alias_pmf, construction_seed, loss_seed,
            block_bytes=4096, degree_scale=12.5, candidate_ns=100)
        control = paired_fixture_measurement(
            native, native.pmf, construction_seed, loss_seed,
            block_bytes=4096, candidate_ns=100,
            started_ns=3_000_000_000, finished_ns=4_000_000_000)
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "trained",
            "coordinates": alias_coordinates,
            "peel_pmf": alias_pmf,
            "goodput": measurement.candidate.goodput(64),
            "trials": 4,
            "block_bytes": 4096,
            "rejected": 0,
            "evaluated_coordinates": alias_coordinates,
            "evaluated_pmf": alias_pmf,
            "paired_measurement": measurement.as_dict(),
            "stock_control_source": peel_codec.STOCK_CONTROL_EMBEDDED,
            "stock_control_measurement": control.as_dict(),
            "shipped_goodput": measurement.identity.goodput(64),
            **measurement.candidate.as_dict(),
        }
        run.return_value = subprocess.CompletedProcess(
            ["funnel"], 0,
            stdout=(
                f"{peel_funnel.FUNNEL_RESULT_PREFIX}"
                f"{json.dumps(receipt)}\n"),
            stderr="")
        profile_mock.return_value = native
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "ineligible candidate"):
            peel_sweep.run_one(
                "bench", 64, None, 1, 2, 4, True,
                target_profile="dispatch-v1", seed_policy="raw",
                loss=0.1, schedule="iid",
                paired_context="/tmp/context.json",
                paired_warmups=0, paired_inner_reps=1,
                max_overhead=16, cache_state="warm",
                evict_bytes=4096, rank_margin=0.0)

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.subprocess.run")
    def test_run_one_replays_shipped_paired_funnel_receipt(
            self, run, profile_mock):
        entry = funnel_entry(64, "shipped")
        search = entry["search_receipt"]
        selected_fields = {
            name: entry[name]
            for name in (
                "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
                "OH_max", "solve_mbps")
        }
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "shipped",
            "coordinates": search["coordinates"],
            "peel_pmf": search["peel_pmf"],
            "goodput": entry["goodput"],
            "trials": search["trials"],
            "block_bytes": search["block_bytes"],
            "rejected": 0,
            "evaluated_coordinates": search["evaluated_coordinates"],
            "evaluated_pmf": search["evaluated_pmf"],
            "paired_measurement": search["paired_measurement"],
            "stock_control_source": search["stock_control_source"],
            "stock_control_measurement":
                search["stock_control_measurement"],
            "shipped_goodput": search["shipped_goodput"],
            **selected_fields,
        }
        run.return_value = subprocess.CompletedProcess(
            ["funnel"], 0,
            stdout=(
                f"{peel_funnel.FUNNEL_RESULT_PREFIX}"
                f"{json.dumps(receipt)}\n"),
            stderr="")
        profile_mock.return_value = fake_profile()
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 1, 2, 4, True,
            target_profile="dispatch-v1", seed_policy="raw",
            loss=0.1, schedule="iid",
            paired_context="/tmp/context.json",
            paired_warmups=0, paired_inner_reps=1,
            max_overhead=16, cache_state="warm",
            evict_bytes=4096, rank_margin=0.0)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(result["scale"], -1.0)
        self.assertEqual(result["peel_pmf"], list(fake_profile().pmf))
        self.assertEqual(
            result["search_receipt"]["paired_measurement"],
            search["paired_measurement"])


if __name__ == "__main__":
    unittest.main()
