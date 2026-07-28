#!/usr/bin/env python3
"""Focused correctness tests for the peel-training Python tools."""

import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
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
        staircase_scale="unset"):
    profile = fake_profile()
    if pmf_digest is None:
        pmf_digest = peel_codec.STOCK_PMF_DIGEST
    target = target_receipt(
        profile, construction_seed, loss_seed, 64, pmf_digest,
        staircase_scale)
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
        f"# compare: N=[64,64] trials/bb=10 loss=0.1 "
        f"seed={loss_seed:#x} "
        "max_message_mib=0 schedule=iid "
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
        "rank_trials": 10,
        "rank_block_bytes": 64,
        "rank_top": 1,
        "target_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "loss": 0.1,
        "schedule": "iid",
    }


def sweep_settings(block_counts=(2,), real_trials_override=10):
    return {
        "proxy_k_ladder": list(block_counts),
        "real_trials_override": real_trials_override,
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


def complete_entry(block_count=64, mode="shipped"):
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


def funnel_entry(block_count=64, mode="shipped", native_shape=False):
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


class PeelCodecTests(unittest.TestCase):
    def setUp(self):
        peel_codec._STOCK_CACHE.clear()

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


class PeelTableTests(unittest.TestCase):
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
            *exact_cli_args(),
        ]), 2)

    def test_shipped_control_is_always_gated_and_ranked(self):
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

    def test_funnel_dedupes_complete_config_but_keeps_scaled_native_pmf(self):
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

    def test_finalists_share_one_true_shipped_rank_control(self):
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

    def test_true_shipped_control_wins_exact_tie(self):
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

    def test_scale_only_candidate_is_injected_and_can_win(self):
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

    @mock.patch("peel_funnel.compare_probe")
    def test_rank_arms_use_one_paired_seed(self, compare):
        compare.side_effect = lambda *args, **kwargs: metrics_for_probe(
            args, kwargs)
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
            trained_result.construction_seed,
            shipped_result.construction_seed)
        self.assertEqual(trained_result.loss_seed, shipped_result.loss_seed)
        self.assertNotEqual(
            trained_result.construction_seed,
            gate_result.construction_seed)
        self.assertNotEqual(trained_result.loss_seed, gate_result.loss_seed)

    @mock.patch("peel_funnel.compare_probe")
    def test_scale_only_uses_stock_hook_even_when_native_sum_is_inexact(
            self, compare):
        compare.side_effect = lambda *args, **kwargs: metrics_for_probe(
            args, kwargs)
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
        )
        funnel.native_profile = replace(fake_profile(), pmf=tuple(pmf))
        funnel.real_probe(
            [1200, 100, 100, 64, 100], 20, 4096, "rank")
        call = compare.call_args
        self.assertNotIn("peel_weights", call.kwargs)
        self.assertEqual(call.kwargs["degree_scale"], 12.0)


class DirectTests(unittest.TestCase):
    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_shipped_control_is_ranked_when_all_trained_arms_fail_gate(
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
    def test_rank_top_one_still_compares_trained_with_shipped(
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
    def test_exact_native_pmf_aliases_are_not_measured_as_trained(
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


class SweepTests(unittest.TestCase):
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
    def test_run_one_parses_complete_recovery_receipt(
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
    def test_run_one_accepts_true_shipped_winner(
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


if __name__ == "__main__":
    unittest.main()
