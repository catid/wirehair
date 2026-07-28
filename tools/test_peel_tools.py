#!/usr/bin/env python3
"""Focused correctness tests for the peel-training Python tools."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
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
        seed=0x1122334455667788):
    completion = " ".join(
        f"{name}={value}"
        for name, value in peel_codec.PRODUCTION_COMPLETION_REGIME.items())
    compare_arm = " ".join(
        f"{name}={value}"
        for name, value in peel_codec.PRODUCTION_COMPARE_ARM.items())
    return (
        f"# compare: N=[64,64] trials/bb=10 loss=0.10 seed={seed:#x} "
        "max_message_mib=0 schedule=iid "
        f"schedule_seed={seed:#x} loss_trace=common-id-v2 "
        f"{completion} {compare_arm}\n"
        "codec bb trials fail N_mean OH_mean OH_sd OH50 OH95 OH99 OH_max "
        "create_MBps encode_MBps decode_MBps recover_MBps\n"
        f"baseline 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} "
        f"{oh95} {oh99} {oh_max} 1 2 {decode_mbps} 4\n"
        f"v2 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} "
        f"{oh95} {oh99} {oh_max} 1 2 {decode_mbps} 4\n"
        f"v2_mixed 64 10 {fail} 64 {oh_mean} {oh_sd} {oh50} {oh95} "
        f"{oh99} {oh_max} 1 2 {decode_mbps} 4\n"
    )


def metrics(seed, fail=0, decode_mbps=1000.0):
    return peel_codec.RecoveryMetrics(
        seed=seed,
        fail=fail,
        oh_mean=0.25,
        oh_sd=0.5,
        oh50=0.0,
        oh95=1.0,
        oh99=2.0,
        oh_max=3.0,
        decode_mbps=decode_mbps,
    )


def fake_profile(block_count=64, shipped_mean=None):
    if shipped_mean is None:
        shipped_mean = block_count * 2.0 / 10.0
    return peel_codec.StockProfile(
        block_count=block_count,
        staircase=10,
        source_hits=2,
        shipped_mean=shipped_mean,
        pmf=(1.0 / 64.0,) * 64,
    )


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
    seed = peel_codec.derive_seed(
        1, "direct-search", block_count, "rank", 10, 64)
    measurement = metrics(seed)
    shipped_measurement = metrics(
        seed, decode_mbps=900.0 if not shipped else measurement.decode_mbps)
    return {
        **coordinates,
        **measurement.as_dict(),
        "goodput": measurement.goodput(block_count),
        "native": profile.as_dict(),
        "peel_pmf": selected_pmf,
        "search_receipt": peel_codec.make_search_receipt(
            measurement,
            mode=mode,
            goodput=measurement.goodput(block_count),
            trials=10,
            block_bytes=64,
            search_kind="direct-real-codec",
            base_seed=1,
            seed_domain="direct-search",
            coordinates={
                name: coordinates[name]
                for name in ("scale", "p1", "tilt", "dmax", "absorb")
            },
            peel_pmf=selected_pmf,
            shipped_control=shipped_measurement,
            shipped_goodput=shipped_measurement.goodput(block_count),
            context={
                "sampling_seed": peel_codec.derive_seed(
                    1, "direct-search", block_count, "sampling"),
            },
        ),
    }


def funnel_entry(block_count=64, mode="shipped", native_shape=False):
    """Build one schema-valid funnel receipt, including shipped controls."""
    entry = complete_entry(block_count, mode)
    profile = fake_profile(block_count)
    seed = peel_codec.derive_seed(
        1, "funnel-search", block_count, "rank", 10, 64)
    receipt = entry["search_receipt"]
    entry["seed"] = receipt["seed"] = seed
    receipt["shipped_control"]["seed"] = seed
    receipt["search_kind"] = "unverified-proxy-funnel"
    receipt["seed_domain"] = "funnel-search"
    receipt["context"] = {
        "proxy_cost_model": peel_funnel.PROXY_COST_MODEL,
        "proxy_measure_regime": dict(peel_funnel.PROXY_MEASURE_REGIME),
        "proxy_ordering": peel_funnel.PROXY_ORDERING_PROTOCOL,
        "search_box": peel_funnel.SEARCH_BOX_PROTOCOL,
        "scale_centi": [0, 5120],
        "warm_start": None,
        "sampling_seed": peel_codec.derive_seed(
            1, "funnel-search", block_count, "sampling"),
        "screen": 1,
        "refine": 0,
        "finals": 2,
        "screen_cells": 1,
        "gate_trials": 1,
        "gate_block_bytes": 64,
        "rank_top": 1,
        "threads": 1,
        "batch": 1,
        "cell_base": 0,
    }
    if mode == "trained":
        entry["scale"] = receipt["coordinates"]["scale"] = 12.0
        if native_shape:
            entry["p1"] = receipt["coordinates"]["p1"] = 100
            entry["peel_pmf"] = list(profile.pmf)
            receipt["peel_pmf"] = list(profile.pmf)
            receipt["peel_pmf_sha256"] = peel_codec.pmf_sha256(profile.pmf)
    return entry


class PeelCodecTests(unittest.TestCase):
    def setUp(self):
        peel_codec._STOCK_CACHE.clear()

    @mock.patch("peel_codec.subprocess.run")
    def test_nonzero_compare_exit_fails_closed(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 7, stdout=compare_stdout(), stderr="fatal")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "exited 7.*fatal"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, seed=0x1234)

    @mock.patch.dict(
        os.environ,
        {
            "WIREHAIR_V2_PEEL_DEGREES": "inherited-bad-value",
            "WIREHAIR_V2_OTHER_TEST_HOOK": "also-bad",
            "SAFE_PARENT_VALUE": "kept",
        },
        clear=True)
    @mock.patch("peel_codec.subprocess.run")
    def test_compare_isolates_environment_and_preserves_metrics(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=compare_stdout(), stderr="")
        result = peel_codec.compare_probe(
            "bench", 64, 10, 64,
            peel_weights=[0.25, 0.75],
            degree_scale=2.5,
            seed=0x1122334455667788)

        command = run.call_args.args[0]
        environment = run.call_args.kwargs["env"]
        self.assertEqual(
            command[command.index("--seed") + 1],
            str(0x1122334455667788))
        self.assertNotIn("--mixed-gf16-rows", command)
        self.assertEqual(
            command[command.index("--max-message-mib") + 1], "0")
        self.assertNotIn("WIREHAIR_V2_OTHER_TEST_HOOK", environment)
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
        self.assertEqual(result.seed, 0x1122334455667788)

    @mock.patch("peel_codec.subprocess.run")
    def test_exact_native_pmf_and_metadata_are_parsed(self, run):
        probabilities = [0.125] + [0.875 / 63.0] * 63
        output = [
            "# peelpmf,N=64,degrees=64,staircase=10,"
            "source_hits=2,shipped_mean=12.8",
            "degree,probability",
        ]
        output.extend(
            f"{degree},{probability:.17g}"
            for degree, probability in enumerate(probabilities, 1))
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout="\n".join(output) + "\n", stderr="")

        profile = peel_codec.stock_profile(sys.executable, 64)
        self.assertEqual(profile.block_count, 64)
        self.assertEqual(profile.staircase, 10)
        self.assertEqual(profile.source_hits, 2)
        self.assertEqual(profile.shipped_mean, 12.8)
        self.assertAlmostEqual(profile.pmf[0], 0.125)
        self.assertNotAlmostEqual(profile.pmf[0], 1.0 / 64.0)
        self.assertEqual(
            run.call_args.args[0],
            [sys.executable, "peelpmf", "--N", "64"])
        self.assertFalse(any(
            key.startswith("WIREHAIR_V2_")
            for key in run.call_args.kwargs["env"]))

    def test_pair_seed_is_shared_but_domains_are_separate(self):
        direct_rank = peel_codec.derive_seed(
            7, "direct-search", 64, "rank", 20, 4096)
        self.assertEqual(
            direct_rank,
            peel_codec.derive_seed(
                7, "direct-search", 64, "rank", 20, 4096))
        self.assertNotEqual(
            direct_rank,
            peel_codec.derive_seed(7, "validation", 64, 20, 4096))

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
                peel_weights=[huge, 1.0], seed=1)
        with self.assertRaisesRegex(
                ValueError, "invalid staircase degree scale"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64, degree_scale=huge, seed=1)
        for trials, block_bytes in (
                (peel_codec._COMPARE_TRIALS_MAX + 1, 64),
                (10, peel_codec._COMPARE_BLOCK_BYTES_MAX + 1),
                (10, 63)):
            with self.subTest(trials=trials, block_bytes=block_bytes):
                with self.assertRaisesRegex(
                        ValueError, "invalid compare K, trial count, or block"):
                    peel_codec.compare_probe(
                        "unused", 64, trials, block_bytes, seed=1)

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_boolean_numeric_aliases(self, run):
        with self.assertRaisesRegex(ValueError, "invalid peel weight vector"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64,
                peel_weights=(value for value in (True, 1.0)), seed=1)
        with self.assertRaisesRegex(
                ValueError, "invalid staircase degree scale"):
            peel_codec.compare_probe(
                "unused", 64, 10, 64, degree_scale=True, seed=1)
        run.assert_not_called()

    @mock.patch("peel_codec.subprocess.run")
    def test_native_pmf_rejects_staircase_outside_production_span(self, run):
        staircase = peel_codec._production_staircase_max(64) + 1
        probabilities = [1.0 / 64.0] * 64
        output = [
            f"# peelpmf,N=64,degrees=64,staircase={staircase},"
            f"source_hits=2,shipped_mean={64.0 * 2.0 / staircase:.17g}",
            "degree,probability",
        ]
        output.extend(
            f"{degree},{probability:.17g}"
            for degree, probability in enumerate(probabilities, 1))
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout="\n".join(output) + "\n", stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid peelpmf metadata"):
            peel_codec.stock_profile(sys.executable, 64)

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_semantically_wrong_n_mean(self, run):
        stdout = compare_stdout().replace(
            "v2_mixed 64 10 0 64 ", "v2_mixed 64 10 0 65 ")
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=stdout, stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "semantically wrong compare row"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, seed=0x1122334455667788)

    @mock.patch("peel_codec.subprocess.run")
    def test_compare_rejects_unstructured_trailing_output(self, run):
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout=compare_stdout() + "unexpected junk\n",
            stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "unexpected compare output"):
            peel_codec.compare_probe(
                "bench", 64, 10, 64, seed=0x1122334455667788)

    @mock.patch("peel_codec.subprocess.run")
    def test_native_pmf_rejects_partial_or_trailing_output(self, run):
        probabilities = [1.0 / 64.0] * 64
        output = [
            "# peelpmf,N=64,degrees=64,staircase=10,"
            "source_hits=2,shipped_mean=12.8",
            "degree,probability",
        ]
        output.extend(
            f"{degree},{probability:.17g}"
            for degree, probability in enumerate(probabilities[:-1], 1))
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout="\n".join(output) + "\n", stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "63 of 64 degrees"):
            peel_codec.stock_profile(sys.executable, 64)

        peel_codec._STOCK_CACHE.clear()
        output.append("64,0.015625")
        output.append("unexpected junk")
        run.return_value = subprocess.CompletedProcess(
            ["bench"], 0, stdout="\n".join(output) + "\n", stderr="")
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "unexpected peelpmf output"):
            peel_codec.stock_profile(sys.executable, 64)

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
        lines = [
            "# peelpmf,N=64,degrees=64,staircase=10,"
            "source_hits=2,shipped_mean=12.8",
            "degree,probability",
        ]
        lines.extend(
            f"{degree},{probability:.17g}"
            for degree, probability in enumerate(probabilities, 1))
        run.return_value = "\n".join(lines) + "\n"
        peel_codec.stock_profile("/bench", 64)
        peel_codec.stock_profile("/bench", 64)
        self.assertEqual(run.call_count, 2)

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
                base_seed=1,
                settings={},
                artifact_identity=before,
            )

    def test_benchmark_identity_mismatch_is_refused(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
        )
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "benchmark identity mismatch"):
            peel_codec.verify_benchmark_identity(document, "/bin/true")


class PeelTableTests(unittest.TestCase):
    def test_schema_requires_artifact_native_and_complete_search_receipts(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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

        contradictions = (
            ("forged goodput", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__(
                 "goodput", 123456.0)),
            ("forged rank seed", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__("seed", 5)),
            ("oversized scale", lambda d:
             d["entries"]["64"].__setitem__("scale", 64000.01)),
            ("wrong seed domain", lambda d:
             d["entries"]["64"]["search_receipt"].__setitem__(
                 "seed_domain", "funnel-search")),
            ("wrong sampling seed", lambda d:
             d["entries"]["64"]["search_receipt"]["context"].__setitem__(
                 "sampling_seed", 0)),
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
                 "shipped_mean", 99.0)),
            ("boolean top-level fail alias", lambda d:
             d["entries"]["64"].__setitem__("fail", False)),
            ("boolean top-level quantile alias", lambda d:
             d["entries"]["64"].__setitem__("OH95", True)),
        )
        for label, damage in contradictions:
            with self.subTest(label=label):
                damaged = json.loads(json.dumps(document))
                damage(damaged)
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
            base_seed=1,
            settings={},
        )
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

    def test_snapshot_hashes_exact_parsed_bytes_and_rejects_duplicate_keys(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
            base_seed=1,
            settings={},
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
             entry["native"].__setitem__("shipped_mean", huge)),
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
                        "out-of-range numeric value"):
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
            base_seed=1,
            settings={},
        )
        native = document["entries"]["64"]["native"]
        native["staircase"] = 10 ** 400
        # Without a native-domain bound the expected mean underflows to zero,
        # and the old absolute tolerance accepted this positive subnormal.
        native["shipped_mean"] = 5e-324
        with self.assertRaisesRegex(
                peel_codec.MeasurementError, "invalid native metadata"):
            peel_codec.validate_table_document(document)

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
                ("shipped", funnel_entry(mode="shipped")),
                ("scaled native shape",
                 funnel_entry(mode="trained", native_shape=True))):
            with self.subTest(label=label):
                document = peel_codec.make_table_document(
                    {64: entry},
                    generator="tools/peel_sweep.py",
                    bench=sys.executable,
                    base_seed=1,
                    settings={"allow_unverified_cost_model": True},
                )
                self.assertEqual(
                    peel_codec.validate_table_document(document)[64]["K"], 64)

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
                peel_params.main(["--table", table, "--K", "64"]), 2)

    def test_validator_refuses_destructive_in_place_output(self):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "input.json")
            peel_codec.write_json_atomic(table, document)
            before = peel_codec.file_sha256(table)
            self.assertEqual(peel_validate.main([
                "--table", table,
                "--out", table,
                "--bench", sys.executable,
            ]), 2)
            self.assertEqual(peel_codec.file_sha256(table), before)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count: fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_failed_validation_does_not_publish_partial_table(
            self, probe, unused_profile):
        entries = {
            64: complete_entry(64),
            96: complete_entry(96),
        }
        document = peel_codec.make_table_document(
            entries,
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
            ])
            self.assertEqual(status, 1)
            first_seed = probe.call_args_list[0].kwargs["seed"]
            second_seed = probe.call_args_list[1].kwargs["seed"]
            self.assertEqual(first_seed, second_seed)
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count:
        fake_profile(block_count, shipped_mean=99.0))
    @mock.patch("peel_validate.compare_probe")
    def test_native_profile_drift_does_not_publish(
            self, probe, unused_profile):
        document = peel_codec.make_table_document(
            {64: complete_entry(64)},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
            ])
            self.assertEqual(status, 1)
            probe.assert_not_called()
            with open(output) as handle:
                self.assertEqual(handle.read(), "sentinel\n")

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count: fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_publishes_complete_paired_receipt(
            self, probe, unused_profile):
        call_count = 0

        def measured(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            return metrics(
                kwargs["seed"],
                decode_mbps=1100.0 if call_count == 1 else 1000.0)

        probe.side_effect = measured
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
                "--seed", "9",
            ]), 0)
            validated, entries = peel_codec.read_table_document(output)
            receipt = entries[64]["validation_receipt"]
            self.assertEqual(receipt["verdict"], "keep")
            self.assertEqual(
                receipt["trained"]["seed"], receipt["shipped"]["seed"])
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
        side_effect=lambda unused_bench, block_count: fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validated_sweep_retains_nested_unverified_opt_in(
            self, probe, unused_profile):
        call_count = 0

        def measured(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            return metrics(
                kwargs["seed"],
                decode_mbps=1100.0 if call_count == 1 else 1000.0)

        probe.side_effect = measured
        source = peel_codec.make_table_document(
            {64: funnel_entry(64, "trained")},
            generator="tools/peel_sweep.py",
            bench=sys.executable,
            base_seed=1,
            settings={"allow_unverified_cost_model": True},
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
                "--seed", "9",
            ]), 0)
            validated, _ = peel_codec.read_table_document(output)
            nested_settings = validated[
                "source_provenance"]["provenance"]["settings"]
            self.assertIs(
                nested_settings["allow_unverified_cost_model"], True)
            del nested_settings["allow_unverified_cost_model"]
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "missing its unverified-cost opt-in"):
                peel_codec.validate_table_document(validated)

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
        side_effect=lambda unused_bench, block_count: fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_revert_uses_fresh_shipped_metrics_and_keeps_search(
            self, probe, unused_profile):
        speeds = iter((900.0, 1000.0))
        probe.side_effect = lambda *args, **kwargs: metrics(
            kwargs["seed"], decode_mbps=next(speeds))
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "trained")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
                "--seed", "9",
            ]), 0)
            validated, entries = peel_codec.read_table_document(output)
            entry = entries[64]
            receipt = entry["validation_receipt"]
            self.assertTrue(entry["reverted_to_shipped"])
            self.assertEqual(entry["decode_mbps"], 1000.0)
            self.assertEqual(entry["seed"], receipt["shipped"]["seed"])
            self.assertEqual(entry["goodput"], receipt["shipped_goodput"])
            self.assertEqual(entry["peel_pmf"], entry["native"]["pmf"])
            self.assertEqual(entry["search_receipt"], original_search)
            self.assertEqual(
                entry["search_receipt"]["coordinates"]["scale"], -1.0)
            damaged = json.loads(json.dumps(validated))
            damaged["entries"]["64"]["gain_pct"] = False
            with self.assertRaisesRegex(
                    peel_codec.MeasurementError,
                    "stale validation summary"):
                peel_codec.validate_table_document(damaged)

    @mock.patch(
        "peel_validate.stock_profile",
        side_effect=lambda unused_bench, block_count: fake_profile(block_count))
    @mock.patch("peel_validate.compare_probe")
    def test_validation_refuses_a_nondecoding_selected_shipped_arm(
            self, probe, unused_profile):
        calls = 0

        def measured(*args, **kwargs):
            nonlocal calls
            calls += 1
            return metrics(kwargs["seed"], fail=1 if calls == 2 else 0)

        probe.side_effect = measured
        document = peel_codec.make_table_document(
            {64: complete_entry(64, "shipped")},
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={},
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
                "--seed", "9",
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
            peel_params.pmf_for_k(64, table, "bench"),
            [2.0 / 3.0, 1.0 / 3.0])

    def test_params_requires_opt_in_for_unverified_proxy_table(self):
        entries = {64: complete_entry(64, "trained")}
        document = peel_codec.make_table_document(
            entries,
            generator="tools/peel_direct.py",
            bench=sys.executable,
            base_seed=1,
            settings={"allow_unverified_cost_model": True},
        )
        with tempfile.TemporaryDirectory() as directory:
            table = os.path.join(directory, "proxy.json")
            peel_codec.write_json_atomic(table, document)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "64",
                ]), 2)
            self.assertEqual(
                peel_params.main([
                    "--table", table,
                    "--bench", sys.executable,
                    "--K", "64",
                    "--allow-unverified-cost-model",
                ]),
                0)


class FunnelTests(unittest.TestCase):
    def test_unverified_proxy_requires_explicit_opt_in(self):
        self.assertEqual(peel_funnel.main(["--K", "64"]), 2)

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

        def probe(vector, trials, block_bytes, tier):
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

        def probe(vector, trials, block_bytes, tier):
            calls.append((vector, tier))
            return metrics(
                len(calls),
                decode_mbps=900.0 if vector is not None else 1000.0)

        funnel.real_probe = probe
        ranked, dead = funnel.real_select(
            [scaled_native, list(scaled_native), None])
        self.assertFalse(dead)
        self.assertEqual(
            [call for call in calls if call[0] is not None],
            [(scaled_native, "gate"), (scaled_native, "rank")])
        self.assertEqual(
            [call for call in calls if call[0] is None],
            [(None, "gate"), (None, "rank")])
        self.assertTrue(any(vector is scaled_native
                            for _, _, vector in ranked))

    def test_large_k_box_contains_native_density_neighborhood(self):
        box = peel_funnel.search_box(
            fake_profile(64000, shipped_mean=554.91329479768785))
        self.assertLessEqual(box[0][1], 55491)
        self.assertGreaterEqual(box[0][2], 55492)

    @staticmethod
    def _measure_funnel():
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.struct = "10:10:12"
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
        compare.side_effect = lambda *args, **kwargs: metrics(kwargs["seed"])
        funnel = peel_funnel.Funnel.__new__(peel_funnel.Funnel)
        funnel.k = 64
        funnel.a = SimpleNamespace(seed=7, bench="bench")
        funnel.native_profile = fake_profile()
        trained = [1200, 100, 100, 2, 100]
        trained_result = funnel.real_probe(
            trained, 20, 4096, "rank")
        shipped_result = funnel.real_probe(
            None, 20, 4096, "rank")
        gate_result = funnel.real_probe(
            None, 5, 64, "gate")
        self.assertEqual(trained_result.seed, shipped_result.seed)
        self.assertNotEqual(trained_result.seed, gate_result.seed)


class DirectTests(unittest.TestCase):
    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_shipped_control_is_ranked_when_all_trained_arms_fail_gate(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            seed=7,
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
        result = direct.solve(64, 7, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(calls[-1], (None, 20, 4096, "rank"))

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_rank_top_one_still_compares_trained_with_shipped(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            seed=7,
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
        result = direct.solve(64, 7, None)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertIsNotNone(calls[-2][0])
        self.assertIsNone(calls[-1][0])

    @mock.patch("peel_direct.stock_profile", return_value=fake_profile())
    def test_exact_native_pmf_aliases_are_not_measured_as_trained(
            self, unused_profile):
        direct = peel_direct.Direct(SimpleNamespace(
            bench="bench",
            seed=7,
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
        result = direct.solve(64, 7, None)
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
        rank_seed = peel_codec.derive_seed(
            7, "funnel-search", 64, "rank", 1, 4096)
        winner_pmf = peel_codec.family(
            fake_profile().pmf, 200, 10, 32, 75)
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
            "shipped_control": metrics(
                rank_seed, decode_mbps=90.0).as_dict(),
            "shipped_goodput": metrics(
                rank_seed, decode_mbps=90.0).goodput(64),
            **metrics(rank_seed, decode_mbps=100.0).as_dict(),
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
        profile.return_value = fake_profile()
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 7, 1, True)
        self.assertIn(
            "--allow-unverified-cost-model", run.call_args.args[0])
        self.assertEqual(result["S"], 10)
        self.assertEqual(result["OH99"], 2.0)
        self.assertEqual(result["seed"], rank_seed)

    @mock.patch("peel_sweep.stock_profile")
    @mock.patch("peel_sweep.budget", return_value=(1, 0, 2, 1))
    @mock.patch("peel_sweep.subprocess.run")
    def test_run_one_accepts_true_shipped_winner(
            self, run, unused_budget, profile):
        rank_seed = peel_codec.derive_seed(
            7, "funnel-search", 64, "rank", 1, 4096)
        shipped = metrics(rank_seed)
        receipt = {
            "schema": peel_funnel.FUNNEL_RESULT_SCHEMA,
            "K": 64,
            "mode": "shipped",
            "coordinates": {
                "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100,
            },
            "peel_pmf": list(fake_profile().pmf),
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
        profile.return_value = fake_profile()
        result, _, _ = peel_sweep.run_one(
            "bench", 64, None, 7, 1, True)
        self.assertTrue(result["reverted_to_shipped"])
        self.assertEqual(result["scale"], -1.0)


if __name__ == "__main__":
    unittest.main()
