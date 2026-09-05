#!/usr/bin/env python3
"""Synthetic protocol tests. No timed worker or fixed namespace is launched."""
import copy
import math
import os
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2PublicOffsetR0 as c
from test_Wh2PublicBorrowedCurrentVsWh1R1 import artifact, health_receipt


def fixture():
    receipt = {"baseline_path": str(c.BASELINE), "baseline_sha256": c.BASELINE_HASH,
               "candidate_path": "/synthetic/candidate.so", "candidate_sha256": "c" * 64}
    target = c.h.synthetic_config(c.h.synthetic_expected(), "a" * 64)["target_identity"]
    config = {"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".config.v1",
              "description_sha256": c.description_hash(), "target_identity": target,
              "bindings": {}, "cells": []}
    for arm, name, base in (("B", "baseline", 0x70000000), ("C", "candidate", 0x80000000)):
        path = receipt[name + "_path"]
        symbols = {name: {"path": path, "address": hex(base + 0x485c0 + index * 128),
                         "elf_offset": hex(0x485c0 + index * 128)}
                   for index, name in enumerate(("wirehair_init_", "wirehair_v2_encoder_create_with_options",
                                                 "wirehair_v2_encode", "wirehair_v2_free"))}
        config["bindings"][arm] = {"path": path, "sha256": receipt[name + "_sha256"], "symbols": symbols,
            "memcpy": {"path": str(c.LIBC), "address": "0x901a14c0", "elf_offset": "0x1a14c0"}}
    for k, b in c.CELLS:
        def oracle(tail):
            raw = b"WHV2" + struct.pack("<HHQQIB3s", 1, 32, 0x4b295bbb47f4f9c9, (k - 1) * b + tail, b, 0, bytes(3))
            return {"profile_hex": raw.hex(), "profile_sha256": c.digest(raw),
                    "source_sha256": c.h.SOURCE_SHA256[k, b, tail],
                    "repair_sha256": "d" * 64, "high_id_sha256": "e" * 64}
        cell = dict(oracle(b), K=k, block_bytes=b, partial=[])
        for tail in (1, b - 1):
            cell["partial"].append(dict(oracle(tail), tail_bytes=tail,
                                        systematic_sha256=c.h.SOURCE_SHA256[k, b, tail]))
        config["cells"].append(cell)
    rows = [config]
    for sequence, (rep, cell, comparison, condition) in enumerate(c.roster()):
        k, b = c.CELLS[cell]
        arena, span = 0x10000000, ((k * b + 8191) // 4096) * 4096
        order = "ABBA" if condition % 2 == 0 else "BAAB"
        addresses = {"arena": hex(arena), "arena_bytes": 2 * span + 1024 + 32768 + 64,
            "span": span, "source": hex(arena + 2048), "output": hex(arena + span + 2336),
            "counters": hex(arena + 2 * span + 1024),
            "handles": ["0x30000000", "0x40000000", "0x50000000", "0x60000000"]}
        slots = []
        for side in c.h.sides_for(order):
            physical = int(side == "right") ^ (condition // 2)
            # Both cross-pairs are genuine10% candidate wins. Mapping changes
            # which logical label receives that time, not which implementation.
            elapsed = 900000 if comparison >= 2 and physical == 0 else 1000000
            slots.append(dict(elapsed_ns=elapsed, logical_lane=side, physical_lane=physical))
        rows.append({"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".panel.v1",
            "description_sha256": c.description_hash(), "sequence": sequence,
            "replicate": rep, "cell_index": cell, "cell_key_sha256": c.cell_key(rep, cell),
            "comparison_index": comparison, "condition": condition, "order": order,
            "mapping": "direct" if condition < 2 else "swapped", "K": k, "block_bytes": b,
            "scope_invocations_per_batch": 8192 // k, "addresses": addresses,
            "source_sha256": c.h.SOURCE_SHA256[k, b, b], "output_sha256": c.h.SOURCE_SHA256[k, b, b],
            "profile_sha256": config["cells"][cell]["profile_sha256"],
            "source_immutable_verified": True, "slots": slots})
    return rows + [c.terminal_value()], receipt


def encoded(rows):
    return b"".join(c.canonical(row) + b"\n" for row in rows)


def evidence_fixture():
    """Self-contained receipts; only their baseline trust-anchor hash is mocked.

    No validator is replaced and no real campaign, compiler, git or host probe
    supplies a fixture. The unpatched anchor rejection has a separate test.
    """
    rows, paths = fixture()
    paths["candidate_sha256"] = c.CANDIDATE_HASH
    rows[0]["bindings"]["C"]["sha256"] = c.CANDIDATE_HASH
    commit, compiler, worker = "1" * 40, "/synthetic/g++-13", "/synthetic/worker"
    base_inputs = {path: c.digest(path.encode()) for path in (
        "CMakeLists.txt", "include/wirehair/wirehair.h", "codec/WirehairV2Profile.cpp",
        "codec/V2BorrowedSourceTest.cpp", "abi/wirehair.map")}
    library_commands = compiler + " -O3 -c codec/WirehairV2Profile.cpp\n"
    basis = {"library": str(c.BASELINE), "library_sha256": c.BASELINE_HASH,
             "source_files": base_inputs, "library_commands": library_commands}
    build = dict(paths, campaign=c.CAMPAIGN, source_commit=commit, baseline_build=basis,
        libc_sha256=c.LIBC_HASH, compiler_path=compiler, compiler_sha256="a" * 64,
        compiler_version="g++ (synthetic) 13.3.0\n", worker_path=worker,
        worker_sha256="b" * 64, source_files=dict(base_inputs), build_output="",
        dynamic="(NEEDED) Shared library: [libc.so.6]\n", ldd="libc.so.6 => /synthetic/libc.so.6\n",
        candidate_build={"directory": "/synthetic", "cache_sha256": "d" * 64,
            "commands": library_commands, "dry_run": "ninja: no work to do.\n",
            "flags": {"BUILD_SHARED_LIBS": "ON", "CMAKE_BUILD_TYPE": "Release",
                "CMAKE_GENERATOR": "Ninja", "MARCH_NATIVE": "ON", "WH_LTO": "OFF", "WH_PGO_MODE": "OFF"}})
    build["source_files"].update({path: c.digest(path.encode()) for path in c.INPUTS})
    build["command"] = [compiler, "-std=c++11", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(commit), "-I" + str(c.ROOT / "include")]
    build["command"] += [str(c.ROOT / path) for path in c.SOURCES] + ["-ldl", "-o", worker]
    interpreter_path = "/synthetic/python3.8"
    interpreter = {"artifact": artifact(interpreter_path), "dont_write_bytecode": True,
        "ignore_environment": True, "implementation": "cpython", "invoked_path": interpreter_path,
        "isolated": True, "no_site": True, "no_user_site": True, "path": interpreter_path,
        "version": "3.8.20 (synthetic) [GCC 13.3.0]", "version_info": [3, 8, 20, "final", 0]}
    git_version = "git version 2.43.0\n"
    git = {"commit": commit, "upstream_commit": commit, "tree": "2" * 40,
        "tracked_clean": True, "git_artifact": artifact("/usr/bin/git"), "git_invoked_path": "/usr/bin/git",
        "git_version_text": git_version, "git_version_sha256": c.digest(git_version.encode()),
        "required_tracked_inputs": {path: c.digest(path.encode())[:40] for path in c.INPUTS}}
    entries = [{"exists": False, "path": str(path)} for path in c.h.R0_ROOTS]
    preserved = {"entries": entries, "snapshot_sha256": c.digest(c.canonical(entries))}
    p = {"build": build, "build_after": None, "controller_affinity": {
        "affinity_before": [0, 56, 120], "affinity_after": [0], "controller_cpu": 0,
        "sibling_available_before": True, "singleton_verified": True, "target_available_before": True},
        "interpreter": interpreter, "interpreter_after": copy.deepcopy(interpreter),
        "git_before": git, "git_after": copy.deepcopy(git),
        "preserved_before": preserved, "preserved_after": copy.deepcopy(preserved),
        "health_before": health_receipt(), "health_after": health_receipt(12), "errors": [],
        "process": {"exit_code": 0, "timed_out": False, "elapsed_seconds": 0.05}}
    claim = {key: build[key] for key in ("baseline_path", "baseline_sha256", "candidate_path",
        "candidate_sha256", "worker_sha256", "source_commit")}
    claim.update(campaign=c.CAMPAIGN, schema=c.PREFIX + ".claim.v1", created_unix_ns=1,
                 description_sha256=c.description_hash(), build_receipt_sha256=c.digest(c.canonical(build)))
    return {"CLAIM": claim, "raw.jsonl": encoded(rows), "worker.stderr": b"",
            "provenance.json": p, "WORKER_STARTED": c.marker_value(claim)}


def summarize(values, elapsed=0.125):
    return c.adjudicate(values["raw.jsonl"], values["worker.stderr"], values["provenance.json"],
                        values["WORKER_STARTED"], values["CLAIM"], elapsed)


def publish_bundle(root, values):
    """Refresh all transport hashes, deliberately not semantic cross-links."""
    files = {name: value if type(value) is bytes else c.canonical(value) + b"\n"
             for name, value in values.items() if value is not None}
    os.chmod(str(root), 0o700)
    directory = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
    try:
        for name, data in files.items():
            c.h.write_exclusive(directory, name, data, 0o400)
        c.h.publish_json(directory, "COMPLETE", {name: c.digest(data) for name, data in files.items()})
    finally:
        os.close(directory)


class OffsetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.receipt = fixture()

    def parse(self, rows):
        return c.parse_transcript(encoded(rows), self.receipt)

    def test_roster_and_independent_counts(self):
        roster = list(c.roster())
        self.assertEqual(len(roster), 1152)
        self.assertEqual(len(set(roster)), 1152)
        self.assertEqual(sum(10 * 8192 for _ in roster), 94371840)
        self.assertEqual(sum(8 * (8192 // c.CELLS[cell][0]) for _, cell, _, _ in roster), 3345408)
        self.assertEqual(c.description()["internal_deadline_seconds"], 60)

    def test_mapping_normalization_preserves_candidate_direction(self):
        stats, failed = self.parse(self.rows)
        self.assertEqual(failed, [])
        self.assertEqual(len(stats["controls"]), 48)
        for item in stats["cells"] + stats["widths"]:
            self.assertAlmostEqual(item["summary"]["geometric_mean"], 0.9, places=12)
            self.assertEqual(item["summary"]["n"], 12)

    def test_control_failure_cannot_be_averaged_away(self):
        rows = copy.deepcopy(self.rows)
        for panel in rows[1:-1]:
            if panel["cell_index"] == 0 and panel["comparison_index"] == 0:
                for slot in panel["slots"]:
                    if slot["physical_lane"] == 0:
                        slot["elapsed_ns"] = 1100000
        _, failed = self.parse(rows)
        self.assertEqual(len([x for x in failed if x.startswith("AA:0:0:")]), 4)

    def test_width_and_cell_gates_remain_required(self):
        samples = {key: 0.0 for key in c.roster()}
        _, failed = c.statistics(samples)
        self.assertEqual(failed, ["width-gain:64"])
        for key in samples:
            rep, cell, comparison, condition = key
            if comparison >= 2:
                samples[key] = math.log(1.04 if cell == 0 else .9) * (1 if condition < 2 else -1)
        _, failed = c.statistics(samples)
        self.assertIn("cell-regression:8:64", failed)

    def test_exact_addresses_and_retention(self):
        for field in ("source", "output", "counters", "span", "arena_bytes", "handles"):
            rows = copy.deepcopy(self.rows)
            addresses = rows[1]["addresses"]
            if field == "handles":
                addresses[field][1] = addresses[field][0]
            elif isinstance(addresses[field], str):
                addresses[field] = hex(int(addresses[field], 16) + 64)
            else:
                addresses[field] += 64
            with self.subTest(field=field), self.assertRaises(c.h.ValidationError):
                self.parse(rows)
        rows = copy.deepcopy(self.rows)
        rows[2]["addresses"]["handles"][0] = "0x35000000"
        with self.assertRaises(c.h.ValidationError):
            self.parse(rows)

    def test_binding_and_oracle_mutations_rejected(self):
        mutations = (
            lambda x: x["bindings"]["C"]["memcpy"].update(address="0xa01a14c0"),
            lambda x: x["bindings"]["C"]["symbols"]["wirehair_v2_encode"].update(path=str(c.BASELINE)),
            lambda x: x["bindings"]["C"].update(sha256="f" * 64),
            lambda x: x["cells"][0].update(profile_hex="00" * 32),
            lambda x: x["cells"][0]["partial"][0].update(tail_bytes=2),
            lambda x: x["cells"][0].update(source_sha256="f" * 64),
        )
        for mutation in mutations:
            rows = copy.deepcopy(self.rows)
            mutation(rows[0])
            with self.subTest(mutation=mutation), self.assertRaises(c.h.ValidationError):
                self.parse(rows)

    def test_roster_type_and_counter_tampering(self):
        for field, value in (("sequence", True), ("comparison_index", 9),
                             ("scope_invocations_per_batch", 0), ("output_sha256", "f" * 64)):
            rows = copy.deepcopy(self.rows)
            rows[1][field] = value
            with self.subTest(field=field), self.assertRaises(c.h.ValidationError):
                self.parse(rows)
        for field, value in (("elapsed_ns", 0), ("elapsed_ns", True), ("physical_lane", True), ("logical_lane", "unknown")):
            rows = copy.deepcopy(self.rows)
            rows[1]["slots"][0][field] = value
            with self.subTest(field=field), self.assertRaises(c.h.ValidationError):
                self.parse(rows)
        rows = copy.deepcopy(self.rows)
        rows[-1]["encode_call_count"] -= 1
        with self.assertRaises(c.h.ValidationError):
            self.parse(rows)

    def test_noncanonical_json_and_missing_records(self):
        raw = encoded(self.rows)
        for bad in (raw[:-1], raw.replace(b'"campaign":', b'"campaign": ', 1),
                    encoded(self.rows[:-2] + self.rows[-1:])):
            with self.subTest(length=len(bad)), self.assertRaises(c.h.ValidationError):
                c.parse_transcript(bad, self.receipt)

    def test_baseline_receipt_is_pinned(self):
        with self.assertRaises(c.h.ValidationError):
            c.validate_basis({"library": str(c.BASELINE), "library_sha256": c.BASELINE_HASH})

    def test_parser_and_statistics_never_probe_live_state(self):
        with mock.patch.object(c.h, "run_checked", side_effect=AssertionError("live probe")), \
                mock.patch.object(c, "file_hash", side_effect=AssertionError("live artifact")):
            self.assertEqual(self.parse(self.rows)[1], [])


class BundleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = evidence_fixture()
        cls.anchor = c.digest(c.canonical(cls.evidence["provenance.json"]["build"]["baseline_build"]))
        with mock.patch.object(c, "BASIS_BUILD_HASH", cls.anchor):
            cls.evidence["summary.json"] = summarize(cls.evidence)
        if cls.evidence["summary.json"]["outcome"] != "pass":
            raise AssertionError(cls.evidence["summary.json"])

    def setUp(self):
        anchor = mock.patch.object(c, "BASIS_BUILD_HASH", self.anchor)
        anchor.start()
        self.addCleanup(anchor.stop)

    def values(self):
        return copy.deepcopy(self.evidence)

    def test_complete_bundle_replays_without_live_probes(self):
        with tempfile.TemporaryDirectory(prefix="offset-replay-test-") as temporary:
            root = Path(temporary)
            publish_bundle(root, self.evidence)
            with mock.patch.object(c, "command", side_effect=AssertionError("live command")), \
                    mock.patch.object(c, "file_hash", side_effect=AssertionError("live artifact")), \
                    mock.patch.object(c.h, "run_worker", side_effect=AssertionError("worker launched")), \
                    mock.patch.object(c.h, "cpu_receipt", side_effect=AssertionError("live health")):
                result = c.replay(root)
            self.assertEqual(result, self.evidence["summary.json"])
            self.assertFalse(result["WH1_compared"])
            self.assertFalse(result["promotion_claimed"])

    def test_invalid_launch_bundle_replays_without_marker(self):
        values = self.values()
        values["raw.jsonl"], values["WORKER_STARTED"] = b"", None
        values["provenance.json"]["process"]["exit_code"] = None
        values["provenance.json"]["errors"] = ["execution:synthetic launch failure"]
        values["summary.json"] = summarize(values)
        self.assertEqual(values["summary.json"]["outcome"], "invalid")
        with tempfile.TemporaryDirectory(prefix="offset-invalid-test-") as temporary:
            root = Path(temporary)
            publish_bundle(root, values)
            self.assertEqual(c.replay(root), values["summary.json"])
            self.assertFalse((root / "WORKER_STARTED").exists())

    def test_process_types_and_time_are_exact(self):
        for field, value in (("exit_code", False), ("exit_code", "0"), ("exit_code", 256),
                             ("timed_out", 0), ("elapsed_seconds", True), ("elapsed_seconds", -1),
                             ("elapsed_seconds", 1), ("elapsed_seconds", float("nan"))):
            values = self.values()
            values["provenance.json"]["process"][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(c.h.ValidationError):
                summarize(values)
        for field, value in (("exit_code", None), ("exit_code", -9), ("exit_code", 1), ("timed_out", True)):
            values = self.values()
            values["provenance.json"]["process"][field] = value
            self.assertIn("worker-failure", summarize(values)["infrastructure_failures"])

    def test_missing_postflight_and_bad_preflight_are_rejected(self):
        mutations = (
            lambda p: p.pop("build_after"),
            lambda p: p["controller_affinity"].update(affinity_after=[120]),
            lambda p: p["controller_affinity"].update(affinity_before=[0, 56, 56, 120]),
            lambda p: p["controller_affinity"].update(affinity_before=[-1, 0, 56, 120]),
            lambda p: p["controller_affinity"].update(affinity_before=[0, 56, 120, 1048576]),
            lambda p: p["interpreter"].update(isolated=False),
            lambda p: p["git_before"].update(tracked_clean=False),
            lambda p: p["git_before"]["required_tracked_inputs"].pop(c.INPUTS[0]),
        )
        for index, mutation in enumerate(mutations):
            values = self.values()
            mutation(values["provenance.json"])
            with self.subTest(index=index), self.assertRaises(c.h.ValidationError):
                summarize(values)

    def test_postflight_drift_and_invalid_health_or_preservation_fail(self):
        mutations = (
            lambda p: p.update(build_after={}),
            lambda p: p["interpreter_after"]["artifact"].update(sha256="f" * 64),
            lambda p: p["git_after"].update(tree="f" * 40),
            lambda p: p["health_before"].update(target_cpu=0),
            lambda p: p["health_after"]["edac"][0].update(value=1),
            lambda p: p["health_after"]["thermal"][0].update(value_millic=1000000),
            lambda p: p.update(health_after=health_receipt(16)),
            lambda p: p.update(preserved_before={}, preserved_after={}),
            lambda p: p["preserved_after"].update(snapshot_sha256="f" * 64),
        )
        for index, mutation in enumerate(mutations):
            values = self.values()
            mutation(values["provenance.json"])
            with self.subTest(index=index):
                result = summarize(values)
                self.assertEqual(result["outcome"], "invalid")
                self.assertTrue(any(x.startswith("postflight:") for x in result["infrastructure_failures"]))

    def test_elapsed_deadline_is_not_a_promotion_gate(self):
        for elapsed in (70, 70.5):
            result = summarize(self.evidence, elapsed)
            self.assertEqual(result["outcome"], "invalid")
            self.assertIn("outer-deadline", result["infrastructure_failures"])
        for elapsed in (-1, True, float("nan"), float("inf")):
            with self.subTest(elapsed=elapsed), self.assertRaises(c.h.ValidationError):
                summarize(self.evidence, elapsed)

    def test_marker_links_and_stderr_are_required(self):
        for field in ("campaign", "schema", "claim_sha256", "description_sha256", "worker_sha256"):
            values = self.values()
            values["WORKER_STARTED"][field] = "wrong"
            self.assertEqual(summarize(values)["outcome"], "invalid")
        values = self.values()
        values["worker.stderr"] = b"unexpected worker diagnostic\n"
        self.assertIn("worker-failure", summarize(values)["infrastructure_failures"])

    def test_rehashed_tamper_cannot_retain_a_passing_summary(self):
        def changed_build(values):
            build = values["provenance.json"]["build"]
            build["source_files"]["include/wirehair/wirehair.h"] = "f" * 64
            values["CLAIM"]["build_receipt_sha256"] = c.digest(c.canonical(build))
            values["WORKER_STARTED"] = c.marker_value(values["CLAIM"])

        mutations = (
            lambda v: v["provenance.json"]["process"].update(exit_code=False),
            lambda v: v["provenance.json"].pop("build_after"),
            lambda v: v["WORKER_STARTED"].update(worker_sha256="f" * 64),
            lambda v: v["summary.json"].update(promotion_claimed=True),
            lambda v: v.update(**{"worker.stderr": b"warning\n"}),
            changed_build,
        )
        for index, mutation in enumerate(mutations):
            values = self.values()
            mutation(values)
            with self.subTest(index=index), tempfile.TemporaryDirectory(prefix="offset-tamper-test-") as temporary:
                root = Path(temporary)
                publish_bundle(root, values)
                with self.assertRaises(c.h.ValidationError):
                    c.replay(root)

    def test_recorded_build_commands_and_source_roster_are_checked(self):
        for field in ("command", "source_files"):
            values = self.values()
            build = values["provenance.json"]["build"]
            if field == "command":
                build[field].insert(1, "-ffast-math")
            else:
                build[field].pop(c.INPUTS[0])
            values["CLAIM"]["build_receipt_sha256"] = c.digest(c.canonical(build))
            values["WORKER_STARTED"] = c.marker_value(values["CLAIM"])
            with self.subTest(field=field), self.assertRaises(c.h.ValidationError):
                summarize(values)


if __name__ == "__main__":
    unittest.main()
