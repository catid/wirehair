#!/usr/bin/env python3
"""Synthetic protocol, parser, replay and one-shot publication tests; no timing."""
import contextlib
import copy
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import shlex
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2PublicAaCalibrationR0 as c
import test_Wh2PublicBorrowedCurrentVsWh1R1 as fixtures


def transcript(ratio=1.0):
    expected = c.h.synthetic_expected()
    original = c.h.synthetic_config(expected, "a" * 64)
    config = {key: original[key] for key in ("compile", "runtime_library_maps_sha256",
                                              "runtime_library_path", "target_identity")}
    config.update(campaign=c.CAMPAIGN, schema=c.PREFIX + ".config.v1",
                  description_sha256=c.description_hash(),
                  cells=[cell for cell in original["cells"] if (cell["K"], cell["block_bytes"]) in c.CELLS])
    oracles = c.validate_config(config, expected, c.description())
    records = [config]
    for sequence, (rep, unit, _, condition) in enumerate(c.roster()):
        k, b, scope = c.TUPLES[unit // 2]
        role = ("C", "L")[unit % 2]
        order = "ABBA" if condition % 2 == 0 else "BAAB"
        output = oracles[k, b, role]["systematic_sha256" if scope.endswith("systematic") else "full_repair_sha256"]
        records.append({"K": k, "block_bytes": b, "condition": condition,
                        "description_sha256": c.description_hash(), "left_output_sha256": output,
                        "right_output_sha256": output, "mapping": "direct" if condition < 2 else "swapped",
                        "order": order, "physical_scratch_addresses": ["0x10000000", "0x20000000"],
                        "replicate": rep, "role": role, "runtime_library_maps_sha256": "a" * 64,
                        "schema": c.PREFIX + ".panel.v1", "scope": scope,
                        "scope_invocations_per_batch": c.batch_invocations(k, scope), "sequence": sequence,
                        "slots": [{"elapsed_ns": round(1000000 * (ratio if side == "left" else 1)),
                                   "logical_lane": side, "physical_lane": (0 if side == "left" else 1) ^ (condition // 2)}
                                  for side in c.h.sides_for(order)], "source_immutable_verified": True,
                        "target_cpu": 120, "tuple_index": unit // 2, "unit_index": unit,
                        "unit_key_sha256": c.digest("aa-cal-r0:rep={}:unit={}".format(rep, unit).encode("ascii"))})
    terminal = {"campaign": c.CAMPAIGN, "encode_call_count": 47431680,
                "measured_batch_count": 6144, "measured_invocation_count": 2411520,
                "panel_count": 768, "record_count": 770, "schema": c.PREFIX + ".terminal.v1",
                "status": "complete", "warmup_batch_count": 1536, "warmup_invocation_count": 602880}
    return records + [terminal], expected


def encoded(records):
    return b"".join(c.canonical(record) + b"\n" for record in records)


def bundle_fixture(ratio=1.0):
    old = fixtures.bundle_fixture()["provenance.json"]
    p = {key: old[key] for key in ("compiler", "interpreter", "source_root", "build_root", "library_path", "controller_affinity")}
    p.update(source_commit=old["expected"]["source_commit"], description=c.description(),
             worker_path=p["build_root"] + "/" + c.TARGET,
             controller_path=p["source_root"] + "/bench/Wh2PublicAaCalibrationR0.py",
             compiler_path=p["compiler"]["path"], expected={key: old["expected"][key]
                 for key in ("compiler_path", "compiler_sha256", "compiler_version", "source_commit", "library_path")})
    p["artifacts_before"] = {key: fixtures.artifact(p[key + "_path"]) for key in ("worker", "library", "controller", "compiler")}
    p["artifacts_before"]["compiler"]["sha256"] = p["compiler"]["sha256"]
    p["git_before"] = old["git"]
    p["git_before"]["required_tracked_inputs"] = {path: c.digest(path.encode("ascii"))[:40] for path in c.INPUTS}
    rows = c.expected_build_tokens(p, p["compiler_path"])
    rows.append(["/usr/bin/cmake", "-E", "cmake_symlink_library", "libwirehair.so.2.0.0", "libwirehair.so.2", "libwirehair.so", "&&", ":"])
    commands = "\n".join(shlex.join(row) for row in rows) + "\n"
    p["build_before"] = {"cache": c.BUILD_CACHE, "compiler_invocation": p["compiler_path"],
                         "commands": commands, "commands_sha256": c.digest(commands.encode("ascii")),
                         "dry_run": "ninja: no work to do.\n"}
    p["dynamic_before"] = fixtures.dynamic_fixture(p["worker_path"], p["library_path"], p["artifacts_before"]["library"])
    entries = [{"exists": False, "path": str(path)} for path in c.h.R0_ROOTS]
    p["preserved_before"] = {"entries": entries, "snapshot_sha256": c.digest(c.canonical(entries))}
    p["health_before"] = fixtures.health_receipt()
    p["postflight_errors"] = dict.fromkeys(c.POSTFLIGHT)
    for key in c.POSTFLIGHT:
        p[key + "_after"] = copy.deepcopy(p[key + "_before"] if key != "interpreter" else p["interpreter"])
    p["process"] = {"exit_code": 0, "elapsed_seconds": 0.1, "timed_out": False, "error": None}
    claim = {"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".claim.v1", "created_unix_ns": 1,
             "description_sha256": c.description_hash(), "source_commit": p["source_commit"],
             "controller_interpreter_sha256": p["interpreter"]["artifact"]["sha256"]}
    claim.update({key + "_sha256": p["artifacts_before"][key]["sha256"] for key in ("worker", "library", "controller")})
    marker = c.marker_value(claim, c.digest(c.canonical(claim) + b"\n"))
    p["marker_sha256"] = c.digest(c.canonical(marker) + b"\n")
    records, _ = transcript(ratio)
    raw = encoded(records)
    summary = c.adjudicate_bundle(raw, b"", p, marker, claim, 0.25)
    return {"CLAIM": claim, "WORKER_STARTED": marker, "raw.jsonl": raw, "worker.stderr": b"",
            "provenance.json": p, "summary.json": summary}


def publish(root, fixture):
    data = {key: value if type(value) is bytes else c.canonical(value) + b"\n"
            for key, value in fixture.items() if value is not None}
    complete = {"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".complete.v1", "status": "complete",
                "outcome": fixture["summary.json"]["outcome"], "files": {name: c.digest(value) for name, value in data.items()}}
    data["COMPLETE"] = c.canonical(complete) + b"\n"
    fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
    try:
        for name, value in data.items():
            c.h.write_exclusive(fd, name, value, 0o400)
    finally:
        os.close(fd)


class CalibrationTest(unittest.TestCase):
    def test_exact_roster_counts(self):
        rows = list(c.roster())
        self.assertEqual(len(rows), 768)
        self.assertEqual(len(set(rows)), 768)
        batches = calls = 0
        for rep, unit, pos, condition in rows:
            self.assertEqual(condition, c.CONDITION_SEQUENCES[rep % 4][pos])
            k, _, scope = c.TUPLES[unit // 2]
            batches += 8 * c.batch_invocations(k, scope)
            calls += 8 * c.batch_invocations(k, scope) * k
        self.assertEqual((batches, calls), (2411520, 37945344))
        self.assertEqual(c.unit_order(0), sorted(range(16), key=lambda unit:
            __import__("hashlib").sha256(("aa-cal-r0:rep=0:unit=" + str(unit)).encode()).hexdigest()))

    def test_warmup_boundaries_and_physical_mapping(self):
        for condition in c.description()["conditions"]:
            previous = condition["warmup"][-1]
            counts = {"left": [0, 0], "right": [0, 0]}
            for lane in c.h.sides_for(condition["order"]):
                counts[lane][0 if lane == previous else 1] += 1
                previous = lane
            self.assertEqual(counts, {"left": [1, 3], "right": [1, 3]})

    def test_pass_fail_no_swapped_negation(self):
        for ratio, fails in ((1.0, 0), (1.03, 80)):
            records, expected = transcript(ratio)
            stats, gates = c.parse_transcript(encoded(records), expected, c.description())
            self.assertEqual(len(gates), fails)
            for result in stats["conditions"] + stats["matched_replicates"]:
                self.assertAlmostEqual(result["summary"]["geometric_mean"], ratio)
        # A changed mapping does not synthesize an opposite-signed observation.
        self.assertAlmostEqual(stats["conditions"][0]["summary"]["log_mean"],
                               stats["conditions"][2]["summary"]["log_mean"])

    def test_matched_means_use_replicates_not_48_independent_values(self):
        samples = {(rep, unit, cond): (rep - 5.5) * .001 + (cond - 1.5) * .0003
                   for rep in range(12) for unit in range(16) for cond in range(4)}
        stats, _ = c.statistics(samples)
        expected = c.h.sample_summary([(rep - 5.5) * .001 for rep in range(12)])
        self.assertAlmostEqual(stats["matched_replicates"][0]["summary"]["log_standard_error"], expected["log_standard_error"])
        self.assertEqual(stats["matched_replicates"][0]["summary"]["n"], 12)

    def test_strict_noise_boundary(self):
        boundary = math.log1p(.02)
        samples = {(rep, unit, cond): boundary for rep in range(12) for unit in range(16) for cond in range(4)}
        _, gates = c.statistics(samples)
        self.assertEqual(len(gates), 80)

    def test_mutated_roster_mapping_identity_and_oracles_invalid(self):
        records, expected = transcript()
        mutations = [lambda r: r.pop(), lambda r: r[1].__setitem__("role", "C/L"),
                     lambda r: r[1].__setitem__("scope_invocations_per_batch", 0),
                     lambda r: r[1].__setitem__("condition", True),
                     lambda r: r[1]["slots"][0].__setitem__("physical_lane", 7),
                     lambda r: r[1]["slots"][0].__setitem__("elapsed_ns", 0),
                     lambda r: r[1].__setitem__("runtime_library_maps_sha256", "b" * 64),
                     lambda r: r[2].__setitem__("physical_scratch_addresses", ["0x30000000", "0x40000000"]),
                     lambda r: r[0]["cells"][0]["oracles"][0].__setitem__("roundtrip_repair_only_verified", False)]
        for mutation in mutations:
            altered = copy.deepcopy(records)
            mutation(altered)
            with self.assertRaises(c.h.ValidationError):
                c.parse_transcript(encoded(altered), expected, c.description())

    def test_parser_limits(self):
        records, expected = transcript()
        for raw in (b"", encoded(records)[:-1], b'{"x":NaN}\n', b'{"x":1,"x":2}\n'):
            with self.assertRaises(c.h.ValidationError):
                c.parse_transcript(raw, expected, c.description())
        with mock.patch.object(c.h, "MAX_RAW_BYTES", 100):
            with self.assertRaises(c.h.ValidationError):
                c.parse_transcript(encoded(records), expected, c.description())

    def test_bundle_replay_pass_fail(self):
        for ratio, code in ((1, 0), (1.03, 2)):
            fixture = bundle_fixture(ratio)
            with tempfile.TemporaryDirectory(prefix="wh2-aa-test-") as root:
                publish(Path(root), fixture)
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(c.replay(Path(root)), code)

    def test_replay_rejects_rehashed_semantic_drift(self):
        mutations = [lambda f: f["summary.json"].__setitem__("outcome", "pass"),
                     lambda f: f["WORKER_STARTED"].__setitem__("description_sha256", "f" * 64),
                     lambda f: f["provenance.json"]["artifacts_before"]["worker"].__setitem__("sha256", "f" * 64),
                     lambda f: f["provenance.json"]["build_before"]["cache"].__setitem__("WH_LTO", "ON")]
        for mutation in mutations:
            fixture = bundle_fixture(1.03)
            mutation(fixture)
            with tempfile.TemporaryDirectory(prefix="wh2-aa-test-") as root:
                publish(Path(root), fixture)
                with self.assertRaises(c.h.ValidationError), contextlib.redirect_stdout(io.StringIO()):
                    c.replay(Path(root))

    def test_infrastructure_invalid_not_aa_fail(self):
        fixture = bundle_fixture()
        p, claim, marker = fixture["provenance.json"], fixture["CLAIM"], fixture["WORKER_STARTED"]
        p["health_after"] = fixtures.health_receipt(16)
        result = c.adjudicate_bundle(fixture["raw.jsonl"], b"", p, marker, claim, 0.25)
        self.assertEqual(result["outcome"], "invalid")
        self.assertEqual(result["failed_aa_gates"], [])
        self.assertTrue(any("SMT sibling" in failure for failure in result["infrastructure_failures"]))

    def run_mock_once(self, mode):
        fixture = bundle_fixture()
        p = fixture["provenance.json"]
        published = []
        original_publish = c.h.publish_json

        def publish_json(fd, name, value):
            published.append(name)
            original_publish(fd, name, value)

        with tempfile.TemporaryDirectory(prefix="wh2-aa-producer-") as parent:
            root = Path(parent) / "bundle"

            def worker(*args):
                if mode == "launch_failure":
                    raise c.h.ValidationError("synthetic launch failure")
                claim = c.h.parse_canonical_document(root / "CLAIM", c.h.MAX_PROBE_BYTES)
                fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
                try:
                    if mode == "bad_marker":
                        c.h.write_exclusive(fd, "WORKER_STARTED", b"malformed\n", 0o400)
                    else:
                        c.h.publish_json(fd, "WORKER_STARTED", c.marker_value(claim, args[3]))
                finally:
                    os.close(fd)
                return (b"malformed\n" if mode == "bad_raw" else fixture["raw.jsonl"]), b"", 0, 0.1, False

            probes = {key: (lambda p, key=key: copy.deepcopy(
                p[key + "_before"] if key != "interpreter" else p["interpreter"])) for key in c.POSTFLIGHT}
            with contextlib.ExitStack() as stack:
                for owner, key, value in ((c, "FIXED_OUTPUT_DIR", root), (c.h, "FIXED_OUTPUT_DIR", root),
                                          (c, "preflight", lambda args: p),
                                          (c.h, "pin_controller", lambda: p["controller_affinity"]),
                                          (c.h, "run_worker", worker), (c, "POSTFLIGHT", probes),
                                          (c.h, "publish_json", publish_json)):
                    stack.enter_context(mock.patch.object(owner, key, value))
                stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                result = c.run_once(SimpleNamespace())
                self.assertEqual(published[-1], "COMPLETE")
                self.assertEqual(c.replay(root), result)
                with self.assertRaises(c.h.ValidationError):
                    c.run_once(SimpleNamespace())
            self.assertEqual(result, 0 if mode == "success" else 1)

    def test_producer_complete_last_replay_and_exclusive_claim(self):
        for mode in ("success", "launch_failure", "bad_raw", "bad_marker"):
            with self.subTest(mode=mode):
                self.run_mock_once(mode)


if __name__ == "__main__":
    unittest.main()
