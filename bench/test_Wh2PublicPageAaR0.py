#!/usr/bin/env python3
"""Synthetic page-layout protocol tests; never launch the timed worker."""
import contextlib
import copy
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
import Wh2PublicPageAaR0 as c
import test_Wh2PublicAaCalibrationR0 as calibration_fixtures


def transcript():
    expected = c.h.synthetic_expected()
    old = c.h.synthetic_config(expected, "a" * 64)
    config = {key: old[key] for key in ("compile", "runtime_library_maps_sha256",
                                       "runtime_library_path", "target_identity")}
    config.update(campaign=c.CAMPAIGN, schema=c.PREFIX + ".config.v1",
                  description_sha256=c.description_hash(),
                  copy_binding={"address": "0x701a14c0", "elf_offset": "0x1a14c0",
                                "path": c.COPY_LIBRARY},
                  cells=[cell for cell in old["cells"]
                         if (cell["K"], cell["block_bytes"]) in c.CELLS])
    c.validate_config(config, expected, c.description())
    rows = [config]
    for sequence, (rep, unit, scenario, _, condition) in enumerate(c.roster()):
        k, b = c.CELLS[unit // 3]
        role, order = c.ROLES[unit % 3], "ABBA" if condition % 2 == 0 else "BAAB"
        arena = 0x10000000
        span = ((k * b + 8191) // 4096) * 4096
        handles = ["0x30000000", "0x40000000"] if role != "M" else ["0x0", "0x0"]
        if scenario == 3:
            handles[1] = handles[0]
        addresses = {"arena": hex(arena), "source": hex(arena + 0x800), "span": span,
                     "outputs": [hex(arena + span * (i + 1) + c.PHASES[scenario][i])
                                 for i in range(2)],
                     "counters": [hex(arena + span * (i + 3) + 0x400) for i in range(2)],
                     "handles": handles}
        source_hash = c.h.SOURCE_SHA256[k, b, b]
        rows.append({
            "K": k, "block_bytes": b, "condition": condition, "scenario": scenario,
            "description_sha256": c.description_hash(), "replicate": rep, "role": role,
            "mapping": "direct" if condition < 2 else "swapped", "order": order,
            "runtime_library_maps_sha256": "a" * 64,
            "copy_binding_sha256": c.digest(c.canonical(config["copy_binding"])),
            "schema": c.PREFIX + ".panel.v1", "scope": "prebuilt-systematic",
            "sequence": sequence, "scope_invocations_per_batch": 8192 // k,
            "source_immutable_verified": True, "target_cpu": 120,
            "tuple_index": unit // 3, "unit_index": unit, "unit_key_sha256": c.unit_key(rep, unit),
            "left_output_sha256": source_hash, "right_output_sha256": source_hash,
            "addresses": addresses,
            "slots": [{"elapsed_ns": 1000000, "logical_lane": side,
                       "physical_lane": (side == "right") ^ (condition // 2)}
                      for side in c.h.sides_for(order)],
        })
    return rows + [c.terminal_value()], expected


def encoded(rows):
    return b"".join(c.canonical(row) + b"\n" for row in rows)


def bundle_fixture():
    fixtures = calibration_fixtures.fixtures
    p = copy.deepcopy(calibration_fixtures.bundle_fixture()["provenance.json"])
    p.update(description=c.description(), worker_path=p["build_root"] + "/" + c.TARGET,
             controller_path=p["source_root"] + "/bench/Wh2PublicPageAaR0.py")
    for key in ("worker", "controller"):
        p["artifacts_before"][key] = fixtures.artifact(p[key + "_path"])
    p["git_before"]["required_tracked_inputs"] = {
        path: c.digest(path.encode("ascii"))[:40] for path in c.INPUTS}
    commands = c.a.expected_build_tokens(p, p["compiler_path"])
    commands.append(["/usr/bin/cmake", "-E", "cmake_symlink_library", "libwirehair.so.2.0.0",
                     "libwirehair.so.2", "libwirehair.so", "&&", ":"])
    commands = "\n".join(shlex.join(row) for row in commands) + "\n"
    p["build_before"].update(commands=commands, commands_sha256=c.digest(commands.encode("ascii")))
    p["dynamic_before"] = fixtures.dynamic_fixture(
        p["worker_path"], p["library_path"], p["artifacts_before"]["library"])
    entries = [{"exists": False, "path": str(path)} for path in c.h.R0_ROOTS]
    p["preserved_before"] = {"entries": entries, "snapshot_sha256": c.digest(c.canonical(entries))}
    p["copy_library_before"] = fixtures.artifact(c.COPY_LIBRARY, c.COPY_LIBRARY_SHA256)
    p["postflight_errors"] = dict.fromkeys(c.a.POSTFLIGHT)
    for key in c.a.POSTFLIGHT:
        p[key + "_after"] = copy.deepcopy(p[key + "_before"] if key != "interpreter" else p["interpreter"])
    claim = {"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".claim.v1", "created_unix_ns": 1,
             "description_sha256": c.description_hash(), "source_commit": p["source_commit"],
             "controller_interpreter_sha256": p["interpreter"]["artifact"]["sha256"]}
    claim.update({key + "_sha256": p["artifacts_before"][key]["sha256"]
                  for key in ("worker", "library", "controller")})
    marker = c.a.marker_value(claim, c.digest(c.canonical(claim) + b"\n"))
    p["marker_sha256"] = c.digest(c.canonical(marker) + b"\n")
    rows, _ = transcript()
    raw = encoded(rows)
    summary = c.a.adjudicate_bundle(raw, b"", p, marker, claim, 0.25)
    c.validate_provenance(p, claim)
    return {"CLAIM": claim, "WORKER_STARTED": marker, "raw.jsonl": raw, "worker.stderr": b"",
            "provenance.json": p, "summary.json": summary}


def publish(root, fixture):
    data = {name: value if type(value) is bytes else c.canonical(value) + b"\n"
            for name, value in fixture.items()}
    data["COMPLETE"] = c.canonical({
        "campaign": c.CAMPAIGN, "schema": c.PREFIX + ".complete.v1", "status": "complete",
        "outcome": fixture["summary.json"]["outcome"],
        "files": {name: c.digest(value) for name, value in data.items()}}) + b"\n"
    fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
    try:
        for name, value in data.items():
            c.h.write_exclusive(fd, name, value, 0o400)
    finally:
        os.close(fd)


class PageAaTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.expected = transcript()

    def parse(self, rows):
        return c.parse_transcript(encoded(rows), self.expected, c.description())

    def test_exact_roster_and_calls(self):
        roster = list(c.roster())
        self.assertEqual(len(roster), 1728)
        self.assertEqual(len(set(roster)), 1728)
        self.assertEqual(sum(10 * 8192 for _ in roster), 141557760)
        self.assertEqual(sum(8 * (8192 // c.CELLS[unit // 3][0])
                             for _, unit, _, _, _ in roster), 9732096)
        for rep, unit, scenario, pos, cond in roster:
            self.assertEqual(cond, c.CONDITION_SEQUENCES[rep % 4][pos])
            self.assertIn(scenario, range(4))
        self.assertEqual(c.a.description(), c.description())

    def test_clean_transcript_and_gate_roster(self):
        stats, failed = self.parse(self.rows)
        self.assertEqual(failed, [])
        self.assertEqual(len(stats["conditions"]), 144)
        self.assertEqual(sum(row["primary"] for row in stats["conditions"]), 24)
        self.assertEqual(len(stats["matched_replicates"]), 36)
        self.assertEqual(sum(row["primary"] for row in stats["matched_replicates"]), 6)

    def test_condition_failures_cannot_average_away(self):
        rows = copy.deepcopy(self.rows)
        for panel in rows[1:-1]:
            if panel["scenario"] != 0 or panel["role"] == "M":
                continue
            for slot in panel["slots"]:
                slot["elapsed_ns"] = 1100000 if slot["physical_lane"] == 0 else 1000000
        stats, failed = self.parse(rows)
        self.assertEqual(len(failed), 24)
        self.assertFalse(any("matched" in item for item in failed))
        for row in stats["matched_replicates"]:
            self.assertAlmostEqual(row["summary"]["log_mean"], 0.0, places=12)

    def test_diagnostics_cannot_decide_primary(self):
        rows = copy.deepcopy(self.rows)
        for panel in rows[1:-1]:
            if panel["scenario"] != 0 or panel["role"] == "M":
                for slot in panel["slots"]:
                    slot["elapsed_ns"] *= 2 if slot["logical_lane"] == "left" else 1
        _, failed = self.parse(rows)
        self.assertEqual(failed, [])

    def test_phase_handle_span_and_canary_receipts_reject(self):
        changes = [
            lambda p: p["addresses"].update(source="0x10000820"),
            lambda p: p["addresses"]["outputs"].__setitem__(0, "0x10000000"),
            lambda p: p["addresses"].update(span=True),
            lambda p: p["addresses"].update(counters=list(reversed(p["addresses"]["counters"]))),
            lambda p: p["addresses"].update(handles=["0x0", "0x0"]),
            lambda p: p.update(source_immutable_verified=False),
            lambda p: p.update(left_output_sha256="0" * 64),
            lambda p: p.update(copy_binding_sha256="0" * 64),
        ]
        target = next(i for i, p in enumerate(self.rows) if p.get("role") == "C" and p.get("scenario") == 0)
        for change in changes:
            with self.subTest(change=change):
                rows = copy.deepcopy(self.rows)
                change(rows[target])
                with self.assertRaises(c.h.ValidationError):
                    self.parse(rows)

    def test_binding_identity_and_roster_reject(self):
        for field, value in (("path", "/tmp/libc.so"), ("elf_offset", "0x1a14c1"),
                             ("address", "0x0")):
            rows = copy.deepcopy(self.rows)
            rows[0]["copy_binding"][field] = value
            with self.assertRaises(c.h.ValidationError):
                self.parse(rows)
        for rows in (self.rows[:-1], self.rows + [self.rows[-1]],
                     [self.rows[0], self.rows[2], self.rows[1]] + self.rows[3:]):
            with self.assertRaises(c.h.ValidationError):
                self.parse(rows)

    def test_noncanonical_and_nonpositive_time_reject(self):
        raw = encoded(self.rows)
        with self.assertRaises(c.h.ValidationError):
            c.parse_transcript(raw.replace(b'"K":8', b'"K":8.0', 1),
                               self.expected, c.description())
        for bad in (0, -1, True, 1.5):
            rows = copy.deepcopy(self.rows)
            rows[1]["slots"][0]["elapsed_ns"] = bad
            with self.assertRaises(c.h.ValidationError):
                self.parse(rows)

    def test_balanced_warmup_predecessors(self):
        for condition in c.description()["conditions"]:
            prior = condition["warmup"][-1]
            same, other = {"left": 0, "right": 0}, {"left": 0, "right": 0}
            for side in c.h.sides_for(condition["order"]):
                (same if prior == side else other)[side] += 1
                prior = side
            self.assertEqual(same, {"left": 1, "right": 1})
            self.assertEqual(other, {"left": 3, "right": 3})

    def test_full_bundle_replay_and_copy_postflight_failure(self):
        fixture = bundle_fixture()
        for changed in (False, True):
            current = copy.deepcopy(fixture)
            if changed:
                p = current["provenance.json"]
                p["copy_library_after"]["sha256"] = "0" * 64
                current["summary.json"] = c.a.adjudicate_bundle(
                    current["raw.jsonl"], b"", p, current["WORKER_STARTED"], current["CLAIM"], 0.25)
                self.assertEqual(current["summary.json"]["outcome"], "invalid")
                self.assertTrue(any(item.startswith("postflight_copy_library:")
                                    for item in current["summary.json"]["infrastructure_failures"]))
            with tempfile.TemporaryDirectory(prefix="wh2-page-replay-") as root:
                publish(Path(root), current)
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(c.a.replay(Path(root)), 1 if changed else 0)

    def test_rehashed_provenance_drift_rejects(self):
        fixture = bundle_fixture()
        mutations = (
            lambda p: p["copy_library_before"].update(sha256="0" * 64),
            lambda p: p["artifacts_before"]["worker"].update(sha256="0" * 64),
            lambda p: p["git_before"]["required_tracked_inputs"].pop("bench/Wh2PublicPageAaR0.cpp"),
            lambda p: p.update(controller_path=p["source_root"] + "/bench/Wh2PublicAaCalibrationR0.py"),
            lambda p: p["build_before"]["cache"].update(WH_LTO="ON"),
        )
        for mutation in mutations:
            current = copy.deepcopy(fixture)
            mutation(current["provenance.json"])
            with tempfile.TemporaryDirectory(prefix="wh2-page-tamper-") as root:
                publish(Path(root), current)  # Recompute all COMPLETE hashes deliberately.
                with self.assertRaises(c.h.ValidationError), contextlib.redirect_stdout(io.StringIO()):
                    c.a.replay(Path(root))

    def test_private_callbacks_and_mock_preflight(self):
        p = bundle_fixture()["provenance.json"]
        args = SimpleNamespace(source_root=p["source_root"], build_root=p["build_root"],
                               worker=p["worker_path"], library=p["library_path"], compiler=p["compiler_path"],
                               expected_source_commit=p["source_commit"], **{
                                   "expected_" + key + "_sha256": p["artifacts_before"][key]["sha256"]
                                   for key in ("worker", "library", "controller")})
        self.assertIs(c.a.preflight, c.preflight)
        self.assertIs(c.a.parse_transcript, c.parse_transcript)
        self.assertIs(c.a.validate_provenance, c.validate_provenance)
        self.assertEqual(c.h.FIXED_OUTPUT_DIR, c.FIXED_OUTPUT_DIR)
        self.assertNotEqual(calibration_fixtures.c.CAMPAIGN, c.a.CAMPAIGN)
        with contextlib.ExitStack() as stack:
            for owner, name, value in (
                (c, "__file__", p["controller_path"]),
                (c.h, "exact_absolute_directory", lambda value, where: Path(value)),
                (c.h, "exact_absolute_file", lambda value, where: Path(value)),
                (c.h, "compiler_receipt", lambda path: copy.deepcopy(p["compiler"])),
                (c.h, "interpreter_receipt", lambda: copy.deepcopy(p["interpreter"])),
                (c, "artifacts", lambda value: copy.deepcopy(p["artifacts_before"])),
                (c.h, "git_receipt", lambda *args: copy.deepcopy(p["git_before"])),
                (c, "build_receipt", lambda value: copy.deepcopy(p["build_before"])),
                (c.h, "dynamic_receipt", lambda *args: copy.deepcopy(p["dynamic_before"])),
                (c.h, "snapshot_r0", lambda: copy.deepcopy(p["preserved_before"])),
                (c.h, "cpu_receipt", lambda: copy.deepcopy(p["health_before"])),
                (c.h, "artifact_receipt", lambda path: copy.deepcopy(p["copy_library_before"])),
            ):
                stack.enter_context(mock.patch.object(owner, name, value))
            probe = stack.enter_context(mock.patch.object(
                c.h, "run_checked", return_value=c.canonical(c.description()) + b"\n"))
            result = c.preflight(args)
            self.assertEqual(result["copy_library_before"], p["copy_library_before"])
            self.assertEqual(result["description"], c.description())
            self.assertEqual(probe.call_args[0][0], [p["worker_path"], "--describe"])
            self.assertEqual(probe.call_args[1]["env"], c.environment(p["library_path"]))
            with mock.patch.object(c, "__file__", "/tmp/controller-copy.py"), \
                    mock.patch.object(c.h, "compiler_receipt") as compiler:
                with self.assertRaisesRegex(c.h.ValidationError, "controller entrypoint"):
                    c.preflight(args)
                compiler.assert_not_called()

    def test_mock_producer_complete_last_replay_and_exclusive_claim(self):
        for launch_failure in (False, True):
            fixture = bundle_fixture()
            p = fixture["provenance.json"]
            published = []
            original_publish = c.h.publish_json

            def publish_json(fd, name, value):
                published.append(name)
                original_publish(fd, name, value)

            with tempfile.TemporaryDirectory(prefix="wh2-page-producer-") as parent:
                root = Path(parent) / "bundle"

                def worker(*args):
                    if launch_failure:
                        raise c.h.ValidationError("synthetic launch failure")
                    claim = c.h.parse_canonical_document(root / "CLAIM", c.h.MAX_PROBE_BYTES)
                    fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
                    try:
                        c.h.publish_json(fd, "WORKER_STARTED", c.a.marker_value(claim, args[3]))
                    finally:
                        os.close(fd)
                    return fixture["raw.jsonl"], b"", 0, 0.1, False

                probes = {key: (lambda value, key=key: copy.deepcopy(
                    value[key + "_before"] if key != "interpreter" else value["interpreter"]))
                          for key in c.a.POSTFLIGHT}
                with contextlib.ExitStack() as stack:
                    for owner, name, value in (
                        (c.a, "FIXED_OUTPUT_DIR", root), (c.h, "FIXED_OUTPUT_DIR", root),
                        (c.a, "preflight", lambda args: p),
                        (c.h, "pin_controller", lambda: p["controller_affinity"]),
                        (c.h, "run_worker", worker), (c.a, "POSTFLIGHT", probes),
                        (c.h, "publish_json", publish_json),
                    ):
                        stack.enter_context(mock.patch.object(owner, name, value))
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                    result = c.a.run_once(SimpleNamespace())
                    self.assertEqual(result, 1 if launch_failure else 0)
                    self.assertEqual(published[-1], "COMPLETE")
                    self.assertEqual(c.a.replay(root), result)
                    with self.assertRaisesRegex(c.h.ValidationError, "namespace already exists"):
                        c.a.run_once(SimpleNamespace())


if __name__ == "__main__":
    unittest.main()
