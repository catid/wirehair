#!/usr/bin/env python3
"""Synthetic crossover tests; never launch a codec or timing worker."""
import copy
import hashlib
import itertools
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2PublicStackCrossoverR0 as c
from test_Wh2PublicOffsetR0 import evidence_fixture as offset_evidence, publish_bundle

CELLS = ((8, 64), (8, 1280), (128, 64), (128, 1280), (8192, 64), (8192, 1280))
CONDITIONS = ((0, 1, 3, 2), (1, 2, 0, 3), (2, 3, 1, 0), (3, 0, 2, 1))


def expected_roster():
    permutations = tuple(itertools.permutations(range(3)))
    for rep in range(12):
        cells = sorted(range(6), key=lambda cell: hashlib.sha256(
            ("stack-crossover-r0:rep=%d:cell=%d" % (rep, cell)).encode("ascii")).digest())
        for cell in cells:
            for target in permutations[(rep + cell) % 6]:
                for index in range(3):
                    comparison = (rep + cell + index) % 3
                    for condition in CONDITIONS[rep % 4]:
                        yield rep, cell, target, comparison, condition


def samples(delta0=0.04, delta1=-0.04, common_drift=True):
    result = {}
    for key in expected_roster():
        rep, cell, target, comparison, condition = key
        value = 0.0
        if comparison == 2:
            # A fixed handle bias need not cross one. Shared replicate drift
            # cancels only when target-minus-null contrasts are paired first.
            value = 0.2 + (0.03 * (rep - 5) if common_drift else 0)
            value += (delta0, delta1, 0.0)[target]
            value *= 1 if condition < 2 else -1
        result[key] = value
    return result


def addresses(phases=(0, 0x10, 0x200, 0x400)):
    arena, span = 0x10000000, 8192
    return {"arena": hex(arena), "arena_bytes": 2 * span + 1024 + 32768 + 64,
            "span": span, "source": hex(arena + 2048), "output": hex(arena + span + 2336),
            "counters": hex(arena + 2 * span + 1024),
            "handles": [hex(0x30000000 + index * 0x10000000 + phase)
                        for index, phase in enumerate(phases)]}


def probe(phase=0x800):
    pre = 0x7fff00000000 + phase
    return {"pre_rsp": hex(pre), "hot_rsp": hex(pre), "final_rsp": hex(pre),
            "frame_address": hex(pre + 96), "correction": 0}


def expected_geometry(addrs, stack):
    def window(target):
        return {(target - delta) % 4096 for delta in range(1, 41)}

    pre = int(stack["pre_rsp"], 16)
    fixed = {(pre + offset) % 4096 for offset in range(256)}
    targets = [(int(handle, 16) + 0x110) % 4096 for handle in addrs["handles"]]
    for i, j in itertools.combinations(range(4), 2):
        left, right = window(targets[i]), window(targets[j])
        if left & right or (left | right) & fixed:
            continue
        for null in range(0, 4096, 16):
            if not window(null) & (left | right | fixed):
                return {"selected_pair": [i, j], "target_phases": [targets[i], targets[j], null],
                        "probe": stack, "fixed_envelope_bytes": 256}
    raise ValueError("no synthetic geometry")


def encoded(rows):
    return b"".join(c.canonical(row) + b"\n" for row in rows)


def transcript_fixture():
    receipt = {"baseline_path": str(c.BASELINE), "baseline_sha256": c.BASELINE_HASH}
    symbols = {name: {"path": str(c.BASELINE), "address": hex(0x70000000 + 0x485c0 + index * 128),
                      "elf_offset": hex(0x485c0 + index * 128)}
               for index, name in enumerate(("wirehair_init_", "wirehair_v2_encoder_create_with_options",
                                              "wirehair_v2_encode", "wirehair_v2_free"))}
    config = {"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".config.v1",
              "description_sha256": c.description_hash(),
              "target_identity": c.h.synthetic_config(c.h.synthetic_expected(), "a" * 64)["target_identity"],
              "bindings": {"baseline": {"path": str(c.BASELINE), "sha256": c.BASELINE_HASH, "symbols": symbols,
                  "memcpy": {"path": str(c.LIBC), "address": "0x901a14c0", "elf_offset": "0x1a14c0"}}}, "cells": []}
    for K, B in CELLS:
        def oracle(tail):
            raw = b"WHV2" + struct.pack("<HHQQIB3s", 1, 32, 0x4b295bbb47f4f9c9, (K - 1) * B + tail, B, 0, bytes(3))
            return {"profile_hex": raw.hex(), "profile_sha256": c.digest(raw),
                    "source_sha256": c.h.SOURCE_SHA256[K, B, tail], "repair_sha256": "d" * 64,
                    "high_id_sha256": "e" * 64}
        partial = [dict(oracle(tail), tail_bytes=tail, systematic_sha256=c.h.SOURCE_SHA256[K, B, tail])
                   for tail in (1, B - 1)]
        config["cells"].append(dict(oracle(B), K=K, block_bytes=B, partial=partial))
    rows = [config]
    for sequence, (rep, cell, target, comparison, condition) in enumerate(expected_roster()):
        K, B = CELLS[cell]
        addrs, initial = addresses(), probe()
        arena, span = int(addrs["arena"], 16), ((K * B + 8191) // 4096) * 4096
        addrs.update(span=span, arena_bytes=2 * span + 1024 + 32768 + 64,
                     output=hex(arena + span + 2336), counters=hex(arena + 2 * span + 1024))
        geometry = expected_geometry(addrs, initial)
        phase, pair = geometry["target_phases"][target], geometry["selected_pair"]
        pre = int(initial["pre_rsp"], 16)
        correction = 4096 + ((pre - phase) % 4096)
        stack = dict(initial, hot_rsp=hex(pre - correction), final_rsp=hex(pre - correction), correction=correction)
        selected = (pair[0], pair[0]) if comparison == 0 else (pair[1], pair[1]) if comparison == 1 else pair

        def slots(bits, warm):
            result = []
            for bit in bits:
                logical = bit ^ (condition % 2)
                physical = logical ^ (condition // 2)
                handle = selected[physical]
                elapsed = 1100000 + (80000, -80000, 0)[target] if handle == pair[0] else 1000000
                result.append({"logical_lane": "right" if logical else "left", "physical_lane": physical,
                               "handle_index": handle, "stack": dict(stack), "elapsed_ns": 0 if warm else elapsed})
            return result

        rows.append({"campaign": c.CAMPAIGN, "schema": c.PREFIX + ".panel.v1",
            "description_sha256": c.description_hash(), "sequence": sequence, "replicate": rep,
            "cell_index": cell, "cell_key_sha256": c.cell_key(rep, cell), "target_index": target,
            "comparison_index": comparison, "condition": condition, "order": "ABBA" if condition % 2 == 0 else "BAAB",
            "mapping": "direct" if condition < 2 else "swapped", "K": K, "block_bytes": B,
            "scope_invocations_per_batch": 8192 // K, "source_sha256": c.h.SOURCE_SHA256[K, B, B],
            "output_sha256": c.h.SOURCE_SHA256[K, B, B], "profile_sha256": config["cells"][cell]["profile_sha256"],
            "source_immutable_verified": True, "addresses": addrs, "geometry": geometry, "target_phase": phase,
            "warmup": slots((0, 1), True), "slots": slots((0, 1, 1, 0, 1, 0, 0, 1), False)})
    return rows + [c.terminal_value()], receipt


def evidence_fixture(raw):
    """Reuse only pure synthetic Offset receipt construction, never its launcher."""
    values = offset_evidence()
    provenance = values["provenance.json"]
    build = provenance["build"]
    for name in ("candidate_path", "candidate_sha256", "candidate_build"):
        del build[name]
    build.update(campaign=c.CAMPAIGN, compiler_path="/usr/bin/x86_64-linux-gnu-g++-13",
                 compiler_sha256=c.COMPILER_HASH,
                 codegen={"disassembly_sha256": "3" * 64, "unwind_sha256": "4" * 64})
    sources = {path: sha for path, sha in build["baseline_build"]["source_files"].items()
               if not path.startswith("bench/")}
    sources.update({path: c.digest(path.encode("ascii")) for path in c.INPUTS})
    sources.update(c.HELPER_HASHES)
    build["source_files"], build["command"] = sources, c.build_command(build)
    for key in ("git_before", "git_after"):
        provenance[key]["required_tracked_inputs"] = {path: c.digest(path.encode("ascii"))[:40] for path in c.INPUTS}
    entries = [{"exists": False, "path": str(path)} for path in c.h.R0_ROOTS]
    for key in ("preserved_before", "preserved_after"):
        provenance[key] = {"entries": entries, "snapshot_sha256": c.digest(c.canonical(entries))}
    claim = {key: build[key] for key in ("baseline_path", "baseline_sha256", "worker_sha256", "source_commit")}
    claim.update(campaign=c.CAMPAIGN, schema=c.PREFIX + ".claim.v1", created_unix_ns=1,
                 description_sha256=c.description_hash(), build_receipt_sha256=c.digest(c.canonical(build)))
    values.update(CLAIM=claim, WORKER_STARTED=c.marker_value(claim))
    values["raw.jsonl"] = raw
    return values


def adjudicate(values, elapsed=0.125):
    return c.adjudicate(values["raw.jsonl"], values["worker.stderr"], values["provenance.json"],
                        values["WORKER_STARTED"], values["CLAIM"], elapsed)


class TranscriptTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.receipt = transcript_fixture()
        cls.raw = encoded(cls.rows)

    def test_complete_native_shaped_transcript_without_live_probes(self):
        with mock.patch.object(c.h, "run_checked", side_effect=AssertionError("live probe")), \
                mock.patch.object(c, "file_hash", side_effect=AssertionError("live artifact")):
            stats, failed = c.parse_transcript(self.raw, self.receipt)
        self.assertEqual(failed, [])
        self.assertEqual(stats["mechanistic_outcome"], "SUPPORTED")
        self.assertEqual(len(self.rows), 2594)

    def test_stack_mapping_warmup_and_retained_geometry_mutations(self):
        for mutation in ("pair", "phase", "warm_stack", "frame", "handle", "warm_time", "measured_time", "sequence"):
            rows = list(self.rows)
            panel = rows[2] = copy.deepcopy(rows[2])
            if mutation == "pair":
                panel["geometry"]["selected_pair"] = [2, 3]
            elif mutation == "phase":
                panel["target_phase"] ^= 16
            elif mutation == "warm_stack":
                stack = panel["warmup"][0]["stack"]
                stack["final_rsp"] = hex(int(stack["final_rsp"], 16) - 16)
            elif mutation == "frame":
                new_frame = hex(int(panel["geometry"]["probe"]["frame_address"], 16) + 16)
                panel["geometry"]["probe"]["frame_address"] = new_frame
                for slot in panel["warmup"] + panel["slots"]:
                    slot["stack"]["frame_address"] = new_frame
            elif mutation == "handle":
                panel["slots"][0]["handle_index"] = 3
            elif mutation == "warm_time":
                panel["warmup"][0]["elapsed_ns"] = 1
            elif mutation == "measured_time":
                panel["slots"][0]["elapsed_ns"] = True
            else:
                panel["sequence"] = True
            with self.subTest(mutation=mutation), self.assertRaises(c.h.ValidationError):
                c.parse_transcript(encoded(rows), self.receipt)

    def test_noncanonical_truncated_binding_oracle_and_terminal_receipts(self):
        for malformed in (self.raw[:-1], self.raw.replace(b'"campaign":', b'"campaign": ', 1), encoded(self.rows[:-2] + self.rows[-1:])):
            with self.subTest(kind="format"), self.assertRaises(c.h.ValidationError):
                c.parse_transcript(malformed, self.receipt)
        for mutation in ("binding", "oracle", "terminal"):
            rows = list(self.rows)
            if mutation == "terminal":
                rows[-1] = dict(rows[-1], encode_call_count=212336639)
            else:
                config = rows[0] = copy.deepcopy(rows[0])
                if mutation == "binding":
                    config["bindings"]["baseline"]["symbols"]["wirehair_v2_encode"]["path"] = "/synthetic/wrong.so"
                else:
                    config["cells"][0]["partial"][0]["source_sha256"] = "f" * 64
            with self.subTest(mutation=mutation), self.assertRaises(c.h.ValidationError):
                c.parse_transcript(encoded(rows), self.receipt)


class BundleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rows, _ = transcript_fixture()
        cls.values = evidence_fixture(encoded(rows))
        cls.anchor = c.digest(c.canonical(cls.values["provenance.json"]["build"]["baseline_build"]))

    def setUp(self):
        # Only the synthetic baseline trust anchor changes; all validators run.
        patch = mock.patch.object(c.o, "BASIS_BUILD_HASH", self.anchor)
        patch.start()
        self.addCleanup(patch.stop)

    def test_complete_bundle_replays_without_live_probes(self):
        values = copy.deepcopy(self.values)
        values["summary.json"] = adjudicate(values)
        self.assertEqual(values["summary.json"]["outcome"], "SUPPORTED")
        with tempfile.TemporaryDirectory(prefix="stack-crossover-replay-test-") as directory:
            root = Path(directory)
            publish_bundle(root, values)
            with mock.patch.object(c, "command", side_effect=AssertionError("live command")), \
                    mock.patch.object(c, "file_hash", side_effect=AssertionError("live artifact")), \
                    mock.patch.object(c.h, "run_worker", side_effect=AssertionError("timing worker")):
                self.assertEqual(c.replay(root), values["summary.json"])
        self.assertFalse(values["summary.json"]["WH1_compared"])
        self.assertFalse(values["summary.json"]["promotion_claimed"])

    def test_worker_outer_deadlines_and_marker_override_supported_statistics(self):
        # Decision precedence only: parsing/statistics have independent tests.
        with mock.patch.object(c, "parse_transcript", return_value=({"mechanistic_outcome": "SUPPORTED"}, ["AA:0:0:0:0"])):
            values = copy.deepcopy(self.values)
            self.assertEqual(adjudicate(values)["outcome"], "CONTROL_FAIL")
            values["worker.stderr"] = b"unexpected diagnostic\n"
            self.assertEqual(adjudicate(values)["outcome"], "INVALID")
        for mutation in ("worker", "outer", "marker", "stderr"):
            values = copy.deepcopy(self.values)
            elapsed = 0.125
            if mutation == "worker":
                values["provenance.json"]["process"]["elapsed_seconds"], elapsed = 60, 61
            elif mutation == "outer":
                elapsed = 70
            elif mutation == "marker":
                values["WORKER_STARTED"]["worker_sha256"] = "f" * 64
            else:
                values["worker.stderr"] = b"unexpected diagnostic\n"
            with self.subTest(mutation=mutation):
                result = adjudicate(values, elapsed)
                self.assertEqual(result["outcome"], "INVALID")
                self.assertTrue(result["infrastructure_failures"])

    def test_rehashed_frozen_helper_mutation_cannot_keep_supported_summary(self):
        values = copy.deepcopy(self.values)
        values["summary.json"] = adjudicate(values)
        build = values["provenance.json"]["build"]
        build["source_files"]["bench/Wh2PublicOffsetR0.py"] = "f" * 64
        values["CLAIM"]["build_receipt_sha256"] = c.digest(c.canonical(build))
        values["WORKER_STARTED"] = c.marker_value(values["CLAIM"])
        with tempfile.TemporaryDirectory(prefix="stack-crossover-tamper-test-") as directory:
            root = Path(directory)
            publish_bundle(root, values)
            with self.assertRaises(c.h.ValidationError):
                c.replay(root)


class GeometryTest(unittest.TestCase):
    def test_first_pair_and_null_use_exact_modular_bytes(self):
        for phases, pre in (((0, 0x10, 0x200, 0x400), 0x800),
                            ((0xef0, 0xf00, 0x110, 0x400), 0x800),
                            ((0, 0x10, 0x200, 0x400), 0xf80),
                            ((0x700, 0x710, 0x200, 0x400), 0x800)):
            addrs, stack = addresses(phases), probe(pre)
            with self.subTest(phases=phases, pre=pre):
                self.assertEqual(c.select_geometry(addrs, stack), expected_geometry(addrs, stack))
        ordinary = c.select_geometry(addresses(), probe())
        self.assertEqual((ordinary["selected_pair"], ordinary["target_phases"]), ([0, 2], [0x110, 0x310, 0]))
        wrapped = c.select_geometry(addresses((0xef0, 0xf00, 0x110, 0x400)), probe())
        self.assertEqual((wrapped["selected_pair"], wrapped["target_phases"]), ([0, 2], [0, 0x220, 0x30]))

    def test_no_pair_or_invalid_handle_and_probe_fail_closed(self):
        with self.assertRaises(c.h.ValidationError):
            c.select_geometry(addresses((0, 0x10, 0x20, 0x10)), probe())
        for mutation in ("duplicate", "zero", "unaligned", "arena", "hot", "final", "correction", "frame"):
            addrs, stack = addresses(), probe()
            if mutation == "duplicate":
                addrs["handles"][1] = addrs["handles"][0]
            elif mutation == "zero":
                addrs["handles"][0] = "0x0"
            elif mutation == "unaligned":
                addrs["handles"][0] = "0x30000001"
            elif mutation == "arena":
                addrs["handles"][0] = addrs["source"]
            elif mutation in ("hot", "final"):
                stack[mutation + "_rsp"] = hex(int(stack["pre_rsp"], 16) - 16)
            elif mutation == "correction":
                stack["correction"] = True
            else:
                stack["frame_address"] = hex(int(stack["pre_rsp"], 16) + 256)
            with self.subTest(mutation=mutation), self.assertRaises(c.h.ValidationError):
                c.select_geometry(addrs, stack)

    def test_measured_stack_matches_probe_frame_and_exact_target(self):
        initial = probe()
        for target in expected_geometry(addresses(), initial)["target_phases"]:
            pre = int(initial["pre_rsp"], 16)
            correction = 4096 + ((pre - target) % 4096)
            stack = dict(initial, hot_rsp=hex(pre - correction), final_rsp=hex(pre - correction), correction=correction)
            c.validate_stack(stack, initial, target)
            for field, value in (("pre_rsp", hex(pre + 16)), ("hot_rsp", hex(pre - correction + 16)),
                                 ("final_rsp", hex(pre - correction - 16)),
                                 ("frame_address", hex(int(initial["frame_address"], 16) + 16)),
                                 ("correction", True), ("correction", correction + 4096)):
                with self.subTest(target=target, field=field), self.assertRaises(c.h.ValidationError):
                    c.validate_stack(dict(stack, **{field: value}), initial, target)


class StatisticsTest(unittest.TestCase):
    def test_exact_roster_balance_and_independent_work_counts(self):
        roster = list(c.roster())
        self.assertEqual(roster, list(expected_roster()))
        self.assertEqual(tuple(map(tuple, c.CELLS)), CELLS)
        self.assertEqual(tuple(map(tuple, c.CONDITIONS)), CONDITIONS)
        self.assertEqual((len(roster), len(set(roster))), (2592, 2592))
        self.assertEqual(8 * 8192 * len(roster), 169869312)
        self.assertEqual(2 * 8192 * len(roster), 42467328)
        self.assertEqual(sum(8 * (8192 // CELLS[cell][0]) for _, cell, _, _, _ in roster), 7527168)
        self.assertEqual(sum(2 * (8192 // CELLS[cell][0]) for _, cell, _, _, _ in roster), 1881792)
        self.assertEqual(str(c.OUTPUT), "/var/tmp/wh2-public-stack-crossover-r0")

    def test_signed_target_minus_null_and_replicate_pairing(self):
        stats, failed = c.statistics(samples())
        self.assertEqual(failed, [])
        self.assertEqual(stats["mechanistic_outcome"], "SUPPORTED")
        self.assertEqual((len(stats["controls"]), len(stats["deltas"])), (144, 48))
        self.assertEqual(sum(row["primary"] for row in stats["deltas"]), 16)
        for row in stats["controls"] + stats["deltas"]:
            self.assertEqual(row["summary"]["n"], 12)
        for row in stats["deltas"]:
            self.assertEqual(row["primary"], row["cell"] in (3, 5))
            self.assertAlmostEqual(row["summary"]["log_mean"], 0.04 if row["target"] == 0 else -0.04, places=14)
            self.assertAlmostEqual(row["summary"]["log_standard_error"], 0, places=14)

    def test_only_frozen_primary_cells_determine_mechanistic_support(self):
        values = samples()
        for key in values:
            rep, cell, target, comparison, condition = key
            if cell not in (3, 5) and comparison == 2 and target < 2:
                values[key] += (-0.2 if target == 0 else 0.2) * (1 if condition < 2 else -1)
        stats, failed = c.statistics(values)
        self.assertEqual(failed, [])
        self.assertEqual(stats["mechanistic_outcome"], "SUPPORTED")
        self.assertTrue(any(row["summary"]["log_mean"] < 0 for row in stats["deltas"]
                            if not row["primary"] and row["target"] == 0))

    def test_conditional_aa_failure_cannot_be_averaged_away(self):
        values = samples()
        for key in values:
            rep, cell, target, comparison, condition = key
            if cell == 3 and target == 0 and comparison == 0 and condition in (0, 2):
                values[key] = 0.03 if condition == 0 else -0.03
        _, failed = c.statistics(values)
        self.assertEqual(set(failed), {"AA:3:0:0:0", "AA:3:0:0:2"})

    def test_within_two_percent_is_distinct_from_inconclusive(self):
        for delta0, delta1, expected in ((0, 0, "WITHIN_2PCT"),
                                         (0.005, 0.005, "WITHIN_2PCT"),
                                         (0.04, 0.04, "INCONCLUSIVE"),
                                         (0, -0.04, "INCONCLUSIVE")):
            with self.subTest(delta0=delta0, delta1=delta1):
                stats, failed = c.statistics(samples(delta0, delta1, common_drift=False))
                self.assertEqual(failed, [])
                self.assertEqual(stats["mechanistic_outcome"], expected)

    def test_nonfinite_or_boolean_samples_fail_closed(self):
        for key in ((0, 0, 0, 0, 0), (0, 0, 0, 2, 0)):
            for value in (True, float("nan"), float("inf")):
                values = samples()
                values[key] = value
                with self.subTest(key=key, value=value), self.assertRaises(c.h.ValidationError):
                    c.statistics(values)
        values = samples()
        del values[0, 0, 0, 2, 0]
        with self.assertRaises(c.h.ValidationError):
            c.statistics(values)

    def test_statistics_do_not_probe_live_state(self):
        with mock.patch.object(c.h, "run_checked", side_effect=AssertionError("live probe")):
            self.assertEqual(c.statistics(samples())[1], [])


if __name__ == "__main__":
    unittest.main()
