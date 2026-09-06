#!/usr/bin/env python3
"""Neutral data packaging/authentication tests; never builds the actual header."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import stat
import struct
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class NativeDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = sibling("_native_data_under_test", "Wh2ThueMorseNativeDataR0.py")
        cls.math = sibling("_native_data_neutral_math", "Wh2ThueMorseR0.py")
        cls.oracle = sibling("_native_data_independent_oracle", "test_Wh2ThueMorseR0.py")

    def setUp(self):
        # Even accidental entry into an actual campaign or recorded artifacts
        # aborts every neutral test. Generation tests replace these with fakes.
        forbidden = ((self.data, "authenticated_inputs"), (self.data, "recovery_bundle"),
                     (self.data.H, "load_inputs"), (self.math, "fixed_feedback"),
                     (self.math, "select_parameter"), (self.math, "run_screen"),
                     (self.math.RowMapper, "row"), (self.math.RowMapper, "reference_row"),
                     (self.math.RADIX, "rank_columns"), (self.math.RADIX, "candidate_matrix"),
                     (self.math.RADIX, "run_screen"), (self.oracle, "geometric_feedback"),
                     (self.oracle, "selection_evidence"), (self.oracle, "verify_result"))
        for module, name in forbidden:
            patch = mock.patch.object(module, name, side_effect=AssertionError("actual corpus/scoring forbidden"))
            patch.start()
            self.addCleanup(patch.stop)
        patch = mock.patch.object(self.data, "sibling", return_value=self.math)
        patch.start()
        self.addCleanup(patch.stop)

    def neutral_report(self):
        pair = self.oracle.neutral_companions()
        packed, lookup = self.data.build_lookup(pair)
        ids = list(range(10))
        trace_hash = hashlib.sha256(struct.pack("<10I", *ids)).hexdigest()
        common = dict(ids=ids, attempted_candidates=10, trace_sha256=trace_hash, ranks=[6] * 5)
        fresh = [dict(root_index=i, B=B, schedule=schedule, **common)
                 for i in range(512) for B in self.data.WIDTHS for schedule in self.data.SCHEDULES]
        hard = [dict(phase=phase, root_index=i, B=2, schedule=schedule, **common)
                for phase in ("training", "validation") for i in range(3) for schedule in self.data.SCHEDULES[1:]]
        prefixes = [dict(K=6, ids=list(range(6 + int(i >= 78))), overhead=int(i >= 78),
                         original_widths=[2, 64] if i < 23 else [2]) for i in range(82)]
        fixtures = [dict(ids=list(range(6 if i < 16 else 7 if i == 16 else 8)),
                         expected_rank=5 if i < 17 else 6) for i in range(33)]
        inputs = dict(pair=dict(A0=pair[0], A1=pair[1]), prefixes=prefixes, fixtures=fixtures,
                      provenance={"neutral_historical": True})
        report = dict(protocol="wirehair.wh2.thue-morse-recovery-r0", outcome="PASS", inputs=inputs,
                      inputs_sha256=self.data.digest(self.data.canonical(inputs)), fresh=fresh, hard=hard,
                      history=[dict(index=i, rank=6, repaired=True) for i in range(82)],
                      fixtures=[dict(index=i, rank=row["expected_rank"]) for i, row in enumerate(fixtures)],
                      evidence=dict(lookup=lookup, coefficient_visit_sha256="1" * 64))
        return report, packed, lookup

    def test_neutral_tables_match_independent_oracle_without_row_scoring(self):
        pair = self.oracle.neutral_companions()
        packed, actual = self.data.build_lookup(pair)
        expected = self.oracle.lookup_evidence(pair)
        self.assertEqual(self.data.canonical(actual), self.data.canonical(expected))
        self.assertEqual(len(packed), 39936)
        self.assertEqual(self.data.digest(packed), actual["tables_sha256"])
        offsets = (0, 6144, 12288, 16896, 21504, 26112, 30720, 39936)
        self.assertEqual([self.data.digest(packed[a:b]) for a, b in zip(offsets, offsets[1:])],
                         [row["sha256"] for row in actual["tables"]])

    def test_complete_projection_counts_tail_order_and_no_duplicate_visit_array(self):
        report, packed, _ = self.neutral_report()
        data = self.data.extract_data(report)
        self.assertEqual(len(data["traces"]), 6162)
        self.assertEqual(len(data["history_cases"]), 105)
        self.assertEqual(len(data["fixture_cases"]), 99)
        self.assertEqual(len(data["partial_cases"]), 90)
        self.assertEqual(data["partial_cases"][:5], [[6144, 2, 1], [6144, 64, 1], [6144, 64, 63],
                                                    [6144, 1280, 1], [6144, 1280, 1279]])
        self.assertEqual(data["partial_cases"][-1], [6161, 1280, 1279])
        self.assertEqual(sum(len(row["ids"]) for group in ("traces", "history", "fixtures") for row in data[group]),
                         62347)
        output = self.data.render_header(data, packed)
        self.assertEqual(output, self.data.render_header(self.data.extract_data(copy.deepcopy(report)), packed))
        self.assertLess(len(output), 4 * 1024 * 1024)
        self.assertIn(b"alignas(64) extern const uint8_t kLookup[39936]", output)
        self.assertIn(b"extern const Trace kTraces[kTraceCount]", output)
        self.assertIn(b"extern const PayloadCase kPartialCases[90]", output)
        self.assertNotIn(b"kVisitIds", output)
        self.assertIn(self.data.digest(self.data.canonical(data)).encode("ascii"), output)

    def test_typed_projection_rejects_incomplete_order_hash_rank_width_mutations(self):
        report, _, _ = self.neutral_report()
        mutations = (
            (("fresh",), report["fresh"][:-1]), (("hard",), report["hard"][:-1]),
            (("history",), report["history"][:-1]), (("fixtures",), report["fixtures"][:-1]),
            (("fresh", 0, "B"), 2.0), (("fresh", 0, "root_index"), False),
            (("fresh", 0, "ids", 0), True), (("fresh", 0, "ids", 1), 0),
            (("fresh", 0, "ranks", 0), 6.0), (("fresh", 0, "ranks"), [5, 4, 6, 6, 6]),
            (("fresh", 0, "trace_sha256"), "0" * 64),
            (("fresh", 0, "attempted_candidates"), True),
            (("hard", 0, "phase"), "validation"),
            (("inputs", "prefixes", 0, "K"), 6.0),
            (("inputs", "prefixes", 0, "original_widths"), [2, 2]),
            (("inputs", "prefixes", 0, "original_widths"), [2.0, 64]),
            (("history", 0, "repaired"), 1), (("fixtures", 0, "rank"), 5.0))
        for path, value in mutations:
            wrong = copy.deepcopy(report)
            target = wrong
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises((ValueError, TypeError)):
                self.data.extract_data(wrong)

    def test_create_only_immutable_neutral_output_and_receipt(self):
        report, packed, lookup = self.neutral_report()
        provenance = {"neutral_recovery": True}
        with tempfile.TemporaryDirectory(prefix="wh2-native-data-neutral-") as directory, \
                mock.patch.object(self.data, "authenticated_inputs", return_value=report), \
                mock.patch.object(self.data, "recovery_bundle", return_value=(report, provenance)), \
                mock.patch.object(self.data.H, "load_inputs", return_value=report["inputs"]), \
                mock.patch.object(self.data, "build_lookup", return_value=(packed, lookup)):
            header, receipt = Path(directory) / "neutral.h", Path(directory) / "neutral.json"
            actual = self.data.generate(str(header), str(receipt))
            for path in (header, receipt):
                info = path.stat()
                self.assertTrue(stat.S_ISREG(info.st_mode))
                self.assertEqual(stat.S_IMODE(info.st_mode), 0o400)
                self.assertEqual(info.st_nlink, 1)
            raw = header.read_bytes()
            self.assertEqual(actual["header_sha256"], self.data.digest(raw))
            self.assertEqual(actual["header_bytes"], len(raw))
            self.assertEqual(receipt.read_bytes(), self.data.canonical(actual))
            self.assertEqual(actual["recovery_provenance"], provenance)
            self.assertEqual(actual["counts"]["coefficient_visits"], 62347)
            self.assertNotIn("timestamp", actual)
            self.assertNotIn("path", actual)
            with self.assertRaises(ValueError):
                self.data.generate(str(header), str(receipt))
            self.assertEqual(header.read_bytes(), raw)

    def test_output_collision_or_symlink_never_enters_authenticated_generation(self):
        with tempfile.TemporaryDirectory(prefix="wh2-native-data-path-") as directory:
            root = Path(directory)
            header, receipt = root / "neutral.h", root / "neutral.json"
            receipt.symlink_to(root / "absent")
            with self.assertRaises(ValueError):
                self.data.generate(str(header), str(receipt))
            self.assertFalse(header.exists())
            with self.assertRaises(ValueError):
                self.data.generate("relative.h", str(root / "other.json"))
            with self.assertRaises(ValueError):
                self.data.generate(str(header), str(header))
            alias = root / "alias"
            alias.symlink_to(root, target_is_directory=True)
            with self.assertRaises(ValueError):
                self.data.generate(str(alias / "new.h"), str(alias / "new.json"))
            self.assertFalse((root / "new.h").exists())

    def test_bundle_authentication_binds_every_member_and_canonical_bytes(self):
        report = {"neutral_recorded": True}
        files = {"raw.json": self.data.canonical(report) + b"\n", "CLAIM.json": b"{}",
                 "summary.json": b"{}", "stderr.txt": b""}
        members = {name: dict(bytes=len(raw), sha256=self.data.digest(raw)) for name, raw in files.items()}
        manifest = self.data.canonical(dict(protocol="wirehair.wh2.thue-morse-recovery-r0", outcome="PASS", files=members))
        files["COMPLETE.json"] = manifest
        original = self.data.recovery_bundle
        # Undo only this neutral test's guard by loading the same function from
        # a separately inert module; no authentic filesystem read is permitted.
        fresh = sibling("_native_data_bundle_neutral", "Wh2ThueMorseNativeDataR0.py")
        with mock.patch.object(fresh, "MANIFEST_SHA", self.data.digest(manifest)), \
                mock.patch.object(fresh, "RAW_SHA", self.data.digest(files["raw.json"])), \
                mock.patch.object(fresh.C, "read_regular", side_effect=lambda path, cap, **kwargs: files[path.name]), \
                mock.patch.object(fresh.os, "scandir") as scan:
            scan.return_value.__enter__.return_value = [SimpleNamespace(name=name) for name in files]
            actual, provenance = fresh.recovery_bundle()
            self.assertEqual(actual, report)
            self.assertEqual(provenance["files"], members)
            files["summary.json"] = b'{"tampered":true}'
            with self.assertRaises(ValueError):
                fresh.recovery_bundle()
            files["summary.json"] = b"{}"
            scan.return_value.__enter__.return_value.append(SimpleNamespace(name="unexpected"))
            with self.assertRaises(ValueError):
                fresh.recovery_bundle()
        self.assertIs(original, self.data.recovery_bundle)

    def test_cmake_neutral_target_has_no_actual_header_dependency(self):
        source = (HERE.parent / "CMakeLists.txt").read_text()
        self.assertIn('option(WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0\n', source)
        self.assertIn('"Build the frozen benchmark-only K6 Thue-Morse native R0 screen" OFF)', source)
        neutral = source.split('add_executable(wh2_thue_morse_native_r0_test EXCLUDE_FROM_ALL', 1)[1]
        neutral = neutral.split('if (BUILD_TESTS)', 1)[0]
        self.assertNotIn('WH2_THUE_NATIVE_HEADER', neutral)
        self.assertNotIn('NativeDataR0', neutral)
        self.assertNotIn('wirehair_v2_policy', neutral)
        self.assertIn('target_link_libraries(wh2_thue_morse_native_r0_test PRIVATE wirehair)', neutral)


if __name__ == "__main__":
    unittest.main()
