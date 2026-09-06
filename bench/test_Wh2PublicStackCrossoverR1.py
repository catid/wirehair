#!/usr/bin/env python3
"""Synthetic R1 adapter regressions; no codec, timing worker or namespace claim."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest import mock

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
import Wh2PublicStackCrossoverR1 as adapter
import Wh2PublicStackCrossoverR0 as frozen

# Private test-module instance: adapt fixtures without modifying either frozen
# source or the normal R0 test module's globals in a combined discovery run.
spec = importlib.util.spec_from_file_location("_stack_r1_frozen_tests",
    BENCH / "test_Wh2PublicStackCrossoverR0.py")
t = importlib.util.module_from_spec(spec)
spec.loader.exec_module(t)
c = adapter.c
t.c = c

TranscriptTest = t.TranscriptTest
BundleTest = t.BundleTest
GeometryTest = t.GeometryTest


class StatisticsTest(t.StatisticsTest):
    def test_exact_roster_balance_and_independent_work_counts(self):
        roster = list(c.roster())
        self.assertEqual(roster, list(t.expected_roster()))
        self.assertEqual(tuple(map(tuple, c.CELLS)), t.CELLS)
        self.assertEqual(tuple(map(tuple, c.CONDITIONS)), t.CONDITIONS)
        self.assertEqual((len(roster), len(set(roster))), (2592, 2592))
        self.assertEqual(8 * 8192 * len(roster), 169869312)
        self.assertEqual(2 * 8192 * len(roster), 42467328)
        self.assertEqual(sum(8 * (8192 // t.CELLS[cell][0])
            for _, cell, _, _, _ in roster), 7527168)
        self.assertEqual(sum(2 * (8192 // t.CELLS[cell][0])
            for _, cell, _, _, _ in roster), 1881792)
        self.assertEqual(str(c.OUTPUT), "/var/tmp/wh2-public-stack-crossover-r1")


class AdapterTest(unittest.TestCase):
    def test_private_namespace_and_unchanged_scientific_description(self):
        self.assertIsNot(c, frozen)
        self.assertIsNot(c.h, frozen.h)
        self.assertEqual(frozen.CAMPAIGN, "wh2-public-stack-crossover-r0")
        self.assertEqual(str(frozen.OUTPUT), "/var/tmp/wh2-public-stack-crossover-r0")
        expected = frozen.description()
        expected.update(campaign=c.CAMPAIGN, schema=c.PREFIX + ".describe.v1")
        self.assertEqual(c.description(), expected)
        self.assertNotEqual(c.description_hash(), frozen.description_hash())
        self.assertEqual(list(c.roster()), list(frozen.roster()))
        self.assertEqual(c.statistics(t.samples()), frozen.statistics(t.samples()))
        self.assertEqual(c.h.CAMPAIGN, c.CAMPAIGN)
        self.assertEqual(c.h.FIXED_OUTPUT_DIR, c.OUTPUT)
        self.assertEqual(c.h.R1_TRACKED_INPUTS, c.INPUTS)
        self.assertIn(frozen.OUTPUT, c.h.R0_ROOTS)
        self.assertNotIn(c.OUTPUT, c.h.R0_ROOTS)
        self.assertEqual(len(c.h.R0_ROOTS), len(set(c.h.R0_ROOTS)))

    def test_literal_frozen_pins_and_complete_additive_source_closure(self):
        for path, sha in c.HELPER_HASHES.items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((c.ROOT / path).read_bytes()).hexdigest(), sha)
        for path in ("bench/Wh2InstalledCodeArtifact.h", "bench/Wh2InstalledCodeArtifact.cpp",
                     "bench/Wh2InstalledCodeArtifactTest.cpp", "bench/Wh2PublicStackCrossoverR1.cpp",
                     "bench/Wh2PublicStackCrossoverR1.py", "bench/test_Wh2PublicStackCrossoverR1.py"):
            self.assertIn(path, c.INPUTS)
        self.assertTrue(set(frozen.INPUTS).issubset(c.INPUTS))
        self.assertEqual(len(c.INPUTS), len(set(c.INPUTS)))
        command = c.build_command({"compiler_path": "/synthetic/g++", "source_commit": "a" * 40,
                                   "worker_path": "/synthetic/worker"})
        self.assertIn(str(c.ROOT / "bench/Wh2InstalledCodeArtifact.cpp"), command)
        self.assertIn(str(c.ROOT / "bench/Wh2PublicStackCrossoverR1.cpp"), command)
        self.assertNotIn(str(c.ROOT / "bench/Wh2PublicStackCrossoverR0.cpp"), command)

    def test_installed_code_probe_is_bound_into_private_live_build_checks(self):
        self.assertIs(c.check_build, adapter.check_build)
        receipt = {"worker_path": "/synthetic/worker"}
        with mock.patch.object(adapter, "_frozen_check_build") as check, \
                mock.patch.object(c, "command", return_value=
                    "stack-crossover-r1 installed libc verified (no codec or clocks)\n") as command:
            c.check_build(receipt)
            check.assert_called_once_with(receipt, live=True)
            command.assert_called_once_with(["/synthetic/worker", "--artifact-preflight"])
            command.reset_mock()
            c.check_build(receipt, live=False)
            command.assert_not_called()
        with mock.patch.object(adapter, "_frozen_check_build"), \
                mock.patch.object(c, "command", return_value="wrong result\n"):
            with self.assertRaises(c.h.ValidationError):
                c.check_build(receipt)
        with mock.patch.object(adapter, "_frozen_check_build", side_effect=c.h.ValidationError("bad build")), \
                mock.patch.object(c, "command") as command:
            with self.assertRaises(c.h.ValidationError):
                c.check_build(receipt)
            command.assert_not_called()

    def test_both_native_libc_checks_are_typed_and_evidence_stays_strict(self):
        source = (BENCH / "Wh2PublicStackCrossoverR1.cpp").read_text()
        self.assertNotIn("old::FileHash(Libc)", source)
        check = "wirehair::wh2_benchmark::VerifyInstalledCodeArtifact(Libc,LibcHash);"
        self.assertEqual(source.count(check), 3)  # entry, terminal, read-only preflight
        authorization = source.split("std::string Authorize(", 1)[1].split("struct MockContext", 1)[0]
        self.assertLess(authorization.index(check), authorization.index("Fd marker("))
        self.assertIn('old::FileHash("/proc/self/exe")==workerHash', authorization)
        self.assertIn("old::ReadFd(claimFd.Value,1024*1024,0400)", authorization)
        terminal = source.split('Require(sequence==2592', 1)[1].split("std::string Authorize(", 1)[0]
        self.assertLess(terminal.index(check), terminal.index("old::Emit(Terminal())"))

    def test_rehashed_stack_r0_helper_mutation_is_rejected(self):
        rows, _ = t.transcript_fixture()
        values = t.evidence_fixture(t.encoded(rows))
        anchor = c.digest(c.canonical(values["provenance.json"]["build"]["baseline_build"]))
        with mock.patch.object(c.o, "BASIS_BUILD_HASH", anchor):
            for path in ("bench/Wh2PublicStackCrossoverR0.cpp", "bench/Wh2PublicStackCrossoverR0.py",
                         "bench/test_Wh2PublicStackCrossoverR0.py"):
                build = copy.deepcopy(values["provenance.json"]["build"])
                build["source_files"][path] = "f" * 64
                with self.subTest(path=path), self.assertRaises(c.h.ValidationError):
                    c.check_build(build, live=False)


if __name__ == "__main__":
    unittest.main()
