#!/usr/bin/env python3
"""Additive broad-pilot adapter tests: synthetic inputs, never codec campaigns."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2ProductionMix3RecoveryR0 as original

ORIGINAL_STATE = (original.PROTOCOL, original.TRAINING_ROOTS, original.HOLDOUT_ROOTS,
                  original.source_receipt, original.check_build, original.frozen_protocol)
import Wh2ProductionMix3RecoveryBroadR0 as broad

c = broad.c
fixture_spec = importlib.util.spec_from_file_location(
    "_wh2_broad_fixture_test_private", Path(__file__).with_name("test_Wh2ProductionMix3RecoveryR0.py"))
fixtures = importlib.util.module_from_spec(fixture_spec)
fixture_spec.loader.exec_module(fixtures)
fixtures.c = c
fixtures.LIBRARY_PATH = str(broad.PINNED_LIBRARY)
OLD_BUNDLE = Path("/var/tmp/wh2-production-mix3-k3k6-r0")


def synthetic_build():
    sources = {path: c.digest(path.encode("ascii")) for path in broad.SOURCE_PATHS}
    basis = {"source_commit": broad.PRODUCTION_COMMIT, "library": str(broad.PINNED_LIBRARY),
             "library_sha256": broad.LIBRARY_SHA256,
             "source_files": {path: value for path, value in sources.items() if "BroadR0" not in path}}
    receipt = {"protocol": broad.PROTOCOL, "source_commit": fixtures.SOURCE_COMMIT,
               "source_files": sources, "production_build": basis,
               "production_source_commit": broad.PRODUCTION_COMMIT,
               "library": str(broad.PINNED_LIBRARY), "library_sha256": broad.LIBRARY_SHA256,
               "worker": "/synthetic/worker", "worker_sha256": "2" * 64,
               "compiler": "/synthetic/compiler", "compiler_sha256": "3" * 64}
    live = {str(broad.ROOT / path): value for path, value in sources.items()}
    live.update({receipt[key]: receipt[key + "_sha256"] for key in ("library", "worker", "compiler")})
    return receipt, live


class BroadRecoveryAdapterTest(unittest.TestCase):
    def test_private_module_and_fixture_isolation(self):
        self.assertIsNot(c, original)
        self.assertEqual((original.PROTOCOL, original.TRAINING_ROOTS, original.HOLDOUT_ROOTS,
                          original.source_receipt, original.check_build, original.frozen_protocol),
                         ORIGINAL_STATE)
        self.assertIs(c.parse_phase.__globals__, c.__dict__)
        self.assertIs(c.validate_cell.__globals__, c.__dict__)
        self.assertIs(fixtures.transcript.__globals__["c"], c)
        self.assertEqual(c.PROTOCOL, broad.PROTOCOL)
        self.assertEqual(c.TRAINING_ROOTS, broad.TRAINING_ROOTS)
        self.assertEqual(c.HOLDOUT_ROOTS, broad.HOLDOUT_ROOTS)

    def test_exact_roots_counts_and_disjointness(self):
        self.assertEqual(broad.PROTOCOL, "wirehair.wh2.production-mix3-k3k6-broad-r0")
        roots = {phase: tuple("0x" + hashlib.sha256(
            (broad.PROTOCOL + ":" + phase + "/" + str(i)).encode("ascii")).hexdigest()[:16]
            for i in range(count)) for phase, count in (("train", 16), ("holdout", 64))}
        self.assertEqual(broad.TRAINING_ROOTS, roots["train"])
        self.assertEqual(broad.HOLDOUT_ROOTS, roots["holdout"])
        self.assertEqual(c.digest(c.canonical(roots)),
                         "2a8ae1333952f4cbd28e4af15cdd689e2314506b18526c377f0a5bb6af79e202")
        all_roots = set(roots["train"] + roots["holdout"])
        self.assertEqual(len(all_roots), 80)
        reserved = {"0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3",
                    "0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b",
                    "0xefd20c982041a46b", "0x8827bc36ed906555", "0x86029f23d6132efa"}
        self.assertFalse(all_roots & (reserved | set(original.TRAINING_ROOTS + original.HOLDOUT_ROOTS)))
        self.assertEqual((c.KS, c.SCHEDULES, c.OVERHEADS),
                         ((3, 6), ("burst", "adversarial", "repair-only"), (0, 1, 2, 4)))
        self.assertEqual((2 * 257 * 16 * 3 + 2 * 2 * 64 * 3) * 4, 101760)

    def test_fresh_roster_and_paired_summary_sizes(self):
        train = fixtures.parse(fixtures.transcript())
        holdout = fixtures.parse(fixtures.transcript("holdout"), "holdout", [1, 2])
        self.assertEqual(len(train["cells"]), 336)
        self.assertEqual(len(holdout["cells"]), 768)
        summary = c.summarize(train, holdout)
        self.assertEqual(summary["outcome"], "PASS")
        self.assertEqual(summary["training"]["cells_per_arm"], 96)
        self.assertEqual(summary["holdout"]["cells_per_arm"], 384)
        self.assertEqual(summary["training"]["repairs"], [2, 0, 0, 0])

    def test_maximum_attempt_roster_fits_unchanged_parser_cap(self):
        # Pure synthetic rows: both first successes occur at255, so every
        # permitted attempt appears without executing any native codec.
        records = fixtures.transcript(selected=(255, 255))
        self.assertEqual(len(records) - 2, 24672)
        text = fixtures.encoded(records)
        self.assertLessEqual(len(text.encode("ascii")), 24279194)
        self.assertLess(len(text), 32 * 1024 * 1024)
        parsed = c.parse_phase(text, "train", fixtures.SOURCE_COMMIT, fixtures.LIBRARY_PATH)
        self.assertEqual(parsed["selected_attempts"], [255, 255])

    def test_inherited_first_success_map_and_construction_invariants(self):
        for name in ("test_first_success_cannot_be_skipped_or_followed_by_more_attempts",
                     "test_candidate_holdout_failure_does_not_reselect",
                     "test_cross_phase_maps_cannot_change",
                     "test_actual_baseline_attempt_is_not_assumed_zero",
                     "test_bad_seed_cannot_vary_by_loss_coordinate",
                     "test_same_attempt_semantics_must_agree_between_arms"):
            with self.subTest(invariant=name):
                getattr(fixtures.ProductionMix3RecoveryTest, name)(self)

    def test_harness_status_is_scoped_but_head_checks_are_unchanged(self):
        status = ["git", "status", "--porcelain", "--untracked-files=no"]
        with mock.patch.object(broad, "_original_command", return_value="") as command:
            broad.harness_command(status)
            command.assert_called_once_with(status + ["--"] + list(broad.SOURCE_PATHS), 30)
            command.reset_mock()
            broad.harness_command(["git", "rev-parse", "HEAD"])
            command.assert_called_once_with(["git", "rev-parse", "HEAD"], 30)
        self.assertNotIn("codec/WirehairV2Profile.cpp", broad.SOURCE_PATHS)

    def test_pinned_build_rejects_artifact_and_coherent_inherited_source_drift(self):
        receipt, live = synthetic_build()
        basis_hash = c.digest(c.canonical(receipt["production_build"]))
        with mock.patch.object(broad, "BASIS_BUILD_SHA256", basis_hash), \
                mock.patch.object(c, "command", return_value=fixtures.SOURCE_COMMIT + "\n"), \
                mock.patch.object(c, "file_digest", side_effect=lambda path: live[str(path)]):
            broad.check_build(receipt)
            # Current wrapper HEAD and the preserved production source basis are
            # deliberately different identities, not interchangeable labels.
            self.assertNotEqual(receipt["source_commit"], receipt["production_source_commit"])
            for key in ("library", "worker", "compiler"):
                path, saved = receipt[key], live[receipt[key]]
                live[path] = "f" * 64
                with self.subTest(artifact=key), self.assertRaises(c.ValidationError):
                    broad.check_build(receipt)
                live[path] = saved
            for path in receipt["production_build"]["source_files"]:
                forged = copy.deepcopy(receipt)
                forged["source_files"][path] = "f" * 64
                absolute = str(broad.ROOT / path)
                saved, live[absolute] = live[absolute], "f" * 64
                with self.subTest(inherited=path), self.assertRaises(c.ValidationError):
                    broad.check_build(forged)
                live[absolute] = saved
            for key in ("production_source_commit", "library_sha256"):
                forged = copy.deepcopy(receipt)
                forged[key] = "f" * len(forged[key])
                with self.subTest(identity=key), self.assertRaises(c.ValidationError):
                    broad.check_build(forged)
            forged = copy.deepcopy(receipt)
            forged["production_build"]["source_commit"] = "f" * 40
            with self.assertRaises(c.ValidationError):
                broad.check_build(forged)

    def test_pinned_library_receipt_checks_freeze_and_library_bytes(self):
        receipt, _ = synthetic_build()
        basis = receipt["production_build"]
        hashes = {str(broad.BASIS_FREEZE): broad.BASIS_FREEZE_SHA256,
                  str(broad.PINNED_LIBRARY): broad.LIBRARY_SHA256}
        with mock.patch.object(broad, "BASIS_BUILD_SHA256", c.digest(c.canonical(basis))), \
                mock.patch.object(Path, "read_bytes", return_value=c.canonical({"build": basis})), \
                mock.patch.object(c, "file_digest", side_effect=lambda path: hashes[str(path)]):
            self.assertEqual(broad.pinned_library_receipt(), basis)
            for path in hashes:
                saved, hashes[path] = hashes[path], "f" * 64
                with self.subTest(path=path), self.assertRaises(c.ValidationError):
                    broad.pinned_library_receipt()
                hashes[path] = saved

    def test_broad_freeze_offline_replay_and_publication_deadline(self):
        receipt, _ = synthetic_build()
        original_fixture = fixtures.bundle_fixture

        def bundle(selected=(1, 2)):
            # Adapt only the fixture's minimal build receipt. Production code
            # continues to validate the full embedded basis during replay.
            with mock.patch.object(c, "frozen_protocol", side_effect=lambda unused: broad.frozen_protocol(receipt)):
                return original_fixture(selected)

        with mock.patch.object(broad, "BASIS_BUILD_SHA256", c.digest(c.canonical(receipt["production_build"]))), \
                mock.patch.object(fixtures, "bundle_fixture", side_effect=bundle):
            freeze = broad.frozen_protocol(receipt)
            self.assertEqual(freeze["max_prefix_decodes"], 101760)
            self.assertEqual(freeze["production_source_commit"], broad.PRODUCTION_COMMIT)
            self.assertEqual(freeze["roster_sha256"], broad.ROSTER_SHA256)
            for name in ("test_offline_replay_never_calls_live_tools",
                         "test_replay_rejects_rehashed_contract_mutations",
                         "test_replay_rejects_rehashed_errors_rosters_maps_and_summaries",
                         "test_publication_hashing_cannot_seal_success_after_deadline"):
                with self.subTest(invariant=name):
                    getattr(fixtures.ProductionMix3RecoveryTest, name)(self)

    @unittest.skipUnless(OLD_BUNDLE.exists(), "immutable .52 bundle unavailable on this host")
    def test_old_bundle_still_replays_without_modification_or_live_tools(self):
        before = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in OLD_BUNDLE.iterdir()}
        with mock.patch.object(original, "command", side_effect=AssertionError("live old replay")):
            result = original.replay(OLD_BUNDLE)
        self.assertEqual(result["outcome"], "FAIL")
        self.assertEqual(result["selected_attempts"], [5, 3])
        self.assertEqual(before, {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in OLD_BUNDLE.iterdir()})


if __name__ == "__main__":
    unittest.main()
