#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

try:
    from . import fixup_holdout_ledger as ledger
except ImportError:
    import fixup_holdout_ledger as ledger


class FixupHoldoutLedgerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def fixups(self, name, rows):
        path = self.root / name
        path.write_text(
            "// generated\n" + "".join(
                f"    {{ {n}, {seed} }},\n" for n, seed in rows
            ),
            encoding="ascii",
        )
        return ledger.read_fixups(path)

    def ohead(
        self, name, loss, rows, policy="disabled", trials=10000,
        metadata_count=None, seed="0x1001", bb=64, start_mode=0,
    ):
        path = self.root / name
        policy_text = f" seed_fixups={policy}" if policy is not None else ""
        seed_text = f" seed={seed}" if seed is not None else ""
        count = len(rows) if metadata_count is None else metadata_count
        path.write_text(
            f"# ohead: threads=4 trials/N={trials} "
            f"loss={loss:.2f}{seed_text}{policy_text} Ns={count} "
            f"bb={bb} startMode={start_mode} range[2,3]\n"
            "N mean p50 p99 p999 max fail\n" +
            "".join(
                f"{n} {mean:.17g} 0 1 1 1 {fail}\n"
                for n, mean, fail in rows
            ),
            encoding="ascii",
        )
        return path

    def test_remove_retain_and_retune_decisions(self):
        peel = self.fixups("peel.inc", [(2, 10)])
        dense = self.fixups("dense.inc", [(3, 20)])
        base10_path = self.ohead("b10", 0.10, [(2, 0.01, 0), (3, 0.08, 0)])
        base30_path = self.ohead("b30", 0.30, [(2, 0.10, 0), (3, 0.20, 0)])
        fixed10_path = self.ohead(
            "f10", 0.10, [(2, 0.01, 0), (3, 0.02, 0)], policy="enabled"
        )
        fixed30_path = self.ohead(
            "f30", 0.30, [(2, 0.10, 0), (3, 0.12, 0)], policy="enabled"
        )
        base10, *_ = ledger.read_ohead(base10_path, 0.10, 10000)
        base30, *_ = ledger.read_ohead(base30_path, 0.30, 10000)
        fixed10, *_ = ledger.read_ohead(fixed10_path, 0.10, 10000)
        fixed30, *_ = ledger.read_ohead(fixed30_path, 0.30, 10000)
        rows = ledger.build_ledger(
            peel, dense, base10, base30, fixed10, fixed30
        )
        self.assertEqual(
            [row["Decision"] for row in rows],
            ["remove-exact-fixup", "retain-exact-fixup"],
        )
        self.assertEqual(rows[1]["FixedMean10"], 0.02)
        self.assertEqual(rows[1]["FixedMean30"], 0.12)
        self.assertEqual(rows[1]["FixedFail10"], 0)
        self.assertEqual(rows[1]["FixedFail30"], 0)
        fixed30[3] = (0.20, 1)
        rows = ledger.build_ledger(
            peel, dense, base10, base30, fixed10, fixed30
        )
        self.assertEqual(rows[1]["Decision"], "retune-required")

    def test_changed_passing_candidate_is_retune(self):
        baseline_peel = {2: 10}
        candidate_peel = {2: 11}
        base10 = {2: (0.08, 0)}
        base30 = {2: (0.20, 0)}
        fixed10 = {2: (0.02, 0)}
        fixed30 = {2: (0.10, 0)}
        rows = ledger.build_ledger(
            candidate_peel,
            {},
            base10,
            base30,
            fixed10,
            fixed30,
            baseline_peel_fixups=baseline_peel,
            baseline_dense_fixups={},
        )
        self.assertEqual(rows[0]["Decision"], "retune-exact-fixup")
        self.assertEqual(rows[0]["PeelFixup"], 11)
        self.assertEqual(rows[0]["BaselinePeelFixup"], 10)

    def test_mean_gates_are_inclusive_and_require_zero_failures(self):
        arguments = ({2: 10}, {}, {2: (0.05, 0)}, {2: (0.15, 0)})
        rows = ledger.build_ledger(*arguments)
        self.assertEqual(rows[0]["Decision"], "remove-exact-fixup")
        rows = ledger.build_ledger(
            {2: 10}, {}, {2: (0.05, 1)}, {2: (0.15, 0)}
        )
        self.assertEqual(rows[0]["Decision"], "retain-pending-enabled-holdout")
        rows = ledger.build_ledger(
            {2: 10}, {}, {2: (0.05004, 0)}, {2: (0.15, 0)}
        )
        self.assertEqual(rows[0]["Decision"], "retain-pending-enabled-holdout")

    def test_coverage_and_metadata_are_strict(self):
        peel = self.fixups("peel.inc", [(2, 10)])
        dense = self.fixups("dense.inc", [(3, 20)])
        path = self.ohead("short", 0.10, [(2, 0.01, 0)], trials=9999)
        with self.assertRaisesRegex(ledger.LedgerError, "trials"):
            ledger.read_ohead(path, 0.10, 10000)
        with self.assertRaisesRegex(ledger.LedgerError, "coverage mismatch"):
            ledger.build_ledger(
                peel, dense, {2: (0.0, 0)}, {2: (0.0, 0)}
            )

        missing_policy = self.ohead(
            "missing-policy", 0.10, [(2, 0.01, 0)], policy=None
        )
        with self.assertRaisesRegex(ledger.LedgerError, "seed_fixups"):
            ledger.read_ohead(missing_policy, 0.10, 10000)

        missing_seed = self.ohead(
            "missing-seed", 0.10, [(2, 0.01, 0)], seed=None
        )
        with self.assertRaisesRegex(ledger.LedgerError, "exactly one"):
            ledger.read_ohead(missing_seed, 0.10, 10000)

        huge_seed = self.ohead(
            "huge-seed", 0.10, [(2, 0.01, 0)], seed=str(1 << 64)
        )
        with self.assertRaisesRegex(ledger.LedgerError, "exceeds uint64"):
            ledger.read_ohead(huge_seed, 0.10, 10000)

        duplicate_policy = self.ohead(
            "duplicate-policy", 0.10, [(2, 0.01, 0)],
            policy="disabled seed_fixups=enabled",
        )
        with self.assertRaisesRegex(ledger.LedgerError, "exactly one"):
            ledger.read_ohead(duplicate_policy, 0.10, 10000)

        bad_count = self.ohead(
            "bad-count", 0.10, [(2, 0.01, 0)], metadata_count=2
        )
        with self.assertRaisesRegex(ledger.LedgerError, "Ns=2"):
            ledger.read_ohead(bad_count, 0.10, 10000)

        excessive_failures = self.ohead(
            "excessive-failures", 0.10, [(2, 0.01, 10001)]
        )
        with self.assertRaisesRegex(ledger.LedgerError, "exceeds trials"):
            ledger.read_ohead(excessive_failures, 0.10, 10000)

        repair_only = self.ohead(
            "repair-only", 0.10, [(2, 0.01, 0)], start_mode=1
        )
        with self.assertRaisesRegex(ledger.LedgerError, "startMode=0"):
            ledger.read_ohead(repair_only, 0.10, 10000)
        matrix_only = self.ohead(
            "matrix-only", 0.10, [(2, 0.01, 0)], bb=0
        )
        with self.assertRaisesRegex(ledger.LedgerError, "real-byte bb=64"):
            ledger.read_ohead(matrix_only, 0.10, 10000)

    def test_holdout_seed_pairing_and_training_disjointness(self):
        ledger.validate_holdout_seeds(1001, 2002, 1001, 2002, {1, 2})
        with self.assertRaisesRegex(ledger.LedgerError, "pair"):
            ledger.validate_holdout_seeds(1001, 2002, 9999, 2002, set())
        with self.assertRaisesRegex(ledger.LedgerError, "distinct"):
            ledger.validate_holdout_seeds(1001, 1001, 1001, 1001, set())
        with self.assertRaisesRegex(ledger.LedgerError, "overlap"):
            ledger.validate_holdout_seeds(1001, 2002, 1001, 2002, {2002})

        manifest = self.root / "manifest.txt"
        manifest.write_text("seed=77\nseeds=1,2,3\n", encoding="ascii")
        self.assertEqual(ledger.read_training_seeds([manifest]), {1, 2, 3, 77})

    def test_empty_generated_fixup_include_round_trips(self):
        try:
            from . import apply_fixup_holdout_ledger as apply
        except ImportError:
            import apply_fixup_holdout_ledger as apply
        path = self.root / "empty.inc"
        path.write_text(apply.render_fixups({}, "peel"), encoding="ascii")
        self.assertEqual(ledger.read_fixups(path), {})
        path.write_text(
            "// retired { 2, 17 },\n    { 3, 19 }, // old { 4, 20 },\n",
            encoding="ascii",
        )
        self.assertEqual(ledger.read_fixups(path), {3: 19})

    def test_atomic_ledger_contains_gates_trials_and_fixed_metrics(self):
        rows = ledger.build_ledger(
            {2: 10}, {}, {2: (0.05, 0)}, {2: (0.15, 0)},
            {2: (0.04, 0)}, {2: (0.14, 0)},
            base_trials10=10001,
            base_trials30=10002,
            fixed_trials10=10003,
            fixed_trials30=10004,
            minimum_trials=10000,
        )
        path = self.root / "ledger.tsv"
        ledger.write_atomic(path, rows)
        fields = path.read_text(encoding="ascii").splitlines()
        self.assertEqual(tuple(fields[0].split("\t")), ledger.LEDGER_HEADER)
        record = dict(zip(ledger.LEDGER_HEADER, fields[1].split("\t")))
        self.assertEqual(record["FixedMean10"], "0.040000000000000001")
        self.assertEqual(record["FixedFail30"], "0")
        self.assertEqual(record["BaseTrials10"], "10001")
        self.assertEqual(record["MeanGate10"], "0.050000000000000003")
        with self.assertRaisesRegex(ledger.LedgerError, "already exists"):
            ledger.write_atomic(path, rows)


if __name__ == "__main__":
    unittest.main()
