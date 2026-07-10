#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

try:
    from . import apply_fixup_holdout_ledger as apply
    from . import fixup_holdout_ledger as ledger
except ImportError:
    import apply_fixup_holdout_ledger as apply
    import fixup_holdout_ledger as ledger


class ApplyFixupHoldoutLedgerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.baseline_peel = {2: 10, 3: 11, 4: 12}
        self.baseline_dense = {3: 20, 4: 21}
        self.candidate_peel = {2: 10, 3: 11, 4: 13}
        self.candidate_dense = dict(self.baseline_dense)
        self.base10 = {2: (0.05, 0), 3: (0.08, 0), 4: (0.09, 0)}
        self.base30 = {2: (0.15, 0), 3: (0.20, 0), 4: (0.21, 0)}
        self.fixed10 = {n: (0.04, 0) for n in (2, 3, 4)}
        self.fixed30 = {n: (0.14, 0) for n in (2, 3, 4)}

    def tearDown(self):
        self.temp.cleanup()

    def build_rows(self):
        return ledger.build_ledger(
            self.candidate_peel,
            self.candidate_dense,
            self.base10,
            self.base30,
            self.fixed10,
            self.fixed30,
            baseline_peel_fixups=self.baseline_peel,
            baseline_dense_fixups=self.baseline_dense,
            base_trials10=10000,
            base_trials30=10000,
            fixed_trials10=10000,
            fixed_trials30=10000,
            minimum_trials=10000,
            base_seed10=1001,
            base_seed30=2002,
            fixed_seed10=1001,
            fixed_seed30=2002,
            block_bytes=64,
            start_mode=0,
            training_seeds={1, 2, 3},
        )

    def fixup_file(self, name, values, kind):
        path = self.root / name
        path.write_text(apply.render_fixups(values, kind), encoding="ascii")
        return path

    def ledger_file(self, name="ledger.tsv", rows=None):
        path = self.root / name
        ledger.write_atomic(path, self.build_rows() if rows is None else rows)
        return path

    def test_remove_retain_and_retune_apply_as_one_pair(self):
        rows = apply.read_ledger(self.ledger_file())
        self.assertEqual(
            [row["Decision"] for row in rows],
            ["remove-exact-fixup", "retain-exact-fixup", "retune-exact-fixup"],
        )
        peel, dense = apply.apply_rows(
            rows, self.baseline_peel, self.baseline_dense
        )
        self.assertEqual(peel, {3: 11, 4: 13})
        self.assertEqual(dense, {3: 20, 4: 21})

        peel_output = self.root / "peel-output.inc"
        dense_output = self.root / "dense-output.inc"
        apply.write_pair_no_replace(
            peel_output,
            apply.render_fixups(peel, "peel"),
            dense_output,
            apply.render_fixups(dense, "dense"),
        )
        self.assertEqual(ledger.read_fixups(peel_output), peel)
        self.assertEqual(ledger.read_fixups(dense_output), dense)

    def test_main_e2e_and_no_clobber(self):
        ledger_path = self.ledger_file()
        peel_input = self.fixup_file("peel.inc", self.baseline_peel, "peel")
        dense_input = self.fixup_file("dense.inc", self.baseline_dense, "dense")
        peel_output = self.root / "new-peel.inc"
        dense_output = self.root / "new-dense.inc"
        arguments = [
            "--ledger", str(ledger_path),
            "--peel-input", str(peel_input),
            "--dense-input", str(dense_input),
            "--peel-output", str(peel_output),
            "--dense-output", str(dense_output),
        ]
        self.assertEqual(apply.main(arguments), 0)
        original_peel = peel_output.read_bytes()
        original_dense = dense_output.read_bytes()
        self.assertEqual(apply.main(arguments), 1)
        self.assertEqual(peel_output.read_bytes(), original_peel)
        self.assertEqual(dense_output.read_bytes(), original_dense)

    def test_unresolved_or_weak_evidence_blocks_all_output(self):
        rows = self.build_rows()
        rows[2]["Decision"] = "retune-required"
        unresolved = self.ledger_file("unresolved.tsv", rows)
        with self.assertRaisesRegex(apply.ApplyError, "unresolved decision"):
            apply.read_ledger(unresolved)

        rows = self.build_rows()
        rows[0]["MinimumTrials"] = 9999
        rows[0]["BaseTrials10"] = 9999
        weak = self.ledger_file("weak.tsv", rows)
        with self.assertRaisesRegex(apply.ApplyError, "too weak"):
            apply.read_ledger(weak)

        rows = self.build_rows()
        rows[0]["BlockBytes"] = 1
        wrong_bytes = self.ledger_file("wrong-bytes.tsv", rows)
        with self.assertRaisesRegex(apply.ApplyError, "bb=64 startMode=0"):
            apply.read_ledger(wrong_bytes)

    def test_baseline_mismatch_and_partial_publication_are_rejected(self):
        rows = apply.read_ledger(self.ledger_file())
        with self.assertRaisesRegex(apply.ApplyError, "peel input"):
            apply.apply_rows(rows, {2: 99, 3: 11, 4: 12}, self.baseline_dense)

        dense_output = self.root / "existing-dense.inc"
        dense_output.write_text("sentinel\n", encoding="ascii")
        peel_output = self.root / "absent-peel.inc"
        with self.assertRaisesRegex(apply.ApplyError, "already exists"):
            apply.write_pair_no_replace(
                peel_output, "peel\n", dense_output, "dense\n"
            )
        self.assertFalse(peel_output.exists())
        self.assertEqual(dense_output.read_text(encoding="ascii"), "sentinel\n")


if __name__ == "__main__":
    unittest.main()
