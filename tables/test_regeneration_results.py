#!/usr/bin/env python3

import csv
import tempfile
import unittest
from pathlib import Path

try:
    from . import regeneration_results as results
except ImportError:
    import regeneration_results as results


class RegenerationResultsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, kind, rows):
        path = self.root / name
        with path.open("w", encoding="ascii", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=results.SCHEMAS[kind],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def small_row(self, n, trials=10, base=123):
        return {
            "N": str(n),
            "DenseCount": "2",
            "DenseSeed": "3",
            "PeelSeed": "4",
            "DenseFailures": "0",
            "PeelScore": "1",
            "CurrentDenseCount": "2",
            "CurrentDenseSeed": "5",
            "CurrentPeelSeed": "6",
            "Trials": str(trials),
            "BaseSeed": str(base),
            "TrialSeed": str(results.derive_trial_seed(base, n, "small")),
            "SmallTablesHash": "234",
            "MethodVersion": "small-v2",
        }

    def dense_count_row(self, n, trials=10, base=123):
        return {
            "N": str(n),
            "DenseCount": "2",
            "SelectedFailures": "1",
            "Trials": str(trials),
            "MaxFailures": "2",
            "LowCountRun": "1",
            "BaseSeed": str(base),
            "TrialSeed": str(
                results.derive_trial_seed(base, n, "dense-count")
            ),
            "QualifyingCount": "1",
            "SourceTablesHash": "1234",
            "MethodVersion": "dense-count-v2",
            "Status": "selected-knee",
        }

    def dense_row(self, index, used=1, trials=10, base=123):
        return {
            "DenseIndex": str(index),
            "Used": str(used),
            "N": "2048" if used else "0",
            "DenseCount": str(index * 4 + 2),
            "DenseSeed": "7",
            "Failures": "1",
            "CurrentDenseSeed": "8",
            "Trials": str(trials),
            "BaseSeed": str(base),
            "TrialSeed": str(results.derive_trial_seed(base, index, "dense-seeds")),
            "SourceTablesHash": "456",
            "MethodVersion": "most-dense-v2",
        }

    def peel_row(self, subdivision, trials=10, base=123):
        first = 2048 + subdivision
        return {
            "Subdivision": str(subdivision),
            "FirstN": str(first),
            "LastN": str(first + 2048),
            "PeelSeed": "9",
            "Score": "2",
            "CurrentPeelSeed": "10",
            "Trials": str(trials),
            "MaxTries": "8",
            "EvaluatedCandidates": "3",
            "BaseSeed": str(base),
            "RequestedNLo": "2048",
            "RequestedNHi": "64000",
            "BaseTablesHash": "789",
            "MethodVersion": "peel-v2",
            "Status": "searched",
        }

    def test_small_shards_merge_like_monolithic(self):
        rows = [self.small_row(2), self.small_row(3)]
        monolithic = self.write("mono.tsv", "small", rows)
        shard_a = self.write("a.tsv", "small", rows[:1])
        shard_b = self.write("b.tsv", "small", rows[1:])
        self.assertEqual(
            results.validate_and_merge("small", [str(monolithic)]),
            results.validate_and_merge("small", [str(shard_b), str(shard_a)]),
        )

    def test_dense_count_accepts_independent_runs(self):
        a = self.write(
            "dense-count-a.tsv", "dense-count", [self.dense_count_row(2, base=11)]
        )
        b = self.write(
            "dense-count-b.tsv", "dense-count", [self.dense_count_row(2, base=12)]
        )
        merged = results.validate_and_merge("dense-count", [str(a), str(b)])
        self.assertEqual([row["BaseSeed"] for row in merged], ["11", "12"])

    def test_dense_count_fallback_and_failure_bounds(self):
        row = self.dense_count_row(2048)
        row["Status"] = "fallback-minimum"
        row["QualifyingCount"] = "0"
        path = self.write("dense-count-fallback.tsv", "dense-count", [row])
        self.assertEqual(
            results.validate_and_merge("dense-count", [str(path)])[0]["Status"],
            "fallback-minimum",
        )
        row["QualifyingCount"] = "1"
        path = self.write("dense-count-inconsistent.tsv", "dense-count", [row])
        with self.assertRaisesRegex(results.ValidationError, "inconsistent fallback"):
            results.validate_and_merge("dense-count", [str(path)])
        row = self.dense_count_row(2)
        row["SelectedFailures"] = "11"
        path = self.write("dense-count-bad.tsv", "dense-count", [row])
        with self.assertRaisesRegex(results.ValidationError, "exceeds Trials"):
            results.validate_and_merge("dense-count", [str(path)])

    def test_duplicate_and_method_mismatch_rejected(self):
        a = self.write("a.tsv", "small", [self.small_row(2)])
        b = self.write("b.tsv", "small", [self.small_row(2)])
        with self.assertRaisesRegex(results.ValidationError, "duplicate"):
            results.validate_and_merge("small", [str(a), str(b)])
        c = self.write("c.tsv", "small", [self.small_row(3, trials=11)])
        with self.assertRaisesRegex(results.ValidationError, "methodology"):
            results.validate_and_merge("small", [str(a), str(c)])

    def test_trial_seed_is_checked(self):
        row = self.small_row(2)
        row["TrialSeed"] = "0"
        path = self.write("bad-seed.tsv", "small", [row])
        with self.assertRaisesRegex(results.ValidationError, "TrialSeed"):
            results.validate_and_merge("small", [str(path)])

    def test_full_coverage_is_exact(self):
        path = self.write("partial.tsv", "small", [self.small_row(2)])
        with self.assertRaisesRegex(results.ValidationError, "incomplete"):
            results.validate_and_merge("small", [str(path)], require_full=True)

    def test_dense_seed_domain_and_usage(self):
        good = self.write("dense.tsv", "dense-seeds", [self.dense_row(13)])
        self.assertEqual(
            results.validate_and_merge("dense-seeds", [str(good)])[0]["N"],
            "2048",
        )
        row = self.dense_row(13, used=0)
        row["N"] = "2048"
        bad = self.write("dense-bad.tsv", "dense-seeds", [row])
        with self.assertRaisesRegex(results.ValidationError, "unused"):
            results.validate_and_merge("dense-seeds", [str(bad)])

    def test_peel_lattice_and_status(self):
        good = self.write("peel.tsv", "peel", [self.peel_row(7)])
        self.assertEqual(
            results.validate_and_merge("peel", [str(good)])[0]["Subdivision"],
            "7",
        )
        row = self.peel_row(7)
        row["LastN"] = str(int(row["LastN"]) + 1)
        bad = self.write("peel-bad.tsv", "peel", [row])
        with self.assertRaisesRegex(results.ValidationError, "lattice"):
            results.validate_and_merge("peel", [str(bad)])

    def test_truncated_and_non_ascii_files_rejected(self):
        truncated = self.root / "truncated.tsv"
        truncated.write_text("\t".join(results.SCHEMAS["small"]) + "\n", encoding="ascii")
        with self.assertRaisesRegex(results.ValidationError, "no result rows"):
            results.validate_and_merge("small", [str(truncated)])
        non_ascii = self.root / "non-ascii.tsv"
        non_ascii.write_bytes(b"N\t\xff\n")
        with self.assertRaisesRegex(results.ValidationError, "malformed ASCII"):
            results.validate_and_merge("small", [str(non_ascii)])

    def test_atomic_output_refuses_overwrite(self):
        rows = [self.small_row(2)]
        output = self.root / "merged.tsv"
        results.write_rows("small", rows, str(output))
        first = output.read_bytes()
        with self.assertRaisesRegex(results.ValidationError, "already exists"):
            results.write_rows("small", rows, str(output))
        self.assertEqual(output.read_bytes(), first)


if __name__ == "__main__":
    unittest.main()
