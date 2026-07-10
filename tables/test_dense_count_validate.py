#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

try:
    from . import dense_count_validate as validate
except ImportError:
    import dense_count_validate as validate


class DenseCountValidateTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def shard(
        self, name, body, seed="1", trials="500", maximum="10",
        low_run="4", selection=None, manifest=True,
    ):
        directory = self.root / name
        directory.mkdir()
        path = directory / "dense_count.out"
        path.write_text(body, encoding="utf-8")
        if manifest:
            selection_text = (
                "" if selection is None else f"dcount_selection={selection}\n"
            )
            (directory / "manifest.txt").write_text(
                "dcount_seed=%s\ndcount_trials=%s\ndcount_max_failures=%s\n"
                "dcount_low_count_run=%s\n%s" % (
                    seed, trials, maximum, low_run, selection_text,
                ),
                encoding="utf-8",
            )
        return path

    def test_compatible_distinct_seeds_and_disjoint_same_seed(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", seed="11")
        b = self.shard("b", header + "2\t2\t0.05\n", seed="12")
        c = self.shard("c", header + "3\t3\t0.0\n", seed="11")
        rows = validate.validate_shards([str(a), str(b), str(c)])
        self.assertEqual([row[0] for row in rows], [2, 2, 3])

    def test_reused_seed_overlap_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", seed="11")
        b = self.shard("b", header + "2\t2\t0.05\n", seed="11")
        with self.assertRaisesRegex(validate.ValidationError, "reused seed"):
            validate.validate_shards([str(a), str(b)])

    def test_incompatible_method_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", trials="500")
        b = self.shard("b", header + "3\t2\t0.1\n", trials="1000")
        with self.assertRaisesRegex(validate.ValidationError, "incompatible"):
            validate.validate_shards([str(a), str(b)])

    def test_seed_and_methodology_metadata_are_canonical(self):
        header = "N\tDenseCount\tLowestFailures\n"
        cases = [
            ("01", "500", "10", "4", "unsigned integer"),
            ("garbage", "500", "10", "4", "unsigned integer"),
            ("1", "0", "0", "4", "outside"),
            ("1", "500", "-1", "4", "unsigned integer"),
            ("1", "500", "10", "0", "outside"),
        ]
        for index, (seed, trials, maximum, low_run, error) in enumerate(cases):
            with self.subTest(index=index):
                path = self.shard(
                    f"metadata-{index}", header + "2\t1\t0\n",
                    seed=seed, trials=trials, maximum=maximum,
                    low_run=low_run,
                )
                with self.assertRaisesRegex(validate.ValidationError, error):
                    validate.read_shard(str(path))

    def test_header_method_cannot_be_overridden_by_manifest(self):
        header = "N\tDenseCount\tLowestFailures\n"
        path = self.shard(
            "contradiction", header + "2\t1\t0\n",
            selection="threshold-run-v1",
        )
        with self.assertRaisesRegex(validate.ValidationError, "contradicts"):
            validate.read_shard(str(path))

    def test_multi_shard_requires_complete_provenance(self):
        header = "N\tDenseCount\tLowestFailures\n"
        proven = self.shard("proven", header + "2\t1\t0\n")
        bare = self.shard(
            "bare", header + "3\t1\t0\n", manifest=False
        )
        with self.assertRaisesRegex(validate.ValidationError, "incomplete"):
            validate.validate_shards([str(proven), str(bare)])

    def test_duplicate_manifest_keys_are_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        path = self.shard("duplicate-key", header + "2\t1\t0\n")
        manifest = path.parent / "manifest.txt"
        manifest.write_text(
            manifest.read_text(encoding="utf-8") + "dcount_seed=2\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(validate.ValidationError, "duplicate"):
            validate.read_shard(str(path))

    def test_shard_schema_and_domains(self):
        cases = [
            ("2\t1\t0.1\n", "expected tab-separated header"),
            ("N\tDenseCount\tLowestFailures\n1\t1\t0\n", "outside"),
            ("N\tDenseCount\tLowestFailures\n2\t401\t0\n", "outside"),
            ("N\tDenseCount\tLowestFailures\n2\t1\tNaN\n", "finite"),
            ("N\tDenseCount\tLowestFailures\n2\t1\t0\n2\t2\t0\n", "duplicate"),
        ]
        for index, (body, error) in enumerate(cases):
            with self.subTest(index=index):
                path = self.shard("bad%d" % index, body)
                with self.assertRaisesRegex(validate.ValidationError, error):
                    validate.read_shard(str(path))

    def test_aggregate_schema_and_domains(self):
        path = self.root / "aggregate.tsv"
        path.write_text(
            "\t".join(validate.AGGREGATE_HEADER) + "\n"
            "2\t2\t1\t2\t0.01\ta,b\n",
            encoding="utf-8",
        )
        self.assertEqual(validate.read_aggregate(str(path))[0][:5], (2, 2, 1, 2, "0.01"))
        path.write_text(
            "\t".join(validate.AGGREGATE_HEADER) + "\n"
            "2\t1\t3\t2\t0.01\ta\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(validate.ValidationError, "min dense"):
            validate.read_aggregate(str(path))


if __name__ == "__main__":
    unittest.main()
