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

    def shard(self, name, body, seed="1", trials="500", maximum="10"):
        directory = self.root / name
        directory.mkdir()
        path = directory / "dense_count.out"
        path.write_text(body, encoding="utf-8")
        (directory / "manifest.txt").write_text(
            "dcount_seed=%s\ndcount_trials=%s\ndcount_max_failures=%s\n"
            "dcount_low_count_run=4\n" % (seed, trials, maximum),
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
