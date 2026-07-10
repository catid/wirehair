#!/usr/bin/env python3

import os
import shutil
import subprocess
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
            nlo="2", nhi="2", manifest=True, omit=(), extra=""):
        directory = self.root / name
        directory.mkdir()
        path = directory / "dense_count.out"
        path.write_text(body, encoding="utf-8")
        if manifest:
            values = {
                "dcount_schema": str(validate.DCOUNT_SCHEMA_VERSION),
                "dcount_methodology": validate.DCOUNT_METHODOLOGY,
                "dcount_seed": str(seed),
                "dcount_trials": str(trials),
                "dcount_max_failures": str(maximum),
                "dcount_low_count_run": "4",
                "dcount_nlo": str(nlo),
                "dcount_nhi": str(nhi),
            }
            (directory / "manifest.txt").write_text(
                "".join(
                    "%s=%s\n" % (key, value)
                    for key, value in values.items() if key not in omit
                ) + extra,
                encoding="utf-8",
            )
        return path

    def aggregate(self, name, rows):
        path = self.root / name
        path.write_text(
            "\t".join(validate.AGGREGATE_HEADER) + "\n" + rows,
            encoding="utf-8",
        )
        return path

    def test_compatible_distinct_seeds_and_disjoint_same_seed(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", seed="11", nlo="2", nhi="2")
        b = self.shard("b", header + "2\t2\t0.05\n", seed="12", nlo="2", nhi="2")
        c = self.shard("c", header + "3\t3\t0.0\n", seed="11", nlo="3", nhi="3")
        rows = validate.validate_shards([str(a), str(b), str(c)])
        self.assertEqual([row[0] for row in rows], [2, 2, 3])

    def test_reused_seed_overlap_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", seed="11", nlo="2", nhi="3")
        b = self.shard("b", header + "3\t2\t0.05\n", seed="11", nlo="3", nhi="4")
        with self.assertRaisesRegex(validate.ValidationError, "reused seed"):
            validate.validate_shards([str(a), str(b)])

    def test_incompatible_method_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("a", header + "2\t1\t0.1\n", trials="500", nlo="2", nhi="2")
        b = self.shard("b", header + "3\t2\t0.1\n", trials="1000", nlo="3", nhi="3")
        with self.assertRaisesRegex(validate.ValidationError, "incompatible"):
            validate.validate_shards([str(a), str(b)])

    def test_multi_shard_requires_every_provenance_field(self):
        header = "N\tDenseCount\tLowestFailures\n"
        required = {
            "schema": "dcount_schema",
            "methodology": "dcount_methodology",
            "seed": "dcount_seed",
            "trials": "dcount_trials",
            "max_failures": "dcount_max_failures",
            "low_count_run": "dcount_low_count_run",
            "nlo": "dcount_nlo",
            "nhi": "dcount_nhi",
        }
        for index, (field, manifest_key) in enumerate(required.items()):
            with self.subTest(field=field):
                a = self.shard(
                    "missing%d" % index, header + "2\t1\t0.1\n",
                    omit=(manifest_key,),
                )
                b = self.shard(
                    "complete%d" % index, header + "3\t1\t0.1\n",
                    seed="2", nlo="3", nhi="3",
                )
                with self.assertRaisesRegex(
                        validate.ValidationError, "missing .*%s" % field):
                    validate.validate_shards([str(a), str(b)])

        bare_a = self.shard("bare_a", header + "2\t1\t0.1\n", manifest=False)
        bare_b = self.shard("bare_b", header + "2\t1\t0.1\n", manifest=False)
        with self.assertRaisesRegex(validate.ValidationError, "complete provenance"):
            validate.validate_shards([str(bare_a), str(bare_b)])

        legacy_seed = self.shard(
            "legacy_seed", header + "2\t1\t0.1\n",
            omit=("dcount_seed",), extra="seed=1\n",
        )
        complete = self.shard(
            "explicit_seed", header + "3\t1\t0.1\n",
            seed="2", nlo="3", nhi="3",
        )
        with self.assertRaisesRegex(validate.ValidationError, "missing seed"):
            validate.validate_shards([str(legacy_seed), str(complete)])

        # A single historical shard remains inspectable without a manifest;
        # provenance becomes mandatory only when evidence is combined.
        self.assertEqual(validate.validate_shards([str(bare_a)])[0][0], 2)

    def test_manifest_types_versions_duplicates_and_ranges(self):
        header = "N\tDenseCount\tLowestFailures\n"
        cases = [
            ({"extra": "dcount_seed=2\n"}, "duplicate manifest key"),
            ({"seed": "-1"}, "unsigned integer"),
            ({"trials": "ten"}, "unsigned integer"),
            ({"extra": "malformed\n"}, "malformed manifest entry"),
            ({"omit": ("dcount_schema",), "extra": "dcount_schema=2\n"}, "unsupported"),
            ({"omit": ("dcount_methodology",),
              "extra": "dcount_methodology=unknown\n"}, "unsupported"),
            ({"nlo": "3", "nhi": "2"}, "range .* reversed"),
            ({"nlo": "3", "nhi": "4"}, "outside declared"),
        ]
        for index, (kwargs, error) in enumerate(cases):
            with self.subTest(index=index):
                path = self.shard(
                    "manifest_bad%d" % index,
                    header + "2\t1\t0.1\n",
                    **kwargs,
                )
                with self.assertRaisesRegex(validate.ValidationError, error):
                    validate.validate_shards([str(path)])

        checksummed = self.shard(
            "checksummed", header + "2\t1\t0.1\n",
            extra="0" * 64 + "  /tmp/a=b/dense_count.out\n",
        )
        self.assertEqual(validate.validate_shards([str(checksummed)])[0][0], 2)

    def test_repeated_input_identity_rejected(self):
        header = "N\tDenseCount\tLowestFailures\n"
        path = self.shard("original", header + "2\t1\t0.1\n")
        alias = self.root / "alias.tsv"
        os.link(path, alias)
        with self.assertRaisesRegex(validate.ValidationError, "repeated dense-count shard"):
            validate.validate_shards([str(path), str(alias)])

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
        path = self.aggregate("aggregate.tsv", "2\t2\t1\t2\t0.01\ta,b\n")
        self.assertEqual(validate.read_aggregate(str(path))[0][:5], (2, 2, 1, 2, "0.01"))
        path.write_text(
            "\t".join(validate.AGGREGATE_HEADER) + "\n"
            "2\t1\t3\t2\t0.01\ta\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(validate.ValidationError, "min dense"):
            validate.read_aggregate(str(path))

    def test_multiple_aggregates_require_unique_paths_and_n(self):
        a = self.aggregate("a.tsv", "2\t2\t1\t2\t0.01\ta,b\n")
        b = self.aggregate("b.tsv", "3\t1\t2\t2\t0.02\tc\n")
        c = self.aggregate("c.tsv", "2\t1\t2\t2\t0.02\td\n")
        self.assertEqual(
            [row[0] for row in validate.validate_aggregates([str(a), str(b)])],
            [2, 3],
        )
        with self.assertRaisesRegex(validate.ValidationError, "repeated dense-count aggregate"):
            validate.validate_aggregates([str(a), str(a)])
        with self.assertRaisesRegex(validate.ValidationError, "duplicate aggregate N=2"):
            validate.validate_aggregates([str(a), str(c)])

    def test_pipeline_callers_enforce_combination_rules(self):
        root = Path(__file__).resolve().parents[1]
        if os.name == "nt":
            self.skipTest("dense-count pipeline scripts are POSIX-only")
        bash = shutil.which("bash")
        cxx = shutil.which("g++") or shutil.which("clang++")
        if bash is None or cxx is None:
            self.skipTest("pipeline caller test requires bash and a GNU-compatible C++ compiler")
        environment = os.environ.copy()
        environment["CXX"] = cxx
        header = "N\tDenseCount\tLowestFailures\n"
        a = self.shard("pipeline_a", header + "2\t1\t0.1\n", seed="11")
        b = self.shard("pipeline_b", header + "2\t2\t0.05\n", seed="12")

        aggregate = subprocess.run(
            [bash, str(root / "tables/aggregate_dense_count_shards.sh"), str(a), str(b)],
            cwd=root, env=environment, text=True, capture_output=True, check=False,
        )
        self.assertEqual(aggregate.returncode, 0, aggregate.stderr)
        self.assertIn("2\t2\t1\t2", aggregate.stdout)
        aggregate_path = self.root / "combined.tsv"
        aggregate_path.write_text(aggregate.stdout, encoding="utf-8")

        derive = subprocess.run(
            [bash, str(root / "tables/derive_dense_count_candidate.sh"), str(aggregate_path)],
            cwd=root, env=environment, text=True, capture_output=True, check=False,
        )
        self.assertEqual(derive.returncode, 0, derive.stderr)
        self.assertIn("2\t2\t0.05\t2\t1", derive.stdout)

        for script, inputs in (
                ("aggregate_dense_count_shards.sh", [a, a]),
                ("compare_dense_count_shards.sh", [a, a]),
                ("derive_dense_count_candidate.sh", [aggregate_path, aggregate_path])):
            with self.subTest(script=script):
                result = subprocess.run(
                    [bash, str(root / "tables" / script)] +
                    [str(path) for path in inputs],
                    cwd=root, env=environment, text=True,
                    capture_output=True, check=False,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")
                self.assertIn("repeated dense-count", result.stderr)


if __name__ == "__main__":
    unittest.main()
