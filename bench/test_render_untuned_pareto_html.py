#!/usr/bin/env python3
"""Focused tests for the historical untuned-Pareto renderer inputs."""

from __future__ import annotations

import copy
import json
import pathlib
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import render_untuned_pareto_html as renderer  # noqa: E402


class HistoricalDiagnosticTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = pathlib.Path(self._tmp.name)

    @staticmethod
    def _cells(start_k: int, outcome: int) -> list[list[int]]:
        # Match the historical seed-major layout while deliberately exposing no
        # seed ID, just like cells_detail in the sealed result JSON.
        return [[k, outcome] for _seed in range(10) for k in range(start_k, 101)]

    @staticmethod
    def _write_doc(root: pathlib.Path, method: str, cells: list[list[int]]) -> pathlib.Path:
        result_dir = root / "results"
        result_dir.mkdir(parents=True, exist_ok=True)
        path = result_dir / f"{method}.json"
        path.write_text(json.dumps({
            "method": method,
            "regions": {
                renderer.HISTORICAL_DIAGNOSTIC_REGION: {
                    "cells_detail": cells,
                },
            },
        }))
        return path

    def _valid_roots(self, tag: str = "valid") -> tuple[pathlib.Path, pathlib.Path]:
        v1 = self.root / tag / "v1"
        v2 = self.root / tag / "v2"
        # Both arms contain K=2..5.  An observed-K intersection would include
        # those cells; the fixed historical diagnostic must not.
        alpha = self._cells(2, 0)
        alpha[4] = [6, -1]
        self._write_doc(v1, "alpha", alpha)
        self._write_doc(v2, "beta", self._cells(2, 1))
        return v1, v2

    @staticmethod
    def _fixture_sources(v1: pathlib.Path, v2: pathlib.Path) -> dict:
        template = {
            "schema": renderer.RESULT_SOURCE_SCHEMA,
            "metadata_status": "synthetic unit-test fixture",
            "round1": {
                "git_tree_oid": "1" * 40,
                "result_file_count": 0,
                "result_file_manifest_sha256": "0" * 64,
            },
            "round2": {
                "git_tree_oid": "2" * 40,
                "result_file_count": 0,
                "result_file_manifest_sha256": "0" * 64,
            },
        }
        return renderer._result_source_identities(v1, v2, template)

    def _build(self, tag: str = "valid") -> dict:
        v1, v2 = self._valid_roots(tag)
        return renderer.historical_diagnostic_overhead(
            v1, v2, renderer.HISTORICAL_DIAGNOSTIC_REGION,
            expected_sources=self._fixture_sources(v1, v2),
        )

    def test_fixed_grid_and_classification(self) -> None:
        cache = self._build()
        self.assertEqual(cache["schema"], renderer.HISTORICAL_DIAGNOSTIC_SCHEMA)
        self.assertEqual(
            cache["classification"],
            {
                "historical": True,
                "diagnostic_only": True,
                "selection_eligible": False,
                "architecture_evidence_eligible": False,
                "seed_ids_present": False,
            },
        )
        self.assertEqual(cache["grid"], {
            "region": "r1_2_100",
            "k_min": 6,
            "k_max": 100,
            "k_count": 95,
            "cells_per_k": 10,
            "cells_per_arm": 950,
        })
        self.assertIn("never architecture evidence", cache["notice"])
        self.assertEqual(cache["sources"]["round1"]["result_file_count"], 1)
        self.assertEqual(cache["sources"]["round2"]["result_file_count"], 1)
        self.assertEqual(cache["arms_sha256"], renderer._canonical_json_sha256(cache["arms"]))
        self.assertEqual(cache["arms"]["alpha"]["nc"], 950)
        self.assertAlmostEqual(cache["arms"]["alpha"]["eohc"], 4 / 950)
        self.assertEqual(cache["arms"]["beta"], {"eohc": 1.0, "nc": 950})
        self.assertEqual(
            renderer.validate_historical_diagnostic_cache(
                cache, {"alpha", "beta"}, expected_sources=cache["sources"],
            ),
            cache["arms"],
        )

    def test_same_shape_source_substitution_is_rejected(self) -> None:
        v1, v2 = self._valid_roots("source-substitution")
        authenticated = self._fixture_sources(v1, v2)
        path = v1 / "results" / "alpha.json"
        doc = json.loads(path.read_text())
        detail = doc["regions"][renderer.HISTORICAL_DIAGNOSTIC_REGION]["cells_detail"]
        detail[50][1] = 1
        path.write_text(json.dumps(doc))
        with self.assertRaisesRegex(SystemExit, "result sources.*manifest_sha256"):
            renderer.historical_diagnostic_overhead(
                v1, v2, renderer.HISTORICAL_DIAGNOSTIC_REGION,
                expected_sources=authenticated,
            )

    def test_authenticated_bytes_are_the_bytes_scored(self) -> None:
        v1, v2 = self._valid_roots("snapshot-toctou")
        authenticated = self._fixture_sources(v1, v2)
        original = renderer._authenticated_result_documents

        def mutate_after_snapshot(*args: object, **kwargs: object) -> tuple[dict, list]:
            receipt, documents = original(*args, **kwargs)
            path = v1 / "results" / "alpha.json"
            doc = json.loads(path.read_text())
            detail = doc["regions"][renderer.HISTORICAL_DIAGNOSTIC_REGION]["cells_detail"]
            for entry in detail:
                if entry[0] >= renderer.HISTORICAL_DIAGNOSTIC_K_MIN:
                    entry[1] = 4
            path.write_text(json.dumps(doc))
            return receipt, documents

        with mock.patch.object(
                renderer, "_authenticated_result_documents",
                side_effect=mutate_after_snapshot):
            cache = renderer.historical_diagnostic_overhead(
                v1, v2, renderer.HISTORICAL_DIAGNOSTIC_REGION,
                expected_sources=authenticated,
            )
        self.assertAlmostEqual(cache["arms"]["alpha"]["eohc"], 4 / 950)

    def test_source_grid_mutations_are_rejected(self) -> None:
        mutations = {
            "missing K": (
                lambda cells: [entry for entry in cells if entry[0] != 42],
                r"missing K=42",
            ),
            "wrong per-K multiplicity": (
                lambda cells: self._remove_one(cells, 42),
                r"K=42 has multiplicity 9",
            ),
            "unexpected K": (
                lambda cells: cells + [[101, 0]],
                r"unexpected K=101",
            ),
        }
        for index, (name, (mutate, message)) in enumerate(mutations.items()):
            with self.subTest(name=name):
                v1, v2 = self._valid_roots(f"mutation-{index}")
                path = v1 / "results" / "alpha.json"
                doc = json.loads(path.read_text())
                detail = doc["regions"][renderer.HISTORICAL_DIAGNOSTIC_REGION]["cells_detail"]
                doc["regions"][renderer.HISTORICAL_DIAGNOSTIC_REGION]["cells_detail"] = mutate(detail)
                path.write_text(json.dumps(doc))
                with self.assertRaisesRegex(SystemExit, message):
                    renderer.historical_diagnostic_overhead(
                        v1, v2, renderer.HISTORICAL_DIAGNOSTIC_REGION,
                        expected_sources=self._fixture_sources(v1, v2),
                    )

    @staticmethod
    def _remove_one(cells: list[list[int]], wanted_k: int) -> list[list[int]]:
        result = list(cells)
        del result[next(i for i, entry in enumerate(result) if entry[0] == wanted_k)]
        return result

    def test_malformed_cache_schema_is_rejected(self) -> None:
        cache = self._build()
        mutations = {
            "legacy unversioned cache": (copy.deepcopy(cache["arms"]), "fields disagree"),
        }

        bad = copy.deepcopy(cache)
        bad["schema"] = "wirehair.untuned-pareto.r1-historical-diagnostic.v0"
        mutations["wrong schema"] = (bad, "schema mismatch")

        bad = copy.deepcopy(cache)
        bad["classification"]["selection_eligible"] = True
        mutations["selectable"] = (bad, "classification.selection_eligible mismatch")

        bad = copy.deepcopy(cache)
        bad["classification"]["historical"] = 1
        mutations["numeric boolean"] = (bad, "classification.historical type mismatch")

        bad = copy.deepcopy(cache)
        bad["grid"]["k_min"] = 5
        mutations["wrong grid"] = (bad, "grid.k_min mismatch")

        bad = copy.deepcopy(cache)
        bad["grid"]["k_min"] = 6.0
        mutations["floating grid count"] = (bad, "grid.k_min type mismatch")

        bad = copy.deepcopy(cache)
        bad["sources"]["round1"]["result_file_manifest_sha256"] = "0" * 64
        mutations["wrong source identity"] = (bad, "cache sources.round1.result_file_manifest_sha256 mismatch")

        bad = copy.deepcopy(cache)
        bad["unexpected"] = True
        mutations["extra top-level field"] = (bad, "fields disagree")

        bad = copy.deepcopy(cache)
        del bad["arms"]["alpha"]["nc"]
        mutations["missing arm field"] = (bad, "fields disagree")

        bad = copy.deepcopy(cache)
        bad["arms"]["alpha"]["nc"] = 949
        mutations["wrong cell count"] = (bad, "must have exactly 950 cells")

        bad = copy.deepcopy(cache)
        bad["arms"]["alpha"]["eohc"] = float("nan")
        mutations["non-finite overhead"] = (bad, "invalid eohc")

        for name, (mutated, message) in mutations.items():
            with self.subTest(name=name):
                with self.assertRaisesRegex(SystemExit, message):
                    renderer.validate_historical_diagnostic_cache(
                        mutated, expected_sources=cache["sources"],
                    )

        with self.assertRaisesRegex(SystemExit, "method sets disagree"):
            renderer.validate_historical_diagnostic_cache(
                cache, {"alpha"}, expected_sources=cache["sources"],
            )

    def test_checked_in_cache_uses_strict_schema(self) -> None:
        repo = pathlib.Path(__file__).resolve().parent.parent
        cache_path = repo / "bench" / "results" / "pareto" / "r1_common.json"
        cache = json.loads(cache_path.read_text())
        arms = renderer.validate_historical_diagnostic_cache(cache)
        self.assertEqual(len(arms), 84)
        self.assertTrue(all(arm["nc"] == 950 for arm in arms.values()))

        page = (repo / "results" / "untuned_pareto_screen.html").read_text()
        self.assertIn("r1: historical diagnostic", page)
        self.assertIn("non-selectable and must never be used as architecture evidence", page)
        self.assertIn("Failure, dead, hard-K and construction metrics hidden", page)
        self.assertNotIn("r1: like-for-like", page)

        rows_cache = json.loads((repo / "bench" / "results" / "pareto" / "rows.json").read_text())
        rows = renderer.validate_rows_cache(rows_cache)
        self.assertEqual(len(rows), 504)


if __name__ == "__main__":
    unittest.main()
