#!/usr/bin/env python3
"""Unit tests for precode result integrity validation."""

import contextlib
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

try:
    from . import validate_results as validator
except ImportError:
    import validate_results as validator


HEADER = validator.BASE + validator.HEAVY + ["runaway_rate"]
NEW_HEADER = HEADER + validator.RUN_METADATA


def valid_row(**changes):
    row = {name: "0" for name in HEADER}
    row.update({
        "K": "1000",
        "scheme": "dense",
        "D": "50",
        "H": "6",
        "oh": "0",
        "trials": "10",
        "fail_rate": "0.2",
        "fail_rate_noheavy": "0.5",
        "def_max": "7",
        "def_pdf": "0:0.5|1:0.3|7:0.2",
    })
    row.update({name: str(value) for name, value in changes.items()})
    return row


def valid_metadata_row(**changes):
    row = valid_row()
    row.update({
        "packet_schedule_exhausted_rate": "0.0", "rowdist": "fixed44",
        "packet_schedule": "iid", "loss": "0.1",
        "identity_systematic": "1", "mix": "3", "base_seed": "12345",
        "paired": "1", "max_inact": "0", "max_row_seconds": "0",
        "requested_trials": "10", "threads": "8", "ge_replay": "0",
        "ge_replay_reverse": "0", "ge_pivot_window": "0",
    })
    row.update({name: str(value) for name, value in changes.items()})
    return row


class ValidateResultsTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write_csv(self, name="results.csv", rows=None, header=None):
        path = self.root / name
        header = HEADER if header is None else header
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, extrasaction="ignore")
            writer.writeheader()
            for row in rows if rows is not None else [valid_row()]:
                writer.writerow(row)
        return path

    def assert_invalid(self, row, message=None):
        path = self.write_csv(rows=[row])
        with self.assertRaises(validator.ValidationError) as caught:
            validator.validate_file(path)
        if message:
            self.assertIn(message, str(caught.exception))

    def run_main(self, *arguments):
        stderr = io.StringIO()
        stdout = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(stdout):
            result = validator.main(list(map(str, arguments)))
        return result, stdout.getvalue(), stderr.getvalue()

    def test_valid_current_schema(self):
        rows = validator.validate_file(self.write_csv())
        self.assertEqual(1, len(rows))
        self.assertEqual("dense", rows[0]["scheme"])

    def test_new_metadata_schema_and_exhaustion_sentinel(self):
        rows = validator.validate_file(self.write_csv(
            "metadata.csv", rows=[valid_metadata_row()], header=NEW_HEADER))
        self.assertEqual(0.0, rows[0]["packet_schedule_exhausted_rate"])

        exhausted = valid_metadata_row(
            fail_rate="1.0", fail_rate_noheavy="1.0", def_mu="0",
            def_max="0", def_pdf="999998:1.0",
            packet_schedule_exhausted_rate="1.0", heavy_divs_mu="0",
        )
        rows = validator.validate_file(self.write_csv(
            "exhausted.csv", rows=[exhausted], header=NEW_HEADER))
        self.assertEqual(1.0, rows[0]["packet_schedule_exhausted_rate"])

        mismatch = dict(exhausted)
        mismatch["packet_schedule_exhausted_rate"] = "0.0"
        with self.assertRaisesRegex(
                validator.ValidationError,
                "packet_schedule_exhausted_rate"):
            validator.validate_file(self.write_csv(
                "exhausted_mismatch.csv", rows=[mismatch],
                header=NEW_HEADER))

    def test_invalid_headers_and_field_counts(self):
        wrong_header = self.write_csv("header.csv", header=HEADER[:-1])
        with self.assertRaisesRegex(validator.ValidationError, "unsupported"):
            validator.validate_file(wrong_header)

        truncated = self.root / "truncated.csv"
        truncated.write_text(",".join(HEADER) + "\n1000,dense\n", encoding="utf-8")
        with self.assertRaisesRegex(validator.ValidationError, "field count"):
            validator.validate_file(truncated)

    def test_invalid_numeric_domains(self):
        cases = [
            (valid_row(trials="x"), "unsigned integer"),
            (valid_row(inact_mu="nan"), "finite range"),
            (valid_row(fail_rate="1.1"), "finite range"),
            (valid_row(inact_mu="-1"), "finite range"),
            (valid_row(K="64000", D="2000"), "model domain"),
            (valid_row(K="64000", oh="2000"), "model domain"),
        ]
        for row, message in cases:
            with self.subTest(message=message, row=row):
                self.assert_invalid(row, message)

    def test_invalid_pdf_and_rescored_rates(self):
        cases = [
            (valid_row(def_pdf="0:0.5|0:0.5"), "duplicate"),
            (valid_row(def_pdf="0:0.8|7:0.1"), "mass"),
            (valid_row(def_pdf="bad"), "malformed"),
            (valid_row(def_pdf="0:0.5|1:nan|7:0.5"), "finite range"),
            (valid_row(def_pdf="0:0.5|1:0.3|70000:0.2", def_max="70000"), "model domain"),
            (valid_row(def_max="8"), "def_max"),
            (valid_row(fail_rate="0.1"), "rescore"),
            (valid_row(fail_rate_noheavy="0.4"), "rescore"),
            (valid_row(
                def_pdf="0:0.5|1:0.3|999999:0.2", def_max="1", runaway_rate="0.1"
            ), "runaway_rate"),
        ]
        for row, message in cases:
            with self.subTest(message=message, row=row):
                self.assert_invalid(row, message)

    def test_duplicate_rows_within_file(self):
        result_file = self.write_csv(rows=[valid_row(), valid_row()])
        with self.assertRaisesRegex(
                validator.ValidationError,
                r"results\.csv:3: duplicate experiment row .*first seen at .*results\.csv:2"):
            validator.validate_file(result_file)

    def test_repeated_path_reports_both_source_locations(self):
        result_file = self.write_csv()
        for manifest in (None, self.write_manifest()):
            with self.subTest(manifest=manifest):
                arguments = [result_file, result_file]
                if manifest is not None:
                    arguments.extend(("--manifest", manifest))
                result, _, error = self.run_main(*arguments)
                self.assertEqual(1, result)
                self.assertIn("duplicate experiment row", error)
                self.assertEqual(2, error.count(f"{result_file}:2"))

    def test_duplicate_cells_across_distinct_files(self):
        first = self.write_csv("first.csv")
        cases = (
            ("identical", valid_row()),
            ("conflicting", valid_row(trials=11)),
        )
        for kind, duplicate in cases:
            for manifest in (None, self.write_manifest()):
                with self.subTest(kind=kind, manifest=manifest):
                    second = self.write_csv("second.csv", rows=[duplicate])
                    arguments = [first, second]
                    if manifest is not None:
                        arguments.extend(("--manifest", manifest))
                    result, _, error = self.run_main(*arguments)
                    self.assertEqual(1, result)
                    self.assertIn("duplicate experiment row", error)
                    self.assertIn(f"{second}:2", error)
                    self.assertIn(f"first seen at {first}:2", error)

    def test_disjoint_shards_with_and_without_manifest(self):
        first = self.write_csv("first.csv")
        second = self.write_csv(
            "second.csv", rows=[valid_row(scheme="ldpc", D=20)]
        )
        for manifest in (None, self.write_manifest(schemes=["dense", "ldpc"])):
            with self.subTest(manifest=manifest):
                arguments = [first, second]
                if manifest is not None:
                    arguments.extend(("--manifest", manifest))
                result, output, error = self.run_main(*arguments)
                self.assertEqual((0, ""), (result, error))
                self.assertIn("2 precode result file(s), 2 row(s)", output)

    def test_grid_trials_and_manifest(self):
        result_file = self.write_csv()
        manifest = self.write_manifest()
        result, _, error = self.run_main(result_file, "--manifest", manifest)
        self.assertEqual((0, ""), (result, error))

        result, _, error = self.run_main(
            result_file, "--manifest", manifest, "--expect-trials", "11"
        )
        self.assertEqual(1, result)
        self.assertIn("incompatible with manifest", error)

        manifest.write_text(json.dumps({
            "K": [1000], "schemes": ["dense", "ldpc"],
            "oh": [0], "trials": 10,
        }), encoding="utf-8")
        result, _, error = self.run_main(result_file, "--manifest", manifest)
        self.assertEqual(1, result)
        self.assertIn("grid mismatch", error)

    def test_duplicate_and_malformed_manifest(self):
        result_file = self.write_csv()
        manifest = self.root / "manifest.json"
        for contents in (
            "not json",
            json.dumps({"K": [1000], "schemes": ["dense"], "oh": [0]}),
            json.dumps({
                "K": [1000, 1000], "schemes": ["dense"],
                "oh": [0], "trials": 10,
            }),
        ):
            with self.subTest(contents=contents):
                manifest.write_text(contents, encoding="utf-8")
                result, _, _ = self.run_main(result_file, "--manifest", manifest)
                self.assertEqual(1, result)

    def write_manifest(self, schemes=None):
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({
            "K": [1000], "schemes": schemes or ["dense"],
            "oh": [0], "trials": 10,
        }), encoding="utf-8")
        return manifest


if __name__ == "__main__":
    unittest.main()
