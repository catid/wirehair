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
        with self.assertRaisesRegex(validator.ValidationError, "duplicate experiment row"):
            validator.validate_file(result_file)

    def test_grid_trials_and_manifest(self):
        result_file = self.write_csv()
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({
            "K": [1000], "schemes": ["dense"], "oh": [0], "trials": 10,
        }), encoding="utf-8")
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


if __name__ == "__main__":
    unittest.main()
