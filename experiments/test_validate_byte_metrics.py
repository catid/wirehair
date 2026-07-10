#!/usr/bin/env python3
"""Unit tests for legacy/v2 byte-metric schema validation."""

import csv
import tempfile
import unittest
from pathlib import Path

try:
    from . import validate_byte_metrics as validator
except ImportError:
    import validate_byte_metrics as validator


class ByteMetricValidatorTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write(self, header, row, name="metrics.csv"):
        path = self.root / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            writer.writerow(row)
        return path

    def xor_row(self, version=2):
        usec = 1.0
        row = {
            "block_bytes": "1280", "working_blocks": "2",
            "ops_per_repeat": "4", "repeats": "1", "median_total_us": "4",
            "median_us_per_xor": str(usec), "median_gib_per_s": "2.384186",
            "checksum": "0x0123456789abcdef",
        }
        if version == 2:
            row.update({
                "schema_version": "2", "logical_bytes_per_op": "1280",
                "estimated_read_bytes_per_op": "2560",
                "estimated_write_bytes_per_op": "1280",
                "logical_gib_per_s": "1.192093",
                "estimated_traffic_gib_per_s": "3.576279",
            })
        return row

    def test_accepts_legacy_and_v2_xor(self):
        legacy = self.write(validator.XOR_V1, self.xor_row(1), "legacy.csv")
        current = self.write(validator.XOR_V2, self.xor_row(2), "current.csv")
        self.assertEqual(("xor", 1, 1), validator.validate_file(legacy))
        self.assertEqual(("xor", 2, 1), validator.validate_file(current))

    def test_rejects_nonfinite_wrong_width_and_bad_ledger(self):
        row = self.xor_row(2)
        row["logical_gib_per_s"] = "nan"
        with self.assertRaises(validator.ValidationError):
            validator.validate_file(self.write(validator.XOR_V2, row, "nan.csv"))

        row = self.xor_row(2)
        row["estimated_read_bytes_per_op"] = "1280"
        with self.assertRaisesRegex(validator.ValidationError, "ledger"):
            validator.validate_file(self.write(validator.XOR_V2, row, "ledger.csv"))

        truncated = self.root / "truncated.csv"
        truncated.write_text(",".join(validator.XOR_V2) + "\n1280,2\n",
                             encoding="utf-8")
        with self.assertRaisesRegex(validator.ValidationError, "field count"):
            validator.validate_file(truncated)

    def test_mixed_trace_uses_explicit_read_write_totals(self):
        row = {
            "case": "fixture", "block_bytes": "10", "recovery_blocks": "2",
            "input_blocks": "1", "ops": "8", "mode": "untiled",
            "tile_bytes": "0", "median_ms": "2", "logical_gib": "0.000000075",
            "memory_gib_3stream": "0.000000225", "gib_per_s": "0.0001125",
            "checksum": "0x0123456789abcdef", "schema_version": "2",
            "logical_work_gib": "0.000000075",
            "estimated_read_gib": "0.000000112",
            "estimated_write_gib": "0.000000075",
            "estimated_traffic_gib": "0.000000187",
            "logical_work_gib_per_s": "0.0000375",
            "estimated_traffic_gib_per_s": "0.0000935",
        }
        path = self.write(validator.TILING_TRACE_V2, row)
        self.assertEqual(("tiling-trace", 2, 1), validator.validate_file(path))


if __name__ == "__main__":
    unittest.main()
