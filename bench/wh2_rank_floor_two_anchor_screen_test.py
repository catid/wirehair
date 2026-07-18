#!/usr/bin/env python3
"""Focused regressions for shared rank-floor thermal parsing."""

from __future__ import annotations

import tempfile
from pathlib import Path
import unittest
from unittest import mock

import wh2_rank_floor_two_anchor_screen as screen


def thermal_row(
    monotonic_s: float,
    *,
    missing_dimms: tuple[int, ...] = (),
    dimm_read_errors: int = 0,
    edac_ce: int = 0,
    edac_ue: int = 0,
) -> bytes:
    dimms = [f"{50.0 + index:.2f}" for index in range(8)]
    for index in missing_dimms:
        dimms[index] = ""
    fields = [
        "2026-07-18T00:00:00.000Z",
        f"{monotonic_s:.6f}",
        "100.0",
        "3600.0",
        "60.0",
        *dimms,
        str(dimm_read_errors),
        "128.0",
        "128.0",
        "128.0",
        str(edac_ce),
        str(edac_ue),
    ]
    if len(fields) != len(screen.THERMAL_FIELDS):
        raise AssertionError("thermal fixture field count changed")
    return (",".join(fields) + "\n").encode("ascii")


class ThermalParserTest(unittest.TestCase):
    def test_stable_readers_reject_same_inode_change_after_final_fstat(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.bin"
            real_stat = screen.os.stat
            for reader in (screen.sha256_file, screen.stable_file_bytes):
                with self.subTest(reader=reader.__name__):
                    path.write_bytes(b"old\n")
                    changed = False

                    def mutate_then_stat(*args, **kwargs):
                        nonlocal changed
                        if not changed:
                            changed = True
                            path.write_bytes(b"new\n")
                        return real_stat(*args, **kwargs)

                    with mock.patch.object(
                            screen.os, "stat", side_effect=mutate_then_stat), \
                         self.assertRaisesRegex(ValueError, "changed"):
                        reader(path)
                    self.assertTrue(changed)

    def test_strict_default_rejects_missing_dimm_temperature(self) -> None:
        row = thermal_row(
            100.0, missing_dimms=(3,), dimm_read_errors=1)
        with self.assertRaisesRegex(ValueError, "missing DIMM temperature"):
            screen.parse_thermal_sample(row, "strict fixture")

    def test_permissive_mode_preserves_real_dimm_error_row(self) -> None:
        row = thermal_row(
            100.0, missing_dimms=(3,), dimm_read_errors=1)
        sample = screen.parse_thermal_sample(
            row, "permissive fixture", allow_dimm_read_errors=True)
        self.assertEqual(sample["raw"], row)
        self.assertEqual(sample["dimm_read_errors"], 1)
        self.assertEqual(len(sample["dimm_temperatures"]), 7)
        self.assertEqual(max(sample["dimm_temperatures"]), 57.0)

    def test_permissive_mode_requires_exact_error_count_and_one_reader(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not match"):
            screen.parse_thermal_sample(
                thermal_row(
                    100.0, missing_dimms=(2,), dimm_read_errors=2),
                "mismatch fixture", allow_dimm_read_errors=True)
        with self.assertRaisesRegex(ValueError, "no readable DIMM"):
            screen.parse_thermal_sample(
                thermal_row(
                    100.0, missing_dimms=tuple(range(8)),
                    dimm_read_errors=8),
                "all missing fixture", allow_dimm_read_errors=True)

    def test_interval_seals_raw_dimm_error_row_byte_exact(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline = thermal_row(100.0)
        interval = thermal_row(
            101.0, missing_dimms=(3,), dimm_read_errors=1)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            sealed = root / "sealed.csv"
            log.write_bytes(header + baseline)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(
                    log, require_zero_edac=False,
                    require_zero_dimm_errors=False)
                with log.open("ab") as output:
                    output.write(interval)
                summary = screen.thermal_finish(
                    log, mark, sealed, require_zero_edac=False,
                    require_zero_dimm_errors=False)
            self.assertEqual(sealed.read_bytes(), header + baseline + interval)
            self.assertEqual(summary["dimm_read_errors_max"], 1)
            self.assertEqual(summary["dimm_max_c"], 57.0)

    def test_edac_rollback_is_negative_and_raw_interval_remains_exact(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline = thermal_row(100.0, edac_ce=5, edac_ue=2)
        interval = (
            thermal_row(101.0, edac_ce=10, edac_ue=3) +
            thermal_row(102.0, edac_ce=7, edac_ue=3)
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            sealed = root / "sealed.csv"
            log.write_bytes(header + baseline)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(
                    log, require_zero_edac=False,
                    require_zero_dimm_errors=False)
                with log.open("ab") as output:
                    output.write(interval)
                summary = screen.thermal_finish(
                    log, mark, sealed, require_zero_edac=False,
                    require_zero_dimm_errors=False)
            self.assertEqual(summary["edac_ce_delta"], -3)
            self.assertEqual(summary["edac_ue_delta"], 1)
            self.assertEqual(sealed.read_bytes(), header + baseline + interval)


if __name__ == "__main__":
    unittest.main()
