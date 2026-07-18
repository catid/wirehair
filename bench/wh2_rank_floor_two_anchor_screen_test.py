#!/usr/bin/env python3
"""Focused regressions for shared rank-floor thermal parsing."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
import unittest
from unittest import mock

import wh2_rank_floor_two_anchor_screen as screen
import wh2_rank_floor_two_anchor_allk as allk


def thermal_row(
    monotonic_s: float,
    *,
    cpu_tctl_c: float = 60.0,
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
        str(cpu_tctl_c),
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

    def test_oversized_thermal_row_is_never_accepted_as_a_suffix(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        row = thermal_row(100.0)
        oversized = b"x" * 20000 + row[row.index(b","):]
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "thermal.csv"
            log.write_bytes(header + oversized)
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(ValueError, "row size limit"):
                screen.thermal_start(log, require_zero_edac=False)

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

    def test_interval_finish_retries_growth_after_read(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline = thermal_row(100.0)
        interval = thermal_row(101.0)
        late = thermal_row(102.0, edac_ce=1)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            sealed = root / "sealed.csv"
            log.write_bytes(header + baseline)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(
                    log, require_zero_edac=False)
            with log.open("ab") as output:
                output.write(interval)
            real_validate = screen.validate_thermal_log_descriptor
            calls = 0

            def append_before_post_read_stat(source, path):
                nonlocal calls
                calls += 1
                if calls == 3:
                    with log.open("ab") as output:
                        output.write(late)
                return real_validate(source, path)

            with mock.patch.object(screen, "validate_thermal_current"), \
                 mock.patch.object(
                     screen, "validate_thermal_log_descriptor",
                     side_effect=append_before_post_read_stat):
                summary = screen.thermal_finish(
                    log, mark, sealed, require_zero_edac=False)
            self.assertGreaterEqual(calls, 6)
            self.assertEqual(summary["edac_ce_delta"], 1)
            self.assertEqual(
                sealed.read_bytes(), header + baseline + interval + late)

    def test_interval_finish_rejects_same_inode_baseline_rewrite(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline = thermal_row(100.0)
        changed = baseline.replace(b",60.0,", b",61.0,", 1)
        self.assertEqual(len(changed), len(baseline))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            sealed = root / "sealed.csv"
            log.write_bytes(header + baseline)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(log, require_zero_edac=False)
            with log.open("r+b") as output:
                output.seek(len(header))
                output.write(changed)
                output.seek(0, os.SEEK_END)
                output.write(thermal_row(101.0))
            with self.assertRaisesRegex(ValueError, "baseline row changed"):
                screen.thermal_finish(
                    log, mark, sealed, require_zero_edac=False)
            self.assertFalse(sealed.exists())

    def test_interval_finish_rejects_inplace_rewrite_at_final_stat(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline = thermal_row(100.0)
        changed = baseline.replace(b",60.0,", b",61.0,", 1)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            sealed = root / "sealed.csv"
            log.write_bytes(header + baseline)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(log, require_zero_edac=False)
            with log.open("ab") as output:
                output.write(thermal_row(101.0))
            real_validate = screen.validate_thermal_log_descriptor
            calls = 0

            def mutate_at_final_stat(source, path):
                nonlocal calls
                calls += 1
                if calls == 3:
                    with log.open("r+b") as output:
                        output.seek(len(header))
                        output.write(changed)
                return real_validate(source, path)

            with mock.patch.object(
                    screen, "validate_thermal_log_descriptor",
                    side_effect=mutate_at_final_stat), \
                 self.assertRaisesRegex(ValueError, "baseline row changed"):
                screen.thermal_finish(
                    log, mark, sealed, require_zero_edac=False)
            self.assertGreaterEqual(calls, 3)
            self.assertFalse(sealed.exists())


class AllKProvenanceTest(unittest.TestCase):
    def test_frozen_controller_sources_must_match_immutable_commit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo = Path(temporary).resolve() / "source"
            repo.mkdir()
            subprocess.run(("git", "init", "-q", str(repo)), check=True)
            subprocess.run(
                ("git", "-C", str(repo), "config", "user.name", "WH2 Test"),
                check=True,
            )
            subprocess.run(
                (
                    "git", "-C", str(repo), "config", "user.email",
                    "wh2-test@example.invalid",
                ),
                check=True,
            )
            source = repo / "controller.py"
            other = repo / "other.py"
            source.write_bytes(b"committed controller\n")
            other.write_bytes(b"other committed controller\n")
            subprocess.run(
                ("git", "-C", str(repo), "add", source.name, other.name),
                check=True)
            subprocess.run(
                ("git", "-C", str(repo), "commit", "-qm", "fixture"),
                check=True,
            )
            commit = subprocess.check_output(
                ("git", "-C", str(repo), "rev-parse", "HEAD"), text=True,
            ).strip()
            frozen = repo.parent / "frozen-controller.py"
            frozen.write_bytes(source.read_bytes())
            git = Path(shutil.which("git") or "/usr/bin/git").resolve(
                strict=True)
            allk.verify_frozen_sources_at_commit(
                repo, commit, ((source, frozen),), git)
            frozen.write_bytes(b"transient worktree bytes\n")
            with self.assertRaisesRegex(
                    allk.CampaignError, "differs from immutable"):
                allk.verify_frozen_sources_at_commit(
                    repo, commit, ((source, frozen),), git)

            source.write_bytes(b"replacement commit controller\n")
            subprocess.run(
                ("git", "-C", str(repo), "add", source.name), check=True)
            subprocess.run(
                ("git", "-C", str(repo), "commit", "-qm", "replacement"),
                check=True,
            )
            replacement = subprocess.check_output(
                ("git", "-C", str(repo), "rev-parse", "HEAD"), text=True,
            ).strip()
            subprocess.run(
                ("git", "-C", str(repo), "replace", commit, replacement),
                check=True,
            )
            frozen.write_bytes(b"committed controller\n")
            allk.verify_frozen_sources_at_commit(
                repo, commit, ((source, frozen),), git)
            frozen.write_bytes(b"replacement commit controller\n")
            with self.assertRaisesRegex(
                    allk.CampaignError, "differs from immutable"):
                allk.verify_frozen_sources_at_commit(
                    repo, commit, ((source, frozen),), git)

            source.unlink()
            source.symlink_to(other.name)
            frozen.write_bytes(other.read_bytes())
            with self.assertRaisesRegex(
                    allk.CampaignError, "differs from immutable"):
                allk.verify_frozen_sources_at_commit(
                    repo, commit, ((source, frozen),), git)

    def test_live_guard_bridges_high_streak_across_its_mark(self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        now = time.monotonic()
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "thermal.csv"
            log.write_bytes(
                header +
                thermal_row(now - 0.2, cpu_tctl_c=91.0) +
                thermal_row(now - 0.1, cpu_tctl_c=91.0)
            )
            abort = allk.threading.Event()
            guard = allk.ThermalGuard(log, abort, limit_c=90.0)
            self.assertEqual(guard.initial_consecutive_high, 2)
            guard.start()
            try:
                with log.open("ab") as output:
                    output.write(thermal_row(now, cpu_tctl_c=91.0))
                self.assertTrue(abort.wait(2.0))
                self.assertIn("3 consecutive samples", guard.error or "")
            finally:
                guard.stop_event.set()
                guard.thread.join(timeout=2.0)

    def test_prepare_anchor_rejects_resealed_staged_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            frozen = root / "frozen"
            frozen.mkdir()
            staged = frozen / "staged.sha256"
            staged.write_bytes(b"first frozen manifest\n")
            schema = "test.prepare.v1"
            fields = ("source_commit", "binary_sha256", "staged_sha256")
            record = {
                "schema": schema, "result_dir": str(root),
                "source_commit": "a" * 40, "binary_sha256": "b" * 64,
                "staged_sha256": allk.sha256_file(staged),
            }
            (root / "prepare.json").write_text(
                json.dumps(record) + "\n", encoding="utf-8")
            self.assertEqual(
                allk.validate_prepare_anchor(
                    root, schema, fields, "staged_sha256"),
                record)
            staged.write_bytes(b"self-consistent replacement manifest\n")
            with self.assertRaisesRegex(allk.CampaignError, "prepare anchor"):
                allk.validate_prepare_anchor(
                    root, schema, fields, "staged_sha256")

    def test_allk_build_policy_rejects_equation_knob_flags(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo = Path(temporary).resolve()
            cache = repo / "CMakeCache.txt"
            entries = {
                "CMAKE_BUILD_TYPE:STRING": "Release",
                "BUILD_TESTS:BOOL": "ON",
                "BUILD_CODEC_V2:BOOL": "ON",
                "WIREHAIR_BUILD_BENCHMARKS:BOOL": "ON",
                "MARCH_NATIVE:BOOL": "ON",
                "WIREHAIR_STRICT_WARNINGS:BOOL": "ON",
                "WIREHAIR_ENABLE_LIBFUZZER:BOOL": "OFF",
                "WH_LTO:STRING": "OFF",
                "WH_PGO_MODE:STRING": "OFF",
                "CMAKE_C_FLAGS:STRING": "",
                "CMAKE_CXX_FLAGS:STRING": "",
                "CMAKE_C_FLAGS_RELEASE:STRING": "-O3 -DNDEBUG",
                "CMAKE_CXX_FLAGS_RELEASE:STRING": "-O3 -DNDEBUG",
                "CMAKE_EXE_LINKER_FLAGS:STRING": "",
                "CMAKE_EXE_LINKER_FLAGS_RELEASE:STRING": "",
                "CMAKE_SHARED_LINKER_FLAGS:STRING": "",
                "CMAKE_SHARED_LINKER_FLAGS_RELEASE:STRING": "",
                "CMAKE_MODULE_LINKER_FLAGS:STRING": "",
                "CMAKE_MODULE_LINKER_FLAGS_RELEASE:STRING": "",
                "CMAKE_HOME_DIRECTORY:INTERNAL": str(repo),
                "CMAKE_C_COMPILER:FILEPATH": "/usr/bin/cc",
                "CMAKE_CXX_COMPILER:FILEPATH": "/usr/bin/c++",
                "CMAKE_GENERATOR:INTERNAL": "Ninja",
            }

            def write_cache(cxx_flags: str) -> None:
                current = dict(entries)
                current["CMAKE_CXX_FLAGS:STRING"] = cxx_flags
                cache.write_text("".join(
                    f"{key}={value}\n" for key, value in current.items()),
                    encoding="utf-8")

            write_cache("")
            policy = allk.validate_candidate_build_policy(cache, repo)
            self.assertTrue(policy["march_native"])
            self.assertTrue(policy["strict_warnings"])
            for flag in ("-DWH_SEED_KNOBS", "-DWH_PEELCAP"):
                with self.subTest(flag=flag):
                    write_cache(flag)
                    with self.assertRaisesRegex(
                            allk.CampaignError, "CMAKE_CXX_FLAGS"):
                        allk.validate_candidate_build_policy(cache, repo)

    def test_frozen_executable_accepts_hardlink_and_rejects_indirection(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            executable = root / "taskset"
            executable.write_bytes(b"mock taskset\n")
            executable.chmod(0o755)
            alias = root / "taskset-hardlink"
            os.link(executable, alias)
            record = {
                "path": str(executable),
                "sha256": allk.sha256_file(
                    executable, require_unique=False),
            }
            self.assertEqual(
                allk.frozen_executable_path(record, "taskset"), executable)
            self.assertEqual(executable.stat().st_nlink, 2)

            indirect = root / "taskset-symlink"
            indirect.symlink_to(executable.name)
            with self.assertRaises(allk.CampaignError):
                allk.frozen_executable_path(
                    {**record, "path": str(indirect)}, "taskset")

            binary = root / "wirehair_v2_bench"
            binary.write_bytes(b"mock binary\n")
            binary.chmod(0o755)
            digest = allk.sha256_file(binary)
            allk.verify_frozen_binary(binary, digest)
            binary.write_bytes(b"mutated binary\n")
            with self.assertRaisesRegex(allk.CampaignError, "binary changed"):
                allk.verify_frozen_binary(binary, digest)

    def test_frozen_thermal_history_rejects_replacement_and_transient_edac(
            self) -> None:
        header = (",".join(screen.THERMAL_FIELDS) + "\n").encode("ascii")
        baseline_row = thermal_row(100.0)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            log.write_bytes(header + baseline_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(log, require_zero_edac=False)
            baseline = {
                "dev": mark["dev"], "ino": mark["ino"],
                "offset": mark["offset"], "edac_ce": mark["edac_ce"],
                "edac_ue": mark["edac_ue"],
                "monotonic_s": mark["monotonic_s"],
                "max_temperature_c": mark["max_temperature_c"],
                "row_sha256": allk.sha256_bytes(baseline_row),
            }
            with log.open("ab") as output:
                output.write(thermal_row(101.0))
            with mock.patch.object(screen, "validate_thermal_current"):
                self.assertEqual(
                    allk.validate_frozen_thermal_source(log, baseline, 5.0),
                    (mark["dev"], mark["ino"], 0, 0))

            with log.open("ab") as output:
                output.write(thermal_row(102.0, edac_ce=1))
                output.write(thermal_row(103.0))
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "history changed"):
                allk.validate_frozen_thermal_source(log, baseline, 5.0)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "thermal.csv"
            log.write_bytes(header + baseline_row)
            with mock.patch.object(screen, "validate_thermal_current"):
                mark = screen.thermal_start(log, require_zero_edac=False)
            baseline = {
                "dev": mark["dev"], "ino": mark["ino"],
                "offset": mark["offset"], "edac_ce": 0, "edac_ue": 0,
                "monotonic_s": mark["monotonic_s"],
                "max_temperature_c": mark["max_temperature_c"],
                "row_sha256": allk.sha256_bytes(baseline_row),
            }
            original = root / "thermal-original.csv"
            log.rename(original)
            log.write_bytes(header + baseline_row + thermal_row(101.0))
            with mock.patch.object(screen, "validate_thermal_current"), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "changed since"):
                allk.validate_frozen_thermal_source(log, baseline, 5.0)


if __name__ == "__main__":
    unittest.main()
