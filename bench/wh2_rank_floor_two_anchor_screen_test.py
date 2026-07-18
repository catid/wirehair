#!/usr/bin/env python3
"""Focused regressions for shared rank-floor thermal parsing."""

from __future__ import annotations

import json
import inspect
import os
import signal
import shutil
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
from pathlib import Path
from contextlib import nullcontext
import unittest
from unittest import mock

import wh2_rank_floor_two_anchor_screen as screen
import wh2_rank_floor_two_anchor_allk as allk
import wh2_two_anchor_phase_screen as phase_screen
import wh2_preferred_attempt_search as preferred_search


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
    def test_all_prepare_paths_use_the_shared_fresh_builder(self) -> None:
        for prepare in (
                allk.prepare, phase_screen.prepare, preferred_search.prepare):
            with self.subTest(module=prepare.__module__):
                source = inspect.getsource(prepare)
                self.assertIn("fresh_pinned_benchmark_build", source)
                self.assertNotIn("--clean-first", source)
                self.assertNotIn("cmake_source_directory", source)

    def test_pinned_build_contract_binding_rejects_split_provenance(
            self) -> None:
        pinned = {
            "build_policy": {"fresh_pinned": True},
            "build_command": ["/usr/bin/cmake", "--build", "/fresh"],
            "source": {"commit": "a" * 40},
            "binary_sha256": "b" * 64,
            "legacy_binary_hint": "/legacy/codec/wirehair_v2_bench",
        }
        contract = {
            "build_policy": pinned["build_policy"],
            "build_command": pinned["build_command"],
            "source_commit": pinned["source"]["commit"],
            "binary_sha256": pinned["binary_sha256"],
            "legacy_binary_hint": pinned["legacy_binary_hint"],
        }
        allk.validate_pinned_build_contract_binding(
            contract, pinned, "fixture", require_legacy_hint=True)
        for field, replacement in (
                ("build_command", ["legacy-build"]),
                ("source_commit", "c" * 40),
                ("binary_sha256", "d" * 64),
                ("legacy_binary_hint", "/different/binary")):
            with self.subTest(field=field):
                mutated = dict(contract)
                mutated[field] = replacement
                with self.assertRaisesRegex(
                        allk.CampaignError, "pinned build"):
                    allk.validate_pinned_build_contract_binding(
                        mutated, pinned, "fixture",
                        require_legacy_hint=True)

    def test_signal_guard_replays_signal_queued_during_entry(self) -> None:
        script = "\n".join((
            "import os, signal",
            "import wh2_rank_floor_two_anchor_allk as subject",
            "guard = subject._PinnedBuildSignalGuard()",
            "real_signal = signal.signal",
            "injected = [False]",
            "def install(signum, handler):",
            "    result = real_signal(signum, handler)",
            "    if (not injected[0] and "
            "getattr(handler, '__self__', None) is guard):",
            "        injected[0] = True",
            "        os.kill(os.getpid(), signal.SIGTERM)",
            "    return result",
            "signal.signal = install",
            "with guard:",
            "    guard.activate()",
            "    raise RuntimeError('queued signal was lost')",
        ))
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(Path(__file__).resolve().parent)
        completed = subprocess.run(
            (sys.executable, "-c", script), env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
        )
        self.assertEqual(
            completed.returncode, -signal.SIGTERM,
            completed.stderr.decode("utf-8", errors="replace"))

    def test_signal_guard_preserves_ignored_signal_semantics(self) -> None:
        previous = signal.signal(signal.SIGTERM, signal.SIG_IGN)
        try:
            guard = allk._PinnedBuildSignalGuard()
            with guard:
                guard.activate()
                signal.raise_signal(signal.SIGTERM)
            self.assertEqual(signal.getsignal(signal.SIGTERM), signal.SIG_IGN)
        finally:
            signal.signal(signal.SIGTERM, previous)

    def test_signal_during_popen_registration_reaps_new_session(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-popen-signal-regression-") as temporary:
            root = Path(temporary).resolve()
            script = "\n".join((
                "import os, signal, sys",
                "from pathlib import Path",
                "import wh2_rank_floor_two_anchor_allk as subject",
                "root = Path(sys.argv[1])",
                "scratch = root / 'tmp'",
                "scratch.mkdir()",
                "real_popen = subject.subprocess.Popen",
                "def launch_then_signal(*_args, **kwargs):",
                "    process = real_popen(('/bin/sleep', '60'), **kwargs)",
                "    (root / 'child.pid').write_text(str(process.pid))",
                "    os.kill(os.getpid(), signal.SIGTERM)",
                "    return process",
                "subject.subprocess.Popen = launch_then_signal",
                "guard = subject._PinnedBuildSignalGuard()",
                "with guard:",
                "    guard.activate()",
                "    subject._run_pinned_command(",
                "        ('/bin/true',), os.environ.copy(), scratch, 10.0,",
                "        'signal registration fixture', guard)",
            ))
            environment = os.environ.copy()
            environment["PYTHONPATH"] = str(Path(__file__).resolve().parent)
            completed = subprocess.run(
                (sys.executable, "-c", script, str(root)), env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
            )
            self.assertEqual(
                completed.returncode, -signal.SIGTERM,
                completed.stderr.decode("utf-8", errors="replace"))
            child_pid = int((root / "child.pid").read_text())
            try:
                os.kill(child_pid, 0)
            except ProcessLookupError:
                pass
            else:
                try:
                    os.killpg(child_pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                self.fail("signal-interrupted pinned child survived cleanup")

    def test_pinned_child_mask_timeout_and_output_cap(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-pinned-command-regression-") as temporary:
            root = Path(temporary).resolve()
            scratch = root / "tmp"
            scratch.mkdir()
            environment = os.environ.copy()
            guard = allk._PinnedBuildSignalGuard()
            with guard:
                guard.activate()
                stdout, _stderr = allk._run_pinned_command(
                    (
                        sys.executable, "-c",
                        "import json, signal; print(json.dumps(sorted("
                        "int(s) for s in signal.pthread_sigmask("
                        "signal.SIG_BLOCK, ()))))",
                    ),
                    environment, scratch, 10.0, "child signal-mask fixture",
                    guard,
                )
                child_mask = json.loads(stdout)
                self.assertNotIn(signal.SIGINT, child_mask)
                self.assertNotIn(signal.SIGTERM, child_mask)

                started = time.monotonic()
                with self.assertRaisesRegex(allk.CampaignError, "timed out"):
                    allk._run_pinned_command(
                        ("/bin/sleep", "60"), environment, scratch, 0.05,
                        "TERM-respecting timeout fixture", guard)
                self.assertLess(time.monotonic() - started, 0.75)

                writer = (
                    "import os, signal; "
                    "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
                    "chunk=b'x'*(1024*1024); "
                    "exec('while True:\\n os.write(1, chunk)')"
                )
                started = time.monotonic()
                with mock.patch.object(
                        allk, "PINNED_OUTPUT_MAX_BYTES", 4096), \
                     self.assertRaisesRegex(
                         allk.CampaignError, "bounded capture limit"):
                    allk._run_pinned_command(
                        (sys.executable, "-c", writer), environment, scratch,
                        10.0, "bounded-output fixture", guard)
                self.assertLess(time.monotonic() - started, 0.75)

    def test_prepare_exit_failure_does_not_publish_partial_result(self) -> None:
        class FailAfterYield:
            def __init__(self, build):
                self.build = build

            def __enter__(self):
                return self.build

            def __exit__(self, _kind, _error, _traceback):
                raise allk.CampaignError("forced post-yield validation failure")

        with tempfile.TemporaryDirectory(
                prefix="wh2-prepare-exit-regression-") as temporary:
            root = Path(temporary).resolve()
            binary = root / "wirehair_v2_bench"
            cache = root / "CMakeCache.txt"
            groups = root / "groups.tsv"
            source_cells = root / "source-cells.json"
            thermal = root / "thermal.csv"
            binary.write_bytes(b"fresh benchmark fixture\n")
            cache.write_bytes(b"fresh cache fixture\n")
            groups.write_bytes(b"groups fixture\n")
            source_cells.write_bytes(b"{}\n")
            thermal.write_bytes(b"thermal fixture\n")
            executable = root / "trusted-tool"
            executable.write_bytes(b"#!/bin/sh\nexit 0\n")
            executable.chmod(0o755)
            tool_record = {
                "path": str(executable),
                "sha256": allk.sha256_file(executable),
            }
            record = {
                "build_policy": {"fresh_pinned": True},
                "build_command": ["/usr/bin/cmake", "--build", "/fresh"],
                "legacy_binary_hint": str(binary),
                "binary_sha256": allk.sha256_file(binary),
            }
            build = SimpleNamespace(
                binary=binary, cache=cache, record=record)
            head = "a" * 40
            thermal_mark = {
                "dev": 1, "ino": 2, "offset": 3,
                "edac_ce": 0, "edac_ue": 0,
                "monotonic_s": 1.0, "max_temperature_c": 50.0,
                "baseline_row": b"thermal fixture\n",
            }
            allk_result = root / "allk-result"
            allk_args = SimpleNamespace(
                groups=groups, thermal=thermal, result_dir=allk_result,
                binary=binary, build_workers=1, workers=1, timeout=1.0,
            )
            with mock.patch.object(
                    allk, "_trusted_executable", return_value=tool_record), \
                 mock.patch.object(
                     allk, "tracked_clean_source",
                     return_value=(root, head)), \
                 mock.patch.object(allk, "load_groups", return_value=[]), \
                 mock.patch.object(
                     allk, "thermal_start", return_value=thermal_mark), \
                 mock.patch.object(
                     allk, "validate_pinned_build_record"), \
                 mock.patch.object(
                     allk, "fresh_pinned_benchmark_build",
                     return_value=FailAfterYield(build)), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "post-yield validation"):
                allk.prepare(allk_args)
            self.assertFalse(allk_result.exists())

            phase_result = root / "phase-result"
            phase_args = SimpleNamespace(
                acknowledge_post_selection_only=True,
                source_cells=source_cells, workers=64, build_workers=1,
                thermal=thermal, result_dir=phase_result, binary=binary,
                timeout=1.0,
            )
            with mock.patch.object(
                    allk, "_trusted_executable", return_value=tool_record), \
                 mock.patch.object(
                     phase_screen, "tracked_clean_sources",
                     return_value=(root, head)), \
                 mock.patch.object(
                     phase_screen, "load_source_cohort", return_value={}), \
                 mock.patch.object(
                     phase_screen, "q0_identity",
                     return_value={"passed": True}), \
                 mock.patch.object(
                     phase_screen, "thermal_start", return_value=thermal_mark), \
                 mock.patch.object(
                     allk, "validate_pinned_build_record"), \
                 mock.patch.object(
                     allk, "fresh_pinned_benchmark_build",
                     return_value=FailAfterYield(build)), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "post-yield validation"):
                phase_screen.prepare(phase_args)
            self.assertFalse(phase_result.exists())

            class SuccessfulBuild(FailAfterYield):
                def __exit__(self, _kind, _error, _traceback):
                    return False

            corrupt_groups_result = root / "allk-corrupt-groups-copy"
            corrupt_groups_args = SimpleNamespace(**vars(allk_args))
            corrupt_groups_args.result_dir = corrupt_groups_result
            real_copyfile = shutil.copyfile

            def corrupt_groups_copy(source, destination, *copy_args, **copy_kwargs):
                copied = real_copyfile(
                    source, destination, *copy_args, **copy_kwargs)
                if Path(source) == groups:
                    Path(destination).write_bytes(b"copy-time corruption\n")
                return copied

            def validate_groups_copy(path, _expected=allk.GROUPS_SHA256):
                if Path(path) == groups:
                    return []
                raise allk.CampaignError(
                    "group ledger SHA256 mismatch after copy")

            with mock.patch.object(
                    allk, "_trusted_executable", return_value=tool_record), \
                 mock.patch.object(
                     allk, "tracked_clean_source",
                     return_value=(root, head)), \
                 mock.patch.object(
                     allk, "load_groups", side_effect=validate_groups_copy), \
                 mock.patch.object(
                     allk, "thermal_start", return_value=thermal_mark), \
                 mock.patch.object(
                     allk, "validate_pinned_build_record"), \
                 mock.patch.object(
                     allk, "fresh_pinned_benchmark_build",
                     return_value=SuccessfulBuild(build)), \
                 mock.patch.object(
                     allk.shutil, "copyfile",
                     side_effect=corrupt_groups_copy), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "group ledger SHA256 mismatch"):
                allk.prepare(corrupt_groups_args)
            self.assertFalse(corrupt_groups_result.exists())
            self.assertEqual(
                list(root.glob(
                    f".{corrupt_groups_result.name}.prepare-*")), [])

            for module, original_args in (
                    (allk, allk_args), (phase_screen, phase_args)):
                result = root / f"{module.__name__}-late-failure"
                late_args = SimpleNamespace(**vars(original_args))
                late_args.result_dir = result
                tracked_name = (
                    "tracked_clean_source" if module is allk
                    else "tracked_clean_sources")
                tracked = mock.patch.object(
                    module, tracked_name, return_value=(root, head))
                cohort = (
                    mock.patch.object(allk, "load_groups", return_value=[])
                    if module is allk else
                    mock.patch.object(
                        phase_screen, "load_source_cohort", return_value={})
                )
                with mock.patch.object(
                        allk, "_trusted_executable",
                        return_value=tool_record), tracked, cohort, \
                     mock.patch.object(
                         module, "thermal_start", return_value=thermal_mark), \
                     mock.patch.object(
                         phase_screen, "q0_identity",
                         return_value={"passed": True}), \
                     mock.patch.object(
                         allk, "validate_pinned_build_record"), \
                     mock.patch.object(
                         allk, "fresh_pinned_benchmark_build",
                         return_value=SuccessfulBuild(build)), \
                     mock.patch.object(
                         allk, "verify_frozen_sources_at_commit",
                         side_effect=allk.CampaignError(
                             "forced late source verification failure")), \
                     self.assertRaisesRegex(
                         allk.CampaignError, "late source verification"):
                    module.prepare(late_args)
                self.assertFalse(result.exists())
                self.assertEqual(
                    list(root.glob(f".{result.name}.prepare-*")), [])

            for module, original_args in (
                    (allk, allk_args), (phase_screen, phase_args)):
                result = root / f"{module.__name__}-published"
                success_args = SimpleNamespace(**vars(original_args))
                success_args.result_dir = result
                tracked_name = (
                    "tracked_clean_source" if module is allk
                    else "tracked_clean_sources")
                tracked = mock.patch.object(
                    module, tracked_name, return_value=(root, head))
                cohort = (
                    mock.patch.object(allk, "load_groups", return_value=[])
                    if module is allk else
                    mock.patch.object(
                        phase_screen, "load_source_cohort", return_value={})
                )
                with mock.patch.object(
                        allk, "_trusted_executable",
                        return_value=tool_record), tracked, cohort, \
                     mock.patch.object(
                         module, "thermal_start", return_value=thermal_mark), \
                     mock.patch.object(
                         phase_screen, "q0_identity",
                         return_value={"passed": True}), \
                     mock.patch.object(
                         allk, "validate_pinned_build_record"), \
                     mock.patch.object(
                         allk, "fresh_pinned_benchmark_build",
                         return_value=SuccessfulBuild(build)), \
                     mock.patch.object(
                         allk, "verify_frozen_sources_at_commit"), \
                     mock.patch("builtins.print"):
                    self.assertEqual(module.prepare(success_args), 0)
                prepare_record = json.loads(
                    (result / "prepare.json").read_bytes())
                self.assertEqual(prepare_record["result_dir"], str(result))
                self.assertEqual(
                    prepare_record["run_command"][1],
                    str(result / "frozen" / Path(module.__file__).name))
                self.assertEqual(
                    list(root.glob(f".{result.name}.prepare-*")), [])

    def test_prepare_result_staging_is_noreplace_and_cleans_failures(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-result-staging-regression-") as temporary:
            root = Path(temporary).resolve()
            final = root / "campaign"
            with allk.prepare_result_staging(final) as (staging, canonical):
                self.assertEqual(canonical, final)
                (staging / "complete").write_bytes(b"sealed\n")
                staged_name = staging.name
                self.assertFalse(final.exists())
            self.assertEqual((final / "complete").read_bytes(), b"sealed\n")
            self.assertFalse((root / staged_name).exists())

            failed = root / "failed"
            with self.assertRaisesRegex(RuntimeError, "fixture failure"):
                with allk.prepare_result_staging(failed) as (staging, _final):
                    (staging / "partial").write_bytes(b"discard\n")
                    failed_staging = staging
                    raise RuntimeError("fixture failure")
            self.assertFalse(failed.exists())
            self.assertFalse(failed_staging.exists())

            raced = root / "raced"
            with self.assertRaisesRegex(
                    allk.CampaignError, "target already exists"):
                with allk.prepare_result_staging(raced) as (staging, _final):
                    (staging / "ours").write_bytes(b"ours\n")
                    raced.mkdir()
                    (raced / "theirs").write_bytes(b"theirs\n")
                    raced_staging = staging
            self.assertEqual((raced / "theirs").read_bytes(), b"theirs\n")
            self.assertFalse((raced / "ours").exists())
            self.assertFalse(raced_staging.exists())

            substituted = root / "substituted"
            displaced = root / "displaced-original-staging"
            replacement = None
            with self.assertRaisesRegex(
                    allk.CampaignError, "staging identity changed"):
                with allk.prepare_result_staging(
                        substituted) as (staging, _final):
                    (staging / "owned").write_bytes(b"original\n")
                    staging.rename(displaced)
                    staging.mkdir()
                    replacement = staging
                    (staging / "replacement").write_bytes(b"unowned\n")
            self.assertFalse(substituted.exists())
            self.assertEqual(
                (displaced / "owned").read_bytes(), b"original\n")
            self.assertEqual(
                (replacement / "replacement").read_bytes(), b"unowned\n")
            shutil.rmtree(displaced)
            shutil.rmtree(replacement)

            durability_failed = root / "durability-failed"
            real_fsync = allk.os.fsync
            with mock.patch.object(
                    allk, "_fsync_prepared_tree",
                    side_effect=OSError(5, "fixture EIO")), \
                 self.assertRaises(OSError):
                with allk.prepare_result_staging(
                        durability_failed) as (staging, _final):
                    (staging / "unflushed").write_bytes(b"not committed\n")
                    durability_staging = staging
            self.assertFalse(durability_failed.exists())
            self.assertFalse(durability_staging.exists())

            committed = root / "committed-before-fsync-error"

            def fail_committed_parent_fsync(descriptor):
                if committed.exists():
                    raise OSError(5, "fixture EIO")
                return real_fsync(descriptor)

            with mock.patch.object(
                    allk, "_fsync_prepared_tree",
                    side_effect=allk._prepared_tree_inventory), \
                 mock.patch.object(
                     allk.os, "fsync",
                     side_effect=fail_committed_parent_fsync), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "publication committed"):
                with allk.prepare_result_staging(
                        committed) as (staging, _final):
                    (staging / "complete").write_bytes(b"committed\n")
                    committed_staging = staging
            self.assertEqual(
                (committed / "complete").read_bytes(), b"committed\n")
            self.assertFalse(committed_staging.exists())

            unreadable = root / "unreadable-subtree"
            real_walk = allk.os.walk
            denied_path = None

            def deny_inventory_walk(*walk_args, **walk_kwargs):
                if walk_kwargs.get("onerror") is not None:
                    error = PermissionError(13, "fixture EACCES", denied_path)
                    walk_kwargs["onerror"](error)
                return real_walk(*walk_args, **walk_kwargs)

            with mock.patch.object(
                    allk.os, "walk", side_effect=deny_inventory_walk), \
                 self.assertRaisesRegex(
                    allk.CampaignError, "cannot inspect complete"):
                with allk.prepare_result_staging(
                        unreadable) as (staging, _final):
                    hidden = staging / "hidden"
                    hidden.mkdir()
                    nested = hidden / "nested"
                    nested.mkdir()
                    (nested / "payload").write_bytes(
                        b"must be inventoried\n")
                    nested.chmod(0)
                    hidden.chmod(0)
                    denied_path = str(hidden)
                    unreadable_staging = staging
            self.assertFalse(unreadable.exists())
            self.assertFalse(unreadable_staging.exists())

            changed_after_flush = root / "changed-after-flush"
            real_rename = allk._rename_directory_noreplace

            def mutate_before_rename(
                    source, target, expected_identity, expected_inventory):
                (source / "payload").write_bytes(b"changed after fsync\n")
                return real_rename(
                    source, target, expected_identity, expected_inventory)

            with mock.patch.object(
                    allk, "_rename_directory_noreplace",
                    side_effect=mutate_before_rename), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "changed after its durability flush"):
                with allk.prepare_result_staging(
                        changed_after_flush) as (staging, _final):
                    (staging / "payload").write_bytes(b"flushed bytes\n")
                    changed_staging = staging
            self.assertFalse(changed_after_flush.exists())
            self.assertFalse(changed_staging.exists())

            def failing_mask_restore(_operation, _signals):
                failing_mask_restore.calls += 1
                if failing_mask_restore.calls == 1:
                    return set()
                raise OSError(5, "fixture mask restore EIO")

            class NoopTransactionGuard:
                def __enter__(self):
                    return self

                def __exit__(self, _kind, _error, _traceback):
                    return False

                def activate(self):
                    return None

                def deferred_acquisition(self):
                    return nullcontext()

                def block_for_cleanup(self):
                    return None

            failing_mask_restore.calls = 0
            mask_failure_after_commit = root / "mask-failure-after-commit"
            with mock.patch.object(
                    allk, "_PinnedBuildSignalGuard",
                    return_value=NoopTransactionGuard()), \
                 mock.patch.object(
                    allk.signal, "pthread_sigmask",
                    side_effect=failing_mask_restore), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "publication committed"):
                with allk.prepare_result_staging(
                        mask_failure_after_commit) as (staging, _final):
                    (staging / "complete").write_bytes(b"committed\n")
            self.assertEqual(
                (mask_failure_after_commit / "complete").read_bytes(),
                b"committed\n",
            )

            failing_mask_restore.calls = 0
            mask_failure_before_commit = root / "mask-failure-before-commit"
            with mock.patch.object(
                    allk, "_PinnedBuildSignalGuard",
                    return_value=NoopTransactionGuard()), \
                 mock.patch.object(
                    allk.signal, "pthread_sigmask",
                    side_effect=failing_mask_restore), \
                 mock.patch.object(
                     allk, "_rename_directory_noreplace",
                     side_effect=allk.CampaignError(
                         "fixture precommit rename failure")), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "precommit rename failure"):
                with allk.prepare_result_staging(
                        mask_failure_before_commit) as (staging, _final):
                    (staging / "partial").write_bytes(b"discard\n")
                    mask_failure_staging = staging
            self.assertFalse(mask_failure_before_commit.exists())
            self.assertFalse(mask_failure_staging.exists())

            def interrupted_mask_restore(_operation, _signals):
                interrupted_mask_restore.calls += 1
                if interrupted_mask_restore.calls == 1:
                    return set()
                raise KeyboardInterrupt("fixture pending signal")

            interrupted_mask_restore.calls = 0
            interrupted_after_commit = root / "interrupted-after-commit"
            with mock.patch.object(
                    allk, "_PinnedBuildSignalGuard",
                    return_value=NoopTransactionGuard()), \
                 mock.patch.object(
                    allk.signal, "pthread_sigmask",
                    side_effect=interrupted_mask_restore), \
                 self.assertRaisesRegex(
                     KeyboardInterrupt, "pending signal"):
                with allk.prepare_result_staging(
                        interrupted_after_commit) as (staging, _final):
                    (staging / "complete").write_bytes(b"committed\n")
            self.assertEqual(
                (interrupted_after_commit / "complete").read_bytes(),
                b"committed\n",
            )

    def test_staging_cleanup_and_inventory_edge_failures(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-staging-edge-regression-") as temporary:
            root = Path(temporary).resolve()

            outside = root / "outside-hardlink-target"
            outside.write_bytes(b"shared inode\n")
            outside.chmod(0o444)
            hardlink_final = root / "hardlink-final"
            with self.assertRaisesRegex(
                    allk.CampaignError, "nonunique non-file"):
                with allk.prepare_result_staging(
                        hardlink_final) as (staging, _final):
                    os.link(outside, staging / "linked")
                    hardlink_staging = staging
            self.assertEqual(outside.stat().st_mode & 0o777, 0o444)
            self.assertEqual(outside.read_bytes(), b"shared inode\n")
            self.assertFalse(hardlink_final.exists())
            self.assertFalse(hardlink_staging.exists())

            setup_final = root / "setup-fsync-failure"
            real_fsync = allk.os.fsync

            def fail_staging_root_fsync(descriptor):
                metadata = os.fstat(descriptor)
                candidates = list(root.glob(
                    f".{setup_final.name}.prepare-*"))
                if (not fail_staging_root_fsync.triggered and candidates and
                        any((candidate.stat().st_dev,
                             candidate.stat().st_ino) ==
                            (metadata.st_dev, metadata.st_ino)
                            for candidate in candidates)):
                    fail_staging_root_fsync.triggered = True
                    raise OSError(5, "fixture setup EIO")
                return real_fsync(descriptor)

            fail_staging_root_fsync.triggered = False
            with mock.patch.object(
                    allk.os, "fsync",
                    side_effect=fail_staging_root_fsync), \
                 self.assertRaisesRegex(OSError, "setup EIO"):
                with allk.prepare_result_staging(setup_final):
                    self.fail("setup fsync failure reached the body")
            self.assertFalse(setup_final.exists())
            self.assertEqual(
                list(root.glob(f".{setup_final.name}.prepare-*")), [])

            capture_final = root / "capture-race"
            capture_displaced = root / "capture-owned-displaced"
            capture_replacement = None
            real_stat = allk.os.stat

            def fail_identity_capture(path, *stat_args, **stat_kwargs):
                nonlocal capture_replacement
                if (not fail_identity_capture.triggered and
                        isinstance(path, str) and
                        path.startswith(f".{capture_final.name}.prepare-") and
                        stat_kwargs.get("dir_fd") is not None and
                        stat_kwargs.get("follow_symlinks") is False):
                    fail_identity_capture.triggered = True
                    candidate = root / path
                    candidate.rename(capture_displaced)
                    candidate.mkdir()
                    capture_replacement = candidate
                    (candidate / "unowned").write_bytes(b"preserve\n")
                    raise OSError(5, "fixture identity capture EIO")
                return real_stat(path, *stat_args, **stat_kwargs)

            fail_identity_capture.triggered = False
            with mock.patch.object(
                    allk.os, "stat", side_effect=fail_identity_capture), \
                 self.assertRaisesRegex(
                     OSError, "identity capture EIO") as caught:
                with allk.prepare_result_staging(capture_final):
                    self.fail("identity-capture failure reached the body")
            self.assertTrue(hasattr(
                caught.exception, "prepared_result_cleanup_error"))
            self.assertEqual(
                (capture_replacement / "unowned").read_bytes(), b"preserve\n")
            self.assertTrue(capture_displaced.is_dir())
            shutil.rmtree(capture_replacement)
            shutil.rmtree(capture_displaced)

            inventory_root = root / "inventory-root"
            inventory_root.mkdir()
            (inventory_root / "payload").write_bytes(b"initial\n")
            metadata = inventory_root.lstat()
            inventory_identity = (metadata.st_dev, metadata.st_ino)
            real_walk = allk.os.walk

            def mutate_after_directory_yield(*walk_args, **walk_kwargs):
                first = True
                for row in real_walk(*walk_args, **walk_kwargs):
                    yield row
                    if first:
                        first = False
                        (inventory_root / "late").write_bytes(b"late\n")

            with mock.patch.object(
                    allk.os, "walk",
                    side_effect=mutate_after_directory_yield), \
                 self.assertRaisesRegex(
                     allk.CampaignError,
                     "directory (entries )?changed during inventory"):
                allk._prepared_tree_inventory(
                    inventory_root, inventory_identity)

            cleanup_race_final = root / "cleanup-race-final"
            displaced = root / "cleanup-race-owned"
            replacement = None
            real_remove = allk._remove_owned_staging_directory

            def swap_at_cleanup(
                    parent_fd, staging_fd, staging_name, expected_identity):
                nonlocal replacement
                candidate = root / staging_name
                candidate.rename(displaced)
                candidate.mkdir()
                replacement = candidate
                (candidate / "unowned").write_bytes(b"preserve\n")
                return real_remove(
                    parent_fd, staging_fd, staging_name, expected_identity)

            with mock.patch.object(
                    allk, "_remove_owned_staging_directory",
                    side_effect=swap_at_cleanup), \
                 self.assertRaisesRegex(RuntimeError, "body failure") as caught:
                with allk.prepare_result_staging(
                        cleanup_race_final) as (staging, _final):
                    (staging / "owned").write_bytes(b"ours\n")
                    raise RuntimeError("body failure")
            self.assertTrue(hasattr(
                caught.exception, "prepared_result_cleanup_error"))
            self.assertEqual(
                (replacement / "unowned").read_bytes(), b"preserve\n")
            self.assertEqual((displaced / "owned").read_bytes(), b"ours\n")
            shutil.rmtree(replacement)
            shutil.rmtree(displaced)

    def test_staging_creation_fsyncs_new_ancestors(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-staging-ancestor-regression-") as temporary:
            root = Path(temporary).resolve()
            final = root / "new-a" / "new-b" / "result"
            real_fsync = allk.os.fsync
            flushed: set[tuple[int, int]] = set()

            def record_fsync(descriptor):
                metadata = os.fstat(descriptor)
                flushed.add((metadata.st_dev, metadata.st_ino))
                return real_fsync(descriptor)

            with mock.patch.object(
                    allk.os, "fsync", side_effect=record_fsync):
                with allk.prepare_result_staging(final) as (staging, _final):
                    (staging / "complete").write_bytes(b"sealed\n")
            expected = {
                (path.stat().st_dev, path.stat().st_ino)
                for path in (root, root / "new-a", root / "new-a" / "new-b")
            }
            self.assertTrue(expected.issubset(flushed))
            self.assertEqual((final / "complete").read_bytes(), b"sealed\n")

    def test_durable_parent_creation_closes_next_fd_on_fsync_error(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-durable-parent-fd-regression-") as temporary:
            root = Path(temporary).resolve()
            before = set(os.listdir("/proc/self/fd"))
            real_fsync = allk.os.fsync

            def fail_new_a_fsync(descriptor):
                metadata = os.fstat(descriptor)
                new_a = root / "new-a"
                if (not fail_new_a_fsync.triggered and new_a.exists() and
                        (metadata.st_dev, metadata.st_ino) ==
                        (new_a.stat().st_dev, new_a.stat().st_ino)):
                    fail_new_a_fsync.triggered = True
                    raise OSError(5, "fixture ancestor fsync EIO")
                return real_fsync(descriptor)

            fail_new_a_fsync.triggered = False
            with mock.patch.object(
                    allk.os, "fsync", side_effect=fail_new_a_fsync), \
                 self.assertRaisesRegex(
                     allk.CampaignError, "ancestor fsync EIO"):
                allk.open_durable_directory(
                    root / "new-a" / "new-b", create=True,
                    reprove_existing=True)
            self.assertEqual(set(os.listdir("/proc/self/fd")), before)
            self.assertTrue((root / "new-a").is_dir())
            self.assertFalse((root / "new-a" / "new-b").exists())

            flushed: set[tuple[int, int]] = set()

            def record_retry_fsync(descriptor):
                metadata = os.fstat(descriptor)
                flushed.add((metadata.st_dev, metadata.st_ino))
                return real_fsync(descriptor)

            with mock.patch.object(
                    allk.os, "fsync", side_effect=record_retry_fsync):
                descriptor = allk.open_durable_directory(
                    root / "new-a" / "new-b", create=True,
                    reprove_existing=True)
                os.close(descriptor)
            expected = {
                (path.stat().st_dev, path.stat().st_ino)
                for path in (root, root / "new-a", root / "new-a" / "new-b")
            }
            self.assertTrue(expected.issubset(flushed))

    def test_sigterm_cleans_outer_result_staging(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-staging-sigterm-regression-") as temporary:
            root = Path(temporary).resolve()
            for nested in (False, True):
                with self.subTest(nested_fresh_guard=nested):
                    final = root / ("nested-result" if nested else "result")
                    if nested:
                        body = (
                            "    with _PinnedBuildSignalGuard() as inner:\n"
                            "        inner.activate()\n"
                            "        print(staging, flush=True)\n"
                            "        time.sleep(60)\n"
                        )
                    else:
                        body = (
                            "    print(staging, flush=True)\n"
                            "    time.sleep(60)\n"
                        )
                    child = (
                        "import sys, time\n"
                        "from pathlib import Path\n"
                        "from wh2_rank_floor_two_anchor_allk import "
                        "prepare_result_staging, _PinnedBuildSignalGuard\n"
                        "final = Path(sys.argv[1])\n"
                        "with prepare_result_staging(final) as "
                        "(staging, _final):\n" + body
                    )
                    environment = os.environ.copy()
                    environment["PYTHONPATH"] = str(
                        Path(__file__).resolve().parent)
                    process = subprocess.Popen(
                        (sys.executable, "-c", child, str(final)),
                        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True, env=environment,
                        start_new_session=True,
                    )
                    try:
                        self.assertIsNotNone(process.stdout)
                        staging_text = process.stdout.readline().strip()
                        self.assertTrue(staging_text)
                        staging = Path(staging_text)
                        os.kill(process.pid, signal.SIGTERM)
                        returncode = process.wait(timeout=5.0)
                        self.assertEqual(returncode, -signal.SIGTERM)
                        self.assertFalse(final.exists())
                        self.assertFalse(staging.exists())
                    finally:
                        if process.poll() is None:
                            os.killpg(process.pid, signal.SIGKILL)
                            process.wait(timeout=5.0)
                        if process.stdout is not None:
                            process.stdout.close()
                        if process.stderr is not None:
                            process.stderr.close()

    def test_frozen_controller_sources_must_match_immutable_commit(self) -> None:
        if shutil.which("git") is None:
            self.skipTest("git is required for immutable-source integration")
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
                "BUILD_SHARED_LIBS:BOOL": "OFF",
                "BUILD_TESTS:BOOL": "ON",
                "BUILD_CODEC_V2:BOOL": "ON",
                "WIREHAIR_BUILD_BOTH:BOOL": "OFF",
                "WIREHAIR_STATIC_PIC:BOOL": "ON",
                "WIREHAIR_BUILD_TOOLS:BOOL": "OFF",
                "WIREHAIR_BUILD_BENCHMARKS:BOOL": "ON",
                "WIREHAIR_ENABLE_SCHEDULED_TESTS:BOOL": "OFF",
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
                "CMAKE_STATIC_LINKER_FLAGS:STRING": "",
                "CMAKE_STATIC_LINKER_FLAGS_RELEASE:STRING": "",
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

    def test_fresh_compiler_probe_rejects_unsupported_driver_family(
            self) -> None:
        compiler = Path(sys.executable).resolve(strict=True)
        clang_predefines = (
            b"#define __GNUC__ 4\n"
            b"#define __clang__ 1\n"
        )
        with tempfile.TemporaryDirectory() as temporary, \
             mock.patch.object(
                 allk, "_run_pinned_command",
                 return_value=(clang_predefines, b"")), \
             self.assertRaisesRegex(
                 allk.CampaignError, "require a native GCC toolchain"):
            allk._compiler_identity(
                compiler, {"PATH": "/usr/bin:/bin"}, Path(temporary),
                native_probe=True, signal_guard=mock.sentinel.guard)

    def test_fresh_pinned_build_ignores_legacy_graph_and_parent_env(
            self) -> None:
        required_tools = (
            "git", "cmake", "taskset", "ninja", "cc", "c++", "ar",
            "ranlib", "ld", "nm", "objcopy", "objdump", "strip",
            "readelf",
        )
        tool_paths = {
            name: shutil.which(name, path=allk.PINNED_BUILD_PATH)
            for name in required_tools
        }
        missing = tuple(
            name for name in required_tools if tool_paths[name] is None)
        if missing:
            self.skipTest(
                "fresh-build integration lacks tools: " + ", ".join(missing))
        for compiler_name, language in (("cc", "c"), ("c++", "c++")):
            compiler = tool_paths[compiler_name]
            self.assertIsNotNone(compiler)
            completed = subprocess.run(
                (str(compiler), "-dM", "-E", "-x", language, "/dev/null"),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=10,
            )
            if (completed.returncode != 0 or completed.stderr or
                    b"#define __clang__" in completed.stdout or
                    b"#define __GNUC__ " not in completed.stdout):
                self.skipTest(
                    "fresh-build integration requires native GCC cc and c++")
        repo = Path(__file__).resolve().parent.parent
        head = subprocess.check_output(
            (str(tool_paths["git"]), "--no-replace-objects", "-C", str(repo),
             "rev-parse", "HEAD"),
            text=True,
        ).strip()
        with tempfile.TemporaryDirectory(
                prefix="wh2-fresh-pinned-regression-") as temporary:
            root = Path(temporary).resolve()
            legacy = root / "legacy"
            (legacy / "codec").mkdir(parents=True)
            hint = legacy / "codec/wirehair_v2_bench"
            hint.write_bytes(b"caller-controlled legacy binary\n")
            hint.chmod(0o755)
            entries = {
                **allk.PINNED_PROJECT_CACHE,
                **allk.PINNED_LANGUAGE_CACHE,
                "CMAKE_HOME_DIRECTORY": ("INTERNAL", str(repo)),
                "CMAKE_C_COMPILER": (
                    "FILEPATH", str(Path(str(tool_paths["cc"])).resolve())),
                "CMAKE_CXX_COMPILER": (
                    "FILEPATH", str(Path(str(tool_paths["c++"])).resolve())),
                "CMAKE_GENERATOR": ("INTERNAL", "Ninja"),
            }
            (legacy / "CMakeCache.txt").write_text("".join(
                f"{key}:{kind}={value}\n"
                for key, (kind, value) in entries.items()),
                encoding="utf-8",
            )
            marker = root / "legacy-graph-ran"
            (legacy / "build.ninja").write_text(
                "rule hostile\n"
                f"  command = touch {marker}\n"
                "build codec/wirehair_v2_bench: hostile\n",
                encoding="utf-8",
            )
            poison = root / "poison-compiler"
            poison.write_text(
                "#!/bin/sh\n"
                f"touch {root / 'poison-compiler-ran'}\n"
                "exit 99\n",
                encoding="utf-8",
            )
            poison.chmod(0o755)
            trusted_tools = {
                name: allk._trusted_executable(name)
                for name in ("git", "cmake", "taskset")
            }
            poisoned_environment = {
                "PATH": str(root),
                "CC": str(poison), "CXX": str(poison),
                "CFLAGS": "-DWH_SEED_KNOBS=1",
                "CXXFLAGS": "-fplugin=/tmp/untrusted.so",
                "LDFLAGS": "-Wl,--wrap=malloc",
                "CMAKE_TOOLCHAIN_FILE": str(root / "untrusted.cmake"),
                "CMAKE_PROJECT_INCLUDE": str(root / "untrusted.cmake"),
                "CMAKE_CXX_COMPILER_LAUNCHER": str(poison),
                "CMAKE_GENERATOR": "Unix Makefiles",
                "COMPILER_PATH": str(root), "CPATH": str(root),
                "LIBRARY_PATH": str(root),
                "MAKEFLAGS": f"-f {root / 'untrusted.mk'}",
                "BASH_ENV": str(root / "untrusted.sh"),
                "ENV": str(root / "untrusted.sh"),
                "LD_PRELOAD": str(root / "untrusted.so"),
            }
            with mock.patch.dict(
                    os.environ, poisoned_environment, clear=False):
                hostile_git = allk._identity_for_executable(
                    poison, "hostile git fixture")
                with self.assertRaisesRegex(
                        allk.CampaignError, "trusted system git"):
                    with allk.fresh_pinned_benchmark_build(
                            repo, head, hint, hostile_git,
                            trusted_tools["cmake"], trusted_tools["taskset"],
                            build_workers=1, cpu_set=None):
                        self.fail("hostile caller-supplied git was accepted")
                self.assertFalse((root / "poison-compiler-ran").exists())
                with allk.fresh_pinned_benchmark_build(
                        repo, head, hint,
                        trusted_tools["git"], trusted_tools["cmake"],
                        trusted_tools["taskset"],
                        build_workers=min(32, os.cpu_count() or 1),
                        cpu_set=None) as build:
                    self.assertFalse(marker.exists())
                    self.assertFalse((root / "poison-compiler-ran").exists())
                    self.assertTrue(build.record["legacy_graph_used"] is False)
                    self.assertNotEqual(
                        build.binary.read_bytes(), hint.read_bytes())
                    self.assertEqual(
                        set(build.record["environment"]),
                        {"PATH", "HOME", "XDG_CONFIG_HOME", "TMPDIR",
                         "LANG", "LC_ALL", "TZ", "SOURCE_DATE_EPOCH"},
                    )
                    command_text = "\0".join(
                        build.record["configure_command"] +
                        build.record["build_command"])
                    self.assertNotIn(str(legacy), command_text)
                    mutated_record = json.loads(json.dumps(build.record))
                    mutated_record["environment"]["PATH"] = str(root)
                    with self.assertRaisesRegex(
                            allk.CampaignError,
                            "environment policy changed"):
                        allk.validate_pinned_build_record(
                            mutated_record, build.cache)
                    materialized = Path(
                        build.record["source"]["materialized_path"])
                    self.assertTrue(materialized.is_dir())
            self.assertFalse(marker.exists())
            self.assertFalse((root / "poison-compiler-ran").exists())
            self.assertFalse(materialized.exists())

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
