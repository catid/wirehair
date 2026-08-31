#!/usr/bin/python3
"""Fake-binary tests for the one-shot fixed-five packet controller."""

from __future__ import annotations

import importlib.util
from dataclasses import replace
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock


MODULE_PATH = Path(__file__).with_name("wh2_fixed5_packet_screen.py")
SPEC = importlib.util.spec_from_file_location("wh2_fixed5_screen", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load fixed5 screen controller")
subject = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = subject
SPEC.loader.exec_module(subject)


FAKE_BENCH = r'''#!/usr/bin/python3
import math
from pathlib import Path
import statistics
import subprocess
import sys
import time

SCENARIO = @SCENARIO@
ROOT = Path(@ROOT@)
COUNTER = Path(@COUNTER@)
CELLS = (
    (0, 10000, 1280, 256), (1, 64000, 32768, 10),
    (2, 64000, 1280, 256), (3, 10000, 32768, 10),
    (4, 10000, 4096, 80), (5, 64000, 4096, 80),
)
PANELS = (
    ("active", "direct_fixed5", 2),
    ("fallback", "direct_fixed5_dispatch", 3),
    ("aa", "legacy_copy", 2),
    ("split", "wide_split2_plus3", 2),
)
T = 2.039513446

def commit():
    return subprocess.check_output(
        ["/usr/bin/git", "-C", str(ROOT), "rev-parse", "HEAD"],
        stderr=subprocess.DEVNULL,
    ).decode("ascii").strip()

def describe():
    source = "0" * 40 if SCENARIO == "source" else commit()
    print("wh2_fixed5_design,version=1,source_commit=%s,rows=256,samples=32,"
          "mix_count=3,peel_seed=4d241359,precode_seed=786f72667573696f,"
          "target_cpu=50,target_full_apic_id=00000064,"
          "target_identity_sha256="
          "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7,"
          "candidate=direct_no_write_fixed5_try,"
          "primary_gate=upper95_lt_0.99,control_band=0.99_to_1.01,"
          "student_t_df31=2.039513446,timed_scope=xor_schedules_only,"
          "timed_counters=0,binary_sha256=controller_required" % source)
    for ordinal, K, B, repetitions in CELLS:
        print("wh2_fixed5_cell,version=1,ordinal=%d,K=%d,block_bytes=%d,"
              "repetitions=%d" % (ordinal, K, B, repetitions))

def stats(pairs):
    logs = [math.log(candidate / legacy) for legacy, candidate in pairs]
    mean = statistics.fmean(logs)
    variance = sum((value - mean) ** 2 for value in logs) / 31
    half = T * math.sqrt(variance / 32)
    return math.exp(mean), math.exp(mean - half), math.exp(mean + half)

def run_cell(arguments):
    ordinal, K, B = map(int, arguments[:3])
    expected_commit = arguments[3]
    target_cpu = int(arguments[4])
    if (CELLS[ordinal][:3] != (ordinal, K, B) or
            expected_commit != commit() or target_cpu != 50):
        return 2
    with COUNTER.open("a", encoding="ascii") as output:
        output.write(str(ordinal) + "\n")
    if SCENARIO == "timeout" and ordinal == 2:
        sys.stdout.write("x" * (2 * 1024 * 1024 + 17))
        sys.stdout.flush()
        time.sleep(3)
        return 1
    if SCENARIO == "exit" and ordinal == 2:
        return 7
    repetitions = CELLS[ordinal][3]
    active_hash = "1" * 16 if K == 10000 else "3" * 16
    if SCENARIO == "row_drift" and ordinal == 2:
        active_hash = "5" * 16
    fallback_hash = "2" * 16 if K == 10000 else "4" * 16
    P = 110 if K == 10000 else 370
    print("wh2_fixed5_meta,version=1,source_commit=%s,ordinal=%d,K=%d,P=%d,"
          "block_bytes=%d,repetitions=%d,rows=256,samples=32,mix_count=3,"
          "peel_seed=4d241359,precode_seed=786f72667573696f,active_degree=2,"
          "fallback_degree=3,active_terms=5,"
          "fallback_terms=6,active_rows_hash=%s,fallback_rows_hash=%s,"
          "input_hash=%016x,active_first_id=0,active_last_id=510,"
          "fallback_first_id=1,fallback_last_id=511,wide_build=1,avx2=1,"
          "avx512=1,candidate=direct_no_write_fixed5_try,direct_wide_available=1,"
          "timed_row_generation=0,timed_pointer_gather=0,timed_counters=0,"
          "panel_order=active_fallback_aa_split" % (
              expected_commit, ordinal, K, P, B, repetitions,
              active_hash, fallback_hash, ordinal + 1))
    if SCENARIO == "partial" and ordinal == 2:
        return 0
    panel_stats = {}
    for panel, candidate_name, degree in PANELS:
        ratio = 1.0
        if panel == "active":
            ratio = 1.0 if SCENARIO in (
                "reject", "reject_missing_target", "exit0_reject") else 0.95
        elif panel == "split":
            ratio = 0.97
        pairs = []
        for sample in range(32):
            legacy = 1000000 + sample * 1009 + ordinal * 17
            candidate = max(1, round(legacy * ratio))
            pairs.append((legacy, candidate))
            order = "AB" if sample % 2 == 0 else "BA"
            print("wh2_fixed5_raw,version=1,cell=%s,candidate=%s,ordinal=%d,"
                  "K=%d,block_bytes=%d,degree=%d,rows=256,repetitions=%d,"
                  "sample=%d,order=%s,legacy_ns=%d,candidate_ns=%d,ratio=%.12f" % (
                      panel, candidate_name, ordinal, K, B, degree, repetitions,
                      sample, order, legacy, candidate, candidate / legacy))
        geomean, lower, upper = stats(pairs)
        panel_stats[panel] = (geomean, lower, upper)
        print("wh2_fixed5_summary,version=1,cell=%s,candidate=%s,ordinal=%d,"
              "K=%d,block_bytes=%d,degree=%d,samples=32,geomean=%.12f,"
              "lower95=%.12f,upper95=%.12f" % (
                  panel, candidate_name, ordinal, K, B, degree,
                  geomean, lower, upper))
    if SCENARIO == "malformed" and ordinal == 2:
        print("unexpected_line")
        return 0
    primary = panel_stats["active"][2] < 0.99
    good_identity = (
        "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7"
    )
    pre_identity = (
        "0" * 64 if SCENARIO == "target" and ordinal == 2 else good_identity
    )
    post_identity = (
        "f" * 64 if SCENARIO == "target_post" and ordinal == 2 else
        pre_identity
    )
    apic = (
        "00000065" if SCENARIO == "target_apic" and ordinal == 2 else
        "00000064"
    )
    target_line = (
        "wh2_fixed5_target,version=1,source_commit=%s,ordinal=%d,K=%d,"
        "block_bytes=%d,target_cpu=50,full_apic_id=%s,"
        "pre_identity_sha256=%s,post_identity_sha256=%s,"
        "canonical_bytes=617,pre_before_cpu=50,pre_after_cpu=50,"
        "post_before_cpu=50,post_after_cpu=50,"
        "pre_before_affinity_count=1,pre_after_affinity_count=1,"
        "post_before_affinity_count=1,post_after_affinity_count=1,"
        "pre_voluntary_delta=0,pre_involuntary_delta=0,"
        "post_voluntary_delta=0,post_involuntary_delta=0,gate=pass" % (
            expected_commit, ordinal, K, B, apic, pre_identity, post_identity)
    )
    omit_target = ordinal == 2 and SCENARIO in (
        "missing_target", "reject_missing_target")
    if not omit_target:
        print(target_line)
        if SCENARIO == "duplicate_target" and ordinal == 2:
            print(target_line)
    print("wh2_fixed5_result,version=1,ordinal=%d,K=%d,block_bytes=%d,"
          "primary=%s,fallback_control=pass,aa_control=pass,controls=pass,"
          "split_promotional=0,mismatch_sink=0,status=%s" % (
              ordinal, K, B, "pass" if primary else "reject",
              "pass" if primary else "reject"))
    if SCENARIO == "exit10_pass" and ordinal == 2:
        return 10
    if SCENARIO == "exit0_reject" and ordinal == 2:
        return 0
    return 0 if primary else 10

if sys.argv[1:] == ["--describe"]:
    describe()
    raise SystemExit(0)
if len(sys.argv) == 7 and sys.argv[1] == "--run-cell":
    raise SystemExit(run_cell(sys.argv[2:]))
raise SystemExit(2)
'''


class Fixture:
    def __init__(self, root: Path, scenario: str = "pass", counters: bool = False):
        self.root = root
        self.source = root / "source"
        self.source.mkdir()
        self.output = root / "evidence"
        self.counter = root / "invocations.txt"
        self.benchmark = self.source / "fake-bench"
        self.taskset = self.source / "fake-taskset"
        self.nm = self.source / "fake-nm"
        bench = FAKE_BENCH.replace("@SCENARIO@", repr(scenario))
        bench = bench.replace("@ROOT@", repr(str(self.source)))
        bench = bench.replace("@COUNTER@", repr(str(self.counter)))
        self.write_executable(self.benchmark, bench)
        self.write_executable(self.taskset, textwrap.dedent(r'''#!/bin/sh
            test "$1" = "-c" || exit 91
            shift 2
            exec "$@"
        '''))
        symbol = "00000000 T gf256_count_calls\n" if counters else ""
        self.write_executable(self.nm, textwrap.dedent('''#!/bin/sh
            printf '%s'
            printf '00000000 T gf256_try_addset5_wide_mem\\n'
        ''' % symbol))
        self.git("init", "-q")
        self.git("config", "user.email", "fake@example.test")
        self.git("config", "user.name", "Fake Screen")
        self.git("add", ".")
        self.git("commit", "-qm", "fixture")
        self.commit = self.git("rev-parse", "HEAD").strip()

    @staticmethod
    def write_executable(path: Path, text: str) -> None:
        path.write_text(text, encoding="ascii")
        path.chmod(path.stat().st_mode | stat.S_IXUSR)

    def git(self, *arguments: str) -> str:
        result = subprocess.run(
            ["/usr/bin/git", "-C", str(self.source), *arguments],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, text=True,
            env={
                "GIT_CONFIG_GLOBAL": "/dev/null",
                "GIT_CONFIG_NOSYSTEM": "1",
                "GIT_NO_REPLACE_OBJECTS": "1",
                "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
            },
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        return result.stdout

    def config(self, binary_hash: str = "") -> subject.Config:
        return subject.Config(
            benchmark=self.benchmark.resolve(), output_dir=self.output,
            source_root=self.source.resolve(), expected_commit=self.commit,
            expected_binary_sha256=(
                binary_hash or subject.sha256_file(self.benchmark)),
            cpu=subject.TARGET_CPU, taskset=self.taskset.resolve(),
            nm=self.nm.resolve(), git=Path("/usr/bin/git"),
        )

    def invocations(self) -> list[int]:
        if not self.counter.exists():
            return []
        return [int(value) for value in self.counter.read_text().splitlines()]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="ascii"))


class ControllerTests(unittest.TestCase):
    def setUp(self) -> None:
        affinity = mock.patch.object(
            subject.os, "sched_getaffinity", return_value={subject.TARGET_CPU})
        affinity.start()
        self.addCleanup(affinity.stop)

    def test_pass_publishes_complete_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            self.assertEqual(subject.run_campaign(fixture.config()), 0)
            self.assertEqual(fixture.invocations(), list(range(6)))
            self.assertEqual(
                read_json(fixture.output / "summary.json")["status"], "pass")
            self.assertTrue((fixture.output / "COMPLETE").is_file())
            complete = read_json(fixture.output / "COMPLETE")
            self.assertEqual(
                complete["inputs_sha256"],
                subject.sha256_file(fixture.output / "inputs.json"))
            self.assertEqual(
                complete["results_sha256"],
                subject.sha256_file(fixture.output / "results.jsonl"))
            self.assertEqual(
                complete["summary_sha256"],
                subject.sha256_file(fixture.output / "summary.json"))
            receipt_hash = complete.pop("receipt_sha256")
            self.assertEqual(
                receipt_hash,
                subject.sha256_bytes(subject.canonical_bytes(complete)))
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(fixture.config())
            self.assertEqual(fixture.invocations(), list(range(6)))

    def test_attempt_claim_durably_contains_complete_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            config = fixture.config()
            preflight = subject.preflight(config)
            inputs = subject.make_inputs(config, preflight)
            claim = subject.claim_attempt(config, preflight, inputs)
            receipt = read_json(claim)
            self.assertEqual(receipt["inputs"], inputs)
            self.assertEqual(
                receipt["inputs_sha256"],
                subject.sha256_bytes(subject.canonical_bytes(inputs)))
            self.assertFalse(fixture.output.exists())
            self.assertEqual(fixture.invocations(), [])

    def test_reject_runs_every_cell_and_returns_ten(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "reject")
            self.assertEqual(subject.run_campaign(fixture.config()), 10)
            self.assertEqual(fixture.invocations(), list(range(6)))
            summary = read_json(fixture.output / "summary.json")
            self.assertEqual(summary["status"], "reject")
            self.assertEqual(summary["outcomes"], ["reject"] * 6)

    def test_malformed_output_is_terminal_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "malformed")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), [0, 1, 2])
            self.assertEqual(read_json(fixture.output / "summary.json")["status"], "invalid")

    def test_partial_output_is_terminal_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "partial")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), [0, 1, 2])
            self.assertTrue((fixture.output / "cell-02.stdout").is_file())

    def test_existing_output_directory_prevents_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            fixture.output.mkdir()
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(fixture.config())
            self.assertEqual(fixture.invocations(), [])

    def test_source_describe_mismatch_fails_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "source")
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(fixture.config())
            self.assertEqual(fixture.invocations(), [])
            self.assertFalse(fixture.output.exists())

    def test_binary_hash_mismatch_fails_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(fixture.config("0" * 64))
            self.assertEqual(fixture.invocations(), [])

    def test_counter_symbol_fails_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), counters=True)
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(fixture.config())
            self.assertEqual(fixture.invocations(), [])

    def test_target_receipt_failures_are_terminal_invalid(self) -> None:
        for scenario in (
                "target", "target_post", "target_apic", "missing_target",
                "duplicate_target"):
            with self.subTest(scenario=scenario), \
                    tempfile.TemporaryDirectory() as directory:
                fixture = Fixture(Path(directory), scenario)
                self.assertEqual(subject.run_campaign(fixture.config()), 1)
                self.assertEqual(fixture.invocations(), [0, 1, 2])
                summary = read_json(fixture.output / "summary.json")
                self.assertEqual(summary["status"], "invalid")

    def test_reject_without_target_receipt_is_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "reject_missing_target")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), [0, 1, 2])
            self.assertEqual(
                read_json(fixture.output / "summary.json")["status"],
                "invalid")

    def test_same_k_row_drift_is_terminal_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "row_drift")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), list(range(6)))
            self.assertEqual(
                read_json(fixture.output / "summary.json")["status"],
                "invalid")

    def test_exit_result_contradictions_are_terminal_invalid(self) -> None:
        for scenario in ("exit10_pass", "exit0_reject"):
            with self.subTest(scenario=scenario), \
                    tempfile.TemporaryDirectory() as directory:
                fixture = Fixture(Path(directory), scenario)
                self.assertEqual(subject.run_campaign(fixture.config()), 1)
                self.assertEqual(fixture.invocations(), [0, 1, 2])
                self.assertEqual(
                    read_json(fixture.output / "summary.json")["status"],
                    "invalid")

    def test_nonfrozen_cpu_fails_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory))
            config = replace(fixture.config(), cpu=subject.TARGET_CPU + 1)
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(config)
            self.assertEqual(fixture.invocations(), [])

    def test_dirty_tracked_or_untracked_source_fails_preflight(self) -> None:
        for tracked in (True, False):
            with self.subTest(tracked=tracked), \
                    tempfile.TemporaryDirectory() as directory:
                fixture = Fixture(Path(directory))
                if tracked:
                    fixture.taskset.write_text(
                        fixture.taskset.read_text(encoding="ascii") +
                        "# dirty\n",
                        encoding="ascii")
                else:
                    (fixture.source / "untracked").write_text(
                        "dirty\n", encoding="ascii")
                with self.assertRaises(subject.ScreenError):
                    subject.run_campaign(fixture.config())
                self.assertEqual(fixture.invocations(), [])

    def test_output_under_source_or_symlink_parent_fails_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = Fixture(root)
            under_source = replace(
                fixture.config(), output_dir=fixture.source / "evidence")
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(under_source)
            real_parent = root / "real-parent"
            real_parent.mkdir()
            linked_parent = root / "linked-parent"
            linked_parent.symlink_to(real_parent, target_is_directory=True)
            through_symlink = replace(
                fixture.config(), output_dir=linked_parent / "evidence")
            with self.assertRaises(subject.ScreenError):
                subject.run_campaign(through_symlink)
            self.assertEqual(fixture.invocations(), [])

    def test_unaccepted_exit_is_captured_and_terminal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = Fixture(Path(directory), "exit")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), [0, 1, 2])
            records = [json.loads(line) for line in
                       (fixture.output / "results.jsonl").read_text().splitlines()]
            self.assertEqual(records[-1]["returncode"], 7)
            self.assertIsNotNone(records[-1]["error"])

    def test_timeout_output_is_bounded_and_marked_truncated(self) -> None:
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(subject, "CELL_TIMEOUT_SECONDS", 1.0):
            fixture = Fixture(Path(directory), "timeout")
            self.assertEqual(subject.run_campaign(fixture.config()), 1)
            self.assertEqual(fixture.invocations(), [0, 1, 2])
            self.assertEqual(
                (fixture.output / "cell-02.stdout").stat().st_size,
                subject.MAX_OUTPUT_BYTES)
            records = [json.loads(line) for line in
                       (fixture.output / "results.jsonl").read_text().splitlines()]
            self.assertIn("truncated oversized stdout", records[-1]["error"])


if __name__ == "__main__":
    unittest.main()
