#!/usr/bin/python3
"""Pure/scratch tests for :mod:`Wh2V2FacadeTimingLaunch`.

These tests never use the fixed campaign namespace, systemd, I2C, the thermal
sampler, or either scientific worker.
"""

from __future__ import annotations

import importlib.util
import ast
import copy
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


MODULE_PATH = Path(__file__).with_name("Wh2V2FacadeTimingLaunch.py")
SPEC = importlib.util.spec_from_file_location("wh2_facade_launcher", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("could not load facade launcher")
launch = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = launch
SPEC.loader.exec_module(launch)


def git(root: Path, *arguments: str) -> bytes:
    result = subprocess.run(
        ["/usr/bin/git", "-C", str(root), *arguments],
        env={
            "GIT_CONFIG_GLOBAL": "/dev/null",
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_NO_REPLACE_OBJECTS": "1",
            "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        },
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.decode("utf-8", "replace"))
    return result.stdout


def sampler_fixture() -> tuple[bytes, bytes, dict]:
    raw = (
        ",".join(launch.SAMPLER_COLUMNS) + "\n" +
        ",".join([
            "2026-08-29T00:00:00.000Z", "1.000000", "1.0", "3000.0",
            "50.0", *(["45.0"] * 8), "0", "0.1", "0.2", "0.3",
            "0", "0",
        ]) + "\n"
    ).encode("ascii")
    sensors = {
        name: {
            "attempt_errors": 0, "hot": False, "hot_streak": 0,
            "jump_c": None, "rate_c_per_s": None, "raw_c": 45.0,
            "reason": "ok", "valid": True,
        }
        for name in launch.SAMPLER_DIMM_FIELDS
    }
    header = {
        "expected_output_owner_uid": launch.CAMPAIGN_UID,
        "raw_columns": list(launch.SAMPLER_COLUMNS),
        "sampler_source_expected_sha256": "a" * 64,
        "sampling": launch.SAMPLER_SAMPLING,
        "schema": launch.SAMPLER_VALIDATION_STREAM_SCHEMA,
        "thresholds": launch.SAMPLER_THRESHOLDS,
    }
    sample = {
        "consecutive_fault_rows": 0, "decision": "continue",
        "edac_ce_delta": 0, "edac_ue_delta": 0, "fault_count": 0,
        "hot_sensors": [], "monotonic_s": 1.0, "read_error_count": 0,
        "sample_index": 0, "schema": launch.SAMPLER_VALIDATION_SCHEMA,
        "sensors": sensors,
    }
    validation = b"".join(
        launch.canonical_bytes(value) + b"\n" for value in (header, sample)
    )
    summary_sensors = {
        name: {
            "attempt_errors": 0, "invalid_samples": 0, "max_hot_streak": 0,
            "max_raw_c": 45.0, "max_valid_c": 45.0, "raw_samples": 1,
            "read_error_samples": 0, "valid_samples": 1,
        }
        for name in launch.SAMPLER_DIMM_FIELDS
    }
    receipt = {
        "cpu_tctl_max_c": 50.0,
        "finished_monotonic_ns": 2_000_000_000,
        "finished_utc": "2026-08-29T00:00:01.000Z",
        "sampler_source": {"expected_sha256": "a" * 64},
        "started_monotonic_ns": 0,
        "started_utc": "2026-08-28T23:59:59.000Z",
        "summary": {
            "consecutive_fault_rows": 0, "decision": "continue",
            "dimm_attempt_errors_total": 0,
            "dimm_invalid_samples_total": 0,
            "dimm_read_error_samples_total": 0,
            "edac_ce_baseline": 0, "edac_ce_delta": 0,
            "edac_ce_last": 0, "edac_ue_baseline": 0,
            "edac_ue_delta": 0, "edac_ue_last": 0,
            "max_consecutive_fault_rows": 0, "sample_count": 1,
            "sensors": summary_sensors,
        },
    }
    return raw, validation, receipt


class CanonicalAndParserTests(unittest.TestCase):
    def test_canonical_document_rejects_duplicate(self) -> None:
        with self.assertRaises(launch.LaunchError):
            launch.parse_canonical_document(b'{"a":1,"a":1}\n', "duplicate")

    def test_exact_public_modes(self) -> None:
        commit = "a" * 40
        execute = launch.parse_args([
            "--execute", "--expected-harness-commit", commit,
        ])
        self.assertTrue(execute.execute)
        with self.assertRaises(launch.LaunchError):
            launch.parse_args([
                "--execute", "--expected-harness-commit=" + commit,
            ])
        with self.assertRaises(launch.LaunchError):
            launch.parse_args(["--selftest", "--execute"])

    def test_systemd_deadline_and_containment_are_exact(self) -> None:
        argv = launch.systemd_run_argv("b" * 40)
        self.assertIn("--property=RuntimeMaxSec=1020s", argv)
        self.assertIn("--property=KillMode=control-group", argv)
        self.assertIn("--property=ExitType=cgroup", argv)
        self.assertIn("--property=AllowedCPUs=120-123", argv)
        self.assertGreater(
            launch.SERVICE_DEADLINE_SECONDS,
            launch.PRE_RELEASE_SECONDS + launch.EXTERNAL_DEADLINE_SECONDS
            + launch.POST_CONTROLLER_SECONDS
            + launch.SERVICE_ACTIVATION_STOP_MARGIN_SECONDS,
        )
        self.assertLess(
            launch.CONTROLLER_ADMISSION_SECONDS
            + launch.INTERNAL_DEADLINE_SECONDS
            + launch.POST_CONTROLLER_SECONDS,
            launch.EXTERNAL_DEADLINE_SECONDS,
        )

    def test_owned_fd_close_return_fault_is_resumable(self) -> None:
        class Owner:
            pass
        owner = Owner()
        owner.fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        target = owner.fd
        real_close = os.close
        fired = False
        def close_then_raise(fd: int) -> None:
            nonlocal fired
            real_close(fd)
            if fd == target and not fired:
                fired = True
                raise RuntimeError("injected owned close return fault")
        with mock.patch.object(launch.os, "close", side_effect=close_then_raise):
            launch.close_object_fd(owner, "fd", "test owned fd")
        self.assertTrue(fired)
        self.assertEqual(owner.fd, -1)
        launch.close_object_fd(owner, "fd", "test owned fd")

    def test_dict_fd_close_return_fault_is_resumable(self) -> None:
        owner = {"evidence": os.open(
            "/dev/null", os.O_RDONLY | os.O_CLOEXEC,
        )}
        target = owner["evidence"]
        real_close = os.close
        fired = False
        def close_then_raise(fd: int) -> None:
            nonlocal fired
            real_close(fd)
            if fd == target and not fired:
                fired = True
                raise RuntimeError("injected dict close return fault")
        with mock.patch.object(launch.os, "close", side_effect=close_then_raise):
            launch.close_dict_fd(owner, "evidence", "test dict fd")
        self.assertTrue(fired)
        self.assertNotIn("evidence", owner)
        launch.close_dict_fd(owner, "evidence", "test dict fd")

    def test_fd_roster_rejects_same_number_different_authority(self) -> None:
        fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        authority = launch.fd_authority_roster((fd,))
        os.close(fd)
        replacement = os.open("/dev/zero", os.O_RDONLY | os.O_CLOEXEC)
        try:
            self.assertEqual(replacement, fd)
            with self.assertRaises(launch.LaunchError):
                launch.verify_fd_authority_roster(
                    authority, "test replaced descriptor",
                )
        finally:
            os.close(replacement)


class GitGateTests(unittest.TestCase):
    def make_repo(self, root: Path) -> str:
        root.mkdir()
        git(root, "init", "--quiet")
        git(root, "config", "user.name", "Facade Test")
        git(root, "config", "user.email", "facade@example.invalid")
        (root / ".gitignore").write_text("ignored.bin\n", encoding="ascii")
        (root / "source.cpp").write_text("int value = 1;\n", encoding="ascii")
        executable = root / "script.sh"
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
        executable.chmod(0o755)
        git(root, "add", ".gitignore", "source.cpp", "script.sh")
        git(root, "commit", "--quiet", "-m", "fixture")
        commit = git(root, "rev-parse", "HEAD").decode("ascii").strip()
        git(root, "checkout", "--quiet", "--detach", commit)
        return commit

    def test_full_tree_gate_and_stat_cache_attack(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-git-") as temporary:
            root = Path(temporary) / "repo"
            commit = self.make_repo(root)
            receipt = launch.verify_git_tree(root, commit)
            self.assertEqual(receipt.commit, commit)
            self.assertEqual(len(receipt.entries), 3)
            source = root / "source.cpp"
            before = source.stat()
            original = source.read_bytes()
            forged = original.replace(b"1", b"2")
            self.assertEqual(len(original), len(forged))
            source.write_bytes(forged)
            os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns))
            with self.assertRaises(launch.LaunchError):
                launch.verify_git_tree(root, commit)

    def test_ignored_file_is_not_clean(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-git-") as temporary:
            root = Path(temporary) / "repo"
            commit = self.make_repo(root)
            (root / "ignored.bin").write_bytes(b"not an input\n")
            with self.assertRaises(launch.LaunchError):
                launch.verify_git_tree(root, commit)


class JournalTests(unittest.TestCase):
    def test_consumed_reservation_never_completes_without_attempt_prefix(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            def injected(point: str, name: str) -> None:
                if point == "artifact-file-fsync" and name == "ATTEMPT":
                    raise RuntimeError("persistent ATTEMPT fault")
            with self.assertRaises(launch.AttemptConsumedError) as caught:
                launch.AttemptJournal.reserve(
                    path, {"attempt": 1}, expected_uid=os.getuid(),
                    expected_gid=os.getgid(), test_parent=True,
                    fault_injector=injected,
                )
            journal = caught.exception.journal
            try:
                self.assertEqual(journal.names, [])
                journal.write_json("terminal.json", {"outcome": "invalid"})
                with self.assertRaises(launch.LaunchError):
                    journal.finish()
                self.assertFalse((path / "COMPLETE").exists())
            finally:
                journal.close()

    def test_consumed_reservation_recovers_one_shot_attempt_publication(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            fired = False
            def injected(point: str, name: str) -> None:
                nonlocal fired
                if (
                    point == "artifact-file-fsync" and name == "ATTEMPT"
                    and not fired
                ):
                    fired = True
                    raise RuntimeError("one-shot ATTEMPT fault")
            with self.assertRaises(launch.AttemptConsumedError) as caught:
                launch.AttemptJournal.reserve(
                    path, {"attempt": 1}, expected_uid=os.getuid(),
                    expected_gid=os.getgid(), test_parent=True,
                    fault_injector=injected,
                )
            journal = caught.exception.journal
            try:
                self.assertEqual(journal.names, ["ATTEMPT"])
                journal.write_json("terminal.json", {"outcome": "invalid"})
                journal.finish()
                self.assertTrue(journal.complete)
            finally:
                journal.close()
                path.chmod(0o700)

    def test_complete_is_last_and_directory_is_sealed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            root = Path(temporary)
            path = root / "attempt"
            journal = launch.AttemptJournal.reserve(
                path,
                launch.self_hashed("attempt.test.v1", {"value": 1}),
                expected_uid=os.getuid(), expected_gid=os.getgid(),
                test_parent=True,
            )
            try:
                journal.write_json("evidence.json", {"value": 2})
                journal.write_json("terminal.json", {"outcome": "invalid"})
                journal.finish()
                self.assertTrue(journal.complete)
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o500)
                self.assertEqual(
                    stat.S_IMODE((path / "COMPLETE").stat().st_mode), 0,
                )
                self.assertEqual(journal.names[-1], "COMPLETE")
            finally:
                journal.close()
                path.chmod(0o700)

    def test_complete_requires_terminal_last(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            try:
                with self.assertRaises(launch.LaunchError):
                    journal.finish()
            finally:
                journal.close()

    def test_unsealed_reserved_stream_forbids_complete(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            stream_fd = -1
            try:
                stream_fd = journal.open_stream("guardian.jsonl")
                journal.write_json("terminal.json", {"outcome": "invalid"})
                with self.assertRaises(launch.LaunchError):
                    journal.finish()
                self.assertFalse(journal.complete)
                self.assertEqual(
                    stat.S_IMODE((path / "guardian.jsonl").stat().st_mode),
                    0o600,
                )
            finally:
                journal.close()

    def test_failed_artifact_publication_rolls_back_exact_namespace(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            def injected(point: str, name: str) -> None:
                if point == "artifact-file-fsync" and name == "evidence.json":
                    raise RuntimeError("injected post-fsync failure")
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
                fault_injector=injected,
            )
            try:
                with self.assertRaises(RuntimeError):
                    journal.write_json("evidence.json", {"value": 2})
                self.assertNotIn("evidence.json", journal.names)
                self.assertFalse((path / "evidence.json").exists())
                self.assertEqual(sorted(os.listdir(journal.directory_fd)), ["ATTEMPT"])
            finally:
                journal.close()

    def test_stream_return_handoff_remains_journal_owned(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            try:
                returned = journal.open_stream("guardian.jsonl")
                self.assertEqual(journal.stream_fds["guardian.jsonl"], returned)
                os.write(returned, b"{}\n")
                journal.seal_stream("guardian.jsonl", returned, 4096)
            finally:
                journal.close()

    def test_empty_stream_discard_recovers_unlink_return_fault(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            real_unlink = os.unlink
            fired = False
            def unlink_then_raise(target: object, *args: object,
                                  **kwargs: object) -> None:
                nonlocal fired
                if target == "guardian.jsonl" and not fired:
                    fired = True
                    real_unlink(target, *args, **kwargs)
                    raise RuntimeError("injected unlink return fault")
                real_unlink(target, *args, **kwargs)
            try:
                journal.open_stream("guardian.jsonl")
                with mock.patch.object(
                    launch.os, "unlink", side_effect=unlink_then_raise,
                ):
                    journal.discard_empty_stream("guardian.jsonl")
                self.assertTrue(fired)
                self.assertNotIn("guardian.jsonl", journal.names)
                self.assertNotIn("guardian.jsonl", journal.stream_fds)
                self.assertFalse((path / "guardian.jsonl").exists())
            finally:
                journal.close()

    def test_resumable_artifact_recovers_exact_fixed_publication(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            try:
                payload = b"sealed evidence\n"
                first = journal.write_bytes("evidence.bin", payload)
                second = journal.write_bytes_resumable("evidence.bin", payload)
                self.assertEqual(first, second)
                with self.assertRaises(launch.LaunchError):
                    journal.write_bytes_resumable("evidence.bin", b"forged\n")
            finally:
                journal.close()

    def test_complete_private_failure_never_publishes_fixed_marker(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            def injected(point: str, name: str) -> None:
                if point == "complete-file-fsync":
                    raise RuntimeError("injected private COMPLETE failure")
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
                fault_injector=injected,
            )
            try:
                journal.write_json("terminal.json", {"outcome": "invalid"})
                with self.assertRaises(RuntimeError):
                    journal.finish()
                self.assertFalse((path / "COMPLETE").exists())
                self.assertFalse(any(name.startswith(".complete.")
                                     for name in os.listdir(path)))
                self.assertFalse(journal.complete)
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o700)
            finally:
                journal.close()

    def test_complete_open_create_then_raise_recovers_name_and_fd(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            try:
                journal.write_json("terminal.json", {"outcome": "invalid"})
                before = launch.open_fd_roster()
                real_open = os.open
                fired = False
                def create_then_raise(target: object, flags: int,
                                      *args: object, **kwargs: object) -> int:
                    nonlocal fired
                    if not fired and str(target).startswith(".complete."):
                        fired = True
                        real_open(target, flags, *args, **kwargs)
                        raise RuntimeError("injected create-then-raise")
                    return real_open(target, flags, *args, **kwargs)
                with mock.patch.object(launch.os, "open", side_effect=create_then_raise):
                    with self.assertRaises(RuntimeError):
                        journal.finish()
                self.assertTrue(fired)
                self.assertEqual(launch.open_fd_roster(), before)
                self.assertFalse((path / "COMPLETE").exists())
                self.assertFalse(any(name.startswith(".complete.")
                                     for name in os.listdir(path)))
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o700)
            finally:
                journal.close()

    def test_complete_retry_after_durability_fault_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-journal-") as temporary:
            path = Path(temporary) / "attempt"
            fired = False
            def injected(point: str, name: str) -> None:
                nonlocal fired
                if point == "complete-file-fsync" and not fired:
                    fired = True
                    raise RuntimeError("one-shot COMPLETE fault")
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
                fault_injector=injected,
            )
            try:
                journal.write_json("terminal.json", {"outcome": "invalid"})
                with self.assertRaises(RuntimeError):
                    journal.finish()
                journal.finish()
                self.assertTrue(journal.complete)
                self.assertEqual(
                    sorted(os.listdir(path)),
                    sorted(["ATTEMPT", "terminal.json", "COMPLETE"]),
                )
            finally:
                journal.close()
                path.chmod(0o700)

    def test_child_bundle_copy_resumes_to_exact_six_artifacts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-copy-fi-") as temporary:
            path = Path(temporary) / "attempt"
            fired = False
            def injected(point: str, name: str) -> None:
                nonlocal fired
                if (
                    point == "artifact-directory-fsync"
                    and name == "screen.summary.json" and not fired
                ):
                    fired = True
                    raise RuntimeError("injected child-copy fault")
            journal = launch.AttemptJournal.reserve(
                path, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
                fault_injector=injected,
            )
            state: dict[str, dict] = {}
            data = {
                name: (name + "\n").encode("ascii")
                for name in launch.CONTROLLER_OUTPUT_NAMES
            }
            bundle = launch.ChildBundle(-1, {}, data, {}, {}, None, None)
            try:
                with self.assertRaises(RuntimeError):
                    launch.copy_child_bundle(bundle, journal, state)
                self.assertNotEqual(
                    set(state), set(launch.CONTROLLER_OUTPUT_NAMES),
                )
                result = launch.copy_child_bundle(bundle, journal, state)
                self.assertIs(result, state)
                self.assertEqual(set(state), set(launch.CONTROLLER_OUTPUT_NAMES))
                for name, root_name in (
                    ("raw.jsonl", "screen.raw.jsonl"),
                    ("summary.json", "screen.summary.json"),
                    ("provenance.json", "screen.provenance.json"),
                    ("current.stderr", "screen.current.stderr"),
                    ("parent.stderr", "screen.parent.stderr"),
                    ("COMPLETE", "screen.child-COMPLETE.json"),
                ):
                    self.assertEqual((path / root_name).read_bytes(), data[name])
            finally:
                journal.close()


class ChildBundleOwnershipTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[launch.BoundDirectory, int]:
        parent_path = root / "controller-parent"
        child_path = parent_path / launch.CONTROLLER_OUTPUT.name
        parent_path.mkdir(mode=0o700)
        child_path.mkdir(mode=0o700)
        values = {
            "raw.jsonl": launch.canonical_bytes({
                "schema":
                    "wirehair.wh2.v2-facade-timing-screen.emergency-invalid.v1",
            }) + b"\n",
            "summary.json": launch.canonical_bytes({
                "schema":
                    "wirehair.wh2.v2-facade-timing-screen.emergency-invalid.v1",
            }) + b"\n",
            "provenance.json": launch.canonical_bytes({
                "schema":
                    "wirehair.wh2.v2-facade-timing-screen.emergency-invalid.v1",
            }) + b"\n",
            "current.stderr": b"",
            "parent.stderr": b"",
        }
        complete = {
            "campaign": launch.CAMPAIGN,
            "files": [
                {
                    "bytes": len(values[name]), "name": name,
                    "sha256": launch.sha256_bytes(values[name]),
                }
                for name in launch.CONTROLLER_OUTPUT_NAMES[:-1]
            ],
            "outcome": "invalid", "preimage_sha256": None,
            "schema": launch.CHILD_COMPLETE_SCHEMA,
        }
        complete["preimage_sha256"] = launch.sha256_bytes(
            launch.canonical_bytes(complete),
        )
        values["COMPLETE"] = launch.canonical_bytes(complete) + b"\n"
        for name, raw in values.items():
            target = child_path / name
            target.write_bytes(raw)
            target.chmod(0o400)
        child_path.chmod(0o500)
        parent_fd = os.open(
            parent_path, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
        )
        info = os.fstat(parent_fd)
        return launch.BoundDirectory(
            parent_path, os.open(
                root, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
            ), os.getuid(), os.getgid(), fd=parent_fd,
            binding=(info.st_dev, info.st_ino), state="bound",
        ), parent_fd

    def test_child_bundle_return_fault_remains_owned_and_closable(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-child-owner-") as temporary:
            parent, _ = self._fixture(Path(temporary))
            owner = launch.OwnershipSlot("test child bundle")
            baseline = launch.open_fd_roster()
            real_adopt = owner.adopt
            def adopt_then_raise(value: object) -> None:
                real_adopt(value)
                raise RuntimeError("injected bundle adopt return fault")
            try:
                with mock.patch.object(owner, "adopt", side_effect=adopt_then_raise):
                    with self.assertRaises(RuntimeError):
                        launch.open_child_bundle(parent, owner=owner)
                self.assertIsInstance(owner.value, launch.ChildBundle)
                self.assertEqual(owner.value.state, "open")
                owner.value.close()
                self.assertEqual(owner.value.state, "closed")
                owner.value = None
                self.assertEqual(launch.open_fd_roster(), baseline)
            finally:
                if owner.value is not None:
                    owner.value.close()
                parent.close()
                launch.close_fd_once(parent.parent_fd, "test parent authority")

    def test_child_file_open_return_fault_closes_exact_fd_delta(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-child-open-") as temporary:
            parent, _ = self._fixture(Path(temporary))
            owner = launch.OwnershipSlot("test child bundle")
            baseline = launch.open_fd_roster()
            real_open = os.open
            fired = False
            def open_then_raise(target: object, flags: int,
                                *args: object, **kwargs: object) -> int:
                nonlocal fired
                if target == "summary.json" and not fired:
                    fired = True
                    real_open(target, flags, *args, **kwargs)
                    raise RuntimeError("injected child open return fault")
                return real_open(target, flags, *args, **kwargs)
            try:
                with mock.patch.object(
                    launch.os, "open", side_effect=open_then_raise,
                ):
                    with self.assertRaises(RuntimeError):
                        launch.open_child_bundle(parent, owner=owner)
                self.assertTrue(fired)
                self.assertIsNone(owner.value)
                self.assertEqual(launch.open_fd_roster(), baseline)
            finally:
                parent.close()
                launch.close_fd_once(parent.parent_fd, "test parent authority")

    def test_child_bundle_close_return_fault_is_idempotent(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-child-close-") as temporary:
            parent, _ = self._fixture(Path(temporary))
            owner = launch.OwnershipSlot("test child bundle")
            baseline = launch.open_fd_roster()
            launch.open_child_bundle(parent, owner=owner)
            bundle = owner.value
            target_fd = bundle.files["summary.json"]
            real_close = os.close
            fired = False
            def close_then_raise(fd: int) -> None:
                nonlocal fired
                real_close(fd)
                if fd == target_fd and not fired:
                    fired = True
                    raise RuntimeError("injected child close return fault")
            try:
                with mock.patch.object(
                    launch.os, "close", side_effect=close_then_raise,
                ):
                    bundle.close()
                self.assertTrue(fired)
                self.assertEqual(bundle.state, "closed")
                bundle.close()
                owner.value = None
                self.assertEqual(launch.open_fd_roster(), baseline)
            finally:
                if owner.value is not None:
                    owner.value.close()
                parent.close()
                launch.close_fd_once(parent.parent_fd, "test parent authority")

    def test_namespace_seal_recovers_closed_bundle_slot_handoff(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-child-seal-") as temporary:
            root = Path(temporary)
            parent, _ = self._fixture(root)
            owner = launch.OwnershipSlot("test child bundle")
            launch.open_child_bundle(parent, owner=owner)
            owner.value.close()
            self.assertEqual(owner.value.state, "closed")
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session.var_tmp_fd = parent.parent_fd
            session.controller_directory = parent
            session._bundle_owner = owner
            session.child_bundle_validated = True
            session.child_copies_complete = True
            try:
                session.seal_namespaces()
                self.assertTrue(session.namespaces_sealed)
                self.assertIsNone(session.bundle)
                self.assertEqual(
                    stat.S_IMODE(os.fstat(parent.fd).st_mode), 0o500,
                )
            finally:
                os.chmod(parent.path, 0o700)
                os.chmod(parent.path / launch.CONTROLLER_OUTPUT.name, 0o700)
                parent.close()
                launch.close_fd_once(parent.parent_fd, "test parent authority")


class GuardianTests(unittest.TestCase):
    def test_guardian_append_recovers_every_short_write(self) -> None:
        armed = {
            "armed_monotonic_ns": 10, "deadline_monotonic_ns": 20,
            "event": "armed", "pid": 7,
            "schema": "wirehair.test.guardian.v1",
        }
        terminal = {
            "deadline_monotonic_ns": 20, "event": "completed",
            "observed_monotonic_ns": 15, "pid": 7,
            "schema": "wirehair.test.guardian.v1",
        }
        armed_raw = launch.canonical_bytes(armed) + b"\n"
        terminal_raw = launch.canonical_bytes(terminal) + b"\n"
        with tempfile.TemporaryDirectory(prefix="wh2-guardian-short-") as temporary:
            directory_fd = os.open(
                temporary, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
            )
            try:
                for record_name, record, prefix, raw in (
                    ("armed", armed, b"", armed_raw),
                    ("terminal", terminal, armed_raw, terminal_raw),
                ):
                    for short_count in range(1, len(raw)):
                        name = "{}-{}".format(record_name, short_count)
                        event_fd = os.open(
                            name, os.O_RDWR | os.O_CREAT | os.O_EXCL
                            | os.O_CLOEXEC, 0o600, dir_fd=directory_fd,
                        )
                        real_write = os.write
                        fired = False
                        def short_once(fd: int, data: bytes) -> int:
                            nonlocal fired
                            if fd == event_fd and not fired:
                                fired = True
                                return real_write(fd, data[:short_count])
                            return real_write(fd, data)
                        try:
                            if prefix:
                                self.assertEqual(
                                    os.write(event_fd, prefix), len(prefix),
                                )
                            with mock.patch.object(
                                launch.os, "write", side_effect=short_once,
                            ):
                                launch.append_guardian_record(
                                    event_fd, directory_fd, record,
                                )
                            self.assertTrue(fired)
                            expected = prefix + raw
                            self.assertEqual(
                                os.pread(event_fd, len(expected), 0), expected,
                            )
                            self.assertEqual(
                                os.fstat(event_fd).st_size, len(expected),
                            )
                        finally:
                            os.close(event_fd)
            finally:
                os.close(directory_fd)

    def test_guardian_append_recovers_write_and_fsync_return_faults(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-guardian-append-") as temporary:
            directory_fd = os.open(
                temporary, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
            )
            event_fd = os.open(
                "events", os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o600, dir_fd=directory_fd,
            )
            record = {
                "event": "armed", "pid": 7,
                "schema": "wirehair.test.guardian.v1",
            }
            expected = launch.canonical_bytes(record) + b"\n"
            real_write = os.write
            real_fsync = os.fsync
            write_fired = False
            fsync_fired = False
            def write_then_raise(fd: int, data: bytes) -> int:
                nonlocal write_fired
                result = real_write(fd, data)
                if fd == event_fd and not write_fired:
                    write_fired = True
                    raise RuntimeError("injected write return fault")
                return result
            def fsync_then_raise(fd: int) -> None:
                nonlocal fsync_fired
                real_fsync(fd)
                if fd == event_fd and not fsync_fired:
                    fsync_fired = True
                    raise RuntimeError("injected fsync return fault")
            try:
                with mock.patch.object(
                    launch.os, "write", side_effect=write_then_raise,
                ), mock.patch.object(
                    launch.os, "fsync", side_effect=fsync_then_raise,
                ):
                    launch.append_guardian_record(
                        event_fd, directory_fd, record,
                    )
                self.assertTrue(write_fired)
                self.assertTrue(fsync_fired)
                self.assertEqual(os.pread(event_fd, len(expected), 0), expected)
                self.assertEqual(os.fstat(event_fd).st_size, len(expected))
            finally:
                os.close(event_fd)
                os.close(directory_fd)

    def _events_fd(self, events: list[dict]) -> tuple[object, int]:
        handle = tempfile.NamedTemporaryFile(
            prefix="wh2-facade-guardian-events-", delete=True,
        )
        handle.write(b"".join(
            launch.canonical_bytes(value) + b"\n" for value in events
        ))
        handle.flush()
        return handle, handle.fileno()

    def test_guardian_event_enum_keys_chronology_and_exit(self) -> None:
        schema = "wirehair.wh2.v2-facade-timing-root-guardian-event.v1"
        armed = {
            "armed_monotonic_ns": 10, "deadline_monotonic_ns": 100,
            "event": "armed", "pid": 7, "schema": schema,
        }
        completed = {
            "deadline_monotonic_ns": 100, "event": "completed",
            "observed_monotonic_ns": 99, "pid": 7, "schema": schema,
        }
        handle, fd = self._events_fd([armed, completed])
        try:
            launch.validate_guardian_events(fd, 7, 100, "completed", 0)
        finally:
            handle.close()
        for mutation in (
            {**completed, "event": "bogus-terminal"},
            {**completed, "extra": 1},
            {**completed, "observed_monotonic_ns": 100},
        ):
            handle, fd = self._events_fd([armed, mutation])
            try:
                with self.assertRaises(launch.LaunchError):
                    launch.validate_guardian_events(fd, 7, 100, None, 0)
            finally:
                handle.close()

    def test_resolved_guardian_receipt_is_idempotent(self) -> None:
        receipt = {"completed": True, "returncode": 0}
        handle = launch.GuardianHandle(
            process=launch.ProcessHandle(
                role="guardian", pid=os.getpid(), pidfd=-1,
                start_ticks=1, argv=("guardian",), stdout_fd=-1,
                stderr_fd=-1, returncode=0,
            ),
            control_fd=-1, event_fd=-1, ready={}, completed=True,
            completion_may_have_been_sent=True, terminal_receipt=receipt,
        )
        self.assertIs(handle.complete(), receipt)


class RuntimePrimitiveTests(unittest.TestCase):
    def test_parent_containment_resumes_after_kill_and_copy_fault(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-containment-fi-") as temporary:
            attempt = Path(temporary) / "attempt"
            fired = False
            def injected(point: str, name: str) -> None:
                nonlocal fired
                if (
                    point == "artifact-file-fsync"
                    and name == "containment-000-result.json"
                    and not fired
                ):
                    fired = True
                    raise RuntimeError("injected result publication fault")
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
                fault_injector=injected,
            )
            kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            root = Path(temporary) / "cgroup"
            session._layout_owner.adopt(launch.CgroupLayout(
                service=root, supervisor=root / "supervisor",
                run=root / "run", sampler=root / "run/sampler",
                experiment=root / "run/experiment",
                verifier=root / "run/verifier", experiment_kill_fd=kill_w,
                sampler_kill_fd=-1, verifier_kill_fd=-1,
                run_kill_fd=-1,
            ))
            sequence = iter(([42], [42], [42], [], []))
            try:
                with mock.patch.object(
                    launch, "cgroup_descendant_pids",
                    side_effect=lambda _path: list(next(sequence)),
                ):
                    with self.assertRaises(RuntimeError):
                        session.contain_cgroup(
                            session.layout.experiment, kill_w, "test-fault",
                        )
                self.assertTrue(fired)
                self.assertEqual(os.read(kill_r, 1), b"1")
                self.assertEqual(len(session.containment_states), 1)
                self.assertIsNone(
                    session.containment_states[0].result_receipt,
                )
                with mock.patch.object(
                    launch, "cgroup_descendant_pids", return_value=[],
                ):
                    roster = session.contain_cgroup(
                        session.layout.experiment, kill_w, "test-fault",
                    )
                self.assertEqual(roster, [42])
                self.assertEqual(len(session.containment_receipts), 1)
                self.assertIsNotNone(
                    session.containment_states[0].result_receipt,
                )
                self.assertTrue(
                    (attempt / "containment-000-intent.json").is_file()
                )
                self.assertTrue(
                    (attempt / "containment-000-result.json").is_file()
                )
            finally:
                journal.close()
                os.close(kill_r)
                os.close(kill_w)

    def test_completed_containment_reconstructs_aggregate_roster(self) -> None:
        session = launch.ExecutionSession(
            launch.LaunchConfig(expected_harness_commit="a" * 40),
        )
        root = Path("/nonexistent/test-cgroup")
        session._layout_owner.adopt(launch.CgroupLayout(
            service=root, supervisor=root / "supervisor",
            run=root / "run", sampler=root / "run/sampler",
            experiment=root / "run/experiment",
            verifier=root / "run/verifier", experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1, run_kill_fd=-1,
        ))
        state = launch.ContainmentState(
            index=0, path=str(session.layout.experiment), reason="test",
            roster=[41, 42], requested_monotonic_ns=1,
            intent_receipt={"sha256": "a" * 64},
            result_receipt={"sha256": "b" * 64},
            completed_monotonic_ns=2,
        )
        session.containment_states.append(state)
        session.containment_receipts.append({
            "intent": state.intent_receipt, "path": state.path,
            "reason": state.reason, "result": state.result_receipt,
            "initial_roster": list(state.roster),
            "observed_roster": list(state.observed_roster),
        })
        session.resume_containments()
        self.assertEqual(session.containment_roster, [41, 42])

    def test_containment_atomic_kill_accepts_roster_growth(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-containment-growth-") as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            root = Path(temporary) / "cgroup"
            session._layout_owner.adopt(launch.CgroupLayout(
                service=root, supervisor=root / "supervisor",
                run=root / "run", sampler=root / "run/sampler",
                experiment=root / "run/experiment",
                verifier=root / "run/verifier", experiment_kill_fd=kill_w,
                sampler_kill_fd=-1, verifier_kill_fd=-1,
                run_kill_fd=-1,
            ))
            sequence = iter(([42], [42, 43], [42, 43], [], []))
            try:
                with mock.patch.object(
                    launch, "cgroup_descendant_pids",
                    side_effect=lambda _path: list(next(sequence)),
                ):
                    roster = session.contain_cgroup(
                        session.layout.experiment, kill_w, "fork-race",
                    )
                self.assertEqual(roster, [42, 43])
                self.assertEqual(os.read(kill_r, 1), b"1")
                state = session.containment_states[0]
                self.assertEqual(state.roster, [42])
                self.assertEqual(state.observed_roster, [42, 43])
                result = json.loads(
                    (attempt / "containment-000-result.json").read_text(
                        encoding="utf-8",
                    )
                )
                self.assertEqual(result["initial_roster"], [42])
                self.assertEqual(result["observed_roster"], [42, 43])
                self.assertEqual(
                    result["kill_scope"],
                    "all-processes-in-exact-cgroup-and-descendants-at-"
                    "atomic-cgroup.kill",
                )
            finally:
                journal.close()
                os.close(kill_r)
                os.close(kill_w)

    def test_containment_kill_return_fault_retains_inner_roster(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-containment-return-") as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            root = Path(temporary) / "cgroup"
            session._layout_owner.adopt(launch.CgroupLayout(
                service=root, supervisor=root / "supervisor",
                run=root / "run", sampler=root / "run/sampler",
                experiment=root / "run/experiment",
                verifier=root / "run/verifier", experiment_kill_fd=kill_w,
                sampler_kill_fd=-1, verifier_kill_fd=-1,
                run_kill_fd=-1,
            ))
            sequence = iter(([42], [42], [42, 43], [], [], []))
            real_kill = launch.kill_cgroup_and_wait
            fired = False
            def returned_then_raised(*arguments: object,
                                     **keywords: object) -> object:
                nonlocal fired
                result = real_kill(*arguments, **keywords)
                if not fired:
                    fired = True
                    raise RuntimeError("injected kill return-handoff fault")
                return result
            try:
                with mock.patch.object(
                    launch, "cgroup_descendant_pids",
                    side_effect=lambda _path: list(next(sequence)),
                ), mock.patch.object(
                    launch, "kill_cgroup_and_wait",
                    side_effect=returned_then_raised,
                ):
                    with self.assertRaises(RuntimeError):
                        session.contain_cgroup(
                            session.layout.experiment, kill_w, "return-fault",
                        )
                    self.assertEqual(
                        session.containment_states[0].observed_roster,
                        [42, 43],
                    )
                    roster = session.contain_cgroup(
                        session.layout.experiment, kill_w, "return-fault",
                    )
                self.assertTrue(fired)
                self.assertEqual(roster, [42, 43])
                self.assertEqual(os.read(kill_r, 1), b"1")
                result = json.loads(
                    (attempt / "containment-000-result.json").read_text(
                        encoding="utf-8",
                    )
                )
                self.assertEqual(result["observed_roster"], [42, 43])
            finally:
                journal.close()
                os.close(kill_r)
                os.close(kill_w)

    def test_final_run_containment_invalidates_scientific_outcome(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-final-quiescence-") as temporary:
            root = Path(temporary) / "cgroup"
            run = root / "run"
            run.mkdir(parents=True)
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._layout_owner.adopt(launch.CgroupLayout(
                service=root, supervisor=root / "supervisor", run=run,
                sampler=run / "sampler", experiment=run / "experiment",
                verifier=run / "verifier", experiment_kill_fd=-1,
                sampler_kill_fd=-1, verifier_kill_fd=-1, run_kill_fd=7,
            ))
            session.child_outcome = "pass"
            sequence = iter(([42], []))
            with mock.patch.object(
                launch, "cgroup_descendant_pids",
                side_effect=lambda _path: list(next(sequence)),
            ), mock.patch.object(
                session, "contain_cgroup", return_value=[42],
            ), mock.patch.object(session, "_wait_tracked_after_kill"):
                session.quiesce()
            self.assertTrue(session.quiesced)
            self.assertTrue(session.unexpected_final_containment)
            self.assertIn(
                "final run containment: LaunchError: unexpected run descendants required "
                "cgroup.kill",
                session.errors,
            )
            self.assertEqual(session.containment_roster, [42])

    def test_adopted_reap_return_fault_retains_pid_status(self) -> None:
        child = os.fork()
        if child == 0:
            os._exit(7)
        try:
            deadline = time.monotonic() + 2.0
            while True:
                observed = os.waitid(
                    os.P_PID, child,
                    os.WEXITED | os.WNOHANG | os.WNOWAIT,
                )
                if observed is not None:
                    break
                if time.monotonic() >= deadline:
                    self.fail("test child did not become waitable")
                time.sleep(0.005)
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            real_reap = launch.reap_adopted_children
            fired = False
            def returned_then_raised(values: list[dict[str, int]]) -> object:
                nonlocal fired
                result = real_reap(values)
                if not fired:
                    fired = True
                    raise RuntimeError("injected adopted-reap return fault")
                return result
            with mock.patch.object(
                launch, "reap_adopted_children",
                side_effect=returned_then_raised,
            ):
                with self.assertRaises(launch.LaunchError):
                    session.close_runtime_authority()
                self.assertEqual(session.adopted_reaps, [{
                    "pid": child, "returncode": 7,
                    "waitid_code": 1, "waitid_status": 7,
                }])
                self.assertIn(child, session.containment_roster)
                session.close_runtime_authority()
            self.assertTrue(fired)
            self.assertTrue(session.runtime_closed)
            with self.assertRaises(ChildProcessError):
                os.waitpid(child, os.WNOHANG)
        finally:
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                os.waitpid(child, 0)
            except ChildProcessError:
                pass

    def test_containment_finishes_after_guardian_deadline_terminal(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-containment-deadline-") as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            root = Path(temporary) / "cgroup"
            session._layout_owner.adopt(launch.CgroupLayout(
                service=root, supervisor=root / "supervisor",
                run=root / "run", sampler=root / "run/sampler",
                experiment=root / "run/experiment",
                verifier=root / "run/verifier", experiment_kill_fd=-1,
                sampler_kill_fd=-1, verifier_kill_fd=-1, run_kill_fd=-1,
            ))
            class TerminalGuardian:
                process = type("P", (), {"pid": 77})()
                def poll(self) -> None:
                    raise launch.LaunchError("deadline terminal")
            session._guardian_owner.adopt(TerminalGuardian())
            session._root_window = (1, 2)
            experiment_calls = 0
            def roster(path: Path) -> list[int]:
                nonlocal experiment_calls
                if path == session.layout.experiment:
                    experiment_calls += 1
                    return [42] if experiment_calls == 1 else []
                return []
            def resolve() -> None:
                session.guardian_authenticated = True
                session.guardian_receipt = {"event": "deadline-fired"}
            try:
                with mock.patch.object(
                    launch, "cgroup_descendant_pids", side_effect=roster,
                ), mock.patch.object(session, "resolve_guardian", side_effect=resolve):
                    result = session.contain_cgroup(
                        session.layout.experiment, -1, "deadline-race",
                    )
                self.assertEqual(result, [42])
                self.assertTrue(session.guardian_authenticated)
                self.assertIsNotNone(
                    session.containment_states[0].result_receipt,
                )
            finally:
                journal.close()

    def test_blocked_spawn_post_fork_fault_reaps_and_closes(self) -> None:
        before = launch.open_fd_roster()
        observed: list[str] = []
        def injected(point: str, role: str) -> None:
            observed.append(point + ":" + role)
            if point == "preamble-post-fork":
                raise RuntimeError("injected post-fork handoff fault")
        with mock.patch.object(launch.CgroupTree, "move_pid"):
            with self.assertRaises(RuntimeError):
                launch.spawn_blocked_role(
                    "fake", Path("/fake/cgroup"), (123,),
                    ("/bin/true",), sampler=False,
                    python_info=os.stat("/usr/bin/python3.12"),
                    fault_injector=injected,
                )
        self.assertIn("preamble-post-fork:fake", observed)
        self.assertEqual(launch.open_fd_roster(), before)

    def test_blocked_fork_return_fault_recovers_unassigned_child(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        real_fork = os.fork
        fired = False
        def fork_then_raise() -> int:
            nonlocal fired
            pid = real_fork()
            if pid > 0 and not fired:
                fired = True
                raise RuntimeError("injected blocked fork return fault")
            return pid
        with mock.patch.object(launch.CgroupTree, "move_pid"), mock.patch.object(
            launch.os, "fork", side_effect=fork_then_raise,
        ):
            with self.assertRaises(RuntimeError):
                launch.fork_blocked_preamble(
                    "fake", Path("/fake/cgroup"), (123,),
                    ("/bin/true",), sampler=False,
                )
        self.assertTrue(fired)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_blocked_child_fork_return_fault_never_escapes_to_caller(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        marker_r, marker_w = os.pipe2(os.O_CLOEXEC)
        parent_pid = os.getpid()
        real_fork = os.fork
        def fork_then_raise_in_child() -> int:
            pid = real_fork()
            if pid == 0:
                raise RuntimeError("injected child-side fork return fault")
            return pid
        def quiet_child_failure(_message: str) -> None:
            os._exit(125)
        preamble = None
        try:
            with mock.patch.object(
                launch.os, "fork", side_effect=fork_then_raise_in_child,
            ), mock.patch.object(
                launch, "child_exec_failure", new=quiet_child_failure,
            ):
                try:
                    preamble = launch.fork_blocked_preamble(
                        "fake", Path("/fake/cgroup"), (123,),
                        ("/bin/true",), sampler=False,
                    )
                except BaseException:
                    if os.getpid() != parent_pid:
                        os.write(marker_w, b"escaped")
                        os._exit(99)
                    raise
            self.assertIsInstance(preamble, launch.BlockedSpawnPreamble)
            preamble.cleanup()
        finally:
            os.close(marker_w)
        try:
            self.assertEqual(os.read(marker_r, 64), b"")
        finally:
            os.close(marker_r)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_blocked_spawn_return_then_raise_uses_adopted_preamble(self) -> None:
        before = launch.open_fd_roster()
        real = launch.fork_blocked_preamble
        def returned_then_raised(*arguments: object,
                                 **keywords: object) -> object:
            real(*arguments, **keywords)
            raise RuntimeError("injected preamble return fault")
        with mock.patch.object(launch.CgroupTree, "move_pid"), mock.patch.object(
            launch, "fork_blocked_preamble", side_effect=returned_then_raised,
        ):
            with self.assertRaises(RuntimeError):
                launch.spawn_blocked_role(
                    "fake", Path("/fake/cgroup"), (123,),
                    ("/bin/true",), sampler=False,
                    python_info=os.stat("/usr/bin/python3.12"),
                )
        self.assertEqual(launch.open_fd_roster(), before)

    @unittest.skipUnless(hasattr(os, "pidfd_open"), "pidfd_open is unavailable")
    def test_blocked_spawn_pidfd_return_fault_closes_reused_number(self) -> None:
        before = launch.open_fd_roster()
        real_pidfd_open = os.pidfd_open
        fired = False
        def pidfd_then_raise(pid: int, flags: int = 0) -> int:
            nonlocal fired
            fd = real_pidfd_open(pid, flags)
            if not fired:
                fired = True
                raise RuntimeError("injected blocked pidfd return fault")
            return fd
        with mock.patch.object(launch.CgroupTree, "move_pid"), mock.patch.object(
            launch.os, "pidfd_open", side_effect=pidfd_then_raise,
        ):
            with self.assertRaises(RuntimeError):
                launch.spawn_blocked_role(
                    "fake", Path("/fake/cgroup"), (123,),
                    ("/bin/true",), sampler=False,
                    python_info=os.stat("/usr/bin/python3.12"),
                )
        self.assertTrue(fired)
        self.assertEqual(launch.open_fd_roster(), before)

    def test_guardian_spawn_post_fork_fault_reaps_and_closes(self) -> None:
        before = launch.open_fd_roster()
        event_fd = os.open("/dev/null", os.O_WRONLY | os.O_CLOEXEC)
        directory_fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
        layout = launch.CgroupLayout(
            service=Path("/fake"), supervisor=Path("/fake/supervisor"),
            run=Path("/fake/run"), sampler=Path("/fake/sampler"),
            experiment=Path("/fake/experiment"),
            verifier=Path("/fake/verifier"), experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1,
            run_kill_fd=kill_w,
        )
        def fake_guardian(control_fd: int, *_arguments: object) -> None:
            os.read(control_fd, 1)
            os._exit(0)
        def injected(point: str, _role: str) -> None:
            if point == "guardian-preamble-post-fork":
                raise RuntimeError("injected guardian post-fork fault")
        try:
            with mock.patch.object(
                launch, "guardian_child", new=fake_guardian,
            ):
                with self.assertRaises(RuntimeError):
                    launch.fork_guardian_preamble(
                        layout, 100, event_fd, directory_fd,
                        fault_injector=injected,
                    )
        finally:
            os.close(kill_r)
            os.close(kill_w)
            os.close(event_fd)
            os.close(directory_fd)
        self.assertEqual(launch.open_fd_roster(), before)

    def test_guardian_fork_return_fault_recovers_unassigned_child(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        event_fd = os.open("/dev/null", os.O_WRONLY | os.O_CLOEXEC)
        directory_fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
        layout = launch.CgroupLayout(
            service=Path("/fake"), supervisor=Path("/fake/supervisor"),
            run=Path("/fake/run"), sampler=Path("/fake/sampler"),
            experiment=Path("/fake/experiment"),
            verifier=Path("/fake/verifier"), experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1,
            run_kill_fd=kill_w,
        )
        real_fork = os.fork
        fired = False
        def fork_then_raise() -> int:
            nonlocal fired
            pid = real_fork()
            if pid > 0 and not fired:
                fired = True
                raise RuntimeError("injected guardian fork return fault")
            return pid
        def fake_guardian(control_fd: int, *_arguments: object) -> None:
            os.read(control_fd, 1)
            os._exit(0)
        try:
            with mock.patch.object(
                launch, "guardian_child", new=fake_guardian,
            ), mock.patch.object(
                launch.os, "fork", side_effect=fork_then_raise,
            ):
                with self.assertRaises(RuntimeError):
                    launch.fork_guardian_preamble(
                        layout, 100, event_fd, directory_fd,
                    )
        finally:
            os.close(kill_r)
            os.close(kill_w)
            os.close(event_fd)
            os.close(directory_fd)
        self.assertTrue(fired)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_guardian_child_fork_return_fault_never_escapes_to_caller(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        event_fd = os.open("/dev/null", os.O_WRONLY | os.O_CLOEXEC)
        directory_fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
        marker_r, marker_w = os.pipe2(os.O_CLOEXEC)
        layout = launch.CgroupLayout(
            service=Path("/fake"), supervisor=Path("/fake/supervisor"),
            run=Path("/fake/run"), sampler=Path("/fake/sampler"),
            experiment=Path("/fake/experiment"),
            verifier=Path("/fake/verifier"), experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1,
            run_kill_fd=kill_w,
        )
        parent_pid = os.getpid()
        real_fork = os.fork
        def fork_then_raise_in_child() -> int:
            pid = real_fork()
            if pid == 0:
                raise RuntimeError("injected guardian child fork return fault")
            return pid
        def quiet_child_failure(_message: str) -> None:
            os._exit(125)
        preamble = None
        try:
            with mock.patch.object(
                launch.os, "fork", side_effect=fork_then_raise_in_child,
            ), mock.patch.object(
                launch, "child_exec_failure", new=quiet_child_failure,
            ):
                try:
                    preamble = launch.fork_guardian_preamble(
                        layout, 100, event_fd, directory_fd,
                    )
                except BaseException:
                    if os.getpid() != parent_pid:
                        os.write(marker_w, b"escaped")
                        os._exit(99)
                    raise
            self.assertIsInstance(preamble, launch.GuardianSpawnPreamble)
            preamble.cleanup()
        finally:
            os.close(marker_w)
            os.close(kill_r)
            os.close(kill_w)
            os.close(event_fd)
            os.close(directory_fd)
        try:
            self.assertEqual(os.read(marker_r, 64), b"")
        finally:
            os.close(marker_r)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_guardian_spawn_return_then_raise_uses_adopted_preamble(self) -> None:
        before = launch.open_fd_roster()
        temporary = tempfile.TemporaryDirectory(prefix="wh2-guardian-spawn-")
        directory_fd = os.open(
            temporary.name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
        )
        event_fd = os.open(
            "events", os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o600, dir_fd=directory_fd,
        )
        kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
        layout = launch.CgroupLayout(
            service=Path("/fake"), supervisor=Path("/fake/supervisor"),
            run=Path("/fake/run"), sampler=Path("/fake/sampler"),
            experiment=Path("/fake/experiment"),
            verifier=Path("/fake/verifier"), experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1,
            run_kill_fd=kill_w,
        )
        real = launch.fork_guardian_preamble
        owner = launch.OwnershipSlot("test guardian preamble")
        fired = False
        def fake_guardian(control_fd: int, *_arguments: object) -> None:
            os.read(control_fd, 1)
            os._exit(0)
        def returned_then_raised(*arguments: object,
                                 **keywords: object) -> object:
            nonlocal fired
            real(*arguments, **keywords)
            fired = True
            raise RuntimeError("injected guardian preamble return fault")
        try:
            with mock.patch.object(
                launch, "guardian_child", new=fake_guardian,
            ):
                with self.assertRaises(RuntimeError):
                    returned_then_raised(
                        layout, 100, event_fd, directory_fd, owner=owner,
                    )
                self.assertIsInstance(
                    owner.value, launch.GuardianSpawnPreamble,
                )
                owner.value.cleanup()
        finally:
            os.close(kill_r)
            os.close(kill_w)
            os.close(event_fd)
            os.close(directory_fd)
            temporary.cleanup()
        self.assertTrue(fired)
        self.assertEqual(launch.open_fd_roster(), before)

    @unittest.skipUnless(hasattr(os, "pidfd_open"), "pidfd_open is unavailable")
    def test_guardian_spawn_pidfd_return_fault_closes_reused_number(self) -> None:
        before = launch.open_fd_roster()
        temporary = tempfile.TemporaryDirectory(prefix="wh2-guardian-pidfd-")
        directory_fd = os.open(
            temporary.name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
        )
        event_fd = os.open(
            "events", os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o600, dir_fd=directory_fd,
        )
        kill_r, kill_w = os.pipe2(os.O_CLOEXEC)
        layout = launch.CgroupLayout(
            service=Path("/fake"), supervisor=Path("/fake/supervisor"),
            run=Path("/fake/run"), sampler=Path("/fake/sampler"),
            experiment=Path("/fake/experiment"),
            verifier=Path("/fake/verifier"), experiment_kill_fd=-1,
            sampler_kill_fd=-1, verifier_kill_fd=-1,
            run_kill_fd=kill_w,
        )
        real_pidfd_open = os.pidfd_open
        real_fstat = os.fstat
        fired = False
        class RootOwnedStat:
            def __init__(self, value: os.stat_result) -> None:
                self.value = value
                self.st_uid = 0
                self.st_gid = 0
            def __getattr__(self, name: str) -> object:
                return getattr(self.value, name)
        def root_event_fstat(fd: int) -> object:
            value = real_fstat(fd)
            return RootOwnedStat(value) if fd == event_fd else value
        def pidfd_then_raise(pid: int, flags: int = 0) -> int:
            nonlocal fired
            fd = real_pidfd_open(pid, flags)
            if not fired:
                fired = True
                raise RuntimeError("injected guardian pidfd return fault")
            return fd
        try:
            with mock.patch.object(
                launch.CgroupTree, "move_pid",
            ), mock.patch.object(
                launch.os, "pidfd_open", side_effect=pidfd_then_raise,
            ), mock.patch.object(
                launch.os, "fstat", side_effect=root_event_fstat,
            ), mock.patch.object(
                launch, "root_process_status_fields", return_value={},
            ):
                with self.assertRaises(RuntimeError) as caught:
                    launch.spawn_guardian(
                        layout, time.monotonic_ns() + 5_000_000_000,
                        event_fd, directory_fd,
                    )
        finally:
            os.close(kill_r)
            os.close(kill_w)
            os.close(event_fd)
            os.close(directory_fd)
            temporary.cleanup()
        self.assertTrue(fired, str(caught.exception))
        self.assertEqual(launch.open_fd_roster(), before)

    @unittest.skipUnless(hasattr(os, "pidfd_open"), "pidfd_open is unavailable")
    def test_sampler_file_return_fault_closes_caller_owned_descriptors(self) -> None:
        before = launch.open_fd_roster()
        pidfd = os.pidfd_open(os.getpid(), 0)
        process = launch.ProcessHandle(
            role="sampler", pid=os.getpid(), pidfd=pidfd,
            start_ticks=launch.process_start_ticks(os.getpid()),
            argv=("/bin/true",), stdout_fd=-1, stderr_fd=-1,
        )
        def returned_then_raised(
            _process: launch.ProcessHandle,
            result: object = None,
        ) -> object:
            self.assertIsNotNone(result)
            assert isinstance(result, dict)
            for name in ("csv", "pid", "validation", "receipt"):
                result[name] = os.open(
                    "/dev/null", os.O_RDONLY | os.O_CLOEXEC,
                )
            raise RuntimeError("injected sampler-file return fault")
        try:
            with mock.patch.object(
                launch, "open_sampler_files", side_effect=returned_then_raised,
            ):
                with self.assertRaises(RuntimeError):
                    launch.wait_sampler_ready(
                        process, mock.Mock(), os.stat("/usr/bin/python3.12"),
                    )
        finally:
            os.close(pidfd)
        self.assertEqual(launch.open_fd_roster(), before)

    def test_fd_delta_return_fault_retains_recovery_evidence(self) -> None:
        baseline = launch.open_fd_roster()
        leaked = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        recovered: list[int] = []
        real = launch.close_fd_delta
        def returned_then_raised(*arguments: object,
                                 **keywords: object) -> object:
            real(*arguments, **keywords)
            raise RuntimeError("injected FD-delta return fault")
        with self.assertRaises(RuntimeError):
            returned_then_raised(baseline, (), closed_sink=recovered)
        self.assertEqual(recovered, [leaked])
        real(baseline, (), closed_sink=recovered)
        self.assertEqual(recovered, [leaked])
        self.assertEqual(launch.open_fd_roster(), baseline)

    def test_close_then_raise_never_reuses_numeric_fd_authority(self) -> None:
        fd = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        process = launch.ProcessHandle(
            role="close-fi", pid=os.getpid(), pidfd=fd,
            start_ticks=1, argv=("/dev/null",), stdout_fd=-1,
            stderr_fd=-1, returncode=0,
        )
        real_close = os.close
        fired = False
        def close_then_raise(candidate: int) -> None:
            nonlocal fired
            if candidate == fd and not fired:
                fired = True
                real_close(candidate)
                raise RuntimeError("injected close return fault")
            real_close(candidate)
        with mock.patch.object(launch.os, "close", side_effect=close_then_raise):
            launch.strict_process_close(process)
        self.assertTrue(fired)
        self.assertEqual(process.pidfd, -1)
        reused = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        try:
            self.assertEqual(reused, fd)
            launch.strict_process_close(process)
            os.fstat(reused)
        finally:
            os.close(reused)

    def test_blocked_release_marks_irreversible_before_gate_and_polls(self) -> None:
        gate_r, gate_w = os.pipe2(os.O_CLOEXEC)
        process = launch.ProcessHandle(
            role="fake", pid=os.getpid(), pidfd=-1,
            start_ticks=123, argv=("/sealed/fake",),
            stdout_fd=-1, stderr_fd=-1,
        )
        blocked = launch.BlockedRole(
            process=process, gate_fd=gate_w,
            bootstrap_security={"bootstrap": True}, expected_cpus=(123,),
            expected_groups=(), python_info=os.stat("/usr/bin/python3.12"),
        )
        order: list[str] = []
        receipt = {
            "pid": process.pid, "process_start_ticks": process.start_ticks,
        }
        try:
            with mock.patch.object(
                launch, "process_security_receipt", return_value=receipt,
            ), mock.patch.object(launch, "process_exited", return_value=False):
                returned, security = blocked.release(
                    lambda: order.append("poll"),
                    lambda: order.append(
                        "marked" if blocked.may_have_released else "unmarked"
                    ),
                )
            self.assertIs(returned, process)
            self.assertEqual(os.read(gate_r, 1), b"R")
            self.assertTrue(blocked.released)
            self.assertTrue(blocked.may_have_released)
            self.assertIs(blocked.final_security, security)
            self.assertEqual(order[0:2], ["poll", "marked"])
            self.assertGreaterEqual(order.count("poll"), 3)
        finally:
            os.close(gate_r)
            if blocked.gate_fd >= 0:
                os.close(blocked.gate_fd)

    def test_sampler_gate_ownership_cannot_erase_release_possibility(self) -> None:
        gate_r, gate_w = os.pipe2(os.O_CLOEXEC)
        process = launch.ProcessHandle(
            role="sampler", pid=os.getpid(), pidfd=-1,
            start_ticks=123, argv=("/sealed/sampler",),
            stdout_fd=-1, stderr_fd=-1,
        )
        blocked = launch.BlockedRole(
            process=process, gate_fd=gate_w,
            bootstrap_security={"bootstrap": True}, expected_cpus=(122,),
            expected_groups=(113,), python_info=os.stat("/usr/bin/python3.12"),
        )
        session = launch.ExecutionSession(
            launch.LaunchConfig(expected_harness_commit="a" * 40),
        )
        session._sampler_blocked_owner.adopt(blocked)
        try:
            with mock.patch.object(
                launch, "close_object_fd",
                side_effect=RuntimeError("injected pre-close async fault"),
            ):
                with self.assertRaises(RuntimeError):
                    blocked.release(
                        release_callback=lambda: setattr(
                            session, "sampler_may_have_started", True,
                        ),
                    )
            self.assertEqual(os.read(gate_r, 1), b"R")
            self.assertGreaterEqual(blocked.gate_fd, 0)
            self.assertTrue(blocked.may_have_released)
            self.assertTrue(session.sampler_release_possible())
            # Even a later erroneous clearing of the caller-side mirror cannot
            # waive authentication after the owned role crossed its gate.
            session.sampler_may_have_started = False
            self.assertTrue(session.sampler_release_possible())
            self.assertFalse(session.completion_licensed())
        finally:
            os.close(gate_r)
            if blocked.gate_fd >= 0:
                os.close(blocked.gate_fd)

    def test_sampler_communicate_return_fault_is_ambiguous(self) -> None:
        process = launch.ProcessHandle(
            role="sampler", pid=os.getpid(), pidfd=-1,
            start_ticks=123, argv=("/sealed/sampler",),
            stdout_fd=-1, stderr_fd=-1,
        )
        authority = launch.SamplerAuthority(
            process=process, held={}, admission={},
        )
        owner = launch.OwnershipSlot("sampler terminal test")
        def communicate(*_arguments: object, **_keywords: object) -> object:
            process.returncode = 0
            return b"UNEXPECTED-SAMPLER-STDOUT", b"", 0
        def inject(point: str) -> None:
            self.assertEqual(point, "sampler-communicate-returned")
            raise RuntimeError("injected communicate return fault")
        with mock.patch.object(
            launch, "process_exited", return_value=False,
        ), mock.patch.object(
            launch, "process_start_ticks", return_value=process.start_ticks,
        ), mock.patch.object(
            launch, "pidfd_signal",
        ), mock.patch.object(
            launch, "communicate_process", side_effect=communicate,
        ):
            with self.assertRaises(RuntimeError):
                launch.stop_and_validate_sampler(
                    authority, object(), object(), object(), owner=owner,
                    fault_injector=inject,
                )
        state = owner.value
        self.assertIsInstance(state, launch.SamplerTerminalState)
        self.assertTrue(state.communication_started)
        self.assertIsNone(state.communication)
        with self.assertRaisesRegex(
            launch.LaunchError, "communication result is ambiguous",
        ):
            launch.stop_and_validate_sampler(
                authority, object(), object(), object(), owner=owner,
            )
        self.assertIsNone(state.proof)
        self.assertIsNone(state.result)

    def test_sampler_mid_drain_failure_is_ambiguous_while_live(self) -> None:
        process = launch.ProcessHandle(
            role="sampler", pid=os.getpid(), pidfd=-1,
            start_ticks=123, argv=("/sealed/sampler",),
            stdout_fd=-1, stderr_fd=-1,
        )
        authority = launch.SamplerAuthority(
            process=process, held={}, admission={},
        )
        owner = launch.OwnershipSlot("live sampler terminal test")
        calls = 0
        def communicate(*_arguments: object, **_keywords: object) -> object:
            nonlocal calls
            calls += 1
            raise RuntimeError("injected failure after draining BAD prefix")
        with mock.patch.object(
            launch, "process_exited", return_value=False,
        ), mock.patch.object(
            launch, "process_start_ticks", return_value=process.start_ticks,
        ), mock.patch.object(
            launch, "pidfd_signal",
        ), mock.patch.object(
            launch, "communicate_process", side_effect=communicate,
        ):
            with self.assertRaisesRegex(RuntimeError, "draining BAD"):
                launch.stop_and_validate_sampler(
                    authority, object(), object(), object(), owner=owner,
                )
            self.assertIsNone(process.returncode)
            with self.assertRaisesRegex(
                launch.LaunchError, "communication result is ambiguous",
            ):
                launch.stop_and_validate_sampler(
                    authority, object(), object(), object(), owner=owner,
                )
        self.assertEqual(calls, 1)
        state = owner.value
        self.assertIsInstance(state, launch.SamplerTerminalState)
        self.assertTrue(state.communication_started)
        self.assertIsNone(state.communication)
        self.assertIsNone(state.proof)
        self.assertIsNone(state.result)

    def test_process_communication_return_fault_replays_exact_transcript(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        stdout_r, stdout_w = os.pipe2(os.O_CLOEXEC)
        stderr_r, stderr_w = os.pipe2(os.O_CLOEXEC)
        child = os.fork()
        if child == 0:
            try:
                os.close(stdout_r)
                os.close(stderr_r)
                os.write(stdout_w, b"exact-stdout")
                os.write(stderr_w, b"exact-stderr")
            finally:
                os._exit(7)
        os.close(stdout_w)
        os.close(stderr_w)
        os.set_blocking(stdout_r, False)
        os.set_blocking(stderr_r, False)
        process = launch.ProcessHandle(
            role="communication", pid=child, pidfd=-1,
            start_ticks=1, argv=("<scratch>",),
            stdout_fd=stdout_r, stderr_fd=stderr_r,
        )
        real = launch.communicate_process
        fired = False
        def returned_then_raised(*arguments: object,
                                 **keywords: object) -> object:
            nonlocal fired
            result = real(*arguments, **keywords)
            if not fired:
                fired = True
                raise RuntimeError("injected communication return fault")
            return result
        try:
            with self.assertRaisesRegex(RuntimeError, "return fault"):
                returned_then_raised(
                    process, time.monotonic() + 2.0, 1024, 1024,
                )
            result = real(process, time.monotonic() + 2.0, 1024, 1024)
            self.assertEqual(result, (b"exact-stdout", b"exact-stderr", 7))
            self.assertEqual(process.communication_result, result)
        finally:
            process.close()
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                os.waitpid(child, 0)
            except ChildProcessError:
                pass
        self.assertTrue(fired)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_process_communication_mid_read_fault_is_ambiguous(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        stdout_r, stdout_w = os.pipe2(os.O_CLOEXEC)
        stderr_r, stderr_w = os.pipe2(os.O_CLOEXEC)
        child = os.fork()
        if child == 0:
            try:
                os.close(stdout_r)
                os.close(stderr_r)
                os.write(stdout_w, b"consumed-but-unpublished")
                time.sleep(1.0)
            finally:
                os._exit(0)
        os.close(stdout_w)
        os.close(stderr_w)
        os.set_blocking(stdout_r, False)
        os.set_blocking(stderr_r, False)
        process = launch.ProcessHandle(
            role="communication", pid=child, pidfd=-1,
            start_ticks=1, argv=("<scratch>",),
            stdout_fd=stdout_r, stderr_fd=stderr_r,
        )
        real_read = os.read
        fired = False
        def read_then_raise(fd: int, maximum: int) -> bytes:
            nonlocal fired
            block = real_read(fd, maximum)
            if fd == stdout_r and block and not fired:
                fired = True
                raise RuntimeError("injected read return fault")
            return block
        try:
            with mock.patch.object(
                launch.os, "read", side_effect=read_then_raise,
            ):
                with self.assertRaisesRegex(RuntimeError, "read return fault"):
                    launch.communicate_process(
                        process, time.monotonic() + 2.0, 1024, 1024,
                    )
            with self.assertRaisesRegex(
                launch.LaunchError, "output read result is ambiguous",
            ):
                launch.communicate_process(
                    process, time.monotonic() + 2.0, 1024, 1024,
                )
        finally:
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                _waited, status_value = os.waitpid(child, 0)
                process.returncode = launch.wait_status_code(status_value)
            except ChildProcessError:
                pass
            process.close()
        self.assertTrue(fired)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_controller_transcript_publication_resumes_after_return_fault(
            self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="wh2-controller-publication-",
        ) as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            session.controller = launch.ProcessHandle(
                role="controller", pid=os.getpid(), pidfd=-1,
                start_ticks=1, argv=("<scratch>",), stdout_fd=-1,
                stderr_fd=-1, returncode=2,
                communication_result=(b"exact-out", b"exact-err", 2),
            )
            real = journal.write_bytes_resumable
            fired = False
            def returned_then_raised(name: str, payload: bytes,
                                     **keywords: object) -> object:
                nonlocal fired
                receipt = real(name, payload, **keywords)
                if name == "controller.stderr" and not fired:
                    fired = True
                    raise RuntimeError("injected controller copy return fault")
                return receipt
            try:
                with mock.patch.object(
                    journal, "write_bytes_resumable",
                    side_effect=returned_then_raised,
                ):
                    with self.assertRaisesRegex(RuntimeError, "return fault"):
                        session.publish_controller_communication()
                self.assertFalse(session.controller_communication_published)
                self.assertEqual(
                    session.publish_controller_communication(),
                    (b"exact-out", b"exact-err"),
                )
                self.assertTrue(session.controller_communication_published)
                self.assertEqual(session.controller_rc, 2)
                self.assertEqual(
                    (attempt / "controller.stdout").read_bytes(), b"exact-out",
                )
                self.assertEqual(
                    (attempt / "controller.stderr").read_bytes(), b"exact-err",
                )
                self.assertEqual(session.phases.count("controller_reaped"), 1)
            finally:
                journal.close()
            self.assertTrue(fired)

    def test_released_controller_requires_published_transcript_for_complete(
            self) -> None:
        session = launch.ExecutionSession(
            launch.LaunchConfig(expected_harness_commit="a" * 40),
        )
        session._journal_owner.adopt(mock.Mock(names=["guardian.jsonl"]))
        session._guardian_owner.adopt(mock.Mock())
        session._root_window = (1, 2)
        session.quiesced = True
        session.postflight_ok = True
        session.namespaces_sealed = True
        session.runtime_closed = True
        session.layout_removed = True
        session.guardian_authenticated = True
        session.controller_released = True
        self.assertFalse(session.completion_licensed())
        session.controller_communication_published = True
        self.assertTrue(session.completion_licensed())

    def test_controller_cancellation_publication_resumes_before_complete(
            self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="wh2-controller-cancellation-",
        ) as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            process = launch.ProcessHandle(
                role="controller", pid=os.getpid(), pidfd=-1,
                start_ticks=1, argv=("<scratch>",), stdout_fd=-1,
                stderr_fd=-1, returncode=125,
                communication_result=(b"", b"", 125),
            )
            cancelled = {
                "returncode": 125, "roster": [process.pid],
                "stderr_sha256": launch.sha256_bytes(b""),
                "stdout_sha256": launch.sha256_bytes(b""),
            }
            blocked = launch.BlockedRole(
                process=process, gate_fd=-1, bootstrap_security={},
                expected_cpus=(123,), expected_groups=(),
                python_info=os.stat("/usr/bin/python3.12"),
                cancellation_roster=[process.pid],
                cancellation_result=dict(cancelled),
            )
            session._controller_blocked_owner.adopt(blocked)
            real = journal.write_bytes_resumable
            fired = False
            def returned_then_raised(name: str, payload: bytes,
                                     **keywords: object) -> object:
                nonlocal fired
                receipt = real(name, payload, **keywords)
                if name == "controller-cancelled.json" and not fired:
                    fired = True
                    raise RuntimeError("injected cancellation copy return fault")
                return receipt
            try:
                with mock.patch.object(
                    journal, "write_bytes_resumable",
                    side_effect=returned_then_raised,
                ):
                    with self.assertRaisesRegex(RuntimeError, "return fault"):
                        session.publish_controller_cancellation(cancelled)
                self.assertFalse(session.controller_cancellation_published)
                session.publish_controller_cancellation(cancelled)
                self.assertTrue(session.controller_cancellation_published)
                self.assertEqual(session.controller_rc, 125)
                self.assertEqual(
                    (attempt / "controller-cancelled.json").read_bytes(),
                    launch.canonical_bytes(cancelled) + b"\n",
                )
                session.quiesced = True
                session.postflight_ok = True
                session.namespaces_sealed = True
                session.runtime_closed = True
                session.layout_removed = True
                self.assertTrue(session.completion_licensed())
                session.controller_cancellation_published = False
                self.assertFalse(session.completion_licensed())
            finally:
                journal.close()
            self.assertTrue(fired)

    def test_sampler_cancellation_is_required_before_complete(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="wh2-sampler-cancellation-",
        ) as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            process = launch.ProcessHandle(
                role="sampler", pid=os.getpid(), pidfd=-1,
                start_ticks=1, argv=("<scratch>",), stdout_fd=-1,
                stderr_fd=-1, returncode=125,
                communication_result=(b"", b"", 125),
            )
            cancelled = {
                "returncode": 125, "roster": [process.pid],
                "stderr_sha256": launch.sha256_bytes(b""),
                "stdout_sha256": launch.sha256_bytes(b""),
            }
            blocked = launch.BlockedRole(
                process=process, gate_fd=-1, bootstrap_security={},
                expected_cpus=(122,), expected_groups=(113,),
                python_info=os.stat("/usr/bin/python3.12"),
                cancellation_roster=[process.pid],
                cancellation_result=dict(cancelled),
            )
            session._sampler_blocked_owner.adopt(blocked)
            session.quiesced = True
            session.postflight_ok = True
            session.namespaces_sealed = True
            session.runtime_closed = True
            session.layout_removed = True
            try:
                self.assertFalse(session.completion_licensed())
                session.publish_sampler_cancellation(cancelled)
                self.assertTrue(session.sampler_cancellation_published)
                self.assertTrue(session.completion_licensed())
                self.assertEqual(
                    (attempt / "sampler-cancelled.json").read_bytes(),
                    launch.canonical_bytes(cancelled) + b"\n",
                )
            finally:
                journal.close()

    def test_released_verifier_publication_resumes_before_resolution(
            self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="wh2-verifier-publication-",
        ) as temporary:
            attempt = Path(temporary) / "attempt"
            journal = launch.AttemptJournal.reserve(
                attempt, {"attempt": 1}, expected_uid=os.getuid(),
                expected_gid=os.getgid(), test_parent=True,
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="a" * 40),
            )
            session._journal_owner.adopt(journal)
            layout = mock.Mock()
            layout.verifier = Path(temporary) / "verifier-cgroup"
            layout.verifier_kill_fd = -1
            session._layout_owner.adopt(layout)
            session._bundle_owner.adopt(mock.Mock())
            session.child_outcome = "pass"
            session.verifier = launch.ProcessHandle(
                role="verifier", pid=os.getpid(), pidfd=-1,
                start_ticks=1, argv=("<scratch>",), stdout_fd=-1,
                stderr_fd=-1, returncode=0,
                communication_result=(b'{"replay":true}\n', b"", 0),
            )
            real = journal.write_bytes_resumable
            fired = False
            def returned_then_raised(name: str, payload: bytes,
                                     **keywords: object) -> object:
                nonlocal fired
                receipt = real(name, payload, **keywords)
                if name == "replay.stderr" and not fired:
                    fired = True
                    raise RuntimeError("injected verifier copy return fault")
                return receipt
            try:
                with mock.patch.object(
                    launch, "validate_replay_result",
                    return_value={"status": "accepted"},
                ), mock.patch.object(
                    launch, "cgroup_descendant_pids", return_value=[],
                ), mock.patch.object(
                    journal, "write_bytes_resumable",
                    side_effect=returned_then_raised,
                ):
                    with self.assertRaisesRegex(RuntimeError, "return fault"):
                        session.publish_verifier_communication()
                self.assertFalse(session.verifier_resolved)
                with mock.patch.object(
                    launch, "validate_replay_result",
                    return_value={"status": "accepted"},
                ), mock.patch.object(
                    launch, "cgroup_descendant_pids", return_value=[],
                ):
                    session.publish_verifier_communication()
                self.assertTrue(session.verifier_resolved)
                self.assertEqual(session.verifier_rc, 0)
                self.assertEqual(session.replay, {"status": "accepted"})
                self.assertEqual(
                    (attempt / "replay.json").read_bytes(),
                    b'{"replay":true}\n',
                )
                self.assertEqual((attempt / "replay.stderr").read_bytes(), b"")
                self.assertEqual(session.phases.count("retained_replay"), 1)
            finally:
                journal.close()
            self.assertTrue(fired)

    def test_cancelled_role_return_fault_retains_output_and_roster(self) -> None:
        process = launch.ProcessHandle(
            role="cancelled", pid=os.getpid(), pidfd=-1,
            start_ticks=1, argv=("<scratch>",), stdout_fd=-1,
            stderr_fd=-1, returncode=125,
            communication_result=(b"retained-out", b"retained-err", 125),
        )
        blocked = launch.BlockedRole(
            process=process, gate_fd=-1, bootstrap_security={},
            expected_cpus=(123,), expected_groups=(),
            python_info=os.stat("/usr/bin/python3.12"),
        )
        real = launch.drain_and_reap_cancelled
        fired = False
        def returned_then_raised(*arguments: object,
                                 **keywords: object) -> object:
            nonlocal fired
            result = real(*arguments, **keywords)
            if not fired:
                fired = True
                raise RuntimeError("injected cancellation return fault")
            return result
        with mock.patch.object(
            launch, "cgroup_descendant_pids", return_value=[],
        ):
            with self.assertRaisesRegex(RuntimeError, "return fault"):
                returned_then_raised(blocked, Path("/scratch/cgroup"), -1)
            recovered = real(blocked, Path("/scratch/cgroup"), -1)
        expected = {
            "returncode": 125, "roster": [os.getpid()],
            "stderr_sha256": launch.sha256_bytes(b"retained-err"),
            "stdout_sha256": launch.sha256_bytes(b"retained-out"),
        }
        self.assertTrue(fired)
        self.assertEqual(recovered, expected)
        self.assertEqual(blocked.cancellation_result, expected)

    def test_maps_parser_rejects_anonymous_and_file_writable_exec(self) -> None:
        samples = (
            b"00400000-00401000 rwxp 00000000 00:00 0\n",
            b"00400000-00401000 rwxp 00000000 08:01 1 /usr/bin/python3.12\n",
        )
        for raw in samples:
            with self.subTest(raw=raw):
                with mock.patch.object(
                    launch, "read_bounded_path", return_value=raw,
                ):
                    with self.assertRaises(launch.LaunchError):
                        launch.mapped_file_paths(999999)

    def test_run_lock_parent_host_policy_is_frozen_read_only(self) -> None:
        info = os.stat("/run/lock", follow_symlinks=False)
        self.assertTrue(stat.S_ISDIR(info.st_mode))
        self.assertEqual((info.st_uid, info.st_gid), (0, 0))
        self.assertEqual(stat.S_IMODE(info.st_mode), 0o1777)


class SamplerReplayTests(unittest.TestCase):
    def test_sampler_schema_and_health_affinity_contracts_are_distinct(self) -> None:
        self.assertEqual(
            launch.SAMPLER_PRODUCER_SCHEMA,
            "wirehair.wh2.thermal_sampler.v2",
        )
        self.assertEqual(
            launch.HEALTH_SAMPLER_ATTESTATION_SCHEMA,
            "wirehair.wh2.direct-systematic-complement-sampler-attestation.v3",
        )
        self.assertEqual(
            launch.HEALTH_CONTROLLER_INITIAL_AFFINITY, (120, 121, 122, 123),
        )

    def test_exact_state_machine_and_summary_replay(self) -> None:
        raw, validation, receipt = sampler_fixture()
        replay = launch.validate_sampler_stream_replay(raw, validation, receipt)
        self.assertEqual(replay["sample_count"], 1)
        self.assertEqual(replay["cpu_tctl_max_c"], 50.0)

    def test_forged_sensors_errors_edac_and_thin_validation_are_rejected(self) -> None:
        raw, validation, receipt = sampler_fixture()
        lines = raw.decode("ascii").splitlines()
        fields = lines[1].split(",")
        mutations = []
        forged = list(fields)
        forged[5] = "999.0"
        mutations.append((lines[0] + "\n" + ",".join(forged) + "\n").encode("ascii"))
        forged = list(fields)
        forged[13] = "123"
        mutations.append((lines[0] + "\n" + ",".join(forged) + "\n").encode("ascii"))
        forged = list(fields)
        forged[17], forged[18] = "456", "789"
        mutations.append((lines[0] + "\n" + ",".join(forged) + "\n").encode("ascii"))
        for forged_raw in mutations:
            with self.assertRaises(launch.LaunchError):
                launch.validate_sampler_stream_replay(
                    forged_raw, validation, receipt,
                )
        values = [json.loads(line) for line in validation.decode("ascii").splitlines()]
        values[1] = {
            "schema": launch.SAMPLER_VALIDATION_SCHEMA,
            "decision": "continue", "sample_index": 0, "monotonic_s": 1.0,
        }
        thin = b"".join(
            launch.canonical_bytes(value) + b"\n" for value in values
        )
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(raw, thin, receipt)

    def test_summary_and_sensor_mutations_are_rejected(self) -> None:
        raw, validation, receipt = sampler_fixture()
        bad_receipt = copy.deepcopy(receipt)
        bad_receipt["summary"]["dimm_attempt_errors_total"] = 1
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(raw, validation, bad_receipt)
        values = [json.loads(line) for line in validation.decode("ascii").splitlines()]
        values[1]["sensors"][launch.SAMPLER_DIMM_FIELDS[0]]["hot"] = True
        bad_validation = b"".join(
            launch.canonical_bytes(value) + b"\n" for value in values
        )
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(raw, bad_validation, receipt)

    def test_validation_timestamp_type_and_raw_domains_are_exact(self) -> None:
        raw, validation, receipt = sampler_fixture()
        values = [json.loads(line) for line in validation.decode("ascii").splitlines()]
        values[1]["monotonic_s"] = 1
        integer_time = b"".join(
            launch.canonical_bytes(value) + b"\n" for value in values
        )
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(raw, integer_time, receipt)

        lines = raw.decode("ascii").splitlines()
        fields = lines[1].split(",")
        fields[5] = "45.1"
        off_grid = (lines[0] + "\n" + ",".join(fields) + "\n").encode("ascii")
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(off_grid, validation, receipt)

        fields = lines[1].split(",")
        fields[0] = "9999-99-99T99:99:99.999Z"
        bad_utc = (lines[0] + "\n" + ",".join(fields) + "\n").encode("ascii")
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(bad_utc, validation, receipt)

        precise_raw, precise_validation, precise_receipt = sampler_fixture()
        precise_lines = precise_raw.decode("ascii").splitlines()
        precise_fields = precise_lines[1].split(",")
        precise_fields[1] = "1.123457"
        precise_fields[0] = "2026-08-29T00:00:00.123Z"
        precise_raw = (
            precise_lines[0] + "\n" + ",".join(precise_fields) + "\n"
        ).encode("ascii")
        precise_values = [
            json.loads(line) for line in precise_validation.decode("ascii").splitlines()
        ]
        precise_values[1]["monotonic_s"] = 1.123456789
        precise_validation = b"".join(
            launch.canonical_bytes(value) + b"\n" for value in precise_values
        )
        replay = launch.validate_sampler_stream_replay(
            precise_raw, precise_validation, precise_receipt,
        )
        self.assertEqual(replay["sample_count"], 1)

    def test_duplicate_formatted_raw_timestamp_is_rejected(self) -> None:
        raw, validation, receipt = sampler_fixture()
        raw_lines = raw.decode("ascii").splitlines()
        second_row = raw_lines[1].split(",")
        second_row[0] = "2026-08-29T00:00:00.001Z"
        # The producer retains full float precision in validation JSON but six
        # decimals in CSV.  Distinct producer times may therefore collide in
        # the formatted stream and must not create an ambiguous pairing key.
        second_row[1] = "1.000000"
        raw = (
            raw_lines[0] + "\n" + raw_lines[1] + "\n"
            + ",".join(second_row) + "\n"
        ).encode("ascii")
        values = [
            json.loads(line) for line in validation.decode("ascii").splitlines()
        ]
        second_sample = copy.deepcopy(values[1])
        second_sample["sample_index"] = 1
        second_sample["monotonic_s"] = 1.0000004
        for sensor in second_sample["sensors"].values():
            sensor["jump_c"] = 0.0
            sensor["rate_c_per_s"] = 0.0
        validation = b"".join(
            launch.canonical_bytes(value) + b"\n"
            for value in (values[0], values[1], second_sample)
        )
        receipt = copy.deepcopy(receipt)
        receipt["summary"]["sample_count"] = 2
        for sensor in receipt["summary"]["sensors"].values():
            sensor["raw_samples"] = 2
            sensor["valid_samples"] = 2
        with self.assertRaises(launch.LaunchError):
            launch.validate_sampler_stream_replay(raw, validation, receipt)

    def test_health_prefixes_require_complete_paired_records(self) -> None:
        raw, validation, _ = sampler_fixture()
        prefix = launch.validate_health_sampler_prefixes(
            raw, validation, "test health prefix",
        )
        self.assertEqual(prefix["sample_ns"], [1_000_000_000])
        for bad_raw, bad_validation in (
            (raw[:-1], validation),
            (raw, validation[:-1]),
            (raw, validation.splitlines(keepends=True)[0]),
            (raw.splitlines(keepends=True)[0], validation),
        ):
            with self.subTest(
                raw_bytes=len(bad_raw), validation_bytes=len(bad_validation),
            ):
                with self.assertRaises(launch.LaunchError):
                    launch.validate_health_sampler_prefixes(
                        bad_raw, bad_validation, "bad health prefix",
                    )

    def test_health_window_is_the_exact_contiguous_terminal_slice(self) -> None:
        raw, validation, _ = sampler_fixture()
        raw_lines = raw.splitlines(keepends=True)
        validation_values = [
            json.loads(line) for line in validation.decode("ascii").splitlines()
        ]
        terminal_raw = [raw_lines[0]]
        terminal_validation = [
            launch.canonical_bytes(validation_values[0]) + b"\n",
        ]
        for index, timestamp in enumerate((1.0, 1.5, 2.0)):
            row = raw_lines[1].decode("ascii").rstrip("\n").split(",")
            row[0] = "2026-08-29T00:00:0{}.000Z".format(index)
            row[1] = "{:.6f}".format(timestamp)
            terminal_raw.append((",".join(row) + "\n").encode("ascii"))
            sample = copy.deepcopy(validation_values[1])
            sample["sample_index"] = index
            sample["monotonic_s"] = timestamp
            terminal_validation.append(launch.canonical_bytes(sample) + b"\n")
        terminal_raw_bytes = b"".join(terminal_raw)
        terminal_validation_bytes = b"".join(terminal_validation)
        expected_raw, expected_validation = launch.select_exact_health_window(
            terminal_raw_bytes, terminal_validation_bytes,
            1_000_000_000, 2_000_000_000, "test health window",
        )
        self.assertEqual(expected_raw, terminal_raw_bytes)
        self.assertEqual(expected_validation, terminal_validation_bytes)
        forged_raw = terminal_raw[0] + terminal_raw[1] + terminal_raw[3]
        forged_validation = (
            terminal_validation[0] + terminal_validation[1]
            + terminal_validation[3]
        )
        self.assertNotEqual(forged_raw, expected_raw)
        self.assertNotEqual(forged_validation, expected_validation)
        with self.assertRaises(launch.LaunchError):
            launch.select_exact_health_window(
                terminal_raw[0] + terminal_raw[1] + terminal_raw[2],
                terminal_validation[0] + terminal_validation[1]
                + terminal_validation[2],
                1_000_000_000, 2_000_000_000,
                "truncated health prefix",
            )


class NamespaceSealTests(unittest.TestCase):
    def test_bound_directory_create_and_seal_are_idempotent(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-bound-") as temporary:
            root = Path(temporary)
            parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            directory = launch.BoundDirectory(
                root / "owned", parent_fd, os.getuid(), os.getgid(),
            )
            try:
                directory.create()
                first_binding = directory.binding
                directory.create()
                self.assertEqual(directory.binding, first_binding)
                directory.seal(())
                first = directory.verify(0o500)
                directory.seal(())
                second = directory.verify(0o500)
                self.assertEqual(launch.full_stat(first), launch.full_stat(second))
            finally:
                directory.close()
                os.chmod(root / "owned", 0o700)
                os.rmdir(root / "owned")
                os.close(parent_fd)

    def test_bound_directory_create_rollback_close_is_idempotent(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-bound-fi-") as temporary:
            root = Path(temporary)
            parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            directory = launch.BoundDirectory(
                root / "owned", parent_fd, os.getuid(), os.getgid(),
            )
            real_close = os.close
            closed_fd = -1
            fired = False
            target_fd = -1
            def rename_fault(*_arguments: object, **_keywords: object) -> None:
                nonlocal target_fd
                target_fd = directory.fd
                raise RuntimeError("injected pre-rename fault")
            def close_then_raise(candidate: int) -> None:
                nonlocal closed_fd, fired
                if candidate == target_fd and candidate >= 0 and not fired:
                    closed_fd = candidate
                    fired = True
                    real_close(candidate)
                    raise RuntimeError("injected rollback close return fault")
                real_close(candidate)
            try:
                with mock.patch.object(
                    launch, "rename_noreplace_at",
                    side_effect=rename_fault,
                ), mock.patch.object(
                    launch.os, "close", side_effect=close_then_raise,
                ):
                    with self.assertRaises(RuntimeError):
                        directory.create()
                self.assertTrue(fired)
                self.assertEqual(directory.fd, -1)
                self.assertEqual(directory.state, "absent")
                self.assertIsNone(directory.binding)
                self.assertEqual(os.listdir(root), [])
                reused = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
                try:
                    self.assertEqual(reused, closed_fd)
                    directory.close()
                    os.fstat(reused)
                finally:
                    os.close(reused)
            finally:
                directory.close()
                os.close(parent_fd)

    def test_absent_precreation_namespaces_are_a_sealed_invalid_state(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-namespaces-") as temporary:
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="d" * 40),
            )
            session.var_tmp_fd = os.open(
                temporary, os.O_RDONLY | os.O_DIRECTORY,
            )
            try:
                session.seal_namespaces()
                session.quiesce()
                self.assertTrue(session.namespaces_sealed)
                self.assertTrue(session.quiesced)
            finally:
                os.close(session.var_tmp_fd)
                session.var_tmp_fd = -1

    def test_mutable_partial_child_forbids_parent_seal(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-facade-namespaces-") as temporary:
            root = Path(temporary)
            parent_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
            controller = launch.BoundDirectory(
                root / "controller", parent_fd, os.getuid(), os.getgid(),
            )
            sampler = launch.BoundDirectory(
                root / "sampler", parent_fd, os.getuid(), os.getgid(),
            )
            session = launch.ExecutionSession(
                launch.LaunchConfig(expected_harness_commit="d" * 40),
            )
            try:
                controller.create()
                sampler.create()
                os.mkdir(launch.CONTROLLER_OUTPUT.name, 0o700,
                         dir_fd=controller.fd)
                session.controller_directory = controller
                session.sampler_directory = sampler
                with self.assertRaises(launch.LaunchError):
                    session.seal_namespaces()
                self.assertFalse(session.namespaces_sealed)
                self.assertEqual(
                    stat.S_IMODE(os.fstat(controller.fd).st_mode), 0o700,
                )
            finally:
                try:
                    os.rmdir(launch.CONTROLLER_OUTPUT.name, dir_fd=controller.fd)
                except OSError:
                    pass
                controller.close()
                sampler.close()
                os.close(parent_fd)


class BuildVectorTests(unittest.TestCase):
    def test_snapshot_clone_environment_binds_exact_campaign_owner(self) -> None:
        if os.geteuid() not in (0, launch.CAMPAIGN_UID):
            self.skipTest("campaign-owner fixture requires UID0 or UID1000")
        with tempfile.TemporaryDirectory(prefix="wh2-clone-environment-") as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            git_directory = source / ".git"
            git_directory.mkdir()
            if os.geteuid() != launch.CAMPAIGN_UID:
                os.chown(source, launch.CAMPAIGN_UID, launch.CAMPAIGN_GID)
                os.chown(git_directory, launch.CAMPAIGN_UID, launch.CAMPAIGN_GID)
            environment = launch.snapshot_clone_environment(source)
            self.assertEqual(environment["SUDO_UID"], str(launch.CAMPAIGN_UID))
            self.assertEqual(
                {key: value for key, value in environment.items()
                 if key != "SUDO_UID"},
                launch.GIT_ENVIRONMENT,
            )

    def test_snapshot_clone_environment_rejects_other_owner(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-clone-owner-") as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            (source / ".git").mkdir()
            wrong_uid = 0 if launch.CAMPAIGN_UID != 0 else 1
            real_stat = launch.os.stat
            def owner_stat(path: str, *, follow_symlinks: bool = True) -> os.stat_result:
                value = real_stat(path, follow_symlinks=follow_symlinks)
                if Path(path) == source:
                    fields = list(value)
                    fields[4] = wrong_uid
                    return os.stat_result(fields)
                return value
            with mock.patch.object(
                launch.os, "stat", side_effect=owner_stat,
            ):
                with self.assertRaises(launch.LaunchError):
                    launch.snapshot_clone_environment(source)

    @unittest.skipUnless(os.geteuid() == 0, "root snapshot clone probe")
    def test_materialize_snapshot_clones_campaign_owned_input_as_root(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        with tempfile.TemporaryDirectory(prefix="wh2-root-snapshot-") as temporary:
            parent = Path(temporary)
            source = parent / "source"
            source.mkdir()
            git(source, "init", "--quiet")
            git(source, "config", "user.name", "Facade Snapshot Test")
            git(source, "config", "user.email", "facade@example.invalid")
            tracked = source / "source.cpp"
            tracked.write_text("int value = 1;\n", encoding="ascii")
            git(source, "add", "source.cpp")
            git(source, "commit", "--quiet", "-m", "snapshot fixture")
            git(source, "checkout", "--quiet", "--detach")
            commit = git(source, "rev-parse", "HEAD").decode("ascii").strip()
            for raw_root, directories, files in os.walk(str(source)):
                for name in directories + files:
                    os.chown(
                        str(Path(raw_root) / name),
                        launch.CAMPAIGN_UID, launch.CAMPAIGN_GID,
                        follow_symlinks=False,
                    )
            os.chown(source, launch.CAMPAIGN_UID, launch.CAMPAIGN_GID)
            destination = parent / "snapshot"
            receipt = launch.materialize_snapshot(source, destination, commit)
            self.assertEqual(receipt.commit, commit)
            self.assertEqual(
                [entry["repository_relative_path"] for entry in receipt.entries],
                ["source.cpp"],
            )
            source_info = os.stat(str(source / "source.cpp"), follow_symlinks=False)
            copied_info = os.stat(str(destination / "source.cpp"), follow_symlinks=False)
            self.assertNotEqual(
                (source_info.st_dev, source_info.st_ino),
                (copied_info.st_dev, copied_info.st_ino),
            )
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_fresh_directory_establishes_mode_under_service_umask(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="wh2-facade-fresh-directory-",
        ) as temporary:
            parent = Path(temporary)
            child = parent / "sealed-root"
            real_stat = launch.os.stat

            class RootProtectedParent:
                def __init__(self, value: os.stat_result) -> None:
                    self._value = value
                    self.st_mode = stat.S_IFDIR | 0o755
                    self.st_uid = 0
                    self.st_gid = 0

                def __getattr__(self, name: str) -> object:
                    return getattr(self._value, name)

            def root_parent_stat(path: object, *arguments: object,
                                 **keywords: object) -> object:
                value = real_stat(path, *arguments, **keywords)
                if Path(path) == parent:
                    return RootProtectedParent(value)
                return value

            previous_umask = os.umask(0o077)
            try:
                with mock.patch.object(
                    launch.os, "stat", side_effect=root_parent_stat,
                ):
                    launch._fresh_directory(
                        child, uid=os.getuid(), gid=os.getgid(), mode=0o755,
                    )
            finally:
                os.umask(previous_umask)
            info = child.stat()
            self.assertEqual(info.st_uid, os.getuid())
            self.assertEqual(info.st_gid, os.getgid())
            self.assertEqual(stat.S_IMODE(info.st_mode), 0o755)

    def test_retained_replay_preserves_invalid_aa_disposition(self) -> None:
        raw = b"synthetic-raw\n"
        failures = ["aa_ci:D/D:prebuilt-systematic:K8:B64"]
        worker_terminals = {
            "current": {"invocation_count": 7},
            "parent": {"invocation_count": 9},
        }
        bundle = mock.Mock()
        bundle.data = {"raw.jsonl": raw}
        bundle.complete = {"outcome": "invalid"}
        bundle.summary = {
            "failures": failures,
            "outcome": "invalid",
            "statistics": {"failed_gates": failures},
            "worker_terminals": worker_terminals,
        }
        replay = {
            "campaign": launch.CAMPAIGN,
            "failed_gates": failures,
            "outcome": "invalid",
            "raw_bytes": len(raw),
            "raw_sha256": launch.sha256_bytes(raw),
            "raw_terminal": {
                "current_invocations": 7,
                "parent_invocations": 9,
            },
            "schema": launch.REPLAY_SCHEMA,
            "statistics": bundle.summary["statistics"],
        }
        encoded = launch.canonical_bytes(replay) + b"\n"
        self.assertEqual(
            launch.validate_replay_result(encoded, b"", 1, bundle), replay,
        )
        with self.assertRaises(launch.LaunchError):
            launch.validate_replay_result(encoded, b"", 2, bundle)

    def test_compile_database_and_ninja_vectors_are_exactly_distinct(self) -> None:
        build = Path("/sealed/build")
        source = Path("/sealed/source/worker.cpp")
        relative_object = "CMakeFiles/worker.dir/worker.cpp.o"
        prefix = ["-O3", "-DNDEBUG", "-std=gnu++11"]
        compile_database, effective_ninja = launch.exact_compile_argv_pair(
            build, source, relative_object, prefix,
        )
        self.assertEqual(compile_database, [
            str(launch.CXX), *prefix,
            "-o", relative_object, "-c", str(source),
        ])
        self.assertEqual(effective_ninja, [
            str(launch.CXX), *prefix,
            "-MD", "-MT", relative_object,
            "-MF", relative_object + ".d",
            "-o", relative_object, "-c", str(source),
        ])
        self.assertNotEqual(compile_database, effective_ninja)

    def test_role_link_vector_preserves_ninja_origin_escape(self) -> None:
        build = Path("/sealed/build")
        output = build / "worker"
        objects = [build / "CMakeFiles/worker.dir/worker.cpp.o"]
        argv = launch.role_link_argv(build, output, objects)
        self.assertIn("-Wl,-rpath,\\$ORIGIN", argv)
        self.assertNotIn("-Wl,-rpath,$ORIGIN", argv)

    def test_build_receipt_link_count_policy_distinguishes_system_dsos(self) -> None:
        role_library = "/sealed/role/libwirehair.so.2"
        self.assertEqual(
            launch.build_receipt_expected_nlink(
                "dynamic_dependencies", role_library, role_library,
            ),
            1,
        )
        self.assertIsNone(
            launch.build_receipt_expected_nlink(
                "dynamic_dependencies", "/usr/lib/libc.so.6", role_library,
            )
        )
        for roster in (
            "source_manifest", "object_closure", "archive_closure",
        ):
            self.assertEqual(
                launch.build_receipt_expected_nlink(
                    roster, "/sealed/artifact", role_library,
                ),
                1,
            )
        with self.assertRaises(launch.LaunchError):
            launch.build_receipt_expected_nlink(
                "unexpected", "/sealed/artifact", role_library,
            )

    @unittest.skipUnless(
        hasattr(os, "pidfd_open")
        and hasattr(launch.signal, "pidfd_send_signal"),
        "pidfd command-containment primitives are unavailable",
    )
    def test_bounded_command_popen_init_return_fault_reaps_exactly(self) -> None:
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        real_init = launch.subprocess.Popen.__init__
        fired = False
        def init_then_raise(process: object, *arguments: object,
                            **keywords: object) -> None:
            nonlocal fired
            real_init(process, *arguments, **keywords)
            if not fired:
                fired = True
                raise RuntimeError("injected Popen init return fault")
        with mock.patch.object(
            launch.subprocess.Popen, "__init__", new=init_then_raise,
        ):
            with self.assertRaises(BaseException):
                launch._bounded_command(
                    ["/bin/sleep", "30"], cwd=None,
                    environment={"PATH": "/usr/bin:/bin"}, timeout=2.0,
                )
        self.assertTrue(fired)
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    @unittest.skipUnless(
        hasattr(os, "fork") and hasattr(os, "pidfd_open")
        and hasattr(launch.signal, "pidfd_send_signal"),
        "subreaper command-containment primitives are unavailable",
    )
    def test_bounded_command_reaps_exited_adopted_child_without_false_leak(
        self,
    ) -> None:
        # The command intentionally exits without waiting for an already-dead
        # child.  Linux reparents that zombie to this process after
        # BoundedCommandOwnership enables subreaper mode.  It must be reaped,
        # but it is not a live descendant and must not invalidate rc=0.
        program = (
            "import os\n"
            "child=os.fork()\n"
            "if child == 0:\n"
            " os._exit(0)\n"
            "os.waitid(os.P_PID,child,os.WEXITED|os.WNOWAIT)\n"
            "os._exit(0)\n"
        )
        before_fds = launch.open_fd_roster()
        before_children = launch.direct_child_pids()
        output, error = launch._bounded_command(
            ["/usr/bin/python3.12", "-c", program], cwd=None,
            environment={"PATH": "/usr/bin:/bin"}, timeout=5.0,
        )
        self.assertEqual(output, b"")
        self.assertEqual(error, b"")
        self.assertEqual(launch.direct_child_pids(), before_children)
        self.assertEqual(launch.open_fd_roster(), before_fds)

    @unittest.skipUnless(
        hasattr(os, "fork") and hasattr(os, "pidfd_open")
        and hasattr(launch.signal, "pidfd_send_signal"),
        "subreaper command-containment primitives are unavailable",
    )
    def test_bounded_command_does_not_signal_unreparented_zombie(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-zombie-parent-") as temporary:
            pid_path = Path(temporary) / "zombie.pid"
            baseline = launch.open_fd_roster()
            child_baseline = tuple(launch.direct_child_pids())
            owner = launch.BoundedCommandOwnership(
                label="unreparented zombie probe", fd_baseline=baseline,
                child_baseline=child_baseline,
                service_baseline=(os.getpid(),),
            )
            parent_pid = os.fork()
            if parent_pid == 0:
                zombie_pid = os.fork()
                if zombie_pid == 0:
                    os._exit(23)
                os.waitid(
                    os.P_PID, zombie_pid, os.WEXITED | os.WNOWAIT,
                )
                pid_path.write_text(str(zombie_pid), encoding="ascii")
                time.sleep(30.0)
                os._exit(0)
            zombie_pid = -1
            try:
                deadline = time.monotonic() + 2.0
                raw_zombie_pid = ""
                while time.monotonic() < deadline:
                    try:
                        raw_zombie_pid = pid_path.read_text(encoding="ascii")
                    except FileNotFoundError:
                        pass
                    if raw_zombie_pid.isdigit():
                        break
                    time.sleep(0.001)
                self.assertTrue(raw_zombie_pid.isdigit())
                zombie_pid = int(raw_zombie_pid)
                def exact_service_roster() -> list:
                    result = [os.getpid()]
                    for candidate in (parent_pid, zombie_pid):
                        if Path("/proc/{}".format(candidate)).is_dir():
                            result.append(candidate)
                    return sorted(result)
                with mock.patch.object(
                    launch, "current_service_pids",
                    side_effect=exact_service_roster,
                ):
                    descendants = owner.quiesce()
                self.assertIn(parent_pid, descendants)
                self.assertNotIn(zombie_pid, descendants)
                for candidate in (parent_pid, zombie_pid):
                    with self.assertRaises(ProcessLookupError):
                        os.kill(candidate, 0)
                    with self.assertRaises(ChildProcessError):
                        os.waitpid(candidate, os.WNOHANG)
            finally:
                for candidate in (parent_pid, zombie_pid):
                    if candidate <= 0:
                        continue
                    try:
                        os.kill(candidate, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    try:
                        os.waitpid(candidate, 0)
                    except ChildProcessError:
                        pass
        self.assertEqual(launch.direct_child_pids(), list(child_baseline))
        self.assertEqual(launch.open_fd_roster(), baseline)

    def test_bounded_command_does_not_classify_disappeared_cgroup_pid_live(
        self,
    ) -> None:
        baseline = launch.open_fd_roster()
        owner = launch.BoundedCommandOwnership(
            label="disappeared descendant probe", fd_baseline=baseline,
            child_baseline=(), service_baseline=(os.getpid(),),
        )
        service_rosters = iter([
            [os.getpid(), 999_999_999], [os.getpid()], [os.getpid()],
        ])
        with mock.patch.object(
            launch, "current_service_pids",
            side_effect=lambda: next(service_rosters),
        ), mock.patch.object(
            launch, "direct_child_pids", return_value=[],
        ), mock.patch.object(
            launch.os, "waitpid", side_effect=ChildProcessError,
        ), mock.patch.object(
            owner, "_signal_pid", return_value=False,
        ):
            self.assertEqual(owner.quiesce(), [])
        self.assertEqual(launch.open_fd_roster(), baseline)

    def test_bounded_command_ignores_shared_cgroup_sibling_arrivals(self) -> None:
        baseline = launch.open_fd_roster()
        owner = launch.BoundedCommandOwnership(
            label="shared service probe", fd_baseline=baseline,
            child_baseline=(), service_baseline=(os.getpid(), 111_111_111),
        )
        with mock.patch.object(
            launch, "current_service_pids",
            side_effect=AssertionError("shared service roster must not be consumed"),
        ), mock.patch.object(
            launch, "direct_child_pids", return_value=[],
        ), mock.patch.object(
            owner, "_signal_pid",
            side_effect=AssertionError("unrelated sibling must not be signaled"),
        ):
            self.assertEqual(owner.quiesce(), [])
        self.assertEqual(launch.open_fd_roster(), baseline)

    def test_bounded_descendant_signal_does_not_hide_roster_failure(self) -> None:
        baseline = launch.open_fd_roster()
        def fake_pidfd_open(pid: int, flags: int = 0) -> int:
            del pid, flags
            return os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        with mock.patch.object(
            launch, "process_start_ticks", return_value=42,
        ), mock.patch.object(
            launch.os, "pidfd_open", side_effect=fake_pidfd_open, create=True,
        ), mock.patch.object(
            launch, "current_service_pids",
            side_effect=FileNotFoundError("injected authority roster failure"),
        ):
            with self.assertRaises(FileNotFoundError):
                launch.BoundedCommandOwnership._signal_pid(999_999_999)
        self.assertEqual(launch.open_fd_roster(), baseline)

    def test_bounded_descendant_natural_exit_preserves_live_observation(self) -> None:
        baseline = launch.open_fd_roster()
        read_fd, write_fd = os.pipe2(os.O_CLOEXEC)
        try:
            def fake_pidfd_open(pid: int, flags: int = 0) -> int:
                del pid, flags
                return os.dup(read_fd)
            with mock.patch.object(
                launch, "process_start_ticks", return_value=42,
            ), mock.patch.object(
                launch.os, "pidfd_open", side_effect=fake_pidfd_open,
                create=True,
            ), mock.patch.object(
                launch, "current_service_pids", return_value=[999_999_999],
            ), mock.patch.object(
                launch.signal, "pidfd_send_signal",
                side_effect=ProcessLookupError, create=True,
            ):
                self.assertTrue(
                    launch.BoundedCommandOwnership._signal_pid(999_999_999)
                )
        finally:
            os.close(read_fd)
            os.close(write_fd)
        self.assertEqual(launch.open_fd_roster(), baseline)

    @unittest.skipUnless(
        hasattr(os, "pidfd_open")
        and hasattr(launch.signal, "pidfd_send_signal"),
        "pidfd command-containment primitives are unavailable",
    )
    def test_bounded_command_rejects_and_reaps_setsid_descendant(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-command-desc-") as temporary:
            pid_path = Path(temporary) / "pid"
            program = (
                "import pathlib,subprocess,sys;"
                "p=subprocess.Popen(['/bin/sleep','30'],"
                "start_new_session=True);"
                "pathlib.Path(sys.argv[1]).write_text(str(p.pid),"
                "encoding='ascii')"
            )
            before_fds = launch.open_fd_roster()
            before_children = launch.direct_child_pids()
            with self.assertRaises(launch.LaunchError):
                launch._bounded_command(
                    ["/usr/bin/python3.12", "-c", program, str(pid_path)],
                    cwd=None, environment={"PATH": "/usr/bin:/bin"},
                    timeout=5.0,
                )
            escaped_pid = int(pid_path.read_text(encoding="ascii"))
            with self.assertRaises(ProcessLookupError):
                os.kill(escaped_pid, 0)
            with self.assertRaises(ChildProcessError):
                os.waitpid(escaped_pid, os.WNOHANG)
            self.assertEqual(launch.direct_child_pids(), before_children)
            self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_build_child_rejects_unserializable_argv_before_acquisition(self) -> None:
        before = launch.open_fd_roster()
        with self.assertRaises(launch.LaunchError):
            launch.bounded_build_child(
                [object()], Path("/"), {"PATH": "/usr/bin:/bin"}, 1.0,
                stdout_limit=1024, stderr_limit=1024,
            )
        self.assertEqual(launch.open_fd_roster(), before)

    @unittest.skipUnless(
        hasattr(os, "pidfd_open")
        and hasattr(launch.signal, "pidfd_send_signal"),
        "pidfd command-containment primitives are unavailable",
    )
    def test_build_child_pidfd_return_fault_reaps_and_closes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-build-pidfd-") as temporary:
            before_fds = launch.open_fd_roster()
            before_children = launch.direct_child_pids()
            real_pidfd_open = os.pidfd_open
            fired = False
            def pidfd_then_raise(pid: int, flags: int = 0) -> int:
                nonlocal fired
                fd = real_pidfd_open(pid, flags)
                if not fired:
                    fired = True
                    raise RuntimeError("injected pidfd return fault")
                return fd
            with mock.patch.object(
                launch.os, "pidfd_open", side_effect=pidfd_then_raise,
            ):
                with self.assertRaises(RuntimeError):
                    launch.bounded_build_child(
                        ["/bin/true"], Path(temporary), {
                            "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
                            "PATH": "/usr/bin:/bin", "TZ": "UTC",
                        }, 2.0, stdout_limit=1024, stderr_limit=1024,
                    )
            self.assertTrue(fired)
            self.assertEqual(launch.direct_child_pids(), before_children)
            self.assertEqual(launch.open_fd_roster(), before_fds)

    def test_frozen_public_facade_source_hashes_are_exact(self) -> None:
        expected = {
            launch.CONTROLLER_RELATIVE: launch.FROZEN_CONTROLLER_SHA256,
            launch.WORKER_RELATIVE: launch.FROZEN_WORKER_SOURCE_SHA256,
            launch.ROLE_CMAKE_RELATIVE: launch.FROZEN_ROLE_CMAKE_SHA256,
            launch.HEALTH_ADAPTER_RELATIVE:
                launch.FROZEN_HEALTH_ADAPTER_SHA256,
        }
        root = MODULE_PATH.parent.parent
        for relative, digest in expected.items():
            with self.subTest(relative=relative):
                self.assertEqual(
                    launch.sha256_bytes((root / relative).read_bytes()), digest,
                )
        screen_path = root / launch.CONTROLLER_RELATIVE
        screen_tree = ast.parse(
            screen_path.read_text(encoding="utf-8"),
            filename=str(screen_path),
        )
        screen_worker_hash = next(
            ast.literal_eval(node.value)
            for node in screen_tree.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "FROZEN_WORKER_SOURCE_SHA256"
                for target in node.targets
            )
        )
        self.assertEqual(
            screen_worker_hash, launch.FROZEN_WORKER_SOURCE_SHA256,
            "launcher and controller frozen worker hashes drifted",
        )

    def test_frozen_health_git_interface_is_cross_file_exact(self) -> None:
        root = MODULE_PATH.parent.parent
        screen_path = root / launch.CONTROLLER_RELATIVE
        adapter_path = root / launch.HEALTH_ADAPTER_RELATIVE
        screen_source = screen_path.read_text(encoding="utf-8")
        adapter_source = adapter_path.read_text(encoding="utf-8")
        compile(screen_source, str(screen_path), "exec", dont_inherit=True)
        compile(adapter_source, str(adapter_path), "exec", dont_inherit=True)
        screen_tree = ast.parse(screen_source, filename=str(screen_path))
        adapter_tree = ast.parse(adapter_source, filename=str(adapter_path))

        adapter_function = next(
            node for node in adapter_tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "git_receipt"
        )
        self.assertEqual(
            [argument.arg for argument in adapter_function.args.args],
            ["root", "expected_commit", "expected_git_sha256", "deadline"],
        )
        self.assertEqual(len(adapter_function.args.defaults), 1)
        self.assertIsInstance(adapter_function.args.defaults[0], ast.Constant)
        self.assertIsNone(adapter_function.args.defaults[0].value)

        nested_git = next(
            node for node in adapter_function.body
            if isinstance(node, ast.FunctionDef) and node.name == "git"
        )
        popen = next(
            node for node in ast.walk(nested_git)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "Popen"
        )
        expected_argv = ast.parse(
            '[git_fd_path,"-c","core.fsmonitor=false","-c",'
            '"safe.directory="+str(root),*arguments]',
            mode="eval",
        ).body
        self.assertEqual(
            ast.dump(popen.args[0], include_attributes=False),
            ast.dump(expected_argv, include_attributes=False),
        )

        load_function = next(
            node for node in screen_tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "load_health_adapter"
        )
        adapter_calls = [
            node for node in ast.walk(load_function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "git_receipt"
        ]
        self.assertEqual(len(adapter_calls), 1)
        self.assertEqual(
            [ast.dump(value, include_attributes=False)
             for value in adapter_calls[0].args],
            [ast.dump(ast.parse(value, mode="eval").body,
                      include_attributes=False) for value in (
                "source_root", "expected_harness_commit", "expected_git_sha256",
            )],
        )

        run_function = next(
            node for node in screen_tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "run_once"
        )
        current_git_expression = ast.dump(ast.parse(
            'prepared.build_receipts["current"]["git_sha256"]', mode="eval",
        ).body, include_attributes=False)
        run_load_calls = [
            node for node in ast.walk(run_function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "load_health_adapter"
        ]
        run_git_calls = [
            node for node in ast.walk(run_function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "git_receipt"
        ]
        self.assertEqual(len(run_load_calls), 1)
        self.assertEqual(len(run_load_calls[0].args), 5)
        self.assertEqual(
            ast.dump(run_load_calls[0].args[4], include_attributes=False),
            current_git_expression,
        )
        self.assertEqual(len(run_git_calls), 1)
        self.assertEqual(len(run_git_calls[0].args), 3)
        self.assertEqual(
            ast.dump(run_git_calls[0].args[2], include_attributes=False),
            current_git_expression,
        )

    def test_role_vectors_bind_role_commits_and_disable_injection_slots(self) -> None:
        argv = launch.role_configure_argv(
            Path("/sealed/harness"), Path("/sealed/current"),
            Path("/sealed/current-build/libwirehair.so.2.0.0"),
            Path("/sealed/current-role"), "current", "c" * 40,
            launch.CURRENT_IMPLEMENTATION_COMMIT,
        )
        for required in (
            "-DWIREHAIR_ROLE=current",
            "-DWIREHAIR_HARNESS_GIT_COMMIT=" + "c" * 40,
            "-DWIREHAIR_IMPLEMENTATION_GIT_COMMIT="
            + launch.CURRENT_IMPLEMENTATION_COMMIT,
            "-DCMAKE_CXX_FLAGS=", "-DCMAKE_EXE_LINKER_FLAGS=",
            "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF",
        ):
            self.assertIn(required, argv)


if __name__ == "__main__":
    unittest.main()
