#!/usr/bin/env python3
"""Offline regression tests for the preferred-attempt holdout controller."""

from __future__ import annotations

import argparse
from contextlib import ExitStack, redirect_stdout
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

import wh2_preferred_attempt_campaign as campaign_module
import wh2_preferred_attempt_search as subject
import wh2_rank_floor_two_anchor_screen as screen_module


class PreferredAttemptSearchTest(unittest.TestCase):
    def fake_repo(self, parent: Path) -> tuple[Path, str]:
        repo = parent / "source"
        remote = parent / "origin.git"
        subprocess.run(("git", "init", "-q", "--bare", str(remote)), check=True)
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
        subprocess.run(
            ("git", "-C", str(repo), "commit", "--allow-empty", "-qm", "fixture"),
            check=True,
        )
        subprocess.run(
            ("git", "-C", str(repo), "remote", "add", "origin", str(remote)),
            check=True,
        )
        commit = subprocess.check_output(
            ("git", "-C", str(repo), "rev-parse", "HEAD"), text=True,
        ).strip()
        return repo, commit

    def fake_contract(self, repo: Path, commit: str) -> dict[str, object]:
        tools: dict[str, dict[str, str]] = {}
        for name in ("git", "ssh", "ssh-add"):
            executable = Path(shutil.which(name) or f"/usr/bin/{name}").resolve(
                strict=True)
            tools[name] = {
                "path": str(executable),
                "sha256": subject.common.sha256_file(executable),
            }
        return {
            "source_repo": str(repo),
            "source_commit": commit,
            "protocol_sha256": "1" * 64,
            "binary_sha256": "2" * 64,
            "drand_wrapper_sha256": "3" * 64,
            "drand_client_bundle_sha256": "4" * 64,
            "drand_package_sha256": "5" * 64,
            "drand_package_lock_sha256": "6" * 64,
            "node": {
                "path": "/usr/bin/node", "version": subject.DRAND_NODE_VERSION,
                "sha256": "7" * 64,
            },
            "npm": {
                "path": "/usr/bin/npm", "version": subject.DRAND_NPM_VERSION,
                "sha256": "8" * 64,
            },
            "system_tools": tools,
        }

    def required_h1_tree(self, result: Path) -> None:
        for relative in subject.HOLDOUT_PHASES["h1"].required_files:
            path = result / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((relative + "\n").encode("ascii"))

    def latest_bound(self, latest_round: int = 1) -> dict[str, object]:
        origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
        return {
            "schema": "wirehair.wh2.verified_latest_bound.v1",
            "latest_round": latest_round,
            "canonical_beacon_sha256": "9" * 64,
            "consensus_origins": origins[:3],
            "raw_response_sha256": {
                origin: f"{index + 1:x}" * 64
                for index, origin in enumerate(origins[:3])
            },
            "wrapper_result_sha256": "a" * 64,
            "verified_ms": subject.now_utc_ms(),
        }

    def complete_prepare_staging(
        self,
        staging: Path,
        final: Path,
        repo: Path,
        commit: str,
        remote: Path,
    ) -> tuple[dict[str, object], dict[str, object]]:
        frozen = staging / "frozen"
        frozen.mkdir(parents=True)
        tools = self.fake_contract(repo, commit)["system_tools"]
        bindings = {
            "source_commit": commit,
            "binary_sha256": "1" * 64,
            "script_sha256": "2" * 64,
            "allk_sha256": "3" * 64,
            "helper_sha256": "4" * 64,
            "campaign_module_sha256": "5" * 64,
            "holdout_module_sha256": "6" * 64,
            "timing_module_sha256": "7" * 64,
            "protocol_sha256": "8" * 64,
            "route_context_sha256": "9" * 64,
            "q0_identity_sha256": "",
            "drand_wrapper_sha256": "a" * 64,
            "drand_package_sha256": "b" * 64,
            "drand_package_lock_sha256": "c" * 64,
            "drand_client_bundle_sha256": "d" * 64,
            "node": {"path": "/node", "version": "test", "sha256": "e" * 64},
            "npm": {"path": "/npm", "version": "test", "sha256": "f" * 64},
            "python": {"path": "/python", "version": "test", "sha256": "0" * 64},
            "system_tools": tools,
            "host_runtime": {"boot_id": "test"},
            "drand_offline_selftest": {"passed": True},
        }
        q0_identity = {"passed": True, "comparisons": 1}
        subject.common.atomic_json(frozen / "q0_identity.json", q0_identity)
        bindings["q0_identity_sha256"] = subject.common.sha256_file(
            frozen / "q0_identity.json")
        contract = {
            "schema": subject.SCHEMA,
            "source_repo": str(repo),
            "holdout_root_files_present": False,
            **bindings,
        }
        subject.common.atomic_json(frozen / "contract.json", contract)
        contract_sha256 = subject.common.sha256_file(frozen / "contract.json")
        subject.common.atomic_write(
            frozen / "staged.sha256",
            subject.common.sha_manifest(
                staging, tuple(frozen.iterdir())),
        )
        prepare = {
            "schema": "wirehair.wh2.h12_preferred_attempt.prepare.v1",
            **bindings,
            "q0_identity": q0_identity,
            "contract_sha256": contract_sha256,
            "staged_seal_sha256": subject.common.sha256_file(
                frozen / "staged.sha256"),
            "development_seed_ledger_sha256": "1" * 64,
            "future_beacon_contract_sha256": "2" * 64,
            "github_remote": str(remote),
            "result_dir": str(final),
            "launch_state": "CLEAN_UNLAUNCHED",
            "holdout_root_files_present": False,
            "prepared_utc_ns": time.time_ns(),
        }
        subject.common.atomic_json(staging / "prepare.json", prepare)
        return prepare, contract

    def fake_timing_recovery_tools(
        self,
        root: Path,
    ) -> dict[str, dict[str, str]]:
        tools: dict[str, dict[str, str]] = {}
        for name in subject.TIMING_HOST_RECOVERY_TOOLS:
            path = root / name
            path.write_bytes((f"mock {name}\n").encode("ascii"))
            path.chmod(0o755)
            tools[name] = {
                "path": str(path.resolve(strict=True)),
                "sha256": subject.common.sha256_file(path),
            }
        return tools

    def fake_timing_host_journal(
        self,
        root: Path,
        boot_id: str,
    ) -> dict[str, object]:
        return {
            "schema": subject.TIMING_HOST_JOURNAL_SCHEMA,
            "boot_id": boot_id,
            "core": 0,
            "siblings": [{"cpu": 64, "online": "1"}],
            "slices": [
                {"unit": "system.slice", "allowed": "",
                 "effective": "0 64"},
                {"unit": "user.slice", "allowed": "",
                 "effective": "0 64"},
                {"unit": "machine.slice", "allowed": "",
                 "effective": ""},
            ],
            "fillers": {
                "command": subject.TimingHostSession._filler_command_record(),
                "count": 2,
                "cpus": [0, 1],
            },
            "recovery_tools": self.fake_timing_recovery_tools(root),
        }

    def test_atomic_write_recovers_dead_partials_and_defers_sigterm(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-atomic-write-") as temporary:
            root = Path(temporary)
            target = root / "wave.json"
            stale = root / ".wave.json.99999999.partial"
            stale.write_bytes(b"stale\n")
            subject.common.atomic_write(target, b"complete\n")
            self.assertEqual(target.read_bytes(), b"complete\n")
            self.assertFalse(stale.exists())

            published = root / "published.json"
            linked_marker = root / ".published.json.99999999.partial"
            published.write_bytes(b"immutable\n")
            os.link(published, linked_marker)
            self.assertEqual(published.stat().st_nlink, 2)
            subject.common.atomic_write_once_or_same(
                published, b"immutable\n")
            self.assertFalse(linked_marker.exists())
            self.assertEqual(published.stat().st_nlink, 1)

            readable = root / "readable.json"
            readable_marker = root / ".readable.json.99999999.partial"
            readable.write_bytes(b"readable\n")
            os.link(readable, readable_marker)
            self.assertEqual(
                subject.common.stable_bytes(readable), b"readable\n")
            self.assertFalse(readable_marker.exists())
            self.assertEqual(readable.stat().st_nlink, 1)

            os.link(readable, readable_marker)
            self.assertEqual(
                subject.common.sha256_file(readable),
                subject.common.sha256_bytes(b"readable\n"))
            self.assertFalse(readable_marker.exists())

            os.link(readable, readable_marker)
            self.assertEqual(
                subject.common.sha256_file(
                    readable, require_unique=False),
                subject.common.sha256_bytes(b"readable\n"))
            self.assertTrue(readable_marker.exists())
            self.assertEqual(readable.stat().st_nlink, 2)
            self.assertEqual(
                subject.common.stable_bytes(readable), b"readable\n")
            self.assertFalse(readable_marker.exists())

            malicious = root / "malicious.json"
            malicious_alias = root / "malicious-alias.json"
            malicious.write_bytes(b"malicious\n")
            os.link(malicious, malicious_alias)
            with self.assertRaisesRegex(
                    subject.common.CampaignError, "nonunique"):
                subject.common.stable_bytes(malicious)
            with self.assertRaisesRegex(
                    subject.common.CampaignError, "nonunique"):
                subject.common.sha256_file(malicious)
            self.assertEqual(
                subject.common.sha256_file(
                    malicious, require_unique=False),
                subject.common.sha256_bytes(b"malicious\n"))
            self.assertTrue(malicious.exists())
            self.assertTrue(malicious_alias.exists())

            racy = root / "racy.json"
            racy.write_bytes(b"old\n")
            real_stat = subject.common.os.stat
            changed = False

            def mutate_before_named_stat(*args, **kwargs):
                nonlocal changed
                if (not changed and kwargs.get("dir_fd") is not None and
                        args and args[0] == racy.name):
                    changed = True
                    racy.write_bytes(b"new\n")
                return real_stat(*args, **kwargs)

            with mock.patch.object(
                    subject.common.os, "stat",
                    side_effect=mutate_before_named_stat), \
                 self.assertRaisesRegex(
                     subject.common.CampaignError, "changed"):
                subject.common.stable_bytes(racy)
            self.assertTrue(changed)

            exclusive = root / "exclusive-directory"
            exclusive_fd = subject.common.create_durable_directory(exclusive)
            os.close(exclusive_fd)
            exclusive_inode = exclusive.stat().st_ino
            with self.assertRaisesRegex(
                    subject.common.CampaignError, "already exists"):
                subject.common.create_durable_directory(exclusive)
            self.assertEqual(exclusive.stat().st_ino, exclusive_inode)

            fault_parent = root / "fault-parent"
            fault_parent.mkdir()
            parent_fd = os.open(
                fault_parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
            opened_child: list[int] = []
            real_open = subject.common.os.open
            real_close = subject.common.os.close
            injected = False

            def track_open(*args: object, **kwargs: object) -> int:
                descriptor = real_open(*args, **kwargs)
                opened_child.append(descriptor)
                return descriptor

            def fail_parent_close(descriptor: int) -> None:
                nonlocal injected
                if descriptor == parent_fd and not injected:
                    injected = True
                    real_close(descriptor)
                    raise OSError("injected parent close failure")
                real_close(descriptor)

            with mock.patch.object(
                    subject.common, "open_durable_directory",
                    return_value=parent_fd), \
                 mock.patch.object(
                     subject.common.os, "open", side_effect=track_open), \
                 mock.patch.object(
                     subject.common.os, "close", side_effect=fail_parent_close), \
                 self.assertRaisesRegex(OSError, "parent close failure"):
                subject.common.create_durable_directory(
                    fault_parent / "created")
            self.assertTrue(injected)
            self.assertEqual(len(opened_child), 1)
            for descriptor in (parent_fd, opened_child[0]):
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

            outside = root / "outside"
            outside.mkdir()
            alias = root / "directory-alias"
            alias.symlink_to(outside, target_is_directory=True)
            with self.assertRaises(subject.CampaignError):
                subject.common.atomic_write_once_or_same(
                    alias / "escaped.json", b"must not escape\n")
            self.assertFalse((outside / "escaped.json").exists())

            live = root / f".wave.json.{os.getpid()}.partial"
            live.write_bytes(b"active\n")
            with self.assertRaisesRegex(
                    subject.common.CampaignError, "writer is still active"):
                subject.common.atomic_write(target, b"replacement\n")
            self.assertEqual(live.read_bytes(), b"active\n")
            live.unlink()

            start_ticks = subject.common._process_start_ticks(os.getpid())
            self.assertIsNotNone(start_ticks)
            reused = root / (
                f".wave.json.{os.getpid()}.{int(start_ticks) + 1}.partial")
            reused.write_bytes(b"stale reused pid\n")
            subject.common.atomic_write(target, b"replacement\n")
            self.assertFalse(reused.exists())
            self.assertEqual(target.read_bytes(), b"replacement\n")

            live_identity = root / (
                f".wave.json.{os.getpid()}.{start_ticks}.partial")
            live_identity.write_bytes(b"active identity\n")
            with self.assertRaisesRegex(
                    subject.common.CampaignError, "writer is still active"):
                subject.common.atomic_write(target, b"must not replace\n")
            self.assertTrue(live_identity.exists())
            live_identity.unlink()

            signaled_target = root / "signaled.json"
            ready = root / "ready"
            script = "\n".join((
                "import pathlib, sys, time",
                f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
                "import wh2_rank_floor_two_anchor_allk as common",
                "original_fsync = common.os.fsync",
                "first = True",
                "def slow_fsync(fd):",
                "    global first",
                "    if first:",
                "        first = False",
                f"        pathlib.Path({str(ready)!r}).write_bytes(b'ready')",
                "        time.sleep(1.0)",
                "    original_fsync(fd)",
                "common.os.fsync = slow_fsync",
                f"common.atomic_write(pathlib.Path({str(signaled_target)!r}), b'x' * 1024)",
            ))
            process = subprocess.Popen((sys.executable, "-c", script))
            for _ in range(500):
                if ready.exists():
                    break
                if process.poll() is not None:
                    break
                time.sleep(0.01)
            self.assertTrue(ready.exists())
            process.terminate()
            self.assertEqual(process.wait(timeout=5), -signal.SIGTERM)
            self.assertEqual(signaled_target.read_bytes(), b"x" * 1024)
            self.assertEqual(
                list(root.glob(".signaled.json.*.partial")), [])

            if hasattr(signal, "pthread_sigmask"):
                before = signal.pthread_sigmask(signal.SIG_BLOCK, set())
                cleanup_target = root / "cleanup-failure.json"
                try:
                    with mock.patch.object(
                            subject.common.os, "replace",
                            side_effect=OSError("replace failed")), \
                         mock.patch.object(
                             subject.common.os, "unlink",
                             side_effect=OSError("cleanup failed")), \
                         self.assertRaisesRegex(OSError, "cleanup failed"):
                        subject.common.atomic_write(cleanup_target, b"value\n")
                    after = signal.pthread_sigmask(signal.SIG_BLOCK, set())
                    self.assertEqual(after, before)
                finally:
                    signal.pthread_sigmask(signal.SIG_SETMASK, before)
                    for partial in root.glob(
                            ".cleanup-failure.json.*.partial"):
                        partial.unlink()

    def test_development_analysis_inventory_discards_dead_manifest_partial(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-analysis-inventory-") as temporary:
            development = Path(temporary)
            artifact = development / "summary.json"
            artifact.write_bytes(b"sealed\n")
            phase = development / "phase_complete.sha256"
            phase.write_bytes(b"phase\n")
            stale = development / \
                ".analysis_complete.sha256.99999999.partial"
            stale.write_bytes(b"interrupted\n")
            self.assertEqual(
                subject.development_analysis_inventory(
                    development, {artifact, phase}),
                {artifact.resolve(strict=True), phase.resolve(strict=True)})
            self.assertFalse(stale.exists())
            alias = development / "summary-alias.json"
            alias.symlink_to(artifact.name)
            with self.assertRaisesRegex(subject.CampaignError, "symlink"):
                subject.development_analysis_inventory(
                    development, {artifact, phase})
            alias.unlink()
            junk = development / "junk.txt"
            junk.write_bytes(b"ambient\n")
            with self.assertRaisesRegex(
                    subject.CampaignError, "artifact set mismatch"):
                subject.development_analysis_inventory(
                    development, {artifact, phase})

    def test_campaign_execution_lock_is_persistent_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-command-lock-") as temporary:
            result = Path(temporary).resolve()
            with subject.campaign_execution_lock(result):
                lock = result / "locks/campaign-execution.lock"
                self.assertTrue(lock.is_file())
                with self.assertRaisesRegex(
                        subject.CampaignError, "controller is active"):
                    with subject.campaign_execution_lock(result):
                        self.fail("a second controller acquired the lock")
            self.assertTrue(lock.is_file())
            with subject.campaign_execution_lock(result):
                pass

    def test_immutable_publication_never_overwrites_racing_target(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-publish-race-") as temporary:
            target = Path(temporary) / "artifact.json"
            real_link = subject.common.os.link

            def competing_link(source: object, destination: object,
                               **kwargs: object) -> None:
                target.write_bytes(b"competitor\n")
                real_link(source, destination, **kwargs)

            with mock.patch.object(
                    subject.common.os, "link", side_effect=competing_link), \
                 self.assertRaises(subject.CampaignError):
                campaign_module.durable_write_once_or_same(
                    target, b"intended\n")
            self.assertEqual(b"competitor\n", target.read_bytes())

    def test_expected_regular_file_readers_do_not_block_on_fifo(self) -> None:
        script = "\n".join((
            "import os, sys, tempfile",
            "from pathlib import Path",
            f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
            "import wh2_preferred_attempt_search as search",
            "import wh2_rank_floor_two_anchor_allk as common",
            "import wh2_rank_floor_two_anchor_screen as screen",
            "with tempfile.TemporaryDirectory() as temporary:",
            "    fifo = Path(temporary) / 'artifact'",
            "    os.mkfifo(fifo)",
            "    for reader in (screen.sha256_file, "
            "                   screen.stable_file_bytes, "
            "                   common.stable_bytes, "
            "                   search.fsync_existing_regular_file):",
            "        try:",
            "            reader(fifo)",
            "        except (ValueError, common.CampaignError, "
            "                search.CampaignError):",
            "            pass",
            "        else:",
            "            raise AssertionError(reader)",
            "    try:",
            "        screen.thermal_start(fifo)",
            "    except ValueError:",
            "        pass",
            "    else:",
            "        raise AssertionError('thermal_start accepted FIFO')",
            "    try:",
            "        screen.thermal_finish(",
            "            fifo, {'dev': 0, 'ino': 0, 'offset': 0},",
            "            fifo.with_name('sealed.csv'))",
            "    except ValueError:",
            "        pass",
            "    else:",
            "        raise AssertionError('thermal_finish accepted FIFO')",
            "    guard = object.__new__(common.ThermalGuard)",
            "    guard.path = fifo",
            "    guard.mark = {'dev': 0, 'ino': 0}",
            "    guard.read_offset = 0",
            "    try:",
            "        guard._new_samples()",
            "    except ValueError:",
            "        pass",
            "    else:",
            "        raise AssertionError('thermal guard accepted FIFO')",
        ))
        completed = subprocess.run(
            (sys.executable, "-c", script), capture_output=True,
            text=True, check=False, timeout=5)
        self.assertEqual(
            completed.returncode, 0,
            msg=completed.stdout + completed.stderr)

    def test_frozen_module_imports_never_create_bytecode_entries(self) -> None:
        module_names = (
            "wh2_preferred_attempt_campaign",
            "wh2_preferred_attempt_holdout",
            "wh2_preferred_attempt_search",
            "wh2_preferred_attempt_timing",
            "wh2_rank_floor_two_anchor_allk",
            "wh2_rank_floor_two_anchor_screen",
        )
        source = Path(__file__).parent
        with tempfile.TemporaryDirectory(
                prefix="wh2-frozen-no-pycache-") as temporary:
            frozen = Path(temporary)
            for module_name in module_names:
                shutil.copyfile(
                    source / f"{module_name}.py",
                    frozen / f"{module_name}.py")
            environment = os.environ.copy()
            environment.pop("PYTHONDONTWRITEBYTECODE", None)
            script = (
                "import pathlib, runpy, sys; "
                "path=pathlib.Path(sys.argv[1]); "
                "sys.path.insert(0, str(path.parent)); "
                "runpy.run_path(str(path), run_name='frozen_probe'); "
                "assert not (path.parent/'__pycache__').exists()"
            )
            for module_name in module_names:
                with self.subTest(module=module_name):
                    shutil.rmtree(frozen / "__pycache__", ignore_errors=True)
                    completed = subprocess.run(
                        (sys.executable, "-c", script,
                         str(frozen / f"{module_name}.py")),
                        env=environment, capture_output=True, text=True,
                        check=False, timeout=10)
                    self.assertEqual(
                        completed.returncode, 0,
                        msg=completed.stdout + completed.stderr)

    def test_prepare_transaction_sigterm_cleanup_and_persistent_flock(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-prepare-transaction-signal-") as temporary:
            root = Path(temporary)
            original_handlers = {
                signum: signal.getsignal(signum)
                for signum in (signal.SIGINT, signal.SIGTERM)
            }
            final = root / "campaign"
            ready = root / "ready"
            script = "\n".join((
                "from pathlib import Path",
                "import signal, sys, time",
                f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
                "import wh2_preferred_attempt_search as subject",
                "signal.signal(signal.SIGINT, signal.default_int_handler)",
                "final, ready = Path(sys.argv[1]), Path(sys.argv[2])",
                "with subject.atomic_result_directory(final) as transaction:",
                "    ready.write_text(str(transaction[0]), encoding='ascii')",
                "    time.sleep(30)",
            ))
            process = subprocess.Popen(
                (sys.executable, "-c", script, str(final), str(ready)))
            for _ in range(500):
                if ready.exists() or process.poll() is not None:
                    break
                time.sleep(0.01)
            self.assertTrue(ready.exists())
            with self.assertRaisesRegex(
                    subject.CampaignError, "lock is held"):
                with subject.atomic_result_directory(final):
                    self.fail("contending prepare unexpectedly acquired the lock")
            process.terminate()
            self.assertEqual(process.wait(timeout=5), -signal.SIGTERM)
            lock = root / ".campaign.prepare.lock"
            self.assertTrue(lock.is_file())
            lock_inode = lock.stat().st_ino
            self.assertEqual(
                [path for path in root.iterdir()
                 if path.is_dir() and path.name.startswith(".campaign.prepare")],
                [],
            )
            self.assertEqual(
                {signum: signal.getsignal(signum)
                 for signum in (signal.SIGINT, signal.SIGTERM)},
                original_handlers,
            )

            interrupt_final = root / "interrupt-campaign"
            interrupt_ready = root / "interrupt-ready"
            interrupted = subprocess.Popen((
                sys.executable, "-c", script, str(interrupt_final),
                str(interrupt_ready),
            ), stderr=subprocess.DEVNULL)
            for _ in range(500):
                if interrupt_ready.exists() or interrupted.poll() is not None:
                    break
                time.sleep(0.01)
            self.assertTrue(interrupt_ready.exists())
            interrupted.send_signal(signal.SIGINT)
            self.assertEqual(interrupted.wait(timeout=5), -signal.SIGINT)
            self.assertEqual(
                [path for path in root.iterdir() if path.is_dir() and
                 path.name.startswith(".interrupt-campaign.prepare")],
                [],
            )

            legacy_final = root / "legacy-campaign"
            legacy_lock = root / ".legacy-campaign.prepare.lock"
            legacy_lock.write_text("pid=99999999\n", encoding="ascii")
            legacy_staging = root / ".legacy-campaign.prepare-99999999"
            legacy_staging.mkdir()
            with self.assertRaisesRegex(RuntimeError, "fixture stop"):
                with subject.atomic_result_directory(legacy_final):
                    self.assertFalse(legacy_staging.exists())
                    raise RuntimeError("fixture stop")
            self.assertTrue(legacy_lock.read_text(encoding="ascii").startswith(
                "schema=wirehair.wh2.prepare_lock.v1\n"))

            live_final = root / "live-legacy-campaign"
            live_lock = root / ".live-legacy-campaign.prepare.lock"
            live_lock.write_text(f"pid={os.getpid()}\n", encoding="ascii")
            with self.assertRaisesRegex(
                    subject.CampaignError, "legacy prepare owner is still active"):
                with subject.atomic_result_directory(live_final):
                    self.fail("live legacy owner was ignored")

            returning_final = root / "returning-handler-campaign"
            previous_term = signal.signal(signal.SIGTERM, lambda *_args: None)

            def returning_delivery(_pid: int, signum: int) -> None:
                handler = signal.getsignal(signum)
                self.assertTrue(callable(handler))
                handler(signum, None)

            try:
                with mock.patch.object(
                        subject.os, "kill", side_effect=returning_delivery), \
                     self.assertRaisesRegex(
                         subject.CampaignError,
                         "handler returned without terminating"):
                    with subject.atomic_result_directory(returning_final):
                        signal.raise_signal(signal.SIGTERM)
                    self.fail("interrupted prepare continued after its body")
            finally:
                signal.signal(signal.SIGTERM, previous_term)

            stale = root / ".campaign.prepare-99999999"
            stale.mkdir()
            (stale / "incomplete").write_bytes(b"discard me\n")
            with self.assertRaisesRegex(RuntimeError, "fixture stop"):
                with subject.atomic_result_directory(
                        final, create=False) as transaction:
                    self.assertFalse(stale.exists())
                    self.assertIsNone(transaction[2])
                    raise RuntimeError("fixture stop")
            self.assertEqual(lock.stat().st_ino, lock_inode)
            self.assertEqual(
                [path for path in root.iterdir()
                 if path.is_dir() and path.name.startswith(".campaign.prepare")],
                [],
            )
            self.assertEqual(
                {signum: signal.getsignal(signum)
                 for signum in (signal.SIGINT, signal.SIGTERM)},
                original_handlers,
            )

    def test_prepare_transaction_reconciles_complete_staging(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-prepare-transaction-recover-") as temporary:
            root = Path(temporary).resolve()
            repo, commit = self.fake_repo(root)
            remote = root / "origin.git"
            final = root / "campaign"
            staging = root / ".campaign.prepare-99999999"
            staging.mkdir()
            prepare, contract = self.complete_prepare_staging(
                staging, final, repo, commit, remote)
            environment = os.environ.copy()
            with mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=environment), \
                 mock.patch.object(
                     subject, "verify_recoverable_prepare_sources"):
                prepare_sha256 = subject.common.sha256_file(
                    staging / "prepare.json")
                published_before_crash = subject.publish_r1_freeze_tag(
                    repo, prepare, prepare_sha256,
                    contract["system_tools"])
                self.assertFalse(
                    (staging / "frozen/r1_freeze_publication.json").exists())
                with subject.atomic_result_directory(final) as transaction:
                    self.assertEqual(transaction[0], final)
                    self.assertEqual(transaction[1], final)
                    self.assertEqual(transaction[2], prepare)
                subject.verify_r1_freeze_publication(
                    final, prepare, contract)
            self.assertFalse(staging.exists())
            self.assertTrue(final.is_dir())
            self.assertTrue(
                (final / "frozen/r1_freeze_publication.json").is_file())
            self.assertEqual(subject.load_fixed_json(
                final / "frozen/r1_freeze_publication.json",
                "recovered R1 receipt")["tag_object"],
                published_before_crash["tag_object"])

            absent = root / "absent-campaign"
            with subject.atomic_result_directory(
                    absent, create=False) as transaction:
                self.assertEqual(transaction, (absent, absent, None))
                self.assertFalse(absent.exists())
                self.assertFalse(
                    (root / ".absent-campaign.prepare.partial").exists())

            broken_final = root / "broken-campaign"
            broken_staging = root / ".broken-campaign.prepare.partial"
            broken_staging.mkdir()
            self.complete_prepare_staging(
                broken_staging, broken_final, repo, commit, remote)
            subject.common.atomic_json(
                broken_staging / "frozen/r1_freeze_publication.json",
                {"status": "malformed"},
            )
            with mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=environment), \
                 mock.patch.object(
                     subject, "verify_recoverable_prepare_sources"), \
                 self.assertRaises(
                        subject.CampaignError):
                with subject.atomic_result_directory(broken_final):
                    self.fail("malformed complete staging was accepted")
            self.assertTrue(broken_staging.is_dir())
            self.assertTrue((broken_staging / "prepare.json").is_file())
            self.assertFalse(broken_final.exists())

    def test_prepare_sigterm_finalizes_complete_staging(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-prepare-complete-signal-") as temporary:
            root = Path(temporary).resolve()
            repo, commit = self.fake_repo(root)
            remote = root / "origin.git"
            final = root / "campaign"
            template = root / "template"
            template.mkdir()
            prepare, contract = self.complete_prepare_staging(
                template, final, repo, commit, remote)
            ready = root / "ready"
            script = "\n".join((
                "from pathlib import Path",
                "import os, shutil, sys, time",
                f"sys.path.insert(0, {str(Path(__file__).parent)!r})",
                "import wh2_preferred_attempt_search as subject",
                "final, template, ready = map(Path, sys.argv[1:])",
                "subject.newest_github_agent_environment = "
                "lambda _tools: os.environ.copy()",
                "subject.verify_recoverable_prepare_sources = "
                "lambda _frozen, _contract: None",
                "with subject.atomic_result_directory(final) as transaction:",
                "    staging = transaction[0]",
                "    shutil.copytree(template / 'frozen', staging / 'frozen')",
                "    shutil.copyfile(template / 'prepare.json', "
                "staging / 'prepare.json')",
                "    ready.write_bytes(b'ready')",
                "    time.sleep(30)",
            ))
            process = subprocess.Popen((
                sys.executable, "-c", script, str(final), str(template),
                str(ready),
            ))
            for _ in range(500):
                if ready.exists() or process.poll() is not None:
                    break
                time.sleep(0.01)
            self.assertTrue(ready.exists())
            process.terminate()
            self.assertEqual(process.wait(timeout=10), -signal.SIGTERM)
            self.assertTrue(final.is_dir())
            self.assertEqual(
                [path for path in root.iterdir()
                 if path.is_dir() and path.name.startswith(".campaign.prepare")],
                [],
            )
            self.assertTrue(
                (final / "frozen/r1_freeze_publication.json").is_file())
            with mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=os.environ.copy()):
                subject.verify_r1_freeze_publication(
                    final, prepare, contract)

    def test_quicknet_formula_and_root_goldens(self) -> None:
        subject.validate_beacon_contract_kats()
        canonical = subject.canonical_beacon_bytes(subject.DRAND_KNOWN_ROUND)
        self.assertEqual(
            hashlib.sha256(canonical).hexdigest(),
            subject.DRAND_CANONICAL_KNOWN_SHA256,
        )
        self.assertEqual(
            subject.holdout_master_root("0" * 64, subject.DRAND_KNOWN_ROUND),
            subject.DRAND_MASTER_ROOT_KNOWN_SHA256,
        )
        mutated = dict(subject.DRAND_KNOWN_ROUND)
        mutated["randomness"] = "0" * 64
        with self.assertRaises(subject.CampaignError):
            subject.canonical_beacon_bytes(mutated)
        extended = dict(subject.DRAND_KNOWN_ROUND)
        extended["unexpected"] = "ignored entropy"
        with self.assertRaises(subject.CampaignError):
            subject.canonical_beacon_bytes(extended)
        seal_ms = subject.DRAND_GENESIS_TIME * 1000 + 1_000_000
        local_round = subject.target_beacon_round(seal_ms)
        latest_round = local_round + 100
        record, _encoded = subject.build_seal_record(
            "h1", "b" * 64,
            self.fake_contract(Path("/tmp"), "c" * 40),
            "git@github.com:fake/wirehair.git", seal_ms,
            self.latest_bound(latest_round),
        )
        self.assertEqual(record["round"], latest_round + 40)
        reseal, _encoded = subject.build_seal_record(
            "h1", "b" * 64,
            self.fake_contract(Path("/tmp"), "c" * 40),
            "git@github.com:fake/wirehair.git", seal_ms,
            self.latest_bound(latest_round), record["round"],
        )
        self.assertEqual(reseal["round"], record["round"] + 1)

    def test_verified_latest_is_only_a_hashed_round_bound(self) -> None:
        beacon = dict(subject.DRAND_KNOWN_ROUND)
        beacon["round"] = 123456
        origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
        canonical_sha256 = hashlib.sha256(
            subject.canonical_beacon_bytes(beacon)).hexdigest()
        observations = [{
            "origin": origin,
            "ok": True,
            "beacon": beacon,
            "canonical_sha256": canonical_sha256,
            "raw_response_sha256": f"{index + 1:x}" * 64,
        } for index, origin in enumerate(origins)]
        wave = {
            "schema": "wirehair.wh2.drand_quicknet_latest_wave.v1",
            "role": "preseal-clock-liveness-lower-bound-only",
            "chain_hash": subject.DRAND_CHAIN_HASH,
            "status": "QUORUM",
            "latest_round": beacon["round"],
            "beacon": beacon,
            "consensus_origins": origins[:3],
            "observations": observations,
        }
        completed = subprocess.CompletedProcess(
            ("node", "wrapper", "bundle", "latest-bound"), 0,
            json.dumps(wave).encode("ascii"), b"",
        )
        contract = {"node": {"path": "/usr/bin/node"}}
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(
                subject.subprocess, "run", return_value=completed))
            stack.enter_context(mock.patch.object(
                subject, "offline_verify_beacon", return_value={}))
            bound = subject.obtain_verified_latest_bound(
                contract, Path("/does/not/matter"))
        self.assertEqual(bound["latest_round"], beacon["round"])
        self.assertNotIn("beacon", bound)
        self.assertNotIn("randomness", json.dumps(bound))
        no_quorum = dict(wave)
        no_quorum["status"] = "TEMPORARY_NO_QUORUM"
        failed = subprocess.CompletedProcess(
            completed.args, 2, json.dumps(no_quorum).encode("ascii"), b"")
        with mock.patch.object(subject.subprocess, "run", return_value=failed):
            with self.assertRaises(subject.CampaignError):
                subject.obtain_verified_latest_bound(
                    contract, Path("/does/not/matter"))
        malformed = subprocess.CompletedProcess(
            completed.args, 0, b"\xff", b"")
        with mock.patch.object(subject.subprocess, "run", return_value=malformed):
            with self.assertRaises(subject.CampaignError):
                subject.obtain_verified_latest_bound(
                    contract, Path("/does/not/matter"))

    def test_agent_refresh_uses_frozen_absolute_ssh_tools(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="ssh-wh2-agent-") as temporary:
            root = Path(temporary)
            repo, commit = self.fake_repo(root)
            tools = self.fake_contract(repo, commit)["system_tools"]
            socket_path = root / "agent.999999"
            listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            listener.bind(str(socket_path))
            try:
                listed = subprocess.CompletedProcess(
                    (), 0, "2048 SHA256:test forwarded-key\n", "")
                authenticated = subprocess.CompletedProcess(
                    (), 1, "", "Hi! You've successfully authenticated.\n")
                with mock.patch.object(
                        subject.Path, "glob", return_value=[socket_path]), \
                     mock.patch.object(
                         subject.subprocess, "run",
                         side_effect=(listed, authenticated)) as run:
                    environment = subject.newest_github_agent_environment(tools)
                ssh = str(tools["ssh"]["path"])
                ssh_add = str(tools["ssh-add"]["path"])
                self.assertEqual(run.call_count, 2)
                self.assertEqual(run.call_args_list[0].args[0][0], ssh_add)
                self.assertEqual(run.call_args_list[1].args[0][0], ssh)
                self.assertEqual(environment["SSH_AUTH_SOCK"], str(socket_path))
                self.assertEqual(
                    environment["GIT_SSH_COMMAND"], subject.shlex.join((
                        ssh, "-o", "BatchMode=yes", "-o",
                        "ConnectTimeout=10")))
                self.assertEqual(environment["GIT_SSH_VARIANT"], "ssh")
                self.assertNotIn("GIT_SSH", environment)
            finally:
                listener.close()

    def test_frozen_executable_records_allow_hardlinks_not_indirection(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-frozen-executable-") as temporary:
            root = Path(temporary).resolve()
            executable = root / "runtime"
            executable.write_bytes(b"mock frozen executable\n")
            executable.chmod(0o755)
            hardlink = root / "runtime-hardlink"
            os.link(executable, hardlink)
            self.assertEqual(executable.stat().st_nlink, 2)
            digest = subject.common.sha256_file(
                executable, require_unique=False)
            tool_record = {"path": str(executable), "sha256": digest}
            runtime_record = {
                **tool_record, "version": subject.DRAND_NODE_VERSION,
            }

            self.assertEqual(
                subject.frozen_tool_path({"taskset": tool_record}, "taskset"),
                executable)
            session = object.__new__(subject.TimingHostSession)
            session.tools = {"taskset": tool_record}
            self.assertEqual(session.tool("taskset"), str(executable))
            with mock.patch.object(
                    subject, "command_version",
                    return_value=subject.DRAND_NODE_VERSION):
                self.assertEqual(subject.frozen_runtime_path(
                    runtime_record, "node", subject.DRAND_NODE_VERSION),
                    executable)
            self.assertTrue(hardlink.samefile(executable))
            self.assertEqual(executable.stat().st_nlink, 2)

            indirect = root / "runtime-symlink"
            indirect.symlink_to(executable.name)
            indirect_tool = {**tool_record, "path": str(indirect)}
            indirect_runtime = {**runtime_record, "path": str(indirect)}
            with self.assertRaises(subject.CampaignError):
                subject.frozen_tool_path(
                    {"taskset": indirect_tool}, "taskset")
            with mock.patch.object(
                    subject, "command_version",
                    return_value=subject.DRAND_NODE_VERSION), \
                 self.assertRaises(subject.CampaignError):
                subject.frozen_runtime_path(
                    indirect_runtime, "node", subject.DRAND_NODE_VERSION)

            loop = root / "runtime-loop"
            loop.symlink_to(loop.name)
            loop_tool = {**tool_record, "path": str(loop)}
            loop_runtime = {**runtime_record, "path": str(loop)}
            with self.assertRaises(subject.CampaignError):
                subject.frozen_tool_path({"taskset": loop_tool}, "taskset")
            with self.assertRaises(subject.CampaignError):
                subject.frozen_runtime_path(
                    loop_runtime, "node", subject.DRAND_NODE_VERSION)

            malformed = {**runtime_record, "sha256": "not-a-digest"}
            with self.assertRaises(subject.CampaignError):
                subject.frozen_runtime_path(
                    malformed, "node", subject.DRAND_NODE_VERSION)

    def test_fake_remote_tag_and_late_ordering(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-tag-test-") as temporary:
            repo, commit = self.fake_repo(Path(temporary))
            contract = self.fake_contract(repo, commit)
            tools = contract["system_tools"]
            environment = os.environ.copy()
            poison = Path(temporary) / "poison-path"
            poison.mkdir()
            fake_git = poison / "git"
            fake_git.write_text("#!/bin/sh\nexit 97\n", encoding="ascii")
            fake_git.chmod(0o755)
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=environment))
                stack.enter_context(mock.patch.dict(os.environ, {
                    "PATH": str(poison) + os.pathsep +
                    environment.get("PATH", ""),
                }))
                seal = "1" * 64
                record = {
                    "round": 123,
                    "round_time_ms": subject.now_utc_ms() + 60_000,
                    "source_commit": commit,
                }
                publication = subject.publish_seal_tag(
                    repo, "h1", seal, "2" * 64, record, tools)
                repeated = subject.publish_seal_tag(
                    repo, "h1", seal, "2" * 64, record, tools)
                self.assertEqual(publication["status"], "CONFIRMED")
                self.assertEqual(publication["tag_object"], repeated["tag_object"])
                late_seal = "3" * 64
                late_record = {
                    "round": 124,
                    "round_time_ms": subject.now_utc_ms() - 1,
                    "source_commit": commit,
                }
                with self.assertRaises(subject.CampaignError):
                    subject.publish_seal_tag(
                        repo, "h2", late_seal, "4" * 64, late_record, tools)
                late_name = (
                    "wh2-h12-kpreferred-h2-seal-" + late_seal[:16]
                )
                _raw, rows = subject.tag_remote_rows(
                    repo, late_name, environment, tools)
                self.assertEqual(rows, {})
                self.assertEqual(subprocess.run(
                    (
                        str(tools["git"]["path"]), "-C", str(repo),
                        "show-ref", "--verify", "--quiet",
                        f"refs/tags/{late_name}",
                    ),
                    check=False,
                ).returncode, 0)

                freeze = Path(temporary) / "freeze"
                (freeze / "frozen").mkdir(parents=True)
                prepare = {
                    "source_commit": commit,
                    "binary_sha256": "1" * 64,
                    "script_sha256": "2" * 64,
                    "contract_sha256": "3" * 64,
                    "protocol_sha256": "4" * 64,
                    "route_context_sha256": "5" * 64,
                    "staged_seal_sha256": "6" * 64,
                    "development_seed_ledger_sha256": "7" * 64,
                    "future_beacon_contract_sha256": "8" * 64,
                    "github_remote": str(Path(temporary) / "origin.git"),
                }
                subject.common.atomic_json(freeze / "prepare.json", prepare)
                prepare_sha = subject.common.sha256_file(freeze / "prepare.json")
                freeze_publication = subject.publish_r1_freeze_tag(
                    repo, prepare, prepare_sha, tools)
                publication_path = freeze / "frozen/r1_freeze_publication.json"
                subject.atomic_state_json(publication_path, freeze_publication)
                freeze_contract = {
                    "source_repo": str(repo), "source_commit": commit,
                    "system_tools": tools,
                }
                subject.verify_r1_freeze_publication(
                    freeze, prepare, freeze_contract)

                git = str(tools["git"]["path"])
                probe_tag = "wh2-replace-ref-probe"
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "tag", "-a", probe_tag, commit, "-m", "wrong annotation",
                ), check=True)
                wrong_object = subprocess.check_output((
                    git, "--no-replace-objects", "-C", str(repo),
                    "rev-parse", f"{probe_tag}^{{tag}}",
                ), text=True).strip()
                good_object = freeze_publication["tag_object"]
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "replace", wrong_object, good_object,
                ), check=True)
                expected_ref = "refs/tags/{}".format(
                    freeze_publication["tag_name"])
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "update-ref", expected_ref, wrong_object,
                ), check=True)
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "push", "--force", "origin",
                    f"{expected_ref}:{expected_ref}",
                ), check=True, capture_output=True)
                forged_raw, _forged_rows = subject.tag_remote_rows(
                    repo, freeze_publication["tag_name"], environment, tools)
                forged = {
                    **freeze_publication,
                    "tag_object": wrong_object,
                    "remote_rows_ascii": forged_raw.decode("ascii"),
                    "remote_rows_sha256": subject.sha256_bytes(forged_raw),
                }
                subject.atomic_state_json(publication_path, forged)
                with self.assertRaisesRegex(
                        subject.CampaignError, "annotation changed"):
                    subject.verify_r1_freeze_publication(
                        freeze, prepare, freeze_contract)
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "update-ref", expected_ref, good_object,
                ), check=True)
                subprocess.run((
                    git, "--no-replace-objects", "-C", str(repo),
                    "push", "--force", "origin",
                    f"{expected_ref}:{expected_ref}",
                ), check=True, capture_output=True)
                subject.atomic_state_json(publication_path, freeze_publication)
                publication_bytes = publication_path.read_bytes()
                publication_path.unlink()
                with self.assertRaises(subject.CampaignError):
                    subject.verify_r1_freeze_publication(
                        freeze, prepare, freeze_contract)
                subject.common.atomic_write(publication_path, publication_bytes)
                mutated_prepare = dict(prepare)
                mutated_prepare["staged_seal_sha256"] = "9" * 64
                subject.common.atomic_json(
                    freeze / "prepare.json", mutated_prepare)
                with self.assertRaises(subject.CampaignError):
                    subject.verify_r1_freeze_publication(
                        freeze, mutated_prepare, freeze_contract)
                mutated_sha = subject.common.sha256_file(freeze / "prepare.json")
                replacement_publication = subject.publish_r1_freeze_tag(
                    repo, mutated_prepare, mutated_sha, tools)
                subject.atomic_state_json(
                    publication_path, replacement_publication)
                subject.verify_r1_freeze_publication(
                    freeze, mutated_prepare, freeze_contract)

    def test_staged_seal_is_immutable_root_of_trust(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-staged-anchor-") as temporary:
            result = Path(temporary)
            frozen = result / "frozen"
            frozen.mkdir()
            subject.common.atomic_json(frozen / "contract.json", {"version": 1})
            subject.common.atomic_write(frozen / "binary", b"original\n")
            staged_files = [frozen / "contract.json", frozen / "binary"]
            subject.common.atomic_write(
                frozen / "staged.sha256",
                subject.common.sha_manifest(result, staged_files),
            )
            prepare = {
                "staged_seal_sha256":
                    subject.common.sha256_file(frozen / "staged.sha256"),
                "contract_sha256":
                    subject.common.sha256_file(frozen / "contract.json"),
            }
            subject.verify_frozen_staged_anchor(result, prepare)
            alias = frozen / "binary-alias"
            (frozen / "binary").rename(alias)
            (frozen / "binary").symlink_to(alias.name)
            with self.assertRaises(subject.CampaignError):
                subject.verify_frozen_staged_anchor(result, prepare)

            (frozen / "binary").unlink()
            alias.rename(frozen / "binary")
            subject.verify_frozen_staged_anchor(result, prepare)
            subject.common.atomic_json(frozen / "contract.json", {"version": 2})
            subject.common.atomic_write(frozen / "binary", b"replacement\n")
            subject.common.atomic_write(
                frozen / "staged.sha256",
                subject.common.sha_manifest(result, staged_files),
            )
            with self.assertRaises(subject.CampaignError):
                subject.verify_frozen_staged_anchor(result, prepare)
            q0 = {"passed": True, "comparisons": 1}
            subject.common.atomic_json(frozen / "q0_identity.json", q0)
            q0_sha = subject.common.sha256_file(frozen / "q0_identity.json")
            contract = {
                "source_commit": "1" * 40,
                "binary_sha256": "2" * 64,
                "script_sha256": "3" * 64,
                "allk_sha256": "4" * 64,
                "helper_sha256": "5" * 64,
                "protocol_sha256": "6" * 64,
                "route_context_sha256": "7" * 64,
                "q0_identity_sha256": q0_sha,
                "drand_wrapper_sha256": "8" * 64,
                "drand_package_sha256": "9" * 64,
                "drand_package_lock_sha256": "a" * 64,
                "drand_client_bundle_sha256": "b" * 64,
                "node": {"version": subject.DRAND_NODE_VERSION},
                "npm": {"version": subject.DRAND_NPM_VERSION},
                "drand_offline_selftest": {"passed": True},
            }
            bound_prepare = {**contract, "q0_identity": q0}
            subject.verify_prepare_contract_bindings(
                frozen, bound_prepare, contract)
            for field in contract:
                mutated = dict(bound_prepare)
                mutated[field] = None
                with self.assertRaises(subject.CampaignError, msg=field):
                    subject.verify_prepare_contract_bindings(
                        frozen, mutated, contract)
            mutated_q0 = dict(bound_prepare)
            mutated_q0["q0_identity"] = {"passed": False}
            with self.assertRaises(subject.CampaignError):
                subject.verify_prepare_contract_bindings(
                    frozen, mutated_q0, contract)

    def test_seal_attempt_never_replaces_a_dangling_final_name(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-seal-name-") as temporary:
            result = Path(temporary).resolve()
            spec = subject.HOLDOUT_PHASES["h1"]
            seal_record = {"schema": "test.seal.v1"}
            seal_bytes = subject.fixed_json_bytes(seal_record)
            seal_sha256 = hashlib.sha256(seal_bytes).hexdigest()
            seal_parent = result / "seals" / spec.seal_tree
            seal_parent.mkdir(parents=True)
            final = seal_parent / seal_sha256
            final.symlink_to("missing-seal-directory", target_is_directory=True)

            with self.assertRaisesRegex(
                    subject.CampaignError, "digest already exists"):
                subject.write_seal_attempt(
                    result, spec, b"manifest\n", seal_record, seal_bytes)
            self.assertTrue(final.is_symlink())

    def test_directory_publication_never_replaces_a_racing_target(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-dir-publish-") as temporary:
            parent = Path(temporary).resolve()
            source = parent / "staging"
            target = parent / "published"
            source.mkdir()
            target.mkdir()
            source_inode = source.stat().st_ino
            target_inode = target.stat().st_ino

            with self.assertRaisesRegex(
                    subject.CampaignError, "target already exists"):
                subject.rename_directory_noreplace(source, target)
            self.assertEqual(source.stat().st_ino, source_inode)
            self.assertEqual(target.stat().st_ino, target_inode)

            target.rmdir()
            subject.rename_directory_noreplace(source, target)
            self.assertFalse(source.exists())
            self.assertEqual(target.stat().st_ino, source_inode)

    def test_candidate_build_policy_is_exact_native_release(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-build-policy-") as temporary:
            repo = Path(temporary).resolve()
            cache = repo / "CMakeCache.txt"
            entries = {
                "CMAKE_BUILD_TYPE": ("STRING", "Release"),
                "BUILD_TESTS": ("BOOL", "ON"),
                "BUILD_CODEC_V2": ("BOOL", "ON"),
                "WIREHAIR_BUILD_BENCHMARKS": ("BOOL", "ON"),
                "MARCH_NATIVE": ("BOOL", "ON"),
                "WIREHAIR_STRICT_WARNINGS": ("BOOL", "ON"),
                "WIREHAIR_ENABLE_LIBFUZZER": ("BOOL", "OFF"),
                "WH_LTO": ("STRING", "OFF"),
                "WH_PGO_MODE": ("STRING", "OFF"),
                "CMAKE_C_FLAGS": ("STRING", ""),
                "CMAKE_CXX_FLAGS": ("STRING", ""),
                "CMAKE_C_FLAGS_RELEASE": ("STRING", "-O3 -DNDEBUG"),
                "CMAKE_CXX_FLAGS_RELEASE": ("STRING", "-O3 -DNDEBUG"),
                "CMAKE_EXE_LINKER_FLAGS": ("STRING", ""),
                "CMAKE_EXE_LINKER_FLAGS_RELEASE": ("STRING", ""),
                "CMAKE_SHARED_LINKER_FLAGS": ("STRING", ""),
                "CMAKE_SHARED_LINKER_FLAGS_RELEASE": ("STRING", ""),
                "CMAKE_MODULE_LINKER_FLAGS": ("STRING", ""),
                "CMAKE_MODULE_LINKER_FLAGS_RELEASE": ("STRING", ""),
                "CMAKE_C_COMPILER": ("FILEPATH", "/usr/bin/cc"),
                "CMAKE_CXX_COMPILER": ("FILEPATH", "/usr/bin/c++"),
                "CMAKE_GENERATOR": ("INTERNAL", "Unix Makefiles"),
                "CMAKE_HOME_DIRECTORY": ("INTERNAL", str(repo)),
            }

            def publish(values: dict[str, tuple[str, str]]) -> None:
                cache.write_text("".join(
                    f"{key}:{kind}={value}\n"
                    for key, (kind, value) in values.items()),
                    encoding="utf-8")

            publish(entries)
            policy = subject.validate_candidate_build_policy(cache, repo)
            self.assertTrue(policy["march_native"])
            self.assertEqual(policy["build_type"], "Release")
            self.assertTrue(policy["tests_enabled"])
            for key, value in (
                    ("MARCH_NATIVE", ("BOOL", "OFF")),
                    ("WH_PGO_MODE", ("STRING", "USE")),
                    ("CMAKE_CONFIGURATION_TYPES",
                     ("STRING", "Debug;Release")),
                    ("CMAKE_GENERATOR",
                     ("INTERNAL", "Ninja Multi-Config")),
                    ("CMAKE_CXX_FLAGS_RELEASE",
                     ("STRING", "-O3 -DNDEBUG -fsanitize=address"))):
                with self.subTest(key=key):
                    mutated = dict(entries)
                    mutated[key] = value
                    publish(mutated)
                    with self.assertRaises(subject.CampaignError):
                        subject.validate_candidate_build_policy(cache, repo)

    def test_prepared_result_root_rejects_relocation_and_aliases(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-result-root-binding-") as temporary:
            parent = Path(temporary).resolve()
            original = parent / "campaign"
            original.mkdir()
            prepare = {"result_dir": str(original)}
            self.assertEqual(
                subject.verify_prepared_result_root(original, prepare),
                original)

            relocated = parent / "campaign-copy"
            relocated.mkdir()
            with self.assertRaisesRegex(
                    subject.CampaignError, "differs from the public freeze"):
                subject.verify_prepared_result_root(relocated, prepare)

            alias = parent / "campaign-alias"
            alias.symlink_to(original, target_is_directory=True)
            with self.assertRaisesRegex(
                    subject.CampaignError, "differs from the public freeze"):
                subject.verify_prepared_result_root(alias, prepare)

            noncanonical = {"result_dir": str(original) + "/."}
            with self.assertRaisesRegex(
                    subject.CampaignError, "differs from the public freeze"):
                subject.verify_prepared_result_root(original, noncanonical)

    def test_phase_manifests_exclude_only_operational_locks(self) -> None:
        for phase, spec in subject.HOLDOUT_PHASES.items():
            with self.subTest(phase=phase), tempfile.TemporaryDirectory(
                    prefix=f"wh2-{phase}-manifest-locks-") as temporary:
                result = Path(temporary)
                for relative in spec.required_files:
                    path = result / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes((relative + "\n").encode("ascii"))
                before = subject.phase_manifest_bytes(result, spec)
                operational = (
                    result / "locks/development-control.lock",
                    result / "locks/h1-paired.lock",
                    result / "beacon/h1.lock",
                    result / "beacon/h2.lock",
                )
                for path in operational:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b"")
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                operational[0].write_bytes(b"operational metadata\n")
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                future_outputs = tuple(
                    relative for relative in spec.forbidden_before_seal
                    if relative == "timing" or relative.startswith("thermal/"))
                self.assertTrue(future_outputs)
                for relative in future_outputs:
                    path = result / relative
                    if relative == "timing":
                        path = path / "sealed-report.json"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b"future phase output\n")
                    self.assertTrue(any(
                        subject.is_relative_prefix(
                            path.relative_to(result).as_posix(), prefix)
                        for prefix in spec.excluded_prefixes))
                self.assertEqual(subject.phase_manifest_bytes(result, spec), before)
                (result / "unexpected.txt").write_bytes(b"not a lock\n")
                with self.assertRaises(subject.CampaignError):
                    subject.phase_manifest_bytes(result, spec)

    def test_development_preseal_requires_exact_terminal_manifests(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-development-readiness-") as temporary:
            result = Path(temporary)
            phase_files: set[Path] = set()
            for top in ("jobs", "ledgers", "completed", "thermal"):
                path = result / top / "artifact.bin"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes((top + "\n").encode("ascii"))
                phase_files.add(path)
            development = result / "development"
            development.mkdir()
            table_path = development / "frozen_table.json"
            table_path.write_bytes(b"{}\n")
            phase_path = development / "phase_complete.sha256"
            phase_path.write_bytes(subject.common.sha_manifest(
                result, phase_files))
            analysis_path = development / "analysis_complete.sha256"
            analysis_path.write_bytes(subject.common.sha_manifest(
                result, {phase_path, table_path}))

            route_context = "7" * 64
            contract = {
                "source_commit": "1" * 40,
                "protocol_sha256": "2" * 64,
                "route_context_sha256": route_context,
            }
            prepare = {"staged_seal_sha256": "3" * 64}
            table_record = {
                "source_commit": contract["source_commit"],
                "staged_seal_sha256": prepare["staged_seal_sha256"],
                "protocol_sha256": contract["protocol_sha256"],
                "route_context_sha256": route_context,
                "phase_seal_sha256": subject.common.sha256_file(phase_path),
            }
            fake_campaign = mock.Mock()
            fake_campaign.load_canonical_object.return_value = table_record

            def verify() -> None:
                with ExitStack() as stack:
                    stack.enter_context(mock.patch.object(
                        subject, "load_frozen_python_module",
                        side_effect=(fake_campaign, object())))
                    stack.enter_context(mock.patch.object(
                        subject, "load_frozen_cohort_bins",
                        return_value=((4096,), {0: (4096,)})))
                    stack.enter_context(mock.patch.object(
                        subject, "load_development_table_and_caches",
                        return_value=({4096: 1}, {0: object()}, "4" * 64)))
                    subject.verify_development_seal_readiness(
                        result, prepare, contract)

            verify()
            (development / "late-artifact.bin").write_bytes(b"late\n")
            with self.assertRaises(subject.CampaignError):
                verify()

    def test_h1_terminal_record_and_h2_required_files_are_exact(self) -> None:
        record = subject.h1_controller_complete_record(
            campaign_module,
            phase_completion_sha256="1" * 64,
            analysis_complete_sha256="2" * 64,
            root_file_sha256="3" * 64,
            parent_table_sha256="4" * 64,
            frozen_table_sha256="5" * 64,
            route_cache_index_sha256="6" * 64,
            route_context_sha256="7" * 64,
            input_routed_K=123,
            retained_routed_K=45,
        )
        campaign_module.verify_sealed_record(
            record,
            "wirehair.wh2.h12_preferred_attempt.h1_controller_complete.v1")
        self.assertEqual(record["cache_count"], subject.BIN_COUNT)
        required = set(subject.HOLDOUT_PHASES["h2"].required_files)
        self.assertTrue({
            "h1/controller_complete.json",
            "h1/controller_complete.json.sha256",
            "h1/route_cache/index.json",
            "h1/route_cache/index.json.sha256",
        }.issubset(required))
        self.assertIn(
            "verify_h1_seal_readiness", subject.run_h1.__code__.co_names)

    def test_holdout_wrapper_failure_replay_is_exact(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-wave-failure-replay-") as temporary:
            waves = Path(temporary)
            seal_sha256 = "1" * 64
            wrapper_failure = {
                "schema": "wirehair.wh2.drand_quicknet_wave.v1",
                "status": "PERMANENT_WRAPPER_FAILURE",
                "stdout_sha256": "2" * 64,
                "stderr_sha256": "3" * 64,
                "returncode": 1,
            }

            def outer(wave: int, completed_ms: int,
                      result: dict[str, object], returncode: int,
                      stderr_sha256: str) -> dict[str, object]:
                return {
                    "schema": "wirehair.wh2.holdout_wave_record.v1",
                    "phase": "h1",
                    "seal_record_sha256": seal_sha256,
                    "wave": wave,
                    "completed_ms": completed_ms,
                    "wrapper_returncode": returncode,
                    "wrapper_stderr_sha256": stderr_sha256,
                    "result": result,
                }

            path = waves / "wave-000001.json"
            valid = outer(1, 10, wrapper_failure, 1, "3" * 64)
            subject.atomic_fixed_json_once(path, valid)
            self.assertEqual(len(subject.load_existing_wave_records(
                waves, "h1", seal_sha256, 123)), 1)
            original = path.read_bytes()
            for mutation in ("extra", "stderr", "returncode", "negative"):
                path.unlink()
                changed = json.loads(original)
                if mutation == "extra":
                    changed["result"]["extra"] = True
                elif mutation == "stderr":
                    changed["result"]["stderr_sha256"] = "4" * 64
                elif mutation == "returncode":
                    changed["result"]["returncode"] = 2
                else:
                    changed["completed_ms"] = -1
                subject.atomic_fixed_json_once(path, changed)
                with self.subTest(mutation=mutation), self.assertRaises(
                        subject.CampaignError):
                    subject.load_existing_wave_records(
                        waves, "h1", seal_sha256, 123)

    def test_distinct_cpu_and_dimm_thermal_limits(self) -> None:
        def row(monotonic: int, cpu_c: int, dimm_c: int) -> bytes:
            values = [
                "2026-07-18T00:00:00Z", str(monotonic), "100", "4000",
                str(cpu_c), *([str(dimm_c)] * 8), "0", "1", "1", "1",
                "0", "0",
            ]
            self.assertEqual(len(values), len(subject.THERMAL_FIELDS))
            return (",".join(values) + "\n").encode("ascii")

        with tempfile.TemporaryDirectory(
                prefix="wh2-separate-thermal-limits-") as temporary:
            root = Path(temporary)
            thermal = root / "thermal.csv"
            header = (",".join(subject.THERMAL_FIELDS) + "\n").encode("ascii")
            baseline_row = row(1, 84, 89)
            interval_row = row(2, 85, 89)
            thermal.write_bytes(header + baseline_row + interval_row)
            metadata = thermal.stat()
            baseline = screen_module.parse_thermal_sample(
                baseline_row, "test baseline")
            mark = {
                "dev": metadata.st_dev,
                "ino": metadata.st_ino,
                "offset": len(header) + len(baseline_row),
                "header": header,
                "baseline_row": baseline_row,
                "baseline": baseline,
            }
            with mock.patch.object(
                    screen_module, "validate_thermal_current",
                    return_value=None):
                separate = subject.thermal_finish(
                    thermal, mark, root / "separate.csv", limit_c=85.0,
                    dimm_limit_c=90.0)
                legacy_equal = subject.thermal_finish(
                    thermal, mark, root / "legacy-equal.csv", limit_c=85.0)
            self.assertEqual(separate["thermal_high_samples"], 1)
            self.assertEqual(
                separate["thermal_high_max_consecutive_samples"], 1)
            self.assertEqual(legacy_equal["thermal_high_samples"], 2)
            self.assertEqual(separate["thermal_limit_c"], 85.0)
            with self.assertRaises(ValueError):
                subject.thermal_finish(
                    thermal, mark, root / "invalid.csv", limit_c=85.0,
                    dimm_limit_c=float("nan"))

    def test_timing_missing_sidecar_repair_requires_host_lock(self) -> None:
        timing = argparse.Namespace(
            PANEL_RESULT_NAME="panel.json",
            PANEL_RESULT_SIDECAR_NAME="panel.json.sha256")
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-lock-repair-") as temporary:
            panel = Path(temporary) / "panel"
            self.assertTrue(
                subject.timing_panel_needs_host_lock(timing, panel))
            panel.mkdir()
            self.assertFalse(
                subject.timing_panel_needs_host_lock(timing, panel))
            self.assertFalse(
                subject.timing_panel_has_complete_seal(timing, panel))
            (panel / "panel.json").write_bytes(b"{}\n")
            self.assertTrue(
                subject.timing_panel_needs_host_lock(timing, panel))
            (panel / "panel.json.sha256").symlink_to("missing")
            self.assertFalse(
                subject.timing_panel_needs_host_lock(timing, panel))
            self.assertFalse(
                subject.timing_panel_has_complete_seal(timing, panel))
            (panel / "panel.json.sha256").unlink()
            (panel / "panel.json.sha256").write_bytes(b"sealed\n")
            self.assertTrue(
                subject.timing_panel_has_complete_seal(timing, panel))
            journal = Path(temporary) / "host-journal.json"
            self.assertFalse(subject.timing_host_lock_required(
                timing, (panel,), journal))
            journal.write_bytes(b"stale host state\n")
            self.assertTrue(subject.timing_host_lock_required(
                timing, (panel,), journal))

    def test_sealed_timing_phase_is_structurally_load_only(self) -> None:
        timing = mock.Mock()
        spec = object()
        config = object()
        probe = object()
        panel = Path("/nonexistent/sealed-panel")
        timing.load_timing_panel_result.return_value = "sealed"
        self.assertEqual(
            subject.load_or_run_timing_panel(
                timing, True, spec, config, probe, panel),
            "sealed")
        timing.load_timing_panel_result.assert_called_once_with(
            panel, spec, config)
        timing.run_or_resume_timing_panel.assert_not_called()

        timing.reset_mock()
        timing.run_or_resume_timing_panel.return_value = "active"
        self.assertEqual(
            subject.load_or_run_timing_panel(
                timing, False, spec, config, probe, panel),
            "active")
        timing.run_or_resume_timing_panel.assert_called_once_with(
            spec, config, probe, panel)
        timing.load_timing_panel_result.assert_not_called()

    def test_timing_thermal_baseline_starts_strict_policy_at_phase(
            self) -> None:
        def row(monotonic_s: float, cpu_c: float) -> bytes:
            fields = [
                "2026-07-18T00:00:00Z", str(monotonic_s), "100", "4000",
                str(cpu_c), *(["50"] * 8), "0", "1", "1", "1", "0", "0",
            ]
            self.assertEqual(len(fields), len(subject.THERMAL_FIELDS))
            return (",".join(fields) + "\n").encode("ascii")

        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-thermal-baseline-") as temporary:
            result = Path(temporary).resolve()
            frozen = result / "frozen"
            frozen.mkdir()
            subject.common.atomic_json(frozen / "contract.json", {})
            timing_dir = result / "timing"
            timing_dir.mkdir()
            thermal = result / "thermal.csv"
            header = (",".join(subject.THERMAL_FIELDS) + "\n").encode("ascii")
            now = time.monotonic()
            prepare_row = row(now - 1.0, 80.0)
            thermal.write_bytes(header + prepare_row)
            prepare_mark = screen_module.thermal_start(
                thermal, require_zero_edac=False)
            contract = {
                "source_commit": "1" * 40,
                "thermal_baseline": {
                    "dev": prepare_mark["dev"], "ino": prepare_mark["ino"],
                    "offset": prepare_mark["offset"],
                    "edac_ce": prepare_mark["edac_ce"],
                    "edac_ue": prepare_mark["edac_ue"],
                    "monotonic_s": prepare_mark["monotonic_s"],
                    "max_temperature_c":
                        prepare_mark["max_temperature_c"],
                    "row_sha256": subject.sha256_bytes(prepare_row),
                },
            }
            policy = {
                "cpu_limit_c": 90.0, "dimm_limit_c": 90.0,
                "timing_cpu_limit_c": 85.0,
                "timing_dimm_limit_c": 90.0,
                "consecutive_samples": 3, "stale_seconds": 5.0,
            }
            with thermal.open("ab") as output:
                for offset in (0.9, 0.8, 0.7):
                    output.write(row(now - offset, 86.0))
            identity = campaign_module.validate_frozen_thermal_source(
                contract, thermal, policy)
            stale = timing_dir / ".thermal_baseline.json.99999999.partial"
            stale.write_bytes(b"interrupted baseline\n")
            baseline, path, digest, timing_identity = \
                subject.load_or_create_timing_thermal_baseline(
                    result, timing_dir, campaign_module, contract,
                    thermal, policy, identity, False)
            self.assertFalse(stale.exists())
            self.assertEqual(timing_identity, identity)
            self.assertGreater(
                baseline["offset"], contract["thermal_baseline"]["offset"])
            self.assertEqual(subject.verify_named_hash_sidecar(path), digest)

            path.with_suffix(".json.sha256").unlink()
            resumed = subject.load_or_create_timing_thermal_baseline(
                result, timing_dir, campaign_module, contract,
                thermal, policy, identity, False)
            self.assertEqual(resumed[0], baseline)
            self.assertEqual(resumed[2], digest)
            path.with_suffix(".json.sha256").unlink()
            with self.assertRaisesRegex(
                    subject.CampaignError, "lacks its thermal baseline"):
                subject.load_or_create_timing_thermal_baseline(
                    result, timing_dir, campaign_module, contract,
                    thermal, policy, identity, True)
            resumed = subject.load_or_create_timing_thermal_baseline(
                result, timing_dir, campaign_module, contract,
                thermal, policy, identity, False)

            with thermal.open("ab") as output:
                output.write(row(now - 0.6, 86.0))
                output.write(row(now - 0.5, 86.0))
            with self.assertRaisesRegex(
                    subject.CampaignError, "consecutive limit"):
                subject.validate_timing_thermal_histories(
                    campaign_module, contract, thermal, policy, baseline)
            # Once the timing phase is durably sealed, later recovery rows are
            # governed by the 90 C recovery policy, not retroactively by 85 C.
            completed = subject.load_or_create_timing_thermal_baseline(
                result, timing_dir, campaign_module, contract,
                thermal, policy, identity, True)
            self.assertEqual(completed[0], baseline)
            self.assertEqual(completed[3], identity)

    def test_timing_host_journal_precedes_first_mocked_mutation(self) -> None:
        boot_id = "11111111-2222-3333-4444-555555555555"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-host-order-") as temporary:
            root = Path(temporary).resolve()
            record = self.fake_timing_host_journal(root, boot_id)
            session = subject.TimingHostSession(
                {"system_tools": {}}, timing, 0)
            session.lock_path = root / "host.lock"
            session.journal_path = root / "host.json"
            events: list[str] = []
            real_write = session._write_host_journal

            def write_journal(value: dict[str, object]) -> None:
                real_write(value)
                self.assertTrue(session.journal_path.exists())
                events.append("durable-journal")

            def first_mutation() -> None:
                self.assertTrue(session.journal_path.exists())
                events.append("first-mutation")
                raise subject.CampaignError("simulated mutation failure")

            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}))
                stack.enter_context(mock.patch.object(
                    session, "_install_signal_handlers"))
                stack.enter_context(mock.patch.object(
                    session, "_quiesce_signal_handlers"))
                stack.enter_context(mock.patch.object(
                    session, "_restore_signal_handlers"))
                stack.enter_context(mock.patch.object(
                    session, "checked", return_value=b""))
                stack.enter_context(mock.patch.object(
                    session, "tool", side_effect=lambda name: f"/mock/{name}"))
                stack.enter_context(mock.patch.object(
                    session, "_snapshot_host_journal", return_value=record))
                stack.enter_context(mock.patch.object(
                    session, "_write_host_journal",
                    side_effect=write_journal))
                stack.enter_context(mock.patch.object(
                    session, "_validate_restored_host_journal"))
                stack.enter_context(mock.patch.object(
                    session, "_stop_fillers", side_effect=first_mutation))
                stack.enter_context(mock.patch.object(
                    session, "_restore_host_state", return_value=[]))
                with self.assertRaisesRegex(
                        subject.CampaignError, "simulated mutation failure"):
                    session.__enter__()
            self.assertEqual(
                events, ["durable-journal", "first-mutation"])
            self.assertFalse(session.journal_path.exists())
            self.assertIsNone(session.lock_fd)

    def test_timing_host_never_signals_a_reused_numeric_pid(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        session = subject.TimingHostSession({}, timing, 0)
        with mock.patch.object(
                session, "_process_start_ticks", return_value=200), \
             mock.patch.object(
                 subject.os, "pidfd_open", return_value=77, create=True), \
             mock.patch.object(
                 subject.signal, "pidfd_send_signal", create=True) as send, \
             mock.patch.object(subject.os, "close") as close:
            session._terminate_process_identities(((12345, 100),))
        send.assert_not_called()
        close.assert_called_once_with(77)

    def test_timing_host_journal_discharge_is_descriptor_bound(self) -> None:
        boot_id = "11111111-2222-3333-4444-555555555555"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-journal-binding-") as temporary:
            root = Path(temporary).resolve()
            record = self.fake_timing_host_journal(root, boot_id)
            session = subject.TimingHostSession({}, timing, 0)
            session.journal_path = root / "host.json"
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}):
                session._write_host_journal(record)
            original = session.journal_path.with_name("original.json")
            session.journal_path.rename(original)
            session.journal_path.write_bytes(session.journal_bytes or b"")
            os.chmod(session.journal_path, 0o600)
            with self.assertRaisesRegex(
                    subject.CampaignError, "changed while being read"):
                session._delete_host_journal()
            self.assertTrue(session.journal_path.exists())
            self.assertTrue(original.exists())
            session._clear_host_snapshot()

    def test_timing_host_recovers_linked_journal_commit_window(self) -> None:
        boot_id = "11111111-2222-3333-4444-555555555555"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-journal-linked-") as temporary:
            root = Path(temporary).resolve()
            record = self.fake_timing_host_journal(root, boot_id)
            writer = subject.TimingHostSession({}, timing, 0)
            writer.journal_path = root / "host.json"
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}):
                writer._write_host_journal(record)
            writer._clear_host_snapshot()
            marker = root / ".host.json.99999999.partial"
            os.link(writer.journal_path, marker)
            self.assertEqual(writer.journal_path.stat().st_nlink, 2)

            recovery = subject.TimingHostSession({}, timing, 0)
            recovery.journal_path = writer.journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}):
                self.assertEqual(recovery._load_host_journal(), record)
            self.assertFalse(marker.exists())
            self.assertEqual(writer.journal_path.stat().st_nlink, 1)
            recovery._clear_host_snapshot()

    def test_timing_host_rejects_nonprivate_global_lock(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-lock-mode-") as temporary:
            root = Path(temporary).resolve()
            lock = root / "host.lock"
            lock.write_bytes(b"")
            lock.chmod(0o644)
            session = subject.TimingHostSession({}, timing, 0)
            session.lock_path = lock
            session.journal_path = root / "host.json"
            with self.assertRaisesRegex(
                    subject.CampaignError, "unique regular file"):
                session.__enter__()
            self.assertIsNone(session.lock_fd)

    def test_timing_host_lock_acquisition_error_closes_descriptor(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-lock-error-") as temporary:
            root = Path(temporary).resolve()
            session = subject.TimingHostSession({}, timing, 0)
            session.lock_path = root / "host.lock"
            session.journal_path = root / "host.json"
            real_open = os.open
            descriptors: list[int] = []

            def opened(*args: object, **kwargs: object) -> int:
                descriptor = real_open(*args, **kwargs)
                descriptors.append(descriptor)
                return descriptor

            with mock.patch.object(subject.os, "open", side_effect=opened), \
                 mock.patch.object(
                     subject.os, "fstat", side_effect=OSError("fstat fault")), \
                 self.assertRaisesRegex(OSError, "fstat fault"):
                session.__enter__()
            self.assertIsNone(session.lock_fd)
            self.assertEqual(len(descriptors), 1)
            with self.assertRaises(OSError):
                os.fstat(descriptors[0])

    def test_timing_host_malformed_journal_closes_bound_descriptor(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        session = subject.TimingHostSession({}, timing, 0)
        descriptor = os.open("/dev/null", os.O_RDONLY)
        with mock.patch.object(
                session, "_open_bound_host_journal",
                return_value=(descriptor, b"[")), \
             mock.patch.object(
                 subject.json, "loads", side_effect=RecursionError("deep")), \
             self.assertRaisesRegex(RecursionError, "deep"):
            session._load_host_journal()
        self.assertIsNone(session.journal_fd)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_timing_host_clear_error_still_releases_global_lock(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-close-error-") as temporary:
            lock = Path(temporary) / "host.lock"
            descriptor = os.open(lock, os.O_CREAT | os.O_RDWR, 0o600)
            fcntl.flock(descriptor, fcntl.LOCK_EX)
            session = subject.TimingHostSession({}, timing, 0)
            session.lock_fd = descriptor
            with mock.patch.object(
                    session, "_quiesce_signal_handlers"), \
                 mock.patch.object(session, "_restore_signal_handlers"), \
                 mock.patch.object(
                     session, "_clear_host_snapshot",
                     side_effect=OSError("journal close fault")), \
                 self.assertRaisesRegex(
                     subject.CampaignError, "journal close fault"):
                session.close()
            self.assertIsNone(session.lock_fd)
            with self.assertRaises(OSError):
                os.fstat(descriptor)

    def test_timing_host_recovers_durable_crash_journal_exactly(self) -> None:
        boot_id = "11111111-2222-3333-4444-555555555555"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-host-journal-") as temporary:
            root = Path(temporary).resolve()
            journal_path = root / "host.json"
            record = self.fake_timing_host_journal(root, boot_id)
            writer = subject.TimingHostSession({}, timing, 0)
            writer.journal_path = journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}):
                writer._write_host_journal(record)
            persisted = journal_path.read_bytes()
            self.assertIn(b'"recovery_tools"', persisted)
            self.assertIn(b'"single_cpu_affinity": true', persisted)

            # Simulate the state left by a hard kill after all mutations: the
            # sibling is offline, cgroups are narrowed, and fillers are gone.
            sibling_state = {64: "0"}
            allowed = {unit: "1" for unit in subject.TIMING_HOST_UNITS}
            effective = {unit: "1" for unit in subject.TIMING_HOST_UNITS}
            desired_effective = {
                item["unit"]: item["effective"]
                for item in record["slices"]
            }
            filler_inventory: list[
                tuple[
                    tuple[int, ...], tuple[int, ...],
                    tuple[tuple[int, int], ...],
                ]
            ] = [((), (), ())]
            events: list[tuple[object, ...]] = []

            def write_online(cpu: int, value: str) -> None:
                self.assertEqual(
                    recovery.tools, record["recovery_tools"])
                events.append(("online", cpu, value))
                sibling_state[cpu] = value

            def set_allowed(unit: str, value: str) -> None:
                self.assertEqual(
                    recovery.tools, record["recovery_tools"])
                events.append(("allowed", unit, value))
                allowed[unit] = value.replace(",", " ")
                effective[unit] = (
                    value.replace(",", " ")
                    if value else desired_effective[unit])

            def restart_fillers() -> None:
                self.assertEqual(
                    recovery.tools, record["recovery_tools"])
                events.append(("fillers",))
                filler_inventory[0] = (
                    (101, 102), (0, 1), ((101, 1), (102, 1)))

            recovery = subject.TimingHostSession({}, timing, 0)
            recovery.journal_path = journal_path
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}))
                stack.enter_context(mock.patch.object(
                    recovery, "_read_sibling_state",
                    side_effect=lambda cpu: sibling_state[cpu]))
                stack.enter_context(mock.patch.object(
                    recovery, "_write_online", side_effect=write_online))
                stack.enter_context(mock.patch.object(
                    recovery, "_set_allowed", side_effect=set_allowed))
                stack.enter_context(mock.patch.object(
                    recovery, "_systemctl_allowed",
                    side_effect=lambda unit: allowed[unit]))
                stack.enter_context(mock.patch.object(
                    recovery, "_systemctl_effective",
                    side_effect=lambda unit: effective[unit]))
                stack.enter_context(mock.patch.object(
                    recovery, "_restart_fillers",
                    side_effect=restart_fillers))
                stack.enter_context(mock.patch.object(
                    recovery, "_filler_inventory",
                    side_effect=lambda: filler_inventory[0]))
                stack.enter_context(mock.patch.object(
                    recovery, "_matching_filler_inventory",
                    side_effect=lambda: filler_inventory[0]))
                recovery._recover_existing_host_journal()

            self.assertFalse(journal_path.exists())
            self.assertEqual(sibling_state, {64: "1"})
            self.assertEqual(allowed, {
                unit: "" for unit in subject.TIMING_HOST_UNITS})
            self.assertEqual(effective, desired_effective)
            self.assertEqual(
                filler_inventory[0],
                ((101, 102), (0, 1), ((101, 1), (102, 1))))
            self.assertEqual(events[0], ("online", 64, "1"))
            self.assertEqual(events[-1], ("fillers",))
            self.assertIsNone(recovery.journal_record)
            self.assertIsNone(recovery.journal_bytes)

    def test_timing_host_retains_cross_boot_journal_without_mutation(self) -> None:
        old_boot = "11111111-2222-3333-4444-555555555555"
        new_boot = "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-host-old-boot-") as temporary:
            root = Path(temporary).resolve()
            journal_path = root / "host.json"
            record = self.fake_timing_host_journal(root, old_boot)
            writer = subject.TimingHostSession({}, timing, 0)
            writer.journal_path = journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": old_boot}):
                writer._write_host_journal(record)
            persisted = journal_path.read_bytes()

            recovery = subject.TimingHostSession({}, timing, 0)
            recovery.journal_path = journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": new_boot}), \
                 mock.patch.object(recovery, "_write_online") as write_online, \
                 mock.patch.object(recovery, "_set_allowed") as set_allowed, \
                 mock.patch.object(
                     recovery, "_restart_fillers") as restart_fillers, \
                 self.assertRaisesRegex(
                     subject.CampaignError, "different boot"):
                recovery._recover_existing_host_journal()
            write_online.assert_not_called()
            set_allowed.assert_not_called()
            restart_fillers.assert_not_called()
            self.assertEqual(journal_path.read_bytes(), persisted)

    def test_timing_host_journal_survives_failed_exact_validation(self) -> None:
        boot_id = "11111111-2222-3333-4444-555555555555"
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        with tempfile.TemporaryDirectory(
                prefix="wh2-timing-host-retained-") as temporary:
            root = Path(temporary).resolve()
            journal_path = root / "host.json"
            record = self.fake_timing_host_journal(root, boot_id)
            writer = subject.TimingHostSession({}, timing, 0)
            writer.journal_path = journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}):
                writer._write_host_journal(record)
            persisted = journal_path.read_bytes()

            recovery = subject.TimingHostSession({}, timing, 0)
            recovery.journal_path = journal_path
            with mock.patch.object(
                    subject, "host_runtime_identity",
                    return_value={"boot_id": boot_id}), \
                 mock.patch.object(
                     recovery, "_restore_host_state",
                     return_value=[subject.CampaignError(
                         "exact validation failed")]), \
                 mock.patch.object(
                     recovery, "_delete_host_journal") as delete, \
                 self.assertRaisesRegex(
                     subject.CampaignError, "exact validation failed"):
                recovery._recover_existing_host_journal()
            delete.assert_not_called()
            self.assertEqual(journal_path.read_bytes(), persisted)

    def test_timing_host_defers_signals_until_mocked_restore_finishes(self) -> None:
        timing = mock.Mock()
        timing.inspect_linux_isolation.return_value = argparse.Namespace(
            sibling_cpus=(64,), numa_node=0)
        session = subject.TimingHostSession({}, timing, 0)
        previous = {
            signum: signal.getsignal(signum)
            for signum in (signal.SIGINT, signal.SIGTERM)
        }
        session._install_signal_handlers()
        try:
            with self.assertRaises(KeyboardInterrupt):
                session._handle_timing_signal(signal.SIGINT, None)
            with self.assertRaises(subject.CampaignError):
                session._handle_timing_signal(signal.SIGTERM, None)
            session.sibling_states = {64: "1"}
            session.slice_states = {"system.slice": ("", "0 64")}
            session.active = True

            def restart_fillers() -> None:
                handler = signal.getsignal(signal.SIGTERM)
                self.assertIs(
                    getattr(handler, "__self__", None), session)
                self.assertIs(
                    getattr(handler, "__func__", None),
                    subject.TimingHostSession._defer_timing_signal)
                handler(signal.SIGTERM, None)

            with ExitStack() as stack:
                read_text = stack.enter_context(mock.patch.object(
                    subject.Path, "read_text", side_effect=("0", "1")))
                write_online = stack.enter_context(mock.patch.object(
                    session, "_write_online"))
                set_allowed = stack.enter_context(mock.patch.object(
                    session, "_set_allowed"))
                effective = stack.enter_context(mock.patch.object(
                    session, "_systemctl_effective",
                    side_effect=("0 64", "0 64")))
                allowed = stack.enter_context(mock.patch.object(
                    session, "_systemctl_allowed", return_value=""))
                stack.enter_context(mock.patch.object(
                    session, "_restart_fillers", side_effect=restart_fillers))
                with self.assertRaisesRegex(
                        subject.CampaignError,
                        "interrupted during restoration by SIGTERM"):
                    session.close()
            self.assertEqual(read_text.call_count, 2)
            write_online.assert_called_once_with(64, "1")
            self.assertEqual(set_allowed.call_args_list, [
                mock.call("system.slice", "0,64"),
                mock.call("system.slice", ""),
            ])
            self.assertEqual(effective.call_count, 2)
            allowed.assert_called_once_with("system.slice")
            self.assertFalse(session.active)
            self.assertFalse(session.signals_installed)
            self.assertEqual(session.sibling_states, {})
            self.assertEqual(session.slice_states, {})
            for signum, handler in previous.items():
                self.assertEqual(signal.getsignal(signum), handler)
        finally:
            if session.signals_installed:
                session.pending_cleanup_signals.clear()
                session._restore_signal_handlers()

    def test_interrupted_turbostat_communication_is_reaped(self) -> None:
        probe = subject.LinuxTimingEvidenceProbe(
            mock.Mock(), mock.Mock(), Path("/thermal"), Path("/evidence"))
        for communicate_effect in (
                [KeyboardInterrupt()],
                [subprocess.TimeoutExpired("turbostat", 10),
                 KeyboardInterrupt()]):
            with self.subTest(calls=len(communicate_effect)):
                process = mock.Mock()
                process.pid = 12345
                process.poll.return_value = None
                process.communicate.side_effect = communicate_effect
                token = subject.TimingEvidenceToken(
                    {}, process, {}, Path("/thermal-out"),
                    Path("/performance-out"), ("turbostat",))
                with ExitStack() as stack:
                    stack.enter_context(mock.patch.object(subject.os, "killpg"))
                    reap = stack.enter_context(mock.patch.object(
                        subject.common, "stop_and_reap_process_group"))
                    with self.assertRaises(KeyboardInterrupt):
                        probe.finish(token)
                reap.assert_called_once_with(process)

    def test_successful_turbostat_leader_cannot_leave_a_child(self) -> None:
        process = subprocess.Popen(
            (
                "/bin/sh", "-c",
                "(sleep 20) >/dev/null 2>&1 & exit 0",
            ),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True)
        process.wait(timeout=2.0)
        probe = subject.LinuxTimingEvidenceProbe(
            mock.Mock(), mock.Mock(), Path("/thermal"), Path("/evidence"))
        token = subject.TimingEvidenceToken(
            {}, process, {}, Path("/thermal-out"),
            Path("/performance-out"), ("turbostat",))
        with self.assertRaisesRegex(
                subject.CampaignError, "unexpected child process"):
            probe.finish(token)
        self.assertFalse(subject.common.process_group_exists(process))

    def test_interrupted_turbostat_startup_is_reaped(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="wh2-turbostat-startup-") as temporary:
            host = mock.Mock()
            host.isolation_is_live.return_value = True
            host.core = 0
            host.isolation.sibling_cpus = ()
            host.tool.side_effect = lambda name: f"/frozen/{name}"
            frozen_baseline = (1, 2, 5, 6)
            probe = subject.LinuxTimingEvidenceProbe(
                mock.Mock(), host, Path("/thermal"), Path(temporary),
                frozen_baseline)
            process = mock.Mock()
            process.pid = 12345
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "thermal_start", return_value={
                        "dev": 1, "ino": 2, "edac_ce": 5, "edac_ue": 6,
                    }))
                stack.enter_context(mock.patch.object(
                    subject.subprocess, "Popen", return_value=process))
                stack.enter_context(mock.patch.object(
                    subject.time, "sleep", side_effect=KeyboardInterrupt))
                reap = stack.enter_context(mock.patch.object(
                    subject.common, "stop_and_reap_process_group"))
                with self.assertRaises(KeyboardInterrupt):
                    probe.start(mock.Mock(), 1, None)
            reap.assert_called_once_with(process)

    def test_state_root_restart_and_manifest_mutation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="wh2-state-test-") as temporary:
            parent = Path(temporary)
            result = parent / "campaign"
            result.mkdir()
            self.required_h1_tree(result)
            repo, commit = self.fake_repo(parent)
            contract = self.fake_contract(repo, commit)
            prepare = {
                "prepared_utc_ns":
                    (subject.now_utc_ms() - 1000) * 1_000_000,
            }
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(
                    subject, "verify_frozen_controller_runtime",
                    return_value=(prepare, contract),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "github_remote_url",
                    return_value="git@github.com:fake/wirehair.git",
                ))
                stack.enter_context(mock.patch.object(
                    subject, "newest_github_agent_environment",
                    return_value=os.environ.copy(),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "obtain_verified_latest_bound",
                    return_value=self.latest_bound(),
                ))
                stack.enter_context(mock.patch.object(
                    subject, "verify_phase_ready_for_seal",
                    return_value=None,
                ))
                stack.enter_context(redirect_stdout(io.StringIO()))
                dangling_state = result / "beacon/h1/state.json"
                dangling_state.parent.mkdir(parents=True, exist_ok=True)
                dangling_state.symlink_to("missing-state.json")
                with self.assertRaises(subject.CampaignError):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertTrue(dangling_state.is_symlink())
                dangling_state.unlink()
                self.assertEqual(subject.seal_holdout(argparse.Namespace(
                    result_dir=result, phase="h1")), 0)
                state = subject.load_fixed_json(
                    result / "beacon/h1/state.json", "test state")
                attempt, _manifest, seal, _seal_sha = subject.load_seal_attempt(
                    result, subject.HOLDOUT_PHASES["h1"], state)
                publication_path = attempt / "publication.json"
                publication_bytes = publication_path.read_bytes()
                stale_publication = attempt / \
                    ".publication.json.99999999.partial"
                stale_publication.write_bytes(b"interrupted receipt\n")
                interrupted_state = dict(state)
                interrupted_state["status"] = "SEALED"
                for key in (
                        "tag_name", "tag_object", "publication_confirmed_ms"):
                    interrupted_state.pop(key)
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)
                with mock.patch.object(
                        subject, "publish_seal_tag",
                        side_effect=AssertionError("must not republish")), \
                     mock.patch.object(
                         subject, "now_utc_ms",
                         return_value=seal["round_time_ms"] - 1):
                    self.assertEqual(subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1")), 0)
                state = subject.load_fixed_json(
                    result / "beacon/h1/state.json", "recovered state")
                self.assertEqual(state["status"], "WAITING_UNTIL_T")
                self.assertEqual(publication_path.read_bytes(), publication_bytes)
                self.assertFalse(stale_publication.exists())
                waiting_state = dict(state)

                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)
                with mock.patch.object(
                        subject, "now_utc_ms",
                        return_value=seal["round_time_ms"] + 1), \
                     self.assertRaisesRegex(
                         subject.CampaignError,
                         "durability could not be reconfirmed"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertEqual(
                    subject.load_fixed_json(
                        result / "beacon/h1/state.json", "unproven recovery")
                    ["status"],
                    "ABANDONED_PUBLICATION_UNPROVEN_LATE")
                self.assertEqual(publication_path.read_bytes(), publication_bytes)
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)

                confirmed_late_start = subject.load_fixed_json(
                    publication_path, "test confirmed publication")
                confirmed_late_start["receipt_write_started_ms"] = \
                    seal["round_time_ms"]
                subject.atomic_state_json(
                    publication_path, confirmed_late_start)
                with self.assertRaisesRegex(
                        subject.CampaignError, "started at or after T"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertEqual(
                    subject.load_fixed_json(
                        result / "beacon/h1/state.json", "late-start recovery")
                    ["status"], "ABANDONED_PUBLICATION_LATE")
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)
                subject.common.atomic_write(publication_path, publication_bytes)

                impossible_late = subject.load_fixed_json(
                    publication_path, "test publication")
                impossible_late["status"] = "ABANDONED_LATE"
                subject.atomic_state_json(publication_path, impossible_late)
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)
                impossible_bytes = publication_path.read_bytes()
                with mock.patch.object(
                        subject, "publish_seal_tag",
                        side_effect=AssertionError("must not republish")), \
                     self.assertRaisesRegex(
                         subject.CampaignError, "publication receipt mismatch"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertEqual(
                    subject.load_fixed_json(
                        result / "beacon/h1/state.json", "failed recovery")
                    ["status"], "SEALED")
                self.assertEqual(publication_path.read_bytes(), impossible_bytes)

                null_late = dict(impossible_late)
                null_late["confirmed_ms"] = seal["round_time_ms"]
                null_late["receipt_write_started_ms"] = \
                    seal["round_time_ms"]
                null_late["receipt_persisted_ms"] = None
                subject.atomic_state_json(publication_path, null_late)
                with self.assertRaisesRegex(
                        subject.CampaignError, "publication receipt mismatch"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertEqual(
                    subject.load_fixed_json(
                        result / "beacon/h1/state.json", "null recovery")
                    ["status"], "SEALED")

                first_write_late = dict(impossible_late)
                first_write_late["confirmed_ms"] = seal["round_time_ms"]
                first_write_late["receipt_write_started_ms"] = \
                    seal["round_time_ms"]
                subject.atomic_state_json(publication_path, first_write_late)
                with mock.patch.object(
                        subject, "publish_seal_tag",
                        side_effect=AssertionError("must not republish")), \
                     self.assertRaisesRegex(
                         subject.CampaignError, "completed at or after T"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                self.assertEqual(
                    subject.load_fixed_json(
                        result / "beacon/h1/state.json", "first-write late")
                    ["status"], "ABANDONED_PUBLICATION_LATE")
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", interrupted_state)

                abandoned_late = dict(impossible_late)
                abandoned_late["receipt_persisted_ms"] = \
                    seal["round_time_ms"]
                subject.atomic_state_json(publication_path, abandoned_late)
                abandoned_bytes = publication_path.read_bytes()
                with mock.patch.object(
                        subject, "publish_seal_tag",
                        side_effect=AssertionError("must not republish")), \
                     self.assertRaisesRegex(
                         subject.CampaignError, "completed at or after T"):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                abandoned_state = subject.load_fixed_json(
                    result / "beacon/h1/state.json", "abandoned recovery")
                self.assertEqual(
                    abandoned_state["status"],
                    "ABANDONED_PUBLICATION_LATE")
                self.assertEqual(publication_path.read_bytes(), abandoned_bytes)

                corrupted_late = dict(abandoned_late)
                corrupted_late["tag_name"] = "wrong-tag"
                subject.atomic_state_json(publication_path, corrupted_late)
                with mock.patch.object(
                        subject, "obtain_verified_latest_bound") as latest, \
                     self.assertRaises(subject.CampaignError):
                    subject.seal_holdout(argparse.Namespace(
                        result_dir=result, phase="h1"))
                latest.assert_not_called()

                subject.common.atomic_write(publication_path, publication_bytes)
                subject.atomic_state_json(
                    result / "beacon/h1/state.json", waiting_state)
                state = waiting_state
                beacon = dict(subject.DRAND_KNOWN_ROUND)
                beacon["round"] = seal["round"]
                origins = [f"https://{host}" for host in subject.DRAND_RELAYS]
                canonical_sha256 = hashlib.sha256(
                    subject.canonical_beacon_bytes(beacon)).hexdigest()
                wave = {
                    "schema": "wirehair.wh2.drand_quicknet_wave.v1",
                    "status": "QUORUM",
                    "chain_hash": subject.DRAND_CHAIN_HASH,
                    "consensus_origins": origins[:3],
                    "observations": [{
                        "origin": origin,
                        "ok": True,
                        "canonical_sha256": canonical_sha256,
                        "beacon": beacon,
                        "fetch_started_ms": 1 + index * 2,
                        "fetch_completed_ms": 2 + index * 2,
                        "raw_response_sha256": f"{index + 1:x}" * 64,
                    } for index, origin in enumerate(origins[:3])] + [{
                        "origin": origins[3],
                        "ok": False,
                        "error": "test timeout",
                    }],
                    "beacon": beacon,
                }
                real_run = subject.subprocess.run
                drand_calls: list[object] = []
                drand_wave = [wave]
                drand_returncode = [0]
                drand_stderr = [b"process fault\n"]
                waves_parent_synced = False

                real_fsync_directory = subject.fsync_directory

                def fsync_directory(path: Path) -> None:
                    nonlocal waves_parent_synced
                    real_fsync_directory(path)
                    waves = result / "beacon/h1/waves"
                    if path == waves.parent and waves.is_dir():
                        waves_parent_synced = True

                def run(command: object, *args: object, **kwargs: object):
                    if (isinstance(command, tuple) and len(command) >= 4 and
                            command[-1] == str(seal["round"]) and
                            "wh2_drand_verify.cjs" in str(command[1])):
                        self.assertTrue(waves_parent_synced)
                        drand_calls.append(command)
                        return subprocess.CompletedProcess(
                            command, drand_returncode[0],
                            json.dumps(drand_wave[0]).encode("ascii"),
                            drand_stderr[0])
                    return real_run(command, *args, **kwargs)

                def offline(
                    _contract: object, _result: object, value: dict[str, object],
                ) -> dict[str, object]:
                    return {
                        "beacon": value,
                        "canonical_sha256": hashlib.sha256(
                            subject.canonical_beacon_bytes(value)).hexdigest(),
                    }

                acquire = argparse.Namespace(
                    result_dir=result, phase="h1", no_wait=False,
                    max_wait_seconds=1.0,
                )
                with ExitStack() as acquire_stack:
                    acquire_stack.enter_context(mock.patch.object(
                        subject.subprocess, "run", side_effect=run))
                    acquire_stack.enter_context(mock.patch.object(
                        subject, "now_utc_ms",
                        return_value=seal["round_time_ms"] + 1,
                    ))
                    acquire_stack.enter_context(mock.patch.object(
                        subject, "offline_verify_beacon", side_effect=offline),
                    )
                    acquire_stack.enter_context(mock.patch.object(
                        subject, "fsync_directory",
                        side_effect=fsync_directory,
                    ))
                    with self.assertRaisesRegex(
                            subject.CampaignError, "process status changed"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 1)
                    failed_wave = \
                        result / "beacon/h1/waves/wave-000001.json"
                    failed_wave.unlink()
                    subject.fsync_directory(failed_wave.parent)
                    subject.atomic_state_json(
                        result / "beacon/h1/state.json", waiting_state)
                    drand_stderr[0] = b""
                    drand_returncode[0] = 2
                    drand_wave[0] = {
                        "status": "TEMPORARY_NO_QUORUM",
                    }
                    with self.assertRaisesRegex(
                            subject.CampaignError, "invalid coverage"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 2)
                    with self.assertRaisesRegex(
                            subject.CampaignError, "invalid coverage"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 2)
                    failed_wave.unlink()
                    subject.fsync_directory(failed_wave.parent)
                    subject.atomic_state_json(
                        result / "beacon/h1/state.json", waiting_state)
                    drand_returncode[0] = 0
                    drand_wave[0] = wave
                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    root = subject.load_fixed_json(
                        result / "holdout/h1_roots.json", "test root")
                    self.assertEqual(root["status"], "ROOTED")
                    self.assertEqual(len(root["roots"]), subject.H1_ROOT_COUNT)
                    self.assertEqual(len(drand_calls), 3)
                    root_path = result / "holdout/h1_roots.json"
                    state_path = result / "beacon/h1/state.json"
                    wave_path = result / "beacon/h1/waves/wave-000001.json"
                    root_bytes = root_path.read_bytes()
                    state_bytes = state_path.read_bytes()
                    wave_bytes = wave_path.read_bytes()

                    root_marker = root_path.parent / \
                        f".{root_path.name}.99999999.partial"
                    os.link(root_path, root_marker)
                    self.assertEqual(root_path.stat().st_nlink, 2)
                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    self.assertFalse(root_marker.exists())
                    self.assertEqual(root_path.stat().st_nlink, 1)
                    self.assertEqual(root_path.read_bytes(), root_bytes)

                    root_path.unlink()
                    subject.atomic_state_json(state_path, waiting_state)
                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    self.assertEqual(len(drand_calls), 3)
                    self.assertEqual(root_path.read_bytes(), root_bytes)

                    bad_stderr = subject.load_fixed_json(
                        wave_path, "test persisted wave")
                    bad_stderr["wrapper_stderr_sha256"] = "a" * 64
                    subject.atomic_state_json(wave_path, bad_stderr)
                    root_path.unlink()
                    subject.atomic_state_json(state_path, waiting_state)
                    with self.assertRaisesRegex(
                            subject.CampaignError, "process status changed"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    subject.common.atomic_write(wave_path, wave_bytes)
                    subject.common.atomic_write(root_path, root_bytes)

                    conflicting_beacon = dict(beacon)
                    conflicting_beacon["signature"] = (
                        ("0" if beacon["signature"][0] != "0" else "1") +
                        beacon["signature"][1:])
                    conflicting_beacon["randomness"] = hashlib.sha256(
                        bytes.fromhex(
                            conflicting_beacon["signature"])).hexdigest()
                    permanent = subject.load_fixed_json(
                        wave_path, "test persisted wave")
                    permanent["wrapper_returncode"] = 3
                    permanent["result"] = {
                        "schema": "wirehair.wh2.drand_quicknet_wave.v1",
                        "status": "PERMANENT_VERIFIED_DISAGREEMENT",
                        "chain_hash": subject.DRAND_CHAIN_HASH,
                        "observations": [
                            dict(wave["observations"][0]),
                            {
                                **dict(wave["observations"][1]),
                                "beacon": conflicting_beacon,
                                "canonical_sha256": hashlib.sha256(
                                    subject.canonical_beacon_bytes(
                                        conflicting_beacon)).hexdigest(),
                            },
                            {
                                "origin": origins[2],
                                "ok": False,
                                "error": "test disagreement",
                            },
                            dict(wave["observations"][3]),
                        ],
                    }
                    subject.atomic_state_json(wave_path, permanent)
                    root_path.unlink()
                    subject.atomic_state_json(state_path, waiting_state)
                    with self.assertRaisesRegex(
                            subject.CampaignError, "permanently blocked"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    self.assertEqual(subject.load_fixed_json(
                        state_path, "blocked replay")["status"],
                        "BLOCKED_PERMANENT")

                    subject.common.atomic_write(wave_path, wave_bytes)
                    alternate = parent / "alternate-wave-ledger"
                    alternate.mkdir()
                    waves_dir = wave_path.parent
                    hidden = waves_dir.with_name("waves-hidden")
                    waves_dir.rename(hidden)
                    waves_dir.symlink_to(alternate, target_is_directory=True)
                    root_path.unlink(missing_ok=True)
                    subject.atomic_state_json(state_path, waiting_state)
                    with mock.patch.object(
                            subject, "verify_manifest_bytes",
                            return_value=None), \
                         self.assertRaisesRegex(
                             subject.CampaignError, "plain directory"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    waves_dir.unlink()
                    hidden.rename(waves_dir)

                    subject.common.atomic_write(root_path, root_bytes)
                    subject.common.atomic_write(state_path, state_bytes)
                    waves_dir.rename(hidden)
                    waves_dir.symlink_to(hidden, target_is_directory=True)
                    with mock.patch.object(
                            subject, "verify_manifest_bytes",
                            return_value=None), \
                         self.assertRaisesRegex(
                             subject.CampaignError, "plain directory"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    waves_dir.unlink()
                    hidden.rename(waves_dir)

                    quorum_two = subject.load_fixed_json(
                        wave_path, "test quorum wave")
                    quorum_two["wave"] = 2
                    wave_two_path = waves_dir / "wave-000002.json"
                    subject.atomic_state_json(wave_two_path, quorum_two)
                    subject.atomic_state_json(wave_path, permanent)
                    rooted_after_terminal = dict(root)
                    rooted_after_terminal["wave"] = 2
                    rooted_after_terminal["wave_record_sha256"] = \
                        subject.common.sha256_file(wave_two_path)
                    subject.atomic_state_json(
                        root_path, rooted_after_terminal)
                    with self.assertRaisesRegex(
                            subject.CampaignError,
                            "continued after a terminal result"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    wave_two_path.unlink()
                    subject.fsync_directory(waves_dir)
                    subject.common.atomic_write(wave_path, wave_bytes)
                    subject.common.atomic_write(root_path, root_bytes)
                    subject.common.atomic_write(state_path, state_bytes)

                    external_root = parent / "external-root.json"
                    external_root.write_bytes(root_bytes)
                    root_path.unlink()
                    root_path.symlink_to(external_root)
                    subject.atomic_state_json(state_path, waiting_state)
                    with mock.patch.object(
                            subject, "verify_manifest_bytes",
                            return_value=None), \
                         self.assertRaisesRegex(
                             subject.CampaignError, "unique regular file"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    root_path.unlink()
                    os.link(external_root, root_path)
                    with mock.patch.object(
                            subject, "verify_manifest_bytes",
                            return_value=None), \
                         self.assertRaisesRegex(
                             subject.CampaignError, "unique regular file"):
                        subject.acquire_holdout(acquire)
                    self.assertEqual(len(drand_calls), 3)
                    root_path.unlink()
                    subject.common.atomic_write(root_path, root_bytes)
                    subject.common.atomic_write(state_path, state_bytes)

                    self.assertEqual(subject.acquire_holdout(acquire), 0)
                    self.assertEqual(len(list(
                        (result / "beacon/h1/waves").glob("wave-*.json"))), 1)
                    for field, value in (
                            ("manifest_sha256", "f" * 64),
                            ("round", root["round"] + 1),
                            ("round_time_ms", root["round_time_ms"] + 1),
                            ("wave", True),
                            ("rooted_ms", "invalid")):
                        mutated_root = dict(root)
                        mutated_root[field] = value
                        subject.atomic_state_json(root_path, mutated_root)
                        with self.assertRaises(subject.CampaignError):
                            subject.acquire_holdout(acquire)
                        subject.common.atomic_write(root_path, root_bytes)
                    mutated_root = dict(root)
                    mutated_root["unrecognized"] = "field"
                    subject.atomic_state_json(root_path, mutated_root)
                    with self.assertRaises(subject.CampaignError):
                        subject.acquire_holdout(acquire)
                    subject.common.atomic_write(root_path, root_bytes)
                    rooted_state = subject.load_fixed_json(
                        state_path, "test rooted state")
                    mutated_state = dict(rooted_state)
                    mutated_state["manifest_sha256"] = "d" * 64
                    subject.atomic_state_json(state_path, mutated_state)
                    with self.assertRaises(subject.CampaignError):
                        subject.seal_holdout(argparse.Namespace(
                            result_dir=result, phase="h1"))
                    subject.common.atomic_write(state_path, state_bytes)
                    for field, value in (
                            ("root_file", "holdout/wrong.json"),
                            ("root_file_sha256", "e" * 64)):
                        mutated_state = dict(rooted_state)
                        mutated_state[field] = value
                        subject.atomic_state_json(state_path, mutated_state)
                        with self.assertRaises(subject.CampaignError):
                            subject.acquire_holdout(acquire)
                        with self.assertRaises(subject.CampaignError):
                            subject.seal_holdout(argparse.Namespace(
                                result_dir=result, phase="h2"))
                        subject.common.atomic_write(state_path, state_bytes)
                    (result / "development/frozen_table.json").write_text(
                        "mutated\n", encoding="ascii")
                    with self.assertRaises(subject.CampaignError):
                        subject.acquire_holdout(acquire)


if __name__ == "__main__":
    unittest.main()
