#!/usr/bin/env python3
"""Warning-strict tests for the source-pinned rv4a verifier."""

import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "wh2_rv4a_parallel_verify_under_test",
    HERE / "wh2_rv4a_parallel_verify.py",
)
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


def _binding(path, byte_limit):
    path = Path(path).absolute()
    descriptor = path.stat()
    payload = path.read_bytes()
    return {
        "path": str(path),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
        "mode": descriptor.st_mode,
        "size": descriptor.st_size,
        "mtime_ns": descriptor.st_mtime_ns,
        "ctime_ns": descriptor.st_ctime_ns,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "byte_limit": byte_limit,
    }


def _manifest(bindings):
    return {
        "schema": verifier.CAMPAIGN_SCHEMA,
        "frozen_roster_sha256": "1" * 64,
        "policy_sha256": "2" * 64,
        "result_free_plan_sha256": "3" * 64,
        "frozen_roster": {
            "native_protocol": verifier.REPAIRTIMING_PROTOCOL,
            "native_schema": verifier.REPAIRTIMING_SCHEMA,
        },
        "runtime_bindings": bindings,
    }


class PinRecordTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        binary = self.root / "benchmark"
        binary.write_bytes(b"binary")
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
        sources = {}
        for name, source_class in verifier.SOURCE_CLASSES.items():
            path = self.root / name
            path.write_bytes(f"{name}\n".encode("ascii"))
            sources[name] = _binding(
                path, verifier.FILE_BYTE_LIMITS[source_class])
        thermal = self.root / "thermal.csv"
        thermal.write_bytes(b"thermal\n")
        thermal_stat = thermal.stat()
        self.bindings = {
            "benchmark": _binding(
                binary, verifier.FILE_BYTE_LIMITS["binary"]),
            "sources": sources,
            "thermal": {
                "path": str(thermal.absolute()),
                "device": thermal_stat.st_dev,
                "inode": thermal_stat.st_ino,
                "mode": thermal_stat.st_mode,
            },
        }
        self.manifest = _manifest(self.bindings)

    def tearDown(self):
        self.temporary.cleanup()

    def test_pin_record_covers_every_runtime_file(self):
        pin = verifier.make_pin_record(self.manifest)
        self.assertEqual(
            set(pin["runtime_files"]),
            set(verifier.RUNTIME_FILE_NAMES),
        )
        self.assertEqual(
            pin["runtime_files"]["benchmark"]["class"], "binary")
        self.assertEqual(
            pin["runtime_files"]["native_repairtiming"]["class"],
            "native",
        )
        self.assertEqual(
            pin["runtime_files"]["parser"]["class"], "python")
        self.assertEqual(
            pin["runtime_files"]["build_configuration"]["class"],
            "configuration",
        )
        verifier._validate_pin_record(pin, self.manifest)

    def test_wrong_runner_native_binary_and_build_pins_fail(self):
        original = verifier.make_pin_record(self.manifest)
        for name in (
            "runner",
            "native_repairtiming",
            "benchmark",
            "codec_build",
            "build_configuration",
        ):
            with self.subTest(name=name):
                changed = copy.deepcopy(original)
                changed["runtime_files"][name]["sha256"] = "f" * 64
                with self.assertRaises(verifier.PinError):
                    verifier._validate_pin_record(
                        changed, self.manifest)

        changed = copy.deepcopy(original)
        changed["thermal"]["inode"] += 1
        with self.assertRaises(verifier.PinError):
            verifier._validate_pin_record(changed, self.manifest)

    def test_runtime_rehash_rejects_binary_and_source_mutation(self):
        verifier._verify_runtime_files(self.bindings)
        for name, path in (
            ("benchmark", Path(self.bindings["benchmark"]["path"])),
            (
                "native",
                Path(self.bindings["sources"][
                    "native_repairtiming"]["path"]),
            ),
            (
                "build",
                Path(self.bindings["sources"]["codec_build"]["path"]),
            ),
        ):
            with self.subTest(name=name):
                verifier._verify_runtime_files(self.bindings)
                original = path.read_bytes()
                path.write_bytes(original + b"x")
                try:
                    with self.assertRaises(verifier.PinError):
                        verifier._verify_runtime_files(self.bindings)
                finally:
                    path.write_bytes(original)
                    if name == "benchmark":
                        self.bindings["benchmark"] = _binding(
                            path,
                            verifier.FILE_BYTE_LIMITS["binary"],
                        )
                    elif name == "native":
                        self.bindings["sources"][
                            "native_repairtiming"] = _binding(
                                path,
                                verifier.FILE_BYTE_LIMITS["native"],
                            )
                    else:
                        self.bindings["sources"]["codec_build"] = _binding(
                            path,
                            verifier.FILE_BYTE_LIMITS["build"],
                        )

    def test_stable_source_snapshot_rejects_symlink(self):
        target = self.root / "target.py"
        target.write_bytes(b"VALUE = 1\n")
        link = self.root / "link.py"
        link.symlink_to(target)
        with self.assertRaises(verifier.PinError):
            verifier._stable_file(
                link.absolute(), byte_limit=1024, retain_payload=True)

    def test_source_forced_cli_executes_exact_verifier(self):
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                str(HERE / "wh2_rv4a_parallel_verify.py"),
                "--help",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertIn(b"preflight-pin", completed.stdout)
        self.assertIn(b"verify", completed.stdout)
        self.assertEqual(completed.stderr, b"")


class PinnedReplayTests(unittest.TestCase):
    def test_bounded_replay_and_exact_synthetic_completion(self):
        runner_path = Path("/tmp/rv4a-synthetic/bench/runner.py")
        parser_path = Path("/tmp/rv4a-synthetic/tools/parser.py")
        context_path = Path("/tmp/rv4a-synthetic/tools/context.py")
        bindings = {
            "benchmark": {
                "path": "/tmp/rv4a-synthetic/bench-bin",
                "device": 1,
                "inode": 1,
                "mode": stat.S_IFREG | 0o755,
                "size": 1,
                "mtime_ns": 1,
                "ctime_ns": 1,
                "sha256": "a" * 64,
                "byte_limit": verifier.FILE_BYTE_LIMITS["binary"],
            },
            "sources": {},
            "thermal": {
                "path": "/tmp/rv4a-synthetic/thermal",
                "device": 1,
                "inode": 2,
                "mode": stat.S_IFREG | 0o644,
            },
        }
        for index, (name, source_class) in enumerate(
                verifier.SOURCE_CLASSES.items(), start=10):
            path = Path(f"/tmp/rv4a-synthetic/{name}")
            if name == "runner":
                path = runner_path
            elif name == "parser":
                path = parser_path
            elif name == "context_tool":
                path = context_path
            elif name == "parallel_verifier":
                path = Path(verifier.__file__).absolute()
            bindings["sources"][name] = {
                "path": str(path),
                "device": 1,
                "inode": index,
                "mode": stat.S_IFREG | 0o644,
                "size": 1,
                "mtime_ns": 1,
                "ctime_ns": 1,
                "sha256": f"{index:064x}",
                "byte_limit": verifier.FILE_BYTE_LIMITS[source_class],
            }
        manifest = _manifest(bindings)
        pin = verifier.make_pin_record(manifest)
        pin_sha256 = verifier.sha256_bytes(
            verifier.canonical_json_bytes(pin))
        manifest["preflight_pin_sha256"] = pin_sha256
        snapshots = {
            name: {
                "path": item["path"],
                "sha256": item["sha256"],
                "payload": b"",
            }
            for name, item in {
                "benchmark": bindings["benchmark"],
                **bindings["sources"],
            }.items()
        }

        campaign = types.ModuleType("synthetic_campaign")
        campaign.__file__ = str(runner_path)
        campaign.__rv4a_source_sha256__ = \
            bindings["sources"]["runner"]["sha256"]
        campaign.CAMPAIGN_SCHEMA = verifier.CAMPAIGN_SCHEMA
        campaign.REQUIRED_REPAIRTIMING_PROTOCOL = \
            verifier.REPAIRTIMING_PROTOCOL
        campaign.REQUIRED_REPAIRTIMING_SCHEMA = verifier.REPAIRTIMING_SCHEMA
        parser = types.ModuleType("repair_timing_codec")
        parser.__file__ = str(parser_path)
        parser.__rv4a_source_sha256__ = \
            bindings["sources"]["parser"]["sha256"]
        context = types.ModuleType("peel_codec")
        context.__file__ = str(context_path)
        context.__rv4a_source_sha256__ = \
            bindings["sources"]["context_tool"]["sha256"]
        campaign.repair = parser
        campaign.peel_codec = context
        calls = []
        campaign._validate_manifest = lambda value: calls.append(
            ("manifest", value))
        campaign.check_runtime_bindings = \
            lambda value, full_hash: calls.append(
                ("bindings", value, full_hash))
        synthetic_completion = {
            "completion_sha256": "4" * 64,
            "receipts_replayed": 1188,
            "exact_reproduction": True,
        }

        def verify_campaign(directory, *, replay_workers):
            calls.append(("verify", Path(directory), replay_workers))
            return synthetic_completion

        campaign.verify_campaign = verify_campaign

        def read_canonical(path, *, byte_limit):
            if Path(path).name == "manifest.json":
                return manifest, "5" * 64
            return pin, pin_sha256

        modules = {
            "repair_timing_codec": parser,
            "peel_codec": context,
        }
        with (
            mock.patch.object(
                verifier, "_read_canonical", side_effect=read_canonical
            ),
            mock.patch.object(
                verifier,
                "_verify_runtime_files",
                side_effect=(snapshots, snapshots),
            ) as rehashed,
            mock.patch.object(
                verifier,
                "_module_from_runner_snapshot",
                return_value=campaign,
            ),
            mock.patch.object(
                verifier,
                "_EXECUTING_VERIFIER_SOURCE_SHA256",
                bindings["sources"]["parallel_verifier"]["sha256"],
            ),
            mock.patch.dict(sys.modules, modules),
        ):
            result = verifier.verify_with_pins(
                "/tmp/rv4a-synthetic/campaign",
                "/tmp/rv4a-synthetic/pins.json",
                pin_sha256,
                replay_workers=7,
            )
        self.assertEqual(result["verified"], synthetic_completion)
        self.assertEqual(rehashed.call_count, 2)
        self.assertIn(
            (
                "verify",
                Path("/tmp/rv4a-synthetic/campaign").absolute(),
                7,
            ),
            calls,
        )
        self.assertEqual(
            sum(item[0] == "bindings" for item in calls), 2)

    def test_worker_bound_and_trusted_digest_fail_closed(self):
        for workers in (0, 33, True):
            with self.subTest(workers=workers):
                with self.assertRaises(verifier.PinError):
                    verifier.verify_with_pins(
                        "/missing",
                        "/missing",
                        "a" * 64,
                        replay_workers=workers,
                    )
        with self.assertRaises(verifier.PinError):
            verifier.verify_with_pins(
                "/missing", "/missing", "not-a-digest")


if __name__ == "__main__":
    unittest.main()
