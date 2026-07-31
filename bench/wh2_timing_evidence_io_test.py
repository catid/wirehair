#!/usr/bin/env python3
"""Adversarial unit tests for WH2 timing evidence filesystem helpers."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import stat
import tempfile
import threading
import unittest
from unittest import mock

import wh2_timing_evidence_io as evidence


class SnapshotTests(unittest.TestCase):
    def test_snapshot_and_streaming_digest_are_exact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "artifact.bin"
            payload = bytes(range(251)) * 10000
            path.write_bytes(payload)
            self.assertEqual(
                evidence.stable_file_snapshot(path, max_bytes=len(payload)),
                payload,
            )
            self.assertEqual(
                evidence.stable_file_sha256(path, max_bytes=len(payload)),
                hashlib.sha256(payload).hexdigest(),
            )

    def test_symlink_fifo_and_directory_are_rejected_without_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            regular = root / "regular"
            regular.write_bytes(b"safe")
            symlink = root / "symlink"
            symlink.symlink_to(regular)
            fifo = root / "fifo"
            os.mkfifo(str(fifo), 0o600)
            directory = root / "directory"
            directory.mkdir()
            for path in (symlink, fifo, directory):
                with self.subTest(path=path.name):
                    with self.assertRaises(evidence.EvidenceIOError):
                        evidence.stable_file_snapshot(path, max_bytes=16)

    def test_symlinked_parent_component_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real = root / "real"
            real.mkdir()
            (real / "artifact").write_bytes(b"safe")
            alias = root / "alias"
            alias.symlink_to(real, target_is_directory=True)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.stable_file_snapshot(alias / "artifact", max_bytes=4)

    def test_hardlink_is_optional_but_rejected_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            first.write_bytes(b"payload")
            os.link(str(first), str(second))
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.stable_file_snapshot(first, max_bytes=7)
            self.assertEqual(
                evidence.stable_file_snapshot(
                    first, max_bytes=7, require_unique=False
                ),
                b"payload",
            )

    def test_oversize_is_rejected_for_snapshot_and_digest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "large"
            path.write_bytes(b"123456789")
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.stable_file_snapshot(path, max_bytes=8)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.stable_file_sha256(path, max_bytes=8)

    def test_path_replacement_between_passes_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "artifact"
            retired = root / "retired"
            replacement = root / "replacement"
            path.write_bytes(b"first")
            replacement.write_bytes(b"other")
            original_read = evidence._read_bounded_pass
            calls = []

            def replacing_read(*args, **kwargs):
                result = original_read(*args, **kwargs)
                calls.append(None)
                if len(calls) == 1:
                    os.rename(str(path), str(retired))
                    os.rename(str(replacement), str(path))
                return result

            with mock.patch.object(
                evidence, "_read_bounded_pass", side_effect=replacing_read
            ):
                with self.assertRaises(evidence.EvidenceIOError):
                    evidence.stable_file_snapshot(path, max_bytes=16)

    def test_supplied_error_type_is_raised(self) -> None:
        class ControllerError(RuntimeError):
            pass

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "directory"
            path.mkdir()
            with self.assertRaises(ControllerError):
                evidence.stable_file_snapshot(
                    path, max_bytes=1, error_type=ControllerError
                )


class CanonicalJSONTests(unittest.TestCase):
    def _write(self, root: Path, raw: bytes) -> Path:
        path = root / "record.json"
        path.write_bytes(raw)
        return path

    def test_exact_canonical_object_and_raw_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = evidence.canonical_json_bytes({"z": 2, "a": [1, True]})
            path = self._write(root, raw)
            value, loaded_raw = evidence.load_canonical_object_snapshot(
                path, max_bytes=1024
            )
            self.assertEqual(value, {"a": [1, True], "z": 2})
            self.assertEqual(loaded_raw, raw)
            self.assertEqual(
                evidence.load_canonical_object(path, max_bytes=1024), value
            )

    def test_duplicate_keys_are_rejected_at_every_depth(self) -> None:
        cases = (
            b'{"a":1,"a":2}\n',
            b'{"a":{"b":1,"b":2}}\n',
        )
        for raw in cases:
            with self.subTest(raw=raw):
                with tempfile.TemporaryDirectory() as temporary:
                    path = self._write(Path(temporary), raw)
                    with self.assertRaises(evidence.EvidenceIOError):
                        evidence.load_canonical_object(path, max_bytes=1024)

    def test_nonfinite_constants_and_overflow_are_rejected(self) -> None:
        for token in (b"NaN", b"Infinity", b"-Infinity", b"1e999"):
            with self.subTest(token=token):
                with tempfile.TemporaryDirectory() as temporary:
                    raw = b'{"value":' + token + b"}\n"
                    path = self._write(Path(temporary), raw)
                    with self.assertRaises(evidence.EvidenceIOError):
                        evidence.load_canonical_object(path, max_bytes=1024)

    def test_non_object_and_noncanonical_encoding_are_rejected(self) -> None:
        for raw in (b"[]\n", b'{"b":2, "a":1}\n', b'{"a":"\xc3\xa9"}\n'):
            with self.subTest(raw=raw):
                with tempfile.TemporaryDirectory() as temporary:
                    path = self._write(Path(temporary), raw)
                    with self.assertRaises(evidence.EvidenceIOError):
                        evidence.load_canonical_object(path, max_bytes=1024)


class PublicationTests(unittest.TestCase):
    def test_immutable_file_is_unique_readonly_and_no_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "artifact"
            evidence.publish_immutable_file(target, b"first")
            self.assertEqual(target.read_bytes(), b"first")
            metadata = target.stat()
            self.assertEqual(metadata.st_nlink, 1)
            self.assertEqual(stat.S_IMODE(metadata.st_mode), 0o400)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.publish_immutable_file(target, b"second")
            self.assertEqual(target.read_bytes(), b"first")

    def test_existing_symlink_is_not_followed_or_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            victim = root / "victim"
            victim.write_bytes(b"victim")
            target = root / "artifact"
            target.symlink_to(victim)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.publish_immutable_file(target, b"attack")
            self.assertTrue(target.is_symlink())
            self.assertEqual(victim.read_bytes(), b"victim")

    def test_concurrent_publishers_have_exactly_one_winner(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "artifact"
            barrier = threading.Barrier(2)
            outcomes = []
            outcome_lock = threading.Lock()

            def publisher(payload: bytes) -> None:
                barrier.wait()
                try:
                    evidence.publish_immutable_file(target, payload)
                except evidence.EvidenceIOError:
                    result = ("error", payload)
                except BaseException as exc:
                    result = ("unexpected", exc)
                else:
                    result = ("success", payload)
                with outcome_lock:
                    outcomes.append(result)

            threads = [
                threading.Thread(target=publisher, args=(b"alpha",)),
                threading.Thread(target=publisher, args=(b"bravo",)),
            ]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(5)
                self.assertFalse(thread.is_alive())
            successes = [payload for outcome, payload in outcomes if outcome == "success"]
            unexpected = [value for outcome, value in outcomes if outcome == "unexpected"]
            self.assertEqual(len(outcomes), 2)
            self.assertEqual(unexpected, [])
            self.assertEqual(len(successes), 1)
            self.assertEqual(target.read_bytes(), successes[0])
            self.assertEqual(target.stat().st_nlink, 1)
            self.assertEqual(list(root.glob(".*.partial")), [])

    def test_directory_publication_is_noreplace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "staging"
            target = root / "published"
            source.mkdir(mode=0o700)
            os.chmod(str(source), 0o500)
            evidence.publish_directory_noreplace(source, target)
            self.assertFalse(source.exists())
            self.assertTrue(target.is_dir())
            self.assertEqual(stat.S_IMODE(target.stat().st_mode), 0o500)

            second_source = root / "second-staging"
            second_source.mkdir(mode=0o700)
            os.chmod(str(second_source), 0o500)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.publish_directory_noreplace(second_source, target)
            self.assertTrue(second_source.is_dir())
            self.assertEqual(stat.S_IMODE(second_source.stat().st_mode), 0o500)

    def test_directory_publication_rejects_writable_or_symlink_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            writable = root / "writable"
            writable.mkdir(mode=0o700)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.publish_directory_noreplace(
                    writable, root / "published-writable"
                )
            real = root / "real"
            real.mkdir(mode=0o700)
            os.chmod(str(real), 0o500)
            alias = root / "alias"
            alias.symlink_to(real, target_is_directory=True)
            with self.assertRaises(evidence.EvidenceIOError):
                evidence.publish_directory_noreplace(alias, root / "published-alias")


class DirectoryAndLockTests(unittest.TestCase):
    def test_private_directory_policy_and_held_descriptor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sealed = root / "sealed"
            sealed.mkdir(mode=0o700)
            os.chmod(str(sealed), 0o500)
            with evidence.held_private_directory(
                sealed, require_writable=False
            ) as descriptor:
                self.assertTrue(stat.S_ISDIR(os.fstat(descriptor).st_mode))

            os.chmod(str(sealed), 0o550)
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.held_private_directory(
                    sealed, require_writable=False
                ):
                    pass
            os.chmod(str(sealed), 0o700)
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.held_private_directory(
                    sealed, require_writable=False
                ):
                    pass

    def test_nonblocking_global_flock_rejects_contention_then_releases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            lock = Path(temporary) / "controller.lock"
            with evidence.nonblocking_global_flock(lock):
                with self.assertRaises(evidence.EvidenceIOError):
                    with evidence.nonblocking_global_flock(lock):
                        pass
            with evidence.nonblocking_global_flock(lock):
                self.assertEqual(stat.S_IMODE(lock.stat().st_mode), 0o600)

    def test_lock_symlink_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            victim = root / "victim"
            victim.write_bytes(b"")
            lock = root / "controller.lock"
            lock.symlink_to(victim)
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.nonblocking_global_flock(lock):
                    pass

    def test_lock_fifo_and_hardlink_are_rejected_without_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fifo = root / "fifo.lock"
            os.mkfifo(str(fifo), 0o600)
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.nonblocking_global_flock(fifo):
                    pass
            first = root / "first.lock"
            second = root / "second.lock"
            first.write_bytes(b"")
            os.chmod(str(first), 0o600)
            os.link(str(first), str(second))
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.nonblocking_global_flock(first):
                    pass

    def test_lock_path_replacement_is_rejected_even_during_unwind(self) -> None:
        class BodyError(RuntimeError):
            pass

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lock = root / "controller.lock"
            retired = root / "retired.lock"
            replacement = root / "replacement.lock"
            replacement.write_bytes(b"")
            os.chmod(str(replacement), 0o600)
            with self.assertRaises(evidence.EvidenceIOError):
                with evidence.nonblocking_global_flock(lock):
                    os.rename(str(lock), str(retired))
                    os.rename(str(replacement), str(lock))
                    raise BodyError("simulated controller failure")


if __name__ == "__main__":
    unittest.main()
