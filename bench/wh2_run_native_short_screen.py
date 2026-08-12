#!/usr/bin/env python3
"""Run the fixed, native WH2 development timing short screen.

This is deliberately a narrow controller for the frozen development domain.
It qualifies, emits, and freezes timing traces before any result job, keeps one
persistent native worker on each selected physical core, traverses the full
timing cohort domain once per independent round with one homogeneous eight-job
wave and barrier at a time, and delegates publication to the fail-closed native
evidence assembler.  Raw-v3 recovery is produced only by the separate recovery
campaign runner; this controller never emits recovery evidence.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import errno
import hashlib
import io
import math
import os
from pathlib import Path
import platform
import re
import selectors
import signal
import socket
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional
from typing import Sequence, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api


DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v1"
DESCRIPTION_FIELDS = frozenset((
    "schema", "source_git_commit", "binary_sha256", "arms",
))
DESCRIPTION_ARM_FIELDS = frozenset((
    "arm", "codec", "arm_descriptor_sha256",
))
EXPECTED_ARMS = (
    ("wirehair2_head", "wirehair2_certified"),
    ("wirehair1", "wirehair1"),
    ("wirehair2_dense_two07_basis_v1", "wirehair2_experiment"),
)
TIMING_CANDIDATE_DESCRIPTOR_SHA256 = (
    "9527f200ad38c7eec6502b2f768fdd67"
    "b92787fb227eed3d7616274ffc2df388"
)
SUMMARY_SCHEMA = "wirehair.wh2.native-short-screen-run.v2"
SUMMARY_FIELDS = frozenset((
    "schema", "status", "output_dir", "source_git_commit",
    "contract_sha256", "worker_binary_sha256", "controller_cpu",
    "worker_cpus", "qualification_worker_cpus",
    "timing_qualification_map_sha256",
    "timing_qualification_execution_receipt_sha256",
    "timing_qualified_domain_sha256", "qualification_attempt_count",
    "qualification_retried_cell_count", "qualification_max_retry_offset",
    "qualification_sum_retry_offsets", "timing_records",
    "timing_trace_manifest_sha256", "timing_freeze_sha256",
    "timing_architecture_artifact_sha256", "timing_result_sha256",
    "timing_execution_receipt_sha256", "timing_validator_summary_sha256",
    "thermal_samples", "cpu_tctl_max_millic", "dimm_max_millic",
    "timing_thermal_samples", "timing_cpu_tctl_max_millic",
    "timing_dimm_max_millic", "qualification_thermal_samples",
    "qualification_cpu_tctl_max_millic", "qualification_dimm_max_millic",
    "overall_cpu_tctl_max_millic", "overall_dimm_max_millic",
    "summary_sha256",
))
MAX_WALL_SECONDS = 7200.0
MAX_RESPONSE_BYTES = 1024 * 1024
TIMING_WORKER_COUNT = 8
ZERO_SHA256 = "0" * 64
CPU_TOKEN = re.compile(r"(?:0|[1-9][0-9]*)\Z")
CACHE_INDEX_TOKEN = re.compile(r"index(?:0|[1-9][0-9]*)\Z")
POSITIVE_TOKEN = re.compile(r"[1-9][0-9]*\Z")
MAX_LINUX_PID = (1 << 31) - 1
MAX_LINUX_PID_TEXT = str(MAX_LINUX_PID)
CPU_SYSFS_ROOT = Path("/sys/devices/system/cpu")

# The producer admits every retry offset in [0,255] independently for all
# 2,304 frozen qualification base cells.  Bound each canonical record before
# publication, then derive the completed-file limits from that exact maximum
# cardinality.  This keeps the producer and strict reopener in agreement even
# for the 589,824-row qualification envelope instead of relying on the much
# smaller common case (one row per base cell).
FROZEN_TIMING_BASE_CELLS = 2304
MAX_QUALIFICATION_ATTEMPTS_PER_CELL = 256
MAX_QUALIFICATION_ATTEMPTS = (
    FROZEN_TIMING_BASE_CELLS * MAX_QUALIFICATION_ATTEMPTS_PER_CELL)
FROZEN_TIMING_PANEL_COUNT = 11
FROZEN_TIMING_RECORDS = FROZEN_TIMING_BASE_CELLS * FROZEN_TIMING_PANEL_COUNT
MAX_QUALIFICATION_NATIVE_RECORD_BYTES = 2048
MAX_QUALIFICATION_AUDIT_RECORD_BYTES = 1024
MAX_TIMING_NATIVE_RECORD_BYTES = 4096
MAX_TIMING_RESULT_RECORD_BYTES = 4096
MAX_SMALL_COMPLETED_ARTIFACT_BYTES = 16 * 1024 * 1024
COMPLETED_ARTIFACT_BYTE_LIMITS = {
    "run-summary.json": MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "timing-freeze.json": MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "timing-traces.jsonl":
        FROZEN_TIMING_BASE_CELLS * MAX_QUALIFICATION_AUDIT_RECORD_BYTES,
    "timing-native-results.jsonl":
        FROZEN_TIMING_RECORDS * MAX_TIMING_NATIVE_RECORD_BYTES,
    "timing-results.jsonl":
        FROZEN_TIMING_RECORDS * MAX_TIMING_RESULT_RECORD_BYTES,
    "timing-execution.json": MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "sampler-attestation.json": MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "timing-qualification-map.json": MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "timing-qualification-audit.jsonl":
        MAX_QUALIFICATION_ATTEMPTS * MAX_QUALIFICATION_AUDIT_RECORD_BYTES,
    "timing-qualification-native.jsonl":
        MAX_QUALIFICATION_ATTEMPTS * MAX_QUALIFICATION_NATIVE_RECORD_BYTES,
    "qualification-sampler-attestation.json":
        MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
    "timing-qualification-execution.json":
        MAX_SMALL_COMPLETED_ARTIFACT_BYTES,
}
MAX_COMPLETED_ARTIFACT_BYTES = max(COMPLETED_ARTIFACT_BYTE_LIMITS.values())
MAX_COMPLETED_BUNDLE_BYTES = sum(COMPLETED_ARTIFACT_BYTE_LIMITS.values())
COMPLETED_SOURCE_NAMES = {
    "summary": "run-summary.json",
    "freeze": "timing-freeze.json",
    "trace": "timing-traces.jsonl",
    "native": "timing-native-results.jsonl",
    "result": "timing-results.jsonl",
    "receipt": "timing-execution.json",
    "sampler": "sampler-attestation.json",
    "qualification_map": "timing-qualification-map.json",
    "qualification_audit": "timing-qualification-audit.jsonl",
    "qualification_native": "timing-qualification-native.jsonl",
    "qualification_sampler": "qualification-sampler-attestation.json",
    "qualification_receipt": "timing-qualification-execution.json",
}
COMPLETED_DEPENDENCY_NAMES = {
    key: name for key, name in COMPLETED_SOURCE_NAMES.items()
    if key != "summary"
}


class RunnerError(RuntimeError):
    """The real native short screen cannot continue safely."""


def fail(message: str) -> None:
    raise RunnerError(message)


CleanupFailure = Tuple[str, BaseException]


def _raise_after_cleanup(
        primary: Optional[BaseException],
        failures: Sequence[CleanupFailure],
        terminal_diagnostics: Sequence[str] = (),
        ) -> None:
    """Report deferred cleanup failures without swallowing control flow."""
    if not failures and not terminal_diagnostics:
        return
    details = []
    if primary is not None:
        details.append("primary failure: {}".format(primary))
    details.extend(
        "{}: {}: {}".format(label, type(failure).__name__, failure)
        for label, failure in failures)
    details.extend(terminal_diagnostics)
    diagnostic = RunnerError("; ".join(details))

    # KeyboardInterrupt, SystemExit, and arbitrary application-defined
    # BaseException subclasses are control flow, not diagnostics.  Cleanup
    # must be exhaustive before one is re-raised, and an already-active
    # control-flow exception remains authoritative.
    control_flow: Optional[BaseException] = None
    if primary is not None and not isinstance(primary, Exception):
        control_flow = primary
    else:
        control_flow = next((
            failure for _label, failure in failures
            if not isinstance(failure, Exception)), None)
    if control_flow is not None:
        raise control_flow from diagnostic
    if primary is not None:
        raise diagnostic from primary
    if failures:
        raise diagnostic from failures[0][1]
    raise diagnostic


def _drain_descriptor_closes(
        descriptors: Iterable[Tuple[str, int]],
        primary: Optional[BaseException],
        ) -> None:
    """Close every descriptor, ignoring only close(2) state uncertainty."""
    failures: List[CleanupFailure] = []
    for label, descriptor in descriptors:
        if descriptor < 0:
            continue
        try:
            os.close(descriptor)
        except OSError:
            # POSIX leaves the descriptor state ambiguous after close(2)
            # reports an error.  Retrying can close a reused descriptor.
            pass
        except BaseException as cleanup:
            failures.append((label, cleanup))
    _raise_after_cleanup(primary, failures)


def _remaining(deadline: float, context: str) -> float:
    value = deadline - time.monotonic()
    if value <= 0.0:
        fail("native short-screen hard wall expired during {}".format(context))
    return value


@contextmanager
def _hard_wall(seconds: float):
    """Interrupt synchronous work and expose one pre-commit teardown hook."""
    if not hasattr(signal, "setitimer") or not hasattr(signal, "ITIMER_REAL"):
        fail("native short-screen hard wall requires POSIX setitimer")
    started = time.monotonic()

    def expired(_signum: int, _frame: Any) -> None:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        fail("{}-second native short-screen hard wall expired".format(
            int(seconds)))

    try:
        old_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, expired)
    except (OSError, ValueError) as exc:
        fail("cannot install native short-screen hard wall: {}".format(exc))
    try:
        old_timer = signal.setitimer(signal.ITIMER_REAL, seconds)
    except (OSError, ValueError) as exc:
        signal.signal(signal.SIGALRM, old_handler)
        fail("cannot install native short-screen hard wall: {}".format(exc))
    finished = False

    def finish() -> None:
        nonlocal finished
        if finished:
            return
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)
        if old_timer[0] > 0.0:
            remaining = max(0.000001,
                            old_timer[0] - (time.monotonic() - started))
            signal.setitimer(signal.ITIMER_REAL, remaining, old_timer[1])
        finished = True

    try:
        yield finish
    finally:
        if not finished:
            finish()


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and contract_api.SHA256.fullmatch(value) \
        is not None


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            while True:
                block = source.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
    except OSError as exc:
        fail("cannot hash {}: {}".format(path, exc))
    return digest.hexdigest()


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | \
        getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(str(path), flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _open_publication_guard(path: Path) -> Tuple[int, os.stat_result]:
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    if nofollow == 0:
        fail("publication requires O_NOFOLLOW")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | nofollow
    descriptor = -1
    try:
        descriptor = os.open(str(path), flags)
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode):
            fail("staged publication artifact must be regular")
        return descriptor, info
    except OSError as exc:
        if descriptor >= 0:
            os.close(descriptor)
        fail("cannot guard staged artifact {}: {}".format(path, exc))
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        raise


def _publish_staged(
        staged: Path, destination: Path,
        expected_identity: Optional[Tuple[int, int]] = None) -> None:
    if destination.exists() or destination.is_symlink():
        fail("refusing to replace existing artifact {}".format(destination))
    guard_fd, staged_info = _open_publication_guard(staged)
    directory_fd = -1
    try:
        directory_fd = os.open(
            str(destination.parent),
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
            getattr(os, "O_CLOEXEC", 0))
        if (expected_identity is not None and expected_identity !=
                (staged_info.st_dev, staged_info.st_ino)):
            fail("staged publication artifact identity changed")
        try:
            os.link(
                "/proc/self/fd/{}".format(guard_fd), destination.name,
                dst_dir_fd=directory_fd,
                follow_symlinks=True)
            published_info = os.stat(
                destination.name, dir_fd=directory_fd,
                follow_symlinks=False)
            if (published_info.st_dev != staged_info.st_dev or
                    published_info.st_ino != staged_info.st_ino):
                fail("published artifact identity changed during link")
            os.fsync(directory_fd)
            if not _unlink_if_identity(
                    staged, staged_info.st_dev, staged_info.st_ino):
                fail("staged publication artifact changed before cleanup")
            os.fsync(directory_fd)
            return
        except FileExistsError:
            fail("refusing to replace existing artifact {}".format(destination))
        except OSError as exc:
            fail("cannot publish {}: {}".format(destination, exc))
        except BaseException:
            # Visible artifacts are immutable dependencies.  A later strict
            # receipt/run-summary is the only completion marker, so never risk
            # deleting a concurrently replaced destination on failure.
            raise
    finally:
        if directory_fd >= 0:
            os.close(directory_fd)
        os.close(guard_fd)


def _unlink_if_identity(path: Path, device: int, inode: int) -> bool:
    """Best-effort cleanup for staging names in the private output directory."""
    try:
        info = os.lstat(str(path))
    except FileNotFoundError:
        return False
    except OSError:
        return False
    if info.st_dev != device or info.st_ino != inode:
        return False
    try:
        os.unlink(str(path))
        return True
    except OSError:
        return False


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    descriptor, raw = tempfile.mkstemp(
        prefix="." + path.name + ".", suffix=".tmp", dir=str(path.parent))
    staged = Path(raw)
    staged_identity: Optional[Tuple[int, int]] = None
    publication_identity: Optional[Tuple[int, int]] = None
    publication_guard = -1
    committed = False
    output = None
    try:
        staged_info = os.fstat(descriptor)
        staged_identity = (staged_info.st_dev, staged_info.st_ino)
        output = os.fdopen(descriptor, "wb")
        descriptor = -1
        with output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        publication_guard, info = _open_publication_guard(staged)
        publication_identity = (info.st_dev, info.st_ino)
        _publish_staged(staged, path, publication_identity)
        if _hash_file(path) != hashlib.sha256(data).hexdigest():
            fail("published artifact changed before commit")
        committed = True
        guard = publication_guard
        try:
            os.close(guard)
        finally:
            publication_guard = -1
        return
    except BaseException:
        if output is not None:
            try:
                output.close()
            except BaseException:
                pass
            finally:
                descriptor = -1
        elif descriptor >= 0:
            try:
                os.close(descriptor)
            except OSError:
                pass
            finally:
                descriptor = -1
        if not committed and publication_identity is not None:
            _unlink_if_identity(
                staged, publication_identity[0], publication_identity[1])
        elif not committed and staged_identity is not None:
            _unlink_if_identity(
                staged, staged_identity[0], staged_identity[1])
        elif not committed:
            try:
                staged.unlink()
            except FileNotFoundError:
                pass
        if publication_guard >= 0:
            guard = publication_guard
            try:
                os.close(guard)
            except OSError:
                pass
            finally:
                publication_guard = -1
        raise


def _atomic_write_object(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_write_bytes(
        path, (contract_api.canonical_json(value) + "\n").encode("utf-8"))


def _open_completion_parent(
        path: Path, descriptor_holder: Optional[List[int]] = None,
        ) -> Tuple[int, Tuple[int, int]]:
    if (descriptor_holder is not None and
            (not isinstance(descriptor_holder, list) or
             descriptor_holder != [-1])):
        fail("completion-marker parent descriptor holder is invalid")
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    directory = getattr(os, "O_DIRECTORY", 0)
    if nofollow == 0 or directory == 0:
        fail("completion-marker publication requires O_NOFOLLOW/O_DIRECTORY")
    descriptor = -1
    try:
        descriptor = os.open(
            str(path.parent), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory)
        if descriptor_holder is not None:
            descriptor_holder[0] = descriptor
        info = os.fstat(descriptor)
        if not stat.S_ISDIR(info.st_mode):
            fail("completion-marker parent must be a real directory")
        identity = (info.st_dev, info.st_ino)
        _verify_completed_directory_path(path.parent, identity)
        return descriptor, identity
    except OSError as exc:
        primary = RunnerError(
            "cannot open completion-marker parent {}: {}".format(
                path.parent, exc))
        if descriptor_holder is not None:
            descriptor_holder[0] = -1
        _drain_descriptor_closes((
            ("completion-parent descriptor cleanup failed", descriptor),
        ), primary)
        raise primary from exc
    except BaseException as primary:
        if descriptor_holder is not None:
            descriptor_holder[0] = -1
        _drain_descriptor_closes((
            ("completion-parent descriptor cleanup failed", descriptor),
        ), primary)
        raise
    raise AssertionError("unreachable")


def _require_completion_parent_descriptor(
        descriptor: int, identity: Tuple[int, int]) -> None:
    try:
        info = os.fstat(descriptor)
    except OSError as exc:
        fail("completion-marker parent descriptor is invalid: {}".format(exc))
    if (not stat.S_ISDIR(info.st_mode) or
            (info.st_dev, info.st_ino) != identity):
        fail("completion-marker parent descriptor identity changed")


def _require_completion_marker_absent(
        path: Path, directory_fd: int) -> None:
    """Require the public completion name to be absent without following it."""
    try:
        os.stat(path.name, dir_fd=directory_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as exc:
        fail("cannot inspect completion marker {}: {}".format(path, exc))
    fail("refusing to replace existing artifact {}".format(path))


def _read_unnamed_completion_marker(
        descriptor: int, byte_limit: int,
        ) -> Tuple[bytes, CompletedFingerprint]:
    """Stably reread a bounded, still-unnamed regular file descriptor."""
    if type(descriptor) is not int or descriptor < 0:
        fail("unnamed completion-marker descriptor is invalid")
    if type(byte_limit) is not int or byte_limit <= 0:
        fail("unnamed completion-marker byte cap is invalid")
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 0 or
                before.st_size > byte_limit):
            fail("unnamed completion marker is not a bounded private file")
        chunks = []
        size = 0
        while True:
            block = os.pread(
                descriptor, min(1024 * 1024, byte_limit + 1 - size), size)
            if not block:
                break
            size += len(block)
            if size > byte_limit:
                fail("unnamed completion marker exceeds its byte cap")
            chunks.append(block)
        after = os.fstat(descriptor)
    except OSError as exc:
        fail("cannot read unnamed completion marker: {}".format(exc))
    before_fingerprint = _completed_fingerprint(before)
    after_fingerprint = _completed_fingerprint(after)
    data = b"".join(chunks)
    if (before_fingerprint != after_fingerprint or
            len(data) != before.st_size):
        fail("unnamed completion marker changed while it was being read")
    return data, after_fingerprint


def _open_unnamed_completion_marker(
        directory_fd: int, data: bytes, descriptor_holder: List[int],
        ) -> CompletedFingerprint:
    """Create and durably fill one O_TMPFILE inode with no directory name."""
    limit = COMPLETED_ARTIFACT_BYTE_LIMITS["run-summary.json"]
    if not isinstance(data, bytes) or not data or len(data) > limit:
        fail("completion marker exceeds its publication byte cap")
    if (not isinstance(descriptor_holder, list) or
            descriptor_holder != [-1]):
        fail("unnamed completion-marker descriptor holder is invalid")
    tmpfile = getattr(os, "O_TMPFILE", 0)
    cloexec = getattr(os, "O_CLOEXEC", 0)
    if tmpfile == 0 or cloexec == 0:
        fail("completion-marker publication requires O_TMPFILE/O_CLOEXEC")
    descriptor = -1
    try:
        descriptor = os.open(
            ".", tmpfile | os.O_RDWR | cloexec, 0o600,
            dir_fd=directory_fd)
        # Populate caller-owned cleanup state before this helper can return.
        # An interrupt in the return/store handoff therefore cannot leak the
        # anonymous inode for the rest of a long-lived controller process.
        descriptor_holder[0] = descriptor
        created = os.fstat(descriptor)
        if not stat.S_ISREG(created.st_mode) or created.st_nlink != 0:
            fail("O_TMPFILE completion marker is not a private regular file")
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                fail("short write while staging unnamed completion marker")
            offset += written
        os.fsync(descriptor)
        observed, fingerprint = _read_unnamed_completion_marker(
            descriptor, limit)
        if (observed != data or
                hashlib.sha256(observed).digest() !=
                hashlib.sha256(data).digest()):
            fail("unnamed completion marker differs from prospective bytes")
        return fingerprint
    except OSError as exc:
        fail("cannot stage unnamed completion marker: {}".format(exc))
    finally:
        if descriptor >= 0 and sys.exc_info()[0] is not None:
            primary = sys.exc_info()[1]
            descriptor_holder[0] = -1
            _drain_descriptor_closes((
                ("unnamed completion-marker descriptor cleanup failed",
                 descriptor),
            ), primary)
    raise AssertionError("unreachable")


def _link_unnamed_completion_marker(
        path: Path, descriptor: int, directory_fd: int,
        directory_identity: Tuple[int, int]) -> Tuple[int, int]:
    """Publish once; after link(2), every outcome is a committed-marker ACK."""
    _require_completion_parent_descriptor(directory_fd, directory_identity)
    _verify_completed_directory_path(path.parent, directory_identity)
    _require_completion_marker_absent(path, directory_fd)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 0:
            fail("unnamed completion marker is no longer private")
        os.link(
            "/proc/self/fd/{}".format(descriptor), path.name,
            dst_dir_fd=directory_fd, follow_symlinks=True)
    except FileExistsError:
        fail("refusing to replace existing artifact {}".format(path))
    except OSError as exc:
        # link(2) may have succeeded before an injected signal or an unusual
        # client-side error became visible.  The public name is never removed:
        # a successful link is the transaction's publish-wins linearization.
        fail("cannot publish completion marker {}: {}".format(path, exc))
    # Directory fsync is an acknowledgement/durability step, not a rollback
    # boundary.  If it fails, the valid linked marker deliberately remains.
    os.fsync(directory_fd)
    return info.st_dev, info.st_ino


class AtomicLineSink:
    """A result stream that is visible only after every row is complete."""

    def __init__(
            self, destination: Path,
            maximum_bytes: Optional[int] = None) -> None:
        if maximum_bytes is not None and (
                type(maximum_bytes) is not int or maximum_bytes <= 0):
            fail("native result stream has an invalid byte cap")
        digest = hashlib.sha256()
        descriptor, raw = tempfile.mkstemp(
            prefix="." + destination.name + ".", suffix=".tmp",
            dir=str(destination.parent))
        staged = Path(raw)
        staged_identity: Optional[Tuple[int, int]] = None
        output = None
        try:
            info = os.fstat(descriptor)
            staged_identity = (info.st_dev, info.st_ino)
            output = os.fdopen(descriptor, "wb")
            descriptor = -1
            self.destination = destination
            self.staged = staged
            self.output = output
            self.staged_identity = staged_identity
            self.sha256 = digest
            self.maximum_bytes = maximum_bytes
            self.bytes_written = 0
            self.published = False
        except BaseException:
            if output is not None:
                try:
                    output.close()
                except BaseException:
                    pass
                finally:
                    descriptor = -1
            elif descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
                finally:
                    descriptor = -1
            if staged_identity is not None:
                _unlink_if_identity(
                    staged, staged_identity[0], staged_identity[1])
            else:
                try:
                    staged.unlink()
                except FileNotFoundError:
                    pass
            raise

    def write(self, line: bytes) -> None:
        if self.output.closed or not line.endswith(b"\n"):
            fail("cannot write an invalid native result line")
        if (self.maximum_bytes is not None and
                len(line) > self.maximum_bytes - self.bytes_written):
            fail("native result stream exceeds its completed-artifact byte cap")
        self.output.write(line)
        self.sha256.update(line)
        self.bytes_written += len(line)

    def publish(self) -> None:
        if self.published:
            fail("native result stream was published twice")
        self.output.flush()
        os.fsync(self.output.fileno())
        info = os.fstat(self.output.fileno())
        if (info.st_dev, info.st_ino) != self.staged_identity:
            fail("native result staging identity changed")
        committed = False
        try:
            _publish_staged(
                self.staged, self.destination, self.staged_identity)
            if _hash_file(self.destination) != self.sha256.hexdigest():
                fail("published native result stream changed before commit")
            self.published = True
            committed = True
            self.output.close()
            return
        except BaseException:
            if not committed:
                _unlink_if_identity(
                    self.staged, info.st_dev, info.st_ino)
                self.published = False
            if not self.output.closed:
                self.output.close()
            raise

    def abort(self) -> None:
        if not self.output.closed:
            self.output.close()
        _unlink_if_identity(
            self.staged, self.staged_identity[0], self.staged_identity[1])


def _parse_canonical_line(data: bytes, context: str) -> Mapping[str, Any]:
    if not data.endswith(b"\n") or data.count(b"\n") != 1:
        fail("{} must be exactly one newline-terminated JSON object".format(
            context))
    logical = data[:-1]
    try:
        value = contract_api._load_json_bytes(logical, context)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if not isinstance(value, dict) or \
            contract_api.canonical_json(value).encode("utf-8") != logical:
        fail("{} is not canonical JSON".format(context))
    return value


def parse_cpu_list(text: str) -> List[int]:
    """Parse a comma-separated, strictly increasing canonical CPU list."""
    if not isinstance(text, str) or not text:
        fail("--cpus must be a nonempty comma-separated CPU list")
    tokens = text.split(",")
    if any(CPU_TOKEN.fullmatch(token) is None for token in tokens):
        fail("--cpus contains a noncanonical CPU ID")
    values = [int(token) for token in tokens]
    if values != sorted(set(values)):
        fail("--cpus must be strictly increasing and unique")
    return values


def parse_pid_list(text: str) -> List[int]:
    """Parse an optional strictly increasing list of controller-owned PIDs."""
    if not isinstance(text, str) or not text:
        fail("--pause-sibling-pids must be a nonempty comma-separated PID list")
    tokens = text.split(",")
    if any(POSITIVE_TOKEN.fullmatch(token) is None for token in tokens):
        fail("--pause-sibling-pids contains a noncanonical PID")
    if any(len(token) > len(MAX_LINUX_PID_TEXT) or
           (len(token) == len(MAX_LINUX_PID_TEXT) and
            token > MAX_LINUX_PID_TEXT)
           for token in tokens):
        fail("--pause-sibling-pids contains a PID outside Linux pid_t")
    values = [int(token) for token in tokens]
    if values != sorted(set(values)):
        fail("--pause-sibling-pids must be strictly increasing and unique")
    return values


def _cpu_topology(cpu: int, sysfs_root: Path) -> Tuple[int, int]:
    topology = sysfs_root / "cpu{}".format(cpu) / "topology"
    values = []
    for name in ("physical_package_id", "core_id"):
        try:
            text = (topology / name).read_text(encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} topology: {}".format(cpu, exc))
        if CPU_TOKEN.fullmatch(text) is None:
            fail("CPU {} has a noncanonical {}".format(cpu, name))
        values.append(int(text))
    return values[0], values[1]


def _cpu_llc_identity(cpu: int, sysfs_root: Path) -> Tuple[int, int]:
    """Read the highest data/unified cache level and ID without assuming index3."""
    cache_root = sysfs_root / "cpu{}".format(cpu) / "cache"
    try:
        entries = sorted(
            (path for path in cache_root.iterdir()
             if CACHE_INDEX_TOKEN.fullmatch(path.name) is not None),
            key=lambda path: int(path.name[5:]))
    except OSError as exc:
        fail("cannot enumerate CPU {} cache topology: {}".format(cpu, exc))
    if not entries:
        fail("CPU {} has no indexed cache topology".format(cpu))
    data_caches: List[Tuple[int, int]] = []
    for entry in entries:
        try:
            level_text = (entry / "level").read_text(
                encoding="ascii").strip()
            cache_type = (entry / "type").read_text(
                encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} cache topology: {}".format(cpu, exc))
        if POSITIVE_TOKEN.fullmatch(level_text) is None or \
                cache_type not in ("Data", "Instruction", "Unified"):
            fail("CPU {} has malformed cache level/type metadata".format(cpu))
        if cache_type == "Instruction":
            continue
        try:
            cache_id_text = (entry / "id").read_text(
                encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} cache ID: {}".format(cpu, exc))
        if CPU_TOKEN.fullmatch(cache_id_text) is None:
            fail("CPU {} has a noncanonical cache ID".format(cpu))
        data_caches.append((int(level_text), int(cache_id_text)))
    if not data_caches:
        fail("CPU {} has no data/unified cache topology".format(cpu))
    highest_level = max(level for level, _ in data_caches)
    ids = {cache_id for level, cache_id in data_caches
           if level == highest_level}
    if len(ids) != 1:
        fail("CPU {} has ambiguous last-level cache IDs".format(cpu))
    return highest_level, next(iter(ids))


def _select_llc_spread_workers(
        candidates: Mapping[Tuple[int, int], int],
        llc_by_core: Mapping[
            Tuple[int, int], Tuple[int, int, int]]) -> List[int]:
    """Pick eight cores with maximum LLC coverage, then round-robin fill."""
    grouped: Dict[Tuple[int, int, int], List[int]] = {}
    for core, cpu in candidates.items():
        grouped.setdefault(llc_by_core[core], []).append(cpu)
    for values in grouped.values():
        values.sort()
    selected: List[int] = []
    group_ids = sorted(grouped)
    while len(selected) < TIMING_WORKER_COUNT:
        progressed = False
        for cache_id in group_ids:
            values = grouped[cache_id]
            if values:
                selected.append(values.pop(0))
                progressed = True
                if len(selected) == TIMING_WORKER_COUNT:
                    break
        if not progressed:
            break
    if len(selected) != TIMING_WORKER_COUNT:
        fail("fewer than eight eligible physical worker cores are available")
    return sorted(selected)


def select_cpu_layout(
        sampler_cpu: int, explicit_workers: Optional[Sequence[int]] = None,
        explicit_controller: Optional[int] = None,
        affinity: Optional[Iterable[int]] = None,
        sysfs_root: Path = CPU_SYSFS_ROOT,
        ) -> Tuple[List[int], int]:
    """Select eight LLC-spread physical cores plus an isolated controller."""
    if type(sampler_cpu) is not int or sampler_cpu < 0:
        fail("sampler CPU must be a nonnegative integer")
    if affinity is None:
        try:
            allowed = sorted(os.sched_getaffinity(0))
        except (AttributeError, OSError) as exc:
            fail("cannot inspect controller CPU affinity: {}".format(exc))
    else:
        affinity_values = list(affinity)
        if any(type(cpu) is not int or cpu < 0 for cpu in affinity_values):
            fail("controller CPU affinity contains an invalid CPU ID")
        allowed = sorted(set(affinity_values))
    if not allowed:
        fail("controller CPU affinity is empty")
    sampler_core = _cpu_topology(sampler_cpu, sysfs_root)
    topology = {cpu: _cpu_topology(cpu, sysfs_root) for cpu in allowed}
    # Cache IDs need not be globally unique across packages or cache levels.
    # Bind both dimensions before comparing last-level cache groups.
    llc_by_cpu = {
        cpu: (topology[cpu][0],) + _cpu_llc_identity(cpu, sysfs_root)
        for cpu in allowed
    }
    allowed_set = set(allowed)

    by_core: Dict[Tuple[int, int], int] = {}
    llc_by_core: Dict[Tuple[int, int], Tuple[int, int, int]] = {}
    for cpu in allowed:
        core = topology[cpu]
        cache_id = llc_by_cpu[cpu]
        if core in llc_by_core and llc_by_core[core] != cache_id:
            fail("SMT siblings disagree on their last-level cache ID")
        llc_by_core[core] = cache_id
        by_core[core] = min(cpu, by_core.get(core, cpu))

    controller_core: Optional[Tuple[int, int]] = None
    if explicit_controller is not None:
        if type(explicit_controller) is not int or explicit_controller < 0 or \
                explicit_controller not in allowed_set:
            fail("controller CPU is outside current sched affinity")
        controller_core = topology[explicit_controller]
        if explicit_controller == sampler_cpu or controller_core == sampler_core:
            fail("controller CPU overlaps the sampler's physical core")

    worker_candidates = {
        core: cpu for core, cpu in by_core.items()
        if core != sampler_core and core != controller_core
    }

    if explicit_workers is not None:
        workers = list(explicit_workers)
        if (workers != sorted(set(workers)) or
                len(workers) != TIMING_WORKER_COUNT):
            fail("explicit worker CPUs must be eight sorted unique IDs")
        if any(type(cpu) is not int or cpu < 0 for cpu in workers):
            fail("explicit worker CPUs contain an invalid CPU ID")
        if not set(workers).issubset(allowed_set):
            fail("explicit worker CPUs are outside current sched affinity")
        worker_cores = [topology[cpu] for cpu in workers]
        if len(worker_cores) != len(set(worker_cores)):
            fail("explicit worker CPUs contain SMT siblings")
        if sampler_cpu in workers or sampler_core in worker_cores:
            fail("worker CPUs overlap the sampler's physical core")
        if controller_core in worker_cores:
            fail("controller CPU shares a physical core with a worker")
    else:
        workers = _select_llc_spread_workers(
            worker_candidates, llc_by_core)
        worker_cores = [topology[cpu] for cpu in workers]

    available_llc_ids = {
        llc_by_core[core] for core in worker_candidates
    }
    required_llc_coverage = min(
        TIMING_WORKER_COUNT, len(available_llc_ids))
    observed_llc_coverage = len({
        llc_by_core[topology[cpu]] for cpu in workers
    })
    if observed_llc_coverage != required_llc_coverage:
        fail("worker CPU roster does not maximize available LLC coverage")

    if explicit_controller is not None:
        controller_cpu = explicit_controller
    else:
        worker_core_set = set(worker_cores)
        controller_candidates = {
            core: cpu for core, cpu in by_core.items()
            if core != sampler_core and core not in worker_core_set
        }
        if controller_candidates:
            controller_cpu = max(controller_candidates.values())
            controller_core = topology[controller_cpu]
        else:
            fail("no dedicated controller core remains outside frozen workers")

    if (len(workers) != TIMING_WORKER_COUNT or
            controller_cpu in workers or controller_core == sampler_core or
            controller_core in {topology[cpu] for cpu in workers}):
        fail("controller/worker/sampler CPU layout is not physically isolated")
    return workers, controller_cpu


def _proc_race(exc: BaseException) -> bool:
    return isinstance(exc, OSError) and exc.errno in (errno.ENOENT, errno.ESRCH)


def _read_proc_task_state(path: Path) -> str:
    try:
        raw = path.read_text(encoding="ascii")
    except (OSError, UnicodeError) as exc:
        if _proc_race(exc):
            raise FileNotFoundError(str(path)) from exc
        fail("cannot read task state {}: {}".format(path, exc))
    close = raw.rfind(")")
    fields = raw[close + 1:].split() if close >= 0 else []
    if len(fields) < 20 or len(fields[0]) != 1:
        fail("task state {} is malformed".format(path))
    return fields[0]


def _process_task_states(pid: int, proc_root: Path = Path("/proc")) \
        -> List[Tuple[int, str]]:
    task_root = proc_root / str(pid) / "task"
    try:
        entries = sorted(
            (entry for entry in task_root.iterdir()
             if entry.name.isascii() and entry.name.isdecimal()),
            key=lambda entry: int(entry.name))
    except OSError as exc:
        if _proc_race(exc):
            return []
        fail("cannot enumerate process {} tasks: {}".format(pid, exc))
    result = []
    for entry in entries:
        try:
            result.append((
                int(entry.name), _read_proc_task_state(entry / "stat")))
        except FileNotFoundError:
            continue
    return result


def _process_start_ticks(pid: int, proc_root: Path = Path("/proc")) -> int:
    stat_path = proc_root / str(pid) / "stat"
    try:
        raw = stat_path.read_text(encoding="ascii")
    except (OSError, UnicodeError) as exc:
        if _proc_race(exc):
            raise ProcessLookupError(pid) from exc
        fail("cannot read process {} identity: {}".format(pid, exc))
    close = raw.rfind(")")
    fields = raw[close + 1:].split() if close >= 0 else []
    if len(fields) < 20 or not fields[19].isascii() or \
            not fields[19].isdecimal():
        fail("process {} identity is malformed".format(pid))
    value = int(fields[19])
    if value <= 0:
        fail("process {} has an invalid start time".format(pid))
    return value


def _scan_exact_affinity_tasks(
        exclusive_cpus: Set[int], proc_root: Path = Path("/proc"),
        ) -> List[Mapping[str, Any]]:
    """Find executable tasks whose entire scheduler affinity is T-owned CPUs."""
    if (not exclusive_cpus or
            any(type(cpu) is not int or cpu < 0 for cpu in exclusive_cpus)):
        fail("SMT isolation scan requires nonempty valid exclusive CPUs")
    try:
        processes = sorted(
            (entry for entry in proc_root.iterdir()
             if entry.name.isascii() and entry.name.isdecimal()),
            key=lambda entry: int(entry.name))
    except OSError as exc:
        fail("cannot enumerate /proc for SMT isolation: {}".format(exc))
    found: List[Mapping[str, Any]] = []
    for process in processes:
        pid = int(process.name)
        try:
            # Kernel threads have no executable link and are not schedulable
            # user workloads.  Every executable process must remain inspectable.
            info = (process / "exe").stat()
            if not stat.S_ISREG(info.st_mode):
                continue
            start_ticks = _process_start_ticks(pid, proc_root)
        except (FileNotFoundError, ProcessLookupError):
            continue
        except OSError as exc:
            if _proc_race(exc):
                continue
            fail("cannot inspect process {} executable: {}".format(pid, exc))
        process_found: List[Mapping[str, Any]] = []
        for tid, state in _process_task_states(pid, proc_root):
            try:
                affinity = set(os.sched_getaffinity(tid))
            except (AttributeError, OSError) as exc:
                if _proc_race(exc):
                    continue
                fail("cannot inspect task {} affinity: {}".format(tid, exc))
            if affinity and affinity.issubset(exclusive_cpus):
                process_found.append({
                    "pid": pid,
                    "tid": tid,
                    "process_start_ticks": start_ticks,
                    "cpu_affinity": sorted(affinity),
                    "state": state,
                })
        try:
            final_start_ticks = _process_start_ticks(pid, proc_root)
        except ProcessLookupError:
            # A task that was exact-pinned to the protected roster was real
            # contamination at the instant we observed it.  Process exit does
            # not erase that evidence; retaining it keeps the boundary scan
            # fail-closed.  Races involving processes with no matching task
            # remain irrelevant and may be skipped.
            found.extend(process_found)
            continue
        if final_start_ticks != start_ticks:
            fail("process {} changed identity during SMT isolation scan".
                 format(pid))
        found.extend(process_found)
    return found


class ManagedPausedSiblingProcesses:
    """Own SIGSTOP/SIGCONT cleanup for explicitly declared load workers."""

    def __init__(
            self, records: Sequence[Mapping[str, Any]] = ()) -> None:
        self._records = [dict(record) for record in records]
        self._pidfds: Dict[int, io.FileIO] = {}
        self._resumed = False

    def _drain_resumed_pidfds(self) -> None:
        """Attempt each idempotent close and retire only proven-closed owners."""
        failures: List[CleanupFailure] = []
        for pid in sorted(list(self._pidfds)):
            handle = self._pidfds.get(pid)
            if not isinstance(handle, io.FileIO):
                failures.append((
                    "managed sibling PID {} pidfd owner is invalid".format(
                        pid),
                    RunnerError("managed sibling pidfd owner type changed")))
                continue
            try:
                try:
                    handle.close()
                finally:
                    # An exception before close leaves the same open object in
                    # the map.  An exception during/after FileIO.close leaves a
                    # closed object whose repeated close is idempotent.  Only
                    # remove the exact object after its C owner says closed.
                    if (handle.closed and
                            self._pidfds.get(pid) is handle):
                        self._pidfds.pop(pid)
            except BaseException as cleanup:
                failures.append((
                    "managed sibling PID {} pidfd cleanup failed".format(pid),
                    cleanup))
        _raise_after_cleanup(None, failures)

    def pause(self, pids: Sequence[int], sibling_cpus: Sequence[int],
              deadline: float) -> None:
        if self._records or self._pidfds or self._resumed:
            fail("managed sibling pause owner is not pristine")
        if (any(type(pid) is not int or pid <= 0 for pid in pids) or
                list(pids) != sorted(set(pids))):
            fail("managed sibling PIDs must be sorted unique positive integers")
        siblings = set(sibling_cpus)
        if pids and not siblings:
            fail("cannot pause sibling load when the timing roster has no SMT")
        if pids and (not callable(getattr(os, "pidfd_open", None)) or
                     not callable(getattr(signal, "pidfd_send_signal", None))):
            fail("managed sibling pause requires Linux pidfd signal support")
        owner_uid = os.geteuid()
        sudo_uid = os.environ.get("SUDO_UID")
        if sudo_uid is not None:
            if (not sudo_uid.isascii() or not sudo_uid.isdecimal() or
                    (len(sudo_uid) > 1 and sudo_uid.startswith("0"))):
                fail("SUDO_UID is not a canonical user ID")
            owner_uid = int(sudo_uid)
        try:
            for pid in pids:
                if pid == os.getpid():
                    fail("controller cannot pause itself as sibling load")
                try:
                    # The only raw-integer interval is pidfd_open's CALL return
                    # to FileIO's CALL entry, before any SIGSTOP ownership is
                    # published.  Once FileIO returns, C-level descriptor state
                    # makes every later close idempotent across Python async
                    # instruction boundaries.  Do not retry a constructor
                    # failure with raw os.close: adoption state is ambiguous.
                    pidfd = io.FileIO(
                        os.pidfd_open(pid, 0), mode="rb", closefd=True)
                    self._pidfds[pid] = pidfd
                    process_info = Path("/proc/{}".format(pid)).stat()
                    executable_info = Path(
                        "/proc/{}/exe".format(pid)).stat()
                    start_ticks = _process_start_ticks(pid)
                    tasks = _process_task_states(pid)
                except (OSError, ProcessLookupError) as exc:
                    fail("cannot inspect managed sibling PID {}: {}".format(
                        pid, exc))
                if (process_info.st_uid != owner_uid or
                        not stat.S_ISREG(executable_info.st_mode) or
                        not tasks):
                    fail("managed sibling PID is not a live controller-owned process")
                process_affinity: Optional[List[int]] = None
                for tid, state in tasks:
                    if state in ("T", "t"):
                        fail("managed sibling PID was already stopped")
                    try:
                        affinity = set(os.sched_getaffinity(tid))
                    except (AttributeError, OSError) as exc:
                        fail("cannot inspect managed sibling PID affinity: {}".
                             format(exc))
                    if not affinity or not affinity.issubset(siblings):
                        fail("managed sibling PID is not pinned only to T siblings")
                    canonical_affinity = sorted(affinity)
                    if process_affinity is None:
                        process_affinity = canonical_affinity
                    elif canonical_affinity != process_affinity:
                        fail("managed sibling PID tasks have different affinities")
                if process_affinity is None:
                    fail("managed sibling PID has no inspectable tasks")
                try:
                    signal.pidfd_send_signal(pidfd.fileno(), 0, None, 0)
                except OSError as exc:
                    fail("managed sibling PID exited before stop: {}".format(
                        exc))
                record = {
                    "pid": pid,
                    "process_start_ticks": start_ticks,
                    "cpu_affinity": process_affinity,
                }
                # Publish provisional ownership before SIGSTOP: every failure
                # after the signal must route this exact identity through the
                # SIGCONT rollback path.
                self._records.append(record)
                signal.pidfd_send_signal(
                    pidfd.fileno(), signal.SIGSTOP, None, 0)
                while True:
                    try:
                        current_ticks = _process_start_ticks(pid)
                    except ProcessLookupError:
                        fail("managed sibling PID exited while stopping")
                    if current_ticks != start_ticks:
                        fail("managed sibling PID changed identity while stopping")
                    current = _process_task_states(pid)
                    if current and all(state in ("T", "t")
                                       for _, state in current):
                        break
                    _remaining(deadline, "stopping managed sibling load")
                    time.sleep(0.01)
                if _process_start_ticks(pid) != start_ticks:
                    fail("managed sibling PID changed identity while stopping")
            return
        except BaseException as primary:
            try:
                # A stop-time hard-wall failure must not donate its expired
                # deadline to rollback.  The manager already owns every PID
                # that may have observed SIGSTOP, so reserve a fresh SIGCONT
                # grace period before propagating the primary failure.
                self.resume(max(time.monotonic() + 5.0, deadline))
            except BaseException as cleanup:
                _raise_after_cleanup(
                    primary,
                    [("managed sibling load rollback failed", cleanup)])
            raise

    @property
    def records(self) -> List[Mapping[str, Any]]:
        return [dict(record) for record in self._records]

    @property
    def cleanup_complete(self) -> bool:
        """Return whether no stopped-process or pidfd cleanup remains."""
        return self._resumed and not self._pidfds

    def resume(self, deadline: float) -> None:
        deferred_control: Optional[BaseException] = None
        if not self._resumed:
            pending = {record["pid"]: record for record in self._records}
            last_failures: Dict[int, BaseException] = {}
            attempt = 0
            while pending:
                attempt += 1
                for pid, record in list(pending.items()):
                    try:
                        pidfd = self._pidfds.get(pid)
                        if (not isinstance(pidfd, io.FileIO) or
                                pidfd.closed):
                            fail("managed sibling PID lost its exact signal handle")
                        try:
                            signal.pidfd_send_signal(
                                pidfd.fileno(), signal.SIGCONT, None, 0)
                        except OSError as exc:
                            if _proc_race(exc):
                                pending.pop(pid)
                                last_failures.pop(pid, None)
                                continue
                            raise
                        while True:
                            try:
                                current_ticks = _process_start_ticks(pid)
                            except ProcessLookupError:
                                break
                            if current_ticks != record["process_start_ticks"]:
                                break
                            states = _process_task_states(pid)
                            if not states or all(state not in ("T", "t")
                                                 for _, state in states):
                                break
                            _remaining(
                                deadline, "resuming managed sibling load")
                            time.sleep(0.01)
                        pending.pop(pid)
                        last_failures.pop(pid, None)
                    except BaseException as cleanup:
                        last_failures[pid] = cleanup
                        if (deferred_control is None and
                                not isinstance(cleanup, Exception)):
                            deferred_control = cleanup
                if not pending:
                    break
                now = time.monotonic()
                # Every exact handle receives at least two attempts, even if
                # the cleanup deadline is already expired.
                if attempt >= 2 and now >= deadline:
                    break
                try:
                    time.sleep(min(0.01, max(0.0, deadline - now)))
                except BaseException as cleanup:
                    if (deferred_control is None and
                            not isinstance(cleanup, Exception)):
                        deferred_control = cleanup

            if pending:
                failures = [
                    ("managed sibling PID {} resume failed".format(pid),
                     last_failures[pid])
                    for pid in sorted(pending)
                ]
                _raise_after_cleanup(deferred_control, failures)
                raise AssertionError("unreachable")

            # Once every stopped record is confirmed continued (or exited),
            # this one state bit switches the same durable handle map from
            # signal ownership to idempotent close ownership.
            transition_attempts = 0
            while not self._resumed:
                try:
                    self._resumed = True
                except BaseException as cleanup:
                    transition_attempts += 1
                    if (deferred_control is None and
                            not isinstance(cleanup, Exception)):
                        deferred_control = cleanup
                    if (not isinstance(cleanup, Exception) and
                            transition_attempts < 8):
                        continue
                    raise

        drain_failures: List[CleanupFailure] = []
        drain_attempts = 0
        while self._pidfds and drain_attempts < 8:
            drain_attempts += 1
            try:
                self._drain_resumed_pidfds()
            except BaseException as cleanup:
                drain_failures.append((
                    "managed sibling pidfd drain attempt {} failed".format(
                        drain_attempts), cleanup))
            if not self._pidfds:
                break
            if drain_attempts >= 2 and time.monotonic() >= deadline:
                break
        if self._pidfds:
            drain_failures.append((
                "managed sibling pidfd drain did not converge",
                RunnerError(
                    "managed sibling cleanup still owns {} pidfd handles "
                    "after {} drain attempts".format(
                        len(self._pidfds), drain_attempts))))
        if drain_failures:
            _raise_after_cleanup(
                deferred_control, drain_failures)
        if deferred_control is not None:
            raise deferred_control


class SiblingIsolationGuard:
    """Fail every T wave boundary if a foreign exact-affinity task appears."""

    def __init__(
            self, worker_cpus: Sequence[int], controller_cpu: int,
            paused_processes: Sequence[Mapping[str, Any]] = (),
            scan: Callable[[Set[int]], List[Mapping[str, Any]]] =
                _scan_exact_affinity_tasks,
            sysfs_root: Path = CPU_SYSFS_ROOT) -> None:
        self.worker_cpus = list(worker_cpus)
        self.controller_cpu = controller_cpu
        self.protected_cpus = sorted(self.worker_cpus + [controller_cpu])
        self.sibling_cpus = native_api.timing_sibling_cpus(
            self.protected_cpus, sysfs_root)
        self.exclusive_cpus = set(
            self.protected_cpus + self.sibling_cpus)
        self.paused_processes = [dict(value) for value in paused_processes]
        self._paused = {
            (value["pid"], value["process_start_ticks"]): value
            for value in self.paused_processes
        }
        self._workers: Dict[Tuple[int, int], int] = {}
        self._checks: List[Mapping[str, Any]] = []
        self._scan = scan

    def bind_workers(self, workers: Sequence["PersistentWorker"]) -> None:
        expected = set(self.worker_cpus)
        if (len(workers) != len(self.worker_cpus) or
                {worker.cpu for worker in workers} != expected or
                any(worker.start_ticks <= 0 for worker in workers)):
            fail("cannot bind an incomplete timing worker roster to SMT guard")
        self._workers = {
            (worker.process.pid, worker.start_ticks): worker.cpu
            for worker in workers
        }
        if len(self._workers) != len(workers):
            fail("timing worker roster reuses a process identity")

    def check(self, stage: str) -> None:
        if not isinstance(stage, str) or not stage or len(stage) > 128:
            fail("SMT isolation check stage is invalid")
        foreign = []
        for task in self._scan(self.exclusive_cpus):
            identity = (task["pid"], task["process_start_ticks"])
            # The controller leader is trusted and separately singleton-pinned
            # before T.  Do not exempt arbitrary threads that happen to share
            # its PID: an unexpected helper pinned to an exclusive sibling is
            # just as capable of contaminating timing as another process.
            if task["pid"] == os.getpid() and task["tid"] == os.getpid():
                if task["cpu_affinity"] == [self.controller_cpu]:
                    continue
                foreign.append(task)
                continue
            worker_cpu = self._workers.get(identity)
            if (worker_cpu is not None and task["tid"] == task["pid"] and
                    task["cpu_affinity"] == [worker_cpu]):
                continue
            paused = self._paused.get(identity)
            if (paused is not None and task["state"] in ("T", "t") and
                    task["cpu_affinity"] == paused["cpu_affinity"]):
                continue
            foreign.append(task)
        now = time.monotonic_ns()
        self._checks.append({
            "ordinal": len(self._checks),
            "stage": stage,
            "monotonic_ns": now,
            "foreign_task_count": len(foreign),
        })
        if foreign:
            task = foreign[0]
            fail("timing SMT isolation found foreign exact-affinity task "
                 "pid={} tid={} cpus={}".format(
                     task["pid"], task["tid"],
                     ",".join(str(cpu) for cpu in task["cpu_affinity"])))

    def attestation(self) -> Mapping[str, Any]:
        if len(self._checks) < 2 or any(
                check["foreign_task_count"] != 0 for check in self._checks):
            fail("cannot attest incomplete timing SMT isolation")
        checks = [dict(value) for value in self._checks]
        return {
            "schema": native_api.SIBLING_ISOLATION_SCHEMA,
            "policy": "exact-affinity-wave-boundary-v1",
            "worker_cpus": list(self.worker_cpus),
            "controller_cpu": self.controller_cpu,
            "protected_cpus": list(self.protected_cpus),
            "sibling_cpus": list(self.sibling_cpus),
            "paused_processes": self.paused_processes,
            "checks": checks,
            "check_count": len(checks),
            "checks_sha256": contract_api.sha256_json(checks),
            "first_check_monotonic_ns": checks[0]["monotonic_ns"],
            "last_check_monotonic_ns": checks[-1]["monotonic_ns"],
            "terminal_status": "clean",
        }


def _pin_controller(cpu: int) -> None:
    try:
        os.sched_setaffinity(0, {cpu})
        observed = os.sched_getaffinity(0)
    except (AttributeError, OSError) as exc:
        fail("cannot pin controller to CPU {}: {}".format(cpu, exc))
    if observed != {cpu}:
        fail("controller did not establish singleton CPU affinity")


def _restore_controller_affinity(affinity: Set[int]) -> None:
    try:
        os.sched_setaffinity(0, affinity)
        observed = os.sched_getaffinity(0)
    except (AttributeError, OSError) as exc:
        fail("cannot restore controller CPU affinity: {}".format(exc))
    if observed != affinity:
        fail("controller CPU affinity did not restore exactly")


def _run_command(argv: Sequence[str], deadline: float,
                 context: str) -> bytes:
    try:
        result = subprocess.run(
            list(argv), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False,
            timeout=_remaining(deadline, context))
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("{} failed: {}".format(context, exc))
    if result.returncode != 0:
        diagnostic = result.stderr.decode("utf-8", "replace").strip()
        fail("{} exited {}{}".format(
            context, result.returncode,
            ": " + diagnostic if diagnostic else ""))
    if result.stderr:
        fail("{} produced unexpected stderr: {}".format(
            context, result.stderr.decode("utf-8", "replace").strip()))
    return result.stdout


def describe_worker(worker_path: Path, deadline: float) -> Mapping[str, Any]:
    try:
        resolved = worker_path.resolve(strict=True)
        info = resolved.stat()
    except OSError as exc:
        fail("cannot resolve native worker {}: {}".format(worker_path, exc))
    if not stat.S_ISREG(info.st_mode) or not os.access(str(resolved), os.X_OK):
        fail("native worker must be an executable regular file")
    binary_sha256 = _hash_file(resolved)
    value = _parse_canonical_line(
        _run_command([str(resolved), "--describe"], deadline,
                     "native worker description"),
        "native worker description")
    if set(value) != DESCRIPTION_FIELDS or \
            value.get("schema") != DESCRIPTION_SCHEMA or \
            not isinstance(value.get("source_git_commit"), str) or \
            re.fullmatch(r"[0-9a-f]{40}",
                         value["source_git_commit"]) is None or \
            value.get("binary_sha256") != binary_sha256:
        fail("native worker description does not bind its executable")
    arms = value.get("arms")
    if not isinstance(arms, list) or len(arms) != len(EXPECTED_ARMS):
        fail("native worker description has an invalid arm roster")
    for index, expected in enumerate(EXPECTED_ARMS):
        arm = arms[index]
        if not isinstance(arm, dict) or set(arm) != DESCRIPTION_ARM_FIELDS or \
                arm.get("arm") != expected[0] or \
                arm.get("codec") != expected[1] or \
                not _is_sha256(arm.get("arm_descriptor_sha256")):
            fail("native worker description has an invalid arm at index {}".format(
                index))
        if index == 2 and arm["arm_descriptor_sha256"] != \
                TIMING_CANDIDATE_DESCRIPTOR_SHA256:
            fail("native worker description has the wrong timing candidate")
    if len({arm["arm_descriptor_sha256"] for arm in arms}) != len(arms):
        fail("native worker description reuses an arm descriptor")
    result = dict(value)
    result["resolved_path"] = str(resolved)
    return result


def _require_worker_source_commit(
        description: Mapping[str, Any], source_commit: str) -> None:
    if description.get("source_git_commit") != source_commit:
        fail("native worker was not built from the clean source HEAD")


def _git_head(deadline: float, repo_root: Optional[Path] = None) -> str:
    """Return HEAD only for a clean, stable codec-source snapshot."""
    root = Path(__file__).resolve().parent.parent if repo_root is None \
        else repo_root
    commands = (
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
        # Any tracked change can invalidate the source attribution (including
        # a build/configuration input outside the obvious codec directories).
        ["git", "status", "--porcelain=v1", "--untracked-files=no"],
        # Untracked experiment/result artifacts are common and harmless.  An
        # untracked build source in one of these locations is not.
        ["git", "ls-files", "--others", "--exclude-standard", "--",
         "CMakeLists.txt", "codec", "bench", "include", "cmake",
         ":(top,glob)*.c", ":(top,glob)*.cc", ":(top,glob)*.cpp",
         ":(top,glob)*.h", ":(top,glob)*.hpp", ":(top,glob)*.inc"],
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
    )
    outputs = []
    try:
        for command in commands:
            result = subprocess.run(
                command, cwd=str(root), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                timeout=_remaining(deadline, "checking source commit"))
            if result.returncode != 0:
                diagnostic = result.stderr.decode("utf-8", "replace").strip()
                fail("cannot establish clean source commit{}".format(
                    ": " + diagnostic if diagnostic else ""))
            outputs.append(result.stdout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot establish clean source commit: {}".format(exc))
    try:
        first = outputs[0].decode("ascii").strip()
        status = outputs[1].decode("utf-8")
        untracked = outputs[2].decode("utf-8")
        second = outputs[3].decode("ascii").strip()
    except UnicodeError as exc:
        fail("git source identity is not valid text: {}".format(exc))
    if re.fullmatch(r"[0-9a-f]{40}", first) is None or first != second:
        fail("git HEAD is not a full lowercase 40-hex commit")
    if status or untracked:
        paths = [line for line in (status + untracked).splitlines() if line]
        fail("codec source tree is dirty or has untracked source: {}".format(
            "; ".join(paths[:8])))
    return first


def _host_identity(controller_cpu: int) -> Mapping[str, Any]:
    return {
        "controller_cpu": controller_cpu,
        "hostname": socket.gethostname(),
        "kernel_release": platform.release(),
        "machine": platform.machine(),
        "system": platform.system(),
    }


def _emit_and_assemble_trace(
        contract: Mapping[str, Any], kind: str,
        description: Mapping[str, Any], output_dir: Path,
        deadline: float) -> Tuple[Path, str]:
    if kind != "recovery":
        fail("timing traces must come from the selected qualification audit")
    worker_path = description["resolved_path"]
    native_path = output_dir / "{}-native-traces.jsonl".format(kind)
    trace_path = output_dir / "{}-traces.jsonl".format(kind)
    data = _run_command(
        [worker_path, "--emit-traces", kind], deadline,
        "{} native trace emission".format(kind))
    _atomic_write_bytes(native_path, data)
    try:
        digest = native_api.assemble_trace_manifest(
            contract, kind, "development", native_path, trace_path)
    except (native_api.NativeEvidenceError,
            contract_api.ContractError) as exc:
        fail(str(exc))
    return trace_path, digest


def _timing_qualification_controls(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        ) -> List[Mapping[str, Any]]:
    described = {
        value["arm"]: value for value in description["arms"]
    }
    controls = []
    for policy in contract["timing"]["qualification"]["controls"]:
        arm = policy["arm"]
        value = described.get(arm)
        if value is None:
            fail("native worker omits a mandatory qualification control")
        controls.append({
            "arm": arm,
            "scope": policy["scope"],
            "binary_sha256": description["binary_sha256"],
            "arm_descriptor_sha256": value["arm_descriptor_sha256"],
            "construction_policy":
                "not_applicable" if arm == "wirehair1" else "raw_base",
            "repair_map_sha256": ZERO_SHA256,
        })
    try:
        contract_api._validate_timing_qualification_controls(
            contract, "development", controls)
    except contract_api.ContractError as exc:
        fail(str(exc))
    return controls


def _freeze_manifest(
        contract: Mapping[str, Any], kind: str, trace_sha256: str,
        source_commit: str, description: Mapping[str, Any],
        cpus: Sequence[int], controller_cpu: int,
        timing_qualification: Optional[
            contract_api.TimingQualification] = None) -> Mapping[str, Any]:
    roster = [value[0] for value in EXPECTED_ARMS]
    arms = []
    for value in description["arms"]:
        arms.append({
            "arm": value["arm"],
            "codec": value["codec"],
            "binary_sha256": description["binary_sha256"],
            "arm_descriptor_sha256": value["arm_descriptor_sha256"],
            "construction_policy": "not_applicable"
                if value["codec"] == "wirehair1" else "raw_base",
            "repair_map_sha256": ZERO_SHA256,
        })
    commands = ([] if kind == "timing" else [[
        description["resolved_path"], "--emit-traces", kind,
    ]]) + [
        [description["resolved_path"], "--worker", str(cpu)] for cpu in cpus
    ]
    if kind == "timing":
        if timing_qualification is None:
            fail("timing freeze requires a validated qualification")
        domain_sha256 = timing_qualification.qualified_domain_sha256
    else:
        domain_sha256 = contract[kind]["domains"]["development"][
            "domain_sha256"]
    return {
        "schema": contract_api.FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": kind,
        "phase": "development",
        "domain_sha256": domain_sha256,
        "source_git_commit": source_commit,
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": commands,
        "cpu_affinity": list(cpus),
        "host_identity": _host_identity(controller_cpu),
        "arms": arms,
    }


def write_development_timing_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        cpus: Sequence[int], controller_cpu: int, source_commit: str,
        trace_sha256: str, output_dir: Path,
        timing_qualification: contract_api.TimingQualification,
        ) -> Mapping[str, Any]:
    """Publish and strictly re-open the sole development timing freeze."""
    path = output_dir / "timing-freeze.json"
    value = _freeze_manifest(
        contract, "timing", trace_sha256, source_commit,
        description, cpus, controller_cpu, timing_qualification)
    _atomic_write_object(path, value)
    try:
        return contract_api.load_freeze_manifest(
            contract, "development", path, "timing", timing_qualification)
    except contract_api.ContractError as exc:
        fail(str(exc))


def _valid_sampler_monotonic_ns(row: Sequence[str]) -> Optional[int]:
    if len(row) != len(native_api.THERMAL_HEADER):
        return None
    try:
        numeric = [float(row[index]) for index in (
            1, 2, 3, 4, 14, 15, 16)]
        dimms = [float(value) for value in row[5:13]]
    except (TypeError, ValueError):
        return None
    if any(not math.isfinite(value) for value in numeric + dimms):
        return None
    canonical_counters = []
    for index in (13, 17, 18):
        value = row[index]
        if not value or not value.isascii() or not value.isdecimal() or \
                (len(value) > 1 and value.startswith("0")):
            return None
        canonical_counters.append(int(value))
    mono, busy, mhz, cpu_c, load1, load5, load15 = numeric
    if (mono < 0.0 or not 0.0 <= busy <= 100.0 or mhz <= 0.0 or
            not 0.0 <= cpu_c <= 120.0 or
            any(not 0.0 < value <= 100.0 for value in dimms) or
            any(value < 0.0 for value in (load1, load5, load15)) or
            any(value != 0 for value in canonical_counters)):
        return None
    return int(round(mono * 1000000000.0))


def _valid_sampler_samples(csv_path: Path) -> List[int]:
    values = []
    try:
        data = csv_path.read_bytes()
        # The live sampler may be in the middle of appending its last row.
        # Never select bytes which the sealed attestation correctly rejects as
        # an unterminated endpoint.
        complete_end = data.rfind(b"\n") + 1
        complete = data[:complete_end]
        text = complete.decode("ascii")
        with io.StringIO(text, newline="") as source:
            reader = csv.reader(source)
            if next(reader, None) != list(native_api.THERMAL_HEADER):
                fail("sampler CSV header differs from the frozen format")
            for row in reader:
                value = _valid_sampler_monotonic_ns(row)
                if value is not None:
                    values.append(value)
    except (OSError, UnicodeError, csv.Error) as exc:
        fail("cannot read sampler CSV: {}".format(exc))
    return values


def _wait_for_sampler_sample(
        csv_path: Path, deadline: float, *, after_ns: Optional[int] = None,
        at_or_after_ns: Optional[int] = None) -> int:
    if (after_ns is None) == (at_or_after_ns is None):
        fail("sampler wait requires exactly one time bound")
    while True:
        samples = _valid_sampler_samples(csv_path)
        if after_ns is not None:
            matches = [value for value in samples if value > after_ns]
        else:
            matches = [value for value in samples
                       if value >= int(at_or_after_ns)]
        if matches:
            return min(matches)
        remaining = _remaining(deadline, "waiting for thermal sample")
        time.sleep(min(0.2, remaining))


def choose_new_sampler_start(csv_path: Path, deadline: float) -> int:
    samples = _valid_sampler_samples(csv_path)
    if not samples:
        fail("sampler CSV has no valid baseline sample")
    return _wait_for_sampler_sample(
        csv_path, deadline, after_ns=max(samples))


def _preflight_sampler(pid: int, cpu: int, script_path: Path,
                       csv_path: Path) -> int:
    descriptors = []
    try:
        for path, context in ((script_path, "sampler script"),
                              (csv_path, "sampler CSV")):
            descriptor, _ = native_api._open_regular_nofollow(path, context)
            descriptors.append(descriptor)
        start_ticks = native_api._parse_proc_start_ticks(pid)
        native_api._verify_live_sampler_process(
            pid, cpu, start_ticks, script_path, csv_path)
        return start_ticks
    except native_api.NativeEvidenceError as exc:
        root_owned = False
        try:
            root_owned = Path("/proc/{}".format(pid)).stat().st_uid == 0
        except OSError:
            pass
        if root_owned and getattr(os, "geteuid", lambda: 1)() != 0:
            fail("{}; sampler PID {} is root-owned and its executable "
                 "provenance requires privilege. Rerun this controller with "
                 "sudo -n python3.12; validation is not weakened".format(
                     exc, pid))
        fail(str(exc))
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
    return 0


class Job:
    def __init__(self, kind: str, cell: int, item: int, ordinal: int,
                 retry_offset: int = 0) -> None:
        if (kind not in ("recovery", "timing", "qualification") or
                any(type(value) is not int or value < 0
                    for value in (cell, item, ordinal, retry_offset)) or
                retry_offset >= 256):
            fail("cannot construct an invalid native job")
        if kind == "recovery" and retry_offset != 0:
            fail("recovery jobs cannot carry a timing retry")
        if kind == "qualification" and (
                item != 0 or ordinal != cell * 256 + retry_offset):
            fail("qualification job identity is inconsistent")
        self.kind = kind
        self.cell = cell
        self.item = item
        self.ordinal = ordinal
        self.retry_offset = retry_offset

    def command(self) -> bytes:
        if self.kind == "recovery":
            prefix = "R"
            item = self.item
        elif self.kind == "qualification":
            prefix = "L"
            item = self.retry_offset
        else:
            prefix = "T"
            item = self.item * 256 + self.retry_offset
        return "{} {} {}\n".format(prefix, self.cell, item).encode("ascii")


class PersistentWorker:
    def __init__(self, cpu: int, process: subprocess.Popen,
                 start_ticks: int) -> None:
        self.cpu = cpu
        self.process = process
        self.start_ticks = start_ticks
        self.pending: Optional[Job] = None
        self.buffer = b""


class _ProvisionalWorker(PersistentWorker):
    """Own a just-spawned child until its public worker handoff completes."""

    def __init__(self, cpu: int) -> None:
        self.cpu = cpu
        self.process: Any = None
        self.start_ticks = 0
        self.pending = None
        self.buffer = b""


def _worker_stderr(worker: PersistentWorker) -> str:
    descriptor = worker.process.stderr.fileno()
    try:
        data = os.read(descriptor, 65536)
    except (BlockingIOError, OSError):
        data = b""
    return data.decode("utf-8", "replace").strip()


def spawn_workers(description: Mapping[str, Any], cpus: Sequence[int],
                  deadline: float) -> List[PersistentWorker]:
    workers: List[PersistentWorker] = []
    provisional_worker: Optional[_ProvisionalWorker] = None
    try:
        for cpu in cpus:
            # Allocate the owner before Popen.  Immediately after Popen returns,
            # its process is stored here so constructor/list-append interrupts
            # cannot leave a live child outside the cleanup pool.  Only the
            # unavoidable syscall CALL-to-Python-STORE sliver precedes it.
            provisional_worker = _ProvisionalWorker(cpu)
            try:
                provisional_worker.process = subprocess.Popen(
                    [description["resolved_path"], "--worker", str(cpu)],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, bufsize=0)
            except OSError as exc:
                fail("cannot spawn native worker on CPU {}: {}".format(cpu, exc))
            workers.append(PersistentWorker(
                cpu, provisional_worker.process, 0))
            provisional_worker = None
        pending = set(range(len(workers)))
        ready_deadline = min(deadline, time.monotonic() + 10.0)
        while pending:
            if time.monotonic() >= ready_deadline:
                fail("native workers did not establish singleton affinity")
            for index in list(pending):
                worker = workers[index]
                if worker.process.poll() is not None:
                    fail("native worker on CPU {} exited during startup: {}".format(
                        worker.cpu, _worker_stderr(worker)))
                try:
                    ticks = native_api._parse_proc_start_ticks(
                        worker.process.pid)
                    affinity = os.sched_getaffinity(worker.process.pid)
                except (native_api.NativeEvidenceError,
                        AttributeError, OSError):
                    continue
                if affinity == {worker.cpu}:
                    worker.start_ticks = ticks
                    os.set_blocking(worker.process.stdout.fileno(), False)
                    os.set_blocking(worker.process.stderr.fileno(), False)
                    pending.remove(index)
            if pending:
                time.sleep(0.01)
        return workers
    except BaseException as primary:
        cleanup_workers = list(workers)
        if (provisional_worker is not None and
                provisional_worker.process is not None and
                all(worker.process is not provisional_worker.process
                    for worker in cleanup_workers)):
            cleanup_workers.append(provisional_worker)
        _finish_worker_pool_cleanup(
            cleanup_workers, False, "partial worker cleanup failed", primary)
        raise


def terminate_workers(workers: Sequence[PersistentWorker]) -> None:
    """Exhaust terminate/kill/reap/close and reject every surviving child."""
    failures: List[CleanupFailure] = []

    def alive(worker: PersistentWorker, stage: str) -> bool:
        try:
            return worker.process.poll() is None
        except BaseException as cleanup:
            failures.append((
                "native worker CPU {} {} poll failed".format(
                    worker.cpu, stage), cleanup))
            # A process whose state cannot be observed must receive every
            # remaining termination/reap attempt.
            return True

    for worker in workers:
        if alive(worker, "pre-terminate"):
            try:
                worker.process.terminate()
            except BaseException as cleanup:
                failures.append((
                    "native worker CPU {} terminate failed".format(
                        worker.cpu), cleanup))
    terminate_limit = time.monotonic() + 2.0
    for worker in workers:
        if alive(worker, "pre-wait"):
            must_kill = False
            try:
                worker.process.wait(timeout=max(
                    0.0, terminate_limit - time.monotonic()))
            except subprocess.TimeoutExpired:
                must_kill = True
            except BaseException as cleanup:
                failures.append((
                    "native worker CPU {} graceful wait failed".format(
                        worker.cpu), cleanup))
                must_kill = True
            if must_kill:
                try:
                    worker.process.kill()
                except BaseException as cleanup:
                    failures.append((
                        "native worker CPU {} kill failed".format(
                            worker.cpu), cleanup))
    reap_limit = time.monotonic() + 2.0
    for worker in workers:
        if alive(worker, "pre-terminal-wait"):
            try:
                worker.process.wait(timeout=max(
                    0.0, reap_limit - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
            except BaseException as cleanup:
                failures.append((
                    "native worker CPU {} terminal wait failed".format(
                        worker.cpu), cleanup))
        for stream in (worker.process.stdin, worker.process.stdout,
                       worker.process.stderr):
            if stream is not None:
                try:
                    stream.close()
                except BaseException as cleanup:
                    failures.append((
                        "native worker CPU {} stream close failed".format(
                            worker.cpu), cleanup))
    survivors = [worker.cpu for worker in workers
                 if alive(worker, "terminal")]
    terminal_diagnostics = []
    if survivors:
        terminal_diagnostics.append(
            "native worker cleanup left unreaped CPUs {}".format(
                ",".join(str(cpu) for cpu in survivors)))
    _raise_after_cleanup(None, failures, terminal_diagnostics)


def _finish_worker_pool_cleanup(
        workers: Sequence[PersistentWorker], clean_shutdown: bool,
        context: str, primary: Optional[BaseException]) -> None:
    """Finish one worker pool and combine cleanup with its active failure."""
    failures: List[CleanupFailure] = []
    if workers and not clean_shutdown:
        try:
            terminate_workers(workers)
        except BaseException as cleanup:
            failures.append((context, cleanup))
    _raise_after_cleanup(primary, failures)


def _write_worker(worker: PersistentWorker, data: bytes) -> None:
    if worker.process.poll() is not None:
        fail("native worker on CPU {} exited before its command".format(
            worker.cpu))
    try:
        written = os.write(worker.process.stdin.fileno(), data)
    except (BrokenPipeError, OSError) as exc:
        fail("cannot command native worker on CPU {}: {}".format(
            worker.cpu, exc))
    if written != len(data):
        fail("short write to native worker on CPU {}".format(worker.cpu))


ResponseValidator = Callable[
    [Mapping[str, Any], bytes, PersistentWorker, Job], int]


def run_job_batch(
        workers: Sequence[PersistentWorker], jobs: Sequence[Job],
        rotation: int, sink: AtomicLineSink, deadline: float,
        validator: ResponseValidator,
        maximum_response_bytes: int = MAX_RESPONSE_BYTES,
        ) -> Tuple[int, Set[int]]:
    """Run one barrier-delimited batch with at most one job per worker."""
    if (type(maximum_response_bytes) is not int or
            not 0 < maximum_response_bytes <= MAX_RESPONSE_BYTES):
        fail("native response byte cap must be in [1,{}]".format(
            MAX_RESPONSE_BYTES))
    if not workers or len(jobs) < len(workers):
        fail("each batch must have enough jobs to use every frozen worker")
    if any(worker.pending is not None or worker.buffer for worker in workers):
        fail("worker state leaked across a job barrier")
    selector = selectors.DefaultSelector()
    next_job = 0
    completed = 0
    maximum_end = 0
    used: Set[int] = set()
    try:
        for worker in workers:
            selector.register(
                worker.process.stdout, selectors.EVENT_READ, worker)
        for offset in range(len(workers)):
            worker = workers[(offset + rotation) % len(workers)]
            worker.pending = jobs[next_job]
            next_job += 1
            _write_worker(worker, worker.pending.command())
        while completed < len(jobs):
            timeout = min(1.0, _remaining(deadline, "native result jobs"))
            events = selector.select(timeout)
            if not events:
                for worker in workers:
                    if worker.process.poll() is not None:
                        fail("native worker on CPU {} exited: {}".format(
                            worker.cpu, _worker_stderr(worker)))
                continue
            for key, _ in events:
                worker = key.data
                try:
                    block = os.read(worker.process.stdout.fileno(), 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    fail("cannot read native worker on CPU {}: {}".format(
                        worker.cpu, exc))
                if not block:
                    fail("native worker on CPU {} reached EOF: {}".format(
                        worker.cpu, _worker_stderr(worker)))
                worker.buffer += block
                if len(worker.buffer) > maximum_response_bytes:
                    fail("native worker response exceeds its exact byte cap")
                if b"\n" not in worker.buffer:
                    continue
                logical, remainder = worker.buffer.split(b"\n", 1)
                if remainder:
                    fail("native worker emitted more than one response per command")
                line = logical + b"\n"
                worker.buffer = b""
                job = worker.pending
                if job is None:
                    fail("native worker responded without an outstanding job")
                value = _parse_canonical_line(line, "native worker response")
                finished_ns = validator(value, line, worker, job)
                if type(finished_ns) is not int or finished_ns < 0:
                    fail("native response validator returned an invalid end time")
                maximum_end = max(maximum_end, finished_ns)
                sink.write(line)
                used.add(worker.cpu)
                completed += 1
                worker.pending = None
                if next_job < len(jobs):
                    worker.pending = jobs[next_job]
                    next_job += 1
                    _write_worker(worker, worker.pending.command())
        if any(worker.pending is not None or worker.buffer for worker in workers):
            fail("worker state is not empty at the job barrier")
        return maximum_end, used
    finally:
        selector.close()


def _strict_qualification_response(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        base_cells: Sequence[Mapping[str, Any]], window_start_ns: int,
        value: Mapping[str, Any], worker: PersistentWorker, job: Job,
        ) -> Tuple[int, bool]:
    payload = value.get("payload")
    integer_fields = (
        "ordinal", "cpu", "worker_pid", "worker_process_start_ticks",
        "started_monotonic_ns", "finished_monotonic_ns",
    )
    if (job.kind != "qualification" or
            set(value) != native_api.NATIVE_RECORD_FIELDS or
            any(type(value.get(field)) is not int for field in integer_fields) or
            value.get("schema") !=
                native_api.TIMING_QUALIFICATION_RECORD_SCHEMA or
            value.get("ordinal") != job.ordinal or
            value.get("cpu") != worker.cpu or
            value.get("worker_pid") != worker.process.pid or
            value.get("worker_process_start_ticks") != worker.start_ticks or
            value.get("worker_binary_sha256") !=
                description["binary_sha256"] or
            not isinstance(payload, dict) or
            set(payload) != native_api.TIMING_QUALIFICATION_PAYLOAD_FIELDS):
        fail("native qualification response does not bind its worker/job")
    start = value["started_monotonic_ns"]
    end = value["finished_monotonic_ns"]
    if not window_start_ns <= start <= end:
        fail("native qualification response has invalid monotonic provenance")
    for field in native_api.SHA256_FIELDS:
        if not _is_sha256(value.get(field)):
            fail("native qualification response {} is not a SHA-256".format(
                field))
    if not 0 <= job.cell < len(base_cells):
        fail("native qualification job is outside the base domain")
    base_cell = base_cells[job.cell]
    base_sha256 = contract_api.sha256_json(base_cell)
    loss_seed = contract_api._qualified_timing_loss_seed(
        base_cell["base_loss_seed"], job.retry_offset)
    qualified_cell = dict(base_cell)
    qualified_cell["base_cell_sha256"] = base_sha256
    qualified_cell["loss_retry_offset"] = job.retry_offset
    qualified_cell["loss_seed"] = loss_seed
    cell_sha256 = contract_api.sha256_json(qualified_cell)
    if (payload.get("ordinal") != job.cell or
            payload.get("loss_retry_offset") != job.retry_offset or
            payload.get("base_cell_sha256") != base_sha256 or
            payload.get("loss_seed") != loss_seed or
            payload.get("cell_sha256") != cell_sha256 or
            value["message_sha256"] !=
                native_api._timing_qualification_message_sha256(base_cell) or
            value["work_sha256"] !=
                native_api._expected_timing_qualification_work_sha256(
                    job.ordinal, cell_sha256)):
        fail("native qualification response substitutes its base/retry identity")
    packet_count = payload.get("packet_count")
    candidate_count = payload.get("candidate_count")
    expected_packets = base_cell["K"] + \
        contract["timing"]["receive_overhead_cap"]
    if (not _is_sha256(payload.get("trace_sha256")) or
            type(packet_count) is not int or
            packet_count != expected_packets or
            type(candidate_count) is not int or
            not packet_count <= candidate_count <=
                256 * expected_packets + 65536):
        fail("native qualification response has invalid trace evidence")
    try:
        head = contract_api._validate_timing_qualification_outcome(
            payload, "wirehair2_head",
            contract["timing"]["receive_overhead_cap"])
        wh1 = contract_api._validate_timing_qualification_outcome(
            payload, "wirehair1",
            contract["timing"]["receive_overhead_cap"])
    except contract_api.ContractError as exc:
        fail(str(exc))
    return end, head == "success" and wh1 == "success"


def run_timing_qualification(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        workers: Sequence[PersistentWorker], output_path: Path,
        window_start_ns: int, deadline: float) -> Tuple[Path, int]:
    """Run lowest, non-speculative retries while keeping base cells parallel."""
    base_cells = list(contract_api.iter_timing_base_cells(
        contract, "development"))
    expected_cells = contract["timing"]["domains"]["development"][
        "expected_cells"]
    if (not workers or len(base_cells) != expected_cells or
            len(workers) > expected_cells):
        fail("qualification worker/base-cell cardinality is invalid")
    if any(worker.pending is not None or worker.buffer for worker in workers):
        fail("worker state leaked into timing qualification")
    sink: Optional[AtomicLineSink] = None
    selector: Optional[selectors.BaseSelector] = None
    next_base = 0
    selected = 0
    maximum_end = 0
    used: Set[int] = set()

    def assign(worker: PersistentWorker, base_ordinal: int,
               retry_offset: int) -> None:
        worker.pending = Job(
            "qualification", base_ordinal, 0,
            base_ordinal * 256 + retry_offset, retry_offset)
        _write_worker(worker, worker.pending.command())

    try:
        sink = AtomicLineSink(
            output_path,
            COMPLETED_ARTIFACT_BYTE_LIMITS[
                "timing-qualification-native.jsonl"])
        selector = selectors.DefaultSelector()
        for worker in workers:
            selector.register(
                worker.process.stdout, selectors.EVENT_READ, worker)
        for worker in workers:
            assign(worker, next_base, 0)
            next_base += 1
        while selected < expected_cells:
            timeout = min(1.0, _remaining(
                deadline, "native timing qualification"))
            events = selector.select(timeout)
            if not events:
                for worker in workers:
                    if worker.process.poll() is not None:
                        fail("qualification worker on CPU {} exited: {}".format(
                            worker.cpu, _worker_stderr(worker)))
                continue
            for key, _ in events:
                worker = key.data
                try:
                    block = os.read(worker.process.stdout.fileno(), 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    fail("cannot read qualification worker on CPU {}: {}".
                         format(worker.cpu, exc))
                if not block:
                    fail("qualification worker on CPU {} reached EOF: {}".
                         format(worker.cpu, _worker_stderr(worker)))
                worker.buffer += block
                if len(worker.buffer) > MAX_RESPONSE_BYTES:
                    fail("qualification response exceeds the bounded line size")
                if b"\n" not in worker.buffer:
                    continue
                logical, remainder = worker.buffer.split(b"\n", 1)
                if remainder:
                    fail("qualification worker emitted multiple responses")
                line = logical + b"\n"
                if len(line) > MAX_QUALIFICATION_NATIVE_RECORD_BYTES:
                    fail("qualification response exceeds its canonical row cap")
                worker.buffer = b""
                job = worker.pending
                if job is None:
                    fail("qualification worker responded without a job")
                value = _parse_canonical_line(
                    line, "native qualification response")
                end, both_success = _strict_qualification_response(
                    contract, description, base_cells, window_start_ns,
                    value, worker, job)
                maximum_end = max(maximum_end, end)
                sink.write(line)
                used.add(worker.cpu)
                worker.pending = None
                if both_success:
                    selected += 1
                    if next_base < expected_cells:
                        assign(worker, next_base, 0)
                        next_base += 1
                else:
                    if job.retry_offset == 255:
                        fail("timing qualification exhausted retry offset 255")
                    # No later attempt for this base is issued until this
                    # exact retryable result has arrived and validated.
                    assign(worker, job.cell, job.retry_offset + 1)
        if next_base != expected_cells or \
                any(worker.pending is not None or worker.buffer
                    for worker in workers):
            fail("timing qualification ended with incomplete worker state")
        if used != {worker.cpu for worker in workers}:
            fail("timing qualification did not use every allowed logical CPU")
        sink.publish()
        return output_path, maximum_end
    finally:
        try:
            if selector is not None:
                selector.close()
        finally:
            if sink is not None:
                sink.abort()


def quit_workers(workers: Sequence[PersistentWorker], deadline: float) -> None:
    for worker in workers:
        if worker.pending is not None or worker.buffer:
            fail("cannot quit a worker with unconsumed result state")
        _write_worker(worker, b"Q\n")
        worker.process.stdin.close()
    for worker in workers:
        try:
            returncode = worker.process.wait(
                timeout=_remaining(deadline, "native worker shutdown"))
        except subprocess.TimeoutExpired:
            fail("native worker on CPU {} did not exit after Q".format(
                worker.cpu))
        stdout = b""
        try:
            while True:
                block = os.read(worker.process.stdout.fileno(), 65536)
                if not block:
                    break
                stdout += block
        except BlockingIOError:
            pass
        stderr = _worker_stderr(worker)
        if returncode != 0 or stdout or stderr:
            fail("native worker on CPU {} failed exact Q shutdown{}".format(
                worker.cpu, ": " + stderr if stderr else ""))
        worker.process.stdout.close()
        worker.process.stderr.close()


def _strict_response_validator(
        contract: Mapping[str, Any], freeze: Mapping[str, Any], kind: str,
        description: Mapping[str, Any], window_start_ns: int,
        timing_qualification: Optional[
            contract_api.TimingQualification] = None,
        ) -> ResponseValidator:
    try:
        ordinal_context = native_api._record_ordinal_context(
            contract, freeze, kind, "development", timing_qualification)
    except (native_api.NativeEvidenceError,
            contract_api.ContractError) as exc:
        fail(str(exc))

    def validate(value: Mapping[str, Any], _line: bytes,
                 worker: PersistentWorker, job: Job) -> int:
        if kind == "timing" and len(_line) > MAX_TIMING_NATIVE_RECORD_BYTES:
            fail("native timing response exceeds its canonical row cap")
        payload = value.get("payload")
        if kind == "recovery":
            if not isinstance(payload, dict):
                fail("native response payload has an unexpected schema")
            try:
                expected_schema, payload_fields = \
                    native_api.recovery_record_schema_fields(freeze, payload)
            except native_api.NativeEvidenceError as exc:
                fail(str(exc))
        else:
            expected_schema = native_api.TIMING_RECORD_SCHEMA
            payload_fields = contract_api.TIMING_RECEIPT_FIELDS
        integer_fields = (
            "ordinal", "cpu", "worker_pid", "worker_process_start_ticks",
            "started_monotonic_ns", "finished_monotonic_ns",
        )
        if job.kind != kind or \
                set(value) != native_api.NATIVE_RECORD_FIELDS or \
                any(type(value.get(field)) is not int
                    for field in integer_fields) or \
                value.get("schema") != expected_schema or \
                value.get("ordinal") != job.ordinal or \
                value.get("cpu") != worker.cpu or \
                value.get("worker_pid") != worker.process.pid or \
                value.get("worker_process_start_ticks") != worker.start_ticks or \
                value.get("worker_binary_sha256") != \
                description["binary_sha256"]:
            fail("native response does not bind its commanded worker/job")
        start = value.get("started_monotonic_ns")
        end = value.get("finished_monotonic_ns")
        if type(start) is not int or type(end) is not int or \
                not window_start_ns <= start <= end:
            fail("native response has invalid monotonic provenance")
        for field in native_api.SHA256_FIELDS:
            if not _is_sha256(value.get(field)):
                fail("native response {} is not a SHA-256".format(field))
        if not isinstance(payload, dict) or set(payload) != payload_fields:
            fail("native response payload has an unexpected schema")
        try:
            ordinal, _ = native_api._expected_record_ordinal(
                contract, freeze, kind, "development", payload,
                ordinal_context)
            work_sha256 = native_api._expected_work_sha256(
                kind, "development", job.ordinal, payload)
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        if ordinal != job.ordinal or value["work_sha256"] != work_sha256:
            fail("native response payload/work hash differs from its command")
        if kind == "recovery":
            frozen_arm = freeze["arms_by_name"].get(payload.get("arm"))
            if (frozen_arm is None or
                    payload.get("binary_sha256") != description["binary_sha256"] or
                    payload.get("arm_descriptor_sha256") !=
                    frozen_arm["arm_descriptor_sha256"]):
                fail("native recovery response differs from its frozen arm")
        else:
            for side in ("left", "right"):
                frozen_arm = freeze["arms_by_name"].get(
                    payload.get(side + "_arm"))
                if (frozen_arm is None or
                        payload.get(side + "_binary_sha256") !=
                        description["binary_sha256"] or
                        payload.get(side + "_arm_descriptor_sha256") !=
                        frozen_arm["arm_descriptor_sha256"]):
                    fail("native timing response differs from its frozen arm")
        return end

    return validate


def _recovery_jobs(contract: Mapping[str, Any], arm_count: int) -> List[Job]:
    """Shared recovery helper retained for the dedicated recovery runner/tests."""
    cells = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    jobs = [
        Job("recovery", cell, arm, cell * arm_count + arm)
        for cell in range(cells) for arm in range(arm_count)
    ]
    if (cells != 360 or arm_count != len(EXPECTED_ARMS) or
            len(jobs) != cells * arm_count):
        fail("development recovery jobs differ from the frozen arm roster")
    return jobs


def _development_timing_shape(
        contract: Mapping[str, Any], panel_count: int,
        ) -> Tuple[int, int, int]:
    timing = contract.get("timing")
    if not isinstance(timing, Mapping):
        fail("development timing policy is missing")
    geometry = timing.get("execution_geometry")
    if geometry != contract_api.EXPECTED_TIMING_EXECUTION_GEOMETRY:
        fail("development timing execution geometry is not frozen")
    if geometry.get("timing_worker_count") != TIMING_WORKER_COUNT or \
            geometry.get("jobs_per_wave") != TIMING_WORKER_COUNT:
        fail("development timing geometry does not use exactly eight workers")
    domains = timing.get("domains")
    domain = domains.get("development") \
        if isinstance(domains, Mapping) else None
    repetitions = domain.get("paired_repetitions") \
        if isinstance(domain, Mapping) else None
    independent_rounds = domain.get("independent_rounds") \
        if isinstance(domain, Mapping) else None
    expected_cells = domain.get("expected_cells") \
        if isinstance(domain, Mapping) else None
    if (type(repetitions) is not int or repetitions <= 0 or
            repetitions % TIMING_WORKER_COUNT != 0):
        fail("development timing repetitions must be a positive multiple of eight")
    if (type(independent_rounds) is not int or independent_rounds <= 0 or
            independent_rounds * TIMING_WORKER_COUNT != repetitions):
        fail("development timing rounds must exactly partition the repetitions")
    if (type(expected_cells) is not int or
            expected_cells != FROZEN_TIMING_BASE_CELLS or
            expected_cells % repetitions != 0):
        fail("development timing cell cardinality differs from the "
             "completed-artifact bound")
    try:
        expected_panels = len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS]))
    except contract_api.ContractError as exc:
        fail(str(exc))
    if (expected_panels != FROZEN_TIMING_PANEL_COUNT or
            type(panel_count) is not int or
            panel_count != FROZEN_TIMING_PANEL_COUNT):
        fail("development timing panel cardinality differs from the "
             "completed-artifact bound")
    return repetitions, independent_rounds, expected_cells


def _timing_job_waves(
        contract: Mapping[str, Any], panel_count: int,
        timing_qualification: contract_api.TimingQualification,
        ) -> List[Tuple[int, List[Job]]]:
    """Build homogeneous, barrier-delimited eight-job timing waves."""
    repetitions, independent_rounds, expected_cells = _development_timing_shape(
        contract, panel_count)
    try:
        cells = list(contract_api.iter_timing_cells(
            contract, "development", timing_qualification))
    except contract_api.ContractError as exc:
        fail(str(exc))
    if len(cells) != expected_cells:
        fail("development timing cell cardinality differs from its contract")

    timing = contract["timing"]
    cell_fields = timing.get("cell_key")
    if (not isinstance(cell_fields, list) or
            any(not isinstance(field, str) for field in cell_fields) or
            len(cell_fields) != len(set(cell_fields)) or
            "replicate" not in cell_fields or "loss_seed" not in cell_fields):
        fail("development timing cell identity is invalid")
    stable_fields = [
        field for field in cell_fields
        if field not in (
            "replicate", "base_loss_seed", "base_cell_sha256",
            "loss_retry_offset", "loss_seed")
    ]
    by_replicate: List[Dict[str, Tuple[int, Mapping[str, Any]]]] = [
        {} for _ in range(repetitions)
    ]
    identity_order: List[str] = []
    expected_field_set = set(cell_fields)
    for cell_ordinal, cell in enumerate(cells):
        if not isinstance(cell, Mapping) or set(cell) != expected_field_set:
            fail("development timing cell has an unexpected identity")
        replicate = cell.get("replicate")
        if (type(replicate) is not int or
                not 0 <= replicate < repetitions):
            fail("development timing cell has an invalid replicate")
        identity = contract_api.canonical_json({
            field: cell[field] for field in stable_fields
        })
        if identity in by_replicate[replicate]:
            fail("duplicate timing cell identity within one replicate")
        by_replicate[replicate][identity] = (cell_ordinal, cell)
        if replicate == 0:
            identity_order.append(identity)

    cells_per_replicate = expected_cells // repetitions
    if (len(identity_order) != cells_per_replicate or
            len(set(identity_order)) != cells_per_replicate):
        fail("development timing replicate zero has invalid cardinality")
    expected_identities = set(identity_order)
    for replicate, values in enumerate(by_replicate):
        if (len(values) != cells_per_replicate or
                set(values) != expected_identities):
            fail("timing cohort identity differs across replicate {}".format(
                replicate))

    waves: List[Tuple[int, List[Job]]] = []
    frozen_arms = [value[0] for value in EXPECTED_ARMS]
    cohorts = [
        (cohort_index, identity, panel)
        for cohort_index, (identity, panel) in enumerate(
            (value, panel)
            for value in identity_order for panel in range(panel_count))
    ]
    # A statistical round is one complete traversal of every frozen cohort.
    # Keeping the round outermost prevents adjacent waves of one K/panel from
    # masquerading as independent observations of host state.
    for round_index in range(independent_rounds):
        first_replicate = round_index * TIMING_WORKER_COUNT
        for cohort_index, identity, panel in cohorts:
            jobs = []
            for replicate in range(
                    first_replicate,
                    first_replicate + TIMING_WORKER_COUNT):
                cell_ordinal, selected_cell = \
                    by_replicate[replicate][identity]
                jobs.append(Job(
                    "timing", cell_ordinal, panel,
                    cell_ordinal * panel_count + panel,
                    selected_cell["loss_retry_offset"]))
            # run_job_batch assigns job position r to worker slot
            # (r + rotation) modulo eight.  Derive that rotation from the
            # bounded contract helper and verify the complete wave before
            # allowing any native command to be issued.
            try:
                rotation = contract_api.timing_worker_slot(
                    contract, "development", frozen_arms, cohort_index,
                    first_replicate)
                for position, replicate in enumerate(range(
                        first_replicate,
                        first_replicate + TIMING_WORKER_COUNT)):
                    expected_slot = contract_api.timing_worker_slot(
                        contract, "development", frozen_arms,
                        cohort_index, replicate)
                    if expected_slot != (
                            position + rotation) % TIMING_WORKER_COUNT:
                        fail("timing worker-slot formula is not a rotation")
            except contract_api.ContractError as exc:
                fail(str(exc))
            waves.append((rotation, jobs))

    expected_cohorts = cells_per_replicate * panel_count
    try:
        frozen_cohorts = contract_api.timing_cohort_count(
            contract, "development", frozen_arms)
    except contract_api.ContractError as exc:
        fail(str(exc))
    expected_waves = expected_cohorts * independent_rounds
    expected_jobs = {
        (cell, panel, cell * panel_count + panel,
         cells[cell]["loss_retry_offset"])
        for cell in range(expected_cells)
        for panel in range(panel_count)
    }
    actual_jobs = [
        (job.cell, job.item, job.ordinal, job.retry_offset)
        for _, jobs in waves for job in jobs
    ]
    if (len(cohorts) != expected_cohorts or
            expected_cohorts != frozen_cohorts or
            len(waves) != expected_waves or
            any(len(jobs) != TIMING_WORKER_COUNT for _, jobs in waves) or
            len(actual_jobs) != len(expected_jobs) or
            len(set(actual_jobs)) != len(actual_jobs) or
            set(actual_jobs) != expected_jobs):
        fail("development timing waves do not cover the exact frozen domain")
    return waves


def _run_timing_jobs(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        timing_qualification: contract_api.TimingQualification,
        description: Mapping[str, Any], workers: Sequence[PersistentWorker],
        output_dir: Path, window_start_ns: int, deadline: float,
        sibling_guard: Optional[SiblingIsolationGuard] = None,
        ) -> Tuple[Path, int]:
    if len(workers) != TIMING_WORKER_COUNT:
        fail("native short screen requires exactly eight timing workers")
    try:
        panel_count = len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS]))
    except contract_api.ContractError as exc:
        fail(str(exc))
    timing_waves = _timing_job_waves(
        contract, panel_count, timing_qualification)
    path = output_dir / "timing-native-results.jsonl"
    timing_sink: Optional[AtomicLineSink] = None
    maximum_end = 0
    try:
        timing_sink = AtomicLineSink(
            path,
            COMPLETED_ARTIFACT_BYTE_LIMITS["timing-native-results.jsonl"])
        expected_timing_cpus = set(worker.cpu for worker in workers)
        timing_validator = _strict_response_validator(
            contract, freeze, "timing", description,
            window_start_ns, timing_qualification)
        for wave_ordinal, (rotation, jobs) in enumerate(timing_waves):
            if sibling_guard is not None:
                sibling_guard.check(
                    "before-timing-wave-{}".format(wave_ordinal))
            end, used = run_job_batch(
                workers, jobs, rotation, timing_sink,
                deadline, timing_validator,
                MAX_TIMING_NATIVE_RECORD_BYTES)
            if sibling_guard is not None:
                sibling_guard.check(
                    "after-timing-wave-{}".format(wave_ordinal))
            if used != expected_timing_cpus:
                fail("timing wave {} did not use every frozen CPU".format(
                    wave_ordinal))
            maximum_end = max(maximum_end, end)
        timing_sink.publish()
        return path, maximum_end
    finally:
        if timing_sink is not None:
            timing_sink.abort()


def _create_output_dir(path: Path) -> Path:
    try:
        path.mkdir(mode=0o700, parents=True, exist_ok=False)
        resolved = path.resolve(strict=True)
        os.chmod(str(resolved), 0o700)
        info = resolved.stat()
        if (not stat.S_ISDIR(info.st_mode) or
                stat.S_IMODE(info.st_mode) != 0o700):
            fail("output path is not a directory")
        _fsync_directory(resolved.parent)
        return resolved
    except FileExistsError:
        fail("output directory already exists: {}".format(path))
    except OSError as exc:
        fail("cannot create output directory {}: {}".format(path, exc))
    return path


CompletedFingerprint = Tuple[int, int, int, int, int, int, int]


def _completed_fingerprint(info: os.stat_result) -> CompletedFingerprint:
    return (
        info.st_dev, info.st_ino, info.st_mode, info.st_nlink, info.st_size,
        info.st_mtime_ns, info.st_ctime_ns,
    )


def _read_completed_descriptor_bytes(
        descriptor: int, context: str, byte_limit: int,
        ) -> Tuple[bytes, CompletedFingerprint]:
    """Stably pread one retained bounded regular-file descriptor."""
    if type(byte_limit) is not int or byte_limit <= 0:
        fail("{} has an invalid completed-artifact byte cap".format(context))
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            fail("{} must be a single-link regular non-symlink file".format(
                context))
        if before.st_size > byte_limit:
            fail("{} exceeds the bounded artifact size".format(context))
        chunks = []
        size = 0
        while True:
            block = os.pread(
                descriptor, min(1024 * 1024, byte_limit + 1 - size), size)
            if not block:
                break
            size += len(block)
            if size > byte_limit:
                fail("{} exceeds the bounded artifact size".format(context))
            chunks.append(block)
        after = os.fstat(descriptor)
    except OSError as exc:
        fail("cannot read {}: {}".format(context, exc))
    stable_before = _completed_fingerprint(before)
    stable_after = _completed_fingerprint(after)
    data = b"".join(chunks)
    if (stable_before != stable_after or len(data) != before.st_size):
        fail("{} changed while it was being read".format(context))
    return data, stable_after


def _open_completed_regular_snapshot(
        path: Path, context: str, directory_fd: int, byte_limit: int,
        descriptor_owner: Dict[str, int], descriptor_key: str,
        ) -> Tuple[bytes, CompletedFingerprint]:
    """Open, snapshot, and retain one named dependency inode."""
    if (not isinstance(descriptor_owner, dict) or
            not isinstance(descriptor_key, str) or not descriptor_key or
            descriptor_key in descriptor_owner):
        fail("completed dependency descriptor owner is invalid")
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(
                context))
        descriptor = os.open(
            str(path), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | getattr(os, "O_NONBLOCK", 0), dir_fd=directory_fd)
        descriptor_owner[descriptor_key] = descriptor
        data, fingerprint = _read_completed_descriptor_bytes(
            descriptor, context, byte_limit)
        named = os.stat(
            str(path), dir_fd=directory_fd, follow_symlinks=False)
        if (not stat.S_ISREG(named.st_mode) or
                _completed_fingerprint(named) != fingerprint):
            fail("{} changed while it was being read".format(context))
        return data, fingerprint
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    finally:
        if descriptor >= 0 and sys.exc_info()[0] is not None:
            primary = sys.exc_info()[1]
            descriptor_owner.pop(descriptor_key, None)
            _drain_descriptor_closes((
                ("completed dependency descriptor cleanup failed",
                 descriptor),
            ), primary)
    raise AssertionError("unreachable")


def _read_completed_regular_bytes(
        path: Path, context: str, directory_fd: int,
        byte_limit: int) -> Tuple[bytes, CompletedFingerprint]:
    """Read one bounded regular file through an already-pinned directory."""
    descriptor_owner: Dict[str, int] = {}
    try:
        data, fingerprint = _open_completed_regular_snapshot(
            path, context, directory_fd, byte_limit,
            descriptor_owner, "artifact")
        return data, fingerprint
    finally:
        descriptor = descriptor_owner.pop("artifact", -1)
        if descriptor >= 0:
            _drain_descriptor_closes((
                ("completed artifact descriptor cleanup failed", descriptor),
            ), sys.exc_info()[1])


def _read_completed_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        ) -> Tuple[Dict[str, bytes], Dict[str, CompletedFingerprint]]:
    snapshots: Dict[str, bytes] = {}
    identities: Dict[str, CompletedFingerprint] = {}
    total_bytes = 0
    for key, name in source_names.items():
        limit = COMPLETED_ARTIFACT_BYTE_LIMITS.get(name)
        if type(limit) is not int or limit <= 0:
            fail("completed timing bundle names an unbounded artifact")
        data, identity = _read_completed_regular_bytes(
            Path(name), "completed timing " + key.replace("_", " "),
            directory_fd, limit)
        total_bytes += len(data)
        if total_bytes > MAX_COMPLETED_BUNDLE_BYTES:
            fail("completed timing evidence exceeds the bounded bundle size")
        snapshots[key] = data
        identities[key] = identity
    return snapshots, identities


def _open_pinned_completed_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Dict[str, int],
        ) -> Tuple[Dict[str, bytes], Dict[str, CompletedFingerprint]]:
    """Snapshot and retain every dependency inode across semantic validation."""
    if not isinstance(descriptors, dict) or descriptors:
        fail("pinned completed timing descriptor holder is invalid")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, CompletedFingerprint] = {}
    total_bytes = 0
    try:
        for key, name in source_names.items():
            limit = COMPLETED_ARTIFACT_BYTE_LIMITS.get(name)
            if type(limit) is not int or limit <= 0:
                fail("completed timing bundle names an unbounded artifact")
            data, fingerprint = _open_completed_regular_snapshot(
                Path(name), "completed timing " + key.replace("_", " "),
                directory_fd, limit, descriptors, key)
            snapshots[key] = data
            fingerprints[key] = fingerprint
            total_bytes += len(data)
            if total_bytes > MAX_COMPLETED_BUNDLE_BYTES:
                fail("completed timing evidence exceeds the bounded bundle size")
        return snapshots, fingerprints
    except BaseException as primary:
        pending = [
            ("completed dependency {} descriptor cleanup failed".format(key),
             descriptor)
            for key, descriptor in descriptors.items()
        ]
        descriptors.clear()
        _drain_descriptor_closes(pending, primary)
        raise


def _reread_pinned_completed_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected_fingerprints: Mapping[str, CompletedFingerprint],
        ) -> Tuple[Dict[str, bytes], Dict[str, CompletedFingerprint]]:
    """Reread retained inodes and prove every public name still points to one."""
    if (set(descriptors) != set(source_names) or
            set(expected_fingerprints) != set(source_names)):
        fail("pinned completed timing bundle is incomplete")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, CompletedFingerprint] = {}
    total_bytes = 0
    for key, name in source_names.items():
        context = "completed timing " + key.replace("_", " ")
        data, fingerprint = _read_completed_descriptor_bytes(
            descriptors[key], context, COMPLETED_ARTIFACT_BYTE_LIMITS[name])
        try:
            named = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            fail("cannot terminally inspect {}: {}".format(context, exc))
        if (not stat.S_ISREG(named.st_mode) or
                fingerprint != expected_fingerprints[key] or
                _completed_fingerprint(named) != fingerprint):
            fail("{} changed during semantic validation".format(context))
        snapshots[key] = data
        fingerprints[key] = fingerprint
        total_bytes += len(data)
        if total_bytes > MAX_COMPLETED_BUNDLE_BYTES:
            fail("completed timing evidence exceeds the bounded bundle size")
    return snapshots, fingerprints


def _require_pinned_completed_bundle_unchanged(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected_fingerprints: Mapping[str, CompletedFingerprint]) -> None:
    """Cheap post-hard-wall identity/metadata check immediately before link."""
    if (set(descriptors) != set(source_names) or
            set(expected_fingerprints) != set(source_names)):
        fail("pinned completed timing bundle is incomplete")
    for key, name in source_names.items():
        try:
            retained = os.fstat(descriptors[key])
            named = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            fail("cannot terminally inspect completed timing {}: {}".format(
                key.replace("_", " "), exc))
        expected = expected_fingerprints[key]
        if (not stat.S_ISREG(retained.st_mode) or
                not stat.S_ISREG(named.st_mode) or
                _completed_fingerprint(retained) != expected or
                _completed_fingerprint(named) != expected):
            fail("completed timing {} changed before publication".format(
                key.replace("_", " ")))


def _verify_completed_directory_path(
        path: Path, expected_identity: Tuple[int, int]) -> None:
    try:
        info = os.stat(str(path), follow_symlinks=False)
    except OSError as exc:
        fail("cannot recheck completed timing directory {}: {}".format(
            path, exc))
    if (not stat.S_ISDIR(info.st_mode) or
            (info.st_dev, info.st_ino) != expected_identity):
        fail("completed timing directory identity changed during validation")


def _commit_completed_timing_screen(
        contract: Mapping[str, Any], output_dir: Path,
        summary: Mapping[str, Any], freeze: Mapping[str, Any],
        timing_summary: Mapping[str, Any],
        execution_receipt: Mapping[str, Any],
        timing_qualification: contract_api.TimingQualification,
        finish_hard_wall: Optional[Callable[[], None]] = None) -> None:
    """Validate privately, finish the hard wall, then publish exactly once."""
    summary_path = output_dir / "run-summary.json"
    summary_bytes = (
        contract_api.canonical_json(summary) + "\n").encode("utf-8")
    parent_fd = -1
    parent_fd_holder = [-1]
    summary_fd_holder = [-1]
    dependency_fds: Dict[str, int] = {}
    parent_identity: Optional[Tuple[int, int]] = None
    try:
        parent_fd, parent_identity = _open_completion_parent(
            summary_path, parent_fd_holder)
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        _verify_completed_directory_path(output_dir, parent_identity)
        _require_completion_marker_absent(summary_path, parent_fd)
        summary_fingerprint = _open_unnamed_completion_marker(
            parent_fd, summary_bytes, summary_fd_holder)
        summary_fd = summary_fd_holder[0]

        dependency_snapshots, dependency_fingerprints = \
            _open_pinned_completed_bundle(
                COMPLETED_DEPENDENCY_NAMES, parent_fd, dependency_fds)
        if (len(summary_bytes) + sum(
                len(value) for value in dependency_snapshots.values()) >
                MAX_COMPLETED_BUNDLE_BYTES):
            fail("completed timing dependencies exceed the bundle byte cap")
        prospective_snapshots = dict(dependency_snapshots)
        prospective_snapshots["summary"] = summary_bytes
        resolved = output_dir.resolve(strict=True)
        resolved_info = resolved.stat()
        if ((resolved_info.st_dev, resolved_info.st_ino) != parent_identity or
                not stat.S_ISDIR(resolved_info.st_mode)):
            fail("completed timing directory identity changed before validation")
        validated = _validate_completed_timing_snapshots(
            contract, resolved, parent_identity, prospective_snapshots)
        if (not isinstance(validated, dict) or set(validated) != {
                "directory", "directory_identity", "run_summary", "freeze",
                "summary", "execution_receipt", "timing_qualification",
                } or
                validated["directory_identity"] != parent_identity or
                contract_api.canonical_json(validated["run_summary"]) !=
                contract_api.canonical_json(summary) or
                contract_api.canonical_json(validated["freeze"]) !=
                contract_api.canonical_json(freeze) or
                contract_api.canonical_json(validated["summary"]) !=
                contract_api.canonical_json(timing_summary) or
                contract_api.canonical_json(validated["execution_receipt"]) !=
                contract_api.canonical_json(execution_receipt) or
                validated["timing_qualification"].map_sha256 !=
                timing_qualification.map_sha256 or
                validated["timing_qualification"].qualified_domain_sha256 !=
                timing_qualification.qualified_domain_sha256):
            fail("prospective timing validation differs from terminal evidence")

        # Semantic validation can be expensive.  Before exposing a marker,
        # reread both the anonymous summary inode and every named dependency,
        # requiring exact bytes and entry fingerprints from the pinned parent.
        final_summary_bytes, final_summary_fingerprint = \
            _read_unnamed_completion_marker(
                summary_fd,
                COMPLETED_ARTIFACT_BYTE_LIMITS["run-summary.json"])
        final_dependencies, final_dependency_fingerprints = \
            _reread_pinned_completed_bundle(
                COMPLETED_DEPENDENCY_NAMES, parent_fd, dependency_fds,
                dependency_fingerprints)
        if (final_summary_bytes != summary_bytes or
                final_summary_fingerprint != summary_fingerprint):
            fail("prospective completion marker changed during validation")
        if (final_dependencies != dependency_snapshots or
                final_dependency_fingerprints != dependency_fingerprints):
            fail("completed timing evidence changed during semantic validation")
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        _verify_completed_directory_path(output_dir, parent_identity)
        _require_completion_marker_absent(summary_path, parent_fd)

        # This is the last potentially unbounded pre-publication hook.  No
        # public marker exists if it raises.  Only bounded identity checks
        # follow it before the publish-wins link(2): after that succeeds, no
        # exception or durability uncertainty may remove the public name.
        if finish_hard_wall is not None:
            finish_hard_wall()
        try:
            post_finish_summary = os.fstat(summary_fd)
        except OSError as exc:
            fail("cannot recheck unnamed completion marker: {}".format(exc))
        if _completed_fingerprint(post_finish_summary) != summary_fingerprint:
            fail("prospective completion marker changed before publication")
        _require_pinned_completed_bundle_unchanged(
            COMPLETED_DEPENDENCY_NAMES, parent_fd, dependency_fds,
            dependency_fingerprints)
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        _verify_completed_directory_path(output_dir, parent_identity)
        _require_completion_marker_absent(summary_path, parent_fd)
        _link_unnamed_completion_marker(
            summary_path, summary_fd, parent_fd, parent_identity)
    finally:
        primary = sys.exc_info()[1]
        pending = [
            ("completed dependency {} descriptor cleanup failed".format(key),
             descriptor)
            for key, descriptor in dependency_fds.items()
        ]
        pending.extend((
            ("completion-marker descriptor cleanup failed",
             summary_fd_holder[0]),
            ("completion-parent descriptor cleanup failed",
             parent_fd_holder[0]),
        ))
        dependency_fds.clear()
        summary_fd_holder[0] = -1
        parent_fd_holder[0] = -1
        # Before link(2), closing the anonymous inode removes it.  After link,
        # descriptor cleanup can never revoke the publish-wins commit.
        _drain_descriptor_closes(pending, primary)


def _write_completed_snapshot(path: Path, data: bytes) -> None:
    try:
        with path.open("xb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except OSError as exc:
        fail("cannot write private timing snapshot {}: {}".format(path, exc))


def _validate_completed_timing_snapshots(
        contract: Mapping[str, Any], resolved: Path,
        directory_identity: Tuple[int, int], snapshots: Mapping[str, bytes],
        ) -> Mapping[str, Any]:
    """Apply the exact loader semantics to already-pinned immutable bytes."""
    if (not isinstance(resolved, Path) or
            not isinstance(directory_identity, tuple) or
            len(directory_identity) != 2 or
            set(snapshots) != set(COMPLETED_SOURCE_NAMES) or
            any(not isinstance(value, bytes) for value in snapshots.values())):
        fail("completed timing snapshot set is invalid")
    total_bytes = 0
    for key, name in COMPLETED_SOURCE_NAMES.items():
        value = snapshots[key]
        limit = COMPLETED_ARTIFACT_BYTE_LIMITS[name]
        if len(value) > limit:
            fail("completed timing snapshot exceeds its artifact byte cap")
        total_bytes += len(value)
    if total_bytes > MAX_COMPLETED_BUNDLE_BYTES:
        fail("completed timing snapshot exceeds its bundle byte cap")

    with tempfile.TemporaryDirectory(prefix="wh2-timing-snapshot-") as raw:
        snapshot_root = Path(raw)
        paths = {key: snapshot_root / name
                 for key, name in COMPLETED_SOURCE_NAMES.items()}
        for key, path in paths.items():
            _write_completed_snapshot(path, snapshots[key])

        run_summary = _parse_canonical_line(
            snapshots["summary"], "completed timing run summary")
        if (set(run_summary) != SUMMARY_FIELDS or
                run_summary.get("schema") != SUMMARY_SCHEMA or
                run_summary.get("status") != "complete" or
                run_summary.get("output_dir") != str(resolved)):
            fail("completed timing run summary is incomplete or misplaced")
        unsigned_summary = {
            key: value for key, value in run_summary.items()
            if key != "summary_sha256"
        }
        if (not _is_sha256(run_summary.get("summary_sha256")) or
                run_summary["summary_sha256"] == ZERO_SHA256 or
                run_summary["summary_sha256"] !=
                    contract_api.sha256_json(unsigned_summary)):
            fail("completed timing run summary self-hash is invalid")
        qualification_map_sha256 = run_summary.get(
            "timing_qualification_map_sha256")
        qualification_execution_sha256 = run_summary.get(
            "timing_qualification_execution_receipt_sha256")
        if (not _is_sha256(qualification_map_sha256) or
                qualification_map_sha256 == ZERO_SHA256 or
                not _is_sha256(qualification_execution_sha256) or
                qualification_execution_sha256 == ZERO_SHA256):
            fail("completed timing qualification hashes are invalid")
        try:
            qualification = contract_api.load_timing_qualification_map(
                contract, "development", paths["qualification_map"],
                paths["qualification_audit"], qualification_map_sha256)
            validated = native_api.validate_execution_receipt(
                contract, "timing", "development", paths["freeze"],
                paths["trace"], paths["native"], paths["result"],
                paths["receipt"], verify_live_sampler=False,
                sampler_path=paths["sampler"],
                timing_qualification_map_path=paths["qualification_map"],
                timing_qualification_audit_path=paths["qualification_audit"],
                timing_qualification_map_sha256=qualification_map_sha256,
                timing_qualification_native_path=
                    paths["qualification_native"],
                timing_qualification_sampler_path=
                    paths["qualification_sampler"],
                timing_qualification_execution_receipt_path=
                    paths["qualification_receipt"],
                timing_qualification_execution_receipt_sha256=
                    qualification_execution_sha256)
            freeze = contract_api.load_freeze_manifest(
                contract, "development", paths["freeze"], "timing",
                qualification)
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))

        if freeze["arm_roster"] != [value[0] for value in EXPECTED_ARMS]:
            fail("completed timing freeze has the wrong exact arm roster")
        for index, (expected_arm, expected_codec) in enumerate(EXPECTED_ARMS):
            arm = freeze["arms"][index]
            if (arm["arm"] != expected_arm or
                    arm["codec"] != expected_codec or
                    arm["binary_sha256"] !=
                        freeze["arms"][0]["binary_sha256"]):
                fail("completed timing freeze substitutes an arm artifact")
        if freeze["arms"][2]["arm_descriptor_sha256"] != \
                TIMING_CANDIDATE_DESCRIPTOR_SHA256:
            fail("completed timing freeze substitutes the Two07 descriptor")

        timing_summary = validated["summary"]
        execution_receipt = validated["execution_receipt"]
        retry_offsets = list(qualification.retry_offsets)
        thermal = execution_receipt["thermal"]
        qualification_thermal = execution_receipt["qualification_thermal"]
        host_identity = freeze.get("host_identity")
        controller_cpu = host_identity.get("controller_cpu") \
            if isinstance(host_identity, dict) else None
        if type(controller_cpu) is not int or controller_cpu < 0:
            fail("completed timing freeze lacks a valid controller CPU")
        expected = {
            "source_git_commit": freeze["source_git_commit"],
            "contract_sha256": contract_api.contract_sha256(contract),
            "worker_binary_sha256": freeze["arms"][0]["binary_sha256"],
            "controller_cpu": controller_cpu,
            "worker_cpus": execution_receipt["worker_cpus"],
            "qualification_worker_cpus":
                execution_receipt["qualification_worker_cpus"],
            "timing_qualification_map_sha256": qualification.map_sha256,
            "timing_qualification_execution_receipt_sha256":
                execution_receipt[
                    "timing_qualification_execution_receipt_sha256"],
            "timing_qualified_domain_sha256":
                qualification.qualified_domain_sha256,
            "qualification_attempt_count":
                len(retry_offsets) + sum(retry_offsets),
            "qualification_retried_cell_count": sum(
                1 for retry in retry_offsets if retry > 0),
            "qualification_max_retry_offset": max(retry_offsets, default=0),
            "qualification_sum_retry_offsets": sum(retry_offsets),
            "timing_records": execution_receipt["record_count"],
            "timing_trace_manifest_sha256":
                freeze["trace_manifest_sha256"],
            "timing_freeze_sha256":
                execution_receipt["freeze_manifest_sha256"],
            "timing_architecture_artifact_sha256":
                contract_api.architecture_artifact_sha256(freeze),
            "timing_result_sha256":
                execution_receipt["result_stream_sha256"],
            "timing_execution_receipt_sha256":
                execution_receipt["receipt_sha256"],
            "timing_validator_summary_sha256":
                execution_receipt["validator_summary_sha256"],
            "thermal_samples": thermal["sample_count"],
            "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "dimm_max_millic": thermal["dimm_max_millic"],
            "timing_thermal_samples": thermal["sample_count"],
            "timing_cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "timing_dimm_max_millic": thermal["dimm_max_millic"],
            "qualification_thermal_samples":
                qualification_thermal["sample_count"],
            "qualification_cpu_tctl_max_millic":
                qualification_thermal["cpu_tctl_max_millic"],
            "qualification_dimm_max_millic":
                qualification_thermal["dimm_max_millic"],
            "overall_cpu_tctl_max_millic": max(
                thermal["cpu_tctl_max_millic"],
                qualification_thermal["cpu_tctl_max_millic"]),
            "overall_dimm_max_millic": max(
                thermal["dimm_max_millic"],
                qualification_thermal["dimm_max_millic"]),
        }
        if (any(contract_api.canonical_json(run_summary.get(field)) !=
                contract_api.canonical_json(value)
                for field, value in expected.items()) or
                run_summary["timing_validator_summary_sha256"] !=
                    contract_api.sha256_json(timing_summary) or
                qualification.source_git_commit !=
                    freeze["source_git_commit"]):
            fail("completed timing summary differs from reopened evidence")

        return {
            "directory": str(resolved),
            "directory_identity": directory_identity,
            "run_summary": run_summary,
            "freeze": freeze,
            "summary": timing_summary,
            "execution_receipt": execution_receipt,
            "timing_qualification": qualification,
        }


def load_completed_timing_screen(
        contract: Mapping[str, Any], timing_screen_dir: Path,
        ) -> Mapping[str, Any]:
    """Strictly reopen one terminal timing-only development evidence bundle."""
    if not isinstance(timing_screen_dir, Path):
        fail("timing screen directory must be a pathlib Path")
    directory_fd = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        directory_flag = getattr(os, "O_DIRECTORY", 0)
        if nofollow == 0 or directory_flag == 0:
            fail("timing screen directory cannot be opened fail-closed")
        directory_fd = os.open(
            str(timing_screen_dir),
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory_flag)
        opened_info = os.fstat(directory_fd)
        if not stat.S_ISDIR(opened_info.st_mode):
            fail("timing screen path must be a real directory, not a symlink")
        resolved = timing_screen_dir.resolve(strict=True)
        resolved_info = resolved.stat()
        if (opened_info.st_dev, opened_info.st_ino) != \
                (resolved_info.st_dev, resolved_info.st_ino):
            fail("timing screen directory identity changed while opening")
    except (OSError, RuntimeError) as exc:
        primary = RunnerError(
            "cannot open timing screen directory {}: {}".format(
                timing_screen_dir, exc))
        descriptor = directory_fd
        directory_fd = -1
        _drain_descriptor_closes((
            ("completed timing directory descriptor cleanup failed",
             descriptor),
        ), primary)
        raise primary from exc
    except BaseException as primary:
        descriptor = directory_fd
        directory_fd = -1
        _drain_descriptor_closes((
            ("completed timing directory descriptor cleanup failed",
             descriptor),
        ), primary)
        raise
    source_names = COMPLETED_SOURCE_NAMES
    try:
        snapshots, snapshot_identities = _read_completed_bundle(
            source_names, directory_fd)
    except BaseException as primary:
        descriptor = directory_fd
        directory_fd = -1
        _drain_descriptor_closes((
            ("completed timing directory descriptor cleanup failed",
             descriptor),
        ), primary)
        raise

    try:
        validated_result = _validate_completed_timing_snapshots(
            contract, resolved, (opened_info.st_dev, opened_info.st_ino),
            snapshots)
    except BaseException as primary:
        descriptor = directory_fd
        directory_fd = -1
        _drain_descriptor_closes((
            ("completed timing directory descriptor cleanup failed",
             descriptor),
        ), primary)
        raise

    # Semantic validation may be expensive.  Reopen the exact directory inode,
    # reread the whole bundle through a new pinned descriptor, and require both
    # bytes and directory-entry fingerprints to equal the initial stable read.
    terminal_directory_fd = -1
    expected_directory_identity = (opened_info.st_dev, opened_info.st_ino)
    try:
        terminal_directory_fd = os.open(
            str(resolved),
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) |
            getattr(os, "O_DIRECTORY", 0))
        terminal_directory_info = os.fstat(terminal_directory_fd)
        if (not stat.S_ISDIR(terminal_directory_info.st_mode) or
                (terminal_directory_info.st_dev,
                 terminal_directory_info.st_ino) !=
                expected_directory_identity):
            fail("completed timing directory changed before terminal reread")
        _verify_completed_directory_path(
            resolved, expected_directory_identity)
        terminal_snapshots, terminal_identities = _read_completed_bundle(
            source_names, terminal_directory_fd)
        _verify_completed_directory_path(
            resolved, expected_directory_identity)
    except OSError as exc:
        fail("cannot terminally reopen completed timing directory: {}".
             format(exc))
    finally:
        primary = sys.exc_info()[1]
        pending = (
            ("terminal completed timing directory descriptor cleanup failed",
             terminal_directory_fd),
            ("completed timing directory descriptor cleanup failed",
             directory_fd),
        )
        terminal_directory_fd = -1
        directory_fd = -1
        _drain_descriptor_closes(pending, primary)
    if (terminal_snapshots != snapshots or
            terminal_identities != snapshot_identities):
        fail("completed timing evidence changed during semantic validation")
    return validated_result


def _resume_managed_sibling_cleanup(
        managed_paused: ManagedPausedSiblingProcesses,
        cleanup_deadline: float) -> List[CleanupFailure]:
    """Converge managed SIGCONT ownership despite one-shot interruptions."""
    failures: List[CleanupFailure] = []
    complete = False
    attempts = 0
    while attempts < 8:
        try:
            complete = managed_paused.cleanup_complete
        except BaseException as cleanup:
            failures.append((
                "managed sibling cleanup-state inspection failed", cleanup))
            complete = False
        if complete:
            break

        attempts += 1
        try:
            managed_paused.resume(cleanup_deadline)
        except BaseException as cleanup:
            failures.append((
                "managed sibling load resume attempt {} failed".format(
                    attempts), cleanup))

        try:
            complete = managed_paused.cleanup_complete
        except BaseException as cleanup:
            failures.append((
                "managed sibling cleanup-state inspection failed", cleanup))
            complete = False
        if complete:
            break

        # Even an expired grace period cannot suppress the second whole call:
        # a one-shot signal may have landed before resume() entered its own
        # retry state machine.  Later attempts remain bounded by both a small
        # fixed cap and the fresh cleanup deadline.
        if attempts >= 2 and time.monotonic() >= cleanup_deadline:
            break

    if not complete:
        failures.append((
            "managed sibling load cleanup did not converge",
            RunnerError(
                "managed sibling load still owns resume state after {} "
                "cleanup attempts".format(attempts))))
    return failures


def _finish_timing_cleanup(
        workers: Sequence[PersistentWorker], clean_shutdown: bool,
        controller_pinned: bool, original_affinity: Set[int],
        primary: Optional[BaseException],
        managed_paused: Optional[ManagedPausedSiblingProcesses] = None,
        paused_resumed: bool = False,
        deadline: Optional[float] = None) -> None:
    """Attempt every cleanup and preserve ordinary and control-flow failures."""
    cleanup_failures: List[Tuple[str, BaseException]] = []
    if workers and not clean_shutdown:
        try:
            terminate_workers(workers)
        except BaseException as cleanup:
            cleanup_failures.append(("timing worker cleanup failed", cleanup))
    if controller_pinned:
        try:
            _restore_controller_affinity(original_affinity)
        except BaseException as cleanup:
            cleanup_failures.append(
                ("controller affinity cleanup failed", cleanup))
    if managed_paused is not None and not paused_resumed:
        try:
            # Cleanup owns a fresh grace period even when the campaign hard
            # wall caused the primary failure.  SIGCONT must never inherit an
            # already-expired benchmark deadline and strand stopped load.
            cleanup_deadline = max(
                time.monotonic() + 5.0,
                deadline if deadline is not None else 0.0)
            cleanup_failures.extend(_resume_managed_sibling_cleanup(
                managed_paused, cleanup_deadline))
        except BaseException as cleanup:
            cleanup_failures.append(
                ("managed sibling load cleanup failed", cleanup))
    _raise_after_cleanup(primary, cleanup_failures)


def run_short_screen(
        args: argparse.Namespace,
        finish_hard_wall: Optional[Callable[[], None]] = None,
        ) -> Mapping[str, Any]:
    if (not math.isfinite(args.deadline_seconds) or
            not 0.0 < args.deadline_seconds <= MAX_WALL_SECONDS):
        fail("--deadline-seconds must be in (0,7200]")
    deadline = time.monotonic() + args.deadline_seconds
    contract = contract_api.load_contract(args.contract)
    description = describe_worker(args.worker, deadline)
    try:
        original_affinity = set(os.sched_getaffinity(0))
    except (AttributeError, OSError) as exc:
        fail("cannot inspect initial controller affinity: {}".format(exc))
    qualification_cpus = sorted(original_affinity)
    if not qualification_cpus:
        fail("qualification requires a nonempty logical CPU affinity")
    sampler_core = _cpu_topology(
        args.sampler_cpu, CPU_SYSFS_ROOT)
    if any(_cpu_topology(
            cpu, CPU_SYSFS_ROOT) == sampler_core
           for cpu in qualification_cpus):
        fail("qualification requires every allowed logical CPU and the sampler "
             "physical core must be outside that affinity")
    _development_timing_shape(
        contract,
        len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS])))
    _preflight_sampler(
        args.sampler_pid, args.sampler_cpu,
        args.sampler_script, args.sampler_csv)
    source_commit = _git_head(deadline)
    _require_worker_source_commit(description, source_commit)
    output_dir = _create_output_dir(args.output_dir)

    controls = _timing_qualification_controls(contract, description)
    qualification_native_path = output_dir / \
        "timing-qualification-native.jsonl"
    qualification_audit_path = output_dir / \
        "timing-qualification-audit.jsonl"
    qualification_map_path = output_dir / "timing-qualification-map.json"
    qualification_execution_receipt_path = output_dir / \
        "timing-qualification-execution.json"
    qualification_window_start_ns = choose_new_sampler_start(
        args.sampler_csv, deadline)

    qualification_workers: List[PersistentWorker] = []
    qualification_clean_shutdown = False
    qualification_maximum_end = 0
    try:
        # Qualification intentionally saturates the complete logical affinity,
        # including SMT siblings.  Timing topology isolation is selected only
        # after this entire pool has exited and been reaped.
        qualification_workers = spawn_workers(
            description, qualification_cpus, deadline)
        _, qualification_maximum_end = run_timing_qualification(
            contract, description, qualification_workers,
            qualification_native_path, qualification_window_start_ns,
            deadline)
        try:
            (timing_qualification, qualification_metadata,
             qualification_execution_receipt_sha256) = \
                native_api.assemble_timing_qualification(
                    contract, "development", qualification_native_path,
                    qualification_audit_path, qualification_map_path,
                    qualification_execution_receipt_path, source_commit,
                    controls, qualification_cpus,
                    verify_live_workers=True)
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        quit_workers(qualification_workers, deadline)
        qualification_clean_shutdown = True
    finally:
        _finish_worker_pool_cleanup(
            qualification_workers, qualification_clean_shutdown,
            "qualification worker cleanup failed", sys.exc_info()[1])
    if any(worker.process.poll() is None for worker in qualification_workers):
        fail("a qualification worker survived the mandatory reap boundary")
    qualification_window_end_ns = _wait_for_sampler_sample(
        args.sampler_csv, deadline,
        at_or_after_ns=qualification_maximum_end)
    qualification_sampler_path = output_dir / \
        "qualification-sampler-attestation.json"
    try:
        native_api.write_sampler_attestation(
            args.sampler_pid, args.sampler_cpu,
            args.sampler_script, args.sampler_csv,
            qualification_window_start_ns, qualification_window_end_ns,
            qualification_sampler_path)
    except native_api.NativeEvidenceError as exc:
        fail(str(exc))

    qualification_map_sha256 = timing_qualification.map_sha256
    try:
        timing_qualification = contract_api.load_timing_qualification_map(
            contract, "development", qualification_map_path,
            qualification_audit_path, qualification_map_sha256)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if _git_head(deadline) != source_commit:
        fail("codec source commit changed during timing qualification")

    # Exact-eight topology selection is deliberately after the qualification
    # process boundary.  No qualification worker/cache can survive into T.
    explicit = parse_cpu_list(args.cpus) if args.cpus is not None else None
    cpus, controller_cpu = select_cpu_layout(
        args.sampler_cpu, explicit, args.controller_cpu,
        affinity=original_affinity)
    if len(cpus) != TIMING_WORKER_COUNT:
        fail("native short screen requires exactly eight physical worker CPUs")
    timing_trace_path = output_dir / "timing-traces.jsonl"
    try:
        timing_trace_sha256 = native_api.publish_timing_trace_manifest(
            contract, "development", timing_qualification,
            timing_trace_path)
    except (native_api.NativeEvidenceError,
            contract_api.ContractError) as exc:
        fail(str(exc))
    if _git_head(deadline) != source_commit:
        fail("codec source commit changed before the result freeze")
    freeze = write_development_timing_freeze(
        contract, description, cpus, controller_cpu, source_commit,
        timing_trace_sha256, output_dir, timing_qualification)
    freeze_path = output_dir / "timing-freeze.json"

    workers: List[PersistentWorker] = []
    clean_shutdown = False
    controller_pinned = False
    managed_paused: Optional[ManagedPausedSiblingProcesses] = \
        ManagedPausedSiblingProcesses()
    paused_resumed = False
    sibling_guard: Optional[SiblingIsolationGuard] = None
    cleanup_primary: Optional[BaseException] = None
    cleanup_capture_complete = False
    try:
        try:
            # Establish caller-visible affinity cleanup ownership before the first
            # isolation scan.  The controller is therefore exact-pinned at every
            # attested stage, including a cpuset containing only T cores/siblings.
            controller_pinned = True
            _pin_controller(controller_cpu)
            pause_text = getattr(args, "pause_sibling_pids", None)
            pause_pids = parse_pid_list(pause_text) \
                if pause_text is not None else []
            protected_cpus = sorted(list(cpus) + [controller_cpu])
            sibling_cpus = native_api.timing_sibling_cpus(protected_cpus)
            managed_paused.pause(pause_pids, sibling_cpus, deadline)
            sibling_guard = SiblingIsolationGuard(
                cpus, controller_cpu, managed_paused.records)
            sibling_guard.check("before-worker-spawn")
            # The validated timing freeze exists before the first T command.
            workers = spawn_workers(description, cpus, deadline)
            sibling_guard.bind_workers(workers)
            sibling_guard.check("before-timing-window")
            timing_window_start_ns = choose_new_sampler_start(
                args.sampler_csv, deadline)
            native_path, maximum_worker_end = _run_timing_jobs(
                contract, freeze, timing_qualification, description, workers,
                output_dir,
                timing_window_start_ns, deadline, sibling_guard)
            sibling_guard.check("after-timing-interval")
            sibling_isolation = sibling_guard.attestation()
            window_end_ns = _wait_for_sampler_sample(
                args.sampler_csv, deadline,
                at_or_after_ns=maximum_worker_end)
            if _git_head(deadline) != source_commit:
                fail("codec source commit changed during the native screen")
            sampler_path = output_dir / "sampler-attestation.json"
            try:
                native_api.write_sampler_attestation(
                    args.sampler_pid, args.sampler_cpu,
                    args.sampler_script, args.sampler_csv,
                    timing_window_start_ns, window_end_ns, sampler_path)
            except native_api.NativeEvidenceError as exc:
                fail(str(exc))

            result_path = output_dir / "timing-results.jsonl"
            receipt_path = output_dir / "timing-execution.json"
            try:
                assembled = native_api.assemble_results(
                    contract, "timing", "development", freeze_path,
                    timing_trace_path, native_path, sampler_path,
                    result_path, receipt_path, verify_live_sampler=True,
                    timing_qualification_map_path=qualification_map_path,
                    timing_qualification_audit_path=qualification_audit_path,
                    timing_qualification_map_sha256=qualification_map_sha256,
                    timing_qualification_native_path=qualification_native_path,
                    timing_qualification_sampler_path=qualification_sampler_path,
                    timing_qualification_execution_receipt_path=
                        qualification_execution_receipt_path,
                    timing_qualification_execution_receipt_sha256=
                        qualification_execution_receipt_sha256,
                    sibling_isolation=sibling_isolation)
                _remaining(deadline, "assembling timing results")
                terminal = native_api.validate_execution_receipt(
                    contract, "timing", "development", freeze_path,
                    timing_trace_path, native_path, result_path,
                    receipt_path, verify_live_sampler=True,
                    sampler_path=sampler_path,
                    timing_qualification_map_path=qualification_map_path,
                    timing_qualification_audit_path=qualification_audit_path,
                    timing_qualification_map_sha256=qualification_map_sha256,
                    timing_qualification_native_path=qualification_native_path,
                    timing_qualification_sampler_path=qualification_sampler_path,
                    timing_qualification_execution_receipt_path=
                        qualification_execution_receipt_path,
                    timing_qualification_execution_receipt_sha256=
                        qualification_execution_receipt_sha256)
                _remaining(deadline, "validating timing execution")
            except (native_api.NativeEvidenceError,
                    contract_api.ContractError) as exc:
                fail(str(exc))
            if contract_api.canonical_json(terminal) != \
                    contract_api.canonical_json(assembled):
                fail("terminal timing execution differs from assembled evidence")
            quit_workers(workers, deadline)
            clean_shutdown = True
            if controller_pinned:
                _restore_controller_affinity(original_affinity)
                controller_pinned = False
            if managed_paused is not None:
                managed_paused.resume(deadline)
                paused_resumed = True

            # Q/reap and affinity restoration are part of the evidence boundary.
            # Revalidate every bound dependency after both have completed so a
            # mutation in that interval cannot be covered by a complete marker.
            if _git_head(deadline) != source_commit:
                fail("codec source commit changed before terminal publication")
            try:
                post_shutdown = native_api.validate_execution_receipt(
                    contract, "timing", "development", freeze_path,
                    timing_trace_path, native_path, result_path,
                    receipt_path, verify_live_sampler=False,
                    sampler_path=sampler_path,
                    timing_qualification_map_path=qualification_map_path,
                    timing_qualification_audit_path=qualification_audit_path,
                    timing_qualification_map_sha256=qualification_map_sha256,
                    timing_qualification_native_path=qualification_native_path,
                    timing_qualification_sampler_path=qualification_sampler_path,
                    timing_qualification_execution_receipt_path=
                        qualification_execution_receipt_path,
                    timing_qualification_execution_receipt_sha256=
                        qualification_execution_receipt_sha256)
                post_shutdown_freeze = contract_api.load_freeze_manifest(
                    contract, "development", freeze_path, "timing",
                    timing_qualification)
                _remaining(deadline, "post-shutdown timing validation")
            except (native_api.NativeEvidenceError,
                    contract_api.ContractError) as exc:
                fail(str(exc))
            if (contract_api.canonical_json(post_shutdown) !=
                    contract_api.canonical_json(assembled) or
                    contract_api.canonical_json(post_shutdown_freeze) !=
                    contract_api.canonical_json(freeze)):
                fail("post-shutdown timing evidence differs from assembled evidence")
            if _git_head(deadline) != source_commit:
                fail("codec source commit changed during terminal validation")

            timing_receipt = post_shutdown["execution_receipt"]
            timing_validator_summary = post_shutdown["summary"]
            thermal = timing_receipt["thermal"]
            qualification_thermal = timing_receipt["qualification_thermal"]
            qualification_retry_offsets = list(
                timing_qualification.retry_offsets)
            summary = {
                "schema": SUMMARY_SCHEMA,
                "status": "complete",
                "output_dir": str(output_dir),
                "source_git_commit": source_commit,
                "contract_sha256": contract_api.contract_sha256(contract),
                "worker_binary_sha256": description["binary_sha256"],
                "controller_cpu": controller_cpu,
                "worker_cpus": list(cpus),
                "qualification_worker_cpus": qualification_cpus,
                "timing_qualification_map_sha256": qualification_map_sha256,
                "timing_qualification_execution_receipt_sha256":
                    qualification_execution_receipt_sha256,
                "timing_qualified_domain_sha256":
                    timing_qualification.qualified_domain_sha256,
                "qualification_attempt_count":
                    qualification_metadata["qualification_attempt_count"],
                "qualification_retried_cell_count": sum(
                    1 for retry in qualification_retry_offsets if retry > 0),
                "qualification_max_retry_offset": max(
                    qualification_retry_offsets, default=0),
                "qualification_sum_retry_offsets": sum(
                    qualification_retry_offsets),
                "timing_records": timing_receipt["record_count"],
                "timing_trace_manifest_sha256":
                    freeze["trace_manifest_sha256"],
                "timing_freeze_sha256": timing_receipt["freeze_manifest_sha256"],
                "timing_architecture_artifact_sha256":
                    contract_api.architecture_artifact_sha256(freeze),
                "timing_result_sha256": timing_receipt["result_stream_sha256"],
                "timing_execution_receipt_sha256":
                    timing_receipt["receipt_sha256"],
                "timing_validator_summary_sha256":
                    timing_receipt["validator_summary_sha256"],
                "thermal_samples": thermal["sample_count"],
                "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
                "dimm_max_millic": thermal["dimm_max_millic"],
                "timing_thermal_samples": thermal["sample_count"],
                "timing_cpu_tctl_max_millic":
                    thermal["cpu_tctl_max_millic"],
                "timing_dimm_max_millic": thermal["dimm_max_millic"],
                "qualification_thermal_samples":
                    qualification_thermal["sample_count"],
                "qualification_cpu_tctl_max_millic":
                    qualification_thermal["cpu_tctl_max_millic"],
                "qualification_dimm_max_millic":
                    qualification_thermal["dimm_max_millic"],
                "overall_cpu_tctl_max_millic": max(
                    thermal["cpu_tctl_max_millic"],
                    qualification_thermal["cpu_tctl_max_millic"]),
                "overall_dimm_max_millic": max(
                    thermal["dimm_max_millic"],
                    qualification_thermal["dimm_max_millic"]),
            }
            summary["summary_sha256"] = contract_api.sha256_json(summary)
            if (set(summary) != SUMMARY_FIELDS or
                    summary["timing_validator_summary_sha256"] !=
                        contract_api.sha256_json(timing_validator_summary)):
                fail("internal timing run summary schema mismatch")
            _commit_completed_timing_screen(
                contract, output_dir, summary, post_shutdown_freeze,
                timing_validator_summary, timing_receipt, timing_qualification,
                finish_hard_wall)
            # The successful transaction is the final fallible action.  Do not add
            # a deadline/provenance check after it without extending its rollback
            # boundary around that check.
            return summary
        finally:
            # Keep the first cleanup suite direct: on CPython 3.11+ a nested
            # try here starts with an uncovered exception-path NOP.  The direct
            # first opcode is protected by the finally unwinder and therefore
            # reaches the outer full-cleanup handler if it is interrupted.
            cleanup_primary = sys.exc_info()[1]
            cleanup_capture_complete = True
            # Retire the campaign's one-shot SIGALRM before any cleanup.  If
            # this call is interrupted, its handler disables the timer and the
            # outer fallback gets an interruption-free retry.
            if finish_hard_wall is not None:
                finish_hard_wall()
            _finish_timing_cleanup(
                workers, clean_shutdown, controller_pinned,
                original_affinity, cleanup_primary, managed_paused,
                paused_resumed, deadline)
    except BaseException as cleanup:
        # Capture can itself be interrupted before its assignment completes.
        # In that case Python chains the displaced body failure as context.
        active_primary = cleanup_primary
        if active_primary is None and not cleanup_capture_complete:
            active_primary = cleanup.__context__

        fallback_failures: List[CleanupFailure] = []
        if cleanup is not active_primary:
            fallback_failures.append((
                "timing cleanup interrupted", cleanup))
        try:
            if finish_hard_wall is not None:
                finish_hard_wall()
        except BaseException as disarm:
            fallback_failures.append((
                "timing hard-wall disarm fallback failed", disarm))
        try:
            # This is deliberately the whole cleanup, not only SIGCONT:
            # interruption may precede worker termination or affinity restore.
            _finish_timing_cleanup(
                workers, clean_shutdown, controller_pinned, original_affinity,
                None, managed_paused, paused_resumed, deadline)
        except BaseException as fallback:
            fallback_failures.append((
                "full timing cleanup fallback failed", fallback))

        if fallback_failures:
            _raise_after_cleanup(active_primary, fallback_failures)
            raise AssertionError("unreachable")
        if cleanup is active_primary:
            raise
        if active_primary is not None:
            raise active_primary
        raise


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument(
        "--worker", type=Path,
        default=Path("build/wirehair_wh2_contract_worker"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sampler-pid", type=int, required=True)
    parser.add_argument("--sampler-cpu", type=int, required=True)
    parser.add_argument("--sampler-script", type=Path, required=True)
    parser.add_argument("--sampler-csv", type=Path, required=True)
    parser.add_argument(
        "--cpus", help="strictly increasing comma-separated logical CPUs")
    parser.add_argument(
        "--controller-cpu", type=int,
        help="dedicated non-worker, non-sampler physical core")
    parser.add_argument(
        "--pause-sibling-pids",
        help="explicit controller-owned exact-affinity load PIDs to stop for T")
    parser.add_argument(
        "--deadline-seconds", type=float, default=MAX_WALL_SECONDS)
    args = parser.parse_args(argv if argv else None)
    try:
        if (not math.isfinite(args.deadline_seconds) or
                not 0.0 < args.deadline_seconds <= MAX_WALL_SECONDS):
            fail("--deadline-seconds must be in (0,7200]")
        with _hard_wall(args.deadline_seconds) as finish_hard_wall:
            summary = run_short_screen(args, finish_hard_wall)
    except (RunnerError, native_api.NativeEvidenceError,
            contract_api.ContractError, OSError, UnicodeError) as exc:
        print("wh2 native short screen: {}".format(exc), file=sys.stderr)
        return 1
    print(contract_api.canonical_json(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
