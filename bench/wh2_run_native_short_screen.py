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
import hashlib
import io
import math
import os
from pathlib import Path
import platform
import re
import selectors
import secrets
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


class RunnerError(RuntimeError):
    """The real native short screen cannot continue safely."""


def fail(message: str) -> None:
    raise RunnerError(message)


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


def _open_completion_parent(path: Path) -> Tuple[int, Tuple[int, int]]:
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    directory = getattr(os, "O_DIRECTORY", 0)
    if nofollow == 0 or directory == 0:
        fail("completion-marker publication requires O_NOFOLLOW/O_DIRECTORY")
    descriptor = -1
    try:
        descriptor = os.open(
            str(path.parent), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory)
        info = os.fstat(descriptor)
        if not stat.S_ISDIR(info.st_mode):
            fail("completion-marker parent must be a real directory")
        identity = (info.st_dev, info.st_ino)
        _verify_completed_directory_path(path.parent, identity)
        return descriptor, identity
    except OSError as exc:
        if descriptor >= 0:
            os.close(descriptor)
        fail("cannot open completion-marker parent {}: {}".format(
            path.parent, exc))
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
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


def _publish_completion_marker(
        path: Path, value: Mapping[str, Any], *,
        parent_fd: Optional[int] = None,
        parent_identity: Optional[Tuple[int, int]] = None,
        ) -> Tuple[int, int]:
    """Install one summary inode, rolling it back on every late failure."""
    path = Path(path)
    data = (contract_api.canonical_json(value) + "\n").encode("utf-8")
    limit = COMPLETED_ARTIFACT_BYTE_LIMITS["run-summary.json"]
    if not data or len(data) > limit:
        fail("completion marker exceeds its publication byte cap")
    owns_parent = parent_fd is None
    if owns_parent:
        if parent_identity is not None:
            fail("completion-marker parent guard is incomplete")
        parent_fd, parent_identity = _open_completion_parent(path)
    elif (type(parent_fd) is not int or parent_fd < 0 or
            not isinstance(parent_identity, tuple) or
            len(parent_identity) != 2):
        fail("completion-marker parent guard is invalid")
    assert parent_fd is not None
    assert parent_identity is not None
    staged_name: Optional[str] = None
    descriptor = -1
    staged_identity: Optional[Tuple[int, int]] = None
    published_identity: Optional[Tuple[int, int]] = None
    committed = False
    try:
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        _verify_completed_directory_path(path.parent, parent_identity)
        for _attempt in range(32):
            candidate = ".{}.{}.tmp".format(
                path.name, secrets.token_hex(12))
            try:
                descriptor = os.open(
                    candidate,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    0o600, dir_fd=parent_fd)
                staged_name = candidate
                break
            except FileExistsError:
                continue
        if descriptor < 0 or staged_name is None:
            fail("cannot allocate a unique staged completion marker")
        created = os.fstat(descriptor)
        if not stat.S_ISREG(created.st_mode):
            fail("staged completion marker is not a regular file")
        staged_identity = (created.st_dev, created.st_ino)
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                fail("short write while staging completion marker")
            offset += written
        os.fsync(descriptor)
        staged = os.fstat(descriptor)
        if (not stat.S_ISREG(staged.st_mode) or
                (staged.st_dev, staged.st_ino) != staged_identity or
                staged.st_size != len(data)):
            fail("staged completion marker is not the written regular file")
        published_identity = staged_identity
        try:
            os.link(
                "/proc/self/fd/{}".format(descriptor), path.name,
                dst_dir_fd=parent_fd, follow_symlinks=True)
        except FileExistsError:
            fail("refusing to replace existing artifact {}".format(path))
        published = os.stat(
            path.name, dir_fd=parent_fd, follow_symlinks=False)
        if (not stat.S_ISREG(published.st_mode) or
                (published.st_dev, published.st_ino) != published_identity):
            fail("published completion-marker identity changed")
        os.fsync(parent_fd)
        staged_entry = os.stat(
            staged_name, dir_fd=parent_fd, follow_symlinks=False)
        if (staged_entry.st_dev, staged_entry.st_ino) != published_identity:
            fail("staged completion-marker name changed during publication")
        os.unlink(staged_name, dir_fd=parent_fd)
        staged_name = None
        os.fsync(parent_fd)
        _verify_completed_directory_path(path.parent, parent_identity)
        committed = True
        return published_identity
    except OSError as exc:
        fail("cannot publish completion marker {}: {}".format(path, exc))
    finally:
        if descriptor >= 0:
            try:
                os.close(descriptor)
            except BaseException:
                pass
        if staged_name is not None:
            try:
                staged_entry = os.stat(
                    staged_name, dir_fd=parent_fd, follow_symlinks=False)
                if (staged_identity is not None and
                        (staged_entry.st_dev, staged_entry.st_ino) ==
                        staged_identity):
                    os.unlink(staged_name, dir_fd=parent_fd)
            except OSError:
                pass
        if not committed and published_identity is not None:
            try:
                visible = os.stat(
                    path.name, dir_fd=parent_fd, follow_symlinks=False)
                if (visible.st_dev, visible.st_ino) == published_identity:
                    os.unlink(path.name, dir_fd=parent_fd)
                    try:
                        os.fsync(parent_fd)
                    except OSError:
                        pass
            except OSError:
                pass
        if owns_parent:
            try:
                os.close(parent_fd)
            except BaseException:
                # All durability and identity checks precede committed=True.
                # A close error cannot invalidate an already committed marker.
                pass
    raise AssertionError("unreachable")


def _rollback_completion_marker(
        path: Path, identity: Tuple[int, int], *,
        parent_fd: Optional[int] = None,
        parent_identity: Optional[Tuple[int, int]] = None) -> None:
    """Remove only the exact summary inode installed by this transaction."""
    path = Path(path)
    owns_parent = parent_fd is None
    if owns_parent:
        if parent_identity is not None:
            fail("completion-marker rollback guard is incomplete")
        parent_fd, parent_identity = _open_completion_parent(path)
    elif (type(parent_fd) is not int or parent_fd < 0 or
            not isinstance(parent_identity, tuple) or
            len(parent_identity) != 2):
        fail("completion-marker rollback guard is invalid")
    assert parent_fd is not None
    assert parent_identity is not None
    try:
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        try:
            visible = os.stat(
                path.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            return
        if (visible.st_dev, visible.st_ino) != identity:
            return
        os.unlink(path.name, dir_fd=parent_fd)
        os.fsync(parent_fd)
        if owns_parent:
            _verify_completed_directory_path(path.parent, parent_identity)
    except OSError as exc:
        fail("cannot roll back completion marker {}: {}".format(path, exc))
    finally:
        if owns_parent:
            try:
                os.close(parent_fd)
            except BaseException:
                pass


def _require_completion_marker(
        path: Path, identity: Tuple[int, int], expected_bytes: bytes, *,
        parent_fd: Optional[int] = None,
        parent_identity: Optional[Tuple[int, int]] = None) -> None:
    path = Path(path)
    owns_parent = parent_fd is None
    if owns_parent:
        if parent_identity is not None:
            fail("completion-marker read guard is incomplete")
        parent_fd, parent_identity = _open_completion_parent(path)
    elif (type(parent_fd) is not int or parent_fd < 0 or
            not isinstance(parent_identity, tuple) or
            len(parent_identity) != 2):
        fail("completion-marker read guard is invalid")
    assert parent_fd is not None
    assert parent_identity is not None
    try:
        _require_completion_parent_descriptor(parent_fd, parent_identity)
        observed, fingerprint = _read_completed_regular_bytes(
            Path(path.name), "published completion marker", parent_fd,
            COMPLETED_ARTIFACT_BYTE_LIMITS["run-summary.json"])
        if ((fingerprint[0], fingerprint[1]) != identity or
                observed != expected_bytes):
            fail("published completion marker changed before commit")
        _verify_completed_directory_path(path.parent, parent_identity)
    finally:
        if owns_parent:
            try:
                os.close(parent_fd)
            except BaseException:
                pass


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
    try:
        for cpu in cpus:
            try:
                process = subprocess.Popen(
                    [description["resolved_path"], "--worker", str(cpu)],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, bufsize=0)
            except OSError as exc:
                fail("cannot spawn native worker on CPU {}: {}".format(cpu, exc))
            workers.append(PersistentWorker(cpu, process, 0))
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
        try:
            terminate_workers(workers)
        except RunnerError as cleanup:
            raise RunnerError(
                "{}; worker cleanup also failed: {}".format(
                    primary, cleanup)) from primary
        raise


def terminate_workers(workers: Sequence[PersistentWorker]) -> None:
    """Exhaust terminate/kill/reap/close and reject every surviving child."""
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.terminate()
            except OSError:
                pass
    limit = time.monotonic() + 2.0
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.wait(timeout=max(0.0, limit - time.monotonic()))
            except subprocess.TimeoutExpired:
                try:
                    worker.process.kill()
                except OSError:
                    pass
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                pass
        for stream in (worker.process.stdin, worker.process.stdout,
                       worker.process.stderr):
            if stream is not None:
                try:
                    stream.close()
                except OSError:
                    pass
    survivors = [
        worker.cpu for worker in workers if worker.process.poll() is None
    ]
    if survivors:
        fail("native worker cleanup left unreaped CPUs {}".format(
            ",".join(str(cpu) for cpu in survivors)))


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
        validator: ResponseValidator) -> Tuple[int, Set[int]]:
    """Run one barrier-delimited batch with at most one job per worker."""
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
                if len(worker.buffer) > MAX_RESPONSE_BYTES:
                    fail("native worker response exceeds the bounded line size")
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
            end, used = run_job_batch(
                workers, jobs, rotation, timing_sink,
                deadline, timing_validator)
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


def _read_completed_regular_bytes(
        path: Path, context: str, directory_fd: int,
        byte_limit: int) -> Tuple[bytes, CompletedFingerprint]:
    """Read one bounded regular file through an already-pinned directory."""
    if type(byte_limit) is not int or byte_limit <= 0:
        fail("{} has an invalid completed-artifact byte cap".format(context))
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(
                context))
        descriptor = os.open(
            str(path), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | getattr(os, "O_NONBLOCK", 0), dir_fd=directory_fd)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} must be a regular non-symlink file".format(context))
        if before.st_size > byte_limit:
            fail("{} exceeds the bounded artifact size".format(context))
        chunks = []
        size = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > byte_limit:
                fail("{} exceeds the bounded artifact size".format(context))
            chunks.append(block)
        after = os.fstat(descriptor)
        stable_before = _completed_fingerprint(before)
        stable_after = _completed_fingerprint(after)
        named = os.stat(
            str(path), dir_fd=directory_fd, follow_symlinks=False)
        stable_named = _completed_fingerprint(named)
        data = b"".join(chunks)
        if (not stat.S_ISREG(named.st_mode) or
                stable_before != stable_after or
                stable_named != stable_after or
                len(data) != before.st_size):
            fail("{} changed while it was being read".format(context))
        return data, stable_after
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return b"", (0, 0, 0, 0, 0, 0, 0)


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


def _require_completed_dependency_bounds(
        output_dir: Path, prospective_summary: bytes) -> None:
    """Apply the strict loader's exact per-file and bundle caps pre-commit."""
    total_bytes = len(prospective_summary)
    if total_bytes > COMPLETED_ARTIFACT_BYTE_LIMITS["run-summary.json"]:
        fail("prospective completion marker exceeds its loader byte cap")
    for name in COMPLETED_SOURCE_NAMES.values():
        if name == "run-summary.json":
            continue
        path = output_dir / name
        try:
            info = os.stat(str(path), follow_symlinks=False)
        except OSError as exc:
            fail("cannot bound completed dependency {}: {}".format(path, exc))
        limit = COMPLETED_ARTIFACT_BYTE_LIMITS[name]
        if not stat.S_ISREG(info.st_mode) or info.st_size > limit:
            fail("completed dependency {} exceeds its producer byte cap".format(
                path))
        total_bytes += info.st_size
        if total_bytes > MAX_COMPLETED_BUNDLE_BYTES:
            fail("completed timing dependencies exceed the bundle byte cap")


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
    """Publish, strictly reopen, and terminally pin the completion marker."""
    summary_path = output_dir / "run-summary.json"
    summary_bytes = (
        contract_api.canonical_json(summary) + "\n").encode("utf-8")
    _require_completed_dependency_bounds(output_dir, summary_bytes)
    parent_fd = -1
    parent_identity: Optional[Tuple[int, int]] = None
    published_identity: Optional[Tuple[int, int]] = None
    committed = False
    try:
        parent_fd, parent_identity = _open_completion_parent(summary_path)
        published_identity = _publish_completion_marker(
            summary_path, summary, parent_fd=parent_fd,
            parent_identity=parent_identity)
        trusted_snapshots, trusted_fingerprints = _read_completed_bundle(
            COMPLETED_SOURCE_NAMES, parent_fd)
        reopened = load_completed_timing_screen(contract, output_dir)
        if (not isinstance(reopened, dict) or set(reopened) != {
                "directory", "directory_identity", "run_summary", "freeze",
                "summary", "execution_receipt", "timing_qualification",
                } or
                reopened["directory_identity"] != parent_identity or
                contract_api.canonical_json(reopened["run_summary"]) !=
                contract_api.canonical_json(summary) or
                contract_api.canonical_json(reopened["freeze"]) !=
                contract_api.canonical_json(freeze) or
                contract_api.canonical_json(reopened["summary"]) !=
                contract_api.canonical_json(timing_summary) or
                contract_api.canonical_json(reopened["execution_receipt"]) !=
                contract_api.canonical_json(execution_receipt) or
                reopened["timing_qualification"].map_sha256 !=
                timing_qualification.map_sha256 or
                reopened["timing_qualification"].qualified_domain_sha256 !=
                timing_qualification.qualified_domain_sha256):
            fail("post-publication timing reopen differs from terminal evidence")
        _require_completion_marker(
            summary_path, published_identity, summary_bytes,
            parent_fd=parent_fd, parent_identity=parent_identity)
        final_snapshots, final_fingerprints = _read_completed_bundle(
            COMPLETED_SOURCE_NAMES, parent_fd)
        if (final_snapshots != trusted_snapshots or
                final_fingerprints != trusted_fingerprints):
            fail("completed timing evidence changed after strict reopen")
        # The path must still name the exact parent inode held for publication.
        # Keep the descriptor live until no fallible commit work remains.  The
        # hard-wall teardown is inside this rollback boundary: afterward no WH2
        # SIGALRM or context-manager teardown can report failure with a marker.
        _verify_completed_directory_path(output_dir, parent_identity)
        if finish_hard_wall is not None:
            finish_hard_wall()
        committed = True
    except BaseException as primary:
        if published_identity is not None:
            try:
                _rollback_completion_marker(
                    summary_path, published_identity, parent_fd=parent_fd,
                    parent_identity=parent_identity)
            except RunnerError as rollback:
                raise RunnerError(
                    "{}; completion-marker rollback also failed: {}".format(
                        primary, rollback)) from primary
        raise
    finally:
        if parent_fd >= 0:
            try:
                os.close(parent_fd)
            except BaseException:
                # A successful close is not part of the on-disk transaction.
                # On success every durability/path check already passed; on
                # failure the exact published inode has already been rolled back.
                pass
    if not committed:
        raise AssertionError("unreachable")


def _write_completed_snapshot(path: Path, data: bytes) -> None:
    try:
        with path.open("xb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except OSError as exc:
        fail("cannot write private timing snapshot {}: {}".format(path, exc))


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
        if directory_fd >= 0:
            os.close(directory_fd)
            directory_fd = -1
        fail("cannot open timing screen directory {}: {}".format(
            timing_screen_dir, exc))
    except BaseException:
        if directory_fd >= 0:
            os.close(directory_fd)
            directory_fd = -1
        raise
    source_names = COMPLETED_SOURCE_NAMES
    try:
        snapshots, snapshot_identities = _read_completed_bundle(
            source_names, directory_fd)
    finally:
        if directory_fd >= 0:
            os.close(directory_fd)

    with tempfile.TemporaryDirectory(prefix="wh2-timing-snapshot-") as raw:
        snapshot_root = Path(raw)
        paths = {key: snapshot_root / name
                 for key, name in source_names.items()}
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

        # Semantic validation may be expensive.  Reopen the exact directory
        # inode, reread every dependency through that pinned descriptor, and
        # require both bytes and directory-entry fingerprints to match the
        # initial snapshot before returning an authoritative bundle.
        terminal_directory_fd = -1
        expected_directory_identity = (
            opened_info.st_dev, opened_info.st_ino)
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
            if terminal_directory_fd >= 0:
                os.close(terminal_directory_fd)
        if (terminal_snapshots != snapshots or
                terminal_identities != snapshot_identities):
            fail("completed timing evidence changed during semantic validation")
        return {
            "directory": str(resolved),
            "directory_identity": expected_directory_identity,
            "run_summary": run_summary,
            "freeze": freeze,
            "summary": timing_summary,
            "execution_receipt": execution_receipt,
            "timing_qualification": qualification,
        }


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
        if qualification_workers and not qualification_clean_shutdown:
            primary = sys.exc_info()[1]
            try:
                terminate_workers(qualification_workers)
            except RunnerError as cleanup:
                if primary is not None:
                    raise RunnerError(
                        "{}; qualification worker cleanup also failed: {}".
                        format(primary, cleanup)) from primary
                raise
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
    try:
        # The validated timing freeze exists before the first T command.
        workers = spawn_workers(description, cpus, deadline)
        controller_pinned = True
        _pin_controller(controller_cpu)
        timing_window_start_ns = choose_new_sampler_start(
            args.sampler_csv, deadline)
        native_path, maximum_worker_end = _run_timing_jobs(
            contract, freeze, timing_qualification, description, workers,
            output_dir,
            timing_window_start_ns, deadline)
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
                    qualification_execution_receipt_sha256)
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
        primary = sys.exc_info()[1]
        cleanup_failures: List[str] = []
        if workers and not clean_shutdown:
            try:
                terminate_workers(workers)
            except RunnerError as cleanup:
                cleanup_failures.append(
                    "timing worker cleanup failed: {}".format(cleanup))
        if controller_pinned:
            try:
                _restore_controller_affinity(original_affinity)
            except RunnerError as cleanup:
                cleanup_failures.append(
                    "controller affinity cleanup failed: {}".format(cleanup))
        if cleanup_failures:
            cleanup_message = "; ".join(cleanup_failures)
            if primary is not None:
                raise RunnerError("{}; {}".format(
                    primary, cleanup_message)) from primary
            fail(cleanup_message)


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
