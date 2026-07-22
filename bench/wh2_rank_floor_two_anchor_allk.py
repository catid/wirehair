#!/usr/bin/env python3
"""Freeze and run an independent production-H12 two-anchor confirmation.

The campaign covers every K in 2..64000 under three fresh seeds and the
burst, adversarial, and repair-only 50% loss schedules.  It compares fixed
D12 with adaptive fixed-D12/two-anchor and adaptive D13 using the exact
mixed10-GF(256)/2-GF(65536), period-244, mix2 production geometry at 64-byte
payloads.  The published R2 group ledger is used only as a full-domain batching
map: its historical salt column is recorded but never applied.  Every cell
uses the named production profile's packet-peel seed xor of zero.

Prepare must run from a clean immutable commit.  Run must use the frozen
script, helper, binary, group ledger, and contract produced by prepare.
Saturated-run timings are provenance only and are never a speed claim.
"""

from __future__ import annotations

import sys

# This module is imported from exact-inventory frozen campaign directories.
sys.dont_write_bytecode = True

import argparse
from array import array
from collections import Counter, defaultdict
import concurrent.futures
from contextlib import contextmanager
import ctypes
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, getcontext
import errno
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import queue
import re
import secrets
import select
import shutil
import signal
import stat
import subprocess
import tempfile
import threading
import time
from typing import Any, Callable, Iterable, Iterator, Sequence

from wh2_rank_floor_two_anchor_screen import (
    THERMAL_FIELDS,
    THERMAL_ROW_MAX_BYTES,
    TIMING_FIELDS,
    binomial_one_sided,
    canonical_json,
    deterministic_equal,
    open_thermal_log,
    parse_bench_output,
    parse_thermal_sample,
    sha256_file as _screen_sha256_file,
    thermal_finish,
    thermal_start,
    validate_thermal_log_descriptor,
    validate_thermal_current,
    validate_thermal_health,
)


GROUPS_SHA256 = (
    "940db61d44fc03462f583c7e2e48bea5c93a830918c8f8e5bf3ba502e840cff4"
)
K_MIN = 2
K_MAX = 64000
K_COUNT = K_MAX - K_MIN + 1
CUTOFF = 4096
SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
ARMS = ("d12", "two_anchor_adaptive", "d13_adaptive")
THERMAL_BASELINE_FIELDS = {
    "dev", "ino", "offset", "edac_ce", "edac_ue", "monotonic_s",
    "max_temperature_c", "row_sha256",
}
BASE_OPTIONS = (
    "precodefail",
    "--bb-list", "64",
    "--overhead", "0",
    "--trials", "1",
    "--threads", "1",
    "--loss", "0.50",
    "--completion", "mixed",
    "--mix-count", "2",
)
JOB_OUTPUT_MAX_BYTES = 64 * 1024 * 1024
PROCESS_POLL_SECONDS = 0.01
PAIRED_HEADER = (
    "K", "group", "group_ledger_salt_unused",
    "active_packet_peel_seed_xor", "schedule", "seed_index", "seed",
    *(
        f"{arm}_{field}"
        for arm in ARMS
        for field in (
            "rank_fail", "error", "heavy_shortfall", "seed_attempt",
            "inact_milli", "binary_def_milli", "heavy_gain_milli",
            "block_xors_milli", "block_muladds_milli",
        )
    ),
)


class CampaignError(RuntimeError):
    pass


class BoundedProcessTimeout(CampaignError):
    pass


class ProcessGroupCleanupError(CampaignError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def sha256_file(path: Path, *, require_unique: bool = True) -> str:
    """Hash one stable, unique regular file using the shared reader."""
    try:
        return _screen_sha256_file(path, require_unique=require_unique)
    except (OSError, ValueError) as exc:
        try:
            metadata = path.lstat()
            recoverable = require_unique and \
                stat.S_ISREG(metadata.st_mode) and \
                metadata.st_nlink == 2
        except OSError:
            recoverable = False
        if recoverable and _reconcile_stale_atomic_publication(path):
            try:
                return _screen_sha256_file(
                    path, require_unique=require_unique)
            except (OSError, ValueError) as retry_exc:
                die(str(retry_exc))
        die(str(exc))


def _canonical_digest(value: Any) -> bool:
    return (
        isinstance(value, str) and len(value) == 64 and
        all(character in "0123456789abcdef" for character in value)
    )


def frozen_executable_path(
    record: Any,
    context: str,
) -> Path:
    """Validate one frozen absolute executable while permitting hardlinks."""
    if (not isinstance(record, dict) or
            set(record) != {"path", "sha256"} or
            not isinstance(record.get("path"), str) or
            not _canonical_digest(record.get("sha256"))):
        die(f"frozen {context} identity record is malformed")
    path = Path(record["path"])
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        die(f"frozen {context} is unavailable: {exc}")
    if (not path.is_absolute() or resolved != path or
            not os.access(path, os.X_OK) or
            sha256_file(path, require_unique=False) != record["sha256"]):
        die(f"frozen {context} changed or is indirect/nonexecutable")
    return path


def executable_identity(name: str) -> dict[str, str]:
    located = shutil.which(name)
    if located is None:
        die(f"required executable is unavailable: {name}")
    try:
        path = Path(located).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        die(f"required executable is unavailable: {name}: {exc}")
    if not path.is_absolute() or not os.access(path, os.X_OK):
        die(f"required executable is not executable: {name}")
    return {
        "path": str(path),
        "sha256": sha256_file(path, require_unique=False),
    }


def python_runtime_identity() -> dict[str, str]:
    try:
        path = Path(sys.executable).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        die(f"Python runtime is unavailable: {exc}")
    if not path.is_absolute() or not os.access(path, os.X_OK):
        die("Python runtime is not an absolute executable")
    return {
        "path": str(path), "version": sys.version,
        "sha256": sha256_file(path, require_unique=False),
    }


def verify_frozen_binary(binary: Path, digest: Any) -> None:
    if (not _canonical_digest(digest) or
            not os.access(binary, os.X_OK) or
            sha256_file(binary) != digest):
        die("frozen benchmark binary changed or is nonexecutable")


def _process_start_ticks(pid: int) -> int | None:
    """Read Linux /proc starttime so PID reuse cannot pin stale markers."""
    try:
        raw = Path(f"/proc/{pid}/stat").read_bytes()
    except (OSError, ValueError):
        return None
    closing = raw.rfind(b") ")
    if closing < 0:
        return None
    fields = raw[closing + 2:].split()
    # The suffix begins at field 3 (state); starttime is field 22.
    if len(fields) <= 19:
        return None
    try:
        value = int(fields[19], 10)
    except ValueError:
        return None
    return value if value >= 0 else None


def _atomic_writer_pid_is_live(
    pid: int,
    start_ticks: int | None = None,
) -> bool:
    if start_ticks is not None:
        return _process_start_ticks(pid) == start_ticks
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def open_durable_directory(
        path: Path, *, create: bool = False,
        reprove_existing: bool = False) -> int:
    """Open an absolute directory without following any component symlink."""
    if reprove_existing and not create:
        die("cannot re-prove directory entries without create=True")
    absolute = path.absolute()
    if (not absolute.is_absolute() or
            any(part in ("", ".", "..") for part in absolute.parts[1:])):
        die(f"noncanonical directory path: {path}")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_DIRECTORY", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    descriptor = os.open("/", flags)
    try:
        for part in absolute.parts[1:]:
            next_descriptor: int | None = None
            try:
                created = False
                if create:
                    try:
                        os.mkdir(part, 0o700, dir_fd=descriptor)
                        created = True
                    except FileExistsError:
                        pass
                next_descriptor = os.open(part, flags, dir_fd=descriptor)
                if created or reprove_existing:
                    # A component can be residue from an earlier call whose
                    # fsync failed after mkdir.  Result-staging setup explicitly
                    # re-proves every traversed edge on retry; ordinary atomic
                    # writes only flush entries created by this invocation.
                    os.fsync(next_descriptor)
                    os.fsync(descriptor)
            except BaseException:
                if next_descriptor is not None:
                    try:
                        os.close(next_descriptor)
                    except OSError:
                        pass
                raise
            previous_descriptor = descriptor
            descriptor = next_descriptor
            os.close(previous_descriptor)
        metadata = os.fstat(descriptor)
        named = os.stat(absolute, follow_symlinks=False)
        if (not stat.S_ISDIR(metadata.st_mode) or
                (metadata.st_dev, metadata.st_ino) !=
                (named.st_dev, named.st_ino)):
            die(f"directory pathname changed while opening: {path}")
        return descriptor
    except OSError as error:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise CampaignError(
            f"cannot open durable directory {path}: {error}") from error
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def create_durable_directory(path: Path) -> int:
    """Exclusively create, fsync, and return one plain final directory."""
    if path.name in ("", ".", ".."):
        die(f"noncanonical directory creation path: {path}")
    parent_fd = open_durable_directory(path.parent, create=True)
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_DIRECTORY", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    descriptor: int | None = None
    try:
        try:
            os.mkdir(path.name, 0o700, dir_fd=parent_fd)
        except FileExistsError as error:
            raise CampaignError(
                f"directory already exists: {path}") from error
        descriptor = os.open(path.name, flags, dir_fd=parent_fd)
        metadata = os.fstat(descriptor)
        named = os.stat(
            path.name, dir_fd=parent_fd, follow_symlinks=False)
        if (not stat.S_ISDIR(metadata.st_mode) or
                (metadata.st_dev, metadata.st_ino) !=
                (named.st_dev, named.st_ino)):
            die(f"directory pathname changed during creation: {path}")
        os.fsync(descriptor)
        os.fsync(parent_fd)
    except BaseException:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except BaseException:
                pass
        try:
            os.close(parent_fd)
        except BaseException:
            pass
        raise
    try:
        os.close(parent_fd)
    except BaseException:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except BaseException:
                pass
        raise
    if descriptor is None:
        die(f"created directory descriptor is unavailable: {path}")
    return descriptor


def _stable_bytes_at(directory_fd: int, name: str, display: Path) -> bytes:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_NONBLOCK", 0))
    try:
        descriptor = os.open(name, flags, dir_fd=directory_fd)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            die(f"refusing symlink input {display}")
        die(f"unable to open unique regular input {display}: {exc}")
    try:
        before = os.fstat(descriptor)
        if stat.S_ISREG(before.st_mode) and before.st_nlink == 2:
            named_before = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
            if ((named_before.st_dev, named_before.st_ino) !=
                    (before.st_dev, before.st_ino)):
                die(f"input pathname changed while being read: {display}")
            if _reconcile_stale_atomic_publication(display, directory_fd):
                before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            die(f"refusing nonunique regular input {display}")
        with os.fdopen(descriptor, "rb", closefd=False) as source:
            data = source.read()
            midpoint = os.fstat(descriptor)
            named_midpoint = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
            source.seek(0)
            confirmation = source.read()
        after = os.fstat(descriptor)
        named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    finally:
        os.close(descriptor)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if (identity(before) != identity(midpoint) or
            len(data) != before.st_size or
            not stat.S_ISREG(named_midpoint.st_mode) or
            named_midpoint.st_nlink != 1 or
            identity(named_midpoint) != identity(midpoint) or
            data != confirmation or
            identity(midpoint) != identity(after) or
            not stat.S_ISREG(named.st_mode) or named.st_nlink != 1 or
            identity(named) != identity(after)):
        die(f"input changed while being read: {display}")
    return data


def _entry_exists_at(directory_fd: int, name: str) -> bool:
    try:
        os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return True


def _discard_stale_atomic_partials(
    path: Path,
    directory_fd: int | None = None,
) -> bool:
    """Remove only exact dead-PID temporaries for one locked target."""
    suffix = ".partial"
    removed = False
    own_descriptor = directory_fd is None
    if directory_fd is None:
        directory_fd = open_durable_directory(path.parent)
    try:
        names = os.listdir(directory_fd)
        for name in names:
            candidate = path.parent / name
            if not name.startswith(".") or not name.endswith(suffix):
                continue
            body = name[1:-len(suffix)]
            prefix = path.name + "."
            if not body.startswith(prefix):
                # Atomic names for dotted siblings can share this target's
                # textual prefix (for example record.json and its
                # record.json.sha256 sidecar).  Parse from the right so each
                # caller touches only its exact target's temporary.
                continue
            identity_fields = body[len(prefix):].split(".")
            if len(identity_fields) == 1:
                pid_text = identity_fields[0]
                start_text = None
            elif (len(identity_fields) == 2 and
                    all(field.isdigit() for field in identity_fields)):
                pid_text, start_text = identity_fields
            else:
                # This is a marker for a dotted sibling target, not this
                # target (for example record.json.sha256.<pid>.partial).
                continue
            if (not pid_text.isdigit() or pid_text == "0" or
                    str(int(pid_text, 10)) != pid_text or
                    (start_text is not None and (
                        start_text == "" or str(int(start_text, 10)) !=
                        start_text))):
                die(f"malformed atomic temporary output: {candidate}")
            pid = int(pid_text, 10)
            start_ticks = (
                int(start_text, 10) if start_text is not None else None)
            if _atomic_writer_pid_is_live(pid, start_ticks):
                die(f"atomic output writer is still active: {candidate}")
            metadata = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
            if (stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 2 and
                    _entry_exists_at(directory_fd, path.name)):
                target_metadata = os.stat(
                    path.name, dir_fd=directory_fd, follow_symlinks=False)
                if (stat.S_ISREG(target_metadata.st_mode) and
                        target_metadata.st_nlink == 2 and
                        (target_metadata.st_dev, target_metadata.st_ino) ==
                        (metadata.st_dev, metadata.st_ino)):
                    # First make the published target name durable while its
                    # second-name marker still exists, then retire the marker.
                    os.fsync(directory_fd)
                    os.unlink(name, dir_fd=directory_fd)
                    removed = True
                    continue
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                die(f"refusing nonunique atomic temporary output: {candidate}")
            os.unlink(name, dir_fd=directory_fd)
            removed = True
    finally:
        if own_descriptor:
            os.close(directory_fd)
    return removed


def _reconcile_stale_atomic_publication(
    path: Path,
    directory_fd: int | None = None,
) -> bool:
    """Retire an exact dead-writer hardlink marker for a published target."""
    own_descriptor = directory_fd is None
    if directory_fd is None:
        directory_fd = open_durable_directory(path.parent)
    try:
        before = os.stat(
            path.name, dir_fd=directory_fd, follow_symlinks=False)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 2:
            return False
        identity = (before.st_dev, before.st_ino)
        if not _discard_stale_atomic_partials(path, directory_fd):
            return False
        os.fsync(directory_fd)
        after = os.stat(
            path.name, dir_fd=directory_fd, follow_symlinks=False)
        if (not stat.S_ISREG(after.st_mode) or after.st_nlink != 1 or
                (after.st_dev, after.st_ino) != identity):
            die(f"refusing nonunique regular input {path}")
        return True
    finally:
        if own_descriptor:
            os.close(directory_fd)


def _fsync_directory(path: Path) -> None:
    descriptor = open_durable_directory(path)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write(path: Path, data: bytes) -> None:
    directory_fd = open_durable_directory(path.parent, create=True)
    writer_pid = os.getpid()
    writer_start_ticks = _process_start_ticks(writer_pid)
    if writer_start_ticks is None:
        os.close(directory_fd)
        die("atomic output writer identity is unavailable")
    temporary_name = (
        f".{path.name}.{writer_pid}.{writer_start_ticks}.partial")
    old_mask = None
    temporary_created = False
    try:
        if hasattr(signal, "pthread_sigmask"):
            old_mask = signal.pthread_sigmask(
                signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
        if _discard_stale_atomic_partials(path, directory_fd):
            os.fsync(directory_fd)
        flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                 getattr(os, "O_NOFOLLOW", 0))
        try:
            descriptor = os.open(
                temporary_name, flags, 0o600, dir_fd=directory_fd)
        except FileExistsError:
            die(f"temporary output already exists: {path.parent / temporary_name}")
        temporary_created = True
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        os.replace(
            temporary_name, path.name,
            src_dir_fd=directory_fd, dst_dir_fd=directory_fd)
        temporary_created = False
        os.fsync(directory_fd)
    finally:
        try:
            if temporary_created and _entry_exists_at(
                    directory_fd, temporary_name):
                os.unlink(temporary_name, dir_fd=directory_fd)
                os.fsync(directory_fd)
        finally:
            # Cleanup failures must never strand SIGINT/SIGTERM blocked in
            # this thread.  In particular, unlink/fsync can fail after the
            # write path has already raised.
            try:
                if old_mask is not None:
                    signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
            finally:
                os.close(directory_fd)


def atomic_write_once_or_same(path: Path, data: bytes) -> None:
    """Durably publish immutable bytes without ever replacing the target."""
    directory_fd = open_durable_directory(path.parent, create=True)
    writer_pid = os.getpid()
    writer_start_ticks = _process_start_ticks(writer_pid)
    if writer_start_ticks is None:
        os.close(directory_fd)
        die("atomic output writer identity is unavailable")
    temporary_name = (
        f".{path.name}.{writer_pid}.{writer_start_ticks}.partial")
    old_mask = None
    temporary_created = False
    target_linked = False
    try:
        if hasattr(signal, "pthread_sigmask"):
            old_mask = signal.pthread_sigmask(
                signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
        if _discard_stale_atomic_partials(path, directory_fd):
            os.fsync(directory_fd)
        if _entry_exists_at(directory_fd, path.name):
            if _stable_bytes_at(directory_fd, path.name, path) != data:
                die(f"immutable output already differs: {path}")
            os.fsync(directory_fd)
            return
        flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                 getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NOFOLLOW", 0))
        descriptor = os.open(
            temporary_name, flags, 0o600, dir_fd=directory_fd)
        temporary_created = True
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        try:
            os.link(
                temporary_name, path.name,
                src_dir_fd=directory_fd, dst_dir_fd=directory_fd,
                follow_symlinks=False)
        except FileExistsError:
            if _stable_bytes_at(directory_fd, path.name, path) != data:
                die(f"immutable output won a conflicting race: {path}")
        else:
            target_linked = True
            # Persist the target link while the temporary remains a durable
            # recovery marker, then persist removal of that second name.
            os.fsync(directory_fd)
            os.unlink(temporary_name, dir_fd=directory_fd)
            temporary_created = False
            target_linked = False
            os.fsync(directory_fd)
        if _stable_bytes_at(directory_fd, path.name, path) != data:
            die(f"immutable output publication was not stable: {path}")
    finally:
        try:
            if (temporary_created and not target_linked and _entry_exists_at(
                    directory_fd, temporary_name)):
                os.unlink(temporary_name, dir_fd=directory_fd)
                os.fsync(directory_fd)
        finally:
            try:
                if old_mask is not None:
                    signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
            finally:
                os.close(directory_fd)


def json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value, indent=2, sort_keys=True, allow_nan=False,
            ensure_ascii=True,
        ) + "\n"
    ).encode("utf-8")


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(path, json_bytes(value))


def stable_bytes(path: Path) -> bytes:
    directory_fd = open_durable_directory(path.parent)
    try:
        return _stable_bytes_at(directory_fd, path.name, path)
    finally:
        os.close(directory_fd)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_prepare_anchor(
    result_dir: Path,
    schema: str,
    fields: Sequence[str],
    staged_field: str,
) -> dict[str, Any]:
    """Bind the mutable staged manifest to its outer prepare record."""
    try:
        record = json.loads(stable_bytes(result_dir / "prepare.json"))
    except json.JSONDecodeError as exc:
        raise CampaignError("prepare anchor is not canonical JSON") from exc
    expected_fields = {"schema", "result_dir", *fields}
    if (not isinstance(record, dict) or set(record) != expected_fields or
            record.get("schema") != schema or
            record.get("result_dir") != str(result_dir) or
            not isinstance(record.get(staged_field), str) or
            not _canonical_digest(record.get(staged_field)) or
            sha256_file(result_dir / "frozen/staged.sha256") !=
            record[staged_field]):
        die("prepare anchor does not bind the frozen staged manifest")
    return record


def strict_uint(text: str, context: str) -> int:
    if not text or (text != "0" and text.startswith("0")) or not text.isdigit():
        die(f"{context}: noncanonical unsigned integer {text!r}")
    return int(text, 10)


def canonical_u64(text: str, context: str) -> int:
    try:
        value = int(text, 0)
    except ValueError:
        die(f"{context}: invalid integer {text!r}")
    if value < 0 or value >= 1 << 64:
        die(f"{context}: integer does not fit u64")
    return value


@dataclass(frozen=True)
class Group:
    group: int
    ledger_salt: str
    ks: tuple[int, ...]


@dataclass(frozen=True)
class Job:
    job: int
    arm: str
    band: str
    seed_index: int
    seed: str
    schedule: str
    group: int
    ledger_salt: str
    ks: tuple[int, ...]

    @property
    def stem(self) -> str:
        return (
            f"job{self.job:05d}.{self.arm}.seed{self.seed_index}."
            f"{self.schedule}.group{self.group:03d}.{self.band}"
        )


def load_groups(path: Path, expected_sha256: str = GROUPS_SHA256) -> list[Group]:
    data = stable_bytes(path)
    if sha256_bytes(data) != expected_sha256:
        die(f"group ledger SHA256 mismatch: {path}")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        die(f"group ledger is not UTF-8: {exc}")
    groups: list[Group] = []
    seen: set[int] = set()
    for number, row in enumerate(csv.reader(text.splitlines(), delimiter="\t"), 1):
        if len(row) != 3:
            die(f"groups:{number}: expected three fields")
        group = strict_uint(row[0], f"groups:{number}:group")
        salt_value = canonical_u64(row[1], f"groups:{number}:salt")
        ledger_salt = f"0x{salt_value:x}"
        values = tuple(
            strict_uint(value, f"groups:{number}:K")
            for value in row[2].split(",")
        )
        if not values or tuple(sorted(values)) != values or len(set(values)) != len(values):
            die(f"groups:{number}: K list is empty, unsorted, or duplicated")
        if any(K < K_MIN or K > K_MAX or K in seen for K in values):
            die(f"groups:{number}: K is outside domain or repeated")
        if group != len(groups):
            die(f"groups:{number}: group IDs are not contiguous from zero")
        seen.update(values)
        groups.append(Group(group, ledger_salt, values))
    if len(groups) != 186 or seen != set(range(K_MIN, K_MAX + 1)):
        die("group ledger is not 186 groups covering exact K=2..64000")
    return groups


def arm_options(arm: str, band: str) -> tuple[str, ...]:
    if arm == "d12":
        return ()
    if arm == "two_anchor_adaptive":
        return ("--binary-dense-two-anchor",) if band == "large" else ()
    if arm == "d13_adaptive":
        return ("--binary-dense-rows", "13") if band == "large" else ()
    die(f"unknown arm {arm}")
    raise AssertionError


def build_jobs(groups: Sequence[Group]) -> list[Job]:
    jobs: list[Job] = []
    cells = 0
    for arm in ARMS:
        for seed_index, seed in enumerate(SEEDS):
            for schedule in SCHEDULES:
                for group in groups:
                    for band, ks in (
                        ("small", tuple(K for K in group.ks if K < CUTOFF)),
                        ("large", tuple(K for K in group.ks if K >= CUTOFF)),
                    ):
                        if not ks:
                            continue
                        jobs.append(Job(
                            len(jobs), arm, band, seed_index, seed, schedule,
                            group.group, group.ledger_salt, ks,
                        ))
                        cells += len(ks)
    expected_cells = len(ARMS) * len(SEEDS) * len(SCHEDULES) * K_COUNT
    if cells != expected_cells:
        die(f"job ledger covers {cells} cells, want {expected_cells}")
    keys = {
        (job.arm, job.seed_index, job.schedule, K)
        for job in jobs for K in job.ks
    }
    if len(keys) != expected_cells:
        die("job ledger has duplicate or missing arm/stratum/K cells")
    return jobs


def make_command(binary: Path, job: Job) -> list[str]:
    return [
        str(binary), *BASE_OPTIONS,
        "--N", ",".join(map(str, job.ks)),
        "--seed", job.seed,
        "--schedule", job.schedule,
        "--packet-peel-seed-xor", "0",
        *arm_options(job.arm, job.band),
    ]


class CpuPool:
    def __init__(self, cpus: Sequence[int], workers: int) -> None:
        self._queue: queue.Queue[int] = queue.Queue()
        for cpu in list(cpus)[:workers]:
            self._queue.put(cpu)

    def acquire(self) -> int:
        return self._queue.get()

    def release(self, cpu: int) -> None:
        self._queue.put(cpu)


def signal_process_group(process: subprocess.Popen[bytes], signum: int) -> None:
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        pass


def process_group_exists(process: subprocess.Popen[bytes]) -> bool:
    try:
        os.killpg(process.pid, 0)
        return True
    except ProcessLookupError:
        return False


def close_process_streams(process: subprocess.Popen[bytes]) -> None:
    for stream in (process.stdout, process.stderr):
        if stream is not None:
            try:
                stream.close()
            except OSError:
                pass


def bounded_capture_exceeded(
    streams: Sequence[Any], output_limit: int,
) -> bool:
    """Check tempfile-backed captures without reading them into memory."""
    if (not isinstance(output_limit, int) or isinstance(output_limit, bool) or
            output_limit <= 0):
        die("bounded capture limit is invalid")
    try:
        return any(
            os.fstat(stream.fileno()).st_size > output_limit
            for stream in streams
        )
    except (AttributeError, OSError, ValueError) as exc:
        raise CampaignError(f"bounded capture inspection failed: {exc}") \
            from exc


def read_bounded_capture(
    stream: Any,
    output_limit: int,
    context: str,
) -> bytes:
    """Read at most one byte past the cap, after the writer group is gone."""
    if (not isinstance(output_limit, int) or isinstance(output_limit, bool) or
            output_limit <= 0 or not isinstance(context, str) or not context):
        die("bounded capture read arguments are invalid")
    try:
        size_before = os.fstat(stream.fileno()).st_size
        if size_before > output_limit:
            die(f"{context} exceeded the bounded capture limit")
        stream.seek(0)
        data = stream.read(output_limit + 1)
        size_after = os.fstat(stream.fileno()).st_size
    except CampaignError:
        raise
    except (AttributeError, OSError, ValueError) as exc:
        raise CampaignError(f"{context} capture read failed: {exc}") from exc
    if (len(data) > output_limit or size_after > output_limit or
            size_after != size_before or len(data) != size_after):
        die(f"{context} exceeded or changed during bounded capture")
    return data


def stop_and_reap_process_group(
    process: subprocess.Popen[bytes], grace_seconds: float = 5.0,
    *,
    graceful: bool = True,
) -> None:
    if (not math.isfinite(grace_seconds) or grace_seconds < 0 or
            not isinstance(graceful, bool)):
        die("process-group cleanup grace is invalid")
    wait_error: BaseException | None = None
    try:
        if graceful:
            signal_process_group(process, signal.SIGTERM)
        # Never drain pipes while cleaning up.  A TERM-ignoring writer can
        # otherwise make communicate() accumulate without a memory bound.
        term_deadline = time.monotonic() + min(grace_seconds, 0.25)
        while (graceful and process_group_exists(process) and
               time.monotonic() < term_deadline):
            process.poll()
            time.sleep(PROCESS_POLL_SECONDS)
        if process_group_exists(process):
            signal_process_group(process, signal.SIGKILL)
        try:
            process.wait(timeout=max(grace_seconds, 0.25))
        except subprocess.TimeoutExpired:
            pass
        except BaseException as error:
            # Finish the kill/proof sequence before replaying an interruption.
            wait_error = error
        deadline = time.monotonic() + grace_seconds
        while process_group_exists(process) and time.monotonic() < deadline:
            process.poll()
            time.sleep(PROCESS_POLL_SECONDS)
        if process_group_exists(process):
            signal_process_group(process, signal.SIGKILL)
            kill_deadline = time.monotonic() + grace_seconds
            while (process_group_exists(process) and
                    time.monotonic() < kill_deadline):
                process.poll()
                time.sleep(PROCESS_POLL_SECONDS)
    finally:
        close_process_streams(process)
    if process.poll() is None or process_group_exists(process):
        cleanup_error = CampaignError(
            f"child process group {process.pid} was not reaped")
        if wait_error is not None:
            raise cleanup_error from wait_error
        raise cleanup_error
    if wait_error is not None:
        raise wait_error


class ProcessRegistry:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.processes: set[subprocess.Popen[bytes]] = set()

    def add(self, process: subprocess.Popen[bytes]) -> None:
        with self.lock:
            self.processes.add(process)

    def remove(self, process: subprocess.Popen[bytes]) -> None:
        with self.lock:
            self.processes.discard(process)

    def signal_all(self, signum: int) -> None:
        with self.lock:
            processes = tuple(self.processes)
        for process in processes:
            signal_process_group(process, signum)

    def count(self) -> int:
        with self.lock:
            return len(self.processes)

    def snapshot(self) -> tuple[subprocess.Popen[bytes], ...]:
        with self.lock:
            return tuple(self.processes)

    def drain(self, timeout: float = 6.0) -> None:
        self.signal_all(signal.SIGTERM)
        deadline = time.monotonic() + timeout
        while self.count() and time.monotonic() < deadline:
            time.sleep(0.05)
        if self.count():
            for process in self.snapshot():
                try:
                    stop_and_reap_process_group(process, timeout)
                except CampaignError:
                    continue
                self.remove(process)
        if self.count():
            die("campaign cleanup could not drain all child process groups")


class CampaignSignalGuard:
    def __init__(self, abort: threading.Event, registry: ProcessRegistry) -> None:
        self.abort = abort
        self.registry = registry
        self.received: int | None = None
        self.previous: dict[int, Any] = {}
        self.read_fd = -1
        self.write_fd = -1
        self.control_stop = False
        self.control_thread: threading.Thread | None = None
        self.control_started = False

    def _handle(self, signum: int, _frame: Any) -> None:
        if self.received is None:
            self.received = signum
        try:
            os.write(self.write_fd, b"s")
        except (BlockingIOError, OSError):
            pass

    def _control_loop(self) -> None:
        while True:
            try:
                readable, _, _ = select.select([self.read_fd], [], [], 1.0)
            except (OSError, ValueError):
                return
            if not readable:
                if self.control_stop:
                    return
                continue
            try:
                os.read(self.read_fd, 4096)
            except (BlockingIOError, OSError):
                pass
            if self.control_stop:
                return
            self.abort.set()
            self.registry.signal_all(signal.SIGTERM)

    def _stop_control(self) -> None:
        self.control_stop = True
        if self.write_fd >= 0:
            try:
                os.write(self.write_fd, b"x")
            except (BlockingIOError, OSError):
                pass
        if self.control_thread is not None and self.control_started:
            self.control_thread.join()
        for signum, handler in self.previous.items():
            signal.signal(signum, handler)
        for fd in (self.read_fd, self.write_fd):
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass
        self.read_fd = self.write_fd = -1

    def __enter__(self) -> "CampaignSignalGuard":
        if threading.current_thread() is not threading.main_thread():
            die("campaign signal guard must run on the main thread")
        self.read_fd, self.write_fd = os.pipe()
        try:
            os.set_blocking(self.write_fd, False)
            self.control_thread = threading.Thread(
                target=self._control_loop, name="campaign-signal-control"
            )
            try:
                self.control_thread.start()
            finally:
                self.control_started = self.control_thread.ident is not None
            for signum in (signal.SIGINT, signal.SIGTERM):
                self.previous[signum] = signal.getsignal(signum)
                signal.signal(signum, self._handle)
        except BaseException:
            self._stop_control()
            raise
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> bool:
        self.abort.set()
        try:
            self.registry.drain()
        finally:
            self._stop_control()
        if self.received is not None:
            name = signal.Signals(self.received).name
            raise CampaignError(f"campaign terminated by external {name}") from exc
        return False


def validate_frozen_thermal_source(
    path: Path,
    baseline: Any,
    stale_seconds: float,
    cpu_limit_c: float = 90.0,
    dimm_limit_c: float = 90.0,
    consecutive_limit: int = 3,
) -> tuple[int, int, int, int]:
    """Revalidate the frozen row and every append-only telemetry row."""
    plain_nonnegative = lambda value: (
        isinstance(value, int) and not isinstance(value, bool) and value >= 0
    )
    if (not math.isfinite(stale_seconds) or stale_seconds <= 0.0 or
            not math.isfinite(cpu_limit_c) or
            not 0.0 <= cpu_limit_c <= 120.0 or
            not math.isfinite(dimm_limit_c) or
            not 0.0 <= dimm_limit_c <= 120.0 or
            not isinstance(consecutive_limit, int) or
            isinstance(consecutive_limit, bool) or consecutive_limit <= 0 or
            not isinstance(baseline, dict) or
            set(baseline) != THERMAL_BASELINE_FIELDS or
            not plain_nonnegative(baseline.get("dev")) or
            baseline.get("dev") == 0 or
            not plain_nonnegative(baseline.get("ino")) or
            baseline.get("ino") == 0 or
            not plain_nonnegative(baseline.get("offset")) or
            baseline.get("offset") == 0 or
            not plain_nonnegative(baseline.get("edac_ce")) or
            not plain_nonnegative(baseline.get("edac_ue")) or
            any(not isinstance(baseline.get(field), (int, float)) or
                isinstance(baseline.get(field), bool) or
                not math.isfinite(float(baseline[field]))
                for field in ("monotonic_s", "max_temperature_c")) or
            float(baseline["monotonic_s"]) < 0.0 or
            not 0.0 <= float(baseline["max_temperature_c"]) <= 120.0 or
            not _canonical_digest(baseline.get("row_sha256"))):
        die("frozen thermal baseline record is malformed")
    try:
        current = thermal_start(
            path, stale_seconds=stale_seconds, require_zero_edac=False)
    except (OSError, ValueError) as exc:
        die(f"frozen thermal source is unavailable: {exc}")
    expected_identity = (baseline["dev"], baseline["ino"])
    if ((current["dev"], current["ino"]) != expected_identity or
            current["offset"] < baseline["offset"] or
            current["edac_ce"] != baseline["edac_ce"] or
            current["edac_ue"] != baseline["edac_ue"] or
            current["monotonic_s"] < float(baseline["monotonic_s"])):
        die("thermal source changed since the public freeze")
    try:
        with open_thermal_log(path) as source:
            metadata = validate_thermal_log_descriptor(source, path)
            if ((metadata.st_dev, metadata.st_ino) != expected_identity or
                    metadata.st_size < baseline["offset"]):
                die("frozen thermal baseline is no longer retained")
            header = source.readline()
            if header != current["header"]:
                die("frozen thermal header changed")
            window_start = max(
                len(header), baseline["offset"] - THERMAL_ROW_MAX_BYTES)
            source.seek(window_start)
            window = source.read(baseline["offset"] - window_start)
            window_rows = window.splitlines(keepends=True)
            if not window_rows:
                die("frozen thermal baseline row is unavailable")
            frozen_row_start = baseline["offset"] - len(window_rows[-1])
            if frozen_row_start < len(header):
                die("frozen thermal baseline overlaps its header")
            if frozen_row_start > len(header):
                source.seek(frozen_row_start - 1)
                if source.read(1) != b"\n":
                    die("frozen thermal baseline exceeds the row size limit")
            source.seek(baseline["offset"])
            interval = source.read(current["offset"] - baseline["offset"])
            if len(interval) != current["offset"] - baseline["offset"]:
                die("frozen thermal history changed while being read")
            validate_thermal_log_descriptor(source, path)
    except (OSError, ValueError) as exc:
        die(f"frozen thermal baseline cannot be reread: {exc}")
    rows = window.splitlines(keepends=True)
    if not rows or not window.endswith(b"\n"):
        die("frozen thermal baseline row is unavailable")
    row = rows[-1]
    if sha256_bytes(row) != baseline["row_sha256"]:
        die("frozen thermal baseline row changed")
    try:
        sample = parse_thermal_sample(row, "frozen thermal baseline")
    except ValueError as exc:
        die(f"frozen thermal baseline row is invalid: {exc}")
    if (sample["edac_ce"] != baseline["edac_ce"] or
            sample["edac_ue"] != baseline["edac_ue"] or
            sample["monotonic_s"] != float(baseline["monotonic_s"]) or
            sample["max_temperature_c"] !=
            float(baseline["max_temperature_c"])):
        die("frozen thermal baseline record differs from its row")
    consecutive_high = int(
        sample["cpu_tctl_c"] >= cpu_limit_c or
        max(sample["dimm_temperatures"]) >= dimm_limit_c)
    if consecutive_high >= consecutive_limit:
        die("frozen thermal history exceeded its consecutive limit")
    if interval and not interval.endswith(b"\n"):
        die("frozen thermal history contains a partial row")
    previous = sample
    for number, history_row in enumerate(
            interval.splitlines(keepends=True), 1):
        try:
            history = parse_thermal_sample(
                history_row, f"frozen thermal history row {number}")
            validate_thermal_health(
                history, f"frozen thermal history row {number}",
                require_zero_edac=False)
        except ValueError as exc:
            die(f"frozen thermal history is invalid: {exc}")
        delta = history["monotonic_s"] - previous["monotonic_s"]
        if delta <= 0.0 or delta > stale_seconds:
            die("frozen thermal history is nonmonotonic or contains a gap")
        if (history["edac_ce"] != baseline["edac_ce"] or
                history["edac_ue"] != baseline["edac_ue"]):
            die("frozen thermal history changed the EDAC counters")
        high = (history["cpu_tctl_c"] >= cpu_limit_c or
                max(history["dimm_temperatures"]) >= dimm_limit_c)
        consecutive_high = consecutive_high + 1 if high else 0
        if consecutive_high >= consecutive_limit:
            die("frozen thermal history exceeded its consecutive limit")
        previous = history
    if (previous["raw"] != current["baseline_row"] or
            previous["monotonic_s"] != current["monotonic_s"]):
        die("frozen thermal history does not reach the current baseline")
    return (
        baseline["dev"], baseline["ino"],
        baseline["edac_ce"], baseline["edac_ue"])


class ThermalGuard:
    def __init__(
        self, path: Path, abort: threading.Event,
        stale_seconds: float = 5.0, limit_c: float = 90.0,
        consecutive_limit: int = 3,
        expected_edac_ce: int = 0, expected_edac_ue: int = 0,
        expected_dev: int | None = None, expected_ino: int | None = None,
    ) -> None:
        self.path = path
        self.abort = abort
        self.stale_seconds = stale_seconds
        self.limit_c = limit_c
        self.consecutive_limit = consecutive_limit
        if (not math.isfinite(limit_c) or not 0.0 <= limit_c <= 120.0 or
                consecutive_limit <= 0 or
                not isinstance(expected_edac_ce, int) or
                isinstance(expected_edac_ce, bool) or expected_edac_ce < 0 or
                not isinstance(expected_edac_ue, int) or
                isinstance(expected_edac_ue, bool) or expected_edac_ue < 0 or
                (expected_dev is None) != (expected_ino is None) or
                (expected_dev is not None and (
                    not isinstance(expected_dev, int) or
                    isinstance(expected_dev, bool) or expected_dev <= 0 or
                    not isinstance(expected_ino, int) or
                    isinstance(expected_ino, bool) or expected_ino <= 0))):
            die("thermal guard policy is outside the canonical range")
        self.expected_edac_ce = expected_edac_ce
        self.expected_edac_ue = expected_edac_ue
        self.mark = thermal_start(
            path, stale_seconds=stale_seconds, require_zero_edac=False)
        if (expected_dev is not None and
                (self.mark["dev"], self.mark["ino"]) !=
                (expected_dev, expected_ino)):
            die("thermal logger identity changed since the frozen baseline")
        self._validate_sample(self.mark["baseline"], "thermal baseline")
        self.initial_consecutive_high = self._trailing_high_samples()
        if self.initial_consecutive_high >= self.consecutive_limit:
            die(
                f"thermal limit {self.limit_c:g} C was already reached for "
                f"{self.initial_consecutive_high} consecutive samples"
            )
        self.read_offset = self.mark["offset"]
        self.last_monotonic = self.mark["monotonic_s"]
        self.stop_event = threading.Event()
        self.error: str | None = None
        self.poll_iterations = 0
        self.samples = 0
        self.high_samples = 0
        self.thread = threading.Thread(target=self._loop, name="thermal-guard")
        self.started = False

    def _trailing_high_samples(self) -> int:
        """Bridge a high-temperature streak across the live guard mark."""
        with open_thermal_log(self.path) as source:
            metadata = validate_thermal_log_descriptor(source, self.path)
            if ((metadata.st_dev, metadata.st_ino) !=
                    (self.mark["dev"], self.mark["ino"]) or
                    metadata.st_size < self.mark["offset"]):
                die("thermal logger changed while seeding the live guard")
            header = source.readline()
            if header != self.mark["header"]:
                die("thermal logger header changed while seeding the live guard")
            start = max(
                len(header),
                self.mark["offset"] -
                self.consecutive_limit * THERMAL_ROW_MAX_BYTES,
            )
            source.seek(start)
            window = source.read(self.mark["offset"] - start)
            validate_thermal_log_descriptor(source, self.path)
        rows = window.splitlines(keepends=True)
        if (not rows or not window.endswith(b"\n") or
                rows[-1] != self.mark["baseline_row"]):
            die("thermal guard mark is not retained byte-exactly")
        retained = rows[-self.consecutive_limit:]
        previous: dict[str, Any] | None = None
        consecutive = 0
        for number, row in enumerate(retained, 1):
            if len(row) > THERMAL_ROW_MAX_BYTES:
                die("thermal guard history contains an oversized row")
            try:
                sample = parse_thermal_sample(
                    row, f"thermal guard history row {number}")
            except ValueError as exc:
                die(f"thermal guard history is invalid: {exc}")
            self._validate_sample(sample, "thermal guard history")
            if previous is not None:
                delta = sample["monotonic_s"] - previous["monotonic_s"]
                if delta <= 0.0 or delta > self.stale_seconds:
                    die("thermal guard history is nonmonotonic or has a gap")
            previous = sample
            high = sample["max_temperature_c"] >= self.limit_c
            consecutive = consecutive + 1 if high else 0
        return consecutive

    def _validate_sample(self, sample: dict[str, Any], context: str) -> None:
        validate_thermal_health(
            sample, context, require_zero_edac=False)
        if (sample["edac_ce"] != self.expected_edac_ce or
                sample["edac_ue"] != self.expected_edac_ue):
            die(f"{context} changed the frozen EDAC counters")

    def start(self) -> None:
        try:
            self.thread.start()
        finally:
            self.started = self.thread.ident is not None

    def _new_samples(self) -> list[dict[str, Any]]:
        with open_thermal_log(self.path) as source:
            file_stat = validate_thermal_log_descriptor(source, self.path)
            if ((file_stat.st_dev, file_stat.st_ino) !=
                    (self.mark["dev"], self.mark["ino"])):
                die("thermal logger inode changed")
            if file_stat.st_size < self.read_offset:
                die("thermal logger shrank")
            source.seek(self.read_offset)
            suffix = source.read()
            validate_thermal_log_descriptor(source, self.path)
        complete_end = suffix.rfind(b"\n")
        if complete_end < 0:
            if len(suffix) > THERMAL_ROW_MAX_BYTES:
                die("thermal logger has an oversized partial row")
            return []
        complete = suffix[:complete_end + 1]
        partial = suffix[complete_end + 1:]
        if len(partial) > THERMAL_ROW_MAX_BYTES:
            die("thermal logger has an oversized partial row")
        lines = complete.splitlines(keepends=True)
        samples = [
            parse_thermal_sample(line, f"live thermal row {number}")
            for number, line in enumerate(lines, 1)
        ]
        self.read_offset += len(complete)
        return samples

    def _loop(self) -> None:
        consecutive = self.initial_consecutive_high
        while not self.stop_event.wait(0.25):
            try:
                samples = self._new_samples()
                self.poll_iterations += 1
                for sample in samples:
                    delta = sample["monotonic_s"] - self.last_monotonic
                    if delta <= 0.0 or delta > self.stale_seconds:
                        die("thermal logger is nonmonotonic or has a sampling gap")
                    self._validate_sample(sample, "live thermal sample")
                    self.last_monotonic = sample["monotonic_s"]
                    self.samples += 1
                    high = sample["max_temperature_c"] >= self.limit_c
                    self.high_samples += high
                    consecutive = consecutive + 1 if high else 0
                    if consecutive >= self.consecutive_limit:
                        die(
                            f"thermal limit {self.limit_c:g} C reached for "
                            f"{consecutive} consecutive samples"
                        )
                current = samples[-1] if samples else {
                    "monotonic_s": self.last_monotonic
                }
                validate_thermal_current(
                    current, self.stale_seconds, "live thermal tail"
                )
            except Exception as exc:  # Preserve the first guard failure for main.
                self.error = str(exc)
                self.abort.set()
                return

    def finish(
        self,
        output: Path,
        *,
        cover_through_monotonic_ns: int | None = None,
    ) -> dict[str, Any]:
        if not self.started:
            die("thermal guard was never started")
        self.stop_event.set()
        self.thread.join()
        if cover_through_monotonic_ns is not None:
            if (type(cover_through_monotonic_ns) is not int or
                    not 0 < cover_through_monotonic_ns < 1 << 64):
                die("thermal coverage target is malformed")
            target_seconds = cover_through_monotonic_ns / 1_000_000_000
            deadline = time.monotonic() + min(
                max(self.stale_seconds, 1.0), 10.0)
            while True:
                current = thermal_start(
                    self.path, stale_seconds=self.stale_seconds,
                    require_zero_edac=False)
                if ((current["dev"], current["ino"]) !=
                        (self.mark["dev"], self.mark["ino"])):
                    die("thermal logger changed before end coverage")
                self._validate_sample(
                    current["baseline"], "thermal end-coverage sample")
                if current["monotonic_s"] >= target_seconds:
                    break
                if time.monotonic() >= deadline:
                    die("thermal logger did not cover the final job end")
                time.sleep(0.05)
        summary = thermal_finish(
            self.path, self.mark, output,
            limit_c=self.limit_c, stale_seconds=self.stale_seconds,
            require_zero_edac=False, dimm_limit_c=self.limit_c,
        )
        if summary["edac_ce_delta"] != 0 or summary["edac_ue_delta"] != 0:
            die("thermal interval changed the frozen EDAC counters")
        summary.update({
            "baseline_row_sha256": sha256_bytes(self.mark["baseline_row"]),
            "guard_poll_iterations": self.poll_iterations,
            "guard_samples": self.samples,
            "guard_high_samples": self.high_samples,
            "guard_limit_c": self.limit_c,
            "guard_error": self.error,
        })
        return summary


def run_job(
    job: Job,
    binary: Path,
    binary_sha256: str,
    taskset: Path,
    taskset_record: dict[str, str],
    result_root: Path,
    pool: CpuPool,
    abort: threading.Event,
    registry: ProcessRegistry,
    timeout: float,
) -> dict[str, Any]:
    if abort.is_set():
        die(f"campaign abort set before job {job.job}")
    cpu = pool.acquire()
    try:
        if abort.is_set():
            die(f"campaign abort set before job {job.job} launch")
        verify_frozen_binary(binary, binary_sha256)
        if frozen_executable_path(taskset_record, "taskset") != taskset:
            die("frozen taskset path changed before job launch")
        command = [str(taskset), "-c", str(cpu), *make_command(binary, job)]
        start_ns = time.time_ns()
        start_monotonic_ns = time.monotonic_ns()
        with tempfile.TemporaryFile() as stdout_file, \
                tempfile.TemporaryFile() as stderr_file:
            process = subprocess.Popen(
                command, stdout=stdout_file, stderr=stderr_file,
                start_new_session=True,
            )
            registered = False
            try:
                registry.add(process)
                registered = True
                deadline = time.monotonic() + timeout
                while True:
                    returncode = process.poll()
                    if bounded_capture_exceeded(
                            (stdout_file, stderr_file),
                            JOB_OUTPUT_MAX_BYTES):
                        stop_and_reap_process_group(
                            process, graceful=False)
                        die(
                            f"job {job.job} output exceeded the bounded "
                            "capture limit")
                    if returncode is not None:
                        break
                    remaining = deadline - time.monotonic()
                    reason = None
                    if abort.is_set():
                        reason = "campaign abort"
                    elif remaining <= 0:
                        reason = f"timeout after {timeout:g}s"
                    if reason is not None:
                        stop_and_reap_process_group(
                            process, graceful=False)
                        die(f"job {job.job} terminated on {reason}")
                    time.sleep(min(PROCESS_POLL_SECONDS, remaining))
            finally:
                leader_live = process.poll() is None
                descendants_live = process_group_exists(process)
                if leader_live or descendants_live:
                    stop_and_reap_process_group(
                        process, graceful=False)
                else:
                    close_process_streams(process)
                if process.poll() is None or process_group_exists(process):
                    die(f"job {job.job} child cleanup was not proven")
                if registered:
                    registry.remove(process)
                if not leader_live and descendants_live:
                    die(f"job {job.job} left an unexpected child process")
            stdout_bytes = read_bounded_capture(
                stdout_file, JOB_OUTPUT_MAX_BYTES,
                f"job {job.job} stdout")
            stderr_bytes = read_bounded_capture(
                stderr_file, JOB_OUTPUT_MAX_BYTES,
                f"job {job.job} stderr")
        end_monotonic_ns = time.monotonic_ns()
        end_ns = time.time_ns()
    finally:
        pool.release(cpu)

    stdout_path = result_root / "stdout" / f"{job.stem}.csv"
    stderr_path = result_root / "stderr" / f"{job.stem}.txt"
    command_path = result_root / "commands" / f"{job.stem}.json"
    atomic_write(stdout_path, stdout_bytes)
    atomic_write(stderr_path, stderr_bytes)
    try:
        stdout = stdout_bytes.decode("utf-8")
        stderr = stderr_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        die(f"job {job.job} emitted non-UTF-8 output: {exc}")
    if process.returncode != 0 or stderr_bytes:
        die(
            f"job {job.job} failed rc={process.returncode}: "
            f"{stderr[-1000:]}"
        )
    options = arm_options(job.arm, job.band)
    rows = parse_bench_output(
        stdout, job.ks, "0x0", job.seed, job.schedule,
        options, True,
    )
    if len(rows) != len(job.ks):
        die(f"job {job.job} result cardinality mismatch")
    record = {
        "job": job.job, "arm": job.arm, "band": job.band,
        "seed_index": job.seed_index, "seed": job.seed,
        "schedule": job.schedule, "group": job.group,
        "group_ledger_salt_unused": job.ledger_salt,
        "active_packet_peel_seed_xor": "0x0",
        "K_count": len(job.ks), "cpu": cpu, "command": command,
        "start_ns": start_ns, "end_ns": end_ns,
        "start_monotonic_ns": start_monotonic_ns,
        "end_monotonic_ns": end_monotonic_ns,
        "elapsed_ns": end_monotonic_ns - start_monotonic_ns,
        "returncode": process.returncode,
        "stdout": str(stdout_path.relative_to(result_root)),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr": str(stderr_path.relative_to(result_root)),
        "stderr_sha256": sha256_file(stderr_path),
        "saturated_timing_speed_claim_valid": False,
    }
    atomic_json(command_path, record)
    return record


def sha_manifest(root: Path, paths: Iterable[Path]) -> bytes:
    records = []
    root = root.resolve(strict=True)
    for path in sorted(paths):
        resolved = path.resolve(strict=True)
        if resolved == root or root not in resolved.parents:
            die(f"manifest target escapes root: {path}")
        records.append(f"{sha256_file(resolved)}  {resolved.relative_to(root)}\n")
    return "".join(records).encode("utf-8")


def verify_sha_manifest(root: Path, manifest: Path) -> dict[Path, str]:
    root = root.resolve(strict=True)
    result: dict[Path, str] = {}
    for number, line in enumerate(stable_bytes(manifest).decode("utf-8").splitlines(), 1):
        digest, separator, relative = line.partition("  ")
        if (
            not separator or len(digest) != 64 or
            any(character not in "0123456789abcdef" for character in digest)
        ):
            die(f"{manifest}:{number}: malformed SHA256 record")
        target = (root / relative).resolve(strict=True)
        if target == root or root not in target.parents or target in result:
            die(f"{manifest}:{number}: unsafe or duplicate target")
        if sha256_file(target) != digest:
            die(f"{manifest}:{number}: target hash mismatch")
        result[target] = digest
    if not result:
        die(f"empty SHA256 manifest: {manifest}")
    return result


PINNED_GIT_ENVIRONMENT = {
    "GIT_ATTR_NOSYSTEM": "1",
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_TERMINAL_PROMPT": "0",
}


def pinned_git_environment() -> dict[str, str]:
    """Return an ambient-copy environment with Git redirects removed."""
    environment = os.environ.copy()
    for key in tuple(environment):
        if key.startswith("GIT_"):
            environment.pop(key)
    environment.update(PINNED_GIT_ENVIRONMENT)
    environment["LANG"] = "C"
    environment["LC_ALL"] = "C"
    return environment


def validate_pinned_git_environment(
    environment: Any,
    extra_git_fields: Iterable[str] = (),
) -> None:
    """Reject ambient Git control variables outside an explicit policy."""
    if not isinstance(environment, dict) or any(
            not isinstance(key, str) or not isinstance(value, str)
            for key, value in environment.items()):
        die("pinned Git environment is not a string mapping")
    extra = set(extra_git_fields)
    if any(not isinstance(key, str) or not key.startswith("GIT_")
           for key in extra):
        die("pinned Git environment extension is malformed")
    allowed = set(PINNED_GIT_ENVIRONMENT) | extra
    unexpected = {key for key in environment
                  if key.startswith("GIT_") and key not in allowed}
    if (unexpected or
            any(environment.get(key) != value
                for key, value in PINNED_GIT_ENVIRONMENT.items()) or
            environment.get("LANG") != "C" or
            environment.get("LC_ALL") != "C"):
        die("pinned Git environment differs from its exact policy")


def tracked_clean_source(
    script: Path,
    helper: Path,
    git: Path,
) -> tuple[Path, str]:
    git_environment = pinned_git_environment()
    validate_pinned_git_environment(git_environment)
    root_result = run_bounded_process_group(
        (
            str(git), "--no-replace-objects", "-C", str(script.parent),
            "rev-parse", "--show-toplevel",
        ),
        timeout=30.0, context="source repository discovery",
        env=git_environment)
    try:
        root_text = root_result.stdout.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise CampaignError("source repository path is not UTF-8") from exc
    if root_result.returncode or root_result.stderr or not root_text:
        die("source repository discovery failed")
    repo = Path(root_text).resolve(strict=True)
    for args, name in (
        (("diff", "--quiet", "HEAD", "--"), "worktree"),
        (("diff", "--cached", "--quiet", "HEAD", "--"), "index"),
    ):
        result = run_bounded_process_group(
            (str(git), "--no-replace-objects", "-C", str(repo), *args),
            timeout=30.0, context=f"source {name} status",
            env=git_environment)
        if result.returncode or result.stdout or result.stderr:
            die(f"tracked source {name} is dirty")
    head_result = run_bounded_process_group(
        (
            str(git), "--no-replace-objects", "-C", str(repo),
            "rev-parse", "HEAD",
        ), timeout=30.0, context="source HEAD discovery",
        env=git_environment)
    try:
        head = head_result.stdout.decode("ascii").strip()
    except UnicodeDecodeError:
        head = ""
    if (head_result.returncode or head_result.stderr or
            re.fullmatch(r"[0-9a-f]{40}", head) is None):
        die("source HEAD discovery failed")
    for path in (script, helper):
        relative = path.relative_to(repo)
        tracked_result = run_bounded_process_group(
            (
                str(git), "--no-replace-objects", "-C", str(repo),
                "cat-file", "blob", f"{head}:{relative}",
            ), timeout=30.0, context=f"immutable source {relative}",
            env=git_environment)
        if (tracked_result.returncode or tracked_result.stderr or
                tracked_result.stdout != stable_bytes(path)):
            die(f"{relative} does not exactly match immutable HEAD")
    return repo, head


def verify_frozen_sources_at_commit(
    repo: Path,
    commit: str,
    sources: Sequence[tuple[Path, Path]],
    git: Path,
) -> None:
    """Bind copied controller bytes to immutable Git objects, not the worktree."""
    if not sources or re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        die("frozen source verification lacks an immutable commit")
    git_environment = pinned_git_environment()
    validate_pinned_git_environment(git_environment)
    for source, frozen in sources:
        try:
            if not source.is_absolute() or ".." in source.parts:
                die("frozen source inventory contains a noncanonical path")
            # The source Path was canonicalized when prepare began.  Keep its
            # lexical repository identity here: resolving it again would let a
            # transient symlink swap redirect this check to another Git path.
            relative = source.relative_to(repo)
            result = run_bounded_process_group(
                (
                    str(git), "--no-replace-objects", "-C", str(repo),
                    "cat-file", "blob", f"{commit}:{relative}",
                ), timeout=30.0,
                context=f"frozen source object {relative}",
                env=git_environment)
            if result.returncode or result.stderr:
                die(f"cannot read immutable source object for {source}")
            expected = result.stdout
        except (OSError, ValueError) as exc:
            raise CampaignError(
                f"cannot read immutable source object for {source}") from exc
        if stable_bytes(frozen) != expected:
            die(
                f"frozen source {frozen.name} differs from immutable "
                f"commit {commit}"
            )


PINNED_BUILD_PATH = "/usr/bin:/bin"
PINNED_BUILD_TIMEOUT_SECONDS = 900.0
PINNED_CONFIGURE_TIMEOUT_SECONDS = 180.0
PINNED_OUTPUT_MAX_BYTES = JOB_OUTPUT_MAX_BYTES
PINNED_PROJECT_CACHE = {
    "BUILD_SHARED_LIBS": ("BOOL", "OFF"),
    "BUILD_TESTS": ("BOOL", "ON"),
    "BUILD_CODEC_V2": ("BOOL", "ON"),
    "MARCH_NATIVE": ("BOOL", "ON"),
    "WIREHAIR_BUILD_BOTH": ("BOOL", "OFF"),
    "WIREHAIR_STATIC_PIC": ("BOOL", "ON"),
    "WIREHAIR_BUILD_TOOLS": ("BOOL", "OFF"),
    "WIREHAIR_BUILD_BENCHMARKS": ("BOOL", "ON"),
    "WIREHAIR_ENABLE_SCHEDULED_TESTS": ("BOOL", "OFF"),
    "WIREHAIR_ENABLE_LIBFUZZER": ("BOOL", "OFF"),
    "WIREHAIR_STRICT_WARNINGS": ("BOOL", "ON"),
    "WH_LTO": ("STRING", "OFF"),
    "WH_PGO_MODE": ("STRING", "OFF"),
}
PINNED_LANGUAGE_CACHE = {
    "CMAKE_BUILD_TYPE": ("STRING", "Release"),
    "CMAKE_CXX_STANDARD": ("STRING", "11"),
    "CMAKE_CXX_STANDARD_REQUIRED": ("BOOL", "ON"),
    "CMAKE_CXX_EXTENSIONS": ("BOOL", "ON"),
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
    "CMAKE_STATIC_LINKER_FLAGS": ("STRING", ""),
    "CMAKE_STATIC_LINKER_FLAGS_RELEASE": ("STRING", ""),
    "CMAKE_INTERPROCEDURAL_OPTIMIZATION": ("BOOL", "OFF"),
    "CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE": ("BOOL", "OFF"),
}
PINNED_EXPLICIT_EMPTY_CACHE = {
    "CMAKE_C_COMPILER_LAUNCHER": "STRING",
    "CMAKE_CXX_COMPILER_LAUNCHER": "STRING",
    "CMAKE_C_LINKER_LAUNCHER": "STRING",
    "CMAKE_CXX_LINKER_LAUNCHER": "STRING",
    "CMAKE_RULE_LAUNCH_COMPILE": "STRING",
    "CMAKE_RULE_LAUNCH_LINK": "STRING",
}
PINNED_FORBIDDEN_CACHE = {
    "CMAKE_TOOLCHAIN_FILE": "FILEPATH",
    "CMAKE_PROJECT_TOP_LEVEL_INCLUDES": "STRING",
    "CMAKE_PROJECT_INCLUDE": "FILEPATH",
    "CMAKE_PROJECT_INCLUDE_BEFORE": "FILEPATH",
    "CMAKE_PROJECT_wirehair_INCLUDE": "FILEPATH",
    "CMAKE_PROJECT_wirehair_INCLUDE_BEFORE": "FILEPATH",
    "CMAKE_USER_MAKE_RULES_OVERRIDE": "FILEPATH",
    "CMAKE_USER_MAKE_RULES_OVERRIDE_C": "FILEPATH",
    "CMAKE_USER_MAKE_RULES_OVERRIDE_CXX": "FILEPATH",
    "CMAKE_MODULE_PATH": "STRING",
    "CMAKE_PREFIX_PATH": "STRING",
    "CMAKE_PROGRAM_PATH": "STRING",
    "CMAKE_INCLUDE_PATH": "STRING",
    "CMAKE_LIBRARY_PATH": "STRING",
    "CMAKE_FIND_ROOT_PATH": "STRING",
    "CMAKE_SYSROOT": "PATH",
    "CMAKE_SYSROOT_COMPILE": "PATH",
    "CMAKE_SYSROOT_LINK": "PATH",
    "CMAKE_C_COMPILER_ARG1": "STRING",
    "CMAKE_CXX_COMPILER_ARG1": "STRING",
    "CMAKE_C_COMPILER_TARGET": "STRING",
    "CMAKE_CXX_COMPILER_TARGET": "STRING",
    "CMAKE_C_COMPILER_EXTERNAL_TOOLCHAIN": "PATH",
    "CMAKE_CXX_COMPILER_EXTERNAL_TOOLCHAIN": "PATH",
}
PINNED_TOOL_CACHE_NAMES = {
    "cmake": None,
    "git": None,
    "taskset": None,
    "ninja": "CMAKE_MAKE_PROGRAM",
    "cc": "CMAKE_C_COMPILER",
    "cxx": "CMAKE_CXX_COMPILER",
    "ar": "CMAKE_AR",
    "ranlib": "CMAKE_RANLIB",
    "ld": "CMAKE_LINKER",
    "nm": "CMAKE_NM",
    "objcopy": "CMAKE_OBJCOPY",
    "objdump": "CMAKE_OBJDUMP",
    "strip": "CMAKE_STRIP",
    "readelf": None,
}


@dataclass(frozen=True)
class FreshPinnedBuild:
    binary: Path
    cache: Path
    record: dict[str, Any]


def _parse_cmake_cache(
    cache: Path,
) -> tuple[bytes, dict[str, tuple[str, str]]]:
    try:
        raw = stable_bytes(cache)
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CampaignError("CMake cache is not UTF-8") from exc
    entries: dict[str, tuple[str, str]] = {}
    for number, line in enumerate(lines, 1):
        if not line or line.startswith(("#", "//")):
            continue
        match = re.fullmatch(r"([^:=]+):([^=]+)=(.*)", line)
        if match is None:
            die(f"CMake cache line {number} is malformed")
        key, kind, value = match.groups()
        if key in entries:
            die(f"CMake cache repeats {key}")
        entries[key] = (kind, value)
    return raw, entries


def _identity_for_executable(path: Path, context: str) -> dict[str, str]:
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        die(f"required {context} executable is unavailable: {exc}")
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        die(f"required {context} path is not an executable regular file")
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved, require_unique=False),
    }


def _trusted_executable(name: str) -> dict[str, str]:
    located = shutil.which(name, path=PINNED_BUILD_PATH)
    if located is None:
        die(f"trusted build PATH lacks {name}")
    return _identity_for_executable(Path(located), name)


def _validate_tool_inventory(tools: Any) -> None:
    if not isinstance(tools, dict) or set(tools) != set(
            PINNED_TOOL_CACHE_NAMES):
        die("pinned build tool inventory changed")
    for name, record in tools.items():
        frozen_executable_path(record, f"pinned-build {name}")


def _process_group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


class _PinnedBuildInterrupted(BaseException):
    pass


class _PinnedBuildSignalGuard:
    """Defer terminal signals until child and workspace cleanup completes."""

    def __init__(self) -> None:
        self.signals = (signal.SIGINT, signal.SIGTERM)
        self.previous: dict[int, Any] = {}
        self.previous_mask: set[signal.Signals] | None = None
        self.pending: list[int] = []
        self.acquisition_depth = 0

    def _block(self) -> None:
        previous_mask = signal.pthread_sigmask(
            signal.SIG_BLOCK, self.signals)
        if self.previous_mask is None:
            self.previous_mask = previous_mask

    def _handler(self, signum: int, _frame: Any) -> None:
        if not self.pending:
            self.pending.append(signum)
        if self.acquisition_depth:
            return
        self._block()
        raise _PinnedBuildInterrupted()

    def __enter__(self) -> "_PinnedBuildSignalGuard":
        if threading.current_thread() is not threading.main_thread():
            die("fresh pinned build must run on the main thread")
        # Keep the signals blocked until activate() runs as the first
        # statement inside the entered `with` body.  That guarantees Python
        # has registered __exit__ before a pending signal can invoke us.
        self._block()
        try:
            for signum in self.signals:
                self.previous[signum] = signal.getsignal(signum)
                if self.previous[signum] != signal.SIG_IGN:
                    signal.signal(signum, self._handler)
        except BaseException:
            for signum, previous in self.previous.items():
                signal.signal(signum, previous)
            if self.previous_mask is not None:
                signal.pthread_sigmask(
                    signal.SIG_SETMASK, self.previous_mask)
            raise
        return self

    def activate(self) -> None:
        if self.previous_mask is None:
            die("pinned-build signal guard was not entered")
        signal.pthread_sigmask(signal.SIG_SETMASK, self.previous_mask)

    @contextmanager
    def deferred_acquisition(self) -> Iterator[None]:
        """Defer raising until a newly created resource is registered."""
        self.acquisition_depth += 1
        try:
            yield
        finally:
            self.acquisition_depth -= 1
            if self.pending:
                # The caller's assignment has completed before context exit,
                # so its encompassing exception handler can now clean up.
                self._block()
                raise _PinnedBuildInterrupted()

    def block_for_cleanup(self) -> None:
        self._block()

    def __exit__(self, _kind: Any, _error: Any, _traceback: Any) -> bool:
        self._block()
        restore_error: BaseException | None = None
        for signum, previous in self.previous.items():
            try:
                signal.signal(signum, previous)
            except BaseException as error:
                if restore_error is None:
                    restore_error = error
        if self.previous_mask is not None:
            signal.pthread_sigmask(signal.SIG_SETMASK, self.previous_mask)
        if self.pending:
            os.kill(os.getpid(), self.pending[0])
            raise CampaignError(
                "pinned-build interruption handler returned")
        if restore_error is not None:
            raise restore_error
        return False


def _terminate_process_group(
    process: subprocess.Popen[Any],
    *,
    graceful: bool = True,
    capture_streams: Sequence[Any] = (),
    output_limit: int | None = None,
) -> None:
    if ((capture_streams and output_limit is None) or
            (output_limit is not None and
             (not isinstance(output_limit, int) or
              isinstance(output_limit, bool) or output_limit <= 0))):
        die("process-group cleanup capture policy is invalid")

    def capture_overflowed() -> bool:
        return output_limit is not None and bounded_capture_exceeded(
            capture_streams, output_limit)

    if graceful and _process_group_exists(process.pid):
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    if graceful:
        deadline = time.monotonic() + 1.0
        while (_process_group_exists(process.pid) and
               time.monotonic() < deadline):
            process.poll()
            if capture_overflowed():
                break
            time.sleep(PROCESS_POLL_SECONDS)
    if _process_group_exists(process.pid):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + 2.0
    while _process_group_exists(process.pid) and time.monotonic() < deadline:
        process.poll()
        time.sleep(PROCESS_POLL_SECONDS)
    try:
        process.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        pass
    if _process_group_exists(process.pid):
        die("pinned build process group survived forced cleanup")


def _run_pinned_command(
    command: Sequence[str],
    env: dict[str, str],
    temporary: Path,
    timeout: float,
    context: str,
    signal_guard: _PinnedBuildSignalGuard,
) -> tuple[bytes, bytes]:
    if (not command or any(not isinstance(part, str) or not part
                           for part in command) or
            not math.isfinite(timeout) or timeout <= 0):
        die(f"invalid {context} command")
    with tempfile.TemporaryFile(dir=str(temporary)) as stdout_file, \
            tempfile.TemporaryFile(dir=str(temporary)) as stderr_file:
        process: subprocess.Popen[Any] | None = None
        try:
            with signal_guard.deferred_acquisition():
                process = subprocess.Popen(
                    tuple(command), stdin=subprocess.DEVNULL,
                    stdout=stdout_file, stderr=stderr_file, env=env,
                    cwd=str(temporary.parent),
                    start_new_session=True,
                )
            deadline = time.monotonic() + timeout
            while True:
                returncode = process.poll()
                sizes = (os.fstat(stdout_file.fileno()).st_size,
                         os.fstat(stderr_file.fileno()).st_size)
                if any(size > PINNED_OUTPUT_MAX_BYTES for size in sizes):
                    _terminate_process_group(
                        process, graceful=False,
                        capture_streams=(stdout_file, stderr_file),
                        output_limit=PINNED_OUTPUT_MAX_BYTES)
                    die(f"{context} output exceeded the bounded capture limit")
                if returncode is not None:
                    break
                if time.monotonic() >= deadline:
                    _terminate_process_group(
                        process,
                        capture_streams=(stdout_file, stderr_file),
                        output_limit=PINNED_OUTPUT_MAX_BYTES)
                    die(f"{context} timed out")
                time.sleep(PROCESS_POLL_SECONDS)
            if _process_group_exists(process.pid):
                _terminate_process_group(
                    process,
                    capture_streams=(stdout_file, stderr_file),
                    output_limit=PINNED_OUTPUT_MAX_BYTES)
                die(f"{context} left a surviving descendant")
            stdout = read_bounded_capture(
                stdout_file, PINNED_OUTPUT_MAX_BYTES, f"{context} stdout")
            stderr = read_bounded_capture(
                stderr_file, PINNED_OUTPUT_MAX_BYTES, f"{context} stderr")
            if returncode:
                detail = stderr[-4096:].decode("utf-8", errors="replace")
                die(f"{context} failed ({returncode}): {detail.strip()}")
            return stdout, stderr
        except BaseException:
            if (process is not None and
                    (process.poll() is None or
                     _process_group_exists(process.pid))):
                _terminate_process_group(
                    process,
                    capture_streams=(stdout_file, stderr_file),
                    output_limit=PINNED_OUTPUT_MAX_BYTES)
            raise


@contextmanager
def deferred_process_group_acquisition() -> Iterator[None]:
    """Defer terminal handlers until a Popen handle has been assigned.

    Blocking signals across ``Popen`` is not safe here: the child inherits the
    caller's signal mask and would then ignore controller TERM/INT requests
    after exec.  Temporarily catching the signals leaves the child mask
    unchanged; caught dispositions reset to default on exec.
    """
    if threading.current_thread() is not threading.main_thread():
        die("process-group acquisition must run on the main thread")
    signals = (signal.SIGINT, signal.SIGTERM)
    previous: dict[int, Any] = {}
    pending: list[tuple[int, Any]] = []

    def defer(signum: int, frame: Any) -> None:
        if not pending:
            pending.append((signum, frame))

    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    try:
        for signum in signals:
            previous[signum] = signal.getsignal(signum)
            if previous[signum] != signal.SIG_IGN:
                signal.signal(signum, defer)
    except BaseException:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        raise
    # The child must inherit the entry mask, never our atomic-transition mask.
    signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
    try:
        yield
    finally:
        # Restore both handlers atomically with respect to terminal delivery.
        # The process handle has already been assigned, so a handler raised by
        # the final unmask is covered by the caller's lifecycle cleanup.
        signal.pthread_sigmask(signal.SIG_BLOCK, signals)
        restore_error: BaseException | None = None
        for signum, handler in previous.items():
            try:
                signal.signal(signum, handler)
            except BaseException as error:
                if restore_error is None:
                    restore_error = error
        unmask_error: BaseException | None = None
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        except BaseException as error:
            unmask_error = error
        if restore_error is not None:
            raise restore_error
        if unmask_error is not None:
            raise unmask_error
        if pending:
            signum, frame = pending[0]
            handler = previous[signum]
            if callable(handler):
                handler(signum, frame)
            elif handler == signal.SIG_DFL:
                # Re-sending a default-fatal signal here would terminate the
                # controller before its encompassing exception handler could
                # reap the just-created group.  Raise only after handle
                # assignment.  The managed wrapper uses its own replaying
                # guard; unmanaged production callers install explicit
                # TimingHost handlers.
                if signum == signal.SIGINT:
                    signal.default_int_handler(signum, frame)
                die("process-group acquisition interrupted by SIGTERM")


@contextmanager
def deferred_terminal_signals() -> Iterator[None]:
    """Block terminal signals only while an existing group is reaped."""
    previous_mask = signal.pthread_sigmask(
        signal.SIG_BLOCK, (signal.SIGINT, signal.SIGTERM))
    try:
        yield
    except BaseException as cleanup_error:
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        except BaseException as replay_error:
            # A failed teardown can mean descendants remain live.  Preserve
            # that safety-critical failure while retaining the deferred
            # interruption as its cause.
            raise cleanup_error from replay_error
        raise
    else:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


def stop_and_reap_process_group_deferred(
    process: subprocess.Popen[bytes], grace_seconds: float = 5.0,
    *,
    on_reaped: Callable[[], None] | None = None,
) -> None:
    """Reap a group before replaying signals, distinguishing cleanup failure."""
    if on_reaped is not None and not callable(on_reaped):
        die("process-group cleanup callback is not callable")
    with deferred_terminal_signals():
        try:
            stop_and_reap_process_group(process, grace_seconds)
        except CampaignError as error:
            raise ProcessGroupCleanupError(
                f"child process group {process.pid} cleanup failed: {error}"
            ) from error
        if on_reaped is not None:
            on_reaped()


def run_bounded_process_group(
    command: Sequence[str],
    *,
    timeout: float,
    context: str,
    input_bytes: bytes | None = None,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    output_limit: int | None = None,
    manage_signals: bool = True,
) -> subprocess.CompletedProcess[bytes]:
    """Run one command in a private session and prove full-group teardown."""
    if output_limit is None:
        output_limit = PINNED_OUTPUT_MAX_BYTES
    if (not command or any(not isinstance(part, str) or not part
                           for part in command) or
            not isinstance(context, str) or not context or
            not math.isfinite(timeout) or timeout <= 0 or
            not isinstance(output_limit, int) or
            isinstance(output_limit, bool) or output_limit <= 0 or
            not isinstance(manage_signals, bool) or
            (input_bytes is not None and
             not isinstance(input_bytes, bytes))):
        die("invalid bounded process-group command")

    def execute(
        signal_guard: _PinnedBuildSignalGuard | None,
    ) -> subprocess.CompletedProcess[bytes]:
        def terminate(
            process: subprocess.Popen[Any], *, graceful: bool = True,
        ) -> None:
            if signal_guard is not None:
                signal_guard.block_for_cleanup()
                _terminate_process_group(
                    process, graceful=graceful,
                    capture_streams=(stdout_file, stderr_file),
                    output_limit=output_limit)
                return
            with deferred_terminal_signals():
                _terminate_process_group(
                    process, graceful=graceful,
                    capture_streams=(stdout_file, stderr_file),
                    output_limit=output_limit)

        with tempfile.TemporaryFile() as stdout_file, \
                tempfile.TemporaryFile() as stderr_file:
            input_file: Any | None = None
            process: subprocess.Popen[Any] | None = None
            try:
                if input_bytes is not None:
                    input_file = tempfile.TemporaryFile()
                    input_file.write(input_bytes)
                    input_file.seek(0)
                acquisition = (
                    signal_guard.deferred_acquisition()
                    if signal_guard is not None else
                    deferred_process_group_acquisition())
                with acquisition:
                    process = subprocess.Popen(
                        tuple(command),
                        stdin=(input_file if input_file is not None
                               else subprocess.DEVNULL),
                        stdout=stdout_file, stderr=stderr_file,
                        cwd=(str(cwd) if cwd is not None else None), env=env,
                        start_new_session=True,
                    )
                deadline = time.monotonic() + timeout
                while True:
                    returncode = process.poll()
                    sizes = (os.fstat(stdout_file.fileno()).st_size,
                             os.fstat(stderr_file.fileno()).st_size)
                    if any(size > output_limit for size in sizes):
                        terminate(process, graceful=False)
                        die(f"{context} output exceeded the bounded capture limit")
                    if returncode is not None:
                        break
                    if time.monotonic() >= deadline:
                        terminate(process)
                        raise BoundedProcessTimeout(f"{context} timed out")
                    time.sleep(PROCESS_POLL_SECONDS)
                if _process_group_exists(process.pid):
                    terminate(process)
                    die(f"{context} left a surviving descendant")
                return subprocess.CompletedProcess(
                    tuple(command), returncode,
                    stdout=read_bounded_capture(
                        stdout_file, output_limit, f"{context} stdout"),
                    stderr=read_bounded_capture(
                        stderr_file, output_limit, f"{context} stderr"))
            except BaseException:
                if (process is not None and
                        (process.poll() is None or
                         _process_group_exists(process.pid))):
                    terminate(process)
                raise
            finally:
                if input_file is not None:
                    input_file.close()

    if not manage_signals:
        return execute(None)
    signal_guard = _PinnedBuildSignalGuard()
    with signal_guard:
        signal_guard.activate()
        return execute(signal_guard)


def _snapshot_path(path_bytes: bytes) -> tuple[str, ...]:
    try:
        text = path_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CampaignError("Git tree contains a non-UTF-8 path") from exc
    path = PurePosixPath(text)
    if (not text or text.startswith("/") or path.is_absolute() or
            "\\" in text or
            any(part in ("", ".", "..") for part in path.parts) or
            any(ord(character) < 32 or ord(character) == 127
                for character in text)):
        die("Git tree contains a noncanonical source path")
    return path.parts


def _write_snapshot_file(path: Path, data: bytes, executable: bool) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o700 if executable else 0o600)
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                die("short write while materializing immutable source")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _verify_source_snapshot(
    source: Path,
    manifest: Sequence[dict[str, Any]],
) -> None:
    expected = {entry["path"]: entry for entry in manifest}
    if len(expected) != len(manifest):
        die("immutable source manifest repeats a path")
    actual: set[str] = set()
    expected_directories = {"."}
    for relative in expected:
        parent = PurePosixPath(relative).parent
        while parent != PurePosixPath("."):
            expected_directories.add(parent.as_posix())
            parent = parent.parent
    actual_directories: set[str] = set()
    for root_text, directories, files in os.walk(source, followlinks=False):
        root = Path(root_text)
        relative_root = root.relative_to(source).as_posix()
        actual_directories.add(
            "." if relative_root == "." else relative_root)
        root_stat = root.lstat()
        if (not stat.S_ISDIR(root_stat.st_mode) or
                stat.S_IMODE(root_stat.st_mode) != 0o555):
            die("immutable source directory permissions changed")
        for name in directories:
            child = root / name
            if child.is_symlink() or not child.is_dir():
                die("immutable source contains an indirect directory")
        for name in files:
            path = root / name
            relative = path.relative_to(source).as_posix()
            entry = expected.get(relative)
            metadata = path.lstat()
            wanted_mode = 0o555 if entry and entry["mode"] == "100755" \
                else 0o444
            if (entry is None or not stat.S_ISREG(metadata.st_mode) or
                    metadata.st_nlink != 1 or
                    stat.S_IMODE(metadata.st_mode) != wanted_mode or
                    metadata.st_size != entry["bytes"] or
                    sha256_file(path) != entry["sha256"]):
                die(f"immutable source snapshot changed at {relative}")
            actual.add(relative)
    if actual != set(expected) or actual_directories != expected_directories:
        die("immutable source snapshot file inventory changed")


def _materialize_commit_snapshot(
    repo: Path,
    commit: str,
    git: Path,
    source: Path,
    env: dict[str, str],
    temporary: Path,
    signal_guard: _PinnedBuildSignalGuard,
) -> dict[str, Any]:
    if source.exists() or os.path.lexists(str(source)):
        die("fresh source snapshot path already exists")
    source.mkdir(mode=0o700)
    head, _stderr = _run_pinned_command(
        (str(git), "--no-replace-objects", "-C", str(repo),
         "rev-parse", "HEAD"),
        env, temporary, 30.0, "immutable source HEAD lookup", signal_guard)
    if head.strip().decode("ascii", errors="strict") != commit:
        die("clean Git environment names a different source commit")
    tree_raw, _stderr = _run_pinned_command(
        (str(git), "--no-replace-objects", "-C", str(repo),
         "rev-parse", f"{commit}^{{tree}}"),
        env, temporary, 30.0, "immutable source tree lookup", signal_guard)
    tree_oid = tree_raw.strip().decode("ascii", errors="strict")
    if re.fullmatch(r"[0-9a-f]{40}", tree_oid) is None:
        die("immutable source tree object is malformed")
    timestamp_raw, _stderr = _run_pinned_command(
        (str(git), "--no-replace-objects", "-C", str(repo),
         "show", "-s", "--format=%ct", commit),
        env, temporary, 30.0, "immutable source timestamp lookup",
        signal_guard)
    try:
        timestamp = int(timestamp_raw.strip(), 10)
    except ValueError:
        die("immutable source commit timestamp is malformed")
    if timestamp < 0:
        die("immutable source commit timestamp is negative")
    tree, _stderr = _run_pinned_command(
        (str(git), "--no-replace-objects", "-C", str(repo),
         "ls-tree", "-r", "-z", "--full-tree", commit),
        env, temporary, 60.0, "immutable source tree inventory",
        signal_guard)
    manifest: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw_entry in tree.split(b"\0"):
        if not raw_entry:
            continue
        try:
            header, path_bytes = raw_entry.split(b"\t", 1)
            mode, kind, oid = header.decode("ascii").split(" ")
        except (UnicodeDecodeError, ValueError) as exc:
            raise CampaignError("Git tree entry is malformed") from exc
        if (mode not in ("100644", "100755") or kind != "blob" or
                re.fullmatch(r"[0-9a-f]{40}", oid) is None):
            die("Git tree contains a symlink, submodule, or nonregular file")
        parts = _snapshot_path(path_bytes)
        relative = PurePosixPath(*parts).as_posix()
        if relative in seen:
            die("Git tree repeats a source path")
        seen.add(relative)
        data, _stderr = _run_pinned_command(
            (str(git), "--no-replace-objects", "-C", str(repo),
             "cat-file", "blob", oid),
            env, temporary, 30.0, f"immutable source blob {relative}",
            signal_guard)
        target = source.joinpath(*parts)
        _write_snapshot_file(target, data, mode == "100755")
        manifest.append({
            "path": relative, "mode": mode, "git_blob_oid": oid,
            "bytes": len(data), "sha256": sha256_bytes(data),
        })
    if not manifest:
        die("immutable source tree is empty")
    for root_text, directories, files in os.walk(source, topdown=False):
        root = Path(root_text)
        for name in files:
            path = root / name
            relative = path.relative_to(source).as_posix()
            entry = next(item for item in manifest
                         if item["path"] == relative)
            path.chmod(0o555 if entry["mode"] == "100755" else 0o444)
        for name in directories:
            (root / name).chmod(0o555)
    source.chmod(0o555)
    manifest.sort(key=lambda entry: entry["path"])
    _verify_source_snapshot(source, manifest)
    return {
        "commit": commit, "tree_oid": tree_oid,
        "materialized_path": str(source),
        "commit_timestamp": timestamp,
        "file_count": len(manifest),
        "manifest_sha256": sha256_bytes(json_bytes(manifest)),
        "manifest": manifest,
    }


def _tree_digest(root: Path) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            continue
        if not stat.S_ISREG(metadata.st_mode):
            die(f"provenance tree contains a nonregular entry: {path}")
        data = stable_bytes(path)
        records.append({
            "path": path.relative_to(root).as_posix(),
            "mode": stat.S_IMODE(metadata.st_mode),
            "bytes": len(data), "sha256": sha256_bytes(data),
        })
    if not records:
        die(f"provenance tree is empty: {root}")
    return {
        "root": str(root), "file_count": len(records),
        "bytes": sum(record["bytes"] for record in records),
        "manifest_sha256": sha256_bytes(json_bytes(records)),
    }


def _cmake_root_identity(
    cmake: Path,
    env: dict[str, str],
    temporary: Path,
    signal_guard: _PinnedBuildSignalGuard,
) -> dict[str, Any]:
    output, _stderr = _run_pinned_command(
        (str(cmake), "--system-information"), env, temporary, 60.0,
        "CMake module-root inventory", signal_guard)
    matches = re.findall(rb'^CMAKE_ROOT "([^"\r\n]+)"$', output,
                         flags=re.MULTILINE)
    if len(matches) != 1:
        die("CMake system information lacks one module root")
    try:
        root = Path(matches[0].decode("utf-8")).resolve(strict=True)
    except (UnicodeDecodeError, OSError, RuntimeError) as exc:
        raise CampaignError("CMake module root is unavailable") from exc
    return _tree_digest(root)


def _compiler_identity(
    compiler: Path,
    env: dict[str, str],
    temporary: Path,
    *,
    native_probe: bool,
    signal_guard: _PinnedBuildSignalGuard,
) -> dict[str, Any]:
    record: dict[str, Any] = _identity_for_executable(
        compiler, compiler.name)
    language = "c++" if native_probe else "c"
    predefines, predefine_stderr = _run_pinned_command(
        (str(compiler), "-dM", "-E", "-x", language, "/dev/null"),
        env, temporary, 30.0, "compiler family probe", signal_guard)
    if (predefine_stderr or b"#define __clang__" in predefines or
            not re.search(rb"^#define __GNUC__ [1-9][0-9]*$", predefines,
                          flags=re.MULTILINE)):
        die(
            "fresh pinned benchmark builds currently require a native GCC "
            "toolchain; Clang and other driver families are unsupported")
    record["family"] = "gcc"
    record["predefines_sha256"] = sha256_bytes(predefines)
    queries = {
        "version_sha256": (str(compiler), "--version"),
        "dumpmachine_sha256": (str(compiler), "-dumpmachine"),
        "dumpspecs_sha256": (str(compiler), "-dumpspecs"),
        "search_dirs_sha256": (str(compiler), "-print-search-dirs"),
        "sysroot_sha256": (str(compiler), "-print-sysroot"),
    }
    for field, command in queries.items():
        stdout, stderr = _run_pinned_command(
            command, env, temporary, 30.0, f"compiler query {field}",
            signal_guard)
        record[field] = sha256_bytes(stdout + b"\0" + stderr)
    if native_probe:
        output = temporary / "native-probe.o"
        stdout, stderr = _run_pinned_command(
            (str(compiler), "-march=native", "-Q", "--help=target",
             "-c", "-x", "c++", "/dev/null", "-o", str(output)),
            env, temporary, 60.0, "native compiler target probe",
            signal_guard)
        try:
            output.unlink()
        except FileNotFoundError:
            die("native compiler target probe did not create its object")
        record["native_target_sha256"] = sha256_bytes(
            stdout + b"\0" + stderr)
    subprograms: dict[str, dict[str, str]] = {}
    for name in ("cc1", "cc1plus", "collect2", "lto-wrapper", "as", "ld"):
        stdout, _stderr = _run_pinned_command(
            (str(compiler), f"-print-prog-name={name}"),
            env, temporary, 30.0, f"compiler subprogram query {name}",
            signal_guard)
        text = stdout.strip().decode("utf-8", errors="strict")
        if not text:
            die(f"compiler did not identify its {name} subprogram")
        path = Path(text)
        if not path.is_absolute():
            located = shutil.which(text, path=env["PATH"])
            if located is None:
                die(f"compiler {name} subprogram is unavailable")
            path = Path(located)
        subprograms[name] = _identity_for_executable(path, name)
    record["subprograms"] = subprograms
    return record


def _validate_compiler_record_paths(record: Any, context: str) -> None:
    required_hashes = {
        "predefines_sha256", "version_sha256", "dumpmachine_sha256",
        "dumpspecs_sha256", "search_dirs_sha256", "sysroot_sha256",
    }
    if (not isinstance(record, dict) or
            record.get("family") != "gcc" or
            not required_hashes.issubset(record) or
            not all(_canonical_digest(record.get(field))
                    for field in required_hashes) or
            not isinstance(record.get("subprograms"), dict) or
            set(record["subprograms"]) != {
                "cc1", "cc1plus", "collect2", "lto-wrapper", "as", "ld"}):
        die(f"pinned {context} compiler record is malformed")
    if ((context == "cxx") != ("native_target_sha256" in record) or
            (context == "cxx" and
             not _canonical_digest(record.get("native_target_sha256")))):
        die(f"pinned {context} native-target record is malformed")
    frozen_executable_path(
        {"path": record.get("path"), "sha256": record.get("sha256")},
        f"pinned {context} compiler")
    for name, identity in record["subprograms"].items():
        frozen_executable_path(identity, f"pinned compiler {name}")


def _current_compiler_records(
    tools: dict[str, dict[str, str]],
    env: dict[str, str],
    temporary: Path,
    signal_guard: _PinnedBuildSignalGuard,
) -> dict[str, Any]:
    return {
        "cc": _compiler_identity(
            Path(tools["cc"]["path"]), env, temporary,
            native_probe=False, signal_guard=signal_guard),
        "cxx": _compiler_identity(
            Path(tools["cxx"]["path"]), env, temporary,
            native_probe=True, signal_guard=signal_guard),
    }


def _pinned_graph_identity(
    build: Path,
    ninja: Path,
    env: dict[str, str],
    temporary: Path,
    signal_guard: _PinnedBuildSignalGuard,
) -> dict[str, Any]:
    relative_files = (
        "CMakeCache.txt", "build.ninja", "CMakeFiles/rules.ninja",
        "compile_commands.json",
    )
    files = []
    for relative in relative_files:
        path = build / relative
        data = stable_bytes(path)
        files.append({
            "path": relative, "bytes": len(data),
            "sha256": sha256_bytes(data),
        })
    commands, stderr = _run_pinned_command(
        (str(ninja), "-C", str(build), "-t", "commands",
         "wirehair_v2_bench"),
        env, temporary, 60.0, "Ninja benchmark command inventory",
        signal_guard)
    if stderr:
        die("Ninja command inventory wrote unexpected stderr")
    normalized = commands.replace(str(build).encode(), b"<BUILD>")
    source = build.parent / "source"
    normalized = normalized.replace(str(source).encode(), b"<SOURCE>")
    forbidden = (
        b"-fsanitize", b"--coverage", b"-fprofile", b" -pg",
        b"-flto", b"-fplugin", b"-specs", b"--sysroot",
        b"-fuse-ld", b" -include", b"-imacros", b" -B",
        b"-DWH_SEED_KNOBS",
        b"-DWH_PEELCAP",
    )
    if any(token in normalized for token in forbidden):
        die("Ninja benchmark command graph contains a forbidden option")
    return {
        "files": files,
        "files_manifest_sha256": sha256_bytes(json_bytes(files)),
        "target_commands_bytes": len(commands),
        "target_commands_sha256": sha256_bytes(commands),
        "normalized_target_commands_sha256": sha256_bytes(normalized),
    }


def _open_cleanup_child_directory(
    parent_fd: int,
    name: str,
    expected_identity: tuple[int, int],
) -> int:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_DIRECTORY", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(name, flags, dir_fd=parent_fd)
    except PermissionError:
        # File modes are irrelevant to unlink and must never be changed: a
        # rejected hardlink may name an inode outside staging.  Directories do
        # need traversal permission, so restore it and then bind the opened
        # descriptor back to the inode observed by the parent.
        path_flags = (getattr(os, "O_PATH", 0) |
                      getattr(os, "O_CLOEXEC", 0) |
                      getattr(os, "O_DIRECTORY", 0) |
                      getattr(os, "O_NOFOLLOW", 0))
        if not getattr(os, "O_PATH", 0):
            die("descriptor-bound cleanup requires Linux O_PATH")
        path_fd = os.open(name, path_flags, dir_fd=parent_fd)
        try:
            path_metadata = os.fstat(path_fd)
            if (not stat.S_ISDIR(path_metadata.st_mode) or
                    (path_metadata.st_dev, path_metadata.st_ino) !=
                    expected_identity):
                die("prepared-result cleanup directory changed before chmod")
            os.chmod(f"/proc/self/fd/{path_fd}", 0o700)
        finally:
            os.close(path_fd)
        descriptor = os.open(name, flags, dir_fd=parent_fd)
    try:
        opened = os.fstat(descriptor)
        named = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (not stat.S_ISDIR(opened.st_mode) or
                (opened.st_dev, opened.st_ino) != expected_identity or
                (named.st_dev, named.st_ino) != expected_identity):
            die("prepared-result cleanup directory identity changed")
        os.fchmod(descriptor, 0o700)
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _make_tree_writable(root: Path) -> None:
    """Restore directory traversal for private-workspace removal only."""
    try:
        root_metadata = root.lstat()
    except FileNotFoundError:
        return
    if not stat.S_ISDIR(root_metadata.st_mode):
        return
    root.chmod(0o700)
    for root_text, directories, _files in os.walk(
            root, topdown=True, followlinks=False):
        directory = Path(root_text)
        directory.chmod(0o700)
        for name in directories:
            child = directory / name
            metadata = child.lstat()
            if stat.S_ISDIR(metadata.st_mode):
                child.chmod(0o700)


def _remove_cleanup_contents(directory_fd: int) -> None:
    for name in sorted(os.listdir(directory_fd)):
        if name in ("", ".", ".."):
            die("prepared-result cleanup found a noncanonical entry")
        metadata = os.stat(
            name, dir_fd=directory_fd, follow_symlinks=False)
        if stat.S_ISDIR(metadata.st_mode):
            identity = (metadata.st_dev, metadata.st_ino)
            child_fd = _open_cleanup_child_directory(
                directory_fd, name, identity)
            try:
                _remove_cleanup_contents(child_fd)
            finally:
                os.close(child_fd)
            confirmed = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
            if (not stat.S_ISDIR(confirmed.st_mode) or
                    (confirmed.st_dev, confirmed.st_ino) != identity):
                die("prepared-result cleanup child identity changed")
            os.rmdir(name, dir_fd=directory_fd)
        else:
            # Unlinking does not require write permission on the inode and is
            # safe for symlinks, FIFOs, and rejected hardlinks.
            os.unlink(name, dir_fd=directory_fd)
    os.fsync(directory_fd)


def _remove_owned_staging_directory(
    parent_fd: int,
    staging_fd: int,
    staging_name: str,
    expected_identity: tuple[int, int],
) -> None:
    opened = os.fstat(staging_fd)
    named = os.stat(
        staging_name, dir_fd=parent_fd, follow_symlinks=False)
    if (not stat.S_ISDIR(opened.st_mode) or
            (opened.st_dev, opened.st_ino) != expected_identity or
            not stat.S_ISDIR(named.st_mode) or
            (named.st_dev, named.st_ino) != expected_identity):
        die("prepared-result staging identity changed before cleanup")
    os.fchmod(staging_fd, 0o700)
    _remove_cleanup_contents(staging_fd)
    confirmed = os.stat(
        staging_name, dir_fd=parent_fd, follow_symlinks=False)
    if (not stat.S_ISDIR(confirmed.st_mode) or
            (confirmed.st_dev, confirmed.st_ino) != expected_identity):
        die("prepared-result staging identity changed during cleanup")
    os.rmdir(staging_name, dir_fd=parent_fd)
    os.fsync(parent_fd)


class _PreparedResultPublicationCommitted(CampaignError):
    """Publication renamed successfully but durability could not be proved."""


def _stat_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )


def _directory_stat_identity(value: os.stat_result) -> tuple[int, ...]:
    # Renaming a directory changes its ctime even though its durable contents
    # and all other relevant metadata are unchanged.
    return (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns,
    )


def _prepared_tree_inventory(
    root: Path,
    expected_identity: tuple[int, int],
) -> tuple[tuple[Any, ...], ...]:
    """Return a stable, content-bound inventory of one private plain tree."""
    root_metadata = root.lstat()
    if (not stat.S_ISDIR(root_metadata.st_mode) or
            (root_metadata.st_dev, root_metadata.st_ino) !=
            expected_identity):
        die("prepared-result staging identity changed")
    inventory: list[tuple[Any, ...]] = []
    expected_directories = {"."}
    visited_directories: set[str] = set()
    directory_records: list[
        tuple[Path, tuple[int, ...], tuple[str, ...]]
    ] = []

    def walk_error(error: OSError) -> None:
        raise CampaignError(
            "cannot inspect complete prepared-result tree: "
            f"{error.filename}: {error}") from error

    for directory_text, names, filenames in os.walk(
            root, topdown=True, onerror=walk_error, followlinks=False):
        names.sort()
        filenames.sort()
        directory = Path(directory_text)
        directory_metadata = directory.lstat()
        if not stat.S_ISDIR(directory_metadata.st_mode):
            die(f"refusing non-directory in prepared result: {directory}")
        relative_directory = directory.relative_to(root)
        relative_text = relative_directory.as_posix()
        if relative_text in visited_directories:
            die(f"prepared-result directory visited twice: {directory}")
        visited_directories.add(relative_text)
        directory_identity = _directory_stat_identity(directory_metadata)
        entry_names = tuple(sorted((*names, *filenames)))
        directory_records.append(
            (directory, directory_identity, entry_names))
        inventory.append((
            "directory", relative_text,
            stat.S_IMODE(directory_metadata.st_mode),
            directory_identity,
        ))
        for name in names:
            child = directory / name
            child_metadata = child.lstat()
            if not stat.S_ISDIR(child_metadata.st_mode):
                die(f"refusing symlink or non-directory in prepared result: {child}")
            expected_directories.add(child.relative_to(root).as_posix())
        for name in filenames:
            child = directory / name
            child_metadata = child.lstat()
            if (not stat.S_ISREG(child_metadata.st_mode) or
                    child_metadata.st_nlink != 1):
                die(f"refusing nonunique non-file in prepared result: {child}")
            data = stable_bytes(child)
            confirmed = child.lstat()
            if _stat_identity(confirmed) != _stat_identity(child_metadata):
                die(f"prepared-result file changed during inventory: {child}")
            inventory.append((
                "file", child.relative_to(root).as_posix(),
                stat.S_IMODE(confirmed.st_mode), confirmed.st_size,
                sha256_bytes(data), _stat_identity(confirmed),
            ))
    if visited_directories != expected_directories:
        die("prepared-result inventory skipped a named directory")
    for (directory, expected_directory_identity,
         expected_entry_names) in reversed(directory_records):
        confirmed_directory = directory.lstat()
        if (not stat.S_ISDIR(confirmed_directory.st_mode) or
                _directory_stat_identity(confirmed_directory) !=
                expected_directory_identity):
            die(f"prepared-result directory changed during inventory: {directory}")
        if tuple(sorted(os.listdir(directory))) != expected_entry_names:
            die(f"prepared-result directory entries changed during inventory: {directory}")
    confirmed_root = root.lstat()
    if ((confirmed_root.st_dev, confirmed_root.st_ino) != expected_identity or
            not stat.S_ISDIR(confirmed_root.st_mode)):
        die("prepared-result staging identity changed during inventory")
    return tuple(inventory)


def _fsync_prepared_tree(
    root: Path,
    expected_identity: tuple[int, int],
) -> tuple[tuple[Any, ...], ...]:
    """Flush a stable, symlink-free prepared tree before its commit rename."""
    before_inventory = _prepared_tree_inventory(root, expected_identity)
    directories: list[tuple[Path, tuple[int, ...]]] = []
    for kind, relative, *_fields in before_inventory:
        path = root if relative == "." else root / relative
        expected_stat_identity = _fields[-1]
        if kind == "directory":
            directories.append((path, expected_stat_identity))
            continue
        parent_fd = open_durable_directory(path.parent)
        flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NOFOLLOW", 0) |
                 getattr(os, "O_NONBLOCK", 0))
        try:
            descriptor = os.open(path.name, flags, dir_fd=parent_fd)
            try:
                descriptor_before = os.fstat(descriptor)
                named_before = os.stat(
                    path.name, dir_fd=parent_fd, follow_symlinks=False)
                if (not stat.S_ISREG(descriptor_before.st_mode) or
                        descriptor_before.st_nlink != 1 or
                        _stat_identity(descriptor_before) !=
                        expected_stat_identity or
                        _stat_identity(named_before) !=
                        _stat_identity(descriptor_before)):
                    die(f"prepared-result file changed before flush: {path}")
                os.fsync(descriptor)
                descriptor_after = os.fstat(descriptor)
                named_after = os.stat(
                    path.name, dir_fd=parent_fd, follow_symlinks=False)
                if (_stat_identity(descriptor_after) !=
                        _stat_identity(descriptor_before) or
                        _stat_identity(named_after) !=
                        _stat_identity(descriptor_after)):
                    die(f"prepared-result file changed during flush: {path}")
            finally:
                os.close(descriptor)
        finally:
            os.close(parent_fd)
    for directory, expected_stat_identity in sorted(
            directories, key=lambda item: len(item[0].parts), reverse=True):
        descriptor = open_durable_directory(directory)
        try:
            before = os.fstat(descriptor)
            if (not stat.S_ISDIR(before.st_mode) or
                    _directory_stat_identity(before) !=
                    expected_stat_identity):
                die(f"prepared-result directory changed before flush: {directory}")
            os.fsync(descriptor)
            after = os.fstat(descriptor)
            named_after = directory.lstat()
            if (_directory_stat_identity(after) !=
                    _directory_stat_identity(before) or
                    _directory_stat_identity(named_after) !=
                    _directory_stat_identity(after)):
                die(f"prepared-result directory changed during flush: {directory}")
        finally:
            os.close(descriptor)
    after_inventory = _prepared_tree_inventory(root, expected_identity)
    if after_inventory != before_inventory:
        die("prepared-result tree changed while it was flushed")
    return after_inventory


def _rename_directory_noreplace(
    source: Path,
    target: Path,
    expected_identity: tuple[int, int],
    expected_inventory: tuple[tuple[Any, ...], ...],
) -> None:
    """Atomically publish a sibling directory without replacing a target."""
    if (not source.is_absolute() or not target.is_absolute() or
            source.parent != target.parent or
            source.name in ("", ".", "..") or
            target.name in ("", ".", "..")):
        die("prepared-result publication paths are noncanonical")
    directory_fd = open_durable_directory(source.parent)
    committed = False
    try:
        before = os.stat(
            source.name, dir_fd=directory_fd, follow_symlinks=False)
        if (not stat.S_ISDIR(before.st_mode) or
                (before.st_dev, before.st_ino) != expected_identity):
            die("prepared-result staging source identity changed")
        if (_prepared_tree_inventory(source, expected_identity) !=
                expected_inventory):
            die("prepared-result tree changed after its durability flush")
        try:
            renameat2 = ctypes.CDLL(None, use_errno=True).renameat2
        except AttributeError as exc:
            raise CampaignError(
                "Linux renameat2 is required for result publication") from exc
        renameat2.argtypes = (
            ctypes.c_int, ctypes.c_char_p,
            ctypes.c_int, ctypes.c_char_p, ctypes.c_uint,
        )
        renameat2.restype = ctypes.c_int
        result = renameat2(
            directory_fd, os.fsencode(source.name),
            directory_fd, os.fsencode(target.name), 1)  # RENAME_NOREPLACE
        if result != 0:
            error_number = ctypes.get_errno()
            if error_number == errno.EEXIST:
                die("prepared-result publication target already exists")
            raise CampaignError(
                "prepared-result no-replace publication failed: "
                f"{os.strerror(error_number)}")
        committed = True
        try:
            after = os.stat(
                target.name, dir_fd=directory_fd, follow_symlinks=False)
            if (not stat.S_ISDIR(after.st_mode) or
                    (after.st_dev, after.st_ino) != expected_identity):
                die("prepared-result publication identity changed")
            if (_prepared_tree_inventory(target, expected_identity) !=
                    expected_inventory):
                die("prepared-result tree changed at its commit rename")
            os.fsync(directory_fd)
            durable = os.stat(
                target.name, dir_fd=directory_fd, follow_symlinks=False)
            if (not stat.S_ISDIR(durable.st_mode) or
                    (durable.st_dev, durable.st_ino) != expected_identity):
                die("prepared-result publication identity changed after flush")
        except Exception as exc:
            raise _PreparedResultPublicationCommitted(
                "prepared-result publication committed, but durability or "
                "identity confirmation failed; preserve and inspect the "
                f"published target {target}: {exc}") from exc
    finally:
        try:
            os.close(directory_fd)
        except OSError:
            # A read-only directory close cannot change the commit outcome.
            # Before rename it is still a normal preparation failure; after
            # rename, parent fsync and identity checks above are authoritative.
            if not committed and sys.exc_info()[0] is None:
                raise


def _attach_secondary_failure(
    active_error: BaseException,
    attribute: str,
    error: BaseException,
    context: str,
) -> None:
    """Keep the primary exception while making a cleanup failure observable."""
    try:
        setattr(active_error, attribute, error)
    except Exception:
        pass
    detail = f"{context}: {type(error).__name__}: {error}"
    add_note = getattr(active_error, "add_note", None)
    if callable(add_note):
        add_note(detail)
        return
    try:
        arguments = active_error.args
        if arguments and isinstance(arguments[0], str):
            active_error.args = (
                f"{arguments[0]} [{detail}]", *arguments[1:])
        else:
            active_error.args = (*arguments, detail)
    except Exception:
        pass


@contextmanager
def prepare_result_staging(requested: Path) -> Iterator[tuple[Path, Path]]:
    """Build a result privately and publish it only after full validation."""
    absolute = requested.absolute()
    parent_fd = open_durable_directory(
        absolute.parent, create=True, reprove_existing=True)
    try:
        parent = absolute.parent.resolve(strict=True)
        opened_parent = os.fstat(parent_fd)
        named_parent = os.stat(parent, follow_symlinks=False)
        if (not stat.S_ISDIR(opened_parent.st_mode) or
                _directory_stat_identity(opened_parent) !=
                _directory_stat_identity(named_parent)):
            die("prepared-result parent path changed after durable creation")
        final = parent / absolute.name
        if (not final.name or
                _entry_exists_at(parent_fd, final.name)):
            die("prepared-result publication target already exists or is empty")

        transaction_guard = _PinnedBuildSignalGuard()
        with transaction_guard:
            staging: Path | None = None
            staging_fd: int | None = None
            identity: tuple[int, int] | None = None
            staging_created = False
            published = False
            try:
                # The cleanup finally is registered before signals are first
                # unblocked.  A nested fresh-build guard may replay into this
                # outer guard, which defers the user's signal until the result
                # staging tree has been removed or committed.
                transaction_guard.activate()
                for _attempt in range(128):
                    name = (
                        f".{final.name}.prepare-{os.getpid()}-"
                        f"{secrets.token_hex(16)}")
                    staging = parent / name
                    try:
                        with transaction_guard.deferred_acquisition():
                            os.mkdir(name, 0o700, dir_fd=parent_fd)
                            staging_created = True
                            metadata = os.stat(
                                name, dir_fd=parent_fd,
                                follow_symlinks=False)
                            if not stat.S_ISDIR(metadata.st_mode):
                                die("prepared-result staging creation changed")
                            identity = (metadata.st_dev, metadata.st_ino)
                            flags = (
                                os.O_RDONLY |
                                getattr(os, "O_CLOEXEC", 0) |
                                getattr(os, "O_DIRECTORY", 0) |
                                getattr(os, "O_NOFOLLOW", 0)
                            )
                            staging_fd = os.open(
                                name, flags, dir_fd=parent_fd)
                            opened = os.fstat(staging_fd)
                            if ((opened.st_dev, opened.st_ino) != identity or
                                    not stat.S_ISDIR(opened.st_mode)):
                                die("prepared-result staging creation raced")
                        break
                    except FileExistsError:
                        staging = None
                        staging_created = False
                        identity = None
                        continue
                else:
                    die("cannot allocate a unique prepared-result staging name")

                if staging is None or staging_fd is None or identity is None:
                    die("prepared-result staging acquisition was incomplete")
                os.fchmod(staging_fd, 0o700)
                path_metadata = staging.lstat()
                if (not stat.S_ISDIR(path_metadata.st_mode) or
                        (path_metadata.st_dev, path_metadata.st_ino) !=
                        identity):
                    die("prepared-result staging path changed after creation")
                os.fsync(staging_fd)
                os.fsync(parent_fd)

                yield staging, final
                inventory = _fsync_prepared_tree(staging, identity)
                old_mask = None
                if hasattr(signal, "pthread_sigmask"):
                    old_mask = signal.pthread_sigmask(
                        signal.SIG_BLOCK, {signal.SIGINT, signal.SIGTERM})
                try:
                    try:
                        _rename_directory_noreplace(
                            staging, final, identity, inventory)
                        published = True
                    except _PreparedResultPublicationCommitted:
                        published = True
                        raise
                finally:
                    if old_mask is not None:
                        restore_error = None
                        try:
                            signal.pthread_sigmask(
                                signal.SIG_SETMASK, old_mask)
                        except OSError as error:
                            restore_error = error
                            try:
                                signal.pthread_sigmask(
                                    signal.SIG_SETMASK, old_mask)
                            except OSError:
                                # We changed only these two signals.  If exact
                                # restore keeps failing, undo just the additions
                                # so the prior mask is reconstructed.
                                added = {
                                    signal.SIGINT, signal.SIGTERM,
                                }.difference(old_mask)
                                try:
                                    if added:
                                        signal.pthread_sigmask(
                                            signal.SIG_UNBLOCK, added)
                                    restore_error = None
                                except OSError:
                                    pass
                            else:
                                restore_error = None
                        if restore_error is not None:
                            active_error = sys.exc_info()[1]
                            if active_error is not None:
                                _attach_secondary_failure(
                                    active_error,
                                    "publication_signal_restore_error",
                                    restore_error,
                                    "publication signal-mask restore failed",
                                )
                            elif published:
                                raise _PreparedResultPublicationCommitted(
                                    "prepared-result publication committed, "
                                    "but its signal mask could not be restored"
                                ) from restore_error
                            else:
                                raise CampaignError(
                                    "prepared-result publication signal mask "
                                    "could not be restored") from restore_error
            finally:
                # Block before any cleanup pathname/descriptor transition.  A
                # pending user signal is replayed by the guard only afterward.
                transaction_guard.block_for_cleanup()
                active_error = sys.exc_info()[1]
                cleanup_error: BaseException | None = None
                if (not published and staging is not None and
                        staging_created):
                    try:
                        if identity is None:
                            die(
                                "prepared-result staging ownership was not "
                                "captured; preserving the unresolved path")
                        if staging_fd is None:
                            staging_fd = _open_cleanup_child_directory(
                                parent_fd, staging.name, identity)
                        _remove_owned_staging_directory(
                            parent_fd, staging_fd, staging.name, identity)
                    except FileNotFoundError as missing_error:
                        try:
                            current_staging = os.stat(
                                staging.name, dir_fd=parent_fd,
                                follow_symlinks=False)
                        except FileNotFoundError:
                            current_staging = None
                        if current_staging is not None:
                            cleanup_error = CampaignError(
                                "prepared-result staging contents changed "
                                "during cleanup")
                            cleanup_error.__cause__ = missing_error
                        else:
                            # A failed pre-commit rename cannot legitimately
                            # remove the root.  If final has our inode, the
                            # commit boundary was crossed; otherwise report the
                            # vanished resource.
                            try:
                                final_metadata = os.stat(
                                    final.name, dir_fd=parent_fd,
                                    follow_symlinks=False)
                            except FileNotFoundError as error:
                                cleanup_error = CampaignError(
                                    "prepared-result staging vanished before "
                                    "cleanup")
                                cleanup_error.__cause__ = error
                            else:
                                if (identity is None or
                                        (final_metadata.st_dev,
                                         final_metadata.st_ino) != identity):
                                    cleanup_error = CampaignError(
                                        "prepared-result staging vanished and "
                                        "final has a different identity")
                                else:
                                    published = True
                    except BaseException as error:
                        cleanup_error = error
                if staging_fd is not None:
                    try:
                        os.close(staging_fd)
                    except OSError as error:
                        if cleanup_error is None:
                            cleanup_error = error
                    staging_fd = None
                if cleanup_error is not None:
                    if active_error is not None:
                        _attach_secondary_failure(
                            active_error, "prepared_result_cleanup_error",
                            cleanup_error,
                            "prepared-result cleanup also failed",
                        )
                    else:
                        raise cleanup_error
    finally:
        active_error = sys.exc_info()[1]
        try:
            os.close(parent_fd)
        except OSError as error:
            if active_error is not None:
                _attach_secondary_failure(
                    active_error, "prepared_result_parent_close_error",
                    error, "prepared-result parent close also failed")
            else:
                raise


def validate_pinned_build_record(
    record: Any,
    cache: Path,
    source_root: Path | None = None,
) -> dict[str, Any]:
    if (not isinstance(record, dict) or
            record.get("schema") != "wirehair.wh2.pinned_build.v1" or
            not isinstance(record.get("source"), dict) or
            not isinstance(record.get("tools"), dict) or
            not isinstance(record.get("environment"), dict) or
            not isinstance(record.get("configure_command"), list) or
            not isinstance(record.get("build_command"), list) or
            not _canonical_digest(record.get("binary_sha256")) or
            not _canonical_digest(record.get("cmake_cache_sha256"))):
        die("pinned build record is malformed")
    if sha256_file(cache) != record["cmake_cache_sha256"]:
        die("pinned build cache differs from its provenance record")
    _raw_cache, cache_entries = _parse_cmake_cache(cache)
    _validate_tool_inventory(record["tools"])
    compiler_records = record.get("compiler_records")
    if not isinstance(compiler_records, dict) or set(compiler_records) != {
            "cc", "cxx"}:
        die("pinned compiler inventory changed")
    for name in ("cc", "cxx"):
        _validate_compiler_record_paths(compiler_records[name], name)
    cmake_root = record.get("cmake_root")
    if (not isinstance(cmake_root, dict) or
            set(cmake_root) != {
                "root", "file_count", "bytes", "manifest_sha256"} or
            not isinstance(cmake_root.get("root"), str) or
            not _canonical_digest(cmake_root.get("manifest_sha256")) or
            _tree_digest(Path(cmake_root["root"])) != cmake_root):
        die("pinned CMake module tree changed")
    source = record["source"]
    if (re.fullmatch(r"[0-9a-f]{40}", str(source.get("commit", ""))) is None or
            re.fullmatch(r"[0-9a-f]{40}",
                         str(source.get("tree_oid", ""))) is None or
            not _canonical_digest(source.get("manifest_sha256")) or
            not isinstance(source.get("materialized_path"), str) or
            not Path(source["materialized_path"]).is_absolute() or
            not isinstance(source.get("commit_timestamp"), int) or
            isinstance(source.get("commit_timestamp"), bool) or
            source["commit_timestamp"] < 0 or
            not isinstance(source.get("manifest"), list) or
            source.get("file_count") != len(source["manifest"]) or
            sha256_bytes(json_bytes(source["manifest"])) !=
            source["manifest_sha256"]):
        die("pinned build source provenance is malformed")
    if cache_entries.get("CMAKE_HOME_DIRECTORY") != (
            "INTERNAL", source["materialized_path"]):
        die("pinned build cache source binding changed")
    environment = record["environment"]
    expected_environment_keys = {
        "PATH", "HOME", "XDG_CONFIG_HOME", "TMPDIR", "LANG", "LC_ALL",
        "TZ", "SOURCE_DATE_EPOCH", *PINNED_GIT_ENVIRONMENT,
    }
    workspace = Path(source["materialized_path"]).parent
    if (set(environment) != expected_environment_keys or
            record.get("working_directory") != str(workspace) or
            environment.get("PATH") != PINNED_BUILD_PATH or
            environment.get("LANG") != "C" or
            environment.get("LC_ALL") != "C" or
            any(environment.get(key) != value
                for key, value in PINNED_GIT_ENVIRONMENT.items()) or
            environment.get("TZ") != "UTC" or
            environment.get("SOURCE_DATE_EPOCH") !=
            str(source.get("commit_timestamp")) or
            environment.get("HOME") != str(workspace / "home") or
            environment.get("TMPDIR") != str(workspace / "tmp") or
            environment.get("XDG_CONFIG_HOME") != str(workspace / "xdg")):
        die("pinned build environment policy changed")
    configure = record["configure_command"]
    build_command = record["build_command"]
    expected_source = source["materialized_path"]
    expected_build = str(workspace / "build")
    if (len(configure) < 7 or
            configure[:7] != [
                record["tools"]["cmake"]["path"], "-S", expected_source,
                "-B", expected_build, "-G", "Ninja"] or
            expected_build not in build_command or
            record["tools"]["cmake"]["path"] not in build_command or
            record.get("legacy_graph_used") is not False or
            not isinstance(record.get("legacy_binary_hint"), str) or
            record.get("legacy_binary_hint") in configure or
            record.get("legacy_binary_hint") in build_command):
        die("pinned build command policy changed")
    policy = record.get("build_policy")
    if (not isinstance(policy, dict) or
            policy.get("cmake_cache_sha256") !=
            record["cmake_cache_sha256"] or
            policy.get("fresh_pinned") is not True):
        die("pinned build policy binding changed")
    graph = record.get("graph")
    expected_graph_files = {
        "CMakeCache.txt", "build.ninja", "CMakeFiles/rules.ninja",
        "compile_commands.json",
    }
    if (not isinstance(graph, dict) or
            not isinstance(graph.get("files"), list) or
            len(graph["files"]) != len(expected_graph_files) or
            not all(isinstance(entry, dict) for entry in graph["files"]) or
            {entry.get("path") for entry in graph["files"]
             if isinstance(entry, dict)} != expected_graph_files or
            not all(_canonical_digest(entry.get("sha256")) and
                    isinstance(entry.get("bytes"), int) and
                    entry["bytes"] >= 0 for entry in graph["files"]) or
            not _canonical_digest(graph.get("files_manifest_sha256")) or
            sha256_bytes(json_bytes(graph["files"])) !=
            graph["files_manifest_sha256"] or
            not all(_canonical_digest(graph.get(field)) for field in (
                "target_commands_sha256",
                "normalized_target_commands_sha256")) or
            not isinstance(graph.get("target_commands_bytes"), int) or
            graph["target_commands_bytes"] <= 0 or
            not _canonical_digest(record.get("legacy_cache_sha256")) or
            re.fullmatch(r"[0-9a-f]+",
                         str(record.get("binary_elf_build_id", ""))) is None or
            not all(_canonical_digest(record.get(field)) for field in (
                "configure_stdout_sha256", "configure_stderr_sha256",
                "build_stdout_sha256", "build_stderr_sha256",
                "elf_notes_sha256", "elf_dynamic_sha256"))):
        die("pinned build graph or output provenance is malformed")
    cache_graph = next(
        entry for entry in graph["files"]
        if entry.get("path") == "CMakeCache.txt")
    if cache_graph.get("sha256") != record["cmake_cache_sha256"]:
        die("pinned build graph cache binding changed")
    if source_root is not None:
        _verify_source_snapshot(source_root, source["manifest"])
    return record


def validate_pinned_build_contract_binding(
    contract: Any,
    pinned_build: Any,
    context: str,
    *,
    require_legacy_hint: bool = False,
) -> None:
    """Cross-bind the human-facing contract to exact build provenance."""
    if (not isinstance(contract, dict) or
            not isinstance(pinned_build, dict) or
            contract.get("build_policy") !=
            pinned_build.get("build_policy") or
            contract.get("build_command") !=
            pinned_build.get("build_command") or
            contract.get("source_commit") !=
            pinned_build.get("source", {}).get("commit") or
            contract.get("binary_sha256") !=
            pinned_build.get("binary_sha256") or
            (require_legacy_hint and
             contract.get("legacy_binary_hint") !=
             pinned_build.get("legacy_binary_hint"))):
        die(f"{context} build provenance differs from its pinned build")


@contextmanager
def _fresh_pinned_benchmark_build_impl(
    repo: Path,
    commit: str,
    binary_hint: Path,
    git_record: dict[str, str],
    cmake_record: dict[str, str],
    taskset_record: dict[str, str] | None,
    *,
    build_workers: int,
    cpu_set: str | None,
    signal_guard: _PinnedBuildSignalGuard,
) -> Iterator[FreshPinnedBuild]:
    """Build the benchmark from immutable Git bytes in a new private tree."""
    if (re.fullmatch(r"[0-9a-f]{40}", commit) is None or
            not isinstance(build_workers, int) or
            isinstance(build_workers, bool) or build_workers <= 0):
        die("fresh pinned build arguments are invalid")
    binary_input = binary_hint.absolute()
    if binary_input.is_symlink():
        die("--binary toolchain hint must not be a symlink")
    hinted_binary = binary_input.resolve(strict=True)
    if not hinted_binary.is_file() or not os.access(hinted_binary, os.X_OK):
        die("--binary toolchain hint must be an executable regular file")
    hinted_build = hinted_binary.parent.parent.resolve(strict=True)
    hinted_cache = hinted_build / "CMakeCache.txt"
    if hinted_binary != (
            hinted_build / "codec/wirehair_v2_bench").resolve(strict=True):
        die("--binary must name the codec/wirehair_v2_bench toolchain hint")
    hint_policy = validate_candidate_build_policy(hinted_cache, repo)
    hinted_cache_sha256 = hint_policy["cmake_cache_sha256"]
    tools = {
        "cmake": cmake_record,
        "git": git_record,
        "taskset": taskset_record,
        "ninja": _trusted_executable("ninja"),
        "cc": _trusted_executable("cc"),
        "cxx": _trusted_executable("c++"),
        "ar": _trusted_executable("ar"),
        "ranlib": _trusted_executable("ranlib"),
        "ld": _trusted_executable("ld"),
        "nm": _trusted_executable("nm"),
        "objcopy": _trusted_executable("objcopy"),
        "objdump": _trusted_executable("objdump"),
        "strip": _trusted_executable("strip"),
        "readelf": _trusted_executable("readelf"),
    }
    git = frozen_executable_path(git_record, "git")
    cmake = frozen_executable_path(cmake_record, "cmake")
    taskset = None
    if cpu_set is not None:
        if taskset_record is None:
            die("pinned build CPU set lacks taskset identity")
        taskset = frozen_executable_path(taskset_record, "taskset")
    elif taskset_record is not None:
        frozen_executable_path(taskset_record, "taskset")
    if taskset_record is None:
        die("fresh pinned build lacks a taskset identity")
    _validate_tool_inventory(tools)
    for name, supplied in (("cmake", cmake_record), ("git", git_record),
                           ("taskset", taskset_record)):
        if supplied != _trusted_executable(name):
            die(f"fresh pinned build requires the trusted system {name}")
    for field, name in (("c_compiler", "cc"), ("cxx_compiler", "cxx")):
        try:
            hinted = Path(str(hint_policy[field])).resolve(strict=True)
        except (OSError, RuntimeError) as exc:
            raise CampaignError("hinted compiler is unavailable") from exc
        if hinted != Path(tools[name]["path"]):
            die("--binary toolchain hint uses a different system compiler")

    workspace: Path | None = None
    try:
        with signal_guard.deferred_acquisition():
            workspace = Path(tempfile.mkdtemp(
                prefix="wirehair-wh2-pinned-build-"))
        workspace.chmod(0o700)
        home = workspace / "home"
        temporary = workspace / "tmp"
        xdg = workspace / "xdg"
        for directory in (home, temporary, xdg):
            directory.mkdir(mode=0o700)
        env = {
            "PATH": PINNED_BUILD_PATH,
            "HOME": str(home),
            "XDG_CONFIG_HOME": str(xdg),
            "TMPDIR": str(temporary),
            "LANG": "C", "LC_ALL": "C", "TZ": "UTC",
            **PINNED_GIT_ENVIRONMENT,
        }
        source = workspace / "source"
        source_record = _materialize_commit_snapshot(
            repo, commit, git, source, env, temporary, signal_guard)
        env["SOURCE_DATE_EPOCH"] = str(source_record["commit_timestamp"])
        cmake_root = _cmake_root_identity(
            cmake, env, temporary, signal_guard)
        compiler_records = _current_compiler_records(
            tools, env, temporary, signal_guard)
        _validate_tool_inventory(tools)
        _verify_source_snapshot(source, source_record["manifest"])
        build = workspace / "build"
        if os.path.lexists(str(build)):
            die("fresh pinned build directory already exists")
        configure_command = [str(cmake), "-S", str(source), "-B", str(build),
                             "-G", "Ninja"]
        configure_values = {
            **PINNED_PROJECT_CACHE, **PINNED_LANGUAGE_CACHE,
            "CMAKE_MAKE_PROGRAM": ("FILEPATH", tools["ninja"]["path"]),
            "CMAKE_C_COMPILER": ("FILEPATH", tools["cc"]["path"]),
            "CMAKE_CXX_COMPILER": ("FILEPATH", tools["cxx"]["path"]),
            "CMAKE_AR": ("FILEPATH", tools["ar"]["path"]),
            "CMAKE_RANLIB": ("FILEPATH", tools["ranlib"]["path"]),
            "CMAKE_LINKER": ("FILEPATH", tools["ld"]["path"]),
            "CMAKE_NM": ("FILEPATH", tools["nm"]["path"]),
            "CMAKE_OBJCOPY": ("FILEPATH", tools["objcopy"]["path"]),
            "CMAKE_OBJDUMP": ("FILEPATH", tools["objdump"]["path"]),
            "CMAKE_STRIP": ("FILEPATH", tools["strip"]["path"]),
            "CMAKE_C_COMPILER_AR": ("FILEPATH", tools["ar"]["path"]),
            "CMAKE_CXX_COMPILER_AR": ("FILEPATH", tools["ar"]["path"]),
            "CMAKE_C_COMPILER_RANLIB":
                ("FILEPATH", tools["ranlib"]["path"]),
            "CMAKE_CXX_COMPILER_RANLIB":
                ("FILEPATH", tools["ranlib"]["path"]),
            "CMAKE_EXPORT_COMPILE_COMMANDS": ("BOOL", "ON"),
            "CMAKE_SUPPRESS_REGENERATION": ("BOOL", "ON"),
            "CMAKE_DISABLE_SOURCE_CHANGES": ("BOOL", "ON"),
            "CMAKE_DISABLE_IN_SOURCE_BUILD": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_OpenMP": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_Python3": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_PkgConfig": ("BOOL", "ON"),
            "CMAKE_FIND_USE_PACKAGE_REGISTRY": ("BOOL", "OFF"),
            "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY": ("BOOL", "OFF"),
        }
        configure_values.update({
            key: (kind, "")
            for key, kind in PINNED_EXPLICIT_EMPTY_CACHE.items()
        })
        for key in sorted(configure_values):
            kind, value = configure_values[key]
            configure_command.append(f"-D{key}:{kind}={value}")
        configure_stdout, configure_stderr = _run_pinned_command(
            configure_command, env, temporary,
            PINNED_CONFIGURE_TIMEOUT_SECONDS, "fresh CMake configure",
            signal_guard)
        if not build.is_dir() or build.is_symlink():
            die("CMake did not create a direct fresh build directory")
        cache = build / "CMakeCache.txt"
        build_policy = validate_candidate_build_policy(
            cache, source, require_pinned=True, expected_tools=tools)
        binary = build / "codec/wirehair_v2_bench"
        if os.path.lexists(str(binary)):
            die("fresh configure unexpectedly created the benchmark target")
        graph = _pinned_graph_identity(
            build, Path(tools["ninja"]["path"]), env, temporary,
            signal_guard)
        prebuild_graph = _pinned_graph_identity(
            build, Path(tools["ninja"]["path"]), env, temporary,
            signal_guard)
        if graph != prebuild_graph:
            die("fresh Ninja graph changed before its build")
        _verify_source_snapshot(source, source_record["manifest"])
        _validate_tool_inventory(tools)
        if _current_compiler_records(
                tools, env, temporary, signal_guard) != compiler_records:
            die("compiler identity changed during fresh configuration")
        if _tree_digest(Path(cmake_root["root"])) != cmake_root:
            die("CMake module tree changed before the build")
        command_prefix = [] if taskset is None else [
            str(taskset), "-c", str(cpu_set)]
        build_command = [
            *command_prefix, str(cmake), "--build", str(build),
            "--target", "wirehair_v2_bench", "--parallel",
            str(build_workers), "--verbose",
        ]
        build_stdout, build_stderr = _run_pinned_command(
            build_command, env, temporary,
            PINNED_BUILD_TIMEOUT_SECONDS, "fresh benchmark build",
            signal_guard)
        if _pinned_graph_identity(
                build, Path(tools["ninja"]["path"]), env,
                temporary, signal_guard) != graph:
            die("fresh Ninja graph changed during its build")
        if validate_candidate_build_policy(
                cache, source, require_pinned=True,
                expected_tools=tools) != build_policy:
            die("fresh CMake policy changed during its build")
        _verify_source_snapshot(source, source_record["manifest"])
        _validate_tool_inventory(tools)
        if _current_compiler_records(
                tools, env, temporary, signal_guard) != compiler_records:
            die("compiler identity changed during the fresh build")
        if _tree_digest(Path(cmake_root["root"])) != cmake_root:
            die("CMake module tree changed during the build")
        binary_sha256 = sha256_file(binary)
        if not os.access(binary, os.X_OK):
            die("fresh benchmark target is not executable")
        readelf = Path(tools["readelf"]["path"])
        notes, notes_stderr = _run_pinned_command(
            (str(readelf), "--notes", str(binary)), env, temporary, 30.0,
            "fresh benchmark ELF build-ID inspection", signal_guard)
        matches = re.findall(
            rb"^[ \t]*Build ID: ([0-9a-f]+)[ \t]*$", notes,
            flags=re.MULTILINE)
        if notes_stderr or len(matches) != 1:
            die("fresh benchmark lacks one lowercase ELF build ID")
        dynamic, dynamic_stderr = _run_pinned_command(
            (str(readelf), "--dynamic", str(binary)), env, temporary, 30.0,
            "fresh benchmark ELF dynamic inspection", signal_guard)
        if dynamic_stderr or re.search(rb"\((?:RPATH|RUNPATH)\)", dynamic):
            die("fresh benchmark contains an unexpected runtime search path")
        record = {
            "schema": "wirehair.wh2.pinned_build.v1",
            "source": source_record,
            "tools": tools,
            "compiler_records": compiler_records,
            "cmake_root": cmake_root,
            "environment": env,
            "working_directory": str(workspace),
            "legacy_binary_hint": str(hinted_binary),
            "legacy_cache_sha256": hinted_cache_sha256,
            "legacy_graph_used": False,
            "configure_command": configure_command,
            "build_command": build_command,
            "configure_stdout_sha256": sha256_bytes(configure_stdout),
            "configure_stderr_sha256": sha256_bytes(configure_stderr),
            "build_stdout_sha256": sha256_bytes(build_stdout),
            "build_stderr_sha256": sha256_bytes(build_stderr),
            "build_policy": build_policy,
            "graph": graph,
            "cmake_cache_sha256": sha256_file(cache),
            "binary_sha256": binary_sha256,
            "binary_elf_build_id": matches[0].decode("ascii"),
            "elf_notes_sha256": sha256_bytes(notes),
            "elf_dynamic_sha256": sha256_bytes(dynamic),
        }
        validate_pinned_build_record(record, cache, source)
        if sha256_file(hinted_cache) != hinted_cache_sha256:
            die("--binary toolchain-hint cache changed during fresh build")
        fresh = FreshPinnedBuild(binary=binary, cache=cache, record=record)
        try:
            yield fresh
        finally:
            if sys.exc_info()[0] is None:
                if (sha256_file(binary) != binary_sha256 or
                        validate_candidate_build_policy(
                            cache, source, require_pinned=True,
                            expected_tools=tools) != build_policy or
                        _pinned_graph_identity(
                            build, Path(tools["ninja"]["path"]), env,
                            temporary, signal_guard) != graph):
                    die("fresh pinned build changed while it was staged")
                validate_pinned_build_record(record, cache, source)
    except BaseException:
        signal_guard.block_for_cleanup()
        raise
    else:
        signal_guard.block_for_cleanup()
    finally:
        if workspace is not None:
            _make_tree_writable(workspace)
            shutil.rmtree(workspace)


@contextmanager
def fresh_pinned_benchmark_build(
    repo: Path,
    commit: str,
    binary_hint: Path,
    git_record: dict[str, str],
    cmake_record: dict[str, str],
    taskset_record: dict[str, str] | None,
    *,
    build_workers: int,
    cpu_set: str | None,
) -> Iterator[FreshPinnedBuild]:
    signal_guard = _PinnedBuildSignalGuard()
    with signal_guard:
        signal_guard.activate()
        try:
            with _fresh_pinned_benchmark_build_impl(
                    repo, commit, binary_hint, git_record, cmake_record,
                    taskset_record, build_workers=build_workers,
                    cpu_set=cpu_set, signal_guard=signal_guard) as build:
                yield build
        except BaseException:
            signal_guard.block_for_cleanup()
            raise
        else:
            signal_guard.block_for_cleanup()


def validate_candidate_build_policy(
    cache: Path,
    repo: Path,
    *,
    require_pinned: bool = False,
    expected_tools: dict[str, dict[str, str]] | None = None,
) -> dict[str, Any]:
    """Require the exact uninstrumented native Release benchmark policy."""
    raw, entries = _parse_cmake_cache(cache)
    required = {
        **PINNED_PROJECT_CACHE,
        **{
            key: value for key, value in PINNED_LANGUAGE_CACHE.items()
            if (require_pinned or key not in (
                "CMAKE_CXX_STANDARD", "CMAKE_CXX_STANDARD_REQUIRED",
                "CMAKE_CXX_EXTENSIONS",
                "CMAKE_INTERPROCEDURAL_OPTIMIZATION",
                "CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE",
            ))
        },
    }
    for key, expected in required.items():
        if entries.get(key) != expected:
            die(f"candidate CMake policy requires {key}={expected[1]!r}")
    for key in (
            "CMAKE_CXX_STANDARD", "CMAKE_CXX_STANDARD_REQUIRED",
            "CMAKE_CXX_EXTENSIONS",
            "CMAKE_INTERPROCEDURAL_OPTIMIZATION",
            "CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE"):
        if key in entries and entries[key] != PINNED_LANGUAGE_CACHE[key]:
            die(f"candidate CMake policy requires {key}="
                f"{PINNED_LANGUAGE_CACHE[key][1]!r}")
    for key in PINNED_EXPLICIT_EMPTY_CACHE:
        if key in entries and entries[key][1]:
            die(f"candidate CMake policy forbids {key}")
        if require_pinned and entries.get(key) != (
                PINNED_EXPLICIT_EMPTY_CACHE[key], ""):
            die(f"fresh CMake policy must explicitly clear {key}")
    for key in PINNED_FORBIDDEN_CACHE:
        if key in entries and (require_pinned or entries[key][1]):
            die(f"candidate CMake policy forbids {key}")
    for key, (_kind, value) in entries.items():
        upper = key.upper()
        if ((upper.endswith("_COMPILER_LAUNCHER") or
             upper.endswith("_LINKER_LAUNCHER") or
             upper.startswith("CMAKE_RULE_LAUNCH_")) and value):
            die(f"candidate CMake policy forbids launcher {key}")
        if (upper.startswith("CMAKE_PROJECT_") and
                upper.endswith(("_INCLUDE", "_INCLUDE_BEFORE")) and value):
            die(f"candidate CMake policy forbids project hook {key}")
    if ("CMAKE_CONFIGURATION_TYPES" in entries and
            entries["CMAKE_CONFIGURATION_TYPES"][1]):
        die("candidate CMake policy requires a single-config generator")
    banned = (
        "-fsanitize", "--coverage", "-fprofile", "-pg", "-flto",
        "-fplugin", "-specs", "--sysroot", "-fuse-ld",
        "-DWH_SEED_KNOBS", "-DWH_PEELCAP", "-Wl,--wrap",
    )
    for key, (_kind, value) in entries.items():
        if "FLAGS" in key and any(token in value for token in banned):
            die(f"candidate CMake policy forbids injected flags in {key}")
        if (key.startswith("CMAKE_INTERPROCEDURAL_OPTIMIZATION") and
                value.upper() not in ("", "OFF", "FALSE", "0")):
            die("candidate CMake policy forbids implicit IPO/LTO")
    home = entries.get("CMAKE_HOME_DIRECTORY")
    if home is None or home[0] != "INTERNAL":
        die("candidate CMake cache lacks its source-tree binding")
    try:
        home_path = Path(home[1]).resolve(strict=True)
    except OSError as exc:
        raise CampaignError("candidate CMake source tree is unavailable") from exc
    if home_path != repo.resolve(strict=True):
        die("candidate CMake cache belongs to another source tree")
    for key in ("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "CMAKE_GENERATOR"):
        if key not in entries or not entries[key][1]:
            die(f"candidate CMake cache lacks {key}")
    allowed_generators = (("INTERNAL", "Ninja"),) if require_pinned else (
        ("INTERNAL", "Unix Makefiles"), ("INTERNAL", "Ninja"))
    if entries["CMAKE_GENERATOR"] not in allowed_generators:
        die("candidate CMake policy has an unapproved generator")
    if require_pinned:
        if expected_tools is None:
            die("fresh CMake policy lacks expected tool identities")
        for name, cache_key in PINNED_TOOL_CACHE_NAMES.items():
            if cache_key is None:
                continue
            entry = entries.get(cache_key)
            if entry is None or not entry[1]:
                die(f"fresh CMake policy lacks {cache_key}")
            try:
                observed = Path(entry[1]).resolve(strict=True)
            except (OSError, RuntimeError) as exc:
                raise CampaignError(
                    f"fresh CMake tool {cache_key} is unavailable") from exc
            if observed != Path(expected_tools[name]["path"]):
                die(f"fresh CMake policy changed {cache_key}")
        for cache_key, tool_name in (
                ("CMAKE_C_COMPILER_AR", "ar"),
                ("CMAKE_CXX_COMPILER_AR", "ar"),
                ("CMAKE_C_COMPILER_RANLIB", "ranlib"),
                ("CMAKE_CXX_COMPILER_RANLIB", "ranlib")):
            entry = entries.get(cache_key)
            if entry is None or not entry[1]:
                die(f"fresh CMake policy lacks {cache_key}")
            try:
                observed = Path(entry[1]).resolve(strict=True)
            except (OSError, RuntimeError) as exc:
                raise CampaignError(
                    f"fresh CMake tool {cache_key} is unavailable") from exc
            if observed != Path(expected_tools[tool_name]["path"]):
                die(f"fresh CMake policy changed {cache_key}")
        for key in (
                "CMAKE_GENERATOR_INSTANCE", "CMAKE_GENERATOR_PLATFORM",
                "CMAKE_GENERATOR_TOOLSET"):
            if key in entries and entries[key][1]:
                die(f"fresh CMake generator binding changed {key}")
        exact_fresh = {
            "CMAKE_EXPORT_COMPILE_COMMANDS": ("BOOL", "ON"),
            "CMAKE_SUPPRESS_REGENERATION": ("BOOL", "ON"),
            "CMAKE_DISABLE_SOURCE_CHANGES": ("BOOL", "ON"),
            "CMAKE_DISABLE_IN_SOURCE_BUILD": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_OpenMP": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_Python3": ("BOOL", "ON"),
            "CMAKE_DISABLE_FIND_PACKAGE_PkgConfig": ("BOOL", "ON"),
            "CMAKE_FIND_USE_PACKAGE_REGISTRY": ("BOOL", "OFF"),
            "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY": ("BOOL", "OFF"),
        }
        for key, expected in exact_fresh.items():
            if entries.get(key) != expected:
                die(f"fresh CMake policy requires {key}={expected[1]!r}")
    return {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.build_policy.v2",
        "cmake_cache_sha256": sha256_bytes(raw),
        "build_type": "Release", "march_native": True,
        "strict_warnings": True, "tests_enabled": True,
        "benchmarks_built": True, "lto": "OFF", "pgo": "OFF",
        "static_library": True, "static_pic": True,
        "cxx_standard": 11, "cxx_extensions": True,
        "c_compiler": entries["CMAKE_C_COMPILER"][1],
        "cxx_compiler": entries["CMAKE_CXX_COMPILER"][1],
        "generator": entries["CMAKE_GENERATOR"][1],
        "fresh_pinned": require_pinned,
    }


def _prepare_with_fresh_pinned_benchmark_build(
    args: argparse.Namespace,
    final_result_dir: Path,
) -> dict[str, Any]:
    script = Path(__file__).resolve(strict=True)
    helper = script.with_name("wh2_rank_floor_two_anchor_screen.py")
    git_record = _trusted_executable("git")
    git = frozen_executable_path(git_record, "git")
    repo, head = tracked_clean_source(script, helper, git)
    taskset_record = _trusted_executable("taskset")
    cmake_record = _trusted_executable("cmake")
    python_runtime = python_runtime_identity()
    taskset = frozen_executable_path(taskset_record, "taskset")
    cmake = frozen_executable_path(cmake_record, "cmake")
    groups_input = args.groups.absolute()
    if groups_input.is_symlink():
        die("--groups must not be a symlink")
    groups = groups_input.resolve(strict=True)
    thermal = args.thermal.resolve(strict=True)
    load_groups(groups)
    result_dir = args.result_dir.resolve()
    frozen = result_dir / "frozen"
    destinations = {
        "script": frozen / script.name,
        "helper": frozen / helper.name,
        "binary": frozen / "wirehair_v2_bench",
        "cmake_cache": frozen / "CMakeCache.txt",
        "groups": frozen / "groups.tsv",
        "pinned_build": frozen / "pinned_build.json",
    }
    with fresh_pinned_benchmark_build(
            repo, head, args.binary, git_record, cmake_record,
            taskset_record, build_workers=args.build_workers,
            cpu_set=None) as fresh_build:
        post_repo, post_head = tracked_clean_source(script, helper, git)
        if post_repo != repo or post_head != head:
            die("source repository changed during the fresh benchmark build")
        if (frozen_executable_path(git_record, "git") != git or
                frozen_executable_path(
                    taskset_record, "taskset") != taskset or
                frozen_executable_path(cmake_record, "cmake") != cmake):
            die("frozen tool path changed during the fresh benchmark build")
        thermal_mark = thermal_start(thermal)
        build_policy = fresh_build.record["build_policy"]
        build_command = fresh_build.record["build_command"]
        legacy_binary_hint = fresh_build.record["legacy_binary_hint"]
        live_binary_sha256 = fresh_build.record["binary_sha256"]
        binary_bytes = stable_bytes(fresh_build.binary)
        cache_bytes = stable_bytes(fresh_build.cache)
        pinned_build_bytes = json_bytes(fresh_build.record)
        if sha256_bytes(binary_bytes) != live_binary_sha256:
            die("fresh benchmark binary changed while it was captured")
        validate_pinned_build_record(fresh_build.record, fresh_build.cache)
    # Do not expose a final result tree until every fresh-build exit check has
    # passed.  A failed context-exit validation must remain safely retryable.
    frozen.mkdir()
    atomic_write(destinations["binary"], binary_bytes)
    atomic_write(destinations["cmake_cache"], cache_bytes)
    atomic_write(destinations["pinned_build"], pinned_build_bytes)
    if sha256_file(destinations["binary"]) != live_binary_sha256:
        die("captured fresh benchmark changed while it was frozen")
    validate_pinned_build_record(
        json.loads(pinned_build_bytes), destinations["cmake_cache"])
    for source, destination in (
        (script, destinations["script"]),
        (helper, destinations["helper"]),
        (groups, destinations["groups"]),
    ):
        shutil.copyfile(source, destination)
    # The source ledger was validated before the build, but the copy itself
    # is a second read and may race a rewrite.  Validate the exact frozen
    # bytes before sealing them into an immutable no-replace result.
    load_groups(destinations["groups"], GROUPS_SHA256)
    destinations["binary"].chmod(0o755)
    destinations["script"].chmod(0o755)
    verify_frozen_sources_at_commit(
        repo, head,
        ((script, destinations["script"]), (helper, destinations["helper"])),
        git,
    )
    if sha256_file(destinations["binary"]) != live_binary_sha256:
        die("fresh benchmark binary changed while inputs were frozen")
    pinned_build_sha256 = sha256_file(destinations["pinned_build"])
    validate_pinned_build_record(
        json.loads(stable_bytes(destinations["pinned_build"])),
        destinations["cmake_cache"],
    )
    contract = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.contract.v6",
        "source_commit": head, "source_repo": str(repo),
        "build_command": build_command,
        "legacy_binary_hint": legacy_binary_hint,
        "binary_sha256": live_binary_sha256,
        "cmake_cache_sha256": sha256_file(destinations["cmake_cache"]),
        "pinned_build_sha256": pinned_build_sha256,
        "build_policy": build_policy,
        "script_sha256": sha256_file(destinations["script"]),
        "helper_sha256": sha256_file(destinations["helper"]),
        "groups_sha256": sha256_file(destinations["groups"]),
        "thermal": str(thermal), "workers": args.workers,
        "taskset": taskset_record,
        "cmake": cmake_record,
        "git": git_record,
        "python": python_runtime,
        "thermal_baseline": {
            "dev": thermal_mark["dev"], "ino": thermal_mark["ino"],
            "offset": thermal_mark["offset"],
            "edac_ce": thermal_mark["edac_ce"],
            "edac_ue": thermal_mark["edac_ue"],
            "monotonic_s": thermal_mark["monotonic_s"],
            "max_temperature_c": thermal_mark["max_temperature_c"],
            "row_sha256": sha256_bytes(thermal_mark["baseline_row"]),
        },
        "timeout_seconds": args.timeout,
        "K_domain": [K_MIN, K_MAX], "cutoff": CUTOFF,
        "arms": ARMS, "seeds": SEEDS, "schedules": SCHEDULES,
        "loss": "0.50", "trials": 1, "block_bytes": 64,
        "independence": (
            "Fresh seeds and every K with packet-peel seed xor fixed to the "
            "named production value zero. Published R2 group membership is "
            "batching/order metadata only; its historical salt column is "
            "sealed for provenance but never applied or analyzed as a scope."
        ),
        "packet_peel_seed_xor": "0x0",
        "group_ledger_salts_applied": False,
        "saturated_timing_speed_claim_valid": False,
        "thermal_fields": THERMAL_FIELDS,
        "thermal_policy": {
            "limit_c": 90.0, "consecutive_samples": 3,
            "stale_seconds": 5.0, "min_cpu_busy_pct": 95.0,
        },
    }
    contract_path = frozen / "contract.json"
    atomic_json(contract_path, contract)
    staged_paths = [*destinations.values(), contract_path]
    staged_path = frozen / "staged.sha256"
    atomic_write(staged_path, sha_manifest(result_dir, staged_paths))
    verify_sha_manifest(result_dir, staged_path)
    prepared = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.prepare.v2",
        "result_dir": str(final_result_dir), "source_commit": head,
        "binary_sha256": contract["binary_sha256"],
        "pinned_build_sha256": pinned_build_sha256,
        "staged_sha256": sha256_file(staged_path),
        "run_command": [
            python_runtime["path"],
            str(final_result_dir / "frozen" / script.name), "run",
            "--result-dir", str(final_result_dir),
        ],
    }
    atomic_json(result_dir / "prepare.json", prepared)
    return prepared


def prepare(args: argparse.Namespace) -> int:
    with prepare_result_staging(args.result_dir) as (staging, final):
        staged_args = argparse.Namespace(**vars(args))
        staged_args.result_dir = staging
        prepared = _prepare_with_fresh_pinned_benchmark_build(
            staged_args, final)
    print(canonical_json(prepared))
    return 0


def load_frozen(
    result_dir: Path,
) -> tuple[
    dict[str, Any], Path, Path, Path, tuple[int, int, int, int], list[Group],
]:
    frozen = result_dir / "frozen"
    prepare_anchor = validate_prepare_anchor(
        result_dir,
        "wirehair.wh2.rank_floor_two_anchor_allk.prepare.v2",
        ("source_commit", "binary_sha256", "pinned_build_sha256",
         "staged_sha256", "run_command"),
        "staged_sha256")
    staged = verify_sha_manifest(result_dir, frozen / "staged.sha256")
    expected = {
        (frozen / name).resolve(strict=True)
        for name in (
            "wh2_rank_floor_two_anchor_allk.py",
            "wh2_rank_floor_two_anchor_screen.py",
            "wirehair_v2_bench", "CMakeCache.txt", "groups.tsv",
            "pinned_build.json", "contract.json",
        )
    }
    if set(staged) != expected:
        die("frozen seal does not cover the exact seven-file staged set")
    expected_names = {path.name for path in expected} | {"staged.sha256"}
    if ({path.name for path in frozen.iterdir()} != expected_names or
            any(path.is_symlink() or not path.is_file()
                for path in frozen.iterdir())):
        die("frozen all-K tree contains an unexpected artifact")
    contract = json.loads(stable_bytes(frozen / "contract.json"))
    try:
        pinned_build = json.loads(stable_bytes(frozen / "pinned_build.json"))
    except json.JSONDecodeError as exc:
        raise CampaignError("pinned build record is not JSON") from exc
    validate_pinned_build_record(
        pinned_build, frozen / "CMakeCache.txt")
    if (
        contract.get("schema") !=
            "wirehair.wh2.rank_floor_two_anchor_allk.contract.v6"
        or tuple(contract.get("arms", ())) != ARMS
        or tuple(contract.get("seeds", ())) != SEEDS
        or tuple(contract.get("schedules", ())) != SCHEDULES
        or contract.get("K_domain") != [K_MIN, K_MAX]
        or contract.get("groups_sha256") != GROUPS_SHA256
        or contract.get("packet_peel_seed_xor") != "0x0"
        or contract.get("group_ledger_salts_applied") is not False
        or tuple(contract.get("thermal_fields", ())) != THERMAL_FIELDS
        or contract.get("thermal_policy") != {
            "limit_c": 90.0, "consecutive_samples": 3,
            "stale_seconds": 5.0, "min_cpu_busy_pct": 95.0,
        }
        or contract.get("cmake_cache_sha256") !=
            sha256_file(frozen / "CMakeCache.txt")
        or contract.get("pinned_build_sha256") !=
            sha256_file(frozen / "pinned_build.json")
    ):
        die("frozen contract does not match this campaign")
    validate_pinned_build_contract_binding(
        contract, pinned_build, "frozen all-K contract",
        require_legacy_hint=True)
    script = (frozen / "wh2_rank_floor_two_anchor_allk.py").resolve(strict=True)
    helper = (frozen / "wh2_rank_floor_two_anchor_screen.py").resolve(strict=True)
    if Path(__file__).resolve(strict=True) != script:
        die("run must be invoked through the frozen campaign script")
    run_command = prepare_anchor.get("run_command")
    python_runtime = contract.get("python")
    if python_runtime != python_runtime_identity():
        die("frozen Python runtime changed")
    if (not isinstance(run_command, list) or len(run_command) != 5 or
            run_command[0] != python_runtime["path"] or
            run_command[1:] != [
                str(script), "run", "--result-dir", str(result_dir)]):
        die("prepare anchor run command changed")
    if (
        sha256_file(script) != contract["script_sha256"] or
        sha256_file(helper) != contract["helper_sha256"]
    ):
        die("frozen script or helper does not match the contract")
    binary = (frozen / "wirehair_v2_bench").resolve(strict=True)
    verify_frozen_binary(binary, contract.get("binary_sha256"))
    if (prepare_anchor.get("source_commit") != contract.get("source_commit") or
            prepare_anchor.get("binary_sha256") !=
            contract.get("binary_sha256") or
            prepare_anchor.get("pinned_build_sha256") !=
            contract.get("pinned_build_sha256")):
        die("prepare anchor differs from the frozen all-K contract")
    taskset = frozen_executable_path(contract.get("taskset"), "taskset")
    frozen_executable_path(contract.get("cmake"), "cmake")
    frozen_executable_path(contract.get("git"), "git")
    groups = load_groups(frozen / "groups.tsv", contract["groups_sha256"])
    thermal = Path(contract["thermal"]).resolve(strict=True)
    thermal_identity = validate_frozen_thermal_source(
        thermal, contract.get("thermal_baseline"),
        float(contract["thermal_policy"]["stale_seconds"]),
        cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
        dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
        consecutive_limit=int(
            contract["thermal_policy"]["consecutive_samples"]))
    return contract, binary, taskset, thermal, thermal_identity, groups


class ArmData:
    def __init__(self) -> None:
        self.seen = bytearray(K_COUNT * len(SEEDS) * len(SCHEDULES))
        size = len(self.seen)
        self.rank = bytearray(size)
        self.error = bytearray(size)
        self.shortfall = array("H", [0]) * size
        self.seed_attempt = array("H", [0]) * size
        self.inact = array("Q", [0]) * size
        self.binary_def = array("Q", [0]) * size
        self.heavy_gain = array("Q", [0]) * size
        self.xors = array("Q", [0]) * size
        self.muladds = array("Q", [0]) * size


ALLK_COMMAND_RECORD_FIELDS = frozenset({
    "job", "arm", "band", "seed_index", "seed", "schedule", "group",
    "group_ledger_salt_unused", "active_packet_peel_seed_xor", "K_count",
    "cpu", "command", "start_ns", "end_ns", "start_monotonic_ns",
    "end_monotonic_ns", "elapsed_ns", "returncode",
    "stdout", "stdout_sha256", "stderr", "stderr_sha256",
    "saturated_timing_speed_claim_valid",
})

ALLK_THERMAL_SUMMARY_FIELDS = frozenset({
    "samples", "sealed_samples_including_baseline", "cpu_busy_min_pct",
    "cpu_tctl_max_c", "dimm_max_c", "dimm_read_errors_max",
    "edac_ce_delta", "edac_ue_delta", "thermal_limit_c",
    "thermal_high_samples", "thermal_high_max_consecutive_samples",
    "guard_poll_iterations", "guard_samples", "guard_high_samples",
    "guard_limit_c", "guard_error",
    "baseline_row_sha256",
})


def _canonical_json_record(path: Path, context: str) -> tuple[bytes, Any]:
    """Load one strict canonical pretty-JSON receipt byte-for-byte."""
    data = stable_bytes(path)
    try:
        record = json.loads(data)
        canonical = json_bytes(record)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as exc:
        raise CampaignError(f"{context} is not valid canonical JSON") from exc
    if data != canonical:
        die(f"{context} is not exact canonical JSON")
    return data, record


def _require_exact_receipt_directory(
    directory: Path,
    expected_names: set[str],
    context: str,
) -> None:
    """Reject missing, extra, indirect, or non-regular receipt entries."""
    try:
        entries = list(directory.iterdir())
    except OSError as exc:
        raise CampaignError(f"cannot inspect {context}: {exc}") from exc
    if {entry.name for entry in entries} != expected_names:
        die(f"{context} does not contain the exact expected artifact set")
    for entry in entries:
        try:
            metadata = entry.lstat()
        except OSError as exc:
            raise CampaignError(f"cannot inspect {context} entry: {exc}") from exc
        if not stat.S_ISREG(metadata.st_mode):
            die(f"{context} contains a non-regular artifact: {entry.name}")


def _thermal_counter_delta(values: Sequence[int]) -> int:
    steps = [current - previous
             for previous, current in zip(values, values[1:])]
    rollback = min(steps, default=0)
    return rollback if rollback < 0 else values[-1] - values[0]


def _verify_allk_thermal_receipt(
    data: bytes,
    summary: dict[str, Any],
    policy: dict[str, Any],
    interval_start_monotonic_ns: int,
    interval_end_monotonic_ns: int,
) -> None:
    """Recompute the sealed thermal interval and exact guard summary."""
    if set(summary) != ALLK_THERMAL_SUMMARY_FIELDS:
        die("all-K thermal summary has the wrong schema")
    lines = data.splitlines(keepends=True)
    expected_header = (",".join(THERMAL_FIELDS) + "\n").encode("ascii")
    if len(lines) < 3 or lines[0] != expected_header:
        die("all-K thermal receipt lacks its exact header/baseline/interval")
    if any(not line.endswith(b"\n") or len(line) > THERMAL_ROW_MAX_BYTES
           for line in lines[1:]):
        die("all-K thermal receipt has a partial or oversized row")
    try:
        samples = [
            parse_thermal_sample(line, f"sealed thermal row {number}")
            for number, line in enumerate(lines[1:], 1)
        ]
        stale_seconds = float(policy["stale_seconds"])
        limit_c = float(policy["limit_c"])
        min_busy = float(policy["min_cpu_busy_pct"])
        consecutive_limit = int(policy["consecutive_samples"])
    except (KeyError, TypeError, ValueError) as exc:
        raise CampaignError(f"all-K thermal receipt is malformed: {exc}") from exc
    if (
        type(interval_start_monotonic_ns) is not int or
        type(interval_end_monotonic_ns) is not int or
        not 0 < interval_start_monotonic_ns <=
            interval_end_monotonic_ns < 1 << 64 or
        not math.isfinite(stale_seconds) or stale_seconds <= 0 or
        not math.isfinite(limit_c) or not 0 <= limit_c <= 120 or
        not math.isfinite(min_busy) or not 0 <= min_busy <= 100 or
        type(policy.get("consecutive_samples")) is not int or
        consecutive_limit <= 0
    ):
        die("all-K thermal policy is malformed")
    if (not _canonical_digest(summary.get("baseline_row_sha256")) or
            summary["baseline_row_sha256"] != sha256_bytes(lines[1])):
        die("all-K thermal receipt does not bind the live guard baseline")
    start_seconds = interval_start_monotonic_ns / 1_000_000_000
    end_seconds = interval_end_monotonic_ns / 1_000_000_000
    if (samples[0]["monotonic_s"] > start_seconds or
            samples[-1]["monotonic_s"] < end_seconds):
        die("all-K thermal receipt does not cover the full job interval")
    for number, (previous, current) in enumerate(
            zip(samples, samples[1:]), 1):
        delta = current["monotonic_s"] - previous["monotonic_s"]
        if delta <= 0 or delta > stale_seconds:
            die(f"all-K thermal row {number} is nonmonotonic or follows a gap")
    interval_samples = samples[1:]
    consecutive_high = 0
    high_samples = 0
    max_consecutive_high = 0
    for sample in samples:
        high = (
            sample["cpu_tctl_c"] >= limit_c or
            max(sample["dimm_temperatures"]) >= limit_c
        )
        if high:
            consecutive_high += 1
            high_samples += 1
            max_consecutive_high = max(
                max_consecutive_high, consecutive_high)
        else:
            consecutive_high = 0
    computed = {
        "samples": len(interval_samples),
        "sealed_samples_including_baseline": len(samples),
        "cpu_busy_min_pct": min(
            sample["cpu_busy_pct"] for sample in interval_samples),
        "cpu_tctl_max_c": max(sample["cpu_tctl_c"] for sample in samples),
        "dimm_max_c": max(
            value for sample in samples
            for value in sample["dimm_temperatures"]),
        "dimm_read_errors_max": max(
            sample["dimm_read_errors"] for sample in samples),
        "edac_ce_delta": _thermal_counter_delta(
            [sample["edac_ce"] for sample in samples]),
        "edac_ue_delta": _thermal_counter_delta(
            [sample["edac_ue"] for sample in samples]),
        "thermal_limit_c": limit_c,
        "thermal_high_samples": high_samples,
        "thermal_high_max_consecutive_samples": max_consecutive_high,
    }
    for field, expected in computed.items():
        if type(summary.get(field)) is not type(expected) or \
                summary[field] != expected:
            die(f"all-K thermal summary field {field} does not replay")
    guard_integer_fields = (
        "guard_poll_iterations", "guard_samples", "guard_high_samples")
    if (
        any(type(summary.get(field)) is not int or summary[field] < 0
            for field in guard_integer_fields) or
        type(summary.get("guard_limit_c")) is not float or
        summary["guard_limit_c"] != limit_c or
        summary.get("guard_error") is not None or
        summary["guard_samples"] > computed["samples"] or
        summary["guard_high_samples"] > summary["guard_samples"] or
        summary["guard_high_samples"] > computed["thermal_high_samples"] or
        (summary["guard_samples"] > 0 and
         summary["guard_poll_iterations"] == 0)
    ):
        die("all-K live thermal guard summary is inconsistent")
    if (
        computed["thermal_high_max_consecutive_samples"] >=
            consecutive_limit or
        computed["dimm_read_errors_max"] != 0 or
        computed["edac_ce_delta"] != 0 or
        computed["edac_ue_delta"] != 0 or
        computed["cpu_busy_min_pct"] < min_busy
    ):
        die("all-K thermal receipt fails the campaign health gate")


def verify_allk_receipts(
    result_dir: Path,
    contract: dict[str, Any],
    jobs: Sequence[Job],
    binary: Path,
    taskset: Path,
    workers: int,
    worker_cpus: Sequence[int],
    thermal_summary: dict[str, Any],
    thermal_baseline_row_sha256: str,
) -> dict[Path, str]:
    """Replay every all-K receipt and return its exact pre-seal hashes."""
    result_root = result_dir / "results"
    if not jobs:
        die("all-K receipt verification has an empty job ledger")
    if any(job.job != ordinal for ordinal, job in enumerate(jobs)):
        die("all-K job ledger is not in exact contiguous job order")
    stems = [job.stem for job in jobs]
    if len(set(stems)) != len(stems):
        die("all-K job ledger contains duplicate artifact stems")
    cpus = tuple(worker_cpus)
    if (
        type(workers) is not int or workers <= 0 or len(cpus) != workers or
        any(type(cpu) is not int or cpu < 0 for cpu in cpus) or
        tuple(sorted(set(cpus))) != cpus
    ):
        die("all-K worker CPU ledger is not canonical")
    if (not isinstance(thermal_summary, dict) or
            not _canonical_digest(thermal_baseline_row_sha256) or
            thermal_summary.get("baseline_row_sha256") !=
                thermal_baseline_row_sha256):
        die("all-K thermal summary is not an object")
    try:
        timeout_seconds = float(contract["timeout_seconds"])
    except (KeyError, TypeError, ValueError) as exc:
        raise CampaignError("all-K timeout contract is malformed") from exc
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        die("all-K timeout contract is malformed")
    verification_wall_ns = time.time_ns()
    verification_monotonic_ns = time.monotonic_ns()
    verify_frozen_binary(binary, contract.get("binary_sha256"))
    if frozen_executable_path(contract.get("taskset"), "taskset") != taskset:
        die("frozen taskset path changed before receipt verification")

    command_names = {f"{stem}.json" for stem in stems}
    stdout_names = {f"{stem}.csv" for stem in stems}
    stderr_names = {f"{stem}.txt" for stem in stems}
    _require_exact_receipt_directory(
        result_root / "commands", command_names, "all-K command directory")
    _require_exact_receipt_directory(
        result_root / "stdout", stdout_names, "all-K stdout directory")
    _require_exact_receipt_directory(
        result_root / "stderr", stderr_names, "all-K stderr directory")
    try:
        root_entries = list(result_root.iterdir())
    except OSError as exc:
        raise CampaignError(f"cannot inspect all-K result directory: {exc}") from exc
    expected_root_names = {
        "commands", "stdout", "stderr", "tasks.jsonl", "run.json",
        "thermal_interval.csv",
    }
    if {entry.name for entry in root_entries} != expected_root_names:
        die("all-K result directory does not have the exact receipt layout")
    for name in ("commands", "stdout", "stderr"):
        try:
            metadata = (result_root / name).lstat()
        except OSError as exc:
            raise CampaignError(f"cannot inspect all-K {name}: {exc}") from exc
        if not stat.S_ISDIR(metadata.st_mode):
            die(f"all-K {name} is not a direct directory")

    hashes: dict[Path, str] = {}
    records: list[dict[str, Any]] = []
    allowed_cpus = set(cpus)
    integer_fields = {
        "job", "seed_index", "group", "K_count", "cpu", "start_ns",
        "end_ns", "start_monotonic_ns", "end_monotonic_ns", "elapsed_ns",
        "returncode",
    }
    intervals_by_cpu: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for job in jobs:
        command_path = result_root / "commands" / f"{job.stem}.json"
        stdout_path = result_root / "stdout" / f"{job.stem}.csv"
        stderr_path = result_root / "stderr" / f"{job.stem}.txt"
        command_data, record = _canonical_json_record(
            command_path, f"job {job.job} command receipt")
        if not isinstance(record, dict) or set(record) != \
                ALLK_COMMAND_RECORD_FIELDS:
            die(f"job {job.job} command receipt has the wrong schema")
        if any(type(record.get(field)) is not int for field in integer_fields):
            die(f"job {job.job} command receipt has a non-integer field")
        cpu = record["cpu"]
        expected_command = [
            str(taskset), "-c", str(cpu), *make_command(binary, job)
        ]
        start_ns = record["start_ns"]
        end_ns = record["end_ns"]
        start_monotonic_ns = record["start_monotonic_ns"]
        end_monotonic_ns = record["end_monotonic_ns"]
        expected_stdout = str(stdout_path.relative_to(result_root))
        expected_stderr = str(stderr_path.relative_to(result_root))
        stdout_data = stable_bytes(stdout_path)
        stderr_data = stable_bytes(stderr_path)
        if (
            record["job"] != job.job or record.get("arm") != job.arm or
            record.get("band") != job.band or
            record["seed_index"] != job.seed_index or
            record.get("seed") != job.seed or
            record.get("schedule") != job.schedule or
            record["group"] != job.group or
            record.get("group_ledger_salt_unused") != job.ledger_salt or
            record.get("active_packet_peel_seed_xor") != "0x0" or
            record["K_count"] != len(job.ks) or cpu not in allowed_cpus or
            record.get("command") != expected_command or
            min(start_ns, end_ns, start_monotonic_ns,
                end_monotonic_ns) <= 0 or
            max(start_ns, end_ns, start_monotonic_ns,
                end_monotonic_ns) >= 1 << 64 or
            end_ns < start_ns or end_monotonic_ns < start_monotonic_ns or
            end_ns > verification_wall_ns + 1_000_000_000 or
            end_monotonic_ns > verification_monotonic_ns or
            abs((start_ns - start_monotonic_ns) -
                (verification_wall_ns - verification_monotonic_ns)) >
                5_000_000_000 or
            abs((end_ns - end_monotonic_ns) -
                (verification_wall_ns - verification_monotonic_ns)) >
                5_000_000_000 or
            record["elapsed_ns"] != end_monotonic_ns - start_monotonic_ns or
            record["elapsed_ns"] >
                math.ceil((timeout_seconds + 1.0) * 1_000_000_000) or
            record["returncode"] != 0 or
            record.get("stdout") != expected_stdout or
            record.get("stderr") != expected_stderr or
            not _canonical_digest(record.get("stdout_sha256")) or
            record["stdout_sha256"] != sha256_bytes(stdout_data) or
            not _canonical_digest(record.get("stderr_sha256")) or
            record["stderr_sha256"] != sha256_bytes(stderr_data) or
            stderr_data != b"" or
            record.get("saturated_timing_speed_claim_valid") is not False
        ):
            die(f"job {job.job} command/provenance receipt mismatch")
        intervals_by_cpu[cpu].append(
            (start_monotonic_ns, end_monotonic_ns, job.job))
        try:
            stdout = stdout_data.decode("utf-8")
            rows = parse_bench_output(
                stdout, job.ks, "0x0", job.seed, job.schedule,
                arm_options(job.arm, job.band), True,
            )
        except (UnicodeDecodeError, ValueError) as exc:
            raise CampaignError(
                f"job {job.job} stdout semantic replay failed: {exc}") from exc
        if len(rows) != len(job.ks):
            die(f"job {job.job} stdout replay cardinality changed")
        records.append(record)
        hashes[command_path.resolve(strict=True)] = sha256_bytes(command_data)
        hashes[stdout_path.resolve(strict=True)] = sha256_bytes(stdout_data)
        hashes[stderr_path.resolve(strict=True)] = sha256_bytes(stderr_data)

    if set(intervals_by_cpu) != allowed_cpus:
        die("all-K receipts do not exercise the exact worker CPU ledger")
    for cpu, intervals in intervals_by_cpu.items():
        ordered = sorted(intervals)
        for previous, current in zip(ordered, ordered[1:]):
            if previous[1] > current[0]:
                die(
                    f"all-K CPU {cpu} has overlapping jobs "
                    f"{previous[2]} and {current[2]}"
                )

    tasks_path = result_root / "tasks.jsonl"
    expected_tasks = b"".join(
        (canonical_json(record) + "\n").encode("utf-8") for record in records
    )
    tasks_data = stable_bytes(tasks_path)
    if tasks_data != expected_tasks:
        die("all-K tasks ledger is not the exact canonical ordered receipt set")
    hashes[tasks_path.resolve(strict=True)] = sha256_bytes(tasks_data)

    thermal_path = result_root / "thermal_interval.csv"
    thermal_data = stable_bytes(thermal_path)
    _verify_allk_thermal_receipt(
        thermal_data, thermal_summary, contract["thermal_policy"],
        min(record["start_monotonic_ns"] for record in records),
        max(record["end_monotonic_ns"] for record in records))
    hashes[thermal_path.resolve(strict=True)] = sha256_bytes(thermal_data)

    expected_run = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.run.v5",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "jobs": len(jobs), "workers": workers,
        "worker_cpus": list(cpus),
        "start_ns": min(record["start_ns"] for record in records),
        "end_ns": max(record["end_ns"] for record in records),
        "start_monotonic_ns": min(
            record["start_monotonic_ns"] for record in records),
        "end_monotonic_ns": max(
            record["end_monotonic_ns"] for record in records),
        "recovery_cells": sum(len(job.ks) for job in jobs),
        "thermal": thermal_summary,
        "saturated_timing_speed_claim_valid": False,
    }
    run_path = result_root / "run.json"
    run_data, run_record = _canonical_json_record(
        run_path, "all-K run receipt")
    if not isinstance(run_record, dict) or run_data != json_bytes(expected_run):
        die("all-K run receipt does not exactly bind the campaign ledger")
    hashes[run_path.resolve(strict=True)] = sha256_bytes(run_data)

    expected_artifacts = 3 * len(jobs) + 3
    if len(hashes) != expected_artifacts:
        die("all-K receipt verification did not cover every artifact exactly")
    return hashes


def cell_index(seed_index: int, schedule: str, K: int) -> int:
    stratum = seed_index * len(SCHEDULES) + SCHEDULES.index(schedule)
    return stratum * K_COUNT + K - K_MIN


def milli(text: str, context: str) -> int:
    try:
        value = Decimal(text)
    except InvalidOperation:
        die(f"{context}: invalid decimal {text!r}")
    scaled = value * Decimal(1000)
    if not value.is_finite() or value < 0 or scaled != scaled.to_integral_value():
        die(f"{context}: expected nonnegative value exact to 0.001")
    result = int(scaled)
    if result >= 1 << 64:
        die(f"{context}: scaled value does not fit u64")
    return result


def sealed_phase_bytes(
    path: Path,
    phase_hashes: dict[Path, str],
    context: str,
) -> bytes:
    """Read one artifact and prove it is still the phase-sealed version."""
    resolved = path.resolve(strict=True)
    expected = phase_hashes.get(resolved)
    data = stable_bytes(path)
    if not _canonical_digest(expected) or sha256_bytes(data) != expected:
        die(f"{context} changed after the phase seal")
    return data


def load_arm_data(
    jobs: Sequence[Job], result_root: Path, phase_hashes: dict[Path, str],
) -> dict[str, ArmData]:
    data = {arm: ArmData() for arm in ARMS}
    for ordinal, job in enumerate(jobs, 1):
        stdout_path = result_root / "stdout" / f"{job.stem}.csv"
        text = sealed_phase_bytes(
            stdout_path, phase_hashes, f"{job.stem} stdout").decode("utf-8")
        rows = parse_bench_output(
            text, job.ks, "0x0", job.seed, job.schedule,
            arm_options(job.arm, job.band), True,
        )
        if len(rows) != len(job.ks):
            die(f"{job.stem} parsed row cardinality changed")
        arm = data[job.arm]
        for K, row in zip(job.ks, rows):
            index = cell_index(job.seed_index, job.schedule, K)
            if arm.seen[index]:
                die(f"duplicate analyzed cell {job.arm}/{job.seed_index}/{job.schedule}/{K}")
            rank = strict_uint(row["rank_fail"], "rank_fail")
            error = strict_uint(row["error"], "error")
            shortfall = strict_uint(row["heavy_shortfall"], "heavy_shortfall")
            attempt = strict_uint(row["seed_attempt"], "seed_attempt")
            if rank > 1 or error > 1 or rank + error > 1:
                die("single-trial row has an invalid outcome")
            if shortfall >= 1 << 16 or attempt >= 1 << 16:
                die("shortfall or seed attempt does not fit u16")
            arm.seen[index] = 1
            arm.rank[index] = rank
            arm.error[index] = error
            arm.shortfall[index] = shortfall
            arm.seed_attempt[index] = attempt
            arm.inact[index] = milli(row["inact_mu"], "inact_mu")
            arm.binary_def[index] = milli(row["binary_def_mu"], "binary_def_mu")
            arm.heavy_gain[index] = milli(row["heavy_gain_mu"], "heavy_gain_mu")
            arm.xors[index] = milli(row["block_xors_mu"], "block_xors_mu")
            arm.muladds[index] = milli(row["block_muladds_mu"], "block_muladds_mu")
        if ordinal % 500 == 0 or ordinal == len(jobs):
            print(f"analysis load {ordinal}/{len(jobs)}", flush=True)
    for arm, values in data.items():
        if values.seen.count(1) != len(values.seen):
            die(f"{arm} does not cover every expected cell exactly once")
    return data


def verify_low_identity(
    jobs: Sequence[Job], result_root: Path, phase_hashes: dict[Path, str],
) -> dict[str, int]:
    low: dict[tuple[int, str, int], dict[str, Job]] = defaultdict(dict)
    for job in jobs:
        if job.band == "small":
            low[(job.seed_index, job.schedule, job.group)][job.arm] = job
    cells = comparisons = 0
    for key, arms in low.items():
        if set(arms) != set(ARMS):
            die(f"low-K identity ledger lacks arms at {key}")
        parsed: dict[str, list[dict[str, str]]] = {}
        for arm, job in arms.items():
            path = result_root / "stdout" / f"{job.stem}.csv"
            parsed[arm] = parse_bench_output(
                sealed_phase_bytes(
                    path, phase_hashes, f"{job.stem} identity stdout"
                ).decode("utf-8"), job.ks, "0x0",
                job.seed, job.schedule, (), True,
            )
        base = parsed["d12"]
        for candidate in ARMS[1:]:
            if len(base) != len(parsed[candidate]):
                die(f"low-K identity row cardinality differs for {candidate}")
            for base_row, candidate_row in zip(base, parsed[candidate]):
                comparisons += len(base_row) - len(TIMING_FIELDS)
                if not deterministic_equal(base_row, candidate_row):
                    die(
                        f"adaptive low-K identity mismatch at {key}, "
                        f"K={base_row['N']}, arm={candidate}"
                    )
        cells += len(base)
    expected = (CUTOFF - K_MIN) * len(SEEDS) * len(SCHEDULES)
    if cells != expected:
        die(f"low-K identity covered {cells} D12 cells, want {expected}")
    return {
        "d12_cells": cells, "candidate_comparisons": cells * 2,
        "non_timing_field_comparisons": comparisons,
        "mismatches": 0,
    }


def failed(data: ArmData, index: int) -> bool:
    return bool(data.rank[index] or data.error[index])


def ratio(numerator: int, denominator: int) -> str | None:
    if denominator == 0:
        return None
    return str(Decimal(numerator) / Decimal(denominator))


def band(K: int) -> str:
    if K < CUTOFF:
        return "K00002-04095"
    if K < 8192:
        return "K04096-08191"
    if K < 16384:
        return "K08192-16383"
    if K < 32768:
        return "K16384-32767"
    return "K32768-64000"


def update_counter(
    counter: Counter[str], data: dict[str, ArmData], index: int,
) -> None:
    counter["cells"] += 1
    for arm in ARMS:
        values = data[arm]
        arm_failed = failed(values, index)
        outcome = "failure_path" if arm_failed else "success"
        counter[f"{arm}_fail"] += arm_failed
        counter[f"{arm}_{outcome}_cells"] += 1
        counter[f"{arm}_error"] += values.error[index]
        counter[f"{arm}_shortfall"] += values.shortfall[index]
        counter[f"{arm}_inact_milli"] += values.inact[index]
        counter[f"{arm}_{outcome}_inact_milli"] += values.inact[index]
        counter[f"{arm}_binary_def_milli"] += values.binary_def[index]
        counter[f"{arm}_heavy_gain_milli"] += values.heavy_gain[index]
        counter[f"{arm}_binary_def_gt15"] += values.binary_def[index] > 15000
        counter[f"{arm}_binary_def_max_milli"] = max(
            counter[f"{arm}_binary_def_max_milli"], values.binary_def[index]
        )
        counter[f"{arm}_xors_milli"] += values.xors[index]
        counter[f"{arm}_{outcome}_xors_milli"] += values.xors[index]
        counter[f"{arm}_muladds_milli"] += values.muladds[index]
        counter[f"{arm}_{outcome}_muladds_milli"] += values.muladds[index]
    base = data["d12"]
    base_fail = failed(base, index)
    for arm in ARMS[1:]:
        candidate = data[arm]
        candidate_fail = failed(candidate, index)
        counter[f"{arm}_repair"] += base_fail and not candidate_fail
        counter[f"{arm}_intro"] += not base_fail and candidate_fail
        counter[f"{arm}_both_fail"] += base_fail and candidate_fail
        counter[f"{arm}_cleared_shortfall"] += (
            bool(base.shortfall[index]) and not candidate.shortfall[index]
        )
        counter[f"{arm}_new_shortfall"] += (
            not base.shortfall[index] and bool(candidate.shortfall[index])
        )
        if not base_fail and not candidate_fail:
            counter[f"{arm}_paired_success"] += 1
            counter[f"{arm}_paired_base_inact_milli"] += base.inact[index]
            counter[f"{arm}_paired_candidate_inact_milli"] += \
                candidate.inact[index]
            counter[f"{arm}_paired_base_xors_milli"] += base.xors[index]
            counter[f"{arm}_paired_candidate_xors_milli"] += candidate.xors[index]
            counter[f"{arm}_paired_base_muladds_milli"] += base.muladds[index]
            counter[f"{arm}_paired_candidate_muladds_milli"] += candidate.muladds[index]


def counter_report(counter: Counter[str]) -> dict[str, Any]:
    arms: dict[str, Any] = {}
    for arm in ARMS:
        success_cells = counter[f"{arm}_success_cells"]
        failure_cells = counter[f"{arm}_failure_path_cells"]
        if success_cells + failure_cells != counter["cells"] or \
                failure_cells != counter[f"{arm}_fail"]:
            die(f"{arm} outcome/cost partition does not cover the exact scope")
        for metric in ("inact", "xors", "muladds"):
            if counter[f"{arm}_{metric}_milli"] != (
                    counter[f"{arm}_success_{metric}_milli"] +
                    counter[f"{arm}_failure_path_{metric}_milli"]):
                die(f"{arm} {metric} cost partition does not replay")
        arms[arm] = {
            "successful_cells": success_cells,
            "failure_path_cells": failure_cells,
            "failures": counter[f"{arm}_fail"],
            "errors": counter[f"{arm}_error"],
            "heavy_shortfall_sum": counter[f"{arm}_shortfall"],
            "inact_all_cells_milli_sum": counter[f"{arm}_inact_milli"],
            "inact_success_milli_sum":
                counter[f"{arm}_success_inact_milli"],
            "inact_failure_path_milli_sum":
                counter[f"{arm}_failure_path_inact_milli"],
            "binary_def_milli_sum": counter[f"{arm}_binary_def_milli"],
            "binary_def_max_milli": counter[f"{arm}_binary_def_max_milli"],
            "binary_def_gt15_cells": counter[f"{arm}_binary_def_gt15"],
            "heavy_gain_milli_sum": counter[f"{arm}_heavy_gain_milli"],
            "block_xors_all_cells_milli_sum": counter[f"{arm}_xors_milli"],
            "block_xors_success_milli_sum":
                counter[f"{arm}_success_xors_milli"],
            "block_xors_failure_path_milli_sum":
                counter[f"{arm}_failure_path_xors_milli"],
            "block_muladds_all_cells_milli_sum":
                counter[f"{arm}_muladds_milli"],
            "block_muladds_success_milli_sum":
                counter[f"{arm}_success_muladds_milli"],
            "block_muladds_failure_path_milli_sum":
                counter[f"{arm}_failure_path_muladds_milli"],
        }
    comparisons: dict[str, Any] = {}
    for arm in ARMS[1:]:
        repairs = counter[f"{arm}_repair"]
        introductions = counter[f"{arm}_intro"]
        both_fail = counter[f"{arm}_both_fail"]
        paired_success = counter[f"{arm}_paired_success"]
        if paired_success + repairs + introductions + both_fail != \
                counter["cells"]:
            die(f"{arm} paired outcome partition does not cover the exact scope")
        comparisons[f"{arm}_vs_d12"] = {
            "repairs": repairs, "introductions": introductions,
            "both_fail": both_fail,
            "cleared_heavy_shortfall": counter[f"{arm}_cleared_shortfall"],
            "new_heavy_shortfall": counter[f"{arm}_new_shortfall"],
            "descriptive_cell_sign_tail":
                binomial_one_sided(repairs, introductions),
            "paired_success_cells": paired_success,
            "paired_success_base_inact_milli_sum":
                counter[f"{arm}_paired_base_inact_milli"],
            "paired_success_candidate_inact_milli_sum":
                counter[f"{arm}_paired_candidate_inact_milli"],
            "inact_ratio_paired_success": ratio(
                counter[f"{arm}_paired_candidate_inact_milli"],
                counter[f"{arm}_paired_base_inact_milli"],
            ),
            "paired_success_base_block_xors_milli_sum":
                counter[f"{arm}_paired_base_xors_milli"],
            "paired_success_candidate_block_xors_milli_sum":
                counter[f"{arm}_paired_candidate_xors_milli"],
            "block_xor_ratio_paired_success": ratio(
                counter[f"{arm}_paired_candidate_xors_milli"],
                counter[f"{arm}_paired_base_xors_milli"],
            ),
            "paired_success_base_block_muladds_milli_sum":
                counter[f"{arm}_paired_base_muladds_milli"],
            "paired_success_candidate_block_muladds_milli_sum":
                counter[f"{arm}_paired_candidate_muladds_milli"],
            "block_muladd_ratio_paired_success": ratio(
                counter[f"{arm}_paired_candidate_muladds_milli"],
                counter[f"{arm}_paired_base_muladds_milli"],
            ),
            "failure_path_costs_excluded_from_ratios": True,
        }
    return {"cells": counter["cells"], "arms": arms, "comparisons": comparisons}


def analyze(
    result_dir: Path, contract: dict[str, Any], groups: Sequence[Group],
    jobs: Sequence[Job], phase_hashes: dict[Path, str],
) -> dict[str, Any]:
    result_root = result_dir / "results"
    expected_result_paths = {
        path.resolve(strict=True)
        for path in result_root.rglob("*") if path.is_file()
        and path.name != "phase_complete.sha256"
    }
    if set(phase_hashes) != expected_result_paths:
        die("result phase seal does not cover the exact result artifact set")
    identity = verify_low_identity(jobs, result_root, phase_hashes)
    data = load_arm_data(jobs, result_root, phase_hashes)
    group_by_k = [-1] * (K_MAX + 1)
    ledger_salt_by_k = [""] * (K_MAX + 1)
    for group in groups:
        for K in group.ks:
            group_by_k[K] = group.group
            ledger_salt_by_k[K] = group.ledger_salt

    analysis = result_dir / "analysis"
    analysis.mkdir(exist_ok=False)
    paired_path = analysis / "paired_cells.csv"
    failures_path = analysis / "failures_and_shortfalls.csv"
    exact_path = analysis / "exact_k_changes.csv"
    scopes_path = analysis / "scopes.csv"
    total = Counter()
    scopes: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    with paired_path.open("w", encoding="utf-8", newline="") as paired_stream, \
            failures_path.open("w", encoding="utf-8", newline="") as failure_stream:
        paired_writer = csv.writer(paired_stream, lineterminator="\n")
        paired_writer.writerow(PAIRED_HEADER)
        failure_writer = csv.writer(failure_stream, lineterminator="\n")
        failure_writer.writerow(PAIRED_HEADER)
        for seed_index, seed in enumerate(SEEDS):
            for schedule in SCHEDULES:
                for K in range(K_MIN, K_MAX + 1):
                    index = cell_index(seed_index, schedule, K)
                    values: list[Any] = [
                        K, group_by_k[K], ledger_salt_by_k[K], "0x0",
                        schedule, seed_index, seed,
                    ]
                    interesting = False
                    for arm in ARMS:
                        arm_data = data[arm]
                        arm_values = [
                            arm_data.rank[index], arm_data.error[index],
                            arm_data.shortfall[index], arm_data.seed_attempt[index],
                            arm_data.inact[index], arm_data.binary_def[index],
                            arm_data.heavy_gain[index], arm_data.xors[index],
                            arm_data.muladds[index],
                        ]
                        values.extend(arm_values)
                        interesting |= bool(
                            arm_data.rank[index] or arm_data.error[index] or
                            arm_data.shortfall[index]
                        )
                    paired_writer.writerow(values)
                    if interesting:
                        failure_writer.writerow(values)
                    update_counter(total, data, index)
                    for scope_type, scope in (
                        ("schedule", schedule),
                        ("seed", f"seed{seed_index}"),
                        ("schedule_seed", f"{schedule}:seed{seed_index}"),
                        ("K_band", band(K)),
                    ):
                        update_counter(scopes[(scope_type, scope)], data, index)

    exact_header = [
        "K", *(
            f"{arm}_{field}" for arm in ARMS
            for field in ("fail_strata", "error_strata", "shortfall_sum")
        ),
        "two_anchor_repairs", "two_anchor_introductions",
        "d13_repairs", "d13_introductions",
    ]
    consistency: dict[str, dict[str, Any]] = {
        arm: {
            "K_better": 0, "K_worse": 0, "K_equal": 0,
            "max_fail_strata": 0, "worse_K": [],
            "new_failure_resonance_K": [],
            "cleared_failure_resonance_K": [],
        }
        for arm in ARMS[1:]
    }
    exact_rows = 0
    with exact_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(exact_header)
        for K in range(K_MIN, K_MAX + 1):
            per_arm: dict[str, tuple[int, int, int]] = {}
            repairs: dict[str, int] = {}
            intros: dict[str, int] = {}
            for arm in ARMS:
                values = data[arm]
                indices = [
                    cell_index(seed_index, schedule, K)
                    for seed_index in range(len(SEEDS)) for schedule in SCHEDULES
                ]
                per_arm[arm] = (
                    sum(failed(values, index) for index in indices),
                    sum(values.error[index] for index in indices),
                    sum(values.shortfall[index] for index in indices),
                )
            base = data["d12"]
            for arm in ARMS[1:]:
                candidate = data[arm]
                indices = [
                    cell_index(seed_index, schedule, K)
                    for seed_index in range(len(SEEDS)) for schedule in SCHEDULES
                ]
                repairs[arm] = sum(
                    failed(base, index) and not failed(candidate, index)
                    for index in indices
                )
                intros[arm] = sum(
                    not failed(base, index) and failed(candidate, index)
                    for index in indices
                )
                base_count, candidate_count = per_arm["d12"][0], per_arm[arm][0]
                report = consistency[arm]
                relation = (
                    "K_better" if candidate_count < base_count else
                    "K_worse" if candidate_count > base_count else "K_equal"
                )
                report[relation] += 1
                report["max_fail_strata"] = max(
                    report["max_fail_strata"], candidate_count
                )
                if candidate_count > base_count:
                    report["worse_K"].append(K)
                if candidate_count >= 2 and base_count < 2:
                    report["new_failure_resonance_K"].append(K)
                if candidate_count < 2 and base_count >= 2:
                    report["cleared_failure_resonance_K"].append(K)
            changed = (
                len(set(per_arm.values())) > 1 or
                any(repairs.values()) or any(intros.values()) or
                any(value[0] >= 2 for value in per_arm.values())
            )
            if changed:
                writer.writerow([
                    K, *(metric for arm in ARMS for metric in per_arm[arm]),
                    repairs["two_anchor_adaptive"],
                    intros["two_anchor_adaptive"],
                    repairs["d13_adaptive"], intros["d13_adaptive"],
                ])
                exact_rows += 1

    for report in consistency.values():
        report["descriptive_K_cluster_sign_tail"] = binomial_one_sided(
            report["K_better"], report["K_worse"]
        )

    with scopes_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("scope_type", "scope", "summary_json"))
        for key in sorted(scopes):
            writer.writerow((*key, canonical_json(counter_report(scopes[key]))))

    summary = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.analysis.v3",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "coverage": {
            "K_min": K_MIN, "K_max": K_MAX, "unique_K": K_COUNT,
            "strata": len(SEEDS) * len(SCHEDULES),
            "cells_per_arm": K_COUNT * len(SEEDS) * len(SCHEDULES),
            "arms": ARMS, "total_recovery_cells": total["cells"] * len(ARMS),
        },
        "independence": contract["independence"],
        "statistical_interpretation": (
            "This is a deterministic census over the sealed seeds and loss "
            "schedules, not an IID population sample. Cell- and K-direction "
            "exact sign tails are descriptive diagnostics only; promotion is "
            "not gated on a p-value and must also use per-seed/schedule "
            "direction, effect size, and new-resonance checks."
        ),
        "cost_interpretation": (
            "Architecture cost ratios and paired sums use only cells where "
            "both arms succeeded. Per-arm all-cell, success, and failure-path "
            "sums are reported separately because rank-first failures may exit "
            "before projection or payload RHS work. Saturated campaign elapsed "
            "times are provenance only; speed claims require pinned full-payload "
            "timing."
        ),
        "adaptive_low_K_identity": identity,
        "overall": counter_report(total),
        "consistency": consistency,
        "failure_resonance_threshold_strata": 2,
        "exact_k_report_rows": exact_rows,
        "saturated_timing_speed_claim_valid": False,
        "phase_seal_sha256": sha256_file(result_root / "phase_complete.sha256"),
    }
    summary_path = analysis / "summary.json"
    atomic_json(summary_path, summary)
    analysis_files = [
        paired_path, failures_path, exact_path, scopes_path, summary_path,
    ]
    analysis_seal = analysis / "analysis_complete.sha256"
    atomic_write(analysis_seal, sha_manifest(result_dir, analysis_files))
    verified = verify_sha_manifest(result_dir, analysis_seal)
    if set(verified) != {path.resolve(strict=True) for path in analysis_files}:
        die("analysis seal does not cover the exact five analysis artifacts")
    return summary


def verify_analysis_publication(
    result_dir: Path,
    summary: dict[str, Any],
) -> dict[Path, str]:
    """Reopen the exact analysis inventory and bind its summary object."""
    analysis = result_dir / "analysis"
    artifact_names = {
        "paired_cells.csv", "failures_and_shortfalls.csv",
        "exact_k_changes.csv", "scopes.csv", "summary.json",
    }
    _require_exact_receipt_directory(
        analysis, artifact_names | {"analysis_complete.sha256"},
        "all-K analysis directory")
    seal = analysis / "analysis_complete.sha256"
    verified = verify_sha_manifest(result_dir, seal)
    expected_paths = {
        (analysis / name).resolve(strict=True) for name in artifact_names
    }
    if set(verified) != expected_paths:
        die("analysis seal does not cover the exact five analysis artifacts")
    summary_path = analysis / "summary.json"
    summary_data = stable_bytes(summary_path)
    if (
        summary_data != json_bytes(summary) or
        verified.get(summary_path.resolve(strict=True)) !=
            sha256_bytes(summary_data)
    ):
        die("analysis summary differs from its sealed in-memory result")
    return verified


def run(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    contract, binary, taskset, thermal, thermal_identity, groups = \
        load_frozen(result_dir)
    jobs = build_jobs(groups)
    result_root = result_dir / "results"
    result_root.mkdir(exist_ok=False)
    for name in ("stdout", "stderr", "commands"):
        (result_root / name).mkdir()
    cpus = sorted(os.sched_getaffinity(0))
    workers = min(int(contract["workers"]), len(cpus), len(jobs))
    if workers <= 0:
        die("campaign has no available worker CPU")
    pool = CpuPool(cpus, workers)
    abort = threading.Event()
    guard = ThermalGuard(
        thermal, abort,
        stale_seconds=float(contract["thermal_policy"]["stale_seconds"]),
        limit_c=float(contract["thermal_policy"]["limit_c"]),
        consecutive_limit=int(contract["thermal_policy"]["consecutive_samples"]),
        expected_dev=thermal_identity[0], expected_ino=thermal_identity[1],
        expected_edac_ce=thermal_identity[2],
        expected_edac_ue=thermal_identity[3],
    )
    if validate_frozen_thermal_source(
            thermal, contract.get("thermal_baseline"),
            float(contract["thermal_policy"]["stale_seconds"]),
            cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
            dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
            consecutive_limit=int(
                contract["thermal_policy"]["consecutive_samples"])) != \
            thermal_identity:
        die("frozen thermal baseline changed before all-K guard start")
    registry = ProcessRegistry()
    records: list[dict[str, Any]] = []
    campaign_error: BaseException | None = None
    thermal_summary: dict[str, Any] | None = None
    with CampaignSignalGuard(abort, registry):
        try:
            guard.start()
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=workers
            ) as executor:
                futures: list[concurrent.futures.Future[dict[str, Any]]] = []
                completed = 0
                try:
                    for job in jobs:
                        if abort.is_set():
                            die("campaign aborted during job submission")
                        futures.append(executor.submit(
                            run_job, job, binary, contract["binary_sha256"],
                            taskset, contract["taskset"], result_root, pool,
                            abort, registry,
                            float(contract["timeout_seconds"]),
                        ))
                    for future in concurrent.futures.as_completed(futures):
                        records.append(future.result())
                        completed += 1
                        if completed % 100 == 0 or completed == len(jobs):
                            print(
                                f"all-K progress {completed}/{len(jobs)}",
                                flush=True,
                            )
                        if guard.error:
                            die(f"thermal guard failed: {guard.error}")
                except BaseException as exc:
                    campaign_error = exc
                    abort.set()
                    registry.signal_all(signal.SIGTERM)
                    for future in futures:
                        future.cancel()
        finally:
            if guard.started:
                coverage_end = (
                    max(record["end_monotonic_ns"] for record in records)
                    if campaign_error is None and len(records) == len(jobs)
                    else None)
                thermal_summary = guard.finish(
                    result_root / "thermal_interval.csv",
                    cover_through_monotonic_ns=coverage_end)
    if guard.error:
        die(f"thermal guard failed: {guard.error}")
    if thermal_summary is None:
        if campaign_error is not None:
            raise campaign_error
        die("thermal guard did not produce a summary")
    if (
        thermal_summary["thermal_high_max_consecutive_samples"] >=
            int(contract["thermal_policy"]["consecutive_samples"]) or
        thermal_summary["dimm_read_errors_max"] != 0 or
        thermal_summary["edac_ce_delta"] != 0 or
        thermal_summary["edac_ue_delta"] != 0
    ):
        die(f"memory telemetry changed during campaign: {thermal_summary}")
    if campaign_error is not None:
        raise campaign_error
    if thermal_summary["cpu_busy_min_pct"] < \
            float(contract["thermal_policy"]["min_cpu_busy_pct"]):
        die(f"CPU utilization fell below campaign gate: {thermal_summary}")
    verify_frozen_binary(binary, contract.get("binary_sha256"))
    if frozen_executable_path(contract.get("taskset"), "taskset") != taskset:
        die("frozen taskset path changed during the all-K campaign")
    if validate_frozen_thermal_source(
            thermal, contract.get("thermal_baseline"),
            float(contract["thermal_policy"]["stale_seconds"]),
            cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
            dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
            consecutive_limit=int(
                contract["thermal_policy"]["consecutive_samples"])) != \
            thermal_identity:
        die("frozen thermal baseline changed during the all-K campaign")
    records.sort(key=lambda record: record["job"])
    if len(records) != len(jobs):
        die(f"completed {len(records)} jobs, want {len(jobs)}")
    tasks_data = b"".join(
        (canonical_json(record) + "\n").encode("utf-8")
        for record in records
    )
    atomic_write(result_root / "tasks.jsonl", tasks_data)
    run_summary = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.run.v5",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "jobs": len(jobs), "workers": workers,
        "worker_cpus": cpus[:workers],
        "start_ns": min(record["start_ns"] for record in records),
        "end_ns": max(record["end_ns"] for record in records),
        "start_monotonic_ns": min(
            record["start_monotonic_ns"] for record in records),
        "end_monotonic_ns": max(
            record["end_monotonic_ns"] for record in records),
        "recovery_cells": sum(len(job.ks) for job in jobs),
        "thermal": thermal_summary,
        "saturated_timing_speed_claim_valid": False,
    }
    atomic_json(result_root / "run.json", run_summary)
    receipt_hashes = verify_allk_receipts(
        result_dir, contract, jobs, binary, taskset, workers,
        cpus[:workers], thermal_summary,
        sha256_bytes(guard.mark["baseline_row"]),
    )
    phase_files = list(receipt_hashes)
    phase_seal = result_root / "phase_complete.sha256"
    atomic_write(phase_seal, sha_manifest(result_dir, phase_files))
    phase_hashes = verify_sha_manifest(result_dir, phase_seal)
    if phase_hashes != receipt_hashes:
        die("all-K receipts changed between semantic replay and phase seal")
    summary = analyze(result_dir, contract, groups, jobs, phase_hashes)
    if verify_sha_manifest(result_dir, phase_seal) != receipt_hashes:
        die("all-K phase artifacts changed during analysis")
    analysis_hashes = verify_analysis_publication(result_dir, summary)
    phase_seal_sha256 = sha256_file(phase_seal)
    analysis_seal = result_dir / "analysis/analysis_complete.sha256"
    analysis_seal_sha256 = sha256_file(analysis_seal)
    if summary.get("phase_seal_sha256") != phase_seal_sha256:
        die("analysis summary does not bind the final phase seal")
    # Recheck after hashing both manifests so a concurrent replacement cannot
    # splice the in-memory summary to a different final complete record.
    if (
        verify_sha_manifest(result_dir, phase_seal) != receipt_hashes or
        verify_analysis_publication(result_dir, summary) != analysis_hashes or
        sha256_file(phase_seal) != phase_seal_sha256 or
        sha256_file(analysis_seal) != analysis_seal_sha256
    ):
        die("all-K evidence changed before final publication")
    complete = {
        "result_dir": str(result_dir),
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "staged_seal_sha256": sha256_file(result_dir / "frozen/staged.sha256"),
        "phase_seal_sha256": phase_seal_sha256,
        "analysis_seal_sha256": analysis_seal_sha256,
        "summary": summary,
    }
    complete_path = result_dir / "complete.json"
    atomic_json(complete_path, complete)
    if stable_bytes(complete_path) != json_bytes(complete):
        die("final all-K publication record changed after creation")
    print(canonical_json(complete))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument(
        "--binary", type=Path, required=True,
        help=("accepted local benchmark used only to confirm the system "
              "toolchain; its binary and generated graph are never reused"))
    prepare_parser.add_argument("--groups", type=Path, required=True)
    prepare_parser.add_argument("--thermal", type=Path, required=True)
    prepare_parser.add_argument("--result-dir", type=Path, required=True)
    prepare_parser.add_argument("--workers", type=int, default=120)
    prepare_parser.add_argument("--build-workers", type=int, default=32)
    prepare_parser.add_argument("--timeout", type=float, default=1200.0)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    getcontext().prec = 100
    if args.mode == "prepare":
        if (args.workers <= 0 or args.build_workers <= 0 or
                not math.isfinite(args.timeout) or args.timeout <= 0):
            die("--workers, --build-workers, and --timeout must be positive")
        return prepare(args)
    return run(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"rank-floor two-anchor all-K: {exc}", file=sys.stderr)
        raise
