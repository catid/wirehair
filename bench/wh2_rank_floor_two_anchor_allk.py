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
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, getcontext
import errno
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import select
import shutil
import signal
import stat
import subprocess
import threading
import time
from typing import Any, Iterable, Sequence

from wh2_rank_floor_two_anchor_screen import (
    THERMAL_FIELDS,
    TIMING_FIELDS,
    binomial_one_sided,
    canonical_json,
    deterministic_equal,
    parse_bench_output,
    parse_thermal_sample,
    sha256_file as _screen_sha256_file,
    thermal_finish,
    thermal_start,
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


def open_durable_directory(path: Path, *, create: bool = False) -> int:
    """Open an absolute directory without following any component symlink."""
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
            created = False
            if create:
                try:
                    os.mkdir(part, 0o700, dir_fd=descriptor)
                    created = True
                except FileExistsError:
                    pass
            next_descriptor = os.open(part, flags, dir_fd=descriptor)
            if created:
                os.fsync(next_descriptor)
                os.fsync(descriptor)
            os.close(descriptor)
            descriptor = next_descriptor
        metadata = os.fstat(descriptor)
        named = os.stat(absolute, follow_symlinks=False)
        if (not stat.S_ISDIR(metadata.st_mode) or
                (metadata.st_dev, metadata.st_ino) !=
                (named.st_dev, named.st_ino)):
            die(f"directory pathname changed while opening: {path}")
        return descriptor
    except OSError as error:
        os.close(descriptor)
        raise CampaignError(
            f"cannot open durable directory {path}: {error}") from error
    except BaseException:
        os.close(descriptor)
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


def stop_and_reap_process_group(
    process: subprocess.Popen[bytes], grace_seconds: float = 5.0,
) -> None:
    try:
        signal_process_group(process, signal.SIGTERM)
        try:
            process.communicate(timeout=grace_seconds)
        except BaseException:
            # Cleanup must continue even if pipe communication itself fails.
            pass
        if process.poll() is None:
            signal_process_group(process, signal.SIGKILL)
            try:
                process.wait(timeout=grace_seconds)
            except BaseException:
                pass
        deadline = time.monotonic() + grace_seconds
        while process_group_exists(process) and time.monotonic() < deadline:
            time.sleep(0.05)
        if process_group_exists(process):
            signal_process_group(process, signal.SIGKILL)
            kill_deadline = time.monotonic() + grace_seconds
            while (process_group_exists(process) and
                    time.monotonic() < kill_deadline):
                time.sleep(0.05)
    finally:
        close_process_streams(process)
    if process.poll() is None or process_group_exists(process):
        raise CampaignError(f"child process group {process.pid} was not reaped")


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


class ThermalGuard:
    def __init__(
        self, path: Path, abort: threading.Event,
        stale_seconds: float = 5.0, limit_c: float = 90.0,
        consecutive_limit: int = 3,
    ) -> None:
        self.path = path
        self.abort = abort
        self.stale_seconds = stale_seconds
        self.limit_c = limit_c
        self.consecutive_limit = consecutive_limit
        if (not math.isfinite(limit_c) or not 0.0 <= limit_c <= 120.0 or
                consecutive_limit <= 0):
            die("thermal guard policy is outside the canonical range")
        self.mark = thermal_start(path, stale_seconds=stale_seconds)
        self.read_offset = self.mark["offset"]
        self.last_monotonic = self.mark["monotonic_s"]
        self.stop_event = threading.Event()
        self.error: str | None = None
        self.poll_iterations = 0
        self.samples = 0
        self.high_samples = 0
        self.thread = threading.Thread(target=self._loop, name="thermal-guard")
        self.started = False

    def start(self) -> None:
        try:
            self.thread.start()
        finally:
            self.started = self.thread.ident is not None

    def _new_samples(self) -> list[dict[str, Any]]:
        with self.path.open("rb") as source:
            file_stat = os.fstat(source.fileno())
            if ((file_stat.st_dev, file_stat.st_ino) !=
                    (self.mark["dev"], self.mark["ino"])):
                die("thermal logger inode changed")
            if file_stat.st_size < self.read_offset:
                die("thermal logger shrank")
            source.seek(self.read_offset)
            suffix = source.read()
        complete_end = suffix.rfind(b"\n")
        if complete_end < 0:
            if len(suffix) > 16384:
                die("thermal logger has an oversized partial row")
            return []
        complete = suffix[:complete_end + 1]
        partial = suffix[complete_end + 1:]
        if len(partial) > 16384:
            die("thermal logger has an oversized partial row")
        lines = complete.splitlines(keepends=True)
        samples = [
            parse_thermal_sample(line, f"live thermal row {number}")
            for number, line in enumerate(lines, 1)
        ]
        self.read_offset += len(complete)
        return samples

    def _loop(self) -> None:
        consecutive = 1 if self.mark["max_temperature_c"] >= self.limit_c else 0
        while not self.stop_event.wait(0.25):
            try:
                samples = self._new_samples()
                self.poll_iterations += 1
                for sample in samples:
                    delta = sample["monotonic_s"] - self.last_monotonic
                    if delta <= 0.0 or delta > self.stale_seconds:
                        die("thermal logger is nonmonotonic or has a sampling gap")
                    validate_thermal_health(sample, "live thermal sample")
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

    def finish(self, output: Path) -> dict[str, Any]:
        if not self.started:
            die("thermal guard was never started")
        self.stop_event.set()
        self.thread.join()
        summary = thermal_finish(
            self.path, self.mark, output,
            limit_c=self.limit_c, stale_seconds=self.stale_seconds,
            dimm_limit_c=self.limit_c,
        )
        summary.update({
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
        command = ["taskset", "-c", str(cpu), *make_command(binary, job)]
        start_ns = time.time_ns()
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
        registered = False
        try:
            registry.add(process)
            registered = True
            deadline = time.monotonic() + timeout
            while True:
                remaining = deadline - time.monotonic()
                try:
                    stdout_bytes, stderr_bytes = process.communicate(
                        timeout=max(0.001, min(0.25, remaining))
                    )
                    break
                except subprocess.TimeoutExpired:
                    reason = None
                    if abort.is_set():
                        reason = "campaign abort"
                    elif remaining <= 0:
                        reason = f"timeout after {timeout:g}s"
                    if reason is None:
                        continue
                    stop_and_reap_process_group(process)
                    die(f"job {job.job} terminated on {reason}")
        finally:
            leader_live = process.poll() is None
            descendants_live = process_group_exists(process)
            if leader_live or descendants_live:
                stop_and_reap_process_group(process)
            else:
                close_process_streams(process)
            if process.poll() is None or process_group_exists(process):
                die(f"job {job.job} child cleanup was not proven")
            if registered:
                registry.remove(process)
            if not leader_live and descendants_live:
                die(f"job {job.job} left an unexpected child process")
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
        "elapsed_ns": end_ns - start_ns, "returncode": process.returncode,
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


def tracked_clean_source(script: Path, helper: Path) -> tuple[Path, str]:
    repo = Path(subprocess.check_output(
        ("git", "-C", str(script.parent), "rev-parse", "--show-toplevel"),
        text=True,
    ).strip()).resolve(strict=True)
    for args, name in (
        (("diff", "--quiet", "HEAD", "--"), "worktree"),
        (("diff", "--cached", "--quiet", "HEAD", "--"), "index"),
    ):
        if subprocess.run(("git", "-C", str(repo), *args), check=False).returncode:
            die(f"tracked source {name} is dirty")
    head = subprocess.check_output(
        ("git", "-C", str(repo), "rev-parse", "HEAD"), text=True,
    ).strip()
    for path in (script, helper):
        relative = path.resolve(strict=True).relative_to(repo)
        tracked = subprocess.check_output(
            ("git", "-C", str(repo), "show", f"HEAD:{relative}"),
        )
        if tracked != stable_bytes(path):
            die(f"{relative} does not exactly match immutable HEAD")
    return repo, head


def cmake_source_directory(cache: Path) -> Path:
    prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
    matches = [
        line[len(prefix):]
        for line in stable_bytes(cache).decode("utf-8").splitlines()
        if line.startswith(prefix)
    ]
    if len(matches) != 1:
        die("CMake cache lacks one canonical source-directory record")
    return Path(matches[0]).resolve(strict=True)


def prepare(args: argparse.Namespace) -> int:
    script = Path(__file__).resolve(strict=True)
    helper = script.with_name("wh2_rank_floor_two_anchor_screen.py")
    repo, head = tracked_clean_source(script, helper)
    binary_input = args.binary.absolute()
    if binary_input.is_symlink():
        die("--binary must not be a symlink")
    binary = binary_input.resolve(strict=True)
    build_dir = binary.parent.parent.resolve(strict=True)
    cache = build_dir / "CMakeCache.txt"
    if binary != (build_dir / "codec/wirehair_v2_bench").resolve(strict=True):
        die("--binary must be the codec/wirehair_v2_bench build target")
    if cmake_source_directory(cache) != repo:
        die("--binary CMake cache belongs to a different source tree")
    build_command = [
        "cmake", "--build", str(build_dir), "--target",
        "wirehair_v2_bench", "--clean-first", "--parallel",
        str(args.build_workers),
    ]
    subprocess.run(build_command, check=True)
    post_repo, post_head = tracked_clean_source(script, helper)
    if post_repo != repo or post_head != head:
        die("source repository changed during the clean benchmark rebuild")
    if binary_input.is_symlink():
        die("clean benchmark rebuild replaced --binary with a symlink")
    binary = binary_input.resolve(strict=True)
    groups_input = args.groups.absolute()
    if groups_input.is_symlink():
        die("--groups must not be a symlink")
    groups = groups_input.resolve(strict=True)
    thermal = args.thermal.resolve(strict=True)
    if binary.is_symlink() or not binary.is_file() or not os.access(binary, os.X_OK):
        die("--binary must be an executable, non-symlink regular file")
    load_groups(groups)
    thermal_start(thermal)
    result_dir = args.result_dir.resolve()
    result_dir.mkdir(parents=True, exist_ok=False)
    frozen = result_dir / "frozen"
    frozen.mkdir()
    destinations = {
        "script": frozen / script.name,
        "helper": frozen / helper.name,
        "binary": frozen / "wirehair_v2_bench",
        "cmake_cache": frozen / "CMakeCache.txt",
        "groups": frozen / "groups.tsv",
    }
    for source, destination in (
        (script, destinations["script"]),
        (helper, destinations["helper"]),
        (binary, destinations["binary"]),
        (cache, destinations["cmake_cache"]),
        (groups, destinations["groups"]),
    ):
        shutil.copyfile(source, destination)
    destinations["binary"].chmod(0o755)
    destinations["script"].chmod(0o755)
    contract = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.contract.v4",
        "source_commit": head, "source_repo": str(repo),
        "build_command": build_command,
        "binary_source": str(binary),
        "binary_sha256": sha256_file(destinations["binary"]),
        "cmake_cache_sha256": sha256_file(destinations["cmake_cache"]),
        "script_sha256": sha256_file(destinations["script"]),
        "helper_sha256": sha256_file(destinations["helper"]),
        "groups_sha256": sha256_file(destinations["groups"]),
        "thermal": str(thermal), "workers": args.workers,
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
        "result_dir": str(result_dir), "source_commit": head,
        "binary_sha256": contract["binary_sha256"],
        "staged_sha256": sha256_file(staged_path),
        "run_command": [
            sys.executable, str(destinations["script"]), "run",
            "--result-dir", str(result_dir),
        ],
    }
    atomic_json(result_dir / "prepare.json", prepared)
    print(canonical_json(prepared))
    return 0


def load_frozen(result_dir: Path) -> tuple[dict[str, Any], Path, Path, list[Group]]:
    frozen = result_dir / "frozen"
    staged = verify_sha_manifest(result_dir, frozen / "staged.sha256")
    expected = {
        (frozen / name).resolve(strict=True)
        for name in (
            "wh2_rank_floor_two_anchor_allk.py",
            "wh2_rank_floor_two_anchor_screen.py",
            "wirehair_v2_bench", "CMakeCache.txt", "groups.tsv",
            "contract.json",
        )
    }
    if set(staged) != expected:
        die("frozen seal does not cover the exact six-file staged set")
    contract = json.loads(stable_bytes(frozen / "contract.json"))
    if (
        contract.get("schema") !=
            "wirehair.wh2.rank_floor_two_anchor_allk.contract.v4"
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
    ):
        die("frozen contract does not match this campaign")
    script = (frozen / "wh2_rank_floor_two_anchor_allk.py").resolve(strict=True)
    helper = (frozen / "wh2_rank_floor_two_anchor_screen.py").resolve(strict=True)
    if Path(__file__).resolve(strict=True) != script:
        die("run must be invoked through the frozen campaign script")
    if (
        sha256_file(script) != contract["script_sha256"] or
        sha256_file(helper) != contract["helper_sha256"]
    ):
        die("frozen script or helper does not match the contract")
    binary = (frozen / "wirehair_v2_bench").resolve(strict=True)
    if not os.access(binary, os.X_OK) or sha256_file(binary) != contract["binary_sha256"]:
        die("frozen binary is missing, non-executable, or changed")
    groups = load_groups(frozen / "groups.tsv", contract["groups_sha256"])
    thermal = Path(contract["thermal"]).resolve(strict=True)
    return contract, binary, thermal, groups


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


def load_arm_data(
    jobs: Sequence[Job], result_root: Path,
) -> dict[str, ArmData]:
    data = {arm: ArmData() for arm in ARMS}
    for ordinal, job in enumerate(jobs, 1):
        stdout_path = result_root / "stdout" / f"{job.stem}.csv"
        text = stable_bytes(stdout_path).decode("utf-8")
        rows = parse_bench_output(
            text, job.ks, "0x0", job.seed, job.schedule,
            arm_options(job.arm, job.band), True,
        )
        arm = data[job.arm]
        for K, row in zip(job.ks, rows, strict=True):
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
    jobs: Sequence[Job], result_root: Path,
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
                stable_bytes(path).decode("utf-8"), job.ks, "0x0",
                job.seed, job.schedule, (), True,
            )
        base = parsed["d12"]
        for candidate in ARMS[1:]:
            for base_row, candidate_row in zip(base, parsed[candidate], strict=True):
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
        counter[f"{arm}_fail"] += failed(values, index)
        counter[f"{arm}_error"] += values.error[index]
        counter[f"{arm}_shortfall"] += values.shortfall[index]
        counter[f"{arm}_inact_milli"] += values.inact[index]
        counter[f"{arm}_binary_def_milli"] += values.binary_def[index]
        counter[f"{arm}_heavy_gain_milli"] += values.heavy_gain[index]
        counter[f"{arm}_binary_def_gt15"] += values.binary_def[index] > 15000
        counter[f"{arm}_binary_def_max_milli"] = max(
            counter[f"{arm}_binary_def_max_milli"], values.binary_def[index]
        )
        counter[f"{arm}_xors_milli"] += values.xors[index]
        counter[f"{arm}_muladds_milli"] += values.muladds[index]
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
            counter[f"{arm}_paired_base_xors_milli"] += base.xors[index]
            counter[f"{arm}_paired_candidate_xors_milli"] += candidate.xors[index]
            counter[f"{arm}_paired_base_muladds_milli"] += base.muladds[index]
            counter[f"{arm}_paired_candidate_muladds_milli"] += candidate.muladds[index]


def counter_report(counter: Counter[str]) -> dict[str, Any]:
    arms: dict[str, Any] = {}
    for arm in ARMS:
        arms[arm] = {
            "failures": counter[f"{arm}_fail"],
            "errors": counter[f"{arm}_error"],
            "heavy_shortfall_sum": counter[f"{arm}_shortfall"],
            "inact_milli_sum": counter[f"{arm}_inact_milli"],
            "binary_def_milli_sum": counter[f"{arm}_binary_def_milli"],
            "binary_def_max_milli": counter[f"{arm}_binary_def_max_milli"],
            "binary_def_gt15_cells": counter[f"{arm}_binary_def_gt15"],
            "heavy_gain_milli_sum": counter[f"{arm}_heavy_gain_milli"],
            "block_xors_milli_sum": counter[f"{arm}_xors_milli"],
            "block_muladds_milli_sum": counter[f"{arm}_muladds_milli"],
        }
    comparisons: dict[str, Any] = {}
    for arm in ARMS[1:]:
        repairs = counter[f"{arm}_repair"]
        introductions = counter[f"{arm}_intro"]
        comparisons[f"{arm}_vs_d12"] = {
            "repairs": repairs, "introductions": introductions,
            "both_fail": counter[f"{arm}_both_fail"],
            "cleared_heavy_shortfall": counter[f"{arm}_cleared_shortfall"],
            "new_heavy_shortfall": counter[f"{arm}_new_shortfall"],
            "descriptive_cell_sign_tail":
                binomial_one_sided(repairs, introductions),
            "paired_success_cells": counter[f"{arm}_paired_success"],
            "block_xor_ratio_paired_success": ratio(
                counter[f"{arm}_paired_candidate_xors_milli"],
                counter[f"{arm}_paired_base_xors_milli"],
            ),
            "block_muladd_ratio_paired_success": ratio(
                counter[f"{arm}_paired_candidate_muladds_milli"],
                counter[f"{arm}_paired_base_muladds_milli"],
            ),
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
    identity = verify_low_identity(jobs, result_root)
    data = load_arm_data(jobs, result_root)
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
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.analysis.v2",
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


def run(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    contract, binary, thermal, groups = load_frozen(result_dir)
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
    )
    registry = ProcessRegistry()
    records: list[dict[str, Any]] = []
    campaign_error: BaseException | None = None
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
                            run_job, job, binary, result_root, pool, abort,
                            registry, float(contract["timeout_seconds"]),
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
                thermal_summary = guard.finish(
                    result_root / "thermal_interval.csv"
                )
    if guard.error:
        die(f"thermal guard failed: {guard.error}")
    if (
        thermal_summary["cpu_busy_min_pct"] <
            float(contract["thermal_policy"]["min_cpu_busy_pct"]) or
        thermal_summary["thermal_high_max_consecutive_samples"] >=
            int(contract["thermal_policy"]["consecutive_samples"]) or
        thermal_summary["dimm_read_errors_max"] != 0 or
        thermal_summary["edac_ce_delta"] != 0 or
        thermal_summary["edac_ue_delta"] != 0
    ):
        die(f"memory telemetry changed during campaign: {thermal_summary}")
    if campaign_error is not None:
        raise campaign_error
    records.sort(key=lambda record: record["job"])
    if len(records) != len(jobs):
        die(f"completed {len(records)} jobs, want {len(jobs)}")
    with (result_root / "tasks.jsonl").open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        for record in records:
            output.write(canonical_json(record) + "\n")
    run_summary = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_allk.run.v3",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "jobs": len(jobs), "workers": workers,
        "recovery_cells": len(ARMS) * len(SEEDS) * len(SCHEDULES) * K_COUNT,
        "thermal": thermal_summary,
        "saturated_timing_speed_claim_valid": False,
    }
    atomic_json(result_root / "run.json", run_summary)
    phase_files = [
        path for path in result_root.rglob("*")
        if path.is_file() and path.name != "phase_complete.sha256"
    ]
    phase_seal = result_root / "phase_complete.sha256"
    atomic_write(phase_seal, sha_manifest(result_dir, phase_files))
    phase_hashes = verify_sha_manifest(result_dir, phase_seal)
    summary = analyze(result_dir, contract, groups, jobs, phase_hashes)
    complete = {
        "result_dir": str(result_dir),
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "staged_seal_sha256": sha256_file(result_dir / "frozen/staged.sha256"),
        "phase_seal_sha256": sha256_file(phase_seal),
        "analysis_seal_sha256": sha256_file(
            result_dir / "analysis/analysis_complete.sha256"
        ),
        "summary": summary,
    }
    atomic_json(result_dir / "complete.json", complete)
    print(canonical_json(complete))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--binary", type=Path, required=True)
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
        if args.workers <= 0 or args.build_workers <= 0 or args.timeout <= 0:
            die("--workers, --build-workers, and --timeout must be positive")
        return prepare(args)
    return run(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"rank-floor two-anchor all-K: {exc}", file=sys.stderr)
        raise
