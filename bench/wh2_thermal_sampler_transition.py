#!/usr/bin/env python3
"""Fail-closed, bounded handoff around the WH2 thermal sampler.

This controller is deliberately specific to the reviewed g8iv.6.22.1 live
dry run.  It stops only the externally sealed old sampler identity, runs the
reviewed candidate through the production P32 ownership/evidence helpers, and
restores the old sampler on every path after the stop has been armed.  It does
not promote the candidate to a long-running sampler.

The live entry point accepts two conspicuous confirmation tokens.  Importing
the module and running its unit tests cannot inspect or mutate live hardware.
"""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import signal
import stat
import subprocess
import sys
import time
from types import ModuleType
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

TRANSITION_ID = "wirehair-g8iv.6.22.1-live-dry-3c291790-v1"
EXECUTE_CONFIRMATION = "STOP-EXACT-OLD-SAMPLER-AND-RESTORE-IT"
TRANSITION_SCHEMA = "wirehair.wh2.thermal_transition.v3"
PHASE_SCHEMA = "wirehair.wh2.thermal_transition.phase.v3"
AUDIT_BINDING_SCHEMA = "wirehair.wh2.thermal_transition.audit_binding.v3"
PLAN_SCHEMA = "wirehair.wh2.thermal_transition.plan.v3"
RENAME_NOREPLACE = 1
SHA256_RE = re.compile(r"[0-9a-f]{64}")

CANDIDATE_SHA256 = (
    "3c291790b6169c07e53dc0924d383455a81addd6da858d8d4d231b9f1263f6a2")
P32_SHA256 = (
    "414806a687d184fc2ec2973fcd0ea4b8ab232d13d8e707404e926d394fefb821")
OLD_SOURCE_SHA256 = (
    "2b84efa91375a96a4a64e09ce5bfd7cba0b85b75028f5a93470cd4ae58aadb01")
OLD_BOOT_ID = "1788608a-7aa1-48de-8f7c-848485be3cc3"
OLD_CMDLINE_SHA256 = (
    "9626d8e11e0e49a6de6bf81c3f3d227f108540a8eb61f40f212c58f0da9fd069")
OLD_LAUNCHER_CMDLINE_SHA256 = (
    "b5cbaa94f51dfe477f0c600e27acdfc0e396b9ec46ab24098bf2a27e1df7c86f")

REPO_CANDIDATE = Path(
    "/tmp/wirehair-wh2-thermal-plausibility-v2/bench/"
    "wirehair_expo_thermal_sampler.py")
REPO_P32 = Path(
    "/tmp/wirehair-wh2-thermal-plausibility-v2/bench/"
    "wh2_p32_dispatch_timing.py")
DRY_ROOT = Path("/dev/shm") / TRANSITION_ID
OLD_ROOT = Path("/dev/shm/wh2-rhs-fusion-batch16-formal-v3")
OLD_SOURCE = OLD_ROOT / "campaign/frozen/wirehair_expo_thermal_sampler.py"
OLD_THERMAL_DIR = OLD_ROOT / "thermal"
OLD_CSV = OLD_THERMAL_DIR / "thermal.csv"
OLD_PID_FILE = OLD_THERMAL_DIR / "sampler.pid"
OLD_ARCHIVE = OLD_THERMAL_DIR / (
    "thermal.pre-g8iv.6.22.1.3c291790.csv")
OLD_UNCLEAN_ARCHIVE = OLD_THERMAL_DIR / (
    "thermal.pre-g8iv.6.22.1.3c291790.unclean.csv")
OLD_STALE_PID_ARCHIVE = OLD_THERMAL_DIR / (
    "sampler.pre-g8iv.6.22.1.3c291790.unclean.pid")
AUDIT_BINDING = OLD_THERMAL_DIR / (
    "future-audit-binding.g8iv.6.22.1.3c291790.json")

OLD_ARGV = (
    "/usr/bin/python3.12", "-I", str(OLD_SOURCE),
    "--csv", str(OLD_CSV), "--pid-file", str(OLD_PID_FILE),
)
OLD_LAUNCH_ARGV = (
    "sudo", "-n", "-b", "/usr/bin/taskset", "--cpu-list", "127",
    *OLD_ARGV,
)
OLD_REPLACEMENT_ARGV = (
    "/usr/bin/python3.12", "-I", "-S", "-B", str(OLD_SOURCE),
    "--csv", str(OLD_CSV), "--pid-file", str(OLD_PID_FILE),
)
OLD_REPLACEMENT_CMDLINE_SHA256 = (
    "a70c9f4bb6e9c13546682b95b13bc4fe8b3f76c1bae1ef01525fcf9f1b9f11e8")

EXECUTE_FLAG_CONTRACT: Dict[str, int | bool] = {
    "bytes_warning": 0,
    "debug": 0,
    "dev_mode": False,
    "dont_write_bytecode": 1,
    "hash_randomization": 1,
    "ignore_environment": 1,
    "inspect": 0,
    "int_max_str_digits": 4300,
    "interactive": 0,
    "isolated": 1,
    "no_site": 1,
    "no_user_site": 1,
    "optimize": 0,
    "quiet": 0,
    "safe_path": True,
    "utf8_mode": 1,
    "verbose": 0,
    "warn_default_encoding": 0,
}

# Exact host executables authorized for the one-off transition.  The plan also
# seals their inode metadata at preparation; the hardcoded digests make a
# same-path replacement fail even before the old sampler can be signalled.
TOOL_CONTRACT: Tuple[Tuple[str, str, str], ...] = (
    ("env", "/usr/bin/env",
     "0aefff8f912fb75716c5d4de3b6acde93edbe8fa280fc8ee895c1226d3e373ef"),
    ("fuser", "/usr/bin/fuser",
     "88aeee250ef0622fb638acd4694a4ffd96ec04d5e25dd2ca880739c4719197f5"),
    ("mpstat", "/usr/bin/mpstat",
     "74695dd7d010730cd1e19efa0664904c1e3b17d746d304c06243e183a5ac9f9c"),
    ("python3", "/usr/bin/python3.12",
     "1643dacd9feaedc58f3cc581e4d22577dfe25c09b10282936186ccf0f2e61118"),
    ("sudo", "/usr/bin/sudo",
     "136f2e48b0295b9fc595b8259cf2411ac43f27ddbfe02b956649ddaa2e92b9fa"),
    ("taskset", "/usr/bin/taskset",
     "a9c851792e54e91fba7b827019380abee54e715b6817899c835e4f221354b260"),
    ("timeout", "/usr/bin/timeout",
     "4fccd5b0192653a2446b745d5385ea547b78e466150e07ade9e2caff2b7f4e08"),
)


class TransitionError(RuntimeError):
    """A substantive transition or evidence failure."""


class TransitionDeadline(TransitionError):
    """The bounded transition exhausted its normal or recovery time."""


def renameat2_function() -> Callable[..., int]:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise TransitionError("renameat2 is unavailable")
    renameat2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint)
    renameat2.restype = ctypes.c_int
    return renameat2


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_fd(fd: int) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    os.lseek(fd, 0, os.SEEK_SET)
    while True:
        block = os.read(fd, 1024 * 1024)
        if not block:
            break
        digest.update(block)
        size += len(block)
    os.lseek(fd, 0, os.SEEK_SET)
    return digest.hexdigest(), size


def stable_file_bytes(path: Path, attempts: int = 20) -> bytes:
    for _attempt in range(attempts):
        before = path.stat()
        raw = path.read_bytes()
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size,
                before.st_mtime_ns) == (
                after.st_dev, after.st_ino, after.st_size,
                after.st_mtime_ns) and len(raw) == after.st_size:
            return raw
        time.sleep(0.02)
    raise TransitionError("file did not become stable: %s" % path)


def sha256_file(path: Path) -> str:
    return sha256_bytes(stable_file_bytes(path))


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if not schema or "schema" in payload or \
            "self_sha256_excluding_field" in payload:
        raise TransitionError("invalid sealed record payload")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        raise TransitionError(name + " schema mismatch")
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA256_RE.fullmatch(claimed) is None:
        raise TransitionError(name + " self hash malformed")
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        raise TransitionError(name + " self hash mismatch")
    return value


def load_canonical(path: Path, schema: str, name: str) -> Dict[str, object]:
    raw = stable_file_bytes(path)
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TransitionError(name + " is not canonical JSON") from exc
    if not isinstance(value, dict) or canonical_json(value) != raw:
        raise TransitionError(name + " is not canonical JSON")
    return verify_sealed(value, schema, name)


def fsync_directory(path: Path) -> None:
    fd = os.open(
        str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def file_binding(path: Path, *, with_hash: bool) -> Dict[str, object]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | \
        getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(path), flags)
    try:
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode):
            raise TransitionError("bound path is not a regular file: %s" % path)
        value: Dict[str, object] = {
            "device": info.st_dev, "gid": info.st_gid,
            "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
            "nlink": info.st_nlink, "size": info.st_size,
            "uid": info.st_uid,
        }
        if with_hash:
            digest, size = sha256_fd(fd)
            if size != info.st_size:
                raise TransitionError("bound file size changed while hashing")
            value["sha256"] = digest
        post = os.stat(str(path), follow_symlinks=False)
        if (post.st_dev, post.st_ino, post.st_size,
                stat.S_IMODE(post.st_mode), post.st_uid, post.st_gid,
                post.st_nlink) != (
                info.st_dev, info.st_ino, info.st_size,
                stat.S_IMODE(info.st_mode), info.st_uid, info.st_gid,
                info.st_nlink):
            raise TransitionError("bound path changed while reading: %s" % path)
        return value
    finally:
        os.close(fd)


def capture_tool_records() -> Dict[str, Dict[str, object]]:
    """Bind every external executable authorized by this transition."""
    records: Dict[str, Dict[str, object]] = {}
    for name, path_text, expected_sha256 in TOOL_CONTRACT:
        if name in records or SHA256_RE.fullmatch(expected_sha256) is None:
            raise TransitionError("tool contract is malformed")
        path = Path(path_text)
        if not path.is_absolute():
            raise TransitionError("tool path is not absolute: " + name)
        binding = file_binding(path, with_hash=True)
        if binding["sha256"] != expected_sha256 or binding["uid"] != 0 or \
                binding["gid"] != 0 or binding["mode"] & 0o022 or \
                not binding["mode"] & 0o111 or binding["nlink"] < 1:
            raise TransitionError("hardcoded tool identity changed: " + name)
        records[name] = {
            "binding": binding,
            "path": path_text,
            "sha256": expected_sha256,
        }
    if set(records) != {
            "env", "fuser", "mpstat", "python3", "sudo", "taskset",
            "timeout"}:
        raise TransitionError("tool contract is incomplete")
    return records


def verify_tool_records(value: object) -> Dict[str, Dict[str, object]]:
    """Rehash and replay the exact preparation-time executable bindings."""
    if not isinstance(value, dict):
        raise TransitionError("sealed tool ledger is malformed")
    current = capture_tool_records()
    if set(value) != set(current):
        raise TransitionError("sealed tool ledger changed")
    for name, record in current.items():
        if value.get(name) != record:
            raise TransitionError("sealed tool identity changed: " + name)
    return current


def verify_running_interpreter(
        tool_records: Mapping[str, Mapping[str, object]], *,
        require_exact_path: bool = False) -> Dict[str, object]:
    """Prove the executing controller is the sealed Python inode."""
    try:
        expected = tool_records["python3"]["binding"]
        if not isinstance(expected, dict):
            raise TypeError("Python binding is not a mapping")
        info = os.stat("/proc/self/exe", follow_symlinks=True)
        observed = {"device": info.st_dev, "inode": info.st_ino}
        expected_path = str(tool_records["python3"]["path"])
        resolved_path = str(Path(sys.executable).resolve(strict=True))
        if observed != {
                "device": expected["device"], "inode": expected["inode"]} or \
                resolved_path != expected_path or \
                (require_exact_path and sys.executable != expected_path):
            raise TransitionError(
                "controller interpreter differs from the sealed Python tool")
        return {**observed, "argv_path": sys.executable,
                "resolved_path": resolved_path}
    except TransitionError:
        raise
    except (KeyError, OSError, RuntimeError, TypeError) as exc:
        raise TransitionError(
            "cannot bind the controller interpreter to the sealed Python tool") \
            from exc


def _directory_binding_fd(fd: int) -> Dict[str, int]:
    info = os.fstat(fd)
    if not stat.S_ISDIR(info.st_mode):
        raise TransitionError("bound parent is not a directory")
    return {
        "device": info.st_dev, "gid": info.st_gid,
        "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink, "uid": info.st_uid,
    }


def directory_binding(path: Path) -> Dict[str, int]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | os.O_DIRECTORY | \
        getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(path), flags)
    try:
        observed = _directory_binding_fd(fd)
        post = os.stat(str(path), follow_symlinks=False)
        if observed != {
                "device": post.st_dev, "gid": post.st_gid,
                "inode": post.st_ino, "mode": stat.S_IMODE(post.st_mode),
                "nlink": post.st_nlink, "uid": post.st_uid}:
            raise TransitionError("bound parent path changed")
        return observed
    finally:
        os.close(fd)


def write_new(path: Path, raw: bytes, *, owner_uid: int,
              mode: int = 0o444) -> Dict[str, object]:
    """Create, fsync, seal, and bind one controller-owned record."""
    parent_binding = directory_binding(path.parent)
    if parent_binding["uid"] != owner_uid or parent_binding["mode"] & 0o022:
        raise TransitionError("record parent is outside the controller trust boundary")
    dir_fd = os.open(
        str(path.parent), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    fd = -1
    try:
        if _directory_binding_fd(dir_fd) != parent_binding:
            raise TransitionError("record parent changed before creation")
        fd = os.open(
            path.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0), mode, dir_fd=dir_fd)
        view = memoryview(raw)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                raise OSError(errno.EIO, "short record write")
            view = view[written:]
        os.fsync(fd)
        os.fchmod(fd, mode)
        os.fsync(fd)
        info = os.fstat(fd)
        if info.st_uid != owner_uid or info.st_nlink != 1 or \
                stat.S_IMODE(info.st_mode) != mode:
            raise TransitionError("new record binding is unsafe")
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent_binding:
            raise TransitionError("record parent changed during creation")
    finally:
        if fd >= 0:
            os.close(fd)
        os.close(dir_fd)
    result = file_binding(path, with_hash=True)
    if directory_binding(path.parent) != parent_binding:
        raise TransitionError("record parent path changed after creation")
    if result["sha256"] != sha256_bytes(raw) or result["size"] != len(raw):
        raise TransitionError("new record content changed")
    return result


def copy_sealed(source: Path, destination: Path, expected_sha256: str,
                *, owner_uid: int) -> Dict[str, object]:
    if SHA256_RE.fullmatch(expected_sha256) is None:
        raise TransitionError("expected copy SHA256 is malformed")
    source_fd = os.open(
        str(source), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK |
        getattr(os, "O_NOFOLLOW", 0))
    destination_fd = -1
    dir_fd = -1
    try:
        source_info = os.fstat(source_fd)
        if not stat.S_ISREG(source_info.st_mode) or source_info.st_nlink != 1:
            raise TransitionError("copy source is not a single-link regular file")
        digest, _size = sha256_fd(source_fd)
        if digest != expected_sha256:
            raise TransitionError("copy source SHA256 changed")
        parent = directory_binding(destination.parent)
        if parent["uid"] != owner_uid or parent["mode"] & 0o022:
            raise TransitionError("copy destination parent is unsafe")
        dir_fd = os.open(
            str(destination.parent), os.O_RDONLY | os.O_CLOEXEC |
            os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0))
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("copy destination parent changed")
        destination_fd = os.open(
            destination.name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0), 0o444, dir_fd=dir_fd)
        os.lseek(source_fd, 0, os.SEEK_SET)
        while True:
            block = os.read(source_fd, 1024 * 1024)
            if not block:
                break
            view = memoryview(block)
            while view:
                written = os.write(destination_fd, view)
                if written <= 0:
                    raise OSError(errno.EIO, "short sealed-copy write")
                view = view[written:]
        os.fsync(destination_fd)
        os.fchmod(destination_fd, 0o444)
        os.fsync(destination_fd)
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("copy destination parent changed during copy")
    finally:
        if destination_fd >= 0:
            os.close(destination_fd)
        if dir_fd >= 0:
            os.close(dir_fd)
        os.close(source_fd)
    binding = file_binding(destination, with_hash=True)
    if directory_binding(destination.parent) != parent:
        raise TransitionError("copy destination parent path changed")
    if binding["sha256"] != expected_sha256 or binding["uid"] != owner_uid or \
            binding["mode"] != 0o444 or binding["nlink"] != 1:
        raise TransitionError("sealed copy binding changed")
    return binding


def rename_noreplace(source: Path, destination: Path,
                     expected_binding: Mapping[str, object], *,
                     parent_uid: int) -> Dict[str, object]:
    """Rename one exact inode with renameat2(RENAME_NOREPLACE)."""
    if source.parent != destination.parent or source.name in ("", ".", "..") or \
            destination.name in ("", ".", ".."):
        raise TransitionError("archive rename must stay within one parent")
    parent = directory_binding(source.parent)
    if parent["uid"] != parent_uid or parent["mode"] & 0o022:
        raise TransitionError("archive parent is outside its trust boundary")
    observed = file_binding(source, with_hash="sha256" in expected_binding)
    if observed != dict(expected_binding):
        raise TransitionError("archive source binding changed")
    dir_fd = os.open(
        str(source.parent), os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY |
        getattr(os, "O_NOFOLLOW", 0))
    try:
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("archive parent changed before rename")
        try:
            os.stat(destination.name, dir_fd=dir_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise FileExistsError(errno.EEXIST, "archive destination exists",
                                  str(destination))
        renameat2 = renameat2_function()
        result = renameat2(
            dir_fd, os.fsencode(source.name), dir_fd,
            os.fsencode(destination.name), RENAME_NOREPLACE)
        if result != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), str(destination))
        os.fsync(dir_fd)
        if _directory_binding_fd(dir_fd) != parent:
            raise TransitionError("archive parent changed during rename")
        try:
            os.stat(source.name, dir_fd=dir_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise TransitionError("archive source name survived rename")
    finally:
        os.close(dir_fd)
    archived = file_binding(destination, with_hash="sha256" in expected_binding)
    if archived != dict(expected_binding):
        raise TransitionError("archive destination binding changed")
    if directory_binding(source.parent) != parent:
        raise TransitionError("archive parent binding changed")
    return archived


@dataclass(frozen=True)
class TransitionPlan:
    transition_id: str = TRANSITION_ID
    root: Path = DRY_ROOT
    controller_uid: int = 1000
    deadline_s: float = 540.0
    recovery_reserve_s: float = 60.0
    emergency_recovery_s: float = 60.0
    # Supplied by the external exact-byte reviewer.  Preparation and execution
    # reject an empty/self-derived value; the controller cannot bless its own
    # arbitrary bytes.
    controller_sha256: str = ""
    candidate_sha256: str = CANDIDATE_SHA256
    p32_sha256: str = P32_SHA256
    old_source_sha256: str = OLD_SOURCE_SHA256
    old_boot_id: str = OLD_BOOT_ID
    old_pid: int = 3320493
    old_start_tick: int = 160912119
    old_cmdline_sha256: str = OLD_CMDLINE_SHA256
    old_launcher_pid: int = 3320490
    old_launcher_start_tick: int = 160912119
    old_launcher_cmdline_sha256: str = OLD_LAUNCHER_CMDLINE_SHA256
    old_launcher_affinity: str = "16-59,81-123"
    old_launcher_uids: Tuple[int, int, int, int] = (1000, 0, 0, 0)
    old_process_group: int = 3320490
    old_session: int = 3320465
    old_csv_device: int = 28
    old_csv_inode: int = 10141
    old_pid_device: int = 28
    old_pid_inode: int = 10140
    old_cpu: int = 127
    candidate_cpu: int = 124
    candidate_sibling: int = 60
    old_source: Path = OLD_SOURCE
    old_csv: Path = OLD_CSV
    old_pid_file: Path = OLD_PID_FILE
    old_archive: Path = OLD_ARCHIVE
    old_unclean_archive: Path = OLD_UNCLEAN_ARCHIVE
    old_stale_pid_archive: Path = OLD_STALE_PID_ARCHIVE
    audit_binding: Path = AUDIT_BINDING
    old_argv: Tuple[str, ...] = OLD_ARGV
    old_launch_argv: Tuple[str, ...] = OLD_LAUNCH_ARGV
    replacement_old_argv: Tuple[str, ...] = OLD_REPLACEMENT_ARGV
    replacement_old_cmdline_sha256: str = OLD_REPLACEMENT_CMDLINE_SHA256

    @property
    def receipts(self) -> Path:
        return self.root / "receipts"

    @property
    def frozen(self) -> Path:
        return self.root / "frozen"

    @property
    def sampler(self) -> Path:
        return self.frozen / "wirehair_expo_thermal_sampler.py"

    @property
    def p32(self) -> Path:
        return self.frozen / "wh2_p32_dispatch_timing.py"

    @property
    def controller(self) -> Path:
        return self.frozen / "wh2_thermal_sampler_transition.py"

    @property
    def plan_receipt(self) -> Path:
        return self.root / "transition-plan.json"


def execute_environment(plan: TransitionPlan) -> Dict[str, str]:
    """The complete environment admitted for the live controller process."""
    return {
        "HOME": str(plan.root / "runtime-home"),
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "TZ": "UTC",
    }


def expected_execute_orig_argv(plan: TransitionPlan) -> Tuple[str, ...]:
    if SHA256_RE.fullmatch(plan.controller_sha256) is None:
        raise TransitionError(
            "externally reviewed controller SHA256 is required")
    return (
        "/usr/bin/python3.12", "-I", "-S", "-B", str(plan.controller),
        "--execute-sealed-transition", plan.transition_id,
        "--expected-controller-sha256", plan.controller_sha256,
        "--confirmation", EXECUTE_CONFIRMATION,
    )


def expected_execute_command(plan: TransitionPlan) -> Tuple[str, ...]:
    environment = execute_environment(plan)
    return (
        "/usr/bin/env", "-i",
        *("%s=%s" % (name, environment[name])
          for name in ("HOME", "PATH", "LANG", "LC_ALL", "TZ")),
        *expected_execute_orig_argv(plan),
    )


def python_isolation_contract(plan: TransitionPlan) -> Dict[str, object]:
    prefix = ["/usr/bin/python3.12", "-I", "-S", "-B"]
    return {
        "candidate_argv_prefix": prefix,
        "controller_orig_argv": list(expected_execute_orig_argv(plan)),
        "helper_argv_prefix": prefix,
        "replacement_old_argv": list(plan.replacement_old_argv),
    }


def verify_execute_runtime(
        plan: TransitionPlan, *,
        observed_environment: Optional[Mapping[str, str]] = None,
        observed_orig_argv: Optional[Sequence[str]] = None,
        observed_flags: Optional[Mapping[str, int | bool]] = None,
) -> Dict[str, object]:
    """Replay env -i, Python isolation flags, and the exact execute argv."""
    environment = dict(os.environ if observed_environment is None
                       else observed_environment)
    orig_argv = tuple(sys.orig_argv if observed_orig_argv is None
                      else observed_orig_argv)
    if observed_flags is None:
        flags = {
            name: getattr(sys.flags, name) for name in EXECUTE_FLAG_CONTRACT}
    else:
        flags = dict(observed_flags)
    expected_environment = execute_environment(plan)
    expected_orig_argv = expected_execute_orig_argv(plan)
    if environment != expected_environment:
        raise TransitionError("execute environment differs from exact env -i contract")
    if orig_argv != expected_orig_argv:
        raise TransitionError("execute sys.orig_argv differs from frozen contract")
    if flags != EXECUTE_FLAG_CONTRACT:
        raise TransitionError("execute sys.flags differ from isolated Python contract")
    return {
        "command": list(expected_execute_command(plan)),
        "environment": expected_environment,
        "flags": dict(EXECUTE_FLAG_CONTRACT),
        "sys_orig_argv": list(expected_orig_argv),
    }


class Deadline:
    def __init__(self, duration_s: float, recovery_reserve_s: float,
                 *, clock: Callable[[], float] = time.monotonic) -> None:
        if not math.isfinite(duration_s) or not math.isfinite(recovery_reserve_s) or \
                duration_s <= recovery_reserve_s or recovery_reserve_s <= 0:
            raise ValueError("transition deadline is malformed")
        self._clock = clock
        self.started = clock()
        self.absolute = self.started + duration_s
        self.normal = self.absolute - recovery_reserve_s

    def remaining(self) -> float:
        return max(0.0, self.absolute - self._clock())

    def now(self) -> float:
        return self._clock()

    def exhausted(self) -> bool:
        return self._clock() >= self.absolute

    def require_normal(self, phase: str) -> None:
        if self._clock() >= self.normal:
            raise TransitionDeadline(
                "normal deadline reached before phase " + phase)

    def start_emergency_recovery(
            self, duration_s: float) -> "EmergencyRecoveryBudget":
        return EmergencyRecoveryBudget(duration_s, clock=self._clock)


class EmergencyRecoveryBudget:
    """Fresh post-stop budget; it bounds waits but never vetoes safe actions."""

    def __init__(self, duration_s: float, *,
                 clock: Callable[[], float] = time.monotonic) -> None:
        if not math.isfinite(duration_s) or duration_s <= 0:
            raise ValueError("emergency recovery budget is malformed")
        self._clock = clock
        self.started = clock()
        self.absolute = self.started + duration_s
        self.safety_extension_count = 0
        self.safety_extension_deadline_max: Optional[float] = None

    def now(self) -> float:
        return self._clock()

    def remaining(self) -> float:
        return max(0.0, self.absolute - self._clock())

    def exhausted(self) -> bool:
        return self._clock() >= self.absolute

    def wait_deadline(self, maximum_wait_s: float, *,
                      minimum_safety_wait_s: float = 0.05) -> float:
        if not math.isfinite(maximum_wait_s) or maximum_wait_s <= 0 or \
                not math.isfinite(minimum_safety_wait_s) or \
                minimum_safety_wait_s <= 0 or \
                minimum_safety_wait_s > maximum_wait_s:
            raise ValueError("recovery wait bound is malformed")
        now = self._clock()
        deadline = max(
            min(self.absolute, now + maximum_wait_s),
            now + minimum_safety_wait_s)
        if deadline > self.absolute:
            self.safety_extension_count += 1
            self.safety_extension_deadline_max = max(
                deadline, self.safety_extension_deadline_max or deadline)
        return deadline

    def receipt(self) -> Dict[str, object]:
        observed = self._clock()
        return {
            "absolute_monotonic_s": self.absolute,
            "exhausted": observed >= self.absolute,
            "observed_monotonic_s": observed,
            "remaining_s": max(0.0, self.absolute - observed),
            "safety_extension_count": self.safety_extension_count,
            "safety_extension_deadline_max":
                self.safety_extension_deadline_max,
            "started_monotonic_s": self.started,
        }


class ReceiptJournal:
    def __init__(self, directory: Path, transition_id: str, owner_uid: int,
                 deadline: Deadline) -> None:
        self.directory = directory
        self.transition_id = transition_id
        self.owner_uid = owner_uid
        self.deadline = deadline
        self.sequence = 0
        self.previous_sha256: Optional[str] = None

    def record(self, phase: str, status_value: str,
               payload: Mapping[str, object]) -> Dict[str, object]:
        if not re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", phase) or \
                status_value not in {"started", "completed", "failed"}:
            raise TransitionError("phase receipt identity is malformed")
        sequence = self.sequence
        record = sealed_record(PHASE_SCHEMA, {
            "absolute_deadline_monotonic_s": self.deadline.absolute,
            "created_utc": utc_now(),
            "payload": dict(payload),
            "phase": phase,
            "previous_receipt_sha256": self.previous_sha256,
            "remaining_s": self.deadline.remaining(),
            "sequence": sequence,
            "status": status_value,
            "transition_id": self.transition_id,
        })
        raw = canonical_json(record)
        path = self.directory / (
            "%04d-%s-%s.json" % (sequence, phase, status_value))
        binding = write_new(path, raw, owner_uid=self.owner_uid)
        self.sequence += 1
        self.previous_sha256 = str(binding["sha256"])
        return record

    def replay(self) -> Dict[str, object]:
        """Re-read the complete durable prefix and verify its hash chain."""
        paths = sorted(self.directory.glob("*.json"))
        if len(paths) != self.sequence:
            raise TransitionError("phase receipt roster is incomplete or enlarged")
        previous: Optional[str] = None
        roster = []
        for sequence, path in enumerate(paths):
            match = re.fullmatch(
                r"([0-9]{4})-([a-z][a-z0-9_-]{0,63})-"
                r"(started|completed|failed)\.json", path.name)
            if match is None or int(match.group(1)) != sequence:
                raise TransitionError("phase receipt filename is noncanonical")
            # Reject FIFOs/devices/symlinks without first opening them through
            # the path-following stable-file reader.  The private parent is the
            # namespace trust boundary; this ordering also makes a malformed
            # receipt roster fail promptly instead of blocking recovery.
            binding = file_binding(path, with_hash=True)
            value = load_canonical(path, PHASE_SCHEMA, "phase receipt")
            if binding["uid"] != self.owner_uid or \
                    binding["mode"] != 0o444 or binding["nlink"] != 1 or \
                    value.get("sequence") != sequence or \
                    value.get("phase") != match.group(2) or \
                    value.get("status") != match.group(3) or \
                    value.get("transition_id") != self.transition_id or \
                    value.get("absolute_deadline_monotonic_s") != \
                        self.deadline.absolute or \
                    value.get("previous_receipt_sha256") != previous:
                raise TransitionError("phase receipt chain changed")
            roster.append({"binding": binding, "path": str(path),
                           "phase": match.group(2),
                           "sequence": sequence, "status": match.group(3)})
            previous = str(binding["sha256"])
        if previous != self.previous_sha256:
            raise TransitionError("phase receipt chain head changed")
        return {"count": len(roster), "head_sha256": previous,
                "roster": roster}


class DeferredSignals:
    def __init__(self) -> None:
        self.previous: Dict[int, object] = {}
        self.requested: list[int] = []

    def __enter__(self) -> "DeferredSignals":
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            try:
                self.previous[int(signum)] = signal.getsignal(signum)
                signal.signal(signum, self._record)
            except (OSError, ValueError) as exc:
                for installed, previous in self.previous.items():
                    signal.signal(installed, previous)
                raise TransitionError(
                    "cannot install transition signal guard") from exc
        return self

    def __exit__(self, _kind: object, _value: object, _trace: object) -> None:
        for signum, previous in self.previous.items():
            signal.signal(signum, previous)
        self.previous.clear()
        if _kind is None:
            self.raise_if_requested()

    def _record(self, signum: int, _frame: object) -> None:
        if not self.requested:
            self.requested.append(int(signum))

    def raise_if_requested(self) -> None:
        if self.requested:
            raise TransitionError(
                "deferred controller signal: " +
                signal.Signals(self.requested[0]).name)


def error_record(exc: BaseException) -> Dict[str, str]:
    return {"message": str(exc), "type": type(exc).__name__}


class TransitionController:
    """State machine independent of the live backend for exhaustive tests."""

    def __init__(self, plan: TransitionPlan, backend: object,
                 journal: ReceiptJournal, deadline: Deadline,
                 signal_guard: Optional[DeferredSignals] = None) -> None:
        self.plan = plan
        self.backend = backend
        self.journal = journal
        self.deadline = deadline
        self.signal_guard = signal_guard
        self.recovery_armed = False
        self.dry_accepted = False
        self.primary_error: Optional[BaseException] = None
        self.recovery_errors: list[BaseException] = []
        self.old_stop: Optional[Dict[str, object]] = None
        self.archive: Optional[Dict[str, object]] = None
        self.candidate_accept: Optional[Dict[str, object]] = None
        self.restored: Optional[Dict[str, object]] = None
        self.audit: Optional[Dict[str, object]] = None
        self.final_replay: Optional[Dict[str, object]] = None
        self.receipt_chain: Optional[Dict[str, object]] = None
        self.recovery_budget: Optional[EmergencyRecoveryBudget] = None
        self.emergency_recovery_completion: Optional[Dict[str, object]] = None
        self.scoring_deadline_error: Optional[TransitionDeadline] = None
        self.scoring_checkpoint_completed_monotonic_s: Optional[float] = None
        self.scoring_evidence_completed_monotonic_s: Optional[float] = None

    def _check_signal(self) -> None:
        if self.signal_guard is not None:
            self.signal_guard.raise_if_requested()

    def _phase(self, name: str, action: Callable[[], Mapping[str, object]]) \
            -> Dict[str, object]:
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " action")
            result = dict(action())
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _best_effort_recovery_receipt(
            self, name: str, status_value: str,
            payload: Mapping[str, object]) -> None:
        try:
            self.journal.record(name, status_value, payload)
        except BaseException as exc:
            self.recovery_errors.append(exc)

    def _record_scoring_deadline_exhaustion(
            self, observed_monotonic_s: Optional[float] = None) -> None:
        observed = self.deadline.now() if observed_monotonic_s is None \
            else observed_monotonic_s
        if not math.isfinite(observed) or observed < self.deadline.absolute or \
                self.scoring_deadline_error is not None:
            return
        error = TransitionDeadline(
            "scoring absolute deadline exhausted; safe recovery continued")
        self.scoring_deadline_error = error
        self.recovery_errors.append(error)
        self._best_effort_recovery_receipt(
            "scoring_deadline", "failed", {
                "error": error_record(error),
                "observed_monotonic_s": observed,
                "recovery_actions_continue": True,
            })

    def _stop_old_phase(self) -> Dict[str, object]:
        """Arm in memory only after the durable stop-start receipt exists."""
        name = "old_stop"
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " signal")
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        # From this assignment onward, every exception and deferred signal must
        # enter recovery.  Before it, the live sampler is provably untouched.
        self.recovery_armed = True
        try:
            result = dict(self.backend.stop_old())
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _candidate_accept_phase(self) -> Dict[str, object]:
        """Retain acceptance evidence even if its completion receipt fails."""
        name = "candidate_accept"
        self.deadline.require_normal(name)
        self._check_signal()
        self.journal.record(name, "started", {})
        self._check_signal()
        try:
            self.deadline.require_normal(name + " action")
            result = dict(self.backend.accept_candidate())
            self.candidate_accept = result
        except BaseException as exc:
            try:
                self.journal.record(name, "failed", {"error": error_record(exc)})
            except BaseException:
                pass
            raise
        self.journal.record(name, "completed", result)
        self.deadline.require_normal(name + " completion")
        self._check_signal()
        return result

    def _recover(self) -> None:
        budget = self.deadline.start_emergency_recovery(
            self.plan.emergency_recovery_s)
        self.recovery_budget = budget
        emergency_setup_error: Optional[BaseException] = None
        emergency_setup: Dict[str, object] = {}
        try:
            emergency_setup = dict(
                self.backend.begin_emergency_recovery(budget))
        except BaseException as exc:
            emergency_setup_error = exc
            self.recovery_errors.append(exc)
        self._best_effort_recovery_receipt(
            "emergency_recovery", "started", {
                **budget.receipt(), "backend": emergency_setup})
        self._record_scoring_deadline_exhaustion()

        self._best_effort_recovery_receipt(
            "candidate_cleanup", "started", {})
        try:
            cleanup = dict(self.backend.cleanup_candidate())
        except BaseException as exc:
            self.recovery_errors.append(exc)
            self._best_effort_recovery_receipt(
                "candidate_cleanup", "failed", {"error": error_record(exc)})
        else:
            self._best_effort_recovery_receipt(
                "candidate_cleanup", "completed", cleanup)

        self._best_effort_recovery_receipt("old_restore", "started", {})
        try:
            self.restored = dict(self.backend.restore_old(self.archive))
        except BaseException as exc:
            self.recovery_errors.append(exc)
            self._best_effort_recovery_receipt(
                "old_restore", "failed", {"error": error_record(exc)})
        else:
            if self.archive is None and isinstance(
                    self.restored.get("archive_record"), dict):
                self.archive = dict(self.restored["archive_record"])
            self._best_effort_recovery_receipt(
                "old_restore", "completed", self.restored)

            self._best_effort_recovery_receipt(
                "audit_binding", "started", {})
            try:
                pre_audit_receipt_chain = self.journal.replay()
                self.audit = dict(self.backend.publish_audit_binding(
                    self.archive, self.restored, self.dry_accepted,
                    self.candidate_accept, pre_audit_receipt_chain))
            except BaseException as exc:
                self.recovery_errors.append(exc)
                self._best_effort_recovery_receipt(
                    "audit_binding", "failed", {"error": error_record(exc)})
            else:
                self._best_effort_recovery_receipt(
                    "audit_binding", "completed", self.audit)

                self._best_effort_recovery_receipt(
                    "final_replay", "started", {})
                try:
                    self.final_replay = dict(self.backend.final_replay(
                        self.candidate_accept, self.archive, self.restored,
                        self.audit))
                except BaseException as exc:
                    self.recovery_errors.append(exc)
                    self._best_effort_recovery_receipt(
                        "final_replay", "failed", {"error": error_record(exc)})
                else:
                    self._best_effort_recovery_receipt(
                        "final_replay", "completed", self.final_replay)

        self._record_scoring_deadline_exhaustion()
        completion = budget.receipt()
        self.emergency_recovery_completion = completion
        emergency_error = emergency_setup_error
        if completion["exhausted"]:
            emergency_error = TransitionDeadline(
                "emergency recovery budget exhausted after safe actions were attempted")
            self.recovery_errors.append(emergency_error)
        self._best_effort_recovery_receipt(
            "emergency_recovery",
            "failed" if emergency_error is not None else "completed",
            {**completion,
             **({"error": error_record(emergency_error)}
                if emergency_error is not None else {})})

    def run(self) -> Dict[str, object]:
        try:
            self._phase("hard_preflight", self.backend.hard_preflight)
            self._phase("recovery_arm", self.backend.arm_recovery)
            self.old_stop = self._stop_old_phase()
            if self.old_stop.get("forced") is True:
                raise TransitionError(
                    "forced old stop vetoes candidate dry run")
            self.archive = self._phase("old_archive", self.backend.archive_old)
            if self.archive.get("forced_stop") is True:
                raise TransitionError(
                    "forced old stop vetoes candidate dry run")
            self._phase("candidate_start", self.backend.start_candidate)
            self._phase("candidate_exercise", self.backend.exercise_candidate)
            self._phase("candidate_stop", self.backend.stop_candidate)
            self._candidate_accept_phase()
            self.dry_accepted = True
        except BaseException as exc:
            self.primary_error = exc
        finally:
            if self.recovery_armed:
                try:
                    self._recover()
                except BaseException as exc:
                    self.recovery_errors.append(exc)

        # Signals received during recovery are deliberately held until cleanup,
        # old restoration, and audit publication have all had their chance.  Do
        # not silently turn such a request into a successful transition.
        if self.signal_guard is not None:
            try:
                self.signal_guard.raise_if_requested()
            except BaseException as exc:
                if self.primary_error is None:
                    self.primary_error = exc
                else:
                    self.recovery_errors.append(exc)

        # Include a scoring overrun that happened during the final recovery
        # receipt/signal replay.  This can only change acceptance; recovery has
        # already had every safe action attempted.
        if self.recovery_armed:
            self._record_scoring_deadline_exhaustion()

        # This is the only terminal journal write.  It deliberately contains
        # the complete evidence candidate but no outcome claim.  Only after the
        # durable write returns can the absolute scoring boundary be evaluated
        # without leaving a serialized success that a late write contradicts.
        try:
            self.journal.record("terminal", "started", {
                "archived_pre_dry": self.archive,
                "candidate_accept": self.candidate_accept,
                "dry_accepted": self.dry_accepted,
                "future_audit_binding": self.audit,
                "final_replay": self.final_replay,
                "old_sampler_restored": self.restored,
                "primary_error": error_record(self.primary_error)
                    if self.primary_error is not None else None,
                "outcome_pending": True,
                "recovery_errors_before_checkpoint": [
                    error_record(exc) for exc in self.recovery_errors],
                "recovery_attempts_finished": self.recovery_armed,
                "transition_id": self.plan.transition_id,
                "transition_success": None,
            })
            self.scoring_checkpoint_completed_monotonic_s = self.deadline.now()
        except BaseException as exc:
            self.recovery_errors.append(exc)
        self._record_scoring_deadline_exhaustion(
            self.scoring_checkpoint_completed_monotonic_s)
        try:
            self.receipt_chain = self.journal.replay()
        except BaseException as exc:
            self.recovery_errors.append(exc)
        # Capture the final scoring observation exactly once.  A deliberately
        # slow replay is scored, while later terminal-object serialization can
        # neither contradict this observation nor reopen the recovery path.
        self.scoring_evidence_completed_monotonic_s = self.deadline.now()
        if self.scoring_evidence_completed_monotonic_s >= \
                self.deadline.absolute and \
                self.scoring_deadline_error is None:
            self._record_scoring_deadline_exhaustion(
                self.scoring_evidence_completed_monotonic_s)
            try:
                self.receipt_chain = self.journal.replay()
            except BaseException as exc:
                self.recovery_errors.append(exc)

        success = self.dry_accepted and self.restored is not None and \
            self.audit is not None and self.final_replay is not None and \
            self.receipt_chain is not None and \
            self.primary_error is None and not self.recovery_errors
        terminal = sealed_record(TRANSITION_SCHEMA, {
            "archived_pre_dry": self.archive,
            "candidate_accept": self.candidate_accept,
            "dry_accepted": self.dry_accepted,
            "finished_utc": utc_now(),
            "future_audit_binding": self.audit,
            "final_replay": self.final_replay,
            "old_sampler_restored": self.restored,
            "primary_error": error_record(self.primary_error)
                if self.primary_error is not None else None,
            "recovery_armed": self.recovery_armed,
            "emergency_recovery": self.emergency_recovery_completion,
            "recovery_errors": [error_record(exc)
                                for exc in self.recovery_errors],
            "receipt_chain_prefix": self.receipt_chain,
            "scoring_deadline": {
                "absolute_monotonic_s": self.deadline.absolute,
                "checkpoint_completed_monotonic_s":
                    self.scoring_checkpoint_completed_monotonic_s,
                "evidence_completed_monotonic_s":
                    self.scoring_evidence_completed_monotonic_s,
                "evidence_replay_defines_scoring_boundary": True,
                "exhausted": self.scoring_evidence_completed_monotonic_s >=
                    self.deadline.absolute,
                "remaining_s": max(
                    0.0, self.deadline.absolute -
                    self.scoring_evidence_completed_monotonic_s),
                "terminal_checkpoint_is_outcome_neutral": True,
            },
            "success": success,
            "transition_id": self.plan.transition_id,
        })
        if not success:
            messages = []
            if self.primary_error is not None:
                messages.append("primary=" + str(self.primary_error))
            messages.extend("recovery=" + str(exc)
                            for exc in self.recovery_errors)
            raise TransitionError(
                "thermal transition failed: " + "; ".join(messages))
        return terminal


PIDFD_SIGNAL_PROGRAM = r"""
import ctypes
import hashlib
import os
import re
import signal
import sys

pid = int(sys.argv[1])
expected_boot = sys.argv[2]
expected_tick = int(sys.argv[3])
expected_cmdline_sha = sys.argv[4]
expected_ppid = int(sys.argv[5])
expected_pgrp = int(sys.argv[6])
expected_session = int(sys.argv[7])
expected_cpu = sys.argv[8]
expected_exe_device = int(sys.argv[9])
expected_exe_inode = int(sys.argv[10])
signum = int(sys.argv[11])
if signum not in (signal.SIGTERM, signal.SIGKILL):
    raise SystemExit(64)
libc = ctypes.CDLL(None, use_errno=True)
pidfd = libc.syscall(434, pid, 0)
if pidfd < 0:
    raise OSError(ctypes.get_errno(), 'pidfd_open')
try:
    with open('/proc/sys/kernel/random/boot_id', 'r', encoding='ascii') as stream:
        boot = stream.read().strip()
    stat_raw = open('/proc/%d/stat' % pid, 'rb').read()
    cmdline_raw = open('/proc/%d/cmdline' % pid, 'rb').read()
    status = open('/proc/%d/status' % pid, 'r', encoding='ascii').read()
    executable = os.stat('/proc/%d/exe' % pid, follow_symlinks=True)
    right = stat_raw.rfind(b')')
    fields = stat_raw[right + 2:].split() if right >= 0 else []
    uid = re.search(r'^Uid:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$',
                    status, re.MULTILINE)
    affinity = re.search(r'^Cpus_allowed_list:\s*(\S+)\s*$',
                         status, re.MULTILINE)
    if (boot != expected_boot or len(fields) < 20 or
            int(fields[1]) not in (expected_ppid, 1) or
            int(fields[2]) != expected_pgrp or
            int(fields[3]) != expected_session or
            int(fields[19]) != expected_tick or
            hashlib.sha256(cmdline_raw).hexdigest() != expected_cmdline_sha or
            (executable.st_dev, executable.st_ino) !=
                (expected_exe_device, expected_exe_inode) or
            uid is None or tuple(map(int, uid.groups())) != (0, 0, 0, 0) or
            affinity is None or affinity.group(1) != expected_cpu):
        raise SystemExit(73)
    if libc.syscall(424, pidfd, signum, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), 'pidfd_send_signal')
finally:
    os.close(pidfd)
""".strip()


LAUNCHER_PIDFD_SIGNAL_PROGRAM = r"""
import ctypes
import hashlib
import os
import re
import signal
import sys

pid = int(sys.argv[1])
expected_boot = sys.argv[2]
expected_tick = int(sys.argv[3])
expected_cmdline_sha = sys.argv[4]
expected_pgrp = int(sys.argv[5])
expected_session = int(sys.argv[6])
expected_cpu = sys.argv[7]
expected_uids = tuple(map(int, sys.argv[8].split(',')))
expected_exe_device = int(sys.argv[9])
expected_exe_inode = int(sys.argv[10])
signum = int(sys.argv[11])
if len(expected_uids) != 4 or signum not in (signal.SIGTERM, signal.SIGKILL):
    raise SystemExit(64)
libc = ctypes.CDLL(None, use_errno=True)
pidfd = libc.syscall(434, pid, 0)
if pidfd < 0:
    raise OSError(ctypes.get_errno(), 'pidfd_open')
try:
    with open('/proc/sys/kernel/random/boot_id', 'r', encoding='ascii') as stream:
        boot = stream.read().strip()
    stat_raw = open('/proc/%d/stat' % pid, 'rb').read()
    cmdline_raw = open('/proc/%d/cmdline' % pid, 'rb').read()
    status = open('/proc/%d/status' % pid, 'r', encoding='ascii').read()
    executable = os.stat('/proc/%d/exe' % pid, follow_symlinks=True)
    right = stat_raw.rfind(b')')
    fields = stat_raw[right + 2:].split() if right >= 0 else []
    uid = re.search(r'^Uid:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$',
                    status, re.MULTILINE)
    affinity = re.search(r'^Cpus_allowed_list:\s*(\S+)\s*$',
                         status, re.MULTILINE)
    if (boot != expected_boot or len(fields) < 20 or int(fields[1]) != 1 or
            int(fields[2]) != expected_pgrp or
            int(fields[3]) != expected_session or
            int(fields[19]) != expected_tick or
            hashlib.sha256(cmdline_raw).hexdigest() != expected_cmdline_sha or
            (executable.st_dev, executable.st_ino) !=
                (expected_exe_device, expected_exe_inode) or
            uid is None or tuple(map(int, uid.groups())) != expected_uids or
            affinity is None or affinity.group(1) != expected_cpu):
        raise SystemExit(73)
    if libc.syscall(424, pidfd, signum, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), 'pidfd_send_signal')
finally:
    os.close(pidfd)
""".strip()


def _parse_proc_stat(raw: bytes) -> Dict[str, int | str]:
    right = raw.rfind(b")")
    fields = raw[right + 2:].split() if right >= 0 else []
    if len(fields) < 37:
        raise TransitionError("process stat is truncated")
    try:
        return {
            "state": fields[0].decode("ascii"),
            "ppid": int(fields[1]), "process_group": int(fields[2]),
            "session": int(fields[3]), "start_tick": int(fields[19]),
            "processor": int(fields[36]),
        }
    except (UnicodeDecodeError, ValueError) as exc:
        raise TransitionError("process stat is malformed") from exc


def capture_owned_session_leader(
        pid: int, *, proc_root: Path = Path("/proc"),
) -> int:
    if type(pid) is not int or pid <= 1:
        raise TransitionError("owned launcher PID is outside its domain")
    try:
        launcher = _parse_proc_stat(
            (proc_root / str(pid) / "stat").read_bytes())
    except FileNotFoundError:
        # sudo -b may have exited after leaving its exact background session.
        # The session cleanup helper treats zero as valid only while the leader
        # PID remains absent, and refuses a live/reused leader at this identity.
        return 0
    if launcher["process_group"] != pid or launcher["session"] != pid:
        raise TransitionError(
            "old restart launcher did not create an exact session")
    return int(launcher["start_tick"])


def _status_field(status: str, name: str) -> str:
    match = re.search(r"^" + re.escape(name) + r":\s*(\S(?:.*\S)?)\s*$",
                      status, re.MULTILINE)
    if match is None:
        raise TransitionError("process status lacks " + name)
    return match.group(1)


def parse_mpstat_idle_receipt(
        raw: bytes, cpu_pair: Sequence[int], *, expected_intervals: int = 3,
) -> Dict[int, float]:
    """Validate the exact finite per-CPU idle receipt requested from mpstat."""
    pair = set(cpu_pair)
    if len(pair) != 2 or any(type(cpu) is not int or cpu < 0 for cpu in pair) or \
            type(expected_intervals) is not int or expected_intervals <= 0:
        raise TransitionError("mpstat receipt contract is malformed")

    def reject_constant(value: str) -> object:
        raise ValueError("non-finite JSON constant: " + value)

    try:
        document = json.loads(
            raw.decode("ascii"), parse_constant=reject_constant)
        hosts = document["sysstat"]["hosts"]
        if not isinstance(hosts, list) or len(hosts) != 1:
            raise TransitionError("mpstat JSON has an unexpected host count")
        statistics = hosts[0]["statistics"]
        if not isinstance(statistics, list) or \
                len(statistics) != expected_intervals:
            raise TransitionError("mpstat JSON has an unexpected interval count")
        idle_by_cpu = {cpu: [] for cpu in pair}
        for interval in statistics:
            loads = interval["cpu-load"]
            if not isinstance(loads, list):
                raise TransitionError("mpstat cpu-load receipt is malformed")
            observed: Dict[int, float] = {}
            for load in loads:
                if not isinstance(load, dict):
                    raise TransitionError("mpstat CPU receipt is malformed")
                cpu_text = str(load["cpu"])
                if not cpu_text.isdecimal() or int(cpu_text) not in pair:
                    continue
                cpu = int(cpu_text)
                if cpu in observed:
                    raise TransitionError("mpstat CPU receipt is duplicated")
                idle_value = load["idle"]
                if type(idle_value) not in (int, float):
                    raise TransitionError("mpstat idle receipt is not numeric")
                idle = float(idle_value)
                if not math.isfinite(idle) or not 0.0 <= idle <= 100.0:
                    raise TransitionError("mpstat idle receipt is outside its domain")
                observed[cpu] = idle
            if set(observed) != pair:
                raise TransitionError("mpstat interval lacks the exact CPU pair")
            for cpu, idle in observed.items():
                idle_by_cpu[cpu].append(idle)
    except TransitionError:
        raise
    except (KeyError, TypeError, ValueError, UnicodeDecodeError,
            json.JSONDecodeError) as exc:
        raise TransitionError("mpstat JSON is malformed") from exc
    if any(len(values) != expected_intervals or
           sum(values) / expected_intervals < 99.0
           for values in idle_by_cpu.values()):
        raise TransitionError("candidate CPU pair is not at least 99% idle")
    return {
        cpu: sum(values) / expected_intervals
        for cpu, values in idle_by_cpu.items()
    }


def _proc_identity(pid: int) -> Dict[str, object]:
    if pid <= 1:
        raise TransitionError("process PID is outside its domain")
    root = Path("/proc") / str(pid)
    stat_raw = (root / "stat").read_bytes()
    cmdline_raw = (root / "cmdline").read_bytes()
    status = (root / "status").read_text(encoding="ascii")
    executable_info = os.stat(str(root / "exe"), follow_symlinks=True)
    if not cmdline_raw.endswith(b"\0"):
        raise TransitionError("process cmdline is not NUL terminated")
    try:
        argv = tuple(
            item.decode("ascii") for item in cmdline_raw[:-1].split(b"\0"))
    except UnicodeDecodeError as exc:
        raise TransitionError("process cmdline is not ASCII") from exc
    parsed = _parse_proc_stat(stat_raw)
    uid_text = _status_field(status, "Uid")
    try:
        uids = tuple(int(value) for value in uid_text.split())
    except ValueError as exc:
        raise TransitionError("process UID receipt is malformed") from exc
    if len(uids) != 4:
        raise TransitionError("process UID receipt is malformed")
    return {
        **parsed, "affinity": _status_field(status, "Cpus_allowed_list"),
        "argv": list(argv), "cmdline_sha256": sha256_bytes(cmdline_raw),
        "executable": {"device": executable_info.st_dev,
                       "inode": executable_info.st_ino},
        "pid": pid, "uids": list(uids),
    }


def _mkdir_exact(path: Path, mode: int, owner_uid: int) -> None:
    path.mkdir(mode=mode)
    os.chmod(path, mode)
    observed = directory_binding(path)
    if observed["uid"] != owner_uid or observed["mode"] != mode:
        raise TransitionError("new directory binding is unsafe: %s" % path)
    fsync_directory(path.parent)


def prepare_transition(plan: TransitionPlan, *, controller_source: Path,
                       candidate_source: Path = REPO_CANDIDATE,
                       p32_source: Path = REPO_P32) -> Dict[str, object]:
    """Freeze all executable bytes without inspecting live state or hardware."""
    if os.geteuid() != plan.controller_uid:
        raise TransitionError("transition preparation controller UID changed")
    if SHA256_RE.fullmatch(plan.controller_sha256) is None or \
            sha256_file(controller_source) != plan.controller_sha256:
        raise TransitionError(
            "controller bytes differ from the externally reviewed SHA256")
    try:
        Deadline(plan.deadline_s, plan.recovery_reserve_s)
        EmergencyRecoveryBudget(plan.emergency_recovery_s)
    except ValueError as exc:
        raise TransitionError("transition deadline plan is malformed") from exc
    tool_records = capture_tool_records()
    if plan.root.exists() or plan.root.is_symlink():
        raise TransitionError("transition root already exists")
    _mkdir_exact(plan.root, 0o700, plan.controller_uid)
    try:
        for directory in (
                plan.frozen, plan.root / "segments", plan.root / "interrupted",
                plan.root / "runtime-home", plan.receipts):
            _mkdir_exact(directory, 0o700, plan.controller_uid)
        sampler_binding = copy_sealed(
            candidate_source, plan.sampler, plan.candidate_sha256,
            owner_uid=plan.controller_uid)
        p32_binding = copy_sealed(
            p32_source, plan.p32, plan.p32_sha256,
            owner_uid=plan.controller_uid)
        controller_binding = copy_sealed(
            controller_source, plan.controller, plan.controller_sha256,
            owner_uid=plan.controller_uid)
        os.chmod(plan.frozen, 0o500)
        fsync_directory(plan.root)
        plan_value = sealed_record(PLAN_SCHEMA, {
            "absolute_deadline_s": plan.deadline_s,
            "candidate": sampler_binding,
            "candidate_cpu": plan.candidate_cpu,
            "candidate_sha256": plan.candidate_sha256,
            "candidate_sibling": plan.candidate_sibling,
            "controller": controller_binding,
            "controller_sha256": plan.controller_sha256,
            "controller_uid": plan.controller_uid,
            "created_utc": utc_now(),
            "destinations": {
                "audit_binding": str(plan.audit_binding),
                "old_archive": str(plan.old_archive),
                "old_stale_pid_archive": str(plan.old_stale_pid_archive),
                "old_unclean_archive": str(plan.old_unclean_archive),
            },
            "emergency_recovery_s": plan.emergency_recovery_s,
            "old": {
                "argv": list(plan.old_argv),
                "boot_id": plan.old_boot_id,
                "cmdline_sha256": plan.old_cmdline_sha256,
                "cpu": plan.old_cpu, "pid": plan.old_pid,
                "csv": {
                    "device": plan.old_csv_device,
                    "inode": plan.old_csv_inode,
                    "path": str(plan.old_csv),
                },
                "launcher": {
                    "affinity": plan.old_launcher_affinity,
                    "cmdline_sha256": plan.old_launcher_cmdline_sha256,
                    "pid": plan.old_launcher_pid,
                    "start_tick": plan.old_launcher_start_tick,
                    "uids": list(plan.old_launcher_uids),
                },
                "process_group": plan.old_process_group,
                "pid_file": {
                    "device": plan.old_pid_device,
                    "inode": plan.old_pid_inode,
                    "path": str(plan.old_pid_file),
                },
                "session": plan.old_session,
                "source_path": str(plan.old_source),
                "source_sha256": plan.old_source_sha256,
                "start_tick": plan.old_start_tick,
                "replacement_argv": list(plan.replacement_old_argv),
                "replacement_cmdline_sha256":
                    plan.replacement_old_cmdline_sha256,
            },
            "p32": p32_binding,
            "p32_sha256": plan.p32_sha256,
            "python_isolation": python_isolation_contract(plan),
            "recovery_reserve_s": plan.recovery_reserve_s,
            "root": str(plan.root),
            "execute_contract": {
                "command": list(expected_execute_command(plan)),
                "environment": execute_environment(plan),
                "flags": dict(EXECUTE_FLAG_CONTRACT),
                "sys_orig_argv": list(expected_execute_orig_argv(plan)),
            },
            "tools": tool_records,
            "transition_id": plan.transition_id,
        })
        plan_binding = write_new(
            plan.plan_receipt, canonical_json(plan_value),
            owner_uid=plan.controller_uid)
        return {"plan": plan_binding, "value": plan_value}
    except BaseException:
        # Preparation happens before any live inspection or stop.  Preserve a
        # partial root for forensic inspection rather than recursively remove.
        raise


def verify_transition_plan(plan: TransitionPlan) -> Dict[str, object]:
    if os.geteuid() != plan.controller_uid:
        raise TransitionError("controller UID differs from the frozen plan")
    value = load_canonical(plan.plan_receipt, PLAN_SCHEMA, "transition plan")
    if value.get("root") != str(plan.root) or \
            value.get("transition_id") != plan.transition_id or \
            value.get("controller_sha256") != plan.controller_sha256 or \
            value.get("controller_uid") != plan.controller_uid or \
            value.get("candidate_sha256") != plan.candidate_sha256 or \
            value.get("p32_sha256") != plan.p32_sha256 or \
            value.get("candidate_cpu") != plan.candidate_cpu or \
            value.get("candidate_sibling") != plan.candidate_sibling or \
            value.get("absolute_deadline_s") != plan.deadline_s or \
            value.get("recovery_reserve_s") != plan.recovery_reserve_s or \
            value.get("emergency_recovery_s") != plan.emergency_recovery_s:
        raise TransitionError("transition plan contract changed")
    try:
        Deadline(plan.deadline_s, plan.recovery_reserve_s)
        EmergencyRecoveryBudget(plan.emergency_recovery_s)
    except ValueError as exc:
        raise TransitionError("transition deadline plan is malformed") from exc
    expected_old = {
        "argv": list(plan.old_argv),
        "boot_id": plan.old_boot_id,
        "cmdline_sha256": plan.old_cmdline_sha256,
        "cpu": plan.old_cpu, "pid": plan.old_pid,
        "csv": {
            "device": plan.old_csv_device,
            "inode": plan.old_csv_inode,
            "path": str(plan.old_csv),
        },
        "launcher": {
            "affinity": plan.old_launcher_affinity,
            "cmdline_sha256": plan.old_launcher_cmdline_sha256,
            "pid": plan.old_launcher_pid,
            "start_tick": plan.old_launcher_start_tick,
            "uids": list(plan.old_launcher_uids),
        },
        "process_group": plan.old_process_group,
        "pid_file": {
            "device": plan.old_pid_device,
            "inode": plan.old_pid_inode,
            "path": str(plan.old_pid_file),
        },
        "session": plan.old_session,
        "source_path": str(plan.old_source),
        "source_sha256": plan.old_source_sha256,
        "start_tick": plan.old_start_tick,
        "replacement_argv": list(plan.replacement_old_argv),
        "replacement_cmdline_sha256": plan.replacement_old_cmdline_sha256,
    }
    if value.get("old") != expected_old:
        raise TransitionError("frozen old-sampler plan changed")
    if value.get("destinations") != {
            "audit_binding": str(plan.audit_binding),
            "old_archive": str(plan.old_archive),
            "old_stale_pid_archive": str(plan.old_stale_pid_archive),
            "old_unclean_archive": str(plan.old_unclean_archive)}:
        raise TransitionError("transition destination plan changed")
    for path, key, digest in (
            (plan.sampler, "candidate", plan.candidate_sha256),
            (plan.p32, "p32", plan.p32_sha256),
            (plan.controller, "controller", plan.controller_sha256)):
        binding = file_binding(path, with_hash=True)
        if binding != value.get(key) or \
                (digest is not None and binding["sha256"] != digest) or \
                binding["uid"] != plan.controller_uid or \
                binding["mode"] != 0o444 or binding["nlink"] != 1:
            raise TransitionError("frozen transition file changed: " + key)
    frozen_parent = directory_binding(plan.frozen)
    if frozen_parent["uid"] != plan.controller_uid or \
            frozen_parent["mode"] != 0o500:
        raise TransitionError("frozen transition parent mode changed")
    for directory in (
            plan.root, plan.root / "segments", plan.root / "interrupted",
            plan.root / "runtime-home", plan.receipts):
        binding = directory_binding(directory)
        if binding["uid"] != plan.controller_uid or binding["mode"] != 0o700:
            raise TransitionError("transition output parent changed")
    tools = verify_tool_records(value.get("tools"))
    if value.get("execute_contract") != {
            "command": list(expected_execute_command(plan)),
            "environment": execute_environment(plan),
            "flags": dict(EXECUTE_FLAG_CONTRACT),
            "sys_orig_argv": list(expected_execute_orig_argv(plan))}:
        raise TransitionError("frozen execute runtime contract changed")
    if value.get("python_isolation") != python_isolation_contract(plan):
        raise TransitionError("frozen Python isolation contract changed")
    if not plan.old_argv or \
            plan.old_argv[0] != tools["python3"]["path"] or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.old_argv) + b"\0") != \
                plan.old_cmdline_sha256 or \
            tuple(plan.replacement_old_argv[:4]) != (
                tools["python3"]["path"], "-I", "-S", "-B") or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.replacement_old_argv) + b"\0") != \
                plan.replacement_old_cmdline_sha256 or \
            len(plan.old_launch_argv) < 4 or \
            tuple(plan.old_launch_argv[:3]) != ("sudo", "-n", "-b") or \
            plan.old_launch_argv[3] != tools["taskset"]["path"] or \
            sha256_bytes(
                b"\0".join(value.encode("ascii")
                            for value in plan.old_launch_argv) + b"\0") != \
                plan.old_launcher_cmdline_sha256:
        raise TransitionError(
            "old sampler argv differs from the sealed executable ledger")
    return value


def load_verified_p32(path: Path, expected_sha256: str) -> ModuleType:
    binding = file_binding(path, with_hash=True)
    if binding["sha256"] != expected_sha256 or binding["mode"] != 0o444 or \
            binding["nlink"] != 1:
        raise TransitionError("P32 helper binding changed")
    name = "_wh2_transition_frozen_p32"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise TransitionError("cannot load frozen P32 helper")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None)
        raise
    if sha256_file(path) != expected_sha256:
        raise TransitionError("P32 helper changed while loading")
    return module


class LiveBackend:
    """Live operations; instantiated only by the confirmed execute entry point."""

    def __init__(self, plan: TransitionPlan, p32: ModuleType,
                 deadline: Deadline,
                 tool_records: Mapping[str, Mapping[str, object]],
                 execute_runtime: Mapping[str, object]) -> None:
        self.plan = plan
        self.p32 = p32
        self.deadline = deadline
        self.old_preflight: Optional[Dict[str, object]] = None
        self.stop_record: Optional[Dict[str, object]] = None
        self.archive_record: Optional[Dict[str, object]] = None
        self.candidate_owner: Optional[object] = None
        self.candidate_timing_start: Optional[float] = None
        self.candidate_benchmark_end: Optional[float] = None
        self.candidate_terminal: Optional[Dict[str, object]] = None
        self.candidate_cleanup_complete = False
        self.recovery_budget: Optional[EmergencyRecoveryBudget] = None
        self.tools = verify_tool_records(tool_records)
        self.execute_runtime = dict(execute_runtime)
        self.design: Dict[str, object] = {
            "controller_uid": plan.controller_uid,
            "immutable_files": {
                "frozen/wirehair_expo_thermal_sampler.py":
                    plan.candidate_sha256,
            },
            "thermal_core": plan.candidate_cpu,
            "tools": self.tools,
            "python_isolation": python_isolation_contract(plan),
            "transition_execute_runtime": self.execute_runtime,
        }

    def begin_emergency_recovery(
            self, budget: EmergencyRecoveryBudget) -> Mapping[str, object]:
        if self.recovery_budget is not None or \
                not isinstance(budget, EmergencyRecoveryBudget):
            raise TransitionError("emergency recovery budget was installed twice")
        self.recovery_budget = budget
        tools = self._verify_tools()
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed before emergency recovery")
        return {"controller_interpreter": verify_running_interpreter(tools),
                "execute_runtime": runtime,
                "tools": tools}

    def _recovery_wait_deadline(
            self, maximum_wait_s: float, *,
            minimum_safety_wait_s: float = 0.05) -> float:
        if self.recovery_budget is None:
            raise TransitionError("recovery action lacks its emergency budget")
        return self.recovery_budget.wait_deadline(
            maximum_wait_s,
            minimum_safety_wait_s=minimum_safety_wait_s)

    def _recovery_now(self) -> float:
        if self.recovery_budget is None:
            raise TransitionError("recovery action lacks its emergency budget")
        return self.recovery_budget.now()

    def _verify_tools(self) -> Dict[str, Dict[str, object]]:
        return verify_tool_records(self.tools)

    def _environment(self) -> Dict[str, str]:
        return self.p32.sanitized_environment(
            self.plan.root / "runtime-home", allocator=False)

    def _i2c_readers(self) -> Tuple[int, ...]:
        return self.p32.sole_i2c_readers(
            Path(str(self.tools["fuser"]["path"])),
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])))

    def _fuser(self, path: Path) -> Tuple[int, ...]:
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])),
            (str(self.tools["fuser"]["path"]), str(path)),
            self._environment())
        if rc == 1 and not stdout and not stderr:
            return ()
        if rc != 0 or re.fullmatch(rb"(?: +[1-9][0-9]*)+", stdout) is None:
            raise TransitionError("fuser result is not canonical for %s" % path)
        expected_label = re.escape(os.fsencode(str(path))) + rb":[ ]*\n"
        if re.fullmatch(expected_label, stderr) is None:
            raise TransitionError("fuser label is not canonical for %s" % path)
        return tuple(sorted(set(int(value) for value in stdout.split())))

    def _verify_topology_and_occupancy(self) -> Dict[str, object]:
        pair = {self.plan.candidate_cpu, self.plan.candidate_sibling}
        expected = "%d,%d" % tuple(sorted(pair))
        for cpu in sorted(pair):
            value = (Path("/sys/devices/system/cpu") / ("cpu%d" % cpu) /
                     "topology/thread_siblings_list").read_text(
                         encoding="ascii").strip()
            if value != expected:
                raise TransitionError("candidate CPU sibling topology changed")
        online = set(self.p32.parse_cpu_list(
            Path("/sys/devices/system/cpu/online").read_text(
                encoding="ascii").strip()))
        if not pair <= online:
            raise TransitionError("candidate CPU pair is offline")

        def scan() -> Tuple[list[Dict[str, object]], list[Dict[str, object]]]:
            pinned = []
            executing = []
            for task in Path("/proc").glob("[0-9]*/task/[0-9]*"):
                try:
                    pid = int(task.parent.parent.name)
                    tid = int(task.name)
                    cmdline = (Path("/proc") / str(pid) / "cmdline").read_bytes()
                    if not cmdline or pid == os.getpid():
                        continue
                    status = (task / "status").read_text(encoding="ascii")
                    allowed = set(self.p32.parse_cpu_list(
                        _status_field(status, "Cpus_allowed_list")))
                    proc_stat = _parse_proc_stat((task / "stat").read_bytes())
                except (OSError, TransitionError, ValueError) as exc:
                    if task.exists():
                        raise TransitionError(
                            "cannot classify a live task on the candidate CPU pair") \
                            from exc
                    continue
                record = {"pid": pid, "tid": tid,
                          "cmdline_sha256": sha256_bytes(cmdline)}
                if allowed and allowed <= pair:
                    pinned.append(record)
                if proc_stat["state"] == "R" and proc_stat["processor"] in pair:
                    executing.append(record)
            return pinned, executing

        pinned_before, executing_before = scan()
        command = [str(self.tools["mpstat"]["path"]), "-o", "JSON", "-P",
                   "%d,%d" % tuple(sorted(pair)), "1", "3"]
        rc, stdout, stderr = self.p32.run_bounded(
            command, self._environment(), 6.0, 128 * 1024, 4096)
        if rc != 0 or stderr:
            raise TransitionError("mpstat occupancy check failed")
        idle_average = parse_mpstat_idle_receipt(stdout, tuple(sorted(pair)))
        pinned_after, executing_after = scan()
        if pinned_before or executing_before or pinned_after or executing_after:
            raise TransitionError("candidate CPU pair has a user-space occupant")
        return {
            "cpu_pair": sorted(pair),
            "idle_average_pct": {
                str(cpu): average for cpu, average in idle_average.items()},
            "mpstat_sha256": sha256_bytes(stdout),
        }

    def _validate_old_child_identity(self) -> Dict[str, object]:
        current = _proc_identity(self.plan.old_pid)
        exact = {
            "affinity": str(self.plan.old_cpu),
            "argv": list(self.plan.old_argv),
            "cmdline_sha256": self.plan.old_cmdline_sha256,
            "pid": self.plan.old_pid,
            "process_group": self.plan.old_process_group,
            "session": self.plan.old_session,
            "start_tick": self.plan.old_start_tick,
            "uids": [0, 0, 0, 0],
        }
        for key, value in exact.items():
            if current.get(key) != value:
                raise TransitionError("old sampler identity changed: " + key)
        python_binding = self.tools["python3"]["binding"]
        if current.get("executable") != {
                "device": python_binding["device"],
                "inode": python_binding["inode"]}:
            raise TransitionError("old sampler identity changed: executable")
        if current.get("ppid") not in (self.plan.old_launcher_pid, 1):
            raise TransitionError("old sampler identity changed: ppid")
        boot_id = Path("/proc/sys/kernel/random/boot_id").read_text(
            encoding="ascii").strip()
        if boot_id != self.plan.old_boot_id:
            raise TransitionError("boot ID changed from the sealed old identity")
        return {"child": current, "boot_id": boot_id}

    def _validate_old_identity(self) -> Dict[str, object]:
        result = self._validate_old_child_identity()
        if result["child"].get("ppid") != self.plan.old_launcher_pid:
            raise TransitionError("old sampler launcher parent changed")
        launcher = _proc_identity(self.plan.old_launcher_pid)
        launcher_exact = {
            "affinity": self.plan.old_launcher_affinity,
            "argv": list(self.plan.old_launch_argv),
            "cmdline_sha256": self.plan.old_launcher_cmdline_sha256,
            "pid": self.plan.old_launcher_pid,
            "ppid": 1,
            "process_group": self.plan.old_process_group,
            "session": self.plan.old_session,
            "start_tick": self.plan.old_launcher_start_tick,
            "uids": list(self.plan.old_launcher_uids),
        }
        sudo_binding = self.tools["sudo"]["binding"]
        launcher_exact["executable"] = {
            "device": sudo_binding["device"],
            "inode": sudo_binding["inode"],
        }
        for key, value in launcher_exact.items():
            if launcher.get(key) != value:
                raise TransitionError("old launcher identity changed: " + key)
        result["launcher"] = launcher
        return result

    def _require_exclusive_old_readers(self) -> Dict[str, object]:
        i2c_readers = self._i2c_readers()
        if i2c_readers != (self.plan.old_pid,):
            raise TransitionError("old sampler is not the sole I2C reader")
        csv_readers = self._fuser(self.plan.old_csv)
        if csv_readers != (self.plan.old_pid,):
            raise TransitionError(
                "fuzz/GFNI timing/auditor or unknown process is reading the live CSV")
        return {"csv_readers": list(csv_readers),
                "i2c_readers": list(i2c_readers)}

    def _validate_thermal_parent(
            self, expected: Optional[Mapping[str, object]] = None,
    ) -> Dict[str, int]:
        paths = (
            self.plan.old_csv, self.plan.old_pid_file, self.plan.old_archive,
            self.plan.old_unclean_archive, self.plan.old_stale_pid_archive,
            self.plan.audit_binding,
        )
        parent = self.plan.old_csv.parent
        if any(path.parent != parent for path in paths):
            raise TransitionError("thermal evidence destinations changed parent")
        observed = directory_binding(parent)
        if observed["uid"] != self.plan.controller_uid or \
                observed["mode"] != 0o700:
            raise TransitionError(
                "thermal evidence parent is outside the controller trust boundary")
        if expected is not None and observed != dict(expected):
            raise TransitionError("thermal evidence parent binding changed")
        return observed

    def hard_preflight(self) -> Mapping[str, object]:
        value = verify_transition_plan(self.plan)
        if value.get("tools") != self.tools:
            raise TransitionError("backend tool ledger differs from the frozen plan")
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed after entry")
        # Discover an unavailable libc/kernel interface before the old sampler
        # can be signalled.  The actual archive still uses RENAME_NOREPLACE and
        # verifies its exact post-rename inode/content binding.
        renameat2_function()
        thermal_parent = self._validate_thermal_parent()
        for path in (
                self.plan.old_archive, self.plan.old_unclean_archive,
                self.plan.old_stale_pid_archive, self.plan.audit_binding):
            if path.exists() or path.is_symlink():
                raise TransitionError("reserved transition destination exists: %s" % path)
        source = file_binding(self.plan.old_source, with_hash=True)
        if source["sha256"] != self.plan.old_source_sha256 or \
                source["uid"] != self.plan.controller_uid or \
                source["mode"] & 0o222 or source["nlink"] != 1:
            raise TransitionError("old sampler source binding changed")
        csv_binding = file_binding(self.plan.old_csv, with_hash=False)
        if (csv_binding["device"], csv_binding["inode"],
                csv_binding["uid"], csv_binding["mode"],
                csv_binding["nlink"]) != (
                self.plan.old_csv_device, self.plan.old_csv_inode,
                0, 0o444, 1):
            raise TransitionError("live old CSV binding changed")
        pid_binding = file_binding(self.plan.old_pid_file, with_hash=True)
        if (pid_binding["device"], pid_binding["inode"],
                pid_binding["uid"], pid_binding["mode"],
                pid_binding["nlink"], pid_binding["size"],
                pid_binding["sha256"]) != (
                self.plan.old_pid_device, self.plan.old_pid_inode,
                0, 0o444, 1, len(str(self.plan.old_pid)) + 1,
                sha256_bytes((str(self.plan.old_pid) + "\n").encode("ascii"))):
            raise TransitionError("live old PID binding changed")
        identity = self._validate_old_identity()
        topology = self._verify_topology_and_occupancy()
        reader_seal = self._require_exclusive_old_readers()
        csv_readers = tuple(reader_seal["csv_readers"])
        self.old_preflight = {
            "csv": csv_binding, "csv_readers": list(csv_readers),
            "identity": identity, "pid_file": pid_binding,
            "source": source, "thermal_parent": thermal_parent,
            "tools": tools, "controller_interpreter": interpreter,
            "execute_runtime": runtime,
            "topology": topology,
        }
        return self.old_preflight

    def arm_recovery(self) -> Mapping[str, object]:
        # The controller writes the completed recovery_arm phase receipt before
        # calling stop_old.  This method performs one last no-side-effect seal.
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed before stop arm")
        identity = self._validate_old_identity()
        if self.old_preflight is None or not isinstance(
                self.old_preflight.get("thermal_parent"), dict):
            raise TransitionError("recovery arm lacks the thermal-parent seal")
        self._validate_thermal_parent(self.old_preflight["thermal_parent"])
        self._require_exclusive_old_readers()
        return {"armed_for_pid": self.plan.old_pid,
                "controller_interpreter": interpreter,
                "execute_runtime": runtime,
                "start_tick": identity["child"]["start_tick"],
                "tools": tools}

    def _old_pidfd_signal(self, signum: int) -> None:
        python_binding = self.tools["python3"]["binding"]
        command = (
            str(self.tools["env"]["path"]), "-i",
            "HOME=" + str(self.plan.root / "runtime-home"),
            "PATH=/usr/bin:/bin", "LC_ALL=C", "LANG=C", "TZ=UTC",
            "PYTHONDONTWRITEBYTECODE=1", str(self.tools["python3"]["path"]),
            "-I", "-S", "-B", "-c",
            PIDFD_SIGNAL_PROGRAM, str(self.plan.old_pid), self.plan.old_boot_id,
            str(self.plan.old_start_tick), self.plan.old_cmdline_sha256,
            str(self.plan.old_launcher_pid), str(self.plan.old_process_group),
            str(self.plan.old_session), str(self.plan.old_cpu),
            str(python_binding["device"]), str(python_binding["inode"]),
            str(signum),
        )
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), command,
            self._environment())
        if rc != 0 or stdout or stderr:
            raise TransitionError(
                "exact old-sampler pidfd helper failed for signal %d" % signum)

    def _old_identity_lives(self) -> bool:
        try:
            self._validate_old_child_identity()
        except (OSError, TransitionError, ValueError):
            return False
        return True

    def _old_launcher_identity_lives(self) -> bool:
        try:
            launcher = _proc_identity(self.plan.old_launcher_pid)
        except (OSError, TransitionError, ValueError):
            return False
        exact = {
            "affinity": self.plan.old_launcher_affinity,
            "argv": list(self.plan.old_launch_argv),
            "cmdline_sha256": self.plan.old_launcher_cmdline_sha256,
            "pid": self.plan.old_launcher_pid,
            "ppid": 1,
            "process_group": self.plan.old_process_group,
            "session": self.plan.old_session,
            "start_tick": self.plan.old_launcher_start_tick,
            "uids": list(self.plan.old_launcher_uids),
        }
        sudo_binding = self.tools["sudo"]["binding"]
        exact["executable"] = {
            "device": sudo_binding["device"],
            "inode": sudo_binding["inode"],
        }
        return all(launcher.get(key) == value for key, value in exact.items())

    def _old_launcher_pidfd_signal(self, signum: int) -> None:
        sudo_binding = self.tools["sudo"]["binding"]
        command = (
            str(self.tools["env"]["path"]), "-i",
            "HOME=" + str(self.plan.root / "runtime-home"),
            "PATH=/usr/bin:/bin", "LC_ALL=C", "LANG=C", "TZ=UTC",
            "PYTHONDONTWRITEBYTECODE=1", str(self.tools["python3"]["path"]),
            "-I", "-S", "-B", "-c", LAUNCHER_PIDFD_SIGNAL_PROGRAM,
            str(self.plan.old_launcher_pid), self.plan.old_boot_id,
            str(self.plan.old_launcher_start_tick),
            self.plan.old_launcher_cmdline_sha256,
            str(self.plan.old_process_group), str(self.plan.old_session),
            self.plan.old_launcher_affinity,
            ",".join(str(value) for value in self.plan.old_launcher_uids),
            str(sudo_binding["device"]), str(sudo_binding["inode"]),
            str(signum),
        )
        rc, stdout, stderr = self.p32.run_privileged_bounded(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), command,
            self._environment())
        if rc != 0 or stdout or stderr:
            raise TransitionError(
                "exact old-launcher pidfd helper failed for signal %d" % signum)

    @staticmethod
    def _raw_proc_pid_present(pid: int) -> bool:
        if type(pid) is not int or pid <= 1:
            raise TransitionError("process PID is outside its domain")
        try:
            Path("/proc/%d" % pid).stat()
        except FileNotFoundError:
            return False
        except OSError as exc:
            raise TransitionError("cannot prove raw process-path absence") from exc
        return True

    def _require_old_pid_absent(self) -> None:
        if self._raw_proc_pid_present(self.plan.old_pid):
            raise TransitionError(
                "old sampler PID remains with an unprovable identity")

    def _require_old_launcher_absent(self) -> None:
        if self._raw_proc_pid_present(self.plan.old_launcher_pid):
            raise TransitionError(
                "original old sudo launcher remains with an unprovable identity")

    def stop_old(self) -> Mapping[str, object]:
        if self.old_preflight is None:
            raise TransitionError("old stop lacks its preflight seal")
        tools_at_stop = self._verify_tools()
        interpreter_at_stop = verify_running_interpreter(tools_at_stop)
        self._require_exclusive_old_readers()
        self._old_pidfd_signal(signal.SIGTERM)
        forced = False
        graceful_deadline = min(self.deadline.absolute - 10.0,
                                time.monotonic() + 15.0)
        while self._old_identity_lives() and time.monotonic() < graceful_deadline:
            time.sleep(0.05)
        if self._old_identity_lives():
            forced = True
            self._old_pidfd_signal(signal.SIGKILL)
            kill_deadline = min(self.deadline.absolute - 5.0,
                                time.monotonic() + 5.0)
            while self._old_identity_lives() and time.monotonic() < kill_deadline:
                time.sleep(0.05)
        if self._old_identity_lives():
            raise TransitionError("exact old sampler survived bounded pidfd stop")
        self._require_old_pid_absent()
        parent_deadline = min(self.deadline.absolute - 5.0,
                              time.monotonic() + 5.0)
        while Path("/proc/%d" % self.plan.old_launcher_pid).exists() and \
                time.monotonic() < parent_deadline:
            time.sleep(0.05)
        if Path("/proc/%d" % self.plan.old_launcher_pid).exists():
            raise TransitionError("old sudo launcher survived child shutdown")
        if self._i2c_readers():
            raise TransitionError("I2C reader remained after old sampler stop")
        if not forced and (self.plan.old_pid_file.exists() or
                           self.plan.old_pid_file.is_symlink()):
            raise TransitionError("graceful old stop retained its PID file")
        self.stop_record = {
            "forced": forced, "i2c_readers_after": [],
            "controller_interpreter": interpreter_at_stop,
            "old_pid": self.plan.old_pid,
            "tools": tools_at_stop,
        }
        return self.stop_record

    def archive_old(self) -> Mapping[str, object]:
        if self.stop_record is None:
            raise TransitionError("old archive lacks a completed stop")
        if self._i2c_readers() or self._fuser(self.plan.old_csv):
            raise TransitionError("old evidence remains open at archive")
        binding = file_binding(self.plan.old_csv, with_hash=True)
        if self._fuser(self.plan.old_csv):
            raise TransitionError("old CSV gained a reader while being archived")
        preflight_csv = self.old_preflight["csv"] if self.old_preflight else None
        if not isinstance(preflight_csv, dict) or any(
                binding[key] != preflight_csv[key]
                for key in ("device", "inode", "uid", "gid", "mode", "nlink")):
            raise TransitionError("old CSV inode changed before archive")
        destination = self.plan.old_unclean_archive \
            if self.stop_record["forced"] else self.plan.old_archive
        rename_noreplace(
            self.plan.old_csv, destination, binding,
            parent_uid=self.plan.controller_uid)
        stale_pid = None
        if self.plan.old_pid_file.exists() or self.plan.old_pid_file.is_symlink():
            if not self.stop_record["forced"]:
                raise TransitionError("graceful stop unexpectedly retained PID evidence")
            stale_binding = file_binding(self.plan.old_pid_file, with_hash=True)
            rename_noreplace(
                self.plan.old_pid_file, self.plan.old_stale_pid_archive,
                stale_binding, parent_uid=self.plan.controller_uid)
            stale_pid = {"binding": stale_binding,
                         "path": str(self.plan.old_stale_pid_archive)}
        self.archive_record = {
            "binding": binding, "forced_stop": self.stop_record["forced"],
            "path": str(destination), "stale_pid": stale_pid,
        }
        return self.archive_record

    def start_candidate(self) -> Mapping[str, object]:
        if self.stop_record is None or self.stop_record.get("forced") is not False:
            raise TransitionError("candidate start is vetoed after forced old stop")
        if self.archive_record is None or \
                self.archive_record.get("forced_stop") is not False or \
                self._i2c_readers():
            raise TransitionError("candidate start lacks an empty-I2C archive seal")
        launch_topology = self._verify_topology_and_occupancy()
        if self._i2c_readers():
            raise TransitionError("I2C ownership appeared during CPU recheck")
        self.candidate_owner = self.p32.start_sampler(
            self.plan.root, self.design, 0)
        self.candidate_timing_start = time.monotonic()
        return {
            "argv": list(self.candidate_owner.identity["cmdline"]),
            "pid": self.candidate_owner.pid,
            "start_tick": self.candidate_owner.identity["start_tick"],
            "cmdline_sha256": self.candidate_owner.identity["cmdline_sha256"],
            "expected_output_owner_uid": self.plan.controller_uid,
            "expected_source_sha256": self.plan.candidate_sha256,
            "launch_topology": launch_topology,
        }

    def exercise_candidate(self) -> Mapping[str, object]:
        if self.candidate_owner is None or self.candidate_timing_start is None:
            raise TransitionError("candidate exercise lacks an owner")
        paths = self.p32._sampler_evidence_paths(
            self.plan.root, 0, final=False)
        last_count = 0
        exercise_deadline = min(self.deadline.normal, time.monotonic() + 12.0)
        while time.monotonic() < exercise_deadline:
            if not self.p32.process_identity_matches(
                    self.candidate_owner.identity, self.plan.candidate_cpu,
                    self.candidate_owner.csv_part):
                raise TransitionError("candidate identity changed during dry exercise")
            rows, samples = self.p32.stable_validated_thermal_rows(
                paths["csv"], paths["validation"], self.plan.candidate_sha256,
                self.plan.controller_uid, attempts=3)
            last_count = len(rows)
            if last_count >= 5:
                self.candidate_benchmark_end = time.monotonic()
                return {"paired_rows_before_stop": last_count,
                        "last_decision": samples[-1]["decision"]}
            time.sleep(0.05)
        raise TransitionDeadline(
            "candidate did not produce five paired rows: %d" % last_count)

    def stop_candidate(self) -> Mapping[str, object]:
        if self.candidate_owner is None or \
                self.candidate_timing_start is None or \
                self.candidate_benchmark_end is None:
            raise TransitionError("candidate stop lacks its exercise interval")
        self.candidate_terminal = self.p32.stop_sampler(
            self.candidate_owner, self.plan.root, self.design,
            self.candidate_timing_start, self.candidate_benchmark_end)
        self.candidate_cleanup_complete = True
        if self._i2c_readers():
            raise TransitionError("candidate stop retained an I2C reader")
        return self.candidate_terminal

    def accept_candidate(self) -> Mapping[str, object]:
        if self.candidate_terminal is None:
            raise TransitionError("candidate acceptance lacks terminal evidence")
        evidence, validation_raw = self.p32.validate_sampler_terminal_evidence(
            self.plan.root, self.design, 0, final=True)
        paths = self.p32._sampler_evidence_paths(self.plan.root, 0, final=True)
        segment_entries = self._candidate_artifact_roster(paths)
        raw_csv = stable_file_bytes(paths["csv"])
        rows = self.p32._parse_thermal_csv(
            raw_csv, allow_raw_dimm_faults=True)
        samples = self.p32._parse_thermal_validation(
            validation_raw, rows, self.plan.candidate_sha256,
            self.plan.controller_uid)
        if self.candidate_timing_start is None or \
                self.candidate_benchmark_end is None:
            raise TransitionError(
                "candidate acceptance lacks its measured interval")
        thermal = self.p32.validate_thermal_interval(
            raw_csv, self.candidate_timing_start, self.candidate_benchmark_end,
            validation_raw=validation_raw,
            expected_source_sha256=self.plan.candidate_sha256,
            expected_output_owner_uid=self.plan.controller_uid)
        if thermal != self.candidate_terminal.get("thermal_summary"):
            raise TransitionError(
                "candidate terminal thermal summary did not replay")
        summary = evidence["sampler_summary"]
        if len(rows) < 6 or any(sample["decision"] != "continue"
                                for sample in samples) or \
                any(sensor["reason"] != "ok" or not sensor["valid"] or
                    sensor["hot"] or sensor["attempt_errors"] != 0
                    for sample in samples
                    for sensor in sample["sensors"].values()) or \
                summary["dimm_invalid_samples_total"] != 0 or \
                summary["dimm_read_error_samples_total"] != 0 or \
                summary["dimm_attempt_errors_total"] != 0 or \
                summary["edac_ce_delta"] != 0 or summary["edac_ue_delta"] != 0:
            raise TransitionError("candidate dry samples are not strictly clean")
        if thermal["cpu_max_c"] >= 85.0 or \
                max(thermal["dimm_max_c"].values()) >= 84.0 or \
                thermal["max_gap_s"] > 2.25:
            raise TransitionError("candidate dry environmental gate failed")
        artifact_bindings = {
            name: {"binding": file_binding(path, with_hash=True),
                   "path": str(path)}
            for name, path in paths.items()
        }
        expected_hashes = {
            "csv": self.candidate_terminal["thermal_csv_sha256"],
            "receipt": self.candidate_terminal["thermal_sampler_evidence"]
                ["terminal_receipt"]["sha256"],
            "validation": evidence["validation_jsonl"]["sha256"],
        }
        if any(record["binding"]["sha256"] != expected_hashes[name] or
               record["binding"]["uid"] != 0 or
               record["binding"]["mode"] != 0o444 or
               record["binding"]["nlink"] != 1
               for name, record in artifact_bindings.items()):
            raise TransitionError("candidate artifact binding changed at acceptance")
        source_binding = file_binding(self.plan.sampler, with_hash=True)
        if source_binding["sha256"] != self.plan.candidate_sha256 or \
                source_binding["uid"] != self.plan.controller_uid or \
                source_binding["mode"] != 0o444 or \
                source_binding["nlink"] != 1:
            raise TransitionError("candidate source binding changed at acceptance")
        return {
            "artifact_bindings": artifact_bindings,
            "artifact_roster": [str(path) for path in segment_entries],
            "candidate_receipt_sha256":
                self.candidate_terminal["thermal_sampler_evidence"]
                    ["terminal_receipt"]["sha256"],
            "raw_sha256": self.candidate_terminal["thermal_csv_sha256"],
            "sample_count": len(rows),
            "source": {"binding": source_binding,
                       "path": str(self.plan.sampler)},
            "validation_sha256": evidence["validation_jsonl"]["sha256"],
        }

    def _candidate_artifact_roster(
            self, paths: Mapping[str, Path]) -> list[Path]:
        if set(paths) != {"csv", "receipt", "validation"} or \
                any(not isinstance(path, Path) for path in paths.values()):
            raise TransitionError("candidate final artifact paths are malformed")
        try:
            segment_entries = sorted(self.plan.root.joinpath("segments").iterdir())
        except OSError as exc:
            raise TransitionError("candidate final artifact roster is unreadable") \
                from exc
        expected_paths = set(paths.values())
        if set(segment_entries) != expected_paths or any(
                path.is_symlink() or not path.is_file()
                for path in segment_entries):
            raise TransitionError("candidate final artifact roster changed")
        return segment_entries

    def _owned_session_members(
            self, session_id: int, *, proc_root: Path = Path("/proc"),
    ) -> Tuple[int, ...]:
        if type(session_id) is not int or session_id <= 1:
            raise TransitionError("candidate session identity is malformed")
        members = []
        for process_path in proc_root.glob("[0-9]*"):
            try:
                pid = int(process_path.name)
                identity = _parse_proc_stat(
                    (process_path / "stat").read_bytes())
            except FileNotFoundError:
                continue
            except (OSError, TransitionError, ValueError) as exc:
                if process_path.exists():
                    raise TransitionError(
                        "cannot classify a live process during candidate cleanup") \
                        from exc
                continue
            if identity["session"] == session_id or \
                    identity["process_group"] == session_id:
                members.append(pid)
        return tuple(sorted(set(members)))

    def _cleanup_proof_pause(self) -> None:
        # This separation is an evidence requirement, not a recovery wait.
        # Even an exhausted emergency budget must not collapse two probes into
        # the same instant and manufacture a vacuous stable-empty proof.
        time.sleep(0.05)

    @staticmethod
    def _candidate_pid_present(pid: int, *,
                               proc_root: Path = Path("/proc")) -> bool:
        if type(pid) is not int or pid <= 1:
            raise TransitionError("candidate PID identity is malformed")
        try:
            (proc_root / str(pid)).stat()
        except FileNotFoundError:
            return False
        except OSError as exc:
            raise TransitionError(
                "cannot prove candidate PID absence during cleanup") from exc
        return True

    def _candidate_cleanup_observation(
            self, owner: object) -> Dict[str, object]:
        returncode = getattr(owner.process, "poll")()
        identity_live = self.p32.process_identity_matches(
            owner.identity, self.plan.candidate_cpu, owner.csv_part)
        session_members = self._owned_session_members(owner.process.pid)
        i2c_readers = self._i2c_readers()
        final_paths = self.p32._sampler_evidence_paths(
            self.plan.root, owner.segment, final=True)
        final_csv = final_paths.get("csv")
        if not isinstance(final_csv, Path):
            raise TransitionError("candidate final CSV path is malformed")
        csv_readers = set()
        for path in {owner.csv_part, final_csv}:
            try:
                info = path.lstat()
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(info.st_mode):
                raise TransitionError(
                    "candidate cleanup CSV path is not a regular file")
            csv_readers.update(self._fuser(path))
        return {
            "candidate_pid_present": self._candidate_pid_present(owner.pid),
            "csv_readers": sorted(csv_readers),
            "exact_identity_live": identity_live,
            "i2c_readers": list(i2c_readers),
            "launcher_returncode": returncode,
            "session_members": list(session_members),
        }

    @staticmethod
    def _cleanup_observation_is_empty(
            observation: Mapping[str, object]) -> bool:
        return observation.get("launcher_returncode") is not None and \
            observation.get("candidate_pid_present") is False and \
            observation.get("exact_identity_live") is False and \
            observation.get("session_members") == [] and \
            observation.get("i2c_readers") == [] and \
            observation.get("csv_readers") == []

    def cleanup_candidate(self) -> Mapping[str, object]:
        owner = self.candidate_owner
        self.candidate_cleanup_complete = False
        kill_error: Optional[BaseException] = None
        if owner is not None:
            try:
                # The privileged launcher may already have exited while the
                # exact root-owned session remains alive.  Always direct the
                # identity-verifying helper at the sealed owned session.
                self.p32._kill_owned_sampler_session(
                    owner, self.plan.root, self.design)
            except BaseException as exc:
                # The helper can report a transport/reap error after it has
                # already cleared the exact session.  Do not infer either
                # success or failure: independently prove stable emptiness.
                kill_error = exc
            try:
                first = self._candidate_cleanup_observation(owner)
                self._cleanup_proof_pause()
                second = self._candidate_cleanup_observation(owner)
                if not self._cleanup_observation_is_empty(first) or \
                        not self._cleanup_observation_is_empty(second):
                    raise TransitionError(
                        "candidate cleanup could not prove stable empty ownership: "
                        "%r then %r" % (first, second))
            except BaseException as proof_error:
                if kill_error is not None:
                    raise TransitionError(
                        "candidate kill helper failed and independent empty-"
                        "ownership proof failed: helper=%s; proof=%s" %
                        (kill_error, proof_error)) from proof_error
                raise
            self.candidate_cleanup_complete = True
            if kill_error is not None:
                raise TransitionError(
                    "candidate kill helper failed after independently proven "
                    "empty ownership: %s" % kill_error) from kill_error
            return {
                "empty_ownership_observations": [first, second],
                "i2c_readers_after": [],
                "owner_created": True,
            }

        first_readers = self._i2c_readers()
        self._cleanup_proof_pause()
        second_readers = self._i2c_readers()
        allowed = ((), (self.plan.old_pid,))
        if first_readers not in allowed or second_readers not in allowed:
            raise TransitionError(
                "candidate cleanup found an unowned I2C reader: "
                "%r then %r" % (first_readers, second_readers))
        self.candidate_cleanup_complete = True
        return {"i2c_readers_after": list(second_readers),
                "owner_created": False,
                "stable_observations": 2}

    def _quiesce_old_for_recovery(self) -> None:
        if self._old_identity_lives():
            # A stop helper can fail after delivering SIGTERM.  Complete the
            # same exact-identity stop before any restart rather than assume the
            # old process is healthy or launch a second reader.
            try:
                self._old_pidfd_signal(signal.SIGTERM)
            except BaseException:
                pass
            wait_deadline = self._recovery_wait_deadline(
                5.0, minimum_safety_wait_s=1.0)
            while self._old_identity_lives() and \
                    self._recovery_now() < wait_deadline:
                time.sleep(0.05)
            if self._old_identity_lives():
                self._old_pidfd_signal(signal.SIGKILL)
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._old_identity_lives() and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        if self._old_identity_lives():
            raise TransitionError("old sampler could not be quiesced for recovery")
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._raw_proc_pid_present(self.plan.old_pid) and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        self._require_old_pid_absent()
        if self._old_launcher_identity_lives():
            try:
                self._old_launcher_pidfd_signal(signal.SIGTERM)
            except BaseException:
                pass
            wait_deadline = self._recovery_wait_deadline(
                5.0, minimum_safety_wait_s=1.0)
            while self._old_launcher_identity_lives() and \
                    self._recovery_now() < wait_deadline:
                time.sleep(0.05)
            if self._old_launcher_identity_lives():
                self._old_launcher_pidfd_signal(signal.SIGKILL)
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._old_launcher_identity_lives() and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        wait_deadline = self._recovery_wait_deadline(
            5.0, minimum_safety_wait_s=1.0)
        while self._raw_proc_pid_present(self.plan.old_launcher_pid) and \
                self._recovery_now() < wait_deadline:
            time.sleep(0.05)
        self._require_old_launcher_absent()
        readers = self._i2c_readers()
        if readers:
            raise TransitionError("I2C is not empty before old recovery")

    def _ensure_recovery_archive(
            self, archive: Optional[Mapping[str, object]]) -> Dict[str, object]:
        def validate_csv_binding(binding: Mapping[str, object]) -> None:
            preflight = getattr(self, "old_preflight", None)
            expected = preflight.get("csv") \
                if isinstance(preflight, dict) else None
            if isinstance(expected, dict) and any(
                    binding.get(key) != expected.get(key)
                    for key in (
                        "device", "inode", "uid", "gid", "mode", "nlink")):
                raise TransitionError(
                    "recovery archive no longer binds the preflight CSV inode")

        def preserve_stale_pid() -> Optional[Dict[str, object]]:
            canonical_exists = self.plan.old_pid_file.exists() or \
                self.plan.old_pid_file.is_symlink()
            archived_exists = self.plan.old_stale_pid_archive.exists() or \
                self.plan.old_stale_pid_archive.is_symlink()
            if canonical_exists and archived_exists:
                raise TransitionError("stale PID archive collision during recovery")
            if canonical_exists:
                pid_binding = file_binding(
                    self.plan.old_pid_file, with_hash=True)
                preflight = getattr(self, "old_preflight", None)
                expected_pid = preflight.get("pid_file") \
                    if isinstance(preflight, dict) else None
                if isinstance(expected_pid, dict) and \
                        pid_binding != expected_pid:
                    raise TransitionError(
                        "recovery stale PID no longer binds the preflight inode")
                rename_noreplace(
                    self.plan.old_pid_file, self.plan.old_stale_pid_archive,
                    pid_binding, parent_uid=self.plan.controller_uid)
                return {"binding": pid_binding,
                        "path": str(self.plan.old_stale_pid_archive)}
            if archived_exists:
                return {
                    "binding": file_binding(
                        self.plan.old_stale_pid_archive, with_hash=True),
                    "path": str(self.plan.old_stale_pid_archive),
                }
            return None

        if archive is not None:
            path = Path(str(archive["path"]))
            forced = archive.get("forced_stop")
            expected_path = self.plan.old_unclean_archive \
                if forced is True else self.plan.old_archive
            if type(forced) is not bool or path != expected_path:
                raise TransitionError("pre-dry archive path changed during recovery")
            other_path = self.plan.old_archive \
                if forced else self.plan.old_unclean_archive
            if other_path.exists() or other_path.is_symlink():
                raise TransitionError(
                    "alternate pre-dry archive appeared during recovery")
            binding = file_binding(path, with_hash=True)
            if binding != archive.get("binding"):
                raise TransitionError("pre-dry archive changed during recovery")
            validate_csv_binding(binding)
            result = dict(archive)
            if not forced and (result.get("stale_pid") is not None or
                               self.plan.old_pid_file.exists() or
                               self.plan.old_pid_file.is_symlink() or
                               self.plan.old_stale_pid_archive.exists() or
                               self.plan.old_stale_pid_archive.is_symlink()):
                raise TransitionError(
                    "graceful archive unexpectedly has stale PID evidence")
            result["stale_pid"] = result.get("stale_pid") or \
                (preserve_stale_pid() if forced else None)
            return result
        discovered = [path for path in (
            self.plan.old_archive, self.plan.old_unclean_archive)
            if path.exists() or path.is_symlink()]
        if len(discovered) > 1:
            raise TransitionError("multiple pre-dry archives exist during recovery")
        if discovered:
            path = discovered[0]
            binding = file_binding(path, with_hash=True)
            validate_csv_binding(binding)
            forced = path == self.plan.old_unclean_archive
            if not forced and (self.plan.old_pid_file.exists() or
                               self.plan.old_pid_file.is_symlink() or
                               self.plan.old_stale_pid_archive.exists() or
                               self.plan.old_stale_pid_archive.is_symlink()):
                raise TransitionError(
                    "graceful archive unexpectedly has stale PID evidence")
            record = {
                "binding": binding,
                "forced_stop": forced,
                "path": str(path),
                "stale_pid": preserve_stale_pid() if forced else None,
            }
            self.archive_record = record
            return record
        if not (self.plan.old_csv.exists() or self.plan.old_csv.is_symlink()):
            raise TransitionError("recovery lacks both old CSV and archive")
        binding = file_binding(self.plan.old_csv, with_hash=True)
        validate_csv_binding(binding)
        if self._fuser(self.plan.old_csv):
            raise TransitionError("old CSV gained a reader during recovery archive")
        destination = self.plan.old_unclean_archive
        rename_noreplace(
            self.plan.old_csv, destination, binding,
            parent_uid=self.plan.controller_uid)
        stale_pid = preserve_stale_pid()
        record = {"binding": binding, "forced_stop": True,
                  "path": str(destination), "stale_pid": stale_pid}
        self.archive_record = record
        return record

    def _replacement_launcher_roster(
            self, restored: Mapping[str, object]) -> list[Dict[str, object]]:
        session_id = restored.get("launcher_session")
        child_pid = restored.get("pid")
        launcher_tick = restored.get("launcher_start_tick")
        expected_command = self._replacement_launch_command()
        expected_command_sha256 = sha256_bytes(
            b"\0".join(value.encode("ascii")
                        for value in expected_command) + b"\0")
        if type(session_id) is not int or session_id <= 1 or \
                type(child_pid) is not int or child_pid <= 1 or \
                type(launcher_tick) is not int or launcher_tick < 0 or \
                restored.get("launcher_command") != list(expected_command) or \
                restored.get("launcher_command_sha256") != \
                    expected_command_sha256:
            raise TransitionError("replacement launcher identity is malformed")
        members = self._owned_session_members(session_id)
        allowed = {session_id, child_pid}
        if child_pid not in members or any(pid not in allowed for pid in members):
            raise TransitionError(
                "replacement launcher session roster has an unknown member")
        roster = []
        for pid in members:
            identity = _proc_identity(pid)
            if identity["session"] != session_id or \
                    identity["process_group"] != session_id:
                raise TransitionError(
                    "replacement launcher roster identity changed")
            if pid == child_pid:
                python_binding = self.tools["python3"]["binding"]
                if identity["argv"] != list(self.plan.replacement_old_argv) or \
                        identity["cmdline_sha256"] != \
                            self.plan.replacement_old_cmdline_sha256 or \
                        identity["start_tick"] != restored.get("start_tick") or \
                        identity["uids"] != [0, 0, 0, 0] or \
                        identity["affinity"] != str(self.plan.old_cpu) or \
                        identity.get("ppid") not in (session_id, 1) or \
                        identity.get("executable") != {
                            "device": python_binding["device"],
                            "inode": python_binding["inode"]}:
                    raise TransitionError(
                        "replacement child changed in launcher roster")
            else:
                initial_launcher = restored.get("launcher_identity")
                if launcher_tick == 0 or not isinstance(initial_launcher, dict) or \
                        identity != initial_launcher:
                    raise TransitionError(
                        "replacement sudo launcher changed after creation")
                self._validate_replacement_launcher_identity(
                    initial_launcher, session_id, launcher_tick)
            roster.append(identity)
        return roster

    def _replacement_launch_command(self) -> Tuple[str, ...]:
        environment = execute_environment(self.plan)
        return (
            str(self.tools["sudo"]["path"]), "-n", "-b",
            str(self.tools["env"]["path"]), "-i",
            *("%s=%s" % (name, environment[name])
              for name in ("HOME", "PATH", "LANG", "LC_ALL", "TZ")),
            str(self.tools["taskset"]["path"]), "--cpu-list",
            str(self.plan.old_cpu), *self.plan.replacement_old_argv,
        )

    def _validate_replacement_launcher_identity(
            self, identity: Mapping[str, object], session_id: int,
            start_tick: int) -> None:
        command = self._replacement_launch_command()
        expected_hash = sha256_bytes(
            b"\0".join(value.encode("ascii") for value in command) + b"\0")
        exact = {
            "argv": list(command),
            "cmdline_sha256": expected_hash,
            "pid": session_id,
            "process_group": session_id,
            "session": session_id,
            "start_tick": start_tick,
            "uids": [self.plan.controller_uid, 0, 0, 0],
        }
        sudo_binding = self.tools["sudo"]["binding"]
        exact["executable"] = {
            "device": sudo_binding["device"],
            "inode": sudo_binding["inode"],
        }
        if any(identity.get(key) != value for key, value in exact.items()) or \
                identity.get("ppid") not in (os.getpid(), 1) or \
                not isinstance(identity.get("affinity"), str) or \
                not identity["affinity"]:
            raise TransitionError(
                "replacement sudo launcher identity is not exact")

    def _launch_old(self) -> Dict[str, object]:
        self._require_old_launcher_absent()
        if self._i2c_readers() or any(
                path.exists() or path.is_symlink()
                for path in (self.plan.old_csv, self.plan.old_pid_file)):
            raise TransitionError("old restart preconditions are not empty")
        source_binding = file_binding(self.plan.old_source, with_hash=True)
        expected_source = self.old_preflight.get("source") \
            if isinstance(self.old_preflight, dict) else None
        if not isinstance(expected_source, dict) or \
                source_binding != expected_source or \
                source_binding["sha256"] != self.plan.old_source_sha256:
            raise TransitionError("old source changed before restart")
        environment = execute_environment(self.plan)
        command = self._replacement_launch_command()
        launcher_start_tick = 0
        launcher_error: Optional[BaseException] = None
        boot_id = Path("/proc/sys/kernel/random/boot_id").read_text(
            encoding="ascii").strip()
        process = subprocess.Popen(
            command, env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            start_new_session=True, close_fds=True)
        last_error: Optional[BaseException] = None
        try:
            try:
                launcher_start_tick = capture_owned_session_leader(process.pid)
            except BaseException as exc:
                launcher_start_tick = 0
                launcher_error = exc
            if launcher_error is not None:
                raise launcher_error
            bootstrap_deadline = self._recovery_wait_deadline(
                20.0, minimum_safety_wait_s=5.0)
            while self._recovery_now() < bootstrap_deadline:
                try:
                    raw_pid = stable_file_bytes(self.plan.old_pid_file)
                    if re.fullmatch(rb"[1-9][0-9]*\n", raw_pid) is None:
                        raise TransitionError("replacement old PID is noncanonical")
                    pid = int(raw_pid)
                    identity = self.p32.capture_process_identity(
                        pid, self.plan.old_cpu, self.plan.old_csv)
                    raw = stable_file_bytes(self.plan.old_csv, attempts=3)
                    rows = self.p32._parse_thermal_csv(raw)
                    if not rows:
                        raise TransitionError(
                            "replacement old CSV has no complete sample")
                    if identity["cmdline"] != \
                            list(self.plan.replacement_old_argv) or \
                            identity["cmdline_sha256"] != \
                                self.plan.replacement_old_cmdline_sha256 or \
                            identity["uid"] != 0 or \
                            identity["boot_id"] != boot_id or \
                            identity["session_id"] != process.pid or \
                            identity["process_group"] != process.pid or \
                            self._i2c_readers() != (pid,) or \
                            self._fuser(self.plan.old_csv) != (pid,):
                        raise TransitionError(
                            "replacement old sampler identity is not exact")
                    csv_binding = file_binding(
                        self.plan.old_csv, with_hash=False)
                    pid_binding = file_binding(
                        self.plan.old_pid_file, with_hash=True)
                    if (csv_binding["uid"], csv_binding["mode"],
                            csv_binding["nlink"]) != (0, 0o444, 1) or \
                            (pid_binding["uid"], pid_binding["mode"],
                             pid_binding["nlink"], pid_binding["sha256"]) != (
                                0, 0o444, 1,
                                sha256_bytes(raw_pid)):
                        raise TransitionError(
                            "replacement old sampler evidence binding is not exact")
                    restored = {
                        "boot_id": identity["boot_id"],
                        "cmdline": identity["cmdline"],
                        "cmdline_sha256": identity["cmdline_sha256"],
                        "csv_initial_size": csv_binding["size"],
                        "csv_live_identity": {
                            key: csv_binding[key] for key in (
                                "device", "gid", "inode", "mode", "nlink",
                                "uid")},
                        "csv_path": str(self.plan.old_csv),
                        "first_sample_monotonic_s": rows[0]["monotonic_s"],
                        "launcher_session": process.pid,
                        "launcher_start_tick": launcher_start_tick,
                        "pid": pid, "pid_binding": pid_binding,
                        "pid_path": str(self.plan.old_pid_file),
                        "source_sha256": self.plan.old_source_sha256,
                        "start_tick": identity["start_tick"],
                    }
                    launcher_returncode = process.poll()
                    if launcher_returncode not in (None, 0):
                        raise TransitionError(
                            "replacement sudo launcher exited nonzero")
                    try:
                        launcher_identity = _proc_identity(process.pid)
                    except (FileNotFoundError, ProcessLookupError):
                        launcher_identity = None
                    except TransitionError:
                        launcher_returncode = process.poll()
                        if launcher_returncode is None:
                            raise
                        launcher_identity = None
                    if launcher_returncode is not None and \
                            launcher_identity is not None:
                        raise TransitionError(
                            "reaped replacement launcher still has an identity")
                    if launcher_identity is not None:
                        if launcher_start_tick == 0:
                            raise TransitionError(
                                "replacement sudo launcher start identity was lost")
                        self._validate_replacement_launcher_identity(
                            launcher_identity, process.pid,
                            launcher_start_tick)
                    restored["launcher_command"] = list(command)
                    restored["launcher_command_sha256"] = sha256_bytes(
                        b"\0".join(value.encode("ascii")
                                    for value in command) + b"\0")
                    restored["launcher_identity"] = launcher_identity
                    restored["launcher_returncode_at_accept"] = \
                        launcher_returncode
                    restored["launcher_roster"] = \
                        self._replacement_launcher_roster(restored)
                    return restored
                except (OSError, TransitionError, ValueError) as exc:
                    last_error = exc
                    time.sleep(0.05)
        except BaseException:
            try:
                self.p32._kill_owned_process_session(
                    process, launcher_start_tick, boot_id,
                    self.plan.root, self.design)
            except BaseException as cleanup:
                raise TransitionError(
                    "old restart bootstrap and exact-session cleanup failed: %s" %
                    cleanup) from cleanup
            raise
        try:
            self.p32._kill_owned_process_session(
                process, launcher_start_tick, boot_id,
                self.plan.root, self.design)
        except BaseException as cleanup:
            raise TransitionError(
                "old restart timeout and exact-session cleanup failed: %s" %
                cleanup) from cleanup
        raise TransitionError(
            "replacement old sampler bootstrap failed: %s" % last_error)

    def restore_old(
            self, archive: Optional[Mapping[str, object]]) -> Mapping[str, object]:
        if not self.candidate_cleanup_complete:
            # Still independently enforce empty I2C below; this flag makes an
            # omitted cleanup attempt visible and fail closed.
            raise TransitionError("old recovery was entered before candidate cleanup")
        self._quiesce_old_for_recovery()
        sealed_archive = self._ensure_recovery_archive(archive)
        restored = self._launch_old()
        restored["archive_record"] = sealed_archive
        restored["archived_pre_dry_sha256"] = \
            sealed_archive["binding"]["sha256"]
        restored["archived_pre_dry_path"] = sealed_archive["path"]
        return restored

    def _revalidate_restored_old(
            self, restored: Mapping[str, object]) -> Dict[str, object]:
        pid = restored.get("pid")
        if type(pid) is not int or pid <= 1:
            raise TransitionError("restored old PID is malformed")
        identity = self.p32.capture_process_identity(
            pid, self.plan.old_cpu, self.plan.old_csv)
        if identity["cmdline"] != list(self.plan.replacement_old_argv) or \
                identity["cmdline_sha256"] != \
                    self.plan.replacement_old_cmdline_sha256 or \
                identity["cmdline_sha256"] != restored.get("cmdline_sha256") or \
                identity["start_tick"] != restored.get("start_tick") or \
                identity["boot_id"] != restored.get("boot_id") or \
                identity["session_id"] != restored.get("launcher_session") or \
                identity["process_group"] != restored.get("launcher_session") or \
                self._i2c_readers() != (pid,) or \
                self._fuser(self.plan.old_csv) != (pid,):
            raise TransitionError("restored old sampler changed before audit binding")
        csv_raw = stable_file_bytes(self.plan.old_csv, attempts=3)
        csv_rows = self.p32._parse_thermal_csv(csv_raw)
        csv_binding = file_binding(self.plan.old_csv, with_hash=False)
        csv_identity = {
            key: csv_binding[key] for key in (
                "device", "gid", "inode", "mode", "nlink", "uid")}
        if csv_identity != restored.get("csv_live_identity"):
            raise TransitionError("restored old CSV inode changed before audit binding")
        initial_size = restored.get("csv_initial_size")
        if type(initial_size) is not int or initial_size <= 0 or \
                csv_binding["size"] < initial_size:
            raise TransitionError("restored old CSV size regressed before audit binding")
        if not csv_rows or len(csv_raw) < initial_size or \
                csv_binding["size"] < len(csv_raw):
            raise TransitionError("restored old CSV content changed before audit binding")
        csv_times = [float(row["monotonic_s"]) for row in csv_rows]
        if csv_times != sorted(csv_times) or \
                any(right <= left or right - left > 2.25
                    for left, right in zip(csv_times, csv_times[1:])) or \
                not 0.0 <= time.monotonic() - csv_times[-1] <= 2.25:
            raise TransitionError("restored old CSV cadence is not live")
        pid_binding = file_binding(self.plan.old_pid_file, with_hash=True)
        source_binding = file_binding(self.plan.old_source, with_hash=True)
        expected_source = self.old_preflight.get("source") \
            if isinstance(self.old_preflight, dict) else None
        if pid_binding != restored.get("pid_binding") or \
                not isinstance(expected_source, dict) or \
                source_binding != expected_source or \
                source_binding["sha256"] != self.plan.old_source_sha256:
            raise TransitionError("restored old PID/source binding changed")
        return {"csv_current_size": csv_binding["size"],
                "csv_stable_sample_count": len(csv_rows),
                "csv_stable_sha256": sha256_bytes(csv_raw),
                "csv_stable_size": len(csv_raw),
                "csv_stable_last_monotonic_s":
                    csv_rows[-1]["monotonic_s"],
                "identity": identity,
                "launcher_roster":
                    self._replacement_launcher_roster(restored)}

    def _replay_external_state(
            self, candidate_accept: Optional[Mapping[str, object]],
            archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object],
    ) -> Dict[str, object]:
        """Re-open every safety-critical artifact and live ownership roster."""
        frozen_plan = verify_transition_plan(self.plan)
        if frozen_plan.get("tools") != self.tools:
            raise TransitionError(
                "frozen transition plan changed during external replay")
        plan_binding = file_binding(self.plan.plan_receipt, with_hash=True)
        if plan_binding["uid"] != self.plan.controller_uid or \
                plan_binding["mode"] != 0o444 or plan_binding["nlink"] != 1:
            raise TransitionError(
                "transition plan receipt binding changed during replay")
        archive_record = archive if archive is not None else \
            restored.get("archive_record")
        if not isinstance(archive_record, dict) or \
                not isinstance(archive_record.get("path"), str) or \
                not isinstance(archive_record.get("binding"), dict):
            raise TransitionError("external replay lacks the pre-dry archive")
        archive_path = Path(str(archive_record["path"]))
        forced_stop = archive_record.get("forced_stop")
        expected_archive_path = self.plan.old_unclean_archive \
            if forced_stop is True else self.plan.old_archive
        if type(forced_stop) is not bool or \
                archive_path != expected_archive_path:
            raise TransitionError("external replay archive destination changed")
        alternate_archive = self.plan.old_archive \
            if forced_stop else self.plan.old_unclean_archive
        if alternate_archive.exists() or alternate_archive.is_symlink():
            raise TransitionError(
                "alternate pre-dry archive appeared during external replay")
        archive_binding = file_binding(archive_path, with_hash=True)
        if archive_binding != archive_record["binding"] or \
                archive_binding["uid"] != 0 or \
                archive_binding["mode"] != 0o444 or \
                archive_binding["nlink"] != 1 or \
                restored.get("archived_pre_dry_path") != str(archive_path) or \
                restored.get("archived_pre_dry_sha256") != \
                    archive_binding["sha256"]:
            raise TransitionError("pre-dry archive changed during external replay")
        stale_result = None
        stale_pid = archive_record.get("stale_pid")
        if stale_pid is not None:
            if not forced_stop:
                raise TransitionError(
                    "graceful archive gained stale PID evidence")
            if not isinstance(stale_pid, dict) or \
                    stale_pid.get("path") != \
                        str(self.plan.old_stale_pid_archive) or \
                    not isinstance(stale_pid.get("binding"), dict):
                raise TransitionError("external replay stale PID receipt changed")
            stale_binding = file_binding(
                self.plan.old_stale_pid_archive, with_hash=True)
            if stale_binding != stale_pid["binding"] or \
                    stale_binding["uid"] != 0 or \
                    stale_binding["mode"] != 0o444 or \
                    stale_binding["nlink"] != 1:
                raise TransitionError("stale PID archive changed during replay")
            stale_result = {"binding": stale_binding,
                            "path": str(self.plan.old_stale_pid_archive)}
        elif self.plan.old_stale_pid_archive.exists() or \
                self.plan.old_stale_pid_archive.is_symlink():
            raise TransitionError(
                "unbound stale PID archive appeared during replay")
        source = file_binding(self.plan.old_source, with_hash=True)
        expected_source = self.old_preflight.get("source") \
            if isinstance(self.old_preflight, dict) else None
        if not isinstance(expected_source, dict) or source != expected_source or \
                source["sha256"] != self.plan.old_source_sha256:
            raise TransitionError("old source changed during external replay")
        candidate_result = self._replay_candidate_accept(
            candidate_accept, restored.get("pid"))
        tools = self._verify_tools()
        interpreter = verify_running_interpreter(tools)
        runtime = verify_execute_runtime(self.plan)
        if runtime != self.execute_runtime:
            raise TransitionError("execute runtime changed during final replay")
        return {
            "archive": {"binding": archive_binding,
                        "path": str(archive_path)},
            "candidate": candidate_result,
            "controller_interpreter": interpreter,
            "execute_runtime": runtime,
            "old_source": {"binding": source,
                           "path": str(self.plan.old_source)},
            "restored_old": self._revalidate_restored_old(restored),
            "stale_pid": stale_result,
            "thermal_parent": self._validate_thermal_parent(
                self.old_preflight.get("thermal_parent")
                if isinstance(self.old_preflight, dict) else None),
            "tools": tools,
            "transition_plan": {
                "binding": plan_binding,
                "path": str(self.plan.plan_receipt),
                "self_sha256_excluding_field":
                    frozen_plan["self_sha256_excluding_field"],
            },
        }

    def _replay_candidate_accept(
            self, candidate_accept: Optional[Mapping[str, object]],
            restored_pid: object,
    ) -> Optional[Dict[str, object]]:
        if self.candidate_owner is None:
            if candidate_accept is not None:
                raise TransitionError(
                    "candidate evidence exists without a launched owner")
            return None
        if not self.candidate_cleanup_complete or \
                type(restored_pid) is not int or restored_pid <= 1:
            raise TransitionError(
                "launched candidate lacks a completed ownership cleanup")
        current: Optional[Dict[str, object]] = None
        if candidate_accept is not None:
            current = dict(self.accept_candidate())
            if current != dict(candidate_accept):
                raise TransitionError(
                    "candidate final evidence changed after acceptance")
        observation = self._candidate_cleanup_observation(self.candidate_owner)
        # This replay runs after the old sampler has been restored, so global
        # I2C ownership must now be exactly that restored PID.  Candidate
        # identity/session/CSV ownership must remain absent independently.
        if observation.get("launcher_returncode") is None or \
                observation.get("candidate_pid_present") is not False or \
                observation.get("exact_identity_live") is not False or \
                observation.get("session_members") != [] or \
                observation.get("csv_readers") != [] or \
                observation.get("i2c_readers") != [restored_pid]:
            raise TransitionError(
                "candidate ownership reappeared during final replay")
        if candidate_accept is None:
            return {"acceptance": None, "ownership": observation}
        assert current is not None
        current["ownership"] = observation
        return current

    def publish_audit_binding(
            self, archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object], dry_accepted: bool,
            candidate_accept: Optional[Mapping[str, object]],
            receipt_chain_prefix: Mapping[str, object],
    ) -> Mapping[str, object]:
        if type(dry_accepted) is not bool or \
                (dry_accepted and not isinstance(candidate_accept, dict)):
            raise TransitionError(
                "audit acceptance claim lacks retained candidate evidence")
        if not isinstance(receipt_chain_prefix, dict) or \
                set(receipt_chain_prefix) != {"count", "head_sha256", "roster"} or \
                type(receipt_chain_prefix.get("count")) is not int or \
                receipt_chain_prefix["count"] <= 0 or \
                not isinstance(receipt_chain_prefix.get("head_sha256"), str) or \
                SHA256_RE.fullmatch(receipt_chain_prefix["head_sha256"]) is None or \
                not isinstance(receipt_chain_prefix.get("roster"), list) or \
                len(receipt_chain_prefix["roster"]) != \
                    receipt_chain_prefix["count"]:
            raise TransitionError("audit receipt-chain prefix is malformed")
        prepublication_replay = self._replay_external_state(
            candidate_accept, archive, restored)
        archive_record = archive if archive is not None else \
            restored.get("archive_record")
        if not isinstance(archive_record, dict) or \
                not isinstance(archive_record.get("path"), str) or \
                not isinstance(archive_record.get("binding"), dict):
            raise TransitionError("audit binding lacks the pre-dry archive")
        archive_path = str(archive_record["path"])
        if archive_path not in {
                str(self.plan.old_archive), str(self.plan.old_unclean_archive)}:
            raise TransitionError("audit archive destination changed")
        archive_binding = file_binding(Path(archive_path), with_hash=True)
        if archive_binding != archive_record["binding"] or \
                archive_binding["uid"] != 0 or \
                archive_binding["mode"] != 0o444 or \
                archive_binding["nlink"] != 1 or \
                restored.get("archived_pre_dry_path") != archive_path or \
                restored.get("archived_pre_dry_sha256") != \
                    archive_binding["sha256"]:
            raise TransitionError("pre-dry archive changed before audit binding")
        stale_pid = archive_record.get("stale_pid")
        if stale_pid is not None:
            if not isinstance(stale_pid, dict) or \
                    stale_pid.get("path") != \
                        str(self.plan.old_stale_pid_archive) or \
                    not isinstance(stale_pid.get("binding"), dict):
                raise TransitionError("stale PID archive receipt changed")
            stale_binding = file_binding(
                self.plan.old_stale_pid_archive, with_hash=True)
            if stale_binding != stale_pid["binding"] or \
                    stale_binding["uid"] != 0 or \
                    stale_binding["mode"] != 0o444 or \
                    stale_binding["nlink"] != 1:
                raise TransitionError("stale PID archive changed before audit binding")
        # Keep the live identity check after potentially lengthy archive hashing,
        # immediately before publishing the future-auditor handoff.
        tools_after_restore = self._verify_tools()
        interpreter_after_restore = verify_running_interpreter(
            tools_after_restore)
        live_revalidation = self._revalidate_restored_old(restored)
        value = sealed_record(AUDIT_BINDING_SCHEMA, {
            "archived_pre_dry": {
                "binding": archive_binding, "path": archive_path},
            "candidate_dry_accepted": dry_accepted,
            "candidate_accept": dict(candidate_accept)
                if candidate_accept is not None else None,
            "created_utc": utc_now(),
            "controller_interpreter": interpreter_after_restore,
            "live_revalidation": live_revalidation,
            "live_old_sampler": dict(restored),
            "prepublication_replay": prepublication_replay,
            "receipt_chain_prefix": dict(receipt_chain_prefix),
            "tools": tools_after_restore,
            "transition_id": self.plan.transition_id,
        })
        binding = write_new(
            self.plan.audit_binding, canonical_json(value),
            owner_uid=self.plan.controller_uid)
        replay = load_canonical(
            self.plan.audit_binding, AUDIT_BINDING_SCHEMA,
            "future audit binding")
        if replay != value:
            raise TransitionError("future audit binding did not replay")
        tools_after_audit = verify_tool_records(replay.get("tools"))
        interpreter_after_audit = verify_running_interpreter(
            tools_after_audit)
        if tools_after_audit != tools_after_restore or \
                tools_after_audit != self.tools or \
                interpreter_after_audit != interpreter_after_restore:
            raise TransitionError(
                "tool identities changed after audit publication")
        return {"binding": binding, "path": str(self.plan.audit_binding),
                "controller_interpreter_after_audit":
                    interpreter_after_audit,
                "tools_after_audit": tools_after_audit, "value": value}

    def final_replay(
            self, candidate_accept: Optional[Mapping[str, object]],
            archive: Optional[Mapping[str, object]],
            restored: Mapping[str, object], audit: Mapping[str, object],
    ) -> Mapping[str, object]:
        """Replay all artifacts and ownership after audit publication."""
        if audit.get("path") != str(self.plan.audit_binding) or \
                not isinstance(audit.get("binding"), dict) or \
                not isinstance(audit.get("value"), dict):
            raise TransitionError("postpublication audit receipt is malformed")
        replay = self._replay_external_state(
            candidate_accept, archive, restored)
        audit_binding = file_binding(self.plan.audit_binding, with_hash=True)
        audit_value = load_canonical(
            self.plan.audit_binding, AUDIT_BINDING_SCHEMA,
            "postpublication future audit binding")
        if audit_binding != audit["binding"] or audit_value != audit["value"]:
            raise TransitionError("future audit binding changed after publication")
        if audit_value.get("prepublication_replay") != \
                audit["value"].get("prepublication_replay") or \
                audit_value.get("candidate_accept") != \
                (dict(candidate_accept)
                 if candidate_accept is not None else None):
            raise TransitionError("audit replay contract changed")
        return {**replay,
                "audit": {"binding": audit_binding,
                          "path": str(self.plan.audit_binding),
                          "value_sha256": audit_binding["sha256"]}}


def execute_transition(plan: TransitionPlan) -> Dict[str, object]:
    execute_runtime = verify_execute_runtime(plan)
    deadline = Deadline(plan.deadline_s, plan.recovery_reserve_s)
    value = verify_transition_plan(plan)
    if Path(__file__).resolve() != plan.controller.resolve():
        raise TransitionError("execute mode requires the frozen controller path")
    tools = value.get("tools")
    if not isinstance(tools, dict):
        raise TransitionError("verified transition plan lacks its tool ledger")
    verify_running_interpreter(tools, require_exact_path=True)
    p32 = load_verified_p32(plan.p32, plan.p32_sha256)
    journal = ReceiptJournal(
        plan.receipts, plan.transition_id, plan.controller_uid, deadline)
    backend = LiveBackend(plan, p32, deadline, tools, execute_runtime)
    lock_path = plan.root / "transition.lock"
    lock_fd = os.open(
        str(lock_path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
        getattr(os, "O_NOFOLLOW", 0),
        0o400)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.fsync(lock_fd)
        fsync_directory(plan.root)
        with DeferredSignals() as guard:
            controller = TransitionController(
                plan, backend, journal, deadline, guard)
            return controller.run()
    finally:
        os.close(lock_fd)


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--prepare-sealed-transition", action="store_true")
    modes.add_argument("--execute-sealed-transition")
    parser.add_argument("--expected-controller-sha256", required=True)
    parser.add_argument("--confirmation")
    args = parser.parse_args(argv)
    if SHA256_RE.fullmatch(args.expected_controller_sha256) is None:
        parser.error("expected controller SHA256 must be canonical lowercase hex")
    if args.prepare_sealed_transition:
        if args.confirmation is not None:
            parser.error("prepare mode does not accept a live confirmation")
    elif args.execute_sealed_transition != TRANSITION_ID or \
            args.confirmation != EXECUTE_CONFIRMATION:
        parser.error("execute mode requires both exact frozen confirmation tokens")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    plan = TransitionPlan(controller_sha256=args.expected_controller_sha256)
    if args.prepare_sealed_transition:
        result = prepare_transition(
            plan, controller_source=Path(__file__).resolve())
    else:
        result = execute_transition(plan)
    print(canonical_json(result).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
