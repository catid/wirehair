#!/usr/bin/env python3
"""Fail-closed Linux filesystem primitives for WH2 timing evidence.

The timing controllers use this module only for small, immutable evidence
artifacts.  Reads are bounded and prove that the descriptor and pathname kept
the same identity for two complete passes.  Publications never replace an
existing name.  Directory publication deliberately has no portable fallback:
Linux ``renameat2(RENAME_NOREPLACE)`` is part of the evidence contract.

All public operations accept an ``error_type`` so a controller can surface its
own substantive-error class without translating exceptions at every call.
"""

from __future__ import annotations

from contextlib import contextmanager
import ctypes
import errno
import fcntl
from functools import wraps
import hashlib
import json
import math
import os
from pathlib import Path
import secrets
import stat
import threading
from typing import Any, Dict, Iterator, Optional, Tuple, Type


class EvidenceIOError(RuntimeError):
    """An evidence path, artifact, lock, or publication is unsafe."""


ErrorType = Type[Exception]
_READ_CHUNK = 1024 * 1024
_RENAME_NOREPLACE = 1


def _fail(error_type: ErrorType, message: str, cause: Optional[BaseException] = None) -> None:
    error = error_type(message)
    if cause is None:
        raise error
    raise error from cause


def _translate_os_errors(message: str) -> Any:
    """Keep ordinary public operations inside the controller's error domain."""
    def decorate(function: Any) -> Any:
        @wraps(function)
        def wrapped(*args: Any, **kwargs: Any) -> Any:
            error_type = kwargs.get("error_type", EvidenceIOError)
            try:
                return function(*args, **kwargs)
            except OSError as exc:
                _fail(error_type, message, exc)
        return wrapped
    return decorate


def _required_flag(name: str, error_type: ErrorType) -> int:
    value = getattr(os, name, None)
    if not isinstance(value, int):
        _fail(error_type, "Linux %s support is required" % name)
    return value


def _directory_flags(error_type: ErrorType) -> int:
    return (
        os.O_RDONLY
        | _required_flag("O_CLOEXEC", error_type)
        | _required_flag("O_DIRECTORY", error_type)
        | _required_flag("O_NOFOLLOW", error_type)
    )


def _file_read_flags(error_type: ErrorType) -> int:
    return (
        os.O_RDONLY
        | _required_flag("O_CLOEXEC", error_type)
        | _required_flag("O_NOFOLLOW", error_type)
        | _required_flag("O_NONBLOCK", error_type)
    )


def _absolute_path(path: Path, error_type: ErrorType) -> Path:
    try:
        raw = os.fspath(path)
    except TypeError as exc:
        _fail(error_type, "evidence path is not path-like", exc)
    if isinstance(raw, bytes):
        _fail(error_type, "evidence paths must be text paths")
    if not raw or "\x00" in raw:
        _fail(error_type, "evidence path is empty or contains NUL")
    return Path(os.path.abspath(raw))


def _directory_binding(value: os.stat_result) -> Tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_mode, value.st_uid, value.st_gid)


def _file_identity(value: os.stat_result) -> Tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _open_directory_path(path: Path, error_type: ErrorType) -> Tuple[int, Path]:
    """Open an absolute directory one non-symlink component at a time."""
    absolute = _absolute_path(path, error_type)
    flags = _directory_flags(error_type)
    try:
        descriptor = os.open("/", flags)
    except OSError as exc:
        _fail(error_type, "cannot open filesystem root", exc)
    try:
        for component in absolute.parts[1:]:
            next_descriptor = -1
            try:
                next_descriptor = os.open(component, flags, dir_fd=descriptor)
                opened = os.fstat(next_descriptor)
                named = os.stat(component, dir_fd=descriptor, follow_symlinks=False)
            except OSError as exc:
                if next_descriptor >= 0:
                    os.close(next_descriptor)
                _fail(error_type, "cannot open plain directory path %s" % absolute, exc)
            if (
                not stat.S_ISDIR(opened.st_mode)
                or _directory_binding(opened) != _directory_binding(named)
            ):
                os.close(next_descriptor)
                _fail(error_type, "directory component identity changed: %s" % absolute)
            previous = descriptor
            descriptor = next_descriptor
            os.close(previous)
        try:
            opened = os.fstat(descriptor)
            named = os.stat(str(absolute), follow_symlinks=False)
        except OSError as exc:
            _fail(error_type, "cannot re-prove directory path %s" % absolute, exc)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or _directory_binding(opened) != _directory_binding(named)
        ):
            _fail(error_type, "directory pathname identity changed: %s" % absolute)
        return descriptor, absolute
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def _verify_directory_binding(
    descriptor: int, absolute: Path, error_type: ErrorType
) -> os.stat_result:
    try:
        opened = os.fstat(descriptor)
        named = os.stat(str(absolute), follow_symlinks=False)
    except OSError as exc:
        _fail(error_type, "directory pathname became unavailable: %s" % absolute, exc)
    if (
        not stat.S_ISDIR(opened.st_mode)
        or _directory_binding(opened) != _directory_binding(named)
    ):
        _fail(error_type, "directory pathname identity changed: %s" % absolute)
    return opened


def _verify_private_directory_stat(
    metadata: os.stat_result,
    display: Path,
    error_type: ErrorType,
    require_writable: Optional[bool],
) -> None:
    mode = stat.S_IMODE(metadata.st_mode)
    if not stat.S_ISDIR(metadata.st_mode):
        _fail(error_type, "not a plain directory: %s" % display)
    if metadata.st_uid != os.geteuid():
        _fail(error_type, "directory is not owned by the effective user: %s" % display)
    if mode & 0o7077:
        _fail(error_type, "directory is not owner-only: %s" % display)
    if mode & (stat.S_IRUSR | stat.S_IXUSR) != (stat.S_IRUSR | stat.S_IXUSR):
        _fail(error_type, "owner cannot read/traverse directory: %s" % display)
    if require_writable is True and not mode & stat.S_IWUSR:
        _fail(error_type, "owner cannot write directory: %s" % display)
    if require_writable is False and mode & 0o222:
        _fail(error_type, "sealed directory is still writable: %s" % display)


@contextmanager
def held_private_directory(
    path: Path,
    *,
    require_writable: Optional[bool] = None,
    error_type: ErrorType = EvidenceIOError,
) -> Iterator[int]:
    """Hold a symlink-free, owner-only directory descriptor.

    ``require_writable=True`` is suitable for private staging/lock parents.
    ``False`` verifies a sealed owner-only directory (normally mode 0500).
    The pathname is re-proved when the context exits.
    """
    if require_writable is not None and not isinstance(require_writable, bool):
        _fail(error_type, "require_writable must be boolean or None")
    descriptor, absolute = _open_directory_path(path, error_type)
    try:
        metadata = _verify_directory_binding(descriptor, absolute, error_type)
        _verify_private_directory_stat(
            metadata, absolute, error_type, require_writable
        )
        try:
            yield descriptor
        finally:
            metadata = _verify_directory_binding(descriptor, absolute, error_type)
            _verify_private_directory_stat(
                metadata, absolute, error_type, require_writable
            )
    finally:
        os.close(descriptor)


def _validate_limit(max_bytes: int, error_type: ErrorType) -> None:
    if (
        not isinstance(max_bytes, int)
        or isinstance(max_bytes, bool)
        or max_bytes < 0
    ):
        _fail(error_type, "max_bytes must be a nonnegative integer")


def _validate_regular(
    metadata: os.stat_result,
    display: Path,
    require_unique: bool,
    max_bytes: int,
    error_type: ErrorType,
) -> None:
    if not stat.S_ISREG(metadata.st_mode):
        _fail(error_type, "refusing non-regular evidence input: %s" % display)
    if require_unique and metadata.st_nlink != 1:
        _fail(error_type, "refusing non-unique evidence input: %s" % display)
    if metadata.st_size < 0 or metadata.st_size > max_bytes:
        _fail(error_type, "evidence input exceeds byte limit: %s" % display)


def _read_bounded_pass(
    descriptor: int,
    max_bytes: int,
    collect: bool,
    display: Path,
    error_type: ErrorType,
) -> Tuple[Any, int]:
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
    except OSError as exc:
        _fail(error_type, "evidence input is not seekable: %s" % display, exc)
    total = 0
    digest = hashlib.sha256()
    chunks = []
    while True:
        allowance = max_bytes - total
        request = min(_READ_CHUNK, allowance + 1)
        try:
            block = os.read(descriptor, request)
        except OSError as exc:
            _fail(error_type, "cannot read evidence input: %s" % display, exc)
        if not block:
            break
        total += len(block)
        if total > max_bytes:
            _fail(error_type, "evidence input exceeds byte limit: %s" % display)
        digest.update(block)
        if collect:
            chunks.append(block)
    return (b"".join(chunks) if collect else digest.digest()), total


def _stable_regular_read(
    path: Path,
    *,
    max_bytes: int,
    require_unique: bool,
    collect: bool,
    error_type: ErrorType,
) -> Any:
    _validate_limit(max_bytes, error_type)
    if not isinstance(require_unique, bool):
        _fail(error_type, "require_unique must be boolean")
    absolute = _absolute_path(path, error_type)
    if absolute.name in ("", ".", ".."):
        _fail(error_type, "evidence input has no filename: %s" % absolute)
    parent_fd, parent = _open_directory_path(absolute.parent, error_type)
    descriptor = -1
    try:
        try:
            descriptor = os.open(
                absolute.name, _file_read_flags(error_type), dir_fd=parent_fd
            )
        except OSError as exc:
            _fail(error_type, "cannot open plain evidence input: %s" % absolute, exc)

        before = os.fstat(descriptor)
        try:
            named_before = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
        except OSError as exc:
            _fail(error_type, "evidence pathname changed while opening: %s" % absolute, exc)
        _validate_regular(before, absolute, require_unique, max_bytes, error_type)
        if _file_identity(before) != _file_identity(named_before):
            _fail(error_type, "evidence descriptor/path identity mismatch: %s" % absolute)

        first, first_size = _read_bounded_pass(
            descriptor, max_bytes, collect, absolute, error_type
        )
        midpoint = os.fstat(descriptor)
        try:
            named_midpoint = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
        except OSError as exc:
            _fail(error_type, "evidence pathname changed while reading: %s" % absolute, exc)
        _validate_regular(midpoint, absolute, require_unique, max_bytes, error_type)

        second, second_size = _read_bounded_pass(
            descriptor, max_bytes, collect, absolute, error_type
        )
        after = os.fstat(descriptor)
        try:
            named_after = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
        except OSError as exc:
            _fail(error_type, "evidence pathname changed after reading: %s" % absolute, exc)
        _validate_regular(after, absolute, require_unique, max_bytes, error_type)

        identities = (
            _file_identity(before),
            _file_identity(named_before),
            _file_identity(midpoint),
            _file_identity(named_midpoint),
            _file_identity(after),
            _file_identity(named_after),
        )
        if any(identity != identities[0] for identity in identities[1:]):
            _fail(error_type, "evidence input changed while reading: %s" % absolute)
        if first_size != before.st_size or second_size != before.st_size:
            _fail(error_type, "evidence size changed while reading: %s" % absolute)
        if first != second:
            _fail(error_type, "evidence contents changed while reading: %s" % absolute)
        _verify_directory_binding(parent_fd, parent, error_type)
        final_descriptor = os.fstat(descriptor)
        try:
            named_final = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
        except OSError as exc:
            _fail(error_type, "evidence pathname changed at confirmation: %s" % absolute, exc)
        if (
            _file_identity(final_descriptor) != identities[0]
            or _file_identity(named_final) != identities[0]
        ):
            _fail(error_type, "evidence input changed at confirmation: %s" % absolute)
        return first
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)


@_translate_os_errors("stable evidence snapshot failed")
def stable_file_snapshot(
    path: Path,
    *,
    max_bytes: int,
    require_unique: bool = True,
    error_type: ErrorType = EvidenceIOError,
) -> bytes:
    """Return a bounded, twice-confirmed snapshot of one regular file."""
    return _stable_regular_read(
        path,
        max_bytes=max_bytes,
        require_unique=require_unique,
        collect=True,
        error_type=error_type,
    )


@_translate_os_errors("stable evidence SHA256 failed")
def stable_file_sha256(
    path: Path,
    *,
    max_bytes: int,
    require_unique: bool = True,
    error_type: ErrorType = EvidenceIOError,
) -> str:
    """Stream and independently confirm a bounded regular-file SHA256."""
    digest = _stable_regular_read(
        path,
        max_bytes=max_bytes,
        require_unique=require_unique,
        collect=False,
        error_type=error_type,
    )
    return digest.hex()


class _StrictJSONError(ValueError):
    pass


def _strict_object(pairs: Any) -> Dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise _StrictJSONError("duplicate JSON key %r" % key)
        result[key] = value
    return result


def _strict_constant(token: str) -> Any:
    raise _StrictJSONError("non-finite JSON constant %s" % token)


def _strict_float(token: str) -> float:
    value = float(token)
    if not math.isfinite(value):
        raise _StrictJSONError("non-finite JSON number %s" % token)
    return value


def canonical_json_bytes(
    value: object, *, error_type: ErrorType = EvidenceIOError
) -> bytes:
    """Encode the compact ASCII canonical-object format used by controllers."""
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("ascii")
    except (
        TypeError,
        ValueError,
        UnicodeError,
        RecursionError,
        OverflowError,
    ) as exc:
        _fail(error_type, "value cannot be encoded as canonical JSON", exc)


def load_canonical_object_snapshot(
    path: Path,
    *,
    max_bytes: int,
    name: str = "JSON evidence",
    require_unique: bool = True,
    error_type: ErrorType = EvidenceIOError,
) -> Tuple[Dict[str, Any], bytes]:
    """Load one strict canonical JSON object and return its exact bytes."""
    raw = stable_file_snapshot(
        path,
        max_bytes=max_bytes,
        require_unique=require_unique,
        error_type=error_type,
    )
    try:
        text = raw.decode("ascii")
        value = json.loads(
            text,
            object_pairs_hook=_strict_object,
            parse_constant=_strict_constant,
            parse_float=_strict_float,
        )
    except (
        UnicodeDecodeError,
        ValueError,
        RecursionError,
        OverflowError,
    ) as exc:
        _fail(error_type, "%s is not strict ASCII JSON" % name, exc)
    if not isinstance(value, dict) or canonical_json_bytes(
        value, error_type=error_type
    ) != raw:
        _fail(error_type, "%s is not an exact canonical JSON object" % name)
    return value, raw


def load_canonical_object(
    path: Path,
    *,
    max_bytes: int,
    name: str = "JSON evidence",
    require_unique: bool = True,
    error_type: ErrorType = EvidenceIOError,
) -> Dict[str, Any]:
    """Load one strict canonical JSON object."""
    value, _raw = load_canonical_object_snapshot(
        path,
        max_bytes=max_bytes,
        name=name,
        require_unique=require_unique,
        error_type=error_type,
    )
    return value


def _validate_publication_path(path: Path, error_type: ErrorType) -> Path:
    absolute = _absolute_path(path, error_type)
    if absolute.name in ("", ".", ".."):
        _fail(error_type, "publication target has no safe filename: %s" % absolute)
    return absolute


@_translate_os_errors("immutable evidence publication failed")
def publish_immutable_file(
    path: Path,
    data: bytes,
    *,
    mode: int = 0o400,
    error_type: ErrorType = EvidenceIOError,
) -> None:
    """Durably hard-link immutable bytes into a previously absent name.

    An existing name is always an error, even when its bytes match.  This
    strict no-clobber rule makes concurrent-publisher mistakes visible.
    """
    if not isinstance(data, bytes):
        _fail(error_type, "immutable publication data must be bytes")
    if (
        not isinstance(mode, int)
        or isinstance(mode, bool)
        or mode < 0
        or mode > 0o777
        or mode & 0o222
        or not mode & stat.S_IRUSR
    ):
        _fail(
            error_type,
            "immutable publication mode must be owner-readable with no write bits",
        )
    absolute = _validate_publication_path(path, error_type)
    temporary_name = ".%s.%d.%d.%s.partial" % (
        absolute.name,
        os.getpid(),
        threading.get_ident(),
        secrets.token_hex(8),
    )
    temporary_exists = False
    with held_private_directory(
        absolute.parent, require_writable=True, error_type=error_type
    ) as parent_fd:
        descriptor = -1
        try:
            flags = (
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | _required_flag("O_CLOEXEC", error_type)
                | _required_flag("O_NOFOLLOW", error_type)
            )
            try:
                descriptor = os.open(
                    temporary_name, flags, 0o600, dir_fd=parent_fd
                )
            except OSError as exc:
                _fail(error_type, "cannot create immutable publication temporary", exc)
            temporary_exists = True
            view = memoryview(data)
            while view:
                try:
                    written = os.write(descriptor, view)
                except OSError as exc:
                    _fail(error_type, "cannot write immutable publication temporary", exc)
                if written <= 0:
                    _fail(error_type, "short write to immutable publication temporary")
                view = view[written:]
            os.fchmod(descriptor, mode)
            os.fsync(descriptor)
            temporary_stat = os.fstat(descriptor)
            if (
                not stat.S_ISREG(temporary_stat.st_mode)
                or temporary_stat.st_nlink != 1
                or temporary_stat.st_size != len(data)
                or stat.S_IMODE(temporary_stat.st_mode) != mode
            ):
                _fail(error_type, "immutable publication temporary changed")
            try:
                os.link(
                    temporary_name,
                    absolute.name,
                    src_dir_fd=parent_fd,
                    dst_dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileExistsError as exc:
                _fail(error_type, "immutable publication target already exists: %s" % absolute, exc)
            except OSError as exc:
                _fail(error_type, "immutable no-clobber link failed: %s" % absolute, exc)
            linked = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
            marker = os.stat(
                temporary_name, dir_fd=parent_fd, follow_symlinks=False
            )
            opened_linked = os.fstat(descriptor)
            if (
                not stat.S_ISREG(linked.st_mode)
                or linked.st_nlink != 2
                or _file_identity(linked) != _file_identity(marker)
                or _file_identity(linked) != _file_identity(opened_linked)
            ):
                _fail(error_type, "immutable publication link identity changed")
            os.fsync(parent_fd)
            os.unlink(temporary_name, dir_fd=parent_fd)
            temporary_exists = False
            os.fsync(parent_fd)
            published = os.stat(
                absolute.name, dir_fd=parent_fd, follow_symlinks=False
            )
            opened_published = os.fstat(descriptor)
            if (
                not stat.S_ISREG(published.st_mode)
                or published.st_nlink != 1
                or published.st_size != len(data)
                or (published.st_dev, published.st_ino)
                != (linked.st_dev, linked.st_ino)
                or stat.S_IMODE(published.st_mode) != mode
                or _file_identity(published) != _file_identity(opened_published)
            ):
                _fail(error_type, "immutable publication was not stable")
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            if temporary_exists:
                try:
                    os.unlink(temporary_name, dir_fd=parent_fd)
                    os.fsync(parent_fd)
                except FileNotFoundError:
                    pass


def _renameat2_function(error_type: ErrorType) -> Any:
    try:
        function = ctypes.CDLL(None, use_errno=True).renameat2
    except (AttributeError, OSError) as exc:
        _fail(error_type, "Linux renameat2(RENAME_NOREPLACE) is required", exc)
    function.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    function.restype = ctypes.c_int
    return function


@_translate_os_errors("no-replace directory publication failed")
def publish_directory_noreplace(
    source: Path,
    target: Path,
    *,
    error_type: ErrorType = EvidenceIOError,
) -> None:
    """Publish one sealed sibling directory with fail-closed renameat2."""
    source_absolute = _validate_publication_path(source, error_type)
    target_absolute = _validate_publication_path(target, error_type)
    if source_absolute == target_absolute or source_absolute.parent != target_absolute.parent:
        _fail(error_type, "directory publication requires distinct sibling paths")
    with held_private_directory(
        source_absolute.parent, require_writable=True, error_type=error_type
    ) as parent_fd:
        source_fd = -1
        try:
            try:
                source_fd = os.open(
                    source_absolute.name,
                    _directory_flags(error_type),
                    dir_fd=parent_fd,
                )
                before = os.fstat(source_fd)
                named_before = os.stat(
                    source_absolute.name,
                    dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except OSError as exc:
                _fail(error_type, "cannot open sealed publication directory", exc)
            if _directory_binding(before) != _directory_binding(named_before):
                _fail(error_type, "sealed publication directory identity changed")
            _verify_private_directory_stat(
                before, source_absolute, error_type, require_writable=False
            )
            os.fsync(source_fd)
            renameat2 = _renameat2_function(error_type)
            result = renameat2(
                parent_fd,
                os.fsencode(source_absolute.name),
                parent_fd,
                os.fsencode(target_absolute.name),
                _RENAME_NOREPLACE,
            )
            if result != 0:
                error_number = ctypes.get_errno()
                if error_number == errno.EEXIST:
                    _fail(
                        error_type,
                        "directory publication target already exists: %s"
                        % target_absolute,
                    )
                _fail(
                    error_type,
                    "renameat2(RENAME_NOREPLACE) failed: %s"
                    % os.strerror(error_number),
                )
            try:
                os.stat(
                    source_absolute.name,
                    dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                _fail(error_type, "directory source name remained after publication")
            try:
                named_after = os.stat(
                    target_absolute.name,
                    dir_fd=parent_fd,
                    follow_symlinks=False,
                )
                opened_after = os.fstat(source_fd)
            except OSError as exc:
                _fail(error_type, "published directory identity is unavailable", exc)
            if (
                _directory_binding(named_after) != _directory_binding(opened_after)
                or (opened_after.st_dev, opened_after.st_ino)
                != (before.st_dev, before.st_ino)
            ):
                _fail(error_type, "published directory identity changed")
            _verify_private_directory_stat(
                opened_after, target_absolute, error_type, require_writable=False
            )
            os.fsync(source_fd)
            os.fsync(parent_fd)
            durable = os.stat(
                target_absolute.name,
                dir_fd=parent_fd,
                follow_symlinks=False,
            )
            durable_opened = os.fstat(source_fd)
            if (
                _directory_binding(durable) != _directory_binding(durable_opened)
                or (durable_opened.st_dev, durable_opened.st_ino)
                != (before.st_dev, before.st_ino)
            ):
                _fail(error_type, "published directory changed after durability flush")
            _verify_private_directory_stat(
                durable_opened,
                target_absolute,
                error_type,
                require_writable=False,
            )
            try:
                os.stat(
                    source_absolute.name,
                    dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                _fail(error_type, "directory source name reappeared after publication")
        finally:
            if source_fd >= 0:
                os.close(source_fd)


@contextmanager
def nonblocking_global_flock(
    path: Path, *, error_type: ErrorType = EvidenceIOError
) -> Iterator[int]:
    """Acquire one private regular-file flock, failing instead of waiting."""
    absolute = _validate_publication_path(path, error_type)
    with held_private_directory(
        absolute.parent, require_writable=True, error_type=error_type
    ) as parent_fd:
        descriptor = -1
        created = False
        try:
            common_flags = (
                os.O_RDWR
                | _required_flag("O_CLOEXEC", error_type)
                | _required_flag("O_NOFOLLOW", error_type)
                | _required_flag("O_NONBLOCK", error_type)
            )
            try:
                descriptor = os.open(
                    absolute.name,
                    common_flags | os.O_CREAT | os.O_EXCL,
                    0o600,
                    dir_fd=parent_fd,
                )
                created = True
            except FileExistsError:
                try:
                    descriptor = os.open(
                        absolute.name, common_flags, dir_fd=parent_fd
                    )
                except OSError as exc:
                    _fail(error_type, "cannot open global evidence lock", exc)
            except OSError as exc:
                _fail(error_type, "cannot create global evidence lock", exc)
            try:
                if created:
                    os.fchmod(descriptor, 0o600)
                    os.fsync(descriptor)
                    os.fsync(parent_fd)
                opened = os.fstat(descriptor)
                named = os.stat(
                    absolute.name, dir_fd=parent_fd, follow_symlinks=False
                )
            except OSError as exc:
                _fail(error_type, "global evidence lock identity is unavailable", exc)
            if (
                not stat.S_ISREG(opened.st_mode)
                or opened.st_nlink != 1
                or opened.st_uid != os.geteuid()
                or stat.S_IMODE(opened.st_mode) != 0o600
                or _file_identity(opened) != _file_identity(named)
            ):
                _fail(error_type, "global evidence lock is not a private unique file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except (BlockingIOError, OSError) as exc:
                _fail(error_type, "global evidence lock is already held", exc)
            try:
                locked_named = os.stat(
                    absolute.name, dir_fd=parent_fd, follow_symlinks=False
                )
                locked_opened = os.fstat(descriptor)
            except OSError as exc:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
                _fail(error_type, "global evidence lock identity is unavailable", exc)
            if _file_identity(locked_named) != _file_identity(locked_opened):
                fcntl.flock(descriptor, fcntl.LOCK_UN)
                _fail(error_type, "global evidence lock pathname changed")
            try:
                try:
                    yield descriptor
                finally:
                    try:
                        confirmed = os.stat(
                            absolute.name,
                            dir_fd=parent_fd,
                            follow_symlinks=False,
                        )
                        current = os.fstat(descriptor)
                    except OSError as exc:
                        _fail(
                            error_type,
                            "global evidence lock pathname became unavailable",
                            exc,
                        )
                    if _file_identity(confirmed) != _file_identity(current):
                        _fail(error_type, "global evidence lock pathname changed")
            finally:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            if descriptor >= 0:
                os.close(descriptor)


__all__ = (
    "EvidenceIOError",
    "canonical_json_bytes",
    "held_private_directory",
    "load_canonical_object",
    "load_canonical_object_snapshot",
    "nonblocking_global_flock",
    "publish_directory_noreplace",
    "publish_immutable_file",
    "stable_file_sha256",
    "stable_file_snapshot",
)
