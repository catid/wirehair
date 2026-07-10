#!/usr/bin/env python3
"""Failure-atomic, no-clobber publication for one or more local files."""

import os
import stat
from pathlib import Path

DIRECTORY_FSYNC_SUPPORTED = os.name != "nt"


class AtomicPublishError(ValueError):
    pass


def _identity(path):
    stat_result = os.lstat(os.fspath(path))
    return stat_result.st_dev, stat_result.st_ino, stat.S_IFMT(stat_result.st_mode)


def _regular_identity(stat_result):
    if not stat.S_ISREG(stat_result.st_mode):
        raise AtomicPublishError("staged output must be a regular file")
    return stat_result.st_dev, stat_result.st_ino, stat.S_IFMT(stat_result.st_mode)


def _link_no_follow(source, destination):
    try:
        os.link(source, destination, follow_symlinks=False)
    except (TypeError, NotImplementedError):
        # Older Windows Python builds do not expose follow_symlinks even though
        # CreateHardLink itself links the named file without replacing it.
        os.link(source, destination)


def _unlink_if_identity(path, identity):
    path = Path(path)
    try:
        if _identity(path) == identity:
            path.unlink()
    except FileNotFoundError:
        pass


def _lexists(path):
    return os.path.lexists(os.fspath(path))


def _sync_directories(parents):
    if not DIRECTORY_FSYNC_SUPPORTED:
        return
    for parent in parents:
        directory_fd = os.open(parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)


def _close_staged_handles(staged):
    errors = []
    for _, _, handle in staged:
        if handle.closed:
            continue
        try:
            handle.close()
        except OSError as exc:
            errors.append(exc)
    return errors


def publish_files_no_replace(entries):
    """Publish all ``(path, bytes-or-str)`` entries, or leave none published.

    Fixed sibling ``.tmp`` names are created exclusively.  Hard-link
    publication is atomic and fails when a destination appears concurrently.
    Any destination linked by this call is rolled back only when its inode still
    matches the staged file, so unrelated concurrent files are never removed.
    """
    normalized = []
    for path, content in entries:
        output = Path(path)
        temporary = Path(str(output) + ".tmp")
        data = content.encode("ascii") if isinstance(content, str) else bytes(content)
        normalized.append((output, temporary, data))
    if not normalized:
        raise AtomicPublishError("no outputs to publish")
    all_paths = [
        os.path.abspath(path)
        for output, temporary, _ in normalized
        for path in (output, temporary)
    ]
    if len(set(all_paths)) != len(all_paths):
        raise AtomicPublishError("output and temporary paths must be distinct")
    for output, temporary, _ in normalized:
        if _lexists(output) or _lexists(temporary):
            raise AtomicPublishError(f"output already exists: {output}")

    staged = []
    published = []
    try:
        for _, temporary, data in normalized:
            handle = temporary.open("xb")
            staged_stat = os.fstat(handle.fileno())
            staged.append((
                temporary,
                _regular_identity(staged_stat),
                handle,
            ))
            try:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            except Exception:
                handle.close()
                raise

        for (output, temporary, data), (_, staged_identity, _) in zip(
            normalized, staged
        ):
            if _identity(temporary) != staged_identity:
                raise AtomicPublishError(
                    f"staged file identity changed before publish: {temporary}"
                )
            _link_no_follow(temporary, output)
            # Record the expected identity before the post-link stat.  If the
            # source was swapped in the narrow pre-link window, replace it with
            # the actual identity before raising so rollback removes our link.
            published.append((output, staged_identity))
            temporary_identity = _identity(temporary)
            output_identity = _identity(output)
            if output_identity != staged_identity:
                if temporary_identity != staged_identity and \
                        output_identity == temporary_identity:
                    published[-1] = (output, output_identity)
                raise AtomicPublishError(
                    f"published inode does not match staged file: {output}"
                )
            if output.read_bytes() != data:
                raise AtomicPublishError(
                    f"published content does not match staged data: {output}"
                )
            if _identity(output) != staged_identity:
                raise AtomicPublishError(
                    f"published inode changed during verification: {output}"
                )

        close_errors = _close_staged_handles(staged)
        if close_errors:
            raise AtomicPublishError(
                f"failed closing {len(close_errors)} staged output file(s)"
            )
        for temporary, staged_identity, _ in staged:
            _unlink_if_identity(temporary, staged_identity)
        _sync_directories({output.parent for output, _, _ in normalized})
    except Exception as primary_error:
        cleanup_errors = _close_staged_handles(staged)
        for output, output_identity in reversed(published):
            try:
                _unlink_if_identity(output, output_identity)
            except OSError as exc:
                cleanup_errors.append(exc)
        try:
            _sync_directories({output.parent for output, _, _ in normalized})
        except OSError as exc:
            cleanup_errors.append(exc)
        if cleanup_errors:
            raise AtomicPublishError(
                f"{primary_error}; rollback had {len(cleanup_errors)} failure(s)"
            ) from primary_error
        raise
    finally:
        _close_staged_handles(staged)
        for temporary, staged_identity, _ in staged:
            try:
                _unlink_if_identity(temporary, staged_identity)
            except OSError:
                pass
