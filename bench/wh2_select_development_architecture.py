#!/usr/bin/env python3
"""Select WH2 architecture only from completed development artifact dirs.

The selector has no summary-only interface.  It strictly reopens the native
four-arm recovery campaign, independently reopens and joins the raw work/rank
sidecar, strictly reopens the timing-only campaign, and publishes a receipt
only when the closed recovery and timing gates produce a winner.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import stat
import sys
import tempfile
from typing import Any, Mapping, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_precodefail_work_screen as work_api
import wh2_run_native_recovery_screen as recovery_api
import wh2_run_native_short_screen as timing_api


MAX_SELECTION_RECEIPT_BYTES = 1024 * 1024
ReceiptIdentity = Tuple[int, int]
ReceiptFingerprint = Tuple[int, int, int, Optional[int], Optional[int]]


class SelectionError(RuntimeError):
    """The authoritative development selection cannot be produced."""


def fail(message: str) -> None:
    raise SelectionError(message)


def _close_fd_quietly(descriptor: int) -> None:
    if descriptor < 0:
        return
    try:
        os.close(descriptor)
    except OSError:
        pass


def _receipt_fingerprint(value: os.stat_result) -> ReceiptFingerprint:
    return (
        value.st_dev, value.st_ino, value.st_size,
        getattr(value, "st_mtime_ns", None),
        getattr(value, "st_ctime_ns", None),
    )


def _receipt_basename(path: Path) -> str:
    name = Path(path).name
    if name in ("", ".", ".."):
        fail("selection receipt path has no regular-file name")
    return name


def _write_private(path: Path, data: bytes) -> None:
    try:
        with path.open("xb") as output:
            output.write(data)
            output.flush()
    except OSError as exc:
        fail("cannot write private selection snapshot {}: {}".format(
            path, exc))


def _verify_real_parent(path: Path, identity: ReceiptIdentity) -> None:
    try:
        current = os.stat(str(path), follow_symlinks=False)
    except OSError as exc:
        fail("selection receipt parent became unavailable: {}".format(exc))
    if (not stat.S_ISDIR(current.st_mode) or
            (current.st_dev, current.st_ino) != identity):
        fail("selection receipt parent identity changed")


def _open_real_parent(path: Path) -> Tuple[int, ReceiptIdentity]:
    parent = path.parent
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    directory = getattr(os, "O_DIRECTORY", 0)
    if nofollow == 0 or directory == 0:
        fail("selection receipt access requires O_NOFOLLOW and O_DIRECTORY")
    descriptor = -1
    try:
        descriptor = os.open(
            str(parent), os.O_RDONLY | directory |
            getattr(os, "O_CLOEXEC", 0) | nofollow)
        opened = os.fstat(descriptor)
        if not stat.S_ISDIR(opened.st_mode):
            fail("selection receipt parent must be a real directory")
        identity = (opened.st_dev, opened.st_ino)
        _verify_real_parent(parent, identity)
        return descriptor, identity
    except OSError as exc:
        _close_fd_quietly(descriptor)
        fail("cannot open selection receipt parent {}: {}".format(
            parent, exc))
    except BaseException:
        _close_fd_quietly(descriptor)
        raise
    raise AssertionError("unreachable")


def _require_parent_descriptor(
        descriptor: int, identity: ReceiptIdentity) -> None:
    """Require one still-open descriptor for the originally named parent."""
    try:
        opened = os.fstat(descriptor)
    except OSError as exc:
        fail("selection receipt parent descriptor is invalid: {}".format(exc))
    if (not stat.S_ISDIR(opened.st_mode) or
            (opened.st_dev, opened.st_ino) != identity):
        fail("selection receipt parent descriptor identity changed")


@dataclass(frozen=True)
class _PreparedReceipt:
    """One fully written, still-unnamed receipt held by an exact inode fd."""

    device: int
    inode: int
    descriptor: int
    byte_count: int
    sha256: str
    canonical_bytes: bytes

    @property
    def identity(self) -> ReceiptIdentity:
        return self.device, self.inode


class _ReceiptPublicationTransaction:
    """Pinned parent and unnamed inode for one publish-wins transaction."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.name = _receipt_basename(self.path)
        self.parent_fd = -1
        self.parent_identity: Optional[ReceiptIdentity] = None
        self.prepared: Optional[_PreparedReceipt] = None

    def open_parent(self) -> None:
        if self.parent_fd >= 0 or self.parent_identity is not None:
            fail("selection receipt publication transaction is already open")
        parent = self.path.parent
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        directory = getattr(os, "O_DIRECTORY", 0)
        if nofollow == 0 or directory == 0:
            fail("selection receipt access requires O_NOFOLLOW and O_DIRECTORY")
        descriptor = -1
        failed = True
        try:
            descriptor = os.open(
                str(parent), os.O_RDONLY | directory |
                getattr(os, "O_CLOEXEC", 0) | nofollow)
            # Transfer ownership before any helper call.  The finally guard
            # recognizes the brief local/self alias interval, so an injected
            # BaseException cannot make both aliases close the same reused fd.
            self.parent_fd = descriptor
            descriptor = -1
            opened = os.fstat(self.parent_fd)
            if not stat.S_ISDIR(opened.st_mode):
                fail("selection receipt parent must be a real directory")
            self.parent_identity = (opened.st_dev, opened.st_ino)
            _verify_real_parent(parent, self.parent_identity)
            failed = False
        except OSError as exc:
            fail("cannot open selection receipt parent {}: {}".format(
                parent, exc))
        finally:
            if descriptor >= 0:
                if self.parent_fd == descriptor:
                    # self is the sole owner; clear only the local alias.
                    descriptor = -1
                else:
                    _close_fd_quietly(descriptor)
                    descriptor = -1
            if failed:
                self.close()

    def require_parent(self) -> Tuple[int, ReceiptIdentity]:
        if self.parent_fd < 0 or self.parent_identity is None:
            fail("selection receipt publication transaction is not open")
        _require_parent_descriptor(self.parent_fd, self.parent_identity)
        return self.parent_fd, self.parent_identity

    def require_prepared(self) -> _PreparedReceipt:
        if self.prepared is None:
            fail("selection receipt publication transaction is not prepared")
        try:
            retained = os.fstat(self.prepared.descriptor)
        except OSError as exc:
            fail("prepared selection receipt descriptor is invalid: {}".format(
                exc))
        if (not stat.S_ISREG(retained.st_mode) or
                (retained.st_dev, retained.st_ino) !=
                self.prepared.identity or
                retained.st_size != self.prepared.byte_count or
                retained.st_nlink != 0):
            fail("prepared selection receipt descriptor changed")
        return self.prepared

    def verify_visible_parent(self) -> None:
        _parent_fd, parent_identity = self.require_parent()
        _verify_real_parent(self.path.parent, parent_identity)

    def close(self) -> None:
        prepared = self.prepared
        self.prepared = None
        descriptor = self.parent_fd
        self.parent_fd = -1
        self.parent_identity = None
        first_async_error: Optional[BaseException] = None
        for owned in (
                prepared.descriptor if prepared is not None else -1,
                descriptor):
            if owned < 0:
                continue
            try:
                os.close(owned)
            except OSError:
                pass
            except BaseException as exc:
                if first_async_error is None:
                    first_async_error = exc
        if first_async_error is not None:
            raise first_async_error


def _read_unnamed_receipt(
        descriptor: int,
        expected_identity: Optional[ReceiptIdentity] = None,
        ) -> Tuple[bytes, ReceiptFingerprint]:
    """Read and fingerprint one exact, regular, still-unnamed descriptor."""
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 0:
            fail("prepared architecture selection is not an unnamed file")
        if (expected_identity is not None and
                (before.st_dev, before.st_ino) != expected_identity):
            fail("prepared architecture selection inode identity changed")
        if before.st_size <= 0 or before.st_size > MAX_SELECTION_RECEIPT_BYTES:
            fail("prepared architecture selection has an invalid byte size")
        offset = 0
        chunks = []
        while offset < before.st_size:
            block = os.pread(
                descriptor, min(64 * 1024, before.st_size - offset), offset)
            if not block:
                fail("prepared architecture selection was truncated")
            chunks.append(block)
            offset += len(block)
        if os.pread(descriptor, 1, before.st_size):
            fail("prepared architecture selection grew while reading")
        after = os.fstat(descriptor)
        if (after.st_nlink != 0 or
                _receipt_fingerprint(before) != _receipt_fingerprint(after)):
            fail("prepared architecture selection changed while reading")
        return b"".join(chunks), _receipt_fingerprint(after)
    except OSError as exc:
        fail("cannot inspect unnamed architecture selection: {}".format(exc))
    raise AssertionError("unreachable")


def _open_unnamed_receipt(transaction: _ReceiptPublicationTransaction) -> int:
    """Allocate an O_TMPFILE inode without ever creating a directory entry."""
    parent_fd, _parent_identity = transaction.require_parent()
    if transaction.prepared is not None:
        fail("selection receipt publication transaction was already used")
    tmpfile = getattr(os, "O_TMPFILE", 0)
    cloexec = getattr(os, "O_CLOEXEC", 0)
    if tmpfile == 0 or cloexec == 0:
        fail("selection receipt publication requires O_TMPFILE and O_CLOEXEC")
    try:
        return os.open(
            ".", tmpfile | os.O_RDWR | cloexec,
            0o600, dir_fd=parent_fd)
    except OSError as exc:
        fail("selection receipt parent does not support O_TMPFILE: {}".format(
            exc))
    raise AssertionError("unreachable")


def _read_stable_receipt_bytes_at(
        path: Path, parent_fd: int, parent_identity: ReceiptIdentity,
        expected_identity: Optional[ReceiptIdentity] = None,
        ) -> Tuple[bytes, ReceiptFingerprint]:
    """Bounded no-follow read through one already-pinned real parent."""
    path = Path(path)
    name = _receipt_basename(path)
    _require_parent_descriptor(parent_fd, parent_identity)
    descriptor = -1
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0) |
            getattr(os, "O_NONBLOCK", 0),
            dir_fd=parent_fd)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("selection receipt must be a regular file")
        if (expected_identity is not None and
                (before.st_dev, before.st_ino) != expected_identity):
            fail("selection receipt is not the published transaction inode")
        if before.st_size <= 0 or before.st_size > MAX_SELECTION_RECEIPT_BYTES:
            fail("selection receipt is empty or exceeds its byte cap")
        remaining = before.st_size
        chunks = []
        while remaining:
            block = os.read(descriptor, min(64 * 1024, remaining))
            if not block:
                fail("selection receipt was truncated while reading")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            fail("selection receipt grew beyond its opened size")
        after = os.fstat(descriptor)
        stable_before = _receipt_fingerprint(before)
        stable_after = _receipt_fingerprint(after)
        entry = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (stable_before != stable_after or
                _receipt_fingerprint(entry) != stable_after):
            fail("selection receipt changed while it was being read")
        data = b"".join(chunks)
        _verify_real_parent(path.parent, parent_identity)
        return data, stable_after
    except OSError as exc:
        fail("cannot read selection receipt {}: {}".format(path, exc))
    finally:
        _close_fd_quietly(descriptor)
    raise AssertionError("unreachable")


def _read_stable_receipt_bytes(
        path: Path,
        expected_identity: Optional[ReceiptIdentity] = None,
        ) -> Tuple[bytes, ReceiptFingerprint]:
    """Bounded, no-follow read pinned to one real parent and regular inode."""
    path = Path(path)
    parent_fd, parent_identity = _open_real_parent(path)
    try:
        return _read_stable_receipt_bytes_at(
            path, parent_fd, parent_identity, expected_identity)
    finally:
        _close_fd_quietly(parent_fd)
    raise AssertionError("unreachable")


def _prepare_receipt(
        transaction: _ReceiptPublicationTransaction,
        contract: Mapping[str, Any], value: Mapping[str, Any],
        ) -> contract_api.ArchitectureSelection:
    """Write, reread, validate, and seal a receipt while it is still unnamed."""
    if not isinstance(transaction, _ReceiptPublicationTransaction):
        fail("selection receipt publication requires transaction state")
    transaction.require_parent()
    data = (contract_api.canonical_json(value) + "\n").encode("utf-8")
    if not data or len(data) > MAX_SELECTION_RECEIPT_BYTES:
        fail("architecture selection receipt exceeds its publication cap")
    descriptor = _open_unnamed_receipt(transaction)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                fail("short write while preparing architecture selection")
            offset += written
        os.fsync(descriptor)
        reopened, fingerprint = _read_unnamed_receipt(descriptor)
        if reopened != data:
            fail("prepared architecture selection differs from written bytes")
        transaction.prepared = _PreparedReceipt(
            device=fingerprint[0], inode=fingerprint[1],
            descriptor=descriptor,
            byte_count=len(data),
            sha256=hashlib.sha256(data).hexdigest(), canonical_bytes=data)
        descriptor = -1
        parsed = contract_api._load_json_bytes(
            reopened, "architecture selection receipt")
        if not isinstance(parsed, dict):
            fail("architecture selection receipt must be a JSON object")
        logical = contract_api.canonical_json(parsed).encode("utf-8")
        if reopened != logical + b"\n":
            fail("architecture selection receipt must be canonical JSON")
        validated = contract_api.validate_selection_receipt(contract, parsed)
        if contract_api.canonical_json(validated) != \
                contract_api.canonical_json(value):
            fail("validated architecture selection differs from candidate")
        selection = contract_api._seal_validated_architecture_selection(
            contract, validated)
        if (type(selection) is not contract_api.ArchitectureSelection or
                selection.contract_sha256 !=
                contract_api.contract_sha256(contract) or
                selection.selection_sha256 != validated["selection_sha256"] or
                selection.canonical_receipt !=
                contract_api.canonical_json(validated)):
            fail("architecture selection sealer returned a malformed handle")
        retained = transaction.require_prepared()
        terminal_bytes, terminal_fingerprint = _read_unnamed_receipt(
            retained.descriptor, expected_identity=retained.identity)
        if (terminal_bytes != retained.canonical_bytes or
                terminal_fingerprint[2] != retained.byte_count or
                hashlib.sha256(terminal_bytes).hexdigest() != retained.sha256):
            fail("prepared architecture selection changed before publication")
        return selection
    except contract_api.ContractError as exc:
        fail(str(exc))
    except OSError as exc:
        fail("cannot prepare architecture selection: {}".format(exc))
    finally:
        if descriptor >= 0:
            prepared = transaction.prepared
            if prepared is not None and prepared.descriptor == descriptor:
                # The transaction owns this exact fd.  Do not close through
                # the stale local alias; its finally will close it exactly once.
                descriptor = -1
            else:
                _close_fd_quietly(descriptor)
                descriptor = -1
    raise AssertionError("unreachable")


def _publish_receipt(transaction: _ReceiptPublicationTransaction) -> None:
    """Link the sealed unnamed inode; a successful link always wins.

    The link is the sole public-name transition and the last operation.  An
    asynchronous exception after the kernel installs the link can lose the
    caller's acknowledgement, but must never trigger removal of a valid
    receipt.
    """
    parent_fd, _parent_identity = transaction.require_parent()
    prepared = transaction.require_prepared()
    retained_bytes, retained_fingerprint = _read_unnamed_receipt(
        prepared.descriptor, expected_identity=prepared.identity)
    if (retained_bytes != prepared.canonical_bytes or
            retained_fingerprint[2] != prepared.byte_count or
            hashlib.sha256(retained_bytes).hexdigest() != prepared.sha256):
        fail("prepared architecture selection changed before final link")
    # This visible-path check is deliberately the final operation before the
    # one-way link transition.
    transaction.verify_visible_parent()
    try:
        return os.link(
            "/proc/self/fd/{}".format(prepared.descriptor), transaction.name,
            dst_dir_fd=parent_fd, follow_symlinks=True)
    except FileExistsError:
        fail("refusing to replace existing selection receipt {}".format(
            transaction.path))
    except OSError as exc:
        fail("cannot link unnamed selection receipt {}: {}".format(
            transaction.path, exc))


def _require_receipt_snapshot(
        path: Path, fingerprint: ReceiptFingerprint,
        expected_bytes: bytes) -> None:
    """Terminally reread the exact stable receipt observed before sealing."""
    expected_identity = (fingerprint[0], fingerprint[1])
    current_bytes, current_fingerprint = _read_stable_receipt_bytes(
        Path(path), expected_identity=expected_identity)
    if (current_fingerprint != fingerprint or
            current_bytes != expected_bytes):
        fail("published architecture selection changed before commit")


def _visible_receipt_identity(path: Path) -> ReceiptIdentity:
    """Snapshot one no-follow regular receipt directory-entry identity."""
    path = Path(path)
    name = _receipt_basename(path)
    parent_fd, parent_identity = _open_real_parent(path)
    try:
        visible = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if not stat.S_ISREG(visible.st_mode):
            fail("selection receipt must be a regular file")
        _verify_real_parent(path.parent, parent_identity)
        return visible.st_dev, visible.st_ino
    except OSError as exc:
        fail("cannot inspect selection receipt {}: {}".format(path, exc))
    finally:
        _close_fd_quietly(parent_fd)
    raise AssertionError("unreachable")


def _revalidate_recovery_summary(
        contract: Mapping[str, Any], campaign: Mapping[str, Any],
        ) -> Mapping[str, Any]:
    """Recompute the ledger summary from the strict campaign snapshot."""
    freeze = campaign.get("freeze")
    rows = campaign.get("rows")
    trace_bytes = campaign.get("trace_bytes")
    receipt = campaign.get("receipt")
    if (not isinstance(freeze, Mapping) or
            not isinstance(rows, list) or
            not isinstance(trace_bytes, bytes) or
            not isinstance(receipt, Mapping)):
        fail("strict recovery loader returned a malformed evidence bundle")
    logical_freeze = {
        key: value for key, value in freeze.items()
        if key != "arms_by_name"
    }
    result_bytes = "".join(
        contract_api.canonical_json(row) + "\n" for row in rows
    ).encode("utf-8")
    with tempfile.TemporaryDirectory(prefix="wh2-selection-recovery-") as raw:
        root = Path(raw)
        freeze_path = root / "recovery-freeze.json"
        trace_path = root / "recovery-traces.jsonl"
        result_path = root / "recovery-results.jsonl"
        _write_private(
            freeze_path,
            (contract_api.canonical_json(logical_freeze) + "\n").encode(
                "utf-8"))
        _write_private(trace_path, trace_bytes)
        _write_private(result_path, result_bytes)
        try:
            summary = contract_api.validate_ledger(
                contract, "development", result_path, freeze_path,
                trace_path)
        except contract_api.ContractError as exc:
            fail(str(exc))
    if (receipt.get("validator_summary_sha256") !=
            contract_api.sha256_json(summary) or
            receipt.get("result_stream_sha256") !=
                recovery_api._hash_bytes(result_bytes)):
        fail("recomputed recovery summary differs from its execution receipt")
    return summary


def _load_recovery_and_work(
        contract: Mapping[str, Any], recovery_campaign_dir: Path,
        work_rank_dir: Path,
        ) -> Mapping[str, Any]:
    try:
        campaign = recovery_api.load_completed_campaign(
            contract, Path(recovery_campaign_dir))
        work = work_api.load_completed_work_screen(
            contract, Path(work_rank_dir))
    except (recovery_api.RecoveryRunnerError, work_api.WorkScreenError,
            contract_api.ContractError) as exc:
        fail(str(exc))
    freeze = campaign.get("freeze")
    rows = campaign.get("rows")
    run_summary = campaign.get("summary")
    receipt = campaign.get("receipt")
    if (not isinstance(freeze, Mapping) or not isinstance(rows, list) or
            not isinstance(run_summary, Mapping) or
            not isinstance(receipt, Mapping)):
        fail("strict recovery loader returned malformed reopened artifacts")
    try:
        work_binding = recovery_api._bind_work_rank_identities(
            rows, work, freeze.get("source_git_commit"))
    except recovery_api.RecoveryRunnerError as exc:
        fail(str(exc))
    if not isinstance(work_binding, Mapping):
        fail("work/rank identity join returned a malformed binding")
    for field in (
            "work_rank_summary_sha256", "work_rank_result_stream_sha256",
            "work_rank_domain_sha256"):
        if (work_binding.get(field) != freeze.get(field) or
                work_binding.get(field) != run_summary.get(field)):
            fail("independent work/rank binding differs in {}".format(field))
    summary = _revalidate_recovery_summary(contract, campaign)
    witness = campaign.get("timing_proxy_witness")
    if not isinstance(witness, Mapping):
        fail("strict recovery loader omitted the timing-proxy witness")
    witness_sha256 = contract_api.sha256_json(witness)
    if (witness_sha256 != freeze.get("timing_proxy_witness_sha256") or
            witness_sha256 != run_summary.get(
                "timing_proxy_witness_sha256")):
        fail("reopened timing-proxy witness differs from recovery evidence")
    return {
        "freeze": freeze,
        "summary": summary,
        "provenance": {
            "recovery_run_summary_sha256": run_summary.get("summary_sha256"),
            "recovery_result_stream_sha256":
                receipt.get("result_stream_sha256"),
            "recovery_execution_receipt_sha256":
                receipt.get("receipt_sha256"),
            "timing_proxy_witness_sha256": witness_sha256,
            "work_rank_summary_sha256":
                work_binding.get("work_rank_summary_sha256"),
            "work_rank_result_stream_sha256":
                work_binding.get("work_rank_result_stream_sha256"),
            "work_rank_domain_sha256":
                work_binding.get("work_rank_domain_sha256"),
        },
    }


def _load_timing(
        contract: Mapping[str, Any], timing_screen_dir: Path,
        ) -> Mapping[str, Any]:
    loader = getattr(timing_api, "load_completed_timing_screen", None)
    if not callable(loader):
        fail("timing runner lacks strict load_completed_timing_screen")
    try:
        loaded = loader(contract, Path(timing_screen_dir))
    except (timing_api.RunnerError, contract_api.ContractError) as exc:
        fail(str(exc))
    expected_keys = {
        "directory", "directory_identity", "run_summary", "freeze",
        "summary", "execution_receipt", "timing_qualification",
    }
    if not isinstance(loaded, dict) or set(loaded) != expected_keys:
        fail("strict timing loader returned the wrong interface")
    run_summary = loaded["run_summary"]
    freeze = loaded["freeze"]
    summary = loaded["summary"]
    receipt = loaded["execution_receipt"]
    qualification = loaded["timing_qualification"]
    if (not isinstance(run_summary, Mapping) or
            not isinstance(freeze, Mapping) or
            not isinstance(summary, Mapping) or
            not isinstance(receipt, Mapping) or
            type(qualification) is not contract_api.TimingQualification):
        fail("strict timing loader returned malformed reopened artifacts")
    expected = {
        "timing_freeze_sha256": contract_api.freeze_manifest_sha256(freeze),
        "timing_architecture_artifact_sha256":
            contract_api.architecture_artifact_sha256(freeze),
        "timing_result_sha256": receipt.get("result_stream_sha256"),
        "timing_execution_receipt_sha256": receipt.get("receipt_sha256"),
        "timing_validator_summary_sha256":
            contract_api.sha256_json(summary),
        "timing_qualification_execution_receipt_sha256": receipt.get(
            "timing_qualification_execution_receipt_sha256"),
    }
    if any(contract_api.canonical_json(run_summary.get(field)) !=
           contract_api.canonical_json(value)
           for field, value in expected.items()):
        fail("timing run summary differs from reopened timing artifacts")
    return {
        "freeze": freeze,
        "summary": summary,
        "timing_qualification": qualification,
        "provenance": {
            "timing_run_summary_sha256":
                run_summary.get("summary_sha256"),
            "timing_result_stream_sha256":
                receipt.get("result_stream_sha256"),
            "timing_execution_receipt_sha256":
                receipt.get("receipt_sha256"),
            "timing_qualification_execution_receipt_sha256": receipt.get(
                "timing_qualification_execution_receipt_sha256"),
        },
    }


def recompute_development_architecture_selection(
        contract: Mapping[str, Any], recovery_campaign_dir: Path,
        work_rank_dir: Path, timing_screen_dir: Path,
        ) -> Optional[Mapping[str, Any]]:
    """Strictly reopen all three artifact directories and recalculate v2."""
    recovery = _load_recovery_and_work(
        contract, recovery_campaign_dir, work_rank_dir)
    timing = _load_timing(contract, timing_screen_dir)
    provenance = {
        **recovery["provenance"],
        **timing["provenance"],
    }
    return contract_api._calculate_development_architecture_selection(
        contract, recovery["summary"], timing["summary"],
        timing["timing_qualification"], recovery["freeze"],
        timing["freeze"], provenance)


def _published_receipt(
        contract: Mapping[str, Any], path: Path,
        expected_identity: Optional[ReceiptIdentity] = None,
        ) -> Tuple[Mapping[str, Any], ReceiptFingerprint, bytes]:
    try:
        data, fingerprint = _read_stable_receipt_bytes(
            Path(path), expected_identity=expected_identity)
        value = contract_api._load_json_bytes(
            data, "architecture selection receipt")
        if not isinstance(value, dict):
            fail("architecture selection receipt must be a JSON object")
        logical = contract_api.canonical_json(value).encode("utf-8")
        if data not in (logical + b"\n", logical + b"\r\n"):
            fail("architecture selection receipt must be canonical JSON")
        return (contract_api.validate_selection_receipt(contract, value),
                fingerprint, data)
    except contract_api.ContractError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def select_development_architecture(
        contract: Mapping[str, Any], recovery_campaign_dir: Path,
        work_rank_dir: Path, timing_screen_dir: Path, output_path: Path,
        ) -> contract_api.ArchitectureSelection:
    """Recompute, seal while unnamed, then atomically publish one winner."""
    output_path = Path(output_path)
    if output_path.exists() or output_path.is_symlink():
        fail("refusing to replace existing selection receipt {}".format(
            output_path))
    result = recompute_development_architecture_selection(
        contract, recovery_campaign_dir, work_rank_dir, timing_screen_dir)
    if result is None:
        fail("development evidence produced no promotable architecture")
    transaction = _ReceiptPublicationTransaction(output_path)
    try:
        transaction.open_parent()
        selection = _prepare_receipt(transaction, contract, result)
        # os.link() below is the publish-wins linearization point.  There is no
        # rollback after it: an asynchronous exception may lose only the
        # acknowledgement, never the fully validated receipt.
        _publish_receipt(transaction)
        return selection
    finally:
        transaction.close()


def load_authoritative_selection(
        contract: Mapping[str, Any], recovery_campaign_dir: Path,
        work_rank_dir: Path, timing_screen_dir: Path,
        selection_receipt_path: Path,
        ) -> contract_api.ArchitectureSelection:
    """Recompute selection and exact-compare one previously published v2."""
    expected = recompute_development_architecture_selection(
        contract, recovery_campaign_dir, work_rank_dir, timing_screen_dir)
    if expected is None:
        fail("development evidence no longer produces a promotable winner")
    selection_receipt_path = Path(selection_receipt_path)
    published_identity = _visible_receipt_identity(selection_receipt_path)
    published, fingerprint, published_bytes = _published_receipt(
        contract, selection_receipt_path,
        expected_identity=published_identity)
    if contract_api.canonical_json(published) != contract_api.canonical_json(
            expected):
        fail("selection receipt differs from recomputed development artifacts")
    selection = contract_api._seal_validated_architecture_selection(
        contract, published)
    _require_receipt_snapshot(
        selection_receipt_path, fingerprint, published_bytes)
    return selection


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument("--recovery-campaign-dir", required=True, type=Path)
    parser.add_argument("--work-rank-dir", required=True, type=Path)
    parser.add_argument("--timing-screen-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        contract = contract_api.load_contract(args.contract)
        selection = select_development_architecture(
            contract, args.recovery_campaign_dir, args.work_rank_dir,
            args.timing_screen_dir, args.output)
        print(selection.selection_sha256)
        return 0
    except (SelectionError, contract_api.ContractError) as exc:
        print("selection error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
