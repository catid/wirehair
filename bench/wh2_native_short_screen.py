#!/usr/bin/env python3
"""Fail-closed assembly for the frozen WH2 native short screen.

The native worker emits canonical JSONL envelopes.  This module checks their
cardinality and runtime provenance, strips the envelopes to the exact frozen
ledger/timing schemas, runs the authoritative contract validator, and only
then publishes the result stream and its terminal execution receipt.

It is intentionally not a campaign scheduler.  Recovery work may be
distributed over the frozen workers, but v6 timing evidence must also prove
the frozen round-major homogeneous-wave CPU placement and non-overlapping
wave barriers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api


TRACE_RECORD_SCHEMA = "wirehair.wh2.native-trace-record.v1"
RECOVERY_RECORD_SCHEMA = "wirehair.wh2.native-recovery-record.v1"
RAW_RECOVERY_RECORD_SCHEMA = "wirehair.wh2.native-recovery-record.v2"
TIMING_RECORD_SCHEMA = "wirehair.wh2.native-timing-record.v4"
TIMING_QUALIFICATION_RECORD_SCHEMA = \
    "wirehair.wh2.native-timing-qualification-record.v1"
TIMING_QUALIFICATION_EXECUTION_SCHEMA = \
    "wirehair.wh2.native-timing-qualification-execution-receipt.v1"
SAMPLER_SCHEMA = "wirehair.wh2.sampler-attestation.v1"
EXECUTION_SCHEMA = "wirehair.wh2.native-execution-receipt.v1"
RAW_EXECUTION_SCHEMA = "wirehair.wh2.native-execution-receipt.v2"
TIMING_EXECUTION_SCHEMA = \
    "wirehair.wh2.native-timing-execution-receipt.v1"

SHA256_FIELDS = frozenset((
    "worker_binary_sha256", "message_sha256", "work_sha256",
))
NATIVE_RECORD_FIELDS = frozenset((
    "schema", "ordinal", "cpu", "worker_pid", "started_monotonic_ns",
    "finished_monotonic_ns", "worker_process_start_ticks",
    "worker_binary_sha256", "message_sha256", "work_sha256", "payload",
))
TRACE_RECORD_FIELDS = frozenset((
    "schema", "ordinal", "cell_sha256", "trace_sha256", "packet_count",
    "candidate_count",
))
SAMPLER_FIELDS = frozenset((
    "schema", "pid", "cpu", "process_start_ticks", "script_path",
    "script_sha256", "csv_path", "csv_device", "csv_inode",
    "window_start_monotonic_ns", "window_end_monotonic_ns",
    "terminal_status",
))
WORKER_FIELDS = frozenset(("cpu", "pid", "process_start_ticks"))
THERMAL_RECEIPT_FIELDS = frozenset((
    "pid", "cpu", "process_start_ticks", "script_path", "script_sha256",
    "csv_path", "csv_device", "csv_inode", "window_csv_sha256",
    "window_start_monotonic_ns", "window_end_monotonic_ns", "sample_count",
    "cpu_tctl_max_millic", "dimm_max_millic", "dimm_read_errors",
    "edac_ce_max", "edac_ue_max", "terminal_status",
))
SAMPLER_IDENTITY_FIELDS = (
    "pid", "cpu", "process_start_ticks", "script_path", "script_sha256",
    "csv_path", "csv_device", "csv_inode", "terminal_status",
)
EXECUTION_FIELDS = frozenset((
    "schema", "contract_sha256", "evidence_kind", "phase",
    "freeze_manifest_sha256", "trace_manifest_sha256",
    "result_stream_sha256", "record_count", "worker_start_monotonic_ns",
    "worker_end_monotonic_ns", "worker_cpus", "workers",
    "worker_binary_sha256s", "message_set_sha256", "work_set_sha256",
    "native_stream_sha256", "arm_descriptor_sha256s", "thermal",
    "validator_summary_sha256", "receipt_sha256",
))
TIMING_EXECUTION_FIELDS = EXECUTION_FIELDS | frozenset((
    "timing_base_domain_sha256", "timing_qualified_domain_sha256",
    "timing_qualification_map_sha256",
    "timing_qualification_audit_sha256",
    "timing_qualification_native_stream_sha256",
    "qualification_worker_start_monotonic_ns",
    "qualification_worker_end_monotonic_ns", "qualification_worker_cpus",
    "qualification_workers", "qualification_thermal",
    "timing_qualification_execution_receipt_sha256",
))
TIMING_QUALIFICATION_EXECUTION_FIELDS = frozenset((
    "schema", "contract_sha256", "phase", "source_git_commit",
    "timing_base_domain_sha256", "timing_qualified_domain_sha256",
    "timing_qualification_map_sha256",
    "timing_qualification_audit_sha256",
    "timing_qualification_native_stream_sha256",
    "qualification_attempt_count", "qualification_allowed_cpus",
    "qualification_worker_start_monotonic_ns",
    "qualification_worker_end_monotonic_ns", "qualification_worker_cpus",
    "qualification_workers", "qualification_worker_binary_sha256s",
    "receipt_sha256",
))
TIMING_QUALIFICATION_EMBEDDED_FIELDS = (
    "timing_qualification_native_stream_sha256",
    "qualification_worker_start_monotonic_ns",
    "qualification_worker_end_monotonic_ns", "qualification_worker_cpus",
    "qualification_workers",
)
TIMING_QUALIFICATION_PAYLOAD_FIELDS = frozenset((
    "ordinal", "base_cell_sha256", "loss_retry_offset", "loss_seed",
    "cell_sha256", "trace_sha256", "packet_count", "candidate_count",
    "wirehair2_head_outcome", "wirehair2_head_decoded_extra",
    "wirehair1_outcome", "wirehair1_decoded_extra",
))
TIMING_QUALIFICATION_MESSAGE_DOMAIN = \
    b"wirehair.wh2.timing-qualification-message.v1\0"
THERMAL_HEADER = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors", "load1",
    "load5", "load15", "edac_ce", "edac_ue",
)
MAX_THERMAL_LINE_BYTES = 65536
CPU_SYSFS_ROOT = Path("/sys/devices/system/cpu")


class NativeEvidenceError(ValueError):
    """The native evidence cannot be selected or published."""


def fail(message: str) -> None:
    raise NativeEvidenceError(message)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and contract_api.SHA256.fullmatch(value) \
        is not None


def _exact_keys(value: Any, fields: Iterable[str], context: str) \
        -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != set(fields):
        fail("{} has an unexpected schema".format(context))
    return value


def _sha256_file(path: Path) -> str:
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


def _open_regular_nofollow(path: Path, context: str) -> Tuple[int, os.stat_result]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    if nofollow == 0:
        fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(context))
    descriptor = -1
    try:
        descriptor = os.open(str(path), flags | nofollow)
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode):
            fail("{} must be a regular non-symlink file".format(context))
        return descriptor, info
    except OSError as exc:
        if descriptor >= 0:
            os.close(descriptor)
        fail("cannot open {} {}: {}".format(context, path, exc))
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        raise


def _sha256_descriptor(descriptor: int) -> str:
    digest = hashlib.sha256()
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            digest.update(block)
        os.lseek(descriptor, 0, os.SEEK_SET)
    except OSError as exc:
        fail("cannot hash an open evidence descriptor: {}".format(exc))
    return digest.hexdigest()


def _sha256_jsonl(values: Iterable[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(contract_api.canonical_json(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _load_canonical_object(path: Path, context: str) -> Mapping[str, Any]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    try:
        value = contract_api._load_json_bytes(data, context)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if not isinstance(value, dict):
        fail("{} must be an object".format(context))
    logical = contract_api.canonical_json(value).encode("utf-8")
    if data not in (logical + b"\n", logical + b"\r\n"):
        fail("{} must be canonical JSON followed by one line ending".format(
            context))
    return value


def _load_canonical_jsonl(path: Path, context: str) \
        -> List[Mapping[str, Any]]:
    try:
        return list(contract_api._parse_canonical_jsonl(path, context))
    except contract_api.ContractError as exc:
        fail(str(exc))
    return []


def _write_logical_jsonl(path: Path, values: Iterable[Mapping[str, Any]]) \
        -> None:
    with path.open("xb") as output:
        for value in values:
            output.write(contract_api.canonical_json(value).encode("utf-8"))
            output.write(b"\n")
        output.flush()
        os.fsync(output.fileno())


def _write_canonical_object(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("xb") as output:
        output.write(contract_api.canonical_json(value).encode("utf-8"))
        output.write(b"\n")
        output.flush()
        os.fsync(output.fileno())


def _temporary_path(destination: Path, suffix: str) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, raw = tempfile.mkstemp(
        prefix="." + destination.name + ".", suffix=suffix,
        dir=str(destination.parent))
    os.close(fd)
    os.unlink(raw)
    return Path(raw)


def _publish_new(
        staged: Path, destination: Path,
        expected_identity: Optional[Tuple[int, int]] = None) -> None:
    guard_fd, staged_info = _open_regular_nofollow(
        staged, "staged publication artifact")
    directory_fd = -1
    try:
        directory_fd = os.open(
            str(destination.parent),
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
            getattr(os, "O_CLOEXEC", 0))
        if (expected_identity is not None and
                expected_identity != (staged_info.st_dev, staged_info.st_ino)):
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
            # A successfully linked destination is already a valid immutable
            # dependency.  Never risk deleting a concurrent path replacement;
            # terminal receipts and run-summary are published last as markers.
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
    except FileNotFoundError:
        return False
    except OSError:
        return False


def _parse_physical_csv_line(data: bytes, context: str) -> Sequence[str]:
    if len(data) > MAX_THERMAL_LINE_BYTES:
        fail("{} exceeds the bounded physical line size".format(context))
    try:
        text = data.decode("ascii")
        rows = list(csv.reader([text], strict=True))
    except (UnicodeError, csv.Error) as exc:
        fail("{} is malformed: {}".format(context, exc))
    if len(rows) != 1:
        fail("{} is not exactly one physical CSV row".format(context))
    return rows[0]


def _trace_cells(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        timing_qualification: Optional[contract_api.TimingQualification] = None,
        ) -> Tuple[List[Mapping[str, Any]], int]:
    if evidence_kind == "recovery":
        cells = list(contract_api.iter_recovery_cells(contract, phase))
        count = contract["recovery"]["domains"][phase][
            "expected_cells_per_arm"]
    elif evidence_kind == "timing":
        if timing_qualification is None:
            fail("timing traces require a validated qualification map")
        cells = list(contract_api.iter_timing_cells(
            contract, phase, timing_qualification))
        count = contract["timing"]["domains"][phase]["expected_cells"]
    else:
        fail("evidence kind must be recovery or timing")
    if len(cells) != count:
        fail("contract cell iterator/cardinality mismatch")
    return cells, count


def assemble_trace_manifest(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        native_path: Path, output_path: Path,
        frozen_sha256: Optional[str] = None,
        timing_qualification: Optional[
            contract_api.TimingQualification] = None) -> str:
    """Validate native trace records and atomically publish a trace manifest."""
    cells, count = _trace_cells(
        contract, evidence_kind, phase, timing_qualification)
    records = _load_canonical_jsonl(native_path, "native trace stream")
    ordered: List[Optional[Mapping[str, Any]]] = [None] * count
    for raw in records:
        row = _exact_keys(raw, TRACE_RECORD_FIELDS, "native trace record")
        if row["schema"] != TRACE_RECORD_SCHEMA:
            fail("native trace record has an unknown schema")
        ordinal = row["ordinal"]
        if type(ordinal) is not int or not 0 <= ordinal < count:
            fail("native trace ordinal is outside the frozen domain")
        if ordered[ordinal] is not None:
            fail("duplicate native trace ordinal {}".format(ordinal))
        for field in ("cell_sha256", "trace_sha256"):
            if not _is_sha256(row[field]):
                fail("native trace {} is not a lowercase SHA-256".format(
                    field))
        if row["cell_sha256"] != contract_api.sha256_json(cells[ordinal]):
            fail("native trace cell hash does not bind its frozen ordinal")
        packet_count = row["packet_count"]
        candidate_count = row["candidate_count"]
        K = cells[ordinal]["K"]
        overhead = 4 if evidence_kind == "recovery" else \
            contract["timing"]["receive_overhead_cap"]
        if type(packet_count) is not int or packet_count != K + overhead:
            fail("native trace packet count differs from its frozen cap")
        candidate_limit = 256 * (K + overhead) + 65536
        if (type(candidate_count) is not int or
                not packet_count <= candidate_count <= candidate_limit):
            fail("native trace candidate count is outside its frozen bound")
        ordered[ordinal] = {
            "ordinal": ordinal,
            "cell_sha256": row["cell_sha256"],
            "trace_sha256": row["trace_sha256"],
        }
    if len(records) != count or any(value is None for value in ordered):
        fail("native trace stream has {} records, expected {}".format(
            len(records), count))
    manifest = [value for value in ordered if value is not None]
    staged = _temporary_path(output_path, ".jsonl")
    publication_identity: Optional[Tuple[int, int]] = None
    publication_guard = -1
    committed = False
    try:
        _write_logical_jsonl(staged, manifest)
        digest = contract_api.trace_manifest_sha256(
            contract, evidence_kind, phase, staged, timing_qualification)
        if frozen_sha256 is not None and digest != frozen_sha256:
            fail("native trace manifest differs from the frozen hash")
        publication_guard, info = _open_regular_nofollow(
            staged, "staged trace manifest")
        publication_identity = (info.st_dev, info.st_ino)
        _publish_new(staged, output_path, publication_identity)
        published_manifest = _load_canonical_jsonl(
            output_path, "published trace manifest")
        published_digest = contract_api.trace_manifest_sha256(
            contract, evidence_kind, phase, output_path,
            timing_qualification)
        if published_manifest != manifest or published_digest != digest:
            fail("published trace manifest changed before commit")
        committed = True
        guard = publication_guard
        try:
            os.close(guard)
        finally:
            publication_guard = -1
        return digest
    except BaseException:
        if not committed and publication_identity is not None:
            _unlink_if_identity(
                staged, publication_identity[0], publication_identity[1])
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


def _parse_proc_start_ticks(pid: int) -> int:
    try:
        data = Path("/proc/{}/stat".format(pid)).read_text(encoding="ascii")
    except OSError as exc:
        fail("cannot inspect live sampler PID {}: {}".format(pid, exc))
    close = data.rfind(")")
    fields = data[close + 2:].split() if close >= 0 else []
    # proc(5): field 22 is starttime.  `fields` starts at field 3.
    if len(fields) <= 19 or not fields[19].isdigit():
        fail("live sampler /proc stat is malformed")
    return int(fields[19])


def _require_process_predates_window(
        process_start_ticks: int, window_start_monotonic_ns: int,
        context: str) -> None:
    try:
        ticks_per_second = os.sysconf("SC_CLK_TCK")
    except (AttributeError, OSError, ValueError) as exc:
        fail("cannot read the process clock tick rate: {}".format(exc))
    if (type(ticks_per_second) is not int or ticks_per_second <= 0 or
            process_start_ticks <= 0 or window_start_monotonic_ns < 0 or
            process_start_ticks * 1000000000 >
                window_start_monotonic_ns * ticks_per_second):
        fail("{} process did not exist at the start of its evidence".format(
            context))


def _verify_live_sampler_process(
        pid: int, cpu: int, process_start_ticks: int,
        script_path: Path, csv_path: Path) -> None:
    if _parse_proc_start_ticks(pid) != process_start_ticks:
        fail("sampler PID was reused or its start time changed")
    try:
        affinity = os.sched_getaffinity(pid)
        cmdline = Path("/proc/{}/cmdline".format(pid)).read_bytes()
        executable = Path("/proc/{}/exe".format(pid)).resolve(strict=True)
    except (AttributeError, OSError) as exc:
        fail("cannot inspect live sampler process: {}".format(exc))
    if affinity != {cpu}:
        fail("sampler does not have singleton declared affinity")
    try:
        arguments = [os.fsdecode(value) for value in cmdline.split(b"\0")
                     if value]
    except UnicodeError as exc:
        fail("live sampler command line is malformed: {}".format(exc))
    try:
        expected_script = str(script_path.resolve(strict=True))
        expected_csv = str(csv_path.resolve(strict=True))
    except OSError as exc:
        fail("cannot resolve live sampler artifacts: {}".format(exc))
    if len(arguments) < 2 or not Path(arguments[0]).is_absolute():
        fail("live sampler command line has no absolute interpreter")
    try:
        command_executable = Path(arguments[0]).resolve(strict=True)
    except OSError as exc:
        fail("cannot resolve live sampler interpreter: {}".format(exc))
    if command_executable != executable or \
            not executable.name.lower().startswith("python"):
        fail("live sampler executable is not its declared Python interpreter")
    allowed_python_flags = frozenset((
        "-B", "-E", "-I", "-O", "-OO", "-q", "-s", "-S", "-u",
        "-v", "-x",
    ))
    script_index = 1
    while (script_index < len(arguments) and
           arguments[script_index] in allowed_python_flags):
        script_index += 1
    if script_index >= len(arguments) or \
            arguments[script_index] in ("-c", "-m"):
        fail("live sampler command line does not execute a script file")
    script_argument = Path(arguments[script_index])
    if not script_argument.is_absolute():
        fail("live sampler script argument must be absolute")
    try:
        executed_script = str(script_argument.resolve(strict=True))
    except OSError as exc:
        fail("cannot resolve executed sampler script: {}".format(exc))
    if executed_script != expected_script:
        fail("live sampler command line does not execute the attested script")
    if arguments.count("--csv") != 1:
        fail("live sampler command line must have exactly one --csv argument")
    try:
        csv_index = arguments.index("--csv")
    except ValueError:
        fail("live sampler command line has no --csv argument")
    if csv_index <= script_index or csv_index + 1 >= len(arguments):
        fail("live sampler command line has no --csv value")
    csv_argument = Path(arguments[csv_index + 1])
    if not csv_argument.is_absolute():
        fail("live sampler --csv argument must be absolute")
    try:
        actual_csv = str(csv_argument.resolve(strict=True))
    except OSError as exc:
        fail("cannot resolve live sampler --csv path: {}".format(exc))
    if actual_csv != expected_csv:
        fail("live sampler command line names a different CSV")
    if _parse_proc_start_ticks(pid) != process_start_ticks:
        fail("sampler exited or changed identity during verification")


def _verify_live_worker_process(
        pid: int, cpu: int, process_start_ticks: int,
        binary_sha256: str) -> None:
    if _parse_proc_start_ticks(pid) != process_start_ticks:
        fail("native worker PID was reused or its start time changed")
    try:
        affinity = os.sched_getaffinity(pid)
        descriptor = os.open(
            "/proc/{}/exe".format(pid),
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    except (AttributeError, OSError) as exc:
        fail("cannot inspect live native worker: {}".format(exc))
    try:
        info = os.fstat(descriptor)
        if affinity != {cpu} or not stat.S_ISREG(info.st_mode):
            fail("native worker is not the declared singleton-pinned executable")
        if _sha256_descriptor(descriptor) != binary_sha256:
            fail("native worker executable bytes differ from its envelope")
    finally:
        os.close(descriptor)
    if _parse_proc_start_ticks(pid) != process_start_ticks:
        fail("native worker exited or changed identity during verification")


def _finite_float(value: str, context: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        fail("{} is not numeric".format(context))
    if not math.isfinite(parsed):
        fail("{} is not finite".format(context))
    return parsed


def _parse_nonnegative_integer(value: str, context: str) -> int:
    if not value or not value.isascii() or not value.isdecimal() or \
            (len(value) > 1 and value.startswith("0")):
        fail("{} is not a canonical nonnegative integer".format(context))
    return int(value)


def _cpu_physical_core(
        cpu: int, sysfs_root: Optional[Path] = None) -> Tuple[int, int]:
    if type(cpu) is not int or cpu < 0:
        fail("CPU topology lookup requires a nonnegative integer")
    root = CPU_SYSFS_ROOT if sysfs_root is None else sysfs_root
    topology = root / "cpu{}".format(cpu) / "topology"
    values = []
    for name in ("physical_package_id", "core_id"):
        try:
            text = (topology / name).read_text(encoding="ascii").strip()
        except (OSError, UnicodeError) as exc:
            fail("cannot read CPU {} physical topology: {}".format(cpu, exc))
        values.append(_parse_nonnegative_integer(
            text, "CPU {} {}".format(cpu, name)))
    return values[0], values[1]


def _thermal_window(
        sampler: Mapping[str, Any], worker_start_ns: int, worker_end_ns: int,
        worker_cpus: Sequence[int], verify_live_sampler: bool,
        controller_cpu: Optional[int] = None,
        sysfs_root: Optional[Path] = None) \
        -> Mapping[str, Any]:
    sampler = _exact_keys(sampler, SAMPLER_FIELDS, "sampler attestation")
    if sampler["schema"] != SAMPLER_SCHEMA or \
            sampler["terminal_status"] != "ok":
        fail("sampler did not report the required terminal status")
    integer_fields = (
        "pid", "cpu", "process_start_ticks", "csv_device", "csv_inode",
        "window_start_monotonic_ns", "window_end_monotonic_ns",
    )
    if any(type(sampler[field]) is not int or sampler[field] < 0
           for field in integer_fields) or sampler["pid"] == 0:
        fail("sampler attestation has an invalid integer field")
    if not _is_sha256(sampler["script_sha256"]):
        fail("sampler script hash is not a lowercase SHA-256")
    for field in ("script_path", "csv_path"):
        if not isinstance(sampler[field], str) or not sampler[field]:
            fail("sampler {} must be a nonempty path".format(field))
    start_ns = sampler["window_start_monotonic_ns"]
    end_ns = sampler["window_end_monotonic_ns"]
    if not start_ns <= worker_start_ns <= worker_end_ns <= end_ns:
        fail("thermal window does not cover the complete worker interval")
    _require_process_predates_window(
        sampler["process_start_ticks"], start_ns, "sampler")
    worker_cores = {
        _cpu_physical_core(cpu, sysfs_root) for cpu in worker_cpus
    }
    sampler_core = _cpu_physical_core(sampler["cpu"], sysfs_root)
    if sampler_core in worker_cores:
        fail("sampler CPU overlaps a frozen worker physical core")
    if controller_cpu is not None:
        if type(controller_cpu) is not int or controller_cpu < 0:
            fail("frozen controller CPU is invalid")
        if sampler_core == _cpu_physical_core(controller_cpu, sysfs_root):
            fail("sampler CPU overlaps the controller physical core")

    script_path = Path(sampler["script_path"])
    csv_path = Path(sampler["csv_path"])
    script_fd = -1
    csv_fd = -1
    csv_read_fd = -1
    samples = []
    window_rows: List[Sequence[str]] = []
    try:
        script_fd, _ = _open_regular_nofollow(
            script_path, "sampler script")
        csv_fd, csv_stat = _open_regular_nofollow(csv_path, "sampler CSV")
        if (csv_stat.st_dev != sampler["csv_device"] or
                csv_stat.st_ino != sampler["csv_inode"]):
            fail("sampler CSV inode differs from its attestation")
        if _sha256_descriptor(script_fd) != sampler["script_sha256"]:
            fail("sampler script bytes differ from their attestation")
        if verify_live_sampler:
            _verify_live_sampler_process(
                sampler["pid"], sampler["cpu"],
                sampler["process_start_ticks"], script_path, csv_path)

        csv_read_fd = os.dup(csv_fd)
        source = os.fdopen(csv_read_fd, "rb")
        csv_read_fd = -1
        with source:
            header_bytes = source.readline(MAX_THERMAL_LINE_BYTES + 1)
            if not header_bytes.endswith(b"\n"):
                fail("sampler CSV header is not newline-terminated")
            header = _parse_physical_csv_line(
                header_bytes, "sampler CSV header")
            if header != list(THERMAL_HEADER):
                fail("sampler CSV header differs from the frozen format")
            prior_ns: Optional[int] = None
            inside_window = False
            row_number = 1
            while True:
                physical_line = source.readline(MAX_THERMAL_LINE_BYTES + 1)
                if not physical_line:
                    break
                row_number += 1
                terminated = physical_line.endswith(b"\n")
                try:
                    row = _parse_physical_csv_line(
                        physical_line,
                        "sampler CSV row {}".format(row_number))
                except NativeEvidenceError:
                    if not inside_window:
                        continue
                    raise
                # Historical rows precede the attested campaign window and
                # are not evidence for this run.  Ignore even malformed old
                # rows until the exact start sample is found; after that,
                # retain the shell runner's fail-closed behavior.
                if not inside_window:
                    if len(row) != len(THERMAL_HEADER):
                        continue
                    try:
                        candidate_s = float(row[1])
                    except (TypeError, ValueError):
                        continue
                    if not math.isfinite(candidate_s) or candidate_s < 0.0:
                        continue
                    candidate_ns = int(round(candidate_s * 1000000000.0))
                    if candidate_ns < start_ns:
                        continue
                    if candidate_ns != start_ns:
                        fail("sampler CSV skipped the exact window start")
                    inside_window = True
                if len(row) != len(THERMAL_HEADER):
                    fail("sampler CSV row {} has wrong width".format(row_number))
                mono_s = _finite_float(row[1], "sampler monotonic time")
                if mono_s < 0.0:
                    fail("sampler monotonic time is negative")
                mono_ns = int(round(mono_s * 1000000000.0))
                if mono_ns < start_ns:
                    fail("sampler window moved backwards before its start")
                if mono_ns > end_ns:
                    fail("sampler CSV skipped the exact window end")
                if prior_ns is not None and (
                        mono_ns <= prior_ns or mono_ns - prior_ns > 5000000000):
                    fail("sampler window is reordered or has a gap over 5 s")
                prior_ns = mono_ns
                busy = _finite_float(row[2], "CPU busy percentage")
                mhz = _finite_float(row[3], "CPU frequency")
                cpu_c = _finite_float(row[4], "CPU temperature")
                dimms = [
                    _finite_float(value, "DIMM temperature")
                    for value in row[5:13]
                ]
                errors = _parse_nonnegative_integer(
                    row[13], "DIMM read-error count")
                for value in row[14:17]:
                    if _finite_float(value, "load average") < 0.0:
                        fail("sampler load average is negative")
                ce = _parse_nonnegative_integer(row[17], "EDAC CE count")
                ue = _parse_nonnegative_integer(row[18], "EDAC UE count")
                if (not 0.0 <= busy <= 100.0 or mhz <= 0.0 or
                        not 0.0 <= cpu_c <= 120.0 or
                        any(not 0.0 < value <= 100.0 for value in dimms)):
                    fail("sampler window contains an out-of-range reading")
                if errors != 0:
                    fail("sampler window contains a DIMM read error")
                if ce != 0 or ue != 0:
                    fail("sampler window contains a nonzero EDAC counter")
                window_rows.append(row)
                samples.append((mono_ns, busy, mhz, cpu_c, max(dimms)))
                # The sampler is intentionally still live while evidence is
                # assembled.  Stop at the sealed endpoint so an append after
                # it (including a temporarily incomplete line) cannot race
                # validation of the already-complete window.
                if mono_ns == end_ns:
                    if not terminated:
                        fail("sampler CSV endpoint row is not newline-terminated")
                    break
    except (OSError, UnicodeError, csv.Error) as exc:
        fail("cannot parse sampler CSV: {}".format(exc))
    finally:
        if csv_read_fd >= 0:
            os.close(csv_read_fd)
        if csv_fd >= 0:
            os.close(csv_fd)
        if script_fd >= 0:
            os.close(script_fd)
    if verify_live_sampler:
        _verify_live_sampler_process(
            sampler["pid"], sampler["cpu"],
            sampler["process_start_ticks"], script_path, csv_path)
    if len(samples) < 2:
        fail("thermal window must contain at least two valid samples")
    if samples[0][0] != start_ns or samples[-1][0] != end_ns:
        fail("sampler window endpoints do not name exact CSV samples")
    window_digest = hashlib.sha256()
    for row in (THERMAL_HEADER,) + tuple(window_rows):
        window_digest.update(",".join(row).encode("ascii"))
        window_digest.update(b"\n")
    return {
        "pid": sampler["pid"],
        "cpu": sampler["cpu"],
        "process_start_ticks": sampler["process_start_ticks"],
        "script_path": sampler["script_path"],
        "script_sha256": sampler["script_sha256"],
        "csv_path": sampler["csv_path"],
        "csv_device": sampler["csv_device"],
        "csv_inode": sampler["csv_inode"],
        # Hash only the exact logical window.  The sole sampler remains live
        # and appends later rows, which must not invalidate a sealed receipt.
        "window_csv_sha256": window_digest.hexdigest(),
        "window_start_monotonic_ns": start_ns,
        "window_end_monotonic_ns": end_ns,
        "sample_count": len(samples),
        "cpu_tctl_max_millic": int(round(max(row[3] for row in samples) * 1000)),
        "dimm_max_millic": int(round(max(row[4] for row in samples) * 1000)),
        "dimm_read_errors": 0,
        "edac_ce_max": 0,
        "edac_ue_max": 0,
        "terminal_status": "ok",
    }


def _require_one_continuous_sampler(
        qualification: Mapping[str, Any], timing: Mapping[str, Any]) -> None:
    if any(qualification.get(field) != timing.get(field)
           for field in SAMPLER_IDENTITY_FIELDS):
        fail("qualification and timing thermal evidence splice samplers")
    if qualification.get("window_end_monotonic_ns") >= \
            timing.get("window_start_monotonic_ns"):
        fail("qualification and timing thermal windows are not sequential")


def write_sampler_attestation(
        pid: int, cpu: int, script_path: Path, csv_path: Path,
        window_start_monotonic_ns: int, window_end_monotonic_ns: int,
        output_path: Path) -> Mapping[str, Any]:
    """Bind a live, singleton-pinned sampler and one exact CSV window."""
    if (type(pid) is not int or pid <= 0 or type(cpu) is not int or cpu < 0 or
            type(window_start_monotonic_ns) is not int or
            type(window_end_monotonic_ns) is not int or
            not 0 <= window_start_monotonic_ns <= window_end_monotonic_ns):
        fail("invalid sampler attestation arguments")
    process_start_ticks = _parse_proc_start_ticks(pid)
    _require_process_predates_window(
        process_start_ticks, window_start_monotonic_ns, "sampler")
    _verify_live_sampler_process(
        pid, cpu, process_start_ticks, script_path, csv_path)
    try:
        resolved_script = str(script_path.resolve(strict=True))
        resolved_csv = str(csv_path.resolve(strict=True))
    except OSError as exc:
        fail("cannot resolve sampler artifacts: {}".format(exc))
    script_fd = -1
    csv_fd = -1
    try:
        script_fd, _ = _open_regular_nofollow(
            script_path, "sampler script")
        csv_fd, csv_info = _open_regular_nofollow(csv_path, "sampler CSV")
        value = {
            "schema": SAMPLER_SCHEMA,
            "pid": pid,
            "cpu": cpu,
            "process_start_ticks": process_start_ticks,
            "script_path": resolved_script,
            "script_sha256": _sha256_descriptor(script_fd),
            "csv_path": resolved_csv,
            "csv_device": csv_info.st_dev,
            "csv_inode": csv_info.st_ino,
            "window_start_monotonic_ns": window_start_monotonic_ns,
            "window_end_monotonic_ns": window_end_monotonic_ns,
            "terminal_status": "ok",
        }
    finally:
        if csv_fd >= 0:
            os.close(csv_fd)
        if script_fd >= 0:
            os.close(script_fd)
    _verify_live_sampler_process(
        pid, cpu, process_start_ticks, script_path, csv_path)
    staged = _temporary_path(output_path, ".json")
    publication_identity: Optional[Tuple[int, int]] = None
    publication_guard = -1
    committed = False
    try:
        _write_canonical_object(staged, value)
        publication_guard, info = _open_regular_nofollow(
            staged, "staged sampler attestation")
        publication_identity = (info.st_dev, info.st_ino)
        _publish_new(staged, output_path, publication_identity)
        published_value = _exact_keys(
            _load_canonical_object(
                output_path, "published sampler attestation"),
            SAMPLER_FIELDS, "published sampler attestation")
        if published_value != value:
            fail("published sampler attestation changed before commit")
        _verify_live_sampler_process(
            pid, cpu, process_start_ticks, script_path, csv_path)
        committed = True
        guard = publication_guard
        try:
            os.close(guard)
        finally:
            publication_guard = -1
        return value
    except BaseException:
        if not committed and publication_identity is not None:
            _unlink_if_identity(
                staged, publication_identity[0], publication_identity[1])
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


def _record_ordinal_context(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        evidence_kind: str, phase: str,
        timing_qualification: Optional[
            contract_api.TimingQualification] = None,
        ) -> Tuple[Any, Mapping[str, int]]:
    """Precompute immutable indexes shared by every row in one receipt."""
    if evidence_kind == "recovery":
        return (
            contract_api._domain_indexes(contract, phase),
            {arm: index for index, arm in enumerate(freeze["arm_roster"])},
        )
    if evidence_kind == "timing":
        if timing_qualification is None:
            fail("timing result indexing requires a validated qualification")
        panels = contract_api.timing_panels(contract, freeze["arm_roster"])
        return (
            contract_api._timing_cell_indexes(
                contract, phase, timing_qualification),
            {
                contract_api.canonical_json(panel): index
                for index, panel in enumerate(panels)
            },
        )
    fail("evidence kind must be recovery or timing")
    return (), {}


def _expected_record_ordinal(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        evidence_kind: str, phase: str, payload: Mapping[str, Any],
        ordinal_context: Tuple[Any, Mapping[str, int]]) \
        -> Tuple[int, str]:
    cell_indexes, item_indexes = ordinal_context
    if evidence_kind == "recovery":
        cell_ordinal, _ = contract_api._cell_ordinal(
            contract, phase, payload, cell_indexes)
        arm = payload.get("arm")
        arm_index = item_indexes.get(arm) if isinstance(arm, str) else None
        if arm_index is None:
            fail("native recovery payload arm is outside the frozen roster")
        return cell_ordinal * len(freeze["arm_roster"]) + arm_index, arm
    cell_ordinal, _ = contract_api._timing_cell_ordinal(
        contract, phase, payload, cell_indexes)
    panel = {
        "panel_kind": payload.get("panel_kind"),
        "scope": payload.get("scope"),
        "left_arm": payload.get("left_arm"),
        "right_arm": payload.get("right_arm"),
    }
    panel_index = item_indexes.get(contract_api.canonical_json(panel))
    if panel_index is None:
        fail("native timing payload has an undeclared or relabeled panel")
    return cell_ordinal * len(item_indexes) + panel_index, str(cell_ordinal)


def _expected_work_sha256(
        evidence_kind: str, phase: str, ordinal: int,
        payload: Mapping[str, Any]) -> str:
    return contract_api.sha256_json({
        "schema": "wirehair.wh2.native-work.v1",
        "evidence_kind": evidence_kind,
        "phase": phase,
        "ordinal": ordinal,
        "cell_sha256": payload["cell_sha256"],
    })


def recovery_record_schema_fields(
        freeze: Mapping[str, Any], payload: Mapping[str, Any],
        ) -> Tuple[str, frozenset]:
    """Select the exact native schema from the frozen per-arm seed policy."""
    if not isinstance(payload, Mapping):
        fail("native result payload has an unexpected schema")
    arm = payload.get("arm")
    if not isinstance(arm, str):
        fail("native recovery payload arm must be a string")
    arms = freeze.get("arms_by_name")
    frozen_arm = arms.get(arm) \
        if isinstance(arms, Mapping) else None
    if frozen_arm is None:
        fail("native recovery payload arm is outside the frozen roster")
    raw_arm = frozen_arm.get("construction_seed_basis") == \
        contract_api.RAW_CONSTRUCTION_SEED_BASIS
    if raw_arm:
        if freeze.get("schema") != contract_api.RAW_FREEZE_SCHEMA:
            fail("uniform raw payload is not bound by a raw v2 freeze")
        return (RAW_RECOVERY_RECORD_SCHEMA,
                contract_api.RAW_RECOVERY_RECORD_FIELDS)
    return RECOVERY_RECORD_SCHEMA, contract_api.LEDGER_FIELDS


def _reverify_native_workers(workers: Sequence[Mapping[str, Any]]) -> None:
    for worker in workers:
        _verify_live_worker_process(
            worker["pid"], worker["cpu"], worker["process_start_ticks"],
            worker["binary_sha256"])


def _timing_qualification_message_sha256(
        base_cell: Mapping[str, Any]) -> str:
    digest = hashlib.sha256()
    digest.update(TIMING_QUALIFICATION_MESSAGE_DOMAIN)
    digest.update(contract_api.canonical_json(base_cell).encode("utf-8"))
    return digest.hexdigest()


def _expected_timing_qualification_work_sha256(
        ordinal: int, cell_sha256: str) -> str:
    return contract_api.sha256_json({
        "schema": "wirehair.wh2.native-work.v1",
        "evidence_kind": "timing_qualification",
        "phase": "development",
        "ordinal": ordinal,
        "cell_sha256": cell_sha256,
    })


def _qualification_controls_as_dicts(
        controls: Sequence[Any]) -> List[Mapping[str, Any]]:
    result: List[Mapping[str, Any]] = []
    for value in controls:
        if isinstance(value, Mapping):
            result.append(dict(value))
        elif (isinstance(value, tuple) and
              len(value) == len(
                  contract_api.TIMING_QUALIFICATION_CONTROL_ORDER)):
            result.append(dict(zip(
                contract_api.TIMING_QUALIFICATION_CONTROL_ORDER, value)))
        else:
            fail("timing qualification control roster is malformed")
    return result


def _validate_timing_qualification_records(
        contract: Mapping[str, Any], phase: str,
        records: Sequence[Mapping[str, Any]], source_git_commit: str,
        controls: Sequence[Any], expected_cpus: Optional[Sequence[int]] = None,
        expected_retry_offsets: Optional[Sequence[int]] = None,
        verify_live_workers: bool = False,
        ) -> Tuple[List[Mapping[str, Any]], List[int], Mapping[str, Any]]:
    """Validate an out-of-order native qualification stream and order it."""
    if phase != "development":
        fail("native qualification currently supports development only")
    if (not isinstance(source_git_commit, str) or
            contract_api.GIT_COMMIT.fullmatch(source_git_commit) is None):
        fail("timing qualification source commit is malformed")
    control_values = _qualification_controls_as_dicts(controls)
    try:
        validated_controls = contract_api._validate_timing_qualification_controls(
            contract, phase, control_values)
    except contract_api.ContractError as exc:
        fail(str(exc))
    control_binaries = {value["binary_sha256"] for value in validated_controls}
    if len(control_binaries) != 1:
        fail("native qualification controls must use one worker binary")
    worker_binary_sha256 = next(iter(control_binaries))

    base_cells = list(contract_api.iter_timing_base_cells(contract, phase))
    expected_count = contract["timing"]["domains"][phase]["expected_cells"]
    if len(base_cells) != expected_count:
        fail("timing base-cell iterator/cardinality mismatch")
    by_attempt: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    attempts_by_base: List[set] = [set() for _ in range(expected_count)]
    cpu_to_pid: Dict[int, int] = {}
    pid_to_cpu: Dict[int, int] = {}
    pid_to_start_ticks: Dict[int, int] = {}
    worker_start: Optional[int] = None
    worker_end: Optional[int] = None
    runtime_workers: Dict[int, Mapping[str, Any]] = {}

    for raw in records:
        row = _exact_keys(
            raw, NATIVE_RECORD_FIELDS, "native timing qualification record")
        if row["schema"] != TIMING_QUALIFICATION_RECORD_SCHEMA:
            fail("native timing qualification record has an unknown schema")
        payload = _exact_keys(
            row["payload"], TIMING_QUALIFICATION_PAYLOAD_FIELDS,
            "native timing qualification payload")
        integer_fields = (
            "ordinal", "cpu", "worker_pid", "worker_process_start_ticks",
            "started_monotonic_ns", "finished_monotonic_ns",
        )
        if any(type(row[field]) is not int for field in integer_fields):
            fail("native timing qualification provenance is noncanonical")
        outer_ordinal = row["ordinal"]
        if not 0 <= outer_ordinal < expected_count * 256:
            fail("native timing qualification ordinal is outside the domain")
        base_ordinal, retry_offset = divmod(outer_ordinal, 256)
        if (type(payload["ordinal"]) is not int or
                payload["ordinal"] != base_ordinal or
                type(payload["loss_retry_offset"]) is not int or
                payload["loss_retry_offset"] != retry_offset):
            fail("native timing qualification ordinal/retry is inconsistent")
        key = (base_ordinal, retry_offset)
        if key in by_attempt:
            fail("duplicate native timing qualification attempt")

        base_cell = base_cells[base_ordinal]
        base_cell_sha256 = contract_api.sha256_json(base_cell)
        loss_seed = contract_api._qualified_timing_loss_seed(
            base_cell["base_loss_seed"], retry_offset)
        qualified_cell = dict(base_cell)
        qualified_cell["base_cell_sha256"] = base_cell_sha256
        qualified_cell["loss_retry_offset"] = retry_offset
        qualified_cell["loss_seed"] = loss_seed
        cell_sha256 = contract_api.sha256_json(qualified_cell)
        if (payload["base_cell_sha256"] != base_cell_sha256 or
                payload["loss_seed"] != loss_seed or
                payload["cell_sha256"] != cell_sha256):
            fail("native timing qualification substitutes its frozen cell")
        for field in ("base_cell_sha256", "cell_sha256", "trace_sha256"):
            if not _is_sha256(payload[field]):
                fail("native timing qualification {} is not a SHA-256".format(
                    field))
        packet_count = payload["packet_count"]
        candidate_count = payload["candidate_count"]
        expected_packets = base_cell["K"] + \
            contract["timing"]["receive_overhead_cap"]
        candidate_limit = 256 * expected_packets + 65536
        if (type(packet_count) is not int or
                packet_count != expected_packets or
                type(candidate_count) is not int or
                not packet_count <= candidate_count <= candidate_limit):
            fail("native timing qualification trace cardinality is invalid")
        try:
            outcomes = (
                contract_api._validate_timing_qualification_outcome(
                    payload, "wirehair2_head",
                    contract["timing"]["receive_overhead_cap"]),
                contract_api._validate_timing_qualification_outcome(
                    payload, "wirehair1",
                    contract["timing"]["receive_overhead_cap"]),
            )
        except contract_api.ContractError as exc:
            fail(str(exc))
        if outcomes not in (
                ("success", "success"),
                ("success", "need_more_at_bound"),
                ("need_more_at_bound", "success"),
                ("need_more_at_bound", "need_more_at_bound")):
            fail("native timing qualification has an invalid disposition")

        cpu = row["cpu"]
        pid = row["worker_pid"]
        start_ticks = row["worker_process_start_ticks"]
        started = row["started_monotonic_ns"]
        finished = row["finished_monotonic_ns"]
        if (cpu < 0 or pid <= 0 or start_ticks <= 0 or started < 0 or
                finished < started):
            fail("native timing qualification provenance is invalid")
        _require_process_predates_window(
            start_ticks, started, "native qualification worker")
        if cpu in cpu_to_pid and cpu_to_pid[cpu] != pid:
            fail("multiple qualification PIDs claim one CPU")
        if pid in pid_to_cpu and pid_to_cpu[pid] != cpu:
            fail("one qualification PID claims multiple CPUs")
        cpu_to_pid[cpu] = pid
        pid_to_cpu[pid] = cpu
        prior_ticks = pid_to_start_ticks.setdefault(pid, start_ticks)
        if prior_ticks != start_ticks:
            fail("qualification worker process start changed")
        if row["worker_binary_sha256"] != worker_binary_sha256:
            fail("qualification worker binary differs from both controls")
        for field in SHA256_FIELDS:
            if not _is_sha256(row[field]):
                fail("native timing qualification {} is not a SHA-256".format(
                    field))
        if row["message_sha256"] != \
                _timing_qualification_message_sha256(base_cell):
            fail("qualification message identity is not the canonical base cell")
        if row["work_sha256"] != \
                _expected_timing_qualification_work_sha256(
                    outer_ordinal, cell_sha256):
            fail("qualification work hash does not bind its exact attempt")
        worker_start = started if worker_start is None else min(
            worker_start, started)
        worker_end = finished if worker_end is None else max(
            worker_end, finished)
        runtime_workers[pid] = {
            "cpu": cpu, "pid": pid, "process_start_ticks": start_ticks,
            "binary_sha256": worker_binary_sha256,
        }
        by_attempt[key] = payload
        attempts_by_base[base_ordinal].add(retry_offset)

    if not records or worker_start is None or worker_end is None:
        fail("native timing qualification stream is empty")
    observed_cpus = sorted(cpu_to_pid)
    allowed_cpus = observed_cpus
    if expected_cpus is not None:
        requested = list(expected_cpus)
        if (requested != sorted(set(requested)) or
                any(type(cpu) is not int or cpu < 0 for cpu in requested)):
            fail("qualification CPU roster is malformed")
        if observed_cpus != requested:
            fail("qualification did not exercise every allowed logical CPU")
        allowed_cpus = requested

    audit: List[Mapping[str, Any]] = []
    retry_offsets: List[int] = []
    for base_ordinal in range(expected_count):
        selected: Optional[int] = None
        for retry_offset in range(256):
            payload = by_attempt.get((base_ordinal, retry_offset))
            if payload is None:
                break
            both_success = (
                payload["wirehair2_head_outcome"] == "success" and
                payload["wirehair1_outcome"] == "success")
            audit.append({
                field: payload[field]
                for field in contract_api.TIMING_QUALIFICATION_AUDIT_FIELDS
            })
            if both_success:
                selected = retry_offset
                break
        if selected is None:
            if (base_ordinal, 255) in by_attempt:
                fail("timing qualification exhausted retry offset 255")
            fail("timing qualification has a missing or unterminated base cell")
        if any(retry > selected for retry in attempts_by_base[base_ordinal]):
            fail("timing qualification speculatively ran a later retry")
        retry_offsets.append(selected)
    if len(by_attempt) != len(audit):
        fail("timing qualification contains an out-of-domain attempt")
    if expected_retry_offsets is not None and \
            tuple(retry_offsets) != tuple(expected_retry_offsets):
        fail("native timing qualification differs from the frozen retry map")

    workers = [runtime_workers[pid] for pid in sorted(
        runtime_workers, key=lambda value: runtime_workers[value]["cpu"])]
    if verify_live_workers:
        _reverify_native_workers(workers)
    metadata = {
        "qualification_attempt_count": len(records),
        "qualification_allowed_cpus": allowed_cpus,
        "qualification_worker_start_monotonic_ns": worker_start,
        "qualification_worker_end_monotonic_ns": worker_end,
        "qualification_worker_cpus": observed_cpus,
        "qualification_workers": [
            {key: worker[key] for key in WORKER_FIELDS}
            for worker in workers
        ],
        "qualification_worker_binary_sha256s": sorted({
            worker["binary_sha256"] for worker in workers
        }),
        "timing_qualification_native_stream_sha256":
            _sha256_jsonl(records),
        "_runtime_workers": workers,
    }
    return audit, retry_offsets, metadata


def load_timing_qualification_execution_receipt(
        contract: Mapping[str, Any], phase: str,
        qualification: contract_api.TimingQualification,
        native_path: Path, receipt_path: Path,
        expected_receipt_sha256: Optional[str] = None,
        expected_cpus: Optional[Sequence[int]] = None,
        verify_live_workers: bool = False,
        ) -> Tuple[Mapping[str, Any], Mapping[str, Any]]:
    """Strictly reopen qualification provenance committed while workers lived."""
    receipt = _exact_keys(
        _load_canonical_object(
            receipt_path, "timing qualification execution receipt"),
        TIMING_QUALIFICATION_EXECUTION_FIELDS,
        "timing qualification execution receipt")
    if receipt["schema"] != TIMING_QUALIFICATION_EXECUTION_SCHEMA:
        fail("timing qualification execution receipt has an unknown schema")
    unsigned = {
        key: value for key, value in receipt.items()
        if key != "receipt_sha256"
    }
    receipt_sha256 = receipt["receipt_sha256"]
    if (not _is_sha256(receipt_sha256) or
            receipt_sha256 != contract_api.sha256_json(unsigned)):
        fail("timing qualification execution receipt self-hash is invalid")
    if (expected_receipt_sha256 is not None and
            receipt_sha256 != expected_receipt_sha256):
        fail("timing qualification execution receipt differs from its identity")
    if (receipt["contract_sha256"] != contract_api.contract_sha256(contract) or
            receipt["phase"] != phase or
            receipt["source_git_commit"] != qualification.source_git_commit or
            receipt["timing_base_domain_sha256"] !=
                qualification.base_domain_sha256 or
            receipt["timing_qualified_domain_sha256"] !=
                qualification.qualified_domain_sha256 or
            receipt["timing_qualification_map_sha256"] !=
                qualification.map_sha256 or
            receipt["timing_qualification_audit_sha256"] !=
                qualification.qualification_audit_sha256):
        fail("timing qualification execution receipt binding is invalid")

    allowed_cpus = receipt["qualification_allowed_cpus"]
    worker_cpus = receipt["qualification_worker_cpus"]
    if (not isinstance(allowed_cpus, list) or not allowed_cpus or
            allowed_cpus != sorted(set(allowed_cpus)) or
            any(type(cpu) is not int or cpu < 0 for cpu in allowed_cpus) or
            worker_cpus != allowed_cpus):
        fail("timing qualification execution receipt CPU roster is invalid")
    if expected_cpus is not None and list(expected_cpus) != allowed_cpus:
        fail("timing qualification execution receipt changes allowed affinity")
    records = _load_canonical_jsonl(
        native_path, "native timing qualification stream")
    _, _, metadata_with_runtime = _validate_timing_qualification_records(
        contract, phase, records, qualification.source_git_commit,
        qualification.controls, expected_cpus=allowed_cpus,
        expected_retry_offsets=qualification.retry_offsets,
        verify_live_workers=verify_live_workers)
    metadata = {
        key: value for key, value in metadata_with_runtime.items()
        if key != "_runtime_workers"
    }
    for field in (
            "timing_qualification_native_stream_sha256",
            "qualification_attempt_count", "qualification_allowed_cpus",
            "qualification_worker_start_monotonic_ns",
            "qualification_worker_end_monotonic_ns",
            "qualification_worker_cpus", "qualification_workers",
            "qualification_worker_binary_sha256s"):
        if receipt[field] != metadata[field]:
            fail("timing qualification execution receipt {} differs from "
                 "native evidence".format(field))
    if (type(receipt["qualification_attempt_count"]) is not int or
            receipt["qualification_attempt_count"] <= 0 or
            type(receipt["qualification_worker_start_monotonic_ns"]) is not int or
            type(receipt["qualification_worker_end_monotonic_ns"]) is not int or
            receipt["qualification_worker_end_monotonic_ns"] <
                receipt["qualification_worker_start_monotonic_ns"] or
            not isinstance(receipt["qualification_workers"], list) or
            len(receipt["qualification_workers"]) != len(allowed_cpus) or
            any(not isinstance(worker, dict) or set(worker) != WORKER_FIELDS
                for worker in receipt["qualification_workers"]) or
            not isinstance(receipt["qualification_worker_binary_sha256s"], list) or
            not receipt["qualification_worker_binary_sha256s"] or
            receipt["qualification_worker_binary_sha256s"] != sorted(set(
                receipt["qualification_worker_binary_sha256s"])) or
            any(not _is_sha256(value) for value in
                receipt["qualification_worker_binary_sha256s"])):
        fail("timing qualification execution receipt provenance is malformed")
    return receipt, metadata


def _require_same_timing_qualification(
        expected: contract_api.TimingQualification,
        actual: contract_api.TimingQualification, context: str) -> None:
    if any(getattr(actual, field) != getattr(expected, field)
           for field in contract_api.TIMING_QUALIFICATION_FIELDS):
        fail("timing qualification changed {}".format(context))


def assemble_timing_qualification(
        contract: Mapping[str, Any], phase: str, native_path: Path,
        audit_path: Path, map_path: Path, execution_receipt_path: Path,
        source_git_commit: str,
        controls: Sequence[Mapping[str, Any]], expected_cpus: Sequence[int],
        verify_live_workers: bool = True,
        ) -> Tuple[contract_api.TimingQualification, Mapping[str, Any], str]:
    """Publish the canonical qualification triple while workers remain live."""
    if (audit_path.exists() or map_path.exists() or
            execution_receipt_path.exists()):
        fail("refusing to replace a timing qualification artifact")
    control_values = _qualification_controls_as_dicts(controls)
    records = _load_canonical_jsonl(
        native_path, "native timing qualification stream")
    audit, retry_offsets, metadata = _validate_timing_qualification_records(
        contract, phase, records, source_git_commit, control_values,
        expected_cpus, verify_live_workers=verify_live_workers)
    selected_traces: List[str] = []
    cursor = 0
    for retry_offset in retry_offsets:
        cursor += retry_offset
        selected_traces.append(audit[cursor]["trace_sha256"])
        cursor += 1
    qualified_domain_sha256 = contract_api._timing_domain_sha256_from_offsets(
        contract, phase, retry_offsets)

    staged_audit = _temporary_path(audit_path, ".jsonl")
    staged_map = _temporary_path(map_path, ".json")
    staged_receipt = _temporary_path(execution_receipt_path, ".json")
    published: List[Tuple[Path, int, int]] = []
    publication_guards: List[int] = []
    committed = False
    try:
        _write_logical_jsonl(staged_audit, audit)
        audit_sha256 = contract_api.timing_qualification_audit_sha256(
            contract, phase, staged_audit)
        map_value = {
            "schema": contract_api.TIMING_QUALIFICATION_MAP_SCHEMA,
            "contract_sha256": contract_api.contract_sha256(contract),
            "phase": phase,
            "source_git_commit": source_git_commit,
            "base_domain_sha256": contract["timing"]["domains"][phase][
                "base_domain_sha256"],
            "qualified_domain_sha256": qualified_domain_sha256,
            "entry_kind": contract_api.TIMING_QUALIFICATION_ENTRY_KIND,
            "controls": control_values,
            "qualification_audit_sha256": audit_sha256,
            "selected_trace_roster_sha256":
                contract_api.timing_selected_trace_roster_sha256(
                    selected_traces),
            "retry_offsets": retry_offsets,
        }
        _write_canonical_object(staged_map, map_value)
        map_sha256 = contract_api.timing_qualification_map_sha256(map_value)
        # Validate the staged triple before making any name visible.
        qualification = contract_api.load_timing_qualification_map(
            contract, phase, staged_map, staged_audit, map_sha256)
        clean_metadata = {
            key: value for key, value in metadata.items()
            if key != "_runtime_workers"
        }
        receipt: Dict[str, Any] = {
            "schema": TIMING_QUALIFICATION_EXECUTION_SCHEMA,
            "contract_sha256": contract_api.contract_sha256(contract),
            "phase": phase,
            "source_git_commit": source_git_commit,
            "timing_base_domain_sha256": qualification.base_domain_sha256,
            "timing_qualified_domain_sha256":
                qualification.qualified_domain_sha256,
            "timing_qualification_map_sha256": map_sha256,
            "timing_qualification_audit_sha256": audit_sha256,
            **clean_metadata,
        }
        receipt["receipt_sha256"] = contract_api.sha256_json(receipt)
        receipt_sha256 = receipt["receipt_sha256"]
        _write_canonical_object(staged_receipt, receipt)
        load_timing_qualification_execution_receipt(
            contract, phase, qualification, native_path, staged_receipt,
            expected_receipt_sha256=receipt_sha256,
            expected_cpus=expected_cpus,
            verify_live_workers=verify_live_workers)
        artifacts = (
            (staged_audit, audit_path), (staged_map, map_path),
            (staged_receipt, execution_receipt_path),
        )
        # Keep every staged inode guarded through the receipt-last commit.
        # Dependencies may remain after a failed run, but without the strict
        # receipt they never constitute terminal qualification evidence.
        published = []
        for staged, destination in artifacts:
            guard, info = _open_regular_nofollow(
                staged, "staged timing qualification artifact")
            publication_guards.append(guard)
            published.append((destination, info.st_dev, info.st_ino))
        for index in (0, 1):
            staged, destination = artifacts[index]
            _, device, inode = published[index]
            _publish_new(staged, destination, (device, inode))
        qualification = contract_api.load_timing_qualification_map(
            contract, phase, map_path, audit_path, map_sha256)
        precommit_receipt, precommit_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, qualification, native_path, staged_receipt,
                expected_receipt_sha256=receipt_sha256,
                expected_cpus=expected_cpus,
                verify_live_workers=verify_live_workers)
        if (precommit_receipt != receipt or
                precommit_metadata != clean_metadata):
            fail("timing qualification changed before receipt commit")
        precommit_qualification = contract_api.load_timing_qualification_map(
            contract, phase, map_path, audit_path, map_sha256)
        _require_same_timing_qualification(
            qualification, precommit_qualification,
            "immediately before receipt commit")
        precommit_receipt, precommit_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, precommit_qualification, native_path,
                staged_receipt, expected_receipt_sha256=receipt_sha256,
                expected_cpus=expected_cpus,
                verify_live_workers=verify_live_workers)
        if (precommit_receipt != receipt or
                precommit_metadata != clean_metadata):
            fail("timing qualification changed immediately before receipt "
                 "commit")
        _, receipt_device, receipt_inode = published[2]
        _publish_new(
            staged_receipt, execution_receipt_path,
            (receipt_device, receipt_inode))
        committed = True
        final_qualification = contract_api.load_timing_qualification_map(
            contract, phase, map_path, audit_path, map_sha256)
        _require_same_timing_qualification(
            precommit_qualification, final_qualification,
            "during receipt commit")
        final_receipt, final_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, final_qualification, native_path,
                execution_receipt_path,
                expected_receipt_sha256=receipt_sha256,
                expected_cpus=expected_cpus,
                verify_live_workers=verify_live_workers)
        if final_receipt != receipt or final_metadata != clean_metadata:
            fail("timing qualification changed across triple publication")
        terminal_qualification = contract_api.load_timing_qualification_map(
            contract, phase, map_path, audit_path, map_sha256)
        _require_same_timing_qualification(
            final_qualification, terminal_qualification,
            "during terminal qualification reopen")
        terminal_receipt, terminal_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, terminal_qualification, native_path,
                execution_receipt_path,
                expected_receipt_sha256=receipt_sha256,
                expected_cpus=expected_cpus,
                verify_live_workers=verify_live_workers)
        if (terminal_receipt != receipt or
                terminal_metadata != clean_metadata):
            fail("timing qualification changed during terminal reopen")
        while publication_guards:
            guard = publication_guards[-1]
            try:
                os.close(guard)
            finally:
                publication_guards.pop()
        return terminal_qualification, terminal_metadata, receipt_sha256
    except BaseException:
        if not committed:
            staged_identities = {
                staged: (device, inode)
                for staged, (_, device, inode) in zip(
                    (staged_audit, staged_map, staged_receipt), published)
            }
            for path in (staged_audit, staged_map, staged_receipt):
                identity = staged_identities.get(path)
                if identity is not None:
                    _unlink_if_identity(path, identity[0], identity[1])
                else:
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
        while publication_guards:
            guard = publication_guards[-1]
            try:
                os.close(guard)
            except OSError:
                pass
            finally:
                publication_guards.pop()
        raise


def publish_timing_trace_manifest(
        contract: Mapping[str, Any], phase: str,
        qualification: contract_api.TimingQualification,
        output_path: Path) -> str:
    """Publish the selected qualification traces as trace-manifest v2."""
    cells = list(contract_api.iter_timing_cells(
        contract, phase, qualification))
    if len(cells) != len(qualification.selected_trace_sha256s):
        fail("timing qualification trace roster has the wrong cardinality")
    rows = [
        {
            "ordinal": ordinal,
            "cell_sha256": contract_api.sha256_json(cell),
            "trace_sha256": qualification.selected_trace_sha256s[ordinal],
        }
        for ordinal, cell in enumerate(cells)
    ]
    staged = _temporary_path(output_path, ".jsonl")
    publication_identity: Optional[Tuple[int, int]] = None
    publication_guard = -1
    committed = False
    try:
        _write_logical_jsonl(staged, rows)
        try:
            digest = contract_api.trace_manifest_sha256(
                contract, "timing", phase, staged, qualification)
            contract_api.load_timing_trace_manifest(
                contract, phase, staged, digest, qualification)
        except contract_api.ContractError as exc:
            fail(str(exc))
        publication_guard, info = _open_regular_nofollow(
            staged, "staged timing trace manifest")
        publication_identity = (info.st_dev, info.st_ino)
        _publish_new(staged, output_path, publication_identity)
        published_rows = _load_canonical_jsonl(
            output_path, "published timing trace manifest")
        try:
            published_digest = contract_api.trace_manifest_sha256(
                contract, "timing", phase, output_path, qualification)
            contract_api.load_timing_trace_manifest(
                contract, phase, output_path, published_digest, qualification)
        except contract_api.ContractError as exc:
            fail(str(exc))
        if published_rows != rows or published_digest != digest:
            fail("published timing trace manifest changed before commit")
        committed = True
        guard = publication_guard
        try:
            os.close(guard)
        finally:
            publication_guard = -1
        return digest
    except BaseException:
        if not committed and publication_identity is not None:
            _unlink_if_identity(
                staged, publication_identity[0], publication_identity[1])
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


def _timing_worker_roster(
        contract: Mapping[str, Any], freeze: Mapping[str, Any]) -> List[int]:
    timing = contract.get("timing") if isinstance(contract, Mapping) else None
    geometry = timing.get("execution_geometry") \
        if isinstance(timing, Mapping) else None
    if geometry != contract_api.EXPECTED_TIMING_EXECUTION_GEOMETRY:
        fail("timing execution geometry is not frozen")
    worker_count = geometry.get("timing_worker_count")
    jobs_per_wave = geometry.get("jobs_per_wave")
    roster = freeze.get("cpu_affinity") if isinstance(freeze, Mapping) else None
    if (type(worker_count) is not int or worker_count != 8 or
            type(jobs_per_wave) is not int or jobs_per_wave != worker_count or
            not isinstance(roster, list) or len(roster) != worker_count or
            any(type(cpu) is not int or cpu < 0 for cpu in roster) or
            roster != sorted(set(roster))):
        fail("timing freeze must bind exactly eight sorted worker CPUs")
    return list(roster)


def _validate_timing_execution_geometry(
        contract: Mapping[str, Any], freeze: Mapping[str, Any], phase: str,
        envelopes: Sequence[Mapping[str, Any]],
        timing_qualification: contract_api.TimingQualification,
        sysfs_root: Optional[Path] = None) -> None:
    """Verify exact v5 round-major placement and every inter-wave barrier."""
    worker_cpus = _timing_worker_roster(contract, freeze)
    worker_cores = [
        _cpu_physical_core(cpu, sysfs_root) for cpu in worker_cpus
    ]
    if len(set(worker_cores)) != len(worker_cores):
        fail("timing worker CPUs are not eight unique physical cores")
    host_identity = freeze.get("host_identity") \
        if isinstance(freeze, Mapping) else None
    controller_cpu = host_identity.get("controller_cpu") \
        if isinstance(host_identity, Mapping) else None
    if type(controller_cpu) is not int or controller_cpu < 0:
        fail("timing freeze does not bind a valid controller CPU")
    if _cpu_physical_core(controller_cpu, sysfs_root) in set(worker_cores):
        fail("timing controller shares a worker physical core")
    geometry = contract["timing"]["execution_geometry"]
    jobs_per_wave = geometry["jobs_per_wave"]
    domain = contract["timing"]["domains"].get(phase)
    repetitions = domain.get("paired_repetitions") \
        if isinstance(domain, Mapping) else None
    independent_rounds = domain.get("independent_rounds") \
        if isinstance(domain, Mapping) else None
    expected_cells = domain.get("expected_cells") \
        if isinstance(domain, Mapping) else None
    if (type(repetitions) is not int or repetitions <= 0 or
            repetitions % jobs_per_wave != 0 or
            type(independent_rounds) is not int or independent_rounds <= 0 or
            independent_rounds * jobs_per_wave != repetitions or
            type(expected_cells) is not int or expected_cells <= 0):
        fail("timing wave domain has invalid cardinality")
    try:
        cells = list(contract_api.iter_timing_cells(
            contract, phase, timing_qualification))
        panels = contract_api.timing_panels(contract, freeze["arm_roster"])
    except (KeyError, contract_api.ContractError) as exc:
        fail(str(exc))
    if len(cells) != expected_cells or not panels:
        fail("timing wave domain differs from its frozen cardinality")

    cell_fields = contract["timing"].get("cell_key")
    if (not isinstance(cell_fields, list) or
            "replicate" not in cell_fields or "loss_seed" not in cell_fields):
        fail("timing cohort identity fields are invalid")
    stable_fields = [
        field for field in cell_fields
        if field not in (
            "replicate", "base_loss_seed", "base_cell_sha256",
            "loss_retry_offset", "loss_seed")
    ]
    identity_order: List[str] = []
    cell_identities: List[str] = []
    for cell in cells:
        if not isinstance(cell, Mapping) or set(cell) != set(cell_fields):
            fail("timing cell has an unexpected cohort identity")
        identity = contract_api.canonical_json({
            field: cell[field] for field in stable_fields
        })
        cell_identities.append(identity)
        if cell.get("replicate") == 0:
            identity_order.append(identity)
    if len(identity_order) != len(set(identity_order)):
        fail("timing replicate zero has duplicate cohort identities")
    identity_indexes = {
        identity: index for index, identity in enumerate(identity_order)
    }
    if (not identity_order or
            any(identity not in identity_indexes for identity in cell_identities)):
        fail("timing cohort identity differs across replicates")

    panel_count = len(panels)
    cohort_count = len(identity_order) * panel_count
    try:
        frozen_cohort_count = contract_api.timing_cohort_count(
            contract, phase, freeze["arm_roster"])
    except (KeyError, contract_api.ContractError) as exc:
        fail(str(exc))
    if cohort_count != frozen_cohort_count:
        fail("timing cohort count differs from frozen execution geometry")
    expected_wave_count = cohort_count * independent_rounds
    expected_record_count = expected_cells * panel_count
    if len(envelopes) != expected_record_count:
        fail("timing execution has an invalid record cardinality")

    waves: List[List[Tuple[int, int, int]]] = [
        [] for _ in range(expected_wave_count)
    ]
    for ordinal, row in enumerate(envelopes):
        if row.get("ordinal") != ordinal:
            fail("timing execution envelopes are not ordinal-complete")
        cell_ordinal, panel_index = divmod(ordinal, panel_count)
        cell = cells[cell_ordinal]
        replicate = cell.get("replicate")
        if type(replicate) is not int or not 0 <= replicate < repetitions:
            fail("timing execution has an invalid replicate")
        stable_index = identity_indexes[cell_identities[cell_ordinal]]
        cohort_index = stable_index * panel_count + panel_index
        round_index = replicate // jobs_per_wave
        # The controller traverses the complete cohort domain once per
        # independent round.  Preserve that chronology in the terminal proof
        # instead of regrouping concurrent lanes as independent samples.
        wave_ordinal = round_index * cohort_count + cohort_index
        try:
            worker_slot = contract_api.timing_worker_slot(
                contract, phase, freeze["arm_roster"], cohort_index,
                replicate)
        except contract_api.ContractError as exc:
            fail(str(exc))
        if (type(worker_slot) is not int or
                not 0 <= worker_slot < len(worker_cpus) or
                row.get("cpu") != worker_cpus[worker_slot]):
            fail("timing result CPU differs from frozen worker-slot placement")
        waves[wave_ordinal].append((
            replicate, row["started_monotonic_ns"],
            row["finished_monotonic_ns"]))

    prior_wave_finish: Optional[int] = None
    for wave_ordinal, values in enumerate(waves):
        round_index = wave_ordinal // cohort_count
        expected_replicates = list(range(
            round_index * jobs_per_wave,
            (round_index + 1) * jobs_per_wave))
        if (len(values) != jobs_per_wave or
                sorted(value[0] for value in values) != expected_replicates):
            fail("timing wave does not contain its exact replicate roster")
        minimum_start = min(value[1] for value in values)
        maximum_finish = max(value[2] for value in values)
        if prior_wave_finish is not None and \
                prior_wave_finish > minimum_start:
            fail("timing waves overlap across a frozen barrier")
        prior_wave_finish = maximum_finish


def _validate_native_records(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        evidence_kind: str, phase: str,
        records: Sequence[Mapping[str, Any]],
        timing_qualification: Optional[
            contract_api.TimingQualification] = None,
        verify_live_workers: bool = False,
        sysfs_root: Optional[Path] = None) \
        -> Tuple[List[Mapping[str, Any]], Mapping[str, Any]]:
    if evidence_kind == "recovery":
        cell_count = contract["recovery"]["domains"][phase][
            "expected_cells_per_arm"]
        expected_count = cell_count * len(freeze["arm_roster"])
    elif evidence_kind == "timing":
        if timing_qualification is None:
            fail("timing native records require a validated qualification")
        _timing_worker_roster(contract, freeze)
        payload_fields = contract_api.TIMING_RECEIPT_FIELDS
        record_schema = TIMING_RECORD_SCHEMA
        cell_count = contract["timing"]["domains"][phase]["expected_cells"]
        expected_count = cell_count * len(contract_api.timing_panels(
            contract, freeze["arm_roster"]))
    else:
        fail("evidence kind must be recovery or timing")
    try:
        ordinal_context = _record_ordinal_context(
            contract, freeze, evidence_kind, phase, timing_qualification)
    except contract_api.ContractError as exc:
        fail(str(exc))
    ordered: List[Optional[Mapping[str, Any]]] = [None] * expected_count
    envelope_ordered: List[Optional[Mapping[str, Any]]] = [None] * expected_count
    worker_start: Optional[int] = None
    worker_end: Optional[int] = None
    cpu_to_pid: Dict[int, int] = {}
    pid_to_cpu: Dict[int, int] = {}
    pid_to_start_ticks: Dict[int, int] = {}
    pid_to_binary: Dict[int, str] = {}
    message_by_cell: Dict[str, str] = {}
    allowed_cpus = set(freeze["cpu_affinity"])
    worker_binaries = set()
    work_by_ordinal: List[Optional[str]] = [None] * expected_count
    for raw in records:
        row = _exact_keys(raw, NATIVE_RECORD_FIELDS, "native result record")
        payload = row["payload"]
        if evidence_kind == "recovery":
            record_schema, payload_fields = recovery_record_schema_fields(
                freeze, payload)
        if row["schema"] != record_schema:
            fail("native result record has an unknown schema")
        _exact_keys(payload, payload_fields, "native result payload")
        try:
            expected_ordinal, message_key = _expected_record_ordinal(
                contract, freeze, evidence_kind, phase, payload,
                ordinal_context)
        except contract_api.ContractError as exc:
            fail(str(exc))
        ordinal = row["ordinal"]
        if type(ordinal) is not int or ordinal != expected_ordinal:
            fail("native result ordinal does not bind its payload identity")
        if ordered[ordinal] is not None:
            fail("duplicate native result ordinal {}".format(ordinal))
        cpu = row["cpu"]
        pid = row["worker_pid"]
        process_start_ticks = row["worker_process_start_ticks"]
        start = row["started_monotonic_ns"]
        end = row["finished_monotonic_ns"]
        if (type(cpu) is not int or cpu not in allowed_cpus or
                type(pid) is not int or pid <= 0 or
                type(process_start_ticks) is not int or
                process_start_ticks <= 0 or
                type(start) is not int or type(end) is not int or
                start < 0 or end < start):
            fail("native result has invalid worker/core/time provenance")
        _require_process_predates_window(
            process_start_ticks, start, "native worker")
        if cpu in cpu_to_pid and cpu_to_pid[cpu] != pid:
            fail("multiple worker PIDs claim one frozen CPU")
        if pid in pid_to_cpu and pid_to_cpu[pid] != cpu:
            fail("one worker PID claims multiple frozen CPUs")
        cpu_to_pid[cpu] = pid
        pid_to_cpu[pid] = cpu
        prior_start_ticks = pid_to_start_ticks.setdefault(
            pid, process_start_ticks)
        if prior_start_ticks != process_start_ticks:
            fail("native worker process start changed within the stream")
        for field in SHA256_FIELDS:
            if not _is_sha256(row[field]):
                fail("native result {} is not a lowercase SHA-256".format(
                    field))
        worker_binary = row["worker_binary_sha256"]
        prior_pid_binary = pid_to_binary.setdefault(pid, worker_binary)
        if prior_pid_binary != worker_binary:
            fail("native worker PID changed executable identity")
        if evidence_kind == "recovery":
            expected_binary = payload["binary_sha256"]
            cell_key = str(ordinal // len(freeze["arm_roster"]))
        else:
            if (payload["left_binary_sha256"] !=
                    payload["right_binary_sha256"]):
                fail("timing sides must run in the same worker executable")
            expected_binary = payload["left_binary_sha256"]
            cell_key = message_key
        if worker_binary != expected_binary:
            fail("native worker binary differs from its result payload")
        prior_message = message_by_cell.setdefault(
            cell_key, row["message_sha256"])
        if prior_message != row["message_sha256"]:
            fail("native arms/panels changed the deterministic cell message")
        worker_binaries.add(worker_binary)
        expected_work = _expected_work_sha256(
            evidence_kind, phase, ordinal, payload)
        if row["work_sha256"] != expected_work:
            fail("native work hash does not bind its exact frozen job")
        work_by_ordinal[ordinal] = row["work_sha256"]
        worker_start = start if worker_start is None else min(worker_start, start)
        worker_end = end if worker_end is None else max(worker_end, end)
        ordered[ordinal] = payload
        envelope_ordered[ordinal] = row
    if len(records) != expected_count or any(value is None for value in ordered):
        fail("native result stream has {} records, expected {}".format(
            len(records), expected_count))
    if set(cpu_to_pid) != allowed_cpus:
        fail("native result stream did not use every frozen worker CPU")
    if worker_start is None or worker_end is None:
        fail("native result stream has no worker interval")
    if any(value is None for value in work_by_ordinal):
        fail("native result stream has an incomplete work-hash roster")
    if evidence_kind == "timing":
        _validate_timing_execution_geometry(
            contract, freeze, phase,
            [value for value in envelope_ordered if value is not None],
            timing_qualification,
            sysfs_root)
    runtime_workers = [
        {
            "cpu": cpu,
            "pid": cpu_to_pid[cpu],
            "process_start_ticks": pid_to_start_ticks[cpu_to_pid[cpu]],
            "binary_sha256": pid_to_binary[cpu_to_pid[cpu]],
        }
        for cpu in sorted(cpu_to_pid)
    ]
    if verify_live_workers:
        _reverify_native_workers(runtime_workers)
    payloads = [value for value in ordered if value is not None]
    envelopes = [value for value in envelope_ordered if value is not None]
    metadata = {
        "record_count": expected_count,
        "worker_start_monotonic_ns": worker_start,
        "worker_end_monotonic_ns": worker_end,
        "worker_cpus": sorted(cpu_to_pid),
        "workers": [
            {
                "cpu": cpu,
                "pid": cpu_to_pid[cpu],
                "process_start_ticks": pid_to_start_ticks[cpu_to_pid[cpu]],
            }
            for cpu in sorted(cpu_to_pid)
        ],
        "worker_binary_sha256s": sorted(worker_binaries),
        "message_set_sha256": contract_api.sha256_json(message_by_cell),
        "work_set_sha256": contract_api.sha256_json(work_by_ordinal),
        "native_stream_sha256": _sha256_jsonl(records),
        "_runtime_workers": runtime_workers,
    }
    return payloads, metadata


def _load_required_timing_qualification(
        contract: Mapping[str, Any], phase: str,
        map_path: Optional[Path], audit_path: Optional[Path],
        map_sha256: Optional[str]) -> contract_api.TimingQualification:
    if map_path is None or audit_path is None or map_sha256 is None:
        fail("timing evidence requires qualification map/audit/hash arguments")
    try:
        return contract_api.load_timing_qualification_map(
            contract, phase, map_path, audit_path, map_sha256)
    except contract_api.ContractError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def assemble_results(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        freeze_path: Path, trace_path: Path, native_path: Path,
        sampler_path: Path, output_path: Path, execution_receipt_path: Path,
        verify_live_sampler: bool = True,
        timing_qualification_map_path: Optional[Path] = None,
        timing_qualification_audit_path: Optional[Path] = None,
        timing_qualification_map_sha256: Optional[str] = None,
        timing_qualification_native_path: Optional[Path] = None,
        timing_qualification_sampler_path: Optional[Path] = None,
        timing_qualification_execution_receipt_path: Optional[Path] = None,
        timing_qualification_execution_receipt_sha256: Optional[str] = None,
        ) -> Mapping[str, Any]:
    """Validate and publish one complete recovery ledger or timing receipt."""
    if output_path.exists() or execution_receipt_path.exists():
        fail("refusing to replace an existing result or execution receipt")
    timing_qualification: Optional[contract_api.TimingQualification] = None
    qualification_metadata: Optional[Mapping[str, Any]] = None
    qualification_sampler: Optional[Mapping[str, Any]] = None
    qualification_execution_receipt: Optional[Mapping[str, Any]] = None
    if evidence_kind == "timing":
        timing_qualification = _load_required_timing_qualification(
            contract, phase, timing_qualification_map_path,
            timing_qualification_audit_path,
            timing_qualification_map_sha256)
        if timing_qualification_native_path is None:
            fail("timing evidence requires the native qualification stream")
        if timing_qualification_sampler_path is None:
            fail("timing evidence requires a qualification thermal attestation")
        if (timing_qualification_execution_receipt_path is None or
                timing_qualification_execution_receipt_sha256 is None):
            fail("timing evidence requires the qualification execution receipt")
        qualification_sampler = _load_canonical_object(
            timing_qualification_sampler_path,
            "qualification sampler attestation")
        qualification_execution_receipt, qualification_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, timing_qualification,
                timing_qualification_native_path,
                timing_qualification_execution_receipt_path,
                expected_receipt_sha256=
                    timing_qualification_execution_receipt_sha256)
    try:
        freeze = contract_api.load_freeze_manifest(
            contract, phase, freeze_path, evidence_kind,
            timing_qualification)
    except contract_api.ContractError as exc:
        fail(str(exc))
    records = _load_canonical_jsonl(native_path, "native result stream")
    payloads, metadata_with_runtime = _validate_native_records(
        contract, freeze, evidence_kind, phase, records,
        timing_qualification,
        verify_live_workers=verify_live_sampler)
    runtime_workers = metadata_with_runtime["_runtime_workers"]
    metadata = {
        key: value for key, value in metadata_with_runtime.items()
        if key != "_runtime_workers"
    }
    if qualification_metadata is not None:
        if qualification_metadata[
                "qualification_worker_end_monotonic_ns"] > \
                metadata["worker_start_monotonic_ns"]:
            fail("qualification workers overlap the exact-eight timing interval")
        qualification_identities = {
            (value["pid"], value["process_start_ticks"])
            for value in qualification_metadata["qualification_workers"]
        }
        timing_identities = {
            (value["pid"], value["process_start_ticks"])
            for value in metadata["workers"]
        }
        if qualification_identities & timing_identities:
            fail("a qualification worker identity survived into timing")
    sampler = _load_canonical_object(sampler_path, "sampler attestation")
    thermal = _thermal_window(
        sampler, metadata["worker_start_monotonic_ns"],
        metadata["worker_end_monotonic_ns"], metadata["worker_cpus"],
        verify_live_sampler,
        controller_cpu=freeze.get("host_identity", {}).get("controller_cpu"))
    qualification_thermal: Optional[Mapping[str, Any]] = None
    if qualification_metadata is not None and \
            qualification_sampler is not None:
        qualification_thermal = _thermal_window(
            qualification_sampler,
            qualification_metadata[
                "qualification_worker_start_monotonic_ns"],
            qualification_metadata[
                "qualification_worker_end_monotonic_ns"],
            qualification_metadata["qualification_worker_cpus"],
            verify_live_sampler)
        _require_one_continuous_sampler(qualification_thermal, thermal)

    staged_result = _temporary_path(output_path, ".jsonl")
    staged_receipt = _temporary_path(execution_receipt_path, ".json")
    publication_identities: List[Tuple[Path, int, int]] = []
    publication_guards: List[int] = []
    committed = False
    try:
        _write_logical_jsonl(staged_result, payloads)
        try:
            if evidence_kind == "recovery":
                summary = contract_api.validate_ledger(
                    contract, phase, staged_result, freeze_path, trace_path)
            elif evidence_kind == "timing":
                summary = contract_api.validate_timing_receipt(
                    contract, phase, staged_result, freeze_path, trace_path,
                    timing_qualification_map_path,
                    timing_qualification_audit_path,
                    timing_qualification_map_sha256)
            else:
                fail("evidence kind must be recovery or timing")
        except contract_api.ContractError as exc:
            fail(str(exc))
        initial_freeze_sha256 = contract_api.freeze_manifest_sha256(freeze)
        if (summary.get("freeze_manifest_sha256") != initial_freeze_sha256 or
                summary.get("architecture_artifact_sha256") !=
                contract_api.architecture_artifact_sha256(freeze)):
            fail("freeze manifest changed during authoritative validation")
        if verify_live_sampler:
            _reverify_native_workers(runtime_workers)
        if timing_qualification is not None:
            final_qualification = _load_required_timing_qualification(
                contract, phase, timing_qualification_map_path,
                timing_qualification_audit_path,
                timing_qualification_map_sha256)
            if (final_qualification.map_sha256 !=
                    timing_qualification.map_sha256 or
                    final_qualification.qualification_audit_sha256 !=
                    timing_qualification.qualification_audit_sha256):
                fail("timing qualification changed during validation")
            final_qualification_receipt, final_qualification_metadata = \
                load_timing_qualification_execution_receipt(
                    contract, phase, final_qualification,
                    timing_qualification_native_path,
                    timing_qualification_execution_receipt_path,
                    expected_receipt_sha256=
                        timing_qualification_execution_receipt_sha256)
            if (final_qualification_receipt !=
                    qualification_execution_receipt or
                    final_qualification_metadata != qualification_metadata):
                fail("timing qualification changed during validation")
        try:
            final_freeze = contract_api.load_freeze_manifest(
                contract, phase, freeze_path, evidence_kind,
                timing_qualification)
            final_trace_sha256 = contract_api.trace_manifest_sha256(
                contract, evidence_kind, phase, trace_path,
                timing_qualification)
        except contract_api.ContractError as exc:
            fail(str(exc))
        if (final_freeze != freeze or
                contract_api.freeze_manifest_sha256(final_freeze) !=
                    initial_freeze_sha256 or
                final_trace_sha256 != freeze["trace_manifest_sha256"]):
            fail("freeze or trace manifest changed during validation")
        final_records = _load_canonical_jsonl(
            native_path, "native result stream after validation")
        if (final_records != records or
                _sha256_jsonl(final_records) != metadata["native_stream_sha256"]):
            fail("native result stream changed during validation")
        final_sampler = _load_canonical_object(
            sampler_path, "sampler attestation after validation")
        final_thermal = _thermal_window(
            final_sampler, metadata["worker_start_monotonic_ns"],
            metadata["worker_end_monotonic_ns"], metadata["worker_cpus"],
            verify_live_sampler,
            controller_cpu=freeze.get("host_identity", {}).get(
                "controller_cpu"))
        if final_sampler != sampler or final_thermal != thermal:
            fail("timing sampler evidence changed during validation")
        if qualification_metadata is not None and \
                qualification_sampler is not None:
            final_qualification_sampler = _load_canonical_object(
                timing_qualification_sampler_path,
                "qualification sampler attestation after validation")
            final_qualification_thermal = _thermal_window(
                final_qualification_sampler,
                qualification_metadata[
                    "qualification_worker_start_monotonic_ns"],
                qualification_metadata[
                    "qualification_worker_end_monotonic_ns"],
                qualification_metadata["qualification_worker_cpus"],
                verify_live_sampler)
            if (final_qualification_sampler != qualification_sampler or
                    final_qualification_thermal != qualification_thermal):
                fail("qualification sampler evidence changed during validation")
            _require_one_continuous_sampler(
                final_qualification_thermal, final_thermal)
        result_sha256 = _sha256_jsonl(payloads)
        if evidence_kind == "timing":
            receipt_schema = TIMING_EXECUTION_SCHEMA
        else:
            receipt_schema = RAW_EXECUTION_SCHEMA \
                if freeze.get("schema") == contract_api.RAW_FREEZE_SCHEMA \
                else EXECUTION_SCHEMA
        receipt: Dict[str, Any] = {
            "schema": receipt_schema,
            "contract_sha256": contract_api.contract_sha256(contract),
            "evidence_kind": evidence_kind,
            "phase": phase,
            "freeze_manifest_sha256":
                initial_freeze_sha256,
            "trace_manifest_sha256": freeze["trace_manifest_sha256"],
            "result_stream_sha256": result_sha256,
            **metadata,
            "arm_descriptor_sha256s": sorted(
                arm["arm_descriptor_sha256"] for arm in freeze["arms"]),
            "thermal": thermal,
            "validator_summary_sha256": contract_api.sha256_json(summary),
        }
        if timing_qualification is not None and \
                qualification_metadata is not None:
            receipt.update({
                "timing_base_domain_sha256":
                    timing_qualification.base_domain_sha256,
                "timing_qualified_domain_sha256":
                    timing_qualification.qualified_domain_sha256,
                "timing_qualification_map_sha256":
                    timing_qualification.map_sha256,
                "timing_qualification_audit_sha256":
                    timing_qualification.qualification_audit_sha256,
                "timing_qualification_execution_receipt_sha256":
                    timing_qualification_execution_receipt_sha256,
                **{
                    field: qualification_metadata[field]
                    for field in TIMING_QUALIFICATION_EMBEDDED_FIELDS
                },
                "qualification_thermal": qualification_thermal,
            })
        receipt["receipt_sha256"] = contract_api.sha256_json(receipt)
        _write_canonical_object(staged_receipt, receipt)
        result_guard, result_info = _open_regular_nofollow(
            staged_result, "staged result stream")
        publication_guards.append(result_guard)
        receipt_guard, receipt_info = _open_regular_nofollow(
            staged_receipt, "staged execution receipt")
        publication_guards.append(receipt_guard)
        publication_identities = [
            (output_path, result_info.st_dev, result_info.st_ino),
            (execution_receipt_path, receipt_info.st_dev, receipt_info.st_ino),
        ]
        _publish_new(
            staged_result, output_path,
            (result_info.st_dev, result_info.st_ino))
        published_payloads = _load_canonical_jsonl(
            output_path, "published result dependency")
        if (published_payloads != payloads or
                _sha256_jsonl(published_payloads) != result_sha256):
            fail("published result dependency changed before receipt commit")
        try:
            if evidence_kind == "recovery":
                published_summary = contract_api.validate_ledger(
                    contract, phase, output_path, freeze_path, trace_path)
            else:
                published_summary = contract_api.validate_timing_receipt(
                    contract, phase, output_path, freeze_path, trace_path,
                    timing_qualification_map_path,
                    timing_qualification_audit_path,
                    timing_qualification_map_sha256)
        except contract_api.ContractError as exc:
            fail(str(exc))
        if published_summary != summary:
            fail("published result dependency changed validator semantics")
        _publish_new(
            staged_receipt, execution_receipt_path,
            (receipt_info.st_dev, receipt_info.st_ino))
        committed = True
        validated = validate_execution_receipt(
            contract, evidence_kind, phase, freeze_path, trace_path,
            native_path, output_path, execution_receipt_path,
            verify_live_sampler=verify_live_sampler,
            sampler_path=sampler_path,
            timing_qualification_map_path=timing_qualification_map_path,
            timing_qualification_audit_path=timing_qualification_audit_path,
            timing_qualification_map_sha256=
                timing_qualification_map_sha256,
            timing_qualification_native_path=timing_qualification_native_path,
            timing_qualification_sampler_path=
                timing_qualification_sampler_path,
            timing_qualification_execution_receipt_path=
                timing_qualification_execution_receipt_path,
            timing_qualification_execution_receipt_sha256=
                timing_qualification_execution_receipt_sha256)
        if (validated["summary"] != summary or
                validated["execution_receipt"] != receipt):
            fail("published result pair differs after terminal validation")
        while publication_guards:
            guard = publication_guards[-1]
            try:
                os.close(guard)
            finally:
                publication_guards.pop()
        return {"summary": summary, "execution_receipt": receipt}
    except BaseException:
        if not committed:
            staged_identities = {
                staged: (device, inode)
                for staged, (_, device, inode) in zip(
                    (staged_result, staged_receipt), publication_identities)
            }
            for path in (staged_result, staged_receipt):
                identity = staged_identities.get(path)
                if identity is not None:
                    _unlink_if_identity(path, identity[0], identity[1])
                else:
                    try:
                        path.unlink()
                    except FileNotFoundError:
                        pass
        while publication_guards:
            guard = publication_guards[-1]
            try:
                os.close(guard)
            except OSError:
                pass
            finally:
                publication_guards.pop()
        raise


def validate_execution_receipt(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        freeze_path: Path, trace_path: Path, native_path: Path,
        result_path: Path, execution_receipt_path: Path,
        verify_live_sampler: bool = True,
        sampler_path: Optional[Path] = None,
        timing_qualification_map_path: Optional[Path] = None,
        timing_qualification_audit_path: Optional[Path] = None,
        timing_qualification_map_sha256: Optional[str] = None,
        timing_qualification_native_path: Optional[Path] = None,
        timing_qualification_sampler_path: Optional[Path] = None,
        timing_qualification_execution_receipt_path: Optional[Path] = None,
        timing_qualification_execution_receipt_sha256: Optional[str] = None,
        ) -> Mapping[str, Any]:
    """Revalidate a terminal receipt against every bound source artifact."""
    if sampler_path is None:
        fail("execution validation requires the sampler attestation")
    receipt_fields = TIMING_EXECUTION_FIELDS \
        if evidence_kind == "timing" else EXECUTION_FIELDS
    receipt = _exact_keys(
        _load_canonical_object(execution_receipt_path, "execution receipt"),
        receipt_fields, "execution receipt")
    timing_qualification: Optional[contract_api.TimingQualification] = None
    qualification_metadata: Optional[Mapping[str, Any]] = None
    qualification_sampler: Optional[Mapping[str, Any]] = None
    qualification_execution_receipt: Optional[Mapping[str, Any]] = None
    if evidence_kind == "timing":
        timing_qualification = _load_required_timing_qualification(
            contract, phase, timing_qualification_map_path,
            timing_qualification_audit_path,
            timing_qualification_map_sha256)
        if timing_qualification_native_path is None:
            fail("timing evidence requires the native qualification stream")
        if timing_qualification_sampler_path is None:
            fail("timing evidence requires a qualification thermal attestation")
        if (timing_qualification_execution_receipt_path is None or
                timing_qualification_execution_receipt_sha256 is None):
            fail("timing evidence requires the qualification execution receipt")
        qualification_sampler = _load_canonical_object(
            timing_qualification_sampler_path,
            "qualification sampler attestation")
        qualification_execution_receipt, qualification_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, timing_qualification,
                timing_qualification_native_path,
                timing_qualification_execution_receipt_path,
                expected_receipt_sha256=
                    timing_qualification_execution_receipt_sha256)
    try:
        freeze = contract_api.load_freeze_manifest(
            contract, phase, freeze_path, evidence_kind,
            timing_qualification)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if evidence_kind == "timing":
        expected_receipt_schema = TIMING_EXECUTION_SCHEMA
    else:
        expected_receipt_schema = RAW_EXECUTION_SCHEMA \
            if freeze.get("schema") == contract_api.RAW_FREEZE_SCHEMA \
            else EXECUTION_SCHEMA
    if receipt["schema"] != expected_receipt_schema:
        fail("execution receipt has an unknown schema")
    unsigned = {key: value for key, value in receipt.items()
                if key != "receipt_sha256"}
    if (not _is_sha256(receipt["receipt_sha256"]) or
            receipt["receipt_sha256"] != contract_api.sha256_json(unsigned)):
        fail("execution receipt self-hash is invalid")
    if (receipt["contract_sha256"] != contract_api.contract_sha256(contract) or
            receipt["evidence_kind"] != evidence_kind or
            receipt["phase"] != phase or
            receipt["freeze_manifest_sha256"] !=
                contract_api.freeze_manifest_sha256(freeze) or
            receipt["trace_manifest_sha256"] !=
                freeze["trace_manifest_sha256"]):
        fail("execution receipt contract/freeze/trace binding is invalid")
    for field in (
            "result_stream_sha256", "native_stream_sha256",
            "message_set_sha256", "work_set_sha256",
            "validator_summary_sha256"):
        if not _is_sha256(receipt[field]):
            fail("execution receipt {} is not a SHA-256".format(field))

    native_records = _load_canonical_jsonl(native_path, "native result stream")
    payloads, metadata_with_runtime = _validate_native_records(
        contract, freeze, evidence_kind, phase, native_records,
        timing_qualification,
        verify_live_workers=verify_live_sampler)
    runtime_workers = metadata_with_runtime["_runtime_workers"]
    metadata = {
        key: value for key, value in metadata_with_runtime.items()
        if key != "_runtime_workers"
    }
    for field in (
            "record_count", "worker_start_monotonic_ns",
            "worker_end_monotonic_ns", "worker_cpus", "workers",
            "worker_binary_sha256s", "message_set_sha256",
            "work_set_sha256", "native_stream_sha256"):
        if receipt[field] != metadata[field]:
            fail("execution receipt {} differs from native evidence".format(
                field))
    if timing_qualification is not None and \
            qualification_metadata is not None:
        qualification_bindings = {
            "timing_base_domain_sha256":
                timing_qualification.base_domain_sha256,
            "timing_qualified_domain_sha256":
                timing_qualification.qualified_domain_sha256,
            "timing_qualification_map_sha256":
                timing_qualification.map_sha256,
            "timing_qualification_audit_sha256":
                timing_qualification.qualification_audit_sha256,
            "timing_qualification_execution_receipt_sha256":
                timing_qualification_execution_receipt_sha256,
            **{
                field: qualification_metadata[field]
                for field in TIMING_QUALIFICATION_EMBEDDED_FIELDS
            },
        }
        for field, expected in qualification_bindings.items():
            if receipt[field] != expected:
                fail("execution receipt {} differs from qualification evidence".
                     format(field))
        if receipt["qualification_worker_end_monotonic_ns"] > \
                receipt["worker_start_monotonic_ns"]:
            fail("qualification workers overlap the exact-eight timing interval")
        qualification_identities = {
            (value["pid"], value["process_start_ticks"])
            for value in receipt["qualification_workers"]
        }
        timing_identities = {
            (value["pid"], value["process_start_ticks"])
            for value in receipt["workers"]
        }
        if qualification_identities & timing_identities:
            fail("a qualification worker identity survived into timing")
        qualification_thermal = _exact_keys(
            receipt["qualification_thermal"], THERMAL_RECEIPT_FIELDS,
            "execution receipt qualification thermal record")
        if qualification_sampler is None:
            fail("qualification sampler attestation is missing")
        actual_qualification_thermal = _thermal_window(
            qualification_sampler,
            receipt["qualification_worker_start_monotonic_ns"],
            receipt["qualification_worker_end_monotonic_ns"],
            receipt["qualification_worker_cpus"], verify_live_sampler)
        if actual_qualification_thermal != qualification_thermal:
            fail("qualification thermal summary differs from sampler bytes")
    if receipt["result_stream_sha256"] != _sha256_jsonl(payloads):
        fail("execution receipt result hash differs from native payloads")
    result_rows = _load_canonical_jsonl(result_path, "published result stream")
    if result_rows != payloads or _sha256_jsonl(result_rows) != \
            receipt["result_stream_sha256"]:
        fail("published result stream differs from the receipted native bytes")
    expected_descriptors = sorted(
        arm["arm_descriptor_sha256"] for arm in freeze["arms"])
    if receipt["arm_descriptor_sha256s"] != expected_descriptors:
        fail("execution receipt descriptor roster differs from the freeze")
    expected_binaries = sorted(set(
        arm["binary_sha256"] for arm in freeze["arms"]))
    if receipt["worker_binary_sha256s"] != expected_binaries:
        fail("execution receipt binary roster differs from the freeze")
    workers = receipt["workers"]
    if (not isinstance(workers, list) or
            any(not isinstance(value, dict) or set(value) != WORKER_FIELDS
                for value in workers)):
        fail("execution receipt worker roster is malformed")

    thermal = _exact_keys(
        receipt["thermal"], THERMAL_RECEIPT_FIELDS,
        "execution receipt thermal record")
    sampler = {
        "schema": SAMPLER_SCHEMA,
        "pid": thermal["pid"],
        "cpu": thermal["cpu"],
        "process_start_ticks": thermal["process_start_ticks"],
        "script_path": thermal["script_path"],
        "script_sha256": thermal["script_sha256"],
        "csv_path": thermal["csv_path"],
        "csv_device": thermal["csv_device"],
        "csv_inode": thermal["csv_inode"],
        "window_start_monotonic_ns": thermal[
            "window_start_monotonic_ns"],
        "window_end_monotonic_ns": thermal["window_end_monotonic_ns"],
        "terminal_status": thermal["terminal_status"],
    }
    sampler_attestation = _exact_keys(
        _load_canonical_object(sampler_path, "sampler attestation"),
        SAMPLER_FIELDS, "sampler attestation")
    if sampler_attestation != sampler:
        fail("execution receipt sampler binding differs from its attestation")
    actual_thermal = _thermal_window(
        sampler_attestation, receipt["worker_start_monotonic_ns"],
        receipt["worker_end_monotonic_ns"], receipt["worker_cpus"],
        verify_live_sampler,
        controller_cpu=freeze.get("host_identity", {}).get("controller_cpu"))
    if actual_thermal != thermal:
        fail("execution receipt thermal summary differs from sampler bytes")
    if evidence_kind == "timing":
        _require_one_continuous_sampler(
            receipt["qualification_thermal"], thermal)

    try:
        if evidence_kind == "recovery":
            summary = contract_api.validate_ledger(
                contract, phase, result_path, freeze_path, trace_path)
        else:
            summary = contract_api.validate_timing_receipt(
                contract, phase, result_path, freeze_path, trace_path,
                timing_qualification_map_path,
                timing_qualification_audit_path,
                timing_qualification_map_sha256)
    except contract_api.ContractError as exc:
        fail(str(exc))
    final_result_rows = _load_canonical_jsonl(
        result_path, "published result stream after validation")
    if (final_result_rows != result_rows or
            _sha256_jsonl(final_result_rows) !=
            receipt["result_stream_sha256"]):
        fail("published result stream changed during authoritative validation")
    if contract_api.sha256_json(summary) != receipt["validator_summary_sha256"]:
        fail("execution receipt validator-summary binding is invalid")
    if timing_qualification is not None and \
            qualification_metadata is not None:
        final_qualification = _load_required_timing_qualification(
            contract, phase, timing_qualification_map_path,
            timing_qualification_audit_path,
            timing_qualification_map_sha256)
        if (final_qualification.map_sha256 !=
                timing_qualification.map_sha256 or
                final_qualification.qualification_audit_sha256 !=
                timing_qualification.qualification_audit_sha256):
            fail("timing qualification changed during receipt validation")
        final_qualification_receipt, final_qualification_metadata = \
            load_timing_qualification_execution_receipt(
                contract, phase, final_qualification,
                timing_qualification_native_path,
                timing_qualification_execution_receipt_path,
                expected_receipt_sha256=
                    timing_qualification_execution_receipt_sha256)
        if (final_qualification_receipt != qualification_execution_receipt or
                final_qualification_metadata != qualification_metadata):
            fail("timing qualification changed during receipt validation")
    final_receipt = _exact_keys(
        _load_canonical_object(
            execution_receipt_path, "execution receipt after validation"),
        receipt_fields, "execution receipt after validation")
    if final_receipt != receipt:
        fail("execution receipt changed during authoritative validation")
    try:
        final_freeze = contract_api.load_freeze_manifest(
            contract, phase, freeze_path, evidence_kind,
            timing_qualification)
        final_trace_sha256 = contract_api.trace_manifest_sha256(
            contract, evidence_kind, phase, trace_path,
            timing_qualification)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if (final_freeze != freeze or
            contract_api.freeze_manifest_sha256(final_freeze) !=
                receipt["freeze_manifest_sha256"] or
            final_trace_sha256 != receipt["trace_manifest_sha256"]):
        fail("freeze or trace manifest changed during receipt validation")
    final_native_records = _load_canonical_jsonl(
        native_path, "native result stream after receipt validation")
    if (final_native_records != native_records or
            _sha256_jsonl(final_native_records) !=
                receipt["native_stream_sha256"]):
        fail("native result stream changed during receipt validation")
    final_sampler_attestation = _exact_keys(
        _load_canonical_object(
            sampler_path, "sampler attestation after receipt validation"),
        SAMPLER_FIELDS, "sampler attestation after receipt validation")
    final_thermal = _thermal_window(
        final_sampler_attestation, receipt["worker_start_monotonic_ns"],
        receipt["worker_end_monotonic_ns"], receipt["worker_cpus"],
        verify_live_sampler,
        controller_cpu=freeze.get("host_identity", {}).get("controller_cpu"))
    if (final_sampler_attestation != sampler_attestation or
            final_thermal != thermal):
        fail("timing sampler evidence changed during receipt validation")
    if timing_qualification is not None and \
            qualification_metadata is not None and \
            qualification_sampler is not None:
        final_qualification_sampler = _load_canonical_object(
            timing_qualification_sampler_path,
            "qualification sampler attestation after receipt validation")
        final_qualification_thermal = _thermal_window(
            final_qualification_sampler,
            qualification_metadata[
                "qualification_worker_start_monotonic_ns"],
            qualification_metadata[
                "qualification_worker_end_monotonic_ns"],
            qualification_metadata["qualification_worker_cpus"],
            verify_live_sampler)
        if (final_qualification_sampler != qualification_sampler or
                final_qualification_thermal != qualification_thermal):
            fail("qualification sampler evidence changed during receipt "
                 "validation")
        _require_one_continuous_sampler(
            final_qualification_thermal, final_thermal)
    if verify_live_sampler:
        _reverify_native_workers(runtime_workers)
    return {"summary": summary, "execution_receipt": receipt}


def _command_trace(args: argparse.Namespace) -> int:
    contract = contract_api.load_contract(args.contract)
    qualification = None
    if args.kind == "timing":
        qualification = _load_required_timing_qualification(
            contract, args.phase, args.timing_qualification_map,
            args.timing_qualification_audit,
            args.timing_qualification_map_sha256)
    frozen_hash = None
    if args.freeze_manifest is not None:
        freeze = contract_api.load_freeze_manifest(
            contract, args.phase, args.freeze_manifest, args.kind,
            qualification)
        frozen_hash = freeze["trace_manifest_sha256"]
    digest = assemble_trace_manifest(
        contract, args.kind, args.phase, args.native, args.output, frozen_hash,
        qualification)
    print(json.dumps({"trace_manifest_sha256": digest}, sort_keys=True))
    return 0


def _command_results(args: argparse.Namespace) -> int:
    contract = contract_api.load_contract(args.contract)
    result = assemble_results(
        contract, args.kind, args.phase, args.freeze_manifest,
        args.trace_manifest, args.native, args.sampler_attestation,
        args.output, args.execution_receipt, verify_live_sampler=True,
        timing_qualification_map_path=args.timing_qualification_map,
        timing_qualification_audit_path=args.timing_qualification_audit,
        timing_qualification_map_sha256=
            args.timing_qualification_map_sha256,
        timing_qualification_native_path=args.timing_qualification_native,
        timing_qualification_sampler_path=args.timing_qualification_sampler,
        timing_qualification_execution_receipt_path=
            args.timing_qualification_execution_receipt,
        timing_qualification_execution_receipt_sha256=
            args.timing_qualification_execution_receipt_sha256)
    print(json.dumps(result["summary"], sort_keys=True, indent=2))
    return 0


def _command_attest_sampler(args: argparse.Namespace) -> int:
    value = write_sampler_attestation(
        args.pid, args.cpu, args.script, args.csv,
        args.window_start_monotonic_ns, args.window_end_monotonic_ns,
        args.output)
    print(json.dumps(value, sort_keys=True, indent=2))
    return 0


def _command_validate_execution(args: argparse.Namespace) -> int:
    contract = contract_api.load_contract(args.contract)
    result = validate_execution_receipt(
        contract, args.kind, args.phase, args.freeze_manifest,
        args.trace_manifest, args.native, args.result,
        args.execution_receipt, verify_live_sampler=True,
        sampler_path=args.sampler_attestation,
        timing_qualification_map_path=args.timing_qualification_map,
        timing_qualification_audit_path=args.timing_qualification_audit,
        timing_qualification_map_sha256=
            args.timing_qualification_map_sha256,
        timing_qualification_native_path=args.timing_qualification_native,
        timing_qualification_sampler_path=args.timing_qualification_sampler,
        timing_qualification_execution_receipt_path=
            args.timing_qualification_execution_receipt,
        timing_qualification_execution_receipt_sha256=
            args.timing_qualification_execution_receipt_sha256)
    print(json.dumps(result["summary"], sort_keys=True, indent=2))
    return 0


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    commands = parser.add_subparsers(dest="command", required=True)
    sampler = commands.add_parser(
        "attest-sampler", help="bind the existing live thermal sampler")
    sampler.add_argument("--pid", type=int, required=True)
    sampler.add_argument("--cpu", type=int, required=True)
    sampler.add_argument("--script", type=Path, required=True)
    sampler.add_argument("--csv", type=Path, required=True)
    sampler.add_argument(
        "--window-start-monotonic-ns", type=int, required=True)
    sampler.add_argument(
        "--window-end-monotonic-ns", type=int, required=True)
    sampler.add_argument("--output", type=Path, required=True)
    sampler.set_defaults(function=_command_attest_sampler)
    trace = commands.add_parser(
        "assemble-traces", help="validate native trace bytes")
    trace.add_argument("--kind", choices=("recovery", "timing"), required=True)
    trace.add_argument("--phase", required=True)
    trace.add_argument("--native", required=True, type=Path)
    trace.add_argument("--output", required=True, type=Path)
    trace.add_argument("--freeze-manifest", type=Path)
    trace.add_argument("--timing-qualification-map", type=Path)
    trace.add_argument("--timing-qualification-audit", type=Path)
    trace.add_argument("--timing-qualification-map-sha256")
    trace.set_defaults(function=_command_trace)
    result = commands.add_parser(
        "assemble-results", help="validate and publish native result bytes")
    result.add_argument("--kind", choices=("recovery", "timing"), required=True)
    result.add_argument("--phase", required=True)
    result.add_argument("--freeze-manifest", required=True, type=Path)
    result.add_argument("--trace-manifest", required=True, type=Path)
    result.add_argument("--native", required=True, type=Path)
    result.add_argument("--sampler-attestation", required=True, type=Path)
    result.add_argument("--output", required=True, type=Path)
    result.add_argument("--execution-receipt", required=True, type=Path)
    result.add_argument("--timing-qualification-map", type=Path)
    result.add_argument("--timing-qualification-audit", type=Path)
    result.add_argument("--timing-qualification-map-sha256")
    result.add_argument("--timing-qualification-native", type=Path)
    result.add_argument("--timing-qualification-sampler", type=Path)
    result.add_argument(
        "--timing-qualification-execution-receipt", type=Path)
    result.add_argument(
        "--timing-qualification-execution-receipt-sha256")
    result.set_defaults(function=_command_results)
    validate = commands.add_parser(
        "validate-execution", help="revalidate a native execution receipt")
    validate.add_argument(
        "--kind", choices=("recovery", "timing"), required=True)
    validate.add_argument("--phase", required=True)
    validate.add_argument("--freeze-manifest", required=True, type=Path)
    validate.add_argument("--trace-manifest", required=True, type=Path)
    validate.add_argument("--native", required=True, type=Path)
    validate.add_argument("--result", required=True, type=Path)
    validate.add_argument("--execution-receipt", required=True, type=Path)
    validate.add_argument("--sampler-attestation", required=True, type=Path)
    validate.add_argument("--timing-qualification-map", type=Path)
    validate.add_argument("--timing-qualification-audit", type=Path)
    validate.add_argument("--timing-qualification-map-sha256")
    validate.add_argument("--timing-qualification-native", type=Path)
    validate.add_argument("--timing-qualification-sampler", type=Path)
    validate.add_argument(
        "--timing-qualification-execution-receipt", type=Path)
    validate.add_argument(
        "--timing-qualification-execution-receipt-sha256")
    validate.set_defaults(function=_command_validate_execution)
    args = parser.parse_args(argv if argv else None)
    try:
        return args.function(args)
    except (NativeEvidenceError, contract_api.ContractError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
