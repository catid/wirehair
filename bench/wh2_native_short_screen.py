#!/usr/bin/env python3
"""Fail-closed assembly for the frozen WH2 native short screen.

The native worker emits canonical JSONL envelopes.  This module checks their
cardinality and runtime provenance, strips the envelopes to the exact frozen
ledger/timing schemas, runs the authoritative contract validator, and only
then publishes the result stream and its terminal execution receipt.

It is intentionally not a campaign scheduler.  A caller may distribute work
however it likes, provided every frozen ordinal is returned exactly once.
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
TIMING_RECORD_SCHEMA = "wirehair.wh2.native-timing-record.v1"
SAMPLER_SCHEMA = "wirehair.wh2.sampler-attestation.v1"
EXECUTION_SCHEMA = "wirehair.wh2.native-execution-receipt.v1"

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
EXECUTION_FIELDS = frozenset((
    "schema", "contract_sha256", "evidence_kind", "phase",
    "freeze_manifest_sha256", "trace_manifest_sha256",
    "result_stream_sha256", "record_count", "worker_start_monotonic_ns",
    "worker_end_monotonic_ns", "worker_cpus", "workers",
    "worker_binary_sha256s", "message_set_sha256", "work_set_sha256",
    "native_stream_sha256", "arm_descriptor_sha256s", "thermal",
    "validator_summary_sha256", "receipt_sha256",
))
THERMAL_HEADER = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors", "load1",
    "load5", "load15", "edac_ce", "edac_ue",
)
MAX_THERMAL_LINE_BYTES = 65536


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
    try:
        descriptor = os.open(str(path), flags | nofollow)
    except OSError as exc:
        fail("cannot open {} {}: {}".format(context, path, exc))
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode):
            fail("{} must be a regular non-symlink file".format(context))
        return descriptor, info
    except BaseException:
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


def _publish_new(staged: Path, destination: Path) -> None:
    try:
        staged_info = os.stat(str(staged), follow_symlinks=False)
    except OSError as exc:
        fail("cannot inspect staged artifact {}: {}".format(staged, exc))
    try:
        os.link(str(staged), str(destination))
    except FileExistsError:
        fail("refusing to replace existing artifact {}".format(destination))
    except OSError as exc:
        fail("cannot publish {}: {}".format(destination, exc))
    try:
        staged.unlink()
    except FileNotFoundError:
        pass
    except BaseException:
        _unlink_if_identity(
            destination, staged_info.st_dev, staged_info.st_ino)
        raise


def _unlink_if_identity(path: Path, device: int, inode: int) -> bool:
    """Remove only our own published hard link, never an EEXIST winner."""
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


def _trace_cells(contract: Mapping[str, Any], evidence_kind: str,
                 phase: str) -> Tuple[List[Mapping[str, Any]], int]:
    if evidence_kind == "recovery":
        cells = list(contract_api.iter_recovery_cells(contract, phase))
        count = contract["recovery"]["domains"][phase][
            "expected_cells_per_arm"]
    elif evidence_kind == "timing":
        cells = list(contract_api.iter_timing_cells(contract, phase))
        count = contract["timing"]["domains"][phase]["expected_cells"]
    else:
        fail("evidence kind must be recovery or timing")
    if len(cells) != count:
        fail("contract cell iterator/cardinality mismatch")
    return cells, count


def assemble_trace_manifest(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        native_path: Path, output_path: Path,
        frozen_sha256: Optional[str] = None) -> str:
    """Validate native trace records and atomically publish a trace manifest."""
    cells, count = _trace_cells(contract, evidence_kind, phase)
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
        if type(packet_count) is not int or packet_count != K + 4:
            fail("native trace packet count must be exactly K+4")
        candidate_limit = 256 * (K + 4) + 65536
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
    try:
        _write_logical_jsonl(staged, manifest)
        digest = contract_api.trace_manifest_sha256(
            contract, evidence_kind, phase, staged)
        if frozen_sha256 is not None and digest != frozen_sha256:
            fail("native trace manifest differs from the frozen hash")
        _publish_new(staged, output_path)
        return digest
    finally:
        try:
            staged.unlink()
        except FileNotFoundError:
            pass


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


def _thermal_window(
        sampler: Mapping[str, Any], worker_start_ns: int, worker_end_ns: int,
        worker_cpus: Sequence[int], verify_live_sampler: bool) \
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
    if sampler["cpu"] in worker_cpus:
        fail("sampler CPU overlaps a frozen worker CPU")

    script_path = Path(sampler["script_path"])
    csv_path = Path(sampler["csv_path"])
    script_fd = -1
    csv_fd = -1
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

        with os.fdopen(os.dup(csv_fd), "rb") as source:
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
    try:
        _write_canonical_object(staged, value)
        _publish_new(staged, output_path)
    finally:
        try:
            staged.unlink()
        except FileNotFoundError:
            pass
    return value


def _expected_record_ordinal(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        evidence_kind: str, phase: str, payload: Mapping[str, Any]) \
        -> Tuple[int, str]:
    if evidence_kind == "recovery":
        indexes = contract_api._domain_indexes(contract, phase)
        cell_ordinal, _ = contract_api._cell_ordinal(
            contract, phase, payload, indexes)
        arm = payload.get("arm")
        if arm not in freeze["arm_roster"]:
            fail("native recovery payload arm is outside the frozen roster")
        arm_index = freeze["arm_roster"].index(arm)
        return cell_ordinal * len(freeze["arm_roster"]) + arm_index, arm
    cells = contract_api._timing_cell_indexes(contract, phase)
    cell_ordinal, _ = contract_api._timing_cell_ordinal(
        contract, phase, payload, cells)
    panel = {
        "panel_kind": payload.get("panel_kind"),
        "scope": payload.get("scope"),
        "left_arm": payload.get("left_arm"),
        "right_arm": payload.get("right_arm"),
    }
    panels = contract_api.timing_panels(contract, freeze["arm_roster"])
    panel_keys = [contract_api.canonical_json(value) for value in panels]
    try:
        panel_index = panel_keys.index(contract_api.canonical_json(panel))
    except ValueError:
        fail("native timing payload has an undeclared or relabeled panel")
    return cell_ordinal * len(panels) + panel_index, str(cell_ordinal)


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


def _reverify_native_workers(workers: Sequence[Mapping[str, Any]]) -> None:
    for worker in workers:
        _verify_live_worker_process(
            worker["pid"], worker["cpu"], worker["process_start_ticks"],
            worker["binary_sha256"])


def _validate_native_records(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        evidence_kind: str, phase: str,
        records: Sequence[Mapping[str, Any]],
        verify_live_workers: bool = False) \
        -> Tuple[List[Mapping[str, Any]], Mapping[str, Any]]:
    if evidence_kind == "recovery":
        payload_fields = contract_api.LEDGER_FIELDS
        record_schema = RECOVERY_RECORD_SCHEMA
        cell_count = contract["recovery"]["domains"][phase][
            "expected_cells_per_arm"]
        expected_count = cell_count * len(freeze["arm_roster"])
    elif evidence_kind == "timing":
        payload_fields = contract_api.TIMING_RECEIPT_FIELDS
        record_schema = TIMING_RECORD_SCHEMA
        cell_count = contract["timing"]["domains"][phase]["expected_cells"]
        expected_count = cell_count * len(contract_api.timing_panels(
            contract, freeze["arm_roster"]))
    else:
        fail("evidence kind must be recovery or timing")
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
        if row["schema"] != record_schema:
            fail("native result record has an unknown schema")
        payload = row["payload"]
        _exact_keys(payload, payload_fields, "native result payload")
        expected_ordinal, message_key = _expected_record_ordinal(
            contract, freeze, evidence_kind, phase, payload)
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
        "native_stream_sha256": _sha256_jsonl(envelopes),
        "_runtime_workers": runtime_workers,
    }
    return payloads, metadata


def assemble_results(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        freeze_path: Path, trace_path: Path, native_path: Path,
        sampler_path: Path, output_path: Path, execution_receipt_path: Path,
        verify_live_sampler: bool = True) -> Mapping[str, Any]:
    """Validate and publish one complete recovery ledger or timing receipt."""
    if output_path.exists() or execution_receipt_path.exists():
        fail("refusing to replace an existing result or execution receipt")
    try:
        freeze = contract_api.load_freeze_manifest(
            contract, phase, freeze_path, evidence_kind)
    except contract_api.ContractError as exc:
        fail(str(exc))
    records = _load_canonical_jsonl(native_path, "native result stream")
    payloads, metadata_with_runtime = _validate_native_records(
        contract, freeze, evidence_kind, phase, records,
        verify_live_workers=verify_live_sampler)
    runtime_workers = metadata_with_runtime["_runtime_workers"]
    metadata = {
        key: value for key, value in metadata_with_runtime.items()
        if key != "_runtime_workers"
    }
    sampler = _load_canonical_object(sampler_path, "sampler attestation")
    thermal = _thermal_window(
        sampler, metadata["worker_start_monotonic_ns"],
        metadata["worker_end_monotonic_ns"], metadata["worker_cpus"],
        verify_live_sampler)

    staged_result = _temporary_path(output_path, ".jsonl")
    staged_receipt = _temporary_path(execution_receipt_path, ".json")
    publication_identities: List[Tuple[Path, int, int]] = []
    try:
        _write_logical_jsonl(staged_result, payloads)
        try:
            if evidence_kind == "recovery":
                summary = contract_api.validate_ledger(
                    contract, phase, staged_result, freeze_path, trace_path)
            elif evidence_kind == "timing":
                summary = contract_api.validate_timing_receipt(
                    contract, phase, staged_result, freeze_path, trace_path)
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
        result_sha256 = _sha256_jsonl(payloads)
        receipt: Dict[str, Any] = {
            "schema": EXECUTION_SCHEMA,
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
        receipt["receipt_sha256"] = contract_api.sha256_json(receipt)
        _write_canonical_object(staged_receipt, receipt)
        result_info = os.stat(str(staged_result), follow_symlinks=False)
        receipt_info = os.stat(str(staged_receipt), follow_symlinks=False)
        publication_identities = [
            (output_path, result_info.st_dev, result_info.st_ino),
            (execution_receipt_path, receipt_info.st_dev, receipt_info.st_ino),
        ]
        _publish_new(staged_result, output_path)
        _publish_new(staged_receipt, execution_receipt_path)
        return {"summary": summary, "execution_receipt": receipt}
    except BaseException:
        for path, device, inode in publication_identities:
            _unlink_if_identity(path, device, inode)
        raise
    finally:
        for path in (staged_result, staged_receipt):
            try:
                path.unlink()
            except FileNotFoundError:
                pass


def validate_execution_receipt(
        contract: Mapping[str, Any], evidence_kind: str, phase: str,
        freeze_path: Path, trace_path: Path, native_path: Path,
        result_path: Path, execution_receipt_path: Path,
        verify_live_sampler: bool = True) -> Mapping[str, Any]:
    """Revalidate a terminal receipt against every bound source artifact."""
    receipt = _exact_keys(
        _load_canonical_object(execution_receipt_path, "execution receipt"),
        EXECUTION_FIELDS, "execution receipt")
    if receipt["schema"] != EXECUTION_SCHEMA:
        fail("execution receipt has an unknown schema")
    unsigned = {key: value for key, value in receipt.items()
                if key != "receipt_sha256"}
    if (not _is_sha256(receipt["receipt_sha256"]) or
            receipt["receipt_sha256"] != contract_api.sha256_json(unsigned)):
        fail("execution receipt self-hash is invalid")
    try:
        freeze = contract_api.load_freeze_manifest(
            contract, phase, freeze_path, evidence_kind)
    except contract_api.ContractError as exc:
        fail(str(exc))
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
    actual_thermal = _thermal_window(
        sampler, receipt["worker_start_monotonic_ns"],
        receipt["worker_end_monotonic_ns"], receipt["worker_cpus"],
        verify_live_sampler)
    if actual_thermal != thermal:
        fail("execution receipt thermal summary differs from sampler bytes")

    try:
        if evidence_kind == "recovery":
            summary = contract_api.validate_ledger(
                contract, phase, result_path, freeze_path, trace_path)
        else:
            summary = contract_api.validate_timing_receipt(
                contract, phase, result_path, freeze_path, trace_path)
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
    if verify_live_sampler:
        _reverify_native_workers(runtime_workers)
    return {"summary": summary, "execution_receipt": receipt}


def _command_trace(args: argparse.Namespace) -> int:
    contract = contract_api.load_contract(args.contract)
    frozen_hash = None
    if args.freeze_manifest is not None:
        freeze = contract_api.load_freeze_manifest(
            contract, args.phase, args.freeze_manifest, args.kind)
        frozen_hash = freeze["trace_manifest_sha256"]
    digest = assemble_trace_manifest(
        contract, args.kind, args.phase, args.native, args.output, frozen_hash)
    print(json.dumps({"trace_manifest_sha256": digest}, sort_keys=True))
    return 0


def _command_results(args: argparse.Namespace) -> int:
    contract = contract_api.load_contract(args.contract)
    result = assemble_results(
        contract, args.kind, args.phase, args.freeze_manifest,
        args.trace_manifest, args.native, args.sampler_attestation,
        args.output, args.execution_receipt, verify_live_sampler=True)
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
        args.execution_receipt, verify_live_sampler=True)
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
    validate.set_defaults(function=_command_validate_execution)
    args = parser.parse_args(argv if argv else None)
    try:
        return args.function(args)
    except (NativeEvidenceError, contract_api.ContractError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
