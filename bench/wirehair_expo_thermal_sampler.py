#!/usr/bin/env python3
"""Low-overhead post-EXPO CPU, DDR5, utilization, and EDAC sampler."""

import argparse
import csv
import ctypes
import fcntl
import glob
import hashlib
import json
import math
import os
import re
import signal
import stat
import sys
import time
from datetime import datetime, timezone


I2C_RDWR = 0x0707
I2C_M_RD = 0x0001
DIMMS = [(1, address) for address in range(0x50, 0x54)] + [
    (2, address) for address in range(0x50, 0x54)
]
# These are safety-contract constants, not command-line tuning knobs.  Every
# terminal receipt binds their exact values.
SAMPLER_SCHEMA = "wirehair.wh2.thermal_sampler.v2"
SAMPLE_VALIDATION_SCHEMA = "wirehair.wh2.thermal_validation_sample.v1"
MIN_PLAUSIBLE_DIMM_C = 0.0
MAX_PLAUSIBLE_DIMM_C = 130.0
DIMM_SAFETY_C = 90.0
MAX_DIMM_JUMP_C = 12.0
MAX_DIMM_RATE_C_PER_S = 6.0
HOT_CONFIRM_SAMPLES = 3
TELEMETRY_FAULT_ABORT_SAMPLES = 8
TERMINAL_INTENT_MARKER = b"!"

EXIT_OK = 0
EXIT_THERMAL_ABORT = 3
EXIT_TELEMETRY_ABORT = 4
EXIT_SAMPLER_ERROR = 5
EXIT_EDAC_ABORT = 6
MAX_UID = (1 << 32) - 2
MAX_TERMINAL_ERRORS = 64
MAX_TERMINAL_ERROR_TEXT = 4000


def bounded_terminal_error_text(value, fallback):
    try:
        text = str(value)
    except BaseException:
        text = fallback
    if not text:
        text = fallback
    return text[:MAX_TERMINAL_ERROR_TEXT]


def append_terminal_error(errors, phase, exc):
    if len(errors) >= MAX_TERMINAL_ERRORS:
        return
    error_type = bounded_terminal_error_text(
        getattr(type(exc), "__name__", "BaseException"), "BaseException")
    errors.append({
        "message": bounded_terminal_error_text(exc, error_type),
        "phase": bounded_terminal_error_text(phase, "unknown"),
        "type": error_type,
    })


class I2CMessage(ctypes.Structure):
    _fields_ = [
        ("addr", ctypes.c_uint16),
        ("flags", ctypes.c_uint16),
        ("len", ctypes.c_uint16),
        ("buf", ctypes.POINTER(ctypes.c_uint8)),
    ]


class I2CTransfer(ctypes.Structure):
    _fields_ = [
        ("msgs", ctypes.POINTER(I2CMessage)),
        ("nmsgs", ctypes.c_uint32),
    ]


def read_text(path):
    try:
        with open(path, "r", encoding="ascii") as stream:
            return stream.read().strip()
    except (OSError, ValueError):
        return None


def find_tctl_path():
    for name_path in glob.glob("/sys/class/hwmon/hwmon*/name"):
        if read_text(name_path) != "k10temp":
            continue
        directory = os.path.dirname(name_path)
        for label_path in glob.glob(os.path.join(directory, "temp*_label")):
            if read_text(label_path) == "Tctl":
                return label_path[:-len("_label")] + "_input"
        fallback = os.path.join(directory, "temp1_input")
        if os.path.exists(fallback):
            return fallback
    return None


def read_cpu_stat():
    raw = read_text("/proc/stat")
    if raw is None or not raw.splitlines():
        raise RuntimeError("/proc/stat is unavailable")
    tokens = raw.splitlines()[0].split()
    if len(tokens) < 9 or tokens[0] != "cpu":
        raise RuntimeError("/proc/stat aggregate CPU row is malformed")
    try:
        ticks = [int(value) for value in tokens[1:]]
    except ValueError as exc:
        raise RuntimeError("/proc/stat aggregate CPU row is malformed") from exc
    if any(value < 0 for value in ticks):
        raise RuntimeError("/proc/stat aggregate CPU counter is negative")
    idle = ticks[3] + (ticks[4] if len(ticks) > 4 else 0)
    # Linux guest and guest_nice (indices 8/9) are already included in user
    # and nice.  Summing them again biases utilization on virtualized hosts.
    return sum(ticks[:8]), idle


def cpu_busy_percent(previous, current):
    total_delta = current[0] - previous[0]
    idle_delta = current[1] - previous[1]
    if total_delta <= 0:
        return ""
    return 100.0 * (total_delta - idle_delta) / total_delta


def average_cpu_mhz():
    values = []
    for path in glob.glob("/sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq"):
        value = read_text(path)
        if value is not None:
            try:
                values.append(float(value) / 1000.0)
            except ValueError:
                pass
    return sum(values) / len(values) if values else ""


def discover_edac_paths(counter):
    paths = tuple(sorted(glob.glob(
        f"/sys/devices/system/edac/mc/mc*/{counter}"
    )))
    if not paths:
        raise RuntimeError(f"EDAC {counter} inventory is empty")
    return paths


def sum_edac(counter, expected_paths):
    paths = discover_edac_paths(counter)
    if paths != tuple(expected_paths):
        raise RuntimeError(f"EDAC {counter} inventory changed")
    total = 0
    for path in paths:
        value = read_text(path)
        if value is None:
            raise RuntimeError(f"EDAC counter became unreadable: {path}")
        if not value.isascii() or not value.isdecimal() or \
                (len(value) > 1 and value.startswith("0")):
            raise RuntimeError(f"EDAC counter is malformed: {path}")
        try:
            count = int(value)
        except ValueError as exc:
            raise RuntimeError(f"EDAC counter is malformed: {path}") from exc
        if count > (1 << 64) - 1:
            raise RuntimeError(f"EDAC counter exceeds uint64: {path}")
        total += count
        if total > (1 << 64) - 1:
            raise RuntimeError("aggregate EDAC counter exceeds uint64")
    return total


def read_spd5118_temperature(fd, address):
    register = (ctypes.c_uint8 * 1)(0x31)
    response = (ctypes.c_uint8 * 2)()
    messages = (I2CMessage * 2)(
        I2CMessage(address, 0, 1, register),
        I2CMessage(address, I2C_M_RD, 2, response),
    )
    transfer = I2CTransfer(messages, 2)
    fcntl.ioctl(fd, I2C_RDWR, transfer)
    raw = response[0] | (response[1] << 8)
    quarter_degrees = (raw >> 2) & 0x7FF
    if quarter_degrees & 0x400:
        quarter_degrees -= 0x800
    return quarter_degrees / 4.0


def read_dimm_temperatures(bus_fds, attempts, retry_delay):
    """Read every DIMM while preserving every final decoded temperature.

    Only failed I2C transactions are retried.  A successfully decoded value is
    raw evidence even when it is physically implausible; plausibility is
    accounted for separately and must never silently clamp or drop that value.
    """
    temperatures = {}
    attempt_errors = {sensor: 0 for sensor in DIMMS}
    pending = list(DIMMS)
    for attempt in range(attempts):
        failed = []
        for bus, address in pending:
            try:
                temperatures[(bus, address)] = read_spd5118_temperature(
                    bus_fds[bus], address)
            except OSError:
                attempt_errors[(bus, address)] += 1
                failed.append((bus, address))
        pending = failed
        if not pending:
            break
        if attempt + 1 < attempts:
            time.sleep(retry_delay)
    return temperatures, pending, attempt_errors


def dimm_name(sensor):
    bus, address = sensor
    return f"dimm_i2c{bus}_{address:02x}_c"


def threshold_metadata():
    return {
        "hot_confirm_samples": HOT_CONFIRM_SAMPLES,
        "max_dimm_jump_c": MAX_DIMM_JUMP_C,
        "max_dimm_rate_c_per_s": MAX_DIMM_RATE_C_PER_S,
        "max_plausible_dimm_c_exclusive": MAX_PLAUSIBLE_DIMM_C,
        "min_plausible_dimm_c_exclusive": MIN_PLAUSIBLE_DIMM_C,
        "telemetry_fault_abort_samples": TELEMETRY_FAULT_ABORT_SAMPLES,
        "dimm_safety_c_inclusive": DIMM_SAFETY_C,
    }


def canonical_json_bytes(value):
    return (json.dumps(
        value, allow_nan=False, separators=(",", ":"), sort_keys=True
    ) + "\n").encode("ascii")


class DimmPlausibilityMonitor:
    """Deterministic per-sensor plausibility and abort state machine."""

    def __init__(self, edac_ce_baseline=0, edac_ue_baseline=0):
        self._check_nonnegative_integer(edac_ce_baseline, "EDAC CE baseline")
        self._check_nonnegative_integer(edac_ue_baseline, "EDAC UE baseline")
        self.edac_ce_baseline = edac_ce_baseline
        self.edac_ue_baseline = edac_ue_baseline
        self.edac_ce_last = edac_ce_baseline
        self.edac_ue_last = edac_ue_baseline
        self.last_sample_time = None
        self.sample_count = 0
        self.consecutive_fault_rows = 0
        self.max_consecutive_fault_rows = 0
        self.decision = "continue"
        self.states = {
            sensor: {
                "attempt_errors": 0,
                "hot_candidate_temp": None,
                "hot_candidate_time": None,
                "hot_streak": 0,
                "invalid_samples": 0,
                "last_valid_temp": None,
                "last_valid_time": None,
                "max_hot_streak": 0,
                "max_raw_c": None,
                "max_valid_c": None,
                "raw_samples": 0,
                "read_error_samples": 0,
                "valid_samples": 0,
            }
            for sensor in DIMMS
        }

    @staticmethod
    def _check_nonnegative_integer(value, label):
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{label} must be a nonnegative integer")

    @staticmethod
    def _coherent_change(previous_temp, previous_time, temp, sample_time):
        elapsed = sample_time - previous_time
        if elapsed <= 0.0:
            return False
        jump = abs(temp - previous_temp)
        return jump <= MAX_DIMM_JUMP_C and \
            jump / elapsed <= MAX_DIMM_RATE_C_PER_S

    def observe(self, sample_time, temperatures, attempt_errors,
                edac_ce, edac_ue):
        if self.decision != "continue":
            raise RuntimeError("plausibility monitor is already terminal")
        if isinstance(sample_time, bool) or \
                not isinstance(sample_time, (int, float)) or \
                not math.isfinite(sample_time) or sample_time < 0.0:
            raise ValueError("sample monotonic time must be finite and nonnegative")
        sample_time = float(sample_time)
        if self.last_sample_time is not None and \
                sample_time <= self.last_sample_time:
            raise ValueError("sample monotonic time must increase")
        if set(temperatures) - set(DIMMS):
            raise ValueError("temperature map contains an unknown DIMM")
        if set(attempt_errors) - set(DIMMS):
            raise ValueError("attempt-error map contains an unknown DIMM")
        self._check_nonnegative_integer(edac_ce, "EDAC CE")
        self._check_nonnegative_integer(edac_ue, "EDAC UE")
        if edac_ce < self.edac_ce_last or edac_ue < self.edac_ue_last:
            raise ValueError("EDAC counters must not decrease")

        # Validate and normalize the entire row before mutating any state.
        # A malformed final sensor must not leave a partially counted sample.
        normalized = {}
        for sensor in DIMMS:
            attempts = attempt_errors.get(sensor, 0)
            self._check_nonnegative_integer(
                attempts, f"attempt errors for {dimm_name(sensor)}")
            temperature = temperatures.get(sensor)
            if temperature is not None:
                if isinstance(temperature, bool) or \
                        not isinstance(temperature, (int, float)) or \
                        not math.isfinite(temperature):
                    raise ValueError(
                        f"raw temperature for {dimm_name(sensor)} must be finite")
                temperature = float(temperature)
            normalized[sensor] = (temperature, attempts)

        self.sample_count += 1
        self.last_sample_time = sample_time
        self.edac_ce_last = edac_ce
        self.edac_ue_last = edac_ue
        sensors = {}
        fault_count = 0
        read_error_count = 0
        hot_sensors = []

        for sensor in DIMMS:
            state = self.states[sensor]
            temperature, attempts = normalized[sensor]
            state["attempt_errors"] += attempts

            if temperature is None:
                state["read_error_samples"] += 1
                state["hot_candidate_temp"] = None
                state["hot_candidate_time"] = None
                state["hot_streak"] = 0
                fault_count += 1
                read_error_count += 1
                sensors[dimm_name(sensor)] = {
                    "attempt_errors": attempts,
                    "hot": False,
                    "hot_streak": 0,
                    "jump_c": None,
                    "rate_c_per_s": None,
                    "raw_c": None,
                    "reason": "read_error",
                    "valid": False,
                }
                continue

            state["raw_samples"] += 1
            if state["max_raw_c"] is None or temperature > state["max_raw_c"]:
                state["max_raw_c"] = temperature

            absolute_valid = (
                MIN_PLAUSIBLE_DIMM_C < temperature < MAX_PLAUSIBLE_DIMM_C)
            # A value beyond the decode-plausibility range is still a safety
            # candidate when it is high.  One such value remains an invalid
            # spike; three coherent high values must not be dismissed for
            # eight rows merely because they also exceed the plausibility cap.
            hot = temperature >= DIMM_SAFETY_C
            if hot:
                if state["hot_candidate_temp"] is not None and \
                        self._coherent_change(
                            state["hot_candidate_temp"],
                            state["hot_candidate_time"],
                            temperature,
                            sample_time):
                    state["hot_streak"] += 1
                else:
                    state["hot_streak"] = 1
                state["hot_candidate_temp"] = temperature
                state["hot_candidate_time"] = sample_time
                state["max_hot_streak"] = max(
                    state["max_hot_streak"], state["hot_streak"])
                if state["hot_streak"] >= HOT_CONFIRM_SAMPLES:
                    hot_sensors.append(dimm_name(sensor))
            else:
                state["hot_candidate_temp"] = None
                state["hot_candidate_time"] = None
                state["hot_streak"] = 0

            jump = None
            rate = None
            reasons = []
            if not absolute_valid:
                reasons.append("absolute_range")
            elif state["last_valid_temp"] is not None:
                elapsed = sample_time - state["last_valid_time"]
                jump = abs(temperature - state["last_valid_temp"])
                rate = jump / elapsed
                if jump > MAX_DIMM_JUMP_C:
                    reasons.append("jump")
                if rate > MAX_DIMM_RATE_C_PER_S:
                    reasons.append("rate")

            valid = not reasons
            if valid:
                state["valid_samples"] += 1
                state["last_valid_temp"] = temperature
                state["last_valid_time"] = sample_time
                if state["max_valid_c"] is None or \
                        temperature > state["max_valid_c"]:
                    state["max_valid_c"] = temperature
                reason = "ok"
            else:
                state["invalid_samples"] += 1
                fault_count += 1
                reason = "+".join(reasons)

            sensors[dimm_name(sensor)] = {
                "attempt_errors": attempts,
                "hot": hot,
                "hot_streak": state["hot_streak"],
                "jump_c": jump,
                "rate_c_per_s": rate,
                "raw_c": temperature,
                "reason": reason,
                "valid": valid,
            }

        if fault_count:
            self.consecutive_fault_rows += 1
            self.max_consecutive_fault_rows = max(
                self.max_consecutive_fault_rows,
                self.consecutive_fault_rows)
        else:
            self.consecutive_fault_rows = 0

        edac_ce_delta = edac_ce - self.edac_ce_baseline
        edac_ue_delta = edac_ue - self.edac_ue_baseline
        if hot_sensors:
            self.decision = "thermal_abort"
        elif edac_ce_delta or edac_ue_delta:
            self.decision = "edac_abort"
        elif self.consecutive_fault_rows >= TELEMETRY_FAULT_ABORT_SAMPLES:
            self.decision = "telemetry_abort"
        else:
            self.decision = "continue"

        return {
            "consecutive_fault_rows": self.consecutive_fault_rows,
            "decision": self.decision,
            "edac_ce_delta": edac_ce_delta,
            "edac_ue_delta": edac_ue_delta,
            "fault_count": fault_count,
            "hot_sensors": hot_sensors,
            "monotonic_s": sample_time,
            "read_error_count": read_error_count,
            "sample_index": self.sample_count - 1,
            "schema": SAMPLE_VALIDATION_SCHEMA,
            "sensors": sensors,
        }

    def summary(self):
        sensors = {}
        for sensor in DIMMS:
            state = self.states[sensor]
            sensors[dimm_name(sensor)] = {
                key: state[key]
                for key in (
                    "attempt_errors", "invalid_samples", "max_hot_streak",
                    "max_raw_c", "max_valid_c", "raw_samples",
                    "read_error_samples", "valid_samples")
            }
        return {
            "consecutive_fault_rows": self.consecutive_fault_rows,
            "decision": self.decision,
            "dimm_attempt_errors_total": sum(
                state["attempt_errors"] for state in self.states.values()),
            "dimm_invalid_samples_total": sum(
                state["invalid_samples"] for state in self.states.values()),
            "dimm_read_error_samples_total": sum(
                state["read_error_samples"] for state in self.states.values()),
            "edac_ce_baseline": self.edac_ce_baseline,
            "edac_ce_delta": self.edac_ce_last - self.edac_ce_baseline,
            "edac_ce_last": self.edac_ce_last,
            "edac_ue_baseline": self.edac_ue_baseline,
            "edac_ue_delta": self.edac_ue_last - self.edac_ue_baseline,
            "edac_ue_last": self.edac_ue_last,
            "max_consecutive_fault_rows": self.max_consecutive_fault_rows,
            "sample_count": self.sample_count,
            "sensors": sensors,
        }


def utc_timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def binding_from_fd(fd):
    """Hash and describe one already-open regular-file inode."""
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode):
        raise RuntimeError("bound evidence descriptor is not a regular file")
    digest = hashlib.sha256()
    offset = 0
    while True:
        chunk = os.pread(fd, 1024 * 1024, offset)
        if not chunk:
            break
        digest.update(chunk)
        offset += len(chunk)
    after = os.fstat(fd)
    if (before.st_dev, before.st_ino, before.st_mode, before.st_nlink,
            before.st_uid, before.st_gid, before.st_size,
            before.st_mtime_ns, before.st_ctime_ns) != \
            (after.st_dev, after.st_ino, after.st_mode, after.st_nlink,
             after.st_uid, after.st_gid, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns) or \
            offset != after.st_size:
        raise RuntimeError("bound evidence metadata changed while hashing")
    return {
        "device": after.st_dev,
        "gid": after.st_gid,
        "inode": after.st_ino,
        "mode": after.st_mode & 0o7777,
        "nlink": after.st_nlink,
        "sha256": digest.hexdigest(),
        "size": after.st_size,
        "uid": after.st_uid,
    }


def _stat_metadata(info):
    return {
        "device": info.st_dev,
        "gid": info.st_gid,
        "inode": info.st_ino,
        "mode": info.st_mode & 0o7777,
        "nlink": info.st_nlink,
        "size": info.st_size,
        "uid": info.st_uid,
    }


def directory_binding_from_fd(fd):
    """Describe one held directory inode without mutable timestamps."""
    info = os.fstat(fd)
    if not stat.S_ISDIR(info.st_mode):
        raise RuntimeError("bound evidence parent is not a directory")
    binding = _stat_metadata(info)
    del binding["size"]
    return binding


def validate_secure_parent(destination):
    """Revalidate the held parent and the pathname which originally named it.

    The parent is the trust boundary for name creation and deletion.  A caller
    with the externally trusted output-owner UID can still defeat POSIX name
    checks, so output parents must be private to that UID rather than shared
    writable space.  Source parents are separately content/path-bound and do
    not use the output-owner contract.
    """
    observed = directory_binding_from_fd(destination["dir_fd"])
    if observed != destination["parent_binding"]:
        raise RuntimeError(
            f"evidence parent binding changed: {destination['parent_path']}")
    expected_owner_uid = destination.get("expected_owner_uid")
    if expected_owner_uid is not None and \
            observed["uid"] != expected_owner_uid:
        raise PermissionError(
            "evidence parent owner does not match the externally trusted UID: "
            f"{destination['parent_path']}")
    if observed["mode"] & 0o022:
        raise PermissionError(
            "evidence parent must not be group/other-writable: "
            f"{destination['parent_path']}")
    try:
        path_info = os.stat(
            destination["parent_path"], follow_symlinks=False)
    except OSError as exc:
        raise RuntimeError(
            f"evidence parent path became unavailable: "
            f"{destination['parent_path']}") from exc
    path_binding = directory_binding_from_stat(path_info)
    if path_binding != observed:
        raise RuntimeError(
            f"evidence parent path binding changed: "
            f"{destination['parent_path']}")


def directory_binding_from_stat(info):
    if not stat.S_ISDIR(info.st_mode):
        raise RuntimeError("evidence parent path is not a directory")
    binding = _stat_metadata(info)
    del binding["size"]
    return binding


def prepare_evidence_destination(path, *, require_absent=True,
                                 expected_owner_uid=None):
    """Open and bind the trusted parent used for all later name operations."""
    if expected_owner_uid is not None and (
            isinstance(expected_owner_uid, bool) or
            not isinstance(expected_owner_uid, int) or
            not 0 <= expected_owner_uid <= MAX_UID):
        raise ValueError("expected output-owner UID is outside its domain")
    canonical = canonical_destination(path)
    parent_path = os.path.dirname(canonical)
    basename = os.path.basename(canonical)
    if not basename or basename in (".", "..") or "/" in basename:
        raise ValueError("evidence basename must be one path component")
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | \
        getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    dir_fd = os.open(parent_path, flags)
    destination = {
        "basename": basename,
        "dir_fd": dir_fd,
        "expected_owner_uid": expected_owner_uid,
        "parent_binding": None,
        "parent_path": parent_path,
        "path": canonical,
    }
    try:
        destination["parent_binding"] = directory_binding_from_fd(dir_fd)
        validate_secure_parent(destination)
        if require_absent:
            try:
                os.stat(basename, dir_fd=dir_fd, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                raise FileExistsError(
                    f"evidence destination already exists: {canonical}")
        validate_secure_parent(destination)
        return destination
    except BaseException:
        try:
            os.close(dir_fd)
        except OSError:
            pass
        raise


def close_evidence_destination(destination):
    if destination is None or destination.get("dir_fd", -1) < 0:
        return
    fd = destination["dir_fd"]
    destination["dir_fd"] = -1
    os.close(fd)


def evidence_destination_metadata(destination):
    if destination is None:
        return None
    return {
        "basename": destination["basename"],
        "expected_owner_uid": destination.get("expected_owner_uid"),
        "parent": {
            "binding": destination["parent_binding"],
            "path": destination["parent_path"],
        },
        "path": destination["path"],
    }


def seal_and_hash_output(output):
    """Seal and hash the exact open inode, immune to pathname replacement."""
    output.flush()
    os.fsync(output.fileno())
    os.fchmod(output.fileno(), 0o444)
    os.fsync(output.fileno())
    binding = binding_from_fd(output.fileno())
    if binding["uid"] != os.geteuid() or binding["nlink"] != 1:
        raise RuntimeError(
            "sealed evidence must be owned by the sampler UID and single-link")
    return binding


def validate_sampling_arguments(interval, dimm_attempts, dimm_retry_delay):
    if not math.isfinite(interval) or interval <= 0.0:
        raise ValueError("--interval must be finite and positive")
    if dimm_attempts <= 0:
        raise ValueError("--dimm-attempts must be positive")
    if not math.isfinite(dimm_retry_delay) or dimm_retry_delay < 0.0:
        raise ValueError("--dimm-retry-delay must be finite and nonnegative")


def open_exclusive_evidence(path, *, newline=None, destination=None,
                            dir_fd=None, basename=None, live_mode=0o444):
    # Request no write bits at creation, so even the open-to-fchmod interval is
    # safe under umask 000.  The creating root descriptor remains writable.
    # O_RDWR permits hashing the bound inode through the descriptor at seal.
    if isinstance(live_mode, bool) or live_mode not in (0o444, 0o600):
        raise ValueError("live evidence mode must be 0444 or 0600")
    owned_destination = None
    if destination is not None:
        if dir_fd is not None or basename is not None:
            raise ValueError("pass either destination or dirfd fields")
        validate_secure_parent(destination)
        dir_fd = destination["dir_fd"]
        basename = destination["basename"]
    elif dir_fd is None:
        owned_destination = prepare_evidence_destination(path)
        destination = owned_destination
        dir_fd = destination["dir_fd"]
        basename = destination["basename"]
    target = basename if dir_fd is not None else path
    if dir_fd is not None and (not basename or "/" in basename):
        raise ValueError("dirfd evidence basename must be one path component")
    fd = -1
    try:
        fd = os.open(
            target, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o444, dir_fd=dir_fd)
        os.fchmod(fd, live_mode)
        if dir_fd is not None:
            os.fsync(dir_fd)
        if destination is not None:
            validate_secure_parent(destination)
        stream = os.fdopen(fd, "w+", newline=newline, encoding="ascii")
        fd = -1  # Ownership transferred to stream.
        return stream
    except BaseException:
        if fd >= 0:
            try:
                os.close(fd)
            except OSError:
                pass
        raise
    finally:
        if owned_destination is not None:
            if sys.exc_info()[0] is None:
                close_evidence_destination(owned_destination)
            else:
                try:
                    close_evidence_destination(owned_destination)
                except OSError:
                    pass


def canonical_destination(path):
    absolute = os.path.abspath(path)
    return os.path.join(
        os.path.realpath(os.path.dirname(absolute) or "."),
        os.path.basename(absolute))


def canonicalize_output_paths(paths):
    canonical = []
    for path in paths:
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        canonical.append(canonical_destination(path))
    if len(set(canonical)) != len(canonical):
        raise ValueError("evidence destinations must be distinct")
    return canonical


def prepare_output_paths(paths):
    canonical = canonicalize_output_paths(paths)
    for destination in canonical:
        prepared = prepare_evidence_destination(destination)
        close_evidence_destination(prepared)
    return canonical


def validate_path_binding(path, binding, *, destination=None,
                          require_sampler_uid=True):
    owned_destination = None
    if binding["nlink"] != 1:
        raise RuntimeError("bound evidence must be single-link")
    if require_sampler_uid and binding["uid"] != os.geteuid():
        raise RuntimeError(
            "bound evidence must be owned by the sampler UID")
    if destination is None:
        owned_destination = prepare_evidence_destination(
            path, require_absent=False)
        destination = owned_destination
    elif canonical_destination(path) != destination["path"]:
        raise ValueError("evidence path does not match its bound destination")
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | \
        getattr(os, "O_NOFOLLOW", 0)
    try:
        validate_secure_parent(destination)
        fd = os.open(
            destination["basename"], flags, dir_fd=destination["dir_fd"])
        try:
            observed = binding_from_fd(fd)
        finally:
            os.close(fd)
        if observed != binding:
            raise RuntimeError(
                f"evidence path binding or content changed: {path}")
        post = os.stat(
            destination["basename"], dir_fd=destination["dir_fd"],
            follow_symlinks=False)
        if _stat_metadata(post) != {
                key: observed[key]
                for key in (
                    "device", "gid", "inode", "mode", "nlink", "size",
                    "uid")
        }:
            raise RuntimeError(
                f"evidence path changed after validation opened it: {path}")
        validate_secure_parent(destination)
    finally:
        if owned_destination is not None:
            close_evidence_destination(owned_destination)


def unlink_bound_destination(destination, binding):
    """Delete only a sealed name inside its held, private parent directory."""
    validate_path_binding(
        destination["path"], binding, destination=destination)
    os.unlink(destination["basename"], dir_fd=destination["dir_fd"])
    os.fsync(destination["dir_fd"])
    try:
        os.stat(
            destination["basename"], dir_fd=destination["dir_fd"],
            follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        raise RuntimeError("deleted evidence name was unexpectedly recreated")
    validate_secure_parent(destination)


def write_terminal_receipt(path, receipt, *, destination=None, output=None):
    without_self_hash = dict(receipt)
    receipt["self_sha256_excluding_field"] = hashlib.sha256(
        canonical_json_bytes(without_self_hash)).hexdigest()
    owned_destination = None
    owned_output = None
    if destination is None:
        owned_destination = prepare_evidence_destination(path)
        destination = owned_destination
    if output is None:
        owned_output = open_exclusive_evidence(
            path, newline="", destination=destination)
        output = owned_output
    try:
        # A running sampler may have published the one-byte intent latch into
        # this already-bound inode.  Replace it in place: creating or renaming
        # a second pathname would discard the controller's held-inode proof.
        output.flush()
        output.seek(0)
        output.write(canonical_json_bytes(receipt).decode("ascii"))
        output.flush()
        output.truncate(output.tell())
        binding = seal_and_hash_output(output)
        validate_path_binding(path, binding, destination=destination)
        return binding
    finally:
        if owned_output is not None:
            owned_output.close()
        if owned_destination is not None:
            close_evidence_destination(owned_destination)


def write_terminal_intent_marker(output):
    """Publish the fixed one-byte pre-cleanup intent latch exactly once."""
    output.flush()
    if output.tell() != 0 or os.fstat(output.fileno()).st_size != 0:
        raise RuntimeError("terminal intent receipt inode is not empty")
    written = output.write(TERMINAL_INTENT_MARKER.decode("ascii"))
    output.flush()
    if written != 1 or os.fstat(output.fileno()).st_size != 1:
        raise RuntimeError("terminal intent marker write was not exactly one byte")


def validate_terminal_intent_state(
        path, output, reservation_binding, marker_written, destination):
    """Bind the empty-or-marked receipt inode before canonical replacement."""
    output.flush()
    observed = binding_from_fd(output.fileno())
    for key in ("device", "gid", "inode", "mode", "nlink", "uid"):
        if observed[key] != reservation_binding[key]:
            raise RuntimeError("terminal intent receipt inode identity changed")
    expected = TERMINAL_INTENT_MARKER if marker_written else b""
    if observed["size"] != len(expected) or observed["sha256"] != \
            hashlib.sha256(expected).hexdigest():
        raise RuntimeError("terminal intent receipt content differs")
    validate_path_binding(path, observed, destination=destination)
    return observed


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--pid-file", required=True)
    parser.add_argument("--validation-jsonl", required=True)
    parser.add_argument("--receipt", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-output-owner-uid", required=True)
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--dimm-attempts", type=int, default=5)
    parser.add_argument("--dimm-retry-delay", type=float, default=0.01)
    args = parser.parse_args(argv)
    try:
        validate_sampling_arguments(
            args.interval, args.dimm_attempts, args.dimm_retry_delay)
    except ValueError as exc:
        parser.error(str(exc))
    if len(args.expected_source_sha256) != 64 or \
            any(character not in "0123456789abcdef"
                for character in args.expected_source_sha256):
        parser.error("--expected-source-sha256 must be 64 lowercase hex digits")
    if re.fullmatch(
            r"(?:0|[1-9][0-9]*)", args.expected_output_owner_uid) is None:
        parser.error("--expected-output-owner-uid must be a canonical uint")
    args.expected_output_owner_uid = int(args.expected_output_owner_uid)
    if args.expected_output_owner_uid > MAX_UID:
        parser.error("--expected-output-owner-uid is outside the UID domain")
    return args


def decision_exit_code(decision):
    return {
        "continue": EXIT_OK,
        "edac_abort": EXIT_EDAC_ABORT,
        "telemetry_abort": EXIT_TELEMETRY_ABORT,
        "thermal_abort": EXIT_THERMAL_ABORT,
    }[decision]


def run_sampler(args):
    # The receipt file is the only prerequisite which cannot report its own
    # failure.  Reserve it first and hold both its inode and trusted parent for
    # the complete run; every later setup and cleanup operation is receipted.
    receipt_destination = None
    receipt_output = None
    try:
        receipt_path, = canonicalize_output_paths([args.receipt])
        receipt_destination = prepare_evidence_destination(
            receipt_path,
            expected_owner_uid=args.expected_output_owner_uid)
        receipt_output = open_exclusive_evidence(
            receipt_path, newline="", destination=receipt_destination)
        receipt_reservation_binding = binding_from_fd(receipt_output.fileno())
        validate_path_binding(
            receipt_path, receipt_reservation_binding,
            destination=receipt_destination)
    except BaseException as exc:
        if receipt_output is not None:
            try:
                receipt_output.close()
            except BaseException:
                pass
        if receipt_destination is not None:
            try:
                close_evidence_destination(receipt_destination)
            except BaseException:
                pass
        print(
            f"failed to reserve terminal thermal receipt: "
            f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_SAMPLER_ERROR

    csv_path = args.csv
    pid_path = args.pid_file
    validation_path = args.validation_jsonl
    source_path = __file__
    destinations = {
        "raw_csv": None,
        "pid_file": None,
        "validation_jsonl": None,
        "sampler_source": None,
    }
    started_utc = None
    started_monotonic_ns = None
    source_fd = -1
    source_binding = None
    source_binding_finished = None
    stop_signal = None
    bus_fds = {}
    pid_file_created = False
    pid_file_removed = False
    pid_binding = None
    raw_output = None
    raw_binding = None
    validation_output = None
    validation_binding = None
    monitor = None
    tctl_path = None
    edac_ce_paths = ()
    edac_ue_paths = ()
    cpu_tctl_max_c = None
    errors = []
    outcome = "sampler_error"
    exit_code = EXIT_SAMPLER_ERROR
    phase = "initialize"
    terminal_intent_attempted = False
    terminal_intent_written = False

    def stop(signum, _frame):
        nonlocal stop_signal
        if stop_signal is None:
            stop_signal = signum

    def record_error(error_phase, exc):
        nonlocal outcome, exit_code
        outcome = "sampler_error"
        exit_code = EXIT_SAMPLER_ERROR
        append_terminal_error(errors, error_phase, exc)

    def publish_terminal_intent():
        nonlocal terminal_intent_attempted, terminal_intent_written
        if terminal_intent_attempted:
            raise RuntimeError("terminal intent marker was attempted twice")
        terminal_intent_attempted = True
        write_terminal_intent_marker(receipt_output)
        terminal_intent_written = True

    try:
        started_utc = utc_timestamp()
        started_monotonic_ns = time.monotonic_ns()
        source_path = canonical_destination(__file__)

        phase = "install_signal_handlers"
        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)

        phase = "prepare_output_paths"
        csv_path, pid_path, validation_path = canonicalize_output_paths([
            args.csv, args.pid_file, args.validation_jsonl])
        if receipt_path in (csv_path, pid_path, validation_path):
            raise ValueError("evidence destinations must be distinct")
        for label, path in (
                ("raw_csv", csv_path),
                ("pid_file", pid_path),
                ("validation_jsonl", validation_path)):
            phase = f"prepare_{label}_destination"
            destinations[label] = prepare_evidence_destination(
                path, expected_owner_uid=args.expected_output_owner_uid)

        phase = "bind_source"
        destinations["sampler_source"] = prepare_evidence_destination(
            source_path, require_absent=False)
        source_fd = os.open(
            destinations["sampler_source"]["basename"],
            os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK |
            getattr(os, "O_NOFOLLOW", 0),
            dir_fd=destinations["sampler_source"]["dir_fd"])
        source_binding = binding_from_fd(source_fd)
        if source_binding["sha256"] != args.expected_source_sha256:
            raise RuntimeError(
                "sampler source does not match externally sealed expected SHA256")
        if source_binding["mode"] & 0o222:
            raise RuntimeError("sampler source is not sealed read-only")
        validate_path_binding(
            source_path, source_binding,
            destination=destinations["sampler_source"],
            require_sampler_uid=False)

        phase = "open_i2c"
        for bus in sorted({bus for bus, _ in DIMMS}):
            bus_fds[bus] = os.open(
                f"/dev/i2c-{bus}", os.O_RDWR | os.O_CLOEXEC)

        phase = "create_pid_file"
        pid_stream = open_exclusive_evidence(
            pid_path, destination=destinations["pid_file"])
        pid_file_created = True
        with pid_stream:
            pid_stream.write(f"{os.getpid()}\n")
            pid_binding = seal_and_hash_output(pid_stream)
            validate_path_binding(
                pid_path, pid_binding, destination=destinations["pid_file"])

        phase = "discover_sensors"
        tctl_path = find_tctl_path()
        if tctl_path is None:
            raise RuntimeError("CPU Tctl sensor is unavailable")
        edac_ce_paths = discover_edac_paths("ce_count")
        edac_ue_paths = discover_edac_paths("ue_count")
        if tuple(os.path.dirname(path) for path in edac_ce_paths) != tuple(
                os.path.dirname(path) for path in edac_ue_paths):
            raise RuntimeError("EDAC CE/UE controller inventories differ")
        edac_ce_baseline = sum_edac("ce_count", edac_ce_paths)
        edac_ue_baseline = sum_edac("ue_count", edac_ue_paths)
        monitor = DimmPlausibilityMonitor(
            edac_ce_baseline, edac_ue_baseline)
        dimm_columns = [dimm_name(sensor) for sensor in DIMMS]
        columns = [
            "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
            *dimm_columns, "dimm_read_errors", "load1", "load5", "load15",
            "edac_ce", "edac_ue",
        ]

        previous_cpu = read_cpu_stat()
        # Give the first utilization sample a real accounting interval.  On a
        # tickless or lightly loaded machine two immediate /proc/stat reads can
        # be identical, which would leave CPU utilization unreceipted.
        next_sample = time.monotonic() + min(args.interval, 0.1)
        phase = "create_streams"
        raw_output = open_exclusive_evidence(
            csv_path, newline="", destination=destinations["raw_csv"],
            live_mode=0o600)
        validation_output = open_exclusive_evidence(
            validation_path, newline="",
            destination=destinations["validation_jsonl"])
        writer = csv.DictWriter(
            raw_output, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        validation_output.write(canonical_json_bytes({
            "expected_output_owner_uid": args.expected_output_owner_uid,
            "raw_columns": columns,
            "sampler_source_expected_sha256": args.expected_source_sha256,
            "sampling": {
                "dimm_attempts": args.dimm_attempts,
                "dimm_retry_delay_s": args.dimm_retry_delay,
                "interval_s": args.interval,
            },
            "schema": "wirehair.wh2.thermal_validation_stream.v1",
            "thresholds": threshold_metadata(),
        }).decode("ascii"))
        raw_output.flush()
        validation_output.flush()
        outcome = "sampling"
        phase = "sampling"

        while stop_signal is None:
            now = time.monotonic()
            if now < next_sample:
                time.sleep(next_sample - now)
            if stop_signal is not None:
                break
            sample_time = time.monotonic()
            current_cpu = read_cpu_stat()
            tctl_text = read_text(tctl_path)
            if tctl_text is None:
                raise RuntimeError("CPU Tctl sensor became unreadable")
            try:
                cpu_tctl_c = float(tctl_text) / 1000.0
            except ValueError as exc:
                raise RuntimeError("CPU Tctl sensor is malformed") from exc
            if not math.isfinite(cpu_tctl_c):
                raise RuntimeError("CPU Tctl sensor is non-finite")
            if cpu_tctl_max_c is None or cpu_tctl_c > cpu_tctl_max_c:
                cpu_tctl_max_c = cpu_tctl_c
            edac_ce = sum_edac("ce_count", edac_ce_paths)
            edac_ue = sum_edac("ue_count", edac_ue_paths)
            temperatures, failed_dimms, attempt_errors = \
                read_dimm_temperatures(
                    bus_fds, args.dimm_attempts, args.dimm_retry_delay)
            row = {
                "utc": utc_timestamp(),
                "monotonic_s": f"{sample_time:.6f}",
                "cpu_busy_pct": cpu_busy_percent(previous_cpu, current_cpu),
                "cpu_avg_mhz": average_cpu_mhz(),
                "cpu_tctl_c": cpu_tctl_c,
                "dimm_read_errors": len(failed_dimms),
                "edac_ce": edac_ce,
                "edac_ue": edac_ue,
            }
            previous_cpu = current_cpu
            load1, load5, load15 = os.getloadavg()
            row.update({"load1": load1, "load5": load5, "load15": load15})
            for sensor, column in zip(DIMMS, dimm_columns):
                row[column] = temperatures.get(sensor, "")
            # The raw row is evidence and is flushed before validation.  A
            # validator failure can therefore never erase the sampled value.
            writer.writerow(row)
            raw_output.flush()
            validation = monitor.observe(
                sample_time, temperatures, attempt_errors, edac_ce, edac_ue)
            validation_output.write(
                canonical_json_bytes(validation).decode("ascii"))
            validation_output.flush()

            if validation["decision"] != "continue":
                publish_terminal_intent()
                outcome = validation["decision"]
                exit_code = decision_exit_code(outcome)
                break
            next_sample += args.interval
            if next_sample <= sample_time:
                next_sample = sample_time + args.interval

        if outcome == "sampling":
            if stop_signal is not None:
                publish_terminal_intent()
            outcome = "stopped"
            exit_code = EXIT_OK
    except BaseException as exc:
        record_error(phase, exc)
        if not terminal_intent_attempted:
            try:
                publish_terminal_intent()
            except BaseException as intent_exc:
                record_error("publish_terminal_intent", intent_exc)
    finally:
        for label, output, path in (
                ("raw_csv", raw_output, csv_path),
                ("validation_jsonl", validation_output, validation_path)):
            if output is None:
                continue
            try:
                binding = seal_and_hash_output(output)
                validate_path_binding(
                    path, binding, destination=destinations[label])
                if label == "raw_csv":
                    raw_binding = binding
                else:
                    validation_binding = binding
            except BaseException as exc:
                record_error(f"seal_{label}", exc)
            finally:
                try:
                    output.close()
                except BaseException as exc:
                    record_error(f"close_{label}", exc)

        for bus, fd in sorted(bus_fds.items()):
            try:
                os.close(fd)
            except BaseException as exc:
                record_error(f"close_i2c_{bus}", exc)

        if pid_file_created:
            try:
                if pid_binding is None:
                    raise RuntimeError("PID evidence was not sealed")
                unlink_bound_destination(
                    destinations["pid_file"], pid_binding)
                pid_file_removed = True
            except BaseException as exc:
                record_error("unlink_pid_file", exc)

    for label, path, binding in (
            ("raw_csv", csv_path, raw_binding),
            ("validation_jsonl", validation_path, validation_binding)):
        if binding is None:
            continue
        try:
            validate_path_binding(
                path, binding, destination=destinations[label])
        except BaseException as exc:
            record_error(f"terminal_revalidate_{label}", exc)

    if source_fd >= 0:
        try:
            source_binding_finished = binding_from_fd(source_fd)
            if source_binding_finished != source_binding:
                raise RuntimeError(
                    "open sampler source inode changed while running")
            if source_binding_finished["sha256"] != \
                    args.expected_source_sha256:
                raise RuntimeError(
                    "sampler source no longer matches expected SHA256")
            validate_path_binding(
                source_path, source_binding,
                destination=destinations["sampler_source"],
                require_sampler_uid=False)
        except BaseException as exc:
            record_error("terminal_revalidate_source", exc)
        finally:
            try:
                os.close(source_fd)
                source_fd = -1
            except BaseException as exc:
                record_error("close_source", exc)

    # Validate and close every non-receipt parent before serializing the
    # terminal record, so directory-binding failures are themselves recorded.
    for label in (
            "raw_csv", "pid_file", "validation_jsonl", "sampler_source"):
        destination = destinations[label]
        if destination is None:
            continue
        try:
            validate_secure_parent(destination)
        except BaseException as exc:
            record_error(f"terminal_revalidate_{label}_parent", exc)
        try:
            close_evidence_destination(destination)
        except BaseException as exc:
            record_error(f"close_{label}_parent", exc)

    try:
        terminal_summary = monitor.summary() if monitor is not None else None
    except BaseException as exc:
        record_error("build_terminal_summary", exc)
        terminal_summary = None
    try:
        terminal_signal = signal.Signals(stop_signal).name \
            if stop_signal is not None else None
    except BaseException as exc:
        record_error("build_terminal_signal", exc)
        terminal_signal = None
    try:
        finished_monotonic_ns = time.monotonic_ns()
        finished_utc = utc_timestamp()
    except BaseException as exc:
        record_error("build_terminal_timestamps", exc)
        finished_monotonic_ns = None
        finished_utc = None

    receipt = {
        "argv": [
            "--csv", csv_path,
            "--pid-file", pid_path,
            "--validation-jsonl", validation_path,
            "--receipt", receipt_path,
            "--expected-source-sha256", args.expected_source_sha256,
            "--expected-output-owner-uid",
            str(args.expected_output_owner_uid),
            "--interval", repr(args.interval),
            "--dimm-attempts", str(args.dimm_attempts),
            "--dimm-retry-delay", repr(args.dimm_retry_delay),
        ],
        "cpu_tctl_max_c": cpu_tctl_max_c,
        "edac_ce_paths": list(edac_ce_paths),
        "edac_ue_paths": list(edac_ue_paths),
        "errors": errors,
        "expected_output_owner_uid": args.expected_output_owner_uid,
        "exit_code": exit_code,
        "finished_monotonic_ns": finished_monotonic_ns,
        "finished_utc": finished_utc,
        "outcome": outcome,
        "pid": os.getpid(),
        "pid_file": {
            "binding": pid_binding,
            "destination": evidence_destination_metadata(
                destinations["pid_file"]),
            "path": pid_path,
            "removed": pid_file_removed,
        },
        "raw_csv": {
            "binding": raw_binding,
            "destination": evidence_destination_metadata(
                destinations["raw_csv"]),
            "path": csv_path,
        },
        "raw_samples_preserved": raw_binding is not None,
        "receipt_file": {
            "destination": evidence_destination_metadata(receipt_destination),
            "path": receipt_path,
            "reservation_binding": receipt_reservation_binding,
        },
        "sampler_source": {
            "binding": source_binding,
            "binding_finished": source_binding_finished,
            "destination": evidence_destination_metadata(
                destinations["sampler_source"]),
            "expected_sha256": args.expected_source_sha256,
            "path": source_path,
            "sha256": source_binding["sha256"]
                if source_binding is not None else None,
        },
        "sampling": {
            "dimm_attempts": args.dimm_attempts,
            "dimm_retry_delay_s": args.dimm_retry_delay,
            "interval_s": args.interval,
        },
        "schema": SAMPLER_SCHEMA,
        "signal": terminal_signal,
        "started_monotonic_ns": started_monotonic_ns,
        "started_utc": started_utc,
        "summary": terminal_summary,
        "thresholds": threshold_metadata(),
        "validation_jsonl": {
            "binding": validation_binding,
            "destination": evidence_destination_metadata(
                destinations["validation_jsonl"]),
            "path": validation_path,
        },
    }
    published_exit_code = exit_code
    try:
        validate_terminal_intent_state(
            receipt_path, receipt_output, receipt_reservation_binding,
            terminal_intent_written, receipt_destination)
        _terminal_receipt_binding = write_terminal_receipt(
            receipt_path, receipt, destination=receipt_destination,
            output=receipt_output)
    except BaseException as exc:
        print(
            f"failed to write terminal thermal receipt: "
            f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_SAMPLER_ERROR
    finally:
        try:
            receipt_output.close()
        except BaseException as exc:
            print(
                f"failed to close terminal thermal receipt: "
                f"{type(exc).__name__}: {exc}", file=sys.stderr)
        try:
            close_evidence_destination(receipt_destination)
        except BaseException as exc:
            print(
                f"failed to close terminal thermal receipt parent: "
                f"{type(exc).__name__}: {exc}", file=sys.stderr)
    return published_exit_code


def main(argv=None):
    return run_sampler(parse_arguments(argv))


if __name__ == "__main__":
    raise SystemExit(main())
