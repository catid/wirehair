#!/usr/bin/env python3
"""Low-overhead post-EXPO CPU, DDR5, utilization, and EDAC sampler."""

import argparse
import csv
import ctypes
import fcntl
import glob
import os
import signal
import time
from datetime import datetime, timezone


I2C_RDWR = 0x0707
I2C_M_RD = 0x0001
DIMMS = [(1, address) for address in range(0x50, 0x54)] + [
    (2, address) for address in range(0x50, 0x54)
]


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
                return label_path.removesuffix("_label") + "_input"
        fallback = os.path.join(directory, "temp1_input")
        if os.path.exists(fallback):
            return fallback
    return None


def read_cpu_stat():
    raw = read_text("/proc/stat")
    if raw is None or not raw.splitlines():
        raise RuntimeError("/proc/stat is unavailable")
    fields = raw.splitlines()[0].split()[1:]
    ticks = [int(value) for value in fields]
    idle = ticks[3] + (ticks[4] if len(ticks) > 4 else 0)
    return sum(ticks), idle


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


def sum_edac(counter):
    total = 0
    found = False
    for path in glob.glob(f"/sys/devices/system/edac/mc/mc*/{counter}"):
        value = read_text(path)
        if value is not None:
            try:
                total += int(value)
                found = True
            except ValueError:
                pass
    return total if found else ""


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
    """Read every DIMM, retrying only transiently failed I2C transactions."""
    temperatures = {}
    pending = list(DIMMS)
    for attempt in range(attempts):
        failed = []
        for bus, address in pending:
            try:
                temperatures[(bus, address)] = read_spd5118_temperature(
                    bus_fds[bus], address
                )
            except OSError:
                failed.append((bus, address))
        pending = failed
        if not pending:
            break
        if attempt + 1 < attempts:
            time.sleep(retry_delay)
    return temperatures, pending


def utc_timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--pid-file", required=True)
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--dimm-attempts", type=int, default=5)
    parser.add_argument("--dimm-retry-delay", type=float, default=0.01)
    args = parser.parse_args()
    if args.interval <= 0.0:
        parser.error("--interval must be positive")
    if args.dimm_attempts <= 0:
        parser.error("--dimm-attempts must be positive")
    if args.dimm_retry_delay < 0.0:
        parser.error("--dimm-retry-delay must be nonnegative")

    running = True

    def stop(_signum, _frame):
        nonlocal running
        running = False

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    bus_fds = {}
    pid_file_created = False
    try:
        for bus in {bus for bus, _ in DIMMS}:
            bus_fds[bus] = os.open(f"/dev/i2c-{bus}", os.O_RDWR)
        os.makedirs(os.path.dirname(args.csv) or ".", exist_ok=True)
        os.makedirs(os.path.dirname(args.pid_file) or ".", exist_ok=True)
        with open(args.pid_file, "x", encoding="ascii") as pid_stream:
            pid_stream.write(f"{os.getpid()}\n")
        pid_file_created = True

        tctl_path = find_tctl_path()
        if tctl_path is None:
            raise RuntimeError("CPU Tctl sensor is unavailable")
        dimm_columns = [
            f"dimm_i2c{bus}_{address:02x}_c" for bus, address in DIMMS
        ]
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
        with open(args.csv, "x", newline="", encoding="ascii") as output:
            writer = csv.DictWriter(output, fieldnames=columns)
            writer.writeheader()
            output.flush()
            while running:
                now = time.monotonic()
                if now < next_sample:
                    time.sleep(next_sample - now)
                if not running:
                    break
                sample_time = time.monotonic()
                current_cpu = read_cpu_stat()
                tctl_text = read_text(tctl_path) if tctl_path else None
                row = {
                    "utc": utc_timestamp(),
                    "monotonic_s": f"{sample_time:.6f}",
                    "cpu_busy_pct": cpu_busy_percent(previous_cpu, current_cpu),
                    "cpu_avg_mhz": average_cpu_mhz(),
                    "cpu_tctl_c": float(tctl_text) / 1000.0
                        if tctl_text is not None else "",
                    "dimm_read_errors": 0,
                    "edac_ce": sum_edac("ce_count"),
                    "edac_ue": sum_edac("ue_count"),
                }
                previous_cpu = current_cpu
                load1, load5, load15 = os.getloadavg()
                row.update({"load1": load1, "load5": load5, "load15": load15})
                temperatures, failed_dimms = read_dimm_temperatures(
                    bus_fds, args.dimm_attempts, args.dimm_retry_delay
                )
                failed_set = set(failed_dimms)
                for (bus, address), column in zip(DIMMS, dimm_columns):
                    if (bus, address) in failed_set:
                        row[column] = ""
                        row["dimm_read_errors"] += 1
                    else:
                        row[column] = temperatures[(bus, address)]
                writer.writerow(row)
                output.flush()
                next_sample += args.interval
                if next_sample <= sample_time:
                    next_sample = sample_time + args.interval
    finally:
        for fd in bus_fds.values():
            os.close(fd)
        if pid_file_created:
            try:
                os.unlink(args.pid_file)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    main()
