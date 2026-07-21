#!/usr/bin/env python3
"""Audit and enforce a single wirehair I2C thermal reader on this host.

Benchmark campaigns launch the DIMM/CPU thermal sampler
(``wirehair_expo_thermal_sampler.py``) and CPU load fillers
(``wirehair_load_fillers.sh``) through ``sudo -n``, so those processes are
root-owned.  A controller that dies without reaching its cleanup path
strands them: a plain ``kill(pid, 0)`` probe from the non-root
orchestrator raises EPERM, and any code that conflates EPERM with ESRCH
concludes the process is dead and leaks it.  On 2026-07-21 exactly that
had happened: one sampler had been orphaned for two days and 128 root
busy-loop fillers had saturated the host for 22 hours.

This standalone tool gives every launcher one audited preflight:

* liveness that never conflates EPERM with ESRCH, using ``sudo -n
  kill -0`` plus ``/proc`` start-tick identity, failing closed when a
  protected PID cannot be classified;
* an exact-cmdline ``/proc`` scan proving how many I2C readers exist;
* ``--mode audit`` reports without side effects; ``--mode enforce``
  fails unless the live reader set is exactly the one identified by
  ``--expect-pid-file`` (or empty when no expectation is given);
* ``--mode stop-orphans`` terminates stranded readers and fillers
  through ``sudo -n``, re-proving each PID's start-tick identity before
  every signal so a reused PID is never killed;
* a JSON receipt recording PID, owner, CPU affinity, start tick and the
  full command line for every match and every action taken.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import stat as stat_module
import subprocess
import sys
import time
from typing import Any, Callable, Mapping, Sequence

SAMPLER_BASENAME = "wirehair_expo_thermal_sampler.py"
FILLER_SIGNATURES = ("wirehair_load_fillers.sh", "while :; do :; done")
TERM_WAIT_SECONDS = 10.0
KILL_WAIT_SECONDS = 5.0
POLL_INTERVAL_SECONDS = 0.25
OWNER_MARKER_PATH = Path("/tmp/wirehair-environment-owner.json")
OWNER_MARKER_SCHEMA = "wirehair.environment_owner.v1"
BOOT_ID_PATH = Path("/proc/sys/kernel/random/boot_id")


class JanitorError(Exception):
    """Raised for environment or verification failures."""


def die(message: str) -> "NoReturn":  # noqa: F821 - documentation type
    raise JanitorError(message)


def _read_cmdline(proc_root: Path, pid: int) -> tuple[str, ...] | None:
    try:
        raw = (proc_root / str(pid) / "cmdline").read_bytes()
    except OSError:
        return None
    if not raw:
        return None
    return tuple(
        token.decode("utf-8", "surrogateescape")
        for token in raw.rstrip(b"\0").split(b"\0"))


def _read_start_ticks(proc_root: Path, pid: int) -> int | None:
    """Return the process start time in clock ticks, or None if gone.

    Field 22 of ``/proc/<pid>/stat`` follows the parenthesised comm,
    which may itself contain spaces and parentheses, so parse from the
    last closing parenthesis.
    """
    try:
        text = (proc_root / str(pid) / "stat").read_text(
            encoding="ascii", errors="replace")
    except OSError:
        return None
    tail = text.rpartition(")")[2].split()
    if len(tail) < 20:
        return None
    try:
        return int(tail[19])
    except ValueError:
        return None


def _read_uid(proc_root: Path, pid: int) -> int | None:
    try:
        text = (proc_root / str(pid) / "status").read_text(
            encoding="ascii", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if line.startswith("Uid:"):
            fields = line.split()
            if len(fields) >= 2 and fields[1].isdigit():
                return int(fields[1])
    return None


def _read_affinity(
    pid: int,
    affinity_fn: Callable[[int], set] = os.sched_getaffinity,
) -> tuple[int, ...] | None:
    try:
        return tuple(sorted(affinity_fn(pid)))
    except OSError:
        return None


WRAPPER_BASENAMES = frozenset((
    "sudo", "env", "taskset", "numactl", "nice", "setpriv", "setsid",
    "sh", "bash", "timeout",
))


def classify_command(cmdline: Sequence[str]) -> str | None:
    """Classify a command line as ``reader``, ``filler`` or None.

    A launcher chain like ``sudo -b env taskset -c 127 python3 sampler``
    leaves the wrapper alive as the sampler's parent with the sampler
    path in its own cmdline; only the process whose executable is the
    interpreter actually holds the I2C bus, so wrapper processes are
    not counted as readers.
    """
    joined = " ".join(cmdline)
    if any(Path(token).name == SAMPLER_BASENAME for token in cmdline):
        if cmdline and Path(cmdline[0]).name in WRAPPER_BASENAMES:
            return None
        return "reader"
    if any(signature in joined for signature in FILLER_SIGNATURES):
        return "filler"
    return None


def pid_liveness(
    pid: int,
    *,
    kill_fn: Callable[[int, int], None] = os.kill,
    run_fn: Callable[..., subprocess.CompletedProcess] = subprocess.run,
    proc_root: Path = Path("/proc"),
) -> str:
    """Return ``alive``, ``dead`` or ``unknown`` without conflating errors.

    EPERM proves only that the caller may not signal the PID; a
    root-owned process is very much alive when it happens.  Escalate
    through ``sudo -n kill -0`` and, if that also cannot answer, report
    ``unknown`` so callers fail closed instead of leaking the process.
    """
    try:
        kill_fn(pid, 0)
    except ProcessLookupError:
        return "dead"
    except PermissionError:
        result = run_fn(
            ("sudo", "-n", "kill", "-0", str(pid)),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            check=False)
        if result.returncode == 0:
            return "alive"
        if (proc_root / str(pid)).exists():
            return "unknown"
        return "dead"
    return "alive"


def scan_processes(
    *,
    proc_root: Path = Path("/proc"),
    affinity_fn: Callable[[int], set] = os.sched_getaffinity,
) -> list[dict[str, Any]]:
    """Inventory every thermal reader and load filler visible in /proc."""
    records: list[dict[str, Any]] = []
    for entry in sorted(proc_root.glob("[0-9]*")):
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        cmdline = _read_cmdline(proc_root, pid)
        if cmdline is None:
            continue
        role = classify_command(cmdline)
        if role is None:
            continue
        start_ticks = _read_start_ticks(proc_root, pid)
        if start_ticks is None:
            continue
        records.append({
            "pid": pid,
            "role": role,
            "uid": _read_uid(proc_root, pid),
            "start_ticks": start_ticks,
            "affinity": _read_affinity(pid, affinity_fn),
            "cmdline": list(cmdline),
        })
    return records


def read_expected_pid(pid_file: Path) -> int:
    """Parse a canonical sampler PID file, refusing surprising content."""
    try:
        metadata = pid_file.lstat()
    except OSError as exc:
        die(f"cannot read expected PID file: {exc}")
    if not stat_module.S_ISREG(metadata.st_mode):
        die(f"expected PID file is not a regular file: {pid_file}")
    text = pid_file.read_text(encoding="ascii")
    if not re.fullmatch(r"[1-9][0-9]*\n?", text):
        die(f"expected PID file is not canonical: {pid_file}")
    return int(text)


def _identity_state(
    record: Mapping[str, Any], proc_root: Path,
) -> str:
    """Return ``present``, ``gone`` or ``reused`` for a scanned record."""
    current = _read_start_ticks(proc_root, int(record["pid"]))
    if current is None:
        return "gone"
    if current != int(record["start_ticks"]):
        return "reused"
    return "present"


def stop_process(
    record: Mapping[str, Any],
    *,
    proc_root: Path = Path("/proc"),
    run_fn: Callable[..., subprocess.CompletedProcess] = subprocess.run,
    sleep_fn: Callable[[float], None] = time.sleep,
    monotonic_fn: Callable[[], float] = time.monotonic,
) -> dict[str, Any]:
    """TERM then KILL one scanned process, re-proving identity each step."""
    pid = int(record["pid"])
    action: dict[str, Any] = {
        "pid": pid,
        "role": record["role"],
        "start_ticks": record["start_ticks"],
        "signals": [],
        "outcome": None,
    }
    for signal_name, wait in (
            ("-TERM", TERM_WAIT_SECONDS), ("-KILL", KILL_WAIT_SECONDS)):
        state = _identity_state(record, proc_root)
        if state != "present":
            action["outcome"] = ("already-exited" if not action["signals"]
                                 else "stopped")
            if state == "reused":
                action["outcome"] = "pid-reused"
            return action
        result = run_fn(
            ("sudo", "-n", "kill", signal_name, str(pid)),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            check=False)
        action["signals"].append(
            {"signal": signal_name, "returncode": result.returncode})
        deadline = monotonic_fn() + wait
        while monotonic_fn() < deadline:
            if _identity_state(record, proc_root) != "present":
                action["outcome"] = "stopped"
                return action
            sleep_fn(POLL_INTERVAL_SECONDS)
    if _identity_state(record, proc_root) == "present":
        action["outcome"] = "unkillable"
    else:
        action["outcome"] = "stopped"
    return action


def load_owner_marker(
    *,
    marker_path: Path = OWNER_MARKER_PATH,
    boot_id_path: Path = BOOT_ID_PATH,
    now_utc: str | None = None,
) -> dict[str, Any] | None:
    """Load the environment-ownership marker, failing closed.

    A live campaign writes this marker so cleanup tooling never kills its
    sealed environment (stopping the sampler or fillers formally
    invalidates the campaign).  A marker that is present but unreadable
    or malformed is treated as ACTIVE with an ``ambiguous`` flag — the
    one direction that can never cost a multi-day run.  Markers from a
    previous boot or past their expiry are inert.
    """
    try:
        raw = marker_path.read_text(encoding="ascii")
    except FileNotFoundError:
        return None
    except OSError:
        return {"ambiguous": True}
    try:
        marker = json.loads(raw)
    except ValueError:
        return {"ambiguous": True}
    if not isinstance(marker, dict) or \
            marker.get("schema") != OWNER_MARKER_SCHEMA:
        return {"ambiguous": True}
    if marker.get("phase") == "complete":
        return None
    try:
        boot_id = boot_id_path.read_text(encoding="ascii").strip()
    except OSError:
        boot_id = None
    if boot_id is not None and marker.get("boot_id") not in (None, boot_id):
        return None
    expires = marker.get("expires_utc")
    current = now_utc
    if current is None:
        result = subprocess.run(
            ("date", "-u", "+%Y-%m-%dT%H:%M:%S+00:00"),
            stdout=subprocess.PIPE, check=False)
        current = result.stdout.decode("ascii", "replace").strip()
    if isinstance(expires, str) and current and expires < current:
        return None
    return marker


def _protected_identities(marker: Mapping[str, Any]) -> set[tuple[int, int]]:
    identities: set[tuple[int, int]] = set()
    for entry in marker.get("protected", ()):
        if isinstance(entry, dict) and \
                isinstance(entry.get("pid"), int) and \
                isinstance(entry.get("start_ticks"), int):
            identities.add((entry["pid"], entry["start_ticks"]))
    return identities


def sudo_available(
    run_fn: Callable[..., subprocess.CompletedProcess] = subprocess.run,
) -> bool:
    result = run_fn(
        ("sudo", "-n", "true"),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    return result.returncode == 0


def run_janitor(
    mode: str,
    *,
    expect_pid_file: Path | None = None,
    override_owner: str | None = None,
    proc_root: Path = Path("/proc"),
    kill_fn: Callable[[int, int], None] = os.kill,
    run_fn: Callable[..., subprocess.CompletedProcess] = subprocess.run,
    affinity_fn: Callable[[int], set] = os.sched_getaffinity,
    sleep_fn: Callable[[float], None] = time.sleep,
    monotonic_fn: Callable[[], float] = time.monotonic,
    marker_path: Path = OWNER_MARKER_PATH,
    boot_id_path: Path = BOOT_ID_PATH,
    now_utc: str | None = None,
) -> dict[str, Any]:
    """Execute one janitor pass and return the receipt."""
    marker = load_owner_marker(
        marker_path=marker_path, boot_id_path=boot_id_path, now_utc=now_utc)
    marker_ambiguous = bool(marker and marker.get("ambiguous"))
    protected = (_protected_identities(marker)
                 if marker and not marker_ambiguous else set())

    filler_leader_protected = bool(marker and not marker_ambiguous and any(
        isinstance(entry, dict) and entry.get("role") == "filler_leader"
        for entry in marker.get("protected", ())))

    def is_protected(record: Mapping[str, Any]) -> bool:
        if marker_ambiguous:
            return True
        if (int(record["pid"]), int(record["start_ticks"])) in protected:
            return True
        # The marker lists only the filler supervisor; its worker children
        # (one busy loop per CPU) belong to the same sealed environment.
        return record["role"] == "filler" and filler_leader_protected

    records = scan_processes(proc_root=proc_root, affinity_fn=affinity_fn)
    for record in records:
        record["liveness"] = pid_liveness(
            int(record["pid"]), kill_fn=kill_fn, run_fn=run_fn,
            proc_root=proc_root)
    readers = [r for r in records if r["role"] == "reader"]
    fillers = [r for r in records if r["role"] == "filler"]
    receipt: dict[str, Any] = {
        "schema": "wh2-thermal-reader-janitor-v1",
        "mode": mode,
        "sampler_basename": SAMPLER_BASENAME,
        "filler_signatures": list(FILLER_SIGNATURES),
        "readers": readers,
        "fillers": fillers,
        "actions": [],
        "violations": [],
        "owner_marker": (None if marker is None else
                         {"ambiguous": True} if marker_ambiguous else
                         {"campaign_root": marker.get("campaign_root"),
                          "phase": marker.get("phase"),
                          "expires_utc": marker.get("expires_utc")}),
        "protected_skipped": [],
    }
    if marker_ambiguous:
        receipt["violations"].append(
            "environment-owner marker is unreadable or malformed; "
            "treating every matched process as protected")
    unknown = [r for r in records if r["liveness"] == "unknown"]
    for record in unknown:
        receipt["violations"].append(
            "cannot determine liveness of protected PID %d" % record["pid"])
    live_readers = [r for r in readers if r["liveness"] == "alive"]
    if mode == "audit":
        return receipt
    if mode == "enforce":
        if expect_pid_file is not None:
            expected_pid = read_expected_pid(expect_pid_file)
            receipt["expected_pid"] = expected_pid
            matching = [
                r for r in live_readers if int(r["pid"]) == expected_pid]
            extra = [
                r for r in live_readers
                if int(r["pid"]) != expected_pid and not is_protected(r)]
            if not matching:
                receipt["violations"].append(
                    "expected sole I2C reader PID %d is not a live sampler"
                    % expected_pid)
            for record in extra:
                receipt["violations"].append(
                    "unexpected extra I2C reader PID %d" % record["pid"])
        else:
            for record in live_readers:
                if is_protected(record):
                    continue
                receipt["violations"].append(
                    "no I2C reader expected but PID %d is live"
                    % record["pid"])
        for record in fillers:
            if record["liveness"] != "dead" and not is_protected(record):
                receipt["violations"].append(
                    "load filler PID %d is still present" % record["pid"])
        return receipt
    if mode == "stop-orphans":
        targets = [r for r in records if r["liveness"] != "dead"]
        if marker is not None:
            owner_root = None if marker_ambiguous else marker.get(
                "campaign_root")
            if override_owner is None or override_owner != owner_root:
                receipt["protected_skipped"] = [
                    {"pid": r["pid"], "role": r["role"]} for r in targets]
                receipt["violations"].append(
                    "environment-owner marker is active; refusing every "
                    "kill (pass --override-owner with the exact campaign "
                    "root to tear down that campaign's own environment)")
                return receipt
        if targets and not sudo_available(run_fn):
            receipt["violations"].append(
                "sudo -n is unavailable; cannot stop root-owned orphans")
            return receipt
        for record in targets:
            action = stop_process(
                record, proc_root=proc_root, run_fn=run_fn,
                sleep_fn=sleep_fn, monotonic_fn=monotonic_fn)
            receipt["actions"].append(action)
            if action["outcome"] not in ("stopped", "already-exited"):
                receipt["violations"].append(
                    "PID %d was not stopped: %s"
                    % (record["pid"], action["outcome"]))
        return receipt
    die(f"unknown mode: {mode}")


def _write_receipt(receipt: Mapping[str, Any], path: Path) -> None:
    serialized = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    scratch = path.with_name(path.name + ".partial")
    scratch.write_text(serialized, encoding="ascii")
    os.replace(scratch, path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--mode", choices=("audit", "enforce", "stop-orphans"),
        default="audit")
    parser.add_argument(
        "--expect-pid-file", type=Path, default=None,
        help="enforce mode: PID file naming the sole allowed I2C reader")
    parser.add_argument(
        "--override-owner", default=None,
        help="stop-orphans: campaign root whose owner marker authorizes "
             "tearing down its own protected environment")
    parser.add_argument("--json-out", type=Path, default=None)
    arguments = parser.parse_args(argv)
    if arguments.expect_pid_file is not None and arguments.mode != "enforce":
        parser.error("--expect-pid-file requires --mode enforce")
    if arguments.override_owner is not None and \
            arguments.mode != "stop-orphans":
        parser.error("--override-owner requires --mode stop-orphans")
    try:
        receipt = run_janitor(
            arguments.mode, expect_pid_file=arguments.expect_pid_file,
            override_owner=arguments.override_owner)
    except JanitorError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if arguments.json_out is not None:
        _write_receipt(receipt, arguments.json_out)
    json.dump(receipt, sys.stdout, indent=2, sort_keys=True)
    print()
    return 1 if receipt["violations"] else 0


if __name__ == "__main__":
    sys.exit(main())
