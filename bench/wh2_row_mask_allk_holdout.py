#!/usr/bin/env python3
"""Prepare, execute, and strictly reduce the WH2 survivor-mask all-K holdout.

The development row-mask screen (wh2_grouped_row_mask_campaign.py) found the
canonical suffix masks pathological and proposed the structural survivors:
prefix mask 0x007 at P48/r3 and prefix mask 0x07f at P32/r7.  This controller
seals the acceptance holdout: each survivor runs the complete all-K domain
(K 2..64000 x 3 seeds x 3 schedules x 1 trial, loss 0.50, bb 64) side by side
with its canonical-suffix comparator (0x380 at P48/r3, 0x3f8 at P32/r7),
cell-for-cell identical to the sealed finalists campaign so its sealed arm
results (p48_r3 685, p32_r7 687, prod244 691, p48_r0 735 failures) are
directly comparable.

Holdout discipline: this is an exact-comparison holdout on the sealed all-K
seeds and schedules of the finalists campaign; NO new independent seeds are
used.  The 16 development hard cells (the finalist field-shortfall union that
selected the survivor masks) were selection material: the reducer reports
them separately as fix-confirmation and EXCLUDES them from the headline
holdout comparison.  Acceptance per survivor requires a lower raw
field-shortfall count and a lower raw total failure count; on headline cells
it also requires a lower total
failure count, at least one actual repair, and zero introductions.  Paired
XOR/GF-muladd cost ratios on common headline successes must remain within the
sealed tolerance, at least one panel source field-shortfall must be genuinely
fixed rather than relabeled, and zero codec errors are required everywhere.

The preparation path is intentionally separate from execution.  It freezes a
freshly built clean-commit test-hook benchmark, the sole CPU/DIMM thermal
sampler, the validated all-K source receipt, and a complete task ledger.
``launch --preflight-only`` performs no benchmark work.  ``launch --execute``
is the only path that starts jobs; it refuses a second SPD/I2C sampler and a
foreign environment-ownership marker, and declares its own ownership marker
at /tmp/wirehair-environment-owner.json before launching jobs.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import hashlib
import io
import json
import math
import os
import re
import select
import shlex
import shutil
import signal
import stat
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


sys.dont_write_bytecode = True
getcontext().prec = 50

SCHEMA = "wirehair.wh2.row_mask_allk_holdout.v1"
CODEC_IMPLEMENTATION_HEAD = "42fb578b3df159970708c2db51684a9ad1abf93c"
SOURCE_HEAD = "0978602bb535712f6136b94c298ef428e4883fbc"
SOURCE_SUMMARY_SHA256 = (
    "c09d1e9b74a0e6cc38fea67d7038a3a77281736d7fd1be4d2867c48d6e3a2e72"
)
SOURCE_SUMMARY_NAME = "validated_summary.json"
SAMPLER_NAME = "wirehair_expo_thermal_sampler.py"
SAMPLER_SHA256 = "2b84efa91375a96a4a64e09ce5bfd7cba0b85b75028f5a93470cd4ae58aadb01"
BINARY_NAME = "wirehair_v2_bench"
CONTROLLER_NAME = "wh2_row_mask_allk_holdout.py"
DIMM_I2C_DEVICES = (Path("/dev/i2c-1"), Path("/dev/i2c-2"))
SUDO_PATH = Path("/usr/bin/sudo")
FUSER_PATH = Path("/usr/bin/fuser")
KILL_PATH = Path("/usr/bin/kill")
CAT_PATH = Path("/usr/bin/cat")
PYTHON_PATH = Path("/usr/bin/python3.12")
I2C_INSPECTION_TIMEOUT_SECONDS = 5.0
I2C_RECHECK_INTERVAL_SECONDS = 1.0
JOB_TIMEOUT_SECONDS = 900.0
PRIVILEGED_PIDFD_SIGNAL_CODE = r"""
import os
import signal
import sys

pid = int(sys.argv[1])
expected_start = int(sys.argv[2])
expected_cmdline = bytes.fromhex(sys.argv[3])
signum = {"TERM": signal.SIGTERM, "KILL": signal.SIGKILL}[sys.argv[4]]
try:
    pidfd = os.pidfd_open(pid, 0)
except ProcessLookupError:
    raise SystemExit(4)
try:
    try:
        with open("/proc/{}/stat".format(pid), "rb") as stream:
            stat_raw = stream.read()
        with open("/proc/{}/cmdline".format(pid), "rb") as stream:
            cmdline_raw = stream.read()
        observed_start = int(stat_raw.rsplit(b") ", 1)[1].split()[19])
    except (FileNotFoundError, ProcessLookupError):
        raise SystemExit(4)
    except (IndexError, OSError, ValueError):
        raise SystemExit(5)
    if observed_start != expected_start or cmdline_raw != expected_cmdline:
        raise SystemExit(3)
    try:
        signal.pidfd_send_signal(pidfd, signum)
    except ProcessLookupError:
        raise SystemExit(4)
finally:
    os.close(pidfd)
"""
PDEATH_EXEC_CODE = r"""
import ctypes
import os
import signal
import sys

expected_parent = int(sys.argv[1])
expected_parent_start = int(sys.argv[2])
command = sys.argv[3:]
if not command:
    raise SystemExit(126)
libc = ctypes.CDLL(None, use_errno=True)
if libc.prctl(1, signal.SIGKILL, 0, 0, 0) != 0:
    raise SystemExit(126)
# PR_SET_PDEATHSIG alone has a fork-to-prctl race.  Comparing the parent after
# arming it closes that race: if the controller already died, never exec the
# benchmark.
try:
    with open("/proc/{}/stat".format(expected_parent), "rb") as stream:
        parent_stat = stream.read()
    parent_start = int(parent_stat.rsplit(b") ", 1)[1].split()[19])
except (FileNotFoundError, IndexError, OSError, ValueError):
    raise SystemExit(125)
if os.getppid() != expected_parent or parent_start != expected_parent_start:
    raise SystemExit(125)
os.execv(command[0], command)
"""
PRIVILEGED_SAMPLER_SUPERVISOR_CODE = r"""
import os
import select
import signal
import subprocess
import sys
import time

controller_pid = int(sys.argv[1])
controller_start = int(sys.argv[2])
pdeath_code = bytes.fromhex(sys.argv[3]).decode("utf-8")
command = sys.argv[4:]
if not command:
    raise SystemExit(126)

try:
    controller_pidfd = os.pidfd_open(controller_pid, 0)
    with open("/proc/{}/stat".format(controller_pid), "rb") as stream:
        raw = stream.read()
    observed_start = int(raw.rsplit(b") ", 1)[1].split()[19])
except (FileNotFoundError, IndexError, OSError, ProcessLookupError, ValueError):
    raise SystemExit(125)
if observed_start != controller_start:
    os.close(controller_pidfd)
    raise SystemExit(125)

poller = select.poll()
poller.register(
    controller_pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
stopping = False

def request_stop(_signum, _frame):
    global stopping
    stopping = True

signal.signal(signal.SIGTERM, request_stop)
signal.signal(signal.SIGINT, request_stop)

def controller_alive():
    return not poller.poll(0)

def group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True

def observe_group_gone(pgid, timeout):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not group_alive(pgid):
            return True
        time.sleep(0.05)
    return not group_alive(pgid)

def leader_ready(poller, timeout):
    return bool(poller.poll(max(0, int(timeout * 1000))))

def stop_child(child, child_poller):
    # start_new_session=True makes child.pid the owned process-group ID.
    # Keep the leader unreaped until the final group signal: its process ID
    # pins the PGID and prevents any unrelated group from reusing that number.
    # The pinned sampler is single-process, while the final pre-reap SIGKILL
    # also cleans up any unexpected same-group descendant.
    if child.returncode is not None:
        if not observe_group_gone(child.pid, 5.0):
            raise RuntimeError("reaped sampler left a live process group")
        return child.returncode
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    leader_ready(child_poller, 15.0)
    # Always issue the last group signal before wait()/poll() can reap the
    # leader.  SIGKILL is harmless to an already-dead unreaped leader.
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    if not leader_ready(child_poller, 5.0):
        raise RuntimeError("sampler leader did not stop")
    returncode = child.wait(timeout=5.0)
    if not observe_group_gone(child.pid, 5.0):
        raise RuntimeError("sampler process group did not stop")
    return returncode

def force_child_without_pidfd(child):
    # The direct child is still unreaped here, so child.pid cannot be reused.
    if child.returncode is not None:
        if not observe_group_gone(child.pid, 5.0):
            raise RuntimeError(
                "reaped unbound sampler left a live process group")
        return child.returncode
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    returncode = child.wait(timeout=5.0)
    if not observe_group_gone(child.pid, 5.0):
        raise RuntimeError("unbound sampler process group did not stop")
    return returncode

child = None
child_pidfd = None
child_poller = None
child_cleaned = False
try:
    if stopping or not controller_alive():
        raise SystemExit(125)
    with open("/proc/self/stat", "rb") as stream:
        own_stat = stream.read()
    own_start = int(own_stat.rsplit(b") ", 1)[1].split()[19])
    if stopping or not controller_alive():
        raise SystemExit(125)
    child = subprocess.Popen([
        sys.executable, "-I", "-c", pdeath_code,
        str(os.getpid()), str(own_start), *command,
    ], start_new_session=True)
    try:
        child_pidfd = os.pidfd_open(child.pid, 0)
    except (OSError, ProcessLookupError):
        force_child_without_pidfd(child)
        child_cleaned = True
        raise SystemExit(126)
    child_poller = select.poll()
    child_poller.register(
        child_pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
    # Close controller-death races on both sides of child creation.
    if stopping or not controller_alive():
        stop_child(child, child_poller)
        child_cleaned = True
        raise SystemExit(125)
    while not leader_ready(child_poller, 0) and not stopping and \
            controller_alive():
        child_poller.poll(200)
    returncode = stop_child(child, child_poller)
    child_cleaned = True
    raise SystemExit(returncode)
finally:
    if child is not None and not child_cleaned:
        if child_poller is None:
            force_child_without_pidfd(child)
        else:
            stop_child(child, child_poller)
    if child_pidfd is not None:
        os.close(child_pidfd)
    os.close(controller_pidfd)
"""
BENCHMARK_EXEC_READY_TIMEOUT_SECONDS = 5.0

SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
SOURCE_ARMS = ("p48_r3", "p32_r7", "p48_r0", "prod244")
FINALIST_ARMS = ("p48_r3", "p32_r7")
# Sealed totals of the source campaign, verified against its
# validated_summary.json before this controller was sealed.
SOURCE_ARM_FAILURES = {
    "p48_r3": 685, "p32_r7": 687, "p48_r0": 735, "prod244": 691,
}
SOURCE_FIELD_SHORTFALLS = {"p48_r3": 7, "p32_r7": 9}
HARD_CELL_COUNT = 16
CANONICAL_SUFFIX_MASKS = {48: 0x380, 32: 0x3F8}
ARMS = (
    {"name": "p48_r3_pfx007", "panel": "p48_r3", "period": 48, "rows": 3,
     "mask": 0x007, "role": "candidate", "comparator": "p48_r3_sfx380"},
    {"name": "p48_r3_sfx380", "panel": "p48_r3", "period": 48, "rows": 3,
     "mask": 0x380, "role": "comparator", "comparator": None},
    {"name": "p32_r7_pfx07f", "panel": "p32_r7", "period": 32, "rows": 7,
     "mask": 0x07F, "role": "candidate", "comparator": "p32_r7_sfx3f8"},
    {"name": "p32_r7_sfx3f8", "panel": "p32_r7", "period": 32, "rows": 7,
     "mask": 0x3F8, "role": "comparator", "comparator": None},
)
SMOKE_CASES = ((48, 3, 0x007), (48, 3, 0x380), (32, 7, 0x07F), (32, 7, 0x3F8))
K_LO = 2
K_HI = 64_000
CHUNK_MAX = 250
DEFAULT_WORKERS = 126
MIN_CAMPAIGN_WORKERS = 64
MAX_CAMPAIGN_WORKERS = 256
# The paired cost gate: the development screen measured XOR ratios of
# exactly 1.0 and muladd ratios within +/-3e-5 for every alternative
# mask, so a sealed 5e-4 ceiling is ~16x the observed noise while still
# rejecting any real cost regression.
COST_RATIO_MAX = Decimal("1.0005")
# The saturation floor is a coarse tripwire for a collapsed or undersized
# worker pool, not a precision gate: at 126 workers on a 128-CPU host the
# sustained ceiling is 98.4% and the measured whole-segment mean lands
# near 95% once the ramp and drain tails are included, so a 95 floor is
# flaky by construction while 90 still trips when saturation is lost.
CPU_BUSY_FLOOR = Decimal("90")
# The floor is a control-stage guarantee: hard segments are sub-minute
# bursts of the 16 development fix-confirmation cells whose interval busy
# mean is dominated by worker-pool spawn latency.  Every segment still
# passes the coverage, cadence, DIMM, EDAC and sampler-exit gates.
CPU_BUSY_FLOOR_STAGES = ("control",)
THERMAL_INTERVAL_SECONDS = Decimal("1")
THERMAL_MIN_GAP_SECONDS = Decimal("0.5")
THERMAL_MAX_GAP_SECONDS = Decimal("2.5")
THERMAL_READY_SAMPLES = 2
THERMAL_READY_TIMEOUT_SECONDS = 12.0
THERMAL_END_TIMEOUT_SECONDS = 5.0
# At one sample per second, 64 MiB accommodates well over 40 hours even at
# 400 bytes per row.  The sealed controller requires at least 64 workers, so
# this also exceeds the one-segment all-jobs-at-timeout bound; at the intended
# 126 workers it has more than 2x margin.
# The cap keeps a corrupt live artifact from turning a freshness check into an
# unbounded allocation; final reduction applies the same limit.
MAX_THERMAL_STREAM_BYTES = 64 * 1024 * 1024
MAX_THERMAL_STDERR_BYTES = 1024 * 1024
MIN_PLAUSIBLE_CPU_C = Decimal("0")
MAX_PLAUSIBLE_CPU_C = Decimal("130")
MIN_PLAUSIBLE_DIMM_C = Decimal("0")
MAX_PLAUSIBLE_DIMM_C = Decimal("130")
OWNER_MARKER_PATH = Path("/tmp/wirehair-environment-owner.json")
OWNER_LOCK_PATH = Path("/tmp/wirehair-environment-owner.lock")
OWNER_MARKER_SCHEMA = "wirehair.environment_owner.v1"
BOOT_ID_PATH = Path("/proc/sys/kernel/random/boot_id")
DEFAULT_OWNER_TTL_HOURS = 168
BUILD_TARGET = "wirehair_v2_bench"
BUILD_CONFIGURE_ARGS = (
    "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DBUILD_TESTS=ON",
    "-DBUILD_CODEC_V2=ON",
    "-DMARCH_NATIVE=ON",
    "-DWIREHAIR_BUILD_BENCHMARKS=ON",
    "-DWIREHAIR_STRICT_WARNINGS=ON",
    "-DWH_LTO=OFF",
    "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
)

CSV_FIELDS = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor",
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
CSV_TIMING_MILLI_FIELDS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
)
CSV_COUNTER_MU_FIELDS = (
    "block_xors_mu", "block_muladds_mu", "mixed_joint_source_xors_mu",
    "mixed_joint_marginal_xors_mu", "mixed_joint_marginal_copies_mu",
    "mixed_joint_active_deltas_mu", "mixed_joint_scratch_bytes_mu",
    "mixed_dual_source_columns_mu",
)
CSV_MILLI_FIELDS = CSV_TIMING_MILLI_FIELDS + CSV_COUNTER_MU_FIELDS
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
RUNTIME_DIRECTORY_NAMES = (
    "raw", "stderr", "exit", "thermal", "segments", "attempts",
    "job_receipts",
)
STOPPED_PROCESS_ACTIONS = {"already_exited", "terminated", "killed"}
# A start-time mismatch proves that the originally captured process object is
# gone even if its numeric PID has already been reused.  Command-identity loss
# with the same start time is not proof of cleanup and is deliberately absent.
DISCOVERY_PROVED_GONE_ACTIONS = STOPPED_PROCESS_ACTIONS | {"identity_reused"}
# Linux exposes process IDs through a signed pid_t to kill(2)/pidfd_open(2).
# The configured pid_max is normally far lower, but this architectural bound
# prevents untrusted receipts from reaching those APIs as unrepresentable
# Python integers.
MAX_PROCESS_PID = (1 << 31) - 1
MAX_SAMPLER_PID_RECEIPT_BYTES = len(str(MAX_PROCESS_PID)) + 1
SamplerPidPublication = Tuple[int, int, int, bytes]


class CampaignError(RuntimeError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_safe_directory(path: Path, label: str) -> None:
    """Reject indirect, foreign-owned, or writable campaign directories."""
    try:
        identity = path.lstat()
    except OSError as exc:
        die("{} directory is unavailable: {}".format(label, exc))
    if not stat.S_ISDIR(identity.st_mode) or identity.st_uid != os.getuid() or \
            identity.st_mode & 0o022:
        die("{} directory is indirect or unsafe: {}".format(label, path))


def verify_runtime_directories(root: Path) -> None:
    verify_safe_directory(root, "campaign root")
    for name in RUNTIME_DIRECTORY_NAMES:
        verify_safe_directory(root / name, name)


def make_private_directory(
    path: Path,
    *,
    parents: bool = False,
    exist_ok: bool = False,
) -> None:
    path.mkdir(mode=0o700, parents=parents, exist_ok=exist_ok)
    # Be explicit even if an unusual platform or pre-existing test fixture
    # ignored the requested creation mode.
    os.chmod(path, 0o700)


def write_once(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        die("immutable artifact is a symlink {}".format(path))
    if path.exists():
        if path.read_bytes() != data:
            die("refusing to replace immutable artifact {}".format(path))
        return
    temporary = Path(
        str(path) + ".part.{}.{}".format(os.getpid(), threading.get_ident())
    )
    try:
        with temporary.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(str(temporary), str(path))
        except FileExistsError:
            if path.read_bytes() != data:
                die("refusing to replace immutable artifact {}".format(path))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def json_lines(values: Iterable[object]) -> bytes:
    return b"".join(canonical_json(value) for value in values)


def strict_json_loads(raw: str) -> object:
    """Decode standards-only JSON while rejecting ambiguous object keys."""
    def reject_constant(token: str) -> object:
        raise ValueError("non-standard JSON constant {!r}".format(token))

    def parse_finite_float(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            raise ValueError("non-finite JSON number {!r}".format(token))
        return value

    def unique_object(pairs: Sequence[Tuple[str, object]]) -> Dict[str, object]:
        result: Dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key {!r}".format(key))
            result[key] = value
        return result

    return json.loads(
        raw,
        parse_constant=reject_constant,
        parse_float=parse_finite_float,
        object_pairs_hook=unique_object,
    )


def load_json(path: Path) -> Dict[str, object]:
    try:
        value = strict_json_loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, ValueError) as exc:
        die("cannot read {}: {}".format(path, exc))
    if not isinstance(value, dict):
        die("{} is not a JSON object".format(path))
    return value


def run_captured(
    argv: Sequence[str],
    *,
    cwd: Path,
    timeout: float,
    environment: Optional[Mapping[str, str]] = None,
    label: str,
) -> Tuple[bytes, bytes]:
    try:
        result = subprocess.run(
            list(argv), cwd=str(cwd), env=None if environment is None else dict(environment),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        die("{} failed to run: {}".format(label, exc))
    if result.returncode:
        detail = result.stderr[-4096:].decode("utf-8", "replace")
        die("{} exited {}: {}".format(label, result.returncode, detail))
    return result.stdout, result.stderr


def git_output(repo: Path, *args: str) -> bytes:
    stdout, stderr = run_captured(
        ("git", *args), cwd=repo, timeout=60.0, label="git " + " ".join(args),
    )
    if stderr:
        die("git {} produced stderr".format(" ".join(args)))
    return stdout


def clean_source_identity(repo: Path) -> Dict[str, object]:
    status = git_output(repo, "status", "--porcelain=v1", "--untracked-files=all")
    if status:
        die("campaign preparation requires a completely clean source worktree")
    head = git_output(repo, "rev-parse", "HEAD").decode("ascii").strip()
    tree = git_output(repo, "rev-parse", "HEAD^{tree}").decode("ascii").strip()
    if re.fullmatch(r"[0-9a-f]{40}", head) is None or \
            re.fullmatch(r"[0-9a-f]{40}", tree) is None:
        die("source commit/tree identity is malformed")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", CODEC_IMPLEMENTATION_HEAD, head],
        cwd=str(repo), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    if ancestor.returncode != 0:
        die("campaign source does not contain the reviewed row-mask implementation")
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", head)
    if not manifest or not manifest.endswith(b"\n"):
        die("source tree manifest is empty or noncanonical")
    return {
        "commit": head,
        "tree_oid": tree,
        "status_porcelain_sha256": sha256_bytes(status),
        "tree_manifest_sha256": sha256_bytes(manifest),
        "tree_manifest_bytes": len(manifest),
    }


def executable_identity(name: str, version_args: Sequence[str]) -> Dict[str, object]:
    located = shutil.which(name)
    if not located:
        die("required build tool is unavailable: {}".format(name))
    path = Path(located).resolve(strict=True)
    stdout, stderr = run_captured(
        (str(path), *version_args), cwd=Path.cwd(), timeout=30.0,
        label="{} version".format(name),
    )
    return {
        "path": str(path), "sha256": sha256_file(path),
        "version_stdout_sha256": sha256_bytes(stdout),
        "version_stderr_sha256": sha256_bytes(stderr),
    }


def parse_cmake_cache_value(cache: Path, key: str) -> str:
    prefix = key + ":"
    matches: List[str] = []
    for line in cache.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix) and "=" in line:
            matches.append(line.split("=", 1)[1])
    if len(matches) != 1:
        die("CMake cache lacks one {} entry".format(key))
    return matches[0]


def fresh_build_binary(repo: Path, frozen: Path, workers: int) -> Dict[str, object]:
    """Build the benchmark from a clean committed tree and freeze provenance."""
    before = clean_source_identity(repo)
    tools = {
        "cmake": executable_identity("cmake", ("--version",)),
        "ninja": executable_identity("ninja", ("--version",)),
        "readelf": executable_identity("readelf", ("--version",)),
    }
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", str(before["commit"]))
    archive = git_output(repo, "archive", "--format=tar", str(before["commit"]))
    source_record = dict(before)
    source_record["source_archive_sha256"] = sha256_bytes(archive)
    build_evidence = frozen / "build"
    make_private_directory(build_evidence)
    (build_evidence / "source_tree.txt").write_bytes(manifest)
    with tempfile.TemporaryDirectory(prefix="wirehair-allk-holdout-build-") as temporary:
        source = Path(temporary) / "source"
        make_private_directory(source)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as stream:
            stream.extractall(path=source, filter="data")
        build = Path(temporary) / "build"
        configure = [
            str(tools["cmake"]["path"]), "-S", str(source), "-B", str(build),
            *BUILD_CONFIGURE_ARGS,
        ]
        environment = dict(os.environ)
        environment.update({"LANG": "C", "LC_ALL": "C", "TZ": "UTC"})
        configure_stdout, configure_stderr = run_captured(
            configure, cwd=Path(temporary), timeout=600.0,
            environment=environment, label="fresh campaign CMake configure",
        )
        build_command = [
            str(tools["cmake"]["path"]), "--build", str(build),
            "--target", BUILD_TARGET, "--parallel", str(workers), "--verbose",
        ]
        build_stdout, build_stderr = run_captured(
            build_command, cwd=Path(temporary), timeout=1800.0,
            environment=environment, label="fresh campaign benchmark build",
        )
        candidate = build / "codec" / BINARY_NAME
        cache = build / "CMakeCache.txt"
        graph = build / "build.ninja"
        compile_commands = build / "compile_commands.json"
        for path, label in ((candidate, "built binary"), (cache, "CMake cache"),
                            (graph, "Ninja graph"),
                            (compile_commands, "compile commands")):
            if not path.is_file() or path.is_symlink():
                die("fresh build lacks a regular {}".format(label))
        if not os.access(str(candidate), os.X_OK):
            die("fresh benchmark is not executable")
        if Path(parse_cmake_cache_value(cache, "CMAKE_HOME_DIRECTORY")).resolve() != source:
            die("fresh CMake cache is not bound to the campaign source")
        compiler = Path(parse_cmake_cache_value(cache, "CMAKE_CXX_COMPILER")).resolve(strict=True)
        compiler_stdout, compiler_stderr = run_captured(
            (str(compiler), "--version"), cwd=Path(temporary), timeout=30.0,
            label="C++ compiler version",
        )
        notes_stdout, notes_stderr = run_captured(
            (str(tools["readelf"]["path"]), "--notes", str(candidate)),
            cwd=Path(temporary), timeout=30.0, label="ELF build-ID inspection",
        )
        build_ids = re.findall(
            rb"^[ \t]*Build ID: ([0-9a-f]+)[ \t]*$", notes_stdout, re.MULTILINE,
        )
        if notes_stderr or len(build_ids) != 1:
            die("fresh benchmark lacks exactly one lowercase ELF build ID")
        evidence_sources = {
            "CMakeCache.txt": cache,
            "build.ninja": graph,
            "compile_commands.json": compile_commands,
        }
        for name, evidence_source in evidence_sources.items():
            shutil.copyfile(str(evidence_source), str(build_evidence / name))
        log_bytes = {
            "configure.stdout": configure_stdout,
            "configure.stderr": configure_stderr,
            "build.stdout": build_stdout,
            "build.stderr": build_stderr,
            "compiler.stdout": compiler_stdout,
            "compiler.stderr": compiler_stderr,
            "readelf-notes.stdout": notes_stdout,
            "readelf-notes.stderr": notes_stderr,
        }
        for name, data in log_bytes.items():
            (build_evidence / name).write_bytes(data)
        shutil.copyfile(str(candidate), str(frozen / BINARY_NAME))
        os.chmod(str(frozen / BINARY_NAME), 0o555)

    after = clean_source_identity(repo)
    if after != before:
        die("source commit/tree/worktree changed during the fresh build")
    evidence_files = sorted(path for path in build_evidence.iterdir() if path.is_file())
    record = {
        "schema": SCHEMA + ".fresh_build",
        "source": source_record,
        "configure_command": configure,
        "build_command": build_command,
        "workers": workers,
        "tools": tools,
        "compiler": {
            "path": str(compiler), "sha256": sha256_file(compiler),
            "version_stdout_sha256": sha256_bytes(compiler_stdout),
            "version_stderr_sha256": sha256_bytes(compiler_stderr),
        },
        "binary_sha256": sha256_file(frozen / BINARY_NAME),
        "binary_elf_build_id": build_ids[0].decode("ascii"),
        "evidence": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in evidence_files
        },
        "fresh_build": True,
    }
    write_once(frozen / "build_receipt.json", canonical_json(record))
    return record


def cell_tuple(record: Mapping[str, object]) -> Tuple[int, int, str]:
    try:
        key = (int(record["K"]), int(record["seed_index"]), str(record["schedule"]))
    except (KeyError, TypeError, ValueError) as exc:
        die("malformed all-K failure record: {}".format(exc))
    if not (K_LO <= key[0] <= K_HI and 0 <= key[1] < len(SEEDS) and
            key[2] in SCHEDULES):
        die("out-of-domain all-K failure key {}".format(key))
    return key


def cell_record(key: Tuple[int, int, str]) -> Dict[str, object]:
    return {"K": key[0], "seed_index": key[1], "schedule": key[2]}


def validate_source_summary(summary: Mapping[str, object]) -> Dict[str, List[dict]]:
    exact = {
        "schema": 3,
        "head": SOURCE_HEAD,
        "K_range": [K_LO, K_HI],
        "K_count": K_HI - K_LO + 1,
        "cells_per_arm": (K_HI - K_LO + 1) * len(SEEDS) * len(SCHEDULES),
        "validation_issue_count": 0,
        "timing_promotional": False,
    }
    for key, expected in exact.items():
        if summary.get(key) != expected:
            die("all-K source receipt {} changed".format(key))
    arms = summary.get("arms")
    if not isinstance(arms, dict) or set(arms) != set(SOURCE_ARMS):
        die("all-K source arm set changed")
    records: Dict[str, List[dict]] = {}
    for arm in SOURCE_ARMS:
        value = arms.get(arm)
        if not isinstance(value, dict) or not isinstance(value.get("failure_records"), list):
            die("all-K source arm {} lacks failure records".format(arm))
        seen: Set[Tuple[int, int, str]] = set()
        parsed: List[dict] = []
        for raw in value["failure_records"]:
            if not isinstance(raw, dict):
                die("all-K source failure record is not an object")
            key = cell_tuple(raw)
            if key in seen:
                die("duplicate all-K source failure {} {}".format(arm, key))
            seen.add(key)
            cause = raw.get("cause")
            if cause not in ("field_shortfall", "q>H"):
                die("unknown all-K failure cause {}".format(cause))
            if cause == "field_shortfall" and (
                    raw.get("binary_def") != 12 or raw.get("heavy_gain") != 11):
                die("field-shortfall receipt changed for {} {}".format(arm, key))
            parsed.append(dict(raw))
        if int(value.get("failures", -1)) != len(parsed):
            die("all-K failure cardinality changed for {}".format(arm))
        if len(parsed) != SOURCE_ARM_FAILURES[arm]:
            die("all-K sealed failure total changed for {}".format(arm))
        records[arm] = parsed
    for arm, count in SOURCE_FIELD_SHORTFALLS.items():
        if sum(row["cause"] == "field_shortfall" for row in records[arm]) != count:
            die("all-K finalist field-shortfall count changed for {}".format(arm))
    return records


def derive_hard_keys(records: Mapping[str, Sequence[dict]]) -> List[dict]:
    sources: Dict[Tuple[int, int, str], List[str]] = defaultdict(list)
    for arm in FINALIST_ARMS:
        for row in records[arm]:
            if row["cause"] == "field_shortfall":
                sources[cell_tuple(row)].append(arm)
    if len(sources) != HARD_CELL_COUNT:
        die("expected exact 16-cell finalist field-shortfall union")
    return [
        {**cell_record(key), "source_arms": sorted(sources[key])}
        for key in sorted(sources)
    ]


def validate_arms(arms: Sequence[Mapping[str, object]]) -> List[dict]:
    """Validate the sealed arm table and return enriched copies."""
    if len(arms) != 4:
        die("arm catalog must contain exactly four arms")
    by_name: Dict[str, dict] = {}
    result: List[dict] = []
    for arm in arms:
        name = arm.get("name")
        if not isinstance(name, str) or not name or name in by_name:
            die("arm catalog contains a missing or duplicate name")
        period = arm.get("period")
        rows = arm.get("rows")
        mask = arm.get("mask")
        if period not in (48, 32) or not isinstance(rows, int) or \
                isinstance(rows, bool) or not isinstance(mask, int) or \
                isinstance(mask, bool):
            die("arm {} has a malformed period/rows/mask".format(name))
        if not 0 < mask < (1 << 10):
            die("arm {} mask is outside the low 10 bits".format(name))
        if bin(mask).count("1") != rows:
            die("arm {} mask popcount does not equal its row count".format(name))
        if arm.get("panel") not in FINALIST_ARMS:
            die("arm {} references an unknown finalist panel".format(name))
        if arm.get("role") not in ("candidate", "comparator"):
            die("arm {} has an unknown role".format(name))
        enriched = dict(arm)
        enriched["mask_hex"] = "0x{:03x}".format(mask)
        enriched["canonical_suffix"] = mask == CANONICAL_SUFFIX_MASKS[period]
        by_name[name] = enriched
        result.append(enriched)
    candidates = [arm for arm in result if arm["role"] == "candidate"]
    comparators = [arm for arm in result if arm["role"] == "comparator"]
    if len(candidates) != 2 or len(comparators) != 2:
        die("arm catalog must contain two candidates and two comparators")
    for arm in comparators:
        if arm.get("comparator") is not None:
            die("comparator arm {} must not name a comparator".format(arm["name"]))
        if not arm["canonical_suffix"]:
            die("comparator arm {} is not the canonical suffix".format(arm["name"]))
    for arm in candidates:
        partner = by_name.get(str(arm.get("comparator")))
        if partner is None or partner["role"] != "comparator":
            die("candidate arm {} lacks a comparator arm".format(arm["name"]))
        if partner["panel"] != arm["panel"] or partner["period"] != arm["period"] \
                or partner["rows"] != arm["rows"]:
            die("candidate arm {} is not paired within its panel".format(arm["name"]))
        if arm["canonical_suffix"] or arm["mask"] == partner["mask"]:
            die("candidate arm {} duplicates the canonical suffix".format(arm["name"]))
    return result


def arms_catalog() -> List[dict]:
    return validate_arms(ARMS)


def hard_strata(hard: Sequence[Mapping[str, object]]) -> Dict[Tuple[int, str], List[int]]:
    strata: Dict[Tuple[int, str], List[int]] = defaultdict(list)
    seen: Set[Tuple[int, int, str]] = set()
    for record in hard:
        key = cell_tuple(record)
        if key in seen:
            die("duplicate development hard cell {}".format(key))
        seen.add(key)
        strata[(key[1], key[2])].append(key[0])
    return {key: sorted(values) for key, values in strata.items()}


def chunked(values: Sequence[int], size: int = CHUNK_MAX) -> Iterator[List[int]]:
    if size <= 0:
        die("chunk size must be positive")
    for offset in range(0, len(values), size):
        yield list(values[offset:offset + size])


def build_tasks(
    binary: str,
    hard: Sequence[Mapping[str, object]],
    k_lo: int = K_LO,
    k_hi: int = K_HI,
) -> List[dict]:
    arms = arms_catalog()
    strata = hard_strata(hard)
    tasks: List[dict] = []
    common = [
        binary, "precodefail", "--bb-list", "64", "--overhead", "0",
        "--trials", "1", "--threads", "1", "--loss", "0.50",
        "--completion", "mixed", "--mix-count", "2",
        "--packet-peel-seed-xor", "0",
    ]
    # Seal the hard fix-confirmation phase before the all-K control phase.
    # The launcher also places a process-pool barrier between these phases,
    # so no control task can start while a hard task is still running.
    for stage in ("hard", "control"):
        for arm in arms:
            for seed_index in range(len(SEEDS)):
                for schedule in SCHEDULES:
                    if stage == "hard":
                        Ks = strata.get((seed_index, schedule), [])
                    else:
                        Ks = list(range(k_lo, k_hi + 1))
                    for chunk_index, chunk in enumerate(chunked(Ks)):
                        job = len(tasks)
                        name = (
                            "{:05d}.{}.{}.seed{}.{}.chunk{:03d}.csv"
                            .format(job, arm["name"], stage, seed_index,
                                    schedule, chunk_index)
                        )
                        task = {
                            "job": job, "arm": arm["name"],
                            "panel": arm["panel"], "role": arm["role"],
                            "period": arm["period"], "rows": arm["rows"],
                            "mask": arm["mask"], "mask_hex": arm["mask_hex"],
                            "canonical_suffix": arm["canonical_suffix"],
                            "stage": stage, "seed_index": seed_index,
                            "seed": SEEDS[seed_index], "schedule": schedule,
                            "chunk_index": chunk_index, "Ks": chunk,
                            "output_name": name,
                        }
                        task["argv"] = [
                            *common, "--N", ",".join(map(str, chunk)),
                            "--seed", SEEDS[seed_index], "--schedule", schedule,
                            "--mixed-period", str(arm["period"]),
                            "--mixed-geometry", "shared-x",
                            "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
                            "--mixed-grouped-gf256-rows", str(arm["rows"]),
                            "--mixed-grouped-gf256-row-mask",
                            "0x{:x}".format(int(arm["mask"])),
                            "--mixed-residue-buckets", "auto",
                            "--binary-dense-two-anchor",
                        ]
                        tasks.append(task)
    k_count = k_hi - k_lo + 1
    expected_cells = len(arms) * (
        len(hard) + k_count * len(SEEDS) * len(SCHEDULES))
    if sum(len(task["Ks"]) for task in tasks) != expected_cells:
        die("task ledger does not cover the exact Cartesian product")
    if any(not task["Ks"] or len(task["Ks"]) > CHUNK_MAX for task in tasks):
        die("empty or oversized task generated")
    if [task["job"] for task in tasks] != list(range(len(tasks))) or \
            len({task["output_name"] for task in tasks}) != len(tasks):
        die("task ledger identity is not canonical")
    return tasks


def plan_from_summary(summary: Mapping[str, object], binary: str) -> Dict[str, object]:
    records = validate_source_summary(summary)
    hard = derive_hard_keys(records)
    tasks = build_tasks(binary, hard)
    return {
        "hard": hard, "arms": arms_catalog(), "tasks": tasks,
        "source_records": records,
    }


def expected_preamble(task: Mapping[str, object]) -> Dict[str, str]:
    return {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": SEEDS[int(task["seed_index"])], "completion": "mixed",
        "mixed_period": str(task["period"]), "mixed_gf256_rows": "10",
        "mixed_gf16_rows": "2", "mixed_geometry": "shared-x",
        "mixed_residue_skew": "0", "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0", "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": str(task["rows"]),
        "mixed_grouped_gf256_row_mask": "0x{:x}".format(int(task["mask"])),
        "mixed_grouped_gf256_hash_seed": "0xb7e15162",
        "mixed_grouped_final_h_a_columns": "12",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "binary_dense_two_anchor": "1", "binary_dense_two_anchor_phase": "0",
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1", "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "salted",
        "full_payload_solve": "0", "schedule": str(task["schedule"]),
    }


def parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die("missing precodefail preamble")
    result: Dict[str, str] = {}
    for token in line[len(prefix):].split():
        if "=" not in token:
            die("malformed preamble token")
        key, value = token.split("=", 1)
        if key in result:
            die("duplicate preamble token {}".format(key))
        result[key] = value
    return result


def decimal_field(row: Mapping[str, str], key: str) -> Decimal:
    try:
        value = Decimal(row[key])
    except (KeyError, InvalidOperation, TypeError) as exc:
        die("invalid decimal field {}: {}".format(key, exc))
    if not value.is_finite():
        die("nonfinite decimal field {}".format(key))
    return value


def benchmark_uint(row: Mapping[str, str], key: str) -> int:
    text = row.get(key)
    if not isinstance(text, str) or \
            re.fullmatch(r"0|[1-9][0-9]*", text) is None:
        die("invalid unsigned benchmark field {}".format(key))
    value = int(text)
    if value >= 1 << 64:
        die("unsigned benchmark field {} exceeds u64".format(key))
    return value


def benchmark_milli(row: Mapping[str, str], key: str) -> Decimal:
    text = row.get(key)
    if not isinstance(text, str) or \
            re.fullmatch(r"(?:0|[1-9][0-9]*)\.[0-9]{3}", text) is None:
        die("invalid fixed-milli benchmark field {}".format(key))
    value = Decimal(text)
    scaled = value * 1000
    if scaled >= 1 << 64:
        die("fixed-milli benchmark field {} exceeds scaled u64".format(key))
    return value


def benchmark_counter_mu(row: Mapping[str, str], key: str) -> int:
    text = row.get(key)
    if not isinstance(text, str) or \
            re.fullmatch(r"(?:0|[1-9][0-9]*)\.000", text) is None:
        die("invalid single-trial counter mean {}".format(key))
    value = int(text[:-4])
    if value >= 1 << 64:
        die("single-trial counter mean {} exceeds u64".format(key))
    return value


def validate_single_trial_semantics(
    row: Mapping[str, str],
) -> Tuple[int, int, int]:
    success, rank_fail, error = (
        benchmark_uint(row, key)
        for key in ("success", "rank_fail", "error")
    )
    if any(value not in (0, 1) for value in (
            success, rank_fail, error)) or \
            success + rank_fail + error != 1 or error:
        die("invalid benchmark outcome")

    inactivated = benchmark_uint(row, "inact_max")
    binary_deficit = benchmark_uint(row, "binary_def_max")
    heavy_gain = benchmark_uint(row, "heavy_gain_min")
    heavy_shortfall = benchmark_uint(row, "heavy_shortfall")
    if any(value >= 1 << 32 for value in (
            inactivated, binary_deficit, heavy_gain)) or \
            heavy_shortfall > 1 or binary_deficit > inactivated or \
            heavy_gain > binary_deficit or \
            (success and (
                binary_deficit > 12 or heavy_gain != binary_deficit)) or \
            (rank_fail and (
                (binary_deficit <= 12 and heavy_gain >= binary_deficit) or
                (binary_deficit > 12 and heavy_gain != 0))):
        die("single-trial rank/deficit/gain semantics mismatch")
    expected_shortfall = (
        1 if rank_fail and binary_deficit <= 12 and
        heavy_gain < binary_deficit else 0)
    expected = {
        "fail_rate": "1.00000000" if rank_fail else "0.00000000",
        "inact_mu": "{}.000".format(inactivated),
        "binary_def_mu": "{}.000".format(binary_deficit),
        "heavy_gain_mu": "{}.000".format(heavy_gain),
        "heavy_shortfall": str(expected_shortfall),
        "first_rank_fail": "0" if rank_fail else "-1",
        "binary_def_hist": "{}:1".format(binary_deficit),
        "heavy_gain_hist": "{}:1".format(heavy_gain),
        "failure_trials": "0" if rank_fail else "",
    }
    if any(row.get(key) != value for key, value in expected.items()):
        die("single-trial aggregate semantics mismatch")
    return success, rank_fail, error


def classify(row: Mapping[str, str]) -> str:
    success, _, _ = validate_single_trial_semantics(row)
    if success:
        return "success"
    return (
        "q>H" if benchmark_uint(row, "binary_def_max") > 12
        else "field_shortfall")


def parse_benchmark_csv(data: bytes, task: Mapping[str, object]) -> List[Dict[str, str]]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        die("benchmark output is not ASCII: {}".format(exc))
    lines = text.splitlines(keepends=True)
    if not lines or not all(line.endswith("\n") for line in lines) or \
            any("\r" in line or line == "\n" for line in lines):
        die("benchmark output is empty or lacks final LF")
    preamble = parse_preamble(lines[0][:-1])
    expected = expected_preamble(task)
    if set(preamble) != set(expected) or preamble != expected:
        die("benchmark preamble differs from the sealed task")
    reader = csv.DictReader(lines[1:])
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        die("benchmark CSV header mismatch")
    rows = list(reader)
    if any(
            None in row or tuple(row) != CSV_FIELDS or
            any(value is None for value in row.values())
            for row in rows):
        die("benchmark CSV row width changed")
    observed_Ks = [benchmark_uint(row, "N") for row in rows]
    if observed_Ks != list(task["Ks"]):
        die("benchmark CSV K ledger mismatch")
    fixed = {
        "bb": "64", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "0", "trials": "1", "active_packet_peel_seed_xor": "0x0",
    }
    for row in rows:
        if any(row.get(key) != value for key, value in fixed.items()):
            die("benchmark fixed CSV field mismatch")
        if benchmark_uint(row, "seed_attempt") >= 1 << 8:
            die("benchmark seed attempt exceeds the codec u8 receipt")
        for key in CSV_TIMING_MILLI_FIELDS:
            benchmark_milli(row, key)
        for key in CSV_COUNTER_MU_FIELDS:
            benchmark_counter_mu(row, key)
        classify(row)
    return rows


def smoke_frozen_binary(binary: Path, frozen: Path) -> Dict[str, object]:
    smoke_dir = frozen / "smoke"
    make_private_directory(smoke_dir)
    records: List[dict] = []
    for index, (period, rows, mask) in enumerate(SMOKE_CASES):
        task = {
            "period": period, "rows": rows, "mask": mask,
            "seed_index": 0, "schedule": "adversarial", "Ks": [257],
        }
        argv = [
            str(binary), "precodefail", "--bb-list", "64", "--overhead", "0",
            "--trials", "1", "--threads", "1", "--loss", "0.50",
            "--completion", "mixed", "--mix-count", "2",
            "--packet-peel-seed-xor", "0", "--N", "257",
            "--seed", SEEDS[0], "--schedule", "adversarial",
            "--mixed-period", str(period), "--mixed-geometry", "shared-x",
            "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
            "--mixed-grouped-gf256-rows", str(rows),
            "--mixed-grouped-gf256-row-mask", "0x{:x}".format(mask),
            "--mixed-residue-buckets", "auto", "--binary-dense-two-anchor",
        ]
        stdout, stderr = run_captured(
            argv, cwd=frozen, timeout=60.0,
            label="bounded holdout smoke P{} mask 0x{:x}".format(period, mask),
        )
        parsed = parse_benchmark_csv(stdout, task)
        if stderr:
            die("bounded holdout smoke produced stderr")
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        stdout_path = smoke_dir / (stem + ".stdout")
        stderr_path = smoke_dir / (stem + ".stderr")
        stdout_path.write_bytes(stdout)
        stderr_path.write_bytes(stderr)
        records.append({
            "case": index, "period": period, "rows": rows,
            "mask": mask, "mask_hex": "0x{:03x}".format(mask), "K": 257,
            "argv": ["frozen/" + BINARY_NAME, *argv[1:]],
            "stdout": str(stdout_path.relative_to(frozen.parent)),
            "stdout_sha256": sha256_bytes(stdout),
            "stderr": str(stderr_path.relative_to(frozen.parent)),
            "stderr_sha256": sha256_bytes(stderr),
            "outcome": classify(parsed[0]),
        })
    record = {
        "schema": SCHEMA + ".binary_smoke", "bounded": True,
        "binary_sha256": sha256_file(binary), "cases": records,
    }
    write_once(frozen / "binary_smoke.json", canonical_json(record))
    return record


def artifact_map(root: Path) -> Dict[str, Path]:
    result = {
        "controller": root / "frozen" / CONTROLLER_NAME,
        "binary": root / "frozen" / BINARY_NAME,
        "sampler": root / "frozen" / SAMPLER_NAME,
        "source_summary": root / "frozen" / SOURCE_SUMMARY_NAME,
        "build_receipt": root / "frozen" / "build_receipt.json",
        "binary_smoke": root / "frozen" / "binary_smoke.json",
        "design": root / "design.json",
        "hard_keys": root / "hard_keys.jsonl",
        "arms": root / "arms.jsonl",
        "tasks": root / "tasks_manifest.jsonl",
    }
    for name in (
        "source_tree.txt", "CMakeCache.txt", "build.ninja",
        "compile_commands.json", "configure.stdout", "configure.stderr",
        "build.stdout", "build.stderr", "compiler.stdout", "compiler.stderr",
        "readelf-notes.stdout", "readelf-notes.stderr",
    ):
        result["build_" + name] = root / "frozen" / "build" / name
    for index, (period, rows, mask) in enumerate(SMOKE_CASES):
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        result["smoke_{}_stdout".format(index)] = \
            root / "frozen" / "smoke" / (stem + ".stdout")
        result["smoke_{}_stderr".format(index)] = \
            root / "frozen" / "smoke" / (stem + ".stderr")
    return result


def owner_time(value: datetime) -> str:
    """Format a UTC instant exactly as the janitor's expiry comparator."""
    if value.tzinfo is None or value.utcoffset() != timedelta(0):
        die("environment-owner timestamps must be UTC")
    return value.strftime("%Y-%m-%dT%H:%M:%S+00:00")


def read_boot_id(boot_id_path: Path = BOOT_ID_PATH) -> Optional[str]:
    try:
        return boot_id_path.read_text(encoding="ascii").strip()
    except OSError:
        return None


def valid_segment_boot_id(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(
        r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-"
        r"[0-9a-f]{4}-[0-9a-f]{12}",
        value,
    ) is not None


def required_boot_id() -> str:
    value = read_boot_id()
    if not valid_segment_boot_id(value):
        die("cannot bind campaign process identities to the current boot")
    assert value is not None
    return value


def owner_protected_entry(role: str, pid: int) -> Dict[str, object]:
    start_ticks = process_start_ticks(pid)
    if start_ticks is None:
        die("cannot bind protected PID {} to its start time".format(pid))
    return {"role": str(role), "pid": int(pid), "start_ticks": int(start_ticks)}


def validate_owner_protected(
    protected: object,
    *,
    require_nonempty: bool,
) -> List[dict]:
    if not isinstance(protected, (list, tuple)) or \
            (require_nonempty and not protected):
        die("environment-owner marker has invalid protected identities")
    entries: List[dict] = []
    seen: Set[Tuple[int, int]] = set()
    for entry in protected:
        if not isinstance(entry, Mapping) or \
                set(entry) != {"role", "pid", "start_ticks"} or \
                not isinstance(entry.get("role"), str) or not entry["role"] or \
                not isinstance(entry.get("pid"), int) or \
                isinstance(entry.get("pid"), bool) or int(entry["pid"]) <= 1 or \
                int(entry["pid"]) > MAX_PROCESS_PID or \
                not isinstance(entry.get("start_ticks"), int) or \
                isinstance(entry.get("start_ticks"), bool) or \
                int(entry["start_ticks"]) < 0:
            die("malformed environment-owner protected entry")
        identity = (int(entry["pid"]), int(entry["start_ticks"]))
        if identity in seen:
            die("environment-owner marker repeats a protected identity")
        seen.add(identity)
        entries.append({
            "role": str(entry["role"]),
            "pid": identity[0],
            "start_ticks": identity[1],
        })
    return entries


def write_owner_marker(
    root: Path,
    phase: str,
    protected: Sequence[Mapping[str, object]],
    ttl_hours: int,
    *,
    marker_path: Path = OWNER_MARKER_PATH,
    boot_id_path: Path = BOOT_ID_PATH,
    now: Optional[datetime] = None,
) -> Dict[str, object]:
    """Atomically publish the wirehair.environment_owner.v1 marker."""
    if phase not in ("executing", "complete"):
        die("unknown environment-owner phase {}".format(phase))
    if not isinstance(ttl_hours, int) or isinstance(ttl_hours, bool) or ttl_hours <= 0:
        die("environment-owner TTL must be a positive integer hour count")
    entries = validate_owner_protected(
        list(protected), require_nonempty=phase == "executing")
    created = datetime.now(timezone.utc) if now is None else now
    marker = {
        "schema": OWNER_MARKER_SCHEMA,
        "campaign_root": str(root),
        "phase": phase,
        "boot_id": read_boot_id(boot_id_path),
        "created_utc": owner_time(created),
        "expires_utc": owner_time(created + timedelta(hours=ttl_hours)),
        "protected": entries,
    }
    scratch = marker_path.with_name(
        marker_path.name + ".part.{}.{}".format(
            os.getpid(), threading.get_ident()))
    scratch_created = False
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | \
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(str(scratch), flags, 0o644)
        scratch_created = True
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(canonical_json(marker))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(str(scratch), str(marker_path))
    except OSError as exc:
        die("cannot publish environment-owner marker: {}".format(exc))
    finally:
        if scratch_created:
            try:
                scratch.unlink()
            except FileNotFoundError:
                pass
    return marker


def load_active_owner_marker(
    *,
    marker_path: Path = OWNER_MARKER_PATH,
    boot_id_path: Path = BOOT_ID_PATH,
    now: Optional[datetime] = None,
) -> Optional[Dict[str, object]]:
    """Return the active ownership marker, mirroring janitor inertness rules.

    A present-but-unreadable or malformed marker fails closed with an error:
    this controller must never launch into an environment whose ownership it
    cannot prove.  Markers that are complete, expired, or from a previous
    boot are inert and return None.
    """
    try:
        identity = marker_path.lstat()
    except FileNotFoundError:
        return None
    except OSError as exc:
        die("environment-owner marker cannot be inspected: {}".format(exc))
    if not stat.S_ISREG(identity.st_mode) or identity.st_uid != os.getuid() or \
            identity.st_nlink != 1 or identity.st_mode & 0o022:
        die("environment-owner marker is not a safe owned regular file")
    try:
        raw = marker_path.read_text(encoding="ascii")
    except (OSError, UnicodeError) as exc:
        die("environment-owner marker is unreadable: {}".format(exc))
    try:
        marker = strict_json_loads(raw)
    except ValueError as exc:
        die("environment-owner marker is malformed JSON: {}".format(exc))
    if not isinstance(marker, dict) or \
            marker.get("schema") != OWNER_MARKER_SCHEMA:
        die("environment-owner marker has an unknown schema")
    if marker.get("phase") not in ("executing", "complete") or \
            not isinstance(marker.get("campaign_root"), str) or \
            not str(marker["campaign_root"]).startswith("/") or \
            not isinstance(marker.get("protected"), list) or \
            marker.get("boot_id") is not None and (
                not isinstance(marker.get("boot_id"), str) or
                not marker["boot_id"]):
        die("environment-owner marker has malformed ownership fields")
    marker["protected"] = validate_owner_protected(
        marker["protected"], require_nonempty=marker["phase"] == "executing")
    created = marker.get("created_utc")
    expires = marker.get("expires_utc")
    timestamp_pattern = (
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T"
        r"[0-9]{2}:[0-9]{2}:[0-9]{2}\+00:00")
    if not isinstance(created, str) or not isinstance(expires, str) or \
            re.fullmatch(timestamp_pattern, created) is None or \
            re.fullmatch(timestamp_pattern, expires) is None:
        die("environment-owner marker has malformed timestamps")
    try:
        created_time = datetime.strptime(created, "%Y-%m-%dT%H:%M:%S%z")
        expires_time = datetime.strptime(expires, "%Y-%m-%dT%H:%M:%S%z")
    except ValueError as exc:
        die("environment-owner marker has invalid timestamps: {}".format(exc))
    if created_time.utcoffset() != timedelta(0) or \
            expires_time.utcoffset() != timedelta(0) or \
            expires_time < created_time:
        die("environment-owner marker timestamp interval is invalid")
    if marker["phase"] == "complete":
        return None
    boot_id = read_boot_id(boot_id_path)
    if boot_id is not None and marker.get("boot_id") not in (None, boot_id):
        return None
    current = owner_time(datetime.now(timezone.utc) if now is None else now)
    if expires < current:
        return None
    return marker


def acquire_environment_lock(
    lock_path: Path = OWNER_LOCK_PATH,
) -> Optional[int]:
    """Acquire the cross-controller launch mutex without following links."""
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(str(lock_path), flags, 0o600)
    except OSError as exc:
        die("cannot open environment-owner lock: {}".format(exc))
    try:
        identity = os.fstat(descriptor)
        path_identity = lock_path.lstat()
        if not stat.S_ISREG(identity.st_mode) or identity.st_nlink != 1 or \
                identity.st_uid != os.getuid() or identity.st_mode & 0o077 or \
                (identity.st_dev, identity.st_ino) != \
                (path_identity.st_dev, path_identity.st_ino):
            die("environment-owner lock is not a safe regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            os.close(descriptor)
            return None
        except OSError as exc:
            die("cannot acquire environment-owner lock: {}".format(exc))
        return descriptor
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def command_prepare(args: argparse.Namespace) -> None:
    root = Path(args.output_root).resolve()
    source_path = Path(args.allk_root).resolve() / SOURCE_SUMMARY_NAME
    sampler = Path(args.sampler).resolve()
    controller = Path(__file__).resolve()
    repo = controller.parent.parent.resolve()
    if root.exists():
        die("output root already exists: {}".format(root))
    staging = root.parent / (root.name + ".prepare.{}".format(os.getpid()))
    if staging.exists():
        die("preparation staging path already exists: {}".format(staging))
    if not source_path.is_file() or source_path.is_symlink() or \
            sha256_file(source_path) != SOURCE_SUMMARY_SHA256:
        die("all-K validated-summary hash mismatch")
    for path, label in ((sampler, "sampler"), (controller, "controller")):
        if not path.is_file() or path.is_symlink():
            die("{} is not a regular file: {}".format(label, path))
    if sha256_file(sampler) != SAMPLER_SHA256:
        die("thermal sampler hash does not match the hardened pinned worker")
    if not isinstance(args.owner_ttl_hours, int) or args.owner_ttl_hours <= 0:
        die("owner TTL must be a positive integer hour count")
    try:
        make_private_directory(staging, parents=True)
        for directory in (
                "frozen", "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            make_private_directory(staging / directory)
        frozen = staging / "frozen"
        shutil.copyfile(str(controller), str(frozen / CONTROLLER_NAME))
        shutil.copyfile(str(sampler), str(frozen / SAMPLER_NAME))
        shutil.copyfile(str(source_path), str(frozen / SOURCE_SUMMARY_NAME))
        if sha256_file(frozen / SOURCE_SUMMARY_NAME) != SOURCE_SUMMARY_SHA256:
            die("copied all-K validated summary changed during preparation")
        for name in (CONTROLLER_NAME, SAMPLER_NAME):
            os.chmod(str(frozen / name), 0o555)

        build = fresh_build_binary(repo, frozen, int(args.workers))
        smoke = smoke_frozen_binary(frozen / BINARY_NAME, frozen)
        if smoke.get("binary_sha256") != build.get("binary_sha256"):
            die("binary smoke is not bound to the fresh build")
        summary = load_json(frozen / SOURCE_SUMMARY_NAME)
        # The task argv is sealed to the final root, not the transient staging
        # directory that is atomically renamed after every receipt is written.
        frozen_binary = str(root / "frozen" / BINARY_NAME)
        plan = plan_from_summary(summary, frozen_binary)
        write_once(staging / "hard_keys.jsonl", json_lines(plan["hard"]))
        write_once(staging / "arms.jsonl", json_lines(plan["arms"]))
        write_once(staging / "tasks_manifest.jsonl", json_lines(plan["tasks"]))

        hard_task_count = sum(
            task["stage"] == "hard" for task in plan["tasks"])
        design = {
            "schema": SCHEMA + ".design", "root": str(root),
            "git_head": build["source"]["commit"],
            "git_tree_oid": build["source"]["tree_oid"],
            "build_receipt_sha256": sha256_file(frozen / "build_receipt.json"),
            "binary_sha256": build["binary_sha256"],
            "sampler_sha256": SAMPLER_SHA256,
            "binary_smoke_sha256": sha256_file(frozen / "binary_smoke.json"),
            "codec_implementation_head": CODEC_IMPLEMENTATION_HEAD,
            "source_head": SOURCE_HEAD,
            "source_summary_sha256": SOURCE_SUMMARY_SHA256,
            "source_arm_failures": dict(SOURCE_ARM_FAILURES),
            "arms": plan["arms"],
            "holdout_policy": (
                "exact-comparison holdout on the sealed all-K seeds and "
                "schedules; no new independent seeds"),
            "dev_hard_cell_policy": (
                "the 16 development hard cells were selection material; they "
                "are reported separately as fix-confirmation and excluded "
                "from the headline holdout comparison"),
            "acceptance_gate": (
                "candidate raw field-shortfall count strictly below its "
                "comparator and raw total failures strictly below its "
                "comparator; on headline cells, total failures strictly "
                "below the comparator with at least one actual repair and "
                "zero introductions; paired XOR/muladd cost ratios on common "
                "headline successes within the sealed tolerance; at least "
                "one panel source field-shortfall genuinely fixed; zero "
                "codec errors"),
            "cost_ratio_max": str(COST_RATIO_MAX),
            "hard_cells": len(plan["hard"]),
            "seeds": list(SEEDS), "schedules": list(SCHEDULES),
            "K_range": [K_LO, K_HI], "bb": 64, "loss": 0.5,
            "overhead": 0, "trials": 1, "mix_count": 2,
            "seed_fixes": "none", "binary_dense_two_anchor": True,
            "cells_per_arm": (
                len(plan["hard"]) +
                (K_HI - K_LO + 1) * len(SEEDS) * len(SCHEDULES)),
            "total_cells": sum(len(task["Ks"]) for task in plan["tasks"]),
            "task_count": len(plan["tasks"]),
            "hard_task_count": hard_task_count,
            "chunk_max": CHUNK_MAX,
            "job_timeout_seconds": JOB_TIMEOUT_SECONDS,
            "benchmark_parent_death_signal": "SIGKILL",
            "thermal_controller_supervision":
                "pidfd+start_ticks+boot_id+ancestry",
            "stage_order": ["hard", "control"],
            "worker_count": args.workers, "thermal_single_sampler_required": True,
            "thermal_cpu_busy_floor_pct": str(CPU_BUSY_FLOOR),
            "thermal_busy_floor_stages": list(CPU_BUSY_FLOOR_STAGES),
            "thermal_interval_seconds": str(THERMAL_INTERVAL_SECONDS),
            "thermal_min_gap_seconds": str(THERMAL_MIN_GAP_SECONDS),
            "thermal_max_gap_seconds": str(THERMAL_MAX_GAP_SECONDS),
            "thermal_ready_samples": THERMAL_READY_SAMPLES,
            "owner_marker_path": str(OWNER_MARKER_PATH),
            "owner_marker_schema": OWNER_MARKER_SCHEMA,
            "owner_ttl_hours": args.owner_ttl_hours,
            "retry_policy": "stage-atomic non-selective retry after failed/interrupted segment",
            "timing_promotional": False,
            "timing_policy": "NONPROMOTIONAL saturated rank-proxy holdout; recovery and work receipts only",
        }
        write_once(staging / "design.json", canonical_json(design))
        paths = artifact_map(staging)
        receipts = {
            "schema": SCHEMA + ".prelaunch", "artifacts": {
                key: {"path": str(path.relative_to(staging)), "sha256": sha256_file(path)}
                for key, path in sorted(paths.items())
            },
        }
        write_once(staging / "prelaunch_receipts.json", canonical_json(receipts))
        os.replace(str(staging), str(root))
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    print(json.dumps({
        "status": "PREPARED_NOT_LAUNCHED", "root": str(root),
        "prelaunch_receipts_sha256": sha256_file(root / "prelaunch_receipts.json"),
        "tasks": len(plan["tasks"]), "cells": design["total_cells"],
        "source_commit": design["git_head"], "source_tree": design["git_tree_oid"],
        "binary_sha256": design["binary_sha256"], "smoke_cases": len(smoke["cases"]),
    }, indent=2, sort_keys=True))


def read_jsonl(path: Path) -> List[dict]:
    result: List[dict] = []
    try:
        raw_lines = path.read_bytes().splitlines(keepends=True)
    except OSError as exc:
        die("cannot read {}: {}".format(path, exc))
    for index, raw in enumerate(raw_lines):
        if not raw.endswith(b"\n"):
            die("{} line {} lacks LF".format(path, index + 1))
        try:
            value = strict_json_loads(raw.decode("ascii"))
        except (UnicodeError, ValueError) as exc:
            die("bad JSONL {} line {}: {}".format(path, index + 1, exc))
        if not isinstance(value, dict) or canonical_json(value) != raw:
            die("noncanonical JSONL {} line {}".format(path, index + 1))
        result.append(value)
    return result


def verify_root(
    root: Path,
    expected_receipts: str,
    *,
    runtime_tolerant: bool = False,
) -> Tuple[dict, List[dict]]:
    verify_runtime_directories(root)
    for path, label in (
            (root / "frozen", "frozen"),
            (root / "frozen" / "build", "frozen build"),
            (root / "frozen" / "smoke", "frozen smoke")):
        verify_safe_directory(path, label)
    receipt_path = root / "prelaunch_receipts.json"
    if not re.fullmatch(r"[0-9a-f]{64}", expected_receipts or ""):
        die("expected receipt hash must be 64 lowercase hex digits")
    if sha256_file(receipt_path) != expected_receipts:
        die("external prelaunch receipt hash mismatch")
    receipts = load_json(receipt_path)
    if receipts.get("schema") != SCHEMA + ".prelaunch":
        die("prelaunch schema mismatch")
    artifacts = receipts.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != set(artifact_map(root)):
        die("prelaunch artifact set mismatch")
    for key, expected_path in artifact_map(root).items():
        record = artifacts[key]
        if not isinstance(record, dict) or record.get("path") != str(expected_path.relative_to(root)):
            die("artifact path mismatch: {}".format(key))
        if expected_path.is_symlink() or not expected_path.is_file():
            die("artifact is not a regular file: {}".format(key))
        if sha256_file(expected_path) != record.get("sha256"):
            die("artifact hash mismatch: {}".format(key))
    controller_record = artifacts["controller"]
    current_controller_path = Path(__file__)
    current_controller = current_controller_path.resolve()
    if current_controller_path.is_symlink() or not current_controller.is_file() or \
            sha256_file(current_controller) != controller_record.get("sha256"):
        die("executing controller differs from the sealed controller artifact")
    expected_frozen_files = {
        CONTROLLER_NAME, BINARY_NAME, SAMPLER_NAME, SOURCE_SUMMARY_NAME,
        "build_receipt.json", "binary_smoke.json",
    }
    frozen_entries = list((root / "frozen").iterdir())
    if {path.name for path in frozen_entries} != \
            expected_frozen_files | {"build", "smoke"} or any(
                path.is_symlink() or (
                    not path.is_file()
                    if path.name in expected_frozen_files else not path.is_dir())
                for path in frozen_entries):
        die("frozen top-level inventory mismatch")
    design = load_json(root / "design.json")
    if design.get("schema") != SCHEMA + ".design" or design.get("root") != str(root):
        die("design identity mismatch")
    if design.get("codec_implementation_head") != CODEC_IMPLEMENTATION_HEAD or \
            design.get("source_head") != SOURCE_HEAD or \
            design.get("source_summary_sha256") != SOURCE_SUMMARY_SHA256 or \
            design.get("sampler_sha256") != SAMPLER_SHA256 or \
            design.get("source_arm_failures") != dict(SOURCE_ARM_FAILURES) or \
            design.get("cost_ratio_max") != str(COST_RATIO_MAX) or \
            design.get("chunk_max") != CHUNK_MAX or \
            design.get("job_timeout_seconds") != JOB_TIMEOUT_SECONDS or \
            design.get("benchmark_parent_death_signal") != "SIGKILL" or \
            design.get("thermal_controller_supervision") != \
                "pidfd+start_ticks+boot_id+ancestry" or \
            design.get("K_range") != [K_LO, K_HI] or \
            design.get("hard_cells") != HARD_CELL_COUNT or \
            design.get("stage_order") != ["hard", "control"] or \
            design.get("seed_fixes") != "none" or \
            design.get("thermal_cpu_busy_floor_pct") != str(CPU_BUSY_FLOOR) or \
            design.get("thermal_busy_floor_stages") != list(CPU_BUSY_FLOOR_STAGES) or \
            design.get("thermal_interval_seconds") != str(THERMAL_INTERVAL_SECONDS) or \
            design.get("thermal_min_gap_seconds") != str(THERMAL_MIN_GAP_SECONDS) or \
            design.get("thermal_max_gap_seconds") != str(THERMAL_MAX_GAP_SECONDS) or \
            design.get("thermal_ready_samples") != THERMAL_READY_SAMPLES or \
            not isinstance(design.get("worker_count"), int) or \
            isinstance(design.get("worker_count"), bool) or \
            not MIN_CAMPAIGN_WORKERS <= design["worker_count"] <= \
                MAX_CAMPAIGN_WORKERS or \
            design.get("owner_marker_path") != str(OWNER_MARKER_PATH) or \
            design.get("owner_marker_schema") != OWNER_MARKER_SCHEMA or \
            not isinstance(design.get("owner_ttl_hours"), int) or \
            isinstance(design.get("owner_ttl_hours"), bool) or \
            design.get("owner_ttl_hours") <= 0 or \
            design.get("timing_promotional") is not False or \
            design.get("retry_policy") != \
                "stage-atomic non-selective retry after failed/interrupted segment":
        die("design trust anchor mismatch")
    if sha256_file(root / "frozen" / SAMPLER_NAME) != SAMPLER_SHA256:
        die("frozen thermal sampler does not match the pinned hardened worker")
    build = load_json(root / "frozen" / "build_receipt.json")
    smoke = load_json(root / "frozen" / "binary_smoke.json")
    source = build.get("source")
    if not isinstance(source, dict) or \
            build.get("schema") != SCHEMA + ".fresh_build" or \
            build.get("fresh_build") is not True or \
            build.get("binary_sha256") != sha256_file(root / "frozen" / BINARY_NAME) or \
            source.get("commit") != design.get("git_head") or \
            source.get("tree_oid") != design.get("git_tree_oid") or \
            design.get("binary_sha256") != build.get("binary_sha256") or \
            design.get("build_receipt_sha256") != sha256_file(
                root / "frozen" / "build_receipt.json"):
        die("fresh-build provenance binding mismatch")
    if re.fullmatch(r"[0-9a-f]{40}", str(source.get("commit", ""))) is None or \
            re.fullmatch(r"[0-9a-f]{40}", str(source.get("tree_oid", ""))) is None or \
            re.fullmatch(r"[0-9a-f]{64}", str(source.get("source_archive_sha256", ""))) is None or \
            source.get("status_porcelain_sha256") != sha256_bytes(b"") or \
            source.get("tree_manifest_sha256") != sha256_file(
                root / "frozen" / "build" / "source_tree.txt") or \
            source.get("tree_manifest_bytes") != (
                root / "frozen" / "build" / "source_tree.txt").stat().st_size:
        die("fresh-build clean source/tree receipt mismatch")
    if smoke.get("schema") != SCHEMA + ".binary_smoke" or \
            smoke.get("bounded") is not True or \
            smoke.get("binary_sha256") != build.get("binary_sha256") or \
            not isinstance(smoke.get("cases"), list) or \
            len(smoke["cases"]) != len(SMOKE_CASES) or \
            design.get("binary_smoke_sha256") != sha256_file(
                root / "frozen" / "binary_smoke.json"):
        die("bounded binary-smoke provenance binding mismatch")
    expected_build_names = {
        path.name for key, path in artifact_map(root).items()
        if key.startswith("build_") and key != "build_receipt"
    }
    expected_smoke_names = {
        path.name for key, path in artifact_map(root).items()
        if key.startswith("smoke_")
    }
    for directory, expected_names in (
            (root / "frozen" / "build", expected_build_names),
            (root / "frozen" / "smoke", expected_smoke_names)):
        entries = list(directory.iterdir())
        actual_names = {path.name for path in entries}
        if actual_names != expected_names or any(
                path.is_symlink() or not path.is_file() for path in entries):
            die("frozen evidence inventory mismatch: {}".format(directory.name))
    evidence = build.get("evidence")
    if not isinstance(evidence, dict) or set(evidence) != expected_build_names:
        die("fresh-build evidence receipt set mismatch")
    for name in expected_build_names:
        path = root / "frozen" / "build" / name
        record = evidence[name]
        if not isinstance(record, dict) or record.get("sha256") != sha256_file(path) or \
                record.get("bytes") != path.stat().st_size:
            die("fresh-build evidence receipt mismatch: {}".format(name))
    for index, ((period, rows, mask), record) in enumerate(zip(SMOKE_CASES, smoke["cases"])):
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        stdout_path = root / "frozen" / "smoke" / (stem + ".stdout")
        stderr_path = root / "frozen" / "smoke" / (stem + ".stderr")
        task = {"period": period, "rows": rows, "mask": mask,
                "seed_index": 0, "schedule": "adversarial", "Ks": [257]}
        if not isinstance(record, dict) or record.get("case") != index or \
                record.get("period") != period or record.get("rows") != rows or \
                record.get("mask") != mask or record.get("K") != 257 or \
                record.get("stdout") != str(stdout_path.relative_to(root)) or \
                record.get("stderr") != str(stderr_path.relative_to(root)) or \
                record.get("stdout_sha256") != sha256_file(stdout_path) or \
                record.get("stderr_sha256") != sha256_file(stderr_path) or \
                stderr_path.stat().st_size:
            die("bounded holdout smoke receipt mismatch")
        parsed = parse_benchmark_csv(stdout_path.read_bytes(), task)
        if record.get("outcome") != classify(parsed[0]):
            die("bounded holdout smoke outcome receipt mismatch")
    summary = load_json(root / "frozen" / SOURCE_SUMMARY_NAME)
    if sha256_file(root / "frozen" / SOURCE_SUMMARY_NAME) != \
            SOURCE_SUMMARY_SHA256:
        die("frozen all-K validated-summary hash mismatch")
    plan = plan_from_summary(summary, str(root / "frozen" / BINARY_NAME))
    expected_files = {
        root / "hard_keys.jsonl": json_lines(plan["hard"]),
        root / "arms.jsonl": json_lines(plan["arms"]),
        root / "tasks_manifest.jsonl": json_lines(plan["tasks"]),
    }
    for path, expected in expected_files.items():
        if path.read_bytes() != expected:
            die("regenerated ledger mismatch: {}".format(path.name))
    tasks = plan["tasks"]
    total_cells = sum(len(task["Ks"]) for task in tasks)
    if design.get("task_count") != len(tasks) or \
            design.get("hard_task_count") != sum(
                task["stage"] == "hard" for task in tasks) or \
            design.get("arms") != plan["arms"] or \
            total_cells % len(plan["arms"]) != 0 or \
            design.get("cells_per_arm") != total_cells // len(plan["arms"]) or \
            design.get("total_cells") != total_cells:
        die("design Cartesian receipt mismatch")
    expected_outputs = {task["output_name"] for task in tasks}
    for directory, suffix in (("raw", ""), ("stderr", ".stderr"), ("exit", ".exit")):
        entries = list((root / directory).iterdir())
        if any(path.is_symlink() or not path.is_file() for path in entries):
            die("invalid {} runtime artifact".format(directory))
        actual = {path.name for path in entries}
        expected = {name + suffix for name in expected_outputs}
        partial = {name for name in actual if ".part." in name}
        if not (actual - partial) <= expected:
            die("unexpected {} output(s)".format(directory))
        if partial and not runtime_tolerant:
            die("partial {} output remains".format(directory))
        if any(not any(name.startswith(prefix + ".part.") for prefix in expected)
               for name in partial):
            die("unexpected partial {} output(s)".format(directory))
    return design, tasks


def job_receipt_path(root: Path, job: int) -> Path:
    return root / "job_receipts" / "job{:05d}.json".format(job)


def completed_jobs(root: Path, tasks: Sequence[dict]) -> Set[int]:
    completed: Set[int] = set()
    for task in tasks:
        name = task["output_name"]
        paths = (
            root / "raw" / name,
            root / "stderr" / (name + ".stderr"),
            root / "exit" / (name + ".exit"),
            job_receipt_path(root, int(task["job"])),
        )
        present = tuple(path.exists() for path in paths)
        if any(present) and not all(present):
            die("partial output triplet for job {}".format(task["job"]))
        if all(present):
            if any(path.is_symlink() or not path.is_file() for path in paths):
                die("completed job contains an indirect artifact")
            receipt = load_json(paths[3])
            if not paths[0].stat().st_size or paths[1].stat().st_size or \
                    paths[2].read_text(encoding="ascii") != "0\n" or \
                    receipt.get("schema") != SCHEMA + ".job" or \
                    not isinstance(receipt.get("job"), int) or \
                    isinstance(receipt.get("job"), bool) or \
                    receipt.get("job") != task["job"] or \
                    receipt.get("stage") != task["stage"] or \
                    receipt.get("output_name") != name or \
                    not isinstance(receipt.get("returncode"), int) or \
                    isinstance(receipt.get("returncode"), bool) or \
                    receipt.get("returncode") != 0 or \
                    receipt.get("stdout_sha256") != sha256_file(paths[0]) or \
                    receipt.get("stderr_sha256") != sha256_file(paths[1]) or \
                    receipt.get("exit_sha256") != sha256_file(paths[2]) or \
                    not isinstance(receipt.get("segment"), int) or \
                    isinstance(receipt.get("segment"), bool) or \
                    not isinstance(receipt.get("started_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("started_monotonic_s"), bool) or \
                    not isinstance(receipt.get("ended_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("ended_monotonic_s"), bool) or \
                    receipt["ended_monotonic_s"] < receipt["started_monotonic_s"]:
                die("unclean completed job {}".format(task["job"]))
            completed.add(int(task["job"]))
    return completed


def process_state(pid: int) -> Optional[str]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(CAT_PATH)
        result = subprocess.run(
            [str(SUDO_PATH), "-n", str(CAT_PATH), str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    try:
        return raw.rsplit(b") ", 1)[1].split(None, 1)[0].decode("ascii")
    except (IndexError, UnicodeError):
        return None


def process_start_ticks(pid: int) -> Optional[int]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(CAT_PATH)
        result = subprocess.run(
            [str(SUDO_PATH), "-n", str(CAT_PATH), str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    try:
        # Field 22 is the process start time; after removing PID and comm it is
        # zero-based field 19 in the remainder beginning with the state.
        return int(raw.rsplit(b") ", 1)[1].split()[19])
    except (IndexError, ValueError):
        return None


def process_parent_pid(pid: int) -> Optional[int]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(CAT_PATH)
        result = subprocess.run(
            [str(SUDO_PATH), "-n", str(CAT_PATH), str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    try:
        # After PID and comm, field 3 (state) is index 0 and PPID is index 1.
        return int(raw.rsplit(b") ", 1)[1].split()[1])
    except (IndexError, ValueError):
        return None


def process_alive(pid: int) -> bool:
    if not isinstance(pid, int) or isinstance(pid, bool):
        die("process PID is malformed")
    if pid <= 1:
        return False
    if pid > MAX_PROCESS_PID:
        die("process PID is out of range")
    try:
        os.kill(pid, 0)
        return process_state(pid) != "Z"
    except ProcessLookupError:
        return False
    except PermissionError:
        # Root-owned samplers return EPERM to a non-root kill(2).  A plain
        # exception must never be interpreted as process death.
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(KILL_PATH)
        try:
            result = subprocess.run(
                [str(SUDO_PATH), "-n", str(KILL_PATH), "-0", str(pid)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            die("cannot inspect permission-protected PID {}: {}".format(
                pid, exc))
        if result.returncode == 0:
            return process_state(pid) != "Z"
        if (Path("/proc") / str(pid)).exists():
            die("cannot determine liveness of permission-protected PID {}".format(pid))
        return False


def process_tokens(pid: int) -> Optional[List[str]]:
    path = Path("/proc") / str(pid) / "cmdline"
    try:
        raw = path.read_bytes()
    except PermissionError:
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(CAT_PATH)
        result = subprocess.run(
            [str(SUDO_PATH), "-n", str(CAT_PATH), str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    return [token.decode("utf-8", "replace") for token in raw.split(b"\0") if token]


def active_owner_live_identities(
    marker: Mapping[str, object],
    *,
    alive_probe: Optional[Any] = None,
    start_ticks_probe: Optional[Any] = None,
) -> List[dict]:
    """Validate an executing marker and return identities still provably live."""
    if marker.get("phase") != "executing":
        die("live-owner inspection requires an executing marker")
    protected = validate_owner_protected(
        marker.get("protected"), require_nonempty=True)
    alive = process_alive if alive_probe is None else alive_probe
    start_ticks = (
        process_start_ticks if start_ticks_probe is None else start_ticks_probe)
    identities: List[dict] = []
    for entry in protected:
        identity = (int(entry["pid"]), int(entry["start_ticks"]))
        if not alive(identity[0]):
            continue
        observed_start = start_ticks(identity[0])
        if observed_start is None:
            if not alive(identity[0]):
                continue
            die("cannot inspect a live protected process identity")
        if int(observed_start) == identity[1]:
            identities.append({
                "role": str(entry["role"]),
                "pid": identity[0],
                "start_ticks": identity[1],
            })
    return identities


def assert_current_controller_owner(root: Path) -> None:
    marker = load_active_owner_marker()
    if marker is None or marker.get("campaign_root") != str(root):
        die("campaign lost its active environment ownership")
    live = active_owner_live_identities(marker)
    current_start = process_start_ticks(os.getpid())
    if current_start is None or not any(
            entry["role"] == "controller" and
            entry["pid"] == os.getpid() and
            entry["start_ticks"] == current_start
            for entry in live):
        die("campaign controller is not the protected environment owner")


def parse_fuser_i2c_result(
    device: Path,
    returncode: int,
    stdout: bytes,
    stderr: bytes,
) -> Tuple[int, ...]:
    """Strictly parse one bounded privileged fuser invocation."""
    if returncode == 1:
        if stdout or stderr:
            die("I2C no-reader result contains diagnostic output")
        return ()
    if returncode != 0:
        die("cannot inspect I2C reader ownership")
    if re.fullmatch(rb"(?: +[1-9][0-9]*)+", stdout) is None:
        die("I2C reader output is noncanonical or empty")
    label = re.escape(os.fsencode(str(device))) + rb":[ ]*\n"
    if re.fullmatch(label, stderr) is None:
        die("I2C reader label output is noncanonical")
    return tuple(sorted(set(int(value) for value in stdout.split())))


def validate_trusted_root_tool(path: Path) -> None:
    try:
        identity = path.lstat()
    except OSError as exc:
        die("trusted root tool is unavailable: {}".format(exc))
    if not stat.S_ISREG(identity.st_mode) or identity.st_uid != 0 or \
            identity.st_mode & 0o022:
        die("trusted root tool is not root-owned and immutable: {}".format(path))


def live_i2c_holders(excluded_pids: Iterable[int] = ()) -> Dict[int, List[str]]:
    """Return processes observed holding a DIMM I2C bus, including renamed readers.

    The sampler opens its I2C buses once and retains those descriptors for its
    lifetime.  Inspecting descriptors is therefore a semantic guard rather
    than a filename convention.  The campaign already requires passwordless
    sudo to launch the sampler; use the same authority to inspect root-owned
    sampler descriptors.  A PID observed by fuser remains evidence even if
    that reader exits before command inspection.
    """
    excluded = {os.getpid(), *excluded_pids}
    validate_trusted_root_tool(SUDO_PATH)
    validate_trusted_root_tool(FUSER_PATH)
    inventory: Dict[int, Set[str]] = defaultdict(set)
    for device in DIMM_I2C_DEVICES:
        if device.is_symlink() or not device.is_char_device():
            die("required I2C device is missing or unsafe: {}".format(device))
        try:
            result = subprocess.run(
                [str(SUDO_PATH), "-n", str(FUSER_PATH), str(device)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            die("cannot inspect live I2C readers: {}".format(exc))
        for pid in parse_fuser_i2c_result(
                device, result.returncode, result.stdout, result.stderr):
            inventory[pid].add(str(device))
    return {
        pid: sorted(devices) for pid, devices in sorted(inventory.items())
        if pid not in excluded
    }


def validate_sole_sampler_inventory(
    inventory: Mapping[int, Sequence[str]],
    sampled_pid: int,
    wrapper_pid: int,
) -> Dict[int, List[str]]:
    """Require the launched sampler on both buses and return every competitor."""
    expected = sorted(str(device) for device in DIMM_I2C_DEVICES)
    observed = list(inventory.get(sampled_pid, ()))
    if observed != expected:
        die("launched sampler does not hold the complete DIMM I2C inventory")
    # Only the process named by the sampler's own PID receipt opens the I2C
    # descriptors.  The sudo wrapper is not an I2C owner and must not be
    # excluded: its PID could have been recycled by a foreign reader.
    owned = {int(sampled_pid)}
    return {
        int(pid): list(devices)
        for pid, devices in sorted(inventory.items())
        if int(pid) not in owned
    }


def sole_sampler_competitors(
    sampled_pid: int,
    wrapper_pid: int,
) -> Dict[int, List[str]]:
    return validate_sole_sampler_inventory(
        live_i2c_holders(), sampled_pid, wrapper_pid)


def is_wirehair_sampler_command(tokens: Sequence[str]) -> bool:
    """Recognize the known sampler family before it has opened its buses."""
    if tokens and Path(tokens[0]).name in ("sh", "bash", "dash"):
        try:
            command_index = list(tokens).index("-c") + 1
            nested = shlex.split(tokens[command_index])
        except (ValueError, IndexError):
            nested = []
        if nested and is_wirehair_sampler_command(nested):
            return True

    def sampler_name(token: str) -> bool:
        name = Path(token).name
        return name == SAMPLER_NAME or (
            name.startswith("wirehair_expo_thermal_sampler") and
            name.endswith(".py"))

    has_sampler_script = False
    for index, token in enumerate(tokens):
        if not sampler_name(token):
            continue
        if index == 0:
            has_sampler_script = True
            break
        python_indices = [
            prior for prior in range(index)
            if re.fullmatch(r"python(?:[0-9]+(?:\.[0-9]+)*)?",
                            Path(tokens[prior]).name)
        ]
        if python_indices:
            interpreter = python_indices[-1]
            if not any(
                    Path(tokens[prior]).name.endswith(".py")
                    for prior in range(interpreter + 1, index)):
                has_sampler_script = True
                break
    # Both options are required by every sampler in this family.  Requiring
    # them avoids treating a short-lived editor, hasher, or reviewer that
    # merely names the sampler source file as an active hardware reader.
    has_csv = any(
        token == "--csv" or token.startswith("--csv=") for token in tokens)
    has_pid_file = any(
        token == "--pid-file" or token.startswith("--pid-file=")
        for token in tokens)
    return has_sampler_script and has_csv and has_pid_file


def other_samplers(excluded_pids: Iterable[int] = ()) -> List[dict]:
    excluded = {os.getpid(), *excluded_pids}
    found: Dict[int, dict] = {}
    for pid, devices in live_i2c_holders(excluded).items():
        tokens = process_tokens(pid)
        found[pid] = {
            "pid": pid,
            "command": (
                " ".join(tokens) if tokens is not None else "<unreadable>"),
            "i2c_devices": devices,
        }
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) in excluded:
            continue
        pid = int(proc.name)
        tokens = process_tokens(pid)
        if tokens is not None and is_wirehair_sampler_command(tokens) and \
                process_alive(pid) and pid not in found:
            found[pid] = {
                "pid": pid, "command": " ".join(tokens), "i2c_devices": [],
            }
    return [found[pid] for pid in sorted(found)]


def atomic_result(path: Path, data: bytes) -> None:
    if path.exists() or path.is_symlink():
        die("refusing to replace runtime artifact {}".format(path))
    temporary = Path(str(path) + ".part.{}.{}".format(
        os.getpid(), threading.get_ident()))
    with temporary.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def segment_path(root: Path, segment: int, kind: str) -> Path:
    return root / "segments" / "segment{:03d}.{}.json".format(segment, kind)


def segment_indices(root: Path) -> List[int]:
    verify_runtime_directories(root)
    indices: Set[int] = set()
    patterns = (
        (root / "segments",
         re.compile(r"segment([0-9]{3})\.(?:intent|ready|final)\.json$"),
         "file"),
        (root / "thermal",
         re.compile(r"segment([0-9]{3})\.(?:csv|pid|stderr)$"),
         "file"),
        (root / "attempts", re.compile(r"segment([0-9]{3})$"), "directory"),
    )
    for directory, pattern, expected_kind in patterns:
        for path in directory.iterdir():
            match = pattern.fullmatch(path.name)
            if match:
                valid = (
                    not path.is_symlink() and
                    (path.is_file() if expected_kind == "file" else path.is_dir())
                )
                if not valid:
                    die("campaign segment evidence has the wrong file type: {}"
                        .format(path))
                indices.add(int(match.group(1)))
    ordered = sorted(indices)
    if ordered and ordered != list(range(ordered[-1] + 1)):
        die("campaign segment numbering has a gap")
    return ordered


def segment_jobs(intent: Mapping[str, object], tasks: Sequence[dict]) -> List[dict]:
    raw = intent.get("jobs")
    if not isinstance(raw, list) or not raw or any(
            not isinstance(job, int) or isinstance(job, bool) for job in raw):
        die("segment intent has an invalid job ledger")
    if raw != sorted(set(raw)):
        die("segment intent job ledger is not canonical")
    task_by_job = {int(task["job"]): task for task in tasks}
    if any(job not in task_by_job for job in raw):
        die("segment intent references an unknown job")
    selected = [task_by_job[job] for job in raw]
    stage = intent.get("stage")
    if stage not in ("hard", "control") or any(task["stage"] != stage for task in selected):
        die("segment intent crosses campaign stages")
    if intent.get("jobs_sha256") != sha256_bytes(json_lines(raw)):
        die("segment intent job hash mismatch")
    return selected


def terminate_verified_process(
    pid: int,
    expected_start_ticks: int,
    expected: Sequence[str],
    *,
    alive_probe: Optional[Any] = None,
    start_ticks_probe: Optional[Any] = None,
    tokens_probe: Optional[Any] = None,
    signal_sender: Optional[Any] = None,
) -> str:
    """Stop one owned identity without ever signaling a reused PID."""
    if not isinstance(expected_start_ticks, int) or \
            isinstance(expected_start_ticks, bool) or expected_start_ticks < 0:
        die("termination requires a valid expected process start time")
    alive = process_alive if alive_probe is None else alive_probe
    start_ticks = (
        process_start_ticks if start_ticks_probe is None else start_ticks_probe)
    read_tokens = process_tokens if tokens_probe is None else tokens_probe

    def identity_state() -> str:
        if not alive(pid):
            return "exited"
        observed_start = start_ticks(pid)
        if observed_start is None:
            if not alive(pid):
                return "exited"
            die("cannot re-prove live process start identity")
        if int(observed_start) != expected_start_ticks:
            return "reused"
        observed_tokens = read_tokens(pid)
        if observed_tokens is None:
            if not alive(pid):
                return "exited"
            die("cannot re-prove live process command identity")
        if observed_tokens != list(expected):
            if not alive(pid):
                return "exited"
            return "changed"
        return "match"

    def send(signal_name: str) -> bool:
        if signal_sender is not None:
            result = signal_sender(signal_name)
            return True if result is None else bool(result)
        validate_trusted_root_tool(SUDO_PATH)
        validate_trusted_root_tool(PYTHON_PATH)
        expected_cmdline = (
            b"\0".join(os.fsencode(token) for token in expected) + b"\0")
        try:
            result = subprocess.run(
                [
                    str(SUDO_PATH), "-n", str(PYTHON_PATH), "-I", "-c",
                    PRIVILEGED_PIDFD_SIGNAL_CODE, str(pid),
                    str(expected_start_ticks), expected_cmdline.hex(),
                    signal_name,
                ],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            die("privileged pidfd signal helper failed: {}".format(exc))
        if result.returncode in (3, 4) and not result.stdout and not result.stderr:
            return False
        if result.returncode or result.stdout or result.stderr:
            die("privileged pidfd signal helper failed closed")
        return True

    def classify_nonmatch(state: str, exited_action: str) -> str:
        if state == "exited":
            return exited_action
        if state == "reused":
            return "identity_reused"
        if state == "changed":
            return "identity_changed"
        die("process identity unexpectedly remained matched")

    initial = identity_state()
    if initial == "exited":
        return "already_exited"
    if initial != "match":
        return "identity_" + initial
    # Re-prove in the caller, then let the privileged helper bind a pidfd and
    # repeat the start/cmdline proof before signaling that process object.
    state = identity_state()
    if state != "match":
        return classify_nonmatch(state, "already_exited")
    if not send("TERM"):
        return classify_nonmatch(identity_state(), "terminated")
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline:
        state = identity_state()
        if state != "match":
            return classify_nonmatch(state, "terminated")
        time.sleep(0.05)
    state = identity_state()
    if state != "match":
        return classify_nonmatch(state, "terminated")
    if not send("KILL"):
        return classify_nonmatch(identity_state(), "killed")
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        state = identity_state()
        if state != "match":
            return classify_nonmatch(state, "killed")
        time.sleep(0.05)
    state = identity_state()
    if state == "match":
        die("verified owned process {} did not stop".format(pid))
    return classify_nonmatch(state, "killed")


def assert_process_identity(
    pid: int,
    expected_start_ticks: int,
    expected_tokens: Sequence[str],
    label: str,
) -> None:
    if not process_alive(pid):
        die("{} exited".format(label))
    observed_start = process_start_ticks(pid)
    observed_tokens = process_tokens(pid)
    if observed_start != expected_start_ticks or \
            observed_tokens != list(expected_tokens):
        die("{} identity changed".format(label))


def terminate_direct_child_by_pidfd(
    process: subprocess.Popen,
    pidfd: int,
) -> str:
    """Boundedly stop a directly spawned user process without PID reuse risk."""
    if process.poll() is not None:
        return "already_exited"
    try:
        signal.pidfd_send_signal(pidfd, signal.SIGTERM, None, 0)
    except ProcessLookupError:
        process.wait(timeout=5.0)
        return "already_exited"
    try:
        process.wait(timeout=10.0)
        return "terminated"
    except subprocess.TimeoutExpired:
        pass
    try:
        signal.pidfd_send_signal(pidfd, signal.SIGKILL, None, 0)
    except ProcessLookupError:
        process.wait(timeout=5.0)
        return "terminated"
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        die("direct child did not stop after pidfd SIGKILL")
    return "killed"


def terminate_direct_unreaped_child(process: subprocess.Popen) -> str:
    """Stop a Popen child before wait/reap allows its PID to be reused."""
    if process.poll() is not None:
        return "already_exited"
    process.terminate()
    try:
        process.wait(timeout=10.0)
        return "terminated"
    except subprocess.TimeoutExpired:
        process.kill()
        try:
            process.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            die("direct unreaped child did not stop after SIGKILL")
        return "killed"


def terminate_privileged_direct_unreaped_child(
    process: subprocess.Popen,
) -> str:
    """Stop our unreaped privileged Popen child with fixed, bounded tools."""
    if process.poll() is not None:
        return "already_exited"
    validate_trusted_root_tool(SUDO_PATH)
    validate_trusted_root_tool(KILL_PATH)

    def send(signal_name: str) -> None:
        try:
            result = subprocess.run(
                [str(SUDO_PATH), "-n", str(KILL_PATH),
                 "-{}".format(signal_name), str(process.pid)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                timeout=I2C_INSPECTION_TIMEOUT_SECONDS,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            die("cannot stop privileged direct child: {}".format(exc))
        if result.returncode and process.poll() is None:
            die("cannot signal privileged direct child")

    send("TERM")
    try:
        process.wait(timeout=10.0)
        return "terminated"
    except subprocess.TimeoutExpired:
        send("KILL")
        try:
            process.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            die("privileged direct child did not stop after SIGKILL")
        return "killed"


def require_stopped_action(action: str, label: str) -> None:
    if action not in STOPPED_PROCESS_ACTIONS:
        die("{} cleanup lost process identity: {}".format(label, action))


def attempt_all_cleanups(
    entries: Sequence[Tuple[int, Any]],
    cleanup: Any,
) -> Tuple[List[Tuple[int, Any, Optional[str]]], List[str]]:
    """Run every cleanup leg even when an earlier identity fails."""
    results: List[Tuple[int, Any, Optional[str]]] = []
    errors: List[str] = []
    for job, process in entries:
        try:
            results.append((job, process, cleanup(job, process)))
        except BaseException as exc:
            errors.append("job{}:{!r}".format(job, exc))
    return results, errors


def campaign_sampled_command(
    sampler_path: Path,
    csv_path: Path,
    pid_file: Path,
) -> List[str]:
    return [
        str(PYTHON_PATH), "-I", str(sampler_path),
        "--csv", str(csv_path), "--pid-file", str(pid_file),
        "--interval", str(THERMAL_INTERVAL_SECONDS),
        "--dimm-attempts", "5", "--dimm-retry-delay", "0.01",
    ]


def campaign_sampler_supervisor_commands(
    sampler_path: Path,
    csv_path: Path,
    pid_file: Path,
    controller_identity: Tuple[int, int],
) -> Tuple[List[str], List[str]]:
    controller_pid, controller_start = controller_identity
    sampled = campaign_sampled_command(sampler_path, csv_path, pid_file)
    supervisor = [
        str(PYTHON_PATH), "-I", "-c", PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
        str(controller_pid), str(controller_start),
        PDEATH_EXEC_CODE.encode("utf-8").hex(), *sampled,
    ]
    return supervisor, [str(SUDO_PATH), "-n", *supervisor]


def campaign_sampler_command_matches(
    tokens: Sequence[str],
    sampler_path: Path,
    csv_path: Path,
    pid_file: Path,
    controller_identity: Tuple[int, int],
) -> bool:
    """Match only the controller-bound sudo/supervisor command.

    The sampled child loses the controller identity from its argv after exec.
    It is therefore signalable only through a durable PID/start receipt, not
    merely because an unrelated process happens to use the same sampler paths.
    """
    supervisor, sudo_wrapper = campaign_sampler_supervisor_commands(
        sampler_path, csv_path, pid_file, controller_identity)
    return list(tokens) in (supervisor, sudo_wrapper)


def launch_thermal_sampler(
    command: Sequence[str],
    sampler_error: Any,
) -> subprocess.Popen:
    """Launch sudo outside any caller terminal so its ancestry is invariant."""
    return subprocess.Popen(
        list(command), stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL, stderr=sampler_error,
        start_new_session=True,
    )


def benchmark_command_matches(
    tokens: Sequence[str],
    expected: Sequence[str],
    controller_identity: Optional[Tuple[int, int]] = None,
) -> bool:
    if list(tokens) == list(expected):
        # Exec removes the controller identity from argv.  Exact bare commands
        # are useful for observation, but are not signal authority during
        # crash reconciliation.
        return controller_identity is None
    prefix = [str(PYTHON_PATH), "-I", "-c", PDEATH_EXEC_CODE]
    if not (len(tokens) == len(prefix) + 2 + len(expected) and
            list(tokens[:len(prefix)]) == prefix and
            str(tokens[len(prefix)]).isdigit() and
            str(tokens[len(prefix) + 1]).isdigit() and
            list(tokens[len(prefix) + 2:]) == list(expected)):
        return False
    return controller_identity is None or (
        int(tokens[len(prefix)]) == controller_identity[0] and
        int(tokens[len(prefix) + 1]) == controller_identity[1])


def discover_owned_processes(
    matchers: Sequence[Tuple[str, Any]],
) -> List[Tuple[str, int, int, List[str]]]:
    """Capture current start/cmdline identities for uniquely matched commands."""
    found: List[Tuple[str, int, int, List[str]]] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        pid = int(proc.name)
        if pid == os.getpid():
            continue
        tokens = process_tokens(pid)
        if tokens is None:
            continue
        labels = [label for label, matcher in matchers if matcher(tokens)]
        if not labels:
            continue
        if len(labels) != 1:
            die("owned process matches more than one sealed command")
        if not process_alive(pid):
            continue
        start_ticks = process_start_ticks(pid)
        rebound_tokens = process_tokens(pid)
        if start_ticks is None or rebound_tokens is None:
            if not process_alive(pid):
                continue
            die("cannot bind discovered owned process identity")
        if rebound_tokens != tokens:
            # A sudo/Python wrapper may legitimately exec between reads.  The
            # bounded caller rescans and binds its post-exec identity.
            continue
        found.append((labels[0], pid, start_ticks, rebound_tokens))
    return sorted(found, key=lambda item: (item[1], item[0]))


def terminate_discovered_owned_processes(
    matchers: Sequence[Tuple[str, Any]],
) -> List[dict]:
    """Reconcile spawn-to-receipt crashes by exact unique command identity."""
    actions: List[dict] = []
    empty_observations = 0
    for _ in range(12):
        identities = discover_owned_processes(matchers)
        if not identities:
            empty_observations += 1
            if empty_observations >= 2:
                return actions
            time.sleep(0.05)
            continue
        empty_observations = 0
        cleanup_errors: List[str] = []
        for label, pid, start_ticks, tokens in identities:
            try:
                action = terminate_verified_process(
                    pid, start_ticks, tokens)
                if action not in DISCOVERY_PROVED_GONE_ACTIONS:
                    die("discovered owned process cleanup lost identity: {}"
                        .format(action))
                actions.append({
                    "kind": label, "pid": pid, "action": action,
                })
            except BaseException as exc:
                # Continue through the whole captured set so one fault cannot
                # strand later children.  Any error still leaves the segment
                # incomplete and prevents receipt publication.
                cleanup_errors.append(
                    "{}:pid{}:{!r}".format(label, pid, exc))
        if cleanup_errors:
            die("discovered owned process cleanup failed closed: {}".format(
                "; ".join(cleanup_errors)))
    die("exact owned command remains live after bounded reconciliation")


def terminate_owned_segment_processes(
    root: Path,
    segment: int,
    tasks: Sequence[dict],
    controller_identity: Optional[Tuple[int, int]] = None,
) -> List[dict]:
    verify_runtime_directories(root)
    if controller_identity is None:
        die("segment reconciliation requires the original controller identity")
    actions: List[dict] = []
    cleanup_errors: List[str] = []
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    sampler_path = root / "frozen" / SAMPLER_NAME
    sampler_matchers: List[Tuple[str, Any]] = [(
        "thermal",
        lambda tokens: campaign_sampler_command_matches(
            tokens, sampler_path, thermal_csv, thermal_pid,
            controller_identity),
    )]
    benchmark_matchers: List[Tuple[str, Any]] = [
        (
            "benchmark_job_{}".format(int(task["job"])),
            lambda tokens, expected=list(task["argv"]):
                benchmark_command_matches(
                    tokens, expected, controller_identity),
        )
        for task in tasks
    ]
    task_by_job = {int(task["job"]): task for task in tasks}
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)

    # Crash cleanup is best-effort exhaustive but fail-closed: try both
    # independently identifiable wrapper categories before consulting any
    # potentially corrupt receipt, then try every remaining safe leg.  One
    # malformed sampler artifact must never strand benchmark children (or
    # vice versa).
    thermal_wrapper_cleanup_proved = False
    for label, matchers in (
            ("thermal_wrapper_discovery", sampler_matchers),
            ("benchmark_wrapper_discovery", benchmark_matchers)):
        try:
            actions.extend(terminate_discovered_owned_processes(matchers))
            if label == "thermal_wrapper_discovery":
                thermal_wrapper_cleanup_proved = True
        except BaseException as exc:
            cleanup_errors.append("{}:{!r}".format(label, exc))

    incomplete_thermal_publication = False
    unbound_thermal_publication = False
    thermal_publication: Optional[SamplerPidPublication] = None
    try:
        pid, thermal_publication = read_sampler_pid_publication(thermal_pid)
        if thermal_publication is not None and pid is None:
            incomplete_thermal_publication = True
        if pid is not None:
            expected_thermal_start: Optional[int] = None
            ready_path = segment_path(root, segment, "ready")
            if ready_path.is_file() and not ready_path.is_symlink():
                ready = load_json(ready_path)
                if ready.get("sampler_pid") == pid and \
                        isinstance(ready.get("sampler_start_ticks"), int) and \
                        not isinstance(
                            ready.get("sampler_start_ticks"), bool):
                    expected_thermal_start = int(
                        ready["sampler_start_ticks"])
            receipt_stale = False
            owned_stopped = False
            initially_alive = process_alive(pid)
            if initially_alive:
                if expected_thermal_start is None:
                    # A complete PID can be published before the ready receipt.
                    # It is not signal authority: defer it to the same exact
                    # wrapper/direct-command/I2C absence proof used for an
                    # incomplete publication.
                    unbound_thermal_publication = True
                else:
                    expected_tokens = campaign_sampled_command(
                        sampler_path, thermal_csv, thermal_pid)
                    action = terminate_verified_process(
                        pid, expected_thermal_start, expected_tokens)
                    if action == "identity_reused":
                        # The receipt's original process object is gone.  Never
                        # signal its numeric replacement.
                        actions.append({
                            "kind": "thermal", "pid": pid,
                            "action": "stale_reused_pid",
                        })
                        receipt_stale = True
                    else:
                        actions.append({
                            "kind": "thermal", "pid": pid, "action": action,
                        })
                        require_stopped_action(
                            action, "receipted thermal reconciliation")
                        owned_stopped = True
            if not unbound_thermal_publication:
                still_alive = process_alive(pid)
                if still_alive and process_start_ticks(pid) is None:
                    if process_alive(pid):
                        die("cannot re-prove a live thermal PID start identity")
                    still_alive = False
                if still_alive and not initially_alive:
                    actions.append({
                        "kind": "thermal", "pid": pid,
                        "action": "stale_reused_pid",
                    })
                    receipt_stale = True
                if not still_alive or receipt_stale or owned_stopped:
                    try:
                        thermal_pid.unlink()
                    except FileNotFoundError:
                        pass
    except BaseException as exc:
        cleanup_errors.append("thermal_receipt:{!r}".format(exc))
    finally:
        try:
            close_sampler_pid_publication(thermal_publication)
        except BaseException as exc:
            cleanup_errors.append(
                "thermal_receipt_close:{!r}".format(exc))

    direct_sampler_blocked = True
    try:
        unreceipted_samplers = discover_owned_processes([(
            "unreceipted_thermal",
            lambda tokens: list(tokens) == campaign_sampled_command(
                sampler_path, thermal_csv, thermal_pid),
        )])
        if unreceipted_samplers:
            die("unreceipted direct sampler command remains live; refusing "
                "command-only signal authority")
        direct_sampler_blocked = False
    except BaseException as exc:
        cleanup_errors.append("thermal_direct_blocker:{!r}".format(exc))
    sampler_inventory_empty = False
    try:
        remaining_samplers = other_samplers()
        if remaining_samplers:
            die("sampler/I2C inventory remains live during reconciliation: {}"
                .format(json.dumps(remaining_samplers, sort_keys=True)))
        sampler_inventory_empty = True
    except BaseException as exc:
        cleanup_errors.append("thermal_inventory:{!r}".format(exc))
    if (incomplete_thermal_publication or
            unbound_thermal_publication) and \
            thermal_wrapper_cleanup_proved and \
            not direct_sampler_blocked and sampler_inventory_empty:
        try:
            thermal_pid.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            cleanup_errors.append(
                "thermal_partial_receipt_cleanup:{!r}".format(exc))

    def cleanup_benchmark_receipt(pid_path: Path) -> None:
        if pid_path.is_symlink() or not pid_path.is_file():
            die("benchmark PID receipt is indirect")
        match = re.fullmatch(r"job([0-9]{5})\.pid", pid_path.name)
        if not match:
            die("malformed benchmark PID receipt")
        job = int(match.group(1))
        if job not in task_by_job:
            die("benchmark PID receipt names an unknown job")
        identity = load_json(pid_path)
        pid = identity.get("pid")
        start_ticks = identity.get("start_ticks")
        if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 1 or \
                pid > MAX_PROCESS_PID or \
                not isinstance(start_ticks, int) or \
                isinstance(start_ticks, bool) or start_ticks < 0:
            die("malformed benchmark PID identity")
        if not process_alive(pid):
            return
        expected_tokens = list(task_by_job[job]["argv"])
        action = terminate_verified_process(
            pid, start_ticks, expected_tokens)
        if action == "identity_reused":
            actions.append({
                "kind": "benchmark", "job": job, "pid": pid,
                "action": "stale_reused_pid",
            })
            return
        actions.append({
            "kind": "benchmark", "job": job, "pid": pid,
            "action": action,
        })
        require_stopped_action(
            action, "benchmark reconciliation job {}".format(job))

    receipt_paths: List[Path] = []
    try:
        if attempt_dir.exists() or attempt_dir.is_symlink():
            if attempt_dir.is_symlink() or not attempt_dir.is_dir():
                die("benchmark attempt directory is indirect")
            receipt_paths = sorted(attempt_dir.glob("*.pid"))
    except BaseException as exc:
        cleanup_errors.append("benchmark_receipts:{!r}".format(exc))
    for pid_path in receipt_paths:
        try:
            cleanup_benchmark_receipt(pid_path)
        except BaseException as exc:
            cleanup_errors.append(
                "benchmark_receipt_{}:{!r}".format(
                    pid_path.name, exc))

    try:
        unreceipted_benchmarks = discover_owned_processes([
            (
                "unreceipted_benchmark_job_{}".format(int(task["job"])),
                lambda tokens, expected=list(task["argv"]):
                    list(tokens) == expected,
            )
            for task in tasks
        ])
        if unreceipted_benchmarks:
            die("unreceipted direct benchmark command remains live; refusing "
                "command-only signal authority")
    except BaseException as exc:
        cleanup_errors.append("benchmark_direct_blocker:{!r}".format(exc))

    if cleanup_errors:
        die("segment process cleanup failed closed: {}".format(
            "; ".join(cleanup_errors)))
    return actions


def rollback_segment_outputs(root: Path, selected: Sequence[dict]) -> List[int]:
    verify_runtime_directories(root)
    rolled_back: List[int] = []
    for task in selected:
        job = int(task["job"])
        name = task["output_name"]
        paths = (
            root / "raw" / name,
            root / "stderr" / (name + ".stderr"),
            root / "exit" / (name + ".exit"),
            job_receipt_path(root, job),
        )
        if any(path.exists() or path.is_symlink() for path in paths):
            rolled_back.append(job)
        for path in paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        for directory, prefix in (
                (root / "raw", name), (root / "stderr", name + ".stderr"),
                (root / "exit", name + ".exit")):
            for partial in directory.glob(prefix + ".part.*"):
                partial.unlink()
    return rolled_back


def seal_attempt_manifest(
    root: Path,
    segment: int,
    allowed_jobs: Set[int],
    *,
    require_complete: bool,
) -> Tuple[str, int]:
    directory = root / "attempts" / "segment{:03d}".format(segment)
    if not directory.exists():
        make_private_directory(directory)
    if directory.is_symlink() or not directory.is_dir():
        die("attempt evidence directory is indirect")
    for partial in directory.glob("*.part.*"):
        partial.unlink()
    files = sorted(directory.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in files):
        die("attempt evidence contains an indirect entry")
    pattern = re.compile(r"job([0-9]{5})\.(pid|stdout|stderr|exit)$")
    matches = [pattern.fullmatch(path.name) for path in files]
    if any(match is None for match in matches):
        die("attempt evidence contains an unexpected file")
    observed = {
        (int(match.group(1)), match.group(2))
        for match in matches if match is not None
    }
    if any(job not in allowed_jobs for job, _ in observed):
        die("attempt evidence references a job outside the segment")
    groups: Dict[int, Set[str]] = defaultdict(set)
    for job, suffix in observed:
        groups[job].add(suffix)
    valid_prefixes = (
        {"pid"},
        {"pid", "stdout"},
        {"pid", "stdout", "stderr"},
        {"pid", "stdout", "stderr", "exit"},
    )
    if any(suffixes not in valid_prefixes for suffixes in groups.values()):
        die("attempt evidence is not an atomic output prefix")
    if require_complete and observed != {
            (job, suffix)
            for job in allowed_jobs
            for suffix in ("pid", "stdout", "stderr", "exit")}:
        die("successful segment attempt evidence is incomplete")
    data = b"".join(
        "{}  {}\n".format(sha256_file(path), path.relative_to(root)).encode("ascii")
        for path in files
    )
    manifest = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
    write_once(manifest, data)
    return sha256_bytes(data), len(files)


def write_segment_final(
    root: Path,
    segment: int,
    intent: Mapping[str, object],
    state: str,
    *,
    failures: Sequence[dict],
    rolled_back: Sequence[int],
    process_actions: Sequence[dict],
    jobs_ended_monotonic_s: Optional[float],
    sampler_returncode: Optional[int],
    validated_thermal_csv_sha256: Optional[str] = None,
) -> Dict[str, object]:
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    thermal_csv_sha256 = bounded_thermal_sha256(thermal_csv)
    thermal_stderr_sha256 = bounded_thermal_stderr_sha256(thermal_stderr)
    if state == "success":
        if not isinstance(validated_thermal_csv_sha256, str) or \
                re.fullmatch(r"[0-9a-f]{64}",
                             validated_thermal_csv_sha256) is None or \
                thermal_csv_sha256 != validated_thermal_csv_sha256:
            die("successful segment thermal bytes changed after validation")
        if thermal_stderr_sha256 != sha256_bytes(b""):
            die("successful segment thermal stderr is not canonically empty")
    elif validated_thermal_csv_sha256 is not None:
        die("non-success segment cannot claim validated thermal bytes")
    raw_jobs = intent.get("jobs")
    if not isinstance(raw_jobs, list) or any(
            not isinstance(job, int) or isinstance(job, bool)
            for job in raw_jobs):
        die("segment final lacks a valid intent job ledger")
    attempts_sha256, attempt_file_count = seal_attempt_manifest(
        root, segment, set(raw_jobs), require_complete=state == "success")
    record = {
        "schema": SCHEMA + ".segment_final", "segment": segment,
        "state": state, "stage": intent.get("stage"),
        "intent_sha256": sha256_file(segment_path(root, segment, "intent")),
        "ready_sha256": (
            sha256_file(segment_path(root, segment, "ready"))
            if segment_path(root, segment, "ready").is_file() else None),
        "ended_utc": datetime.now(timezone.utc).isoformat(),
        "jobs_ended_monotonic_s": jobs_ended_monotonic_s,
        "jobs": intent.get("jobs"), "jobs_sha256": intent.get("jobs_sha256"),
        "failures": list(failures), "rolled_back_jobs": sorted(rolled_back),
        "published_jobs": list(intent.get("jobs", [])) if state == "success" else [],
        "process_actions": list(process_actions),
        "sampler_returncode": sampler_returncode,
        "thermal_csv": thermal_csv.name,
        "thermal_csv_sha256": thermal_csv_sha256,
        "thermal_stderr": thermal_stderr.name,
        "thermal_stderr_sha256": thermal_stderr_sha256,
        "attempt_manifest": "segment{:03d}.attempts.sha256".format(segment),
        "attempt_manifest_sha256": attempts_sha256,
        "attempt_file_count": attempt_file_count,
        "retry_policy": "entire stage segment; no successful survivor is retained on failure",
    }
    write_once(segment_path(root, segment, "final"), canonical_json(record))
    return record


def reconcile_incomplete_segments(root: Path, tasks: Sequence[dict]) -> List[dict]:
    reconciled: List[dict] = []
    current_boot_id = required_boot_id()
    for segment in segment_indices(root):
        intent_path = segment_path(root, segment, "intent")
        final_path = segment_path(root, segment, "final")
        if final_path.exists():
            if not intent_path.exists():
                die("finalized segment lacks its intent")
            continue
        if not intent_path.exists():
            die("campaign segment artifacts lack an immutable intent")
        intent = load_json(intent_path)
        controller_pid = intent.get("controller_pid")
        controller_start = intent.get("controller_start_ticks")
        segment_boot_id = intent.get("boot_id")
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                not isinstance(intent.get("segment"), int) or \
                isinstance(intent.get("segment"), bool) or \
                intent.get("segment") != segment or \
                not valid_segment_boot_id(segment_boot_id) or \
                not isinstance(controller_pid, int) or \
                isinstance(controller_pid, bool) or controller_pid <= 1 or \
                controller_pid > MAX_PROCESS_PID or \
                not isinstance(controller_start, int) or \
                isinstance(controller_start, bool) or controller_start < 0:
            die("incomplete segment intent is malformed")
        selected = segment_jobs(intent, tasks)
        if segment_boot_id == current_boot_id:
            actions = terminate_owned_segment_processes(
                root, segment, selected,
                (controller_pid, controller_start))
        else:
            # Linux start ticks are boot-relative.  A boot change proves every
            # old process object is gone, but forbids interpreting any current
            # PID/cmdline collision as signal authority.
            current_samplers = other_samplers()
            if current_samplers:
                die("prior-boot reconciliation is blocked by current "
                    "sampler(s): {}".format(json.dumps(
                        current_samplers, sort_keys=True)))
            actions = [{
                "kind": "segment", "action": "boot_changed",
                "from_boot_id": segment_boot_id,
                "to_boot_id": current_boot_id,
            }]
            try:
                (root / "thermal" /
                 "segment{:03d}.pid".format(segment)).unlink()
            except FileNotFoundError:
                pass
            except OSError as exc:
                die("cannot remove stale prior-boot thermal PID receipt: {}"
                    .format(exc))
        rolled_back = rollback_segment_outputs(root, selected)
        final = write_segment_final(
            root, segment, intent, "interrupted", failures=[],
            rolled_back=rolled_back, process_actions=actions,
            jobs_ended_monotonic_s=None, sampler_returncode=None,
        )
        reconciled.append(final)
    return reconciled


def read_bounded_regular_bytes(
    path: Path,
    *,
    allow_missing: bool,
    max_bytes: int,
    label: str,
) -> Optional[Tuple[bytes, bool]]:
    """Read at most one bounded, path-bound regular-file snapshot."""
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
        getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        descriptor = os.open(str(path), flags)
    except FileNotFoundError:
        if allow_missing:
            return None
        die("{} is missing: {}".format(label, path))
    except OSError as exc:
        die("cannot read {} {}: {}".format(label, path, exc))
    try:
        identity = os.fstat(descriptor)
        if not stat.S_ISREG(identity.st_mode) or identity.st_nlink != 1:
            die("{} is not a direct regular file: {}".format(label, path))
        chunks: List[bytes] = []
        remaining = max_bytes + 1
        while remaining > 0:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        try:
            path_identity = path.lstat()
        except FileNotFoundError:
            die("{} disappeared during read: {}".format(label, path))
        expected_metadata = (
            identity.st_uid, identity.st_gid,
            stat.S_IMODE(identity.st_mode), identity.st_nlink)
        if len(raw) > max_bytes:
            die("{} exceeds the bounded size policy: {}".format(label, path))
        if not stat.S_ISREG(after.st_mode) or after.st_nlink != 1 or \
                not stat.S_ISREG(path_identity.st_mode) or \
                path_identity.st_nlink != 1 or \
                (identity.st_dev, identity.st_ino) != \
                (after.st_dev, after.st_ino) or \
                (after.st_dev, after.st_ino) != \
                (path_identity.st_dev, path_identity.st_ino) or \
                (
                    after.st_uid, after.st_gid,
                    stat.S_IMODE(after.st_mode), after.st_nlink) != \
                expected_metadata or \
                (
                    path_identity.st_uid, path_identity.st_gid,
                    stat.S_IMODE(path_identity.st_mode),
                    path_identity.st_nlink) != expected_metadata:
            die("{} identity changed during read: {}".format(label, path))
        stable = after.st_size == len(raw)
    except OSError as exc:
        die("cannot read {} {}: {}".format(label, path, exc))
    finally:
        try:
            os.close(descriptor)
        except OSError as exc:
            die("cannot close {} {}: {}".format(label, path, exc))
    return raw, stable


def read_thermal_snapshot_with_sha256(
    path: Path,
    *,
    allow_incomplete: bool,
) -> Optional[Tuple[List[Dict[str, str]], str]]:
    snapshot = read_bounded_regular_bytes(
        path, allow_missing=allow_incomplete,
        max_bytes=MAX_THERMAL_STREAM_BYTES, label="thermal stream")
    if snapshot is None:
        return None
    raw, stable = snapshot
    if not stable or not raw or not raw.endswith(b"\n"):
        if allow_incomplete:
            return None
        die("thermal stream has an incomplete trailing record: {}".format(path))
    try:
        text = raw.decode("ascii")
        reader = csv.DictReader(io.StringIO(text, newline=""))
        if tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
            die("thermal header mismatch: {}".format(path))
        rows = list(reader)
    except (UnicodeError, csv.Error) as exc:
        die("cannot parse thermal stream {}: {}".format(path, exc))
    if any(set(row) != set(THERMAL_FIELDS) or any(value is None for value in row.values())
           for row in rows):
        die("thermal stream contains a malformed row")
    return rows, sha256_bytes(raw)


def read_thermal_snapshot(
    path: Path,
    *,
    allow_incomplete: bool,
) -> Optional[List[Dict[str, str]]]:
    snapshot = read_thermal_snapshot_with_sha256(
        path, allow_incomplete=allow_incomplete)
    return None if snapshot is None else snapshot[0]


def read_thermal_rows(path: Path) -> List[Dict[str, str]]:
    return read_thermal_rows_with_sha256(path)[0]


def read_thermal_rows_with_sha256(
    path: Path,
) -> Tuple[List[Dict[str, str]], str]:
    snapshot = read_thermal_snapshot_with_sha256(
        path, allow_incomplete=False)
    assert snapshot is not None
    return snapshot


def bounded_thermal_sha256(path: Path) -> Optional[str]:
    snapshot = read_bounded_regular_bytes(
        path, allow_missing=True, max_bytes=MAX_THERMAL_STREAM_BYTES,
        label="thermal stream")
    if snapshot is None:
        return None
    raw, stable = snapshot
    if not stable:
        die("thermal stream changed during terminal hashing: {}".format(path))
    return sha256_bytes(raw)


def bounded_thermal_stderr_sha256(path: Path) -> Optional[str]:
    snapshot = read_bounded_regular_bytes(
        path, allow_missing=True, max_bytes=MAX_THERMAL_STDERR_BYTES,
        label="thermal sampler stderr")
    if snapshot is None:
        return None
    raw, stable = snapshot
    if not stable:
        die("thermal sampler stderr changed during terminal hashing: {}"
            .format(path))
    return sha256_bytes(raw)


def read_live_thermal_rows(path: Path) -> List[Dict[str, str]]:
    deadline = time.monotonic() + 0.5
    while True:
        rows = read_thermal_snapshot(path, allow_incomplete=True)
        if rows is not None:
            return rows
        if time.monotonic() >= deadline:
            die("live thermal stream remained incomplete: {}".format(path))
        time.sleep(0.005)


def read_sampler_pid_publication(
    pid_file: Path,
    previous: Optional[SamplerPidPublication] = None,
) -> Tuple[Optional[int], Optional[SamplerPidPublication]]:
    """Read one append-only sampler PID publication without parsing a prefix."""
    opened_here = False
    if previous is not None:
        if not isinstance(previous, tuple) or len(previous) != 4 or \
                any(not isinstance(value, int) for value in previous[:3]) or \
                not isinstance(previous[3], bytes):
            die("thermal sampler PID publication tracker is malformed")
        descriptor = previous[0]
    else:
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        try:
            descriptor = os.open(str(pid_file), flags)
        except FileNotFoundError:
            return None, None
        except OSError as exc:
            die("cannot open thermal sampler PID publication: {}".format(exc))
        opened_here = True
    try:
        identity = os.fstat(descriptor)
        if not stat.S_ISREG(identity.st_mode) or identity.st_nlink != 1 or \
                identity.st_uid not in (0, os.geteuid()) or \
                identity.st_mode & (stat.S_IWGRP | stat.S_IWOTH):
            die("thermal sampler PID publication is not a safe regular file")
        os.lseek(descriptor, 0, os.SEEK_SET)
        chunks: List[bytes] = []
        remaining = MAX_SAMPLER_PID_RECEIPT_BYTES + 1
        while remaining > 0:
            chunk = os.read(descriptor, remaining)
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
        try:
            path_identity = pid_file.lstat()
        except FileNotFoundError:
            die("thermal sampler PID publication disappeared during read")
        def safe_publication_identity(value: os.stat_result) -> bool:
            return stat.S_ISREG(value.st_mode) and value.st_nlink == 1 and \
                value.st_uid in (0, os.geteuid()) and not (
                    value.st_mode & (stat.S_IWGRP | stat.S_IWOTH))

        if not safe_publication_identity(after) or \
                not safe_publication_identity(path_identity) or \
                (identity.st_dev, identity.st_ino) != \
                (after.st_dev, after.st_ino) or \
                (after.st_dev, after.st_ino) != \
                (path_identity.st_dev, path_identity.st_ino) or \
                (
                    identity.st_uid, identity.st_gid,
                    stat.S_IMODE(identity.st_mode), identity.st_nlink) != (
                    after.st_uid, after.st_gid,
                    stat.S_IMODE(after.st_mode), after.st_nlink) or \
                (
                    after.st_uid, after.st_gid,
                    stat.S_IMODE(after.st_mode), after.st_nlink) != (
                    path_identity.st_uid, path_identity.st_gid,
                    stat.S_IMODE(path_identity.st_mode),
                    path_identity.st_nlink):
            die("thermal sampler PID publication identity changed")
        if len(raw) > MAX_SAMPLER_PID_RECEIPT_BYTES:
            die("thermal sampler PID publication is oversized")
        publication = (
            descriptor, int(after.st_dev), int(after.st_ino), raw)
        if previous is not None:
            if publication[:3] != previous[:3]:
                die("thermal sampler PID publication was replaced")
            if len(raw) < len(previous[3]) or \
                    not raw.startswith(previous[3]):
                die("thermal sampler PID publication did not grow append-only")

        stable_size = after.st_size == len(raw)
        if raw == b"":
            return None, publication
        if re.fullmatch(rb"[1-9][0-9]*", raw) is not None:
            if int(raw) > MAX_PROCESS_PID:
                die("thermal sampler PID publication prefix is out of range")
            return None, publication
        if re.fullmatch(rb"[1-9][0-9]*\n", raw) is None or not stable_size:
            die("thermal sampler PID file is malformed")
        pid = int(raw[:-1])
        if pid <= 1 or pid > MAX_PROCESS_PID:
            die("thermal sampler PID is out of range")
        return pid, publication
    except OSError as exc:
        if opened_here:
            try:
                os.close(descriptor)
            except OSError:
                pass
        die("cannot read thermal sampler PID publication: {}".format(exc))
    except BaseException:
        if opened_here:
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def close_sampler_pid_publication(
    publication: Optional[SamplerPidPublication],
) -> None:
    if publication is None:
        return
    try:
        os.close(publication[0])
    except OSError as exc:
        die("cannot close thermal sampler PID publication: {}".format(exc))


def sampler_identity(pid_file: Path, sampler_path: Path, csv_path: Path) -> Optional[int]:
    publication: Optional[SamplerPidPublication] = None
    try:
        pid, publication = read_sampler_pid_publication(pid_file)
        if pid is None:
            return None
        if not process_alive(pid):
            return None
        tokens = process_tokens(pid)
        if tokens != campaign_sampled_command(sampler_path, csv_path, pid_file):
            if not process_alive(pid):
                return None
            die("thermal sampler PID is not bound to this segment")
        return pid
    finally:
        close_sampler_pid_publication(publication)


def pidfd_has_exited(pidfd: int) -> bool:
    poller = select.poll()
    poller.register(pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
    return bool(poller.poll(0))


def bind_benchmark_exec_or_terminal(
    process: subprocess.Popen,
    pidfd: int,
    wrapper_argv: Sequence[str],
    benchmark_argv: Sequence[str],
    abort: threading.Event,
) -> Tuple[int, bool]:
    """Bind a direct child until its benchmark exec or pinned terminal state.

    The benchmark can legitimately exec and finish before /proc exposes its
    post-exec argv to this thread.  Its unreaped direct-child object and pidfd
    still pin one PID/start identity, and its exact terminal stdout/stderr/exit
    receipts are validated later.  Treat that terminal transition as safe
    evidence rather than inventing a live-identity failure.  A process that is
    still live must continue to prove either the exact controller-bound wrapper
    or the exact frozen benchmark command.

    Returns ``(start_ticks, exec_observed)``.  ``exec_observed`` is diagnostic;
    both outcomes are safe for PID-receipt publication because the latter is
    already terminal and therefore cannot become signal authority.
    """
    start_ticks = process_start_ticks(process.pid)
    if start_ticks is None:
        die("cannot bind benchmark direct-child start identity")
    deadline = time.monotonic() + BENCHMARK_EXEC_READY_TIMEOUT_SECONDS
    while time.monotonic() < deadline:
        observed_start = process_start_ticks(process.pid)
        tokens = process_tokens(process.pid)
        if observed_start != start_ticks:
            # A direct child remains unreaped throughout this function, so its
            # PID cannot be reused.  Losing or changing /proc start identity is
            # never a benign completion shortcut.
            die("benchmark direct-child start identity changed")
        if tokens == list(benchmark_argv):
            # Repeat both live identity fields.  At either boundary, pidfd
            # readability proves the same pinned object became terminal and
            # makes a transient/empty cmdline harmless.
            if pidfd_has_exited(pidfd):
                return start_ticks, True
            rebound_start = process_start_ticks(process.pid)
            rebound_tokens = process_tokens(process.pid)
            if rebound_start != start_ticks:
                die("benchmark direct-child start identity changed")
            if rebound_tokens == list(benchmark_argv):
                if pidfd_has_exited(pidfd):
                    return start_ticks, True
                return start_ticks, True
            if pidfd_has_exited(pidfd):
                return start_ticks, True
            if rebound_tokens:
                die("benchmark identity changed during live exec binding")
            # /proc/PID/cmdline can be transiently empty while exec replaces
            # the argument pages.  The pidfd still says this object is live;
            # retry rather than misclassifying that transition.
            time.sleep(0.005)
            continue
        if pidfd_has_exited(pidfd):
            return start_ticks, False
        if tokens and tokens != list(wrapper_argv):
            die("benchmark wrapper changed command unexpectedly")
        if abort.is_set():
            die("benchmark launch aborted before exec")
        time.sleep(0.005)
    die("benchmark wrapper did not exec before timeout")


def bind_sampler_identity(
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
    controller_identity: Tuple[int, int],
    wrapper_identity: Tuple[int, int, Sequence[str]],
    publication: Optional[SamplerPidPublication] = None,
) -> Optional[Tuple[int, int, List[str]]]:
    """Bind sampled child -> supervisor -> this unreaped sudo Popen."""
    owns_publication = publication is None
    publication_closed = False
    observed_publication: Optional[SamplerPidPublication] = None

    def close_owned_publication() -> None:
        nonlocal publication_closed
        if owns_publication and not publication_closed:
            publication_closed = True
            close_sampler_pid_publication(observed_publication)

    try:
        pid, observed_publication = read_sampler_pid_publication(
            pid_file, publication)
        if pid is None:
            close_owned_publication()
            return None
        if publication is not None and observed_publication != publication:
            die("thermal sampler PID publication changed before binding")
        assert observed_publication is not None
        wrapper_pid, expected_wrapper_start, expected_wrapper_tokens = \
            wrapper_identity
        exact_supervisor, exact_sudo_wrapper = \
            campaign_sampler_supervisor_commands(
                sampler_path, csv_path, pid_file, controller_identity)
        expected_wrapper_tokens = list(expected_wrapper_tokens)
        if not isinstance(wrapper_pid, int) or \
                isinstance(wrapper_pid, bool) or \
                wrapper_pid <= 1 or wrapper_pid > MAX_PROCESS_PID or \
                not isinstance(expected_wrapper_start, int) or \
                isinstance(expected_wrapper_start, bool) or \
                expected_wrapper_start < 0 or \
                expected_wrapper_tokens not in (
                    exact_supervisor, exact_sudo_wrapper):
            die("thermal sampler wrapper identity is malformed")
    except BaseException:
        close_owned_publication()
        raise
    child_pidfd: Optional[int] = None
    supervisor_pidfd: Optional[int] = None
    wrapper_pidfd: Optional[int] = None
    try:
        try:
            wrapper_pidfd = os.pidfd_open(wrapper_pid, 0)
        except ProcessLookupError:
            die("thermal sampler wrapper exited during identity binding")
        if pidfd_has_exited(wrapper_pidfd):
            die("thermal sampler wrapper exited during identity binding")
        wrapper_start = process_start_ticks(wrapper_pid)
        wrapper_tokens = process_tokens(wrapper_pid)
        if wrapper_start != expected_wrapper_start or \
                wrapper_tokens not in (
                    exact_supervisor, exact_sudo_wrapper):
            die("thermal sampler wrapper identity changed")
        try:
            child_pidfd = os.pidfd_open(pid, 0)
        except ProcessLookupError:
            return None
        if pidfd_has_exited(child_pidfd):
            return None
        expected_child = campaign_sampled_command(
            sampler_path, csv_path, pid_file)
        child_start = process_start_ticks(pid)
        child_parent = process_parent_pid(pid)
        child_tokens = process_tokens(pid)
        if child_start is None or child_parent is None or child_parent <= 1 or \
                child_tokens != expected_child:
            die("cannot bind exact live thermal sampler identity")
        supervisor_pid = child_parent
        if supervisor_pid == wrapper_pid:
            supervisor_pidfd = wrapper_pidfd
        else:
            try:
                supervisor_pidfd = os.pidfd_open(supervisor_pid, 0)
            except ProcessLookupError:
                die("thermal sampler supervisor exited during identity binding")
        supervisor_start = process_start_ticks(supervisor_pid)
        supervisor_tokens = process_tokens(supervisor_pid)
        supervisor_parent = process_parent_pid(supervisor_pid)
        if supervisor_start is None or supervisor_tokens != exact_supervisor:
            die("thermal sampler parent is not the controller-bound supervisor")
        if supervisor_pid == wrapper_pid:
            # sudo may exec directly.  Permit only its exact one-way argv
            # transition to the supervisor command in the same PID/start.
            if wrapper_tokens != exact_supervisor or \
                    expected_wrapper_tokens not in (
                        exact_sudo_wrapper, exact_supervisor):
                die("thermal sampler direct supervisor topology is invalid")
        elif wrapper_tokens != exact_sudo_wrapper or \
                expected_wrapper_tokens != exact_sudo_wrapper or \
                supervisor_parent != wrapper_pid:
            die("thermal sampler forked supervisor topology is invalid")
        # Repeat every numeric-PID and ancestry proof while all object handles
        # are held.  The sampler launch has no controlling TTY, so the
        # supported shape is wrapper -> supervisor -> sampled child; direct
        # sudo exec remains safe when wrapper == supervisor.
        stable_pid, stable_publication = read_sampler_pid_publication(
            pid_file, observed_publication)
        if stable_pid != pid or \
                stable_publication != observed_publication or \
                process_start_ticks(pid) != child_start or \
                process_parent_pid(pid) != supervisor_pid or \
                process_tokens(pid) != expected_child or \
                process_start_ticks(supervisor_pid) != supervisor_start or \
                process_tokens(supervisor_pid) != exact_supervisor or \
                (supervisor_pid != wrapper_pid and
                 process_parent_pid(supervisor_pid) != wrapper_pid) or \
                process_start_ticks(wrapper_pid) != expected_wrapper_start or \
                process_tokens(wrapper_pid) != wrapper_tokens or \
                pidfd_has_exited(child_pidfd) or \
                pidfd_has_exited(supervisor_pidfd) or \
                pidfd_has_exited(wrapper_pidfd):
            die("thermal sampler identity changed during stable binding")
        return pid, child_start, expected_child
    except OSError as exc:
        die("cannot bind thermal sampler identity: {}".format(exc))
    finally:
        descriptors: List[int] = []
        for descriptor in (
                supervisor_pidfd, child_pidfd, wrapper_pidfd):
            if descriptor is not None and descriptor not in descriptors:
                descriptors.append(descriptor)
        close_errors: List[str] = []
        for descriptor in descriptors:
            try:
                os.close(descriptor)
            except OSError as exc:
                close_errors.append("{}:{!r}".format(descriptor, exc))
        try:
            close_owned_publication()
        except BaseException as exc:
            close_errors.append("publication:{!r}".format(exc))
        if close_errors:
            die("thermal sampler pidfd cleanup failed: {}".format(
                "; ".join(close_errors)))


def wait_sampler_ready(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
    controller_identity: Tuple[int, int],
    wrapper_identity: Tuple[int, int, Sequence[str]],
    identity_sink: List[Tuple[int, int, List[str]]],
) -> Tuple[int, List[Dict[str, str]]]:
    deadline = time.monotonic() + THERMAL_READY_TIMEOUT_SECONDS
    publication: Optional[SamplerPidPublication] = None
    try:
        while time.monotonic() < deadline:
            if sampler.poll() is not None:
                die("thermal sampler exited before readiness")
            published_pid, publication = read_sampler_pid_publication(
                pid_file, publication)
            identity = None
            if published_pid is not None:
                identity = bind_sampler_identity(
                    pid_file, sampler_path, csv_path, controller_identity,
                    wrapper_identity, publication)
            pid = None if identity is None else identity[0]
            if identity is not None:
                if identity_sink and identity_sink[0] != identity:
                    die("launched thermal sampler identity changed before "
                        "readiness")
                if not identity_sink:
                    identity_sink.append(identity)
            if pid is not None and csv_path.is_file():
                rows = read_thermal_snapshot(
                    csv_path, allow_incomplete=True)
                if rows is not None and len(rows) >= THERMAL_READY_SAMPLES:
                    last = decimal_field(rows[-1], "monotonic_s")
                    age = Decimal(str(time.monotonic())) - last
                    if 0 <= age <= THERMAL_MAX_GAP_SECONDS:
                        return pid, rows
            time.sleep(0.05)
        die("thermal sampler did not provide two live data rows")
    finally:
        close_sampler_pid_publication(publication)


def wait_thermal_coverage(
    sampler: subprocess.Popen,
    csv_path: Path,
    end_monotonic_s: float,
) -> bool:
    deadline = time.monotonic() + THERMAL_END_TIMEOUT_SECONDS
    target = Decimal(str(end_monotonic_s))
    while time.monotonic() < deadline:
        if sampler.poll() is not None:
            return False
        rows = read_live_thermal_rows(csv_path)
        if rows and decimal_field(rows[-1], "monotonic_s") >= target:
            return True
        time.sleep(0.05)
    return False


def stop_launched_sampler(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
    sampled_identity: Optional[Tuple[int, int, Sequence[str]]],
) -> Tuple[Optional[int], List[dict]]:
    actions: List[dict] = []
    cleanup_errors: List[str] = []
    if sampled_identity is not None:
        pid, start_ticks, tokens = sampled_identity
        try:
            action = terminate_verified_process(pid, start_ticks, tokens)
            actions.append({
                "kind": "thermal", "pid": pid,
                "action": action,
            })
            if action not in STOPPED_PROCESS_ACTIONS:
                cleanup_errors.append("thermal:{}".format(action))
        except BaseException as exc:
            cleanup_errors.append("thermal_exception:{!r}".format(exc))
    returncode: Optional[int] = None
    try:
        returncode = sampler.wait(timeout=15.0)
    except subprocess.TimeoutExpired:
        try:
            action = terminate_privileged_direct_unreaped_child(sampler)
            actions.append({
                "kind": "thermal_wrapper", "pid": sampler.pid,
                "action": action,
            })
            if action not in STOPPED_PROCESS_ACTIONS:
                cleanup_errors.append(
                    "thermal_wrapper:{}".format(action))
        except BaseException as exc:
            cleanup_errors.append(
                "thermal_wrapper_exception:{!r}".format(exc))
        try:
            returncode = sampler.wait(timeout=5.0)
        except BaseException as exc:
            cleanup_errors.append(
                "thermal_wrapper_wait_exception:{!r}".format(exc))
    except BaseException as exc:
        cleanup_errors.append("thermal_wrapper_wait_exception:{!r}".format(exc))
        try:
            action = terminate_privileged_direct_unreaped_child(sampler)
            actions.append({
                "kind": "thermal_wrapper", "pid": sampler.pid,
                "action": action,
            })
            if action not in STOPPED_PROCESS_ACTIONS:
                cleanup_errors.append(
                    "thermal_wrapper:{}".format(action))
        except BaseException as cleanup_exc:
            cleanup_errors.append(
                "thermal_wrapper_exception:{!r}".format(cleanup_exc))
    # If readiness never captured the child, refuse to bless a PID read only
    # during cleanup.  It may already have been reused; never signal it.
    if sampled_identity is None:
        try:
            unknown_pid = sampler_identity(pid_file, sampler_path, csv_path)
            if unknown_pid is not None:
                cleanup_errors.append(
                    "uncaptured_live_sampler:{}".format(unknown_pid))
        except BaseException as exc:
            cleanup_errors.append(
                "uncaptured_sampler_inspection:{!r}".format(exc))
    remaining_readers: Optional[List[dict]] = None
    try:
        remaining_readers = other_samplers()
        if remaining_readers:
            cleanup_errors.append(
                "remaining_i2c_readers:{}".format(json.dumps(
                    remaining_readers, sort_keys=True)))
    except BaseException as exc:
        cleanup_errors.append("post_shutdown_i2c_check:{!r}".format(exc))
    if remaining_readers == [] and not cleanup_errors:
        try:
            pid_file.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            cleanup_errors.append("pid_file_cleanup:{!r}".format(exc))
    if cleanup_errors:
        die("campaign sampler shutdown failed closed: {}".format(
            "; ".join(cleanup_errors)))
    return returncode, actions


def cleanup_failed_sampler_binding(sampler: subprocess.Popen) -> None:
    """Best-effort all cleanup legs after wrapper identity capture fails."""
    errors: List[str] = []
    try:
        action = terminate_privileged_direct_unreaped_child(sampler)
        require_stopped_action(action, "sampler wrapper binding")
    except BaseException as exc:
        errors.append("wrapper_stop:{!r}".format(exc))
    try:
        remaining_readers = other_samplers()
        if remaining_readers:
            errors.append("remaining_i2c_readers:{}".format(json.dumps(
                remaining_readers, sort_keys=True)))
    except BaseException as exc:
        errors.append("post_shutdown_i2c_check:{!r}".format(exc))
    if errors:
        die("failed sampler binding cleanup was not proved: {}".format(
            "; ".join(errors)))


def benchmark_future_result(
    future: concurrent.futures.Future,
    job: int,
) -> dict:
    """Return one job-bound worker result without collapsing equal failures."""
    if future.cancelled():
        return {"job": job, "status": "future_cancelled"}
    try:
        result = future.result()
    except BaseException as exc:
        return {
            "job": job, "status": "worker_exception",
            "detail": repr(exc),
        }
    if not isinstance(result, dict) or result.get("job") != job:
        return {
            "job": job, "status": "worker_exception",
            "detail": "worker returned a mismatched job identity",
        }
    return result


def run_segment(
    root: Path,
    design: Mapping[str, object],
    tasks: Sequence[dict],
    stage: str,
    previously_complete: int,
    resume: bool,
) -> Dict[str, object]:
    indices = segment_indices(root)
    segment = indices[-1] + 1 if indices else 0
    controller_start_ticks = process_start_ticks(os.getpid())
    if controller_start_ticks is None:
        die("cannot bind campaign controller start identity")
    jobs = sorted(int(task["job"]) for task in tasks)
    intent = {
        "schema": SCHEMA + ".segment_intent", "segment": segment,
        "boot_id": required_boot_id(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "created_monotonic_s": time.monotonic(),
        "controller_pid": os.getpid(),
        "controller_start_ticks": controller_start_ticks,
        "stage": stage, "jobs": jobs,
        "jobs_sha256": sha256_bytes(json_lines(jobs)),
        "workers": design["worker_count"], "previously_complete": previously_complete,
        "resume": bool(resume),
        "retry_policy": "stage-atomic non-selective retry",
    }
    write_once(segment_path(root, segment, "intent"), canonical_json(intent))
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    make_private_directory(attempt_dir)
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    sampler_path = root / "frozen" / SAMPLER_NAME
    sampler_error = thermal_stderr.open("xb")
    sampled_command = campaign_sampled_command(
        sampler_path, thermal_csv, thermal_pid)
    sampler_command = [
        str(SUDO_PATH), "-n", str(PYTHON_PATH), "-I", "-c",
        PRIVILEGED_SAMPLER_SUPERVISOR_CODE,
        str(os.getpid()), str(controller_start_ticks),
        PDEATH_EXEC_CODE.encode("utf-8").hex(),
        *sampled_command,
    ]
    validate_trusted_root_tool(SUDO_PATH)
    validate_trusted_root_tool(PYTHON_PATH)
    try:
        sampler = launch_thermal_sampler(sampler_command, sampler_error)
    except OSError as exc:
        sampler_error.close()
        return write_segment_final(
            root, segment, intent, "failed",
            failures=[{"status": "sampler_spawn_error", "detail": str(exc)}],
            rolled_back=[], process_actions=[], jobs_ended_monotonic_s=None,
            sampler_returncode=None,
        )
    binding_error: Optional[BaseException] = None
    sampler_wrapper_start: Optional[int] = None
    sampler_wrapper_tokens: Optional[List[str]] = None
    try:
        sampler_wrapper_start = process_start_ticks(sampler.pid)
        sampler_wrapper_tokens = process_tokens(sampler.pid)
    except BaseException as exc:
        binding_error = exc
    if binding_error is None and (
            sampler_wrapper_start is None or not sampler_wrapper_tokens):
        binding_error = CampaignError(
            "cannot bind launched sampler wrapper identity")
    if binding_error is not None:
        try:
            cleanup_failed_sampler_binding(sampler)
        finally:
            sampler_error.close()
        raise binding_error
    assert sampler_wrapper_start is not None
    assert sampler_wrapper_tokens is not None
    sampled_identities: List[Tuple[int, int, List[str]]] = []
    abort = threading.Event()
    # The registry owns every direct benchmark child's pidfd from immediately
    # after spawn until wait/poll proves that child exited.  A worker exception
    # must never discard the coordinator's last safe cleanup handle.
    active: Dict[
        int, Tuple[subprocess.Popen, Optional[int], threading.Lock]
    ] = {}
    emergency_active: List[Tuple[int, subprocess.Popen]] = []
    active_lock = threading.Lock()
    failures: List[dict] = []
    process_actions: List[dict] = []
    jobs_started: Optional[float] = None
    jobs_ended: Optional[float] = None
    sampler_returncode: Optional[int] = None
    sealed_thermal_sha256: Optional[str] = None
    ready_written = False
    sampler_cleanup_proved = False

    def register_active(
        job: int,
        process: subprocess.Popen,
    ) -> threading.Lock:
        lease = threading.Lock()
        with active_lock:
            if job in active:
                die("benchmark job is already registered active")
            active[job] = (process, None, lease)
        return lease

    def attach_active_pidfd(
        job: int,
        process: subprocess.Popen,
        pidfd: int,
    ) -> None:
        with active_lock:
            current = active.get(job)
            if current is None or current[0] is not process or \
                    current[1] is not None:
                die("benchmark active-registry identity changed")
            active[job] = (process, pidfd, current[2])

    def release_stopped_active(
        job: int,
        process: subprocess.Popen,
        lease: threading.Lock,
        *,
        lease_held: bool = False,
    ) -> bool:
        """Drop one registry entry only after its child has been reaped."""
        def release() -> bool:
            owned_pidfd: Optional[int] = None
            with active_lock:
                current = active.get(job)
                if current is None:
                    # Registry removal is allowed only after a successful
                    # cleanup/reap while holding this same lease.
                    return process.returncode is not None
                if current[0] is not process or current[2] is not lease:
                    die("benchmark active-registry process changed")
                if process.returncode is None:
                    return False
                _, owned_pidfd, _ = active.pop(job)
            if owned_pidfd is not None:
                os.close(owned_pidfd)
            return True

        if lease_held:
            return release()
        with lease:
            return release()

    def terminate_registered(
        job: int,
        process: subprocess.Popen,
    ) -> Optional[str]:
        """Attempt one lease without removing it on any cleanup failure."""
        duplicate: Optional[int] = None
        with active_lock:
            current = active.get(job)
            if current is None:
                if process.returncode is None:
                    die("live benchmark lost its active-registry lease")
                return None
            if current[0] is not process:
                die("benchmark active-registry process changed")
            lease = current[2]
        with lease:
            with active_lock:
                current = active.get(job)
                if current is None:
                    if process.returncode is None:
                        die("live benchmark lost its active-registry lease")
                    return None
                if current[0] is not process or current[2] is not lease:
                    die("benchmark active-registry process changed")
                if current[1] is None:
                    # This window exists only between registry publication and
                    # pidfd attachment (or after pidfd_open failed).  Hold
                    # both locks across PID signaling/reap so the worker
                    # cannot attach or publish a reused numeric identity.
                    action = terminate_direct_unreaped_child(process)
                    require_stopped_action(
                        action, "active benchmark job {}".format(job))
                    if process.returncode is None:
                        die("benchmark cleanup action did not prove child exit")
                    active.pop(job)
                    return action
                duplicate = os.dup(current[1])
            try:
                action = terminate_direct_child_by_pidfd(process, duplicate)
                require_stopped_action(
                    action, "active benchmark job {}".format(job))
                if process.returncode is None:
                    die("benchmark cleanup action did not prove child exit")
                if not release_stopped_active(
                        job, process, lease, lease_held=True):
                    die("benchmark cleanup could not release stopped lease")
                return action
            finally:
                os.close(duplicate)

    def terminate_active() -> None:
        """Try every current child; retain failed leases and aggregate errors."""
        with active_lock:
            snapshot = [
                (job, process) for job, (process, _, _) in active.items()
            ]
            emergency_snapshot = list(emergency_active)
        results, cleanup_errors = attempt_all_cleanups(
            snapshot, terminate_registered)
        for job, process, action in results:
            if action is not None:
                process_actions.append({
                    "kind": "benchmark", "job": job, "pid": process.pid,
                    "action": action,
                })

        def terminate_emergency(
            job: int,
            process: subprocess.Popen,
        ) -> str:
            action = terminate_direct_unreaped_child(process)
            require_stopped_action(
                action, "emergency benchmark job {}".format(job))
            if process.returncode is None:
                die("emergency benchmark cleanup did not prove child exit")
            with active_lock:
                try:
                    emergency_active.remove((job, process))
                except ValueError:
                    die("emergency benchmark registry changed")
            return action

        emergency_results, emergency_errors = attempt_all_cleanups(
            emergency_snapshot, terminate_emergency)
        for job, process, action in emergency_results:
            process_actions.append({
                "kind": "benchmark_emergency", "job": job,
                "pid": process.pid, "action": action,
            })
        cleanup_errors.extend(
            "emergency_{}".format(error) for error in emergency_errors)
        if cleanup_errors:
            die("benchmark cleanup failed closed: {}".format(
                "; ".join(cleanup_errors)))

    def run_one(task: dict) -> dict:
        job = int(task["job"])
        if abort.is_set():
            return {"job": job, "status": "cancelled_before_start"}
        started_monotonic = time.monotonic()
        wrapper_argv = [
            str(PYTHON_PATH), "-I", "-c", PDEATH_EXEC_CODE,
            str(os.getpid()), str(controller_start_ticks), *task["argv"],
        ]
        try:
            process = subprocess.Popen(
                wrapper_argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                start_new_session=True,
            )
        except OSError as exc:
            return {"job": job, "status": "spawn_error", "detail": str(exc)}
        try:
            lease = register_active(job, process)
        except BaseException as registration_error:
            try:
                action = terminate_direct_unreaped_child(process)
                require_stopped_action(
                    action, "failed benchmark registration job {}"
                    .format(job))
            except BaseException as cleanup_error:
                with active_lock:
                    emergency_active.append((job, process))
                raise CampaignError(
                    "benchmark registration and fallback cleanup failed: "
                    "registration={!r}; cleanup={!r}".format(
                        registration_error, cleanup_error)
                ) from cleanup_error
            raise
        timed_out = False
        stdout = stderr = b""
        try:
            # The lease serializes exec/start/cmdline proof and receipt
            # publication against every cleanup path.  The child stays
            # unreaped while the pidfd and repeated /proc proofs bind the
            # numeric PID to the exact exec'd benchmark object.
            with lease:
                with active_lock:
                    current = active.get(job)
                    if current is None or current[0] is not process or \
                            current[1] is not None or \
                            current[2] is not lease or \
                            process.returncode is not None:
                        die("benchmark stopped before pidfd binding")
                try:
                    pidfd = os.pidfd_open(process.pid, 0)
                except (OSError, ProcessLookupError) as exc:
                    die("cannot bind benchmark pidfd: {}".format(exc))
                try:
                    attach_active_pidfd(job, process, pidfd)
                except BaseException:
                    os.close(pidfd)
                    raise
                start_ticks, exec_observed = bind_benchmark_exec_or_terminal(
                    process, pidfd, wrapper_argv, task["argv"], abort)
                with active_lock:
                    current = active.get(job)
                    if current is None or current[0] is not process or \
                            current[1] != pidfd or current[2] is not lease:
                        die("benchmark active identity changed before receipt")
                if not pidfd_has_exited(pidfd):
                    if not exec_observed:
                        die("live benchmark lacks an observed exec identity")
                    receipt_start = process_start_ticks(process.pid)
                    receipt_tokens = process_tokens(process.pid)
                    if receipt_start != start_ticks:
                        die("benchmark start identity changed before receipt")
                    if receipt_tokens != task["argv"] and \
                            not pidfd_has_exited(pidfd):
                        die("benchmark command identity changed before receipt")
                atomic_result(
                    attempt_dir / "job{:05d}.pid".format(job),
                    canonical_json({
                        "pid": process.pid, "start_ticks": start_ticks,
                    }),
                )
            if abort.is_set():
                action = terminate_registered(job, process)
                if action is not None:
                    process_actions.append({
                        "kind": "benchmark", "job": job, "pid": process.pid,
                        "action": action,
                    })
            try:
                stdout, stderr = process.communicate(
                    timeout=JOB_TIMEOUT_SECONDS)
            except subprocess.TimeoutExpired:
                timed_out = True
                action = terminate_registered(job, process)
                if action is not None:
                    process_actions.append({
                        "kind": "benchmark", "job": job, "pid": process.pid,
                        "action": action,
                    })
                stdout, stderr = process.communicate(timeout=15.0)
        except BaseException as original:
            cleanup_error: Optional[BaseException] = None
            try:
                action = terminate_registered(job, process)
                if action is not None:
                    process_actions.append({
                        "kind": "benchmark", "job": job,
                        "pid": process.pid, "action": action,
                    })
            except BaseException as exc:
                cleanup_error = exc
            if process.returncode is not None:
                try:
                    process.communicate(timeout=15.0)
                except (OSError, subprocess.TimeoutExpired) as exc:
                    if cleanup_error is None:
                        cleanup_error = exc
            if cleanup_error is not None:
                raise CampaignError(
                    "benchmark job {} failed and cleanup was not proved: "
                    "original={!r}; cleanup={!r}".format(
                        job, original, cleanup_error)) from cleanup_error
            raise
        finally:
            release_stopped_active(job, process, lease)
        ended_monotonic = time.monotonic()
        exit_bytes = (str(process.returncode) + "\n").encode("ascii")
        for suffix, data in (("stdout", stdout), ("stderr", stderr), ("exit", exit_bytes)):
            atomic_result(attempt_dir / "job{:05d}.{}".format(job, suffix), data)
        if timed_out:
            return {
                "job": job, "status": "timeout",
                "timeout_seconds": JOB_TIMEOUT_SECONDS,
            }
        if process.returncode != 0:
            return {"job": job, "status": "exit", "returncode": process.returncode}
        try:
            parse_benchmark_csv(stdout, task)
            if stderr:
                die("successful benchmark wrote stderr")
        except CampaignError as exc:
            return {"job": job, "status": "validation_error", "detail": str(exc)}
        if abort.is_set():
            return {"job": job, "status": "cancelled_after_run"}
        name = task["output_name"]
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        atomic_result(raw_path, stdout)
        atomic_result(stderr_path, stderr)
        atomic_result(exit_path, exit_bytes)
        receipt = {
            "schema": SCHEMA + ".job", "job": job, "segment": segment,
            "stage": stage, "output_name": name, "returncode": 0,
            "started_monotonic_s": started_monotonic,
            "ended_monotonic_s": ended_monotonic,
            "stdout_sha256": sha256_bytes(stdout),
            "stderr_sha256": sha256_bytes(stderr),
            "exit_sha256": sha256_bytes(exit_bytes),
        }
        write_once(job_receipt_path(root, job), canonical_json(receipt))
        return {"job": job, "status": "success"}

    state = "failed"
    caught: Optional[BaseException] = None
    try:
        sampled_pid, ready_rows = wait_sampler_ready(
            sampler, thermal_pid, sampler_path, thermal_csv,
            (os.getpid(), controller_start_ticks),
            (sampler.pid, sampler_wrapper_start, sampler_wrapper_tokens),
            sampled_identities,
        )
        if len(sampled_identities) != 1 or \
                sampled_identities[0][0] != sampled_pid:
            die("thermal sampler readiness lacks one captured identity")
        _, sampled_start_ticks, sampled_tokens = sampled_identities[0]
        validate_thermal_rows(ready_rows, segment=segment)

        def check_owned_sampler() -> Dict[int, List[str]]:
            assert_process_identity(
                sampled_pid, sampled_start_ticks, sampled_tokens,
                "launched thermal sampler")
            return sole_sampler_competitors(sampled_pid, sampler.pid)

        competitors = check_owned_sampler()
        if competitors:
            die("another I2C reader appeared during sampler startup: {}".format(
                json.dumps(competitors, sort_keys=True)))
        assert_current_controller_owner(root)
        jobs_started = time.monotonic()
        ready = {
            "schema": SCHEMA + ".segment_ready", "segment": segment,
            "intent_sha256": sha256_file(segment_path(root, segment, "intent")),
            "ready_utc": datetime.now(timezone.utc).isoformat(),
            "jobs_started_monotonic_s": jobs_started,
            "sampler_pid": sampled_pid, "samples_at_ready": len(ready_rows),
            "sampler_start_ticks": sampled_start_ticks,
            "last_sample_monotonic_s": ready_rows[-1]["monotonic_s"],
        }
        write_once(segment_path(root, segment, "ready"), canonical_json(ready))
        ready_written = True
        # Refresh the environment-ownership marker so cleanup tooling can
        # prove the sole I2C sampler of this segment belongs to a live
        # campaign (janitor schema wirehair.environment_owner.v1).
        write_owner_marker(
            root, "executing",
            [owner_protected_entry("controller", os.getpid()),
             {"role": "reader", "pid": sampled_pid,
              "start_ticks": sampled_start_ticks}],
            int(design["owner_ttl_hours"]),
        )
        next_i2c_check = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=int(design["worker_count"])) as pool:
            pending = {
                pool.submit(run_one, task): int(task["job"])
                for task in tasks
            }
            while pending:
                now = time.monotonic()
                if now >= next_i2c_check and not abort.is_set():
                    readers = check_owned_sampler()
                    validate_thermal_rows(
                        read_live_thermal_rows(thermal_csv), segment=segment)
                    if readers:
                        failures.append({
                            "status": "concurrent_i2c_reader",
                            "readers": readers,
                        })
                        abort.set()
                        terminate_active()
                    assert_current_controller_owner(root)
                    next_i2c_check = now + I2C_RECHECK_INTERVAL_SECONDS
                if sampler.poll() is not None and not abort.is_set():
                    failures.append({"status": "sampler_exited", "returncode": sampler.returncode})
                    abort.set()
                    terminate_active()
                done, _ = concurrent.futures.wait(
                    tuple(pending), timeout=0.2,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    job = pending.pop(future)
                    result = benchmark_future_result(future, job)
                    if result.get("status") != "success":
                        failures.append(result)
                        if not abort.is_set():
                            abort.set()
                            terminate_active()
                if abort.is_set():
                    terminate_active()
                    for future in pending:
                        future.cancel()
        jobs_ended = time.monotonic()
        if not failures:
            check_owned_sampler()
        if not failures and not wait_thermal_coverage(sampler, thermal_csv, jobs_ended):
            failures.append({"status": "thermal_end_coverage_failed"})
        if not failures:
            readers = check_owned_sampler()
            final_rows = read_live_thermal_rows(thermal_csv)
            final_thermal = validate_thermal_rows(
                final_rows, segment=segment)
            assert jobs_started is not None and jobs_ended is not None
            validate_successful_thermal_coverage(
                final_rows, final_thermal,
                start=Decimal(str(jobs_started)),
                end=Decimal(str(jobs_ended)),
                stage=stage, segment=segment)
            if readers:
                failures.append({
                    "status": "concurrent_i2c_reader",
                    "readers": readers,
                })
            assert_current_controller_owner(root)
        if not failures:
            state = "success"
    except BaseException as exc:
        caught = exc
        abort.set()
        try:
            terminate_active()
        except BaseException as cleanup_exc:
            failures.append({
                "status": "benchmark_stop_error",
                "detail": repr(cleanup_exc),
            })
        failures.append({"status": "launcher_exception", "detail": repr(exc)})
        jobs_ended = time.monotonic() if ready_written else None
    finally:
        try:
            sampler_returncode, sampler_actions = stop_launched_sampler(
                sampler, thermal_pid, sampler_path, thermal_csv,
                sampled_identities[0] if sampled_identities else None,
            )
            process_actions.extend(sampler_actions)
            sampler_cleanup_proved = True
        except BaseException as exc:
            if caught is None:
                caught = exc
            failures.append({"status": "sampler_stop_error", "detail": repr(exc)})
            state = "failed"
        sampler_error.close()
    try:
        terminate_active()
    except BaseException as exc:
        if caught is None:
            caught = exc
        failures.append({
            "status": "benchmark_stop_error", "detail": repr(exc),
        })
        state = "failed"
    with active_lock:
        unresolved_benchmarks = [
            {"job": job, "pid": process.pid}
            for job, (process, _, _) in active.items()
        ]
        unresolved_benchmarks.extend(
            {"job": job, "pid": process.pid, "emergency": True}
            for job, process in emergency_active
        )
    if unresolved_benchmarks:
        rollback_segment_outputs(root, tasks)
        die("benchmark cleanup was not proved; segment remains incomplete: {}"
            .format(json.dumps(unresolved_benchmarks, sort_keys=True)))
    if not sampler_cleanup_proved:
        rollback_segment_outputs(root, tasks)
        if caught is not None:
            raise caught
        die("campaign sampler cleanup was not proved")
    if sampler_returncode != 0:
        failures.append({"status": "sampler_returncode", "returncode": sampler_returncode})
        state = "failed"
    if not thermal_stderr.is_file() or thermal_stderr.stat().st_size:
        failures.append({"status": "sampler_stderr_nonempty_or_missing"})
        state = "failed"
    if state == "success":
        assert jobs_started is not None and jobs_ended is not None
        try:
            # The live checks above can race the sampler's final write.  Only
            # the bytes observed after the writer has exited are immutable and
            # therefore safe to bless in the segment-final hash.
            sealed_rows, candidate_thermal_sha256 = \
                read_thermal_rows_with_sha256(thermal_csv)
            sealed_thermal = validate_thermal_rows(
                sealed_rows, segment=segment)
            validate_successful_thermal_coverage(
                sealed_rows, sealed_thermal,
                start=Decimal(str(jobs_started)),
                end=Decimal(str(jobs_ended)),
                stage=stage, segment=segment)
            sealed_thermal_sha256 = candidate_thermal_sha256
        except Exception as exc:
            sealed_thermal_sha256 = None
            failures.append({
                "status": "sealed_thermal_validation_error",
                "detail": repr(exc),
            })
            state = "failed"
    rolled_back: List[int] = []
    if state != "success":
        rolled_back = rollback_segment_outputs(root, tasks)
    final = write_segment_final(
        root, segment, intent, state, failures=failures,
        rolled_back=rolled_back, process_actions=process_actions,
        jobs_ended_monotonic_s=jobs_ended,
        sampler_returncode=sampler_returncode,
        validated_thermal_csv_sha256=sealed_thermal_sha256,
    )
    if caught is not None:
        raise caught
    return final


def command_launch(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    lock_descriptor = acquire_environment_lock()
    if lock_descriptor is None:
        if args.preflight_only:
            print(json.dumps({
                "status": "BLOCKED_ACTIVE_CAMPAIGN_CONTROLLER",
                "root": str(root),
                "requested_stage": args.stage,
            }, indent=2, sort_keys=True))
            return
        die("another campaign controller holds the environment-owner lock")
    try:
        command_launch_locked(args, root)
    finally:
        os.close(lock_descriptor)


def command_launch_locked(args: argparse.Namespace, root: Path) -> None:
    design, tasks = verify_root(
        root, args.expected_receipts_sha256, runtime_tolerant=True,
    )
    owner = load_active_owner_marker()
    owner_live = (
        active_owner_live_identities(owner) if owner is not None else [])
    same_root_owner = (
        owner if owner is not None and owner.get("campaign_root") == str(root)
        else None)
    foreign_owner = (
        owner if owner is not None and owner.get("campaign_root") != str(root)
        else None)
    receipt_partials = sorted(
        path for directory in (
            root / "segments", root / "job_receipts", root / "raw",
            root / "stderr", root / "exit", root / "attempts",
        )
        for path in directory.rglob("*.part.*")
    )
    indices = segment_indices(root)
    incomplete = [
        segment for segment in indices
        if not segment_path(root, segment, "final").is_file()
    ]
    needs_reconciliation = bool(incomplete or receipt_partials)
    live_controller = any(
        entry["role"] == "controller" for entry in owner_live)
    orphan_reader_resume = bool(
        same_root_owner is not None and args.resume and needs_reconciliation and
        owner_live and not live_controller and
        all(entry["role"] == "reader" for entry in owner_live))
    if args.preflight_only:
        samplers = other_samplers()
        complete_count: Optional[int] = None
        if not needs_reconciliation:
            complete_count = len(completed_jobs(root, tasks))
        owner_scope = (
            "same_root" if same_root_owner is not None
            else "foreign" if foreign_owner is not None else None)
        print(json.dumps({
            "status": (
                "RESUME_ORPHAN_RECONCILIATION_REQUIRED"
                if orphan_reader_resume
                else "BLOCKED_SAME_ROOT_LIVE_OWNER"
                if same_root_owner is not None and owner_live
                else "BLOCKED_FOREIGN_LIVE_OWNER"
                if foreign_owner is not None and owner_live
                else "BLOCKED_FOREIGN_ENVIRONMENT_OWNER" if foreign_owner
                else "RESUME_STALE_OWNER_REQUIRED"
                if same_root_owner is not None and not args.resume
                else "RESUME_RECONCILIATION_REQUIRED"
                if needs_reconciliation and not args.resume
                else "RESUME_EXISTING_HISTORY_REQUIRED"
                if indices and not args.resume
                else "BLOCKED_CONCURRENT_SAMPLER" if samplers
                else "READY_NOT_LAUNCHED"),
            "root": str(root), "tasks": len(tasks), "complete": complete_count,
            "remaining": None if complete_count is None else len(tasks) - complete_count,
            "workers": design["worker_count"], "requested_stage": args.stage,
            "single_sampler": not bool(samplers), "other_samplers": samplers,
            "environment_owner": (
                None if owner is None else {
                    "campaign_root": owner.get("campaign_root"),
                    "phase": owner.get("phase"),
                    "expires_utc": owner.get("expires_utc"),
                    "scope": owner_scope,
                    "live_protected": owner_live,
                }),
            "incomplete_segments": incomplete,
            "orphan_receipt_partials": [str(path.relative_to(root)) for path in receipt_partials],
        }, indent=2, sort_keys=True))
        return
    if owner_live and not orphan_reader_resume:
        die("environment is owned by live protected process(es): {}".format(
            json.dumps({
                "campaign_root": owner.get("campaign_root") if owner else None,
                "live_protected": owner_live,
            }, sort_keys=True)))
    if foreign_owner is not None:
        die("environment is owned by another campaign: {}".format(
            json.dumps({
                "campaign_root": foreign_owner.get("campaign_root"),
                "expires_utc": foreign_owner.get("expires_utc"),
            }, sort_keys=True)))
    if same_root_owner is not None and not args.resume:
        die("stale same-root owner requires explicit --resume")
    if needs_reconciliation and not args.resume:
        die("interrupted segment(s) require --resume reconciliation")
    if indices and not args.resume:
        die("existing campaign segment history requires --resume")
    reconciled: List[dict] = []
    if args.resume:
        reconciled = reconcile_incomplete_segments(root, tasks)
        # Reconciliation must consume every live-process identity before any
        # fsynced partial receipt is discarded.  Attempt PID partials are
        # removed by the segment sealer only after exact-command cleanup.
        for path in receipt_partials:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    complete = completed_jobs(root, tasks)
    remaining = [task for task in tasks if int(task["job"]) not in complete]
    by_stage = {
        stage: [task for task in remaining if task["stage"] == stage]
        for stage in ("hard", "control")
    }
    samplers = other_samplers()
    if samplers:
        die("refusing concurrent SPD/I2C reader(s): {}".format(
            json.dumps(samplers, sort_keys=True)))
    if complete and not args.resume:
        die("outputs already exist; inspect then pass --resume")
    if not remaining:
        terminal_owner = load_active_owner_marker()
        terminal_owner_live = (
            active_owner_live_identities(terminal_owner)
            if terminal_owner is not None else [])
        if terminal_owner is not None and \
                terminal_owner.get("campaign_root") != str(root):
            die("environment owner changed during complete-ledger recovery")
        if terminal_owner is not None and terminal_owner_live:
            die("complete-ledger recovery still has live protected processes")
        if terminal_owner is not None and args.resume:
            # The controller may have died after the last successful segment
            # final but before replacing its executing owner marker.  Re-prove
            # the complete terminal ledger and the empty sampler inventory,
            # then finish that one missing lifecycle transition.
            validate_thermal(root, design, tasks)
            write_owner_marker(
                root, "complete", [], int(design["owner_ttl_hours"]))
            print(json.dumps({
                "status": "ALREADY_COMPLETE_OWNER_FINALIZED",
                "root": str(root), "tasks": len(tasks),
                "reconciled_segments": len(reconciled),
            }, indent=2, sort_keys=True))
            return
        die("campaign ledger is already complete")
    hard_incomplete = bool(by_stage["hard"])
    if args.stage == "control" and hard_incomplete:
        die("control stage requires the complete hard-stage Cartesian product")
    requested_stages = ("hard", "control") if args.stage == "all" else (args.stage,)
    if not any(by_stage[stage] for stage in requested_stages):
        die("requested campaign stage is already complete")
    # Declare environment ownership before any job or sampler starts, so
    # cleanup tooling never tears down this campaign's sealed environment.
    write_owner_marker(
        root, "executing",
        [owner_protected_entry("controller", os.getpid())],
        int(design["owner_ttl_hours"]),
    )
    assert_current_controller_owner(root)
    launched: List[dict] = []
    for stage in requested_stages:
        stage_tasks = by_stage[stage]
        if not stage_tasks:
            continue
        assert_current_controller_owner(root)
        # Recheck immediately before each privileged I2C sampler launch.
        samplers = other_samplers()
        if samplers:
            die("refusing concurrent SPD/I2C reader(s): {}".format(
                json.dumps(samplers, sort_keys=True)))
        final = run_segment(
            root, design, stage_tasks, stage, len(complete), bool(args.resume),
        )
        launched.append(final)
        print(json.dumps(final, indent=2, sort_keys=True))
        if final.get("state") != "success":
            # The marker deliberately stays in phase "executing" (bounded by
            # its TTL): a failed segment fails closed for cleanup tooling.
            raise SystemExit(1)
        complete.update(int(task["job"]) for task in stage_tasks)
    # Seal: no owned processes remain after a fully successful launch pass.
    assert_current_controller_owner(root)
    write_owner_marker(root, "complete", [], int(design["owner_ttl_hours"]))
    print(json.dumps({
        "status": "LAUNCH_OK", "reconciled_segments": len(reconciled),
        "launched_segments": [record["segment"] for record in launched],
        "complete": len(complete), "remaining": len(tasks) - len(complete),
    }, indent=2, sort_keys=True))


def validate_attempt_manifest(
    root: Path,
    segment: int,
    final: Mapping[str, object],
    allowed_jobs: Set[int],
    *,
    require_complete: bool,
) -> None:
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    if not attempt_dir.is_dir() or attempt_dir.is_symlink():
        die("segment attempt directory is missing or indirect")
    files = sorted(path for path in attempt_dir.iterdir() if path.is_file())
    entries = list(attempt_dir.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in entries):
        die("segment attempt inventory contains a symlink")
    allowed = re.compile(r"job([0-9]{5})\.(?:pid|stdout|stderr|exit)$")
    matches = [allowed.fullmatch(path.name) for path in files]
    if any(match is None for match in matches):
        die("segment attempt inventory contains an unexpected artifact")
    groups: Dict[int, Set[str]] = defaultdict(set)
    for path, match in zip(files, matches):
        assert match is not None
        job = int(match.group(1))
        if job not in allowed_jobs:
            die("segment attempt inventory names a job outside its intent")
        suffix = path.suffix[1:]
        groups[job].add(suffix)
        if suffix == "pid":
            identity = load_json(path)
            value = identity.get("pid")
            start_ticks = identity.get("start_ticks")
            if set(identity) != {"pid", "start_ticks"} or \
                    path.read_bytes() != canonical_json(identity) or \
                    not isinstance(value, int) or isinstance(value, bool) or value <= 1 or \
                    value > MAX_PROCESS_PID or \
                    not isinstance(start_ticks, int) or isinstance(start_ticks, bool) or \
                    start_ticks < 0:
                die("segment attempt PID receipt is out of range")
    complete_suffixes = {"pid", "stdout", "stderr", "exit"}
    valid_prefixes = (
        {"pid"},
        {"pid", "stdout"},
        {"pid", "stdout", "stderr"},
        complete_suffixes,
    )
    for job, suffixes in groups.items():
        # A controller crash can interrupt the three atomic output writes
        # after the PID receipt.  Preserve and authenticate any such prefix;
        # only a successful segment requires the complete four-file set.
        if suffixes not in valid_prefixes:
            die("segment attempt is not an atomic output prefix")
        if "exit" in suffixes:
            exit_path = attempt_dir / "job{:05d}.exit".format(job)
            try:
                exit_bytes = exit_path.read_bytes()
            except OSError:
                die("segment attempt exit receipt is malformed")
            if re.fullmatch(rb"(?:0|[1-9][0-9]*|-[1-9][0-9]*)\n",
                            exit_bytes) is None:
                die("segment attempt exit receipt is malformed")
            exit_value = int(exit_bytes[:-1])
            if exit_value < -(signal.NSIG - 1) or exit_value > 255:
                die("segment attempt exit receipt is out of range")
    if require_complete and (set(groups) != allowed_jobs or
                             any(value != complete_suffixes for value in groups.values())):
        die("successful segment lacks complete attempt evidence for every job")
    data = b"".join(
        "{}  {}\n".format(sha256_file(path), path.relative_to(root)).encode("ascii")
        for path in files
    )
    manifest = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
    if not manifest.is_file() or manifest.is_symlink() or manifest.read_bytes() != data or \
            final.get("attempt_manifest") != manifest.name or \
            final.get("attempt_manifest_sha256") != sha256_bytes(data) or \
            not isinstance(final.get("attempt_file_count"), int) or \
            isinstance(final.get("attempt_file_count"), bool) or \
            final.get("attempt_file_count") != len(files):
        die("segment attempt manifest mismatch")


def validate_thermal_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    segment: int,
) -> Dict[str, object]:
    if len(rows) < THERMAL_READY_SAMPLES:
        die("ready segment {} has fewer than two thermal rows".format(segment))
    monotonic: List[Decimal] = []
    busy: List[Decimal] = []
    cpu: List[Decimal] = []
    dimm_by_field: Dict[str, List[Decimal]] = {
        key: [] for key in THERMAL_FIELDS[5:13]
    }
    utc_previous: Optional[datetime] = None
    for row in rows:
        utc_text = row.get("utc")
        if not isinstance(utc_text, str) or re.fullmatch(
                r"[0-9]{4}-[0-9]{2}-[0-9]{2}T"
                r"[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}Z",
                utc_text) is None:
            die("thermal UTC timestamp is not in sampler format")
        try:
            utc = datetime.fromisoformat(utc_text.replace("Z", "+00:00"))
        except (KeyError, ValueError) as exc:
            die("invalid thermal UTC timestamp: {}".format(exc))
        if utc.utcoffset() != timedelta(0) or \
                (utc_previous is not None and utc <= utc_previous):
            die("thermal UTC timestamps are not strictly increasing")
        utc_previous = utc
        monotonic_text = row.get("monotonic_s")
        if not isinstance(monotonic_text, str) or re.fullmatch(
                r"(?:0|[1-9][0-9]*)\.[0-9]{6}",
                monotonic_text) is None:
            die("thermal monotonic timestamp is not in sampler format")
        current = decimal_field(row, "monotonic_s")
        if monotonic and current <= monotonic[-1]:
            die("thermal monotonic timestamps are not strictly increasing")
        if monotonic and not (
                THERMAL_MIN_GAP_SECONDS <= current - monotonic[-1] <=
                THERMAL_MAX_GAP_SECONDS):
            die("thermal cadence gap is outside policy in segment {}".format(segment))
        monotonic.append(current)
        if not row["cpu_busy_pct"]:
            die("thermal row lacks CPU busy utilization")
        value = decimal_field(row, "cpu_busy_pct")
        if value.is_signed() or value > 100:
            die("thermal CPU busy percentage is out of range")
        busy.append(value)
        for key in ("cpu_avg_mhz", "load1", "load5", "load15"):
            if not row[key]:
                die("thermal row lacks {}".format(key))
            auxiliary = decimal_field(row, key)
            if auxiliary.is_signed() or \
                    (key == "cpu_avg_mhz" and auxiliary == 0):
                die("thermal auxiliary metric {} is out of range".format(key))
        if not row["cpu_tctl_c"]:
            die("thermal row lacks CPU Tctl")
        cpu_temperature = decimal_field(row, "cpu_tctl_c")
        if not MIN_PLAUSIBLE_CPU_C < cpu_temperature < MAX_PLAUSIBLE_CPU_C:
            die("thermal CPU Tctl is physically implausible")
        cpu.append(cpu_temperature)
        for key in THERMAL_FIELDS[5:13]:
            if not row[key]:
                die("thermal row lacks DIMM field {}".format(key))
            dimm_temperature = decimal_field(row, key)
            if not MIN_PLAUSIBLE_DIMM_C < dimm_temperature < \
                    MAX_PLAUSIBLE_DIMM_C:
                die("thermal DIMM temperature is physically implausible")
            dimm_by_field[key].append(dimm_temperature)
        if any(row.get(key) != "0" for key in (
                "dimm_read_errors", "edac_ce", "edac_ue")):
            die("thermal hardware error receipt is nonzero")
    return {
        "samples": len(rows), "monotonic": monotonic, "busy": busy,
        "cpu": cpu, "dimm_by_field": dimm_by_field,
    }


def validate_successful_thermal_coverage(
    rows: Sequence[Mapping[str, str]],
    thermal: Mapping[str, object],
    *,
    start: Decimal,
    end: Decimal,
    stage: str,
    segment: int,
) -> List[Tuple[Decimal, Decimal]]:
    monotonic = thermal["monotonic"]
    if not isinstance(monotonic, list) or not monotonic or \
            end <= start or monotonic[0] > start or monotonic[-1] < end:
        die("successful thermal segment lacks full launch start/end coverage")
    # Each sampler row's cpu_busy_pct is the /proc/stat delta for the interval
    # ending at that row.  Weight only the overlap of interval (i-1, i] with
    # the benchmark window; an arithmetic row mean lets long low-utilization
    # gaps hide behind short high-utilization intervals.
    coverage_busy: List[Tuple[Decimal, Decimal]] = []
    for index in range(1, len(monotonic)):
        left = max(monotonic[index - 1], start)
        right = min(monotonic[index], end)
        if right <= left:
            continue
        coverage_busy.append((
            decimal_field(rows[index], "cpu_busy_pct"),
            right - left,
        ))
    covered_duration = sum(
        (duration for _, duration in coverage_busy), Decimal(0))
    if not coverage_busy or covered_duration != end - start:
        die("successful thermal segment has an interval coverage gap")
    weighted_mean = sum(
        (value * duration for value, duration in coverage_busy), Decimal(0)
    ) / covered_duration
    if stage in CPU_BUSY_FLOOR_STAGES and weighted_mean < CPU_BUSY_FLOOR:
        die("successful segment {} CPU busy mean is below the sealed floor"
            .format(segment))
    return coverage_busy


def validate_thermal(root: Path, design: Mapping[str, object], tasks: Sequence[dict]) -> Dict[str, object]:
    indices = segment_indices(root)
    if not indices:
        die("campaign has no launch segments")
    expected_segment_files: Set[str] = set()
    successful_jobs: Set[int] = set()
    successful_segments = failed_segments = interrupted_segments = 0
    samples = 0
    all_busy: List[Tuple[Decimal, Decimal]] = []
    all_cpu: List[Decimal] = []
    dimm_maxima: Dict[str, Decimal] = {}
    seen_control_success = False
    expected_thermal_files: Set[str] = set()
    task_by_job = {int(task["job"]): task for task in tasks}
    for segment in indices:
        intent_path = segment_path(root, segment, "intent")
        ready_path = segment_path(root, segment, "ready")
        final_path = segment_path(root, segment, "final")
        manifest_path = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
        expected_segment_files.update((intent_path.name, final_path.name, manifest_path.name))
        if not intent_path.is_file() or intent_path.is_symlink() or \
                not final_path.is_file() or final_path.is_symlink():
            die("campaign contains an unreconciled segment")
        intent = load_json(intent_path)
        final = load_json(final_path)
        selected = segment_jobs(intent, tasks)
        jobs = [int(task["job"]) for task in selected]
        stage = str(intent.get("stage"))
        expected_jobs = [
            int(task["job"]) for task in tasks
            if task["stage"] == stage and int(task["job"]) not in successful_jobs
        ]
        hard_remaining = any(
            task["stage"] == "hard" and int(task["job"]) not in successful_jobs
            for task in tasks
        )
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                not isinstance(intent.get("segment"), int) or \
                isinstance(intent.get("segment"), bool) or \
                intent.get("segment") != segment or \
                not valid_segment_boot_id(intent.get("boot_id")) or \
                not isinstance(intent.get("controller_pid"), int) or \
                isinstance(intent.get("controller_pid"), bool) or \
                intent["controller_pid"] <= 1 or \
                intent["controller_pid"] > MAX_PROCESS_PID or \
                not isinstance(intent.get("controller_start_ticks"), int) or \
                isinstance(intent.get("controller_start_ticks"), bool) or \
                intent["controller_start_ticks"] < 0 or \
                not isinstance(intent.get("workers"), int) or \
                isinstance(intent.get("workers"), bool) or \
                intent.get("workers") != design.get("worker_count") or \
                not isinstance(intent.get("previously_complete"), int) or \
                isinstance(intent.get("previously_complete"), bool) or \
                intent.get("previously_complete") != len(successful_jobs) or \
                jobs != expected_jobs or (stage == "control" and hard_remaining) or \
                intent.get("retry_policy") != "stage-atomic non-selective retry" or \
                final.get("schema") != SCHEMA + ".segment_final" or \
                not isinstance(final.get("segment"), int) or \
                isinstance(final.get("segment"), bool) or \
                final.get("segment") != segment or \
                final.get("stage") != stage or \
                final.get("intent_sha256") != sha256_file(intent_path) or \
                not isinstance(final.get("jobs"), list) or \
                any(not isinstance(job, int) or isinstance(job, bool)
                    for job in final.get("jobs", [])) or \
                final.get("jobs") != jobs or \
                final.get("jobs_sha256") != intent.get("jobs_sha256") or \
                final.get("retry_policy") != (
                    "entire stage segment; no successful survivor is retained "
                    "on failure"):
            die("segment intent/final binding mismatch")
        state = final.get("state")
        if state not in ("success", "failed", "interrupted"):
            die("unknown segment terminal state")
        rolled_back = final.get("rolled_back_jobs")
        published = final.get("published_jobs")
        failures = final.get("failures")
        if not isinstance(rolled_back, list) or \
                any(not isinstance(job, int) or isinstance(job, bool) or job not in jobs
                    for job in rolled_back) or rolled_back != sorted(set(rolled_back)) or \
                not isinstance(published, list) or \
                any(not isinstance(job, int) or isinstance(job, bool)
                    for job in published) or \
                not isinstance(failures, list):
            die("segment terminal retry ledger is malformed")
        validate_attempt_manifest(
            root, segment, final, set(jobs), require_complete=state == "success",
        )
        ready: Optional[Dict[str, object]] = None
        if ready_path.exists():
            expected_segment_files.add(ready_path.name)
            if ready_path.is_symlink() or not ready_path.is_file():
                die("segment ready receipt is indirect")
            ready = load_json(ready_path)
            samples_at_ready = ready.get("samples_at_ready")
            if ready.get("schema") != SCHEMA + ".segment_ready" or \
                    not isinstance(ready.get("segment"), int) or \
                    isinstance(ready.get("segment"), bool) or \
                    ready.get("segment") != segment or \
                    ready.get("intent_sha256") != sha256_file(intent_path) or \
                    final.get("ready_sha256") != sha256_file(ready_path) or \
                    not isinstance(samples_at_ready, int) or isinstance(samples_at_ready, bool) or \
                    samples_at_ready < THERMAL_READY_SAMPLES or \
                    not isinstance(ready.get("sampler_pid"), int) or \
                    isinstance(ready.get("sampler_pid"), bool) or ready["sampler_pid"] <= 1 or \
                    ready["sampler_pid"] > MAX_PROCESS_PID or \
                    not isinstance(ready.get("sampler_start_ticks"), int) or \
                    isinstance(ready.get("sampler_start_ticks"), bool) or \
                    ready["sampler_start_ticks"] < 0 or \
                    not isinstance(ready.get("jobs_started_monotonic_s"), (int, float)) or \
                    isinstance(ready.get("jobs_started_monotonic_s"), bool):
                die("segment readiness receipt mismatch")
        elif final.get("ready_sha256") is not None:
            die("segment final references a missing readiness receipt")
        csv_path = root / "thermal" / "segment{:03d}.csv".format(segment)
        stderr_path = root / "thermal" / "segment{:03d}.stderr".format(segment)
        pid_path = root / "thermal" / "segment{:03d}.pid".format(segment)
        if pid_path.exists() or pid_path.is_symlink():
            die("finalized segment retains a thermal PID file")
        csv_exists = csv_path.is_file() and not csv_path.is_symlink()
        stderr_exists = stderr_path.is_file() and not stderr_path.is_symlink()
        if csv_exists != (final.get("thermal_csv_sha256") is not None) or \
                stderr_exists != (final.get("thermal_stderr_sha256") is not None):
            die("segment thermal artifact cardinality mismatch")
        if ready is not None and (not csv_exists or not stderr_exists):
            die("ready segment lacks complete thermal artifacts")
        if csv_exists:
            expected_thermal_files.add(csv_path.name)
        if stderr_exists:
            expected_thermal_files.add(stderr_path.name)
            observed_stderr_sha256 = bounded_thermal_stderr_sha256(
                stderr_path)
            if final.get("thermal_stderr") != stderr_path.name or \
                    final.get("thermal_stderr_sha256") != \
                    observed_stderr_sha256 or \
                    state == "success" and observed_stderr_sha256 != \
                    sha256_bytes(b""):
                die("segment thermal stderr hash binding mismatch")
        thermal: Optional[Dict[str, object]] = None
        rows: List[Dict[str, str]] = []
        if csv_exists:
            if final.get("thermal_csv") != csv_path.name or \
                    (stderr_exists and final.get("thermal_stderr") != stderr_path.name):
                die("segment thermal hash binding mismatch")
            if state == "success":
                rows, observed_csv_sha256 = \
                    read_thermal_rows_with_sha256(csv_path)
            else:
                observed_csv_sha256 = bounded_thermal_sha256(csv_path)
            if final.get("thermal_csv_sha256") != observed_csv_sha256:
                die("segment thermal hash binding mismatch")
            if state == "success":
                thermal = validate_thermal_rows(rows, segment=segment)
                if ready is not None:
                    count = int(ready["samples_at_ready"])
                    if count > len(rows) or \
                            ready.get("last_sample_monotonic_s") != \
                            rows[count - 1]["monotonic_s"] or \
                            Decimal(str(
                                ready["jobs_started_monotonic_s"])) < \
                            decimal_field(
                                rows[count - 1], "monotonic_s") or \
                            Decimal(str(
                                ready["jobs_started_monotonic_s"])) - \
                            decimal_field(
                                rows[count - 1], "monotonic_s") > \
                            THERMAL_MAX_GAP_SECONDS:
                        die("segment readiness is not bound to its live "
                            "thermal prefix")
                samples += int(thermal["samples"])
                all_cpu.extend(thermal["cpu"])
                for key, values in thermal["dimm_by_field"].items():
                    dimm_maxima[key] = max(
                        dimm_maxima.get(key, values[0]), max(values))
        if state == "success":
            successful_segments += 1
            end_value = final.get("jobs_ended_monotonic_s")
            if ready is None or thermal is None or \
                    not isinstance(final.get("sampler_returncode"), int) or \
                    isinstance(final.get("sampler_returncode"), bool) or \
                    final.get("sampler_returncode") != 0 or \
                    final.get("failures") != [] or final.get("rolled_back_jobs") != [] or \
                    final.get("published_jobs") != jobs or \
                    not isinstance(end_value, (int, float)) or isinstance(end_value, bool):
                die("successful segment terminal receipt is incomplete")
            start = Decimal(str(ready["jobs_started_monotonic_s"]))
            end = Decimal(str(final["jobs_ended_monotonic_s"]))
            coverage_busy = validate_successful_thermal_coverage(
                rows, thermal, start=start, end=end, stage=stage,
                segment=segment)
            all_busy.extend(coverage_busy)
            if intent.get("stage") == "control":
                seen_control_success = True
            elif seen_control_success:
                die("successful hard segment appears after a successful control segment")
            overlap = successful_jobs.intersection(jobs)
            if overlap:
                die("successful job appears in more than one segment")
            for job in jobs:
                receipt = load_json(job_receipt_path(root, job))
                task = task_by_job[job]
                name = task["output_name"]
                if receipt.get("segment") != segment or \
                        Decimal(str(receipt["started_monotonic_s"])) < start or \
                        Decimal(str(receipt["ended_monotonic_s"])) > end or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.stdout".format(job)) != \
                            sha256_file(root / "raw" / name) or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.stderr".format(job)) != \
                            sha256_file(root / "stderr" / (name + ".stderr")) or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.exit".format(job)) != \
                            sha256_file(root / "exit" / (name + ".exit")):
                    die("job receipt is bound to the wrong successful segment")
            successful_jobs.update(jobs)
        else:
            if state == "failed":
                failed_segments += 1
            else:
                interrupted_segments += 1
            if final.get("published_jobs") != []:
                die("non-success segment retained published jobs")
            if state == "failed" and not failures:
                die("failed segment lacks a failure receipt")
            for job in jobs:
                receipt_path = job_receipt_path(root, job)
                if receipt_path.is_file() and load_json(receipt_path).get("segment") == segment:
                    die("non-success segment retained a job receipt")
    actual_segment_files = {path.name for path in (root / "segments").iterdir()}
    if actual_segment_files != expected_segment_files:
        die("segment receipt inventory contains missing or unexpected files")
    actual_thermal_files = {path.name for path in (root / "thermal").iterdir()}
    if actual_thermal_files != expected_thermal_files:
        die("thermal segment inventory contains missing or unexpected files")
    expected_attempt_dirs = {"segment{:03d}".format(index) for index in indices}
    attempt_entries = list((root / "attempts").iterdir())
    if {path.name for path in attempt_entries} != expected_attempt_dirs or \
            any(not path.is_dir() or path.is_symlink() for path in attempt_entries):
        die("attempt segment directory inventory mismatch")
    completed = completed_jobs(root, tasks)
    if successful_jobs != completed:
        die("successful segment/job receipt coverage mismatch")
    if completed != set(task_by_job):
        die("successful segment receipts do not cover the complete campaign")
    expected_job_receipts = {
        job_receipt_path(root, job).name for job in completed
    }
    job_receipt_entries = list((root / "job_receipts").iterdir())
    if {path.name for path in job_receipt_entries} != expected_job_receipts or \
            any(not path.is_file() or path.is_symlink() for path in job_receipt_entries):
        die("job receipt inventory mismatch")
    if not all_busy:
        die("campaign has no successful thermal load samples")
    busy_duration = sum(
        (duration for _, duration in all_busy), Decimal(0))
    if busy_duration <= 0:
        die("campaign has no positive successful thermal load duration")
    successful_busy_mean = sum(
        (value * duration for value, duration in all_busy), Decimal(0)
    ) / busy_duration
    return {
        "segments": len(indices), "successful_segments": successful_segments,
        "failed_segments": failed_segments,
        "interrupted_segments": interrupted_segments, "samples": samples,
        "successful_busy_mean": str(successful_busy_mean),
        "cpu_max_c": str(max(all_cpu)) if all_cpu else None,
        "dimm_max_c_by_field": {
            key: str(dimm_maxima[key]) for key in sorted(dimm_maxima)
        },
    }


def domain_iter(k_lo: int, k_hi: int) -> Iterator[Tuple[int, int, str]]:
    for K in range(k_lo, k_hi + 1):
        for seed_index in range(len(SEEDS)):
            for schedule in SCHEDULES:
                yield (K, seed_index, schedule)


def reduce_outcomes(
    outcomes: Mapping[Tuple[str, str, int, int, str], Mapping[str, object]],
    hard: Sequence[Mapping[str, object]],
    source_records: Mapping[str, Sequence[dict]],
    *,
    k_lo: int = K_LO,
    k_hi: int = K_HI,
) -> Dict[str, object]:
    """Score both survivors against their sealed canonical comparators.

    ``outcomes`` maps ``(arm, stage, K, seed_index, schedule)`` to a parsed
    cell record.  The 16 development hard cells are scored separately as
    fix-confirmation and excluded from the headline holdout comparison; the
    acceptance gate requires raw field-shortfall and raw total-recovery
    improvement, strict headline total-recovery improvement, at least one
    actual headline repair, zero headline introductions, a genuine panel
    source field-shortfall fix, and paired cost ratios within tolerance.
    """
    arms = arms_catalog()
    by_name = {arm["name"]: arm for arm in arms}
    hard_keys = [cell_tuple(cell) for cell in hard]
    hard_set = set(hard_keys)
    if len(hard_set) != len(hard):
        die("duplicate development hard cell in reduction")
    if any(not (k_lo <= key[0] <= k_hi) for key in hard_keys):
        die("development hard cell outside the reduced K domain")
    domain_count = (k_hi - k_lo + 1) * len(SEEDS) * len(SCHEDULES)
    if len(outcomes) != len(arms) * (domain_count + len(hard)):
        die("reduction outcome cardinality mismatch")
    panel_fs: Dict[str, Set[Tuple[int, int, str]]] = {}
    for panel in FINALIST_ARMS:
        rows = source_records.get(panel)
        if not isinstance(rows, (list, tuple)):
            die("source records lack finalist panel {}".format(panel))
        panel_fs[panel] = {
            cell_tuple(row) for row in rows
            if row.get("cause") == "field_shortfall"
        }

    def validate_reduced_cell(record: object, label: str) -> Mapping[str, object]:
        if not isinstance(record, Mapping):
            die("{} is not a reduced cell record".format(label))
        if record.get("outcome") not in ("success", "field_shortfall", "q>H"):
            die("{} has an unknown reduced outcome".format(label))
        for field in ("xors", "muladds"):
            value = record.get(field)
            if not isinstance(value, Decimal) or not value.is_finite() or \
                    value.is_signed():
                die("{} has an invalid nonnegative work metric {}".format(
                    label, field))
        for field in ("seed_attempt", "inactivated", "binary_deficit"):
            value = record.get(field)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                die("{} has an invalid nonnegative integer receipt {}".format(
                    label, field))
        return record

    arm_failures: Dict[str, Dict[Tuple[int, int, str], str]] = {}
    for arm in arms:
        failure_map: Dict[Tuple[int, int, str], str] = {}
        for key in domain_iter(k_lo, k_hi):
            record = outcomes.get((arm["name"], "control", *key))
            if record is None:
                die("missing control cell {} for arm {}".format(key, arm["name"]))
            record = validate_reduced_cell(
                record, "{} control {}".format(arm["name"], key))
            outcome = record["outcome"]
            if outcome != "success":
                failure_map[key] = str(outcome)
        arm_failures[arm["name"]] = failure_map
    # The canonical comparators must reproduce the sealed all-K source
    # results exactly (same failing cells with the same causes) before any
    # survivor is scored; this also re-proves the sealed totals.
    for arm in arms:
        if arm["role"] != "comparator":
            continue
        source_map = {
            cell_tuple(row): str(row["cause"])
            for row in source_records[arm["panel"]]
            if k_lo <= cell_tuple(row)[0] <= k_hi
        }
        if arm_failures[arm["name"]] != source_map:
            die("comparator arm {} failed exact all-K source replay".format(arm["name"]))
    # Hard-stage fix-confirmation replays must agree with the all-K sweep
    # cell-for-cell (chunking must never change an outcome), and comparator
    # hard replays must reproduce the sealed source expectations.
    for arm in arms:
        failure_map = arm_failures[arm["name"]]
        for key in hard_keys:
            record = outcomes.get((arm["name"], "hard", *key))
            if record is None:
                die("missing hard replay cell {} for arm {}".format(key, arm["name"]))
            record = validate_reduced_cell(
                record, "{} hard {}".format(arm["name"], key))
            if record["outcome"] != failure_map.get(key, "success"):
                die("hard-stage replay diverged from the all-K sweep for {} {}"
                    .format(arm["name"], key))
        if arm["role"] == "comparator":
            for key in hard_keys:
                expected = ("field_shortfall" if key in panel_fs[arm["panel"]]
                            else "success")
                if outcomes[(arm["name"], "hard", *key)]["outcome"] != expected:
                    die("comparator {} failed sealed hard-cell replay at {}"
                        .format(arm["name"], key))

    def summarize(name: str) -> Dict[str, object]:
        failure_map = arm_failures[name]
        causes: Dict[str, int] = defaultdict(int)
        by_seed: Dict[int, int] = defaultdict(int)
        by_schedule: Dict[str, int] = defaultdict(int)
        multiplicity: Dict[int, int] = defaultdict(int)
        failure_records: List[dict] = []
        for key in sorted(failure_map):
            cause = failure_map[key]
            causes[cause] += 1
            by_seed[key[1]] += 1
            by_schedule[key[2]] += 1
            multiplicity[key[0]] += 1
            failure_records.append({**cell_record(key), "cause": cause})
        return {
            "role": by_name[name]["role"], "panel": by_name[name]["panel"],
            "mask": by_name[name]["mask"], "mask_hex": by_name[name]["mask_hex"],
            "cells": domain_count, "failures": len(failure_map),
            "field_shortfalls": causes.get("field_shortfall", 0),
            "causes": dict(sorted(causes.items())),
            "failures_by_seed": {str(key): by_seed[key] for key in sorted(by_seed)},
            "failures_by_schedule": {
                key: by_schedule[key] for key in sorted(by_schedule)
            },
            "weak_K": len(multiplicity),
            "weak_K_multiplicity": {
                str(K): multiplicity[K] for K in sorted(multiplicity)
            },
            "max_K_multiplicity": max(multiplicity.values()) if multiplicity else 0,
            "failure_records": failure_records,
        }

    arm_summaries = {arm["name"]: summarize(arm["name"]) for arm in arms}
    candidates: Dict[str, dict] = {}
    for arm in arms:
        if arm["role"] != "candidate":
            continue
        comparator = by_name[str(arm["comparator"])]
        candidate_map = arm_failures[arm["name"]]
        comparator_map = arm_failures[comparator["name"]]
        repair_keys = sorted(key for key in comparator_map if key not in candidate_map)
        introduction_keys = sorted(key for key in candidate_map if key not in comparator_map)
        headline_repairs = [key for key in repair_keys if key not in hard_set]
        headline_introductions = [
            key for key in introduction_keys if key not in hard_set
        ]
        receipt_deltas = {
            field: {"cells": 0, "candidate_minus_comparator_sum": 0, "max_abs": 0}
            for field in ("seed_attempt", "inactivated", "binary_deficit")
        }
        candidate_work = [Decimal(0), Decimal(0)]
        comparator_work = [Decimal(0), Decimal(0)]
        common_success = 0
        for key in domain_iter(k_lo, k_hi):
            candidate_record = outcomes[(arm["name"], "control", *key)]
            comparator_record = outcomes[(comparator["name"], "control", *key)]
            # Configuration selection is inherited from the codec, with no
            # experiment-specific seed patches.  A mask can therefore alter
            # the selected attempt or solve dimensions; preserve those
            # paired deltas as results instead of rejecting them.
            for field, stats in receipt_deltas.items():
                delta = int(candidate_record[field]) - int(comparator_record[field])
                if delta:
                    stats["cells"] += 1
                    stats["candidate_minus_comparator_sum"] += delta
                    stats["max_abs"] = max(stats["max_abs"], abs(delta))
            if key in hard_set:
                continue
            if key not in candidate_map and key not in comparator_map:
                common_success += 1
                candidate_work[0] += candidate_record["xors"]
                candidate_work[1] += candidate_record["muladds"]
                comparator_work[0] += comparator_record["xors"]
                comparator_work[1] += comparator_record["muladds"]
        xor_ratio = (
            candidate_work[0] / comparator_work[0] if comparator_work[0] else None)
        muladd_ratio = (
            candidate_work[1] / comparator_work[1] if comparator_work[1] else None)
        candidate_fs = arm_summaries[arm["name"]]["field_shortfalls"]
        comparator_fs = arm_summaries[comparator["name"]]["field_shortfalls"]
        candidate_raw_failures = len(candidate_map)
        comparator_raw_failures = len(comparator_map)
        candidate_headline_fs = sum(
            cause == "field_shortfall" and key not in hard_set
            for key, cause in candidate_map.items())
        comparator_headline_fs = sum(
            cause == "field_shortfall" and key not in hard_set
            for key, cause in comparator_map.items())
        candidate_headline_failures = sum(
            key not in hard_set for key in candidate_map)
        comparator_headline_failures = sum(
            key not in hard_set for key in comparator_map)
        fix_rows: List[dict] = []
        fixed = residual = 0
        for key in hard_keys:
            candidate_outcome = candidate_map.get(key, "success")
            comparator_outcome = comparator_map.get(key, "success")
            in_panel = key in panel_fs[arm["panel"]]
            if in_panel:
                if candidate_outcome == "success":
                    fixed += 1
                else:
                    residual += 1
            fix_rows.append({
                **cell_record(key),
                "source_arms": sorted(
                    panel for panel in FINALIST_ARMS if key in panel_fs[panel]),
                "panel_source_failure": in_panel,
                "candidate": candidate_outcome,
                "comparator": comparator_outcome,
            })
        cost_within = (
            xor_ratio is not None and muladd_ratio is not None and
            xor_ratio <= COST_RATIO_MAX and muladd_ratio <= COST_RATIO_MAX)
        fix_confirmed = fixed > 0
        improved = candidate_fs < comparator_fs
        raw_recovery_improved = (
            candidate_raw_failures < comparator_raw_failures)
        recovery_improved = (
            candidate_headline_failures < comparator_headline_failures and
            bool(headline_repairs))
        no_new = not headline_introductions
        acceptance = {
            "raw_field_shortfalls": {
                "candidate": candidate_fs, "comparator": comparator_fs,
            },
            "headline_field_shortfalls": {
                "candidate": candidate_headline_fs,
                "comparator": comparator_headline_fs,
            },
            "gate_field_shortfall_improved": improved,
            "raw_total_failures": {
                "candidate": candidate_raw_failures,
                "comparator": comparator_raw_failures,
            },
            "gate_raw_total_recovery_improved": raw_recovery_improved,
            "headline_total_failures": {
                "candidate": candidate_headline_failures,
                "comparator": comparator_headline_failures,
            },
            "gate_actual_recovery_improved": recovery_improved,
            "headline_introductions": len(headline_introductions),
            "gate_zero_headline_introductions": no_new,
            "cost_ratio_max": str(COST_RATIO_MAX),
            "common_success_xor_ratio": (
                None if xor_ratio is None else str(xor_ratio)),
            "common_success_muladd_ratio": (
                None if muladd_ratio is None else str(muladd_ratio)),
            "gate_cost_within_tolerance": cost_within,
            "panel_source_failures_fixed": fixed,
            "gate_dev_fix_confirmed": fix_confirmed,
            "accepted": (
                improved and raw_recovery_improved and recovery_improved and
                no_new and cost_within and fix_confirmed),
        }
        candidates[arm["name"]] = {
            "comparator": comparator["name"],
            "repairs_vs_comparator": len(repair_keys),
            "introductions_vs_comparator": len(introduction_keys),
            "repair_keys_vs_comparator": [cell_record(key) for key in repair_keys],
            "introduction_keys_vs_comparator": [
                cell_record(key) for key in introduction_keys
            ],
            "headline_repairs": len(headline_repairs),
            "headline_introductions": len(headline_introductions),
            "headline_introduction_keys": [
                cell_record(key) for key in headline_introductions
            ],
            "paired_receipt_deltas_vs_comparator": receipt_deltas,
            "headline_common_success": common_success,
            "fix_confirmation": {
                "policy": (
                    "the development hard cells were selection material; "
                    "fix-confirmation only, excluded from the headline "
                    "holdout comparison"),
                "panel_source_failures_fixed": fixed,
                "panel_source_failures_residual": residual,
                "cells": fix_rows,
            },
            "acceptance": acceptance,
        }
    return {
        "schema": SCHEMA + ".reduction",
        "codec_errors": 0,
        "K_range": [k_lo, k_hi],
        "cells_per_arm_control": domain_count,
        "hard_cells": [dict(cell) for cell in hard],
        "arms": arm_summaries,
        "candidates": candidates,
        "accepted": all(
            block["acceptance"]["accepted"] for block in candidates.values()),
    }


def command_reduce(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    if len(completed_jobs(root, tasks)) != len(tasks):
        die("campaign is incomplete")
    thermal = validate_thermal(root, design, tasks)
    outcomes: Dict[Tuple[str, str, int, int, str], dict] = {}
    manifest_lines = bytearray()
    runtime_paths: List[Path] = []
    for directory_name in ("segments", "thermal", "attempts", "job_receipts"):
        directory = root / directory_name
        for path in directory.rglob("*"):
            if path.is_symlink():
                die("runtime evidence contains a symlink")
            if path.is_file():
                runtime_paths.append(path)
    for path in sorted(runtime_paths, key=lambda value: str(value.relative_to(root))):
        manifest_lines.extend("{}  {}\n".format(
            sha256_file(path), path.relative_to(root)).encode("ascii"))
    for task in tasks:
        name = task["output_name"]
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        for path in (raw_path, stderr_path, exit_path):
            manifest_lines.extend("{}  {}\n".format(
                sha256_file(path), path.relative_to(root)).encode("ascii"))
        rows = parse_benchmark_csv(raw_path.read_bytes(), task)
        for row in rows:
            try:
                K = int(row["N"])
                seed_attempt = int(row["seed_attempt"])
                inactivated = int(row["inact_max"])
                binary_deficit = int(row["binary_def_max"])
            except (KeyError, TypeError, ValueError) as exc:
                die("benchmark integer result field is malformed: {}".format(
                    exc))
            key = (str(task["arm"]), str(task["stage"]), K,
                   int(task["seed_index"]), str(task["schedule"]))
            if key in outcomes:
                die("duplicate Cartesian cell {}".format(key))
            outcomes[key] = {
                "outcome": classify(row),
                "xors": decimal_field(row, "block_xors_mu"),
                "muladds": decimal_field(row, "block_muladds_mu"),
                "seed_attempt": seed_attempt,
                "inactivated": inactivated,
                "binary_deficit": binary_deficit,
            }
    if len(outcomes) != int(design["total_cells"]):
        die("Cartesian result count mismatch")
    hard = read_jsonl(root / "hard_keys.jsonl")
    source_records = validate_source_summary(
        load_json(root / "frozen" / SOURCE_SUMMARY_NAME))
    reduction = reduce_outcomes(outcomes, hard, source_records)
    data_manifest = bytes(manifest_lines)
    data_manifest_sha256 = sha256_bytes(data_manifest)
    write_once(root / "data_manifest.sha256", data_manifest)
    summary = {
        "schema": SCHEMA + ".validated", "validation_issue_count": 0,
        "root": str(root), "prelaunch_receipts_sha256": args.expected_receipts_sha256,
        "data_manifest_sha256": data_manifest_sha256, "design": design,
        "thermal": thermal, "reduction": reduction,
    }
    write_once(root / "validated_summary.json", canonical_json(summary))
    print(json.dumps({
        "status": "VALIDATION_OK", "cells": len(outcomes),
        "accepted": reduction["accepted"],
        "acceptance_by_candidate": {
            name: block["acceptance"]["accepted"]
            for name, block in sorted(reduction["candidates"].items())
        },
        "validated_summary_sha256": sha256_file(root / "validated_summary.json"),
        "data_manifest_sha256": data_manifest_sha256,
    }, indent=2, sort_keys=True))


def make_thermal_bytes(
    monotonic_values: Sequence[str],
    *,
    edac_ce: str = "0",
    blank_dimm: bool = False,
    busy: str = "99.5",
    cpu_temperature: str = "60",
    dimm_temperature: str = "45",
    blank_busy_indices: Sequence[int] = (),
    row_overrides: Optional[Mapping[str, str]] = None,
) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=THERMAL_FIELDS, lineterminator="\n")
    writer.writeheader()
    for index, monotonic_value in enumerate(monotonic_values):
        row = {key: "0" for key in THERMAL_FIELDS}
        row.update({
            "utc": "2026-07-21T00:00:{:02d}.000Z".format(index),
            "monotonic_s": "{:.6f}".format(Decimal(monotonic_value)),
            "cpu_busy_pct": busy,
            "cpu_avg_mhz": "5000", "cpu_tctl_c": cpu_temperature,
            "dimm_read_errors": "0", "load1": "128", "load5": "128",
            "load15": "128", "edac_ce": edac_ce, "edac_ue": "0",
        })
        if index in blank_busy_indices:
            row["cpu_busy_pct"] = ""
        for key in THERMAL_FIELDS[5:13]:
            row[key] = dimm_temperature
        if blank_dimm:
            row[THERMAL_FIELDS[5]] = ""
        if row_overrides:
            row.update(row_overrides)
        writer.writerow(row)
    return output.getvalue().encode("ascii")


def make_success_fixture(
    root: Path,
    *,
    stage: str = "hard",
    busy: str = "99.5",
) -> Tuple[List[dict], dict]:
    """Synthesize one finalized successful segment for validator selftests."""
    make_private_directory(root, parents=True, exist_ok=True)
    for directory in RUNTIME_DIRECTORY_NAMES:
        make_private_directory(root / directory, parents=True, exist_ok=True)
    task = {"job": 0, "stage": stage, "output_name": "job.csv", "argv": ["/bin/true"]}
    intent = {
        "schema": SCHEMA + ".segment_intent", "segment": 0,
        "boot_id": required_boot_id(),
        "created_utc": "2026-07-21T00:00:00+00:00",
        "created_monotonic_s": 99.5,
        "controller_pid": 99999998, "controller_start_ticks": 1,
        "stage": stage, "jobs": [0],
        "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
        "previously_complete": 0, "resume": False,
        "retry_policy": "stage-atomic non-selective retry",
    }
    write_once(segment_path(root, 0, "intent"), canonical_json(intent))
    make_private_directory(root / "attempts" / "segment000")
    (root / "attempts" / "segment000" / "job00000.pid").write_bytes(
        canonical_json({"pid": 99999999, "start_ticks": 1}))
    ready = {
        "schema": SCHEMA + ".segment_ready", "segment": 0,
        "intent_sha256": sha256_file(segment_path(root, 0, "intent")),
        "ready_utc": "2026-07-21T00:00:00+00:00",
        "jobs_started_monotonic_s": 101.1, "sampler_pid": 99999999,
        "sampler_start_ticks": 1,
        "samples_at_ready": 2, "last_sample_monotonic_s": "101.000000",
    }
    write_once(segment_path(root, 0, "ready"), canonical_json(ready))
    thermal_csv = root / "thermal" / "segment000.csv"
    thermal_stderr = root / "thermal" / "segment000.stderr"
    thermal_csv.write_bytes(make_thermal_bytes(("100.0", "101.0", "102.0"), busy=busy))
    thermal_stderr.write_bytes(b"")
    raw = root / "raw" / "job.csv"
    stderr = root / "stderr" / "job.csv.stderr"
    exit_path = root / "exit" / "job.csv.exit"
    raw.write_bytes(b"synthetic output\n")
    stderr.write_bytes(b"")
    exit_path.write_bytes(b"0\n")
    (root / "attempts" / "segment000" / "job00000.stdout").write_bytes(raw.read_bytes())
    (root / "attempts" / "segment000" / "job00000.stderr").write_bytes(b"")
    (root / "attempts" / "segment000" / "job00000.exit").write_bytes(b"0\n")
    write_once(job_receipt_path(root, 0), canonical_json({
        "schema": SCHEMA + ".job", "job": 0, "segment": 0,
        "stage": stage, "output_name": "job.csv", "returncode": 0,
        "started_monotonic_s": 101.2, "ended_monotonic_s": 101.4,
        "stdout_sha256": sha256_file(raw), "stderr_sha256": sha256_file(stderr),
        "exit_sha256": sha256_file(exit_path),
    }))
    write_segment_final(
        root, 0, intent, "success", failures=[], rolled_back=[],
        process_actions=[], jobs_ended_monotonic_s=101.5,
        sampler_returncode=0,
        validated_thermal_csv_sha256=bounded_thermal_sha256(thermal_csv),
    )
    return [task], {"worker_count": 1}


def synth_cell(
    outcome: str,
    *,
    xors: Decimal = Decimal("1000"),
    muladds: Decimal = Decimal("100"),
    seed_attempt: int = 1,
    inactivated: int = 5,
    binary_deficit: int = 0,
) -> Dict[str, object]:
    return {
        "outcome": outcome, "xors": xors, "muladds": muladds,
        "seed_attempt": seed_attempt, "inactivated": inactivated,
        "binary_deficit": binary_deficit,
    }


def synthetic_reduction_fixture(
    k_lo: int = 2,
    k_hi: int = 13,
) -> Tuple[Dict[Tuple[str, str, int, int, str], dict], List[dict], Dict[str, List[dict]]]:
    """Build a small internally consistent accept-path reduction fixture.

    Every comparator field shortfall is one of its panel's development cells,
    matching the real sealed source.  Each comparator also has an independent
    headline q>H failure which its candidate repairs, so the fixture proves
    raw field-shortfall fix-confirmation and independent headline recovery.
    """
    h1 = (3, 0, "burst")
    h2 = (5, 1, "adversarial")
    h3 = (7, 2, "repair-only")
    h4 = (9, 0, "burst")
    q1 = (11, 1, "burst")
    q2 = (12, 2, "repair-only")
    source_records = {
        "p48_r3": [
            {**cell_record(h1), "cause": "field_shortfall"},
            {**cell_record(h2), "cause": "field_shortfall"},
            {**cell_record(q1), "cause": "q>H"},
        ],
        "p32_r7": [
            {**cell_record(h3), "cause": "field_shortfall"},
            {**cell_record(h4), "cause": "field_shortfall"},
            {**cell_record(q2), "cause": "q>H"},
        ],
    }
    hard = [
        {**cell_record(key), "source_arms": arms}
        for key, arms in sorted({
            h1: ["p48_r3"], h2: ["p48_r3"],
            h3: ["p32_r7"], h4: ["p32_r7"],
        }.items())
    ]
    failure_maps = {
        "p48_r3_sfx380": {
            h1: "field_shortfall", h2: "field_shortfall",
            q1: "q>H"},
        "p48_r3_pfx007": {h2: "field_shortfall"},
        "p32_r7_sfx3f8": {
            h3: "field_shortfall", h4: "field_shortfall",
            q2: "q>H"},
        "p32_r7_pfx07f": {},
    }
    outcomes: Dict[Tuple[str, str, int, int, str], dict] = {}
    hard_keys = [cell_tuple(cell) for cell in hard]
    for arm in arms_catalog():
        failure_map = failure_maps[arm["name"]]
        for key in domain_iter(k_lo, k_hi):
            outcomes[(arm["name"], "control", *key)] = synth_cell(
                failure_map.get(key, "success"))
        for key in hard_keys:
            outcomes[(arm["name"], "hard", *key)] = synth_cell(
                failure_map.get(key, "success"))
    return outcomes, hard, source_records


def expect_reject(function: Any, label: str) -> None:
    try:
        function()
    except CampaignError:
        return
    die("{} tamper selftest was accepted".format(label))


def command_selftest(_: argparse.Namespace) -> None:
    # Canonical-JSON round trips and immutability semantics.
    sample = {"b": [1, 2], "a": {"z": "x", "k": None}, "n": 7}
    encoded = canonical_json(sample)
    if json.loads(encoded.decode("ascii")) != sample or \
            not encoded.endswith(b"\n") or b" " in encoded:
        die("canonical JSON round-trip selftest failed")
    if canonical_json(json.loads(encoded.decode("ascii"))) != encoded:
        die("canonical JSON stability selftest failed")
    try:
        canonical_json({"bad": float("nan")})
    except ValueError:
        pass
    else:
        die("canonical JSON NaN rejection selftest failed")

    # Arm-catalog and mask validation.
    arms = arms_catalog()
    if [arm["name"] for arm in arms] != [arm["name"] for arm in ARMS] or \
            {arm["mask_hex"] for arm in arms} != {"0x007", "0x380", "0x07f", "0x3f8"} or \
            any(bin(arm["mask"]).count("1") != arm["rows"] for arm in arms) or \
            sum(arm["canonical_suffix"] for arm in arms) != 2:
        die("arm catalog selftest failed")
    bad_popcount = [dict(arm) for arm in ARMS]
    bad_popcount[0]["mask"] = 0x00F
    expect_reject(lambda: validate_arms(bad_popcount), "mask popcount")
    bad_range = [dict(arm) for arm in ARMS]
    bad_range[0]["mask"] = 1 << 10
    expect_reject(lambda: validate_arms(bad_range), "mask range")
    bad_suffix = [dict(arm) for arm in ARMS]
    bad_suffix[1]["mask"] = 0x1C0
    expect_reject(lambda: validate_arms(bad_suffix), "noncanonical comparator")
    bad_pair = [dict(arm) for arm in ARMS]
    bad_pair[0]["comparator"] = "p32_r7_sfx3f8"
    expect_reject(lambda: validate_arms(bad_pair), "cross-panel pairing")

    # Preamble parsing.
    if parse_preamble("# precodefail: a=1 b=two") != {"a": "1", "b": "two"}:
        die("preamble selftest failed")
    expect_reject(
        lambda: parse_preamble("# precodefail: duplicate=1 duplicate=2"),
        "duplicate preamble",
    )
    single_trial = {
        "success": "1", "rank_fail": "0", "error": "0",
        "fail_rate": "0.00000000",
        "inact_mu": "12.000", "inact_max": "12",
        "binary_def_mu": "12.000", "binary_def_max": "12",
        "heavy_gain_mu": "12.000", "heavy_gain_min": "12",
        "heavy_shortfall": "0", "first_rank_fail": "-1",
        "binary_def_hist": "12:1", "heavy_gain_hist": "12:1",
        "failure_trials": "",
    }
    if classify(single_trial) != "success":
        die("single-trial success semantics selftest failed")
    forged_trial = dict(single_trial)
    forged_trial.update({
        "inact_mu": "13.000", "inact_max": "13",
        "binary_def_mu": "13.000", "binary_def_max": "13",
        "binary_def_hist": "13:1",
    })
    expect_reject(
        lambda: classify(forged_trial),
        "single-trial contradictory success",
    )

    # The launch guard must identify renamed I2C readers from their retained
    # device descriptors, and recognize wrapped known samplers before their
    # descriptors are open without classifying unrelated Python commands.
    if parse_fuser_i2c_result(
            Path("/dev/i2c-1"), 0, b" 41 17 41",
            b"/dev/i2c-1:         \n") != (17, 41) or \
            parse_fuser_i2c_result(
                Path("/dev/i2c-2"), 1, b"", b"") != ():
        die("I2C reader parser selftest failed")
    expect_reject(
        lambda: parse_fuser_i2c_result(
            Path("/dev/i2c-1"), 0, b"41\n",
            b"/dev/i2c-1:         \n"),
        "malformed I2C reader output",
    )
    expect_reject(
        lambda: parse_fuser_i2c_result(
            Path("/dev/i2c-1"), 1, b"", b"warning\n"),
        "diagnostic I2C no-reader output",
    )
    expected_inventory = {
        71: ["/dev/i2c-1", "/dev/i2c-2"],
    }
    if validate_sole_sampler_inventory(
            expected_inventory, 71, 70) != {} or \
            validate_sole_sampler_inventory({
                **expected_inventory, 91: ["/dev/i2c-1"]}, 71, 70) != \
                {91: ["/dev/i2c-1"]}:
        die("sole sampler inventory selftest failed")
    expect_reject(
        lambda: validate_sole_sampler_inventory(
            {71: ["/dev/i2c-1"]}, 71, 70),
        "incomplete owned I2C inventory",
    )
    if not is_wirehair_sampler_command([
            "sudo", "-n", "python3",
            "/tmp/wirehair_expo_thermal_sampler_hardened.py",
            "--csv", "/tmp/thermal.csv", "--pid-file", "/tmp/thermal.pid",
    ]) or not is_wirehair_sampler_command([
            "python3", "/tmp/wirehair_expo_thermal_sampler.py",
            "--csv=/tmp/thermal.csv", "--pid-file=/tmp/thermal.pid",
    ]) or not is_wirehair_sampler_command([
            "bash", "-c",
            "python3 /tmp/wirehair_expo_thermal_sampler_hardened.py "
            "--csv /tmp/thermal.csv --pid-file /tmp/thermal.pid",
    ]) or is_wirehair_sampler_command([
            "python3", "/tmp/unrelated_csv_worker.py",
            "--csv", "/tmp/output.csv", "--pid-file", "/tmp/worker.pid",
    ]) or is_wirehair_sampler_command([
            "python3", "/tmp/review.py",
            "/tmp/wirehair_expo_thermal_sampler.py",
            "--csv", "notes", "--pid-file", "fixture",
    ]) or is_wirehair_sampler_command([
            "rg", "read_spd5118_temperature",
            "/tmp/wirehair_expo_thermal_sampler.py",
    ]):
        die("sampler command recognition selftest failed")

    # A process that exits and has its PID reused after TERM must never receive
    # the later KILL.  The identity transition is injected deterministically.
    termination_state = {"start_ticks": 99}
    sent_signals: List[str] = []

    def simulated_signal(signal_name: str) -> None:
        sent_signals.append(signal_name)
        if signal_name == "TERM":
            termination_state["start_ticks"] = 100

    action = terminate_verified_process(
        1234, 99, ["owned", "worker"],
        alive_probe=lambda _: True,
        start_ticks_probe=lambda _: termination_state["start_ticks"],
        tokens_probe=lambda _: ["owned", "worker"],
        signal_sender=simulated_signal,
    )
    if action != "identity_reused" or sent_signals != ["TERM"]:
        die("PID-reuse-safe termination selftest failed")
    sent_signals.clear()
    action = terminate_verified_process(
        1234, 99, ["owned", "worker"],
        alive_probe=lambda _: True,
        start_ticks_probe=lambda _: 99,
        tokens_probe=lambda _: ["foreign", "worker"],
        signal_sender=simulated_signal,
    )
    if action != "identity_changed" or sent_signals:
        die("changed-command termination selftest failed")

    # Task-ledger identity, cardinality, chunking, and argv shape.
    _, hard, _ = synthetic_reduction_fixture()
    tasks = build_tasks("/frozen/wirehair_v2_bench", hard, k_lo=2, k_hi=501)
    if [task["job"] for task in tasks] != list(range(len(tasks))) or \
            len({task["output_name"] for task in tasks}) != len(tasks) or \
            sum(len(task["Ks"]) for task in tasks) != 4 * (4 + 500 * 9):
        die("task-ledger identity/cardinality selftest failed")
    stages = [task["stage"] for task in tasks]
    if stages != sorted(stages, key=("hard", "control").index):
        die("task-ledger stage ordering selftest failed")
    if any(len(task["Ks"]) > CHUNK_MAX or not task["Ks"] for task in tasks):
        die("task-ledger chunk bound selftest failed")
    cartesian = {
        (task["arm"], task["stage"], K, task["seed_index"], task["schedule"])
        for task in tasks for K in task["Ks"]
    }
    if len(cartesian) != 4 * (4 + 500 * 9):
        die("task-ledger Cartesian uniqueness selftest failed")
    control_chunks = [
        task for task in tasks
        if task["stage"] == "control" and task["arm"] == "p48_r3_pfx007" and
        task["seed_index"] == 0 and task["schedule"] == "burst"
    ]
    if [len(task["Ks"]) for task in control_chunks] != [250, 250] or \
            control_chunks[0]["Ks"][0] != 2 or control_chunks[1]["Ks"][-1] != 501:
        die("task-ledger chunk layout selftest failed")
    probe = tasks[0]
    argv = probe["argv"]
    for flag, value in (
            ("--mixed-grouped-gf256-rows", str(probe["rows"])),
            ("--mixed-grouped-gf256-row-mask", "0x{:x}".format(probe["mask"])),
            ("--mixed-period", str(probe["period"])),
            ("--loss", "0.50"), ("--bb-list", "64"), ("--trials", "1")):
        if value != argv[argv.index(flag) + 1]:
            die("task argv flag selftest failed: {}".format(flag))

    # Reducer accept path on synthetic receipts.
    outcomes, hard, source_records = synthetic_reduction_fixture()
    reduction = reduce_outcomes(outcomes, hard, source_records, k_lo=2, k_hi=13)
    p48 = reduction["candidates"]["p48_r3_pfx007"]
    p32 = reduction["candidates"]["p32_r7_pfx07f"]
    if not reduction["accepted"] or reduction["codec_errors"] != 0 or \
            p48["acceptance"]["raw_field_shortfalls"] != {"candidate": 1, "comparator": 2} or \
            p32["acceptance"]["raw_field_shortfalls"] != {"candidate": 0, "comparator": 2} or \
            p48["acceptance"]["headline_field_shortfalls"] != \
                {"candidate": 0, "comparator": 0} or \
            not p48["acceptance"]["gate_actual_recovery_improved"] or \
            not p48["acceptance"]["gate_raw_total_recovery_improved"] or \
            p48["repairs_vs_comparator"] != 2 or p48["introductions_vs_comparator"] != 0 or \
            p32["repairs_vs_comparator"] != 3 or \
            p48["fix_confirmation"]["panel_source_failures_fixed"] != 1 or \
            p48["fix_confirmation"]["panel_source_failures_residual"] != 1 or \
            p48["acceptance"]["common_success_xor_ratio"] != "1" or \
            reduction["arms"]["p48_r3_sfx380"]["weak_K_multiplicity"] != \
                {"3": 1, "5": 1, "11": 1}:
        die("reducer accept selftest failed")

    # Reducer rejection paths: headline introduction, cost, tie, source
    # replay, hard-stage divergence, and cardinality.
    introduced = dict(outcomes)
    introduced[("p48_r3_pfx007", "control", 4, 0, "burst")] = synth_cell("q>H")
    verdict = reduce_outcomes(introduced, hard, source_records, k_lo=2, k_hi=13)
    if verdict["accepted"] or \
            verdict["candidates"]["p48_r3_pfx007"]["acceptance"]["gate_zero_headline_introductions"]:
        die("reducer headline-introduction rejection selftest failed")
    costly = dict(outcomes)
    costly[("p48_r3_pfx007", "control", 6, 0, "burst")] = synth_cell(
        "success", xors=Decimal("1100"))
    verdict = reduce_outcomes(costly, hard, source_records, k_lo=2, k_hi=13)
    if verdict["accepted"] or \
            verdict["candidates"]["p48_r3_pfx007"]["acceptance"]["gate_cost_within_tolerance"]:
        die("reducer cost rejection selftest failed")
    relabeled = dict(outcomes)
    hard_key_set = {cell_tuple(cell) for cell in hard}
    for key in domain_iter(2, 13):
        comparator = relabeled[("p48_r3_sfx380", "control", *key)]
        relabeled[("p48_r3_pfx007", "control", *key)] = dict(comparator)
        if key in hard_key_set:
            relabeled[("p48_r3_pfx007", "hard", *key)] = dict(comparator)
        if comparator["outcome"] != "field_shortfall":
            continue
        relabeled[("p48_r3_pfx007", "control", *key)] = synth_cell("q>H")
        if key in hard_key_set:
            relabeled[("p48_r3_pfx007", "hard", *key)] = synth_cell("q>H")
    verdict = reduce_outcomes(relabeled, hard, source_records, k_lo=2, k_hi=13)
    if verdict["accepted"] or \
            verdict["candidates"]["p48_r3_pfx007"]["acceptance"][
                "gate_actual_recovery_improved"]:
        die("reducer cause-relabel rejection selftest failed")
    relabeled_with_repair = dict(outcomes)
    panel_keys = {
        cell_tuple(record)
        for record in source_records["p48_r3"]
        if record["cause"] == "field_shortfall"
    }
    for key in panel_keys:
        relabeled_with_repair[
            ("p48_r3_pfx007", "control", *key)] = synth_cell("q>H")
        relabeled_with_repair[
            ("p48_r3_pfx007", "hard", *key)] = synth_cell("q>H")
    verdict = reduce_outcomes(
        relabeled_with_repair, hard, source_records, k_lo=2, k_hi=13)
    relabel_acceptance = verdict["candidates"]["p48_r3_pfx007"]["acceptance"]
    if verdict["accepted"] or \
            not relabel_acceptance["gate_actual_recovery_improved"] or \
            relabel_acceptance["gate_dev_fix_confirmed"]:
        die("reducer hard-relabel fix-confirmation selftest failed")
    raw_gain_erased = dict(outcomes)
    for key in ((7, 2, "repair-only"), (9, 0, "burst")):
        raw_gain_erased[
            ("p48_r3_pfx007", "control", *key)] = synth_cell("q>H")
        raw_gain_erased[
            ("p48_r3_pfx007", "hard", *key)] = synth_cell("q>H")
    verdict = reduce_outcomes(
        raw_gain_erased, hard, source_records, k_lo=2, k_hi=13)
    raw_acceptance = verdict["candidates"]["p48_r3_pfx007"]["acceptance"]
    if verdict["accepted"] or \
            not raw_acceptance["gate_actual_recovery_improved"] or \
            not raw_acceptance["gate_dev_fix_confirmed"] or \
            raw_acceptance["gate_raw_total_recovery_improved"]:
        die("reducer raw-total recovery rejection selftest failed")
    tied = dict(outcomes)
    tied[("p48_r3_pfx007", "control", 4, 2, "adversarial")] = \
        synth_cell("field_shortfall")
    verdict = reduce_outcomes(tied, hard, source_records, k_lo=2, k_hi=13)
    if verdict["accepted"] or \
            verdict["candidates"]["p48_r3_pfx007"]["acceptance"]["gate_field_shortfall_improved"]:
        die("reducer tie rejection selftest failed")
    unfaithful = dict(outcomes)
    unfaithful[("p48_r3_sfx380", "control", 3, 0, "burst")] = synth_cell("success")
    unfaithful[("p48_r3_sfx380", "hard", 3, 0, "burst")] = synth_cell("success")
    expect_reject(
        lambda: reduce_outcomes(unfaithful, hard, source_records, k_lo=2, k_hi=13),
        "comparator source replay",
    )
    diverged = dict(outcomes)
    diverged[("p48_r3_pfx007", "hard", 3, 0, "burst")] = synth_cell("field_shortfall")
    expect_reject(
        lambda: reduce_outcomes(diverged, hard, source_records, k_lo=2, k_hi=13),
        "hard-stage replay divergence",
    )
    missing = dict(outcomes)
    del missing[("p32_r7_pfx07f", "control", 2, 0, "burst")]
    expect_reject(
        lambda: reduce_outcomes(missing, hard, source_records, k_lo=2, k_hi=13),
        "missing cell",
    )

    with tempfile.TemporaryDirectory(prefix="wirehair-allk-holdout-selftest-") as temporary:
        # Write-once immutability.
        immutable = Path(temporary) / "immutable.json"
        write_once(immutable, canonical_json({"v": 1}))
        write_once(immutable, canonical_json({"v": 1}))
        expect_reject(
            lambda: write_once(immutable, canonical_json({"v": 2})),
            "write-once replacement",
        )

        # Environment-ownership marker lifecycle.
        marker_path = Path(temporary) / "owner.json"
        boot_path = Path(temporary) / "boot_id"
        boot_path.write_text("boot-1\n", encoding="ascii")
        moment = datetime(2026, 7, 21, 0, 0, 0, tzinfo=timezone.utc)
        entry = {"role": "controller", "pid": 1234, "start_ticks": 99}
        marker = write_owner_marker(
            Path("/tmp/example-root"), "executing", [entry], 48,
            marker_path=marker_path, boot_id_path=boot_path, now=moment)
        if marker_path.read_bytes() != canonical_json(marker) or \
                marker["schema"] != OWNER_MARKER_SCHEMA or \
                marker["boot_id"] != "boot-1" or \
                marker["created_utc"] != "2026-07-21T00:00:00+00:00" or \
                marker["expires_utc"] != "2026-07-23T00:00:00+00:00" or \
                marker["protected"] != [entry] or \
                marker_path.with_name(marker_path.name + ".partial").exists():
            die("owner marker publish selftest failed")
        active = load_active_owner_marker(
            marker_path=marker_path, boot_id_path=boot_path, now=moment)
        if active is None or active.get("campaign_root") != "/tmp/example-root":
            die("owner marker active-load selftest failed")
        if active_owner_live_identities(
                active, alive_probe=lambda pid: pid == 1234,
                start_ticks_probe=lambda pid: 99) != [entry] or \
                active_owner_live_identities(
                    active, alive_probe=lambda pid: pid == 1234,
                    start_ticks_probe=lambda pid: 100):
            die("owner protected-identity selftest failed")
        if load_active_owner_marker(
                marker_path=marker_path, boot_id_path=boot_path,
                now=moment + timedelta(hours=49)) is not None:
            die("owner marker expiry selftest failed")
        boot_path.write_text("boot-2\n", encoding="ascii")
        if load_active_owner_marker(
                marker_path=marker_path, boot_id_path=boot_path, now=moment) is not None:
            die("owner marker boot inertness selftest failed")
        boot_path.write_text("boot-1\n", encoding="ascii")
        write_owner_marker(
            Path("/tmp/example-root"), "complete", [], 48,
            marker_path=marker_path, boot_id_path=boot_path, now=moment)
        if load_active_owner_marker(
                marker_path=marker_path, boot_id_path=boot_path, now=moment) is not None:
            die("owner marker completion selftest failed")
        marker_path.write_bytes(canonical_json({
            "schema": OWNER_MARKER_SCHEMA,
            "campaign_root": "/tmp/example-root",
            "phase": "complete",
            "boot_id": "boot-1",
            "created_utc": "2026-07-21T00:00:00+00:00",
            "expires_utc": "2026-07-23T00:00:00+00:00",
            "protected": [{"role": "reader", "pid": True, "start_ticks": 1}],
        }))
        expect_reject(
            lambda: load_active_owner_marker(
                marker_path=marker_path, boot_id_path=boot_path, now=moment),
            "malformed complete owner marker",
        )
        marker_path.write_text("{malformed", encoding="ascii")
        expect_reject(
            lambda: load_active_owner_marker(
                marker_path=marker_path, boot_id_path=boot_path, now=moment),
            "malformed owner marker",
        )
        expect_reject(
            lambda: write_owner_marker(
                Path("/tmp/example-root"), "executing",
                [{"role": "controller", "pid": 1234}], 48,
                marker_path=marker_path, boot_id_path=boot_path, now=moment),
            "malformed protected entry",
        )

        lock_path = Path(temporary) / "owner.lock"
        first_lock = acquire_environment_lock(lock_path)
        if first_lock is None:
            die("environment lock acquisition selftest failed")
        second_lock = acquire_environment_lock(lock_path)
        if second_lock is not None:
            os.close(second_lock)
            die("environment lock exclusion selftest failed")
        os.close(first_lock)
        third_lock = acquire_environment_lock(lock_path)
        if third_lock is None:
            die("environment lock release selftest failed")
        os.close(third_lock)

        # Thermal receipt validation: acceptance, tampering, and the
        # control-stage-only CPU busy floor scoping.
        fixture = Path(temporary) / "success"
        fixture.mkdir()
        fixture_tasks, fixture_design = make_success_fixture(fixture)
        receipt = validate_thermal(fixture, fixture_design, fixture_tasks)
        if receipt["successful_segments"] != 1 or receipt["samples"] != 3:
            die("synthetic segment validation selftest failed")

        def tampered_copy(name: str) -> Path:
            target = Path(temporary) / name
            shutil.copytree(fixture, target)
            return target

        cadence = tampered_copy("cadence")
        cadence_csv = cadence / "thermal" / "segment000.csv"
        cadence_csv.write_bytes(make_thermal_bytes(("100.0", "101.0", "104.0")))
        cadence_final = load_json(segment_path(cadence, 0, "final"))
        cadence_final["thermal_csv_sha256"] = sha256_file(cadence_csv)
        segment_path(cadence, 0, "final").write_bytes(canonical_json(cadence_final))
        expect_reject(
            lambda: validate_thermal(cadence, fixture_design, fixture_tasks),
            "thermal cadence",
        )

        dimm = tampered_copy("dimm")
        dimm_csv = dimm / "thermal" / "segment000.csv"
        dimm_csv.write_bytes(make_thermal_bytes(("100.0", "101.0", "102.0"), blank_dimm=True))
        dimm_final = load_json(segment_path(dimm, 0, "final"))
        dimm_final["thermal_csv_sha256"] = sha256_file(dimm_csv)
        segment_path(dimm, 0, "final").write_bytes(canonical_json(dimm_final))
        expect_reject(
            lambda: validate_thermal(dimm, fixture_design, fixture_tasks),
            "missing DIMM",
        )

        impossible_cpu = tampered_copy("impossible-cpu")
        impossible_cpu_csv = impossible_cpu / "thermal" / "segment000.csv"
        impossible_cpu_csv.write_bytes(make_thermal_bytes(
            ("100.0", "101.0", "102.0"), cpu_temperature="9999"))
        impossible_cpu_final = load_json(
            segment_path(impossible_cpu, 0, "final"))
        impossible_cpu_final["thermal_csv_sha256"] = sha256_file(
            impossible_cpu_csv)
        segment_path(impossible_cpu, 0, "final").write_bytes(
            canonical_json(impossible_cpu_final))
        expect_reject(
            lambda: validate_thermal(
                impossible_cpu, fixture_design, fixture_tasks),
            "implausible CPU temperature",
        )

        impossible_dimm = tampered_copy("impossible-dimm")
        impossible_dimm_csv = impossible_dimm / "thermal" / "segment000.csv"
        impossible_dimm_csv.write_bytes(make_thermal_bytes(
            ("100.0", "101.0", "102.0"), dimm_temperature="-9999"))
        impossible_dimm_final = load_json(
            segment_path(impossible_dimm, 0, "final"))
        impossible_dimm_final["thermal_csv_sha256"] = sha256_file(
            impossible_dimm_csv)
        segment_path(impossible_dimm, 0, "final").write_bytes(
            canonical_json(impossible_dimm_final))
        expect_reject(
            lambda: validate_thermal(
                impossible_dimm, fixture_design, fixture_tasks),
            "implausible DIMM temperature",
        )

        edac = tampered_copy("edac")
        edac_csv = edac / "thermal" / "segment000.csv"
        edac_csv.write_bytes(make_thermal_bytes(("100.0", "101.0", "102.0"), edac_ce="1"))
        edac_final = load_json(segment_path(edac, 0, "final"))
        edac_final["thermal_csv_sha256"] = sha256_file(edac_csv)
        segment_path(edac, 0, "final").write_bytes(canonical_json(edac_final))
        expect_reject(
            lambda: validate_thermal(edac, fixture_design, fixture_tasks),
            "EDAC",
        )

        coverage = tampered_copy("coverage")
        coverage_ready = load_json(segment_path(coverage, 0, "ready"))
        coverage_ready["jobs_started_monotonic_s"] = 99.0
        segment_path(coverage, 0, "ready").write_bytes(canonical_json(coverage_ready))
        coverage_final = load_json(segment_path(coverage, 0, "final"))
        coverage_final["ready_sha256"] = sha256_file(segment_path(coverage, 0, "ready"))
        segment_path(coverage, 0, "final").write_bytes(canonical_json(coverage_final))
        expect_reject(
            lambda: validate_thermal(coverage, fixture_design, fixture_tasks),
            "thermal start coverage",
        )

        sampler_exit = tampered_copy("sampler-exit")
        sampler_final = load_json(segment_path(sampler_exit, 0, "final"))
        sampler_final["sampler_returncode"] = 9
        segment_path(sampler_exit, 0, "final").write_bytes(canonical_json(sampler_final))
        expect_reject(
            lambda: validate_thermal(sampler_exit, fixture_design, fixture_tasks),
            "sampler exit",
        )

        hard_burst = Path(temporary) / "hard-burst"
        hard_burst.mkdir()
        hard_tasks, hard_design = make_success_fixture(hard_burst, busy="20.0")
        hard_receipt = validate_thermal(hard_burst, hard_design, hard_tasks)
        if hard_receipt["successful_segments"] != 1:
            die("hard burst segment below the control floor was not accepted")

        control_floor = Path(temporary) / "control-floor"
        control_floor.mkdir()
        control_tasks, control_design = make_success_fixture(
            control_floor, stage="control", busy=str(CPU_BUSY_FLOOR))
        control_receipt = validate_thermal(
            control_floor, control_design, control_tasks)
        if control_receipt["successful_segments"] != 1:
            die("control segment at the sealed floor was not accepted")

        control_low = Path(temporary) / "control-low"
        control_low.mkdir()
        low_tasks, low_design = make_success_fixture(
            control_low, stage="control", busy="89.9")
        expect_reject(
            lambda: validate_thermal(control_low, low_design, low_tasks),
            "control busy floor",
        )

        interrupted = Path(temporary) / "interrupted"
        make_private_directory(interrupted)
        for directory in RUNTIME_DIRECTORY_NAMES:
            make_private_directory(interrupted / directory)
        retry_task = {"job": 0, "stage": "hard", "output_name": "retry.csv", "argv": ["/bin/true"]}
        retry_intent = {
            "schema": SCHEMA + ".segment_intent", "segment": 0,
            "boot_id": required_boot_id(),
            "controller_pid": 99999998, "controller_start_ticks": 1,
            "stage": "hard", "jobs": [0],
            "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
            "retry_policy": "stage-atomic non-selective retry",
        }
        write_once(segment_path(interrupted, 0, "intent"), canonical_json(retry_intent))
        make_private_directory(interrupted / "attempts" / "segment000")
        (interrupted / "raw" / "retry.csv").write_bytes(b"failed output\n")
        (interrupted / "stderr" / "retry.csv.stderr").write_bytes(b"failure\n")
        (interrupted / "exit" / "retry.csv.exit").write_bytes(b"1\n")
        reconciled = reconcile_incomplete_segments(interrupted, [retry_task])
        if len(reconciled) != 1 or reconciled[0]["state"] != "interrupted" or \
                completed_jobs(interrupted, [retry_task]) or \
                any((interrupted / directory / name).exists() for directory, name in (
                    ("raw", "retry.csv"), ("stderr", "retry.csv.stderr"),
                    ("exit", "retry.csv.exit"))) or \
                reconcile_incomplete_segments(interrupted, [retry_task]):
            die("interrupted failed-triplet resume selftest failed")
        expect_reject(lambda: validate_thermal_rows([], segment=0), "header-only thermal")
    print("wh2_row_mask_allk_holdout selftest: PASS")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="freeze a new holdout; never launches jobs")
    prepare.add_argument("--output-root", required=True)
    prepare.add_argument("--allk-root", required=True)
    prepare.add_argument("--sampler", required=True)
    prepare.add_argument(
        "--workers", type=int, default=DEFAULT_WORKERS,
        choices=range(MIN_CAMPAIGN_WORKERS, MAX_CAMPAIGN_WORKERS + 1))
    prepare.add_argument("--owner-ttl-hours", type=int, default=DEFAULT_OWNER_TTL_HOURS)
    prepare.set_defaults(function=command_prepare)
    launch = sub.add_parser("launch", help="preflight or execute a frozen holdout")
    launch.add_argument("--root", required=True)
    launch.add_argument("--expected-receipts-sha256", required=True)
    mode = launch.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight-only", action="store_true")
    mode.add_argument("--execute", action="store_true")
    launch.add_argument("--resume", action="store_true")
    launch.add_argument("--stage", choices=("hard", "control", "all"), default="all")
    launch.set_defaults(function=command_launch)
    reduce_parser = sub.add_parser("reduce", help="strictly validate and reduce all outputs")
    reduce_parser.add_argument("--root", required=True)
    reduce_parser.add_argument("--expected-receipts-sha256", required=True)
    reduce_parser.set_defaults(function=command_reduce)
    selftest = sub.add_parser("selftest", help="run bounded pure-Python invariants")
    selftest.set_defaults(function=command_selftest)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "resume", False) and not getattr(args, "execute", False):
        parser.error("--resume requires --execute")
    try:
        args.function(args)
    except CampaignError as exc:
        print("campaign error: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
