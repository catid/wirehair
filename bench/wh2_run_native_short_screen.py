#!/usr/bin/env python3
"""Run the fixed, native WH2 development short screen.

This is deliberately a narrow controller for the frozen development domain.
It emits and freezes arm-free traces before starting any result job, keeps one
persistent native worker on each selected physical core, executes timing one
homogeneous eight-job wave and barrier at a time, and delegates publication to
the fail-closed native evidence assembler.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import hashlib
import io
import math
import os
from pathlib import Path
import platform
import re
import selectors
import signal
import socket
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional
from typing import Sequence, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api


DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v1"
DESCRIPTION_FIELDS = frozenset((
    "schema", "source_git_commit", "binary_sha256", "arms",
))
DESCRIPTION_ARM_FIELDS = frozenset((
    "arm", "codec", "arm_descriptor_sha256",
))
EXPECTED_ARMS = (
    ("wirehair2_head", "wirehair2_certified"),
    ("wirehair1", "wirehair1"),
    ("wirehair2_identity", "wirehair2_experiment"),
)
SUMMARY_SCHEMA = "wirehair.wh2.native-short-screen-run.v1"
MAX_WALL_SECONDS = 7200.0
MAX_RESPONSE_BYTES = 1024 * 1024
TIMING_WORKER_COUNT = 8
ZERO_SHA256 = "0" * 64
CPU_TOKEN = re.compile(r"(?:0|[1-9][0-9]*)\Z")
CACHE_INDEX_TOKEN = re.compile(r"index(?:0|[1-9][0-9]*)\Z")
POSITIVE_TOKEN = re.compile(r"[1-9][0-9]*\Z")


class RunnerError(RuntimeError):
    """The real native short screen cannot continue safely."""


def fail(message: str) -> None:
    raise RunnerError(message)


def _remaining(deadline: float, context: str) -> float:
    value = deadline - time.monotonic()
    if value <= 0.0:
        fail("native short-screen hard wall expired during {}".format(context))
    return value


@contextmanager
def _hard_wall(seconds: float):
    """Interrupt synchronous validation too, not just subprocess/select waits."""
    if not hasattr(signal, "setitimer") or not hasattr(signal, "ITIMER_REAL"):
        fail("native short-screen hard wall requires POSIX setitimer")
    started = time.monotonic()

    def expired(_signum: int, _frame: Any) -> None:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        fail("{}-second native short-screen hard wall expired".format(
            int(seconds)))

    try:
        old_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, expired)
    except (OSError, ValueError) as exc:
        fail("cannot install native short-screen hard wall: {}".format(exc))
    try:
        old_timer = signal.setitimer(signal.ITIMER_REAL, seconds)
    except (OSError, ValueError) as exc:
        signal.signal(signal.SIGALRM, old_handler)
        fail("cannot install native short-screen hard wall: {}".format(exc))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)
        if old_timer[0] > 0.0:
            remaining = max(0.000001,
                            old_timer[0] - (time.monotonic() - started))
            signal.setitimer(signal.ITIMER_REAL, remaining, old_timer[1])


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and contract_api.SHA256.fullmatch(value) \
        is not None


def _hash_file(path: Path) -> str:
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


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | \
        getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(str(path), flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_staged(staged: Path, destination: Path) -> None:
    if destination.exists() or destination.is_symlink():
        fail("refusing to replace existing artifact {}".format(destination))
    try:
        staged_info = os.stat(str(staged), follow_symlinks=False)
    except OSError as exc:
        fail("cannot inspect staged artifact {}: {}".format(staged, exc))
    linked = False
    try:
        os.link(str(staged), str(destination))
        linked = True
        _fsync_directory(destination.parent)
        staged.unlink()
    except FileExistsError:
        fail("refusing to replace existing artifact {}".format(destination))
    except OSError as exc:
        if linked:
            _unlink_if_identity(
                destination, staged_info.st_dev, staged_info.st_ino)
        fail("cannot publish {}: {}".format(destination, exc))
    except BaseException:
        # The hard-wall signal can arrive after link(2) but before the
        # directory fsync.  Do not leave a result that was never durably
        # published, and never unlink a path another process replaced.
        if linked:
            _unlink_if_identity(
                destination, staged_info.st_dev, staged_info.st_ino)
        raise


def _unlink_if_identity(path: Path, device: int, inode: int) -> bool:
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
    except OSError:
        return False


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    descriptor, raw = tempfile.mkstemp(
        prefix="." + path.name + ".", suffix=".tmp", dir=str(path.parent))
    staged = Path(raw)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        _publish_staged(staged, path)
    finally:
        try:
            staged.unlink()
        except FileNotFoundError:
            pass


def _atomic_write_object(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_write_bytes(
        path, (contract_api.canonical_json(value) + "\n").encode("utf-8"))


class AtomicLineSink:
    """A result stream that is visible only after every row is complete."""

    def __init__(self, destination: Path) -> None:
        self.destination = destination
        descriptor, raw = tempfile.mkstemp(
            prefix="." + destination.name + ".", suffix=".tmp",
            dir=str(destination.parent))
        self.staged = Path(raw)
        self.output = os.fdopen(descriptor, "wb")
        self.published = False

    def write(self, line: bytes) -> None:
        if self.output.closed or not line.endswith(b"\n"):
            fail("cannot write an invalid native result line")
        self.output.write(line)

    def publish(self) -> None:
        if self.published:
            fail("native result stream was published twice")
        self.output.flush()
        os.fsync(self.output.fileno())
        self.output.close()
        _publish_staged(self.staged, self.destination)
        self.published = True

    def abort(self) -> None:
        if not self.output.closed:
            self.output.close()
        try:
            self.staged.unlink()
        except FileNotFoundError:
            pass


def _parse_canonical_line(data: bytes, context: str) -> Mapping[str, Any]:
    if not data.endswith(b"\n") or data.count(b"\n") != 1:
        fail("{} must be exactly one newline-terminated JSON object".format(
            context))
    logical = data[:-1]
    try:
        value = contract_api._load_json_bytes(logical, context)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if not isinstance(value, dict) or \
            contract_api.canonical_json(value).encode("utf-8") != logical:
        fail("{} is not canonical JSON".format(context))
    return value


def parse_cpu_list(text: str) -> List[int]:
    """Parse a comma-separated, strictly increasing canonical CPU list."""
    if not isinstance(text, str) or not text:
        fail("--cpus must be a nonempty comma-separated CPU list")
    tokens = text.split(",")
    if any(CPU_TOKEN.fullmatch(token) is None for token in tokens):
        fail("--cpus contains a noncanonical CPU ID")
    values = [int(token) for token in tokens]
    if values != sorted(set(values)):
        fail("--cpus must be strictly increasing and unique")
    return values


def _cpu_topology(cpu: int, sysfs_root: Path) -> Tuple[int, int]:
    topology = sysfs_root / "cpu{}".format(cpu) / "topology"
    values = []
    for name in ("physical_package_id", "core_id"):
        try:
            text = (topology / name).read_text(encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} topology: {}".format(cpu, exc))
        if CPU_TOKEN.fullmatch(text) is None:
            fail("CPU {} has a noncanonical {}".format(cpu, name))
        values.append(int(text))
    return values[0], values[1]


def _cpu_llc_identity(cpu: int, sysfs_root: Path) -> Tuple[int, int]:
    """Read the highest data/unified cache level and ID without assuming index3."""
    cache_root = sysfs_root / "cpu{}".format(cpu) / "cache"
    try:
        entries = sorted(
            (path for path in cache_root.iterdir()
             if CACHE_INDEX_TOKEN.fullmatch(path.name) is not None),
            key=lambda path: int(path.name[5:]))
    except OSError as exc:
        fail("cannot enumerate CPU {} cache topology: {}".format(cpu, exc))
    if not entries:
        fail("CPU {} has no indexed cache topology".format(cpu))
    data_caches: List[Tuple[int, int]] = []
    for entry in entries:
        try:
            level_text = (entry / "level").read_text(
                encoding="ascii").strip()
            cache_type = (entry / "type").read_text(
                encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} cache topology: {}".format(cpu, exc))
        if POSITIVE_TOKEN.fullmatch(level_text) is None or \
                cache_type not in ("Data", "Instruction", "Unified"):
            fail("CPU {} has malformed cache level/type metadata".format(cpu))
        if cache_type == "Instruction":
            continue
        try:
            cache_id_text = (entry / "id").read_text(
                encoding="ascii").strip()
        except OSError as exc:
            fail("cannot read CPU {} cache ID: {}".format(cpu, exc))
        if CPU_TOKEN.fullmatch(cache_id_text) is None:
            fail("CPU {} has a noncanonical cache ID".format(cpu))
        data_caches.append((int(level_text), int(cache_id_text)))
    if not data_caches:
        fail("CPU {} has no data/unified cache topology".format(cpu))
    highest_level = max(level for level, _ in data_caches)
    ids = {cache_id for level, cache_id in data_caches
           if level == highest_level}
    if len(ids) != 1:
        fail("CPU {} has ambiguous last-level cache IDs".format(cpu))
    return highest_level, next(iter(ids))


def _select_llc_spread_workers(
        candidates: Mapping[Tuple[int, int], int],
        llc_by_core: Mapping[
            Tuple[int, int], Tuple[int, int, int]]) -> List[int]:
    """Pick eight cores with maximum LLC coverage, then round-robin fill."""
    grouped: Dict[Tuple[int, int, int], List[int]] = {}
    for core, cpu in candidates.items():
        grouped.setdefault(llc_by_core[core], []).append(cpu)
    for values in grouped.values():
        values.sort()
    selected: List[int] = []
    group_ids = sorted(grouped)
    while len(selected) < TIMING_WORKER_COUNT:
        progressed = False
        for cache_id in group_ids:
            values = grouped[cache_id]
            if values:
                selected.append(values.pop(0))
                progressed = True
                if len(selected) == TIMING_WORKER_COUNT:
                    break
        if not progressed:
            break
    if len(selected) != TIMING_WORKER_COUNT:
        fail("fewer than eight eligible physical worker cores are available")
    return sorted(selected)


def select_cpu_layout(
        sampler_cpu: int, explicit_workers: Optional[Sequence[int]] = None,
        explicit_controller: Optional[int] = None,
        affinity: Optional[Iterable[int]] = None,
        sysfs_root: Path = Path("/sys/devices/system/cpu"),
        ) -> Tuple[List[int], int]:
    """Select eight LLC-spread physical cores plus an isolated controller."""
    if type(sampler_cpu) is not int or sampler_cpu < 0:
        fail("sampler CPU must be a nonnegative integer")
    if affinity is None:
        try:
            allowed = sorted(os.sched_getaffinity(0))
        except (AttributeError, OSError) as exc:
            fail("cannot inspect controller CPU affinity: {}".format(exc))
    else:
        affinity_values = list(affinity)
        if any(type(cpu) is not int or cpu < 0 for cpu in affinity_values):
            fail("controller CPU affinity contains an invalid CPU ID")
        allowed = sorted(set(affinity_values))
    if not allowed:
        fail("controller CPU affinity is empty")
    sampler_core = _cpu_topology(sampler_cpu, sysfs_root)
    topology = {cpu: _cpu_topology(cpu, sysfs_root) for cpu in allowed}
    # Cache IDs need not be globally unique across packages or cache levels.
    # Bind both dimensions before comparing last-level cache groups.
    llc_by_cpu = {
        cpu: (topology[cpu][0],) + _cpu_llc_identity(cpu, sysfs_root)
        for cpu in allowed
    }
    allowed_set = set(allowed)

    by_core: Dict[Tuple[int, int], int] = {}
    llc_by_core: Dict[Tuple[int, int], Tuple[int, int, int]] = {}
    for cpu in allowed:
        core = topology[cpu]
        cache_id = llc_by_cpu[cpu]
        if core in llc_by_core and llc_by_core[core] != cache_id:
            fail("SMT siblings disagree on their last-level cache ID")
        llc_by_core[core] = cache_id
        by_core[core] = min(cpu, by_core.get(core, cpu))

    controller_core: Optional[Tuple[int, int]] = None
    if explicit_controller is not None:
        if type(explicit_controller) is not int or explicit_controller < 0 or \
                explicit_controller not in allowed_set:
            fail("controller CPU is outside current sched affinity")
        controller_core = topology[explicit_controller]
        if explicit_controller == sampler_cpu or controller_core == sampler_core:
            fail("controller CPU overlaps the sampler's physical core")

    worker_candidates = {
        core: cpu for core, cpu in by_core.items()
        if core != sampler_core and core != controller_core
    }

    if explicit_workers is not None:
        workers = list(explicit_workers)
        if (workers != sorted(set(workers)) or
                len(workers) != TIMING_WORKER_COUNT):
            fail("explicit worker CPUs must be eight sorted unique IDs")
        if any(type(cpu) is not int or cpu < 0 for cpu in workers):
            fail("explicit worker CPUs contain an invalid CPU ID")
        if not set(workers).issubset(allowed_set):
            fail("explicit worker CPUs are outside current sched affinity")
        worker_cores = [topology[cpu] for cpu in workers]
        if len(worker_cores) != len(set(worker_cores)):
            fail("explicit worker CPUs contain SMT siblings")
        if sampler_cpu in workers or sampler_core in worker_cores:
            fail("worker CPUs overlap the sampler's physical core")
        if controller_core in worker_cores:
            fail("controller CPU shares a physical core with a worker")
    else:
        workers = _select_llc_spread_workers(
            worker_candidates, llc_by_core)
        worker_cores = [topology[cpu] for cpu in workers]

    available_llc_ids = {
        llc_by_core[core] for core in worker_candidates
    }
    required_llc_coverage = min(
        TIMING_WORKER_COUNT, len(available_llc_ids))
    observed_llc_coverage = len({
        llc_by_core[topology[cpu]] for cpu in workers
    })
    if observed_llc_coverage != required_llc_coverage:
        fail("worker CPU roster does not maximize available LLC coverage")

    if explicit_controller is not None:
        controller_cpu = explicit_controller
    else:
        worker_core_set = set(worker_cores)
        controller_candidates = {
            core: cpu for core, cpu in by_core.items()
            if core != sampler_core and core not in worker_core_set
        }
        if controller_candidates:
            controller_cpu = max(controller_candidates.values())
            controller_core = topology[controller_cpu]
        else:
            fail("no dedicated controller core remains outside frozen workers")

    if (len(workers) != TIMING_WORKER_COUNT or
            controller_cpu in workers or controller_core == sampler_core or
            controller_core in {topology[cpu] for cpu in workers}):
        fail("controller/worker/sampler CPU layout is not physically isolated")
    return workers, controller_cpu


def _pin_controller(cpu: int) -> None:
    try:
        os.sched_setaffinity(0, {cpu})
        observed = os.sched_getaffinity(0)
    except (AttributeError, OSError) as exc:
        fail("cannot pin controller to CPU {}: {}".format(cpu, exc))
    if observed != {cpu}:
        fail("controller did not establish singleton CPU affinity")


def _restore_controller_affinity(affinity: Set[int]) -> None:
    try:
        os.sched_setaffinity(0, affinity)
        observed = os.sched_getaffinity(0)
    except (AttributeError, OSError) as exc:
        fail("cannot restore controller CPU affinity: {}".format(exc))
    if observed != affinity:
        fail("controller CPU affinity did not restore exactly")


def _run_command(argv: Sequence[str], deadline: float,
                 context: str) -> bytes:
    try:
        result = subprocess.run(
            list(argv), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False,
            timeout=_remaining(deadline, context))
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("{} failed: {}".format(context, exc))
    if result.returncode != 0:
        diagnostic = result.stderr.decode("utf-8", "replace").strip()
        fail("{} exited {}{}".format(
            context, result.returncode,
            ": " + diagnostic if diagnostic else ""))
    if result.stderr:
        fail("{} produced unexpected stderr: {}".format(
            context, result.stderr.decode("utf-8", "replace").strip()))
    return result.stdout


def describe_worker(worker_path: Path, deadline: float) -> Mapping[str, Any]:
    try:
        resolved = worker_path.resolve(strict=True)
        info = resolved.stat()
    except OSError as exc:
        fail("cannot resolve native worker {}: {}".format(worker_path, exc))
    if not stat.S_ISREG(info.st_mode) or not os.access(str(resolved), os.X_OK):
        fail("native worker must be an executable regular file")
    binary_sha256 = _hash_file(resolved)
    value = _parse_canonical_line(
        _run_command([str(resolved), "--describe"], deadline,
                     "native worker description"),
        "native worker description")
    if set(value) != DESCRIPTION_FIELDS or \
            value.get("schema") != DESCRIPTION_SCHEMA or \
            not isinstance(value.get("source_git_commit"), str) or \
            re.fullmatch(r"[0-9a-f]{40}",
                         value["source_git_commit"]) is None or \
            value.get("binary_sha256") != binary_sha256:
        fail("native worker description does not bind its executable")
    arms = value.get("arms")
    if not isinstance(arms, list) or len(arms) != len(EXPECTED_ARMS):
        fail("native worker description has an invalid arm roster")
    for index, expected in enumerate(EXPECTED_ARMS):
        arm = arms[index]
        if not isinstance(arm, dict) or set(arm) != DESCRIPTION_ARM_FIELDS or \
                arm.get("arm") != expected[0] or \
                arm.get("codec") != expected[1] or \
                not _is_sha256(arm.get("arm_descriptor_sha256")):
            fail("native worker description has an invalid arm at index {}".format(
                index))
    if len({arm["arm_descriptor_sha256"] for arm in arms}) != len(arms):
        fail("native worker description reuses an arm descriptor")
    result = dict(value)
    result["resolved_path"] = str(resolved)
    return result


def _require_worker_source_commit(
        description: Mapping[str, Any], source_commit: str) -> None:
    if description.get("source_git_commit") != source_commit:
        fail("native worker was not built from the clean source HEAD")


def _git_head(deadline: float, repo_root: Optional[Path] = None) -> str:
    """Return HEAD only for a clean, stable codec-source snapshot."""
    root = Path(__file__).resolve().parent.parent if repo_root is None \
        else repo_root
    commands = (
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
        # Any tracked change can invalidate the source attribution (including
        # a build/configuration input outside the obvious codec directories).
        ["git", "status", "--porcelain=v1", "--untracked-files=no"],
        # Untracked experiment/result artifacts are common and harmless.  An
        # untracked build source in one of these locations is not.
        ["git", "ls-files", "--others", "--exclude-standard", "--",
         "CMakeLists.txt", "codec", "bench", "include", "cmake",
         ":(top,glob)*.c", ":(top,glob)*.cc", ":(top,glob)*.cpp",
         ":(top,glob)*.h", ":(top,glob)*.hpp", ":(top,glob)*.inc"],
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
    )
    outputs = []
    try:
        for command in commands:
            result = subprocess.run(
                command, cwd=str(root), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                timeout=_remaining(deadline, "checking source commit"))
            if result.returncode != 0:
                diagnostic = result.stderr.decode("utf-8", "replace").strip()
                fail("cannot establish clean source commit{}".format(
                    ": " + diagnostic if diagnostic else ""))
            outputs.append(result.stdout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot establish clean source commit: {}".format(exc))
    try:
        first = outputs[0].decode("ascii").strip()
        status = outputs[1].decode("utf-8")
        untracked = outputs[2].decode("utf-8")
        second = outputs[3].decode("ascii").strip()
    except UnicodeError as exc:
        fail("git source identity is not valid text: {}".format(exc))
    if re.fullmatch(r"[0-9a-f]{40}", first) is None or first != second:
        fail("git HEAD is not a full lowercase 40-hex commit")
    if status or untracked:
        paths = [line for line in (status + untracked).splitlines() if line]
        fail("codec source tree is dirty or has untracked source: {}".format(
            "; ".join(paths[:8])))
    return first


def _host_identity(controller_cpu: int) -> Mapping[str, Any]:
    return {
        "controller_cpu": controller_cpu,
        "hostname": socket.gethostname(),
        "kernel_release": platform.release(),
        "machine": platform.machine(),
        "system": platform.system(),
    }


def _emit_and_assemble_trace(
        contract: Mapping[str, Any], kind: str,
        description: Mapping[str, Any], output_dir: Path,
        deadline: float) -> Tuple[Path, str]:
    worker_path = description["resolved_path"]
    native_path = output_dir / "{}-native-traces.jsonl".format(kind)
    trace_path = output_dir / "{}-traces.jsonl".format(kind)
    data = _run_command(
        [worker_path, "--emit-traces", kind], deadline,
        "{} native trace emission".format(kind))
    _atomic_write_bytes(native_path, data)
    try:
        digest = native_api.assemble_trace_manifest(
            contract, kind, "development", native_path, trace_path)
    except (native_api.NativeEvidenceError,
            contract_api.ContractError) as exc:
        fail(str(exc))
    return trace_path, digest


def _freeze_manifest(
        contract: Mapping[str, Any], kind: str, trace_sha256: str,
        source_commit: str, description: Mapping[str, Any],
        cpus: Sequence[int], controller_cpu: int) -> Mapping[str, Any]:
    roster = [value[0] for value in EXPECTED_ARMS]
    arms = []
    for value in description["arms"]:
        arms.append({
            "arm": value["arm"],
            "codec": value["codec"],
            "binary_sha256": description["binary_sha256"],
            "arm_descriptor_sha256": value["arm_descriptor_sha256"],
            "construction_policy": "not_applicable"
                if value["codec"] == "wirehair1" else "raw_base",
            "repair_map_sha256": ZERO_SHA256,
        })
    commands = [
        [description["resolved_path"], "--emit-traces", kind],
    ] + [
        [description["resolved_path"], "--worker", str(cpu)] for cpu in cpus
    ]
    return {
        "schema": contract_api.FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": kind,
        "phase": "development",
        "domain_sha256": contract[kind]["domains"]["development"][
            "domain_sha256"],
        "source_git_commit": source_commit,
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": commands,
        "cpu_affinity": list(cpus),
        "host_identity": _host_identity(controller_cpu),
        "arms": arms,
    }


def write_development_freezes(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        cpus: Sequence[int], controller_cpu: int, source_commit: str,
        trace_sha256s: Mapping[str, str], output_dir: Path,
        ) -> Mapping[str, Mapping[str, Any]]:
    """Publish and re-open both freezes before a result process may start."""
    loaded: Dict[str, Mapping[str, Any]] = {}
    for kind in ("recovery", "timing"):
        path = output_dir / "{}-freeze.json".format(kind)
        value = _freeze_manifest(
            contract, kind, trace_sha256s[kind], source_commit,
            description, cpus, controller_cpu)
        _atomic_write_object(path, value)
        try:
            loaded[kind] = contract_api.load_freeze_manifest(
                contract, "development", path, kind)
        except contract_api.ContractError as exc:
            fail(str(exc))
    return loaded


def _valid_sampler_monotonic_ns(row: Sequence[str]) -> Optional[int]:
    if len(row) != len(native_api.THERMAL_HEADER):
        return None
    try:
        numeric = [float(row[index]) for index in (
            1, 2, 3, 4, 14, 15, 16)]
        dimms = [float(value) for value in row[5:13]]
    except (TypeError, ValueError):
        return None
    if any(not math.isfinite(value) for value in numeric + dimms):
        return None
    canonical_counters = []
    for index in (13, 17, 18):
        value = row[index]
        if not value or not value.isascii() or not value.isdecimal() or \
                (len(value) > 1 and value.startswith("0")):
            return None
        canonical_counters.append(int(value))
    mono, busy, mhz, cpu_c, load1, load5, load15 = numeric
    if (mono < 0.0 or not 0.0 <= busy <= 100.0 or mhz <= 0.0 or
            not 0.0 <= cpu_c <= 120.0 or
            any(not 0.0 < value <= 100.0 for value in dimms) or
            any(value < 0.0 for value in (load1, load5, load15)) or
            any(value != 0 for value in canonical_counters)):
        return None
    return int(round(mono * 1000000000.0))


def _valid_sampler_samples(csv_path: Path) -> List[int]:
    values = []
    try:
        data = csv_path.read_bytes()
        # The live sampler may be in the middle of appending its last row.
        # Never select bytes which the sealed attestation correctly rejects as
        # an unterminated endpoint.
        complete_end = data.rfind(b"\n") + 1
        complete = data[:complete_end]
        text = complete.decode("ascii")
        with io.StringIO(text, newline="") as source:
            reader = csv.reader(source)
            if next(reader, None) != list(native_api.THERMAL_HEADER):
                fail("sampler CSV header differs from the frozen format")
            for row in reader:
                value = _valid_sampler_monotonic_ns(row)
                if value is not None:
                    values.append(value)
    except (OSError, UnicodeError, csv.Error) as exc:
        fail("cannot read sampler CSV: {}".format(exc))
    return values


def _wait_for_sampler_sample(
        csv_path: Path, deadline: float, *, after_ns: Optional[int] = None,
        at_or_after_ns: Optional[int] = None) -> int:
    if (after_ns is None) == (at_or_after_ns is None):
        fail("sampler wait requires exactly one time bound")
    while True:
        samples = _valid_sampler_samples(csv_path)
        if after_ns is not None:
            matches = [value for value in samples if value > after_ns]
        else:
            matches = [value for value in samples
                       if value >= int(at_or_after_ns)]
        if matches:
            return min(matches)
        remaining = _remaining(deadline, "waiting for thermal sample")
        time.sleep(min(0.2, remaining))


def choose_new_sampler_start(csv_path: Path, deadline: float) -> int:
    samples = _valid_sampler_samples(csv_path)
    if not samples:
        fail("sampler CSV has no valid baseline sample")
    return _wait_for_sampler_sample(
        csv_path, deadline, after_ns=max(samples))


def _preflight_sampler(pid: int, cpu: int, script_path: Path,
                       csv_path: Path) -> int:
    descriptors = []
    try:
        for path, context in ((script_path, "sampler script"),
                              (csv_path, "sampler CSV")):
            descriptor, _ = native_api._open_regular_nofollow(path, context)
            descriptors.append(descriptor)
        start_ticks = native_api._parse_proc_start_ticks(pid)
        native_api._verify_live_sampler_process(
            pid, cpu, start_ticks, script_path, csv_path)
        return start_ticks
    except native_api.NativeEvidenceError as exc:
        root_owned = False
        try:
            root_owned = Path("/proc/{}".format(pid)).stat().st_uid == 0
        except OSError:
            pass
        if root_owned and getattr(os, "geteuid", lambda: 1)() != 0:
            fail("{}; sampler PID {} is root-owned and its executable "
                 "provenance requires privilege. Rerun this controller with "
                 "sudo -n python3.12; validation is not weakened".format(
                     exc, pid))
        fail(str(exc))
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
    return 0


class Job:
    def __init__(self, kind: str, cell: int, item: int, ordinal: int) -> None:
        if (kind not in ("recovery", "timing") or
                any(type(value) is not int or value < 0
                    for value in (cell, item, ordinal))):
            fail("cannot construct an invalid native job")
        self.kind = kind
        self.cell = cell
        self.item = item
        self.ordinal = ordinal

    def command(self) -> bytes:
        prefix = "R" if self.kind == "recovery" else "T"
        return "{} {} {}\n".format(prefix, self.cell, self.item).encode("ascii")


class PersistentWorker:
    def __init__(self, cpu: int, process: subprocess.Popen,
                 start_ticks: int) -> None:
        self.cpu = cpu
        self.process = process
        self.start_ticks = start_ticks
        self.pending: Optional[Job] = None
        self.buffer = b""


def _worker_stderr(worker: PersistentWorker) -> str:
    descriptor = worker.process.stderr.fileno()
    try:
        data = os.read(descriptor, 65536)
    except (BlockingIOError, OSError):
        data = b""
    return data.decode("utf-8", "replace").strip()


def spawn_workers(description: Mapping[str, Any], cpus: Sequence[int],
                  deadline: float) -> List[PersistentWorker]:
    workers: List[PersistentWorker] = []
    try:
        for cpu in cpus:
            try:
                process = subprocess.Popen(
                    [description["resolved_path"], "--worker", str(cpu)],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, bufsize=0)
            except OSError as exc:
                fail("cannot spawn native worker on CPU {}: {}".format(cpu, exc))
            workers.append(PersistentWorker(cpu, process, 0))
        pending = set(range(len(workers)))
        ready_deadline = min(deadline, time.monotonic() + 10.0)
        while pending:
            if time.monotonic() >= ready_deadline:
                fail("native workers did not establish singleton affinity")
            for index in list(pending):
                worker = workers[index]
                if worker.process.poll() is not None:
                    fail("native worker on CPU {} exited during startup: {}".format(
                        worker.cpu, _worker_stderr(worker)))
                try:
                    ticks = native_api._parse_proc_start_ticks(
                        worker.process.pid)
                    affinity = os.sched_getaffinity(worker.process.pid)
                except (native_api.NativeEvidenceError,
                        AttributeError, OSError):
                    continue
                if affinity == {worker.cpu}:
                    worker.start_ticks = ticks
                    os.set_blocking(worker.process.stdout.fileno(), False)
                    os.set_blocking(worker.process.stderr.fileno(), False)
                    pending.remove(index)
            if pending:
                time.sleep(0.01)
        return workers
    except BaseException:
        terminate_workers(workers)
        raise


def terminate_workers(workers: Sequence[PersistentWorker]) -> None:
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.terminate()
            except OSError:
                pass
    limit = time.monotonic() + 2.0
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.wait(timeout=max(0.0, limit - time.monotonic()))
            except subprocess.TimeoutExpired:
                try:
                    worker.process.kill()
                except OSError:
                    pass
    for worker in workers:
        if worker.process.poll() is None:
            try:
                worker.process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                pass
        for stream in (worker.process.stdin, worker.process.stdout,
                       worker.process.stderr):
            if stream is not None:
                try:
                    stream.close()
                except OSError:
                    pass


def _write_worker(worker: PersistentWorker, data: bytes) -> None:
    if worker.process.poll() is not None:
        fail("native worker on CPU {} exited before its command".format(
            worker.cpu))
    try:
        written = os.write(worker.process.stdin.fileno(), data)
    except (BrokenPipeError, OSError) as exc:
        fail("cannot command native worker on CPU {}: {}".format(
            worker.cpu, exc))
    if written != len(data):
        fail("short write to native worker on CPU {}".format(worker.cpu))


ResponseValidator = Callable[
    [Mapping[str, Any], bytes, PersistentWorker, Job], int]


def run_job_batch(
        workers: Sequence[PersistentWorker], jobs: Sequence[Job],
        rotation: int, sink: AtomicLineSink, deadline: float,
        validator: ResponseValidator) -> Tuple[int, Set[int]]:
    """Run one barrier-delimited batch with at most one job per worker."""
    if not workers or len(jobs) < len(workers):
        fail("each batch must have enough jobs to use every frozen worker")
    if any(worker.pending is not None or worker.buffer for worker in workers):
        fail("worker state leaked across a job barrier")
    selector = selectors.DefaultSelector()
    for worker in workers:
        selector.register(worker.process.stdout, selectors.EVENT_READ, worker)
    next_job = 0
    completed = 0
    maximum_end = 0
    used: Set[int] = set()
    try:
        for offset in range(len(workers)):
            worker = workers[(offset + rotation) % len(workers)]
            worker.pending = jobs[next_job]
            next_job += 1
            _write_worker(worker, worker.pending.command())
        while completed < len(jobs):
            timeout = min(1.0, _remaining(deadline, "native result jobs"))
            events = selector.select(timeout)
            if not events:
                for worker in workers:
                    if worker.process.poll() is not None:
                        fail("native worker on CPU {} exited: {}".format(
                            worker.cpu, _worker_stderr(worker)))
                continue
            for key, _ in events:
                worker = key.data
                try:
                    block = os.read(worker.process.stdout.fileno(), 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    fail("cannot read native worker on CPU {}: {}".format(
                        worker.cpu, exc))
                if not block:
                    fail("native worker on CPU {} reached EOF: {}".format(
                        worker.cpu, _worker_stderr(worker)))
                worker.buffer += block
                if len(worker.buffer) > MAX_RESPONSE_BYTES:
                    fail("native worker response exceeds the bounded line size")
                if b"\n" not in worker.buffer:
                    continue
                logical, remainder = worker.buffer.split(b"\n", 1)
                if remainder:
                    fail("native worker emitted more than one response per command")
                line = logical + b"\n"
                worker.buffer = b""
                job = worker.pending
                if job is None:
                    fail("native worker responded without an outstanding job")
                value = _parse_canonical_line(line, "native worker response")
                finished_ns = validator(value, line, worker, job)
                if type(finished_ns) is not int or finished_ns < 0:
                    fail("native response validator returned an invalid end time")
                maximum_end = max(maximum_end, finished_ns)
                sink.write(line)
                used.add(worker.cpu)
                completed += 1
                worker.pending = None
                if next_job < len(jobs):
                    worker.pending = jobs[next_job]
                    next_job += 1
                    _write_worker(worker, worker.pending.command())
        if any(worker.pending is not None or worker.buffer for worker in workers):
            fail("worker state is not empty at the job barrier")
        return maximum_end, used
    finally:
        selector.close()


def quit_workers(workers: Sequence[PersistentWorker], deadline: float) -> None:
    for worker in workers:
        if worker.pending is not None or worker.buffer:
            fail("cannot quit a worker with unconsumed result state")
        _write_worker(worker, b"Q\n")
        worker.process.stdin.close()
    for worker in workers:
        try:
            returncode = worker.process.wait(
                timeout=_remaining(deadline, "native worker shutdown"))
        except subprocess.TimeoutExpired:
            fail("native worker on CPU {} did not exit after Q".format(
                worker.cpu))
        stdout = b""
        try:
            while True:
                block = os.read(worker.process.stdout.fileno(), 65536)
                if not block:
                    break
                stdout += block
        except BlockingIOError:
            pass
        stderr = _worker_stderr(worker)
        if returncode != 0 or stdout or stderr:
            fail("native worker on CPU {} failed exact Q shutdown{}".format(
                worker.cpu, ": " + stderr if stderr else ""))
        worker.process.stdout.close()
        worker.process.stderr.close()


def _strict_response_validator(
        contract: Mapping[str, Any], freeze: Mapping[str, Any], kind: str,
        description: Mapping[str, Any], window_start_ns: int,
        ) -> ResponseValidator:
    expected_schema = native_api.RECOVERY_RECORD_SCHEMA if kind == "recovery" \
        else native_api.TIMING_RECORD_SCHEMA
    payload_fields = contract_api.LEDGER_FIELDS if kind == "recovery" \
        else contract_api.TIMING_RECEIPT_FIELDS

    def validate(value: Mapping[str, Any], _line: bytes,
                 worker: PersistentWorker, job: Job) -> int:
        integer_fields = (
            "ordinal", "cpu", "worker_pid", "worker_process_start_ticks",
            "started_monotonic_ns", "finished_monotonic_ns",
        )
        if job.kind != kind or \
                set(value) != native_api.NATIVE_RECORD_FIELDS or \
                any(type(value.get(field)) is not int
                    for field in integer_fields) or \
                value.get("schema") != expected_schema or \
                value.get("ordinal") != job.ordinal or \
                value.get("cpu") != worker.cpu or \
                value.get("worker_pid") != worker.process.pid or \
                value.get("worker_process_start_ticks") != worker.start_ticks or \
                value.get("worker_binary_sha256") != \
                description["binary_sha256"]:
            fail("native response does not bind its commanded worker/job")
        start = value.get("started_monotonic_ns")
        end = value.get("finished_monotonic_ns")
        if type(start) is not int or type(end) is not int or \
                not window_start_ns <= start <= end:
            fail("native response has invalid monotonic provenance")
        for field in native_api.SHA256_FIELDS:
            if not _is_sha256(value.get(field)):
                fail("native response {} is not a SHA-256".format(field))
        payload = value.get("payload")
        if not isinstance(payload, dict) or set(payload) != payload_fields:
            fail("native response payload has an unexpected schema")
        try:
            ordinal, _ = native_api._expected_record_ordinal(
                contract, freeze, kind, "development", payload)
            work_sha256 = native_api._expected_work_sha256(
                kind, "development", job.ordinal, payload)
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        if ordinal != job.ordinal or value["work_sha256"] != work_sha256:
            fail("native response payload/work hash differs from its command")
        if kind == "recovery":
            frozen_arm = freeze["arms_by_name"].get(payload.get("arm"))
            if (frozen_arm is None or
                    payload.get("binary_sha256") != description["binary_sha256"] or
                    payload.get("arm_descriptor_sha256") !=
                    frozen_arm["arm_descriptor_sha256"]):
                fail("native recovery response differs from its frozen arm")
        else:
            for side in ("left", "right"):
                frozen_arm = freeze["arms_by_name"].get(
                    payload.get(side + "_arm"))
                if (frozen_arm is None or
                        payload.get(side + "_binary_sha256") !=
                        description["binary_sha256"] or
                        payload.get(side + "_arm_descriptor_sha256") !=
                        frozen_arm["arm_descriptor_sha256"]):
                    fail("native timing response differs from its frozen arm")
        return end

    return validate


def _recovery_jobs(contract: Mapping[str, Any], arm_count: int) -> List[Job]:
    cells = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    jobs = [
        Job("recovery", cell, arm, cell * arm_count + arm)
        for cell in range(cells) for arm in range(arm_count)
    ]
    if (cells != 360 or arm_count != len(EXPECTED_ARMS) or
            len(jobs) != cells * arm_count):
        fail("development recovery jobs differ from the frozen arm roster")
    return jobs


def _development_timing_shape(
        contract: Mapping[str, Any], panel_count: int,
        ) -> Tuple[int, int]:
    timing = contract.get("timing")
    if not isinstance(timing, Mapping):
        fail("development timing policy is missing")
    geometry = timing.get("execution_geometry")
    if geometry != contract_api.EXPECTED_TIMING_EXECUTION_GEOMETRY:
        fail("development timing execution geometry is not frozen")
    if geometry.get("timing_worker_count") != TIMING_WORKER_COUNT or \
            geometry.get("jobs_per_wave") != TIMING_WORKER_COUNT:
        fail("development timing geometry does not use exactly eight workers")
    domains = timing.get("domains")
    domain = domains.get("development") \
        if isinstance(domains, Mapping) else None
    repetitions = domain.get("paired_repetitions") \
        if isinstance(domain, Mapping) else None
    expected_cells = domain.get("expected_cells") \
        if isinstance(domain, Mapping) else None
    if (type(repetitions) is not int or repetitions <= 0 or
            repetitions % TIMING_WORKER_COUNT != 0):
        fail("development timing repetitions must be a positive multiple of eight")
    if (type(expected_cells) is not int or expected_cells <= 0 or
            expected_cells % repetitions != 0):
        fail("development timing cell cardinality is invalid")
    try:
        expected_panels = len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS]))
    except contract_api.ContractError as exc:
        fail(str(exc))
    if type(panel_count) is not int or panel_count != expected_panels:
        fail("development timing panel cardinality is invalid")
    return repetitions, expected_cells


def _timing_job_waves(
        contract: Mapping[str, Any], panel_count: int,
        ) -> List[Tuple[int, List[Job]]]:
    """Build homogeneous, barrier-delimited eight-job timing waves."""
    repetitions, expected_cells = _development_timing_shape(
        contract, panel_count)
    try:
        cells = list(contract_api.iter_timing_cells(
            contract, "development"))
    except contract_api.ContractError as exc:
        fail(str(exc))
    if len(cells) != expected_cells:
        fail("development timing cell cardinality differs from its contract")

    timing = contract["timing"]
    cell_fields = timing.get("cell_key")
    if (not isinstance(cell_fields, list) or
            any(not isinstance(field, str) for field in cell_fields) or
            len(cell_fields) != len(set(cell_fields)) or
            "replicate" not in cell_fields or "loss_seed" not in cell_fields):
        fail("development timing cell identity is invalid")
    stable_fields = [
        field for field in cell_fields
        if field not in ("replicate", "loss_seed")
    ]
    by_replicate: List[Dict[str, Tuple[int, Mapping[str, Any]]]] = [
        {} for _ in range(repetitions)
    ]
    identity_order: List[str] = []
    expected_field_set = set(cell_fields)
    for cell_ordinal, cell in enumerate(cells):
        if not isinstance(cell, Mapping) or set(cell) != expected_field_set:
            fail("development timing cell has an unexpected identity")
        replicate = cell.get("replicate")
        if (type(replicate) is not int or
                not 0 <= replicate < repetitions):
            fail("development timing cell has an invalid replicate")
        identity = contract_api.canonical_json({
            field: cell[field] for field in stable_fields
        })
        if identity in by_replicate[replicate]:
            fail("duplicate timing cell identity within one replicate")
        by_replicate[replicate][identity] = (cell_ordinal, cell)
        if replicate == 0:
            identity_order.append(identity)

    cells_per_replicate = expected_cells // repetitions
    if (len(identity_order) != cells_per_replicate or
            len(set(identity_order)) != cells_per_replicate):
        fail("development timing replicate zero has invalid cardinality")
    expected_identities = set(identity_order)
    for replicate, values in enumerate(by_replicate):
        if (len(values) != cells_per_replicate or
                set(values) != expected_identities):
            fail("timing cohort identity differs across replicate {}".format(
                replicate))

    waves: List[Tuple[int, List[Job]]] = []
    cohort_index = 0
    frozen_arms = [value[0] for value in EXPECTED_ARMS]
    for identity in identity_order:
        for panel in range(panel_count):
            for wave_index, first_replicate in enumerate(
                    range(0, repetitions, TIMING_WORKER_COUNT)):
                jobs = []
                for replicate in range(
                        first_replicate,
                        first_replicate + TIMING_WORKER_COUNT):
                    cell_ordinal, _ = by_replicate[replicate][identity]
                    jobs.append(Job(
                        "timing", cell_ordinal, panel,
                        cell_ordinal * panel_count + panel))
                # run_job_batch assigns job position r to worker slot
                # (r + rotation) modulo eight.  Derive that rotation from the
                # bounded contract helper and verify the complete wave before
                # allowing any native command to be issued.
                try:
                    rotation = contract_api.timing_worker_slot(
                        contract, "development", frozen_arms, cohort_index,
                        first_replicate)
                    for position, replicate in enumerate(range(
                            first_replicate,
                            first_replicate + TIMING_WORKER_COUNT)):
                        expected_slot = contract_api.timing_worker_slot(
                            contract, "development", frozen_arms,
                            cohort_index, replicate)
                        if expected_slot != (
                                position + rotation) % TIMING_WORKER_COUNT:
                            fail("timing worker-slot formula is not a rotation")
                except contract_api.ContractError as exc:
                    fail(str(exc))
                waves.append((rotation, jobs))
            cohort_index += 1

    expected_cohorts = cells_per_replicate * panel_count
    try:
        frozen_cohorts = contract_api.timing_cohort_count(
            contract, "development", frozen_arms)
    except contract_api.ContractError as exc:
        fail(str(exc))
    expected_waves = expected_cohorts * (
        repetitions // TIMING_WORKER_COUNT)
    expected_jobs = {
        (cell, panel, cell * panel_count + panel)
        for cell in range(expected_cells)
        for panel in range(panel_count)
    }
    actual_jobs = [
        (job.cell, job.item, job.ordinal)
        for _, jobs in waves for job in jobs
    ]
    if (cohort_index != expected_cohorts or
            expected_cohorts != frozen_cohorts or
            len(waves) != expected_waves or
            any(len(jobs) != TIMING_WORKER_COUNT for _, jobs in waves) or
            len(actual_jobs) != len(expected_jobs) or
            len(set(actual_jobs)) != len(actual_jobs) or
            set(actual_jobs) != expected_jobs):
        fail("development timing waves do not cover the exact frozen domain")
    return waves


def _run_native_jobs(
        contract: Mapping[str, Any], freezes: Mapping[str, Mapping[str, Any]],
        description: Mapping[str, Any], workers: Sequence[PersistentWorker],
        output_dir: Path, window_start_ns: int, deadline: float,
        ) -> Tuple[Mapping[str, Path], int]:
    if len(workers) != TIMING_WORKER_COUNT:
        fail("native short screen requires exactly eight timing workers")
    try:
        panel_count = len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS]))
    except contract_api.ContractError as exc:
        fail(str(exc))
    timing_waves = _timing_job_waves(contract, panel_count)
    paths = {
        "recovery": output_dir / "recovery-native-results.jsonl",
        "timing": output_dir / "timing-native-results.jsonl",
    }
    recovery_sink = AtomicLineSink(paths["recovery"])
    timing_sink = AtomicLineSink(paths["timing"])
    maximum_end = 0
    try:
        recovery_end, recovery_cpus = run_job_batch(
            workers, _recovery_jobs(contract, len(EXPECTED_ARMS)), 0,
            recovery_sink, deadline, _strict_response_validator(
                contract, freezes["recovery"], "recovery", description,
                window_start_ns))
        if recovery_cpus != set(worker.cpu for worker in workers):
            fail("recovery did not exercise every frozen CPU")
        recovery_sink.publish()
        maximum_end = max(maximum_end, recovery_end)

        expected_timing_cpus = set(worker.cpu for worker in workers)
        timing_validator = _strict_response_validator(
            contract, freezes["timing"], "timing", description,
            window_start_ns)
        for wave_ordinal, (rotation, jobs) in enumerate(timing_waves):
            end, used = run_job_batch(
                workers, jobs, rotation, timing_sink,
                deadline, timing_validator)
            if used != expected_timing_cpus:
                fail("timing wave {} did not use every frozen CPU".format(
                    wave_ordinal))
            maximum_end = max(maximum_end, end)
        timing_sink.publish()
        return paths, maximum_end
    finally:
        recovery_sink.abort()
        timing_sink.abort()


def _create_output_dir(path: Path) -> Path:
    try:
        path.mkdir(parents=True, exist_ok=False)
        resolved = path.resolve(strict=True)
        if not stat.S_ISDIR(resolved.stat().st_mode):
            fail("output path is not a directory")
        _fsync_directory(resolved.parent)
        return resolved
    except FileExistsError:
        fail("output directory already exists: {}".format(path))
    except OSError as exc:
        fail("cannot create output directory {}: {}".format(path, exc))
    return path


def run_short_screen(args: argparse.Namespace) -> Mapping[str, Any]:
    if (not math.isfinite(args.deadline_seconds) or
            not 0.0 < args.deadline_seconds <= MAX_WALL_SECONDS):
        fail("--deadline-seconds must be in (0,7200]")
    deadline = time.monotonic() + args.deadline_seconds
    contract = contract_api.load_contract(args.contract)
    description = describe_worker(args.worker, deadline)
    try:
        original_affinity = set(os.sched_getaffinity(0))
    except (AttributeError, OSError) as exc:
        fail("cannot inspect initial controller affinity: {}".format(exc))
    explicit = parse_cpu_list(args.cpus) if args.cpus is not None else None
    cpus, controller_cpu = select_cpu_layout(
        args.sampler_cpu, explicit, args.controller_cpu,
        affinity=original_affinity)
    # Recovery remains work-conserving on this same frozen roster.  Timing
    # uses one exact job per worker and a barrier after every homogeneous wave.
    if len(cpus) != TIMING_WORKER_COUNT:
        fail("native short screen requires exactly eight physical worker CPUs")
    _development_timing_shape(
        contract,
        len(contract_api.timing_panels(
            contract, [value[0] for value in EXPECTED_ARMS])))
    _preflight_sampler(
        args.sampler_pid, args.sampler_cpu,
        args.sampler_script, args.sampler_csv)
    source_commit = _git_head(deadline)
    _require_worker_source_commit(description, source_commit)
    output_dir = _create_output_dir(args.output_dir)

    trace_paths: Dict[str, Path] = {}
    trace_sha256s: Dict[str, str] = {}
    for kind in ("recovery", "timing"):
        path, digest = _emit_and_assemble_trace(
            contract, kind, description, output_dir, deadline)
        trace_paths[kind] = path
        trace_sha256s[kind] = digest
    if _git_head(deadline) != source_commit:
        fail("codec source commit changed before the result freeze")
    freezes = write_development_freezes(
        contract, description, cpus, controller_cpu, source_commit,
        trace_sha256s, output_dir)
    freeze_paths = {
        kind: output_dir / "{}-freeze.json".format(kind)
        for kind in ("recovery", "timing")
    }

    workers: List[PersistentWorker] = []
    clean_shutdown = False
    controller_pinned = False
    try:
        # Both validated freezes exist before the first R/T command.
        workers = spawn_workers(description, cpus, deadline)
        controller_pinned = True
        _pin_controller(controller_cpu)
        window_start_ns = choose_new_sampler_start(args.sampler_csv, deadline)
        native_paths, maximum_worker_end = _run_native_jobs(
            contract, freezes, description, workers, output_dir,
            window_start_ns, deadline)
        window_end_ns = _wait_for_sampler_sample(
            args.sampler_csv, deadline,
            at_or_after_ns=maximum_worker_end)
        if _git_head(deadline) != source_commit:
            fail("codec source commit changed during the native screen")
        sampler_path = output_dir / "sampler-attestation.json"
        try:
            native_api.write_sampler_attestation(
                args.sampler_pid, args.sampler_cpu,
                args.sampler_script, args.sampler_csv,
                window_start_ns, window_end_ns, sampler_path)
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))

        assembled: Dict[str, Mapping[str, Any]] = {}
        for kind in ("recovery", "timing"):
            result_path = output_dir / "{}-results.jsonl".format(kind)
            receipt_path = output_dir / "{}-execution.json".format(kind)
            try:
                assembled[kind] = native_api.assemble_results(
                    contract, kind, "development", freeze_paths[kind],
                    trace_paths[kind], native_paths[kind], sampler_path,
                    result_path, receipt_path, verify_live_sampler=True)
                _remaining(deadline, "assembling {} results".format(kind))
                native_api.validate_execution_receipt(
                    contract, kind, "development", freeze_paths[kind],
                    trace_paths[kind], native_paths[kind], result_path,
                    receipt_path, verify_live_sampler=True)
                _remaining(deadline, "validating {} execution".format(kind))
            except (native_api.NativeEvidenceError,
                    contract_api.ContractError) as exc:
                fail(str(exc))
        quit_workers(workers, deadline)
        clean_shutdown = True

        recovery_receipt = assembled["recovery"]["execution_receipt"]
        timing_receipt = assembled["timing"]["execution_receipt"]
        thermal = timing_receipt["thermal"]
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "complete",
            "output_dir": str(output_dir),
            "source_git_commit": source_commit,
            "contract_sha256": contract_api.contract_sha256(contract),
            "worker_binary_sha256": description["binary_sha256"],
            "controller_cpu": controller_cpu,
            "worker_cpus": list(cpus),
            "recovery_records": recovery_receipt["record_count"],
            "timing_records": timing_receipt["record_count"],
            "recovery_freeze_sha256":
                recovery_receipt["freeze_manifest_sha256"],
            "timing_freeze_sha256": timing_receipt["freeze_manifest_sha256"],
            "recovery_result_sha256":
                recovery_receipt["result_stream_sha256"],
            "timing_result_sha256": timing_receipt["result_stream_sha256"],
            "thermal_samples": thermal["sample_count"],
            "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "dimm_max_millic": thermal["dimm_max_millic"],
        }
        _atomic_write_object(output_dir / "run-summary.json", summary)
        return summary
    finally:
        if workers and not clean_shutdown:
            terminate_workers(workers)
        if controller_pinned:
            _restore_controller_affinity(original_affinity)


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument(
        "--worker", type=Path,
        default=Path("build/wirehair_wh2_contract_worker"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sampler-pid", type=int, required=True)
    parser.add_argument("--sampler-cpu", type=int, required=True)
    parser.add_argument("--sampler-script", type=Path, required=True)
    parser.add_argument("--sampler-csv", type=Path, required=True)
    parser.add_argument(
        "--cpus", help="strictly increasing comma-separated logical CPUs")
    parser.add_argument(
        "--controller-cpu", type=int,
        help="dedicated non-worker, non-sampler physical core")
    parser.add_argument(
        "--deadline-seconds", type=float, default=MAX_WALL_SECONDS)
    args = parser.parse_args(argv if argv else None)
    try:
        if (not math.isfinite(args.deadline_seconds) or
                not 0.0 < args.deadline_seconds <= MAX_WALL_SECONDS):
            fail("--deadline-seconds must be in (0,7200]")
        with _hard_wall(args.deadline_seconds):
            summary = run_short_screen(args)
    except (RunnerError, native_api.NativeEvidenceError,
            contract_api.ContractError, OSError, UnicodeError) as exc:
        print("wh2 native short screen: {}".format(exc), file=sys.stderr)
        return 1
    print(contract_api.canonical_json(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
