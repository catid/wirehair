#!/usr/bin/env python3
"""Seal, execute, and strictly reduce the WH2 degree-balanced all-K holdout.

Preparation and execution are deliberately separate.  Preparation builds two
fresh binaries from exact commits: the production base itself, and the
reviewed degree-balanced candidate.  It then proves that the candidate binary
with the experiment flag disabled matches the base binary on every
deterministic preamble, outcome, dimension, seed-attempt, and work field over a
bounded matrix (excluding only six measured timing fields), freezes independent
holdout seeds and an exact paired all-K ledger, and emits a hash-locked launch
receipt.  Only ``launch --execute`` can start benchmark work; preparation,
selftests, and preflight never do.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


sys.dont_write_bytecode = True
getcontext().prec = 50

SCHEMA = "wirehair.wh2.degree_balanced_allk_holdout.v1"
BASE_COMMIT = "3c8a7a31e0e0b592f906d52bf7fef0bafebbf2f9"
SAMPLER_NAME = "wirehair_expo_thermal_sampler.py"
BASE_BINARY_NAME = "wirehair_v2_bench.base"
CANDIDATE_BINARY_NAME = "wirehair_v2_bench.degree_balanced"
CONTROLLER_NAME = "wh2_degree_balanced_allk_campaign.py"

# These seeds are an independent sealed holdout.  They are the first 64 bits
# of SHA-256("wirehair.wh2.degree-balanced.allk.holdout.v1|seed|i") and are
# intentionally disjoint from the three development seeds.
SEEDS = (
    "0x475b1dd3864534a3",
    "0x65eda454ed6be8ec",
    "0x6fd9d0c558cad523",
)
DEVELOPMENT_SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
ARMS = ("base", "degree_balanced")
K_LO = 2
K_HI = 64_000
CHUNK_MAX = 256
DEFAULT_WORKERS = os.cpu_count() or 120
STAGES = tuple(
    "seed{}_{}".format(seed_index, schedule)
    for seed_index in range(len(SEEDS)) for schedule in SCHEDULES
)
MAX_STDOUT_BYTES_PER_TASK = 96 * 1024
MAX_STDERR_BYTES_PER_TASK = 16 * 1024
TASK_TIMEOUT_SECONDS = 300
MAX_STAGE_SECONDS = 1800
MAX_THERMAL_BYTES_PER_SEGMENT = 2 * 1024 * 1024
MAX_ATTEMPTS_PER_STAGE = 3
MAX_RUNTIME_OUTPUT_BYTES = 3 * 1024 * 1024 * 1024
CPU_BUSY_FLOOR = Decimal("95")
THERMAL_INTERVAL_SECONDS = Decimal("1")
THERMAL_MIN_GAP_SECONDS = Decimal("0.5")
THERMAL_MAX_GAP_SECONDS = Decimal("2.5")
THERMAL_READY_SAMPLES = 2
THERMAL_READY_TIMEOUT_SECONDS = 12.0
THERMAL_END_TIMEOUT_SECONDS = 5.0
BUILD_TARGET = "wirehair_v2_bench"
EXPECTED_CANDIDATE_PATHS = (
    "codec/PrecodeTest.cpp",
    "codec/V2BenchCliTest.cmake",
    "codec/WirehairV2Bench.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h",
    "codec/fuzz/V2ProfileFuzz.cpp",
)
ALLOWED_CANDIDATE_UNTRACKED = (
    "build-balanced-asan/",
    "build-balanced-clang-strict/",
    "build-balanced-strict/",
    "build-balanced/",
)
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
NONDETERMINISTIC_REPLAY_FIELDS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
)
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)


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


def load_json(path: Path) -> Dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
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
    tracked_status = git_output(
        repo, "status", "--porcelain=v1", "--untracked-files=no")
    if tracked_status:
        die("candidate preparation requires no tracked source changes")
    status = git_output(
        repo, "status", "--porcelain=v1", "--untracked-files=normal")
    untracked: List[str] = []
    for raw in status.decode("utf-8").splitlines():
        if not raw.startswith("?? "):
            die("candidate status contains a non-untracked entry")
        untracked.append(raw[3:])
    if any(path not in ALLOWED_CANDIDATE_UNTRACKED for path in untracked):
        die("candidate worktree contains an unreviewed untracked path")
    head = git_output(repo, "rev-parse", "HEAD").decode("ascii").strip()
    tree = git_output(repo, "rev-parse", "HEAD^{tree}").decode("ascii").strip()
    if re.fullmatch(r"[0-9a-f]{40}", head) is None or \
            re.fullmatch(r"[0-9a-f]{40}", tree) is None:
        die("source commit/tree identity is malformed")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", BASE_COMMIT, head],
        cwd=str(repo), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    if ancestor.returncode != 0:
        die("candidate source does not descend from the sealed production base")
    changed = tuple(filter(None, git_output(
        repo, "diff", "--name-only", "{}..{}".format(BASE_COMMIT, head),
    ).decode("utf-8").splitlines()))
    if changed != EXPECTED_CANDIDATE_PATHS:
        die("candidate/base changed-file set is not the reviewed six-file patch")
    precode = git_output(repo, "show", "{}:codec/WirehairV2Precode.cpp".format(head))
    bench = git_output(repo, "show", "{}:codec/WirehairV2Bench.cpp".format(head))
    if b"AppendDegreeBalancedStaircaseEdges" not in precode or \
            b"--degree-balanced-staircase" not in bench:
        die("candidate commit lacks the reviewed degree-balanced implementation/CLI")
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", head)
    if not manifest or not manifest.endswith(b"\n"):
        die("source tree manifest is empty or noncanonical")
    patch = git_output(repo, "diff", "--binary", BASE_COMMIT, head, "--", *changed)
    return {
        "commit": head,
        "tree_oid": tree,
        "tracked_status_porcelain_sha256": sha256_bytes(tracked_status),
        "untracked_status_porcelain_sha256": sha256_bytes(status),
        "allowed_untracked_paths": untracked,
        "tree_manifest_sha256": sha256_bytes(manifest),
        "tree_manifest_bytes": len(manifest),
        "base_commit": BASE_COMMIT,
        "changed_paths": list(changed),
        "base_to_candidate_patch_sha256": sha256_bytes(patch),
        "base_to_candidate_patch_bytes": len(patch),
    }


def commit_identity(repo: Path, commit: str) -> Dict[str, object]:
    resolved = git_output(repo, "rev-parse", "{}^{{commit}}".format(commit)).decode(
        "ascii").strip()
    tree = git_output(repo, "rev-parse", "{}^{{tree}}".format(resolved)).decode(
        "ascii").strip()
    if re.fullmatch(r"[0-9a-f]{40}", resolved) is None or \
            re.fullmatch(r"[0-9a-f]{40}", tree) is None:
        die("build commit/tree identity is malformed")
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", resolved)
    archive = git_output(repo, "archive", "--format=tar", resolved)
    return {
        "commit": resolved, "tree_oid": tree,
        "tree_manifest_sha256": sha256_bytes(manifest),
        "tree_manifest_bytes": len(manifest),
        "source_archive_sha256": sha256_bytes(archive),
        "manifest": manifest, "archive": archive,
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


def fresh_build_binary(
    repo: Path,
    frozen: Path,
    workers: int,
    *,
    commit: str,
    role: str,
    binary_name: str,
) -> Dict[str, object]:
    """Build one exact commit and freeze independently checkable provenance."""
    before = clean_source_identity(repo)
    source_record = commit_identity(repo, commit)
    manifest = source_record.pop("manifest")
    archive = source_record.pop("archive")
    tools = {
        "cmake": executable_identity("cmake", ("--version",)),
        "ninja": executable_identity("ninja", ("--version",)),
        "readelf": executable_identity("readelf", ("--version",)),
    }
    source_record["repository_head_at_build"] = before["commit"]
    source_record["repository_clean_tracked_status_sha256"] = before[
        "tracked_status_porcelain_sha256"]
    build_evidence = frozen / ("build_" + role)
    build_evidence.mkdir()
    (build_evidence / "source_tree.txt").write_bytes(manifest)
    with tempfile.TemporaryDirectory(
            prefix="wirehair-degree-balanced-{}-build-".format(role)) as temporary:
        source = Path(temporary) / "source"
        source.mkdir()
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
        candidate = build / "codec" / "wirehair_v2_bench"
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
        for name, source in evidence_sources.items():
            shutil.copyfile(str(source), str(build_evidence / name))
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
        shutil.copyfile(str(candidate), str(frozen / binary_name))
        os.chmod(str(frozen / binary_name), 0o555)

    after = clean_source_identity(repo)
    if after != before:
        die("source commit/tree/worktree changed during the fresh build")
    evidence_files = sorted(path for path in build_evidence.iterdir() if path.is_file())
    record = {
        "schema": SCHEMA + ".fresh_build",
        "role": role,
        "source": source_record,
        "candidate_repository_snapshot": before,
        "configure_command": configure,
        "build_command": build_command,
        "workers": workers,
        "tools": tools,
        "compiler": {
            "path": str(compiler), "sha256": sha256_file(compiler),
            "version_stdout_sha256": sha256_bytes(compiler_stdout),
            "version_stderr_sha256": sha256_bytes(compiler_stderr),
        },
        "binary_name": binary_name,
        "binary_sha256": sha256_file(frozen / binary_name),
        "binary_elf_build_id": build_ids[0].decode("ascii"),
        "evidence": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in evidence_files
        },
        "fresh_build": True,
    }
    write_once(frozen / ("build_{}_receipt.json".format(role)), canonical_json(record))
    return record


def cell_tuple(record: Mapping[str, object]) -> Tuple[int, int, str]:
    try:
        key = (int(record["K"]), int(record["seed_index"]), str(record["schedule"]))
    except (KeyError, TypeError, ValueError) as exc:
        die("malformed all-K cell record: {}".format(exc))
    if not (K_LO <= key[0] <= K_HI and 0 <= key[1] < len(SEEDS) and
            key[2] in SCHEDULES):
        die("out-of-domain all-K cell key {}".format(key))
    return key


def cell_record(key: Tuple[int, int, str]) -> Dict[str, object]:
    return {"K": key[0], "seed_index": key[1], "schedule": key[2]}


def expected_preamble(task: Mapping[str, object]) -> Dict[str, str]:
    result = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": SEEDS[int(task["seed_index"])], "completion": "mixed",
        "mixed_period": "244", "mixed_gf256_rows": "10",
        "mixed_gf16_rows": "2", "mixed_geometry": "frozen",
        "mixed_residue_skew": "0", "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0", "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": "0",
        "mixed_grouped_gf256_hash_seed": "0x0",
        "mixed_grouped_final_h_a_columns": "0",
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
    if task["arm"] == "degree_balanced":
        result["degree_balanced_staircase"] = "1"
    return result


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


def parse_benchmark_csv(data: bytes, task: Mapping[str, object]) -> List[Dict[str, str]]:
    if len(data) > int(task.get("stdout_max_bytes", MAX_STDOUT_BYTES_PER_TASK)):
        die("benchmark stdout exceeds the sealed per-task bound")
    try:
        text_value = data.decode("ascii")
    except UnicodeDecodeError as exc:
        die("benchmark output is not ASCII: {}".format(exc))
    lines = text_value.splitlines(keepends=True)
    if not lines or not all(line.endswith("\n") for line in lines) or \
            any("\r" in line or line == "\n" for line in lines):
        die("benchmark output is empty or lacks final LF")
    preamble = parse_preamble(lines[0][:-1])
    if preamble != expected_preamble(task):
        die("benchmark preamble differs from the sealed task")
    reader = csv.DictReader(lines[1:])
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        die("benchmark CSV header mismatch")
    rows = list(reader)
    if any(None in row or set(row) != set(CSV_FIELDS) or
           any(value is None for value in row.values()) for row in rows):
        die("benchmark CSV contains missing or extra fields")
    try:
        observed_Ks = [int(row["N"]) for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        die("benchmark CSV has an invalid K ledger: {}".format(exc))
    if observed_Ks != list(task["Ks"]):
        die("benchmark CSV K ledger mismatch")
    fixed = {
        "bb": "64", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "0", "trials": "1", "active_packet_peel_seed_xor": "0x0",
    }
    for row in rows:
        if any(row.get(key) != value for key, value in fixed.items()):
            die("benchmark fixed CSV field mismatch")
        validate_record(row)
    return rows


def task_argv(
    binary: str,
    arm: str,
    seed_index: int,
    schedule: str,
    Ks: Sequence[int],
) -> List[str]:
    if arm not in ARMS or schedule not in SCHEDULES or not Ks:
        die("invalid task argv request")
    argv = [
        binary, "precodefail", "--bb-list", "64", "--overhead", "0",
        "--trials", "1", "--threads", "1", "--loss", "0.50",
        "--completion", "mixed", "--mix-count", "2",
        "--packet-peel-seed-xor", "0", "--N", ",".join(map(str, Ks)),
        "--seed", SEEDS[seed_index], "--schedule", schedule,
        "--mixed-period", "244", "--mixed-geometry", "frozen",
        "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
        "--mixed-grouped-gf256-rows", "0", "--mixed-residue-buckets", "auto",
        "--binary-dense-two-anchor",
    ]
    if arm == "degree_balanced":
        argv.append("--degree-balanced-staircase")
    return argv


def smoke_frozen_binaries(
    base_binary: Path,
    candidate_binary: Path,
    frozen: Path,
) -> Dict[str, object]:
    smoke_dir = frozen / "smoke"
    smoke_dir.mkdir()
    smoke_Ks = [2, 257, 64_000]
    records: List[dict] = []
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            case = len(records)
            stem = "case{:02d}.seed{}.{}".format(case, seed_index, schedule)
            task_base = {
                "arm": "base", "seed_index": seed_index, "schedule": schedule,
                "Ks": smoke_Ks, "stdout_max_bytes": MAX_STDOUT_BYTES_PER_TASK,
            }
            task_candidate = dict(task_base, arm="degree_balanced")
            invocations = (
                ("base", base_binary, task_base),
                ("candidate_base_replay", candidate_binary, task_base),
                ("degree_balanced", candidate_binary, task_candidate),
            )
            outputs: Dict[str, bytes] = {}
            semantic_rows: Dict[str, List[dict]] = {}
            invocation_receipts: Dict[str, dict] = {}
            for name, binary, task in invocations:
                argv = task_argv(
                    str(binary), str(task["arm"]), seed_index, schedule, smoke_Ks)
                stdout, stderr = run_captured(
                    argv, cwd=frozen, timeout=180.0,
                    label="bounded {} smoke seed{} {}".format(name, seed_index, schedule),
                )
                if stderr or len(stderr) > MAX_STDERR_BYTES_PER_TASK:
                    die("bounded binary smoke produced stderr")
                parsed = parse_benchmark_csv(stdout, task)
                stdout_path = smoke_dir / (stem + "." + name + ".stdout")
                stderr_path = smoke_dir / (stem + "." + name + ".stderr")
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(stderr)
                outputs[name] = stdout
                semantic_rows[name] = [
                    {
                        field: row[field] for field in CSV_FIELDS
                        if field not in NONDETERMINISTIC_REPLAY_FIELDS
                    }
                    for row in parsed
                ]
                invocation_receipts[name] = {
                    "argv": ["frozen/" + binary.name, *argv[1:]],
                    "stdout": str(stdout_path.relative_to(frozen.parent)),
                    "stdout_sha256": sha256_bytes(stdout),
                    "stdout_bytes": len(stdout),
                    "stderr": str(stderr_path.relative_to(frozen.parent)),
                    "stderr_sha256": sha256_bytes(stderr),
                    "outcomes": [classify(row) for row in parsed],
                }
            if semantic_rows["base"] != semantic_rows["candidate_base_replay"]:
                die("candidate flag-off binary does not semantically replay the production base")
            records.append({
                "case": case, "seed_index": seed_index, "seed": SEEDS[seed_index],
                "schedule": schedule, "Ks": smoke_Ks,
                "base_replay_semantic_identical": True,
                "base_replay_byte_identical": (
                    outputs["base"] == outputs["candidate_base_replay"]),
                "replay_excluded_fields": list(NONDETERMINISTIC_REPLAY_FIELDS),
                "base_semantic_sha256": sha256_bytes(canonical_json(
                    semantic_rows["base"])),
                "candidate_base_semantic_sha256": sha256_bytes(canonical_json(
                    semantic_rows["candidate_base_replay"])),
                "invocations": invocation_receipts,
            })
    record = {
        "schema": SCHEMA + ".binary_smoke", "bounded": True,
        "base_commit": BASE_COMMIT,
        "base_binary_sha256": sha256_file(base_binary),
        "candidate_binary_sha256": sha256_file(candidate_binary),
        "cases": records,
    }
    write_once(frozen / "binary_smoke.json", canonical_json(record))
    return record


def build_tasks(
    base_binary: str,
    candidate_binary: str,
    *,
    k_lo: int = K_LO,
    k_hi: int = K_HI,
    chunk_max: int = CHUNK_MAX,
) -> List[dict]:
    if not (K_LO <= k_lo <= k_hi <= K_HI and 1 <= chunk_max <= CHUNK_MAX):
        die("invalid task-ledger bounds")
    binary_by_arm = {
        "base": base_binary,
        "degree_balanced": candidate_binary,
    }
    tasks: List[dict] = []
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            stage = "seed{}_{}".format(seed_index, schedule)
            for chunk_index, lo in enumerate(range(k_lo, k_hi + 1, chunk_max)):
                Ks = list(range(lo, min(k_hi, lo + chunk_max - 1) + 1))
                for arm in ARMS:
                    job = len(tasks)
                    output_name = (
                        "{:05d}.{}.seed{}.{}.chunk{:03d}.csv"
                        .format(job, arm, seed_index, schedule, chunk_index)
                    )
                    task = {
                        "job": job, "arm": arm, "stage": stage,
                        "seed_index": seed_index, "seed": SEEDS[seed_index],
                        "schedule": schedule, "chunk_index": chunk_index,
                        "K_lo": Ks[0], "K_hi": Ks[-1], "Ks": Ks,
                        "output_name": output_name,
                        "stdout_max_bytes": MAX_STDOUT_BYTES_PER_TASK,
                        "stderr_max_bytes": MAX_STDERR_BYTES_PER_TASK,
                        "timeout_seconds": TASK_TIMEOUT_SECONDS,
                    }
                    task["argv"] = task_argv(
                        binary_by_arm[arm], arm, seed_index, schedule, Ks)
                    tasks.append(task)
    return tasks


def independent_task_audit(
    tasks: Sequence[dict],
    *,
    k_lo: int = K_LO,
    k_hi: int = K_HI,
) -> Dict[str, object]:
    expected_Ks = list(range(k_lo, k_hi + 1))
    strata: Dict[Tuple[str, int, str], List[int]] = {
        (arm, seed_index, schedule): []
        for arm in ARMS for seed_index in range(len(SEEDS))
        for schedule in SCHEDULES
    }
    logical = hashlib.sha256()
    outputs: Set[str] = set()
    for expected_job, task in enumerate(tasks):
        if task.get("job") != expected_job:
            die("task job sequence mismatch")
        arm = str(task.get("arm"))
        seed_index = int(task.get("seed_index", -1))
        schedule = str(task.get("schedule"))
        stage = "seed{}_{}".format(seed_index, schedule)
        Ks = [int(value) for value in task.get("Ks", [])]
        if arm not in ARMS or not 0 <= seed_index < len(SEEDS) or \
                schedule not in SCHEDULES or task.get("seed") != SEEDS[seed_index] or \
                task.get("stage") != stage or not Ks or len(Ks) > CHUNK_MAX or \
                Ks != list(range(Ks[0], Ks[-1] + 1)) or \
                task.get("K_lo") != Ks[0] or task.get("K_hi") != Ks[-1]:
            die("task semantic receipt mismatch at job {}".format(expected_job))
        expected_binary = task["argv"][0]
        if task["argv"] != task_argv(
                expected_binary, arm, seed_index, schedule, Ks):
            die("task argv differs from regenerated architecture")
        if task.get("stdout_max_bytes") != MAX_STDOUT_BYTES_PER_TASK or \
                task.get("stderr_max_bytes") != MAX_STDERR_BYTES_PER_TASK or \
                task.get("timeout_seconds") != TASK_TIMEOUT_SECONDS:
            die("task output bound differs")
        name = str(task.get("output_name"))
        if name in outputs:
            die("duplicate task output name")
        outputs.add(name)
        strata[(arm, seed_index, schedule)].extend(Ks)
        for K in Ks:
            logical.update("{},{},{},{}\n".format(
                arm, K, seed_index, schedule).encode("ascii"))
    if len(tasks) % 2:
        die("paired task ledger has odd cardinality")
    for offset in range(0, len(tasks), 2):
        base, candidate = tasks[offset:offset + 2]
        if base["arm"] != "base" or candidate["arm"] != "degree_balanced" or \
                any(base[key] != candidate[key] for key in (
                    "stage", "seed_index", "seed", "schedule", "chunk_index",
                    "K_lo", "K_hi", "Ks")):
            die("adjacent task pair is not an exact architecture A/B")
    expected_digest = hashlib.sha256()
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            for K in expected_Ks:
                for arm in ARMS:
                    expected_digest.update("{},{},{},{}\n".format(
                        arm, K, seed_index, schedule).encode("ascii"))
    # The generated task order is chunk then arm, while the direct audit order
    # above is cell then arm.  Compare sets through independent stratum replay
    # and receipt both deterministic orderings rather than conflating them.
    if any(values != expected_Ks for values in strata.values()):
        die("task ledger does not exactly cover every all-K stratum")
    cells_per_arm = len(expected_Ks) * len(SEEDS) * len(SCHEDULES)
    total_cells = cells_per_arm * len(ARMS)
    per_job_bound = (
        MAX_STDOUT_BYTES_PER_TASK + MAX_STDERR_BYTES_PER_TASK + 4096)
    max_stage_tasks = max(
        sum(task["stage"] == stage for task in tasks) for stage in STAGES)
    max_runtime = (
        len(STAGES) * MAX_ATTEMPTS_PER_STAGE * max_stage_tasks * per_job_bound +
        len(tasks) * per_job_bound +
        len(STAGES) * MAX_ATTEMPTS_PER_STAGE *
            MAX_THERMAL_BYTES_PER_SEGMENT
    )
    if max_runtime > MAX_RUNTIME_OUTPUT_BYTES:
        die("sealed task ledger exceeds the campaign output budget")
    return {
        "schema": SCHEMA + ".task_audit", "audit": "OK",
        "audit_method": "independent per-arm stratum replay plus adjacent A/B pairing",
        "K_range": [k_lo, k_hi], "K_count": len(expected_Ks),
        "arms": list(ARMS), "seeds": list(SEEDS), "schedules": list(SCHEDULES),
        "stages": list(STAGES), "task_count": len(tasks),
        "task_pairs": len(tasks) // 2, "cells_per_arm": cells_per_arm,
        "total_cells": total_cells, "chunk_max": CHUNK_MAX,
        "generated_logical_cells_sha256": logical.hexdigest(),
        "direct_cell_major_sha256": expected_digest.hexdigest(),
        "max_stdout_bytes_per_task": MAX_STDOUT_BYTES_PER_TASK,
        "max_stderr_bytes_per_task": MAX_STDERR_BYTES_PER_TASK,
        "max_attempts_per_stage": MAX_ATTEMPTS_PER_STAGE,
        "max_stage_tasks": max_stage_tasks,
        "task_timeout_seconds": TASK_TIMEOUT_SECONDS,
        "max_stage_seconds": MAX_STAGE_SECONDS,
        "max_thermal_bytes_per_segment": MAX_THERMAL_BYTES_PER_SEGMENT,
        "max_runtime_output_bytes": max_runtime,
        "campaign_output_budget_bytes": MAX_RUNTIME_OUTPUT_BYTES,
    }


def artifact_map(root: Path) -> Dict[str, Path]:
    result = {
        "controller": root / "frozen" / CONTROLLER_NAME,
        "sampler": root / "frozen" / SAMPLER_NAME,
        "base_binary": root / "frozen" / BASE_BINARY_NAME,
        "candidate_binary": root / "frozen" / CANDIDATE_BINARY_NAME,
        "base_build_receipt": root / "frozen" / "build_base_receipt.json",
        "candidate_build_receipt": root / "frozen" / "build_candidate_receipt.json",
        "candidate_patch": root / "frozen" / "base_to_candidate.patch",
        "binary_smoke": root / "frozen" / "binary_smoke.json",
        "task_audit": root / "independent_task_audit.json",
        "tasks": root / "tasks_manifest.jsonl",
        "design": root / "design.json",
    }
    evidence_names = (
        "source_tree.txt", "CMakeCache.txt", "build.ninja",
        "compile_commands.json", "configure.stdout", "configure.stderr",
        "build.stdout", "build.stderr", "compiler.stdout", "compiler.stderr",
        "readelf-notes.stdout", "readelf-notes.stderr",
    )
    for role in ("base", "candidate"):
        for name in evidence_names:
            result["build_{}_{}".format(role, name)] = (
                root / "frozen" / ("build_" + role) / name)
    for case in range(len(SEEDS) * len(SCHEDULES)):
        seed_index = case // len(SCHEDULES)
        schedule = SCHEDULES[case % len(SCHEDULES)]
        stem = "case{:02d}.seed{}.{}".format(case, seed_index, schedule)
        for invocation in ("base", "candidate_base_replay", "degree_balanced"):
            for stream in ("stdout", "stderr"):
                key = "smoke_{}_{}_{}".format(case, invocation, stream)
                result[key] = root / "frozen" / (
                    "smoke/{0}.{1}.{2}".format(stem, invocation, stream))
    return result


def command_prepare(args: argparse.Namespace) -> None:
    root = Path(args.output_root).resolve()
    candidate_repo = Path(args.candidate_repo).resolve()
    sampler = Path(args.sampler).resolve()
    controller = Path(__file__).resolve()
    if root.exists():
        die("output root already exists: {}".format(root))
    staging = root.parent / (root.name + ".prepare.{}".format(os.getpid()))
    if staging.exists():
        die("preparation staging path already exists: {}".format(staging))
    if not candidate_repo.is_dir():
        die("candidate repository is not a directory")
    for path, label in ((sampler, "sampler"), (controller, "controller")):
        if not path.is_file() or path.is_symlink():
            die("{} is not a regular file: {}".format(label, path))
    source_snapshot = clean_source_identity(candidate_repo)
    candidate_commit = str(source_snapshot["commit"])
    if candidate_commit == BASE_COMMIT:
        die("candidate commit is identical to the production base")
    patch = git_output(
        candidate_repo, "diff", "--binary", BASE_COMMIT, candidate_commit,
        "--", *EXPECTED_CANDIDATE_PATHS,
    )
    if sha256_bytes(patch) != source_snapshot["base_to_candidate_patch_sha256"]:
        die("candidate patch changed between source audit steps")
    try:
        staging.mkdir(parents=True)
        for directory in (
                "frozen", "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (staging / directory).mkdir()
        frozen = staging / "frozen"
        shutil.copyfile(str(controller), str(frozen / CONTROLLER_NAME))
        shutil.copyfile(str(sampler), str(frozen / SAMPLER_NAME))
        (frozen / "base_to_candidate.patch").write_bytes(patch)
        os.chmod(str(frozen / CONTROLLER_NAME), 0o555)
        os.chmod(str(frozen / SAMPLER_NAME), 0o555)

        base_build = fresh_build_binary(
            candidate_repo, frozen, int(args.workers), commit=BASE_COMMIT,
            role="base", binary_name=BASE_BINARY_NAME,
        )
        candidate_build = fresh_build_binary(
            candidate_repo, frozen, int(args.workers), commit=candidate_commit,
            role="candidate", binary_name=CANDIDATE_BINARY_NAME,
        )
        if base_build["source"]["commit"] != BASE_COMMIT or \
                candidate_build["source"]["commit"] != candidate_commit:
            die("fresh builds are not bound to the intended commits")
        if candidate_build["candidate_repository_snapshot"] != source_snapshot:
            die("candidate repository identity changed during preparation")
        smoke = smoke_frozen_binaries(
            frozen / BASE_BINARY_NAME, frozen / CANDIDATE_BINARY_NAME, frozen)
        if smoke["base_binary_sha256"] != base_build["binary_sha256"] or \
                smoke["candidate_binary_sha256"] != candidate_build["binary_sha256"]:
            die("binary smoke is not bound to both fresh builds")

        tasks = build_tasks(
            str(root / "frozen" / BASE_BINARY_NAME),
            str(root / "frozen" / CANDIDATE_BINARY_NAME),
        )
        task_audit = independent_task_audit(tasks)
        write_once(staging / "tasks_manifest.jsonl", json_lines(tasks))
        write_once(staging / "independent_task_audit.json", canonical_json(task_audit))

        design = {
            "schema": SCHEMA + ".design", "root": str(root),
            "base_commit": BASE_COMMIT,
            "base_tree_oid": base_build["source"]["tree_oid"],
            "candidate_commit": candidate_commit,
            "candidate_tree_oid": candidate_build["source"]["tree_oid"],
            "base_to_candidate_patch_sha256": sha256_bytes(patch),
            "base_to_candidate_patch_bytes": len(patch),
            "candidate_changed_paths": list(EXPECTED_CANDIDATE_PATHS),
            "base_build_receipt_sha256": sha256_file(
                frozen / "build_base_receipt.json"),
            "candidate_build_receipt_sha256": sha256_file(
                frozen / "build_candidate_receipt.json"),
            "base_binary_sha256": base_build["binary_sha256"],
            "candidate_binary_sha256": candidate_build["binary_sha256"],
            "binary_smoke_sha256": sha256_file(frozen / "binary_smoke.json"),
            "base_replay_policy": (
                "candidate flag-off output must match the exact base-commit "
                "binary on every deterministic preamble, outcome, dimension, "
                "seed-attempt, and work field; only six measured timing fields "
                "are excluded"),
            "arms": {
                "base": {
                    "binary": "frozen/" + BASE_BINARY_NAME,
                    "degree_balanced_staircase": False,
                    "role": "exact production-base commit replay",
                },
                "degree_balanced": {
                    "binary": "frozen/" + CANDIDATE_BINARY_NAME,
                    "degree_balanced_staircase": True,
                    "role": "raw corrected degree-balanced staircase",
                },
            },
            "architecture_common": {
                "completion": "mixed", "mix_count": 2, "period": 244,
                "geometry": "frozen", "gf256_rows": 10, "gf16_rows": 2,
                "grouped_gf256_rows": 0, "binary_dense_two_anchor": True,
                "binary_dense_two_anchor_phase": 0,
                "packet_peel_seed_xor": 0, "packet_peel_seed_table": "none",
            },
            "seed_policy": (
                "independent SHA-256-derived holdout seeds; no experiment-specific "
                "K or seed fixes; inherited production profile selection is paired"),
            "seed_derivation_domain": (
                "wirehair.wh2.degree-balanced.allk.holdout.v1|seed|i"),
            "seeds": list(SEEDS), "development_seeds": list(DEVELOPMENT_SEEDS),
            "development_seed_overlap": sorted(set(SEEDS) & set(DEVELOPMENT_SEEDS)),
            "schedules": list(SCHEDULES), "K_range": [K_LO, K_HI],
            "K_count": K_HI - K_LO + 1,
            "bb": 64, "loss": 0.5, "overhead": 0, "trials": 1,
            "hard_loss_policy": (
                "50 percent loss at zero overhead under burst, adversarial, and "
                "repair-only schedules"),
            "cells_per_arm": (K_HI - K_LO + 1) * len(SEEDS) * len(SCHEDULES),
            "total_cells": task_audit["total_cells"],
            "task_count": len(tasks), "task_pairs": len(tasks) // 2,
            "chunk_max": CHUNK_MAX, "stage_order": list(STAGES),
            "worker_count": args.workers, "all_cpu_policy": (
                "one worker per available logical CPU unless explicitly sealed lower"),
            "thermal_single_sampler_required": True,
            "thermal_cpu_busy_floor_pct": str(CPU_BUSY_FLOOR),
            "thermal_interval_seconds": str(THERMAL_INTERVAL_SECONDS),
            "thermal_min_gap_seconds": str(THERMAL_MIN_GAP_SECONDS),
            "thermal_max_gap_seconds": str(THERMAL_MAX_GAP_SECONDS),
            "thermal_ready_samples": THERMAL_READY_SAMPLES,
            "retry_policy": (
                "stage-atomic non-selective retry after failed/interrupted segment"),
            "max_attempts_per_stage": MAX_ATTEMPTS_PER_STAGE,
            "task_timeout_seconds": TASK_TIMEOUT_SECONDS,
            "max_stage_seconds": MAX_STAGE_SECONDS,
            "max_thermal_bytes_per_segment": MAX_THERMAL_BYTES_PER_SEGMENT,
            "output_policy": {
                "stdout_max_bytes_per_task": MAX_STDOUT_BYTES_PER_TASK,
                "stderr_max_bytes_per_task": MAX_STDERR_BYTES_PER_TASK,
                "max_runtime_output_bytes": task_audit["max_runtime_output_bytes"],
                "campaign_output_budget_bytes": MAX_RUNTIME_OUTPUT_BYTES,
            },
            "timing_promotional": False,
            "timing_policy": (
                "NONPROMOTIONAL saturated rank-proxy recovery campaign; use paired "
                "XOR, muladd, and inactivation work only, with success/failure "
                "strata reported separately"),
            "reducer_policy": (
                "weak-K count and multiplicity, q histograms, repairs, introductions, "
                "and paired work for all four outcome quadrants"),
            "tasks_manifest_sha256": sha256_file(staging / "tasks_manifest.jsonl"),
            "independent_task_audit_sha256": sha256_file(
                staging / "independent_task_audit.json"),
        }
        if design["development_seed_overlap"]:
            die("sealed holdout seeds overlap development")
        if args.workers != (os.cpu_count() or args.workers):
            die("--workers must equal the available logical CPU count")
        write_once(staging / "design.json", canonical_json(design))
        paths = artifact_map(staging)
        receipts = {
            "schema": SCHEMA + ".prelaunch",
            "artifacts": {
                key: {"path": str(path.relative_to(staging)),
                      "sha256": sha256_file(path), "bytes": path.stat().st_size}
                for key, path in sorted(paths.items())
            },
        }
        write_once(staging / "prelaunch_receipts.json", canonical_json(receipts))

        executable = {
            frozen / CONTROLLER_NAME, frozen / SAMPLER_NAME,
            frozen / BASE_BINARY_NAME, frozen / CANDIDATE_BINARY_NAME,
        }
        for path in frozen.rglob("*"):
            if path.is_file():
                path.chmod(0o555 if path in executable else 0o444)
        for path in (
                staging / "tasks_manifest.jsonl",
                staging / "independent_task_audit.json",
                staging / "design.json", staging / "prelaunch_receipts.json"):
            path.chmod(0o444)
        os.replace(str(staging), str(root))
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    receipt_hash = sha256_file(root / "prelaunch_receipts.json")
    print(json.dumps({
        "status": "PREPARED_NOT_LAUNCHED", "root": str(root),
        "prelaunch_receipts_sha256": receipt_hash,
        "base_commit": BASE_COMMIT, "candidate_commit": candidate_commit,
        "base_binary_sha256": design["base_binary_sha256"],
        "candidate_binary_sha256": design["candidate_binary_sha256"],
        "base_replay_smoke_cases": len(smoke["cases"]),
        "tasks": len(tasks), "task_pairs": len(tasks) // 2,
        "cells": design["total_cells"], "workers": design["worker_count"],
        "preflight_command": (
            "{}/frozen/{} launch --root {} --expected-receipts-sha256 {} "
            "--preflight-only".format(
                root, CONTROLLER_NAME, root, receipt_hash)),
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
            value = json.loads(raw.decode("ascii"))
        except (UnicodeError, json.JSONDecodeError) as exc:
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
    receipt_path = root / "prelaunch_receipts.json"
    if not re.fullmatch(r"[0-9a-f]{64}", expected_receipts or ""):
        die("expected receipt hash must be 64 lowercase hex digits")
    if sha256_file(receipt_path) != expected_receipts:
        die("external prelaunch receipt hash mismatch")
    receipts = load_json(receipt_path)
    if receipts.get("schema") != SCHEMA + ".prelaunch":
        die("prelaunch schema mismatch")
    artifacts = receipts.get("artifacts")
    paths = artifact_map(root)
    if not isinstance(artifacts, dict) or set(artifacts) != set(paths):
        die("prelaunch artifact set mismatch")
    for key, expected_path in paths.items():
        record = artifacts[key]
        if not isinstance(record, dict) or \
                record.get("path") != str(expected_path.relative_to(root)) or \
                expected_path.is_symlink() or not expected_path.is_file() or \
                record.get("sha256") != sha256_file(expected_path) or \
                record.get("bytes") != expected_path.stat().st_size:
            die("prelaunch artifact mismatch: {}".format(key))

    design = load_json(root / "design.json")
    if design.get("schema") != SCHEMA + ".design" or \
            design.get("root") != str(root) or \
            design.get("base_commit") != BASE_COMMIT or \
            design.get("candidate_changed_paths") != list(EXPECTED_CANDIDATE_PATHS) or \
            design.get("seeds") != list(SEEDS) or \
            design.get("development_seeds") != list(DEVELOPMENT_SEEDS) or \
            design.get("development_seed_overlap") != [] or \
            design.get("schedules") != list(SCHEDULES) or \
            design.get("K_range") != [K_LO, K_HI] or \
            design.get("stage_order") != list(STAGES) or \
            design.get("thermal_cpu_busy_floor_pct") != str(CPU_BUSY_FLOOR) or \
            design.get("thermal_interval_seconds") != str(THERMAL_INTERVAL_SECONDS) or \
            design.get("thermal_min_gap_seconds") != str(THERMAL_MIN_GAP_SECONDS) or \
            design.get("thermal_max_gap_seconds") != str(THERMAL_MAX_GAP_SECONDS) or \
            design.get("thermal_ready_samples") != THERMAL_READY_SAMPLES or \
            design.get("max_attempts_per_stage") != MAX_ATTEMPTS_PER_STAGE or \
            design.get("task_timeout_seconds") != TASK_TIMEOUT_SECONDS or \
            design.get("max_stage_seconds") != MAX_STAGE_SECONDS or \
            design.get("max_thermal_bytes_per_segment") != \
                MAX_THERMAL_BYTES_PER_SEGMENT or \
            design.get("retry_policy") != \
                "stage-atomic non-selective retry after failed/interrupted segment" or \
            design.get("worker_count") != (os.cpu_count() or DEFAULT_WORKERS):
        die("design trust anchor mismatch")
    output_policy = design.get("output_policy")
    if not isinstance(output_policy, dict) or \
            output_policy.get("stdout_max_bytes_per_task") != MAX_STDOUT_BYTES_PER_TASK or \
            output_policy.get("stderr_max_bytes_per_task") != MAX_STDERR_BYTES_PER_TASK or \
            output_policy.get("campaign_output_budget_bytes") != MAX_RUNTIME_OUTPUT_BYTES:
        die("design output-bound policy mismatch")
    if design.get("base_to_candidate_patch_sha256") != sha256_file(
            root / "frozen" / "base_to_candidate.patch"):
        die("candidate patch hash binding mismatch")

    build_specs = (
        ("base", BASE_BINARY_NAME, BASE_COMMIT, design.get("base_tree_oid")),
        ("candidate", CANDIDATE_BINARY_NAME, design.get("candidate_commit"),
         design.get("candidate_tree_oid")),
    )
    for role, binary_name, commit, tree_oid in build_specs:
        receipt_name = "build_{}_receipt.json".format(role)
        build = load_json(root / "frozen" / receipt_name)
        binary = root / "frozen" / binary_name
        source = build.get("source")
        snapshot = build.get("candidate_repository_snapshot")
        evidence_dir = root / "frozen" / ("build_" + role)
        if not isinstance(source, dict) or not isinstance(snapshot, dict) or \
                build.get("schema") != SCHEMA + ".fresh_build" or \
                build.get("role") != role or build.get("fresh_build") is not True or \
                build.get("binary_name") != binary_name or \
                build.get("binary_sha256") != sha256_file(binary) or \
                source.get("commit") != commit or source.get("tree_oid") != tree_oid or \
                source.get("repository_head_at_build") != design.get("candidate_commit") or \
                source.get("repository_clean_tracked_status_sha256") != sha256_bytes(b"") or \
                snapshot.get("commit") != design.get("candidate_commit") or \
                snapshot.get("tracked_status_porcelain_sha256") != sha256_bytes(b""):
            die("{} fresh-build provenance binding mismatch".format(role))
        tree_path = evidence_dir / "source_tree.txt"
        if source.get("tree_manifest_sha256") != sha256_file(tree_path) or \
                source.get("tree_manifest_bytes") != tree_path.stat().st_size or \
                re.fullmatch(r"[0-9a-f]{64}",
                             str(source.get("source_archive_sha256", ""))) is None:
            die("{} source archive/tree receipt mismatch".format(role))
        if snapshot.get("base_commit") != BASE_COMMIT or \
                snapshot.get("changed_paths") != list(EXPECTED_CANDIDATE_PATHS) or \
                snapshot.get("base_to_candidate_patch_sha256") != \
                    design.get("base_to_candidate_patch_sha256") or \
                any(path not in ALLOWED_CANDIDATE_UNTRACKED
                    for path in snapshot.get("allowed_untracked_paths", [])):
            die("{} candidate repository snapshot mismatch".format(role))
        evidence = build.get("evidence")
        expected_names = {
            path.name for key, path in paths.items()
            if key.startswith("build_{}_".format(role)) and
            key not in ("build_{}_receipt".format(role),)
        }
        entries = list(evidence_dir.iterdir())
        if {path.name for path in entries} != expected_names or \
                any(path.is_symlink() or not path.is_file() for path in entries) or \
                not isinstance(evidence, dict) or set(evidence) != expected_names:
            die("{} build evidence inventory mismatch".format(role))
        for name in expected_names:
            path = evidence_dir / name
            record = evidence[name]
            if not isinstance(record, dict) or \
                    record.get("sha256") != sha256_file(path) or \
                    record.get("bytes") != path.stat().st_size:
                die("{} build evidence receipt mismatch: {}".format(role, name))
        design_hash_key = "{}_binary_sha256".format(role)
        design_receipt_key = "{}_build_receipt_sha256".format(role)
        if design.get(design_hash_key) != build.get("binary_sha256") or \
                design.get(design_receipt_key) != sha256_file(
                    root / "frozen" / receipt_name):
            die("{} design/build receipt mismatch".format(role))

    smoke = load_json(root / "frozen" / "binary_smoke.json")
    if smoke.get("schema") != SCHEMA + ".binary_smoke" or \
            smoke.get("bounded") is not True or \
            smoke.get("base_commit") != BASE_COMMIT or \
            smoke.get("base_binary_sha256") != design.get("base_binary_sha256") or \
            smoke.get("candidate_binary_sha256") != \
                design.get("candidate_binary_sha256") or \
            design.get("binary_smoke_sha256") != sha256_file(
                root / "frozen" / "binary_smoke.json") or \
            not isinstance(smoke.get("cases"), list) or \
            len(smoke["cases"]) != len(SEEDS) * len(SCHEDULES):
        die("bounded binary-smoke provenance binding mismatch")
    for case, record in enumerate(smoke["cases"]):
        seed_index = case // len(SCHEDULES)
        schedule = SCHEDULES[case % len(SCHEDULES)]
        Ks = [2, 257, 64_000]
        stem = "case{:02d}.seed{}.{}".format(case, seed_index, schedule)
        if not isinstance(record, dict) or record.get("case") != case or \
                record.get("seed_index") != seed_index or \
                record.get("seed") != SEEDS[seed_index] or \
                record.get("schedule") != schedule or record.get("Ks") != Ks or \
                record.get("base_replay_semantic_identical") is not True or \
                record.get("replay_excluded_fields") != \
                    list(NONDETERMINISTIC_REPLAY_FIELDS):
            die("bounded smoke case receipt mismatch")
        invocations = record.get("invocations")
        if not isinstance(invocations, dict) or set(invocations) != {
                "base", "candidate_base_replay", "degree_balanced"}:
            die("bounded smoke invocation set mismatch")
        bytes_by_name: Dict[str, bytes] = {}
        semantic_by_name: Dict[str, List[dict]] = {}
        for name, binary_name, arm in (
                ("base", BASE_BINARY_NAME, "base"),
                ("candidate_base_replay", CANDIDATE_BINARY_NAME, "base"),
                ("degree_balanced", CANDIDATE_BINARY_NAME, "degree_balanced")):
            stdout_path = root / "frozen" / "smoke" / (
                stem + "." + name + ".stdout")
            stderr_path = root / "frozen" / "smoke" / (
                stem + "." + name + ".stderr")
            invocation = invocations[name]
            task = {
                "arm": arm, "seed_index": seed_index, "schedule": schedule,
                "Ks": Ks, "stdout_max_bytes": MAX_STDOUT_BYTES_PER_TASK,
            }
            parsed = parse_benchmark_csv(stdout_path.read_bytes(), task)
            expected_argv = task_argv(
                str(root / "frozen" / binary_name), arm, seed_index, schedule, Ks)
            if not isinstance(invocation, dict) or \
                    invocation.get("argv") != [
                        "frozen/" + binary_name, *expected_argv[1:]] or \
                    invocation.get("stdout") != str(stdout_path.relative_to(root)) or \
                    invocation.get("stderr") != str(stderr_path.relative_to(root)) or \
                    invocation.get("stdout_sha256") != sha256_file(stdout_path) or \
                    invocation.get("stdout_bytes") != stdout_path.stat().st_size or \
                    invocation.get("stderr_sha256") != sha256_file(stderr_path) or \
                    stderr_path.stat().st_size or \
                    invocation.get("outcomes") != [classify(row) for row in parsed]:
                die("bounded smoke invocation receipt mismatch")
            bytes_by_name[name] = stdout_path.read_bytes()
            semantic_by_name[name] = [
                {
                    field: row[field] for field in CSV_FIELDS
                    if field not in NONDETERMINISTIC_REPLAY_FIELDS
                }
                for row in parsed
            ]
        if semantic_by_name["base"] != \
                semantic_by_name["candidate_base_replay"] or \
                record.get("base_replay_byte_identical") != (
                    bytes_by_name["base"] ==
                    bytes_by_name["candidate_base_replay"]) or \
                record.get("base_semantic_sha256") != sha256_bytes(canonical_json(
                    semantic_by_name["base"])) or \
                record.get("candidate_base_semantic_sha256") != sha256_bytes(
                    canonical_json(semantic_by_name["candidate_base_replay"])):
            die("bounded exact deterministic base replay differs")

    tasks = build_tasks(
        str(root / "frozen" / BASE_BINARY_NAME),
        str(root / "frozen" / CANDIDATE_BINARY_NAME),
    )
    audit = independent_task_audit(tasks)
    if (root / "tasks_manifest.jsonl").read_bytes() != json_lines(tasks) or \
            (root / "independent_task_audit.json").read_bytes() != canonical_json(audit) or \
            design.get("tasks_manifest_sha256") != sha256_file(
                root / "tasks_manifest.jsonl") or \
            design.get("independent_task_audit_sha256") != sha256_file(
                root / "independent_task_audit.json") or \
            design.get("task_count") != len(tasks) or \
            design.get("task_pairs") != len(tasks) // 2 or \
            design.get("total_cells") != audit["total_cells"] or \
            audit.get("max_attempts_per_stage") != MAX_ATTEMPTS_PER_STAGE or \
            audit.get("task_timeout_seconds") != TASK_TIMEOUT_SECONDS or \
            audit.get("max_stage_seconds") != MAX_STAGE_SECONDS or \
            audit.get("max_thermal_bytes_per_segment") != \
                MAX_THERMAL_BYTES_PER_SEGMENT or \
            output_policy.get("max_runtime_output_bytes") != \
                audit["max_runtime_output_bytes"]:
        die("regenerated paired all-K ledger mismatch")

    expected_outputs = {str(task["output_name"]): task for task in tasks}
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
        if any(not any(name.startswith(prefix + suffix + ".part.")
                       for prefix in expected_outputs) for name in partial):
            die("unexpected partial {} output(s)".format(directory))
        if directory in ("raw", "stderr"):
            for path in entries:
                if ".part." in path.name:
                    continue
                base_name = path.name[:-len(suffix)] if suffix else path.name
                limit_key = (
                    "stdout_max_bytes" if directory == "raw"
                    else "stderr_max_bytes")
                if path.stat().st_size > int(expected_outputs[base_name][limit_key]):
                    die("{} runtime artifact exceeds its sealed bound".format(directory))
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
                    paths[0].stat().st_size > int(task["stdout_max_bytes"]) or \
                    paths[1].stat().st_size > int(task["stderr_max_bytes"]) or \
                    paths[2].read_text(encoding="ascii") != "0\n" or \
                    receipt.get("schema") != SCHEMA + ".job" or \
                    receipt.get("job") != task["job"] or \
                    receipt.get("stage") != task["stage"] or \
                    receipt.get("output_name") != name or \
                    receipt.get("returncode") != 0 or \
                    receipt.get("stdout_sha256") != sha256_file(paths[0]) or \
                    receipt.get("stdout_bytes") != paths[0].stat().st_size or \
                    receipt.get("stderr_sha256") != sha256_file(paths[1]) or \
                    receipt.get("stderr_bytes") != paths[1].stat().st_size or \
                    receipt.get("exit_sha256") != sha256_file(paths[2]) or \
                    not isinstance(receipt.get("segment"), int) or \
                    isinstance(receipt.get("segment"), bool) or \
                    not isinstance(receipt.get("started_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("started_monotonic_s"), bool) or \
                    not isinstance(receipt.get("ended_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("ended_monotonic_s"), bool) or \
                    receipt["ended_monotonic_s"] < receipt["started_monotonic_s"]:
                die("unclean completed job {}".format(task["job"]))
            parse_benchmark_csv(paths[0].read_bytes(), task)
            completed.add(int(task["job"]))
    return completed


def process_state(pid: int) -> Optional[str]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
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
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
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


def process_alive(pid: int) -> bool:
    if pid <= 1:
        return False
    try:
        os.kill(pid, 0)
        return process_state(pid) != "Z"
    except ProcessLookupError:
        return False
    except PermissionError:
        # Root-owned samplers return EPERM to a non-root kill(2).  A plain
        # exception must never be interpreted as process death.
        result = subprocess.run(
            ["sudo", "-n", "kill", "-0", str(pid)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
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
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    return [token.decode("utf-8", "replace") for token in raw.split(b"\0") if token]


def other_samplers(excluded_pids: Iterable[int] = ()) -> List[dict]:
    excluded = {os.getpid(), *excluded_pids}
    found: List[dict] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) in excluded:
            continue
        pid = int(proc.name)
        tokens = process_tokens(pid)
        if tokens is not None and any(
                Path(token).name == SAMPLER_NAME for token in tokens) and \
                process_alive(pid):
            found.append({"pid": pid, "command": " ".join(tokens)})
    return sorted(found, key=lambda item: item["pid"])


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
    indices: Set[int] = set()
    patterns = (
        (root / "segments", re.compile(r"segment([0-9]{3})\.(?:intent|ready|final)\.json$")),
        (root / "thermal", re.compile(r"segment([0-9]{3})\.(?:csv|pid|stderr)$")),
        (root / "attempts", re.compile(r"segment([0-9]{3})$")),
    )
    for directory, pattern in patterns:
        for path in directory.iterdir():
            match = pattern.fullmatch(path.name)
            if match:
                indices.add(int(match.group(1)))
    ordered = sorted(indices)
    if ordered and ordered != list(range(ordered[-1] + 1)):
        die("campaign segment numbering has a gap")
    return ordered


def stage_attempt_history(root: Path) -> Dict[str, Dict[str, int]]:
    history = {
        stage: {"attempts": 0, "successful": 0} for stage in STAGES
    }
    for segment in segment_indices(root):
        intent_path = segment_path(root, segment, "intent")
        if intent_path.is_symlink() or not intent_path.is_file():
            die("segment history lacks a regular immutable intent")
        stage = load_json(intent_path).get("stage")
        if stage not in STAGES:
            die("segment intent has an unknown stage")
        history[str(stage)]["attempts"] += 1
        final_path = segment_path(root, segment, "final")
        if final_path.exists() or final_path.is_symlink():
            if final_path.is_symlink() or not final_path.is_file():
                die("segment history has an indirect final receipt")
            state = load_json(final_path).get("state")
            if state not in ("success", "failed", "interrupted"):
                die("segment history has an unknown terminal state")
            if state == "success":
                history[str(stage)]["successful"] += 1
    if any(value["successful"] > 1 for value in history.values()):
        die("a holdout stage has more than one successful attempt")
    return history


def stage_attempt_counts(root: Path) -> Dict[str, int]:
    return {
        stage: value["attempts"]
        for stage, value in stage_attempt_history(root).items()
    }


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
    if stage not in STAGES or any(task["stage"] != stage for task in selected):
        die("segment intent crosses campaign stages")
    if intent.get("jobs_sha256") != sha256_bytes(json_lines(raw)):
        die("segment intent job hash mismatch")
    return selected


def terminate_verified_process(pid: int, expected: Sequence[str]) -> str:
    if not process_alive(pid):
        return "already_exited"
    tokens = process_tokens(pid)
    if tokens != list(expected):
        die("refusing to terminate PID {} with changed command identity".format(pid))
    subprocess.run(
        ["sudo", "-n", "kill", "-TERM", str(pid)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline and process_alive(pid):
        time.sleep(0.05)
    if process_alive(pid):
        subprocess.run(
            ["sudo", "-n", "kill", "-KILL", str(pid)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
        deadline = time.monotonic() + 5.0
        while time.monotonic() < deadline and process_alive(pid):
            time.sleep(0.05)
        if process_alive(pid):
            die("verified owned process {} did not stop".format(pid))
        return "killed"
    return "terminated"


def terminate_owned_segment_processes(root: Path, segment: int, tasks: Sequence[dict]) -> List[dict]:
    actions: List[dict] = []
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    if thermal_pid.is_symlink():
        die("owned thermal PID file is a symlink")
    if thermal_pid.exists():
        try:
            pid = int(thermal_pid.read_text(encoding="ascii").strip())
        except (OSError, ValueError):
            die("malformed owned thermal PID file for segment {}".format(segment))
        ready_path = segment_path(root, segment, "ready")
        expected_start: Optional[int] = None
        if ready_path.is_file() and not ready_path.is_symlink():
            ready = load_json(ready_path)
            if ready.get("sampler_pid") == pid and \
                    isinstance(ready.get("sampler_start_ticks"), int) and \
                    not isinstance(ready.get("sampler_start_ticks"), bool):
                expected_start = int(ready["sampler_start_ticks"])
        if process_alive(pid) and expected_start is not None and \
                process_start_ticks(pid) != expected_start:
            actions.append({"kind": "thermal", "pid": pid,
                            "action": "stale_reused_pid"})
            try:
                thermal_pid.unlink()
            except FileNotFoundError:
                pass
            pid = -1
        tokens = process_tokens(pid) if process_alive(pid) else None
        if process_alive(pid) and tokens is None:
            die("cannot inspect the live owned thermal sampler")
        if tokens is not None:
            frozen_sampler = str(root / "frozen" / SAMPLER_NAME)
            if frozen_sampler not in tokens or str(thermal_csv) not in tokens:
                actions.append({"kind": "thermal", "pid": pid,
                                "action": "stale_changed_command"})
                try:
                    thermal_pid.unlink()
                except FileNotFoundError:
                    pass
                tokens = None
            if tokens is not None:
                actions.append({"kind": "thermal", "pid": pid,
                                "action": terminate_verified_process(pid, tokens)})
        if not process_alive(pid):
            try:
                thermal_pid.unlink()
            except FileNotFoundError:
                pass
    task_by_job = {int(task["job"]): task for task in tasks}
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    if attempt_dir.is_dir():
        for pid_path in sorted(attempt_dir.glob("job*.pid")):
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
                    not isinstance(start_ticks, int) or isinstance(start_ticks, bool) or \
                    start_ticks < 0:
                die("malformed benchmark PID identity")
            if not process_alive(pid):
                continue
            if process_start_ticks(pid) != start_ticks:
                actions.append({"kind": "benchmark", "job": job, "pid": pid,
                                "action": "stale_reused_pid"})
                continue
            tokens = process_tokens(pid)
            if tokens is None:
                die("cannot inspect a live owned benchmark PID")
            if tokens != task_by_job[job]["argv"]:
                actions.append({"kind": "benchmark", "job": job, "pid": pid,
                                "action": "stale_changed_command"})
                continue
            actions.append({"kind": "benchmark", "job": job, "pid": pid,
                            "action": terminate_verified_process(pid, tokens)})
    return actions


def rollback_segment_outputs(root: Path, selected: Sequence[dict]) -> List[int]:
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


def seal_attempt_manifest(root: Path, segment: int) -> Tuple[str, int]:
    directory = root / "attempts" / "segment{:03d}".format(segment)
    if not directory.exists():
        directory.mkdir()
    if directory.is_symlink() or not directory.is_dir():
        die("attempt evidence directory is indirect")
    for partial in directory.glob("*.part.*"):
        partial.unlink()
    files = sorted(path for path in directory.iterdir() if path.is_file())
    if any(path.is_symlink() for path in files):
        die("attempt evidence contains a symlink")
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
) -> Dict[str, object]:
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    attempts_sha256, attempt_file_count = seal_attempt_manifest(root, segment)
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
        "thermal_csv_sha256": sha256_file(thermal_csv) if thermal_csv.is_file() else None,
        "thermal_stderr": thermal_stderr.name,
        "thermal_stderr_sha256": (
            sha256_file(thermal_stderr) if thermal_stderr.is_file() else None),
        "attempt_manifest": "segment{:03d}.attempts.sha256".format(segment),
        "attempt_manifest_sha256": attempts_sha256,
        "attempt_file_count": attempt_file_count,
        "retry_policy": "entire stage segment; no successful survivor is retained on failure",
    }
    write_once(segment_path(root, segment, "final"), canonical_json(record))
    return record


def reconcile_incomplete_segments(root: Path, tasks: Sequence[dict]) -> List[dict]:
    reconciled: List[dict] = []
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
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                intent.get("segment") != segment:
            die("incomplete segment intent is malformed")
        selected = segment_jobs(intent, tasks)
        actions = terminate_owned_segment_processes(root, segment, selected)
        rolled_back = rollback_segment_outputs(root, selected)
        final = write_segment_final(
            root, segment, intent, "interrupted", failures=[],
            rolled_back=rolled_back, process_actions=actions,
            jobs_ended_monotonic_s=None, sampler_returncode=None,
        )
        reconciled.append(final)
    return reconciled


def read_thermal_rows(path: Path) -> List[Dict[str, str]]:
    try:
        if path.stat().st_size > MAX_THERMAL_BYTES_PER_SEGMENT:
            die("thermal stream exceeds its sealed segment bound")
    except OSError as exc:
        die("cannot stat thermal stream {}: {}".format(path, exc))
    try:
        with path.open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
                die("thermal header mismatch: {}".format(path))
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        die("cannot parse thermal stream {}: {}".format(path, exc))
    if any(set(row) != set(THERMAL_FIELDS) or any(value is None for value in row.values())
           for row in rows):
        die("thermal stream contains a malformed row")
    return rows


def sampler_identity(pid_file: Path, sampler_path: Path, csv_path: Path) -> Optional[int]:
    if pid_file.is_symlink():
        die("thermal sampler PID file is a symlink")
    if not pid_file.is_file():
        return None
    try:
        pid = int(pid_file.read_text(encoding="ascii").strip())
    except (OSError, ValueError):
        die("thermal sampler PID file is malformed")
    if not process_alive(pid):
        return None
    tokens = process_tokens(pid)
    if tokens is None or str(sampler_path) not in tokens or str(csv_path) not in tokens:
        die("thermal sampler PID is not bound to this segment")
    return pid


def wait_sampler_ready(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
) -> Tuple[int, List[Dict[str, str]]]:
    deadline = time.monotonic() + THERMAL_READY_TIMEOUT_SECONDS
    while time.monotonic() < deadline:
        if sampler.poll() is not None:
            die("thermal sampler exited before readiness")
        pid = sampler_identity(pid_file, sampler_path, csv_path)
        if pid is not None and csv_path.is_file():
            rows = read_thermal_rows(csv_path)
            if len(rows) >= THERMAL_READY_SAMPLES:
                last = decimal_field(rows[-1], "monotonic_s")
                age = Decimal(str(time.monotonic())) - last
                if 0 <= age <= THERMAL_MAX_GAP_SECONDS:
                    return pid, rows
        time.sleep(0.05)
    die("thermal sampler did not provide two live data rows")


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
        rows = read_thermal_rows(csv_path)
        if rows and decimal_field(rows[-1], "monotonic_s") >= target:
            return True
        time.sleep(0.05)
    return False


def stop_launched_sampler(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
) -> Tuple[Optional[int], List[dict]]:
    actions: List[dict] = []
    pid = sampler_identity(pid_file, sampler_path, csv_path)
    if pid is None and sampler.poll() is None:
        tokens = process_tokens(sampler.pid)
        if tokens is None or str(sampler_path) not in tokens or str(csv_path) not in tokens:
            die("cannot prove ownership of launched thermal sampler")
        pid = sampler.pid
    if pid is not None:
        tokens = process_tokens(pid) if process_alive(pid) else None
        if tokens is not None:
            actions.append({"kind": "thermal", "pid": pid,
                            "action": terminate_verified_process(pid, tokens)})
    try:
        returncode = sampler.wait(timeout=15.0)
    except subprocess.TimeoutExpired:
        tokens = process_tokens(sampler.pid)
        if tokens is None:
            die("launched sampler wrapper did not exit and cannot be inspected")
        actions.append({"kind": "thermal_wrapper", "pid": sampler.pid,
                        "action": terminate_verified_process(sampler.pid, tokens)})
        returncode = sampler.wait(timeout=5.0)
    try:
        pid_file.unlink()
    except FileNotFoundError:
        pass
    return returncode, actions


def run_segment(
    root: Path,
    design: Mapping[str, object],
    tasks: Sequence[dict],
    stage: str,
    previously_complete: int,
    resume: bool,
) -> Dict[str, object]:
    history = stage_attempt_history(root)
    if history.get(stage, {}).get("successful", 0):
        die("stage {} already has a successful attempt".format(stage))
    if history.get(stage, {}).get("attempts", 0) >= MAX_ATTEMPTS_PER_STAGE:
        die("stage {} exhausted its sealed attempt budget".format(stage))
    indices = segment_indices(root)
    segment = indices[-1] + 1 if indices else 0
    jobs = sorted(int(task["job"]) for task in tasks)
    intent = {
        "schema": SCHEMA + ".segment_intent", "segment": segment,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "created_monotonic_s": time.monotonic(),
        "stage": stage, "jobs": jobs,
        "jobs_sha256": sha256_bytes(json_lines(jobs)),
        "workers": design["worker_count"], "previously_complete": previously_complete,
        "resume": bool(resume),
        "retry_policy": "stage-atomic non-selective retry",
    }
    write_once(segment_path(root, segment, "intent"), canonical_json(intent))
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    attempt_dir.mkdir()
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    sampler_path = root / "frozen" / SAMPLER_NAME
    sampler_error = thermal_stderr.open("xb")
    sampler_command = [
        "sudo", "-n", sys.executable, str(sampler_path),
        "--csv", str(thermal_csv), "--pid-file", str(thermal_pid),
        "--interval", str(THERMAL_INTERVAL_SECONDS),
        "--dimm-attempts", "5", "--dimm-retry-delay", "0.01",
    ]
    try:
        sampler = subprocess.Popen(
            sampler_command, stdout=subprocess.DEVNULL, stderr=sampler_error,
        )
    except OSError as exc:
        sampler_error.close()
        return write_segment_final(
            root, segment, intent, "failed",
            failures=[{"status": "sampler_spawn_error", "detail": str(exc)}],
            rolled_back=[], process_actions=[], jobs_ended_monotonic_s=None,
            sampler_returncode=None,
        )
    abort = threading.Event()
    active: Dict[int, subprocess.Popen] = {}
    active_lock = threading.Lock()
    failures: List[dict] = []
    process_actions: List[dict] = []
    jobs_ended: Optional[float] = None
    sampler_returncode: Optional[int] = None
    ready_written = False

    def terminate_active() -> None:
        with active_lock:
            snapshot = list(active.items())
        for job, process in snapshot:
            if process.poll() is not None:
                continue
            tokens = process_tokens(process.pid)
            if tokens is None:
                try:
                    process.terminate()
                    action = "terminated_via_owned_popen"
                except ProcessLookupError:
                    action = "already_exited"
                process_actions.append({
                    "kind": "benchmark", "job": job, "pid": process.pid,
                    "action": action,
                })
                continue
            process_actions.append({
                "kind": "benchmark", "job": job, "pid": process.pid,
                "action": terminate_verified_process(process.pid, tokens),
            })

    def run_one(task: dict) -> dict:
        job = int(task["job"])
        if abort.is_set():
            return {"job": job, "status": "cancelled_before_start"}
        started_monotonic = time.monotonic()
        try:
            process = subprocess.Popen(
                task["argv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                start_new_session=True,
            )
        except OSError as exc:
            return {"job": job, "status": "spawn_error", "detail": str(exc)}
        start_ticks = process_start_ticks(process.pid)
        if start_ticks is None:
            process.terminate()
            process.wait()
            return {"job": job, "status": "pid_identity_error"}
        atomic_result(
            attempt_dir / "job{:05d}.pid".format(job),
            canonical_json({"pid": process.pid, "start_ticks": start_ticks}),
        )
        with active_lock:
            active[job] = process
        if abort.is_set() and process.poll() is None:
            tokens = process_tokens(process.pid)
            if tokens is not None:
                terminate_verified_process(process.pid, tokens)
            else:
                process.terminate()
        timed_out = False
        try:
            stdout, stderr = process.communicate(
                timeout=float(task["timeout_seconds"]))
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                process.terminate()
                process_actions.append({
                    "kind": "benchmark", "job": job, "pid": process.pid,
                    "action": "terminated_after_task_timeout",
                })
            except ProcessLookupError:
                pass
            try:
                stdout, stderr = process.communicate(timeout=10.0)
            except subprocess.TimeoutExpired:
                process.kill()
                process_actions.append({
                    "kind": "benchmark", "job": job, "pid": process.pid,
                    "action": "killed_after_task_timeout",
                })
                stdout, stderr = process.communicate(timeout=5.0)
        with active_lock:
            active.pop(job, None)
        ended_monotonic = time.monotonic()
        exit_bytes = (str(process.returncode) + "\n").encode("ascii")
        for suffix, data in (("stdout", stdout), ("stderr", stderr), ("exit", exit_bytes)):
            atomic_result(attempt_dir / "job{:05d}.{}".format(job, suffix), data)
        if len(stdout) > int(task["stdout_max_bytes"]):
            return {"job": job, "status": "stdout_bound_exceeded",
                    "bytes": len(stdout), "bound": task["stdout_max_bytes"]}
        if len(stderr) > int(task["stderr_max_bytes"]):
            return {"job": job, "status": "stderr_bound_exceeded",
                    "bytes": len(stderr), "bound": task["stderr_max_bytes"]}
        if timed_out:
            return {"job": job, "status": "task_timeout",
                    "timeout_seconds": task["timeout_seconds"]}
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
            "stdout_bytes": len(stdout),
            "stderr_sha256": sha256_bytes(stderr),
            "stderr_bytes": len(stderr),
            "exit_sha256": sha256_bytes(exit_bytes),
        }
        write_once(job_receipt_path(root, job), canonical_json(receipt))
        return {"job": job, "status": "success"}

    state = "failed"
    caught: Optional[BaseException] = None
    try:
        sampled_pid, ready_rows = wait_sampler_ready(
            sampler, thermal_pid, sampler_path, thermal_csv,
        )
        sampled_start_ticks = process_start_ticks(sampled_pid)
        if sampled_start_ticks is None:
            die("cannot bind thermal sampler to its process start time")
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
        future_results: List[dict] = []
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=int(design["worker_count"])) as pool:
            pending = {pool.submit(run_one, task) for task in tasks}
            while pending:
                if time.monotonic() - jobs_started > MAX_STAGE_SECONDS and \
                        not abort.is_set():
                    failures.append({
                        "status": "stage_timeout",
                        "timeout_seconds": MAX_STAGE_SECONDS,
                    })
                    abort.set()
                    terminate_active()
                if sampler.poll() is not None and not abort.is_set():
                    failures.append({"status": "sampler_exited", "returncode": sampler.returncode})
                    abort.set()
                    terminate_active()
                done, pending = concurrent.futures.wait(
                    pending, timeout=0.2,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    if future.cancelled():
                        result = {"status": "future_cancelled"}
                    else:
                        try:
                            result = future.result()
                        except BaseException as exc:
                            result = {"status": "worker_exception", "detail": repr(exc)}
                    future_results.append(result)
                    if result.get("status") != "success" and not abort.is_set():
                        failures.append(result)
                        abort.set()
                        terminate_active()
                if abort.is_set():
                    terminate_active()
                    for future in pending:
                        future.cancel()
            for result in future_results:
                if result.get("status") != "success" and result not in failures:
                    failures.append(result)
        jobs_ended = time.monotonic()
        if not failures and not wait_thermal_coverage(sampler, thermal_csv, jobs_ended):
            failures.append({"status": "thermal_end_coverage_failed"})
        if not failures:
            state = "success"
    except BaseException as exc:
        caught = exc
        abort.set()
        terminate_active()
        failures.append({"status": "launcher_exception", "detail": repr(exc)})
        jobs_ended = time.monotonic() if ready_written else None
    finally:
        try:
            sampler_returncode, sampler_actions = stop_launched_sampler(
                sampler, thermal_pid, sampler_path, thermal_csv,
            )
            process_actions.extend(sampler_actions)
        except BaseException as exc:
            if caught is None:
                caught = exc
            failures.append({"status": "sampler_stop_error", "detail": repr(exc)})
            state = "failed"
        sampler_error.close()
    if sampler_returncode != 0:
        failures.append({"status": "sampler_returncode", "returncode": sampler_returncode})
        state = "failed"
    if not thermal_stderr.is_file() or thermal_stderr.stat().st_size:
        failures.append({"status": "sampler_stderr_nonempty_or_missing"})
        state = "failed"
    rolled_back: List[int] = []
    if state != "success":
        rolled_back = rollback_segment_outputs(root, tasks)
    final = write_segment_final(
        root, segment, intent, state, failures=failures,
        rolled_back=rolled_back, process_actions=process_actions,
        jobs_ended_monotonic_s=jobs_ended,
        sampler_returncode=sampler_returncode,
    )
    if caught is not None:
        raise caught
    return final


def command_launch(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(
        root, args.expected_receipts_sha256, runtime_tolerant=True,
    )
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
    if args.preflight_only:
        samplers = other_samplers()
        complete_count: Optional[int] = None
        if not needs_reconciliation:
            complete_count = len(completed_jobs(root, tasks))
        print(json.dumps({
            "status": (
                "RESUME_RECONCILIATION_REQUIRED" if needs_reconciliation
                else "BLOCKED_CONCURRENT_SAMPLER" if samplers
                else "READY_NOT_LAUNCHED"),
            "root": str(root), "tasks": len(tasks), "complete": complete_count,
            "remaining": None if complete_count is None else len(tasks) - complete_count,
            "workers": design["worker_count"], "requested_stage": args.stage,
            "stage_order": list(STAGES),
            "stage_attempt_counts": stage_attempt_counts(root),
            "max_attempts_per_stage": MAX_ATTEMPTS_PER_STAGE,
            "single_sampler": not bool(samplers), "other_samplers": samplers,
            "incomplete_segments": incomplete,
            "orphan_receipt_partials": [
                str(path.relative_to(root)) for path in receipt_partials],
        }, indent=2, sort_keys=True))
        return
    if needs_reconciliation and not args.resume:
        die("interrupted segment(s) require --resume reconciliation")
    if indices and not args.resume:
        die("existing campaign segment history requires --resume")
    reconciled: List[dict] = []
    if args.resume:
        for path in receipt_partials:
            path.unlink()
        reconciled = reconcile_incomplete_segments(root, tasks)
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    complete = completed_jobs(root, tasks)
    stage_all = {
        stage: [task for task in tasks if task["stage"] == stage]
        for stage in STAGES
    }
    by_stage = {
        stage: [task for task in stage_all[stage]
                if int(task["job"]) not in complete]
        for stage in STAGES
    }
    for stage in STAGES:
        done = len(stage_all[stage]) - len(by_stage[stage])
        if done not in (0, len(stage_all[stage])):
            die("stage {} is selectively complete; reconciliation is unsafe".format(stage))
    samplers = other_samplers()
    if samplers:
        die("refusing concurrent SPD/I2C reader(s): {}".format(
            json.dumps(samplers, sort_keys=True)))
    if complete and not args.resume:
        die("outputs already exist; inspect then pass --resume")
    if not any(by_stage.values()):
        die("campaign ledger is already complete")
    requested_stages = STAGES if args.stage == "all" else (args.stage,)
    first_requested = STAGES.index(requested_stages[0])
    if any(by_stage[stage] for stage in STAGES[:first_requested]):
        die("requested stage requires every earlier holdout stratum")
    if not any(by_stage[stage] for stage in requested_stages):
        die("requested campaign stage is already complete")
    launched: List[dict] = []
    for stage in requested_stages:
        stage_tasks = by_stage[stage]
        if not stage_tasks:
            continue
        if any(by_stage[prior] for prior in STAGES[:STAGES.index(stage)]):
            die("cannot cross an incomplete earlier holdout stratum")
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
            raise SystemExit(1)
        complete.update(int(task["job"]) for task in stage_tasks)
        by_stage[stage] = []
    print(json.dumps({
        "status": "LAUNCH_OK", "reconciled_segments": len(reconciled),
        "launched_segments": [record["segment"] for record in launched],
        "complete": len(complete), "remaining": len(tasks) - len(complete),
    }, indent=2, sort_keys=True))


def decimal_field(row: Mapping[str, str], key: str) -> Decimal:
    try:
        value = Decimal(row[key])
    except (KeyError, InvalidOperation, TypeError) as exc:
        die("invalid decimal field {}: {}".format(key, exc))
    if not value.is_finite():
        die("nonfinite decimal field {}".format(key))
    return value


def validate_record(row: Mapping[str, str]) -> Dict[str, object]:
    try:
        success, rank_fail, error = (
            int(row[key]) for key in ("success", "rank_fail", "error"))
        inactivated = int(row["inact_max"])
        binary_deficit = int(row["binary_def_max"])
        heavy_gain = int(row["heavy_gain_min"])
        heavy_shortfall = int(row["heavy_shortfall"])
        seed_attempt = int(row["seed_attempt"])
        first_rank_fail = int(row["first_rank_fail"])
    except (KeyError, TypeError, ValueError) as exc:
        die("invalid benchmark integer field: {}".format(exc))
    if any(value not in (0, 1) for value in (success, rank_fail, error)) or \
            success + rank_fail + error != 1 or error:
        die("invalid benchmark outcome")
    if min(inactivated, binary_deficit, heavy_gain, heavy_shortfall, seed_attempt) < 0:
        die("negative benchmark counter")
    if heavy_shortfall not in (0, 1):
        die("invalid heavy-shortfall receipt")
    decimal_metrics = {
        key: decimal_field(row, key)
        for key in (
            "inact_mu", "binary_def_mu", "heavy_gain_mu",
            "block_xors_mu", "block_muladds_mu", "solve_ms_mu",
            "build_ms_mu", "peel_ms_mu", "project_ms_mu",
            "residual_ms_mu", "backsub_ms_mu",
        )
    }
    if any(value < 0 for value in decimal_metrics.values()):
        die("negative benchmark metric")
    if decimal_metrics["inact_mu"] != inactivated or \
            decimal_metrics["binary_def_mu"] != binary_deficit or \
            decimal_metrics["heavy_gain_mu"] != heavy_gain:
        die("one-trial mean/max/min receipt mismatch")
    expected_shortfall = int(
        rank_fail == 1 and binary_deficit <= 12 and heavy_gain < binary_deficit)
    if heavy_shortfall != expected_shortfall:
        die("heavy-shortfall receipt is inconsistent")
    outcome = (
        "success" if success else
        "q>H" if binary_deficit > 12 else
        "field_shortfall" if heavy_shortfall else
        "other_rank_fail"
    )
    if decimal_field(row, "fail_rate") != (
            Decimal(0) if success else Decimal(1)):
        die("one-trial fail-rate receipt mismatch")
    if row["binary_def_hist"] != "{}:1".format(binary_deficit) or \
            row["heavy_gain_hist"] != "{}:1".format(heavy_gain) or \
            first_rank_fail != (-1 if success else 0) or \
            row["failure_trials"] != ("" if success else "0"):
        die("one-trial histogram/failure receipt mismatch")
    for key in (
            "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
            "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
            "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu"):
        if decimal_field(row, key) != 0:
            die("production separate-bucket receipt {} is nonzero".format(key))
    return {
        "outcome": outcome, "success": success, "rank_fail": rank_fail,
        "binary_deficit": binary_deficit, "q": binary_deficit,
        "heavy_gain": heavy_gain, "heavy_shortfall": heavy_shortfall,
        "seed_attempt": seed_attempt, "inactivated": Decimal(inactivated),
        "xors": decimal_metrics["block_xors_mu"],
        "muladds": decimal_metrics["block_muladds_mu"],
    }


def classify(row: Mapping[str, str]) -> str:
    return str(validate_record(row)["outcome"])


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
        if match is None:
            die("segment attempt inventory changed during validation")
        job = int(match.group(1))
        if job not in allowed_jobs:
            die("segment attempt inventory names a job outside its intent")
        suffix = path.suffix[1:]
        groups[job].add(suffix)
        if suffix == "pid":
            identity = load_json(path)
            value = identity.get("pid")
            start_ticks = identity.get("start_ticks")
            if not isinstance(value, int) or isinstance(value, bool) or value <= 1 or \
                    not isinstance(start_ticks, int) or isinstance(start_ticks, bool) or \
                    start_ticks < 0:
                die("segment attempt PID receipt is out of range")
    complete_suffixes = {"pid", "stdout", "stderr", "exit"}
    for job, suffixes in groups.items():
        if suffixes != {"pid"} and suffixes != complete_suffixes:
            die("segment attempt contains a partial completed artifact set")
        if suffixes == complete_suffixes:
            exit_path = attempt_dir / "job{:05d}.exit".format(job)
            try:
                int(exit_path.read_text(encoding="ascii").strip())
            except (OSError, ValueError):
                die("segment attempt exit receipt is malformed")
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
        try:
            utc = datetime.fromisoformat(row["utc"].replace("Z", "+00:00"))
        except (KeyError, ValueError) as exc:
            die("invalid thermal UTC timestamp: {}".format(exc))
        if utc.tzinfo is None or (utc_previous is not None and utc <= utc_previous):
            die("thermal UTC timestamps are not strictly increasing")
        utc_previous = utc
        current = decimal_field(row, "monotonic_s")
        if monotonic and current <= monotonic[-1]:
            die("thermal monotonic timestamps are not strictly increasing")
        if monotonic and not (
                THERMAL_MIN_GAP_SECONDS <= current - monotonic[-1] <=
                THERMAL_MAX_GAP_SECONDS):
            die("thermal cadence gap is outside policy in segment {}".format(segment))
        monotonic.append(current)
        if row["cpu_busy_pct"]:
            value = decimal_field(row, "cpu_busy_pct")
            if value < 0 or value > 100:
                die("thermal CPU busy percentage is out of range")
            busy.append(value)
        if not row["cpu_tctl_c"]:
            die("thermal row lacks CPU Tctl")
        cpu.append(decimal_field(row, "cpu_tctl_c"))
        for key in THERMAL_FIELDS[5:13]:
            if not row[key]:
                die("thermal row lacks DIMM field {}".format(key))
            dimm_by_field[key].append(decimal_field(row, key))
        try:
            dimm_errors = int(row["dimm_read_errors"])
            edac_ce = int(row["edac_ce"])
            edac_ue = int(row["edac_ue"])
        except (KeyError, ValueError) as exc:
            die("invalid thermal hardware counter: {}".format(exc))
        if dimm_errors != 0 or edac_ce != 0 or edac_ue != 0:
            die("thermal hardware error receipt is nonzero")
    return {
        "samples": len(rows), "monotonic": monotonic, "busy": busy,
        "cpu": cpu, "dimm_by_field": dimm_by_field,
    }


def validate_thermal(root: Path, design: Mapping[str, object], tasks: Sequence[dict]) -> Dict[str, object]:
    indices = segment_indices(root)
    if not indices:
        die("campaign has no launch segments")
    expected_segment_files: Set[str] = set()
    successful_jobs: Set[int] = set()
    successful_segments = failed_segments = interrupted_segments = 0
    samples = 0
    all_busy: List[Decimal] = []
    all_cpu: List[Decimal] = []
    dimm_maxima: Dict[str, Decimal] = {}
    last_success_stage_index = -1
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
        if stage not in STAGES:
            die("segment intent has an unknown holdout stratum")
        stage_index = STAGES.index(stage)
        expected_jobs = [
            int(task["job"]) for task in tasks
            if task["stage"] == stage and int(task["job"]) not in successful_jobs
        ]
        earlier_remaining = any(
            STAGES.index(str(task["stage"])) < stage_index and
            int(task["job"]) not in successful_jobs for task in tasks
        )
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                intent.get("segment") != segment or \
                intent.get("workers") != design.get("worker_count") or \
                intent.get("previously_complete") != len(successful_jobs) or \
                jobs != expected_jobs or earlier_remaining or \
                intent.get("retry_policy") != "stage-atomic non-selective retry" or \
                final.get("schema") != SCHEMA + ".segment_final" or \
                final.get("segment") != segment or \
                final.get("intent_sha256") != sha256_file(intent_path) or \
                final.get("jobs") != jobs or final.get("jobs_sha256") != intent.get("jobs_sha256"):
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
                    ready.get("segment") != segment or \
                    ready.get("intent_sha256") != sha256_file(intent_path) or \
                    final.get("ready_sha256") != sha256_file(ready_path) or \
                    not isinstance(samples_at_ready, int) or isinstance(samples_at_ready, bool) or \
                    samples_at_ready < THERMAL_READY_SAMPLES or \
                    not isinstance(ready.get("sampler_pid"), int) or \
                    isinstance(ready.get("sampler_pid"), bool) or ready["sampler_pid"] <= 1 or \
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
        thermal: Optional[Dict[str, object]] = None
        rows: List[Dict[str, str]] = []
        if csv_exists:
            if final.get("thermal_csv") != csv_path.name or \
                    final.get("thermal_csv_sha256") != sha256_file(csv_path) or \
                    (stderr_exists and final.get("thermal_stderr") != stderr_path.name):
                die("segment thermal hash binding mismatch")
            if ready is not None and stderr_path.stat().st_size:
                die("thermal sampler stderr is nonempty")
            rows = read_thermal_rows(csv_path)
            if ready is not None:
                thermal = validate_thermal_rows(rows, segment=segment)
                count = int(ready["samples_at_ready"])
                if count > len(rows) or \
                        ready.get("last_sample_monotonic_s") != rows[count - 1]["monotonic_s"] or \
                        Decimal(str(ready["jobs_started_monotonic_s"])) < \
                            decimal_field(rows[count - 1], "monotonic_s") or \
                        Decimal(str(ready["jobs_started_monotonic_s"])) - \
                            decimal_field(rows[count - 1], "monotonic_s") > THERMAL_MAX_GAP_SECONDS:
                    die("segment readiness is not bound to its live thermal prefix")
                samples += int(thermal["samples"])
                all_cpu.extend(thermal["cpu"])
                for key, values in thermal["dimm_by_field"].items():
                    dimm_maxima[key] = max(dimm_maxima.get(key, values[0]), max(values))
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
            monotonic = thermal["monotonic"]
            if end < start or monotonic[0] > start or monotonic[-1] < end:
                die("successful thermal segment lacks full launch start/end coverage")
            # Include the first sample at/after the end boundary and the sample
            # immediately before start, since CPU busy is interval-derived.
            first_after_end = next(index for index, value in enumerate(monotonic) if value >= end)
            first_after_start = next(index for index, value in enumerate(monotonic) if value >= start)
            coverage_set = set(range(max(0, first_after_start - 1), first_after_end + 1))
            coverage_busy = [
                decimal_field(rows[index], "cpu_busy_pct")
                for index in sorted(coverage_set)
                if rows[index]["cpu_busy_pct"]
            ]
            if not coverage_busy or sum(coverage_busy) / len(coverage_busy) < CPU_BUSY_FLOOR:
                die("successful segment CPU busy mean is below the sealed floor")
            all_busy.extend(coverage_busy)
            if stage_index != last_success_stage_index + 1:
                die("successful holdout strata are not in the sealed order")
            last_success_stage_index = stage_index
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
        if stderr_exists and final.get("thermal_stderr") != stderr_path.name or \
                (stderr_exists and final.get("thermal_stderr_sha256") != sha256_file(stderr_path)):
            die("segment thermal stderr hash binding mismatch")
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
    return {
        "segments": len(indices), "successful_segments": successful_segments,
        "failed_segments": failed_segments,
        "interrupted_segments": interrupted_segments, "samples": samples,
        "successful_busy_mean": str(sum(all_busy) / len(all_busy)),
        "cpu_max_c": str(max(all_cpu)) if all_cpu else None,
        "dimm_max_c_by_field": {
            key: str(dimm_maxima[key]) for key in sorted(dimm_maxima)
        },
    }


def command_reduce(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    if len(completed_jobs(root, tasks)) != len(tasks):
        die("campaign is incomplete")
    thermal = validate_thermal(root, design, tasks)
    records: Dict[str, Dict[Tuple[int, int, str], dict]] = {
        arm: {} for arm in ARMS
    }
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
        name = str(task["output_name"])
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        for path in (raw_path, stderr_path, exit_path):
            manifest_lines.extend("{}  {}\n".format(
                sha256_file(path), path.relative_to(root)).encode("ascii"))
        rows = parse_benchmark_csv(raw_path.read_bytes(), task)
        arm = str(task["arm"])
        for row in rows:
            key = (int(row["N"]), int(task["seed_index"]), str(task["schedule"]))
            if key in records[arm]:
                die("duplicate Cartesian cell {}/{}".format(arm, key))
            records[arm][key] = validate_record(row)

    expected_keys = {
        (K, seed_index, schedule)
        for K in range(K_LO, K_HI + 1)
        for seed_index in range(len(SEEDS))
        for schedule in SCHEDULES
    }
    if any(set(records[arm]) != expected_keys for arm in ARMS):
        die("result cube does not exactly cover both all-K arms")

    metric_names = ("xors", "muladds", "inactivated")

    def work_summary(values: Sequence[Mapping[str, object]]) -> Dict[str, object]:
        result: Dict[str, object] = {"cells": len(values)}
        for metric in metric_names:
            total = sum((Decimal(value[metric]) for value in values), Decimal(0))
            result[metric + "_sum"] = str(total)
            result[metric + "_mean"] = (
                str(total / len(values)) if values else None)
        return result

    def arm_summary(arm: str) -> Dict[str, object]:
        values = records[arm]
        failure_records: List[dict] = []
        weak_multiplicity: Dict[int, int] = defaultdict(int)
        causes: Dict[str, int] = defaultdict(int)
        by_seed: Dict[int, int] = defaultdict(int)
        by_schedule: Dict[str, int] = defaultdict(int)
        q_all: Dict[int, int] = defaultdict(int)
        q_failure: Dict[int, int] = defaultdict(int)
        seed_attempt_hist: Dict[int, int] = defaultdict(int)
        success_values: List[dict] = []
        failure_values: List[dict] = []
        for key in sorted(values):
            record = values[key]
            K, seed_index, schedule = key
            outcome = str(record["outcome"])
            causes[outcome] += 1
            q = int(record["q"])
            q_all[q] += 1
            seed_attempt_hist[int(record["seed_attempt"])] += 1
            if outcome == "success":
                success_values.append(record)
            else:
                failure_values.append(record)
                weak_multiplicity[K] += 1
                by_seed[seed_index] += 1
                by_schedule[schedule] += 1
                q_failure[q] += 1
                failure_records.append({
                    **cell_record(key), "cause": outcome, "q": q,
                    "heavy_gain": int(record["heavy_gain"]),
                    "seed_attempt": int(record["seed_attempt"]),
                })
        multiplicity_hist: Dict[int, int] = defaultdict(int)
        for multiplicity in weak_multiplicity.values():
            multiplicity_hist[multiplicity] += 1
        return {
            "cells": len(values), "outcomes": dict(sorted(causes.items())),
            "failures": len(failure_values), "weak_K": len(weak_multiplicity),
            "weak_K_multiplicity": {
                str(K): weak_multiplicity[K] for K in sorted(weak_multiplicity)
            },
            "weak_K_multiplicity_histogram": {
                str(value): multiplicity_hist[value]
                for value in sorted(multiplicity_hist)
            },
            "failures_by_seed": {
                str(seed): by_seed[seed] for seed in sorted(by_seed)
            },
            "failures_by_schedule": {
                schedule: by_schedule[schedule] for schedule in sorted(by_schedule)
            },
            "q_histogram_all": {str(q): q_all[q] for q in sorted(q_all)},
            "q_histogram_failures": {
                str(q): q_failure[q] for q in sorted(q_failure)
            },
            "seed_attempt_histogram": {
                str(value): seed_attempt_hist[value]
                for value in sorted(seed_attempt_hist)
            },
            "success_work": work_summary(success_values),
            "failure_work": work_summary(failure_values),
            "failure_records": failure_records,
        }

    quadrants: Dict[str, List[Tuple[Tuple[int, int, str], dict, dict]]] = {
        "both_success": [], "repair": [], "introduction": [], "both_fail": [],
    }
    q_transitions: Dict[Tuple[int, int], int] = defaultdict(int)
    receipt_deltas = {
        field: {"cells": 0, "candidate_minus_base_sum": 0, "max_abs": 0}
        for field in ("seed_attempt", "binary_deficit")
    }
    inactivation_delta_cells = 0
    inactivation_delta_sum = Decimal(0)
    inactivation_delta_max = Decimal(0)
    for key in sorted(expected_keys):
        base = records["base"][key]
        candidate = records["degree_balanced"][key]
        base_success = base["outcome"] == "success"
        candidate_success = candidate["outcome"] == "success"
        quadrant = (
            "both_success" if base_success and candidate_success else
            "repair" if not base_success and candidate_success else
            "introduction" if base_success and not candidate_success else
            "both_fail"
        )
        quadrants[quadrant].append((key, base, candidate))
        q_transitions[(int(base["q"]), int(candidate["q"]))] += 1
        for field, stats in receipt_deltas.items():
            delta = int(candidate[field]) - int(base[field])
            if delta:
                stats["cells"] += 1
                stats["candidate_minus_base_sum"] += delta
                stats["max_abs"] = max(stats["max_abs"], abs(delta))
        delta_inact = Decimal(candidate["inactivated"]) - Decimal(base["inactivated"])
        if delta_inact:
            inactivation_delta_cells += 1
            inactivation_delta_sum += delta_inact
            inactivation_delta_max = max(inactivation_delta_max, abs(delta_inact))

    paired_work: Dict[str, object] = {}
    for quadrant, values in quadrants.items():
        base_values = [value[1] for value in values]
        candidate_values = [value[2] for value in values]
        paired_work[quadrant] = {
            "cells": len(values),
            "base": work_summary(base_values),
            "degree_balanced": work_summary(candidate_values),
            "candidate_minus_base_sums": {
                metric: str(sum(
                    (Decimal(candidate[metric]) - Decimal(base[metric])
                     for _, base, candidate in values), Decimal(0)))
                for metric in metric_names
            },
        }
    common = paired_work["both_success"]
    for metric in metric_names:
        base_total = Decimal(common["base"][metric + "_sum"])
        candidate_total = Decimal(
            common["degree_balanced"][metric + "_sum"])
        common[metric + "_ratio"] = (
            str(candidate_total / base_total) if base_total else None)

    def key_receipts(
        values: Sequence[Tuple[Tuple[int, int, str], dict, dict]]
    ) -> List[dict]:
        return [
            {
                **cell_record(key), "base_cause": base["outcome"],
                "candidate_cause": candidate["outcome"],
                "base_q": int(base["q"]), "candidate_q": int(candidate["q"]),
            }
            for key, base, candidate in values
        ]

    comparison = {
        "repairs": len(quadrants["repair"]),
        "introductions": len(quadrants["introduction"]),
        "net_failure_change": (
            len(quadrants["introduction"]) - len(quadrants["repair"])),
        "both_success": len(quadrants["both_success"]),
        "both_fail": len(quadrants["both_fail"]),
        "repair_keys": key_receipts(quadrants["repair"]),
        "introduction_keys": key_receipts(quadrants["introduction"]),
        "both_fail_keys": key_receipts(quadrants["both_fail"]),
        "q_transition_histogram": {
            "{}->{}".format(before, after): q_transitions[(before, after)]
            for before, after in sorted(q_transitions)
        },
        "paired_work_by_outcome_quadrant": paired_work,
        "paired_integer_receipt_deltas": receipt_deltas,
        "paired_inactivation_receipt_delta": {
            "cells": inactivation_delta_cells,
            "candidate_minus_base_sum": str(inactivation_delta_sum),
            "max_abs": str(inactivation_delta_max),
        },
    }
    data_manifest = bytes(manifest_lines)
    write_once(root / "data_manifest.sha256", data_manifest)
    summary = {
        "schema": SCHEMA + ".validated", "validation_issue_count": 0,
        "root": str(root),
        "prelaunch_receipts_sha256": args.expected_receipts_sha256,
        "data_manifest_sha256": sha256_bytes(data_manifest),
        "design": design, "thermal": thermal,
        "arms": {arm: arm_summary(arm) for arm in ARMS},
        "comparison": comparison,
    }
    encoded = canonical_json(summary)
    if len(encoded) > 64 * 1024 * 1024:
        die("validated summary exceeds the sealed 64 MiB bound")
    write_once(root / "validated_summary.json", encoded)
    print(json.dumps({
        "status": "VALIDATION_OK",
        "cells": sum(len(values) for values in records.values()),
        "repairs": comparison["repairs"],
        "introductions": comparison["introductions"],
        "validated_summary_sha256": sha256_file(
            root / "validated_summary.json"),
        "data_manifest_sha256": sha256_bytes(data_manifest),
    }, indent=2, sort_keys=True))


def command_selftest(_: argparse.Namespace) -> None:
    expected_holdout = tuple(
        "0x" + hashlib.sha256(
            "wirehair.wh2.degree-balanced.allk.holdout.v1|seed|{}".format(
                index).encode("ascii")).hexdigest()[:16]
        for index in range(3)
    )
    if SEEDS != expected_holdout or set(SEEDS) & set(DEVELOPMENT_SEEDS):
        die("independent holdout seed derivation selftest failed")
    tasks = build_tasks(
        "/frozen/base", "/frozen/candidate", k_lo=2, k_hi=19, chunk_max=5)
    audit = independent_task_audit(tasks, k_lo=2, k_hi=19)
    if [task["job"] for task in tasks] != list(range(len(tasks))) or \
            len({task["output_name"] for task in tasks}) != len(tasks) or \
            audit["task_count"] != 72 or audit["task_pairs"] != 36 or \
            audit["total_cells"] != 324:
        die("paired all-K task cardinality selftest failed")
    for offset in range(0, len(tasks), 2):
        base, candidate = tasks[offset:offset + 2]
        if "--degree-balanced-staircase" in base["argv"] or \
                candidate["argv"][-1] != "--degree-balanced-staircase":
            die("A/B architecture isolation selftest failed")
    if parse_preamble("# precodefail: a=1 b=two") != {"a": "1", "b": "two"}:
        die("preamble selftest failed")
    try:
        parse_preamble("# precodefail: duplicate=1 duplicate=2")
    except CampaignError:
        pass
    else:
        die("duplicate preamble rejection selftest failed")

    def expect_reject(function: Any, label: str) -> None:
        try:
            function()
        except CampaignError:
            return
        die("{} tamper selftest was accepted".format(label))

    def thermal_bytes(monotonic_values: Sequence[str], *, edac_ce: str = "0",
                      blank_dimm: bool = False) -> bytes:
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=THERMAL_FIELDS, lineterminator="\n")
        writer.writeheader()
        for index, monotonic_value in enumerate(monotonic_values):
            row = {key: "0" for key in THERMAL_FIELDS}
            row.update({
                "utc": "2026-07-18T00:00:{:02d}.000Z".format(index),
                "monotonic_s": monotonic_value, "cpu_busy_pct": "99.5",
                "cpu_avg_mhz": "5000", "cpu_tctl_c": "60",
                "dimm_read_errors": "0", "load1": "128", "load5": "128",
                "load15": "128", "edac_ce": edac_ce, "edac_ue": "0",
            })
            for key in THERMAL_FIELDS[5:13]:
                row[key] = "45"
            if blank_dimm:
                row[THERMAL_FIELDS[5]] = ""
            writer.writerow(row)
        return output.getvalue().encode("ascii")

    def benchmark_bytes(task: Mapping[str, object]) -> bytes:
        row = {key: "0" for key in CSV_FIELDS}
        row.update({
            "N": str(task["Ks"][0]), "bb": "64", "heavy_family": "periodic",
            "mix_count": "2", "overhead": "0", "trials": "1",
            "success": "1", "rank_fail": "0", "error": "0", "fail_rate": "0",
            "inact_mu": "12", "inact_max": "12",
            "binary_def_mu": "12", "binary_def_max": "12",
            "heavy_gain_mu": "12", "heavy_gain_min": "12",
            "heavy_shortfall": "0", "seed_attempt": "0",
            "block_xors_mu": "10", "block_muladds_mu": "2",
            "first_rank_fail": "-1", "binary_def_hist": "12:1",
            "heavy_gain_hist": "12:1", "failure_trials": "",
            "active_packet_peel_seed_xor": "0x0",
        })
        output = io.StringIO(newline="")
        output.write("# precodefail: " + " ".join(
            "{}={}".format(key, value)
            for key, value in expected_preamble(task).items()) + "\n")
        writer = csv.DictWriter(output, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
        return output.getvalue().encode("ascii")

    def make_success_fixture(root: Path) -> Tuple[List[dict], dict]:
        for directory in (
                "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (root / directory).mkdir(parents=True, exist_ok=True)
        task = {
            "job": 0, "arm": "base", "stage": STAGES[0],
            "seed_index": 0, "seed": SEEDS[0], "schedule": SCHEDULES[0],
            "Ks": [2], "output_name": "job.csv", "argv": ["/bin/true"],
            "stdout_max_bytes": MAX_STDOUT_BYTES_PER_TASK,
            "stderr_max_bytes": MAX_STDERR_BYTES_PER_TASK,
        }
        intent = {
            "schema": SCHEMA + ".segment_intent", "segment": 0,
            "created_utc": "2026-07-18T00:00:00+00:00",
            "created_monotonic_s": 99.5, "stage": STAGES[0], "jobs": [0],
            "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
            "previously_complete": 0, "resume": False,
            "retry_policy": "stage-atomic non-selective retry",
        }
        write_once(segment_path(root, 0, "intent"), canonical_json(intent))
        (root / "attempts" / "segment000").mkdir()
        (root / "attempts" / "segment000" / "job00000.pid").write_bytes(
            canonical_json({"pid": 99999999, "start_ticks": 1}))
        ready = {
            "schema": SCHEMA + ".segment_ready", "segment": 0,
            "intent_sha256": sha256_file(segment_path(root, 0, "intent")),
            "ready_utc": "2026-07-18T00:00:00+00:00",
            "jobs_started_monotonic_s": 101.1, "sampler_pid": 99999999,
            "sampler_start_ticks": 1,
            "samples_at_ready": 2, "last_sample_monotonic_s": "101.0",
        }
        write_once(segment_path(root, 0, "ready"), canonical_json(ready))
        thermal_csv = root / "thermal" / "segment000.csv"
        thermal_stderr = root / "thermal" / "segment000.stderr"
        thermal_csv.write_bytes(thermal_bytes(("100.0", "101.0", "102.0")))
        thermal_stderr.write_bytes(b"")
        raw = root / "raw" / "job.csv"
        stderr = root / "stderr" / "job.csv.stderr"
        exit_path = root / "exit" / "job.csv.exit"
        raw.write_bytes(benchmark_bytes(task))
        stderr.write_bytes(b"")
        exit_path.write_bytes(b"0\n")
        (root / "attempts" / "segment000" / "job00000.stdout").write_bytes(
            raw.read_bytes())
        (root / "attempts" / "segment000" / "job00000.stderr").write_bytes(b"")
        (root / "attempts" / "segment000" / "job00000.exit").write_bytes(b"0\n")
        write_once(job_receipt_path(root, 0), canonical_json({
            "schema": SCHEMA + ".job", "job": 0, "segment": 0,
            "stage": STAGES[0], "output_name": "job.csv", "returncode": 0,
            "started_monotonic_s": 101.2, "ended_monotonic_s": 101.4,
            "stdout_sha256": sha256_file(raw), "stdout_bytes": raw.stat().st_size,
            "stderr_sha256": sha256_file(stderr), "stderr_bytes": 0,
            "exit_sha256": sha256_file(exit_path),
        }))
        write_segment_final(
            root, 0, intent, "success", failures=[], rolled_back=[],
            process_actions=[], jobs_ended_monotonic_s=101.5,
            sampler_returncode=0,
        )
        return [task], {"worker_count": 1}

    with tempfile.TemporaryDirectory(
            prefix="wirehair-degree-balanced-allk-selftest-") as temporary:
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
        cadence_csv.write_bytes(thermal_bytes(("100.0", "101.0", "104.0")))
        cadence_final = load_json(segment_path(cadence, 0, "final"))
        cadence_final["thermal_csv_sha256"] = sha256_file(cadence_csv)
        segment_path(cadence, 0, "final").write_bytes(canonical_json(cadence_final))
        expect_reject(
            lambda: validate_thermal(cadence, fixture_design, fixture_tasks),
            "thermal cadence",
        )

        dimm = tampered_copy("dimm")
        dimm_csv = dimm / "thermal" / "segment000.csv"
        dimm_csv.write_bytes(thermal_bytes(
            ("100.0", "101.0", "102.0"), blank_dimm=True))
        dimm_final = load_json(segment_path(dimm, 0, "final"))
        dimm_final["thermal_csv_sha256"] = sha256_file(dimm_csv)
        segment_path(dimm, 0, "final").write_bytes(canonical_json(dimm_final))
        expect_reject(
            lambda: validate_thermal(dimm, fixture_design, fixture_tasks),
            "missing DIMM",
        )

        edac = tampered_copy("edac")
        edac_csv = edac / "thermal" / "segment000.csv"
        edac_csv.write_bytes(thermal_bytes(
            ("100.0", "101.0", "102.0"), edac_ce="1"))
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
        segment_path(coverage, 0, "ready").write_bytes(
            canonical_json(coverage_ready))
        coverage_final = load_json(segment_path(coverage, 0, "final"))
        coverage_final["ready_sha256"] = sha256_file(
            segment_path(coverage, 0, "ready"))
        segment_path(coverage, 0, "final").write_bytes(
            canonical_json(coverage_final))
        expect_reject(
            lambda: validate_thermal(coverage, fixture_design, fixture_tasks),
            "thermal start coverage",
        )

        sampler_exit = tampered_copy("sampler-exit")
        sampler_final = load_json(segment_path(sampler_exit, 0, "final"))
        sampler_final["sampler_returncode"] = 9
        segment_path(sampler_exit, 0, "final").write_bytes(
            canonical_json(sampler_final))
        expect_reject(
            lambda: validate_thermal(
                sampler_exit, fixture_design, fixture_tasks),
            "sampler exit",
        )

        interrupted = Path(temporary) / "interrupted"
        for directory in (
                "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (interrupted / directory).mkdir(parents=True, exist_ok=True)
        retry_task = {
            "job": 0, "stage": STAGES[0], "output_name": "retry.csv",
            "argv": ["/bin/true"], "stdout_max_bytes": 1024,
            "stderr_max_bytes": 1024,
        }
        retry_intent = {
            "schema": SCHEMA + ".segment_intent", "segment": 0,
            "stage": STAGES[0], "jobs": [0],
            "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
            "retry_policy": "stage-atomic non-selective retry",
        }
        write_once(
            segment_path(interrupted, 0, "intent"), canonical_json(retry_intent))
        (interrupted / "attempts" / "segment000").mkdir()
        (interrupted / "raw" / "retry.csv").write_bytes(b"failed output\n")
        (interrupted / "stderr" / "retry.csv.stderr").write_bytes(b"failure\n")
        (interrupted / "exit" / "retry.csv.exit").write_bytes(b"1\n")
        reconciled = reconcile_incomplete_segments(interrupted, [retry_task])
        if len(reconciled) != 1 or reconciled[0]["state"] != "interrupted" or \
                completed_jobs(interrupted, [retry_task]) or \
                any((interrupted / directory / name).exists()
                    for directory, name in (
                        ("raw", "retry.csv"), ("stderr", "retry.csv.stderr"),
                        ("exit", "retry.csv.exit"))) or \
                reconcile_incomplete_segments(interrupted, [retry_task]):
            die("interrupted failed-triplet resume selftest failed")
        expect_reject(
            lambda: validate_thermal_rows([], segment=0),
            "header-only thermal",
        )
        bound_task = {
            "arm": "base", "seed_index": 0, "schedule": SCHEDULES[0],
            "Ks": [2], "stdout_max_bytes": 8,
        }
        expect_reject(
            lambda: parse_benchmark_csv(b"0123456789", bound_task),
            "stdout bound",
        )

        foreign_pid = Path(temporary) / "foreign-pid"
        for directory in (
                "segments", "thermal", "attempts", "frozen"):
            (foreign_pid / directory).mkdir(parents=True, exist_ok=True)
        (foreign_pid / "attempts" / "segment000").mkdir()
        sleeper = subprocess.Popen(["sleep", "30"])
        try:
            (foreign_pid / "thermal" / "segment000.pid").write_text(
                "{}\n".format(sleeper.pid), encoding="ascii")
            actions = terminate_owned_segment_processes(
                foreign_pid, 0, [])
            if actions != [{
                    "kind": "thermal", "pid": sleeper.pid,
                    "action": "stale_changed_command",
            }] or sleeper.poll() is not None or \
                    (foreign_pid / "thermal" / "segment000.pid").exists():
                die("foreign thermal PID recovery selftest failed")
        finally:
            if sleeper.poll() is None:
                sleeper.terminate()
            sleeper.wait(timeout=5.0)

        def history_fixture(name: str, states: Sequence[Optional[str]]) -> Path:
            history_root = Path(temporary) / name
            for directory in ("segments", "thermal", "attempts"):
                (history_root / directory).mkdir(parents=True, exist_ok=True)
            for segment, state in enumerate(states):
                write_once(segment_path(history_root, segment, "intent"), canonical_json({
                    "schema": SCHEMA + ".segment_intent", "segment": segment,
                    "stage": STAGES[0], "jobs": [0],
                }))
                if state is not None:
                    write_once(segment_path(history_root, segment, "final"), canonical_json({
                        "schema": SCHEMA + ".segment_final", "segment": segment,
                        "state": state,
                    }))
            return history_root

        capped = history_fixture(
            "attempt-cap", ("failed", "interrupted", None))
        if stage_attempt_counts(capped)[STAGES[0]] != MAX_ATTEMPTS_PER_STAGE:
            die("stage attempt history did not count every intent")
        expect_reject(
            lambda: run_segment(
                capped, {"worker_count": 1}, [], STAGES[0], 0, True),
            "stage attempt cap",
        )
        if segment_path(capped, MAX_ATTEMPTS_PER_STAGE, "intent").exists():
            die("attempt-cap rejection wrote a new segment intent")

        succeeded = history_fixture("successful-stage", ("success",))
        expect_reject(
            lambda: run_segment(
                succeeded, {"worker_count": 1}, [], STAGES[0], 0, True),
            "post-success stage attempt",
        )
        if segment_path(succeeded, 1, "intent").exists():
            die("post-success rejection wrote a new segment intent")
    print("wh2_degree_balanced_allk_campaign selftest: PASS")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="freeze a new campaign; never launches jobs")
    prepare.add_argument("--output-root", required=True)
    prepare.add_argument("--candidate-repo", required=True)
    prepare.add_argument("--sampler", required=True)
    prepare.add_argument("--workers", type=int, default=DEFAULT_WORKERS, choices=range(1, 257))
    prepare.set_defaults(function=command_prepare)
    launch = sub.add_parser("launch", help="preflight or execute a frozen campaign")
    launch.add_argument("--root", required=True)
    launch.add_argument("--expected-receipts-sha256", required=True)
    mode = launch.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight-only", action="store_true")
    mode.add_argument("--execute", action="store_true")
    launch.add_argument("--resume", action="store_true")
    launch.add_argument("--stage", choices=(*STAGES, "all"), default="all")
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
