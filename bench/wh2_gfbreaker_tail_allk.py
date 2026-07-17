#!/usr/bin/env python3
"""Freeze, run, and analyze a fresh all-K tail-anchored R2 confirmation.

The evidentiary domain is deterministically derived only after the candidate
source and the successful 540-cell screen artifact are immutable.  It never
reads any H15-v5 commitment or holdout material.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Iterable, Sequence

import wh2_gfbreaker_tail_anchored_screen as common


getcontext().prec = 60

SCHEMA = "wirehair.wh2.gfbreaker_tail_allk.v1"
DOMAIN = b"wirehair.wh2.gfbreaker.tail-anchored.all-k.v1\0"
EXPECTED_SCREEN_SEAL_SHA256 = (
    "6f00bfe1471f91fb23e0b51540da3e2761888edb96edd457f08f294425cfa38d"
)
DEFAULT_SCREEN_SEAL = Path(
    "/tmp/wh2-gfbreaker-tail-anchored-screen-20260717T160000Z-v3/"
    "artifact.sha256"
)
DEFAULT_CURRENT_BINARY = common.DEFAULT_CURRENT_BINARY
EXPECTED_CURRENT_BINARY_SHA256 = common.EXPECTED_CURRENT_BINARY_SHA256
DEFAULT_THERMAL = Path("/tmp/wirehair-enoq-thermal.csv")
GROUP_COUNT = 120
K_FIRST = 2
K_LAST = 64000
K_COUNT = K_LAST - K_FIRST + 1
SCHEDULES = common.SCHEDULES
LOSSES = ("0.35", "0.50", "0.65")
OUTPUT_LOSSES = {
    "0.35": "0.34999999999999998",
    "0.50": "0.5",
    "0.65": "0.65000000000000002",
}
ARMS = (
    ("r0", "current", 0),
    ("current_r1", "current", 1),
    ("nested_r1", "nested", 1),
    ("current_r2", "current", 2),
    ("nested_r2", "nested", 2),
)
ARM_NAMES = tuple(value[0] for value in ARMS)
EXPECTED_CELLS = K_COUNT * len(SCHEDULES) * len(LOSSES)
EXPECTED_JOBS = GROUP_COUNT * len(SCHEDULES) * len(LOSSES) * len(ARMS)
THERMAL_HEADER = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
BANNED_TRAINING_SEEDS = {
    0xA6B527B9AE8DE8D7,
    0x75446E1619E81D8A,
    0x007E9DD892A319F5,
    0xA11CE520F84D877E,
}
SEED_DERIVATION_NOTE = (
    "seeds use u64le(SHA256(domain || commit || NUL || "
    "screen_seal_sha256 || NUL || label)[0:8]); 0 falls back to "
    "bytes[8:16] then 1; salts use the nonzero low uint32 of the same "
    "derivation"
)
COMPLEMENT_NOTE = (
    "0.35/0.65 edge strata bitwise-complement the loss seed and packet-row "
    "salt parameters; SplitMix schedules are not claimed to be packetwise "
    "complements"
)


@dataclass(frozen=True)
class Stratum:
    index: int
    schedule: str
    loss: str
    seed: int
    salt: int
    relation: str


@dataclass(frozen=True)
class Job:
    job: int
    arm: str
    binary_name: str
    rows: int
    stratum: int
    group: int

    @property
    def stem(self) -> str:
        return f"job{self.job:04d}"


@dataclass(frozen=True)
class ParsedJob:
    preamble: dict[str, str]
    diagnostics: dict[int, tuple[str, str]]
    rows: dict[int, tuple[str, ...]]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def derive_u64(commit: str, screen_seal: str, label: str) -> int:
    digest = hashlib.sha256(
        DOMAIN + commit.encode("ascii") + b"\0" +
        screen_seal.encode("ascii") + b"\0" + label.encode("ascii")
    ).digest()
    value = int.from_bytes(digest[:8], "little")
    return value or int.from_bytes(digest[8:16], "little") or 1


def derive_u32(commit: str, screen_seal: str, label: str) -> int:
    value = derive_u64(commit, screen_seal, label) & 0xFFFFFFFF
    return value or 1


def make_strata(commit: str, screen_seal: str) -> tuple[Stratum, ...]:
    result: list[Stratum] = []
    for schedule in SCHEDULES:
        edge_seed = derive_u64(commit, screen_seal, f"{schedule}/edge/seed")
        edge_salt = derive_u32(commit, screen_seal, f"{schedule}/edge/salt")
        middle_seed = derive_u64(commit, screen_seal, f"{schedule}/middle/seed")
        middle_salt = derive_u32(commit, screen_seal, f"{schedule}/middle/salt")
        triples = (
            ("0.35", edge_seed, edge_salt, "edge"),
            ("0.50", middle_seed, middle_salt, "independent-middle"),
            (
                "0.65", edge_seed ^ ((1 << 64) - 1),
                edge_salt ^ ((1 << 32) - 1), "seed-salt-complement",
            ),
        )
        for loss, seed, salt, relation in triples:
            result.append(Stratum(
                len(result), schedule, loss, seed, salt, relation,
            ))
    seeds = [value.seed for value in result]
    salts = [value.salt for value in result]
    if (len(result) != 9 or len(set(seeds)) != 9 or len(set(salts)) != 9 or
        any(value == 0 for value in seeds + salts) or
        BANNED_TRAINING_SEEDS.intersection(seeds)):
        common.die("fresh stratum derivation violated uniqueness/freshness")
    for offset in range(0, 9, 3):
        if (result[offset].seed ^ result[offset + 2].seed != (1 << 64) - 1 or
            result[offset].salt ^ result[offset + 2].salt != (1 << 32) - 1):
            common.die("edge complement relation is not exact")
    return tuple(result)


def groups() -> tuple[tuple[int, ...], ...]:
    result = tuple(
        tuple(range(K_FIRST + group, K_LAST + 1, GROUP_COUNT))
        for group in range(GROUP_COUNT)
    )
    flattened = [k for group in result for k in group]
    if len(flattened) != K_COUNT or set(flattened) != set(range(K_FIRST, K_LAST + 1)):
        common.die("all-K groups do not partition the exact domain")
    return result


def jobs() -> tuple[Job, ...]:
    result: list[Job] = []
    for stratum in range(9):
        for group in range(GROUP_COUNT):
            for arm, binary_name, rows in ARMS:
                result.append(Job(
                    len(result), arm, binary_name, rows, stratum, group,
                ))
    if len(result) != EXPECTED_JOBS:
        common.die("job ledger cardinality mismatch")
    return tuple(result)


def write_tsv(path: Path, rows: Iterable[Sequence[object]]) -> None:
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists() or temporary.is_symlink():
        common.die(f"stale temporary output {temporary}")
    with temporary.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def write_csv(path: Path, rows: Iterable[Sequence[object]]) -> None:
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists() or temporary.is_symlink():
        common.die(f"stale temporary output {temporary}")
    with temporary.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def seal_files(root: Path, paths: Sequence[Path], destination: Path) -> str:
    unique: set[Path] = set()
    lines: list[str] = []
    for path in sorted(paths, key=lambda item: str(item.relative_to(root))):
        resolved = common.regular_file(path, "sealed artifact")
        if root.resolve() not in resolved.parents or resolved in unique:
            common.die(f"duplicate/out-of-root sealed artifact {resolved}")
        unique.add(resolved)
        lines.append(f"{common.sha256(resolved)}  {resolved.relative_to(root.resolve())}\n")
    common.atomic_text(destination, "".join(lines))
    return common.sha256(destination)


def verify_seal(root: Path, seal: Path) -> dict[Path, str]:
    root = root.resolve()
    seal = common.regular_file(seal, "seal")
    result: dict[Path, str] = {}
    for number, line in enumerate(seal.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^/].*)", line)
        if not match:
            common.die(f"{seal}:{number}: malformed seal line")
        relative = Path(match.group(2))
        if relative.is_absolute() or ".." in relative.parts:
            common.die(f"{seal}:{number}: escaping path")
        path = common.regular_file(root / relative, "sealed file")
        if root not in path.parents or path in result:
            common.die(f"{seal}:{number}: duplicate/out-of-root path")
        digest = common.sha256(path)
        if digest != match.group(1):
            common.die(f"{seal}:{number}: SHA256 mismatch")
        result[path] = digest
    return result


def git_capture(worktree: Path, *args: str) -> bytes:
    try:
        return subprocess.check_output(
            ["git", "-C", str(worktree), *args], stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        common.die(f"git {' '.join(args)} failed: {exc}")


def require_head_blob(worktree: Path, path: Path, head: str) -> str:
    try:
        relative = path.resolve().relative_to(worktree.resolve())
    except ValueError:
        common.die(f"tracked campaign input is outside worktree: {path}")
    blob = git_capture(worktree, "show", f"{head}:{relative.as_posix()}")
    blob_digest = hashlib.sha256(blob).hexdigest()
    file_digest = common.sha256(path)
    if blob_digest != file_digest:
        common.die(f"campaign input does not exactly match HEAD: {relative}")
    return blob_digest


def require_frozen_executor(out: Path) -> None:
    frozen = out.resolve() / "frozen"
    expected_script = (frozen / Path(__file__).name).resolve()
    expected_helper = (
        frozen / "wh2_gfbreaker_tail_anchored_screen.py"
    ).resolve()
    if Path(__file__).resolve() != expected_script:
        common.die(f"run this phase with the frozen harness {expected_script}")
    if Path(common.__file__).resolve() != expected_helper:
        common.die(f"campaign helper is not the frozen copy {expected_helper}")


def rebuild_nested_binary(worktree: Path, binary: Path) -> tuple[Path, list[str]]:
    binary = binary.resolve()
    build_dir = worktree.resolve() / "build-tail"
    expected_binary = build_dir / "codec/wirehair_v2_bench"
    if binary != expected_binary:
        common.die(
            f"nested binary must use the audited build tree {expected_binary}"
        )
    cache = common.regular_file(build_dir / "CMakeCache.txt", "nested build cache")
    home_lines = [
        line.split("=", 1)[1]
        for line in cache.read_text(encoding="utf-8").splitlines()
        if line.startswith("CMAKE_HOME_DIRECTORY:INTERNAL=")
    ]
    if len(home_lines) != 1 or Path(home_lines[0]).resolve() != worktree.resolve():
        common.die("nested CMake build tree is not configured from this worktree")
    command = [
        "cmake", "--build", str(build_dir), "--target", "wirehair_v2_bench",
        "--clean-first", "-j", str(min(120, os.cpu_count() or 1)),
    ]
    try:
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        common.die(f"nested benchmark rebuild failed: {exc}")
    return common.regular_file(binary, "rebuilt nested binary"), command


def copy_regular(source: Path, destination: Path, executable: bool = False) -> None:
    source = common.regular_file(source, "freeze source")
    with source.open("rb") as input_stream, destination.open("xb") as output_stream:
        shutil.copyfileobj(input_stream, output_stream, 1024 * 1024)
        output_stream.flush()
        os.fsync(output_stream.fileno())
    destination.chmod(0o755 if executable else 0o644)
    if common.sha256(source) != common.sha256(destination):
        common.die(f"freeze copy mismatch for {source}")


def prepare(args: argparse.Namespace) -> None:
    out = args.out.resolve()
    if out.exists() or out.is_symlink():
        common.die(f"prepare refuses existing output {out}")
    worktree = Path(__file__).resolve().parents[1]
    if git_capture(worktree, "diff", "--quiet") or git_capture(
        worktree, "diff", "--cached", "--quiet"
    ):
        common.die("prepare requires no tracked or staged source changes")
    head = git_capture(worktree, "rev-parse", "HEAD").decode("ascii").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", head):
        common.die("invalid source commit")
    campaign_script_sha = require_head_blob(worktree, Path(__file__), head)
    helper = Path(__file__).with_name("wh2_gfbreaker_tail_anchored_screen.py")
    helper_sha = require_head_blob(worktree, helper, head)
    screen_seal = common.regular_file(args.screen_seal, "screen artifact seal")
    screen_digest = common.sha256(screen_seal)
    if screen_digest != EXPECTED_SCREEN_SEAL_SHA256:
        common.die(f"screen seal SHA256 {screen_digest}, want {EXPECTED_SCREEN_SEAL_SHA256}")
    current = common.regular_file(args.current_binary, "current binary")
    nested, nested_build_command = rebuild_nested_binary(
        worktree, args.nested_binary,
    )
    if git_capture(worktree, "rev-parse", "HEAD").decode("ascii").strip() != head:
        common.die("source HEAD changed during nested benchmark rebuild")
    if git_capture(worktree, "diff", "--quiet") or git_capture(
        worktree, "diff", "--cached", "--quiet"
    ):
        common.die("tracked source changed during nested benchmark rebuild")
    if (require_head_blob(worktree, Path(__file__), head) != campaign_script_sha or
        require_head_blob(worktree, helper, head) != helper_sha):
        common.die("campaign source identity changed during nested rebuild")
    current_sha = common.sha256(current)
    nested_sha = common.sha256(nested)
    if current_sha != EXPECTED_CURRENT_BINARY_SHA256 or current_sha == nested_sha:
        common.die("current/nested binary identity invariant failed")

    frozen = out / "frozen"
    frozen.mkdir(parents=True)
    copy_regular(current, frozen / "current_bench", executable=True)
    copy_regular(nested, frozen / "nested_bench", executable=True)
    copy_regular(Path(__file__), frozen / Path(__file__).name)
    copy_regular(helper, frozen / helper.name)

    strata = make_strata(head, screen_digest)
    all_groups = groups()
    all_jobs = jobs()
    write_tsv(frozen / "strata.tsv", [[
        "stratum", "schedule", "loss", "seed", "salt", "relation",
    ]] + [[
        value.index, value.schedule, value.loss, hex(value.seed),
        hex(value.salt), value.relation,
    ] for value in strata])
    write_tsv(frozen / "groups.tsv", [["group", "ks"]] + [[
        index, ",".join(map(str, values)),
    ] for index, values in enumerate(all_groups)])
    write_tsv(frozen / "jobs.tsv", [[
        "job", "arm", "binary", "rows", "stratum", "group",
    ]] + [[
        value.job, value.arm, value.binary_name, value.rows,
        value.stratum, value.group,
    ] for value in all_jobs])
    source = {
        "source_commit": head,
        "tracked_status": git_capture(
            worktree, "status", "--short", "--untracked-files=no",
        ).decode("utf-8"),
        "screen_artifact_seal": str(screen_seal),
        "screen_artifact_seal_sha256": screen_digest,
        "current_binary_source": str(current),
        "nested_binary_source": str(nested),
        "nested_build_directory": str(worktree / "build-tail"),
        "nested_build_source_directory": str(worktree),
        "nested_build_command": nested_build_command,
        "campaign_script_head_blob_sha256": campaign_script_sha,
        "campaign_helper_head_blob_sha256": helper_sha,
    }
    common.atomic_json(frozen / "source.json", source)
    contract = {
        "schema": SCHEMA,
        "prepared_utc": utc_now(),
        "source_commit": head,
        "screen_artifact_seal_sha256": screen_digest,
        "seed_derivation": SEED_DERIVATION_NOTE,
        "domain_hex": DOMAIN.hex(),
        "freshness_banned_training_seeds": [
            hex(value) for value in sorted(BANNED_TRAINING_SEEDS)
        ],
        "K_first": K_FIRST, "K_last": K_LAST, "K_count": K_COUNT,
        "strata": len(strata), "cells": EXPECTED_CELLS,
        "groups": GROUP_COUNT, "arms": list(ARM_NAMES),
        "jobs": EXPECTED_JOBS,
        "block_bytes": 64, "seed_block_bytes": 1280,
        "overhead": 0, "trials": 1,
        "schedules": list(SCHEDULES), "losses": list(LOSSES),
        "complement_note": COMPLEMENT_NOTE,
        "current_binary_sha256": common.sha256(frozen / "current_bench"),
        "nested_binary_sha256": common.sha256(frozen / "nested_bench"),
        "timing_claim_valid": False,
        "thermal_cpu_abort_c": 90, "thermal_dimm_abort_c": 90,
        "thermal_cpu_busy_min_pct": 95, "thermal_stale_seconds": 5,
        "edac_policy": "abort-on-any-increase",
    }
    common.atomic_json(frozen / "contract.json", contract)
    freeze_paths = [
        frozen / "current_bench", frozen / "nested_bench",
        frozen / Path(__file__).name, frozen / helper.name,
        frozen / "strata.tsv", frozen / "groups.tsv", frozen / "jobs.tsv",
        frozen / "source.json", frozen / "contract.json",
    ]
    seal = seal_files(out, freeze_paths, frozen / "frozen.sha256")
    verify_frozen(out)
    print(json.dumps({
        "output": str(out), "source_commit": head,
        "frozen_seal_sha256": seal,
        "strata": [
            {
                "index": value.index, "schedule": value.schedule,
                "loss": value.loss, "seed": hex(value.seed),
                "salt": hex(value.salt), "relation": value.relation,
            }
            for value in strata
        ],
        "cells": EXPECTED_CELLS, "jobs": EXPECTED_JOBS,
    }, sort_keys=True))


def read_tsv(path: Path, header: Sequence[str]) -> list[list[str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        if next(reader, None) != list(header):
            common.die(f"{path}: header mismatch")
        return list(reader)


def verify_frozen(out: Path) -> tuple[dict[str, object], tuple[Stratum, ...], tuple[tuple[int, ...], ...], tuple[Job, ...]]:
    out = out.resolve()
    frozen = out / "frozen"
    hashes = verify_seal(out, frozen / "frozen.sha256")
    expected = {
        (frozen / name).resolve()
        for name in (
            "current_bench", "nested_bench", Path(__file__).name,
            "wh2_gfbreaker_tail_anchored_screen.py", "strata.tsv",
            "groups.tsv", "jobs.tsv", "source.json", "contract.json",
        )
    }
    if set(hashes) != expected:
        common.die("frozen seal does not cover the exact nine-file set")
    enumerated_frozen = {
        common.regular_file(path, "frozen artifact")
        for path in frozen.iterdir() if path.is_file()
    }
    frozen_directories = {path.name for path in frozen.iterdir() if path.is_dir()}
    if (enumerated_frozen != expected | {(frozen / "frozen.sha256").resolve()} or
        not frozen_directories.issubset({"__pycache__"})):
        common.die("frozen directory contains an unsealed or missing artifact")
    contract = json.loads((frozen / "contract.json").read_text(encoding="utf-8"))
    source = json.loads((frozen / "source.json").read_text(encoding="utf-8"))
    if not isinstance(contract, dict) or not isinstance(source, dict):
        common.die("frozen contract/source must be JSON objects")
    source_keys = {
        "source_commit", "tracked_status", "screen_artifact_seal",
        "screen_artifact_seal_sha256", "current_binary_source",
        "nested_binary_source", "campaign_script_head_blob_sha256",
        "campaign_helper_head_blob_sha256", "nested_build_directory",
        "nested_build_source_directory", "nested_build_command",
    }
    if set(source) != source_keys or source.get("tracked_status") != "":
        common.die("frozen source metadata invariant mismatch")
    head = source.get("source_commit")
    if not isinstance(head, str) or not re.fullmatch(r"[0-9a-f]{40}", head):
        common.die("frozen source commit is not canonical")
    build_dir = Path(str(source.get("nested_build_directory")))
    build_source = Path(str(source.get("nested_build_source_directory")))
    build_command = source.get("nested_build_command")
    if (not build_dir.is_absolute() or not build_source.is_absolute() or
        build_dir != build_source / "build-tail" or
        Path(str(source.get("nested_binary_source"))) !=
            build_dir / "codec/wirehair_v2_bench" or
        not isinstance(build_command, list) or len(build_command) != 8 or
        build_command[:7] != [
            "cmake", "--build", str(build_dir), "--target",
            "wirehair_v2_bench", "--clean-first", "-j",
        ] or
        not isinstance(build_command[7], str) or
        not build_command[7].isdigit() or not 1 <= int(build_command[7]) <= 120):
        common.die("frozen nested build provenance is malformed")
    screen_digest = source.get("screen_artifact_seal_sha256")
    if screen_digest != EXPECTED_SCREEN_SEAL_SHA256:
        common.die("frozen source screen seal is not the pinned screen")
    if (source.get("campaign_script_head_blob_sha256") !=
            common.sha256(frozen / Path(__file__).name) or
        source.get("campaign_helper_head_blob_sha256") != common.sha256(
            frozen / "wh2_gfbreaker_tail_anchored_screen.py"
        )):
        common.die("frozen harness/helper do not match recorded HEAD blobs")
    current_digest = common.sha256(frozen / "current_bench")
    nested_digest = common.sha256(frozen / "nested_bench")
    if current_digest != EXPECTED_CURRENT_BINARY_SHA256 or current_digest == nested_digest:
        common.die("frozen current/nested binary identity invariant failed")
    expected_contract = {
        "schema": SCHEMA, "source_commit": head,
        "screen_artifact_seal_sha256": screen_digest,
        "seed_derivation": SEED_DERIVATION_NOTE, "domain_hex": DOMAIN.hex(),
        "freshness_banned_training_seeds": [
            hex(value) for value in sorted(BANNED_TRAINING_SEEDS)
        ],
        "K_first": K_FIRST, "K_last": K_LAST, "K_count": K_COUNT,
        "strata": 9, "cells": EXPECTED_CELLS, "groups": GROUP_COUNT,
        "arms": list(ARM_NAMES), "jobs": EXPECTED_JOBS,
        "block_bytes": 64, "seed_block_bytes": 1280,
        "overhead": 0, "trials": 1, "schedules": list(SCHEDULES),
        "losses": list(LOSSES), "complement_note": COMPLEMENT_NOTE,
        "current_binary_sha256": current_digest,
        "nested_binary_sha256": nested_digest,
        "timing_claim_valid": False,
        "thermal_cpu_abort_c": 90, "thermal_dimm_abort_c": 90,
        "thermal_cpu_busy_min_pct": 95, "thermal_stale_seconds": 5,
        "edac_policy": "abort-on-any-increase",
    }
    if set(contract) != set(expected_contract) | {"prepared_utc"}:
        common.die("frozen contract key set mismatch")
    for key, value in expected_contract.items():
        if contract.get(key) != value:
            common.die(f"frozen contract invariant mismatch for {key}")
    if not isinstance(contract.get("prepared_utc"), str) or not re.fullmatch(
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z",
        str(contract["prepared_utc"]),
    ):
        common.die("frozen preparation timestamp is not canonical")
    stratum_rows = read_tsv(
        frozen / "strata.tsv",
        ("stratum", "schedule", "loss", "seed", "salt", "relation"),
    )
    if any(len(row) != 6 for row in stratum_rows):
        common.die("frozen stratum row field count mismatch")
    strata = tuple(Stratum(
        common.strict_uint(row[0], "stratum index"), row[1], row[2],
        common.canonical_u64(row[3], "stratum seed"),
        common.canonical_u64(row[4], "stratum salt"), row[5],
    ) for row in stratum_rows)
    expected_strata = make_strata(
        head, str(screen_digest)
    )
    if strata != expected_strata:
        common.die("frozen strata differ from deterministic derivation")
    group_rows = read_tsv(frozen / "groups.tsv", ("group", "ks"))
    parsed_groups: list[tuple[int, ...]] = []
    for index, row in enumerate(group_rows):
        if len(row) != 2 or common.strict_uint(row[0], "group") != index:
            common.die("frozen group index mismatch")
        values = tuple(common.strict_uint(item, "group K") for item in row[1].split(","))
        parsed_groups.append(values)
    if tuple(parsed_groups) != groups():
        common.die("frozen groups differ from exact all-K partition")
    job_rows = read_tsv(
        frozen / "jobs.tsv", ("job", "arm", "binary", "rows", "stratum", "group"),
    )
    if any(len(row) != 6 for row in job_rows):
        common.die("frozen job row field count mismatch")
    parsed_jobs = tuple(Job(
        common.strict_uint(row[0], "job"), row[1], row[2],
        common.strict_uint(row[3], "job rows"),
        common.strict_uint(row[4], "job stratum"),
        common.strict_uint(row[5], "job group"),
    ) for row in job_rows)
    if parsed_jobs != jobs():
        common.die("frozen jobs differ from deterministic ledger")
    return contract, strata, tuple(parsed_groups), parsed_jobs


def output_loss(loss: str) -> str:
    try:
        return OUTPUT_LOSSES[loss]
    except KeyError:
        common.die(f"unsupported frozen loss {loss}")


def make_command(out: Path, job: Job, stratum: Stratum, ks: Sequence[int]) -> list[str]:
    binary = out / "frozen" / f"{job.binary_name}_bench"
    result = [
        str(binary), "precodefail", "--N", ",".join(map(str, ks)),
        "--bb-list", "64", "--seed-block-bytes", "1280",
        "--overhead", "0", "--trials", "1", "--threads", "1",
        "--loss", stratum.loss, "--seed", hex(stratum.seed),
        "--schedule", stratum.schedule, "--completion", "mixed",
        "--mix-count", "2", "--packet-peel-seed-xor", str(stratum.salt),
        "--mixed-gf256-rows", "11", "--mixed-gf16-rows", "4",
        "--mixed-period", "32", "--mixed-geometry", "shared-x",
        "--mixed-residue-schedule", "hashed",
        "--mixed-residue-hash-seed", "68", "--mixed-residue-hash-keyed",
        "--mixed-independent-extension-residues",
        "--mixed-extension-residue-seed-xor", "78",
        "--mixed-fused-schedule-buckets", "auto",
    ]
    if job.rows:
        result.extend(("--mixed-independent-gf256-breaker-rows", str(job.rows)))
    return result


def validate_row(row: tuple[str, ...], k: int, salt: int, context: str) -> None:
    if len(row) != len(common.CSV_HEADER):
        common.die(f"{context}: CSV field count")
    if row[:6] != (str(k), "64", "periodic", "2", "0", "1"):
        common.die(f"{context}: fixed geometry mismatch")
    for index in (1, 3, 4, 5, 6, 7, 8, 11, 13, 15, 16, 23):
        common.strict_uint(row[index], f"{context}:{common.CSV_HEADER[index]}")
    for index in (9, 10, 12, 14, 17, 18, 19, 20, 21, 22, 24, 25):
        common.strict_decimal(row[index], f"{context}:{common.CSV_HEADER[index]}")
    success, rank_fail, error = map(int, row[6:9])
    if success + rank_fail + error != 1 or int(row[16]) > rank_fail:
        common.die(f"{context}: outcome mismatch")
    if row[9] != ("1.00000000" if rank_fail + error else "0.00000000"):
        common.die(f"{context}: fail rate mismatch")
    binary_count, binary_sum = common.parse_histogram(row[27], context + ":binary")
    heavy_count, heavy_sum = common.parse_histogram(row[28], context + ":heavy")
    if (binary_count != 1 or heavy_count != 1 or
        Decimal(row[10]) != int(row[11]) or
        Decimal(row[12]) != binary_sum or int(row[13]) != binary_sum or
        Decimal(row[14]) != heavy_sum or int(row[15]) != heavy_sum):
        common.die(f"{context}: histogram summary mismatch")
    if (rank_fail == 0 and row[26] != "-1") or (rank_fail == 1 and row[26] != "0"):
        common.die(f"{context}: first failure mismatch")
    failures = () if not row[29] else tuple(row[29].split("|"))
    if len(failures) != rank_fail + error or any(item != "0" for item in failures):
        common.die(f"{context}: failure trial ledger mismatch")
    if common.canonical_u64(row[30], context) != salt:
        common.die(f"{context}: active salt mismatch")


def parse_job_output(data: bytes, job: Job, stratum: Stratum, ks: Sequence[int]) -> ParsedJob:
    context = job.stem
    try:
        lines = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        common.die(f"{context}: invalid UTF-8: {exc}")
    expected_lines = 2 + len(ks) * (2 if job.rows else 1)
    if len(lines) != expected_lines:
        common.die(f"{context}: lines={len(lines)}, want {expected_lines}")
    preamble = common.parse_preamble(lines[0], context)
    expected = {
        "trials": "1", "threads": "1", "loss": output_loss(stratum.loss),
        "seed": hex(stratum.seed), "completion": "mixed", "mixed_period": "32",
        "mixed_gf256_rows": "11", "mixed_gf16_rows": "4",
        "mixed_geometry": "shared-x", "mixed_residue_skew": "0",
        "mixed_residue_schedule": "hashed", "mixed_residue_hash_seed": "0x44",
        "mixed_residue_hash_keyed": "1",
        "mixed_independent_extension_residues": "1",
        "mixed_extension_residue_seed_xor": "0x4e",
        "mixed_fused_schedule_buckets_mode": "-1",
        "source_hits_override": "0", "packet_peel_seed_xor": hex(stratum.salt),
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1", "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "1280", "overhead_stream": "salted",
        "full_payload_solve": "0", "schedule": stratum.schedule,
    }
    if job.rows:
        expected.update({
            "mixed_independent_gf256_breaker_residues": "1",
            "mixed_independent_gf256_breaker_rows": str(job.rows),
            "mixed_gf256_breaker_residue_seed_xor": "0xb7e15162",
        })
    if preamble != expected:
        common.die(f"{context}: exact preamble mismatch")
    if tuple(next(csv.reader([lines[1]]))) != common.CSV_HEADER:
        common.die(f"{context}: header mismatch")
    diagnostics: dict[int, tuple[str, str]] = {}
    rows: dict[int, tuple[str, ...]] = {}
    cursor = 2
    for k in ks:
        if job.rows:
            match = common.DIAGNOSTIC_RE.fullmatch(lines[cursor])
            if not match or int(match.group(1)) != k or int(match.group(3)) != job.rows:
                common.die(f"{context}: K={k} diagnostic mismatch")
            diagnostics[k] = (match.group(2), match.group(4))
            cursor += 1
        row = tuple(next(csv.reader([lines[cursor]])))
        validate_row(row, k, stratum.salt, f"{context}:K={k}")
        rows[k] = row
        cursor += 1
    if cursor != len(lines) or tuple(rows) != tuple(ks):
        common.die(f"{context}: output order mismatch")
    return ParsedJob(preamble, diagnostics, rows)


def validate_thermal_row(
    values: Sequence[str], baseline_edac: tuple[int, int],
) -> tuple[float, float, float, float]:
    if len(values) != len(THERMAL_HEADER):
        common.die("thermal row field count")
    if not re.fullmatch(
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}Z",
        values[0],
    ):
        common.die("thermal row UTC timestamp is not canonical")
    monotonic = float(values[1])
    busy = float(values[2])
    mhz = float(values[3])
    cpu = float(values[4])
    dimms = [float(value) for value in values[5:13]]
    loads = [float(value) for value in values[14:17]]
    if not all(math.isfinite(value) for value in [
        monotonic, busy, mhz, cpu, *dimms, *loads,
    ]):
        common.die("thermal row contains non-finite value")
    if (monotonic <= 0 or not 0 <= busy <= 100 or mhz <= 0 or cpu < 0 or
        any(value < 0 for value in [*dimms, *loads])):
        common.die("thermal row contains an out-of-range measurement")
    if common.strict_uint(values[13], "thermal DIMM errors") != 0:
        common.die("thermal sampler reported a DIMM read error")
    edac = (
        common.strict_uint(values[17], "EDAC CE"),
        common.strict_uint(values[18], "EDAC UE"),
    )
    if edac != baseline_edac:
        common.die(f"EDAC changed from {baseline_edac} to {edac}")
    if cpu >= 90 or max(dimms) >= 90:
        common.die(f"thermal limit reached CPU={cpu} DIMM={max(dimms)}")
    if busy < 95:
        common.die(f"CPU utilization fell below campaign floor: {busy}%")
    return monotonic, busy, cpu, max(dimms)


class ThermalMonitor:
    def __init__(self, source: Path, destination: Path, abort: threading.Event):
        self.source = common.regular_file(source, "thermal source")
        self.destination = destination
        self.abort = abort
        self.stop_event = threading.Event()
        self.error: BaseException | None = None
        self.thread: threading.Thread | None = None
        self.samples = 0
        self.cpu_busy_min = 100.0
        self.cpu_max = 0.0
        self.dimm_max = 0.0

    def start(self) -> None:
        data = self.source.read_bytes()
        if not data.endswith(b"\n"):
            common.die("thermal source does not end at a row boundary")
        lines = data.decode("utf-8").splitlines()
        if not lines or tuple(next(csv.reader([lines[0]]))) != THERMAL_HEADER:
            common.die("thermal source header mismatch")
        baseline_values = next(csv.reader([lines[-1]]))
        if len(baseline_values) != len(THERMAL_HEADER):
            common.die("thermal source baseline row field count")
        baseline_edac = (
            common.strict_uint(baseline_values[17], "baseline CE"),
            common.strict_uint(baseline_values[18], "baseline UE"),
        )
        mono, busy, cpu, dimm = validate_thermal_row(baseline_values, baseline_edac)
        if abs(time.monotonic() - mono) > 5:
            common.die("thermal source is stale before launch")
        source_status = self.source.stat()
        source_inode = (source_status.st_dev, source_status.st_ino)
        offset = len(data)

        def worker() -> None:
            last_seen = mono
            try:
                with self.source.open("rb") as input_stream, self.destination.open(
                    "x", encoding="utf-8", newline=""
                ) as output_stream:
                    writer = csv.writer(output_stream, lineterminator="\n")
                    writer.writerow(THERMAL_HEADER)
                    writer.writerow(baseline_values)
                    self.samples = 1
                    self.cpu_busy_min = busy
                    self.cpu_max = cpu
                    self.dimm_max = dimm
                    input_stream.seek(offset)
                    while not self.stop_event.is_set():
                        position = input_stream.tell()
                        line = input_stream.readline()
                        if not line:
                            if time.monotonic() - last_seen > 5:
                                common.die("thermal source became stale under load")
                            time.sleep(0.1)
                            continue
                        if not line.endswith(b"\n"):
                            input_stream.seek(position)
                            if time.monotonic() - last_seen > 5:
                                common.die("thermal source ended in a stale partial row")
                            time.sleep(0.05)
                            continue
                        values = next(csv.reader([line.decode("utf-8").rstrip("\n")]))
                        next_mono, next_busy, next_cpu, next_dimm = validate_thermal_row(
                            values, baseline_edac,
                        )
                        if next_mono <= last_seen or abs(time.monotonic() - next_mono) > 5:
                            common.die("thermal source row is stale or nonmonotonic")
                        writer.writerow(values)
                        output_stream.flush()
                        self.samples += 1
                        self.cpu_busy_min = min(self.cpu_busy_min, next_busy)
                        self.cpu_max = max(self.cpu_max, next_cpu)
                        self.dimm_max = max(self.dimm_max, next_dimm)
                        last_seen = next_mono
                    while True:
                        line = input_stream.readline()
                        if not line:
                            break
                        if not line.endswith(b"\n"):
                            break
                        values = next(csv.reader([line.decode("utf-8").rstrip("\n")]))
                        next_mono, next_busy, next_cpu, next_dimm = validate_thermal_row(
                            values, baseline_edac,
                        )
                        if next_mono <= last_seen or abs(time.monotonic() - next_mono) > 5:
                            common.die("final thermal row is stale or nonmonotonic")
                        writer.writerow(values)
                        self.samples += 1
                        self.cpu_busy_min = min(self.cpu_busy_min, next_busy)
                        self.cpu_max = max(self.cpu_max, next_cpu)
                        self.dimm_max = max(self.dimm_max, next_dimm)
                        last_seen = next_mono
                    output_stream.flush()
                    os.fsync(output_stream.fileno())
                final_status = self.source.stat()
                if (final_status.st_dev, final_status.st_ino) != source_inode:
                    common.die("thermal source identity changed during run")
            except BaseException as exc:
                self.error = exc
                self.abort.set()
                kill_active(signal.SIGTERM)

        self.thread = threading.Thread(target=worker, name="thermal-monitor", daemon=True)
        self.thread.start()

    def stop(self) -> None:
        deadline = time.monotonic() + 5
        while (self.samples < 2 and self.error is None and
               not self.abort.is_set() and self.thread and self.thread.is_alive() and
               time.monotonic() < deadline):
            time.sleep(0.05)
        self.stop_event.set()
        if self.thread:
            self.thread.join(timeout=10)
            if self.thread.is_alive():
                common.die("thermal monitor did not stop")
        if self.error:
            raise self.error
        if self.samples < 2:
            common.die("thermal monitor captured fewer than two samples")


ACTIVE_LOCK = threading.Lock()
ACTIVE: dict[int, subprocess.Popen[bytes]] = {}


def kill_active(signum: int = signal.SIGTERM) -> None:
    with ACTIVE_LOCK:
        processes = list(ACTIVE.values())
    for process in processes:
        if process.poll() is not None:
            continue
        try:
            os.killpg(process.pid, signum)
        except ProcessLookupError:
            pass


def run_one(
    out: Path, raw: Path, job: Job, stratum: Stratum, ks: Sequence[int],
    cpu: int, taskset: str, timeout: float, environment: dict[str, str],
    abort: threading.Event,
) -> bool:
    stdout_path = raw / f"{job.stem}.stdout"
    stderr_path = raw / f"{job.stem}.stderr"
    if stdout_path.exists() or stderr_path.exists():
        if not stdout_path.is_file() or not stderr_path.is_file() or stderr_path.stat().st_size:
            common.die(f"{job.stem}: invalid partial/existing raw pair")
        parse_job_output(stdout_path.read_bytes(), job, stratum, ks)
        return False
    command = [taskset, "-c", str(cpu), *make_command(out, job, stratum, ks)]
    with ACTIVE_LOCK:
        if abort.is_set():
            common.die("campaign aborted before job launch")
        process = subprocess.Popen(
            command, cwd=out, env=environment, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True,
        )
        ACTIVE[job.job] = process
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        stdout, stderr = process.communicate()
        common.die(
            f"{job.stem}: timeout after {timeout} seconds; "
            f"discarded stdout={len(stdout)} stderr={len(stderr)} bytes"
        )
    finally:
        with ACTIVE_LOCK:
            ACTIVE.pop(job.job, None)
    if process.returncode != 0 or stderr:
        common.die(
            f"{job.stem}: exit={process.returncode}, discarded "
            f"stdout={len(stdout)} stderr={len(stderr)} bytes"
        )
    parse_job_output(stdout, job, stratum, ks)
    common.atomic_bytes(stdout_path, stdout)
    common.atomic_bytes(stderr_path, stderr)
    return True


def run_campaign(args: argparse.Namespace) -> None:
    out = args.out.resolve()
    require_frozen_executor(out)
    contract, strata, all_groups, all_jobs = verify_frozen(out)
    if (args.workers <= 0 or args.workers > GROUP_COUNT or
        args.timeout_seconds <= 0):
        common.die(f"workers must be in [1,{GROUP_COUNT}] and timeout positive")
    taskset = shutil.which("taskset")
    if not taskset:
        common.die("taskset is required")
    cpus = sorted(os.sched_getaffinity(0))[:args.workers]
    if not cpus:
        common.die("empty CPU affinity")
    results = out / "results"
    raw = results / "raw"
    thermal_dir = results / "thermal"
    raw.mkdir(parents=True, exist_ok=True)
    thermal_dir.mkdir(parents=True, exist_ok=True)
    if (results / "phase_complete.sha256").exists():
        verify_results(out)
        print(json.dumps({"already_complete": True, "output": str(out)}, sort_keys=True))
        return
    existing_segments: list[int] = []
    for path in thermal_dir.iterdir():
        match = re.fullmatch(r"segment([0-9]{3})\.csv", path.name)
        if not match or not path.is_file():
            common.die(f"unexpected thermal segment artifact {path}")
        existing_segments.append(int(match.group(1)))
    if sorted(existing_segments) != list(range(len(existing_segments))):
        common.die("thermal segment indices are not contiguous from zero")
    segment = len(existing_segments)
    if segment > 999:
        common.die("thermal segment index exhausted")
    thermal_path = thermal_dir / f"segment{segment:03d}.csv"
    abort = threading.Event()
    monitor = ThermalMonitor(args.thermal, thermal_path, abort)
    environment = dict(os.environ)
    for key in list(environment):
        if key.startswith("WIREHAIR_") or key in {
            "ASAN_OPTIONS", "UBSAN_OPTIONS", "LSAN_OPTIONS", "TSAN_OPTIONS",
        }:
            environment.pop(key, None)
    environment.update({"LC_ALL": "C", "TZ": "UTC"})
    by_group: list[list[Job]] = [[] for _ in range(GROUP_COUNT)]
    for job in all_jobs:
        by_group[job.group].append(job)
    for batch in by_group:
        batch.sort(key=lambda value: (value.stratum, ARM_NAMES.index(value.arm)))
    bins: list[list[int]] = [[] for _ in cpus]
    for group in range(GROUP_COUNT):
        bins[group % len(cpus)].append(group)
    total = len(all_jobs)
    already = sum(
        (raw / f"{job.stem}.stdout").is_file() and
        (raw / f"{job.stem}.stderr").is_file()
        for job in all_jobs
    )
    completed = already
    progress_lock = threading.Lock()
    start = time.monotonic()

    def run_bin(slot: int, group_indices: Sequence[int]) -> None:
        nonlocal completed
        cpu = cpus[slot]
        try:
            for group in group_indices:
                for job in by_group[group]:
                    fresh_job = run_one(
                        out, raw, job, strata[job.stratum], all_groups[job.group],
                        cpu, taskset, args.timeout_seconds, environment, abort,
                    )
                    with progress_lock:
                        if fresh_job:
                            completed += 1
                        fresh = completed - already
                        if fresh_job and (
                            fresh % 25 == 0 or completed == total
                        ):
                            elapsed = time.monotonic() - start
                            rate = fresh / elapsed
                            eta = (total - completed) / rate if rate else 0
                            print(
                                f"progress={completed}/{total} fresh={fresh} "
                                f"elapsed_s={elapsed:.1f} eta_s={eta:.1f}",
                                flush=True,
                            )
        except BaseException:
            abort.set()
            kill_active(signal.SIGTERM)
            raise

    previous_signals = {
        signum: signal.getsignal(signum)
        for signum in (signal.SIGINT, signal.SIGTERM)
    }

    def handle_signal(signum: int, _frame: object) -> None:
        abort.set()
        kill_active(signal.SIGTERM)
        raise KeyboardInterrupt(f"campaign interrupted by signal {signum}")

    for signum in previous_signals:
        signal.signal(signum, handle_signal)
    try:
        monitor.start()
        try:
            with ThreadPoolExecutor(max_workers=len(cpus)) as executor:
                futures = [
                    executor.submit(run_bin, slot, group_indices)
                    for slot, group_indices in enumerate(bins) if group_indices
                ]
                for future in as_completed(futures):
                    future.result()
                    if abort.is_set():
                        kill_active(signal.SIGTERM)
                        common.die("thermal monitor aborted campaign")
        except BaseException:
            abort.set()
            kill_active(signal.SIGTERM)
            time.sleep(0.2)
            kill_active(signal.SIGKILL)
            raise
        finally:
            monitor.stop()
    finally:
        for signum, previous in previous_signals.items():
            signal.signal(signum, previous)
    if abort.is_set():
        common.die("campaign aborted")
    for job in all_jobs:
        stdout_path = raw / f"{job.stem}.stdout"
        stderr_path = raw / f"{job.stem}.stderr"
        if not stdout_path.is_file() or not stderr_path.is_file() or stderr_path.stat().st_size:
            common.die(f"{job.stem}: missing/nonempty final raw artifact")
        parse_job_output(
            stdout_path.read_bytes(), job, strata[job.stratum], all_groups[job.group],
        )
    run_record = results / f"run_segment{segment:03d}.json"
    common.atomic_json(run_record, {
        "schema": SCHEMA + ".run",
        "segment": segment, "completed_utc": utc_now(),
        "elapsed_seconds": str(Decimal(str(time.monotonic() - start))),
        "workers": len(cpus), "cpus": cpus,
        "jobs_already_present": already, "jobs_completed": total,
        "thermal_samples": monitor.samples,
        "thermal_cpu_busy_min_pct": monitor.cpu_busy_min,
        "thermal_cpu_max_c": monitor.cpu_max,
        "thermal_dimm_max_c": monitor.dimm_max,
        "timing_claim_valid": False,
    })
    result_paths = [
        path for path in raw.iterdir() if path.is_file()
    ] + [
        path for path in thermal_dir.iterdir() if path.is_file()
    ] + [
        path for path in results.glob("run_segment*.json") if path.is_file()
    ]
    seal = seal_files(out, result_paths, results / "phase_complete.sha256")
    verify_results(out)
    print(json.dumps({
        "output": str(out), "jobs": total, "cells": EXPECTED_CELLS,
        "phase_seal_sha256": seal, "thermal_samples": monitor.samples,
        "thermal_cpu_busy_min_pct": monitor.cpu_busy_min,
        "thermal_cpu_max_c": monitor.cpu_max,
        "thermal_dimm_max_c": monitor.dimm_max,
    }, sort_keys=True))


def verify_results(out: Path) -> dict[Path, str]:
    out = out.resolve()
    _, strata, all_groups, all_jobs = verify_frozen(out)
    results = out / "results"
    hashes = verify_seal(out, results / "phase_complete.sha256")
    raw = results / "raw"
    expected_raw = {
        (raw / f"{job.stem}.{extension}").resolve()
        for job in all_jobs for extension in ("stdout", "stderr")
    }
    actual_raw = {path for path in hashes if path.parent == raw.resolve()}
    enumerated_raw = {
        common.regular_file(path, "raw result")
        for path in raw.iterdir()
    }
    if actual_raw != expected_raw or enumerated_raw != expected_raw:
        common.die("results seal raw file set mismatch")
    thermal = [path for path in hashes if path.parent == (results / "thermal").resolve()]
    run_records = [
        path for path in hashes
        if path.parent == results.resolve() and path.name.startswith("run_segment")
    ]
    enumerated_thermal = {
        common.regular_file(path, "thermal result")
        for path in (results / "thermal").iterdir()
    }
    enumerated_root = {
        common.regular_file(path, "result root file")
        for path in results.iterdir() if path.is_file() and
        path.name != "phase_complete.sha256"
    }
    result_directories = {
        path.name for path in results.iterdir() if path.is_dir()
    }
    if (not thermal or not run_records or
        set(hashes) != expected_raw | set(thermal) | set(run_records) or
        enumerated_thermal != set(thermal) or enumerated_root != set(run_records) or
        result_directories != {"raw", "thermal"}):
        common.die("results seal non-raw file set mismatch")
    thermal_by_segment: dict[int, Path] = {}
    for path in thermal:
        match = re.fullmatch(r"segment([0-9]{3})\.csv", path.name)
        if not match or int(match.group(1)) in thermal_by_segment:
            common.die("sealed thermal segment name/index mismatch")
        thermal_by_segment[int(match.group(1))] = path
    if sorted(thermal_by_segment) != list(range(len(thermal_by_segment))):
        common.die("sealed thermal segment indices are not contiguous")
    records_by_segment: dict[int, Path] = {}
    for path in run_records:
        match = re.fullmatch(r"run_segment([0-9]{3})\.json", path.name)
        if not match or int(match.group(1)) in records_by_segment:
            common.die("sealed run record name/index mismatch")
        records_by_segment[int(match.group(1))] = path
    if (not set(records_by_segment).issubset(thermal_by_segment) or
        max(thermal_by_segment) not in records_by_segment):
        common.die("sealed run records do not cover the completed segment")
    for job in all_jobs:
        stderr = raw / f"{job.stem}.stderr"
        if stderr.stat().st_size:
            common.die(f"{job.stem}: sealed stderr is nonempty")
        parse_job_output(
            (raw / f"{job.stem}.stdout").read_bytes(),
            job, strata[job.stratum], all_groups[job.group],
        )
    baseline_edac: tuple[int, int] | None = None
    thermal_stats: dict[int, tuple[int, float, float, float]] = {}
    previous_mono: float | None = None
    previous_values: tuple[str, ...] | None = None
    for segment, path in sorted(thermal_by_segment.items()):
        count = 0
        busy_min = 100.0
        cpu_max = 0.0
        dimm_max = 0.0
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.reader(stream)
            if tuple(next(reader, ())) != THERMAL_HEADER:
                common.die(f"{path}: thermal header mismatch")
            for values in reader:
                if baseline_edac is None:
                    if len(values) != len(THERMAL_HEADER):
                        common.die(f"{path}: thermal row field count")
                    baseline_edac = (
                        common.strict_uint(values[17], "thermal baseline CE"),
                        common.strict_uint(values[18], "thermal baseline UE"),
                    )
                mono, busy, cpu, dimm = validate_thermal_row(values, baseline_edac)
                if previous_mono is not None and mono <= previous_mono:
                    boundary_duplicate = (
                        count == 0 and mono == previous_mono and
                        tuple(values) == previous_values
                    )
                    if not boundary_duplicate:
                        common.die(
                            "sealed thermal samples regress or duplicate "
                            "within a segment"
                        )
                previous_mono = mono
                previous_values = tuple(values)
                count += 1
                busy_min = min(busy_min, busy)
                cpu_max = max(cpu_max, cpu)
                dimm_max = max(dimm_max, dimm)
        if count < 2:
            common.die(f"sealed thermal segment {segment} has fewer than two samples")
        thermal_stats[segment] = (count, busy_min, cpu_max, dimm_max)
    run_keys = {
        "schema", "segment", "completed_utc", "elapsed_seconds", "workers",
        "cpus", "jobs_already_present", "jobs_completed", "thermal_samples",
        "thermal_cpu_busy_min_pct", "thermal_cpu_max_c", "thermal_dimm_max_c",
        "timing_claim_valid",
    }
    for segment, path in sorted(records_by_segment.items()):
        record = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(record, dict) or set(record) != run_keys:
            common.die(f"{path}: run record key set mismatch")
        elapsed = common.strict_decimal(
            str(record.get("elapsed_seconds")), f"{path}:elapsed_seconds",
        )
        cpus = record.get("cpus")
        if (record.get("schema") != SCHEMA + ".run" or
            record.get("segment") != segment or elapsed <= 0 or
            not isinstance(record.get("elapsed_seconds"), str) or
            not isinstance(record.get("workers"), int) or
            not isinstance(cpus, list) or
            record["workers"] != len(cpus) or
            not 1 <= record["workers"] <= GROUP_COUNT or
            any(not isinstance(cpu, int) or cpu < 0 for cpu in cpus) or
            len(set(cpus)) != len(cpus) or
            not isinstance(record.get("jobs_already_present"), int) or
            not 0 <= record["jobs_already_present"] <= EXPECTED_JOBS or
            record.get("jobs_completed") != EXPECTED_JOBS or
            record.get("timing_claim_valid") is not False or
            not isinstance(record.get("completed_utc"), str) or
            not re.fullmatch(
                r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z",
                record["completed_utc"],
            )):
            common.die(f"{path}: run record invariant mismatch")
        count, busy_min, cpu_max, dimm_max = thermal_stats[segment]
        if (record.get("thermal_samples") != count or
            record.get("thermal_cpu_busy_min_pct") != busy_min or
            record.get("thermal_cpu_max_c") != cpu_max or
            record.get("thermal_dimm_max_c") != dimm_max):
            common.die(f"{path}: run record/thermal statistics mismatch")
    return hashes


def failed(row: Sequence[str]) -> bool:
    return row[7] != "0" or row[8] != "0"


def milli(text: str, context: str) -> int:
    value = common.strict_decimal(text, context) * 1000
    if value != value.to_integral_value():
        common.die(f"{context}: non-integral milli-unit")
    return int(value)


def update_comparison(counter: Counter[str], label: str, base: bool, candidate: bool) -> None:
    if base and not candidate:
        counter[f"nested_repair_{label}"] += 1
    elif not base and candidate:
        counter[f"nested_intro_{label}"] += 1


def binomial_one_sided(repairs: int, introductions: int) -> dict[str, object]:
    discordant = repairs + introductions
    if discordant == 0:
        return {"numerator": 1, "denominator": 1, "decimal": "1"}
    numerator = sum(math.comb(discordant, k) for k in range(introductions + 1))
    denominator = 1 << discordant
    return {
        "numerator": numerator, "denominator": denominator,
        "decimal": str(Decimal(numerator) / Decimal(denominator)),
    }


CELL_HEADER = (
    "K", "stratum", "schedule", "loss", "seed", "salt",
    *tuple(
        value
        for arm in ARM_NAMES
        for value in (
            f"{arm}_failed", f"{arm}_error", f"{arm}_heavy_shortfall",
            f"{arm}_block_xors_milli", f"{arm}_block_muladds_milli",
        )
    ),
    "r1_exact_identity",
)
TOTAL_KEYS = (
    "cells",
    *tuple(
        value
        for arm in ARM_NAMES
        for value in (
            f"{arm}_fail", f"{arm}_error", f"{arm}_heavy_shortfall",
        )
    ),
    "nested_repair_r0", "nested_intro_r0",
    "nested_repair_r1", "nested_intro_r1",
    "nested_repair_current_r2", "nested_intro_current_r2",
    "r1_identity_mismatch", "current_nested_common_success",
    "current_r2_block_xors_milli_common",
    "nested_r2_block_xors_milli_common",
    "current_r2_block_muladds_milli_common",
    "nested_r2_block_muladds_milli_common",
)
PER_K_HEADER = (
    "K", "cells",
    *tuple(
        value
        for arm in ARM_NAMES
        for value in (
            f"{arm}_fail", f"{arm}_error", f"{arm}_heavy_shortfall",
        )
    ),
    "nested_repair_r0", "nested_intro_r0", "nested_repair_r1",
    "nested_intro_r1", "nested_repair_current_r2",
    "nested_intro_current_r2", "r1_identity_mismatch",
)


def cell_values(
    k: int, stratum: Stratum, outputs: dict[str, ParsedJob],
) -> list[object]:
    result: list[object] = [
        k, stratum.index, stratum.schedule, stratum.loss,
        hex(stratum.seed), hex(stratum.salt),
    ]
    for arm in ARM_NAMES:
        row = outputs[arm].rows[k]
        result.extend((
            int(failed(row)), int(row[8]), int(row[16]),
            milli(row[24], f"K={k}/{arm}/xor"),
            milli(row[25], f"K={k}/{arm}/muladd"),
        ))
    current = outputs["current_r1"]
    nested = outputs["nested_r1"]
    identity = (
        current.preamble == nested.preamble and
        current.diagnostics[k] == nested.diagnostics[k] and
        tuple(current.rows[k][index] for index in common.NON_TIMING_INDICES) ==
        tuple(nested.rows[k][index] for index in common.NON_TIMING_INDICES)
    )
    result.append(int(identity))
    return result


def update_totals(counter: Counter[str], values: Sequence[object]) -> None:
    row = dict(zip(CELL_HEADER, values))
    counter["cells"] += 1
    for arm in ARM_NAMES:
        counter[f"{arm}_fail"] += int(row[f"{arm}_failed"])
        counter[f"{arm}_error"] += int(row[f"{arm}_error"])
        counter[f"{arm}_heavy_shortfall"] += int(row[f"{arm}_heavy_shortfall"])
    nested = bool(row["nested_r2_failed"])
    for label, arm in (
        ("r0", "r0"), ("r1", "current_r1"), ("current_r2", "current_r2"),
    ):
        update_comparison(counter, label, bool(row[f"{arm}_failed"]), nested)
    counter["r1_identity_mismatch"] += 1 - int(row["r1_exact_identity"])
    if not bool(row["current_r2_failed"]) and not nested:
        counter["current_nested_common_success"] += 1
        for metric in ("block_xors_milli", "block_muladds_milli"):
            counter[f"current_r2_{metric}_common"] += int(row[f"current_r2_{metric}"])
            counter[f"nested_r2_{metric}_common"] += int(row[f"nested_r2_{metric}"])


def ratio(numerator: int, denominator: int) -> str | None:
    return str(Decimal(numerator) / Decimal(denominator)) if denominator else None


def normalized_totals(counter: Counter[str]) -> dict[str, int]:
    return {key: counter[key] for key in TOTAL_KEYS}


def work_summary(totals: Counter[str]) -> dict[str, object]:
    return {
        "nested_over_current_r2_xor_common_success": ratio(
            totals["nested_r2_block_xors_milli_common"],
            totals["current_r2_block_xors_milli_common"],
        ),
        "nested_over_current_r2_muladd_common_success": ratio(
            totals["nested_r2_block_muladds_milli_common"],
            totals["current_r2_block_muladds_milli_common"],
        ),
        "common_success_cells": totals["current_nested_common_success"],
    }


def exact_p_values(totals: Counter[str]) -> dict[str, dict[str, object]]:
    return {
        label: binomial_one_sided(
            totals[f"nested_repair_{label}"],
            totals[f"nested_intro_{label}"],
        )
        for label in ("r0", "r1", "current_r2")
    }


def analysis_gates(
    totals: Counter[str], work: dict[str, object],
    p_values: dict[str, dict[str, object]],
    scopes: dict[tuple[str, str], Counter[str]],
) -> dict[str, bool]:
    errors = sum(totals[f"{arm}_error"] for arm in ARM_NAMES)
    expected_scope_keys = (
        {("schedule", value) for value in SCHEDULES} |
        {("loss", value) for value in LOSSES} |
        {("stratum", str(value)) for value in range(9)}
    )
    exact_scopes = set(scopes) == expected_scope_keys and all(
        counter["cells"] == (K_COUNT if scope_type == "stratum" else K_COUNT * 3)
        for (scope_type, _), counter in scopes.items()
    )
    return {
        "exact_cells": totals["cells"] == EXPECTED_CELLS,
        "exact_scope_partition": exact_scopes,
        "zero_errors": errors == 0,
        "exact_r1_identity": totals["r1_identity_mismatch"] == 0,
        "nested_better_than_r0": totals["nested_r2_fail"] < totals["r0_fail"],
        "nested_better_than_r1":
            totals["nested_r2_fail"] < totals["current_r1_fail"],
        "nested_better_than_current_r2":
            totals["nested_r2_fail"] < totals["current_r2_fail"],
        "nested_net_repairs_vs_r0":
            totals["nested_repair_r0"] > totals["nested_intro_r0"],
        "nested_net_repairs_vs_r1":
            totals["nested_repair_r1"] > totals["nested_intro_r1"],
        "nested_net_repairs_vs_current_r2":
            totals["nested_repair_current_r2"] >
            totals["nested_intro_current_r2"],
        "nested_exact_p_le_0_05_vs_r0":
            Decimal(str(p_values["r0"]["decimal"])) <= Decimal("0.05"),
        "nested_exact_p_le_0_05_vs_r1":
            Decimal(str(p_values["r1"]["decimal"])) <= Decimal("0.05"),
        "nested_exact_p_le_0_05_vs_current_r2":
            Decimal(str(p_values["current_r2"]["decimal"])) <= Decimal("0.05"),
        "no_adverse_schedule_slice_vs_current_r2": all(
            counter["nested_repair_current_r2"] >=
                counter["nested_intro_current_r2"]
            for (scope_type, _), counter in scopes.items()
            if scope_type == "schedule"
        ),
        "no_adverse_loss_slice_vs_current_r2": all(
            counter["nested_repair_current_r2"] >=
                counter["nested_intro_current_r2"]
            for (scope_type, _), counter in scopes.items()
            if scope_type == "loss"
        ),
        "no_adverse_stratum_vs_current_r2": all(
            counter["nested_repair_current_r2"] >=
                counter["nested_intro_current_r2"]
            for (scope_type, _), counter in scopes.items()
            if scope_type == "stratum"
        ),
        "xor_work_le_1_04":
            work["nested_over_current_r2_xor_common_success"] is not None and
            Decimal(str(work["nested_over_current_r2_xor_common_success"])) <=
            Decimal("1.04"),
        "muladd_work_le_1_001":
            work["nested_over_current_r2_muladd_common_success"] is not None and
            Decimal(str(work["nested_over_current_r2_muladd_common_success"])) <=
            Decimal("1.001"),
    }


def summary_text(totals: Counter[str], work: dict[str, object], passed: bool) -> str:
    lines = [
        f"passed={int(passed)}",
        f"cells={totals['cells']}",
        *(f"{arm}_failures={totals[f'{arm}_fail']}" for arm in ARM_NAMES),
        *(f"{arm}_errors={totals[f'{arm}_error']}" for arm in ARM_NAMES),
        *(
            f"{arm}_heavy_shortfalls={totals[f'{arm}_heavy_shortfall']}"
            for arm in ARM_NAMES
        ),
        f"r1_identity_mismatches={totals['r1_identity_mismatch']}",
        f"nested_vs_r0_repairs={totals['nested_repair_r0']}",
        f"nested_vs_r0_introductions={totals['nested_intro_r0']}",
        f"nested_vs_r1_repairs={totals['nested_repair_r1']}",
        f"nested_vs_r1_introductions={totals['nested_intro_r1']}",
        f"nested_vs_current_r2_repairs={totals['nested_repair_current_r2']}",
        f"nested_vs_current_r2_introductions={totals['nested_intro_current_r2']}",
        f"xor_nested_over_current_common="
        f"{work['nested_over_current_r2_xor_common_success']}",
        f"muladd_nested_over_current_common="
        f"{work['nested_over_current_r2_muladd_common_success']}",
        "timing_claim_valid=0",
    ]
    return "\n".join(lines) + "\n"


def analyze(args: argparse.Namespace) -> None:
    out = args.out.resolve()
    require_frozen_executor(out)
    _, strata, all_groups, all_jobs = verify_frozen(out)
    verify_results(out)
    analysis = out / "analysis"
    if analysis.exists() or analysis.is_symlink():
        common.die(f"analysis refuses existing {analysis}")
    analysis.mkdir()
    raw = out / "results/raw"
    job_map = {(job.stratum, job.group, job.arm): job for job in all_jobs}
    totals: Counter[str] = Counter()
    per_k = [Counter() for _ in range(K_LAST + 1)]
    scopes: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    cells_temp = analysis / "cells.csv.partial"
    with cells_temp.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(CELL_HEADER)
        for stratum in strata:
            for group_index, ks in enumerate(all_groups):
                parsed: dict[str, ParsedJob] = {}
                for arm in ARM_NAMES:
                    job = job_map[(stratum.index, group_index, arm)]
                    parsed[arm] = parse_job_output(
                        (raw / f"{job.stem}.stdout").read_bytes(), job, stratum, ks,
                    )
                for k in ks:
                    values = cell_values(k, stratum, parsed)
                    writer.writerow(values)
                    update_totals(totals, values)
                    update_totals(per_k[k], values)
                    for key in (
                        ("schedule", stratum.schedule),
                        ("loss", stratum.loss),
                        ("stratum", str(stratum.index)),
                    ):
                        update_totals(scopes[key], values)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(cells_temp, analysis / "cells.csv")
    if totals["cells"] != EXPECTED_CELLS:
        common.die(f"analyzed cells={totals['cells']}, want {EXPECTED_CELLS}")

    per_k_rows: list[list[object]] = [list(PER_K_HEADER)]
    exact_rows: list[list[object]] = [list(PER_K_HEADER)]
    for k in range(K_FIRST, K_LAST + 1):
        counter = per_k[k]
        row = [k] + [counter[name] for name in PER_K_HEADER[1:]]
        if counter["cells"] != 9:
            common.die(f"K={k} cells={counter['cells']}, want 9")
        per_k_rows.append(row)
        if any(counter[name] for name in (
            "nested_repair_r0", "nested_intro_r0", "nested_repair_r1",
            "nested_intro_r1", "nested_repair_current_r2",
            "nested_intro_current_r2",
        )):
            exact_rows.append(row)
    write_csv(analysis / "per_k.csv", per_k_rows)
    write_csv(analysis / "exact_k_changes.csv", exact_rows)
    scope_header = ("scope_type", "scope", *PER_K_HEADER[1:])
    write_csv(analysis / "scopes.csv", [scope_header] + [[
        scope_type, scope, *(counter[name] for name in PER_K_HEADER[1:])
    ] for (scope_type, scope), counter in sorted(scopes.items())])

    work = work_summary(totals)
    p_values = exact_p_values(totals)
    gates = analysis_gates(totals, work, p_values, scopes)
    expected_pass = all(gates.values())
    summary = {
        "schema": SCHEMA + ".analysis",
        "analyzed_utc": utc_now(), "passed": expected_pass,
        "gates": gates, "totals": normalized_totals(totals), "work": work,
        "one_sided_exact": p_values, "timing_claim_valid": False,
        "timing_note": "Saturated all-K confirmation; structural work only.",
    }
    common.atomic_json(analysis / "summary.json", summary)
    common.atomic_text(
        analysis / "summary.txt",
        summary_text(totals, work, bool(summary["passed"])),
    )
    analysis_paths = [
        analysis / "cells.csv", analysis / "per_k.csv",
        analysis / "exact_k_changes.csv", analysis / "scopes.csv",
        analysis / "summary.json", analysis / "summary.txt",
    ]
    seal = seal_files(out, analysis_paths, analysis / "analysis_complete.sha256")
    verify_analysis(out)
    print(json.dumps({
        "passed": summary["passed"], "totals": summary["totals"],
        "work": work, "one_sided_exact": p_values,
        "analysis_seal_sha256": seal,
    }, sort_keys=True))


def verify_csv_rows(
    path: Path, header: Sequence[str], rows: Iterable[Sequence[object]],
) -> None:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        if tuple(next(reader, ())) != tuple(header):
            common.die(f"{path}: header mismatch")
        for number, expected in enumerate(rows, 2):
            actual = next(reader, None)
            canonical = [str(value) for value in expected]
            if actual is None:
                common.die(f"{path}: ended before row {number}")
            if actual != canonical:
                common.die(f"{path}:{number}: derived row mismatch")
        if next(reader, None) is not None:
            common.die(f"{path}: contains extra rows")


def verify_analysis(out: Path) -> None:
    out = out.resolve()
    _, strata, all_groups, all_jobs = verify_frozen(out)
    verify_results(out)
    analysis = out / "analysis"
    hashes = verify_seal(out, analysis / "analysis_complete.sha256")
    expected = {
        (analysis / name).resolve()
        for name in (
            "cells.csv", "per_k.csv", "exact_k_changes.csv", "scopes.csv",
            "summary.json", "summary.txt",
        )
    }
    if set(hashes) != expected:
        common.die("analysis seal file set mismatch")
    enumerated_analysis = {
        common.regular_file(path, "analysis artifact")
        for path in analysis.iterdir()
    }
    if enumerated_analysis != expected | {
        (analysis / "analysis_complete.sha256").resolve()
    }:
        common.die("analysis directory contains an unsealed or missing artifact")
    raw = out / "results/raw"
    job_map = {(job.stratum, job.group, job.arm): job for job in all_jobs}
    totals: Counter[str] = Counter()
    per_k = [Counter() for _ in range(K_LAST + 1)]
    scopes: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    with (analysis / "cells.csv").open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        if tuple(next(reader, ())) != CELL_HEADER:
            common.die("analysis cells header mismatch")
        number = 1
        for stratum in strata:
            for group_index, ks in enumerate(all_groups):
                parsed: dict[str, ParsedJob] = {}
                for arm in ARM_NAMES:
                    job = job_map[(stratum.index, group_index, arm)]
                    parsed[arm] = parse_job_output(
                        (raw / f"{job.stem}.stdout").read_bytes(),
                        job, stratum, ks,
                    )
                for k in ks:
                    number += 1
                    expected_values = cell_values(k, stratum, parsed)
                    actual = next(reader, None)
                    if actual != [str(value) for value in expected_values]:
                        common.die(
                            f"analysis cells:{number}: sealed-raw reproduction mismatch"
                        )
                    update_totals(totals, expected_values)
                    update_totals(per_k[k], expected_values)
                    for key in (
                        ("schedule", stratum.schedule),
                        ("loss", stratum.loss),
                        ("stratum", str(stratum.index)),
                    ):
                        update_totals(scopes[key], expected_values)
        if next(reader, None) is not None:
            common.die("analysis cells contain rows beyond the exact raw domain")
    if totals["cells"] != EXPECTED_CELLS:
        common.die("independent analysis cell count mismatch")

    def per_k_rows(changes_only: bool) -> Iterable[Sequence[object]]:
        for k in range(K_FIRST, K_LAST + 1):
            counter = per_k[k]
            if counter["cells"] != 9:
                common.die(f"independent K={k} cell count mismatch")
            if changes_only and not any(counter[name] for name in (
                "nested_repair_r0", "nested_intro_r0", "nested_repair_r1",
                "nested_intro_r1", "nested_repair_current_r2",
                "nested_intro_current_r2",
            )):
                continue
            yield [k, *(counter[name] for name in PER_K_HEADER[1:])]

    verify_csv_rows(
        analysis / "per_k.csv", PER_K_HEADER, per_k_rows(False),
    )
    verify_csv_rows(
        analysis / "exact_k_changes.csv", PER_K_HEADER, per_k_rows(True),
    )
    scope_header = ("scope_type", "scope", *PER_K_HEADER[1:])
    verify_csv_rows(analysis / "scopes.csv", scope_header, (
        [scope_type, scope, *(counter[name] for name in PER_K_HEADER[1:])]
        for (scope_type, scope), counter in sorted(scopes.items())
    ))

    summary = json.loads((analysis / "summary.json").read_text(encoding="utf-8"))
    if not isinstance(summary, dict) or set(summary) != {
        "schema", "analyzed_utc", "passed", "gates", "totals", "work",
        "one_sided_exact", "timing_claim_valid", "timing_note",
    }:
        common.die("analysis summary key set mismatch")
    work = work_summary(totals)
    p_values = exact_p_values(totals)
    gates = analysis_gates(totals, work, p_values, scopes)
    expected_pass = all(gates.values())
    if (summary.get("schema") != SCHEMA + ".analysis" or
        summary.get("totals") != normalized_totals(totals) or
        summary.get("work") != work or
        summary.get("one_sided_exact") != p_values or
        summary.get("gates") != gates or
        summary.get("passed") is not expected_pass or
        summary.get("timing_claim_valid") is not False or
        summary.get("timing_note") !=
            "Saturated all-K confirmation; structural work only." or
        not isinstance(summary.get("analyzed_utc"), str) or
        not re.fullmatch(
            r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z",
            summary["analyzed_utc"],
        )):
        common.die("analysis summary independent reproduction mismatch")
    expected_text = summary_text(totals, work, expected_pass)
    if (analysis / "summary.txt").read_text(encoding="utf-8") != expected_text:
        common.die("analysis text summary independent reproduction mismatch")


def verify_command(args: argparse.Namespace) -> None:
    require_frozen_executor(args.out)
    verify_frozen(args.out)
    verify_results(args.out)
    verify_analysis(args.out)
    print(json.dumps({"verified": True, "output": str(args.out.resolve())}, sort_keys=True))


def self_test() -> None:
    commit = "1" * 40
    seal = "2" * 64
    first = make_strata(commit, seal)
    second = make_strata(commit, seal)
    if first != second or len(groups()) != GROUP_COUNT or len(jobs()) != EXPECTED_JOBS:
        common.die("self-test deterministic geometry")
    if binomial_one_sided(11, 0) != {
        "numerator": 1, "denominator": 2048,
        "decimal": str(Decimal(1) / Decimal(2048)),
    } or binomial_one_sided(1, 1) != {
        "numerator": 3, "denominator": 4,
        "decimal": "0.75",
    }:
        common.die("self-test exact binomial")
    print("self-test: PASS")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--out", type=Path, required=True)
    prepare_parser.add_argument("--screen-seal", type=Path, default=DEFAULT_SCREEN_SEAL)
    prepare_parser.add_argument("--current-binary", type=Path, default=DEFAULT_CURRENT_BINARY)
    prepare_parser.add_argument(
        "--nested-binary", type=Path,
        default=Path(__file__).resolve().parents[1] /
            "build-tail/codec/wirehair_v2_bench",
    )
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--out", type=Path, required=True)
    run_parser.add_argument("--thermal", type=Path, default=DEFAULT_THERMAL)
    run_parser.add_argument("--workers", type=int, default=120)
    run_parser.add_argument("--timeout-seconds", type=float, default=900)
    analyze_parser = subparsers.add_parser("analyze")
    analyze_parser.add_argument("--out", type=Path, required=True)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--out", type=Path, required=True)
    subparsers.add_parser("self-test")
    return result


def main(argv: Sequence[str]) -> int:
    args = parser().parse_args(argv)
    if args.command == "prepare":
        prepare(args)
    elif args.command == "run":
        run_campaign(args)
    elif args.command == "analyze":
        analyze(args)
    elif args.command == "verify":
        verify_command(args)
    else:
        self_test()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except common.ScreenError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
