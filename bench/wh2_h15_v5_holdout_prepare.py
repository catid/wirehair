#!/usr/bin/env python3
"""Prepare and seal the seed-isolated normalized-H15-v5 holdout campaign.

The real commitment is deliberately not a command-line input.  The only
post-freeze interface that can consume commitment bytes is ``reveal`` and it
reads stdin exactly once, after verifying that this harness commit is present
on a remote ref.  All other commands accept only the canonical three-row
holdout reveal snapshot.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from collections.abc import Callable, Iterable
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath


TABLE_SHA256 = "ef1729d961f71bc604e972194495e9e508459d62db786d0d9df56c3ce4198f3c"
SELECTION_COMPLETION_SHA256 = "4fdc34d4ff18bbbe5e037bd159bbb8c04d0e58e762cc065e5823ce5cf43495e8"
SOURCE_COMMITMENT_SHA256 = "f76eedbf4813597fa84af4b21a1236c957d11ff4a745585ea0442d9676645728"
PAYLOAD_K_SHA256 = "ae59d0cbeee5124b776b82ebaf107fe837d4a286946c8cbc270927fb8ec0035f"
PAYLOAD_PHASE_SEEDS = (
    "0x3321e92f548584b0", "0xa25642abac616007",
    "0x1b3a569b7011f225", "0xccb2312c22b50c26",
    "0xd8acf21ed37ad3a3", "0x989fa1566b3ae7bd",
    "0x46c04a5dceb47370", "0x206b83942cd89aab",
)
PAYLOAD_FIRST_ARMS = ("v4", "v5", "v5", "v4", "v5", "v4", "v4", "v5")
LOOKUP_FIRST_ARMS = ("v4", "v5", "v5", "v4") * 3
LOOKUP_CHECKSUMS = {
    "normalized-h15-v4": "0x5f244ad7ce774dfc",
    "normalized-h15-v5": "0x447d503555521b60",
}
SCHEDULES = ("iid", "burst", "permutation", "systematic-first",
             "repair-only", "adversarial")
LOSSES = ("0.10", "0.35", "0.50", "0.65")
ARMS = (("v4", "normalized-h15-v4"), ("v5", "normalized-h15-v5"))
STRUCTURAL_K = (
    2, 3, 63, 64, 65, 319, 320, 321, 944, 945, 946, 1279, 1280,
    1281, 3199, 3200, 3201, 9999, 10000, 10001, 19999, 20000,
    20001, 29999, 30000, 30001, 63999, 64000,
)
REVEAL_HEADER = ("partition", "index", "domain", "sha256", "derived_u64", "usage")
TABLE_HEADER = ("K", "v4_salt", "v5_salt", "changed", "decision")
GROUP_HEADER = ("group_id", "k_count", "k_csv")
PAYLOAD_K_HEADER = ("K", "class", "v4_salt", "v5_salt")
RECOVERY_JOB_HEADER = (
    "job_id", "stem", "phase", "group_id", "seed_index", "seed",
    "schedule", "loss", "arm", "profile", "trials", "k_count", "k_csv",
)
PAYLOAD_JOB_HEADER = (
    "job_id", "stem", "phase", "width", "arm", "profile", "seed",
    "first_arm", "overhead", "loss", "schedule", "trials", "k_count", "k_csv",
)
LOOKUP_JOB_HEADER = (
    "job_id", "stem", "phase", "cpu", "arm", "profile", "first_arm",
    "lookups", "permutation_step", "domain_repeats", "expected_checksum",
)
EXPECTED_RECOVERY_JOBS = 17424
EXPECTED_PAIRED_ROWS = 4720824
EXPECTED_RECOVERY_ARTIFACTS = 87123
EXPECTED_PAYLOAD_ARTIFACTS = 163
EXPECTED_LOOKUP_ARTIFACTS = 123
PRE_REVEAL_ROLES = (
    ("contract", "freeze/contract.json"),
    ("gates", "freeze/gates.json"),
    ("build_provenance", "freeze/build_provenance.tsv"),
    ("selection_table", "freeze/selection_table.tsv"),
    ("table_dump", "freeze/table_dump.tsv"),
    ("groups", "freeze/groups.tsv"),
    ("payload_k", "freeze/payload_k.tsv"),
    ("payload_jobs", "freeze/payload_jobs.tsv"),
    ("lookup_jobs", "freeze/lookup_jobs.tsv"),
    ("prepare", "freeze/tools/wh2_h15_v5_holdout_prepare.py"),
    ("recovery_job_runner", "freeze/tools/wh2_h15_v5_holdout_job.sh"),
    ("recovery_launcher", "freeze/tools/wh2_h15_v5_holdout_launch.sh"),
    ("analyzer", "freeze/tools/wh2_h15_v5_holdout_analyze.py"),
    ("performance", "freeze/tools/wh2_h15_v5_holdout_performance.py"),
    ("payload_artifact_manifest", "performance/payload_artifact_manifest.sha256"),
    ("payload_summary", "performance/payload_summary.json"),
    ("lookup_artifact_manifest", "performance/lookup_artifact_manifest.sha256"),
    ("lookup_summary", "performance/lookup_summary.json"),
    ("static_size_summary", "performance/static_size_summary.json"),
)
RECOVERY_INPUT_ROLES = (
    ("pre_reveal_seal", "meta/pre_reveal_seal.sha256"),
    ("holdout_reveal", "recovery/holdout_reveal.tsv"),
    ("reveal_provenance", "recovery/reveal_provenance.json"),
    ("recovery_manifest", "recovery/manifest.json"),
    ("recovery_jobs", "recovery/jobs.tsv"),
    ("frozen_binary", "freeze/wirehair_v2_bench"),
)
COMPLETION_ROLES = (
    ("pre_reveal_seal", "meta/pre_reveal_seal.sha256"),
    ("recovery_input_seal", "recovery/input_seal.sha256"),
    ("recovery_artifact_manifest", "recovery/artifact_manifest.sha256"),
    ("summary", "analysis/summary.json"),
    ("recovery_scopes", "analysis/recovery_scopes.csv"),
    ("per_k_recovery", "analysis/per_k_recovery.csv"),
    ("work_scopes", "analysis/work_scopes.csv"),
    ("semantic_identity", "analysis/semantic_identity.json"),
    ("environment", "analysis/environment.json"),
)
SEAL_HEADER = "role\tsha256\tpath\n"


class AuditError(RuntimeError):
    """A deterministic campaign-contract violation."""


def fail(message: str) -> None:
    raise AuditError(message)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_relative(text: str) -> str:
    if not text or "\\" in text or "\n" in text or "\r" in text:
        fail(f"unsafe relative path {text!r}")
    path = PurePosixPath(text)
    if path.is_absolute() or any(part in ("", ".", "..") for part in path.parts):
        fail(f"unsafe relative path {text!r}")
    normalized = path.as_posix()
    if normalized != text:
        fail(f"noncanonical relative path {text!r}")
    return normalized


PORTABLE_COMPONENT_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")


def canonical_remote_ref(text: str) -> str:
    prefix = "refs/remotes/"
    if not isinstance(text, str) or not text.startswith(prefix):
        fail(f"unsafe remote-tracking ref {text!r}")
    suffix = text[len(prefix):]
    components = suffix.split("/")
    if (len(components) < 2 or ".." in text or
            any(not PORTABLE_COMPONENT_RE.fullmatch(component) or
                component.endswith((".", ".lock")) for component in components)):
        fail(f"unsafe remote-tracking ref {text!r}")
    return text


def canonical_tree_path(text: str) -> str:
    relative = canonical_relative(text)
    if any(not PORTABLE_COMPONENT_RE.fullmatch(component) or component.endswith(".")
           for component in relative.split("/")):
        fail(f"unsafe portable Git tree path {text!r}")
    return relative


def _read_fd(fd: int, label: str) -> bytes:
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode):
        fail(f"{label}: expected regular file")
    chunks: list[bytes] = []
    while True:
        chunk = os.read(fd, 1024 * 1024)
        if not chunk:
            break
        chunks.append(chunk)
    after = os.fstat(fd)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    if identity(before) != identity(after) or sum(map(len, chunks)) != before.st_size:
        fail(f"{label}: file changed while it was read")
    return b"".join(chunks)


def require_no_symlink_components(path: Path, include_leaf: bool = True) -> Path:
    absolute = Path(os.path.abspath(path))
    parts = absolute.parts
    current = Path(parts[0])
    stop = len(parts) if include_leaf else len(parts) - 1
    for component in parts[1:stop]:
        current /= component
        st = current.lstat()
        if stat.S_ISLNK(st.st_mode):
            fail(f"{path}: symlink path component is forbidden: {current}")
    return absolute


def read_once(path: Path) -> bytes:
    require_no_symlink_components(path)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    try:
        return _read_fd(fd, str(path))
    finally:
        os.close(fd)


def open_root(root: Path) -> int:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(root, flags)
    st = os.fstat(fd)
    if st.st_uid != os.geteuid() or st.st_mode & 0o077:
        os.close(fd)
        fail(f"{root}: snapshot root must be owned by this user and mode 0700")
    return fd


def read_relative(root: Path, relative: str) -> bytes:
    parts = canonical_relative(relative).split("/")
    current = open_root(root)
    try:
        for component in parts[:-1]:
            flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
            if hasattr(os, "O_NOFOLLOW"):
                flags |= os.O_NOFOLLOW
            nxt = os.open(component, flags, dir_fd=current)
            st = os.fstat(nxt)
            if st.st_uid != os.geteuid() or st.st_mode & 0o077:
                os.close(nxt)
                fail(f"{root}/{relative}: snapshot directory is not private")
            os.close(current)
            current = nxt
        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        fd = os.open(parts[-1], flags, dir_fd=current)
        try:
            return _read_fd(fd, f"{root}/{relative}")
        finally:
            os.close(fd)
    finally:
        os.close(current)


def open_relative_directory(root: Path, relative: str) -> int:
    if relative == ".":
        return open_root(root)
    parts = canonical_relative(relative).split("/")
    current = open_root(root)
    try:
        for component in parts:
            flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY
            if hasattr(os, "O_NOFOLLOW"):
                flags |= os.O_NOFOLLOW
            nxt = os.open(component, flags, dir_fd=current)
            st = os.fstat(nxt)
            if st.st_uid != os.geteuid() or st.st_mode & 0o077:
                os.close(nxt)
                fail(f"{root}/{relative}: snapshot directory is not private")
            os.close(current)
            current = nxt
        return current
    except BaseException:
        os.close(current)
        raise


def verify_exact_directory_entries(
        root: Path, relative: str, expected_files: Iterable[str],
        expected_directories: Iterable[str] = ()) -> None:
    def names(values: Iterable[str], label: str) -> set[str]:
        result: set[str] = set()
        for value in values:
            if canonical_relative(value) != value or "/" in value or value in result:
                fail(f"{relative}: invalid/duplicate expected {label} entry {value!r}")
            result.add(value)
        return result

    files = names(expected_files, "file")
    directories = names(expected_directories, "directory")
    if files & directories:
        fail(f"{relative}: an entry cannot be both a file and directory")
    fd = open_relative_directory(root, relative)
    try:
        before = os.listdir(fd)
        actual = set(before)
        expected = files | directories
        if len(actual) != len(before) or actual != expected:
            missing = sorted(expected - actual)[:3]
            extra = sorted(actual - expected)[:3]
            fail(f"{relative}: exact entry set mismatch: missing={missing} extra={extra}")
        for name in files:
            st = os.stat(name, dir_fd=fd, follow_symlinks=False)
            if not stat.S_ISREG(st.st_mode):
                fail(f"{relative}/{name}: expected a regular nonsymlink file")
        for name in directories:
            st = os.stat(name, dir_fd=fd, follow_symlinks=False)
            if (not stat.S_ISDIR(st.st_mode) or st.st_uid != os.geteuid() or
                    st.st_mode & 0o077):
                fail(f"{relative}/{name}: expected a private nonsymlink directory")
        if set(os.listdir(fd)) != actual:
            fail(f"{relative}: entry set changed while it was enumerated")
    finally:
        os.close(fd)


def ensure_private_root(root: Path, fresh: bool) -> None:
    if fresh:
        require_no_symlink_components(root, include_leaf=False)
        root.mkdir(mode=0o700)
    require_no_symlink_components(root)
    st = root.lstat()
    if not stat.S_ISDIR(st.st_mode) or stat.S_ISLNK(st.st_mode):
        fail(f"{root}: campaign root must be a nonsymlink directory")
    if st.st_uid != os.geteuid() or st.st_mode & 0o077:
        fail(f"{root}: campaign root must be owned by this user and private")


def ensure_dir(root: Path, relative: str) -> None:
    current = root
    for component in canonical_relative(relative).split("/"):
        current = current / component
        try:
            current.mkdir(mode=0o700)
        except FileExistsError:
            st = current.lstat()
            if (not stat.S_ISDIR(st.st_mode) or stat.S_ISLNK(st.st_mode) or
                    st.st_uid != os.geteuid() or st.st_mode & 0o077):
                fail(f"{current}: expected nonsymlink directory")


def write_exclusive(root: Path, relative: str, data: bytes, executable: bool = False) -> None:
    relative = canonical_relative(relative)
    parent = str(PurePosixPath(relative).parent)
    if parent != ".":
        ensure_dir(root, parent)
    path = root / relative
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags, 0o700 if executable else 0o600)
    try:
        offset = 0
        while offset < len(data):
            offset += os.write(fd, data[offset:])
        os.fsync(fd)
    finally:
        os.close(fd)
    if read_relative(root, relative) != data:
        fail(f"{path}: exclusive write did not preserve exact bytes")


def copy_snapshot(root: Path, relative: str, source: Path, executable: bool = False) -> str:
    data = read_once(source)
    write_exclusive(root, relative, data, executable)
    return sha256_bytes(data)


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("ascii")


def parse_json(data: bytes, label: str) -> object:
    try:
        text = data.decode("ascii")
        value = json.loads(text, object_pairs_hook=_unique_object,
                           parse_constant=lambda token: fail(
                               f"{label}: nonfinite JSON constant {token!r}"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        fail(f"{label}: invalid JSON: {exc}")
    if canonical_json(value) != data:
        fail(f"{label}: JSON is not canonical")
    return value


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def canonical_tsv(header: tuple[str, ...], rows: Iterable[Iterable[object]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    writer.writerow(header)
    writer.writerows(rows)
    data = output.getvalue().encode("ascii")
    if b'"' in data or b"\r" in data:
        fail("TSV serialization unexpectedly required quoting")
    return data


def parse_tsv(data: bytes, header: tuple[str, ...], label: str) -> list[dict[str, str]]:
    if b"\r" in data or not data.endswith(b"\n") or b"\n\n" in data:
        fail(f"{label}: requires LF-only framing, final LF, and no blank rows")
    try:
        text = data.decode("ascii")
    except UnicodeError:
        fail(f"{label}: expected ASCII")
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if tuple(reader.fieldnames or ()) != header:
        fail(f"{label}: header mismatch")
    rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        fail(f"{label}: malformed TSV row")
    if canonical_tsv(header, ([row[field] for field in header] for row in rows)) != data:
        fail(f"{label}: noncanonical TSV serialization")
    return rows


def uint(text: str, label: str, maximum: int | None = None) -> int:
    if not re.fullmatch(r"(?:0|[1-9][0-9]*)", text):
        fail(f"{label}: noncanonical unsigned integer {text!r}")
    value = int(text)
    if maximum is not None and value > maximum:
        fail(f"{label}: exceeds {maximum}")
    return value


def salt(text: str, label: str) -> int:
    if not re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]*)", text):
        fail(f"{label}: noncanonical salt {text!r}")
    value = int(text, 16)
    if value > 255:
        fail(f"{label}: salt exceeds uint8")
    return value


def parse_table(data: bytes, require_pinned_hash: bool = True) -> dict[int, dict[str, str]]:
    if require_pinned_hash and sha256_bytes(data) != TABLE_SHA256:
        fail("selection table does not match the independently sealed SHA-256")
    parsed = parse_tsv(data, TABLE_HEADER, "selection table")
    result: dict[int, dict[str, str]] = {}
    counts = {0x27: 0, 0x79: 0, 0x6F: 0, 0xA8: 0}
    inherited = 0
    changed = 0
    for line, row in enumerate(parsed, 2):
        k = uint(row["K"], f"selection table:{line}:K", 64000)
        v4 = salt(row["v4_salt"], f"selection table:{line}:v4")
        v5 = salt(row["v5_salt"], f"selection table:{line}:v5")
        bit = uint(row["changed"], f"selection table:{line}:changed", 1)
        if k < 2 or k in result or bit != (v4 != v5):
            fail(f"selection table:{line}: K/order/change mismatch")
        if v4:
            inherited += 1
            if v5 != v4 or row["decision"] != "immutable-nonzero-v4":
                fail(f"selection table:{line}: changed inherited v4 entry")
        elif v5:
            changed += 1
            if v5 not in counts or row["decision"] != "selected":
                fail(f"selection table:{line}: unapproved selected salt")
            counts[v5] += 1
        elif row["decision"] != "retain-zero":
            fail(f"selection table:{line}: zero decision mismatch")
        result[k] = row
    if (list(result) != list(range(2, 64001)) or inherited != 70 or changed != 98 or
            counts != {0x27: 75, 0x79: 17, 0x6F: 6, 0xA8: 0}):
        fail("selection table invariant/count mismatch")
    return result


def make_groups() -> tuple[bytes, list[tuple[int, tuple[int, ...]]]]:
    groups = [(group, tuple(range(2 + group, 64001, 120))) for group in range(120)]
    if {k for _, ks in groups for k in ks} != set(range(2, 64001)):
        fail("stride groups do not partition K=2..64000")
    data = canonical_tsv(GROUP_HEADER, (
        (group, len(ks), ",".join(map(str, ks))) for group, ks in groups
    ))
    return data, groups


def make_payload_k(table: dict[int, dict[str, str]], enforce_hash: bool = True) -> bytes:
    nonzero = {k for k, row in table.items() if row["v5_salt"] != "0x0"}
    structural = set(STRUCTURAL_K)
    if len(nonzero) != 168 or nonzero & structural:
        fail("payload nonzero/structural support contract mismatch")
    selected = nonzero | structural
    eligible = set(table) - selected
    distances = {k: min(abs(k - chosen) for chosen in selected) for k in eligible}
    controls: list[int] = []
    for _ in range(60):
        chosen = max(eligible, key=lambda k: (distances[k], -k))
        controls.append(chosen)
        eligible.remove(chosen)
        for k in eligible:
            distances[k] = min(distances[k], abs(k - chosen))
    classes = ({k: ("changed-v5" if table[k]["changed"] == "1" else "immutable-v4")
                for k in nonzero} |
               {k: "structural-control" for k in structural} |
               {k: "maximin-control" for k in controls})
    data = canonical_tsv(PAYLOAD_K_HEADER, (
        (k, classes[k], table[k]["v4_salt"], table[k]["v5_salt"])
        for k in sorted(classes)
    ))
    if len(classes) != 256:
        fail("payload panel must contain exactly 256 unique K")
    if enforce_hash and sha256_bytes(data) != PAYLOAD_K_SHA256:
        fail("deterministic payload panel SHA-256 mismatch")
    return data


def make_payload_jobs(payload_data: bytes) -> bytes:
    panel = parse_tsv(payload_data, PAYLOAD_K_HEADER, "payload K panel")
    ks = tuple(uint(row["K"], "payload K") for row in panel)
    k_csv = ",".join(map(str, ks))
    rows: list[tuple[object, ...]] = []
    for phase, (seed, first) in enumerate(zip(PAYLOAD_PHASE_SEEDS, PAYLOAD_FIRST_ARMS)):
        widths = (1280, 4096) if phase % 2 == 0 else (4096, 1280)
        arms = ARMS if first == "v4" else tuple(reversed(ARMS))
        for width in widths:
            for arm, profile in arms:
                job_id = len(rows)
                rows.append((job_id, f"payload{job_id:02d}_p{phase}_{width}_{arm}", phase,
                             width, arm, profile, seed, first, 4, "0.10", "repair-only",
                             1, len(ks), k_csv))
    if len(rows) != 32:
        fail("payload ledger must contain exactly 32 jobs")
    return canonical_tsv(PAYLOAD_JOB_HEADER, rows)


def make_lookup_jobs() -> bytes:
    rows: list[tuple[object, ...]] = []
    for phase, first in enumerate(LOOKUP_FIRST_ARMS):
        arms = ARMS if first == "v4" else tuple(reversed(ARMS))
        for arm, profile in arms:
            job_id = len(rows)
            rows.append((job_id, f"lookup{job_id:02d}_p{phase}_{arm}", phase,
                         0 if phase % 2 == 0 else 64, arm, profile, first,
                         63999 * 2048, 40501, 2048, LOOKUP_CHECKSUMS[profile]))
    if len(rows) != 24:
        fail("lookup ledger must contain exactly 24 jobs")
    return canonical_tsv(LOOKUP_JOB_HEADER, rows)


def make_gates() -> dict[str, object]:
    return {
        "environment": {
            "all_dimm_blackout_forbidden": True,
            "cpu_busy_mean_min_pct": 98.0,
            "cpu_busy_sample_min_pct": 90.0,
            "dimm_coverage_min_fraction": 0.99,
            "edac_delta_max": 0,
            "endpoint_slack_s": 2.5,
            "live_abort_consecutive_hot_samples": 3,
            "live_abort_no_complete_sample_s": 5.0,
            "max_consecutive_dimm_misses": 5,
            "max_sample_gap_s": 2.5,
            "temperature_abort_c": 90.0,
        },
        "lookup": {
            "absolute_median_delta_ns_max": 2.0,
            "lookup_specific_code_data_bytes_max": 4096,
            "logical_table_bytes_max": 2048,
            "median_ratio_max": 1.25,
            "per_phase_ratio_max": 1.50,
            "phases": 12,
            "reference_checksum_identity": True,
            "whole_binary_section_deltas_report_only": True,
        },
        "payload": {
            "all_256_solve_ratio_max": 1.01,
            "changed_98_combined_ratio_max": 1.02,
            "errors_max": 0,
            "median_phase_ratio_max": 1.01,
            "muladd_ratio_max": 1.0025,
            "paired_both_success_only": True,
            "paired_success_min_per_phase_width": 231,
            "v5_only_rank_failures_max": 0,
            "width_solve_ratio_max": 1.02,
            "xor_ratio_max": 1.01,
        },
        "recovery": {
            "changed_baseline_failures_min": 40,
            "changed_failure_reduction_min_fraction": 0.25,
            "changed_improvement_p_max": 0.01,
            "errors_max": 0,
            "holm_alpha": 0.05,
            "holm_changed_k_family_size": 98,
            "holm_scope_family_size": 57,
            "material_min_discordant_delta": 3,
            "material_paired_risk_increase_min": 0.0005,
            "material_v5_failure_ratio_min": 1.10,
            "muladd_ratio_max": 1.0025,
            "overall_v5_failures_not_above_v4": True,
            "semantic_identity_unchanged_k": True,
            "xor_ratio_max": 1.01,
        },
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.gates.v1",
    }


def build_provenance(binary_data: bytes, reference_data: bytes, source_commit: str) -> bytes:
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        fail("source commit must be exactly 40 lowercase hexadecimal digits")
    return canonical_tsv(("key", "value"), (
        ("binary_sha256", sha256_bytes(binary_data)),
        ("reference_binary_sha256", sha256_bytes(reference_data)),
        ("selection_completion_seal_sha256", SELECTION_COMPLETION_SHA256),
        ("selection_table_sha256", TABLE_SHA256),
        ("source_commit", source_commit),
    ))


def prepare_freeze(args: argparse.Namespace) -> None:
    root = args.output
    ensure_private_root(root, fresh=True)
    table_data = read_once(args.selection_table)
    dump_data = read_once(args.table_dump)
    if table_data != dump_data:
        fail("selection table and frozen binary table dump differ byte-for-byte")
    table = parse_table(table_data)
    binary_data = read_once(args.binary)
    reference_data = read_once(args.reference_binary)
    groups_data, _ = make_groups()
    payload_data = make_payload_k(table)
    payload_jobs = make_payload_jobs(payload_data)
    lookup_jobs = make_lookup_jobs()
    tool_specs = (
        ("prepare", "freeze/tools/wh2_h15_v5_holdout_prepare.py", args.prepare, True),
        ("recovery_job_runner", "freeze/tools/wh2_h15_v5_holdout_job.sh", args.job_runner, True),
        ("recovery_launcher", "freeze/tools/wh2_h15_v5_holdout_launch.sh", args.launcher, True),
        ("analyzer", "freeze/tools/wh2_h15_v5_holdout_analyze.py", args.analyzer, True),
        ("performance", "freeze/tools/wh2_h15_v5_holdout_performance.py", args.performance, True),
    )
    tool_data = {role: read_once(path) for role, _, path, _ in tool_specs}
    tool_hashes = {role: sha256_bytes(data) for role, data in tool_data.items()}
    contract = {
        "binary_sha256": sha256_bytes(binary_data),
        "completion_roles": [role for role, _ in COMPLETION_ROLES],
        "geometry": {
            "breadth_invocations": 17280,
            "breadth_paired_rows": 4607928,
            "changed_k": 98,
            "depth_invocations": 144,
            "depth_paired_rows": 112896,
            "holdout_rows": 3,
            "losses": list(LOSSES),
            "paired_rows": EXPECTED_PAIRED_ROWS,
            "recovery_invocations": EXPECTED_RECOVERY_JOBS,
            "schedules": list(SCHEDULES),
            "stride_groups": 120,
        },
        "lookup_artifacts": EXPECTED_LOOKUP_ARTIFACTS,
        "payload_artifacts": EXPECTED_PAYLOAD_ARTIFACTS,
        "payload_k_sha256": PAYLOAD_K_SHA256,
        "pre_reveal_roles": [role for role, _ in PRE_REVEAL_ROLES],
        "recovery_artifacts": EXPECTED_RECOVERY_ARTIFACTS,
        "recovery_input_roles": [role for role, _ in RECOVERY_INPUT_ROLES],
        "reference_binary_sha256": sha256_bytes(reference_data),
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.contract.v1",
        "selection_completion_seal_sha256": SELECTION_COMPLETION_SHA256,
        "selection_table_sha256": TABLE_SHA256,
        "source_commitment_sha256": SOURCE_COMMITMENT_SHA256,
        "source_commitment_usage": "post-push stdin-only one-read provenance; never a path input",
        "source_revision": args.source_commit,
        "tool_sha256": tool_hashes,
    }
    outputs = (
        ("freeze/contract.json", canonical_json(contract), False),
        ("freeze/gates.json", canonical_json(make_gates()), False),
        ("freeze/build_provenance.tsv", build_provenance(binary_data, reference_data, args.source_commit), False),
        ("freeze/selection_table.tsv", table_data, False),
        ("freeze/table_dump.tsv", dump_data, False),
        ("freeze/groups.tsv", groups_data, False),
        ("freeze/payload_k.tsv", payload_data, False),
        ("freeze/payload_jobs.tsv", payload_jobs, False),
        ("freeze/lookup_jobs.tsv", lookup_jobs, False),
        ("freeze/wirehair_v2_bench", binary_data, True),
        ("freeze/reference_wirehair_v2_bench", reference_data, True),
    )
    for relative, data, executable in outputs:
        write_exclusive(root, relative, data, executable)
    for role, relative, _source, executable in tool_specs:
        write_exclusive(root, relative, tool_data[role], executable)
    ensure_dir(root, "performance")
    ensure_dir(root, "meta")
    print("prepared seed-free holdout freeze inputs (performance and pre-reveal seal pending)")


def parse_reveal(data: bytes, partition: str = "holdout") -> list[dict[str, str]]:
    rows = parse_tsv(data, REVEAL_HEADER, f"{partition} reveal")
    seen_digest: set[str] = set()
    seen_seed: set[str] = set()
    ordered: dict[int, dict[str, str]] = {}
    for line, row in enumerate(rows, 2):
        if row["partition"] != partition:
            fail(f"reveal:{line}: expected only partition {partition!r}")
        index = uint(row["index"], f"reveal:{line}:index", 2)
        if index in ordered:
            fail(f"reveal:{line}: duplicate index")
        if not re.fullmatch(r"[\x20-\x7e]+", row["domain"]) or not re.fullmatch(r"[\x20-\x7e]+", row["usage"]):
            fail(f"reveal:{line}: domain/usage must be canonical printable ASCII")
        lowered = (row["domain"] + " " + row["usage"]).lower()
        if partition not in lowered or (partition == "holdout" and "selection" in lowered):
            fail(f"reveal:{line}: domain/usage partition mismatch")
        digest = row["sha256"]
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            fail(f"reveal:{line}: invalid SHA-256")
        derived = "0x" + digest[:16]
        if row["derived_u64"] != derived or digest in seen_digest or derived in seen_seed:
            fail(f"reveal:{line}: derived seed or uniqueness mismatch")
        seen_digest.add(digest)
        seen_seed.add(derived)
        ordered[index] = row
    if set(ordered) != {0, 1, 2}:
        fail("reveal must contain exact indices 0..2")
    result = [ordered[index] for index in range(3)]
    if canonical_tsv(REVEAL_HEADER, ([row[key] for key in REVEAL_HEADER] for row in result)) != data:
        fail("reveal row order is not canonical 0..2")
    return result


def parse_whole_commitment_once(data: bytes, expected_sha: str) -> bytes:
    if sha256_bytes(data) != expected_sha:
        fail("whole commitment SHA-256 differs from the frozen procedural provenance")
    rows = parse_tsv(data, REVEAL_HEADER, "whole commitment")
    holdout = [row for row in rows if row["partition"] == "holdout"]
    # Validate only the emitted subset through the same canonical parser.  No
    # selection value is written, printed, returned, or included in provenance.
    subset = canonical_tsv(REVEAL_HEADER, ([row[key] for key in REVEAL_HEADER] for row in holdout))
    parse_reveal(subset)
    if len(rows) != 6:
        fail("whole commitment must contain exactly three selection and three holdout rows")
    return subset


def validate_seal_commit_state(
        source_revision: str, source_resolved: str | None, head: str,
        parents: tuple[str, ...], remote_ref: str, upstream_ref: str,
        tracking_commit: str, live_commit: str, live_ref: str,
        expected_live_ref: str, changes: tuple[tuple[str, str], ...],
        sealed_tree_path: str, sealed_mode: str) -> None:
    commit_pattern = r"[0-9a-f]{40}"
    canonical_remote_ref(remote_ref)
    canonical_tree_path(sealed_tree_path)
    if not re.fullmatch(commit_pattern, source_revision):
        fail("frozen source revision is not an exact commit ID")
    if source_resolved is None or not re.fullmatch(commit_pattern, source_resolved):
        fail("frozen source revision does not resolve to an existing commit")
    if source_resolved != source_revision:
        fail("frozen source revision did not resolve to its exact commit ID")
    if not re.fullmatch(commit_pattern, head):
        fail("executing repository HEAD is not an exact commit ID")
    if parents != (source_revision,):
        fail("final HEAD must be exactly one seal-only commit atop the frozen source revision")
    if upstream_ref != remote_ref:
        fail("configured upstream differs from the required remote-tracking ref")
    if tracking_commit != head:
        fail("local remote-tracking ref is ahead of or behind final HEAD")
    if live_commit != head or live_ref != expected_live_ref:
        fail("live remote branch tip differs from final HEAD")
    if changes != (("A", sealed_tree_path),):
        fail("seal commit must add exactly the canonical sealed-tree path and no other change")
    if sealed_mode != "100644":
        fail("committed pre-reveal seal must be one non-executable regular file")


def verify_git_push(repo: Path, remote_ref: str, sealed_tree_path: str,
                    local_seal: bytes, source_revision: str) -> tuple[str, str]:
    repo = require_no_symlink_components(repo)
    remote_ref = canonical_remote_ref(remote_ref)
    tree_path = canonical_tree_path(sealed_tree_path)
    if not re.fullmatch(r"[0-9a-f]{40}", source_revision):
        fail("source revision must be an exact lowercase SHA-1 commit ID")
    git_environment = {key: value for key, value in os.environ.items()
                       if not key.startswith("GIT_")}
    git_environment.update({
        "GIT_NO_REPLACE_OBJECTS": "1",
        "GIT_SSH_COMMAND": "ssh -o BatchMode=yes",
        "GIT_TERMINAL_PROMPT": "0",
    })
    def git(*arguments: str) -> str:
        result = subprocess.run(
            ("git", "--no-replace-objects", "-C", str(repo), *arguments), check=False,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=git_environment)
        if result.returncode != 0:
            fail(f"git {' '.join(arguments)} failed: {result.stderr.strip()}")
        return result.stdout.strip()

    def raw_commit_parents(commit: str) -> tuple[str, ...]:
        if git("cat-file", "-t", commit) != "commit":
            fail(f"{commit}: object is not a raw commit")
        raw = git("cat-file", "commit", commit)
        headers, separator, _message = raw.partition("\n\n")
        lines = headers.splitlines()
        if (not separator or not lines or
                not re.fullmatch(r"tree [0-9a-f]{40}", lines[0])):
            fail(f"{commit}: malformed raw commit headers")
        parents: list[str] = []
        for line in lines[1:]:
            if line.startswith("parent "):
                if not re.fullmatch(r"parent [0-9a-f]{40}", line):
                    fail(f"{commit}: malformed raw parent header")
                parents.append(line[7:])
            elif parents:
                break
        return tuple(parents)

    if (git("rev-parse", "--is-inside-work-tree") != "true" or
            Path(git("rev-parse", "--show-toplevel")) != repo):
        fail("Git repository selection does not resolve to the requested worktree")
    if git("status", "--porcelain", "--untracked-files=no"):
        fail("revealer requires a clean tracked worktree")
    head = git("rev-parse", "--verify", "HEAD^{commit}")
    if git("cat-file", "-t", source_revision) != "commit":
        fail("frozen source revision is not a raw commit object")
    source_resolved = source_revision
    parents = raw_commit_parents(head)
    upstream = git("rev-parse", "--symbolic-full-name", "@{upstream}")
    tracking = git("rev-parse", "--verify", f"{remote_ref}^{{commit}}")
    match = re.fullmatch(r"refs/remotes/([^/]+)/(.+)", remote_ref)
    if not match:
        fail("remote ref cannot be mapped to a remote/branch for live verification")
    remote_name, branch_name = match.groups()
    expected_live_ref = f"refs/heads/{branch_name}"
    live = subprocess.run(("git", "--no-replace-objects", "-C", str(repo),
                           "ls-remote", "--exit-code", remote_name, expected_live_ref), check=False,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                          env=git_environment)
    live_lines = live.stdout.splitlines()
    if live.returncode != 0 or len(live_lines) != 1:
        fail("live Git remote branch could not be resolved uniquely")
    live_fields = live_lines[0].split("\t")
    if len(live_fields) != 2:
        fail("live Git remote returned a malformed branch record")
    diff = subprocess.run(
        ("git", "--no-replace-objects", "-C", str(repo), "diff-tree",
         "--no-commit-id", "--name-status",
         "-r", "-z", source_revision, head), check=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=git_environment)
    if diff.returncode != 0:
        fail(f"could not inspect seal-only commit: {diff.stderr.decode('utf-8', 'replace').strip()}")
    diff_fields = diff.stdout.rstrip(b"\0").split(b"\0") if diff.stdout else []
    if len(diff_fields) % 2:
        fail("seal-only commit returned a malformed changed-path record")
    try:
        changes = tuple((diff_fields[index].decode("ascii"),
                         diff_fields[index + 1].decode("ascii"))
                        for index in range(0, len(diff_fields), 2))
    except UnicodeError:
        fail("seal-only commit changed a non-ASCII path")
    tree = subprocess.run(
        ("git", "--no-replace-objects", "-C", str(repo), "ls-tree", "-z",
         head, "--", tree_path), check=False, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env=git_environment)
    if tree.returncode != 0:
        fail(f"could not inspect committed seal mode: {tree.stderr.decode('utf-8', 'replace').strip()}")
    tree_match = re.fullmatch(
        rb"(100[0-7]{3}) blob [0-9a-f]{40}\t" +
        re.escape(tree_path.encode("ascii")) + rb"\x00", tree.stdout)
    sealed_mode = tree_match.group(1).decode("ascii") if tree_match else ""
    validate_seal_commit_state(
        source_revision, source_resolved, head, parents, remote_ref, upstream,
        tracking, live_fields[0], live_fields[1], expected_live_ref, changes,
        tree_path, sealed_mode)
    blob = subprocess.run(
        ("git", "--no-replace-objects", "-C", str(repo), "show", f"{head}:{tree_path}"),
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=git_environment)
    if blob.returncode != 0 or blob.stdout != local_seal:
        fail("the exact pre-reveal seal is not committed in the pushed harness revision")
    return head, live_fields[0]


def seal_data(root: Path, roles: tuple[tuple[str, str], ...]) -> bytes:
    rows = []
    for role, relative in roles:
        data = read_relative(root, relative)
        rows.append((role, sha256_bytes(data), canonical_relative(relative)))
    return canonical_tsv(("role", "sha256", "path"), rows)


def verify_seal(root: Path, seal_relative: str,
                roles: tuple[tuple[str, str], ...]) -> dict[str, str]:
    data = read_relative(root, seal_relative)
    rows = parse_tsv(data, ("role", "sha256", "path"), seal_relative)
    expected = [(role, path) for role, path in roles]
    actual = [(row["role"], row["path"]) for row in rows]
    if actual != expected:
        fail(f"{seal_relative}: role/path order mismatch")
    result: dict[str, str] = {}
    for row in rows:
        if not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
            fail(f"{seal_relative}: invalid digest")
        actual_digest = sha256_bytes(read_relative(root, row["path"]))
        if actual_digest != row["sha256"]:
            fail(f"{seal_relative}: {row['role']} digest mismatch")
        result[row["role"]] = row["sha256"]
    if seal_data(root, roles) != data:
        fail(f"{seal_relative}: noncanonical or stale seal")
    return result


def verify_hash_manifest(root: Path, relative: str, expected_count: int) -> None:
    rows = parse_tsv(read_relative(root, relative), ("sha256", "path"), relative)
    if len(rows) != expected_count:
        fail(f"{relative}: expected exactly {expected_count} records")
    paths: set[str] = set()
    for row in rows:
        path = canonical_relative(row["path"])
        if path in paths or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
            fail(f"{relative}: duplicate path or malformed digest")
        paths.add(path)
        if sha256_bytes(read_relative(root, path)) != row["sha256"]:
            fail(f"{relative}: artifact digest mismatch for {path}")


def validate_freeze(root: Path) -> dict[str, object]:
    contract_data = read_relative(root, "freeze/contract.json")
    manifest = parse_json(contract_data, "freeze contract")
    if not isinstance(manifest, dict):
        fail("freeze contract root must be an object")
    expected_keys = {
        "binary_sha256", "completion_roles", "geometry", "lookup_artifacts",
        "payload_artifacts", "payload_k_sha256", "pre_reveal_roles",
        "recovery_artifacts", "recovery_input_roles", "reference_binary_sha256",
        "schema", "selection_completion_seal_sha256", "selection_table_sha256",
        "source_commitment_sha256", "source_commitment_usage", "source_revision",
        "tool_sha256",
    }
    expected_geometry = {
        "breadth_invocations": 17280, "breadth_paired_rows": 4607928,
        "changed_k": 98, "depth_invocations": 144, "depth_paired_rows": 112896,
        "holdout_rows": 3, "losses": list(LOSSES), "paired_rows": EXPECTED_PAIRED_ROWS,
        "recovery_invocations": EXPECTED_RECOVERY_JOBS, "schedules": list(SCHEDULES),
        "stride_groups": 120,
    }
    if (set(manifest) != expected_keys or
            manifest.get("schema") != "wirehair.wh2.normalized_h15_v5.holdout.contract.v1" or
            manifest.get("geometry") != expected_geometry or
            manifest.get("payload_artifacts") != EXPECTED_PAYLOAD_ARTIFACTS or
            manifest.get("lookup_artifacts") != EXPECTED_LOOKUP_ARTIFACTS or
            manifest.get("recovery_artifacts") != EXPECTED_RECOVERY_ARTIFACTS or
            manifest.get("payload_k_sha256") != PAYLOAD_K_SHA256 or
            manifest.get("selection_table_sha256") != TABLE_SHA256 or
            manifest.get("selection_completion_seal_sha256") != SELECTION_COMPLETION_SHA256 or
            manifest.get("source_commitment_sha256") != SOURCE_COMMITMENT_SHA256 or
            manifest.get("source_commitment_usage") !=
            "post-push stdin-only one-read provenance; never a path input" or
            manifest.get("pre_reveal_roles") != [role for role, _ in PRE_REVEAL_ROLES] or
            manifest.get("recovery_input_roles") != [role for role, _ in RECOVERY_INPUT_ROLES] or
            manifest.get("completion_roles") != [role for role, _ in COMPLETION_ROLES] or
            not re.fullmatch(r"[0-9a-f]{40}", str(manifest.get("source_revision", "")))):
        fail("freeze contract fields differ from the exact campaign contract")
    binary = read_relative(root, "freeze/wirehair_v2_bench")
    reference = read_relative(root, "freeze/reference_wirehair_v2_bench")
    if (sha256_bytes(binary) != manifest.get("binary_sha256") or
            sha256_bytes(reference) != manifest.get("reference_binary_sha256")):
        fail("frozen candidate/reference binary hash mismatch")
    table_data = read_relative(root, "freeze/selection_table.tsv")
    dump_data = read_relative(root, "freeze/table_dump.tsv")
    if table_data != dump_data:
        fail("frozen table dump differs from the independently sealed selection table")
    table = parse_table(table_data)
    groups_data, _ = make_groups()
    if read_relative(root, "freeze/groups.tsv") != groups_data:
        fail("frozen groups differ from exact stride reconstruction")
    payload = make_payload_k(table)
    if read_relative(root, "freeze/payload_k.tsv") != payload:
        fail("frozen payload panel differs from deterministic reconstruction")
    if read_relative(root, "freeze/payload_jobs.tsv") != make_payload_jobs(payload):
        fail("frozen payload jobs differ from exact reconstruction")
    if read_relative(root, "freeze/lookup_jobs.tsv") != make_lookup_jobs():
        fail("frozen lookup jobs differ from exact reconstruction")
    if read_relative(root, "freeze/gates.json") != canonical_json(make_gates()):
        fail("frozen gate file differs from exact predeclared thresholds")
    tool_hashes = manifest.get("tool_sha256")
    expected_tool_roles = {role for role, _ in PRE_REVEAL_ROLES[9:14]}
    if not isinstance(tool_hashes, dict) or set(tool_hashes) != expected_tool_roles:
        fail("freeze contract tool hash map is malformed")
    for role, relative in PRE_REVEAL_ROLES[9:14]:
        if tool_hashes.get(role) != sha256_bytes(read_relative(root, relative)):
            fail(f"frozen tool hash mismatch: {role}")
    expected_provenance = build_provenance(binary, reference, str(manifest["source_revision"]))
    if read_relative(root, "freeze/build_provenance.tsv") != expected_provenance:
        fail("build provenance differs from frozen binary/revision hashes")
    return manifest


def expected_performance_paths(root: Path, panel: str) -> list[str]:
    if panel == "payload":
        header = PAYLOAD_JOB_HEADER
        jobs_path = "freeze/payload_jobs.tsv"
        expected_jobs = 32
        expected_artifacts = EXPECTED_PAYLOAD_ARTIFACTS
    elif panel == "lookup":
        header = LOOKUP_JOB_HEADER
        jobs_path = "freeze/lookup_jobs.tsv"
        expected_jobs = 24
        expected_artifacts = EXPECTED_LOOKUP_ARTIFACTS
    else:
        fail(f"unknown performance panel {panel!r}")
    jobs = parse_tsv(read_relative(root, jobs_path), header, jobs_path)
    if len(jobs) != expected_jobs:
        fail(f"{panel}: performance job count mismatch")
    paths: list[str] = []
    for row in jobs:
        stem = row["stem"]
        paths.extend((f"performance/{panel}/raw/{stem}.stdout",
                      f"performance/{panel}/raw/{stem}.stderr",
                      f"performance/{panel}/raw/{stem}.time",
                      f"performance/{panel}/status/{stem}.ok",
                      f"performance/{panel}/commands/{stem}.txt"))
    paths.extend((f"performance/{panel}/run_start.json",
                  f"performance/{panel}/run_finish.json",
                  f"performance/{panel}/thermal_interval.csv"))
    if len(paths) != expected_artifacts or len(set(paths)) != len(paths):
        fail(f"{panel}: performance artifact path contract mismatch")
    return paths


def verify_performance_artifact_sets(root: Path) -> None:
    for panel in ("payload", "lookup"):
        paths = expected_performance_paths(root, panel)
        panel_root = f"performance/{panel}"
        for leaf in ("raw", "status", "commands"):
            parent = f"{panel_root}/{leaf}"
            expected = [PurePosixPath(path).name for path in paths
                        if PurePosixPath(path).parent.as_posix() == parent]
            verify_exact_directory_entries(root, parent, expected)
        panel_files = [PurePosixPath(path).name for path in paths
                       if PurePosixPath(path).parent.as_posix() == panel_root]
        verify_exact_directory_entries(
            root, panel_root, panel_files, ("raw", "status", "commands"))
    verify_exact_directory_entries(
        root, "performance",
        ("payload_artifact_manifest.sha256", "payload_summary.json",
         "lookup_artifact_manifest.sha256", "lookup_summary.json",
         "static_size_summary.json"),
        ("payload", "lookup"))


def validate_performance(root: Path) -> None:
    verify_performance_artifact_sets(root)
    for panel, manifest_path, count in (
        ("payload", "performance/payload_artifact_manifest.sha256", EXPECTED_PAYLOAD_ARTIFACTS),
        ("lookup", "performance/lookup_artifact_manifest.sha256", EXPECTED_LOOKUP_ARTIFACTS),
    ):
        verify_hash_manifest(root, manifest_path, count)
        rows = parse_tsv(read_relative(root, manifest_path), ("sha256", "path"), manifest_path)
        if [row["path"] for row in rows] != expected_performance_paths(root, panel):
            fail(f"{manifest_path}: paths/order differ from exact artifact contract")
        finish = parse_json(read_relative(root, f"performance/{panel}/run_finish.json"),
                            f"{panel} run finish")
        if (not isinstance(finish, dict) or not isinstance(finish.get("environment"), dict) or
                finish["environment"].get("gate") != "PASS"):
            fail(f"{panel}: sealed environment report is absent or failed")
        analyzer = root / "freeze/tools/wh2_h15_v5_holdout_analyze.py"
        environment_audit = subprocess.run(
            (sys.executable, str(analyzer), "audit-environment", "--root", str(root),
             "--panel", panel), check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        if environment_audit.returncode != 0:
            fail(f"{panel}: independent environment audit failed: "
                 f"{environment_audit.stderr.strip()}")
    def metric(value: object, label: str, allow_negative: bool = False) -> Decimal:
        try:
            parsed = Decimal(str(value))
        except InvalidOperation:
            fail(f"{label}: invalid decimal")
        if not parsed.is_finite() or (not allow_negative and parsed < 0):
            fail(f"{label}: expected a finite {'signed' if allow_negative else 'nonnegative'} decimal")
        return parsed

    payload = parse_json(read_relative(root, "performance/payload_summary.json"), "payload summary")
    payload_keys = {
            "all_256_solve_ratio", "changed_98_solve_ratio", "class_coverage",
            "coverage_by_phase_width", "coverage_valid", "errors", "excluded_observations",
            "gate", "included_observations", "median_phase_ratio", "muladd_ratio",
            "phase_ratios", "schema", "v4_only_repairs", "v5_only_rank_failures",
            "width_ratios", "xor_ratio"}
    if (not isinstance(payload, dict) or set(payload) != payload_keys or
            payload.get("schema") != "wirehair.wh2.holdout.payload_summary.v1"):
        fail("payload summary schema mismatch")
    integer_fields = ("errors", "excluded_observations", "included_observations",
                      "v4_only_repairs", "v5_only_rank_failures")
    if any(type(payload.get(key)) is not int or payload[key] < 0 for key in integer_fields):
        fail("payload summary count field is malformed")
    class_totals = {"changed-v5": 98 * 16, "immutable-v4": 70 * 16,
                    "structural-control": 28 * 16, "maximin-control": 60 * 16}
    class_coverage = payload.get("class_coverage")
    cells = payload.get("coverage_by_phase_width")
    if (not isinstance(class_coverage, dict) or set(class_coverage) != set(class_totals) or
            any(not isinstance(value, dict) or set(value) != {"included", "excluded"} or
                any(type(value[key]) is not int or value[key] < 0 for key in value)
                for value in class_coverage.values()) or
            any(class_coverage[name]["included"] + class_coverage[name]["excluded"] != total
                for name, total in class_totals.items()) or
            not isinstance(cells, dict) or set(cells) != {
                f"phase{phase}:width{width}" for phase in range(8) for width in (1280, 4096)} or
            any(not isinstance(value, dict) or set(value) != {"included", "excluded"} or
                any(type(value[key]) is not int or value[key] < 0 for key in value) or
                value["included"] < 231 or value["included"] + value["excluded"] != 256
                for value in cells.values())):
        fail("payload coverage ledger is malformed or below its frozen minimum")
    phases = payload.get("phase_ratios")
    widths = payload.get("width_ratios")
    phase_values = ([metric(value, "phase payload ratio") for value in phases]
                    if isinstance(phases, list) and len(phases) == 8 else [])
    phase_median = ((sorted(phase_values)[3] + sorted(phase_values)[4]) / 2
                    if phase_values else Decimal(-1))
    if (payload.get("errors") != 0 or payload.get("v5_only_rank_failures") != 0 or
            payload.get("coverage_valid") is not True or
            payload["included_observations"] + payload["excluded_observations"] != 4096 or
            sum(value["included"] for value in cells.values()) != payload["included_observations"] or
            sum(value["excluded"] for value in cells.values()) != payload["excluded_observations"] or
            sum(value["included"] for value in class_coverage.values()) != payload["included_observations"] or
            sum(value["excluded"] for value in class_coverage.values()) != payload["excluded_observations"] or
            payload["v4_only_repairs"] > payload["excluded_observations"] or
            not isinstance(phases, list) or len(phases) != 8 or
            not isinstance(widths, dict) or set(widths) != {"1280", "4096"} or
            metric(payload["all_256_solve_ratio"], "all payload ratio") > Decimal("1.01") or
            metric(payload["changed_98_solve_ratio"], "changed payload ratio") > Decimal("1.02") or
            any(metric(value, "width payload ratio") > Decimal("1.02") for value in widths.values()) or
            any(value <= 0 for value in phase_values) or
            metric(payload["median_phase_ratio"], "median payload ratio") != phase_median or
            metric(payload["median_phase_ratio"], "median payload ratio") > Decimal("1.01") or
            metric(payload["xor_ratio"], "payload XOR ratio") > Decimal("1.01") or
            metric(payload["muladd_ratio"], "payload muladd ratio") > Decimal("1.0025") or
            payload.get("gate") != "PASS"):
        fail("payload summary violates its exact schema or frozen gates")
    lookup = parse_json(read_relative(root, "performance/lookup_summary.json"), "lookup summary")
    if (not isinstance(lookup, dict) or set(lookup) != {
            "absolute_median_delta_ns", "checksums", "gate", "median_ratio",
            "per_phase_ratios", "schema"} or
            lookup.get("schema") != "wirehair.wh2.holdout.lookup_summary.v1" or
            lookup.get("checksums") != LOOKUP_CHECKSUMS or
            not isinstance(lookup.get("per_phase_ratios"), list) or
            len(lookup["per_phase_ratios"]) != 12):
        fail("lookup summary schema mismatch")
    lookup_ratios = [metric(value, "lookup phase ratio") for value in lookup["per_phase_ratios"]]
    lookup_median = (sorted(lookup_ratios)[5] + sorted(lookup_ratios)[6]) / 2
    if (any(value <= 0 or value > Decimal("1.50") for value in lookup_ratios) or
            metric(lookup["median_ratio"], "lookup median ratio") != lookup_median or
            lookup_median > Decimal("1.25") or
            metric(lookup["absolute_median_delta_ns"], "lookup absolute median delta") > Decimal("2") or
            lookup.get("gate") != "PASS"):
        fail("lookup summary violates its exact schema or frozen gates")
    static = parse_json(read_relative(root, "performance/static_size_summary.json"), "static summary")
    size_keys = ("logical_table_bytes", "lookup_symbol_bytes")
    symbol_sizes = static.get("lookup_symbols") if isinstance(static, dict) else None
    if (not isinstance(static, dict) or set(static) != {
            "candidate_sections", "gate", "logical_table_bytes", "lookup_symbol_bytes",
            "lookup_symbols", "reference_sections", "schema", "whole_binary_section_deltas"} or
            static.get("schema") != "wirehair.wh2.holdout.static_size.v1" or
            any(type(static.get(key)) is not int for key in size_keys) or
            any(not isinstance(static.get(key), dict) or
                set(static[key]) != {"text", "data", "bss"} or
                any(type(value) is not int for value in static[key].values())
                for key in ("candidate_sections", "reference_sections", "whole_binary_section_deltas")) or
            any(static["candidate_sections"][key] - static["reference_sections"][key] !=
                static["whole_binary_section_deltas"][key] for key in ("text", "data", "bss")) or
            static.get("logical_table_bytes") != 756 or static.get("logical_table_bytes") > 2048 or
            not isinstance(symbol_sizes, dict) or set(symbol_sizes) != {
                "benchmark_wrapper", "helper", "keys", "offsets", "salts"} or
            any(type(value) is not int or value <= 0 for value in symbol_sizes.values()) or
            symbol_sizes.get("keys") != 336 or symbol_sizes.get("offsets") != 252 or
            symbol_sizes.get("salts") != 168 or
            static.get("lookup_symbol_bytes") != sum(symbol_sizes.values()) or
            static.get("lookup_symbol_bytes") > 4096 or
            static.get("gate") != "PASS"):
        fail("static-size summary violates its exact schema or frozen gates")
    performance = root / "freeze/tools/wh2_h15_v5_holdout_performance.py"
    performance_audit = subprocess.run(
        (sys.executable, str(performance), "audit-results", "--root", str(root)),
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    if performance_audit.returncode != 0:
        fail(f"independent performance-result audit failed: {performance_audit.stderr.strip()}")
    verify_performance_artifact_sets(root)


def revalidate_nested_performance(
        root: Path, expected_seal: bytes | None = None,
        expected_records: dict[str, str] | None = None) -> None:
    seal_before = read_relative(root, "meta/pre_reveal_seal.sha256")
    records_before = verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    if ((expected_seal is not None and seal_before != expected_seal) or
            (expected_records is not None and records_before != expected_records)):
        fail("pre-reveal seal differs from the caller's validated snapshot")
    validate_performance(root)
    records_after = verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    if (records_after != records_before or
            read_relative(root, "meta/pre_reveal_seal.sha256") != seal_before or
            (expected_seal is not None and seal_before != expected_seal) or
            (expected_records is not None and records_after != expected_records)):
        fail("pre-reveal seal changed while nested performance targets were revalidated")


def seal_pre_reveal(args: argparse.Namespace) -> None:
    root = args.root
    ensure_private_root(root, fresh=False)
    validate_freeze(root)
    validate_performance(root)
    data = seal_data(root, PRE_REVEAL_ROLES)
    write_exclusive(root, "meta/pre_reveal_seal.sha256", data)
    # Revalidate the semantics after materializing the digest list, closing
    # the validate-then-seal mutation window.
    validate_freeze(root)
    validate_performance(root)
    verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    print(f"sealed {len(PRE_REVEAL_ROLES)} pre-reveal roles sha256={sha256_bytes(data)}")


def reveal(args: argparse.Namespace) -> None:
    root = args.root
    ensure_private_root(root, fresh=False)
    remote_ref = canonical_remote_ref(args.remote_ref)
    sealed_tree_path = canonical_tree_path(args.sealed_tree_path)
    freeze_contract = validate_freeze(root)
    seal_records = verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    local_seal = read_relative(root, "meta/pre_reveal_seal.sha256")
    head, remote = verify_git_push(
        args.repo, remote_ref, sealed_tree_path, local_seal,
        str(freeze_contract["source_revision"]))
    verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    if validate_freeze(root) != freeze_contract:
        fail("freeze contract changed during pushed-seal verification")
    if read_relative(root, "meta/pre_reveal_seal.sha256") != local_seal:
        fail("pre-reveal seal changed during pushed-freeze verification")
    if sha256_bytes(read_once(Path(__file__))) != seal_records["prepare"]:
        fail("running revealer bytes differ from the sealed prepare tool")
    revalidate_nested_performance(root, local_seal, seal_records)
    # This is intentionally the process's sole read of commitment bytes.
    source = sys.stdin.buffer.read()
    subset = parse_whole_commitment_once(source, SOURCE_COMMITMENT_SHA256)
    if (read_relative(root, "meta/pre_reveal_seal.sha256") != local_seal or
            verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES) !=
            seal_records):
        fail("pre-reveal seal changed during the one-read commitment input")
    provenance = {
        "emitted_partition": "holdout",
        "emitted_rows": 3,
        "harness_head": head,
        "holdout_reveal_sha256": sha256_bytes(subset),
        "pre_reveal_seal_sha256": sha256_bytes(local_seal),
        "remote_commit": remote,
        "remote_ref": remote_ref,
        "sealed_tree_path": sealed_tree_path,
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.reveal_provenance.v1",
        "source_commitment_reads": 1,
        "source_commitment_sha256": SOURCE_COMMITMENT_SHA256,
        "source_path_accepted": False,
        "sealed_prepare_sha256": seal_records["prepare"],
    }
    write_exclusive(root, "recovery/holdout_reveal.tsv", subset)
    write_exclusive(root, "recovery/reveal_provenance.json", canonical_json(provenance))
    print("emitted canonical three-row holdout reveal snapshot")


def depth_seed(row: dict[str, str]) -> str:
    material = ("wirehair.wh2.normalized-h15-v5.holdout-depth.v1\t" +
                row["index"] + "\t" + row["domain"] + "\t" + row["sha256"] + "\n")
    return "0x" + hashlib.sha256(material.encode("ascii")).hexdigest()[:16]


def read_groups(data: bytes) -> list[tuple[int, tuple[int, ...]]]:
    rows = parse_tsv(data, GROUP_HEADER, "groups")
    result = []
    for line, row in enumerate(rows, 2):
        group = uint(row["group_id"], f"groups:{line}:id", 119)
        ks = tuple(uint(value, f"groups:{line}:K", 64000) for value in row["k_csv"].split(","))
        if group != len(result) or len(ks) != uint(row["k_count"], f"groups:{line}:count"):
            fail(f"groups:{line}: group ordering/count mismatch")
        if ks != tuple(range(2 + group, 64001, 120)):
            fail(f"groups:{line}: not the frozen stride partition")
        result.append((group, ks))
    expected, _ = make_groups()
    if data != expected:
        fail("groups differ from exact reconstructed stride groups")
    return result


def make_recovery_jobs(groups: list[tuple[int, tuple[int, ...]]],
                       reveals: list[dict[str, str]], changed: tuple[int, ...]) -> bytes:
    rows: list[tuple[object, ...]] = []
    def add(phase: str, group: str, seed_index: int, seed: str, schedule: str,
            loss: str, arm: str, profile: str, trials: int, ks: tuple[int, ...]) -> None:
        job_id = len(rows)
        rows.append((job_id, f"recovery{job_id:05d}_{phase}_{group}_s{seed_index}_{schedule}_l{loss[2:]}_{arm}",
                     phase, group, seed_index, seed, schedule, loss, arm, profile,
                     trials, len(ks), ",".join(map(str, ks))))
    for group, ks in groups:
        for seed_index, reveal_row in enumerate(reveals):
            for schedule in SCHEDULES:
                for loss in LOSSES:
                    for arm, profile in ARMS:
                        add("breadth", str(group), seed_index, reveal_row["derived_u64"],
                            schedule, loss, arm, profile, 1, ks)
    for seed_index, reveal_row in enumerate(reveals):
        for schedule in SCHEDULES:
            for loss in LOSSES:
                for arm, profile in ARMS:
                    add("depth", "changed98", seed_index, depth_seed(reveal_row),
                        schedule, loss, arm, profile, 16, changed)
    paired = sum(int(row[10]) * int(row[11]) for row in rows if row[8] == "v4")
    if len(rows) != EXPECTED_RECOVERY_JOBS or paired != EXPECTED_PAIRED_ROWS:
        fail(f"recovery grid mismatch: jobs={len(rows)} paired_rows={paired}")
    return canonical_tsv(RECOVERY_JOB_HEADER, rows)


def validate_reveal_provenance(
        provenance: object, reveal_data: bytes, pre_seal: bytes,
        sealed_prepare_sha256: str) -> None:
    expected_keys = {
        "emitted_partition", "emitted_rows", "harness_head", "holdout_reveal_sha256",
        "pre_reveal_seal_sha256", "remote_commit", "remote_ref", "schema",
        "sealed_prepare_sha256", "sealed_tree_path", "source_commitment_reads",
        "source_commitment_sha256", "source_path_accepted",
    }
    if not isinstance(provenance, dict) or set(provenance) != expected_keys:
        fail("reveal provenance has an unexpected schema")
    remote_ref = provenance.get("remote_ref")
    tree_path = provenance.get("sealed_tree_path")
    harness_head = provenance.get("harness_head")
    remote_commit = provenance.get("remote_commit")
    if (not isinstance(remote_ref, str) or not isinstance(tree_path, str) or
            not isinstance(harness_head, str) or not isinstance(remote_commit, str)):
        fail("reveal provenance Git fields must be strings")
    canonical_remote_ref(remote_ref)
    canonical_tree_path(tree_path)
    if (provenance.get("schema") !=
            "wirehair.wh2.normalized_h15_v5.holdout.reveal_provenance.v1" or
            provenance.get("emitted_partition") != "holdout" or
            type(provenance.get("emitted_rows")) is not int or
            provenance.get("emitted_rows") != 3 or
            type(provenance.get("source_commitment_reads")) is not int or
            provenance.get("source_commitment_reads") != 1 or
            provenance.get("source_path_accepted") is not False or
            provenance.get("source_commitment_sha256") != SOURCE_COMMITMENT_SHA256 or
            provenance.get("holdout_reveal_sha256") != sha256_bytes(reveal_data) or
            provenance.get("pre_reveal_seal_sha256") != sha256_bytes(pre_seal) or
            provenance.get("sealed_prepare_sha256") != sealed_prepare_sha256 or
            not re.fullmatch(r"[0-9a-f]{40}", harness_head) or
            not re.fullmatch(r"[0-9a-f]{40}", remote_commit) or
            harness_head != remote_commit):
        fail("reveal provenance does not prove the one-read/path-free Git contract")


def prepare_recovery(args: argparse.Namespace) -> None:
    root = args.root
    ensure_private_root(root, fresh=False)
    freeze_contract = validate_freeze(root)
    pre_seal = read_relative(root, "meta/pre_reveal_seal.sha256")
    pre_records = verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    revalidate_nested_performance(root, pre_seal, pre_records)
    reveal_data = read_relative(root, "recovery/holdout_reveal.tsv")
    reveals = parse_reveal(reveal_data)
    provenance = parse_json(read_relative(root, "recovery/reveal_provenance.json"), "reveal provenance")
    validate_reveal_provenance(
        provenance, reveal_data, pre_seal, pre_records["prepare"])
    groups_data = read_relative(root, "freeze/groups.tsv")
    groups = read_groups(groups_data)
    table_data = read_relative(root, "freeze/selection_table.tsv")
    table = parse_table(table_data)
    changed = tuple(k for k, row in table.items() if row["changed"] == "1")
    jobs = make_recovery_jobs(groups, reveals, changed)
    if (read_relative(root, "meta/pre_reveal_seal.sha256") != pre_seal or
            verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES) !=
            pre_records):
        fail("pre-reveal seal changed while recovery inputs were prepared")
    revalidate_nested_performance(root, pre_seal, pre_records)
    if validate_freeze(root) != freeze_contract:
        fail("freeze contract changed while recovery inputs were prepared")
    manifest = {
        "artifact_records": EXPECTED_RECOVERY_ARTIFACTS,
        "binary_sha256": freeze_contract["binary_sha256"],
        "holdout_reveal_sha256": sha256_bytes(reveal_data),
        "jobs": EXPECTED_RECOVERY_JOBS,
        "jobs_sha256": sha256_bytes(jobs),
        "paired_rows": EXPECTED_PAIRED_ROWS,
        "pre_reveal_seal_sha256": sha256_bytes(pre_seal),
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.recovery.v1",
        "seed_roles": {"breadth": "raw holdout derived_u64", "depth": "domain-separated SHA-256 derivative"},
    }
    write_exclusive(root, "recovery/jobs.tsv", jobs)
    write_exclusive(root, "recovery/manifest.json", canonical_json(manifest))
    input_seal = seal_data(root, RECOVERY_INPUT_ROLES)
    write_exclusive(root, "recovery/input_seal.sha256", input_seal)
    input_records = verify_seal(root, "recovery/input_seal.sha256", RECOVERY_INPUT_ROLES)
    if (input_records["pre_reveal_seal"] != sha256_bytes(pre_seal) or
            input_records["frozen_binary"] != freeze_contract["binary_sha256"] or
            validate_freeze(root) != freeze_contract or
            verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES) !=
            pre_records):
        fail("recovery input seal does not bind the validated freeze/pre-reveal state")
    print(f"prepared {EXPECTED_RECOVERY_JOBS} jobs / {EXPECTED_PAIRED_ROWS} paired rows")


def run_final_analysis_audit(root: Path) -> None:
    analyzer = root / "freeze/tools/wh2_h15_v5_holdout_analyze.py"
    analysis_audit = subprocess.run(
        (sys.executable, str(analyzer), "audit-analysis", "--root", str(root)),
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    if analysis_audit.returncode != 0:
        fail(f"independent final analysis audit failed: {analysis_audit.stderr.strip()}")


def run_recovery_artifact_audit(root: Path) -> None:
    analyzer = root / "freeze/tools/wh2_h15_v5_holdout_analyze.py"
    artifact_audit = subprocess.run(
        (sys.executable, str(analyzer), "verify-artifacts", "--root", str(root)),
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    if artifact_audit.returncode != 0:
        fail(f"independent recovery-artifact audit failed: {artifact_audit.stderr.strip()}")


def verify_completion_layout(root: Path, sealed: bool) -> None:
    verify_exact_directory_entries(
        root, "freeze/tools",
        tuple(PurePosixPath(path).name for _role, path in PRE_REVEAL_ROLES[9:14]))
    verify_exact_directory_entries(
        root, "freeze",
        ("contract.json", "gates.json", "build_provenance.tsv", "selection_table.tsv",
         "table_dump.tsv", "groups.tsv", "payload_k.tsv", "payload_jobs.tsv",
         "lookup_jobs.tsv", "wirehair_v2_bench", "reference_wirehair_v2_bench"),
        ("tools",))
    verify_exact_directory_entries(
        root, "analysis",
        tuple(PurePosixPath(path).name for _role, path in COMPLETION_ROLES[3:]))
    meta_files = ["pre_reveal_seal.sha256"]
    if sealed:
        meta_files.append("completion_seal.sha256")
    verify_exact_directory_entries(root, "meta", meta_files)
    verify_exact_directory_entries(
        root, ".", (), ("freeze", "performance", "meta", "recovery", "analysis"))


def validate_completion_inputs(root: Path, completion_sealed: bool) -> dict[str, object]:
    verify_completion_layout(root, completion_sealed)
    validate_freeze(root)
    pre_seal = read_relative(root, "meta/pre_reveal_seal.sha256")
    pre_records = verify_seal(root, "meta/pre_reveal_seal.sha256", PRE_REVEAL_ROLES)
    verify_seal(root, "recovery/input_seal.sha256", RECOVERY_INPUT_ROLES)
    verify_hash_manifest(root, "recovery/artifact_manifest.sha256", EXPECTED_RECOVERY_ARTIFACTS)
    run_recovery_artifact_audit(root)
    run_final_analysis_audit(root)
    # Re-expand both nested artifact families after the potentially long
    # recovery analysis, then reassert the direct recovery-input seal.
    revalidate_nested_performance(root, pre_seal, pre_records)
    verify_hash_manifest(root, "recovery/artifact_manifest.sha256", EXPECTED_RECOVERY_ARTIFACTS)
    run_recovery_artifact_audit(root)
    verify_seal(root, "recovery/input_seal.sha256", RECOVERY_INPUT_ROLES)
    verify_completion_layout(root, completion_sealed)
    summary = parse_json(read_relative(root, "analysis/summary.json"), "analysis summary")
    if not isinstance(summary, dict) or summary.get("gate") not in ("PASS", "INCONCLUSIVE"):
        fail("completion refuses a failed or malformed analysis summary")
    return summary


def seal_completion(args: argparse.Namespace) -> None:
    root = args.root
    ensure_private_root(root, fresh=False)
    validate_completion_inputs(root, completion_sealed=False)
    audited_hashes = {role: sha256_bytes(read_relative(root, path))
                      for role, path in COMPLETION_ROLES}
    data = seal_data(root, COMPLETION_ROLES)
    seal_rows = parse_tsv(data, ("role", "sha256", "path"), "completion seal candidate")
    if {row["role"]: row["sha256"] for row in seal_rows} != audited_hashes:
        fail("completion inputs changed after independent analysis audit")
    write_exclusive(root, "meta/completion_seal.sha256", data)
    # Acceptance is post-materialization and expanded: a direct outer-seal
    # check would not notice a raw target changing behind either manifest.
    verify_completion(argparse.Namespace(root=root))
    print(f"sealed {len(COMPLETION_ROLES)} completion roles sha256={sha256_bytes(data)}")


def verify_completion(args: argparse.Namespace) -> None:
    root = args.root
    ensure_private_root(root, fresh=False)
    seal_before = read_relative(root, "meta/completion_seal.sha256")
    records_before = verify_seal(root, "meta/completion_seal.sha256", COMPLETION_ROLES)
    validate_completion_inputs(root, completion_sealed=True)
    records_after = verify_seal(root, "meta/completion_seal.sha256", COMPLETION_ROLES)
    verify_completion_layout(root, sealed=True)
    if (records_after != records_before or
            read_relative(root, "meta/completion_seal.sha256") != seal_before):
        fail("completion seal changed while nested targets were revalidated")
    print("completion verification: PASS (performance + recovery + analysis expanded)")


def expect_failure(action: Callable[[], object], label: str) -> None:
    try:
        action()
    except (AuditError, OSError, UnicodeError, ValueError, csv.Error):
        return
    fail(f"selftest mutation unexpectedly passed: {label}")


def synthetic_table() -> bytes:
    inherited = set(k for k in range(2, 1000) if k not in STRUCTURAL_K)
    inherited = set(sorted(inherited)[:70])
    candidates = [k for k in range(1000, 64001) if k not in STRUCTURAL_K][:98]
    rows = []
    for k in range(2, 64001):
        if k in inherited:
            v4 = v5 = "0x1"
            decision = "immutable-nonzero-v4"
        elif k in candidates:
            v4, v5 = "0x0", "0x27"
            decision = "selected"
        else:
            v4 = v5 = "0x0"
            decision = "retain-zero"
        rows.append((k, v4, v5, int(v4 != v5), decision))
    return canonical_tsv(TABLE_HEADER, rows)


def selftest(_args: argparse.Namespace) -> None:
    for token in (b"NaN", b"Infinity", b"-Infinity"):
        expect_failure(lambda token=token: parse_json(
            b'{"value": ' + token + b"}\n", "nonfinite JSON"),
            f"JSON {token.decode('ascii')}")
    expect_failure(lambda: canonical_json({"value": float("nan")}),
                   "nonfinite JSON serialization")
    seal_path = "campaign/pre_reveal_seal.sha256"
    source_commit = "1" * 40
    head_commit = "2" * 40
    remote_ref = "refs/remotes/origin/holdout"
    binding: dict[str, object] = {
        "source_revision": source_commit, "source_resolved": source_commit,
        "head": head_commit, "parents": (source_commit,), "remote_ref": remote_ref,
        "upstream_ref": remote_ref, "tracking_commit": head_commit,
        "live_commit": head_commit, "live_ref": "refs/heads/holdout",
        "expected_live_ref": "refs/heads/holdout",
        "changes": (("A", seal_path),), "sealed_tree_path": seal_path,
        "sealed_mode": "100644",
    }
    validate_seal_commit_state(**binding)  # type: ignore[arg-type]
    binding_mutations: list[tuple[str, dict[str, object]]] = []
    bad = dict(binding)
    bad.update({"source_revision": "3" * 40, "source_resolved": None})
    binding_mutations.append(("nonexistent arbitrary source revision", bad))
    bad = dict(binding)
    bad.update({"source_revision": "3" * 40, "source_resolved": "3" * 40})
    binding_mutations.append(("wrong existing source revision", bad))
    bad = dict(binding)
    bad["changes"] = (("A", seal_path), ("A", "campaign/extra"))
    binding_mutations.append(("extra seal-commit path", bad))
    bad = dict(binding)
    bad["parents"] = (source_commit, "3" * 40)
    binding_mutations.append(("merge seal commit", bad))
    bad = dict(binding)
    bad.update({"tracking_commit": "3" * 40, "live_commit": "3" * 40})
    binding_mutations.append(("remote ahead/local behind", bad))
    bad = dict(binding)
    bad["sealed_mode"] = "120000"
    binding_mutations.append(("symlink seal blob", bad))
    # This pure gate is called by verify_git_push before reveal's sole stdin
    # read; every provenance/remote mutation must therefore fail pre-input.
    for label, state in binding_mutations:
        expect_failure(
            lambda state=state: validate_seal_commit_state(
                **state),  # type: ignore[arg-type]
            label)
    with tempfile.TemporaryDirectory(prefix="wh2-h15-v5-git-selftest-") as text:
        git_root = Path(text)
        repo = git_root / "repo"
        remote = git_root / "remote.git"
        selftest_git_environment = {
            key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
        def test_git(*arguments: str, bare: bool = False) -> str:
            command = (("git", "--git-dir", str(remote)) if bare else
                       ("git", "-C", str(repo)))
            result = subprocess.run(
                (*command, *arguments), check=False, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, env=selftest_git_environment)
            if result.returncode != 0:
                fail(f"selftest git {' '.join(arguments)} failed: {result.stderr.strip()}")
            return result.stdout.strip()
        subprocess.run(("git", "init", "--bare", str(remote)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=selftest_git_environment)
        subprocess.run(("git", "init", "-b", "source", str(repo)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=selftest_git_environment)
        test_git("config", "user.name", "Holdout Selftest")
        test_git("config", "user.email", "holdout-selftest@example.invalid")
        (repo / "source.txt").write_bytes(b"frozen source\n")
        test_git("add", "source.txt")
        test_git("commit", "-m", "frozen source")
        source_revision = test_git("rev-parse", "HEAD")
        test_git("remote", "add", "origin", str(remote))
        test_git("push", "-u", "origin", "HEAD:refs/heads/holdout")
        seal_relative = "campaign/pre_reveal_seal.sha256"
        seal_bytes = b"role\tsha256\tpath\n"
        (repo / "campaign").mkdir()
        (repo / seal_relative).write_bytes(seal_bytes)
        test_git("add", seal_relative)
        test_git("commit", "-m", "seal only")
        test_git("push", "origin", "HEAD:refs/heads/holdout")
        verified_head, verified_remote = verify_git_push(
            repo, "refs/remotes/origin/holdout", seal_relative,
            seal_bytes, source_revision)
        if verified_head != verified_remote or verified_head != test_git("rev-parse", "HEAD"):
            fail("selftest seal-only Git verification returned inconsistent tips")
        redirect = git_root / "redirect.git"
        subprocess.run(("git", "init", "--bare", str(redirect)), check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=selftest_git_environment)
        injected = {
            "GIT_CONFIG_COUNT": "1", "GIT_CONFIG_KEY_0": "core.useReplaceRefs",
            "GIT_CONFIG_VALUE_0": "true", "GIT_DIR": str(redirect),
            "GIT_REPLACE_REF_BASE": "refs/replace/adversarial",
            "GIT_WORK_TREE": str(git_root / "wrong-worktree"),
        }
        saved_environment = {key: os.environ.get(key) for key in injected}
        os.environ.update(injected)
        try:
            env_head, env_remote = verify_git_push(
                repo, "refs/remotes/origin/holdout", seal_relative,
                seal_bytes, source_revision)
        finally:
            for key, value in saved_environment.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
        if env_head != verified_head or env_remote != verified_remote:
            fail("selftest sanitized Git environment changed verified tips")
        expect_failure(lambda: verify_git_push(
            repo, "refs/remotes/origin/holdout", seal_relative,
            seal_bytes, "3" * 40), "nonexistent source revision in real Git gate")
        valid_tree = test_git("rev-parse", f"{verified_head}^{{tree}}")
        bad_head = test_git(
            "commit-tree", valid_tree, "-p", verified_head, "-m", "bad raw lineage")
        replacement = test_git(
            "commit-tree", valid_tree, "-p", source_revision, "-m", "spoofed lineage")
        test_git("update-ref", "refs/heads/source", bad_head, verified_head)
        test_git("push", "origin", "HEAD:refs/heads/holdout")
        test_git("replace", bad_head, replacement)
        expect_failure(lambda: verify_git_push(
            repo, "refs/remotes/origin/holdout", seal_relative,
            seal_bytes, source_revision), "Git replacement-object parent spoof")
        test_git("replace", "-d", bad_head)
        info = repo / ".git/info"
        info.mkdir(exist_ok=True)
        (info / "grafts").write_text(
            f"{bad_head} {source_revision}\n", encoding="ascii")
        expect_failure(lambda: verify_git_push(
            repo, "refs/remotes/origin/holdout", seal_relative,
            seal_bytes, source_revision), "Git graft parent spoof")
    groups_data, groups = make_groups()
    if len(groups) != 120 or sum(len(ks) for _, ks in groups) != 63999:
        fail("selftest stride geometry mismatch")
    # Use a minimal valid shape for generic geometry; the production path pins
    # and parses the independently sealed table before it reaches this code.
    fake_rows = parse_tsv(synthetic_table(), TABLE_HEADER, "synthetic table")
    table = {int(row["K"]): row for row in fake_rows}
    payload = make_payload_k(table, enforce_hash=False)
    if len(parse_tsv(payload, PAYLOAD_K_HEADER, "synthetic payload")) != 256:
        fail("selftest payload count mismatch")
    if len(parse_tsv(make_payload_jobs(payload), PAYLOAD_JOB_HEADER, "payload jobs")) != 32:
        fail("selftest payload job count mismatch")
    if len(parse_tsv(make_lookup_jobs(), LOOKUP_JOB_HEADER, "lookup jobs")) != 24:
        fail("selftest lookup job count mismatch")
    reveal_rows = []
    selection_rows = []
    for partition, target in (("selection", selection_rows), ("holdout", reveal_rows)):
        for index in range(3):
            digest = hashlib.sha256(f"synthetic-{partition}-{index}".encode()).hexdigest()
            target.append((partition, index, f"wirehair.h15-v5.{partition}.{index}", digest,
                           "0x" + digest[:16], f"{partition} campaign row {index}"))
    reveal_data = canonical_tsv(REVEAL_HEADER, reveal_rows)
    reveals = parse_reveal(reveal_data)
    changed = tuple(int(row["K"]) for row in fake_rows if row["changed"] == "1")
    jobs = make_recovery_jobs(groups, reveals, changed)
    if len(parse_tsv(jobs, RECOVERY_JOB_HEADER, "recovery jobs")) != EXPECTED_RECOVERY_JOBS:
        fail("selftest recovery job count mismatch")
    whole = canonical_tsv(REVEAL_HEADER, selection_rows + reveal_rows)
    if parse_whole_commitment_once(whole, sha256_bytes(whole)) != reveal_data:
        fail("selftest one-read reveal subset mismatch")
    synthetic_pre_seal = b"synthetic pre-reveal seal\n"
    synthetic_prepare_sha = "4" * 64
    synthetic_head = "5" * 40
    good_provenance = {
        "emitted_partition": "holdout", "emitted_rows": 3,
        "harness_head": synthetic_head,
        "holdout_reveal_sha256": sha256_bytes(reveal_data),
        "pre_reveal_seal_sha256": sha256_bytes(synthetic_pre_seal),
        "remote_commit": synthetic_head,
        "remote_ref": "refs/remotes/origin/holdout",
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.reveal_provenance.v1",
        "sealed_prepare_sha256": synthetic_prepare_sha,
        "sealed_tree_path": "campaign/pre_reveal_seal.sha256",
        "source_commitment_reads": 1,
        "source_commitment_sha256": SOURCE_COMMITMENT_SHA256,
        "source_path_accepted": False,
    }
    validate_reveal_provenance(
        good_provenance, reveal_data, synthetic_pre_seal, synthetic_prepare_sha)
    provenance_mutations = (
        ("mismatched pushed commits", "remote_commit", "6" * 40),
        ("boolean commitment-read count", "source_commitment_reads", True),
        ("double-dot remote ref", "remote_ref", "refs/remotes/origin/bad..branch"),
        ("empty remote-ref component", "remote_ref", "refs/remotes/origin//holdout"),
        ("nonportable remote-ref character", "remote_ref", "refs/remotes/origin/bad branch"),
        ("tree path traversal", "sealed_tree_path", "campaign/../seal.sha256"),
        ("noncanonical tree separators", "sealed_tree_path", "campaign//seal.sha256"),
        ("nonportable tree-path character", "sealed_tree_path", "campaign/bad seal.sha256"),
    )
    for label, field, value in provenance_mutations:
        bad_provenance = dict(good_provenance)
        bad_provenance[field] = value
        expect_failure(lambda bad_provenance=bad_provenance: validate_reveal_provenance(
            bad_provenance, reveal_data, synthetic_pre_seal, synthetic_prepare_sha), label)
    mutations = []
    mutations.append(("CRLF reveal", reveal_data.replace(b"\n", b"\r\n")))
    mutations.append(("no final LF", reveal_data[:-1]))
    mutations.append(("reordered", canonical_tsv(REVEAL_HEADER, reversed(reveal_rows))))
    bad = [list(row) for row in reveal_rows]
    bad[2][0] = "selection"
    mutations.append(("selection row", canonical_tsv(REVEAL_HEADER, bad)))
    bad = [list(row) for row in reveal_rows]
    bad[2][4] = "0x0000000000000000"
    mutations.append(("bad derived", canonical_tsv(REVEAL_HEADER, bad)))
    for label, data in mutations:
        expect_failure(lambda data=data: parse_reveal(data), label)
    with tempfile.TemporaryDirectory(prefix="wh2-h15-v5-holdout-selftest-") as text:
        root = Path(text) / "private"
        ensure_private_root(root, fresh=True)
        write_exclusive(root, "a/b", b"stable\n")
        if read_relative(root, "a/b") != b"stable\n":
            fail("selftest relative snapshot mismatch")
        os.symlink("b", root / "a/link")
        expect_failure(lambda: read_relative(root, "a/link"), "leaf symlink")
        os.symlink("a", root / "parent-link")
        expect_failure(lambda: read_relative(root, "parent-link/b"), "parent symlink")
        expect_failure(lambda: read_once(root / "parent-link/b"), "live-source parent symlink")
        expect_failure(lambda: read_relative(root, "../escape"), "path traversal")
        write_exclusive(root, "freeze/payload_jobs.tsv", make_payload_jobs(payload))
        write_exclusive(root, "freeze/lookup_jobs.tsv", make_lookup_jobs())
        for panel in ("payload", "lookup"):
            for path in expected_performance_paths(root, panel):
                write_exclusive(root, path, b"fixture\n")
        for path in (
                "performance/payload_artifact_manifest.sha256",
                "performance/payload_summary.json",
                "performance/lookup_artifact_manifest.sha256",
                "performance/lookup_summary.json",
                "performance/static_size_summary.json"):
            write_exclusive(root, path, b"fixture\n")
        verify_performance_artifact_sets(root)
        write_exclusive(root, "performance/payload/raw/unmanifested.tmp", b"extra\n")
        expect_failure(lambda: verify_performance_artifact_sets(root),
                       "post-performance-seal unmanifested extra file")
        (root / "performance/payload/raw/unmanifested.tmp").unlink()
        verify_performance_artifact_sets(root)
        write_exclusive(root, "sealed/item", b"sealed\n")
        roles = (("item", "sealed/item"),)
        write_exclusive(root, "sealed/seal.tsv", seal_data(root, roles))
        verify_seal(root, "sealed/seal.tsv", roles)
        (root / "sealed/item").write_bytes(b"mutated\n")
        expect_failure(lambda: verify_seal(root, "sealed/seal.tsv", roles),
                       "post-seal artifact mutation")
        # A completion seal intentionally binds the nested artifact manifest,
        # not every large raw file directly.  Prove that the outer seal alone
        # remains valid after a target mutation/deletion and that the expanded
        # verifier closes that otherwise-invisible gap.
        target_relative = "nested/raw/result.tsv"
        manifest_relative = "nested/recovery_artifact_manifest.sha256"
        summary_relative = "nested/analysis_summary.json"
        completion_relative = "nested/completion_seal.sha256"
        original = b"raw result\n"
        write_exclusive(root, target_relative, original)
        write_exclusive(root, manifest_relative, canonical_tsv(
            ("sha256", "path"), ((sha256_bytes(original), target_relative),)))
        write_exclusive(root, summary_relative, b'{"gate":"PASS"}\n')
        completion_roles = (
            ("recovery_artifact_manifest", manifest_relative),
            ("summary", summary_relative),
        )
        write_exclusive(root, completion_relative, seal_data(root, completion_roles))
        verify_seal(root, completion_relative, completion_roles)
        verify_hash_manifest(root, manifest_relative, 1)
        verify_exact_directory_entries(root, "nested/raw", ("result.tsv",))
        write_exclusive(root, "nested/raw/unmanifested.tmp", b"extra\n")
        verify_seal(root, completion_relative, completion_roles)
        verify_hash_manifest(root, manifest_relative, 1)
        expect_failure(lambda: verify_exact_directory_entries(
            root, "nested/raw", ("result.tsv",)),
            "post-completion unmanifested extra file")
        (root / "nested/raw/unmanifested.tmp").unlink()
        verify_exact_directory_entries(root, "nested/raw", ("result.tsv",))
        (root / target_relative).write_bytes(b"mutated raw result\n")
        verify_seal(root, completion_relative, completion_roles)
        expect_failure(lambda: verify_hash_manifest(root, manifest_relative, 1),
                       "nested target mutation after completion seal")
        (root / target_relative).write_bytes(original)
        verify_hash_manifest(root, manifest_relative, 1)
        (root / target_relative).unlink()
        verify_seal(root, completion_relative, completion_roles)
        expect_failure(lambda: verify_hash_manifest(root, manifest_relative, 1),
                       "nested target deletion after completion seal")
    with tempfile.TemporaryDirectory(prefix="wh2-completion-layout-selftest-") as text:
        root = Path(text) / "private"
        ensure_private_root(root, fresh=True)
        for relative in ("performance", "recovery"):
            ensure_dir(root, relative)
        for relative in (
                "freeze/contract.json", "freeze/gates.json", "freeze/build_provenance.tsv",
                "freeze/selection_table.tsv", "freeze/table_dump.tsv", "freeze/groups.tsv",
                "freeze/payload_k.tsv", "freeze/payload_jobs.tsv", "freeze/lookup_jobs.tsv",
                "freeze/wirehair_v2_bench", "freeze/reference_wirehair_v2_bench"):
            write_exclusive(root, relative, b"fixture\n")
        for _role, relative in PRE_REVEAL_ROLES[9:14]:
            write_exclusive(root, relative, b"fixture\n")
        for _role, relative in COMPLETION_ROLES[3:]:
            write_exclusive(root, relative, b"fixture\n")
        write_exclusive(root, "meta/pre_reveal_seal.sha256", b"fixture\n")
        write_exclusive(root, "meta/completion_seal.sha256", b"fixture\n")
        verify_completion_layout(root, sealed=True)
        write_exclusive(root, "meta/unmanifested.tmp", b"extra\n")
        expect_failure(lambda: verify_completion_layout(root, sealed=True),
                       "post-completion extra metadata file")
        (root / "meta/unmanifested.tmp").unlink()
        write_exclusive(root, "unmanifested.tmp", b"extra\n")
        expect_failure(lambda: verify_completion_layout(root, sealed=True),
                       "post-completion extra campaign-root file")
    print("selftest: ok (geometry, raw Git provenance, exact completion sets, nested seals)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    freeze = commands.add_parser("prepare-freeze")
    freeze.add_argument("--output", type=Path, required=True)
    freeze.add_argument("--binary", type=Path, required=True)
    freeze.add_argument("--reference-binary", type=Path, required=True)
    freeze.add_argument("--selection-table", type=Path, required=True)
    freeze.add_argument("--table-dump", type=Path, required=True)
    freeze.add_argument("--source-commit", required=True)
    freeze.add_argument("--prepare", type=Path, required=True)
    freeze.add_argument("--job-runner", type=Path, required=True)
    freeze.add_argument("--launcher", type=Path, required=True)
    freeze.add_argument("--analyzer", type=Path, required=True)
    freeze.add_argument("--performance", type=Path, required=True)
    freeze.set_defaults(func=prepare_freeze)
    preseal = commands.add_parser("seal-pre-reveal")
    preseal.add_argument("--root", type=Path, required=True)
    preseal.set_defaults(func=seal_pre_reveal)
    reveal_parser = commands.add_parser("reveal")
    reveal_parser.add_argument("--root", type=Path, required=True)
    reveal_parser.add_argument("--repo", type=Path, required=True)
    reveal_parser.add_argument("--remote-ref", required=True)
    reveal_parser.add_argument("--sealed-tree-path", required=True,
                               help="repo-relative path containing the exact committed pre-reveal seal")
    reveal_parser.set_defaults(func=reveal)
    recovery = commands.add_parser("prepare-recovery")
    recovery.add_argument("--root", type=Path, required=True)
    recovery.set_defaults(func=prepare_recovery)
    complete = commands.add_parser("seal-completion")
    complete.add_argument("--root", type=Path, required=True)
    complete.set_defaults(func=seal_completion)
    verify_complete = commands.add_parser("verify-completion")
    verify_complete.add_argument("--root", type=Path, required=True)
    verify_complete.set_defaults(func=verify_completion)
    test = commands.add_parser("selftest")
    test.set_defaults(func=selftest)
    args = parser.parse_args()
    try:
        args.func(args)
    except (AuditError, OSError, UnicodeError, ValueError, InvalidOperation, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
