#!/usr/bin/env python3
"""Prepare the sealed WH2 H15-v5 Stage G/Stage T selection campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

GROUPS_SHA256 = "940db61d44fc03462f583c7e2e48bea5c93a830918c8f8e5bf3ba502e840cff4"
DISCOVERY_SHA256 = "ecc722c5e44aee2ed468948cc4c7ee903ec4943b49ea4b67778b37a598e0b125"
IMMUTABLE_SHA256 = (
    ("discovery_failures.csv", "c9fe1730f632bc406c229391f137a30db5f933a5d21e405d822658b74e2d5dd7"),
    ("exact_salt_scores.csv", "ff238d4a21a340c866408ada53d22da474a880f054c1d50ecbba07007bf128e4"),
    ("requested_k_salt.tsv", "93fe6d38b090d9054e0c293276d3caef7e0fd1dc2306105fdfd9351755b2c130"),
)
SOURCE_COMMITMENT_SHA256 = "f76eedbf4813597fa84af4b21a1236c957d11ff4a745585ea0442d9676645728"
SOURCE_COMMITMENT_BINDING = ("opaque whole-file procedural provenance only; "
                             "no membership claim; no source commitment path was accepted or read")
ARMS = (("baseline", 0x00), ("c27", 0x27), ("c79", 0x79),
        ("c6f", 0x6F), ("ca8", 0xA8))
G_SCHEDULES = ("burst", "adversarial", "repair-only")
T_SCHEDULES = ("iid", "burst", "permutation", "systematic-first",
               "repair-only", "adversarial")
T_LOSSES = ("0.35", "0.50", "0.65")
REVEAL_HEADER = ("partition", "index", "domain", "sha256", "derived_u64", "usage")
DISCOVERY_HEADER = ("K", "certified_v4_salt", "witness_count", "witness_ids")
JOB_HEADER = (
    "job_id", "stem", "stage", "schedule", "seed_index", "seed", "loss",
    "arm", "salt", "trials", "k_count", "k_csv",
)
EXPECTED_G_GROUPS = 120
EXPECTED_G_K = 63929
EXPECTED_DISCOVERY_K = 466
EXPECTED_G_LOGICAL_TRIALS = 958935
EXPECTED_T_LOGICAL_TRIALS = 167760
EXPECTED_JOBS = 1980
EXPECTED_LOGICAL_TRIALS = 1126695


class AuditError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise AuditError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    actual = file_sha256(path)
    if actual != expected:
        fail(f"{path}: SHA-256 mismatch: expected {expected}, found {actual}")


def uint(text: str, context: str, maximum: int | None = None) -> int:
    if not re.fullmatch(r"(?:0|[1-9][0-9]*)", text):
        fail(f"{context}: noncanonical unsigned integer {text!r}")
    value = int(text)
    if maximum is not None and value > maximum:
        fail(f"{context}: value exceeds {maximum}")
    return value


def salt(text: str, context: str) -> int:
    if not re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]*)", text):
        fail(f"{context}: noncanonical lowercase hexadecimal salt {text!r}")
    value = int(text, 16)
    if value > 255:
        fail(f"{context}: salt exceeds uint8")
    return value


def canonical_text(text: str, context: str) -> None:
    if not text or text != text.strip() or not re.fullmatch(r"[\x20-\x7e]+", text):
        fail(f"{context}: expected nonempty canonical printable-ASCII text")


def read_groups(path: Path, expected_sha: str = GROUPS_SHA256) -> list[tuple[int, tuple[int, ...]]]:
    require_hash(path, expected_sha)
    groups: list[tuple[int, int, tuple[int, ...]]] = []
    covered: set[int] = set()
    with path.open(newline="") as stream:
        for line_number, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if len(fields) != 3:
                fail(f"{path}:{line_number}: expected three TSV fields")
            group_id = uint(fields[0], f"{path}:{line_number}:group_id", 185)
            if group_id != line_number - 1:
                fail(f"{path}:{line_number}: group IDs must be ordered 0..185")
            group_salt = salt(fields[1], f"{path}:{line_number}:salt")
            parts = fields[2].split(",")
            ks = tuple(uint(part, f"{path}:{line_number}:K", 64000) for part in parts)
            if not ks or ks != tuple(sorted(set(ks))) or ks[0] < 2:
                fail(f"{path}:{line_number}: K list must be nonempty, unique, and sorted")
            if covered.intersection(ks):
                fail(f"{path}:{line_number}: duplicate K across groups")
            covered.update(ks)
            groups.append((group_id, group_salt, ks))
    if len(groups) != 186 or covered != set(range(2, 64001)):
        fail(f"{path}: expected 186 groups partitioning exact K range 2..64000")
    selected = [(group_id, ks) for group_id, group_salt, ks in groups if group_salt == 0]
    if len(selected) != EXPECTED_G_GROUPS or sum(len(ks) for _, ks in selected) != EXPECTED_G_K:
        fail(f"{path}: expected exactly 120 salt-0 groups covering 63929 K")
    return selected


def read_discovery(path: Path, expected_sha: str = DISCOVERY_SHA256) -> tuple[int, ...]:
    require_hash(path, expected_sha)
    ks: list[int] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != DISCOVERY_HEADER:
            fail(f"{path}: discovery header mismatch")
        for line_number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                fail(f"{path}:{line_number}: malformed CSV row")
            k = uint(row["K"], f"{path}:{line_number}:K", 64000)
            if k < 2:
                fail(f"{path}:{line_number}: K is below 2")
            if salt(row["certified_v4_salt"], f"{path}:{line_number}:salt") != 0:
                fail(f"{path}:{line_number}: discovery K is not certified salt 0")
            if uint(row["witness_count"], f"{path}:{line_number}:witness_count") == 0:
                fail(f"{path}:{line_number}: empty witness count")
            canonical_text(row["witness_ids"], f"{path}:{line_number}:witness_ids")
            ks.append(k)
    if len(ks) != EXPECTED_DISCOVERY_K or ks != sorted(set(ks)):
        fail(f"{path}: expected exactly 466 unique, sorted salt-0 K")
    return tuple(ks)


def read_reveal(path: Path) -> list[dict[str, str]]:
    rows: dict[int, dict[str, str]] = {}
    seen_sha: set[str] = set()
    seen_derived: set[str] = set()
    seen_domains: set[str] = set()
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        fail(f"{path}: reveal must be valid UTF-8")
    if b"\r" in data or not data.endswith(b"\n") or b"\n\n" in data:
        fail(f"{path}: reveal requires LF-only lines, final newline, and no blank lines")
    with io.StringIO(text, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REVEAL_HEADER:
            fail(f"{path}: reveal header mismatch")
        for line_number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                fail(f"{path}:{line_number}: malformed reveal row")
            if row["partition"] != "selection":
                fail(f"{path}:{line_number}: only the selection partition is permitted")
            index = uint(row["index"], f"{path}:{line_number}:index", 2)
            if index in rows:
                fail(f"{path}:{line_number}: duplicate selection index {index}")
            canonical_text(row["domain"], f"{path}:{line_number}:domain")
            canonical_text(row["usage"], f"{path}:{line_number}:usage")
            domain_lower = row["domain"].lower()
            if "selection" not in domain_lower or "holdout" in domain_lower:
                fail(f"{path}:{line_number}: domain must identify selection, never holdout")
            if domain_lower in seen_domains:
                fail(f"{path}:{line_number}: duplicate selection domain")
            if "holdout" in row["usage"].lower():
                fail(f"{path}:{line_number}: usage must not identify holdout")
            digest = row["sha256"]
            if not re.fullmatch(r"[0-9a-f]{64}", digest):
                fail(f"{path}:{line_number}: SHA-256 must be 64 lowercase hex digits")
            if digest in seen_sha:
                fail(f"{path}:{line_number}: duplicate revealed SHA-256")
            expected_derived = "0x" + digest[:16]
            if row["derived_u64"] != expected_derived:
                fail(f"{path}:{line_number}: derived_u64 must equal {expected_derived}")
            if expected_derived in seen_derived:
                fail(f"{path}:{line_number}: duplicate derived selection seed")
            seen_sha.add(digest)
            seen_derived.add(expected_derived)
            seen_domains.add(domain_lower)
            rows[index] = dict(row)
    if set(rows) != {0, 1, 2}:
        fail(f"{path}: expected exact selection indices 0..2")
    ordered = [rows[index] for index in range(3)]
    canonical = "\t".join(REVEAL_HEADER) + "\n"
    canonical += "".join(
        "\t".join(row[field] for field in REVEAL_HEADER) + "\n"
        for row in ordered
    )
    if data != canonical.encode("ascii"):
        fail(f"{path}: reveal must use exact canonical unquoted TSV serialization")
    return ordered


def make_jobs(
    groups: list[tuple[int, tuple[int, ...]]], discovery: tuple[int, ...],
    reveals: list[dict[str, str]],
) -> list[dict[str, object]]:
    jobs: list[dict[str, object]] = []

    def add(stage: str, stem_tail: str, schedule: str, seed_index: int,
            loss: str, arm: str, arm_salt: int, trials: int, ks: tuple[int, ...]) -> None:
        job_id = len(jobs)
        jobs.append({
            "job_id": job_id,
            "stem": f"job{job_id:04d}_{stage}_{stem_tail}",
            "stage": stage,
            "schedule": schedule,
            "seed_index": seed_index,
            "seed": reveals[seed_index]["derived_u64"],
            "loss": loss,
            "arm": arm,
            "salt": arm_salt,
            "trials": trials,
            "k_count": len(ks),
            "k_csv": ",".join(map(str, ks)),
        })
    for group_id, ks in groups:
        for schedule in G_SCHEDULES:
            for arm, arm_salt in ARMS:
                add("G", f"group{group_id:03d}_{schedule}_seed0_{arm}_l050",
                    schedule, 0, "0.50", arm, arm_salt, 1, ks)
    for seed_index in (1, 2):
        for schedule in T_SCHEDULES:
            for loss in T_LOSSES:
                for arm, arm_salt in ARMS:
                    add("T", f"{schedule}_seed{seed_index}_{arm}_l{loss.replace('.', '')}",
                        schedule, seed_index, loss, arm, arm_salt, 2, discovery)
    g_jobs = [row for row in jobs if row["stage"] == "G"]
    t_jobs = [row for row in jobs if row["stage"] == "T"]
    g_trials = sum(int(row["trials"]) * int(row["k_count"]) for row in g_jobs)
    t_trials = sum(int(row["trials"]) * int(row["k_count"]) for row in t_jobs)
    logical_trials = g_trials + t_trials
    if (len(g_jobs), len(t_jobs), g_trials, t_trials) != (
        1800, 180, EXPECTED_G_LOGICAL_TRIALS, EXPECTED_T_LOGICAL_TRIALS,
    ) or len(jobs) != EXPECTED_JOBS or logical_trials != EXPECTED_LOGICAL_TRIALS:
        fail(f"constructed grid mismatch: jobs={len(jobs)} logical_trials={logical_trials}")
    return jobs


def write_tsv(path: Path, header: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_exclusive(source: Path, destination: Path) -> None:
    with source.open("rb") as reader, destination.open("xb") as writer:
        for chunk in iter(lambda: reader.read(1024 * 1024), b""):
            writer.write(chunk)


def prepare_paths(
    groups_path: Path, discovery_path: Path, reveal_path: Path, output: Path,
    immutable: tuple[tuple[str, Path, str], ...],
    groups_sha: str = GROUPS_SHA256, discovery_sha: str = DISCOVERY_SHA256,
) -> dict[str, object]:
    groups = read_groups(groups_path, groups_sha)
    discovery = read_discovery(discovery_path, discovery_sha)
    reveals = read_reveal(reveal_path)
    group_k = {k for _, ks in groups for k in ks}
    if not set(discovery).issubset(group_k):
        fail("discovery K are not a subset of the salt-0 source groups")
    if tuple(name for name, _, _ in immutable) != tuple(name for name, _ in IMMUTABLE_SHA256):
        fail("immutable metadata names/order mismatch")
    for _, path, expected in immutable:
        require_hash(path, expected)
    meta_inputs = (("groups.tsv", groups_path, groups_sha),
                   ("discovery_k.csv", discovery_path, discovery_sha)) + immutable
    jobs = make_jobs(groups, discovery, reveals)
    if output.exists():
        if not output.is_dir() or any(output.iterdir()):
            fail(f"output must be a fresh or empty directory: {output}")
    else:
        output.mkdir()
    meta = output / "meta"
    meta.mkdir()
    reveal_rows: list[dict[str, object]] = [dict(row) for row in reveals]
    write_tsv(meta / "selection_seeds.tsv", REVEAL_HEADER, reveal_rows)
    write_tsv(meta / "jobs.tsv", JOB_HEADER, jobs)
    selection_seeds_sha = file_sha256(meta / "selection_seeds.tsv")
    for name, source, expected in meta_inputs:
        copy_exclusive(source, meta / name)
        require_hash(meta / name, expected)
    manifest: dict[str, object] = {
        "schema": "wirehair.wh2.normalized_h15_v5.selection.v1",
        "source_commitment_sha256": SOURCE_COMMITMENT_SHA256,
        "source_commitment_binding": SOURCE_COMMITMENT_BINDING,
        # Bind the already validated snapshots copied above.  Do not reopen
        # live source paths after the campaign inputs have been materialized.
        "groups_sha256": groups_sha,
        "discovery_k_sha256": discovery_sha,
        # read_reveal accepts only this exact canonical serialization, so the
        # copied snapshot digest is also the digest of the bytes that were
        # parsed.  Never re-hash a potentially changed live source path here.
        "selection_reveal_sha256": selection_seeds_sha,
        "selection_seeds_sha256": selection_seeds_sha,
        "jobs_sha256": file_sha256(meta / "jobs.tsv"),
        "immutable_meta_sha256": {name: expected for name, _, expected in meta_inputs},
        "selection_seed_roles": {"0": "StageG", "1": "StageT", "2": "StageT"},
        "arms": [{"arm": arm, "salt": value, "salt_hex": f"0x{value:x}"}
                 for arm, value in ARMS],
        "stage_g": {"groups": len(groups), "k": EXPECTED_G_K, "jobs": 1800,
                    "logical_trials": EXPECTED_G_LOGICAL_TRIALS},
        "stage_t": {"seeds": 2, "k": len(discovery), "jobs": 180,
                    "logical_trials": EXPECTED_T_LOGICAL_TRIALS},
        "jobs": len(jobs),
        "logical_trials": EXPECTED_LOGICAL_TRIALS,
    }
    with (meta / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


def prepare(args: argparse.Namespace) -> None:
    paths = {name: getattr(args, Path(name).stem.replace("-", "_"))
             for name, _ in IMMUTABLE_SHA256}
    immutable = tuple((name, paths[name], digest) for name, digest in IMMUTABLE_SHA256)
    manifest = prepare_paths(args.groups, args.discovery_k, args.reveal, args.output, immutable)
    print(f"prepared {manifest['jobs']} jobs / {manifest['logical_trials']} logical trials")


def expect_failure(action: Callable[[], object], label: str) -> None:
    try:
        action()
    except AuditError:
        return
    fail(f"selftest mutation unexpectedly passed: {label}")


def split_evenly(values: list[int], count: int) -> list[list[int]]:
    quotient, remainder = divmod(len(values), count)
    result: list[list[int]] = []
    offset = 0
    for index in range(count):
        size = quotient + (index < remainder)
        result.append(values[offset:offset + size])
        offset += size
    return result


def selftest(_args: argparse.Namespace) -> None:
    with tempfile.TemporaryDirectory(prefix="wh2-h15-v5-selection-selftest-") as temp_text:
        temp = Path(temp_text)
        groups_path = temp / "groups.tsv"
        zero_chunks = split_evenly(list(range(2, 63931)), 120)
        other_chunks = split_evenly(list(range(63931, 64001)), 66)
        with groups_path.open("x", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            for group_id, ks in enumerate(zero_chunks + other_chunks):
                writer.writerow((group_id, "0x0" if group_id < 120 else "0x1",
                                 ",".join(map(str, ks))))
        discovery_path = temp / "discovery_k.csv"
        with discovery_path.open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=DISCOVERY_HEADER, lineterminator="\n")
            writer.writeheader()
            for witness, k in enumerate(range(2, 468)):
                writer.writerow({"K": k, "certified_v4_salt": "0x0", "witness_count": 1,
                                 "witness_ids": f"w{witness:04d}"})
        base_rows: list[dict[str, object]] = []
        for index in range(3):
            digest = hashlib.sha256(f"synthetic-selection-{index}".encode()).hexdigest()
            base_rows.append({
                "partition": "selection", "index": index,
                "domain": f"wirehair.wh2.h15-v5.selection.{index}", "sha256": digest,
                "derived_u64": "0x" + digest[:16],
                "usage": "StageG" if index == 0 else "StageT",
            })

        def reveal_file(name: str, rows: list[dict[str, object]]) -> Path:
            path = temp / name
            write_tsv(path, REVEAL_HEADER, rows)
            return path

        reveal = reveal_file("reveal.tsv", base_rows)
        immutable_list = []
        for name, _ in IMMUTABLE_SHA256:
            path = temp / name
            path.write_text(f"synthetic {name}\n")
            immutable_list.append((name, path, file_sha256(path)))
        immutable = tuple(immutable_list)
        output = temp / "output"
        manifest = prepare_paths(groups_path, discovery_path, reveal, output, immutable,
                                 file_sha256(groups_path), file_sha256(discovery_path))
        if manifest["jobs"] != EXPECTED_JOBS or manifest["logical_trials"] != EXPECTED_LOGICAL_TRIALS:
            fail("selftest generated counts mismatch")
        expected_contract = {
            "source_commitment_sha256": SOURCE_COMMITMENT_SHA256,
            "source_commitment_binding": SOURCE_COMMITMENT_BINDING,
            "selection_seed_roles": {"0": "StageG", "1": "StageT", "2": "StageT"},
            "arms": [{"arm": arm, "salt": value, "salt_hex": f"0x{value:x}"} for arm, value in ARMS],
            "stage_g": {"groups": 120, "k": 63929, "jobs": 1800, "logical_trials": 958935},
            "stage_t": {"seeds": 2, "k": 466, "jobs": 180, "logical_trials": 167760},
        }
        if any(manifest.get(key) != value for key, value in expected_contract.items()):
            fail("selftest manifest contract mismatch")
        for name, digest in manifest["immutable_meta_sha256"].items():
            require_hash(output / "meta" / name, digest)
        with (output / "meta/jobs.tsv").open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if tuple(reader.fieldnames or ()) != JOB_HEADER:
                fail("selftest jobs.tsv header mismatch")
            jobs = list(reader)
        stage_counts = {stage: sum(row["stage"] == stage for row in jobs) for stage in ("G", "T")}
        if (len(jobs) != EXPECTED_JOBS or stage_counts != {"G": 1800, "T": 180} or
                sum(int(row["trials"]) * int(row["k_count"]) for row in jobs)
                != EXPECTED_LOGICAL_TRIALS):
            fail("selftest jobs.tsv ledger mismatch")
        expect_failure(lambda: read_groups(groups_path, "0" * 64), "bad hash")
        mutations: list[tuple[str, Path]] = []
        reveal_data = reveal.read_bytes()
        for label, payload in (("CRLF", reveal_data.replace(b"\n", b"\r\n")),
                               ("no final LF", reveal_data[:-1]),
                               ("blank line", reveal_data.replace(b"\n", b"\n\n", 1)),
                               ("invalid UTF-8", b"\xff" + reveal_data),
                               ("quoted field", reveal_data.replace(
                                   str(base_rows[0]["domain"]).encode(),
                                   b'"' + str(base_rows[0]["domain"]).encode() + b'"', 1)),
                               ("reordered rows", b"\n".join(
                                   reveal_data.splitlines()[:1] +
                                   list(reversed(reveal_data.splitlines()[1:]))) + b"\n")):
            path = temp / f"bad-{label.replace(' ', '-')}.tsv"
            path.write_bytes(payload)
            mutations.append((label, path))
        for label, field, value in (
            ("bad index", "index", "03"),
            ("holdout partition", "partition", "holdout"),
            ("holdout domain", "domain", "wirehair.selection.holdout.2"),
            ("holdout usage", "usage", "selection-holdout-provenance"),
            ("noncanonical usage", "usage", " selection provenance"),
            ("bad derived", "derived_u64", "0x0000000000000000"),
            ("noncanonical derived", "derived_u64", "0X" + str(base_rows[2]["sha256"])[:16]),
            ("duplicate index", "index", "1"),
            ("duplicate domain", "domain", str(base_rows[1]["domain"]).upper()),
        ):
            rows = [dict(row) for row in base_rows]
            rows[2][field] = value
            mutations.append((label, reveal_file(f"bad-{label.replace(' ', '-')}.tsv", rows)))
        duplicate_sha = [dict(row) for row in base_rows]
        duplicate_sha[2]["sha256"] = duplicate_sha[1]["sha256"]
        duplicate_sha[2]["derived_u64"] = duplicate_sha[1]["derived_u64"]
        mutations.append(("duplicate SHA", reveal_file("bad-duplicate-sha.tsv", duplicate_sha)))
        duplicate_derived = [dict(row) for row in base_rows]
        prefix = str(duplicate_derived[1]["sha256"])[:16]
        duplicate_derived[2]["sha256"] = prefix + str(duplicate_derived[2]["sha256"])[16:]
        duplicate_derived[2]["derived_u64"] = "0x" + prefix
        mutations.append(("duplicate derived", reveal_file("bad-duplicate-derived.tsv", duplicate_derived)))
        bad_header = temp / "bad-header.tsv"
        write_tsv(bad_header, tuple(reversed(REVEAL_HEADER)), base_rows)
        mutations.append(("bad header", bad_header))
        for label, path in mutations:
            expect_failure(lambda path=path: read_reveal(path), label)
        provenance_rows = [dict(row) for row in base_rows]
        for index, row in enumerate(provenance_rows):
            row["usage"] = f"arbitrary provenance {index}"
        provenance = read_reveal(reveal_file("provenance-usage.tsv", provenance_rows))
        provenance_jobs = make_jobs(
            [(group_id, tuple(ks)) for group_id, ks in enumerate(zero_chunks)],
            tuple(range(2, 468)), provenance,
        )
        if any(row["seed_index"] != 0 for row in provenance_jobs if row["stage"] == "G"):
            fail("usage provenance incorrectly changed the Stage G role")
        if any(row["seed_index"] not in (1, 2) for row in provenance_jobs if row["stage"] == "T"):
            fail("usage provenance incorrectly changed a Stage T role")
    print(f"selftest: ok (1980 jobs, 1126695 logical trials; {1 + len(mutations)} mutations rejected)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--groups", type=Path, required=True)
    prepare_parser.add_argument("--discovery-k", type=Path, required=True)
    prepare_parser.add_argument("--reveal", type=Path, required=True)
    prepare_parser.add_argument("--output", type=Path, required=True)
    prepare_parser.add_argument("--discovery-failures", type=Path, required=True)
    prepare_parser.add_argument("--exact-salt-scores", type=Path, required=True)
    prepare_parser.add_argument("--requested-k-salt", type=Path, required=True)
    prepare_parser.set_defaults(func=prepare)
    selftest_parser = commands.add_parser("selftest")
    selftest_parser.set_defaults(func=selftest)
    args = parser.parse_args()
    try:
        args.func(args)
    except (AuditError, OSError, UnicodeError, csv.Error, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
