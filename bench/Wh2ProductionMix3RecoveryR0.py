#!/usr/bin/env python3
"""Bounded public-MIX3 finite-map feasibility pilot, never a promotion gate."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import struct
import subprocess
import sys
import time

PROTOCOL = "wirehair.wh2.production-mix3-k3k6-r0"
KS = (3, 6)
OVERHEADS = (0, 1, 2, 4)
SCHEDULES = ("burst", "adversarial", "repair-only")
TRAINING_ROOTS = (
    "0x7ccd510f122fc160", "0xb889883a79549774", "0xb5666de0987896af",
    "0x30a0ac4ab2e861cc", "0x20e6b10b1cc3838e", "0x894e0a8fedcd6cb5")
HOLDOUT_ROOTS = (
    "0x688e4a7ca826b448", "0xa2c6fb6d887f8efe", "0xe07d30ba3a9d921f",
    "0x958b941967c7fbba", "0x6f00ea1b271de4eb", "0x754a9eb98cd69323",
    "0x3f3f217a15779be6", "0x8307ffc9acf675e0", "0xb7af2e5f6ccb44fa")
ROOT = Path(__file__).resolve().parent.parent
WORKER_SOURCE = ROOT / "bench/Wh2ProductionMix3RecoveryR0.cpp"
CELL_KEYS = {"type", "phase", "arm", "K", "attempt", "root_index", "root", "schedule",
             "profile_hex", "profile_sha256", "source_sha256", "attempted_candidates",
             "trace_sha256", "ids", "packets_hex", "outcomes", "recovered_sha256"}
MASK64 = (1 << 64) - 1


class ValidationError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")


def digest(value):
    return hashlib.sha256(value).hexdigest()


def file_digest(path):
    return digest(Path(path).read_bytes())


def integer(value, lower, upper, label):
    require(type(value) is int and lower <= value <= upper, label)
    return value


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def parse_json(value):
    try:
        return json.loads(value, object_pairs_hook=no_duplicates,
                          parse_constant=lambda value: (_ for _ in ()).throw(ValidationError(value)))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValidationError(str(error)) from error


def source_bytes(K):
    return bytes((73 * i + 19 * K + 11) & 255 for i in range(K * 2))


def profile_bytes(K, attempt):
    # Exact public 32-byte WHV2 wire format, not the native C struct layout.
    return b"WHV2" + struct.pack("<HHQQIB3s", 1, 32, 0x4b295bbb47f4f9c9,
                                K * 2, 2, attempt, bytes(3))


def trace(K, root_hex, schedule):
    require(K in KS and schedule in SCHEDULES, "trace shape")
    state = (int(root_hex, 16) ^ (K * 0x9e3779b97f4a7c15) ^ (2 * 0xbf58476d1ce4e5b9)) & MASK64
    state ^= 0x10fade

    def unit():
        nonlocal state
        state = (state + 0x9e3779b97f4a7c15) & MASK64
        value = state
        value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        value ^= value >> 31
        return (value >> 11) / float(1 << 53)

    ids, attempted, burst_left = [], 0, 0
    while len(ids) < K + 4:
        require(attempted < 256 * (K + 4) + 65536, "trace candidate limit")
        candidate = attempted
        attempted += 1
        if schedule == "burst":
            if burst_left:
                burst_left -= 1
                drop = True
            elif unit() < 0.5 / (8.0 - 7.0 * 0.5):
                burst_left, drop = 7, True
            else:
                drop = False
        else:
            drop = unit() < 0.5
        if not drop:
            packet_id = (0xffffffff - 2 * candidate if schedule == "adversarial"
                         else K + candidate if schedule == "repair-only" else candidate)
            ids.append(packet_id & 0xffffffff)
    return ids, attempted


def validate_cell(row, phase, arm, K, attempt, root_index, root, schedule):
    require(type(row) is dict and set(row) == CELL_KEYS, "cell schema")
    for key, expected in (("type", "cell"), ("phase", phase), ("arm", arm), ("K", K),
                          ("attempt", attempt), ("root_index", root_index),
                          ("root", root), ("schedule", schedule)):
        require(type(row[key]) is type(expected) and row[key] == expected, "cell " + key)
    profile = profile_bytes(K, attempt)
    source_hash = digest(source_bytes(K))
    require(row["profile_hex"] == profile.hex() and row["profile_sha256"] == digest(profile),
            "profile bytes/hash")
    require(row["source_sha256"] == source_hash, "source hash")
    ids, attempted = trace(K, root, schedule)
    require(type(row["ids"]) is list and all(type(x) is int for x in row["ids"]) and
            row["ids"] == ids, "raw packet IDs")
    require(type(row["attempted_candidates"]) is int and row["attempted_candidates"] == attempted,
            "attempted candidates")
    require(row["trace_sha256"] == digest(b"".join(struct.pack("<I", x) for x in ids)), "trace hash")
    outcomes = row["outcomes"]
    require(type(outcomes) is list and len(outcomes) == 4 and
            all(type(x) is int and x in (0, 1, 10) for x in outcomes), "outcome schema")
    if 10 in outcomes:
        require(arm == "candidate" and outcomes == [10] * 4 and row["packets_hex"] is None,
                "construction failure")
    else:
        packets = row["packets_hex"]
        require(type(packets) is str and len(packets) == (K + 4) * 4 and
                all(c in "0123456789abcdef" for c in packets), "packet payload receipt")
        require(all(outcomes[i] >= outcomes[i + 1] for i in range(3)), "non-nested outcomes")
        payload = bytes.fromhex(packets)
        for i, packet_id in enumerate(ids):
            if packet_id < K:
                require(payload[2 * i:2 * i + 2] == source_bytes(K)[2 * packet_id:2 * packet_id + 2],
                        "systematic packet receipt")
    require(row["recovered_sha256"] == [source_hash if x == 0 else None for x in outcomes],
            "recovered byte oracle")


def parse_phase(text, phase, source_commit=None, library_path=None, selected_attempts=None):
    require(phase in ("train", "holdout"), "phase")
    require(type(text) is str and text.endswith("\n") and len(text) <= 32 * 1024 * 1024,
            "raw stream length/termination")
    lines = text.splitlines()
    require(len(lines) >= 2 and all(lines), "raw stream lines")
    records = [parse_json(line) for line in lines]
    begin, terminal = records[0], records[-1]
    require(type(begin) is dict and set(begin) == {"type", "protocol", "phase", "source_commit", "library_path"},
            "begin schema")
    require(begin["type"] == "begin" and begin["phase"] == phase and begin["protocol"] == PROTOCOL,
            "begin identity")
    require(type(begin["source_commit"]) is str and len(begin["source_commit"]) == 40 and
            all(c in "0123456789abcdef" for c in begin["source_commit"]), "source commit format")
    require(type(begin["library_path"]) is str and Path(begin["library_path"]).is_absolute(), "library path")
    if source_commit is not None:
        require(begin["source_commit"] == source_commit, "source commit mismatch")
    if library_path is not None:
        require(begin["library_path"] == str(library_path), "library path mismatch")
    require(type(terminal) is dict and set(terminal) == {"type", "phase", "rows", "selected_attempts"},
            "terminal schema")
    require(terminal["type"] == "terminal" and terminal["phase"] == phase and
            type(terminal["rows"]) is int and terminal["rows"] == len(records) - 2, "terminal count")
    published = terminal["selected_attempts"]
    require(type(published) is list and len(published) == 2 and
            all(type(x) is int and -1 <= x <= 255 for x in published), "terminal map")
    if phase == "holdout":
        require(selected_attempts is not None and list(selected_attempts) == published and
                all(type(x) is int and 0 <= x <= 255 for x in selected_attempts), "sealed holdout map")
    roots = TRAINING_ROOTS if phase == "train" else HOLDOUT_ROOTS
    cells, pos, baselines, selected, duplicates = records[1:-1], 0, [], [], {}

    def arm_rows(K, arm, attempt):
        nonlocal pos
        group = []
        for root_index, root in enumerate(roots):
            for schedule in SCHEDULES:
                require(pos < len(cells), "missing cell")
                row = cells[pos]
                validate_cell(row, phase, arm, K, attempt, root_index, root, schedule)
                identity = (K, attempt, root, schedule)
                semantic = {key: row[key] for key in CELL_KEYS - {"arm"}}
                if identity in duplicates:
                    require(duplicates[identity] == semantic, "same-attempt semantic mismatch")
                duplicates[identity] = semantic
                group.append(row)
                pos += 1
        require(len({10 in row["outcomes"] for row in group}) == 1, "inconsistent constructor result")
        return all(row["outcomes"][0] == 0 for row in group)

    for k_index, K in enumerate(KS):
        require(pos < len(cells) and type(cells[pos]) is dict, "missing baseline")
        baseline = integer(cells[pos].get("attempt"), 0, 255, "baseline attempt")
        baselines.append(baseline)
        arm_rows(K, "baseline", baseline)
        if phase == "train":
            choice = -1
            for attempt in range(256):
                if arm_rows(K, "candidate", attempt):
                    choice = attempt
                    break
            selected.append(choice)
        else:
            selected.append(published[k_index])
            arm_rows(K, "candidate", published[k_index])
    require(pos == len(cells), "extra cells/late selection")
    require(selected == published, "non-first-success map")
    return {"begin": begin, "cells": cells, "selected_attempts": selected,
            "baseline_attempts": baselines}


def paired_counts(parsed, selected):
    baseline, candidate = {}, {}
    for row in parsed["cells"]:
        key = (row["K"], row["root"], row["schedule"])
        if row["arm"] == "baseline":
            baseline[key] = row
        elif row["attempt"] == selected[KS.index(row["K"])]:
            candidate[key] = row
    require(set(baseline) == set(candidate), "paired summary roster")
    result = {"cells_per_arm": len(baseline)}
    for name in ("baseline_failures", "candidate_failures", "repairs", "introductions"):
        result[name] = [0] * 4
    for key, row in baseline.items():
        for h in range(4):
            before, after = row["outcomes"][h] != 0, candidate[key]["outcomes"][h] != 0
            result["baseline_failures"][h] += int(before)
            result["candidate_failures"][h] += int(after)
            result["repairs"][h] += int(before and not after)
            result["introductions"][h] += int(after and not before)
    return result


def summarize(train, holdout):
    selected, baseline = train["selected_attempts"], train["baseline_attempts"]
    result = {"selected_attempts": selected, "baseline_attempts": baseline,
              "map_changed": selected != baseline, "training": None, "holdout": None}
    if -1 in selected:
        require(holdout is None, "holdout after exhausted training")
        result["outcome"] = "EXHAUSTED"
        return result
    require(holdout is not None and holdout["selected_attempts"] == selected and
            holdout["baseline_attempts"] == baseline and
            holdout["begin"]["source_commit"] == train["begin"]["source_commit"] and
            holdout["begin"]["library_path"] == train["begin"]["library_path"], "phase continuity")
    result["training"] = paired_counts(train, selected)
    result["holdout"] = paired_counts(holdout, selected)
    result["outcome"] = ("FAIL" if result["holdout"]["candidate_failures"][0]
                         else "PASS" if result["map_changed"] else "NO_CHANGE")
    return result


def command(argv, timeout=30):
    return subprocess.check_output(argv, cwd=str(ROOT), env={"PATH": "/usr/bin:/bin", "LC_ALL": "C"},
                                   timeout=timeout, stderr=subprocess.STDOUT).decode("utf-8")


def source_receipt():
    paths = command(["git", "ls-files", "--", "CMakeLists.txt", "cmake", "codec", "include", "tables",
                     "wirehair.cpp", "Wirehair*", "gf256.cpp", "gf256.h", "bench/Wh2FrozenTrace.*",
                     "bench/Wh2ProductionMix3RecoveryR0.*", "bench/test_Wh2ProductionMix3RecoveryR0.py"]).splitlines()
    # Include the additive files during permitted precommit selftests as well.
    paths = sorted(set(paths) | {str(WORKER_SOURCE.relative_to(ROOT)),
                                "bench/Wh2ProductionMix3RecoveryR0.py"})
    return {path: file_digest(ROOT / path) for path in paths}


def write_json(path, value):
    with Path(path).open("xb") as output:
        output.write(canonical(value) + b"\n")


def build(output, library, compiler, library_build):
    output, library = Path(output).resolve(), Path(library).resolve(strict=True)
    compiler, library_build = Path(compiler).resolve(strict=True), Path(library_build).resolve(strict=True)
    require(not output.exists(), "build output already exists")
    dry_run = command(["/usr/bin/ninja", "-C", str(library_build), "-n", "wirehair"])
    require(dry_run.splitlines()[-1] == "ninja: no work to do.", "production library build is stale")
    require(library.parent == library_build, "library/build directory mismatch")
    sources = source_receipt()
    commit = command(["git", "rev-parse", "HEAD"]).strip()
    output.mkdir()
    worker = output / "worker"
    argv = [str(compiler), "-std=c++11", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(commit),
            "-I" + str(ROOT / "include"), str(WORKER_SOURCE), str(ROOT / "bench/Wh2FrozenTrace.cpp"),
            str(library), "-Wl,-rpath," + str(library.parent), "-ldl", "-o", str(worker)]
    build_output = command(argv, timeout=60)
    require(source_receipt() == sources, "sources changed during build")
    description = parse_json(command([str(worker), "--describe"]))
    require(description == {"type": "begin", "protocol": PROTOCOL, "phase": "describe",
                            "source_commit": commit, "library_path": str(library)}, "worker binding")
    receipt = {"protocol": PROTOCOL, "source_commit": commit, "source_files": sources,
               "compiler": str(compiler), "compiler_sha256": file_digest(compiler),
               "compiler_version": command([str(compiler), "--version"]), "command": argv,
               "build_output": build_output, "library": str(library), "library_sha256": file_digest(library),
               "library_build": str(library_build), "library_dry_run": dry_run,
               "library_cache_sha256": file_digest(library_build / "CMakeCache.txt"),
               "library_commands": command(["/usr/bin/ninja", "-C", str(library_build), "-t", "commands", "wirehair"]),
               "worker": str(worker), "worker_sha256": file_digest(worker)}
    write_json(output / "build.json", receipt)
    return receipt


def check_build(receipt):
    require(receipt["protocol"] == PROTOCOL, "build protocol")
    require(command(["git", "rev-parse", "HEAD"]).strip() == receipt["source_commit"], "build commit drift")
    require(source_receipt() == receipt["source_files"], "build source drift")
    for path_key, hash_key in (("worker", "worker_sha256"), ("library", "library_sha256"),
                              ("compiler", "compiler_sha256")):
        require(file_digest(receipt[path_key]) == receipt[hash_key], path_key + " changed")


def execute_worker(argv, raw_path, error_path, timeout):
    with raw_path.open("xb") as stdout, error_path.open("xb") as stderr:
        process = subprocess.Popen(argv, cwd=str(ROOT), env={"PATH": "/usr/bin:/bin", "LC_ALL": "C"},
                                   stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            result = process.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise
    require(result == 0 and error_path.stat().st_size == 0, "worker failed")


def frozen_protocol(receipt):
    return {"protocol": PROTOCOL, "build": receipt, "training_roots": list(TRAINING_ROOTS),
            "holdout_roots": list(HOLDOUT_ROOTS), "K": list(KS), "block_bytes": 2,
            "schedules": list(SCHEDULES), "loss_ppm": 500000, "overheads": list(OVERHEADS),
            "worker_budget_seconds": 60, "outer_budget_seconds": 70, "max_prefix_decodes": 37440}


def publish_result(bundle, summary, started):
    write_json(bundle / "summary.json", summary)
    hashes = {p.name: file_digest(p) for p in sorted(bundle.iterdir()) if p.is_file()}
    elapsed = time.monotonic() - started
    if elapsed >= 70:
        # This summary is ours and is not yet sealed by COMPLETE.  Never publish
        # an earlier scientific result after hashing/publication exhausts time.
        summary = {"outcome": "INVALID", "error": "outer deadline during publication",
                   "protocol": PROTOCOL, "elapsed_seconds": elapsed,
                   "promotion_claimed": False, "all_K_claimed": False, "timing_claimed": False}
        (bundle / "summary.json").write_bytes(canonical(summary) + b"\n")
        hashes["summary.json"] = file_digest(bundle / "summary.json")
    write_json(bundle / "COMPLETE", hashes)
    return summary


def run_once(build_path, bundle_path):
    receipt = parse_json(Path(build_path).read_bytes())
    check_build(receipt)
    require(command(["git", "status", "--porcelain", "--untracked-files=no"]).strip() == "", "tracked tree dirty")
    for source in receipt["source_files"]:
        command(["git", "ls-files", "--error-unmatch", "--", source])
    # The source must already be published before an experiment consumes its roster.
    require(command(["git", "rev-parse", "@{upstream}"]).strip() == receipt["source_commit"], "source not pushed")
    bundle = Path(bundle_path).resolve()
    bundle.mkdir(exist_ok=False)
    freeze = frozen_protocol(receipt)
    write_json(bundle / "freeze.json", freeze)
    started = time.monotonic()
    summary = None
    try:
        execute_worker([receipt["worker"], "--train", "60"], bundle / "train.jsonl", bundle / "train.stderr", 60)
        train = parse_phase((bundle / "train.jsonl").read_text(), "train", receipt["source_commit"], receipt["library"])
        selected = train["selected_attempts"]
        selection = {"protocol": PROTOCOL, "selected_attempts": selected,
                     "train_sha256": file_digest(bundle / "train.jsonl")}
        write_json(bundle / "selection.json", selection)
        selection_hash = file_digest(bundle / "selection.json")
        holdout = None
        if -1 not in selected:
            remaining = 60.0 - (time.monotonic() - started)
            require(remaining >= 1.0, "aggregate worker budget exhausted")
            seconds = min(60, int(math.floor(remaining)))
            execute_worker([receipt["worker"], "--holdout", str(seconds)] + [str(x) for x in selected],
                           bundle / "holdout.jsonl", bundle / "holdout.stderr", remaining)
            holdout = parse_phase((bundle / "holdout.jsonl").read_text(), "holdout", receipt["source_commit"],
                                  receipt["library"], selected)
        require(file_digest(bundle / "selection.json") == selection_hash, "sealed map changed")
        summary = summarize(train, holdout)
        check_build(receipt)
        require(time.monotonic() - started < 70, "outer deadline")
    except Exception as error:
        summary = {"outcome": "INVALID", "error": str(error)}
    summary.update({"protocol": PROTOCOL, "elapsed_seconds": time.monotonic() - started,
                    "promotion_claimed": False, "all_K_claimed": False, "timing_claimed": False})
    return publish_result(bundle, summary, started)


def replay(bundle_path):
    bundle = Path(bundle_path)
    complete = parse_json((bundle / "COMPLETE").read_bytes())
    require(type(complete) is dict and set(complete) == {p.name for p in bundle.iterdir() if p.name != "COMPLETE"},
            "bundle files")
    for name, expected in complete.items():
        require(Path(name).name == name and file_digest(bundle / name) == expected, "bundle digest")
    freeze = parse_json((bundle / "freeze.json").read_bytes())
    require(type(freeze) is dict and type(freeze.get("build")) is dict and
            canonical(freeze) == canonical(frozen_protocol(freeze["build"])), "frozen protocol")
    summary = parse_json((bundle / "summary.json").read_bytes())
    require(summary["protocol"] == PROTOCOL and summary["promotion_claimed"] is False and
            summary["all_K_claimed"] is False and summary["timing_claimed"] is False, "summary scope")
    if summary["outcome"] == "INVALID":
        return summary  # Partial infrastructure receipts make no scientific claim.
    require(summary["outcome"] in ("PASS", "FAIL", "NO_CHANGE", "EXHAUSTED"), "unknown outcome")
    elapsed = summary.get("elapsed_seconds")
    require(type(elapsed) in (int, float) and math.isfinite(elapsed) and 0 <= elapsed < 70,
            "elapsed bound")
    expected_files = {"freeze.json", "train.jsonl", "train.stderr", "selection.json", "summary.json"}
    if summary["outcome"] != "EXHAUSTED":
        expected_files.update(("holdout.jsonl", "holdout.stderr"))
    require(set(complete) == expected_files, "scientific bundle file roster")
    for name in expected_files:
        if name.endswith(".stderr"):
            require((bundle / name).read_bytes() == b"", "worker stderr")
    build = freeze["build"]
    train = parse_phase((bundle / "train.jsonl").read_text(), "train", build["source_commit"], build["library"])
    selection = parse_json((bundle / "selection.json").read_bytes())
    require(selection == {"protocol": PROTOCOL, "selected_attempts": train["selected_attempts"],
                          "train_sha256": file_digest(bundle / "train.jsonl")}, "selection receipt")
    holdout = None
    if -1 not in train["selected_attempts"]:
        holdout = parse_phase((bundle / "holdout.jsonl").read_text(), "holdout", build["source_commit"],
                              build["library"], train["selected_attempts"])
    expected = summarize(train, holdout)
    expected_keys = set(expected) | {"protocol", "elapsed_seconds", "promotion_claimed", "all_K_claimed", "timing_claimed"}
    require(set(summary) == expected_keys and
            all(canonical(summary.get(key)) == canonical(value) for key, value in expected.items()),
            "summary mismatch")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build", metavar="NEW_DIRECTORY")
    modes.add_argument("--run", metavar="NEW_BUNDLE")
    modes.add_argument("--replay", metavar="BUNDLE")
    parser.add_argument("--library")
    parser.add_argument("--library-build")
    parser.add_argument("--compiler", default="/usr/bin/g++-13")
    parser.add_argument("--build-receipt")
    args = parser.parse_args()
    if args.build:
        require(args.library and args.library_build, "build needs library and library-build")
        result = build(args.build, args.library, args.compiler, args.library_build)
    elif args.run:
        require(args.build_receipt, "run needs build-receipt")
        result = run_once(args.build_receipt, args.run)
    else:
        result = replay(args.replay)
    print(canonical(result).decode("ascii"))


if __name__ == "__main__":
    try:
        main()
    except (ValidationError, OSError, subprocess.SubprocessError) as error:
        print("production MIX3 recovery r0: " + str(error), file=sys.stderr)
        sys.exit(2)
