#!/usr/bin/env python3
"""Bounded joint-width MIX3 first-success screen; development evidence only.

Replay checks a trusted producer's receipts, not authenticity or omitted codec
calls. Rejected attempts retain their first failure; selected attempts are fully
revalidated before a separately invoked, never-retuned fresh holdout.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import math
import os
from pathlib import Path
import selectors
import signal
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
WORKER_SOURCE = Path(__file__).resolve()
PROTOCOL = "wirehair.wh2.production-mix3-k3k6-joint-r0"
NAMESPACE = Path("/var/tmp/wh2-production-mix3-k3k6-joint-r0")
spec = importlib.util.spec_from_file_location("_wh2_joint_portability_private",
                                             ROOT / "bench/Wh2ProductionMix3RecoveryPortabilityR0.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
b, c = p.b, p.c
c.PROTOCOL = PROTOCOL
KS, WIDTHS = c.KS, p.WIDTHS
STREAM_LIMIT, TOTAL_LIMIT = 33554432, 67108864
ROSTER_SHA256 = "8245376b0c91593e099f35207b436a0ed6343fd212bc6119f6357d74c1fd2470"
LEDGER_SHA256 = "e12d079b994aa04d7dd47656cf42895c1910f4239434ee8f3a2aae9bb29e0b1b"
ORIGINS_SHA256 = "5362613e8ade02cec8932fad0fc1fd0e290eaa343b2eb11471ebd79a094034ff"
BAD_ORIGINS_SHA256 = "617c0db46b2a719d14f780f3ed44872c5ecf65d6f499e2500d3438d89dcf4613"
PRIOR_ROOTS_SHA256 = "94a2b0733dca059c3ff439799c81a8e5513244031d1720d550db98e2ed9090ea"
BAD_PROFILES = [[3, 2, x] for x in range(1, 5)] + [[6, 64, 3]]
HISTORY = (
    (52, "wh2-production-mix3-k3k6-r0", "2963be04d328988a9f7d2e363295cbea57fe76ca2d5585125341f2cceb46f993"),
    (54, "wh2-production-mix3-k3k6-broad-r0", "105243293b25c76b749fa4df7f0f3970c5cf3c0378f88d8c5a767c3b2c6e6734"),
    (55, "wh2-production-mix3-k3k6-portability-r0", "f1dc2656c26337edfcd7994b9d6204ef2ee80af70599e69a85d604a19e1b3996"))
RESERVED = ("0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3",
            "0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b",
            "0xefd20c982041a46b", "0x8827bc36ed906555", "0x86029f23d6132efa")
TRAINING_ROOTS, HOLDOUT_ROOTS = (tuple("0x" + c.digest((PROTOCOL + ":" + phase + "/" + str(i)).encode("ascii"))[:16]
                                    for i in range(count)) for phase, count in (("train", 16), ("holdout", 64)))
SOURCE_PATHS = tuple(sorted(set(p.SOURCE_PATHS) - {"bench/test_Wh2ProductionMix3RecoveryPortabilityR0.py"} |
                            {"bench/Wh2ProductionMix3RecoveryJointR0.py",
                             "bench/test_Wh2ProductionMix3RecoveryJointR0.py"}))
IMMUTABLE_HELPERS = {
    "bench/Wh2ProductionMix3RecoveryPortabilityR0.py": "559a377ce8d2cc0c8bdef1f8ec5d4b2198b7b9657ce18df780672066aab6afa1",
    "bench/Wh2ProductionMix3RecoveryBroadR0.py": "dfecce6141de2371d5fa818ab0eabaf1659ea557db71f06febf597a5b57ddacc",
    "bench/Wh2ProductionMix3RecoveryR0.py": "ea6768b7162919fa9546b1e87626fad3e89b644151947910f77c505cb7b1c7d0"}


def exact(actual, expected, label):
    c.require(c.canonical(actual) == c.canonical(expected), label)


def hex_value(value, length, label):
    c.require(type(value) is str and len(value) == length and
              all(x in "0123456789abcdef" for x in value), label)


def validate_history(history):
    c.require(type(history) is dict and set(history) == {"prefixes", "origins", "bad_origins", "prior_roots"},
              "history schema")
    for key, count, digest in (("prefixes", 73, LEDGER_SHA256), ("origins", 126, ORIGINS_SHA256),
                               ("bad_origins", 360, BAD_ORIGINS_SHA256)):
        c.require(type(history[key]) is list and len(history[key]) == count and
                  c.digest(c.canonical(history[key])) == digest, key + " ledger drift")
    exact(sorted({(row["K"], row["B"], row["attempt"]) for row in history["bad_origins"]}),
          BAD_PROFILES, "known construction profiles")
    prior = history["prior_roots"]
    c.require(type(prior) is list and len(prior) == 127 and prior == sorted(set(prior)) and
              c.digest(c.canonical(prior)) == PRIOR_ROOTS_SHA256, "prior root roster")
    for root in prior:
        c.require(type(root) is str and root.startswith("0x"), "prior root format")
        hex_value(root[2:], 16, "prior root format")
    roots = TRAINING_ROOTS + HOLDOUT_ROOTS
    c.require(c.digest(c.canonical({"train": TRAINING_ROOTS, "holdout": HOLDOUT_ROOTS})) == ROSTER_SHA256 and
              len(set(roots)) == 80 and not set(roots) & (set(prior) | set(RESERVED)), "root disjointness/hash")


def load_history():
    origins, bad_origins, prefixes, prior = [], [], {}, set()
    for issue, name, manifest_sha in HISTORY:
        bundle = Path("/var/tmp") / name
        c.require(c.file_digest(bundle / "COMPLETE") == manifest_sha, "historical manifest drift")
        manifest = c.parse_json((bundle / "COMPLETE").read_bytes())
        for filename, digest in manifest.items():
            c.require(Path(filename).name == filename and c.file_digest(bundle / filename) == digest,
                      "historical artifact drift")
        freeze = c.parse_json((bundle / "freeze.json").read_bytes())
        prior.update(freeze.get("training_roots", []))
        prior.update(freeze["holdout_roots"])
        for filename in sorted(x for x in manifest if x.endswith(".jsonl")):
            B = int(filename[1:-6]) if issue == 55 else 2
            for line, raw in enumerate((bundle / filename).read_bytes().splitlines(), 1):
                row = c.parse_json(raw)
                if row["type"] != "cell":
                    continue
                origin = dict(issue=issue, file=filename, line=line, K=row["K"], B=B,
                              arm=row["arm"], attempt=row["attempt"], root=row["root"], schedule=row["schedule"])
                if 10 in row["outcomes"]:
                    bad_origins.append(origin)
                for overhead, outcome in zip(c.OVERHEADS, row["outcomes"]):
                    if outcome == 1:
                        ids = row["ids"][:row["K"] + overhead]
                        origins.append(dict(origin, overhead=overhead, ids=ids))
                        prefixes.setdefault((row["K"], tuple(ids)), set()).add(B)
    ledger = [dict(K=K, ids=list(ids), original_widths=sorted(widths), overhead=len(ids) - K)
              for (K, ids), widths in sorted(prefixes.items())]
    history = dict(prefixes=ledger, origins=origins, bad_origins=bad_origins, prior_roots=sorted(prior))
    validate_history(history)
    return history


def source_receipt():
    return {path: c.file_digest(ROOT / path) for path in SOURCE_PATHS}


def validate_receipt(receipt):
    keys = {"protocol", "source_commit", "source_files", "worker", "worker_sha256", "interpreter",
            "interpreter_sha256", "interpreter_version", "interpreter_flags", "library", "library_sha256",
            "production_source_commit", "production_build"}
    c.require(type(receipt) is dict and set(receipt) == keys, "build schema")
    b.validate_basis(receipt["production_build"])
    c.require(receipt["protocol"] == PROTOCOL and receipt["production_source_commit"] == b.PRODUCTION_COMMIT and
              receipt["library"] == str(b.PINNED_LIBRARY) and receipt["library_sha256"] == b.LIBRARY_SHA256,
              "production identity")
    hex_value(receipt["source_commit"], 40, "source commit")
    c.require(type(receipt["source_files"]) is dict and set(receipt["source_files"]) == set(SOURCE_PATHS),
              "source closure")
    for digest in receipt["source_files"].values():
        hex_value(digest, 64, "source digest")
    for path, digest in IMMUTABLE_HELPERS.items():
        c.require(receipt["source_files"][path] == digest, "immutable helper drift: " + path)
    c.require(receipt["worker"] == str(WORKER_SOURCE) and
              receipt["worker_sha256"] == receipt["source_files"][str(WORKER_SOURCE.relative_to(ROOT))] and
              type(receipt["interpreter"]) is str and Path(receipt["interpreter"]).is_absolute() and
              type(receipt["interpreter_version"]) is str and bool(receipt["interpreter_version"]), "worker/interpreter")
    hex_value(receipt["interpreter_sha256"], 64, "interpreter digest")
    exact(receipt["interpreter_flags"], ["-I", "-B", "-S"], "interpreter flags")
    for path in set(SOURCE_PATHS) & set(receipt["production_build"]["source_files"]):
        c.require(receipt["source_files"][path] == receipt["production_build"]["source_files"][path],
                  "inherited source drift")


def check_build(receipt):
    validate_receipt(receipt)
    c.require(c.command(["git", "rev-parse", "HEAD"]).strip() == receipt["source_commit"], "build commit drift")
    exact(source_receipt(), receipt["source_files"], "build source drift")
    for key in ("worker", "library", "interpreter"):
        c.require(c.file_digest(receipt[key]) == receipt[key + "_sha256"], key + " changed")


def build(output):
    output = Path(output).resolve()
    c.require(not output.exists(), "build output exists")
    basis, history = b.pinned_library_receipt(), load_history()
    interpreter = str(Path(sys.executable).resolve())
    receipt = dict(protocol=PROTOCOL, source_commit=c.command(["git", "rev-parse", "HEAD"]).strip(),
                   source_files=source_receipt(), worker=str(WORKER_SOURCE), worker_sha256=c.file_digest(WORKER_SOURCE),
                   interpreter=interpreter, interpreter_sha256=c.file_digest(interpreter),
                   interpreter_version=sys.version, interpreter_flags=["-I", "-B", "-S"],
                   library=str(b.PINNED_LIBRARY), library_sha256=b.LIBRARY_SHA256,
                   production_source_commit=b.PRODUCTION_COMMIT, production_build=basis)
    check_build(receipt)
    output.mkdir()
    c.write_json(output / "build.json", receipt)
    c.write_json(output / "inputs.json", history)
    return receipt


def frozen_protocol(receipt, history):
    validate_receipt(receipt)
    validate_history(history)
    return dict(protocol=PROTOCOL, build=receipt, history=history, K=list(KS), widths=list(WIDTHS),
                training_roots=list(TRAINING_ROOTS), holdout_roots=list(HOLDOUT_ROOTS), roster_sha256=ROSTER_SHA256,
                ledger_sha256=LEDGER_SHA256, origins_sha256=ORIGINS_SHA256, bad_origins_sha256=BAD_ORIGINS_SHA256,
                bad_profiles=BAD_PROFILES, schedules=list(c.SCHEDULES), loss_ppm=500000, overheads=list(c.OVERHEADS),
                worker_budget_seconds=60, outer_budget_seconds=70, stream_limit_exclusive=STREAM_LIMIT,
                aggregate_limit_inclusive=TOTAL_LIMIT, max_prefix_cases=141750, max_constructors=1560)


def validate_freeze(freeze):
    c.require(type(freeze) is dict and "build" in freeze and "history" in freeze, "freeze schema")
    exact(freeze, frozen_protocol(freeze["build"], freeze["history"]), "frozen protocol")


def case_specs(K, kind, roots, history):
    for B in WIDTHS:
        if kind == "regression":
            for index, row in enumerate(history["prefixes"]):
                if row["K"] == K:
                    yield dict(B=B, kind=kind, index=index, root=None, schedule=None,
                               ids=row["ids"], attempted_candidates=None, overheads=[row["overhead"]])
        else:
            c.require(kind == "fresh", "case kind")
            for root_index, root in enumerate(roots):
                for schedule_index, schedule in enumerate(c.SCHEDULES):
                    ids, attempted = p.trace(K, B, root, schedule)
                    yield dict(B=B, kind=kind, index=root_index * 3 + schedule_index, root=root, schedule=schedule,
                               ids=ids, attempted_candidates=attempted, overheads=list(c.OVERHEADS))


def specs_for(K, phase, history):
    roots = HOLDOUT_ROOTS if phase == "holdout" else TRAINING_ROOTS
    result = [] if phase == "holdout" else list(case_specs(K, "regression", roots, history))
    return result + list(case_specs(K, "fresh", roots, history))


def prefix_cell(row):
    result = dict(row)
    count = row["K"] + row["overheads"][0]
    result.update(ids=row["ids"][:count], overheads=row["overheads"][:1], outcomes=row["outcomes"][:1],
                  recovered_sha256=row["recovered_sha256"][:1],
                  packets_hex=None if row["packets_hex"] is None else row["packets_hex"][:count * row["B"] * 2])
    result["trace_sha256"] = c.digest(struct.pack("<" + "I" * count, *result["ids"]))
    return result


def success_token(row):
    result = prefix_cell(row)
    result["packets_sha256"] = c.digest(bytes.fromhex(result.pop("packets_hex")))
    return c.canonical(result) + b"\n"


def check_time(deadline):
    c.require(time.monotonic() < deadline, "worker deadline")


def make_cell(native, state, K, B, arm, spec, overheads, deadline):
    handle, profile, message, source, bad, actual = state
    ids = spec["ids"][:K + max(overheads)]
    check_time(deadline)
    packets = None if bad else b"".join(native.encode(handle, packet_id, B) for packet_id in ids)
    outcomes = []
    for overhead in overheads:
        check_time(deadline)
        outcomes.append(10 if bad else native.decode(profile, source, B, ids, packets, K + overhead))
    c.require(message.raw == source, "source mutated during packet/decode work")
    source_hash = c.digest(source)
    return dict(spec, type="cell", K=K, arm=arm, attempt=actual, ids=ids, overheads=list(overheads),
                profile_hex=profile.hex(), profile_sha256=c.digest(profile), source_sha256=source_hash,
                trace_sha256=c.digest(struct.pack("<" + "I" * len(ids), *ids)),
                packets_hex=None if bad else packets.hex(), outcomes=outcomes,
                recovered_sha256=[source_hash if outcome == 0 else None for outcome in outcomes])


def validate_cell(row, K, arm, attempt, spec, search=False):
    overheads = spec["overheads"][:1] if search else spec["overheads"]
    ids = spec["ids"][:K + max(overheads)]
    expected = dict(spec, type="cell", K=K, arm=arm, attempt=attempt, ids=ids, overheads=overheads)
    c.require(type(row) is dict and set(row) == set(expected) | {
        "profile_hex", "profile_sha256", "source_sha256", "trace_sha256", "packets_hex", "outcomes", "recovered_sha256"},
        "cell schema")
    for key, value in expected.items():
        exact(row[key], value, "cell " + key)
    profile, source = p.profile_bytes(K, spec["B"], attempt), p.source_bytes(K, spec["B"])
    c.require(row["profile_hex"] == profile.hex() and row["profile_sha256"] == c.digest(profile) and
              row["source_sha256"] == c.digest(source) and
              row["trace_sha256"] == c.digest(struct.pack("<" + "I" * len(ids), *ids)), "cell byte identities")
    outcomes = row["outcomes"]
    c.require(type(outcomes) is list and len(outcomes) == len(overheads) and
              all(type(x) is int and x in (0, 1, 10) for x in outcomes), "outcomes")
    if [K, spec["B"], attempt] in BAD_PROFILES:
        c.require(outcomes == [10] * len(overheads), "known BadSeed profile contradicted")
    if 10 in outcomes:
        c.require(arm == "candidate" and outcomes == [10] * len(overheads) and row["packets_hex"] is None,
                  "BadSeed shape")
    else:
        hex_value(row["packets_hex"], len(ids) * spec["B"] * 2, "packet bytes")
        payload, B = bytes.fromhex(row["packets_hex"]), spec["B"]
        for i, packet_id in enumerate(ids):
            if packet_id < K:
                c.require(payload[i * B:(i + 1) * B] == source[packet_id * B:(packet_id + 1) * B],
                          "systematic packet bytes")
        c.require(all(outcomes[i] >= outcomes[i + 1] for i in range(len(outcomes) - 1)), "nested outcomes")
    exact(row["recovered_sha256"], [c.digest(source) if x == 0 else None for x in outcomes], "recovered oracle")


def search_attempt(native, K, attempt, specs, deadline):
    constructors, passed, witness, digest = [], 0, None, hashlib.sha256()
    with contextlib.ExitStack() as stack:
        states = {}
        for B in WIDTHS:
            check_time(deadline)
            state = stack.enter_context(native.encoder(K, B, attempt))
            states[B] = state
            constructors.append(dict(B=B, profile_hex=state[1].hex(), outcome=10 if state[4] else 0))
            if state[4]:
                break
        if constructors[-1]["outcome"] == 0:
            for spec in specs:
                row = make_cell(native, states[spec["B"]], K, spec["B"], "candidate", spec,
                                spec["overheads"][:1], deadline)
                if row["outcomes"] != [0]:
                    witness = row
                    break
                digest.update(success_token(row))
                passed += 1
    return dict(type="attempt", K=K, attempt=attempt, constructors=constructors, passed=passed,
                pass_sha256=digest.hexdigest(), witness=witness, accepted=passed == len(specs))


def begin_record(phase, freeze, selection_sha256):
    return dict(type="begin", protocol=PROTOCOL, phase=phase, source_commit=freeze["build"]["source_commit"],
                library_path=freeze["build"]["library"], selection_sha256=selection_sha256)


def worker(phase, freeze_path, seconds, selection_path=None, selection_sha256=None):
    c.require(phase in ("search", "validate", "holdout"), "worker phase")
    c.require(type(seconds) in (int, float) and math.isfinite(seconds) and 0 < seconds <= 60, "worker seconds")
    deadline = time.monotonic() + seconds
    freeze = c.parse_json(Path(freeze_path).read_bytes())
    validate_freeze(freeze)
    check_build(freeze["build"])
    selected = [-1, -1]
    if phase != "search":
        c.require(selection_path is not None and c.file_digest(selection_path) == selection_sha256, "sealed map hash")
        selected = validate_selection(c.parse_json(Path(selection_path).read_bytes()))
        c.require(-1 not in selected, "exhausted worker map")
    else:
        c.require(selection_path is None and selection_sha256 is None, "search map supplied")
    native, rows = p.Native(), 0
    p.emit(begin_record(phase, freeze, selection_sha256))
    for k_index, K in enumerate(KS):
        specs = specs_for(K, phase, freeze["history"])
        if phase == "search":
            for attempt in range(256):
                row = search_attempt(native, K, attempt, specs, deadline)
                p.emit(row)
                rows += 1
                if row["accepted"]:
                    selected[k_index] = attempt
                    break
        else:
            for arm, attempt in (("baseline", None), ("candidate", selected[k_index])):
                with contextlib.ExitStack() as stack:
                    states = {}
                    for B in WIDTHS:
                        check_time(deadline)
                        states[B] = stack.enter_context(native.encoder(K, B, attempt))
                    for spec in specs:
                        row = make_cell(native, states[spec["B"]], K, spec["B"], arm, spec, spec["overheads"], deadline)
                        p.emit(row)
                        rows += 1
    check_time(deadline)
    p.emit(dict(type="terminal", phase=phase, rows=rows, selected_attempts=selected))


def records(text, phase, freeze, selection_sha256):
    c.require(type(text) is str and text.endswith("\n") and len(text.encode("ascii")) < STREAM_LIMIT,
              "stream length/termination")
    rows = [c.parse_json(line) for line in text.splitlines()]
    c.require(len(rows) >= 2, "stream records")
    exact(rows[0], begin_record(phase, freeze, selection_sha256), "begin")
    terminal = rows[-1]
    c.require(type(terminal) is dict and set(terminal) == {"type", "phase", "rows", "selected_attempts"}, "terminal schema")
    exact({key: terminal[key] for key in ("type", "phase", "rows")},
          dict(type="terminal", phase=phase, rows=len(rows) - 2), "terminal count")
    selected = terminal["selected_attempts"]
    c.require(type(selected) is list and len(selected) == 2 and all(type(x) is int and -1 <= x <= 255 for x in selected),
              "terminal map")
    return rows[1:-1], selected


def parse_search(text, freeze):
    rows, selected = records(text, "search", freeze, None)
    position, choices, prefix_cases, constructors = 0, [], 0, 0
    for K in KS:
        specs, choice = specs_for(K, "search", freeze["history"]), -1
        for attempt in range(256):
            c.require(position < len(rows), "missing attempt")
            row = rows[position]
            position += 1
            c.require(type(row) is dict and set(row) == {"type", "K", "attempt", "constructors", "passed",
                                                        "pass_sha256", "witness", "accepted"}, "attempt schema")
            exact([row["type"], row["K"], row["attempt"]], ["attempt", K, attempt], "attempt order")
            cons = row["constructors"]
            c.require(type(cons) is list and 1 <= len(cons) <= 3, "constructor roster")
            for index, con in enumerate(cons):
                c.require(type(con) is dict and set(con) == {"B", "profile_hex", "outcome"} and
                          type(con["outcome"]) is int and con["outcome"] in (0, 10), "constructor schema")
                exact([con["B"], con["profile_hex"]], [WIDTHS[index], p.profile_bytes(K, WIDTHS[index], attempt).hex()],
                      "constructor profile/order")
                if [K, WIDTHS[index], attempt] in BAD_PROFILES:
                    c.require(con["outcome"] == 10, "known BadSeed constructor contradicted")
                c.require(con["outcome"] == 0 or index == len(cons) - 1, "constructor after failure")
            constructors += len(cons)
            bad = cons[-1]["outcome"] == 10
            c.require(bad or len(cons) == 3, "missing successful constructor")
            passed = c.integer(row["passed"], 0, len(specs), "passed constraints")
            hex_value(row["pass_sha256"], 64, "passed hash")
            if passed == 0:
                c.require(row["pass_sha256"] == c.digest(b""), "empty passed hash")
            c.require(type(row["accepted"]) is bool, "accepted type")
            if bad:
                c.require(passed == 0 and row["witness"] is None and not row["accepted"], "construction witness")
            elif passed == len(specs):
                c.require(row["witness"] is None and row["accepted"], "accepted attempt")
                choice = attempt
            else:
                validate_cell(row["witness"], K, "candidate", attempt, specs[passed], search=True)
                c.require(row["witness"]["outcomes"] == [1] and not row["accepted"], "first failure witness")
                prefix_cases += 1
            prefix_cases += passed
            if choice != -1:
                break
        choices.append(choice)
    c.require(position == len(rows), "extra attempts/late selection")
    exact(selected, choices, "first-success map")
    c.require(prefix_cases <= 129792 and constructors <= 1536, "search bounds")
    return dict(selected_attempts=choices, records=rows, prefix_cases=prefix_cases, constructors=constructors)


def validate_selection(selection, search=None, search_sha256=None):
    c.require(type(selection) is dict and set(selection) == {"protocol", "selected_attempts", "search_sha256"} and
              selection["protocol"] == PROTOCOL, "selection schema")
    selected = selection["selected_attempts"]
    c.require(type(selected) is list and len(selected) == 2 and
              all(type(x) is int and -1 <= x <= 255 for x in selected), "selection attempts")
    hex_value(selection["search_sha256"], 64, "search digest")
    if search is not None:
        exact(selected, search["selected_attempts"], "sealed search map")
        c.require(selection["search_sha256"] == search_sha256, "sealed search hash")
    return selected


def parse_validation(text, phase, freeze, selection, search):
    c.require(phase in ("validate", "holdout"), "validation phase")
    selection_sha = c.digest(c.canonical(selection) + b"\n")
    rows, selected = records(text, phase, freeze, selection_sha)
    exact(selected, selection["selected_attempts"], "validation map")
    c.require(-1 not in selected, "validation after exhaustion")
    position, baselines, semantic = 0, {}, {}
    for k_index, K in enumerate(KS):
        specs = specs_for(K, phase, freeze["history"])
        for arm in ("baseline", "candidate"):
            digest, attempts, bad_by_width = hashlib.sha256(), {}, {}
            for spec in specs:
                c.require(position < len(rows), "missing validation cell")
                row, B = rows[position], spec["B"]
                position += 1
                c.require(type(row) is dict, "validation row")
                attempt = selected[k_index] if arm == "candidate" else c.integer(row.get("attempt"), 0, 255, "baseline attempt")
                validate_cell(row, K, arm, attempt, spec)
                c.require(attempts.setdefault(B, attempt) == attempt, "constructor attempt drift")
                bad = row["outcomes"][0] == 10
                c.require(bad_by_width.setdefault(B, bad) == bad, "constructor status drift")
                identity = (K, B, attempt, spec["kind"], spec["index"])
                value = dict(row, arm="candidate")
                if identity in semantic:
                    exact(value, semantic[identity], "same attempt semantic mismatch")
                semantic[identity] = value
                if arm == "candidate" and phase == "validate":
                    c.require(all(x == 0 for x in row["outcomes"]), "selected revalidation failure")
                    digest.update(success_token(row))
            if arm == "baseline":
                baselines[str(K)] = [attempts[B] for B in WIDTHS]
            elif phase == "validate":
                accepted = next(row for row in search["records"] if row["K"] == K and row["accepted"])
                c.require(digest.hexdigest() == accepted["pass_sha256"], "selected search/revalidation byte mismatch")
    c.require(position == len(rows), "extra validation cells")
    return dict(rows=rows, baseline_attempts=baselines)


def paired(rows, B, kind):
    arms = {arm: [row for row in rows if row["B"] == B and row["kind"] == kind and row["arm"] == arm]
            for arm in ("baseline", "candidate")}
    pairs = list(zip(arms["baseline"], arms["candidate"]))
    c.require(len(arms["baseline"]) == len(arms["candidate"]), "unpaired cells")
    overheads = sorted({oh for row in arms["baseline"] for oh in row["overheads"]})
    counts = []
    for overhead in overheads:
        outcomes = [(base["outcomes"][base["overheads"].index(overhead)], cand["outcomes"][cand["overheads"].index(overhead)])
                    for base, cand in pairs if overhead in base["overheads"]]
        counts.append(dict(overhead=overhead, pairs=len(outcomes), baseline_failures=sum(a != 0 for a, _ in outcomes),
                           candidate_failures=sum(z != 0 for _, z in outcomes), repairs=sum(a != 0 and z == 0 for a, z in outcomes),
                           introductions=sum(a == 0 and z != 0 for a, z in outcomes)))
    return counts


def summarize(search, validation=None, holdout=None):
    selected = search["selected_attempts"]
    result = dict(selected_attempts=selected, search_prefix_cases=search["prefix_cases"],
                  search_constructors=search["constructors"], total_prefix_cases=search["prefix_cases"],
                  total_constructors=search["constructors"], widths=[])
    if -1 in selected:
        c.require(validation is None and holdout is None, "validation after exhaustion")
        result["outcome"] = "EXHAUSTED"
        return result
    c.require(validation is not None and holdout is not None, "missing validation/holdout")
    exact(validation["baseline_attempts"], holdout["baseline_attempts"], "baseline phase drift")
    result["baseline_attempts"] = validation["baseline_attempts"]
    result["total_prefix_cases"] += 2742 + 9216
    result["total_constructors"] += 24
    for B in WIDTHS:
        result["widths"].append(dict(B=B, regression=paired(validation["rows"], B, "regression"),
                                     training=paired(validation["rows"], B, "fresh"),
                                     holdout=paired(holdout["rows"], B, "fresh")))
    result["outcome"] = "FAIL" if any(row["holdout"][0]["candidate_failures"] for row in result["widths"]) else "PASS"
    c.require(result["total_prefix_cases"] <= 141750, "prefix case cap")
    return result


def execute_worker(argv, raw_path, error_path, timeout, prior_bytes=0):
    """Count pipe bytes before writing: a child cannot transiently exceed caps."""
    c.require(0 <= prior_bytes <= TOTAL_LIMIT and 0 < timeout <= 60, "supervisor budget")
    started, total = time.monotonic(), prior_bytes
    with contextlib.ExitStack() as stack:
        outputs = [stack.enter_context(Path(path).open("xb")) for path in (raw_path, error_path)]
        selector = stack.enter_context(selectors.DefaultSelector())
        process = subprocess.Popen(argv, cwd=str(ROOT), env={"PATH": "/usr/bin:/bin", "LC_ALL": "C"},
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
        streams = (process.stdout, process.stderr)
        sizes = [0, 0]
        try:
            for index, stream in enumerate(streams):
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ, index)
            while selector.get_map():
                remaining = timeout - (time.monotonic() - started)
                c.require(remaining > 0, "worker deadline")
                for key, _ in selector.select(min(remaining, 0.05)):
                    data = os.read(key.fd, 65536)
                    if not data:
                        selector.unregister(key.fileobj)
                        continue
                    index = key.data
                    # Preserve bounded bytes already produced, then fail closed.
                    allowed = min(STREAM_LIMIT - 1 - sizes[index], TOTAL_LIMIT - total)
                    kept = data[:allowed]
                    outputs[index].write(kept)
                    sizes[index] += len(kept)
                    total += len(kept)
                    c.require(len(kept) == len(data), "raw byte cap")
            remaining = timeout - (time.monotonic() - started)
            c.require(remaining > 0, "worker deadline")
            code = process.wait(timeout=remaining)
            c.require(code == 0 and sizes[1] == 0, "worker failed/stderr")
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise
        finally:
            for stream in streams:
                stream.close()
    return total


def run_once(build_path, bundle_path):
    receipt = c.parse_json(Path(build_path).read_bytes())
    check_build(receipt)
    history = c.parse_json(Path(build_path).with_name("inputs.json").read_bytes())
    exact(history, load_history(), "frozen historical input drift")
    c.require(not c.command(["git", "status", "--porcelain", "--untracked-files=no", "--"] + list(SOURCE_PATHS)).strip(),
              "harness tree dirty")
    for path in SOURCE_PATHS:
        c.command(["git", "ls-files", "--error-unmatch", "--", path])
    c.require(c.command(["git", "rev-parse", "@{upstream}"]).strip() == receipt["source_commit"], "source not pushed")
    bundle = Path(bundle_path).resolve()
    c.require(bundle == NAMESPACE, "frozen one-shot namespace")
    bundle.mkdir(exist_ok=False)
    started = time.monotonic()
    freeze = frozen_protocol(receipt, history)
    c.write_json(bundle / "freeze.json", freeze)
    completed, workers, worker_elapsed, raw_bytes, selection_sha = [], [], 0.0, 0, None
    try:
        search = validation = holdout = None
        for phase in ("search", "validate", "holdout"):
            remaining = min(60 - worker_elapsed, 70 - (time.monotonic() - started))
            c.require(remaining > 0, "aggregate worker deadline")
            argv = [receipt["interpreter"], "-I", "-B", "-S", receipt["worker"], "--worker", phase,
                    "--freeze", str(bundle / "freeze.json"), "--seconds", str(remaining)]
            if selection_sha is not None:
                c.require(c.file_digest(bundle / "selection.json") == selection_sha, "sealed map changed")
                argv += ["--selection", str(bundle / "selection.json"), "--selection-sha256", selection_sha]
            worker_started = time.monotonic()
            try:
                raw_bytes = execute_worker(argv, bundle / (phase + ".jsonl"), bundle / (phase + ".stderr"), remaining, raw_bytes)
            finally:
                elapsed = time.monotonic() - worker_started
                workers.append(dict(phase=phase, elapsed_seconds=elapsed))
                worker_elapsed += elapsed
            c.require(worker_elapsed < 60, "aggregate worker deadline")
            text = (bundle / (phase + ".jsonl")).read_text(encoding="ascii")
            if phase == "search":
                search = parse_search(text, freeze)
                completed.append(phase)
                selection = dict(protocol=PROTOCOL, selected_attempts=search["selected_attempts"],
                                 search_sha256=c.file_digest(bundle / "search.jsonl"))
                c.write_json(bundle / "selection.json", selection)
                selection_sha = c.file_digest(bundle / "selection.json")
                if -1 in search["selected_attempts"]:
                    break
            else:
                parsed = parse_validation(text, phase, freeze, selection, search)
                if phase == "validate":
                    validation = parsed
                else:
                    holdout = parsed
                completed.append(phase)
        result = summarize(search, validation, holdout)
        c.require(c.file_digest(bundle / "selection.json") == selection_sha, "sealed map changed")
        check_build(receipt)
        c.require(time.monotonic() - started < 70, "outer deadline")
    except Exception as error:
        result = dict(outcome="INVALID", error=str(error), completed_phases=completed)
    result.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic() - started, workers=workers,
                  promotion_claimed=False, all_K_claimed=False, timing_claimed=False)
    return publish_result(bundle, result, started)


def publish_result(bundle, result, started):
    c.write_json(bundle / "summary.json", result)
    hashes = {path.name: c.file_digest(path) for path in sorted(bundle.iterdir()) if path.is_file()}
    completed = result.get("completed_phases", ["search"] if result["outcome"] == "EXHAUSTED"
                           else ["search", "validate", "holdout"])

    def refresh():
        nonlocal result
        elapsed = time.monotonic() - started
        if elapsed >= 70:
            result = dict(outcome="INVALID", error="outer deadline during publication", completed_phases=completed,
                          protocol=PROTOCOL, elapsed_seconds=elapsed, workers=result["workers"],
                          promotion_claimed=False, all_K_claimed=False, timing_claimed=False)
        else:
            result = dict(result, elapsed_seconds=elapsed)
        data = c.canonical(result) + b"\n"
        # Only our new, unsealed summary is replaced; old bundles are untouched.
        (bundle / "summary.json").write_bytes(data)
        hashes["summary.json"] = c.digest(data)

    refresh()  # Account for all bulk artifact hashing in the final receipt.
    if time.monotonic() - started >= 70:
        refresh()  # Final pre-seal boundary, including the final summary write.
    # The terminal small manifest write is outside this self-referential timer.
    # The coordinating process also checks complete CLI elapsed <70 seconds.
    c.write_json(bundle / "COMPLETE", hashes)
    return result


def replay(bundle_path):
    bundle = Path(bundle_path)
    complete = c.parse_json((bundle / "COMPLETE").read_bytes())
    c.require(type(complete) is dict and set(complete) == {path.name for path in bundle.iterdir() if path.name != "COMPLETE"},
              "bundle roster")
    for name, digest in complete.items():
        c.require(Path(name).name == name and c.file_digest(bundle / name) == digest, "bundle digest")
    freeze = c.parse_json((bundle / "freeze.json").read_bytes())
    validate_freeze(freeze)
    result = c.parse_json((bundle / "summary.json").read_bytes())
    common = {"protocol", "elapsed_seconds", "workers", "promotion_claimed", "all_K_claimed", "timing_claimed"}
    c.require(type(result) is dict and result.get("protocol") == PROTOCOL and
              all(result.get(key) is False for key in ("promotion_claimed", "all_K_claimed", "timing_claimed")), "summary scope")
    elapsed = result.get("elapsed_seconds")
    c.require(type(elapsed) in (int, float) and math.isfinite(elapsed) and elapsed >= 0, "elapsed")
    raw = [name for name in complete if name.endswith((".jsonl", ".stderr"))]
    c.require(all((bundle / name).stat().st_size < STREAM_LIMIT for name in raw) and
              sum((bundle / name).stat().st_size for name in raw) <= TOTAL_LIMIT, "raw byte bounds")
    invalid = result.get("outcome") == "INVALID"
    if invalid:
        c.require(set(result) == common | {"outcome", "error", "completed_phases"} and
                  type(result["error"]) is str and bool(result["error"]), "INVALID schema")
        completed = result["completed_phases"]
        c.require(type(completed) is list and
                  completed in ([], ["search"], ["search", "validate"], ["search", "validate", "holdout"]), "phase progress")
    else:
        c.require(result.get("outcome") in ("PASS", "FAIL", "EXHAUSTED") and elapsed < 70, "scientific outcome/time")
        completed = ["search"] if result["outcome"] == "EXHAUSTED" else ["search", "validate", "holdout"]
    workers = result.get("workers")
    c.require(type(workers) is list and len(completed) <= len(workers) <= min(3, len(completed) + int(invalid)),
              "worker execution roster")
    for index, worker in enumerate(workers):
        c.require(type(worker) is dict and set(worker) == {"phase", "elapsed_seconds"} and
                  worker["phase"] == ("search", "validate", "holdout")[index], "worker execution schema/order")
        duration = worker["elapsed_seconds"]
        c.require(type(duration) in (int, float) and math.isfinite(duration) and duration >= 0, "worker elapsed")
    worker_elapsed = sum(worker["elapsed_seconds"] for worker in workers)
    c.require(worker_elapsed <= elapsed and (invalid or worker_elapsed < 60), "aggregate worker elapsed")
    required = {"freeze.json", "summary.json"}
    search = validation = holdout = selection = None
    for phase in completed:
        required.update((phase + ".jsonl", phase + ".stderr"))
        c.require((bundle / (phase + ".stderr")).read_bytes() == b"", "completed worker stderr")
        text = (bundle / (phase + ".jsonl")).read_text(encoding="ascii")
        if phase == "search":
            search = parse_search(text, freeze)
            required.add("selection.json")
            selection = c.parse_json((bundle / "selection.json").read_bytes())
            validate_selection(selection, search, c.file_digest(bundle / "search.jsonl"))
        else:
            parsed = parse_validation(text, phase, freeze, selection, search)
            if phase == "validate":
                validation = parsed
            else:
                holdout = parsed
    extras = set(complete) - required
    if invalid:
        next_phase = ("search", "validate", "holdout")[len(completed)] if len(completed) < 3 else None
        allowed = {next_phase + suffix for suffix in (".jsonl", ".stderr")} if next_phase else set()
        c.require(required <= set(complete) and extras <= allowed, "INVALID partial phase roster")
        if search is not None and -1 in search["selected_attempts"]:
            c.require(not extras and validation is None, "work after exhausted search")
        return result
    c.require(set(complete) == required, "scientific file roster")
    expected = summarize(search, validation, holdout)
    c.require(set(result) == set(expected) | common, "summary schema")
    for key, value in expected.items():
        exact(result[key], value, "summary " + key)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build")
    modes.add_argument("--run")
    modes.add_argument("--replay")
    modes.add_argument("--worker", choices=("search", "validate", "holdout"))
    modes.add_argument("--selftest", action="store_true")
    parser.add_argument("--build-receipt")
    parser.add_argument("--freeze")
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--selection")
    parser.add_argument("--selection-sha256")
    args = parser.parse_args()
    if args.worker:
        c.require(args.freeze and args.seconds is not None, "worker arguments")
        worker(args.worker, args.freeze, args.seconds, args.selection, args.selection_sha256)
    elif args.selftest:
        # This inherited selftest uses only the neutral root and each actual
        # baseline's own attempt, never a joint candidate or frozen fresh root.
        p.selftest()
    elif args.build:
        p.emit(build(args.build))
    elif args.run:
        c.require(args.build_receipt, "run needs build receipt")
        p.emit(run_once(args.build_receipt, args.run))
    else:
        p.emit(replay(args.replay))


if __name__ == "__main__":
    try:
        main()
    except (c.ValidationError, OSError, subprocess.SubprocessError) as error:
        print("production MIX3 joint r0: " + str(error), file=sys.stderr)
        sys.exit(2)
