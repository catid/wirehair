#!/usr/bin/env python3
"""Inert, bounded authentication of existing .58/.62 recovery evidence.

No historical campaign module, coefficient mapper, trace generator, or codec
is imported. Loading copies recorded bytes; it never scores the fixed pair.
"""
import importlib.util
import io
import os
from pathlib import Path


_SPEC = importlib.util.spec_from_file_location(
    "_thue_recovery_historical_reader", Path(__file__).resolve().with_name(
        "Wh2NoncommutingRadixRunR0.py"))
C = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(C)

WIDTH_LOCAL = Path("/var/tmp/wh2-production-mix3-k3k6-width-local-r0")
THUE_MORSE = Path("/var/tmp/wh2-thue-morse-r0")
WIDTH_MANIFEST = "135302425b619f70c62d0fd4a6f0b3eec642314fdf19741c77986336be4b78b7"
THUE_MANIFEST = "ef1ae96b32686560548ce75a039bd370dd61d3d4eb5ae3ffed345036c248f675"
THUE_RAW = "0af037929a52e31263d4b5681cba761a91880ac2e0f3b13dfd9a89de60a01975"
PAIR_SHA256 = "d28da7ebd5ab6b9589589ffc8a4f146fe1ce3a90ed3b37fcff77a8fc2efd9bfb"
PREFIX_SHA256 = "e9a7e15d5eec740255b467dc8b39dfc5b31f2d0b08441b2a0b3bb930c20b5c86"
ORIGIN_SHA256 = "653df33ce21259d6cb4c2d0721598b7f9484f6927f66593260858b8bdc66e77b"
INPUT_CAP = 32 * 1024 * 1024
AGGREGATE_CAP = 64 * 1024 * 1024
WIDTH_FILES = ("freeze.json", "holdout.jsonl", "holdout.stderr", "search.jsonl",
               "search.stderr", "selection.json", "summary.json", "validate.jsonl", "validate.stderr")
THUE_FILES = ("CLAIM.json", "raw.json", "stderr.txt", "summary.json")
# Exact main-contract root literals; no retired profile-selection code import.
# Projection: wh2_benchmark_contract_v4.json (wh2-pure-gf256-1pct-v7).
# Extra exclusion-only literals: wh2_mix2_seed_repair_contract.py:132-147.
HARD_TRAINING_ROOTS = ("0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3")
HARD_VALIDATION_ROOTS = ("0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b")
MAIN_ROOTS = HARD_TRAINING_ROOTS + HARD_VALIDATION_ROOTS + (
    "0xb889883a79549774", "0xb5666de0987896af", "0x8bfca269b0bc01e0",
    "0xefd20c982041a46b", "0x8827bc36ed906555", "0x86029f23d6132efa")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, lower, upper):
    return type(value) is int and lower <= value <= upper


def root_value(value):
    return (type(value) is str and len(value) == 18 and value.startswith("0x") and
            all(c in "0123456789abcdef" for c in value[2:]))


def packet_ids(values):
    require(type(values) is list and all(integer(x, 0, (1 << 32) - 1) for x in values)
            and len(values) == len(set(values)), "historical packet IDs")
    return values


def read_bundle(directory, manifest_name, pin, filenames, structured, budget, deadline=None):
    """Authenticate every exact member, retaining no paths from untrusted JSON."""
    def read(name):
        remaining = AGGREGATE_CAP - budget[0]
        require(remaining >= 0, "historical aggregate byte cap")
        raw = C.read_regular(directory / name, min(INPUT_CAP, remaining), deadline=deadline)
        budget[0] += len(raw)
        require(budget[0] <= AGGREGATE_CAP, "historical aggregate byte cap")
        return raw

    raw = read(manifest_name)
    require(C.sha(raw) == pin, "historical manifest pin")
    manifest = C.strict_json(raw)
    require(raw == C.canonical(manifest) + (b"" if structured else b"\n"),
            "historical canonical manifest")
    if structured:
        require(type(manifest) is dict and set(manifest) == {"protocol", "outcome", "files"}
                and manifest["protocol"] == "wirehair.wh2.thue-morse-r0"
                and manifest["outcome"] == "PASS", "Thue-Morse manifest verdict")
        members = manifest["files"]
    else:
        members = manifest
    require(type(members) is dict and set(members) == set(filenames), "historical member roster")
    expected_names, observed_names = set(filenames) | {manifest_name}, set()
    with os.scandir(str(directory)) as paths:
        for path in paths:
            C.time_left(deadline)
            require(path.name in expected_names and path.name not in observed_names,
                    "historical directory roster")
            observed_names.add(path.name)
    require(observed_names == expected_names, "historical directory roster")
    content, entries = {}, {}
    for name in sorted(filenames):
        data = read(name)
        actual = dict(bytes=len(data), sha256=C.sha(data))
        expected = actual if structured else actual["sha256"]
        require(C.canonical(members[name]) == C.canonical(expected), "historical member pin: " + name)
        content[name], entries[name] = data, actual
    return content, dict(path=str(directory), manifest_sha256=pin,
                         manifest_bytes=len(raw), files=entries)


def validate_origin(origin):
    keys = {"issue", "file", "line", "K", "B", "arm", "attempt", "root", "schedule", "overhead", "ids"}
    require(type(origin) is dict and set(origin) == keys, "historical origin schema")
    require(integer(origin["issue"], 0, 100) and origin["issue"] in (52, 54, 55, 56, 58)
            and type(origin["file"]) is str and Path(origin["file"]).name == origin["file"]
            and origin["file"].endswith(".jsonl") and integer(origin["line"], 1, 1000000)
            and type(origin["K"]) is int and origin["K"] == 6
            and type(origin["B"]) is int and origin["B"] in (2, 64, 1280)
            and origin["arm"] in ("baseline", "candidate") and integer(origin["attempt"], 0, 255)
            and integer(origin["overhead"], 0, 4), "historical origin coordinates")
    require((origin["root"] is None and origin["schedule"] is None) or
            (root_value(origin["root"]) and origin["schedule"] in ("burst", "adversarial", "repair-only")),
            "historical origin trace provenance")
    require(len(packet_ids(origin["ids"])) == 6 + origin["overhead"], "historical prefix length")


def extract_origins(previous, streams, deadline=None):
    """Prior order, then lexicographic filenames, raw line order, overhead order."""
    require(type(previous) is list, "prior origins schema")
    origins = []
    for origin in previous:
        require(type(origin) is dict and type(origin.get("K")) is int, "prior origin K")
        if origin["K"] == 6:
            validate_origin(origin)
            origins.append(origin)
    for filename in sorted(streams):
        for line, raw in enumerate(io.BytesIO(streams[filename]), 1):
            C.time_left(deadline)
            require(len(raw) <= 1024 * 1024, "historical JSONL line cap")
            row = C.strict_json(raw)
            require(type(row) is dict and type(row.get("type")) is str, "historical JSONL record")
            if row["type"] == "attempt":
                require(type(row.get("accepted")) is bool, "historical attempt acceptance")
                witness = row.get("witness")
                if witness is None:
                    continue  # Accepted or typed BadSeed; neither is a recovery origin.
                require(not row["accepted"] and type(witness) is dict and
                        witness.get("outcomes") == [1] and
                        C.canonical([witness.get("K"), witness.get("B"), witness.get("attempt")]) ==
                        C.canonical([row.get("K"), row.get("B"), row.get("attempt")]),
                        "historical first-failure witness")
                row = witness
            elif row["type"] != "cell":
                continue
            require(type(row.get("K")) is int, "historical cell K")
            if row["K"] != 6:
                continue
            overheads, outcomes = row.get("overheads"), row.get("outcomes")
            require(type(overheads) is list and type(outcomes) is list and len(overheads) == len(outcomes)
                    and all(integer(x, 0, 4) for x in overheads)
                    and all(integer(x, 0, 10) for x in outcomes), "historical outcomes schema")
            ids = packet_ids(row.get("ids"))
            for overhead, outcome in zip(overheads, outcomes):
                if outcome == 1:
                    origin = dict(issue=58, file=filename, line=line, K=6, B=row.get("B"),
                                  arm=row.get("arm"), attempt=row.get("attempt"), root=row.get("root"),
                                  schedule=row.get("schedule"), overhead=overhead, ids=ids[:6 + overhead])
                    validate_origin(origin)
                    origins.append(origin)
    return origins


def prefix_ledger(origins):
    prefixes = {}
    for origin in origins:
        validate_origin(origin)
        prefixes.setdefault(tuple(origin["ids"]), set()).add(origin["B"])
    return [dict(K=6, ids=list(ids), original_widths=sorted(widths), overhead=len(ids) - 6)
            for ids, widths in sorted(prefixes.items())]


def extract_fixtures(final):
    require(type(final) is list, "Thue-Morse final evidence")
    deficient, repaired, seen = [], [], set()
    for row in final:
        require(type(row) is dict and type(row.get("family")) is str and
                integer(row.get("index"), 0, 100000), "Thue-Morse final coordinates")
        ids = packet_ids(row.get("ids"))
        ranks = row.get("prefix_ranks")
        require(type(ranks) is dict and all(type(k) is str and k in ("6", "7", "8", "9", "10", "17")
                and integer(v, 0, 6) for k, v in ranks.items()), "Thue-Morse recorded ranks")
        failures = [int(k) for k in sorted(ranks, key=int) if ranks[k] < 6]
        if not failures:
            continue
        require(len(ids) >= 8 and ranks.get("8") == 6 and all(k in (6, 7) for k in failures)
                and all(ranks[str(k)] == 5 for k in failures), "Thue-Morse expected deficient pattern")
        for length in failures:
            deficient.append(dict(family=row["family"], index=row["index"], ids=ids[:length], expected_rank=5))
        key = tuple(ids[:8])
        if key not in seen:
            seen.add(key)
            repaired.append(dict(family=row["family"], index=row["index"], ids=list(key), expected_rank=6))
    return deficient + repaired


def load_inputs(deadline=None):
    budget = [0]
    width, width_provenance = read_bundle(WIDTH_LOCAL, "COMPLETE", WIDTH_MANIFEST,
                                         WIDTH_FILES, False, budget, deadline)
    freeze = C.strict_json(width["freeze.json"])
    require(freeze["protocol"] == "wirehair.wh2.production-mix3-k3k6-width-local-r0", "historical freeze protocol")
    history = freeze["history"]
    require(len(history["origins"]) == 212 and C.sha(C.canonical(history["origins"])) ==
            "89e514a5b6a673c0fe5a5188503c0f125717eb35962f48cd7bfa55a3646d8d39", "embedded prior origins")
    origins = extract_origins(history["origins"], {name: width[name] for name in
                              ("search.jsonl", "validate.jsonl", "holdout.jsonl")}, deadline)
    prefixes = prefix_ledger(origins)
    require(len(origins) == 187 and C.sha(C.canonical(origins)) == ORIGIN_SHA256, "complete origin ledger")
    require(len(prefixes) == 82 and C.sha(C.canonical(prefixes)) == PREFIX_SHA256, "complete prefix ledger")
    roots = history["prior_roots"] + freeze["training_roots"] + freeze["holdout_roots"] + list(MAIN_ROOTS)
    require(all(root_value(root) for root in roots), "excluded root format")
    excluded_roots = sorted(set(roots))
    require(len(excluded_roots) == 297, "excluded root accounting")
    del width, freeze, history
    thue, thue_provenance = read_bundle(THUE_MORSE, "COMPLETE.json", THUE_MANIFEST,
                                       THUE_FILES, True, budget, deadline)
    require(C.sha(thue["raw.json"]) == THUE_RAW, "Thue-Morse raw pin")
    report = C.strict_json(thue["raw.json"])
    require(report["protocol"] == "wirehair.wh2.thue-morse-r0" and report["outcome"] == "PASS",
            "Thue-Morse recorded verdict")
    selection = report["selection"]
    require(C.sha(C.canonical(selection)) == report["selection_sha256"], "Thue-Morse selection seal")
    matrices = selection["A0"], selection["A1"]
    require(all(type(matrix) is list and len(matrix) == 6 and all(type(row) is list and len(row) == 6
                and all(integer(x, 0, 255) for x in row) for row in matrix) for matrix in matrices),
            "Thue-Morse recorded matrix shape")
    pair_hash = C.sha(bytes(x for matrix in matrices for row in matrix for x in row))
    require(pair_hash == selection["pair_sha256"] == PAIR_SHA256, "fixed pair pin")
    fixtures = extract_fixtures(report["final"])
    require(len(fixtures) == 33 and sum(row["expected_rank"] == 5 for row in fixtures) == 17,
            "Thue-Morse fixture accounting")
    C.time_left(deadline)
    return dict(provenance=dict(width_local=width_provenance, thue_morse=thue_provenance),
                prefixes=prefixes, origins=origins, excluded_roots=excluded_roots,
                pair=dict(A0=matrices[0], A1=matrices[1], pair_sha256=pair_hash,
                          selection_sha256=report["selection_sha256"]), fixtures=fixtures)


def input_receipt(deadline=None):
    inputs = load_inputs(deadline)
    return dict(provenance=inputs["provenance"], inputs_sha256=C.sha(C.canonical(inputs)),
                counts={key: len(inputs[key]) for key in ("prefixes", "origins", "excluded_roots", "fixtures")},
                hashes={key: C.sha(C.canonical(inputs[key])) for key in
                        ("prefixes", "origins", "excluded_roots", "fixtures")},
                pair_sha256=inputs["pair"]["pair_sha256"],
                selection_sha256=inputs["pair"]["selection_sha256"])
