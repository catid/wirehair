#!/usr/bin/env python3
"""Six-cell width-local public MIX3 feasibility screen, never certification.

Build requires the four immutable, machine-local historical bundles (.52/.54/
.55/.56) and the pinned original production DSO. The resulting inputs embed a
hash-pinned .56 receipt and complete origin ledgers; replay needs none of those
old paths and never invokes the codec. Rejected first-failure observations are
not a census of unvisited widths or independent loss trials.
"""
import contextlib
import hashlib
import importlib.util
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
WORKER_SOURCE = Path(__file__).resolve()
PROTOCOL = "wirehair.wh2.production-mix3-k3k6-width-local-r0"
NAMESPACE = Path("/var/tmp/wh2-production-mix3-k3k6-width-local-r0")
spec = importlib.util.spec_from_file_location("_wh2_width_local_joint_private",
                                             ROOT / "bench/Wh2ProductionMix3RecoveryJointR0.py")
j = importlib.util.module_from_spec(spec)
spec.loader.exec_module(j)
p, b, c = j.p, j.b, j.c
CELLS = ((3, 2), (3, 64), (3, 1280), (6, 2), (6, 64), (6, 1280))
ROSTER_SHA256 = "2297d69d6b7a5a83f715589f981deea0d48a676599e02e5697d8ecff6dcfa89f"
LEDGER_SHA256 = "6f1fac82e97e72c124b83c495316eca69acb4180097f8ea60e05674404fca33b"
ORIGINS_SHA256 = "89e514a5b6a673c0fe5a5188503c0f125717eb35962f48cd7bfa55a3646d8d39"
BAD_ORIGINS_SHA256 = "2cca4560caa79b35f7d5103a694231e22d2797a4f466109fb1b3f9dce3241a54"
KNOWN_BAD_SHA256 = "1438930c6732e27b8f9a34ce7a8b788f441d9a4f25400fb109ce9b38cfb218dd"
PRIOR_ROOTS_SHA256 = "5b89e3a62123156bd8a4918f34d1364e52dd9e0d091300a302d944c344c3997d"
CONSUMED_ROOTS_SHA256 = "ada764d2cb8afa51e9ea13c130e8536ae0b2d23a4ad1209b532987739c1b291e"
UNCONSUMED_ROOTS_SHA256 = "85e942184cb461bd91914c9f97e293d13daf8564413af97a9440ff86af641a65"
JOINT_BUNDLE = Path("/var/tmp/wh2-production-mix3-k3k6-joint-r0")
JOINT_COMPLETE_SHA256 = "3b939cc601edf2b75e0ce62a0d06744b4b9213cb738f0e2f23e5164fe317a524"
JOINT_SEARCH_SHA256 = "88f3f8515167e869d5c4cfcbdb07172c7f75172dbf5e4f16ade4844aed15a2cf"
TRAINING_ROOTS, HOLDOUT_ROOTS = (tuple("0x" + c.digest((PROTOCOL + ":" + phase + "/" + str(i)).encode("ascii"))[:16]
                                    for i in range(count)) for phase, count in (("train", 16), ("holdout", 64)))
SOURCE_PATHS = tuple(sorted(set(j.SOURCE_PATHS) - {"bench/test_Wh2ProductionMix3RecoveryJointR0.py"} |
                            {"bench/Wh2ProductionMix3RecoveryWidthLocalR0.py",
                             "bench/test_Wh2ProductionMix3RecoveryWidthLocalR0.py"}))
IMMUTABLE_HELPERS = dict(j.IMMUTABLE_HELPERS, **{
    "bench/Wh2ProductionMix3RecoveryJointR0.py": "c67769b4a0bf0c4511d183704378e6a3486b0cc9f210a67a189fae5498371858"})
_old_load_history = j.load_history


def history_from_basis(basis):
    c.require(type(basis) is dict and set(basis) == {"complete", "freeze", "search", "stderr", "selection", "summary"},
              "historical basis schema")
    c.require(c.digest(c.canonical(basis["complete"]) + b"\n") == JOINT_COMPLETE_SHA256, "joint manifest drift")
    files = {"freeze.json": c.canonical(basis["freeze"]) + b"\n",
             "search.jsonl": b"".join(c.canonical(row) + b"\n" for row in basis["search"]),
             "search.stderr": basis["stderr"].encode("ascii"),
             "selection.json": c.canonical(basis["selection"]) + b"\n",
             "summary.json": c.canonical(basis["summary"]) + b"\n"}
    j.exact({name: c.digest(value) for name, value in files.items()}, basis["complete"], "joint artifact identities")
    c.require(c.digest(files["search.jsonl"]) == JOINT_SEARCH_SHA256, "joint raw identity")
    previous = basis["freeze"]["history"]
    j.validate_history(previous)  # Immutable .52/.54/.55 ledgers and old roots.
    origins, bad_origins = list(previous["origins"]), list(previous["bad_origins"])
    for line, row in enumerate(basis["search"], 1):
        if row["type"] != "attempt":
            continue
        base = dict(issue=56, file="search.jsonl", line=line, K=row["K"], arm="candidate", attempt=row["attempt"])
        witness = row["witness"]
        if witness is not None:
            c.require(witness["outcomes"] == [1], "historical recovery witness")
            origins.append(dict(base, B=witness["B"], root=witness["root"], schedule=witness["schedule"],
                                overhead=witness["overheads"][0], ids=witness["ids"]))
        else:
            con = row["constructors"][-1]
            c.require(con["outcome"] == 10, "historical construction witness")
            bad_origins.append(dict(base, B=con["B"], root=None, schedule=None))
    prefixes = {}
    for row in origins:
        prefixes.setdefault((row["K"], tuple(row["ids"])), set()).add(row["B"])
    ledger = [dict(K=K, ids=list(ids), original_widths=sorted(widths), overhead=len(ids) - K)
              for (K, ids), widths in sorted(prefixes.items())]
    consumed = sorted(set(previous["prior_roots"]) | set(basis["freeze"]["training_roots"]))
    unconsumed = sorted(basis["freeze"]["holdout_roots"])
    known_bad = [list(profile) for profile in sorted({(row["K"], row["B"], row["attempt"]) for row in bad_origins})]
    history = dict(basis=basis, prefixes=ledger, origins=origins, bad_origins=bad_origins, known_bad=known_bad,
                   consumed_roots=consumed, unconsumed_holdout_roots=unconsumed,
                   prior_roots=sorted(set(consumed) | set(unconsumed)))
    for key, count, digest in (("prefixes", 75, LEDGER_SHA256), ("origins", 212, ORIGINS_SHA256),
                               ("bad_origins", 786, BAD_ORIGINS_SHA256), ("known_bad", 426, KNOWN_BAD_SHA256),
                               ("prior_roots", 207, PRIOR_ROOTS_SHA256), ("consumed_roots", 143, CONSUMED_ROOTS_SHA256),
                               ("unconsumed_holdout_roots", 64, UNCONSUMED_ROOTS_SHA256)):
        c.require(len(history[key]) == count and c.digest(c.canonical(history[key])) == digest, key + " history drift")
    roots = TRAINING_ROOTS + HOLDOUT_ROOTS
    c.require(c.digest(c.canonical({"train": TRAINING_ROOTS, "holdout": HOLDOUT_ROOTS})) == ROSTER_SHA256 and
              len(set(roots)) == 80 and not set(roots) & (set(history["prior_roots"]) | set(j.RESERVED)),
              "fresh root disjointness/hash")
    return history


def validate_history(history):
    c.require(type(history) is dict and "basis" in history, "history schema")
    j.exact(history, history_from_basis(history["basis"]), "complete no-omission history")


def load_history():
    previous = _old_load_history()  # Requires and verifies the three earlier bundles.
    c.require(c.file_digest(JOINT_BUNDLE / "COMPLETE") == JOINT_COMPLETE_SHA256, "joint manifest drift")
    basis = {key: c.parse_json((JOINT_BUNDLE / filename).read_bytes()) for key, filename in
             (("complete", "COMPLETE"), ("freeze", "freeze.json"), ("selection", "selection.json"), ("summary", "summary.json"))}
    for name, digest in basis["complete"].items():
        c.require(Path(name).name == name and c.file_digest(JOINT_BUNDLE / name) == digest, "joint artifact drift")
    basis["search"] = [c.parse_json(line) for line in (JOINT_BUNDLE / "search.jsonl").read_bytes().splitlines()]
    basis["stderr"] = (JOINT_BUNDLE / "search.stderr").read_text(encoding="ascii")
    j.exact(previous, basis["freeze"]["history"], "embedded earlier history disagrees with original bundles")
    return history_from_basis(basis)


def frozen_protocol(receipt, history):
    j.validate_receipt(receipt)
    validate_history(history)
    return dict(protocol=PROTOCOL, build=receipt, history=history, cells=[list(cell) for cell in CELLS],
                training_roots=list(TRAINING_ROOTS), holdout_roots=list(HOLDOUT_ROOTS), roster_sha256=ROSTER_SHA256,
                schedules=list(c.SCHEDULES), loss_ppm=500000, overheads=list(c.OVERHEADS),
                worker_budget_seconds=60, outer_budget_seconds=70, stream_limit_exclusive=j.STREAM_LIMIT,
                aggregate_limit_inclusive=j.TOTAL_LIMIT, max_prefix_cases=143298, max_constructors=1560)


@contextlib.contextmanager
def at_width(B, history):
    # These are private module globals, restored even after an oracle exception.
    original = j.WIDTHS, j.TRAINING_ROOTS, j.HOLDOUT_ROOTS, j.BAD_PROFILES
    j.WIDTHS, j.TRAINING_ROOTS, j.HOLDOUT_ROOTS, j.BAD_PROFILES = (B,), TRAINING_ROOTS, HOLDOUT_ROOTS, history["known_bad"]
    try:
        yield
    finally:
        j.WIDTHS, j.TRAINING_ROOTS, j.HOLDOUT_ROOTS, j.BAD_PROFILES = original


def specs_for(K, B, phase, history):
    with at_width(B, history):
        return j.specs_for(K, phase, history)


def validate_selection(selection, search=None, search_sha256=None):
    c.require(type(selection) is dict and set(selection) == {"protocol", "selected_attempts", "search_sha256"} and
              selection["protocol"] == PROTOCOL, "selection schema")
    selected = selection["selected_attempts"]
    c.require(type(selected) is list and len(selected) == 6 and
              all(type(x) is int and -1 <= x <= 255 for x in selected), "six-cell selection")
    j.hex_value(selection["search_sha256"], 64, "search digest")
    if search is not None:
        j.exact(selected, search["selected_attempts"], "sealed search map")
        c.require(selection["search_sha256"] == search_sha256, "sealed search hash")
    return selected


def worker(phase, freeze_path, seconds, selection_path=None, selection_sha256=None):
    c.require(phase in ("search", "validate", "holdout") and type(seconds) in (int, float) and
              math.isfinite(seconds) and 0 < seconds <= 60, "worker arguments")
    deadline = j.time.monotonic() + seconds
    freeze = c.parse_json(Path(freeze_path).read_bytes())
    j.validate_freeze(freeze)
    j.check_build(freeze["build"])
    selected = [-1] * 6
    if phase != "search":
        c.require(selection_path is not None and c.file_digest(selection_path) == selection_sha256, "sealed map hash")
        selected = validate_selection(c.parse_json(Path(selection_path).read_bytes()))
        c.require(-1 not in selected, "work after exhausted cell")
    else:
        c.require(selection_path is None and selection_sha256 is None, "search map supplied")
    native, rows = p.Native(), 0
    p.emit(j.begin_record(phase, freeze, selection_sha256))
    for cell_index, (K, B) in enumerate(CELLS):
        with at_width(B, freeze["history"]):
            specs = j.specs_for(K, phase, freeze["history"])
            if phase == "search":
                for attempt in range(256):
                    row = dict(j.search_attempt(native, K, attempt, specs, deadline), B=B)
                    p.emit(row)
                    rows += 1
                    if row["accepted"]:
                        selected[cell_index] = attempt
                        break
            else:
                for arm, attempt in (("baseline", None), ("candidate", selected[cell_index])):
                    j.check_time(deadline)
                    with native.encoder(K, B, attempt) as state:
                        for case in specs:
                            p.emit(j.make_cell(native, state, K, B, arm, case, case["overheads"], deadline))
                            rows += 1
    j.check_time(deadline)
    p.emit(dict(type="terminal", phase=phase, rows=rows, selected_attempts=selected))


def records(text, phase, freeze, selection_sha256):
    c.require(type(text) is str and text.endswith("\n") and len(text.encode("ascii")) < j.STREAM_LIMIT,
              "stream length/termination")
    rows = [c.parse_json(line) for line in text.splitlines()]
    c.require(len(rows) >= 2, "stream records")
    j.exact(rows[0], j.begin_record(phase, freeze, selection_sha256), "begin")
    terminal = rows[-1]
    c.require(type(terminal) is dict and set(terminal) == {"type", "phase", "rows", "selected_attempts"}, "terminal schema")
    j.exact({key: terminal[key] for key in ("type", "phase", "rows")},
            dict(type="terminal", phase=phase, rows=len(rows) - 2), "terminal count")
    selected = validate_selection(dict(protocol=PROTOCOL, selected_attempts=terminal["selected_attempts"], search_sha256="0" * 64))
    return rows[1:-1], selected


def parse_search(text, freeze):
    rows, selected = records(text, "search", freeze, None)
    position, choices, prefix_cases = 0, [], 0
    for K, B in CELLS:
        with at_width(B, freeze["history"]):
            specs, choice = j.specs_for(K, "search", freeze["history"]), -1
            for attempt in range(256):
                c.require(position < len(rows), "missing attempt")
                row = rows[position]
                position += 1
                c.require(type(row) is dict and set(row) == {"type", "K", "B", "attempt", "constructors", "passed",
                                                            "pass_sha256", "witness", "accepted"}, "attempt schema")
                j.exact([row["type"], row["K"], row["B"], row["attempt"]], ["attempt", K, B, attempt], "attempt order")
                cons = row["constructors"]
                c.require(type(cons) is list and len(cons) == 1 and type(cons[0]) is dict and
                          set(cons[0]) == {"B", "profile_hex", "outcome"}, "one local constructor")
                con = cons[0]
                j.exact([con["B"], con["profile_hex"]], [B, p.profile_bytes(K, B, attempt).hex()], "constructor profile")
                c.require(type(con["outcome"]) is int and con["outcome"] in (0, 10), "constructor outcome")
                if [K, B, attempt] in j.BAD_PROFILES:
                    c.require(con["outcome"] == 10, "known BadSeed constructor contradicted")
                passed = c.integer(row["passed"], 0, len(specs), "passed constraints")
                j.hex_value(row["pass_sha256"], 64, "passed hash")
                c.require(passed != 0 or row["pass_sha256"] == c.digest(b""), "empty passed hash")
                c.require(type(row["accepted"]) is bool, "accepted type")
                if con["outcome"] == 10:
                    c.require(passed == 0 and row["witness"] is None and not row["accepted"], "construction witness")
                elif passed == len(specs):
                    c.require(row["witness"] is None and row["accepted"], "accepted attempt")
                    choice = attempt
                else:
                    j.validate_cell(row["witness"], K, "candidate", attempt, specs[passed], search=True)
                    c.require(row["witness"]["outcomes"] == [1] and not row["accepted"], "first failure witness")
                    prefix_cases += 1
                prefix_cases += passed
                if choice != -1:
                    break
            choices.append(choice)
    c.require(position == len(rows), "extra attempts/late selection")
    j.exact(selected, choices, "six-cell first-success map")
    c.require(prefix_cases <= 131328 and len(rows) <= 1536, "search bounds")
    return dict(selected_attempts=choices, records=rows, prefix_cases=prefix_cases, constructors=len(rows))


def parse_validation(text, phase, freeze, selection, search):
    c.require(phase in ("validate", "holdout"), "validation phase")
    rows, selected = records(text, phase, freeze, c.digest(c.canonical(selection) + b"\n"))
    j.exact(selected, selection["selected_attempts"], "validation map")
    c.require(-1 not in selected, "validation after exhaustion")
    position, baselines = 0, []
    for cell_index, (K, B) in enumerate(CELLS):
        with at_width(B, freeze["history"]):
            specs, semantic = j.specs_for(K, phase, freeze["history"]), {}
            for arm in ("baseline", "candidate"):
                digest, actual, bad = hashlib.sha256(), None, None
                for case in specs:
                    c.require(position < len(rows) and type(rows[position]) is dict, "missing validation cell")
                    row = rows[position]
                    position += 1
                    attempt = selected[cell_index] if arm == "candidate" else c.integer(row.get("attempt"), 0, 255, "baseline attempt")
                    j.validate_cell(row, K, arm, attempt, case)
                    current_bad = row["outcomes"][0] == 10
                    if actual is None:
                        actual, bad = attempt, current_bad
                    c.require(actual == attempt and bad == current_bad, "constructor drift")
                    key, value = (attempt, case["kind"], case["index"]), dict(row, arm="candidate")
                    if key in semantic:
                        j.exact(value, semantic[key], "same attempt semantic mismatch")
                    semantic[key] = value
                    if arm == "candidate" and phase == "validate":
                        c.require(all(x == 0 for x in row["outcomes"]), "selected revalidation failure")
                        digest.update(j.success_token(row))
                if arm == "baseline":
                    baselines.append(actual)
                elif phase == "validate":
                    accepted = next(row for row in search["records"] if row["K"] == K and row["B"] == B and row["accepted"])
                    c.require(digest.hexdigest() == accepted["pass_sha256"], "selected search/revalidation byte mismatch")
    c.require(position == len(rows), "extra validation cells")
    return dict(rows=rows, baseline_attempts=baselines)


def summarize(search, validation=None, holdout=None):
    selected = search["selected_attempts"]
    result = dict(selected_attempts=selected, search_prefix_cases=search["prefix_cases"], search_constructors=search["constructors"],
                  total_prefix_cases=search["prefix_cases"], total_constructors=search["constructors"], cells=[])
    if -1 in selected:
        c.require(validation is None and holdout is None, "validation after exhaustion")
        result["outcome"] = "EXHAUSTED"
        return result
    c.require(validation is not None and holdout is not None, "missing validation/holdout")
    j.exact(validation["baseline_attempts"], holdout["baseline_attempts"], "baseline phase drift")
    result["baseline_attempts"] = validation["baseline_attempts"]
    result["total_prefix_cases"] += 2754 + 9216
    result["total_constructors"] += 24
    for K, B in CELLS:
        training = [row for row in validation["rows"] if row["K"] == K]
        held = [row for row in holdout["rows"] if row["K"] == K]
        result["cells"].append(dict(K=K, B=B, regression=j.paired(training, B, "regression"),
                                    training=j.paired(training, B, "fresh"), holdout=j.paired(held, B, "fresh")))
    result["outcome"] = "FAIL" if any(cell["holdout"][0]["candidate_failures"] for cell in result["cells"]) else "PASS"
    c.require(result["total_prefix_cases"] <= 143298 and result["total_constructors"] <= 1560, "case/constructor cap")
    return result


# Dispatch only inside this private imported module. The immutable controller
# already implements source/receipt pins, live pipe caps, separate worker/outer
# budgets, sealed-map phases, no-rerun namespace, and strict offline replay.
j.PROTOCOL = c.PROTOCOL = PROTOCOL
j.NAMESPACE, j.WORKER_SOURCE, j.SOURCE_PATHS = NAMESPACE, WORKER_SOURCE, SOURCE_PATHS
j.IMMUTABLE_HELPERS = IMMUTABLE_HELPERS
j.load_history, j.frozen_protocol = load_history, frozen_protocol
j.worker, j.parse_search, j.parse_validation = worker, parse_search, parse_validation
j.validate_selection, j.summarize = validate_selection, summarize
j.__doc__ = __doc__
build, check_build, validate_receipt = j.build, j.check_build, j.validate_receipt
run_once, replay, validate_freeze = j.run_once, j.replay, j.validate_freeze


if __name__ == "__main__":
    try:
        j.main()
    except (c.ValidationError, OSError, subprocess.SubprocessError) as error:
        print("production MIX3 width-local r0: " + str(error), file=sys.stderr)
        sys.exit(2)
