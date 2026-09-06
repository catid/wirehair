#!/usr/bin/env python3
"""Frozen .63 K6 coefficient-rank recovery screen, never a codec benchmark.

Imports do not derive fresh roots, load historical artifacts, or score the
sealed candidate. Only run_screen()/--worker binds the actual campaign.
The helpers accept explicit neutral inputs for pre-campaign unit tests.
"""

import hashlib
import importlib.util
from pathlib import Path
import resource
import struct
import sys


def _sibling(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, str(Path(__file__).resolve().with_name(filename)))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


TM = _sibling("_wh2_thue_recovery_math", "Wh2ThueMorseR0.py")
RADIX = TM.RADIX
require = RADIX.require
canonical = RADIX.canonical
PROTOCOL = "wirehair.wh2.thue-morse-recovery-r0"
PAIR_SHA256 = "d28da7ebd5ab6b9589589ffc8a4f146fe1ce3a90ed3b37fcff77a8fc2efd9bfb"
WIDTHS = (2, 64, 1280)
SCHEDULES = ("iid", "burst", "adversarial", "repair-only")
OVERHEADS = (0, 1, 2, 4)
HARD_ROOTS = (
    ("training", ("0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3")),
    ("validation", ("0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b")),
)
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
CANDIDATE_LIMIT = 68096


def _integer(value, lower, upper, label):
    require(type(value) is int and lower <= value <= upper, label)
    return value


def _root(value):
    require(type(value) is str and len(value) == 18 and value.startswith("0x") and
            all(x in "0123456789abcdef" for x in value[2:]), "loss root encoding")
    return int(value, 16)


def _check(budget):
    if budget is not None:
        budget.check()


def load_inputs(deadline=None):
    """Read/authenticate history only on explicit request, never at import."""
    history = _sibling("_wh2_thue_recovery_history", "Wh2ThueMorseRecoveryHistoryR0.py")
    return history.load_inputs(deadline=deadline)


def fresh_roots():
    """CAMPAIGN ONLY: derive the frozen 512 roots, without retry or selection."""
    stem = PROTOCOL + ":fresh/"
    return ["0x" + hashlib.sha256((stem + str(i)).encode("ascii")).hexdigest()[:16]
            for i in range(512)]


def fresh_domain(root_count):
    _integer(root_count, 1, 512, "fresh root count")
    for root_index in range(root_count):
        for width in WIDTHS:
            for schedule in SCHEDULES:
                yield root_index, width, schedule


def hard_domain():
    # These six roots are the v4/v7 benchmark contract projection, not the
    # later MIX2 selection/validation partition used only for exclusions.
    for phase, roots in HARD_ROOTS:
        for root_index, root in enumerate(roots):
            for schedule in SCHEDULES[1:]:
                yield phase, root_index, 2, schedule, root


def trace(B, root_hex, schedule, budget=None):
    """Exact K6 Wh2FrozenTrace law, with only its full-ten candidate cap.

    In particular, a burst's seven continuation drops consume no RNG values.
    OH0/1/2/4 do not reseed or generate separate traces; callers slice this one.
    """
    require(type(B) is int and B in WIDTHS, "trace width")
    require(type(schedule) is str and schedule in SCHEDULES, "trace schedule")
    state = (_root(root_hex) ^ (6 * 0x9e3779b97f4a7c15) ^
             (B * 0xbf58476d1ce4e5b9)) & MASK64
    if schedule != "iid":
        state ^= 0x10fade

    def unit():
        nonlocal state
        state = (state + 0x9e3779b97f4a7c15) & MASK64
        value = state
        value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        value ^= value >> 31
        return float(value >> 11) * (1.0 / 9007199254740992.0)

    loss = float(100000 if schedule == "iid" else 500000) / 1000000.0
    ids, attempted, remaining = [], 0, 0
    while len(ids) < 10 and attempted < CANDIDATE_LIMIT:
        if attempted % 256 == 0:
            _check(budget)
        candidate = attempted
        attempted += 1
        if schedule == "burst":
            if remaining:
                remaining -= 1
                drop = True
            else:
                drop = unit() < loss / (8.0 - 7.0 * loss)
                if drop:
                    remaining = 7
        else:
            drop = unit() < loss
        if not drop:
            packet_id = (MASK32 - 2 * candidate if schedule == "adversarial" else
                         6 + candidate if schedule == "repair-only" else candidate)
            ids.append(packet_id & MASK32)
    require(len(ids) == 10, "frozen trace candidate cap")
    require(len(set(ids)) == 10, "duplicate delivered packet ID")
    return dict(ids=ids, attempted_candidates=attempted,
                trace_sha256=hashlib.sha256(b"".join(struct.pack("<I", x) for x in ids)).hexdigest())


class CoefficientLedger:
    """Bind every ordered ID/row visit, including repeated IDs across traces."""
    def __init__(self, mapper, budget=None):
        self.mapper = mapper
        self.budget = budget
        self.visits = 0
        self.digest = hashlib.sha256()

    def columns(self, ids):
        require(type(ids) in (list, tuple) and 6 <= len(ids) <= 10, "coefficient ID roster")
        for packet_id in ids:
            _integer(packet_id, 0, MASK32, "packet ID range")
        require(len(set(ids)) == len(ids), "duplicate coefficient packet ID")
        result = []
        for packet_id in ids:
            _check(self.budget)
            row = self.mapper.row(packet_id)
            require(type(row) in (list, tuple, bytes) and len(row) == 6 and
                    all(type(x) is int and 0 <= x < 256 for x in row) and any(row),
                    "invalid or zero coefficient row")
            row = tuple(row)
            self.digest.update(struct.pack("<I", packet_id) + bytes(row))
            self.visits += 1
            result.append(row)
        return tuple(result)


def _rank(ranker, columns):
    return _integer(ranker(columns), 1, 6, "rank must be an exact nonzero integer in 1..6")


def _valid_ranks(ranks):
    require(type(ranks) is list and len(ranks) == 5, "five nested prefix ranks")
    for value in ranks:
        _integer(value, 1, 6, "prefix rank range")
    require(all(0 <= right - left <= 1 for left, right in zip(ranks, ranks[1:])),
            "nonmonotone or impossible nested rank")


def _trace_record(B, root, schedule, ledger, ranker, budget):
    record = trace(B, root, schedule, budget)
    columns = ledger.columns(record["ids"])
    record["ranks"] = [_rank(ranker, columns[:length]) for length in range(6, 11)]
    _valid_ranks(record["ranks"])
    return record


def _trace_summary(records, zero_required=False):
    failure = {str(overhead): 0 for overhead in OVERHEADS}
    first = {str(overhead): 0 for overhead in range(5)}
    first["beyond4"] = 0
    count = 0
    for record in records:
        ranks = record["ranks"]
        _valid_ranks(ranks)
        count += 1
        for overhead in OVERHEADS:
            failure[str(overhead)] += int(ranks[overhead] < 6)
        success = next((str(i) for i, rank in enumerate(ranks) if rank == 6), "beyond4")
        first[success] += 1
    passed = failure["0"] == 0 if zero_required else count > 0 and 100 * failure["0"] <= count
    return dict(traces=count, failure_by_overhead=failure, first_success=first, oh0_pass=passed)


def summarize(fresh, hard, historical, fixtures, expected_fixtures):
    """Empirical per-cell gates; no confidence claim or historical-failure veto."""
    summary = _trace_summary(fresh)
    cells = []
    for B in WIDTHS:
        for schedule in SCHEDULES:
            rows = [row for row in fresh if row["B"] == B and row["schedule"] == schedule]
            cell = dict(B=B, schedule=schedule, **_trace_summary(rows))
            cells.append(cell)
    summary["cells"] = cells
    require(len(fixtures) == len(expected_fixtures), "fixture summary coverage")
    matched = 0
    for i, (actual, expected) in enumerate(zip(fixtures, expected_fixtures)):
        require(type(actual["index"]) is int and actual["index"] == i, "fixture summary order")
        _integer(actual["rank"], 1, 6, "fixture summary rank")
        _integer(expected["expected_rank"], 1, 6, "fixture expected rank")
        matched += int(actual["rank"] == expected["expected_rank"])
    repaired = 0
    for i, row in enumerate(historical):
        require(type(row["index"]) is int and row["index"] == i, "history summary order")
        _integer(row["rank"], 1, 6, "historical summary rank")
        require(type(row["repaired"]) is bool and row["repaired"] == (row["rank"] == 6),
                "history repaired flag")
        repaired += int(row["repaired"])
    return dict(overheads=list(OVERHEADS), fresh=summary, hard=_trace_summary(hard, True),
                history=dict(prefixes=len(historical), repaired=repaired, unrepaired=len(historical) - repaired),
                fixtures=dict(checks=len(fixtures), matched=matched))


def _empty_result():
    return dict(protocol=PROTOCOL, outcome="INVALID", inputs=None, inputs_sha256=None,
                roots=[], fresh=[], hard=[], history=[], fixtures=[], evidence={}, counts={}, summary={})


def evaluate_recovery(inputs, roots, fresh_specs, hard_specs, mapper, ranker=None,
                      budget=None, record=None):
    """Evaluate explicitly injected data; this helper never selects roots/pairs.

    Small unrelated domains and mappers are allowed for neutral tests. The
    actual run_screen entry point separately enforces every frozen count.
    """
    result = _empty_result() if record is None else record
    require(result == _empty_result(), "nonempty recovery output record")
    ranker = RADIX.rank_columns if ranker is None else ranker
    require(type(inputs) is dict and type(roots) is list and len(roots) > 0, "recovery inputs/roots")
    for root in roots:
        _root(root)
    require(len(set(roots)) == len(roots), "duplicate fresh root")
    require(not set(roots) & set(inputs["excluded_roots"]), "fresh root intersects excluded domain")
    result["inputs"], result["roots"] = inputs, list(roots)
    result["inputs_sha256"] = RADIX.digest(canonical(inputs))
    result["evidence"]["field"] = RADIX.init_field()
    result["evidence"]["lookup"] = mapper.evidence()
    ledger = CoefficientLedger(mapper, budget)
    seen = set()
    for root_index, B, schedule in fresh_specs:
        _check(budget)
        _integer(root_index, 0, len(roots) - 1, "fresh root index")
        key = root_index, B, schedule
        require(key not in seen, "duplicate fresh trace cell")
        seen.add(key)
        row = dict(root_index=root_index, B=B, schedule=schedule,
                   **_trace_record(B, roots[root_index], schedule, ledger, ranker, budget))
        result["fresh"].append(row)
    seen = set()
    for phase, root_index, B, schedule, root in hard_specs:
        _check(budget)
        require(type(phase) is str and phase in ("training", "validation"), "hard phase")
        _integer(root_index, 0, 2, "hard root index")
        require(type(B) is int and B == 2 and schedule in SCHEDULES[1:], "hard trace coordinates")
        key = phase, root_index, B, schedule
        require(key not in seen, "duplicate hard trace cell")
        seen.add(key)
        row = dict(phase=phase, root_index=root_index, B=B, schedule=schedule,
                   **_trace_record(B, root, schedule, ledger, ranker, budget))
        result["hard"].append(row)
    for index, prefix in enumerate(inputs["prefixes"]):
        _check(budget)
        require(type(prefix["K"]) is int and prefix["K"] == 6, "historical K")
        overhead = _integer(prefix["overhead"], 0, 1, "historical overhead")
        require(len(prefix["ids"]) == 6 + overhead, "historical prefix length")
        rank = _rank(ranker, ledger.columns(prefix["ids"]))
        result["history"].append(dict(index=index, rank=rank, repaired=rank == 6))
    for index, fixture in enumerate(inputs["fixtures"]):
        _check(budget)
        expected = _integer(fixture["expected_rank"], 5, 6, "fixture expected rank")
        rank = _rank(ranker, ledger.columns(fixture["ids"]))
        result["fixtures"].append(dict(index=index, rank=rank))
        require(rank == expected, "sealed .62 fixture discrepancy at index " + str(index))
    result["evidence"].update(coefficient_visits=ledger.visits,
                               coefficient_visit_sha256=ledger.digest.hexdigest())
    result["counts"] = dict(fresh_roots=len(roots), fresh_traces=len(result["fresh"]),
                            hard_traces=len(result["hard"]), history_prefixes=len(result["history"]),
                            history_origins=len(inputs["origins"]), fixture_checks=len(result["fixtures"]),
                            rank_checks=5 * (len(result["fresh"]) + len(result["hard"])) +
                            len(result["history"]) + len(result["fixtures"]), coefficient_visits=ledger.visits)
    result["summary"] = summarize(result["fresh"], result["hard"], result["history"],
                                  result["fixtures"], inputs["fixtures"])
    require(RADIX.digest(canonical(inputs)) == result["inputs_sha256"], "sealed recovery inputs changed")
    _check(budget)
    result["outcome"] = "PASS" if (result["summary"]["fresh"]["oh0_pass"] and
                                   all(cell["oh0_pass"] for cell in result["summary"]["fresh"]["cells"]) and
                                   result["summary"]["hard"]["oh0_pass"]) else "FAIL"
    return result


def run_screen():
    """Sole actual-campaign entry: bind the already-selected pair, never choose."""
    result = _empty_result()
    budget = RADIX.Budget()
    try:
        inputs = load_inputs(deadline=budget.deadline)
        pair_data = inputs["pair"]
        pair = tuple(tuple(tuple(row) for row in pair_data[key]) for key in ("A0", "A1"))
        require(pair_data["pair_sha256"] == PAIR_SHA256 and
                RADIX.digest(b"".join(RADIX.matrix_bytes(matrix) for matrix in pair)) == PAIR_SHA256,
                "sealed pair bytes/hash")
        require(len(inputs["excluded_roots"]) == 297 and len(set(inputs["excluded_roots"])) == 297 and
                len(inputs["prefixes"]) == 82 and len(inputs["origins"]) == 187 and
                len(inputs["fixtures"]) == 33, "frozen history/fixture counts")
        roots = fresh_roots()
        require(len(roots) == 512, "frozen fresh root count")
        for root in roots:
            _root(root)
        require(len(set(roots)) == 512 and not set(roots) & set(inputs["excluded_roots"]),
                "frozen fresh root collision")
        mapper = TM.RowMapper(pair, budget)
        evaluate_recovery(inputs, roots, fresh_domain(512), hard_domain(), mapper,
                          budget=budget, record=result)
        require(result["counts"] == dict(fresh_roots=512, fresh_traces=6144, hard_traces=18,
                history_prefixes=82, history_origins=187, fixture_checks=33, rank_checks=30925,
                coefficient_visits=62347), "frozen complete recovery accounting")
        require(all(cell["traces"] == 512 for cell in result["summary"]["fresh"]["cells"]),
                "frozen per-cell denominator")
        budget.check()
    except Exception as error:
        result["outcome"] = "INVALID"
        result["error"] = (type(error).__name__ + ": " + str(error))[:1024]
    return result


def main(argv):
    if argv != ["--worker"]:
        sys.stderr.write("usage: Wh2ThueMorseRecoveryR0.py --worker\n")
        return 2
    try:
        resource.setrlimit(resource.RLIMIT_AS, (512 * 1024 * 1024, 512 * 1024 * 1024))
        result = run_screen()
        raw = canonical(result) + b"\n"
        require(len(raw) <= RADIX.STDOUT_LIMIT, "worker JSON byte cap")
        sys.stdout.buffer.write(raw)
        sys.stdout.buffer.flush()
        return 1 if result["outcome"] == "INVALID" else 0
    except Exception as error:
        sys.stderr.write(("thue-morse-recovery-r0: " + type(error).__name__ + ": " + str(error))[:1024] + "\n")
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
