#!/usr/bin/env python3
"""Neutral .63 tests and explicit, post-result independent recovery replay.

No construction or fresh campaign root is evaluated on import. The replay
uses the independently reviewed 8/8/8/8 mapper and polynomial-long-division
GF256 arithmetic, not the worker's 10/7/7/8 mapper or rank implementation.
"""
import hashlib
import importlib.util
import io
import itertools
import json
from pathlib import Path
import struct
import sys
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
PROTOCOL = "wirehair.wh2.thue-morse-recovery-r0"
MASK64 = (1 << 64) - 1
MAX_ID = (1 << 32) - 1
WIDTHS = (2, 64, 1280)
SCHEDULES = ("iid", "burst", "adversarial", "repair-only")
OVERHEADS = (0, 1, 2, 4)
TRAINING = ("0xd1b54a32d192ed03", "0x94d049bb133111eb", "0x8538ecb5bd456ea3")
# The original frozen benchmark contract's v7 validation roots, not v9's
# additional reserved roots (which belong only in the exclusion ledger).
VALIDATION = ("0xc0ac29b7c97c50dd", "0x3f84d5b5b5470917", "0x9216d5d98979fb1b")


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


old = sibling("_recovery_independent_thue", "test_Wh2ThueMorseR0.py")
arithmetic = old.oracle
canonical, require, equal = old.canonical, old.require, old.equal
OracleMapper, rank_columns = old.OracleMapper, old.rank_columns


def exact(actual, expected, message):
    # Python container equality alone accepts True==1 and 6.0==6.
    equal(canonical(actual), canonical(expected), message)


def integer(value, low, high, message):
    require(type(value) is int and low <= value <= high, message)
    return value


def root_value(value):
    require(type(value) is str and len(value) == 18 and value[:2] == "0x" and
            all(c in "0123456789abcdef" for c in value[2:]), "root encoding")
    return int(value[2:], 16)


def uniform_stream(seed):
    """Literal C++ SplitMix64 Next/Unit, including unsigned overflow."""
    integer(seed, 0, MASK64, "RNG seed")
    while True:
        seed = (seed + 0x9e3779b97f4a7c15) & MASK64
        value = ((seed ^ (seed >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        value ^= value >> 31
        yield (value >> 11) * (1.0 / 9007199254740992.0)


def id_hash(ids):
    for packet_id in ids:
        integer(packet_id, 0, MAX_ID, "packet ID")
    return hashlib.sha256(b"".join(struct.pack("<I", x) for x in ids)).hexdigest()


def frozen_trace(B, root, schedule, K=6, draws=None):
    """Independent loop; optional neutral draws also test exact burst length.

    The generic K argument exists only to replay unchanged C++ golden traces.
    The scientific replay calls this at K6 and the three frozen widths only.
    """
    integer(B, 1, MAX_ID, "trace width")
    integer(K, 2, 64000, "trace K")
    require(type(schedule) is str and schedule in SCHEDULES, "trace schedule")
    state = (root_value(root) ^ (K * 0x9e3779b97f4a7c15) ^
             (B * 0xbf58476d1ce4e5b9)) & MASK64
    if schedule != "iid":
        state ^= 0x10fade
    random = uniform_stream(state) if draws is None else iter(draws)
    probability = (100000 if schedule == "iid" else 500000) / 1000000.0
    start_probability = probability / (8.0 - 7.0 * probability)
    delivered, burst_end = [], -1
    for candidate in range(256 * (K + 4) + 65536):
        if schedule == "burst" and candidate < burst_end:
            continue  # No random draw while a previously started burst lasts.
        draw = next(random)
        require(type(draw) in (float, int) and 0 <= draw < 1, "RNG unit range")
        if schedule == "burst":
            if draw < start_probability:
                burst_end = candidate + 8
                continue
        elif draw < probability:
            continue
        packet_id = (MAX_ID - 2 * candidate if schedule == "adversarial" else
                     K + candidate if schedule == "repair-only" else candidate)
        delivered.append(packet_id & MAX_ID)
        if len(delivered) == K + 4:
            return dict(ids=delivered, attempted_candidates=candidate + 1,
                        trace_sha256=id_hash(delivered))
    raise ValueError("independent trace exceeded frozen candidate limit")


def replay_roots():
    """CAMPAIGN NAMESPACE: called only by explicit post-result replay."""
    return ["0x" + hashlib.sha256((PROTOCOL + ":fresh/" + str(index)).encode("ascii"))
            .hexdigest()[:16] for index in range(512)]


def hard_roster():
    for phase, roots in (("training", TRAINING), ("validation", VALIDATION)):
        for index, root in enumerate(roots):
            for schedule in SCHEDULES[1:]:
                yield phase, index, 2, schedule, root


def load_inputs():
    """Authenticate historical evidence without importing retired searches."""
    return sibling("_recovery_replay_history", "Wh2ThueMorseRecoveryHistoryR0.py").load_inputs()


def check_ranks(ranks):
    require(type(ranks) is list and len(ranks) == 5, "five prefix ranks required")
    for value in ranks:
        integer(value, 1, 6, "prefix rank must be an exact nonzero integer")
    require(all(a <= b <= a + 1 for a, b in zip(ranks, ranks[1:])),
            "prefix rank is not a monotone one-row extension")


def statistics(records, hard=False):
    failures = {str(h): 0 for h in OVERHEADS}
    first = {str(h): 0 for h in range(5)}
    first["beyond4"] = 0
    for record in records:
        ranks = record["ranks"]
        check_ranks(ranks)
        for h in OVERHEADS:
            failures[str(h)] += int(ranks[h] < 6)
        success = next((str(h) for h, rank in enumerate(ranks) if rank == 6), "beyond4")
        first[success] += 1
    return dict(traces=len(records), failure_by_overhead=failures, first_success=first,
                oh0_pass=(failures["0"] == 0 if hard else bool(records) and 100 * failures["0"] <= len(records)))


def summarize(fresh, hard, history, fixtures):
    cells = []
    for B in WIDTHS:
        for schedule in SCHEDULES:
            records = [row for row in fresh if row["B"] == B and row["schedule"] == schedule]
            cells.append(dict(B=B, schedule=schedule, **statistics(records)))
    overall = statistics(fresh)
    overall["cells"] = cells
    return dict(overheads=list(OVERHEADS), fresh=overall, hard=statistics(hard, hard=True),
                history=dict(prefixes=len(history), repaired=sum(row["repaired"] for row in history),
                             unrepaired=sum(not row["repaired"] for row in history)),
                fixtures=dict(checks=len(fixtures), matched=len(fixtures)))


def verify_result(report):
    """Explicit replay of complete PASS/FAIL evidence; never reselect a pair."""
    require(type(report) is dict, "recovery result object")
    equal(set(report), {"protocol", "outcome", "inputs", "inputs_sha256", "roots", "fresh",
                        "hard", "history", "fixtures", "evidence", "counts", "summary"}, "result keys")
    equal(report["protocol"], PROTOCOL, "recovery protocol")
    require(report["outcome"] in ("PASS", "FAIL"), "only complete scientific verdicts replay")
    inputs = load_inputs()
    exact(report["inputs"], inputs, "authenticated historical inputs differ")
    equal(report["inputs_sha256"], hashlib.sha256(canonical(inputs)).hexdigest(), "input seal")
    for key, count in (("prefixes", 82), ("origins", 187), ("excluded_roots", 297), ("fixtures", 33)):
        require(type(inputs[key]) is list and len(inputs[key]) == count, "input roster: " + key)
    roots = replay_roots()
    exact(report["roots"], roots, "fresh root derivation/order")
    require(len(roots) == 512 and len(set(roots)) == 512, "fresh root collision")
    for root in roots + inputs["excluded_roots"]:
        root_value(root)
    require(len(set(inputs["excluded_roots"])) == 297 and
            not set(roots).intersection(inputs["excluded_roots"]), "freshness exclusion collision")
    for key, count in (("fresh", 6144), ("hard", 18), ("history", 82), ("fixtures", 33)):
        require(type(report[key]) is list and len(report[key]) == count, "incomplete roster: " + key)

    pair = tuple(tuple(tuple(old.byte(x) for x in row) for row in inputs["pair"][key])
                 for key in ("A0", "A1"))
    require(all(len(matrix) == 6 and all(len(row) == 6 for row in matrix) for matrix in pair),
            "authenticated matrix shape")
    mapper = OracleMapper(pair)
    visit_digest, visits, rank_checks = hashlib.sha256(), 0, 0

    def columns(ids):
        nonlocal visits
        result = []
        for packet_id in ids:
            integer(packet_id, 0, MAX_ID, "coefficient packet ID")
            column = mapper.row(packet_id)
            require(type(column) in (tuple, list) and len(column) == 6 and
                    all(type(x) is int and 0 <= x <= 255 for x in column) and any(column),
                    "zero or invalid coefficient")
            visits += 1
            visit_digest.update(struct.pack("<I", packet_id) + bytes(column))
            result.append(column)
        return result

    def trace_record(root_index, B, schedule, root, recorded, phase=None):
        nonlocal rank_checks
        expected = dict(root_index=root_index, B=B, schedule=schedule,
                        **frozen_trace(B, root, schedule))
        if phase is not None:
            expected["phase"] = phase
        cols = columns(expected["ids"])
        expected["ranks"] = [rank_columns(cols[:length]) for length in range(6, 11)]
        check_ranks(expected["ranks"])
        rank_checks += 5
        exact(recorded, expected, "independent trace/rank differs")
        return expected

    fresh, cursor = [], 0
    for root_index, root in enumerate(roots):
        for B in WIDTHS:
            for schedule in SCHEDULES:
                fresh.append(trace_record(root_index, B, schedule, root, report["fresh"][cursor]))
                cursor += 1
    hard = [trace_record(index, B, schedule, root, recorded, phase)
            for (phase, index, B, schedule, root), recorded in zip(hard_roster(), report["hard"])]
    equal(len(hard), 18, "hard projection roster")
    history, fixtures = [], []
    for index, prefix in enumerate(inputs["prefixes"]):
        rank = rank_columns(columns(prefix["ids"]))
        integer(rank, 1, 6, "history rank")
        rank_checks += 1
        row = dict(index=index, rank=rank, repaired=(rank == 6))
        exact(report["history"][index], row, "history prefix differs")
        history.append(row)
    for index, fixture in enumerate(inputs["fixtures"]):
        rank = rank_columns(columns(fixture["ids"]))
        integer(rank, 1, 6, "fixture rank")
        rank_checks += 1
        equal(rank, fixture["expected_rank"], "known fixed-pair fixture differs")
        row = dict(index=index, rank=rank)
        exact(report["fixtures"][index], row, "fixture record differs")
        fixtures.append(row)
    evidence = dict(field=arithmetic.field_evidence(), lookup=old.lookup_evidence(pair),
                    coefficient_visits=visits, coefficient_visit_sha256=visit_digest.hexdigest())
    exact(report["evidence"], evidence, "arithmetic/lookup/coefficient evidence differs")
    counts = dict(fresh_roots=512, fresh_traces=6144, hard_traces=18, history_prefixes=82,
                  history_origins=187, fixture_checks=33, rank_checks=30925, coefficient_visits=62347)
    equal(visits, 62347, "coefficient visit accounting")
    equal(rank_checks, 30925, "rank check accounting")
    exact(report["counts"], counts, "result counts differ")
    summary = summarize(fresh, hard, history, fixtures)
    exact(report["summary"], summary, "aggregate recovery evidence differs")
    passed = (summary["fresh"]["oh0_pass"] and summary["hard"]["oh0_pass"] and
              all(cell["oh0_pass"] for cell in summary["fresh"]["cells"]))
    equal(report["outcome"], "PASS" if passed else "FAIL", "scientific gate differs")
    return dict(protocol=PROTOCOL, verified=True, outcome=report["outcome"], counts=counts, summary=summary)


def replay_file(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(4 * 1024 * 1024 + 1)
    require(len(raw) <= 4 * 1024 * 1024, "replay byte cap")
    report = json.loads(raw, object_pairs_hook=arithmetic.reject_duplicate_keys,
                        parse_constant=lambda value: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    equal(raw, canonical(report) + b"\n", "noncanonical result encoding")
    return verify_result(report)


class NeutralMapper:
    """Six nonzero unit columns; never a campaign equation law."""
    def row(self, packet_id):
        return tuple(int(i == packet_id % 6) for i in range(6))

    def evidence(self):
        return {"neutral_lookup": True}


def neutral_inputs():
    pair = old.neutral_companions()
    prefixes = [dict(K=6, ids=list(range(1000 + i * 10, 1006 + i * 10 + int(i >= 78))),
                     original_widths=[2], overhead=int(i >= 78)) for i in range(82)]
    fixtures = [dict(family="neutral", index=i,
                     ids=list(range(10000 + 10 * i, 10000 + 10 * i + (6 if i < 16 else 7 if i == 16 else 8))),
                     expected_rank=5 if i < 17 else 6) for i in range(33)]
    return dict(provenance={"neutral": True}, prefixes=prefixes,
                origins=[{"neutral_origin": i} for i in range(187)],
                excluded_roots=["0x{:016x}".format(0x9000000000000000 + i) for i in range(297)],
                pair=dict(A0=pair[0], A1=pair[1], pair_sha256="neutral", selection_sha256="neutral"),
                fixtures=fixtures)


def neutral_ranker(fresh_failures=0, hard_failure=False):
    """Scripted ranks for full-roster accounting, not coefficient scoring."""
    calls = [0]

    def rank(columns):
        index = calls[0]
        calls[0] += 1
        if index < 6144 * 5:
            trace, prefix = divmod(index, 5)
            # First B2/IID cell only: six failures pass pooled but fail its cell.
            return 5 if trace % 12 == 0 and trace // 12 < fresh_failures and prefix == 0 else 6
        if index < 6162 * 5:
            return 5 if hard_failure and index == 6144 * 5 else 6
        if index < 6162 * 5 + 82:
            return 5 if index == 6162 * 5 else 6
        return 5 if index < 6162 * 5 + 82 + 17 else 6
    return calls, rank


class RecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = sibling("_recovery_worker_under_test", "Wh2ThueMorseRecoveryR0.py")

    def setUp(self):
        self.module = sys.modules[__name__]
        # All tests are neutral. Nested mock scopes may replace these with
        # explicit unrelated roots/data, never with the frozen construction.
        forbidden = ((self.module, "replay_roots"), (self.module, "load_inputs"),
                     (self.worker, "fresh_roots"), (self.worker, "run_screen"),
                     (self.worker, "load_inputs"), (self.worker.TM, "fixed_feedback"),
                     (self.worker.TM, "select_parameter"), (self.worker.TM, "run_screen"),
                     (self.worker.TM, "RowMapper"),
                     (self.worker.RADIX, "run_screen"), (self.worker.RADIX, "candidate_matrix"),
                     (old, "geometric_feedback"), (old, "selection_evidence"),
                     (old, "verify_result"), (arithmetic, "derive_recorded_candidate"))
        for module, name in forbidden:
            patch = mock.patch.object(module, name,
                                      side_effect=AssertionError("actual roots/construction/campaign forbidden"))
            patch.start()
            self.addCleanup(patch.stop)

    def test_cpp_trace_goldens_and_le32_hash(self):
        goldens = (
            ("iid", [0, 1, 3, 4, 5, 6], 7,
             "37391e3501b8cb36bf90f089771535ab6824a2ac518e0c00ac7c600746b5c3ab"),
            ("burst", list(range(6)), 6,
             "cd9a54ed1f18bf97db08914e280ea7349e11ca2c4885a4d8052552ceba84208d"),
            ("adversarial", [4294967295, 4294967291, 4294967287, 4294967283, 4294967281, 4294967271], 13,
             "e5ce1cd184f9ff951dcfcf196f1656baa15fadb745b4dd94f05fb5d1e10375ae"),
            ("repair-only", [2, 4, 6, 8, 9, 14], 13,
             "ff58d29f143b0d1d52f6a6f25b0c2ef53ed4a123898639918cf72091ea6d3b5a"))
        for schedule, ids, attempted, digest in goldens:
            self.assertEqual(frozen_trace(2, TRAINING[0], schedule, K=2),
                             dict(ids=ids, attempted_candidates=attempted, trace_sha256=digest))
        self.assertEqual(id_hash([0, 0x12345678, MAX_ID]),
                         "cdf516d8944dafbd0e8a427e7517d767a63667051927d23e2d5380d6a6d86921")

    def test_neutral_width_schedule_traces_match_worker(self):
        root = "0x0123456789abcdef"
        for B in WIDTHS:
            for schedule in SCHEDULES:
                expected = frozen_trace(B, root, schedule)
                self.assertEqual(self.worker.trace(B, root, schedule), expected)
                self.assertEqual(len(set(expected["ids"])), 10)
                self.assertLessEqual(expected["attempted_candidates"], 68096)
        for invalid in (True, 2.0, 0, -1):
            with self.assertRaises(ValueError):
                frozen_trace(invalid, root, "iid")
        for invalid in ("0X0123456789abcdef", "0x123", 1, True):
            with self.assertRaises(ValueError):
                frozen_trace(2, invalid, "iid")

    def test_exact_width_seed_and_iid_salt_boundary(self):
        root = "0x0123456789abcdef"
        for B in WIDTHS:
            mixed = (int(root, 16) ^ (6 * 0x9e3779b97f4a7c15) ^ (B * 0xbf58476d1ce4e5b9)) & MASK64
            for schedule in SCHEDULES:
                with mock.patch.object(self.module, "uniform_stream", return_value=itertools.repeat(0.99)) as rng:
                    frozen_trace(B, root, schedule)
                rng.assert_called_once_with(mixed if schedule == "iid" else mixed ^ 0x10fade)
        first = frozen_trace(2, root, "adversarial", draws=itertools.repeat(0.9))
        self.assertEqual(first["ids"], [MAX_ID - 2 * i for i in range(10)])

    def test_burst_drop_length_draw_consumption_and_single_full_cap(self):
        root = "0x0123456789abcdef"
        draws = mock.Mock(side_effect=[0.0] + [0.9] * 10)
        stream = iter(draws, None)
        actual = frozen_trace(2, root, "burst", draws=stream)
        self.assertEqual(actual["ids"], list(range(8, 18)))
        self.assertEqual(actual["attempted_candidates"], 18)
        self.assertEqual(draws.call_count, 11)
        consecutive = frozen_trace(2, root, "burst", draws=[0.0, 0.0] + [0.9] * 10)
        self.assertEqual(consecutive["ids"], list(range(16, 26)))
        boundary = frozen_trace(2, root, "burst", draws=itertools.repeat(0.5 / (8.0 - 7.0 * 0.5)))
        self.assertEqual(boundary["ids"], list(range(10)))  # Strict less-than, not <=.
        late = frozen_trace(2, root, "iid", draws=itertools.chain(itertools.repeat(0.0, 68086), itertools.repeat(0.9)))
        self.assertEqual(late["attempted_candidates"], 68096)
        self.assertEqual(late["ids"][-1], 68095)
        with self.assertRaises(ValueError):
            frozen_trace(2, root, "iid", draws=itertools.repeat(0.0))

    def test_rank_arithmetic_and_neutral_mapper(self):
        pair = old.neutral_companions()
        mapper = OracleMapper(pair)
        self.assertEqual([mapper.row(i) for i in range(6)], list(old.identity()))
        for packet_id in (255, 256, 1023, 1024, 65535, 65536, 1 << 31, MAX_ID):
            self.assertEqual(mapper.row(packet_id), old.direct_dyadic_row(pair, packet_id))
        columns = [mapper.row(i) for i in range(10)]
        self.assertEqual([rank_columns(columns[:n]) for n in range(6, 11)], [6] * 5)
        self.assertEqual(rank_columns([columns[0]] * 10), 1)
        self.assertEqual(rank_columns(old.identity()[:5]), 5)
        for wrong in ([5, 4, 6, 6, 6], [4, 6, 6, 6, 6], [True, 6, 6, 6, 6],
                      [6.0] * 5, [6] * 4, [7] * 5, [0] * 5):
            with self.assertRaises(ValueError):
                check_ranks(wrong)

    def test_first_success_keeps_overhead_three_and_beyond_four_separate(self):
        rows = [dict(ranks=ranks) for ranks in ([6] * 5, [5, 6, 6, 6, 6], [4, 5, 6, 6, 6],
                                              [3, 4, 5, 6, 6], [2, 3, 4, 5, 6], [5] * 5)]
        value = statistics(rows)
        self.assertEqual(value["first_success"], dict({str(i): 1 for i in range(5)}, beyond4=1))
        self.assertEqual(value["failure_by_overhead"], {"0": 5, "1": 4, "2": 3, "4": 1})
        self.assertFalse(value["oh0_pass"])
        self.assertTrue(statistics([dict(ranks=[5, 6, 6, 6, 6])] * 5 + [dict(ranks=[6] * 5)] * 507)["oh0_pass"])
        self.assertFalse(statistics([dict(ranks=[5, 6, 6, 6, 6])] * 6 + [dict(ranks=[6] * 5)] * 506)["oh0_pass"])

    def test_complete_domain_and_exact_hard_projection_roots(self):
        expected = [(i, B, schedule) for i in range(512) for B in WIDTHS for schedule in SCHEDULES]
        self.assertEqual(list(self.worker.fresh_domain(512)), expected)
        self.assertEqual(len(expected), 6144)
        self.assertEqual(list(self.worker.hard_domain()), list(hard_roster()))
        self.assertEqual(len(list(hard_roster())), 18)
        self.assertEqual(list(hard_roster())[-1], ("validation", 2, 2, "repair-only", VALIDATION[2]))

    def make_neutral_report(self, fresh_failures=6, hard_failure=False):
        inputs = neutral_inputs()
        roots = ["0x{:016x}".format(0x7000000000000000 + i) for i in range(512)]
        hard = [(phase, i, 2, schedule, "0x{:016x}".format(0xb000000000000000 + p * 3 + i))
                for p, phase in enumerate(("training", "validation"))
                for i in range(3) for schedule in SCHEDULES[1:]]
        calls, ranker = neutral_ranker(fresh_failures, hard_failure)
        result = self.worker.evaluate_recovery(inputs, roots, self.worker.fresh_domain(512), hard,
                                               NeutralMapper(), ranker=ranker)
        self.assertEqual(calls[0], 30925)
        return json.loads(canonical(result)), inputs, roots, hard

    def replay_neutral(self, report, inputs, roots, hard, fresh_failures=6, hard_failure=False,
                       mapper=None):
        calls, ranker = neutral_ranker(fresh_failures, hard_failure)
        with mock.patch.object(self.module, "load_inputs", return_value=inputs), \
                mock.patch.object(self.module, "replay_roots", return_value=roots), \
                mock.patch.object(self.module, "hard_roster", return_value=iter(hard)), \
                mock.patch.object(self.module, "OracleMapper", return_value=mapper or NeutralMapper()), \
                mock.patch.object(self.module, "rank_columns", side_effect=ranker), \
                mock.patch.object(old, "lookup_evidence", return_value=NeutralMapper().evidence()):
            result = verify_result(report)
        self.assertEqual(calls[0], 30925)
        return result

    def test_complete_roundtrip_rejects_cell_masking_and_finishes_all_work(self):
        report, inputs, roots, hard = self.make_neutral_report()
        self.assertEqual(report["outcome"], "FAIL")
        self.assertTrue(report["summary"]["fresh"]["oh0_pass"])
        self.assertFalse(report["summary"]["fresh"]["cells"][0]["oh0_pass"])
        self.assertEqual(report["summary"]["fresh"]["failure_by_overhead"]["0"], 6)
        self.assertEqual(report["summary"]["fresh"]["cells"][0]["traces"], 512)
        self.assertTrue(all(cell["oh0_pass"] for cell in report["summary"]["fresh"]["cells"][1:]))
        self.assertEqual(report["counts"]["coefficient_visits"], 62347)
        self.assertEqual(report["history"][0], dict(index=0, rank=5, repaired=False))
        self.assertEqual(report["summary"]["history"]["unrepaired"], 1)
        self.assertEqual(report["summary"]["fixtures"], dict(checks=33, matched=33))
        self.assertEqual(self.replay_neutral(report, inputs, roots, hard)["outcome"], "FAIL")
        self.assertLess(len(canonical(report)) + 1, 4 * 1024 * 1024)

        def altered(path, value):
            changed = json.loads(canonical(report))
            target = changed
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            return changed

        mutations = (
            (("protocol",), "other"), (("outcome",), "PASS"),
            (("fresh",), report["fresh"][:-1]), (("hard",), report["hard"][:-1]),
            (("history",), report["history"][:-1]), (("fixtures",), report["fixtures"][:-1]),
            (("fresh", 0, "root_index"), False), (("fresh", 0, "B"), 2.0),
            (("fresh", 0, "ranks", 0), 5.0), (("fresh", 0, "ids", 0), True),
            (("fresh", 0, "trace_sha256"), "0" * 64),
            (("fresh", 0, "attempted_candidates"), report["fresh"][0]["attempted_candidates"] + 1),
            (("fresh", 1), report["fresh"][0]), (("fresh", 0, "ranks"), [5, 4, 6, 6, 6]),
            (("inputs", "origins", 0, "neutral_origin"), False),
            (("inputs_sha256",), "0" * 64), (("roots", 1), roots[0]),
            (("counts", "rank_checks"), 30925.0),
            (("evidence", "coefficient_visits"), 62347.0),
            (("summary", "fresh", "cells", 0, "oh0_pass"), True),
            (("summary", "fresh", "first_success", "0"), 6138.0))
        for path, value in mutations:
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.replay_neutral(altered(path, value), inputs, roots, hard)

    def test_complete_pass_and_independent_hard_zero_gate(self):
        report, inputs, roots, hard = self.make_neutral_report(fresh_failures=5)
        self.assertEqual(report["outcome"], "PASS")
        self.assertEqual(self.replay_neutral(report, inputs, roots, hard, fresh_failures=5)["outcome"], "PASS")
        # The descriptive historical failure remains, and never enters a fresh
        # denominator or vetoes this otherwise passing neutral result.
        self.assertEqual(report["summary"]["history"]["unrepaired"], 1)
        report, inputs, roots, hard = self.make_neutral_report(fresh_failures=0, hard_failure=True)
        self.assertEqual(report["outcome"], "FAIL")
        self.assertTrue(report["summary"]["fresh"]["oh0_pass"])
        self.assertEqual(report["summary"]["hard"]["failure_by_overhead"]["0"], 1)
        self.assertEqual(self.replay_neutral(report, inputs, roots, hard,
                                           fresh_failures=0, hard_failure=True)["outcome"], "FAIL")

    def test_worker_invalid_coefficients_ranks_duplicates_and_fixture(self):
        inputs = neutral_inputs()
        inputs.update(prefixes=[], origins=[], fixtures=[])
        roots = ["0x0123456789abcdef"]
        specs = [(0, 2, "iid")]
        for wrong in ((0,) * 6, (True, 0, 0, 0, 0, 0), (1.0, 0, 0, 0, 0, 0), (1,) * 5):
            mapper = NeutralMapper()
            with mock.patch.object(mapper, "row", return_value=wrong), self.assertRaises(ValueError):
                self.worker.evaluate_recovery(inputs, roots, specs, [], mapper, ranker=lambda cols: 6)
        for wrong in (True, 6.0, 0, 7):
            with self.assertRaises(ValueError):
                self.worker.evaluate_recovery(inputs, roots, specs, [], NeutralMapper(), ranker=lambda cols: wrong)
        for ranks in ([5, 4, 6, 6, 6], [4, 6, 6, 6, 6]):
            with self.assertRaises(ValueError):
                self.worker.evaluate_recovery(inputs, roots, specs, [], NeutralMapper(), ranker=mock.Mock(side_effect=ranks))
        with self.assertRaises(ValueError):
            self.worker.evaluate_recovery(inputs, roots, specs * 2, [], NeutralMapper(), ranker=lambda cols: 6)
        with self.assertRaises(ValueError):
            self.worker.evaluate_recovery(inputs, roots * 2, specs, [], NeutralMapper(), ranker=lambda cols: 6)
        inputs["fixtures"] = [dict(ids=list(range(6)), expected_rank=5)]
        with self.assertRaises(ValueError):
            self.worker.evaluate_recovery(inputs, roots, specs, [], NeutralMapper(), ranker=lambda cols: 6)

    def test_replay_parser_rejects_noncanonical_duplicate_and_truncated_json(self):
        report = {"neutral": 1}
        with mock.patch.object(Path, "open", return_value=io.BytesIO(canonical(report) + b"\n")), \
                mock.patch.object(self.module, "verify_result", return_value={"verified": True}) as verify:
            self.assertEqual(replay_file("unused"), {"verified": True})
            verify.assert_called_once_with(report)
        for raw in (b'{"neutral":1}', b'{ "neutral":1}\n', b'{"x":1,"x":2}\n',
                    b'{"x":NaN}\n', b'{"neutral":', b'{}\ntrailer', b'x' * (4 * 1024 * 1024 + 1)):
            with self.subTest(raw=raw[:40]), mock.patch.object(Path, "open", return_value=io.BytesIO(raw)):
                with self.assertRaises(ValueError):
                    replay_file("unused")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--replay":
        try:
            print(canonical(replay_file(sys.argv[2])).decode("ascii"))
        except (ValueError, OSError, TypeError, KeyError, IndexError) as error:
            print("independent recovery replay: " + str(error)[:1000], file=sys.stderr)
            sys.exit(1)
    else:
        unittest.main()
