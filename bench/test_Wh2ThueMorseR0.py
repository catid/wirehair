#!/usr/bin/env python3
"""Independent neutral tests and explicit post-campaign Thue--Morse replay.

The fixed geometric polynomial is evaluated only by explicit replay. Unit
tests prohibit that entry point and never construct a campaign candidate.
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


ROOT = Path(__file__).resolve().parent
PROTOCOL = "wirehair.wh2.thue-morse-r0"
UINT32_MAX = (1 << 32) - 1


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


oracle = sibling("_thue_independent_arithmetic", "test_Wh2NoncommutingRadixR0.py")
PRODUCTS, INVERSES = oracle.PRODUCTS, oracle.INVERSES
identity, multiply, vector_product = oracle.identity, oracle.multiply, oracle.vector_product
rank_columns, inverse = oracle.rank_columns, oracle.inverse
canonical, require, equal = oracle.canonical, oracle.require, oracle.equal


def byte(value):
    require(type(value) is int and 0 <= value < 256, "not a GF256 byte")
    return value


def parity(value):
    return bin(value).count("1") % 2


def geometric_feedback():
    """Replay-only fixed polynomial; forbidden in all neutral unit tests."""
    coefficients = [1]
    for r in range(6):
        root = oracle.scalar_power(2, r)
        next_coefficients = [0] * (len(coefficients) + 1)
        for i, value in enumerate(coefficients):
            next_coefficients[i] ^= PRODUCTS[value][root]
            next_coefficients[i + 1] ^= value
        coefficients = next_coefficients
    require(coefficients[-1] == 1, "nonmonic independent polynomial")
    return tuple(coefficients[:-1])


def companion(feedback):
    require(len(feedback) == 6, "companion feedback length")
    feedback = tuple(byte(x) for x in feedback)
    return tuple(tuple(int(r == c + 1) if c < 5 else feedback[r]
                       for c in range(6)) for r in range(6))


def determinant(matrix):
    """Scalar determinant for interpolation, not polynomial expansion."""
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "determinant shape")
    value = 0
    for permutation in itertools.permutations(range(n)):
        term = 1
        for r, c in enumerate(permutation):
            term = PRODUCTS[term][byte(matrix[r][c])]
        value ^= term  # Signs coincide in characteristic two.
    return value


def polynomial_roster():
    for m in range(2, 5):
        words = tuple(word for word in itertools.product((0, 1), repeat=m - 1)
                      if m < 4 or len(set(word)) == 2)
        for columns in itertools.combinations(range(1, m + 5), m):
            masks = {tuple(t for t in range(1, m)
                           if t in columns and word[t - 1]) for word in words}
            for active in sorted(masks):
                if active:
                    yield m, columns, active


def complementary_minor(feedback, m, columns, active, parameter):
    coefficients = tuple(feedback) + (1,)
    rows = []
    for r in range(m):
        row = []
        for c in columns:
            value = coefficients[c - r] if 0 <= c - r <= 6 else 0
            if c == r and r in active:
                value ^= parameter
            row.append(value)
        rows.append(tuple(row))
    return tuple(rows)


def determinant_polynomial(feedback, m, columns, active):
    samples = tuple(determinant(complementary_minor(feedback, m, columns, active, x))
                    for x in (0, 1, 2))
    constant, at_one, at_two = samples
    linear_plus_quadratic = at_one ^ constant
    quadratic = PRODUCTS[at_two ^ constant ^ PRODUCTS[2][linear_plus_quadratic]][INVERSES[6]]
    coefficients = [constant, linear_plus_quadratic ^ quadratic, quadratic]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def polynomial_value(coefficients, x):
    value = 0
    for coefficient in reversed(coefficients):
        value = PRODUCTS[value][x] ^ coefficient
    return value


def neutral_companions():
    a = companion((1, 1, 0, 1, 0, 1))
    b = companion((3, 1, 0, 1, 0, 1))
    return a, b


def sequential_columns(generators, length, phase=0):
    product = identity()
    result = []
    for n in range(length):
        result.append(tuple(row[0] for row in product))
        product = multiply(product, generators[phase ^ parity(n)])
    return tuple(result)


class OracleMapper:
    """Independent 8/8/8/8 decomposition, not worker's 10/7/7/8 lookup."""
    def __init__(self, generators):
        blocks = [tuple(generators)]
        for _ in range(31):
            left, right = blocks[-1]
            blocks.append((multiply(left, right), multiply(right, left)))
        self.blocks = tuple(blocks)
        self.low = tuple(sequential_columns(generators, 256, phase) for phase in (0, 1))
        self.tables = {}
        for offset in (8, 16, 24):
            self.tables[offset] = tuple(tuple(self.segment(value, offset, phase)
                                              for value in range(256)) for phase in (0, 1))
        self.cache = {}

    def segment(self, value, offset, phase):
        result = identity()
        for bit in range(7, -1, -1):
            if value & (1 << bit):
                result = multiply(result, self.blocks[offset + bit][phase])
                phase ^= 1
        return result

    def row(self, packet_id):
        require(type(packet_id) is int and 0 <= packet_id <= UINT32_MAX,
                "oracle packet ID outside uint32")
        if packet_id not in self.cache:
            high, middle, lower, low = ((packet_id >> shift) & 255 for shift in (24, 16, 8, 0))
            phase_middle = parity(high)
            phase_lower = phase_middle ^ parity(middle)
            phase_low = phase_lower ^ parity(lower)
            value = self.low[phase_low][low]
            value = vector_product(self.tables[8][phase_lower][lower], value)
            value = vector_product(self.tables[16][phase_middle][middle], value)
            self.cache[packet_id] = vector_product(self.tables[24][0][high], value)
        return self.cache[packet_id]


def direct_dyadic_row(generators, packet_id):
    """Literal high-to-low right product for neutral arbitrary-ID checks."""
    blocks = [tuple(generators)]
    for _ in range(31):
        a, b = blocks[-1]
        blocks.append((multiply(a, b), multiply(b, a)))
    value = identity()
    phase = 0
    for bit in range(31, -1, -1):
        if packet_id & (1 << bit):
            value = multiply(value, blocks[bit][phase])
            phase ^= 1
    return tuple(row[0] for row in value)


def final_roster():
    for index, (parameter, ids) in enumerate(oracle.ordinary_roster()):
        yield "ordinary", index, parameter, ids, (6, 7, 8, 9)
    for index, (parameter, ids) in enumerate(oracle.modular_roster()):
        yield "modular", index, parameter, ids, (6, 7, 8, 10, 17)
    yield "old-companion", 0, None, tuple(6 + t * 16777217 for t in range(10)), (6, 7, 8, 10)
    for r in range(29):
        stride = 1 << r
        for index, start in enumerate((0, 6, UINT32_MAX - 9 * stride)):
            yield "dyadic", r * 3 + index, dict(exponent=r, start_index=index), tuple(start + t * stride for t in range(10)), (6, 7, 8, 10)


def selection_evidence(feedback):
    feedback = tuple(feedback)
    records, excluded = [], {0, feedback[0]}
    for m, columns, active in polynomial_roster():
        coefficients = determinant_polynomial(feedback, m, columns, active)
        require(coefficients[0] != 0 and len(coefficients) <= 3,
                "independent local polynomial contradicts proof")
        roots = [x for x in range(256) if polynomial_value(coefficients, x) == 0]
        excluded.update(roots)
        records.append(dict(m=m, columns=columns, active=active,
                            degree_bound=len(active), coefficients=coefficients, roots=roots))
    equal(len(records), 190, "independent polynomial count")
    equal(sum(row["degree_bound"] for row in records), 240, "independent degree bound")
    eligible = [x for x in range(256) if x not in excluded]
    require(len(eligible) >= 14, "independent root-union existence contradiction")
    parameter = eligible[0]
    changed = (feedback[0] ^ parameter,) + feedback[1:]
    pair = companion(feedback), companion(changed)
    pair_bytes = bytes(x for matrix in pair for row in matrix for x in row)
    return dict(feedback=feedback, polynomials=records, polynomial_count=190,
                degree_bound_sum=240, excluded=sorted(excluded), eligible=eligible,
                **{"lambda": parameter, "A0": pair[0], "A1": pair[1],
                   "pair_sha256": hashlib.sha256(pair_bytes).hexdigest()})


def local_recurrence_columns(pair, word):
    columns = list(identity())
    for t in range(4):
        feedback = tuple(row[-1] for row in pair[int(word[t])])
        columns.append(tuple(oracle.dot(feedback, tuple(columns[t + k][r] for k in range(6)))
                             for r in range(6)))
    return tuple(columns)


def local_evidence(pair):
    # Derive the four-bit factors from aligned substitution rather than
    # treating a literal worker roster as the oracle.
    words = {"".join("01" if bit == 0 else "10" for bit in pair_bits)
             for pair_bits in itertools.product((0, 1), repeat=2)}
    words.update("{}{}{}{}".format(1 - a, b, 1 - b, c)
                 for a, b, c in itertools.product((0, 1), repeat=3)
                 if len({a, b, c}) == 2)
    equal(len(words), 10, "independent four-bit factor count")
    result = []
    for word in sorted(words):
        columns = local_recurrence_columns(pair, word)
        for offsets in itertools.combinations(range(10), 6):
            require(rank_columns(tuple(columns[i] for i in offsets)) == 6,
                    "independent local certificate disagrees with selector")
        result.append(dict(word=word, passed=210))
    return result


def cyclic_evidence(pair):
    # The replay reconstructs companions with D=lambda*e0*e5^T. Its left
    # Krylov columns A0^j*(lambda*e0)=lambda*e_j are already a basis; checking
    # the right Krylov rows below therefore proves all 36 outer products span
    # Mat_6. This is not a general algebra test for arbitrary matrix pairs.
    first, second = pair
    difference = tuple(tuple(a ^ b for a, b in zip(left, right))
                       for left, right in zip(first, second))
    rows, current = [], difference
    for _ in range(6):
        rows.append(current[0])
        current = multiply(current, first)
    rank = rank_columns(zip(*rows))
    invertible = [inverse(matrix) is not None for matrix in pair]
    noncommuting = multiply(first, second) != multiply(second, first)
    require(rank == 6 and all(invertible) and noncommuting,
            "independent cyclic/full-algebra proof disagreement")
    return dict(cyclic_difference_rows=rows, cyclic_difference_rank=rank,
                invertible=invertible, noncommuting=noncommuting)


def lookup_evidence(pair):
    blocks = [tuple(pair)]
    for _ in range(31):
        left, right = blocks[-1]
        blocks.append((multiply(left, right), multiply(right, left)))
    raw_matrix = lambda matrix: bytes(x for row in matrix for x in row)
    block_bytes = b"".join(raw_matrix(blocks[level][phase])
                           for phase in (0, 1) for level in range(32))
    tables, payload = [], []
    for name, start, width, phases in (("low10", 0, 10, (0, 1)),
                                       ("middle10", 10, 7, (0, 1)),
                                       ("middle17", 17, 7, (0, 1)),
                                       ("high24", 24, 8, (0,))):
        for phase in phases:
            if start == 0:
                raw = b"".join(bytes(column) for column in
                               sequential_columns(pair, 1 << width, phase))
            else:
                values = []
                for value in range(1 << width):
                    current, state = identity(), phase
                    for bit in range(width - 1, -1, -1):
                        if value & (1 << bit):
                            current = multiply(current, blocks[start + bit][state])
                            state ^= 1
                    values.append(raw_matrix(current))
                raw = b"".join(values)
            tables.append(dict(name=name, start_bit=start, width=width, phase=phase,
                               bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
            payload.append(raw)
    packed = b"".join(payload)
    equal(len(packed), 39936, "independent packed lookup bytes")
    return dict(packed_bytes=39936, max_matvecs=3, max_gf_products=108, max_xors=90,
                dyadic_blocks_sha256=hashlib.sha256(block_bytes).hexdigest(),
                tables=tables, tables_sha256=hashlib.sha256(packed).hexdigest())


def final_evidence(pair, recorded):
    require(type(recorded) is list and len(recorded) == 5620, "final roster length differs")
    mapper = OracleMapper(pair)
    seen_ids, seen_columns = set(), {}
    first_collision, visits = None, 0
    visit_hash = hashlib.sha256()

    def visit(packet_id):
        nonlocal first_collision, visits
        column = mapper.row(packet_id)
        visits += 1
        visit_hash.update(struct.pack("<I", packet_id) + bytes(column))
        seen_ids.add(packet_id)
        key = oracle.normalized(column)
        if key is None:
            if first_collision is None:
                first_collision = dict(kind="zero", ids=[packet_id], normalized=list(column))
        else:
            previous = seen_columns.setdefault(key, packet_id)
            if previous != packet_id and first_collision is None:
                first_collision = dict(kind="proportional", ids=[previous, packet_id],
                                       normalized=list(key))
        return column

    systematic = tuple(visit(i) for i in range(6))
    equal(systematic, identity(), "independent systematic identity differs")
    maximum = visit(UINT32_MAX)
    seams = []
    for exponent in (10, 17, 24, 31):
        start = (1 << exponent) - 5
        ids = tuple(range(start, start + 10))
        columns = tuple(visit(packet_id) for packet_id in ids)
        for offsets in itertools.combinations(range(10), 6):
            require(rank_columns(tuple(columns[i] for i in offsets)) == 6,
                    "independent dyadic seam contradicts local certificate")
        seams.append(dict(exponent=exponent, start=start, ids=ids,
                          coefficients=columns, passed=210))

    cursor, passed = 0, True
    for family, index, parameter, ids, prefixes in final_roster():
        columns = tuple(visit(packet_id) for packet_id in ids)
        ranks = {str(length): rank_columns(columns[:length]) for length in prefixes}
        expected = dict(family=family, index=index, parameter=parameter, ids=ids,
                        prefix_ranks=ranks, rank=ranks[str(len(ids))])
        equal(canonical(recorded[cursor]), canonical(expected),
              "independent final row differs at {}:{}".format(family, index))
        passed &= expected["rank"] == 6
        cursor += 1
    equal(cursor, 5620, "independent final roster accounting")
    equal(visits, 92475, "independent coefficient visit accounting")
    evidence = dict(systematic=dict(ids=list(range(6)), coefficients=systematic),
                    maximum=dict(id=UINT32_MAX, coefficients=maximum), seams=seams,
                    aliases=dict(visits=visits, unique_ids=len(seen_ids),
                                 normalized_rows=len(seen_columns), first_collision=first_collision,
                                 coefficient_visit_sha256=visit_hash.hexdigest()))
    return evidence, passed and first_collision is None


def verify_result(report):
    """Reconstruct the sole frozen choice; never select a final-pattern rescue."""
    require(type(report) is dict, "result must be an object")
    equal(set(report), {"protocol", "outcome", "selection", "selection_sha256", "local",
                        "final", "evidence", "counts"}, "result keys differ")
    equal(report["protocol"], PROTOCOL, "protocol differs")
    require(report["outcome"] in ("PASS", "FAIL"),
            "incomplete INVALID result cannot establish complete mathematical replay")
    require(type(report["selection"]) is dict and
            type(report["selection"].get("polynomials")) is list and
            len(report["selection"]["polynomials"]) == 190, "selection polynomial roster")
    selection = selection_evidence(geometric_feedback())
    equal(canonical(report["selection"]), canonical(selection), "selection evidence differs")
    equal(report["selection_sha256"], hashlib.sha256(canonical(selection)).hexdigest(),
          "selected pair seal differs")
    pair = selection["A0"], selection["A1"]
    equal(canonical(report["local"]), canonical(local_evidence(pair)), "local certificate differs")
    evidence = dict(field=oracle.field_evidence(), pair=cyclic_evidence(pair),
                    lookup=lookup_evidence(pair))
    final, passed = final_evidence(pair, report["final"])
    evidence.update(final)
    equal(canonical(report["evidence"]), canonical(evidence), "final evidence differs")
    outcome = "PASS" if passed else "FAIL"
    equal(report["outcome"], outcome, "result decision differs")
    counts = dict(polynomial_count=190, degree_bound_sum=240,
                  eligible_count=len(selection["eligible"]), local_rank_checks=2100,
                  seam_rank_checks=840, final_patterns=5620, final_visits=92428, total_visits=92475)
    equal(canonical(report["counts"]), canonical(counts), "work counts differ")
    return dict(protocol=PROTOCOL, verified=True, outcome=outcome, counts=counts)


def replay_file(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(4 * 1024 * 1024 + 1)
    require(len(raw) <= 4 * 1024 * 1024, "replay byte cap")
    report = json.loads(raw, object_pairs_hook=oracle.reject_duplicate_keys,
                        parse_constant=lambda value: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    equal(raw, canonical(report) + b"\n", "noncanonical result encoding")
    return verify_result(report)


class ThueMorseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = sibling("_thue_worker_under_test", "Wh2ThueMorseR0.py")
        cls.worker.RADIX.init_field()

    def setUp(self):
        forbidden = ((sys.modules[__name__], "geometric_feedback"),
                     (self.worker, "fixed_feedback"), (self.worker, "run_screen"),
                     (self.worker.RADIX, "candidate_matrix"),
                     (self.worker.RADIX, "run_screen"),
                     (oracle, "derive_recorded_candidate"))
        for module, name in forbidden:
            patch = mock.patch.object(module, name,
                                      side_effect=AssertionError("actual construction/campaign forbidden"))
            patch.start()
            self.addCleanup(patch.stop)

    def test_field_and_neutral_companion_orientation(self):
        equal(canonical(self.worker.RADIX.init_field()), canonical(oracle.field_evidence()), "field")
        self.assertEqual(len({oracle.scalar_power(2, i) for i in range(255)}), 255)
        pair = neutral_companions()
        for matrix, feedback in zip(pair, ((1, 1, 0, 1, 0, 1), (3, 1, 0, 1, 0, 1))):
            self.assertEqual(self.worker.companion(feedback), matrix)
            self.assertEqual(multiply(matrix, inverse(matrix)), identity())
        self.assertNotEqual(multiply(*pair), multiply(*reversed(pair)))
        self.assertEqual(cyclic_evidence(pair)["cyclic_difference_rank"], 6)
        for invalid in (True, 1.0, -1, 256):
            with self.assertRaises(ValueError):
                companion((invalid, 0, 0, 0, 0, 0))
            with self.assertRaises(ValueError):
                self.worker.companion((invalid, 0, 0, 0, 0, 0))

    def test_190_polynomial_roster_and_interpolation_on_neutral_feedback(self):
        roster = tuple(polynomial_roster())
        self.assertEqual(roster, tuple(self.worker.constraint_specs()))
        self.assertEqual(len(roster), 190)
        self.assertEqual(sum(len(active) for _, _, active in roster), 240)
        self.assertEqual(tuple(sum(len(active) for m, _, active in roster if m == wanted)
                               for wanted in (2, 3, 4)), (5, 40, 195))
        feedback = (1, 1, 0, 1, 0, 1)
        for m, columns, active in roster:
            coefficients = determinant_polynomial(feedback, m, columns, active)
            self.assertEqual(coefficients, self.worker.constraint_polynomial(feedback, m, columns, active))
            self.assertLessEqual(len(coefficients) - 1, len(active))
            for x in (0, 1, 2, 3, 7, 255):
                self.assertEqual(polynomial_value(coefficients, x),
                                 determinant(complementary_minor(feedback, m, columns, active, x)))
        for coefficients in ((1,), (1, 1), (1, 0, 1), (0, 1, 1), (7, 11, 23)):
            roots = [x for x in range(256) if polynomial_value(coefficients, x) == 0]
            self.assertEqual(self.worker.polynomial_roots(coefficients), roots)

    def test_mocked_polynomials_select_first_eligible_once(self):
        feedback = (1, 1, 0, 1, 0, 1)
        with mock.patch.object(self.worker, "constraint_polynomial", return_value=(1, 1)) as worker_poly, \
                mock.patch.object(sys.modules[__name__], "determinant_polynomial", return_value=(1, 1)) as oracle_poly:
            actual = self.worker.select_parameter(feedback)
            expected = selection_evidence(feedback)
        self.assertEqual(canonical(actual), canonical(expected))
        self.assertEqual(worker_poly.call_count, 190)
        self.assertEqual(oracle_poly.call_count, 190)
        self.assertEqual(actual["lambda"], 2)
        self.assertEqual(actual["excluded"], [0, 1])
        with mock.patch.object(self.worker, "constraint_polynomial", return_value=(0, 1)):
            with self.assertRaises(ValueError):
                self.worker.select_parameter(feedback)

    def test_local_recurrence_uses_right_products_and_all_four_bit_factors(self):
        pair = neutral_companions()
        sequence = sequential_columns(pair, 64)
        self.assertEqual(sequence[:6], identity())
        self.assertEqual(len(self.worker.WORDS), 10)
        observed = {"".join(str(parity(n + t)) for t in range(4)) for n in range(60)}
        self.assertEqual(observed, set(self.worker.WORDS))
        for word in self.worker.WORDS:
            self.assertEqual(local_recurrence_columns(pair, word), self.worker.local_columns(pair, word))
        product = identity()
        wrong = identity()[0]
        different = False
        for n in range(63):
            self.assertEqual(sequence[n], tuple(row[0] for row in product))
            wrong = vector_product(pair[parity(n)], wrong)
            product = multiply(product, pair[parity(n)])
            different |= wrong != sequence[n + 1]
        self.assertTrue(different, "left-update orientation negative control was vacuous")
        with mock.patch.object(sys.modules[__name__], "rank_columns", return_value=6) as rank:
            result = local_evidence(pair)
        self.assertEqual(rank.call_count, 2100)
        self.assertEqual(result, [dict(word=word, passed=210) for word in self.worker.WORDS])

    def test_independent_mapper_at_dyadic_seams_and_uint32_extremes(self):
        pair = neutral_companions()
        mapper, worker_mapper = OracleMapper(pair), self.worker.RowMapper(pair)
        sequence = sequential_columns(pair, 1030)
        for packet_id in range(1030):
            self.assertEqual(mapper.row(packet_id), sequence[packet_id])
            self.assertEqual(worker_mapper.row(packet_id), sequence[packet_id])
        ids = {UINT32_MAX, UINT32_MAX - 1, 0x12345678, 0xaaaaaaaa, 0x55555555}
        for bit in (8, 10, 16, 17, 24, 31):
            ids.update((1 << bit) + offset for offset in range(-5, 5))
        for packet_id in sorted(ids, reverse=True):
            expected = direct_dyadic_row(pair, packet_id)
            self.assertEqual(mapper.row(packet_id), expected)
            self.assertEqual(worker_mapper.row(packet_id), expected)
        for packet_id in (True, 1.0, -1, 1 << 32):
            with self.assertRaises(ValueError):
                mapper.row(packet_id)
            with self.assertRaises(ValueError):
                worker_mapper.row(packet_id)

    def test_lookup_and_cyclic_evidence_hashes_on_neutral_pair(self):
        pair = neutral_companions()
        mapper = self.worker.RowMapper(pair)
        expected = lookup_evidence(pair)
        self.assertEqual(canonical(mapper.evidence()), canonical(expected))
        self.assertEqual(sum(item["bytes"] for item in expected["tables"]), 39936)
        self.assertEqual(len(expected["tables"]), 7)
        self.assertEqual(canonical(self.worker.pair_evidence(pair, (1, 1, 0, 1, 0, 1), 2)),
                         canonical(cyclic_evidence(pair)))
        for chosen in (0, 1, True, 2.0, -1, 256):
            with self.assertRaises(ValueError):
                self.worker.pair_evidence(pair, (1, 1, 0, 1, 0, 1), chosen)

    def test_full_mocked_final_failure_keeps_every_roster_and_visit(self):
        class FakeMapper:
            def __init__(self, *args):
                pass

            def row(self, packet_id):
                return tuple(int(i == packet_id % 6) for i in range(6))

            reference_row = row

            def evidence(self):
                return {"neutral_lookup": True}

        def rank_fixture():
            calls = [0]

            def rank(columns):
                calls[0] += 1
                # All seams pass; the first ordinary full-nine row fails.
                return 5 if calls[0] == 844 else 6
            return calls, rank

        worker_calls, worker_rank = rank_fixture()
        with mock.patch.object(self.worker, "RowMapper", FakeMapper), \
                mock.patch.object(self.worker.RADIX, "rank_columns", side_effect=worker_rank):
            recorded, evidence, passed = self.worker.evaluate_final(neutral_companions(), mock.Mock())
        self.assertFalse(passed)
        self.assertEqual(worker_calls[0], 28540)
        self.assertEqual(len(recorded), 5620)
        self.assertEqual(sum(len(row["ids"]) for row in recorded), 92428)
        self.assertEqual(evidence["aliases"]["visits"], 92475)
        self.assertEqual(recorded[0]["rank"], 5)
        self.assertEqual(recorded[-1]["index"], 86)
        self.assertEqual(evidence["aliases"]["first_collision"]["ids"], [3, UINT32_MAX])

        independent_calls, independent_rank = rank_fixture()
        with mock.patch.object(sys.modules[__name__], "OracleMapper", FakeMapper), \
                mock.patch.object(sys.modules[__name__], "rank_columns", side_effect=independent_rank):
            independent, independent_passed = final_evidence(neutral_companions(), recorded)
        self.assertFalse(independent_passed)
        self.assertEqual(independent_calls[0], 28540)
        expected = {key: value for key, value in evidence.items() if key != "lookup"}
        self.assertEqual(canonical(independent), canonical(expected))

        changed = json.loads(canonical(recorded))
        changed[0]["ids"][0] = False
        changed_rank = json.loads(canonical(recorded))
        changed_rank[0]["rank"] = 5.0
        for wrong in (recorded[:-1], changed, changed_rank):
            _, ranks = rank_fixture()
            with mock.patch.object(sys.modules[__name__], "OracleMapper", FakeMapper), \
                    mock.patch.object(sys.modules[__name__], "rank_columns", side_effect=ranks):
                with self.assertRaises(ValueError):
                    final_evidence(neutral_companions(), wrong)

    def test_mock_replay_selection_seal_no_reselection_and_nested_types(self):
        feedback = (1, 1, 0, 1, 0, 1)
        with mock.patch.object(sys.modules[__name__], "determinant_polynomial", return_value=(1, 1)):
            selection = selection_evidence(feedback)
        local = [dict(word=word, passed=210) for word in self.worker.WORDS]
        lookup = {"neutral_lookup": True}
        final_observation = dict(aliases={"visits": 92475})
        report = dict(protocol=PROTOCOL, outcome="FAIL", selection=selection,
                      selection_sha256=hashlib.sha256(canonical(selection)).hexdigest(),
                      local=local, final=[{"neutral_final": True}],
                      evidence=dict(field=oracle.field_evidence(), pair=cyclic_evidence((selection["A0"], selection["A1"])),
                                    lookup=lookup, **final_observation),
                      counts=dict(polynomial_count=190, degree_bound_sum=240, eligible_count=254,
                                  local_rank_checks=2100, seam_rank_checks=840,
                                  final_patterns=5620, final_visits=92428, total_visits=92475))
        report = json.loads(canonical(report))
        with mock.patch.object(sys.modules[__name__], "geometric_feedback", return_value=feedback) as fixed, \
                mock.patch.object(sys.modules[__name__], "determinant_polynomial", return_value=(1, 1)) as polynomials, \
                mock.patch.object(sys.modules[__name__], "local_evidence", return_value=local), \
                mock.patch.object(sys.modules[__name__], "lookup_evidence", return_value=lookup), \
                mock.patch.object(sys.modules[__name__], "final_evidence", return_value=(final_observation, False)) as final:
            self.assertEqual(verify_result(report)["outcome"], "FAIL")
            self.assertEqual(fixed.call_count, 1)
            self.assertEqual(polynomials.call_count, 190)
            self.assertEqual(final.call_count, 1)

            def altered(path, value):
                changed = json.loads(canonical(report))
                target = changed
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                return changed

            mutations = [
                altered(("protocol",), "other"),
                altered(("outcome",), "PASS"),
                altered(("selection", "polynomials"), report["selection"]["polynomials"][:-1]),
                altered(("selection", "polynomials", 0, "m"), 2.0),
                altered(("selection", "polynomials", 0, "active"), [True]),
                altered(("selection", "polynomials", 0, "roots"), []),
                altered(("selection", "eligible"), list(range(3, 256))),
                altered(("selection", "lambda"), 3),
                altered(("selection", "A0", 0, 0), False),
                altered(("selection_sha256",), "0" * 64),
                altered(("local", 0, "passed"), 210.0),
                altered(("evidence", "pair", "cyclic_difference_rank"), 6.0),
                altered(("counts", "total_visits"), 92474),
            ]
            resealed = altered(("selection", "lambda"), 3)
            resealed["selection_sha256"] = hashlib.sha256(canonical(resealed["selection"])).hexdigest()
            mutations.append(resealed)
            missing = json.loads(canonical(report))
            del missing["local"]
            mutations.append(missing)
            for changed in mutations:
                with self.assertRaises(ValueError):
                    verify_result(changed)

    def test_final_roster_count_order_and_extreme_dyadic_bounds(self):
        expected = tuple(final_roster())
        self.assertEqual(canonical(expected), canonical(tuple(self.worker.final_patterns())))
        self.assertEqual(len(expected), 5620)
        self.assertEqual(sum(len(row[3]) for row in expected), 92428)
        self.assertEqual(sum(len(row[4]) for row in expected), 27700)
        for family, index, parameter, ids, _ in expected[5533:]:
            self.assertEqual(family, "dyadic")
            self.assertEqual(index, 3 * parameter["exponent"] + parameter["start_index"])
            self.assertEqual(len(set(ids)), 10)
            self.assertTrue(all(type(packet_id) is int and 0 <= packet_id <= UINT32_MAX for packet_id in ids))
            if parameter["start_index"] == 2:
                self.assertEqual(ids[-1], UINT32_MAX)

    def test_replay_file_rejects_noncanonical_duplicate_and_truncated_json(self):
        report = {"neutral": 1}
        with mock.patch.object(Path, "open", return_value=io.BytesIO(canonical(report) + b"\n")), \
                mock.patch.object(sys.modules[__name__], "verify_result", return_value={"verified": True}) as verify:
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
            print("independent Thue-Morse replay: " + str(error), file=sys.stderr)
            sys.exit(1)
    else:
        unittest.main()
