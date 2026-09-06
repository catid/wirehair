#!/usr/bin/env python3
"""Independent arithmetic tests and post-campaign replay for radix R0.

Unit tests never consume the frozen SHA candidate pool.  Real candidate
derivation belongs only to explicit replay of already published evidence.
"""
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import struct
import sys
import unittest
from unittest import mock


PROTOCOL = "wirehair.wh2.noncommuting-radix-r0"
ROOT = Path(__file__).resolve().parent


def polynomial_product(a, b):
    """Carryless polynomial product followed by long division, not xtime."""
    value = 0
    for bit in range(8):
        if b & (1 << bit):
            value ^= a << bit
    for bit in range(14, 7, -1):
        if value & (1 << bit):
            value ^= 0x14d << (bit - 8)
    return value


PRODUCTS = tuple(bytes(polynomial_product(a, b) for b in range(256))
                 for a in range(256))


def scalar_power(a, n):
    value = 1
    while n:
        if n & 1:
            value = PRODUCTS[value][a]
        a = PRODUCTS[a][a]
        n >>= 1
    return value


INVERSES = (0,) + tuple(scalar_power(a, 254) for a in range(1, 256))


def identity(size=6):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def dot(left, right):
    value = 0
    for a, b in zip(left, right):
        value ^= PRODUCTS[a][b]
    return value


def multiply(a, b):
    columns = tuple(zip(*b))
    return tuple(tuple(dot(row, column) for column in columns) for row in a)


def vector_product(matrix, column):
    return tuple(dot(row, column) for row in matrix)


def power(matrix, exponent):
    result = identity(len(matrix))
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix)
        matrix = multiply(matrix, matrix)
        exponent >>= 1
    return result


def insert_column(column, basis):
    """Highest-coordinate-first column-space elimination."""
    value = list(column)
    for pivot in range(len(value) - 1, -1, -1):
        if not value[pivot]:
            continue
        if pivot in basis:
            scale = PRODUCTS[value[pivot]]
            value = [a ^ scale[b] for a, b in zip(value, basis[pivot])]
        else:
            scale = PRODUCTS[INVERSES[value[pivot]]]
            basis[pivot] = tuple(scale[a] for a in value)
            return True
    return False


def rank_columns(columns):
    basis = {}
    for column in columns:
        insert_column(column, basis)
    return len(basis)


def inverse(matrix):
    """Track column combinations, then express each unit vector in that basis."""
    size = len(matrix)
    basis = {}
    for column_index, column in enumerate(zip(*matrix)):
        value = list(column)
        coefficients = [int(i == column_index) for i in range(size)]
        for pivot in range(size - 1, -1, -1):
            if not value[pivot]:
                continue
            if pivot in basis:
                scale = PRODUCTS[value[pivot]]
                old_value, old_coefficients = basis[pivot]
                value = [a ^ scale[b] for a, b in zip(value, old_value)]
                coefficients = [a ^ scale[b] for a, b in
                                zip(coefficients, old_coefficients)]
            else:
                scale = PRODUCTS[INVERSES[value[pivot]]]
                basis[pivot] = (tuple(scale[a] for a in value),
                                tuple(scale[a] for a in coefficients))
                break
        else:
            return None
    result_columns = []
    for column_index in range(size):
        target = [int(i == column_index) for i in range(size)]
        coefficients = [0] * size
        for pivot in range(size - 1, -1, -1):
            if target[pivot]:
                scale = PRODUCTS[target[pivot]]
                value, expression = basis[pivot]
                target = [a ^ scale[b] for a, b in zip(target, value)]
                coefficients = [a ^ scale[b] for a, b in
                                zip(coefficients, expression)]
        if any(target):
            raise AssertionError("independent inverse did not span a unit vector")
        result_columns.append(tuple(coefficients))
    return tuple(zip(*result_columns))


def vandermonde_column(phase):
    return tuple(scalar_power(2, (degree * phase) % 255) for degree in range(6))


def systematic_transform():
    return inverse(tuple(zip(*(vandermonde_column(j) for j in range(6)))))


def crossing_roster():
    # Enumerate whole six-subsets independently of the worker's leftmost loop.
    values = [row for row in itertools.combinations(range(-9, 9), 6)
              if row[0] < 0 <= row[-1] and row[-1] - row[0] <= 9]
    return tuple(sorted(values, key=lambda row: (-row[0], row[1:])))


def radix_value(digits):
    value = 0
    for digit in reversed(digits):
        value = value * 255 + digit
    return value


def ordinary_roster():
    for z in itertools.product(range(-2, 3), repeat=4):
        if radix_value(z) <= 0:
            continue
        initial = tuple(-120 * min(component, 0) for component in z)
        ids = tuple(radix_value(tuple(initial[i] + 15 * t * z[i]
                                      for i in range(4))) for t in range(9))
        yield z, ids


def modular_roster():
    # The first nonzero coordinate is fixed, avoiding worker-style filtering
    # over all 17**4 tuples while preserving exact lexicographic order.
    for leading in range(3, -1, -1):
        for tail in itertools.product(range(17), repeat=3 - leading):
            z = (0,) * leading + (1,) + tail
            ids = tuple(radix_value(tuple(15 * ((t * value) % 17)
                                           for value in z)) for t in range(17))
            yield z, ids


def independent_row(generators, packet_id, transform=None):
    phase = packet_id % 255
    block = packet_id // 255
    value = vandermonde_column(phase)
    for generator in generators:
        value = vector_product(power(generator, block % 255), value)
        block //= 255
    if block:
        raise ValueError("packet ID exceeds the four-digit block domain")
    return vector_product(systematic_transform() if transform is None else transform, value)


def load_worker():
    spec = importlib.util.spec_from_file_location("_noncommuting_radix_worker",
                                                ROOT / "Wh2NoncommutingRadixR0.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def neutral_generators():
    diagonal = tuple(tuple(scalar_power(2, i) if i == j else 0
                           for j in range(6)) for i in range(6))
    result = []
    for i in range(4):
        a = [list(row) for row in identity()]
        a[i][i + 1] = 1
        a = tuple(map(tuple, a))
        result.append(multiply(multiply(a, diagonal), inverse(a)))
    return tuple(result)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("ascii")


def matrix_hash(matrix):
    return hashlib.sha256(bytes(value for row in matrix for value in row)).hexdigest()


def derive_recorded_candidate(stage, attempt):
    """Replay-only derivation. Never called by unmocked unit tests."""
    prefix = PROTOCOL + ":stage/{}/candidate/{}/chunk/".format(stage, attempt)
    values = b"".join(hashlib.sha256((prefix + str(k)).encode("ascii")).digest()
                      for k in range(2))[:36]
    return tuple(tuple(values[6 * i:6 * i + 6]) for i in range(6))


class OracleMapper:
    def __init__(self, generators):
        self.transform = systematic_transform()
        self.tables = []
        for generator in generators:
            table = [identity()]
            for _ in range(254):
                table.append(multiply(table[-1], generator))
            self.tables.append(tuple(table))
        self.cache = {}

    def row(self, packet_id):
        if type(packet_id) is not int or not 0 <= packet_id <= 0xffffffff:
            raise ValueError("oracle packet ID outside uint32")
        if packet_id not in self.cache:
            value = vandermonde_column(packet_id % 255)
            block = packet_id // 255
            for table in self.tables:
                value = vector_product(table[block % 255], value)
                block //= 255
            if block:
                raise ValueError("oracle block digits not exhausted")
            self.cache[packet_id] = vector_product(self.transform, value)
        return self.cache[packet_id]


def normalized(column):
    for value in column:
        if value:
            scale = PRODUCTS[INVERSES[value]]
            return tuple(scale[x] for x in column)
    return None


def algebra_evidence(generators):
    queue = [identity()]
    basis = {}
    insert_column(tuple(value for row in queue[0] for value in row), basis)
    cursor = 0
    products = 0
    while cursor < len(queue) and len(basis) < 36:
        matrix = queue[cursor]
        cursor += 1
        for generator in generators:
            product = multiply(matrix, generator)
            products += 1
            if insert_column(tuple(value for row in product for value in row), basis):
                queue.append(product)
            if len(basis) == 36:
                break
    return dict(dimension=len(basis), examined_products=products,
                basis_sha256=hashlib.sha256(b"".join(bytes(x for row in matrix for x in row)
                                                      for matrix in queue)).hexdigest())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def equal(actual, expected, message):
    require(type(actual) is type(expected) and actual == expected, message)


def reject_duplicate_keys(pairs):
    value = {}
    for key, item in pairs:
        require(key not in value, "duplicate JSON key")
        value[key] = item
    return value


def field_evidence():
    return dict(polynomial=0x14d, alpha=2, alpha_order=255,
                multiplication_sha256=hashlib.sha256(b"".join(PRODUCTS)).hexdigest(),
                inverse_sha256=hashlib.sha256(bytes(INVERSES)).hexdigest())


def oracle_candidate(matrix, cumulative, crossing):
    matrix_rank = rank_columns(zip(*matrix))
    record = dict(A_rank=matrix_rank, A_sha256=matrix_hash(matrix),
                  G_sha256=None, R_sha256=None, witness=None, passed=0)
    if matrix_rank != 6:
        return record, None, None
    diagonal = tuple(tuple(scalar_power(2, i) if i == j else 0
                           for j in range(6)) for i in range(6))
    generator = multiply(multiply(matrix, diagonal), inverse(matrix))
    ratio = multiply(cumulative, generator)
    record.update(G_sha256=matrix_hash(generator), R_sha256=matrix_hash(ratio))
    columns = {i: (vandermonde_column(i + 255) if i < 0 else
                   vector_product(ratio, vandermonde_column(i))) for i in range(-9, 9)}
    for offsets in crossing:
        rank = rank_columns(tuple(columns[i] for i in offsets))
        if rank < 6:
            record["witness"] = dict(offsets=list(offsets), rank=rank)
            break
        record["passed"] += 1
    return record, generator, ratio


def oracle_final(selected, recorded):
    generators = tuple(item["G"] for item in selected)
    mapper = OracleMapper(generators)
    seen_ids, seen_columns = set(), {}
    first_collision = None
    visit_hash = hashlib.sha256()
    visits = 0

    def visit(packet_id):
        nonlocal visits, first_collision
        column = mapper.row(packet_id)
        visits += 1
        visit_hash.update(struct.pack("<I", packet_id) + bytes(column))
        seen_ids.add(packet_id)
        key = normalized(column)
        if key is None:
            if first_collision is None:
                first_collision = dict(kind="zero", ids=[packet_id], normalized=list(column))
        else:
            previous = seen_columns.setdefault(key, packet_id)
            if previous != packet_id and first_collision is None:
                first_collision = dict(kind="proportional", ids=[previous, packet_id],
                                       normalized=list(key))
        return column

    powers = []
    for stage, generator in enumerate(generators):
        checks = {}
        for exponent in (255, 85, 51, 15):
            value = power(generator, exponent)
            scalar = all(value[i][j] == (value[0][0] if i == j else 0)
                         for i in range(6) for j in range(6))
            checks[str(exponent)] = value == identity() if exponent == 255 else not scalar
        require(all(checks.values()), "independent conjugate order disagreement")
        powers.append(dict(stage=stage, checks=checks))
    algebra = algebra_evidence(generators)
    systematic = tuple(visit(i) for i in range(6))
    equal(systematic, identity(), "independent systematic identity differs")
    maximum = visit(0xffffffff)
    seams = []
    crossing = crossing_roster()
    for stage in range(4):
        boundary = 255 ** (stage + 1)
        ids = tuple(boundary + offset for offset in range(-9, 9))
        columns = {offset: visit(boundary + offset) for offset in range(-9, 9)}
        for offsets in crossing:
            require(rank_columns(tuple(columns[i] for i in offsets)) == 6,
                    "independent actual seam contradicts accepted certificate")
        seams.append(dict(stage=stage, boundary=boundary, ids=list(ids), passed=1050,
                          coefficients=[list(columns[i]) for i in range(-9, 9)]))

    require(type(recorded) is list and len(recorded) == 5533, "final roster length differs")
    cursor, passed = 0, True
    families = (("ordinary", ordinary_roster(), (6, 7, 8, 9)),
                ("modular", modular_roster(), (6, 7, 8, 10, 17)),
                ("old-companion", ((None, tuple(6 + t * 16777217 for t in range(10))),),
                 (6, 7, 8, 10)))
    for family, roster, prefixes in families:
        for index, (parameter, ids) in enumerate(roster):
            columns = tuple(visit(packet_id) for packet_id in ids)
            ranks = {str(length): rank_columns(columns[:length]) for length in prefixes}
            expected = dict(family=family, index=index, parameter=parameter, ids=ids,
                            prefix_ranks=ranks, rank=ranks[str(len(ids))])
            equal(canonical(recorded[cursor]), canonical(expected),
                  "independent final row differs at {}:{}".format(family, index))
            passed &= expected["rank"] == 6
            cursor += 1
    equal(cursor, 5533, "independent final roster accounting")
    equal(visits, 91637, "independent coefficient visit accounting")
    evidence = dict(powers=powers, algebra=algebra,
                    systematic=dict(ids=list(range(6)), coefficients=systematic),
                    maximum=dict(id=0xffffffff, coefficients=maximum), seams=seams,
                    aliases=dict(visits=visits, unique_ids=len(seen_ids),
                                 normalized_rows=len(seen_columns), first_collision=first_collision,
                                 coefficient_visit_sha256=visit_hash.hexdigest()))
    return evidence, passed and algebra["dimension"] == 36 and first_collision is None


def verify_result(report):
    """Verify recorded attempts only; never search beyond a published roster."""
    require(type(report) is dict, "result must be an object")
    equal(set(report), {"protocol", "outcome", "search", "selected", "exhausted_stage",
                        "selected_sha256", "final", "evidence", "counts"}, "result keys differ")
    equal(report["protocol"], PROTOCOL, "protocol differs")
    require(report["outcome"] in ("PASS", "FAIL", "EXHAUSTED"),
            "incomplete INVALID result cannot establish a complete mathematical replay")
    search = report["search"]
    require(type(search) is list and 0 < len(search) <= 512, "search roster size")
    require(type(report["selected"]) is list, "selected roster type")
    crossing = crossing_roster()
    selected, cursor, checks, exhausted = [], 0, 0, None
    cumulative = identity()
    for stage in range(4):
        accepted = False
        for attempt in range(128):
            require(cursor < len(search), "missing recorded candidate; no implicit retry")
            row = search[cursor]
            require(type(row) is dict, "candidate record type")
            equal(row.get("stage"), stage, "candidate stage differs")
            equal(row.get("attempt"), attempt, "candidate attempt differs")
            matrix = derive_recorded_candidate(stage, attempt)
            expected, generator, ratio = oracle_candidate(matrix, cumulative, crossing)
            accepted = expected["A_rank"] == 6 and expected["witness"] is None
            expected.update(stage=stage, attempt=attempt,
                            status="ACCEPTED" if accepted else
                            "SINGULAR" if expected["A_rank"] < 6 else "REJECTED")
            equal(canonical(row), canonical(expected), "candidate certificate differs")
            checks += expected["passed"] + int(expected["witness"] is not None)
            cursor += 1
            if accepted:
                require(expected["passed"] == 1050, "accepted certificate incomplete")
                item = dict(stage=stage, attempt=attempt, A=matrix, G=generator, R=ratio,
                            A_sha256=expected["A_sha256"], G_sha256=expected["G_sha256"],
                            R_sha256=expected["R_sha256"])
                selected.append(item)
                cumulative = ratio
                break
        if not accepted:
            exhausted = stage
            break
    equal(cursor, len(search), "extra attempts after first acceptance/exhaustion")
    equal(canonical(report["selected"]), canonical(selected), "selected matrix receipts differ")
    equal(report["selected_sha256"], hashlib.sha256(canonical(selected)).hexdigest(),
          "selected seal differs")
    equal(report["exhausted_stage"], exhausted, "exhausted stage differs")
    evidence = dict(field=field_evidence())
    if exhausted is not None:
        equal(report["final"], [], "exhausted search evaluated final roster")
        outcome = "EXHAUSTED"
    else:
        final_evidence, passed = oracle_final(selected, report["final"])
        evidence.update(final_evidence)
        outcome = "PASS" if passed else "FAIL"
    equal(canonical(report["evidence"]), canonical(evidence), "final evidence differs")
    equal(report["outcome"], outcome, "result decision differs")
    counts = dict(attempts=len(search), selected=len(selected), candidate_rank_checks=checks,
                  final_patterns=0 if exhausted is not None else 5533,
                  final_visits=0 if exhausted is not None else 91558)
    equal(canonical(report["counts"]), canonical(counts), "work counts differ")
    require(checks <= 537600, "candidate certificate rank-check cap")
    return dict(protocol=PROTOCOL, verified=True, outcome=outcome, counts=counts)


class ArithmeticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.worker = load_worker()
        cls.worker.init_field()
        cls.candidate_function = staticmethod(cls.worker.candidate_matrix)

    def setUp(self):
        # Any accidental unmocked candidate-pool use is a unit-test failure.
        for module, name in ((self.worker, "candidate_matrix"),
                             (sys.modules[__name__], "derive_recorded_candidate")):
            patch = mock.patch.object(module, name,
                                      side_effect=AssertionError("actual SHA pool forbidden in unit tests"))
            patch.start()
            self.addCleanup(patch.stop)

    def test_field_against_independent_polynomial_long_division(self):
        for a in range(256):
            for b in range(256):
                self.assertEqual(self.worker.gf_mul_raw(a, b), PRODUCTS[a][b])
        self.assertEqual(len({scalar_power(2, i) for i in range(255)}), 255)
        self.assertEqual(scalar_power(2, 255), 1)
        for divisor in (3, 5, 17):
            self.assertNotEqual(scalar_power(2, 255 // divisor), 1)
        for value in range(1, 256):
            self.assertEqual(PRODUCTS[value][INVERSES[value]], 1)

    def test_neutral_matrix_arithmetic_and_inverse(self):
        matrices = (identity(),) + neutral_generators()
        for a in matrices:
            self.assertEqual(self.worker.matrix_inverse(a), inverse(a))
            self.assertEqual(multiply(a, inverse(a)), identity())
            self.assertEqual(self.worker.matrix_rank(a), rank_columns(zip(*a)))
            for b in matrices:
                self.assertEqual(self.worker.matrix_multiply(a, b), multiply(a, b))
            for exponent in (0, 1, 15, 51, 85, 254, 255):
                self.assertEqual(self.worker.matrix_power(a, exponent), power(a, exponent))
        singular = tuple(identity()[i // 2] for i in range(6))
        self.assertIsNone(inverse(singular))
        self.assertEqual(self.worker.matrix_rank(singular), 3)

    def test_roster_counts_and_exact_order(self):
        crossing = crossing_roster()
        self.assertEqual(len(crossing), 1050)
        self.assertEqual(tuple(self.worker.crossing_patterns()), crossing)
        ordinary = tuple(ordinary_roster())
        modular = tuple(modular_roster())
        self.assertEqual(len(ordinary), 312)
        self.assertEqual(len(modular), 5220)
        self.assertEqual(tuple(self.worker.ordinary_patterns()), ordinary)
        self.assertEqual(tuple(self.worker.modular_patterns()), modular)
        self.assertEqual(sum(len(ids) for _, ids in ordinary + modular) + 10, 91558)
        self.assertEqual(len(ordinary) + len(modular) + 1, 5533)
        self.assertEqual(len({packet for _, ids in modular for packet in ids}), 17 ** 4)
        for _, ids in ordinary:
            self.assertEqual(len(set(ids)), 9)
            self.assertEqual(len({b - a for a, b in zip(ids, ids[1:])}), 1)
            self.assertTrue(0 <= ids[0] < ids[-1] <= 0xffffffff)

    def test_neutral_noncommuting_carry_and_arbitrary_lookup(self):
        generators = neutral_generators()
        self.assertNotEqual(multiply(generators[0], generators[1]),
                            multiply(generators[1], generators[0]))
        mapper = self.worker.RowMapper(generators)
        packet_ids = (0xffffffff, 0, 5, 255 ** 4, 254, 255, 255 ** 2 - 1,
                      255 ** 2, 255 ** 3 - 1, 255 ** 3, 4, 3, 2, 1)
        for packet_id in packet_ids:
            self.assertEqual(mapper.row(packet_id), independent_row(generators, packet_id))
        for j in range(6):
            self.assertEqual(mapper.row(j), identity()[j])
        cumulative = identity()
        for m in range(4):
            cumulative = multiply(cumulative, generators[m])
            before = identity()
            for i in range(m):
                before = multiply(power(generators[i], 254), before)
            relative = multiply(inverse(before), generators[m])
            self.assertEqual(relative, cumulative)

    def test_candidate_bytes_use_exact_two_mocked_hash_messages(self):
        messages = []

        class FakeHash:
            def __init__(self, message):
                messages.append(message)
                self.value = bytes([len(messages)]) * 32

            def digest(self):
                return self.value

        with mock.patch.object(self.worker.hashlib, "sha256", side_effect=FakeHash):
            matrix = self.candidate_function(3, 127)
        self.assertEqual(messages, [(PROTOCOL + ":stage/3/candidate/127/chunk/" + str(k)).encode("ascii")
                                    for k in range(2)])
        self.assertEqual(bytes(x for row in matrix for x in row), bytes([1]) * 32 + bytes([2]) * 4)

    def test_negative_inputs_and_zero_leading_projective_coordinates(self):
        for value in (-1, 256, True, 1.0):
            with self.assertRaises(ValueError):
                self.worker.gf_mul_raw(value, 1)
        with self.assertRaises(ValueError):
            self.worker.matrix_inverse(tuple((0,) * 6 for _ in range(6)))
        with self.assertRaises(ValueError):
            self.worker.matrix_power(identity(), -1)
        self.assertEqual(self.worker.algebra_dimension((identity(),)), algebra_evidence((identity(),)))
        value = (0, 0, 2, 4, 6, 8)
        scaled = tuple(PRODUCTS[7][x] for x in value)

        class FakeMapper:
            def row(self, packet_id):
                return value if packet_id == 0 else scaled if packet_id == 1 else (0,) * 6

        ledger = self.worker.AliasLedger(FakeMapper())
        ledger.visit(0)
        ledger.visit(0)
        self.assertIsNone(ledger.first_collision)
        ledger.visit(1)
        self.assertEqual(ledger.first_collision,
                         dict(kind="proportional", ids=[0, 1], normalized=list(normalized(value))))
        zero_ledger = self.worker.AliasLedger(FakeMapper())
        zero_ledger.visit(2)
        self.assertEqual(zero_ledger.first_collision["kind"], "zero")

    def test_neutral_full_matrix_algebra(self):
        diagonal = tuple(tuple(scalar_power(2, i) if i == j else 0
                               for j in range(6)) for i in range(6))
        cycle = tuple(tuple(int(j == (i + 1) % 6) for j in range(6)) for i in range(6))
        generators = (diagonal, cycle, identity(), identity())
        expected = algebra_evidence(generators)
        self.assertEqual(expected["dimension"], 36)
        self.assertEqual(self.worker.algebra_dimension(generators), expected)

    def mock_checked(self, accepted):
        return dict(A_rank=6, A_sha256=matrix_hash(identity()), G=identity(), R=identity(),
                    G_sha256=matrix_hash(identity()), R_sha256=matrix_hash(identity()),
                    witness=None if accepted else dict(offsets=[-1, 0, 1, 2, 3, 4], rank=5),
                    passed=1050 if accepted else 0)

    def test_mocked_search_first_pass_and_exhaustion(self):
        class NoBudget:
            def check(self):
                pass

        expected_attempts = [2, 1, 0, 3]
        coordinates = []

        def source(stage, attempt):
            coordinates.append((stage, attempt))
            return identity()

        def checker(matrix, cumulative, budget):
            stage, attempt = coordinates[-1]
            return self.mock_checked(attempt == expected_attempts[stage])

        search, selected, exhausted = self.worker.search_generators(source, checker, NoBudget())
        self.assertIsNone(exhausted)
        self.assertEqual([item["attempt"] for item in selected], expected_attempts)
        self.assertEqual(coordinates, [(s, a) for s in range(4) for a in range(expected_attempts[s] + 1)])
        self.assertEqual(len(search), 10)
        coordinates.clear()
        search, selected, exhausted = self.worker.search_generators(
            source, lambda *args: self.mock_checked(False), NoBudget())
        self.assertEqual(exhausted, 0)
        self.assertEqual(selected, [])
        self.assertEqual(coordinates, [(0, a) for a in range(128)])
        with mock.patch.object(self.worker, "search_generators", return_value=(search, selected, exhausted)), \
                mock.patch.object(self.worker, "evaluate_final", side_effect=AssertionError("exhausted final work")):
            result = self.worker.run_screen()
        self.assertEqual(result["outcome"], "EXHAUSTED")
        self.assertEqual(result["final"], [])

    def test_mocked_incomplete_acceptance_is_invalid(self):
        checked = self.mock_checked(True)
        checked["passed"] = 1049
        with self.assertRaises(ValueError):
            self.worker.search_generators(lambda *args: identity(), lambda *args: checked)

    def test_mocked_final_failure_still_visits_every_frozen_pattern(self):
        class NoBudget:
            def check(self):
                pass

        class FakeMapper:
            def __init__(self, *args, **kwargs):
                pass

            def row(self, packet_id):
                return identity()[packet_id % 6]

        selected = [dict(G=identity()) for _ in range(4)]
        nonscalar = neutral_generators()[0]
        rank_calls = []

        def rank_fake(columns):
            rank_calls.append(len(columns))
            # First ordinary pattern's full-nine check fails; later rows stay.
            return 5 if len(rank_calls) == 4204 else 6

        with mock.patch.object(self.worker, "matrix_power", side_effect=lambda a, n: identity() if n == 255 else nonscalar), \
                mock.patch.object(self.worker, "algebra_dimension", return_value=dict(dimension=1)), \
                mock.patch.object(self.worker, "RowMapper", FakeMapper), \
                mock.patch.object(self.worker, "rank_columns", side_effect=rank_fake):
            final, evidence, passed = self.worker.evaluate_final(selected, NoBudget())
        self.assertFalse(passed)
        self.assertEqual(final[0]["rank"], 5)
        self.assertEqual(len(final), 5533)
        self.assertEqual(sum(len(row["ids"]) for row in final), 91558)
        self.assertEqual(evidence["aliases"]["visits"], 91637)
        self.assertEqual(len(rank_calls), 31552)
        self.assertEqual(final[-1]["family"], "old-companion")
        # Exercise the complete independent final-record/alias plumbing using
        # only mocked neutral rows and ranks, never a candidate evaluation.
        current = sys.modules[__name__]
        with mock.patch.object(current, "OracleMapper", FakeMapper), \
                mock.patch.object(current, "power", side_effect=lambda a, n: identity() if n == 255 else nonscalar), \
                mock.patch.object(current, "algebra_evidence", return_value=dict(dimension=1)), \
                mock.patch.object(current, "rank_columns", side_effect=rank_fake):
            rank_calls.clear()
            replayed, replay_passed = oracle_final(selected, json.loads(canonical(final)))
            self.assertEqual(canonical(replayed), canonical(evidence))
            self.assertFalse(replay_passed)
            self.assertEqual(len(rank_calls), 31552)
            for mutation in ("missing", "boolean_id", "float_rank"):
                changed = json.loads(canonical(final))
                if mutation == "missing":
                    changed.pop()
                elif mutation == "boolean_id":
                    changed[0]["ids"][0] = False
                else:
                    changed[0]["prefix_ranks"]["6"] = 6.0
                rank_calls.clear()
                with self.assertRaises(ValueError):
                    oracle_final(selected, changed)
        search = [dict(self.mock_checked(True), stage=stage, attempt=0) for stage in range(4)]
        with mock.patch.object(self.worker, "search_generators", return_value=(search, selected, None)) as search_mock, \
                mock.patch.object(self.worker, "evaluate_final", return_value=(final, evidence, False)) as final_mock:
            result = self.worker.run_screen()
        self.assertEqual(result["outcome"], "FAIL")
        self.assertEqual(result["final"], final)
        self.assertEqual(search_mock.call_count, 1)
        self.assertEqual(final_mock.call_count, 1)

    def test_independent_mock_replay_exhaustion_and_receipt_mutations(self):
        checked = self.mock_checked(False)
        base_record = {key: value for key, value in checked.items() if key not in ("G", "R")}
        records = [dict(base_record, stage=0, attempt=i, status="REJECTED") for i in range(128)]
        report = dict(protocol=PROTOCOL, outcome="EXHAUSTED", search=records, selected=[],
                      exhausted_stage=0, selected_sha256=hashlib.sha256(canonical([])).hexdigest(),
                      final=[], evidence=dict(field=field_evidence()),
                      counts=dict(attempts=128, selected=0, candidate_rank_checks=128,
                                  final_patterns=0, final_visits=0))
        changes = (
            lambda r: r["search"][0].update(stage=True),
            lambda r: r["search"][1].update(attempt=2),
            lambda r: r["search"][0].update(passed=1),
            lambda r: r["search"][0]["witness"].update(rank=True),
            lambda r: r["search"].pop(),
            lambda r: r["search"].append(dict(r["search"][-1])),
            lambda r: r.update(selected_sha256="0" * 64),
            lambda r: r.update(exhausted_stage=False),
            lambda r: r["counts"].update(candidate_rank_checks=127),
            lambda r: r["evidence"]["field"].update(alpha_order=254),
            lambda r: r.update(final=[{}]),
        )
        current = sys.modules[__name__]
        with mock.patch.object(current, "derive_recorded_candidate", return_value=identity()) as source, \
                mock.patch.object(current, "oracle_candidate",
                                  side_effect=lambda *args: (dict(base_record), identity(), identity())):
            self.assertEqual(verify_result(report)["outcome"], "EXHAUSTED")
            self.assertEqual(source.call_args_list[:128], [mock.call(0, a) for a in range(128)])
            for change in changes:
                changed = json.loads(canonical(report))
                change(changed)
                with self.assertRaises(ValueError):
                    verify_result(changed)

    def test_independent_mock_replay_selected_seal_and_final_failure(self):
        checked = self.mock_checked(True)
        base_record = {key: value for key, value in checked.items() if key not in ("G", "R")}
        selected = [dict(stage=i, attempt=0, A=identity(), G=identity(), R=identity(),
                         A_sha256=matrix_hash(identity()), G_sha256=matrix_hash(identity()),
                         R_sha256=matrix_hash(identity())) for i in range(4)]
        report = dict(protocol=PROTOCOL, outcome="FAIL", selected=selected,
                      search=[dict(base_record, stage=i, attempt=0, status="ACCEPTED") for i in range(4)],
                      exhausted_stage=None, selected_sha256=hashlib.sha256(canonical(selected)).hexdigest(),
                      final=[{}] * 5533, evidence=dict(field=field_evidence(), algebra=dict(dimension=1)),
                      counts=dict(attempts=4, selected=4, candidate_rank_checks=4200,
                                  final_patterns=5533, final_visits=91558))
        current = sys.modules[__name__]
        with mock.patch.object(current, "derive_recorded_candidate", return_value=identity()) as source, \
                mock.patch.object(current, "oracle_candidate", side_effect=lambda *args: (dict(base_record), identity(), identity())), \
                mock.patch.object(current, "oracle_final", return_value=(dict(algebra=dict(dimension=1)), False)) as final:
            self.assertEqual(verify_result(report)["outcome"], "FAIL")
            self.assertEqual(source.call_count, 4)
            self.assertEqual(final.call_count, 1)
            changed = json.loads(canonical(report))
            changed["outcome"] = "PASS"
            with self.assertRaises(ValueError):
                verify_result(changed)
            changed = json.loads(canonical(report))
            changed["selected"][0]["A"][0][0] = True
            changed["selected_sha256"] = hashlib.sha256(canonical(changed["selected"])).hexdigest()
            with self.assertRaises(ValueError):
                verify_result(changed)


def replay_file(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(8 * 1024 * 1024 + 1)
    require(len(raw) <= 8 * 1024 * 1024, "replay byte cap")
    report = json.loads(raw, object_pairs_hook=reject_duplicate_keys,
                        parse_constant=lambda value: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    equal(raw, canonical(report) + b"\n", "noncanonical result encoding")
    return verify_result(report)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--replay":
        try:
            print(canonical(replay_file(sys.argv[2])).decode("ascii"))
        except (ValueError, OSError, TypeError, KeyError, IndexError) as error:
            print("independent radix replay: " + str(error), file=sys.stderr)
            sys.exit(1)
    else:
        unittest.main()
