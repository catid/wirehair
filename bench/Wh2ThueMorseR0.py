#!/usr/bin/env python3
"""Frozen .62 constructive GF256 structural screen; not a codec or benchmark.

Imports are inert. Only run_screen()/--worker constructs the frozen geometric
feedback and chooses its sole parameter. Historical campaigns are never run.
"""

import importlib.util
import itertools
from pathlib import Path
import resource
import sys


_SPEC = importlib.util.spec_from_file_location(
    "_wh2_thue_radix_arithmetic", str(Path(__file__).resolve().with_name(
        "Wh2NoncommutingRadixR0.py")))
RADIX = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(RADIX)

PROTOCOL = "wirehair.wh2.thue-morse-r0"
WORDS = ("0010", "0011", "0100", "0101", "0110",
         "1001", "1010", "1011", "1100", "1101")
UINT32_MAX = (1 << 32) - 1
require = RADIX.require
canonical = RADIX.canonical


def _field():
    if RADIX.MUL is None:
        RADIX.init_field()


def _bytes(values):
    values = tuple(values)
    require(all(type(x) is int and 0 <= x < 256 for x in values), "GF256 byte range")
    return values


def _trim(values):
    values = list(values)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values) if values else (0,)


def polynomial_add(left, right):
    left, right = _bytes(left), _bytes(right)
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] ^= value
    for index, value in enumerate(right):
        result[index] ^= value
    return _trim(result)


def polynomial_multiply(left, right):
    _field()
    left, right = _bytes(left), _bytes(right)
    if not left or not right:
        return (0,)
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] ^= RADIX.MUL[x][y]
    return _trim(result)


def polynomial_evaluate(coefficients, value):
    _field()
    coefficients = _bytes(coefficients)
    require(type(value) is int and 0 <= value < 256, "polynomial argument")
    result = 0
    for coefficient in reversed(coefficients):
        result = RADIX.MUL[result][value] ^ coefficient
    return result


def polynomial_roots(coefficients):
    return [value for value in range(256)
            if polynomial_evaluate(coefficients, value) == 0]


def determinant(rows):
    """Scalar forward elimination, including det(empty)=1."""
    _field()
    rows = [list(_bytes(row)) for row in rows]
    size = len(rows)
    require(size <= 6 and all(len(row) == size for row in rows), "determinant shape")
    result = 1
    for column in range(size):
        pivot = next((r for r in range(column, size) if rows[r][column]), None)
        if pivot is None:
            return 0
        rows[column], rows[pivot] = rows[pivot], rows[column]
        value = rows[column][column]
        result = RADIX.MUL[result][value]
        inverse = RADIX.INV[value]
        for r in range(column + 1, size):
            factor = RADIX.MUL[rows[r][column]][inverse]
            for c in range(column + 1, size):
                rows[r][c] ^= RADIX.MUL[factor][rows[column][c]]
            rows[r][column] = 0
    return result


def companion(feedback):
    feedback = _bytes(feedback)
    require(len(feedback) == 6, "feedback shape")
    return tuple(tuple(feedback[r] if c == 5 else int(r == c + 1)
                       for c in range(6)) for r in range(6))


def constraint_specs():
    for m in (2, 3, 4):
        suffixes = {word[-(m - 1):] for word in WORDS}
        for columns in itertools.combinations(range(1, m + 5), m):
            active_sets = {tuple(t for t in columns if t < m and suffix[t - 1] == "1")
                           for suffix in suffixes}
            for active in sorted(active_sets):
                if active:
                    yield m, columns, active


def constraint_polynomial(feedback, m, columns, active):
    """Laplace expansion in the selected diagonal lambda entries, in char 2."""
    feedback = _bytes(feedback)
    columns, active = tuple(columns), tuple(active)
    require(len(feedback) == 6 and type(m) is int and m in (2, 3, 4), "constraint shape")
    require(len(columns) == m and tuple(sorted(set(columns))) == columns and
            all(type(c) is int and 1 <= c <= m + 4 for c in columns), "constraint columns")
    require(tuple(sorted(set(active))) == active and len(active) <= 2 and
            all(type(t) is int and t in columns and 1 <= t < m for t in active),
            "constraint active diagonal")
    coefficients = feedback + (1,)
    h = tuple(tuple(coefficients[c - r] if 0 <= c - r <= 6 else 0
                    for c in columns) for r in range(m))
    result = [0] * (len(active) + 1)
    for degree in range(len(active) + 1):
        for removed in itertools.combinations(active, degree):
            submatrix = tuple(tuple(h[r][j] for j, c in enumerate(columns) if c not in removed)
                              for r in range(m) if r not in removed)
            result[degree] ^= determinant(submatrix)
    return _trim(result)


def fixed_feedback():
    """Campaign-only: do not call from pre-campaign neutral tests."""
    _field()
    polynomial = (1,)
    for r in range(6):
        polynomial = polynomial_multiply(polynomial, (RADIX.EXP[r], 1))
    require(len(polynomial) == 7 and polynomial[-1] == 1, "geometric polynomial degree")
    return polynomial[:-1]


def select_parameter(feedback, budget=None, record=None):
    """Deterministic root exclusion; injectable feedback supports neutral tests."""
    feedback = _bytes(feedback)
    require(len(feedback) == 6 and feedback[0] != 0, "invertible baseline feedback")
    result = {} if record is None else record
    result.update(feedback=list(feedback), polynomials=[])
    excluded = {0, feedback[0]}
    degree_sum = 0
    for m, columns, active in constraint_specs():
        if budget is not None:
            budget.check()
        polynomial = constraint_polynomial(feedback, m, columns, active)
        item = dict(m=m, columns=list(columns), active=list(active), degree_bound=len(active),
                    coefficients=list(polynomial), roots=[])
        result["polynomials"].append(item)
        require(polynomial[0] != 0, "zero determinant constant: " + str((m, columns, active)))
        require(len(polynomial) - 1 <= len(active) <= 2, "determinant degree contradiction")
        item["roots"] = polynomial_roots(polynomial)
        require(len(item["roots"]) <= len(active), "polynomial root bound contradiction")
        excluded.update(item["roots"])
        degree_sum += len(active)
    require(len(result["polynomials"]) == 190 and degree_sum == 240, "polynomial roster accounting")
    eligible = sorted(set(range(256)) - excluded)
    require(len(eligible) >= 14, "constructive existence bound contradicted")
    chosen = eligible[0]
    other = (feedback[0] ^ chosen,) + feedback[1:]
    a0, a1 = companion(feedback), companion(other)
    result.update(polynomial_count=190, degree_bound_sum=degree_sum,
                  excluded=sorted(excluded), eligible=eligible, **{"lambda": chosen},
                  A0=a0, A1=a1,
                  pair_sha256=RADIX.digest(RADIX.matrix_bytes(a0) + RADIX.matrix_bytes(a1)))
    return result


def local_columns(pair, word):
    require(word in WORDS and len(pair) == 2, "local word/pair")
    columns = list(RADIX.matrix_identity())
    product = RADIX.matrix_identity()
    for bit in word:
        matrix = pair[int(bit)]
        feedback = tuple(row[5] for row in matrix)
        columns.append(RADIX.matrix_vector(product, feedback))
        product = RADIX.matrix_multiply(product, matrix)
    return tuple(columns)


def pair_evidence(pair, feedback, chosen):
    feedback = _bytes(feedback)
    require(len(feedback) == 6 and type(chosen) is int and
            0 < chosen < 256 and chosen != feedback[0], "companion parameter range")
    a0, a1 = pair
    require(a0 == companion(feedback) and
            a1 == companion((feedback[0] ^ chosen,) + tuple(feedback[1:])),
            "companion reconstruction disagreement")
    identity = RADIX.matrix_identity()
    for matrix in pair:
        inverse = RADIX.matrix_inverse(matrix)
        require(RADIX.matrix_multiply(matrix, inverse) == identity and
                RADIX.matrix_multiply(inverse, matrix) == identity, "companion inverse disagreement")
    noncommuting = RADIX.matrix_multiply(a0, a1) != RADIX.matrix_multiply(a1, a0)
    require(a0 != a1 and noncommuting, "distinct/noncommuting companion contradiction")
    difference = tuple(tuple(x ^ y for x, y in zip(left, right)) for left, right in zip(a0, a1))
    rows = []
    product = difference
    for unused in range(6):
        rows.append(product[0])
        product = RADIX.matrix_multiply(product, a0)
    rank = RADIX.matrix_rank(rows)
    require(rank == 6, "cyclic difference algebra proof contradicted")
    return dict(invertible=[True, True], noncommuting=True,
                cyclic_difference_rows=rows, cyclic_difference_rank=rank)


def _parity(value):
    return bin(value).count("1") & 1


class RowMapper:
    """Packed shared-lookup model; diagnostic block/cache storage is not a profile claim."""
    def __init__(self, pair, budget=None):
        require(len(pair) == 2, "companion pair count")
        for matrix in pair:
            RADIX.matrix_bytes(matrix)
        self.blocks = [[pair[0]], [pair[1]]]
        for level in range(31):
            if budget is not None:
                budget.check()
            left, right = self.blocks[0][level], self.blocks[1][level]
            self.blocks[0].append(RADIX.matrix_multiply(left, right))
            self.blocks[1].append(RADIX.matrix_multiply(right, left))
        descriptions, payloads = [], []

        def table(name, start_bit, width, phase, vectors=False):
            product = RADIX.matrix_identity()
            data = bytearray()
            for value in range(1 << width):
                if budget is not None and value % 32 == 0:
                    budget.check()
                data.extend(bytes(row[0] for row in product) if vectors else RADIX.matrix_bytes(product))
                if value + 1 < (1 << width):
                    product = RADIX.matrix_multiply(
                        product, self.blocks[phase ^ _parity(value)][start_bit])
            data = bytes(data)
            payloads.append(data)
            descriptions.append(dict(name=name, start_bit=start_bit, width=width, phase=phase,
                                     bytes=len(data), sha256=RADIX.digest(data)))
            return data

        self.low = tuple(table("low10", 0, 10, phase, True) for phase in range(2))
        self.middle10 = tuple(table("middle10", 10, 7, phase) for phase in range(2))
        self.middle17 = tuple(table("middle17", 17, 7, phase) for phase in range(2))
        self.high = table("high24", 24, 8, 0)
        payload = b"".join(payloads)
        require(len(payload) == 39936, "packed lookup payload accounting")
        self._evidence = dict(packed_bytes=len(payload), max_matvecs=3, max_gf_products=108,
                              max_xors=90, tables=descriptions, tables_sha256=RADIX.digest(payload),
                              dyadic_blocks_sha256=RADIX.digest(b"".join(
                                  RADIX.matrix_bytes(matrix) for phase in self.blocks for matrix in phase)))
        self.cache = {}

    @staticmethod
    def _apply(table, index, vector):
        offset = index * 36
        output = []
        for r in range(6):
            value = RADIX.MUL[table[offset + r * 6]][vector[0]]
            for c in range(1, 6):
                value ^= RADIX.MUL[table[offset + r * 6 + c]][vector[c]]
            output.append(value)
        return tuple(output)

    def row(self, packet_id):
        require(type(packet_id) is int and 0 <= packet_id <= UINT32_MAX, "packet ID range")
        if packet_id not in self.cache:
            high, mid17, mid10, low = packet_id >> 24, (packet_id >> 17) & 127, (packet_id >> 10) & 127, packet_id & 1023
            phase17 = _parity(high)
            phase10 = phase17 ^ _parity(mid17)
            phase0 = phase10 ^ _parity(mid10)
            value = self.low[phase0][low * 6:low * 6 + 6]
            value = self._apply(self.middle10[phase10], mid10, value)
            value = self._apply(self.middle17[phase17], mid17, value)
            self.cache[packet_id] = self._apply(self.high, high, value)
        return self.cache[packet_id]

    def reference_row(self, packet_id):
        """Separate bit decomposition for mandatory-ID lookup concordance."""
        require(type(packet_id) is int and 0 <= packet_id <= UINT32_MAX, "packet ID range")
        value = (1, 0, 0, 0, 0, 0)
        for level in range(32):
            if (packet_id >> level) & 1:
                value = RADIX.matrix_vector(self.blocks[_parity(packet_id >> (level + 1))][level], value)
        return value

    def evidence(self):
        return self._evidence


def final_patterns():
    yield from RADIX.final_patterns()
    for exponent in range(29):
        stride = 1 << exponent
        for start_index, start in enumerate((0, 6, UINT32_MAX - 9 * stride)):
            ids = tuple(start + j * stride for j in range(10))
            yield "dyadic", 3 * exponent + start_index, dict(exponent=exponent, start_index=start_index), ids, (6, 7, 8, 10)


def evaluate_final(pair, budget):
    mapper = RowMapper(pair, budget)
    ledger = RADIX.AliasLedger(mapper)

    def checked_visit(packet_id):
        value = ledger.visit(packet_id)
        require(value == mapper.reference_row(packet_id), "packed/bitwise row disagreement")
        return value

    systematic = [checked_visit(j) for j in range(6)]
    require(tuple(systematic) == RADIX.matrix_identity(), "systematic row disagreement")
    maximum = checked_visit(UINT32_MAX)
    seams = []
    for exponent in (10, 17, 24, 31):
        start = (1 << exponent) - 5
        ids = tuple(start + j for j in range(10))
        columns = tuple(checked_visit(packet_id) for packet_id in ids)
        passed = 0
        for selected in itertools.combinations(range(10), 6):
            if passed % 32 == 0:
                budget.check()
            require(RADIX.rank_columns(tuple(columns[j] for j in selected)) == 6,
                    "actual seam local certificate disagreement")
            passed += 1
        require(passed == 210, "seam rank accounting")
        seams.append(dict(exponent=exponent, start=start, ids=ids, coefficients=columns, passed=passed))
    final, visits = [], 0
    for family, index, parameter, ids, prefixes in final_patterns():
        budget.check()
        require(len(set(ids)) == len(ids), "duplicate diagnostic ID")
        columns = tuple(ledger.visit(packet_id) for packet_id in ids)
        visits += len(ids)
        ranks = {str(length): RADIX.rank_columns(columns[:length]) for length in prefixes}
        final.append(dict(family=family, index=index, parameter=parameter, ids=ids,
                          prefix_ranks=ranks, rank=ranks[str(len(ids))]))
    require(len(final) == 5620 and visits == 92428 and ledger.visits == 92475,
            "final pattern/visit accounting")
    evidence = dict(lookup=mapper.evidence(), systematic=dict(ids=list(range(6)), coefficients=systematic),
                    maximum=dict(id=UINT32_MAX, coefficients=maximum), seams=seams, aliases=ledger.evidence())
    passed = ledger.first_collision is None and all(item["rank"] == 6 for item in final)
    return final, evidence, passed


def run_screen():
    budget = RADIX.Budget()
    result = dict(protocol=PROTOCOL, outcome="INVALID", selection=None, selection_sha256=None,
                  local=[], final=[], evidence={}, counts={})
    try:
        result["evidence"]["field"] = RADIX.init_field()
        feedback = fixed_feedback()
        result["selection"] = {}
        selection = select_parameter(feedback, budget, result["selection"])
        result["selection_sha256"] = RADIX.digest(canonical(selection))
        pair = (selection["A0"], selection["A1"])
        result["evidence"]["pair"] = pair_evidence(pair, feedback, selection["lambda"])
        for word in WORDS:
            columns = local_columns(pair, word)
            passed = 0
            for selected in itertools.combinations(range(10), 6):
                if passed % 32 == 0:
                    budget.check()
                require(RADIX.rank_columns(tuple(columns[j] for j in selected)) == 6,
                        "constructive local certificate disagreement: " + word + ":" + str(selected))
                passed += 1
            require(passed == 210, "local rank accounting")
            result["local"].append(dict(word=word, passed=passed))
        final, evidence, passed = evaluate_final(pair, budget)
        result["final"] = final
        result["evidence"].update(evidence)
        result["counts"] = dict(polynomial_count=selection["polynomial_count"],
                                degree_bound_sum=selection["degree_bound_sum"], eligible_count=len(selection["eligible"]),
                                local_rank_checks=sum(item["passed"] for item in result["local"]),
                                seam_rank_checks=sum(item["passed"] for item in evidence["seams"]),
                                final_patterns=len(final), final_visits=sum(len(item["ids"]) for item in final),
                                total_visits=evidence["aliases"]["visits"])
        require(result["counts"]["local_rank_checks"] == 2100 and
                result["counts"]["seam_rank_checks"] == 840, "local/seam accounting")
        require(RADIX.digest(canonical(selection)) == result["selection_sha256"], "selected evidence changed")
        budget.check()
        result["outcome"] = "PASS" if passed else "FAIL"
    except RADIX.ScreenInvalid as error:
        result["outcome"] = "INVALID"
        result["error"] = str(error)[:1024]
    return result


def main(argv):
    if argv != ["--worker"]:
        sys.stderr.write("usage: Wh2ThueMorseR0.py --worker\n")
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
        sys.stderr.write(("thue-morse-r0: " + type(error).__name__ + ": " + str(error))[:1024] + "\n")
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
