#!/usr/bin/env python3
"""Frozen K6 ordinary-GF256 structural screen; not a codec or benchmark.

Importing this module never generates a candidate or evaluates the campaign.
The only campaign entry point is run_screen(), or the isolated --worker CLI.
"""
import collections
import hashlib
import itertools
import json
import resource
import struct
import sys
import time


PROTOCOL = "wirehair.wh2.noncommuting-radix-r0"
POLYNOMIAL = 0x14D
UINT32_MAX = (1 << 32) - 1
STDOUT_LIMIT = 4 * 1024 * 1024
MUL = None
INV = None
EXP = None
FIELD_EVIDENCE = None


class ScreenInvalid(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ScreenInvalid(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("ascii")


def digest(value):
    return hashlib.sha256(value).hexdigest()


def gf_mul_raw(x, y):
    """Polynomial shift/reduction oracle used to construct scalar tables."""
    require(type(x) is int and type(y) is int and 0 <= x < 256 and 0 <= y < 256,
            "GF256 byte range")
    result = 0
    while y:
        if y & 1:
            result ^= x
        x <<= 1
        if x & 256:
            x ^= POLYNOMIAL
        y >>= 1
    return result


def init_field():
    global MUL, INV, EXP, FIELD_EVIDENCE
    if FIELD_EVIDENCE is not None:
        return dict(FIELD_EVIDENCE)
    # Establish the full alpha cycle using polynomial arithmetic, not log tables.
    exp = [1]
    for unused in range(254):
        exp.append(gf_mul_raw(exp[-1], 2))
    require(len(set(exp)) == 255 and set(exp) == set(range(1, 256)) and
            gf_mul_raw(exp[-1], 2) == 1, "alpha does not have order 255")
    mul = tuple(bytes(gf_mul_raw(x, y) for y in range(256)) for x in range(256))
    inverse = [0] * 256
    for exponent, value in enumerate(exp):
        inverse[value] = exp[(-exponent) % 255]
        require(mul[value][inverse[value]] == 1, "GF256 inverse disagreement")
    MUL, INV, EXP = mul, tuple(inverse), tuple(exp)
    FIELD_EVIDENCE = dict(polynomial=POLYNOMIAL, alpha=2, alpha_order=255,
                          multiplication_sha256=digest(b"".join(mul)),
                          inverse_sha256=digest(bytes(inverse)))
    return dict(FIELD_EVIDENCE)


def matrix_identity():
    return tuple(tuple(int(r == c) for c in range(6)) for r in range(6))


def matrix_bytes(matrix):
    require(len(matrix) == 6 and all(len(row) == 6 for row in matrix), "matrix shape")
    require(all(type(x) is int and 0 <= x < 256 for row in matrix for x in row),
            "matrix element")
    return bytes(x for row in matrix for x in row)


def matrix_multiply(left, right):
    init_field()
    columns = tuple(zip(*right))
    return tuple(tuple(_dot(row, column) for column in columns) for row in left)


def _dot(left, right):
    value = 0
    for x, y in zip(left, right):
        value ^= MUL[x][y]
    return value


def matrix_vector(matrix, vector):
    init_field()
    return tuple(_dot(row, vector) for row in matrix)


def matrix_rank(rows):
    """Forward row elimination, with ascending column and row pivot order."""
    init_field()
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    width = len(rows[0])
    require(all(len(row) == width for row in rows), "rank rectangle")
    pivot_row = 0
    for column in range(width):
        pivot = next((r for r in range(pivot_row, len(rows)) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        source = rows[pivot_row]
        scale = MUL[INV[source[column]]]
        for c in range(column, width):
            source[c] = scale[source[c]]
        for r in range(pivot_row + 1, len(rows)):
            if rows[r][column]:
                table = MUL[rows[r][column]]
                for c in range(column, width):
                    rows[r][c] ^= table[source[c]]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def rank_columns(columns):
    return matrix_rank(tuple(zip(*columns))) if columns else 0


def matrix_inverse(matrix):
    init_field()
    matrix_bytes(matrix)
    rows = [list(row) + list(unit) for row, unit in zip(matrix, matrix_identity())]
    for column in range(6):
        pivot = next((r for r in range(column, 6) if rows[r][column]), None)
        require(pivot is not None, "singular inverse")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        source = rows[column]
        table = MUL[INV[source[column]]]
        rows[column] = source = [table[value] for value in source]
        for r in range(6):
            if r != column and rows[r][column]:
                table = MUL[rows[r][column]]
                rows[r] = [x ^ table[y] for x, y in zip(rows[r], source)]
    return tuple(tuple(row[6:]) for row in rows)


def matrix_power(matrix, exponent):
    require(type(exponent) is int and exponent >= 0, "matrix exponent")
    result = matrix_identity()
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, matrix)
        matrix = matrix_multiply(matrix, matrix)
        exponent >>= 1
    return result


def is_scalar(matrix):
    value = matrix[0][0]
    return all(matrix[r][c] == (value if r == c else 0)
               for r in range(6) for c in range(6))


def vandermonde(phase):
    init_field()
    require(type(phase) is int and 0 <= phase < 255, "phase range")
    return tuple(EXP[(r * phase) % 255] for r in range(6))


def systematic_basis():
    return matrix_inverse(tuple(zip(*(vandermonde(j) for j in range(6)))))


def candidate_matrix(stage, attempt):
    require(type(stage) is int and 0 <= stage < 4 and type(attempt) is int and
            0 <= attempt < 128, "candidate coordinates")
    stem = PROTOCOL + ":stage/" + str(stage) + "/candidate/" + str(attempt) + "/chunk/"
    raw = b"".join(hashlib.sha256((stem + str(k)).encode("ascii")).digest()
                   for k in range(2))[:36]
    return tuple(tuple(raw[r * 6:r * 6 + 6]) for r in range(6))


def crossing_patterns():
    for s in range(1, 10):
        for rest in itertools.combinations(range(-s + 1, -s + 10), 5):
            if rest[-1] >= 0:
                yield (-s,) + rest


class Budget:
    def __init__(self, seconds=60):
        self.deadline = time.monotonic() + seconds

    def check(self):
        require(time.monotonic() < self.deadline, "worker deadline")


def check_conjugator(matrix, cumulative, budget):
    rank = matrix_rank(matrix)
    result = dict(A_rank=rank, A_sha256=digest(matrix_bytes(matrix)), G=None, R=None,
                  G_sha256=None, R_sha256=None, witness=None, passed=0)
    if rank != 6:
        return result
    diagonal = tuple(tuple(EXP[r] if r == c else 0 for c in range(6)) for r in range(6))
    generator = matrix_multiply(matrix_multiply(matrix, diagonal), matrix_inverse(matrix))
    ratio = matrix_multiply(cumulative, generator)
    result.update(G=generator, R=ratio, G_sha256=digest(matrix_bytes(generator)),
                  R_sha256=digest(matrix_bytes(ratio)))
    columns = {offset: (vandermonde(255 + offset) if offset < 0 else
                        matrix_vector(ratio, vandermonde(offset))) for offset in range(-9, 9)}
    for offsets in crossing_patterns():
        if result["passed"] % 32 == 0:
            budget.check()
        rank = rank_columns(tuple(columns[x] for x in offsets))
        if rank != 6:
            result["witness"] = dict(offsets=list(offsets), rank=rank)
            return result
        result["passed"] += 1
    require(result["passed"] == 1050, "crossing roster size")
    return result


def search_generators(candidate_source=None, checker=None, budget=None):
    """Injectable only for neutral unit tests; campaign uses the fixed defaults."""
    candidate_source = candidate_matrix if candidate_source is None else candidate_source
    checker = check_conjugator if checker is None else checker
    budget = Budget() if budget is None else budget
    search, selected = [], []
    cumulative = matrix_identity()
    for stage in range(4):
        for attempt in range(128):
            budget.check()
            matrix = candidate_source(stage, attempt)
            checked = checker(matrix, cumulative, budget)
            accepted = checked["A_rank"] == 6 and checked["witness"] is None
            require(not accepted or checked["passed"] == 1050, "incomplete candidate certificate")
            record = {key: value for key, value in checked.items() if key not in ("G", "R")}
            record.update(stage=stage, attempt=attempt,
                          status="ACCEPTED" if accepted else
                          ("SINGULAR" if checked["A_rank"] < 6 else "REJECTED"))
            search.append(record)
            if accepted:
                selected.append(dict(stage=stage, attempt=attempt, A=matrix, G=checked["G"],
                                     R=checked["R"], A_sha256=checked["A_sha256"],
                                     G_sha256=checked["G_sha256"], R_sha256=checked["R_sha256"]))
                cumulative = checked["R"]
                break
        else:
            return search, selected, stage
    return search, selected, None


def ordinary_patterns():
    weights = tuple(255 ** r for r in range(4))
    for z in itertools.product(range(-2, 3), repeat=4):
        steps = tuple(15 * x for x in z)
        if sum(x * weight for x, weight in zip(steps, weights)) <= 0:
            continue
        starts = tuple(-8 * min(x, 0) for x in steps)
        yield z, tuple(sum((start + t * step) * weight for start, step, weight in
                           zip(starts, steps, weights)) for t in range(9))


def modular_patterns():
    weights = tuple(255 ** r for r in range(4))
    for z in itertools.product(range(17), repeat=4):
        first = next((x for x in z if x), None)
        if first != 1:
            continue
        yield z, tuple(sum(15 * ((t * x) % 17) * weight for x, weight in
                           zip(z, weights)) for t in range(17))


def final_patterns():
    for family, iterator, prefixes in (("ordinary", ordinary_patterns(), (6, 7, 8, 9)),
                                        ("modular", modular_patterns(), (6, 7, 8, 10, 17))):
        for index, (parameter, ids) in enumerate(iterator):
            yield family, index, parameter, ids, prefixes
    yield "old-companion", 0, None, tuple(6 + t * 16777217 for t in range(10)), (6, 7, 8, 10)


class RowMapper:
    def __init__(self, generators, S=None, budget=None):
        require(len(generators) == 4, "generator count")
        self.S = systematic_basis() if S is None else S
        self.powers = []
        for generator in generators:
            matrix_bytes(generator)
            powers = [matrix_identity()]
            for unused in range(254):
                if budget is not None:
                    budget.check()
                powers.append(matrix_multiply(powers[-1], generator))
            self.powers.append(tuple(powers))
        self.final_powers = tuple(matrix_multiply(self.S, power) for power in self.powers[3])
        self.phases = tuple(vandermonde(j) for j in range(255))
        self.cache = {j: matrix_vector(self.S, self.phases[j]) for j in range(255)}

    def row(self, packet_id):
        require(type(packet_id) is int and 0 <= packet_id <= UINT32_MAX, "packet ID range")
        if packet_id not in self.cache:
            phase, block = packet_id % 255, packet_id // 255
            value = self.phases[phase]
            for r in range(3):
                digit, block = block % 255, block // 255
                if digit:
                    value = matrix_vector(self.powers[r][digit], value)
            require(block <= 1, "highest packet digit")
            value = matrix_vector(self.final_powers[block], value)
            self.cache[packet_id] = value
        return self.cache[packet_id]


class AliasLedger:
    def __init__(self, mapper):
        self.mapper = mapper
        self.ids = set()
        self.normalized = {}
        self.first_collision = None
        self.visits = 0
        self.sha = hashlib.sha256()

    def visit(self, packet_id):
        value = self.mapper.row(packet_id)
        self.visits += 1
        self.sha.update(struct.pack("<I", packet_id) + bytes(value))
        self.ids.add(packet_id)
        first = next((x for x in value if x), 0)
        if not first:
            if self.first_collision is None:
                self.first_collision = dict(kind="zero", ids=[packet_id], normalized=list(value))
            return value
        key = tuple(MUL[INV[first]][x] for x in value)
        previous = self.normalized.get(key)
        if previous is not None and previous != packet_id and self.first_collision is None:
            self.first_collision = dict(kind="proportional", ids=[previous, packet_id], normalized=list(key))
        self.normalized.setdefault(key, packet_id)
        return value

    def evidence(self):
        return dict(visits=self.visits, unique_ids=len(self.ids),
                    normalized_rows=len(self.normalized), first_collision=self.first_collision,
                    coefficient_visit_sha256=self.sha.hexdigest())


def algebra_dimension(generators, budget=None):
    """BFS on independent original matrices; elimination uses row-major pivots."""
    pivots = {}
    basis = []
    queue = collections.deque()
    products = 0

    def insert(matrix):
        vector = list(matrix_bytes(matrix))
        for position in range(36):
            if not vector[position]:
                continue
            if position in pivots:
                table = MUL[vector[position]]
                vector = [x ^ table[y] for x, y in zip(vector, pivots[position])]
            else:
                table = MUL[INV[vector[position]]]
                pivots[position] = tuple(table[x] for x in vector)
                basis.append(matrix)
                queue.append(matrix)
                return

    init_field()
    insert(matrix_identity())
    while queue and len(basis) < 36:
        if budget is not None:
            budget.check()
        parent = queue.popleft()
        for generator in generators:
            insert(matrix_multiply(parent, generator))
            products += 1
            if len(basis) == 36:
                break
    return dict(dimension=len(basis), examined_products=products,
                basis_sha256=digest(b"".join(matrix_bytes(matrix) for matrix in basis)))


def evaluate_final(selected, budget):
    generators = tuple(item["G"] for item in selected)
    powers = []
    for stage, generator in enumerate(generators):
        checks = {str(exponent): (matrix_power(generator, exponent) == matrix_identity()
                                 if exponent == 255 else not is_scalar(matrix_power(generator, exponent)))
                  for exponent in (255, 85, 51, 15)}
        powers.append(dict(stage=stage, checks=checks))
    require(all(all(item["checks"].values()) for item in powers), "conjugate power disagreement")
    algebra = algebra_dimension(generators, budget)
    mapper = RowMapper(generators, budget=budget)
    ledger = AliasLedger(mapper)
    systematic = [ledger.visit(j) for j in range(6)]
    require(tuple(systematic) == matrix_identity(), "identity systematic disagreement")
    maximum = ledger.visit(UINT32_MAX)
    seams = []
    for stage in range(4):
        boundary = 255 ** (stage + 1)
        ids = tuple(boundary + offset for offset in range(-9, 9))
        columns = {offset: ledger.visit(boundary + offset) for offset in range(-9, 9)}
        passed = 0
        for offsets in crossing_patterns():
            if passed % 32 == 0:
                budget.check()
            require(rank_columns(tuple(columns[x] for x in offsets)) == 6,
                    "actual seam disagrees with accepted relative certificate")
            passed += 1
        require(passed == 1050, "actual seam roster size")
        seams.append(dict(stage=stage, boundary=boundary, ids=list(ids), passed=passed,
                          coefficients=[list(columns[offset]) for offset in range(-9, 9)]))
    final = []
    visits = 0
    for family, index, parameter, ids, prefixes in final_patterns():
        budget.check()
        require(len(set(ids)) == len(ids), "duplicate pattern ID")
        columns = tuple(ledger.visit(packet_id) for packet_id in ids)
        visits += len(ids)
        ranks = {str(length): rank_columns(columns[:length]) for length in prefixes}
        final.append(dict(family=family, index=index, parameter=parameter, ids=ids,
                          prefix_ranks=ranks, rank=ranks[str(len(ids))]))
    require(len(final) == 5533 and visits == 91558, "final roster accounting")
    evidence = dict(powers=powers, algebra=algebra,
                    systematic=dict(ids=list(range(6)), coefficients=systematic),
                    maximum=dict(id=UINT32_MAX, coefficients=maximum), seams=seams,
                    aliases=ledger.evidence())
    passed = algebra["dimension"] == 36 and ledger.first_collision is None and all(
        item["rank"] == 6 for item in final)
    return final, evidence, passed


def run_screen():
    budget = Budget()
    result = dict(protocol=PROTOCOL, outcome="INVALID", search=[], selected=[],
                  exhausted_stage=None, final=[], evidence={}, counts={})
    try:
        result["evidence"]["field"] = init_field()
        search, selected, exhausted = search_generators(budget=budget)
        result.update(search=search, selected=selected, exhausted_stage=exhausted)
        result["selected_sha256"] = digest(canonical(selected))
        if exhausted is not None:
            result["outcome"] = "EXHAUSTED"
        else:
            final, evidence, passed = evaluate_final(selected, budget)
            result["final"] = final
            result["evidence"].update(evidence)
            result["outcome"] = "PASS" if passed else "FAIL"
        result["counts"] = dict(attempts=len(search), selected=len(selected),
                                candidate_rank_checks=sum(item["passed"] + int(item["witness"] is not None)
                                                          for item in search),
                                final_patterns=len(result["final"]),
                                final_visits=sum(len(item["ids"]) for item in result["final"]))
        require(result["counts"]["attempts"] <= 512 and
                result["counts"]["candidate_rank_checks"] <= 537600, "search accounting cap")
        budget.check()
    except ScreenInvalid as error:
        result["outcome"] = "INVALID"
        result["error"] = str(error)[:1024]
    return result


def main(argv):
    if argv != ["--worker"]:
        sys.stderr.write("usage: Wh2NoncommutingRadixR0.py --worker\n")
        return 2
    try:
        resource.setrlimit(resource.RLIMIT_AS, (512 * 1024 * 1024, 512 * 1024 * 1024))
        result = run_screen()
        raw = canonical(result) + b"\n"
        require(len(raw) <= STDOUT_LIMIT, "worker JSON byte cap")
        sys.stdout.buffer.write(raw)
        sys.stdout.buffer.flush()
        return 1 if result["outcome"] == "INVALID" else 0
    except Exception as error:
        sys.stderr.write(("noncommuting-radix-r0: " + type(error).__name__ + ": " + str(error))[:1024] + "\n")
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
