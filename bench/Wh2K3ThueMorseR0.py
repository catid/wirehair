#!/usr/bin/env python3
"""One-shot .82 K3 structural/recovery screen, not a codec or speed test.

Import is inert. Only --worker selects the frozen pair or scores its domains.
The scientific protocol is preregistered in wirehair-sxvz.16.1.20.82.
"""
import hashlib
import importlib.util
import itertools
from pathlib import Path
import resource
import struct
import sys


HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


F = sibling("_k3_tm_field", "Wh2NoncommutingRadixR0.py")
H = sibling("_k3_tm_history_reader", "Wh2ThueMorseRecoveryHistoryR0.py")
require, canonical = F.require, F.canonical
PROTOCOL = "wirehair.wh2.k3-thue-morse-r0"
WORDS = ("0010", "0011", "0100", "0101", "0110",
         "1001", "1010", "1011", "1100", "1101")
WIDTHS = (2, 64, 1280)
SCHEDULES = ("iid", "burst", "adversarial", "repair-only")
TRIPLES = tuple(itertools.combinations(range(7), 3))
MAX_ID = (1 << 32) - 1
MASK64 = (1 << 64) - 1
DEVELOPMENT_TRACE_SHA = "7cb9e6b5e6951eb84e65fa75faa279e927be1c455c6ce5f2cfd458f94aef9d93"


def identity(n):
    return tuple(tuple(int(r == c) for c in range(n)) for r in range(n))


def companion(feedback):
    require(len(feedback) >= 2 and all(type(v) is int and 0 <= v < 256 for v in feedback),
            "companion feedback")
    n = len(feedback)
    return tuple(tuple(feedback[r] if c == n - 1 else int(r == c + 1)
                       for c in range(n)) for r in range(n))


def multiply_polynomial(x, y):
    """Independent carryless product followed by long division, no GF tables."""
    product = 0
    for bit in range(8):
        if y & (1 << bit):
            product ^= x << bit
    for bit in range(14, 7, -1):
        if product & (1 << bit):
            product ^= 0x14d << (bit - 8)
    return product


def determinant3(rows, mul=None):
    require(len(rows) == 3 and all(len(row) == 3 for row in rows), "3x3 determinant")
    if mul is None:
        F.init_field()
        mul = lambda x, y: F.MUL[x][y]
    a, b, c = rows
    return (mul(a[0], mul(b[1], c[2]) ^ mul(b[2], c[1])) ^
            mul(a[1], mul(b[0], c[2]) ^ mul(b[2], c[0])) ^
            mul(a[2], mul(b[0], c[1]) ^ mul(b[1], c[0])))


def local_columns(pair, word):
    n = len(pair[0])
    product = identity(n)
    columns = list(identity(n))
    for bit in word:
        matrix = pair[int(bit)]
        columns.append(F.matrix_vector(product, tuple(row[-1] for row in matrix)))
        product = F.matrix_multiply(product, matrix)
    return columns


def choose_pair(feedback, words, triples, budget, record):
    """Select only by local algebra; injected inputs support neutral tests."""
    require(len(feedback) == 3 and feedback[0] != 0, "invertible K3 feedback")
    for value in range(1, 256):
        if value == feedback[0]:
            continue
        budget.check()
        other = (feedback[0] ^ value,) + tuple(feedback[1:])
        pair = (companion(feedback), companion(other))
        item = dict(parameter=value, checked=0, first_failure=None)
        record.append(item)
        for word in words:
            columns = local_columns(pair, word)
            for triple in triples:
                item["checked"] += 1
                if determinant3(tuple(columns[i] for i in triple)) == 0:
                    item["first_failure"] = dict(word=word, columns=list(triple))
                    break
            if item["first_failure"] is not None:
                break
        if item["first_failure"] is None:
            return pair
    return None


def parity(value):
    return bin(value).count("1") & 1


class Mapper:
    """Dimension-parametric version of the K6 packed dyadic lookup."""
    def __init__(self, pair, budget):
        self.n = n = len(pair[0])
        require(len(pair) == 2 and all(len(m) == n and all(len(r) == n for r in m)
                and all(type(v) is int and 0 <= v < 256 for r in m for v in r) for m in pair),
                "lookup pair dimensions/bytes")
        self.blocks = [[pair[0]], [pair[1]]]
        for level in range(31):
            budget.check()
            left, right = self.blocks[0][level], self.blocks[1][level]
            self.blocks[0].append(F.matrix_multiply(left, right))
            self.blocks[1].append(F.matrix_multiply(right, left))
        payloads = []

        def table(bit, width, phase, vectors=False):
            product, data = identity(n), bytearray()
            for value in range(1 << width):
                if value % 32 == 0:
                    budget.check()
                data.extend(bytes(row[0] for row in product) if vectors else
                            bytes(v for row in product for v in row))
                if value + 1 < (1 << width):
                    product = F.matrix_multiply(product, self.blocks[phase ^ parity(value)][bit])
            payloads.append(bytes(data))
            return bytes(data)

        self.low = tuple(table(0, 10, phase, True) for phase in range(2))
        self.mid10 = tuple(table(10, 7, phase) for phase in range(2))
        self.mid17 = tuple(table(17, 7, phase) for phase in range(2))
        self.high = table(24, 8, 0)
        self.payload = b"".join(payloads)
        self.cache = {}

    def apply(self, table, index, vector):
        n = self.n
        offset = index * n * n
        return F.matrix_vector(tuple(tuple(table[offset + r * n + c] for c in range(n))
                                     for r in range(n)), vector)

    def reference_row(self, packet_id):
        # Separately multiply dyadic blocks high-bit first as full matrices.
        product = identity(self.n)
        for bit in range(31, -1, -1):
            if packet_id & (1 << bit):
                product = F.matrix_multiply(product, self.blocks[parity(packet_id >> (bit + 1))][bit])
        return tuple(row[0] for row in product)

    def row(self, packet_id):
        require(type(packet_id) is int and 0 <= packet_id <= MAX_ID, "packet ID")
        if packet_id not in self.cache:
            high, mid17, mid10, low = packet_id >> 24, (packet_id >> 17) & 127, (packet_id >> 10) & 127, packet_id & 1023
            p17 = parity(high)
            p10 = p17 ^ parity(mid17)
            p0 = p10 ^ parity(mid10)
            n = self.n
            vector = self.low[p0][low * n:(low + 1) * n]
            vector = self.apply(self.mid10[p10], mid10, vector)
            vector = self.apply(self.mid17[p17], mid17, vector)
            vector = self.apply(self.high, high, vector)
            require(vector == self.reference_row(packet_id), "lookup/reference disagreement")
            require(any(vector), "zero packet equation")
            self.cache[packet_id] = vector
        return self.cache[packet_id]


def history_inputs(deadline=None):
    """Authenticate and project K3 evidence, never execute historical workers."""
    width, provenance = H.read_bundle(H.WIDTH_LOCAL, "COMPLETE", H.WIDTH_MANIFEST,
                                      H.WIDTH_FILES, False, [0], deadline)
    freeze = H.C.strict_json(width["freeze.json"])
    previous = freeze["history"]
    require(len(previous["origins"]) == 212 and F.digest(canonical(previous["origins"])) ==
            "89e514a5b6a673c0fe5a5188503c0f125717eb35962f48cd7bfa55a3646d8d39", "prior origin pin")
    origins = [row for row in previous["origins"] if row["K"] == 3]
    for name in ("holdout.jsonl", "search.jsonl", "validate.jsonl"):
        for line, data in enumerate(width[name].splitlines(), 1):
            H.C.time_left(deadline)
            row = H.C.strict_json(data)
            if row["type"] == "attempt":
                row = row.get("witness")
                if row is None:
                    continue
            elif row["type"] != "cell":
                continue
            if row["K"] != 3:
                continue
            require(len(row["overheads"]) == len(row["outcomes"]), "history outcomes")
            for overhead, outcome in zip(row["overheads"], row["outcomes"]):
                if outcome == 1:
                    origins.append(dict(issue=58, file=name, line=line, K=3, B=row["B"],
                                        arm=row["arm"], attempt=row["attempt"], root=row["root"],
                                        schedule=row["schedule"], overhead=overhead, ids=row["ids"][:3 + overhead]))
    prefixes = {}
    for row in origins:
        require(row["K"] == 3 and row["B"] in WIDTHS and 0 <= row["overhead"] <= 4 and
                len(H.packet_ids(row["ids"])) == 3 + row["overhead"], "K3 historical origin")
        prefixes.setdefault(tuple(row["ids"]), set()).add(row["B"])
    ledger = [dict(ids=list(ids), original_widths=sorted(widths)) for ids, widths in sorted(prefixes.items())]
    excluded = sorted(set(previous["prior_roots"] + freeze["training_roots"] +
                          freeze["holdout_roots"] + list(H.MAIN_ROOTS) + [
                              "0x" + F.digest(("wirehair.wh2.thue-morse-recovery-r0:fresh/" + str(i)).encode("ascii"))[:16]
                              for i in range(512)]))
    result = dict(provenance=provenance, origins=origins, prefixes=ledger, excluded_roots=excluded)
    # These are authenticated input projections, fixed before candidate scoring.
    for key, count, digest in (
            ("origins", 154, "ff02648eb7cced955312ab56b4e027b13afe471744a234e67d288b5a27057e0d"),
            ("prefixes", 45, "ee13f6d8795ab552705ba026a8fe92dddee65f50ba67ab9085d15b5ce9cce1b0"),
            ("excluded_roots", 809, "c2d46fcdb1c358fa050fea9911e8a1a28a6cb6dc199b5077fbd1100e46467144")):
        require(len(result[key]) == count and F.digest(canonical(result[key])) == digest, "K3 input ledger: " + key)
    return result


def trace(B, root, schedule):
    require(B in WIDTHS and H.root_value(root) and schedule in SCHEDULES, "trace coordinates")
    state = (int(root, 16) ^ 3 * 0x9e3779b97f4a7c15 ^ B * 0xbf58476d1ce4e5b9) & MASK64
    if schedule != "iid":
        state ^= 0x10fade

    def uniform():
        nonlocal state
        state = (state + 0x9e3779b97f4a7c15) & MASK64
        value = ((state ^ (state >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        return ((value ^ (value >> 31)) >> 11) * (2.0 ** -53)

    loss, burst, ids = (0.1 if schedule == "iid" else 0.5), 0, []
    for candidate in range(67328):
        if schedule == "burst" and burst:
            burst -= 1
            continue
        if uniform() < (loss / (8 - 7 * loss) if schedule == "burst" else loss):
            if schedule == "burst":
                burst = 7
            continue
        ids.append(MAX_ID - 2 * candidate if schedule == "adversarial" else
                   3 + candidate if schedule == "repair-only" else candidate)
        if len(ids) == 7:
            return ids
    raise F.ScreenInvalid("trace candidate cap")


def checked_rank(rows):
    rank = F.matrix_rank(rows)
    full = any(determinant3(triple, multiply_polynomial) != 0
               for triple in itertools.combinations(rows, 3))
    if full:
        reference = 3
    else:
        two = any(multiply_polynomial(a[x], b[y]) != multiply_polynomial(a[y], b[x])
                  for a, b in itertools.combinations(rows, 2)
                  for x, y in itertools.combinations(range(3), 2))
        reference = 2 if two else int(any(any(row) for row in rows))
    require(reference == rank, "independent determinant/rank disagreement")
    return rank


def check_window(ids, mapper):
    rows = [mapper.row(packet_id) for packet_id in ids]
    deficient = []
    for selected in TRIPLES:
        triple = tuple(rows[i] for i in selected)
        det = determinant3(triple)
        require(det == determinant3(triple, multiply_polynomial), "determinant oracle disagreement")
        require((det != 0) == (F.matrix_rank(triple) == 3), "rank oracle disagreement")
        if not det:
            deficient.append(list(selected))
    return dict(ids=list(ids), deficient=deficient)


def trace_result(B, root, schedule, mapper):
    ids = trace(B, root, schedule)
    rows = [mapper.row(i) for i in ids]
    ranks = [checked_rank(rows[:length]) for length in range(3, 8)]
    require(all(0 <= b - a <= 1 for a, b in zip(ranks, ranks[1:])), "nested ranks")
    return dict(B=B, root=root, schedule=schedule, ids=ids, ranks=ranks)


def summarize_fresh(rows, per_cell=512):
    require(len(rows) == per_cell * len(WIDTHS) * len(SCHEDULES), "fresh complete denominator")
    cells = []
    for B, schedule in itertools.product(WIDTHS, SCHEDULES):
        selected = [row for row in rows if row["B"] == B and row["schedule"] == schedule]
        require(len(selected) == per_cell and len({row["root"] for row in selected}) == per_cell,
                "fresh cell denominator/duplicate")
        for row in selected:
            require(len(row["ranks"]) == 5 and all(type(rank) is int and 0 <= rank <= 3 for rank in row["ranks"])
                    and all(0 <= b - a <= 1 for a, b in zip(row["ranks"], row["ranks"][1:])), "fresh rank values")
        failures = [sum(row["ranks"][oh] < 3 for row in selected) for oh in range(5)]
        first = [per_cell - failures[0]] + [failures[i - 1] - failures[i] for i in range(1, 5)] + [failures[4]]
        cells.append(dict(B=B, schedule=schedule, traces=per_cell, failures=failures, first_success=first))
    return dict(cells=cells, failures=[sum(cell["failures"][i] for cell in cells) for i in range(5)],
                fresh_pass=all(cell["failures"][0] * 100 <= per_cell for cell in cells))


def run_screen():
    budget = F.Budget()
    result = dict(protocol=PROTOCOL, outcome="INVALID", selection=[], pair=None,
                  local=[], seams=[], development=[], hard=[], history=[], fresh=[],
                  inputs=None, evidence={}, summary={})
    try:
        result["inputs"] = history_inputs(budget.deadline)
        result["evidence"]["field"] = F.init_field()
        pair = choose_pair((8, 14, 7), WORDS, TRIPLES, budget, result["selection"])
        if pair is None:
            result["outcome"] = "EXHAUSTED"
            return result
        result["pair"] = pair
        for word in WORDS:
            columns = local_columns(pair, word)
            for triple in TRIPLES:
                require(checked_rank(tuple(columns[i] for i in triple)) == 3, "selected local certificate")
            result["local"].append(dict(word=word, columns=columns, checked=35))
        mapper = Mapper(pair, budget)
        require(len(mapper.payload) == 13056, "lookup byte count")
        require(tuple(mapper.row(i) for i in range(3)) == identity(3), "systematic rows")
        # Literal sequential products independently validate the dyadic recurrence.
        product = identity(3)
        for packet_id in range(2049):
            require(mapper.row(packet_id) == tuple(row[0] for row in product), "sequential product oracle")
            product = F.matrix_multiply(product, pair[parity(packet_id)])
        result["evidence"].update(lookup_bytes=len(mapper.payload), lookup_sha256=F.digest(mapper.payload),
                                  pair_sha256=F.digest(bytes(v for m in pair for r in m for v in r)))
        starts = [(1 << exponent) - 3 for exponent in range(2, 32)] + [MAX_ID - 6]
        for start in starts:
            budget.check()
            result["seams"].append(check_window(range(start, start + 7), mapper))
        development_hash = hashlib.sha256()
        for trial, root in enumerate(H.HARD_TRAINING_ROOTS):
            for schedule_index, schedule in enumerate(SCHEDULES):
                ids = trace(2, root, schedule)
                ordinal = trial * 120 + schedule_index * 30 + 1
                development_hash.update(struct.pack("<Q7I", ordinal, *ids))
                result["development"].append(check_window(ids, mapper))
        require(development_hash.hexdigest() == DEVELOPMENT_TRACE_SHA, "original K3 trace identity")
        for root, B, schedule in itertools.product(H.HARD_TRAINING_ROOTS + H.HARD_VALIDATION_ROOTS, WIDTHS, SCHEDULES):
            budget.check()
            result["hard"].append(trace_result(B, root, schedule, mapper))
        for prefix in result["inputs"]["prefixes"]:
            budget.check()
            result["history"].append(dict(ids=prefix["ids"], rank=checked_rank([mapper.row(i) for i in prefix["ids"]])))
        structural = (all(not row["deficient"] for row in result["seams"] + result["development"]) and
                      all(row["ranks"][0] == 3 for row in result["hard"]) and
                      all(row["rank"] == 3 for row in result["history"]))
        result["summary"] = dict(structural_pass=structural, fresh_entered=False)
        if structural:
            roots = ["0x" + F.digest((PROTOCOL + ":fresh/" + str(i)).encode("ascii"))[:16] for i in range(512)]
            require(len(set(roots)) == 512 and not set(roots) & set(result["inputs"]["excluded_roots"]), "fresh root collision")
            result["summary"]["fresh_entered"] = True
            for root, B, schedule in itertools.product(roots, WIDTHS, SCHEDULES):
                budget.check()
                result["fresh"].append(trace_result(B, root, schedule, mapper))
            result["summary"].update(summarize_fresh(result["fresh"]))
        result["evidence"]["unique_rows"] = [dict(id=i, row=mapper.cache[i]) for i in sorted(mapper.cache)]
        result["evidence"]["inputs_sha256"] = F.digest(canonical(result["inputs"]))
        result["counts"] = dict(local_triples=350, seam_triples=31 * 35, development_triples=420,
                                hard_traces=len(result["hard"]), history_prefixes=len(result["history"]),
                                fresh_traces=len(result["fresh"]), unique_rows=len(mapper.cache))
        require(len(result["local"]) == 10 and len(result["seams"]) == 31 and len(result["development"]) == 12
                and len(result["hard"]) == 72 and len(result["history"]) == len(result["inputs"]["prefixes"])
                and len(result["fresh"]) == (6144 if structural else 0), "complete screen accounting")
        budget.check()
        result["outcome"] = "PASS" if structural and result["summary"]["fresh_pass"] else "FAIL"
    except Exception as error:
        result["outcome"] = "INVALID"
        result["error"] = (type(error).__name__ + ": " + str(error))[:1024]
    return result


def main(argv):
    if argv != ["--worker"]:
        sys.stderr.write("usage: Wh2K3ThueMorseR0.py --worker\n")
        return 2
    resource.setrlimit(resource.RLIMIT_AS, (512 * 1024 * 1024, 512 * 1024 * 1024))
    result = run_screen()
    raw = canonical(result) + b"\n"
    require(len(raw) <= F.STDOUT_LIMIT, "worker output cap")
    sys.stdout.buffer.write(raw)
    sys.stdout.buffer.flush()
    return int(result["outcome"] == "INVALID")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
