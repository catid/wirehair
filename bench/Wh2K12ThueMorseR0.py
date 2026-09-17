#!/usr/bin/env python3
"""Finite K12 GF(256) structural screen; never a codec or timing claim.

The worker is launched only by Wh2K12ThueMorseRunR0.py after an exact claim has
been sealed.  Importing this module is inert and never selects a candidate.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import struct
import sys
import time

PROTOCOL = 'wirehair.wh2.k12-thue-morse-r1'
POLYNOMIAL = 0x14d
K = 12
MAX_ID = (1 << 32) - 1
MASK64 = (1 << 64) - 1
WIDTHS = (2, 64, 1280)
SCHEDULES = ('iid', 'burst', 'adversarial', 'repair-only')
WORDS = ('0010', '0011', '0100', '0101', '0110',
         '1001', '1010', '1011', '1100', '1101')
MINORS = tuple(itertools.combinations(range(K + 4), K))
ROOTS = 512
CANDIDATES = tuple(range(1, 256))
OUTPUT_LIMIT = 4 * 1024 * 1024
DEADLINE_SECONDS = 60
INVENTORY = Path('/var/tmp/wh2-uncovered-band-inventory-r0')
OUTPUT = Path('/var/tmp/wh2-k12-thue-morse-r1')
INVENTORY_COMPLETE_SHA = '16fcf13214cd25362fe66ee35f62c1c9616ed082e343fa5e2009d81d61359ea0'


class ScreenInvalid(ValueError):
    pass


def require(ok, reason):
    if not ok:
        raise ScreenInvalid(reason)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True, allow_nan=False) + '\n').encode('ascii')


def digest(value):
    return hashlib.sha256(value).hexdigest()


class Budget:
    def __init__(self, seconds=DEADLINE_SECONDS):
        self.deadline = time.monotonic() + seconds

    def check(self):
        require(time.monotonic() < self.deadline, 'worker deadline')


def gf_mul(x, y):
    require(type(x) is int and type(y) is int and 0 <= x < 256 and 0 <= y < 256,
            'GF256 byte range')
    result = 0
    while y:
        if y & 1:
            result ^= x
        x <<= 1
        if x & 256:
            x ^= POLYNOMIAL
        y >>= 1
    return result


MUL = tuple(bytes(gf_mul(x, y) for y in range(256)) for x in range(256))
INV = (0,) + tuple(next(b for b in range(1, 256) if MUL[a][b] == 1)
                  for a in range(1, 256))


def identity(n=K):
    return tuple(tuple(int(r == c) for c in range(n)) for r in range(n))


def matrix_bytes(matrix):
    require(len(matrix) == K and all(len(row) == K for row in matrix), 'matrix shape')
    require(all(type(x) is int and 0 <= x < 256 for row in matrix for x in row),
            'matrix element range')
    return bytes(x for row in matrix for x in row)


def matrix_multiply(left, right):
    require(len(left) == K and len(right) == K, 'matrix dimensions')
    return tuple(tuple(_dot(row, column) for column in zip(*right)) for row in left)


def _dot(left, right):
    value = 0
    for x, y in zip(left, right):
        value ^= MUL[x][y]
    return value


def matrix_vector(matrix, vector):
    return tuple(_dot(row, vector) for row in matrix)


def matrix_rank(rows):
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    width = len(rows[0])
    require(all(len(row) == width and
                all(type(x) is int and 0 <= x < 256 for x in row) for row in rows),
            'rank rectangle')
    pivot = 0
    for column in range(width):
        chosen = next((r for r in range(pivot, len(rows)) if rows[r][column]), None)
        if chosen is None:
            continue
        rows[pivot], rows[chosen] = rows[chosen], rows[pivot]
        scale = MUL[INV[rows[pivot][column]]]
        rows[pivot] = [scale[x] for x in rows[pivot]]
        for r in range(pivot + 1, len(rows)):
            if rows[r][column]:
                table = MUL[rows[r][column]]
                rows[r] = [x ^ table[y] for x, y in zip(rows[r], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def matrix_inverse(matrix):
    matrix_bytes(matrix)
    rows = [list(row) + list(unit) for row, unit in zip(matrix, identity())]
    for column in range(K):
        chosen = next((r for r in range(column, K) if rows[r][column]), None)
        require(chosen is not None, 'singular inverse')
        rows[column], rows[chosen] = rows[chosen], rows[column]
        scale = MUL[INV[rows[column][column]]]
        rows[column] = [scale[x] for x in rows[column]]
        for r in range(K):
            if r != column and rows[r][column]:
                table = MUL[rows[r][column]]
                rows[r] = [x ^ table[y] for x, y in zip(rows[r], rows[column])]
    return tuple(tuple(row[K:]) for row in rows)


def companion(feedback):
    require(len(feedback) == K and feedback[0] != 0 and
            all(type(x) is int and 0 <= x < 256 for x in feedback), 'feedback')
    return tuple(tuple(feedback[r] if c == K - 1 else int(r == c + 1)
                       for c in range(K)) for r in range(K))


def fixed_feedback():
    polynomial = [1]
    value = 1
    roots = []
    for _ in range(K):
        roots.append(value)
        value = gf_mul(value, 2)
    require(len(set(roots)) == K, 'distinct feedback roots')
    for root in roots:
        output = [0] * (len(polynomial) + 1)
        for i, coefficient in enumerate(polynomial):
            output[i] ^= MUL[root][coefficient]
            output[i + 1] ^= coefficient
        polynomial = output
    require(len(polynomial) == K + 1 and polynomial[-1] == 1 and polynomial[0],
            'monic invertible feedback')
    return tuple(polynomial[:-1])


def parity(value):
    return bin(value).count('1') & 1


def local_columns(pair, word):
    product = identity()
    columns = list(identity())
    for bit in word:
        matrix = pair[int(bit)]
        columns.append(matrix_vector(product, tuple(row[-1] for row in matrix)))
        product = matrix_multiply(product, matrix)
    require(len(columns) == K + 4, 'local column count')
    return columns


def choose_pair(feedback, budget):
    records = []
    for value in CANDIDATES:
        budget.check()
        if value == feedback[0]:
            continue
        pair = (companion(feedback), companion((feedback[0] ^ value,) + feedback[1:]))
        entry = dict(parameter=value, checked=0, first_failure=None)
        records.append(entry)
        for word in WORDS:
            columns = local_columns(pair, word)
            for selected in MINORS:
                entry['checked'] += 1
                if matrix_rank([columns[i] for i in selected]) != K:
                    entry['first_failure'] = dict(word=word, columns=list(selected))
                    break
            if entry['first_failure'] is not None:
                break
        if entry['first_failure'] is None:
            return pair, records
    return None, records


class Mapper:
    """Packed dyadic K12 mapping, independently checked against full products."""
    def __init__(self, pair, budget):
        self.pair = pair
        self.blocks = [[pair[0]], [pair[1]]]
        for level in range(31):
            budget.check()
            self.blocks[0].append(matrix_multiply(self.blocks[0][level], self.blocks[1][level]))
            self.blocks[1].append(matrix_multiply(self.blocks[1][level], self.blocks[0][level]))
        payloads = []

        def table(bit, width, phase, vectors=False):
            product = identity()
            data = bytearray()
            for value in range(1 << width):
                if value % 32 == 0:
                    budget.check()
                data.extend(bytes(row[0] for row in product) if vectors else matrix_bytes(product))
                if value + 1 < (1 << width):
                    product = matrix_multiply(product, self.blocks[phase ^ parity(value)][bit])
            payloads.append(bytes(data))
            return bytes(data)

        self.low = tuple(table(0, 10, phase, True) for phase in range(2))
        self.mid10 = tuple(table(10, 7, phase) for phase in range(2))
        self.mid17 = tuple(table(17, 7, phase) for phase in range(2))
        self.high = table(24, 8, 0)
        self.payload = b''.join(payloads)
        require(len(self.payload) == 135168, 'exact K12 lookup geometry')
        self.cache = {}

    def apply(self, table, index, vector):
        offset = index * K * K
        matrix = tuple(tuple(table[offset + r * K + c] for c in range(K)) for r in range(K))
        return matrix_vector(matrix, vector)

    def reference_row(self, packet_id):
        product = identity()
        for bit in range(31, -1, -1):
            if packet_id & (1 << bit):
                product = matrix_multiply(product,
                    self.blocks[parity(packet_id >> (bit + 1))][bit])
        return tuple(row[0] for row in product)

    def row(self, packet_id):
        require(type(packet_id) is int and 0 <= packet_id <= MAX_ID, 'packet ID')
        if packet_id not in self.cache:
            high, mid17, mid10 = packet_id >> 24, (packet_id >> 17) & 127, (packet_id >> 10) & 127
            low = packet_id & 1023
            phase17 = parity(high)
            phase10 = phase17 ^ parity(mid17)
            phase0 = phase10 ^ parity(mid10)
            vector = self.low[phase0][low * K:(low + 1) * K]
            vector = self.apply(self.mid10[phase10], mid10, vector)
            vector = self.apply(self.mid17[phase17], mid17, vector)
            vector = self.apply(self.high, high, vector)
            require(vector == self.reference_row(packet_id), 'lookup/reference disagreement')
            require(any(vector), 'zero packet equation')
            self.cache[packet_id] = vector
        return self.cache[packet_id]


def trace(width, root, schedule):
    state = (int(root, 16) ^ K * 0x9e3779b97f4a7c15 ^ width * 0xbf58476d1ce4e5b9) & MASK64
    if schedule != 'iid':
        state ^= 0x10fade
    loss = .1 if schedule == 'iid' else .5
    threshold = (loss / (8 - 7 * loss) if schedule == 'burst' else loss) * 2.0**53
    result, burst = [], 0
    for candidate in range(65536):
        if schedule == 'burst' and burst:
            burst -= 1
            continue
        state = (state + 0x9e3779b97f4a7c15) & MASK64
        value = ((state ^ (state >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        accepted = ((value ^ (value >> 31)) >> 11) >= threshold
        if not accepted:
            if schedule == 'burst':
                burst = 7
            continue
        result.append(MAX_ID - 2 * candidate if schedule == 'adversarial' else
                      K + candidate if schedule == 'repair-only' else candidate)
        if len(result) == K + 4:
            return result
    raise ScreenInvalid('trace candidate cap')


def history_inputs(deadline=None):
    if deadline is not None:
        require(time.monotonic() < deadline, 'history deadline')
    complete = INVENTORY / 'complete.json'
    complete_raw = complete.read_bytes()
    require(complete.is_file() and digest(complete_raw) == INVENTORY_COMPLETE_SHA,
            'sealed baseline inventory')
    manifest = json.loads(complete_raw, object_pairs_hook=dict,
                          parse_constant=lambda value: (_ for _ in ()).throw(
                              ValueError('nonfinite baseline manifest')))
    require(canonical(manifest) == complete_raw, 'canonical baseline manifest')
    require(set(manifest) == {'analysis.json', 'claim.json', 'process.json', 'raw.jsonl', 'stderr.txt'},
            'baseline complete members')
    require({path.name for path in INVENTORY.iterdir()} ==
            set(manifest) | {'complete.json'}, 'exact baseline directory members')
    raw_path = INVENTORY / 'raw.jsonl'
    raw = raw_path.read_bytes()
    require(len(raw) == 86272424, 'baseline raw byte count')
    require(digest(raw) == '3d26cab298d0bebb7f293893aaf69e3ea9fb19b72dd347333b022cdc923c97c8',
            'baseline raw content')
    for name, expected in manifest.items():
        require(name in {'analysis.json', 'claim.json', 'process.json', 'raw.jsonl', 'stderr.txt'},
                'baseline manifest name')
        path = INVENTORY / name
        require(path.is_file() and path.stat().st_mode & 0o777 == 0o400,
                'sealed baseline member')
        data = raw if name == 'raw.jsonl' else path.read_bytes()
        require(digest(data) == expected, 'baseline member hash')
        if deadline is not None:
            require(time.monotonic() < deadline, 'history deadline')
    origins, roots = [], set()
    lines = raw.splitlines()
    case_ordinal = 0
    for line in lines:
        if deadline is not None:
            require(time.monotonic() < deadline, 'history deadline')
        row = json.loads(line, object_pairs_hook=dict,
                         parse_constant=lambda value: (_ for _ in ()).throw(
                             ValueError('nonfinite baseline record')))
        require(canonical(row) == line + b'\n', 'canonical baseline record')
        if row.get('type') != 'case':
            continue
        require(row['ordinal'] == case_ordinal, 'baseline chronology')
        case_ordinal += 1
        case = row['case']
        if case['root']:
            roots.add(case['root'])
        if case['k'] != K:
            continue
        for arm, result in enumerate(row['arms']):
            for overhead in range(5):
                if not result['first'] or result['first'] > K + overhead:
                    origins.append(dict(ordinal=row['ordinal'], arm=arm, b=case['b'], tail=case['tail'],
                                        group=case['group'], root=case['root'], schedule=case['schedule'],
                                        overhead=overhead, ids=case['ids'][:K + overhead]))
    prefixes = {}
    for origin in origins:
        prefixes.setdefault(tuple(origin['ids']), set()).add(origin['b'])
    result = dict(origins=origins,
                  prefixes=[dict(ids=list(ids), original_widths=sorted(widths))
                            for ids, widths in sorted(prefixes.items())],
                  roots=sorted(roots))
    require(len(origins) == 63 and len(result['prefixes']) == 60 and len(result['roots']) == 64,
            'complete K12 history projection')
    require(case_ordinal == 3096, 'complete baseline case chronology')
    require(digest(canonical(origins)) == '21c1ffe5cf89d8bc99846b51ff50d7e9d7b9a7bd7c723a2ab73cfc4a7f38b2ec',
            'K12 history origins pin')
    require(digest(canonical(result['prefixes'])) == '9ba913961b79f21992b01c7f9750d1604307dca655258a989c4342601c3a2026',
            'K12 history prefix pin')
    require(digest(canonical(result['roots'])) == 'a3bc3986862980d0698a9d5b2850879d113ea2a2c94bc237392926b607707e03',
            'K12 history roots pin')
    return result


def fresh_roots(history):
    excluded = set(history['roots'])
    for prior in ('wirehair.wh2.k2-thue-morse-r0', 'wirehair.wh2.k3-thue-morse-r0',
                  'wirehair.wh2.k5-thue-morse-r0', 'wirehair.wh2.k8-thue-morse-r0',
                  'wirehair.wh2.thue-morse-recovery-r0'):
        excluded.update('0x' + digest((prior + ':fresh/' + str(i)).encode())[:16]
                        for i in range(ROOTS))
    roots = ['0x' + digest((PROTOCOL + ':fresh/' + str(i)).encode())[:16]
             for i in range(ROOTS)]
    require(len(set(roots)) == ROOTS and not set(roots) & excluded, 'fresh root collision')
    return roots


def claimed_inputs():
    """Authenticate the controller's immutable claim before candidate work."""
    path = OUTPUT / 'CLAIM.json'
    raw = path.read_bytes()
    require(len(raw) <= 1024 * 1024, 'claim size')
    claim = json.loads(raw, object_pairs_hook=dict,
                       parse_constant=lambda value: (_ for _ in ()).throw(
                           ValueError('nonfinite claim')))
    compact_json = lambda value: json.dumps(value, sort_keys=True, separators=(',', ':'),
                                             ensure_ascii=True, allow_nan=False).encode('ascii')
    compact = compact_json(claim)
    require(compact == raw and
            set(claim) == {'protocol', 'receipt_sha256', 'receipt'} and
            claim['protocol'] == PROTOCOL and
            claim['receipt_sha256'] == digest(compact_json(claim['receipt'])),
            'claim identity')
    return digest(raw)


def check_mapper(mapper, budget):
    require(matrix_rank(mapper.pair[0]) == K and matrix_rank(mapper.pair[1]) == K,
            'invertible candidate pair')
    require(matrix_multiply(mapper.pair[0], mapper.pair[1]) !=
            matrix_multiply(mapper.pair[1], mapper.pair[0]),
            'noncommuting candidate pair')
    require(tuple(mapper.row(i) for i in range(K)) == identity(), 'systematic K12 rows')
    product = identity()
    for i in range(2049):
        budget.check()
        require(mapper.row(i) == tuple(row[0] for row in product), 'sequential row oracle')
        product = matrix_multiply(product, mapper.pair[parity(i)])
    seams = []
    for exponent in range(3, 32):
        for offset in (-1, 0, 1):
            packet_id = (1 << exponent) + offset
            seams.append(dict(id=packet_id, row=list(mapper.row(packet_id))))
    seams.append(dict(id=MAX_ID, row=list(mapper.row(MAX_ID))))
    return seams


def check_history(mapper, history, budget):
    checked = []
    for origin in history['origins']:
        budget.check()
        rows = [mapper.row(packet_id) for packet_id in origin['ids']]
        checked.append(dict(origin=origin, rank=matrix_rank(rows)))
        require(checked[-1]['rank'] == K, 'retained history full rank')
    return checked


def trace_result(mapper, width, root, schedule, budget):
    ids = trace(width, root, schedule)
    rows = [mapper.row(packet_id) for packet_id in ids]
    ranks = [matrix_rank(rows[:count]) for count in range(K, K + 5)]
    require(all(0 <= later - earlier <= 1 for earlier, later in zip(ranks, ranks[1:])),
            'nested fresh ranks')
    budget.check()
    return dict(b=width, root=root, schedule=schedule, ids=ids, ranks=ranks)


def summarize_fresh(rows, per_cell=ROOTS):
    require(len(rows) == per_cell * len(WIDTHS) * len(SCHEDULES), 'fresh denominator')
    cells = []
    for width, schedule in itertools.product(WIDTHS, SCHEDULES):
        selected = [r for r in rows if r['b'] == width and r['schedule'] == schedule]
        require(len(selected) == per_cell and len({r['root'] for r in selected}) == per_cell,
                'fresh cell denominator')
        failures = [sum(r['ranks'][overhead] < K for r in selected) for overhead in range(5)]
        cells.append(dict(b=width, schedule=schedule, traces=per_cell, failures=failures,
                          first_success=[per_cell - failures[0]] +
                          [failures[i - 1] - failures[i] for i in range(1, 5)] + [failures[4]]))
    failures = [sum(cell['failures'][i] for cell in cells) for i in range(5)]
    return dict(cells=cells, failures=failures,
                fresh_pass=failures[0] * 100 <= len(rows) and
                all(cell['failures'][0] * 100 <= per_cell for cell in cells))


def run_screen(claim):
    budget = Budget()
    history = history_inputs()
    result = dict(protocol=PROTOCOL, claim_sha256=claim, outcome='INVALID', selection=[],
                  pair=None, history_projection=dict(origins=len(history['origins']),
                  prefixes=len(history['prefixes']), roots=len(history['roots'])),
                  seams=[], history=[], fresh=[], summary={})
    feedback = fixed_feedback()
    pair, selection = choose_pair(feedback, budget)
    result['selection'] = selection
    if pair is None:
        result['outcome'] = 'EXHAUSTED'
        return result
    result['pair'] = pair
    result['feedback'] = feedback
    mapper = Mapper(pair, budget)
    result['lookup'] = dict(bytes=len(mapper.payload), sha256=digest(mapper.payload))
    result['seams'] = check_mapper(mapper, budget)
    result['history'] = check_history(mapper, history, budget)
    if any(row['rank'] != K for row in result['history']):
        result['outcome'] = 'FAIL'
        return result
    fresh = []
    for root, width, schedule in itertools.product(fresh_roots(history), WIDTHS, SCHEDULES):
        fresh.append(trace_result(mapper, width, root, schedule, budget))
    result['fresh'] = fresh
    result['summary'] = summarize_fresh(fresh)
    result['outcome'] = 'PASS' if result['summary']['fresh_pass'] else 'FAIL'
    result['counts'] = dict(local_words=len(WORDS), local_minors=len(WORDS) * len(MINORS),
                            seams=len(result['seams']), history=len(result['history']),
                            fresh=len(fresh), unique_rows=len(mapper.cache))
    return result


def set_limits():
    resource.setrlimit(resource.RLIMIT_CPU, (DEADLINE_SECONDS, DEADLINE_SECONDS))
    resource.setrlimit(resource.RLIMIT_AS, (512 * 1024**2, 512 * 1024**2))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (OUTPUT_LIMIT, OUTPUT_LIMIT))


def main(argv):
    require(argv == ['--worker'], 'usage')
    set_limits()
    claim = claimed_inputs()
    result = run_screen(claim)
    raw = canonical(result)
    require(len(raw) <= OUTPUT_LIMIT, 'worker output cap')
    sys.stdout.buffer.write(raw)
    sys.stdout.buffer.flush()
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except BaseException as error:
        print(type(error).__name__ + ': ' + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
