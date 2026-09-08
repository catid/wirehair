#!/usr/bin/env python3
"""Sealed .82 evidence to native build data; never select or score candidates."""
import argparse
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('_k3_native_reader', HERE / 'Wh2NoncommutingRadixRunR0.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
ROOT = Path('/var/tmp/wh2-k3-thue-morse-r0')
MANIFEST_SHA = '96b50c60249b0dd3a010d6550bfa0a22c7ae2d3d86746aa153177a6a0a19ffc5'
RAW_SHA = '28ea54fadc3474e411066d76b8f1abe040b586064a5e9195899ee57aac3ef3df'
LOOKUP_SHA = 'c78d6f350767bc5336f36eae30424914347f42b4314465b590fee6c1612e9d15'
MEMBERS = ('CLAIM.json', 'raw.json', 'stderr.txt', 'summary.json')
WIDTHS = (2, 64, 1280)
CAP = 4 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_report():
    raw = C.read_regular(ROOT / 'COMPLETE.json', CAP)
    require(C.sha(raw) == MANIFEST_SHA, 'manifest identity')
    manifest = C.strict_json(raw)
    require(raw == C.canonical(manifest) and manifest['protocol'] == 'wirehair.wh2.k3-thue-morse-r0'
            and manifest['outcome'] == 'PASS' and set(manifest['files']) == set(MEMBERS), 'manifest schema')
    require({p.name for p in ROOT.iterdir()} == set(MEMBERS) | {'COMPLETE.json'}, 'bundle roster')
    content = {}
    for name in MEMBERS:
        content[name] = C.read_regular(ROOT / name, CAP)
        require(manifest['files'][name] == dict(bytes=len(content[name]), sha256=C.sha(content[name])), 'member identity')
    require(not content['stderr.txt'] and C.sha(content['raw.json']) == RAW_SHA, 'raw/stderr identity')
    report = C.strict_json(content['raw.json'])
    require(content['raw.json'] == C.canonical(report) + b'\n' and report['outcome'] == 'PASS', 'report encoding/outcome')
    require(report['evidence']['lookup_sha256'] == LOOKUP_SHA, 'recorded lookup identity')
    return report


def integer(value, low, high):
    require(type(value) is int and low <= value <= high, 'integer range/type')
    return value


def ids(values, count):
    require(type(values) is list and len(values) == count and len(set(values)) == count, 'ID count/uniqueness')
    return [integer(v, 0, (1 << 32) - 1) for v in values]


def extract(report):
    """Copy recorded IDs/ranks only; no trace generation or matrix scoring."""
    require(len(report['fresh']) == 6144 and len(report['hard']) == 72 and len(report['history']) == 45,
            'trace/history roster')
    traces = []
    for row in report['fresh'] + report['hard']:
        require(type(row['B']) is int and row['B'] in WIDTHS and row['ranks'] == [3] * 5
                and all(type(v) is int for v in row['ranks']), 'retained full-rank trace')
        traces.append(dict(B=row['B'], ids=ids(row['ids'], 7), ranks=list(row['ranks'])))
    require(len(report['inputs']['prefixes']) == 45, 'history input roster')
    history = []
    for recorded, prefix in zip(report['history'], report['inputs']['prefixes']):
        require(recorded['ids'] == prefix['ids'] and recorded['rank'] == 3, 'history concordance')
        count = integer(len(recorded['ids']), 3, 7)
        widths = prefix['original_widths']
        require(widths and widths == sorted(set(widths)) and all(type(B) is int and B in WIDTHS for B in widths),
                'historical width roster')
        history.append(dict(ids=ids(prefix['ids'], count), widths=widths))
    require(len(report['development']) == 12 and len(report['seams']) == 31, 'window roster')
    windows = []
    for row in report['development'] + report['seams']:
        require(row['deficient'] == [], 'window certificate')
        windows.append(ids(row['ids'], 7))
    rows = report['evidence']['unique_rows']
    require(len(rows) == 2226 and [row['id'] for row in rows] == sorted({row['id'] for row in rows}), 'row roster')
    for row in rows:
        integer(row['id'], 0, (1 << 32) - 1)
        require(len(row['row']) == 3, 'row shape')
        for value in row['row']: integer(value, 0, 255)
    pair = report['pair']
    require(pair == [[[0,0,8],[1,0,14],[0,1,7]], [[0,0,9],[1,0,14],[0,1,7]]], 'sealed companion pair')
    return dict(pair=pair, traces=traces, history=history, windows=windows, rows=rows)


def multiply(a, b):
    product = 0
    for bit in range(8):
        if b & (1 << bit): product ^= a << bit
    for bit in range(14, 7, -1):
        if product & (1 << bit): product ^= 0x14d << (bit - 8)
    return product


def build_lookup(pair):
    """Reconstruct fixed byte tables, not rows, ranks or recovery scores."""
    require(len(pair) == 2 and all(len(m) == 3 and all(len(r) == 3 and
            all(type(v) is int and 0 <= v < 256 for v in r) for r in m) for m in pair), 'lookup pair shape')
    table = [[multiply(a, b) for b in range(256)] for a in range(256)]

    def product(a, b):
        return [[table[a[r][0]][b[0][c]] ^ table[a[r][1]][b[1][c]] ^ table[a[r][2]][b[2][c]]
                 for c in range(3)] for r in range(3)]

    blocks = [pair]
    for level in range(31):
        a, b = blocks[level]
        blocks.append([product(a, b), product(b, a)])
    output = bytearray()
    for start, width, phase, vectors in ((0,10,0,True), (0,10,1,True), (10,7,0,False),
            (10,7,1,False), (17,7,0,False), (17,7,1,False), (24,8,0,False)):
        p = [[1,0,0],[0,1,0],[0,0,1]]
        for index in range(1 << width):
            output.extend([r[0] for r in p] if vectors else [v for r in p for v in r])
            p = product(p, blocks[start][phase ^ (bin(index).count('1') & 1)])
    return bytes(output)


def render(data, lookup):
    require(len(lookup) == 13056 and C.sha(lookup) == LOOKUP_SHA, 'lookup bytes/hash')
    lines = ['// Generated only from authenticated sealed K3 evidence.', '#include <cstdint>',
             'namespace wh2_k3_data {', 'struct Trace { unsigned B; std::uint32_t ids[7]; unsigned ranks[5]; };',
             'struct Prefix { unsigned count, widths; std::uint32_t ids[7]; };',
             'struct Row { std::uint32_t id; std::uint8_t values[3]; };',
             'static const char kRawSha[] = "' + RAW_SHA + '";',
             'static const char kLookupSha[] = "' + LOOKUP_SHA + '";',
             'alignas(64) static const std::uint8_t kLookup[13056] = {']
    for offset in range(0, len(lookup), 24):
        lines.append(','.join(str(v) for v in lookup[offset:offset + 24]) + ',')
    lines.append('};\nstatic const Trace kTraces[] = {')
    for row in data['traces']:
        lines.append('{%d,{%s},{%s}},' % (row['B'], ','.join(str(v) + 'u' for v in row['ids']),
                                        ','.join(str(v) for v in row['ranks'])))
    lines.append('};\nstatic const Prefix kHistory[] = {')
    for row in data['history']:
        mask = sum(1 << WIDTHS.index(B) for B in row['widths'])
        lines.append('{%d,%d,{%s}},' % (len(row['ids']), mask, ','.join(str(v) + 'u' for v in row['ids'])))
    lines.append('};\nstatic const std::uint32_t kWindows[][7] = {')
    lines.extend('{' + ','.join(str(v) + 'u' for v in row) + '},' for row in data['windows'])
    lines.append('};\nstatic const Row kRows[] = {')
    lines.extend('{%du,{%s}},' % (r['id'], ','.join(str(v) for v in r['row'])) for r in data['rows'])
    lines.append('};\n} // namespace wh2_k3_data\n')
    raw = '\n'.join(lines).encode('ascii')
    require(len(raw) <= CAP, 'generated header cap')
    return raw


def write_header(path, raw):
    # A build may request unchanged output again; never overwrite different bytes.
    try:
        C.write_new(path, raw)
    except FileExistsError:
        require(C.read_regular(path, CAP) == raw, 'existing generated header differs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--header', required=True, type=Path)
    args = parser.parse_args()
    report = load_report()
    data = extract(report)
    raw = render(data, build_lookup(data['pair']))
    write_header(args.header, raw)
    print('K3_DATA bytes=%d sha256=%s history_width_cases=%d' %
          (len(raw), C.sha(raw), sum(len(r['widths']) for r in data['history'])))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(type(error).__name__ + ': ' + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
