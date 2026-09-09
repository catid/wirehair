#!/usr/bin/env python3
"""Reconstruct native K5 fixtures from sealed evidence; never select or score."""
import argparse
import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('_k5_native_reader', HERE / 'Wh2K3NativeDataR0.py')
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)
C, require, integer, ids = R.C, R.require, R.integer, R.ids
ROOT = Path('/var/tmp/wh2-k5-thue-morse-r0')
PROTOCOL = 'wirehair.wh2.k5-thue-morse-r0'
MANIFEST_SHA = '04047ef22709785fd94be51372359d1735e09a10cb6f081763d909760baed2d2'
RAW_SHA = 'acddcd2d6c8dbfa7284aa900980b2845a03130b5f2f721515ac8d45bb580af76'
LOOKUP_SHA = '4ac8059aba3b5797c8789c4258a1bda52e5fdb005592d601591705940cfe76c9'
WIDTHS, CAP = R.WIDTHS, R.CAP
FEEDBACK = (121, 110, 207, 198, 31)
PAIR = [[[int(r == c + 1) if c < 4 else FEEDBACK[r] ^ (phase if r == 0 else 0)
          for c in range(5)] for r in range(5)] for phase in range(2)]


def load_report():
    raw = C.read_regular(ROOT / 'COMPLETE.json', CAP)
    require(C.sha(raw) == MANIFEST_SHA, 'manifest identity')
    manifest = C.strict_json(raw)
    require(raw == C.canonical(manifest) and manifest['protocol'] == PROTOCOL
            and manifest['outcome'] == 'PASS' and set(manifest['files']) == set(R.MEMBERS), 'manifest schema')
    require({p.name for p in ROOT.iterdir()} == set(R.MEMBERS) | {'COMPLETE.json'}, 'bundle roster')
    content = {}
    for name in R.MEMBERS:
        content[name] = C.read_regular(ROOT / name, CAP)
        require(manifest['files'][name] == dict(bytes=len(content[name]), sha256=C.sha(content[name])), 'member identity')
    require(not content['stderr.txt'] and C.sha(content['raw.json']) == RAW_SHA, 'raw/stderr identity')
    report = C.strict_json(content['raw.json'])
    require(content['raw.json'] == C.canonical(report) + b'\n' and report['outcome'] == 'PASS'
            and report['protocol'] == PROTOCOL, 'report encoding/outcome')
    require(report['evidence']['lookup_sha256'] == LOOKUP_SHA and
            report['evidence']['lookup_bytes'] == 29440, 'recorded lookup identity')
    return report


def extract(report):
    """Copy exact recorded IDs/ranks and original widths, without matrix scoring."""
    require(report['pair'] == PAIR, 'sealed companion pair')
    require(len(report['fresh']) == 6144 and len(report['hard']) == 72, 'trace roster')
    traces = []
    for row in report['fresh'] + report['hard']:
        require(type(row['B']) is int and row['B'] in WIDTHS and len(row['ranks']) == 5, 'trace shape')
        ranks = [integer(v, 0, 5) for v in row['ranks']]
        require(ranks == sorted(ranks) and ranks[1:] == [5] * 4, 'retained prefix ranks')
        traces.append(dict(B=row['B'], ids=ids(row['ids'], 9), ranks=ranks))
    require(sum(t['ranks'][0] != 5 for t in traces[:6144]) == 11 and
            all(t['ranks'][0] == 5 for t in traces[6144:]), 'fresh/hard rank counts')
    prefixes = report['inputs']['prefixes']
    require(len(prefixes) == 54 and len(report['history']) == 54 and
            len(report['inputs']['origins']) == 79, 'history roster')
    widths = {}
    for origin in report['inputs']['origins']:
        p = tuple(ids(origin['ids'], integer(len(origin['ids']), 5, 9)))
        require(type(origin['b']) is int and origin['b'] in WIDTHS, 'original width')
        widths.setdefault(p, set()).add(origin['b'])
    history = []
    for recorded, prefix in zip(report['history'], prefixes):
        p = ids(prefix, integer(len(prefix), 5, 9))
        require(recorded['ids'] == p and type(recorded['rank']) is int and recorded['rank'] == 5
                and tuple(p) in widths, 'history concordance')
        history.append(dict(ids=p, widths=sorted(widths[tuple(p)])))
    require(len({tuple(r['ids']) for r in history}) == 54 and len(widths) == 54 and
            sum(len(r['widths']) for r in history) == 57 and
            sum(len(r['ids']) * len(r['widths']) for r in history) == 290, 'history width coverage')
    require(len(report['seams']) == 30, 'window roster')
    windows = []
    for window in report['seams']:
        require(window['deficient'] == [], 'window certificate')
        windows.append(ids(window['ids'], 9))
    rows = report['evidence']['unique_rows']
    require(len(rows) == 2270 and [r['id'] for r in rows] == sorted({r['id'] for r in rows}), 'row roster')
    for row in rows:
        integer(row['id'], 0, (1 << 32) - 1)
        require(len(row['row']) == 5, 'row shape')
        for value in row['row']: integer(value, 0, 255)
    return dict(traces=traces, history=history, windows=windows, rows=rows)


def build_lookup():
    """Fixed pair only. Independent carryless polynomial multiplication."""
    table = [[R.multiply(a, b) for b in range(256)] for a in range(256)]

    def product(a, b):
        result = [[0] * 5 for _ in range(5)]
        for r in range(5):
            for c in range(5):
                for k in range(5): result[r][c] ^= table[a[r][k]][b[k][c]]
        return result

    blocks = [PAIR]
    for level in range(31):
        a, b = blocks[level]
        blocks.append([product(a, b), product(b, a)])
    output = bytearray()
    for start, width, phase, vectors in ((0,10,0,True), (0,10,1,True), (10,7,0,False),
            (10,7,1,False), (17,7,0,False), (17,7,1,False), (24,8,0,False)):
        p = [[int(r == c) for c in range(5)] for r in range(5)]
        for index in range(1 << width):
            output.extend([r[0] for r in p] if vectors else [v for r in p for v in r])
            p = product(p, blocks[start][phase ^ (bin(index).count('1') & 1)])
    require(len(output) == 29440 and C.sha(output) == LOOKUP_SHA, 'lookup bytes/hash')
    return bytes(output)


def render(data, lookup):
    require(len(lookup) == 29440 and C.sha(lookup) == LOOKUP_SHA, 'lookup bytes/hash')
    lines = ['// Generated only from authenticated sealed K5 evidence.', '#include <cstdint>',
             'namespace wh2_k5_data {', 'struct Trace { unsigned B; std::uint32_t ids[9]; unsigned ranks[5]; };',
             'struct Prefix { unsigned count, widths; std::uint32_t ids[9]; };',
             'struct Row { std::uint32_t id; std::uint8_t values[5]; };',
             'static const char kRawSha[] = "' + RAW_SHA + '";',
             'static const char kLookupSha[] = "' + LOOKUP_SHA + '";',
             'alignas(64) static const std::uint8_t kLookup[29440] = {']
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
    lines.append('};\nstatic const std::uint32_t kWindows[][9] = {')
    lines.extend('{' + ','.join(str(v) + 'u' for v in row) + '},' for row in data['windows'])
    lines.append('};\nstatic const Row kRows[] = {')
    lines.extend('{%du,{%s}},' % (r['id'], ','.join(str(v) for v in r['row'])) for r in data['rows'])
    lines.append('};\n} // namespace wh2_k5_data\n')
    raw = '\n'.join(lines).encode('ascii')
    require(len(raw) <= CAP, 'generated header cap')
    return raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--header', required=True, type=Path)
    args = parser.parse_args()
    data = extract(load_report())
    raw = render(data, build_lookup())
    R.write_header(args.header, raw)
    print('K5_DATA bytes=%d sha256=%s history_width_cases=57' % (len(raw), C.sha(raw)))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(type(error).__name__ + ': ' + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
