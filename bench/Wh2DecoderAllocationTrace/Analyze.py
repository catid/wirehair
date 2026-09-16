"""Strict reconstruction of the bounded allocation trace, never timing analysis."""
import argparse
import csv
import gzip
import itertools
import json

HEADER = ('fixture,k,width,tail,policy,family,history,order,position,arm,kind,event,'
          'allocate,array,bytes,pointer,caller_owner,caller_offset,stack,dso_base,phase,arg1,arg2,packet').split(',')
KEY = ('fixture', 'family', 'history', 'order', 'position')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def shape(fixture):
    require(0 <= fixture < 42, 'fixture range')
    k = (3, 5, 8)[fixture // 14]
    variant = fixture // 2 % 7
    width = (2, 2, 64, 256, 257, 1280, 1280)[variant]
    return k, width, 1 if variant in (1, 6) else width, fixture % 2


def disjoint(a, na, b, nb):
    return na <= b-a if a <= b else nb <= a-b


def read_rows(lines):
    for index, line in enumerate(lines):
        require(len(line) <= 1024 and line.endswith('\n'), 'row extent')
        values = next(csv.reader([line]))
        if index == 0:
            require(values == HEADER, 'header')
            continue
        require(len(values) == len(HEADER), 'column count')
        require(all(v and v.isascii() and v.isdigit() for v in values), 'unsigned decimal')
        numbers = list(map(int, values))
        require(all(n < 2**64 for n in numbers), 'integer range')
        yield dict(zip(HEADER, numbers))


def validate(rows, key):
    fixture, family, history, order, position = key
    k, width, tail, policy = shape(fixture)
    arm = (0, 1, 1, 0)[position] ^ order
    require(192 < len(rows) <= 1312, 'bounded complete observation')
    for r in rows:
        require(tuple(r[n] for n in KEY) == key and
                (r['k'], r['width'], r['tail'], r['policy'], r['arm']) ==
                (k, width, tail, policy, arm), 'roster')
        require(r['dso_base'] and r['dso_base'] == rows[0]['dso_base'], 'DSO base')
    events, calls = rows[:192], rows[192:]
    require(len(calls) % 32 == 0, 'call cycles')
    stride = len(calls) // 32
    feeds = stride - 3
    require(k <= feeds <= 32, 'first-success feed bound')
    for i, e in enumerate(events):
        slot = i % 6
        sizes = (296, {3: 104, 5: 136, 8: 200}[k], (k+1)*width, 0, 0, 0)
        require((e['kind'], e['event'], e['allocate'], e['array'], e['bytes'], e['phase']) ==
                (0, i, int(slot < 3), int(slot in (2, 3)), sizes[slot], 0 if slot < 3 else 3),
                'allocation/free sequence')
        require(e['pointer'] and e['stack'] and e['caller_offset'] and
                e['pointer'] + e['bytes'] < 2**64 and e['caller_owner'] in (0, 1) and
                (not e['caller_owner'] or slot == 5), 'event addresses/owner')
        require(not (e['arg1'] or e['arg2'] or e['packet']), 'unused event fields')
        if slot >= 3:
            require(e['pointer'] == events[i-slot+5-slot]['pointer'], 'original-pointer free')
    for i, c in enumerate(calls):
        cycle, slot = divmod(i, stride)
        phase = 0 if slot == 0 else 1 if slot <= feeds else 2 if slot == feeds+1 else 3
        require((c['kind'], c['event'], c['phase']) == (1, i, phase), 'call sequence')
        require(not any(c[n] for n in ('allocate', 'array', 'caller_owner', 'caller_offset', 'stack')),
                'unused call fields')
        require(c['pointer'] == (0 if phase == 0 else events[cycle*6]['pointer']), 'call handle')
        packet = (2**32-1-2*(slot-1) if family else k+slot-1) if phase == 1 else 0
        require(c['packet'] == packet and bool(c['arg1']) == (phase != 3) and
                bool(c['arg2']) == (phase in (0, 2)), 'call arguments')
        require(c['arg1'] + c['bytes'] < 2**64 and c['arg2'] + 8 < 2**64, 'argument range')
        require(c['bytes'] == (32 if phase == 0 else width if phase == 1 else
                              (k-1)*width+tail if phase == 2 else 0), 'call extent')
        if phase != 3:
            if phase in (0, 2):
                require(disjoint(c['arg1'], c['bytes'], c['arg2'], 8), 'public output overlap')
            if phase == 0:
                require(c['arg1'] == calls[0]['arg1'], 'staged descriptor address')
            elif phase == 1:
                require(c['arg1'] == calls[1]['arg1'] + (slot-1)*width, 'staged packet stride')
            else:
                require(c['arg1'] == calls[feeds+1]['arg1'] + cycle*3*k*width, 'recover cycle stride')
            for allocation in events[cycle*6:cycle*6+3]:
                require(disjoint(c['arg1'], c['bytes'], allocation['pointer'], allocation['bytes']),
                        'public/private overlap')
                if c['arg2']:
                    require(disjoint(c['arg2'], 8, allocation['pointer'], allocation['bytes']),
                            'pointer/count output overlap')
    for cycle in range(32):
        live = events[cycle*6:cycle*6+3]
        for a, b in itertools.combinations(live, 2):
            require(disjoint(a['pointer'], a['bytes'], b['pointer'], b['bytes']), 'live overlap')
    return events, calls, feeds


def comparison(a, b):
    ea, ca, _ = a
    eb, cb, _ = b
    return dict(pairs=1,
                same_private_addresses=int([r['pointer'] for r in ea] == [r['pointer'] for r in eb]),
                same_allocator_stack=int([r['stack'] for r in ea] == [r['stack'] for r in eb]),
                same_public_arguments=int([(r['phase'], r['arg1'], r['arg2'], r['bytes'], r['packet']) for r in ca] ==
                                          [(r['phase'], r['arg1'], r['arg2'], r['bytes'], r['packet']) for r in cb]))


def analyze(lines):
    groups = iter(itertools.groupby(read_rows(lines), lambda r: tuple(r[n] for n in KEY)))
    counts, comparisons, placements, dsos, steps = {}, {}, {}, {}, {}
    for fixture in range(42):
        for family in range(2):
            for history in range(3):
                for order in range(2):
                    panel = []
                    for position in range(4):
                        key, source = next(groups, (None, ()))
                        expected = (fixture, family, history, order, position)
                        require(key == expected, 'complete ordered roster')
                        rows = list(itertools.islice(source, 1313))
                        result = validate(rows, expected)
                        panel.append(result)
                        events, calls, feeds = result
                        arm = rows[0]['arm']
                        steps.setdefault((fixture, family), feeds)
                        require(steps[fixture, family] == feeds, 'same first-success count')
                        dsos.setdefault(arm, rows[0]['dso_base'])
                        require(dsos[arm] == rows[0]['dso_base'], 'fixed DSO base')
                        counts['observations'] = counts.get('observations', 0) + 1
                        counts['events'] = counts.get('events', 0) + len(events)
                        counts['calls'] = counts.get('calls', 0) + len(calls)
                        for e in events:
                            if not e['allocate']:
                                continue
                            pk = '{}/history{}/arm{}/allocation{}'.format(fixture, history, arm, e['event'] % 6)
                            places = placements.setdefault(pk, dict(mod64=set(), mod4096=set()))
                            places['mod64'].add(e['pointer'] % 64)
                            places['mod4096'].add(e['pointer'] % 4096)
                    for name, x, y in (('B/B', 0 if not order else 1, 3 if not order else 2),
                                       ('C/C', 1 if not order else 0, 2 if not order else 3),
                                       ('C/B', 0, 1), ('C/B', 2, 3)):
                        value = comparison(panel[x], panel[y])
                        for label in (name, '{}/history{}/{}'.format(fixture, history, name)):
                            total = comparisons.setdefault(label, {n: 0 for n in value})
                            for n in value:
                                total[n] += value[n]
    require(next(groups, None) is None, 'trailing observations')
    require(counts['observations'] == 2016 and counts['events'] == 387072, 'complete trace')
    require(len(dsos) == 2 and dsos[0] != dsos[1], 'distinct DSO mappings')
    return dict(scope='Non-timing small-decoder observation; not historical or speed attribution',
                counts=counts, comparisons=comparisons, dso_bases=dsos,
                placements={k: {n: sorted(v) for n, v in values.items()} for k, values in placements.items()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace')
    args = parser.parse_args()
    opener = gzip.open if args.trace.endswith('.gz') else open
    with opener(args.trace, 'rt', encoding='ascii', newline='') as stream:
        print(json.dumps(analyze(stream), sort_keys=True, indent=2))
