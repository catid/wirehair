#!/usr/bin/env python3
"""Frozen ordinary K3 recovery, both source policies; retained cohort only."""
import argparse
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import struct
import sys
import time

HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


U = sibling('_k3_recovery_build_io', 'Wh2K3OrdinaryCostR0.py')
D = sibling('_k3_recovery_retained_data', 'Wh2K3NativeDataR0.py')
A = U.A
ROOT = U.ROOT
PROTOCOL = 'wirehair.wh2.k3-ordinary-recovery-r0'
OUTPUT = Path('/var/tmp/wh2-k3-ordinary-recovery-r0')
MODES = ('native', 'scalar', 'asan')
WIDTHS = (2, 64, 1280)
ARMS = ('certified_wh2', 'ordinary_independent', 'wh1', 'ordinary_borrowed')
ENV_KEYS = U.ENV_KEYS + ('ASAN_OPTIONS', 'UBSAN_OPTIONS')
SANITIZERS = dict(ASAN_OPTIONS='detect_leaks=1:detect_stack_use_after_return=1',
                  UBSAN_OPTIONS='halt_on_error=1')
SOURCES = ('bench/Wh2K3OrdinaryRecoveryR0.cpp', 'bench/Wh2K3OrdinaryRecoveryR0.py',
           'bench/test_Wh2K3OrdinaryRecoveryR0.py', 'bench/Wh2K3OrdinaryCostR0.py',
           'bench/Wh2AlignedIntermediateCostR0.py', 'bench/Wh2K3NativeDataR0.py',
           'bench/Wh2NoncommutingRadixRunR0.py', 'SMALL_WIRE_PROFILES.md',
           'V2_WIRE_PROFILE.md', 'WH2_BORROWED_SOURCE_API.md',
           'bench/test_Wh2K3OrdinaryCostR0.py')


@lru_cache(maxsize=1)
def retained():
    report = D.load_report()
    data = D.extract(report)
    return report, data


@lru_cache(maxsize=1)
def roster():
    report, data = retained()
    result = []
    for group, name in enumerate(('fresh', 'hard')):
        for index, row in enumerate(report[name]):
            result.append(dict(group=group, index=index, B=row['B'], ids=row['ids'],
                               schedule=row['schedule'], root=row['root']))
    for index, row in enumerate(data['history']):
        for width in row['widths']:
            result.append(dict(group=2, index=index, B=width, ids=row['ids']))
    A.exact(len(result), 6269, 'retained trace/width roster')
    return result


@lru_cache(maxsize=1)
def coefficient_rows():
    return {r['id']: r['row'] for r in retained()[1]['rows']}


@lru_cache(maxsize=256)
def products(coefficient):
    # Polynomial arithmetic only; never import candidate payload tables.
    return bytes(D.multiply(coefficient, b) for b in range(256))


@lru_cache(maxsize=3)
def message(width):
    return bytes((37*i+i//11) % 256 for i in range(3*width))


@lru_cache(maxsize=8192)
def candidate_packet(width, packet_id):
    row = coefficient_rows()[packet_id]
    a, b, c = (products(v) for v in row)
    source = message(width)
    packet = bytes(a[source[j]] ^ b[source[width+j]] ^ c[source[2*width+j]] for j in range(width))
    return A.sha(packet)


@lru_cache(maxsize=8192)
def candidate_first(ids):
    rows = coefficient_rows()
    return next((n for n in range(3, len(ids)+1) if U.scalar_rank([rows[i] for i in ids[:n]]) == 3), 0)


def valid_hex(value, size):
    A.require(type(value) is str and len(value) == 2*size and re.fullmatch('[0-9a-f]+', value) is not None,
              'hex bytes')


def summarize(records, expected):
    A.exact(len(records), len(expected), 'summary whole cohort')
    groups = []
    for group, name in enumerate(('retained_fresh', 'hard', 'historical_original_width')):
        selected = [(r, e) for r, e in zip(records, expected) if e['group'] == group]
        curves = [[sum(a['first'] == 0 or a['first'] > 3+oh for r, _ in selected for a in [r['arms'][arm]])
                   for oh in range(5)] for arm in range(4)] if group != 2 else None
        pairs = []
        for candidate, control in ((1,0), (1,2), (3,0), (3,2)):
            pair = dict(control=ARMS[control], candidate=ARMS[candidate], fixed=[], introduced=[])
            for oh in range(5):
                fixed, introduced = [], []
                for r, e in selected:
                    if group == 2 and 3+oh > len(e['ids']):
                        continue
                    x, y = r['arms'][control]['first'], r['arms'][candidate]['first']
                    control_failed, candidate_failed = x == 0 or x > 3+oh, y == 0 or y > 3+oh
                    key = [e['index'], e['B']]
                    if control_failed and not candidate_failed: fixed.append(key)
                    if candidate_failed and not control_failed: introduced.append(key)
                pair['fixed'].append(fixed); pair['introduced'].append(introduced)
            pairs.append(pair)
        groups.append(dict(group=name, cases=len(selected), failure_counts_oh0_to_4=curves,
                           unresolved=[sum(r['arms'][a]['first'] == 0 for r, _ in selected) for a in range(4)],
                           paired=pairs))
    cells = []
    schedules = sorted({e['schedule'] for e in expected if e['group'] == 0})
    A.exact(len(schedules), 4, 'four retained loss schedules')
    for width in WIDTHS:
        for schedule in schedules:
            selected = [r for r, e in zip(records, expected) if e['group'] == 0 and e['B'] == width and e['schedule'] == schedule]
            A.exact(len(selected), 512, 'whole retained cell')
            curves = [[sum(r['arms'][a]['first'] == 0 or r['arms'][a]['first'] > 3+oh for r in selected)
                       for oh in range(5)] for a in range(4)]
            cells.append(dict(B=width, schedule=schedule, cases=512, failure_counts_oh0_to_4=curves))
    fresh = groups[0]['failure_counts_oh0_to_4']
    wins = {ARMS[c]: fresh[c][0] < fresh[0][0] and fresh[c][0] < fresh[2][0] for c in (1,3)}
    return dict(outcome='PASS', arms=list(ARMS), groups=groups, cells=cells,
                comparative_oh0_win=all(wins.values()), policy_oh0_wins=wins,
                fresh_sample=False, speed_claimed=False, all_K_claimed=False,
                production_promotion_claimed=False)


def verify(raw, claim, mode):
    A.require(0 < len(raw) <= U.RAW_CAP and raw.endswith(b'\n'), 'bounded complete stream')
    rows = [A.decode(line) for line in raw.splitlines()]
    expected = roster()
    A.exact(len(rows), len(expected)+2, 'whole recovery cohort')
    header, records, footer = rows[0], rows[1:-1], rows[-1]
    A.exact(header, dict(type='header', protocol=PROTOCOL, claim=claim, backend=mode,
                       retained_raw_sha256=D.RAW_SHA, features=[0]*4 if mode == 'scalar' else [1]*4,
                       sources=[A.sha(message(b)) for b in WIDTHS]), 'exact replay header')
    profiles = {}
    for row, e in zip(records, expected):
        A.exact(set(row), {'type', 'group', 'index', 'B', 'ids', 'arms'}, 'record schema')
        A.exact({k: row[k] for k in ('group', 'index', 'B', 'ids')},
                {k: e[k] for k in ('group', 'index', 'B', 'ids')}, 'retained exact chronology')
        A.exact(row['type'], 'record', 'record type')
        A.exact(len(row['arms']), 4, 'four actual API routes')
        width, ids = e['B'], e['ids']
        for arm, r in enumerate(row['arms']):
            A.exact(set(r), {'profile', 'packets', 'feed', 'first', 'recoveries', 'checked'}, 'arm schema')
            valid_hex(r['profile'], 32)
            A.exact(r['profile'], profiles.setdefault((width, arm), r['profile']), 'stable full descriptor')
            profile = bytes.fromhex(r['profile'])
            if arm == 0:
                A.exact(profile[:16], struct.pack('<4sHHQ', b'WHV2', 1, 32, 0x4b295bbb47f4f9c9),
                        'explicit certified WH2 descriptor')
                A.exact(profile[16:28], struct.pack('<QI', 3*width, width), 'certified WH2 dimensions')
                A.exact(profile[29:], bytes(3), 'certified reserved bytes')
            elif arm in (1,3):
                A.exact(profile, struct.pack('<4sHHQQII', b'WHV2', 1, 32, 0x67c1043ecaa9e184,
                                            3*width, width, 0), 'ordinary WHV2 K3 identity')
            else:
                A.exact(profile, b'\0'*32, 'WH1 has no descriptor')
            A.exact(len(r['packets']), len(ids), 'all retained packets')
            for packet in r['packets']:
                valid_hex(packet, 32)
            first = A.integer(r['first'], 0, len(ids))
            A.require(first == 0 or first >= 3, 'possible first success')
            A.exact(r['feed'], [1]*(first-1)+[0] if first else [1]*len(ids), 'all prefix feed statuses')
            A.exact(r['recoveries'], 2 if first else 0, 'two byte-exact recoveries')
            A.exact(r['checked'], True, 'payload/guard validation')
            if arm in (1,3):
                A.exact(r['packets'], [candidate_packet(width, i) for i in ids], 'independent K3 payload hashes')
                A.exact(first, candidate_first(tuple(ids)), 'independent first-success rank')
                A.require(first != 0, 'retained candidate must recover')
                if e['group'] != 2:
                    A.exact(first, 3, 'retained zero-overhead candidate')
        A.exact(row['arms'][1], row['arms'][3], 'policy-independent complete result')
    A.exact(footer, dict(type='footer', records=6269, checked=True), 'complete footer')
    return summarize(records, expected), records


def build(mode, output):
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT not in output.parents and output != ROOT and not output.exists() and not output.is_symlink(),
              'fresh external build')
    archives, dependencies = U.qualified_inputs(mode)
    output.mkdir(mode=0o700)
    U.command([sys.executable, '-I', '-B', '-S', HERE/'Wh2K3NativeDataR0.py', '--header', output/'Wh2K3NativeData.inc'])
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_K3_ORDINARY_RECOVERY_BACKEND='+str(MODES.index(mode)),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer'] if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar': flags.append('-DANDROID=1')
    dependencies.update(ROOT/s for s in SOURCES)
    commands, objects = [], []
    for source in (ROOT/SOURCES[0], HERE/'Wh2FrozenTrace.cpp'):
        obj, dep = output/(source.stem+'.o'), output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)]
        U.command(args); commands.append(args); objects.append(obj)
        dependencies.update(Path(s).resolve(strict=True) for s in shlex.split(dep.read_text().replace('\\\n', '').split(': ', 1)[1]))
    executable = output/'recovery_worker'
    args = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan': args.append('-fsanitize=address,undefined')
    args += list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    U.command(args); commands.append(args)
    names = [line.split()[-1] for line in U.command(['/usr/bin/nm', '-g', '--defined-only', executable]).decode().splitlines() if line.split()]
    A.require(names.count('GF256Ctx') == names.count('gf256_init_') == 1, 'single GF runtime')
    A.require(not any(n.startswith('wh2_small_') for n in names), 'no prototype codec')
    for prefix in ('wirehair_', 'wirehair_v2_', 'wirehair_small_'):
        for name in ('encode', 'decode', 'free'):
            A.exact(names.count(prefix+name), 1, 'actual installed API')
    for name in ('c++', 'cc', 'as', 'ld', 'nm', 'ar', 'make'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for name in ('cc1', 'cc1plus', 'collect2'):
        dependencies.add(Path(U.command(['/usr/bin/c++', '-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    for target in (executable, Path(sys.executable)):
        dependencies.update(Path(w).resolve(strict=True) for w in U.command(['/usr/bin/ldd', target]).decode().split() if w.startswith('/'))
    dependencies.add(Path(sys.executable).resolve(strict=True))
    neutral = U.command([executable, '--neutral'])
    A.require(neutral.startswith(b'PASS 12 neutral cases,'), 'neutral public API replay')
    A.publish(output/'neutral.txt', neutral)
    manifest = dict(protocol=PROTOCOL, mode=mode, commands=commands,
                    inputs=[U.pin(p) for p in sorted(dependencies)],
                    artifacts=[U.pin(p) for p in sorted(output.iterdir())])
    A.publish(output/'manifest.json', A.canonical(manifest))
    print(json.dumps(dict(mode=mode, neutral=neutral.decode().strip(), scientific_launch=False)))


def current(receipt):
    A.exact(set(receipt), {'protocol', 'head', 'executables', 'environment', 'pins'}, 'receipt schema')
    A.exact(receipt['protocol'], PROTOCOL, 'receipt protocol')
    A.exact(receipt['head'], U.command(['git', 'rev-parse', 'HEAD']).decode().strip(), 'source HEAD')
    A.exact(receipt['environment'], {k: os.environ.get(k) for k in ENV_KEYS}, 'runtime environment')
    A.exact(receipt['environment'], dict({k: None for k in U.ENV_KEYS}, **SANITIZERS), 'full sanitizer policy')
    A.exact(set(receipt['executables']), set(MODES), 'all qualified executables')
    declared = {p['path']: p for p in receipt['pins']}
    A.exact(len(declared), len(receipt['pins']), 'unique receipt pins')
    for pin in receipt['pins']:
        A.exact(U.pin(pin['path']), pin, 'unchanged input')
    closure = {}
    for mode in MODES:
        executable = Path(receipt['executables'][mode])
        A.require(executable.is_absolute() and executable.name == 'recovery_worker' and
                  executable.parent.name == mode and str(executable) in declared, 'pinned executable binding')
        path = executable.parent/'manifest.json'
        A.require(str(path) in declared, 'pinned build manifest')
        manifest = A.decode(A.read_regular(path, 1024*1024))
        A.exact((manifest['protocol'], manifest['mode']), (PROTOCOL, mode), 'bound backend build')
        A.require(str(executable) in {p['path'] for p in manifest['artifacts']}, 'built executable artifact')
        for pin in manifest['inputs']+manifest['artifacts']+[declared[str(path)]]:
            if pin['path'] in closure: A.exact(closure[pin['path']], pin, 'shared manifest input')
            closure[pin['path']] = pin
    A.exact(receipt['pins'], sorted(closure.values(), key=lambda p: p['path']), 'complete manifest closure')
    retained()


def receipt(build_dir):
    pins, executables = {}, {}
    for mode in MODES:
        path = build_dir/mode/'manifest.json'
        manifest = A.decode(A.read_regular(path, 1024*1024))
        A.exact((manifest['protocol'], manifest['mode']), (PROTOCOL, mode), 'qualified backend')
        for pin in manifest['inputs']+manifest['artifacts']+[U.pin(path)]:
            if pin['path'] in pins: A.exact(pins[pin['path']], pin, 'shared input identity')
            pins[pin['path']] = pin
        executables[mode] = str(build_dir/mode/'recovery_worker')
    for pin in pins.values():
        A.exact(U.pin(pin['path']), pin, 'current qualification')
        path = Path(pin['path'])
        if ROOT in path.parents:
            A.exact(path.read_bytes(), U.command(['git', 'cat-file', 'blob', 'HEAD:'+str(path.relative_to(ROOT))]), 'committed input')
    environment = {k: os.environ.get(k) for k in ENV_KEYS}
    A.exact(environment, dict({k: None for k in U.ENV_KEYS}, **SANITIZERS), 'ordinary allocator and full sanitizer policy')
    result = dict(protocol=PROTOCOL, head=U.command(['git', 'rev-parse', 'HEAD']).decode().strip(),
                  executables=executables, environment=environment,
                  pins=sorted(pins.values(), key=lambda p: p['path']))
    current(result)
    return result


def run(path):
    begin = time.monotonic()
    frozen = A.read_regular(path, 1024*1024); claim = A.decode(frozen)
    A.exact(frozen, A.canonical(claim), 'canonical receipt'); current(claim)
    OUTPUT.mkdir(mode=0o700); A.publish(OUTPUT/'CLAIM.json', frozen)
    analysis, reference = dict(outcome='INVALID'), None
    try:
        for mode in MODES:
            raw, error, code, failure = U.capture(claim['executables'][mode], A.sha(frozen),
                min(begin+300, time.monotonic()+90), [OUTPUT/(mode+'.raw.jsonl'), OUTPUT/(mode+'.stderr.txt')])
            A.require(failure is None and code == 0 and error == b'', mode+' worker/observer: '+str(failure))
            summary, records = verify(raw, A.sha(frozen), mode)
            if reference is None:
                reference = records; analysis = summary
            else:
                A.exact(records, reference, 'cross-backend complete byte/status parity')
                A.exact(summary, analysis, 'cross-backend summary parity')
        current(claim)
        A.require(time.monotonic()-begin < 300, 'whole replay deadline')
    except Exception as error:
        analysis = dict(outcome='INVALID', failure=str(error))
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'analysis.json', A.canonical(analysis))
    members = [U.pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members) < 256*1024**2-65536, 'bundle cap')
    A.publish(OUTPUT/'COMPLETE.json', A.canonical(dict(protocol=PROTOCOL, outcome=analysis['outcome'], files=members)))
    print(json.dumps({k: v for k, v in analysis.items() if k not in ('groups', 'cells')}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('mode', choices=MODES); b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    args = parser.parse_args()
    if args.command == 'build': build(args.mode, args.output)
    elif args.command == 'receipt': A.publish(args.output, A.canonical(receipt(args.build_dir.resolve(strict=True))))
    else: run(args.receipt)
