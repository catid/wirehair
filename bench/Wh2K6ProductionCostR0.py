#!/usr/bin/env python3
"""Fixed .81 actual-library K6 lifecycle gate; no old experiment is rerun or rescored."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import selectors
import re
import struct
from functools import lru_cache
import shlex
import subprocess
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = sibling('k6_production_cost_common', 'Wh2AlignedIntermediateCostR0.py')
ROOT = A.ROOT
PROTOCOL = 'wirehair.wh2.k6-production-cost-r0'
OUTPUT = Path('/var/tmp/wh2-k6-production-cost-r0')
QUALIFIED = Path('/tmp/wh2-k6-production.JQkPgH')
PRIOR = Path('/var/tmp/wh2-thue-decoder-cost-r0/raw.jsonl')
PRIOR_SHA = '149096c9f4caa6f5303cc8eb9062bd776b10702f1bda3c55a7241066ef19494f'
ARCHIVES = dict(
    native='ef289be958fc945aead207c9838237dce684290f8574a50647bbf49a65410576',
    scalar='f0829b419cfd858cf57acc0bff0e1d299d400bd87f9dd9a9db08f3fe10aab7fc',
    asan='9cb2c6af58acd6b9222707aabd5c4d132c9e987f48275e469f6a7e2aaa2d5478')
PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (2, 1))
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0)
WIDTHS = (2, 64, 1280)
CALLBACKS, BATCH = 19440, 128
RAW_CAP, ERR_CAP = 64*1024**2, 65536
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
NEW = ('bench/Wh2K6ProductionCostR0.cpp', 'bench/Wh2K6ProductionCostR0.py',
       'bench/test_Wh2K6ProductionCostR0.py')


def pin(path):
    p = Path(path)
    owned = any(base == p or base in p.parents for base in (ROOT, Path('/tmp'), Path('/var/tmp')))
    return A.pin(p, installed=not owned)


def command(args):
    return subprocess.check_output(list(map(str, args)), timeout=60, cwd=ROOT)


def roster():
    index = 0
    for r in range(12):
        for s in range(2):
            for ws in range(3):
                for ms in range(3):
                    for cs in range(5):
                        order, width, metric = (r+s) % 2, (r+s+ws) % 3, (r+s+ws+ms) % 3
                        comparison = (2*r+s+ws+metric+cs) % 5
                        for p in range(18):
                            bin_number = r+12*(r % 4) if p < 2 else (r+6*((p-2)//8)) % 12+12*(((p-2) % 8)//2)
                            yield [index, r, order, width, metric, comparison, p,
                                   PAIRS[comparison][SIDES[p] ^ order], (2*bin_number+1)*1000000//96]
                            index += 1


def statistics(records):
    A.exact(len(records), CALLBACKS, 'complete statistical cohort')
    groups = {}
    for start in range(0, CALLBACKS, 18):
        panel, contrasts = records[start:start+18], []
        for j in range(8):
            values = {}
            for row in panel[2+2*j:4+2*j]:
                c, clocks = row['coordinate'], row['observation']['clocks']
                values[SIDES[c[6]] ^ c[2]] = clocks[3]-clocks[2]
            A.require(set(values) == {0, 1} and min(values.values()) > 0, 'paired positive durations')
            contrasts.append(math.log(values[1])-math.log(values[0]))
        c = panel[0]['coordinate']
        groups.setdefault((c[3], c[4], c[5], c[2]), []).append(math.fsum(contrasts)/8)
    A.exact(len(groups), 90, 'separate full-cost cells')
    results, failed_controls, failed_treatments = [], [], []
    bound = math.log1p(.02)
    for key, values in sorted(groups.items()):
        width, metric, comparison, order = key
        estimate = A.confidence(values)
        item = dict(width=WIDTHS[width], metric=metric, comparison=comparison, order=order,
                    estimate=estimate, replicate_logs=values)
        if comparison < 3:
            passed = -bound < estimate['lower95_log'] and estimate['upper95_log'] < bound
            item['control_pass'] = passed
            if not passed:
                failed_controls.append(list(key))
        else:
            item['treatment_pass'] = estimate['upper95_log'] < 0
            item['upper_ratio_limit'] = 1.0
            if not item['treatment_pass']:
                failed_treatments.append(list(key))
        results.append(item)
    return dict(outcome='CONTROL_FAIL' if failed_controls else 'FAIL' if failed_treatments else 'PASS',
                statistics=results, failed_controls=failed_controls, failed_treatments=failed_treatments,
                speed_qualified=not failed_controls and not failed_treatments,
                all_K_claimed=False, recovery_rate_claimed=False, production_promotion_claimed=False)



def clocks(observation, previous):
    A.exact(set(observation), {'clocks', 'before', 'after'}, 'observation keys')
    c = observation['clocks']
    A.require(len(c) == 6, 'clock shape')
    for value in c:
        A.integer(value)
    A.require(c[0] <= c[2] < c[3] <= c[5] and c[1] <= c[4] and c[4]-c[1] <= c[5]-c[0], 'clock ordering')
    for key in ('before', 'after'):
        A.exact(len(observation[key]), 4, 'counter shape')
        for value in observation[key]:
            A.integer(value)
    A.require(all(a <= b for a, b in zip(observation['before'], observation['after'])), 'counter ordering')
    if previous is not None:
        A.require(previous['clocks'][5] <= c[0] and previous['clocks'][4] <= c[1] and
                  all(a <= b for a, b in zip(previous['after'], observation['before'])), 'cross-record ordering')


def multiply(a, b):
    product = 0
    for bit in range(8):
        if b & (1 << bit):
            product ^= a << bit
    for bit in range(14, 7, -1):
        if product & (1 << bit):
            product ^= 0x14d << (bit-8)
    return product


@lru_cache(maxsize=1)
def selected_rows():
    # Independent polynomial arithmetic and binary prefix decomposition;
    # never consult the candidate's embedded packed table.
    def product(a, b):
        out = [0]*36
        for r in range(6):
            for c in range(6):
                for k in range(6):
                    out[6*r+c] ^= multiply(a[6*r+k], b[6*k+c])
        return out
    pair = []
    for feedback in ((124,127,152,84,241,63), (125,127,152,84,241,63)):
        matrix = [0]*36
        for i in range(5):
            matrix[6*(i+1)+i] = 1
        for i, value in enumerate(feedback):
            matrix[6*i+5] = value
        pair.append(matrix)
    levels = [pair]
    for _ in range(31):
        a, b = levels[-1]
        levels.append([product(a,b), product(b,a)])
    rows = []
    for packet_id in tuple(range(12))+tuple(0xffffffff-2*j for j in range(6)):
        vector = [1,0,0,0,0,0]
        for bit in range(32):
            if packet_id & (1 << bit):
                matrix = levels[bit][bin(packet_id >> (bit+1)).count('1') % 2]
                out = [0]*6
                for r in range(6):
                    for c in range(6):
                        out[r] ^= multiply(matrix[6*r+c], vector[c])
                vector = out
        rows.append(tuple(vector))
    return tuple(rows)


@lru_cache(maxsize=3)
def candidate_fixture(width):
    A.require(width in WIDTHS, 'frozen width')
    source = bytes((37*i+i//11) % 256 for i in range(6*width))
    packets = bytearray()
    for row in selected_rows():
        for j in range(width):
            value = 0
            for k in range(6):
                value ^= multiply(row[k], source[k*width+j])
            packets.append(value)
    profile = struct.pack('<4sHHQQII', b'WHK6', 1, 32, 0x5748324b36544d31, 6*width, width, 0)
    return profile.hex(), packets.hex()


def prior_header():
    raw = A.read_regular(PRIOR, 64*1024**2)
    A.exact(A.sha(raw), PRIOR_SHA, 'immutable old correctness corpus; no old timing analysis')
    return A.decode(raw.splitlines()[0])


def verify_fixtures(header, old):
    A.exact(len(header['fixtures']), 3, 'fixture count')
    A.exact(len(old['fixtures']), 3, 'prior fixture count')
    for f, prior, width in zip(header['fixtures'], old['fixtures'], WIDTHS):
        A.exact(set(f), {'width', 'source', 'arms'}, 'fixture fields')
        A.exact(f['width'], width, 'fixture width')
        A.exact(f['source'], bytes((37*i+i//11) % 256 for i in range(6*width)).hex(), 'fixed source')
        A.exact(len(f['arms']), 3, 'three fixtures')
        public_profile = next(h['profile_hex'] for h in prior['handles'] if h['arm'] == 3)
        candidate_profile, candidate_packets = candidate_fixture(width)
        for a, arm in enumerate(f['arms']):
            A.exact(set(arm), {'profile', 'packets', 'steps'}, 'arm fields')
            A.exact(arm['profile'], (public_profile, candidate_profile, '00'*32)[a], 'distinct descriptor')
            A.exact(arm['packets'], prior['packets_hex'][(3,1,2)[a]], 'same qualified packet corpus')
            A.exact(len(arm['steps']), 2, 'decoder families')
            prior_arm = next(h for h in prior['handles'] if h['arm'] == (3,1,2)[a])
            for family, count in enumerate(arm['steps']):
                A.integer(count, 6, 12)
                A.exact(count, len(prior_arm['decode'][family]['feeds']), 'same first-success endpoint')
        A.exact(f['arms'][1]['packets'], candidate_packets, 'independent selected-pair payload oracle')


def verify(raw, claim, old):
    A.require(0 < len(raw) <= RAW_CAP and raw.endswith(b'\n'), 'bounded complete stream')
    rows = [A.decode(line) for line in raw.splitlines()]
    A.exact(len(rows), CALLBACKS+2, 'whole raw cohort')
    header, footer = rows[0], rows[-1]
    A.exact((header['type'], header['protocol'], header['claim'], header['batch']),
            ('header', PROTOCOL, claim, BATCH), 'header identity')
    A.exact(header['identity_hex'], old['identity_before']['canonical_hex'], 'frozen captured CPU identity')
    verify_fixtures(header, old)
    previous = header['prelude']; clocks(previous, None)
    total_work = 0
    for row, coordinate in zip(rows[1:-1], roster()):
        A.exact(set(row), {'type', 'coordinate', 'ready', 'target', 'wait', 'observation', 'counts',
                          'addresses', 'address_count', 'complete', 'checked'}, 'record fields')
        A.exact(row['type'], 'record', 'record type')
        A.exact(row['coordinate'], coordinate, 'fixed chronology')
        index, rep, order, width, metric, comparison, position, arm, q = coordinate
        del index, rep, comparison, position
        A.integer(row['ready']); A.integer(row['target'])
        A.require(previous['clocks'][5] <= row['ready'] and row['target'] == row['ready']+q, 'relative target')
        A.exact(len(row['wait']), 4, 'wait shape')
        w0, c0, w1, c1 = [A.integer(v) for v in row['wait']]
        A.require(row['ready'] <= w0 <= w1 and w1 >= row['target'] and c0 <= c1 and
                  previous['clocks'][4] <= c0, 'wait chronology')
        clocks(row['observation'], previous)
        observed = row['observation']['clocks']
        A.require(observed[0] >= w1 and observed[1] >= c1 and observed[2] >= row['target'], 'retained delayed start')
        previous = row['observation']; total_work += observed[3]-observed[2]
        steps = header['fixtures'][width]['arms'][arm]['steps'][metric-1] if metric else 0
        expected = [0 if metric else BATCH, 0 if metric else 18*BATCH,
                    BATCH if metric else 0, BATCH*steps, BATCH if metric else 0, BATCH]
        A.exact(row['counts'], expected, 'whole lifecycle ledger')
        A.exact(row['address_count'], BATCH, 'all handles retained')
        A.exact(len(row['addresses']), BATCH, 'address roster')
        for address in row['addresses']:
            A.integer(address, 1)
        A.exact(row['complete'], True, 'every call completed')
        A.exact(row['checked'], True, 'every output checked')
    A.exact(footer, dict(type='footer', complete=True, records=CALLBACKS, work_ns=total_work), 'complete footer')
    A.require(total_work <= 80000000000, 'inner work cap')
    return statistics(rows[1:-1])


def qualified_inputs(mode):
    base = QUALIFIED/mode
    archives = [base/'libwirehair.a']
    A.exact(pin(archives[0])['sha256'], ARCHIVES[mode], 'neutral-qualified actual library')
    dependencies = set(archives)
    commands = A.decode(A.read_regular(base/'compile_commands.json', 1024*1024))
    built = {str(Path(c['output'])): c for c in commands if (base/c['output']).is_file()}
    prefix = 'CMakeFiles/wirehair_objects.dir/' if mode == 'native' else 'CMakeFiles/wirehair.dir/'
    producing = {output for output in built if output.startswith(prefix)}
    A.require(producing, 'actual producing library target')
    # Validate every archive member against its actual compiled object, not
    # merely a nearby source list. Offline generator commands are not built.
    for archive in archives:
        members = command(['/usr/bin/ar', 't', archive]).decode().splitlines()
        A.require(len(members) == len(set(members)) and members, 'unique archive members')
        A.exact(set(members), {Path(output).name for output in producing}, 'complete actual library members')
        for member in members:
            matches = [base/output for output in producing if Path(output).name == member]
            A.exact(len(matches), 1, 'one actual object for archive member')
            A.exact(A.sha(command(['/usr/bin/ar', 'p', archive, member])),
                    pin(matches[0])['sha256'], 'archive/object bytes')
    dep_text = command(['/usr/bin/ninja', '-C', base, '-t', 'deps']).decode()
    seen = set()
    for block in dep_text.strip().split('\n\n'):
        lines = block.splitlines()
        if not lines:
            continue
        match = re.fullmatch(r'(.+): #deps ([0-9]+), deps mtime [0-9]+ \(VALID\)', lines[0])
        A.require(match is not None, 'valid recorded compile dependencies')
        output, count = match[1], int(match[2])
        A.require(output in built and output not in seen, 'actual dependency object')
        seen.add(output)
        paths = [Path(line.strip()).resolve(strict=True) for line in lines[1:]]
        A.exact(len(paths), count, 'whole dependency record')
        A.require(Path(built[output]['file']).resolve(strict=True) in paths, 'actual TU in dependency closure')
        dependencies.update(paths)
        dependencies.add(base/output)
    A.exact(seen, set(built), 'all compiled objects have dependency closure')
    for relative in ('compile_commands.json', 'CMakeCache.txt', 'build.ninja', 'CMakeFiles/rules.ninja',
                     'k6_codec_test', 'k6_c_consumer', 'k6_payload_test', 'Testing/Temporary/LastTest.log'):
        dependencies.add(base/relative)
    dependencies.update((ROOT/'CMakeLists.txt', ROOT/'abi/wirehair.map', ROOT/'K6_WIRE_PROFILE.md', PRIOR))
    prior_header()
    return archives, dependencies


def build(mode, output):
    output = output.parent.resolve(strict=True)/output.name
    A.require(output.is_absolute() and ROOT not in output.parents and output != ROOT and
              not output.exists() and not output.is_symlink(), 'fresh external build')
    archives, dependencies = qualified_inputs(mode)
    output.mkdir(mode=0o700)
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_K6_COST_NEUTRAL='+str(int(mode != 'native')),
             '-I'+str(ROOT), '-I'+str(ROOT/'bench'), '-I'+str(ROOT/'include')]
    flags += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer'] if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags += ['-DANDROID=1']
    sources = [ROOT/NEW[0]]+[ROOT/'bench'/name for name in
               ('Wh2FrozenTrace.cpp', 'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp')]
    commands, objects = [], []
    dependencies.update(ROOT/n for n in NEW)
    dependencies.add(ROOT/'bench/Wh2AlignedIntermediateCostR0.py')
    for source in sources:
        obj = output/(source.stem+'.o'); dep = output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)]
        command(args); commands.append(args); objects.append(obj)
        dependencies.update(Path(s).resolve(strict=True) for s in shlex.split(dep.read_text().replace('\\\n', '').split(': ', 1)[1]))
    executable = output/'cost_worker'
    args = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan':
        args += ['-fsanitize=address,undefined']
    args += list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    command(args); commands.append(args)
    symbols = command(['/usr/bin/nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    A.require(names.count('GF256Ctx') == names.count('gf256_init_') == 1, 'single GF runtime')
    A.require(all(names.count(prefix+name) == 1 for prefix in ('wirehair_v2_', 'wirehair_k6_') for name in
                  ('encode', 'decode', 'free')), 'actual public WH2 and K6 APIs')
    A.require(all(names.count(n) == 1 for n in ('wirehair_encode', 'wirehair_decode', 'wirehair_free')), 'actual WH1 API')
    text = command(['/usr/bin/nm', '-C', executable]).decode().splitlines()
    A.require(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line for line in text) == 1, 'one common WORK body')
    for name in ('c++', 'cc', 'as', 'ld', 'nm', 'ar', 'ninja'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for compiler in ('cc1', 'cc1plus', 'collect2'):
        dependencies.add(Path(command(['/usr/bin/c++', '-print-prog-name='+compiler]).decode().strip()).resolve(strict=True))
    for target in (executable, Path(sys.executable)):
        linked = command(['/usr/bin/ldd', target]).decode()
        dependencies.update(Path(word).resolve(strict=True) for word in linked.split() if word.startswith('/'))
    dependencies.add(Path(sys.executable).resolve(strict=True))
    manifest = dict(protocol=PROTOCOL, mode=mode, commands=commands,
                    inputs=[pin(p) for p in sorted(dependencies)], artifacts=[pin(p) for p in sorted(output.iterdir())])
    A.publish(output/'manifest.json', A.canonical(manifest))
    print(json.dumps(dict(mode=mode, executable=str(executable), scientific_launch=False)))


def current(receipt):
    A.exact(receipt['protocol'], PROTOCOL, 'receipt protocol')
    A.exact(receipt['environment'], {k: os.environ.get(k) for k in ENV_KEYS}, 'allocator environment')
    A.exact(command(['git', 'rev-parse', 'HEAD']).decode().strip(), receipt['head'], 'source HEAD')
    for p in receipt['pins']:
        A.exact(pin(Path(p['path'])), p, 'unchanged receipt input')


def receipt(build_dir):
    manifest = A.decode(A.read_regular(build_dir/'manifest.json', 1024*1024))
    A.exact((manifest['protocol'], manifest['mode']), (PROTOCOL, 'native'), 'native full cost build')
    pins = {p['path']: p for p in manifest['inputs']+manifest['artifacts']}
    pins[str(build_dir/'manifest.json')] = pin(build_dir/'manifest.json')
    for p in pins.values():
        A.exact(pin(Path(p['path'])), p, 'current build/input')
        path = Path(p['path'])
        if ROOT in path.parents:
            A.exact(path.read_bytes(), command(['git', 'cat-file', 'blob', 'HEAD:'+str(path.relative_to(ROOT))]), 'committed source')
    result = dict(protocol=PROTOCOL, head=command(['git', 'rev-parse', 'HEAD']).decode().strip(),
                  executable=str(build_dir/'cost_worker'), environment={k: os.environ.get(k) for k in ENV_KEYS},
                  pins=sorted(pins.values(), key=lambda p: p['path']))
    A.require(all(v is None for v in result['environment'].values()), 'ordinary allocator policy')
    current(result)
    return result


def capture(executable, claim, deadline, spools):
    buffers = [bytearray(), bytearray()]; files = []; child = None; failure = None
    selector = selectors.DefaultSelector()
    try:
        for p in spools:
            files.append(os.open(str(p), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600))
        child = subprocess.Popen([str(executable), '--worker', claim], stdin=subprocess.DEVNULL,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        for i, stream in enumerate((child.stdout, child.stderr)):
            os.set_blocking(stream.fileno(), False); selector.register(stream, selectors.EVENT_READ, i)
        while selector.get_map():
            for key, _ in selector.select(min(.05, A.time_left(deadline))):
                block = os.read(key.fileobj.fileno(), 65536)
                if not block:
                    selector.unregister(key.fileobj); continue
                i = key.data; available = (RAW_CAP, ERR_CAP)[i]-len(buffers[i])
                pending = memoryview(block[:available])
                while pending:
                    n = os.write(files[i], pending); A.require(n > 0, 'spool progress')
                    buffers[i].extend(pending[:n]); pending = pending[n:]
                A.require(len(block) <= available, 'worker output cap')
        child.wait(timeout=A.time_left(deadline))
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        failure = str(error)
    finally:
        if child is not None:
            if child.poll() is None:
                try:
                    child.kill()
                except ProcessLookupError:
                    pass
            child.wait()
            for stream in (child.stdout, child.stderr):
                try:
                    stream.close()
                except OSError as error:
                    failure = failure or 'pipe cleanup: '+str(error)
        for fd in files:
            for action in (lambda: os.fsync(fd), lambda: os.fchmod(fd, 0o400), lambda: os.close(fd)):
                try:
                    action()
                except OSError as error:
                    failure = failure or 'spool cleanup: '+str(error)
        try:
            selector.close()
        except OSError as error:
            failure = failure or 'selector cleanup: '+str(error)
    return bytes(buffers[0]), bytes(buffers[1]), None if child is None else child.returncode, failure


def run(receipt_path):
    begin = time.monotonic()
    frozen = A.read_regular(receipt_path, 1024*1024); claimed = A.decode(frozen)
    A.exact(frozen, A.canonical(claimed), 'canonical receipt'); current(claimed)
    os.mkdir(str(OUTPUT), 0o700); A.publish(OUTPUT/'CLAIM.json', frozen)
    analysis = dict(outcome='INVALID', failure=None)
    try:
        raw, error, code, failure = capture(claimed['executable'], A.sha(frozen),
                                          min(begin+180, time.monotonic()+120), [OUTPUT/'raw.jsonl', OUTPUT/'stderr.txt'])
        A.require(failure is None and code == 0 and error == b'', 'worker/observer failure: '+str(failure))
        prior = prior_header()
        analysis = verify(raw, A.sha(frozen), prior)
        current(claimed); A.require(time.monotonic()-begin < 180, 'whole controller deadline')
    except Exception as error:
        analysis = dict(outcome='INVALID', failure=str(error))
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'analysis.json', A.canonical(analysis))
    members = [pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members) < 128*1024**2-65536, 'bundle cap')
    A.publish(OUTPUT/'COMPLETE.json', A.canonical(dict(protocol=PROTOCOL, outcome=analysis['outcome'], files=members)))
    print(json.dumps({k: v for k, v in analysis.items() if k not in ('statistics',)}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('mode', choices=tuple(ARCHIVES)); b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    args = parser.parse_args()
    if args.command == 'build':
        build(args.mode, args.output)
    elif args.command == 'receipt':
        A.publish(args.output, A.canonical(receipt(args.build_dir.resolve(strict=True))))
    else:
        run(args.receipt)
