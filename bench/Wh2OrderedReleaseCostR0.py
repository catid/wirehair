#!/usr/bin/env python3
"""Fixed .78 full-cost gate. Build/neutral work never consumes its namespace."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import selectors
import shlex
import subprocess
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = sibling('ordered_cost_common', 'Wh2AlignedIntermediateCostR0.py')
I = sibling('ordered_cost_isolation', 'Wh2OrderedRhsReleaseIsolation.py')
ROOT = I.ROOT
PROTOCOL = 'wirehair.wh2.ordered-release-cost-r0'
OUTPUT = Path('/var/tmp/wh2-ordered-release-cost-r0')
ISOLATED = Path('/tmp/wh2-release-isolation.KAG4qU')
MANIFESTS = dict(native='f2c634eadfc55ca9b5bebbdc40c47e290ef35f0b84e02df18475726fcbc68e49',
                 asan='b254b72e61cebc4aa7d429eee06f5c56916eaaa29227ddfc09d310a5b92a9b5f',
                 scalar='26e89927b1aab43a0d041ddec101efb0f287dff52bf3c74c35a2b193ec3b267e')
PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (2, 1))
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0)
WIDTHS = (2, 64, 1280)
CALLBACKS, BATCH = 19440, 128
RAW_CAP, ERR_CAP = 64*1024**2, 65536
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
NEW = ('bench/Wh2OrderedReleaseCostR0.cpp', 'bench/Wh2OrderedReleaseCostR0.py',
       'bench/test_Wh2OrderedReleaseCostR0.py')


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
    results, failed_controls, failed_local, wh1_not_faster = [], [], [], []
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
        elif comparison == 3:
            limit = 0.0 if width == 2 and metric > 0 else bound
            item['local_pass'] = estimate['upper95_log'] < limit
            item['upper_ratio_limit'] = math.exp(limit)
            if not item['local_pass']:
                failed_local.append(list(key))
        else:
            if estimate['upper95_log'] >= 0:
                wh1_not_faster.append(list(key))
        results.append(item)
    return dict(outcome='CONTROL_FAIL' if failed_controls else 'FAIL' if failed_local else 'PASS',
                statistics=results, failed_controls=failed_controls, failed_local=failed_local,
                WH1_not_faster=wh1_not_faster, WH1_qualified=not failed_controls,
                all_WH1_faster=not failed_controls and not wh1_not_faster,
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


def verify_fixtures(header, old):
    A.exact(len(header['fixtures']), 3, 'fixture count')
    for f, prior, width in zip(header['fixtures'], old['fixtures'], WIDTHS):
        A.exact(f['width'], width, 'fixture width')
        A.exact(f['source'], bytes((37*i+i//11) % 256 for i in range(6*width)).hex(), 'fixed source')
        A.exact(len(f['arms']), 3, 'three fixtures')
        public_profile = next(h['profile_hex'] for h in prior['handles'] if h['arm'] == 3)
        for a, arm in enumerate(f['arms']):
            A.exact(arm['profile'], public_profile if a < 2 else '00'*32, 'actual public descriptor')
            A.exact(arm['packets'], prior['packets_hex'][3 if a < 2 else 2], 'unchanged packet corpus')
            A.exact(len(arm['steps']), 2, 'decoder families')
            for count in arm['steps']:
                A.integer(count, 6, 12)
        A.exact(f['arms'][0], f['arms'][1], 'B/C fixture parity')


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


def build(mode, output):
    output = output.parent.resolve(strict=True)/output.name
    A.require(output.is_absolute() and ROOT not in output.parents and output != ROOT and
              not output.exists() and not output.is_symlink(), 'fresh external build')
    isolation = ISOLATED/(mode+'-final')
    raw = A.read_regular(isolation/'manifest.json', 1024*1024)
    A.exact(A.sha(raw), MANIFESTS[mode], 'qualified isolation manifest')
    selected = A.decode(raw)
    for p in selected['inputs']+selected['artifacts']:
        A.exact(pin(Path(p['path'])), p, 'unchanged isolated artifacts')
    output.mkdir(mode=0o700)
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_RELEASE_COST_NEUTRAL='+str(int(mode != 'native')),
             '-I'+str(ROOT), '-I'+str(ROOT/'bench'), '-I'+str(ROOT/'include'), '-I'+str(isolation)]
    flags += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer'] if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags += ['-DANDROID=1']
    sources = [ROOT/NEW[0]]+[ROOT/'bench'/name for name in
               ('Wh2FrozenTrace.cpp', 'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp')]
    commands, objects, dependencies = [], [], {ROOT/n for n in NEW}
    dependencies.update(Path(p['path']) for p in selected['inputs']+selected['artifacts'])
    dependencies.update((isolation/'manifest.json', ROOT/'bench/Wh2AlignedIntermediateCostR0.py'))
    for source in sources:
        obj = output/(source.stem+'.o'); dep = output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)]
        command(args); commands.append(args); objects.append(obj)
        dependencies.update(Path(s).resolve(strict=True) for s in shlex.split(dep.read_text().replace('\\\n', '').split(': ', 1)[1]))
    archive = I.QUALIFIED/mode/'production/libwirehair.a'
    objects.extend(isolation/name for name in I.MEMBERS)
    executable = output/'cost_worker'
    args = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan':
        args += ['-fsanitize=address,undefined']
    args += list(map(str, objects))+[str(archive), '-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    command(args); commands.append(args)
    symbols = command(['/usr/bin/nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    A.require(names.count('GF256Ctx') == names.count('gf256_init_') == 1, 'single GF runtime')
    A.require(all(names.count(prefix+name) == 1 for prefix in ('', I.PREFIX) for name in
                  ('wirehair_v2_encode', 'wirehair_v2_decode', 'wirehair_v2_free')), 'both actual V2 APIs')
    A.require(all(names.count(n) == 1 for n in ('wirehair_encode', 'wirehair_decode', 'wirehair_free')), 'actual WH1 API')
    text = command(['/usr/bin/nm', '-C', executable]).decode().splitlines()
    A.require(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line for line in text) == 1, 'one common WORK body')
    for name in ('c++', 'as', 'ld', 'nm', 'objcopy'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
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
    old = A.decode(A.read_regular(I.QUALIFIED/'receipt.json', 1024*1024))
    pins = {p['path']: p for p in old['pins']+manifest['inputs']+manifest['artifacts']}
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
        with Path('/var/tmp/wh2-thue-decoder-cost-r0/raw.jsonl').open('rb') as stream:
            prior = A.decode(stream.readline(1024*1024))
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
    b = sub.add_parser('build'); b.add_argument('mode', choices=tuple(MANIFESTS)); b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    args = parser.parse_args()
    if args.command == 'build':
        build(args.mode, args.output)
    elif args.command == 'receipt':
        A.publish(args.output, A.canonical(receipt(args.build_dir.resolve(strict=True))))
    else:
        run(args.receipt)
