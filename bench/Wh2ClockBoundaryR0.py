#!/usr/bin/env python3
"""Build/qualify/run a one-shot codec-free clock/computation diagnostic."""
import argparse
import array
import functools
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = Path('/var/tmp/wh2-clock-boundary-r0')
PROTOCOL = 'wirehair.wh2.clock-boundary-r0'
SAMPLES, NEUTRAL, ITERATIONS = 524288, 256, 65536
SEED = 0x9e3779b97f4a7c15
RAW_CAP, ERR_CAP = 192 * 1024**2, 65536
SOURCES = ('bench/Wh2ClockBoundaryR0.cpp', 'bench/Wh2ClockBoundaryR0.py',
           'bench/test_Wh2ClockBoundaryR0.py', 'bench/Wh2ClockBoundaryR0.md',
           'bench/Wh2FrozenTrace.h', 'bench/Wh2FrozenTrace.cpp')
ENV_KEYS = ('LD_PRELOAD', 'LD_LIBRARY_PATH', 'GLIBC_TUNABLES', 'MALLOC_PERTURB_',
            'MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def pin(path):
    path = Path(path).resolve(strict=True)
    h, size = hashlib.sha256(), 0
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024**2), b''):
            h.update(data)
            size += len(data)
    return dict(path=str(path), bytes=size, sha256=h.hexdigest())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, sort_keys=True, separators=(',', ':'), allow_nan=False)
        stream.write('\n')


def reject_constant(value):
    raise ValueError('nonfinite JSON: ' + value)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def decode(line):
    return json.loads(line, object_pairs_hook=unique_object, parse_constant=reject_constant)


@functools.lru_cache(None)
def answer():
    value = SEED
    for _ in range(ITERATIONS):
        value ^= (value << 13) & ((1 << 64) - 1)
        value ^= value >> 7
        value ^= (value << 17) & ((1 << 64) - 1)
    return value


def identity(cpuid):
    require(isinstance(cpuid, list) and len(cpuid) == 7, 'CPUID shape')
    ids = [0, 1, 0x80000000, 0x80000001, 0x80000007, 0x80000008, 0x80000021]
    require(all(isinstance(r, list) and len(r) == 5 and all(type(x) is int and 0 <= x < 2**32 for x in r) for r in cpuid), 'CPUID words')
    require([r[0] for r in cpuid] == ids, 'CPUID leaves')
    require([cpuid[0][j] for j in (2, 4, 3)] == [0x68747541, 0x69746e65, 0x444d4163], 'AMD vendor')
    require(cpuid[1][1] == 0x00b00f81 and cpuid[1][2] >> 24 == 100 and not cpuid[1][3] & (1 << 31), 'physical CPU50')
    require(cpuid[3][4] & (1 << 27) and cpuid[4][4] & (1 << 8) and cpuid[6][1] & 4, 'ordered invariant RDTSCP')


def records(path):
    with Path(path).open('rb') as stream:
        while True:
            line = stream.readline(8193)
            if not line:
                return
            require(len(line) <= 8192 and line.endswith(b'\n'), 'bounded complete line')
            yield decode(line)


def analyze(path, count, claim, neutral=False, failure=None):
    raw_pin = pin(path)
    stream = iter(records(path))
    h = next(stream)
    require(set(h) == {'type','protocol','claim','neutral','samples','iterations','seed','answer','aux','start','cpuid'}, 'header schema')
    require(h['type'] == 'header' and h['protocol'] == PROTOCOL and h['claim'] == claim, 'header identity')
    require(h['neutral'] is neutral and h['samples'] == count and h['iterations'] == ITERATIONS and
            h['seed'] == SEED and h['answer'] == answer() and h['aux'] == 50, 'fixed workload')
    require(len(h['start']) == 2 and all(type(x) is int and x >= 0 for x in h['start']), 'start clocks')
    identity(h['cpuid'])
    durations = [array.array('Q') for _ in range(3)]
    phase_durations = [array.array('Q') for _ in range(48)]
    counter_totals, between_totals = [0]*4, [0]*4
    first, previous = None, None
    for index in range(count):
        r = next(stream)
        require(isinstance(r, list) and len(r) == 23 and all(type(x) is int and 0 <= x < 2**64 for x in r), 'record words')
        require(r[0] == index, 'complete ordinal sequence')
        bad = failure if index + 1 == count else None
        if bad == 'clock':
            require(r[22] == 7 and r[4] == r[8] == r[2] == r[21] == 0 and r[9] == answer(), 'retained partial clock failure')
            require(r[5] < r[6] < r[7] and r[18:21] == [50]*3, 'partial ordered stamps')
            require(all(r[11+2*j] == 0 for j in range(4)), 'partial counters absent')
            require(previous is not None and previous[2] <= r[1] and previous[4] <= r[3] and previous[8] < r[5], 'partial inter-record clocks')
            require(all(previous[11+2*j] <= r[10+2*j] for j in range(4)), 'partial inter-record counters')
            continue
        require(r[22] == 11 and r[9] == (answer() ^ (bad == 'result')), 'computation result and stage')
        require(r[1] <= r[2] and r[3] < r[4] and r[5] < r[6] < r[7] < r[8], 'clock/TSC ordering')
        require(r[18:22] == [h['aux']]*4, 'AUX migration')
        for j in range(4):
            require(r[10+2*j] <= r[11+2*j], 'counter ordering')
            counter_totals[j] += r[11+2*j] - r[10+2*j]
        if previous:
            require(previous[2] <= r[1] and previous[4] <= r[3] and previous[8] < r[5], 'inter-record clocks')
            for j in range(4):
                require(previous[11+2*j] <= r[10+2*j], 'inter-record counters')
                between_totals[j] += r[10+2*j] - previous[11+2*j]
        else:
            require(h['start'][0] <= r[3] and h['start'][1] <= r[1], 'startup clocks')
            first = r
        for j in range(3):
            durations[j].append(r[6+j] - r[5+j])
        phase_durations[(r[3] % 1000000)*48//1000000].append(r[7]-r[6])
        previous = r
    footer = next(stream)
    require(set(footer) == {'type','complete','records','codec_calls','end'}, 'footer schema')
    require(footer['type'] == 'footer' and footer['complete'] is (failure is None) and
            footer['records'] == count and footer['codec_calls'] == 0, 'footer completion')
    require(len(footer['end']) == 2 and all(type(x) is int and x >= 0 for x in footer['end']) and
            footer['end'][0] >= previous[4] and footer['end'][1] >= previous[2], 'footer clocks')
    require(next(stream, None) is None, 'trailing records')
    if neutral:
        require(pin(path) == raw_pin, 'raw changed during neutral analysis')
        return dict(outcome='EXPECTED_FAILURE' if failure else 'NEUTRAL_PASS', records=count,
                    cpuid=h['cpuid'], answer=answer(), codec_calls=0)
    require(failure is None and count == SAMPLES, 'scientific fixed count')
    require(all(phase_durations), 'all 48 passive absolute phases represented')
    scale = (previous[8]-first[5]) / (previous[4]-first[3])
    medians = [statistics.median(d) for d in durations]
    long_counts = [sum(t >= 50000*scale for t in d) for d in durations]
    events = []
    for r in records(path):
        if not isinstance(r, list):
            continue
        ticks = [r[6+j]-r[5+j] for j in range(3)]
        flags = [ticks[0] >= 50000*scale, ticks[1] >= 1.5*medians[1], ticks[2] >= 50000*scale]
        if any(flags):
            events.append(dict(record=r, ticks=ticks, ns_equivalent=[t/scale for t in ticks], flags=flags,
                               counter_delta=[r[11+2*j]-r[10+2*j] for j in range(4)]))
    require(pin(path) == raw_pin, 'raw changed during analysis')
    return dict(outcome='DIAGNOSTIC_COMPLETE', speed_qualified=False, codec_calls=0, records=count,
                iterations=count*ITERATIONS, tick_per_ns=scale, median_ticks=medians,
                maxima_ticks=[max(d) for d in durations], counts_ge_50us=long_counts,
                counter_totals=counter_totals, between_counter_totals=between_totals,
                phase_counts=[len(d) for d in phase_durations],
                phase_compute_median_ticks=[statistics.median(d) for d in phase_durations],
                phase_compute_max_ticks=[max(d) for d in phase_durations], events=events)


def command(args, **kwargs):
    return subprocess.run(list(map(str,args)), cwd=ROOT, check=True, timeout=60,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs).stdout


def snapshot_env():
    return {k: os.environ.get(k) for k in ENV_KEYS}


def check_pins(manifest):
    require(snapshot_env() == manifest['environment'], 'environment changed')
    for expected in manifest['pins']:
        require(pin(expected['path']) == expected, 'pin changed: '+expected['path'])


def build(directory):
    directory = Path(directory).resolve()
    directory.mkdir(parents=False, exist_ok=False)
    compiler = Path(shutil.which('g++')).resolve()
    files = {ROOT/name for name in SOURCES}
    files.update(Path(command([compiler, '-print-prog-name='+name]).decode().strip()).resolve() if name == 'cc1plus'
                 else Path(shutil.which(name)).resolve() for name in ('cc1plus','as','ld'))
    files.add(compiler)
    commands, neutral_results = [], []
    for interpreter in ('python3','python3.8'):
        executable = Path(command([interpreter,'-c','import sys; print(sys.executable)']).decode().strip()).resolve()
        files.add(executable)
        result = subprocess.run([str(executable),'-m','unittest','bench.test_Wh2ClockBoundaryR0'],
                                cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
        require(result.returncode == 0, 'reader unit tests')
        write(directory/(interpreter+'-tests.json'),dict(command=[str(executable),'-m','unittest','bench.test_Wh2ClockBoundaryR0'],
                                                       stdout=result.stdout.decode(),stderr=result.stderr.decode(),exit=result.returncode))
        files.add(directory/(interpreter+'-tests.json'))
    for mode in ('native','asan'):
        out = directory/mode
        out.mkdir()
        flags = ['-std=c++11','-Wall','-Wextra','-Wconversion','-Wsign-conversion','-Werror',
                 '-fno-lto','-fno-pie','-ffunction-sections','-fdata-sections']
        flags += ['-O3','-DWH2_CLOCK_NEUTRAL_ONLY=0'] if mode == 'native' else [
            '-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer','-DWH2_CLOCK_NEUTRAL_ONLY=1']
        objects = []
        for name in ('Wh2ClockBoundaryR0','Wh2FrozenTrace'):
            obj, dep = out/(name+'.o'), out/(name+'.d')
            args = [compiler]+flags+['-MD','-MF',dep,'-c',ROOT/'bench'/(name+'.cpp'),'-o',obj]
            command(args); commands.append(list(map(str,args)))
            deps = shlex.split(dep.read_text().replace('\\\n',' ').split(':',1)[1])
            files.update((ROOT/p).resolve() for p in deps)
            files.update((obj,dep)); objects.append(obj)
        worker = out/'worker'
        args = [compiler]+flags+['-no-pie']+objects+['-Wl,--gc-sections','-o',worker]
        command(args); commands.append(list(map(str,args))); files.add(worker)
        ldd = command(['ldd',worker]).decode()
        files.update(Path(p).resolve() for p in re.findall(r'(?:=>\s+|^\s*)(/\S+)',ldd,re.M))
        contract = decode(command([worker,'--contract']))
        require(contract == dict(samples=SAMPLES,neutral_samples=NEUTRAL,iterations=ITERATIONS,
                                 cpu_seconds=100,wall_seconds=120,address_space_mib=128,
                                 claim_path=str(OUTPUT/'CLAIM.json'),neutral_only=mode == 'asan'), 'native frozen contract')
        write(out/'runtime.json',dict(ldd=ldd,contract=contract))
        files.add(out/'runtime.json')
        env = dict(os.environ, ASAN_OPTIONS='detect_leaks=1:detect_stack_use_after_return=1:abort_on_error=1',
                   UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
        write(out/'neutral-claim.json',dict(neutral=True))
        fixture = pin(out/'neutral-claim.json')
        auth_results = []
        for path,sha,expected in [(fixture['path'],fixture['sha256'],0),
                                  (fixture['path'],'0'*64,1),(str(out/'missing-claim'),'0'*64,1)]:
            r = subprocess.run([str(worker),'--neutral-claim',path,sha],env=env,cwd=ROOT,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20)
            require(r.returncode == expected and bool(r.stderr) == bool(expected), 'claim authentication check')
            auth_results.append(dict(path=path,sha256=sha,exit=r.returncode,stdout=r.stdout.decode(),stderr=r.stderr.decode()))
        write(out/'claim-tests.json',auth_results)
        files.update((out/'neutral-claim.json',out/'claim-tests.json'))
        if mode == 'asan':
            r = subprocess.run([str(worker),'--worker','0'*64],env=env,cwd=ROOT,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20)
            require(r.returncode == 1 and b'neutral worker cannot measure' in r.stderr and not r.stdout, 'ASan scientific mode rejection')
            write(out/'mode-rejection.json',dict(exit=r.returncode,stdout=r.stdout.decode(),stderr=r.stderr.decode()))
            files.add(out/'mode-rejection.json')
        for suffix,failure in (('neutral',None),('neutral-bad-result','result'),('neutral-fail-clock','clock')):
            raw, err = out/(suffix+'.jsonl'), out/(suffix+'.stderr')
            with raw.open('xb') as stdout, err.open('xb') as stderr:
                r = subprocess.run([str(worker),'--'+suffix],stdout=stdout,stderr=stderr,env=env,timeout=20,cwd=ROOT)
            require(r.returncode == (1 if failure else 0), 'neutral exit')
            require(bool(err.stat().st_size) == bool(failure), 'neutral stderr')
            neutral_results.append(analyze(raw,NEUTRAL,'0'*64,True,failure))
            files.update((raw,err))
    require(all(r['cpuid'] == neutral_results[0]['cpuid'] for r in neutral_results), 'backend identity')
    write(directory/'build.json', dict(commands=commands,neutral=neutral_results))
    files.add(directory/'build.json')
    files.update(Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None) and Path(m.__file__).is_file())
    manifest = dict(protocol=PROTOCOL,environment=snapshot_env(),pins=[pin(p) for p in sorted(files)],cpuid=neutral_results[0]['cpuid'],
                    native=str(directory/'native/worker'),neutral_pass=True)
    write(directory/'manifest.json',manifest)
    print(json.dumps(dict(build=str(directory),pins=len(files),neutral_runs=len(neutral_results))))


def run(manifest_path):
    manifest_path = Path(manifest_path).resolve()
    manifest = decode(manifest_path.read_bytes())
    require(manifest['protocol'] == PROTOCOL and manifest['neutral_pass'] is True, 'qualified build')
    check_pins(manifest)
    for source in SOURCES:
        require(command(['git','show','HEAD:'+source]) == (ROOT/source).read_bytes(), 'uncommitted producing source: '+source)
    head = command(['git','rev-parse','HEAD']).decode().strip()
    claim = dict(protocol=PROTOCOL,head=head,manifest=pin(manifest_path),pins=manifest['pins'],environment=snapshot_env(),
                 interpreter=pin(sys.executable),native=manifest['native'])
    OUTPUT.mkdir(exist_ok=False)
    write(OUTPUT/'CLAIM.json',claim)
    claim_sha = pin(OUTPUT/'CLAIM.json')['sha256']
    start = time.monotonic()
    failure, process = None, None
    try:
        with (OUTPUT/'raw.jsonl').open('xb') as stdout, (OUTPUT/'stderr.txt').open('xb') as stderr:
            process = subprocess.Popen([manifest['native'],'--worker',claim_sha],stdout=stdout,stderr=stderr,cwd=ROOT)
            while process.poll() is None:
                require(time.monotonic()-start < 150, 'observer wall cap')
                require((OUTPUT/'raw.jsonl').stat().st_size <= RAW_CAP and (OUTPUT/'stderr.txt').stat().st_size <= ERR_CAP, 'output cap')
                time.sleep(.05)
            require((OUTPUT/'raw.jsonl').stat().st_size <= RAW_CAP and (OUTPUT/'stderr.txt').stat().st_size <= ERR_CAP, 'terminal output cap')
            require(process.returncode == 0 and (OUTPUT/'stderr.txt').stat().st_size == 0, 'worker success')
        check_pins(manifest)
        require(pin(manifest_path) == claim['manifest'], 'manifest changed')
        require(pin(OUTPUT/'CLAIM.json')['sha256'] == claim_sha, 'claim changed')
        require(pin(sys.executable) == claim['interpreter'], 'interpreter changed')
        require(command(['git','rev-parse','HEAD']).decode().strip() == head, 'source HEAD changed')
        require(next(records(OUTPUT/'raw.jsonl'))['cpuid'] == manifest['cpuid'], 'qualified target changed')
        analysis = analyze(OUTPUT/'raw.jsonl',SAMPLES,claim_sha)
    except Exception as error:
        failure = str(error)
        if process is not None and process.poll() is None:
            process.kill(); process.wait(timeout=10)
        analysis = dict(outcome='INVALID',error=failure,speed_qualified=False)
    analysis['elapsed_seconds'] = time.monotonic()-start
    write(OUTPUT/'analysis.json',analysis)
    members = [p for p in OUTPUT.iterdir() if p.is_file()]
    write(OUTPUT/'COMPLETE.json',dict(protocol=PROTOCOL,outcome=analysis['outcome'],files=[pin(p) for p in sorted(members)]))
    for p in OUTPUT.iterdir():
        p.chmod(0o400)
    print(json.dumps(dict(outcome=analysis['outcome'],directory=str(OUTPUT))))
    require(failure is None,'diagnostic invalid: '+str(failure))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('build','run','analyze'))
    p.add_argument('path')
    args = p.parse_args()
    if args.mode == 'build':
        build(args.path)
    elif args.mode == 'run':
        run(args.path)
    else:
        print(json.dumps(analyze(args.path,SAMPLES,pin(OUTPUT/'CLAIM.json')['sha256']),sort_keys=True))


if __name__ == '__main__':
    main()
