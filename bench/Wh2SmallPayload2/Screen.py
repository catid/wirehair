#!/usr/bin/env python3
"""Prospective bounded diagnostic. Never a production/all-K promotion receipt."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUTPUT = Path('/var/tmp/wh2-small-payload2-screen-r0')
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
PAIRS = ((0,0),(1,1),(2,2),(0,1),(2,1))
HEADER = ['rep','cell','pair','order','position','arm','ns']
ENV = ('LD_PRELOAD','LD_LIBRARY_PATH','GLIBC_TUNABLES','MALLOC_PERTURB_',
       'MALLOC_TRIM_THRESHOLD_','MALLOC_MMAP_THRESHOLD_','MALLOC_TOP_PAD_')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def roster():
    for rep in range(12):
        for slot in range(96):
            cell = (slot + rep*17) % 96
            for ps in range(5):
                pair = (ps+rep+cell) % 5
                for order in range(2):
                    for pos in range(18):
                        yield [rep, cell, pair, order, pos, PAIRS[pair][SIDES[pos]^order]]


def shape(cell):
    metric, fixture = cell % 4, cell // 4
    policy, index = fixture % 2, fixture // 2
    k, variant = (3,5,8)[index//4], index%4
    width = (2,2,64,1280)[variant]
    return dict(k=k, width=width, tail=1 if variant == 1 else width,
                policy=('independent','borrowed')[policy],
                metric=('prebuilt-repair','full-encoder','low-decoder','distant-decoder')[metric])


def confidence(logs):
    require(len(logs) == 12 and all(math.isfinite(x) for x in logs), '12 finite replicates')
    mean = math.fsum(logs)/12
    se = math.sqrt(math.fsum((x-mean)**2 for x in logs)/(11*12))
    half = 2.200985160082949*se
    return dict(ratio=math.exp(mean), lower95=math.exp(mean-half), upper95=math.exp(mean+half))


def analyze(stream):
    reader = csv.reader(stream)
    require(next(reader, None) == HEADER, 'header')
    groups, panel, total_ns = {}, [], 0
    for expected in roster():
        raw = next(reader, None)
        require(raw is not None and len(raw) == 7 and all(x.isascii() and x.isdecimal() for x in raw), 'complete numeric row')
        row = list(map(int, raw))
        require(row[:6] == expected and 0 < row[6] < 180_000_000_000, 'chronology/positive duration')
        total_ns += row[6]
        require(total_ns < 180_000_000_000, 'aggregate WORK cap')
        panel.append(row)
        if expected[4] == 17:
            contrasts = []
            for pos in range(2,18,2):
                times = {SIDES[r[4]]^r[3]:r[6] for r in panel[pos:pos+2]}
                require(set(times) == {0,1}, 'paired sides')
                contrasts.append(math.log(times[1])-math.log(times[0]))
            groups.setdefault(tuple(expected[1:4]), []).append(math.fsum(contrasts)/8)
            panel = []
    require(next(reader, None) is None and not panel and len(groups) == 960, 'exact complete cohort')
    stats, bad_controls, bad_candidate, wh1_not_faster = [], [], [], []
    for (cell,pair,order), logs in sorted(groups.items()):
        s, ci = shape(cell), confidence(logs)
        item = dict(cell=cell, pair=pair, order=order, shape=s, **ci)
        if pair < 3:
            passed = 1/1.02 < ci['lower95'] and ci['upper95'] < 1.02
            if not passed:
                bad_controls.append([cell,pair,order])
        elif pair == 3:
            # Strict B2 encoder/repair win. All wider or decoder paths must
            # pass the same prospectively declared 2% screening envelope.
            limit = 1.0 if s['width'] == 2 and cell%4 < 2 else 1.02
            passed = ci['upper95'] < limit
            if not passed:
                bad_candidate.append([cell,pair,order])
        else:
            passed = ci['upper95'] < 1.0
            if not passed:
                wh1_not_faster.append([cell,pair,order])
        item['pass'] = passed
        stats.append(item)
    return dict(outcome='CONTROL_FAIL' if bad_controls else 'FAIL' if bad_candidate else 'PASS',
                failed_controls=bad_controls, failed_candidate=bad_candidate,
                candidate_not_proven_faster_than_WH1=wh1_not_faster, statistics=stats,
                rows=207360, production_promotion_claimed=False, recovery_rate_claimed=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write('\n')


def run(build, output=None, here=None, protocol='wh2-small-payload2-screen-r0'):
    # Additional experiments must supply their own fixed namespace and source
    # directory. The original default namespace remains spent, never reusable.
    output = OUTPUT if output is None else output
    here = HERE if here is None else here
    build = build.resolve(strict=True)
    require(not any(os.environ.get(k) is not None for k in ENV), 'unmodified allocator/loader environment')
    cache = (build/'CMakeCache.txt').read_text()
    require(all(f'{name}:BOOL=OFF' in cache for name in ('SANITIZE','PORTABLE','COUNT')), 'native uncounted build')
    require((build/'libbaseline.so').resolve() != (build/'libcandidate.so').resolve() and
            digest(build/'libbaseline.so') != digest(build/'libcandidate.so'), 'distinct baseline/candidate files')
    # Exact content pins for this diagnostic; not the older formal campaign's
    # transitive toolchain/loader provenance contract.
    sources = subprocess.check_output(['git','ls-files'], cwd=ROOT, text=True).splitlines()
    paths = [ROOT/p for p in sources if (ROOT/p).suffix in ('.cpp','.c','.h','.inc','.cmake') or p == 'CMakeLists.txt']
    paths += list(HERE.glob('*')) + list(here.glob('*'))
    paths += [build/p for p in ('screen','libbaseline.so','libcandidate.so','build.ninja','compile_commands.json','CMakeCache.txt')]
    paths += list((build/'candidate').glob('*'))
    pins = {str(p):digest(p) for p in paths if p.is_file()}
    output.mkdir(mode=0o700)  # A spent namespace is never overwritten/retried.
    write(output/'claim.json', dict(protocol=protocol, pins=pins,
        head=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        batch=32, replicates=12, rows_per_order=207360, cpu=50, limits='2% AA/retention; strict B2 encoder/repair gain',
        scope='diagnostic only; not all-K, installed-package, or production promotion'))
    results = {}
    for mode in ('run','run-reverse'):
        require(all(digest(p) == h for p,h in pins.items()), 'inputs changed before launch')
        with (output/(mode+'.csv')).open('xb') as raw, (output/(mode+'.stderr')).open('xb') as error:
            try:
                completed = subprocess.run([str(build/'screen'), str(build/'libbaseline.so'),
                    str(build/'libcandidate.so'), mode], stdout=raw, stderr=error, timeout=240, check=False)
                code = completed.returncode
            except subprocess.TimeoutExpired:
                code = 'TIMEOUT'
        try:
            require(code == 0 and (output/(mode+'.stderr')).stat().st_size == 0, 'worker failed: '+str(code))
            require(all(digest(p) == h for p,h in pins.items()), 'inputs changed during launch')
            with (output/(mode+'.csv')).open() as f:
                result = analyze(f)
        except (ValueError, OSError) as error:
            result = dict(outcome='INVALID', failure=str(error))
        result.update(exit=code, raw_sha256=digest(output/(mode+'.csv')))
        results[mode] = result
        write(output/(mode+'.json'), result)
        print(mode, result['outcome'], 'controls', len(result.get('failed_controls',[])),
              'candidate', len(result.get('failed_candidate',[])), flush=True)
    write(output/'complete.json', {p.name:digest(p) for p in output.iterdir() if p.is_file()})
    for p in output.iterdir():
        p.chmod(0o400)


def replay(bundle):
    members = json.loads((bundle/'complete.json').read_text())
    expected = {'claim.json','run.csv','run.stderr','run.json',
                'run-reverse.csv','run-reverse.stderr','run-reverse.json'}
    require(set(members) == expected, 'bundle members')
    require({p.name for p in bundle.iterdir()} == expected | {'complete.json'}, 'exact bundle')
    require(all(digest(bundle/name) == value for name,value in members.items()), 'bundle hashes')
    result = {}
    for mode in ('run','run-reverse'):
        recorded = json.loads((bundle/(mode+'.json')).read_text())
        require(recorded['exit'] == 0 and (bundle/(mode+'.stderr')).stat().st_size == 0,
                'original worker must have completed successfully')
        with (bundle/(mode+'.csv')).open() as f:
            fresh = analyze(f)
        fresh.update(exit=0, raw_sha256=digest(bundle/(mode+'.csv')))
        require(fresh == recorded, 'exact recorded analysis')
        result[mode] = fresh['outcome']
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', type=Path)
    group.add_argument('--analyze', type=Path)
    args = parser.parse_args()
    if args.run:
        run(args.run)
    else:
        print(json.dumps(replay(args.analyze), sort_keys=True))
