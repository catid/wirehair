#!/usr/bin/env python3
"""Frozen baseline-only alignment diagnosis. No production speed qualification."""
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
OUTPUT = Path('/var/tmp/wh2-small-basis-alignment-r0')
BASELINE_SHA = '82f9cf8985ea0873263593dd1dd945bdfcdeb993e32f473977225b6aa24090ed'
CANDIDATE_SHA = 'd10c490244d0415c6daeb035264b6a6e9484185a0476d19685e3319cece883b3'
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
HEADER = 'rep,cell,pair,order,position,side,carrier,offset,ns,handle,evaluator,basis,source,output'.split(',')
ENV = ('LD_PRELOAD','LD_LIBRARY_PATH','GLIBC_TUNABLES','MALLOC_PERTURB_',
       'MALLOC_TRIM_THRESHOLD_','MALLOC_MMAP_THRESHOLD_','MALLOC_TOP_PAD_')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def roster():
    for rep in range(12):
        for slot in range(8):
            cell = (slot+rep*3)%8
            for ps in range(10):
                pair = (ps+rep+cell)%10
                for order in range(2):
                    for pos in range(18):
                        side = SIDES[pos]^order
                        carrier = side if pair < 4 else (pair-4)//3
                        offset = pair*16 if pair < 4 else (((pair-4)%3+1)*16 if side else 0)
                        yield [rep,cell,pair,order,pos,side,carrier,offset]


def confidence(logs):
    require(len(logs) == 12 and all(math.isfinite(x) for x in logs), '12 finite replicates')
    mean = math.fsum(logs)/12
    half = 2.200985160082949*math.sqrt(math.fsum((x-mean)**2 for x in logs)/(11*12))
    return dict(ratio=math.exp(mean), lower95=math.exp(mean-half), upper95=math.exp(mean+half))


def analyze(lines):
    reader = csv.reader(lines)
    require(next(reader,None) == HEADER,'header')
    panel, groups, carriers, fixed, total = [], {}, {}, None, 0
    for expected in roster():
        raw = next(reader,None)
        require(raw is not None and len(raw) == 14 and all(x.isascii() and x.isdecimal() for x in raw),'complete numeric row')
        row = list(map(int,raw))
        require(row[:8] == expected and 0 < row[8] < 120_000_000_000,'chronology/duration')
        total += row[8]
        require(total < 120_000_000_000,'aggregate WORK cap')
        pointers = row[9:]
        require(all(0 < p < 2**64 for p in pointers),'valid pointers')
        current = tuple(row[i] for i in (9,10,12,13))
        require(all(p%64 == 0 for p in current) and len(set(current)) == 4,'aligned distinct fixed objects')
        if fixed is None:
            fixed = current
        require(current == fixed,'fixed handle/evaluator/source/output addresses')
        base = row[11]-row[7]
        require(base%4096 == 0,'page-aligned carrier origin')
        carriers.setdefault(row[6],base)
        require(carriers[row[6]] == base,'fixed carrier addresses')
        if len(carriers) == 2:
            require(carriers[1]-carriers[0] == 16384,'distinct contiguous carriers')
            ranges = [(fixed[0],512),(fixed[1],512),(fixed[2],8*1280),
                      (fixed[3]-64,64*2*8*1280+128)] + [(base-4096,16384) for base in carriers.values()]
            ranges.sort()
            require(all(0 < start < start+length < 2**64 for start,length in ranges) and
                    all(a+n <= b for (a,n),(b,_) in zip(ranges,ranges[1:])),'disjoint fixed storage')
        panel.append(row)
        if row[4] == 17:
            contrasts = []
            for pos in range(2,18,2):
                times = {r[5]:r[8] for r in panel[pos:pos+2]}
                require(set(times) == {0,1},'paired sides')
                contrasts.append(math.log(times[1])-math.log(times[0]))
            groups.setdefault(tuple(row[1:4]),[]).append(math.fsum(contrasts)/8)
            panel = []
    require(next(reader,None) is None and not panel and len(groups) == 160,'exact complete cohort')
    statistics, failed_controls, missing_effects = [], [], []
    for (cell,pair,order), logs in sorted(groups.items()):
        ci = confidence(logs)
        control = pair < 4
        primary = cell >= 6 and not control  # K8/B1280, both policies, every offset/carrier/order.
        passed = 1/1.02 < ci['lower95'] and ci['upper95'] < 1.02 if control else ci['lower95'] > 1
        if control and not passed:
            failed_controls.append([cell,pair,order])
        if primary and not passed:
            missing_effects.append([cell,pair,order])
        statistics.append(dict(cell=cell,pair=pair,order=order,primary=primary,control=control,passed=passed,**ci))
    return dict(outcome='CONTROL_FAIL' if failed_controls else 'NO_UNIFORM_PRIMARY_EFFECT' if missing_effects else 'ALIGNMENT_EFFECT_DETECTED',
                failed_controls=failed_controls, missing_primary_effects=missing_effects, statistics=statistics,
                rows=34560, fixed_addresses=fixed, carriers=carriers,
                production_speed_claimed=False, historical_slowdown_explained=False)


def trace_analyze(lines):
    reader = csv.reader(lines)
    require(next(reader,None) == 'k,width,policy,creation_order,arm,allocation,bytes,array,address,mod64,mod4096'.split(','),'trace header')
    bases = []
    for k in (3,5,8):
        for width in (64,1280):
            for policy in range(2):
                for order in range(2):
                    pair = []
                    for slot in range(2):
                        arm = slot^order
                        for allocation in range(3):
                            raw = next(reader,None)
                            require(raw is not None and len(raw) == 11 and all(x.isascii() and x.isdecimal() for x in raw),'trace row')
                            row = list(map(int,raw))
                            require(row[:6] == [k,width,policy,order,arm,allocation] and row[6]>0 and row[7] == (allocation == 1),'trace chronology')
                            require(0 < row[8] < 2**64 and row[9:] == [row[8]%64,row[8]%4096],'trace address')
                            if allocation == 1:
                                require(row[6] == k*width,'basis size')
                                pair.append(row[8]); bases.append(row)
                    require(pair[0] != pair[1],'private bases differ')
    require(next(reader,None) is None,'trace extent')
    return dict(rows=144,bases=bases,timing_claimed=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True)
        f.write('\n')


def run(build,baseline,candidate):
    build,baseline,candidate = (p.resolve(strict=True) for p in (build,baseline,candidate))
    cache = (build/'CMakeCache.txt').read_text().splitlines()
    require('SANITIZE:BOOL=OFF' in cache and 'CMAKE_HOME_DIRECTORY:INTERNAL='+str(HERE) in cache,'native diagnostic build')
    require(not any(k in os.environ for k in ENV),'unmodified allocator/loader environment')
    require(digest(baseline) == BASELINE_SHA and digest(candidate) == CANDIDATE_SHA,'exact previously recorded native DSOs')
    paths = list(HERE.glob('*')) + [baseline,candidate] + [build/p for p in ('trace','controlled','build.ninja','CMakeCache.txt')]
    pins = {str(p):digest(p) for p in paths if p.is_file()}
    OUTPUT.mkdir(mode=0o700)
    write(OUTPUT/'claim.json',dict(protocol='wh2-small-basis-alignment-r0',pins=pins,
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        cpu=50,batch=64,replicates=12,rows=34560,scope='baseline-only hot prebuilt repair alignment diagnosis; no speed promotion'))
    for name,command in (
        ('trace-normal',[str(build/'trace'),str(baseline),str(candidate),'normal']),
        ('trace-reverse',[str(build/'trace'),str(baseline),str(candidate),'reverse']),
        ('controlled',[str(build/'controlled'),str(baseline),'run'])):
        require(all(digest(p) == h for p,h in pins.items()),'input pins before launch')
        with (OUTPUT/(name+'.csv')).open('xb') as out, (OUTPUT/(name+'.stderr')).open('xb') as err:
            try:
                code = subprocess.run(command,stdout=out,stderr=err,timeout=150,check=False).returncode
            except subprocess.TimeoutExpired:
                code = 'TIMEOUT'
        require(code == 0 and (OUTPUT/(name+'.stderr')).stat().st_size == 0,'worker failed: '+str(code))
        require(all(digest(p) == h for p,h in pins.items()),'input pins after launch')
        with (OUTPUT/(name+'.csv')).open() as f:
            result = analyze(f) if name == 'controlled' else trace_analyze(f)
        result.update(exit=code,raw_sha256=digest(OUTPUT/(name+'.csv')))
        write(OUTPUT/(name+'.json'),result)
        print(name,result.get('outcome','OBSERVED'),flush=True)
    write(OUTPUT/'complete.json',{p.name:digest(p) for p in OUTPUT.iterdir() if p.is_file()})
    for p in OUTPUT.iterdir():
        p.chmod(0o400)


def replay(bundle):
    manifest = json.loads((bundle/'complete.json').read_text())
    expected = {'claim.json'} | {name+suffix for name in ('trace-normal','trace-reverse','controlled') for suffix in ('.csv','.stderr','.json')}
    require(set(manifest) == expected and {p.name for p in bundle.iterdir()} == expected | {'complete.json'},'exact bundle members')
    require(all(digest(bundle/p) == h for p,h in manifest.items()),'bundle hashes')
    for name in ('trace-normal','trace-reverse','controlled'):
        old = json.loads((bundle/(name+'.json')).read_text())
        require(old['exit'] == 0 and not (bundle/(name+'.stderr')).stat().st_size,'successful worker')
        with (bundle/(name+'.csv')).open() as f:
            fresh = analyze(f) if name == 'controlled' else trace_analyze(f)
        fresh.update(exit=0,raw_sha256=digest(bundle/(name+'.csv')))
        require(json.loads(json.dumps(fresh)) == old,'exact recorded analysis')
    return old['outcome']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run',type=Path,metavar='BUILD')
    group.add_argument('--analyze',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--candidate',type=Path)
    args = parser.parse_args()
    if args.run:
        require(args.baseline is not None and args.candidate is not None,'DSO paths required')
        run(args.run,args.baseline,args.candidate)
    else:
        print(replay(args.analyze))
