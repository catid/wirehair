#!/usr/bin/env python3
"""One-shot alignment candidate screen; not a production/all-K qualification."""
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
OUTPUT = Path('/var/tmp/wh2-aligned-small-basis-screen-r0')
LIBRARY_SHA = ('bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe',
               '59a98ab34920cf40c2cfba6aa26b63d8ff929e9ec0823050ab07d21bb2ee05c7')
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
PAIRS = ((0,0),(1,1),(2,2),(0,1),(2,1))
NATURAL = 'rep,cell,pair,order,position,arm,ns,source,output,profile,low,distant'.split(',')
CONTROLLED = 'rep,cell,pair,order,position,side,arm,carrier,offset,ns,handle,evaluator,raw_basis,basis,source,output'.split(',')
TRACE = 'k,width,policy,creation_order,arm,allocation,bytes,array,address,mod64,mod4096,view,view_mod64,view_mod4096'.split(',')
ENV = ('LD_PRELOAD','LD_LIBRARY_PATH','GLIBC_TUNABLES','MALLOC_PERTURB_',
       'MALLOC_TRIM_THRESHOLD_','MALLOC_MMAP_THRESHOLD_','MALLOC_TOP_PAD_')
WORKERS = ('trace-normal','trace-reverse','controlled','controlled-reverse','natural','natural-reverse')


def require(ok, why):
    if not ok:
        raise ValueError(why)


def shape(cell):
    fixture = cell//3
    if fixture < 42:
        k, variant = (3,5,8)[fixture//14], (fixture//2)%7
        width = (2,2,64,256,257,1280,1280)[variant]
        tail = 1 if variant in (1,6) else width
    else:
        k, width = 16, (64,1280)[(fixture-42)//2]
        tail = width
    return dict(k=k,width=width,tail=tail,policy=fixture%2,metric=cell%3,
                route='ordinary' if k == 3 else 'certified' if k == 16 else 'explicit-small')


def roster(controlled):
    cells, pairs, shift = (6,16,5) if controlled else (138,5,17)
    for rep in range(12):
        for slot in range(cells):
            cell = (slot+rep*shift)%cells
            for ps in range(pairs):
                pair = (ps+rep+cell)%pairs
                for order in range(2):
                    for pos in range(18):
                        side = SIDES[pos]^order
                        if controlled:
                            arm = pair//4 if pair < 8 else side
                            carrier = side if pair < 8 else (pair-8)//4
                            offset = pair%4*16
                            yield [rep,cell,pair,order,pos,side,arm,carrier,offset]
                        else:
                            yield [rep,cell,pair,order,pos,PAIRS[pair][side]]


def numeric(raw, length):
    require(raw is not None and len(raw) == length and all(x.isascii() and x.isdecimal() for x in raw),'complete numeric row')
    return list(map(int,raw))


def disjoint(ranges):
    ranges = sorted(ranges)
    require(all(0 < a < a+n < 2**64 for a,n in ranges) and
            all(a+n <= b for (a,n),(b,_) in zip(ranges,ranges[1:])),'valid disjoint storage')


def confidence(logs):
    require(len(logs) == 12 and all(math.isfinite(x) for x in logs),'12 finite replicates')
    mean = math.fsum(logs)/12
    half = 2.200985160082949*math.sqrt(math.fsum((x-mean)**2 for x in logs)/(11*12))
    return dict(ratio=math.exp(mean),lower95=math.exp(mean-half),upper95=math.exp(mean+half))


def analyze(lines, controlled):
    reader = csv.reader(lines)
    header = CONTROLLED if controlled else NATURAL
    require(next(reader,None) == header,'header')
    panel, groups, fixed, carriers, total, count = [], {}, {}, {}, 0, 0
    duration = 9 if controlled else 6
    for expected in roster(controlled):
        row = numeric(next(reader,None),len(header))
        require(row[:duration] == expected and 0 < row[duration] < 120_000_000_000,'chronology/duration')
        count += 1
        total += row[duration]
        require(total < 120_000_000_000,'aggregate WORK cap including warmups')
        require(all(0 < p < 2**64 for p in row[duration+1:]),'valid pointers')
        if controlled:
            current = tuple(row[i] for i in (10,11,14,15))
            require(all(p%64 == 0 for p in current),'fixed object alignment')
            fixed.setdefault('process',current)
            require(fixed['process'] == current,'fixed process addresses')
            base = row[12]-row[8]
            require(base%4096 == 0 and row[13] == row[12]+(64-row[8] if row[6] and row[8] else 0),'raw/effective view geometry')
            carriers.setdefault(row[7],base)
            require(carriers[row[7]] == base,'fixed carriers')
            if len(carriers) == 2:
                require(carriers[1]-carriers[0] == 16384,'carrier separation')
                disjoint([(current[0],512),(current[1],512),(current[2],8*1280),
                          (current[3]-64,64*2*8*1280+128)]+[(b-4096,16384) for b in carriers.values()])
        else:
            s = shape(row[1]); k,b,tail = (s[x] for x in ('k','width','tail'))
            current = tuple(row[7:])
            key = str(row[0])+':'+str(row[1])
            fixed.setdefault(key,current)
            require(fixed[key] == current,'fixed lifecycle-cell workspaces')
            disjoint([(current[0],(k-1)*b+tail),(current[1]-16,32*3*k*b+32),
                      (current[2],32),(current[3],32*b),(current[4],32*b)])
        panel.append(row)
        if row[4] == 17:
            contrasts = []
            for pos in range(2,18,2):
                times = {SIDES[r[4]]^r[3]:r[duration] for r in panel[pos:pos+2]}
                require(set(times) == {0,1},'paired sides')
                contrasts.append(math.log(times[1])-math.log(times[0]))
            groups.setdefault(tuple(row[1:4]),[]).append(math.fsum(contrasts)/8)
            panel = []
    require(next(reader,None) is None and not panel and len(groups) == (192 if controlled else 1380),'exact cohort')
    stats, bad_aa, bad_cb, bad_wh1 = [], [], [], []
    for (cell,pair,order), logs in sorted(groups.items()):
        ci = confidence(logs)
        control = pair < (8 if controlled else 3)
        primary = (pair >= 8 and pair%4 != 0) if controlled else (pair == 3 and shape(cell)['k'] == 8 and
            shape(cell)['width'] == shape(cell)['tail'] == 1280 and cell%3 == 0)
        wh1 = not controlled and pair == 4
        passed = (1/1.02 < ci['lower95'] and ci['upper95'] < 1.02) if control else ci['upper95'] < (1 if primary or wh1 else 1.02)
        if not passed:
            (bad_aa if control else bad_wh1 if wh1 else bad_cb).append([cell,pair,order])
        stats.append(dict(cell=cell,pair=pair,order=order,control=control,primary=primary,wh1=wh1,passed=passed,**ci))
    return dict(outcome='CONTROL_FAIL' if bad_aa else 'FAIL' if bad_cb else 'PASS',
                failed_controls=bad_aa,failed_candidate=bad_cb,candidate_not_proven_faster_than_WH1=bad_wh1,
                statistics=stats,rows=count,total_work_ns=total,fixed_addresses=fixed,carriers=carriers,
                production_promotion_claimed=False,recovery_rate_claimed=False)


def trace_analyze(lines):
    reader = csv.reader(lines)
    require(next(reader,None) == TRACE,'trace header')
    bases = []
    for k in (3,5,8):
        for width in (64,1280):
            for policy in range(2):
                for order in range(2):
                    spans = []
                    for slot in range(2):
                        arm = slot^order
                        for allocation in range(3):
                            r = numeric(next(reader,None),14)
                            extra = 63 if arm and width == 1280 else 0
                            size = (296,k*width+extra,48+8*k)[allocation]
                            view = r[8]+((64-r[8]%64)%64 if allocation == 1 and extra else 0)
                            require(r[:8] == [k,width,policy,order,arm,allocation,size,int(allocation == 1)],'trace allocation contract')
                            require(r[9:] == [r[8]%64,r[8]%4096,view,view%64,view%4096],'trace raw/effective addresses')
                            spans.append((r[8],size))
                            if allocation == 1:
                                require(r[8] <= view and view+k*width <= r[8]+size,'basis within owner')
                                bases.append(r)
                    disjoint(spans)
    require(next(reader,None) is None,'trace extent')
    return dict(rows=144,bases=bases,timing_claimed=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True); f.write('\n')


def decisions(results):
    require(set(results) == set(WORKERS),'all worker results required')
    timed = {k:r['outcome'] for k,r in results.items() if not k.startswith('trace')}
    require(all(x in ('PASS','FAIL','CONTROL_FAIL') for x in timed.values()),'valid layer outcomes')
    return dict(outcome='CONTROL_FAIL' if 'CONTROL_FAIL' in timed.values() else 'FAIL' if 'FAIL' in timed.values() else 'PASS',
                layers=timed,production_promotion_claimed=False,recovery_rate_claimed=False)


def run(build):
    build = build.resolve(strict=True)
    cache = (build/'CMakeCache.txt').read_text().splitlines()
    require('SANITIZE:BOOL=OFF' in cache and 'CMAKE_HOME_DIRECTORY:INTERNAL='+str(HERE) in cache,'native observer build')
    libraries = [Path(x.split('=',1)[1]).resolve(strict=True) for x in cache if x.startswith('DSO_DIR:PATH=')]
    require(len(libraries) == 1,'one DSO directory')
    libraries = libraries[0]
    dsos = [libraries/('lib'+arm+'.so') for arm in ('baseline','candidate')]
    require(tuple(map(digest,dsos)) == LIBRARY_SHA,'exact neutral-qualified native DSOs')
    require(not any(k in os.environ for k in ENV),'unmodified loader/allocator environment')
    sources = [ROOT/p for p in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
               if Path(p).suffix in ('.h','.cpp','.c','.inc','.cmake') or Path(p).name == 'CMakeLists.txt']
    sources += [p for p in HERE.iterdir() if p.is_file()]
    for p in HERE.iterdir():
        if p.is_file():
            committed = subprocess.check_output(['git','show','HEAD:'+str(p.relative_to(ROOT))],cwd=ROOT)
            require(committed == p.read_bytes(),'commit final screen before launch')
    sources += list((HERE.parent/'Wh2AlignedSmallBasis').glob('*'))
    required = [HERE/p for p in ('CMakeLists.txt','NaturalMain.inc','ControlledMain.inc','Screen.py','test_Screen.py','README.md')]
    required += dsos+[libraries/p for p in ('candidate/WirehairV2Profile.cpp','CMakeCache.txt','build.ninja','compile_commands.json','Testing/Temporary/LastTest.log')]
    required += [build/p for p in ('natural','controlled','trace','Natural.cpp','Controlled.cpp','Trace.cpp','CMakeCache.txt','build.ninja','compile_commands.json')]
    require(all(p.is_file() for p in required),'all required source/build inputs exist')
    sources += required
    pins = {str(p):digest(p) for p in sources if p.is_file()}
    OUTPUT.mkdir(mode=0o700)
    write(OUTPUT/'claim.json',dict(protocol='wh2-aligned-small-basis-screen-r0',pins=pins,
          head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          cpu=50,replicates=12,natural_batch=32,controlled_batch=64,scope='bounded screen only; no production/all-K qualification'))
    results = {}
    for name in WORKERS:
        trace = name.startswith('trace')
        executable = 'trace' if trace else 'controlled' if name.startswith('controlled') else 'natural'
        mode = ('reverse' if name.endswith('reverse') else 'normal') if trace else ('run-reverse' if name.endswith('reverse') else 'run')
        require(all(digest(p) == h for p,h in pins.items()),'input pins before launch')
        with (OUTPUT/(name+'.csv')).open('xb') as out, (OUTPUT/(name+'.stderr')).open('xb') as err:
            try:
                code = subprocess.run([str(build/executable)]+list(map(str,dsos))+[mode],stdout=out,stderr=err,timeout=150,check=False).returncode
            except subprocess.TimeoutExpired:
                code = 'TIMEOUT'
        require(code == 0 and not (OUTPUT/(name+'.stderr')).stat().st_size,'worker failed: '+name+': '+str(code))
        require(all(digest(p) == h for p,h in pins.items()),'input pins after launch')
        with (OUTPUT/(name+'.csv')).open() as f:
            r = trace_analyze(f) if trace else analyze(f,executable == 'controlled')
        r.update(exit=code,raw_sha256=digest(OUTPUT/(name+'.csv')))
        results[name] = r
        write(OUTPUT/(name+'.json'),r)
        print(name,r.get('outcome','OBSERVED'),'controls',len(r.get('failed_controls',[])),'candidate',len(r.get('failed_candidate',[])),flush=True)
    write(OUTPUT/'decision.json',decisions(results))
    write(OUTPUT/'complete.json',{p.name:digest(p) for p in OUTPUT.iterdir() if p.is_file()})
    for p in OUTPUT.iterdir():
        p.chmod(0o400)


def replay(bundle):
    manifest = json.loads((bundle/'complete.json').read_text())
    expected = {'claim.json','decision.json'} | {name+suffix for name in WORKERS for suffix in ('.csv','.stderr','.json')}
    require(set(manifest) == expected and {p.name for p in bundle.iterdir()} == expected | {'complete.json'},'exact bundle members')
    require(all(digest(bundle/p) == h for p,h in manifest.items()),'bundle hashes')
    results = {}
    for name in WORKERS:
        old = json.loads((bundle/(name+'.json')).read_text())
        require(old['exit'] == 0 and not (bundle/(name+'.stderr')).stat().st_size,'successful worker')
        with (bundle/(name+'.csv')).open() as f:
            fresh = trace_analyze(f) if name.startswith('trace') else analyze(f,name.startswith('controlled'))
        fresh.update(exit=0,raw_sha256=digest(bundle/(name+'.csv')))
        require(json.loads(json.dumps(fresh)) == old,'exact recorded analysis')
        results[name] = fresh
    result = decisions(results)
    require(result == json.loads((bundle/'decision.json').read_text()),'aggregate decision')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run',type=Path,metavar='BUILD')
    group.add_argument('--analyze',type=Path)
    args = parser.parse_args()
    if args.run:
        run(args.run)
    else:
        print(json.dumps(replay(args.analyze),sort_keys=True))
