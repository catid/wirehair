#!/usr/bin/env python3
"""One prospective unit-diagonal decoder lifecycle screen; no all-K qualification."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import Launch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUTPUT = Path('/var/tmp/wh2-small-unit-diagonal-screen-r0')
QUALIFIED = Path('/tmp/wh2-small-unit-diagonal.MfPjFJaS')
NEUTRAL_SHA = '93bfff5491f451c15ff8be50ec4453bef8770c0dbbc24e290ad48caa1822f532'
OBSERVER_MANIFEST = Path('/tmp/wh2-small-unit-diagonal-screen.tSZm35BW/OBSERVER.sha256')
OBSERVER_SHA = 'ca5fa7ff12828b0c52f03470b95cf4a9914950b03ac0923fcb934848511fc3e3'
LIBRARY_SHA = ('bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe',
               'af37e1286af35fe7695e9d423aafd3d4cba2ea7127ef432ce89f807ef990d0ba')
HEADER = 'rep,cell,pair,order,position,arm,ns,steps,source,output,descriptor,descriptor_out,low,distant'.split(',')
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
PAIRS = ((0,0),(1,1),(2,2),(0,1),(2,1))
WORKERS = ('normal','reverse')
ENV = ('LD_PRELOAD','LD_LIBRARY_PATH','LD_AUDIT','GLIBC_TUNABLES','MALLOC_PERTURB_',
       'MALLOC_TRIM_THRESHOLD_','MALLOC_MMAP_THRESHOLD_','MALLOC_TOP_PAD_',
       'ASAN_OPTIONS','UBSAN_OPTIONS')
BATCH, REPS, CELLS = 64, 12, 210
CAP = 240_000_000_000


def require(ok, why):
    if not ok:
        raise ValueError(why)


def shape(cell):
    require(type(cell) is int and 0 <= cell < CELLS,'valid cell')
    fixture = cell//3
    route, variant = fixture//10, (fixture//2)%5
    width = (2,2,64,1280,1280)[variant]
    return dict(k=(3,3,5,8,3,16,3)[route],width=width,tail=1 if variant in (1,4) else width,
                policy=fixture%2,route=route,metric=cell%3,small=route < 4 or route == 6)


def roster():
    for rep in range(REPS):
        for slot in range(CELLS):
            cell = (slot+rep*17)%CELLS
            for ps in range(5):
                pair = (ps+rep+cell)%5
                for order in range(2):
                    for pos in range(18):
                        yield [rep,cell,pair,order,pos,PAIRS[pair][SIDES[pos]^order]]


def numeric(raw):
    require(raw is not None and len(raw) == len(HEADER) and
            all(x.isascii() and x.isdecimal() for x in raw),'complete numeric row')
    return list(map(int,raw))


def disjoint(ranges):
    ranges = sorted(ranges)
    require(all(0 < a < a+n < 2**64 for a,n in ranges) and
            all(a+n <= b for (a,n),(b,_) in zip(ranges,ranges[1:])),'valid disjoint public storage')


def confidence(logs):
    require(len(logs) == REPS and all(math.isfinite(x) for x in logs),'12 finite replicates')
    mean = math.fsum(logs)/REPS
    half = 2.200985160082949*math.sqrt(math.fsum((x-mean)**2 for x in logs)/(11*12))
    return dict(ratio=math.exp(mean),lower95=math.exp(mean-half),upper95=math.exp(mean+half))


def summarize(groups):
    require(set(groups) == {(c,p,o) for c in range(CELLS) for p in range(5) for o in range(2)},'exact interval roster')
    stats, aa, retention, wh1 = [], [], [], []
    for (cell,pair,order), logs in sorted(groups.items()):
        ci = confidence(logs)
        passed = (1/1.02 < ci['lower95'] and ci['upper95'] < 1.02) if pair < 3 else ci['upper95'] < (1.02 if pair == 3 else 1)
        if not passed:
            (aa if pair < 3 else retention if pair == 3 else wh1).append([cell,pair,order])
        stats.append(dict(cell=cell,pair=pair,order=order,passed=passed,**ci))
    required_wh1 = [key for key in wh1 if shape(key[0])['small']]
    primary = []
    for metric in (1,2):
        cells = [c for c in range(CELLS) if shape(c)['small'] and c%3 == metric]
        require(len(cells) == 50,'all 50 affected small fixtures in decoder primary')
        for order in range(2):
            # Keep covariance: one equal-fixture aggregate per replicate.
            logs = [math.fsum(groups[c,3,order][r] for c in cells)/50 for r in range(REPS)]
            ci = confidence(logs)
            primary.append(dict(metric=metric,order=order,passed=ci['upper95'] < 1,**ci))
    return dict(outcome='CONTROL_FAIL' if aa else 'FAIL' if retention or required_wh1 or any(not p['passed'] for p in primary) else 'PASS',
                failed_controls=aa,failed_retention=retention,primary=primary,statistics=stats,
                candidate_not_proven_faster_than_WH1=wh1,failed_required_wh1=required_wh1,
                production_promotion_claimed=False,recovery_rate_claimed=False)


def analyze(lines):
    reader = csv.reader(lines)
    require(next(reader,None) == HEADER,'header')
    panel, groups, fixed, steps, total, count = [], {}, {}, {}, 0, 0
    for expected in roster():
        row = numeric(next(reader,None))
        require(row[:6] == expected and 0 < row[6] < CAP,'chronology/duration')
        count += 1
        total += row[6]
        require(total < CAP,'aggregate WORK cap including warmups')
        s = shape(row[1]); k,b,tail = (s[x] for x in ('k','width','tail'))
        require(row[7] == 0 if s['metric'] == 0 else k <= row[7] <= 32,'first-success bound')
        # C/B must agree; WH1 keeps its own workload and may require more rows.
        key = str(row[1])+':'+str(int(row[5] == 2))
        steps.setdefault(key,row[7])
        require(steps[key] == row[7],'fixed arm first-success count')
        current = tuple(row[8:])
        key = str(row[0])+':'+str(row[1])
        fixed.setdefault(key,current)
        require(fixed[key] == current,'fixed lifecycle-cell public workspaces')
        disjoint([(current[0],(k-1)*b+tail),(current[1]-16,BATCH*3*k*b+32),
                  (current[2],32),(current[3],32),(current[4],32*b),(current[5],32*b)])
        panel.append(row)
        if row[4] == 17:
            contrasts = []
            for pos in range(2,18,2):
                times = {SIDES[r[4]]^r[3]:r[6] for r in panel[pos:pos+2]}
                require(set(times) == {0,1},'paired sides')
                contrasts.append(math.log(times[1])-math.log(times[0]))
            groups.setdefault(tuple(row[1:4]),[]).append(math.fsum(contrasts)/8)
            panel = []
    require(next(reader,None) is None and not panel,'exact row cohort')
    result = summarize(groups)
    result.update(rows=count,total_work_ns=total,fixed_addresses=fixed,first_success=steps)
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True); f.write('\n')


def checked_manifest(path, expected, relative_root=None):
    require(digest(path) == expected,'qualification manifest unchanged: '+str(path))
    pins = {str(path):expected}
    for line in path.read_text().splitlines():
        sha, separator, filename = line.partition('  ')
        require(separator and len(sha) == 64 and all(c in '0123456789abcdef' for c in sha) and filename,
                'canonical qualification manifest entry')
        member = Path(filename)
        if not member.is_absolute():
            require(relative_root is not None and '..' not in member.parts,'explicit manifest-relative root')
            member = relative_root/member
        filename = str(member.resolve(strict=True))
        require(filename not in pins,'unique qualification manifest member')
        require(digest(filename) == sha,'qualified artifact unchanged: '+filename)
        pins[filename] = sha
    require(len(pins) > 1,'nonempty qualification manifest')
    return pins


def decisions(results):
    require(set(results) == set(WORKERS),'both load orders required')
    outcomes = {k:r['outcome'] for k,r in results.items()}
    require(all(x in ('PASS','FAIL','CONTROL_FAIL') for x in outcomes.values()),'valid worker outcomes')
    return dict(outcome='CONTROL_FAIL' if 'CONTROL_FAIL' in outcomes.values() else 'FAIL' if 'FAIL' in outcomes.values() else 'PASS',
                load_orders=outcomes,production_promotion_claimed=False,recovery_rate_claimed=False)


def inputs(build):
    require(build == OBSERVER_MANIFEST.parent/'native','exact neutral-qualified observer build')
    require(subprocess.check_output(['ninja','-n'],cwd=build,text=True) == 'ninja: no work to do.\n',
            'observer build is current; rebuild only the observer before claiming a namespace')
    cache = (build/'CMakeCache.txt').read_text().splitlines()
    require('SANITIZE:BOOL=OFF' in cache and 'CMAKE_HOME_DIRECTORY:INTERNAL='+str(HERE) in cache and
            'DSO_DIR:PATH='+str(QUALIFIED/'native') in cache,'native observer and exact qualified library directory')
    dsos = [QUALIFIED/'native'/('lib'+arm+'.so') for arm in ('baseline','candidate')]
    require(tuple(map(digest,dsos)) == LIBRARY_SHA,'exact qualified native DSOs')
    manifest = QUALIFIED/'NEUTRAL.sha256'
    qualified_pins = checked_manifest(manifest,NEUTRAL_SHA,ROOT)
    observer_pins = checked_manifest(OBSERVER_MANIFEST,OBSERVER_SHA)
    require(not any(k in os.environ for k in ENV),'unmodified loader/allocator/sanitizer environment')
    sources = [ROOT/p for p in subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
               if Path(p).suffix in ('.h','.cpp','.c','.inc','.cmake') or Path(p).name == 'CMakeLists.txt']
    own = [HERE/p for p in ('CMakeLists.txt','Generate.py','Screen.py','Launch.py','test_Screen.py','test_Launch.py','test_Generate.py','README.md')]
    for p in own:
        require(subprocess.check_output(['git','show','HEAD:'+str(p.relative_to(ROOT))],cwd=ROOT) == p.read_bytes(),'commit final screen before launch')
    shared = [ROOT/p for p in ('bench/Wh2SmallPayload2/Screen.cpp','bench/Wh2SmallDormantCoreScreen/Main.inc')]
    required = own + shared + dsos + [manifest]
    required += [build/p for p in ('worker','Worker.cpp','CMakeCache.txt','build.ninja','compile_commands.json','Testing/Temporary/LastTest.log')]
    required += [QUALIFIED/'native'/p for p in ('candidate/WirehairSmallCore.h','CMakeCache.txt','build.ninja','compile_commands.json','Testing/Temporary/LastTest.log')]
    require(all(p.is_file() for p in required),'all required source/build/test inputs exist')
    pins = {str(p):digest(p) for p in sources+required}
    pins.update(qualified_pins)
    pins.update(observer_pins)
    return dsos,pins


def run(build):
    # Check the spent namespace before any build inspection or subprocess.
    require(not OUTPUT.exists(),'namespace already spent')
    build = build.resolve(strict=True)
    dsos,pins = inputs(build)
    OUTPUT.mkdir(mode=0o700)
    write(OUTPUT/'claim.json',dict(protocol=OUTPUT.name,pins=pins,
          head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          cpu=50,replicates=REPS,batch=BATCH,cells=CELLS,worker_cap_seconds=240,controller_timeout_seconds=270,
          launcher_limits=Launch.limits(),launcher_environment=Launch.ENVIRONMENT,
          scope='bounded lifecycle screen only; no production/all-K or recovery-rate qualification'))
    results = {}
    name = 'prelaunch'
    try:
        for name in WORKERS:
            require(all(digest(p) == h for p,h in pins.items()),'input pins before launch')
            mode = 'run' if name == 'normal' else 'run-reverse'
            captured = Launch.capture([str(build/'worker')]+list(map(str,dsos))+[mode],
                                      OUTPUT/(name+'.csv'),OUTPUT/(name+'.stderr'))
            write(OUTPUT/(name+'.exit.json'),captured)
            successful_capture(OUTPUT,name,captured)
            code = captured['exit']
            require(all(digest(p) == h for p,h in pins.items()),'input pins after launch')
            with (OUTPUT/(name+'.csv')).open() as f:
                r = analyze(f)
            r.update(exit=code,raw_sha256=digest(OUTPUT/(name+'.csv')))
            results[name] = r
            write(OUTPUT/(name+'.json'),r)
            print(name,r['outcome'],'controls',len(r['failed_controls']),'retention',len(r['failed_retention']),
                  'required WH1',len(r['failed_required_wh1']),flush=True)
        write(OUTPUT/'decision.json',decisions(results))
        write(OUTPUT/'complete.json',{p.name:digest(p) for p in OUTPUT.iterdir() if p.is_file()})
    except Exception as e:
        write(OUTPUT/'failed.json',dict(worker=name,error=str(e),namespace_spent=True))
        raise
    finally:
        for p in OUTPUT.iterdir():
            p.chmod(0o400)


def successful_capture(bundle, name, captured):
    require(captured['exit'] == 0 and captured['failure'] is None and
            0 < captured['wall_seconds'] <= Launch.WALL_SECONDS and
            captured['stdout_bytes'] == (bundle/(name+'.csv')).stat().st_size <= Launch.STDOUT_BYTES and
            captured['stderr_bytes'] == (bundle/(name+'.stderr')).stat().st_size == 0,
            'worker failed or capture invalid: '+name)


def replay(bundle):
    manifest = json.loads((bundle/'complete.json').read_text())
    expected = {'claim.json','decision.json'} | {name+suffix for name in WORKERS for suffix in ('.csv','.stderr','.exit.json','.json')}
    require(set(manifest) == expected and {p.name for p in bundle.iterdir()} == expected | {'complete.json'},'exact bundle members')
    require(all(digest(bundle/p) == h for p,h in manifest.items()),'bundle hashes')
    claim = json.loads((bundle/'claim.json').read_text())
    require(claim['protocol'] == OUTPUT.name and (claim['cpu'],claim['replicates'],claim['batch'],claim['cells']) == (50,12,64,210) and
            (claim['worker_cap_seconds'],claim['controller_timeout_seconds']) == (240,270) and claim['launcher_limits'] == Launch.limits() and
            claim['launcher_environment'] == Launch.ENVIRONMENT,'frozen protocol')
    results = {}
    for name in WORKERS:
        old = json.loads((bundle/(name+'.json')).read_text())
        require(old['exit'] == 0,'successful worker')
        successful_capture(bundle,name,json.loads((bundle/(name+'.exit.json')).read_text()))
        with (bundle/(name+'.csv')).open() as f:
            fresh = analyze(f)
        fresh.update(exit=0,raw_sha256=digest(bundle/(name+'.csv')))
        require(json.loads(json.dumps(fresh)) == old,'exact recorded analysis')
        results[name] = fresh
    result = decisions(results)
    require(result == json.loads((bundle/'decision.json').read_text()),'aggregate decision')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run',type=Path,metavar='NATIVE_OBSERVER_BUILD')
    group.add_argument('--analyze',type=Path,metavar='BUNDLE')
    args = parser.parse_args()
    if args.run:
        run(args.run)
    else:
        print(json.dumps(replay(args.analyze),sort_keys=True))
