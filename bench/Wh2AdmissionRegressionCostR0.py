#!/usr/bin/env python3
"""Prospectively frozen admission pre/post DSO screen.

Build/neutral qualification and raw verifier only at this stage. No scientific
launch command: closed library provenance, receipt and controller are required.
Native DSO code is never claimed to be sanitizer-instrumented.
"""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import sys


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


N = sibling('admission_cost_neutral', 'Wh2AdmissionRegressionNeutral.py')
O = sibling('admission_cost_capture', 'Wh2K3OrdinaryCostR0.py')
A = N.A
ROOT = A.ROOT
PROTOCOL = 'wirehair.wh2.admission-regression-cost-r0'
OUTPUT = Path('/var/tmp/wh2-admission-regression-cost-r0')
CALLBACKS = 51840
PAIRS = ((0,0), (1,1), (0,1))
SIDES = (0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0)
ENV_KEYS = O.ENV_KEYS+('LD_AUDIT','LD_DEBUG')
FAMILIES = ('certified','small','wh1','k6')
CASES = tuple((0,k,b,2) for k in (2,3,4,6,128) for b in (2,1280))
CASES += tuple((0,3,b,1) for b in (2,1280))
CASES += tuple((1,3,b,p) for b in (2,1280) for p in (1,2))
CASES += tuple((2,3,b,2) for b in (2,1280))
CASES += tuple((3,6,b,2) for b in (2,1280))
NEUTRAL_DIR = Path('/tmp/wh2-admission-regression-neutral.DfWbHH')
NEUTRAL_SHA = ('a6a05ef144072dac836513781facb952a37753aab6309e2102e4605b0a8ba8f3',
               '2e03b86213d966536947349da28e7f8d252ca15334cbd93687679406672ab0b2')
NEW = ('bench/Wh2AdmissionRegressionCostR0.cpp','bench/Wh2AdmissionRegressionCostR0.py',
       'bench/test_Wh2AdmissionRegressionCostR0.py')
RAW_CAP, ERR_CAP = 192*1024**2, 65536


def batch(case):
    return 4 if case[1] == 128 else 128


def command(args):
    return subprocess.check_output(list(map(str,args)), cwd=ROOT, timeout=60)


def metadata():
    result = []
    for path, digest in N.LIBRARIES:
        raw = A.read_regular(path,4*1024**2)
        A.exact(A.sha(raw),digest,'qualified versioned DSO')
        elf = N.Elf(raw)
        A.exact(len(elf.slots),6,'exact internal public GOT roster')
        runtime_slots = []
        for i,section in enumerate(elf.sections):
            if section[1] != 4:
                continue
            symbol_section = elf.sections[section[6]]
            names = elf.body(symbol_section[6])
            symbols = list(struct.iter_unpack('<IBBHQQ',elf.body(section[6])))
            for offset,info,addend in struct.iter_unpack('<QQq',elf.body(i)):
                if info & 0xffffffff != 7:
                    continue
                name_offset,_,_,symbol_index,_,_ = symbols[info>>32]
                name = names[name_offset:names.index(b'\0',name_offset)].decode('ascii')
                if name.startswith('wirehair_'):
                    continue
                A.require(symbol_index==0 and addend==0 and offset%8==0 and elf.allocated(offset,8,3),
                          'allocated external runtime GOT')
                runtime_slots.append(dict(name=name,offset=offset))
        runtime_slots.sort(key=lambda s:s['name'])
        A.exact(len(runtime_slots),37,'complete external runtime GOT roster')
        A.exact(len({s['name'] for s in runtime_slots}),37,'unique external runtime imports')
        result.append(dict(path=str(path),sha256=digest,
                           exports=[dict(name=n,offset=elf.symbol(n,2)[0]) for n in sorted(elf.exports)],
                           slots=[dict(name=n,offset=p) for n,p in sorted(elf.slots)],
                           runtime_slots=runtime_slots,
                           context=elf.symbol('GF256Ctx',1)[0],
                           context_bytes=elf.symbol('GF256Ctx',1)[1],
                           getter=elf.symbol('gf256_get_active_x86_cpu_features',2)[0]))
    A.exact([s['name'] for s in result[0]['runtime_slots']],
            [s['name'] for s in result[1]['runtime_slots']],'identical runtime import names')
    return result


def bindings_header(meta):
    def entries(values):
        return '{{'+','.join('{'+json.dumps(v['name'])+','+str(v['offset'])+'}' for v in values)+'}}'
    rows = []
    for lib in meta:
        rows.append('{'+','.join((json.dumps(lib['path']),json.dumps(lib['sha256']),
                    entries(lib['exports']),entries(lib['slots']),entries(lib['runtime_slots']),str(lib['context']),
                    str(lib['context_bytes']),str(lib['getter'])))+'}')
    return ('// Generated from authenticated exact native ELF metadata.\n'
            'const LibrarySpec library_specs[2]={'+',\n'.join(rows)+'};\n').encode()


def roster():
    index = 0
    for r in range(12):
        for s in range(2):
            for ws in range(len(CASES)):
                for ms in range(2):
                    for cs in range(3):
                        order, which, metric = (r+s)%2, (r+s+ws)%len(CASES), (r+s+ws+ms)%2
                        comparison = (2*r+s+ws+metric+cs)%3
                        for p in range(18):
                            phase = r+12*(r%4) if p<2 else (r+6*((p-2)//8))%12+12*(((p-2)%8)//2)
                            yield [index,r,order,which,metric,comparison,p,
                                   PAIRS[comparison][SIDES[p]^order],(2*phase+1)*1000000//96]
                            index += 1


def statistics(records):
    A.exact(len(records),CALLBACKS,'complete statistical cohort')
    for row,coordinate in zip(records,roster()):
        A.exact(row['coordinate'],coordinate,'complete fixed statistical chronology')
    groups = {}
    for start in range(0,CALLBACKS,18):
        panel, contrasts = records[start:start+18], []
        for j in range(8):
            values = {}
            for row in panel[2+2*j:4+2*j]:
                c, clocks = row['coordinate'], row['observation']['clocks']
                side = SIDES[c[6]]^c[2]
                A.require(side not in values,'unique paired logical side')
                values[side] = A.integer(clocks[3])-A.integer(clocks[2])
            A.require(set(values)=={0,1} and min(values.values())>0,'positive paired durations')
            contrasts.append(math.log(values[1])-math.log(values[0]))
        c = panel[0]['coordinate']
        groups.setdefault((c[3],c[4],c[5],c[2]),[]).append(math.fsum(contrasts)/8)
    A.exact(len(groups),240,'all separate cost cells')
    results, controls, regressions, uncertain = [], [], [], []
    bound = math.log1p(.02)
    for key, values in sorted(groups.items()):
        A.exact(len(values),12,'t11 replicate count')
        which, metric, comparison, order = key
        estimate = A.confidence(values)
        lo, hi = estimate['lower95_log'], estimate['upper95_log']
        item = dict(case=list(CASES[which]),metric=metric,comparison=comparison,order=order,
                    estimate=estimate,replicate_logs=values)
        if comparison<2:
            passed = -bound<lo and hi<bound
            item['control_pass'] = passed
            if not passed: controls.append(list(key))
        else:
            item.update(resolved_slowdown=lo>0,uncertainty_bounded=hi<bound,
                        resolved_improvement=hi<0,screen_pass=lo<=0 and hi<bound)
            if lo>0: regressions.append(list(key))
            if hi>=bound: uncertain.append(list(key))
        results.append(item)
    return dict(outcome='CONTROL_FAIL' if controls else 'REGRESSION' if regressions else
                'INCONCLUSIVE' if uncertain else 'PASS', statistics=results,
                failed_controls=controls, resolved_regressions=regressions, uncertain=uncertain,
                all_K_claimed=False, WH1_speed_qualified=False, static_speed_qualified=False,
                recovery_rate_claimed=False, production_promotion_claimed=False)


def prior_records():
    records = None
    for order, digest in zip(('old-new','new-old'),NEUTRAL_SHA):
        raw = A.read_regular(NEUTRAL_DIR/('qualified-'+order+'.json'),1024*1024)
        A.exact(A.sha(raw),digest,'prior independently exercised neutral DSO records')
        result = A.decode(raw)
        A.exact(result['order'],order,'prior load order')
        if records is not None:
            A.exact(result['records'],records,'prior complete load-order identity')
        records = result['records']
    A.exact(len(records),216,'prior neutral roster')
    mapped = {tuple(r['case']):r for r in records}
    A.exact(len(mapped),216,'unique prior neutral cases')
    return mapped


def verify_bindings(bindings, meta):
    A.exact(len(bindings),2,'two DSO bindings')
    for bound, lib in zip(bindings,meta):
        A.exact(set(bound),{'base','context','exports','slots','providers','runtime_targets'},'binding schema')
        base = A.integer(bound['base'],1)
        A.exact(bound['context'],base+lib['context'],'private GF context')
        A.exact(bound['exports'],[base+s['offset'] for s in lib['exports']],'all public symbol addresses')
        public = {s['name']:base+s['offset'] for s in lib['exports']}
        A.exact(bound['slots'],[public[s['name']] for s in lib['slots']],'all internal GOT targets')
        A.exact(len(bound['providers']),5,'runtime provider roster')
        for address in bound['providers']:
            A.integer(address,1)
        A.exact(len(bound['runtime_targets']),37,'actual runtime GOT roster')
        for address in bound['runtime_targets']:
            A.integer(address,1)
    A.exact(bindings[0]['providers'],bindings[1]['providers'],'common runtime providers')
    A.exact(bindings[0]['runtime_targets'],bindings[1]['runtime_targets'],'common actual runtime GOT providers')
    x,y = (b['context'] for b in bindings)
    A.require(x+meta[0]['context_bytes']<=y or y+meta[1]['context_bytes']<=x,'distinct GF storage')


def verify_header(header, order, claim, meta):
    A.exact(set(header),{'type','protocol','claim','load_order','identity_hex','prelude','bindings','fixtures'},
            'header schema')
    A.exact((header['type'],header['protocol'],header['claim'],header['load_order']),
            ('header',PROTOCOL,claim,order),'header identity')
    encoded = header['identity_hex']
    A.require(type(encoded) is str and len(encoded)==2*O.TARGET_BYTES and
              re.fullmatch('[0-9a-f]+',encoded) is not None,'binary target identity')
    A.exact(A.sha(bytes.fromhex(encoded)),O.TARGET_SHA,'frozen CPU50 identity')
    O.clocks(header['prelude'],None)
    verify_bindings(header['bindings'],meta)
    prior = prior_records()
    A.exact(len(header['fixtures']),len(CASES),'full frozen case roster')
    for f, c in zip(header['fixtures'],CASES):
        family,k,b,policy = c
        A.exact(set(f),{'case','batch','source','arms'},'fixture schema')
        A.exact(f['case'],list(c),'case identity'); A.exact(f['batch'],batch(c),'fixed batch')
        A.exact(f['source'],bytes((37*i+i//11)%256 for i in range(k*b)).hex(),'fixed source')
        A.exact(len(f['arms']),2,'two arms'); A.exact(f['arms'][0],f['arms'][1],'exact pre/post fixture')
        old = prior[(FAMILIES[family],k,b,b,policy)]
        for arm in f['arms']:
            A.exact(set(arm),{'profile','packets','steps'},'arm schema')
            A.exact(arm['profile'],old['profile'],'prior exact descriptor')
            encoded = arm['packets']
            A.require(type(encoded) is str and len(encoded)==2*(k+14)*b and
                      re.fullmatch('[0-9a-f]+',encoded) is not None,'complete packet bytes')
            packet_bytes = bytes.fromhex(encoded)
            A.exact([A.sha(packet_bytes[j*b:(j+1)*b]) for j in range(k+14)],
                    old['packets'],'prior independent C ABI packet checks')
            A.exact(arm['steps'],old['first'],'prior own first-success endpoint')
            A.exact(old['ids'],list(range(k+8))+[0xffffffff-2*j for j in range(6)],'frozen packet IDs')



def verify_rows(rows, claim, order, meta):
    A.exact(len(rows),CALLBACKS+2,'whole raw cohort')
    header,footer = rows[0],rows[-1]
    verify_header(header,order,claim,meta)
    previous = header['prelude']; work = 0
    for row,coordinate in zip(rows[1:-1],roster()):
        A.exact(set(row),{'type','coordinate','ready','target','wait','observation','counts',
                         'addresses','address_count','complete','checked'},'record schema')
        A.exact(row['type'],'record','record type')
        A.exact(row['coordinate'],coordinate,'fixed chronology')
        _,_,_,which,metric,_,_,arm,q = coordinate
        c = CASES[which]; cycles = batch(c)
        ready,target = A.integer(row['ready']),A.integer(row['target'])
        A.require(previous['clocks'][5]<=ready and target==ready+q,'relative completion delay')
        A.exact(len(row['wait']),4,'wait shape')
        w0,c0,w1,c1 = [A.integer(v) for v in row['wait']]
        A.require(ready<=w0<=w1 and w1>=target and previous['clocks'][4]<=c0<=c1,'wait chronology')
        O.clocks(row['observation'],previous)
        observed = row['observation']['clocks']
        A.require(observed[0]>=w1 and observed[1]>=c1 and observed[2]>=target,'delayed WORK start')
        previous = row['observation']; work += observed[3]-observed[2]
        steps = header['fixtures'][which]['arms'][arm]['steps'] if metric else 0
        A.exact(row['counts'],[0 if metric else cycles,0 if metric else cycles*(c[1]+14),
                              cycles if metric else 0,cycles*steps,cycles if metric else 0,cycles],
                'every attempted API call')
        A.exact(row['address_count'],cycles,'full handle roster')
        A.exact(len(row['addresses']),128,'fixed address transport')
        for j,address in enumerate(row['addresses']):
            if j<cycles:
                A.integer(address,1)
            else:
                A.exact(address,0,'unused address slot')
        A.exact(row['complete'],True,'all WORK complete')
        A.exact(row['checked'],True,'all outputs checked')
    A.exact(footer,dict(type='footer',complete=True,records=CALLBACKS,work_ns=work),'terminal complete footer')
    A.require(work<=150000000000,'inner WORK cap')
    return statistics(rows[1:-1])


def verify(raw, claim, order, meta):
    A.require(0<len(raw)<=RAW_CAP and raw.endswith(b'\n'),'complete bounded raw stream')
    return verify_rows([A.decode(line) for line in raw.splitlines()],claim,order,meta)


def combine(results):
    A.exact(len(results),2,'both separate load orders')
    outcomes = [r['outcome'] for r in results]
    A.require(all(o in ('PASS','CONTROL_FAIL','REGRESSION','INCONCLUSIVE') for o in outcomes),'terminal screen outcomes')
    outcome = next((o for o in ('CONTROL_FAIL','REGRESSION','INCONCLUSIVE') if o in outcomes),'PASS')
    return dict(protocol=PROTOCOL,outcome=outcome,load_orders=results,
                shared_preserved_path_screen_pass=outcome=='PASS',
                WH1_speed_qualified=False,static_speed_qualified=False,all_K_claimed=False,
                recovery_rate_claimed=False,production_promotion_claimed=False)


def build(mode, output):
    A.require(mode in ('native','asan-driver'),'explicit native or instrumented driver build')
    A.require(output.is_absolute(),'absolute external build')
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT not in output.parents and output!=ROOT and not output.exists() and
              not output.is_symlink(),'fresh external build')
    A.require(all(os.environ.get(k) is None for k in ENV_KEYS),'clean allocator/loader environment')
    A.exact(command(['/usr/bin/c++','-dumpfullversion']).decode().strip(),'13.3.0','frozen GNU compiler')
    if mode=='asan-driver':
        A.exact(os.environ.get('ASAN_OPTIONS'),'detect_leaks=1:detect_stack_use_after_return=1','instrumented-driver ASAN policy')
        A.exact(os.environ.get('UBSAN_OPTIONS'),'halt_on_error=1','instrumented-driver UBSAN policy')
    meta = metadata()
    output.mkdir(mode=0o700)
    A.publish(output/'AdmissionLibraryBindings.h',bindings_header(meta))
    A.publish(output/'library-metadata.json',A.canonical(meta))
    flags = ['-std=c++11','-Wall','-Wextra','-Wpedantic','-Werror','-fno-lto','-fPIC',
             '-DWIREHAIR_STATIC=1','-DWH2_ADMISSION_REGRESSION_NEUTRAL='+str(int(mode!='native')),
             '-I'+str(ROOT),'-I'+str(ROOT/'bench'),'-I'+str(ROOT/'include'),'-I'+str(output)]
    flags += ['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer'] if mode=='asan-driver' else ['-O3','-g1']
    sources = [ROOT/NEW[0]]+[ROOT/'bench'/n for n in
              ('Wh2FrozenTrace.cpp','Wh2PublicBorrowedTargetIdentity.cpp','Wh2RdpruTargetIdentityV2.cpp')]
    dependencies = {ROOT/n for n in NEW}
    dependencies.update((ROOT/'bench/Wh2AdmissionRegressionNeutral.py', ROOT/'bench/Wh2K3OrdinaryCostR0.py',
                         ROOT/'bench/Wh2AlignedIntermediateCostR0.py'))
    dependencies.update(p for p,_ in N.LIBRARIES)
    dependencies.update(NEUTRAL_DIR/('qualified-'+o+'.json') for o in ('old-new','new-old'))
    objects, commands = [], []
    for source in sources:
        obj = output/(source.stem+'.o'); dep = output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-MD','-MF',str(dep),'-c',str(source),'-o',str(obj)]
        command(args); commands.append(args); objects.append(obj)
        dependencies.update(Path(p).resolve(strict=True) for p in
                            shlex.split(dep.read_text().replace('\\\n','').split(': ',1)[1]))
    exe = output/'cost_worker'
    args = ['/usr/bin/c++','-fno-lto','-no-pie']
    if mode=='asan-driver': args += ['-fsanitize=address,undefined']
    args += list(map(str,objects))+['-ldl','-pthread','-Wl,-Map,'+str(output/'link.map'),'-o',str(exe)]
    command(args); commands.append(args)
    names = [line.split()[-1] for line in command(['/usr/bin/nm','-g',exe]).decode().splitlines() if line.split()]
    A.require(not any(n.startswith(('wirehair_','gf256_')) or n=='GF256Ctx' for n in names),'no linked Wirehair/GF runtime')
    text = command(['/usr/bin/nm','-C',exe]).decode().splitlines()
    A.exact(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line for line in text),1,'single common WORK')
    for name in ('c++','cc','as','ld','nm'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for name in ('cc1','cc1plus','collect2'):
        dependencies.add(Path(command(['/usr/bin/c++','-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    for target in [exe,Path(sys.executable)]+[p for p,_ in N.LIBRARIES]:
        dependencies.update(Path(w).resolve(strict=True) for w in command(['/usr/bin/ldd',target]).decode().split()
                            if w.startswith('/'))
    dependencies.add(Path(sys.executable).resolve(strict=True))
    for order,name in enumerate(('old-new','new-old')):
        A.publish(output/('neutral-'+name+'.txt'),command([exe,'--neutral',name]))
        raw = command([exe,'--neutral-fixtures',name])
        verify_header(A.decode(raw),order,'0'*64,meta)
        A.publish(output/('fixtures-'+name+'.json'),raw)
    rejections = []
    for tail in ([],['--worker'],['--worker','0'*64,'bad-order'],['--worker','0'*64,'old-new'],
                 ['--neutral','bad-order'],['--neutral-fixtures','bad-order'],
                 ['--neutral','old-new','extra'],['--unexpected'],['--neutral-target','old-new']):
        result = subprocess.run([str(exe)]+tail,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,timeout=10)
        A.require(result.returncode==1 and result.stdout==b'' and result.stderr.startswith(b'INVALID:'),'negative worker CLI')
        rejections.append(dict(arguments=tail,returncode=result.returncode,stderr=result.stderr.decode()))
    A.publish(output/'negative-cli.json',A.canonical(rejections))
    manifest = dict(protocol=PROTOCOL,mode=mode,commands=commands,
                    environment={k:os.environ.get(k) for k in ENV_KEYS+('ASAN_OPTIONS','UBSAN_OPTIONS')},
                    inputs=[O.pin(p) for p in sorted(dependencies)],
                    artifacts=[O.pin(p) for p in sorted(output.iterdir())],
                    scientific_launch=False,library_source_provenance_closed=False,
                    sanitized_library_code=False)
    A.publish(output/'manifest.json',A.canonical(manifest))
    print(json.dumps(dict(mode=mode,executable=str(exe),scientific_launch=False,
                          library_source_provenance_closed=False)))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    b = sub.add_parser('build')
    b.add_argument('mode',choices=('native','asan-driver')); b.add_argument('output',type=Path)
    args = parser.parse_args()
    build(args.mode,args.output)
