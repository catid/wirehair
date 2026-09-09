#!/usr/bin/env python3
"""Installed WHV2 K5 recovery, six ownership-matched arms; retained cohort only."""
import argparse
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path
import re
import selectors
import subprocess
import struct
import sys
import time

HERE = Path(__file__).resolve().parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


U = sibling('_k5_recovery_build_io', 'Wh2K5PublicCostR0.py')
D = sibling('_k5_recovery_retained_data', 'Wh2K5NativeDataR0.py')
A = U.A
ROOT = U.ROOT
PROTOCOL = 'wirehair.wh2.k5-public-recovery-r0'
OUTPUT = Path('/var/tmp/wh2-k5-public-recovery-r0')
MODES = ('native', 'scalar', 'asan')
PROTOTYPE_RAW = Path('/var/tmp/wh2-k5-serialized-recovery-r0/native.raw.jsonl')
PROTOTYPE_SHA = 'd6951f2594434c739776d0c37e1e4ed3b44317d91a49819dc832c28fee837303'
WIDTHS = (2, 64, 1280)
ARMS = U.ARM_NAMES
K = 5
RAW_CAP, ERR_CAP = 64*1024**2, 65536
PAIRED = ((1,0),(1,2),(4,3),(4,5))
ENV_KEYS = U.ENV_KEYS + ('ASAN_OPTIONS', 'UBSAN_OPTIONS')
SANITIZERS = dict(ASAN_OPTIONS='detect_leaks=1:detect_stack_use_after_return=1:halt_on_error=1',
                  UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
SOURCES = ('bench/Wh2K5PublicRecoveryR0.cpp', 'bench/Wh2K5PublicRecoveryR0.py',
           'bench/test_Wh2K5PublicRecoveryR0.py', 'bench/Wh2K5PublicCostR0.py',
           'bench/Wh2AlignedIntermediateCostR0.py', 'bench/Wh2K5NativeDataR0.py',
           'bench/Wh2K3NativeDataR0.py', 'bench/Wh2NoncommutingRadixRunR0.py',
           'bench/Wh2K5PublicRecoveryR0.md', 'V2_WIRE_PROFILE.md', 'WH2_BORROWED_SOURCE_API.md')



@lru_cache(maxsize=1)
def retained():
    report = D.load_report()
    data = D.extract(report)
    return report, data


@lru_cache(maxsize=2)
def roster(neutral=False):
    report, data = retained()
    result = []
    if neutral:
        streams = [list(range(9)),list(range(5,14)),[0xffffffff-2*j for j in range(9)],[0]*9]
        for width in WIDTHS:
            for index, ids in enumerate(streams):
                result.append(dict(group=3,index=index,B=width,ids=ids))
        for index in list(range(12))+list(range(6132,6144)):
            row = report['fresh'][index]
            result.append(dict(group=4,index=index,B=row['B'],ids=row['ids']))
        A.exact(len(result),36,'bounded first/last-cell smoke roster')
    else:
        for group, name in enumerate(('fresh','hard')):
            for index, row in enumerate(report[name]):
                result.append(dict(group=group,index=index,B=row['B'],ids=row['ids'],
                                   schedule=row['schedule'],root=row['root']))
        for index, row in enumerate(data['history']):
            for width in row['widths']:
                result.append(dict(group=2,index=index,B=width,ids=row['ids']))
        A.exact(len(result),6273,'retained trace/width roster')
        A.exact(sum(len(e['ids']) for e in result),56234,'exact original packet horizons')
    return result


@lru_cache(maxsize=256)
def products(coefficient):
    return bytes(D.R.multiply(coefficient,b) for b in range(256))


@lru_cache(maxsize=1)
def powers():
    def product(a,b):
        result = [0]*25
        for i in range(K):
            for j in range(K):
                for k in range(K): result[K*i+j] ^= products(a[K*i+k])[b[K*k+j]]
        return tuple(result)
    matrices = []
    for feedback in ((121,110,207,198,31),(120,110,207,198,31)):
        m = [0]*25
        for i in range(K-1): m[K*(i+1)+i] = 1
        for i in range(K): m[K*i+K-1] = feedback[i]
        matrices.append(tuple(m))
    levels = [tuple(matrices)]
    for _ in range(31):
        a,b = levels[-1]
        levels.append((product(a,b),product(b,a)))
    return tuple(levels)


@lru_cache(maxsize=8192)
def coefficient(packet_id):
    row = (1,0,0,0,0)
    for bit in range(32):
        if packet_id & (1 << bit):
            matrix = powers()[bit][bin(packet_id >> (bit+1)).count('1') % 2]
            next_row = [0]*K
            for i in range(K):
                for j in range(K): next_row[i] ^= products(matrix[K*i+j])[row[j]]
            row = tuple(next_row)
    return row


@lru_cache(maxsize=3)
def message(width):
    A.require(width in WIDTHS,'frozen width')
    return bytes((37*i+i//11) % 256 for i in range(K*width))


@lru_cache(maxsize=32768)
def packet_hash(width,row):
    a,b,c,d,e = (products(v) for v in row)
    source = message(width)
    packet = bytes(a[source[j]] ^ b[source[width+j]] ^ c[source[2*width+j]] ^
                   d[source[3*width+j]] ^ e[source[4*width+j]] for j in range(width))
    return A.sha(packet)


@lru_cache(maxsize=32768)
def first_success(rows):
    # Independent inverse-free incremental elimination, no decoder state.
    basis = {}
    for index,incoming in enumerate(rows):
        row = list(incoming)
        for column in range(K):
            if not row[column]: continue
            if column not in basis:
                basis[column] = row
                if len(basis)==K: return index+1
                break
            pivot = basis[column]
            a,b = products(pivot[column]),products(row[column])
            row = [a[x] ^ b[y] for x,y in zip(row,pivot)]
    return 0


def candidate_first(ids):
    return first_success(tuple(coefficient(i) for i in ids))


def valid_hex(value, size):
    A.require(type(value) is str and len(value) == 2*size and re.fullmatch('[0-9a-f]+', value) is not None,
              'hex bytes')


def summarize(records,expected):
    A.exact(len(records),len(expected),'summary whole cohort')
    def paired(selected,historical=False):
        pairs = []
        for candidate,control in PAIRED:
            pair = dict(control=ARMS[control],candidate=ARMS[candidate],fixed=[],introduced=[])
            for overhead in range(5):
                fixed,introduced = [],[]
                for r,e in selected:
                    if historical and K+overhead>len(e['ids']): continue
                    x,y = r['arms'][control]['first'],r['arms'][candidate]['first']
                    control_failed,candidate_failed = x==0 or x>K+overhead,y==0 or y>K+overhead
                    key = [e['index'],e['B']]
                    if control_failed and not candidate_failed: fixed.append(key)
                    if candidate_failed and not control_failed: introduced.append(key)
                pair['fixed'].append(fixed); pair['introduced'].append(introduced)
            pairs.append(pair)
        return pairs
    def curves(selected):
        return [[sum(r['arms'][a]['first']==0 or r['arms'][a]['first']>K+overhead for r,_ in selected)
                 for overhead in range(5)] for a in range(6)]
    groups = []
    for group,name in enumerate(('retained_fresh','hard','historical_original_width')):
        selected = [(r,e) for r,e in zip(records,expected) if e['group']==group]
        groups.append(dict(group=name,cases=len(selected),
            failure_counts_oh0_to_4=curves(selected) if group!=2 else None,
            unresolved=[sum(r['arms'][a]['first']==0 for r,_ in selected) for a in range(6)],
            paired=paired(selected,group==2)))
    cells = []
    schedules = sorted({e['schedule'] for e in expected if e['group']==0})
    A.exact(len(schedules),4,'four retained loss schedules')
    for width in WIDTHS:
        for schedule in schedules:
            selected = [(r,e) for r,e in zip(records,expected) if e['group']==0 and e['B']==width and e['schedule']==schedule]
            A.exact(len(selected),512,'whole retained cell')
            cells.append(dict(B=width,schedule=schedule,cases=512,
                              failure_counts_oh0_to_4=curves(selected),paired=paired(selected)))
    fresh = groups[0]['failure_counts_oh0_to_4']
    wins = {ARMS[c]:all(fresh[c][0]<fresh[a][0] for cc,a in PAIRED if cc==c) for c in (1,4)}
    return dict(outcome='PASS',arms=list(ARMS),groups=groups,cells=cells,
                comparative_oh0_win=all(wins.values()),policy_oh0_wins=wins,
                candidate_cells_at_most_one_percent=all(100*cell['failure_counts_oh0_to_4'][a][0]<=512
                                                       for cell in cells for a in (1,4)),
                fresh_sample=False,speed_claimed=False,all_K_claimed=False,
                production_promotion_claimed=False)


def verify(raw,claim,mode,neutral=False):
    A.require(0<len(raw)<=RAW_CAP and raw.endswith(b'\n'),'bounded complete stream')
    rows = [A.decode(line) for line in raw.splitlines()]
    expected = roster(neutral)
    A.exact(len(rows),len(expected)+2,'whole recovery cohort')
    header,records,footer = rows[0],rows[1:-1],rows[-1]
    A.exact(header,dict(type='header',protocol=PROTOCOL,claim=claim,backend=mode,
                       scope='neutral' if neutral else 'retained',retained_raw_sha256=D.RAW_SHA,
                       features=[0]*4 if mode=='scalar' else [1]*4,
                       sources=[A.sha(message(b)) for b in WIDTHS]),'exact replay header')
    profiles = {}
    for row,e in zip(records,expected):
        A.exact(set(row),{'type','group','index','B','ids','arms'},'record schema')
        A.exact({k:row[k] for k in ('group','index','B','ids')},
                {k:e[k] for k in ('group','index','B','ids')},'retained exact chronology')
        A.exact(row['type'],'record','record type'); A.exact(len(row['arms']),6,'six ownership-matched API routes')
        width,ids = e['B'],e['ids']
        for arm,r in enumerate(row['arms']):
            A.exact(set(r),{'profile','packets','rows','feed','first','recoveries','counts','checked'},'arm schema')
            valid_hex(r['profile'],32)
            A.exact(r['profile'],profiles.setdefault((width,arm),r['profile']),'stable full descriptor')
            profile = bytes.fromhex(r['profile'])
            if arm in (0,3):
                A.exact(profile[:16],struct.pack('<4sHHQ',b'WHV2',1,32,0x4b295bbb47f4f9c9),'explicit certified WH2 descriptor')
                A.exact(profile[16:28],struct.pack('<QI',K*width,width),'certified WH2 dimensions')
                A.exact(profile[29:],bytes(3),'certified reserved bytes')
            elif arm in (1,4):
                A.exact(profile,struct.pack('<4sHHQQII',b'WHV2',1,32,0x80070c81bfe375f1,K*width,width,0),'installed sealed WHV2 K5 identity')
            else: A.exact(profile,bytes(32),'WH1 has no descriptor')
            valid_hex(r['rows'],K*len(ids))
            encoded_rows = bytes.fromhex(r['rows'])
            observed = tuple(tuple(encoded_rows[i:i+K]) for i in range(0,len(encoded_rows),K))
            if arm in (1,4):
                A.exact(observed,tuple(coefficient(i) for i in ids),'independent sealed K5 rows')
            A.exact(len(r['packets']),len(ids),'all retained packets')
            for packet in r['packets']: valid_hex(packet,32)
            A.exact(r['packets'],[packet_hash(width,rr) for rr in observed],'independent every-arm packet hashes')
            first = A.integer(r['first'],0,len(ids))
            A.require(first==0 or first>=K,'possible first success')
            A.exact(first,first_success(observed),'independent every-arm first-success rank')
            A.exact(r['feed'],[1]*(first-1)+[0] if first else [1]*len(ids),'all prefix feed statuses')
            A.exact(r['recoveries'],2 if first else 0,'two byte-exact recoveries')
            A.exact(r['counts'],[6,6*len(ids),1,first if first else len(ids),2 if first else 0,7],'whole attempted API ledger')
            A.exact(r['checked'],True,'payload/guard validation')
            if arm in (1,4) and not neutral:
                A.require(first!=0,'retained candidate must recover')
                if e['group']!=2:
                    ranks = retained()[0]['fresh' if e['group']==0 else 'hard'][e['index']]['ranks']
                    A.exact(first,next(K+j for j,v in enumerate(ranks) if v==K),'sealed retained prefix outcome')
        for a in range(3): A.exact(row['arms'][a],row['arms'][a+3],'policy-independent complete result')
    A.exact(footer,dict(type='footer',records=len(expected),checked=True),'complete footer')
    return (dict(outcome='PASS',neutral_cases=36) if neutral else summarize(records,expected)),records


def retained_parity(records):
    """Exact old outcomes, permitting only the explicitly versioned descriptor."""
    raw = A.read_regular(PROTOTYPE_RAW,RAW_CAP)
    A.exact(A.sha(raw),PROTOTYPE_SHA,'authenticated prototype recovery evidence')
    rows = [A.decode(line) for line in raw.splitlines()]
    A.exact(len(rows),6275,'whole prototype recovery roster')
    expected = rows[1:-1]
    for row in expected:
        for arm in (1,4):
            A.exact(row['arms'][arm]['profile'],struct.pack('<4sHHQQII',b'WHK5',1,32,
                    0x5748324b35544d31,K*row['B'],row['B'],0).hex(),'original prototype descriptor')
            row['arms'][arm]['profile'] = struct.pack('<4sHHQQII',b'WHV2',1,32,
                    0x80070c81bfe375f1,K*row['B'],row['B'],0).hex()
    A.exact(records,expected,'complete installed versus retained prototype parity')


def context_size(symbols):
    rows = [line.split() for line in symbols.decode().splitlines() if line.split() and line.split()[-1]=='GF256Ctx']
    A.require(len(rows)==1 and len(rows[0])==4 and rows[0][2]=='B','one sized GF context symbol')
    A.require(re.fullmatch('[0-9a-fA-F]+',rows[0][1]) is not None,'GF context size encoding')
    size = int(rows[0][1],16)
    A.require(0<size<1024*1024,'bounded GF context size')
    return size


def build(mode, output):
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT not in output.parents and output != ROOT and not output.exists() and not output.is_symlink(),
              'fresh external build')
    archives, dependencies, provenance = U.qualified_inputs(mode)
    frozen = {}
    dependencies.update(ROOT/s for s in SOURCES)
    dependencies.add(PROTOTYPE_RAW)
    dependencies.add(Path(sys.executable).resolve(strict=True))
    for name in ('c++','cc','as','ld','nm','ar','ldd','bash'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for compiler in ('cc1','cc1plus','collect2'):
        dependencies.add(Path(U.command(['/usr/bin/c++','-print-prog-name='+compiler]).decode().strip()).resolve(strict=True))
    U.freeze_inputs(dependencies,frozen)
    for record in provenance['members']:
        A.exact(frozen[Path(record['object']['path'])],record['object'],'unchanged producing object')
    for key in ('archive','qualification_log'):
        A.exact(frozen[Path(provenance[key]['path'])],provenance[key],'unchanged qualified evidence')
    output.mkdir(mode=0o700)
    A.publish(output/'qualified-library.json',A.canonical(provenance))
    U.command([sys.executable, '-I', '-B', '-S', HERE/'Wh2K5NativeDataR0.py', '--header', output/'Wh2K5NativeData.inc'])
    U.freeze_inputs(dependencies,frozen)
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_K5_PUBLIC_RECOVERY_BACKEND='+str(MODES.index(mode)),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    # gf256.h has a compile-ISA-dependent private context. ASAN's preserved
    # archive was built with -march=native; every direct table reader must match.
    flags += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-march=native'] if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar': flags.append('-DANDROID=1')
    commands, objects = [], []
    for source in (ROOT/SOURCES[0], HERE/'Wh2FrozenTrace.cpp'):
        obj, dep = output/(source.stem+'.o'), output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-M','-MT',str(obj),'-MF',str(dep),str(source)]
        U.command(args); commands.append(args)
        before = U.preprocessor_dependencies(dep.read_bytes(),obj)
        dependencies.update(before); U.freeze_inputs(dependencies,frozen)
        args = ['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)]
        U.command(args); commands.append(args); objects.append(obj)
        after = U.preprocessor_dependencies(dep.read_bytes(),obj)
        A.exact(sorted(after),sorted(before),'actual compiler dependency closure')
        U.freeze_inputs(dependencies,frozen)
    executable = output/'recovery_worker'
    args = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan': args.append('-fsanitize=address,undefined')
    args += list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    U.command(args); commands.append(args)
    # Include linker-loaded startup objects, archives and linker scripts, not
    # only the final executable and its dynamic runtime libraries.
    loads = [line[5:] for line in (output/'link.map').read_text().splitlines() if line.startswith('LOAD ')]
    A.require(loads and all(Path(p).is_absolute() for p in loads),'absolute linker input roster')
    for path in loads:
        resolved = Path(path).resolve(strict=True)
        if output not in resolved.parents:
            dependencies.add(resolved)
    symbols = U.command(['/usr/bin/nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    A.require(names.count('GF256Ctx') == names.count('gf256_init_') == 1, 'single GF runtime')
    A.require(all(names.count('wirehair_v2_'+name)==1 for name in
                  ('encoder_create_profile_id_with_options','decoder_create','encode','decode','recover','free')),
              'one actual installed WHV2 API for both profiles')
    A.require(not any(n.startswith('wh2_small_') for n in names),'no benchmark K5 boundary linked')
    A.require(all(names.count(n) == 1 for n in ('wirehair_encode', 'wirehair_decode', 'wirehair_free')), 'actual WH1 API')
    tool_binaries = [p for p in dependencies if os.access(str(p),os.X_OK) and
                     A.read_regular(p,256*1024*1024,installed=True)[:4]==b'\x7fELF']
    for target in tool_binaries+[executable, Path(sys.executable).resolve(strict=True)]:
        dependencies.update(Path(w).resolve(strict=True) for w in U.command(['/usr/bin/ldd', target]).decode().split() if w.startswith('/'))
    dependencies.add(Path(sys.executable).resolve(strict=True))
    contract = U.command([executable,'--contract'])
    context_bytes = context_size(U.command(['/usr/bin/nm','-S','--defined-only',executable]))
    A.exact(A.decode(contract),dict(K=5,arms=6,records=6273,cpu_seconds=120,wall_seconds=150,
            address_space_mib=384,gf_context_bytes=context_bytes,asan_shadow_exempt=mode=='asan',backend=mode,
            claim_path=str(OUTPUT/'CLAIM.json')),'compiled resource/launch contract')
    A.publish(output/'contract.json',contract)
    neutral_claim = output/'neutral-claim.json'; claimed = A.canonical(dict(protocol=PROTOCOL,neutral=True))
    A.publish(neutral_claim,claimed)
    A.exact(U.command([executable,'--neutral-claim',neutral_claim,A.sha(claimed)]),
            b'PASS claim authentication\n','positive claim authentication')
    negatives = []
    for argv in ([executable,'--neutral-claim',neutral_claim,'0'*64],
                 [executable,'--neutral-claim',output/'absent',A.sha(claimed)],
                 [executable,'--neutral-claim',neutral_claim,'g'*64],
                 [executable],[executable,'--unknown'],[executable,'--worker','0']):
        p = subprocess.run(list(map(str,argv)),cwd=ROOT,stdin=subprocess.DEVNULL,
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10)
        A.require(p.returncode==1 and not p.stdout and p.stderr.startswith(b'INVALID:'),'negative CLI/authentication')
        negatives.append(dict(argv=list(map(str,argv)),code=p.returncode,stderr=p.stderr.decode()))
    A.publish(output/'negative-cli.json',A.canonical(negatives))
    neutral = U.command([executable,'--neutral'])
    A.exact(neutral,b'PASS 36 neutral cases, six ownership-matched APIs, all-arm packet/rank oracle, four late-call cleanup checks\n','neutral API replay')
    A.publish(output/'neutral.txt',neutral)
    fixture = U.command([executable,'--neutral-fixtures'])
    verify(fixture,'0'*64,mode,neutral=True)
    A.publish(output/'fixtures.jsonl',fixture)
    dependencies.update(Path(m.__file__).resolve(strict=True) for m in list(sys.modules.values())
                        if getattr(m,'__file__',None) and Path(m.__file__).is_file())
    U.freeze_inputs(dependencies,frozen)
    manifest = dict(protocol=PROTOCOL, mode=mode, commands=commands,
                    inputs=[frozen[p] for p in sorted(dependencies)],
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
    reference = None
    for mode in MODES:
        _,records = verify(A.read_regular(Path(receipt['executables'][mode]).parent/'fixtures.jsonl',RAW_CAP),'0'*64,mode,True)
        if reference is None: reference = records
        else: A.exact(records,reference,'complete neutral backend parity')


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


def run(path):
    begin = time.monotonic()
    frozen = A.read_regular(path, 1024*1024); claim = A.decode(frozen)
    A.exact(frozen, A.canonical(claim), 'canonical receipt'); current(claim)
    OUTPUT.mkdir(mode=0o700); A.publish(OUTPUT/'CLAIM.json', frozen)
    analysis, reference = dict(outcome='INVALID'), None
    try:
        for mode in MODES:
            raw, error, code, failure = capture(claim['executables'][mode], A.sha(frozen),
                min(begin+600, time.monotonic()+180), [OUTPUT/(mode+'.raw.jsonl'), OUTPUT/(mode+'.stderr.txt')])
            A.require(failure is None and code == 0 and error == b'', mode+' worker/observer: '+str(failure))
            summary, records = verify(raw, A.sha(frozen), mode)
            if reference is None:
                reference = records; analysis = summary
            else:
                A.exact(records, reference, 'cross-backend complete byte/status parity')
                A.exact(summary, analysis, 'cross-backend summary parity')
        retained_parity(reference)
        current(claim)
        A.require(time.monotonic()-begin < 600, 'whole replay deadline')
    except Exception as error:
        analysis = dict(outcome='INVALID', failure=str(error))
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'analysis.json', A.canonical(analysis))
    members = [U.pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members) < 256*1024**2-65536, 'bundle cap')
    A.publish(OUTPUT/'COMPLETE.json', A.canonical(dict(protocol=PROTOCOL, outcome=analysis['outcome'], files=members)))
    print(json.dumps({k: v for k, v in analysis.items() if k not in ('groups', 'cells')}, sort_keys=True))


def replay():
    complete = A.decode(A.read_regular(OUTPUT/'COMPLETE.json',1024*1024))
    A.exact(set(complete),{'protocol','outcome','files'},'sealed complete schema')
    A.exact(complete['protocol'],PROTOCOL,'sealed protocol')
    names = {'CLAIM.json','analysis.json'} | {m+s for m in MODES for s in ('.raw.jsonl','.stderr.txt')}
    actual_names = {p.name for p in OUTPUT.iterdir()}-{'COMPLETE.json'}
    if complete['outcome']!='INVALID': A.exact(actual_names,names,'whole sealed backend roster')
    else: A.require({'CLAIM.json','analysis.json'}<=actual_names<=names,'bounded partial roster')
    A.exact(len(complete['files']),len(actual_names),'unique sealed members')
    A.exact({p['path'] for p in complete['files']},{str(OUTPUT/n) for n in actual_names},'exact sealed member paths')
    for p in complete['files']: A.exact(U.pin(p['path']),p,'sealed member identity')
    frozen = A.read_regular(OUTPUT/'CLAIM.json',1024*1024)
    claimed = A.decode(frozen); A.exact(frozen,A.canonical(claimed),'canonical sealed receipt'); current(claimed)
    stored = A.decode(A.read_regular(OUTPUT/'analysis.json',16*1024*1024))
    A.exact(stored['protocol'],PROTOCOL,'analysis protocol'); A.exact(stored['outcome'],complete['outcome'],'sealed outcome')
    if complete['outcome']!='INVALID':
        reference = None
        for mode in MODES:
            A.exact(A.read_regular(OUTPUT/(mode+'.stderr.txt'),ERR_CAP),b'','backend stderr')
            actual,records = verify(A.read_regular(OUTPUT/(mode+'.raw.jsonl'),RAW_CAP),A.sha(frozen),mode)
            A.exact({k:v for k,v in stored.items() if k not in ('protocol','elapsed_seconds')},actual,'complete backend replay')
            if reference is None: reference = records
            else: A.exact(records,reference,'complete backend byte/status parity')
        retained_parity(reference)
    return stored


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('mode', choices=MODES); b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    sub.add_parser('replay')
    args = parser.parse_args()
    if args.command == 'build': build(args.mode, args.output)
    elif args.command == 'receipt': A.publish(args.output, A.canonical(receipt(args.build_dir.resolve(strict=True))))
    elif args.command == 'run': run(args.receipt)
    else: print(json.dumps({k:v for k,v in replay().items() if k not in ('groups','cells')},sort_keys=True))
