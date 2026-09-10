#!/usr/bin/env python3
"""Frozen K8 serialized recovery; six ownership-matched arms, retained cohort only."""
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


U = sibling('_k8_recovery_build_io', 'Wh2K8RecoveryBuildR0.py')
D = sibling('_k8_recovery_retained_data', 'Wh2K8NativeDataR0.py')
A = U.A
ROOT = U.ROOT
PROTOCOL = 'wirehair.wh2.k8-serialized-recovery-r0'
OUTPUT = Path('/var/tmp/wh2-k8-serialized-recovery-r0')
MODES = ('native', 'scalar', 'asan')
WIDTHS = (2, 64, 1280)
ARMS = ('ordinary_wh2_independent','k8_independent','wh1_owned',
        'ordinary_wh2_borrowed','k8_borrowed','wh1_borrowed')
K = 8
RAW_CAP, ERR_CAP = 96*1024**2, 65536
PAIRED = ((1,0),(1,2),(4,3),(4,5))
ENV_KEYS = U.ENV_KEYS + ('ASAN_OPTIONS', 'UBSAN_OPTIONS')
SANITIZERS = dict(ASAN_OPTIONS='detect_leaks=1:detect_stack_use_after_return=1:halt_on_error=1',
                  UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
SOURCES = ('bench/Wh2K8SerializedRecoveryR0.cpp', 'bench/Wh2SmallRecoveryWorkerR0.h',
           'bench/Wh2K8SerializedRecoveryR0.py', 'bench/test_Wh2K8SerializedRecoveryR0.py',
           'bench/Wh2K8RecoveryBuildR0.py', 'bench/test_Wh2K8RecoveryBuildR0.py',
           'bench/Wh2K8SerializedRecoveryR0.md', 'bench/Wh2K8NativeDataR0.py',
           'bench/Wh2K3NativeDataR0.py', 'bench/Wh2NoncommutingRadixRunR0.py',
           'bench/Wh2AlignedIntermediateCostR0.py', 'V2_WIRE_PROFILE.md',
           'WH2_BORROWED_SOURCE_API.md')


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
        streams = [list(range(12)),list(range(8,20)),[0xffffffff-2*j for j in range(12)],[0]*12]
        for width in WIDTHS:
            for shape, tail in enumerate((width,1)):
                for index, ids in enumerate(streams):
                    result.append(dict(group=3,index=shape*4+index,B=width,tail=tail,ids=ids))
        for index in list(range(12))+list(range(6132,6144)):
            row = report['fresh'][index]
            result.append(dict(group=4,index=index,B=row['B'],tail=row['B'],ids=row['ids']))
        A.exact(len(result),48,'bounded first/last-cell smoke roster')
    else:
        for group, name in enumerate(('fresh','hard')):
            for index, row in enumerate(report[name]):
                result.append(dict(group=group,index=index,B=row['B'],tail=row['B'],ids=row['ids'],
                                   schedule=row['schedule'],root=row['root']))
        for index, row in enumerate(data['history']):
            result.append(dict(group=2,index=index,B=row['B'],tail=row['tail'],ids=row['ids'],
                               origin=report['inputs']['origins'][index]))
        A.exact(len(result),6260,'retained trace/origin roster')
        A.exact(sum(len(e['ids']) for e in result),74946,'exact original packet horizons')
    return result


@lru_cache(maxsize=256)
def products(coefficient):
    return bytes(D.R.multiply(coefficient,b) for b in range(256))


@lru_cache(maxsize=1)
def powers():
    def product(a,b):
        result = [0]*(K*K)
        for i in range(K):
            for j in range(K):
                for k in range(K): result[K*i+j] ^= products(a[K*i+k])[b[K*k+j]]
        return tuple(result)
    matrices = []
    for feedback in ((96,19,186,153,85,252,7,255),(98,19,186,153,85,252,7,255)):
        m = [0]*(K*K)
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
    row = (1,)+(0,)*(K-1)
    for bit in range(32):
        if packet_id & (1 << bit):
            matrix = powers()[bit][bin(packet_id >> (bit+1)).count('1') % 2]
            next_row = [0]*K
            for i in range(K):
                for j in range(K): next_row[i] ^= products(matrix[K*i+j])[row[j]]
            row = tuple(next_row)
    return row


@lru_cache(maxsize=6)
def message(width,tail):
    A.require(type(width) is int and width in WIDTHS,'frozen width')
    A.integer(tail,1,width)
    return bytes((37*i+i//11) % 256 for i in range((K-1)*width+tail))


@lru_cache(maxsize=32768)
def packet_hash(width,tail,row,length):
    A.require(len(row)==K and length in (width,tail),'packet shape')
    source = message(width,tail)+bytes(width-tail)
    packet = bytearray(width)
    for k,value in enumerate(row):
        table = products(value)
        for j in range(width): packet[j] ^= table[source[k*width+j]]
    return A.sha(packet[:length])


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
                    key = [e['index'],e['B'],e['tail']]
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
    superior = all(len(pair['fixed'][0])>len(pair['introduced'][0]) for pair in groups[0]['paired'])
    target = (all(100*fresh[a][0]<=6144 for a in (1,4)) and
              all(100*cell['failure_counts_oh0_to_4'][a][0]<=512 for cell in cells for a in (1,4)) and
              all(groups[1]['failure_counts_oh0_to_4'][a][0]==0 and
                  groups[2]['unresolved'][a]==0 for a in (1,4)))
    return dict(outcome='PASS' if target and superior else 'FAIL',execution_valid=True,
                retained_target_qualified=target,retained_paired_oh0_superior=superior,
                arms=list(ARMS),groups=groups,cells=cells,
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
                       sources=[A.sha(message(b,b)) for b in WIDTHS]),'exact replay header')
    profiles = {}
    for row,e in zip(records,expected):
        A.exact(set(row),{'type','group','index','B','tail','ids','arms'},'record schema')
        A.exact({k:row[k] for k in ('group','index','B','tail','ids')},
                {k:e[k] for k in ('group','index','B','tail','ids')},'retained exact chronology')
        A.exact(row['type'],'record','record type'); A.exact(len(row['arms']),6,'six ownership-matched API routes')
        width,tail,ids = e['B'],e['tail'],e['ids']
        size = len(message(width,tail))
        for arm,r in enumerate(row['arms']):
            A.exact(set(r),{'profile','packets','rows','feed','first','recoveries','counts','checked'},'arm schema')
            valid_hex(r['profile'],32)
            A.exact(r['profile'],profiles.setdefault((width,tail,arm),r['profile']),'stable full descriptor')
            profile = bytes.fromhex(r['profile'])
            if arm in (0,3):
                A.exact(profile[:16],struct.pack('<4sHHQ',b'WHV2',1,32,0x4b295bbb47f4f9c9),'ordinary certified WH2 descriptor')
                A.exact(profile[16:28],struct.pack('<QI',size,width),'certified WH2 dimensions')
                A.exact(profile[29:],bytes(3),'certified reserved bytes')
            elif arm in (1,4):
                A.exact(profile,struct.pack('<4sHHQQII',b'WHK8',1,32,0x5748324b38544d31,size,width,0),'sealed WHK8 identity')
            else: A.exact(profile,bytes(32),'WH1 has no descriptor')
            valid_hex(r['rows'],K*len(ids))
            encoded_rows = bytes.fromhex(r['rows'])
            observed = tuple(tuple(encoded_rows[i:i+K]) for i in range(0,len(encoded_rows),K))
            if arm in (1,4):
                A.exact(observed,tuple(coefficient(i) for i in ids),'independent sealed K8 rows')
            A.exact(len(r['packets']),len(ids),'all retained packets')
            for packet in r['packets']: valid_hex(packet,32)
            A.exact(r['packets'],[packet_hash(width,tail,rr,tail if packet_id==K-1 else width)
                                       for packet_id,rr in zip(ids,observed)],'independent every-arm packet hashes')
            first = A.integer(r['first'],0,len(ids))
            A.require(first==0 or first>=K,'possible first success')
            A.exact(first,first_success(observed),'independent every-arm first-success rank')
            A.exact(r['feed'],[1]*(first-1)+[0] if first else [1]*len(ids),'all prefix feed statuses')
            A.exact(r['recoveries'],2 if first else 0,'two byte-exact recoveries')
            A.exact(r['counts'],[9,9*len(ids),1,first if first else len(ids),2 if first else 0,10],'whole attempted API ledger')
            A.exact(r['checked'],True,'payload/guard validation')
            if arm in (1,4) and not neutral:
                if e['group']!=2:
                    ranks = retained()[0]['fresh' if e['group']==0 else 'hard'][e['index']]['ranks']
                    A.exact(first,next((K+j for j,v in enumerate(ranks) if v==K),0),'sealed retained prefix outcome')
        for a in range(3): A.exact(row['arms'][a],row['arms'][a+3],'policy-independent complete result')
    A.exact(footer,dict(type='footer',records=len(expected),checked=True),'complete footer')
    return (dict(outcome='PASS',neutral_cases=48) if neutral else summarize(records,expected)),records


def build(mode,output):
    return U.build(mode,output)


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
        manifest = A.decode(A.read_regular(path, 4*1024*1024))
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
        manifest = A.decode(A.read_regular(path, 4*1024*1024))
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
    analysis, reference = dict(outcome='INVALID',execution_valid=False), None
    try:
        for mode in MODES:
            raw, error, code, failure = capture(claim['executables'][mode], A.sha(frozen),
                min(begin+900, time.monotonic()+245), [OUTPUT/(mode+'.raw.jsonl'), OUTPUT/(mode+'.stderr.txt')])
            A.require(failure is None and code == 0 and error == b'', mode+' worker/observer: '+str(failure))
            summary, records = verify(raw, A.sha(frozen), mode)
            if reference is None:
                reference = records; analysis = summary
            else:
                A.exact(records, reference, 'cross-backend complete byte/status parity')
                A.exact(summary, analysis, 'cross-backend summary parity')
        current(claim)
        A.require(time.monotonic()-begin < 900, 'whole replay deadline')
    except Exception as error:
        analysis = dict(outcome='INVALID',execution_valid=False, failure=str(error))
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'analysis.json', A.canonical(analysis))
    members = [U.pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members) < 384*1024**2-65536, 'bundle cap')
    A.publish(OUTPUT/'COMPLETE.json', A.canonical(dict(protocol=PROTOCOL, outcome=analysis['outcome'], files=members)))
    print(json.dumps({k: v for k, v in analysis.items() if k not in ('groups', 'cells')}, sort_keys=True))


def replay():
    complete = A.decode(A.read_regular(OUTPUT/'COMPLETE.json',4*1024*1024))
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
