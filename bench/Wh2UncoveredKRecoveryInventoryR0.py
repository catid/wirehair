#!/usr/bin/env python3
"""Frozen baseline-only recovery inventory. No candidate, timings or promotion.

Native generator rows are observed with basis messages, not inferred from the
decoder. Independent GF(256) arithmetic checks packet linearity and rank.
Actual decoder success, never rank alone, defines the reported failure counts.
"""
import argparse
import ctypes as T
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path
import resource
import struct
import subprocess
import sys
import time

SPEC = importlib.util.spec_from_file_location('inventory_provenance',
    Path(__file__).with_name('Wh2SmallIsolationPreservedCostR0.py'))
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)
R, A, N, ROOT = P.R, P.A, P.R.N, P.ROOT
PROTOCOL = 'wirehair.wh2.uncovered-k-recovery-inventory-r0'
OUTPUT = Path('/var/tmp/wh2-uncovered-k-recovery-inventory-r0')
LIBRARY = P.LIBRARIES[1]
KS, WIDTHS = (2, 4, 5, 8), (2, 64, 1280)
SCHEDULES = ('iid', 'burst', 'adversarial', 'repair-only')
MASK = (1 << 64) - 1
RAW_CAP, ERR_CAP = 128*1024**2, 65536
ENV_KEYS = R.ENV_KEYS + ('ASAN_OPTIONS', 'UBSAN_OPTIONS')
SOURCES = ('bench/Wh2UncoveredKRecoveryInventoryR0.py',
           'bench/test_Wh2UncoveredKRecoveryInventoryR0.py',
           'bench/Wh2SmallIsolationPreservedCostR0.py',
           'bench/Wh2AdmissionRegressionCostR0.py',
           'bench/Wh2AdmissionRegressionNeutral.py',
           'bench/Wh2K3OrdinaryCostR0.py',
           'bench/Wh2AlignedIntermediateCostR0.py')


def roots():
    return ['0x'+A.sha((PROTOCOL+':inventory/'+str(i)).encode())[:16] for i in range(64)]


def trace(k, b, root, schedule):
    A.require(k in KS and b in WIDTHS and schedule in SCHEDULES, 'trace dimensions')
    state = (int(root, 16) ^ k*0x9e3779b97f4a7c15 ^ b*0xbf58476d1ce4e5b9) & MASK
    if schedule != 'iid': state ^= 0x10fade
    def uniform():
        nonlocal state
        state = (state+0x9e3779b97f4a7c15) & MASK
        x = ((state ^ (state >> 30))*0xbf58476d1ce4e5b9) & MASK
        x = ((x ^ (x >> 27))*0x94d049bb133111eb) & MASK
        return ((x ^ (x >> 31)) >> 11)*2.0**-53
    ids, burst = [], 0
    loss = .1 if schedule == 'iid' else .5
    for candidate in range(65536):
        if schedule == 'burst' and burst:
            burst -= 1
            continue
        if uniform() < (loss/(8-7*loss) if schedule == 'burst' else loss):
            if schedule == 'burst': burst = 7
            continue
        ids.append(0xffffffff-2*candidate if schedule == 'adversarial' else
                   k+candidate if schedule == 'repair-only' else candidate)
        if len(ids) == k+4: return ids
    raise ValueError('trace candidate bound')


@lru_cache(maxsize=1)
def roster():
    rows = []
    for k in KS:
        for b in WIDTHS:
            for schedule in SCHEDULES:
                for root in roots():
                    rows.append(dict(group='inventory', k=k, b=b, tail=b,
                        schedule=schedule, root=root, ids=trace(k,b,root,schedule)))
    for k in KS:
        for b in WIDTHS:
            for tail in (1,b):
                rows.append(dict(group='hard', k=k, b=b, tail=tail,
                    schedule='low-repair', root=None, ids=list(range(k,2*k+4))))
    A.exact(len(rows),3096,'complete prospective roster')
    return rows


def dimensions(row): return row['k'], row['b'], row['tail']


def source(k,b,tail):
    return bytes((37*i+i//11)%256 for i in range((k-1)*b+tail))


def multiply(x,y):
    value = 0
    for bit in range(8):
        if (y >> bit) & 1: value ^= x << bit
    for bit in range(14,7,-1):
        if (value >> bit) & 1: value ^= 0x14d << (bit-8)
    return value


@lru_cache(maxsize=256)
def products(x): return bytes(multiply(x,y) for y in range(256))


def packet(row, data, b, packet_id):
    value = 0
    for i,c in enumerate(row):
        part = data[i*b:(i+1)*b].ljust(b,b'\0')
        value ^= int.from_bytes(part.translate(products(c)), 'little')
    out = value.to_bytes(b,'little')
    return out[:len(data)-(len(row)-1)*b] if packet_id == len(row)-1 else out


def rank(rows, k):
    basis = {}
    for row in rows:
        v = list(row)
        A.require(len(v)==k and all(type(c) is int and 0<=c<=255 for c in v),'GF256 row')
        for i in range(k):
            if not v[i]: continue
            if i in basis:
                table = products(v[i])
                v = [a ^ table[c] for a,c in zip(v,basis[i])]
            else:
                table = products(products(v[i]).index(1))
                basis[i] = [table[c] for c in v]
                break
    return len(basis)


def first_rank(rows, k):
    return next((n for n in range(k,len(rows)+1) if rank(rows[:n],k)==k),0)


def guard(data):
    return (T.c_ubyte*(len(data)+2)).from_buffer_copy(b'\xa5'+data+b'\xa5')


class Api:
    def __init__(self, library, arm):
        A.require(arm in ('wh2','wh1'),'actual baseline arm')
        self.arm = arm
        self.calls = [0]*6  # encoder create, encode, decoder create, feed, recover, free
        prefix = 'wirehair_v2_' if arm=='wh2' else 'wirehair_'
        call, V, U32, U64, INT = library.call, N.V, N.U32, N.U64, N.INT
        self.free = call(prefix+'free',None,V)
        self.encode = call(prefix+'encode',INT,V,U32,V,U32,V)
        self.decode = call(prefix+'decode',INT,V,U32,V,U32)
        if arm=='wh2':
            self.create = call(prefix+'encoder_create',INT,V,U64,U32,V,U32,V,V)
            self.receiver = call(prefix+'decoder_create',INT,V,U32,V)
            self.recover = call(prefix+'recover',INT,V,V,U64,V)
        else:
            self.create = call(prefix+'encoder_create_owned_ex',INT,V,V,U64,U32,V)
            self.receiver = call(prefix+'decoder_create_ex',INT,V,U64,U32,V)
            self.recover = call(prefix+'recover',INT,V,V,U64)

    def packets(self,data,b,ids):
        encoder, written = N.V(), N.U32()
        incoming = T.create_string_buffer(data,len(data))
        profile = (T.c_ubyte*32)()
        try:
            self.calls[0] += 1
            status = (self.create(incoming,len(data),b,profile,32,T.byref(written),T.byref(encoder))
                      if self.arm=='wh2' else self.create(None,incoming,len(data),b,T.byref(encoder)))
            A.exact(status,0,'owned ordinary encoder creation')
            A.require(encoder.value is not None,'live encoder')
            if self.arm=='wh2':
                A.exact(written.value,32,'descriptor length')
                check_profile(bytes(profile).hex(),len(data),b)
            # Both constructors promise independence from caller-owned input.
            T.memset(incoming,0xcc,len(data))
            result = []
            k = (len(data)+b-1)//b
            for i in ids:
                out = guard(bytes([0xa5])*b); written.value = 0xffffffff
                self.calls[1] += 1
                A.exact(self.encode(encoder,i,T.byref(out,1),b,T.byref(written)),0,'native encode')
                size = len(data)-(k-1)*b if i==k-1 else b
                A.exact(written.value,size,'exact packet length')
                A.exact(bytes(out)[:1]+bytes(out)[1+size:],bytes([0xa5])*(b+2-size),'encode guards')
                result.append(bytes(out)[1:1+size])
            A.exact(incoming.raw,bytes([0xcc])*len(data),'owned source remains untouched')
            return bytes(profile).hex(),result
        finally:
            if encoder.value is not None:
                self.calls[5] += 1; self.free(encoder)

    def receive(self,profile,data,b,ids,packets):
        decoder = N.V()
        encoded = T.create_string_buffer(bytes.fromhex(profile),32)
        try:
            self.calls[2] += 1
            status = (self.receiver(encoded,32,T.byref(decoder)) if self.arm=='wh2' else
                      self.receiver(None,len(data),b,T.byref(decoder)))
            A.exact(status,0,'ordinary decoder creation')
            A.require(decoder.value is not None,'live decoder')
            statuses = []
            for i,payload in zip(ids,packets):
                incoming = guard(payload); before = bytes(incoming)
                self.calls[3] += 1
                status = self.decode(decoder,i,T.byref(incoming,1),len(payload))
                A.require(status in (0,1),'decode status')
                A.exact(bytes(incoming),before,'immutable packet and guards')
                statuses.append(status)
                if status==0: break
            first = len(statuses) if statuses[-1]==0 else 0
            recovered = []
            if first:
                for _ in range(2):
                    out = guard(bytes([0xa5])*len(data)); written = N.U64(0xffffffffffffffff)
                    self.calls[4] += 1
                    status = (self.recover(decoder,T.byref(out,1),len(data),T.byref(written))
                              if self.arm=='wh2' else self.recover(decoder,T.byref(out,1),len(data)))
                    A.exact(status,0,'recover success')
                    if self.arm=='wh2': A.exact(written.value,len(data),'recover length')
                    A.exact(bytes(out),b'\xa5'+data+b'\xa5','exact recovered bytes and guards')
                    recovered.append(A.sha(data))
            return dict(feed=statuses,first=first,recovered=recovered)
        finally:
            if decoder.value is not None:
                self.calls[5] += 1; self.free(decoder)


def check_profile(profile,m,b):
    raw = bytes.fromhex(profile)
    A.exact(len(raw),32,'descriptor bytes')
    A.exact(raw[:28],struct.pack('<4sHHQQI',b'WHV2',1,32,0x4b295bbb47f4f9c9,m,b),
            'ordinary uncovered K selects existing GF256 equations')
    A.exact(raw[29:],bytes(3),'descriptor reserved bytes')


def id_sets(cases):
    result = {}
    for r in cases: result.setdefault(dimensions(r),set()).update(r['ids'])
    return result


def observe_rows(apis,cases):
    result = []
    for (k,b,tail),ids in sorted(id_sets(cases).items()):
        ids = sorted(ids); arms = []
        for api in apis:
            columns, profiles = [], []
            for column in range(k):
                data = bytes(int(i//b==column) for i in range((k-1)*b+tail))
                profile, packets = api.packets(data,b,ids)
                profiles.append(profile)
                coefficients = [p[0] for p in packets]
                for i,c,payload in zip(ids,coefficients,packets):
                    row = [0]*k; row[column] = c
                    A.exact(payload,packet(row,data,b,i),'native basis payload across all lanes')
                columns.append(coefficients)
            A.exact(len(set(profiles)),1,'content-independent construction descriptor')
            arms.append(dict(profile=profiles[0],rows=[list(r) for r in zip(*columns)]))
        result.append(dict(k=k,b=b,tail=tail,ids=ids,arms=arms))
    return result


def row_index(observed,cases):
    expected = id_sets(cases)
    A.exact(len(observed),len(expected),'every observed dimension')
    result = {}
    for r,key in zip(observed,sorted(expected)):
        A.exact(set(r),{'k','b','tail','ids','arms'},'basis observation schema')
        A.exact(dimensions(r),key,'basis dimension chronology')
        A.exact(r['ids'],sorted(expected[key]),'complete unique packet-ID roster')
        A.exact(len(r['arms']),2,'both native generators')
        k,b,tail = key
        for arm,x in enumerate(r['arms']):
            A.exact(set(x),{'profile','rows'},'basis arm schema')
            if arm==0: check_profile(x['profile'],(k-1)*b+tail,b)
            else: A.exact(x['profile'],'00'*32,'WH1 descriptor placeholder')
            A.exact(len(x['rows']),len(r['ids']),'every native row')
            for packet_id,row in zip(r['ids'],x['rows']):
                A.require(len(row)==k and all(type(v) is int and 0<=v<256 for v in row),'row bytes')
                if packet_id<k: A.exact(row,[int(i==packet_id) for i in range(k)],'systematic equation')
            result[(key,arm)] = x['profile'],dict(zip(r['ids'],x['rows']))
    return result


def evaluate(apis,cases,observed,emit):
    index = row_index(observed,cases)
    for ordinal,case in enumerate(cases):
        k,b,tail = dimensions(case); data = source(k,b,tail); arms = []
        for arm,api in enumerate(apis):
            profile,packets = api.packets(data,b,case['ids'])
            expected,rows = index[(dimensions(case),arm)]
            A.exact(profile,expected,'same native profile for basis and message')
            for i,payload in zip(case['ids'],packets):
                A.exact(payload,packet(rows[i],data,b,i),'independent GF256 packet oracle')
            # packets() has already freed the encoder before receiver creation.
            r = api.receive(profile,data,b,case['ids'],packets)
            r.update(profile=profile,packets=[p.hex() for p in packets],
                     rank_first=first_rank([rows[i] for i in case['ids']],k))
            if r['first']: A.require(0<r['rank_first']<=r['first'],'success requires full rank')
            arms.append(r)
        emit(dict(type='case',ordinal=ordinal,case=case,arms=arms))


def summarize(records):
    def counts(selected,arm,field):
        return [sum(r['arms'][arm][field]==0 or r['arms'][arm][field]>r['case']['k']+oh
                    for r in selected) for oh in range(5)]
    cells,totals = [],[]
    for k in KS:
        selected = [r for r in records if r['case']['group']=='inventory' and r['case']['k']==k]
        A.exact(len(selected),768,'complete per-K inventory denominator')
        actual = [counts(selected,a,'first') for a in (0,1)]
        totals.append(dict(k=k,traces=768,failures=actual,rank_failures=[counts(selected,a,'rank_first') for a in (0,1)],
            worse=[r['ordinal'] for r in selected if (r['arms'][0]['first']==0 or r['arms'][0]['first']>k) and r['arms'][1]['first']==k],
            better=[r['ordinal'] for r in selected if (r['arms'][1]['first']==0 or r['arms'][1]['first']>k) and r['arms'][0]['first']==k]))
        for b in WIDTHS:
            for schedule in SCHEDULES:
                cell = [r for r in selected if r['case']['b']==b and r['case']['schedule']==schedule]
                A.exact(len(cell),64,'complete cell denominator')
                cells.append(dict(k=k,b=b,schedule=schedule,traces=64,failures=[counts(cell,a,'first') for a in (0,1)]))
    ordered = sorted(totals,key=lambda r:(r['failures'][0][0]-r['failures'][1][0],r['failures'][0][4],r['failures'][0][0],-r['k']),reverse=True)
    top = ordered[0]
    recommend = top['k'] if top['failures'][0][0]>top['failures'][1][0] and top['failures'][0][0]*100>top['traces'] else None
    hard = [r for r in records if r['case']['group']=='hard']
    A.exact(len(hard),24,'separate full/partial hard cases')
    return dict(protocol=PROTOCOL,outcome='DIAGNOSTIC_COMPLETE',arms=['ordinary_wh2','wh1'],totals=totals,cells=cells,
        hard_failures=[counts(hard,a,'first') for a in (0,1)],priority_order=[r['k'] for r in ordered],recommended_k=recommend,
        decoder_lag=[sum(r['arms'][a]['rank_first']>0 and (r['arms'][a]['first']==0 or r['arms'][a]['first']>r['arms'][a]['rank_first']) for r in records) for a in (0,1)],
        candidate_tested=False,speed_claimed=False,holdout=False,all_K_claimed=False,promotion_claimed=False)


def check_library(report):
    A.exact(set(report),{'path','sha256','base','exports','slots','providers','context','context_bytes','features'},'library report schema')
    A.exact((report['path'],report['sha256']),(str(LIBRARY[0]),LIBRARY[1]),'native library identity')
    raw = A.read_regular(LIBRARY[0],4*1024**2)
    A.exact(A.sha(raw),LIBRARY[1],'actual ELF bytes')
    elf = N.Elf(raw); base = report['base']
    A.require(type(base) is int and base>0,'loaded base')
    exports = {name:dict(address=base+elf.symbol(name,2)[0],offset=elf.symbol(name,2)[0],size=elf.symbol(name,2)[1]) for name in elf.exports}
    A.exact(report['exports'],exports,'all owned exports')
    A.exact(report['slots'],[dict(name=name,offset=offset,target=exports[name]['address']) for name,offset in elf.slots],'own internal calls')
    context,size = elf.symbol('GF256Ctx',1)
    A.exact((report['context'],report['context_bytes']),(base+context,size),'native GF context')
    A.exact(report['features'],[1]*4,'native GF features')
    A.exact(set(report['providers']),set(N.PROVIDERS),'runtime provider roster')
    A.require(all(type(x) is int and x>0 for x in report['providers'].values()),'resolved providers')


def expected_calls(records,cases):
    ids = id_sets(cases)
    creates = sum(key[0] for key in ids)+len(cases)
    encodes = sum(key[0]*len(values) for key,values in ids.items())+sum(len(r['ids']) for r in cases)
    return [[creates,encodes,len(cases),sum(len(r['arms'][a]['feed']) for r in records),
             sum(len(r['arms'][a]['recovered']) for r in records),creates+len(cases)] for a in (0,1)]


def verify(path,claim):
    records = []; cases = roster()
    A.require(path.stat().st_size<=RAW_CAP,'raw bound')
    with path.open('rb') as stream:
        def next_record():
            line = stream.readline(8*1024**2+1)
            A.require(0<len(line)<=8*1024**2,'bounded complete JSONL record')
            value = A.decode(line)
            A.exact(line,A.canonical(value),'canonical JSONL record')
            return value
        header = next_record()
        A.exact(set(header),{'type','protocol','claim','library','coefficients'},'header schema')
        A.exact((header['type'],header['protocol'],header['claim']),('header',PROTOCOL,claim),'header identity')
        check_library(header['library'])
        index = row_index(header['coefficients'],cases)
        for ordinal,case in enumerate(cases):
            r = next_record()
            A.exact(set(r),{'type','ordinal','case','arms'},'case schema')
            A.exact((r['type'],r['ordinal'],r['case']),('case',ordinal,case),'whole case chronology')
            A.exact(len(r['arms']),2,'paired actual decoders')
            k,b,tail = dimensions(case); data = source(k,b,tail)
            for arm,x in enumerate(r['arms']):
                A.exact(set(x),{'profile','packets','feed','first','recovered','rank_first'},'arm schema')
                profile,rows = index[(dimensions(case),arm)]
                A.exact(x['profile'],profile,'actual descriptor')
                A.exact(x['packets'],[packet(rows[i],data,b,i).hex() for i in case['ids']],'every payload byte')
                A.require(type(x['first']) is int and (x['first']==0 or k<=x['first']<=k+4),'first success range')
                A.exact(x['feed'],[1]*(x['first']-1)+[0] if x['first'] else [1]*(k+4),'actual first success trace')
                A.exact(x['recovered'],[A.sha(data)]*2 if x['first'] else [],'two exact recoveries')
                A.exact(x['rank_first'],first_rank([rows[i] for i in case['ids']],k),'independent row rank')
                if x['first']: A.require(0<x['rank_first']<=x['first'],'successful decode has full information')
                del x['packets']
            records.append(r)
        A.exact(next_record(),dict(type='footer',complete=True,records=len(cases),calls=expected_calls(records,cases)),'complete footer and API ledger')
        A.exact(stream.read(),b'','no trailing records')
    return summarize(records)


def environment():
    result = {k:os.environ.get(k) for k in ENV_KEYS}
    A.exact(result,{k:None for k in ENV_KEYS},'clean native environment')
    return result


def receipt(neutral):
    A.require(neutral.is_absolute(),'absolute external neutral qualification')
    A.exact(A.decode(A.read_regular(neutral,65536)),neutral_expected(),'neutral qualification')
    proof,inputs = P.provenance()
    inputs.update(ROOT/p for p in SOURCES)
    inputs.update((Path(__file__).resolve(),Path(sys.executable).resolve(strict=True),neutral))
    for target in (LIBRARY[0],Path(sys.executable)):
        inputs.update(R.runtime_dependencies(target))
    frozen = dict(protocol=PROTOCOL,head=R.command(['git','rev-parse','HEAD']).decode().strip(),
        source_head=P.SOURCE_HEAD,interpreter=str(Path(sys.executable).resolve(strict=True)),neutral=str(neutral),environment=environment(),
        library_provenance=proof,pins=[R.O.pin(p) for p in sorted(inputs)])
    current(frozen)
    return frozen


def current(frozen):
    A.exact(set(frozen),{'protocol','head','source_head','interpreter','neutral','environment','library_provenance','pins'},'receipt schema')
    A.exact((frozen['protocol'],frozen['source_head']),(PROTOCOL,P.SOURCE_HEAD),'receipt identities')
    A.exact(frozen['environment'],environment(),'receipt environment')
    A.exact(R.command(['git','rev-parse','HEAD']).decode().strip(),frozen['head'],'exact source HEAD')
    paths = {p['path'] for p in frozen['pins']}
    A.exact(len(paths),len(frozen['pins']),'unique receipt paths')
    for p in frozen['pins']:
        path = Path(p['path']); A.exact(R.O.pin(path),p,'current input')
        if ROOT in path.parents:
            data = R.command(['git','cat-file','blob',frozen['head']+':'+str(path.relative_to(ROOT))])
            A.exact((len(data),A.sha(data)),(p['bytes'],p['sha256']),'committed source')
    A.require(all(str(ROOT/p) in paths for p in SOURCES),'mandatory sources')
    A.require(str(LIBRARY[0]) in paths and frozen['interpreter'] in paths,'mandatory runtime')
    neutral = Path(frozen['neutral'])
    A.require(neutral.is_absolute() and ROOT not in neutral.parents and str(neutral) in paths,'mandatory external neutral qualification')
    A.exact(A.decode(A.read_regular(neutral,65536)),neutral_expected(),'positive neutral qualification')
    for target in (LIBRARY[0],Path(frozen['interpreter'])):
        A.require(all(str(p) in paths for p in R.runtime_dependencies(target)),'complete loader dependencies')
    proof,inputs = P.provenance()
    A.exact(proof,frozen['library_provenance'],'exact current-library producing chain')
    A.require(all(str(p) in paths for p in inputs),'complete producer inputs')


def emit(row):
    sys.stdout.buffer.write(A.canonical(row)); sys.stdout.buffer.flush()


def worker(claim):
    A.require(len(claim)==64 and all(c in '0123456789abcdef' for c in claim),'claim hex')
    raw = A.read_regular(OUTPUT/'CLAIM.json',1024*1024)
    A.exact(A.sha(raw),claim,'actual inventory claim')
    current(A.decode(raw))
    resource.setrlimit(resource.RLIMIT_CPU,(120,120))
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    resource.setrlimit(resource.RLIMIT_FSIZE,(RAW_CAP,RAW_CAP))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    library = N.Library(0,(LIBRARY,)); apis = [Api(library,a) for a in ('wh2','wh1')]
    observed = observe_rows(apis,roster())
    emit(dict(type='header',protocol=PROTOCOL,claim=claim,library=library.report,coefficients=observed))
    evaluate(apis,roster(),observed,emit)
    library.check_bindings()
    emit(dict(type='footer',complete=True,records=len(roster()),calls=[api.calls for api in apis]))


def neutral_expected():
    return dict(protocol=PROTOCOL,library_sha256=LIBRARY[1],cases=24,arms=2,checked=True,scientific_launch=False)


def neutral(path):
    A.require(path.is_absolute() and ROOT not in path.parents and path.parent.exists(),'external neutral report')
    A.require(not path.exists() and not path.is_symlink(),'new neutral report')
    environment(); library = N.Library(0,(LIBRARY,))
    apis = [Api(library,a) for a in ('wh2','wh1')]
    cases = [r for r in roster() if r['group']=='hard']; observed = observe_rows(apis,cases)
    rows = []; evaluate(apis,cases,observed,rows.append)
    A.exact(len(rows),24,'neutral hard roster')
    library.check_bindings(); A.publish(path,A.canonical(neutral_expected()))


def check_process(process):
    A.exact(set(process),{'returncode','error','elapsed_seconds'},'process schema')
    A.require(type(process['returncode']) is int and process['returncode']==0 and
              process['error'] is None,'worker completion')
    A.require(type(process['elapsed_seconds']) in (int,float) and
              0<process['elapsed_seconds']<150,'observer deadline')


def run(path):
    raw = A.read_regular(path,1024*1024); frozen = A.decode(raw)
    A.exact(raw,A.canonical(frozen),'canonical receipt'); current(frozen)
    OUTPUT.mkdir(mode=0o700); A.publish(OUTPUT/'CLAIM.json',raw)
    files = []; child = None; error = None; begin = time.monotonic()
    try:
        for name in ('raw.jsonl','stderr.txt'):
            files.append(os.open(str(OUTPUT/name),os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_NOFOLLOW,0o600))
        child = subprocess.Popen([frozen['interpreter'],str(Path(__file__).resolve()),'worker',A.sha(raw)],
            stdin=subprocess.DEVNULL,stdout=files[0],stderr=files[1],close_fds=True)
        while child.poll() is None:
            A.require(time.monotonic()-begin<150,'observer deadline')
            A.require(os.fstat(files[0]).st_size<=RAW_CAP and os.fstat(files[1]).st_size<=ERR_CAP,'output bounds')
            try: child.wait(timeout=.05)
            except subprocess.TimeoutExpired: pass
    except BaseException as problem:
        error = type(problem).__name__+': '+str(problem)
    finally:
        if child is not None:
            if child.poll() is None: child.kill()
            child.wait()
        for fd in files:
            os.fsync(fd); os.fchmod(fd,0o400); os.close(fd)
    process = dict(returncode=None if child is None else child.returncode,error=error,elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'process.json',A.canonical(process))
    try:
        # Check again after reaping: exit may race the last observer poll.
        check_process(process)
        A.exact(A.read_regular(OUTPUT/'stderr.txt',ERR_CAP),b'','empty worker stderr')
        current(frozen); analysis = verify(OUTPUT/'raw.jsonl',A.sha(raw))
    except Exception as problem:
        analysis = dict(protocol=PROTOCOL,outcome='INVALID',error=str(problem),promotion_claimed=False)
    A.publish(OUTPUT/'analysis.json',A.canonical(analysis))
    members = [R.O.pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members)<256*1024**2,'bundle bound')
    A.publish(OUTPUT/'COMPLETE.json',A.canonical(dict(protocol=PROTOCOL,outcome=analysis['outcome'],files=members)))
    print(json.dumps(dict(outcome=analysis['outcome'],recommended_k=analysis.get('recommended_k'),process=process)))


def replay():
    complete = A.decode(A.read_regular(OUTPUT/'COMPLETE.json',65536))
    A.exact(set(complete),{'protocol','outcome','files'},'terminal schema')
    A.exact(complete['protocol'],PROTOCOL,'terminal identity')
    names = {'CLAIM.json','raw.jsonl','stderr.txt','process.json','analysis.json'}
    A.exact({p.name for p in OUTPUT.iterdir()},names|{'COMPLETE.json'},'whole bundle')
    A.exact(len(complete['files']),len(names),'unique terminal members')
    A.exact({Path(p['path']).name for p in complete['files']},names,'terminal members')
    for p in complete['files']:
        path = Path(p['path']); A.exact(path.parent,OUTPUT,'member parent')
        A.exact(path.stat().st_mode&0o777,0o400,'sealed member')
        A.exact(R.O.pin(path),p,'unchanged member')
    raw = A.read_regular(OUTPUT/'CLAIM.json',1024*1024); current(A.decode(raw))
    process = A.decode(A.read_regular(OUTPUT/'process.json',65536))
    check_process(process)
    A.exact(A.read_regular(OUTPUT/'stderr.txt',ERR_CAP),b'','worker stderr')
    result = verify(OUTPUT/'raw.jsonl',A.sha(raw))
    A.exact(result,A.decode(A.read_regular(OUTPUT/'analysis.json',1024*1024)),'exact outcome replay')
    A.exact(complete['outcome'],result['outcome'],'terminal outcome')
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    n = sub.add_parser('neutral'); n.add_argument('output',type=Path)
    r = sub.add_parser('receipt'); r.add_argument('neutral',type=Path); r.add_argument('output',type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt',type=Path)
    w = sub.add_parser('worker'); w.add_argument('claim')
    sub.add_parser('replay')
    args = parser.parse_args()
    if args.command=='neutral': neutral(args.output)
    elif args.command=='receipt': A.publish(args.output,A.canonical(receipt(args.neutral)))
    elif args.command=='run': run(args.receipt)
    elif args.command=='worker': worker(args.claim)
    else: print(json.dumps(replay(),sort_keys=True))
