"""Baseline-only uncovered-size inventory; no candidate, timing or promotion.

Numeric, packet-oracle and public-API logic is carried from the reviewed
Wh2UncoveredKRecoveryInventoryR0 without its historical producer/receipt chain.
"""
import ctypes as T
from functools import lru_cache
import hashlib
from pathlib import Path
import struct
import Support as A
N = A

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = 'wirehair.wh2.uncovered-band-inventory-r0'
QUALIFIED = Path('/tmp/wh2-small-unit-diagonal.MfPjFJaS')
LIBRARIES = {
    'native': (QUALIFIED/'native/libbaseline.so',
               'bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe'),
    'portable': (QUALIFIED/'portable/libbaseline.so',
                 '5f161aaa2e2356ef64f39dcbd913451101e9058cf26df91d39fd8b0431dc2452'),
}
KS, WIDTHS = (7, 9, 12, 16), (2, 64, 1280)
SCHEDULES = ('iid', 'burst', 'adversarial', 'repair-only')
MASK = (1 << 64) - 1
RAW_CAP = 128*1024**2


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



def neutral_roster():
    return [dict(group='neutral', k=k, b=b, tail=tail, schedule='systematic',
                 root=None, ids=list(range(k+4)))
            for k in KS for b in (17,65) for tail in (1,b)]


def check_library(report, mode):
    A.require(mode in LIBRARIES, 'declared baseline backend')
    path, expected = LIBRARIES[mode]
    A.exact(set(report), {'path','sha256','base','exports'}, 'public binding report schema')
    A.exact((report['path'],report['sha256']), (str(path),expected), 'baseline DSO identity')
    raw = A.read_regular(path,16*1024**2)
    A.exact(A.sha(raw),expected,'retained baseline DSO bytes')
    names = {'wirehair_init_', 'wirehair_encoder_create_owned_ex', 'wirehair_decoder_create_ex'}
    names.update('wirehair_'+name for name in ('encode','decode','recover','free'))
    names.update('wirehair_v2_'+name for name in
                 ('encoder_create','decoder_create','encode','decode','recover','free'))
    A.exact(set(report['exports']),names,'complete public binding roster')
    base = report['base']
    A.require(type(base) is int and 0 < base < 2**64,'reported load base')
    for record in report['exports'].values():
        A.exact(set(record),{'address','offset'},'public binding schema')
        A.require(type(record['offset']) is int and 0 < record['offset'] < len(raw) and
                  type(record['address']) is int and
                  record['address'] == base+record['offset'] < 2**64,'public binding bounds')
    # Actual ownership is checked with dladdr before and after native calls.
    # This retained-data check does not claim to rerun dynamic loading.


def expected_calls(records,cases):
    ids = id_sets(cases)
    creates = sum(key[0] for key in ids)+len(cases)
    encodes = sum(key[0]*len(values) for key,values in ids.items())+sum(len(r['ids']) for r in cases)
    return [[creates,encodes,len(cases),sum(len(r['arms'][a]['feed']) for r in records),
             sum(len(r['arms'][a]['recovered']) for r in records),creates+len(cases)] for a in (0,1)]


def verify(path,claim,mode='native',neutral=False):
    records = []; cases = neutral_roster() if neutral else roster()
    parity = hashlib.sha256()
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
        check_library(header['library'],mode)
        index = row_index(header['coefficients'],cases)
        parity.update(A.canonical({key:header[key] for key in ('type','protocol','coefficients')}))
        for ordinal,case in enumerate(cases):
            r = next_record()
            parity.update(A.canonical(r))
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
                if neutral: A.exact((x['first'],x['rank_first']),(k,k),'systematic neutral succeeds at K')
                del x['packets']
            records.append(r)
        footer = next_record()
        A.exact(footer,dict(type='footer',complete=True,records=len(cases),calls=expected_calls(records,cases)),'complete footer and API ledger')
        parity.update(A.canonical(footer))
        A.exact(stream.read(),b'','no trailing records')
    return dict(cases=len(cases), checked=True, scientific_launch=False, mode=mode,
                parity_sha256=parity.hexdigest()) if neutral else summarize(records)



def run_inventory(mode, cases, claim, emit):
    A.require(mode in LIBRARIES,'known inventory backend')
    library = A.Library(*LIBRARIES[mode])
    apis = [Api(library,arm) for arm in ('wh2','wh1')]
    observed = observe_rows(apis,cases)
    library.check_bindings()
    emit(dict(type='header',protocol=PROTOCOL,claim=claim,
              library=library.report,coefficients=observed))
    evaluate(apis,cases,observed,emit)
    library.check_bindings()
    emit(dict(type='footer',complete=True,records=len(cases),calls=[api.calls for api in apis]))
