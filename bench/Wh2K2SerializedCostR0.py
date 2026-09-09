#!/usr/bin/env python3
"""K2 lifecycle gate with deferred result I/O; one separately frozen launch."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import selectors
import re
import struct
from functools import lru_cache
import shlex
import subprocess
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = sibling('k2_serialized_cost_common', 'Wh2AlignedIntermediateCostR0.py')
ROOT = A.ROOT
PROTOCOL = 'wirehair.wh2.k2-serialized-cost-r0'
OUTPUT = Path('/var/tmp/wh2-k2-serialized-cost-r0')
U = sibling('k2_control_producer', 'Wh2K5PublicCostR0.py')
QUALIFIED = Path('/tmp/wh2-k5-public-recovery-qualified.bksSezzu')
SMALL = Path('/tmp/wh2-k2-serialized.rnxDJH4G')
K = 2
SEALED = Path('/var/tmp/wh2-k2-thue-morse-r0')
LIB_MANIFEST_SHA = dict(
    native='bb89d5d4551989dbd02d7128ca8bff12b02c135e08ea9e6d3f4667ae19eb4bff',
    scalar='2eaaa3df1dc30239fa970076b3cc63c8b4ac32996522f89920c8c8efbcedee5b',
    asan='4ae747bcc734699d704a871e9c3a58a186ea58fcdccdc6e23bfa6f56d01c704b')
# Existing CPU50 identity constants from Wh2RdpruTargetIdentityV2.cpp, not a
# post-run retarget. The binary serializer has always emitted 617 bytes.
TARGET_BYTES = 617
TARGET_SHA = '3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7'
ARCHIVES = U.ARCHIVES
PAIRS = ((0,0),(1,1),(2,2),(3,3),(4,4),(5,5),(0,1),(2,1),(3,4),(5,4))
ARM_NAMES = ('ordinary_independent','k2_independent','wh1_owned',
             'ordinary_borrowed','k2_borrowed','wh1_borrowed')
SIDES = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0)
WIDTHS = (2, 64, 1280)
SHAPES = ((2,2),(2,1),(64,64),(64,1),(1280,1280),(1280,1))
PACKETS = tuple(range(8))+tuple(0xffffffff-2*j for j in range(6))
CALLBACKS, BATCH = 77760, 128
RAW_CAP, ERR_CAP = 192*1024**2, 65536
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
NEW = ('bench/Wh2K2SerializedCostR0.cpp', 'bench/Wh2K2SerializedCostR0.py',
       'bench/test_Wh2K2SerializedCostR0.py',
       'bench/Wh2K2SerializedCostR0.md',
       'bench/Wh2AlignedIntermediateCostR0.py')



def pin(path):
    p = Path(path)
    owned = any(base == p or base in p.parents for base in (ROOT, Path('/tmp'), Path('/var/tmp')))
    return A.pin(p, installed=not owned)


def command(args):
    return subprocess.check_output(list(map(str, args)), timeout=60, cwd=ROOT)


def roster():
    index = 0
    for r in range(12):
        for s in range(2):
            for ws in range(6):
                for ms in range(3):
                    for cs in range(10):
                        order, width, metric = (r+s) % 2, (r+s+ws) % 6, (r+s+ws+ms) % 3
                        comparison = (2*r+s+ws+metric+cs) % 10
                        for p in range(18):
                            bin_number = r+12*(r % 4) if p < 2 else (r+6*((p-2)//8)) % 12+12*(((p-2) % 8)//2)
                            yield [index, r, order, width, metric, comparison, p,
                                   PAIRS[comparison][SIDES[p] ^ order], (2*bin_number+1)*1000000//96]
                            index += 1


def statistics(records):
    A.exact(len(records), CALLBACKS, 'complete statistical cohort')
    groups = {}
    for start in range(0, CALLBACKS, 18):
        panel, contrasts = records[start:start+18], []
        for j in range(8):
            values = {}
            for row in panel[2+2*j:4+2*j]:
                c, clocks = row['coordinate'], row['observation']['clocks']
                values[SIDES[c[6]] ^ c[2]] = clocks[3]-clocks[2]
            A.require(set(values) == {0, 1} and min(values.values()) > 0, 'paired positive durations')
            contrasts.append(math.log(values[1])-math.log(values[0]))
        c = panel[0]['coordinate']
        groups.setdefault((c[3], c[4], c[5], c[2]), []).append(math.fsum(contrasts)/8)
    A.exact(len(groups), 360, 'separate full-cost cells')
    results, failed_controls, failed_treatments = [], [], []
    bound = math.log1p(.02)
    for key, values in sorted(groups.items()):
        width, metric, comparison, order = key
        estimate = A.confidence(values)
        item = dict(width=SHAPES[width][0], tail=SHAPES[width][1], metric=metric, comparison=comparison, order=order,
                    estimate=estimate, replicate_logs=values,
                    comparison_arms=[ARM_NAMES[a] for a in PAIRS[comparison]])
        if comparison < 6:
            passed = -bound < estimate['lower95_log'] and estimate['upper95_log'] < bound
            item['control_pass'] = passed
            if not passed:
                failed_controls.append(list(key))
        else:
            item['treatment_pass'] = estimate['upper95_log'] < 0
            item['upper_ratio_limit'] = 1.0
            if not item['treatment_pass']:
                failed_treatments.append(list(key))
        results.append(item)
    return dict(outcome='CONTROL_FAIL' if failed_controls else 'FAIL' if failed_treatments else 'PASS',
                statistics=results, failed_controls=failed_controls, failed_treatments=failed_treatments,
                speed_qualified=not failed_controls and not failed_treatments,
                all_K_claimed=False, recovery_rate_claimed=False, production_promotion_claimed=False)



def clocks(observation, previous):
    A.exact(set(observation), {'clocks', 'before', 'after'}, 'observation keys')
    c = observation['clocks']
    A.require(len(c) == 6, 'clock shape')
    for value in c:
        A.integer(value)
    A.require(c[0] <= c[2] < c[3] <= c[5] and c[1] <= c[4] and c[4]-c[1] <= c[5]-c[0], 'clock ordering')
    for key in ('before', 'after'):
        A.exact(len(observation[key]), 4, 'counter shape')
        for value in observation[key]:
            A.integer(value)
    A.require(all(a <= b for a, b in zip(observation['before'], observation['after'])), 'counter ordering')
    if previous is not None:
        A.require(previous['clocks'][5] <= c[0] and previous['clocks'][4] <= c[1] and
                  all(a <= b for a, b in zip(previous['after'], observation['before'])), 'cross-record ordering')


def multiply(a, b):
    product = 0
    for bit in range(8):
        if b & (1 << bit):
            product ^= a << bit
    for bit in range(14, 7, -1):
        if product & (1 << bit):
            product ^= 0x14d << (bit-8)
    return product


@lru_cache(maxsize=1)
def selected_rows():
    # No candidate lookup or native generator is consulted.
    def product(a, b):
        out = [0]*(K*K)
        for r in range(K):
            for c in range(K):
                for j in range(K): out[K*r+c] ^= multiply(a[K*r+j], b[K*j+c])
        return out
    pair = []
    for feedback in ((2,3),(3,3)):
        matrix = [0]*(K*K)
        for i in range(K-1): matrix[K*(i+1)+i] = 1
        for i, value in enumerate(feedback): matrix[K*i+K-1] = value
        pair.append(matrix)
    levels = [pair]
    for _ in range(31):
        a, b = levels[-1]
        levels.append([product(a,b),product(b,a)])
    rows = []
    for packet_id in PACKETS:
        vector = [1]+[0]*(K-1)
        for bit in range(32):
            if packet_id & (1 << bit):
                matrix = levels[bit][bin(packet_id >> (bit+1)).count('1') % 2]
                out = [0]*K
                for r in range(K):
                    for c in range(K): out[r] ^= multiply(matrix[K*r+c],vector[c])
                vector = out
        rows.append(tuple(vector))
    return tuple(rows)


def packet_lengths(width, tail):
    return [tail if packet_id==K-1 else width for packet_id in PACKETS]


def payload(rows, width, tail):
    source = bytes((37*i+i//11) % 256 for i in range(width+tail))
    output = bytearray()
    for row, length in zip(rows,packet_lengths(width,tail)):
        for j in range(length):
            value = 0
            for k in range(K):
                offset = k*width+j
                if offset<len(source): value ^= multiply(row[k],source[offset])
            output.append(value)
        output.extend([0xa5]*(width-length))
    return output.hex()


@lru_cache(maxsize=6)
def candidate_fixture(width,tail):
    A.require((width,tail) in SHAPES, 'frozen shape')
    profile = struct.pack('<4sHHQQII',b'WHK2',1,32,0x5748324b32544d31,width+tail,width,0)
    return profile.hex(),payload(selected_rows(),width,tail)


def prior_header(build_dir):
    # Both captures use the exact binary-safe HeaderJson/Hex transport.
    header = A.decode(A.read_regular(build_dir/'fixtures.json', 1024*1024))
    target = A.decode(A.read_regular(build_dir/'target.json', 1024*1024))
    fields = {'type','protocol','claim','batch','identity_hex','prelude','fixtures'}
    for captured in (header, target):
        A.exact(set(captured), fields, 'neutral reference header fields')
        A.exact((captured['type'], captured['protocol'], captured['claim'], captured['batch']),
                ('header', PROTOCOL, '0'*64, BATCH), 'neutral reference identity')
        clocks(captured['prelude'], None)
    verify_fixtures(target, header)
    encoded = target['identity_hex']
    A.require(type(encoded) is str and len(encoded) == 2*TARGET_BYTES and
              re.fullmatch('[0-9a-f]+', encoded) is not None, 'complete binary target encoding')
    binary = bytes.fromhex(encoded)
    A.exact(A.sha(binary), TARGET_SHA, 'existing frozen CPU50 identity digest')
    header['identity_before'] = dict(canonical_hex=encoded)
    return header


def scalar_rank(rows):
    rows = [list(r) for r in rows]
    rank = 0
    for column in range(K):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = next(i for i in range(1, 256) if multiply(rows[rank][column], i) == 1)
        rows[rank] = [multiply(v, inverse) for v in rows[rank]]
        for i in range(len(rows)):
            if i != rank:
                factor = rows[i][column]
                rows[i] = [a ^ multiply(factor, b) for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def verify_fixtures(header, old):
    A.exact(len(header['fixtures']),6,'fixture count')
    A.exact(len(old['fixtures']),6,'reference fixture count')
    selected = selected_rows()
    for f, prior, (width,tail) in zip(header['fixtures'],old['fixtures'],SHAPES):
        A.exact(set(f),{'width','message','packet_bytes','source','arms'},'fixture fields')
        A.exact(f['width'],width,'fixture width')
        A.exact(f['message'],width+tail,'fixture message length')
        A.exact(f['packet_bytes'],packet_lengths(width,tail),'meaningful packet lengths')
        A.exact(f['source'],bytes((37*i+i//11) % 256 for i in range(width+tail)).hex(),'fixed K2 source')
        A.exact(len(f['arms']),6,'six ownership-matched fixtures')
        A.exact(f,prior,'pinned native neutral fixture')
        for a, arm in enumerate(f['arms']):
            A.exact(set(arm),{'profile','packets','steps','rows','fallback'},'arm fields')
            for key, size in (('profile',64),('packets',28*width),('rows',28*K)):
                A.require(type(arm[key]) is str and len(arm[key])==size and
                          re.fullmatch('[0-9a-f]+',arm[key]) is not None,'exact '+key+' hex')
            A.exact(len(arm['steps']),2,'decoder families')
            for count in arm['steps']: A.integer(count,K,K+6)
            if a in (0,3):
                profile = bytes.fromhex(arm['profile'])
                A.exact(profile[:16],struct.pack('<4sHHQ',b'WHV2',1,32,0x4b295bbb47f4f9c9),
                        'explicit certified WH2 descriptor')
                A.exact(profile[16:28],struct.pack('<QI',width+tail,width),'actual WH2 dimensions')
                A.exact(profile[29:],bytes(3),'certified reserved bytes')
            if a in (2,5): A.exact(arm['profile'],'00'*32,'WH1 has no descriptor')
            raw_rows = bytes.fromhex(arm['rows'])
            rows = tuple(tuple(raw_rows[i:i+K]) for i in range(0,14*K,K))
            A.exact(arm['packets'],payload(rows,width,tail),'independent native-row packet arithmetic')
            for family in range(2):
                slots = list(range(8,14) if family else range(K,K+6))+list(range(K))
                endpoint = next(i for i in range(K,K+7) if scalar_rank([rows[j] for j in slots[:i]])==K)
                A.exact(arm['steps'][family],endpoint,'independent first-success rank')
            A.exact(arm['fallback'],[n>6 for n in arm['steps']],'actual systematic fallback')
            if a in (1,4):
                p, packets = candidate_fixture(width,tail)
                A.exact(arm['profile'],p,'sealed K2 descriptor')
                A.exact(rows,selected,'independent sealed K2 rows')
                A.exact(arm['packets'],packets,'independent selected-pair payload oracle')
        for a in range(3): A.exact(f['arms'][a],f['arms'][a+3],'policy-independent complete fixture')


def verify(raw, claim, old):
    A.require(0 < len(raw) <= RAW_CAP and raw.endswith(b'\n'), 'bounded complete stream')
    rows = [A.decode(line) for line in raw.splitlines()]
    A.exact(len(rows), CALLBACKS+2, 'whole raw cohort')
    header, footer = rows[0], rows[-1]
    A.exact(set(header), {'type','protocol','claim','batch','identity_hex','prelude','fixtures'}, 'header schema')
    A.exact((header['type'], header['protocol'], header['claim'], header['batch']),
            ('header', PROTOCOL, claim, BATCH), 'header identity')
    A.exact(header['identity_hex'], old['identity_before']['canonical_hex'], 'frozen captured CPU identity')
    verify_fixtures(header, old)
    previous = header['prelude']; clocks(previous, None)
    total_work = 0
    for row, coordinate in zip(rows[1:-1], roster()):
        A.exact(set(row), {'type', 'coordinate', 'ready', 'target', 'wait', 'observation', 'counts',
                          'addresses', 'address_count', 'complete', 'checked'}, 'record fields')
        A.exact(row['type'], 'record', 'record type')
        A.exact(row['coordinate'], coordinate, 'fixed chronology')
        index, rep, order, width, metric, comparison, position, arm, q = coordinate
        del index, rep, comparison, position
        A.integer(row['ready']); A.integer(row['target'])
        A.require(previous['clocks'][5] <= row['ready'] and row['target'] == row['ready']+q, 'relative target')
        A.exact(len(row['wait']), 4, 'wait shape')
        w0, c0, w1, c1 = [A.integer(v) for v in row['wait']]
        A.require(row['ready'] <= w0 <= w1 and w1 >= row['target'] and c0 <= c1 and
                  previous['clocks'][4] <= c0, 'wait chronology')
        clocks(row['observation'], previous)
        observed = row['observation']['clocks']
        A.require(observed[0] >= w1 and observed[1] >= c1 and observed[2] >= row['target'], 'retained delayed start')
        previous = row['observation']; total_work += observed[3]-observed[2]
        steps = header['fixtures'][width]['arms'][arm]['steps'][metric-1] if metric else 0
        expected = [0 if metric else BATCH, 0 if metric else 14*BATCH,
                    BATCH if metric else 0, BATCH*steps, BATCH if metric else 0, BATCH]
        A.exact(row['counts'], expected, 'whole lifecycle ledger')
        A.exact(row['address_count'], BATCH, 'all handles retained')
        A.exact(len(row['addresses']), BATCH, 'address roster')
        for address in row['addresses']:
            A.integer(address, 1)
        A.exact(row['complete'], True, 'every call completed')
        A.exact(row['checked'], True, 'every output checked')
    A.exact(footer, dict(type='footer', complete=True, records=CALLBACKS, work_ns=total_work), 'complete footer')
    A.require(total_work <= 180000000000, 'inner work cap')
    return statistics(rows[1:-1])


def verify_publication(raw, mode, old):
    """Three neutral real-codec captures; never a scientific timing cohort."""
    A.require(mode in ('success','last-recover','throw-recover','last-clock','last-source'),'neutral mode')
    A.require(0 < len(raw) < 1024*1024 and raw.endswith(b'\n'),'neutral stream cap')
    rows = [A.decode(line) for line in raw.splitlines()]
    A.exact(len(rows),5,'all neutral captures and footer')
    header,footer = rows[0],rows[-1]
    A.exact(set(header),{'type','protocol','claim','batch','identity_hex','prelude','fixtures'},'neutral header fields')
    A.exact((header['type'],header['protocol'],header['claim'],header['batch'],header['identity_hex']),
            ('header',PROTOCOL,'0'*64,BATCH,b'neutral-deferred'.hex()),'neutral identity')
    verify_fixtures(header,old)
    previous = header['prelude']; clocks(previous,None); total = 0
    coordinates = list(roster())
    for i,row in enumerate(rows[1:-1]):
        A.exact(set(row),{'type','coordinate','ready','target','wait','observation','counts',
                         'addresses','address_count','complete','checked'},'neutral record fields')
        coordinate = coordinates[i*180]
        A.exact(row['type'],'record','neutral record type')
        A.exact(row['coordinate'],coordinate,'neutral coordinate')
        A.exact((row['ready'],row['target'],row['wait']),(0,0,[0]*4),'no neutral phase targeting')
        metric,width,arm = coordinate[4],coordinate[3],coordinate[7]
        steps = header['fixtures'][width]['arms'][arm]['steps'][metric-1] if metric else 0
        A.exact(row['counts'],[0 if metric else BATCH,0 if metric else 14*BATCH,
                              BATCH if metric else 0,BATCH*steps,BATCH if metric else 0,BATCH],
                'complete attempted last-call ledger')
        A.exact(row['address_count'],BATCH,'all neutral handles')
        A.exact(len(row['addresses']),BATCH,'neutral address shape')
        for address in row['addresses']: A.integer(address,1)
        failed = i==2 and mode!='success'
        partial = i==2 and mode=='last-clock'
        A.exact(row['checked'],not failed,'neutral byte check status')
        A.exact(row['complete'],not failed or partial or mode=='last-source','neutral completion status')
        o = row['observation']
        A.exact(o['before'],[0]*4,'neutral before counters')
        A.exact(o['after'],[0]*4,'neutral after counters')
        if partial:
            A.exact(set(o),{'clocks','before','after'},'partial observation schema')
            c = o['clocks']; A.exact(len(c),6,'partial clock shape')
            for value in c: A.integer(value)
            A.require(previous['clocks'][5]<=c[0]<=c[2]<c[3] and
                      previous['clocks'][4]<=c[1]<=c[4] and c[5]==0,'partial final clock retention')
            A.exact(o['before'],[0]*4,'partial before counters')
            A.exact(o['after'],[0]*4,'partial after counters')
        else:
            clocks(o,previous); total+=o['clocks'][3]-o['clocks'][2]
        previous=o
    A.exact(footer,dict(type='footer',complete=mode=='success',records=3,work_ns=total),'neutral failure footer')
    return dict(neutral=True,records=3,mode=mode,speed_qualified=False)


def qualified_inputs(mode):
    """Reuse closed producer evidence; never invoke a spent verifier or build."""
    path = QUALIFIED/mode/'manifest.json'
    raw = A.read_regular(path,2*1024**2)
    A.exact(A.sha(raw),LIB_MANIFEST_SHA[mode],'retained control producer manifest')
    manifest = A.decode(raw)
    declared = {p['path']:p for p in manifest['inputs']+manifest['artifacts']}
    proof = QUALIFIED/mode/'qualified-library.json'
    A.exact(pin(proof),declared[str(proof)],'bound library provenance')
    provenance = A.decode(A.read_regular(proof,2*1024**2))
    A.exact((provenance['source_head'],provenance['mode'],provenance['library_source_provenance_closed']),
            (U.SOURCE_HEAD,mode,True),'closed original producer')
    archive = Path(provenance['archive']['path'])
    A.exact(pin(archive),provenance['archive'],'archive identity')
    A.exact(pin(archive)['sha256'],ARCHIVES[mode],'qualified actual WH1/ordinary WH2')
    objects,sources = U.library_commands(provenance['commands'].encode(),mode)
    A.exact(command(['/usr/bin/ar','t',archive]).decode().splitlines(),[p.name for p in objects],'complete archive roster')
    paths = {path,proof,archive,Path('/usr/bin/git').resolve(strict=True)}
    for record,obj,source in zip(provenance['members'],objects,sources):
        A.exact(pin(obj),record['object'],'bound producer object')
        A.exact(A.sha(command(['/usr/bin/ar','p',archive,obj.name])),pin(obj)['sha256'],'archive/object agreement')
        paths.add(obj)
        paths.update(U.library_dependencies(record['dependencies'].encode(),
            obj.relative_to(U.QUALIFIED/U.LIB_DIRS[mode]),source,obj.stat().st_mtime_ns))
    A.exact(len(provenance['members']),18,'all control translation units')
    snapshots = {}
    for p in sorted(paths-{path,proof}):
        expected = declared[str(p)]
        if pin(p)!=expected:
            A.exact(str(p.relative_to(ROOT)),'codec/WirehairSmallCore.h','sole changed producing source')
            old = command(['git','cat-file','blob',U.SOURCE_HEAD+':codec/WirehairSmallCore.h'])
            A.exact(dict(path=str(p),bytes=len(old),sha256=A.sha(old)),expected,'historical source bytes')
            A.exact(pin(p)['sha256'],'11fe03c80d705e29bd694ab9795ae21ee6379a5320dfab12f689f3aff9806381','qualified K2 guards only')
            snapshots[p.name] = old
        else:
            A.exact(pin(p),expected,'unchanged producer input')
    log = Path(provenance['qualification_log']['path'])
    A.exact(pin(log),provenance['qualification_log'],'preserved control qualification log'); paths.add(log)
    # The native three production users were independently recompiled at K2
    # admission. Their bytes also equal this exact control archive's members.
    identity_path = Path('/tmp/wh2-k2-native.jtQp2mfD/object-identity.json')
    A.exact(pin(identity_path)['sha256'],'539349f52389f7d30521a9b2d8f6d353938b4a4816038715fa281239ad899115','K2 code-generation proof')
    paths.add(identity_path)
    identity = A.decode(A.read_regular(identity_path,1024*1024))
    if mode=='native':
        for name in ('WirehairSmall.cpp','WirehairSmallK5.cpp','WirehairV2Profile.cpp'):
            record = next(r for r in identity['objects'] if Path(r['path']).name=='shared-'+name+'-normal.o')
            p = Path(record['path']); paths.add(p)
            A.exact(pin(p)['sha256'],record['sha256'],'retained normal production object')
            A.exact(A.sha(command(['/usr/bin/ar','p',archive,name+'.o'])),record['sha256'],'current production object equals archive')
    qualification = SMALL/'QUALIFIED.json'
    A.exact(pin(qualification)['sha256'],'616ecfb9727a4ad48990d3fa00cf8b559ca2f6fd5e767c83c906269b2fbbdf6c','qualified K2 boundary')
    paths.add(qualification)
    q = A.decode(A.read_regular(qualification,1024*1024))
    for name,item in q['files'].items():
        p = SMALL/name
        A.exact(dict(bytes=pin(p)['bytes'],sha256=pin(p)['sha256']),item,'qualified boundary artifact')
        paths.add(p)
    report = A.decode(A.read_regular(SMALL/(mode+'-qualification.json'),1024*1024))
    for name,digest in report['source_hashes'].items():
        A.exact(pin(ROOT/name)['sha256'],digest,'unchanged qualified boundary source'); paths.add(ROOT/name)
    paths.update(SEALED/name for name in ('COMPLETE.json','CLAIM.json','raw.json','summary.json','stderr.txt'))
    return [archive],paths,provenance,snapshots


def build(mode, output):
    output = output.parent.resolve(strict=True)/output.name
    A.require(output.is_absolute() and ROOT not in output.parents and output != ROOT and
              not output.exists() and not output.is_symlink(), 'fresh external build')
    archives, dependencies, provenance, snapshots = qualified_inputs(mode)
    frozen = {}; U.freeze_inputs(dependencies,frozen)
    output.mkdir(mode=0o700)
    A.publish(output/'qualified-library.json',A.canonical(provenance))
    for name,raw in snapshots.items(): A.publish(output/('qualified-'+name),raw)
    boundary = output/'boundary'
    args = ['cmake','-S',str(ROOT/'bench/Wh2SmallNative'),'-B',str(boundary),
            '-DWH2_SMALL_LIBRARY='+str(archives[0]),'-DWH2_SMALL_TEST_DIMENSION=2',
            '-DCMAKE_BUILD_TYPE='+('Debug' if mode=='asan' else 'Release')]
    if mode=='scalar': args.append('-DWH2_SMALL_PORTABLE=ON')
    if mode=='asan':
        flags = '-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -march=native'
        args += ['-DCMAKE_C_FLAGS='+flags,'-DCMAKE_CXX_FLAGS='+flags,
                 '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined']
    for tag,argv in (('configure',args),('build',['cmake','--build',str(boundary),'-j4']),
                     ('test',['ctest','--test-dir',str(boundary),'-V','--timeout','120'])):
        raw = command(argv); A.publish(output/('boundary-'+tag+'.txt'),raw)
        if tag=='test': A.require(b'100% tests passed, 0 tests failed out of 7' in raw,'current-GF boundary qualification')
    archives.insert(0,boundary/'libwh2_small_serialized.a')
    dependencies.update(p for p in boundary.rglob('*') if p.is_file())
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_SMALL_CODEC_K=2', '-DWH2_K2_SERIALIZED_COST_NEUTRAL='+str(int(mode != 'native')),
             '-I'+str(ROOT), '-I'+str(ROOT/'bench'), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += ['-O1', '-g', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-march=native'] if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags += ['-DANDROID=1']
    sources = [ROOT/NEW[0]]+[ROOT/'bench'/name for name in
               ('Wh2FrozenTrace.cpp', 'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp')]
    commands, objects = [], []
    dependencies.update(ROOT/n for n in NEW)
    dependencies.add(ROOT/'bench/Wh2K5PublicCostR0.py')
    for name in ('c++','cc','as','ld','nm','ar','ldd','bash','cmake','make'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for compiler in ('cc1','cc1plus','collect2'):
        dependencies.add(Path(command(['/usr/bin/c++','-print-prog-name='+compiler]).decode().strip()).resolve(strict=True))
    U.freeze_inputs(dependencies,frozen)
    dependencies.add(ROOT/'bench/Wh2AlignedIntermediateCostR0.py')
    for source in sources:
        obj = output/(source.stem+'.o'); dep = output/(source.stem+'.d')
        args = ['/usr/bin/c++']+flags+['-M','-MT',str(obj),'-MF',str(dep),str(source)]
        command(args); commands.append(args)
        before = U.preprocessor_dependencies(dep.read_bytes(),obj)
        dependencies.update(before); U.freeze_inputs(dependencies,frozen)
        args = ['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)]
        command(args); commands.append(args); objects.append(obj)
        after = U.preprocessor_dependencies(dep.read_bytes(),obj)
        A.exact(sorted(after),sorted(before),'actual compiler dependency closure')
        U.freeze_inputs(dependencies,frozen)
    executable = output/'cost_worker'
    args = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan':
        args += ['-fsanitize=address,undefined']
    args += list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)]
    command(args); commands.append(args)
    loads = [line[5:] for line in (output/'link.map').read_text().splitlines() if line.startswith('LOAD ')]
    A.require(loads and all(Path(p).is_absolute() for p in loads),'absolute linker input roster')
    dependencies.update(Path(p).resolve(strict=True) for p in loads if output not in Path(p).resolve(strict=True).parents)
    symbols = command(['/usr/bin/nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    A.require(names.count('GF256Ctx') == names.count('gf256_init_') == 1, 'single GF runtime')
    A.require(all(names.count(prefix+name) == 1 for prefix in ('wirehair_v2_', 'wh2_small_') for name in
                  ('encode', 'decode', 'free')), 'actual WH2 and qualified serialized K2 APIs')
    A.require(all(names.count(n) == 1 for n in ('wirehair_encode', 'wirehair_decode', 'wirehair_free')), 'actual WH1 API')
    # Check the compile-ISA-dependent private context before any GF access.
    sizes = command(['/usr/bin/nm','-S','--defined-only',executable]).decode().splitlines()
    context_size = [int(line.split()[1],16) for line in sizes if line.split()[-1]=='GF256Ctx']
    A.exact(len(context_size),1,'one sized GF context')
    contract = command([executable,'--contract'])
    A.exact(A.decode(contract),dict(K=K,batch=BATCH,callbacks=CALLBACKS,cpu_seconds=240,
            wall_seconds=300,work_seconds=180,address_space_mib=384,
            claim_path=str(OUTPUT/'CLAIM.json'),neutral_only=mode!='native',
            gf_context_bytes=context_size[0]),'compiled resource/launch/GF ABI contract')
    A.publish(output/'contract.json',contract)
    text = command(['/usr/bin/nm', '-C', executable]).decode().splitlines()
    A.require(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line for line in text) == 1, 'one common WORK body')
    for name in ('c++', 'cc', 'as', 'ld', 'nm', 'ar', 'make'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for compiler in ('cc1', 'cc1plus', 'collect2'):
        dependencies.add(Path(command(['/usr/bin/c++', '-print-prog-name='+compiler]).decode().strip()).resolve(strict=True))
    tool_binaries = [p for p in dependencies if os.access(str(p),os.X_OK) and
                     A.read_regular(p,256*1024*1024,installed=True)[:4]==b'\x7fELF']
    for target in tool_binaries+[executable, Path(sys.executable).resolve(strict=True)]:
        linked = command(['/usr/bin/ldd', target]).decode()
        dependencies.update(Path(word).resolve(strict=True) for word in linked.split() if word.startswith('/'))
    dependencies.add(Path(sys.executable).resolve(strict=True))
    fixture_raw = command([executable, '--neutral-fixtures'])
    header = A.decode(fixture_raw)
    verify_fixtures(header, header)
    A.publish(output/'fixtures.json', fixture_raw)
    if mode == 'native':
        A.publish(output/'target.json', command([executable, '--neutral-target']))
        prior_header(output)
    # The very same Authenticate helper gates Worker and this positive test.
    # Do not create a scientific namespace merely to test launch binding.
    A.exact(command([executable,'--claim-path']),str(OUTPUT/'CLAIM.json').encode()+b'\n','compiled claim path')
    neutral_claim = output/'neutral-claim.json'
    neutral_bytes = A.canonical(dict(protocol=PROTOCOL,neutral=True))
    A.publish(neutral_claim,neutral_bytes)
    A.exact(command([executable,'--neutral-claim',neutral_claim,A.sha(neutral_bytes)]),
            b'PASS claim authentication\n','positive claim authentication')
    negatives = []
    for argv in ([executable,'--neutral-claim',neutral_claim,'0'*64],
                 [executable,'--neutral-claim',output/'absent',A.sha(neutral_bytes)],
                 [executable,'--neutral-claim',neutral_claim,'g'*64],
                 [executable],[executable,'--unknown'],[executable,'--worker','0']):
        p = subprocess.run(list(map(str,argv)),cwd=ROOT,stdin=subprocess.DEVNULL,
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10)
        A.require(p.returncode==1 and not p.stdout and p.stderr.startswith(b'INVALID:'),'negative CLI/authentication')
        negatives.append(dict(argv=list(map(str,argv)),code=p.returncode,stderr=p.stderr.decode()))
    A.publish(output/'negative-cli.json',A.canonical(negatives))
    A.publish(output/'neutral.txt',command([executable,'--neutral']))
    for case in ('success','last-recover','throw-recover','last-clock','last-source'):
        raw = command([executable,'--neutral-publication',case])
        verify_publication(raw,case,header)
        A.publish(output/('publication-'+case+'.jsonl'),raw)
    # Real stdio failures, not only a mock sink. Never launch scientific WORK.
    failures = []
    with open('/dev/full','wb') as full:
        p = subprocess.run([str(executable),'--neutral-publication','success'],
                           stdout=full,stderr=subprocess.PIPE,stdin=subprocess.DEVNULL,timeout=60)
    A.require(p.returncode==1 and b'output stream' in p.stderr,'full output device rejects publication')
    failures.append(dict(sink='full',code=p.returncode,stderr=p.stderr.decode()))
    read_fd, write_fd = os.pipe()
    os.close(read_fd)
    try:
        p = subprocess.run([str(executable),'--neutral-publication','success'],
                           stdout=write_fd,stderr=subprocess.PIPE,stdin=subprocess.DEVNULL,timeout=60)
    finally:
        os.close(write_fd)
    A.require(p.returncode==1 and b'output stream' in p.stderr,'broken pipe is an explicit output error')
    failures.append(dict(sink='broken-pipe',code=p.returncode,stderr=p.stderr.decode()))
    A.publish(output/'publication-output-errors.json',A.canonical(failures))
    # Include the actual interpreter imports used by this builder/controller.
    dependencies.update(Path(m.__file__).resolve(strict=True) for m in list(sys.modules.values())
                        if getattr(m,'__file__',None) and Path(m.__file__).is_file())
    U.freeze_inputs(dependencies,frozen)
    manifest = dict(protocol=PROTOCOL, mode=mode, commands=commands,
                    inputs=[pin(p) for p in sorted(dependencies)], artifacts=[pin(p) for p in sorted(output.iterdir()) if p.is_file()])
    A.publish(output/'manifest.json', A.canonical(manifest))
    print(json.dumps(dict(mode=mode, executable=str(executable), scientific_launch=False)))


def current(receipt):
    A.exact(set(receipt), {'protocol','head','executable','environment','pins'}, 'receipt schema')
    A.exact(receipt['protocol'], PROTOCOL, 'receipt protocol')
    A.exact(receipt['environment'], {k: os.environ.get(k) for k in ENV_KEYS}, 'allocator environment')
    A.exact(receipt['environment'], {k: None for k in ENV_KEYS}, 'ordinary allocator policy')
    A.exact(command(['git', 'rev-parse', 'HEAD']).decode().strip(), receipt['head'], 'source HEAD')
    declared = {p['path']: p for p in receipt['pins']}
    A.exact(len(declared), len(receipt['pins']), 'unique receipt pins')
    for p in receipt['pins']:
        A.exact(pin(Path(p['path'])), p, 'unchanged receipt input')
    executable = Path(receipt['executable'])
    A.require(executable.is_absolute() and executable.name == 'cost_worker' and
              executable.parent.name == 'native' and str(executable) in declared, 'pinned executable binding')
    path = executable.parent/'manifest.json'
    A.require(str(path) in declared, 'pinned build manifest')
    manifest = A.decode(A.read_regular(path, 1024*1024))
    A.exact((manifest['protocol'], manifest['mode']), (PROTOCOL, 'native'), 'bound native build')
    A.require(str(executable) in {p['path'] for p in manifest['artifacts']}, 'built executable artifact')
    closure = {}
    for p in manifest['inputs']+manifest['artifacts']+[declared[str(path)]]:
        if p['path'] in closure:
            A.exact(closure[p['path']], p, 'shared manifest input')
        closure[p['path']] = p
    A.exact(receipt['pins'], sorted(closure.values(), key=lambda p: p['path']), 'complete manifest closure')
    prior_header(executable.parent)


def receipt(build_dir):
    manifest = A.decode(A.read_regular(build_dir/'manifest.json', 1024*1024))
    A.exact((manifest['protocol'], manifest['mode']), (PROTOCOL, 'native'), 'native full cost build')
    pins = {p['path']: p for p in manifest['inputs']+manifest['artifacts']}
    pins[str(build_dir/'manifest.json')] = pin(build_dir/'manifest.json')
    for p in pins.values():
        A.exact(pin(Path(p['path'])), p, 'current build/input')
        path = Path(p['path'])
        if ROOT in path.parents:
            A.exact(path.read_bytes(), command(['git', 'cat-file', 'blob', 'HEAD:'+str(path.relative_to(ROOT))]), 'committed source')
    result = dict(protocol=PROTOCOL, head=command(['git', 'rev-parse', 'HEAD']).decode().strip(),
                  executable=str(build_dir/'cost_worker'), environment={k: os.environ.get(k) for k in ENV_KEYS},
                  pins=sorted(pins.values(), key=lambda p: p['path']))
    A.require(all(v is None for v in result['environment'].values()), 'ordinary allocator policy')
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


def run(receipt_path):
    begin = time.monotonic()
    frozen = A.read_regular(receipt_path, 1024*1024); claimed = A.decode(frozen)
    A.exact(frozen, A.canonical(claimed), 'canonical receipt'); current(claimed)
    os.mkdir(str(OUTPUT), 0o700); A.publish(OUTPUT/'CLAIM.json', frozen)
    analysis = dict(outcome='INVALID', failure=None)
    try:
        raw, error, code, failure = capture(claimed['executable'], A.sha(frozen),
                                          min(begin+480, time.monotonic()+330), [OUTPUT/'raw.jsonl', OUTPUT/'stderr.txt'])
        A.require(failure is None and code == 0 and error == b'', 'worker/observer failure: '+str(failure))
        prior = prior_header(Path(claimed['executable']).parent)
        analysis = verify(raw, A.sha(frozen), prior)
        current(claimed); A.require(time.monotonic()-begin < 480, 'whole controller deadline')
    except Exception as error:
        analysis = dict(outcome='INVALID', failure=str(error))
    analysis.update(protocol=PROTOCOL, elapsed_seconds=time.monotonic()-begin)
    A.publish(OUTPUT/'analysis.json', A.canonical(analysis))
    members = [pin(p) for p in sorted(OUTPUT.iterdir())]
    A.require(sum(p['bytes'] for p in members) < 256*1024**2-65536, 'bundle cap')
    A.publish(OUTPUT/'COMPLETE.json', A.canonical(dict(protocol=PROTOCOL, outcome=analysis['outcome'], files=members)))
    print(json.dumps({k: v for k, v in analysis.items() if k not in ('statistics',)}, sort_keys=True))


def replay():
    complete = A.decode(A.read_regular(OUTPUT/'COMPLETE.json',1024*1024))
    A.exact(complete['protocol'],PROTOCOL,'sealed protocol')
    names = {'CLAIM.json','raw.jsonl','stderr.txt','analysis.json'}
    A.exact({p.name for p in OUTPUT.iterdir()},names|{'COMPLETE.json'},'sealed complete roster')
    A.exact({Path(p['path']).name for p in complete['files']},names,'sealed manifest roster')
    A.exact(len(complete['files']),4,'sealed manifest uniqueness')
    for p in complete['files']:
        A.exact(str(OUTPUT/Path(p['path']).name),p['path'],'sealed member path')
        A.exact(pin(Path(p['path'])),p,'sealed member identity')
    frozen = A.read_regular(OUTPUT/'CLAIM.json',1024*1024)
    receipt = A.decode(frozen); A.exact(frozen,A.canonical(receipt),'canonical sealed receipt'); current(receipt)
    stored = A.decode(A.read_regular(OUTPUT/'analysis.json',4*1024*1024))
    A.exact(stored['protocol'],PROTOCOL,'analysis protocol')
    A.exact(stored['outcome'],complete['outcome'],'sealed outcome')
    if stored['outcome']!='INVALID':
        A.exact(A.read_regular(OUTPUT/'stderr.txt',ERR_CAP),b'','successful worker stderr')
        actual = verify(A.read_regular(OUTPUT/'raw.jsonl',RAW_CAP),A.sha(frozen),
                        prior_header(Path(receipt['executable']).parent))
        A.exact({k:v for k,v in stored.items() if k not in ('protocol','elapsed_seconds')},
                actual,'complete independent-language replay')
    return stored


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('mode', choices=tuple(ARCHIVES)); b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    sub.add_parser('replay')
    args = parser.parse_args()
    if args.command == 'build':
        build(args.mode, args.output)
    elif args.command == 'receipt':
        A.publish(args.output, A.canonical(receipt(args.build_dir.resolve(strict=True))))
    elif args.command == 'run':
        run(args.receipt)
    else:
        print(json.dumps({k:v for k,v in replay().items() if k!='statistics'},sort_keys=True))
