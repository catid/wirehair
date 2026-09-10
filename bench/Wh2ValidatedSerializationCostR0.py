#!/usr/bin/env python3
"""One frozen lifecycle comparison of the validated-local serialization overlay."""
import argparse
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path
import struct
import sys

def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name,Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

B = sibling('_serialization_cost_build','Wh2ValidatedSerializationCostBuildR0.py')
D = sibling('_serialization_cost_deferred','Wh2CurrentPreservedDeferredCostR0.py')
R = B.T.runtime(B.PREPARED/'Wh2ValidatedSerializationDerivedR0.py')
A, ROOT = B.A, B.ROOT
PROTOCOL = 'wirehair.wh2.validated-serialization-cost-r0'
OUTPUT = Path('/var/tmp/wh2-validated-serialization-cost-r0')
# Filled from the sole byte-verified preparation, before any worker qualification.
CANDIDATE_SHA = '91aaaad2e752b304d8b0c1fdca7113051651eed4b9a6344ba6fe27be1c07298f'
PREPARATION_SHA = '5c789653460989147291d1fc02c958beffb156cb29100987ead9cac2f4dee4f5'
CASES = R.CASES + tuple((4,k,b,p) for k in (3,5,8) for b in (2,64,1280) for p in (1,2))
SOURCES = D.NEW + ('bench/Wh2ValidatedSerializationCostR0.py',
    'bench/Wh2ValidatedSerializationCostBuildR0.py','bench/test_Wh2ValidatedSerializationCostR0.py',
    'bench/Wh2ValidatedSerializationCostR0.md','bench/Wh2K8PublicCostBuildR0.py',
    'bench/Wh2ValidatedSerialization/CMakeLists.txt','bench/Wh2ValidatedSerialization/Overlay.cmake',
    'bench/Wh2ProfileClassification/CheckParity.cmake',
    'bench/Wh2ProfileClassification/ValidationTest.cpp',
    'bench/Wh2ValidatedSerializationRuntimeR0.py','bench/Wh2ValidatedSerialization/README.md',
    'test/V2SmallCodecTest.cpp','codec/V2ProfileTest.cpp','test/V2CertifiedCompatibility.cpp',
    str(B.COMMON),str(B.WORKER),str(B.PREPARED/'Wh2ValidatedSerializationDerivedR0.py'))
PROFILE_IDS = {3:0x67c1043ecaa9e184,5:0x80070c81bfe375f1,8:0x7a9276b85c730ae0}

def packet_ids(k): return tuple(range(k+8))+tuple(0xffffffff-2*j for j in range(6))

@lru_cache(maxsize=256)
def products(a): return bytes(R.O.multiply(a,b) for b in range(256))

@lru_cache(maxsize=3)
def rows(k):
    A.require(k in (3,5,8),'explicit installed small dimension')
    def multiply(a,b):
        result = [0]*(k*k)
        for i in range(k):
            for j in range(k):
                for c in range(k): result[k*i+j] ^= products(a[k*i+c])[b[k*c+j]]
        return tuple(result)
    taps = {3:(8,14,7),5:(121,110,207,198,31),8:(96,19,186,153,85,252,7,255)}[k]
    matrices = []
    for phase in range(2):
        matrix = [0]*(k*k)
        for i in range(k-1): matrix[k*(i+1)+i] = 1
        for i,v in enumerate(taps): matrix[k*i+k-1] = v ^ ((2 if k==8 else 1)*phase if i==0 else 0)
        matrices.append(tuple(matrix))
    levels = [tuple(matrices)]
    for _ in range(31):
        a,b = levels[-1]; levels.append((multiply(a,b),multiply(b,a)))
    result = []
    for packet in packet_ids(k):
        row = (1,)+(0,)*(k-1)
        for bit in range(32):
            if packet & (1<<bit):
                matrix = levels[bit][bin(packet>>(bit+1)).count('1')%2]
                next_row = [0]*k
                for i in range(k):
                    for j in range(k): next_row[i] ^= products(matrix[k*i+j])[row[j]]
                row = tuple(next_row)
        result.append(row)
    return tuple(result)

def first_success(incoming,k):
    basis = {}
    for index,original in enumerate(incoming):
        row = list(original)
        for column in range(k):
            if not row[column]: continue
            if column not in basis:
                basis[column] = row
                if len(basis)==k: return index+1
                break
            pivot = basis[column]; a,b = products(pivot[column]),products(row[column])
            row = [a[x]^b[y] for x,y in zip(row,pivot)]
    raise ValueError('small fixture never reaches full rank')

@lru_cache(maxsize=9)
def small_fixture(k,b):
    source = bytes((37*i+i//11)%256 for i in range(k*b))
    packets = bytearray()
    for row in rows(k):
        for j in range(b):
            value = 0
            for c in range(k): value ^= products(row[c])[source[c*b+j]]
            packets.append(value)
    profile = struct.pack('<4sHHQQII',b'WHV2',1,32,PROFILE_IDS[k],k*b,b,0)
    r = rows(k)
    return dict(profile=profile.hex(),packets=packets.hex(),steps=first_success(r[k:]+r[:k],k))

def verify_header(header,order,claim,meta,protocol=PROTOCOL):
    A.exact(len(header['fixtures']),38,'complete current-plus-small fixture roster')
    # The preserved 20 equations still match their independent historical data.
    R.verify_header(dict(header,fixtures=header['fixtures'][:20]),order,claim,meta,protocol)
    for fixture,case in zip(header['fixtures'][20:],CASES[20:]):
        _,k,b,_ = case
        A.exact(set(fixture),{'case','batch','source','arms'},'small fixture schema')
        A.exact(fixture['case'],list(case),'actual installed small route')
        A.exact(fixture['batch'],128,'frozen small batch')
        A.exact(fixture['source'],bytes((37*i+i//11)%256 for i in range(k*b)).hex(),'small source')
        A.exact(fixture['arms'],[small_fixture(k,b)]*2,'independent profile, full packet menu and prefix rank')

def combine(results,protocol=PROTOCOL):
    result = R.combine(results,protocol); improvements = []
    for case in R.CASES:
        if case[0]!=0: continue
        for metric in (0,1):
            cells = [s for load in results for s in load['statistics']
                     if s['case']==list(case) and s['metric']==metric and s['comparison']==2]
            A.exact(len(cells),4,'four separate benefit cells')
            A.exact([s['order'] for s in cells],[0,1,0,1],'both measurement orders in both load orders')
            if all(s['estimate']['upper95_log']<0 for s in cells): improvements.append(dict(case=list(case),metric=metric))
    result.update(certified_four_cell_improvements=improvements,
        candidate_retained=result['outcome']=='PASS' and bool(improvements),pre_admission_restoration_qualified=False)
    return result

def provenance(proof_dir=None):
    result,inputs = B.provenance(PREPARATION_SHA,proof_dir,SOURCES)
    A.exact(result[1]['original']['sha256'],CANDIDATE_SHA,'prospectively fixed candidate DSO')
    return result,inputs

def qualify(executable,output,meta,mode):
    for name in ('old-new','new-old'):
        A.exact(A.read_regular(output/('neutral-'+name+'.txt'),65536),
                b'PASS neutral98496-coordinate roster, 304 native WORK cases; no timing\n','complete neutral roster')
    D.qualify(executable,output,meta,mode,PROTOCOL,verify_header)

def settings():
    A.require(len(CANDIDATE_SHA)==64,'pin prepared candidate before building observer')
    return R.Configuration(PROTOCOL,OUTPUT,((B.DSO,B.DSO_SHA),(B.PREPARED/'libwirehair.so.2.0.0',CANDIDATE_SHA)),
        SOURCES,provenance,str(B.WORKER),qualify,CASES,verify_header,combine)

def enter_clean_environment():
    """Delegated subprocesses inherit only this explicit, checked environment."""
    environment = B.P.process_environment()
    if dict(os.environ)!=environment:
        os.execve(sys.executable,[sys.executable,str(Path(__file__).resolve())]+sys.argv[1:],environment)
        raise RuntimeError('execve unexpectedly returned')
    A.exact(dict(os.environ),environment,'exact controller and inherited worker environment')

def main():
    enter_clean_environment()
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='command',required=True)
    sub.add_parser('prepare')
    b=sub.add_parser('build'); b.add_argument('mode',choices=('native','asan-driver')); b.add_argument('output',type=Path)
    r=sub.add_parser('receipt'); r.add_argument('build_dir',type=Path); r.add_argument('output',type=Path)
    r=sub.add_parser('run'); r.add_argument('receipt',type=Path)
    sub.add_parser('replay'); args=parser.parse_args()
    if args.command=='prepare': B.prepare(); return
    cfg=settings()
    if args.command=='build': B.build(args.mode,args.output,cfg,R)
    elif args.command=='receipt': A.publish(args.output,A.canonical(R.receipt(args.build_dir,cfg)))
    elif args.command=='run': R.run(args.receipt,cfg)
    else:
        result = R.replay(cfg)
        print(json.dumps(dict(outcome=result['outcome'],candidate_retained=result['candidate_retained'],exact_replay=True)))

if __name__=='__main__': main()
