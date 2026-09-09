#!/usr/bin/env python3
"""One prospective current/certified-fall-through cost and small-path retention screen."""
import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import shlex
import struct

SPEC = importlib.util.spec_from_file_location('fallthrough_deferred',
    Path(__file__).with_name('Wh2CurrentPreservedDeferredCostR0.py'))
D = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(D)
R, U, H, A, ROOT = D.R, D.U, D.H, D.A, D.ROOT
PROTOCOL = 'wirehair.wh2.certified-fallthrough-cost-r0'
OUTPUT = Path('/var/tmp/wh2-certified-fallthrough-cost-r0')
PREPARED = Path('/tmp/wh2-fallthrough-cost-inputs.rwdQ1q8v')
NEUTRAL = Path('/tmp/wh2-certified-fallthrough-neutral.NdEtjlM4')
OVERLAY = NEUTRAL/'WirehairV2ProfileFallthrough.cpp'
OBJECT = NEUTRAL/'CMakeFiles/fallthrough_candidate.dir/WirehairV2ProfileFallthrough.cpp.o'
OBJECT_SHA = '473b113cc96ca1e8056c78b473af4c47aa8f5227ce6a7a07e8d97a10ae737bde'
CANDIDATE_SHA = 'be7514b66cad509a42f39b8d4f458ebd24208c0b53d2748aaea7260f3af9786e'
CASES = R.CASES + tuple((4,k,b,p) for k in (3,5) for b in (2,64,1280) for p in (1,2))
SOURCES = D.NEW + ('bench/Wh2CertifiedFallthroughCostR0.cpp',
    'bench/Wh2CertifiedFallthroughCostR0.py','bench/test_Wh2CertifiedFallthroughCostR0.py',
    'bench/Wh2CertifiedFallthroughCostR0.md',
    'bench/Wh2CertifiedFallthrough/CMakeLists.txt','bench/Wh2CertifiedFallthrough/README.md',
    'test/V2SmallCodecTest.cpp','test/V2CertifiedCompatibility.cpp')


def packet_ids(k):
    return tuple(range(k+8))+tuple(0xffffffff-2*j for j in range(6))


@lru_cache(maxsize=256)
def products(a):
    # Carryless polynomial arithmetic, independent of all native lookup tables.
    return bytes(R.O.multiply(a,b) for b in range(256))


@lru_cache(maxsize=2)
def rows(k):
    A.require(k in (3,5),'explicit small dimension')
    def matmul(a,b):
        out = [0]*(k*k)
        for i in range(k):
            for j in range(k):
                for c in range(k): out[k*i+j] ^= products(a[k*i+c])[b[k*c+j]]
        return tuple(out)
    feedback = ((8,14,7),(9,14,7)) if k==3 else ((121,110,207,198,31),(120,110,207,198,31))
    pair = []
    for taps in feedback:
        m = [0]*(k*k)
        for i in range(k-1): m[k*(i+1)+i] = 1
        for i,v in enumerate(taps): m[k*i+k-1] = v
        pair.append(tuple(m))
    levels = [pair]
    for _ in range(31):
        a,b = levels[-1]; levels.append((matmul(a,b),matmul(b,a)))
    result = []
    for packet in packet_ids(k):
        row = (1,)+(0,)*(k-1)
        for bit in range(32):
            if packet & (1<<bit):
                matrix = levels[bit][bin(packet>>(bit+1)).count('1')%2]
                out = [0]*k
                for i in range(k):
                    for j in range(k): out[i] ^= products(matrix[k*i+j])[row[j]]
                row = tuple(out)
        result.append(row)
    return tuple(result)


def first_success(incoming, k):
    basis = {}
    for index,original in enumerate(incoming):
        row = list(original)
        for c in range(k):
            if not row[c]: continue
            if c not in basis:
                basis[c] = row
                if len(basis)==k: return index+1
                break
            pivot = basis[c]; a,b = products(pivot[c]),products(row[c])
            row = [a[x]^b[y] for x,y in zip(row,pivot)]
    raise ValueError('small fixture never reaches full rank')


@lru_cache(maxsize=6)
def small_fixture(k,b):
    source = bytes((37*i+i//11)%256 for i in range(k*b))
    packets = bytearray()
    for row in rows(k):
        for j in range(b):
            v = 0
            for c in range(k): v ^= products(row[c])[source[c*b+j]]
            packets.append(v)
    profile = struct.pack('<4sHHQQII',b'WHV2',1,32,
        0x67c1043ecaa9e184 if k==3 else 0x80070c81bfe375f1,k*b,b,0)
    # Decoder sees eight low repairs, six distant repairs, then systematic.
    r = rows(k)
    steps = first_success(r[k:]+r[:k],k)
    return dict(profile=profile.hex(),packets=packets.hex(),steps=steps)


def verify_header(header,order,claim,meta,protocol=PROTOCOL):
    A.exact(len(header['fixtures']),32,'complete extended fixture roster')
    # Preserve the exact old 20-fixture oracle, including its historical pins.
    R.verify_header(dict(header,fixtures=header['fixtures'][:20]),order,claim,meta,protocol)
    for fixture,c in zip(header['fixtures'][20:],CASES[20:]):
        _,k,b,_ = c
        A.exact(set(fixture),{'case','batch','source','arms'},'public-small fixture schema')
        A.exact(fixture['case'],list(c),'exact public-small route')
        A.exact(fixture['batch'],128,'small lifecycle batch')
        A.exact(fixture['source'],bytes((37*i+i//11)%256 for i in range(k*b)).hex(),'small source')
        A.exact(fixture['arms'],[small_fixture(k,b)]*2,'independent small profile/payload/first-success oracle')


def combine(results,protocol=PROTOCOL):
    result = R.combine(results,protocol)
    improvements = []
    for c in R.CASES:
        if c[0]!=0: continue
        for metric in (0,1):
            cells = [r for load in results for r in load['statistics']
                if r['case']==list(c) and r['metric']==metric and r['comparison']==2]
            A.exact(len(cells),4,'both measurement orders in both load orders')
            A.exact([r['order'] for r in cells],[0,1,0,1],'complete four-cell improvement gate')
            if all(r['estimate']['upper95_log']<0 for r in cells):
                improvements.append(dict(case=list(c),metric=metric))
    result.update(certified_four_cell_improvements=improvements,
        candidate_retained=result['outcome']=='PASS' and bool(improvements),
        pre_admission_restoration_qualified=False)
    return result


def candidate_compile(output,dep):
    return ['/usr/bin/c++','-DNDEBUG','-DWIREHAIR_BUILDING=1',
        '-I'+str(ROOT/'include'),'-I'+str(ROOT/'codec'),'-I'+str(ROOT),
        '-std=gnu++11','-fPIC','-O3','-Wall','-Wextra','-Wpedantic','-Werror',
        '-MD','-MF',str(dep),'-o',str(output),'-c',str(OVERLAY)]


def candidate_inputs():
    source = ROOT/'codec/WirehairV2Profile.cpp'
    A.exact(U.pin(source)['sha256'],'3e3fd6e0dbf1f0d5d12cc2409f6b0ea71b2853d210457f76a1e3f21fc09d68b6','frozen production facade')
    text = source.read_text()
    for old,new,count in (
        ('if (codec && codec->SmallK)','if (codec && CAT_UNLIKELY(codec->SmallK))',1),
        ('if (impl->SmallK)','if (CAT_UNLIKELY(impl->SmallK))',2),
        ('const WirehairV2Result result = impl->SmallK ?',
         'const WirehairV2Result result = CAT_UNLIKELY(impl->SmallK) ?',1)):
        A.exact(text.count(old),count,'unambiguous four-hint overlay'); text=text.replace(old,new)
    A.exact(OVERLAY.read_bytes(),text.encode(),'exact equation-neutral generated overlay')
    A.exact(U.pin(OVERLAY)['sha256'],'c72422d514d8ec565332f803b0985b8db30ef3cc481e5381ab91d81f216e15bd','qualified overlay')
    A.exact(U.pin(OBJECT)['sha256'],OBJECT_SHA,'qualified candidate object')
    for name,digest in (
        ('Testing/Temporary/LastTest.log','1a25fc6a2cffe56093854d8e2dc803912763d0c334ced6e302f792715f5a6739'),
        ('certified-parity.json','fd863d086ebd964af91d098d5318b7e605db40342dcdbe278c753db237351f88'),
        ('link-proof.json','df3e902103b289218835d1fa70e855360f0c6d94bf9c7a4a689ccaa309fec095')):
        A.exact(U.pin(NEUTRAL/name)['sha256'],digest,'prior neutral correctness evidence')
    # Pin preserved objects, binaries, maps and build recipes, never rebuild or
    # perform CTest discovery there. Production inputs are checked by U as well.
    inputs = {p for p in NEUTRAL.rglob('*') if p.is_file()}
    inputs.update(ROOT/n for n in SOURCES)
    for name in ('bench/Wh2CertifiedFallthrough/CMakeLists.txt','test/V2SmallCodecTest.cpp',
                 'test/V2CertifiedCompatibility.cpp'):
        A.exact(A.read_regular(ROOT/name,1024*1024),R.command(['git','cat-file','blob','3af52cb:'+name]),
            'unchanged producing neutral recipe/test source')
    records = A.decode(A.read_regular(NEUTRAL/'compile_commands.json',65536))
    record = [r for r in records if r['output']==str(OBJECT.relative_to(NEUTRAL))]
    A.exact(len(record),1,'one original candidate compile')
    expected = candidate_compile(OBJECT.relative_to(NEUTRAL),'unused')
    del expected[expected.index('-MD'):expected.index('-o')]
    A.exact(shlex.split(record[0]['command']),expected,'qualified candidate compile recipe')
    # Reconstruct the entire compiler input set without writing preserved files.
    args = expected[:expected.index('-o')]+['-M','-MT','candidate',str(OVERLAY)]
    raw = R.command(args)
    dependencies = U.preprocessor_dependencies(raw,'candidate')
    inputs.update(dependencies)
    return inputs,dict(overlay=U.pin(OVERLAY),object=U.pin(OBJECT),
        compile=record[0],preprocessor_command=args,preprocessor_dependencies=raw.decode(),
        neutral_evidence=[U.pin(p) for p in sorted(inputs) if NEUTRAL in p.parents])


def preparation_inputs():
    previous,inputs = D.provenance()
    installed = previous[1]
    candidates,proof = candidate_inputs(); inputs.update(candidates)
    original = Path(installed['original']['path'])
    linked = installed['shared_link'][len(H.native_link(original)):-1]
    objects = [original.parent/p for p in linked]
    A.exact(len(objects),18,'original ordered 18-object DSO')
    selected = [i for i,p in enumerate(objects) if p.name=='WirehairV2Profile.cpp.o']
    A.exact(len(selected),1,'one replaced facade object')
    objects[selected[0]] = OBJECT
    return installed,objects,inputs,proof


def prepare():
    A.require(PREPARED.is_dir() and not PREPARED.is_symlink() and not list(PREPARED.iterdir()),
        'sole empty external preparation directory')
    installed,objects,inputs,proof = preparation_inputs()
    H.relink(dict(installed,proof_name='proof-old.so'),PREPARED)
    obj,dep = PREPARED/'proof-profile.o',PREPARED/'proof-profile.d'
    args = candidate_compile(obj,dep); R.command(args)
    A.exact(U.pin(obj)['sha256'],OBJECT_SHA,'reproduced candidate from exact generated source')
    target = PREPARED/'libwirehair.so.2.0.0'
    link = H.native_link(target)[:-1]+[str(target)]+list(map(str,objects))+['-lm']
    R.command(link)
    A.publish(PREPARED/'preparation.json',A.canonical(dict(
        original=installed,overlay_provenance=proof,compile=args,link=link,
        inputs=[U.pin(p) for p in sorted(inputs)],
        artifacts=[U.pin(p) for p in sorted(PREPARED.iterdir())])))
    print(json.dumps(dict(candidate=U.pin(target),scientific_launch=False)))


def provenance(proof_dir=None):
    installed,objects,inputs,proof = preparation_inputs()
    prepared = A.decode(A.read_regular(PREPARED/'preparation.json',1024*1024))
    A.exact(prepared['original'],installed,'original preparation chain')
    A.exact(prepared['overlay_provenance'],proof,'qualified overlay preparation chain')
    # Harness sources may advance during neutral qualification, never during
    # timing. All preparation-time non-harness inputs must still be identical.
    for pin in prepared['inputs']+prepared['artifacts']:
        path = Path(pin['path'])
        if path not in {ROOT/n for n in SOURCES}:
            A.exact(U.pin(path),pin,'unchanged prepared input/artifact')
        inputs.add(path)
    inputs.add(PREPARED/'preparation.json')
    original = Path(installed['original']['path'])
    target = PREPARED/'libwirehair.so.2.0.0'
    A.exact(U.pin(target)['sha256'],CANDIDATE_SHA,'prospectively pinned candidate DSO')
    link = H.native_link(target)[:-1]+[str(target)]+list(map(str,objects))+['-lm']
    A.exact(prepared['link'],link,'single-object original-order link')
    A.exact(prepared['compile'],candidate_compile(PREPARED/'proof-profile.o',PREPARED/'proof-profile.d'),
        'exact candidate recompile')
    reports = [dict(installed,proof_name='proof-old.so'),dict(source_head=U.SOURCE_HEAD,
        original=U.pin(target),proof_name='proof-new.so',proof_sha256=CANDIDATE_SHA,
        shared_link=H.native_link(target)+list(map(str,objects))+['-lm'],
        original_dso=U.pin(original),overlay_provenance=proof,
        proof_object=dict(name='proof-profile.o',sha256=OBJECT_SHA))]
    if proof_dir is not None:
        for report in reports: H.relink(report,proof_dir)
        obj,dep = proof_dir/'proof-profile.o',proof_dir/'proof-profile.d'
        R.command(candidate_compile(obj,dep))
        A.exact(U.pin(obj)['sha256'],OBJECT_SHA,'fresh qualified candidate object proof')
    return reports,inputs


def qualify(executable,output,meta,mode):
    for name in ('old-new','new-old'):
        A.exact(A.read_regular(output/('neutral-'+name+'.txt'),65536),
            b'PASS neutral82944-coordinate roster, 256 native WORK cases; no timing\n',
            'extended neutral lifecycle roster')
    D.qualify(executable,output,meta,mode,PROTOCOL,verify_header)


def settings():
    A.require(len(CANDIDATE_SHA)==64,'candidate must be pinned before qualification')
    return R.Configuration(PROTOCOL,OUTPUT,(D.LIBRARIES[1],(PREPARED/'libwirehair.so.2.0.0',CANDIDATE_SHA)),
        SOURCES,provenance,'bench/Wh2CertifiedFallthroughCostR0.cpp',qualify,CASES,verify_header,combine)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    sub.add_parser('prepare')
    b = sub.add_parser('build'); b.add_argument('mode',choices=('native','asan-driver')); b.add_argument('output',type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir',type=Path); r.add_argument('output',type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt',type=Path)
    sub.add_parser('replay')
    args = parser.parse_args()
    if args.command=='prepare': prepare(); return
    cfg = settings()
    if args.command=='build': R.build(args.mode,args.output,cfg)
    elif args.command=='receipt': A.publish(args.output,A.canonical(R.receipt(args.build_dir,cfg)))
    elif args.command=='run': R.run(args.receipt,cfg)
    else: print(json.dumps(dict(outcome=R.replay(cfg)['outcome'],exact_replay=True)))


if __name__=='__main__': main()
