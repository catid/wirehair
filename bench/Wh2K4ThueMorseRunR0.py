#!/usr/bin/env python3
"""K4 one-shot controller and read-only sealed-inventory projection."""
import importlib.util
import io
from pathlib import Path
import sys

SPEC=importlib.util.spec_from_file_location('k4_tm_controller',
    Path(__file__).with_name('Wh2NoncommutingRadixRunR0.py'))
C=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(C)
C.PROTOCOL='wirehair.wh2.k4-thue-morse-r0'
C.OUTPUT=Path('/var/tmp/wh2-k4-thue-morse-r0')
C.SOURCES=('bench/Wh2K4ThueMorseR0.py','bench/Wh2K4ThueMorseRunR0.py',
    'bench/test_Wh2K4ThueMorseR0.py','bench/Wh2K4ThueMorseR0.md',
    'bench/Wh2K3ThueMorseR0.py','bench/Wh2K5ThueMorseR0.py',
    'bench/Wh2K5ThueMorseRunR0.py','bench/Wh2NoncommutingRadixR0.py',
    'bench/Wh2NoncommutingRadixRunR0.py','bench/test_Wh2NoncommutingRadixRunR0.py',
    'bench/Wh2ThueMorseRecoveryHistoryR0.py','bench/Wh2FrozenTrace.cpp',
    'bench/Wh2FrozenTrace.h','bench/wh2_benchmark_contract_v4.json')
INVENTORY=Path('/var/tmp/wh2-uncovered-k-recovery-inventory-r0')
MANIFEST='93a9068986dc24e1caa2931f5be897022e34f13b8f5ca1170a294d87365c5f2c'
PROJECTIONS={
    'origins':(38,'2efbce06e3db3928635764c51c74dd0f79966b7050b4751f592c60f14d75b578'),
    'prefixes':(37,'96cb80f4458adca9f5cc6ed773f9cca27ba9e0b5a048dc0cd8c9f1ca11be42ff'),
    'roots':(64,'6a03f629cae261c18128ba4707b2656dda7833a9f836948442f9fb92e8c20bbf')}


def require(value,message):
    if not value: raise ValueError(message)


def project(raw,deadline=None):
    origins=[]; roots=set(); cases=0; expected_ordinal=0
    for line in io.BytesIO(raw):
        C.time_left(deadline); require(len(line)<=8*1024**2,'historical line cap')
        row=C.strict_json(line)
        if row.get('type')!='case': continue
        require(row['ordinal']==expected_ordinal,'historical chronology'); expected_ordinal+=1
        case=row['case']
        if case['root']: roots.add(case['root'])
        if case['k']!=4: continue
        cases+=1
        require(len(case['ids'])==8 and len(set(case['ids']))==8,'original K4 horizon')
        for arm,result in enumerate(row['arms']):
            for oh in range(5):
                if not result['first'] or result['first']>4+oh:
                    origins.append(dict(ordinal=row['ordinal'],arm=arm,b=case['b'],tail=case['tail'],
                        group=case['group'],root=case['root'],schedule=case['schedule'],
                        overhead=oh,ids=case['ids'][:4+oh]))
    require(cases==774 and expected_ordinal==3096,'entire historical K4 roster')
    widths={}
    for o in origins: widths.setdefault(tuple(o['ids']),set()).add(o['b'])
    result=dict(origins=origins,prefixes=[dict(ids=list(ids),original_widths=sorted(bs))
        for ids,bs in sorted(widths.items())],roots=sorted(roots))
    for key,(count,digest) in PROJECTIONS.items():
        require(len(result[key])==count and C.sha(C.canonical(result[key]))==digest,'historical projection '+key)
    return result


def inventory_inputs(deadline=None):
    names={'CLAIM.json','raw.jsonl','stderr.txt','process.json','analysis.json','COMPLETE.json'}
    require({p.name for p in INVENTORY.iterdir()}==names,'exact historical roster')
    raw=C.read_regular(INVENTORY/'COMPLETE.json',65536,deadline=deadline)
    require(C.sha(raw)==MANIFEST,'historical manifest SHA')
    manifest=C.strict_json(raw)
    require(manifest['outcome']=='DIAGNOSTIC_COMPLETE' and len(manifest['files'])==5 and
        {Path(p['path']).name for p in manifest['files']}==names-{'COMPLETE.json'},'historical manifest members')
    total=len(raw); corpus=None
    for p in manifest['files']:
        path=Path(p['path'])
        require(path.parent==INVENTORY and path.stat().st_mode&0o777==0o400,'sealed historical member')
        data=C.read_regular(path,64*1024**2-total,deadline=deadline); total+=len(data)
        require(len(data)==p['bytes'] and C.sha(data)==p['sha256'],'historical member content')
        if path.name=='raw.jsonl': corpus=data
    require(corpus is not None,'complete original corpus')
    result=project(corpus,deadline)
    result['provenance']=dict(path=str(INVENTORY),manifest_sha256=MANIFEST,files=manifest['files'])
    return result


base_receipt=C.current_receipt


def current_receipt(deadline=None):
    receipt=base_receipt(deadline)
    history=inventory_inputs(deadline)
    receipt['inventory']=dict(provenance=history['provenance'],projection_sha256=C.sha(C.canonical(history)))
    return receipt


C.current_receipt=current_receipt


def claimed_inputs():
    raw=C.read_regular(C.OUTPUT/'CLAIM.json',1024*1024)
    claim=C.strict_json(raw)
    require(raw==C.canonical(claim) and set(claim)=={'protocol','receipt_sha256','receipt'},'canonical claim')
    require(claim['protocol']==C.PROTOCOL and claim['receipt_sha256']==C.sha(C.canonical(claim['receipt'])),
        'claim identity')
    require(claim['receipt']==current_receipt(),'all current claimed inputs')
    return C.sha(raw)


if __name__=='__main__':
    try: sys.exit(C.main())
    except Exception as error:
        print(type(error).__name__+': '+str(error)[:1000],file=sys.stderr); sys.exit(1)
