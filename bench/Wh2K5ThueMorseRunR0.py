#!/usr/bin/env python3
"""K5 one-shot launcher and read-only authenticated inventory projection."""
import importlib.util
import io
from pathlib import Path
import sys


SPEC = importlib.util.spec_from_file_location('k5_capture',
    Path(__file__).with_name('Wh2NoncommutingRadixRunR0.py'))
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
C.PROTOCOL = 'wirehair.wh2.k5-thue-morse-r0'
C.OUTPUT = Path('/var/tmp/wh2-k5-thue-morse-r0')
C.SOURCES = ('bench/Wh2K5ThueMorseR0.py', 'bench/Wh2K5ThueMorseRunR0.py',
    'bench/test_Wh2K5ThueMorseR0.py', 'bench/Wh2K3ThueMorseR0.py',
    'bench/Wh2NoncommutingRadixR0.py', 'bench/Wh2NoncommutingRadixRunR0.py',
    'bench/test_Wh2NoncommutingRadixRunR0.py', 'bench/Wh2ThueMorseRecoveryHistoryR0.py',
    'bench/Wh2FrozenTrace.cpp', 'bench/Wh2FrozenTrace.h', 'bench/wh2_benchmark_contract_v4.json')
INVENTORY = Path('/var/tmp/wh2-uncovered-k-recovery-inventory-r0')
MANIFEST = '93a9068986dc24e1caa2931f5be897022e34f13b8f5ca1170a294d87365c5f2c'
PROJECTIONS = {
    'origins': (79, '72d639fbd329c57faab90cec0c85877c816d7be77c53fc3565a0acbf404240f9'),
    'prefixes': (54, 'c22b92bf2a6271c4e5264a9d8c97d1e8081e5c705c757434a2e0ccc7aa4ea168'),
    'roots': (64, '6a03f629cae261c18128ba4707b2656dda7833a9f836948442f9fb92e8c20bbf')}


def require(value, message):
    if not value: raise ValueError(message)


def project(raw, deadline=None):
    origins, roots = [], set()
    for line in io.BytesIO(raw):
        C.time_left(deadline)
        require(len(line)<=8*1024**2, 'historical line bound')
        row = C.strict_json(line)
        if row.get('type')!='case': continue
        case = row['case']
        if case['root']: roots.add(case['root'])
        if case['k']!=5: continue
        for arm, result in enumerate(row['arms']):
            for overhead in range(5):
                if not result['first'] or result['first']>5+overhead:
                    origins.append(dict(ordinal=row['ordinal'],arm=arm,b=case['b'],
                        overhead=overhead,ids=case['ids'][:5+overhead]))
    result = dict(origins=origins, prefixes=[list(p) for p in sorted({tuple(o['ids']) for o in origins})],
                  roots=sorted(roots))
    for name, (count, digest) in PROJECTIONS.items():
        require(len(result[name])==count and C.sha(C.canonical(result[name]))==digest,
                'historical projection: '+name)
    return result


def inventory_inputs(deadline=None):
    names = {'CLAIM.json','raw.jsonl','stderr.txt','process.json','analysis.json','COMPLETE.json'}
    require({p.name for p in INVENTORY.iterdir()}==names,'historical member roster')
    manifest = C.read_regular(INVENTORY/'COMPLETE.json',65536,deadline=deadline)
    require(C.sha(manifest)==MANIFEST,'historical manifest pin')
    record = C.strict_json(manifest)
    require(record['outcome']=='DIAGNOSTIC_COMPLETE','historical outcome')
    require(len(record['files'])==5 and {Path(p['path']).name for p in record['files']}==names-{'COMPLETE.json'},
            'historical unique members')
    total = len(manifest); data = None
    for p in record['files']:
        path = Path(p['path'])
        require(path.parent==INVENTORY and path.stat().st_mode&0o777==0o400,'historical member path/mode')
        raw = C.read_regular(path,64*1024**2-total,deadline=deadline); total+=len(raw)
        require(len(raw)==p['bytes'] and C.sha(raw)==p['sha256'],'historical member pin')
        if path.name=='raw.jsonl': data=raw
    require(data is not None,'historical raw stream')
    result = project(data,deadline)
    result['provenance'] = dict(path=str(INVENTORY),manifest_sha256=MANIFEST,
        manifest_bytes=len(manifest),files=record['files'])
    return result


_base_receipt = C.current_receipt


def current_receipt(deadline=None):
    receipt = _base_receipt(deadline)
    history = inventory_inputs(deadline)
    receipt['inventory'] = dict(provenance=history['provenance'],projection_sha256=C.sha(C.canonical(history)))
    return receipt


C.current_receipt = current_receipt


def claimed_inputs():
    raw = C.read_regular(C.OUTPUT/'CLAIM.json',1024*1024)
    claim = C.strict_json(raw)
    require(raw==C.canonical(claim) and set(claim)=={'protocol','receipt_sha256','receipt'},'canonical claim')
    require(claim['protocol']==C.PROTOCOL and claim['receipt_sha256']==C.sha(C.canonical(claim['receipt'])),
            'claim identity')
    require(claim['receipt']==current_receipt(),'current authenticated worker inputs')
    return C.sha(raw)


if __name__=='__main__':
    try: sys.exit(C.main())
    except Exception as error:
        print(type(error).__name__+': '+str(error)[:1000],file=sys.stderr)
        sys.exit(1)
