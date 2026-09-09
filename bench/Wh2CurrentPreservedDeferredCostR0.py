#!/usr/bin/env python3
"""Current installed preserved paths, fixed cohort and deferred publication."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess

SPEC = importlib.util.spec_from_file_location('current_preserved_shared_gate',
    Path(__file__).with_name('Wh2AdmissionRegressionCostR0.py'))
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)
U = R.sibling('current_preserved_installed_proof','Wh2K5PublicCostR0.py')
H = R.sibling('current_preserved_historical_proof','Wh2SmallIsolationPreservedCostR0.py')

A, ROOT = R.A, R.ROOT
PROTOCOL = 'wirehair.wh2.current-preserved-deferred-cost-r0'
OUTPUT = Path('/var/tmp/wh2-current-preserved-deferred-cost-r0')
LIBRARIES = (R.N.LIBRARIES[0], (U.QUALIFIED/'native-default/libwirehair.so.2.0.0',
    '0a3b1b964ced8a549b16dabfe04a0f35562f5eb268167703abacf0f6f8eacd88'))
NEW = R.NEW + ('bench/Wh2CurrentPreservedDeferredCostR0.cpp',
    'bench/Wh2CurrentPreservedDeferredCostR0.py',
    'bench/test_Wh2CurrentPreservedDeferredCostR0.py',
    'bench/Wh2CurrentPreservedDeferredCostR0.md',
    'bench/Wh2K5PublicCostR0.py','bench/Wh2SmallIsolationPreservedCostR0.py')
PUBLICATIONS = ('success','last-recover','throw-recover','last-clock','last-source')


def provenance(proof_dir=None):
    historical, inputs = R.library_provenance(
        diagnostic_exclusions={H.LOST_DIAGNOSTIC['path']:H.LOST_DIAGNOSTIC})
    archives, current_inputs, archive_proof = U.qualified_inputs('native')
    inputs.update(current_inputs)
    dso, digest = LIBRARIES[1]
    A.exact(U.pin(dso)['sha256'],digest,'actual current installed shared library')
    lines = R.command(['/usr/bin/ninja','-C',dso.parent,'-t','commands',dso.name]).decode().splitlines()
    A.exact(lines[:-1],archive_proof['commands'].splitlines()[:-1],'same installed shared producer compiles')
    framed = shlex.split(lines[-1])
    A.exact((framed[:2],framed[-2:]),([':','&&'],['&&',':']),'shared linker framing')
    shared = framed[2:-2]
    expected = H.native_link(dso)+[str(Path(m['object']['path']).relative_to(dso.parent))
                                 for m in archive_proof['members']]+['-lm']
    A.exact(shared,expected,'exact current 18-object shared link')
    diagnostic = dict(claimed=H.LOST_DIAGNOSTIC,
        observed=U.pin(H.LOST_DIAGNOSTIC['path']),
        reason='Prior CTest discovery overwrote a non-producing historical diagnostic; no old receipt is rebound')
    inputs.update((dso,Path(H.LOST_DIAGNOSTIC['path'])))
    report = dict(source_head=U.SOURCE_HEAD,archive=U.pin(archives[0]),
        installed_archive_provenance=archive_proof,shared_link=shared,
        original=U.pin(dso),proof_name='proof-new.so',proof_sha256=digest,
        excluded_historical_diagnostics=[diagnostic])
    result = [historical[0],report]
    if proof_dir is not None:
        # Reconstruct into fresh files only; neither proof DSO is ever loaded.
        for lib in result:
            H.relink(lib,proof_dir)
    return result,inputs


def verify_publication(raw,mode,order,meta,reference_header):
    A.require(mode in PUBLICATIONS,'explicit neutral publication case')
    A.require(raw.endswith(b'\n') and len(raw)<4*1024*1024,'complete bounded neutral publication')
    rows = [A.decode(line) for line in raw.splitlines()]
    A.exact(len(rows),5,'three retained neutral records with header/footer')
    header,records,footer = rows[0],rows[1:-1],rows[-1]
    R.verify_header(header,order,'0'*64,meta,PROTOCOL)
    # DSO bases legitimately vary between these separate neutral processes.
    # Their complete binding graphs are independently checked above.
    A.exact({k:v for k,v in header.items() if k!='bindings'},
            {k:v for k,v in reference_header.items() if k!='bindings'},'immutable initial header content')
    previous=header['prelude']; work=0
    for i,(row,coordinate) in enumerate(zip(records,list(R.roster())[54:57])):
        A.exact(set(row),{'type','coordinate','ready','target','wait','observation','counts',
                         'addresses','address_count','complete','checked'},'neutral record schema')
        A.exact(row['type'],'record','neutral record type')
        A.exact(row['coordinate'],coordinate,'unchanged decoder coordinates')
        A.exact(coordinate[3:6],[0,1,1],'neutral certified candidate decoder')
        A.exact(row['ready'],previous['clocks'][5],'neutral ready chronology')
        A.exact(row['target'],row['ready'],'neutral omits scientific delay')
        A.exact(row['wait'],[row['ready'],previous['clocks'][4],row['ready'],previous['clocks'][4]],
                'neutral fake-clock wait')
        final = i==2
        A.exact(row['complete'],not(final and mode in ('last-recover','throw-recover')),'retained WORK completion')
        A.exact(row['checked'],not final or mode=='success','retained validation completion')
        steps=header['fixtures'][0]['arms'][1]['steps']
        A.exact(row['counts'],[0,0,128,128*steps,128,128],'whole neutral final-call ledger')
        A.exact(row['address_count'],128,'whole neutral handle roster')
        A.exact(len(row['addresses']),128,'all neutral handles retained')
        for address in row['addresses']: A.integer(address,1)
        if final and mode=='last-clock':
            observed=row['observation']
            A.exact(set(observed),{'clocks','before','after'},'partial clock schema')
            c=observed['clocks']; A.exact(len(c),6,'partial clock shape')
            for value in c: A.integer(value)
            A.require(c[0]>previous['clocks'][5] and c[1]>=previous['clocks'][4] and c[2]>c[0],
                      'partial clock prefix')
            A.exact(c[3:],[0,0,0],'unfinished final clock remains zero')
            A.exact(observed['before'],[0]*4,'neutral initial counters')
            A.exact(observed['after'],[0]*4,'unfinished neutral counters')
        else:
            R.O.clocks(row['observation'],previous)
            c=row['observation']['clocks']; work+=c[3]-c[2]; previous=row['observation']
    A.exact(footer,dict(type='footer',complete=mode=='success',records=3,work_ns=work),'neutral failure footer')
    return dict(mode=mode,load_order=order,records=3,scientific_launch=False)


def qualify(executable,output,meta,mode):
    for order,name in enumerate(('old-new','new-old')):
        reference=A.decode(A.read_regular(output/('fixtures-'+name+'.json'),4*1024*1024))
        for case in PUBLICATIONS:
            raw=R.command([executable,'--neutral-publication',name,case])
            verify_publication(raw,case,order,meta,reference)
            A.publish(output/('publication-'+name+'-'+case+'.jsonl'),raw)
    failures=[]
    with open('/dev/full','wb') as full:
        p=subprocess.run([str(executable),'--neutral-publication','old-new','success'],
            stdin=subprocess.DEVNULL,stdout=full,stderr=subprocess.PIPE,timeout=60)
    A.require(p.returncode==1 and b'output stream' in p.stderr,'full output device is terminal')
    failures.append(dict(sink='full',code=p.returncode,stderr=p.stderr.decode()))
    read_fd,write_fd=os.pipe(); os.close(read_fd)
    try:
        p=subprocess.run([str(executable),'--neutral-publication','old-new','success'],
            stdin=subprocess.DEVNULL,stdout=write_fd,stderr=subprocess.PIPE,timeout=60)
    finally:
        os.close(write_fd)
    A.require(p.returncode==1 and b'output stream' in p.stderr,'broken pipe is a reported terminal output error')
    failures.append(dict(sink='broken-pipe',code=p.returncode,stderr=p.stderr.decode()))
    A.publish(output/'publication-output-errors.json',A.canonical(failures))


SETTINGS = R.Configuration(PROTOCOL,OUTPUT,LIBRARIES,NEW,provenance,
    'bench/Wh2CurrentPreservedDeferredCostR0.cpp',qualify)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build'); b.add_argument('mode',choices=('native','asan-driver')); b.add_argument('output',type=Path)
    r=sub.add_parser('receipt'); r.add_argument('build_dir',type=Path); r.add_argument('output',type=Path)
    r=sub.add_parser('run'); r.add_argument('receipt',type=Path)
    sub.add_parser('replay')
    args=parser.parse_args()
    if args.command=='build': R.build(args.mode,args.output,SETTINGS)
    elif args.command=='receipt': A.publish(args.output,A.canonical(R.receipt(args.build_dir,SETTINGS)))
    elif args.command=='run': R.run(args.receipt,SETTINGS)
    else: print(json.dumps(dict(outcome=R.replay(SETTINGS)['outcome'],exact_replay=True)))


if __name__=='__main__': main()
