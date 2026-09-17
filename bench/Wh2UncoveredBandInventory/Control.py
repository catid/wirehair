"""One-shot baseline diagnostic with bounded neutral qualification and capture."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

import Inventory as I
import Launch
import Support as A

HERE = Path(__file__).resolve().parent
ROOT = I.ROOT
OUTPUT = Path('/var/tmp/wh2-uncovered-band-inventory-r0')
QUALIFICATION = Path('/var/tmp/wh2-uncovered-band-qualification.C0lOL3VL')
CODEC_MANIFEST = I.QUALIFIED/'NEUTRAL.sha256'
CODEC_SHA = '93bfff5491f451c15ff8be50ec4453bef8770c0dbbc24e290ad48caa1822f532'
OWN = ('Inventory.py','Support.py','Control.py','Launch.py','test_Inventory.py',
       'test_Support.py','test_Launch.py','test_Control.py','README.md')
RUNTIME = ('Inventory.py','Support.py','Control.py','Launch.py')
ENV_KEYS = ('LD_PRELOAD','LD_LIBRARY_PATH','LD_AUDIT','GLIBC_TUNABLES',
            'MALLOC_PERTURB_','MALLOC_TRIM_THRESHOLD_','MALLOC_MMAP_THRESHOLD_',
            'MALLOC_TOP_PAD_','ASAN_OPTIONS','UBSAN_OPTIONS','PYTHONPATH')


def git(*arguments):
    return subprocess.check_output(['git']+list(arguments),cwd=ROOT,timeout=30)


def head():
    return git('rev-parse','HEAD').decode().strip()


def environment():
    A.require(not any(key in os.environ for key in ENV_KEYS),'clean diagnostic environment')


def codec_pins():
    raw = A.read_regular(CODEC_MANIFEST,1024*1024)
    A.exact(A.sha(raw),CODEC_SHA,'retained baseline engineering manifest')
    pins = {str(CODEC_MANIFEST):CODEC_SHA}
    for line in raw.decode().splitlines():
        sha, separator, name = line.partition('  ')
        A.require(separator and len(sha)==64 and set(sha)<=set('0123456789abcdef') and
                  Path(name).is_absolute(),'canonical absolute codec pin')
        path = Path(name).resolve(strict=True)
        A.require(str(path) not in pins,'unique codec pin')
        A.exact(A.pin(path)['sha256'],sha,'unchanged qualified artifact')
        pins[str(path)] = sha
    A.exact(len(pins),1040,'all 1039 engineering members plus manifest')
    return pins


def neutral_inputs(mode, interpreter=None):
    A.require(mode in I.LIBRARIES,'neutral backend')
    environment()
    return dict(protocol=I.PROTOCOL,mode=mode,
                sources={str(HERE/name):A.pin(HERE/name)['sha256'] for name in RUNTIME},
                interpreter=A.pin(sys.executable if interpreter is None else interpreter),
                library=A.pin(I.LIBRARIES[mode][0]),limits=Launch.limits(),
                environment=Launch.ENVIRONMENT)


def neutral_paths(mode):
    A.require(mode in I.LIBRARIES,'neutral backend')
    return {name:QUALIFICATION/(mode+suffix) for name,suffix in
            (('claim','.claim.json'),('raw','.jsonl'),('stderr','.stderr.txt'),
             ('process','.process.json'),('proof','.json'))}


def neutral_worker(mode, claim_sha):
    paths = neutral_paths(mode)
    raw = A.read_regular(paths['claim'],65536)
    A.exact(A.sha(raw),claim_sha,'actual neutral claim')
    claim = A.decode(raw)
    A.exact(claim,neutral_inputs(mode),'unchanged neutral producer and interpreter')
    I.run_inventory(mode,I.neutral_roster(),claim_sha,
                    lambda record: sys.stdout.buffer.write(A.canonical(record)))
    sys.stdout.buffer.flush()
    A.exact(claim,neutral_inputs(mode),'unchanged neutral inputs after execution')


def neutral_evidence(mode):
    paths = neutral_paths(mode)
    raw_claim = A.read_regular(paths['claim'],65536)
    claim = A.decode(raw_claim)
    A.exact(claim,neutral_inputs(mode,claim['interpreter']['path']),
            'neutral evidence produced by the current runtime sources')
    process = A.decode(A.read_regular(paths['process'],65536))
    check_capture(process,paths['raw'],paths['stderr'])
    result = I.verify(paths['raw'],A.sha(raw_claim),mode,True)
    return dict(protocol=I.PROTOCOL,mode=mode,result=result,
                artifacts={name:A.pin(path) for name,path in paths.items() if name!='proof'})


def neutral(mode):
    """Bounded systematic check on different widths; not a scientific cohort."""
    paths = neutral_paths(mode)
    A.require(not any(path.exists() or path.is_symlink() for path in paths.values()),
              'fresh neutral outputs')
    codec_pins()
    claim = neutral_inputs(mode)
    A.publish(paths['claim'],claim)
    try:
        result = Launch.capture([claim['interpreter']['path'],'-B',str(Path(__file__).resolve()),
                                 'neutral-worker',mode,A.sha(paths['claim'].read_bytes())],
                                paths['raw'],paths['stderr'])
        A.publish(paths['process'],result)
        A.exact(claim,neutral_inputs(mode),'unchanged neutral producer after capture')
        proof = neutral_evidence(mode)
        A.publish(paths['proof'],proof)
        print(mode,proof['result'],flush=True)
    finally:
        for path in paths.values():
            if path.exists(): path.chmod(0o400)


def inputs(interpreter=None):
    environment()
    pins = codec_pins()
    parity, producers = [], []
    for mode in I.LIBRARIES:
        paths = neutral_paths(mode)
        proof = A.decode(A.read_regular(paths['proof'],65536))
        expected = neutral_evidence(mode)
        A.exact(proof,expected,'exact complete backend neutral evidence')
        parity.append(proof['result']['parity_sha256'])
        for path in paths.values():
            pins[str(path)] = A.pin(path)['sha256']
        neutral_claim = A.decode(A.read_regular(paths['claim'],65536))
        for name,pin in neutral_claim['sources'].items(): pins[name] = pin
        pin = neutral_claim['interpreter']
        producers.append(pin)
        pins[pin['path']] = pin['sha256']
    A.exact(parity[0],parity[1],'all native/portable neutral records agree')
    for name in OWN:
        path = HERE/name
        raw = path.read_bytes()
        A.exact(git('show','HEAD:'+str(path.relative_to(ROOT))),raw,'commit source before diagnostic')
        pins[str(path)] = A.sha(raw)
    for name in ('python312.tests.log','python38.tests.log'):
        path = QUALIFICATION/name
        pins[str(path)] = A.pin(path)['sha256']
    interpreter = Path(sys.executable if interpreter is None else interpreter).resolve(strict=True)
    producer = A.pin(interpreter)
    for qualified in producers:
        A.exact(producer,qualified,'scientific producer passed both neutral backends')
    pins[str(interpreter)] = producer['sha256']
    return dict(protocol=I.PROTOCOL,head=head(),pins=pins,interpreter=str(interpreter),
                library_mode='native',limits=Launch.limits(),environment=Launch.ENVIRONMENT,
                scope='baseline recovery diagnostic only; not speed, holdout or complete toolchain provenance')


def current(claim):
    # Reanalysis may use another Python version while preserving the actual
    # producing interpreter identity instead of silently substituting it.
    A.exact(claim,inputs(claim['interpreter']),'exact current diagnostic inputs')


def worker(claim_sha):
    A.require(len(claim_sha)==64 and set(claim_sha)<=set('0123456789abcdef'),'claim hex')
    raw = A.read_regular(OUTPUT/'claim.json',1024*1024)
    A.exact(A.sha(raw),claim_sha,'actual inventory claim')
    claim = A.decode(raw)
    A.exact(str(Path(sys.executable).resolve(strict=True)),claim['interpreter'],
            'worker uses the claimed producing interpreter')
    current(claim)
    I.run_inventory('native',I.roster(),claim_sha,
                    lambda record: sys.stdout.buffer.write(A.canonical(record)))
    sys.stdout.buffer.flush()
    current(claim)


def check_capture(result,raw,stderr):
    A.exact(set(result),{'exit','failure','wall_seconds','stdout_bytes','stderr_bytes'},
            'complete process capture schema')
    A.require(type(result['exit']) is int and result['exit']==0 and result['failure'] is None,
              'successful bounded worker')
    A.require(type(result['wall_seconds']) in (int,float) and
              0 < result['wall_seconds'] <= Launch.WALL_SECONDS,'worker wall cap')
    A.exact(result['stdout_bytes'],raw.stat().st_size,'complete raw capture')
    A.require(0 < result['stdout_bytes'] <= I.RAW_CAP,'raw size cap')
    A.exact(result['stderr_bytes'],0,'empty captured stderr')
    A.exact(stderr.stat().st_size,0,'empty retained stderr')


def run():
    A.require(not OUTPUT.exists() and not OUTPUT.is_symlink(),'scientific namespace already spent')
    claim = inputs()
    OUTPUT.mkdir(mode=0o700)
    try:
        A.publish(OUTPUT/'claim.json',claim)
        claim_sha = A.sha((OUTPUT/'claim.json').read_bytes())
        result = Launch.capture([claim['interpreter'],'-B',str(Path(__file__).resolve()),
                                 'worker',claim_sha],OUTPUT/'raw.jsonl',OUTPUT/'stderr.txt')
        A.publish(OUTPUT/'process.json',result)
        check_capture(result,OUTPUT/'raw.jsonl',OUTPUT/'stderr.txt')
        current(claim)
        analysis = I.verify(OUTPUT/'raw.jsonl',claim_sha)
        A.publish(OUTPUT/'analysis.json',analysis)
        A.publish(OUTPUT/'complete.json',{path.name:A.pin(path)['sha256'] for path in OUTPUT.iterdir()})
        print(analysis['outcome'],'priority',analysis['priority_order'],
              'recommended',analysis['recommended_k'],flush=True)
    except BaseException as error:
        A.publish(OUTPUT/'failed.json',dict(error=type(error).__name__+': '+str(error),
                                          outcome='INVALID',namespace_spent=True))
        raise
    finally:
        for path in OUTPUT.iterdir():
            path.chmod(0o400)


def replay():
    names = {'claim.json','raw.jsonl','stderr.txt','process.json','analysis.json'}
    complete = A.decode(A.read_regular(OUTPUT/'complete.json',65536))
    A.exact(set(complete),names,'exact complete inventory')
    A.exact({path.name for path in OUTPUT.iterdir()},names|{'complete.json'},'exact bundle members')
    for name,sha in complete.items():
        A.exact(A.pin(OUTPUT/name)['sha256'],sha,'sealed retained member')
    claim_raw = A.read_regular(OUTPUT/'claim.json',1024*1024)
    claim = A.decode(claim_raw)
    current(claim)
    check_capture(A.decode(A.read_regular(OUTPUT/'process.json',65536)),OUTPUT/'raw.jsonl',OUTPUT/'stderr.txt')
    analysis = I.verify(OUTPUT/'raw.jsonl',A.sha(claim_raw))
    A.exact(analysis,A.decode(A.read_regular(OUTPUT/'analysis.json',1024*1024)),'complete retained reanalysis')
    print(analysis['outcome'],'priority',analysis['priority_order'],
          'recommended',analysis['recommended_k'],flush=True)
    return analysis


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    n = sub.add_parser('neutral'); n.add_argument('mode',choices=tuple(I.LIBRARIES))
    n = sub.add_parser('neutral-worker'); n.add_argument('mode',choices=tuple(I.LIBRARIES)); n.add_argument('claim')
    w = sub.add_parser('worker'); w.add_argument('claim')
    sub.add_parser('run'); sub.add_parser('replay'); sub.add_parser('preflight')
    args = parser.parse_args()
    if args.command=='neutral': neutral(args.mode)
    elif args.command=='neutral-worker': neutral_worker(args.mode,args.claim)
    elif args.command=='worker': worker(args.claim)
    elif args.command=='run': run()
    elif args.command=='replay': replay()
    else: print('pins',len(inputs()['pins']))
