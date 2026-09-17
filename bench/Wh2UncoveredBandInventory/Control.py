"""Disabled draft controller; qualification is incomplete (see README.md)."""
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
       'test_Support.py','test_Launch.py','test_Draft.py','README.md')
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


def neutral(mode):
    """Different widths/systematic prefixes; no scientific namespace or timing."""
    I.require_qualified()
    A.require(mode in I.LIBRARIES,'neutral backend')
    environment()
    raw_path, proof_path = (QUALIFICATION/(mode+suffix) for suffix in ('.jsonl','.json'))
    A.require(not raw_path.exists() and not proof_path.exists(),'fresh neutral outputs')
    with raw_path.open('xb') as stream:
        I.run_inventory(mode,I.neutral_roster(),'neutral',
                        lambda record: stream.write(A.canonical(record)))
    result = I.verify(raw_path,'neutral',mode,True)
    raw_path.chmod(0o400)
    A.publish(proof_path,dict(protocol=I.PROTOCOL,mode=mode,raw=A.pin(raw_path),result=result,
                             library=A.pin(I.LIBRARIES[mode][0])))
    print(mode,result,flush=True)


def inputs(interpreter=None):
    I.require_qualified()
    environment()
    pins = codec_pins()
    for mode in I.LIBRARIES:
        proof_path = QUALIFICATION/(mode+'.json')
        proof = A.decode(A.read_regular(proof_path,65536))
        raw_path = QUALIFICATION/(mode+'.jsonl')
        expected = dict(protocol=I.PROTOCOL,mode=mode,raw=A.pin(raw_path),
                        result=I.verify(raw_path,'neutral',mode,True),
                        library=A.pin(I.LIBRARIES[mode][0]))
        A.exact(proof,expected,'exact complete backend neutral evidence')
        for path in (proof_path,raw_path):
            pins[str(path)] = A.pin(path)['sha256']
    for name in OWN:
        path = HERE/name
        raw = path.read_bytes()
        A.exact(git('show','HEAD:'+str(path.relative_to(ROOT))),raw,'commit source before diagnostic')
        pins[str(path)] = A.sha(raw)
    interpreter = Path(sys.executable if interpreter is None else interpreter).resolve(strict=True)
    pins[str(interpreter)] = A.pin(interpreter)['sha256']
    return dict(protocol=I.PROTOCOL,head=head(),pins=pins,interpreter=str(interpreter),
                library_mode='native',limits=Launch.limits(),environment=Launch.ENVIRONMENT,
                scope='baseline recovery diagnostic only; not speed, holdout or complete toolchain provenance')


def current(claim):
    # Reanalysis may use another Python version while preserving the actual
    # producing interpreter identity instead of silently substituting it.
    A.exact(claim,inputs(claim['interpreter']),'exact current diagnostic inputs')


def worker(claim_sha):
    I.require_qualified()
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


def check_capture(result,root):
    A.require(type(result['exit']) is int and result['exit']==0 and result['failure'] is None,
              'successful bounded worker')
    A.require(type(result['wall_seconds']) in (int,float) and
              0 < result['wall_seconds'] <= Launch.WALL_SECONDS,'worker wall cap')
    A.exact(result['stdout_bytes'],(root/'raw.jsonl').stat().st_size,'complete raw capture')
    A.require(0 < result['stdout_bytes'] <= I.RAW_CAP,'raw size cap')
    A.exact(result['stderr_bytes'],0,'empty captured stderr')
    A.exact((root/'stderr.txt').stat().st_size,0,'empty retained stderr')


def run():
    I.require_qualified()
    A.require(not OUTPUT.exists(),'scientific namespace already spent')
    claim = inputs()
    OUTPUT.mkdir(mode=0o700)
    A.publish(OUTPUT/'claim.json',claim)
    try:
        claim_sha = A.sha((OUTPUT/'claim.json').read_bytes())
        result = Launch.capture([claim['interpreter'],'-B',str(Path(__file__).resolve()),
                                 'worker',claim_sha],OUTPUT/'raw.jsonl',OUTPUT/'stderr.txt')
        A.publish(OUTPUT/'process.json',result)
        check_capture(result,OUTPUT)
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
    I.require_qualified()
    names = {'claim.json','raw.jsonl','stderr.txt','process.json','analysis.json'}
    complete = A.decode(A.read_regular(OUTPUT/'complete.json',65536))
    A.exact(set(complete),names,'exact complete inventory')
    A.exact({path.name for path in OUTPUT.iterdir()},names|{'complete.json'},'exact bundle members')
    for name,sha in complete.items():
        A.exact(A.pin(OUTPUT/name)['sha256'],sha,'sealed retained member')
    claim_raw = A.read_regular(OUTPUT/'claim.json',1024*1024)
    claim = A.decode(claim_raw)
    current(claim)
    check_capture(A.decode(A.read_regular(OUTPUT/'process.json',65536)),OUTPUT)
    analysis = I.verify(OUTPUT/'raw.jsonl',A.sha(claim_raw))
    A.exact(analysis,A.decode(A.read_regular(OUTPUT/'analysis.json',1024*1024)),'complete retained reanalysis')
    print(analysis['outcome'],'priority',analysis['priority_order'],
          'recommended',analysis['recommended_k'],flush=True)
    return analysis


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    n = sub.add_parser('neutral'); n.add_argument('mode',choices=tuple(I.LIBRARIES))
    w = sub.add_parser('worker'); w.add_argument('claim')
    sub.add_parser('run'); sub.add_parser('replay'); sub.add_parser('preflight')
    args = parser.parse_args()
    if args.command=='neutral': neutral(args.mode)
    elif args.command=='worker': worker(args.claim)
    elif args.command=='run': run()
    elif args.command=='replay': replay()
    else: print('pins',len(inputs()['pins']))
