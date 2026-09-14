#!/usr/bin/env python3
"""Fresh current-source engineering checks after loss of /tmp build artifacts.

No historical library identity, scientific receipt or timing qualification is
claimed. This retains fresh compiler/test outputs for subsequent independent
producing-closure qualification. It never runs a scientific worker.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
MODES = ('native','scalar','asan')
SANITIZERS = dict(ASAN_OPTIONS='detect_leaks=1:detect_stack_use_after_return=1:halt_on_error=1',
                  UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
TARGETS = ('small_codec_test','v2_small_codec_test','v2_small_k5_codec_test','v2_small_k8_codec_test',
           'small_c_consumer','v2_small_c_consumer','v2_small_k5_c_consumer','v2_small_k8_c_consumer',
           'k6_codec_test','k6_c_consumer','k6_payload_test','v2_borrowed_facade_fault_test',
           'gf256_inplace_test','portability_roundtrip_test')
SHARED = ('small_c_consumer_shared','v2_small_c_consumer_shared','v2_small_k5_c_consumer_shared',
          'v2_small_k8_c_consumer_shared','k6_c_consumer_shared')


def require(ok, why):
    if not ok: raise ValueError(why)


def canonical(value):
    return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()


def pin(path):
    path=Path(path)
    require(path.is_file() and not path.is_symlink(), 'regular evidence file')
    raw=path.read_bytes()
    return dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())


def publish(path, raw):
    with Path(path).open('xb') as stream:
        stream.write(raw); stream.flush(); os.fsync(stream.fileno())
    Path(path).chmod(0o400)


def environment():
    return dict(PATH='/usr/bin:/bin',LANG='C',LC_ALL='C',TZ='UTC',**SANITIZERS)


def build(mode, base):
    require(mode in MODES,'backend')
    base=Path(base)
    require(base.is_absolute() and base.parent==Path('/var/tmp') and
            base.name.startswith('wh2-k4-fresh-neutral-r0.') and base==base.resolve(strict=True) and
            not base.is_symlink(),'fresh durable qualification root')
    output=base/mode
    output.mkdir(mode=0o700)
    commands=[]
    frozen={}

    def stable():
        for path,expected in frozen.items(): require(pin(path)==expected,'source changed: '+str(path))

    def run(argv, cwd=ROOT, timeout=120):
        argv=list(map(str,argv)); index=len(commands)
        invocation=dict(argv=argv,cwd=str(cwd),environment=environment())
        commands.append(invocation); stable()
        publish(output/('command-%03d.argv.json'%index),canonical(invocation))
        try:
            p=subprocess.run(argv,cwd=cwd,env=environment(),stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=timeout)
            stdout,stderr=p.stdout,p.stderr; result=dict(returncode=p.returncode)
        except subprocess.TimeoutExpired as error:
            stdout,stderr=error.stdout or b'',error.stderr or b''; result=dict(timeout=True)
        except OSError as error:
            stdout,stderr=b'',b''; result=dict(error=str(error))
        publish(output/('command-%03d.stdout'%index),stdout)
        publish(output/('command-%03d.stderr'%index),stderr)
        publish(output/('command-%03d.result.json'%index),canonical(result))
        require(result==dict(returncode=0),'fresh neutral command failed: '+str(argv))
        require(len(stdout)<=32*1024**2 and len(stderr)<=1024**2,'neutral output cap')
        stable(); return stdout

    try:
        head=run(['/usr/bin/git','rev-parse','HEAD']).decode().strip()
        names=run(['/usr/bin/git','ls-files','-z']).decode().split('\0')
        # Freeze all repository build/source inputs, including this new runner.
        selected={ROOT/name for name in names if name and
                  (Path(name).suffix in ('.cpp','.c','.h','.hpp','.inc','.cmake','.py') or
                   Path(name).name=='CMakeLists.txt')}
        selected.add(Path(__file__).resolve())
        frozen.update((path,pin(path)) for path in sorted(selected))
        publish(output/'SOURCE.json',canonical(dict(head=head,source_pins=list(frozen.values()))))
        library=output/'library'
        flags={'native':'','scalar':'-DANDROID',
               'asan':'-fsanitize=address,undefined -fno-omit-frame-pointer -march=native'}[mode]
        options=['-G','Ninja','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
                 '-DCMAKE_BUILD_TYPE='+('Debug' if mode=='asan' else 'Release'),
                 '-DCMAKE_C_FLAGS='+flags,'-DCMAKE_CXX_FLAGS='+flags,
                 '-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF']
        run(['/usr/bin/cmake','-S',ROOT,'-B',library]+options+[
            '-DBUILD_TESTS=ON','-DBUILD_CODEC_V2=OFF','-DMARCH_NATIVE=OFF',
            '-DWIREHAIR_STRICT_WARNINGS=ON','-DWH_LTO=OFF','-DWH_PGO_MODE=OFF',
            '-DWIREHAIR_BUILD_BOTH='+('ON' if mode=='native' else 'OFF')])
        targets=list(TARGETS)+(list(SHARED) if mode=='native' else [])
        run(['/usr/bin/ninja','-C',library,'-d','keepdepfile','-j','8','wirehair']+targets,timeout=300)
        tests='^('+'|'.join(targets)+')$'
        test_output=run(['/usr/bin/ctest','--test-dir',library,'--output-on-failure','-j','4','-R',tests],timeout=240)
        require(('100% tests passed, 0 tests failed out of '+str(len(targets))).encode() in test_output,
                'exact selected library test count')
        boundary=output/'k4'
        run(['/usr/bin/cmake','-S',ROOT/'bench/Wh2SmallNative','-B',boundary]+options+[
            '-DWH2_SMALL_TEST_DIMENSION=4','-DWH2_SMALL_LIBRARY='+str(library/'libwirehair.a'),
            '-DWH2_SMALL_PORTABLE='+('ON' if mode=='scalar' else 'OFF'),
            '-DPython3_EXECUTABLE=/usr/bin/python3'])
        run(['/usr/bin/ninja','-C',boundary,'-d','keepdepfile','-j','8'],timeout=240)
        test_output=run(['/usr/bin/ctest','--test-dir',boundary,'--output-on-failure','-j','4'],timeout=240)
        require(b'100% tests passed, 0 tests failed out of 7' in test_output,'all seven K4 boundary tests')
        require(pin(boundary/'Wh2K4NativeData.inc')['sha256']==
                '608bafbe37ac0ba3aa94f5d030390af6cc623f9cf0791e86c5bcb8d72f277763',
                'exact retained fixture bytes, no new selection')
        require(run(['/usr/bin/git','rev-parse','HEAD']).decode().strip()==head,'HEAD changed')
        stable()
        artifacts=[pin(path) for path in sorted(output.rglob('*')) if path.is_file() and not path.is_symlink()]
        result=dict(schema='wh2-k4-fresh-neutral-r0',mode=mode,status='PASS',head=head,
                    library_tests=len(targets),boundary_tests=7,commands=commands,
                    source_pins=list(frozen.values()),artifacts=artifacts,
                    producing_source_closure=False,scientific_launch=False,
                    scope='fresh neutral checks; independent producing/tool/runtime qualification still required')
        publish(output/'RESULT.json',canonical(result))
        print(json.dumps(dict(mode=mode,status='PASS',library_tests=len(targets),boundary_tests=7,
                              result=pin(output/'RESULT.json'),scientific_launch=False)),flush=True)
    except Exception as error:
        publish(output/'FAILED.json',canonical(dict(mode=mode,status='FAILED',failure=str(error),
                    commands=commands,scientific_launch=False)))
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=MODES); parser.add_argument('output',type=Path)
    args=parser.parse_args(); build(args.mode,args.output)
