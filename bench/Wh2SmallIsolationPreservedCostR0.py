#!/usr/bin/env python3
"""Frozen pre-admission versus small-state-isolation preserved-path cost gate.

Reuse the original 20-case native WORK and complete observer/statistical gates.
No old namespace, sample, protocol identity or library binding is overwritten.
Ordinary K3 retention/WH1 speed is a separate gate, not a claim of this screen.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shlex

SPEC = importlib.util.spec_from_file_location('small_isolation_shared_gate',
    Path(__file__).with_name('Wh2AdmissionRegressionCostR0.py'))
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)
A, O, ROOT = R.A, R.O, R.ROOT
PROTOCOL = 'wirehair.wh2.small-isolation-preserved-cost-r0'
OUTPUT = Path('/var/tmp/wh2-small-isolation-preserved-cost-r0')
SOURCE_HEAD = '5e902dfad74152af55fa3a404608e9a6d57f53d8'
BASE = Path('/tmp/wh2-small-isolation.b24lX5/native')
ARCHIVE_SHA = '21fa78cc3cccbacd01be13bdd66eb383d54f479f7fdc7de26bb98487380991ac'
LIBRARIES = (R.N.LIBRARIES[0], (BASE/'libwirehair.so.2.0.0',
    '6772f3818ba040ba1d939852e285e3717dc1710e9f9e3b880a929a194a82a71a'))
PREFIX = 'CMakeFiles/wirehair_objects.dir/'
NEW = R.NEW + ('bench/Wh2SmallIsolationPreservedCostR0.py',
               'bench/test_Wh2SmallIsolationPreservedCostR0.py')
# CTest -N unexpectedly overwrote this non-producing historical diagnostic.
# Never claim the old receipt remains wholly intact or alter its original pins.
# The original R0 entrypoint still rejects the mismatch. This new gate records
# the loss and independently proves every code-producing input and artifact.
LOST_DIAGNOSTIC = dict(path='/tmp/wh2-v2-k3-admission.TECbfU/native/Testing/Temporary/LastTest.log',
    bytes=756592,sha256='07b574cbe29cdbf8d1d98f8dc19a0b71a7350490c3c102e70c42916060ec2b54')


def dependency_roster(recorded):
    result = {}
    for block in recorded.strip().split('\n\n'):
        lines = block.splitlines()
        if not lines or not lines[0].startswith(PREFIX):
            continue
        match = re.fullmatch(r'(.+): #deps ([0-9]+), deps mtime [0-9]+ \(VALID\)',lines[0])
        A.require(match is not None and match[1] not in result,'valid unique candidate dependency record')
        paths = [str(Path(line.strip()).resolve(strict=True)) for line in lines[1:]]
        A.exact(len(paths),int(match[2]),'complete candidate compiler dependencies')
        # GCC may record a header twice (gf256.h does so here). Preserve the
        # emitted count and sequence rather than rejecting valid dependencies.
        result[match[1]] = paths
    return result


def native_compile(relative):
    return ['/usr/bin/c++','-DWIREHAIR_BUILDING=1','-I'+str(ROOT/'include'),
            '-O3','-DNDEBUG','-std=gnu++11','-fPIC','-Wall','-Wextra','-Wpedantic','-Werror',
            '-MD','-MT',relative,'-MF',relative+'.d','-o',relative,'-c',str(ROOT/relative[len(PREFIX):-2])]


def native_link(dso):
    return ['/usr/bin/c++','-fPIC','-O3','-DNDEBUG','-Wl,--version-script='+str(ROOT/'abi/wirehair.map'),
            '-shared','-Wl,-soname,libwirehair.so.2','-o',dso.name]


def relink(report, output):
    target = output/report['proof_name']
    A.require(not target.exists() and not target.is_symlink(),'fresh proof DSO')
    original = Path(report['original']['path'])
    prefix = native_link(original)
    linked = report['shared_link'][len(prefix):-1]
    args = prefix[:-1]+[str(target)]+[str(original.parent/p) for p in linked]+['-lm']
    R.command(args)
    A.exact(A.read_regular(target,4*1024**2),A.read_regular(original,4*1024**2),
            'exact untouched candidate/baseline DSO relink')
    target.chmod(0o400)


def provenance(proof_dir=None):
    # The historical source/object chains remain tied to their original Git
    # snapshots, not to current production files. This function never runs codecs.
    historical, inputs = R.library_provenance(diagnostic_exclusions={LOST_DIAGNOSTIC['path']:LOST_DIAGNOSTIC})
    previous = historical[1]
    previous_objects = {Path(p['path']).name:p for p in previous['objects']}
    previous_sources = {p['path']:p for p in previous['historical_source_blobs']}
    archive = BASE/'libwirehair.a'
    A.exact(O.pin(archive)['sha256'],ARCHIVE_SHA,'exact isolation archive')
    dso,digest = LIBRARIES[1]
    A.exact(O.pin(dso)['sha256'],digest,'exact isolation native DSO')
    objects = sorted((BASE/PREFIX).rglob('*.o'))
    members = R.command(['/usr/bin/ar','t',archive]).decode().splitlines()
    A.exact(len(objects),17,'candidate object count')
    A.exact(len(members),17,'candidate archive count')
    A.exact(set(members),{p.name for p in objects},'complete candidate object/archive roster')
    deps = dependency_roster(R.command(['/usr/bin/ninja','-C',BASE,'-t','deps']).decode())
    A.exact(set(deps),{str(p.relative_to(BASE)) for p in objects},'all producing dependency records')
    commands = R.command(['/usr/bin/ninja','-C',BASE,'-t','commands',dso.name]).decode().splitlines()
    A.exact(len(commands),18,'17 candidate compiles and one link')
    compiles = [shlex.split(line) for line in commands[:-1]]
    A.exact(sorted(compiles),sorted(native_compile(p) for p in deps),'exact native producing compiles')
    shared = shlex.split(commands[-1])
    A.exact((shared[:2],shared[-2:]),([':','&&'],['&&',':']),'Ninja linker framing')
    shared = shared[2:-2]; prefix = native_link(dso)
    A.exact(shared[:len(prefix)],prefix,'exact candidate linker flags')
    A.exact(shared[-1],'-lm','candidate link library')
    A.exact(len(shared[len(prefix):-1]),17,'all candidate linked objects')
    A.exact(set(shared[len(prefix):-1]),set(deps),'shared/archive same objects')
    current_sources = {}
    for relative, paths in deps.items():
        A.require(str(ROOT/relative[len(PREFIX):-2]) in paths,'actual candidate source dependency')
        for path in paths:
            p = Path(path); inputs.add(p)
            if ROOT in p.parents and path not in current_sources:
                blob = R.command(['git','cat-file','blob',SOURCE_HEAD+':'+str(p.relative_to(ROOT))])
                declared = dict(path=path,bytes=len(blob),sha256=A.sha(blob))
                A.exact(O.pin(p),declared,'candidate source bound to committed production')
                current_sources[path] = declared
    changed, pins = [], []
    for obj in objects:
        pin = O.pin(obj); pins.append(pin); inputs.add(obj)
        A.exact(A.sha(R.command(['/usr/bin/ar','p',archive,obj.name])),pin['sha256'],'actual archive member')
        if pin['sha256'] != previous_objects[obj.name]['sha256']:
            changed.append(obj.name)
        else:
            relative = str(obj.relative_to(BASE))
            A.exact(set(deps[relative]),set(previous['dependencies'][relative]),'identical object compiler inputs')
            for path in deps[relative]:
                if path in current_sources:
                    A.exact(current_sources[path],previous_sources[path],'unchanged object committed source')
    A.exact(changed,['WirehairV2Profile.cpp.o'],'one changed production object')
    changed_object = BASE/PREFIX/'codec/WirehairV2Profile.cpp.o'
    relative = str(changed_object.relative_to(BASE))
    proof_object = dict(name='proof-profile.o',sha256=O.pin(changed_object)['sha256'])
    if proof_dir is not None:
        target,dep = proof_dir/proof_object['name'],proof_dir/'proof-profile.d'
        A.require(not target.exists() and not target.is_symlink() and not dep.exists() and not dep.is_symlink(),
                  'fresh changed-object proof')
        args = native_compile(relative)
        args[args.index('-o')+1] = str(target)
        args[args.index('-MF')+1] = str(dep)
        R.command(args)
        A.exact(O.pin(target)['sha256'],proof_object['sha256'],'changed object reproduced from committed sources')
        text = A.read_regular(dep,1024*1024).decode().replace('\\\n','')
        name,paths = text.split(':',1)
        A.exact(name,relative,'proof compile target')
        paths = [str(Path(p).resolve(strict=True)) for p in shlex.split(paths)]
        A.exact(set(paths),set(deps[relative]),'same proof compiler dependency closure')
        target.chmod(0o400); dep.chmod(0o400)
    inputs.update((archive,dso,BASE/'build.ninja',BASE/'.ninja_deps',BASE/'CMakeCache.txt',
                   BASE/'CMakeFiles/rules.ninja'))
    diagnostic = dict(claimed=LOST_DIAGNOSTIC,observed=O.pin(Path(LOST_DIAGNOSTIC['path'])),
                      reason='Prior ctest -N overwrote non-producing LastTest.log; historical receipt not wholly intact')
    inputs.add(Path(LOST_DIAGNOSTIC['path']))
    report = dict(source_head=SOURCE_HEAD,historical_source_blobs=sorted(current_sources.values(),key=lambda p:p['path']),
                  archive=O.pin(archive),objects=pins,dependencies=deps,compile_commands=commands[:-1],
                  shared_link=shared,original=O.pin(dso),proof_name='proof-new.so',proof_sha256=digest,
                  proof_object=proof_object,changed_from_admission=changed,
                  admission_claim=previous['historical_claim'],admission_archive=previous['archive'],
                  excluded_historical_diagnostics=[diagnostic])
    result = [historical[0],report]
    if proof_dir is not None:
        for library in result: relink(library,proof_dir)
    return result,inputs


SETTINGS = R.Configuration(PROTOCOL,OUTPUT,LIBRARIES,NEW,provenance)


def main(settings=SETTINGS):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    b = sub.add_parser('build'); b.add_argument('mode',choices=('native','asan-driver')); b.add_argument('output',type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir',type=Path); r.add_argument('output',type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt',type=Path)
    sub.add_parser('replay')
    args = parser.parse_args()
    if args.command == 'build': R.build(args.mode,args.output,settings)
    elif args.command == 'receipt': A.publish(args.output,A.canonical(R.receipt(args.build_dir,settings)))
    elif args.command == 'run': R.run(args.receipt,settings)
    else:
        result = R.replay(settings)
        print(json.dumps(dict(outcome=result['outcome'],exact_replay=True)))


if __name__ == '__main__': main()
