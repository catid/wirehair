#!/usr/bin/env python3
"""Build and check direct native DSO consumers in a fresh directory; no science."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys
import types
import xml.etree.ElementTree as ET

import Prepare as P
import Bindings as L

A, B = P.A, P.B


def helper_files():
    """Include nested module_from_spec helpers omitted from sys.modules."""
    pending = [P, L]
    seen, paths = set(), set()
    while pending:
        current = pending.pop()
        if id(current) in seen:
            continue
        seen.add(id(current))
        path = getattr(current, '__file__', None)
        if path and Path(path).is_file():
            paths.add(Path(path).resolve(strict=True))
        pending.extend(value for value in vars(current).values() if isinstance(value, types.ModuleType))
    return paths | B.imported_files() | {Path(__file__).resolve(strict=True)}


def consumer_names():
    common = ('k3', 'k5', 'k8', 'borrowed', 'metadata', 'legacy_c', 'legacy_cpp',
              'small_c', 'k6_c', 'k3_c', 'k5_c', 'k8_explicit_c', 'certified')
    return {arm+'_'+name for arm in ('baseline', 'candidate') for name in common} | {
        'baseline_profile', 'candidate_contracts', 'candidate_borrowed_k8', 'candidate_ordinary_c'}


def test_names():
    return (consumer_names()-{'baseline_certified', 'candidate_certified'}) | {
        'baseline_exports', 'candidate_exports', 'certified_parity'}


def runtime_elf(path):
    raw = A.read_regular(path, 256*1024**2, installed=True)
    # Executable permission is not required for dlopen: Python extension DSOs
    # are commonly mode 0644. Relocatable objects are not ldd inputs.
    return len(raw) >= 20 and raw[:4] == b'\x7fELF' and raw[5] == 1 and raw[16:18] in (b'\x02\x00', b'\x03\x00')


def check_test_results(raw):
    suite = ET.fromstring(raw)
    A.exact(suite.tag, 'testsuite', 'CTest JUnit root')
    A.exact(suite.attrib['tests'], '31', 'complete native shared test count')
    A.exact(suite.attrib['failures'], '0', 'no failed native shared tests')
    A.exact(suite.attrib['disabled'], '0', 'no disabled native shared tests')
    A.exact(suite.attrib.get('skipped', '0'), '0', 'no skipped native shared tests')
    cases = list(suite.findall('testcase'))
    A.exact(len(cases), 31, 'exact testcase record count')
    A.exact({case.attrib['name'] for case in cases}, test_names(), 'exact passed native shared roster')
    for case in cases:
        A.require(case.attrib.get('status') == 'run' and case.find('failure') is None and
                  case.find('skipped') is None and case.find('error') is None, 'every required test actually passed')


def qualify(prepared, output):
    prepared = Path(prepared).resolve(strict=True)
    proof = P.verify_prepared(prepared)
    output = Path(output)
    output = output.parent.resolve(strict=True)/output.name
    A.require(P.ROOT != output and P.ROOT not in output.parents and not output.exists() and
              not output.is_symlink(), 'fresh external neutral test directory')
    output.mkdir(mode=0o700)
    dependencies = helper_files() | {prepared/'PREPARED.json', P.HERE/'CMakeLists.txt',
                                    P.HERE/'ProviderCheck.cpp', Path(sys.executable).resolve(strict=True),
                                    Path('/usr/bin/ctest').resolve(strict=True)}
    dependencies.update(P.ROOT/name for name in (
        'bench/Wh2K8OrdinarySelector/TestOverlay.cmake',
        'test/cmake/CheckElfExports.cmake', 'test/cmake/ParseElfExports.cmake',
        'bench/Wh2ProfileClassification/CheckParity.cmake'))
    # Producer/tool/runtime pins are inherited only after full current verification.
    dependencies.update(Path(row['path']) for row in proof['files'])
    for path in tuple(dependencies):
        if runtime_elf(path):
            raw = B.command(['ldd', path])
            A.require(b'not found' not in raw, 'resolved qualifier/ctypes runtime')
            dependencies.update(Path(word).resolve(strict=True) for word in raw.decode().split() if word.startswith('/'))
    frozen, commands = {}, []
    B.freeze_inputs(dependencies, frozen)

    def run(argv, expected=0, clean_stderr=True, cwd=P.ROOT):
        argv = list(map(str, argv))
        index = len(commands)
        record = dict(argv=argv, cwd=str(cwd), expected=expected)
        A.publish(output/('command-%03d.request.json'%index), A.canonical(record))
        B.freeze_inputs(dependencies, frozen)
        A.require(Path(argv[0]).resolve(strict=True) in frozen, 'invoked tool is frozen before execution')
        result = subprocess.run(argv, cwd=cwd, env=B.process_environment(), stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
        A.publish(output/('command-%03d.stdout'%index), result.stdout)
        A.publish(output/('command-%03d.stderr'%index), result.stderr)
        record['returncode'] = result.returncode
        commands.append(record)
        A.publish(output/('command-%03d.result.json'%index), A.canonical(record))
        A.require(result.returncode == expected and (not clean_stderr or not result.stderr),
                  'neutral qualification command failed: '+str(argv))
        A.require(len(result.stdout) <= 16*1024**2 and len(result.stderr) <= 65536, 'neutral output cap')
        B.freeze_inputs(dependencies, frozen)
        return result.stdout

    build = output/'build'
    run(['/usr/bin/cmake', '-S', P.HERE, '-B', build, '-G', 'Ninja', '-DCMAKE_BUILD_TYPE=Release',
         '-DPython3_EXECUTABLE='+str(Path(sys.executable).resolve(strict=True)),
         '-DWH2_K8_SHARED_PREPARED='+str(prepared)])
    database = build/'compile_commands.json'
    dependencies.update((database, build/'build.ninja', build/'CMakeFiles/rules.ninja',
                         build/'CTestTestfile.cmake', build/'CMakeCache.txt'))
    entries = A.decode(A.read_regular(database, 4*1024**2))
    A.exact(len(entries), 60, 'two translation units per direct consumer')
    targets = [row['output'].split('/')[1] for row in entries]
    A.require(all(target.endswith('.dir') for target in targets), 'consumer object directory shape')
    A.exact({target[:-4] for target in targets}, consumer_names(), 'exact consumer target roster')
    A.require(all(targets.count(name+'.dir') == 2 for name in consumer_names()), 'two TUs for every consumer')
    roster = A.decode(run(['/usr/bin/ctest', '--test-dir', build, '--show-only=json-v1']))
    A.exact(len(roster['tests']), 31, 'complete declared CTest count')
    A.exact({row['name'] for row in roster['tests']}, test_names(), 'exact declared CTest roster')
    plans = []
    for index, entry in enumerate(entries):
        argv = shlex.split(entry['command'])
        A.exact(set(entry), {'command', 'directory', 'file', 'output'}, 'CMake consumer recipe schema')
        A.exact(entry['directory'], str(build), 'new consumer build directory')
        A.exact(argv[-4:], ['-o', entry['output'], '-c', entry['file']], 'consumer compiler recipe')
        A.require(argv[0] in ('/usr/bin/cc', '/usr/bin/c++') and
                  not any('WIREHAIR_STATIC' in arg for arg in argv), 'shared-only consumer definitions')
        obj, dep = build/entry['output'], output/('consumer-%03d.d'%index)
        run(argv[:-4]+['-M', '-MT', obj, '-MF', dep, entry['file']], cwd=build)
        before = B.preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj)
        dependencies.update(before)
        B.freeze_inputs(dependencies, frozen)
        plans.append(dict(recipe=entry, object=obj, before=before, dependency=dep))
    # All consumer sources and headers are pinned before the first compilation.
    run(['/usr/bin/cmake', '--build', build, '-j', '8'])
    maps = sorted(build.glob('*.map'))
    A.exact({p.stem for p in maps}, consumer_names(), 'all direct shared-consumer link maps')
    consumer_proofs, negatives, artifacts = [], [], set()
    for path in maps:
        name = path.stem
        arm = name.split('_', 1)[0]
        A.require(arm in ('baseline', 'candidate'), 'explicit consumer arm')
        dso = prepared/arm/'libwirehair.so.2.0.0'
        loads = [line[5:] for line in path.read_text().splitlines() if line.startswith('LOAD ')]
        project_libraries = [p for p in loads if ('wirehair' in Path(p).name and
                             Path(p).suffix in ('.a', '.o')) or 'libwirehair.so' in Path(p).name]
        A.exact(project_libraries, [str(dso)], 'one actual DSO and no static codec link')
        executable = build/name
        dependencies.add(executable)
        symbols = run(['/usr/bin/nm', '--defined-only', executable]).decode()
        names = [line.split()[-1] for line in symbols.splitlines() if line.split()]
        A.require(not any(name.startswith('wirehair_') or name == 'GF256Ctx' for name in names),
                  'consumer defines no codec API/private GF context')
        runtime = run(['/usr/bin/ldd', executable]).decode()
        resolved = [line.split()[2] for line in runtime.splitlines() if line.strip().startswith('libwirehair.so.2 =>')]
        A.require(len(resolved) == 1 and Path(resolved[0]).resolve(strict=True) == dso,
                  'SONAME resolves exact intended arm')
        # Explicitly force the wrong DSO in a fresh process. ProviderCheck must
        # reject it before main, even for tests whose bytes happen to agree.
        wrong = 'candidate' if arm == 'baseline' else 'baseline'
        env = B.process_environment()
        env['LD_LIBRARY_PATH'] = str(prepared/wrong)
        negative = subprocess.run([str(executable)], cwd=P.ROOT, env=env, stdin=subprocess.DEVNULL,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
        A.publish(output/(name+'-wrong-provider.stdout'), negative.stdout)
        A.publish(output/(name+'-wrong-provider.stderr'), negative.stderr)
        A.publish(output/(name+'-wrong-provider.result.json'), A.canonical(dict(returncode=negative.returncode)))
        A.require(negative.returncode == 1 and not negative.stdout and
                  negative.stderr == b'Shared test resolved the wrong Wirehair provider\n',
                  'positive wrong-provider detector')
        captures = [output/(name+'-wrong-provider.'+suffix) for suffix in ('stdout', 'stderr', 'result.json')]
        artifacts.update(captures)
        negatives.append(dict(name=name, returncode=negative.returncode, wrong_provider=str(prepared/wrong),
                              captures=[B.pin(p) for p in captures]))
        dependencies.add(executable)
        consumer_proofs.append(dict(name=name, executable=B.pin(executable), map=B.pin(path), provider=B.pin(dso)))
    for reverse in (False, True):
        raw = run([sys.executable, P.HERE/'Bindings.py', prepared]+(['--reverse'] if reverse else []))
        report = A.decode(raw)
        A.exact(report['cases'], [48, 48], 'both actual DSO API rosters')
        A.publish(output/('bindings-reverse.json' if reverse else 'bindings-forward.json'), raw)
        artifacts.add(output/('bindings-reverse.json' if reverse else 'bindings-forward.json'))
    junit = output/'ctest-results.xml'
    run(['/usr/bin/ctest', '--test-dir', build, '--output-on-failure', '--no-tests=error',
         '--output-junit', junit, '-j', '4'])
    check_test_results(A.read_regular(junit, 4*1024**2))
    artifacts.add(junit)
    # Compare compiler dependency closure after building, without overwriting
    # CMake's own dependency records or any retained historical artifact.
    for index, plan in enumerate(plans):
        entry = plan['recipe']
        argv = shlex.split(entry['command'])
        after_path = output/('consumer-%03d-after.d'%index)
        run(argv[:-4]+['-M', '-MT', plan['object'], '-MF', after_path, entry['file']], cwd=build)
        A.exact(B.preprocessor_dependencies(A.read_regular(after_path, 2*1024**2), plan['object']),
                plan['before'], 'unchanged consumer compiler closure')
        dependencies.add(plan['object'])
    P.verify_prepared(prepared)
    artifacts.update(output/('command-%03d.%s'%(index, suffix)) for index in range(len(commands))
                     for suffix in ('request.json', 'stdout', 'stderr', 'result.json'))
    artifacts.update(maps)
    dependencies.update(artifacts)
    B.freeze_inputs(dependencies, frozen)
    result = dict(schema=1, scope='native actual shared API/ownership/provider qualification; no timing',
                  prepared=B.pin(prepared/'PREPARED.json'), commands=commands, consumers=consumer_proofs,
                  negatives=negatives, junit=B.pin(junit),
                  tests=B.pin(build/'Testing/Temporary/LastTest.log'),
                  files=[frozen[p] for p in sorted(frozen)])
    A.publish(output/'QUALIFIED.json', A.canonical(result))
    print('PASS direct native shared consumers, both loader orders, provider negatives and certified bytes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    qualify(args.prepared, args.output)
