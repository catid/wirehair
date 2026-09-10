#!/usr/bin/env python3
"""Neutral relocated installed-shared consumers; never launches timing or loss science."""
import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import Qualify as Q

P, A, B = Q.P, Q.A, Q.B
ORIGINAL = P.CURRENT/'package-shared/relocated-prefix'
LIBRARY = Path('lib64/wirehair-e2e/libwirehair.so.2.0.0')
DIRECT = {'package_c_consumer', 'package_cpp_consumer', 'package_k6_consumer',
          'package_small_consumer', 'package_v2_small_consumer', 'package_v2_k5_consumer',
          'package_v2_k8_consumer', 'package_ordinary_cpp'}


def names(arm):
    A.require(arm in ('baseline', 'candidate'), 'explicit package arm')
    return DIRECT | ({'package_ordinary_c'} if arm == 'candidate' else set())


def check_tests(raw, arm):
    suite = ET.fromstring(raw)
    expected = names(arm) | {'package_plugin_round_trip'}
    A.exact(suite.tag, 'testsuite', 'package JUnit root')
    A.exact(suite.attrib['tests'], str(len(expected)), 'package test count')
    for key in ('failures', 'disabled', 'skipped', 'errors'):
        A.exact(suite.attrib.get(key, '0'), '0', 'no unsuccessful package tests')
    cases = list(suite.findall('testcase'))
    A.exact(len(cases), len(expected), 'unique package test count')
    A.exact({case.attrib['name'] for case in cases}, expected, 'exact package test roster')
    for case in cases:
        A.require(case.attrib.get('status') == 'run' and not any(
            case.find(tag) is not None for tag in ('failure', 'skipped', 'error')),
            'every package test actually passed')


def inventory(prefix):
    """Record bytes and symlink spelling; reject special files or escaped links."""
    result = {}
    for directory, directories, files in os.walk(prefix, followlinks=False):
        for name in directories:
            A.require(not (Path(directory)/name).is_symlink(), 'no symlinked package directories')
        for name in files:
            path = Path(directory)/name
            relative = str(path.relative_to(prefix))
            if path.is_symlink():
                target = os.readlink(path)
                A.require(not Path(target).is_absolute() and prefix in path.resolve(strict=True).parents,
                          'package symlink stays in its relocated prefix')
                result[relative] = dict(kind='symlink', target=target)
            else:
                pin = B.pin(path)
                result[relative] = dict(kind='file', bytes=pin['bytes'], sha256=pin['sha256'])
    A.require(result, 'nonempty package')
    return result


def check_overlay(original, observed, candidate=None):
    expected = {name: dict(row) for name, row in original.items()}
    if candidate is not None:
        expected[str(LIBRARY)].update(bytes=candidate['bytes'], sha256=candidate['sha256'])
    A.exact(observed, expected, 'package preserves every byte/link except the selected DSO')


def check_defined_symbols(raw, target):
    allowed = set()
    if target in ('package_c_consumer', 'package_cpp_consumer', 'package_plugin'):
        allowed.update(('wirehair_package_round_trip', 'wirehair_package_v2_round_trip',
                        'wirehair_package_v2_selector_failures'))
    if target == 'package_plugin': allowed.add('wirehair_plugin_round_trip')
    for line in raw.decode().splitlines():
        if not line.split(): continue
        name = line.split()[-1]
        if name.startswith(('wirehair_', 'gf256_')) or name == 'GF256Ctx':
            # Local compiler clones of the three known test helpers are not
            # codec implementations. No public codec symbol shares these names.
            A.require(name.split('.', 1)[0] in allowed, 'no statically compiled codec API/GF context')


def qualify(prepared, output):
    prepared = Path(prepared).resolve(strict=True)
    P.verify_prepared(prepared)
    candidate_pin = B.pin(prepared/'candidate'/LIBRARY.name)
    A.exact(candidate_pin['sha256'], P.CANDIDATE_DSO_SHA, 'exact prepared candidate bytes')
    output = Path(output)
    output = output.parent.resolve(strict=True)/output.name
    A.require(output != P.ROOT and P.ROOT not in output.parents and not output.exists()
              and not output.is_symlink(), 'fresh external package directory')
    original = inventory(ORIGINAL)
    A.exact(sum(row['kind'] == 'file' for row in original.values()), 19, 'retained package file count')
    A.exact(sum(row['kind'] == 'symlink' for row in original.values()), 2, 'retained SONAME links')
    A.exact(original[str(LIBRARY)]['sha256'], P.DSO_SHA, 'original package is baseline DSO')
    output.mkdir(mode=0o700)
    dependencies = Q.helper_files() | P.compiler_inputs() | {
        Path(__file__).resolve(), Path(sys.executable).resolve(strict=True), Path('/usr/bin/ctest'),
        prepared/'PREPARED.json', P.HERE/'package/CMakeLists.txt',
        P.HERE/'package/OrdinaryConsumer.cpp', P.HERE/'ProviderCheck.cpp',
        P.ROOT/'test/package/CMakeLists.txt', P.ROOT/'test/V2SmallK8CConsumer.c'}
    dependencies.update(ORIGINAL/name for name, row in original.items() if row['kind'] == 'file')
    for path in tuple(dependencies):
        if Q.runtime_elf(path):
            raw = B.command(['ldd', path])
            A.require(b'not found' not in raw, 'package qualifier runtime resolves')
            dependencies.update(Path(word).resolve(strict=True) for word in raw.decode().split() if word.startswith('/'))
    frozen, commands, results = {}, [], []
    B.freeze_inputs(dependencies, frozen)

    def run(argv, cwd=P.ROOT, expected=0, environment=None):
        argv = list(map(str, argv))
        dependencies.add(Path(argv[0]).resolve(strict=True))
        B.freeze_inputs(dependencies, frozen)
        env = B.process_environment() if environment is None else environment
        index = len(commands)
        record = dict(argv=argv, cwd=str(cwd), environment=env, expected=expected)
        stem = output/('command-%03d'%index)
        A.publish(str(stem)+'.request.json', A.canonical(record))
        process = subprocess.run(argv, cwd=cwd, env=env, stdin=subprocess.DEVNULL,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
        record['returncode'] = process.returncode
        for suffix, raw in (('stdout', process.stdout), ('stderr', process.stderr),
                            ('result.json', A.canonical(record))):
            A.publish(str(stem)+'.'+suffix, raw)
        commands.append(record)
        A.require(len(process.stdout) <= 16*1024**2 and len(process.stderr) <= 65536, 'package command output caps')
        A.require(process.returncode == expected and (expected != 0 or not process.stderr),
                  'package command failed: '+str(argv))
        B.freeze_inputs(dependencies, frozen)
        return process

    prefixes = {}
    for arm in ('baseline', 'candidate'):
        parent = output/arm
        parent.mkdir()
        staged = parent/'staging-prefix'
        shutil.copytree(ORIGINAL, staged, symlinks=True)
        check_overlay(original, inventory(staged))
        source = prepared/arm/LIBRARY.name
        dependencies.add(source)
        B.freeze_inputs(dependencies, frozen)
        if arm == 'candidate': shutil.copyfile(source, staged/LIBRARY)
        prefix = parent/'relocated-prefix'
        staged.rename(prefix)
        check_overlay(original, inventory(prefix), candidate_pin if arm == 'candidate' else None)
        A.exact(B.pin(prefix/LIBRARY)['sha256'], B.pin(source)['sha256'], 'exact prepared DSO in relocated package')
        dependencies.update(prefix/name for name, row in original.items() if row['kind'] == 'file')
        prefixes[arm] = prefix

    for arm in ('baseline', 'candidate'):
        prefix = prefixes[arm]
        build = output/arm/'build'
        dso = prefix/LIBRARY
        run(['/usr/bin/cmake', '-S', P.HERE/'package', '-B', build, '-G', 'Ninja',
             '-DCMAKE_BUILD_TYPE=Release', '-DWH2_PACKAGE_ARM='+arm,
             '-DWH2_RELOCATED_PREFIX='+str(prefix)])
        dependencies.update(build/name for name in ('compile_commands.json', 'CMakeCache.txt', 'build.ninja',
                            'CMakeFiles/rules.ninja', 'CTestTestfile.cmake', 'upstream/CTestTestfile.cmake',
                            'imported-target.txt'))
        # Ninja's configure dependency roster closes the actual CMake inputs,
        # including installed package exports and system CMake modules.
        configure = [line for line in (build/'build.ninja').read_text().splitlines()
                     if line.startswith('build build.ninja: RERUN_CMAKE | ')]
        A.exact(len(configure), 1, 'one configure input roster')
        for name in shlex.split(configure[0].split(' | ', 1)[1]):
            path = Path(name)
            dependencies.add((path if path.is_absolute() else build/path).resolve(strict=True))
        entries = A.decode(A.read_regular(build/'compile_commands.json', 4*1024**2))
        targets = names(arm) | {'package_plugin', 'package_plugin_host'}
        counts = {name: 2 for name in targets}
        counts.update(package_c_consumer=3, package_cpp_consumer=3, package_plugin=3, package_plugin_host=1)
        observed, headers, plans = {}, {}, []
        for index, entry in enumerate(entries):
            target = next(part[:-4] for part in Path(entry['output']).parts if part.endswith('.dir'))
            observed[target] = observed.get(target, 0) + 1
            argv = shlex.split(entry['command'])
            A.exact(argv[-4:], ['-o', entry['output'], '-c', entry['file']], 'package compiler recipe')
            A.exact(entry['directory'], str(build), 'package compiler working directory')
            obj, dep = build/entry['output'], output/arm/('consumer-%03d.d'%index)
            run(argv[:-4]+['-M', '-MT', obj, '-MF', dep, entry['file']], cwd=build)
            inputs = B.preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj)
            headers.setdefault(target, set()).update(inputs)
            A.require(not any(P.ROOT/'include' in path.parents for path in inputs), 'no repository public-header fallback')
            if target != 'package_plugin_host':
                A.require('-DWIREHAIR_DLL=1' in argv, 'actual installed shared interface')
            dependencies.update(inputs)
            dependencies.add(dep)
            plans.append((argv, obj, inputs))
        A.exact(observed, counts, 'exact package consumer/TU roster')
        for target in targets - {'package_plugin_host'}:
            A.require(prefix/'include/wirehair/wirehair.h' in headers[target], 'target consumes installed public header')
        A.require(prefix/'include/wirehair/wirehair.hpp' in headers['package_ordinary_cpp'],
                  'ordinary C++ consumer uses installed wrapper')
        run(['/usr/bin/cmake', '--build', build, '-j', '8'])
        maps = sorted(build.glob('*.map'))
        A.exact({path.stem for path in maps}, targets, 'complete package link maps')
        negatives = []
        for path in maps:
            name = path.stem
            binary = build/'upstream'/('libpackage_plugin.so' if name == 'package_plugin' else name)
            if name in ('package_ordinary_cpp', 'package_ordinary_c'): binary = build/name
            dependencies.update((path, binary))
            loads = [line[5:] for line in path.read_text().splitlines() if line.startswith('LOAD ')]
            resolved_loads = [(Path(p) if Path(p).is_absolute() else build/p).resolve(strict=True) for p in loads]
            dependencies.update(resolved_loads)
            project = [str(p) for p in resolved_loads if 'libwirehair' in p.name]
            A.exact(project, [] if name == 'package_plugin_host' else [str(dso)], 'only relocated shared codec linked')
            check_defined_symbols(run(['/usr/bin/nm', '--defined-only', binary]).stdout, name)
            raw = run(['/usr/bin/ldd', binary]).stdout
            A.require(b'not found' not in raw, 'package consumer runtime resolves')
            dependencies.update(Path(word).resolve(strict=True) for word in raw.decode().split() if word.startswith('/'))
            providers = [line.split()[2] for line in raw.decode().splitlines()
                         if line.strip().startswith('libwirehair.so.2 =>')]
            A.exact([str(Path(p).resolve(strict=True)) for p in providers],
                    [] if name == 'package_plugin_host' else [str(dso)], 'exact relocated runtime provider')
            if name == 'package_plugin_host': continue
            env = B.process_environment()
            env['LD_LIBRARY_PATH'] = str(prefixes['candidate' if arm == 'baseline' else 'baseline']/LIBRARY.parent)
            argv = [binary] if name != 'package_plugin' else [build/'upstream/package_plugin_host', binary]
            negative = run(argv, expected=1, environment=env)
            A.exact(negative.stdout, b'', 'wrong-provider test never reaches application output')
            A.exact(negative.stderr, b'Shared test resolved the wrong Wirehair provider\n', 'wrong relocated provider rejected')
            negatives.append(name)
        junit = output/arm/'ctest-results.xml'
        run(['/usr/bin/ctest', '--test-dir', build, '--output-on-failure', '--no-tests=error',
             '--timeout', '30', '--output-junit', junit, '-j', '4'])
        check_tests(A.read_regular(junit, 4*1024**2), arm)
        for index, (argv, obj, before) in enumerate(plans):
            dep = output/arm/('consumer-%03d-after.d'%index)
            run(argv[:-4]+['-M', '-MT', obj, '-MF', dep, argv[-1]], cwd=build)
            A.exact(B.preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj), before, 'unchanged package compiler closure')
            dependencies.update((dep, obj))
        dependencies.update((junit, build/'Testing/Temporary/LastTest.log'))
        results.append(dict(arm=arm, prefix=str(prefix), files=inventory(prefix), junit=B.pin(junit),
                            consumers=sorted(targets), negatives=sorted(negatives)))
    A.exact(inventory(ORIGINAL), original, 'historical package remains unchanged')
    for row in results:
        check_overlay(original, inventory(Path(row['prefix'])), candidate_pin if row['arm'] == 'candidate' else None)
    P.verify_prepared(prepared)
    dependencies.update(output/('command-%03d.%s'%(index, suffix)) for index in range(len(commands))
                        for suffix in ('request.json', 'stdout', 'stderr', 'result.json'))
    B.freeze_inputs(dependencies, frozen)
    A.publish(output/'PACKAGE_QUALIFIED.json', A.canonical(dict(schema=1,
        scope='native relocated C/C++/plugin shared consumers only; no timing or default promotion',
        prepared=B.pin(prepared/'PREPARED.json'), original=dict(prefix=str(ORIGINAL), files=original),
        packages=results, commands=commands, files=[frozen[p] for p in sorted(frozen)])))
    print('PASS relocated baseline/candidate packages, 19 CTests and 19 wrong-provider rejections')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    qualify(args.prepared, args.output)
