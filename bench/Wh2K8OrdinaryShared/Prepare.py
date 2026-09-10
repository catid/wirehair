#!/usr/bin/env python3
"""Reproduce the native DSO and replace only its Profile object; never time code.

All commands write into a fresh external directory. Historical recipes, objects,
and qualification logs are read-only inputs, not old verifiers to rerun.
"""
import argparse
import importlib.util
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


B = module('_k8_shared_build_helpers', ROOT/'bench/Wh2K8PublicCostBuildR0.py')
C = module('_k8_shared_candidate', ROOT/'bench/Wh2K8OrdinaryCost/Candidate.py')
A = B.A
CURRENT = B.PRODUCTION/'native'
DSO_SHA = 'ef684abac606667e30e3b0de1204b0897aa0fa93a4f5ec2268d1ba951dce03b2'
CANDIDATE_DSO_SHA = 'fe7527a7470e7ea8e1790eb84b3761689e2bb61c883d6df905a1b3889699df3b'
PROFILE_INDEX = B.PRODUCERS.index('codec/WirehairV2Profile.cpp')
LINK_INPUTS = B.LINK_INPUTS + ('crtbeginS.o', 'crtendS.o')


def compiler_inputs():
    result = {(Path('/usr/bin')/name).resolve(strict=True) for name in
              ('c++', 'cc', 'as', 'ld', 'ar', 'ranlib', 'ninja', 'nm', 'ldd', 'cmake', 'readelf')}
    for name in ('cc1', 'cc1plus', 'collect2'):
        result.add(Path(B.command(['c++', '-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    for name in LINK_INPUTS:
        path = Path(B.command(['c++', '-print-file-name='+name]).decode().strip())
        A.require(path.is_absolute(), 'resolved link input '+name)
        result.add(path.resolve(strict=True))
    for path in tuple(result):
        raw = A.read_regular(path, 256*1024**2, installed=True)
        if len(raw) >= 20 and raw[:4] == b'\x7fELF' and raw[5] == 1 and raw[16:18] in (b'\x02\x00', b'\x03\x00'):
            linked = B.command(['ldd', path])
            A.require(b'not found' not in linked, 'resolved compiler runtime')
            result.update(Path(word).resolve(strict=True) for word in linked.decode().split() if word.startswith('/'))
    return result


def dependency_paths(output):
    return [output/(Path(name).name+'.o.d') for name in B.PRODUCERS] + [output/(C.SOURCE_NAME+'.d')]


def check_dependencies(raw, target, declared, pins):
    actual = B.preprocessor_dependencies(raw, target)
    A.exact({row['path'] for row in declared}, {str(path) for path in actual},
            'complete actual compiler dependencies')
    A.exact(len(declared), len(actual), 'unique compiler dependencies')
    for row in declared:
        A.exact(row, pins[row['path']], 'compiler dependency pin cross-link')


def link_recipe(raw, recipes):
    lines = raw.decode('utf-8').splitlines()
    A.exact(len(lines), 20, 'nineteen producer commands and one shared link')
    for line, (source, obj, flags) in zip(lines[:-1], recipes):
        relative = str(obj.relative_to(CURRENT))
        A.exact(shlex.split(line), ['/usr/bin/c++'] + flags +
                ['-MD', '-MT', relative, '-MF', relative+'.d', '-o', relative, '-c', str(source)],
                'shared build uses the same original producer recipe')
    framed = shlex.split(lines[-1])
    A.exact((framed[:2], framed[-2:]), ([':', '&&'], ['&&', ':']), 'Ninja link framing')
    expected = ['/usr/bin/c++', '-fPIC', '-O3', '-DNDEBUG',
                '-Wl,--version-script='+str(ROOT/'abi/wirehair.map'), '-shared',
                '-Wl,-soname,libwirehair.so.2', '-o', 'libwirehair.so.2.0.0']
    expected += [str(obj.relative_to(CURRENT)) for _, obj, _ in recipes] + ['-lm']
    A.exact(framed[2:-2], expected, 'exact original native shared link')
    return expected


def replacement_objects(objects, replacement):
    A.exact(len(objects), 19, 'complete shared object roster')
    A.exact(objects[PROFILE_INDEX].name, 'WirehairV2Profile.cpp.o', 'original Profile position')
    result = list(objects)
    result[PROFILE_INDEX] = replacement
    A.exact(sum(a != b for a, b in zip(objects, result)), 1, 'only Profile is replaced')
    return result


def proof_shape(proof, output):
    A.exact(set(proof), {'schema', 'scope', 'original_dso', 'original_link', 'members',
                        'candidate', 'links', 'commands', 'files', 'captures', 'interpreter'}, 'preparation schema')
    A.exact(proof['schema'], 3, 'preparation version')
    A.exact([row['arm'] for row in proof['links']], ['baseline', 'candidate'], 'prepared arm order')
    for row in proof['links']:
        A.exact(row['dso']['path'], str(output/row['arm']/'libwirehair.so.2.0.0'), 'exact prepared DSO path')
    A.exact(proof['links'][0]['dso']['sha256'], DSO_SHA, 'exact baseline identity')
    A.exact(proof['links'][1]['dso']['sha256'], CANDIDATE_DSO_SHA, 'independently inspected candidate DSO')
    A.exact(proof['candidate']['object']['sha256'], C.OBJECTS['native'][1], 'exact qualified selector object')
    A.exact([row['source']['path'] for row in proof['members']],
            [str(ROOT/name) for name in B.PRODUCERS], 'nineteen producing source order')
    A.exact(len(proof['commands']), 62, 'complete preparation command count')
    A.require(proof['original_link'] and proof['files'], 'nonempty producing graph')
    pins = B.pin_map(proof['files'])
    mandatory = {HERE/'Prepare.py', HERE/'CMakeLists.txt', HERE/'Bindings.py',
                 Path(A.__file__).resolve(), Path(B.__file__).resolve(), Path(C.__file__).resolve(),
                 CURRENT/'build.ninja', CURRENT/'CMakeFiles/rules.ninja', CURRENT/'compile_commands.json',
                 CURRENT/'libwirehair.so.2.0.0', CURRENT/'libwirehair.a', ROOT/'abi/wirehair.map',
                 B.PRODUCTION/'NEUTRAL_QUALIFIED.json'}
    mandatory.add(Path(proof['interpreter']['path']))
    mandatory.update(compiler_inputs())
    for row in proof['members']:
        mandatory.update(Path(row[key]['path']) for key in ('source', 'original', 'reproduced'))
        A.require(row['dependencies'], 'nonempty producer dependency closure')
        mandatory.update(Path(p['path']) for p in row['dependencies'])
    for key in ('source', 'object', 'qualified_object', 'compile_database', 'qualification_log'):
        mandatory.add(Path(proof['candidate'][key]['path']))
    mandatory.update(Path(p['path']) for p in proof['candidate']['dependencies'])
    mandatory.update(Path(row['dso']['path']) for row in proof['links'])
    A.require({str(p) for p in mandatory} <= set(pins), 'all mandatory producing inputs and outputs pinned')
    captures = B.pin_map(proof['captures'])
    expected = {str(output/('command-%03d.%s'%(i, suffix))) for i in range(62)
                for suffix in ('request.json', 'stdout', 'stderr', 'result.json')}
    expected.update(str(path) for path in dependency_paths(output))
    A.exact(set(captures), expected, 'all command captures pinned exactly once')
    A.exact(len(proof['captures']), len(expected), 'unique capture roster')


def expected_commands(output, recipes, original_link):
    plan = [(['/usr/bin/ninja', '-C', str(CURRENT), '-t', 'commands', 'libwirehair.so.2.0.0'], ROOT)]
    for source, old, flags in recipes:
        obj, dep = output/old.name, output/(old.name+'.d')
        plan.append((['/usr/bin/c++']+flags+['-M', '-MT', str(obj), '-MF', str(dep), str(source)], CURRENT))
    row, flags = C.recipe(A.decode(A.read_regular(C.QUALIFIED/'compile_commands.json', 4*1024**2)), 'native')
    obj, dep, source = output/(C.SOURCE_NAME+'.o'), output/(C.SOURCE_NAME+'.d'), C.QUALIFIED/C.SOURCE_NAME
    plan.append((['/usr/bin/c++']+flags+['-M', '-MT', str(obj), '-MF', str(dep), str(source)], C.QUALIFIED))
    for current_source, old, current_flags in recipes:
        current_obj, current_dep = output/old.name, output/(old.name+'.d')
        plan.append((['/usr/bin/c++']+current_flags+['-MD', '-MT', str(current_obj), '-MF', str(current_dep),
                    '-o', str(current_obj), '-c', str(current_source)], CURRENT))
        plan.append((['/usr/bin/ar', 'p', str(CURRENT/'libwirehair.a'), old.name], ROOT))
    plan.append((['/usr/bin/c++']+flags+['-MD', '-MT', str(obj), '-MF', str(dep),
                '-o', str(obj), '-c', str(source)], C.QUALIFIED))
    objects = [output/old.name for _, old, _ in recipes]
    for arm, inputs in (('baseline', objects), ('candidate', replacement_objects(objects, obj))):
        directory = output/arm
        plan.append((original_link[:8]+[str(directory/'libwirehair.so.2.0.0')]+list(map(str, inputs))+
                     ['-lm', '-Wl,-Map,'+str(directory/'link.map')], ROOT))
    return [dict(argv=argv, cwd=str(cwd), returncode=0, expected=0) for argv, cwd in plan]


def verify_prepared(output):
    output = Path(output).resolve(strict=True)
    proof = A.decode(A.read_regular(output/'PREPARED.json', 4*1024**2))
    proof_shape(proof, output)
    A.exact(proof['original_dso'], B.pin(CURRENT/'libwirehair.so.2.0.0'), 'fixed original DSO cross-link')
    A.exact(proof['interpreter'], B.pin(Path(proof['interpreter']['path'])), 'producer interpreter pin')
    C.source_identity(A.read_regular(ROOT/'codec/WirehairV2Profile.cpp', 1024**2),
                      A.read_regular(C.QUALIFIED/C.SOURCE_NAME, 1024**2))
    for record in proof['files']+proof['captures']:
        A.exact(B.pin(Path(record['path'])), record, 'unchanged preparation input/output')
    recipes = B.producer_recipes(A.decode(A.read_regular(CURRENT/'compile_commands.json', 4*1024**2)), 'native')
    original_link = link_recipe(A.read_regular(output/'command-000.stdout', 1024**2), recipes)
    A.exact(proof['original_link'], original_link, 'authenticated original shared link')
    A.exact(proof['commands'], expected_commands(output, recipes, original_link), 'exact complete producing chronology')
    pins = B.pin_map(proof['files'])
    candidate = proof['candidate']
    candidate_database = C.QUALIFIED/'compile_commands.json'
    A.exact(B.pin(candidate_database)['sha256'], C.DATABASE_SHA, 'fixed candidate compile database')
    candidate_row, _ = C.recipe(A.decode(A.read_regular(candidate_database, 4*1024**2)), 'native')
    candidate_fields = dict(source=C.QUALIFIED/C.SOURCE_NAME, object=output/(C.SOURCE_NAME+'.o'),
                            qualified_object=C.QUALIFIED/candidate_row['output'],
                            compile_database=candidate_database,
                            qualification_log=C.QUALIFIED/'Testing/Temporary/LastTest.log')
    for key, path in candidate_fields.items():
        A.exact(candidate[key], B.pin(path), 'prescribed candidate '+key)
        A.exact(candidate[key], pins[str(path)], 'candidate producing pin cross-link')
    A.exact(candidate['recipe'], candidate_row, 'candidate recipe cross-link')
    check_dependencies(A.read_regular(output/(C.SOURCE_NAME+'.d'), 2*1024**2),
                       output/(C.SOURCE_NAME+'.o'), candidate['dependencies'], pins)
    for index, record in enumerate(proof['commands']):
        A.exact(A.decode(A.read_regular(output/('command-%03d.result.json'%index), 65536)), record,
                'captured actual terminal command status')
        A.exact(A.decode(A.read_regular(output/('command-%03d.request.json'%index), 65536)),
                {key: record[key] for key in ('argv', 'cwd')}, 'captured command request')
        A.exact(A.read_regular(output/('command-%03d.stderr'%index), 65536), b'', 'empty successful command stderr')
    objects = [output/old.name for _, old, _ in recipes]
    for index, ((source, old, _), member) in enumerate(zip(recipes, proof['members'])):
        obj = objects[index]
        A.exact(member['original'], B.pin(old), 'original object cross-link')
        A.exact(member['reproduced'], B.pin(obj), 'reproduced object cross-link')
        A.exact(member['source'], B.pin(source), 'source pin cross-link')
        check_dependencies(A.read_regular(output/(old.name+'.d'), 2*1024**2), obj, member['dependencies'], pins)
        A.exact(A.read_regular(obj, 4*1024**2), A.read_regular(old, 4*1024**2), 'reproduced producer bytes')
        A.exact(A.read_regular(output/('command-%03d.stdout'%(22+2*index)), 4*1024**2),
                A.read_regular(obj, 4*1024**2), 'captured original archive member')
    for row in proof['links']:
        A.exact(B.pin(Path(row['dso']['path'])), row['dso'], 'unchanged actual DSO')
        A.exact(row['map']['path'], str(output/row['arm']/'link.map'), 'exact link-map path')
        A.exact(B.pin(Path(row['map']['path'])), row['map'], 'unchanged actual link map')
        expected_index = 60 if row['arm'] == 'baseline' else 61
        A.exact(row['argv'], proof['commands'][expected_index]['argv'], 'shared link command cross-link')
        loads = [str(Path(line[5:]).resolve(strict=True))
                 for line in A.read_regular(row['map']['path'], 4*1024**2).decode().splitlines() if line.startswith('LOAD ')]
        A.require(set(loads) <= set(pins), 'closed shared linker LOAD inputs')
        candidate_obj = output/(C.SOURCE_NAME+'.o')
        wanted = objects if row['arm'] == 'baseline' else replacement_objects(objects, candidate_obj)
        observed = [Path(p) for p in loads if Path(p) in objects or Path(p) == candidate_obj]
        A.exact(observed, wanted, 'exact loaded nineteen-object order')
    return proof


def prepare(output):
    output = Path(output)
    output = output.parent.resolve(strict=True)/output.name
    A.require(output.is_absolute() and ROOT != output and ROOT not in output.parents and
              not output.exists() and not output.is_symlink(), 'fresh external output')
    A.require(not any(os.environ.get(k) for k in B.ENV_KEYS + ('LD_AUDIT', 'LD_DEBUG')),
              'no allocator or dynamic-loader overrides')
    inventory_path = B.PRODUCTION/'NEUTRAL_QUALIFIED.json'
    A.exact(B.pin(inventory_path)['sha256'], B.INVENTORY_SHA, 'retained K8 qualification')
    inventory = B.pin_map(A.decode(A.read_regular(inventory_path, 4*1024**2))['files'])
    database = CURRENT/'compile_commands.json'
    original_dso = CURRENT/'libwirehair.so.2.0.0'
    original_archive = CURRENT/'libwirehair.a'
    initial = {inventory_path, database, original_dso, original_archive,
               CURRENT/'build.ninja', CURRENT/'CMakeFiles/rules.ninja', ROOT/'CMakeLists.txt', ROOT/'abi/wirehair.map',
               HERE/'Prepare.py', HERE/'CMakeLists.txt', HERE/'Bindings.py'}
    for path in (database, original_dso, original_archive, ROOT/'CMakeLists.txt'):
        A.exact(B.pin(path), inventory[str(path)], 'retained original producing input')
    A.exact(B.pin(ROOT/'abi/wirehair.map')['sha256'],
            'eb15977acc8a0fba143f56e70333c6f55562081aeee2eee4b98f18062a3752dd', 'unchanged export map')
    A.exact(B.pin(database)['sha256'], B.DATABASES['native'], 'original compile database')
    A.exact(B.pin(original_dso)['sha256'], DSO_SHA, 'original native DSO')
    recipes = B.producer_recipes(A.decode(A.read_regular(database, 8*1024**2)), 'native')
    initial.add(Path(sys.executable).resolve(strict=True))
    initial.update(compiler_inputs())
    for path in tuple(initial):
        if os.access(str(path), os.X_OK) and A.read_regular(path, 256*1024**2, installed=True)[:4] == b'\x7fELF':
            raw = B.command(['ldd', path])
            A.require(b'not found' not in raw, 'resolved tool runtime')
            initial.update(Path(word).resolve(strict=True) for word in raw.decode().split() if word.startswith('/'))
    dependencies = initial | B.imported_files() | {Path(A.__file__).resolve(strict=True)}
    frozen, commands = {}, []
    B.freeze_inputs(dependencies, frozen)
    output.mkdir(mode=0o700)

    def run(argv, expected=0, clean_stderr=True, cwd=ROOT):
        argv = list(map(str, argv))
        record = dict(argv=argv, cwd=str(cwd))
        index = len(commands)
        commands.append(record)
        B.freeze_inputs(dependencies, frozen)
        A.publish(output/('command-%03d.request.json'%index), A.canonical(record))
        result = subprocess.run(argv, cwd=cwd, env=B.process_environment(), stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        A.publish(output/('command-%03d.stdout'%index), result.stdout)
        A.publish(output/('command-%03d.stderr'%index), result.stderr)
        record.update(returncode=result.returncode, expected=expected)
        A.publish(output/('command-%03d.result.json'%index), A.canonical(record))
        A.require(len(result.stdout) <= 16*1024**2 and len(result.stderr) <= 65536, 'command output cap')
        A.require(result.returncode == expected and (not clean_stderr or not result.stderr),
                  'preparation command failed: '+str(argv))
        B.freeze_inputs(dependencies, frozen)
        return result

    raw = run(['/usr/bin/ninja', '-C', CURRENT, '-t', 'commands', original_dso.name]).stdout
    original_link = link_recipe(raw, recipes)
    plans = []
    for source, old, flags in recipes:
        obj, dep = output/old.name, output/(old.name+'.d')
        run(['/usr/bin/c++']+flags+['-M', '-MT', obj, '-MF', dep, source], cwd=CURRENT)
        before = B.preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj)
        dependencies.update(before | {old})
        B.freeze_inputs(dependencies, frozen)
        plans.append((source, old, flags, obj, dep, before))
    # The candidate planner authenticates the same source, recipe, object and
    # original forty-test log as the completed static qualification, read-only.
    candidate = C.plan('native', output, dependencies, frozen, run, B.freeze_inputs,
                       B.preprocessor_dependencies, B.pin)
    members = []
    for source, old, flags, obj, dep, before in plans:
        run(['/usr/bin/c++']+flags+['-MD', '-MT', obj, '-MF', dep, '-o', obj, '-c', source], cwd=CURRENT)
        A.exact(B.preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj), before, 'compiler closure')
        A.exact(A.read_regular(obj, 64*1024**2), A.read_regular(old, 64*1024**2), 'byte-identical producer')
        A.exact(run(['/usr/bin/ar', 'p', original_archive, old.name]).stdout,
                A.read_regular(obj, 64*1024**2), 'matching original archive member')
        dependencies.add(obj)
        B.freeze_inputs(dependencies, frozen)
        members.append(dict(source=B.pin(source), original=B.pin(old), reproduced=B.pin(obj),
                            dependencies=[frozen[p] for p in sorted(before)]))
    candidate_proof = C.compile_candidate(candidate, dependencies, frozen, run, B.freeze_inputs,
                                          B.preprocessor_dependencies, B.pin)
    objects = [p[3] for p in plans]
    links = []
    for arm, inputs in (('baseline', objects), ('candidate', replacement_objects(objects, candidate['object']))):
        directory = output/arm
        directory.mkdir(mode=0o700)
        dso, link_map = directory/original_dso.name, directory/'link.map'
        argv = original_link[:8] + [str(dso)] + list(map(str, inputs)) + ['-lm', '-Wl,-Map,'+str(link_map)]
        run(argv)
        loads = [Path(line[5:]).resolve(strict=True) for line in link_map.read_text().splitlines()
                 if line.startswith('LOAD ')]
        A.require(set(loads) <= dependencies, 'all actual linker inputs frozen before link')
        A.exact([p for p in loads if p in objects or p == candidate['object']], inputs,
                'original nineteen-object order and exact one-object replacement')
        if arm == 'baseline':
            A.exact(A.read_regular(dso, 4*1024**2), A.read_regular(original_dso, 4*1024**2),
                    'byte-identical complete native DSO')
        else:
            A.exact(B.pin(dso)['sha256'], CANDIDATE_DSO_SHA, 'unchanged independently inspected shared candidate')
        dependencies.add(dso)
        B.freeze_inputs(dependencies, frozen)
        # Preserve the actual SONAME without changing it to aid side-by-side tests.
        (directory/'libwirehair.so.2').symlink_to(dso.name)
        links.append(dict(arm=arm, dso=B.pin(dso), map=B.pin(link_map), argv=argv))
    B.freeze_inputs(dependencies, frozen)
    captures = [B.pin(output/('command-%03d.%s'%(i, suffix))) for i in range(len(commands))
                for suffix in ('request.json', 'stdout', 'stderr', 'result.json')]
    captures.extend(B.pin(path) for path in dependency_paths(output))
    proof = dict(schema=3, scope='native shared producer qualification only; no timing or promotion',
                 original_dso=B.pin(original_dso), original_link=original_link, members=members,
                 candidate=candidate_proof, links=links, commands=commands, captures=captures,
                 interpreter=B.pin(Path(sys.executable).resolve(strict=True)),
                 files=[frozen[p] for p in sorted(frozen)])
    A.publish(output/'PREPARED.json', A.canonical(proof))
    verify_prepared(output)
    print('PASS exact nineteen producers, original DSO and selector; shared tests remain separate', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--verify', action='store_true', help='Check new neutral preparation only; never executes a codec')
    args = parser.parse_args()
    if args.verify:
        verify_prepared(args.output)
        print('PASS unchanged native shared preparation')
    else:
        prepare(args.output)
