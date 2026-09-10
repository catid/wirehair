"""Closed ordinary-K8 selector reproduction for the freshly derived builder."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shlex

ROOT = Path(__file__).resolve().parents[2]
QUALIFIED = Path('/tmp/wh2-k8-ordinary-qualified.97E2f2VQ')
DATABASE_SHA = '1bf3cc82076fb157b89a0cf94851ca6494546a1c3d35879d579cdbfc7819719b'
SOURCE_SHA = '172324fd048daf154267f60a6e76df1e1b9a2545b9eac49059d0be21cfe99a3f'
ORIGINAL_SHA = '975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15'
OBJECTS = {
    'native': ('k8_ordinary_candidate', '4c5e748f1fa618c88d8d9b8cd58045a1d4bf28211bceea63db4d726362ebf784'),
    'scalar': ('k8_ordinary_scalar', '78e0eb810723ed89cd7f18cbaa61984265578f5e8f242e06359a7d5f44bf493b'),
    'asan': ('k8_ordinary_asan', '59404f5d4f1572df89ed7d32425bf2c741cdea12759fd097f988deb6853b750a'),
}
SOURCE_NAME = 'WirehairV2K8OrdinarySelector.cpp'


def require(value, why):
    if not value:
        raise ValueError(why)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def verify_derivation(directory):
    spec = importlib.util.spec_from_file_location('_ordinary_k8_derivation', Path(__file__).with_name('Derive.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    directory = module.verify(directory)
    return {directory/name for name in module.encoded_outputs()}


def source_identity(original, candidate):
    require(digest(original) == ORIGINAL_SHA and digest(candidate) == SOURCE_SHA, 'exact selector source pins')
    anchor = b"""    if (SmallShape<3>(messageBytes, blockBytes)) {
        return wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_SMALL_K3_2026_09, message, messageBytes,
            blockBytes, serializedProfileOut, serializedProfileCapacity,
            serializedProfileBytesOut, codecOut);
    }
"""
    extra = anchor.replace(b'SmallShape<3>', b'SmallShape<8>').replace(b'SMALL_K3_', b'SMALL_K8_')
    require(original.count(anchor) == 1 and original.replace(anchor, anchor+extra) == candidate,
            'exact six-line selector derivation')


def recipe(database, mode):
    require(mode in OBJECTS, 'candidate backend')
    target = 'CMakeFiles/'+OBJECTS[mode][0]+'.dir/'+SOURCE_NAME+'.o'
    entries = [row for row in database if row.get('output') == target]
    require(len(entries) == 1, 'unique nonfault candidate recipe')
    row = entries[0]
    require(set(row) == {'directory', 'command', 'file', 'output'} and
            row['directory'] == str(QUALIFIED) and row['file'] == str(QUALIFIED/SOURCE_NAME),
            'candidate recipe binding')
    argv = shlex.split(row['command'])
    require(argv[0] == '/usr/bin/c++' and argv[-4:] == ['-o', target, '-c', row['file']] and
            not any('WIREHAIR_TESTING' in a or 'ENABLE_TEST_HOOKS' in a for a in argv),
            'actual production-style candidate recipe')
    flags = argv[1:-4]
    require(('-DANDROID=1' in flags) == (mode == 'scalar') and
            ('-fsanitize=address,undefined' in flags) == (mode == 'asan') and
            (mode != 'asan' or '-march=native' in flags), 'candidate private ABI/flags')
    return row, flags


def plan(mode, output, dependencies, frozen, run, freeze, parse_dependencies, pin):
    source = QUALIFIED/SOURCE_NAME
    original = ROOT/'codec/WirehairV2Profile.cpp'
    database_path = QUALIFIED/'compile_commands.json'
    raw_database = database_path.read_bytes()
    require(digest(raw_database) == DATABASE_SHA, 'qualified candidate compile database')
    source_identity(original.read_bytes(), source.read_bytes())
    row, flags = recipe(json.loads(raw_database), mode)
    old_object = QUALIFIED/row['output']
    require(pin(old_object)['sha256'] == OBJECTS[mode][1], 'qualified candidate object')
    log = QUALIFIED/'Testing/Temporary/LastTest.log'
    require(pin(log)['sha256'] == '0665ec2927728dc0f890acdbca5a9ffdd3d762a8de99e53a2b5295a547cb4eef',
            'complete forty-test qualification log')
    dependencies.update((source, original, old_object, database_path, log, Path(__file__).resolve()))
    freeze(dependencies, frozen)
    obj, dep = output/(SOURCE_NAME+'.o'), output/(SOURCE_NAME+'.d')
    run(['/usr/bin/c++']+flags+['-M', '-MT', str(obj), '-MF', str(dep), str(source)], cwd=QUALIFIED)
    before = parse_dependencies(dep.read_bytes(), obj)
    dependencies.update(before)
    freeze(dependencies, frozen)
    return dict(source=source, original=old_object, object=obj, dependency=dep,
                before=before, flags=flags, recipe=row, mode=mode, database=database_path, log=log)


def compile_candidate(plan, dependencies, frozen, run, freeze, parse_dependencies, pin):
    obj, dep, source = plan['object'], plan['dependency'], plan['source']
    run(['/usr/bin/c++']+plan['flags']+['-MD', '-MT', str(obj), '-MF', str(dep),
                                    '-o', str(obj), '-c', str(source)], cwd=QUALIFIED)
    require(parse_dependencies(dep.read_bytes(), obj) == plan['before'], 'candidate compiler closure')
    require(obj.read_bytes() == plan['original'].read_bytes(), 'byte-identical qualified candidate object')
    dependencies.add(obj)
    freeze(dependencies, frozen)
    return dict(source=pin(source), object=pin(obj), qualified_object=pin(plan['original']),
                recipe=plan['recipe'], compile_database=pin(plan['database']),
                qualification_log=pin(plan['log']), dependencies=[pin(p) for p in sorted(plan['before'])])


def check_link(text, archive, candidate, members):
    loads = [line[5:] for line in text.splitlines() if line.startswith('LOAD ')]
    require(loads.count(str(candidate)) == 1 and loads.count(str(archive)) == 1 and
            loads.index(str(candidate)) < loads.index(str(archive)), 'candidate precedes exact archive')
    selected = set(re.findall(re.escape(str(archive))+r'\(([^()\n]+)\)', text))
    require(selected and selected <= set(members) and 'WirehairV2Profile.cpp.o' not in selected,
            'old facade must not be extracted')
    require({'wirehair.cpp.o', 'gf256.cpp.o', 'WirehairSmall.cpp.o', 'WirehairSmallK8.cpp.o'} <= selected,
            'actual library owns WH1/core/K8/GF')
    return sorted(selected)
