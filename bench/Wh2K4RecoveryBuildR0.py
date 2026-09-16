#!/usr/bin/env python3
"""Build only the new K4 recovery observer; reuse authenticated production archives.

No historical controller/verifier is imported or invoked. Production recipes
are parsed as data; only the K4 boundary recipe is reproduced with fresh outputs.
The sole scientific
launch belongs to Wh2K4SerializedRecoveryR0.py, not this module.
"""
import importlib.util
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SPEC = importlib.util.spec_from_file_location('_k4_recovery_build_io', HERE/'Wh2AlignedIntermediateCostR0.py')
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
MODES = ('native', 'scalar', 'asan')
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
PROTOCOL = 'wirehair.wh2.k4-serialized-recovery-r0'
PRODUCTION = Path('/tmp/wh2-v2-k8-admission.YNN2XmAx')
PROOFS = Path('/tmp/wh2-k8-public-recovery-qualified.i7v7bPnq')
PRIOR = Path('/var/tmp/wh2-k8-public-recovery-r0')
PRIOR_PROTOCOL = 'wirehair.wh2.k8-public-recovery-r0'
PRIOR_COMPLETE_SHA = '699415eef17fcc486aeb3e5487feef0a9dcd4812110497df1370014c47b7af02'
PRIOR_CLAIM_SHA = 'b5335bb0aa2d0e8f354aebc2df5b3d8d51f65f7070e45b6e22b44499c5819ab1'
SOURCE_HEAD = '5a300505da3754dbdbb8cbfcba4b8a79bad1b0cb'
LIB_MANIFEST_SHA = dict(
    native='c5460ad57bc3328be535d0eb67f7f3139965e85d607e12d5b4991a8f0151ce01',
    scalar='b16e5ca24e87d611702cd318578278f75f8cb031d1636edaaaa0ccce71312f5a',
    asan='5838a98584c5805e905a1dac80f3a2bd7c9ce44c349aa82923f82cdc2da0ece9')
SMALL = Path('/tmp/wh2-k4-serialized.yegpP9w9')
BOUNDARY_AUDIT = Path('/tmp/wh2-k4-serialized-independent.VQMU5xiC/report-python312.json')
BOUNDARY_AUDIT_SHA = 'c0a5afa2e2250da819bde1b9eeefb9e13e895ce7a81666564fb143ff766054d9'
BOUNDARY_HEAD = '84eeea90c34332716c8ae4219989aca6b08ba2a8'
CORE_OLD_SHA = '26167117230258275cc0d522c3039564f5ec5abe8009cbc48faae0e7b9f04d51'
CORE_SHA = '5b0acdd096d24b76351bacd1718c44ca5b37d4df587fe7334822a9f61f1e0b8c'
FIXTURE_SHA = '608bafbe37ac0ba3aa94f5d030390af6cc623f9cf0791e86c5bcb8d72f277763'
SEALED = Path('/var/tmp/wh2-k4-thue-morse-r0')
COMPLETE_SHA = '8f1c9357b125edffc282df4cde29592130e27566f86946a7a85b0e92f175d281'
ARCHIVES = dict(native='0d4a0a7ec3cfa4b8e0f08df75e2c12935b51f28653c8c953a32bb49a5658eff5',
                scalar='bd7eb2eeb0eeaea15e11d131c6fba2921bf44cb1526dea8dff86171b496de867',
                asan='a7043eb65eb864aefd210f51b360d4591a365c60710c3aec87baf0abac3d4f1a')
GF_BYTES = dict(native=141328, scalar=137232, asan=157728)
ASAN_FLAGS = '-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -march=native'
LINK_INPUTS = ('crt1.o','crti.o','crtbegin.o','crtend.o','crtn.o','libstdc++.so',
               'libm.so','libm.so.6','libmvec.so','libgcc.a','libgcc_s.so','libgcc_s.so.1',
               'libc.so','libc.so.6','libc_nonshared.a','libpthread.a','liblto_plugin.so',
               'libasan_preinit.o','libasan.so','libubsan.so')
PRODUCERS = (
    'wirehair.cpp', 'codec/WirehairSmall.cpp', 'codec/WirehairSmallK5.cpp',
    'codec/WirehairSmallK8.cpp', 'codec/WirehairK6.cpp', 'codec/WirehairK6Core.cpp',
    'gf256.cpp', 'WirehairCodec.cpp', 'WirehairTools.cpp', 'codec/WirehairV2Codec.cpp',
    'codec/WirehairV2Peel.cpp', 'codec/WirehairV2Plan.cpp', 'codec/WirehairV2Policy.cpp',
    'codec/WirehairV2Precode.cpp', 'codec/WirehairV2PrecodeDecode.cpp',
    'codec/WirehairV2PrecodeEncode.cpp', 'codec/WirehairV2Profile.cpp',
    'codec/WirehairV2Seeds.cpp', 'codec/WirehairV2Solve.cpp')

def pin(path):
    path = Path(path)
    owned = any(base == path or base in path.parents for base in (ROOT, Path('/tmp'), Path('/var/tmp')))
    return A.pin(path, installed=not owned)


def command(argv):
    """Small read-only inspections; compilation uses the logged build runner."""
    argv = list(map(str, argv))
    if argv[0] in ('git', 'ar', 'nm', 'ldd', 'c++'):
        argv[0] = '/usr/bin/'+argv[0]
    result = subprocess.run(argv, cwd=ROOT, env=process_environment(), stdin=subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
    A.require(result.returncode == 0 and not result.stderr, 'inspection failed: '+str(argv))
    A.require(len(result.stdout) <= 32*1024**2, 'inspection output cap')
    return result.stdout


def process_environment():
    """Explicit compiler/observer environment excludes implicit search flags."""
    result = dict(PATH='/usr/bin:/bin', LANG='C', LC_ALL='C', TZ='UTC')
    result.update((key, os.environ[key]) for key in ('ASAN_OPTIONS', 'UBSAN_OPTIONS') if key in os.environ)
    return result


def freeze_inputs(paths, frozen):
    for path in sorted(paths):
        value = pin(path)
        if path in frozen:
            A.exact(value, frozen[path], 'input changed during qualification')
        else:
            frozen[path] = value


def pin_map(records):
    """Input/artifact overlap is permitted only for byte-identical records."""
    result = {}
    for record in records:
        A.exact(set(record), {'path', 'bytes', 'sha256'}, 'pin schema')
        path = Path(record['path'])
        A.require(path.is_absolute() and '..' not in path.parts, 'absolute pin path')
        A.integer(record['bytes'], 0, 256*1024**2)
        A.require(type(record['sha256']) is str and re.fullmatch('[0-9a-f]{64}', record['sha256']), 'pin digest')
        if str(path) in result:
            A.exact(result[str(path)], record, 'conflicting duplicate pin')
        result[str(path)] = record
    return result


def preprocessor_dependencies(raw, target):
    text = raw.decode().replace('\\\n', '')
    prefix = str(target)+': '
    A.require(text.startswith(prefix), 'exact dependency target')
    paths = shlex.split(text[len(prefix):])
    A.require(paths and all(Path(p).is_absolute() for p in paths), 'absolute compiler inputs')
    resolved = {Path(p).resolve(strict=True) for p in paths}
    for p in paths:
        original = Path(p)
        # Normalize harmless lexical '..' before checking for symlink redirects.
        lexical = Path(os.path.normpath(str(original)))
        if ROOT in lexical.parents:
            A.exact(lexical, original.resolve(strict=True), 'repository dependency redirect')
    return resolved


def producer_recipes(database, mode):
    """Validate the exact current CMake producer commands, not test-support TUs."""
    A.require(mode in MODES, 'producer backend')
    base = PRODUCTION/mode
    prefix = 'CMakeFiles/wirehair_objects.dir/' if mode=='native' else 'CMakeFiles/wirehair.dir/'
    entries = [entry for entry in database if entry['output'].startswith(prefix)]
    A.exact(len(entries),19,'complete production TU roster')
    by_source = {entry['file']:entry for entry in entries}
    A.exact(len(by_source),19,'unique producing sources')
    A.exact(set(by_source),{str(ROOT/name) for name in PRODUCERS},'exact producing source set')
    flags = ['-DWIREHAIR_BUILDING=1']+([] if mode=='native' else ['-DWIREHAIR_STATIC=1'])
    flags += ['-I'+str(ROOT/'include')]
    flags += {'native':['-O3','-DNDEBUG'], 'scalar':['-DANDROID','-O3','-DNDEBUG'],
              'asan':['-fsanitize=address,undefined','-fno-omit-frame-pointer','-g']}[mode]
    flags += ['-std=gnu++11','-fPIC','-Wall','-Wextra','-Wpedantic','-Werror']
    if mode=='asan': flags.append('-march=native')
    result = []
    for name in PRODUCERS:
        source = ROOT/name
        entry = by_source[str(source)]
        output = prefix+name+'.o'
        expected = ['/usr/bin/c++']+flags+['-o',output,'-c',str(source)]
        A.exact(set(entry),{'directory','command','file','output'},'CMake compile entry schema')
        A.exact((entry['directory'],entry['output']),(str(base),output),'actual producer directory/target')
        A.exact(shlex.split(entry['command']),expected,'exact production compiler recipe')
        result.append((source,base/output,flags))
    return result


def boundary_recipe(database, build, mode):
    A.exact(len(database), 8, 'eight qualified boundary/test translation units')
    selected = [r for r in database if r['output'].startswith('CMakeFiles/wh2_small_serialized.dir/')]
    A.exact(len(selected), 1, 'one producing boundary translation unit')
    entry = selected[0]
    A.exact(set(entry), {'directory', 'command', 'file', 'output'}, 'boundary compile schema')
    target = 'CMakeFiles/wh2_small_serialized.dir/home/catid/wirehair/bench/Wh2SmallSerialized.cpp.o'
    A.exact((entry['directory'], entry['file'], entry['output']),
            (str(build), str(HERE/'Wh2SmallSerialized.cpp'), target), 'exact boundary producer')
    flags = ['-DWH2_SMALL_CODEC_K=4']
    if mode == 'scalar':
        flags = ['-DANDROID']+flags+['-DWH2_SMALL_EXPECT_PORTABLE=1']
    flags += ['-I'+str(ROOT), '-I'+str(ROOT/'include'), '-I'+str(build)]
    flags += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-march=native', '-g'] if mode == 'asan' else ['-O3', '-DNDEBUG']
    flags += ['-std=c++11', '-fPIC', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-strict-aliasing', '-fno-lto']
    argv = ['/usr/bin/c++']+flags+['-o', target, '-c', str(HERE/'Wh2SmallSerialized.cpp')]
    A.exact(shlex.split(entry['command']), argv, 'exact qualified boundary flags')
    return build/target, argv


def context_size(raw):
    matches = re.findall(rb'^\S+\s+([0-9a-fA-F]+)\s+[Bb]\s+GF256Ctx$', raw, re.M)
    A.exact(len(matches), 1, 'unique sized GF256 context')
    size = int(matches[0], 16)
    A.require(size in GF_BYTES.values(), 'qualified GF context size')
    return size


def external_output(output):
    """Never add files inside the repository or any reused evidence root."""
    output = Path(output)
    protected = (ROOT, PRODUCTION, PROOFS, PRIOR, SMALL, SEALED, BOUNDARY_AUDIT.parent, Path('/var/tmp'))
    A.require(output.is_absolute() and '..' not in output.parts and not output.is_symlink() and
              output == output.resolve() and output not in (Path('/'), Path('/tmp')) and
              not any(base == output or base in output.parents for base in protected), 'external snapshot directory')
    return output


def qualified_inputs(mode, output):
    """Authenticate old producing evidence as data; never rebind changed paths."""
    A.require(mode in MODES, 'backend')
    output = external_output(output)
    dependencies, frozen, snapshots, historical = set(), {}, {}, []

    def retain(path, expected=None):
        path = Path(path)
        value = pin(path)
        if expected is not None:
            A.exact(value, expected, 'retained input identity: '+str(path))
        if path in frozen:
            A.exact(value, frozen[path], 'input changed during provenance inspection')
        dependencies.add(path); frozen[path] = value
        return value

    # Authenticate the closed library-producing proof through its original
    # scientific COMPLETE -> CLAIM -> manifest. Do not execute its old reader.
    complete_path, claim_path = PRIOR/'COMPLETE.json', PRIOR/'CLAIM.json'
    A.exact(retain(complete_path)['sha256'], PRIOR_COMPLETE_SHA, 'original library COMPLETE')
    complete = A.decode(A.read_regular(complete_path, 65536))
    A.exact((complete['protocol'], complete['outcome']), (PRIOR_PROTOCOL, 'PASS'), 'original library qualification')
    members = pin_map(complete['files'])
    A.exact(retain(claim_path, members[str(claim_path)])['sha256'], PRIOR_CLAIM_SHA, 'original library CLAIM')
    claim = A.decode(A.read_regular(claim_path, 1024**2))
    A.exact((claim['protocol'], claim['head']), (PRIOR_PROTOCOL, SOURCE_HEAD), 'original producing HEAD')
    claimed = pin_map(claim['pins'])
    manifest_path = PROOFS/mode/'manifest.json'
    A.exact(retain(manifest_path, claimed[str(manifest_path)])['sha256'],
            LIB_MANIFEST_SHA[mode], 'closed original producer manifest')
    manifest = A.decode(A.read_regular(manifest_path, 2*1024**2))
    A.exact((manifest['protocol'], manifest['mode']), (PRIOR_PROTOCOL, mode), 'original manifest backend')
    declared = pin_map(manifest['inputs']+manifest['artifacts'])
    for name, record in declared.items():
        A.exact(record, claimed[name], 'original manifest/receipt closure')
    # All historical installed tools and runtimes are checked BEFORE ar/git/
    # compiler inspection. Old non-producing repository docs need not be current.
    for name, record in declared.items():
        path = Path(name)
        if not any(base == path or base in path.parents for base in (ROOT, Path('/tmp'), Path('/var/tmp'))):
            retain(path, record)
    for name in ('c++', 'ar', 'nm', 'git', 'ldd'):
        A.require((Path('/usr/bin')/name).resolve(strict=True) in frozen, 'authenticated inspection tool')
    proof_path = PROOFS/mode/'qualified-library.json'
    retain(proof_path, declared[str(proof_path)])
    producer = A.decode(A.read_regular(proof_path, 2*1024**2))
    A.exact((producer['protocol'], producer['mode'], producer['producing_source_closure']),
            (PRIOR_PROTOCOL, mode, True), 'historically closed library producer')
    archive = PRODUCTION/mode/'libwirehair.a'
    for key in ('archive', 'reproduced', 'compile_database', 'qualification_log', 'neutral_inventory'):
        record = producer[key]
        A.exact(record, declared[record['path']], 'producer artifact binding')
        retain(Path(record['path']), record)
    A.exact(producer['archive']['path'], str(archive), 'actual production archive path')
    A.exact(frozen[archive]['sha256'], ARCHIVES[mode], 'actual WH1/WH2 archive identity')
    A.exact((producer['reproduced']['bytes'], producer['reproduced']['sha256']),
            (frozen[archive]['bytes'], ARCHIVES[mode]), 'whole reproduced archive identity')
    database_path = PRODUCTION/mode/'compile_commands.json'
    A.exact(producer['compile_database']['path'], str(database_path), 'original compile database path')
    recipes = producer_recipes(A.decode(A.read_regular(database_path, 2*1024**2)), mode)
    A.exact(command(['ar', 't', archive]).decode().splitlines(),
            [obj.name for _, obj, _ in recipes], 'complete nineteen-member archive order')
    A.exact(len(producer['members']), 19, 'all original producer records')
    producing = {}
    for (source, obj, _), record in zip(recipes, producer['members']):
        A.exact(set(record), {'source', 'object', 'reproduced', 'dependencies'}, 'original member schema')
        A.exact((record['source'], record['object']['path']), (str(source), str(obj)), 'original member source/object')
        for key in ('object', 'reproduced'):
            original = record[key]
            A.exact(original, declared[original['path']], 'original member pin binding')
            retain(Path(original['path']), original)
        A.exact((record['object']['bytes'], record['object']['sha256']),
                (record['reproduced']['bytes'], record['reproduced']['sha256']), 'original/reproduced object bytes')
        member_raw = command(['ar', 'p', archive, obj.name])
        A.exact((len(member_raw), A.sha(member_raw)),
                (record['object']['bytes'], record['object']['sha256']), 'actual archive member bytes')
        deps = pin_map(record['dependencies'])
        A.require(str(source) in deps, 'producing translation unit pinned')
        for name, value in deps.items():
            A.exact(value, declared[name], 'authenticated producing dependency')
            if name in producing: A.exact(producing[name], value, 'shared producing dependency')
            producing[name] = value
    for name, expected in sorted(producing.items()):
        path = Path(name)
        if path == ROOT/'codec/WirehairSmallCore.h':
            A.exact(expected['sha256'], CORE_OLD_SHA, 'original producing core identity')
            raw = command(['git', 'cat-file', 'blob', SOURCE_HEAD+':codec/WirehairSmallCore.h'])
            A.exact(dict(path=name, bytes=len(raw), sha256=A.sha(raw)), expected, 'exact original core snapshot')
            A.exact(retain(path)['sha256'], CORE_SHA, 'separately qualified current K4 core')
            snapshot = output/'historical-production-WirehairSmallCore.h'
            snapshots[snapshot] = raw
            historical.append(dict(original=expected, snapshot=dict(path=str(snapshot), bytes=len(raw), sha256=A.sha(raw))))
        else:
            retain(path, expected)
    # Authenticate the original preprocessor/compile/archive invocation slice
    # and its outputs. Later historical observer commands are out of scope.
    A.require(len(manifest['commands']) >= 40, 'complete original producing log')
    producing_commands = []
    for phase in range(2):
        for source, obj, flags in recipes:
            replica = PROOFS/mode/obj.name
            dep = Path(str(replica)+'.d')
            retain(dep, declared[str(dep)])
            recorded = preprocessor_dependencies(A.read_regular(dep, 2*1024**2), replica)
            member = producer['members'][len(producing_commands) % 19]
            A.exact(recorded, {Path(p['path']) for p in member['dependencies']}, 'original compiler dependency closure')
            tail = (['-M', '-MT', str(replica), '-MF', str(dep), str(source)] if phase == 0 else
                    ['-MD', '-MT', str(replica), '-MF', str(dep), '-o', str(replica), '-c', str(source)])
            producing_commands.append(dict(argv=['/usr/bin/c++']+flags+tail, cwd=str(PRODUCTION/mode)))
    reproduced_archive = producer['reproduced']['path']
    producing_commands += [dict(argv=['/usr/bin/ar', 'qc', reproduced_archive]+
                                [m['reproduced']['path'] for m in producer['members']], cwd=str(ROOT)),
                           dict(argv=['/usr/bin/ranlib', reproduced_archive], cwd=str(ROOT))]
    A.exact(manifest['commands'][:40], producing_commands, 'exact original producing commands')
    for index in range(40):
        for suffix in ('stdout', 'stderr', 'result.json'):
            path = PROOFS/mode/('command-%03d.%s' % (index, suffix))
            retain(path, declared[str(path)])
            if suffix == 'result.json':
                A.exact(A.decode(A.read_regular(path, 65536)), dict(returncode=0), 'successful original producing command')
            elif suffix == 'stderr':
                A.exact(A.read_regular(path, 65536), b'', 'clean original producing stderr')

    A.exact(retain(BOUNDARY_AUDIT)['sha256'], BOUNDARY_AUDIT_SHA, 'independent K4 boundary audit')
    audit = A.decode(A.read_regular(BOUNDARY_AUDIT, 65536))
    A.exact((audit['status'], audit['tests']), ('PASS', 105), 'qualified K4 boundary milestone')
    report_path = SMALL/mode/'RESULT.json'
    retain(report_path, audit['modes'][mode]['result'])
    report = A.decode(A.read_regular(report_path, 2*1024**2))
    A.exact((report['status'], report['tests']), ('PASS', 35), 'qualified backend tests')
    frozen_path = SMALL/mode/'FROZEN.json'
    retain(frozen_path, report['frozen'])
    original = A.decode(A.read_regular(frozen_path, 2*1024**2))
    A.exact((original['schema'], original['mode'], original['base_head']),
            ('wh2-k4-serialized-neutral-build-v1', mode, BOUNDARY_HEAD), 'original K4 boundary inputs')
    A.exact(original['library'], frozen[archive], 'matching production archive')
    source_pins = pin_map(original['source_pins'])
    build = SMALL/mode/'k4'
    artifacts = pin_map(report['artifacts'])
    for name, record in artifacts.items():
        if build in Path(name).parents: retain(Path(name), record)
    for step in report['steps']:
        if step['name'].startswith('k4-'):
            A.exact(step['returncode'], 0, 'retained K4 qualification command')
            for stream in ('stdout', 'stderr'):
                retain(Path(step[stream]['path']), step[stream])
    small_archive = build/'libwh2_small_serialized.a'
    fixture = build/'Wh2K4NativeData.inc'
    A.exact(frozen[fixture]['sha256'], FIXTURE_SHA, 'unchanged qualified K4 fixture')
    obj, recipe = boundary_recipe(A.decode(A.read_regular(build/'compile_commands.json', 1024**2)), build, mode)
    retain(obj)
    A.exact(command(['ar', 't', small_archive]).decode().splitlines(), [obj.name], 'sole K4 archive member')
    member_raw = command(['ar', 'p', small_archive, obj.name])
    A.exact((len(member_raw), A.sha(member_raw)), (frozen[obj]['bytes'], frozen[obj]['sha256']), 'K4 archive/member equality')
    depfile = Path(str(obj)+'.d'); retain(depfile)
    depraw = A.read_regular(depfile, 2*1024**2)
    boundary_deps = preprocessor_dependencies(depraw, obj.relative_to(build))
    A.require(HERE/'Wh2SmallSerialized.cpp' in boundary_deps and fixture in boundary_deps, 'boundary source/fixture dependency')
    for path in sorted(boundary_deps):
        if ROOT in path.parents:
            retain(path, source_pins[str(path)])
        elif path == fixture:
            retain(path, artifacts[str(path)])
        else:
            # Every installed boundary header is covered by the authenticated
            # production manifest; reproduction additionally proves this closure.
            retain(path, declared[str(path)])
    complete_path = SEALED/'COMPLETE.json'
    A.exact(retain(complete_path)['sha256'], COMPLETE_SHA, 'sealed retained K4 evidence')
    sealed = A.decode(A.read_regular(complete_path, 65536))
    A.exact(set(sealed['files']), {'CLAIM.json', 'raw.json', 'summary.json', 'stderr.txt'}, 'sealed structural evidence roster')
    for name, record in sealed['files'].items():
        retain(SEALED/name, dict(path=str(SEALED/name), **record))
    freeze_inputs(dependencies, frozen)
    provenance = dict(mode=mode, historical_production=producer, production_manifest=frozen[manifest_path],
        production_complete=frozen[PRIOR/'COMPLETE.json'], production_claim=frozen[PRIOR/'CLAIM.json'],
        serialized_audit=frozen[BOUNDARY_AUDIT], serialized_report=frozen[report_path],
        serialized_frozen=frozen[frozen_path], serialized_archive=frozen[small_archive],
        serialized_object=frozen[obj], serialized_command=recipe, serialized_dependencies=depraw.decode(),
        historical_snapshots=historical, input_pins=[frozen[p] for p in sorted(frozen)],
        scope='authenticated original production; K4 boundary still requires exact reproduction')
    return [small_archive, archive], dependencies, provenance, snapshots


def imported_files():
    return {Path(m.__file__).resolve(strict=True) for m in list(sys.modules.values())
            if getattr(m, '__file__', None) and Path(m.__file__).is_file()}


def build(mode, output):
    """Fresh observer build and bounded neutral checks, never scientific work."""
    A.require(mode in MODES, 'backend')
    output = Path(output)
    A.require(output.name == mode, 'mode-named build directory')
    output = output.parent.resolve(strict=True)/output.name
    external_output(output)
    A.require(ROOT not in output.parents and output != ROOT and not output.exists() and not output.is_symlink(), 'fresh external build')
    dependencies = {Path(__file__).resolve(), HERE/'Wh2K4SerializedRecoveryR0.py'}
    frozen = {}
    freeze_inputs(dependencies, frozen)
    spec = importlib.util.spec_from_file_location('_k4_recovery_reader_for_build', HERE/'Wh2K4SerializedRecoveryR0.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    A.exact((reader.PROTOCOL, reader.MODES), (PROTOCOL, MODES), 'new recovery reader contract')
    A.exact({k: os.environ.get(k) for k in reader.ENV_KEYS},
            dict({k: None for k in reader.ENV_KEYS}, **reader.SANITIZERS), 'ordinary allocator/sanitizer environment')
    dependencies.update(ROOT/name for name in reader.SOURCES)
    dependencies.update(imported_files())
    dependencies.add(Path(sys.executable).resolve(strict=True))
    freeze_inputs(dependencies, frozen)
    archives, reused, provenance, snapshots = qualified_inputs(mode, output)
    for record in provenance['input_pins']:
        path = Path(record['path'])
        A.exact(pin(path), record, 'validated producer input unchanged before use')
        if path in frozen:
            A.exact(frozen[path], record, 'shared input unchanged during provenance inspection')
        frozen[path] = record
    dependencies.update(reused)
    output.mkdir(mode=0o700)
    for path, raw in snapshots.items():
        A.publish(path, raw)
        dependencies.add(path)
    for record in provenance['historical_snapshots']:
        A.exact(pin(Path(record['snapshot']['path'])), record['snapshot'], 'published original bytes')
    fixture = SMALL/mode/'k4'/'Wh2K4NativeData.inc'
    A.publish(output/'Wh2K4NativeData.inc', A.read_regular(fixture, 2*1024**2))
    A.exact(pin(output/'Wh2K4NativeData.inc')['sha256'], FIXTURE_SHA, 'copied retained fixture')
    dependencies.add(output/'Wh2K4NativeData.inc')
    for name in ('c++', 'cc', 'as', 'ld', 'nm', 'ar', 'ranlib', 'ldd', 'bash', 'git'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    freeze_inputs(dependencies, frozen)
    for name in ('cc1', 'cc1plus', 'collect2'):
        dependencies.add(Path(command(['c++', '-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    for name in LINK_INPUTS:
        path = Path(command(['c++', '-print-file-name='+name]).decode().strip())
        A.require(path.is_absolute(), 'resolved compiler/linker input')
        dependencies.add(path.resolve(strict=True))
    elf = [p for p in dependencies if os.access(str(p), os.X_OK) and
           A.read_regular(p, 256*1024**2, installed=not any(base in p.parents for base in
                          (ROOT, Path('/tmp'), Path('/var/tmp'))))[:4] == b'\x7fELF']
    for path in elf:
        linked = command(['ldd', path])
        A.require(b'not found' not in linked, 'resolved tool runtimes')
        dependencies.update(Path(p).resolve(strict=True) for p in linked.decode().split() if p.startswith('/'))
    freeze_inputs(dependencies, frozen)
    commands = []

    def run(argv, expected=0, clean_stderr=True, cwd=ROOT):
        argv = list(map(str, argv))
        index = len(commands)
        commands.append(dict(argv=argv, cwd=str(cwd)))
        freeze_inputs(dependencies, frozen)
        try:
            result = subprocess.run(argv, cwd=cwd, env=process_environment(), stdin=subprocess.DEVNULL,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        except subprocess.TimeoutExpired as error:
            A.publish(output/('command-%03d.stdout' % index), error.stdout or b'')
            A.publish(output/('command-%03d.stderr' % index), error.stderr or b'')
            A.publish(output/('command-%03d.result.json' % index), A.canonical(dict(timeout=True)))
            raise
        A.publish(output/('command-%03d.stdout' % index), result.stdout)
        A.publish(output/('command-%03d.stderr' % index), result.stderr)
        A.publish(output/('command-%03d.result.json' % index), A.canonical(dict(returncode=result.returncode)))
        A.require(len(result.stdout) <= 16*1024**2 and len(result.stderr) <= 65536, 'build/neutral output cap')
        A.require(result.returncode == expected and (not clean_stderr or not result.stderr), 'build/neutral command failed: '+str(argv))
        freeze_inputs(dependencies, frozen)
        return result

    # Close the boundary's installed-header gap by reproducing ONLY its one
    # object and archive. Preserve original compiler cwd and object basename;
    # all writes go to this fresh output, never the qualification directory.
    boundary_dir = SMALL/mode/'k4'
    original_obj, argv = boundary_recipe(
        A.decode(A.read_regular(boundary_dir/'compile_commands.json', 1024**2)), boundary_dir, mode)
    source, boundary_flags = HERE/'Wh2SmallSerialized.cpp', argv[1:-4]
    obj, dep = output/original_obj.name, output/(original_obj.name+'.d')
    run(['/usr/bin/c++']+boundary_flags+['-M', '-MT', str(obj), '-MF', str(dep), str(source)], cwd=boundary_dir)
    before = preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj)
    A.exact(before, preprocessor_dependencies(provenance['serialized_dependencies'].encode(),
            original_obj.relative_to(boundary_dir)), 'exact original/reproducing boundary dependencies')
    dependencies.update(before); freeze_inputs(dependencies, frozen)
    run(['/usr/bin/c++']+boundary_flags+['-MD', '-MT', str(obj), '-MF', str(dep),
                                     '-o', str(obj), '-c', str(source)], cwd=boundary_dir)
    A.exact(preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj), before, 'boundary compiler closure')
    A.exact(A.read_regular(obj, 64*1024**2), A.read_regular(original_obj, 64*1024**2), 'byte-identical boundary object')
    dependencies.update((obj, dep)); freeze_inputs(dependencies, frozen)
    reproduced = output/'libwh2_small_serialized-reproduced.a'
    run(['/usr/bin/ar', 'qc', reproduced, obj]); run(['/usr/bin/ranlib', reproduced])
    A.exact(A.read_regular(reproduced, 128*1024**2), A.read_regular(archives[0], 128*1024**2), 'byte-identical boundary archive')
    dependencies.update((obj, dep, reproduced)); freeze_inputs(dependencies, frozen)
    provenance.update(serialized_reproduced_object=pin(obj), serialized_reproduced_archive=pin(reproduced),
        serialized_reproducing_dependencies=[frozen[p] for p in sorted(before)],
        producing_source_closure=True,
        scope='authenticated original nineteen-object production; exact newly reproduced K4 boundary')
    A.publish(output/'qualified-library.json', A.canonical(provenance))
    dependencies.add(output/'qualified-library.json'); freeze_inputs(dependencies, frozen)

    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_SMALL_CODEC_K=4', '-DWH2_SMALL_RECOVERY_BACKEND='+str(MODES.index(mode)),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += shlex.split(ASAN_FLAGS) if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags.append('-DANDROID=1')
    objects = []
    for source in (HERE/'Wh2K4SerializedRecoveryR0.cpp', HERE/'Wh2FrozenTrace.cpp'):
        obj, dep = output/(source.stem+'.o'), output/(source.stem+'.d')
        run(['/usr/bin/c++']+flags+['-M', '-MT', str(obj), '-MF', str(dep), str(source)])
        before = preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj)
        dependencies.update(before)
        freeze_inputs(dependencies, frozen)
        run(['/usr/bin/c++']+flags+['-MD', '-MF', str(dep), '-c', str(source), '-o', str(obj)])
        A.exact(preprocessor_dependencies(A.read_regular(dep, 2*1024**2), obj), before, 'actual compile dependency closure')
        objects.append(obj)
        dependencies.add(obj)
        freeze_inputs(dependencies, frozen)
    executable = output/'recovery_worker'
    argv = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan':
        argv.append('-fsanitize=address,undefined')
    run(argv+list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)])
    loads = [line[5:] for line in A.read_regular(output/'link.map', 4*1024**2).decode().splitlines() if line.startswith('LOAD ')]
    A.require(loads and all(Path(p).is_absolute() for p in loads), 'absolute linker input roster')
    loaded = {Path(p).resolve(strict=True) for p in loads}
    A.require(loaded <= dependencies, 'all linker inputs pinned before link')
    symbols = command(['nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    expected = ['GF256Ctx', 'gf256_init_', 'wirehair_encoder_create', 'wirehair_encoder_create_ex',
                'wirehair_encoder_create_owned', 'wirehair_encoder_create_owned_ex',
                'wirehair_decoder_create', 'wirehair_encode', 'wirehair_decode', 'wirehair_recover', 'wirehair_free']
    expected += ['wh2_small_'+name for name in ('profile_validate', 'encoder_create', 'encoder_create_profile',
                 'decoder_create', 'encoder_detach_input', 'encode', 'decode', 'recover', 'free')]
    expected += ['wirehair_v2_'+name for name in ('encoder_create', 'encoder_create_with_options',
                 'encoder_create_profile_id_with_options', 'decoder_create', 'encode', 'decode', 'recover', 'free')]
    A.require(all(names.count(name) == 1 for name in expected), 'one exact private/ordinary/WH1 API and GF runtime')
    size = context_size(command(['nm', '-S', '--defined-only', executable]))
    A.exact(size, GF_BYTES[mode], 'matching private GF ABI')
    dependencies.add(executable)
    elf = [p for p in dependencies if os.access(str(p), os.X_OK) and
           A.read_regular(p, 256*1024**2, installed=ROOT not in p.parents and output not in p.parents)[:4] == b'\x7fELF']
    for path in elf:
        linked = command(['ldd', path])
        A.require(b'not found' not in linked, 'resolved runtime libraries')
        dependencies.update(Path(p).resolve(strict=True) for p in linked.decode().split() if p.startswith('/'))
    dependencies.update(imported_files())
    freeze_inputs(dependencies, frozen)
    contract = run([executable, '--contract']).stdout
    A.exact(A.decode(contract), dict(K=4, arms=6, records=6254, cpu_seconds=180, wall_seconds=210,
        address_space_mib=384, asan_shadow_exempt=mode == 'asan', backend=mode,
        claim_path='/var/tmp/wh2-k4-serialized-recovery-r0/CLAIM.json'), 'compiled neutral/resource contract')
    A.publish(output/'contract.json', contract)
    claim_path = output/'neutral-claim.json'
    claim = A.canonical(dict(protocol=PROTOCOL, neutral=True))
    A.publish(claim_path, claim)
    A.exact(run([executable, '--neutral-claim', claim_path, A.sha(claim)]).stdout,
            b'PASS claim authentication\n', 'positive neutral claim authentication')
    negatives = []
    for argv in ([executable, '--neutral-claim', claim_path, '0'*64],
                 [executable, '--neutral-claim', output/'absent', A.sha(claim)],
                 [executable, '--neutral-claim', claim_path, 'g'*64],
                 [executable], [executable, '--unknown'], [executable, '--worker', '0']):
        result = run(argv, expected=1, clean_stderr=False)
        A.require(not result.stdout and result.stderr.startswith(b'INVALID:'), 'negative CLI/claim rejection')
        negatives.append(dict(argv=list(map(str, argv)), code=result.returncode, stderr=result.stderr.decode()))
    A.publish(output/'negative-cli.json', A.canonical(negatives))
    neutral = run([executable, '--neutral']).stdout
    A.exact(neutral, b'PASS 48 neutral cases, six ownership-matched APIs, all-arm packet/rank oracle, twelve late-call cleanup checks\n', 'neutral cleanup/API gates')
    A.publish(output/'neutral.txt', neutral)
    fixture_raw = run([executable, '--neutral-fixtures']).stdout
    reader.verify(fixture_raw, '0'*64, mode, neutral=True)
    A.publish(output/'fixtures.jsonl', fixture_raw)
    dependencies.update(imported_files())
    freeze_inputs(dependencies, frozen)
    # Re-check the exact archives and all current/historical inputs after the
    # final observer execution, before publication of the new build receipt.
    artifacts = [pin(p) for p in sorted(output.iterdir())]
    freeze_inputs(dependencies, frozen)
    manifest = dict(protocol=PROTOCOL, mode=mode, commands=commands,
                    inputs=[frozen[p] for p in sorted(dependencies)], artifacts=artifacts)
    A.publish(output/'manifest.json', A.canonical(manifest))
    return manifest
