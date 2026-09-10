#!/usr/bin/env python3
"""Build only the new K8 recovery observer; reuse closed, unchanged archives.

No historical controller/verifier is imported or invoked. Producing recipes
and dependency records are parsed as data, never executed. The sole scientific
launch belongs to Wh2K8SerializedRecoveryR0.py, not this module.
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
SPEC = importlib.util.spec_from_file_location('_k8_recovery_build_io', HERE/'Wh2AlignedIntermediateCostR0.py')
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
MODES = ('native', 'scalar', 'asan')
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
PROTOCOL = 'wirehair.wh2.k8-serialized-recovery-r0'
PRODUCTION = Path('/tmp/wh2-v2-k5-admission.VR7a0QAY')
PROOFS = Path('/tmp/wh2-k5-public-recovery-qualified.bksSezzu')
SMALL = Path('/tmp/wh2-k8-serialized.qQHyZqaa')
SEALED = Path('/var/tmp/wh2-k8-thue-morse-r0')
SOURCE_HEAD = '421626d0883e11560dac1f635a2ccf68710b30e4'
QUALIFIED_HEAD = 'ee160033345894585db6c5312b4567d96126d9d7'
LIB_DIRS = dict(native='native-default', scalar='scalar', asan='asan')
LIB_MANIFEST_SHA = dict(
    native='bb89d5d4551989dbd02d7128ca8bff12b02c135e08ea9e6d3f4667ae19eb4bff',
    scalar='2eaaa3df1dc30239fa970076b3cc63c8b4ac32996522f89920c8c8efbcedee5b',
    asan='4ae747bcc734699d704a871e9c3a58a186ea58fcdccdc6e23bfa6f56d01c704b')
ARCHIVES = dict(
    native='ddf322a6798676cdfbf29bcf05744a2d70dc7a0449d81947eda60301e4e62464',
    scalar='a4961ebe7a056ccb0b82655e32fe696773f52865993214e08d54b2eac8ea4a91',
    asan='44086d6961a27578c1f15ef3f18ed03df9ba182768576ce0795e3597d9dc52ab')
SMALL_ARCHIVES = dict(
    native='0ad76d0471b5b3fda81f40bec87e52e328f4cb6c0dca5f6dbdd90f4e84d583c6',
    scalar='4313edd59c522fd11b0a24b9b5042fe0df036cf45eaae204ad463b9b80dd1a10',
    asan='6f1e29b12999a867fb909c3d017b95c97a3978955bfe5f75e2415bb791341c6b')
QUALIFIED_SHA = '991c5a69b5fa22d34ad0d51922e81dd64a678a867f38add7d121c05b342bd1a4'
CORE_OLD_SHA = '00c299e130df9fbe3acbc920e3f7c94dc72dcd0b4e581dd6dd337735c242da40'
CORE_SHA = '3f534e7b2b127860e5497d9af3d090fe97d300346c116b37cbffea9f7b0d6769'
README_SHA = 'd1b4b923e27adb54e777dda2bb22c870d1b3ef58e9ee05cdc70a26fd1295b3c0'
FIXTURE_SHA = '99e2b20da1405755136311d58bd7391474866bdee1a10ccf22a2de7482b4aaea'
COMPLETE_SHA = '4911bedbb288c8a7e39577c35dbe6c31a0a2990e0d4dfab1ee4466a32f99fdb0'
GF_BYTES = dict(native=141328, scalar=137232, asan=157728)
ASAN_FLAGS = '-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -march=native'
PRODUCER_SOURCES = (
    'wirehair.cpp', 'codec/WirehairSmall.cpp', 'codec/WirehairSmallK5.cpp',
    'codec/WirehairK6.cpp', 'codec/WirehairK6Core.cpp', 'gf256.cpp',
    'WirehairCodec.cpp', 'WirehairTools.cpp', 'codec/WirehairV2Codec.cpp',
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


def library_commands(raw, mode):
    A.require(mode in MODES, 'backend')
    base = PRODUCTION/LIB_DIRS[mode]
    prefix = Path('CMakeFiles/wirehair_objects.dir' if mode == 'native' else 'CMakeFiles/wirehair.dir')
    defines = ['-DWIREHAIR_BUILDING=1']+([] if mode == 'native' else ['-DWIREHAIR_STATIC=1'])
    flags = {
        'native': '-O3 -DNDEBUG -std=gnu++11 -fPIC -Wall -Wextra -Wpedantic -Werror',
        'scalar': '-DANDROID=1 -DWIREHAIR_SMALL_EXPECT_PORTABLE=1 -O3 -DNDEBUG -std=gnu++11 -fPIC -Wall -Wextra -Wpedantic -Werror',
        'asan': '-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -g -std=gnu++11 -fPIC -Wall -Wextra -Wpedantic -Werror -march=native'}
    lines = raw.decode().splitlines()
    A.exact(len(lines), 19, 'complete eighteen-object producer commands')
    objects = []
    for line, name in zip(lines, PRODUCER_SOURCES):
        obj = prefix/(name+'.o')
        expected = ['/usr/bin/c++']+defines+['-I'+str(ROOT/'include')]+shlex.split(flags[mode])
        expected += ['-MD', '-MT', str(obj), '-MF', str(obj)+'.d', '-o', str(obj), '-c', str(ROOT/name)]
        A.exact(shlex.split(line), expected, 'exact historical producer command')
        objects.append(base/obj)
    archive = [':', '&&', '/usr/bin/cmake', '-E', 'rm', '-f', 'libwirehair.a', '&&',
               '/usr/bin/ar', 'qc', 'libwirehair.a']+[str(o.relative_to(base)) for o in objects]
    archive += ['&&', '/usr/bin/ranlib', 'libwirehair.a', '&&', ':']
    A.exact(shlex.split(lines[-1]), archive, 'exact historical archive recipe')
    return objects


def library_dependencies(raw, relative, source, object_mtime):
    lines = raw.decode().splitlines()
    A.require(lines, 'nonempty original dependency record')
    match = re.fullmatch(re.escape(str(relative))+r': #deps ([0-9]+), deps mtime ([0-9]+) \(VALID\)', lines[0])
    A.require(match is not None and lines[-1] == '', 'complete original dependency record')
    A.exact(int(match[2]), object_mtime, 'producing object mtime')
    A.exact(len(lines)-2, int(match[1]), 'original dependency multiplicity')
    A.require(all(p.startswith('    /') for p in lines[1:-1]), 'original absolute dependencies')
    originals = [Path(p[4:]) for p in lines[1:-1]]
    result = {p.resolve(strict=True) for p in originals}
    A.require(source in result, 'producing source in dependency record')
    A.require(all(p == p.resolve(strict=True) for p in originals if ROOT in p.parents), 'historical source redirect')
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
    flags = ['-DWH2_SMALL_CODEC_K=8']
    if mode == 'scalar':
        flags = ['-DANDROID']+flags+['-DWH2_SMALL_EXPECT_PORTABLE=1']
    flags += ['-I'+str(ROOT), '-I'+str(ROOT/'include'), '-I'+str(build)]
    flags += shlex.split(ASAN_FLAGS)+['-g'] if mode == 'asan' else ['-O3', '-DNDEBUG']
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


def qualified_inputs(mode, output):
    """Read-only reuse: archives, current paths, proof, output-local snapshots.

    Snapshots are returned as bytes, not written here. In particular the old
    production header is never rebound to the current K8 header at its path.
    """
    A.require(mode in MODES, 'backend')
    output = Path(output)
    A.require(output.is_absolute() and ROOT not in output.parents and output != ROOT and
              '..' not in output.parts and not output.is_symlink(), 'external snapshot directory')
    dependencies, frozen, snapshots = set(), {}, {}

    def retain(path, expected=None):
        path = Path(path)
        value = pin(path)
        if expected is not None:
            A.exact(value, expected, 'retained input identity: '+str(path))
        if path in frozen:
            A.exact(value, frozen[path], 'input changed during provenance inspection')
        dependencies.add(path)
        frozen[path] = value
        return value

    manifest_path = PROOFS/mode/'manifest.json'
    A.exact(retain(manifest_path)['sha256'], LIB_MANIFEST_SHA[mode], 'closed original producer manifest')
    manifest = A.decode(A.read_regular(manifest_path, 2*1024**2))
    declared = pin_map(manifest['inputs']+manifest['artifacts'])
    # Authenticate inspection tools/runtime before using ar or git, not merely
    # after they have supplied the archive/member/source evidence.
    for path_string, expected in declared.items():
        path = Path(path_string)
        if ROOT not in path.parents and Path('/tmp') not in path.parents and Path('/var/tmp') not in path.parents:
            retain(path, expected)
    proof_path = PROOFS/mode/'qualified-library.json'
    retain(proof_path, declared[str(proof_path)])
    producer = A.decode(A.read_regular(proof_path, 2*1024**2))
    A.exact((producer['source_head'], producer['mode'], producer['library_source_provenance_closed']),
            (SOURCE_HEAD, mode, True), 'historically closed library producer')
    base = PRODUCTION/LIB_DIRS[mode]
    archive = base/'libwirehair.a'
    A.exact(retain(archive, producer['archive'])['sha256'], ARCHIVES[mode], 'exact actual WH1/WH2 archive')
    objects = library_commands(producer['commands'].encode(), mode)
    A.exact(command(['ar', 't', archive]).decode().splitlines(), [o.name for o in objects], 'complete archive roster')
    A.exact(len(producer['members']), 18, 'all original producer records')
    producing = set()
    for obj, name, record in zip(objects, PRODUCER_SOURCES, producer['members']):
        retain(obj, record['object'])
        A.exact(record['source'], str(ROOT/name), 'original member source')
        A.exact(A.sha(command(['ar', 'p', archive, obj.name])), frozen[obj]['sha256'], 'member/object equality')
        relative = obj.relative_to(base)
        paths = library_dependencies(record['dependencies'].encode(), relative, ROOT/name, obj.stat().st_mtime_ns)
        A.exact(paths, preprocessor_dependencies(record['preprocessor_dependencies'].encode(), relative),
                'retained independent original dependency agreement')
        producing.update(paths)
    historical = []
    for path in sorted(producing):
        expected = declared[str(path)]
        if path == ROOT/'codec/WirehairSmallCore.h':
            A.exact(expected['sha256'], CORE_OLD_SHA, 'original producing core identity')
            raw = command(['git', 'cat-file', 'blob', SOURCE_HEAD+':codec/WirehairSmallCore.h'])
            A.exact(dict(path=str(path), bytes=len(raw), sha256=A.sha(raw)), expected, 'exact original core snapshot')
            A.exact(retain(path)['sha256'], CORE_SHA, 'separately qualified current K8 core')
            snapshot = output/'historical-production-WirehairSmallCore.h'
            snapshots[snapshot] = raw
            historical.append(dict(original=expected, snapshot=dict(path=str(snapshot), bytes=len(raw), sha256=A.sha(raw))))
        else:
            retain(path, expected)
    # Preserve producer metadata in addition to the installed pins above.
    # Do not import unrelated old scientific source/controllers into new closure.
    for path_string, expected in declared.items():
        path = Path(path_string)
        if base in path.parents or path in (ROOT/'CMakeLists.txt', ROOT/'abi/wirehair.map'):
            retain(path, expected)
    retain(Path(producer['qualification_log']['path']), producer['qualification_log'])

    qualified_path = SMALL/'QUALIFIED.json'
    A.exact(retain(qualified_path)['sha256'], QUALIFIED_SHA, 'K8 serialized qualification')
    qualification = A.decode(A.read_regular(qualified_path, 2*1024**2))
    A.exact((qualification['outcome'], qualification['runtime_tests']), ('QUALIFIED_SERIALIZED', 84), 'K8 correctness milestone')
    for path, record in qualification['files'].items():
        retain(Path(path), dict(path=path, **record))
    report_path = SMALL/(mode+'-qualification.json')
    report = A.decode(A.read_regular(report_path, 2*1024**2))
    A.exact((report['source_head'], report['mode'], report['tests'], report['outcome']),
            (QUALIFIED_HEAD, mode, 28, 'PASS'), 'qualified K8 backend')
    A.exact(report['source_hashes'], qualification['source_hashes'], 'same frozen source roster')
    A.exact((report['library'], report['library_sha256']), (str(archive), ARCHIVES[mode]), 'matching production input')
    A.exact(report['driver_sha256'], frozen[SMALL/'qualify.py']['sha256'], 'qualified build driver')
    for key in ('compiler', 'interpreter'):
        A.exact(retain(Path(report[key]))['sha256'], report[key+'_sha256'], 'qualified '+key)
    for name, digest in report['source_hashes'].items():
        path = ROOT/name
        if name == 'bench/Wh2SmallNative/README.md':
            old = SMALL/'README-pre-result.md'
            A.exact(retain(old)['sha256'], README_SHA, 'qualified nonproducing documentation snapshot')
            A.exact(digest, README_SHA, 'original README hash remains bound')
            snapshot = output/'historical-qualification-README.md'
            snapshots[snapshot] = A.read_regular(old, 2*1024**2)
            historical.append(dict(original=dict(path=str(path), bytes=len(snapshots[snapshot]), sha256=digest),
                                   snapshot=dict(path=str(snapshot), bytes=len(snapshots[snapshot]), sha256=digest)))
        else:
            A.exact(retain(path)['sha256'], digest, 'qualified current boundary source')
    for path, record in report['artifacts'].items():
        retain(Path(path), dict(path=path, **record))
    for record in report['commands']:
        A.exact(record['returncode'], 0, 'retained successful qualification command')
        for stream in ('stdout', 'stderr'):
            path = SMALL/(mode+'-'+record['name']+'.'+stream)
            A.exact(retain(path)['sha256'], record[stream+'_sha256'], 'qualified command stream')
    build = SMALL/(mode+'-k8')
    small_archive = build/'libwh2_small_serialized.a'
    A.exact(frozen[small_archive]['sha256'], SMALL_ARCHIVES[mode], 'exact K8 boundary archive')
    fixture = build/'Wh2K8NativeData.inc'
    A.exact(frozen[fixture]['sha256'], FIXTURE_SHA, 'unchanged generated K8 fixture')
    obj, recipe = boundary_recipe(A.decode(A.read_regular(build/'compile_commands.json', 1024**2)), build, mode)
    retain(obj)
    A.exact(command(['ar', 't', small_archive]).decode().splitlines(), [obj.name], 'sole K8 archive member')
    A.exact(A.sha(command(['ar', 'p', small_archive, obj.name])), frozen[obj]['sha256'], 'K8 archive/member equality')
    depfile = Path(str(obj)+'.d')
    depraw = A.read_regular(depfile, 2*1024**2)
    boundary_deps = preprocessor_dependencies(depraw, obj.relative_to(build))
    A.require(HERE/'Wh2SmallSerialized.cpp' in boundary_deps and fixture in boundary_deps, 'boundary source/fixture dependency')
    for path in sorted(boundary_deps):
        if ROOT in path.parents or path == fixture:
            A.require(path in frozen, 'every boundary producing repository/fixture pin')
        else:
            retain(path, declared[str(path)])
    # The old producer manifest covers every producing system dependency of
    # the K8 archive, not every system dependency of its separate test binaries.
    complete = SEALED/'COMPLETE.json'
    A.exact(retain(complete)['sha256'], COMPLETE_SHA, 'sealed retained K8 evidence')
    sealed = A.decode(A.read_regular(complete, 65536))
    A.exact(set(sealed['files']), {'CLAIM.json', 'raw.json', 'summary.json', 'stderr.txt'}, 'sealed evidence roster')
    for name, record in sealed['files'].items():
        retain(SEALED/name, dict(path=str(SEALED/name), **record))
    freeze_inputs(dependencies, frozen)
    provenance = dict(mode=mode, historical_production=producer, production_manifest=frozen[manifest_path],
        serialized_qualification=frozen[qualified_path], serialized_report=frozen[report_path],
        serialized_archive=frozen[small_archive], serialized_object=frozen[obj],
        serialized_command=recipe, serialized_dependencies=depraw.decode(),
        historical_snapshots=historical, input_pins=[frozen[p] for p in sorted(frozen)],
        scope='exact retained archives; original source identities preserved, no old verifier rerun')
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
    A.require(ROOT not in output.parents and output != ROOT and not output.exists() and not output.is_symlink(), 'fresh external build')
    dependencies = {Path(__file__).resolve(), HERE/'Wh2K8SerializedRecoveryR0.py'}
    frozen = {}
    freeze_inputs(dependencies, frozen)
    spec = importlib.util.spec_from_file_location('_k8_recovery_reader_for_build', HERE/'Wh2K8SerializedRecoveryR0.py')
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
    A.publish(output/'qualified-library.json', A.canonical(provenance))
    fixture = SMALL/(mode+'-k8')/'Wh2K8NativeData.inc'
    A.publish(output/'Wh2K8NativeData.inc', A.read_regular(fixture, 2*1024**2))
    A.exact(pin(output/'Wh2K8NativeData.inc')['sha256'], FIXTURE_SHA, 'copied retained fixture')
    dependencies.update((output/'Wh2K8NativeData.inc', output/'qualified-library.json'))
    for name in ('c++', 'cc', 'as', 'ld', 'nm', 'ar', 'ranlib', 'ldd', 'bash', 'git'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    for name in ('cc1', 'cc1plus', 'collect2'):
        dependencies.add(Path(command(['c++', '-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    freeze_inputs(dependencies, frozen)
    commands = []

    def run(argv, expected=0, clean_stderr=True):
        argv = list(map(str, argv))
        index = len(commands)
        commands.append(argv)
        freeze_inputs(dependencies, frozen)
        result = subprocess.run(argv, cwd=ROOT, env=process_environment(), stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
        A.publish(output/('command-%03d.stdout' % index), result.stdout)
        A.publish(output/('command-%03d.stderr' % index), result.stderr)
        A.publish(output/('command-%03d.result.json' % index), A.canonical(dict(returncode=result.returncode)))
        A.require(len(result.stdout) <= 16*1024**2 and len(result.stderr) <= 65536, 'build/neutral output cap')
        A.require(result.returncode == expected and (not clean_stderr or not result.stderr), 'build/neutral command failed: '+str(argv))
        freeze_inputs(dependencies, frozen)
        return result

    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_SMALL_CODEC_K=8', '-DWH2_SMALL_RECOVERY_BACKEND='+str(MODES.index(mode)),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += shlex.split(ASAN_FLAGS) if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags.append('-DANDROID=1')
    objects = []
    for source in (HERE/'Wh2K8SerializedRecoveryR0.cpp', HERE/'Wh2FrozenTrace.cpp'):
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
    A.exact(A.decode(contract), dict(K=8, arms=6, records=6260, cpu_seconds=180, wall_seconds=210,
        address_space_mib=384, asan_shadow_exempt=mode == 'asan', backend=mode,
        claim_path='/var/tmp/wh2-k8-serialized-recovery-r0/CLAIM.json'), 'compiled neutral/resource contract')
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
