#!/usr/bin/env python3
"""Build the installed WHV2 K8 timing observer; reproduce all archive objects first.

Existing qualification directories and spent controllers are never modified or
executed. Original CMake compile recipes are data; their output paths alone are
redirected into a fresh directory. Byte identity is required before archive reuse.
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
SPEC = importlib.util.spec_from_file_location('_k8_public_build_io', HERE/'Wh2AlignedIntermediateCostR0.py')
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
MODES = ('native', 'scalar', 'asan')
PROTOCOL = 'wirehair.wh2.k8-public-cost-r0'
ENV_KEYS = ('MALLOC_TRIM_THRESHOLD_', 'MALLOC_MMAP_THRESHOLD_', 'MALLOC_TOP_PAD_',
            'MALLOC_PERTURB_', 'GLIBC_TUNABLES', 'LD_PRELOAD', 'LD_LIBRARY_PATH')
PRODUCTION = Path('/tmp/wh2-v2-k8-admission.YNN2XmAx')
INVENTORY_SHA = '61535c7d6611de03095e55b906c8e4b4096148bbe743327b9039c554957b7285'
ARCHIVES = dict(native='0d4a0a7ec3cfa4b8e0f08df75e2c12935b51f28653c8c953a32bb49a5658eff5',
                scalar='bd7eb2eeb0eeaea15e11d131c6fba2921bf44cb1526dea8dff86171b496de867',
                asan='a7043eb65eb864aefd210f51b360d4591a365c60710c3aec87baf0abac3d4f1a')
DATABASES = dict(native='1d51beb1861ea024c72c1ea81600cdb00192ffcc2e1986781bb9399fae27f5a2',
                 scalar='8162ca0b09d255514818274f48cc08d8b710977fb89228ba3051247d53307c23',
                 asan='84d8efd846114d4b65b89077918553d55d4b3d2e17d062f33a3037cb2d3753a1')
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
OBSERVERS = tuple(HERE/name for name in ('Wh2K8PublicCostR0.cpp', 'Wh2FrozenTrace.cpp',
    'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp'))

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



def imported_files():
    return {Path(m.__file__).resolve(strict=True) for m in list(sys.modules.values())
            if getattr(m, '__file__', None) and Path(m.__file__).is_file()}


def context_size(raw):
    matches = re.findall(rb'^\S+\s+([0-9a-fA-F]+)\s+[Bb]\s+GF256Ctx$', raw, re.M)
    A.exact(len(matches), 1, 'unique sized GF256 context')
    return int(matches[0], 16)


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


def build(mode, output):
    """Reproduce installed archive and qualify bounded new observer; no science."""
    A.require(mode in MODES, 'backend')
    output = Path(output)
    A.require(output.name==mode,'mode-named build directory')
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT not in output.parents and output!=ROOT and not output.exists() and
              not output.is_symlink(),'fresh external build')
    spec = importlib.util.spec_from_file_location('_k8_public_reader_for_build',HERE/'Wh2K8PublicCostR0.py')
    reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)
    A.exact((reader.PROTOCOL,reader.MODES),(PROTOCOL,MODES),'reader identity')
    A.exact({k:os.environ.get(k) for k in reader.ENV_KEYS},
            dict({k:None for k in reader.ENV_KEYS},**reader.SANITIZERS),'allocator/sanitizer environment')
    dependencies = {ROOT/name for name in reader.SOURCES}|set(OBSERVERS)|imported_files()
    frozen, commands = {}, []
    inventory_path = PRODUCTION/'NEUTRAL_QUALIFIED.json'
    inventory_bytes = A.read_regular(inventory_path,4*1024**2)
    A.exact(A.sha(inventory_bytes),INVENTORY_SHA,'original neutral inventory')
    inventory = pin_map(A.decode(inventory_bytes)['files'])
    database_path = PRODUCTION/mode/'compile_commands.json'
    original_archive = PRODUCTION/mode/'libwirehair.a'
    qualification_log = PRODUCTION/mode/'final-qualification.log'
    for path in (ROOT/'CMakeLists.txt',database_path,original_archive,qualification_log):
        A.exact(pin(path),inventory[str(path)],'retained neutral input')
        dependencies.add(path)
    A.exact(pin(database_path)['sha256'],DATABASES[mode],'fixed compile database')
    A.exact(pin(original_archive)['sha256'],ARCHIVES[mode],'actual installed archive')
    recipes = producer_recipes(A.decode(A.read_regular(database_path,8*1024**2)),mode)
    dependencies.add(inventory_path)
    dependencies.add(Path(sys.executable).resolve(strict=True))
    for name in ('c++','cc','as','ld','nm','ar','ranlib','ldd','bash','git'):
        dependencies.add((Path('/usr/bin')/name).resolve(strict=True))
    freeze_inputs(dependencies,frozen)
    for name in ('cc1','cc1plus','collect2'):
        dependencies.add(Path(command(['c++','-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    # Resolve all expected startup/link inputs before any archive/tool link.
    for name in LINK_INPUTS:
        path = Path(command(['c++','-print-file-name='+name]).decode().strip())
        A.require(path.is_absolute(),'resolved compiler/linker input')
        dependencies.add(path.resolve(strict=True))
    elf = [p for p in dependencies if os.access(str(p),os.X_OK) and
           A.read_regular(p,256*1024**2,installed=ROOT not in p.parents and Path('/tmp') not in p.parents
                          and Path('/var/tmp') not in p.parents)[:4]==b'\x7fELF']
    for path in elf:
        linked = command(['ldd',path])
        A.require(b'not found' not in linked,'resolved tool runtimes')
        dependencies.update(Path(p).resolve(strict=True) for p in linked.decode().split() if p.startswith('/'))
    dependencies.update(imported_files())
    freeze_inputs(dependencies,frozen)
    output.mkdir(mode=0o700)

    def run(argv, expected=0, clean_stderr=True, cwd=ROOT):
        argv = list(map(str,argv)); index = len(commands)
        commands.append(dict(argv=argv,cwd=str(cwd)))
        freeze_inputs(dependencies,frozen)
        result = subprocess.run(argv,cwd=cwd,env=process_environment(),stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
        A.publish(output/('command-%03d.stdout'%index),result.stdout)
        A.publish(output/('command-%03d.stderr'%index),result.stderr)
        A.publish(output/('command-%03d.result.json'%index),A.canonical(dict(returncode=result.returncode)))
        A.require(len(result.stdout)<=16*1024**2 and len(result.stderr)<=65536,'command output cap')
        A.require(result.returncode==expected and (not clean_stderr or not result.stderr),
                  'build/neutral command failed: '+str(argv))
        freeze_inputs(dependencies,frozen)
        return result

    plans = []
    for source,original_object,flags in recipes:
        obj = output/original_object.name
        dep = output/(original_object.name+'.d')
        run(['/usr/bin/c++']+flags+['-M','-MT',str(obj),'-MF',str(dep),str(source)],cwd=PRODUCTION/mode)
        before = preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj)
        dependencies.update(before); dependencies.add(original_object)
        freeze_inputs(dependencies,frozen)
        plans.append((source,original_object,flags,obj,dep,before))
    # Every production source/header is frozen before the first compilation.
    members = []
    for source,original_object,flags,obj,dep,before in plans:
        run(['/usr/bin/c++']+flags+['-MD','-MT',str(obj),'-MF',str(dep),
                                  '-o',str(obj),'-c',str(source)],cwd=PRODUCTION/mode)
        A.exact(preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj),before,'production compiler closure')
        original_bytes = A.read_regular(original_object,64*1024**2)
        A.exact(A.read_regular(obj,64*1024**2),original_bytes,'reproduced complete production object')
        A.exact(command(['ar','p',original_archive,original_object.name]),original_bytes,'actual archive member bytes')
        members.append(dict(source=str(source),object=pin(original_object),reproduced=pin(obj),
                            dependencies=[frozen[p] for p in sorted(before)]))
        dependencies.add(obj); freeze_inputs(dependencies,frozen)
    A.exact(command(['ar','t',original_archive]).decode().splitlines(),
            [old.name for _,old,_ in recipes],'nineteen-member archive order')
    reproduced_archive = output/'libwirehair-reproduced.a'
    run(['/usr/bin/ar','qc',reproduced_archive]+[p[3] for p in plans])
    run(['/usr/bin/ranlib',reproduced_archive])
    A.exact(A.read_regular(reproduced_archive,128*1024**2),
            A.read_regular(original_archive,128*1024**2),'byte-identical complete archive')
    proof = dict(protocol=PROTOCOL,mode=mode,archive=pin(original_archive),
                 reproduced=pin(reproduced_archive),members=members,
                 producing_source_closure=True,neutral_inventory=pin(inventory_path),
                 compile_database=pin(database_path),qualification_log=pin(qualification_log))
    A.publish(output/'qualified-library.json',A.canonical(proof))
    archives = [original_archive]
    dependencies.add(output/'qualified-library.json')
    freeze_inputs(dependencies,frozen)
    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_SMALL_COST_NEUTRAL='+str(int(mode!='native')),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += shlex.split(ASAN_FLAGS) if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags.append('-DANDROID=1')
    objects = []
    for source in OBSERVERS:
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
    executable = output/'cost_worker'
    argv = ['/usr/bin/c++', '-fno-lto', '-no-pie']
    if mode == 'asan':
        argv.append('-fsanitize=address,undefined')
    run(argv+list(map(str, objects+archives))+['-pthread', '-Wl,-Map,'+str(output/'link.map'), '-o', str(executable)])
    dependencies.add(executable)
    freeze_inputs(dependencies, frozen)
    loads = [line[5:] for line in A.read_regular(output/'link.map', 4*1024**2).decode().splitlines() if line.startswith('LOAD ')]
    A.require(loads and all(Path(p).is_absolute() for p in loads), 'absolute linker input roster')
    loaded = {Path(p).resolve(strict=True) for p in loads}
    A.require(loaded <= dependencies, 'all linker inputs pinned before link')
    symbols = command(['nm', '-g', '--defined-only', executable]).decode().splitlines()
    names = [line.split()[-1] for line in symbols if line.split()]
    expected = ['GF256Ctx', 'gf256_init_', 'wirehair_encoder_create', 'wirehair_encoder_create_ex',
                'wirehair_encoder_create_owned', 'wirehair_encoder_create_owned_ex',
                'wirehair_decoder_create', 'wirehair_encode', 'wirehair_decode', 'wirehair_recover', 'wirehair_free']
    A.require(not any(name.startswith('wh2_small_') for name in names),'no benchmark boundary linked')
    expected += ['wirehair_v2_'+name for name in ('encoder_create', 'encoder_create_with_options',
                 'encoder_create_profile_id_with_options', 'decoder_create', 'encode', 'decode', 'recover', 'free')]
    A.require(all(names.count(name) == 1 for name in expected), 'one actual WHV2/WH1 API and GF runtime')
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
    A.exact(A.decode(contract), dict(K=8, batch=128, callbacks=77760, cpu_seconds=240, wall_seconds=300,
        work_seconds=180, address_space_mib=384, neutral_only=mode!='native', gf_context_bytes=size,
        claim_path='/var/tmp/wh2-k8-public-cost-r0/CLAIM.json'), 'compiled neutral/resource contract')
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
    text = command(['nm', '-C', executable]).decode().splitlines()
    A.require(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line
                  for line in text)==1, 'one common WORK body')
    fixture_raw = run([executable, '--neutral-fixtures']).stdout
    header = A.decode(fixture_raw)
    reader.verify_fixtures(header, header)
    A.publish(output/'fixtures.json', fixture_raw)
    if mode=='native':
        A.publish(output/'target.json', run([executable, '--neutral-target']).stdout)
        reader.prior_header(output)
    A.exact(run([executable, '--claim-path']).stdout,
            str(reader.OUTPUT/'CLAIM.json').encode()+b'\n', 'compiled claim path')
    neutral = run([executable, '--neutral']).stdout
    A.exact(neutral, b'PASS neutral77760-coordinate roster,216 actual WORK cases, last-call cleanup/capture (216 checked)\n',
            'complete actual WORK qualification')
    A.publish(output/'neutral.txt', neutral)
    rejected = run([executable, '--neutral-profile-rejection']).stdout
    A.exact(rejected, b'PASS installed K8 descriptor rejection\n', 'actual descriptor rejection')
    A.publish(output/'profile-rejection.txt', rejected)
    for case in ('success','last-recover','throw-recover','last-clock','last-source'):
        raw = run([executable, '--neutral-publication', case]).stdout
        reader.verify_publication(raw,case,header)
        A.publish(output/('publication-'+case+'.jsonl'),raw)
    failures = []
    with open('/dev/full','wb') as full:
        result = subprocess.run([str(executable),'--neutral-publication','success'],
            cwd=ROOT, env=process_environment(), stdout=full, stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL, timeout=60)
    A.require(result.returncode==1 and b'output stream' in result.stderr, 'full output device rejects publication')
    failures.append(dict(sink='full',code=result.returncode,stderr=result.stderr.decode()))
    read_fd, write_fd = os.pipe()
    os.close(read_fd)
    try:
        result = subprocess.run([str(executable),'--neutral-publication','success'],
            cwd=ROOT, env=process_environment(), stdout=write_fd, stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL, timeout=60)
    finally:
        os.close(write_fd)
    A.require(result.returncode==1 and b'output stream' in result.stderr, 'broken pipe rejects publication')
    failures.append(dict(sink='broken-pipe',code=result.returncode,stderr=result.stderr.decode()))
    A.publish(output/'publication-output-errors.json',A.canonical(failures))
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
