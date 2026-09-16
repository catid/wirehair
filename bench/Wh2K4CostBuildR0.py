#!/usr/bin/env python3
"""New K4 cost observer build and neutral qualification; never launches science.

The imported recovery module supplies only read-only producing-data parsers
and filesystem helpers. Its build/controller/verifier is never invoked.
"""
import importlib.util
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SPEC = importlib.util.spec_from_file_location('_k4_cost_archive_data', HERE/'Wh2K4RecoveryBuildR0.py')
U = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(U)
A = U.A
MODES, ENV_KEYS, ARCHIVES = U.MODES, U.ENV_KEYS, U.ARCHIVES
SMALL, FIXTURE_SHA, GF_BYTES = U.SMALL, U.FIXTURE_SHA, U.GF_BYTES
ASAN_FLAGS, LINK_INPUTS = U.ASAN_FLAGS, U.LINK_INPUTS
PROTOCOL = 'wirehair.wh2.k4-serialized-cost-r0'
pin, pin_map, command = U.pin, U.pin_map, U.command
process_environment, imported_files = U.process_environment, U.imported_files
freeze_inputs, external_output = U.freeze_inputs, U.external_output
qualified_inputs, boundary_recipe = U.qualified_inputs, U.boundary_recipe
preprocessor_dependencies, context_size = U.preprocessor_dependencies, U.context_size


def verify_qualified_library(provenance):
    """R0 has no reusable qualified producer proof; fail closed for every input.

    A numerical PASS cannot repair missing producer qualification.  This check
    also runs during receipt validation so the old fresh receipt cannot be
    accepted for a new launch, including after deleting its fresh marker.
    A separately reviewed explicit qualifier is required before replacing
    this unconditional block. The historical reader remains unchanged.
    """
    raise ValueError('K4 R0 producing closure is unqualified; preserve spent R0 evidence')


def output_device_check(executable, sink, output, index, commands):
    """Preserve each output-error attempt before checking its result."""
    argv = [str(executable), '--neutral-publication', 'success']
    A.require(sink in ('full','broken-pipe'), 'known output device')
    commands.append(dict(argv=argv, cwd=str(ROOT), stdout_sink=sink))
    prefix = output/('output-device-%d' % index)
    full = None
    write_fd = None
    result = None
    try:
        if sink == 'full':
            full = open('/dev/full','wb')
            target = full
        else:
            read_fd, write_fd = os.pipe()
            os.close(read_fd)
            target = write_fd
        try:
            result = subprocess.run(argv, stdout=target, stderr=subprocess.PIPE,
                                    stdin=subprocess.DEVNULL, env=process_environment(),
                                    cwd=ROOT, timeout=60)
            status = dict(returncode=result.returncode)
            stderr = result.stderr
        except subprocess.TimeoutExpired as error:
            status = dict(timeout=True)
            stderr = error.stderr or b''
        except OSError as error:
            status = dict(error=str(error))
            stderr = b''
    except OSError as error:
        status = dict(error=str(error))
        stderr = b''
    finally:
        try:
            if full is not None: full.close()
            if write_fd is not None: os.close(write_fd)
        except OSError as error:
            status['cleanup_error'] = str(error)
    A.publish(Path(str(prefix)+'.stderr'),stderr)
    A.publish(Path(str(prefix)+'.json'),A.canonical(dict(argv=argv,cwd=str(ROOT),sink=sink,**status)))
    A.require(result is not None and result.returncode==1 and 'cleanup_error' not in status and
              len(stderr)<=65536 and b'output stream' in stderr,
              'output device rejects publication: '+sink)
    return dict(sink=sink,code=result.returncode,stderr=stderr.decode())


def build(mode, output):
    """Fresh observer build and bounded neutral checks, never scientific work."""
    A.require(mode in MODES, 'backend')
    output = Path(output)
    A.require(output.name == mode, 'mode-named build directory')
    output = output.parent.resolve(strict=True)/output.name
    external_output(output)
    A.require(ROOT not in output.parents and output != ROOT and not output.exists() and not output.is_symlink(), 'fresh external build')
    dependencies = {Path(__file__).resolve(), HERE/'Wh2K4SerializedCostR0.py'}
    frozen = {}
    freeze_inputs(dependencies, frozen)
    spec = importlib.util.spec_from_file_location('_k4_cost_reader_for_build', HERE/'Wh2K4SerializedCostR0.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    A.exact((reader.PROTOCOL, reader.MODES), (PROTOCOL, MODES), 'new cost reader contract')
    A.exact({k: os.environ.get(k) for k in reader.ENV_KEYS},
            dict({k: None for k in reader.ENV_KEYS}, **reader.SANITIZERS), 'ordinary allocator/sanitizer environment')
    dependencies.update(ROOT/name for name in reader.NEW)
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
    verify_qualified_library(provenance)
    A.publish(output/'qualified-library.json', A.canonical(provenance))
    dependencies.add(output/'qualified-library.json'); freeze_inputs(dependencies, frozen)

    flags = ['-std=c++11', '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-fno-lto', '-fPIC',
             '-DWIREHAIR_STATIC=1', '-DWH2_SMALL_CODEC_K=4', '-DWH2_SMALL_COST_NEUTRAL='+str(int(mode != 'native')),
             '-I'+str(ROOT), '-I'+str(HERE), '-I'+str(ROOT/'include'), '-I'+str(output)]
    flags += shlex.split(ASAN_FLAGS) if mode == 'asan' else ['-O3', '-g1']
    if mode == 'scalar':
        flags.append('-DANDROID=1')
    objects = []
    for source in (HERE/'Wh2K4SerializedCostR0.cpp', HERE/'Wh2FrozenTrace.cpp',
                   HERE/'Wh2PublicBorrowedTargetIdentity.cpp', HERE/'Wh2RdpruTargetIdentityV2.cpp'):
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
    contract = run([executable,'--contract']).stdout
    A.exact(A.decode(contract),dict(K=reader.K,batch=reader.BATCH,callbacks=reader.CALLBACKS,cpu_seconds=240,
            wall_seconds=300,work_seconds=180,address_space_mib=384,
            claim_path=str(reader.OUTPUT/'CLAIM.json'),neutral_only=mode!='native',
            gf_context_bytes=size),'compiled resource/launch/GF ABI contract')
    A.publish(output/'contract.json',contract)
    text = command(['/usr/bin/nm', '-C', executable]).decode().splitlines()
    A.require(sum('RunWork(' in line and '.cold' not in line and '[clone' not in line for line in text) == 1, 'one common WORK body')
    fixture_raw = run([executable, '--neutral-fixtures']).stdout
    header = A.decode(fixture_raw)
    reader.verify_fixtures(header, header)
    A.publish(output/'fixtures.json', fixture_raw)
    if mode == 'native':
        A.publish(output/'target.json', run([executable, '--neutral-target']).stdout)
        reader.prior_header(output)
    # The very same Authenticate helper gates Worker and this positive test.
    # Do not create a scientific namespace merely to test launch binding.
    A.exact(run([executable,'--claim-path']).stdout,str(reader.OUTPUT/'CLAIM.json').encode()+b'\n','compiled claim path')
    neutral_claim = output/'neutral-claim.json'
    neutral_bytes = A.canonical(dict(protocol=PROTOCOL,neutral=True))
    A.publish(neutral_claim,neutral_bytes)
    A.exact(run([executable,'--neutral-claim',neutral_claim,A.sha(neutral_bytes)]).stdout,
            b'PASS claim authentication\n','positive claim authentication')
    negatives = []
    for argv in ([executable,'--neutral-claim',neutral_claim,'0'*64],
                 [executable,'--neutral-claim',output/'absent',A.sha(neutral_bytes)],
                 [executable,'--neutral-claim',neutral_claim,'g'*64],
                 [executable],[executable,'--unknown'],[executable,'--worker','0']):
        p = run(argv,expected=1,clean_stderr=False)
        A.require(p.returncode==1 and not p.stdout and p.stderr.startswith(b'INVALID:'),'negative CLI/authentication')
        negatives.append(dict(argv=list(map(str,argv)),code=p.returncode,stderr=p.stderr.decode()))
    A.publish(output/'negative-cli.json',A.canonical(negatives))
    neutral=run([executable,'--neutral']).stdout
    A.exact(neutral,b'PASS neutral77760-coordinate roster,216 actual WORK cases, last-call cleanup/capture (216 checked)\n','complete neutral WORK qualification')
    A.publish(output/'neutral.txt',neutral)
    for case in ('success','last-recover','throw-recover','last-clock','last-source'):
        raw = run([executable,'--neutral-publication',case]).stdout
        reader.verify_publication(raw,case,header)
        A.publish(output/('publication-'+case+'.jsonl'),raw)
    # Real stdio failures, not only a mock sink. Never launch scientific WORK.
    failures = []
    for index, sink in enumerate(('full','broken-pipe')):
        freeze_inputs(dependencies,frozen)
        failures.append(output_device_check(executable,sink,output,index,commands))
        freeze_inputs(dependencies,frozen)
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
