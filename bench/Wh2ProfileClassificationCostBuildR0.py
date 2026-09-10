#!/usr/bin/env python3
"""Closed, isolated current-library and worker inputs; never launch timing here."""
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import types

HERE = Path(__file__).resolve().parent
def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

P = sibling('_classification_archive_recipes', 'Wh2K8PublicCostBuildR0.py')
A, ROOT = P.A, P.ROOT
PREPARED = Path('/tmp/wh2-profile-classification-inputs.0OlVQ33g')
NEUTRAL = Path('/tmp/wh2-profile-classification-qualified.lu7R3FrP')
BASE = P.PRODUCTION/'native'
DSO = BASE/'libwirehair.so.2.0.0'
DSO_SHA = 'ef684abac606667e30e3b0de1204b0897aa0fa93a4f5ec2268d1ba951dce03b2'
OVERLAY = NEUTRAL/'WirehairV2ProfileClassification.cpp'
OBJECT = NEUTRAL/'CMakeFiles/classification_candidate.dir/WirehairV2ProfileClassification.cpp.o'
OBJECT_SHA = '0abb85d6e4092c5dbbe6a47bc8ae42f44188bc3abb033e631ef5a09e6784f146'
OVERLAY_SHA = '745707ed4e5337965e842a5a17d3541458beb9bfc5d0a5c92d14c4ff4093edb3'
CORE_SHA = '1ae362b8d1296aefccc7f05d29d56a8af1a818b9ae3de06f4b917a223f0cba63'
DEFERRED_SHA = 'fe79e356433f5786ff724955d45e17309c9e24e4f0edfe7ed800c896cda73e34'
WORKER = PREPARED/'Wh2ProfileClassificationWorker.cpp'
COMMON = PREPARED/'Wh2ProfileClassificationCommon.h'
LINK_INPUTS = P.LINK_INPUTS+('crtbeginS.o','crtendS.o','libdl.a')
MUTABLE_HARNESS = ('bench/Wh2ProfileClassificationCostR0.py',
                   'bench/Wh2ProfileClassificationCostBuildR0.py',
                   'bench/test_Wh2ProfileClassificationCostR0.py',
                   'bench/Wh2ProfileClassificationCostR0.md')
pin = P.pin

def imported_files():
    """Include spec-loaded repository modules absent from sys.modules too."""
    pending,seen,result = list(sys.modules.values())+[P],set(),set()
    while pending:
        module = pending.pop()
        if not isinstance(module,types.ModuleType) or id(module) in seen: continue
        seen.add(id(module))
        name = getattr(module,'__file__',None)
        if not name or not Path(name).is_file(): continue
        path = Path(name).resolve(strict=True); result.add(path)
        if ROOT in path.parents:
            pending.extend(value for value in vars(module).values() if isinstance(value,types.ModuleType))
    return result

def verify_link_inputs(path, frozen_inputs):
    loads = [line[5:] for line in A.read_regular(path,4*1024**2).decode().splitlines()
             if line.startswith('LOAD ')]
    A.require(loads and all(Path(p).is_absolute() for p in loads),'absolute linker LOAD inputs')
    resolved = {Path(p).resolve(strict=True) for p in loads}
    A.require(resolved<=set(frozen_inputs),'all actual link inputs frozen before link')
    return resolved

def exact_replace(text, old, new, count=1):
    A.exact(text.count(old),count,'unique generated-source anchor: '+old)
    return text.replace(old,new)

def worker_sources():
    original = A.read_regular(HERE/'Wh2AdmissionRegressionCostR0.cpp',1024**2)
    A.exact(A.sha(original),CORE_SHA,'unchanged historical WORK source')
    common = original.decode()
    common = exact_replace(common,'case_count=32,max_batch=128,callbacks=82944',
                           'case_count=38,max_batch=128,callbacks=98496')
    common = exact_replace(common,'{PublicSmall,5,1280,1},{PublicSmall,5,1280,2}',
        '{PublicSmall,5,1280,1},{PublicSmall,5,1280,2},\n'
        '    {PublicSmall,8,2,1},{PublicSmall,8,2,2},{PublicSmall,8,64,1},{PublicSmall,8,64,2},\n'
        '    {PublicSmall,8,1280,1},{PublicSmall,8,1280,2}')
    common = exact_replace(common,'WIREHAIR_V2_PROFILE_SMALL_K5_2026_09',
        '(c.k==5?WIREHAIR_V2_PROFILE_SMALL_K5_2026_09:WIREHAIR_V2_PROFILE_SMALL_K8_2026_09)',2)
    common = exact_replace(common,'actual ordinary K3/explicit K5 constructor route',
                           'actual ordinary K3/explicit K5/K8 constructor route')
    common = exact_replace(common,'max_batch*(20*1280+128)','max_batch*(22*1280+128)')
    # Only the CPU cap grows prospectively. Wall, work and observer caps stay fixed.
    common = exact_replace(common,'cpu={180,180}','cpu={210,210}')
    common = exact_replace(common,'UINT64_C(180000000000)','UINT64_C(210000000000)',2)
    deferred = A.read_regular(HERE/'Wh2CurrentPreservedDeferredCostR0.cpp',1024**2)
    A.exact(A.sha(deferred),DEFERRED_SHA,'unchanged historical deferred source')
    worker = exact_replace(deferred.decode(),'#include "Wh2AdmissionRegressionCostR0.cpp"',
                           '#include "Wh2ProfileClassificationCommon.h"')
    worker = '#define WH2_ADMISSION_PUBLIC_SMALL 1\n'+worker
    worker = exact_replace(worker,'cpu={180,180}','cpu={210,210}')
    worker = exact_replace(worker,'UINT64_C(180000000000)','UINT64_C(210000000000)',3)
    # Assert literal byte identity of the measured body, not just a source claim.
    start,end = 'NOINLINE void RunWork(', 'void Prepare('
    A.exact(common[common.index(start):common.index(end)],
            original.decode()[original.decode().index(start):original.decode().index(end)],'unchanged WORK body')
    return common.encode(),worker.encode()

def candidate_compile(output, dep):
    return ['/usr/bin/c++','-DNDEBUG','-DWIREHAIR_BUILDING=1',
        '-I'+str(ROOT/'include'),'-I'+str(ROOT/'codec'),'-I'+str(ROOT),
        '-std=gnu++11','-fPIC','-O3','-Wall','-Wextra','-Wpedantic','-Werror',
        '-MD','-MF',str(dep),'-o',str(output),'-c',str(OVERLAY)]

def link(target, objects, map_path=None):
    return ['/usr/bin/c++','-fPIC','-O3','-DNDEBUG','-Wl,--version-script='+str(ROOT/'abi/wirehair.map'),
        '-shared','-Wl,-soname,libwirehair.so.2','-o',str(target)]+list(map(str,objects))+['-lm']+(
        [] if map_path is None else ['-Wl,-Map,'+str(map_path)])

def overlay_inputs():
    source = ROOT/'codec/WirehairV2Profile.cpp'
    raw = A.read_regular(source,1024**2)
    A.exact(A.sha(raw),'975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15','unchanged production facade')
    text = exact_replace(raw.decode(),
        'if (IsSmallProfileId(profile.profile_id)) {\n        const bool shape',
        'if (profile.profile_id != WIREHAIR_V2_PROFILE_CERTIFIED_2026_07) {\n        const bool shape')
    A.exact(A.read_regular(OVERLAY,1024**2),text.encode(),'exact one-predicate overlay')
    A.exact(pin(OVERLAY)['sha256'],OVERLAY_SHA,'qualified overlay')
    A.exact(pin(OBJECT)['sha256'],OBJECT_SHA,'qualified native candidate')
    checks = {'Testing/Temporary/LastTest.log':'bcdf6c1110b26286f5d8b345b637d238d62602a5c1853f8e8825a3df6eb86eba',
        'certified-BASELINE.bin':'2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b',
        'certified-CANDIDATE.bin':'2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b'}
    for name,digest in checks.items(): A.exact(pin(NEUTRAL/name)['sha256'],digest,'retained candidate correctness')
    return {source,OVERLAY,OBJECT}|{p for p in NEUTRAL.rglob('*') if p.is_file()}

def prepare():
    A.require(PREPARED.is_dir() and not PREPARED.is_symlink() and not list(PREPARED.iterdir()),'sole empty preparation directory')
    inventory = P.PRODUCTION/'NEUTRAL_QUALIFIED.json'
    A.exact(pin(inventory)['sha256'],P.INVENTORY_SHA,'retained installed inventory')
    database = BASE/'compile_commands.json'
    A.exact(pin(database)['sha256'],P.DATABASES['native'],'retained compile recipes')
    archive = BASE/'libwirehair.a'
    A.exact(pin(archive)['sha256'],P.ARCHIVES['native'],'actual installed archive')
    A.exact(pin(DSO)['sha256'],DSO_SHA,'actual installed DSO')
    recipes = P.producer_recipes(A.decode(A.read_regular(database,8*1024**2)),'native')
    inputs = overlay_inputs()|{inventory,database,archive,DSO,ROOT/'abi/wirehair.map',Path(__file__).resolve()}
    retained = P.pin_map(A.decode(A.read_regular(inventory,4*1024**2))['files'])
    for path in (ROOT/'CMakeLists.txt',database,archive,BASE/'final-qualification.log'):
        A.exact(pin(path),retained[str(path)],'retained installed qualification input')
        inputs.add(path)
    inputs.update((BASE/'build.ninja',BASE/'CMakeFiles/rules.ninja'))
    inputs.update(p for source,original,_ in recipes for p in (source,original))
    inputs.update(HERE/n for n in ('Wh2ProfileClassificationCostR0.py','Wh2AdmissionRegressionCostR0.cpp',
                                  'Wh2CurrentPreservedDeferredCostR0.cpp','Wh2K8PublicCostBuildR0.py'))
    inputs.add(Path(sys.executable).resolve(strict=True)); inputs.update(imported_files())
    for name in ('c++','cc','as','ld','nm','ar','ranlib','ninja','ldd','bash','git'):
        inputs.add((Path('/usr/bin')/name).resolve(strict=True))
    for name in ('cc1','cc1plus','collect2'):
        inputs.add(Path(P.command(['c++','-print-prog-name='+name]).decode().strip()).resolve(strict=True))
    for name in LINK_INPUTS:
        path = Path(P.command(['c++','-print-file-name='+name]).decode().strip())
        A.require(path.is_absolute(),'resolved linker input'); inputs.add(path.resolve(strict=True))
    for p in list(inputs):
        if os.access(str(p),os.X_OK) and A.read_regular(p,256*1024**2,installed=ROOT not in p.parents and Path('/tmp') not in p.parents)[:4]==b'\x7fELF':
            raw = P.command(['ldd',p]); A.require(b'not found' not in raw,'resolved tool runtimes')
            inputs.update(Path(w).resolve(strict=True) for w in raw.decode().split() if w.startswith('/'))
    frozen,commands = {},[]
    P.freeze_inputs(inputs,frozen)
    def run(argv, cwd=ROOT):
        argv = list(map(str,argv)); P.freeze_inputs(inputs,frozen)
        result = subprocess.run(argv,cwd=cwd,env=P.process_environment(),stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
        index = len(commands); commands.append(dict(argv=argv,cwd=str(cwd),returncode=result.returncode))
        A.publish(PREPARED/('prepare-%03d.stdout'%index),result.stdout)
        A.publish(PREPARED/('prepare-%03d.stderr'%index),result.stderr)
        A.require(len(result.stdout)<=16*1024**2 and len(result.stderr)<=65536,'bounded preparation output')
        A.require(result.returncode==0 and not result.stderr,'preparation command failed: '+str(argv))
        P.freeze_inputs(inputs,frozen); return result.stdout
    plans = []
    for source,original,flags in recipes:
        obj,dep = PREPARED/original.name,PREPARED/(original.name+'.d')
        run(['/usr/bin/c++']+flags+['-M','-MT',str(obj),'-MF',str(dep),str(source)],BASE)
        before = P.preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj)
        inputs.update(before); inputs.add(original); P.freeze_inputs(inputs,frozen)
        plans.append((source,original,flags,obj,dep,before))
    candidate = PREPARED/'proof-profile.o'; candidate_dep = PREPARED/'proof-profile.d'
    args = candidate_compile(candidate,candidate_dep)
    prefix = args[:args.index('-MD')]
    run(prefix+['-M','-MT',str(candidate),'-MF',str(candidate_dep),str(OVERLAY)])
    candidate_before = P.preprocessor_dependencies(A.read_regular(candidate_dep,2*1024**2),candidate)
    inputs.update(candidate_before); P.freeze_inputs(inputs,frozen)
    for source,original,flags,obj,dep,before in plans:
        run(['/usr/bin/c++']+flags+['-MD','-MT',str(obj),'-MF',str(dep),'-o',str(obj),'-c',str(source)],BASE)
        A.exact(P.preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj),before,'complete producer dependencies')
        A.exact(pin(obj)['sha256'],pin(original)['sha256'],'exact original object reproduction')
        A.exact(A.sha(run(['/usr/bin/ar','p',archive,original.name])),pin(original)['sha256'],'actual archive member')
        inputs.add(obj); P.freeze_inputs(inputs,frozen)
    A.exact(run(['/usr/bin/ar','t',archive]).decode().splitlines(),[p[1].name for p in plans],'complete archive order')
    reproduced = PREPARED/'libwirehair-reproduced.a'
    run(['/usr/bin/ar','qc',reproduced]+[p[3] for p in plans]); run(['/usr/bin/ranlib',reproduced])
    A.exact(pin(reproduced)['sha256'],P.ARCHIVES['native'],'whole archive reproduction')
    run(args)
    A.exact(pin(candidate)['sha256'],OBJECT_SHA,'exact frozen candidate reproduction')
    A.exact(P.preprocessor_dependencies(A.read_regular(candidate_dep,2*1024**2),candidate),candidate_before,'candidate compiler closure')
    inputs.add(candidate); P.freeze_inputs(inputs,frozen)
    recorded = run(['/usr/bin/ninja','-C',BASE,'-t','commands',DSO.name]).decode().splitlines()
    A.exact(len(recorded),20,'nineteen compiles and shared link')
    framed = shlex.split(recorded[-1]); A.exact((framed[:2],framed[-2:]),([':','&&'],['&&',':']),'shared link framing')
    original_objects = [p[1] for p in plans]
    expected = link(Path(DSO.name),[p.relative_to(BASE) for p in original_objects])
    A.exact(framed[2:-2],expected,'original ordered shared link')
    baseline_map = PREPARED/'proof-old.map'
    run(link(PREPARED/'proof-old.so',[p[3] for p in plans],baseline_map))
    verify_link_inputs(baseline_map,inputs)
    A.exact(pin(PREPARED/'proof-old.so')['sha256'],DSO_SHA,'complete baseline DSO reproduction')
    replaced = [candidate if p.name=='WirehairV2Profile.cpp.o' else p for p in original_objects]
    A.exact(sum(p==candidate for p in replaced),1,'one replaced object')
    target = PREPARED/'libwirehair.so.2.0.0'; candidate_map = PREPARED/'candidate.map'
    run(link(target,replaced,candidate_map)); verify_link_inputs(candidate_map,inputs)
    common,worker = worker_sources(); A.publish(COMMON,common); A.publish(WORKER,worker)
    inputs.update(imported_files()); P.freeze_inputs(inputs,frozen)
    report = dict(original=pin(DSO),candidate=pin(target),original_objects=list(map(str,original_objects)),
        candidate_objects=list(map(str,replaced)),commands=commands,
        inputs=[frozen[p] for p in sorted(inputs)],artifacts=[pin(p) for p in sorted(PREPARED.iterdir())],
        scientific_launch=False,producing_source_closure=True)
    A.publish(PREPARED/'preparation.json',A.canonical(report))
    print(A.canonical(dict(candidate=pin(target),scientific_launch=False)).decode(),end='')

def prepared_inputs(preparation_sha, harness_sources=()):
    report_path = PREPARED/'preparation.json'
    A.require(type(preparation_sha) is str and len(preparation_sha)==64,'prospectively pinned preparation')
    A.exact(pin(report_path)['sha256'],preparation_sha,'exact independently reviewed preparation')
    report = A.decode(A.read_regular(report_path,4*1024**2))
    A.exact(set(report),{'original','candidate','original_objects','candidate_objects','commands',
                        'inputs','artifacts','scientific_launch','producing_source_closure'},'preparation schema')
    A.exact(report['scientific_launch'],False,'neutral preparation only')
    A.exact(report['producing_source_closure'],True,'closed production inputs')
    # New harness/test/docs may advance before science; never production inputs.
    mutable = {ROOT/name for name in MUTABLE_HARNESS if name in harness_sources}
    P.pin_map(report['inputs']+report['artifacts'])
    expected = [str(BASE/'CMakeFiles/wirehair_objects.dir'/(name+'.o')) for name in P.PRODUCERS]
    A.exact(report['original_objects'],expected,'exact nineteen original objects in link order')
    replacement = [str(PREPARED/'proof-profile.o') if Path(p).name=='WirehairV2Profile.cpp.o' else p for p in expected]
    A.exact(report['candidate_objects'],replacement,'sole fixed Profile replacement')
    A.exact(report['original'],dict(path=str(DSO),bytes=686992,sha256=DSO_SHA),'fixed original DSO')
    A.exact(report['candidate']['path'],str(PREPARED/'libwirehair.so.2.0.0'),'fixed candidate DSO path')
    A.exact(pin(PREPARED/'libwirehair.so.2.0.0'),report['candidate'],'actual prepared candidate DSO')
    inputs = {report_path}
    for record in report['inputs']+report['artifacts']:
        path = Path(record['path']); inputs.add(path)
        if path not in mutable: A.exact(pin(path),record,'unchanged prepared input/artifact')
    for path,raw in zip((COMMON,WORKER),worker_sources()):
        A.exact(A.read_regular(path,1024**2),raw,'exact generated worker')
    inputs.update(overlay_inputs())
    return report,inputs

def provenance(preparation_sha, proof_dir=None, harness_sources=()):
    report,inputs = prepared_inputs(preparation_sha,harness_sources)
    original,candidate = report['original'],report['candidate']
    A.exact(original,dict(path=str(DSO),bytes=686992,sha256=DSO_SHA),'fixed original library')
    reports = [dict(original=p,proof_name=name,proof_sha256=p['sha256'],
                    preparation=pin(PREPARED/'preparation.json'))
               for p,name in ((original,'proof-old.so'),(candidate,'proof-new.so'))]
    reports[1]['proof_object'] = dict(name='proof-profile.o',sha256=OBJECT_SHA)
    if proof_dir is not None:
        for p,name,key in ((original,'proof-old.so','original_objects'),(candidate,'proof-new.so','candidate_objects')):
            target = proof_dir/name; A.require(not target.exists(),'fresh DSO proof')
            map_path = proof_dir/(name+'.map')
            P.command(link(target,[Path(s) for s in report[key]],map_path))
            verify_link_inputs(map_path,inputs)
            A.exact(pin(target)['sha256'],p['sha256'],'exact single-object relink')
        obj,dep = proof_dir/'proof-profile.o',proof_dir/'proof-profile.d'
        P.command(candidate_compile(obj,dep)); A.exact(pin(obj)['sha256'],OBJECT_SHA,'candidate proof recompile')
    return reports,inputs

def build(mode, output, cfg, R):
    """Original observer qualification with pre/post compiler/link input closure."""
    A.require(mode in ('native','asan-driver'),'explicit observer backend')
    A.require(output.is_absolute() and output.name==('native' if mode=='native' else 'asan-driver'),
              'absolute mode-named external build')
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT not in output.parents and output!=ROOT and not output.exists() and
              not output.is_symlink(),'fresh observer directory')
    A.require(all(os.environ.get(k) is None for k in R.ENV_KEYS),'ordinary allocator/loader environment')
    A.exact(P.command(['/usr/bin/c++','-dumpfullversion']).decode().strip(),'13.3.0','qualified compiler')
    expected_san = ('detect_leaks=1:detect_stack_use_after_return=1','halt_on_error=1') if mode=='asan-driver' else (None,None)
    A.exact(tuple(os.environ.get(k) for k in ('ASAN_OPTIONS','UBSAN_OPTIONS')),expected_san,'observer sanitizer environment')
    meta = R.metadata(cfg.libraries); output.mkdir(mode=0o700)
    reports,dependencies = cfg.provenance(output)
    A.publish(output/'library-provenance.json',A.canonical(reports))
    A.publish(output/'AdmissionLibraryBindings.h',R.bindings_header(meta))
    A.publish(output/'library-metadata.json',A.canonical(meta))
    flags = ['-std=c++11','-Wall','-Wextra','-Wpedantic','-Werror','-fno-lto','-fPIC',
        '-DWIREHAIR_STATIC=1','-DWH2_ADMISSION_REGRESSION_NEUTRAL='+str(int(mode!='native')),
        '-I'+str(ROOT),'-I'+str(HERE),'-I'+str(ROOT/'include'),'-I'+str(output)]
    flags += ['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer'] if mode=='asan-driver' else ['-O3','-g1']
    flags += ['-DWH2_ADMISSION_PROTOCOL='+json.dumps(cfg.protocol),
              '-DWH2_ADMISSION_CLAIM_PATH='+json.dumps(str(cfg.output/'CLAIM.json'))]
    sources = [ROOT/cfg.worker_source]+[HERE/n for n in
        ('Wh2FrozenTrace.cpp','Wh2PublicBorrowedTargetIdentity.cpp','Wh2RdpruTargetIdentityV2.cpp')]
    dependencies.update(sources)
    dependencies.update(ROOT/n for n in cfg.sources)
    dependencies.update(output/n for n in ('library-provenance.json','AdmissionLibraryBindings.h','library-metadata.json'))
    dependencies.update(R.NEUTRAL_DIR/('qualified-'+o+'.json') for o in ('old-new','new-old'))
    dependencies.update(imported_files()); dependencies.add(Path(sys.executable).resolve(strict=True))
    # Preparation has already closed compiler/linker/runtime tools. Add the
    # explicit observer startup inputs too, then verify actual LOAD records.
    for name in LINK_INPUTS:
        path = Path(P.command(['c++','-print-file-name='+name]).decode().strip())
        A.require(path.is_absolute(),'resolved observer linker input'); dependencies.add(path.resolve(strict=True))
    frozen,commands = {},[]; P.freeze_inputs(dependencies,frozen)
    def run(argv, expected=0, clean=True):
        argv = list(map(str,argv)); P.freeze_inputs(dependencies,frozen)
        result = subprocess.run(argv,cwd=ROOT,env=P.process_environment(),stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
        index=len(commands); commands.append(dict(argv=argv,returncode=result.returncode))
        A.publish(output/('command-%03d.stdout'%index),result.stdout)
        A.publish(output/('command-%03d.stderr'%index),result.stderr)
        A.require(len(result.stdout)<=16*1024**2 and len(result.stderr)<=65536,'bounded observer command output')
        A.require(result.returncode==expected and (not clean or not result.stderr),'observer command failed: '+str(argv))
        P.freeze_inputs(dependencies,frozen); return result
    plans=[]
    for source in sources:
        obj,dep = output/(source.stem+'.o'),output/(source.stem+'.d')
        run(['/usr/bin/c++']+flags+['-M','-MT',str(obj),'-MF',str(dep),str(source)])
        before=P.preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj)
        dependencies.update(before); P.freeze_inputs(dependencies,frozen)
        plans.append((source,obj,dep,before))
    objects=[]
    for source,obj,dep,before in plans:
        run(['/usr/bin/c++']+flags+['-MD','-MT',str(obj),'-MF',str(dep),'-c',str(source),'-o',str(obj)])
        A.exact(P.preprocessor_dependencies(A.read_regular(dep,2*1024**2),obj),before,'observer compiler input closure')
        objects.append(obj); dependencies.add(obj); P.freeze_inputs(dependencies,frozen)
    executable=output/'cost_worker'
    args=['/usr/bin/c++','-fno-lto','-no-pie']+(['-fsanitize=address,undefined'] if mode=='asan-driver' else [])
    run(args+list(map(str,objects))+['-ldl','-pthread','-Wl,-Map,'+str(output/'link.map'),'-o',str(executable)])
    dependencies.add(executable); P.freeze_inputs(dependencies,frozen)
    verify_link_inputs(output/'link.map',dependencies)
    names=[line.split()[-1] for line in run(['/usr/bin/nm','-g',executable]).stdout.decode().splitlines() if line.split()]
    A.require(not any(n.startswith(('wirehair_','gf256_')) or n=='GF256Ctx' for n in names),'no linked codec runtime')
    R.verify_work_symbol(run(['/usr/bin/nm','-C',executable]).stdout.decode())
    for target in [executable,Path(sys.executable)]+[p for p,_ in cfg.libraries]: dependencies.update(R.runtime_dependencies(target))
    dependencies.update(imported_files()); P.freeze_inputs(dependencies,frozen)
    A.exact(run([executable,'--binding']).stdout,(cfg.protocol+'\n'+str(cfg.output/'CLAIM.json')+'\n').encode(),'compiled claim binding')
    neutral=A.canonical(dict(protocol=cfg.protocol,purpose='neutral claim authentication only'))
    neutral_path=output/'neutral-claim.json'; A.publish(neutral_path,neutral)
    A.exact(run([executable,'--neutral-claim',A.sha(neutral),neutral_path]).stdout,b'','positive claim check')
    rejected=run([executable,'--neutral-claim','0'*64,neutral_path],1,False)
    A.exact((rejected.stdout,rejected.stderr),(b'',b'INVALID: claim bytes\n'),'negative claim check')
    A.publish(output/'claim-binding.json',A.canonical(R.claim_binding(cfg)))
    for order,name in enumerate(('old-new','new-old')):
        A.publish(output/('neutral-'+name+'.txt'),run([executable,'--neutral',name]).stdout)
        raw=run([executable,'--neutral-fixtures',name]).stdout
        cfg.header_checker(A.decode(raw),order,'0'*64,meta,cfg.protocol)
        A.publish(output/('fixtures-'+name+'.json'),raw)
    rejections=[]
    for tail in ([],['--worker'],['--worker','0'*64,'bad-order'],['--worker','0'*64,'old-new'],
                 ['--neutral','bad-order'],['--neutral-fixtures','bad-order'],
                 ['--neutral','old-new','extra'],['--unexpected'],['--neutral-target','old-new']):
        result=run([executable]+tail,1,False)
        A.require(not result.stdout and result.stderr.startswith(b'INVALID:'),'negative observer CLI')
        rejections.append(dict(arguments=tail,returncode=result.returncode,stderr=result.stderr.decode()))
    A.publish(output/'negative-cli.json',A.canonical(rejections))
    P.freeze_inputs(dependencies,frozen)
    cfg.qualify(executable,output,meta,mode)
    P.freeze_inputs(dependencies,frozen)
    dependencies.update(imported_files()); P.freeze_inputs(dependencies,frozen)
    manifest=dict(protocol=cfg.protocol,mode=mode,commands=commands,
        environment={k:os.environ.get(k) for k in R.ENV_KEYS+('ASAN_OPTIONS','UBSAN_OPTIONS')},
        inputs=[frozen[p] for p in sorted(dependencies)],artifacts=[pin(p) for p in sorted(output.iterdir())],
        scientific_launch=False,library_source_provenance_closed=True,sanitized_library_code=False)
    A.publish(output/'manifest.json',A.canonical(manifest))
    print(A.canonical(dict(mode=mode,executable=str(executable),scientific_launch=False)).decode(),end='')
