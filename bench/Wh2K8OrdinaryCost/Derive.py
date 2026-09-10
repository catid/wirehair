#!/usr/bin/env python3
"""Prepare the ordinary-K8 lifecycle instrument; never launch scientific timing."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / 'bench'
PROTOCOL = 'wirehair.wh2.k8-ordinary-cost-r0'
TEMPLATES = {
    'Wh2K8PublicCostR0.cpp': ('4bd7b3904f1caf0f49bac390139a132959e32605a117ffb8dbee586eff21cbdc', (0, 1, 1)),
    'Wh2K8PublicCostR0.py': ('5b1961ebcf5fca25804844c4dbe5b537c8dc1988dcb8b2a5ef7434c7be1c216a', (7, 1, 1)),
    'Wh2K8PublicCostBuildR0.py': ('2672d0d9e1ef2aa20c26cd839b9a49a560ec22e70332f41c8e2998d5cc1b58b6', (2, 1, 1)),
    'test_Wh2K8PublicCostR0.py': ('ea3ec4cb70b0e34e6d0c400f376c0834a45b53632afcf79a07b3a49e4cf162be', (2, 0, 1)),
    'test_Wh2K8PublicCostBuildR0.py': ('2721458c24f039acdebd9917aa35bbd98cc49f72c88ca2729085c78fd4a438ea', (2, 0, 0)),
}
RENAMES = (
    ('Wh2K8PublicCost', 'Wh2K8OrdinaryCost'),
    ('wirehair.wh2.k8-public-cost-r0', PROTOCOL),
    ('/var/tmp/wh2-k8-public-cost-r0', '/var/tmp/wh2-k8-ordinary-cost-r0'),
)
SOURCES = tuple('bench/' + p for p in TEMPLATES) + (
    'bench/Wh2K8PublicCostR0.md', 'bench/Wh2AlignedIntermediateCostR0.py',
    'bench/Wh2K8OrdinaryCost/Derive.py', 'bench/Wh2K8OrdinaryCost/Candidate.py',
    'bench/Wh2K8OrdinaryCost/test_Derive.py', 'bench/Wh2K8OrdinaryCost/README.md',
    'bench/Wh2K8OrdinarySelector/CMakeLists.txt', 'bench/Wh2K8OrdinarySelector/Overlay.cmake',
    'bench/Wh2K8OrdinarySelector/TestOverlay.cmake', 'bench/Wh2K8OrdinarySelector/SelectorChecks.inc',
    'bench/Wh2K8OrdinarySelector/ContractTest.cpp', 'bench/Wh2K8OrdinarySelector/README.md',
)


def require(value, why):
    if not value:
        raise ValueError(why)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def replace(text, old, new, count=1):
    require(old and text.count(old) == count, 'missing/ambiguous derivation anchor: ' + old[:80])
    return text.replace(old, new)


def templates():
    result = {}
    for name, (digest, _) in TEMPLATES.items():
        raw = (HERE / name).read_bytes()
        require(sha(raw) == digest, 'historical template changed: ' + name)
        result[name] = raw.decode('utf-8')
    return result


def derive(inputs):
    """Pure text derivation; historical controllers are never imported/executed."""
    require(set(inputs) == set(TEMPLATES), 'complete template roster')
    output = {}
    for name, text in inputs.items():
        require(sha(text.encode()) == TEMPLATES[name][0], 'unauthenticated template')
        for (old, new), count in zip(RENAMES, TEMPLATES[name][1]):
            text = replace(text, old, new, count)
        output[name.replace('PublicCost', 'OrdinaryCost')] = text

    name = 'Wh2K8OrdinaryCostR0.cpp'
    text = output[name]
    old = """        wirehair_v2_encoder_create(source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_with_options(source,message,b,&o,p.data(),32,&n,&handle);"""
    text = replace(text, old, """        wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,message,b,&o,p.data(),32,&n,&handle);""")
    old = """    const int r=wirehair_v2_encoder_create_profile_id_with_options(
        WIREHAIR_V2_PROFILE_SMALL_K8_2026_09,source,message,b,&o,p.data(),32,&n,&handle);"""
    text = replace(text, old, """    const int r=policy==WirehairV2EncoderSource_Independent ?
        wirehair_v2_encoder_create(source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_with_options(source,message,b,&o,p.data(),32,&n,&handle);""")
    text = replace(text, 'ordinary control still selects certified equations',
                   'explicit certified control selects certified equations')
    text = replace(text, '// Same installed candidate constructor routes as the retained recovery gate.',
                   '// Actual ordinary candidate routes; no explicit-profile timing inheritance.')
    text = replace(text, '// Actual installed WHV2 K8 lifecycle gate. Spent prototype sources stay unchanged.',
                   '// Ordinary K8 selector lifecycle gate. All historical instruments stay unchanged.')
    output[name] = text

    name = 'Wh2K8OrdinaryCostR0.py'
    text = output[name]
    text = replace(text, 'Path(__file__).with_name(filename)',
                   "(Path(filename) if Path(filename).is_absolute() else Path(__file__).with_name(filename))")
    text = replace(text, "'Wh2AlignedIntermediateCostR0.py')",
                   repr(str(HERE / 'Wh2AlignedIntermediateCostR0.py')) + ')')
    start = text.index('SOURCES = (')
    end = text.index('\n\n', start)
    text = replace(text, text[start:end], 'SOURCES = ' + repr(SOURCES))
    text = replace(text, "('ordinary_independent','k8_independent','wh1_owned',\n"
                        "             'ordinary_borrowed','k8_borrowed','wh1_borrowed')",
                        "('certified_independent','ordinary_k8_independent','wh1_owned',\n"
                        "             'certified_borrowed','ordinary_k8_borrowed','wh1_borrowed')")
    output[name] = text

    name = 'Wh2K8OrdinaryCostBuildR0.py'
    text = output[name]
    text = replace(text, 'HERE = Path(__file__).resolve().parent\nROOT = HERE.parent',
                   'GENERATED = Path(__file__).resolve().parent\nROOT = Path(' + repr(str(ROOT)) +
                   ')\nHERE = ROOT / "bench"')
    anchor = 'SPEC.loader.exec_module(A)'
    text = replace(text, anchor, anchor + '\n'
        'SPEC = importlib.util.spec_from_file_location("_ordinary_k8_candidate", HERE/"Wh2K8OrdinaryCost/Candidate.py")\n'
        'C = importlib.util.module_from_spec(SPEC)\nSPEC.loader.exec_module(C)')
    old = """OBSERVERS = tuple(HERE/name for name in ('Wh2K8OrdinaryCostR0.cpp', 'Wh2FrozenTrace.cpp',
    'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp'))"""
    text = replace(text, old, """OBSERVERS = (GENERATED/'Wh2K8OrdinaryCostR0.cpp',) + tuple(HERE/name for name in (
    'Wh2FrozenTrace.cpp', 'Wh2PublicBorrowedTargetIdentity.cpp', 'Wh2RdpruTargetIdentityV2.cpp'))""")
    text = replace(text, "HERE/'Wh2K8OrdinaryCostR0.py'", "GENERATED/'Wh2K8OrdinaryCostR0.py'")
    anchor = "dependencies = {ROOT/name for name in reader.SOURCES}|set(OBSERVERS)|imported_files()"
    text = replace(text, anchor, anchor + '\n'
                   '    dependencies.add(GENERATED/"DERIVATION.json")\n'
                   '    dependencies.update(C.verify_derivation(GENERATED))')
    anchor = '    # Every production source/header is frozen before the first compilation.'
    text = replace(text, anchor, """    candidate = C.plan(mode, output, dependencies, frozen, run, freeze_inputs,
                       preprocessor_dependencies, pin)
""" + anchor)
    anchor = '    proof = dict(protocol=PROTOCOL,mode=mode,archive=pin(original_archive),'
    text = replace(text, anchor, """    candidate_proof = C.compile_candidate(candidate, dependencies, frozen, run,
                                          freeze_inputs, preprocessor_dependencies, pin)
""" + anchor)
    text = replace(text, 'compile_database=pin(database_path),qualification_log=pin(qualification_log))',
                   'compile_database=pin(database_path),qualification_log=pin(qualification_log),\n'
                   '                 ordinary_selector=candidate_proof)')
    text = replace(text, '    archives = [original_archive]',
                   '    archives = [candidate["object"], original_archive]')
    anchor = "    A.require(loaded <= dependencies, 'all linker inputs pinned before link')"
    text = replace(text, anchor, anchor + '\n'
                   '    C.check_link(A.read_regular(output/"link.map", 4*1024**2).decode(),\n'
                   '                 original_archive, candidate["object"], [p[1].name for p in recipes])')
    output[name] = text

    name = 'test_Wh2K8OrdinaryCostR0.py'
    output[name] = replace(output[name], "(C.ROOT/'bench/Wh2K8OrdinaryCostR0.cpp')",
                           "(Path(__file__).with_name('Wh2K8OrdinaryCostR0.cpp'))")
    output[name] = replace(output[name],
        "('ordinary_independent','k8_independent','wh1_owned',\n"
        "                                     'ordinary_borrowed','k8_borrowed','wh1_borrowed')",
        "('certified_independent','ordinary_k8_independent','wh1_owned',\n"
        "                                     'certified_borrowed','ordinary_k8_borrowed','wh1_borrowed')")
    for name, text in output.items():
        # Old source paths in the controller SOURCES tuple are deliberate pins.
        if name != 'Wh2K8OrdinaryCostR0.py':
            require('/var/tmp/wh2-k8-public-cost-r0' not in text, 'stale scientific path')
        if name.endswith('.py'):
            compile(text, name, 'exec')
    return output


def encoded_outputs():
    return {name: text.encode() for name, text in derive(templates()).items()}


def record(outputs):
    return dict(protocol=PROTOCOL, root=str(ROOT), generator_sha=sha(Path(__file__).read_bytes()),
                templates={name: digest for name, (digest, _) in TEMPLATES.items()},
                outputs={name: sha(raw) for name, raw in outputs.items()})


def verify(directory):
    directory = Path(directory).resolve(strict=True)
    outputs = encoded_outputs()
    actual = json.loads((directory/'DERIVATION.json').read_bytes())
    require(actual == record(outputs), 'derivation receipt changed')
    for name, raw in outputs.items():
        path = directory/name
        require(not path.is_symlink() and path.is_file() and path.read_bytes() == raw,
                'generated source changed: ' + name)
    return directory


def prepare(directory):
    directory = Path(directory)
    require(directory.is_absolute() and directory.name not in ('', '.', '..'), 'absolute fresh directory')
    parent = directory.parent.resolve(strict=True)
    directory = parent/directory.name
    require(directory != ROOT and ROOT not in directory.parents and not directory.exists() and
            not directory.is_symlink(), 'fresh external derivation only')
    outputs = encoded_outputs()
    directory.mkdir(mode=0o700)
    for name, raw in outputs.items():
        with (directory/name).open('xb') as stream:
            stream.write(raw)
    with (directory/'DERIVATION.json').open('xb') as stream:
        stream.write((json.dumps(record(outputs), sort_keys=True, separators=(',', ':'))+'\n').encode())
    verify(directory)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'verify', 'build'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--mode', choices=('native', 'scalar', 'asan'))
    args = parser.parse_args()
    require((args.command == 'build') == (args.mode is not None), 'mode required only for build')
    directory = prepare(args.directory) if args.command == 'prepare' else verify(args.directory)
    if args.command == 'build':
        # No run/receipt/replay forwarding before independent qualification.
        return subprocess.call([sys.executable, str(directory/'Wh2K8OrdinaryCostR0.py'),
                                'build', args.mode, str(directory/args.mode)])
    print(json.dumps(dict(protocol=PROTOCOL, directory=str(directory), timing_launched=False)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
