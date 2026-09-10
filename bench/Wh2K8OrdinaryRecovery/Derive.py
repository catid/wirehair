#!/usr/bin/env python3
"""Prepare an ordinary-K8 retained recovery instrument; never launch science."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT/'bench'
PROTOCOL = 'wirehair.wh2.k8-ordinary-recovery-r0'
TEMPLATES = {
    'Wh2K8PublicRecoveryR0.cpp': ('61a4fd6e431c77ee64a7036d3e5fe2f9cb1dda82bfb7ee7c15597e6e894a9a66', (0, 1, 1)),
    'Wh2K8PublicRecoveryR0.py': ('143247c2ad286d1a007e6e4d69c7afcceee572a8e83d67807321e63ff18f5014', (7, 1, 1)),
    'Wh2K8PublicRecoveryBuildR0.py': ('604ddcb1eba0032bacfd091d25e12194b3a8342414afc8ba8a8cfe3922926fb0', (2, 1, 1)),
    'test_Wh2K8PublicRecoveryR0.py': ('737b880ef9ae353c2f69d536f20bf9eccdf28278c8a60b774351406763acb5aa', (1, 0, 0)),
    'test_Wh2K8PublicRecoveryBuildR0.py': ('03c4d68333e35d1bd60de94c136a31afcec331e81bd479f28db73bd126b4f465', (1, 0, 0)),
}
RENAMES = (
    ('Wh2K8PublicRecovery', 'Wh2K8OrdinaryRecovery'),
    ('wirehair.wh2.k8-public-recovery-r0', PROTOCOL),
    ('/var/tmp/wh2-k8-public-recovery-r0', '/var/tmp/wh2-k8-ordinary-recovery-r0'),
)
SOURCES = tuple('bench/'+name for name in TEMPLATES) + (
    'bench/Wh2K8PublicRecoveryR0.md', 'bench/Wh2K8NativeDataR0.py',
    'bench/Wh2K3NativeDataR0.py', 'bench/Wh2NoncommutingRadixRunR0.py',
    'bench/Wh2AlignedIntermediateCostR0.py', 'V2_WIRE_PROFILE.md',
    'WH2_BORROWED_SOURCE_API.md', 'bench/Wh2K8OrdinaryCost/Candidate.py',
    'bench/Wh2K8OrdinaryRecovery/Derive.py', 'bench/Wh2K8OrdinaryRecovery/Retained.py',
    'bench/Wh2K8OrdinaryRecovery/test_Derive.py', 'bench/Wh2K8OrdinaryRecovery/README.md',
    'bench/Wh2K8OrdinarySelector/CMakeLists.txt', 'bench/Wh2K8OrdinarySelector/Overlay.cmake',
    'bench/Wh2K8OrdinarySelector/TestOverlay.cmake', 'bench/Wh2K8OrdinarySelector/SelectorChecks.inc',
    'bench/Wh2K8OrdinarySelector/ContractTest.cpp', 'bench/Wh2K8OrdinarySelector/README.md',
)


def require(value, why):
    if not value:
        raise ValueError(why)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def replace(text, old, new, count=1):
    require(old and text.count(old) == count, 'missing/ambiguous derivation anchor: '+old[:80])
    return text.replace(old, new)


def templates():
    result = {}
    for name, (digest, _) in TEMPLATES.items():
        raw = (HERE/name).read_bytes()
        require(sha(raw) == digest, 'historical template changed: '+name)
        result[name] = raw.decode('utf-8')
    return result


def derive(inputs):
    """Old scientific instruments are authenticated inert text, never executed."""
    require(set(inputs) == set(TEMPLATES), 'complete template roster')
    output = {}
    for name, text in inputs.items():
        require(sha(text.encode()) == TEMPLATES[name][0], 'unauthenticated template')
        for (old, new), count in zip(RENAMES, TEMPLATES[name][1]):
            text = replace(text, old, new, count)
        output[name.replace('PublicRecovery', 'OrdinaryRecovery')] = text

    name = 'Wh2K8OrdinaryRecoveryR0.cpp'
    text = output[name]
    text = replace(text, '''        wirehair_v2_encoder_create(source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_with_options(source,message,b,&o,p.data(),32,&n,&handle);''',
        '''        wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,source,message,b,&o,p.data(),32,&n,&handle);''')
    text = replace(text, '''    const int r=wirehair_v2_encoder_create_profile_id_with_options(
        WIREHAIR_V2_PROFILE_SMALL_K8_2026_09,source,message,b,&o,p.data(),32,&n,&handle);''',
        '''    const int r=policy==WirehairV2EncoderSource_Independent ?
        wirehair_v2_encoder_create(source,message,b,p.data(),32,&n,&handle) :
        wirehair_v2_encoder_create_with_options(source,message,b,&o,p.data(),32,&n,&handle);''')
    text = replace(text, 'ordinary selects certified profile', 'explicit control selects certified profile')
    text = replace(text, '// Actual WHV2 candidate boundary. No benchmark facade or descriptor translation.',
                   '// Actual ordinary WHV2 candidate constructors; no prototype facade or descriptor translation.')
    output[name] = text

    name = 'Wh2K8OrdinaryRecoveryR0.py'
    text = output[name]
    text = replace(text, "'Wh2K8NativeDataR0.py')", repr(str(HERE/'Wh2K8NativeDataR0.py'))+')')
    anchor = "D = sibling('_k8_recovery_retained_data', "
    end = text.index('\n', text.index(anchor))
    text = text[:end] + "\nI = sibling('_ordinary_k8_installed_parity', "+repr(str(HERE/'Wh2K8OrdinaryRecovery/Retained.py'))+")" + text[end:]
    text = replace(text, "('ordinary_wh2_independent','k8_independent','wh1_owned',\n"
                         "        'ordinary_wh2_borrowed','k8_borrowed','wh1_borrowed')",
                         "('certified_independent','ordinary_k8_independent','wh1_owned',\n"
                         "        'certified_borrowed','ordinary_k8_borrowed','wh1_borrowed')")
    start = text.index('SOURCES = ('); end = text.index('\n\n', start)
    text = replace(text, text[start:end], 'SOURCES = '+repr(SOURCES))
    text = replace(text, 'ordinary certified WH2 descriptor', 'explicit certified WH2 descriptor')
    anchor = "    A.exact(len(declared), len(receipt['pins']), 'unique receipt pins')"
    text = replace(text, anchor, anchor+"\n    A.require({str(p) for p in I.PINS} <= set(declared), 'installed evidence closure')")
    anchor = "        _,records = verify(A.read_regular(Path(receipt['executables'][mode]).parent/'fixtures.jsonl',RAW_CAP),'0'*64,mode,True)"
    text = replace(text, anchor, anchor+'\n        I.neutral(records, mode, A)')
    text = replace(text, '        retained_parity(reference)',
                   '        retained_parity(reference)\n        I.verify(reference, A, U.pin)', 2)
    output[name] = text

    name = 'Wh2K8OrdinaryRecoveryBuildR0.py'
    text = output[name]
    text = replace(text, 'HERE = Path(__file__).resolve().parent\nROOT = HERE.parent',
                   'GENERATED = Path(__file__).resolve().parent\nROOT = Path('+repr(str(ROOT))+')\nHERE = ROOT/"bench"')
    anchor = 'SPEC.loader.exec_module(A)'
    text = replace(text, anchor, anchor+'\n'
        'SPEC = importlib.util.spec_from_file_location("_ordinary_k8_candidate", HERE/"Wh2K8OrdinaryCost/Candidate.py")\n'
        'C = importlib.util.module_from_spec(SPEC)\nSPEC.loader.exec_module(C)\n'
        'SPEC = importlib.util.spec_from_file_location("_ordinary_k8_recovery_derivation", HERE/"Wh2K8OrdinaryRecovery/Derive.py")\n'
        'G = importlib.util.module_from_spec(SPEC)\nSPEC.loader.exec_module(G)')
    text = replace(text, "HERE/'Wh2K8OrdinaryRecoveryR0.py'", "GENERATED/'Wh2K8OrdinaryRecoveryR0.py'")
    anchor = '    dependencies = {ROOT/name for name in reader.SOURCES}|imported_files()'
    text = replace(text, anchor, anchor+'\n'
        '    G.verify(GENERATED)\n'
        '    dependencies.update(GENERATED/name for name in G.encoded_outputs())\n'
        '    dependencies.add(GENERATED/"DERIVATION.json")\n'
        '    dependencies.update(reader.I.PINS)')
    anchor = '    # Every production source/header is frozen before the first compilation.'
    text = replace(text, anchor, '    candidate = C.plan(mode, output, dependencies, frozen, run, freeze_inputs,\n'
                   '                       preprocessor_dependencies, pin)\n'+anchor)
    anchor = '    proof = dict(protocol=PROTOCOL,mode=mode,archive=pin(original_archive),'
    text = replace(text, anchor, '    candidate_proof = C.compile_candidate(candidate, dependencies, frozen, run,\n'
                   '                                          freeze_inputs, preprocessor_dependencies, pin)\n'+anchor)
    text = replace(text, 'compile_database=pin(database_path),qualification_log=pin(qualification_log))',
                   'compile_database=pin(database_path),qualification_log=pin(qualification_log),\n'
                   '                 ordinary_selector=candidate_proof)')
    text = replace(text, '    archives = [original_archive]', '    archives = [candidate["object"], original_archive]')
    text = replace(text, "HERE/'Wh2K8OrdinaryRecoveryR0.cpp'", "GENERATED/'Wh2K8OrdinaryRecoveryR0.cpp'")
    anchor = "    A.require(loaded <= dependencies, 'all linker inputs pinned before link')"
    text = replace(text, anchor, anchor+'\n'
                   '    C.check_link(A.read_regular(output/"link.map", 4*1024**2).decode(),\n'
                   '                 original_archive, candidate["object"], [p[1].name for p in recipes])')
    text = replace(text, "'encoder_create_profile_id_with_options', 'decoder_create'",
                   "'encoder_create_profile_id', 'encoder_create_profile_id_with_options', 'decoder_create'")
    text = replace(text, "    reader.verify(fixture_raw, '0'*64, mode, neutral=True)",
                   "    _, fixture_records = reader.verify(fixture_raw, '0'*64, mode, neutral=True)\n"
                   '    reader.I.neutral(fixture_records, mode, A)')
    output[name] = text

    name = 'test_Wh2K8OrdinaryRecoveryR0.py'
    text = replace(output[name], "'ordinary certified'", "'explicit certified'")
    anchor = "        source=add('/synthetic/source.cpp',b'source')"
    text = replace(text, anchor, anchor+'\n        installed=[add(str(p), b"installed") for p in C.I.PINS]')
    text = replace(text, 'inputs=[source],artifacts=[artifact,fixture]', 'inputs=[source]+installed,artifacts=[artifact,fixture]')
    text = replace(text, "patch.object(C,'verify',return_value=({},[])):",
                   "patch.object(C,'verify',return_value=({},[])),patch.object(C.I,'neutral'):")
    text = replace(text, "lambda r:r['pins'].remove(source),lambda r:r['pins'].append(source),",
                   "lambda r:r['pins'].remove(source),lambda r:r['pins'].append(source),\n"
                   "                           lambda r:r['pins'].remove(installed[0]),")
    output[name] = text
    for name, text in output.items():
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
    require(json.loads((directory/'DERIVATION.json').read_bytes()) == record(outputs), 'derivation receipt changed')
    for name, raw in outputs.items():
        path = directory/name
        require(not path.is_symlink() and path.is_file() and path.read_bytes() == raw,
                'generated source changed: '+name)
    return directory


def prepare(directory):
    directory = Path(directory)
    require(directory.is_absolute() and directory.name not in ('', '.', '..'), 'absolute fresh directory')
    directory = directory.parent.resolve(strict=True)/directory.name
    require(directory != ROOT and ROOT not in directory.parents and not directory.exists() and
            not directory.is_symlink(), 'fresh external derivation only')
    outputs = encoded_outputs()
    directory.mkdir(mode=0o700)
    for name, raw in outputs.items():
        with (directory/name).open('xb') as stream:
            stream.write(raw)
    with (directory/'DERIVATION.json').open('xb') as stream:
        stream.write((json.dumps(record(outputs), sort_keys=True, separators=(',', ':'))+'\n').encode())
    return verify(directory)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'verify', 'build'))
    parser.add_argument('directory', type=Path)
    parser.add_argument('--mode', choices=('native', 'scalar', 'asan'))
    args = parser.parse_args()
    require((args.command == 'build') == (args.mode is not None), 'mode required only for build')
    directory = prepare(args.directory) if args.command == 'prepare' else verify(args.directory)
    if args.command == 'build':
        return subprocess.call([sys.executable, str(directory/'Wh2K8OrdinaryRecoveryR0.py'),
                                'build', args.mode, str(directory/args.mode)])
    print(json.dumps(dict(protocol=PROTOCOL, directory=str(directory), science_launched=False)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
