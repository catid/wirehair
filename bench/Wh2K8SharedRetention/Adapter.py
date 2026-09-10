#!/usr/bin/env python3
"""Prepare/qualify gate A's actual shared observer. No scientific launch command."""
import argparse
import importlib.util
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


P = module('_k8_retention_shared_proof', ROOT/'bench/Wh2K8OrdinaryShared/Prepare.py')
F = module('_k8_retention_fixture_helpers', ROOT/'bench/Wh2ProfileClassificationCostR0.py')
S, R, D, A, B = F.B, F.R, F.D, P.A, P.B
PROTOCOL = 'wirehair.wh2.k8-ordinary-shared-retention-r0'
OUTPUT = Path('/var/tmp/wh2-k8-ordinary-shared-retention-r0')
PREPARED = Path('/tmp/wh2-k8-ordinary-shared-final.Iy6tBm1H/prepared')
PREPARATION_SHA = 'd98b53db95464544d73b4f1cc5f68eaf5b00d8c6062b23e5a91112e7a55b7379'
QUALIFICATIONS = (
    (PREPARED.parent/'neutral-closed/QUALIFIED.json',
     'db047c9c3d7eb1f526cb99596630a255bc8717b042498317af0784c24decce82'),
    (Path('/tmp/wh2-k8-ordinary-package.8vNek285/neutral-helpers/PACKAGE_QUALIFIED.json'),
     '625e5236976fea2d315320a67f926572b8f136badc909cd9d4644200b68eaeaf'))
LIBRARIES = ((PREPARED/'baseline/libwirehair.so.2.0.0', P.DSO_SHA),
             (PREPARED/'candidate/libwirehair.so.2.0.0', P.CANDIDATE_DSO_SHA))
CASES = R.CASES + tuple((4,k,b,p) for k in (3,5,8) for b in (2,64,1280) for p in (1,2))
COMMON = 'Wh2K8RetentionCommon.h'
WORKER = 'Wh2K8RetentionWorker.cpp'
HELPERS = {
    'Wh2ProfileClassificationCostBuildR0.py': 'd8ac4d6094ef4cd06ecee6f7cb260882928101f82af1bc8657e94005e3fbd54f',
    'Wh2ProfileClassificationCostR0.py': '57e0d699c604fd888c93041e76dfb660736e9123eedaee093171951b73c67567',
    'Wh2AdmissionRegressionCostR0.py': '167e08af526f7527e6afa12276c465bbb9d844b126efd6bbdbfb78e923921dcb',
    'Wh2CurrentPreservedDeferredCostR0.py': 'ae1f2207c541518a4a2cb6a652120ae51de6e80565970f4a15dff9d1264d2f7f'}
SLOTS = ('wirehair_recover_block_ex', 'wirehair_v2_encoder_create',
         'wirehair_v2_encoder_create_profile', 'wirehair_v2_encoder_create_profile_id',
         'wirehair_v2_profile_deserialize', 'wirehair_v2_profile_serialize')


def fixed_helpers():
    for name, digest in HELPERS.items():
        A.exact(B.pin(ROOT/'bench'/name)['sha256'], digest, 'unchanged reusable observer helper')


def worker_sources():
    fixed_helpers()
    # This helper generates ONLY the shared observer, not its rejected codec.
    # Keep the existing 38-case adapter, six/six bindings, WORK and deferred IO.
    common, worker = S.worker_sources()
    worker = S.exact_replace(worker.decode(), '#include "Wh2ProfileClassificationCommon.h"',
                            '#include "'+COMMON+'"').encode()
    A.require(b'slot_count' not in common and b'std::array<SymbolSpec,6> slots;' in common,
              'original six/six public binding graph, no serializer variant')
    return common, worker


def prepare(output):
    output = Path(output)
    output = output.parent.resolve(strict=True)/output.name
    A.require(ROOT != output and ROOT not in output.parents and not output.exists() and
              not output.is_symlink(), 'fresh external observer source directory')
    common, worker = worker_sources()
    output.mkdir(mode=0o700)
    A.publish(output/COMMON, common)
    A.publish(output/WORKER, worker)
    A.publish(output/'ADAPTER.json', A.canonical(dict(schema=1, protocol=PROTOCOL,
        sources=[B.pin(output/name) for name in (COMMON, WORKER)],
        helpers=[B.pin(ROOT/'bench'/name) for name in sorted(HELPERS)],
        scientific_launch=False)))
    print('PASS unchanged 38-case WORK and six/six binding observer prepared; no timing')


def verify_adapter(adapter):
    adapter = Path(adapter).resolve(strict=True)
    proof = A.decode(A.read_regular(adapter/'ADAPTER.json', 65536))
    A.exact(set(proof), {'schema','protocol','sources','helpers','scientific_launch'}, 'adapter schema')
    A.exact((proof['schema'], proof['protocol'], proof['scientific_launch']), (1, PROTOCOL, False), 'neutral adapter identity')
    A.exact(proof['sources'], [B.pin(adapter/name) for name in (COMMON, WORKER)], 'exact generated source roster')
    A.exact(proof['helpers'], [B.pin(ROOT/'bench'/name) for name in sorted(HELPERS)], 'exact helper roster')
    for name, expected in zip((COMMON, WORKER), worker_sources()):
        A.exact(A.read_regular(adapter/name, 1024**2), expected, 'complete generated observer source')
    return {adapter/name for name in (COMMON, WORKER, 'ADAPTER.json')}


def library_inputs():
    A.exact(B.pin(PREPARED/'PREPARED.json')['sha256'], PREPARATION_SHA, 'accepted exact native preparation')
    proof = P.verify_prepared(PREPARED)
    inputs = {PREPARED/'PREPARED.json'}
    pins = proof['files'] + proof['captures']
    for path, digest in QUALIFICATIONS:
        A.exact(B.pin(path)['sha256'], digest, 'accepted native/package qualification')
        inputs.add(path)
        qualified = A.decode(A.read_regular(path, 4*1024**2))
        pins += qualified['files']
        if 'tests' in qualified: pins += [qualified['tests'], qualified['junit']]
    for pin in pins:
        path = Path(pin['path'])
        A.exact(B.pin(path), pin, 'unchanged neutral source/build/consumer evidence')
        inputs.add(path)
    return proof, inputs


def check_metadata(meta):
    A.exact(len(meta), 2, 'two actual DSOs')
    for row, (path, digest) in zip(meta, LIBRARIES):
        A.exact((row['path'], row['sha256'], row['context_bytes']), (str(path), digest, 141328),
                'actual prepared provider and native GF context')
        A.exact(tuple(s['name'] for s in row['slots']), SLOTS, 'exact six internal public bindings')
        A.exact(len(row['exports']), 53, 'complete export count')
        A.exact(len(row['runtime_slots']), 37, 'complete runtime count')
    A.exact([s['name'] for s in meta[0]['exports']], [s['name'] for s in meta[1]['exports']],
            'unchanged exact public exports')


def provenance(adapter, proof_dir=None):
    proof, inputs = library_inputs()
    inputs.update(verify_adapter(adapter))
    inputs.update(S.imported_files())
    inputs.update(P.compiler_inputs())
    check_metadata(R.metadata(LIBRARIES))
    result = []
    for index, row in enumerate(proof['links']):
        name = 'proof-'+('old' if index == 0 else 'new')+'.so'
        result.append(dict(original=row['dso'], proof_name=name, proof_sha256=row['dso']['sha256'],
            preparation=B.pin(PREPARED/'PREPARED.json'), adapter=B.pin(Path(adapter)/'ADAPTER.json'),
            qualifications=[B.pin(path) for path, _ in QUALIFICATIONS],
            producing_objects=[arg for arg in row['argv'] if arg.endswith('.o')]))
        if proof_dir is not None:
            # The native producer proof already reproduced all nineteen objects.
            # Relink the exact same roster into this fresh observer build only.
            target, link_map = proof_dir/name, proof_dir/(name+'.map')
            A.require(not target.exists() and not link_map.exists(), 'fresh observer relink proofs')
            argv = list(row['argv'])
            argv[argv.index('-o')+1] = str(target)
            A.require(argv[-1].startswith('-Wl,-Map,'), 'single original map argument')
            argv[-1] = '-Wl,-Map,'+str(link_map)
            B.command(argv)
            S.verify_link_inputs(link_map, inputs)
            A.exact(B.pin(target)['sha256'], row['dso']['sha256'], 'byte-exact nineteen-object relink')
    return result, inputs


def verify_header(header, order, claim, meta, protocol=PROTOCOL):
    check_metadata(meta)
    # First twenty fixtures use pinned immutable neutral data, not an old run
    # or old receipt verifier. The remaining eighteen use independent GF rows.
    F.verify_header(header, order, claim, meta, protocol)


def combine(results, protocol=PROTOCOL):
    result = R.combine(results, protocol)
    result.update(current_path_retention_qualified=result['outcome'] == 'PASS',
                  pre_admission_restoration_qualified=False, ordinary_K8_WH1_speed_qualified=False,
                  historical_K3_workload_retention_qualified=False)
    return result


def qualify(executable, output, meta, mode):
    for order in ('old-new', 'new-old'):
        A.exact(A.read_regular(output/('neutral-'+order+'.txt'), 65536),
                b'PASS neutral98496-coordinate roster, 304 native WORK cases; no timing\n', 'complete gate-A neutral roster')
    D.qualify(executable, output, meta, mode, PROTOCOL, verify_header)


def settings(adapter):
    adapter = Path(adapter).resolve(strict=True)
    verify_adapter(adapter)
    sources = tuple('bench/'+name for name in HELPERS) + R.NEW + (
        'bench/Wh2CurrentPreservedDeferredCostR0.cpp', 'bench/Wh2K8SharedRetention/Adapter.py',
        'bench/Wh2K8SharedRetention/test_Adapter.py', 'bench/Wh2K8OrdinaryShared/Gates.md',
        str(adapter/COMMON), str(adapter/WORKER), str(adapter/'ADAPTER.json'))
    return R.Configuration(PROTOCOL, OUTPUT, LIBRARIES, sources,
        lambda proof_dir=None: provenance(adapter, proof_dir), str(adapter/WORKER), qualify, CASES, verify_header, combine)


def main():
    environment = B.process_environment()
    if dict(os.environ) != environment:
        os.execve(sys.executable, [sys.executable, str(Path(__file__).resolve())]+sys.argv[1:], environment)
        raise RuntimeError('execve returned')
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare'); p.add_argument('output', type=Path)
    b = sub.add_parser('build'); b.add_argument('mode', choices=('native','asan-driver'))
    b.add_argument('adapter', type=Path); b.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.output)
    else: S.build(args.mode, args.output, settings(args.adapter), R)


if __name__ == '__main__': main()
