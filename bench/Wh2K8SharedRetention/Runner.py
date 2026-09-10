#!/usr/bin/env python3
"""Explicit one-shot gate-A controller; never select a historical experiment."""
import argparse
import os
from pathlib import Path
import sys

import Adapter as K

QUALIFIED = Path('/tmp/wh2-k8-shared-retention-final.4O8Fkmso')
ADAPTER = QUALIFIED/'adapter'
MANIFESTS = (
    ('native', 'a7a15094c55f8d3f8e0ca6d9cef90edc3f0841ba3ae209f63fd5b3dfcab444b1'),
    ('asan-driver', '1441cf2bce3ab4ab0f4cbb051d528808d2cc304247006b3b8f490bb9fc546906'))
NEW = ('bench/Wh2K8SharedRetention/Runner.py',
       'bench/Wh2K8SharedRetention/test_Runner.py',
       'bench/Wh2K8SharedRetention/Freeze.md')


def observer_inputs():
    """Bind both accepted neutral observers, never their scientific controls."""
    inputs, pins = set(), {}
    for mode, digest in MANIFESTS:
        path = QUALIFIED/mode/'manifest.json'
        K.A.exact(K.B.pin(path)['sha256'], digest, 'accepted neutral observer manifest')
        manifest = K.A.decode(K.A.read_regular(path, 1024**2))
        K.A.exact((manifest['protocol'], manifest['mode'], manifest['scientific_launch'],
                   manifest['library_source_provenance_closed'], manifest['sanitized_library_code']),
                  (K.PROTOCOL, mode, False, True, False), 'accepted observer scope')
        inputs.add(path)
        for pin in manifest['inputs'] + manifest['artifacts']:
            if pin['path'] in pins:
                K.A.exact(pins[pin['path']], pin, 'consistent neutral prerequisite overlap')
            pins[pin['path']] = pin
            K.A.exact(K.B.pin(pin['path']), pin, 'unchanged neutral prerequisite')
            inputs.add(Path(pin['path']))
    return inputs


def provenance(proof_dir=None):
    result, inputs = K.provenance(ADAPTER, proof_dir)
    inputs.update(observer_inputs())
    return result, inputs


def settings():
    cfg = K.settings(ADAPTER)
    return cfg._replace(sources=cfg.sources + NEW, provenance=provenance)


def enter_clean_environment():
    environment = K.B.process_environment()
    if dict(os.environ) != environment:
        os.execve(sys.executable, [sys.executable, str(Path(__file__).resolve())]+sys.argv[1:], environment)
        raise RuntimeError('execve returned')


def main(argv=None):
    enter_clean_environment()
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('mode', choices=('native','asan-driver'))
    b.add_argument('output', type=Path)
    r = sub.add_parser('receipt'); r.add_argument('build_dir', type=Path); r.add_argument('output', type=Path)
    r = sub.add_parser('run'); r.add_argument('receipt', type=Path)
    sub.add_parser('replay')
    args = parser.parse_args(argv)
    cfg = settings()
    if args.command == 'build': K.S.build(args.mode, args.output, cfg, K.R)
    elif args.command == 'receipt': K.A.publish(args.output, K.A.canonical(K.R.receipt(args.build_dir, cfg)))
    elif args.command == 'run': K.R.run(args.receipt, cfg)
    else:
        result = K.R.replay(cfg)
        print(K.A.canonical(dict(outcome=result['outcome'], exact_replay=True)).decode(), end='')


if __name__ == '__main__': main()
