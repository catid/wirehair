#!/usr/bin/env python3
"""K12 one-shot launcher and authenticated structural-screen projection."""
import importlib.util
import hashlib
from pathlib import Path
import sys

SPEC = importlib.util.spec_from_file_location('k12_tm_controller',
    Path(__file__).with_name('Wh2NoncommutingRadixRunR0.py'))
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
C.PROTOCOL = 'wirehair.wh2.k12-thue-morse-r1'
C.OUTPUT = Path('/var/tmp/wh2-k12-thue-morse-r1')
C.SOURCES = ('bench/Wh2K12ThueMorseR0.py', 'bench/Wh2K12ThueMorseRunR0.py',
             'bench/test_Wh2K12ThueMorseR0.py', 'bench/Wh2K12ThueMorseR0.md')
C.INVENTORY = Path('/var/tmp/wh2-uncovered-band-inventory-r0')
C.MANIFEST = '16fcf13214cd25362fe66ee35f62c1c9616ed082e343fa5e2009d81d61359ea0'
C.PROJECTIONS = {
    'origins': (63, '21c1ffe5cf89d8bc99846b51ff50d7e9d7b9a7bd7c723a2ab73cfc4a7f38b2ec'),
    'prefixes': (60, '9ba913961b79f21992b01c7f9750d1604307dca655258a989c4342601c3a2026'),
    'roots': (64, 'a3bc3986862980d0698a9d5b2850879d113ea2a2c94bc237392926b607707e03')}

K12_SPEC = importlib.util.spec_from_file_location('k12_screen_inputs',
    Path(__file__).with_name('Wh2K12ThueMorseR0.py'))
K12 = importlib.util.module_from_spec(K12_SPEC)
K12_SPEC.loader.exec_module(K12)


BASE = C.current_receipt


def current_receipt(deadline=None):
    receipt = BASE(deadline)
    projection = K12.history_inputs(deadline)
    history = dict(provenance=dict(path=str(C.INVENTORY), manifest_sha256=C.MANIFEST),
                   projection_sha256=hashlib.sha256(K12.canonical(projection)).hexdigest())
    receipt['history'] = history
    return receipt


C.current_receipt = current_receipt


def main():
    return C.main()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ': ' + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
