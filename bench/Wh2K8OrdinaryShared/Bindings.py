#!/usr/bin/env python3
"""Neutral shared-loader binding checks; no timing or loss cohort."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import sys

import Prepare as P

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('_k8_shared_binding_helpers', HERE.parent/'Wh2AdmissionRegressionCostR0.py')
R = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = R
spec.loader.exec_module(R)
A = R.A


def exercise(lib, candidate):
    """Small deterministic API checks after both DSOs are loaded, not loss tests."""
    N, T = R.N, R.N.T
    bare = lib.call('wirehair_v2_encoder_create', N.INT, N.V, N.U64, N.U32, N.V, N.U32,
                    T.POINTER(N.U32), T.POINTER(N.V))
    options = lib.call('wirehair_v2_encoder_create_with_options', N.INT, N.V, N.U64, N.U32,
                       T.POINTER(N.Options), N.V, N.U32, T.POINTER(N.U32), T.POINTER(N.V))
    encode = lib.call('wirehair_v2_encode', N.INT, N.V, N.U32, N.V, N.U32, T.POINTER(N.U32))
    create_decoder = lib.call('wirehair_v2_decoder_create', N.INT, N.V, N.U32, T.POINTER(N.V))
    decode = lib.call('wirehair_v2_decode', N.INT, N.V, N.U32, N.V, N.U32)
    recover = lib.call('wirehair_v2_recover', N.INT, N.V, N.V, N.U64, T.POINTER(N.U64))
    free = lib.call('wirehair_v2_free', None, N.V)
    count = 0
    for k in (3, 5, 8, 16):
        for block in (2, 64):
            for tail in (1, block):
                for policy in (0, 1, 2):
                    message = bytes((37*i+i//11) % 256 for i in range((k-1)*block+tail))
                    source = T.create_string_buffer(message)
                    descriptor = (T.c_ubyte*32)()
                    written, encoder, decoder = N.U32(), N.V(), N.V()
                    try:
                        if policy:
                            opts = N.Options(16, 1, policy, 0)
                            status = options(source, len(message), block, T.byref(opts), descriptor,
                                             32, T.byref(written), T.byref(encoder))
                        else:
                            status = bare(source, len(message), block, descriptor, 32,
                                          T.byref(written), T.byref(encoder))
                        A.require(status == 0 and encoder.value and written.value == 32, 'actual shared ordinary creation')
                        expected = 0x67c1043ecaa9e184 if k == 3 else (
                            0x7a9276b85c730ae0 if candidate and k == 8 else 0x4b295bbb47f4f9c9)
                        A.exact(int.from_bytes(bytes(descriptor)[8:16], 'little'), expected,
                                'bare/options internal binding selects this DSO profile')
                        packets = []
                        for packet_id in range(k):
                            packet = (T.c_ubyte*block)()
                            A.exact(encode(encoder, packet_id, packet, block, T.byref(written)), 0, 'shared systematic encode')
                            expected_bytes = tail if packet_id == k-1 else block
                            A.exact(written.value, expected_bytes, 'systematic payload length')
                            A.exact(bytes(packet)[:expected_bytes], message[packet_id*block:packet_id*block+expected_bytes],
                                    'independent systematic payload oracle')
                            packets.append((packet, expected_bytes))
                        free(encoder); encoder = N.V()
                        T.memset(source, 0xcc, len(message))
                        A.require(create_decoder(descriptor, 32, T.byref(decoder)) == 0 and decoder.value,
                                  'receiver created after sender destruction')
                        for packet_id, (packet, size) in enumerate(packets):
                            A.exact(decode(decoder, packet_id, packet, size), 0 if packet_id == k-1 else 1,
                                    'shared systematic first success')
                        recovered = (T.c_ubyte*(len(message)+2))()
                        recovered_bytes = N.U64()
                        for _ in range(2):
                            T.memset(recovered, 0xa5, len(recovered))
                            A.exact(recover(decoder, T.byref(recovered, 1), len(message), T.byref(recovered_bytes)),
                                    0, 'shared repeated recovery')
                            A.exact(recovered_bytes.value, len(message), 'shared recovered length')
                            A.exact(bytes(recovered), b'\xa5'+message+b'\xa5', 'guarded shared recovery bytes')
                        count += 1
                    finally:
                        free(encoder); free(decoder)
    A.exact(count, 48, 'complete two-library API roster')
    return count


def check(prepared, reverse):
    A.require(not any(key in os.environ for key in R.ENV_KEYS), 'no dynamic-loader/allocator overrides')
    proof = P.verify_prepared(prepared)
    libraries = tuple((Path(row['dso']['path']), row['dso']['sha256']) for row in proof['links'])
    A.exact(len(libraries), 2, 'two prepared native DSOs')
    metadata = R.metadata(libraries)
    loaded = [None, None]
    for index in ((1, 0) if reverse else (0, 1)):
        loaded[index] = R.N.Library(index, libraries)
    A.require(loaded[0].base != loaded[1].base and
              loaded[0].report['context'] != loaded[1].report['context'], 'distinct libraries/GF contexts')
    runtime = []
    for lib, meta in zip(loaded, metadata):
        targets = []
        for row in meta['runtime_slots']:
            observed = R.N.V.from_address(lib.base+row['offset']).value
            expected = R.N.address(getattr(lib.dso, row['name']))
            A.exact(observed, expected, 'actual external runtime binding '+row['name'])
            targets.append(dict(name=row['name'], address=observed))
        runtime.append(targets)
        lib.check_bindings()
        A.exact(lib.report['context_bytes'], 141328, 'native GF context size')
    A.exact(runtime[0], runtime[1], 'common actual runtime providers')
    cases = [exercise(lib, index == 1) for index, lib in enumerate(loaded)]
    for lib in loaded:
        lib.check_bindings()
    P.verify_prepared(prepared)
    return dict(scope='neutral native shared bindings only; no timing', reverse=reverse,
                metadata=metadata, bindings=[lib.report for lib in loaded], runtime=runtime, cases=cases)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared', type=Path)
    parser.add_argument('--reverse', action='store_true')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2, 512*1024**2))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    signal.alarm(30)
    print(json.dumps(check(args.prepared.resolve(strict=True), args.reverse), sort_keys=True))
