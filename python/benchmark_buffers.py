#!/usr/bin/env python3
"""Measure Python output allocation and throughput for reusable buffers."""

import argparse
import gc
import os
from pathlib import Path
import sys
import time
import tracemalloc


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--module-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--message-bytes", type=int, default=262144)
    parser.add_argument("--block-bytes", type=int, default=1200)
    parser.add_argument("--iterations", type=int, default=2000)
    return parser.parse_args()


def measure(name, method, calls, logical_bytes, operation):
    # Warm ctypes' array-type cache so peak numbers describe steady-state calls.
    operation(0)
    gc.collect()
    tracemalloc.start()
    start = time.perf_counter()
    checksum = 0
    for index in range(calls):
        checksum ^= operation(index)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    mib_per_second = calls * logical_bytes / elapsed / (1024 * 1024)
    return (name, method, calls, elapsed, calls / elapsed,
            mib_per_second, peak, checksum)


def main():
    args = parse_args()
    if args.message_bytes < 2 or args.block_bytes < 1 or args.iterations < 1:
        raise ValueError("message, block, and iteration counts must be positive")
    block_count = (
        args.message_bytes + args.block_bytes - 1) // args.block_bytes
    if block_count < 2 or block_count > 64000:
        raise ValueError("benchmark requires 2..64000 source blocks")

    sys.path.insert(0, str(args.module_dir.resolve()))
    if args.library is not None:
        os.environ["WIREHAIR_LIBRARY"] = str(args.library.resolve())
    import whirehair as wh

    library = wh.initialize()
    message = bytes(
        (index * 73 + index // 11 + 19) & 0xff
        for index in range(args.message_bytes))
    results = []

    with wh.Encoder.create(
            message, args.block_bytes, library=library, owned=True) as encoder:
        decoder = wh.Decoder.create(
            len(message), args.block_bytes, library=library)
        packet = bytearray(args.block_bytes)
        complete = False
        for block_id in range(block_count):
            written = encoder.encode_into(block_id, packet)
            complete = decoder.decode(block_id, memoryview(packet)[:written])
        if not complete:
            decoder.close()
            raise RuntimeError("systematic benchmark decode did not complete")

        repair_id = block_count + 100
        encode_output = bytearray(args.block_bytes)

        def encode_alloc(index):
            value = encoder.encode(repair_id + index)
            return value[0] ^ value[-1]

        def encode_into(index):
            written = encoder.encode_into(repair_id + index, encode_output)
            return written ^ encode_output[0] ^ encode_output[written - 1]

        results.append(measure(
            "encode", "bytes", args.iterations, args.block_bytes,
            encode_alloc))
        results.append(measure(
            "encode", "into", args.iterations, args.block_bytes,
            encode_into))

        recovery_calls = max(8, min(200, args.iterations // 20))
        recovery_output = bytearray(len(message))

        def recover_alloc(_index):
            value = decoder.recover()
            return value[0] ^ value[-1]

        def recover_into(_index):
            written = decoder.recover_into(recovery_output)
            return written ^ recovery_output[0] ^ recovery_output[-1]

        results.append(measure(
            "recover", "bytes", recovery_calls, len(message),
            recover_alloc))
        results.append(measure(
            "recover", "into", recovery_calls, len(message),
            recover_into))

        block_output = bytearray(args.block_bytes)
        block_calls = args.iterations

        def block_alloc(index):
            value = decoder.recover_block(index % block_count)
            return value[0] ^ value[-1]

        def block_into(index):
            block_id = index % block_count
            written = decoder.recover_block_into(block_id, block_output)
            return written ^ block_output[0] ^ block_output[written - 1]

        results.append(measure(
            "recover_block", "bytes", block_calls, args.block_bytes,
            block_alloc))
        results.append(measure(
            "recover_block", "into", block_calls, args.block_bytes,
            block_into))
        decoder.close()

    print(
        "operation       method calls   seconds    calls/s      MiB/s "
        "traced_peak_B checksum")
    for result in results:
        print("%-15s %-6s %5d %9.4f %10.0f %10.1f %13d %8d" % result)
    print(
        "traced_peak_B excludes reusable buffers allocated before tracing; "
        "it includes steady-state Python/ctypes call scaffolding.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
