#!/usr/bin/env python3
"""Native-library E2E for the installed ctypes binding."""

import argparse
import os
from pathlib import Path
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    return parser.parse_args()


def expect_wirehair_error(wh, operation, result, function):
    try:
        function()
    except wh.WirehairError as error:
        if error.operation != operation or error.result != result:
            raise AssertionError(
                "expected %s result %d, got %s result %d" %
                (operation, result, error.operation, error.result))
    else:
        raise AssertionError("expected %s to fail" % operation)


def expect_value_error(function, required_text):
    try:
        function()
    except ValueError as error:
        if required_text not in str(error):
            raise AssertionError(
                "validation error did not contain %r: %s" %
                (required_text, error))
    else:
        raise AssertionError("expected wrapper validation to fail")


def main():
    args = parse_args()
    module_dir = args.module_dir.resolve()
    library_path = args.library.resolve()
    if not (module_dir / "whirehair.py").is_file():
        raise RuntimeError("installed binding not found in " + str(module_dir))
    if not library_path.is_file():
        raise RuntimeError("native library not found: " + str(library_path))

    sys.path.insert(0, str(module_dir))
    os.environ["WIREHAIR_LIBRARY"] = str(library_path)
    import whirehair as wh

    library = wh.initialize()
    original = bytes((index * 73 + index // 11 + 19) & 0xff for index in range(4099))
    mutable = bytearray(original)
    block_bytes = 113
    block_count = (len(original) + block_bytes - 1) // block_bytes

    expect_wirehair_error(
        wh,
        "wirehair_encoder_create",
        wh.Wirehair_BadInput_SmallN,
        lambda: wh.Encoder.create(b"x" * block_bytes, block_bytes, library=library),
    )
    expect_wirehair_error(
        wh,
        "wirehair_decoder_create",
        wh.Wirehair_BadInput_SmallN,
        lambda: wh.Decoder.create(block_bytes, block_bytes, library=library),
    )

    with wh.Encoder.create(mutable, block_bytes, library=library, owned=True) as encoder:
        mutable[:] = b"\0" * len(mutable)
        with wh.Decoder.create(len(original), block_bytes, library=library) as decoder:
            expect_wirehair_error(
                wh,
                "wirehair_encode",
                wh.Wirehair_InvalidInput,
                lambda: encoder.encode(0, capacity=block_bytes - 1),
            )
            expect_wirehair_error(
                wh,
                "wirehair_recover",
                wh.Wirehair_NeedMore,
                decoder.recover,
            )
            expect_wirehair_error(
                wh,
                "wirehair_recover_block",
                wh.Wirehair_NeedMore,
                lambda: decoder.recover_block(0),
            )
            expect_value_error(
                lambda: decoder.decode(0, b"x" * (block_bytes + 1)),
                "exceeds block_bytes",
            )
            ordered_ids = []
            for first in range(0, block_count + 128, 5):
                ordered_ids.extend(reversed(range(first, min(first + 5, block_count + 128))))
            packet_ids = [block_id for block_id in ordered_ids
                          if (block_id * 7 + 3) % 13 >= 4]
            if packet_ids == sorted(packet_ids) or len(packet_ids) == len(ordered_ids):
                raise AssertionError("loss/reorder schedule did not exercise both behaviors")

            complete = False
            delivered = 0
            for block_id in packet_ids:
                delivered += 1
                if decoder.decode(block_id, encoder.encode(block_id)):
                    complete = True
                    break
            if not complete:
                raise AssertionError("native decoder did not complete within the bounded schedule")
            if decoder.recover() != original:
                raise AssertionError("native full-message recovery mismatch")

            for block_id in (0, block_count // 2, block_count - 1):
                expected = original[block_id * block_bytes:(block_id + 1) * block_bytes]
                actual = decoder.recover_block(block_id, capacity=len(expected))
                if actual != expected:
                    raise AssertionError("native block recovery mismatch at %d" % block_id)

            converted = decoder.become_encoder()
            try:
                if not decoder.closed:
                    raise AssertionError("source decoder remained open after conversion")
                try:
                    decoder.recover()
                except RuntimeError as error:
                    if "closed" not in str(error):
                        raise
                else:
                    raise AssertionError("closed source decoder remained usable")
                for block_id in (0, block_count - 1, block_count, block_count + 17):
                    if converted.encode(block_id) != encoder.encode(block_id):
                        raise AssertionError(
                            "decoder-to-encoder conversion mismatch at %d" % block_id)
            finally:
                converted.close()

    print(
        "native Python E2E passed: bytes=%d blocks=%d delivered=%d library=%s" %
        (len(original), block_count, delivered, library_path)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
