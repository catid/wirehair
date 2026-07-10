#!/usr/bin/env python3
"""Native-library E2E for the installed ctypes binding."""

import argparse
import gc
import os
from pathlib import Path
import sys
from unittest import mock
import weakref


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


def expect_memory_error(function, required_text):
    try:
        function()
    except MemoryError as error:
        if required_text not in str(error):
            raise AssertionError(
                "allocation error did not contain %r: %s" %
                (required_text, error))
    else:
        raise AssertionError("expected injected allocation failure")


def complete_systematic_decoder(wh, library, encoder, original, block_bytes):
    decoder = wh.Decoder.create(
        len(original), block_bytes, library=library,
        profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT)
    block_count = (len(original) + block_bytes - 1) // block_bytes
    complete = False
    for block_id in range(block_count):
        complete = decoder.decode(block_id, encoder.encode(block_id))
    if not complete or decoder.recover() != original:
        decoder.close()
        raise AssertionError("systematic native decoder did not recover")
    return decoder


def exercise_conversion_failure_atomicity(
        wh, library, encoder, original, block_bytes):
    """Inject Python allocation/commit failures around the real native call."""
    native_become = library.wirehair_decoder_becomes_encoder
    native_release = wh._release_codec
    releases = {}

    def tracked_release(native_library, pointer_value):
        releases[pointer_value] = releases.get(pointer_value, 0) + 1
        return native_release(native_library, pointer_value)

    original_finish = wh._Codec._finish_construction

    def fail_after_finalizer(instance):
        original_finish(instance)
        raise MemoryError("injected after-finalizer failure")

    failure_patches = (
        ("new failure", lambda: mock.patch.object(
            wh.Encoder, "__new__",
            side_effect=MemoryError("injected new failure"))),
        ("init failure", lambda: mock.patch.object(
            wh.Encoder, "__init__",
            side_effect=MemoryError("injected init failure"))),
        ("finalizer failure", lambda: mock.patch.object(
            wh.weakref, "finalize",
            side_effect=MemoryError("injected finalizer failure"))),
        ("after-finalizer failure", lambda: mock.patch.object(
            wh._Codec, "_finish_construction", new=fail_after_finalizer)),
    )

    with mock.patch.object(wh, "_release_codec", new=tracked_release):
        for required_text, failure_patch in failure_patches:
            decoder = complete_systematic_decoder(
                wh, library, encoder, original, block_bytes)
            pointer_value = decoder.pointer.value
            release_before = releases.get(pointer_value, 0)
            conversion_calls = []

            def counting_become(pointer):
                conversion_calls.append(pointer.value)
                return native_become(pointer)

            with mock.patch.object(
                    library, "wirehair_decoder_becomes_encoder",
                    new=counting_become), failure_patch():
                expect_memory_error(decoder.become_encoder, required_text)
            gc.collect()
            if decoder.closed or decoder.recover() != original:
                raise AssertionError(
                    "%s did not preserve a usable native decoder" %
                    required_text)
            if conversion_calls:
                raise AssertionError(
                    "%s reached irreversible native conversion" %
                    required_text)
            if releases.get(pointer_value, 0) != release_before:
                raise AssertionError("pre-conversion failure freed the decoder")
            decoder.close()
            gc.collect()
            if releases.get(pointer_value, 0) != release_before + 1:
                raise AssertionError(
                    "%s did not free exactly once" % required_text)

        decoder = wh.Decoder.create(
            len(original), block_bytes, library=library,
            profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT)
        pointer_value = decoder.pointer.value
        release_before = releases.get(pointer_value, 0)
        expect_wirehair_error(
            wh, "wirehair_decoder_becomes_encoder", wh.Wirehair_InvalidInput,
            decoder.become_encoder)
        if decoder.closed:
            raise AssertionError("precondition failure closed the native decoder")
        block_count = (len(original) + block_bytes - 1) // block_bytes
        complete = False
        for block_id in range(block_count):
            complete = decoder.decode(block_id, encoder.encode(block_id))
        if not complete or decoder.recover() != original:
            raise AssertionError(
                "decoder was not usable after native precondition failure")
        decoder.close()
        if releases.get(pointer_value, 0) != release_before + 1:
            raise AssertionError("precondition failure cleanup was not exact")

        decoder = complete_systematic_decoder(
            wh, library, encoder, original, block_bytes)
        pointer_value = decoder.pointer.value
        release_before = releases.get(pointer_value, 0)

        def convert_then_raise(pointer):
            result = native_become(pointer)
            if result != wh.Wirehair_Success:
                raise AssertionError(
                    "native conversion failed before injected exception: %d" %
                    result)
            raise MemoryError("injected post-conversion failure")

        with mock.patch.object(
                library, "wirehair_decoder_becomes_encoder",
                new=convert_then_raise):
            expect_memory_error(decoder.become_encoder, "post-conversion")
        if (not decoder.closed or
                releases.get(pointer_value, 0) != release_before + 1):
            raise AssertionError(
                "post-conversion exception left a live or multiply-owned codec")
        try:
            decoder.recover()
        except RuntimeError:
            pass
        else:
            raise AssertionError("post-conversion source remained callable")
        decoder.close()
        gc.collect()
        if releases.get(pointer_value, 0) != release_before + 1:
            raise AssertionError("post-conversion cleanup freed more than once")

        decoder = complete_systematic_decoder(
            wh, library, encoder, original, block_bytes)
        pointer_value = decoder.pointer.value
        release_before = releases.get(pointer_value, 0)

        def fail_commit():
            raise MemoryError("injected ownership commit failure")

        decoder._relinquish = fail_commit
        expect_memory_error(decoder.become_encoder, "ownership commit")
        if (not decoder.closed or
                releases.get(pointer_value, 0) != release_before + 1):
            raise AssertionError(
                "commit failure left the converted source apparently usable")
        try:
            decoder.recover()
        except RuntimeError:
            pass
        else:
            raise AssertionError("commit-failure source remained callable")
        decoder.close()
        gc.collect()
        if releases.get(pointer_value, 0) != release_before + 1:
            raise AssertionError("commit-failure cleanup freed more than once")

        decoder = complete_systematic_decoder(
            wh, library, encoder, original, block_bytes)
        pointer_value = decoder.pointer.value
        release_before = releases.get(pointer_value, 0)
        converted = decoder.become_encoder()
        converted_reference = weakref.ref(converted)
        del converted
        gc.collect()
        if converted_reference() is not None:
            raise AssertionError("converted encoder was not finalized")
        if (not decoder.closed or
                releases.get(pointer_value, 0) != release_before + 1):
            raise AssertionError("converted encoder GC did not free exactly once")
        decoder.close()


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
    expect_wirehair_error(
        wh,
        "wirehair_wire_profile_init",
        wh.Wirehair_InvalidInput,
        lambda: wh.Encoder.create(
            b"profile-check", 4, library=library,
            profile_id=0x123456789abcdef0),
    )
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

    with wh.Encoder.create(
            mutable, block_bytes, library=library, owned=True,
            profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT) as encoder:
        mutable[:] = b"\0" * len(mutable)
        with wh.Decoder.create(
                len(original), block_bytes, library=library,
                profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT) as decoder:
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

        exercise_conversion_failure_atomicity(
            wh, library, encoder, original, block_bytes)

    # Exercise the native borrowed-buffer path and the historical equation
    # profile with repair packets only; systematic packets would not
    # distinguish profiles.
    pre_source = bytearray(
        (index * 41 + index // 7 + 5) & 0xff for index in range(257))
    pre_expected = bytes(pre_source)
    pre_block_bytes = 31
    pre_block_count = (
        len(pre_source) + pre_block_bytes - 1) // pre_block_bytes
    with wh.Encoder.create(
            pre_source, pre_block_bytes, library=library, owned=False,
            profile_id=wh.WIREHAIR_LEGACY_PROFILE_PRE_FIXUP) as pre_encoder:
        with wh.Decoder.create(
                len(pre_source), pre_block_bytes, library=library,
                profile_id=wh.WIREHAIR_LEGACY_PROFILE_PRE_FIXUP) as pre_decoder:
            pre_complete = False
            for packet_id in range(
                    pre_block_count, pre_block_count * 2 + 128):
                if pre_decoder.decode(
                        packet_id, pre_encoder.encode(packet_id)):
                    pre_complete = True
                    break
            if not pre_complete or pre_decoder.recover() != pre_expected:
                raise AssertionError(
                    "borrowed PRE_FIXUP repair-only round trip failed")

    print(
        "native Python E2E passed: bytes=%d blocks=%d delivered=%d library=%s" %
        (len(original), block_count, delivered, library_path)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
