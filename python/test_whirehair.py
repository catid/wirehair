#!/usr/bin/env python3
"""Unit tests for the ctypes Wirehair binding."""

import ctypes
import gc
import io
import os
import sys
import tempfile
import threading
import unittest
import weakref
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))
import whirehair as wh
import wirehair as canonical


class FakeFunction:
    def __init__(self, function):
        self.function = function
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.function(*args)


def number(value):
    return value.value if hasattr(value, "value") else int(value)


class FakeLibrary:
    def __init__(self):
        self.next_pointer = 100
        self.codecs = {}
        self.freed = []
        self.init_result = wh.Wirehair_Success
        self.encoder_result = wh.Wirehair_Success
        self.decoder_result = wh.Wirehair_Success
        self.encode_result = wh.Wirehair_Success
        self.detach_result = wh.Wirehair_Success
        self.detach_calls = 0
        self.decode_terminal = None
        self.always_need_more = False
        self.recover_result = wh.Wirehair_Success
        self.become_encoder_result = wh.Wirehair_Success
        self.become_encoder_calls = 0
        self.raise_after_become_encoder = False
        self.corrupt_recovery = False
        self.raise_encode = False
        self.raise_recover = False
        self.last_recover_capacity = None
        self.encode_calls = 0
        self.recover_calls = 0
        self.recover_block_calls = 0
        self.before_encode = None

        for name in (
            "wirehair_result_string", "wirehair_init_", "wirehair_encoder_create",
            "wirehair_wire_profile_init",
            "wirehair_encoder_create_ex", "wirehair_encoder_create_owned",
            "wirehair_encoder_create_owned_ex",
            "wirehair_encoder_create_profile_ex", "wirehair_encode",
            "wirehair_encoder_detach_input",
            "wirehair_decoder_create", "wirehair_decoder_create_ex",
            "wirehair_decoder_create_profile_ex",
            "wirehair_decode", "wirehair_recover", "wirehair_recover_block",
            "wirehair_recover_block_ex", "wirehair_decoder_becomes_encoder",
            "wirehair_free",
        ):
            setattr(self, name, FakeFunction(getattr(self, "_" + name)))

    def _new(self, state, output=None):
        pointer = self.next_pointer
        self.next_pointer += 1
        self.codecs[pointer] = state
        if output is not None:
            ctypes.cast(output, ctypes.POINTER(ctypes.c_void_p))[0] = pointer
        return pointer

    def _state(self, pointer):
        return self.codecs[number(pointer)]

    def _wirehair_result_string(self, result):
        return ("fake result %d" % number(result)).encode("ascii")

    def _wirehair_init_(self, _version):
        return self.init_result

    def _wirehair_wire_profile_init(self, profile_id, output):
        profile_id = number(profile_id)
        if profile_id not in (
                wh.WIREHAIR_LEGACY_PROFILE_PRE_FIXUP,
                wh.WIREHAIR_LEGACY_PROFILE_CURRENT):
            return wh.Wirehair_InvalidInput
        profile = ctypes.cast(
            output, ctypes.POINTER(wh.WirehairWireProfile))[0]
        profile.struct_bytes = ctypes.sizeof(wh.WirehairWireProfile)
        profile.profile_version = wh.WIREHAIR_WIRE_PROFILE_VERSION
        profile.profile_id = profile_id
        return wh.Wirehair_Success

    def _encoder_state(self, data, size, block, owned):
        size = number(size)
        state = {"kind": "encoder", "size": size, "block": number(block)}
        if owned:
            state["message"] = ctypes.string_at(data, size)
        else:
            state["source"] = number(data)
        return state

    def _wirehair_encoder_create(self, _reuse, data, size, block):
        if self.encoder_result != wh.Wirehair_Success:
            return None
        return self._new(self._encoder_state(data, size, block, False))

    def _wirehair_encoder_create_ex(self, _reuse, data, size, block, output):
        if self.encoder_result != wh.Wirehair_Success:
            return self.encoder_result
        self._new(self._encoder_state(data, size, block, False), output)
        return wh.Wirehair_Success

    def _wirehair_encoder_create_owned(self, _reuse, data, size, block):
        if self.encoder_result != wh.Wirehair_Success:
            return None
        return self._new(self._encoder_state(data, size, block, True))

    def _wirehair_encoder_create_owned_ex(self, _reuse, data, size, block, output):
        if self.encoder_result != wh.Wirehair_Success:
            return self.encoder_result
        self._new(self._encoder_state(data, size, block, True), output)
        return wh.Wirehair_Success

    def _wirehair_encoder_create_profile_ex(
            self, _reuse, data, size, block, _profile, flags, output):
        if self.encoder_result != wh.Wirehair_Success:
            return self.encoder_result
        self._new(self._encoder_state(
            data, size, block,
            bool(number(flags) & wh.WIREHAIR_ENCODER_OWN_INPUT)), output)
        return wh.Wirehair_Success

    def _encoder_message(self, state):
        if "message" in state:
            return state["message"]
        return ctypes.string_at(state["source"], state["size"])

    def _wirehair_encode(self, pointer, block_id, output, capacity, written):
        self.encode_calls += 1
        if self.before_encode is not None:
            self.before_encode()
        if self.raise_encode:
            raise RuntimeError("injected encode exception")
        if self.encode_result != wh.Wirehair_Success:
            return self.encode_result
        state = self._state(pointer)
        message = self._encoder_message(state)
        block = state["block"]
        block_id = number(block_id)
        start = block_id * block
        if start < len(message):
            packet = message[start:start + block]
        else:
            packet = (message + b"\0" * block)[:block]
        packet = packet[:number(capacity)]
        ctypes.memmove(output, packet, len(packet))
        ctypes.cast(written, ctypes.POINTER(ctypes.c_uint32))[0] = len(packet)
        return wh.Wirehair_Success

    def _wirehair_encoder_detach_input(self, pointer):
        self.detach_calls += 1
        if self.detach_result != wh.Wirehair_Success:
            return self.detach_result
        state = self._state(pointer)
        if state.get("kind") != "encoder":
            return wh.Wirehair_InvalidInput
        state["message"] = self._encoder_message(state)
        state.pop("source", None)
        return wh.Wirehair_Success

    def _wirehair_decoder_create(self, _reuse, size, block):
        if self.decoder_result != wh.Wirehair_Success:
            return None
        return self._new({"kind": "decoder", "size": number(size),
                          "block": number(block), "packets": {}})

    def _wirehair_decoder_create_ex(self, _reuse, size, block, output):
        if self.decoder_result != wh.Wirehair_Success:
            return self.decoder_result
        self._new({"kind": "decoder", "size": number(size),
                   "block": number(block), "packets": {}}, output)
        return wh.Wirehair_Success

    def _wirehair_decoder_create_profile_ex(
            self, reuse, size, block, _profile, output):
        return self._wirehair_decoder_create_ex(reuse, size, block, output)

    def _wirehair_decode(self, pointer, block_id, data, size):
        if self.decode_terminal is not None:
            return self.decode_terminal
        state = self._state(pointer)
        state["packets"][number(block_id)] = ctypes.string_at(data, number(size))
        if self.always_need_more:
            return wh.Wirehair_NeedMore
        needed = (state["size"] + state["block"] - 1) // state["block"]
        return (wh.Wirehair_Success if len(state["packets"]) >= needed else
                wh.Wirehair_NeedMore)

    def _original_message(self, state):
        encoders = [codec for codec in self.codecs.values()
                    if codec.get("kind") == "encoder"]
        if encoders:
            return self._encoder_message(encoders[0])
        chunks = [state["packets"][key] for key in sorted(state["packets"])]
        return b"".join(chunks)[:state["size"]]

    def _wirehair_recover(self, pointer, output, size):
        self.recover_calls += 1
        if self.raise_recover:
            raise RuntimeError("injected recover exception")
        if self.recover_result != wh.Wirehair_Success:
            return self.recover_result
        state = self._state(pointer)
        message = self._original_message(state)[:number(size)]
        if self.corrupt_recovery and message:
            message = bytes([message[0] ^ 1]) + message[1:]
        ctypes.memmove(output, message, len(message))
        return wh.Wirehair_Success

    def _wirehair_recover_block(self, pointer, block_id, output, written):
        self.recover_block_calls += 1
        state = self._state(pointer)
        return self._recover_block(state, block_id, output, state["block"], written)

    def _wirehair_recover_block_ex(self, pointer, block_id, output, capacity, written):
        self.recover_block_calls += 1
        state = self._state(pointer)
        self.last_recover_capacity = number(capacity)
        return self._recover_block(state, block_id, output, number(capacity), written)

    def _recover_block(self, state, block_id, output, capacity, written):
        start = number(block_id) * state["block"]
        packet = self._original_message(state)[start:start + state["block"]]
        if capacity < len(packet):
            ctypes.cast(written, ctypes.POINTER(ctypes.c_uint32))[0] = 0
            return wh.Wirehair_InvalidInput
        ctypes.memmove(output, packet, len(packet))
        ctypes.cast(written, ctypes.POINTER(ctypes.c_uint32))[0] = len(packet)
        return wh.Wirehair_Success

    def _wirehair_decoder_becomes_encoder(self, pointer):
        self.become_encoder_calls += 1
        if self.become_encoder_result != wh.Wirehair_Success:
            return self.become_encoder_result
        state = self._state(pointer)
        state["message"] = self._original_message(state)
        state["kind"] = "encoder"
        if self.raise_after_become_encoder:
            raise MemoryError("injected post-conversion failure")
        return wh.Wirehair_Success

    def _wirehair_free(self, pointer):
        value = number(pointer)
        self.freed.append(value)
        self.codecs.pop(value, None)


class BindingTests(unittest.TestCase):
    def setUp(self):
        self.library = FakeLibrary()

    def completed_decoder(self, message=b"abcdefgh", block_bytes=4):
        decoder = wh.Decoder.create(len(message), block_bytes, self.library)
        block_count = (len(message) + block_bytes - 1) // block_bytes
        for block_id in range(block_count):
            start = block_id * block_bytes
            decoder.decode(block_id, message[start:start + block_bytes])
        self.assertEqual(message, decoder.recover())
        return decoder

    def assert_decoder_usable(self, decoder, message=b"abcdefgh"):
        self.assertFalse(decoder.closed)
        self.assertEqual("decoder", self.library._state(decoder.pointer)["kind"])
        self.assertEqual(message, decoder.recover())

    def test_import_is_lazy_and_prototypes_preserve_pointer_width(self):
        self.assertIsNone(wh._default_library)
        wh.configure_library(self.library)
        self.assertIs(self.library.wirehair_encoder_create.restype, ctypes.c_void_p)
        self.assertIs(self.library.wirehair_decoder_create.restype, ctypes.c_void_p)
        self.assertEqual(ctypes.c_void_p,
                         self.library.wirehair_encoder_create_ex.argtypes[-1]._type_)
        self.assertEqual(
            [ctypes.c_void_p],
            self.library.wirehair_encoder_detach_input.argtypes)
        self.assertEqual("fake result 9", wh.result_string(9, self.library))
        self.assertEqual(16, ctypes.sizeof(wh.WirehairWireProfile))

    def test_canonical_and_compatibility_imports_share_public_objects(self):
        self.assertEqual("2.0.0", wh.__version__)
        self.assertEqual(wh.__version__, canonical.__version__)
        self.assertIs(wh.Encoder, canonical.Encoder)
        self.assertIs(wh.Decoder, canonical.Decoder)
        self.assertIs(wh.initialize, canonical.initialize)
        self.assertEqual(wh.__all__, canonical.__all__)
        self.assertNotIn("ctypes", wh.__all__)

    def test_explicit_wire_profiles(self):
        message = b"abcdefgh"
        with wh.Encoder.create(
                message, 4, self.library, owned=True,
                profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT) as encoder:
            with wh.Decoder.create(
                    len(message), 4, self.library,
                    profile_id=wh.WIREHAIR_LEGACY_PROFILE_CURRENT) as decoder:
                self.assertFalse(decoder.decode(0, encoder.encode(0)))
                self.assertTrue(decoder.decode(1, encoder.encode(1)))
                self.assertEqual(message, decoder.recover())

        with self.assertRaises(wh.WirehairError) as caught:
            wh.Encoder.create(
                message, 4, self.library, profile_id=0x123456789abcdef0)
        self.assertEqual("wirehair_wire_profile_init", caught.exception.operation)
        self.assertEqual(wh.Wirehair_InvalidInput, caught.exception.result)

    def test_windows_install_discovers_prefix_bin(self):
        prefix = os.path.abspath(os.path.join(os.path.dirname(wh.__file__), os.pardir))
        expected = os.path.join(prefix, "bin", "wirehair.dll")
        loaded = object()
        attempted = []

        def load(candidate):
            attempted.append(candidate)
            if candidate == expected:
                return loaded
            raise OSError("not this candidate")

        with mock.patch.object(sys, "platform", "win32"), \
                mock.patch.dict(
                    os.environ,
                    {"WIREHAIR_LIBRARY": "", "WIREHAIR_PREFIX": ""}), \
                mock.patch.object(wh.ctypes.util, "find_library", return_value=None), \
                mock.patch.object(wh.os.path, "isdir", return_value=False), \
                mock.patch.object(wh.ctypes, "CDLL", side_effect=load):
            self.assertIs(loaded, wh._load_wirehair())
        self.assertIn(expected, attempted)

    def test_platform_library_names_and_authoritative_overrides(self):
        self.assertEqual(
            ("wirehair.dll", "libwirehair.dll"),
            wh._library_names("win32"))
        self.assertEqual(
            ("libwirehair.dylib", "libwirehair.2.dylib"),
            wh._library_names("darwin"))
        self.assertEqual(
            ("libwirehair.so", "libwirehair.so.2"),
            wh._library_names("linux"))

        intended = os.path.abspath("missing-intended-wirehair.so")
        with mock.patch.dict(
                os.environ,
                {"WIREHAIR_LIBRARY": intended,
                 "WIREHAIR_PREFIX": "/must/not/fallback"},
                clear=True), \
                mock.patch.object(wh.ctypes, "CDLL", side_effect=OSError("missing")) \
                as load, \
                mock.patch.object(wh.ctypes.util, "find_library") as find, \
                self.assertRaisesRegex(OSError, "WIREHAIR_LIBRARY"):
            wh._load_wirehair()
        load.assert_called_once_with(intended)
        find.assert_not_called()

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.dict(
                    os.environ,
                    {"WIREHAIR_LIBRARY": "",
                     "WIREHAIR_PREFIX": temporary},
                    clear=True), \
                mock.patch.object(
                    wh.ctypes, "CDLL", side_effect=OSError("missing")), \
                mock.patch.object(wh.ctypes.util, "find_library") as find, \
                self.assertRaisesRegex(OSError, "WIREHAIR_PREFIX"):
            wh._load_wirehair()
        find.assert_not_called()

    def test_authoritative_prefix_finds_nested_platform_library_only(self):
        loaded = object()
        with tempfile.TemporaryDirectory() as temporary:
            prefix = os.path.abspath(temporary)
            nested = os.path.join(
                prefix, "lib64", "custom", "libwirehair.so.2")
            os.makedirs(os.path.dirname(nested))
            with open(nested, "wb"):
                pass
            attempted = []

            def load(candidate):
                attempted.append(candidate)
                if candidate == nested:
                    return loaded
                raise OSError("not this candidate")

            with mock.patch.object(sys, "platform", "linux"), \
                    mock.patch.dict(
                        os.environ,
                        {"WIREHAIR_LIBRARY": "",
                         "WIREHAIR_PREFIX": prefix},
                        clear=True), \
                    mock.patch.object(wh.ctypes, "CDLL", side_effect=load), \
                    mock.patch.object(wh.ctypes.util, "find_library") as find:
                self.assertIs(loaded, wh._load_wirehair())
            self.assertIn(nested, attempted)
            find.assert_not_called()

    def test_macos_authoritative_prefix_discovers_versioned_dylib(self):
        loaded = object()
        with tempfile.TemporaryDirectory() as temporary:
            prefix = os.path.abspath(temporary)
            expected = os.path.join(
                prefix, "lib", "libwirehair.2.dylib")
            os.makedirs(os.path.dirname(expected))
            with open(expected, "wb"):
                pass

            def load(candidate):
                if candidate == expected:
                    return loaded
                raise OSError("not this candidate")

            with mock.patch.object(sys, "platform", "darwin"), \
                    mock.patch.dict(
                        os.environ,
                        {"WIREHAIR_LIBRARY": "",
                         "WIREHAIR_PREFIX": prefix},
                        clear=True), \
                    mock.patch.object(wh.ctypes, "CDLL", side_effect=load), \
                    mock.patch.object(wh.os, "walk") as walk, \
                    mock.patch.object(wh.ctypes.util, "find_library") as find:
                self.assertIs(loaded, wh._load_wirehair())
            walk.assert_not_called()
            find.assert_not_called()

    def test_context_close_double_close_and_finalizer(self):
        encoder = wh.Encoder.create(b"abcdefgh", 4, self.library)
        pointer = encoder.pointer.value
        encoder.close()
        encoder.close()
        self.assertEqual([pointer], self.library.freed)
        with self.assertRaises(RuntimeError):
            encoder.encode(0)

        decoder = wh.Decoder.create(8, 4, self.library)
        pointer = decoder.pointer.value
        reference = weakref.ref(decoder)
        del decoder
        gc.collect()
        self.assertIsNone(reference())
        self.assertEqual(1, self.library.freed.count(pointer))

    def test_owned_and_borrowed_encoder_lifetimes(self):
        source = bytearray(b"abcdefgh")
        encoder = wh.Encoder.create(source, 4, self.library, owned=False)
        self.assertEqual(b"abcd", encoder.encode(0))
        self.assertTrue(hasattr(encoder, "_borrowed_source"))
        with self.assertRaises(BufferError):
            source.extend(b"!")
        encoder.close()
        self.assertFalse(hasattr(encoder, "_borrowed_source"))
        source.extend(b"!")
        self.assertEqual(b"abcdefgh!", bytes(source))

        source = bytearray(b"abcdefgh")
        encoder = wh.Encoder.create(source, 4, self.library, owned=True)
        source[:] = b"XXXXXXXX"
        self.assertEqual(b"abcd", encoder.encode(0))
        encoder.close()

    def test_encoder_detach_input_releases_borrow_and_is_idempotent(self):
        source = bytearray(b"abcdefghij")
        encoder = wh.Encoder.create(source, 4, self.library, owned=False)
        self.assertTrue(hasattr(encoder, "_borrowed_source"))
        self.assertIsNone(encoder.detach_input())
        self.assertFalse(hasattr(encoder, "_borrowed_source"))
        source[:] = b"XXXXXXXXXX"
        self.assertEqual(b"abcd", encoder.encode(0))
        self.assertEqual(b"ij", encoder.encode(2))
        self.assertIsNone(encoder.detach_input())
        self.assertEqual(2, self.library.detach_calls)

        self.library.detach_result = wh.Wirehair_InvalidInput
        with self.assertRaises(wh.WirehairError):
            encoder.detach_input()
        encoder.close()
        with self.assertRaises(RuntimeError):
            encoder.detach_input()

    def test_encoder_detach_failure_retains_borrowed_source_pin(self):
        source = bytearray(b"abcdefghij")
        encoder = wh.Encoder.create(source, 4, self.library, owned=False)
        self.library.detach_result = wh.Wirehair_InvalidInput
        with self.assertRaises(wh.WirehairError):
            encoder.detach_input()
        self.assertTrue(hasattr(encoder, "_borrowed_source"))
        with self.assertRaises(BufferError):
            source.extend(b"!")

        self.library.detach_result = wh.Wirehair_Success
        self.assertIsNone(encoder.detach_input())
        source.extend(b"!")
        self.assertEqual(b"abcdefghij!", bytes(source))
        encoder.close()

    def test_reusable_output_buffers_exact_counts_and_tail_preservation(self):
        message = b"abcdefghij"
        with wh.Encoder.create(message, 4, self.library) as encoder:
            packet = bytearray(b"\xa5" * 9)
            self.assertEqual(4, encoder.encode_into(0, memoryview(packet)[2:]))
            self.assertEqual(b"\xa5\xa5abcd\xa5\xa5\xa5", bytes(packet))

            final = bytearray(2)
            self.assertEqual(2, encoder.encode_into(2, final))
            self.assertEqual(b"ij", bytes(final))
            final.extend(b"!")
            self.assertEqual(b"ij!", bytes(final))

            repair = bytearray(b"\xa5" * 8)
            self.assertEqual(4, encoder.encode_into(3, memoryview(repair)[1:5]))
            self.assertEqual(b"\xa5", bytes(repair[:1]))
            self.assertEqual(b"\xa5\xa5\xa5", bytes(repair[5:]))

            with wh.Decoder.create(len(message), 4, self.library) as decoder:
                self.assertFalse(decoder.decode(0, encoder.encode(0)))
                self.assertFalse(decoder.decode(1, encoder.encode(1)))
                self.assertTrue(decoder.decode(2, encoder.encode(2)))

                recovered = bytearray(b"\xa5" * (len(message) + 5))
                self.assertEqual(len(message), decoder.recover_into(recovered))
                self.assertEqual(message, bytes(recovered[:len(message)]))
                self.assertEqual(b"\xa5" * 5, bytes(recovered[len(message):]))
                recovered.extend(b"!")

                block = bytearray(b"\xa5" * 6)
                self.assertEqual(
                    4, decoder.recover_block_into(1, memoryview(block)[1:5]))
                self.assertEqual(b"\xa5efgh\xa5", bytes(block))

                final_block = bytearray(2)
                self.assertEqual(
                    2, decoder.recover_block_into(2, final_block))
                self.assertEqual(b"ij", bytes(final_block))

                self.assertEqual(message, decoder.recover())
                self.assertEqual(b"ij", decoder.recover_block(2, capacity=2))

    def test_output_buffer_rejections_happen_before_native_calls(self):
        message = b"abcdefghij"
        encoder = wh.Encoder.create(message, 4, self.library)
        decoder = wh.Decoder.create(len(message), 4, self.library)

        def assert_encode_rejected(output, exception):
            calls = self.library.encode_calls
            with self.assertRaises(exception):
                encoder.encode_into(0, output)
            self.assertEqual(calls, self.library.encode_calls)

        readonly = memoryview(b"readonly")
        noncontiguous = memoryview(bytearray(8))[::2]
        released = memoryview(bytearray(4))
        released.release()
        assert_encode_rejected(readonly, TypeError)
        assert_encode_rejected(noncontiguous, ValueError)
        assert_encode_rejected(released, ValueError)
        undersized = bytearray(b"xyz")
        assert_encode_rejected(undersized, ValueError)
        self.assertEqual(b"xyz", bytes(undersized))
        assert_encode_rejected(None, TypeError)

        class HugeView:
            readonly = False
            c_contiguous = True

            def __init__(self, size):
                self.nbytes = size

            def cast(self, _format):
                return self

        with mock.patch.object(
                wh, "memoryview", return_value=HugeView(0x100000000),
                create=True):
            assert_encode_rejected(object(), OverflowError)

        recover_calls = self.library.recover_calls
        with self.assertRaises(ValueError):
            decoder.recover_into(bytearray(len(message) - 1))
        self.assertEqual(recover_calls, self.library.recover_calls)
        with mock.patch.object(
                wh, "memoryview", return_value=HugeView(0x10000000000000000),
                create=True):
            with self.assertRaises(OverflowError):
                decoder.recover_into(object())
        self.assertEqual(recover_calls, self.library.recover_calls)

        block_calls = self.library.recover_block_calls
        with self.assertRaises(ValueError):
            decoder.recover_block_into(0, bytearray(3))
        with self.assertRaises(ValueError):
            decoder.recover_block_into(3, bytearray(4))
        self.assertEqual(block_calls, self.library.recover_block_calls)

        sentinel = bytearray(b"KEEP")
        encode_calls = self.library.encode_calls
        encoder.close()
        with self.assertRaisesRegex(RuntimeError, "closed"):
            encoder.encode_into(0, sentinel)
        self.assertEqual(b"KEEP", bytes(sentinel))
        self.assertEqual(encode_calls, self.library.encode_calls)

        decoder.close()
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover_into(bytearray(len(message)))
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover_block_into(0, sentinel)

    def test_output_buffer_is_pinned_only_for_native_call(self):
        entered = threading.Event()
        release = threading.Event()
        output = bytearray(4)
        errors = []

        def block_native_call():
            entered.set()
            if not release.wait(10):
                raise RuntimeError("test did not release fake native call")

        self.library.before_encode = block_native_call
        with wh.Encoder.create(b"abcdefgh", 4, self.library) as encoder:
            def run_encode():
                try:
                    encoder.encode_into(0, output)
                except BaseException as exc:
                    errors.append(exc)

            worker = threading.Thread(target=run_encode)
            worker.start()
            try:
                self.assertTrue(entered.wait(10))
                with self.assertRaises(BufferError):
                    output.extend(b"!")
            finally:
                release.set()
                worker.join(10)
            self.assertFalse(worker.is_alive())
            self.assertEqual([], errors)
            self.assertEqual(b"abcd", bytes(output))
            output.extend(b"!")
            self.assertEqual(b"abcd!", bytes(output))

    def test_output_buffer_pin_is_released_after_call_exception(self):
        output = bytearray(4)
        with wh.Encoder.create(b"abcdefgh", 4, self.library) as encoder:
            self.library.raise_encode = True
            with self.assertRaisesRegex(RuntimeError, "injected encode"):
                encoder.encode_into(0, output)
            output.extend(b"!")
            self.assertEqual(5, len(output))

        self.library = FakeLibrary()
        decoder = self.completed_decoder()
        recovered = bytearray(8)
        self.library.raise_recover = True
        with self.assertRaisesRegex(RuntimeError, "injected recover"):
            decoder.recover_into(recovered)
        recovered.extend(b"!")
        self.assertEqual(9, len(recovered))
        decoder.close()

    def test_decode_recover_block_and_conversion_transfer_ownership(self):
        message = b"abcdefgh"
        with wh.Encoder.create(message, 4, self.library) as source:
            decoder = wh.Decoder.create(len(message), 4, self.library)
            self.assertFalse(decoder.decode(0, source.encode(0)))
            self.assertTrue(decoder.decode(1, source.encode(1)))
            self.assertEqual(message, decoder.recover())
            self.assertEqual(b"efgh", decoder.recover_block(1, capacity=4))
            self.assertEqual(4, self.library.last_recover_capacity)
            calls = self.library.recover_block_calls
            with self.assertRaises(ValueError):
                decoder.recover_block(0, capacity=3)
            self.assertEqual(calls, self.library.recover_block_calls)
            converted_pointer = decoder.pointer.value
            converted = decoder.become_encoder()
            self.assertTrue(decoder.closed)
            self.assertEqual(b"abcd", converted.encode(0))
            converted_output = bytearray(4)
            self.assertEqual(4, converted.encode_into(1, converted_output))
            self.assertEqual(b"efgh", bytes(converted_output))
            self.assertIsNone(converted.detach_input())
            self.assertEqual(b"abcd", converted.encode(0))
            decoder.close()
            converted.close()
            self.assertEqual(1, self.library.freed.count(converted_pointer))

    def test_conversion_wrapper_failures_have_one_cleanup_owner(self):
        original_init = wh.Encoder.__init__
        for failure_point in (
                "before", "after", "close", "detach_finalizer",
                "replace_finalizer"):
            with self.subTest(failure_point=failure_point):
                self.library = FakeLibrary()
                decoder = self.completed_decoder()
                pointer = decoder.pointer.value

                def fail(instance, *arguments, **keywords):
                    if failure_point == "before":
                        raise RuntimeError("injected wrapper failure")
                    original_init(instance, *arguments, **keywords)
                    if failure_point == "close":
                        instance.close()
                    elif failure_point == "detach_finalizer":
                        replacement = weakref.finalize(instance, lambda: None)
                        replacement.detach()
                        instance._finalizer = replacement
                    elif failure_point == "replace_finalizer":
                        instance._finalizer = object()
                    raise RuntimeError("injected wrapper failure")

                with mock.patch.object(wh.Encoder, "__init__", new=fail):
                    with self.assertRaisesRegex(
                            RuntimeError, "injected wrapper failure"):
                        decoder.become_encoder()
                gc.collect()
                if failure_point == "close":
                    self.assertTrue(decoder.closed)
                    self.assertEqual([pointer], self.library.freed)
                else:
                    self.assert_decoder_usable(decoder)
                    self.assertEqual([], self.library.freed)
                self.assertEqual(0, self.library.become_encoder_calls)
                decoder.close()
                gc.collect()
                self.assertEqual(1, self.library.freed.count(pointer))

    def test_conversion_finalizer_failures_happen_before_native_mutation(self):
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        with mock.patch.object(
                wh.weakref, "finalize",
                side_effect=MemoryError("injected finalizer allocation failure")):
            with self.assertRaisesRegex(MemoryError, "finalizer allocation"):
                decoder.become_encoder()
        gc.collect()
        self.assert_decoder_usable(decoder)
        self.assertEqual(0, self.library.become_encoder_calls)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

        self.library = FakeLibrary()
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        original_finish = wh._Codec._finish_construction

        def fail_after_finalizer(instance):
            original_finish(instance)
            raise MemoryError("injected failure after finalizer allocation")

        with mock.patch.object(
                wh._Codec, "_finish_construction", new=fail_after_finalizer):
            with self.assertRaisesRegex(MemoryError, "after finalizer"):
                decoder.become_encoder()
        gc.collect()
        self.assert_decoder_usable(decoder)
        self.assertEqual(0, self.library.become_encoder_calls)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

    def test_conversion_failed_finalizer_detach_closes_shared_owner(self):
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value

        class UndetachableFinalizer:
            alive = True

            def detach(self):
                raise MemoryError("injected finalizer detach failure")

        def fail_with_undetachable_finalizer(instance):
            finalizer = UndetachableFinalizer()
            instance._owner_finalizer = finalizer
            instance._finalizer = finalizer
            raise MemoryError("injected wrapper finalization failure")

        with mock.patch.object(
                wh._Codec, "_finish_construction",
                new=fail_with_undetachable_finalizer):
            with self.assertRaisesRegex(MemoryError, "wrapper finalization"):
                decoder.become_encoder()
        self.assertEqual(0, self.library.become_encoder_calls)
        self.assertTrue(decoder.closed)
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover()
        self.assertEqual([pointer], self.library.freed)
        decoder.close()
        gc.collect()
        self.assertEqual(1, self.library.freed.count(pointer))

    def test_conversion_native_failure_leaves_truthful_source_state(self):
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        self.library.become_encoder_result = wh.Wirehair_InvalidInput
        with self.assertRaises(wh.WirehairError) as caught:
            decoder.become_encoder()
        self.assertEqual(wh.Wirehair_InvalidInput, caught.exception.result)
        self.assert_decoder_usable(decoder)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

        self.library = FakeLibrary()
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        self.library.become_encoder_result = wh.Wirehair_Error
        with self.assertRaises(wh.WirehairError) as caught:
            decoder.become_encoder()
        self.assertEqual(wh.Wirehair_Error, caught.exception.result)
        self.assertTrue(decoder.closed)
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover()
        self.assertEqual([pointer], self.library.freed)
        decoder.close()
        self.assertEqual(1, self.library.freed.count(pointer))

    def test_post_native_conversion_exceptions_close_exactly_once(self):
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        self.library.raise_after_become_encoder = True
        with self.assertRaisesRegex(MemoryError, "post-conversion"):
            decoder.become_encoder()
        self.assertTrue(decoder.closed)
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover()
        self.assertEqual([pointer], self.library.freed)
        decoder.close()
        gc.collect()
        self.assertEqual(1, self.library.freed.count(pointer))

        self.library = FakeLibrary()
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value

        def fail_commit():
            raise MemoryError("injected ownership commit failure")

        decoder._relinquish = fail_commit
        with self.assertRaisesRegex(MemoryError, "ownership commit"):
            decoder.become_encoder()
        self.assertTrue(decoder.closed)
        with self.assertRaisesRegex(RuntimeError, "closed"):
            decoder.recover()
        self.assertEqual([pointer], self.library.freed)
        decoder.close()
        gc.collect()
        self.assertEqual(1, self.library.freed.count(pointer))

    def test_wrapper_new_failures_preserve_exact_ownership(self):
        for wrapper, create in (
                (wh.Encoder,
                 lambda: wh.Encoder.create(b"abcdefgh", 4, self.library)),
                (wh.Decoder,
                 lambda: wh.Decoder.create(8, 4, self.library))):
            with self.subTest(wrapper=wrapper.__name__):
                self.library = FakeLibrary()
                pointer = self.library.next_pointer
                with mock.patch.object(
                        wrapper, "__new__",
                        side_effect=RuntimeError("injected new failure")):
                    with self.assertRaisesRegex(
                            RuntimeError, "injected new failure"):
                        create()
                self.assertEqual([pointer], self.library.freed)

        self.library = FakeLibrary()
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        with mock.patch.object(
                wh.Encoder, "__new__",
                side_effect=RuntimeError("injected new failure")):
            with self.assertRaisesRegex(RuntimeError, "injected new failure"):
                decoder.become_encoder()
        self.assert_decoder_usable(decoder)
        self.assertEqual(0, self.library.become_encoder_calls)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

    def test_converted_encoder_gc_frees_shared_owner_once(self):
        decoder = self.completed_decoder()
        pointer = decoder.pointer.value
        converted = decoder.become_encoder()
        self.assertTrue(decoder.closed)
        reference = weakref.ref(converted)
        del converted
        gc.collect()
        self.assertIsNone(reference())
        self.assertEqual(1, self.library.freed.count(pointer))
        decoder.close()
        self.assertEqual(1, self.library.freed.count(pointer))

    def test_creation_wrapper_failure_releases_native_handle(self):
        for wrapper, create in (
                (wh.Encoder,
                 lambda: wh.Encoder.create(b"abcdefgh", 4, self.library)),
                (wh.Decoder,
                 lambda: wh.Decoder.create(8, 4, self.library))):
            original_init = wrapper.__init__
            for failure_point in (
                    "before", "after", "close", "detach_finalizer",
                    "replace_finalizer"):
                with self.subTest(
                        wrapper=wrapper.__name__, failure_point=failure_point):
                    self.library = FakeLibrary()
                    pointer = self.library.next_pointer

                    def fail(instance, *arguments, **keywords):
                        if failure_point == "before":
                            raise RuntimeError("injected wrapper failure")
                        original_init(instance, *arguments, **keywords)
                        if failure_point == "close":
                            instance.close()
                        elif failure_point == "detach_finalizer":
                            replacement = weakref.finalize(instance, lambda: None)
                            replacement.detach()
                            instance._finalizer = replacement
                        elif failure_point == "replace_finalizer":
                            instance._finalizer = object()
                        raise RuntimeError("injected wrapper failure")

                    with mock.patch.object(wrapper, "__init__", new=fail):
                        with self.assertRaisesRegex(
                                RuntimeError, "injected wrapper failure"):
                            create()
                    gc.collect()
                    self.assertEqual([pointer], self.library.freed)
                    self.assertNotIn(pointer, self.library.codecs)

    def test_validation_and_structured_creation_errors(self):
        self.library.encoder_result = wh.Wirehair_OOM
        with self.assertRaises(wh.WirehairError) as caught:
            wh.Encoder.create(b"abcdefgh", 4, self.library)
        self.assertEqual(wh.Wirehair_OOM, caught.exception.result)
        self.assertEqual([], self.library.freed)
        self.library.encoder_result = wh.Wirehair_Success
        with self.assertRaises((TypeError, ValueError)):
            wh.Encoder.create("not bytes", 4, self.library)
        with self.assertRaises(ValueError):
            wh.Decoder.create(8, 0, self.library)
        with self.assertRaises(ValueError):
            wh.Decoder.create(8, 0x80000000, self.library)

    def run_example(self, **options):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            result = wh._run_example(library=self.library, **options)
        return result, stdout.getvalue(), stderr.getvalue()

    def test_example_non_ascii_success_and_exact_cleanup(self):
        message = "héllø 世界".encode("utf-8") * 4
        result, output, error = self.run_example(
            message=message, block_size=8, max_packets=80)
        self.assertEqual((0, ""), (result, error))
        self.assertIn("recovery succeeded", output)
        self.assertEqual(2, len(self.library.freed))
        self.assertEqual(2, len(set(self.library.freed)))

    def test_example_perpetual_need_more_is_bounded(self):
        self.library.always_need_more = True
        result, _, error = self.run_example(
            message=b"abcdefgh", block_size=4, max_packets=7)
        self.assertEqual(1, result)
        self.assertIn("still needs data", error)
        self.assertEqual(2, len(self.library.freed))

    def test_example_terminal_and_corrupt_recovery_fail(self):
        self.library.decode_terminal = wh.Wirehair_ExtraInsufficient
        result, _, error = self.run_example(message=b"abcdefgh", block_size=4)
        self.assertEqual(1, result)
        self.assertIn("wirehair_decode failed", error)
        self.assertEqual(2, len(self.library.freed))

        self.library = FakeLibrary()
        self.library.corrupt_recovery = True
        result, _, error = self.run_example(message=b"abcdefgh", block_size=4)
        self.assertEqual(1, result)
        self.assertIn("digest verification", error)
        self.assertEqual(2, len(self.library.freed))

    def test_example_exceptions_and_decoder_creation_failure_cleanup(self):
        for attribute, expected_frees in (("raise_encode", 2), ("raise_recover", 2)):
            with self.subTest(attribute=attribute):
                self.library = FakeLibrary()
                setattr(self.library, attribute, True)
                result, _, error = self.run_example(message=b"abcdefgh", block_size=4)
                self.assertEqual(1, result)
                self.assertIn("injected", error)
                self.assertEqual(expected_frees, len(self.library.freed))

        self.library = FakeLibrary()
        self.library.decoder_result = wh.Wirehair_OOM
        result, _, error = self.run_example(message=b"abcdefgh", block_size=4)
        self.assertEqual(1, result)
        self.assertIn("decoder_create", error)
        self.assertEqual(1, len(self.library.freed))


if __name__ == "__main__":
    unittest.main()
