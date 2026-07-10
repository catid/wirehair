#!/usr/bin/env python3
"""Unit tests for the ctypes Wirehair binding."""

import ctypes
import gc
import io
import os
import sys
import unittest
import weakref
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))
import whirehair as wh


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
        self.decode_terminal = None
        self.always_need_more = False
        self.recover_result = wh.Wirehair_Success
        self.become_encoder_result = wh.Wirehair_Success
        self.corrupt_recovery = False
        self.raise_encode = False
        self.raise_recover = False
        self.last_recover_capacity = None

        for name in (
            "wirehair_result_string", "wirehair_init_", "wirehair_encoder_create",
            "wirehair_encoder_create_ex", "wirehair_encoder_create_owned",
            "wirehair_encoder_create_owned_ex", "wirehair_encode",
            "wirehair_decoder_create", "wirehair_decoder_create_ex",
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

    def _encoder_message(self, state):
        if "message" in state:
            return state["message"]
        return ctypes.string_at(state["source"], state["size"])

    def _wirehair_encode(self, pointer, block_id, output, capacity, written):
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
        state = self._state(pointer)
        return self._recover_block(state, block_id, output, state["block"], written)

    def _wirehair_recover_block_ex(self, pointer, block_id, output, capacity, written):
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
        if self.become_encoder_result != wh.Wirehair_Success:
            return self.become_encoder_result
        state = self._state(pointer)
        state["message"] = self._original_message(state)
        state["kind"] = "encoder"
        return wh.Wirehair_Success

    def _wirehair_free(self, pointer):
        value = number(pointer)
        self.freed.append(value)
        self.codecs.pop(value, None)


class BindingTests(unittest.TestCase):
    def setUp(self):
        self.library = FakeLibrary()

    def test_import_is_lazy_and_prototypes_preserve_pointer_width(self):
        self.assertIsNone(wh._default_library)
        wh.configure_library(self.library)
        self.assertIs(self.library.wirehair_encoder_create.restype, ctypes.c_void_p)
        self.assertIs(self.library.wirehair_decoder_create.restype, ctypes.c_void_p)
        self.assertEqual(ctypes.c_void_p,
                         self.library.wirehair_encoder_create_ex.argtypes[-1]._type_)
        self.assertEqual("fake result 9", wh.result_string(9, self.library))

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
                mock.patch.dict(os.environ, {"WIREHAIR_LIBRARY": ""}), \
                mock.patch.object(wh.ctypes.util, "find_library", return_value=None), \
                mock.patch.object(wh.os.path, "isdir", return_value=False), \
                mock.patch.object(wh.ctypes, "CDLL", side_effect=load):
            self.assertIs(loaded, wh._load_wirehair())
        self.assertIn(expected, attempted)

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
        del source
        gc.collect()
        self.assertEqual(b"abcd", encoder.encode(0))
        self.assertTrue(hasattr(encoder, "_borrowed_source"))
        encoder.close()

        source = bytearray(b"abcdefgh")
        encoder = wh.Encoder.create(source, 4, self.library, owned=True)
        source[:] = b"XXXXXXXX"
        self.assertEqual(b"abcd", encoder.encode(0))
        encoder.close()

    def test_decode_recover_block_and_conversion_transfer_ownership(self):
        message = b"abcdefgh"
        with wh.Encoder.create(message, 4, self.library) as source:
            decoder = wh.Decoder.create(len(message), 4, self.library)
            self.assertFalse(decoder.decode(0, source.encode(0)))
            self.assertTrue(decoder.decode(1, source.encode(1)))
            self.assertEqual(message, decoder.recover())
            self.assertEqual(b"efgh", decoder.recover_block(1, capacity=4))
            self.assertEqual(4, self.library.last_recover_capacity)
            with self.assertRaises(wh.WirehairError):
                decoder.recover_block(0, capacity=3)
            converted_pointer = decoder.pointer.value
            converted = decoder.become_encoder()
            self.assertTrue(decoder.closed)
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
                decoder = wh.Decoder.create(8, 4, self.library)
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
                    self.assertFalse(decoder.closed)
                    self.assertEqual([], self.library.freed)
                decoder.close()
                gc.collect()
                self.assertEqual(1, self.library.freed.count(pointer))

    def test_conversion_native_failure_retains_cleanup_ownership(self):
        decoder = wh.Decoder.create(8, 4, self.library)
        pointer = decoder.pointer.value
        self.library.become_encoder_result = wh.Wirehair_Error
        with self.assertRaises(wh.WirehairError) as caught:
            decoder.become_encoder()
        self.assertEqual(wh.Wirehair_Error, caught.exception.result)
        self.assertFalse(decoder.closed)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

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
        decoder = wh.Decoder.create(8, 4, self.library)
        pointer = decoder.pointer.value
        with mock.patch.object(
                wh.Encoder, "__new__",
                side_effect=RuntimeError("injected new failure")):
            with self.assertRaisesRegex(RuntimeError, "injected new failure"):
                decoder.become_encoder()
        self.assertFalse(decoder.closed)
        self.assertEqual([], self.library.freed)
        decoder.close()
        self.assertEqual([pointer], self.library.freed)

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
        self.assertIn("differs", error)
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
