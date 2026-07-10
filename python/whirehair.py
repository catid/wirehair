# Copyleft (c) 2019 Daniel Norte de Moraes <danielcheagle@gmail.com>.
# This code is placed in the public domain and provided without warranty.

"""Safe ctypes bindings and a bounded Wirehair round-trip example."""

import ctypes
import ctypes.util
import os
import sys
import weakref


WIREHAIR_VERSION = 2
WIREHAIR_WIRE_PROFILE_VERSION = 1
WIREHAIR_LEGACY_PROFILE_PRE_FIXUP = 0xe1b9f77f1c90f680
WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07 = 0x4d241359db07bb07
WIREHAIR_LEGACY_PROFILE_CURRENT = WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07
WIREHAIR_ENCODER_OWN_INPUT = 1

Wirehair_Success = 0
Wirehair_NeedMore = 1
Wirehair_InvalidInput = 2
Wirehair_BadDenseSeed = 3
Wirehair_BadPeelSeed = 4
Wirehair_BadInput_SmallN = 5
Wirehair_BadInput_LargeN = 6
Wirehair_ExtraInsufficient = 7
Wirehair_Error = 8
Wirehair_OOM = 9
Wirehair_UnsupportedPlatform = 10
WirehairResult_Count = 11
WirehairResult_Padding = 0x7fffffff

_default_library = None


class WirehairWireProfile(ctypes.Structure):
    """In-process descriptor for a trusted legacy equation profile."""

    _fields_ = [
        ("struct_bytes", ctypes.c_uint32),
        ("profile_version", ctypes.c_uint32),
        ("profile_id", ctypes.c_uint64),
    ]


class WirehairError(RuntimeError):
    """A Wirehair C API operation returned a terminal result."""

    def __init__(self, operation, result, library=None):
        self.operation = operation
        self.result = int(result)
        try:
            description = result_string(self.result, library)
        except Exception:
            description = "result %d" % self.result
        super().__init__("%s failed: %s" % (operation, description))


def _load_wirehair():
    if sys.platform.startswith("win"):
        names = ["wirehair.dll", "libwirehair.dll"]
    elif sys.platform == "darwin":
        names = ["libwirehair.dylib", "libwirehair.2.dylib"]
    else:
        names = ["libwirehair.so", "libwirehair.so.2"]

    candidates = []
    env_path = os.environ.get("WIREHAIR_LIBRARY")
    if env_path:
        candidates.append(env_path)

    here = os.path.abspath(os.path.dirname(__file__))
    prefix = os.path.abspath(os.path.join(here, os.pardir))
    roots = (os.path.join(prefix, "bin"), os.path.join(prefix, "lib"),
             os.path.join(prefix, "lib64"))
    for directory in (here,) + roots:
        for name in names:
            candidates.append(os.path.join(directory, name))
    for root in roots:
        if os.path.isdir(root):
            for directory, _, files in os.walk(root):
                for name in names:
                    if name in files:
                        candidates.append(os.path.join(directory, name))

    found = ctypes.util.find_library("wirehair")
    if found:
        candidates.append(found)
    candidates.extend(names)

    errors = []
    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            return ctypes.CDLL(candidate)
        except OSError as exc:
            errors.append("%s: %s" % (candidate, exc))
    raise OSError("Unable to load Wirehair shared library. Tried:\n" +
                  "\n".join(errors))


def _prototype(library, name, argtypes, restype):
    function = getattr(library, name)
    function.argtypes = argtypes
    function.restype = restype


def configure_library(library):
    """Declare every public C API prototype on an injected CDLL-like object."""
    if getattr(library, "_wirehair_configured", False):
        return library
    codec = ctypes.c_void_p
    result = ctypes.c_int32
    u32p = ctypes.POINTER(ctypes.c_uint32)
    codecpp = ctypes.POINTER(codec)
    profilep = ctypes.POINTER(WirehairWireProfile)
    _prototype(library, "wirehair_result_string", [result], ctypes.c_char_p)
    _prototype(library, "wirehair_init_", [ctypes.c_int32], result)
    _prototype(library, "wirehair_wire_profile_init",
               [ctypes.c_uint64, profilep], result)
    _prototype(library, "wirehair_encoder_create",
               [codec, ctypes.c_void_p, ctypes.c_uint64, ctypes.c_uint32], codec)
    _prototype(library, "wirehair_encoder_create_ex",
               [codec, ctypes.c_void_p, ctypes.c_uint64, ctypes.c_uint32, codecpp], result)
    _prototype(library, "wirehair_encoder_create_owned",
               [codec, ctypes.c_void_p, ctypes.c_uint64, ctypes.c_uint32], codec)
    _prototype(library, "wirehair_encoder_create_owned_ex",
               [codec, ctypes.c_void_p, ctypes.c_uint64, ctypes.c_uint32, codecpp], result)
    _prototype(library, "wirehair_encoder_create_profile_ex",
               [codec, ctypes.c_void_p, ctypes.c_uint64, ctypes.c_uint32,
                profilep, ctypes.c_uint32, codecpp], result)
    _prototype(library, "wirehair_encode",
               [codec, ctypes.c_uint32, ctypes.c_void_p, ctypes.c_uint32, u32p], result)
    _prototype(library, "wirehair_decoder_create",
               [codec, ctypes.c_uint64, ctypes.c_uint32], codec)
    _prototype(library, "wirehair_decoder_create_ex",
               [codec, ctypes.c_uint64, ctypes.c_uint32, codecpp], result)
    _prototype(library, "wirehair_decoder_create_profile_ex",
               [codec, ctypes.c_uint64, ctypes.c_uint32, profilep, codecpp], result)
    _prototype(library, "wirehair_decode",
               [codec, ctypes.c_uint32, ctypes.c_void_p, ctypes.c_uint32], result)
    _prototype(library, "wirehair_recover",
               [codec, ctypes.c_void_p, ctypes.c_uint64], result)
    _prototype(library, "wirehair_recover_block",
               [codec, ctypes.c_uint32, ctypes.c_void_p, u32p], result)
    _prototype(library, "wirehair_recover_block_ex",
               [codec, ctypes.c_uint32, ctypes.c_void_p, ctypes.c_uint32, u32p], result)
    _prototype(library, "wirehair_decoder_becomes_encoder", [codec], result)
    _prototype(library, "wirehair_free", [codec], None)
    library._wirehair_configured = True
    return library


def get_library(library=None):
    """Return a configured injected or lazily loaded shared library."""
    global _default_library
    if library is not None:
        return configure_library(library)
    if _default_library is None:
        _default_library = configure_library(_load_wirehair())
    return _default_library


def result_string(result, library=None):
    value = get_library(library).wirehair_result_string(int(result))
    if not value:
        return "Wirehair result %d" % int(result)
    if isinstance(value, str):
        return value
    return value.decode("utf-8", "replace")


def initialize(library=None, expected_version=WIREHAIR_VERSION):
    library = get_library(library)
    result = library.wirehair_init_(int(expected_version))
    if result != Wirehair_Success:
        raise WirehairError("wirehair_init", result, library)
    return library


def _bounded_integer(value, name, maximum):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("%s must be an integer" % name)
    if value < 0 or value > maximum:
        raise ValueError("%s is outside the supported range" % name)
    return value


def _input_buffer(data):
    try:
        view = memoryview(data)
    except TypeError as exc:
        raise TypeError("data must be bytes-like") from exc
    if not view.contiguous:
        raise ValueError("data must be contiguous")
    view = view.cast("B")
    if not view:
        raise ValueError("data must not be empty")
    array_type = ctypes.c_uint8 * len(view)
    if view.readonly:
        array = array_type.from_buffer_copy(view)
    else:
        array = array_type.from_buffer(view)
    return view, array


def _check(operation, result, library):
    if result != Wirehair_Success:
        raise WirehairError(operation, result, library)


def _profile_descriptor(profile_id, library):
    profile_id = _bounded_integer(profile_id, "profile_id", 0xffffffffffffffff)
    profile = WirehairWireProfile()
    result = library.wirehair_wire_profile_init(
        profile_id, ctypes.byref(profile))
    _check("wirehair_wire_profile_init", result, library)
    return profile


def _release_codec(library, pointer_value):
    library.wirehair_free(ctypes.c_void_p(pointer_value))


class _NativeCodecOwner:
    """Idempotent owner shared while a native handle changes wrappers."""

    def __init__(self, library, pointer):
        value = (pointer.value if isinstance(pointer, ctypes.c_void_p)
                 else int(pointer))
        if not value:
            raise RuntimeError("Wirehair returned a null codec after success")
        self.library = library
        self.pointer_value = value

    @property
    def alive(self):
        return self.pointer_value is not None

    def release(self):
        value = self.pointer_value
        if value is None:
            return
        self.pointer_value = None
        _release_codec(self.library, value)


def _construct_owned_codec(wrapper_type, library, pointer, message_bytes,
                           block_bytes, owner=None):
    release_on_failure = owner is None
    if owner is None:
        owner = _NativeCodecOwner(library, pointer)
    elif (owner.library is not library or
          owner.pointer_value != pointer.value):
        raise RuntimeError("native codec owner does not match the handle")
    try:
        wrapper = wrapper_type.__new__(wrapper_type)
    except BaseException:
        if release_on_failure:
            owner.release()
        raise
    try:
        wrapper_type.__init__(
            wrapper, library, owner, message_bytes, block_bytes)
        _Codec._finish_construction(wrapper)
    except BaseException:
        if release_on_failure:
            owner.release()
        else:
            try:
                _Codec._relinquish(wrapper)
            except BaseException:
                # A prepared transfer wrapper must never retain a finalizer
                # that can later free the decoder's still-shared handle.  If
                # detaching that finalizer itself fails, close the shared
                # owner so neither wrapper can advertise a usable codec.
                owner.release()
        raise
    return wrapper


class _Codec:
    def __init__(self, library, owner, message_bytes, block_bytes):
        deferred = isinstance(owner, _NativeCodecOwner)
        if not deferred:
            owner = _NativeCodecOwner(library, owner)
        if owner.library is not library or not owner.alive:
            raise RuntimeError("native codec owner is not available")
        self._library = library
        self._owner = owner
        self._owns_codec = True
        self._pointer_value = owner.pointer_value
        self.message_bytes = message_bytes
        self.block_bytes = block_bytes
        self._owner_finalizer = None
        self._finalizer = None
        if not deferred:
            _Codec._finish_construction(self)

    def _finish_construction(self):
        if (not getattr(self, "_owns_codec", False) or
                not self._owner.alive):
            raise RuntimeError("codec was closed during wrapper construction")
        finalizer = weakref.finalize(self, self._owner.release)
        self._owner_finalizer = finalizer
        self._finalizer = finalizer

    def _relinquish(self):
        self._owns_codec = False
        self._pointer_value = None
        finalizer = getattr(self, "_owner_finalizer", None)
        if finalizer is not None and finalizer.alive:
            finalizer.detach()

    @property
    def closed(self):
        return (not getattr(self, "_owns_codec", False) or
                not self._owner.alive)

    @property
    def pointer(self):
        if self.closed:
            raise RuntimeError("codec is closed")
        return ctypes.c_void_p(self._owner.pointer_value)

    def close(self):
        if getattr(self, "_owns_codec", False):
            self._owns_codec = False
            self._pointer_value = None
            try:
                self._owner.release()
            finally:
                finalizer = self._owner_finalizer
                if finalizer is not None and finalizer.alive:
                    finalizer.detach()

    def __enter__(self):
        if self.closed:
            raise RuntimeError("codec is closed")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False


class Encoder(_Codec):
    """Context-managed Wirehair encoder."""

    @classmethod
    def create(cls, message, block_bytes, library=None, owned=True,
               profile_id=None):
        """Create an encoder.

        ``owned=True`` asks the C library to copy the message.  With
        ``owned=False`` the binding retains a writable view, or a stable
        snapshot for read-only inputs, for the complete encoder lifetime.
        ``profile_id`` selects a trusted legacy equation profile; ``None``
        uses the raw API's frozen current profile.
        """
        library = initialize(library)
        block_bytes = _bounded_integer(block_bytes, "block_bytes", 0x7fffffff)
        if block_bytes == 0:
            raise ValueError("block_bytes must be positive")
        view, buffer = _input_buffer(message)
        output = ctypes.c_void_p()
        if profile_id is None:
            create = (library.wirehair_encoder_create_owned_ex if owned else
                      library.wirehair_encoder_create_ex)
            result = create(None, ctypes.cast(buffer, ctypes.c_void_p),
                            len(view), block_bytes, ctypes.byref(output))
        else:
            profile = _profile_descriptor(profile_id, library)
            flags = WIREHAIR_ENCODER_OWN_INPUT if owned else 0
            result = library.wirehair_encoder_create_profile_ex(
                None, ctypes.cast(buffer, ctypes.c_void_p), len(view),
                block_bytes, ctypes.byref(profile), flags,
                ctypes.byref(output))
        _check("wirehair_encoder_create", result, library)
        encoder = _construct_owned_codec(
            cls, library, output, len(view), block_bytes)
        if not owned:
            encoder._borrowed_source = (view, buffer)
        return encoder

    @classmethod
    def _from_pointer(cls, library, pointer, message_bytes, block_bytes,
                      owner=None):
        return _construct_owned_codec(
            cls, library, pointer, message_bytes, block_bytes, owner=owner)

    def encode(self, block_id, capacity=None):
        block_id = _bounded_integer(block_id, "block_id", 0xffffffff)
        if capacity is None:
            capacity = self.block_bytes
        capacity = _bounded_integer(capacity, "capacity", 0xffffffff)
        if capacity == 0:
            raise ValueError("capacity must be positive")
        output = (ctypes.c_uint8 * capacity)()
        written = ctypes.c_uint32()
        result = self._library.wirehair_encode(
            self.pointer, block_id, ctypes.cast(output, ctypes.c_void_p),
            capacity, ctypes.byref(written))
        _check("wirehair_encode", result, self._library)
        if written.value > capacity:
            raise RuntimeError("wirehair_encode returned an invalid length")
        return bytes(output[:written.value])


class Decoder(_Codec):
    """Context-managed Wirehair decoder."""

    @classmethod
    def create(cls, message_bytes, block_bytes, library=None, profile_id=None):
        """Create a decoder, optionally for a trusted legacy profile id."""
        library = initialize(library)
        message_bytes = _bounded_integer(message_bytes, "message_bytes", 0xffffffffffffffff)
        block_bytes = _bounded_integer(block_bytes, "block_bytes", 0x7fffffff)
        if message_bytes == 0 or block_bytes == 0:
            raise ValueError("message_bytes and block_bytes must be positive")
        output = ctypes.c_void_p()
        if profile_id is None:
            result = library.wirehair_decoder_create_ex(
                None, message_bytes, block_bytes, ctypes.byref(output))
        else:
            profile = _profile_descriptor(profile_id, library)
            result = library.wirehair_decoder_create_profile_ex(
                None, message_bytes, block_bytes, ctypes.byref(profile),
                ctypes.byref(output))
        _check("wirehair_decoder_create", result, library)
        return _construct_owned_codec(
            cls, library, output, message_bytes, block_bytes)

    def decode(self, block_id, data):
        block_id = _bounded_integer(block_id, "block_id", 0xffffffff)
        view, buffer = _input_buffer(data)
        if len(view) > self.block_bytes:
            raise ValueError("encoded block exceeds block_bytes")
        result = self._library.wirehair_decode(
            self.pointer, block_id, ctypes.cast(buffer, ctypes.c_void_p), len(view))
        if result == Wirehair_NeedMore:
            return False
        _check("wirehair_decode", result, self._library)
        return True

    def recover(self):
        output = (ctypes.c_uint8 * self.message_bytes)()
        result = self._library.wirehair_recover(
            self.pointer, ctypes.cast(output, ctypes.c_void_p), self.message_bytes)
        _check("wirehair_recover", result, self._library)
        return bytes(output)

    def recover_block(self, block_id, capacity=None):
        block_id = _bounded_integer(block_id, "block_id", 0xffffffff)
        if capacity is None:
            capacity = self.block_bytes
        capacity = _bounded_integer(capacity, "capacity", 0xffffffff)
        if capacity == 0:
            raise ValueError("capacity must be positive")
        output = (ctypes.c_uint8 * capacity)()
        written = ctypes.c_uint32()
        result = self._library.wirehair_recover_block_ex(
            self.pointer, block_id, ctypes.cast(output, ctypes.c_void_p),
            capacity, ctypes.byref(written))
        _check("wirehair_recover_block", result, self._library)
        if written.value > capacity:
            raise RuntimeError("wirehair_recover_block returned an invalid length")
        return bytes(output[:written.value])

    def become_encoder(self):
        """Convert a completed decoder, transferring its ownership.

        The replacement Python wrapper is fully allocated before the native
        codec is irreversibly converted.  Ordinary allocation failures leave
        this decoder usable; if even failure cleanup cannot safely detach the
        candidate, the shared owner is closed instead.  Once native conversion
        succeeds, any exceptional transfer failure also closes the shared
        native owner instead of leaving a decoder wrapper that points at an
        encoder-mode codec.
        """
        pointer = self.pointer
        owner = self._owner
        encoder = Encoder._from_pointer(
            self._library, pointer, self.message_bytes, self.block_bytes,
            owner=owner)
        if encoder._owner is not owner or encoder.closed:
            encoder.close()
            raise RuntimeError("converted encoder did not adopt native ownership")

        try:
            result = self._library.wirehair_decoder_becomes_encoder(pointer)
        except BaseException:
            # A foreign-function call that raises has unknown native state.
            # Conservatively close both wrappers through their shared owner.
            encoder.close()
            raise

        if result != Wirehair_Success:
            if result == Wirehair_InvalidInput:
                # The native entry point validates this precondition before
                # mutating the decoder, so discard only the prepared wrapper.
                try:
                    _Codec._relinquish(encoder)
                except BaseException:
                    owner.release()
            else:
                # Other native failures may leave a failed/partially converted
                # codec and cannot be rolled back safely.
                encoder.close()
            _check("wirehair_decoder_becomes_encoder", result, self._library)

        try:
            self._relinquish()
        except BaseException:
            encoder.close()
            raise
        return encoder


def _run_example(message=None, block_size=32, library=None, max_packets=None):
    """Run a bounded loss/reorder recovery example; return a process status."""
    if message is None:
        message = ("Wirehair Python example: loss recovery with UTF-8 data. "
                   "Zażółć gęślą jaźń. ").encode("utf-8") * 3
    try:
        message = bytes(memoryview(message))
        if not message:
            raise ValueError("message must not be empty")
        block_size = _bounded_integer(block_size, "block_size", 0x7fffffff)
        if block_size == 0:
            raise ValueError("block_size must be positive")
        block_count = (len(message) + block_size - 1) // block_size
        if max_packets is None:
            max_packets = block_count + max(64, block_count // 2)
        max_packets = _bounded_integer(max_packets, "max_packets", 0xffffffff)
        if max_packets == 0:
            raise ValueError("max_packets must be positive")

        with Encoder.create(message, block_size, library=library) as encoder:
            with Decoder.create(len(message), block_size, library=library) as decoder:
                complete = False
                for block_id in range(max_packets):
                    if (block_id + 1) % 10 == 0:
                        continue
                    packet = encoder.encode(block_id)
                    if decoder.decode(block_id, packet):
                        complete = True
                        break
                if not complete:
                    raise RuntimeError(
                        "decoder still needs data after %d packet identifiers" %
                        max_packets)
                recovered = decoder.recover()
                if recovered != message:
                    raise RuntimeError("recovered message differs from input")
        print("Wirehair recovery succeeded (%d bytes)" % len(message))
        return 0
    except Exception as exc:
        print("Wirehair example failed: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(_run_example())
