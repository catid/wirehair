# Wirehair Python binding

`wirehair-fec` is the standard Python distribution for Wirehair. New code
should use the correctly spelled import:

```python
import wirehair

message = b"application payload"
with wirehair.Encoder.create(message, block_bytes=1200) as encoder:
    first_packet = encoder.encode(0)
```

The historical `import whirehair` spelling remains supported. Both names
export the same classes, constants, and functions.

## Native-library policy

The wheel is pure Python (`py3-none-any`) and deliberately does **not** bundle
`wirehair.dll`, `libwirehair.so`, or `libwirehair.dylib`. Install a matching
Wirehair shared library from the same release through CMake, a system package,
or your application's deployment system. The binding checks major API version
2 and requires every symbol used by that release; an older 2.x library missing
newer entry points is rejected. Static archives cannot be loaded by this binding.

Discovery is deterministic and follows this order:

1. A non-empty `WIREHAIR_LIBRARY` is authoritative and names the exact shared
   library file to load. Failure does not fall back to another installation.
2. A non-empty `WIREHAIR_PREFIX` is authoritative and searches only its
   `bin`, `lib`, and `lib64` trees using platform library names.
3. CMake-prefix and active Python-prefix locations are searched.
4. The operating-system loader is queried (`find_library` and normal loader
   names).

Linux uses `libwirehair.so`/`libwirehair.so.2`, macOS uses
`libwirehair.dylib`/`libwirehair.2.dylib`, and Windows uses
`wirehair.dll`/`libwirehair.dll`. `wirehair.initialize()` validates the native
API version before a codec is created.

## Installation modes

Build and install the Python distribution with any PEP 517 frontend:

```sh
python -m pip install .
python -m pip install -e .  # editable source checkout
```

A shared or dual CMake installation also installs both import names under its
`python` directory for prefix-oriented deployments. Static-only CMake installs
omit all Python artifacts because they do not contain a loadable library.

The reusable `encode_into`, `recover_into`, and `recover_block_into` methods
write to caller-owned `bytearray` or writable contiguous `memoryview` storage.
Decode/recovery success is not authentication: verify a digest or MAC obtained
from trusted application metadata before using recovered bytes.

`Encoder.detach_input()` severs the native encoder from its source after the
recovery columns created during construction are ready. For borrowed input it
also releases the binding's retained buffer view, so the caller may resize,
modify, or destroy the source immediately after success. For owned input it
releases the native private message copy. The method is idempotent and packets
stay byte-identical, but systematic packet generation becomes slower because
the encoder must reconstruct those packets from recovery columns. Do not call
it concurrently with encode, reuse, conversion, or close.
