# Wirehair CI commands

The workflow entry points are ordinary local commands. Build directories and
install prefixes are explicit so static, shared, and instrumented trees do not
contaminate one another.

Fast GCC/Clang or MSVC lanes run the registered CTests, install Wirehair, build
and run the C-language consumer in `test/package`, and run the Python fake
binding tests. Shared lanes additionally load the installed library and binding
for a real loss, reorder, block-recovery, and decoder-to-encoder round trip.
The workflow passes `--strict`, which applies
`-Wall -Wextra -Wpedantic -Werror` with GCC/Clang or `/W4 /WX` with MSVC to the
library and installed consumer builds. Native shared-library coverage includes
structured creation, encode, and premature-recovery failures in addition to
the successful loss/reorder/conversion round trip. The strict Visual Studio
shared lane also runs the Visual Studio 2022 `dumpbin /exports` tool against
the installed DLL and requires exactly all 15 public `wirehair_*` C entry
points. Its complete export table is retained in `ci-logs` for failure
diagnostics. These VS17 jobs are pinned to the `windows-2022` runner image so
the image and requested Visual Studio generator cannot drift apart.
The hosted static GCC and Clang lanes additionally pass `--tool-coverage` to
build the excluded offline generators, V2 benchmark CLI, and `wirehair_whx`
validation harness. Their CTest run includes the heavy-matrix generator
regression, V2 benchmark CLI contract, and whx selftest with seed fixups
explicitly enabled.

```sh
python3 ci/run_ci.py matrix --linkage static --strict \
  --tool-coverage \
  --build-dir build/ci-gcc-static --install-dir build/install-gcc-static \
  --generator Ninja --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++

python3 ci/run_ci.py matrix --linkage shared --strict \
  --build-dir build/ci-gcc-shared --install-dir build/install-gcc-shared \
  --generator Ninja --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++
```

The sanitizer commands use compiler and linker instrumentation on the full
target graph. Leak detection is enabled by the workflow environment. Current
Linux kernels can place PIE mappings inside GCC TSan's fixed shadow region;
the TSan command launches CTest through `setarch -R` to prevent that runtime
collision, then runs both the normal legacy core test and the 24-thread cold
initialization test.

```sh
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
python3 ci/run_ci.py sanitizer --config Debug --generator Ninja \
  --cmake-arg=-DCMAKE_C_COMPILER=clang \
  --cmake-arg=-DCMAKE_CXX_COMPILER=clang++

TSAN_OPTIONS=halt_on_error=1:second_deadlock_stack=1 \
python3 ci/run_ci.py tsan --config Debug --generator Ninja \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++
```

Scheduled coverage fixes all codec seeds and runs small, medium, maximum-N,
large-block, deterministic-loss, and bounded fuzz recovery. It enables the
registered scheduled CTest profiles and runs `ctest -L scheduled` before the
additional direct and fuzz cases. Failures from `unit_test`,
`large_message_test`, or `whx` remain process failures and their logs are
written below the selected build directory.

```sh
CXX=g++ python3 ci/run_ci.py scheduled --generator Ninja \
  --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++
```

The QEMU lane requires Clang plus `qemu-x86_64` from `qemu-user` (or the static
equivalent). It runs the portable artifact on Conroe and a caller-selected
`-mavx2` artifact on Haswell. It does not claim that the AVX2 artifact is
portable to older CPUs. The gate requires exact dispatch reports: every SIMD
kernel disabled on Conroe, and only SSSE3 plus AVX2 enabled on Haswell.

```sh
CXX=clang++ python3 ci/run_ci.py qemu --build-dir build/ci-qemu
```

This command requires MinGW-w64 and cross-builds static and shared packages,
then builds the pure-C package consumer. It asserts that only the shared
consumer imports `libwirehair.dll` and that the DLL exports exactly the 15
public C entry points.

```sh
python3 ci/run_ci.py mingw --generator Ninja \
  --build-dir build/ci-mingw --install-dir build/install-mingw
```
