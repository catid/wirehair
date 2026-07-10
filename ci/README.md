# Wirehair CI commands

The workflow entry points are ordinary local commands. Build directories and
install prefixes are explicit so static, shared, and instrumented trees do not
contaminate one another.

## Architecture portability budget

The per-change workflow spends at most 45 additional runner-minutes on the two
failure modes most likely to be affected by an ordinary patch:

- `ubuntu-24.04-arm` gets 25 minutes for a native strict ARM64 shared build,
  the NEON GF256 differential/self-test, a byte-exact golden round trip, and
  installed C/C++/plugin consumers.
- `ubuntu-24.04` gets 20 minutes for an emulated ELF32 i686 strict build, the
  checked-allocation boundary/failure tests, and a golden round trip.

The weekly/manual scheduled workflow reserves at most 60 additional
runner-minutes for lower-frequency platform boundaries: 30 minutes for a
single Ubuntu job that runs the s390x big-endian heavy-window oracle and the
RISC-V scalar smoke round trip, plus 30 minutes for native macOS ARM64 package
and Mach-O ABI validation. Keeping the two cross targets in one scheduled job
avoids paying for the same QEMU/toolchain setup twice.

Linux portability jobs select GCC 13 executables and fail unless the compiler
reports the Ubuntu 24.04 `13.3` family. Emulated jobs likewise fail outside the
Ubuntu 24.04 QEMU `8.2` family. Security patch revisions may move within those
families. The macOS job uses the architecture-specific `macos-15` runner,
selects `/Applications/Xcode_16.4.app`, and requires Xcode `16.4` build `16F6`.
Every lane prints and validates the actual object class, machine, byte order,
pointer width, and selected GF256 backend; a mislabeled runner or accidentally
native artifact cannot silently satisfy a cross lane.

Developers on Ubuntu 24.04 can reproduce every Linux target after installing
the corresponding GCC 13 cross compiler and `qemu-user`:

```sh
python3 ci/run_portability.py cross --target aarch64 \
  --build-dir build/portability-aarch64 \
  --install-dir build/install-portability-aarch64
python3 ci/run_portability.py cross --target i686 \
  --build-dir build/portability-i686
python3 ci/run_portability.py cross --target s390x \
  --build-dir build/portability-s390x
python3 ci/run_portability.py cross --target riscv64 \
  --build-dir build/portability-riscv64
```

The QEMU ARM64 variant also installs the library and runs its external package
consumers, making the native ARM acceptance path reproducible without access
to an ARM host. On a `macos-15` ARM runner with Xcode 16.4 selected, the native
equivalent is:

```sh
python3 ci/run_portability.py macos-arm64 \
  --build-dir build/portability-macos-arm64 \
  --install-dir build/install-portability-macos-arm64
```

The Mach-O linker allowlist is generated from the same versioned
`abi/wirehair.map` manifest used by ELF, with Darwin's ABI underscore added at
configure time. The macOS lane then checks the installed dylib with `nm` and
rejects both missing API symbols and leaked C++ implementation symbols.

Repository data integrity runs once on Linux, outside the compiler/linkage
matrix. It runs the dense-count, precode-result/ranking, and byte-metric
validator unit tests; lints the heavy-matrix documentation; validates every
tracked precode CSV and error-log invariant; and exercises the byte-ledger
tools with bounded 1--10 MiB working sets and three timing repeats. Any
malformed artifact stops the lane at the validator that reports its file and
row. The complete hosted gate is also available locally:

```sh
python3 ci/run_ci.py data-integrity
```

## Hosted coverage ownership

Every compiler/linkage cell still builds the complete target graph, runs the
compiler-sensitive CTests and quick large-message profiles, installs its own
artifact, and builds and executes an external C package consumer. The following
map assigns gates that do not gain coverage from matrix repetition to one
explicit owner:

| Gate | Hosted owner | Other matrix cells |
| --- | --- | --- |
| Table/result/byte-metric validators | `data-integrity` | Never repeated |
| CI-runner and fake Python-binding unit tests | Linux GCC shared (`--repository-gates`) | Skipped |
| Core C/C++ CTests and quick large-message profile | All six build cells | Run for their own compiler/linkage |
| Nested static/shared/dual package producers and build-policy suite | Linux GCC shared (`ci-once` CTest label) | Excluded with `ctest -LE` |
| Installed external C consumer | All six build cells | Run against each cell's install |
| Installed native Python auto-discovery plus loss/reorder/recovery/conversion | Each shared cell | Static cells have no loadable library |
| Reproducible Python wheel/sdist, isolated/editable install, and uninstall | Linux GCC shared (`--repository-gates`) | Pure wheel remains platform-neutral |
| ELF ABI/export checks | Linux GCC/Clang shared | Non-ELF cells skip them |
| MSVC `dumpbin` ABI/export check | Windows shared | Other cells use their native ABI gates |
| Offline generator/V2 CLI smoke | Windows static (`--tool-smoke`) | Built and run once with MSVC |
| ASan+UBSan+LSan and TSan | Dedicated sanitizer jobs | Nested uninstrumented producers excluded |
| MinGW package/import/export checks | `mingw-cross` | Run once for static and shared PE artifacts |
| Native ARM64/NEON and emulated i686 boundaries | Dedicated per-change portability jobs | Not repeated in the six-cell matrix |
| High-N, large-block, coverage-guided fuzz soaks, and lower-frequency QEMU portability | Scheduled workflow | Ordinary CTest retains the fixed-corpus 10,000-mutation fuzz smoke |

The `ci-once` label belongs only to the three nested package producer tests and
the multi-build policy test. `ci/run_ci.py matrix` excludes that label unless
`--repository-gates` is present, and the sanitizer lane always excludes it
because those child configurations do not inherit sanitizer instrumentation.
CTest uses the lane's `--jobs` value while retaining `--output-on-failure` and
the workflow uploads top-level plus nested configure/test diagnostics on every
failure.

Fast GCC/Clang or MSVC lanes install Wirehair and run the coverage assigned
above. Shared lanes additionally load the installed library and binding for a
real loss, reorder, block-recovery, and decoder-to-encoder round trip. The
single Linux package owner also relocates a nested-libdir install, compiles and
runs pure-C consumers from both normal and `--static` pkg-config queries, and
cross-checks pkg-config/CMake versions plus installed license/manifests.
The workflow passes `--strict`, which applies
`-Wall -Wextra -Wpedantic -Werror` with GCC/Clang or `/W4 /WX` with MSVC to the
library and installed consumer builds. Native shared-library coverage includes
structured creation, encode, and premature-recovery failures in addition to
the successful loss/reorder/conversion round trip. The strict Visual Studio
shared lane also runs the Visual Studio 2022 `dumpbin /exports` tool against
the installed DLL and requires exactly the public C entry points in the
versioned `abi/wirehair.map` allowlist. Its complete export table is retained
in `ci-logs` for failure diagnostics. These VS17 jobs are pinned to the
`windows-2022` runner image so the image and requested Visual Studio generator
cannot drift apart.

```sh
python3 ci/run_ci.py matrix --linkage static --strict \
  --build-dir build/ci-gcc-static --install-dir build/install-gcc-static \
  --generator Ninja --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++

python3 ci/run_ci.py matrix --linkage shared --strict \
  --repository-gates \
  --build-dir build/ci-gcc-shared --install-dir build/install-gcc-shared \
  --generator Ninja --cmake-arg=-DCMAKE_C_COMPILER=gcc \
  --cmake-arg=-DCMAKE_CXX_COMPILER=g++
```

## CI-time comparison

The successful hosted baseline
[`29113415903`](https://github.com/catid/wirehair/actions/runs/29113415903)
used 18.67 summed job-minutes across the nine jobs present in that run: 12.28
minutes in the six build matrix cells, with a 2:50 matrix critical path, and a
5:33 baseline-set critical path owned by the sanitizer job. Those numbers are
wall-clock differences from the public job timestamps and exclude billing-
rounding policy. Data-integrity and architecture-portability jobs added by
separate changes after that run are outside this optimization comparison.

On the same 128-CPU development host and the fully enumerated 32-test tree at
the time of the optimization, the old serial full-CTest schedule took 70.26
seconds per cell. The new two-worker owner schedule took 46.50 seconds, while
a non-owner's 28-test compiler-sensitive subset took 31.93 seconds. Across six
matrix cells, the isolated CTest portion dropped from 7.03 to 3.44 summed
minutes (51%).

The authoritative after run
[`29128713010`](https://github.com/catid/wirehair/actions/runs/29128713010)
at commit `06d80d1` passed all 12 jobs. Comparing the same six build matrix
cells, summed time fell from 12.28 to 11.47 minutes (6.6%), while their critical
path moved from 2:50 to 2:52. The same nine-job baseline set used 20.47 summed
minutes and retained a sanitizer-owned 5:35 critical path, versus 18.67 and
5:33 before. That raw total is not a fixed-workload comparison: the audit grew
the registered graph from 32 to 41 tests and lengthened the TSan, package, and
ABI gates after the baseline. The separately added repository-integrity,
ARM64, and i686 jobs contributed another 1.77 summed minutes, for 22.23 across
all 12 jobs, without extending the 5:35 workflow wall time. The observed
hosted result replaces the earlier linear planning estimate; configure, build,
install, and expanded-test costs dominate more than the isolated CTest timing
suggested.

For an internal pull request, feature-branch pushes previously triggered both
`push` and `pull_request` copies of the entire workflow. `push` is now limited
to `master` and release tags; the pull request owns feature commits, and the
eventual merge commit owns the one `master` push. Thus the same commit no longer
runs two full workflows. Direct `master`, release-tag, and manual runs remain
available.

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
consumer imports `libwirehair.dll` and that the DLL exports exactly the public
C entry points in `abi/wirehair.map`.

```sh
python3 ci/run_ci.py mingw --generator Ninja \
  --build-dir build/ci-mingw --install-dir build/install-mingw
```
