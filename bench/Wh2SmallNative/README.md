# K3/K6 native correctness harness

This standalone build tests the private compile-time small-block core. It does
not change the public library, select a wire profile, or measure performance.
The K3 fixture generator authenticates the retained `.82` recovery bundle and
reconstructs its sealed table without rerunning the selector or loss campaign.

```sh
cmake -S bench/Wh2SmallNative -B /tmp/wh2-small-native-build \
  -DWH2_SMALL_LIBRARY=/absolute/path/to/qualified/libwirehair.a
cmake --build /tmp/wh2-small-native-build
ctest --test-dir /tmp/wh2-small-native-build --output-on-failure
```

Use a fresh, task-owned build directory. This harness needs the original local
`/var/tmp/wh2-k3-thue-morse-r0` evidence; missing or changed evidence fails closed.
The ordinary library and installed package do not depend on this directory.

The private `gf256.h` compilation settings must match the linked library's
settings. For an `ANDROID`-selected portable-GF library, pass
`-DWH2_SMALL_PORTABLE=ON`. This checks the portable backend, **not** an Android
platform or real non-GFNI-host performance. Sanitizer compile and link flags
must likewise match the selected instrumented library.

`--neutral` checks K3/K6 widths and partial tails, field products, allocation
failures, aliases, dependent/conflicting packets, and K6 parity with the
untouched original core. `--corpus` checks all 7,774 retained K3 payload cases,
including the 6,144 fresh and 72 hard traces, historical original-width
prefixes, and every triple in the development/seam windows. Neither mode
compares speed with WH1 or establishes a new failure-rate sample.
