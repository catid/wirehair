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

The build also produces a separate, non-LTO serialized K3 C boundary and five
additional tests. Its shared `Facade<K,Traits>` implementation retains K6's
independent/borrowed input, transactional allocating detach, descriptor
validation and permanent conflict-poison behavior. Only the benchmark K3
wrapper instantiates that boundary for external callers; the installed K6
implementation remains untouched. K3 uses `WHK3` and profile ID
`0x5748324b33544d31`, not an existing or retired WH2/K6 identity.

`small_serialized` tests K3 ownership/error/byte behavior; `small_serialized_k6`
instantiates the same template at K6 and checks installed K6 descriptor,
packet, feed, recovery and detach parity. `small_serialized_c` verifies the C
ABI, including decoder creation before any encoder exists. The remaining two
tests replay the neutral and retained K3 corpus through the external serialized
functions. All seven tests are correctness checks, not speed measurements.
When enabling sanitizers, instrument both C and C++ compile flags and the link
flags; the supplied archive must use the same sanitizers/backend configuration.

## Separate K3 serialized performance qualification

At source `5caa4a9`, the separate `Wh2K3SerializedCostR1.py` gate passed all
54 same-code timing controls and all 36 comparisons against actual WH1 and
current public WH2. It retained 2,488,320 fresh codec lifecycles, checking every
output. Time reductions versus WH1, spanning both measurement orders:

| Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---|---|---|
| 2 | 65.8-65.9% | 95.4% | 91.0% |
| 64 | 72.7% | 95.1% | 91.3% |
| 1280 | 75.2-75.5% | 93.9% | 92.2% |

These are serialized-prototype, static-call results on one GFNI-capable host,
using borrowed immutable input and three full source blocks. Encoder work
includes creation, descriptor output, 18 packets and free; decoder work
includes creation, feed through first success, recovery and free. Both decoder
streams reached success after three packets for every arm. This is not an
integrated-library, cold-start, shared-call, partial-tail, non-GFNI or all-K
speed result. This prototype gate did not change the normal library or defaults.

The immutable outcome is `/var/tmp/wh2-k3-serialized-cost-r1`. Its raw stream
SHA-256 is `35d3ce567924ae21418f0b16986e6cf89a02a162e0b2c64bab6d89036a250604`.
A Python 3.8 full replay verified all 601 receipt pins; a separately written
raw chronology, API ledger and confidence-interval audit reproduced all 90
decisions. The R0 run remains invalid because its neutral diagnostic truncated
the binary CPU identity at an embedded NUL. R1 corrected that transport and
added prelaunch validation without changing the candidate or performance gate.
Neither namespace may be rerun or its observations filtered or rescored.

The retained K3 recovery screen had no zero-overhead failures in 6,144 traces
and passed all 72 hard traces. Native and serialized payload replay confirmed
those same cases; this timing run adds no recovery-rate sample or paired WH1
recovery-rate comparison. Production integration must preserve the sealed
equations and explicit descriptor, then qualify the actual library separately.

The later opt-in integration is documented in
[`SMALL_WIRE_PROFILES.md`](../../SMALL_WIRE_PROFILES.md). Its separately frozen
`Wh2K3ProductionCostR0.py` actual-library gate at `ffa7739` also passed all
54 same-code controls and 36 WH1/current WH2 comparisons; see that document
for its distinct measurements and scope. `Wh2SmallProductionParity.cpp` reuses
this correctness corpus
through both the old serialized prototype archive and the new library,
comparing descriptors, packet bytes, prefix status and recovered bytes. It is
an engineering replay, not a timing or fresh recovery experiment, and does
not alter the frozen source harness or old outcome bundles.
