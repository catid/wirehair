# Large-message end-to-end profiles

`large_message_test` exercises the public API as a transport would: it encodes
all original packets, drops a deterministic subset, delivers the remainder in
a seeded permutation, emits bounded repair packets, decodes, recovers the full
message, and verifies every recovered byte.  The final original packet can be
short, so partial-final handling is covered rather than inferred.

## Resource policy

Every profile declares `WIREHAIR_LARGE_MAX_PAYLOAD_BYTES`.  The executable
rejects the profile before allocation unless
`N * block_bytes <= MAX_PAYLOAD_BYTES`.  This quantity is the reproducible
payload-pressure policy, not an estimate of process RSS: the test keeps input
and recovered buffers, the codec has its own working state, and sanitizers add
redzones.  Successful and failing runs print elapsed milliseconds and peak RSS
on platforms with `getrusage`; the shell runner also records `/usr/bin/time`
elapsed seconds and maximum RSS when available.

The executable's safe standalone default is `N=256`, `block_bytes=257`, a
113-byte final block, five losses, a 64-repair cap, a 64 MiB payload-policy cap,
and a 60-second internal deadline.  Invalid numbers, resource-policy failures,
terminal decode errors, repair-cap exhaustion, and deadline failures all exit
nonzero.

## Profiles

| Profile | N | Block bytes | Final bytes | Payload policy | Purpose |
| --- | ---: | ---: | ---: | ---: | --- |
| `fast_partial` | 256 | 257 | 113 | 1 MiB | Default CTest encode/loss/reorder/repair/recover verification |
| `n64000_byte` | 64000 | 1 | 1 | 64 KiB | Maximum supported N and one-byte blocks |
| `packet_partial` | 4096 | 1280 | 733 | 8 MiB | Packet-sized blocks and partial final packet |
| `large_block_partial` | 32 | 1048576 | 1048559 | 32 MiB | Large blocks and bounded payload memory pressure |

The quick suite repeats the success profile and compares its packet-order and
loss line, then asserts four failure paths: an injected API decode error,
repair-cap exhaustion, an injected internal timeout, and pre-allocation
resource rejection.  It requires the exact status and expected diagnostic for
each and rejects sanitizer diagnostics, so an unrelated crash cannot satisfy
an expected-failure contract.  CTest's own
`TIMEOUT` remains the outer bound for a call that stops making progress inside
the codec.

## Running

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON
cmake --build build --target large_message_test -j
test/run_large_message_profiles.sh build/large_message_test quick
test/run_large_message_profiles.sh build/large_message_test scheduled
```

The default CTest registration contains only the quick resource-bounded
profile and its expected-failure contracts.  Configure with
`-DWIREHAIR_ENABLE_SCHEDULED_TESTS=ON` to register the three scheduled
profiles, then select them with `ctest --test-dir build -L scheduled
--output-on-failure`.

For sanitizer validation, configure a separate build with the project's
supported ASan/UBSan flags.  Run the quick suite and the high-N/packet profiles;
the 32 MiB large-block case is intentionally scheduled because sanitizer RSS
is substantially larger than its payload-policy value.

## Reference bounds

Reference measurements on 2026-07-10 with GCC 13.3 Release and
ASan+UBSan/RelWithDebInfo builds were:

| Profile | Release time / RSS | ASan+UBSan time / RSS |
| --- | ---: | ---: |
| `fast_partial` | 0.00 s / 3 MiB | 0.00 s / 10 MiB |
| `n64000_byte` | 0.06 s / 27 MiB | 0.12 s / 40 MiB |
| `packet_partial` | 0.01 s / 29 MiB | 0.04 s / 40 MiB |
| `large_block_partial` | 0.09 s / 203 MiB | 0.48 s / 239 MiB |

These are observations, not pass thresholds; CTest's declared timeout and the
payload-policy cap are the portable bounds.  CI hosts should retain runner
output so time/RSS trends remain visible.
