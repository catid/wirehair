# Out-of-line two-byte payload diagnostic

Benchmark-only follow-up `wirehair-6llt` to the rejected
`Wh2SmallPayload2` inline experiment. No production, profile, recovery, or
default-selection changes. Neutral qualification alone is not speed evidence.

The inline candidate enlarged the shared encoder by 830 bytes and shifted
later functions. This variant keeps the original five-argument call boundary
and moves the B2 decision and source-major two-byte dot product into a separate
translation unit, compiled without LTO. Wider blocks tail-dispatch to the
existing kernel; other dimensions and WH_COUNT retain their original routes.
The source-major loop reuses each coefficient's table base for both bytes.
This tests a concrete alternative to duplicating the fast path in the hot
encoder, not a retry of the spent inline timing experiment.

Use `cmake -S bench/Wh2SmallPayload2Outlined -B <fresh-external-build> -G Ninja
-DCMAKE_BUILD_TYPE=Release`, then build and run CTest. Also qualify separate
`PORTABLE`, `SANITIZE`, and `COUNT` builds. ASan/UBSan gates use leak detection,
stack-use-after-return detection, and halt-on-UB. The exhaustive independent
polynomial oracle and public-API parity fixtures are reused from the inline
experiment; its candidate and its controller are not reused.

Inspect the native `EncodeSmall` and surrounding symbols against the paired
baseline. A smaller or identical encoder is only code-shape evidence. A
future speed decision requires a new prospectively frozen bounded screen,
including wider-path retention, both policies/load orders, and valid A/A
controls. Do not run the old controller against this build.

## Native code-shape result

`CodeShape.py <baseline.so> <candidate.so>` verifies the complete machine bytes
of `EncodeSmall`, `DecodeSmall`, `RecoverSmall`, and the public WHV2 encode,
decode and recover functions. All six retain their addresses and sizes; only
the three expected payload call destinations in `EncodeSmall` change.
`EncodeSmall` returns from the inline candidate's 6,298 bytes to 5,468 bytes.
The outlined helper is 482 bytes and adds 497 bytes including text padding.

Whole-text isolation **failed**: 27 other bytes differ after accounting for
the four payload call targets. They include reordered constructor-cleanup
helpers, register choices in the standalone small facade, and calls to the
direct-systematic recovery function, whose custom section shifts by 496 bytes.
The diagnostic reports these differences; it does not normalize them away.
The new helper also adds a branch/tail jump for wider repairs. Neither an
unchanged hot encoder nor neutral tests establish retention or speed.

Native build: `/tmp/wh2-small-payload2-outlined-native.XqDmenwR`; other modes:
`/tmp/wh2-small-payload2-outlined-neutral.S6Us0ju0/{scalar,asan,count}`.
All four modes passed 11/11 CTests before timing, including the independent
1,835,008-product oracle and 1,248 mixed/fallback cases. Additional direct
dispatch tests cover 36 null no-ops, 272 cases with counts 1–17, and exact
equality of all six operation call/byte counters in WH_COUNT (including a
nonzero witness). Five code-shape synthetic tests pass on Python 3.8/3.12.

Native baseline SHA256:
`82f9cf8985ea0873263593dd1dd945bdfcdeb993e32f473977225b6aa24090ed`.
Candidate SHA256:
`d10c490244d0415c6daeb035264b6a6e9484185a0476d19685e3319cece883b3`.

## Prospective matched-workspace screen

Issue `wirehair-viwn`. The prior observer assigned both arm and output buffer
from `side`; reversing observation order did not exchange the buffer mapping.
This was a structural confound, **not a proven cause** of its failures. The
old worker and spent results remain unchanged and rejected.

This build generates a distinct worker from the hash-pinned original source.
Within a workload cell, all arms now use one source allocation and one output
workspace. Their profiles, decoder repair bytes, and own-first-success counts
are copied into common staging addresses before WORK. Prebuilt borrowed handles
outlive no staging allocation. Output reset, staging, and independent output
verification are outside WORK. Staging copies still read the arm-specific
fixture allocations, so cache preparation is not claimed address-identical.
These changes remove fixed arm-to-workspace mapping inside WORK; they do not
guarantee unbiased timings or cure all noise.

After final neutral tests and independent review, run only
`python3 -B bench/Wh2SmallPayload2Outlined/Screen.py --run <native-build>`.
Sole namespace: `/var/tmp/wh2-small-payload2-outlined-screen-r0`. No retry,
trimming, resizing, pooling, or rescuing failed controls. The new wrapper pins
both controller modules, all experiment files, generated worker and libraries.
Like the original diagnostic, it is not full transitive build/runtime provenance.

Freeze the original 96-cell roster (K3/K5/K8, B2 full/partial, B64/B1280,
two source policies, four metrics), CPU50, batch32, 12 replicates, five
comparison pairs, two side orders, and separate DSO load-order processes.
All A/A t11 95% intervals must lie inside the reciprocal 2% envelope. Every
B2 encoder/repair C/B upper bound must be below 1; other C/B bounds must be
below 1.02, only a screening tolerance. Report C/WH1 separately. Both load
orders must pass for further qualification; no production, ordinary-default,
preserved-certified-path, universal-speed or recovery-rate claim follows.
The explicit K5/K8 routes are still not ordinary defaults.
