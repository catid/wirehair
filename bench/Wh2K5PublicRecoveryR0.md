# Installed K5 retained recovery R0

Prospectively frozen in `wirehair-sxvz.16.1.20.82.3.1.6.2` before source edits.
Protocol `wirehair.wh2.k5-public-recovery-r0`; sole output directory
`/var/tmp/wh2-k5-public-recovery-r0`. This is recovery qualification, not timing,
a new holdout, a default change, or an all-K claim.

The candidate is the installed explicit WHV2 profile `0x80070c81bfe375f1`,
seed zero, with the unchanged sealed K5 GF(256) companion pair. Six arms use
actual public APIs: certified WH2 independent, K5 independent, WH1 owned,
certified WH2 borrowed, K5 borrowed, and WH1 borrowed. Both WHV2 profiles use
the same public call adapters. No prototype codec is linked.

The roster is unchanged: 6,144 retained loss traces, 72 hard cases, and 57
original-width historical cases, containing 56,234 packet IDs. Widths are
2, 64, and 1280 bytes. Native, portable arithmetic, and ASAN/UBSAN backends
repeat the same cases; policies and backends do not increase sample size.
No historical prefix is extended.

Each arm uses one real and five basis encoders, all destroyed before any
receiver. Independent polynomial coefficients, basis-derived every-arm
payloads, and full-rank prefixes must agree with actual packet bytes and
first successful decode. Each successful message is recovered twice with
byte/length/guard checks. Complete attempted API ledgers, immutable inputs,
policy parity and backend parity are mandatory.

Full records must reproduce the authenticated earlier prototype native raw
(`d6951f2594434c739776d0c37e1e4ed3b44317d91a49819dc832c28fee837303`)
after only replacing its candidate descriptor with the installed WHV2 identity.
Expected retained OH0 failures are K5 11/6144, WH1 122/6144, certified WH2
439/6144. All K5 traces recover with one extra packet; every 512-trace cell
is at most 1%. K5 fixes all 122 WH1 failures but introduces 11 different
failures. This is not per-trace or per-cell dominance or a population bound.

Build qualification reuses `Wh2K5PublicCostR0.py`'s exact installed archive,
18 producing objects, compiler recipes/dependencies, producing Git blobs,
and preserved test logs. Inputs are frozen before use and checked after
compilation. Driver ISA flags match each archive; its compiled GF context size
must equal the linked symbol size before any codec call. Original sources,
builds, and sealed namespaces are untouched.
Neutral checks cover 36 cases, exhaustive field arithmetic, four late-call
cleanup failures, positive claim authentication, and invalid CLI/claim cases.
Python 3.8 and 3.12 tests and repeated source-reading passes precede launch.

Workers have CPU 120s, wall 150s, and address-space 384MiB caps; only ASAN's
shadow mapping is exempt from the address-space limit. The controller allows
180s per worker and 600s overall. Raw output is limited to 64MiB per backend,
stderr to 64KiB, and the sealed bundle to 256MiB. Allocator/loader overrides
must be absent. Exact sanitizer settings are:

```
ASAN_OPTIONS=detect_leaks=1:detect_stack_use_after_return=1:halt_on_error=1
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
```

Commit and push qualified sources before creating the receipt and claiming
the sole run. Require exact replay and a separately written complete
raw/arithmetic/provenance audit before advancing HEAD or pinned reports.
No rerun, retuning, sample filtering, or speed promotion follows from this gate.
