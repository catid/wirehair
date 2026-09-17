# One-allocation small-decoder storage

Issue `wirehair-sxvz.16.1.20.85`. Benchmark-only candidate, **neutral checks
passed; subsequent speed screen CONTROL_FAIL; candidate not retained**.
Production is unchanged. This is a new allocation
reduction, not a retry of dormant-core lifetime, alignment or serializer changes.

## Change

Current installed small WHV2 decoders allocate a public facade, private decoder,
and `(K+1)*block_bytes` packet slab separately. The generated candidate combines
the latter two allocations. Public construction goes from three allocations to
two; private core construction goes from two to one. The public facade, unused
general core, encoder storage, row algebra, tables and profile defaults are
unchanged. Both source policies retain their existing semantics.

`Generate.py` checks the exact original core SHA before changing only decoder
storage management in an external copy. It checks addition overflow, requests
one raw scalar allocation, placement-constructs a nonthrowing decoder and a
trailing uninitialized byte array, then sets the original RHS/scratch pointers.
The decoder has a raw slab pointer and a class-specific unsized delete that
releases the original scalar allocation. It never frees the interior slab
separately or supplies `sizeof(Decoder)` as the full deallocation size.

The explicit byte-array placement starts its lifetime without relying on
implicit-lifetime wording. Supported modern GCC/Clang nonallocating placement
array new has no cookie. A null result is checked before pointer setup. The
existing dimension check bounds multiplication; the additional subtraction
check bounds object-plus-slab addition. No 32-bit execution is claimed.

Decoder size is 104/136/200 bytes for K3/K5/K8 on the tested ABI. Its fundamental
alignment is sufficient for ordinary allocation. Consequently the trailing slab
starts eight bytes past a 16-byte boundary on an aligned base, unlike an
independently allocated slab. This is an inherent placement effect of coallocation,
not an extra alignment optimization or evidence of faster execution. Existing
unaligned GF256 kernels and odd-width/tail tests remain essential.

All four small-profile translation units see the same generated private header.
Baseline and candidate tests link separate archives; shared parity uses separate
hidden-internal-symbol DSOs. No baseline/candidate template definitions are mixed
within a static executable. Test-only allocation replacement and forced-OOM
facades never enter the performance libraries.

## Neutral evidence

Final fresh builds: `/tmp/wh2-small-decoder-storage-final.6dWpqN24/`.

| Configuration | CTests | Result |
| --- | ---: | --- |
| GCC 13 native | 21 | Pass |
| GCC 13 forced portable | 21 | Pass |
| GCC 13 ASan + UBSan, leaks and fake stacks | 21 | Pass |
| Clang 18 C++, GCC 13 C consumers | 21 | Pass |

The suite includes both arms' complete WHV2 K3/K5/K8 tests, WHK3 standalone C++
facade tests, actual shared-library C consumers, profile and borrowed-source
tests, separately hook-enabled facade faults, allocation/lifetime probes,
two DSO load-order parity tests and the 48-case certified compatibility stream.
All four configurations preserve the 2,180,292-byte certified stream, SHA256
`2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.

The new storage probe checks 114 private-core cases and 456 public decoder
lifecycles per arm/configuration, across K3/K5/K8, 19 widths from 1 to 4096,
full/one-byte tails, low/distant repairs and both source policies. Block1's two
tail entries coincide; these are coverage cases, not independent rate samples.
It checks exact allocation bases/sizes/kinds, all matching frees, no packet-time
allocations, object/slab and cross-boundary aliases, reverse pivots, partial
padding, dependent conflicts, duplicates and repeated recovery into freshly
reset guarded buffers. Public repairs retain systematic fallback and their own
first-success endpoint. Every constructor allocation is fault-injected:
1,596 baseline and 1,026 candidate OOM cases, with complete unwind accounting.

Existing tests are not weakened to ignore allocation checks. Generated copies
change WHV2's exact decoder count from three to two, and WHK3's exact count plus
two exhaustive decoder loops from three to two. All other bytes are preserved;
source hashes and four Python synthetic tests enforce the generation boundaries.
Those tests pass under Python 3.8 and 3.12.

Initial neutral artifacts remain under
`/tmp/wh2-small-decoder-storage.lknR9awQ/`. The first Clang direct-core allocation
probe failed its count/base assertion while candidate and public tests passed.
A non-inlined factory boundary fixes the test's allocation-observation assumption;
public library constructors already have a separate translation-unit boundary.
The original failure log is retained as `clang-initial-tests.log`. This is not a
demonstrated codec defect or a timing result. Independent review also added the
placement-array null guard, fresh recovery canaries, adjacency aliases and
standalone-facade coverage before final qualification.

The initial library list placed Profile before the small facade objects. Before
any timing, the final build restores the established screening-library link
order. The final native baseline exactly reproduces the retained unchanged DSO:
`bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`.

Generated core SHA256:
`d7dd701e28f315febe70883c0e618525fb5ab5384eccbb4edd04541352d6f906`.

| Configuration | Candidate DSO SHA256 |
| --- | --- |
| Native | `15d27877bd0afa7178f1f410f13acf1114f45b44eec4f512b388eaa0b434b237` |
| Portable | `46dc509c85da5408b313d514e70ac871c7f9b624211cb86f23f4f68448aa37aa` |
| Sanitized | `e3a30257646fb01963d2ef5ccc8302f36262c24cc8e03b97dda2f50c3fc36bce` |
| Clang | `564b586ac357e43dfd73b43fde99fd38ab0fcc68a31fd28ef9db0d80844fc3eb` |

`NEUTRAL.sha256` in the final build root records 1,006 files: build/generated
artifacts, compiled source/header dependencies, relevant build/test definitions,
final and preliminary logs, and the initial Clang failure log. Its SHA256 is
`8b619caafe0c14366bfde96ce371740d19cd7adf8fb4e56fe8e23d55acbf1b2e`.
All entries were rehashed successfully after publication. This document is not
part of that manifest, and the manifest is not a complete compiler/linker/runtime
input seal. Repeated independent source and final artifact review is clean.

These are neutral compatibility and ownership results, not new recovery samples,
speed qualification, complete transitive toolchain provenance or an all-K claim.
The prospectively frozen decoder lifecycle/retention/WH1 screen
(`wirehair-sxvz.16.1.20.85.1`) is now terminal CONTROL_FAIL in both loading
orders: 10/16 failed same-code bounds and 141/168 failed retention bounds,
respectively. Nominal decoder means improved 0.71–1.77%, but failed controls
and retention prevent promotion. All 907,200 rows and both decisions were
independently reconstructed before result-documentation advancement; see
[the full screen result](../Wh2SmallDecoderStorageScreen/README.md).
This candidate is not retained; no retry, alignment rescue or production
adoption. The allocation reduction is correctness-qualified only, not a
demonstrated acceptable speedup or new recovery-rate result.

Reproduce correctness only in a **new external** directory:
`cmake -S bench/Wh2SmallDecoderStorage -B <fresh-dir> -G Ninja`, build and run
CTest. Use `-DPORTABLE=ON` or `-DSANITIZE=ON` for those modes. Never rebuild or
rewrite retained qualified inputs or historical screen artifacts.
