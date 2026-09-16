# Alignment-only small-profile encoder candidate

Issue `wirehair-hgiu`. Benchmark-only, pure GF(256); production unchanged.
The passing `Wh2SmallBasisAlignment` intervention motivates this distinct
candidate. It does not revive either rejected two-byte payload helper.

Only `CreateSmallEncoder` changes in a hash-pinned external source copy.
For block widths at least 256 bytes and divisible by 64, allocate 63 extra
bytes and pass the next 64-byte boundary to the existing evaluator. Keep the
original pointer in the same unique_ptr for deletion. Narrow/irregular-width
allocation sizes, codec layouts, equations, decoder, source-policy semantics,
and partial-tail padding are unchanged. Check size overflow before allocation.
At most 63 extra owned bytes and no extra allocation are required per encoder.

Build externally with CMake/Ninja; qualify native, `PORTABLE=ON`, and
`SANITIZE=ON` builds. Run sanitizer tests with leak/fake-stack checks enabled
and halt on undefined behavior. The shared existing WHV2 tests provide an
independent polynomial packet oracle, ownership, alias and lifecycle checks.
`AllocationTest.cpp` forces all four malloc-compatible raw mod64 offsets,
verifies the exact copied interior view and surrounding guards, original-owner
deletion, every OOM boundary, and full/partial packet parity across K3/K5/K8,
13 boundary widths, and both ownership policies. It manually poisons unused
storage/lifetimes under ASan. Shared-library comparisons run in both load orders.

No timing controller or promotion decision is included yet. A separate bounded
screen must test lifecycle cost, aligned-baseline retention, internal placement,
existing paths and WH1 comparisons. Neutral identity alone is not speed or
recovery-rate evidence; this candidate does not change recovery equations.

## Neutral qualification (2026-09-16)

GNU 13.3 native, portable-arithmetic and ASan/UBSan builds each passed 10/10
CTest cases. Sanitizer tests used `detect_leaks=1`,
`detect_stack_use_after_return=1` and `UBSAN_OPTIONS=halt_on_error=1`.
Each arm's forced-placement test passed 936 constructions, 3,384 injected OOM
boundaries and 14,976 packets. Tiny tail duplicates are retained in that roster;
these are engineering cases, not independent recovery samples. The existing
WHV2 tests cover 69 distinct width/tail shapes and all constructor routes,
including ownership detachment and an independent polynomial packet oracle.

The native pair retains the exact addresses, sizes and machine bytes of
`EncodeSmall`, `DecodeSmall`, `RecoverSmall` and public WHV2 encode/decode/
recover. `EncodeSmall` remains 5,468 bytes. Total `.text` grows by 128 bytes;
whole-text identity and performance retention are **not** claimed.

Builds: `/tmp/wh2-aligned-small-basis.Zu0KNnhY/{native,scalar,asan}`.
Native baseline SHA256:
`bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`.
Native candidate SHA256:
`59a98ab34920cf40c2cfba6aa26b63d8ff929e9ec0823050ab07d21bb2ee05c7`.
Generated candidate profile source SHA256:
`2189df2596002626418c8f812849d5f3b10bf9c7e410187a58e8d5a81c29623c`.

Retained `Testing/Temporary/LastTest.log` SHA256 values:

- native: `2aa80c0c09b1ca6936a3f4ecfc0800b9fb224ccfa22554e6b02d62c79e3181f7`
- scalar: `7bcfab2b8f5df5d1a508acef282b5f5c101e411c6657782d88cac607ea66127b`
- asan: `3fbf014d98edceae77bcf8f2557ac8212665e97fb48cfbbfd3da0591865f11db`

This checkpoint qualifies construction/packet safety, not promotion. Do not
run a previous candidate's spent timing controller against these libraries.

## Speed decision

The separate `Wh2AlignedSmallBasisScreen` at `047d20b` is terminal
**CONTROL_FAIL** in both natural-lifecycle load orders. Controlled repairs
passed, but lifecycle retention/control failures prevent promotion. See that
screen's README for independently audited results and limitations. This
candidate remains benchmark-only; no production or recovery change was made.
Follow-up `wirehair-6juk` is attribution work, not permission to rerun this
candidate or recycle its spent namespace.
