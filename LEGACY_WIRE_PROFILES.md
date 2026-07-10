# Legacy wire profiles

Wirehair's legacy packet ID is an input to a deterministic equation generator.
Changing a dense or peel seed therefore changes the wire contract even when the
C ABI and `WIREHAIR_VERSION` do not change.  Raw packets contain only an ID and
payload, so they cannot identify or authenticate that contract themselves.

Applications that persist or exchange packets across library builds should put
the numeric `profile_id` in authenticated or otherwise trusted framing, call
`wirehair_wire_profile_init()`, and use the explicit profile create APIs.  They
must also authenticate the original message with an application digest or MAC.
A profile identifier is not an integrity check.

## Supported profiles

| Profile | Numeric ID | Canonical SHA-256 input | Seed-table fingerprint | N=2373 repair fingerprint | N=5550 sample fingerprint |
| --- | --- | --- | --- | --- | --- |
| `WIREHAIR_LEGACY_PROFILE_PRE_FIXUP` | `e1b9f77f1c90f680` | `wirehair:legacy-v2:pre-fixup:557c00c707a4b6a51db312c113b8036dadbe132e` (`e1b9f77f1c90f680dd4dc6f94167e72b98a2e766e01df5975f9ef680b957e638`) | `ea803c84c587877a` | `a6c1a4d75faec4d8` | `18497889fd2e5fd7` |
| `WIREHAIR_LEGACY_PROFILE_FIXUPS_2026_07` (`CURRENT`) | `4d241359db07bb07` | `wirehair:legacy-v2:fixups:63f759171e904285d5e9661125b03942d6a396be` (`4d241359db07bb07523d234e280b528d9066de741625a5fa40b6bee70c005479`) | `4294485edd1cc4dd` | `a9020b085ecd5dfa` | `a87074dc238666c0` |

The anchor commits identify the intended, verified packet bytes, originally
reproduced on little-endian hosts.  Those revisions contained a broken
big-endian heavy-window branch and are not supported as separate big-endian
profiles.  The named profiles now produce the same pinned bytes on both
little- and big-endian hosts; strict s390x emulation runs the golden/profile
and randomized recovery suites.

The seed fingerprint is FNV-1a over, for every `N` from 2 through 64000, the
little-endian tuple `(N, dense_count, dense_seed, peel_seed)`.  The repair
fingerprint is FNV-1a over `(little-endian packet_id, one-byte packet)` for the
2373 repair-only packets with IDs 2373 through 4745 and the message generator
pinned by `LegacyCoreTest`.  These are compatibility regression values, not
cryptographic hashes.  The N=5550 fingerprint covers four pinned repair IDs at
a point where both peel and dense exact-N seeds differ.

`LegacyCoreTest` also pins a wider corpus of five repair IDs at 301 boundary
and deterministic pseudo-random block counts.  Its pre-fixup fingerprint is
`d69a6563d81cf247` and its current-profile fingerprint is `26898d9bc72d6f10`.
That exact-multiple, one-byte-block pre-fixup corpus was reproduced byte-for-
byte with the `557c00c` anchor, covering later row-generator optimizations as
well as the seed-table split.

The profiles deliberately use the corrected full-width packet-ID comparison
for final-block sizing.  The `557c00c` encoder truncated a repair packet when
its ID was at least 65536, its low 16 bits equaled `N - 1`, and the message had
a short final block; its decoder already used the full-width rule.  Reproducing
that encoder bug would make profile output depend on an accidental ID alias and
would not match the anchor decoder's contract, so those truncated packets are
not a supported profile.  A golden fixture pins a seven-byte repair packet at
ID 65543 for `N=8`, where the anchor bug returned only four bytes.

Both profiles freeze all equation-affecting choices:

- the PCG32 transition/output function and row seeding in `PCGRandom`;
- peel weight, column iteration, mix-column generation, and packet-ID mapping;
- the complete dense-count, dense-seed, and peel-seed tables;
- the pre-fixup profile's absence of exact-N overrides, or the current profile's
  768 peel and 163 dense exact-N overrides;
- Shuffle-2 dense construction, the shipped 6x18 heavy matrix, and field
  polynomial `0x14D`.

Builds compiled with the private `WH_SEED_KNOBS` experiment macro deliberately
permit seed, dense-count, and `WH_PEELCAP` peel-weight overrides.  Nondefault
`WH_HEAVY_ROWS` or `WH_HEAVY_COLS` likewise change the equation geometry.
Defining the legacy `CAT_IDENTITY_LOWER_RIGHT` experiment replaces part of the
shipped heavy matrix and also changes packet equations.
These builds do not implement the named wire profiles: explicit profile
initialization/creation returns `Wirehair_UnsupportedPlatform`, while the raw
experiment APIs remain available but non-conforming.  Do not use such binaries
to persist or exchange packets under a standard profile identifier.

The unframed APIs (`wirehair_encoder_create*` and `wirehair_decoder_create*`)
are fixed to `WIREHAIR_LEGACY_PROFILE_CURRENT`.  They cannot detect packets
created with another profile.  Explicit APIs reject malformed or unknown
trusted descriptors before creating a codec.

## History inventory and support policy

There were no Git release tags at the 2026-07-10 audit point.  The initial
2018 snapshot at `058c0e6` differs from the supported pre-fixup anchor
`557c00c` and is not assigned a supported profile by this library.  Exact-N
correction tables were then accumulated on `master` through the 45
equation-table commits below.  Every intermediate
snapshot is a distinct equation contract and may have been consumed directly
from `master`, but none was a tagged release.  The current library intentionally
does not pretend it can reproduce those 45 partial tables: their identifiers
are unsupported and must be rejected.  A deployment that used one must retain
that exact binary/source snapshot, recover under it, verify an application
digest, and re-encode with a supported profile.

```text
ec226e3 2026-06-03 Add overhead seed fixups
ff739c7 2026-06-03 Add optimization experiment knobs and fixups
d616766 2026-06-11 Extend verified overhead fixups
f61231f 2026-07-03 Apply triple-gated [32000,64000] seed fixups
4165c5c 2026-07-07 Add matrix-only seed screening
f58dfbf 2026-07-07 Add bb0-screened overhead fixups
9496ac4 2026-07-07 Add top-severe overhead fixups
d6e55fd 2026-07-07 Add next severe overhead fixups
c094128 2026-07-07 Add additional severe overhead fixups
c967448 2026-07-07 Add 39304 severe overhead fixup
cece00d 2026-07-07 Add 56052 severe overhead fixup
e1e03da 2026-07-07 Add 49065 severe overhead fixup
ec8c76c 2026-07-07 Add 40915 severe overhead fixup
37e345e 2026-07-07 Add more severe overhead fixups
3635d3f 2026-07-07 Add 44538 severe overhead fixup
b52b851 2026-07-07 Add 32204 severe overhead fixup
78a89d1 2026-07-07 Add top5q severe overhead fixups
1a5115c 2026-07-07 Add top5s severe overhead fixups
5285bc5 2026-07-07 Add top5t severe overhead fixups
41661b3 2026-07-08 Add top5u severe overhead fixup
1ce5e25 2026-07-08 Add top5v severe overhead fixups
34b9141 2026-07-08 Add top5w severe overhead fixups
7b35be7 2026-07-08 Add top5x severe overhead fixup
491bf03 2026-07-08 Add top5y severe overhead fixup
5ca30b0 2026-07-08 Add top5z severe overhead fixups
d167712 2026-07-08 Add top5ab severe overhead fixup
97d9895 2026-07-08 Add top5ac severe overhead fixup
e902500 2026-07-08 Add top5af severe overhead fixups
a41dfe0 2026-07-08 Add top5ag severe overhead fixups
210d990 2026-07-08 Add top5ah severe overhead fixups
9ac4bea 2026-07-08 Add top5ai severe overhead fixups
857f49e 2026-07-08 Add top5al severe overhead fixup
542eb67 2026-07-08 Add top5am severe overhead fixup
c8e63c0 2026-07-08 Add top5an severe overhead fixup
b39b0fb 2026-07-08 Add top5ao severe overhead fixup
769605c 2026-07-08 Add top5ap severe overhead fixups
17f550e 2026-07-08 Add top5ar severe overhead fixup
0594587 2026-07-08 Add top5as severe overhead fixups
e984c51 2026-07-08 Add top5at severe overhead fixup
cc38f57 2026-07-08 Add top5au severe overhead fixup
fe566c8 2026-07-08 Add top5av severe overhead fixups
5fd4e09 2026-07-08 Add top5aw severe overhead fixups
81ed748 2026-07-08 Add top5ax severe overhead fixup
0c0039e 2026-07-08 Add top5ay severe overhead fixup
0affe73 2026-07-08 Add top5az severe overhead fixups
```

Any future change to an item in the frozen list requires a new profile ID,
new fingerprints/golden vectors, and an explicit migration note.  Reusing an
existing ID with changed equations is a compatibility failure.
