# Dormant certified core in small-profile handles

Issue `wirehair-fcsn`; performance screening is `wirehair-fcsn.1`.
**Neutral correctness/lifetime qualification passed; no speed qualification or
production promotion.** No timing namespace has been launched for this candidate.

## Change and rationale

The production small-profile facade contains a 240-byte general
`wirehair_v2::Codec` that it never uses. Its constructor initializes three null
owners and a 216-byte profile; its destructor follows the null-owner cleanup
path. The generated candidate replaces only this member's declaration and
lifetime management with a same-layout anonymous union. It placement-constructs
and explicitly destroys the core only for immutable `SmallK == 0`. Copies and
moves are explicitly deleted. No handle is retagged or reuses the union storage.

All nine direct `Impl` expressions are either on newly allocated certified
handles or guarded by `SmallK`. Small encoder/decoder allocation failures still
unwind typed `unique_ptr`s, and descriptor-publication failure still goes through
`FreePublicCodec`. The evaluator/basis destruction order is unchanged. The new
conditional lifetime is the only codec change: no alignment modification,
allocation removal, row/seed/default change, new field, or recovery equation.

This is independent of the rejected aligned-basis candidate. It does not explain
or rescue that candidate's historical timing results. Production source remains
unchanged at SHA256 `975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15`.

## Qualification

Fresh builds are retained under `/tmp/wh2-small-dormant-core.qARBuvHF/`.

| Build | CTests | Scope |
| --- | ---: | --- |
| Native | 21/21 pass | Both arms, real SIMD dispatch |
| Portable | 21/21 pass | Both arms, forced portable arithmetic |
| ASan+UBSan | 21/21 pass | All codec/test objects instrumented; leaks and fake stacks enabled |

Each arm runs the existing complete K3/K5/K8 public-codec suites, profile and
borrowed-source tests, unchanged-basis allocation test, a new lifetime test,
layout comparison and separately hook-enabled facade-OOM tests. Small-codec
suites cover independent/borrowed/legacy constructor routes, transactional OOM,
aliases, partial tails, conflicting/duplicate packets, repeated recovery,
protected borrowed sources, allocation-free detach and C++ moves/replacement.

The reused allocation test checks 936 forced placements, 3,384 OOM boundaries,
and 14,976 packets per arm. Neither arm enables its alignment-candidate macro:
both retain the original basis allocation sizes and copy positions.

The new lifetime executable tests all applicable ordinary/explicit-ID/serialized
constructor families, legacy/independent/borrowed policies, widths 2/64/1280 and
full/one-byte tails. This includes ordinary small K3, explicit/serialized small
K3/K5/K8, ordinary certified K16, and explicit/serialized certified K3/K16.
Ordinary certified K3 is not a route: its ordinary constructor selects small K3.

Per arm/backend it checks 432 successful encoder/decoder lifetimes and 1,179
injected OOM boundaries: all 819 small-profile boundaries plus 360 certified
positive-control failures (only the first two constructor allocations). This
is **not** exhaustive certified OOM coverage. The existing facade-fault tests
provide additional transactional certified controls.

Only the lifetime executable wraps the actual external general-core C1/C2/D1/D2
symbols. The tracked facade address must match the core address, whose offset
zero is separately verified. Small baseline lifetimes execute exactly one
constructor/destructor; small candidate lifetimes execute zero of each.
Certified positive controls execute one of each. Facade-allocation failure
executes neither; later tested OOM failures clean up the appropriate active core.

For small-profile handles, the candidate lifetime allocator fills the first
240 facade bytes with `0xa5`
**before outer construction**, and under ASan additionally user-poisons that
storage. It remains poisoned through all API operations and full destruction,
including OOM unwinds. Replacement delete alone unpoisons/checks the canary.
Native canaries detect writes, not reads; the fully instrumented sanitizer run
adds read detection. No poison or wrapping is linked into the actual candidate
DSOs or static archives.

Layout tests construct separate original/current models (never reinterpret a
live public handle). Core size 240, certified size 272, small size 296, all
alignments, every header/owner offset and noncopyability match. Two DSO load-order
neutral parity runs cover the retained 24-fixture explicit-small worker. Certified
compatibility independently checks both producers' successful exits, empty stderr
and all 2,180,292 bytes over 48 cases before comparing streams. Every backend
matches SHA256 `2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.
These are compatibility checks, not new recovery-rate samples or all-K evidence.

Source review caught and fixed a test-only descriptor assumption before runs:
the test now retains the successful encoder's selected descriptor, never assumes
certified attempt zero. The first native lifetime probes failed a 256-entry
registry bound: total allocation calls were incorrectly used as live-slot
indices. Certified construction reaches 268 calls while freeing temporaries.
The fixed registry reuses live slots and keeps failure ordinals separate. Small
OOM enumeration remains exhaustive; certified enumeration is explicitly bounded
to the two lifetime transitions above. The failed log is retained as
`native-initial-lifetime-bound-failure.log`. Candidate source/DSO did not change
while correcting the probe. Repeated final source passes found no remaining bugs.

## Actual machine-code effects

Native `CreateEncoderForProfile` and `wirehair_v2_decoder_create` each go from
four general-core constructor callsites to one; only the certified call remains.
Small initialization stores begin at offset `0xf0`, with tag `0x101` and owner
offsets `0x110/0x118/0x120` preserved. No store initializes the dormant 240 bytes.

Destructor callsites still exist behind `SmallK == 0` guards. Valid small tags
skip them, as the independent call-count probe verifies. These guards add a tag
load/branch: `FreeSmallCodec` grows 292 to 324 bytes and each small exception
cleanup helper grows 135 to 166 bytes. Certified header initialization now
precedes its core constructor, and certified error cleanup gains tag guards.
Ordinary encoder code grows 1,386 to 1,441 bytes; explicit-ID code grows 1,910
to 1,947. Total ELF `.text` shrinks 80 bytes, but function placement changes.

Public encode/decode/recover/detach/free, serialized constructors, and the small
encode/decode/recover helpers retain normalized instructions. This is not a
claim of identical addresses, equal timing, or zero retained-path regression.
The candidate removes real executed work but does not merely delete machine
code. Actual lifecycle and certified retention gates remain necessary.

Independent linked-data inspection finds 152 changed 32-bit entries across 12
code-relative jump tables. Each resolves to the same function and intra-function
offset in the two libraries; all other `.rodata` bytes, including coefficient
tables, match. Two relative relocation addends change only to follow relocated
initialization/finalization functions. Do not carry forward the older alignment
experiment's whole-`.rodata`/relocation-payload identity claim to this candidate.

## Qualified inputs and next gate

Generated profile source, identical across all three modes:
`9c8ccd1c8ae40b3d075061a9e8237c2ad431e46339d6ff047cd9f257bb3a4a59`.

| Mode | Candidate DSO SHA256 |
| --- | --- |
| Native | `f9c2e226f180fb885561e65b0556df5b39c43e49bb53141ba9f4724cc5a8f9be` |
| Portable | `1b37aaece6f39fff9b0c14f790dde8bcc177ea37659d2207c482bf3027d0eb6f` |
| Sanitized | `44001f634412d3bcd0bd7aedc66bb1bb8606c053c4cf9c53ee95661f2aa7d1bb` |

Native baseline is byte-identical to the retained production baseline:
`bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`.
Do not rebuild these qualified libraries, generated sources or test logs.
`NEUTRAL.sha256` in the retained build root pins sources, generated candidates,
all six DSOs and archives, lifetime/layout executables, compile/test definitions,
final test logs, compatibility streams and the initial failed log. Its SHA256 is
`dc84b6d96b4d5994c717a339a11326229a9108a860d730e61a7aca5900efa4f4`.

Next: freeze a short actual-lifecycle screen using these exact libraries, with
ordinary K3, explicit small routes, both source policies, full/partial payloads,
preserved certified K3/K16, both load orders, A/A controls and separate WH1
comparisons. Primary wins, retention bounds, roster, statistics and stop rule
must be fixed before timing. No production adoption, broad speed claim or all-K
campaign is justified until that screen passes. Existing spent namespaces remain
closed; this candidate has not been timed.

To reproduce correctness, use a **different fresh external** GNU/Linux build:
`cmake -S bench/Wh2SmallDormantCore -B <fresh-dir> -G Ninja`, then build and run
CTest. Use `-DPORTABLE=ON` or `-DSANITIZE=ON` for those respective modes. Test-only
layout/fault objects and linker wrapping must stay out of any performance DSO.
