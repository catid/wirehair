# Validated-local WHV2 serialization candidate

This is a single benchmark-only extraction of the descriptor byte writer.
Production source, profiles, defaults and qualified archives are not modified.
There is no timing entry point or performance result here.

## Mechanism and scope

At production source SHA-256
`975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15`,
the descriptor-producing constructors call the public serializer after their
private profile has already been validated. On successful paths:

| Constructor | Host validations before candidate | After candidate |
| --- | ---: | ---: |
| Ordinary certified | 2 | 1 |
| Explicit certified profile ID | 3 | 2 |
| Small K3/K5/K8 profile ID (including ordinary K3 dispatch) | 2 | 1 |

These counts concern facade validation, not validation inside the codec.
Options wrappers inherit their selected constructor's count. Each removed
call repeats dimension validation, including division, and public serializer
alias/capacity checks on distinct, fixed-size stack locals.

`Overlay.cmake` extracts the existing little-endian byte stores into
`SerializeValidatedProfile`, a private helper that writes all 32 bytes,
including zero reserved bytes. Public serialization retains every preflight,
its distinct temporary, and its final copy, so exact/partial host-output alias
support and error precedence remain unchanged. The three constructors use the
helper only after their prior validations succeed. Earlier cleanup, the
explicit-ID identity check, and descriptor/handle publication order remain.
`MakePublicProfile` still expands and compares the canonical equation profile;
that safety check is not redundant and is not removed.

The decoder never serializes a descriptor. This candidate cannot directly
remove decoder work or explain unchanged-WH1 timing shifts in earlier screens.
The separate duplicate parse in the serialized-profile options constructor is
not changed. Rejected validator predicates, branch hints, division-removal,
copied alignment and buffer-release variants are not part of this candidate.

The three constructor calls no longer pass through the exported serializer,
so intentionally interposing that function no longer intercepts construction.
The public API does not promise this internal call composition. A future DSO
cost gate must nevertheless authenticate the changed binding graph explicitly;
the old six-slot internal GOT roster must not be blindly reused.

## Neutral build

Use a fresh external build directory, native GNU 13.3/Linux, and
`-DWH2_VALIDATED_SERIALIZATION_NEUTRAL=ON`. The default is compile-only.
The baseline Profile object must reproduce the qualified K8 object exactly.
The native, portable-arithmetic and fully ASan/UBSan candidate objects link
against the existing, hash-checked K8 archives without rebuilding them.

The checks reuse the existing K3/K5/K8 public suites, full profile suite,
borrowed-source and borrowed-fault suites, independent host/wire metadata
oracle, and the 48-case certified byte-stream comparison. The previous
metadata and parity harnesses are reused read-only, not their spent timing
instruments or receipts. Full lifecycle timing, source-policy/old-path
retention and actual WH1 comparisons require separate prospective gates.

The borrowed-fault suite specifically uses separate baseline/candidate facade
objects with `WIREHAIR_TESTING=1`, enabling its allocation-failure hook. All
other targets use production-style facade objects. An initial incorrect build
omitted that hook and failed the three borrowed-fault tests (22/25 passed);
the disabled injection also left the test's unexpectedly successful handles
unfreed. That build and its failure/leak logs are retained at
`/tmp/wh2-validated-serialization-neutral.yGhFOWxg`, not qualification evidence.
The correction changes only test wiring, not the candidate source transform.

## Neutral qualification result

The corrected build at
`/tmp/wh2-validated-serialization-qualified.B0LfQJ69` passes all 28 checks:
K3/K5/K8, profile, borrowed-source, separate baseline/candidate borrowed-fault,
and baseline/candidate metadata suites in native, portable and sanitizer modes,
plus native certified byte-stream parity. ASan/UBSan checks enable leak and
fake-stack detection; no sanitizer report remains in this qualified run.
Each metadata arm covers 66,928 host and 66,820 wire cases. Both certified
emitters exit successfully with empty stderr and exactly 2,180,292 identical
bytes, matching the retained 48-case hash
`2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.

The native baseline object reproduces byte-for-byte. The candidate's generated
source hash is
`b51ee03b2e7454d561d6e42b6abeb10b2038812b5acfeb9bb5c26ab09fe73ac6`;
its production-style native object hash is
`5291f174f65d974712f92175434b90a2b0039b88ccc947fa848871c66c4284b1`.
The portable object is
`018ef15722afc645878bc80079a8dac27bc284318d8fa6f39f97d49dc3ed3a07`,
and the fully instrumented object is
`a57ac11fc340526a259a7256960c2e786335c0791c29daff0526e5f05b8f2e5e`.
The qualification `Testing/Temporary/LastTest.log` hash is
`533212a82a14ac5750439b5764b62e639d9b6b933732cae031c4c873b4b6a02f`.
Do not run test discovery or reconfigure in this archived directory: either
can overwrite evidence. Use a fresh build for later checks.

Native assembly confirms all three constructor relocations to exported
`wirehair_v2_profile_serialize` are removed. All existing public definitions
remain. The byte writer is one 56-byte internal helper, called from the public
serializer and the three constructor sites. The compiler also creates a
partial validator clone; validation counts alone do not describe every
code-layout effect. Observed function sizes, excluding cold fragments:

| Function | Baseline bytes | Candidate bytes |
| --- | ---: | ---: |
| Public serializer | 670 | 446 |
| Ordinary encoder constructor | 1,386 | 1,306 |
| Explicit profile-ID constructor | 1,910 | 1,517 |

GNU `size` reports total text 38,136 to 37,768 bytes for the Profile object.
The canonical `MakePublicProfile` source and all code from the serialized
options constructor through decode/recover/free remain textually unchanged.
These are correctness and static-mechanism results, not elapsed-time evidence,
restored preserved-path performance, improved recovery, or promotion.
