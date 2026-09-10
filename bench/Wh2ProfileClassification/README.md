# One-comparison host-profile classification

Tracked in `wirehair-sxvz.16.1.20.84.3.5`. This is an isolated compile and
correctness experiment, not a codec promotion or timing result.

The subsequent [full-lifecycle cost screen](../Wh2ProfileClassificationCostR0.md)
did not retain this candidate: overall CONTROL_FAIL, with a separately
controlled reverse-load REGRESSION and no qualifying certified-path benefit.
The compile/correctness evidence below remains valid, but is not a speed win.

After `IsSupportedProfileId` accepts the ID, the supported set contains exactly
the certified profile and K3/K5/K8. In `ValidateHostProfile` only, replace its
second three-ID small-family predicate with `profile_id != CERTIFIED`. The
earlier membership and platform checks, dimensions, small-shape validation,
seed attempts and their error priority stay unchanged. Unsupported and retired
IDs must still fail before this predicate. No extra state or other edits.

The baseline installed native object still contains redundant K5/K3/K8 tests
after dimension validation on the certified-valid path. That is a concrete
instruction-removal hypothesis, not proof of the observed lifecycle regression's
cause. The separately rejected four branch hints, cached block count, division
removal, copied alignment and buffer-release variants are not included.

Configure/build in a fresh external directory with GNU 13.3. The generated
translation unit changes exactly one uniquely anchored predicate in the
hash-pinned source. The baseline must reproduce the preserved K8 object
byte-for-byte. Inspect the candidate assembly before proceeding: no useful
work elimination means stop. Do not run or rebind any historical verifier.

A surviving compiler mechanism requires native, portable-arithmetic and full
sanitizer correctness checks, all installed small-profile contracts and
certified-byte parity before a separately frozen lifecycle timing experiment.
Neither this target nor its build launches timing or a recovery cohort.

## Compile and neutral result

The one-predicate candidate passes the compiler mechanism screen. The native
baseline reproduces the installed object exactly. `ValidateHostProfile` shrinks
from 386 to 368 bytes: the certified ID is checked first and the post-dimension
small-family classification uses one certified-ID comparison instead of the
three small-ID tests. However, serialization now inlines the validator and
grows from 670 to 1,038 bytes. Other function locations also move. These are
both potential benefits and costs to measure, not evidence of net speed.

Final native GNU 13.3/Linux qualification is preserved at
`/tmp/wh2-profile-classification-qualified.lu7R3FrP`. Configure it with
`WH2_CLASSIFICATION_NEUTRAL=ON` only after the compile screen. All 19 CTests
pass: K3/K5/K8 lifecycle/ownership/OOM/conflict/alias checks and the existing
full public-profile suite in native, portable-arithmetic and ASan/UBSan modes;
baseline and candidate metadata probes in all three modes; plus native
certified-byte parity. ASan/UBSan instruments both the replacement object and
the underlying retained archive, not just the test driver. Leak and fake-stack
detection are enabled. Source reviews by the main worker and an independent
bug worker found no remaining issue after repeated passes.

Each metadata probe checks 66,928 host cases and 66,820 wire cases, including
all four supported IDs, all 256 attempt values, 65 geometry slots (some
deliberately repeat), generic and small-profile dimension boundaries, reserved
bytes, unknown/tombstoned/prototype IDs and combined-error precedence. Wire
inputs are independently packed; a different floor-based dimension oracle
checks status, complete guarded outputs, failed-parse zeroing and unchanged
inputs. These operations do not create codecs or allocate large messages.

The link maps extract the original archive's Profile object only for the three
metadata baseline arms. All other targets use the selected replacement (or
the byte-identical reconstructed baseline for its parity emitter). The two
48-case certified streams both contain exactly 2,180,292 bytes and match the
preserved pre-admission digest. Both emitters exit zero with empty stderr.

| Artifact | SHA-256 |
| --- | --- |
| Original native Profile object | `229da2af98169c140774fff3bc970a23b6688bba1586ae49fae862554aa64071` |
| Exact generated source | `745707ed4e5337965e842a5a17d3541458beb9bfc5d0a5c92d14c4ff4093edb3` |
| Native candidate object | `0abb85d6e4092c5dbbe6a47bc8ae42f44188bc3abb033e631ef5a09e6784f146` |
| Portable candidate object | `546e43a53b0a8438a58babf18a15168434a1b6c3c4212bd000ecb6cfe3b438a8` |
| ASan/UBSan candidate object | `85a63c76ec068d595ccb17fa429367e6ef23930a4696810ee21f112af973205a` |
| Final CTest log | `bcdf6c1110b26286f5d8b345b637d238d62602a5c1853f8e8825a3df6eb86eba` |
| Either complete certified byte stream | `2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b` |

Production source, defaults, equations and qualified installed libraries are
unchanged. No performance screen or recovery sample was launched. This neutral
result permits preparing a separately frozen current-versus-candidate lifecycle
screen, including certified, WH1 and K3/K5/K8 retention controls. It does not
restore preserved-path speed, establish strict WH1 superiority or qualify
promotion. Preserve these build directories; do not run CTest discovery or
reconfigure them after treating these log/object hashes as retained evidence.
