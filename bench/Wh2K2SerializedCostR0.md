# K2 serialized lifecycle cost R0

Protocol `wirehair.wh2.k2-serialized-cost-r0`; sole scientific namespace
`/var/tmp/wh2-k2-serialized-cost-r0`. Building and neutral testing do not launch
a timing cohort. This benchmark does not install or promote K2.

## Frozen scope

Six actual C API arms compare ordinary WH2, the qualified WHK2 boundary and
WH1, separately for independent/owned and borrowed input. Ordinary independent
WH2 uses `wirehair_v2_encoder_create`; borrowed uses its options constructor.
At K2 both must still select the certified profile. WH1 uses its actual owned
or borrowed constructor. K2 retains the sealed `(2,3)` / `(3,3)` equations.

Shapes are `(block bytes, tail bytes)`:
`(2,2), (2,1), (64,64), (64,1), (1280,1280), (1280,1)`.
Message length is block plus tail. All arms emit the same 14 IDs: 0 through 7,
then `UINT32_MAX - 2*j` for j=0..5. ID1 emits only the meaningful tail; all
repair packets remain full width. Unwritten packet-slot bytes are guarded,
not counted as emitted data.

Encoder work is create, descriptor, all 14 packets and free. Decoder work is
create, every feed through its own first success, recover and free. Low-ID
receivers see IDs2..7 then0,1; distant receivers see the six distant IDs then0,1.
First-success counts are frozen from the neutral fixtures and independently
checked against generator-row rank. Every timed status is checked; use of
systematic fallback is recorded. Neither packet counts nor recovery endpoints
are equalized between arms.

Reuse the deferred-output K5 instrument: 12 replicates, both orders, ten pairs
(six same-code controls and four candidate comparisons), three metrics, six
shapes, 18 positions and batch128: 77,760 callbacks. Each panel retains its
two prelude observations and eight paired contrasts. All 48 relative phase
bins remain covered. All 216 control confidence intervals must fit strictly
inside +/-log(1.02); all 144 candidate upper95 time-ratio bounds must be below
one. No 5% minimum, filtering, subtraction, pooling, resizing or rescue.

Worker caps: CPU240s, wall300s, timed work180s, address space384MiB.
Observer330s; total controller480s; raw192MiB, stderr64KiB, bundle256MiB.
Native CPU50 and its existing 617-byte identity/dispatch are fixed. Portable
and sanitizer binaries reject scientific execution. Result formatting/output
occurs only after measurement; failed and incomplete publication is rejected.

## Qualification and provenance

Reuse the exact current-production-equivalent archives retained at
`/tmp/wh2-v2-k5-admission.VR7a0QAY`, with authenticated producer records from
the installed K5 recovery qualification. Preserve original source bytes when
K2's compile-time guards differ; do not relax or rerun a historical verifier.
The retained native object-identity proof must match the actual archive's
three SmallCore users. Rebuild only the benchmark boundary in a fresh directory
against the selected GF runtime; ASAN direct-table users must match its
`-march=native` private context layout. Run all seven existing boundary checks.

Before science: qualify 216 real WORK cases, the complete roster, independent
polynomial packet and rank oracles, all meaningful tail lengths/guards,
ownership lifetime, late errors/exceptions/clock/source failures and output
failures. Require backend fixture agreement and synthetic gate mutation tests.
Commit and pin qualified source before the sole launch, and audit the retained
result before advancing any pinned source or documentation.

Passing this prototype gate alone does not prove installed K2 speed, a new
recovery rate, preserved-path regression repair, all-K superiority, shared-call
or cold-start performance, or real non-GFNI-host performance.
