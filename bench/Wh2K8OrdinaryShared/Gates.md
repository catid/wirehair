# Prospective ordinary K8 shared-call qualification

Preparation under `wirehair-sxvz.16.1.20.82.5.3.4.3`. This specifies distinct
questions before implementing or launching their instruments. It is not a
timing result, a claimed scientific namespace, or default promotion. Existing
spent experiments and their decisions remain unchanged.

Use the exact native DSOs qualified by `Prepare.py`: current `ef684aba…` and
six-line-selector candidate `fe7527a7…`. Replace only Profile in its original
nineteen-object position. No serialization, validation, hint, field, ownership,
SONAME, export-map or linker-binding changes accompany the selector.

## Actual APIs and ownership

| Route | Encoder construction |
| --- | --- |
| Ordinary WH2 independent | `wirehair_v2_encoder_create` (bare constructor) |
| Ordinary WH2 borrowed | `wirehair_v2_encoder_create_with_options`, borrowed policy |
| Explicit certified/K5/K8 | `wirehair_v2_encoder_create_profile_id_with_options`, matched policy |
| WH1 independent | `wirehair_encoder_create_owned_ex(nullptr, ...)` |
| WH1 borrowed | `wirehair_encoder_create_ex(nullptr, ...)` |
| Opt-in K3/K6 | Existing `wirehair_small_encoder_create` / `wirehair_k6_encoder_create`, original roster policy |

Explicit controls use actual profile IDs: certified `4b295bbb47f4f9c9`, K5
`80070c81bfe375f1`, K8 `7a9276b85c730ae0`. Never substitute an ordinary
constructor for a certified control. Every WHV2 receiver consumes its own
actual serialized descriptor. WH1 receivers use `wirehair_decoder_create_ex`.
Count all actual preparation/copies in creation and all cleanup in teardown.

Policies share receiver implementations; repeated decoder observations do
not create independent recovery samples. The old ordinary K3 speed screen
compared both WH2 policies with borrowed WH1. Gate C corrects that ownership
scope prospectively; it does not rewrite the old measurements.

## A: candidate-versus-current retention

Reuse the 38-case roster and WORK from the prior shared lifecycle instrument,
not its rejected candidate or spent controller. The cases are full messages:

- Certified borrowed K2/3/4/6/128, each at B2 and B1280: ten cases.
- Certified independent K3 at B2 and B1280: two cases.
- Opt-in K3 at B2/B1280 under both policies: four cases.
- Borrowed WH1 K3 and borrowed opt-in K6 at B2/B1280: four cases.
- Ordinary K3 and explicit K5/K8 at B2/B64/B1280 under both policies:
  eighteen cases.

The encoder lifecycle creates, produces the descriptor, encodes IDs `0..K+7`
and six IDs `UINT32_MAX-2*j`, then frees. The decoder creates, consumes eight
low repairs, six distant repairs, then systematics, stopping at its actual
first success, recovers, and frees. Use batch128 except K128 batch4.

Each case/metric has three pairs: current/current, candidate/candidate,
current/candidate. Analyze both measurement orders in each separate DSO load
order. Across both load orders this is 196,992 callbacks, 608 same-code
controls and 304 retention decisions.

This gate asks for retention, not an incidental certified-profile speedup.
Unlike the rejected serializer experiment, it does not require a certified
metric to improve in all four order cells. It also does not test the historical
K3 workload completely: this encoder emits seventeen K3 packets and the
decoder has one combined stream. Gate C retains the eighteen-packet and
separate-stream workload explicitly.

Reference WORK: `bench/Wh2AdmissionRegressionCostR0.cpp`, `RunWork`.
Reference roster: `bench/Wh2ValidatedSerializationCostR0.md`.

## B: ordinary K8 versus ownership-matched WH1

Four arms: candidate ordinary independent/borrowed, and current-DSO WH1
owned/borrowed. The WH1 provider is fixed to the current DSO before measuring;
it is not chosen after inspecting timings.

Reuse the ordinary K8 six-shape workload: B2/B64/B1280, full and one-byte
tails. Encoder creation/descriptor, 24 packets (IDs0..15 and eight distant
IDs), and free are timed. Low-ID and distant-ID decoder lifecycles are separate;
each supplies eight repairs then eight systematics, stopping at its own first
success, recovering and freeing. Retain any systematic fallback. Batch128.

Four A/A pairs and two ownership-matched superiority pairs produce 93,312
callbacks, 288 controls and 144 superiority decisions across both load orders.
This gate deliberately makes no new shared-speed claim against certified WH2.

Reference WORK: `bench/Wh2K8PublicCostR0.cpp`; ordinary route mapping and
retained workload: `bench/Wh2K8OrdinaryCost/README.md`.

## C: ownership-matched K3 speed and complete K3 retention

Six arms: current ordinary K3, candidate ordinary K3, and current-DSO WH1,
each under independent/owned and borrowed input. Three full-message shapes:
B2/B64/B1280. Reuse the historical K3 WORK rather than gate A's shorter work.

Encoder creation/descriptor, eighteen packets (IDs0..11 and six distant IDs),
and free are timed. Separate low/distant decoder streams supply six repairs
then systematic IDs0..2, at most nine feeds, stopping at the arm's actual first
success, recovering and freeing. Batch128.

Use twelve pairs: six A/A, two current-K3/WH1 superiority, two candidate-K3/WH1
superiority, and two candidate/current-K3 retention. Direct candidate/WH1
comparisons avoid inferring superiority by combining different confidence
intervals. Across both load orders: 93,312 callbacks, 216 controls,
144 superiority decisions and 72 retention decisions.

Reference WORK: `bench/Wh2K3OrdinaryCostR0.cpp`, `RunWork`; replace its old
borrowed-only WH1 construction mapping explicitly for the independent arm.

## Common decisions, bounds and launch conditions

Keep twelve replicates, both measurement orders, eighteen positions, eight
paired contrasts per panel, all 48 completion-relative conditioning phases,
and deferred formatting/publication. CPU50 and its existing exact target
identity are required. Both DSO load orders run separately; never pool them.
Treatment log ratios are `log(candidate_time/current_time)` for retention and
`log(WH2_time/WH1_time)` for superiority. The reused reader computes logical
side1 minus side0, so freeze pair tuples in control-then-treatment order;
measurement order changes execution order, never the ratio's meaning.

- A/A: every 95% log-ratio interval lies strictly inside `+/-log(1.02)`.
- Retention: any `lower95_log > 0` is a resolved slowdown; any
  `upper95_log >= log(1.02)` leaves retention uncertain. Either prevents PASS.
- Superiority: every `upper95_log < 0`; there is no 5% minimum improvement.
- A failed control dominates treatment estimates. Keep every observation;
  no filtering, pooling, timing subtraction, resizing, rerunning, or retuning.

Retention means no resolved slowdown within this instrument's declared
precision. It is not proof of exact zero regression or strict WH1 superiority.

Per gate: worker WORK150s/CPU210s/wall210s/address512MiB; observer240s per
load, controller600s; raw192MiB per load, stderr64KiB, complete bundle448MiB.
The complete proposed roster is 383,616 callbacks across three distinct gates.
Run A first. If A does not pass, stop this timing sequence before B or C;
retain its terminal result and diagnose before proposing different work.

Before any launch, implement and independently qualify the actual shared
adapter, complete arm/pair rosters, descriptor/payload/rank fixtures, ledger,
late-call cleanup/exception/clock failures, and publication failures. Freeze
the complete source/build/runtime/target/command inputs and exact gate-specific
contract, commit/push, then claim its one-shot namespace. This document alone
does not authorize an unqualified instrument. Authenticate all 53 exports,
six internal public GOT targets, 37 runtime targets and separate GF contexts
before and after calls in both load orders. No static codec fallback.

## Separate unresolved requirement: pre-admission restoration

None of A/B/C compares against the pre-admission DSO or repairs its recorded
regressions. Preserved explicit-certified constructors do not enter the
ordinary selector; ordinary K3 returns before its new K8 branch. Consequently
there is no direct reduced-work restoration mechanism in this six-line diff.
Do not launch a selector-only rerun of old restoration workloads or revive
rejected hint/validator/serializer/division/alignment/release-order families.

Even passing A/B/C would establish only their stated shared-call scopes.
Pre-admission restoration, broader K/width/hardware and construction-seed
validation, installed K5 speed, and the full WH2-versus-WH1 objective remain
separate requirements. No default promotion follows from preparation alone.
