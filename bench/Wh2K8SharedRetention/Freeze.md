# Ordinary K8 shared retention gate A: frozen controller

Protocol `wirehair.wh2.k8-ordinary-shared-retention-r0`; sole namespace
`/var/tmp/wh2-k8-ordinary-shared-retention-r0`. Tracked in
`wirehair-sxvz.16.1.20.82.5.3.4.5`. The complete prospective WORK, roster,
ownership, order, target, statistics and caps are gate A plus the common rules
in [Gates.md](../Wh2K8OrdinaryShared/Gates.md). This file does not change them.

`Runner.py` passes the exact new Adapter configuration to the generic closed
builder, receipt verifier, bounded controller and full replay. Every such call
has explicit configuration; no old default, historical receipt/run/replay,
candidate or scientific namespace is selected. Both accepted native and
ASan-driver observer manifests and all their inputs/artifacts become mandatory
new build inputs. The Adapter, worker derivation and codec DSOs are unchanged.
The Runner, its tests and this contract are also mandatory committed inputs.

Native build flags are unchanged. ASan-driver builds remain neutral-only and
load native DSOs. Native receipts reject instrumented workers, missing pins,
source drift, incorrect protocol, wrong provider and changed prerequisites.
Before any scientific launch, independently review exact source/build/receipt,
qualify the new controller with synthetic-only tests, and commit/push its exact
source. The separately prepared independent auditor must also be frozen before
launch. No tuning decision may use the new outcome before that freeze.

The native worker performs 98,496 callbacks per load order, 304 same-code
controls and 152 retention decisions. Keep both load orders separate and all
196,992 observations. Each run retains actual packet/descriptor/recovery
validation, own-first-success endpoints, complete attempted-call ledgers,
addresses, clock/fault records and all 48 delay phases. No timing subtraction,
pooling, filtering, sample resizing or rerun. Formatting/output are deferred.

The controller always attempts both orders after an ordinary per-order capture
or validation failure, while preserving the overall resource deadline. It
seals raw data, stderr, process records, analysis and the original claim. An
infrastructure/validation failure gives INVALID, never PASS. Otherwise outcome
precedence is CONTROL_FAIL, REGRESSION, INCONCLUSIVE, PASS. Gate A requires no
incidental speedup. There is no five-percent minimum improvement. Any resolved
slowdown prevents retention; a failed same-code precision control prevents
treatment qualification. An already-existing namespace cannot be overwritten.

WORK150s/CPU210s/wall210s/address512MiB per worker; observer240s per load,
controller600s; raw192MiB per load, stderr64KiB and bundle448MiB. Native
library bindings, source/producer/build/runtime pins and target identity must
remain exact. The controller cleans the allocator/loader environment using the
existing explicit allowlist. No allocator or loader policy is tuned.

Use fresh external mode-named build children:

```sh
python3 bench/Wh2K8SharedRetention/Runner.py build native /absolute/fresh/parent/native
env ASAN_OPTIONS=detect_leaks=1:detect_stack_use_after_return=1 UBSAN_OPTIONS=halt_on_error=1 python3 bench/Wh2K8SharedRetention/Runner.py build asan-driver /absolute/fresh/parent/asan-driver
python3 -m unittest discover -s bench/Wh2K8SharedRetention -p 'test_*.py' -v
python3 bench/Wh2K8SharedRetention/Runner.py receipt /absolute/fresh/parent/native /absolute/fresh/parent/receipt.json
```

Only after all prelaunch requirements pass, invoke `Runner.py run` once with
that exact receipt. After terminal completion, exact `Runner.py replay` and
the independent full raw/fixture/statistical/provenance audit must both finish
before producing HEAD or any pinned source advances. If A does not pass, stop
the B/C timing sequence and retain its complete outcome. No default promotion,
pre-admission restoration, new recovery rate, full historical K3 retention,
ordinary K8-versus-WH1 shared speed or all-K guarantee follows from gate A.
