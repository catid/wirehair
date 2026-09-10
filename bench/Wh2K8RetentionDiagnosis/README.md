# Retained ordinary-K8 shared-retention diagnosis

Read-only follow-up to the permanently spent
`wirehair.wh2.k8-ordinary-shared-retention-r0` gate, tracked in
`wirehair-sxvz.16.1.20.82.5.3.4.6`. The original result remains
**CONTROL_FAIL**; the reverse library load independently passes its controls
and reports seventeen resolved slowdown cells. Retention is unqualified,
gates B/C remain unlaunched, and ordinary production K8 still selects certified
WH2. This diagnosis changes no codec, equation, measurement or decision.

## Scope and reproducibility

`Analyze.py` uses only the Python standard library. It authenticates the sealed
independent `POSTRUN.json` audit anchor, hashes all eight scientific bundle
members, streams both complete raw files, and recomputes the original
per-replicate paired logs exactly. It does not import the scientific controller,
invoke an old strict current-HEAD verifier, execute a codec, or recalculate a
qualification verdict.

The complete report visits 196,992 records and 912 statistical cells. All
preludes, observations and cells remain present; the 432 expanded records are
the twelve complete eighteen-position panels for each of the two failed
controls, not a selected replacement sample. Each cell includes descriptive
WORK/thread/wall/wait/gap distributions, fault/switch deltas and complete
ordered public-handle address-vector identities. Original statistics are
copied unchanged after their paired-log equality check.

Neutral tests, with no codec execution:

```sh
python3 -B -m unittest discover -s bench/Wh2K8RetentionDiagnosis -v
python3.8 -B -m unittest discover -s bench/Wh2K8RetentionDiagnosis -v
```

Both interpreters pass all eight tests. Read-only report regeneration is
permitted with a fresh external output path, for example:

```sh
python3 -B bench/Wh2K8RetentionDiagnosis/Analyze.py /absolute/fresh/external/report.json
```

The parent directory must already exist. The output is created exclusively,
flushed and sealed mode 0400; an existing output is never overwritten. The
script rejects repository and source-bundle output paths. Do not put a new
report in any other spent scientific namespace either. Historical controller
and independent-auditor strict current-HEAD checks must not be rerun now that
their producing HEAD and pinned documents have advanced.

## Findings

### Added K8 checks do not execute on valid K3 constructor paths

Independent disassembly review compared the exact measured DSOs. For full K3
messages with two- and 64-byte blocks, each bare ordinary constructor body
executes 48 non-NOP instructions, counting its call/return but excluding
callees and options-wrapper work. Each has three comparisons, four conditional
branches and one call, with the same `0xc8` stack frame. These are equal
counts, not identical instructions. The baseline has one additional NOP. The
candidate skips every added K8 shape check.

The constructor remains at DSO offset `0x4ccd0` but grows from 1,386 to 1,485
bytes. Its K3 branch direction and code placement change. The borrowed-options
wrapper moves by `0x60`; its 1,282-byte body changes only the relative
displacement of the existing ordinary-constructor call. Downstream K3
construction, the lookup, serializer and relevant PLT entries have identical
bytes at identical offsets. Both independent and borrowed K3 use this same
bare constructor path; no borrowed-specific reduction in executed work was
identified.

WH1 create/encode/free symbols also retain identical offsets, sizes and bytes.
These facts rule out added executed K8 checks as the K3 explanation, but not
physical effects of branch direction, placement, internal memory or machine
state. They establish neither a timing cause nor a new optimization.

### Public handle placement does not explain the selected encoder differences

For all 32 encoder treatment cells covering WH1 K3 B2/B1280 and ordinary
public K3 B2/B64/B1280, available source policies, both measurement orders and
both load orders, every one of 108 paired complete address vectors matches
across logical sides. Different public-handle alignment is not the explanation
for these cells. Internal allocation addresses and physical pages are not
recorded, so equality of those cannot be inferred.

31 of these 32 cells have zero capture fault/switch deltas, including every
WH1 B2 and ordinary K3 cell. The exception is current-first WH1 B1280,
key `[17,0,2,0]`: candidate-side prelude record 46405 has one involuntary
switch. It remains in the report. Inter-capture counters are a separate scope
and are not uniformly zero.

Thread-bracket medians move with WORK medians in the investigated encoder
cells. The original control-passing reverse-load borrowed K3 encoder estimates
remain 2.008-2.066% slower at B2 and 2.033-2.304% slower at B64. WH1 B2
remains 2.835-3.044% slower. Descriptive medians are not new speed estimates.

### Failed-control excursions resemble the prior boundary association

The largest positive-log panels of the two failed controls include these
retained records:

| Cell key | Replicate | Record | WORK ns | Thread bracket ns | Start modulo 1 ms, ns |
| --- | ---: | ---: | ---: | ---: | ---: |
| `[13,1,0,0]` | 2 | 17643 | 11167 | 11267 | 995314 |
| `[13,1,0,0]` | 2 | 17653 | 10496 | 10606 | 995034 |
| `[21,1,1,0]` | 6 | 50877 | 20440 | 20571 | 988714 |

All three WORK intervals cross an absolute millisecond boundary. Their whole
eighteen-record panels have zero capture and gap counter changes. Each failed
cell has one involuntary switch elsewhere in its inter-capture gaps, not in
these largest panels. This resembles the already documented
[clock-boundary association](../Wh2ClockBoundaryR0.md), not a new cause,
proof of a clock bug, or reason to remove an observation or repeat a gate.

Thread brackets include capture work, and counter brackets are wider than
inner WORK. Gap and wait durations overlap and must not be added. Gap counters
cannot attribute events specifically to waiting. Internal allocations,
physical pages, cache/predictor/frequency state and some interruptions remain
unobserved. The cause of the regressions and control excursions is unresolved.

## Independent review and retained attempts

The independent bugworker read the source/tests, compared the exact ELF
paths, and used a separate streaming standard-library calculation without
importing this analyzer. It checked all 196,992 raw records, every aggregate
in all 912 cells, every copied original-statistics object and all fields of
the 432 focused records. Complete address vectors were compared directly, not
only by hash. The review found no remaining implementation/accounting bug and
corrected the counter interpretation to 31/32 zero-capture cells as above.
Repeated main-agent source-reading passes are clean.

Diagnostic-only development failures remain disclosed: the initial side tuple
omitted two entries and failed synthetic tests; the first actual read rejected
the legitimate 2,681,854-byte raw header at an overly small 2-MiB line bound
and produced no report. The corrected 4-MiB bounded reader has a 3-MiB
synthetic test. Later additions clarify overlapping gap/wait scope, reject
undeclared failed cells, check complete publication, and test both measurement
orders. No codec or scientific timing was rerun during these fixes.

The final reports are byte-identical on Python 3.12 and 3.8, each 3,303,965
bytes, in `/tmp/wh2-k8-retention-diagnosis.SDeYBhpA`:

| Artifact | SHA256 |
| --- | --- |
| `final-python312.json`, `final-python38.json` | `70a0bd9aa322b5c918cd39fc547d3a5637ea432f5f53101c210b3e8948296764` |
| `INDEPENDENT_REVIEW.md` (10,649 bytes, mode 0400) | `868d9a47f8210edddb026d6076eae9d7763e5f3ab080224dbeed2b1b885d74e1` |
| Earlier `report-python312.json`, `report-python38.json` | `15c44f3248da298ab55a8ed4f3c1f9e26f920acda3b66b0366c35ee3fa46a514` |

The fixed audit anchor is
`/tmp/wh2-k8-retention-auditor.jY1sGYtx/POSTRUN.json`, SHA256
`f030a6f8e27005dacdf818fd0d9e402f03f872c298263a6ae94ea5185ff8e4ec`.
It authenticates the original complete scientific file roster; the
[gate-A result](../Wh2K8SharedRetention/README.md#independently-audited-gate-a-result-not-retained)
records those identities and both frozen load-order decisions.

No selector/hint/validator/serializer/alignment/release-order rescue is
justified by this diagnosis. A next structural experiment for an uncovered K
must have its own prospective bounded plan and independent gates. Existing
scoped K3/K5/K6/K8 results remain separate evidence. This report establishes
neither pre-admission restoration nor current-path retention, default
promotion, shared WH1 superiority or the full all-K speed/recovery objective.
