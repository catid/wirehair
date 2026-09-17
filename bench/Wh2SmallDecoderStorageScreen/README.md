# Coallocated small-decoder lifecycle screen

Prospective experiment for `wirehair-sxvz.16.1.20.85.1`. Protocol was recorded
in Beads before this observer was implemented. Production is unchanged.
No timing result is claimed by this directory's existence or neutral tests.

The sole candidate from `bench/Wh2SmallDecoderStorage` combines the small
decoder object and packet slab in one allocation. Public allocations fall
from three to two, private core allocations from two to one. This is not
the rejected dormant-core candidate. Equations, encoders, public facades,
certified cores and defaults are unchanged. The naturally different private
slab alignment is part of the measured cost; it is not adjusted or tuned.

## Immutable codec inputs

Use existing qualified libraries only; this observer never compiles codecs:

`/tmp/wh2-small-decoder-storage-final.6dWpqN24/`

- Native baseline SHA256: `bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`
- Native candidate SHA256: `15d27877bd0afa7178f1f410f13acf1114f45b44eec4f512b388eaa0b434b237`
- Qualified manifest SHA256: `8b619caafe0c14366bfde96ce371740d19cd7adf8fb4e56fe8e23d55acbf1b2e`

Do not rerun builds/tests into those directories. Observer-only builds use
`/tmp/wh2-small-decoder-storage-screen.r5PFefBn/{native,portable,asan}` and
load corresponding already-qualified DSOs. Native and portable observations
use the same non-ISA-specialized observer source; ASan/UBSan instruments both
observer and its matching DSOs, with leak and fake-stack checks enabled.
Sanitized workers reject timing modes.

`Generate.py` pins and transforms two previously reviewed observer components
into a new external `Worker.cpp`. It adds explicit WHK3 handle/descriptor/API
dispatch, ordinary/certified WHV2 routing, matched descriptor output, and the
new roster. A WHK3 descriptor is checked against independently constructed
literal bytes, never reinterpreted as WHV2. All public function addresses
must resolve to the requested DSO. Both DSO loading orders are checked.

## Frozen workload and decision

Seven routes, in order: ordinary WHV2 K3; explicit small WHV2 K3, K5, K8;
explicit certified WHV2 K3, K16; separate WHK3 API. Each route includes
B2/full, B2/tail1, B64/full, B1280/full, B1280/tail1, with independent and
borrowed source policies: 70 fixtures and 210 metric cells.

Metrics are complete encoder construction + 3K packets + destruction,
complete low-repair decoder lifecycle, and complete distant-repair decoder
lifecycle. Encoders emit K systematic, K low and K distant packets. Each
decoder uses its own verified first-success endpoint, recovers, then frees.
WH1 uses real ownership-matched APIs and its own packets/endpoint. No
prebuilt codec or preparation is omitted from a measured lifecycle.

Per DSO load order: CPU50 singleton, batch64, 12 replicates, two observation
orders, five comparisons (BB, CC, WH1/WH1, CB, CWH1), 18 positions per panel
(two retained warmups followed by eight adjacent pairs). All 453,600 rows
are retained, including 50,400 warmups. Both load orders total 907,200 rows.
Public source, output, descriptor input/output and both repair-array addresses
must match within every cell/replicate. Output guards, padding, complete
packets and recovered bytes are verified after every WORK callback.

Paired log ratios are averaged across eight pairs within each replicate.
Two-sided t11 95% intervals use 12 replicate means, separately per
cell/comparison/observation order. Ratios are candidate time / reference
time; lower is better. Each DSO load order must independently satisfy:

- All 1,260 same-code intervals strictly inside `[1/1.02,1.02]`.
- All 420 candidate/baseline upper bounds strictly below 1.02.
- All 300 candidate/WH1 upper bounds for the 50 affected small fixtures
  (routes 0,1,2,3,6), three metrics and two orders strictly below 1.
  The other 120 certified WH1 bounds are reported but not required; their
  pre-existing deficits remain unresolved.
- Four decoder-only primary bounds (low/distant x two orders) strictly below
  1. For each, average the 50 affected fixture log ratios *within each
  replicate* before computing the interval, preserving covariance. Encoder
  retention is required; encoder improvement is not a decoder primary.

Any failed same-code gate yields CONTROL_FAIL, overriding nominal wins.
Otherwise failed retention, required WH1 or decoder primary yields FAIL.
Both loading orders must pass; no pooling, subsets or reweighting. Mean
improvement is not a claim of improvement in every cell or simultaneous
confidence across intervals. No all-K or recovery-rate qualification follows
from this bounded screen.

## Bounded one-shot execution

The only scientific namespace is
`/var/tmp/wh2-small-decoder-storage-screen-r0`.
Existence permanently spends it; never retry, rename, trim, tune or rescue.
Complete both orders after statistical failure; stop after infrastructure,
payload or identity failure. No timing launch before neutral qualification,
independent review, manifest sealing, and committing/pushing exact sources.

`Launch.py` is the actual pinned launcher, not a documentary wrapper. It
uses a fresh process group, explicit PATH/LANG/LC_ALL/TZ environment,
CPU250s, address-space384MiB, core0 and file-size128MiB limits. A selector
drains both pipes with stdout128MiB, stderr64KiB and wall270s limits. On
timeout/cap/error it kills the worker process group, reaps the direct worker,
and preserves bounded prefixes. The worker separately enforces elapsed240s
including preparation, checks and final publication/flush; summed measured
WORK, including warmups, must be below240s.
Nonzero exit, stderr, missing rows, invalid endpoints or address mismatch
invalidates the run. No scientific retry occurs on failure.

`Screen.py` checks the neutral manifests, exact DSOs, generated worker,
observer sources/binary/build inputs/test logs and committed own files before
claiming the namespace. It pins these and repository codec sources in
`claim.json`, rechecks all pins before/after each worker, retains captures and
analyses, seals bundle hashes, and supports exact reanalysis. Observer seal
excludes Screen.py to avoid self-reference; the committed controller itself
is independently pinned in the launch claim. Provenance is artifact/neutral
consistency, not a complete transitive toolchain/runtime qualification.

Before any HEAD or result-documentation advancement, perform exact replay
and independent complete raw/fixture/statistics/provenance audits. Preserve
all outcomes. Only a passing screen can justify evaluating production adoption.

## Observer checks

Final observer qualification passed both neutral loading orders in each of
native, portable and full ASan/UBSan builds (6/6 CTests), plus 24/24 synthetic
tests on each Python3.12 and Python3.8. Each neutral run exercises all
70 fixtures x three metrics x three APIs x 64 lifecycles. Original codec
manifest entries were rehashed unchanged. Independent source review caught
and fixed a missing elapsed-cap check after final output publication; earlier
observer builds/logs remain as `*-prepublication`, not timing evidence.
Final generated worker SHA256:
`d395c62ecac4819adfc6125aa50f0d7efb997d5a3d0e43e7533e07a3b0608565`.
Native/portable observer binary SHA256:
`b63fcbc9cece6c072b76a3f06e2994d2a6e4ea72714a807b12e06e9a01482ff5`.

Run `python3 -B -m unittest discover -s bench/Wh2SmallDecoderStorageScreen -v`
and the equivalent Python3.8 invocation. Tests cover exact chronology,
warmups, first-success parity, addresses, covariance, thresholds, required
WH1, standalone inclusion, encoder exclusion from primary, spent namespace,
manifests, resource/environment limits and bounded subprocess failure paths.

Configure each fresh observer build with this directory as CMake source,
Ninja, `DSO_DIR` pointing to the matching qualified library directory, and
`SANITIZE=ON` only for ASan. Build its worker and run both neutral CTests.
Seal the final observer inputs/results before the one native timing launch:

```sh
python3 -B bench/Wh2SmallDecoderStorageScreen/Screen.py --run /tmp/wh2-small-decoder-storage-screen.r5PFefBn/native
```
