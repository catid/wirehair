# Small-decoder unit-diagonal experiment

Issue: `wirehair-sxvz.16.1.20.86`. Experiment only; production code, wire
equations, profiles and defaults are unchanged. No timing has been run and
this is not a performance promotion or a recovery-rate improvement.

`Generate.py` pins the original `WirehairSmallCore.h` and changes exactly two
coefficient loops in `Decoder::Feed`. An occupied pivot row has zeros before
its pivot and a diagonal coefficient of one. Eliminating its pivot therefore
sets the incoming coefficient to zero without multiplying by one. Normalizing
a new nonzero pivot sets its coefficient to one after computing the inverse.
Both remaining loops start at the next column.

This invariant holds for reverse pivot insertion and arbitrary accepted repair
rows. Successful recovery only clears coefficients above the diagonal;
subsequent full-rank feeds cannot insert a new pivot. A dependent contradiction
changes scratch storage, not the retained basis. Payload arithmetic, memory
allocation, ownership and layout are unchanged. The independent candidate is
not combined with the rejected coallocated-decoder or payload experiments.

## Correctness and mechanism evidence

The final tests on 2026-09-17 used the production sources at `309bf03` and the
files in this directory. Artifacts remain at
`/tmp/wh2-small-unit-diagonal.MfPjFJaS/`.

| Configuration | Final CTests |
| --- | ---: |
| GCC 13 native | 23/23 passed |
| GCC 13 portable arithmetic selector | 23/23 passed |
| GCC 13 full ASan/UBSan, leak detection and fake stacks | 23/23 passed |
| Clang 18 C++, GCC 13 C consumers | 23/23 passed |

Python 3.8 and 3.12 each passed both generator tests. The portable build uses
`ANDROID=1` to select arithmetic; it is not performance evidence from a real
non-GFNI host. Logs are the four `*-final-tests.log` files and each build's
`Testing/Temporary/LastTest.log`.

The namespace-isolated mechanism test uses an independent polynomial GF(256)
oracle, exhausts all 65,536 products and every nonzero pivot value, and tests
all six template dimensions (2, 3, 4, 5, 6, 8). Per arm and configuration:

| Mechanism counter | Result |
| --- | ---: |
| Cases / feeds / recoveries | 3,312 / 127,032 / 82,480 |
| Active eliminations / new pivots | 213,871 / 15,456 |
| Baseline / candidate scalar coefficient multiplications | 1,905,038 / 1,675,711 |
| Multiplications removed | 229,327 |
| Full-echelon checks / nonidentity bases | 3,312 / 1,698 |

The exact saving is one multiplication per active elimination plus one per
new pivot. These are operation counts, not measured speedups. Coverage includes
dense dependent equations before first recovery, reverse insertion, conflicts,
zero rows, distant IDs, partial tails, repeated recovery and guarded buffers.
K2/K4/K6 template coverage does not install those routes; the existing production
K6 implementation is separate.

Both arms also pass the unchanged public K3/K5/K8 and standalone small-codec
tests, C shared-library consumers, profile/borrowed/fault tests, allocation and
lifetime checks, and parity in both DSO load orders. Private/public decoder
allocation counts remain 2/3. All eight certified streams are 2,180,292 bytes
with SHA-256
`2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.

Native baseline and candidate DSO SHA-256 values, respectively:

```text
bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe
af37e1286af35fe7695e9d423aafd3d4cba2ea7127ef432ce89f807ef990d0ba
```

These identities are not a complete transitive build-provenance seal. That
seal and a separately frozen, bounded lifecycle/retention/WH1 timing protocol
remain tracked in the issue before any timing or production adoption.

## Reproduce in a fresh external directory

From the repository root, on 64-bit Linux:

```bash
wh2_diagonal_build=$(mktemp -d /tmp/wh2-unit-diagonal-check.XXXXXXXX)
python3 -B -m unittest discover -s bench/Wh2SmallUnitDiagonal -v
cmake -S bench/Wh2SmallUnitDiagonal -B "$wh2_diagonal_build" -G Ninja
cmake --build "$wh2_diagonal_build" -j 4
ctest --test-dir "$wh2_diagonal_build" --output-on-failure
```

Use separate fresh directories with `-DPORTABLE=ON`, `-DSANITIZE=ON`, or
`-DCMAKE_CXX_COMPILER=clang++-18` for the other configurations. Preserve the
original evidence directory; do not overwrite it with reproduction runs.
