# Precode and table-regeneration handoff (2026-07-10)

## Git coordinates

- Remote branch: `origin/codex/precode-handoff-20260710`
- Branch point: `6a6d3d9dc8a33a169298f7565b067e7ede048cd3`
- Implementation checkpoint: `90fed75e072976f273f38781e90f7fd18aca3227`
- The branch point was also `origin/master` when this handoff was prepared. The
  checkpoint was intentionally kept separate so that newer work can merge it
  and resolve conflicts in context.

The branch contains one large implementation checkpoint followed by this
handoff note. It is meant to be merged as a branch; the implementation commit
can also be reviewed independently as `6a6d3d9..90fed75`.

## Work represented by the branch

### `wirehair-axd`: version-3 V2 precode codec

The isolated V2 path now has a real packet-equation contract, intermediate
value solver, precode decoder, message-level encoder/decoder adapters, facade
entry points, deterministic seed-attempt selection, and focused solve,
round-trip, and seed-selection tests. The main implementation is in
`codec/WirehairV2Solve.*`, `codec/WirehairV2PrecodeDecode.*`, and the expanded
`codec/WirehairV2PrecodeEncode.*` files. `codec/README.md` describes the packet
contract and the important integrity limitation: this is an erasure codec, not
packet authentication.

The V1-compatible facade path remains distinct. V2 packets use the explicit
precode initializers and bind all matrix-affecting profile fields. Do not weaken
the decoder's profile/options validation while resolving merge conflicts.

Benchmark support adds real precode compare/failure modes, and the CI changes
exercise explicit tools in hosted GCC/Clang lanes. The Beads issue remains
`in_progress`; its full acceptance criteria and accumulated design notes are
the source of truth:

```sh
bd show wirehair-axd
```

### `wirehair-jwg`: offline table regeneration hardening

The table tools gain deterministic selection policies, atomic/no-clobber result
publication, bounded input handling, strict holdout/fixup ledgers, regeneration
result validation/application, a full-regeneration driver, and Python/C++
tests. The generator and policy changes are deliberately strict about
provenance and partial output.

This checkpoint does not assert that the expensive full regeneration campaign
finished. Do not promote partial campaign output or alter manifest-frozen
inputs based only on this branch. Current campaign status, deferred launch
gates, and exact result hashes live in:

```sh
bd show wirehair-jwg
```

Beads data was pushed with `bd dolt push`. `.beads/issues.jsonl` is only a
passive export; if it conflicts during integration, do not use `bd import` as a
merge mechanism.

## Validation completed for the checkpoint

The handoff pass completed these checks successfully:

```sh
git diff --cached --check
python3 -m unittest discover -s tables -p 'test_*.py'  # 47 passed
python3 -m unittest ci.test_run_ci                     # 19 passed
bash -n bench/test_whx.sh tables/aggregate_dense_count_shards.sh \
  tables/derive_dense_count_candidate.sh \
  tables/run_regeneration_slice.sh tables/run_full_regeneration.sh
python3 -m py_compile ci/run_ci.py ci/test_run_ci.py tables/*.py
```

A fresh C++ rebuild and complete CTest run were not performed during the short
branch-publication pass. After merging onto the receiving branch, use a fresh
build directory so older generated state cannot hide integration problems:

```sh
cmake -S . -B build/handoff -G Ninja \
  -DBUILD_TESTS=ON -DBUILD_CODEC_V2=ON \
  -DWIREHAIR_BUILD_TOOLS=ON -DWIREHAIR_BUILD_BENCHMARKS=ON \
  -DWIREHAIR_BUILD_VALIDATION_HARNESS=ON
cmake --build build/handoff -j
ctest --test-dir build/handoff --output-on-failure
```

The most relevant focused CTest names contain `wirehair_v2_precode`,
`wirehair_v2_policy`, `atomic_result_file`, or `peel_selection_policy`.

## Deliberately excluded local files

Four protected scratch/result files were not committed or deleted:

```text
experiments/precode/results/codecport_cert_K64000_20260703.csv
experiments/precode/results/codecport_n13_ic_cert_K10000_32000_64000_20260707.csv.err.tmp.GzXKp4
experiments/precode/results/codecport_n13_ic_cert_K10000_32000_64000_20260707.csv.tmp.j2SLv0
probe1.cpp
```

The first CSV is partial, and the two temporary files are interrupted campaign
artifacts. They are not certification evidence and should not be swept into a
merge commit. The files remain only in the originating worktree and are not on
the remote handoff branch.
