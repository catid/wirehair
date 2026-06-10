# Row-Op Tiling Experiments

Standalone cache-tiling experiments for block-data row operations.  These do
not modify the production codec.

Build:

```bash
bash experiments/tiling/build.sh
```

Run:

```bash
./experiments/tiling/rowop_tiling
```

Quick correctness-only check:

```bash
./experiments/tiling/rowop_tiling --verify-only
```

The benchmark replays the same `dst ^= src` row-op schedule in two ways:

- untiled: each operation processes the whole symbol before the next operation;
- tiled: the full operation schedule is replayed independently for byte tiles.

The tiled order is algebraically valid for byte-wise RHS data because each tile
contains independent byte lanes and preserves operation order inside the tile.

## Real Codec Schedule Capture

Build the `bench/whx` harness with value-stage row-op logging enabled:

```bash
TAG=oplog OUT=/tmp/whx_oplog EXTRA='-DWH_OPLOG=1 -DWH_COUNT' bash bench/build.sh
```

Set `WH_OPLOG_PATH` while running a decode workload:

```bash
WH_OPLOG_PATH=/tmp/wirehair_oplog.csv /tmp/whx_oplog count --N 1000 --bb 102400 --loss 0.10
```

The CSV columns are:

```text
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
```

Logged stages are `InitializeColumnValues`, `MultiplyDenseValues`,
`AddSubdiagonalValues`, `BackSubstituteAboveDiagonal`, and `Substitute`.
`op_type` includes `zero`, `memcpy`, `xor`, `addset`, `add2`, `muladd`, and
`div`.  The offset fields are needed for final-block padding and partial final
block operations.
