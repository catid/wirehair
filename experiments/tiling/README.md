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

Replay the captured schedule in this harness:

```bash
./experiments/tiling/rowop_tiling --schedule /tmp/wirehair_oplog.csv --block-bytes 102400 --tiles 4096,8192,16384,32768,65536 --repeats 5
```

Trace replay validates its peak simultaneous storage before allocating replay
buffers or printing a CSV header.  The default `--max-memory-mib 256` aggregate
policy includes the operation-vector capacity, timing samples, input, and base,
untiled, and tiled recovery copies.  It also reserves 1 MiB for the bounded CSV
parser, stream buffer, and allocator metadata.  Adjust the policy explicitly
for a known larger workload.  Trace lines above 4096 bytes and traces above one
million operations are rejected so malformed input cannot grow parser storage,
the replay schedule, or the byte ledger without bound; `--max-operations` can
lower the operation ceiling for constrained validation runs.

The CSV columns are:

```text
stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes
```

Logged stages are `InitializeColumnValues`, `MultiplyDenseValues`,
`AddSubdiagonalValues`, `BackSubstituteAboveDiagonal`, and `Substitute`.
`op_type` includes `zero`, `memcpy`, `xor`, `addset`, `add2`, `muladd`, and
`div`.  The offset fields are needed for final-block padding and partial final
block operations.

## Byte-accounting schema

Schema-v2 result rows preserve `logical_gib`, `memory_gib_3stream`, and
`gib_per_s` for historical readers, then append explicit logical-work,
estimated-read, estimated-write, and estimated-total-traffic columns.  The
per-operation ledger distinguishes zero, memcpy, XOR, addset, add2, muladd,
and divide instead of applying three streams to every trace operation.

Logical work is the destination span processed.  Estimated traffic counts
explicit operands and destination writes but excludes write-allocate,
cache-line rounding, prefetching, and cache effects; it is not measured DRAM
bandwidth.  Validate old and new outputs with:

```bash
python3 experiments/validate_byte_metrics.py result.csv
```
