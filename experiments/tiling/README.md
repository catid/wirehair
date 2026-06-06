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
