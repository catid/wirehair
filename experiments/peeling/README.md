# Peeling Experiments

This directory is reserved for standalone peeling-schedule experiments.
Experimental heuristics should live here instead of behind compile-time branches
in `WirehairCodec.cpp`.

Keep production codec changes separate from sweep harnesses:
- reuse public headers or copied experiment-only model code here;
- report residual unpeeled rows/cols and peel-side runtime separately from dense
  matrix behavior;
- remember that `N` is the number of blocks, while bytes are `N * block_bytes`.

Build:

```bash
bash experiments/peeling/build.sh
```

Example sweep:

```bash
./experiments/peeling/peel_sweep --N 128,2048,32000 --trials 40 --source loss --loss 0.10
```

The harness models only the peel graph. It feeds generated Wirehair peel rows
sequentially, allows singleton rows to avalanche immediately, then applies a
selected greedy inactivation schedule. Dense matrix success/failure is outside
the measurement target.
