# Peeling Experiments

This directory is reserved for standalone peeling-schedule experiments.
Experimental heuristics should live here instead of behind compile-time branches
in `WirehairCodec.cpp`.

Keep production codec changes separate from sweep harnesses:
- reuse public headers or copied experiment-only model code here;
- report residual unpeeled rows/cols and peel-side runtime separately from dense
  matrix behavior;
- remember that `N` is the number of blocks, while bytes are `N * block_bytes`.
