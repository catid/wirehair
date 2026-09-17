# WH1 vs WH2 timing snapshot

The [CSV](wh1-vs-wh2.csv) transcribes the rounded timings recorded on
2026-09-17 in Beads issue `wirehair-za1o`. It compares the legacy Wirehair 1
public C API with the default, pure-GF(256) Wirehair 2 public C API in the
same process. It does **not** substitute an opt-in WH2 profile such as K6.

The recorded workload is an all-systematic, no-loss lifecycle probe:
encoder creation, encoding the K original packets, decoder creation,
decode and recovery. The notes identify WH1's `wirehair_encoder_create_ex`
and WH2's `wirehair_v2_encoder_create` / `wirehair_v2_decoder_create` paths.
K is the number of original blocks, B is bytes per block, and the full-block
message size is K × B. The plotted metric is elapsed lifecycle time in
milliseconds, not separate encoder or decoder throughput. Lower is better.

| Blocks (K) | Block bytes (B) | Message bytes | WH1 (ms) | WH2 (ms) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 64 | 512 | 0.003 | 0.007 |
| 128 | 64 | 8,192 | 0.029 | 0.038 |
| 512 | 64 | 32,768 | 0.132 | 0.170 |
| 1024 | 64 | 65,536 | 0.334 | 0.395 |
| 8 | 1280 | 10,240 | 0.012 | 0.013 |
| 128 | 1280 | 163,840 | 0.163 | 0.112 |
| 512 | 1280 | 655,360 | 0.639 | 0.595 |
| 1024 | 1280 | 1,310,720 | 1.288 | 1.182 |

## Scope and limitations

These are one-host diagnostic observations, not a speed qualification gate.
The retained issue summary has no raw samples, confidence intervals, exact
source/build identity, CPU identity, repetition count, or complete timing-boundary
specification. Do not infer these from the current checkout or machine.
In particular, the 0.001 ms rounding is coarse at K=8. Both plot panels use
the same logarithmic axes; lines connect measured points only as visual guides
and are not measurements of intervening sizes.

WH2 took less time at B=1280 for K=128, 512, and 1024, but more time at all
four B=64 points and at K=8/B=1280. This snapshot does not establish that WH2
is always faster, and it says nothing about performance under loss, recovery
failure probability, or required repair overhead. It is separate from the
older WH1-only measurements in the root README, which use a different workload
and machine/build context.

## Regenerate the plots

From the repository root, using only Python 3.8+ and its standard library:

```sh
python3 docs/benchmarks/plot_wh1_vs_wh2.py
python3 docs/benchmarks/plot_wh1_vs_wh2.py --check
python3 -m unittest discover -s docs/benchmarks -p 'test_*.py'
```

This reproduces the SVG from the checked-in CSV; it does **not** rerun the
timing experiment. A future measurement refresh should retain the benchmark
source, exact build/host metadata, timing boundaries, and raw repeated samples
alongside its results. Keep the default WH2 and WH1 workloads matched, and
report lossy recovery measurements separately from this no-loss lifecycle.
