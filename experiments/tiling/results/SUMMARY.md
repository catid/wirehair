# Row-Op Tiling Results

These are standalone row-op replay experiments, not production codec changes.
The schedule is synthetic but intentionally models repeated source-block fanout:
for each source block, eight `dst ^= src` operations are emitted before moving
to the next source.

Validation:

- `bash experiments/tiling/build.sh`
- `./experiments/tiling/rowop_tiling --verify-only`
- ASan/UBSan compile and `--verify-only`
- Full benchmark checksums matched for untiled and every tile size in each
  case.

Aggregate CSV:
`experiments/tiling/results/rowop_tiling_1280_100k_1m.csv`.

This persisted file is schema v1.  Its GiB/s field is a three-stream
operation estimate, not measured DRAM bandwidth.  All synthetic operations are
XORs, so the estimate corresponds to two explicit reads plus one write; timing,
checksums, and within-case speedups remain valid.

Best tile by block size:

| case | block bytes | untiled estimated-traffic GiB/s | best tile | tiled estimated-traffic GiB/s | speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| mtu1280 | 1280 | 57.614 | 16 KiB | 57.254 | 0.994x |
| kib100 | 102400 | 70.086 | 16 KiB | 90.032 | 1.285x |
| mib1 | 1048576 | 74.728 | 16 KiB | 156.202 | 2.090x |

**Stale-data caveat (2026-06-10):** the CSV and the mtu1280 row above predate
the e26a6cc harness fix.  Every published mtu1280 "tiled" row used tile sizes
(16 KiB+) at or above the 1280-byte block size, which the fixed harness now
skips as degenerate untiled replays; no genuine sub-block (256/512 B) tile was
ever measured for mtu1280.  Treat the mtu1280 readout below as "no evidence",
not "measured neutral", until the fixed binary regenerates this CSV.  The
kib100/mib1 rows used tiles below the block size and are unaffected.

Readout:

- Tiling is neutral/slightly negative for MTU-sized 1280-byte symbols.
  (See stale-data caveat: this row was a degenerate measurement.)
- For 100 KiB symbols, 16 KiB tiles improved synthetic row-op throughput by
  about 28%.
- For 1 MiB symbols, 16-64 KiB tiles roughly doubled synthetic row-op
  throughput, with 16 KiB best in this run.
- Larger tiles quickly lose the benefit: at 1 MiB, 256 KiB tiles were only
  1.20x over untiled, and 512 KiB was close to untiled.
- This supports `wirehair-2sc` as a real large-block experiment, but it does
  not prove end-to-end codec speedup.

## Actual codec WH_OPLOG replay (2026-07-03)

The replay loader was also run against real WH_OPLOG schedules captured from
the `bench/whx count` and `bench/whx repro` paths.  Captured schedules are
stored next to this summary:

- `wh_oplog_N128_bb1280_full_20260703.csv`
- `wh_oplog_N128_bb102400_full_20260703.csv`
- `wh_oplog_N128_bb1048576_full_20260703.csv`
- `wh_oplog_repro_N64_bb102400_partial_20260703.csv`

Aggregate replay CSV:
`experiments/tiling/results/rowop_tiling_real_oplog_20260703.csv`.

This is also schema v1.  Its hard-coded three-stream rate is only a legacy
normalization for mixed operation traces: zero, memcpy, addset, add2, muladd,
and divide have different read/write ledgers.  The recorded times, checksums,
and speedup ratios remain comparable because both replay arms use the same
schedule and denominator.  Schema v2 reports the mixed ledger explicitly.

Validation:

- `bash experiments/tiling/build.sh`
- `TAG=oplog OUT=/tmp/whx_oplog EXTRA='-DWH_OPLOG=1 -DWH_COUNT' bash bench/build.sh`
- `./experiments/tiling/rowop_tiling --verify-only`
- Full-block captures:
  `WH_OPLOG_PATH=/tmp/wirehair_oplog_N128_bb*.csv /tmp/whx_oplog count --N 128 --bb <1280|102400|1048576> --loss 0.10`
- Partial-final-block capture:
  `WH_OPLOG_PATH=/tmp/wirehair_oplog_repro_N64_bb102400_partial.csv /tmp/whx_oplog repro 0x12345678 64 102400 0 0.10`
- Every captured schedule passed `rowop_tiling --schedule ... --verify-only`.
  The partial trace had a final block of 85050/102400 bytes and included
  nonzero-offset padding rows.

Best real-schedule replay results:

| case | block bytes | ops | untiled legacy normalized GiB/s | best tile | tiled legacy normalized GiB/s | speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| WH_OPLOG N=128 full | 1280 | 4125 | 254.062 | 512 B | 180.402 | 0.710x |
| WH_OPLOG N=128 full | 102400 | 4172 | 193.461 | 16 KiB | 238.516 | 1.233x |
| WH_OPLOG N=128 full | 1048576 | 4157 | 95.699 | 16 KiB | 181.863 | 1.900x |
| WH_OPLOG N=64 partial final block | 102400 | 1910 | 185.963 | 8 KiB | 206.941 | 1.113x |

Readout:

- Actual codec schedules confirm the synthetic shape: sub-block tiling is a
  loss for MTU-sized symbols, modestly helpful around 100 KiB, and strongly
  helpful at 1 MiB.
- 16 KiB is the best full-block tile in these runs for both 100 KiB and 1 MiB
  symbols; 4-32 KiB is the useful region.
- The partial-final-block replay matched untiled output across all tested tile
  sizes, covering the padding-sensitive schedule rows.
- These experiments validate the schedule transformation and throughput
  direction.  Any production optimization should still be gated by block size
  and measured end-to-end in the composed decoder path.
