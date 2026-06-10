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

Best tile by block size:

| case | block bytes | untiled GiB/s | best tile | tiled GiB/s | speedup |
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
  not prove end-to-end codec speedup.  A production prototype still needs to
  log actual `AddSubdiagonalValues`/`BackSubstituteAboveDiagonal`/`Substitute`
  row-op schedules and replay those schedules by tiles while preserving final
  block padding behavior.
