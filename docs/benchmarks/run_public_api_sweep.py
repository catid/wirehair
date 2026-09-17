#!/usr/bin/env python3
"""Run the bounded public-API WH1/WH2 size sweep used by the README plots."""

import argparse
import csv
import subprocess
from pathlib import Path


K_POINTS = (8, 128, 512, 1024) + tuple(range(1000, 64001, 1000))
BLOCKS = (64, 1280)
PROFILE_CERTIFIED = 0x4B295BBB47F4F9C9


def run(binary, mode, k, block_bytes, trials, cpu):
    command = [str(binary), mode, str(k), str(block_bytes), str(trials)]
    if cpu is not None:
        command = ["taskset", "-c", str(cpu)] + command
    completed = subprocess.run(command, check=True, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != trials:
        raise RuntimeError("{} K={} B={} returned {} rows, expected {}".format(
            mode, k, block_bytes, len(rows), trials))
    for row in rows:
        if int(row["blocks"]) != k or int(row["block_bytes"]) != block_bytes:
            raise RuntimeError("dimension mismatch in {} K={} B={}".format(
                mode, k, block_bytes))
        if int(row["wh2_profile"]) != PROFILE_CERTIFIED:
            raise RuntimeError("default WH2 profile identity changed")
    return rows


def write_rows(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True,
                        help="compiled PublicApiSweep executable")
    parser.add_argument("--speed-trials", type=int, default=5)
    parser.add_argument("--recovery-trials", type=int, default=16)
    parser.add_argument("--cpu", type=int,
                        help="optional CPU for taskset pinning")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    if not args.binary.is_file() or args.speed_trials < 1 or args.recovery_trials < 1:
        parser.error("binary must exist and trial counts must be positive")

    speed_raw, recovery_raw = [], []
    for block_bytes in BLOCKS:
        for k in K_POINTS:
            speed_raw.extend(run(args.binary, "timing", k, block_bytes,
                                 args.speed_trials, args.cpu))
            recovery_raw.extend(run(args.binary, "recovery", k, block_bytes,
                                    args.recovery_trials, args.cpu))
            print("completed K={} B={}".format(k, block_bytes), flush=True)

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    speed_fields = ["blocks", "block_bytes", "speed_trials", "loss_rate",
                    "wh1_lifecycle_ms", "wh2_lifecycle_ms"]
    speed_rows = []
    for block_bytes in BLOCKS:
        for k in K_POINTS:
            rows = [r for r in speed_raw if int(r["blocks"]) == k and
                    int(r["block_bytes"]) == block_bytes]
            speed_rows.append({
                "blocks": k, "block_bytes": block_bytes,
                "speed_trials": len(rows), "loss_rate": 0.0,
                "wh1_lifecycle_ms": sum(int(r["wh1_ns"]) for r in rows) /
                len(rows) / 1e6,
                "wh2_lifecycle_ms": sum(int(r["wh2_ns"]) for r in rows) /
                len(rows) / 1e6,
            })
    recovery_fields = ["blocks", "block_bytes", "trials", "loss_rate",
                       "wh1_exact_k_rate", "wh2_exact_k_rate",
                       "wh1_eventual_rate", "wh2_eventual_rate",
                       "wh1_terminal_failures", "wh2_terminal_failures"]
    recovery_rows = []
    for block_bytes in BLOCKS:
        for k in K_POINTS:
            rows = [r for r in recovery_raw if int(r["blocks"]) == k and
                    int(r["block_bytes"]) == block_bytes]
            n = len(rows)
            wh1 = [int(r["wh1_first"]) for r in rows]
            wh2 = [int(r["wh2_first"]) for r in rows]
            for first in wh1 + wh2:
                if first and not k <= first <= k + 4:
                    raise RuntimeError("recovery horizon violation K={}".format(k))
            recovery_rows.append({
                "blocks": k, "block_bytes": block_bytes, "trials": n,
                "loss_rate": 0.10,
                "wh1_exact_k_rate": sum(first == k for first in wh1) / n,
                "wh2_exact_k_rate": sum(first == k for first in wh2) / n,
                "wh1_eventual_rate": sum(first >= k for first in wh1) / n,
                "wh2_eventual_rate": sum(first >= k for first in wh2) / n,
                "wh1_terminal_failures": sum(first == 0 for first in wh1),
                "wh2_terminal_failures": sum(first == 0 for first in wh2),
            })

    write_rows(output / "wh1-vs-wh2-large-k-speed-trials.csv", speed_raw,
               list(speed_raw[0]))
    write_rows(output / "wh1-vs-wh2-large-k-recovery-trials.csv", recovery_raw,
               list(recovery_raw[0]))
    write_rows(output / "wh1-vs-wh2-large-k-speed.csv", speed_rows, speed_fields)
    write_rows(output / "wh1-vs-wh2-large-k-recovery.csv", recovery_rows,
               recovery_fields)


if __name__ == "__main__":
    main()
