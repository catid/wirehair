#!/usr/bin/env python3
"""Serial ABBA timing panels for the fixed-D12 two-anchor experiment."""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Sequence

from wh2_rank_floor_two_anchor_screen import (
    arm_options,
    canonical_json,
    parse_bench_output,
    sha256_bytes,
    sha256_file,
)


BASE_OPTIONS = (
    "precodefail",
    "--overhead", "0",
    "--threads", "1",
    "--loss", "0.50",
    "--schedule", "adversarial",
    "--completion", "mixed",
    "--mix-count", "2",
    "--packet-peel-seed-xor", "0",
    "--mixed-gf256-rows", "11",
    "--mixed-gf16-rows", "4",
    "--mixed-period", "32",
    "--mixed-geometry", "shared-x",
    "--mixed-residue-schedule", "hashed",
    "--mixed-residue-hash-seed", "68",
    "--mixed-residue-hash-keyed",
    "--mixed-independent-extension-residues",
    "--mixed-extension-residue-seed-xor", "78",
    "--mixed-independent-gf256-breaker-rows", "0",
    "--mixed-fused-schedule-buckets", "auto",
    "--seed-block-bytes", "1280",
    "--full-payload-solve",
)
METRICS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu", "block_xors_mu",
    "block_muladds_mu",
)
PANELS = (
    ("two_anchor_vs_d12", "d12", "two_anchor_adaptive"),
    ("d13_vs_d12", "d12", "d13_adaptive"),
    ("d14_vs_d12", "d12", "d14"),
)


def parse_u32_list(text: str) -> list[int]:
    result = [int(value, 10) for value in text.split(",")]
    if not result or any(value <= 0 or value > 0xFFFFFFFF for value in result):
        raise ValueError("expected a comma-separated positive u32 list")
    if len(result) != len(set(result)):
        raise ValueError("list values must be unique")
    return result


def parse_trial_map(text: str, widths: list[int]) -> dict[int, int]:
    result: dict[int, int] = {}
    for item in text.split(","):
        key, separator, value = item.partition(":")
        if not separator:
            raise ValueError("trial map entries must use width:trials")
        width = int(key, 10)
        trials = int(value, 10)
        if width <= 0 or trials <= 0 or width in result:
            raise ValueError("trial map requires unique positive values")
        result[width] = trials
    if set(result) != set(widths):
        raise ValueError("trial map must specify every requested width once")
    return result


def median(values: list[Decimal]) -> Decimal:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot take an empty median")
    middle = len(ordered) // 2
    if len(ordered) & 1:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / Decimal(2)


def run_phase(
    task_id: str,
    arm: str,
    binary: Path,
    K: int,
    width: int,
    trials: int,
    seed: str,
    cpu: int,
    result_dir: Path,
) -> dict[str, Any]:
    options = arm_options(arm, "large")
    command = [
        "taskset", "-c", str(cpu), str(binary), *BASE_OPTIONS,
        "--N", str(K), "--bb-list", str(width),
        "--trials", str(trials), "--seed", seed, *options,
    ]
    start_ns = time.time_ns()
    process = subprocess.run(
        command, check=False, text=True, encoding="utf-8",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    end_ns = time.time_ns()
    stdout_path = result_dir / "stdout" / f"{task_id}.csv"
    stderr_path = result_dir / "stderr" / f"{task_id}.txt"
    stdout_path.write_text(process.stdout, encoding="utf-8", newline="\n")
    stderr_path.write_text(process.stderr, encoding="utf-8", newline="\n")
    if process.returncode != 0 or process.stderr:
        raise RuntimeError(
            f"{task_id} failed rc={process.returncode}: {process.stderr[-1000:]}"
        )
    rows = parse_bench_output(
        process.stdout, (K,), "0x0", seed, "adversarial", options, True,
        trials, width
    )
    if len(rows) != 1:
        raise ValueError(f"{task_id} emitted {len(rows)} result rows")
    row = rows[0]
    if (int(row["N"]) != K or int(row["bb"]) != width or
            int(row["trials"]) != trials or int(row["success"]) != trials or
            int(row["rank_fail"]) != 0 or int(row["error"]) != 0):
        raise ValueError(f"{task_id} did not complete every solve successfully")
    return {
        "task_id": task_id, "arm": arm, "command": command,
        "start_ns": start_ns, "end_ns": end_ns,
        "returncode": process.returncode,
        "stdout": str(stdout_path.relative_to(result_dir)),
        "stdout_sha256": sha256_file(stdout_path),
        "stderr": str(stderr_path.relative_to(result_dir)),
        "stderr_sha256": sha256_file(stderr_path), "row": row,
    }


def summarize_panel(phases: Sequence[dict[str, Any]]) -> dict[str, Any]:
    if len(phases) != 4:
        raise ValueError("ABBA panel must contain four phases")
    reference = (phases[0], phases[3])
    candidate = (phases[1], phases[2])
    if (reference[0]["arm"] != reference[1]["arm"] or
            candidate[0]["arm"] != candidate[1]["arm"]):
        raise ValueError("panel is not ABBA ordered")
    metrics: dict[str, Any] = {}
    for metric in METRICS:
        reference_values = [Decimal(phase["row"][metric]) for phase in reference]
        candidate_values = [Decimal(phase["row"][metric]) for phase in candidate]
        reference_median = median(reference_values)
        candidate_median = median(candidate_values)
        metrics[metric] = {
            "reference_samples": [str(value) for value in reference_values],
            "candidate_samples": [str(value) for value in candidate_values],
            "reference_median": str(reference_median),
            "candidate_median": str(candidate_median),
            "candidate_over_reference": str(candidate_median / reference_median),
            "candidate_minus_reference": str(candidate_median - reference_median),
        }
    return {
        "reference": reference[0]["arm"],
        "candidate": candidate[0]["arm"], "metrics": metrics,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--K", type=int, default=48466)
    parser.add_argument("--widths", default="64,1280,4096")
    parser.add_argument("--trials", default="64:64,1280:24,4096:8")
    parser.add_argument("--seed", default="0x7a11ce0000000000")
    parser.add_argument("--cpu", type=int, default=-1)
    args = parser.parse_args()

    binary = args.binary.resolve(strict=True)
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise ValueError("--binary is not an executable regular file")
    if args.K < 4096 or args.K > 64000:
        raise ValueError("--K must be in [4096,64000] for active adaptive arms")
    widths = parse_u32_list(args.widths)
    trials_by_width = parse_trial_map(args.trials, widths)
    seed_value = int(args.seed, 0)
    if not 0 <= seed_value <= 0xFFFFFFFFFFFFFFFF:
        raise ValueError("--seed must fit u64")
    seed = f"0x{seed_value:016x}"
    cpus = sorted(os.sched_getaffinity(0))
    cpu = cpus[0] if args.cpu < 0 else args.cpu
    if cpu not in cpus:
        raise ValueError("--cpu is outside this process affinity")

    result_dir = args.result_dir.resolve()
    result_dir.mkdir(parents=True, exist_ok=False)
    (result_dir / "stdout").mkdir()
    (result_dir / "stderr").mkdir()
    getcontext().prec = 80

    all_phases: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for width in widths:
        width_summary: dict[str, Any] = {}
        trials = trials_by_width[width]
        for panel_name, reference, candidate in PANELS:
            sequence = (reference, candidate, candidate, reference)
            phases: list[dict[str, Any]] = []
            for phase_index, arm in enumerate(sequence):
                task_id = f"bb{width}.{panel_name}.phase{phase_index + 1}.{arm}"
                phase = run_phase(
                    task_id, arm, binary, args.K, width, trials,
                    seed, cpu, result_dir,
                )
                phases.append(phase)
                all_phases.append(phase)
            width_summary[panel_name] = summarize_panel(phases)
        summary[str(width)] = width_summary

    aggregate_reference = Decimal(0)
    aggregate_candidate = Decimal(0)
    equal_width_ratios: list[Decimal] = []
    for width in widths:
        solve = summary[str(width)]["two_anchor_vs_d12"]["metrics"][
            "solve_ms_mu"
        ]
        reference = Decimal(solve["reference_median"])
        candidate = Decimal(solve["candidate_median"])
        aggregate_reference += reference
        aggregate_candidate += candidate
        equal_width_ratios.append(candidate / reference)
    summary["aggregate_two_anchor_vs_d12"] = {
        "sum_reference_solve_ms": str(aggregate_reference),
        "sum_candidate_solve_ms": str(aggregate_candidate),
        "ratio_of_summed_solve_medians": str(
            aggregate_candidate / aggregate_reference
        ),
        "equal_width_mean_solve_ratio": str(
            sum(equal_width_ratios) / Decimal(len(equal_width_ratios))
        ),
    }

    with (result_dir / "commands.jsonl").open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        for phase in all_phases:
            output.write(canonical_json({
                key: value for key, value in phase.items() if key != "row"
            }) + "\n")
    with (result_dir / "results.jsonl").open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        for phase in all_phases:
            output.write(canonical_json({
                "task_id": phase["task_id"], "arm": phase["arm"],
                "row": phase["row"],
            }) + "\n")
    (result_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )

    repo = Path(__file__).resolve().parents[1]
    head = subprocess.check_output(("git", "rev-parse", "HEAD"), cwd=repo)
    diff = subprocess.check_output(
        ("git", "diff", "--binary", "--no-ext-diff"), cwd=repo
    )
    status = subprocess.check_output(
        ("git", "status", "--short"), cwd=repo
    ).decode("utf-8")
    helper = Path(__file__).resolve().with_name(
        "wh2_rank_floor_two_anchor_screen.py"
    )
    artifacts = {
        str(path.relative_to(result_dir)): sha256_file(path)
        for path in sorted(result_dir.rglob("*")) if path.is_file()
    }
    manifest = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_timing.v1",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "helper_script": str(helper), "helper_script_sha256": sha256_file(helper),
        "binary": str(binary), "binary_sha256": sha256_file(binary),
        "git_head": head.strip().decode("ascii"),
        "git_diff_sha256": sha256_bytes(diff), "git_status": status,
        "K": args.K, "widths": widths, "trials_by_width": trials_by_width,
        "seed": seed, "cpu": cpu, "panels": PANELS,
        "artifacts": artifacts,
    }
    manifest_path = result_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(canonical_json({
        "result_dir": str(result_dir),
        "manifest_sha256": sha256_file(manifest_path), "summary": summary,
    }))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"rank-floor two-anchor timing: {exc}", file=sys.stderr)
        raise
