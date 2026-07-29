#!/usr/bin/env python3
"""Run the pre-result-frozen za5v survivor S/D2-neighbor screen."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import signal
import subprocess
import sys
import time


REPO = Path(
    "/home/catid/wirehair/.claude/worktrees/"
    "wf_1b4ec124-567-2"
).resolve()
sys.path.insert(0, str(REPO / "tools"))

import band_timing_codec as band
import peel_codec


SCHEMA = "wirehair.wh2.za5v.local-neighbors.v1"
K_VALUES = (2, 4, 5, 9, 12, 40, 56, 92, 94, 100)
CONSTRUCTION_SEED = 0x6789CAFE
LOSS_SEED = 0x243F6A8885A308D3
ORDER_SEED = 0x5A5A4C4F43414C31
THERMAL_SOURCE = Path(
    "/tmp/wh2-za5v-bandtiming-thermal-a1e1997.csv"
)
BENCH = REPO / "build-audit/codec/wirehair_v2_bench"


def candidate_roster():
    roster = [
        ("dispatch10p2_s0_d4_frozen", 10, 2, "frozen", 0, 4),
    ]
    for gf256 in (8, 9):
        for staircase_offset in (-1, 0, 1):
            for dense_rows in (3, 4, 5):
                roster.append((
                    f"pure{gf256}_s{staircase_offset:+d}_d{dense_rows}_"
                    "frozen",
                    gf256, 0, "frozen", staircase_offset, dense_rows,
                ))
    return tuple(roster)


CANDIDATES = candidate_roster()


def canonical_json(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    )


def sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def runtime_bindings():
    thermal_stat = THERMAL_SOURCE.stat()
    return {
        "bench_sha256": sha256_bytes(BENCH.read_bytes()),
        "runner_sha256": sha256_bytes(
            Path(__file__).resolve().read_bytes()
        ),
        "parser_sha256": sha256_bytes(
            Path(band.__file__).resolve().read_bytes()
        ),
        "context_tool_sha256": sha256_bytes(
            Path(peel_codec.__file__).resolve().read_bytes()
        ),
        "thermal_source_device": thermal_stat.st_dev,
        "thermal_source_inode": thermal_stat.st_ino,
    }


def atomic_json(path, value):
    payload = (canonical_json(value) + "\n").encode("ascii")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with open(temporary, "xb") as output:
        output.write(payload)
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary, path)
    return sha256_bytes(payload)


def descriptor_for(job):
    dispatch = band.dispatch_band_descriptor(job["K"])
    return band.BandDescriptor(
        staircase=max(1, dispatch.staircase + job["staircase_offset"]),
        dense_rows=job["dense_rows"],
        gf256_rows=job["gf256"],
        gf16_rows=job["gf16"],
        period=244,
        x_mode=job["x"],
    )


def expected_request(job):
    return {
        "block_count": job["K"],
        "block_bytes": 2,
        "candidate": descriptor_for(job),
        "construction_seed": CONSTRUCTION_SEED,
        "loss": 0.1,
        "loss_seed": LOSS_SEED,
        "schedule": "iid",
        "warmup_replicates": 0,
        "replicates": 3,
        "inner_reps": 1024,
        "max_overhead": 16,
        "cache_state": "warm",
        "systematic_cache": job["cache"],
        "evict_bytes": 4096,
        "required_margin": 0.0,
    }


def run_worker(job_path, output_path):
    with open(job_path, encoding="ascii") as source:
        job = json.load(source)
    affinity = sorted(os.sched_getaffinity(0))
    if affinity != [job["cpu"]]:
        raise RuntimeError(
            f"worker affinity {affinity!r} != requested {[job['cpu']]!r}"
        )
    before_bindings = runtime_bindings()
    if before_bindings != job["bindings"]:
        raise RuntimeError("worker runtime binding changed before launch")
    request = expected_request(job)
    context = peel_codec.make_paired_context_config(
        str(THERMAL_SOURCE),
        1000,
        cache_state=request["cache_state"],
        evict_bytes=request["evict_bytes"],
    )
    started = time.time_ns()
    measurement = band.bandtiming_probe(
        str(BENCH),
        request["block_count"],
        request["block_bytes"],
        request["candidate"],
        construction_seed=request["construction_seed"],
        loss=request["loss"],
        loss_seed=request["loss_seed"],
        schedule=request["schedule"],
        warmup_replicates=request["warmup_replicates"],
        replicates=request["replicates"],
        inner_reps=request["inner_reps"],
        max_overhead=request["max_overhead"],
        cache_state=request["cache_state"],
        systematic_cache=request["systematic_cache"],
        evict_bytes=request["evict_bytes"],
        context=context,
        required_margin=request["required_margin"],
    )
    finished = time.time_ns()
    receipt = measurement.as_dict()
    band.replay_bandtiming_receipt(
        receipt,
        expected_request=request,
    )
    if (receipt["context"]["bound"]["thermal_device"] !=
            job["bindings"]["thermal_source_device"] or
            receipt["context"]["bound"]["thermal_inode"] !=
            job["bindings"]["thermal_source_inode"]):
        raise RuntimeError("worker receipt used a different thermal source")
    after_bindings = runtime_bindings()
    if after_bindings != before_bindings:
        raise RuntimeError("worker runtime binding changed during measurement")
    envelope = {
        "schema": SCHEMA,
        "job": job,
        "wall_started_unix_ns": started,
        "wall_finished_unix_ns": finished,
        "bindings_before": before_bindings,
        "bindings_after": after_bindings,
        "receipt": receipt,
    }
    atomic_json(output_path, envelope)


def jobs():
    result = []
    sequence = 0
    for (name, gf256, gf16, x_mode, staircase_offset,
         dense_rows) in CANDIDATES:
        for block_count in K_VALUES:
            for cache in ("off", "on"):
                result.append({
                    "sequence": sequence,
                    "candidate": name,
                    "K": block_count,
                    "gf256": gf256,
                    "gf16": gf16,
                    "x": x_mode,
                    "staircase_offset": staircase_offset,
                    "dense_rows": dense_rows,
                    "cache": cache,
                })
                sequence += 1
    random.Random(ORDER_SEED).shuffle(result)
    for order, job in enumerate(result):
        job["order"] = order
    return result


def contrast(receipt, name):
    return next(
        item for item in receipt["contrasts"] if item["name"] == name
    )


def contrast_summary(receipt, name):
    item = contrast(receipt, name)
    eligible = item["eligible_replicates"]
    return {
        f"{name}_eligible_replicates": eligible,
        f"{name}_cost_ratio":
            math.exp(item["log_cost_mean"]) if eligible else "",
        f"{name}_recovery_regressions": item["recovery_regressions"],
        f"{name}_recovery_improvements": item["recovery_improvements"],
        f"{name}_both_failures": item["both_failures"],
    }


def summarize(directory, manifest, completed):
    fieldnames = (
        "order", "candidate", "K", "cache", "cpu", "stream_sha256",
        "weak_cells", "censored_panel_observations",
        "fault_contaminated_panel_observations",
        "encoder_candidate_dispatch_eligible_replicates",
        "encoder_candidate_dispatch_cost_ratio",
        "encoder_candidate_dispatch_recovery_regressions",
        "encoder_candidate_dispatch_recovery_improvements",
        "encoder_candidate_dispatch_both_failures",
        "encoder_candidate_wh1_eligible_replicates",
        "encoder_candidate_wh1_cost_ratio",
        "encoder_candidate_wh1_recovery_regressions",
        "encoder_candidate_wh1_recovery_improvements",
        "encoder_candidate_wh1_both_failures",
        "decoder_candidate_dispatch_eligible_replicates",
        "decoder_candidate_dispatch_cost_ratio",
        "decoder_candidate_dispatch_recovery_regressions",
        "decoder_candidate_dispatch_recovery_improvements",
        "decoder_candidate_dispatch_both_failures",
        "decoder_candidate_wh1_eligible_replicates",
        "decoder_candidate_wh1_cost_ratio",
        "decoder_candidate_wh1_recovery_regressions",
        "decoder_candidate_wh1_recovery_improvements",
        "decoder_candidate_wh1_both_failures",
        "direct_candidate_dispatch_eligible_replicates",
        "direct_candidate_dispatch_cost_ratio",
        "direct_candidate_dispatch_recovery_regressions",
        "direct_candidate_dispatch_recovery_improvements",
        "direct_candidate_dispatch_both_failures",
        "candidate_decoder_successes", "dispatch_decoder_successes",
        "wh1_decoder_successes", "cpu_max_c", "dimm_max_c",
    )
    rows = []
    file_hashes = []
    for job in manifest["jobs"]:
        path = directory / job["output"]
        payload = path.read_bytes()
        file_hashes.append({
            "output": job["output"],
            "sha256": sha256_bytes(payload),
        })
        envelope = json.loads(payload)
        if envelope["schema"] != SCHEMA or envelope["job"] != job:
            raise RuntimeError(f"job envelope mismatch: {path}")
        if (envelope["bindings_before"] != manifest["bindings"] or
                envelope["bindings_after"] != manifest["bindings"]):
            raise RuntimeError(f"job runtime binding mismatch: {path}")
        receipt = envelope["receipt"]
        band.replay_bandtiming_receipt(
            receipt,
            expected_request=expected_request(job),
        )
        if receipt["candidate_descriptor"] != descriptor_for(job).as_dict():
            raise RuntimeError(f"candidate descriptor mismatch: {path}")
        if (receipt["context"]["bound"]["thermal_device"] !=
                manifest["bindings"]["thermal_source_device"] or
                receipt["context"]["bound"]["thermal_inode"] !=
                manifest["bindings"]["thermal_source_inode"]):
            raise RuntimeError(f"job thermal source mismatch: {path}")
        arms = {
            (item["scope"], item["arm"]): item
            for item in receipt["arms"]
        }
        thermal = receipt["context"]["thermal"]
        rows.append({
            "order": job["order"],
            "candidate": job["candidate"],
            "K": job["K"],
            "cache": job["cache"],
            "cpu": job["cpu"],
            "stream_sha256": receipt["stream_sha256"],
            "weak_cells": len(receipt["weak_cells"]),
            "censored_panel_observations": sum(
                len(item["censored_panels"])
                for item in receipt["replicates"]
            ),
            "fault_contaminated_panel_observations": sum(
                len(item["fault_contaminated_panels"])
                for item in receipt["replicates"]
            ),
            **contrast_summary(receipt, "encoder_candidate_dispatch"),
            **contrast_summary(receipt, "encoder_candidate_wh1"),
            **contrast_summary(receipt, "decoder_candidate_dispatch"),
            **contrast_summary(receipt, "decoder_candidate_wh1"),
            **contrast_summary(receipt, "direct_candidate_dispatch"),
            "candidate_decoder_successes":
                arms[("decoder", "candidate")]["successes"],
            "dispatch_decoder_successes":
                arms[("decoder", "dispatch")]["successes"],
            "wh1_decoder_successes":
                arms[("decoder", "wh1")]["successes"],
            "cpu_max_c": thermal["cpu_tctl_max_c"],
            "dimm_max_c": thermal["dimm_max_c"],
        })
    summary_path = directory / "summary.csv"
    temporary = summary_path.with_suffix(".csv.tmp")
    with open(temporary, "x", newline="", encoding="ascii") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        output.flush()
        os.fsync(output.fileno())
    os.replace(temporary, summary_path)
    if runtime_bindings() != manifest["bindings"]:
        raise RuntimeError("campaign runtime binding changed before completion")
    completed.update({
        "job_files": file_hashes,
        "summary_sha256": sha256_bytes(summary_path.read_bytes()),
        "completed_unix_ns": time.time_ns(),
    })
    atomic_json(directory / "complete.json", completed)


def run_campaign(directory):
    if directory.exists():
        raise RuntimeError(f"refusing existing campaign directory {directory}")
    if not BENCH.is_file() or not os.access(BENCH, os.X_OK):
        raise RuntimeError(f"benchmark is not executable: {BENCH}")
    if not THERMAL_SOURCE.is_file():
        raise RuntimeError(f"thermal source does not exist: {THERMAL_SOURCE}")
    cpus = sorted(os.sched_getaffinity(0))
    if not cpus:
        raise RuntimeError("campaign has no allowed CPUs")
    directory.mkdir(parents=True)
    job_dir = directory / "jobs"
    result_dir = directory / "results"
    job_dir.mkdir()
    result_dir.mkdir()
    campaign_jobs = jobs()
    bindings = runtime_bindings()
    for index, job in enumerate(campaign_jobs):
        job["cpu"] = cpus[index % len(cpus)]
        job["output"] = f"results/{index:03d}.json"
        job["bindings"] = bindings
    manifest = {
        "schema": SCHEMA,
        "created_unix_ns": time.time_ns(),
        "repository": str(REPO),
        "bench": str(BENCH),
        "bindings": bindings,
        "thermal_source": str(THERMAL_SOURCE),
        "K_values": list(K_VALUES),
        "candidates": [
            {
                "name": name, "gf256": gf256, "gf16": gf16,
                "x": x_mode, "staircase_offset": staircase_offset,
                "dense_rows": dense_rows,
            }
            for (name, gf256, gf16, x_mode, staircase_offset,
                 dense_rows) in CANDIDATES
        ],
        "construction_seed": CONSTRUCTION_SEED,
        "loss_seed": LOSS_SEED,
        "order_seed": ORDER_SEED,
        "cpus": cpus,
        "jobs": campaign_jobs,
    }
    manifest_sha = atomic_json(directory / "manifest.json", manifest)
    for index, job in enumerate(campaign_jobs):
        atomic_json(job_dir / f"{index:03d}.json", job)

    script = Path(__file__).resolve()
    for wave_start in range(0, len(campaign_jobs), len(cpus)):
        if runtime_bindings() != bindings:
            raise RuntimeError("campaign runtime binding changed between waves")
        wave = campaign_jobs[wave_start:wave_start + len(cpus)]
        processes = []
        try:
            for offset, job in enumerate(wave):
                index = wave_start + offset
                command = [
                    "taskset", "-c", str(job["cpu"]),
                    sys.executable, str(script), "--worker",
                    "--job", str(job_dir / f"{index:03d}.json"),
                    "--output", str(directory / job["output"]),
                ]
                log_path = directory / f"worker-{index:03d}.log"
                log = open(log_path, "xb")
                try:
                    process = subprocess.Popen(
                        command, cwd=REPO, stdout=log,
                        stderr=subprocess.STDOUT, start_new_session=True,
                    )
                except BaseException:
                    log.close()
                    raise
                processes.append((index, process, log))
            failures = []
            for index, process, unused_log in processes:
                status = process.wait()
                if status != 0:
                    failures.append((index, status))
            if failures:
                raise RuntimeError(f"campaign wave failed: {failures!r}")
        except BaseException:
            for unused_index, process, unused_log in processes:
                if process.poll() is None:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
            for unused_index, process, unused_log in processes:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
            raise
        finally:
            for unused_index, unused_process, log in processes:
                if not log.closed:
                    log.close()
        if runtime_bindings() != bindings:
            raise RuntimeError("campaign runtime binding changed during wave")
    summarize(
        directory,
        manifest,
        {
            "schema": SCHEMA,
            "manifest_sha256": manifest_sha,
            "jobs": len(campaign_jobs),
        },
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--job", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--directory", type=Path)
    args = parser.parse_args()
    if args.worker:
        if args.job is None or args.output is None or args.directory:
            parser.error("--worker requires --job and --output only")
        run_worker(args.job, args.output)
    else:
        if args.directory is None or args.job or args.output:
            parser.error("campaign mode requires --directory only")
        run_campaign(args.directory.resolve())


if __name__ == "__main__":
    main()
