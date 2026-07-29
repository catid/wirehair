#!/usr/bin/env python3
"""Frozen all-K reliability/timing campaign for bead wirehair-za5v.

This runner deliberately separates the immutable statistical population from
runtime CPU assignment.  The selection population and every possible holdout
population can therefore be hashed before a CPU is chosen.  Each native
receipt is replayed with :mod:`band_timing_codec` both immediately after the
probe and again while the final summaries are built.

The large campaign is intentionally not started by importing this module.
Use ``plan`` to inspect the frozen logical population and ``run`` to execute
one phase into a directory that must not already exist.
"""

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
import time
import zlib


REPOSITORY = Path(__file__).resolve().parents[1]
TOOLS = REPOSITORY / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import band_timing_codec as band
import peel_codec


CAMPAIGN_SCHEMA = "wirehair.wh2.za5v.all-k-campaign.v1"
JOB_SCHEMA = "wirehair.wh2.za5v.all-k-job.v1"
RESULT_SCHEMA = "wirehair.wh2.za5v.all-k-result.v1"
SUMMARY_SCHEMA = "wirehair.wh2.za5v.all-k-summary.v1"
COMPLETION_SCHEMA = "wirehair.wh2.za5v.all-k-completion.v1"
CELL_LEDGER_SCHEMA = "wirehair.wh2.za5v.all-k-cell.v1"
DECISION_SCHEMA = "wirehair.wh2.za5v.selection-decision.v1"

K_VALUES = tuple(range(2, 101))
SCHEDULES = (
    "iid",
    "burst",
    "permutation",
    "systematic-first",
    "repair-only",
    "adversarial",
)

BLOCK_BYTES = 2
LOSS = 0.10
WARMUP_REPLICATES = 0
INNER_REPS = 1
MAX_OVERHEAD = 64
CACHE_STATE = "warm"
SYSTEMATIC_CACHE = "off"
EVICT_BYTES = 4096
REQUIRED_MARGIN = 0.02
THERMAL_INTERVAL_MS = 1000
MAX_WORKERS = 32
RECOVERY_PAIRS = (
    ("candidate", "dispatch"),
    ("candidate", "wh1"),
    ("dispatch", "wh1"),
)

_PARENT_SIGNAL_STATE = None


def _parent_signal_error(signum):
    return CampaignError(
        "campaign parent received "
        f"{signal.Signals(signum).name}")


def _raise_pending_parent_signal(state):
    for item in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(item, signal.SIG_IGN)
    raise _parent_signal_error(state["pending"])


def _parent_signal_safe_point():
    """Raise a recorded TERM/HUP only where cleanup is guaranteed."""
    state = _PARENT_SIGNAL_STATE
    if state is not None and state["pending"] is not None:
        _raise_pending_parent_signal(state)

SELECTION_REPLICATES = 768
HOLDOUT_REPLICATES = 256
SELECTION_CONSTRUCTION_BASE = 0x13579BDF
SELECTION_LOSS_BASE = 0x0123456789ABCDEF
HOLDOUT_CONSTRUCTION_BASE = 0x2468ACE0
HOLDOUT_LOSS_BASE = 0xFEDCBA9876543210
SELECTION_ORDER_SEED = 0x5E1EC71020260728
HOLDOUT_ORDER_SEED = 0xA01D0E7020260728

RESULT_CLASSES = ("success", "need_more", "weak", "error")
TIMING_ARM_KEYS = (
    ("encoder", "candidate"),
    ("encoder", "dispatch"),
    ("encoder", "wh1"),
    ("decoder", "candidate"),
    ("decoder", "dispatch"),
    ("decoder", "wh1"),
    ("direct", "candidate"),
    ("direct", "dispatch"),
)
DECISION_TIMING_PANELS = (
    (
        "direct_candidate_dispatch_log_cost",
        "direct_candidate_dispatch",
        "direct",
        "candidate",
        "dispatch",
    ),
    (
        "encoder_candidate_wh1_log_cost",
        "encoder_candidate_wh1",
        "encoder",
        "candidate",
        "wh1",
    ),
    (
        "decoder_candidate_wh1_log_cost",
        "decoder_candidate_wh1",
        "decoder",
        "candidate",
        "wh1",
    ),
)
CANDIDATE_CONTROL_PANELS = (
    ("encoder_candidate_dispatch", "encoder", "dispatch"),
    ("encoder_candidate_wh1", "encoder", "wh1"),
    ("decoder_candidate_dispatch", "decoder", "dispatch"),
    ("decoder_candidate_wh1", "decoder", "wh1"),
    ("direct_candidate_dispatch", "direct", "dispatch"),
)
PRODUCTION_BLOCKER_SCHEMA = "wirehair.wh2.za5v.production-blockers.v1"

# The names are part of the frozen experiment contract.  ``staircase_delta``
# is applied to the exact dispatch-v1 staircase S0 at each K.
CANDIDATES = (
    {
        "name": "pure8_s0m1_d3",
        "staircase_delta": -1,
        "dense_rows": 3,
        "gf256_rows": 8,
        "gf16_rows": 0,
        "period": 244,
        "x_mode": "frozen",
    },
    {
        "name": "pure8_s0_d3",
        "staircase_delta": 0,
        "dense_rows": 3,
        "gf256_rows": 8,
        "gf16_rows": 0,
        "period": 244,
        "x_mode": "frozen",
    },
    {
        "name": "pure8_s0m1_d5",
        "staircase_delta": -1,
        "dense_rows": 5,
        "gf256_rows": 8,
        "gf16_rows": 0,
        "period": 244,
        "x_mode": "frozen",
    },
    {
        "name": "pure9_s0m1_d3",
        "staircase_delta": -1,
        "dense_rows": 3,
        "gf256_rows": 9,
        "gf16_rows": 0,
        "period": 244,
        "x_mode": "frozen",
    },
)
CANDIDATE_BY_NAME = {item["name"]: item for item in CANDIDATES}


class CampaignError(RuntimeError):
    """The frozen campaign contract or its evidence was violated."""


def _is_frozen_candidate_name(value):
    return type(value) is str and value in CANDIDATE_BY_NAME


def canonical_json_bytes(value):
    """Return the one accepted JSON representation for campaign evidence."""
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ) + "\n"
    ).encode("ascii")


def sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def canonical_sha256(value):
    return sha256_bytes(canonical_json_bytes(value))


def _fsync_directory(directory):
    descriptor = os.open(str(directory), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _commit_no_replace(temporary, path):
    """Atomically publish ``temporary`` without ever replacing ``path``."""
    try:
        os.link(temporary, path, follow_symlinks=False)
    except FileExistsError:
        raise CampaignError(f"refusing to replace existing file {path}")
    os.unlink(temporary)
    _fsync_directory(path.parent)


def atomic_json(path, value):
    """Create one canonical JSON artifact atomically, refusing replacement."""
    path = Path(path)
    payload = canonical_json_bytes(value)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with open(temporary, "xb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        _commit_no_replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return sha256_bytes(payload)


def atomic_gzip_json(path, value):
    """Create deterministic compact gzip JSON atomically."""
    path = Path(path)
    payload = canonical_json_bytes(value)
    compressed = gzip.compress(payload, compresslevel=6, mtime=0)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with open(temporary, "xb") as output:
            output.write(compressed)
            output.flush()
            os.fsync(output.fileno())
        _commit_no_replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return {
        "compressed_sha256": sha256_bytes(compressed),
        "uncompressed_sha256": sha256_bytes(payload),
        "compressed_bytes": len(compressed),
        "uncompressed_bytes": len(payload),
    }


def _strict_json_payload(payload, label):
    try:
        value = peel_codec.strict_json_loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise CampaignError(f"{label} is not strict JSON: {error}")
    if canonical_json_bytes(value) != payload:
        raise CampaignError(f"{label} is not canonical JSON")
    return value


def read_canonical_json(path):
    path = Path(path)
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise CampaignError(f"could not read {path}: {error}")
    return _strict_json_payload(payload, str(path)), sha256_bytes(payload)


def read_canonical_gzip_json(path):
    path = Path(path)
    try:
        compressed = path.read_bytes()
        payload = gzip.decompress(compressed)
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as error:
        raise CampaignError(f"could not read gzip evidence {path}: {error}")
    return (
        _strict_json_payload(payload, str(path)),
        {
            "compressed_sha256": sha256_bytes(compressed),
            "uncompressed_sha256": sha256_bytes(payload),
            "compressed_bytes": len(compressed),
            "uncompressed_bytes": len(payload),
        },
    )


def construction_seed(base, replicate):
    return (
        base ^ (((replicate + 1) * 0xD192ED03) & 0xFFFFFFFF)
    ) & 0xFFFFFFFF


def loss_seed(base, replicate):
    return (
        base ^
        (((replicate + 1) * 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF)
    ) & 0xFFFFFFFFFFFFFFFF


def _seed_set(base, replicates, derivation):
    return {derivation(base, replicate) for replicate in range(replicates)}


def seed_disjointness_proof():
    """Exhaustively prove the frozen selection and holdout seed sets differ."""
    selection_construction = _seed_set(
        SELECTION_CONSTRUCTION_BASE,
        SELECTION_REPLICATES,
        construction_seed,
    )
    holdout_construction = _seed_set(
        HOLDOUT_CONSTRUCTION_BASE,
        HOLDOUT_REPLICATES,
        construction_seed,
    )
    selection_loss = _seed_set(
        SELECTION_LOSS_BASE,
        SELECTION_REPLICATES,
        loss_seed,
    )
    holdout_loss = _seed_set(
        HOLDOUT_LOSS_BASE,
        HOLDOUT_REPLICATES,
        loss_seed,
    )
    if (
        len(selection_construction) != SELECTION_REPLICATES
        or len(holdout_construction) != HOLDOUT_REPLICATES
        or len(selection_loss) != SELECTION_REPLICATES
        or len(holdout_loss) != HOLDOUT_REPLICATES
        or selection_construction & holdout_construction
        or selection_loss & holdout_loss
    ):
        raise CampaignError("selection and holdout seed sets are not disjoint")

    def digest(values):
        return canonical_sha256(sorted(values))

    return {
        "construction": {
            "derivation":
                band.BANDTIMING_CONSTRUCTION_SEED_DERIVATION,
            "selection_base": SELECTION_CONSTRUCTION_BASE,
            "selection_count": len(selection_construction),
            "selection_set_sha256": digest(selection_construction),
            "holdout_base": HOLDOUT_CONSTRUCTION_BASE,
            "holdout_count": len(holdout_construction),
            "holdout_set_sha256": digest(holdout_construction),
            "intersection_count": 0,
        },
        "loss": {
            "derivation": band.BANDTIMING_LOSS_SEED_DERIVATION,
            "selection_base": SELECTION_LOSS_BASE,
            "selection_count": len(selection_loss),
            "selection_set_sha256": digest(selection_loss),
            "holdout_base": HOLDOUT_LOSS_BASE,
            "holdout_count": len(holdout_loss),
            "holdout_set_sha256": digest(holdout_loss),
            "intersection_count": 0,
        },
    }


def candidate_descriptor(name, block_count):
    if (
        not _is_frozen_candidate_name(name)
        or type(block_count) is not int
        or block_count not in K_VALUES
    ):
        raise CampaignError(f"unknown frozen candidate {name!r}")
    try:
        spec = CANDIDATE_BY_NAME[name]
    except KeyError:
        raise CampaignError(f"unknown frozen candidate {name!r}")
    staircase = (
        band.dispatch_band_descriptor(block_count).staircase
        + spec["staircase_delta"]
    )
    return band.BandDescriptor(
        staircase=staircase,
        dense_rows=spec["dense_rows"],
        gf256_rows=spec["gf256_rows"],
        gf16_rows=spec["gf16_rows"],
        period=spec["period"],
        x_mode=spec["x_mode"],
    )


def _phase_parameters(phase):
    if phase == "selection":
        return {
            "replicates": SELECTION_REPLICATES,
            "construction_seed_base": SELECTION_CONSTRUCTION_BASE,
            "loss_seed_base": SELECTION_LOSS_BASE,
            "order_seed": SELECTION_ORDER_SEED,
        }
    if phase == "holdout":
        return {
            "replicates": HOLDOUT_REPLICATES,
            "construction_seed_base": HOLDOUT_CONSTRUCTION_BASE,
            "loss_seed_base": HOLDOUT_LOSS_BASE,
            "order_seed": HOLDOUT_ORDER_SEED,
        }
    raise CampaignError(f"unknown campaign phase {phase!r}")


def _job_without_identity(phase, candidate_name, block_count, schedule):
    parameters = _phase_parameters(phase)
    descriptor = candidate_descriptor(candidate_name, block_count)
    return {
        "schema": JOB_SCHEMA,
        "phase": phase,
        "candidate": candidate_name,
        "K": block_count,
        "schedule": schedule,
        "candidate_descriptor": descriptor.as_dict(),
        "dispatch_descriptor":
            band.dispatch_band_descriptor(block_count).as_dict(),
        "block_bytes": BLOCK_BYTES,
        "loss": LOSS,
        "construction_seed_base":
            parameters["construction_seed_base"],
        "loss_seed_base": parameters["loss_seed_base"],
        "warmup_replicates": WARMUP_REPLICATES,
        "replicates": parameters["replicates"],
        "inner_reps": INNER_REPS,
        "max_overhead": MAX_OVERHEAD,
        "cache_state": CACHE_STATE,
        "systematic_cache": SYSTEMATIC_CACHE,
        "evict_bytes": EVICT_BYTES,
        "required_margin": REQUIRED_MARGIN,
        "thermal_sampling_interval_ms": THERMAL_INTERVAL_MS,
    }


def build_pre_cpu_jobs(phase, survivor=None):
    """Return the immutable, deterministically ordered logical job list."""
    if phase == "selection":
        if survivor is not None:
            raise CampaignError("selection does not accept a survivor")
        candidate_names = tuple(item["name"] for item in CANDIDATES)
    elif phase == "holdout":
        if not _is_frozen_candidate_name(survivor):
            raise CampaignError(
                "holdout requires one explicit survivor from the frozen roster"
            )
        candidate_names = (survivor,)
    else:
        raise CampaignError(f"unknown campaign phase {phase!r}")

    jobs = []
    for candidate_name in candidate_names:
        for block_count in K_VALUES:
            _parent_signal_safe_point()
            for schedule in SCHEDULES:
                base = _job_without_identity(
                    phase, candidate_name, block_count, schedule)
                jobs.append({
                    **base,
                    "job_id": canonical_sha256(base),
                })
    if len({job["job_id"] for job in jobs}) != len(jobs):
        raise CampaignError("frozen job identities are not unique")
    order_seed = _phase_parameters(phase)["order_seed"]

    def order_key(job):
        seed = order_seed.to_bytes(8, "big")
        return hashlib.sha256(seed + job["job_id"].encode("ascii")).digest()

    jobs.sort(key=lambda job: (order_key(job), job["job_id"]))
    return tuple({**job, "order": order} for order, job in enumerate(jobs))


def pre_cpu_job_list_hash(phase, survivor=None):
    return canonical_sha256(list(build_pre_cpu_jobs(phase, survivor)))


def selection_policy():
    """Return the frozen, result-free selection and promotion rule."""
    timing_mean = (
        "arithmetic-mean-of-per-K-schedule-job-log_cost_mean-with-equal-"
        "job-weight"
    )
    aggregate_upper = (
        "max(two-sided-Student-t-95-percent-upper-bound-of-equal-job-"
        "log-cost-means,arithmetic-mean-of-authenticated-per-job-"
        "log-cost-CI-highs)"
    )
    return {
        "schema": "wirehair.wh2.za5v.selection-policy.v1",
        "raw_evidence": {
            "seed_treatment": "no-retries-no-fixes-no-repairs",
            "candidate_population": "all-K2-to-100-all-six-schedules",
            "recovery_weighting": "equal-weight-per-matched-raw-cell",
            "timing_weighting": timing_mean,
            "timing_completeness": {
                "cross_contrast":
                    "every-equal-weight-job-must-have-valid-finite-CI-"
                    "effective-floor-and-boolean-left-faster-evidence",
                "both_corresponding_aa_arms":
                    "every-equal-weight-job-must-have-valid-finite-AA-CI-"
                    "and-floor-evidence",
                "minimum_eligible_replicates_per_job": {
                    "selection": SELECTION_REPLICATES // 2,
                    "holdout": HOLDOUT_REPLICATES // 2,
                },
                "coverage_applies_to":
                    "cross-contrast-and-each-corresponding-AA-arm",
            },
        },
        "metrics": {
            "final_unrecovered_at_h64":
                "sum-final-unrecovered-over-all-recovery-groups",
            "overhead_tail_auc_h0_to_h64":
                "discrete-sum-of-unrecovered-counts-at-each-integer-h-0-to-64",
            "unique_weak_constructions":
                "count-unique-K-plus-construction-seed",
            "direct_candidate_dispatch_log_cost": timing_mean,
            "encoder_candidate_wh1_log_cost": timing_mean,
            "decoder_candidate_wh1_log_cost": timing_mean,
            "log_cost_sign": "negative-means-left-arm-faster",
        },
        "aggregate_timing_confirmation": {
            "job_mean_interval":
                "two-sided-Student-t-95-percent-interval-across-equal-"
                "weighted-job-log-cost-means",
            "aggregate_upper": aggregate_upper,
            "mean_effective_floor":
                "arithmetic-mean-of-authenticated-per-job-effective-"
                "floors-which-include-the-frozen-practical-margin-and-AA-"
                "noise",
            "pass_rule":
                "aggregate-upper-strictly-less-than-negative-mean-"
                "effective-floor",
            "per_job_left_faster":
                "reported-count-not-an-every-job-veto",
        },
        "raw_pareto_surface": {
            "dimensions": [
                "final_unrecovered_at_h64",
                "overhead_tail_auc_h0_to_h64",
                "unique_weak_constructions",
                "direct_candidate_dispatch_log_cost",
                "encoder_candidate_wh1_log_cost",
                "decoder_candidate_wh1_log_cost",
            ],
            "direction": "minimize-every-dimension",
            "dominance":
                "all-dimensions-less-than-or-equal-and-at-least-one-strict",
            "timing_incomplete_candidate":
                "excluded-from-pareto-and-selection",
        },
        "selection": {
            "eligibility":
                "candidate-final-unrecovered-and-tail-AUC-each-less-than-or-"
                "equal-to-both-matched-dispatch-and-matched-WH1",
            "ranking":
                "minimum-direct-candidate-dispatch-log-cost-on-eligible-"
                "raw-Pareto-surface",
            "tie_break": "lexicographically-smallest-candidate-name",
            "weak_seed_rule":
                "descriptive-Pareto-metric-not-an-eligibility-veto-repairs-"
                "deferred-until-after-architecture-selection",
            "no_eligible_candidate":
                "no-survivor-retain-dispatch-and-do-not-run-holdout",
        },
        "holdout": {
            "population": "predeclared-disjoint-seeds-selected-survivor-only",
            "recovery_confirmation":
                "same-final-unrecovered-and-tail-AUC-no-worse-than-both-"
                "matched-dispatch-and-matched-WH1",
            "direct_speed_confirmation": {
                "metric": "direct_candidate_dispatch_log_cost",
                "rule": aggregate_upper
                    + "-strictly-less-than-negative-mean-effective-floor",
            },
            "raw_architecture_rule":
                "constructor-weakness-and-candidate-specific-regressions-"
                "are-reported-but-do-not-erase-raw-architecture-selection",
            "failure_rule":
                "retain-dispatch-no-fallback-candidate-selection",
        },
        "production_reliability": {
            "candidate_control_panels": [
                panel for panel, unused_scope, unused_control
                in CANDIDATE_CONTROL_PANELS
            ],
            "repairable_constructor_codes": [3, 4],
            "repairable_rule":
                "any-candidate-weak-construction-in-selection-or-holdout-"
                "requires-seed-repair-and-independent-retest",
            "nonrepairable_rule":
                "any-candidate-error-outcome-even-shared-or-any-candidate-"
                "only-nonweak-regression-in-any-candidate-control-panel-"
                "requires-reliability-investigation",
            "shared_need_more_rule":
                "not-a-production-blocker-recovery-comparison-remains-"
                "authoritative",
            "cumulative_evidence":
                "authenticated-selected-candidate-selection-plus-disjoint-"
                "holdout",
        },
        "production_contract_promotion": {
            "requires_holdout_confirmation": True,
            "wh1_throughput_objectives": {
                "encoder_candidate_wh1_log_cost": aggregate_upper,
                "decoder_candidate_wh1_log_cost": aggregate_upper,
            },
            "requires_no_cumulative_reliability_blockers": True,
            "actions": {
                "promote-new-contract":
                    "all-recovery-speed-and-reliability-gates-pass",
                "repair-seeds-and-retest":
                    "raw-architecture-confirmed-but-known-candidate-weak-"
                    "construction-remains",
                "investigate-reliability":
                    "raw-architecture-confirmed-but-nonrepairable-blocker-"
                    "remains",
                "retain-dispatch":
                    "raw-architecture-recovery-or-speed-or-WH1-throughput-"
                    "gate-fails",
            },
            "failure_rule":
                "follow-prebound-action-routing-and-record-each-remaining-"
                "speed-or-reliability-gap",
        },
    }


def frozen_roster():
    return {
        "K_values": list(K_VALUES),
        "schedules": list(SCHEDULES),
        "candidates": [dict(item) for item in CANDIDATES],
        "dispatch": {
            "profile": peel_codec.TARGET_PROFILE,
            "seed_policy": peel_codec.TARGET_SEED_POLICY,
            "descriptor": "dispatch_band_descriptor(K)",
        },
        "request": {
            "block_bytes": BLOCK_BYTES,
            "loss": LOSS,
            "warmup_replicates": WARMUP_REPLICATES,
            "inner_reps": INNER_REPS,
            "max_overhead": MAX_OVERHEAD,
            "cache_state": CACHE_STATE,
            "systematic_cache": SYSTEMATIC_CACHE,
            "evict_bytes": EVICT_BYTES,
            "required_margin": REQUIRED_MARGIN,
            "thermal_sampling_interval_ms": THERMAL_INTERVAL_MS,
        },
        "phases": {
            "selection": _phase_parameters("selection"),
            "holdout": _phase_parameters("holdout"),
        },
        "selection_policy": selection_policy(),
    }


def frozen_roster_sha256():
    return canonical_sha256(frozen_roster())


def holdout_hashes():
    return {
        name: pre_cpu_job_list_hash("holdout", name)
        for name in CANDIDATE_BY_NAME
    }


def validate_frozen_contract():
    """Validate every native request without running the benchmark."""
    proof = seed_disjointness_proof()
    selection = build_pre_cpu_jobs("selection")
    if len(selection) != 2376:
        raise CampaignError("selection population is not 2,376 jobs")
    matched = {}
    for job in selection:
        _parent_signal_safe_point()
        key = (
            job["K"],
            job["schedule"],
            job["construction_seed_base"],
            job["loss_seed_base"],
            job["replicates"],
        )
        matched[key] = matched.get(key, 0) + 1
        validate_job(job)
    if (
        len(matched) != len(K_VALUES) * len(SCHEDULES)
        or set(matched.values()) != {len(CANDIDATES)}
    ):
        raise CampaignError("candidate jobs do not share an identical sample")
    # Each matched cell is the Cartesian product of the frozen K/schedule
    # roster and this ordered list of paired derived seeds.  Since every job
    # above was checked to use the exact same bases and replicate count, this
    # compact contract proves equality across candidates without serializing
    # 1.8 million cells in the dry planner.
    paired_seeds = [
        {
            "replicate": replicate,
            "construction_seed": construction_seed(
                SELECTION_CONSTRUCTION_BASE, replicate),
            "loss_seed": loss_seed(SELECTION_LOSS_BASE, replicate),
        }
        for replicate in range(SELECTION_REPLICATES)
    ]
    expected_derived = (
        len(K_VALUES) * len(SCHEDULES) * SELECTION_REPLICATES)
    matched_cell_hash = canonical_sha256({
        "K_values": list(K_VALUES),
        "schedules": list(SCHEDULES),
        "paired_seeds": paired_seeds,
        "cells_per_candidate": expected_derived,
    })
    for survivor in CANDIDATE_BY_NAME:
        holdout = build_pre_cpu_jobs("holdout", survivor)
        if len(holdout) != 594:
            raise CampaignError("holdout population is not 594 jobs")
        for job in holdout:
            _parent_signal_safe_point()
            validate_job(job)
    return {
        "frozen_roster_sha256": frozen_roster_sha256(),
        "selection_jobs": len(selection),
        "selection_job_list_sha256":
            canonical_sha256(list(selection)),
        "holdout_jobs_per_survivor": 594,
        "holdout_job_list_sha256": holdout_hashes(),
        "matched_selection_cells_per_candidate": expected_derived,
        "matched_selection_cell_set_sha256": matched_cell_hash,
        "seed_disjointness": proof,
    }


def descriptor_from_job(job):
    descriptor = job.get("candidate_descriptor")
    if (
        not isinstance(descriptor, dict)
        or set(descriptor) != {"S", "D2", "gf256", "gf16", "P", "x"}
    ):
        raise CampaignError("job candidate descriptor is malformed")
    return band.BandDescriptor(
        staircase=descriptor["S"],
        dense_rows=descriptor["D2"],
        gf256_rows=descriptor["gf256"],
        gf16_rows=descriptor["gf16"],
        period=descriptor["P"],
        x_mode=descriptor["x"],
    )


def expected_request(job):
    return {
        "block_count": job["K"],
        "block_bytes": job["block_bytes"],
        "candidate": descriptor_from_job(job),
        "dispatch_profile": peel_codec.TARGET_PROFILE,
        "seed_policy": peel_codec.TARGET_SEED_POLICY,
        "construction_seed": job["construction_seed_base"],
        "loss": job["loss"],
        "loss_seed": job["loss_seed_base"],
        "schedule": job["schedule"],
        "warmup_replicates": job["warmup_replicates"],
        "replicates": job["replicates"],
        "inner_reps": job["inner_reps"],
        "max_overhead": job["max_overhead"],
        "cache_state": job["cache_state"],
        "systematic_cache": job["systematic_cache"],
        "evict_bytes": job["evict_bytes"],
        "required_margin": job["required_margin"],
    }


def validate_job(job):
    if not isinstance(job, dict):
        raise CampaignError("job is not a dictionary")
    phase = job.get("phase")
    candidate_name = job.get("candidate")
    block_count = job.get("K")
    schedule = job.get("schedule")
    if (
        phase not in ("selection", "holdout")
        or not _is_frozen_candidate_name(candidate_name)
        or block_count not in K_VALUES
        or schedule not in SCHEDULES
    ):
        raise CampaignError("job is outside the frozen roster")
    expected = _job_without_identity(
        phase, candidate_name, block_count, schedule)
    identity = canonical_sha256(expected)
    if job != {
        **expected,
        "job_id": identity,
        "order": job.get("order"),
    }:
        raise CampaignError("job differs from its frozen request")
    if type(job["order"]) is not int or job["order"] < 0:
        raise CampaignError("job order is invalid")
    request = expected_request(job)
    band.validate_bandtiming_dimensions(**request)
    return request


def _stable_file_binding(path, *, executable=False):
    path = Path(path).resolve()
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
        try:
            before = os.fstat(descriptor)
            digest = hashlib.sha256()
            payload_bytes = 0
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
                payload_bytes += len(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        path_stat = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError(f"could not bind runtime file {path}: {error}")
    stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
              "st_ctime_ns")
    if (
        any(getattr(before, name) != getattr(after, name) for name in stable)
        or any(
            getattr(after, name) != getattr(path_stat, name)
            for name in ("st_dev", "st_ino", "st_mode", "st_size")
        )
        or not stat.S_ISREG(after.st_mode)
        or payload_bytes != after.st_size
        or (executable and not os.access(path, os.X_OK))
    ):
        raise CampaignError(f"runtime file is unstable or invalid: {path}")
    return {
        "path": str(path),
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": after.st_mode,
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
        "sha256": digest.hexdigest(),
    }


def _thermal_binding(path):
    path = Path(path).resolve()
    try:
        descriptor = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignError(f"could not bind thermal source {path}: {error}")
    if not stat.S_ISREG(descriptor.st_mode):
        raise CampaignError(f"thermal source is not a regular file: {path}")
    return {
        "path": str(path),
        "device": descriptor.st_dev,
        "inode": descriptor.st_ino,
        "mode": descriptor.st_mode,
    }


def runtime_bindings(bench_path, thermal_source):
    return {
        "benchmark": _stable_file_binding(bench_path, executable=True),
        "parser": _stable_file_binding(Path(band.__file__)),
        "context_tool": _stable_file_binding(Path(peel_codec.__file__)),
        "runner": _stable_file_binding(Path(__file__)),
        "thermal": _thermal_binding(thermal_source),
    }


def _validate_runtime_bindings_schema(expected):
    file_fields = {
        "path", "device", "inode", "mode", "size",
        "mtime_ns", "ctime_ns", "sha256",
    }
    thermal_fields = {"path", "device", "inode", "mode"}
    if (
        not isinstance(expected, dict)
        or set(expected) != {
            "benchmark", "parser", "context_tool", "runner", "thermal"
        }
    ):
        raise CampaignError("runtime binding schema is invalid")
    for name in ("benchmark", "parser", "context_tool", "runner"):
        item = expected[name]
        if (
            not isinstance(item, dict)
            or set(item) != file_fields
            or not isinstance(item.get("path"), str)
            or not Path(item["path"]).is_absolute()
            or any(type(item.get(field)) is not int or item[field] < 0
                   for field in (
                       "device", "inode", "mode", "size",
                       "mtime_ns", "ctime_ns"))
            or not peel_codec._is_sha256(item.get("sha256"))
        ):
            raise CampaignError(
                f"runtime {name} binding schema is invalid")
    thermal = expected["thermal"]
    if (
        not isinstance(thermal, dict)
        or set(thermal) != thermal_fields
        or not isinstance(thermal.get("path"), str)
        or not Path(thermal["path"]).is_absolute()
        or any(type(thermal.get(field)) is not int or thermal[field] < 0
               for field in ("device", "inode", "mode"))
    ):
        raise CampaignError("runtime thermal binding schema is invalid")
    return expected


def check_runtime_bindings(expected, *, full_hash):
    _validate_runtime_bindings_schema(expected)
    thermal = _thermal_binding(expected["thermal"]["path"])
    if thermal != expected["thermal"]:
        raise CampaignError("thermal source identity changed")
    for name in ("benchmark", "parser", "context_tool", "runner"):
        item = expected[name]
        path = Path(item["path"])
        try:
            current = os.stat(path, follow_symlinks=False)
        except OSError as error:
            raise CampaignError(
                f"runtime binding disappeared for {name}: {error}")
        for field, attribute in (
            ("device", "st_dev"),
            ("inode", "st_ino"),
            ("mode", "st_mode"),
            ("size", "st_size"),
            ("mtime_ns", "st_mtime_ns"),
            ("ctime_ns", "st_ctime_ns"),
        ):
            if getattr(current, attribute) != item[field]:
                raise CampaignError(f"runtime binding changed for {name}")
        if full_hash:
            rebound = _stable_file_binding(
                path, executable=(name == "benchmark"))
            if rebound != item:
                raise CampaignError(f"runtime content changed for {name}")


def _require_normal_priority():
    if not hasattr(os, "getpriority"):
        raise CampaignError("normal-priority validation requires os.getpriority")
    priority = os.getpriority(os.PRIO_PROCESS, 0)
    if priority != 0:
        raise CampaignError(
            f"campaign workers must run at normal priority, got nice {priority}")


def _create_fresh_directory(directory):
    directory = Path(directory)
    if directory.exists():
        raise CampaignError(
            f"refusing existing campaign directory {directory}")
    try:
        directory.mkdir(parents=True)
    except FileExistsError:
        raise CampaignError(
            f"refusing existing campaign directory {directory}")
    _fsync_directory(directory.parent)
    return directory


def _load_selection_contract(path, requested_survivor):
    _parent_signal_safe_point()
    path = Path(path).resolve()
    manifest, manifest_sha256 = read_canonical_json(path)
    _validate_manifest(manifest)
    if (
        manifest.get("schema") != CAMPAIGN_SCHEMA
        or manifest.get("phase") != "selection"
        or manifest.get("frozen_roster_sha256") != frozen_roster_sha256()
        or manifest.get("pre_cpu_job_list_sha256")
            != pre_cpu_job_list_hash("selection")
        or manifest.get("predeclared_holdout_job_list_sha256")
            != holdout_hashes()
    ):
        raise CampaignError("selection manifest is not the frozen contract")
    completion_path = path.parent / "completion.json"
    completion, completion_sha256 = read_canonical_json(completion_path)
    if not isinstance(completion, dict):
        raise CampaignError("selection completion is not an object")
    checksum_path = path.parent / "completion.sha256"
    try:
        checksum = checksum_path.read_text(encoding="ascii")
    except OSError as error:
        raise CampaignError(f"could not read {checksum_path}: {error}")
    expected_checksum = f"{completion_sha256}  completion.json\n"
    if (
        checksum != expected_checksum
        or completion.get("schema") != COMPLETION_SCHEMA
        or completion.get("phase") != "selection"
        or completion.get("manifest_sha256") != manifest_sha256
        or completion.get("jobs") != 2376
    ):
        raise CampaignError("selection campaign is not complete")
    verified = verify_campaign(path.parent)
    _parent_signal_safe_point()
    if (
        verified["manifest_sha256"] != manifest_sha256
        or verified["completion_sha256"] != completion_sha256
        or verified["jobs"] != 2376
    ):
        raise CampaignError(
            "selection campaign failed strict parent verification")
    summary, summary_evidence = read_canonical_gzip_json(
        path.parent / "summary.json.gz")
    if not isinstance(summary, dict):
        raise CampaignError("selection summary is not an object")
    if completion.get("summary") != {
        "path": "summary.json.gz",
        **summary_evidence,
    }:
        raise CampaignError(
            "selection decision is not bound by completion summary hash")
    decision = summary.get("decision")
    selected = _require_selected_survivor(
        decision, requested_survivor)
    blockers = _validate_production_blockers(
        decision["selected_candidate_production_blockers"], selected)
    contract = {
        "path": str(path),
        "manifest_sha256": manifest_sha256,
        "completion_sha256": completion_sha256,
        "selection_decision_sha256": canonical_sha256(decision),
        "selected_survivor": selected,
        "selected_candidate_production_blockers": blockers,
    }
    return contract


def make_plan(phase, survivor=None):
    jobs = build_pre_cpu_jobs(phase, survivor)
    validation = validate_frozen_contract()
    return {
        "schema": CAMPAIGN_SCHEMA,
        "phase": phase,
        "survivor": survivor,
        "frozen_roster": frozen_roster(),
        "frozen_roster_sha256": validation["frozen_roster_sha256"],
        "pre_cpu_job_count": len(jobs),
        "pre_cpu_job_list_sha256": canonical_sha256(list(jobs)),
        "predeclared_holdout_job_list_sha256":
            validation["holdout_job_list_sha256"],
        "matched_selection_cells_per_candidate":
            validation["matched_selection_cells_per_candidate"],
        "matched_selection_cell_set_sha256":
            validation["matched_selection_cell_set_sha256"],
        "seed_disjointness": validation["seed_disjointness"],
    }


def _build_manifest(
        phase, survivor, jobs, cpus, workers, bindings,
        selection_contract):
    assignments = []
    for index, job in enumerate(jobs):
        assignments.append({
            "order": index,
            "job_id": job["job_id"],
            "cpu": cpus[index % workers],
            "job_file": f"jobs/{index:04d}.json",
            "output": f"results/{index:04d}.json.gz",
            "log": f"logs/{index:04d}.log",
        })
    return {
        **make_plan(phase, survivor),
        "created_unix_ns": time.time_ns(),
        "repository": str(REPOSITORY),
        "workers": workers,
        "worker_cpus": cpus[:workers],
        "priority": "normal-nice-0",
        "runtime_bindings": bindings,
        "selection_contract": selection_contract,
        "pre_cpu_jobs": list(jobs),
        "assignments": assignments,
    }


def _write_job_files(directory, manifest, manifest_sha256):
    files = []
    for assignment, job in zip(
            manifest["assignments"], manifest["pre_cpu_jobs"]):
        _parent_signal_safe_point()
        record = _expected_job_record(
            manifest, manifest_sha256, assignment, job)
        path = directory / assignment["job_file"]
        digest = atomic_json(path, record)
        files.append({
            "path": assignment["job_file"],
            "sha256": digest,
        })
    return files


def _expected_job_record(manifest, manifest_sha256, assignment, job):
    return {
        "schema": JOB_SCHEMA,
        "manifest_sha256": manifest_sha256,
        "runtime_bindings": manifest["runtime_bindings"],
        "job": job,
        "assignment": assignment,
    }


def _verify_job_files(directory, manifest, manifest_sha256):
    evidence = []
    for assignment, job in zip(
            manifest["assignments"], manifest["pre_cpu_jobs"]):
        _parent_signal_safe_point()
        path = directory / assignment["job_file"]
        record, digest = read_canonical_json(path)
        if record != _expected_job_record(
                manifest, manifest_sha256, assignment, job):
            raise CampaignError(f"job file changed: {path}")
        evidence.append({
            "path": assignment["job_file"],
            "sha256": digest,
        })
    return evidence


def _verify_empty_logs(directory, manifest):
    evidence = []
    for assignment in manifest["assignments"]:
        _parent_signal_safe_point()
        path = directory / assignment["log"]
        try:
            payload = path.read_bytes()
        except OSError as error:
            raise CampaignError(f"could not read worker log {path}: {error}")
        if payload:
            raise CampaignError(f"successful worker log is not empty: {path}")
        evidence.append({
            "path": assignment["log"],
            "sha256": sha256_bytes(payload),
            "bytes": len(payload),
        })
    return evidence


def _validate_worker_record(record):
    if (
        not isinstance(record, dict)
        or set(record) != {
            "schema", "manifest_sha256", "runtime_bindings",
            "job", "assignment",
        }
        or record["schema"] != JOB_SCHEMA
        or not peel_codec._is_sha256(record["manifest_sha256"])
    ):
        raise CampaignError("worker job record is malformed")
    _validate_runtime_bindings_schema(record["runtime_bindings"])
    job = record["job"]
    validate_job(job)
    assignment = record["assignment"]
    if (
        not isinstance(assignment, dict)
        or set(assignment) != {
            "order", "job_id", "cpu", "job_file", "output", "log"
        }
        or assignment["order"] != job["order"]
        or assignment["job_id"] != job["job_id"]
        or type(assignment["cpu"]) is not int
        or assignment["cpu"] < 0
        or any(
            type(assignment[name]) is not str
            for name in ("job_file", "output", "log")
        )
    ):
        raise CampaignError("worker CPU assignment is malformed")
    return job, assignment


def run_worker(job_file, output_path):
    """Run and strictly replay one frozen native job."""
    _require_normal_priority()
    record, unused_sha = read_canonical_json(job_file)
    job, assignment = _validate_worker_record(record)
    job_file = Path(job_file).resolve()
    output_path = Path(output_path).resolve()
    campaign_directory = job_file.parent.parent
    if (
        job_file != (campaign_directory / assignment["job_file"]).resolve()
        or output_path
            != (campaign_directory / assignment["output"]).resolve()
    ):
        raise CampaignError("worker output path does not match assignment")
    try:
        os.sched_setaffinity(0, {assignment["cpu"]})
    except (AttributeError, OSError) as error:
        raise CampaignError(f"could not pin worker CPU: {error}")
    affinity = sorted(os.sched_getaffinity(0))
    if affinity != [assignment["cpu"]]:
        raise CampaignError("worker CPU affinity does not match assignment")

    bindings = record["runtime_bindings"]
    before_bindings = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if before_bindings != bindings:
        raise CampaignError("worker runtime binding changed before launch")
    request = expected_request(job)
    context = peel_codec.make_paired_context_config(
        bindings["thermal"]["path"],
        job["thermal_sampling_interval_ms"],
        cache_state=request["cache_state"],
        evict_bytes=request["evict_bytes"],
    )
    wall_started_ns = time.time_ns()
    measurement = band.bandtiming_probe(
        bindings["benchmark"]["path"],
        request["block_count"],
        request["block_bytes"],
        request["candidate"],
        dispatch_profile=request["dispatch_profile"],
        seed_policy=request["seed_policy"],
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
    wall_finished_ns = time.time_ns()
    receipt = measurement.as_dict()
    band.replay_bandtiming_receipt(
        receipt, expected_request=request)
    bound = receipt["context"]["bound"]
    if (
        bound["cpu_affinity"] != [assignment["cpu"]]
        or bound["thermal_device"] != bindings["thermal"]["device"]
        or bound["thermal_inode"] != bindings["thermal"]["inode"]
    ):
        raise CampaignError(
            "native receipt does not bind the assigned CPU and thermal source")
    after_bindings = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if after_bindings != before_bindings:
        raise CampaignError("worker runtime binding changed during probe")
    envelope = {
        "schema": RESULT_SCHEMA,
        "manifest_sha256": record["manifest_sha256"],
        "job": job,
        "assignment": assignment,
        "wall_started_unix_ns": wall_started_ns,
        "wall_finished_unix_ns": wall_finished_ns,
        "runtime_bindings_before": before_bindings,
        "runtime_bindings_after": after_bindings,
        "receipt": receipt,
    }
    atomic_gzip_json(output_path, envelope)


def _kill_process_groups(processes):
    _kill_process_groups_with_grace(processes, 5.0)


def _process_group_exists(process_group):
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    return True


def _kill_process_groups_with_grace(processes, grace_seconds):
    """Terminate every recorded group, even after its leader has exited."""
    _kill_process_groups_with_grace_impl(processes, grace_seconds)


def _kill_process_groups_with_grace_impl(processes, grace_seconds):
    processes = tuple(processes)
    process_groups = tuple(sorted({
        process.pid for process, unused_log in processes
    }))
    for process_group in process_groups:
        try:
            os.killpg(process_group, signal.SIGTERM)
        except ProcessLookupError:
            pass

    deadline = time.monotonic() + grace_seconds
    surviving = set(process_groups)
    while surviving:
        # ``poll`` reaps leaders that exited after TERM, but group liveness is
        # checked independently because an orphaned descendant may remain.
        for process, unused_log in processes:
            process.poll()
        surviving = {
            process_group for process_group in surviving
            if _process_group_exists(process_group)
        }
        remaining = deadline - time.monotonic()
        if not surviving or remaining <= 0.0:
            break
        time.sleep(min(0.01, remaining))

    for process_group in sorted(surviving):
        try:
            os.killpg(process_group, signal.SIGKILL)
        except ProcessLookupError:
            pass

    # Reap leaders separately from process-group cleanup.  In particular, an
    # already reaped leader says nothing about descendants in the same pgid.
    for process, unused_log in processes:
        try:
            process.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()

    confirmation_deadline = time.monotonic() + max(1.0, grace_seconds)
    surviving = set(process_groups)
    while surviving:
        surviving = {
            process_group for process_group in surviving
            if _process_group_exists(process_group)
        }
        remaining = confirmation_deadline - time.monotonic()
        if not surviving or remaining <= 0.0:
            break
        time.sleep(min(0.01, remaining))
    if surviving:
        raise CampaignError(
            "worker process groups survived SIGKILL: "
            + ",".join(str(item) for item in sorted(surviving)))


def _register_process(processes, process, log):
    """Register a launched process group before any signal safe point."""
    processes.append((process, log))


def _run_wave(directory, manifest, assignments):
    processes = []
    script = Path(__file__).resolve()
    try:
        for assignment in assignments:
            # TERM/HUP handlers only record intent.  Check immediately before
            # launch and after the returned process group is registered, so
            # there is no asynchronous exception window that can orphan it.
            _parent_signal_safe_point()
            log = None
            process = None
            appended = False
            try:
                log_path = directory / assignment["log"]
                log = open(log_path, "xb")
                command = [
                    sys.executable,
                    str(script),
                    "_worker",
                    "--job-file",
                    str(directory / assignment["job_file"]),
                    "--output",
                    str(directory / assignment["output"]),
                ]
                process = subprocess.Popen(
                    command,
                    cwd=str(REPOSITORY),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                _register_process(processes, process, log)
                appended = True
            finally:
                try:
                    if process is not None and not appended:
                        _kill_process_groups(((process, log),))
                finally:
                    if not appended and log is not None:
                        log.close()
            _parent_signal_safe_point()

        checks = 0
        checks_while_running = 0
        first_check_ns = None
        last_check_ns = None
        failures = []
        while True:
            _parent_signal_safe_point()
            check_runtime_bindings(
                manifest["runtime_bindings"], full_hash=False)
            checked_ns = time.monotonic_ns()
            if first_check_ns is None:
                first_check_ns = checked_ns
            last_check_ns = checked_ns
            checks += 1
            running = False
            for process, unused_log in processes:
                status = process.poll()
                if status is None:
                    running = True
                elif status != 0 and process.pid not in {
                        item[0] for item in failures}:
                    failures.append((process.pid, status))
            if running:
                checks_while_running += 1
            if failures:
                raise CampaignError(f"campaign wave failed: {failures}")
            if not running:
                break
            _parent_signal_safe_point()
            time.sleep(0.25)
        check_runtime_bindings(
            manifest["runtime_bindings"], full_hash=True)
        _parent_signal_safe_point()
        return {
            "first_order": assignments[0]["order"],
            "last_order": assignments[-1]["order"],
            "checks": checks,
            "checks_while_running": checks_while_running,
            "first_check_monotonic_ns": first_check_ns,
            "last_check_monotonic_ns": last_check_ns,
            "runtime_bindings_sha256":
                canonical_sha256(manifest["runtime_bindings"]),
        }
    except BaseException:
        _kill_process_groups(processes)
        raise
    finally:
        for unused_process, log in processes:
            if not log.closed:
                log.close()


def _new_group(job):
    def arm():
        return {
            "result_classes": {name: 0 for name in RESULT_CLASSES},
            "byte_recovery_failures": 0,
            "overhead_successes": 0,
            "overhead_exact": [0] * (MAX_OVERHEAD + 1),
        }

    def pair(left, right):
        return {
            "left_arm": left,
            "right_arm": right,
            "final": {
                "both_recovered": 0,
                "left_only_recovered": 0,
                "right_only_recovered": 0,
                "neither_recovered": 0,
            },
            "_threshold_deltas": {
                "both_recovered": [0] * (MAX_OVERHEAD + 2),
                "left_only_recovered": [0] * (MAX_OVERHEAD + 2),
                "right_only_recovered": [0] * (MAX_OVERHEAD + 2),
                "neither_recovered": [0] * (MAX_OVERHEAD + 2),
            },
        }

    return {
        "candidate_name": job["candidate"],
        "K": job["K"],
        "schedule": job["schedule"],
        "loss": job["loss"],
        "raw_cells": 0,
        "recovery_categories": {
            "both_success": 0,
            "candidate_only_failure": 0,
            "dispatch_only_failure": 0,
            "shared_failure": 0,
        },
        "candidate": arm(),
        "dispatch": arm(),
        "wh1": arm(),
        "censored_replicates": 0,
        "censored_panel_observations": {
            name: 0 for name in band.BANDTIMING_PANELS
        },
        "fault_contaminated_replicates": 0,
        "fault_contaminated_panel_observations": {
            name: 0 for name in band.BANDTIMING_PANELS
        },
        "_paired_recovery": {
            f"{left}_vs_{right}": pair(left, right)
            for left, right in RECOVERY_PAIRS
        },
    }


class RecoveryAggregator:
    """Aggregate matched raw cells without collapsing failures or censoring."""

    def __init__(self, phase, jobs):
        self.phase = phase
        self.jobs = tuple(jobs)
        self.groups = {}
        self.seen_cells = set()
        tested_candidates = {job["candidate"] for job in self.jobs}
        if (
            (phase == "selection"
             and tested_candidates != set(CANDIDATE_BY_NAME))
            or (phase == "holdout" and len(tested_candidates) != 1)
        ):
            raise CampaignError(
                "aggregator candidate population does not match its phase")
        self.candidate_weak = {
            name: set() for name in sorted(tested_candidates)
        }
        self.dispatch_weak = set()
        self.candidate_constructor_classes = {}
        self.candidate_constructor_schedules = {}
        self.dispatch_constructor_classes = {}
        self.dispatch_constructor_schedules = {}
        self.wh1_constructor_classes = {}
        self.wh1_constructor_schedules = {}
        self.wh1_weak = set()
        self.control_signatures = {}
        self.control_candidates = {}

    @staticmethod
    def _validate_class(value, label):
        if value not in RESULT_CLASSES:
            raise CampaignError(f"{label} has invalid result class {value!r}")

    @classmethod
    def _construction_result(cls, replicate, arm):
        """Consume the parser-authenticated constructor result, never infer it."""
        result = replicate[f"{arm}_construct_result"]
        outcome = replicate[f"{arm}_construct_class"]
        cls._validate_class(outcome, f"{arm} construction")
        if (
            type(result) is not int
            or result == 9
            or not 0 <= result <= 10
        ):
            raise CampaignError(
                f"{arm} construction result code is invalid")
        try:
            expected = band._generic_result_class(result)
        except band.MeasurementError as error:
            raise CampaignError(
                f"{arm} construction result code is invalid") from error
        if outcome != expected:
            raise CampaignError(
                f"{arm} construction class contradicts its result code")
        return result, outcome

    @staticmethod
    def _update_arm(arm, outcome, overhead):
        arm["result_classes"][outcome] += 1
        if outcome == "success":
            if type(overhead) is not int or not 0 <= overhead <= MAX_OVERHEAD:
                raise CampaignError("successful recovery overhead is invalid")
            arm["overhead_successes"] += 1
            arm["overhead_exact"][overhead] += 1
        else:
            if overhead != -1:
                raise CampaignError("failed recovery reported an overhead")
            arm["byte_recovery_failures"] += 1

    @staticmethod
    def _increment_range(deltas, start, finish):
        """Add one to inclusive threshold range [start, finish]."""
        if start > finish:
            return
        deltas[start] += 1
        deltas[finish + 1] -= 1

    @classmethod
    def _update_paired_recovery(
            cls, pair, left_class, left_overhead,
            right_class, right_overhead):
        left_ok = left_class == "success"
        right_ok = right_class == "success"
        if left_ok and right_ok:
            final_category = "both_recovered"
        elif left_ok:
            final_category = "left_only_recovered"
        elif right_ok:
            final_category = "right_only_recovered"
        else:
            final_category = "neither_recovered"
        pair["final"][final_category] += 1

        left_at = (
            left_overhead
            if left_ok else MAX_OVERHEAD + 1
        )
        right_at = (
            right_overhead
            if right_ok else MAX_OVERHEAD + 1
        )
        deltas = pair["_threshold_deltas"]
        earliest = min(left_at, right_at)
        cls._increment_range(
            deltas["neither_recovered"], 0, earliest - 1)
        if left_at < right_at:
            cls._increment_range(
                deltas["left_only_recovered"],
                left_at,
                min(right_at - 1, MAX_OVERHEAD),
            )
        elif right_at < left_at:
            cls._increment_range(
                deltas["right_only_recovered"],
                right_at,
                min(left_at - 1, MAX_OVERHEAD),
            )
        latest = max(left_at, right_at)
        if latest <= MAX_OVERHEAD:
            cls._increment_range(
                deltas["both_recovered"], latest, MAX_OVERHEAD)

    def add_replicate(self, job, replicate):
        required = {
            "replicate", "construction_seed", "loss_seed", "trace_sha256",
            "candidate_construct_result", "candidate_construct_class",
            "dispatch_construct_result", "dispatch_construct_class",
            "wh1_construct_result", "wh1_construct_class",
            "encoder_candidate_result_class",
            "encoder_candidate_overhead",
            "encoder_dispatch_result_class",
            "encoder_dispatch_overhead",
            "encoder_wh1_result_class", "encoder_wh1_overhead",
            "decoder_candidate_result_class",
            "decoder_candidate_overhead",
            "decoder_dispatch_result_class",
            "decoder_dispatch_overhead",
            "decoder_wh1_result_class", "decoder_wh1_overhead",
            "direct_candidate_result_class", "direct_candidate_overhead",
            "direct_dispatch_result_class", "direct_dispatch_overhead",
            "censored_panels", "fault_contaminated_panels",
        }
        if not isinstance(replicate, dict) or set(replicate) != required:
            raise CampaignError("replicate receipt schema is incomplete")
        index = replicate["replicate"]
        if (
            type(index) is not int
            or not 0 <= index < job["replicates"]
            or replicate["construction_seed"] != construction_seed(
                job["construction_seed_base"], index)
            or replicate["loss_seed"] != loss_seed(
                job["loss_seed_base"], index)
            or not peel_codec._is_sha256(replicate["trace_sha256"])
        ):
            raise CampaignError("replicate does not match its frozen seeds")

        for scope, arm in TIMING_ARM_KEYS:
            self._validate_class(
                replicate[f"{scope}_{arm}_result_class"],
                f"{scope}/{arm}",
            )
        candidate_class = replicate["decoder_candidate_result_class"]
        dispatch_class = replicate["decoder_dispatch_result_class"]
        wh1_class = replicate["decoder_wh1_result_class"]
        candidate_construct_result, candidate_constructor = \
            self._construction_result(replicate, "candidate")
        dispatch_construct_result, dispatch_constructor = \
            self._construction_result(replicate, "dispatch")
        wh1_construct_result, wh1_constructor = \
            self._construction_result(replicate, "wh1")
        candidate_constructor_key = (
            job["candidate"], job["K"], replicate["construction_seed"])
        dispatch_constructor_key = (
            job["K"], replicate["construction_seed"])
        candidate_constructor_signature = (
            candidate_construct_result, candidate_constructor)
        dispatch_constructor_signature = (
            dispatch_construct_result, dispatch_constructor)
        wh1_constructor_signature = (
            wh1_construct_result, wh1_constructor)
        prior = self.candidate_constructor_classes.setdefault(
            candidate_constructor_key, candidate_constructor_signature)
        if prior != candidate_constructor_signature:
            raise CampaignError(
                "candidate constructor result/class changed across schedules")
        self.candidate_constructor_schedules.setdefault(
            candidate_constructor_key, set()).add(job["schedule"])
        prior = self.dispatch_constructor_classes.setdefault(
            dispatch_constructor_key, dispatch_constructor_signature)
        if prior != dispatch_constructor_signature:
            raise CampaignError(
                "dispatch constructor result/class changed across duplicated "
                "jobs")
        self.dispatch_constructor_schedules.setdefault(
            dispatch_constructor_key, set()).add(job["schedule"])
        prior = self.wh1_constructor_classes.setdefault(
            dispatch_constructor_key, wh1_constructor_signature)
        if prior != wh1_constructor_signature:
            raise CampaignError(
                "WH1 constructor result/class changed across duplicated jobs")
        self.wh1_constructor_schedules.setdefault(
            dispatch_constructor_key, set()).add(job["schedule"])

        cell_key = (
            job["candidate"],
            job["K"],
            job["schedule"],
            replicate["construction_seed"],
            replicate["loss_seed"],
        )
        if cell_key in self.seen_cells:
            raise CampaignError("duplicate matched recovery cell")
        self.seen_cells.add(cell_key)
        control_key = (
            job["K"],
            job["schedule"],
            replicate["construction_seed"],
            replicate["loss_seed"],
        )
        control_signature = (
            replicate["trace_sha256"],
            dispatch_construct_result,
            dispatch_constructor,
            replicate["encoder_dispatch_result_class"],
            replicate["encoder_dispatch_overhead"],
            replicate["decoder_dispatch_result_class"],
            replicate["decoder_dispatch_overhead"],
            replicate["direct_dispatch_result_class"],
            wh1_construct_result,
            wh1_constructor,
            replicate["encoder_wh1_result_class"],
            replicate["encoder_wh1_overhead"],
            replicate["decoder_wh1_result_class"],
            replicate["decoder_wh1_overhead"],
        )
        # Direct candidate/dispatch observations intentionally share a
        # pair-local fixed prefix: max(candidate recovery, dispatch recovery).
        # Its recorded overhead therefore belongs to the candidate comparison,
        # not to the candidate-independent control signature.  The result
        # class remains invariant: a recovered dispatch must also solve at any
        # larger pair prefix, while an unrecovered dispatch always uses K+64.
        # The native receipt parser authenticates that pair-local relationship
        # before this aggregator consumes the replicate.
        prior = self.control_signatures.setdefault(
            control_key, control_signature)
        if prior != control_signature:
            raise CampaignError(
                "dispatch/WH1 control changed across candidate duplicates")
        candidates = self.control_candidates.setdefault(control_key, set())
        if job["candidate"] in candidates:
            raise CampaignError("duplicate candidate copy of a control cell")
        candidates.add(job["candidate"])
        group_key = (job["candidate"], job["K"], job["schedule"])
        group = self.groups.setdefault(group_key, _new_group(job))
        group["raw_cells"] += 1

        candidate_ok = candidate_class == "success"
        dispatch_ok = dispatch_class == "success"
        if candidate_ok and dispatch_ok:
            category = "both_success"
        elif not candidate_ok and dispatch_ok:
            category = "candidate_only_failure"
        elif candidate_ok and not dispatch_ok:
            category = "dispatch_only_failure"
        else:
            category = "shared_failure"
        group["recovery_categories"][category] += 1

        self._update_arm(
            group["candidate"],
            candidate_class,
            replicate["decoder_candidate_overhead"],
        )
        self._update_arm(
            group["dispatch"],
            dispatch_class,
            replicate["decoder_dispatch_overhead"],
        )
        self._update_arm(
            group["wh1"],
            wh1_class,
            replicate["decoder_wh1_overhead"],
        )
        recovery = {
            "candidate": (
                candidate_class,
                replicate["decoder_candidate_overhead"],
            ),
            "dispatch": (
                dispatch_class,
                replicate["decoder_dispatch_overhead"],
            ),
            "wh1": (
                wh1_class,
                replicate["decoder_wh1_overhead"],
            ),
        }
        for left, right in RECOVERY_PAIRS:
            self._update_paired_recovery(
                group["_paired_recovery"][f"{left}_vs_{right}"],
                *recovery[left],
                *recovery[right],
            )

        censored = replicate["censored_panels"]
        contaminated = replicate["fault_contaminated_panels"]
        for values, count_name, observation_name in (
            (
                censored,
                "censored_replicates",
                "censored_panel_observations",
            ),
            (
                contaminated,
                "fault_contaminated_replicates",
                "fault_contaminated_panel_observations",
            ),
        ):
            if (
                not isinstance(values, list)
                or any(
                    type(value) is not str
                    or value not in band.BANDTIMING_PANELS
                    for value in values
                )
                or len(values) != len(set(values))
            ):
                raise CampaignError("replicate panel censor list is invalid")
            if values:
                group[count_name] += 1
            for value in values:
                group[observation_name][value] += 1

        weak_key = (job["K"], replicate["construction_seed"])
        if candidate_constructor == "weak":
            self.candidate_weak[job["candidate"]].add(weak_key)
        if dispatch_constructor == "weak":
            self.dispatch_weak.add(weak_key)
        if wh1_constructor == "weak":
            self.wh1_weak.add(weak_key)

        return {
            "schema": CELL_LEDGER_SCHEMA,
            "phase": self.phase,
            "candidate": job["candidate"],
            "K": job["K"],
            "schedule": job["schedule"],
            "loss": job["loss"],
            "replicate": index,
            "construction_seed": replicate["construction_seed"],
            "loss_seed": replicate["loss_seed"],
            "trace_sha256": replicate["trace_sha256"],
            "candidate_construct_result": candidate_construct_result,
            "candidate_construction_class": candidate_constructor,
            "dispatch_construct_result": dispatch_construct_result,
            "dispatch_construction_class": dispatch_constructor,
            "wh1_construct_result": wh1_construct_result,
            "wh1_construction_class": wh1_constructor,
            "candidate_recovery_class": candidate_class,
            "dispatch_recovery_class": dispatch_class,
            "wh1_recovery_class": wh1_class,
            "candidate_overhead":
                replicate["decoder_candidate_overhead"],
            "dispatch_overhead":
                replicate["decoder_dispatch_overhead"],
            "wh1_overhead": replicate["decoder_wh1_overhead"],
            "recovery_category": category,
            "candidate_byte_recovery_failure": not candidate_ok,
            "dispatch_byte_recovery_failure": not dispatch_ok,
            "wh1_byte_recovery_failure": wh1_class != "success",
            "censored_panels": list(censored),
            "fault_contaminated_panels": list(contaminated),
        }

    @staticmethod
    def _finish_arm(arm, raw_cells):
        successes = arm["overhead_successes"]
        failures = raw_cells - successes
        if (
            sum(arm["result_classes"].values()) != raw_cells
            or sum(arm["overhead_exact"]) != successes
            or arm["byte_recovery_failures"] != failures
        ):
            raise CampaignError("recovery arm totals are inconsistent")
        tail = []
        running = failures
        suffix = sum(arm["overhead_exact"][1:])
        for overhead in range(MAX_OVERHEAD + 1):
            if overhead == 0:
                running = failures + suffix
            else:
                suffix -= arm["overhead_exact"][overhead]
                running = failures + suffix
            tail.append(running)
        return {
            **arm,
            "final_unrecovered": failures,
            "unrecovered_by_overhead_h_0_to_64": tail,
        }

    @staticmethod
    def _finish_pair(pair, raw_cells):
        if sum(pair["final"].values()) != raw_cells:
            raise CampaignError(
                "paired final recovery counts do not partition cells")
        threshold_counts = {
            name: []
            for name in pair["_threshold_deltas"]
        }
        for name, deltas in pair["_threshold_deltas"].items():
            running = 0
            for overhead in range(MAX_OVERHEAD + 1):
                running += deltas[overhead]
                threshold_counts[name].append(running)
        thresholds = []
        for overhead in range(MAX_OVERHEAD + 1):
            counts = {
                name: values[overhead]
                for name, values in threshold_counts.items()
            }
            if sum(counts.values()) != raw_cells:
                raise CampaignError(
                    "paired threshold counts do not partition cells")
            thresholds.append({
                "h": overhead,
                "denominator": raw_cells,
                **counts,
            })
        if {
            name: thresholds[-1][name]
            for name in pair["final"]
        } != pair["final"]:
            raise CampaignError(
                "paired final recovery counts disagree with h=64")
        return {
            "left_arm": pair["left_arm"],
            "right_arm": pair["right_arm"],
            "final": dict(pair["final"]),
            "by_overhead_h_0_to_64": thresholds,
        }

    def finish(self):
        expected_cells = sum(job["replicates"] for job in self.jobs)
        if len(self.seen_cells) != expected_cells:
            raise CampaignError(
                f"recovery ledger has {len(self.seen_cells)} cells; "
                f"expected {expected_cells}")
        expected_control_candidates = (
            set(CANDIDATE_BY_NAME)
            if self.phase == "selection"
            else {job["candidate"] for job in self.jobs}
        )
        if any(
                candidates != expected_control_candidates
                for candidates in self.control_candidates.values()
        ):
            raise CampaignError(
                "control cells lack their exact candidate copy population")
        expected_schedules = set(SCHEDULES)
        if any(
                schedules != expected_schedules
                for schedules in self.candidate_constructor_schedules.values()
        ):
            raise CampaignError(
                "candidate constructor class lacks all six schedule duplicates")
        if any(
                schedules != expected_schedules
                for schedules in self.dispatch_constructor_schedules.values()
        ):
            raise CampaignError(
                "dispatch constructor class lacks all six schedule duplicates")
        if any(
                schedules != expected_schedules
                for schedules in self.wh1_constructor_schedules.values()
        ):
            raise CampaignError(
                "WH1 constructor class lacks all six schedule duplicates")
        groups = []
        for key in sorted(self.groups):
            group = self.groups[key]
            if group["raw_cells"] != _phase_parameters(
                    self.phase)["replicates"]:
                raise CampaignError("recovery group has missing replicates")
            if sum(group["recovery_categories"].values()) != group["raw_cells"]:
                raise CampaignError("recovery categories do not partition cells")
            finished = dict(group)
            del finished["_paired_recovery"]
            finished["candidate"] = self._finish_arm(
                group["candidate"], group["raw_cells"])
            finished["dispatch"] = self._finish_arm(
                group["dispatch"], group["raw_cells"])
            finished["wh1"] = self._finish_arm(
                group["wh1"], group["raw_cells"])
            finished["paired_recovery"] = [
                self._finish_pair(
                    group["_paired_recovery"][f"{left}_vs_{right}"],
                    group["raw_cells"],
                )
                for left, right in RECOVERY_PAIRS
            ]
            candidate_dispatch = finished["paired_recovery"][0][
                "by_overhead_h_0_to_64"]
            finished[
                "paired_candidate_dispatch_by_overhead_h_0_to_64"
            ] = [
                {
                    "h": item["h"],
                    "denominator": item["denominator"],
                    "both_recovered": item["both_recovered"],
                    "candidate_only_unrecovered":
                        item["right_only_recovered"],
                    "dispatch_only_unrecovered":
                        item["left_only_recovered"],
                    "both_unrecovered": item["neither_recovered"],
                }
                for item in candidate_dispatch
            ]
            groups.append(finished)

        candidate_weak = {}
        dispatch_comparisons = {}
        wh1_comparisons = {}
        for name, values in self.candidate_weak.items():
            ordered = [
                {"K": block_count, "construction_seed": seed}
                for block_count, seed in sorted(values)
            ]
            candidate_weak[name] = {
                "count": len(ordered),
                "cells": ordered,
            }
            def weak_cells(cells):
                return [
                    {"K": block_count, "construction_seed": seed}
                    for block_count, seed in sorted(cells)
                ]

            def comparison(control):
                candidate_only = values - control
                shared = values & control
                control_only = control - values
                return {
                    "candidate_only": {
                        "count": len(candidate_only),
                        "cells": weak_cells(candidate_only),
                    },
                    "shared": {
                        "count": len(shared),
                        "cells": weak_cells(shared),
                    },
                    "control_only": {
                        "count": len(control_only),
                        "cells": weak_cells(control_only),
                    },
                }

            dispatch_comparisons[name] = comparison(self.dispatch_weak)
            wh1_comparisons[name] = comparison(self.wh1_weak)
        dispatch_ordered = [
            {"K": block_count, "construction_seed": seed}
            for block_count, seed in sorted(self.dispatch_weak)
        ]
        wh1_ordered = [
            {"K": block_count, "construction_seed": seed}
            for block_count, seed in sorted(self.wh1_weak)
        ]
        return {
            "raw_cells": len(self.seen_cells),
            "matched_control_cells": {
                "unique_cells": len(self.control_signatures),
                "copies_per_cell": len(expected_control_candidates),
                "candidate_copies": sorted(expected_control_candidates),
                "identity":
                    "K-schedule-construction-seed-loss-seed",
                "signature":
                    "trace-dispatch-construct-result-class-encoder-decoder-"
                    "result-class-overhead-direct-result-class-"
                    "wh1-construct-result-class-"
                    "encoder-decoder-result-class-overhead",
            },
            "groups": groups,
            "independent_construction_weaknesses": {
                "identity": "unique-K-plus-construction-seed",
                "candidate": candidate_weak,
                "dispatch": {
                    "count": len(dispatch_ordered),
                    "cells": dispatch_ordered,
                },
                "wh1": {
                    "count": len(wh1_ordered),
                    "cells": wh1_ordered,
                },
                "candidate_dispatch_intersections": dispatch_comparisons,
                "candidate_wh1_intersections": wh1_comparisons,
            },
        }


def _classify_candidate_reliability(contrasts, replicates):
    """Enrich authenticated panel totals from exact replicate result classes."""
    enriched = [dict(item) for item in contrasts]
    by_name = {item.get("name"): item for item in enriched}
    if len(by_name) != len(enriched):
        raise CampaignError("candidate reliability panels are duplicated")
    panel_counts = {
        panel: {
            "candidate_constructor_weak_regressions": 0,
            "candidate_nonrepairable_regressions": 0,
        }
        for panel, unused_scope, unused_control
        in CANDIDATE_CONTROL_PANELS
    }
    candidate_errors = {
        "constructor": 0,
        "encoder": 0,
        "decoder": 0,
        "direct": 0,
    }
    for replicate in replicates:
        if not isinstance(replicate, dict):
            raise CampaignError(
                "candidate reliability replicate is malformed")
        construct_result = replicate.get("candidate_construct_result")
        construct_class = replicate.get("candidate_construct_class")
        if (
            type(construct_result) is not int
            or construct_result < 0
            or construct_class not in RESULT_CLASSES
        ):
            raise CampaignError(
                "candidate constructor reliability evidence is malformed")
        if construct_class == "weak" and construct_result not in (3, 4):
            raise CampaignError(
                "only BadDense/BadPeel are repairable constructor weakness")
        if construct_class == "error":
            candidate_errors["constructor"] += 1

        for scope in ("encoder", "decoder", "direct"):
            candidate_class = replicate.get(
                f"{scope}_candidate_result_class")
            if candidate_class not in RESULT_CLASSES:
                raise CampaignError(
                    "candidate result reliability evidence is malformed")
            if candidate_class == "error":
                candidate_errors[scope] += 1

        for panel, scope, control in CANDIDATE_CONTROL_PANELS:
            candidate_class = replicate[
                f"{scope}_candidate_result_class"]
            control_class = replicate.get(
                f"{scope}_{control}_result_class")
            if control_class not in RESULT_CLASSES:
                raise CampaignError(
                    "control result reliability evidence is malformed")
            if candidate_class == "success" or control_class != "success":
                continue
            if construct_class == "weak" and construct_result in (3, 4):
                panel_counts[panel][
                    "candidate_constructor_weak_regressions"] += 1
            else:
                panel_counts[panel][
                    "candidate_nonrepairable_regressions"] += 1

    for panel, counts in panel_counts.items():
        contrast = by_name.get(panel)
        regressions = (
            None if contrast is None
            else contrast.get("recovery_regressions")
        )
        if (
            type(regressions) is not int
            or regressions < 0
            or sum(counts.values()) != regressions
        ):
            raise CampaignError(
                f"candidate reliability classification disagrees with {panel}")
        contrast.update(counts)
    candidate_errors["total"] = sum(candidate_errors.values())
    return enriched, candidate_errors


class JobEvidenceAggregator:
    """Preserve timing/CIs and aggregate authenticated thermal evidence."""

    def __init__(self, phase, jobs):
        self.phase = phase
        self.jobs = tuple(jobs)
        if any(
                not isinstance(job, dict)
                or type(job.get("job_id")) is not str
                for job in self.jobs):
            raise CampaignError("timing job population is malformed")
        self.expected_job_ids = {job["job_id"] for job in self.jobs}
        if len(self.expected_job_ids) != len(self.jobs):
            raise CampaignError("timing job population has duplicate ids")
        self.seen_job_ids = set()
        self.timing_jobs = []
        self.thermal = {
            "jobs": 0,
            "rows": 0,
            "valid_rows": 0,
            "missing_read_rows": 0,
            "invalid_rows": 0,
            "cpu_tctl_max_c": None,
            "dimm_max_c": None,
            "edac_ce_max": 0,
            "edac_ue_max": 0,
            "cpu_assignments": set(),
            "cache_states": set(),
            "thermal_source_identities": set(),
            "cpu_models": set(),
        }

    @staticmethod
    def _panel_counts(replicates, field):
        counts = {name: 0 for name in band.BANDTIMING_PANELS}
        replicate_count = 0
        for replicate in replicates:
            values = (
                replicate.get(field)
                if isinstance(replicate, dict) else None
            )
            if (
                not isinstance(values, list)
                or any(
                    type(panel) is not str or panel not in counts
                    for panel in values
                )
            ):
                raise CampaignError(
                    "timing replicate panel evidence is malformed")
            if values:
                replicate_count += 1
            for panel in values:
                counts[panel] += 1
        return replicate_count, counts

    def add(self, job, receipt):
        job_id = job["job_id"]
        if job_id not in self.expected_job_ids or job_id in self.seen_job_ids:
            raise CampaignError("timing receipt has an unexpected duplicate job")
        self.seen_job_ids.add(job_id)
        arms = receipt.get("arms")
        contrasts = receipt.get("contrasts")
        if (
            not isinstance(arms, list)
            or not all(isinstance(item, dict) for item in arms)
            or not isinstance(contrasts, list)
            or not all(isinstance(item, dict) for item in contrasts)
            or any(
                type(item.get("scope")) is not str
                or type(item.get("arm")) is not str
                for item in arms
            )
            or any(
                type(item.get("name")) is not str
                for item in contrasts
            )
        ):
            raise CampaignError("timing receipt summaries are malformed")
        arm_keys = (
            [(item.get("scope"), item.get("arm")) for item in arms]
        )
        contrast_names = (
            [item.get("name") for item in contrasts]
        )
        if (
            len(arm_keys) != len(TIMING_ARM_KEYS)
            or set(arm_keys) != set(TIMING_ARM_KEYS)
            or len(set(arm_keys)) != len(arm_keys)
            or len(contrast_names) != len(band.BANDTIMING_CROSS_PANELS)
            or set(contrast_names) != set(band.BANDTIMING_CROSS_PANELS)
            or len(set(contrast_names)) != len(contrast_names)
            or not peel_codec._is_sha256(receipt.get("stream_sha256"))
        ):
            raise CampaignError("timing receipt cardinality is invalid")
        replicates = receipt.get("replicates")
        if not isinstance(replicates, list):
            raise CampaignError("timing receipt replicates are missing")
        enriched_contrasts, candidate_errors = \
            _classify_candidate_reliability(contrasts, replicates)
        censored_replicates, censored_panels = self._panel_counts(
            replicates, "censored_panels")
        contaminated_replicates, contaminated_panels = self._panel_counts(
            replicates, "fault_contaminated_panels")

        context = receipt.get("context")
        if (
            not isinstance(context, dict)
            or not isinstance(context.get("bound"), dict)
            or not isinstance(context.get("thermal"), dict)
        ):
            raise CampaignError("timing receipt context is incomplete")
        bound = context["bound"]
        thermal = context["thermal"]
        required_thermal = (
            "rows", "valid_rows", "missing_read_rows", "invalid_rows",
            "cpu_tctl_max_c", "dimm_max_c", "edac_ce_max", "edac_ue_max",
        )
        evidence = receipt.get("evidence")
        if (
            any(name not in thermal for name in required_thermal)
            or any(
                type(thermal.get(name)) is not int
                or thermal[name] < 0
                for name in (
                    "rows", "valid_rows", "missing_read_rows",
                    "invalid_rows", "edac_ce_max", "edac_ue_max",
                )
            )
            or not _finite_number(thermal.get("cpu_tctl_max_c"))
            or not _finite_number(thermal.get("dimm_max_c"))
            or not isinstance(evidence, dict)
            or not peel_codec._is_sha256(
                evidence.get("context_sha256"))
            or not peel_codec._is_sha256(
                evidence.get("final_context_sha256"))
        ):
            raise CampaignError("thermal summary is incomplete")
        if (
            thermal["missing_read_rows"] != 0
            or thermal["invalid_rows"] != 0
            or thermal["edac_ce_max"] != 0
            or thermal["edac_ue_max"] != 0
            or thermal["rows"] != thermal["valid_rows"]
            or thermal["rows"] < 1
        ):
            raise CampaignError(
                "thermal evidence contains a sensor or EDAC error")
        affinity = bound.get("cpu_affinity")
        if (
            not isinstance(affinity, list)
            or len(affinity) != 1
            or type(affinity[0]) is not int
            or affinity[0] < 0
            or bound.get("cache_state") != CACHE_STATE
            or type(bound.get("thermal_device")) is not int
            or bound["thermal_device"] < 0
            or type(bound.get("thermal_inode")) is not int
            or bound["thermal_inode"] < 0
            or type(bound.get("cpu_model")) is not str
            or not bound["cpu_model"]
        ):
            raise CampaignError("timing context policy is invalid")

        self.thermal["jobs"] += 1
        for name in (
            "rows", "valid_rows", "missing_read_rows", "invalid_rows"):
            self.thermal[name] += thermal[name]
        self.thermal["cpu_tctl_max_c"] = (
            thermal["cpu_tctl_max_c"]
            if self.thermal["cpu_tctl_max_c"] is None
            else max(
                self.thermal["cpu_tctl_max_c"],
                thermal["cpu_tctl_max_c"],
            )
        )
        self.thermal["dimm_max_c"] = (
            thermal["dimm_max_c"]
            if self.thermal["dimm_max_c"] is None
            else max(self.thermal["dimm_max_c"], thermal["dimm_max_c"])
        )
        self.thermal["edac_ce_max"] = max(
            self.thermal["edac_ce_max"], thermal["edac_ce_max"])
        self.thermal["edac_ue_max"] = max(
            self.thermal["edac_ue_max"], thermal["edac_ue_max"])
        self.thermal["cpu_assignments"].add(affinity[0])
        self.thermal["cache_states"].add(bound["cache_state"])
        self.thermal["thermal_source_identities"].add((
            bound["thermal_device"], bound["thermal_inode"]))
        self.thermal["cpu_models"].add(bound.get("cpu_model"))

        self.timing_jobs.append({
            "job_id": job_id,
            "order": job["order"],
            "candidate": job["candidate"],
            "K": job["K"],
            "schedule": job["schedule"],
            "stream_sha256": receipt["stream_sha256"],
            "arms": [dict(item) for item in arms],
            "contrasts": enriched_contrasts,
            "candidate_error_outcomes": candidate_errors,
            "censoring": {
                "replicates": len(replicates),
                "censored_replicates": censored_replicates,
                "censored_panel_observations": censored_panels,
                "fault_contaminated_replicates": contaminated_replicates,
                "fault_contaminated_panel_observations":
                    contaminated_panels,
            },
            "context": {
                "context_sha256":
                    receipt["evidence"]["context_sha256"],
                "final_context_sha256":
                    receipt["evidence"]["final_context_sha256"],
                "cpu_affinity": list(affinity),
                "cache_state": bound["cache_state"],
                "thermal_device": bound["thermal_device"],
                "thermal_inode": bound["thermal_inode"],
                "thermal_rows": thermal["rows"],
                "thermal_valid_rows": thermal["valid_rows"],
                "thermal_missing_read_rows": thermal["missing_read_rows"],
                "thermal_invalid_rows": thermal["invalid_rows"],
                "cpu_tctl_max_c": thermal["cpu_tctl_max_c"],
                "dimm_max_c": thermal["dimm_max_c"],
                "edac_ce_max": thermal["edac_ce_max"],
                "edac_ue_max": thermal["edac_ue_max"],
            },
        })

    def finish(self):
        if self.seen_job_ids != self.expected_job_ids:
            raise CampaignError("timing summary has missing jobs")
        if (
            self.thermal["jobs"] != len(self.jobs)
            or self.thermal["missing_read_rows"] != 0
            or self.thermal["invalid_rows"] != 0
            or self.thermal["edac_ce_max"] != 0
            or self.thermal["edac_ue_max"] != 0
            or len(self.thermal["thermal_source_identities"]) != 1
            or self.thermal["cache_states"] != {CACHE_STATE}
            or len(self.thermal["cpu_models"]) != 1
        ):
            raise CampaignError("aggregate thermal/context evidence is invalid")
        timing_jobs = sorted(
            self.timing_jobs, key=lambda item: item["order"])
        if [item["order"] for item in timing_jobs] != list(
                range(len(self.jobs))):
            raise CampaignError("timing jobs are not a complete ordering")
        thermal = dict(self.thermal)
        thermal["cpu_assignments"] = sorted(thermal["cpu_assignments"])
        thermal["cache_states"] = sorted(thermal["cache_states"])
        thermal["thermal_source_identities"] = [
            {"device": device, "inode": inode}
            for device, inode in sorted(
                thermal["thermal_source_identities"])
        ]
        thermal["cpu_models"] = sorted(thermal["cpu_models"])
        return {
            "timing_analysis_unit":
                "one-authenticated-bandtiming-job-per-"
                "candidate-K-schedule",
            "timing_job_count": len(timing_jobs),
            "timing_jobs": timing_jobs,
            "thermal_context": thermal,
        }


def _require_selected_survivor(decision, requested_survivor):
    """Reject caller discretion after the frozen selection evidence exists."""
    if (
        not isinstance(decision, dict)
        or decision.get("schema") != DECISION_SCHEMA
        or decision.get("phase") != "selection"
        or decision.get("policy_sha256")
            != canonical_sha256(selection_policy())
        or "selected_survivor" not in decision
    ):
        raise CampaignError("selection decision evidence is malformed")
    selected = decision["selected_survivor"]
    if selected is None:
        raise CampaignError(
            "selection policy produced no survivor; holdout is forbidden")
    if not _is_frozen_candidate_name(selected):
        raise CampaignError("selection decision named an unknown survivor")
    if requested_survivor != selected:
        raise CampaignError(
            "requested holdout survivor differs from frozen selection")
    rows = decision.get("candidate_metrics")
    selected_rows = (
        [
            row for row in rows
            if isinstance(row, dict)
            and row.get("candidate") == selected
        ]
        if isinstance(rows, list) else []
    )
    blockers = decision.get("selected_candidate_production_blockers")
    if (
        len(selected_rows) != 1
        or _validate_production_blockers(blockers, selected)
            != selected_rows[0].get("production_reliability")
    ):
        raise CampaignError(
            "selected candidate production blockers are malformed")
    return selected


def _arm_recovery_metric(groups, arm):
    final = 0
    tail_auc = 0
    for group in groups:
        summary = group.get(arm)
        if not isinstance(summary, dict):
            raise CampaignError("decision recovery arm is missing")
        arm_final = summary.get("final_unrecovered")
        tail = summary.get("unrecovered_by_overhead_h_0_to_64")
        if (
            type(arm_final) is not int
            or arm_final < 0
            or not isinstance(tail, list)
            or len(tail) != MAX_OVERHEAD + 1
            or any(type(value) is not int or value < 0 for value in tail)
            or tail[-1] != arm_final
        ):
            raise CampaignError("decision recovery tail is malformed")
        final += arm_final
        tail_auc += sum(tail)
    return {
        "final_unrecovered_at_h64": final,
        "overhead_tail_auc_h0_to_h64": tail_auc,
    }


def _finite_number(value):
    return (
        type(value) in (int, float)
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def _student_t_mean_interval(values):
    if (
        not isinstance(values, list)
        or len(values) < 4
        or any(type(value) is not float or not math.isfinite(value)
               for value in values)
    ):
        raise CampaignError(
            "aggregate timing needs at least four finite job means")
    count = len(values)
    mean = math.fsum(values) / count
    variance = math.fsum(
        (value - mean) ** 2 for value in values) / (count - 1)
    half = (
        peel_codec._student_t_critical_95(count - 1)
        * math.sqrt(variance / count)
    )
    low = mean - half
    high = mean + half
    if not all(math.isfinite(value) for value in (mean, low, high)):
        raise CampaignError("aggregate timing interval is not finite")
    return mean, low, high


def _validate_aa_timing_arm(arm, raw_replicates, label):
    eligible = arm.get("aa_eligible_replicates")
    mean = arm.get("aa_log_cost_mean")
    low = arm.get("aa_log_cost_ci_low")
    high = arm.get("aa_log_cost_ci_high")
    floor = arm.get("aa_floor_log")
    if (
        type(eligible) is not int
        or not 0 <= eligible <= raw_replicates
        or not _finite_number(mean)
        or not _finite_number(floor)
        or floor < 0.0
    ):
        raise CampaignError(f"{label} AA timing evidence is malformed")
    available = eligible >= 4
    if available:
        if (
            not _finite_number(low)
            or not _finite_number(high)
            or not low <= mean <= high
            or floor != max(abs(low), abs(high))
        ):
            raise CampaignError(f"{label} AA confidence interval is invalid")
    elif low is not None or high is not None or floor != 0.0:
        raise CampaignError(
            f"{label} unavailable AA confidence interval is contradictory")
    return {
        "eligible_replicates": eligible,
        "ci_available": available,
        "log_cost_mean": float(mean),
        "log_cost_ci_low": None if low is None else float(low),
        "log_cost_ci_high": None if high is None else float(high),
        "floor_log": float(floor),
    }


def _validate_decision_contrast(
        contrast, left_aa, right_aa, raw_replicates, label):
    required_lower, unused_upper = \
        peel_codec._paired_practical_log_bounds(REQUIRED_MARGIN)
    required_floor = -required_lower
    eligible = contrast.get("eligible_replicates")
    mean = contrast.get("log_cost_mean")
    low = contrast.get("log_cost_ci_low")
    high = contrast.get("log_cost_ci_high")
    floor = contrast.get("effective_floor_log")
    available = contrast.get("timing_ci_available")
    left_faster = contrast.get("left_faster")
    counts = {
        name: contrast.get(name)
        for name in (
            "recovery_regressions",
            "recovery_improvements",
            "both_failures",
            "candidate_constructor_weak_regressions",
            "candidate_nonrepairable_regressions",
        )
    }
    if (
        type(eligible) is not int
        or not 0 <= eligible <= raw_replicates
        or not _finite_number(mean)
        or not _finite_number(floor)
        or floor < 0.0
        or type(available) is not bool
        or type(left_faster) is not bool
        or any(type(value) is not int or value < 0
               for value in counts.values())
        or (
            counts["recovery_regressions"]
            + counts["recovery_improvements"]
            + counts["both_failures"]
        ) > raw_replicates
        or (
            eligible
            + counts["recovery_regressions"]
            + counts["recovery_improvements"]
            + counts["both_failures"]
        ) > raw_replicates
        or (
            counts["candidate_constructor_weak_regressions"]
            + counts["candidate_nonrepairable_regressions"]
        ) != counts["recovery_regressions"]
        or contrast.get("required_margin") != REQUIRED_MARGIN
        or contrast.get("left_aa_floor_log") != left_aa["floor_log"]
        or contrast.get("right_aa_floor_log") != right_aa["floor_log"]
        or floor != max(
            required_floor,
            left_aa["floor_log"],
            right_aa["floor_log"],
        )
    ):
        raise CampaignError(f"{label} timing evidence is malformed")
    if available != (eligible >= 4):
        raise CampaignError(f"{label} timing CI availability is contradictory")
    if available:
        if (
            not _finite_number(low)
            or not _finite_number(high)
            or not low <= mean <= high
        ):
            raise CampaignError(f"{label} confidence interval is invalid")
    elif low is not None or high is not None:
        raise CampaignError(
            f"{label} unavailable confidence interval is contradictory")
    if left_faster and (
        not available
        or not left_aa["ci_available"]
        or not right_aa["ci_available"]
        or counts["recovery_regressions"] != 0
        or high >= -floor
    ):
        raise CampaignError(f"{label} left-faster claim is contradictory")
    return {
        "eligible_replicates": eligible,
        "ci_available": available,
        "log_cost_mean": float(mean),
        "log_cost_ci_low": None if low is None else float(low),
        "log_cost_ci_high": None if high is None else float(high),
        "effective_floor_log": float(floor),
        "left_faster": left_faster,
        **counts,
    }


def _timing_metrics(phase, timing_jobs):
    raw_replicates = _phase_parameters(phase)["replicates"]
    coverage_floor = raw_replicates // 2
    records = {
        metric: []
        for metric, unused_panel, unused_scope, unused_left, unused_right
        in DECISION_TIMING_PANELS
    }
    for job in timing_jobs:
        contrasts = job.get("contrasts")
        arms = job.get("arms")
        if (
            not isinstance(contrasts, list)
            or not all(isinstance(item, dict) for item in contrasts)
            or not isinstance(arms, list)
            or not all(isinstance(item, dict) for item in arms)
            or any(
                type(item.get("name")) is not str
                for item in contrasts
            )
            or any(
                type(item.get("scope")) is not str
                or type(item.get("arm")) is not str
                for item in arms
            )
        ):
            raise CampaignError("decision timing evidence is missing")
        by_name = {item.get("name"): item for item in contrasts}
        by_arm = {
            (item.get("scope"), item.get("arm")): item
            for item in arms
        }
        if (
            len(by_name) != len(contrasts)
            or len(by_arm) != len(arms)
            or set(by_arm) != set(TIMING_ARM_KEYS)
        ):
            raise CampaignError("decision timing evidence is duplicated")
        for metric, panel, scope, left, right in DECISION_TIMING_PANELS:
            contrast = by_name.get(panel)
            left_arm = by_arm.get((scope, left))
            right_arm = by_arm.get((scope, right))
            if (
                not isinstance(contrast, dict)
                or not isinstance(left_arm, dict)
                or not isinstance(right_arm, dict)
            ):
                raise CampaignError(
                    f"decision timing contrast {panel} is missing")
            left_aa = _validate_aa_timing_arm(
                left_arm, raw_replicates, f"{panel} left")
            right_aa = _validate_aa_timing_arm(
                right_arm, raw_replicates, f"{panel} right")
            records[metric].append({
                "cross": _validate_decision_contrast(
                    contrast, left_aa, right_aa,
                    raw_replicates, panel),
                "left_aa": left_aa,
                "right_aa": right_aa,
            })

    if not timing_jobs or any(
            len(items) != len(timing_jobs)
            for items in records.values()):
        raise CampaignError("decision timing population is incomplete")

    point_metrics = {}
    confirmation = {}
    timing_complete = True
    for metric, items in records.items():
        means = [item["cross"]["log_cost_mean"] for item in items]
        mean, t_low, t_high = _student_t_mean_interval(means)
        cross_complete = all(
            item["cross"]["ci_available"]
            and item["cross"]["eligible_replicates"] >= coverage_floor
            for item in items
        )
        aa_complete = all(
            item[side]["ci_available"]
            and item[side]["eligible_replicates"] >= coverage_floor
            for item in items
            for side in ("left_aa", "right_aa")
        )
        evidence_complete = cross_complete and aa_complete
        timing_complete = timing_complete and evidence_complete
        mean_ci_high = (
            math.fsum(
                item["cross"]["log_cost_ci_high"] for item in items)
            / len(items)
            if all(item["cross"]["ci_available"] for item in items)
            else None
        )
        mean_floor = (
            math.fsum(
                item["cross"]["effective_floor_log"] for item in items)
            / len(items)
        )
        aggregate_upper = (
            max(t_high, mean_ci_high)
            if mean_ci_high is not None else None
        )
        aggregate_pass = (
            evidence_complete
            and aggregate_upper < -mean_floor
        )
        point_metrics[metric] = mean
        confirmation[metric] = {
            "job_count": len(items),
            "coverage_floor_per_job": coverage_floor,
            "cross_minimum_eligible_replicates": min(
                item["cross"]["eligible_replicates"] for item in items),
            "cross_ci_available_job_count": sum(
                item["cross"]["ci_available"] for item in items),
            "cross_jobs_meeting_coverage_floor": sum(
                item["cross"]["eligible_replicates"] >= coverage_floor
                for item in items),
            "aa_minimum_eligible_replicates": min(
                item[side]["eligible_replicates"]
                for item in items
                for side in ("left_aa", "right_aa")),
            "left_aa_minimum_eligible_replicates": min(
                item["left_aa"]["eligible_replicates"]
                for item in items),
            "right_aa_minimum_eligible_replicates": min(
                item["right_aa"]["eligible_replicates"]
                for item in items),
            "aa_ci_available_arm_count": sum(
                item[side]["ci_available"]
                for item in items
                for side in ("left_aa", "right_aa")),
            "aa_arms_meeting_coverage_floor": sum(
                item[side]["eligible_replicates"] >= coverage_floor
                for item in items
                for side in ("left_aa", "right_aa")),
            "left_aa_jobs_meeting_coverage_floor": sum(
                item["left_aa"]["eligible_replicates"] >= coverage_floor
                for item in items),
            "right_aa_jobs_meeting_coverage_floor": sum(
                item["right_aa"]["eligible_replicates"] >= coverage_floor
                for item in items),
            "aa_arm_observations": 2 * len(items),
            "log_cost_mean": mean,
            "job_mean_t95_ci_low": t_low,
            "job_mean_t95_ci_high": t_high,
            "mean_per_job_ci_high": mean_ci_high,
            "mean_effective_floor_log": mean_floor,
            "aggregate_upper_log_cost": aggregate_upper,
            "left_faster_job_count": sum(
                item["cross"]["left_faster"] for item in items),
            "recovery_regressions": sum(
                item["cross"]["recovery_regressions"] for item in items),
            "recovery_improvements": sum(
                item["cross"]["recovery_improvements"] for item in items),
            "both_failures": sum(
                item["cross"]["both_failures"] for item in items),
            "candidate_constructor_weak_regressions": sum(
                item["cross"][
                    "candidate_constructor_weak_regressions"]
                for item in items),
            "candidate_nonrepairable_regressions": sum(
                item["cross"]["candidate_nonrepairable_regressions"]
                for item in items),
            "evidence_complete": evidence_complete,
            "aggregate_confirmation_pass": aggregate_pass,
        }
    return point_metrics, confirmation, timing_complete


def _validate_production_blockers(value, candidate):
    panel_names = {
        panel for panel, unused_scope, unused_control
        in CANDIDATE_CONTROL_PANELS
    }
    if (
        not isinstance(value, dict)
        or set(value) != {
            "schema", "candidate", "candidate_weak_constructions",
            "candidate_error_outcomes", "candidate_control_panels",
            "repairable_seed_blocker",
            "nonrepairable_reliability_blocker",
        }
        or value.get("schema") != PRODUCTION_BLOCKER_SCHEMA
        or value.get("candidate") != candidate
        or type(value.get("candidate_weak_constructions")) is not int
        or value["candidate_weak_constructions"] < 0
        or type(value.get("repairable_seed_blocker")) is not bool
        or type(value.get("nonrepairable_reliability_blocker")) is not bool
    ):
        raise CampaignError("production blocker evidence is malformed")
    errors = value["candidate_error_outcomes"]
    if (
        not isinstance(errors, dict)
        or set(errors) != {
            "constructor", "encoder", "decoder", "direct", "total"
        }
        or any(type(item) is not int or item < 0
               for item in errors.values())
        or errors["total"] != sum(
            errors[name]
            for name in ("constructor", "encoder", "decoder", "direct"))
    ):
        raise CampaignError("candidate error evidence is malformed")
    panels = value["candidate_control_panels"]
    if (
        not isinstance(panels, dict)
        or set(panels) != panel_names
    ):
        raise CampaignError("candidate panel blocker evidence is malformed")
    nonrepairable = 0
    weak_regressions = 0
    for panel in panel_names:
        counts = panels[panel]
        if (
            not isinstance(counts, dict)
            or set(counts) != {
                "recovery_regressions", "recovery_improvements",
                "both_failures",
                "candidate_constructor_weak_regressions",
                "candidate_nonrepairable_regressions",
            }
            or any(type(item) is not int or item < 0
                   for item in counts.values())
            or (
                counts["candidate_constructor_weak_regressions"]
                + counts["candidate_nonrepairable_regressions"]
            ) != counts["recovery_regressions"]
        ):
            raise CampaignError(
                "candidate panel blocker counts are malformed")
        weak_regressions += \
            counts["candidate_constructor_weak_regressions"]
        nonrepairable += counts["candidate_nonrepairable_regressions"]
    if (
        (weak_regressions > 0
         and value["candidate_weak_constructions"] == 0)
        or value["repairable_seed_blocker"]
            != (value["candidate_weak_constructions"] > 0)
        or value["nonrepairable_reliability_blocker"]
            != (errors["total"] > 0 or nonrepairable > 0)
    ):
        raise CampaignError("production blocker flags are inconsistent")
    return value


def _production_reliability_metrics(
        phase, candidate, timing_jobs, weak_count):
    raw_replicates = _phase_parameters(phase)["replicates"]
    panels = {
        panel: {
            "recovery_regressions": 0,
            "recovery_improvements": 0,
            "both_failures": 0,
            "candidate_constructor_weak_regressions": 0,
            "candidate_nonrepairable_regressions": 0,
        }
        for panel, unused_scope, unused_control
        in CANDIDATE_CONTROL_PANELS
    }
    errors = {
        "constructor": 0,
        "encoder": 0,
        "decoder": 0,
        "direct": 0,
        "total": 0,
    }
    for job in timing_jobs:
        job_errors = job.get("candidate_error_outcomes")
        if (
            not isinstance(job_errors, dict)
            or set(job_errors) != set(errors)
            or any(type(item) is not int or item < 0
                   for item in job_errors.values())
            or any(
                job_errors[field] > raw_replicates
                for field in ("constructor", "encoder", "decoder", "direct"))
            or job_errors["total"] != sum(
                job_errors[name]
                for name in ("constructor", "encoder", "decoder", "direct"))
        ):
            raise CampaignError(
                "decision candidate error evidence is malformed")
        for field in errors:
            errors[field] += job_errors[field]
        contrasts = job.get("contrasts")
        if (
            not isinstance(contrasts, list)
            or not all(isinstance(item, dict) for item in contrasts)
            or any(
                type(item.get("name")) is not str
                for item in contrasts
            )
        ):
            raise CampaignError(
                "decision production contrasts are missing")
        by_name = {item.get("name"): item for item in contrasts}
        if len(by_name) != len(contrasts):
            raise CampaignError(
                "decision production contrasts are duplicated")
        for panel in panels:
            contrast = by_name.get(panel)
            if not isinstance(contrast, dict):
                raise CampaignError(
                    f"decision production contrast {panel} is missing")
            for field in panels[panel]:
                value = contrast.get(field)
                if type(value) is not int or value < 0:
                    raise CampaignError(
                        "decision production count is malformed")
                panels[panel][field] += value
            if (
                contrast["candidate_constructor_weak_regressions"]
                + contrast["candidate_nonrepairable_regressions"]
                != contrast["recovery_regressions"]
                or (
                    contrast["recovery_regressions"]
                    + contrast["recovery_improvements"]
                    + contrast["both_failures"]
                ) > raw_replicates
            ):
                raise CampaignError(
                    "decision production regression split is inconsistent")
    value = {
        "schema": PRODUCTION_BLOCKER_SCHEMA,
        "candidate": candidate,
        "candidate_weak_constructions": weak_count,
        "candidate_error_outcomes": errors,
        "candidate_control_panels": panels,
        "repairable_seed_blocker": weak_count > 0,
        "nonrepairable_reliability_blocker": (
            errors["total"] > 0
            or any(
                counts["candidate_nonrepairable_regressions"] > 0
                for counts in panels.values())
        ),
    }
    return _validate_production_blockers(value, candidate)


def _combine_production_blockers(
        selection_blockers, holdout_blockers, candidate):
    selection_blockers = _validate_production_blockers(
        selection_blockers, candidate)
    holdout_blockers = _validate_production_blockers(
        holdout_blockers, candidate)
    errors = {
        field: (
            selection_blockers["candidate_error_outcomes"][field]
            + holdout_blockers["candidate_error_outcomes"][field]
        )
        for field in (
            "constructor", "encoder", "decoder", "direct", "total")
    }
    panels = {}
    for panel, unused_scope, unused_control in CANDIDATE_CONTROL_PANELS:
        panels[panel] = {
            field: (
                selection_blockers["candidate_control_panels"][
                    panel][field]
                + holdout_blockers["candidate_control_panels"][
                    panel][field]
            )
            for field in (
                "recovery_regressions", "recovery_improvements",
                "both_failures",
                "candidate_constructor_weak_regressions",
                "candidate_nonrepairable_regressions",
            )
        }
    weak_count = (
        selection_blockers["candidate_weak_constructions"]
        + holdout_blockers["candidate_weak_constructions"]
    )
    combined = {
        "schema": PRODUCTION_BLOCKER_SCHEMA,
        "candidate": candidate,
        "candidate_weak_constructions": weak_count,
        "candidate_error_outcomes": errors,
        "candidate_control_panels": panels,
        "repairable_seed_blocker": weak_count > 0,
        "nonrepairable_reliability_blocker": (
            errors["total"] > 0
            or any(
                counts["candidate_nonrepairable_regressions"] > 0
                for counts in panels.values())
        ),
    }
    return _validate_production_blockers(combined, candidate)


def _candidate_decision_metrics(phase, aggregates, timing_and_thermal):
    if (
        not isinstance(aggregates, dict)
        or not isinstance(timing_and_thermal, dict)
    ):
        raise CampaignError("decision aggregate evidence is malformed")
    groups_evidence = aggregates.get("groups")
    timing_jobs_evidence = timing_and_thermal.get("timing_jobs")
    if (
        not isinstance(groups_evidence, list)
        or not all(isinstance(group, dict) for group in groups_evidence)
        or not isinstance(timing_jobs_evidence, list)
        or not all(isinstance(job, dict) for job in timing_jobs_evidence)
    ):
        raise CampaignError("decision aggregate population is malformed")
    if phase == "selection":
        expected_candidates = set(CANDIDATE_BY_NAME)
    elif phase == "holdout":
        if any(
                not _is_frozen_candidate_name(group.get("candidate_name"))
                for group in groups_evidence):
            raise CampaignError("decision candidate population is invalid")
        expected_candidates = {
            group["candidate_name"] for group in groups_evidence
        }
    else:
        raise CampaignError("decision candidate population is invalid")
    if (
        not expected_candidates
        or (phase == "holdout" and len(expected_candidates) != 1)
    ):
        raise CampaignError("decision candidate population is invalid")
    expected_cells = {
        (block_count, schedule)
        for block_count in K_VALUES
        for schedule in SCHEDULES
    }
    groups_by_candidate = {name: [] for name in expected_candidates}
    for group in groups_evidence:
        candidate_name = group.get("candidate_name")
        block_count = group.get("K")
        schedule = group.get("schedule")
        if (
            not isinstance(group, dict)
            or not _is_frozen_candidate_name(candidate_name)
            or candidate_name not in groups_by_candidate
            or type(block_count) is not int
            or block_count not in K_VALUES
            or type(schedule) is not str
            or schedule not in SCHEDULES
        ):
            raise CampaignError("decision recovery group is unexpected")
        groups_by_candidate[candidate_name].append(group)
    timing_by_candidate = {name: [] for name in expected_candidates}
    for job in timing_jobs_evidence:
        candidate_name = job.get("candidate")
        block_count = job.get("K")
        schedule = job.get("schedule")
        if (
            not isinstance(job, dict)
            or not _is_frozen_candidate_name(candidate_name)
            or candidate_name not in timing_by_candidate
            or type(block_count) is not int
            or block_count not in K_VALUES
            or type(schedule) is not str
            or schedule not in SCHEDULES
        ):
            raise CampaignError("decision timing job is unexpected")
        timing_by_candidate[candidate_name].append(job)
    weakness = aggregates.get("independent_construction_weaknesses")
    weak = (
        weakness.get("candidate")
        if isinstance(weakness, dict) else None
    )
    if not isinstance(weak, dict) or set(weak) != expected_candidates:
        raise CampaignError("decision weak-construction evidence is invalid")

    result = {}
    control_reference = None
    for name in sorted(expected_candidates):
        groups = groups_by_candidate[name]
        timing_jobs = timing_by_candidate[name]
        if (
            {
                (group.get("K"), group.get("schedule"))
                for group in groups
            } != expected_cells
            or len(groups) != len(expected_cells)
            or {
                (job.get("K"), job.get("schedule"))
                for job in timing_jobs
            } != expected_cells
            or len(timing_jobs) != len(expected_cells)
        ):
            raise CampaignError(
                "decision evidence lacks the exact K/schedule population")
        recovery = {
            arm: _arm_recovery_metric(groups, arm)
            for arm in ("candidate", "dispatch", "wh1")
        }
        controls = (recovery["dispatch"], recovery["wh1"])
        if control_reference is None:
            control_reference = controls
        elif controls != control_reference:
            raise CampaignError(
                "decision controls changed across candidate copies")
        timing, timing_confirmation, timing_complete = \
            _timing_metrics(phase, timing_jobs)
        weak_item = weak[name]
        weak_count = (
            weak_item.get("count")
            if isinstance(weak_item, dict) else None
        )
        if type(weak_count) is not int or weak_count < 0:
            raise CampaignError(
                "decision weak-construction count is invalid")
        candidate = recovery["candidate"]
        recovery_pass = all(
            candidate["final_unrecovered_at_h64"]
                <= recovery[control]["final_unrecovered_at_h64"]
            and candidate["overhead_tail_auc_h0_to_h64"]
                <= recovery[control]["overhead_tail_auc_h0_to_h64"]
            for control in ("dispatch", "wh1")
        )
        result[name] = {
            "candidate": name,
            "recovery": recovery,
            "unique_weak_constructions": weak_count,
            "timing": timing,
            "timing_confirmation": timing_confirmation,
            "timing_complete": timing_complete,
            "recovery_constraints_pass": recovery_pass,
            "production_reliability": _production_reliability_metrics(
                phase, name, timing_jobs, weak_count),
        }
    return result


def derive_campaign_decision(
        phase, aggregates, timing_and_thermal, *, survivor=None,
        selection_production_blockers=None):
    """Execute the prebound selection/holdout policy on replayed evidence."""
    metrics = _candidate_decision_metrics(
        phase, aggregates, timing_and_thermal)
    policy_sha256 = canonical_sha256(selection_policy())
    dimensions = (
        ("recovery", "candidate", "final_unrecovered_at_h64"),
        ("recovery", "candidate", "overhead_tail_auc_h0_to_h64"),
        ("unique_weak_constructions",),
        ("timing", "direct_candidate_dispatch_log_cost"),
        ("timing", "encoder_candidate_wh1_log_cost"),
        ("timing", "decoder_candidate_wh1_log_cost"),
    )

    def value(row, path):
        for item in path:
            row = row[item]
        return row

    if phase == "selection":
        if (
            survivor is not None
            or selection_production_blockers is not None
            or set(metrics) != set(CANDIDATE_BY_NAME)
        ):
            raise CampaignError("selection decision population is invalid")
        complete = {
            name for name, row in metrics.items()
            if row["timing_complete"]
        }
        pareto = set()
        for name in complete:
            row = metrics[name]
            dominated = any(
                all(
                    value(other, path) <= value(row, path)
                    for path in dimensions
                )
                and any(
                    value(other, path) < value(row, path)
                    for path in dimensions
                )
                for other_name, other in metrics.items()
                if other_name in complete and other_name != name
            )
            if not dominated:
                pareto.add(name)
        eligible = sorted(
            name for name in pareto
            if metrics[name]["recovery_constraints_pass"]
        )
        selected = (
            min(
                eligible,
                key=lambda name: (
                    metrics[name]["timing"][
                        "direct_candidate_dispatch_log_cost"],
                    name,
                ),
            )
            if eligible else None
        )
        rows = []
        for name in sorted(metrics):
            rows.append({
                **metrics[name],
                "on_raw_pareto_surface": name in pareto,
                "selection_eligible":
                    name in pareto
                    and metrics[name]["timing_complete"]
                    and metrics[name]["recovery_constraints_pass"],
            })
        return {
            "schema": DECISION_SCHEMA,
            "phase": "selection",
            "policy_sha256": policy_sha256,
            "candidate_metrics": rows,
            "raw_pareto_candidates": sorted(pareto),
            "eligible_pareto_candidates": eligible,
            "selected_survivor": selected,
            "selected_candidate_production_blockers": (
                metrics[selected]["production_reliability"]
                if selected is not None else None
            ),
            "status": "selected" if selected is not None else "no-survivor",
        }

    if phase != "holdout" or not _is_frozen_candidate_name(survivor):
        raise CampaignError("holdout decision population is invalid")
    if set(metrics) != {survivor}:
        raise CampaignError("holdout decision tested the wrong survivor")
    row = metrics[survivor]
    selection_production_blockers = _validate_production_blockers(
        selection_production_blockers, survivor)
    cumulative_blockers = _combine_production_blockers(
        selection_production_blockers,
        row["production_reliability"],
        survivor,
    )
    timing = row["timing"]
    confirmation = row["timing_confirmation"]
    direct_metric = "direct_candidate_dispatch_log_cost"
    encoder_metric = "encoder_candidate_wh1_log_cost"
    decoder_metric = "decoder_candidate_wh1_log_cost"
    direct_pass = confirmation[direct_metric][
        "aggregate_confirmation_pass"]
    holdout_confirmed = row["recovery_constraints_pass"] and direct_pass
    metric_panels = {
        direct_metric: "direct_candidate_dispatch",
        encoder_metric: "encoder_candidate_wh1",
        decoder_metric: "decoder_candidate_wh1",
    }
    production_timing_evidence = {}
    for metric, panel in metric_panels.items():
        panel_clean = (
            cumulative_blockers["candidate_control_panels"][panel][
                "candidate_nonrepairable_regressions"] == 0
        )
        production_timing_evidence[metric] = {
            **confirmation[metric],
            "cumulative_candidate_nonrepairable_regressions":
                cumulative_blockers["candidate_control_panels"][panel][
                    "candidate_nonrepairable_regressions"],
            "production_objective_pass":
                confirmation[metric]["aggregate_confirmation_pass"]
                and panel_clean,
        }
    direct_production_pass = production_timing_evidence[direct_metric][
        "production_objective_pass"]
    encoder_timing_pass = confirmation[encoder_metric][
        "aggregate_confirmation_pass"]
    decoder_timing_pass = confirmation[decoder_metric][
        "aggregate_confirmation_pass"]
    encoder_pass = production_timing_evidence[encoder_metric][
        "production_objective_pass"]
    decoder_pass = production_timing_evidence[decoder_metric][
        "production_objective_pass"]
    reliability_clean = (
        not cumulative_blockers["repairable_seed_blocker"]
        and not cumulative_blockers[
            "nonrepairable_reliability_blocker"]
    )
    promotion_ready = (
        holdout_confirmed
        and direct_production_pass
        and encoder_pass
        and decoder_pass
        and reliability_clean
    )
    if promotion_ready:
        production_action = "promote-new-contract"
    elif (
        not holdout_confirmed
        or not encoder_timing_pass
        or not decoder_timing_pass
    ):
        production_action = "retain-dispatch"
    elif (
        cumulative_blockers["nonrepairable_reliability_blocker"]
    ):
        production_action = "investigate-reliability"
    elif cumulative_blockers["repairable_seed_blocker"]:
        production_action = "repair-seeds-and-retest"
    else:
        production_action = "retain-dispatch"
    return {
        "schema": DECISION_SCHEMA,
        "phase": "holdout",
        "policy_sha256": policy_sha256,
        "tested_survivor": survivor,
        "candidate_metrics": row,
        "recovery_confirmed": row["recovery_constraints_pass"],
        "direct_speed_confirmed": direct_pass,
        "direct_production_objective_confirmed":
            direct_production_pass,
        "direct_speed_evidence":
            production_timing_evidence[direct_metric],
        "holdout_confirmed": holdout_confirmed,
        "wh1_throughput_objectives": {
            "encoder_confirmed": encoder_pass,
            "decoder_confirmed": decoder_pass,
            "encoder_candidate_wh1_log_cost":
                timing["encoder_candidate_wh1_log_cost"],
            "decoder_candidate_wh1_log_cost":
                timing["decoder_candidate_wh1_log_cost"],
            "encoder_evidence":
                production_timing_evidence[encoder_metric],
            "decoder_evidence":
                production_timing_evidence[decoder_metric],
        },
        "production_reliability": {
            "selection":
                selection_production_blockers,
            "holdout": row["production_reliability"],
            "cumulative": cumulative_blockers,
        },
        "production_promotion_ready": promotion_ready,
        "production_action": production_action,
    }


class AtomicGzipJsonLines:
    """Stream canonical JSON lines to one atomic deterministic gzip member."""

    def __init__(self, path):
        self.path = Path(path)
        self.temporary = self.path.with_name(self.path.name + ".tmp")
        self.raw_sha256 = hashlib.sha256()
        self.raw_bytes = 0
        self.rows = 0
        self.file = None
        self.gzip_file = None

    def __enter__(self):
        self.file = open(self.temporary, "xb")
        self.gzip_file = gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=6,
            fileobj=self.file,
            mtime=0,
        )
        return self

    def write(self, value):
        payload = canonical_json_bytes(value)
        self.raw_sha256.update(payload)
        self.raw_bytes += len(payload)
        self.rows += 1
        self.gzip_file.write(payload)

    def __exit__(self, exception_type, exception, traceback):
        try:
            if self.gzip_file is not None:
                self.gzip_file.close()
            if self.file is not None:
                self.file.flush()
                os.fsync(self.file.fileno())
                self.file.close()
            if exception_type is None:
                _commit_no_replace(self.temporary, self.path)
            else:
                try:
                    os.unlink(self.temporary)
                except FileNotFoundError:
                    pass
        except BaseException:
            try:
                os.unlink(self.temporary)
            except FileNotFoundError:
                pass
            raise

    def evidence(self):
        compressed = self.path.read_bytes()
        return {
            "rows": self.rows,
            "compressed_sha256": sha256_bytes(compressed),
            "uncompressed_sha256": self.raw_sha256.hexdigest(),
            "compressed_bytes": len(compressed),
            "uncompressed_bytes": self.raw_bytes,
        }


def _replay_result(
        directory, manifest, manifest_sha256, assignment, job):
    path = directory / assignment["output"]
    envelope, evidence = read_canonical_gzip_json(path)
    required = {
        "schema", "manifest_sha256", "job", "assignment",
        "wall_started_unix_ns", "wall_finished_unix_ns",
        "runtime_bindings_before", "runtime_bindings_after", "receipt",
    }
    if (
        not isinstance(envelope, dict)
        or set(envelope) != required
        or envelope["schema"] != RESULT_SCHEMA
        or envelope["manifest_sha256"] != manifest_sha256
        or envelope["job"] != job
        or envelope["assignment"] != assignment
        or envelope["runtime_bindings_before"]
            != manifest["runtime_bindings"]
        or envelope["runtime_bindings_after"]
            != manifest["runtime_bindings"]
        or type(envelope["wall_started_unix_ns"]) is not int
        or type(envelope["wall_finished_unix_ns"]) is not int
        or envelope["wall_finished_unix_ns"]
            <= envelope["wall_started_unix_ns"]
    ):
        raise CampaignError(f"result envelope is inconsistent: {path}")
    receipt = envelope["receipt"]
    request = expected_request(job)
    band.replay_bandtiming_receipt(receipt, expected_request=request)
    if (
        receipt["candidate_descriptor"] != job["candidate_descriptor"]
        or receipt["dispatch_descriptor"] != job["dispatch_descriptor"]
        or len(receipt["replicates"]) != job["replicates"]
    ):
        raise CampaignError(f"native receipt changed its job: {path}")
    bound = receipt["context"]["bound"]
    if (
        bound["cpu_affinity"] != [assignment["cpu"]]
        or bound["thermal_device"]
            != manifest["runtime_bindings"]["thermal"]["device"]
        or bound["thermal_inode"]
            != manifest["runtime_bindings"]["thermal"]["inode"]
    ):
        raise CampaignError(f"native runtime binding is inconsistent: {path}")
    return receipt, evidence


def build_summaries(directory, manifest, output_directory):
    """Strictly replay all results and build deterministic compact summaries."""
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    ledger_path = output_directory / "cell-ledger.jsonl.gz"
    summary_path = output_directory / "summary.json.gz"
    jobs = tuple(manifest["pre_cpu_jobs"])
    manifest_sha256 = canonical_sha256(manifest)
    job_files = _verify_job_files(
        directory, manifest, manifest_sha256)
    aggregator = RecoveryAggregator(manifest["phase"], jobs)
    job_evidence = JobEvidenceAggregator(manifest["phase"], jobs)
    result_files = []
    stream_hashes = set()
    with AtomicGzipJsonLines(ledger_path) as ledger:
        for assignment, job in zip(manifest["assignments"], jobs):
            _parent_signal_safe_point()
            receipt, evidence = _replay_result(
                directory, manifest, manifest_sha256, assignment, job)
            stream = receipt["stream_sha256"]
            if stream in stream_hashes:
                raise CampaignError("duplicate native stream hash")
            stream_hashes.add(stream)
            job_evidence.add(job, receipt)
            result_files.append({
                "path": assignment["output"],
                **evidence,
                "stream_sha256": stream,
            })
            for replicate in receipt["replicates"]:
                ledger.write(aggregator.add_replicate(job, replicate))
    ledger_evidence = ledger.evidence()
    aggregates = aggregator.finish()
    timing_and_thermal = job_evidence.finish()
    decision = derive_campaign_decision(
        manifest["phase"],
        aggregates,
        timing_and_thermal,
        survivor=manifest["survivor"],
        selection_production_blockers=(
            manifest["selection_contract"][
                "selected_candidate_production_blockers"]
            if manifest["phase"] == "holdout" else None
        ),
    )
    if timing_and_thermal["thermal_context"]["cpu_assignments"] != \
            manifest["worker_cpus"]:
        raise CampaignError(
            "thermal evidence did not cover every assigned worker CPU")
    summary = {
        "schema": SUMMARY_SCHEMA,
        "phase": manifest["phase"],
        "survivor": manifest["survivor"],
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "pre_cpu_job_list_sha256":
            manifest["pre_cpu_job_list_sha256"],
        "job_count": len(jobs),
        "cell_ledger": {
            "path": "cell-ledger.jsonl.gz",
            **ledger_evidence,
        },
        "decision": decision,
        **aggregates,
        **timing_and_thermal,
    }
    summary_evidence = atomic_gzip_json(summary_path, summary)
    return {
        "job_files": job_files,
        "result_files": result_files,
        "cell_ledger": {
            "path": "cell-ledger.jsonl.gz",
            **ledger_evidence,
        },
        "summary": {
            "path": "summary.json.gz",
            **summary_evidence,
        },
    }


def _write_completion_checksum(directory, completion_sha256):
    path = directory / "completion.sha256"
    temporary = path.with_name(path.name + ".tmp")
    payload = f"{completion_sha256}  completion.json\n".encode("ascii")
    try:
        with open(temporary, "xb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        _commit_no_replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


@contextmanager
def _parent_termination_handlers():
    """Record TERM/HUP and raise only at explicit cleanup-safe points."""
    global _PARENT_SIGNAL_STATE
    handled = (signal.SIGTERM, signal.SIGHUP)
    previous = {
        item: signal.getsignal(item)
        for item in handled
    }
    if _PARENT_SIGNAL_STATE is not None:
        raise CampaignError("campaign parent signal handler is already active")
    state = {"pending": None}

    def interrupt(signum, unused_frame):
        # Never raise asynchronously.  In particular, a Python handler may
        # run between entry into an ``except`` clause and its first cleanup
        # statement, where raising would skip that same handler entirely.
        if state["pending"] is None:
            state["pending"] = signum

    completed = False
    try:
        _PARENT_SIGNAL_STATE = state
        for item in handled:
            signal.signal(item, interrupt)
    except BaseException:
        for item in handled:
            signal.signal(item, previous[item])
        _PARENT_SIGNAL_STATE = None
        raise
    try:
        yield
        completed = True
    finally:
        # Block delivery while handing the dispositions back.  A signal
        # arriving before the block is recorded; one arriving while blocked
        # is delivered under the restored prior disposition after cleanup.
        old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, handled)
        try:
            pending = state["pending"]
            for item in handled:
                signal.signal(item, previous[item])
            _PARENT_SIGNAL_STATE = None
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        if completed and pending is not None:
            raise _parent_signal_error(pending)


def _run_campaign_impl(
        phase, directory, bench_path, thermal_source, *,
        survivor=None, selection_manifest=None, workers=MAX_WORKERS):
    """Run one complete frozen phase and finalize only after strict replay."""
    _require_normal_priority()
    if type(workers) is not int or not 1 <= workers <= MAX_WORKERS:
        raise CampaignError("worker count must be in [1, 32]")
    if phase == "selection":
        if survivor is not None or selection_manifest is not None:
            raise CampaignError(
                "selection forbids survivor and selection-manifest")
        selection_contract = None
    elif phase == "holdout":
        if (
            not _is_frozen_candidate_name(survivor)
            or selection_manifest is None
        ):
            raise CampaignError(
                "holdout requires an explicit frozen survivor and "
                "completed selection manifest")
        selection_contract = _load_selection_contract(
            selection_manifest, survivor)
    else:
        raise CampaignError(f"unknown campaign phase {phase!r}")

    _parent_signal_safe_point()
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < workers:
        raise CampaignError(
            f"requested {workers} workers but only {len(cpus)} CPUs allowed")
    cpus = cpus[:workers]
    jobs = build_pre_cpu_jobs(phase, survivor)
    _parent_signal_safe_point()
    bindings = runtime_bindings(bench_path, thermal_source)
    directory = _create_fresh_directory(Path(directory).resolve())
    for name in ("jobs", "results", "logs"):
        (directory / name).mkdir()
    manifest = _build_manifest(
        phase, survivor, jobs, cpus, workers, bindings, selection_contract)
    manifest_sha256 = atomic_json(directory / "manifest.json", manifest)
    job_files = _write_job_files(directory, manifest, manifest_sha256)
    _parent_signal_safe_point()

    runtime_monitor = []
    assignments = manifest["assignments"]
    for start in range(0, len(assignments), workers):
        _parent_signal_safe_point()
        check_runtime_bindings(bindings, full_hash=True)
        runtime_monitor.append(_run_wave(
            directory, manifest, assignments[start:start + workers]))
        check_runtime_bindings(bindings, full_hash=True)
        _parent_signal_safe_point()
    log_files = _verify_empty_logs(directory, manifest)
    summary_evidence = build_summaries(directory, manifest, directory)
    _parent_signal_safe_point()
    if summary_evidence["job_files"] != job_files:
        raise CampaignError("job file hashes changed during campaign")
    final_bindings = runtime_bindings(
        bindings["benchmark"]["path"], bindings["thermal"]["path"])
    if final_bindings != bindings:
        raise CampaignError("campaign runtime binding changed before completion")
    completion = {
        "schema": COMPLETION_SCHEMA,
        "phase": phase,
        "survivor": survivor,
        "manifest_sha256": manifest_sha256,
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "pre_cpu_job_list_sha256":
            manifest["pre_cpu_job_list_sha256"],
        "jobs": len(jobs),
        "runtime_binding_monitor": {
            "total_checks": sum(
                wave["checks"] for wave in runtime_monitor),
            "total_checks_while_running": sum(
                wave["checks_while_running"] for wave in runtime_monitor),
            "waves": runtime_monitor,
        },
        "runtime_bindings_final": final_bindings,
        "log_files": log_files,
        **summary_evidence,
        "completed_unix_ns": time.time_ns(),
    }
    completion_sha256 = atomic_json(
        directory / "completion.json", completion)
    _write_completion_checksum(directory, completion_sha256)
    _parent_signal_safe_point()
    return completion


def run_campaign(
        phase, directory, bench_path, thermal_source, *,
        survivor=None, selection_manifest=None, workers=MAX_WORKERS):
    """Run a phase while TERM/HUP deterministically reap worker groups."""
    # Only the campaign parent enters this wrapper.  ``_worker`` dispatches
    # directly to ``run_worker`` and retains default signal dispositions, so
    # the parent's process-group TERM/KILL cleanup remains effective.
    with _parent_termination_handlers():
        completion = _run_campaign_impl(
            phase,
            directory,
            bench_path,
            thermal_source,
            survivor=survivor,
            selection_manifest=selection_manifest,
            workers=workers,
        )
        _parent_signal_safe_point()
        return completion


def _validate_manifest(manifest):
    common_runtime_fields = {
        "created_unix_ns", "repository", "workers", "worker_cpus",
        "priority", "runtime_bindings", "selection_contract",
        "pre_cpu_jobs", "assignments",
    }
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema") != CAMPAIGN_SCHEMA
        or manifest.get("phase") not in ("selection", "holdout")
    ):
        raise CampaignError("campaign manifest is not the frozen contract")
    phase = manifest["phase"]
    survivor = manifest.get("survivor")
    if phase == "holdout" and not isinstance(survivor, str):
        raise CampaignError("holdout survivor is malformed")
    expected_plan = make_plan(phase, survivor)
    worker_cpus = manifest.get("worker_cpus")
    if (
        set(manifest) != set(expected_plan) | common_runtime_fields
        or any(manifest.get(name) != value
               for name, value in expected_plan.items())
        or type(manifest.get("created_unix_ns")) is not int
        or manifest["created_unix_ns"] <= 0
        or manifest.get("repository") != str(REPOSITORY)
        or type(manifest.get("workers")) is not int
        or not 1 <= manifest["workers"] <= MAX_WORKERS
        or not isinstance(worker_cpus, list)
        or len(worker_cpus) != manifest["workers"]
        or any(type(cpu) is not int or cpu < 0
               for cpu in worker_cpus)
        or manifest.get("priority") != "normal-nice-0"
        or not isinstance(manifest.get("runtime_bindings"), dict)
    ):
        raise CampaignError("campaign manifest runtime contract is invalid")
    if worker_cpus != sorted(set(worker_cpus)):
        raise CampaignError("campaign manifest runtime contract is invalid")
    _validate_runtime_bindings_schema(manifest["runtime_bindings"])
    selection_contract = manifest.get("selection_contract")
    if phase == "selection":
        if selection_contract is not None:
            raise CampaignError("selection manifest fabricated a parent")
    elif (
        not isinstance(selection_contract, dict)
        or set(selection_contract) != {
            "path", "manifest_sha256", "completion_sha256",
            "selection_decision_sha256", "selected_survivor",
            "selected_candidate_production_blockers",
        }
        or not isinstance(selection_contract["path"], str)
        or not Path(selection_contract["path"]).is_absolute()
        or not peel_codec._is_sha256(
            selection_contract["manifest_sha256"])
        or not peel_codec._is_sha256(
            selection_contract["completion_sha256"])
        or not peel_codec._is_sha256(
            selection_contract["selection_decision_sha256"])
        or selection_contract["selected_survivor"] != survivor
        or _validate_production_blockers(
            selection_contract[
                "selected_candidate_production_blockers"],
            survivor,
        ) is not selection_contract[
            "selected_candidate_production_blockers"]
    ):
        raise CampaignError("holdout selection contract is invalid")

    jobs = build_pre_cpu_jobs(phase, survivor)
    if manifest.get("pre_cpu_jobs") != list(jobs):
        raise CampaignError("campaign manifest population changed")
    assignments = manifest.get("assignments")
    if (
        not isinstance(assignments, list)
        or len(assignments) != len(jobs)
        or any(not isinstance(item, dict) for item in assignments)
        or any(type(item.get("job_id")) is not str for item in assignments)
        or len({item.get("job_id") for item in assignments}) != len(jobs)
    ):
        raise CampaignError("campaign assignments are incomplete")
    for index, (assignment, job) in enumerate(zip(assignments, jobs)):
        _parent_signal_safe_point()
        if (
            not isinstance(assignment, dict)
            or set(assignment) != {
                "order", "job_id", "cpu", "job_file", "output", "log"
            }
            or assignment.get("order") != index
            or assignment.get("job_id") != job["job_id"]
            or assignment.get("cpu")
                != manifest["worker_cpus"][index % manifest["workers"]]
            or assignment.get("job_file") != f"jobs/{index:04d}.json"
            or assignment.get("output") != f"results/{index:04d}.json.gz"
            or assignment.get("log") != f"logs/{index:04d}.log"
        ):
            raise CampaignError("campaign CPU assignment changed")
    return jobs


def _stream_file_sha256(path):
    digest = hashlib.sha256()
    size = 0
    try:
        with open(path, "rb") as source:
            while True:
                chunk = source.read(1024 * 1024)
                if not chunk:
                    break
                _parent_signal_safe_point()
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise CampaignError(f"could not hash {path}: {error}")
    return digest.hexdigest(), size


def _verify_stored_summary_artifacts(directory, completion):
    summary_value, summary_evidence = read_canonical_gzip_json(
        directory / "summary.json.gz")
    expected_summary = completion.get("summary")
    if (
        not isinstance(expected_summary, dict)
        or not isinstance(summary_value, dict)
        or expected_summary != {
            "path": "summary.json.gz",
            **summary_evidence,
        }
        or summary_value.get("schema") != SUMMARY_SCHEMA
        or summary_value.get("phase") != completion["phase"]
        or summary_value.get("pre_cpu_job_list_sha256")
            != completion["pre_cpu_job_list_sha256"]
    ):
        raise CampaignError("stored summary does not match completion")
    ledger = completion.get("cell_ledger")
    if (
        not isinstance(ledger, dict)
        or set(ledger) != {
            "path", "rows", "compressed_sha256", "uncompressed_sha256",
            "compressed_bytes", "uncompressed_bytes",
        }
        or ledger.get("path") != "cell-ledger.jsonl.gz"
    ):
        raise CampaignError("cell ledger evidence is malformed")
    digest, size = _stream_file_sha256(
        directory / "cell-ledger.jsonl.gz")
    if (
        digest != ledger["compressed_sha256"]
        or size != ledger["compressed_bytes"]
    ):
        raise CampaignError("stored cell ledger does not match completion")


def _validate_runtime_monitor(monitor, manifest):
    if (
        not isinstance(monitor, dict)
        or set(monitor) != {
            "total_checks", "total_checks_while_running", "waves"
        }
        or type(monitor.get("total_checks")) is not int
        or type(monitor.get("total_checks_while_running")) is not int
        or not isinstance(monitor["waves"], list)
    ):
        raise CampaignError("runtime binding monitor is malformed")
    expected_wave_count = (
        len(manifest["assignments"]) + manifest["workers"] - 1
    ) // manifest["workers"]
    if len(monitor["waves"]) != expected_wave_count:
        raise CampaignError("runtime binding monitor lost a wave")
    check_total = running_total = 0
    next_order = 0
    expected_binding = canonical_sha256(manifest["runtime_bindings"])
    for wave in monitor["waves"]:
        _parent_signal_safe_point()
        if (
            not isinstance(wave, dict)
            or set(wave) != {
                "first_order", "last_order", "checks",
                "checks_while_running", "first_check_monotonic_ns",
                "last_check_monotonic_ns", "runtime_bindings_sha256",
            }
            or type(wave["first_order"]) is not int
            or wave["first_order"] != next_order
            or type(wave["last_order"]) is not int
            or wave["last_order"] < wave["first_order"]
            or wave["last_order"] >= len(manifest["assignments"])
            or wave["last_order"] - wave["first_order"] + 1
                > manifest["workers"]
            or type(wave["checks"]) is not int
            or wave["checks"] < 1
            or type(wave["checks_while_running"]) is not int
            or not 1 <= wave["checks_while_running"] <= wave["checks"]
            or type(wave["first_check_monotonic_ns"]) is not int
            or type(wave["last_check_monotonic_ns"]) is not int
            or not 0 < wave["first_check_monotonic_ns"]
                <= wave["last_check_monotonic_ns"]
            or wave["runtime_bindings_sha256"] != expected_binding
        ):
            raise CampaignError("runtime binding wave evidence is invalid")
        next_order = wave["last_order"] + 1
        check_total += wave["checks"]
        running_total += wave["checks_while_running"]
    if (
        next_order != len(manifest["assignments"])
        or monitor["total_checks"] != check_total
        or monitor["total_checks_while_running"] != running_total
    ):
        raise CampaignError("runtime binding monitor totals are invalid")


def verify_campaign(directory):
    """Replay a completed campaign and reproduce both summary hashes."""
    directory = Path(directory).resolve()
    manifest, manifest_sha256 = read_canonical_json(
        directory / "manifest.json")
    _validate_manifest(manifest)
    if manifest["phase"] == "holdout":
        reopened = _load_selection_contract(
            manifest["selection_contract"]["path"],
            manifest["survivor"],
        )
        if reopened != manifest["selection_contract"]:
            raise CampaignError(
                "holdout parent selection evidence changed")
    completion, completion_sha256 = read_canonical_json(
        directory / "completion.json")
    try:
        checksum = (directory / "completion.sha256").read_text(
            encoding="ascii")
    except OSError as error:
        raise CampaignError(f"could not read completion checksum: {error}")
    if checksum != f"{completion_sha256}  completion.json\n":
        raise CampaignError("completion checksum does not match")
    completion_fields = {
        "schema", "phase", "survivor", "manifest_sha256",
        "frozen_roster_sha256", "pre_cpu_job_list_sha256", "jobs",
        "runtime_binding_monitor", "runtime_bindings_final",
        "job_files", "log_files", "result_files", "cell_ledger", "summary",
        "completed_unix_ns",
    }
    if (
        not isinstance(completion, dict)
        or set(completion) != completion_fields
        or completion.get("schema") != COMPLETION_SCHEMA
        or completion.get("manifest_sha256") != manifest_sha256
        or completion.get("phase") != manifest["phase"]
        or completion.get("survivor") != manifest["survivor"]
        or completion.get("frozen_roster_sha256")
            != manifest["frozen_roster_sha256"]
        or completion.get("pre_cpu_job_list_sha256")
            != manifest["pre_cpu_job_list_sha256"]
        or completion.get("jobs") != len(manifest["pre_cpu_jobs"])
        or completion.get("runtime_bindings_final")
            != manifest["runtime_bindings"]
        or type(completion.get("completed_unix_ns")) is not int
        or completion["completed_unix_ns"] < manifest["created_unix_ns"]
    ):
        raise CampaignError("completion does not match manifest")
    _validate_runtime_monitor(
        completion.get("runtime_binding_monitor"), manifest)
    _verify_stored_summary_artifacts(directory, completion)
    check_runtime_bindings(manifest["runtime_bindings"], full_hash=True)
    with tempfile.TemporaryDirectory(
            prefix="wh2-za5v-verify-") as temporary:
        reproduced = build_summaries(
            directory, manifest, Path(temporary))
    if (
        reproduced["job_files"] != completion["job_files"]
        or reproduced["result_files"] != completion["result_files"]
        or reproduced["cell_ledger"] != completion["cell_ledger"]
        or reproduced["summary"] != completion["summary"]
        or _verify_empty_logs(directory, manifest) != completion["log_files"]
    ):
        raise CampaignError("strict replay did not reproduce summary hashes")
    return {
        "manifest_sha256": manifest_sha256,
        "completion_sha256": completion_sha256,
        "jobs": len(manifest["pre_cpu_jobs"]),
        "raw_cells":
            len(manifest["pre_cpu_jobs"])
            * _phase_parameters(manifest["phase"])["replicates"],
    }


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    plan_parser = subparsers.add_parser(
        "plan", help="print the frozen pre-CPU phase plan")
    plan_parser.add_argument(
        "--phase", required=True, choices=("selection", "holdout"))
    plan_parser.add_argument(
        "--survivor", choices=tuple(CANDIDATE_BY_NAME))

    run_parser = subparsers.add_parser(
        "run", help="execute one frozen campaign phase")
    run_parser.add_argument(
        "--phase", required=True, choices=("selection", "holdout"))
    run_parser.add_argument("--directory", required=True, type=Path)
    run_parser.add_argument("--bench", required=True, type=Path)
    run_parser.add_argument("--thermal-source", required=True, type=Path)
    run_parser.add_argument(
        "--workers", type=int, default=MAX_WORKERS)
    run_parser.add_argument(
        "--survivor", choices=tuple(CANDIDATE_BY_NAME))
    run_parser.add_argument("--selection-manifest", type=Path)

    verify_parser = subparsers.add_parser(
        "verify", help="strictly replay a completed campaign")
    verify_parser.add_argument("--directory", required=True, type=Path)

    worker_parser = subparsers.add_parser("_worker")
    worker_parser.add_argument("--job-file", required=True, type=Path)
    worker_parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if args.command == "plan":
        if args.phase == "selection" and args.survivor is not None:
            raise CampaignError("selection plan forbids --survivor")
        if args.phase == "holdout" and args.survivor is None:
            raise CampaignError("holdout plan requires --survivor")
        print(canonical_json_bytes(
            make_plan(args.phase, args.survivor)
        ).decode("ascii"), end="")
    elif args.command == "run":
        run_campaign(
            args.phase,
            args.directory,
            args.bench,
            args.thermal_source,
            survivor=args.survivor,
            selection_manifest=args.selection_manifest,
            workers=args.workers,
        )
    elif args.command == "verify":
        print(canonical_json_bytes(
            verify_campaign(args.directory)
        ).decode("ascii"), end="")
    else:
        run_worker(args.job_file, args.output)


if __name__ == "__main__":
    try:
        main()
    except (CampaignError, peel_codec.MeasurementError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
